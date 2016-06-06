/*

    ________/\\\\\\\\\__/\\\\\\\\\\\\\____/\\\_____________        
     _____/\\\////////__\/\\\/////////\\\_\/\\\_____________       
      ___/\\\/___________\/\\\_______\/\\\_\/\\\_____________      
       __/\\\_____________\/\\\\\\\\\\\\\/__\/\\\_____________     
        _\/\\\_____________\/\\\/////////____\/\\\_____________    
         _\//\\\____________\/\\\_____________\/\\\_____________   
          __\///\\\__________\/\\\_____________\/\\\_____________  
           ____\////\\\\\\\\\_\/\\\_____________\/\\\\\\\\\\\\\\\_ 
            _______\/////////__\///______________\///////////////__


                         C P L  -  L I B R A R Y 

           Copyright (C) 2012-2015 Edward Smith & David Trevelyan

License

    This file is part of CPL-Library.

    CPL-Library is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    CPL-Library is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with CPL-Library.  If not, see <http://www.gnu.org/licenses/>.

Description

    See CPLSocket.H 

*/
#include "CPLSocket.H"
#include "blockMesh.H"
#include "PstreamGlobals.H"


// Initialise CFD realm communicator
void CPLSocket::initComms (int& argc, char**& argv)
{

    MPI_Init(&argc, &argv);

    // Split MPI_COMM_WORLD to md and cfd realm comms,
    // store cfd realm comm in global variable for Pstream and
    // store rank in realm

    CPL::init (CPL::cfd_realm, realmComm);
    Foam::PstreamGlobals::CPLRealmComm = realmComm;
    MPI_Comm_rank (realmComm, &rankRealm);

}



// Analyse mesh topology and perform CFD-side CPL_init.
void CPLSocket::\
    initCFD (const Foam::Time &runTime, const Foam::fvMesh &mesh) {

    Foam::Info << "CPLSocket: Analysing processor and mesh topology"
               << Foam::endl;

    // Read from decomposePar dictionary the number of processors in each
    // direction. Must be decomposed with the "simple" method.
    Foam::IOdictionary decomposeDict
    (
        Foam::IOobject ("decomposeParDict", runTime.time().system(), runTime,
                        IOobject::MUST_READ, IOobject::NO_WRITE, false)
    );
    Foam::dictionary simpleCoeffs = decomposeDict.subDict ("simpleCoeffs");
    Foam::Vector<int> np = simpleCoeffs.lookup ("n");
    nprocs = np.x() * np.y() * np.z();

    // Define arrays needed by MPI & CPL cart create routines 
    npxyz[0] = np.x();
    npxyz[1] = np.y();
    npxyz[2] = np.z();
    periods[0] = 1;
    periods[1] = 0;
    periods[2] = 1;

    // Get info needed to calculate this processor's spatial coordinates
    Foam::pointField points = mesh.points();
    Foam::Vector<double> procMaxPoints (maxPoints (points));
    Foam::Vector<double> procMinPoints (minPoints (points));

    // Reduce all with min/max operator to find global min/max
    Foam::Vector<double> globMaxPoints (procMaxPoints);
    Foam::Vector<double> globMinPoints (procMinPoints);
    Foam::reduce (globMaxPoints, maxOp<Foam::Vector<double>>());
    Foam::reduce (globMinPoints, minOp<Foam::Vector<double>>());

    Foam::Vector<double> domainLength (globMaxPoints - globMinPoints);

//    std::cout << "domainlength " << domainLength.x() << " " << domainLength.y() << " " << domainLength.z() << std::endl;
//    std::cout << "proc_min" << procMinPoints.x() << " " << procMinPoints.y() << " " << procMinPoints.z() \
    << "glob_min" << globMinPoints.x() << " " << globMinPoints.y() << " " << globMinPoints.z() << std::endl;
   
    // Store this processor's coordinates in space
    myCoords.push_back
    (
        nint (np.x() * ((procMinPoints.x() - globMinPoints.x())
                       / domainLength.x()))
    );
    myCoords.push_back
    (
        nint (np.y() * ((procMinPoints.y() - globMinPoints.y())
                       / domainLength.y()))
    );
    myCoords.push_back
    (
        nint (np.z() * ((procMinPoints.z() - globMinPoints.z())
                       / domainLength.z()))
    );


    Foam::Info << "CPLSocket: Defining new MPI Cartesian communicator"
               << Foam::endl;

    // Create custom cartesian communicator (cartComm) based on myCoords
    //std::cout << "npxyz: " << npxyz[0] << " " << npxyz[1] << " " << npxyz[2] << "my_coords: " << myCoords[0] << " " << myCoords[1] << " " << myCoords[2] << std::endl;
    CPL::Cart_create (realmComm, 3, npxyz, periods, 
                      myCoords.data(), &cartComm);

    MPI_Comm_rank (cartComm, &rankCart);


    // Prepare inputs for CPL::cfd_init 
    double dt_cfd = runTime.deltaTValue();
    int nsteps = nint ((runTime.endTime().value() - \
                        runTime.startTime().value()) / dt_cfd);

    // Domain dimensions
    xyzL[0] = domainLength.x();
    xyzL[1] = domainLength.y();
    xyzL[2] = domainLength.z();
  
    // dummy density for now TODO(djt06@ic.ac.uk) remove density from coupler
    // cfd_init function input
    double dummyDensity = -666.0;

    Foam::IOdictionary blockMeshDict
    (
        Foam::IOobject ("blockMeshDict", "../constant/polyMesh", runTime,
                        IOobject::MUST_READ, IOobject::NO_WRITE, false)
    );

    Foam::word dummyRegionName("dummy");
    Foam::blockMesh blocks(blockMeshDict, dummyRegionName);
    Foam::Vector<int> meshDensity = blocks[0].meshDensity();
   
    // Global number of cells
    ncxyz[0] = meshDensity.x();
    ncxyz[1] = meshDensity.y();
    ncxyz[2] = meshDensity.z();

    // Origin of the domain
    double xyz_orig[3] = {0.0, 0.0, 0.0};

    // Initialise CPL library
    CPL::setup_cfd (nsteps, dt_cfd, cartComm, xyzL, xyz_orig, ncxyz, dummyDensity);

    getCellTopology();
    allocateBuffers();

    // Store some values from CPL that are useful later
    CPLDensity = CPL::density_cfd();

    Foam::Info << "OpenFOAM CPL topology initialisation complete" << Foam::endl;
    return;

}

void CPLSocket::\                                                                                                                                                                                                                      
getCellTopology() {                                                                                                                                                                                                                       


    // Cell sizes
    dx = CPL::get<double> ("dx");
    dy = CPL::get<double> ("dy");
    dz = CPL::get<double> ("dz");
   
    // Cell bounds for the overlap region
    CPL::get_olap_limits(olapRegion.data());
    
    // Cell bounds for velocity BCs region
    velBCRegion = olapRegion;
    velBCRegion[3] = velBCRegion[2];
    CPL::my_proc_portion(velBCRegion.data(), velBCPortion.data());

    CPL::get_no_cells(velBCPortion.data(), velBCCells);

    // Cell bounds for the constrained region
    CPL::get_cnst_limits(cnstFRegion.data());
    CPL::my_proc_portion(cnstFRegion.data(), cnstFPortion.data());

    CPL::get_no_cells(cnstFPortion.data(), cnstFCells);
}

void CPLSocket::\
allocateBuffers() {
    // Received stress field
    int zeroShapeStress[4] = {9, 0, 0, 0};
    int sendShape[4] = {9, cnstFCells[0], cnstFCells[1], cnstFCells[2]};
    sendStressBuff.resize(4, sendShape);
    recvStressBuff.resize(4, zeroShapeStress);

    // LAMMPS computed velocity field
    int recvShape[4] = {4, velBCCells[0], velBCCells[1], velBCCells[2]};
    int zeroShapeVel[4] = {4, 0, 0, 0};
    recvVelocityBuff.resize (4, recvShape);
    sendVelocityBuff.resize (4, zeroShapeVel);
}  
    
// Packs the 9 components of the stress-tensor to the socket's CPL::ndArray
// storage.
void CPLSocket::packStress(volVectorField &U, dimensionedScalar &nu, fvMesh &mesh)
{

    // Evaluate the stress tensor sigma at all local cells (forget pressure for
    // now)
    Foam::dimensionedScalar mu(CPLDensity*nu);
    Foam::volTensorField gradU(fvc::grad(U));
    Foam::volTensorField sigma(mu*(gradU + gradU.T()));

    // Loop over socket cells
    int icmin = cnstFPortion[0];
    int jcmin = cnstFPortion[2];
    int kcmin = cnstFPortion[4];

    for (int ix = 0; ix < sendStressBuff.shape(1); ix++) {
        for (int iy = 0; iy < sendStressBuff.shape(2); iy++) {
            for (int iz = 0; iz < sendStressBuff.shape(3); iz++) {

                // Global position at cell center
                Foam::point globalPos
                (
                    (static_cast<double>(ix + icmin) + 0.5) * dx,
                    (static_cast<double>(iy + jcmin) + 0.5) * dy,
                    (static_cast<double>(iz + kcmin) + 0.5) * dz
                );

                // Get value of stress 9D by interpolating
                Foam::label cell = mesh.findCell(globalPos);
                sendStressBuff(0,ix,iy,iz) = sigma[cell].xx();
                sendStressBuff(1,ix,iy,iz) = sigma[cell].xy();
                sendStressBuff(2,ix,iy,iz) = sigma[cell].xz();
                sendStressBuff(3,ix,iy,iz) = sigma[cell].yx();
                sendStressBuff(4,ix,iy,iz) = sigma[cell].yy();
                sendStressBuff(5,ix,iy,iz) = sigma[cell].yz();
                sendStressBuff(6,ix,iy,iz) = sigma[cell].zx();
                sendStressBuff(7,ix,iy,iz) = sigma[cell].zy();
                sendStressBuff(8,ix,iy,iz) = sigma[cell].zz();

            }
        }
    }

}

// Unpacks the 3 components of the velocity-tensor from the socket's
// recvVelocity (cpl::ndArray) storage.
double CPLSocket::\
unpackVelocity(volVectorField &U, fvMesh &mesh) {

    // Take the mean across spanwise direction if specified by coupler
    // y direction is axis #2
    if (CPL::get<int>("cpl_cfd_bc_slice")) {

        Foam::Info << "CPL_CFD_BC_SLICE is on: averaging CFD recvVelocity "
                      "in the x-z plane" << Foam::endl;

        // Number of cells in the local processor in x-z plane
        int N = recvVelocityBuff.shape(1) * recvVelocityBuff.shape(3);

        // For every component and y-value 
        for (int j = 0; j < recvVelocityBuff.shape(2); ++j) {
            for (int c = 0; c < recvVelocityBuff.shape(0); ++c) {
            // Sum across the x-z plane 
            double total = 0.0;
                for (int k = 0; k < recvVelocityBuff.shape(3); ++k)
                    for (int i = 0; i < recvVelocityBuff.shape(1); ++i)
                        total += recvVelocityBuff(c, i, j, k);

            // Find mean by dividing sum by number of cells 
                for (int k = 0; k < recvVelocityBuff.shape(3); ++k)
                    for (int i = 0; i < recvVelocityBuff.shape(1); ++i)
                        recvVelocityBuff(c, i, j, k) = total / static_cast<double> (N);
            }
        }
    } 

    int applyBCx = CPL::get<int> ("cpl_cfd_bc_x");
    int applyBCy = CPL::get<int> ("cpl_cfd_bc_y");
    int applyBCz = CPL::get<int> ("cpl_cfd_bc_z");

    Foam::string receivePatchName ("CPLReceiveMD");
    Foam::label rvPatchID = mesh.boundary().findPatchID(receivePatchName);

    if (rvPatchID == -1) {
        FatalErrorIn ( "CPLSocket::unpack()")
            << " Could not find patch ID " << receivePatchName << ". "
               " Aborting."
            << exit(FatalError);
    }

    Foam::fvPatchVectorField& rvPatch = U.boundaryField()[rvPatchID];
    const Foam::vectorField faceCenters = mesh.boundary()[rvPatchID].Cf();

    for (int faceI = 0; faceI != faceCenters.size(); ++faceI) {

        double facex = faceCenters[faceI].x();
        double facey = faceCenters[faceI].y();
        double facez = faceCenters[faceI].z();

        // Find the cell indices for this position recvVelocity(:, ix, iy, iz)
        int glob_cell[3];
        CPL::map_coord2cell(facex, facey, facez, glob_cell);
        int loc_cell[3];
        bool valid_cell = CPL::map_glob2loc_cell(velBCPortion.data(), glob_cell, loc_cell);

        // For now, the coupler passes momentum and mass in velocity array
        // This is due for review. 
        if (valid_cell) {
            //std::cout << "cell: " << loc_cell[0] << " " << loc_cell[1] << " " << loc_cell[2] \
               << "limits: " << velBCPortion[0] << " " << velBCPortion[1] << " " << velBCPortion[2] << " " << velBCPortion[3] << " " << velBCPortion[4] << " " << velBCPortion[5] << std::endl;

            double m = recvVelocityBuff(3, loc_cell[0], loc_cell[1], loc_cell[2]); 
            double recvvy = recvVelocityBuff(1, loc_cell[0], loc_cell[1], loc_cell[2]) /m;
            double recvvx = recvVelocityBuff(0, loc_cell[0], loc_cell[1], loc_cell[2])/m;
            double recvvz = recvVelocityBuff(2, loc_cell[0], loc_cell[1], loc_cell[2])/m;

            // Received velocity was averaged in the cell region BELOW 
            // the domain. So, to get the velocity required at the BOUNDARY, 
            // we need to take the average of the velocities ABOVE (i.e. the 
            // bottom OpenFOAM cell velocity) AND BELOW (i.e. the received 
            // velocity) the domain face.
            Foam::point closestCellCentre (facex, facey + 0.5*dy, facez);
            Foam::label cell = mesh.findCell (closestCellCentre);
            double vx = (recvvx + U[cell].x()) / 2.0;
            double vy = (recvvy + U[cell].y()) / 2.0;
            double vz = (recvvz + U[cell].z()) / 2.0;

            if (applyBCx) rvPatch[faceI].x() = vx;
            if (applyBCy) rvPatch[faceI].y() = vy;
            if (applyBCz) rvPatch[faceI].z() = vz;
            //return recvvy + m + recvvx + recvvz;
         }
        }

    if (applyBCx) Foam::Info << "MD->CFD BC x-velocity applied." << Foam::endl;
    if (applyBCy) Foam::Info << "MD->CFD BC y-velocity applied." << Foam::endl;
    if (applyBCz) Foam::Info << "MD->CFD BC z-velocity applied." << Foam::endl;

}

// Sends 9 components of the stress-tensor to overlapping MD processes.
void CPLSocket::\
sendStress() {

    const int comm_style = CPL::get<int> ("comm_style"); 
    const int gath_scat = CPL::get<int> ("comm_style_gath_scat"); 
    const int send_recv = CPL::get<int> ("comm_style_send_recv"); 


    if (comm_style == gath_scat) {
        // Send stress from CFD to MD processes
     /**
    for (int ix = 0; ix < sendStressBuff.shape(1); ix++) {
        for (int iy = 0; iy < sendStressBuff.shape(2); iy++) {
            for (int iz = 0; iz < sendStressBuff.shape(3); iz++) {

                if (cnstFPortion[2] >= 0)
                    std::cout << "Stress:" << sendStressBuff(1,ix,iy,iz) << " " << sendStressBuff(4,ix,iy,iz) << " " <<sendStressBuff(7,ix,iy,iz) << std::endl;
             }
         }
     }**/

        CPL::scatter (sendStressBuff.data(), sendStressBuff.shapeData(),
                      cnstFRegion.data(), recvStressBuff.data(), 
                      recvStressBuff.shapeData());
    }
    else if (comm_style == send_recv) {
        FatalErrorIn ("CPLSocket::send()")
            << " invalid comm_style."
               " Currently only gather/scatter supported by OpenFOAM."
            << exit(FatalError);
    }
    else {
        FatalErrorIn ("CPLSocket::send()")
            << " unrecognised comm_style."
               " Currently only gather/scatter supported by OpenFOAM."
            << exit(FatalError);
    }
  
}

// Receives 3 components of the velocity vector from overlapping MD processes.
void CPLSocket::recvVelocity()
{

    const int comm_style = CPL::get<int> ("comm_style"); 
    const int gath_scat = CPL::get<int> ("comm_style_gath_scat"); 
    const int send_recv = CPL::get<int> ("comm_style_send_recv"); 

    if (comm_style == gath_scat) {
        // Receive velocity from CFD to MD processes
        CPL::gather (sendVelocityBuff.data(), sendVelocityBuff.shapeData(),
                     velBCRegion.data(), recvVelocityBuff.data(), 
                     recvVelocityBuff.shapeData());
   /**
    for (int ix = 0; ix < recvVelocityBuff.shape(1); ix++) {
        for (int iy = 0; iy < recvVelocityBuff.shape(2); iy++) {
            for (int iz = 0; iz < recvVelocityBuff.shape(3); iz++) {
                if (cnstFPortion[2] >= 0)
                    std::cout << "Stress:" << recvVelocityBuff(0,ix,iy,iz) << " " << recvVelocityBuff(1,ix,iy,iz) << " " <<recvVelocityBuff(2,ix,iy,iz) << std::endl;
             }
         }
    }
**/ 
    }
    else if (comm_style == send_recv) {
        FatalErrorIn ("CPLSocket::receive()")
            << " invalid comm_style."
               " Currently only gather/scatter supported by OpenFOAM."
            << exit(FatalError);
    }
    else {
        FatalErrorIn ("CPLSocket::receive()")
            << " unrecognised comm_style."
               " Currently only gather/scatter supported by OpenFOAM."
            << exit(FatalError);
    }
  
}
