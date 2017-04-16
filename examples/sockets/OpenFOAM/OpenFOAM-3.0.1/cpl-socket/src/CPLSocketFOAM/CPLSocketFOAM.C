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

           Copyright (C) 2012-2017 Edward Smith & David Trevelyan

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

    See CPLSocketFOAM.H 

*/
#include "CPLSocketFOAM.H"
#include "blockMesh.H"
#include "PstreamGlobals.H"


// Initialise CFD realm communicator
void CPLSocketFOAM::initComms (int& argc, char**& argv)
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
void CPLSocketFOAM::initCFD (const Foam::Time &runTime, 
                             const Foam::fvMesh &mesh) {

    Foam::Info << "CPLSocketFOAM: Analysing processor and mesh topology"
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


    Foam::Info << "CPLSocketFOAM: Defining new MPI Cartesian communicator"
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
    CPL::setup_cfd (cartComm, xyzL, xyz_orig, ncxyz);

    getCellTopology();
    //allocateBuffers(1);

    // Store some values from CPL that are useful later
    CPLDensity = CPL::density_cfd();

    Foam::Info << "OpenFOAM CPL topology initialisation complete" << Foam::endl;
    return;

}

void CPLSocketFOAM::getCellTopology() {                                                                                                                                                                                                                       

    // Cell sizes
    dx = CPL::get<double> ("dx");
    dy = CPL::get<double> ("dy");
    dz = CPL::get<double> ("dz");

    // Get overlap extents
    CPL::get_olap_limits(olapRegion.data());
   
    // Processor cell bounds for the overlap region
    CPL::my_proc_portion(olapRegion.data(), olapPortion.data());
    CPL::get_no_cells(olapPortion.data(), olapCells);
    
    // Processor cell bounds for velocity BCs region
    velBCRegion = olapRegion;
    //velBCRegion[3] = velBCRegion[2];
    CPL::my_proc_portion(velBCRegion.data(), velBCPortion.data());
    CPL::get_no_cells(velBCPortion.data(), velBCCells);

    // Processor cell bounds for the constrained region
    CPL::get_cnst_limits(cnstFRegion.data());
    CPL::my_proc_portion(cnstFRegion.data(), cnstFPortion.data());
    CPL::get_no_cells(cnstFPortion.data(), cnstFCells);
}

void CPLSocketFOAM::allocateBuffers(int sendtype) {

    //Check what is to be packed and sent

    // 3,9 or 12-component for every local cell
    int sendShape[4] = {0, cnstFCells[0], cnstFCells[1], cnstFCells[2]};

    if (sendtype == PACKVELONLY)
    {
        sendShape[0] = 3;
    }
    else if (sendtype == PACKSTRESSONLY)
    {
        sendShape[0] = 9;
    }
    else if (sendtype == PACKVELSTRESS)
    {
        sendShape[0] = 12;
    }
    else if (sendtype == DEBUG)
        sendShape[0] = 3;
    else
    {
        FatalErrorIn
        (
            "CPLSocketFOAM::pack()"
        )
            << " sendtype flag not of known type " << sendtype << ". "
               " Aborting."
            << exit(FatalError);
    }

    sendBuf.resize(4, sendShape);

    // LAMMPS computed velocity field
    int recvVelocityShape[4] = {4, velBCCells[0], velBCCells[1], velBCCells[2]};
    recvVelocityBuff.resize (4, recvVelocityShape);

    // LAMMPS olap size field
    int recvShape[4] = {4, olapCells[0], olapCells[1], olapCells[2]};
    recvBuf.resize (4, recvShape);

}


// Packs the components of velocity and stress-tensor to the socket's CPL::ndArray
// storage.
void CPLSocketFOAM::pack(volVectorField &U, 
                     dimensionedScalar &nu, 
                     fvMesh &mesh, 
                     int sendtype)
{

    // Evaluate the stress tensor sigma at all local cells 
    // (forget pressure for now)
    Foam::dimensionedScalar mu(CPLDensity*nu);
    Foam::volTensorField gradU(fvc::grad(U));
    Foam::volTensorField sigma(mu*(gradU + gradU.T()));

    //Reallocate buffer depending on send type
    allocateBuffers(sendtype);

    //printf("pack %d %d %d %d %d %d %d %d \n", sendShape[0], sendShape[1], sendShape[2], sendShape[3], 
    //                                          sendBuf.shape(0), sendBuf.shape(1), sendBuf.shape(2), sendBuf.shape(3));
    
    // Loop over socket cells, -1 for Fortran to C++ indexing
    int icmin = cnstFPortion[0];
    int jcmin = cnstFPortion[2];
    int kcmin = cnstFPortion[4];

    for (int ix=0; ix<sendBuf.shape(1); ix++) {
        for (int iy=0; iy<sendBuf.shape(2); iy++) {
            for (int iz=0; iz<sendBuf.shape(3); iz++) {

                // Global position at cell center
                Foam::point globalPos
                (
                    (static_cast<double>(ix + icmin) + 0.5) * dx,
                    (static_cast<double>(iy + jcmin) + 0.5) * dy,
                    (static_cast<double>(iz + kcmin) + 0.5) * dz

                );

                Foam::label cell = mesh.findCell(globalPos);
                if (sendtype == PACKVELONLY)
                {
                    // Get value of velocity 3D
                    sendBuf(0,ix,iy,iz) = U[cell].x();
                    sendBuf(1,ix,iy,iz) = U[cell].y();
                    sendBuf(2,ix,iy,iz) = U[cell].z();

//                    printf("vel %d %d %d %4.2f %4.2f %4.2f  \n", ix,iy,iz, 
//                            sendBuf(0,ix,iy,iz), sendBuf(1,ix,iy,iz), sendBuf(2,ix,iy,iz));
                }
                else if (sendtype == PACKSTRESSONLY)
                {
                    // Get value of stress 9D by interpolating
                    sendBuf(0,ix,iy,iz) = sigma[cell].xx();
                    sendBuf(1,ix,iy,iz) = sigma[cell].xy();
                    sendBuf(2,ix,iy,iz) = sigma[cell].xz();
                    sendBuf(3,ix,iy,iz) = sigma[cell].yx();
                    sendBuf(4,ix,iy,iz) = sigma[cell].yy();
                    sendBuf(5,ix,iy,iz) = sigma[cell].yz();
                    sendBuf(6,ix,iy,iz) = sigma[cell].zx();
                    sendBuf(7,ix,iy,iz) = sigma[cell].zy();
                    sendBuf(8,ix,iy,iz) = sigma[cell].zz();

//                    printf("stress %d %d %d %4.2f %4.2f %4.2f %4.2f %4.2f %4.2f %4.2f %4.2f %4.2f %4.2f \n", ix,iy,iz, 
//                            sendBuf(0,ix,iy,iz), sendBuf(1,ix,iy,iz), sendBuf(2,ix,iy,iz), 
//                            sendBuf(3,ix,iy,iz), sendBuf(4,ix,iy,iz), sendBuf(5,ix,iy,iz), 
//                            sendBuf(6,ix,iy,iz), sendBuf(7,ix,iy,iz), sendBuf(8,ix,iy,iz));

                }
                else if (sendtype == PACKVELSTRESS)
                {
                    // Get value of velocity 3D
                    sendBuf(0,ix,iy,iz) = U[cell].x();
                    sendBuf(1,ix,iy,iz) = U[cell].y();
                    sendBuf(2,ix,iy,iz) = U[cell].z();

                    // Get value of stress 9D by interpolating
                    sendBuf(3,ix,iy,iz) = sigma[cell].xx();
                    sendBuf(4,ix,iy,iz) = sigma[cell].xy();
                    sendBuf(5,ix,iy,iz) = sigma[cell].xz();
                    sendBuf(6,ix,iy,iz) = sigma[cell].yx();
                    sendBuf(7,ix,iy,iz) = sigma[cell].yy();
                    sendBuf(8,ix,iy,iz) = sigma[cell].yz();
                    sendBuf(9,ix,iy,iz) = sigma[cell].zx();
                    sendBuf(10,ix,iy,iz) = sigma[cell].zy();
                    sendBuf(11,ix,iy,iz) = sigma[cell].zz();
                }
                else if (sendtype == DEBUG)
                {
                    sendBuf(0,ix,iy,iz) = globalPos[0];
                    sendBuf(1,ix,iy,iz) = globalPos[1];
                    sendBuf(2,ix,iy,iz) = globalPos[2];
                }

            }
        }
    }

}


//// Unpacks scalar from socket
//double CPLSocketFOAM::unpack(volScalarField &s, fvMesh &mesh) {

//    for (int ix=0; ix<recvBuf.shape(1); ix++) {
//        for (int iy=0; iy<recvBuf.shape(2); iy++) {
//            for (int iz=0; iz<recvBuf.shape(3); iz++) {

//                double e  = recvBuf(0, ix, iy, iz);
//                double glob_pos[3];
//                CPL::map_cell2coord(ix, iy, iz, glob_pos);
//                Foam::point closestCellCentre (glob_pos[0], glob_pos[1], glob_pos[2]);
//                Foam::label cell = mesh.findCell (closestCellCentre);

//                s[cell] = e;

//            }
//        }
//    }

//}


//// Unpacks vector from socket
//double CPLSocketFOAM::unpack(volVectorField &F, fvMesh &mesh) {

//    for (int ix=0; ix<recvBuf.shape(1); ix++) {
//        for (int iy=0; iy<recvBuf.shape(2); iy++) {
//            for (int iz=0; iz<recvBuf.shape(3); iz++) {

//                double Fx = recvBuf(0, ix, iy, iz);
//                double Fy = recvBuf(1, ix, iy, iz);
//                double Fz = recvBuf(2, ix, iy, iz);

//                double glob_pos[3];
//                CPL::map_cell2coord(ix, iy, iz, glob_pos);
//                Foam::point closestCellCentre (glob_pos[0], glob_pos[1], glob_pos[2]);
//                Foam::label cell = mesh.findCell (closestCellCentre);

//                F[cell].x() = Fx;
//                F[cell].y() = Fy;
//                F[cell].z() = Fz;

//            }
//        }
//    }

//}


// Unpacks the components from the socket's
double CPLSocketFOAM::unpackPorousForce(volVectorField &F, volScalarField &eps, fvMesh &mesh) {

    for (int ix=0; ix<recvBuf.shape(1); ix++) {
        for (int iy=0; iy<recvBuf.shape(2); iy++) {
            for (int iz=0; iz<recvBuf.shape(3); iz++) {

                double Fx = recvBuf(0, ix, iy, iz);
                double Fy = recvBuf(1, ix, iy, iz);
                double Fz = recvBuf(2, ix, iy, iz);
                double e  = recvBuf(3, ix, iy, iz);

                double glob_pos[3];
                CPL::map_cell2coord(ix, iy, iz, glob_pos);
                Foam::point closestCellCentre(glob_pos[0]+0.5*dx, glob_pos[1]+0.5*dy, glob_pos[2]+0.5*dz);
                Foam::label cell = mesh.findCell(closestCellCentre);

//                Foam::Info << " cell " << ix << " " << iy << " " << iz << " " 
//                           << glob_pos[0]+0.5*dx << " " 
//                           << glob_pos[1]+0.5*dy << " " 
//                           << glob_pos[2]+0.5*dz << " " 
//                           << mesh.C()[cell] << " " << Foam::endl;

                eps[cell] = e;
                F[cell].x() = Fx;
                F[cell].y() = Fy;
                F[cell].z() = Fz;

            }
        }
    }


//    forAll(mesh.C(), cell)
//    {
//        Foam::Info << rankRealm << " " << cell << " eps " 
//               << eps[cell] << " "
//          << mesh.C()[cell] << " Fx,y,z " 
//                 << F[cell].x() << " " << 
//                    F[cell].y() << " " << 
//                    F[cell].z() << Foam::endl;

//    }

}


// Unpacks the 3 components of the velocity-tensor from the socket's
// recvVelocity (cpl::ndArray) storage into a boundary condition.
double CPLSocketFOAM::unpackVelocity(volVectorField &U, fvMesh &mesh) {

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
        FatalErrorIn ( "CPLSocketFOAM::unpack()")
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
        int glob_cell[3]; int loc_cell[3];
        CPL::map_coord2cell(facex, facey, facez, glob_cell);
        glob_cell[1] += 1; //E.S. Added this line otherwise out of domain error from map_glob2loc_cell!!
        bool valid_cell = CPL::map_glob2loc_cell(velBCPortion.data(), glob_cell, loc_cell);

//        printf("cells %d %d %d %d %d %d %d %d  \n", 
//               rankRealm, valid_cell, loc_cell[0], loc_cell[1], loc_cell[2], 
//               glob_cell[0], glob_cell[1], glob_cell[2]);

        // For now, the coupler passes momentum and mass in velocity array
        // This is due for review.
        if (valid_cell) {
            //std::cout << "cell: " << loc_cell[0] << " " << loc_cell[1] << " " << loc_cell[2] \
               << "limits: " << velBCPortion[0] << " " << velBCPortion[1] << " " << velBCPortion[2] << " " << velBCPortion[3] << " " << velBCPortion[4] << " " << velBCPortion[5] << std::endl;

            double m = recvVelocityBuff(3, loc_cell[0], loc_cell[1], loc_cell[2]); 
            double recvvx = recvVelocityBuff(0, loc_cell[0], loc_cell[1], loc_cell[2])/m;
            double recvvy = recvVelocityBuff(1, loc_cell[0], loc_cell[1], loc_cell[2])/m;
            double recvvz = recvVelocityBuff(2, loc_cell[0], loc_cell[1], loc_cell[2])/m;

            // Received velocity was averaged in the cell region BELOW 
            // the domain. So, to get the velocity required at the BOUNDARY, 
            // we need to take the average of the velocities ABOVE (i.e. the 
            // bottom OpenFOAM cell velocity) AND BELOW (i.e. the received 
            // velocity) the domain face.
            Foam::point closestCellCentre (facex, facey + 0.5*dy, facez);
            Foam::label cell = mesh.findCell (closestCellCentre);
            double vx = recvvx; 
            double vy = recvvy; 
            double vz = recvvz; 
//            double vx = (recvvx + U[cell].x()) / 2.0;
//            double vy = (recvvy + U[cell].y()) / 2.0;
//            double vz = (recvvz + U[cell].z()) / 2.0;

//            printf("vel %d %d %d %d %4.2f %4.2f %4.2f  \n", 
//                   rankRealm, loc_cell[0], loc_cell[1], loc_cell[2], vx, vy, vz);

            if (applyBCx) rvPatch[faceI].x() = vx;
            if (applyBCy) rvPatch[faceI].y() = vy;
            if (applyBCz) rvPatch[faceI].z() = vz;
            //return recvvy + m + recvvx + recvvz;
         }
    }

//    if (applyBCx) Foam::Info << "MD->CFD BC x-velocity applied." << Foam::endl;
//    if (applyBCy) Foam::Info << "MD->CFD BC y-velocity applied." << Foam::endl;
//    if (applyBCz) Foam::Info << "MD->CFD BC z-velocity applied." << Foam::endl;

}


// Sends buffer to overlapping MD processes.
void CPLSocketFOAM::send() {
    CPL::send(sendBuf.data(), sendBuf.shapeData(), cnstFRegion.data());
}

// Receives buffer from overlapping MD processes.
void CPLSocketFOAM::recv()
{
    CPL::recv(recvBuf.data(), recvBuf.shapeData(), olapPortion.data());
}


// Receives 3 components of the velocity vector from overlapping MD processes.
void CPLSocketFOAM::recvVelocity()
{
    CPL::recv(recvVelocityBuff.data(), 
              recvVelocityBuff.shapeData(), 
              velBCPortion.data());
}



