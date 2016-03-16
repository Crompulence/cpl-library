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

// Default constructor, no arguments, does nothing
CPLSocket::CPLSocket() {}


// Initialise CFD realm communicator
void CPLSocket::initComms (int& argc, char**& argv)
{

    MPI_Init(&argc, &argv);

    // Split MPI_COMM_WORLD to md and cfd realm comms,
    // store cfd realm comm in global variable for Pstream and
    // store rank in realm

    CPL::create_comm (CPL::cfd_realm, realmComm);
    Foam::PstreamGlobals::CPLRealmComm = realmComm;
    MPI_Comm_rank (realmComm, &rankRealm);

}



// Analyse mesh topology and perform CFD-side CPL_init.
void CPLSocket::initCFD
(
    const Foam::Time &time,
    const Foam::fvMesh &mesh
)
{

    Foam::Info << "CPLSocket: Analysing processor and mesh topology"
               << Foam::endl;

    // Read from decomposePar dictionary the number of processors in each
    // direction. Must be decomposed with the "simple" method.
    Foam::IOdictionary decomposeDict
    (
        Foam::IOobject
        (
            "decomposeParDict",
            time.time().system(),
            time,
            IOobject::MUST_READ,
            IOobject::NO_WRITE,
            false
        )
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
   
    // Store this processor's coordinates in space (0 indexed) 
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
    CPL::Cart_create
    (
        realmComm,
        ndims,
        npxyz,
        periods,
        myCoords.data(),
        &cartComm
    );
    MPI_Comm_rank (cartComm, &rankCart);


    // Prepare inputs for CPL::cfd_init 
    double dt_cfd = time.deltaTValue();
    int nsteps = nint
    (
        (time.endTime().value() - time.startTime().value()) / dt_cfd
    );

    this->allReduceProcCoords();

    double xyzL[ndims] =
    {
        domainLength.x(),
        domainLength.y(),
        domainLength.z()
    };

    Foam::IOdictionary blockMeshDict
    (
        Foam::IOobject
        (
            "blockMeshDict",
            "../constant/polyMesh",
            time,
            IOobject::MUST_READ,
            IOobject::NO_WRITE,
            false
        )
    );
    Foam::word dummyRegionName("dummy");
    Foam::blockMesh blocks(blockMeshDict, dummyRegionName);
    Foam::Vector<int> meshDensity = blocks[0].meshDensity();
   
    // Global number of cells
    ncxg = meshDensity.x();
    ncyg = meshDensity.y();
    nczg = meshDensity.z();
    int ncxyzg[ndims] = {ncxg, ncyg, nczg};

    int ijkcmin[ndims]; for (int i=0; i<ndims; i++) ijkcmin[i] = 1;
    int ijkcmax[ndims]; for (int i=0; i<ndims; i++) ijkcmax[i] = ncxyzg[i];

    // Number of local cells per processor
    ncxl = ncxg / np.x();
    ncyl = ncyg / np.y();
    nczl = nczg / np.z();
    iTmax.clear(); jTmax.clear(); kTmax.clear();
    iTmin.clear(); jTmin.clear(); kTmin.clear();
    for (int i=0; i<np.x(); i++) iTmax.push_back((i+1)*ncxl);
    for (int j=0; j<np.y(); j++) jTmax.push_back((j+1)*ncyl);
    for (int k=0; k<np.z(); k++) kTmax.push_back((k+1)*nczl);
    for (int i=0; i<np.x(); i++) iTmin.push_back(iTmax[i] - ncxl + 1);
    for (int j=0; j<np.y(); j++) jTmin.push_back(jTmax[j] - ncyl + 1);
    for (int k=0; k<np.z(); k++) kTmin.push_back(kTmax[k] - nczl + 1);

    int tempSizexy[2] = {ncxg+1, ncyg+1};
    CPL::ndArray<double> xg (2, tempSizexy), yg (2, tempSizexy); 
    std::vector<double> zg; // 1D only in z
    
    // CFD grid cell sizes
    dx = domainLength.x() / static_cast<double> (ncxg);
    dy = domainLength.y() / static_cast<double> (ncyg);
    dz = domainLength.z() / static_cast<double> (nczg);
   
    // Create a structured grid to pass to CPL
    for (int i=0; i<ncxg+1; i++)
    {
        for (int j=0; j<ncyg+1; j++)
        {
            xg(i, j) = i*dx;
            yg(i, j) = j*dy;
        }
    }
    
    for (int k=0; k<nczg+1; k++)
    {
        zg.push_back(k*dz);
    }

    // Store cell indices for later
    /*int shape[3] = {ncxl, ncyl, nczl};
    cellIndex.resize(3, shape);
    double x, y, z;
    int celli;
    Foam::Vector<double> cellCentre;
    for (int i=0; i<ncxl; i++)
    {
        x = (static_cast<double>(i+iTmin[myCoords[0]]) + 0.5)*dx;

        for (int j=0; j<ncyl; j++)
        {
            y = (static_cast<double>(j+jTmin[myCoords[1]]) + 0.5)*dy;

            for (int k=0; k<nczl; k++)
            {
                z = (static_cast<double>(k+kTmin[myCoords[2]]) + 0.5)*dz;
               
                cellCentre.x() = x;
                cellCentre.y() = y;
                cellCentre.z() = z;

                Foam::Info << cellCentre << Foam::endl;
                celli = 1;//mesh.findCell(cellCentre);
                if (celli != -1)
                {
                    cellIndex(i, j, k) = celli;
                }
            }
        }
    }*/

    // dummy density for now TODO(djt06@ic.ac.uk) remove density from coupler
    // cfd_init function input
    double dummyDensity = -666.0;

    CPL::cfd_init
    (
        nsteps,
        dt_cfd,
        cartComm,
        allCoords.data(),
        npxyz,
        xyzL,
        ncxyzg,
        dummyDensity,
        ijkcmax,
        ijkcmin,
        iTmin.data(),
        iTmax.data(),
        jTmin.data(),
        jTmax.data(),
        kTmin.data(),
        kTmax.data(),
        xg.data(),
        yg.data(),
        zg.data()
    );

    // Store some values from CPL that are useful later
    CPLDensity = CPL::density_cfd();

    olap_limits.resize(6);
    olap_limits[0] = CPL::get<int> ("icmin_olap");
    olap_limits[1] = CPL::get<int> ("icmax_olap");
    olap_limits[2] = CPL::get<int> ("jcmin_olap");
    olap_limits[3] = CPL::get<int> ("jcmax_olap");
    olap_limits[4] = CPL::get<int> ("kcmin_olap");
    olap_limits[5] = CPL::get<int> ("kcmax_olap");

    Foam::Info << "OpenFOAM CPL topology initialisation complete" << Foam::endl;
    return;

}

void CPLSocket::allReduceProcCoords()
{

    // Resize allCoords to handle nproc by ndims shape and set all elements to
    // zero
    int arrayDims = 2;
    int shape[2] = {ndims, nprocs};
    allCoords.resize (arrayDims, shape);
    allCoords = 0;

    // Set my coordinates in the right place
    // int myRank;
    // MPI_Comm_rank(cartComm, &myRank);
    for (int coord = 0; coord < ndims; ++coord)
    {
        // Fortran counts from 1!
        allCoords (coord, rankCart) = myCoords[coord] + 1;
    }
    
    // Reduce on cartComm
    CPL::ndArray<int> buf = allCoords;
    MPI_Allreduce
    (
        allCoords.data(),
        buf.data(),
        allCoords.size(),
        MPI_INT,
        MPI_SUM,
        cartComm
    );
    allCoords = buf;
    
    return;

}


// Packs the 9 components of the stress-tensor to the socket's CPL::ndArray
// storage.
void CPLSocket::pack(volVectorField &U, dimensionedScalar &nu, fvMesh &mesh)
{

    // Evaluate the stress tensor sigma at all local cells (forget pressure for
    // now)
    Foam::dimensionedScalar mu(CPLDensity*nu);
    Foam::volTensorField gradU(fvc::grad(U));
    Foam::volTensorField sigma(mu*(gradU + gradU.T()));

    // 9-component stress tensor for every local cell, CPL::scatter requires
    // that the data on the "receiving" process is size 0
    int sendShape[4] = {9, ncxl, ncyl, nczl};
    int zeroShape[4] = {0, 0, 0, 0};

    // Clear and reallocate packed data
    sendStress.resize(4, sendShape);
    recvStress.resize(4, zeroShape);
    
    // Loop over socket cells, -1 for Fortran to C++ indexing
    int myiTmin = iTmin[myCoords[0]] - 1;
    int myjTmin = jTmin[myCoords[1]] - 1;
    int mykTmin = kTmin[myCoords[2]] - 1;
    for (int ix=0; ix<sendStress.shape(1); ix++)
    {
        for (int iy=0; iy<sendStress.shape(2); iy++)
        {
            for (int iz=0; iz<sendStress.shape(3); iz++)
            {

                // Global position at cell center
                Foam::point globalPos
                (
                    (static_cast<double>(ix+myiTmin) + 0.5)*dx,
                    (static_cast<double>(iy+myjTmin) + 0.5)*dy,
                    (static_cast<double>(iz+mykTmin) + 0.5)*dz
                );

                // Get value of stress 9D by interpolating
                Foam::label cell = mesh.findCell(globalPos);
                sendStress(0,ix,iy,iz) = sigma[cell].xx();
                sendStress(1,ix,iy,iz) = sigma[cell].xy();
                sendStress(2,ix,iy,iz) = sigma[cell].xz();
                sendStress(3,ix,iy,iz) = sigma[cell].yx();
                sendStress(4,ix,iy,iz) = sigma[cell].yy();
                sendStress(5,ix,iy,iz) = sigma[cell].yz();
                sendStress(6,ix,iy,iz) = sigma[cell].zx();
                sendStress(7,ix,iy,iz) = sigma[cell].zy();
                sendStress(8,ix,iy,iz) = sigma[cell].zz();

            }
        }
    }

}

// Unpacks the 3 components of the velocity-tensor from the socket's
// recvVelocity (cpl::ndArray) storage.
void CPLSocket::unpack
(
    volVectorField &U,
    fvMesh &mesh
)
{


    // Take the mean across spanwise direction if specified by coupler
    // y direction is axis #2
    if (CPL::get<int>("cpl_cfd_bc_slice"))
    {

        Foam::Info << "CPL_CFD_BC_SLICE is on: averaging CFD recvVelocity "
                      "in the x-z plane" << Foam::endl;

        // Number of cells in x-z plane
        int N = recvVelocity.shape(1) * recvVelocity.shape(3);

        // For every component and y-value 
        for (int j = 0; j < recvVelocity.shape(2); ++j)
        {
        for (int c = 0; c < recvVelocity.shape(0); ++c)
        {
            // Sum across the x-z plane 
            double total = 0.0;
            for (int k = 0; k < recvVelocity.shape(3); ++k)
            {
            for (int i = 0; i < recvVelocity.shape(1); ++i)
            {
                total += recvVelocity(c, i, j, k);
            }
            }

            // Find mean by dividing sum by number of cells 
            for (int k = 0; k < recvVelocity.shape(3); ++k)
            {
            for (int i = 0; i < recvVelocity.shape(1); ++i)
            {
                recvVelocity(c, i, j, k) = total / static_cast<double> (N);
            }
            }

        }
        }

    } 

    int applyBCx = CPL::get<int> ("cpl_cfd_bc_x");
    int applyBCy = CPL::get<int> ("cpl_cfd_bc_y");
    int applyBCz = CPL::get<int> ("cpl_cfd_bc_z");

    Foam::string receivePatchName ("CPLReceiveMD");
    Foam::label rvPatchID = mesh.boundary().findPatchID(receivePatchName);
    if (rvPatchID == -1)
    {
        FatalErrorIn
        (
            "CPLSocket::unpack()"
        )
            << " Could not find patch ID " << receivePatchName << ". "
               " Aborting."
            << exit(FatalError);
    }
    Foam::fvPatchVectorField& rvPatch = U.boundaryField()[rvPatchID];
    const Foam::vectorField faceCenters = mesh.boundary()[rvPatchID].Cf();

    //TODO (djt06@ic.ac.uk)
    int myiTmin = iTmin[myCoords[0]] - 1;
    int myjTmin = jTmin[myCoords[1]] - 1;
    int mykTmin = kTmin[myCoords[2]] - 1;

    for (int faceI = 0; faceI != faceCenters.size(); ++faceI)
    {

        double facex = faceCenters[faceI].x();
        double facey = faceCenters[faceI].y();
        double facez = faceCenters[faceI].z();

        // Find the cell indices for this position recvVelocity(:, ix, iy, iz)
        int ix = static_cast<int>(facex / dx) - myiTmin;
        int iy = static_cast<int>(facey / dy) - myjTmin;
        int iz = static_cast<int>(facez / dz) - mykTmin;

        // For now, the coupler passes momentum and mass in velocity array
        // This is due for review. 
        double m = recvVelocity(3, ix, iy, iz); 
        double recvvx = recvVelocity(0, ix, iy, iz) / m;
        double recvvy = recvVelocity(1, ix, iy, iz) / m;
        double recvvz = recvVelocity(2, ix, iy, iz) / m;

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
    }

    if (applyBCx) Foam::Info << "MD->CFD BC x-velocity applied." << Foam::endl;
    if (applyBCy) Foam::Info << "MD->CFD BC y-velocity applied." << Foam::endl;
    if (applyBCz) Foam::Info << "MD->CFD BC z-velocity applied." << Foam::endl;

}

// Sends 9 components of the stress-tensor to overlapping MD processes.
void CPLSocket::send()
{
    const int comm_style = CPL::get<int> ("comm_style"); 
    const int gath_scat = CPL::get<int> ("comm_style_gath_scat"); 
    const int send_recv = CPL::get<int> ("comm_style_send_recv"); 

    std::vector<int> limits (olap_limits);

    if (comm_style == gath_scat)
    {
        // Send stress from CFD to MD processes
        CPL::scatter
        (
            sendStress.data(),
            sendStress.shapeData(),
            limits.data(),
            recvStress.data(),
            recvStress.shapeData()
        );
    }
    else if (comm_style == send_recv)
    {
        FatalErrorIn
        (
            "CPLSocket::send()"
        )
            << " invalid comm_style."
               " Currently only gather/scatter supported by OpenFOAM."
            << exit(FatalError);
    }
    else
    {
        FatalErrorIn
        (
            "CPLSocket::send()"
        )
            << " unrecognised comm_style."
               " Currently only gather/scatter supported by OpenFOAM."
            << exit(FatalError);
    }
  
}

// Receives 3 components of the velocity vector from overlapping MD processes.
void CPLSocket::receive()
{
    // 3-component velocity vector for every local cell +1 for cell mass,
    // CPL::gather requires that the data on the "sending" process is size 0
    int zeroShape[4] = {4, 0, 0, 0};
    int recvShape[4] = {4, ncxl, 1, nczl};

    // Clear and reallocate packed data
    sendVelocity.resize(4, zeroShape);
    recvVelocity.resize(4, recvShape);

    // MD actually sends a slab of cells 1 below the bottom of the CFD domain,
    // so set receive limits region to be 1 cell thick in vertical direction. 
    std::vector<int> limits(olap_limits);
    limits[3] = limits[2];

    const int comm_style = CPL::get<int> ("comm_style"); 
    const int gath_scat = CPL::get<int> ("comm_style_gath_scat"); 
    const int send_recv = CPL::get<int> ("comm_style_send_recv"); 

    if (comm_style == gath_scat)
    {
        // Receive velocity from CFD to MD processes
        CPL::gather
        (
            sendVelocity.data(),
            sendVelocity.shapeData(),
            limits.data(),
            recvVelocity.data(),
            recvVelocity.shapeData()
        );
    }
    else if (comm_style == send_recv)
    {
        FatalErrorIn
        (
            "CPLSocket::receive()"
        )
            << " invalid comm_style."
               " Currently only gather/scatter supported by OpenFOAM."
            << exit(FatalError);
    }
    else 
    {
        FatalErrorIn
        (
            "CPLSocket::receive()"
        )
            << " unrecognised comm_style."
               " Currently only gather/scatter supported by OpenFOAM."
            << exit(FatalError);

    }
  
}
