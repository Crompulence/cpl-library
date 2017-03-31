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

    See CPLSocketFOAM.H 

*/
#include "CPLSocketFOAM.H"
#include "blockMesh.H"
#include <sstream>
#include <unistd.h>


// Initialise CFD realm communicator
void CPLSocketFOAM::initComms (int& argc, char**& argv) {
    CPL::init (CPL::cfd_realm, realmComm);
    MPI_Comm_rank (realmComm, &rankRealm);

}

// Analyse mesh topology and perform CFD-side CPL_init.
void CPLSocketFOAM::initCFD (const Foam::Time &runTime, const Foam::fvMesh &mesh) {

	meshSearcher = new Foam::meshSearch(mesh);

    Foam::Info << "CPLSocketFOAM: Analysing processor and mesh topology"
               << Foam::endl;

    // Read from decomposePar dictionary the number of processors in each
    // direction. Must be decomposed with the "simple" method.
	//TODO: THis should not be necessary if the cell grid is defined in the config
	//file for the overlap region
    Foam::IOdictionary decomposeDict(Foam::IOobject ("decomposeParDict", 
									 runTime.time().system(), runTime,
                        			 IOobject::MUST_READ, IOobject::NO_WRITE,
									 false));
    Foam::dictionary simpleCoeffs = decomposeDict.subDict ("simpleCoeffs");
    Foam::Vector<int> np = simpleCoeffs.lookup ("n");
    nprocs = np.x() * np.y() * np.z();

    // Define arrays needed by MPI & CPL cart create routines 
    npxyz[0] = np.x();
    npxyz[1] = np.y();
    npxyz[2] = np.z();
    periods[0] = 0;
    periods[1] = 0;
    periods[2] = 0;

    Foam::Info << "CPLSocketFOAM: Defining new MPI Cartesian communicator"
               << Foam::endl;

    // Create custom cartesian communicator (cartComm) based on myCoords	
	// Assume xyz ordering (Simple decomposition)
	myCoords[2] = Pstream::myProcNo() / (np.x()*np.y());
	int mod_aux = Pstream::myProcNo() % (np.x()*np.y());
	myCoords[1] = mod_aux / np.x();
	myCoords[0] = mod_aux % np.x();
	
    CPL::Cart_create (realmComm, 3, npxyz, periods, 
                      myCoords.data(), &cartComm);

    MPI_Comm_rank (cartComm, &rankCart);


    // Prepare inputs for CPL::cfd_init 
    double dt_cfd = runTime.deltaTValue();
    int nsteps = nint ((runTime.endTime().value() - \
                        runTime.startTime().value()) / dt_cfd);

    Foam::IOdictionary blockMeshDict(Foam::IOobject ("polyMesh/blockMeshDict", 
									 runTime.time().constant(), runTime,
                        			 IOobject::MUST_READ, 
									 IOobject::NO_WRITE, false));

	Foam::List<Foam::Vector<double>> vertices(blockMeshDict.lookup("vertices"));

    // Domain dimensions
    xyzL[0] = vertices[1][0] - vertices[0][0];
    xyzL[1] = vertices[3][1] - vertices[0][1];
    xyzL[2] = vertices[4][2] - vertices[0][2];
  
    double dummyDensity = -666.0;


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
	CPL::set_timing(0, nsteps, dt_cfd);
    CPL::setup_cfd (cartComm, xyzL, xyz_orig, ncxyz);


    getCellTopology();
    allocateBuffers();

    // Store some values from CPL that are useful later
    CPLDensity = CPL::density_cfd();

    Foam::Info << "OpenFOAM CPL topology initialisation complete" << Foam::endl;

}

void CPLSocketFOAM::getCellTopology() {                                                                                                                                                                                                                       
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

void CPLSocketFOAM::allocateBuffers() {
    // Received stress field
    int sendShape[4] = {9, cnstFCells[0], cnstFCells[1], cnstFCells[2]};
    sendStressBuff.resize(4, sendShape);

    // LAMMPS computed velocity field
    int recvShape[4] = {4, velBCCells[0], velBCCells[1], velBCCells[2]};
    recvVelocityBuff.resize (4, recvShape);
}  
    
// Packs the 9 components of the stress-tensor to the socket's CPL::ndArray
// storage.
void CPLSocketFOAM::packStress(volVectorField &U, dimensionedScalar &nu, fvMesh &mesh)
{
    Foam::dimensionedScalar mu(CPLDensity*nu);
	Foam::volSymmTensorField sigma(nu*2*dev(symm(fvc::grad(U))));
	if (CPL::is_proc_inside(cnstFPortion.data())) {

		// Loop over socket cells
		Foam::label cell;
		Foam::point globalPos;
		int glob_cell[3], loc_cell[3];
		for (int ix = cnstFPortion[0]; ix <= cnstFPortion[1]; ix++) {
			for (int iy = cnstFPortion[2]; iy <= cnstFPortion[3]; iy++) {
				for (int iz = cnstFPortion[4]; iz <= cnstFPortion[5]; iz++) {
					// Global position at cell center
					glob_cell[0] = ix;
					glob_cell[1] = iy;
					glob_cell[2] = iz;
					CPL::map_glob2loc_cell(cnstFPortion.data(), glob_cell, loc_cell);
					globalPos = Foam::point((glob_cell[0] + 0.5) * dx,
											(glob_cell[1] + 0.5) * dy, 
											(glob_cell[2] + 0.5) * dz);
					cell = meshSearcher->findNearestCell(globalPos);
					if (cell != -1) {
						sendStressBuff(0,loc_cell[0],loc_cell[1],loc_cell[2]) = sigma[cell].xx();
						sendStressBuff(1,loc_cell[0],loc_cell[1],loc_cell[2]) = sigma[cell].xy();
						sendStressBuff(2,loc_cell[0],loc_cell[1],loc_cell[2]) = sigma[cell].xz();
						sendStressBuff(3,loc_cell[0],loc_cell[1],loc_cell[2]) = sigma[cell].xy();
						sendStressBuff(4,loc_cell[0],loc_cell[1],loc_cell[2]) = sigma[cell].yy();
						sendStressBuff(5,loc_cell[0],loc_cell[1],loc_cell[2]) = sigma[cell].yz();
						sendStressBuff(6,loc_cell[0],loc_cell[1],loc_cell[2]) = sigma[cell].xz();
						sendStressBuff(7,loc_cell[0],loc_cell[1],loc_cell[2]) = sigma[cell].yz();
						sendStressBuff(8,loc_cell[0],loc_cell[1],loc_cell[2]) = sigma[cell].zz();
					}
					else
						std::cerr << "Warning: The point (" << (glob_cell[0]+0.5)*dx << "," << 
									(glob_cell[1]+0.5)*dy << "," << (glob_cell[2]+0.5)*dz << 
									") is outside the mesh for whatever reason! - Cell: " << 
									cell << std::endl;

				}
			}
		}
	}
}

// Unpacks the 3 components of the velocity-tensor from the socket's
// recvVelocity (cpl::ndArray) storage.
double CPLSocketFOAM::\
unpackVelocity(volVectorField &U, fvMesh &mesh) {

	if (CPL::is_proc_inside(velBCPortion.data())) { 

		// TODO: Make this a utility general function that can be used on buffers
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

		// Apply BCs only in certain directions.
		int applyBCx = CPL::get<int> ("cpl_cfd_bc_x");
		int applyBCy = CPL::get<int> ("cpl_cfd_bc_y");
		int applyBCz = CPL::get<int> ("cpl_cfd_bc_z");

		// Patch receiving B.Cs
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

		Foam::label cell;
		Foam::point closestCellCentre;
		for (int faceI = 0; faceI != faceCenters.size(); ++faceI) {
			double facex = faceCenters[faceI].x();
			double facey = faceCenters[faceI].y();
			double facez = faceCenters[faceI].z();

			// Find the cell indices for this position recvVelocity(:, ix, iy, iz)
			int glob_cell[3];
			CPL::map_coord2cell(facex, facey, facez, glob_cell);
			int loc_cell[3];
			CPL::map_glob2loc_cell(velBCPortion.data(), glob_cell, loc_cell);

			//double m = recvVelocityBuff(3, loc_cell[0], loc_cell[1], loc_cell[2]); 
			double m = 1.0;
			double recvvx = recvVelocityBuff(0, loc_cell[0], loc_cell[1], loc_cell[2])/m;
			double recvvy = recvVelocityBuff(1, loc_cell[0], loc_cell[1], loc_cell[2]) /m;
			double recvvz = recvVelocityBuff(2, loc_cell[0], loc_cell[1], loc_cell[2])/m;

			if (applyBCx) rvPatch[faceI].x() = recvvx;
			if (applyBCy) rvPatch[faceI].y() = recvvy;
			if (applyBCz) rvPatch[faceI].z() = recvvz;
		}

		if (applyBCx) Foam::Info << "MD->CFD BC x-velocity applied." << Foam::endl;
		if (applyBCy) Foam::Info << "MD->CFD BC y-velocity applied." << Foam::endl;
		if (applyBCz) Foam::Info << "MD->CFD BC z-velocity applied." << Foam::endl;
	}
}

// Sends 9 components of the stress-tensor to overlapping MD processes.
void CPLSocketFOAM::sendStress() {
        CPL::send (sendStressBuff.data(), sendStressBuff.shapeData(),
                   cnstFRegion.data());
}

// Receives 3 components of the velocity vector from overlapping MD processes.
void CPLSocketFOAM::recvVelocity() {
    CPL::recv(recvVelocityBuff.data(), recvVelocityBuff.shapeData(), 
			  velBCRegion.data());
}
  
