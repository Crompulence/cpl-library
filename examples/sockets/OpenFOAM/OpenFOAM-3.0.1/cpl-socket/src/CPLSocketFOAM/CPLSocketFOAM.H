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

    "Socket" class for interfacing with CPL-Library. 

*/
#ifndef CPLSocketFOAM_H
#define CPLSocketFOAM_H

#include<vector>
#include<memory>
#include "mpi.h"
#include "fvCFD.H"
#include "error.H"
#include "meshSearch.H"

#include "cpl/CPL_ndArray.h"
#include "cpl/CPL.h"

// Nearest integer
int nint(double d)
{
    return std::floor(d+0.5);
}

class CPLSocketFOAM
{   

public:
    
    // Construct from no arguments 
    CPLSocketFOAM() : olapRegion(6), velBCRegion(6), cnstFRegion(6),
                   olapPortion(6), velBCPortion(6), cnstFPortion(6), myCoords(3) {}
    
    // Invoke CPL_Library and store CFD realm communicator
    void initComms(int& argc, char**& argv);

    // Analyse topology and perform CFD-side initialisation
    void initCFD (const Foam::Time &runTime, const Foam::fvMesh &mesh);

    void getCellTopology();
    void allocateBuffers(int);

    // Pack and unpack data into CPL-Library formats
    void packStress(volVectorField &U, dimensionedScalar &nu, fvMesh &mesh);
    double unpackVelocity(volVectorField &U, fvMesh &mesh);
    void pack(volVectorField &U, dimensionedScalar &nu, fvMesh &mesh, int);
    double unpackPorousForce(volVectorField &F, volScalarField &eps, fvMesh &mesh);
    volVectorField divideFieldsVectorbyScalar(volVectorField &F, volScalarField &eps, fvMesh &mesh);

    // Send and receive the data with CPL-Library
    void sendStress();
    void recvVelocity();
    void send();
    void recv();

    int const DEBUG = -1;
    int const PACKVELONLY = 0;
    int const PACKSTRESSONLY = 1;
    int const PACKVELSTRESS = 2;

	~CPLSocketFOAM() {
		delete meshSearcher;
	}

private:


    // Number of processors in x, y and z directions, and periodicity
    int npxyz[3]; 
    int nprocs;
    int periods[3];

    // Domain dimensions
    double xyzL[3];

    // Number of cells globally in CFD domain 
    int ncxyz[3];

    // Size of the cells in each direction
    double dx, dy, dz;

    // Density of the fluid assumed by CPL
    double CPLDensity;

    // Communicators for use with CPL_Library 
    MPI_Comm realmComm; 
    MPI_Comm cartComm; 

    // Rank of this processor in cartComm 
    int rankRealm;
    int rankCart;

    // Coordinates of the local processor
    std::vector<int> myCoords;

    // Data to be sent/received with CPL-Library
    CPL::ndArray<double> sendBuf;
    CPL::ndArray<double> recvBuf;
    CPL::ndArray<double> sendStressBuff;
    CPL::ndArray<double> recvVelocityBuff;

    // Region cell limits
    std::vector<int> olapRegion;
    std::vector<int> velBCRegion;
    std::vector<int> cnstFRegion;

    // Portions of the region in the local processor
    std::vector<int> olapPortion;
    std::vector<int> velBCPortion;
    std::vector<int> cnstFPortion;
    
    // Number of cells in the portion in a region
    int olapCells[3];
    int velBCCells[3];
    int cnstFCells[3];

	//meshSearch local search for points in a mesh
	Foam::meshSearch* meshSearcher;

};


#endif
