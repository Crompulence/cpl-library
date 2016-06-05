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
#ifndef CPLSocket_H
#define CPLSocket_H

#include<vector>
#include<memory>
#include "mpi.h"
#include "fvCFD.H"
#include "error.H"
#include "CPL.h"

// Nearest integer
int nint(double d)
{
    return std::floor(d+0.5);
}

class CPLSocket
{   
public:
    
    // Construct from no arguments 
    CPLSocket();
    
    // Invoke CPL_Library and store CFD realm communicator
    void initComms(int& argc, char**& argv);

    // Analyse topology and perform CFD-side initialisation
    void initCFD
    (
        const Foam::Time &time, 
        const Foam::fvMesh &mesh
    );

    // Pack and unpack data into CPL-Library formats
    void pack(volVectorField &U, dimensionedScalar &nu, fvMesh &mesh);
    void unpack(volVectorField &U, fvMesh &mesh);

    // Send and receive the data with CPL-Library
    void send();
    void receive();

private:

    // Number of spatial dimensions (fixed to three for CPL_Library)
    static const int ndims=3;

    // Number of processors in x, y and z directions, and periodicity
    int npxyz[ndims]; 
    int nprocs;
    int periods[ndims];

    // Number of cells on a local processor
    int ncxl;
    int ncyl;
    int nczl;

    // Number of cells globally in CFD domain 
    int ncxg;
    int ncyg;
    int nczg;

    // Size of the cells in each direction
    double dx;
    double dy;
    double dz;

    // Density of the fluid assumed by CPL
    double CPLDensity;

    // Cell extents in each direction local to processor
    std::vector<int> iTmax, jTmax, kTmax;
    std::vector<int> iTmin, jTmin, kTmin;
    
    // Communicators for use with CPL_Library 
    MPI_Comm realmComm; 
    MPI_Comm cartComm; 

    // Rank of this processor in cartComm 
    int rankRealm;
    int rankCart;

    // Coordinates of own and all processor(s)
    std::vector<int> myCoords;        
    CPL::ndArray<int> allCoords;        
    void allReduceProcCoords();
    
    // Indices of cells by x, y, z coordinates (3D) (ncxl, ncyl, nczl)
    CPL::ndArray<int> cellIndex;

    // Data to be sent/received with CPL-Library
    CPL::ndArray<double> sendVelocity;
    CPL::ndArray<double> sendStress;
    CPL::ndArray<double> recvVelocity;
    CPL::ndArray<double> recvStress;

    // Overlap limits
    std::vector<int> olap_limits;

    // Some utility functions
    Foam::Vector<double> maxPoints (Foam::pointField points)
    {
        double xmax=-std::numeric_limits<double>::max();
        double ymax=-std::numeric_limits<double>::max();
        double zmax=-std::numeric_limits<double>::max();
        double x, y, z;
        forAll(points, i)
        {
            x = points[i].x();
            y = points[i].y();
            z = points[i].z();
            if (x > xmax) xmax = x;
            if (y > ymax) ymax = y;
            if (z > zmax) zmax = z;
        }

        Foam::Vector<double> coords(xmax, ymax, zmax);

        return coords;
    }

    Foam::Vector<double> minPoints (Foam::pointField points)
    {
        double xmin=std::numeric_limits<double>::max();
        double ymin=std::numeric_limits<double>::max();
        double zmin=std::numeric_limits<double>::max();
        double x, y, z;
        forAll(points, i)
        {
            x = points[i].x();
            y = points[i].y();
            z = points[i].z();
            if (x < xmin) xmin = x;
            if (y < ymin) ymin = y;
            if (z < zmin) zmin = z;
        }

        Foam::Vector<double> coords (xmin, ymin, zmin);

        return coords;
    }

};


#endif
