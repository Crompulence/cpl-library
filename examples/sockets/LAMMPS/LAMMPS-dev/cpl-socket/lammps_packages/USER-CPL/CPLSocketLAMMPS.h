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

    "CPLSocketLAMMPS" class for interfacing with CPL-Library.

Author(s)

    David Trevelyan

*/
#ifndef CPL_SOCKET_H_INCLUDED
#define CPL_SOCKET_H_INCLUDED

#include<vector>
#include<memory>
#include "mpi.h"
#include "lammps.h"
#include "CPL.h"
#include "fix_cpl_force.h"

#include "fix_ave_chunk.h"
#include "region.h"
typedef CPL::ndArray<double> arrayDoub;
typedef std::shared_ptr<arrayDoub> shaPtrArrayDoub;


class CPLSocketLAMMPS
{

public:
    
    // Construct from no arguments
    CPLSocketLAMMPS() : myCoords(3), olapRegion(6), velBCRegion(6), cnstFRegion(6),
                     velBCPortion(6), cnstFPortion(6) {}
    //~CPLSocketLAMMPS();

    // Timesteps and timestep ratio
    int nsteps;
    int timestep_ratio;
    
    // Initialisation routines 
    void initComms ();
    void initMD (LAMMPS_NS::LAMMPS *lammps);

    // Data preparation and communication 
    void packVelocity(const LAMMPS_NS::LAMMPS *lammps);
    void sendVelocity();
    void unpackStress(const LAMMPS_NS::LAMMPS *lammps);
    void recvStress();

    // Useful information for main level program
    const MPI_Comm realmCommunicator() {return realmComm;}
    const MPI_Comm cartCommunicator() {return cartComm;}
    const bool isRootProcess() {return (rankRealm == 0);}

    // Clean up MPI/CPL communicators
    void finalizeComms();

    void setupFixMDtoCFD(LAMMPS_NS::LAMMPS *lammps); 
    void setupFixCFDtoMD(LAMMPS_NS::LAMMPS *lammps); 


private:
    
    double VELBC_BELOW = 0.0;
    double VELBC_ABOVE = 0.0;

    // Cartesian coordinates of the processor
    std::vector<int> myCoords;

    // Communicators for use with CPL_Library
    MPI_Comm realmComm;
    MPI_Comm cartComm;

    // Rank of this processor in realm and cartComm
    int rankRealm;
    int rankCart;

    // Communication regions in the overlap
    std::vector<int> olapRegion;
    std::vector<int> velBCRegion;
    std::vector<int> cnstFRegion;

    // Portions of the regions in the local processor
    std::vector<int> velBCPortion;
    std::vector<int> cnstFPortion;

    // Number of cells in the portion for each region
    int velBCCells[3];
    int cnstFCells[3];

    // Data to be sent/received with CPL-Library
    arrayDoub sendVelocityBuff;
    arrayDoub sendStressBuff;
    arrayDoub recvVelocityBuff;
    arrayDoub recvStressBuff;


    // Cell sizes
    double dx, dy, dz;
    //Appropriate region, compute and fix    
    class LAMMPS_NS::Region *cfdbcregion, *cplforceregion;
    class LAMMPS_NS::Compute *cfdbccompute;
    class LAMMPS_NS::Fix *cfdbcfix, *cplforcefix;
    class LAMMPS_NS::Group *cplforcegroup;
    
    // Fix that applies the momentum constrain
    FixCPLForce* cplfix;

    // Internal grid
    arrayDoub cfd_xg; 
    arrayDoub cfd_yg;
    arrayDoub cfd_zg;

    // Internal routines
    void getCellTopology();
    void allocateBuffers();
};

#endif // CPL_SOCKET_H_INCLUDED
