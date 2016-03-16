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

    "SocketLAMMPS" class for interfacing with CPL-Library.

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

namespace CPL
{

    class SocketLAMMPS
    {

    public:
        
        // Construct from no arguments
        SocketLAMMPS();
        //~SocketLAMMPS();

        // Timesteps and timestep ratio
        int nsteps;
        int timestep_ratio;
        
        // Initialisation routines 
        void initComms (int argc, char **argv);
        void initMD (LAMMPS_NS::LAMMPS *lammps);

        // Data preparation and communication 
        void pack (const LAMMPS_NS::LAMMPS *lammps);
        void unpack (const LAMMPS_NS::LAMMPS *lammps);
        void send();
        void receive();

        // Useful information for main level program
        const MPI_Comm realmCommunicator() {return realmComm;}
        const MPI_Comm cartCommunicator() {return cartComm;}
        const bool isRootProcess() {return (rankRealm == 0);}

        // Clean up MPI/CPL communicators
        void finalizeComms();

    private:

        // Communicators for use with CPL_Library
        MPI_Comm realmComm;
        MPI_Comm cartComm;

        // Rank of this processor in realm and cartComm
        int rankRealm;
        int rankCart;

        // Coordinates of own and all processor(s)
        std::vector<int> myCoords;
        CPL::ndArray<int> allCoords;
        // Also store with Fortran indexing convention (start at 1, not 0)
        std::vector<int> myFCoords;
        CPL::ndArray<int> allFCoords;
        
        // Data to be sent/received with CPL-Library
        CPL::ndArray<double> sendVelocity;
        CPL::ndArray<double> sendStress;
        CPL::ndArray<double> recvVelocity;
        CPL::ndArray<double> recvStress;
    
        // Minimum and maximum CPL cell indices on this processor
        std::vector<int> extents;
        // Minimum and maximum CPL cell indices of global overlap region
        std::vector<int> olap_limits;

        // Number of CPL cells on this processor in each direction
        int ncxl;
        int ncyl;
        int nczl;

        // Internal routines
        void getTimestepRatio();
        void getCellTopology();
        void setupFixMDtoCFD(LAMMPS_NS::LAMMPS *lammps); 
        void setupFixCFDtoMD(LAMMPS_NS::LAMMPS *lammps); 
        void map_cfd_to_lammps(double *coord, LAMMPS_NS::LAMMPS *lammps);

    };

}
#endif // CPL_SOCKET_H_INCLUDED
