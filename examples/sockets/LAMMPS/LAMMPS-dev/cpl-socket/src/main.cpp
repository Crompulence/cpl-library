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

    Simple program that couples LAMMPS to an external solver via
    CPL-Library. This program takes a single command-line argument,
    which is the path to a setup script (arg[1]) that is performed
    line-by-line before the coupled run.

Author(s)

    David Trevelyan

*/
#include <memory>
#include <fstream>
#include <string>
#include <iostream>
#include "mpi.h"

// LAMMPS
#include "lammps.h"
#include "input.h"

// CPL
#include "CPL.h"
#include "socket.h"

int main(int narg, char **arg)
{
    // Create socket and initialise MPI/CPL communicators
    SocketLAMMPS cpl;
    cpl.initComms (narg, arg);

    // Create lammps instance on socket realm communicator
//    auto lmp = std::make_unique <LAMMPS_NS::LAMMPS> (1, arg, 
//                                                     cpl.realmCommunicator());
    LAMMPS_NS::LAMMPS *lmp = new LAMMPS_NS::LAMMPS(1, arg, cpl.realmCommunicator());

    // Open LAMMPS input script and perform setup line-by-line
    std::ifstream inputFile (arg[1]);
    if (!inputFile.is_open())
    {
        std::cerr << "ERROR: Could not open LAMMPS setup script "
                  << arg[1] << std::endl;
        MPI_Abort (MPI_COMM_WORLD, 1);
    }
    std::string line;
    while (!inputFile.eof())
    {
        std::getline (inputFile, line);
        lmp->input->one (line.c_str());
    }

    cpl.initMD (lmp);

    for (int step = 0; step < cpl.nsteps; step += cpl.timestep_ratio)
    {
        // Communications
        cpl.recvStress();
        cpl.unpackStress(lmp);
        cpl.packVelocity(lmp);
        cpl.sendVelocity();

        // Continue LAMMPS 
        line = "run " + std::to_string(cpl.timestep_ratio);
        lmp->input->one (line.c_str());

    }

    // Finalize MPI
    delete lmp;
    cpl.finalizeComms();

}
