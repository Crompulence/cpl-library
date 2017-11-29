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

    See CPL_cartCreate.h

Author(s)
    
    David Trevelyan

*/
#include<vector>
#include<iostream>
#include "mpi.h"
#include "CPL_cartCreate.h"

void CPL::Cart_create
(

    const MPI_Comm oldComm,           // same as MPI_Cart_create
    int ndims,                        // same as MPI_Cart_create
    const CPL::IntVector& dims,       // same as MPI_Cart_create
    const CPL::IntVector& periods,    // same as MPI_Cart_create
    const CPL::IntVector& coords,    // desired input coordinates for this process
    MPI_Comm* newCartComm             // returned MPI Cart comm with desired coords
)
{

    MPI_Comm dummyCartComm, tempComm;
    int oldCommRank, desiredRank, newCartCommRank, newCoordsArray[3];
    CPL::IntVector newCoords(3);

    // Save old communicator rank for debugging purposes
    MPI_Comm_rank(oldComm, &oldCommRank);

    // New dummy cart comm, save what the rank of desired coords would be
    MPI_Cart_create(oldComm, ndims, &dims[0], &periods[0], false, &dummyCartComm);
    MPI_Cart_rank(dummyCartComm, &coords[0], &desiredRank);
 
    // Force this rank in new comm with comm_split trick
    MPI_Comm_split(oldComm, 0, desiredRank, &tempComm);

    // Create cartesian topology based on desied ranks
    MPI_Cart_create(tempComm, ndims, &dims[0], &periods[0], false, newCartComm);

    // Check coordinates match desired values
    MPI_Comm_rank(*newCartComm, &newCartCommRank);
    MPI_Cart_coords(*newCartComm, newCartCommRank, 3, &newCoords[0]);
    
    if ((coords != newCoords).max())
    {
        std::cerr
            << "New coordinates in cartesian communicator "
            << " do not match the desired coordinates for rank "
            << oldCommRank
            << " on communicator "
            << oldComm
            << std::endl;
        throw;
    }

}
