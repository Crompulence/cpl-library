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

    Function that creates a new MPI Cartesian communicator, where the
    coordinates of each process are specified by the user. All inputs
    to the function CPL::Cart_create are the same as required by
    MPI_Cart_create, except for the addition of the user-specified coordinates
    for this process "coords[]".

Author(s)
    
    David Trevelyan

*/

#ifndef CPL_CART_CREATE_H_INCLUDED
#define CPL_CART_CREATE_H_INCLUDED

#include<vector>
#include "mpi.h"

namespace CPL
{

    void Cart_create
    (
        const MPI_Comm oldComm, // same as MPI_Cart_create
        const int ndims,        // same as MPI_Cart_create
        const int dims[],       // same as MPI_Cart_create
        const int periods[],    // same as MPI_Cart_create
        const int coords[],     // desired input coordinates for this process
        MPI_Comm *newCartComm   // returned MPI Cart comm with desired coords
    );

}

#endif // CPL_CART_CREATE_H_INCLUDED
