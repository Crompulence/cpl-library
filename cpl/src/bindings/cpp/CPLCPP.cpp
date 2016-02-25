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

    Defines functions declared in namespace CPL, that wrap the C language
    bindings provided by cpl/src/bindings/c/CPLC.[f90&&h]

*/
#include "CPLCPP.h"

void CPL::create_comm
(
    int  calling_realm,
    int& returned_realm_comm
)
{
    CPLC_create_comm
    (
    		calling_realm,
     		&returned_realm_comm
       //returned_realm_comm

    );
}

void CPL::cfd_init
(
    int nsteps,
    double dt,
    int icomm_grid,
    int icoord[],
    int npxyz_cfd[],
    double xyzL[],
    int ncxyz[],
    double density,
    int ijkcmax[],
    int ijkcmin[],
    int iTmin[],
    int iTmax[],
    int jTmin[],
    int jTmax[],
    int kTmin[],
    int kTmax[],
    double xgrid[],
    double ygrid[],
    double zgrid[]
)
{
    CPLC_cfd_init
    (
        nsteps,
        dt,
        icomm_grid,
        icoord,
        npxyz_cfd,
        xyzL,
        ncxyz,
        density,
        ijkcmax,
        ijkcmin,
        iTmin,
        iTmax,
        jTmin,
        jTmax,
        kTmin,
        kTmax,
        xgrid,
        ygrid,
        zgrid
    );
}

void CPL::md_init
(
    int& nsteps,
    int& initialstep,
    double dt,
    int icomm_grid,
    int icoord[],
    int npxyz_md[],
    double globaldomain[],
    double density
)
{
    CPLC_md_init
    (
        &nsteps,
        &initialstep,
        dt,
        icomm_grid,
        icoord,
        npxyz_md,
        globaldomain,
        density
    );
}


void CPL::send
(
		double* asend,	
		double* asend_shape,
		int ndims,
		int icmin, 
		int icmax,
		int jcmin, 
		int jcmax,
		int kcmin, 
		int kcmax,
		bool& send_flag

)
{
		CPLC_send
		(
				asend,	
				asend_shape,
				ndims,
				icmin, 
				icmax,
				jcmin, 
				jcmax,
				kcmin, 
				kcmax,
				&send_flag
		);
}


void CPL::recv
(
		double* arecv,	
		double* arecv_shape,
		int ndims,
		int icmin, 
		int icmax,
		int jcmin, 
		int jcmax,
		int kcmin, 
		int kcmax,
		bool& recv_flag

)
{
		CPLC_recv
		(
				arecv,	
				arecv_shape,
				ndims,
				icmin, 
				icmax,
				jcmin, 
				jcmax,
				kcmin, 
				kcmax,
				&recv_flag
		);
}

void CPL::scatter
(
    double* scatterarray,
    int* scatter_shape,
    int* limits,
    double* recvarray,
    int* recv_shape
)
{
    CPLC_scatter
    (
        scatterarray,
        scatter_shape,
        limits,
        recvarray,
        recv_shape
    );
}

void CPL::gather
(
    double* gatherarray,
    int* gather_shape,
    int* limits,
    double* recvarray,
    int* recv_shape
)
{
    CPLC_gather
    (
        gatherarray,
        gather_shape,
        limits,
        recvarray,
        recv_shape
    );
}

void CPL::proc_extents
(
    int coord[],
    int realm,
    int extents[]
)
{
    CPLC_proc_extents
    (
        coord,
        realm,
        extents
    );
}

void CPL::proc_portion
(
    int coord[],
    int realm,
    int limits[],
    int portion[]
)
{
    CPLC_proc_portion
    (
        coord,
        realm,
        limits,
        portion
    );
}

double* CPL::map_cfd2md_global
(
    double r_cfd[]
)
{
    return CPLC_map_cfd2md_global (r_cfd);
}

// Setters
void CPL::set_output_mode
(
    int mode
)
{
    CPLC_set_output_mode(mode);
}

// Getters
double CPL::density_cfd()
{
    return CPLC_density_cfd();
}
