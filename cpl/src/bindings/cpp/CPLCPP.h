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

    Declares functions that wrap the C-bindings from cpl/src/bindings/c/CPLC.h
    in the namespace CPL, such that CPL::function_name corresponds to
    CPLC_function_name.

Author(s)

    David Trevelyan

*/

#ifndef CPLCPP_H_INCLUDED
#define CPLCPP_H_INCLUDED
#include <string>
#include "CPLC.h"

// CPL namespace for C++ bindings
namespace CPL
{

    // Realm identifiers, also found in coupler_module.f90
    static const int cfd_realm = 1;
    static const int md_realm = 2;

    void create_comm
    (
        int  calling_realm,
        int& returned_realm_comm
    );

    // All the info (plus some, need to streamline!) from cfd for mapping
    void cfd_init
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
    );

    void md_init
    (
        int& nsteps,
        int& initialstep,
        double dt,
        int icomm_grid,
        int icoord[],
        int npxyz_md[],
        double globaldomain[],
        double density
    );

	void send
	(
			double* asend,	
			int* asend_shape,
			int ndims,
			int icmin, 
			int icmax,
			int jcmin, 
			int jcmax,
			int kcmin, 
			int kcmax,
			bool send_flag
	);

	void recv
	(
			double* arecv,	
			int* arecv_shape,
			int ndims,
			int icmin, 
			int icmax,
			int jcmin, 
			int jcmax,
			int kcmin, 
			int kcmax,
			bool recv_flag
	);

    void scatter
    (
        double* scatterarray,
        int* scatter_shape,
        int* limits,
        double* recvarray,
        int* recv_shape
    );

    void gather
    (
        double* gatherarray,
        int* gather_shape,
        int* limits,
        double* recvarray,
        int* recv_shape
    );

    void proc_extents
    (
        int coord[],
        int realm,
        int extents[]
    );

    void proc_portion
    (
        int coord[],
        int realm,
        int limits[],
        int portion[]
    );


    double* map_cfd2md_global
    (
        double r_cfd[]
    );

    void set_output_mode
    (
        int mode
    );

    // Getters
    template<class T>
    T get(std::string name)
    {
        T (*fp)();
 
        if (name == "ncx") fp = reinterpret_cast<T(*)()> (&CPLC_ncx);
        if (name == "ncy") fp = reinterpret_cast<T(*)()> (&CPLC_ncy);
        if (name == "ncz") fp = reinterpret_cast<T(*)()> (&CPLC_ncz);
        if (name == "icmin_olap") fp = reinterpret_cast<T(*)()> (&CPLC_icmin_olap);
        if (name == "jcmin_olap") fp = reinterpret_cast<T(*)()> (&CPLC_jcmin_olap);
        if (name == "kcmin_olap") fp = reinterpret_cast<T(*)()> (&CPLC_kcmin_olap);
        if (name == "icmax_olap") fp = reinterpret_cast<T(*)()> (&CPLC_icmax_olap);
        if (name == "jcmax_olap") fp = reinterpret_cast<T(*)()> (&CPLC_jcmax_olap);
        if (name == "kcmax_olap") fp = reinterpret_cast<T(*)()> (&CPLC_kcmax_olap);
        if (name == "icmin_cnst") fp = reinterpret_cast<T(*)()> (&CPLC_icmin_cnst);
        if (name == "jcmin_cnst") fp = reinterpret_cast<T(*)()> (&CPLC_jcmin_cnst);
        if (name == "kcmin_cnst") fp = reinterpret_cast<T(*)()> (&CPLC_kcmin_cnst);
        if (name == "icmax_cnst") fp = reinterpret_cast<T(*)()> (&CPLC_icmax_cnst);
        if (name == "jcmax_cnst") fp = reinterpret_cast<T(*)()> (&CPLC_jcmax_cnst);
        if (name == "kcmax_cnst") fp = reinterpret_cast<T(*)()> (&CPLC_kcmax_cnst);
        if (name == "timestep_ratio") fp = reinterpret_cast<T(*)()> (&CPLC_timestep_ratio);
        if (name == "dx") fp = reinterpret_cast<T(*)()> (&CPLC_dx);
        if (name == "dy") fp = reinterpret_cast<T(*)()> (&CPLC_dy);
        if (name == "dz") fp = reinterpret_cast<T(*)()> (&CPLC_dz);
        if (name == "xg") fp = reinterpret_cast<T(*)()> (&CPLC_xg);
        if (name == "yg") fp = reinterpret_cast<T(*)()> (&CPLC_yg);
        if (name == "zg") fp = reinterpret_cast<T(*)()> (&CPLC_zg);
        if (name == "comm_style") fp = reinterpret_cast<T(*)()> (&CPLC_comm_style); 
        if (name == "comm_style_gath_scat") fp = reinterpret_cast<T(*)()> (&CPLC_comm_style_gath_scat); 
        if (name == "comm_style_send_recv") fp = reinterpret_cast<T(*)()> (&CPLC_comm_style_send_recv); 
        if (name == "cpl_cfd_bc_x") fp = reinterpret_cast<T(*)()> (&CPLC_cpl_cfd_bc_x);
        if (name == "cpl_cfd_bc_y") fp = reinterpret_cast<T(*)()> (&CPLC_cpl_cfd_bc_y);
        if (name == "cpl_cfd_bc_z") fp = reinterpret_cast<T(*)()> (&CPLC_cpl_cfd_bc_z);
        if (name == "cpl_md_bc_slice") fp = reinterpret_cast<T(*)()> (&CPLC_cpl_md_bc_slice);
        if (name == "cpl_cfd_bc_slice") fp = reinterpret_cast<T(*)()> (&CPLC_cpl_cfd_bc_slice);

        T got = fp();

        return got;

    }

    double density_cfd();

}

#endif // CPLCPP_H_INCLUDED
