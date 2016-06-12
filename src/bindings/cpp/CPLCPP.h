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

    void init ( int  calling_realm, int& returned_realm_comm); 
    void finalize(); 

    // All the info (plus some, need to streamline!) from cfd for mapping
    void setup_cfd
    (
        int icomm_grid,
        double xyzL[],
        double xyz_orig[],
        int ncxyz[]
    );

    void setup_md
    (
        int icomm_grid,
        double xyzL[],
        double xyz_orig[]
    );


    bool send(double* asend, int* asend_shape, int* limits);
    bool recv(double* arecv, int* arecv_shape, int* limits);

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

    void my_proc_extents
    (
        int extents[]
    );

    void proc_portion
    (
        int coord[],
        int realm,
        int limits[],
        int portion[]
    );


    void my_proc_portion
    (
        int limits[],
        int portion[]
    );

    void map_cell2coord
    (
        int i,
        int j,
        int k,
        double coord_xyz[]
    );

    bool map_coord2cell
    (
        double x,
        double y,
        double z,
        int cell_ijk[]
    );

    void get_no_cells
    (
        int limits[],
        int no_cells[]
    );

    bool map_glob2loc_cell
    (
        int limits[],
        int glob_cell[],
        int loc_cell[]
    );

    void get_olap_limits
    (
        int limits[]
    );

    void get_cnst_limits
    (
        int limits[]
    );

    bool map_cfd2md_coord
    (
        double cfd_coord[],
        double md_coord[]
    );

    bool map_md2cfd_coord
    (
        double md_coord[],
        double cfd_coord[]
    );

    bool overlap();

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
        if (name == "comm_style") fp = reinterpret_cast<T(*)()> (&CPLC_comm_style); 
        if (name == "comm_style_gath_scat") fp = reinterpret_cast<T(*)()> (&CPLC_comm_style_gath_scat); 
        if (name == "comm_style_send_recv") fp = reinterpret_cast<T(*)()> (&CPLC_comm_style_send_recv); 
        if (name == "cpl_cfd_bc_x") fp = reinterpret_cast<T(*)()> (&CPLC_cpl_cfd_bc_x);
        if (name == "cpl_cfd_bc_y") fp = reinterpret_cast<T(*)()> (&CPLC_cpl_cfd_bc_y);
        if (name == "cpl_cfd_bc_z") fp = reinterpret_cast<T(*)()> (&CPLC_cpl_cfd_bc_z);
        if (name == "cpl_md_bc_slice") fp = reinterpret_cast<T(*)()> (&CPLC_cpl_md_bc_slice);
        if (name == "cpl_cfd_bc_slice") fp = reinterpret_cast<T(*)()> (&CPLC_cpl_cfd_bc_slice);
        if (name == "cpl_x_orig_cfd") fp = reinterpret_cast<T(*)()> (&CPLC_x_orig_cfd);
        if (name == "cpl_y_orig_cfd") fp = reinterpret_cast<T(*)()> (&CPLC_y_orig_cfd);
        if (name == "cpl_z_orig_cfd") fp = reinterpret_cast<T(*)()> (&CPLC_z_orig_cfd);
        if (name == "cpl_x_orig_md") fp = reinterpret_cast<T(*)()> (&CPLC_x_orig_md);
        if (name == "cpl_y_orig_md") fp = reinterpret_cast<T(*)()> (&CPLC_y_orig_md);
        if (name == "cpl_z_orig_md") fp = reinterpret_cast<T(*)()> (&CPLC_z_orig_md);

        T got = fp();

        return got;

    }

    double density_cfd();

}

#endif // CPLCPP_H_INCLUDED
