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

    Declares the functions defined in cpl/src/bindings/CPLC.f90 for inclusion
    in C language programs.
    
Author(s)

    David Trevelyan

*/

extern "C" void CPLC_create_comm
(
    int  calling_realm,
    int* returned_realm_comm
);

extern "C" void CPLC_cfd_init
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

extern "C" void CPLC_test_python
(
 		int int_p, 
		double doub_p, 
		bool bool_p, 
		int* int_pptr,
		double* doub_pptr, 
		int* int_pptr_dims,
		int* doub_pptr_dims
);

extern "C" void CPLC_md_init
(
    int* nsteps,
    int* initialstep,
    double dt,
    int icomm_grid,
    int icoord[],
    int npxyz_md[],
    double globaldomain[],
		double density
);

extern "C" void CPLC_send
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
		bool* send_flag

);

extern "C" void CPLC_recv
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
		bool* recv_flag

);



extern "C" void CPLC_scatter
(
    double* scatterarray,
    int* scatter_shape,
    int* limits,
    double* recvarray,
    int* recv_shape
);

extern "C" void CPLC_gather
(
    double* gatherarray,
    int* gather_shape,
    int* limits,
    double* recvarray,
    int* recv_shape
);

extern "C" void CPLC_proc_extents
(
    int coord[],
    int realm,
    int extents[]
);

extern "C" void CPLC_proc_portion
(
    int coord[],
    int realm,
    int limits[],
    int portion[]
);

extern "C" double* CPLC_map_cfd2md_global
(
    double r_cfd[]
);

// Setters
extern "C" void CPLC_set_output_mode
(
    int mode
);

// Getters
extern "C" int CPLC_icmin_olap();
extern "C" int CPLC_jcmin_olap();
extern "C" int CPLC_kcmin_olap();
extern "C" int CPLC_icmax_olap();
extern "C" int CPLC_jcmax_olap();
extern "C" int CPLC_kcmax_olap();
extern "C" int CPLC_icmin_cnst();
extern "C" int CPLC_jcmin_cnst();
extern "C" int CPLC_kcmin_cnst();
extern "C" int CPLC_icmax_cnst();
extern "C" int CPLC_jcmax_cnst();
extern "C" int CPLC_kcmax_cnst();
extern "C" int CPLC_timestep_ratio();
extern "C" int CPLC_ncx();
extern "C" int CPLC_ncy();
extern "C" int CPLC_ncz();
extern "C" int CPLC_npx_md();
extern "C" int CPLC_npy_md();
extern "C" int CPLC_npz_md();
extern "C" int CPLC_npx_cfd();
extern "C" int CPLC_npy_cfd();
extern "C" int CPLC_npz_cfd();
extern "C" int CPLC_comm_style();
extern "C" int CPLC_comm_style_send_recv();
extern "C" int CPLC_comm_style_gath_scat();
extern "C" int CPLC_cpl_cfd_bc_x();
extern "C" int CPLC_cpl_cfd_bc_y();
extern "C" int CPLC_cpl_cfd_bc_z();
extern "C" int CPLC_cpl_md_bc_slice();
extern "C" int CPLC_cpl_cfd_bc_slice();
extern "C" double CPLC_density_cfd();
extern "C" double CPLC_dx();
extern "C" double CPLC_dy();
extern "C" double CPLC_dz();
extern "C" double CPLC_xl_md();
extern "C" double CPLC_yl_md();
extern "C" double CPLC_zl_md();
extern "C" double CPLC_xl_cfd();
extern "C" double CPLC_yl_cfd();
extern "C" double CPLC_zl_cfd();
extern "C" double* CPLC_xg();
extern "C" double* CPLC_yg();
extern "C" double* CPLC_zg();
