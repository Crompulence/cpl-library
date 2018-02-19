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

    David Trevelyan, Edward Smith, Eduardo Ramos Fernandez

*/

/** MODULE IO **/

#ifdef JSON_SUPPORT
extern "C" void CPLC_load_param_file(const char* fname);
extern "C" void CPLC_close_param_file();
extern "C" void CPLC_get_real_param(const char* section, const char* param_name, double* real_param);
extern "C" void CPLC_get_real_array_param(const char* section, const char* param_name, double** real_param_array, int* array_len);
extern "C" void CPLC_get_int_param(const char* section, const char* param_name, int* int_param);
extern "C" void CPLC_get_int_array_param(const char* section, const char* param_name, int** int_param_array, int* array_len);
extern "C" void CPLC_get_boolean_param(const char* section, const char* param_name, bool* boolean_param);
extern "C" void CPLC_get_boolean_array_param(const char* section, const char* param_name, bool** boolean_param_array, int* array_len);
extern "C" void CPLC_get_string_param(const char* section, const char* param_name, char** string_param);
extern "C" void CPLC_get_string_array_param(const char* section, const char* param_name, char** string_param_array, int* array_len);
#endif 

extern "C" void CPLC_init(int  calling_realm, int* returned_realm_comm);
extern "C" void CPLC_finalize();
extern "C" void CPLC_setup_cfd(int icomm_grid, double xyzL[], 
                               double xyz_orig[], int ncxyz[]);

extern "C" void CPLC_test_python(int int_p, double doub_p, bool bool_p, 
                                 int* int_pptr, double* doub_pptr, 
                                 int* int_pptr_dims, int* doub_pptr_dims);

extern "C" void CPLC_setup_md(int icomm_grid, double xyzL[], 
                              double xyz_orig[]);

extern "C" void CPLC_set_timing(int initialstep, int nsteps, double dt);
extern "C" void CPLC_send(double* asend, int* asend_shape, int* limits, bool *send_flag);
extern "C" void CPLC_recv(double* arecv, int* arecv_shape, int* limits, bool *recv_flag);

extern "C" void CPLC_scatter(double* scatterarray, int* scatter_shape, int* limits, 
                             double* recvarray, int* recv_shape);

extern "C" void CPLC_gather(double* gatherarray, int* gather_shape, int* limits, 
                            double* recvarray, int* recv_shape);

extern "C" void CPLC_swaphalos(double* A, int* A_shape);

extern "C" void CPLC_proc_extents(int coord[], int realm, int extents[]);
extern "C" void CPLC_my_proc_extents(int extents[]);
extern "C" void CPLC_proc_portion(int coord[], int realm, int limits[], int portion[]);
extern "C" void CPLC_my_proc_portion(int limits[], int portion[]);
extern "C" void CPLC_map_cell2coord(int i, int j, int k, double coord_xyz[]);
extern "C" bool CPLC_map_coord2cell(double x,  double y,  double z,  int cell_ijk[]);
extern "C" void CPLC_get_no_cells(int limits[], int no_cells[]);
extern "C" bool CPLC_map_glob2loc_cell(int limits[], int glob_cell[], int loc_cell[]);
extern "C" void CPLC_get_olap_limits(int limits[]);
extern "C" void CPLC_get_cnst_limits(int limits[]);
extern "C" void CPLC_get_bnry_limits(int limits[]);
extern "C" bool CPLC_map_cfd2md_coord(double coord_cfd[], double coord_md[]);
extern "C" bool CPLC_map_md2cfd_coord(double coord_md[], double coord_cfd[]);
extern "C" bool CPLC_overlap();
extern "C" bool CPLC_is_proc_inside(int region[]);

// Setters
extern "C" void CPLC_set_output_mode
( int mode
);

// Getters
extern "C" int CPLC_nsteps_md();
extern "C" int CPLC_nsteps_cfd();
extern "C" int CPLC_nsteps_coupled();
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
extern "C" int CPLC_icmin_bnry();
extern "C" int CPLC_jcmin_bnry();
extern "C" int CPLC_kcmin_bnry();
extern "C" int CPLC_icmax_bnry();
extern "C" int CPLC_jcmax_bnry();
extern "C" int CPLC_kcmax_bnry();
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
extern "C" double CPLC_x_orig_cfd();
extern "C" double CPLC_y_orig_cfd();
extern "C" double CPLC_z_orig_cfd();
extern "C" double CPLC_x_orig_md();
extern "C" double CPLC_y_orig_md();
extern "C" double CPLC_z_orig_md();
extern "C" double CPLC_xl_cfd();
extern "C" double CPLC_yl_cfd();
extern "C" double CPLC_zl_cfd();
extern "C" double* CPLC_xg();
extern "C" double* CPLC_yg();
extern "C" double* CPLC_zg();
