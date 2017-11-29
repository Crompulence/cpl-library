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
#include "cpl.h"
#include <iostream>
#include <sstream>

// Macro definitions
#define CONSTINT(X) (const_cast<int*>(X))
#define CONSTDOUBLE(X) (const_cast<double*>(X))

#ifdef JSON_SUPPORT
void CPL::load_param_file(std::string fname) {
	CPLC_load_param_file(fname.c_str());
}

void CPL::close_param_file() {
	CPLC_close_param_file();
}

void CPL::get_file_param(const std::string section, const std::string param_name,
					 	 double& real_param) {
	CPLC_get_real_param(section.c_str(), param_name.c_str(), &real_param);
}

void CPL::get_file_param(const std::string section, const std::string param_name,
		 				 std::vector<double> &real_param_array){
	int array_len = 0;
	double* real_param_array_c;
	CPLC_get_real_array_param(section.c_str(), param_name.c_str(),
							  &real_param_array_c, &array_len);
	real_param_array.assign(real_param_array_c, real_param_array_c + array_len);
}

void CPL::get_file_param(const std::string section, const std::string param_name,
						 int& int_param) {
	CPLC_get_int_param(section.c_str(), param_name.c_str(), &int_param);
}

void CPL::get_file_param(const std::string section, const std::string param_name,
						 std::vector<int> &int_param_array){
	int array_len = 0;
	int* int_param_array_c;
	CPLC_get_int_array_param(section.c_str(), param_name.c_str(),
							 &int_param_array_c, &array_len);
	int_param_array.assign(int_param_array_c, int_param_array_c + array_len);
}

void CPL::get_file_param(const std::string section, const std::string param_name,
						 bool& boolean_param) {
	CPLC_get_boolean_param(section.c_str(), param_name.c_str(), &boolean_param);
}

void CPL::get_file_param(const std::string section, const std::string param_name,
						 std::vector<bool> &boolean_param_array){
	int array_len = 0;
	bool* boolean_param_array_c;
	CPLC_get_boolean_array_param(section.c_str(), param_name.c_str(),
								 &boolean_param_array_c, &array_len);
	boolean_param_array.assign(boolean_param_array_c,
							   boolean_param_array_c + array_len);
}

void CPL::get_file_param(const std::string section, const std::string param_name,
						 std::string& string_param) {
	char* sp;
	CPLC_get_string_param(section.c_str(), param_name.c_str(), &sp);
	string_param = std::string(sp);
}

void CPL::get_file_param(const std::string section, const std::string param_name,
						 std::vector<std::string> &string_param_array){
	int array_len = 0;
	char* string_param_array_c;
	CPLC_get_string_array_param(section.c_str(), param_name.c_str(),
								&string_param_array_c, &array_len);
	string_param_array.resize(array_len);
	for (int s = 0; s < array_len; s++) {
		std::stringstream myStreamString;
		myStreamString << ((char**)string_param_array_c)[s];
		string_param_array[s] = myStreamString.str();
	}
	// Free allocated array inside Fortran bindings. A bit dirty but it does the job.
	free(string_param_array_c);

}
#endif 


void CPL::init(int calling_realm, int& returned_realm_comm) {
    CPLC_init(calling_realm, &returned_realm_comm);
}


void CPL::finalize() {
    CPLC_finalize();
}

void CPL::setup_cfd(int icomm_grid, const CPL::DoubVector& xyzL, const CPL::DoubVector& xyz_orig, 
               const CPL::IntVector& ncxyz) {
    CPLC_setup_cfd(icomm_grid, CONSTDOUBLE(&xyzL[0]),
                   CONSTDOUBLE(&xyz_orig[0]), CONSTINT(&ncxyz[0]));
}

void CPL::setup_md(int icomm_grid, const CPL::DoubVector& xyzL, 
                   const CPL::DoubVector& xyz_orig) {
    CPLC_setup_md(icomm_grid, CONSTDOUBLE(&xyzL[0]), CONSTDOUBLE(&xyz_orig[0]));
}

void CPL::set_timing(int initialstep, int nsteps, double dt) {
    CPLC_set_timing(initialstep, nsteps, dt);
}

bool CPL::send (const CPL::DoubNdArray& asend, const CPL::IntVector& asend_shape, 
                const CPL::IntVector& limits) {
    bool send_flag;
    CPLC_send(CONSTDOUBLE(asend.data()), CONSTINT(&asend_shape[0]), CONSTINT(&limits[0]), &send_flag);
    return send_flag;
}

bool CPL::recv (CPL::DoubNdArray& arecv, const CPL::IntVector& arecv_shape, const CPL::IntVector& limits) {
    bool recv_flag; 
    CPLC_recv(arecv.data(), CONSTINT(&arecv_shape[0]), CONSTINT(&limits[0]), &recv_flag);
    return recv_flag;
}

void CPL::scatter (CPL::DoubVector& scatterarray, const CPL::IntVector& scatter_shape, 
                   const CPL::IntVector& limits, CPL::DoubVector& recvarray, const CPL::IntVector& recv_shape) {
    CPLC_scatter (&scatterarray[0], CONSTINT(&scatter_shape[0]), CONSTINT(&limits[0]), &recvarray[0], 
                  CONSTINT(&recv_shape[0]));
}

void CPL::gather (CPL::DoubVector& gatherarray, const CPL::IntVector& gather_shape,
                  const CPL::IntVector& limits, CPL::DoubVector& recvarray, const CPL::IntVector& recv_shape) {
    CPLC_gather (&gatherarray[0], CONSTINT(&gather_shape[0]), CONSTINT(&limits[0]), &recvarray[0], 
                 CONSTINT(&recv_shape[0]));
}

void CPL::proc_extents (const CPL::IntVector& coord, int realm, CPL::IntVector& extents) {
   CPLC_proc_extents (CONSTINT(&coord[0]), realm, &extents[0]);
}

void CPL::my_proc_extents (CPL::IntVector& extents) {
    CPLC_my_proc_extents (&extents[0]);
}

void CPL::proc_portion (const CPL::IntVector& coord, int realm,
                        const CPL::IntVector& limits, CPL::IntVector& portion) {
    CPLC_proc_portion (CONSTINT(&coord[0]), realm, CONSTINT(&limits[0]), &portion[0]);
}

void CPL::my_proc_portion (const CPL::IntVector& limits, CPL::IntVector& portion) {
    CPLC_my_proc_portion (CONSTINT(&limits[0]), &portion[0]);
}

void CPL::map_cell2coord (int i, int j, int k, CPL::DoubVector& coord_xyz) {
    CPLC_map_cell2coord (i, j, k, &coord_xyz[0]);
}

bool CPL::map_coord2cell (double x, double y, double z, CPL::IntVector& cell_ijk) {
    return CPLC_map_coord2cell (x, y, z, &cell_ijk[0]);
}

void CPL::get_no_cells (const CPL::IntVector& limits, CPL::IntVector& no_cells) {
    CPLC_get_no_cells(CONSTINT(&limits[0]), &no_cells[0]);
}

bool CPL::map_glob2loc_cell (const CPL::IntVector& limits, const CPL::IntVector& glob_cell,
                             CPL::IntVector& loc_cell) {
    return CPLC_map_glob2loc_cell(CONSTINT(&limits[0]), CONSTINT(&glob_cell[0]), &loc_cell[0]);
}

void CPL::get_olap_limits (CPL::IntVector& limits) {
    CPLC_get_olap_limits(&limits[0]);
}


void CPL::get_cnst_limits (CPL::IntVector& limits) {
    CPLC_get_cnst_limits(&limits[0]);
}

void CPL::get_bnry_limits (CPL::IntVector& limits) {
    CPLC_get_bnry_limits(&limits[0]);
}

bool CPL::map_cfd2md_coord (const CPL::DoubVector& cfd_coord, CPL::DoubVector& md_coord) {
    return CPLC_map_cfd2md_coord (CONSTDOUBLE(&cfd_coord[0]), &md_coord[0]);
}

bool CPL::map_md2cfd_coord (const CPL::DoubVector& md_coord, CPL::DoubVector& cfd_coord) {
    return CPLC_map_md2cfd_coord (CONSTDOUBLE(&md_coord[0]), &cfd_coord[0]);
}

bool CPL::overlap() {
    return CPLC_overlap();
}

bool CPL::is_proc_inside(const CPL::IntVector& region) {
	return CPLC_is_proc_inside(CONSTINT(&region[0]));
}

// Setters
void CPL::set_output_mode (int mode) {
    CPLC_set_output_mode(mode);
}

// Getters
double CPL::density_cfd() {
    return CPLC_density_cfd();
}
