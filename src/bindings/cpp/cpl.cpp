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
						 std::string& string_param) {
	char* sp;
	CPLC_get_string_param(section.c_str(), param_name.c_str(), &sp);
	string_param = std::string(sp);
}
//TODO NOTE: These are not working for the moment

// void CPL::get_file_param(const std::string section, const std::string param_name,
// 						 std::vector<bool> &boolean_param_array){
// 	int array_len = 0;
// 	bool* boolean_param_array_c;
// 	CPLC_get_boolean_array_param(section.c_str(), param_name.c_str(),
// 								 &boolean_param_array_c, &array_len);
// 	boolean_param_array.assign(boolean_param_array_c,
// 							   boolean_param_array_c + array_len);
// }

// void CPL::get_file_param(const std::string section, const std::string param_name,
// 						 std::vector<std::string> &string_param_array){
// 	int array_len = 0;
// 	char* string_param_array_c;
// 	CPLC_get_string_array_param(section.c_str(), param_name.c_str(),
// 								&string_param_array_c, &array_len);
// 	string_param_array.resize(array_len);
// 	for (int s = 0; s < array_len; s++) {
// 		std::stringstream myStreamString;
// 		myStreamString << ((char**)string_param_array_c)[s];
// 		string_param_array[s] = myStreamString.str();
// 	}
// 	// Free allocated array inside Fortran bindings. A bit dirty but it does the job.
// 	free(string_param_array_c);
//
// }
#endif 

void CPL::init(int calling_realm, int& returned_realm_comm) {
    CPLC_init(calling_realm, &returned_realm_comm);
}


void CPL::finalize() {
    CPLC_finalize();
}

void CPL::setup_cfd(int icomm_grid, double xyzL[],
                    double xyz_orig[], int ncxyz[]) {
    CPLC_setup_cfd(icomm_grid, xyzL, xyz_orig, ncxyz);
}

void CPL::setup_md(int icomm_grid, double xyzL[], 
                   double xyz_orig[]) {
    CPLC_setup_md(icomm_grid, xyzL, xyz_orig);
}

void CPL::set_timing(int initialstep, int nsteps, double dt) {
    CPLC_set_timing(initialstep, nsteps, dt);
}

bool CPL::send (double* asend, int* asend_shape, int* limits) {
    bool send_flag;
    CPLC_send(asend, asend_shape, limits, &send_flag);
    return send_flag;
}

bool CPL::recv (double* arecv,	int* arecv_shape, int* limits) {
    bool recv_flag;
    CPLC_recv(arecv, arecv_shape, limits, &recv_flag);
    return recv_flag;
}

void CPL::scatter (double* scatterarray, int* scatter_shape, int* limits,
                   double* recvarray, int* recv_shape) {
    CPLC_scatter (scatterarray, scatter_shape, limits, recvarray, recv_shape);
}

void CPL::gather (double* gatherarray, int* gather_shape, int* limits,
                  double* recvarray, int* recv_shape) {
    CPLC_gather (gatherarray, gather_shape, limits, recvarray, recv_shape);
}

void CPL::proc_extents (int coord[], int realm, int extents[]) {
   CPLC_proc_extents (coord, realm, extents);
}

void CPL::my_proc_extents (int extents[]) {
    CPLC_my_proc_extents (extents);
}

void CPL::proc_portion (int coord[], int realm,
                        int limits[], int portion[]) {
    CPLC_proc_portion (coord, realm, limits, portion);
}

void CPL::my_proc_portion (int limits[], int portion[])
{
    CPLC_my_proc_portion (limits, portion);
}

void CPL::map_cell2coord (int i, int j, int k, double coord_xyz[]) {
    CPLC_map_cell2coord (i, j, k, coord_xyz);
}

bool CPL::map_coord2cell (double x, double y, double z, int cell_ijk[]) {
    return CPLC_map_coord2cell (x, y, z, cell_ijk);
}

void CPL::get_no_cells (int limits[], int no_cells[]) {
    CPLC_get_no_cells(limits, no_cells);
}

bool CPL::map_glob2loc_cell (int limits[], int glob_cell[], int loc_cell[]) {
    return CPLC_map_glob2loc_cell(limits, glob_cell, loc_cell);
}

void CPL::get_olap_limits (int limits[]) {
    CPLC_get_olap_limits(limits);
}


void CPL::get_cnst_limits (int limits[]) {
    CPLC_get_cnst_limits(limits);
}

void CPL::get_bnry_limits (int limits[]) {
    CPLC_get_bnry_limits(limits);
}


bool CPL::map_cfd2md_coord (double cfd_coord[], double md_coord[]) {
    return CPLC_map_cfd2md_coord (cfd_coord, md_coord);
}

bool CPL::map_md2cfd_coord (double md_coord[], double cfd_coord[]) {
    return CPLC_map_md2cfd_coord (md_coord, cfd_coord);
}

bool CPL::overlap() {
    return CPLC_overlap();
}

bool CPL::is_proc_inside(int region[]) {
	return CPLC_is_proc_inside(region);
}

// Setters
void CPL::set_output_mode (int mode) {
    CPLC_set_output_mode(mode);
}

// Getters
double CPL::density_cfd() {
    return CPLC_density_cfd();
}
