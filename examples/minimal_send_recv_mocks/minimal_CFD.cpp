#include "mpi.h"
#include <iostream>

#include "cpl.h"
#include "CPL_ndArray.h"

using namespace std;

void get_arrays(CPL::ndArray<double>* recv_array, int recv_size, 
                CPL::ndArray<double>* send_array, int send_size)
{
    int Ncells[3]; int olap_limits[6], portion[6];
    CPL::get_olap_limits(olap_limits);
    CPL::my_proc_portion(olap_limits, portion);
    CPL::get_no_cells(portion, Ncells);

    int send_shape[4] = {send_size, Ncells[0], Ncells[1], Ncells[2]};
    send_array->resize(4, send_shape);
    int recv_shape[4] = {recv_size, Ncells[0], Ncells[1], Ncells[2]};
    recv_array->resize(4, recv_shape);

}

int main() {

    bool flag;
    int CFD_realm = 1, CFD_COMM;
    CPL::ndArray<double> send_array, recv_array;

    MPI_Init(NULL, NULL); 
    CPL::init(CFD_realm, CFD_COMM);

    int npxyz[3] = {1, 1, 1}; int periods[3] = {1, 1, 1}; int CART_COMM;
    MPI_Cart_create(CFD_COMM, 3, npxyz, periods, true, &CART_COMM);
    double xyzL[3] = {1.0, 1.0, 1.0}; double xyz_orig[3] = {0.0, 0.0, 0.0};
    int ncxyz[3] = {32, 32, 32};
    CPL::setup_cfd(CART_COMM, xyzL, xyz_orig, ncxyz);
    get_arrays(&recv_array, 4, &send_array, 1);

	for (int time = 0; time < 5; time++) {
        flag = CPL::recv(&recv_array);
        std::cout << "CFD " << time << " " << recv_array(0,0,0,0) << std::endl;
        send_array(0,0,0,0) = 2.*time;
        flag = CPL::send(&send_array);
    }

   // Release all coupler comms 
    CPL::finalize();
    MPI_Finalize();
   
}
