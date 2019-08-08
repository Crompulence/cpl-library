#include "mpi.h"
#include <iostream>

#include "cpl.h"
#include "CPL_ndArray.h"

using namespace std;

int main() {

    bool flag;
    int CFD_realm = 1;
    MPI_Comm CFD_COMM, CART_COMM;
    CPL::ndArray<double> send_array, recv_array;

    MPI_Init(NULL, NULL); 
    CPL::init(CFD_realm, CFD_COMM);

    int npxyz[3] = {1, 1, 1}; int periods[3] = {1, 1, 1};
    MPI_Cart_create(CFD_COMM, 3, npxyz, periods, 1, &CART_COMM);

    //MPI_Cart_get(icomm_grid, 3, npxyz_cfd, cart_periods, cart_coords);

    double xyzL[3] = {1.0, 1.0, 1.0}; double xyz_orig[3] = {0.0, 0.0, 0.0};
    int ncxyz[3] = {32, 32, 32};
    CPL::setup_cfd(CART_COMM, xyzL, xyz_orig, ncxyz);
    CPL::get_arrays(&recv_array, 4, &send_array, 1);

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
