#include "mpi.h"
#include <iostream>

#include "cpl.h"
#include "CPL_ndArray.h"

using namespace std;

int main() {

    bool flag;
    int MD_realm = 2, MD_COMM;
    CPL::ndArray<double> send_array, recv_array;

    MPI_Init(NULL, NULL); 
    CPL::init(MD_realm, MD_COMM);

    int npxyz[3] = {1, 1, 1}; int periods[3] = {1, 1, 1}; int CART_COMM;
    MPI_Cart_create(MD_COMM, 3, npxyz, periods, true, &CART_COMM);
    double xyzL[3] = {1.0, 1.0, 1.0}; double xyz_orig[3] = {0.0, 0.0, 0.0};
    CPL::setup_md(CART_COMM, xyzL, xyz_orig);
    CPL::get_arrays(&recv_array, 1, &send_array, 4);

	for (int time = 0; time < 5; time++) {
        send_array(0,0,0,0) = 5.*time;
        flag = CPL::send(&send_array);
        flag = CPL::recv(&recv_array);
        std::cout << "MD " << time << " " << recv_array(0,0,0,0) << std::endl;
    }

   // Release all coupler comms 
    CPL::finalize();
    MPI_Finalize();
   
}
