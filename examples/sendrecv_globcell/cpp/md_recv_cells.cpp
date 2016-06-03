#include "CPL.h"
#include "mpi.h"
#include <iostream>
#include <stdio.h>
#include <string>

using namespace std;

int main() {
   MPI_Init(NULL, NULL); 

   int MD_realm = 2, MD_COMM;
   CPL::init(MD_realm, MD_COMM);

   // Parameters
   int dt = 0.1;
   double density = 0.8;
   int nsteps, initialstep;

   // Parameters of the cpu topology (cartesian grid)
   double xyzL[3] = {10.0, 10.0, 10.0};
   double xyz_orig[3] = {0.0, 0.0, 0.0};
   int npxyz[3] = {4, 2, 2};

   int nprocs_realm;
   MPI_Comm_size(MD_COMM, &nprocs_realm);

   // Create communicators and check that number of processors is consistent
   if (nprocs_realm != (npxyz[0] * npxyz[1] * npxyz[2])) {
      cout << "Non-coherent number of processes." << endl;
      MPI_Abort(MPI_COMM_WORLD, 1);
   }

   // Setup cartesian topology
   int rank;
   MPI_Comm_rank(MD_COMM, &rank);
   int periods[3] = {1, 1, 1};
   int CART_COMM;
   MPI_Cart_create(MD_COMM, 3, npxyz, periods, true, &CART_COMM);

   // Coupler setup
   CPL::setup_md(nsteps, initialstep, dt, CART_COMM, xyzL, xyz_orig, density);

   // Get detail for grid
   int Ncells[3];
   int olap_limits[6], portion[6];
   CPL::get_olap_limits(olap_limits);
   CPL::my_proc_portion(olap_limits, portion);
   CPL::get_no_cells(portion, Ncells);

   // Pack recv_array with cell coordinates. Each cell in the array carries
   // its global cell number within the overlap region.
   int recv_shape[4] = {3, Ncells[0], Ncells[1], Ncells[2]};
   CPL::ndArray<double> recv_array(4, recv_shape);
   CPL::recv(recv_array.data(), recv_array.shapeData(), olap_limits);
   bool no_error = true;
   if (CPL::overlap()) {
       int ii, jj, kk;
       for (int i = 0; i < Ncells[0]; i++) {
           ii = i + portion[0];
           for (int j = 0; j < Ncells[1]; j++) {
               jj = j + portion[2];
               for (int k = 0; k < Ncells[2]; k++) {
                   kk = k + portion[4];
                   if (((double) ii - recv_array(0, i, j, k)) > 1e-8) {
                       printf("Error -- portion in x: %d %d MD rank: %d cell i: %d recv_array: %f\n",\
                               portion[0], portion[1], rank, ii, recv_array(0, i, j, k));
                       no_error = false;
                   }
                   if (((double) jj - recv_array(1, i, j, k)) > 1e-8) {
                       printf("Error -- portion in y: %d %d MD rank: %d cell j: %d recv_array: %f\n",\
                               portion[2], portion[3], rank, jj, recv_array(1, i, j, k));
                       no_error = false;
                   }
                   if (((double) kk - recv_array(2, i, j, k)) > 1e-8) {
                       printf("Error -- portion in z: %d %d MD rank: %d cell k: %d recv_array: %f\n",\
                               portion[1], portion[5], rank, kk, recv_array(2, i, j, k));
                       no_error = false;
                   }
               }
            }
        }
    }

   // Block before checking if successful
   MPI_Barrier(MD_COMM);
   if (CPL::overlap() && no_error)
       printf("MD -- (rank=%2d) CELLS HAVE BEEN RECEIVED CORRECTLY.\n", rank);
   MPI_Barrier(MPI_COMM_WORLD);
}
