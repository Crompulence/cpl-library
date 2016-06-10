#include "CPL.h"
#include "mpi.h"
#include <iostream>

using namespace std;

int main() {
   MPI_Init(NULL, NULL); 

   int CFD_realm = 1, CFD_COMM;
   CPL::init(CFD_realm, CFD_COMM);

   // Parameters of the cpu topology (cartesian grid)
   double xyzL[3] = {10.0, 10.0, 10.0};
   double xyz_orig[3] = {0.0, 0.0, 0.0};
   int npxyz[3] = {2, 2, 1};
   int ncxyz[3] = {64, 18, 64};

   int nprocs_realm;
   MPI_Comm_size(CFD_COMM, &nprocs_realm);

   // Create communicators and check that number of processors is consistent
   if (nprocs_realm != (npxyz[0] * npxyz[1] * npxyz[2])) {
      cout << "Non-coherent number of processes." << endl;
      MPI_Abort(MPI_COMM_WORLD, 1);
   }

   // Setup cartesian topology
   int rank;
   MPI_Comm_rank(CFD_COMM, &rank);
   int periods[3] = {1, 1, 1};
   int CART_COMM;
   MPI_Cart_create(CFD_COMM, 3, npxyz, periods, true, &CART_COMM);

   // Coupler setup
   CPL::setup_cfd(CART_COMM, xyzL, xyz_orig, ncxyz);

   // Get detail for grid
   int Ncells[3];
   int olap_limits[6], portion[6];
   CPL::get_olap_limits(olap_limits);
   CPL::my_proc_portion(olap_limits, portion);
   CPL::get_no_cells(portion, Ncells);

   // Pack send_array with cell coordinates. Each cell in the array carries
   // its global cell number within the overlap region.
   int send_shape[4] = {3, Ncells[0], Ncells[1], Ncells[2]};
   CPL::ndArray<double> send_array(4, send_shape);
   int ii, jj, kk;
   for (int i = 0; i < Ncells[0]; i++) {
       ii = i + portion[0];
       for (int j = 0; j < Ncells[1]; j++) {
           jj = j + portion[2];
           for (int k = 0; k < Ncells[2]; k++) {
               kk = k + portion[4];

               send_array(0, i, j, k) = (double) ii;
               send_array(1, i, j, k) = (double) jj;
               send_array(2, i, j, k) = (double) kk;
           }
        }
    }

   CPL::send(send_array.data(), send_array.shapeData(), olap_limits);

   // Block before checking if successful
   MPI_Barrier(MPI_COMM_WORLD);

   // Release all coupler comms 
   CPL::finalize();
   MPI_free_comm(CFD_COMM)
   MPI_free_comm(CART_COMM)
   MPI_Finalize();
   
}
