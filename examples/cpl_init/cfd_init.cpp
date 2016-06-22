#include "CPL.h"
#include "mpi.h"
#include <iostream>

using namespace std;

int main() {
   MPI_Init(NULL, NULL); 

   int CFD_realm = 1, CFD_COMM;
   CPL::init(CFD_realm, CFD_COMM);

   int nprocs; int rank;
   MPI_Comm_size(CFD_COMM, &nprocs);   
   MPI_Comm_rank(CFD_COMM, &rank);

   cout << "MD code processor " << rank+1 << " of " << nprocs << endl;

   CPL::finalize();
   MPI_Finalize();
}



