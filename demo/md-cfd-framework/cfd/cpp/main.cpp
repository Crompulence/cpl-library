#include "mpi.h"
#include "CPL.h"

int main (int argc, char** argv)
{
    MPI_Init (&argc, &argv);

    MPI_Comm realmComm;
    CPL::create_comm (CPL::cfd_realm, realmComm);

    MPI_Finalize();
}
