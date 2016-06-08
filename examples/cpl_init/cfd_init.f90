program main_CFD
   use cpl, only : CPL_init, CPL_finalize
   use mpi
   implicit none

   integer :: rank, nprocs, ierr
   integer :: CFD_COMM
   integer, parameter :: CFD_realm=1

   !Initialise MPI
   call MPI_Init(ierr)

   !Create MD Comm by spliting world
   call CPL_init(CFD_realm, CFD_COMM, ierr)

   !get local processor rank and print
   call MPI_comm_size(CFD_COMM, nprocs, ierr)
   call MPI_comm_rank(CFD_COMM, rank, ierr)

   print*, "CFD code processor ", rank+1, " of ", nprocs

   !No need for seperate CPL finalise as MPI finalise takes care of this
   call CPL_finalize(ierr)
   call MPI_finalize(ierr)

end program main_CFD
