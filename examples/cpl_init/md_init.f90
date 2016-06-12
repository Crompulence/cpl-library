program main_MD
   use cpl, only : CPL_init, CPL_finalize
   use mpi
   implicit none

   integer :: rank, nprocs, ierr
   integer :: MD_COMM
   integer, parameter :: MD_realm=2

   !Initialise MPI
   call MPI_Init(ierr)

   !Create MD Comm by spliting world
   call CPL_init(MD_realm, MD_COMM, ierr)

   call MPI_comm_size(MD_COMM, nprocs, ierr)
   call MPI_comm_rank(MD_COMM, rank, ierr)

   print*, "MD code processor ", rank+1, " of ", nprocs

   !No need for seperate CPL finalise as MPI finalise takes care of this
   call CPL_finalize(ierr)
   call MPI_comm_free(MD_COMM,ierr)
   call MPI_finalize(ierr)

end program main_MD


