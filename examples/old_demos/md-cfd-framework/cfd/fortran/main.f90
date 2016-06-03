program main 
    use mpi
    use cpl, only : cfd_realm, CPL_init
    implicit none

    integer :: rank, ierror, CFD_COMM

    call MPI_Init(ierror)
    call CPL_init(cfd_realm, CFD_COMM, ierror)
    call MPI_Finalize(ierror)

end program main 
