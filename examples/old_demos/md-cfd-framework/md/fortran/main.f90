program main
    use mpi
    use cpl, only : md_realm, CPL_init
    implicit none

    integer :: rank, ierror, MD_COMM

    call MPI_Init(ierror)
    call CPL_init(md_realm, MD_COMM, ierror)
    call MPI_Finalize(ierror)

end program main 
