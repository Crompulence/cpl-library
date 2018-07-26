
program main_CFD
    use cpl, only : CPL_Init, CPL_setup_md, CPL_send, &
                    CPL_recv, CPL_get_arrays, CPL_finalize
    use mpi
    implicit none

    integer :: time, MD_COMM, CART_COMM, ierr, MD_realm=2
    double precision, dimension(:,:,:,:), & 
         allocatable  :: send_array, recv_array

    call MPI_Init(ierr)
    call CPL_init(MD_realm, MD_COMM, ierr)
    call MPI_Cart_create(MD_COMM, 3, (/1, 1, 1/), & 
                        (/.true.,.true.,.true./), & 
                         .true., CART_COMM, ierr)
    call CPL_setup_md(CART_COMM, (/1.d0, 1.d0, 1.d0/), &
                       (/0.d0, 0.d0, 0.d0/))

    call CPL_get_arrays(recv_array, 1, send_array, 4)

    do time=1,5
        send_array(1,:,:,:) = 5.*time
        call CPL_send(send_array)
        call CPL_recv(recv_array)
        print*, "MD", time, recv_array(1,1,1,1)
    enddo

   call CPL_finalize(ierr)
   call MPI_comm_free(MD_COMM,ierr)
   call MPI_finalize(ierr)

end program main_CFD
