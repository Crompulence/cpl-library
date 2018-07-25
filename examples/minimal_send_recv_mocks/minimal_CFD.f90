
program main_CFD
    use cpl, only : CPL_Init, CPL_setup_cfd, CPL_send, &
                    CPL_recv, CPL_get_arrays, CPL_finalize
    use mpi
    implicit none

    integer :: time, CFD_COMM, CART_COMM, ierr, CFD_realm=1
    double precision, dimension(:,:,:,:), & 
         allocatable  :: send_array, recv_array

    call MPI_Init(ierr)
    call CPL_init(CFD_realm, CFD_COMM, ierr)
    call MPI_Cart_create(CFD_COMM, 3, (/1, 1, 1/), & 
                        (/.true.,.true.,.true./), & 
                         .true., CART_COMM, ierr)
    call CPL_setup_cfd(CART_COMM, (/1.d0, 1.d0, 1.d0/), &
                       (/0.d0, 0.d0, 0.d0/), &
                       (/32, 32, 32/))

    call CPL_get_arrays(recv_array, 4, send_array, 1)

    do time=1,5
        call CPL_recv(recv_array)
        print*, "CFD", time, recv_array(1,1,1,1)
        send_array(1,:,:,:) = 2.*time
        call CPL_send(send_array)
    enddo

   call CPL_finalize(ierr)
   call MPI_comm_free(CFD_COMM,ierr)
   call MPI_finalize(ierr)

end program main_CFD
