program MD_cpl_example
    use cpl, only : CPL_init, CPL_setup_MD, & 
                    CPL_get_olap_limits, CPL_my_proc_portion, &
                    CPL_get_no_cells, CPL_send, CPL_recv, &
					CPL_finalize
    use mpi
    implicit none

    logical :: recv_flag, send_flag, NO_ERROR
    integer :: i,j,k,ii,jj,kk,ierr
    integer :: NPx, NPy, NPz, NProcs, rank
    integer :: nprocs_realm
    integer :: CART_COMM, MD_COMM
    integer, parameter :: md_realm=2
    integer, dimension(3) :: npxyz, Ncells
    integer, dimension(6) :: portion, limits
    double precision, parameter :: pi = 3.14159265359
    double precision, dimension(3)  :: xyzL, xyz_orig
    double precision, dimension(:,:,:,:), & 
         allocatable  :: send_array

    !Initialise MPI
    call MPI_Init(ierr)

    !Create MD Comm by spliting world
    call CPL_init(MD_realm, MD_COMM, ierr)

    ! Parameters of the cpu topology (cartesian grid)
    xyzL = (/10.d0, 10.d0, 10.d0/)
    xyz_orig = (/0.d0, 0.d0, 0.d0/)
    npxyz = (/ 2, 2, 1/)

    ! Create communicators and check that number of processors is consistent
    call MPI_Comm_size(MD_COMM, nprocs_realm, ierr) 
    if (nprocs_realm .ne. product(npxyz)) then
        print'(4(a,i6))', "Non-coherent number of processes in MD ", nprocs_realm, & 
                " no equal to ",  npxyz(1), " X ", npxyz(2), " X ", npxyz(3)
        call MPI_Abort(MPI_COMM_WORLD, 1, ierr)
    endif

    !Setup cartesian topology
    call MPI_comm_rank(MD_COMM, rank, ierr)
    call MPI_Cart_create(MD_COMM, 3, npxyz, (/1, 1, 1/), & 
                         .true., CART_COMM, ierr)

    !Coupler setup
    call CPL_setup_md(CART_COMM, xyzL, xyz_orig)

    !Get detail for grid
    call CPL_get_olap_limits(limits)
    call CPL_my_proc_portion(limits, portion)
    call CPL_get_no_cells(portion, Ncells)

    ! Pack send_array with cell coordinates. Each cell in the array carries
    ! its global cell number within the overlap region.
    allocate(send_array(1, Ncells(1), Ncells(2), Ncells(3)))
    do i = 1,Ncells(1)
    do j = 1,Ncells(2)
    do k = 1,Ncells(3)
        ! -2 indices to match c++ and python indexing in portion and i,j,k
        ii = i + portion(1) - 2
        jj = j + portion(3) - 2
        kk = k + portion(5) - 2

        send_array(1,i,j,k) = sin(2.d0*pi*ii/Ncells(1)-0.25d0*jj*pi)*cos(2.d0*pi*kk/Ncells(3))
    enddo
    enddo
    enddo
    call CPL_send(send_array, limits, send_flag)

    !Release all coupler comms 
    call CPL_finalize(ierr)

    !Deallocate arrays and finalise MPI
    deallocate(send_array)
    call MPI_finalize(ierr)

end program MD_cpl_example
