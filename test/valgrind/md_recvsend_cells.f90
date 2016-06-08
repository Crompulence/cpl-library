program md_cpl_example
    use cpl, only : CPL_init, CPL_setup_md, & 
                    CPL_get_olap_limits, CPL_my_proc_portion, &
                    CPL_get_no_cells, CPL_send, CPL_recv, &
					CPL_overlap, CPL_finalize, CPL_overlap
    use array_stuff
    use mpi
    implicit none

    logical :: recv_flag,send_flag, NO_ERROR
    integer :: i,j,k,ii,jj,kk,ierr,errorcode
    integer :: rank, nprocs_realm, nsteps, initialstep
    integer :: CART_COMM, MD_COMM
    integer, parameter :: md_realm=2
    integer, dimension(3) :: npxyz, Ncells
    integer, dimension(6) :: portion, limits
    double precision :: dt, density
    double precision, dimension(3)  :: xyzL, xyz_orig
    double precision, dimension(:,:,:,:), & 
        allocatable  :: recv_array, send_array

    !Initialise MPI
    call MPI_Init(ierr)

    !Create MD Comm by spliting world
    call CPL_init(md_realm, MD_COMM, ierr)

    !Parameters
    dt = 0.002
    density = 0.8
	nsteps = 5000
	initialstep = 0

    ! Parameters of the cpu topology (cartesian grid)
    call read_input(xyzL=xyzL, xyz_orig=xyz_orig, & 
                    npxyz_MD=npxyz)

    ! Create communicators and check that number of processors is consistent
    call MPI_Comm_size(MD_COMM, nprocs_realm, ierr) 

    if (nprocs_realm .ne. product(npxyz)) then
        print'(4(a,i6))', "Non-coherent number of processes in MD ", nprocs_realm, & 
                " no equal to ",  npxyz(1), " X ", npxyz(2), " X ", npxyz(3)
        call MPI_Abort(MPI_COMM_WORLD, errorcode, ierr)
    endif

    !Setup cartesian topology
    call MPI_comm_rank(MD_COMM, rank, ierr)
    call MPI_Cart_create(MD_COMM, 3, npxyz, (/.true.,.true.,.true./), & 
                         .true., CART_COMM, ierr)

    !Coupler setup
    call CPL_setup_md(nsteps, initialstep, dt, CART_COMM, & 
                      xyzL, xyz_orig, density)

    !Get detail for grid
    call CPL_get_olap_limits(limits)
    call CPL_my_proc_portion(limits, portion)
    call CPL_get_no_cells(portion, Ncells)

    !Coupled Recieve and print
    allocate(recv_array(3, Ncells(1), Ncells(2), Ncells(3)))
    recv_array = 0.d0
    call CPL_recv(recv_array, limits, recv_flag)
    call print_array(recv_array, rank, no_error)

    !Block before checking if successful
    call MPI_Barrier(MPI_COMM_WORLD, ierr)
    if (CPL_overlap() .and. no_error) then
        print'(a,a,i2,a)', "MD -- ", "(rank=", rank, ") CELLS HAVE BEEN RECEIVED CORRECTLY."
    endif

    !Coupled Send
    call fill_array(send_array)
    call CPL_send(send_array, limits, send_flag)

    !Release all coupler comms 
    call MPI_Barrier(MPI_COMM_WORLD, ierr)
    call CPL_finalize(ierr)

    !Deallocate arrays and finalise MPI
    deallocate(recv_array, send_array)
	call MPI_COMM_FREE(CART_COMM, ierr)
	call MPI_COMM_FREE(MD_COMM, ierr)
    call MPI_finalize(ierr)

end program
