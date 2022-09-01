program MD_cpl_example
    use cpl, only : CPL_init, CPL_setup_MD, CPL_setup_CFD,  & 
                    CPL_get_olap_limits, CPL_my_proc_portion, &
                    CPL_get_no_cells, CPL_send, CPL_recv, &
                    CPL_finalize
    use mpi
    use stuff
    implicit none

    logical :: recv_flag, send_flag, NO_ERROR, error, CFD
    character(len=200)    :: arg, scaling
    integer :: i,j,k,ic,jc,kc,ii,jj,kk,nt,ierr, N, P_CFD, P_MD
    integer :: NPx, NPy, NPz, NProcs, rank
    integer :: nprocs_realm, time
    integer :: CART_COMM, COMM, realm
    integer, parameter :: cfd_realm=1, md_realm=2
    integer, dimension(3) :: npxyz, npxyz_CFD, npxyz_MD, Ncells, ncxyz
    integer, dimension(6) :: portion, limits
    double precision :: t1, t2
    double precision, parameter :: pi = 3.14159265359
    double precision, dimension(3)  :: xyzL, xyz_orig, temp
    double precision, dimension(:,:,:,:), & 
         allocatable  :: send_array, recv_array

    scaling="weak"

    !We cannot specify command arguments in MPMD mode!?
#if TYPE
    arg="MD"
#else
    arg="CFD"
#endif
    !call get_cmd_args(arg)
    print*, "Run of type", arg

    if (arg .eq. "CFD") then
        !print*, "Input is CFD"
        CFD = .true.
        realm = CFD_realm
    elseif (arg .eq. "MD") then
        !print*, "Input is MD"
        CFD = .false.
        realm = MD_realm
    else
        print*, arg
        stop "Error, input flag should be CFD or MD"
    endif


    !Initialise MPI
    call MPI_Init(ierr)

    !Create MD Comm by spliting world
    call CPL_init(realm, COMM, ierr)

    ! Parameters of the cpu topology (cartesian grid)
    xyzL = (/1.d0, 1.d0, 1.d0/)
    xyz_orig = (/0.d0, 0.d0, 0.d0/)
    call read_input('./input', ncxyz, npxyz_CFD, npxyz_MD)
    if (CFD) then
        npxyz = npxyz_CFD
    else
        npxyz = npxyz_MD
    endif

    if (scaling .eq. "weak") then
        ncxyz = ncxyz * npxyz_CFD
    endif


    ! Create communicators and check that number of processors is consistent
    call MPI_Comm_size(COMM, nprocs_realm, ierr) 
    if (nprocs_realm .ne. product(npxyz)) then
        print'(4(a,i6))', "Non-coherent number of processes in MD ", nprocs_realm, & 
                " not equal to ",  npxyz(1), " X ", npxyz(2), " X ", npxyz(3)
        print'(a, i6)', "Attempting to find factors that work with ", nprocs_realm
    call find3factors(nprocs_realm, npxyz(1), npxyz(2), npxyz(3))
        print'(a, 3i6)', "Using: ", npxyz
    endif

    !Setup cartesian topology
    call MPI_comm_rank(COMM, rank, ierr)
    call MPI_Cart_create(COMM, 3, npxyz, (/1, 1, 1/), & 
                         .true., CART_COMM, ierr)

    !Coupler setup
    if (CFD) then
        call CPL_setup_cfd(CART_COMM, xyzL, xyz_orig, ncxyz)
    else
        call CPL_setup_md(CART_COMM, xyzL, xyz_orig)
    endif
    !Get detail for grid
    call CPL_get_olap_limits(limits)
    call CPL_my_proc_portion(limits, portion)
    call CPL_get_no_cells(portion, Ncells)

    allocate(send_array(3, Ncells(1), Ncells(2), Ncells(3)))
    allocate(recv_array(3, Ncells(1), Ncells(2), Ncells(3)))
    send_array = 0.d0; recv_array = 0.d0


    t1 = MPI_wtime()

    !Use 400 steps and 0.05 to match OpenFOAM scaling tests with ~50k cells
    !which is the value stated to give good scaling in
    !(http://www.mcs.anl.gov/~fischer/nek5000/sprague_nek5000_dec2010.pdf)
    do time = 1,20

        !print*, time, Ncells

        if (CFD) then
            !Create array packed with cell numbers 
            call create_array(send_array)
            call CPL_send(send_array, limits, send_flag)
            call check_array(send_array, no_error)

            !Recv array of data
            call CPL_recv(recv_array, limits, recv_flag)
            call check_array(recv_array, no_error)
            if (.not. no_error) stop "Error - CFD recv not correct"
        else
            !Recv array of data
            call CPL_recv(recv_array, limits, recv_flag)
            call check_array(recv_array, no_error)
            if (.not. no_error) stop "Error - MD recv not correct"

            !Create array packed with cell numbers 
            call create_array(send_array)
            call CPL_send(send_array, limits, send_flag)

        endif

        !Strong scaling -- wait the same time per processor
    if (scaling .eq. "strong") then
        do i=1,ncxyz(1)
        do j=1,ncxyz(2)
        do k=1,ncxyz(3)
            call milisleep(0.005d0)
        enddo
        enddo
        enddo
        
        !Weak Scaling -- wait is proportional to ncells on processor
    elseif (scaling .eq. "weak") then
        print*, "Ncells = ", Ncells
        
        do i=1,Ncells(1)
        do j=1,Ncells(2)
        do k=1,Ncells(3)
            !All of the remaining domain needs to include some work
            call milisleep(0.005d0)
        !    call milisleep(0.05d0/dble(Ncells(1)*Ncells(2)*Ncells(3)))
            !Something arbitrary here, access some data with a sqrt
            do ic=-1,1
            do jc=-1,1
            do kc=-1,1
                ii = i + ic
                jj = j + jc
                kk = k + kc
                if (ii .lt. 1) ii = 1
                if (jj .lt. 1) jj = 1
                if (kk .lt. 1) kk = 1
                if (ii .gt. Ncells(1)) ii = Ncells(1)
                if (jj .gt. Ncells(2)) jj = Ncells(2)
                if (kk .gt. Ncells(3)) kk = Ncells(3)
                temp = 0.d0
                do nt=1,5
                        temp=temp+send_array(:,i,j,k)+sqrt(recv_array(:,ii,jj,kk))
                enddo
            enddo
            enddo
            enddo
        enddo
        enddo
        enddo
    else
        stop "Scaling should be set to weak or strong"
    endif

    enddo

    N = product(ncxyz)
    P_CFD = product(npxyz_CFD)
    P_MD = product(npxyz_MD)
    t2 = MPI_wtime()
    if (CFD .and. rank .eq. 0) then
        call write_output("./output"//trim(scaling), N, npxyz_CFD, npxyz_MD, t2-t1)
    endif

    !Release all coupler comms 
    call CPL_finalize(ierr)

    !Deallocate arrays and finalise MPI
    deallocate(send_array)
    call MPI_finalize(ierr)

end program MD_cpl_example
