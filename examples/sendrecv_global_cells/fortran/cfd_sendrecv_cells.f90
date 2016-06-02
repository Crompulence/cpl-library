

module array_stuff
    implicit none

contains

    subroutine read_input(xyzL, xyz_orig, npxyz_MD, ncxyz, npxyz_CFD)

        integer, dimension(3), optional, intent(out) :: npxyz_MD, npxyz_CFD, ncxyz
        double precision, dimension(3), intent(out)  :: xyzL, xyz_orig

        integer :: dummy

        open(unit=10, file='./fortran.in')
        read(10, *) xyzL(1)
        read(10, *) xyzL(2)
        read(10, *) xyzL(3)
        read(10, *) xyz_orig(1)
        read(10, *) xyz_orig(2)
        read(10, *) xyz_orig(3)
        if (present(npxyz_CFD)) then
            read(10, *) npxyz_CFD(1)
            read(10, *) npxyz_CFD(2)
            read(10, *) npxyz_CFD(3)
        else
            read(10, *) dummy
            read(10, *) dummy
            read(10, *) dummy
        endif
        if (present(ncxyz)) then
            read(10, *) ncxyz(1)
            read(10, *) ncxyz(2)
            read(10, *) ncxyz(3)
        else
            read(10, *) dummy
            read(10, *) dummy
            read(10, *) dummy
        endif
        if (present(npxyz_MD)) then
            read(10, *) npxyz_MD(1)
            read(10, *) npxyz_MD(2)
            read(10, *) npxyz_MD(3)
        else
            read(10, *) dummy
            read(10, *) dummy
            read(10, *) dummy
        endif
        close(10)

        !print*, xyzL, xyz_orig, npxyz_MD, ncxyz, npxyz_CFD

    end subroutine 

    subroutine fill_array(A)
        use cpl, only : CPL_get_olap_limits, CPL_my_proc_portion, &
                        CPL_get_no_cells

        double precision, dimension(:,:,:,:), &
            allocatable, intent(out)  :: A

        logical :: no_error
        integer :: i,j,k,ii,jj,kk,ierr
        integer, dimension(3) :: Ncells
        integer, dimension(6) :: portion, limits

        call CPL_get_olap_limits(limits)
        call CPL_my_proc_portion(limits, portion)
        call CPL_get_no_cells(portion, Ncells)

        allocate(A(3, Ncells(1), Ncells(2), Ncells(3)))
        do i = 1,Ncells(1)
        do j = 1,Ncells(2)
        do k = 1,Ncells(3)
            ii = i + portion(1)
            jj = j + portion(3)
            kk = k + portion(5)

            A(1,i,j,k) = ii
            A(2,i,j,k) = jj
            A(3,i,j,k) = kk
        enddo
        enddo
        enddo

    end subroutine fill_array


    subroutine print_array(A, rank, no_error)
        use cpl, only : CPL_get_olap_limits, CPL_my_proc_portion, &
                        CPL_get_no_cells

        logical, intent(out), optional :: no_error
        integer, intent(in), optional :: rank
        double precision, dimension(:,:,:,:), &
            allocatable, intent(in)  :: A

        integer :: i,j,k,ii,jj,kk,ierr
        integer, dimension(3) :: Ncells
        integer, dimension(6) :: portion, limits

        call CPL_get_olap_limits(limits)
        call CPL_my_proc_portion(limits, portion)
        call CPL_get_no_cells(portion, Ncells)

        if (any(portion .ne. -666)) then
            no_error = .true.
            do i = 1,Ncells(1)
            do j = 1,Ncells(2)
            do k = 1,Ncells(3)
                ii = i + portion(1)
                jj = j + portion(3)
                kk = k + portion(5)

                if ((dble(ii) - A(1,i,j,k)) .gt. 1e-8) then 
                    print'(a,2i5,a,i5,a,i6,a,f10.5)', "ERROR -- portion: ", portion(1:2), & 
                           " MD rank: ", rank, " cell i: ",ii, & 
                           " recv_array: ", A(1,i,j,k)
                    no_error = .false.
                endif
                if ((dble(jj) - A(2,i,j,k)) .gt. 1e-8) then 
                    print'(a,2i5,a,i5,a,i6,a,f10.5)', "ERROR -- portion: ", portion(3:4), & 
                           " MD rank: ", rank, " cell j: ", jj , & 
                           " recv_array: ", A(2,i,j,k)
                    no_error = .false.  
                endif
                if ((dble(kk) - A(3,i,j,k)) .gt. 1e-8) then 
                    print'(a,2i5,a,i5,a,i6,a,f10.5)', "ERROR -- portion: ", portion(5:6), & 
                           " MD rank: ", rank, " cell k: ", kk , & 
                           " recv_array: ", A(3,i,j,k)
                    no_error = .false.
                endif
            enddo
            enddo
            enddo

        endif

        if (no_error) print*, "CFD -- NO PROBLEMS HAVE BEEN DETECTED!"

    end subroutine print_array


end module array_stuff



program cfd_cpl_example
    use cpl, only : CPL_init, CPL_setup_cfd, & 
                    CPL_get_olap_limits, CPL_my_proc_portion, &
                    CPL_get_no_cells, CPL_send, CPL_recv
    use array_stuff, only : fill_array, print_array
    use mpi
    implicit none

    logical :: recv_flag, send_flag, NO_ERROR
    integer :: i,j,k,ii,jj,kk,ierr,errorcode
    integer :: NPx, NPy, NPz, NProcs, rank
    integer :: nprocs_realm, nsteps, initialstep
    integer :: CART_COMM, CFD_COMM
    integer, parameter :: cfd_realm=1
    integer, dimension(3) :: npxyz, Ncells, ncxyz
    integer, dimension(6) :: portion, limits
    double precision :: dt, density
    double precision, dimension(3)  :: xyzL, xyz_orig
    double precision, dimension(:,:,:,:), & 
         allocatable  :: send_array, recv_array

    !Initialise MPI
    call MPI_Init(ierr)

    !Create MD Comm by spliting world
    call CPL_init(CFD_realm, CFD_COMM, ierr)

    !Parameters
    dt = 0.1
    density = 0.8
    nsteps = 100

    ! Parameters of the cpu topology (cartesian grid)
    xyzL = (/10.d0, 10.d0, 10.d0/)
    xyz_orig = (/0.d0, 0.d0, 0.d0/)
    npxyz = (/ 2, 2, 1/)
    ncxyz = (/ 64, 18, 64 /)

    !call read_input(xyzL=xyzL, xyz_orig=xyz_orig, & 
    !                npxyz_CFD=npxyz, ncxyz=ncxyz)

    ! Create communicators and check that number of processors is consistent
	call MPI_Comm_size(CFD_COMM, nprocs_realm, ierr) 
    if (nprocs_realm .ne. product(npxyz)) then
        print*, "Non-coherent number of processes"
        call MPI_Abort(MPI_COMM_WORLD, errorcode, ierr)
    endif

    !Setup cartesian topology
    call MPI_comm_rank(CFD_COMM, rank, ierr)
	call MPI_Cart_create(CFD_COMM, 3, npxyz, (/.true.,.true.,.true./), & 
                         .true., CART_COMM, ierr)

    !Coupler setup
    call CPL_setup_cfd(nsteps, dt, CART_COMM, xyzL, xyz_orig, ncxyz, density)

    !Get detail for grid
    call CPL_get_olap_limits(limits)
    call CPL_my_proc_portion(limits, portion)
    call CPL_get_no_cells(portion, Ncells)

    !Pack array with cell data
    call fill_array(send_array)
    call CPL_send(send_array, limits, send_flag)

    !Block before checking if successful
    call MPI_Barrier(MPI_COMM_WORLD, ierr)

    !Coupled Recieve and print
    allocate(recv_array(3, Ncells(1), Ncells(2), Ncells(3)))
    recv_array = 0.d0
    call CPL_recv(recv_array, limits, recv_flag)
    call print_array(recv_array, rank, no_error)

    !Block before checking if successful
    call MPI_Barrier(MPI_COMM_WORLD, ierr)
    if (no_error) print*, "CFD -- NO PROBLEMS HAVE BEEN DETECTED!"

end program cfd_cpl_example
