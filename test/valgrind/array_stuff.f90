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

    end subroutine print_array


end module array_stuff
