
module stuff
    implicit none

    contains

    subroutine read_input(filename, ncxyz, npxyz_CFD, npxyz_MD)

        character(*), intent(in) :: filename
        integer :: i
        integer, dimension(3), intent(out) :: ncxyz, npxyz_CFD, npxyz_MD

        open(200,file=trim(filename))
        do i = 1,3
            read(200,*) ncxyz(i)
        enddo
        do i = 1,3
            read(200,*) npxyz_CFD(i)
        enddo
        !Remote realm number of processes
        do i = 1,3
            read(200,*) npxyz_MD(i)
        enddo
        close(200)

    end subroutine read_input

    subroutine write_output(filename, N, P_CFD, P_MD, dt)

        logical                      :: exists
        character(*), intent(in)     :: filename
        integer, intent(in)          :: N
        integer, dimension(3), intent(in) :: P_CFD, P_MD
        double precision, intent(in) :: dt

        inquire(file="output", exist=exists)
        if (exists) then 
            open(100, file='output', status="old", & 
                      position="append", action="write")
        else
            open(100, file='output', status="new", action="write")
        endif
        write(100,'(i12,6i8,f18.5)') N, P_CFD, P_MD, dt
        close(100)

    end subroutine write_output


    subroutine get_cmd_args(arg)
        implicit none

        integer :: i,argcount
        character(len=200), intent(out) :: arg

        argcount = command_argument_count()
        if (argcount .gt. 1) stop "Only one arg needed"

        if (argcount.gt.0) then         !If more than 0 arguments   
            do i=1,argcount         !Loop through all arguments...
                call get_command_argument(i,arg) 
            enddo
        endif

    end subroutine get_cmd_args


    subroutine create_array(send_array)
        use cpl, only : CPL_overlap, CPL_get_olap_limits, CPL_my_proc_portion, CPL_get_no_cells

        double precision, dimension(:,:,:,:), intent(out) :: send_array
        integer :: i,j,k,ii,jj,kk

        integer, dimension(3) :: Ncells
        integer, dimension(6) :: limits, portion

        call CPL_get_olap_limits(limits)
        call CPL_my_proc_portion(limits, portion)
        call CPL_get_no_cells(portion, Ncells)

        ! Pack send_array with cell coordinates. Each cell in the array carries
        ! its global cell number within the overlap region.
        do i = 1,Ncells(1)
        do j = 1,Ncells(2)
        do k = 1,Ncells(3)
            ! -2 indices to match c++ and python indexing in portion and i,j,k
            ii = i + portion(1) - 2
            jj = j + portion(3) - 2
            kk = k + portion(5) - 2

            send_array(1,i,j,k) = ii
            send_array(2,i,j,k) = jj
            send_array(3,i,j,k) = kk
                                  
        enddo
        enddo
        enddo

    end subroutine create_array


    subroutine check_array(recv_array, no_error)
        use cpl, only : CPL_overlap, CPL_get_olap_limits, CPL_my_proc_portion, CPL_get_no_cells

        logical, intent(out) :: no_error
        double precision, dimension(:,:,:,:), intent(in) :: recv_array
        integer :: i,j,k,ii,jj,kk
        integer, dimension(3) :: Ncells
        integer, dimension(6) :: limits, portion

        no_error = .false.

        call CPL_get_olap_limits(limits)
        call CPL_my_proc_portion(limits, portion)
        call CPL_get_no_cells(portion, Ncells)

        ! Check that every processor inside the overlap region receives correctly the cell
        ! number.  
        if (CPL_overlap()) then
            no_error = .true.
            do i = 1, Ncells(1)
            do j = 1, Ncells(2)
            do k = 1, Ncells(3)
                ! -2 indices to match c++ and python indexing in portion and i,j,k
                ii = i + portion(1) - 2
                jj = j + portion(3) - 2
                kk = k + portion(5) - 2

                if ((dble(ii) - recv_array(1,i,j,k)) .gt. 1e-8) then 
                    print'(a,i5,a,i6,a,f10.5)', "ERROR -- portion in x: ", portion(1:2), & 
                           " cell i: ",ii, & 
                           " recv_array: ", recv_array(1,i,j,k)
                    no_error = .false.
                endif
                if ((dble(jj) - recv_array(2,i,j,k)) .gt. 1e-8) then 
                    print'(a,i5,a,i6,a,f10.5)', "ERROR -- portion in y: ", portion(3:4), & 
                          " cell j: ", jj , & 
                           " recv_array: ", recv_array(2,i,j,k)
                    no_error = .false.  
                endif
                if ((dble(kk) - recv_array(3,i,j,k)) .gt. 1e-8) then 
                    print'(a,i5,a,i6,a,f10.5)', "ERROR -- portion in z: ", portion(5:6), & 
                           " cell k: ", kk , & 
                           " recv_array: ", recv_array(3,i,j,k)
                    no_error = .false.
                endif
            enddo
            enddo
            enddo
        endif

    end subroutine check_array

    subroutine milisleep(dt)

        ! desired sleep interval
        double precision, intent(in) :: dt                
        integer,dimension(8) :: t ! arguments for date_and_time
        integer :: ms1,ms2        ! start and end times [ms]

        ! Get start time:
        call date_and_time(values=t)
        ms1=(t(5)*3600+t(6)*60+t(7))*1000+t(8)

        do ! check time:
          call date_and_time(values=t)
          ms2=(t(5)*3600+t(6)*60+t(7))*1000+t(8)
          if(ms2-ms1>=dt) exit
        enddo

    end subroutine milisleep



!--------------------------------------------------------------------------------------
! Split an integer into three factors minimising value of each

subroutine find3factors(n,nx,ny,nz)
	implicit none

	integer, intent(in)		:: n
	integer, intent(out)	:: nx,ny,nz

	integer :: nfactors,  minval1, minval2
	integer,dimension(1)	:: minlocation
	integer, allocatable, dimension(:) :: factors, nonzerof
	integer,parameter :: large=99999999

	!Check for sensible configuration
	nx = ceiling( dble(n)**(1.d0/3.d0))
	ny = ceiling((dble(n)/dble(nx))**(1.d0/2.d0))
	nz = ceiling( dble(n)/(dble(nx)*dble(ny)))
	if (nx * ny * nz .eq. n) return

	!Otherwise find it out the hard way
	if (n .ge. 4) then
		! Find all prime factors
		allocate(factors(n/2))
		factors = 0
		call primefactors(n,factors,nfactors)
	else
		! First 3 numbers are primes
		allocate(factors(1))
		factors = n
		nfactors = 1
	endif

	! Reduce/increase number of factors to three
	if (nfactors .eq. 1) then
		nx = factors(1); ny = 1         ; nz = 1
	elseif (nfactors .eq. 2) then
		nx = factors(1); ny = factors(2); nz = 1
	elseif (nfactors .eq. 3) then
		nx = factors(1); ny = factors(2); nz = factors(3)		
	elseif (nfactors .gt. 3) then
		allocate(nonzerof(nfactors))
		nonzerof = factors(1:nfactors)

		do while (nfactors .gt. 3)
			!Multiple two minimum values and store in 
			minlocation = minloc(nonzerof)
			minval1 = nonzerof(minlocation(1))
			nonzerof(minlocation(1)) = LARGE
			minlocation = minloc(nonzerof)
			minval2 = nonzerof(minlocation(1))
			nonzerof(minlocation(1)) = LARGE

			nonzerof(minlocation(1))=minval1*minval2
			nfactors = nfactors - 1
		enddo

		minlocation = minloc(nonzerof)
		nz = nonzerof(minlocation(1)); nonzerof(minlocation(1))=LARGE
		minlocation = minloc(nonzerof)
		ny = nonzerof(minlocation(1)); nonzerof(minlocation(1))=LARGE
		minlocation = minloc(nonzerof)
		nx = nonzerof(minlocation(1)); nonzerof(minlocation(1))=LARGE
 		
	endif

	if (n - nx*ny*nz .ne. 0) stop "ERROR in find3factors"

end subroutine find3factors

!--------------------------------------------------------------------------------------
!		SUBROUTINE TO FIND THE PRIME FACTORS OF A NUMBER
! Start with 2, check whether 2 is a factor by seeing if MOD(<input_number>,2)
! is zero. If it is zero, then 2 becomes a factor. If not, check with the next number.
! When a factor is found, divide the given number with the factor found. However,
! do not move to the next possible factor - a number can occur more than once as a factor

subroutine primefactors(num, factors, f)
	implicit none

	integer, intent(in) :: num  !input number
	integer,intent(out), dimension((num/2))::factors !array to store factors
	integer, intent(inout) :: f
	integer :: i, n

	i = 2  !eligible factor
	f = 1  !number of factors
	n = num !store input number into a temporary variable
	do
		if (mod(n,i) == 0) then !if i divides 2, it is a factor
			factors(f) = i
			f = f+1
			n = n/i
		else
			i = i+1     !not a factor. move to next number
		end if
		if (n == 1) then		
			f = f-1		!its value will be one more than the number of factors
			exit
		end if
	end do

end subroutine primefactors

end module stuff

