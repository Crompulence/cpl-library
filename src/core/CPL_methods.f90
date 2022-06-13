!=============================================================================
!
!    ________/\\\\\\\\\__/\\\\\\\\\\\\\____/\\\_____________
!     _____/\\\////////__\/\\\/////////\\\_\/\\\_____________
!      ___/\\\/___________\/\\\_______\/\\\_\/\\\_____________
!       __/\\\_____________\/\\\\\\\\\\\\\/__\/\\\_____________
!        _\/\\\_____________\/\\\/////////____\/\\\_____________
!         _\//\\\____________\/\\\_____________\/\\\_____________
!          __\///\\\__________\/\\\_____________\/\\\_____________
!           ____\////\\\\\\\\\_\/\\\_____________\/\\\\\\\\\\\\\\\_
!            _______\/////////__\///______________\///////////////__
!
!
!                         C P L  -  L I B R A R Y
!
!           Copyright (C) 2012-2022 Edward Smith
!
!=============================================================================

module coupler
!
!Routines accessible from application ( molecular or continuum ) after
!the name, in parenthesis, is the realm in which each routine must be called
!
!  - CPL_send            (cfd+md)   sends grid data exchanged between realms ( generic interface)
!  - CPL_recv            (cfd+md)   receives data exchanged between realms ( generic interface)
!  - CPL_cfd_get         (cfd)    returns coupler internal parameters for CFD realm
!  - CPL_get             (md)    returns coupler internal parameters
!  - CPL_md_get_save_period     (md)    auxiliary used for testing [depricated]
!  - CPL_md_get_average_period  (md)    returns average period of BC [depricated]
!  - CPL_md_get_md_per_cfd_dt   (md)    returns the number of step MD does for each CFD step [depricated]
!
!Also see `coupler module <fortran_api.html#coupler-module>`_ which contains all
!routines to setup the simulation and all the variables defining the setup topology.
!
!.. code-block:: guess
!
!    ________/\\\\\\\\\__/\\\\\\\\\\\\\____/\\\_____________
!     _____/\\\////////__\/\\\/////////\\\_\/\\\_____________
!      ___/\\\/___________\/\\\_______\/\\\_\/\\\_____________
!       __/\\\_____________\/\\\\\\\\\\\\\/__\/\\\_____________
!        _\/\\\_____________\/\\\/////////____\/\\\_____________
!         _\//\\\____________\/\\\_____________\/\\\_____________
!          __\///\\\__________\/\\\_____________\/\\\_____________
!           ____\////\\\\\\\\\_\/\\\_____________\/\\\\\\\\\\\\\\\_
!            _______\/////////__\///______________\///////////////__
!
!
!                         C P L  -  L I B R A R Y
!
!           Copyright (C) 2012-2022 Edward Smith
!
!**Author(s)**
!
! - Edward Smith Novemeber 2011 to present
! - Eduardo Ramos Fernandez 2015 to 2018
! - David Trevelyan September 2012 to December 2015
! - Lucian Anton, November 2011  
!
!**License**
!
!    This file is part of CPL-Library.
!
!    CPL-Library is free software: you can redistribute it and/or modify
!    it under the terms of the GNU General Public License as published by
!    the Free Software Foundation, either version 3 of the License, or
!    (at your option) any later version.
!
!    CPL-Library is distributed in the hope that it will be useful,
!    but WITHOUT ANY WARRANTY; without even the implied warranty of
!    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
!    GNU General Public License for more details.
!
!    You should have received a copy of the GNU General Public License
!    along with CPL-Library.  If not, see <http://www.gnu.org/licenses/>.
!
    implicit none

    private CPL_send_full, CPL_send_min, CPL_recv_full, CPL_recv_min

    public CPL_gather, CPL_scatter, CPL_get, CPL_proc_extents,&
           CPL_my_proc_extents, CPL_proc_portion, CPL_my_proc_portion, &
           CPL_map_cell2coord, CPL_map_coord2cell, CPL_get_no_cells, &
           CPL_map_glob2loc_cell, CPL_get_olap_limits, CPL_get_cnst_limits, &
           CPL_get_bnry_limits, CPL_map_cfd2md_coord, CPL_map_md2cfd_coord, & 
           CPL_overlap, CPL_realm, CPL_cart_coords, CPL_cfd_dt, CPL_md2cfd, &
           CPL_cfd2md, CPL_is_proc_inside, CPL_swaphalos

    interface CPL_send
		!
		!Interface to allow single call to CPL send to be used for different inputs
		!
        module procedure CPL_send_full, CPL_send_min
    end interface CPL_send

    interface CPL_recv
		!
		!Interface to allow single call to CPL recv to be used for different inputs
		!
        module procedure CPL_recv_full, CPL_recv_min
    end interface CPL_recv


contains


!=============================================================================
!                _____ _                 _       _   _
!               /  ___(_)               | |     | | (_)
!               \ `--. _ _ __ ___  _   _| | __ _| |_ _  ___  _ __
!                `--. \ | '_ ` _ \| | | | |/ _` | __| |/ _ \| '_ \
!               /\__/ / | | | | | | |_| | | (_| | |_| | (_) | | | |
!               \____/|_|_| |_| |_|\__,_|_|\__,_|\__|_|\___/|_| |_|
!
!------------------------------------------------------------------------------
!                              CPL_gather                                     -
!------------------------------------------------------------------------------
subroutine CPL_gather(gatherarray, npercell, limits, recvarray)
!
!Perform gather operation on CPL_OLAP_COMM communicator. The CFD processor
!is the root process. The gathered data is effectively "slotted" into the
!correct part of the recvarray, and is intented for use in providing the
!CFD simulation boundary conditions with data obtained from the MD
!simulation.
!
!**Synopsis**
!
!  - CPL_gather(gatherarray,npercell,limits,recvarray)
!
!**Inputs**
!
! - gatherarray
!   - Assumed shape array of data to be gathered from each MD processor
!     in the overlap communicator. Must be size 0 on CFD processor.
!
! - limits
!   - Integer array of length 6, specifying the global cell extents of the
!     region to be gathered, is the same on ALL processors.
!
! - npercell
!   - number of data points per cell to be gathered (integer)
!     Note: should be the same as size(gatherarray(1)) for MD
!     processor. E.G. npercell = 3 for gathering 3D velocities.
!
!**Inputs/Output**
!
! - recvarray
!   - The array in which the gathered values are to be stored on the CFD
!     processor. The only values to be changed in recvarray are:
!     recvarray(limits(1):limits(2),limits(3):limits(4),limits(5):limits(6))
!
!**Output**
!  - NONE
!
! .. sectionauthor:: David Trevelyan
!
    use mpi
    use coupler_module
    implicit none

    integer, intent(in) :: npercell
    integer, intent(in) :: limits(6)
    real(kind(0.d0)), dimension(:,:,:,:), intent(in)    :: gatherarray
    real(kind(0.d0)), dimension(:,:,:,:), intent(inout) :: recvarray

    integer :: sendcount
    integer, dimension(:), allocatable :: recvcounts, displs
    real(kind(0.d0)), dimension(:), allocatable :: sendbuf
    real(kind(0.d0)), dimension(:), allocatable :: recvbuf
   
    if (.not. CPL_overlap()) return

    call check_limits_consistency

    call prepare_gatherv_parameters
    
    if (realm .eq. md_realm) call pack_sendbuf
    
    call MPI_gatherv(sendbuf, sendcount, MPI_DOUBLE_PRECISION, recvbuf, &
                     recvcounts, displs, MPI_DOUBLE_PRECISION, CFDid_olap, &
                     CPL_OLAP_COMM, ierr)

!    if (realm .eq.md_realm) then
!        print'(a,5i8,3f18.12)', "md gather", sendcount, shape(gatherarray), sum(sendbuf), minval(sendbuf), maxval(sendbuf)
!    elseif( realm .eq. cfd_realm) then
!        do sendcount=1,size(recvcounts)
!            print'(a,6i8,3f18.12)', "cfd gather", sendcount, recvcounts(sendcount), & 
!                                    shape(recvarray), sum(recvbuf), minval(recvbuf), maxval(recvbuf)
!        enddo
!    endif


    if (realm .eq. cfd_realm) call unpack_recvbuf

    call deallocate_gather_u
    
contains

    subroutine check_limits_consistency
        implicit none

        logical :: consistent = .true.

        integer :: i, cfd_limits(6), md_limits(6)
        integer :: md_limits_send(6), cfd_limits_send(6)

        select case (realm)
        case (md_realm)

            md_limits = limits
            cfd_limits = 0

        case (cfd_realm)

            md_limits = 0
            cfd_limits = limits

        end select

        md_limits_send = md_limits
        cfd_limits_send = cfd_limits

        ! Since we set md/cfd limits to zero in opposing realms, MPI_MAX
        ! will result in the correct limits on both sides
        call MPI_Allreduce(md_limits_send, md_limits, 6, MPI_INTEGER, MPI_MAX, &
                           CPL_OLAP_COMM, ierr)
        call MPI_Allreduce(cfd_limits_send, cfd_limits, 6, MPI_INTEGER, MPI_MAX, &
                           CPL_OLAP_COMM, ierr)

        do i = 1, 6
            if (md_limits(i) .ne. cfd_limits(i)) consistent = .false.
        end do
        
        if (.not. consistent) then
            print'(2(a,6i8))', 'MD limits=', md_limits, ' CFD limits=', cfd_limits
            call error_abort("CPL_gather error - MD and CFD limits not consistent in CPL_gather")
        else
            return
        end if

    end subroutine check_limits_consistency

    subroutine prepare_gatherv_parameters
        implicit none

        integer :: coord(3),portion(6)
        integer :: ncells,bufsize
        integer :: trank_olap,tid_olap
        
        ! Check send limits are inside overlap region
        if (limits(1) .lt. icmin .or. &
            limits(2) .gt. icmax .or. &
            limits(3) .lt. jcmin .or. &
            limits(4) .gt. jcmax .or. &
            limits(5) .lt. kcmin .or. &
            limits(6) .gt. kcmax) then

            if (limits(1) .lt. icmin) then
                print'(2(a,6i8))', "CPL_gather x minimum limit = ", limits(1), & 
                                   " is less than domain limit ", icmin
            endif
            if (limits(2) .gt. icmax) then
                print'(2(a,6i8))', "CPL_gather x maximum limit = ", limits(2), &
                                   " is greater than domain limit ", icmax
            endif
            if (limits(3) .lt. jcmin) then
                print'(2(a,6i8))', "CPL_gather y minimum limit = ", limits(3),  &
                                   " is less than domain limit ", jcmin
            endif
            if (limits(4) .gt. jcmax) then
                print'(2(a,6i8))', "CPL_gather y maximum limit = ", limits(4),  &
                                   " is greater than domain limit ", jcmax
            endif
            if (limits(5) .lt. kcmin) then
                print'(2(a,6i8))', "CPL_gather z minimum limit = ", limits(5),  &
                                   " is less than domain limit ", kcmin
            endif
            if (limits(6) .gt. kcmax) then
                print'(2(a,6i8))', "CPL_gather z maximum limit = ", limits(6),  &
                                   " is greater than domain limit ", kcmax
            endif
            
        call error_abort("CPL_gather error - Gather limits are outside global domain. " // &
                             "Aborting simulation.")
            
        end if
        

        ! Check if CFD processor has tried to "send" anything
        if (myid_olap.eq.CFDid_olap .and. product(shape(gatherarray)).ne.0) then
        call error_abort('CPL_gather error - CFD proc input to CPL_gather: '          // &
                             'gatherarray has nonzero size. Aborting ' // &
                             'from prepare_gatherv_parameters' )
        end if

        ! Allocate send buffer
        if (realm.eq.md_realm) then
            call CPL_Cart_coords(CPL_CART_COMM,rank_cart,realm,3,coord,ierr)
            call CPL_proc_portion(coord,md_realm,limits,portion,ncells)
            bufsize = npercell*ncells
        else
            bufsize = 0
        end if
        allocate(sendbuf(bufsize))

        ! Allocate array of sendbuffer sizes and populate it
        allocate(recvcounts(nproc_olap))
        do trank_olap = 1,nproc_olap
            tid_olap = trank_olap - 1
            if (tid_olap .eq. CFDid_olap) then
                recvcounts(trank_olap) = 0
            else
                call CPL_Cart_coords(CPL_OLAP_COMM,trank_olap,md_realm,3, &
                                     coord,ierr)
                call CPL_proc_portion(coord,md_realm,limits,portion,ncells)
                recvcounts(trank_olap) = npercell*ncells
            end if
        end do
    
        ! Own sendbuffer size
        sendcount = size(sendbuf)
        ! Sanity check
        if (sendcount .ne. recvcounts(rank_olap)) then
            call error_abort('Send buffer sizes calculated incorrectly '//&
                             'in prepare_gatherv_parameters. Aborting.')
        end if

        ! Allocate recvbuffer on CFD proc
        allocate(recvbuf(sum(recvcounts)))

        ! Calculate displacements for each proc in array recvbuffer
        allocate(displs(nproc_olap))
        displs(1) = 0
        do trank_olap=2,nproc_olap
            displs(trank_olap) = sum(recvcounts(1:trank_olap-1))    
        end do

    end subroutine prepare_gatherv_parameters

    subroutine pack_sendbuf
        implicit none

        integer :: pos, ixyz, icell, jcell, kcell
        integer :: i,j,k
        integer :: coord(3),portion(6)

        ! Get MD processor extents and cells portion of send region
        call CPL_Cart_coords(CPL_CART_COMM,rank_cart,realm,3,coord,ierr) 
        call CPL_proc_portion(coord,realm,limits,portion)

        ! If MD proc has nothing to send, exit
        if (any(portion.eq.VOID)) return

        pos = 1
        do ixyz  = 1,npercell 
        do icell = portion(1),portion(2) 
        do jcell = portion(3),portion(4) 
        do kcell = portion(5),portion(6) 
            ! Map to local coords (assumed shape array has
            ! lower bound 1 by default when input to subroutine) 
            i = icell - portion(1) + 1
            j = jcell - portion(3) + 1
            k = kcell - portion(5) + 1
            sendbuf(pos) = gatherarray(ixyz,i,j,k)
            pos = pos + 1
        end do
        end do
        end do
        end do
    
    end subroutine pack_sendbuf
    
    subroutine unpack_recvbuf
        implicit none

        integer :: coord(3),md_portion(6),cfd_portion(6)
        integer :: trank_olap, tid_olap
        integer :: pos,ixyz,icell,jcell,kcell
        integer :: i,j,k

        ! Get CFD proc coords and extents, allocate suitable array
        call CPL_cart_coords(CPL_OLAP_COMM,rank_olap,cfd_realm,3,coord,ierr)
        call CPL_proc_portion(coord,cfd_realm,limits, cfd_portion)

        ! Loop over all processors in overlap comm
        do trank_olap = 1,nproc_olap

            tid_olap = trank_olap - 1
            if (tid_olap .eq. CFDid_olap) cycle

            call CPL_Cart_coords(CPL_OLAP_COMM,trank_olap,md_realm,3,coord,ierr)
            call CPL_proc_portion(coord,md_realm,limits,md_portion)

            if (any(md_portion.eq.VOID)) cycle

            ! Set position and unpack MD proc's part of recvbuf to
            ! correct region of recvarray   
            pos = displs(trank_olap) + 1    
            do ixyz = 1,npercell
            do icell = md_portion(1),md_portion(2)
            do jcell = md_portion(3),md_portion(4)
            do kcell = md_portion(5),md_portion(6)

                i = icell - cfd_portion(1) + 1 
                j = jcell - cfd_portion(3) + 1 
                k = kcell - cfd_portion(5) + 1 

                recvarray(ixyz,i,j,k) = recvbuf(pos)
                pos = pos + 1

!               write(8000+myid_world,'(a,i4,a,i4,a,i4,a,i4,a,f20.1)'),   &
!                     'recvarray(',ixyz,',',icell,',',jcell,',',kcell,') =', &
!                      recvarray(ixyz,i,j,k)

            end do  
            end do  
            end do
            end do
                    
        end do

    end subroutine unpack_recvbuf
    
    subroutine deallocate_gather_u
        implicit none

        if(allocated(recvcounts)) deallocate(recvcounts)
        if(allocated(displs))     deallocate(displs)
        if(allocated(sendbuf))    deallocate(sendbuf)
        if(allocated(recvbuf))    deallocate(recvbuf)

    end subroutine deallocate_gather_u

end subroutine CPL_gather

!------------------------------------------------------------------------------
!                              CPL_scatter                                    -
!------------------------------------------------------------------------------
subroutine CPL_scatter(scatterarray,npercell,limits,recvarray)
!
!Scatter cell-wise data from CFD processor to corresponding MD processors
!on the overlap communicator CPL_OLAP_COMM.
!
!**Synopsis**
!
!  - CPL_scatter(scatterarray,npercell,limits,recvarray)
!
!**Inputs**
!
!  - scatterarray
!   - assumed shape array of data to be scattered (real(kind(0.d0)))
!
!  - limits
!   - integer array of length 6, specifying the global cell extents of the
!     region to be scattered, is the same on all processors.
!
!  - npercell
!   - number of data points per cell to be scattered (integer).
!     Note: should be the same as size(scatterarray(1)) for CFD proc
!
!**Inputs/Output**
!  - recvarray
!   - the array in which the scattered values are stored on the MD
!     processors.
!
!**Output**
!  - NONE
!
! .. sectionauthor:: David Trevelyan
!
    use coupler_module
    use mpi
    implicit none

    integer,intent(in) :: npercell
    integer,intent(in) :: limits(6)
    real(kind(0.d0)),dimension(:,:,:,:),intent(in) :: scatterarray
    real(kind(0.d0)),dimension(:,:,:,:),intent(inout) :: recvarray

    integer :: recvcount
    integer, dimension(:), allocatable :: displs, sendcounts
    real(kind(0.d0)), dimension(:), allocatable :: recvbuf
    real(kind(0.d0)), dimension(:), allocatable :: scatterbuf

    if (.not. CPL_overlap()) return

    call prepare_scatterv_parameters

    if (realm .eq. cfd_realm) call pack_scatterbuf

    call MPI_scatterv(scatterbuf, sendcounts, displs, MPI_DOUBLE_PRECISION, &
                      recvbuf, recvcount, MPI_DOUBLE_PRECISION, CFDid_olap, &
                      CPL_OLAP_COMM, ierr)

!    if (realm .eq.md_realm) then
!        print'(a,3f18.12)', "md scatter", sum(recvbuf),minval(recvbuf),maxval(recvbuf)
!    elseif( realm .eq. cfd_realm) then
!        !print'(a,3f18.12)', "cfd scatter", sum(scatterbuf),minval(scatterbuf),maxval(scatterbuf)
!    endif


    if (realm .eq. md_realm)  call unpack_scatterbuf

    call deallocate_scatter_s

contains

    subroutine prepare_scatterv_parameters
        implicit none

        integer :: ncells
        integer :: coord(3),portion(6)
        integer :: bufsize
        integer :: trank_olap, tid_olap

        ! Check send limits are inside overlap region
        if (limits(1) .lt. icmin .or. &
            limits(2) .gt. icmax .or. &
            limits(3) .lt. jcmin .or. &
            limits(4) .gt. jcmax .or. &
            limits(5) .lt. kcmin .or. &
            limits(6) .gt. kcmax) then

            if (limits(1) .lt. icmin) then
                print'(2(a,6i8))', "CPL_scatter x minimum limit = ", limits(1), " is less than domain limit ", icmin
            endif
            if (limits(2) .gt. icmax) then
                print'(2(a,6i8))', "CPL_scatter x maximum limit = ", limits(2), " is greater than domain limit ", icmax
            endif
            if (limits(3) .lt. jcmin) then
                print'(2(a,6i8))', "CPL_scatter y minimum limit = ", limits(3), " is less than domain limit ", jcmin
            endif
            if (limits(4) .gt. jcmax) then
                print'(2(a,6i8))', "CPL_scatter y maximum limit = ", limits(4), " is greater than domain limit ", jcmax
            endif
            if (limits(5) .lt. kcmin) then
                print'(2(a,6i8))', "CPL_scatter z minimum limit = ", limits(5), " is less than domain limit ", kcmin
            endif
            if (limits(6) .gt. kcmax) then
                print'(2(a,6i8))', "CPL_scatter z maximum limit = ", limits(6), " is greater than domain limit ", kcmax
            endif
            
        call error_abort("CPL_scatter error - Scatter limits are outside global domain. " // &
                             "Aborting simulation.")
            
        end if

        if (realm.eq.cfd_realm) then
            call CPL_Cart_coords(CPL_CART_COMM,rank_cart,cfd_realm,3,coord,ierr)
            call CPL_proc_portion(coord,cfd_realm,limits,portion,ncells)
            bufsize = npercell*ncells
        else
            bufsize = 0
        end if

        allocate(scatterbuf(bufsize))
        allocate(sendcounts(nproc_olap))
        allocate(displs(nproc_olap))

        ! Loop over all procs in overlap comm
        do trank_olap = 1,nproc_olap
            tid_olap = trank_olap - 1
            ! Calculate number of data points to scatter to each proc
            if (tid_olap .eq. CFDid_olap) then
                sendcounts(trank_olap) = 0
            else
                call CPL_Cart_coords(CPL_OLAP_COMM,trank_olap,md_realm,3, &
                                     coord, ierr)
                call CPL_proc_portion(coord,md_realm,limits,portion,ncells)
                sendcounts(trank_olap) = npercell*ncells
            end if
        end do

        ! Get number of data points this MD proc will receive
        recvcount = sendcounts(rank_olap)

        ! Calculate starting positions of each MD proc region in
        ! scatterbuf array
        displs(1) = 0
        do trank_olap=2,nproc_olap
            displs(trank_olap) = sum(sendcounts(1:trank_olap-1))
        end do

        ! Allocate space to receive data
        allocate(recvbuf(sum(sendcounts)))

    end subroutine prepare_scatterv_parameters

    subroutine pack_scatterbuf
        implicit none
        
        integer :: pos, n
        integer :: tid_olap
        integer :: cfd_coord(3), md_coord(3),cfdextents(6),portion(6)
        integer :: ixyz, icell, jcell, kcell
        integer :: i,j,k

        ! Grab CFD proc extents to be used for local cell mapping (when an 
        ! array is an input to subroutine it has lower bound 1 by default)
        call CPL_cart_coords(CPL_CART_COMM,rank_cart,cfd_realm,3,cfd_coord,ierr)
        !call CPL_proc_extents(cfd_coord,cfd_realm,cfdextents)
        call CPL_proc_portion(cfd_coord,cfd_realm,limits,cfdextents)

        ! Loop over procs in olap comm and pack scatter buffer 
        ! in separate regions for each MD proc
        pos = 1
        do n = 1,nproc_olap

            tid_olap = n - 1
            if (tid_olap.eq.CFDid_olap) cycle

            call CPL_Cart_coords(CPL_OLAP_COMM,n,md_realm,3,md_coord,ierr)
            call CPL_proc_portion(md_coord,md_realm,limits,portion)
            if (any(portion.eq.VOID)) cycle

            do ixyz = 1,npercell
            do icell= portion(1),portion(2)
            do jcell= portion(3),portion(4)
            do kcell= portion(5),portion(6)
                i = icell - cfdextents(1) + 1
                j = jcell - cfdextents(3) + 1
                !print*, "CPL: ", j, jcell, "|",cfdextents, "Coord: ", cfd_coord
                k = kcell - cfdextents(5) + 1
                scatterbuf(pos) = scatterarray(ixyz,i,j,k)
                pos = pos + 1
            end do
            end do
            end do
            end do
            
        end do
    
    end subroutine pack_scatterbuf  

    subroutine unpack_scatterbuf
        implicit none

        integer :: pos, ierr
        integer :: coord(3),portion(6),extents(6)
        integer :: ixyz,icell,jcell,kcell
        integer :: i,j,k

        call CPL_cart_coords(CPL_OLAP_COMM,rank_olap,md_realm,3,coord,ierr)
        call CPL_proc_portion(coord,realm,limits, extents)
        call CPL_proc_portion(coord,realm,limits,portion)
        if (any(portion.eq.VOID)) return

        pos = 1
        do ixyz = 1,npercell
        do icell= portion(1),portion(2)
        do jcell= portion(3),portion(4)
        do kcell= portion(5),portion(6)
            i = icell - portion(1) + 1
            j = jcell - portion(3) + 1
            k = kcell - portion(5) + 1
            recvarray(ixyz,i,j,k) = recvbuf(pos)
!           write(7000+myid_realm,'(i4,a,i4,a,i4,a,i4,a,i4,a,f20.1)'),        &
!                 rank_cart,' recvarray(',ixyz,',',icell,',',jcell,',',kcell, &
!                 ') =',recvarray(ixyz,i,j,k)
            pos = pos + 1
        end do  
        end do  
        end do
        end do

    end subroutine unpack_scatterbuf
    
    subroutine deallocate_scatter_s
        implicit none

        if(allocated(displs))     deallocate(displs)
        if(allocated(scatterbuf)) deallocate(scatterbuf)
        if(allocated(recvbuf))    deallocate(recvbuf)
        if(allocated(sendcounts)) deallocate(sendcounts)

    end subroutine deallocate_scatter_s

end subroutine CPL_scatter

! ----------------------------------------------------------------------------
subroutine CPL_send_full(asend, limits, send_flag)
!
!Send four dimensional array *asend* of data from all processors in the 
!current realm with data between global cell array *limits* to the 
!corresponding processors from the other realm.
!
!**Remarks**
!
!Assumes the coupler has been initialised with `CPL_init <#f/_/cpl_init>`_ and 
!topological mapping has been setup using either `CPL_setup_md <#f/_/cpl_setup_md>`_ 
!or `CPL_setup_cfd <#f/_/cpl_setup_cfd>`_ as appropriate.
!
!**Synopsis**
!
!.. code-block:: fortran
!
!  CPL_send(asend, limits, send_flag)    
!
!**Inputs**
!
! - asend
!
!   - Array of data to send. Should be a four dimensional array allocated using the number of cells on the current processor between the limits. Size should be be obtained from `CPL_my_proc_portion(limits, portion) <#f/coupler/cpl_my_proc_portion>`_.
! - limits 
!
!   - Limits in global cell coordinates, must be the same as corresponding recieve
!
!**Outputs**
!
! - send_flag
!
!   - Returned flag which indicates success or failure of send process
!
!**Example**
!
!.. code-block:: fortran
!
!  call CPL_get_olap_limits(limits)
!  call CPL_my_proc_portion(limits, portion)
!  call CPL_get_no_cells(portion, Ncells)
!
!  !Coupled Send array
!  allocate(A(3, Ncells(1), Ncells(2), Ncells(3)))
!
!  do i =portion(1),portion(2)
!  do j =portion(3),portion(4)
!  do k =portion(5),portion(6)
!     ii = i-portion(1)+1; jj = j-portion(3)+1; kk = k-portion(5)+1
!     A(1,ii,jj,kk) = dble(i)
!     A(2,ii,jj,kk) = dble(j)
!     A(3,ii,jj,kk) = dble(k)
!  enddo
!  enddo
!  enddo
!  call CPL_send(A, limits, send_flag)
!
! .. sectionauthor::Edward Smith
! ----------------------------------------------------------------------------
    use mpi
    use coupler_module, only : md_realm, cfd_realm, error_abort, & 
                               CPL_GRAPH_COMM, myid_graph,olap_mask, &
                               rank_world, realm, ierr, VOID, &
                               CPL_setup_complete, REALM_NAME, &
                               myid_realm, rootid_realm, & 
                               rank_intersect, CPL_REALM_INTERSECTION_COMM
    implicit none

    
    logical, intent(out), optional  :: send_flag !Flag set if processor has passed data   
    integer, dimension(6), intent(in) :: limits ! Global cell indices with minimum and maximum values to send
    real(kind=kind(0.d0)),dimension(:,:,:,:), intent(in):: asend ! Array containing data to send
   
    !Neighbours
    integer                             :: nneighbors   
    integer,dimension(:),allocatable    :: id_neighbors, req
    integer,dimension(:,:),allocatable :: status

    ! local indices 
    integer :: icell,jcell,kcell
    integer :: n,pos,iclmin,iclmax,jclmin,jclmax,kclmin,kclmax
    integer :: npercell

    ! auxiliaries 
    integer                             :: nbr, ndata, itag, destid, Ncells
    integer,dimension(3)                :: pcoords, Ncell
    integer,dimension(6)                :: portion, myportion
    real(kind=kind(0.d0)), allocatable  :: vbuf(:)

    !Check counter
    integer       :: source, check_send
    integer, save :: first_counter = 4, ncalls=0

    !Check setup is complete
    if (CPL_setup_complete .ne. 1) then
       call error_abort("Error CPL_send called before CPL_setup_"//REALM_NAME(realm))
    endif

    ! This local CFD domain is outside MD overlap zone 
    if (olap_mask(rank_world) .eqv. .false.) return

    ! Number of components at each grid point
    npercell = size(asend,1)

    !Get current processors portion
    call CPL_my_proc_portion(limits, myportion)

    !Check size is consistent with array
    if (myportion(2)-myportion(1)+1 .ne. size(asend,2)) then
        print*, myportion(2)-myportion(1)+1, size(asend,2)
        call error_abort("Error in CPL send -- x size of asend must be the same as portion(2)-portion(1)+1")
    endif
    if (myportion(4)-myportion(3)+1 .ne. size(asend,3)) then
        print*, myportion(4)-myportion(3)+1, size(asend,3)
        call error_abort("Error in CPL send -- y size of asend must be the same as portion(4)-portion(3)+1")
    endif
    if (myportion(6)-myportion(5)+1 .ne. size(asend,4)) then
        print*, myportion(6)-myportion(5)+1, size(asend,4)
        call error_abort("Error in CPL send -- z size of asend must be the same as portion(6)-portion(5)+1")
    endif

    !Check send/recv are consistent size
    if (first_counter .ne. 0) then
        first_counter = first_counter - 1
        check_send = 0
        call MPI_AllReduce(npercell, check_send, 1, MPI_INTEGER, & 
                           MPI_SUM, CPL_REALM_INTERSECTION_COMM, ierr)
        call MPI_comm_size(CPL_REALM_INTERSECTION_COMM, n, ierr)
        !print*, "Send npercell=", npercell, check_send, n
        if (check_send/n .ne. npercell) then
            print*, "Error in CPL send in ", realm_name(realm), & 
                    " realm, expected send size: ", npercell, " but recv is expecting ", check_send/n
            call error_abort("Error in CPL send -- first index of send array does not match recv data size")
        endif
    endif

    !Get neighbours
    call MPI_Graph_neighbors_count(CPL_GRAPH_COMM,myid_graph,nneighbors,ierr)
    allocate(id_neighbors(nneighbors))
    call MPI_Graph_neighbors(CPL_GRAPH_COMM,myid_graph,nneighbors,id_neighbors,ierr )

    !Set sendflag to false and only change if anything is sent
    if (present(send_flag)) send_flag = .false.

    !Keep track of number of sends
    ncalls = ncalls + 1



    ! Send to all attached processors
    allocate(req(nneighbors)); req = MPI_REQUEST_NULL
    allocate(status(MPI_STATUS_SIZE, nneighbors));   
    ! loop through the maps and send the corresponding sections of asend
    do nbr = 1, nneighbors

        !Get taget processor from mapping
        destid = id_neighbors(nbr)

        ! ----------------- pack data for destid-----------------------------
        if (realm .eq. cfd_realm) then
            !Get portion of nbr MD processor to send to
            call CPL_Cart_coords(CPL_GRAPH_COMM, destid+1,  md_realm, 3, pcoords, ierr)
            call CPL_proc_portion(pcoords, md_realm, limits, portion, Ncells)
        elseif (realm .eq. md_realm) then
            !Data to send is based on current processor as MD proc size < CFD proc size
            call CPL_my_proc_portion(limits, portion)
            call CPL_get_no_cells(portion, Ncell)
            Ncells = product(Ncell)
        endif

        !Get index in local cell coordinates for data
        iclmin = portion(1)-myportion(1)+1;   iclmax = portion(2)-myportion(1)+1
        jclmin = portion(3)-myportion(3)+1;   jclmax = portion(4)-myportion(3)+1
        kclmin = portion(5)-myportion(5)+1;   kclmax = portion(6)-myportion(5)+1

        ! Amount of data to be sent
        if (any(portion.eq.VOID)) then
            !print*, 'VOID send qqqq',realm_name(realm),rank_world,rank_realm
            ndata = 0
        else
            ndata = npercell * Ncells
            if (allocated(vbuf)) deallocate(vbuf); allocate(vbuf(ndata))
            if (present(send_flag)) send_flag = .true.

            ! Pack array into buffer
            pos = 1
!            print'(a,7i8,i4,21i3)', 'send qqqq', rank_world, rank_realm,            & 
!                                   rank_olap, ndata, nbr, destid, size(asend), pos, &
!                                   iclmin, iclmax, jclmin, jclmax, kclmin, kclmax,  &
!                                   limits, portion
            do kcell=kclmin,kclmax
            do jcell=jclmin,jclmax
            do icell=iclmin,iclmax
            do n = 1,npercell
                vbuf(pos) = asend(n, icell, jcell, kcell)
                !write(98000+destid+1+10*rank_world,'(3i8,f20.5)') rank_world,destid+1,n, vbuf(pos)
                pos = pos + 1
            end do
            end do
            end do
            end do
            ! ----------------- pack data for destid -----------------------------

            ! Send data 
            !Using MPI_TAG_UB directly causes error with OpenMPI
            !call MPI_COMM_GET_ATTR(CPL_GRAPH_COMM, MPI_TAG_UB, MPI_TAG_UB_int, ATTR_FLAG, ierr)
            !itag = mod(ncalls, MPI_TAG_UB_int) !ncall could go over max tag value for long runs
            itag = ncalls
            call MPI_sSend(vbuf, ndata, MPI_DOUBLE_PRECISION, destid, itag, CPL_GRAPH_COMM, ierr)
            !call MPI_iSend(vbuf, ndata, MPI_DOUBLE_PRECISION, destid, itag, CPL_GRAPH_COMM, req(nbr), ierr)

        endif

    enddo

    !call MPI_waitall(nneighbors, req, status, ierr)

    !Barrier for CPL_isend version
    !call MPI_barrier(CPL_GRAPH_COMM, ierr)

end subroutine CPL_send_full

subroutine CPL_send_min(asend, send_flag)
    use coupler_module, only : md_realm, cfd_realm, realm

    real(kind=kind(0.d0)),dimension(:,:,:,:), intent(in):: asend ! Array containing data to send
    logical, intent(out), optional  :: send_flag !Flag set if processor has passed data   

    integer, dimension(6) :: limits ! Global cell indices with minimum and maximum values to send

    if (realm .eq. cfd_realm) then
        call CPL_get_cnst_limits(limits)
    elseif (realm .eq. md_realm) then
        call CPL_get_bnry_limits(limits)
    endif

    call CPL_send_full(asend, limits, send_flag)

end subroutine CPL_send_min


! ----------------------------------------------------------------------------
subroutine CPL_recv_full(arecv, limits, recv_flag)
! ----------------------------------------------------------------------------
!
!Receive data from to local grid from the associated ranks from the other 
!realm
!
!**Remarks**
!
!Assumes the coupler has been initialised with `CPL_init <#f/_/cpl_init>`_ and 
!topological mapping has been setup using either `CPL_setup_md <#f/_/cpl_setup_md>`_ 
!or `CPL_setup_cfd <#f/_/cpl_setup_cfd>`_ as appropriate.
!
!**Synopsis**
!
!.. code-block:: fortran
!
!  CPL_recv(arecv, limits, recv_flag)    
!
!**Inputs**
!
! - arecv
!
!   - Array of data to recv. Should be a four dimensional array allocated using the number of cells on the current processor between the limits. Size should be be obtained from `CPL_my_proc_portion(limits, portion) <#f/coupler/cpl_my_proc_portion>`_.
! - limits 
!
!   - Limits in global cell coordinates, must be the same as corresponding send
!
!**Outputs**
!
! - recv_flag
!
!   - Returned flag which indicates success or failure of recv process
!
!**Example**
!
!.. code-block:: guess
!
!  call CPL_get_olap_limits(limits)
!  call CPL_my_proc_portion(limits, portion)
!  call CPL_get_no_cells(portion, Ncells)
!
!  !Coupled Recieve
!  allocate(A(3, Ncells(1), Ncells(2), Ncells(3)))
!  call CPL_recv(A, limits, recv_flag)
!
! .. sectionauthor::Edward Smith
!
    use mpi
    use coupler_module, only : md_realm, cfd_realm, rank_graph, &
                               error_abort,CPL_GRAPH_COMM,myid_graph,olap_mask, &
                               rank_world, realm, realm_name, & 
                               iblock_realm,jblock_realm,kblock_realm,VOID,ierr, &
                               REALM_NAME, realm, & 
                               rank_intersect, CPL_REALM_INTERSECTION_COMM
    implicit none

    logical, intent(out), optional  :: recv_flag  !Flag set if processor has received data
    integer, intent(in), dimension(6) :: limits  ! Global cell indices with minimum and maximum values to recv
    real(kind(0.d0)), dimension(:,:,:,:),intent(inout):: arecv ! Pre allocated array that recieves data 

    !Neighbours
    integer                             :: nneighbors   
    integer,dimension(:),allocatable    :: id_neighbors
                                                         
    ! local indices 
    integer :: n,nbr,icell,jcell,kcell
    integer :: pos,iclmin,iclmax,jclmin,jclmax,kclmin,kclmax
    integer :: pcoords(3),npercell,ndata,ncells
    integer,dimension(6) :: portion, myportion

    ! auxiliaries 
    logical :: ATTR_FLAG
    integer :: itag, sourceid,start_address, MPI_TAG_UB_int
    integer,dimension(:),allocatable   :: req
    integer,dimension(:,:),allocatable :: status
    real(kind(0.d0)),dimension(:), allocatable ::  vbuf

    !Check counter
    integer, save :: first_counter = 4, check_recv, ncalls=0

    ! This local CFD domain is outside MD overlap zone 
    if (olap_mask(rank_world).eqv. .false.) return

    ! Number of components at each grid point
    npercell = size(arecv,1)

    !Check send/recv are consistent size
    if (first_counter .ne. 0) then
        check_recv = 0
        first_counter = first_counter - 1
        call MPI_AllReduce(size(arecv,1), check_recv, 1, MPI_INTEGER, & 
                        MPI_SUM, CPL_REALM_INTERSECTION_COMM, ierr)
        call MPI_comm_size(CPL_REALM_INTERSECTION_COMM, n, ierr)

        !print*, "Recv npercell=", npercell, check_recv, n
        if (check_recv/n .ne. npercell) then
            print*, "Error in CPL recv in ", realm_name(realm), & 
                    " realm, expected recv size: ", npercell, " but got ", check_recv/n
            call error_abort("Error in CPL recv -- first index of recv array does not match sent data size")
        endif
    endif

    ! Get local grid box ranges seen by this rank for CFD and allocate buffer
    if (realm .eq. cfd_realm) then 

        !Load CFD cells per processor
        call CPL_Cart_coords(CPL_GRAPH_COMM, rank_graph, cfd_realm, 3, pcoords, ierr) 

        ! If limits passed to recv routine, use these instead
        ! of overlap/processor limits
        call CPL_proc_portion(pcoords, cfd_realm, limits, portion, ncells)

        ! Amount of data to receive from all MD processors
        ndata = npercell * ncells
        allocate(vbuf(ndata)); vbuf = 0.d0

    endif    

    !Keep track of number of recvs
    ncalls = ncalls + 1
  
    !Get neighbours
    call MPI_Graph_neighbors_count(CPL_GRAPH_COMM,myid_graph,nneighbors,ierr)
    allocate(id_neighbors(nneighbors))
    call MPI_Graph_neighbors(CPL_GRAPH_COMM,myid_graph,nneighbors,id_neighbors,ierr )

    ! Receive from all attached processors
    allocate(req(nneighbors)); req = MPI_REQUEST_NULL
    allocate(status(MPI_STATUS_SIZE, nneighbors)); 
    start_address = 1 
    do nbr = 1, nneighbors

        !Get source processor from topology graph
        sourceid =  id_neighbors(nbr)

        ! Get size of data to receive from source processors
        if (realm .eq. cfd_realm) then
            !CFD realm receives data based on size of MD processor domain
            call CPL_Cart_coords(CPL_GRAPH_COMM, sourceid+1, md_realm, 3, pcoords, ierr) 
        elseif (realm .eq. md_realm) then
            !MD realm receives data as big as own processor domain
            pcoords = (/iblock_realm, jblock_realm, kblock_realm /)
        endif

        ! If limits passed to recv routine, use these instead
        ! of overlap/processor limits
        call CPL_proc_portion(pcoords, md_realm, limits, portion, ncells)

        !Only receive if overlapping
        if (any(portion.eq.VOID)) then
            ndata = 0
            if (req(nbr) .ne. MPI_REQUEST_NULL) then
                call MPI_Request_free(req(nbr), ierr)
            endif   
            if (present(recv_flag)) recv_flag = .false.
        else
            ! Amount of data to receive
            ndata = npercell * ncells
            if (present(recv_flag)) recv_flag = .true.
            
            if (realm .eq. md_realm) then
                !Allocate receive buffer for MD processors
                if (allocated(vbuf)) deallocate(vbuf); allocate(vbuf(ndata)); vbuf=0.d0
            endif

            ! Receive section of data
            !Using MPI_TAG_UB directly causes error with OpenMPI
            !call MPI_COMM_GET_ATTR(CPL_GRAPH_COMM, MPI_TAG_UB, MPI_TAG_UB_int, ATTR_FLAG, ierr)
            !itag = mod(ncalls, MPI_TAG_UB_int) !ncall could go over max tag value for long runs
            itag = ncalls 
            call MPI_irecv(vbuf(start_address), ndata, MPI_DOUBLE_PRECISION, sourceid, itag, &
                                    CPL_GRAPH_COMM, req(nbr), ierr)

        endif

        !Increment pointer ready to receive next piece of data      
        start_address = start_address + ndata

    enddo
    call MPI_waitall(nneighbors, req, status, ierr)

    !free all requests
!    do nbr = 1, nneighbors
!    print*, req(nbr)
!        if (req(nbr) .ne. MPI_REQUEST_NULL) then
!            call MPI_Request_free(req(nbr), ierr)
!        endif
!    enddo
!    deallocate(req)

    !if (rank_world .eq. 33) then
    !   do n = 1,size(vbuf)
    !       write(98000+rank_world,*) rank_world,n, vbuf(n)
    !   enddo
    !endif

    !Get current processors portion
    call CPL_my_proc_portion(limits, myportion)

    !Check size is consistent with array
    if (myportion(2)-myportion(1)+1 .ne. size(arecv,2)) then
        print*, myportion(2)-myportion(1)+1, size(arecv,2)
        call error_abort("Error in CPL recv -- x size of arecv must be the same as portion(2)-portion(1)+1")
    endif
    if (myportion(4)-myportion(3)+1 .ne. size(arecv,3)) then
        print*, myportion(4)-myportion(3)+1, size(arecv,3)
        call error_abort("Error in CPL recv -- y size of arecv must be the same as portion(4)-portion(3)+1")
    endif
    if (myportion(6)-myportion(5)+1 .ne. size(arecv,4)) then
        print*, myportion(6)-myportion(5)+1, size(arecv,4)
        call error_abort("Error in CPL recv -- z size of arecv must be the same as portion(6)-portion(5)+1")
    endif

    ! ----------------- Unpack data -----------------------------
    start_address = 1
    do nbr = 1, nneighbors

        !Get source processor from topology graph
        sourceid =  id_neighbors(nbr)

        ! Get size of data to receive from source processors
        if (realm .eq. cfd_realm) then
            !CFD always receives data
            if (present(recv_flag)) recv_flag = .true.
            !CFD realm receives data based on size of MD processor domain
            call CPL_Cart_coords(CPL_GRAPH_COMM, sourceid+1, md_realm, 3, pcoords, ierr)
            call CPL_proc_portion(pcoords, md_realm, limits, portion, ncells)
        elseif (realm .eq. md_realm) then
            !MD realm receives data as big as own processor domain
            call CPL_my_proc_portion(limits, portion)
        endif

        iclmin = portion(1)-myportion(1)+1;   iclmax = portion(2)-myportion(1)+1
        jclmin = portion(3)-myportion(3)+1;   jclmax = portion(4)-myportion(3)+1
        kclmin = portion(5)-myportion(5)+1;   kclmax = portion(6)-myportion(5)+1
                
        ! Unpack array into buffer
        if (any(portion.eq.VOID)) then
            !print*, 'VOID recv qqqq',realm_name(realm),rank_world,rank_realm,rank_graph2rank_world(sourceid+1),recv_flag
            ndata = 0
        else
            ! Get local extents in received region
            pos = start_address; ndata = npercell * ncells
!            print'(a,3i4,4i10,12i4,l)', 'recv qqqq', rank_world, rank_realm, rank_olap, &
!                                       ndata, nbr, size(arecv), start_address,          &
!                                       iclmin, iclmax, jclmin, jclmax, kclmin, kclmax,  &
!                                       portion, recv_flag

            do kcell=kclmin,kclmax
            do jcell=jclmin,jclmax
            do icell=iclmin,iclmax
            do n = 1,npercell
                arecv(n,icell,jcell,kcell) = vbuf(pos)
                !print'(6i5,f15.5)', rank_world,pos,n,icell,jcell,kcell,vbuf(pos)
                pos = pos + 1
            end do
            end do
            end do
            end do


        endif

        !Increment reading position of packed data
        start_address = start_address + ndata

    enddo
    ! ----------------- Unpack data -----------------------------

    !Barrier for CPL_isend version
    !call MPI_barrier(CPL_GRAPH_COMM, ierr)
           
end subroutine CPL_recv_full


subroutine CPL_recv_min(arecv, recv_flag)
    use coupler_module, only : md_realm, cfd_realm, realm

    ! Array containing data to recv
    real(kind=kind(0.d0)),dimension(:,:,:,:), intent(out):: arecv 
    logical, intent(out), optional  :: recv_flag !Flag set if processor has passed data   

    integer, dimension(6) :: limits ! Global cell indices with minimum and maximum values to recv

    if (realm .eq. cfd_realm) then
        call CPL_get_bnry_limits(limits)
    elseif (realm .eq. md_realm) then
        call CPL_get_cnst_limits(limits)
    endif

    call CPL_recv_full(arecv, limits, recv_flag)

end subroutine CPL_recv_min

! ----------------------------------------------------------------------------
subroutine CPL_swaphalos(A)
! ----------------------------------------------------------------------------
!
!Swap halos with adjacent processors on current realm cart comm
!
!**Remarks**
!
!Assumes the coupler has been initialised with `CPL_init <#f/_/cpl_init>`_ and 
!topological mapping has been setup using either `CPL_setup_md <#f/_/cpl_setup_md>`_ 
!or `CPL_setup_cfd <#f/_/cpl_setup_cfd>`_ as appropriate.
!
!
!**Synopsis**
!
!.. code-block:: fortran
!
!  CPL_swaphalos(A)    
!
!**Inputs**
!
! - A
!
!   - Array of data to swaphalo. Should be a four dimensional array allocated using the number of cells including one halos on the current processor. Size should be be obtained from `CPL_my_proc_portion(limits, portion) <#f/coupler/cpl_my_proc_portion>`_ with an extra halo on each side. This halo should have been used to collect overflow on current processor.
! 
!
!**Outputs**
!
! - A
!
!   - Array of data with outer halo obtained from adjacent processes.
!
!**Example**
!
!.. code-block:: fortran
!
!  call CPL_get(icmax_olap=icmax_olap)
!  call CPL_get_olap_limits(limits)
!  call CPL_my_proc_portion(limits, portion)
!  call CPL_get_no_cells(portion, Ncells)
!
!  Setup an increasing number in each cell in 2D
!  allocate(A(1, Ncells(1)+2, Ncells(2)+2, Ncells(3)+2))
!  A = 0.0d0
!  do i=2,Ncells(1)+1
!  do j=2,Ncells(2)+1
!  do k=2,Ncells(3)+1
!      ii = i + portion(1)
!      jj = j + portion(3)
!      A(1,i,j,k) = ii + jj*icmax_olap
!  enddo
!  enddo
!  enddo
!  call CPL_swaphalo(A)
!
! .. sectionauthor::Edward Smith
!
    use coupler_module, only : CPL_CART_COMM
    implicit none

    real(kind=kind(0.d0)), intent(inout) :: A(:,:,:,:) !4D Array of data with empty outer halo.

    integer                                                 :: nresults
    integer                                                 :: n,i,j,k,ic,jc,kc
    real(kind=kind(0.d0)),dimension(:,:,:,:),allocatable    :: buf

    integer                                                 :: nhalo,na(3),nxyz(3),nhb(3)
    integer,dimension(:,:),allocatable                      :: halo

    nresults = size(A,1)
    na = (/ size(A,2), size(A,3), size(A,4) /)
    nhb = (/1,1,1/)
    nxyz = na - 2*nhb

    !Get halo cells, input is cell limits minus halo
    call establish_halo_cells(nxyz, halo)
    nhalo = size(halo,1)

    !Pack bins into array of cells
    allocate(buf(na(1),na(2),na(3),nresults))
    buf = 0.d0
    do i=1,size(A,2)
    do j=1,size(A,3)
    do k=1,size(A,4)
        do n=1,nresults
            buf(i, j, k, n) = A(n, i, j, k)
        enddo
    enddo
    enddo
    enddo

    !Copy halo values to adjacent cells in array ready to send
    do n = 1, nhalo
        i = halo(n,1); j = halo(n,2); k = halo(n,3)

        !Change in number of Molecules in halo cells
        ic = i + int(heaviside(na(1)-1-i)-heaviside(i-2))
        jc = j + int(heaviside(na(2)-1-j)-heaviside(j-2))
        kc = k + int(heaviside(na(3)-1-k)-heaviside(k-2))

        buf(i,j,k,:) = buf(i,j,k,:) + buf(ic,jc,kc,:)

    enddo

    !Exchange faces with adjacent processors
    call updatefaces(buf,na(1),na(2),na(3),nresults,1,CPL_CART_COMM)
    call updatefaces(buf,na(1),na(2),na(3),nresults,2,CPL_CART_COMM)
    call updatefaces(buf,na(1),na(2),na(3),nresults,3,CPL_CART_COMM)

    !halo values to correct cells in array
!    do n = 1, nhalo
!        i = halo(n,1); j = halo(n,2); k = halo(n,3)

!        !Change in number of Molecules in halo cells
!        ic = i + int(heaviside(na(1)-1-i)-heaviside(i-2))
!        jc = j + int(heaviside(na(2)-1-j)-heaviside(j-2))
!        kc = k + int(heaviside(na(3)-1-k)-heaviside(k-2))

!        buf(ic,jc,kc,:) = buf(ic,jc,kc,:) + buf(i,j,k,:)

!    enddo

    !Unpack array of cells
    do i=1,size(A,2)
    do j=1,size(A,3)
    do k=1,size(A,4)
        do n=1,nresults
            A(n, i, j, k) = buf(i, j, k, n)
        enddo
    enddo
    enddo
    enddo
    deallocate(buf)

end subroutine CPL_swaphalos


subroutine updatefaces(A, n1, n2, n3, nresults, ixyz, icomm_grid)
!
! updatefaces --  Facilitate the MPI based exchange of data
! Update face halo cells by passing to neighbours
!
    use mpi
    implicit none

    integer,intent(in)                                :: n1,n2,n3,nresults,icomm_grid
    real(kind(0.d0)),intent(inout)                    :: A(:,:,:,:)

    integer                                         :: ixyz, sendrecv_tag, ierr
    integer                                         :: icount,isource,idest
    real(kind(0.d0)),dimension(:,:,:,:),allocatable    :: buf1, buf2

    !Choose a unique sendrecv_tag for this operation
    sendrecv_tag = 1002

    !Determine size of send buffer and copy to buffer data to pass to lower neighbour
    select case (ixyz)
    case (1)
        allocate(buf1(1,n2,n3,nresults), buf2(1,n2,n3,nresults))
        icount = 1*n2*n3*nresults
        buf2 = 0.d0
        buf1(1,:,:,:) = A(1,:,:,:)
    case (2)
        allocate(buf1(n1,1,n3,nresults), buf2(n1,1,n3,nresults))
        icount = n1*1*n3*nresults
        buf2 = 0.d0
        buf1(:,1,:,:) = A(:,1,:,:)
    case (3)
        allocate(buf1(n1,n2,1,nresults), buf2(n1,n2,1,nresults))
        icount = n1*n2*1*nresults
        buf2 = 0.d0
        buf1(:,:,1,:) = A(:,:,1,:)
    case default
        stop "updateBorder: invalid value for ixyz"
    end select

    ! Send to lower neighbor
    call MPI_Cart_shift(icomm_grid, ixyz-1, -1, isource, idest, ierr)
    call MPI_sendrecv(buf1, icount, MPI_double_precision, idest,   sendrecv_tag, &
                      buf2, icount, MPI_double_precision, isource, sendrecv_tag, &
                      icomm_grid, MPI_STATUS_IGNORE, ierr)

    !Save recieved data from upper neighbour and copy to buffer data to pass to upper neighbour
    select case (ixyz)
    case (1)
        buf1(1,:,:,:)= A(n1,:,:,:)
        A(n1,:,:,:)  = buf2(1,:,:,:)
    case (2)
        buf1(:,1,:,:)= A(:,n2,:,:)
        A(:,n2,:,:)  = buf2(:,1,:,:)

    case (3)
        buf1(:,:,1,:)= A(:,:,n3,:)
        A(:,:,n3,:)  = buf2(:,:,1,:)
    end select

    ! Send to upper neighbor
    call MPI_Cart_shift(icomm_grid, ixyz-1, +1, isource, idest, ierr)
    call MPI_sendrecv(buf1, icount, MPI_double_precision, idest,   sendrecv_tag, &
                      buf2, icount, MPI_double_precision, isource, sendrecv_tag, &
                      icomm_grid, MPI_STATUS_IGNORE, ierr)


    !Save recieved data from lower neighbour
    select case (ixyz)
    case (1)
        A(1,:,:,:) = buf2(1,:,:,:)
    case (2)
        A(:,1,:,:) = buf2(:,1,:,:) 
    case (3)
        A(:,:,1,:) = buf2(:,:,1,:)
    end select

    deallocate(buf1, buf2)

end subroutine updatefaces

!-------------------------------------------------------------------
subroutine establish_surface_cells(ncells, nhb, surfacecells)
!
!Establish and store indices of cells which are on the outer domain
!but not in the halo
!
    implicit none

    integer, dimension(3), intent(in) :: ncells !Number of cells in x,y and z
    integer, dimension(3), intent(in) :: nhb    !Number of halos in x,y and z
    integer, dimension(:,:), allocatable, intent(out) :: surfacecells !N by 3 array of Surface cells

    integer        :: n, nsurfacecells
    integer        :: icell, jcell, kcell

    nsurfacecells=    2*( ncells(1)   * ncells(2) &
                    +  (ncells(3)-2)* ncells(2) &
                    +  (ncells(3)-2)*(ncells(1)-2))

    allocate(surfacecells(nsurfacecells,3))

    n = 1
    do kcell=1, ncells(3)+2
    do jcell=1, ncells(2)+2
    do icell=1, ncells(1)+2

        !Remove inner part of domain
        if((icell .gt. 1+nhb(1) .and. icell .lt. ncells(1)+nhb(1)) .and. &
           (jcell .gt. 1+nhb(2) .and. jcell .lt. ncells(2)+nhb(2)) .and. &
           (kcell .gt. 1+nhb(3) .and. kcell .lt. ncells(3)+nhb(2))) cycle
        !Remove outer cells leaving only 1 layer of surface cells
        if((icell .lt. 1+nhb(1) .or. icell .gt. ncells(1)+nhb(1)) .or. &
           (jcell .lt. 1+nhb(2) .or. jcell .gt. ncells(2)+nhb(2)) .or. &
           (kcell .lt. 1+nhb(3) .or. kcell .gt. ncells(3)+nhb(3))) cycle

        surfacecells(n,1)=icell
        surfacecells(n,2)=jcell
        surfacecells(n,3)=kcell
        n = n + 1

    enddo
    enddo
    enddo

end subroutine establish_surface_cells


!-------------------------------------------------------------------
subroutine establish_halo_cells(ncells, halocells)
!
!Establish and store indices of cells which are in the halo
!
    implicit none

    integer, dimension(3), intent(in) :: ncells  !Number of cells in x,y and z
    integer, dimension(:,:), allocatable, intent(out) :: halocells !N by 3 array of halo cells

    integer        :: n, nhalocells
    integer        :: icell, jcell, kcell

    nhalocells  =    2*((ncells(1)+2)*(ncells(2)+2) &
                    +  (ncells(3)  )*(ncells(2)+2) &
                    +  (ncells(3)  )*(ncells(1)  ))

    allocate(halocells(nhalocells,3))

    n = 1
    do kcell=1, ncells(3)+2
    do jcell=1, ncells(2)+2
    do icell=1, ncells(1)+2

        !Remove inner part of domain
        if((icell .gt. (1) .and. icell .lt. (ncells(1)+2)) .and. &
           (jcell .gt. (1) .and. jcell .lt. (ncells(2)+2)) .and. &
           (kcell .gt. (1) .and. kcell .lt. (ncells(3)+2))) cycle

        halocells(n,1)=icell
        halocells(n,2)=jcell
        halocells(n,3)=kcell
        n = n + 1

    enddo
    enddo
    enddo

end subroutine establish_halo_cells


function heaviside(x)
    implicit none

    real(kind(0.d0))                :: heaviside
    integer    ,intent(in)                :: x

    heaviside = nint(0.5*sign(1,x)+1)

end function



!-------------------------------------------------------------------

!subroutine CPL_pack(unpacked,packed,realm,icmax_pack,icmin_pack,jcmax_pack, & 
!                                          jcmin_pack,kcmax_pack,kcmin_pack    )
!    use coupler_module, only : CPL_CART_COMM,rank_cart,md_realm,cfd_realm, & 
!                               error_abort,CPL_GRAPH_COMM,myid_graph,realm_name
!    implicit none

!    integer, intent(in)                                         :: realm
!    real(kind=kind(0.d0)),dimension(:,:,:,:), intent(in)        :: unpacked
!    real(kind=kind(0.d0)),dimension(:),allocatable, intent(out) :: packed

!    ! Optional minimum and maximum values to pack
!    integer, intent(in), optional    :: icmax_pack,icmin_pack
!    integer, intent(in), optional    :: jcmax_pack,jcmin_pack
!    integer, intent(in), optional    :: kcmax_pack,kcmin_pack

!    integer                          :: pos,n,nbr,id_nbr,icell,jcell,kcell,ierr
!    integer                          :: npercell,ncells,nneighbors
!    integer,dimension(3)             :: coord
!    integer,dimension(6)             :: extents,gextents
!    integer,dimension(:),allocatable :: id_neighbors

!    !Amount of data per cell
!    npercell = size(unpacked,1)

!    !Allocate packing buffer
!    if (allocated(packed)) deallocate(packed)
!    allocate(packed(size(unpacked)))

!    ! Get neighbour topology to determine ordering of packed data
!    call MPI_Graph_neighbors_count(CPL_GRAPH_COMM,myid_graph,nneighbors,ierr)
!    allocate(id_neighbors(nneighbors))
!    call MPI_Graph_neighbors(CPL_GRAPH_COMM,myid_graph,nneighbors,id_neighbors,ierr)

!    ! Loop through all neighbours which will be sent to and order the data 
!    ! appropriately to send each correctly
!    do nbr = 1,nneighbors

!        if (realm .eq. cfd_realm) then

!            ! Get MD neighbour
!            id_nbr = id_neighbors(nbr)
!            call CPL_Cart_coords(CPL_GRAPH_COMM,id_nbr+1,md_realm,3,coord,ierr) 
!            call CPL_olap_extents(coord,md_realm,extents,ncells)

!            ! Get offset of neighbouring processor
!            pos = id_nbr * npercell * ncells

!            !print*,'Pack',rank_cart,realm_name(realm),coord,nbr,id_nbr,extents,coord

!        elseif (realm .eq. md_realm) then
!            !Get own processor coordinates 
!            call CPL_Cart_coords(CPL_CART_COMM,rank_cart,realm,3,coord,ierr) 
!            call CPL_olap_extents(coord,realm,gextents,ncells)
!            ! Get local extents
!            extents(1) = 1; extents(2) = gextents(2)-gextents(1)
!            extents(3) = 1; extents(4) = gextents(4)-gextents(3)
!            extents(5) = 1; extents(6) = gextents(6)-gextents(5)
!            pos = 1
!        endif

!        ! Pack array into buffer
!        do kcell=extents(5),extents(6)
!        do jcell=extents(3),extents(4)
!        do icell=extents(1),extents(2)
!        do n = 1,npercell

!            packed(pos) = unpacked(n,icell,jcell,kcell)
!            pos = pos + 1

!        end do
!        end do
!        end do
!        end do

!    end do

!    !Sanity check
!    if (size(packed) .ne. npercell*ncells) then
!        !print*, 'data size', size(packed), 'expected size', npercell*ncells
!        call error_abort("CPL_pack error - cell array does not match expected extents")
!    endif

!end subroutine CPL_pack


!!-------------------------------------------------------------------

!subroutine CPL_unpack(packed,unpacked,realm)
!    use coupler_module, only : CPL_CART_COMM,rank_cart,md_realm,cfd_realm, & 
!                               error_abort,CPL_GRAPH_COMM,myid_graph
!    implicit none

!    integer, intent(in)                                               :: realm
!    real(kind=kind(0.d0)),dimension(:,:,:,:),allocatable, intent(out) :: unpacked
!    real(kind=kind(0.d0)),dimension(:),allocatable, intent(inout)     :: packed

!    integer                          :: pos,n,nbr,id_nbr,icell,jcell,kcell,ierr
!    integer                          :: npercell,ncells,nneighbors
!    integer,dimension(3)             :: coord
!    integer,dimension(6)             :: extents,gextents
!    integer,dimension(:),allocatable :: id_neighbors

!    call CPL_Cart_coords(CPL_CART_COMM,rank_cart,realm,3,coord,ierr) 
!    call CPL_proc_extents(coord,realm,extents,ncells)

!    !Amount of data per cell
!    npercell = size(packed)/ncells

!    !Allocate packing buffer
!    if (allocated(unpacked)) deallocate(unpacked)
!    allocate(unpacked(npercell,1:extents(2)-extents(1), &
!                               1:extents(4)-extents(3), &
!                               1:extents(6)-extents(5)))

!    ! Get neighbour topology to determine ordering of packed data
!    call MPI_Graph_neighbors_count(CPL_GRAPH_COMM,myid_graph,nneighbors,ierr)
!    allocate(id_neighbors(nneighbors))
!    call MPI_Graph_neighbors(CPL_GRAPH_COMM,myid_graph,nneighbors,id_neighbors,ierr)

!    ! Loop through all neighbours which will be sent to and order the data 
!    ! appropriately to send each correctly
!    do nbr = 1,nneighbors

!        if (realm .eq. cfd_realm) then
!            ! Get MD neighbour
!            id_nbr = id_neighbors(nbr)
!            call CPL_Cart_coords(CPL_GRAPH_COMM,id_nbr,md_realm,3,coord,ierr) 
!            call CPL_proc_extents(coord,md_realm,extents,ncells)
!            ! Get offset of neighbouring processor
!            pos = id_nbr * npercell * ncells    !ASSUMES all ncell the same!!
!        elseif (realm .eq. md_realm) then
!            !Get own processor coordinates 
!            call CPL_Cart_coords(CPL_CART_COMM,rank_cart,realm,3,coord,ierr) 
!            call CPL_proc_extents(coord,realm,gextents,ncells)
!            ! Get local extents
!            extents(1) = 1; extents(2) = gextents(2)-gextents(1)
!            extents(3) = 1; extents(4) = gextents(4)-gextents(3)
!            extents(5) = 1; extents(6) = gextents(6)-gextents(5)
!            pos = 1
!        endif

!        ! Unpack buffer into array
!        do kcell=extents(5),extents(6)
!        do jcell=extents(3),extents(4)
!        do icell=extents(1),extents(2)
!        do n = 1,npercell

!            unpacked(n,icell,jcell,kcell) = packed(pos)
!            pos = pos + 1

!        end do
!        end do
!        end do
!        end do

!    end do

!    !Deallocate packed buffer
!    deallocate(packed)

!end subroutine CPL_unpack


!-------------------------------------------------------------------
!                   CPL_proc_extents                  -
!-------------------------------------------------------------------
subroutine CPL_proc_extents(coord, realm, extents, ncells)
!
!Gets maximum and minimum cells for processor coordinates
!
!**Synopsis**
!
!  - CPL_proc_extents(coord,realm,extents,ncells)
!
!**Inputs**
!
!  - coord
!   - processor cartesian coordinate (3 x integer) 
!
!  - realm
!   - cfd_realm (1) or md_realm (2) (integer) 
!
!**Inputs/Output**
!  - NONE
!
!**Output**
!
! - extents
!   - Six components array which defines processor extents
!     xmin,xmax,ymin,ymax,zmin,zmax (6 x integer) 
!
! - ncells (optional)
!   - number of cells on processor (integer) 
!
! .. sectionauthor:: David Trevelyan
! .. sectionauthor:: Edward Smith
!
    use mpi
    use coupler_module, only: md_realm,      cfd_realm,      &
                              icPmin_md,     icPmax_md,      &
                              jcPmin_md,     jcPmax_md,      &
                              kcPmin_md,     kcPmax_md,      &
                              icPmin_cfd,    icPmax_cfd,     &
                              jcPmin_cfd,    jcPmax_cfd,     &
                              kcPmin_cfd,    kcPmax_cfd,     &
                              error_abort!, & 
                              !CPL_setup_complete, REALM_NAME
    implicit none

    character(250)       :: strng
    integer, intent(in)  :: coord(3), realm
    integer, intent(out) :: extents(6)
    integer, optional, intent(out) :: ncells

    !Check setup is complete
    !if (CPL_setup_complete .ne. 1) then
    !    call error_abort("Error CPL_extents/portion called before CPL_setup_"//REALM_NAME(realm))
    !endif

    select case(realm)
    case(md_realm)
        if (size(icPmax_md) .lt. coord(1)) then
            write(strng,'(a,i8,a,i8,a)') 'CPL_proc_extents error - MD proc ',       & 
                                           coord(1), ' extents requested but only ', & 
                                           size(icPmax_md) , ' in x'
            call error_abort(trim(strng))
        endif
        if (size(jcPmax_md) .lt. coord(2)) then
            write(strng,'(a,i8,a,i8,a)') 'CPL_proc_extents error - MD proc ',       & 
                                           coord(2), ' extents requested but only ', & 
                                           size(icPmax_md) , ' in y'
            call error_abort(trim(strng))
        endif
        if (size(kcPmax_md) .lt. coord(3)) then
            write(strng,'(a,i8,a,i8,a)') 'CPL_proc_extents error - MD proc ',       & 
                                           coord(3), ' extents requested but only ', & 
                                           size(icPmax_md) , ' in z'
            call error_abort(trim(strng))
        endif
        extents = (/icPmin_md(coord(1)),icPmax_md(coord(1)), & 
                    jcPmin_md(coord(2)),jcPmax_md(coord(2)), & 
                    kcPmin_md(coord(3)),kcPmax_md(coord(3))/)
    case(cfd_realm)
        if (size(icPmax_cfd) .lt. coord(1)) then
            write(strng,'(a,i8,a,i8,a)') 'CPL_proc_extents error - CFD proc ',      & 
                                           coord(1), ' extents requested but only ', & 
                                           size(icPmax_cfd) , ' in x'
            call error_abort(trim(strng))
        endif
        if (size(jcPmax_cfd) .lt. coord(2)) then
            write(strng,'(a,i8,a,i8,a)') 'CPL_proc_extents error - CFD proc ',      & 
                                           coord(2), ' extents requested but only ', & 
                                           size(icPmax_cfd) , ' in y'
            call error_abort(trim(strng))
        endif
        if (size(kcPmax_cfd) .lt. coord(3)) then
            write(strng,'(a,i8,a,i8,a)') 'CPL_proc_extents error - CFD proc ',      & 
                                           coord(3), ' extents requested but only ', & 
                                           size(icPmax_cfd) , ' in z'
            call error_abort(trim(strng))
        endif
        extents = (/icPmin_cfd(coord(1)),icPmax_cfd(coord(1)), & 
                    jcPmin_cfd(coord(2)),jcPmax_cfd(coord(2)), & 
                    kcPmin_cfd(coord(3)),kcPmax_cfd(coord(3))/)

    case default
        call error_abort('CPL_proc_extents error - Wrong realm in CPL_proc_extents')
    end select

    if (present(ncells)) then
        ncells = (extents(2) - extents(1) + 1) * &
                 (extents(4) - extents(3) + 1) * &
                 (extents(6) - extents(5) + 1)
    end if

end subroutine CPL_proc_extents

!-------------------------------------------------------------------
!                   CPL_olap_extents                  -
!-------------------------------------------------------------------
subroutine CPL_olap_extents(coord,realm,extents,ncells)
!
!Get maximum and minimum cells for current communicator within
!the overlapping region only
!
!**Synopsis** 
!
!  - CPL_olap_extents(coord,realm,extents,ncells)
!
!**Inputs** 
!
!  - coord
!   - processor cartesian coordinate (3 x integer) 
!
!  - realm
!   - cfd_realm (1) or md_realm (2) (integer) 
!
!**Inputs/Output**
!  - NONE
!
!**Output**
!
!  - extents
!   - Six components array which defines processor extents within
!     the overlap region only: xmin,xmax,ymin,ymax,zmin,zmax (6 x integer) 
!
!  - ncells (optional)
!   - number of cells on processor (integer) 
!
! .. sectionauthor:: David Trevelyan
!
    use mpi
    use coupler_module, only: md_realm,      cfd_realm,      &
                              icPmin_md,     icPmax_md,      &
                              jcPmin_md,     jcPmax_md,      &
                              kcPmin_md,     kcPmax_md,      &
                              icPmin_cfd,    icPmax_cfd,     &
                              jcPmin_cfd,    jcPmax_cfd,     &
                              kcPmin_cfd,    kcPmax_cfd,     &
                              icmin_olap,    icmax_olap,     & 
                              jcmin_olap,    jcmax_olap,     & 
                              kcmin_olap,    kcmax_olap,     & 
                              error_abort
    implicit none

    integer, intent(in)  :: coord(3), realm
    integer, intent(out) :: extents(6)
    integer, optional, intent(out) :: ncells

    select case(realm)
    case(md_realm)
        extents(1) = max(icPmin_md(coord(1)),icmin_olap)
        extents(2) = min(icPmax_md(coord(1)),icmax_olap) 
        extents(3) = max(jcPmin_md(coord(2)),jcmin_olap) 
        extents(4) = min(jcPmax_md(coord(2)),jcmax_olap) 
        extents(5) = max(kcPmin_md(coord(3)),kcmin_olap) 
        extents(6) = min(kcPmax_md(coord(3)),kcmax_olap) 
    case(cfd_realm)
        extents(1) = max(icPmin_cfd(coord(1)),icmin_olap)
        extents(2) = min(icPmax_cfd(coord(1)),icmax_olap) 
        extents(3) = max(jcPmin_cfd(coord(2)),jcmin_olap) 
        extents(4) = min(jcPmax_cfd(coord(2)),jcmax_olap) 
        extents(5) = max(kcPmin_cfd(coord(3)),kcmin_olap) 
        extents(6) = min(kcPmax_cfd(coord(3)),kcmax_olap) 
    case default
        call error_abort('CPL_olap_extents error - Wrong realm in CPL_olap_extents')
    end select

    if (present(ncells)) then
        ncells = (extents(2) - extents(1) + 1) * &
                 (extents(4) - extents(3) + 1) * &
                 (extents(6) - extents(5) + 1)
    end if

end subroutine CPL_olap_extents

!-------------------------------------------------------------------
!                   CPL_proc_portion                  -
!-------------------------------------------------------------------
subroutine CPL_proc_portion(coord, realm, limits, portion, ncells)
!
!Get maximum and minimum cell indices, i.e. the 'portion', of the
!input cell extents 'limits' that is contributed by the processor
!specified by processor coord. 
!
!**Remarks**
!
!Assumes the coupler has been initialised with `CPL_init <#f/_/cpl_init>`_ and 
!topological mapping has been setup using either `CPL_setup_md <#f/_/cpl_setup_md>`_ 
!or `CPL_setup_cfd <#f/_/cpl_setup_cfd>`_ as appropriate.
! - Note: limits(6) and portion(6) are of the form: (xmin,xmax,ymin,ymax,zmin,zmax)
!
!**Synopsis**
!
!.. code-block:: fortran
!
!  CPL_proc_portion(coord, realm, limits, portion, ncells)
!
!**Inputs**
!
! - coord
!
!   - processor cartesian coordinate (3 x integer) 
! - realm
!
!   - cfd_realm (1) or md_realm (2) (integer) 
! - limits
!
!   - Array of cell extents that specify the input region. 
!
!
!**Outputs**
!
! - portion
!   - Array of cell extents that define the local processor's
!     contribution to the input region 'limits'.
!
! - ncells (optional)
!    - number of cells in portion (integer) 
!
!
! .. sectionauthor:: David Trevelyan
!
    use mpi
    use coupler_module, only: VOID
    implicit none

    integer, intent(in)  :: coord(3), limits(6),realm
    integer, intent(out) :: portion(6)
    integer, optional, intent(out) :: ncells
    integer :: extents(6)

    call CPL_proc_extents(coord, realm, extents)

    if (extents(1).gt.limits(2) .or. &
        extents(2).lt.limits(1) .or. &
        extents(3).gt.limits(4) .or. &
        extents(4).lt.limits(3) .or. &
        extents(5).gt.limits(6) .or. &
        extents(6).lt.limits(5)) then

        portion(:) = VOID
        if(present(ncells)) ncells = 0
        return
        
    end if

    portion(1) = max(extents(1),limits(1))                          
    portion(2) = min(extents(2),limits(2))                          
    portion(3) = max(extents(3),limits(3))                          
    portion(4) = min(extents(4),limits(4))                          
    portion(5) = max(extents(5),limits(5))                          
    portion(6) = min(extents(6),limits(6))                          

    if (present(ncells)) then
        ncells = (portion(2) - portion(1) + 1) * &
                 (portion(4) - portion(3) + 1) * &
                 (portion(6) - portion(5) + 1)
    end if

end subroutine CPL_proc_portion 


subroutine CPL_my_proc_portion(limits, portion)
!
!Get maximum and minimum cell indices, i.e. the 'portion', of the
!input cell extents 'limits' that is contributed by calling processor.
!
!**Remarks**
!
!Assumes the coupler has been initialised with `CPL_init <#f/_/cpl_init>`_ and 
!topological mapping has been setup using either `CPL_setup_md <#f/_/cpl_setup_md>`_ 
!or `CPL_setup_cfd <#f/_/cpl_setup_cfd>`_ as appropriate.
! - Note: limits(6) and portion(6) are of the form: (xmin,xmax,ymin,ymax,zmin,zmax)
!
!**Synopsis**
!
!.. code-block:: fortran
!
!  CPL_my_proc_portion(limits, portion)
!
!**Inputs**
!
! - limits
!
!   - Array of cell extents that specify the input region. 
!
!**Outputs**
!
! - portion
!   - Array of cell extents that define the local processor's
!     part of the input region 'limits'.
!
! .. sectionauthor:: Eduardo Ramos Fernandez
!
    use coupler_module, only: rank_cart, realm, md_realm, &
                              cfd_realm, rank2coord_cfd, &
                              rank2coord_md
    implicit none

    integer, intent(in)  :: limits(6)
    integer, intent(out) :: portion(6)

    integer :: mycoords(3)

    if (realm .eq. md_realm) then
        mycoords = rank2coord_md(1:3, rank_cart)
    else
        mycoords = rank2coord_cfd(1:3, rank_cart)
    
    end if
    call CPL_proc_portion(mycoords, realm, limits, portion) 

end subroutine CPL_my_proc_portion


subroutine CPL_my_proc_extents(extents)
    use coupler_module, only: rank_cart, realm, md_realm, &
                              cfd_realm, rank2coord_cfd, &
                              rank2coord_md
    implicit none

    integer, intent(out) :: extents(6)

    integer :: mycoords(3)

    if (realm .eq. md_realm) then
        mycoords = rank2coord_md(1:3, rank_cart)
    else
        mycoords = rank2coord_cfd(1:3, rank_cart)
    
    end if
    call CPL_proc_extents(mycoords, realm, extents) 

end subroutine CPL_my_proc_extents


!-------------------------------------------------------------------
!                   CPL_Cart_coords                                -
!-------------------------------------------------------------------

subroutine CPL_Cart_coords(COMM, rank, realm, maxdims, coords, ierr)
!
!Determines process coords in appropriate realm's cartesian topology 
!given a rank in any communicator
!
!**Remarks**
!
!Assumes the coupler has been initialised with `CPL_init <#f/_/cpl_init>`_ and 
!topological mapping has been setup using either `CPL_setup_md <#f/_/cpl_setup_md>`_ 
!or `CPL_setup_cfd <#f/_/cpl_setup_cfd>`_ as appropriate.
!
!
!**Synopsis**
!
!.. code-block:: fortran
!
!  CPL_Cart_coords(COMM, rank, realm, maxdims, coords, ierr)
!
!**Inputs**
!
! - COMM
!
!   - communicator with cartesian structure (handle) 
!
! - rank
!
!   - rank of a process within group of comm (integer) NOTE fortran convention rank=1 to nproc
!
! - realm
!
!   - cfd_realm (1) or md_realm (2) (integer) 
!
! - maxdims
!
!   - length of vector coords in the calling program (integer) 
!
!**Outputs**
!
! - coords
!
!   - integer array (of size ndims) containing the Cartesian coordinates 
!     of specified process (integer) 
! - ierr
!
!   - error flag
! .. sectionauthor:: Edward Smith
!
    use coupler_module, only :  CPL_WORLD_COMM, CPL_REALM_COMM, CPL_INTER_COMM, & 
                                CPL_CART_COMM, CPL_OLAP_COMM, CPL_GRAPH_COMM,   &
                                CPL_REALM_INTERSECTION_COMM, md_realm,cfd_realm, &
                                rank_world2rank_mdcart,     &
                                rank_world2rank_cfdcart,    &
                                rank_mdrealm2rank_world,    &
                                rank_cfdrealm2rank_world,    &
                                rank_olap2rank_world,    &
                                rank_graph2rank_world, &
                                rank2coord_cfd, rank2coord_md, &
                                COUPLER_ERROR_CART_COMM, VOID, nproc_cfd,nproc_md, & 
                                error_abort
    implicit none

    integer, intent(in)     :: COMM, rank, realm, maxdims
    integer, intent(out)    :: coords(maxdims), ierr

    integer                 :: worldrank, cartrank

    !Get rank in world COMM from current COMM
    if (COMM .eq. CPL_WORLD_COMM) then
        ! -  -  -  -  -  -  -  -  -  -  -  -  -
        worldrank = rank
        ! -  -  -  -  -  -  -  -  -  -  -  -  -
    elseif(COMM .eq. CPL_REALM_COMM) then
        ! -  -  -  -  -  -  -  -  -  -  -  -  -
        if (realm .eq. cfd_realm) then
            if (allocated(rank_cfdrealm2rank_world) .eqv. .false.) then
                call error_abort("CPL_Cart_coords Error - " // &
                                 "Setup not complete for CFD CPL_REALM_COMM")
            elseif (rank .gt. size(rank_cfdrealm2rank_world)) then
                print*, 'rank = ', rank, 'comm size = ', size(rank_cfdrealm2rank_world)
                call error_abort("CPL_Cart_coords Error -Specified rank is not in CFD realm")
            endif
            worldrank = rank_cfdrealm2rank_world(rank)
        elseif (realm .eq. md_realm) then
            if (allocated(rank_mdrealm2rank_world) .eqv. .false.) then
                call error_abort("CPL_Cart_coords Error - Setup not " // &
                                 "complete for MD CPL_REALM_COMM")
            elseif (rank .gt. size(rank_mdrealm2rank_world)) then
                print*, 'rank = ', rank, 'comm size = ', size(rank_cfdrealm2rank_world)
                call error_abort("CPL_Cart_coords Error -Specified rank " // &
                                 "is not in MD realm")
            endif
            worldrank = rank_mdrealm2rank_world(rank)
        endif
        ! -  -  -  -  -  -  -  -  -  -  -  -  -
    elseif(COMM .eq. CPL_CART_COMM) then
        ! -  -  -  -  -  -  -  -  -  -  -  -  -
        if (realm .eq. cfd_realm) then
            if (allocated(rank2coord_cfd) .eqv. .false.) &
                call error_abort("CPL_Cart_coords Error - Setup not complete" // &
                                 " for CFD CPL_CART_COMM")
            coords = rank2coord_cfd(:,rank)
        elseif (realm .eq. md_realm) then
            if (allocated(rank2coord_md) .eqv. .false.) &
                call error_abort("CPL_Cart_coords Error - Setup not complete " // &
                                 "for MD CPL_CART_COMM")
            coords = rank2coord_md(:,rank)
        endif
        return
        ! -  -  -  -  -  -  -  -  -  -  -  -  -
    elseif(COMM .eq. CPL_OLAP_COMM) then
        ! -  -  -  -  -  -  -  -  -  -  -  -  -
        if (allocated(rank_olap2rank_world) .eqv. .false.) then
            call error_abort("CPL_Cart_coords Error - Setup not complete for CPL_OLAP_COMM")
        elseif (rank .gt. size(rank_olap2rank_world)) then
            print*, 'rank = ', rank, 'comm size = ', size(rank_olap2rank_world)
            call error_abort("CPL_Cart_coords Error - Specified rank is not in overlap")
        endif
        worldrank = rank_olap2rank_world(rank)
        ! -  -  -  -  -  -  -  -  -  -  -  -  -
    elseif(COMM .eq. CPL_GRAPH_COMM) then
        ! -  -  -  -  -  -  -  -  -  -  -  -  -
        if (allocated(rank_graph2rank_world) .eqv. .false.) then
            call error_abort("CPL_Cart_coords Error - Setup not complete for CPL_GRAPH_COMM")
        elseif (rank .gt. size(rank_graph2rank_world)) then
            call error_abort("CPL_Cart_coords Error - Specified rank is not in graph")
        endif
        worldrank = rank_graph2rank_world(rank)
        ! -  -  -  -  -  -  -  -  -  -  -  -  -
    elseif(COMM .eq. CPL_REALM_INTERSECTION_COMM) then
        ! -  -  -  -  -  -  -  -  -  -  -  -  -
        call error_abort("CPL_Cart_coords Error - Intersection COMM not programmed")
        ! -  -  -  -  -  -  -  -  -  -  -  -  -

    elseif(COMM .eq. CPL_INTER_COMM) then
        ! -  -  -  -  -  -  -  -  -  -  -  -  -
        call error_abort("CPL_Cart_coords Error - Intercomm between realms id" // &
                                 " - use realm comms instead")
        ierr = COUPLER_ERROR_CART_COMM 
        return
        ! -  -  -  -  -  -  -  -  -  -  -  -  -
    else 
        ! -  -  -  -  -  -  -  -  -  -  -  -  -
        call error_abort("CPL_Cart_coords Error - Unknown COMM")
        ierr = COUPLER_ERROR_CART_COMM 
        return
        ! -  -  -  -  -  -  -  -  -  -  -  -  -
    endif
    
    !Get rank in realm cartesian communicator
    if (realm .eq. cfd_realm) then
        if (allocated(rank_world2rank_cfdcart) .eqv. .false.) then
            call error_abort("CPL_Cart_coords Error - world to cart mapping " // &
                                 "not initialised correctly")
        endif
        cartrank = rank_world2rank_cfdcart(worldrank)
        if (cartrank .eq. VOID) call error_abort("CPL_Cart_coords Error - void element in mapping")
        if (cartrank .gt. nproc_cfd) call error_abort("CPL_Cart_coords Error - rank not in cfd realm")
    elseif (realm .eq. md_realm) then
        if (allocated(rank_world2rank_mdcart) .eqv. .false.) then
            call error_abort("CPL_Cart_coords Error - world to cart " // &
                                 "mapping not initialised correctly")
        endif
        cartrank = rank_world2rank_mdcart(worldrank)
        if (cartrank .eq. VOID) call error_abort("CPL_Cart_coords Error - void element in mapping")
        if (cartrank .gt. nproc_md) call error_abort("CPL_Cart_coords Error - rank not in md realm")
    endif

    !Get cartesian coordinate in appropriate realm
    if (realm .eq. cfd_realm) then
        coords = rank2coord_cfd(:,cartrank)
    elseif (realm .eq. md_realm) then
        coords = rank2coord_md(:,cartrank)
    endif

    !Success
    ierr = 0
 
end subroutine CPL_Cart_coords

!-------------------------------------------------------------------
!                   CPL_get_rank                       -
!-------------------------------------------------------------------
subroutine CPL_get_rank(COMM, rank)
!
!Return rank of current processor in specified COMM 
!
!**Remarks**
!
!Assumes the coupler has been initialised with `CPL_init <#f/_/cpl_init>`_ and 
!topological mapping has been setup using either `CPL_setup_md <#f/_/cpl_setup_md>`_ 
!or `CPL_setup_cfd <#f/_/cpl_setup_cfd>`_ as appropriate.
!
!
!**Synopsis**
!
!.. code-block:: fortran
!
!  CPL_get_rank(COMM, rank)
!
!**Inputs**
!
! - COMM
!   - communicator with cartesian structure (handle) 
!
!**Outputs**
!
! - rank
!   - rank of a process within group of comm (integer) 
!      NOTE fortran convention rank=1 to nproc
!
! .. sectionauthor:: Edward Smith
!
    use coupler_module, only :  CPL_WORLD_COMM, CPL_REALM_COMM, CPL_INTER_COMM, & 
                                CPL_CART_COMM,  CPL_OLAP_COMM,  CPL_GRAPH_COMM, &
                                CPL_REALM_INTERSECTION_COMM,rank_world,         &
                                rank_realm,rank_cart,rank_olap,rank_graph,error_abort

    integer, intent(in)     :: COMM
    integer, intent(out)    :: rank

    !Get rank in world COMM from current COMM
    if (COMM .eq. CPL_WORLD_COMM) then
        rank = rank_world
        return
    elseif(COMM .eq. CPL_REALM_COMM) then
        rank = rank_realm
        return
    elseif(COMM .eq. CPL_CART_COMM) then
        rank = rank_cart
        return
    elseif(COMM .eq. CPL_OLAP_COMM) then
        rank = rank_olap
        return
    elseif(COMM .eq. CPL_GRAPH_COMM) then
        rank = rank_graph
        return
    elseif(COMM .eq. CPL_REALM_INTERSECTION_COMM) then
        call error_abort("CPL_get_rank Error - Intersection COMM not programmed")
    elseif(COMM .eq. CPL_INTER_COMM) then
        call error_abort("CPL_get_rank Error - No rank in Intercomm" // &
                                 " - use realm comms instead")
    else 
        call error_abort("CPL_get_rank Error - Unknown COMM")
    endif
    

end subroutine CPL_get_rank


!-------------------------------------------------------------------
!                   CPL_olap_check                                         -
!-------------------------------------------------------------------
function CPL_overlap() result(p)
!
!Check if current processor is in the overlap region
!
!**Synopsis**
!
!  - CPL_olap_check()
!
!**Inputs** Parameters
!
!  - NONE
!
! - Returns
!
!  - CPL_olap_check
!   - True if calling processor is in the overlap region
!     and false otherwise
!
! .. sectionauthor:: Edward Smith
!
    use coupler_module, only : olap_mask, rank_world
    implicit none

    logical :: p

    p = olap_mask(rank_world)

end function CPL_overlap


!============================================================================
!
! Utility functions and subroutines that extract parameters from internal modules 
!
!-----------------------------------------------------------------------------    

!-------------------------------------------------------------------
!                   CPL_get                                        -
!-------------------------------------------------------------------
subroutine CPL_get(icmax_olap,icmin_olap,jcmax_olap,jcmin_olap,  & 
                   kcmax_olap,kcmin_olap,density_cfd,density_md, &
                   dt_cfd,dt_MD,dx,dy,dz,ncx,ncy,ncz,xg,yg,zg,   &
                   xL_md,xL_cfd,yL_md,yL_cfd,zL_md,zL_cfd,       &
                   npx_md,npy_md,npz_md,npx_cfd,npy_cfd,npz_cfd, &
                   constraint_algo, constraint_CVflag,           &
                   constraint_OT, constraint_NCER,               &
                   constraint_Flekkoy, constraint_off,           &
                   constraint_CV,                                &
                   icmin_cnst, icmax_cnst,                       &
                   jcmin_cnst, jcmax_cnst,                       &
                   kcmin_cnst, kcmax_cnst,                       &
                   md_cfd_match_cellsize, staggered_averages,    &
                   cpl_cfd_bc_slice, cpl_md_bc_slice,            &
                   nsteps_md, nsteps_cfd, nsteps_coupled,        &
                   cpl_cfd_bc_x, cpl_cfd_bc_y, cpl_cfd_bc_z,     &
                   timestep_ratio, comm_style,                   &
                   sendtype_cfd_to_md, sendtype_md_to_cfd)
!Wrapper to retrieve (read only) parameters from the coupler_module 
!Note - this ensures all variable in the coupler are protected
!from corruption by either CFD or MD codes
!
!**Synopsis**
!
!  - CPL_get([see coupler_module])
!
!**Input**
!
!  - NONE
!
!**Output**
!
!  - See below
!
! .. sectionauthor:: Edward Smith

    use coupler_module, only :  icmax_olap_=>icmax_olap,         &
                                icmin_olap_=>icmin_olap,         &
                                jcmax_olap_=>jcmax_olap,         &
                                jcmin_olap_=>jcmin_olap,         &
                                kcmax_olap_=>kcmax_olap,         &
                                kcmin_olap_=>kcmin_olap,         &
                                density_cfd_=>density_cfd,       &             
                                density_md_=>density_md,         &
                                dt_cfd_=>dt_cfd,dt_MD_=>dt_MD,   &
                                dx_=>dx,dy_=>dy,dz_=>dz,         &
                                ncx_=>ncx,ncy_=> ncy,ncz_=> ncz, &
                                xL_md_ =>xL_md,                  &
                                yL_md_ =>yL_md,                  &
                                zL_md_ =>zL_md,                  &
                                xL_cfd_=> xL_cfd,                &
                                yL_cfd_=> yL_cfd,                &
                                zL_cfd_=> zL_cfd,                &
                                npx_md_=> npx_md,                &
                                npy_md_=> npy_md,                &
                                npz_md_=> npz_md,                &
                                npx_cfd_=> npx_cfd,                &
                                npy_cfd_=> npy_cfd,                &
                                npz_cfd_=> npz_cfd,                &
                                xg_=>xg, yg_=> yg, zg_=> zg,     &
                                timestep_ratio_ => timestep_ratio, &
                                md_cfd_match_cellsize_ => md_cfd_match_cellsize, &
                                staggered_averages_ => staggered_averages, &
                                constraint_algo_ => constraint_algo,       &
                                constraint_CVflag_ => constraint_CVflag,   &
                                constraint_OT_ => constraint_OT,           &
                                constraint_NCER_ => constraint_NCER,       & 
                                constraint_Flekkoy_ => constraint_Flekkoy, &
                                constraint_CV_ => constraint_CV,           &
                                constraint_off_ => constraint_off,         &
                                icmin_cnst_ => icmin_cnst, &
                                icmax_cnst_ => icmax_cnst, &
                                jcmin_cnst_ => jcmin_cnst, &
                                jcmax_cnst_ => jcmax_cnst, &
                                kcmin_cnst_ => kcmin_cnst, &
                                kcmax_cnst_ => kcmax_cnst, &
                                cpl_cfd_bc_slice_ => cpl_cfd_bc_slice, &
                                cpl_md_bc_slice_ => cpl_md_bc_slice, &
                                nsteps_md_ => nsteps_md, &
                                nsteps_cfd_ => nsteps_cfd, &
                                nsteps_coupled_ => nsteps_coupled, &
                                cpl_cfd_bc_x_ => cpl_cfd_bc_x, &
                                cpl_cfd_bc_y_ => cpl_cfd_bc_y, &
                                cpl_cfd_bc_z_ => cpl_cfd_bc_z, &
                                comm_style_ => comm_style, &
                                comm_style_send_recv_ => comm_style_send_recv,&
                                comm_style_gath_scat_ => comm_style_gath_scat,&
                                sendtype_cfd_to_md_ => sendtype_cfd_to_md,    &
                                sendtype_md_to_cfd_ =>  sendtype_md_to_cfd
    implicit none

    logical,dimension(3),optional,intent(out) :: staggered_averages

    integer, optional, intent(out)          :: icmax_olap
    integer, optional, intent(out)          :: icmin_olap
    integer, optional, intent(out)          :: jcmax_olap
    integer, optional, intent(out)          :: jcmin_olap
    integer, optional, intent(out)          :: kcmax_olap
    integer, optional, intent(out)          :: kcmin_olap
    integer, optional, intent(out)          :: ncx,ncy,ncz
    integer, optional, intent(out)          :: npx_md
    integer, optional, intent(out)          :: npy_md
    integer, optional, intent(out)          :: npz_md
    integer, optional, intent(out)          :: npx_cfd
    integer, optional, intent(out)          :: npy_cfd
    integer, optional, intent(out)          :: npz_cfd
    integer, optional, intent(out)          :: md_cfd_match_cellsize
    integer, optional, intent(out)          :: timestep_ratio

    integer, optional, intent(out)          :: constraint_algo
    integer, optional, intent(out)          :: constraint_CVflag
    integer, optional, intent(out)          :: constraint_OT
    integer, optional, intent(out)          :: constraint_NCER
    integer, optional, intent(out)          :: constraint_Flekkoy
    integer, optional, intent(out)          :: constraint_CV
    integer, optional, intent(out)          :: constraint_off
    integer, optional, intent(out)          :: comm_style 
    integer, optional, intent(out)          :: icmax_cnst
    integer, optional, intent(out)          :: icmin_cnst
    integer, optional, intent(out)          :: jcmax_cnst
    integer, optional, intent(out)          :: jcmin_cnst
    integer, optional, intent(out)          :: kcmax_cnst
    integer, optional, intent(out)          :: kcmin_cnst
    integer, optional, intent(out)          :: cpl_cfd_bc_slice 
    integer, optional, intent(out)          :: cpl_cfd_bc_x 
    integer, optional, intent(out)          :: cpl_cfd_bc_y 
    integer, optional, intent(out)          :: cpl_cfd_bc_z 
    integer, optional, intent(out)          :: cpl_md_bc_slice 
    integer, optional, intent(out)          :: nsteps_md, nsteps_cfd
    integer, optional, intent(out)          :: nsteps_coupled
    integer, optional, intent(out)          :: sendtype_cfd_to_md
    integer, optional, intent(out)          :: sendtype_md_to_cfd

    real(kind(0.d0)), optional, intent(out) :: density_cfd
    real(kind(0.d0)), optional, intent(out) :: density_md
    real(kind(0.d0)), optional, intent(out) :: dt_cfd
    real(kind(0.d0)), optional, intent(out) :: dt_MD
    real(kind(0.d0)), optional, intent(out) :: dx
    real(kind(0.d0)), optional, intent(out) :: dy
    real(kind(0.d0)), optional, intent(out) :: dz
    real(kind(0.d0)), optional, intent(out) :: xL_md
    real(kind(0.d0)), optional, intent(out) :: xL_cfd
    real(kind(0.d0)), optional, intent(out) :: yL_md
    real(kind(0.d0)), optional, intent(out) :: yL_cfd
    real(kind(0.d0)), optional, intent(out) :: zL_md
    real(kind(0.d0)), optional, intent(out) :: zL_cfd

    real(kind(0.d0)), dimension(:,:,:),allocatable,optional,intent(out) :: xg
    real(kind(0.d0)), dimension(:,:,:),allocatable,optional,intent(out) :: yg
    real(kind(0.d0)), dimension(:,:,:),allocatable,optional,intent(out) :: zg

    !Overlap extents
    if (present(icmax_olap)) icmax_olap = icmax_olap_
    if (present(icmin_olap)) icmin_olap = icmin_olap_
    if (present(jcmax_olap)) jcmax_olap = jcmax_olap_
    if (present(jcmin_olap)) jcmin_olap = jcmin_olap_
    if (present(kcmax_olap)) kcmax_olap = kcmax_olap_
    if (present(kcmin_olap)) kcmin_olap = kcmin_olap_

    !Density
    if (present(density_cfd)) density_cfd = density_cfd_
    if (present(density_md )) density_md  = density_md_

    !Cells
    if (present(ncx)) ncx = ncx_
    if (present(ncy)) ncy = ncy_
    if (present(ncz)) ncz = ncz_
            
    !Temporal and spatial steps
    if (present(dt_cfd)) dt_cfd= dt_cfd_
    if (present(dt_MD )) dt_MD = dt_MD_
    if (present(nsteps_md)) nsteps_md = nsteps_md_
    if (present(nsteps_cfd)) nsteps_cfd = nsteps_cfd_
    if (present(nsteps_coupled)) nsteps_coupled= nsteps_coupled_
    if (present(dx)) dx = dx_
    if (present(dy)) dy = dy_   
    if (present(dz)) dz = dz_

    !Domain sizes
    if (present(xL_md )) xL_md = xL_md_
    if (present(xL_cfd)) xL_cfd= xL_cfd_
    if (present(yL_md )) yL_md = yL_md_
    if (present(yL_cfd)) yL_cfd= yL_cfd_
    if (present(zL_md )) zL_md = zL_md_
    if (present(zL_cfd)) zL_cfd= zL_cfd_

    !Number of processors
    if (present(npx_md)) npx_md = npx_md_
    if (present(npy_md)) npy_md = npy_md_
    if (present(npz_md)) npz_md = npz_md_
    if (present(npx_cfd)) npx_cfd = npx_cfd_
    if (present(npy_cfd)) npy_cfd = npy_cfd_
    if (present(npz_cfd)) npz_cfd = npz_cfd_

    !The mesh
    if (present(xg)) then
        allocate(xg(size(xg_,1),size(xg_,2),size(xg_,3)))
        xg = xg_
    endif
    if (present(yg)) then
        allocate(yg(size(yg_,1),size(yg_,2),size(yg_,3)))
        yg = yg_
    endif
    if (present(zg)) then
        allocate(zg(size(zg_,1),size(zg_,2),size(zg_,3)))
        zg = zg_
    endif

    !Coupler input parameters
    if (present(staggered_averages)) staggered_averages = staggered_averages_
    if (present(timestep_ratio)) timestep_ratio = timestep_ratio_
    if (present(md_cfd_match_cellsize)) then
        md_cfd_match_cellsize = md_cfd_match_cellsize_
    end if

    ! Constraint information
    if (present(constraint_algo)) constraint_algo = constraint_algo_
    if (present(constraint_CVflag)) constraint_CVflag = constraint_CVflag_
    if (present(constraint_OT)) constraint_OT = constraint_OT_
    if (present(constraint_NCER)) constraint_NCER = constraint_NCER_
    if (present(constraint_Flekkoy)) constraint_Flekkoy = constraint_Flekkoy_
    if (present(constraint_CV)) constraint_CV = constraint_CV_
    if (present(constraint_off)) constraint_off = constraint_off_
    if (present(icmin_cnst)) icmin_cnst = icmin_cnst_
    if (present(icmax_cnst)) icmax_cnst = icmax_cnst_
    if (present(jcmin_cnst)) jcmin_cnst = jcmin_cnst_
    if (present(jcmax_cnst)) jcmax_cnst = jcmax_cnst_
    if (present(kcmin_cnst)) kcmin_cnst = kcmin_cnst_
    if (present(kcmax_cnst)) kcmax_cnst = kcmax_cnst_

    ! Communication style
    if (present(comm_style)) comm_style = comm_style_
    if (present(sendtype_cfd_to_md)) sendtype_cfd_to_md = sendtype_cfd_to_md_
    if (present(sendtype_md_to_cfd)) sendtype_md_to_cfd = sendtype_md_to_cfd_

    ! Coupling styles
    if (present(cpl_cfd_bc_slice)) cpl_cfd_bc_slice = cpl_cfd_bc_slice_
    if (present(cpl_md_bc_slice)) cpl_md_bc_slice = cpl_md_bc_slice_
    if (present(cpl_cfd_bc_x)) cpl_cfd_bc_x = cpl_cfd_bc_x_
    if (present(cpl_cfd_bc_y)) cpl_cfd_bc_y = cpl_cfd_bc_y_
    if (present(cpl_cfd_bc_z)) cpl_cfd_bc_z = cpl_cfd_bc_z_

end subroutine CPL_get


function CPL_comm_style()
    use coupler_module, only: comm_style
    implicit none

    integer :: CPL_comm_style
    
    CPL_comm_style = comm_style

end function

! Function to get current realm
function CPL_realm()
    use coupler_module, only : realm
    implicit none
    
    integer :: CPL_realm

    CPL_realm = realm

end function CPL_realm


subroutine CPL_meshgrid(xg, yg, zg)
    implicit none

    real(kind(0.d0)), dimension(:,:,:), & 
        allocatable,intent(out) :: xg, yg, zg

    call CPL_get(xg=xg, yg=yg, zg=zg)

end subroutine CPL_meshgrid


!=============================================================================
function CPL_map_md2cfd_coord(coord_md, coord_cfd) result(valid_coord)
!
! Map global MD position to global CFD coordinate frame
!
    use coupler_module, only :  xL_md,xg,icmin_olap,icmax_olap, & 
                                yL_md,yg,jcmin_olap,jcmax_olap, & 
                                zL_md,zg,kcmin_olap,kcmax_olap, &
                                x_orig_md, y_orig_md, z_orig_md,&
                                x_orig_cfd, y_orig_cfd, z_orig_cfd,&
                                xL_cfd, yL_cfd, zL_cfd

    implicit none

    real(kind(0.d0)), intent(out)  :: coord_cfd(3)
    real(kind(0.d0)), intent(in)   :: coord_md(3)
    real(kind(0.d0))               :: md_only(3), coord_limits_cfd(6), &
                                      coord_limits_md(6)
    logical                        :: valid_coord


    coord_limits_md(1) = x_orig_md
    coord_limits_md(2) = x_orig_md + xL_md
    coord_limits_md(3) = y_orig_md
    coord_limits_md(4) = y_orig_md + yL_md
    coord_limits_md(5) = z_orig_md
    coord_limits_md(6) = z_orig_md + zL_md

    valid_coord = is_coord_inside(coord_md, coord_limits_md)
    
    if (valid_coord) then
        !Get size of MD domain which has no CFD cells overlapping
        !This should be general enough to include grid stretching
        !and total overlap in any directions 
        md_only(1) = xL_md-abs(xg(icmax_olap+1,1,1) - xg(icmin_olap,1,1))
        md_only(2) = yL_md-abs(yg(1,jcmax_olap+1,1) - yg(1,jcmin_olap,1))
        md_only(3) = zL_md-abs(zg(1,1,kcmax_olap+1) - zg(1,1,kcmin_olap))

        ! CFD has origin at bottom left while MD origin at centre
        coord_limits_cfd(1) = x_orig_cfd
        coord_limits_cfd(2) = x_orig_cfd + xL_cfd
        coord_limits_cfd(3) = y_orig_cfd
        coord_limits_cfd(4) = y_orig_cfd + yL_cfd
        coord_limits_cfd(5) = z_orig_cfd
        coord_limits_cfd(6) = z_orig_cfd + zL_cfd
        !coord_md = coord_cfd + md_only + md_xyz_orig
        coord_cfd(1) = abs(coord_md(1) - x_orig_md - md_only(1)) + x_orig_cfd
        coord_cfd(2) = abs(coord_md(2) - y_orig_md - md_only(2)) + y_orig_cfd
        coord_cfd(3) = abs(coord_md(3) - z_orig_md - md_only(3)) + z_orig_cfd
        valid_coord = is_coord_inside(coord_cfd, coord_limits_cfd)
    endif

end function CPL_map_md2cfd_coord



!=============================================================================
!-----------------------------------------------------------------------------
function CPL_map_cfd2md_coord(coord_cfd, coord_md) result(valid_coord)
!
!Map global CFD position in global MD coordinate frame
!
    use coupler_module, only :  xL_md,xg,icmin_olap,icmax_olap, & 
                                yL_md,yg,jcmin_olap,jcmax_olap, & 
                                zL_md,zg,kcmin_olap,kcmax_olap, &
                                x_orig_md, y_orig_md, z_orig_md,&
                                x_orig_cfd, y_orig_cfd, z_orig_cfd,&
                                xL_cfd, yL_cfd, zL_cfd

    implicit none

    real(kind(0.d0)), intent(in)  :: coord_cfd(3)
    real(kind(0.d0)), intent(out) :: coord_md(3)
    real(kind(0.d0))              :: md_only(3), coord_limits_cfd(6), &
                                     coord_limits_md(6)
    logical                       :: valid_coord

    coord_limits_cfd(1) = x_orig_cfd
    coord_limits_cfd(2) = x_orig_cfd + xL_cfd
    coord_limits_cfd(3) = y_orig_cfd
    coord_limits_cfd(4) = y_orig_cfd + yL_cfd
    coord_limits_cfd(5) = z_orig_cfd
    coord_limits_cfd(6) = z_orig_cfd + zL_cfd

    coord_limits_md(1) = x_orig_md
    coord_limits_md(2) = x_orig_md + xL_md
    coord_limits_md(3) = y_orig_md
    coord_limits_md(4) = y_orig_md + yL_md
    coord_limits_md(5) = z_orig_md
    coord_limits_md(6) = z_orig_md + zL_md

    valid_coord = is_coord_inside(coord_cfd, coord_limits_cfd)
    
    if (valid_coord) then
        !Get size of MD domain which has no CFD cells overlapping
        !This should be general enough to include grid stretching
        !and total overlap in any directions 
        md_only(1) = xL_md-abs(xg(icmax_olap+1,1,1) - xg(icmin_olap,1,1))
        md_only(2) = yL_md-abs(yg(1,jcmax_olap+1,1) - yg(1,jcmin_olap,1))
        md_only(3) = zL_md-abs(zg(1,1,kcmax_olap+1) - zg(1,1,kcmin_olap))

        ! CFD has origin at bottom left while MD origin at centre
        !coord_md = coord_cfd + md_only + md_xyz_orig
        coord_md(1) = abs(coord_cfd(1) - x_orig_cfd) + md_only(1) + x_orig_md
        coord_md(2) = abs(coord_cfd(2) - y_orig_cfd) + md_only(2) + y_orig_md
        coord_md(3) = abs(coord_cfd(3) - z_orig_cfd) + md_only(3) + z_orig_md
        valid_coord = is_coord_inside(coord_md, coord_limits_md)
    endif

end function CPL_map_cfd2md_coord

function CPL_md2cfd(coord_md) result(coord_cfd)
    use coupler_module, only : error_abort
    implicit none

    real(kind(0.d0)), intent(in)  :: coord_md(3)
    real(kind(0.d0))              :: coord_cfd(3)

    logical :: valid

    valid = CPL_map_md2cfd_coord(coord_md, coord_cfd)
!    if (.not. valid) then
!        print'(a,l,2(a,3f10.5))', "Valid coords= ",valid, " CFD coords= ", coord_cfd, " MD coords= ", coord_md
!!        call error_abort("CPL_md2cfd error - coord outside limits") 
!    endif    

end function CPL_md2cfd


function CPL_cfd2md(coord_cfd) result(coord_md)
    use coupler_module, only : error_abort
    implicit none

    real(kind(0.d0)), intent(in)  :: coord_cfd(3)
    real(kind(0.d0))              :: coord_md(3)

    logical :: valid

    valid = CPL_map_cfd2md_coord(coord_cfd, coord_md)
!    if (.not. valid) then
!        print'(a,l,2(a,3f10.5))', "Valid coords= ",valid, " CFD coords= ", coord_cfd, " MD coords= ", coord_md
!!        call error_abort("CPL_cfd2md error - coord outside limits") 
!    endif    

end function CPL_cfd2md
!-----------------------------------------------------------------------------

subroutine CPL_map_cell2coord(i, j, k, coord_xyz)
    use coupler_module, only: xg, yg, zg, realm, &
                              ncx, ncy, ncz, &
                              md_realm, cfd_realm, & 
                              maxgridsize, error_abort

    integer, intent(in)  :: i, j, k
    real(kind(0.d0)), intent(out) :: coord_xyz(3)

    real(kind(0.d0)) :: coord_md(3)
    integer :: olap_limits(6)
    logical :: aux_ret

    call CPL_get_olap_limits(olap_limits)

    if (.not. is_cell_inside((/i, j, k/), olap_limits)) then
        print*, "cell:", (/i,j,k/)
        call error_abort("CPL_map_cell2coord error - Cell is outside overlap region. " // &
                         "Aborting simulation.") 
    end if

    if (ncx*ncy*ncz .lt. maxgridsize) then
        coord_xyz(1) = xg(i, j, k)
        coord_xyz(2) = yg(i, j, k) 
        coord_xyz(3) = zg(i, j, k)
    else
        coord_xyz(1) = xg(i, 1, 1)
        coord_xyz(2) = yg(1, j, 1) 
        coord_xyz(3) = zg(1, 1, k)
    endif

    if (realm .eq. md_realm) then
        aux_ret = CPL_map_cfd2md_coord(coord_xyz, coord_md)
        coord_xyz = coord_md
    end if


end subroutine CPL_map_cell2coord

!-----------------------------------------------------------------------------

 !function CPL_map_coord2cell(x, y, z, cell_ijk) result(ret)
 !
 !   use coupler_module, only: dx, dy, dz, &
 !                              icmin_olap, jcmin_olap, kcmin_olap
 !
 !   real(kind(0.d0)), intent(in)  :: x, y, z
 !   integer, intent(out)         :: cell_ijk(3)
 ! 
 !   integer          :: ixyz
 !   real(kind(0.d0)) :: olap_lo(3)
 !   integer          :: olap_limits(6)
 !    logical          :: ret
 !
 !   call CPL_map_cell2coord(icmin_olap, jcmin_olap, kcmin_olap, olap_lo)
 !
 !   cell_ijk(1) = ceiling((x - olap_lo(1)) / dx)
 !    cell_ijk(2) = ceiling((y - olap_lo(2)) / dy)
 !   cell_ijk(3) = ceiling((z - olap_lo(3)) / dz)
 !
 !   call CPL_get_olap_limits(olap_limits)
 !
 !   ret = .true.
 !    if (.not. is_cell_inside(cell_ijk, olap_limits)) then
 !       ret = .false.
 !   end if
 !   
 !end function CPL_map_coord2cell
 !
 !
 
!!-----------------------------------------------------------------------------
!! This is the new version introduced by Edu which causes an error in my code
function CPL_map_coord2cell(x, y, z, cell_ijk) result(ret)
 
    use coupler_module, only: dx, dy, dz, &
                            icmin_olap, jcmin_olap, kcmin_olap, &
                              icmax_olap, jcmax_olap, kcmax_olap

    real(kind(0.d0)), intent(in)  :: x, y, z
    integer, intent(out)         :: cell_ijk(3)

    real(kind(0.d0)) :: olap_lo(3), olap_hi(3)
    integer          :: olap_limits(6)
    logical          :: ret

    ! Returns coordinates in the correct realm already
    call CPL_map_cell2coord(icmin_olap, jcmin_olap, kcmin_olap, olap_lo)
    cell_ijk(1) = floor(abs(x - olap_lo(1)) / dx) + 1
    cell_ijk(2) = floor(abs(y - olap_lo(2)) / dy) + 1
    cell_ijk(3) = floor(abs(z - olap_lo(3)) / dz) + 1

    ! Returns coordinates in the correct realm already
    call CPL_map_cell2coord(icmax_olap, jcmax_olap, kcmax_olap, olap_hi)
    ! Get the maximun coordinate in the domain
    olap_hi = olap_hi + (/dx, dy, dz/)

    !NOTE: The highest coordinate would give a cell number outside the 
    !      grid so one has to be subtracted from the cell number.

    if (x == olap_hi(1)) then
        cell_ijk(1) = cell_ijk(1) - 1
    end if
    if (y == olap_hi(2)) then
        cell_ijk(2) = cell_ijk(2) - 1
    end if
    if (z == olap_hi(3)) then
        cell_ijk(3) = cell_ijk(3) - 1
    end if

    call CPL_get_olap_limits(olap_limits)

    ret = .true.
    if (.not. is_cell_inside(cell_ijk, olap_limits)) then
        ret = .false.
    end if
    
end function CPL_map_coord2cell



!-----------------------------------------------------------------------------

subroutine CPL_get_no_cells(limits, no_cells)

   integer, intent(in)  :: limits(6)
   integer, intent(out) :: no_cells(3)

   ! TODO: Check for limits
   
   no_cells(1) = limits(2) - limits(1) + 1
   no_cells(2) = limits(4) - limits(3) + 1
   no_cells(3) = limits(6) - limits(5) + 1

end subroutine CPL_get_no_cells

!-----------------------------------------------------------------------------
!-----------------------------------------------------------------------------
function CPL_map_glob2loc_cell(limits, glob_cell, loc_cell) result(ret)

    use coupler_module, only :  VOID, error_abort
                                

    integer, intent(in)  :: limits(6)
    integer, intent(in)  :: glob_cell(3)
    integer, intent(out) :: loc_cell(3)

    logical :: ret

    ! Check if global cell is within the limits of the region specified
    if (is_cell_inside(glob_cell, limits)) then
        loc_cell(1) = glob_cell(1) - limits(1) + 1
        loc_cell(2) = glob_cell(2) - limits(3) + 1
        loc_cell(3) = glob_cell(3) - limits(5) + 1
        ret = .true.
    else
        loc_cell = VOID
        ret = .false.
    end if

end function CPL_map_glob2loc_cell

!-----------------------------------------------------------------------------

subroutine CPL_get_olap_limits(limits)
! ----------------------------------------------------------------------------
!Get limits of overlap region as cell indices in global coordinate system, 
!these are as specified in the input file (COUPLER.in).
!
!**Remarks**
!
!Assumes the coupler has been initialised with `CPL_init <#f/_/cpl_init>`_ and 
!topological mapping has been setup using either `CPL_setup_md <#f/_/cpl_setup_md>`_ 
!or `CPL_setup_cfd <#f/_/cpl_setup_cfd>`_ as appropriate.
! - Note: limits(6) are of the form: (xmin,xmax,ymin,ymax,zmin,zmax)
!
!**Synopsis**
!
!.. code-block:: fortran
!
!  CPL_get_olap_limits(limits)
!
!**Outputs**
!
! - limits
!
!   - Array of cell extents that specify the overlap region. 
!
! .. sectionauthor:: Edward Smith
!
    use coupler_module, only :  icmin_olap, icmax_olap, & 
                                jcmin_olap, jcmax_olap, & 
                                kcmin_olap, kcmax_olap

   integer, intent(out) :: limits(6)

   limits(1) = icmin_olap
   limits(2) = icmax_olap
   limits(3) = jcmin_olap
   limits(4) = jcmax_olap
   limits(5) = kcmin_olap
   limits(6) = kcmax_olap

end subroutine CPL_get_olap_limits

! ----------------------------------------------------------------------------
subroutine CPL_get_cnst_limits(limits)
!
!Get limits of constraint region as cell indices in global coordinate system, 
!these are as specified in the input file (COUPLER.in) and must be in the
!overlap region.
!
!**Remarks**
!
!Assumes the coupler has been initialised with `CPL_init <#f/_/cpl_init>`_ and 
!topological mapping has been setup using either `CPL_setup_md <#f/_/cpl_setup_md>`_ 
!or `CPL_setup_cfd <#f/_/cpl_setup_cfd>`_ as appropriate.
! - Note: limits(6) are of the form: (xmin,xmax,ymin,ymax,zmin,zmax)
!
!**Synopsis**
!
!.. code-block:: fortran
!
!  CPL_get_cnst_limits(limits)
!
!**Outputs**
!
! - limits
!
!   - Array of cell extents that specify the constrained region. 
!
! .. sectionauthor:: Edward Smith
!
    use coupler_module, only :  icmin_cnst, icmax_cnst, & 
                                jcmin_cnst, jcmax_cnst, & 
                                kcmin_cnst, kcmax_cnst

   integer, intent(out) :: limits(6)

   limits(1) = icmin_cnst
   limits(2) = icmax_cnst
   limits(3) = jcmin_cnst
   limits(4) = jcmax_cnst
   limits(5) = kcmin_cnst
   limits(6) = kcmax_cnst

end subroutine CPL_get_cnst_limits

!-----------------------------------------------------------------------------

subroutine CPL_get_bnry_limits(limits)
!
!Get limits of boundary region as cell indices in global coordinate system, 
!these are as specified in the input file (COUPLER.in) and must be in the
!overlap region.
!
!**Remarks**
!
!Assumes the coupler has been initialised with `CPL_init <#f/_/cpl_init>`_ and 
!topological mapping has been setup using either `CPL_setup_md <#f/_/cpl_setup_md>`_ 
!or `CPL_setup_cfd <#f/_/cpl_setup_cfd>`_ as appropriate.
! - Note: limits(6) are of the form: (xmin,xmax,ymin,ymax,zmin,zmax)
!
!**Synopsis**
!
!.. code-block:: fortran
!
!  CPL_get_bnry_limits(limits)
!
!**Outputs**
!
! - limits
!
!   - Array of cell extents that specify the boundary region. 
!
! .. sectionauthor:: Edward Smith
!
    use coupler_module, only :  icmin_bnry, icmax_bnry, & 
                                jcmin_bnry, jcmax_bnry, & 
                                kcmin_bnry, kcmax_bnry

   integer, intent(out) :: limits(6)

   limits(1) = icmin_bnry
   limits(2) = icmax_bnry
   limits(3) = jcmin_bnry
   limits(4) = jcmax_bnry
   limits(5) = kcmin_bnry
   limits(6) = kcmax_bnry

end subroutine CPL_get_bnry_limits

subroutine CPL_get_arrays(recv_array, recv_size, &
                          send_array, send_size)
!
!A helper function to get arrays of the required size for cells local to
!the current processor
!
!**Example**
!
!The first example shows the CFD side of the exchange, with send/recv arrays 
!obtained from CPL_get_arrays
!
!.. literalinclude:: ../../../examples/minimal_send_recv_mocks/minimal_CFD.f90
!
!The corresponding MD code which matches this, note CPL_get_arrays send and recv
!arrays are swapped over (as send from CFD in recv on MD and vice versa)
!
!.. literalinclude:: ../../../examples/minimal_send_recv_mocks/minimal_MD.f90
!
!These are then both run together, either MPMD or port connect modes
!
! .. sectionauthor:: Edward Smith
!
    integer :: recv_size, send_size
    integer :: limits(6), portion(6), Ncells(3)
    double precision, dimension(:,:,:,:), & 
         allocatable  :: send_array, recv_array

    !Get detail for grid
    call CPL_get_olap_limits(limits)
    call CPL_my_proc_portion(limits, portion)
    call CPL_get_no_cells(portion, Ncells)

    if (allocated(recv_array)) deallocate(recv_array)
    if (allocated(send_array)) deallocate(send_array)

    allocate(recv_array(recv_size, Ncells(1), Ncells(2), Ncells(3)))
    allocate(send_array(send_size, Ncells(1), Ncells(2), Ncells(3)))

end subroutine CPL_get_arrays

!-----------------------------------------------------------------------------

function CPL_is_proc_inside(region) result(res)
    use coupler_module, only: VOID

    integer, intent(in) :: region(6)
    logical :: res

    res = .true.
    if (any(region.eq.VOID)) then
        res = .false.
    end if
end function CPL_is_proc_inside

function is_cell_inside(cell, limits) result(res)
    use coupler_module, only :  VOID

    integer, intent(in) :: cell(3)
    integer, intent(in) :: limits(6)
    logical :: res
   
    res = .true.
    ! Check send limits are inside a certain region
    if (any(limits.eq.VOID)) then
        res =.false.
    end if
    if (cell(1) .lt. limits(1) .or. &
        cell(1) .gt. limits(2) .or. &
        cell(2) .lt. limits(3) .or. &
        cell(2) .gt. limits(4) .or. &
        cell(3) .lt. limits(5) .or. &
        cell(3) .gt. limits(6)) then

        res = .false.
    end if

end function is_cell_inside

!PRIVATE FUNCTION

function is_coord_inside(coord, coord_limits) result(res)

    real(kind(0.d0)), intent(in)   :: coord(3)
    real(kind(0.d0)), intent(in)   :: coord_limits(6)
    logical :: res

    res = .true.
    ! Check send limits are inside a certain region
    if (coord(1) .lt. coord_limits(1) .or. &
        coord(1) .gt. coord_limits(2) .or. &
        coord(2) .lt. coord_limits(3) .or. &
        coord(2) .gt. coord_limits(4) .or. &
        coord(3) .lt. coord_limits(5) .or. &
        coord(3) .gt. coord_limits(6)) then

        res = .false.
    end if
end function is_coord_inside


!-----------------------------------------------------------------------------
 
function coupler_md_get_save_period() result(p)
    use coupler_module, only : save_period
    implicit none

    integer p
    p = save_period

end function coupler_md_get_save_period

!----------------------------------------------------------------------------- 

function coupler_md_get_average_period() result(p)
    use coupler_module, only : average_period
    implicit none

    integer p
    p = average_period

end function coupler_md_get_average_period

!----------------------------------------------------------------------------- 

function coupler_md_get_md_steps_per_cfd_dt() result(p)
    use coupler_module, only : timestep_ratio 
    implicit none

    integer p
    p = timestep_ratio 

end function coupler_md_get_md_steps_per_cfd_dt

!-----------------------------------------------------------------------------

!-----------------------------------------------------------------------------

function CPL_cfd_dt() result(p)
    use coupler_module, only : dt_CFD  
    implicit none

    real(kind=kind(0.d0)) p
    p = dt_CFD

end function CPL_cfd_dt

!------------------------------------------------------------------------------

subroutine test_python (integer_p, double_p, bool_p, integer_pptr, double_pptr)
  integer, intent(in) :: integer_p
  real(kind(0.d0)), intent(in) :: double_p
  logical, intent(in) :: bool_p
  integer, dimension(:) :: integer_pptr
  real(kind(0.d0)), dimension(:) :: double_pptr

  print *, 'Integer: ', integer_p
  print *, 'Double: ', double_p
  print *, 'bool: ', bool_p
  print *, 'integer_pptr: ', integer_pptr
  print *, 'double_pptr: ', double_pptr
 end subroutine

!------------------------------------------------------------------------------

subroutine MPI_errorcheck(ierr)
    use mpi

    integer, intent(in) :: ierr

    integer             :: resultlen, newierr
    character(12)       :: err_buffer

    call MPI_Error_string(ierr, err_buffer, resultlen, newierr)
    print*, err_buffer

end subroutine MPI_errorcheck

end module coupler
