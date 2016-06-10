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
!           Copyright (C) 2012-2015 Edward Smith & David Trevelyan
!
!License
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
!
!Description
!
!
!Author(s)
! .. codeauthor:: Edward Smith 
! .. codeauthor:: David Trevelyan September 2012 to De
! .. codeauthor:: Lucian Anton, November 2011  
!
!! Routines accessible from application ( molecular or continuum ) after
!! the name, in parenthesis, is the realm in which each routine must be called
!!
! SIMULATION SIMULATION SIMULATION SIMULATION SIMULATION SIMULATION SIMULATION
!!
!! - CPL_send_data            (cfd+md)   sends grid data exchanged between
!!                                      realms ( generic interface)
!!
!! - CPL_recv_data            (cfd+md)   receives data exchanged between realms
!!                                      ( generic interface)
!!
!! - CPL_cfd_get               (cfd)    returns coupler internal parameters
!!                                      for CFD realm
!!
!! - CPL_md_get                 (md)    returns coupler internal parameters
!!                                      for MD realm
!!
!! - CPL_md_get_save_period     (md)    auxiliary used for testing
!!
!! - CPL_md_get_average_period  (md)    returns average period of BC
!!
!! - CPL_md_get_md_per_cfd_dt   (md)    returns the number of step MD does for 
!!                                      each CFD step
!!
!! - CPL_md_get_nsteps          (md)    returm CFD nsteps  
!!
!! - CPL_md_get_dt_cfd          (md)    returns MD dt
!!
!! - CPL_md_set                 (md)    sets zL if CFD is 2D
!!
!! - CPL_md_get_density         (md)    gets CFD density
!!
!! - CPL_md_get_cfd_id          (md)    id for CFD code, possible values set 
!!                                      in coupler_parameters
!!
!! @see coupler_module
!=============================================================================

module coupler
    implicit none

    private

    public CPL_send, CPL_recv, CPL_gather, CPL_scatter, CPL_get, CPL_proc_extents,&
           CPL_my_proc_extents, CPL_proc_portion, CPL_my_proc_portion, &
           CPL_map_cell2coord, CPL_map_coord2cell, CPL_get_no_cells, &
           CPL_map_glob2loc_cell, CPL_get_olap_limits, CPL_get_cnst_limits, & 
           CPL_map_cfd2md_coord, CPL_map_md2cfd_coord, CPL_overlap

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
!>
!! Perform gather operation on CPL_OLAP_COMM communicator. The CFD processor
!! is the root process. The gathered data is effectively "slotted" into the
!! correct part of the recvarray, and is intented for use in providing the
!! CFD simulation boundary conditions with data obtained from the MD
!! simulation.
!!
!! - Synopsis
!!
!!  - CPL_gather(gatherarray,npercell,limits,recvarray)
!!
!! - Input
!!
!!  - gatherarray
!!   - Assumed shape array of data to be gathered from each MD processor
!!     in the overlap communicator. Must be size 0 on CFD processor.
!!
!!  - limits
!!   - Integer array of length 6, specifying the global cell extents of the
!!     region to be gathered, is the same on ALL processors.
!!
!!  - npercell
!!   - number of data points per cell to be gathered (integer)
!!     Note: should be the same as size(gatherarray(1)) for MD
!!     processor. E.G. npercell = 3 for gathering 3D velocities.
!!
!! - Input/Output
!!  - recvarray
!!   - The array in which the gathered values are to be stored on the CFD
!!     processor. The only values to be changed in recvarray are:
!!     recvarray(limits(1):limits(2),limits(3):limits(4),limits(5):limits(6))
!!
!! - Output Parameters
!!  - NONE
!!
!! .. sectionauthor:: David Trevelyan

subroutine CPL_gather(gatherarray, npercell, limits, recvarray)!todo better name than recvarray
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
        integer :: coord(3),portion(6),mdextents(6)

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
!>
!! Scatter cell-wise data from CFD processor to corresponding MD processors
!! on the overlap communicator CPL_OLAP_COMM.
!!
!! - Synopsis
!!
!!  - CPL_scatter(scatterarray,npercell,limits,recvarray)
!!
!! - Input
!!
!!  - scatterarray
!!   - assumed shape array of data to be scattered (real(kind(0.d0)))
!!
!!  - limits
!!   - integer array of length 6, specifying the global cell extents of the
!!     region to be scattered, is the same on all processors.
!!
!!  - npercell
!!   - number of data points per cell to be scattered (integer).
!!     Note: should be the same as size(scatterarray(1)) for CFD proc
!!
!! - Input/Output
!!  - recvarray
!!   - the array in which the scattered values are stored on the MD
!!     processors.
!!
!! - Output
!!  - NONE
!!
!! .. sectionauthor:: David Trevelyan
subroutine CPL_scatter(scatterarray,npercell,limits,recvarray)
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
subroutine CPL_send(asend, limits, send_flag)
! ----------------------------------------------------------------------------
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
!.. code-block:: c
!
!  CPL_send(
!           asend, 
!           limits, 
!           send_flag
!           )    
!
!**Inputs**
!
! - *asend*
!
!   - Array of data to send. Should be a four dimensional array allocated using the number of cells on the current processor between the limits. Size should be be obtained from `CPL_my_proc_portion(limits, portion) <#f/_/cpl_my_proc_portion>`_.
! 
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
!.. code-block:: guess
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
    use coupler_module, only : md_realm,cfd_realm, & 
                               error_abort,CPL_GRAPH_COMM,myid_graph,olap_mask, &
                               rank_world, realm, rank_realm,rank_olap, & 
                               iblock_realm,jblock_realm,kblock_realm,ierr, VOID, &
							   CPL_setup_complete, REALM_NAME, realm
    implicit none

    
    logical, intent(out), optional  :: send_flag !Flag set if processor has passed data   
    integer, dimension(6), intent(in) :: limits ! Global cell indices with minimum and maximum values to send
    real(kind=kind(0.d0)),dimension(:,:,:,:), intent(in):: asend ! Array containing data to send
   
    !Neighbours
    integer                             :: nneighbors   
    integer,dimension(:),allocatable    :: id_neighbors

    ! local indices 
    integer :: icell,jcell,kcell
    integer :: n,pos,iclmin,iclmax,jclmin,jclmax,kclmin,kclmax
    integer :: npercell

    ! auxiliaries 
    integer                             :: nbr, ndata, itag, destid, Ncells
    integer,dimension(3)                :: pcoords, Ncell
    integer,dimension(6)                :: portion, myportion, portion_CFD
    real(kind=kind(0.d0)), allocatable  :: vbuf(:)

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

    !Get neighbours
    call MPI_Graph_neighbors_count(CPL_GRAPH_COMM,myid_graph,nneighbors,ierr)
    allocate(id_neighbors(nneighbors))
    call MPI_Graph_neighbors(CPL_GRAPH_COMM,myid_graph,nneighbors,id_neighbors,ierr )

    !Set sendflag to false and only change if anything is sent
    if (present(send_flag)) send_flag = .false.
  
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
            itag = 0 !mod( ncalls, MPI_TAG_UB) !Attention ncall could go over max tag value for long runs!!
            call MPI_sSend(vbuf, ndata, MPI_DOUBLE_PRECISION, destid, itag, CPL_GRAPH_COMM, ierr)

        endif

    enddo

end subroutine CPL_send


! ----------------------------------------------------------------------------
subroutine CPL_recv(arecv, limits, recv_flag)
! ----------------------------------------------------------------------------
!
! Receive data from to local grid from the associated ranks from the other 
! realm
!
!**Remarks**
!
!Assumes the coupler has been initialised with `CPL_init <#f/_/cpl_init>`_ and 
!topological mapping has been setup using either `CPL_setup_md <#f/_/cpl_setup_md>`_ 
!or `CPL_setup_cfd <#f/_/cpl_setup_cfd>`_ as appropriate.
!
!**Synopsis**
!
!.. code-block:: c
!
!  CPL_send(
!           arecv, 
!           limits, 
!           recv_flag
!           )    
!
!**Inputs**
!
! - *arecv*
!
!   - Array of data to recv. Should be a four dimensional array allocated using the number of cells on the current processor between the limits. Size should be be obtained from `CPL_my_proc_portion(limits, portion) <#f/_/cpl_my_proc_portion>`_.
! 
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
! ----------------------------------------------------------------------------
    use mpi
    use coupler_module, only : md_realm,cfd_realm, & 
                               rank_graph, &
                               error_abort,CPL_GRAPH_COMM,myid_graph,olap_mask, &
                               rank_world, realm, rank_realm, rank_olap, & 
                               iblock_realm,jblock_realm,kblock_realm,VOID,ierr, &
							   CPL_setup_complete, REALM_NAME, realm
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
    integer,dimension(6) :: portion, myportion, portion_CFD

    ! auxiliaries 
    integer :: itag, sourceid,start_address
    integer,dimension(:),allocatable   :: req
    integer,dimension(:,:),allocatable :: status
    real(kind(0.d0)),dimension(:), allocatable ::  vbuf

	!Check setup is complete
	if (CPL_setup_complete .ne. 1) then
		call error_abort("Error CPL_recv called before CPL_setup_"//REALM_NAME(realm))
	endif
 
    ! This local CFD domain is outside MD overlap zone 
    if (olap_mask(rank_world).eqv. .false.) return

    ! Number of components at each grid point
    npercell = size(arecv,1)

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

    !Get neighbours
    call MPI_Graph_neighbors_count(CPL_GRAPH_COMM,myid_graph,nneighbors,ierr)
    allocate(id_neighbors(nneighbors))
    call MPI_Graph_neighbors(CPL_GRAPH_COMM,myid_graph,nneighbors,id_neighbors,ierr )

    ! Receive from all attached processors
    allocate(req(nneighbors))
    allocate(status(MPI_STATUS_SIZE,nneighbors))
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
            pcoords = (/iblock_realm,jblock_realm,kblock_realm /)
        endif

        ! If limits passed to recv routine, use these instead
        ! of overlap/processor limits
        call CPL_proc_portion(pcoords, md_realm, limits, portion, ncells)

        !Only receive if overlapping
        if (any(portion.eq.VOID)) then
            ndata = 0
            req(nbr) = MPI_REQUEST_NULL
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
            itag = 0
            call MPI_irecv(vbuf(start_address),ndata,MPI_DOUBLE_PRECISION,sourceid,itag,&
                                    CPL_GRAPH_COMM,req(nbr),ierr)

        endif

        !Increment pointer ready to receive next piece of data      
        start_address = start_address + ndata

    enddo
    call MPI_waitall(nneighbors, req, status, ierr)
    call MPI_errorcheck(ierr)

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
        call error_abort("Error in CPL send -- x size of arecv must be the same as portion(2)-portion(1)+1")
    endif
    if (myportion(4)-myportion(3)+1 .ne. size(arecv,3)) then
        print*, myportion(4)-myportion(3)+1, size(arecv,3)
        call error_abort("Error in CPL send -- y size of arecv must be the same as portion(4)-portion(3)+1")
    endif
    if (myportion(6)-myportion(5)+1 .ne. size(arecv,4)) then
        print*, myportion(6)-myportion(5)+1, size(arecv,4)
        call error_abort("Error in CPL send -- z size of arecv must be the same as portion(6)-portion(5)+1")
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
           
end subroutine CPL_recv

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
!>
!!
!! Gets maximum and minimum cells for processor coordinates
!!
!! - Synopsis
!!
!!  - CPL_proc_extents(coord,realm,extents,ncells)
!!
!! - Input
!!
!!  - coord
!!   - processor cartesian coordinate (3 x integer) 
!!
!!  - realm
!!   - cfd_realm (1) or md_realm (2) (integer) 
!!
!! - Input/Output
!!  - NONE
!!
!! - Output
!!
!!  - extents
!!   - Six components array which defines processor extents
!!     xmin,xmax,ymin,ymax,zmin,zmax (6 x integer) 
!!
!!  - ncells (optional)
!!   - number of cells on processor (integer) 
!!
!! .. sectionauthor:: David Trevelyan
!! .. sectionauthor:: Edward Smith
subroutine CPL_proc_extents(coord, realm, extents, ncells)
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
	!	call error_abort("Error CPL_extents/portion called before CPL_setup_"//REALM_NAME(realm))
	!endif

    select case(realm)
    case(md_realm)
        if (size(icPmax_md) .lt. coord(1)) then
            write(strng,'(a,i8,a,i8,a)'), 'CPL_proc_extents error - MD proc ',       & 
                                           coord(1), ' extents requested but only ', & 
                                           size(icPmax_md) , ' in x'
            call error_abort(trim(strng))
        endif
        if (size(jcPmax_md) .lt. coord(2)) then
            write(strng,'(a,i8,a,i8,a)'), 'CPL_proc_extents error - MD proc ',       & 
                                           coord(2), ' extents requested but only ', & 
                                           size(icPmax_md) , ' in y'
            call error_abort(trim(strng))
        endif
        if (size(kcPmax_md) .lt. coord(3)) then
            write(strng,'(a,i8,a,i8,a)'), 'CPL_proc_extents error - MD proc ',       & 
                                           coord(3), ' extents requested but only ', & 
                                           size(icPmax_md) , ' in z'
            call error_abort(trim(strng))
        endif
        extents = (/icPmin_md(coord(1)),icPmax_md(coord(1)), & 
                    jcPmin_md(coord(2)),jcPmax_md(coord(2)), & 
                    kcPmin_md(coord(3)),kcPmax_md(coord(3))/)
    case(cfd_realm)
        if (size(icPmax_cfd) .lt. coord(1)) then
            write(strng,'(a,i8,a,i8,a)'), 'CPL_proc_extents error - CFD proc ',      & 
                                           coord(1), ' extents requested but only ', & 
                                           size(icPmax_cfd) , ' in x'
            call error_abort(trim(strng))
        endif
        if (size(jcPmax_cfd) .lt. coord(2)) then
            write(strng,'(a,i8,a,i8,a)'), 'CPL_proc_extents error - CFD proc ',      & 
                                           coord(2), ' extents requested but only ', & 
                                           size(icPmax_cfd) , ' in y'
            call error_abort(trim(strng))
        endif
        if (size(kcPmax_cfd) .lt. coord(3)) then
            write(strng,'(a,i8,a,i8,a)'), 'CPL_proc_extents error - CFD proc ',      & 
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
!>
!!
!! Get maximum and minimum cells for current communicator within
!! the overlapping region only
!!
!! - Synopsis 
!!
!!  - CPL_olap_extents(coord,realm,extents,ncells)
!!
!! - Input 
!!
!!  - coord
!!   - processor cartesian coordinate (3 x integer) 
!!
!!  - realm
!!   - cfd_realm (1) or md_realm (2) (integer) 
!!
!! - Input/Output
!!  - NONE
!!
!! - Output 
!!
!!  - extents
!!   - Six components array which defines processor extents within
!!     the overlap region only: xmin,xmax,ymin,ymax,zmin,zmax (6 x integer) 
!!
!!  - ncells (optional)
!!   - number of cells on processor (integer) 
!!
!! .. sectionauthor:: David Trevelyan
subroutine CPL_olap_extents(coord,realm,extents,ncells)
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
!.. code-block:: c
!
!  CPL_proc_portion(
!                   coord,
!                   realm,
!                   limits,
!                   portion,
!                   ncells
!                   )
!
!**Inputs**
!
! - *coord*
!
!   - processor cartesian coordinate (3 x integer) 
! - *realm*
!
!   - cfd_realm (1) or md_realm (2) (integer) 
! - *limits*
!
!   - Array of cell extents that specify the input region. 
!
!
!**Outputs**
!
! - *portion*
!   - Array of cell extents that define the local processor's
!     contribution to the input region 'limits'.
!
! - *ncells (optional)*
!    - number of cells in portion (integer) 
!
!
! .. sectionauthor:: David Trevelyan

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
!.. code-block:: c
!
!  CPL_proc_portion(
!                   limits,
!                   portion,
!                   )
!
!**Inputs**
!
! - *limits*
!
!   - Array of cell extents that specify the input region. 
!
!**Outputs**
!
! - *portion*
!   - Array of cell extents that define the local processor's
!     contribution to the input region 'limits'.
!
! .. sectionauthor:: Eduardo Ramos Fernandez
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
! Determines process coords in appropriate realm's cartesian topology 
! given a rank in any communicator
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
!.. code-block:: c
!
!  CPL_Cart_coords(
!                  COMM, 
!                  rank, 
!                  realm, 
!                  maxdims, 
!                  coords, 
!                  ierr
!                  )
!
!**Inputs**
!
! - *comm*
!
!   - communicator with cartesian structure (handle) 
! - *realm*
!
!   - cfd_realm (1) or md_realm (2) (integer) 
! - *rank*
!
!   - rank of a process within group of comm (integer) 
!      NOTE fortran convention rank=1 to nproc
! - *maxdims*
!
!   - length of vector coords in the calling program (integer) 
!
!**Outputs**
!
! - *coords*
!
!   - integer array (of size ndims) containing the Cartesian coordinates 
!     of specified process (integer) 
! - *ierr*
!
!   - error flag
! .. sectionauthor:: Edward Smith
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

    integer, intent(in)     :: COMM, realm, rank, maxdims
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
! Return rank of current processor in specified COMM 
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
!.. code-block:: c
!
!  CPL_get_rank(
!               COMM, 
!               rank
!               )
!
!**Inputs**
!
! - *comm*
!   - communicator with cartesian structure (handle) 
!
!**Outputs**
!
! - *rank*
!   - rank of a process within group of comm (integer) 
!      NOTE fortran convention rank=1 to nproc
!
! .. sectionauthor:: Edward Smith

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
!>
!! Check if current processor is in the overlap region
!!
!! - Synopsis
!!
!!  - CPL_olap_check()
!!
!! - Input Parameters
!!
!!  - NONE
!!
!! - Returns
!!
!!  - CPL_olap_check
!!   - True if calling processor is in the overlap region
!!     and false otherwise
!!
!! .. sectionauthor:: Edward Smith


function CPL_overlap() result(p)
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
                   cpl_cfd_bc_x, cpl_cfd_bc_y, cpl_cfd_bc_z,     &
                   timestep_ratio, comm_style)
! Wrapper to retrieve (read only) parameters from the coupler_module 
! Note - this ensures all variable in the coupler are protected
! from corruption by either CFD or MD codes
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
                                cpl_cfd_bc_x_ => cpl_cfd_bc_x, &
                                cpl_cfd_bc_y_ => cpl_cfd_bc_y, &
                                cpl_cfd_bc_z_ => cpl_cfd_bc_z, &
                                comm_style_ => comm_style, &
                                comm_style_send_recv_ => comm_style_send_recv, &
                                comm_style_gath_scat_ => comm_style_gath_scat
    implicit none

    logical,dimension(3),optional,intent(out) :: staggered_averages

    integer, optional, intent(out)          :: icmax_olap ,icmin_olap
    integer, optional, intent(out)          :: jcmax_olap ,jcmin_olap
    integer, optional, intent(out)          :: kcmax_olap ,kcmin_olap
    integer, optional, intent(out)          :: ncx,ncy,ncz
    integer, optional, intent(out)          :: npx_md,npy_md,npz_md
    integer, optional, intent(out)          :: npx_cfd,npy_cfd,npz_cfd
    integer, optional, intent(out)          :: md_cfd_match_cellsize,timestep_ratio

    integer, optional, intent(out)          :: constraint_algo
    integer, optional, intent(out)          :: constraint_CVflag
    integer, optional, intent(out)          :: constraint_OT
    integer, optional, intent(out)          :: constraint_NCER
    integer, optional, intent(out)          :: constraint_Flekkoy
    integer, optional, intent(out)          :: constraint_CV
    integer, optional, intent(out)          :: constraint_off
    integer, optional, intent(out)          :: comm_style 
    integer, optional, intent(out)          :: icmax_cnst, icmin_cnst
    integer, optional, intent(out)          :: jcmax_cnst, jcmin_cnst
    integer, optional, intent(out)          :: kcmax_cnst, kcmin_cnst
    integer, optional, intent(out)          :: cpl_cfd_bc_slice 
    integer, optional, intent(out)          :: cpl_cfd_bc_x 
    integer, optional, intent(out)          :: cpl_cfd_bc_y 
    integer, optional, intent(out)          :: cpl_cfd_bc_z 
    integer, optional, intent(out)          :: cpl_md_bc_slice 

    real(kind(0.d0)), optional, intent(out) :: density_cfd,density_md
    real(kind(0.d0)), optional, intent(out) :: dt_cfd,dt_MD
    real(kind(0.d0)), optional, intent(out) :: dx,dy,dz
    real(kind(0.d0)), optional, intent(out) :: xL_md,xL_cfd
    real(kind(0.d0)), optional, intent(out) :: yL_md,yL_cfd
    real(kind(0.d0)), optional, intent(out) :: zL_md,zL_cfd

    real(kind(0.d0)), dimension(:,:),allocatable,optional,intent(out) :: xg,yg
    real(kind(0.d0)), dimension(:)  ,allocatable,optional,intent(out) :: zg

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
        allocate(xg(size(xg_,1),size(xg_,2)))
        xg = xg_
    endif
    if (present(yg)) then
        allocate(yg(size(yg_,1),size(yg_,2)))
        yg = yg_
    endif
    if (present(zg)) then
        allocate(zg(size(zg_,1)))
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

end function


!=============================================================================
!> Map global MD position to global CFD coordinate frame
!-----------------------------------------------------------------------------
function CPL_map_md2cfd_coord(coord_md, coord_cfd) result(valid_coord)
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

    valid_coord = is_coord_inside(coord_md, coord_limits_md)
    
    if (valid_coord) then
        !Get size of MD domain which has no CFD cells overlapping
        !This should be general enough to include grid stretching
        !and total overlap in any directions 
        md_only(1) = xL_md-abs(xg(icmax_olap+1,1) - xg(icmin_olap,1))
        md_only(2) = yL_md-abs(yg(1,jcmax_olap+1) - yg(1,jcmin_olap))
        md_only(3) = zL_md-abs(zg( kcmax_olap+1 ) - zg( kcmin_olap ))

        ! CFD has origin at bottom left while MD origin at centre
        !coord_md = coord_cfd + md_only + md_xyz_orig
        coord_cfd(1) = abs(coord_md(1) - x_orig_md - md_only(1)) + x_orig_cfd
        coord_cfd(2) = abs(coord_md(2) - y_orig_md - md_only(2)) + y_orig_cfd
        coord_cfd(3) = abs(coord_md(3) - z_orig_md - md_only(3)) + z_orig_cfd
        valid_coord = is_coord_inside(coord_cfd, coord_limits_cfd)
    endif

end function CPL_map_md2cfd_coord



!=============================================================================
!> Map global CFD position in global MD coordinate frame
!-----------------------------------------------------------------------------
function CPL_map_cfd2md_coord(coord_cfd, coord_md) result(valid_coord)
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
        md_only(1) = xL_md-abs(xg(icmax_olap+1,1) - xg(icmin_olap,1))
        md_only(2) = yL_md-abs(yg(1,jcmax_olap+1) - yg(1,jcmin_olap))
        md_only(3) = zL_md-abs(zg( kcmax_olap+1 ) - zg( kcmin_olap ))

        ! CFD has origin at bottom left while MD origin at centre
        !coord_md = coord_cfd + md_only + md_xyz_orig
        coord_md(1) = abs(coord_cfd(1) - x_orig_cfd) + md_only(1) + x_orig_md
        coord_md(2) = abs(coord_cfd(2) - y_orig_cfd) + md_only(2) + y_orig_md
        coord_md(3) = abs(coord_cfd(3) - z_orig_cfd) + md_only(3) + z_orig_md
        valid_coord = is_coord_inside(coord_md, coord_limits_md)
    endif

end function CPL_map_cfd2md_coord

!-----------------------------------------------------------------------------

subroutine CPL_map_cell2coord(i, j, k, coord_xyz)

    use coupler_module, only: xg, yg, zg, realm, &
                              md_realm, cfd_realm, error_abort
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

    coord_xyz(1) = xg(i, j)
    coord_xyz(2) = yg(i, j) 
    coord_xyz(3) = zg(k)

    if (realm .eq. md_realm) then
        aux_ret = CPL_map_cfd2md_coord(coord_xyz, coord_md)
        coord_xyz = coord_md
    end if


end subroutine CPL_map_cell2coord

!-----------------------------------------------------------------------------

function CPL_map_coord2cell(x, y, z, cell_ijk) result(ret)

    use coupler_module, only: dx, dy, dz, &
                              icmin_olap, jcmin_olap, kcmin_olap

    real(kind(0.d0)), intent(in)  :: x, y, z
    integer, intent(out)         :: cell_ijk(3)

    real(kind(0.d0)) :: olap_lo(3)
    integer          :: olap_limits(6)
    logical          :: ret

    call CPL_map_cell2coord(icmin_olap, jcmin_olap, kcmin_olap, olap_lo)

    cell_ijk(1) = int(x - olap_lo(1)) / dx + 1
    cell_ijk(2) = int(y - olap_lo(2)) / dy + 1
    cell_ijk(3) = int(z - olap_lo(3)) / dz + 1

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
    integer :: olap_limits(6)

    call CPL_get_olap_limits(olap_limits)

    ! Check if cell is inside the overlap region
    if (.not. is_cell_inside(glob_cell, olap_limits)) then 
        print*, "cell:" , glob_cell
        call error_abort("CPL_map_glob2loc_cell error - Cell is outside overlap region. " // &
                         "Aborting simulation.") 
    end if

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

!-----------------------------------------------------------------------------

subroutine CPL_get_cnst_limits(limits)


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

function is_cell_inside(cell, limits) result(res)
    use coupler_module, only :  icmin_olap, icmax_olap, & 
                                jcmin_olap, jcmax_olap, & 
                                kcmin_olap, kcmax_olap

    integer, intent(in) :: cell(3)
    integer, intent(in) :: limits(6)
    logical :: res
    
    res = .true.
    ! Check send limits are inside overlap region
    if (cell(1) .lt. limits(1) .or. &
        cell(1) .gt. limits(2) .or. &
        cell(2) .lt. limits(3) .or. &
        cell(2) .gt. limits(4) .or. &
        cell(3) .lt. limits(5) .or. &
        cell(3) .gt. limits(6)) then

        res = .false.
    end if
end function is_cell_inside

function is_coord_inside(coord, coord_limits) result(res)
    use coupler_module, only :  icmin_olap, icmax_olap, & 
                                jcmin_olap, jcmax_olap, & 
                                kcmin_olap, kcmax_olap, &
                                dx, dy, dz

    real(kind(0.d0)), intent(in)   :: coord(3)
    real(kind(0.d0)), intent(in)   :: coord_limits(6)
    logical :: res

    res = .true.
    ! Check send limits are inside overlap region
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

function coupler_md_get_nsteps() result(p)
    use coupler_module, only :  nsteps_md
    implicit none 

     integer p
     p = nsteps_md

end function coupler_md_get_nsteps

!-----------------------------------------------------------------------------

function coupler_md_get_dt_cfd() result(p)
    use coupler_module, only : dt_CFD  
    implicit none

    real(kind=kind(0.d0)) p
    p = dt_CFD

end function coupler_md_get_dt_cfd

!------------------------------------------------------------------------------

subroutine test_python (integer_p, double_p, bool_p, integer_pptr, double_pptr)
  integer, intent(in) :: integer_p
  real(kind(0.d0)), intent(in) :: double_p
  logical, intent(in) :: bool_p
  integer, dimension(:) :: integer_pptr
  real(kind(0.d0)), dimension(:) :: double_pptr

  print *, 'Integer: ', integer_p
  print *, 'Double: ', double_p
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
