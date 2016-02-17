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
!
!   Edward Smith
!   David Trevelyan
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
!! @author  Lucian Anton, November 2011  
!! @author Edward Smith, Dave Trevelyan September 2012
!! @see coupler_module
!=============================================================================

module coupler
    implicit none

    interface CPL_send
        module procedure CPL_send_3d, CPL_send_4d
    end interface

    interface CPL_recv
        module procedure CPL_recv_3d, CPL_recv_4d
    end interface

    private CPL_send_3d, CPL_send_4d, &
        CPL_send_xd, CPL_recv_3d, CPL_recv_4d,&
        CPL_recv_xd

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
!! @author David Trevelyan
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
            limits(6) .lt. kcmax) then
            
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
        call CPL_proc_extents(coord,realm,mdextents)
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
            i = icell - mdextents(1) + 1
            j = jcell - mdextents(3) + 1
            k = kcell - mdextents(5) + 1
            sendbuf(pos) = gatherarray(ixyz,i,j,k)
            pos = pos + 1
        end do
        end do
        end do
        end do
    
    end subroutine pack_sendbuf
    
    subroutine unpack_recvbuf
        implicit none

        integer :: coord(3),portion(6),cfdextents(6)
        integer :: trank_olap, tid_olap
        integer :: pos,ixyz,icell,jcell,kcell
        integer :: i,j,k

        integer :: tempextents(6)

        ! Get CFD proc coords and extents, allocate suitable array
        call CPL_cart_coords(CPL_OLAP_COMM,rank_olap,cfd_realm,3,coord,ierr)
        call CPL_proc_extents(coord,cfd_realm,cfdextents)
        call CPL_proc_portion(coord,cfd_realm,limits,portion)

        ! Loop over all processors in overlap comm
        do trank_olap = 1,nproc_olap

            tid_olap = trank_olap - 1
            if (tid_olap .eq. CFDid_olap) cycle

            call CPL_Cart_coords(CPL_OLAP_COMM,trank_olap,md_realm,3,coord,ierr)
            call CPL_proc_portion(coord,md_realm,limits,portion)

            call CPL_proc_extents(coord,md_realm,tempextents)
            
            if (any(portion.eq.VOID)) cycle

            ! Set position and unpack MD proc's part of recvbuf to
            ! correct region of recvarray   
            pos = displs(trank_olap) + 1    
            do ixyz = 1,npercell
            do icell = portion(1),portion(2)
            do jcell = portion(3),portion(4)
            do kcell = portion(5),portion(6)

                i = icell - cfdextents(1) + 1 
                j = jcell - cfdextents(3) + 1 
                k = kcell - cfdextents(5) + 1 

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
!! @author David Trevelyan
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
            limits(6) .lt. kcmax) then
            
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
        integer :: coord(3),cfdextents(6),portion(6)
        integer :: ixyz, icell, jcell, kcell
        integer :: i,j,k

        ! Grab CFD proc extents to be used for local cell mapping (when an 
        ! array is an input to subroutine it has lower bound 1 by default)
        call CPL_cart_coords(CPL_CART_COMM,rank_cart,cfd_realm,3,coord,ierr)
        call CPL_proc_extents(coord,cfd_realm,cfdextents)

        ! Loop over procs in olap comm and pack scatter buffer 
        ! in separate regions for each MD proc
        pos = 1
        do n = 1,nproc_olap

            tid_olap = n - 1
            if (tid_olap.eq.CFDid_olap) cycle

            call CPL_Cart_coords(CPL_OLAP_COMM,n,md_realm,3,coord,ierr)
            call CPL_proc_portion(coord,md_realm,limits,portion)
            if (any(portion.eq.VOID)) cycle

            do ixyz = 1,npercell
            do icell= portion(1),portion(2)
            do jcell= portion(3),portion(4)
            do kcell= portion(5),portion(6)
                i = icell - cfdextents(1) + 1
                j = jcell - cfdextents(3) + 1
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
        call CPL_proc_extents(coord,realm,extents)
        call CPL_proc_portion(coord,realm,limits,portion)
        if (any(portion.eq.VOID)) return

        pos = 1
        do ixyz = 1,npercell
        do icell= portion(1),portion(2)
        do jcell= portion(3),portion(4)
        do kcell= portion(5),portion(6)
            i = icell - extents(1) + 1
            j = jcell - extents(3) + 1
            k = kcell - extents(5) + 1
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


!=============================================================================
!>
!! CPL_send_data wrapper for 3d arrays
!! see CPL_send_xd for input description
!! @see coupler#subroutine_CPL_send_xd
!-----------------------------------------------------------------------------
subroutine CPL_send_3d(temp,icmin_send,icmax_send,jcmin_send, & 
                            jcmax_send,kcmin_send,kcmax_send,send_flag)
    use coupler_module, only : icmin_olap,icmax_olap, & 
                               jcmin_olap,jcmax_olap, &
                               kcmin_olap,kcmax_olap,error_abort
    implicit none
 
    logical, intent(out), optional                      :: send_flag
    integer, intent(in), optional                       :: icmax_send,icmin_send
    integer, intent(in), optional                       :: jcmax_send,jcmin_send
    integer, intent(in), optional                       :: kcmax_send,kcmin_send
    real(kind=kind(0.d0)),dimension(:,:,:), intent(in)  :: temp
    
    integer :: n1,n2,n3,n4
    integer :: icmin,icmax,jcmin,jcmax,kcmin,kcmax
    real(kind=kind(0.d0)),dimension(:,:,:,:),allocatable :: asend


    !Revert to default i domain sending - top of overlap to bottom of overlap
    if ((present(icmax_send)) .and. (present(icmin_send))) then
            icmax = icmax_send; icmin = icmin_send
    elseif ((.not. present(icmax_send)).and.(.not. present(icmin_send))) then
            icmax = icmax_olap; icmin = icmin_olap
    else
        call error_abort("CPL_send error - both max and min i limits " // &
                         "required and only one supplied")
    endif

    !Revert to default j domain sending - top of overlap to bottom of overlap
    if ((present(jcmax_send)) .and. (present(jcmin_send))) then
            jcmax = jcmax_send; jcmin = jcmin_send
    elseif ((.not. present(jcmax_send)).and.(.not. present(jcmin_send))) then
            jcmax = jcmax_olap; jcmin = jcmin_olap
    else
        call error_abort("CPL_send error - both max and min j limits " // &
                         "required and only one supplied")
    endif

    !Revert to default k domain sending - top of overlap to bottom of overlap
    if ((present(kcmax_send)) .and. (present(kcmin_send))) then
            kcmax = kcmax_send; kcmin = kcmin_send
    elseif ((.not. present(kcmax_send)).and.(.not. present(kcmin_send))) then
            kcmax = kcmax_olap; kcmin = kcmin_olap
    else
        call error_abort("CPL_send error - both max and min k limits " // &
                         "required and only one supplied")
    endif

   
    n1 = 1
    n2 = size(temp,1)
    n3 = size(temp,2)
    n4 = size(temp,3)

    !Add padding column to 3D array to make it 4D
    allocate(asend(n1,n2,n3,n4))
    asend(1,:,:,:) = temp(:,:,:)

    call CPL_send_xd(asend,icmax,icmin,jcmax, & 
                           jcmin,kcmax,kcmin,send_flag )

end subroutine CPL_send_3d

!=============================================================================
!>
!! CPL_send_data wrapper for 4d arrays
!! see CPL_send_xd for input description
!! @see coupler#subroutine_CPL_send_xd
!-----------------------------------------------------------------------------
subroutine CPL_send_4d(asend,icmin_send,icmax_send,jcmin_send, & 
                             jcmax_send,kcmin_send,kcmax_send,send_flag)
    use coupler_module, only : icmin_olap,icmax_olap, & 
                               jcmin_olap,jcmax_olap, &
                               kcmin_olap,kcmax_olap,error_abort
    implicit none
 
    logical, intent(out), optional                      :: send_flag
    integer, intent(in), optional                       :: icmax_send,icmin_send
    integer, intent(in), optional                       :: jcmax_send,jcmin_send
    integer, intent(in), optional                       :: kcmax_send,kcmin_send
    real(kind=kind(0.d0)),dimension(:,:,:,:), intent(in) :: asend
    
    integer :: npercell
    integer :: icmin,icmax,jcmin,jcmax,kcmin,kcmax

    npercell = size(asend,1)
    !if ((present(icmax_send))) print*, 'icmax_send', icmax_send
    !if ((present(icmin_send))) print*, 'icmin_send', icmin_send
    !if ((present(jcmax_send))) print*, 'jcmax_send', jcmax_send
    !if ((present(jcmin_send))) print*, 'jcmin_send', jcmin_send
    !if ((present(kcmax_send))) print*, 'kcmax_send', kcmax_send
    !if ((present(kcmin_send))) print*, 'kcmin_send', kcmin_send

    !Revert to default i domain sending - top of overlap to bottom of overlap
    if ((present(icmax_send)) .and. (present(icmin_send))) then
            icmax = icmax_send; icmin = icmin_send
    elseif ((.not. present(icmax_send)).and.(.not. present(icmin_send))) then
            icmax = icmax_olap; icmin = icmin_olap
    else
        call error_abort("CPL_send error - both max and min i limits " // &
                         "required and only one supplied")
    endif

    !Revert to default j domain sending - top of overlap to bottom of overlap
    if ((present(jcmax_send)) .and. (present(jcmin_send))) then
            jcmax = jcmax_send; jcmin = jcmin_send
    elseif ((.not. present(jcmax_send)).and.(.not. present(jcmin_send))) then
            jcmax = jcmax_olap; jcmin = jcmin_olap
    else
        call error_abort("CPL_send error - both max and min j limits " // &
                         "required and only one supplied")
    endif

    !Revert to default k domain sending - top of overlap to bottom of overlap
    if ((present(kcmax_send)) .and. (present(kcmin_send))) then
            kcmax = kcmax_send; kcmin = kcmin_send
    elseif ((.not. present(kcmax_send)).and.(.not. present(kcmin_send))) then
            kcmax = kcmax_olap; kcmin = kcmin_olap
    else
        call error_abort("CPL_send error - both max and min k limits " // &
                         "required and only one supplied")
    endif

    !send_extents = (/npercell,icmax-icmin+1,jcmax-jcmin+1,kcmax-kcmin+1 /) 
    !if (any(shape(asend) .lt. send_extents)) then
    !   print'(2(a,4i5))', '  Shape of input array = ', shape(asend), & 
    !                     '  Passed range = ',send_extents 
    !   call error_abort("CPL_send error - Specified send range is greater" // &
    !                    "than the number of cells on processor")
    !endif
 
    call CPL_send_xd(asend,icmin,icmax,jcmin, & 
                           jcmax,kcmin,kcmax,send_flag )

end subroutine CPL_send_4d

!=============================================================================
!                       CPL_send_xd
!>
!! Send data from the local grid to the associated ranks from the other 
!! realm
!!
!! - Synopsis
!!
!!  - CPL_send_xd(asend,icmin_send,icmax_send,jcmin_send,  
!!                           jcmax_send,kcmin_send,kcmax_send,send_flag)
!!
!! - Input Parameters
!!
!!   - asend
!!
!!   - icmin_send
!!
!!   - icmax_send
!!
!!   - jcmin_send
!!
!!   - jcmax_send
!!
!!   - kcmin_send
!!
!!   - kcmax_send
!!
!! - Output Parameter
!!
!!   - send_flag
!!
!! @author Edward Smith
! ----------------------------------------------------------------------------
subroutine CPL_send_xd(asend,icmin_send,icmax_send,jcmin_send, & 
                             jcmax_send,kcmin_send,kcmax_send,send_flag)
    use mpi
    use coupler_module, only : md_realm,cfd_realm, & 
                               error_abort,CPL_GRAPH_COMM,myid_graph,olap_mask, &
                               rank_world, realm, & 
                               iblock_realm,jblock_realm,kblock_realm,ierr, VOID
    implicit none

    !Flag set if processor has passed data
    logical, intent(out), optional  :: send_flag

    ! Minimum and maximum values to send
    integer, intent(in) :: icmax_send,icmin_send, & 
                           jcmax_send,jcmin_send, & 
                           kcmax_send,kcmin_send

   ! Array containing data distributed on the grid
    real(kind=kind(0.d0)),dimension(:,:,:,:), intent(in):: asend
   
    !Number of halos
    !integer :: nh = 1 !?? todo needed?

    !Neighbours
    integer                             :: nneighbors   
    integer,dimension(:),allocatable    :: id_neighbors

    ! local indices 
    integer :: icell,jcell,kcell
    integer :: n,pos,iclmin,iclmax,jclmin,jclmax,kclmin,kclmax
    integer :: npercell

    ! auxiliaries 
    integer                             :: nbr,ndata,itag,destid,ncells
    integer                             :: pcoords(3)
    integer,dimension(6)                :: extents,limits,portion
    real(kind=kind(0.d0)), allocatable  :: vbuf(:)

    ! This local CFD domain is outside MD overlap zone 
    if (olap_mask(rank_world) .eqv. .false.) return

    ! Save limits array of Minimum and maximum values to send
    limits = (/ icmin_send,icmax_send,jcmin_send,jcmax_send,kcmin_send,kcmax_send /)

    ! Get local grid box ranges seen by this rank for CFD
    if (realm .eq. cfd_realm) then 
        !Load extents of CFD processor
        pcoords = (/iblock_realm,jblock_realm,kblock_realm /)
        call CPL_proc_extents(pcoords,cfd_realm,extents)
    endif

    ! Number of components at each grid point
    npercell = size(asend,1)

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
            !Get extents of nbr MD processor to send to
            call CPL_Cart_coords(CPL_GRAPH_COMM, destid+1,  md_realm, 3, pcoords, ierr) 
        elseif (realm .eq. md_realm) then
            !Data to send is based on current processor as MD proc size < CFD proc size
            pcoords = (/iblock_realm,jblock_realm,kblock_realm /)
            !Get extents of current processor
            call CPL_proc_extents(pcoords,md_realm,extents)
        endif

        ! If limits passed to send routine, use these instead
        ! of overlap/processor limits
        call CPL_proc_portion(pcoords,md_realm,limits,portion,ncells)

        ! Amount of data to be sent
        if (any(portion.eq.VOID)) then
            !print*, 'VOID send qqqq',realm_name(realm),rank_world,rank_realm
            ndata = 0
        else

            ! Get data range on processor's local extents
            iclmin = portion(1)-extents(1)+1;   iclmax = portion(2)-extents(1)+1
            jclmin = portion(3)-extents(3)+1;   jclmax = portion(4)-extents(3)+1
            kclmin = portion(5)-extents(5)+1;   kclmax = portion(6)-extents(5)+1

            ndata = npercell * ncells
            if (allocated(vbuf)) deallocate(vbuf); allocate(vbuf(ndata))
            if (present(send_flag)) send_flag = .true.
            ! Pack array into buffer
            pos = 1
            !print'(a,5i4,2i6,i4,24i4)', 'send qqqq',rank_world,rank_realm,rank_olap,ndata,nbr,destid, & 
            !                       size(asend),pos,&
            !                       iclmin,   iclmax,   jclmin,   jclmax,   kclmin,   kclmax,     &
            !                       icmin_send,icmax_send,jcmin_send,jcmax_send,kcmin_send,kcmax_send, & 
            !                       portion, extents
            do kcell=kclmin,kclmax
            do jcell=jclmin,jclmax
            do icell=iclmin,iclmax
            do n = 1,npercell
                vbuf(pos) = asend(n,icell,jcell,kcell)
                !write(98000+destid+1+10*rank_world,'(3i8,f20.5)') rank_world,destid+1,n, vbuf(pos)
                pos = pos + 1
            end do
            end do
            end do
            end do
            ! ----------------- pack data for destid -----------------------------

            ! Send data 
            itag = 0 !mod( ncalls, MPI_TAG_UB) !Attention ncall could go over max tag value for long runs!!
            call MPI_send(vbuf, ndata, MPI_DOUBLE_PRECISION, destid, itag, CPL_GRAPH_COMM, ierr)

        endif

    enddo

end subroutine CPL_send_xd

!=============================================================================
!>
!! CPL_recv_xd wrapper for 3d arrays
!! see CPL_recv_xd for input description
!! @see coupler#subroutine_CPL_recv_xd
!-----------------------------------------------------------------------------
subroutine CPL_recv_3d(temp,icmin_recv,icmax_recv,jcmin_recv, & 
                            jcmax_recv,kcmin_recv,kcmax_recv,recv_flag)
    use coupler_module, only : icmin_olap,icmax_olap, & 
                               jcmin_olap,jcmax_olap, &
                               kcmin_olap,kcmax_olap,error_abort
    implicit none

    logical, intent(out), optional                  :: recv_flag
    integer, intent(in), optional                   :: icmax_recv,icmin_recv
    integer, intent(in), optional                   :: jcmax_recv,jcmin_recv
    integer, intent(in), optional                   :: kcmax_recv,kcmin_recv
    real(kind(0.d0)),dimension(:,:,:),intent(inout) :: temp 
                                                          
    integer :: n1,n2,n3,n4
    integer :: icmin,icmax,jcmin,jcmax,kcmin,kcmax
    real(kind(0.d0)),dimension(:,:,:,:),allocatable  :: arecv

    !if ((present(icmax_recv))) print*, 'icmax_recv', icmax_recv
    !if ((present(icmin_recv))) print*, 'icmin_recv', icmin_recv
    !if ((present(jcmax_recv))) print*, 'jcmax_recv', jcmax_recv
    !if ((present(jcmin_recv))) print*, 'jcmin_recv', jcmin_recv
    !if ((present(kcmax_recv))) print*, 'kcmax_recv', kcmax_recv
    !if ((present(kcmin_recv))) print*, 'kcmin_recv', kcmin_recv

    !Revert to default i domain sending - top of overlap to bottom of overlap
    if ((present(icmax_recv)) .and. (present(icmin_recv))) then
            icmax = icmax_recv; icmin = icmin_recv
    elseif ((.not. present(icmax_recv)).and.(.not. present(icmin_recv))) then
            icmax = icmax_olap; icmin = icmin_olap
    else
        call error_abort("CPL_recv error - both max and min i limits " // &
                         "required and only one supplied")
    endif

    !Revert to default j domain sending - top of overlap to bottom of overlap
    if ((present(jcmax_recv)) .and. (present(jcmin_recv))) then
            jcmax = jcmax_recv; jcmin = jcmin_recv
    elseif ((.not. present(jcmax_recv)).and.(.not. present(jcmin_recv))) then
            jcmax = jcmax_olap; jcmin = jcmin_olap
    else
        call error_abort("CPL_recv error - both max and min j limits " // &
                         "required and only one supplied")
    endif

    !Revert to default k domain sending - top of overlap to bottom of overlap
    if ((present(kcmax_recv)) .and. (present(kcmin_recv))) then
            kcmax = kcmax_recv; kcmin = kcmin_recv
    elseif ((.not. present(kcmax_recv)).and.(.not. present(kcmin_recv))) then
            kcmax = kcmax_olap; kcmin = kcmin_olap
    else
        call error_abort("CPL_recv error - both max and min k limits " // &
                         "required and only one supplied")
    endif
 
    n1 = 1 
    n2 = size(temp,1)
    n3 = size(temp,2)
    n4 = size(temp,3)

    !Add padding column to 3D array to make it 4D
    allocate(arecv(n1,n2,n3,n4))
    call CPL_recv_xd(arecv,icmin,icmax,jcmin, & 
                           jcmax,kcmin,kcmax,recv_flag)
    temp(:,:,:) =   arecv(1,:,:,:) 

end subroutine CPL_recv_3d

!=============================================================================
!>
!! CPL_recv_xd  wrapper for 4d arrays
!! See CPL_recv_xd for input description
!! @see coupler#subroutine_CPL_recv_xd
!-----------------------------------------------------------------------------
subroutine CPL_recv_4d(arecv,icmin_recv,icmax_recv,jcmin_recv, & 
                             jcmax_recv,kcmin_recv,kcmax_recv,recv_flag)
    use coupler_module, only : icmin_olap,icmax_olap, & 
                               jcmin_olap,jcmax_olap, &
                               kcmin_olap,kcmax_olap,error_abort
    implicit none
 
    logical, intent(out), optional                        :: recv_flag
    integer, intent(in), optional                         :: icmax_recv,icmin_recv
    integer, intent(in), optional                         :: jcmax_recv,jcmin_recv
    integer, intent(in), optional                         :: kcmax_recv,kcmin_recv
    real(kind=kind(0.d0)),dimension(:,:,:,:), intent(out) :: arecv
    
    integer :: icmin,icmax,jcmin,jcmax,kcmin,kcmax

    !if ((present(icmax_recv))) print*, 'icmax_recv', icmax_recv
    !if ((present(icmin_recv))) print*, 'icmin_recv', icmin_recv
    !if ((present(jcmax_recv))) print*, 'jcmax_recv', jcmax_recv
    !if ((present(jcmin_recv))) print*, 'jcmin_recv', jcmin_recv
    !if ((present(kcmax_recv))) print*, 'kcmax_recv', kcmax_recv
    !if ((present(kcmin_recv))) print*, 'kcmin_recv', kcmin_recv

    !Revert to default i domain sending - top of overlap to bottom of overlap
    if ((present(icmax_recv)) .and. (present(icmin_recv))) then
            icmax = icmax_recv; icmin = icmin_recv
    elseif ((.not. present(icmax_recv)).and.(.not. present(icmin_recv))) then
            icmax = icmax_olap; icmin = icmin_olap
    else
        call error_abort("CPL_recv error - both max and min i limits " // &
                         "required and only one supplied")
    endif

    !Revert to default j domain sending - top of overlap to bottom of overlap
    if ((present(jcmax_recv)) .and. (present(jcmin_recv))) then
            jcmax = jcmax_recv; jcmin = jcmin_recv
    elseif ((.not. present(jcmax_recv)).and.(.not. present(jcmin_recv))) then
            jcmax = jcmax_olap; jcmin = jcmin_olap
    else
        call error_abort("CPL_recv error - both max and min j limits " // &
                         "required and only one supplied")
    endif

    !Revert to default k domain sending - top of overlap to bottom of overlap
    if ((present(kcmax_recv)) .and. (present(kcmin_recv))) then
            kcmax = kcmax_recv; kcmin = kcmin_recv
    elseif ((.not. present(kcmax_recv)).and.(.not. present(kcmin_recv))) then
            kcmax = kcmax_olap; kcmin = kcmin_olap
    else
        call error_abort("CPL_recv error - both max and min k limits " // &
                         "required and only one supplied")
    endif
 
    call CPL_recv_xd(arecv,icmin,icmax,jcmin, & 
                           jcmax,kcmin,kcmax,recv_flag)

end subroutine CPL_recv_4d


!=============================================================================
!                       CPL_recv_xd
!>
!! Receive data from to local grid from the associated ranks from the other 
!! realm
!!
!! - Synopsis
!!
!!  - CPL_recv_xd(arecv,icmin_recv,icmax_recv,jcmin_recv,  
!!                           jcmax_recv,kcmin_recv,kcmax_recv,recv_flag)
!!
!! - Input Parameters
!!
!!   - arecv
!!
!!   - icmin_recv
!!
!!   - icmax_recv
!!
!!   - jcmin_recv
!!
!!   - jcmax_recv
!!
!!   - kcmin_recv
!!
!!   - kcmax_recv
!!
!! - Output Parameter
!!
!!   - recv_flag
!!
!! @author Edward Smith
! ----------------------------------------------------------------------------
!-----------------------------------------------------------------------------
subroutine CPL_recv_xd(arecv,icmin_recv,icmax_recv,jcmin_recv, & 
                             jcmax_recv,kcmin_recv,kcmax_recv,recv_flag)
    use mpi
    use coupler_module, only : md_realm,cfd_realm, & 
                               rank_graph, &
                               error_abort,CPL_GRAPH_COMM,myid_graph,olap_mask, &
                               rank_world, realm, & 
                               iblock_realm,jblock_realm,kblock_realm,VOID,ierr
    implicit none

    !Flag set if processor has received data
    logical, intent(out), optional                  :: recv_flag

    ! Minimum and maximum values of j to send
    integer, intent(in) :: icmax_recv,icmin_recv, & 
                           jcmax_recv,jcmin_recv, & 
                           kcmax_recv,kcmin_recv

    ! Array that recieves grid distributed data 
    real(kind(0.d0)), dimension(:,:,:,:),intent(inout):: arecv     

    !Neighbours
    integer                             :: nneighbors   
    integer,dimension(:),allocatable    :: id_neighbors
                                                         
    ! local indices 
    integer :: n,nbr,icell,jcell,kcell
    integer :: pos,iclmin,iclmax,jclmin,jclmax,kclmin,kclmax
    integer :: pcoords(3),npercell,ndata,ncells
    integer,dimension(6) :: extents,portion,limits

    ! auxiliaries 
    integer :: itag, sourceid,start_address
    integer,dimension(:),allocatable   :: req
    integer,dimension(:,:),allocatable :: status
    real(kind(0.d0)),dimension(:), allocatable ::  vbuf
 
    ! This local CFD domain is outside MD overlap zone 
    if (olap_mask(rank_world).eqv. .false.) return

    ! Save limits array of Minimum and maximum values to recv
    limits = (/ icmin_recv,icmax_recv,jcmin_recv,jcmax_recv,kcmin_recv,kcmax_recv /)

    ! Number of components at each grid point
    npercell = size(arecv,1)

    ! Get local grid box ranges seen by this rank for CFD and allocate buffer
    if (realm .eq. cfd_realm) then 

        !Load CFD cells per processor
        call CPL_Cart_coords(CPL_GRAPH_COMM, rank_graph, cfd_realm, 3, pcoords, ierr) 
        call CPL_proc_extents(pcoords,cfd_realm,extents)

        ! If limits passed to recv routine, use these instead
        ! of overlap/processor limits
        call CPL_proc_portion(pcoords,cfd_realm,limits,portion,ncells)

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
        call CPL_proc_portion(pcoords,md_realm,limits,portion,ncells)

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

    !if (rank_world .eq. 33) then
    !   do n = 1,size(vbuf)
    !       write(98000+rank_world,*) rank_world,n, vbuf(n)
    !   enddo
    !endif

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
        elseif (realm .eq. md_realm) then
            !MD realm receives data as big as own processor domain
            pcoords = (/iblock_realm,jblock_realm,kblock_realm /)
            !Get extents of current processor/overlap region
            call CPL_proc_extents(pcoords,md_realm,extents)
        endif

        ! If limits passed to recv routine, use these instead
        ! of overlap/processor limits
        call CPL_proc_portion(pcoords,md_realm,limits,portion,ncells)
                
        ! Unpack array into buffer
        if (any(portion.eq.VOID)) then
            !print*, 'VOID recv qqqq',realm_name(realm),rank_world,rank_realm,rank_graph2rank_world(sourceid+1),recv_flag
            ndata = 0
        else
            ! Get local extents in received region
            pos = start_address; ndata = npercell * ncells
            iclmin = portion(1)-extents(1)+1;   iclmax = portion(2)-extents(1)+1
            jclmin = portion(3)-extents(3)+1;   jclmax = portion(4)-extents(3)+1
            kclmin = portion(5)-extents(5)+1;   kclmax = portion(6)-extents(5)+1
            !print'(a,5i4,2i6,i4,18i4,l)', 'recv qqqq',rank_world,rank_realm,rank_olap,ndata,nbr, & 
            !                           rank_graph2rank_world(sourceid+1),size(arecv),start_address,&
            !                           iclmin,iclmax,jclmin,jclmax,kclmin,kclmax, & 
            !                           portion,extents,recv_flag
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
           
end subroutine CPL_recv_xd

!-------------------------------------------------------------------

subroutine CPL_pack(unpacked,packed,realm,icmax_pack,icmin_pack,jcmax_pack, & 
                                          jcmin_pack,kcmax_pack,kcmin_pack    )
    use coupler_module, only : CPL_CART_COMM,rank_cart,md_realm,cfd_realm, & 
                               error_abort,CPL_GRAPH_COMM,myid_graph,realm_name
    implicit none

    integer, intent(in)                                         :: realm
    real(kind=kind(0.d0)),dimension(:,:,:,:), intent(in)        :: unpacked
    real(kind=kind(0.d0)),dimension(:),allocatable, intent(out) :: packed

    ! Optional minimum and maximum values to pack
    integer, intent(in), optional    :: icmax_pack,icmin_pack
    integer, intent(in), optional    :: jcmax_pack,jcmin_pack
    integer, intent(in), optional    :: kcmax_pack,kcmin_pack

    integer                          :: pos,n,nbr,id_nbr,icell,jcell,kcell,ierr
    integer                          :: npercell,ncells,nneighbors
    integer,dimension(3)             :: coord
    integer,dimension(6)             :: extents,gextents
    integer,dimension(:),allocatable :: id_neighbors

    !Amount of data per cell
    npercell = size(unpacked,1)

    !Allocate packing buffer
    if (allocated(packed)) deallocate(packed)
    allocate(packed(size(unpacked)))

    ! Get neighbour topology to determine ordering of packed data
    call MPI_Graph_neighbors_count(CPL_GRAPH_COMM,myid_graph,nneighbors,ierr)
    allocate(id_neighbors(nneighbors))
    call MPI_Graph_neighbors(CPL_GRAPH_COMM,myid_graph,nneighbors,id_neighbors,ierr)

    ! Loop through all neighbours which will be sent to and order the data 
    ! appropriately to send each correctly
    do nbr = 1,nneighbors

        if (realm .eq. cfd_realm) then

            ! Get MD neighbour
            id_nbr = id_neighbors(nbr)
            call CPL_Cart_coords(CPL_GRAPH_COMM,id_nbr+1,md_realm,3,coord,ierr) 
            call CPL_olap_extents(coord,md_realm,extents,ncells)

            ! Get offset of neighbouring processor
            pos = id_nbr * npercell * ncells

            !print*,'Pack',rank_cart,realm_name(realm),coord,nbr,id_nbr,extents,coord

        elseif (realm .eq. md_realm) then
            !Get own processor coordinates 
            call CPL_Cart_coords(CPL_CART_COMM,rank_cart,realm,3,coord,ierr) 
            call CPL_olap_extents(coord,realm,gextents,ncells)
            ! Get local extents
            extents(1) = 1; extents(2) = gextents(2)-gextents(1)
            extents(3) = 1; extents(4) = gextents(4)-gextents(3)
            extents(5) = 1; extents(6) = gextents(6)-gextents(5)
            pos = 1
        endif

        ! Pack array into buffer
        do kcell=extents(5),extents(6)
        do jcell=extents(3),extents(4)
        do icell=extents(1),extents(2)
        do n = 1,npercell

            packed(pos) = unpacked(n,icell,jcell,kcell)
            pos = pos + 1

        end do
        end do
        end do
        end do

    end do

    !Sanity check
    if (size(packed) .ne. npercell*ncells) then
        !print*, 'data size', size(packed), 'expected size', npercell*ncells
        call error_abort("CPL_pack error - cell array does not match expected extents")
    endif

end subroutine CPL_pack


!-------------------------------------------------------------------

subroutine CPL_unpack(packed,unpacked,realm)
    use coupler_module, only : CPL_CART_COMM,rank_cart,md_realm,cfd_realm, & 
                               error_abort,CPL_GRAPH_COMM,myid_graph
    implicit none

    integer, intent(in)                                               :: realm
    real(kind=kind(0.d0)),dimension(:,:,:,:),allocatable, intent(out) :: unpacked
    real(kind=kind(0.d0)),dimension(:),allocatable, intent(inout)     :: packed

    integer                          :: pos,n,nbr,id_nbr,icell,jcell,kcell,ierr
    integer                          :: npercell,ncells,nneighbors
    integer,dimension(3)             :: coord
    integer,dimension(6)             :: extents,gextents
    integer,dimension(:),allocatable :: id_neighbors

    call CPL_Cart_coords(CPL_CART_COMM,rank_cart,realm,3,coord,ierr) 
    call CPL_proc_extents(coord,realm,extents,ncells)

    !Amount of data per cell
    npercell = size(packed)/ncells

    !Allocate packing buffer
    if (allocated(unpacked)) deallocate(unpacked)
    allocate(unpacked(npercell,1:extents(2)-extents(1), &
                               1:extents(4)-extents(3), &
                               1:extents(6)-extents(5)))

    ! Get neighbour topology to determine ordering of packed data
    call MPI_Graph_neighbors_count(CPL_GRAPH_COMM,myid_graph,nneighbors,ierr)
    allocate(id_neighbors(nneighbors))
    call MPI_Graph_neighbors(CPL_GRAPH_COMM,myid_graph,nneighbors,id_neighbors,ierr)

    ! Loop through all neighbours which will be sent to and order the data 
    ! appropriately to send each correctly
    do nbr = 1,nneighbors

        if (realm .eq. cfd_realm) then
            ! Get MD neighbour
            id_nbr = id_neighbors(nbr)
            call CPL_Cart_coords(CPL_GRAPH_COMM,id_nbr,md_realm,3,coord,ierr) 
            call CPL_proc_extents(coord,md_realm,extents,ncells)
            ! Get offset of neighbouring processor
            pos = id_nbr * npercell * ncells    !ASSUMES all ncell the same!!
        elseif (realm .eq. md_realm) then
            !Get own processor coordinates 
            call CPL_Cart_coords(CPL_CART_COMM,rank_cart,realm,3,coord,ierr) 
            call CPL_proc_extents(coord,realm,gextents,ncells)
            ! Get local extents
            extents(1) = 1; extents(2) = gextents(2)-gextents(1)
            extents(3) = 1; extents(4) = gextents(4)-gextents(3)
            extents(5) = 1; extents(6) = gextents(6)-gextents(5)
            pos = 1
        endif

        ! Unpack buffer into array
        do kcell=extents(5),extents(6)
        do jcell=extents(3),extents(4)
        do icell=extents(1),extents(2)
        do n = 1,npercell

            unpacked(n,icell,jcell,kcell) = packed(pos)
            pos = pos + 1

        end do
        end do
        end do
        end do

    end do

    !Deallocate packed buffer
    deallocate(packed)

end subroutine CPL_unpack


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
!! @author David Trevelyan
subroutine CPL_proc_extents(coord,realm,extents,ncells)
    use mpi
    use coupler_module, only: md_realm,      cfd_realm,      &
                              icPmin_md,     icPmax_md,      &
                              jcPmin_md,     jcPmax_md,      &
                              kcPmin_md,     kcPmax_md,      &
                              icPmin_cfd,    icPmax_cfd,     &
                              jcPmin_cfd,    jcPmax_cfd,     &
                              kcPmin_cfd,    kcPmax_cfd,     &
                              error_abort
    implicit none

    integer, intent(in)  :: coord(3), realm
    integer, intent(out) :: extents(6)
    integer, optional, intent(out) :: ncells

    select case(realm)
    case(md_realm)
        extents = (/icPmin_md(coord(1)),icPmax_md(coord(1)), & 
                    jcPmin_md(coord(2)),jcPmax_md(coord(2)), & 
                    kcPmin_md(coord(3)),kcPmax_md(coord(3))/)
    case(cfd_realm)
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
!! @author David Trevelyan
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
!>
!! Get maximum and minimum cell indices, i.e. the 'portion', of the
!! input cell extents 'limits' that is contributed by the current
!! overlapping processor. 
!!
!! - Synopsis
!!  - CPL_proc_portion(coord,realm,limits,portion,ncells)
!!
!! - Input
!!
!!  - coord
!!   - processor cartesian coordinate (3 x integer) 
!!
!!  - realm
!!   - cfd_realm (1) or md_realm (2) (integer) 
!!
!!  - limits(6)
!!   - Array of cell extents that specify the input region. 
!!
!! - Input/Output
!!  - NONE
!!
!! - Output
!!
!!  - portion(6) 
!!   - Array of cell extents that define the local processor's
!!     contribution to the input region 'limits'.
!!
!!   - ncells (optional)
!!    - number of cells in portion (integer) 
!!
!! - Note: limits(6) and portion(6) are of the form:
!!   (xmin,xmax,ymin,ymax,zmin,zmax)
!!
!! @author David Trevelyan

subroutine CPL_proc_portion(coord,realm,limits,portion,ncells)
    use mpi
    use coupler_module, only: VOID
    implicit none

    integer, intent(in)  :: coord(3), limits(6),realm
    integer, intent(out) :: portion(6)
    integer, optional, intent(out) :: ncells
    integer :: extents(6)

    call CPL_proc_extents(coord,realm,extents)

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

!-------------------------------------------------------------------
!                   CPL_Cart_coords                                -
!-------------------------------------------------------------------

!>
!! Determines process coords in appropriate realm's cartesian topology 
!! given a rank in any communicator
!!
!! - Synopsis
!!
!!  - CPL_Cart_coords(COMM, rank, realm, maxdims, coords, ierr)
!!
!! - Input Parameters
!!
!!  - comm
!!   - communicator with cartesian structure (handle) 
!!
!!  - realm
!!   - cfd_realm (1) or md_realm (2) (integer) 
!!
!!  - rank
!!   - rank of a process within group of comm (integer) 
!!      NOTE fortran convention rank=1 to nproc
!!
!!  - maxdims
!!   - length of vector coords in the calling program (integer) 
!!
!! - Output Parameter
!!
!!  - coords
!!   - integer array (of size ndims) containing the Cartesian coordinates 
!!     of specified process (integer) 
!!
!!  - ierr
!!   - error flag
!! @author Edward Smith

subroutine CPL_Cart_coords(COMM, rank, realm, maxdims, coords, ierr)
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
!>
!! Return rank of current processor in specified COMM 
!!
!! - Synopsis
!!
!!  - CPL_get_rank(COMM, rank)
!!
!! - Input Parameters
!!
!!  - comm
!!   - communicator with cartesian structure (handle) 
!!
!! - Output Parameter
!!
!!  - rank
!!   - rank of a process within group of comm (integer) 
!!      NOTE fortran convention rank=1 to nproc
!!
!! @author Edward Smith

subroutine CPL_get_rank(COMM,rank)
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
!! @author Edward Smith


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
!>
!! Wrapper to retrieve (read only) parameters from the coupler_module 
!! Note - this ensures all variable in the coupler are protected
!! from corruption by either CFD or MD codes
!!
!! - Synopsis
!!
!!  - CPL_get([see coupler_module])
!!
!! - Input Parameters
!!
!!  - NONE
!!
!! - Output Parameter
!!
!!  - @see coupler_module
!!
!! @author Edward Smith

subroutine CPL_get(icmax_olap,icmin_olap,jcmax_olap,jcmin_olap,  & 
                   kcmax_olap,kcmin_olap,density_cfd,density_md, &
                   dt_cfd,dt_MD,dx,dy,dz,ncx,ncy,ncz,xg,yg,zg,   &
                   xL_md,xL_cfd,yL_md,yL_cfd,zL_md,zL_cfd,       &
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
!> Get molecule's global position from position local to processor.
!-----------------------------------------------------------------------------
function globalise(r) result(rg)
    use coupler_module, only :  xLl,iblock_realm,npx_md, & 
                                yLl,jblock_realm,npy_md, & 
                                zLl,kblock_realm,npz_md
    implicit none

    real(kind(0.d0)),intent(in) :: r(3)
    real(kind(0.d0))            :: rg(3)

    rg(1) = r(1) - 0.5d0*xLl*(npx_md-1) + xLl*(iblock_realm-1)
    rg(2) = r(2) - 0.5d0*yLl*(npy_md-1) + yLl*(jblock_realm-1)
    rg(3) = r(3) - 0.5d0*zLl*(npz_md-1) + zLl*(kblock_realm-1)

end function globalise

!=============================================================================
!> Get local position on processor from molecule's global position.
!-----------------------------------------------------------------------------
function localise(r) result(rg)
    use coupler_module, only :  xLl,iblock_realm,npx_md, & 
                                yLl,jblock_realm,npy_md, & 
                                zLl,kblock_realm,npz_md
    implicit none

    real(kind(0.d0)),intent(in) :: r(3)
    real(kind(0.d0)) rg(3)

    !Global domain has origin at centre
    rg(1) = r(1) - xLl*(iblock_realm-1)+0.5d0*xLl*(npx_md-1)
    rg(2) = r(2) - yLl*(jblock_realm-1)+0.5d0*yLl*(npy_md-1)
    rg(3) = r(3) - zLl*(kblock_realm-1)+0.5d0*zLl*(npz_md-1)

end function localise

!=============================================================================
!> Map global MD position to global CFD coordinate frame
!-----------------------------------------------------------------------------
function map_md2cfd_global(r) result(rg)
    use coupler_module, only :  xL_md,xg,icmin_olap,icmax_olap, & 
                                yL_md,yg,jcmin_olap,jcmax_olap, & 
                                zL_md,zg,kcmin_olap,kcmax_olap
    implicit none

    real(kind(0.d0)),intent(in) :: r(3)
    real(kind(0.d0)):: md_only(3), rg(3)

    !Get size of MD domain which has no CFD cells overlapping
    !This should be general enough to include grid stretching
    !and total overlap in any directions 
    md_only(1) = xL_md-(xg(icmax_olap+1,1) - xg(icmin_olap,1))
    md_only(2) = yL_md-(yg(1,jcmax_olap+1) - yg(1,jcmin_olap))
    md_only(3) = zL_md-(zg( kcmax_olap+1 ) - zg( kcmin_olap ))

    ! CFD has origin at bottom left while MD origin at centre
    rg(1) = r(1) + 0.5d0*xL_md - md_only(1)
    rg(2) = r(2) + 0.5d0*yL_md - md_only(2)
    rg(3) = r(3) + 0.5d0*zL_md - md_only(3)

end function map_md2cfd_global


!=============================================================================
!> Map global CFD position in global MD coordinate frame
!-----------------------------------------------------------------------------
function map_cfd2md_global(r) result(rg)
    use coupler_module, only :  xL_md,xg,icmin_olap,icmax_olap, & 
                                yL_md,yg,jcmin_olap,jcmax_olap, & 
                                zL_md,zg,kcmin_olap,kcmax_olap
    implicit none

    real(kind(0.d0)),intent(in) :: r(3)
    real(kind(0.d0))            :: md_only(3), rg(3)

    !Get size of MD domain which has no CFD cells overlapping
    !This should be general enough to include grid stretching
    !and total overlap in any directions 
    md_only(1) = xL_md-(xg(icmax_olap+1,1) - xg(icmin_olap,1))
    md_only(2) = yL_md-(yg(1,jcmax_olap+1) - yg(1,jcmin_olap))
    md_only(3) = zL_md-(zg( kcmax_olap+1 ) - zg( kcmin_olap ))

    ! CFD has origin at bottom left while MD origin at centre
    rg(1) = r(1) - 0.5d0*xL_md + md_only(1)
    rg(2) = r(2) - 0.5d0*yL_md + md_only(2)
    rg(3) = r(3) - 0.5d0*zL_md + md_only(3)

end function map_cfd2md_global

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

end module coupler
