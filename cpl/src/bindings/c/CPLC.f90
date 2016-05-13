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
!   This file is part of CPL-Library.
!
!   CPL-Library is free software: you can redistribute it and/or modify
!   it under the terms of the GNU General Public License as published by
!   the Free Software Foundation, either version 3 of the License, or
!   (at your option) any later version.
!
!   CPL-Library is distributed in the hope that it will be useful,
!   but WITHOUT ANY WARRANTY; without even the implied warranty of
!   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
!   GNU General Public License for more details.
!
!   You should have received a copy of the GNU General Public License
!   along with CPL-Library.  If not, see <http://www.gnu.org/licenses/>.
!
!
!Description
!
!   Wrappers for the functions found in cpl/src/libsrc/coupler.f90 and
!   cpl/src/libsrc/coupler_module.f90, bound for use with C-type languages.
!
!   The call-signatures for the C-bound functions sometimes require extra
!   arguments to specify shapes of arrays that were otherwise not required by
!   Fortran module interfaces. In light of this difference, the prefixes for
!   each function have been changed from "CPL_" to "CPLC_". So, for example,
!   CPLC_scatter takes extra arguments that define the shape of the array
!   inputs, and wraps CPL_scatter from the original library functions.
!
!   Also included are a number of "get"-style functions for coupler variables.
!
!Author(s)
!
!   David Trevelyan
!

module CPLC

    use iso_c_binding
    implicit none

contains



    subroutine CPLC_create_comm(calling_realm, returned_realm_comm) &
        bind (C, name="CPLC_create_comm")
        use CPL, only: CPL_create_comm
        implicit none
        
        integer(C_INT), value :: calling_realm
        integer(C_INT)        :: returned_realm_comm
        
        integer :: ierror

        call CPL_create_comm(calling_realm, returned_realm_comm, ierror)

    end subroutine CPLC_create_comm




    subroutine CPLC_cfd_init(nsteps, dt, icomm_grid, icoord, npxyz_cfd, xyzL, &
                             ncxyz, density, ijkcmax, ijkcmin, iTmin, iTmax,  &
                             jTmin, jTmax, kTmin, kTmax, xgrid, ygrid, zgrid) &
        bind (C, name="CPLC_cfd_init")
        use CPL, only: coupler_cfd_init
        implicit none

        ! Integers
        integer(C_INT), value :: nsteps
        integer(C_INT), value :: icomm_grid
        type(C_PTR), value :: ijkcmin ! (3)
        type(C_PTR), value :: ijkcmax ! (3)
        type(C_PTR), value :: npxyz_cfd ! (3)
        type(C_PTR), value :: ncxyz ! (3)
        type(C_PTR), value :: iTmin, iTmax ! (npx)
        type(C_PTR), value :: jTmin, jTmax ! (npy)
        type(C_PTR), value :: kTmin, kTmax ! (npz)
        type(C_PTR), value :: icoord ! (3, nprocs)

        ! Reals
        real(C_DOUBLE), value :: dt
        real(C_DOUBLE), value :: density
        type(C_PTR), value :: xyzL ! (3)
        type(C_PTR), value :: xgrid, ygrid ! (ncx, ncy)
        type(C_PTR), value :: zgrid ! (ncz)

        ! Fortran equivalent array pointers
        integer, dimension(:), pointer :: ijkcmin_f, ijkcmax_f
        integer, dimension(:), pointer :: npxyz_cfd_f, ncxyz_f
        integer, dimension(:), pointer :: iTmin_f, iTmax_f
        integer, dimension(:), pointer :: jTmin_f, jTmax_f
        integer, dimension(:), pointer :: kTmin_f, kTmax_f
        integer, dimension(:,:), pointer :: icoord_f
        real(kind(0.d0)), dimension(:), pointer :: xyzL_f
        real(kind(0.d0)), dimension(:), pointer :: zgrid_f
        real(kind(0.d0)), dimension(:,:), pointer :: xgrid_f, ygrid_f

        ! Convenient local variables to be stored
        integer :: npx, npy, npz, nprocs, ncx, ncy, ncz

        ! Convert arrays with convenience local variables first
        call C_F_POINTER(npxyz_cfd, npxyz_cfd_f, [3])
        call C_F_POINTER(ncxyz, ncxyz_f, [3])

        ! Store convenience variables
        npx = npxyz_cfd_f(1)
        npy = npxyz_cfd_f(2)
        npz = npxyz_cfd_f(3)
        nprocs = npx*npy*npz
        ncx = ncxyz_f(1)
        ncy = ncxyz_f(2)
        ncz = ncxyz_f(3)

        ! Convert remaining array pointers
        ! Integers
        call C_F_POINTER(ijkcmin, ijkcmin_f, [3])
        call C_F_POINTER(ijkcmax, ijkcmax_f, [3])
        call C_F_POINTER(npxyz_cfd, npxyz_cfd_f, [3])
        call C_F_POINTER(ncxyz, ncxyz_f, [3])
        call C_F_POINTER(iTmin, iTmin_f, [npx])
        call C_F_POINTER(iTmax, iTmax_f, [npx])
        call C_F_POINTER(jTmin, jTmin_f, [npy])
        call C_F_POINTER(jTmax, jTmax_f, [npy])
        call C_F_POINTER(kTmin, kTmin_f, [npz])
        call C_F_POINTER(kTmax, kTmax_f, [npz])
        call C_F_POINTER(icoord, icoord_f, [3, nprocs])

        ! Reals
        call C_F_POINTER(xyzL, xyzL_f, [3])
        call C_F_POINTER(xgrid, xgrid_f, [ncx+1, ncy+1])
        call C_F_POINTER(ygrid, ygrid_f, [ncx+1, ncy+1])
        call C_F_POINTER(zgrid, zgrid_f, [ncz+1])

        ! Call the old routine
        call coupler_cfd_init(nsteps, dt, icomm_grid, icoord_f, npxyz_cfd_f,  &
                              xyzL_f, ncxyz_f, density, ijkcmax_f, ijkcmin_f, &
                              iTmin_f, iTmax_f, jTmin_f, jTmax_f, kTmin_f,    &
                              kTmax_f, xgrid_f, ygrid_f, zgrid_f)

    end subroutine CPLC_cfd_init



    subroutine CPLC_test_python(int_p, doub_p, bool_p, int_pptr, doub_pptr, int_pptr_dims, doub_pptr_dims) &
       bind (C, name="CPLC_test_python")
       !use CPL, only: CPL_test_python
       implicit none

       integer(C_INT), value :: int_p
       real(C_DOUBLE), value :: doub_p
       logical(C_BOOL), value :: bool_p
       type(C_PTR), value :: int_pptr
       type(C_PTR), value :: doub_pptr
       type(C_PTR), value :: int_pptr_dims
       type(C_PTR), value :: doub_pptr_dims

       !Fortran equivalent array pointers
       integer, dimension(:,:), pointer :: int_pptr_f
       real(kind(0.d0)), dimension(:,:), pointer :: doub_pptr_f
       integer, dimension(:), pointer :: int_pptr_dims_f, doub_pptr_dims_f

       !Aux variables
       integer :: int_dimx, int_dimy
       integer :: doub_dimx, doub_dimy
       integer :: i, j

       call C_F_POINTER(int_pptr_dims, int_pptr_dims_f, [2])
       call C_F_POINTER(doub_pptr_dims, doub_pptr_dims_f, [2])


       int_dimx = int_pptr_dims_f(1)
       int_dimy = int_pptr_dims_f(2)
       print * , 'int_dimx: ', int_dimx
       print * , 'int_dimy: ', int_dimy

       doub_dimx = doub_pptr_dims_f(1)
       doub_dimy = doub_pptr_dims_f(2)
       print * , 'doub_dimx: ', doub_dimx
       print * , 'doub_dimy: ', doub_dimy

       

       call C_F_POINTER(int_pptr, int_pptr_f, [int_dimx, int_dimy])
       call C_F_POINTER(doub_pptr, doub_pptr_f, [doub_dimx, doub_dimy])

       print *, 'int_p: ', int_p
       print *, 'doub_p: ', doub_p
       if (bool_p) then
         print *, 'bool_p: True'
       else
         print *, 'bool_p: False'
       end if
       print *, 'int_pptr: ', int_pptr_f
       !do i=1,int_dimx
       ! do j=1,int_dimy
       !   print *, int_pptr_f(i,j)
       ! end do
       !end do
       print *, 'doub_pptr: ', doub_pptr_f(5, 2)

       !integer(C_INT) :: returned_realm_comm
       !call CPL_test_python()
    end subroutine CPLC_test_python

    subroutine CPLC_md_init(nsteps, initialstep, dt, icomm_grid, icoord, &
                            npxyz_md, globaldomain, density) &
        bind(C, name="CPLC_md_init")
        use CPL, only: coupler_md_init
        implicit none
       
        ! Integers
        integer(C_INT) :: nsteps ! NOT by value, will be modified
        integer(C_INT) :: initialstep ! NOT by value, will be modified
        integer(C_INT), value :: icomm_grid
        type(C_PTR), value :: icoord ! (3, nprocs)
        type(C_PTR), value :: npxyz_md ! (3)
       
        ! Reals
        real(C_DOUBLE), value :: dt
        real(C_DOUBLE), value :: density
        type(C_PTR), value :: globaldomain ! (3)

        ! Fortran equivalent arrays
        integer, dimension(:), pointer :: npxyz_md_f
        integer, dimension(:,:), pointer :: icoord_f
        real(kind(0.d0)), dimension(:), pointer :: globaldomain_f

        ! Other function internals
        integer :: nprocs

        ! C -> Fortran array pointer type conversions
        call C_F_POINTER(npxyz_md, npxyz_md_f, [3])
        call C_F_POINTER(globaldomain, globaldomain_f, [3])

        ! Store convenience variables
        nprocs = npxyz_md_f(1)*npxyz_md_f(2)*npxyz_md_f(3)
   
        ! Remaining conversions
        call C_F_POINTER(icoord, icoord_f, [3, nprocs])

        ! Call routine
        call coupler_md_init(nsteps, initialstep, dt, icomm_grid, icoord_f, &
                             npxyz_md_f, globaldomain_f, density)
    
    end subroutine CPLC_md_init




    subroutine CPLC_send(asend, asend_shape, ndims, icmax, icmin, jcmax, &
                         jcmin, kcmax, kcmin, send_flag) &
        bind(C, name="CPLC_send")
        use CPL, only: CPL_send!, CPL_send_4d, CPL_send_3d
        implicit none

        ! Integers
        integer(C_INT), value :: ndims
        integer(C_INT), value :: icmax, icmin
        integer(C_INT), value :: jcmax, jcmin
        integer(C_INT), value :: kcmax, kcmin

        ! Boolean
        logical(C_BOOL) :: send_flag

        ! Inputs
        type(C_PTR), value :: asend
        type(C_PTR), value :: asend_shape

        ! Fortran equivalent array pointers
        real(kind(0.d0)), dimension(:,:,:,:), pointer :: asend4D_f
        real(kind(0.d0)), dimension(:,:,:), pointer :: asend3D_f
        integer, dimension(:), pointer :: asend_shape_f

        ! Other useful variables internal to this subroutine
        integer :: s1, s2, s3, s4

        send_flag  = .true.

        if (ndims .eq. 3) then
          call C_F_POINTER(asend_shape, asend_shape_f, [3])
          s1 = asend_shape_f(1) 
          s2 = asend_shape_f(2) 
          s3 = asend_shape_f(3) 
          call C_F_POINTER(asend, asend3D_f, [s1, s2, s3])
          call CPL_send(asend3D_f, icmax, icmin, jcmax, jcmin, kcmax, kcmin)
        elseif (ndims .eq. 4) then
          call C_F_POINTER(asend_shape, asend_shape_f, [4])
          s1 = asend_shape_f(1) 
          s2 = asend_shape_f(2) 
          s3 = asend_shape_f(3) 
          s4 = asend_shape_f(4) 
          call C_F_POINTER(asend, asend4D_f, [s1, s2, s3, s4])
          call CPL_send(asend4D_f, icmax, icmin, jcmax, jcmin, kcmax, kcmin)
        else
          print*, "Error in CPLC_send -- only 3 and 4 is allowed for ndims."
        end if


    end subroutine CPLC_send


    subroutine CPLC_recv(arecv, arecv_shape, ndims, icmax, icmin, jcmax, &
                         jcmin, kcmax, kcmin, recv_flag) &
        bind(C, name="CPLC_recv")
        use CPL, only: CPL_recv
        implicit none

        ! Integers
        integer(C_INT), value :: ndims
        integer(C_INT), value :: icmax, icmin
        integer(C_INT), value :: jcmax, jcmin
        integer(C_INT), value :: kcmax, kcmin


        ! Boolean
        logical(C_BOOL) :: recv_flag

        ! Inputs
        type(C_PTR), value :: arecv
        type(C_PTR), value :: arecv_shape

        ! Fortran equivalent array pointers
        real(kind(0.d0)), dimension(:,:,:,:), pointer :: arecv4D_f
        real(kind(0.d0)), dimension(:,:,:), pointer :: arecv3D_f
        integer, dimension(:), pointer :: arecv_shape_f

        ! Other useful variables internal to this subroutine
        integer :: s1, s2, s3, s4

        recv_flag = .true.
        if (ndims .eq. 3) then
          call C_F_POINTER(arecv_shape, arecv_shape_f, [3])
          s1 = arecv_shape_f(1) 
          s2 = arecv_shape_f(2) 
          s3 = arecv_shape_f(3) 
          call C_F_POINTER(arecv, arecv3D_f, [s1, s2, s3])
          call CPL_recv(arecv3D_f, icmax, icmin, jcmax, jcmin, kcmax, kcmin)
        elseif (ndims .eq. 4) then
          call C_F_POINTER(arecv_shape, arecv_shape_f, [4])
          s1 = arecv_shape_f(1) 
          s2 = arecv_shape_f(2) 
          s3 = arecv_shape_f(3) 
          s4 = arecv_shape_f(4) 
          call C_F_POINTER(arecv, arecv4D_f, [s1, s2, s3, s4])
          call CPL_recv(arecv4D_f, icmax, icmin, jcmax, jcmin, kcmax, kcmin)
        else
          print*, "Error only 3 and 4 is allowed for ndims."
        end if


    end subroutine CPLC_recv




    subroutine CPLC_scatter(scatterarray, scatter_shape, limits, recvarray, &
                            recv_shape) &
        bind(C, name="CPLC_scatter")
        use CPL, only: CPL_scatter
        implicit none
   
        ! Inputs
        type(C_PTR), value :: scatterarray ! (sn, sx, sy, sz) (see below)
        type(C_PTR), value :: scatter_shape ! (4)
        type(C_PTR), value :: limits ! (6)
        type(C_PTR), value :: recv_shape ! (4)
        type(C_PTR), value :: recvarray ! (rn, rx, ry, rz) (see below)

        ! Fortran equivalent array pointers
        integer, dimension(:), pointer :: scatter_shape_f
        integer, dimension(:), pointer :: limits_f
        integer, dimension(:), pointer :: recv_shape_f
        real(kind(0.d0)), dimension(:,:,:,:), pointer :: scatterarray_f
        real(kind(0.d0)), dimension(:,:,:,:), pointer :: recvarray_f

        ! Other useful variables internal to this subroutine
        integer :: sx, sy, sz ! Number of cells in (s)cattered array
        integer :: rx, ry, rz ! Number of cells in (r)eceived array
        integer :: sn ! Number of values per cell in sent array
        integer :: rn ! Number of values per cell in recv array

        ! Store limits array first
        call C_F_POINTER(scatter_shape, scatter_shape_f, [4])
        call C_F_POINTER(recv_shape, recv_shape_f, [4])
        call C_F_POINTER(limits, limits_f, [6])

        ! Store shape of scattered array
        sn = scatter_shape_f(1)
        sx = scatter_shape_f(2)
        sy = scatter_shape_f(3)
        sz = scatter_shape_f(4)

        ! Store shape of received array
        !if (recv_shape_f(1) .ne. sn) then
        !    call error_abort("Error: npercell in scatter_shape and "// &
        !                     "recv_shape must are not the same. Aborting.")
        !end if
        rn = recv_shape_f(1)
        rx = recv_shape_f(2)
        ry = recv_shape_f(3)
        rz = recv_shape_f(4)

        ! Remaining conversions
        call C_F_POINTER(scatterarray, scatterarray_f, [sn, sx, sy, sz])
        call C_F_POINTER(recvarray, recvarray_f, [rn, rx, ry, rz])

        ! Library call
        call CPL_scatter(scatterarray_f, sn, limits_f, recvarray_f)

    end subroutine CPLC_scatter




    subroutine CPLC_gather(gatherarray, gather_shape, limits, recvarray, &
                            recv_shape) &
        bind(C, name="CPLC_gather")
        use CPL, only: CPL_gather
        implicit none
   
        ! Inputs
        type(C_PTR), value :: gatherarray ! (gn, gx, gy, gz) (see below)
        type(C_PTR), value :: gather_shape ! (4)
        type(C_PTR), value :: limits ! (6)
        type(C_PTR), value :: recv_shape ! (4)
        type(C_PTR), value :: recvarray ! (rn, rx, ry, rz) (see below)

        ! Fortran equivalent array pointers
        integer, dimension(:), pointer :: gather_shape_f
        integer, dimension(:), pointer :: limits_f
        integer, dimension(:), pointer :: recv_shape_f
        real(kind(0.d0)), dimension(:,:,:,:), pointer :: gatherarray_f
        real(kind(0.d0)), dimension(:,:,:,:), pointer :: recvarray_f

        ! Other useful variables internal to this subroutine
        integer :: gx, gy, gz ! Number of cells in (s)cattered array
        integer :: rx, ry, rz ! Number of cells in (r)eceived array
        integer :: gn, rn ! Number of values per cell in arrays

        ! Store limits array first
        call C_F_POINTER(gather_shape, gather_shape_f, [4])
        call C_F_POINTER(recv_shape, recv_shape_f, [4])
        call C_F_POINTER(limits, limits_f, [6])

        ! Store shape of gathered array
        gn = gather_shape_f(1)
        gx = gather_shape_f(2)
        gy = gather_shape_f(3)
        gz = gather_shape_f(4)

        ! Received array shape
        rn = recv_shape_f(1)
        rx = recv_shape_f(2)
        ry = recv_shape_f(3)
        rz = recv_shape_f(4)

        ! Remaining conversions
        call C_F_POINTER(gatherarray, gatherarray_f, [gn, gx, gy, gz])
        call C_F_POINTER(recvarray, recvarray_f, [rn, rx, ry, rz])

        ! Library call
        call CPL_gather(gatherarray_f, rn, limits_f, recvarray_f)


    end subroutine CPLC_gather




    subroutine CPLC_proc_extents(coord, realm, extents) &
        bind(C, name="CPLC_proc_extents")
        use CPL, only: CPL_proc_extents
        
        ! Inputs
        type(C_PTR), value :: coord ! (3)
        type(C_PTR), value :: extents ! (6)
        integer(C_INT), value :: realm

        ! Fortran equivalent arrays
        integer, dimension(:), pointer :: coord_f
        integer, dimension(:), pointer :: extents_f

        call C_F_POINTER(coord, coord_f, [3])
        call C_F_POINTER(extents, extents_f, [6])

        call CPL_proc_extents(coord_f, realm, extents_f)

    end subroutine CPLC_proc_extents




    subroutine CPLC_proc_portion(coord, realm, limits, portion) &
        bind(C, name="CPLC_proc_portion")
        use CPL, only: CPL_proc_portion
        
        ! Inputs
        type(C_PTR), value :: coord ! (3)
        integer(C_INT), value :: realm
        type(C_PTR), value :: limits ! (6)
        type(C_PTR), value :: portion ! (6)

        ! Fortran equivalent arrays
        integer, dimension(:), pointer :: coord_f
        integer, dimension(:), pointer :: limits_f
        integer, dimension(:), pointer :: portion_f

        call C_F_POINTER(coord, coord_f, [3])
        call C_F_POINTER(limits, limits_f, [6])
        call C_F_POINTER(portion, portion_f, [6])

        call CPL_proc_portion(coord_f, realm, limits_f, portion_f)

    end subroutine CPLC_proc_portion




    type(C_PTR) function CPLC_map_cfd2md_global(r_cfd) &
        bind(C, name="CPLC_map_cfd2md_global")
        use CPL, only: map_cfd2md_global
    
        ! Input position
        type(C_PTR), value :: r_cfd ! (3)
       
        ! Fortran equivalent array
        real(kind(0.d0)), dimension(:), pointer :: r_cfd_f
       
        ! Output position
        real(C_DOUBLE), dimension(:), pointer :: r_md_f

        call C_F_POINTER(r_cfd, r_cfd_f, [3])

        ! Allocate array but DON'T deallocate, otherwise data
        ! will be nullified and the calling C program will just
        ! get undefined values 
        allocate(r_md_f(3)) 

        r_md_f = map_cfd2md_global(r_cfd_f)

        CPLC_map_cfd2md_global = C_LOC(r_md_f(1))

    end function CPLC_map_cfd2md_global
     

    ! Setters:

    subroutine CPLC_set_output_mode(mode) &
        bind(C, name="CPLC_set_output_mode")
        use CPL, only: set_output_mode
        
        integer(C_INT), intent(in), value :: mode

        call set_output_mode(mode)
        
    end subroutine CPLC_set_output_mode

    ! Getters: integers

    integer(C_INT) function CPLC_ncx() &
        bind(C, name="CPLC_ncx")
        use CPL, only: ncx
        implicit none

        CPLC_ncx = ncx
        
    end function CPLC_ncx

    integer(C_INT) function CPLC_ncy() &
        bind(C, name="CPLC_ncy")
        use CPL, only: ncy
        implicit none

        CPLC_ncy = ncy
        
    end function CPLC_ncy

    integer(C_INT) function CPLC_ncz() &
        bind(C, name="CPLC_ncz")
        use CPL, only: ncz
        implicit none

        CPLC_ncz = ncz
        
    end function CPLC_ncz

    integer(C_INT) function CPLC_npx_md() &
        bind(C, name="CPLC_npx_md")
        use CPL, only: npx_md
        implicit none

        CPLC_npx_md = npx_md
        
    end function CPLC_npx_md

    integer(C_INT) function CPLC_npy_md() &
        bind(C, name="CPLC_npy_md")
        use CPL, only: npy_md
        implicit none

        CPLC_npy_md = npy_md
        
    end function CPLC_npy_md

    integer(C_INT) function CPLC_npz_md() &
        bind(C, name="CPLC_npz_md")
        use CPL, only: npz_md
        implicit none

        CPLC_npz_md = npz_md
        
    end function CPLC_npz_md


    integer(C_INT) function CPLC_npx_cfd() &
        bind(C, name="CPLC_npx_cfd")
        use CPL, only: npx_cfd
        implicit none

        CPLC_npx_cfd = npx_cfd
        
    end function CPLC_npx_cfd

    integer(C_INT) function CPLC_npy_cfd() &
        bind(C, name="CPLC_npy_cfd")
        use CPL, only: npy_cfd
        implicit none

        CPLC_npy_cfd = npy_cfd
        
    end function CPLC_npy_cfd

    integer(C_INT) function CPLC_npz_cfd() &
        bind(C, name="CPLC_npz_cfd")
        use CPL, only: npz_cfd
        implicit none

        CPLC_npz_cfd = npz_cfd
        
    end function CPLC_npz_cfd


    integer(C_INT) function CPLC_icmin_olap() &
        bind(C, name="CPLC_icmin_olap")
        use CPL, only: icmin_olap
        implicit none

        CPLC_icmin_olap = icmin_olap
        
    end function CPLC_icmin_olap

    integer(C_INT) function CPLC_jcmin_olap() &
        bind(C, name="CPLC_jcmin_olap")
        use CPL, only: jcmin_olap
        implicit none

        CPLC_jcmin_olap = jcmin_olap
        
    end function CPLC_jcmin_olap

    integer(C_INT) function CPLC_kcmin_olap() &
        bind(C, name="CPLC_kcmin_olap")
        use CPL, only: kcmin_olap
        implicit none

        CPLC_kcmin_olap = kcmin_olap
        
    end function CPLC_kcmin_olap

    integer(C_INT) function CPLC_icmax_olap() &
        bind(C, name="CPLC_icmax_olap")
        use CPL, only: icmax_olap
        implicit none

        CPLC_icmax_olap = icmax_olap
        
    end function CPLC_icmax_olap

    integer(C_INT) function CPLC_jcmax_olap() &
        bind(C, name="CPLC_jcmax_olap")
        use CPL, only: jcmax_olap
        implicit none

        CPLC_jcmax_olap = jcmax_olap
        
    end function CPLC_jcmax_olap

    integer(C_INT) function CPLC_kcmax_olap() &
        bind(C, name="CPLC_kcmax_olap")
        use CPL, only: kcmax_olap
        implicit none

        CPLC_kcmax_olap = kcmax_olap
        
    end function CPLC_kcmax_olap

    integer(C_INT) function CPLC_icmin_cnst() &
        bind(C, name="CPLC_icmin_cnst")
        use CPL, only: icmin_cnst
        implicit none

        CPLC_icmin_cnst = icmin_cnst
        
    end function CPLC_icmin_cnst

    integer(C_INT) function CPLC_jcmin_cnst() &
        bind(C, name="CPLC_jcmin_cnst")
        use CPL, only: jcmin_cnst
        implicit none

        CPLC_jcmin_cnst = jcmin_cnst
        
    end function CPLC_jcmin_cnst

    integer(C_INT) function CPLC_kcmin_cnst() &
        bind(C, name="CPLC_kcmin_cnst")
        use CPL, only: kcmin_cnst
        implicit none

        CPLC_kcmin_cnst = kcmin_cnst
        
    end function CPLC_kcmin_cnst

    integer(C_INT) function CPLC_icmax_cnst() &
        bind(C, name="CPLC_icmax_cnst")
        use CPL, only: icmax_cnst
        implicit none

        CPLC_icmax_cnst = icmax_cnst
        
    end function CPLC_icmax_cnst

    integer(C_INT) function CPLC_jcmax_cnst() &
        bind(C, name="CPLC_jcmax_cnst")
        use CPL, only: jcmax_cnst
        implicit none

        CPLC_jcmax_cnst = jcmax_cnst
        
    end function CPLC_jcmax_cnst

    integer(C_INT) function CPLC_kcmax_cnst() &
        bind(C, name="CPLC_kcmax_cnst")
        use CPL, only: kcmax_cnst
        implicit none

        CPLC_kcmax_cnst = kcmax_cnst
        
    end function CPLC_kcmax_cnst

    integer(C_INT) function CPLC_timestep_ratio() &
        bind(C, name="CPLC_timestep_ratio")
        use CPL, only: timestep_ratio
        implicit none

        CPLC_timestep_ratio = timestep_ratio
        
    end function CPLC_timestep_ratio

    integer(C_INT) function CPLC_cpl_md_bc_slice() &
        bind(C, name="CPLC_cpl_md_bc_slice")
        use CPL, only: cpl_md_bc_slice
        implicit none

        CPLC_cpl_md_bc_slice = cpl_md_bc_slice

    end function CPLC_cpl_md_bc_slice

    integer(C_INT) function CPLC_cpl_cfd_bc_slice() &
        bind(C, name="CPLC_cpl_cfd_bc_slice")
        use CPL, only: cpl_cfd_bc_slice
        implicit none

        CPLC_cpl_cfd_bc_slice = cpl_cfd_bc_slice

    end function CPLC_cpl_cfd_bc_slice

    integer(C_INT) function CPLC_cpl_cfd_bc_x() &
        bind(C, name="CPLC_cpl_cfd_bc_x")
        use CPL, only: cpl_cfd_bc_x
        implicit none

        CPLC_cpl_cfd_bc_x = cpl_cfd_bc_x 

    end function CPLC_cpl_cfd_bc_x 

    integer(C_INT) function CPLC_cpl_cfd_bc_y() &
        bind(C, name="CPLC_cpl_cfd_bc_y")
        use CPL, only: cpl_cfd_bc_y
        implicit none

        CPLC_cpl_cfd_bc_y = cpl_cfd_bc_y 

    end function CPLC_cpl_cfd_bc_y 

    integer(C_INT) function CPLC_cpl_cfd_bc_z() &
        bind(C, name="CPLC_cpl_cfd_bc_z")
        use CPL, only: cpl_cfd_bc_z
        implicit none

        CPLC_cpl_cfd_bc_z = cpl_cfd_bc_z

    end function CPLC_cpl_cfd_bc_z 

    integer(C_INT) function CPLC_comm_style() &
        bind(C, name="CPLC_comm_style")
        use CPL, only: comm_style
        implicit none

        CPLC_comm_style = comm_style
        
    end function CPLC_comm_style

    integer(C_INT) function CPLC_comm_style_gath_scat() &
        bind(C, name="CPLC_comm_style_gath_scat")
        use CPL, only: comm_style_gath_scat
        implicit none

        CPLC_comm_style_gath_scat = comm_style_gath_scat
        
    end function CPLC_comm_style_gath_scat

    integer(C_INT) function CPLC_comm_style_send_recv() &
        bind(C, name="CPLC_comm_style_send_recv")
        use CPL, only: comm_style_send_recv
        implicit none

        CPLC_comm_style_send_recv = comm_style_send_recv
        
    end function CPLC_comm_style_send_recv


    integer(C_INT) function CPLC_overlap() &
        bind(C, name="CPLC_overlap")
        use CPL, only: CPL_overlap
        implicit none

        if (CPL_overlap()) then 
            CPLC_overlap = 1
        else
            CPLC_overlap = 0
        endif
        
    end function CPLC_overlap

    !Getters: doubles

    real(C_DOUBLE) function CPLC_density_cfd() &
        bind(C, name="CPLC_density_cfd")
        use CPL, only: density_cfd
        implicit none

        CPLC_density_cfd = density_cfd
        
    end function CPLC_density_cfd

    real(C_DOUBLE) function CPLC_dx() &
        bind(C, name="CPLC_dx")
        use CPL, only: dx
        implicit none

        CPLC_dx = dx
        
    end function CPLC_dx

    real(C_DOUBLE) function CPLC_dy() &
        bind(C, name="CPLC_dy")
        use CPL, only: dy
        implicit none

        CPLC_dy = dy
        
    end function CPLC_dy

    real(C_DOUBLE) function CPLC_dz() &
        bind(C, name="CPLC_dz")
        use CPL, only: dz
        implicit none

        CPLC_dz = dz
        
    end function CPLC_dz

    real(C_DOUBLE) function CPLC_xl_md() &
        bind(C, name="CPLC_xl_md")
        use CPL, only: xl_md
        implicit none

        CPLC_xl_md = xl_md
        
    end function CPLC_xl_md

    real(C_DOUBLE) function CPLC_yl_md() &
        bind(C, name="CPLC_yl_md")
        use CPL, only: yl_md
        implicit none

        CPLC_yl_md = yl_md
        
    end function CPLC_yl_md

    real(C_DOUBLE) function CPLC_zl_md() &
        bind(C, name="CPLC_zl_md")
        use CPL, only: zl_md
        implicit none

        CPLC_zl_md = zl_md
        
    end function CPLC_zl_md


    real(C_DOUBLE) function CPLC_xl_cfd() &
        bind(C, name="CPLC_xl_cfd")
        use CPL, only: xl_cfd
        implicit none

        CPLC_xl_cfd = xl_cfd
        
    end function CPLC_xl_cfd

    real(C_DOUBLE) function CPLC_yl_cfd() &
        bind(C, name="CPLC_yl_cfd")
        use CPL, only: yl_cfd
        implicit none

        CPLC_yl_cfd = yl_cfd
        
    end function CPLC_yl_cfd

    real(C_DOUBLE) function CPLC_zl_cfd() &
        bind(C, name="CPLC_zl_cfd")
        use CPL, only: zl_cfd
        implicit none

        CPLC_zl_cfd = zl_cfd
        
    end function CPLC_zl_cfd

    !Getters: pointers
    
    type(C_PTR) function CPLC_xg() &
        bind(C, name="CPLC_xg")
        use CPL, only: xg
        implicit none

        CPLC_xg = C_LOC(xg)
    
    end function CPLC_xg

    type(C_PTR) function CPLC_yg() &
        bind(C, name="CPLC_yg")
        use CPL, only: yg
        implicit none

        CPLC_yg = C_LOC(yg)
    
    end function CPLC_yg

    type(C_PTR) function CPLC_zg() &
        bind(C, name="CPLC_zg")
        use CPL, only: zg
        implicit none

        CPLC_zg = C_LOC(zg)
    
    end function CPLC_zg




end module CPLC
