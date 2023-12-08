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
!   David Trevelyan, Eduardo Ramos Fernandez
!

module CPLC

    use iso_c_binding
    implicit none

contains

    !Convert c bool to fortran
    subroutine C_F_BOOL(cbool, fbool)
        logical(C_BOOL),intent(in) :: cbool
        logical,intent(out) :: fbool
        if (cbool) then
            fbool = .true.
        else
            fbool = .false.
        endif
    end subroutine C_F_BOOL

    !Convert fortran to c bool
    subroutine F_C_BOOL(fbool, cbool)
        logical,intent(in) :: fbool
        logical(C_BOOL),intent(out) :: cbool
        if (fbool) then
            cbool = .true.
        else
            cbool = .false.
        endif
    end subroutine F_C_BOOL

    !Convert fortran to C array of booleans
    subroutine F_C_BOOL_ARRAY(fbool_array, cbool_array)
        implicit none

        logical, dimension(:), intent(in) :: fbool_array
        logical(C_BOOL), dimension(:), intent(out) :: cbool_array
        integer :: i

        do i = 1, size(fbool_array)
            call F_C_BOOL(fbool_array(i), cbool_array(i))
        end do
    end subroutine F_C_BOOL_ARRAY

    !Convert C to fortran array of booleans
    subroutine C_F_BOOL_ARRAY(cbool_array, fbool_array)
        implicit none

        logical(C_BOOL), dimension(:), intent(in) :: cbool_array
        logical, dimension(:), intent(out) :: fbool_array
        integer :: i

        do i = 1, size(cbool_array)
            call C_F_BOOL(cbool_array(i), fbool_array(i))
        end do
    end subroutine C_F_BOOL_ARRAY
 
    ! CPL_io module routines

    subroutine c_f_string(c_string, f_string)
        implicit none 
        type(C_PTR), intent(in) :: c_string
        character(:), intent(out), allocatable :: f_string
        ! Array for accessing string pointed at by C pointer 
        character, pointer :: string_ptr(:)
        integer :: i    ! string index
        interface
          ! Steal std C library function rather than writing our own.
          function strlen(s) bind(c, name='strlen')
            use iso_c_binding 
            implicit none
            !----
            type(C_PTR), intent(in), value :: s
            integer(c_size_t) :: strlen
          end function strlen
        end interface
        !****
        ! Map C pointer to fortran character array
        call c_f_pointer(c_string, string_ptr, [strlen(c_string)])
        ! Allocate fortran character variable to the c string's length
        allocate(character(size(string_ptr)) :: f_string)
        ! Copy across (with possible kind conversion) characters
        do i = 1, size(string_ptr)
            f_string(i:i) = string_ptr(i)
        end do
    end subroutine c_f_string      

    subroutine f_c_string(f_string, c_string)
        implicit none

        type(C_PTR), intent(out) :: c_string
        character(len=:), intent(inout), allocatable, target :: f_string

        f_string = f_string//C_NULL_CHAR
        c_string = c_loc(f_string)
    end subroutine f_c_string

    subroutine f_c_string_array(f_string_array, c_string_array)
        use commondefs, only: CPL_STRING_MAX_LEN
        implicit none

        type(C_PTR), dimension(:), allocatable, intent(out) :: c_string_array
        integer :: i, s_len
        character(len=CPL_STRING_MAX_LEN), dimension(:), target, allocatable,intent(inout):: f_string_array

        print*, "INSIDE0"
        allocate(c_string_array(size(f_string_array)))
        print*, "INSIDE1"
        do i = 1, size(f_string_array)
            s_len = len_trim(f_string_array(i))
            print*, "loop:",i
            f_string_array(i)(s_len+1:s_len+1) = C_NULL_CHAR
            c_string_array(i) = c_loc(f_string_array(i))
        end do

        print*, "INSIDE2"

    end subroutine f_c_string_array

#ifdef JSON_SUPPORT
    subroutine CPLC_load_param_file(fname) &
        bind (C, name="CPLC_load_param_file")
        use CPL, only: CPL_load_param_file
        implicit none

        type(C_PTR), intent(in), value :: fname
        character(:), allocatable :: fname_f

        call c_f_string(fname, fname_f)

        call CPL_load_param_file(fname_f)

    end subroutine  CPLC_load_param_file

    subroutine CPLC_close_param_file() &
        bind (C, name="CPLC_close_param_file")
        use CPL, only: CPL_close_param_file
        implicit none

        call CPL_close_param_file()
    end subroutine

    subroutine CPLC_get_real_param(section, param_name, real_param) &
        bind (C, name="CPLC_get_real_param")
        use CPL, only: get_real_param 
        implicit none
        
        real(C_DOUBLE), intent(out) :: real_param
        type(C_PTR), intent(in), value :: param_name
        type(C_PTR), intent(in), value :: section
        character(:), allocatable :: param_name_f
        character(:), allocatable :: section_f

        call c_f_string(param_name, param_name_f)
        call c_f_string(section, section_f)

        call get_real_param(section_f, param_name_f, real_param)

    end subroutine CPLC_get_real_param

    subroutine CPLC_get_real_array_param(section, param_name, real_array_param, array_len) &
        bind (C, name="CPLC_get_real_array_param")
        use CPL, only: get_real_array_param
        implicit none
        
        integer(C_INT), intent(out) :: array_len
        type(C_PTR), intent(out) :: real_array_param

        type(C_PTR), intent(in), value :: param_name
        type(C_PTR), intent(in), value :: section

        character(:), allocatable :: param_name_f
        character(:), allocatable :: section_f
        real(kind(0.d0)), dimension(:), allocatable, target, save:: real_array_param_f
        !real(kind(0.d0)), dimension(:), pointer :: r_f

        call c_f_string(param_name, param_name_f)
        call c_f_string(section, section_f)

        call get_real_array_param(section_f, param_name_f, real_array_param_f)

        real_array_param = c_loc(real_array_param_f(1))
        array_len = size(real_array_param_f)
    end subroutine CPLC_get_real_array_param


    subroutine CPLC_get_boolean_param(section, param_name, boolean_param) &
        bind (C, name="CPLC_get_boolean_param")
        use CPL, only: get_boolean_param
        implicit none
        
        logical(C_BOOL), intent(out) :: boolean_param
        type(C_PTR), intent(in), value :: param_name
        type(C_PTR), intent(in), value :: section
        character(:), allocatable :: param_name_f
        character(:), allocatable :: section_f
        logical :: boolean_param_f

        call c_f_string(param_name, param_name_f)
        call c_f_string(section, section_f)

        call C_F_BOOL(boolean_param, boolean_param_f) 
        call get_boolean_param(section_f, param_name_f, boolean_param_f)
        call F_C_BOOL(boolean_param_f, boolean_param) 

    end subroutine CPLC_get_boolean_param

    subroutine CPLC_get_boolean_array_param(section, param_name, boolean_array_param, array_len) &
        bind (C, name="CPLC_get_boolean_array_param")
        use CPL, only: get_boolean_array_param
        implicit none
        
        integer(C_INT), intent(out) :: array_len
        type(C_PTR), intent(out) :: boolean_array_param

        type(C_PTR), intent(in), value :: param_name
        type(C_PTR), intent(in), value :: section

        character(:), allocatable :: param_name_f
        character(:), allocatable :: section_f
        logical, dimension(:), allocatable :: boolean_array_param_f
        logical(C_BOOL), dimension(:), allocatable, target, save:: boolean_array_param_c
        !boolean(kind(0.d0)), dimension(:), pointer :: r_f

        call c_f_string(param_name, param_name_f)
        call c_f_string(section, section_f)

        call get_boolean_array_param(section_f, param_name_f, boolean_array_param_f)
        allocate(boolean_array_param_c(size(boolean_array_param_f)))
        call F_C_BOOL_ARRAY(boolean_array_param_f, boolean_array_param_c)

        boolean_array_param = c_loc(boolean_array_param_c(1))
        array_len = size(boolean_array_param_c)
    end subroutine CPLC_get_boolean_array_param
 
    
   subroutine CPLC_get_int_param(section, param_name, int_param) &
        bind (C, name="CPLC_get_int_param")
        use CPL, only: get_int_param
        implicit none
        
        integer(C_INT), intent(out) :: int_param
        type(C_PTR), intent(in), value :: param_name
        type(C_PTR), intent(in), value :: section
        character(:), allocatable :: param_name_f
        character(:), allocatable :: section_f

        call c_f_string(param_name, param_name_f)
        call c_f_string(section, section_f)

        call get_int_param(section_f, param_name_f, int_param)

    end subroutine CPLC_get_int_param


    subroutine CPLC_get_int_array_param(section, param_name, int_array_param, array_len) &
        bind (C, name="CPLC_get_int_array_param")
        use CPL, only: get_int_array_param
        implicit none
        
        integer(C_INT), intent(out) :: array_len
        type(C_PTR), intent(out) :: int_array_param

        type(C_PTR), intent(in), value :: param_name
        type(C_PTR), intent(in), value :: section

        character(:), allocatable :: param_name_f
        character(:), allocatable :: section_f
        integer, dimension(:), allocatable, target, save :: int_array_param_f

        call c_f_string(param_name, param_name_f)
        call c_f_string(section, section_f)

        call get_int_array_param(section_f, param_name_f, int_array_param_f)
        int_array_param = c_loc(int_array_param_f)
        array_len = size(int_array_param_f)
    end subroutine CPLC_get_int_array_param

   subroutine CPLC_get_string_param(section, param_name, string_param) &
        bind (C, name="CPLC_get_string_param")
        use CPL, only: get_string_param
        implicit none

        type(C_PTR), intent(out) :: string_param
        type(C_PTR), intent(in), value :: param_name
        type(C_PTR), intent(in), value :: section
        character(len=:), allocatable, save:: string_param_f
        character(:), allocatable :: section_f
        character(:), allocatable :: param_name_f


        call c_f_string(param_name, param_name_f)
        call c_f_string(section, section_f)

        call get_string_param(section_f, param_name_f, string_param_f)
        call f_c_string(string_param_f, string_param)

    end subroutine CPLC_get_string_param


    subroutine CPLC_get_string_array_param(section, param_name, string_array_param, array_len) &
        bind (C, name="CPLC_get_string_array_param")
        use CPL, only: get_string_array_param, CPL_STRING_MAX_LEN
        implicit none
        
        integer(C_INT), intent(out) :: array_len
        type(C_PTR), dimension(:), allocatable, intent(out) :: string_array_param

        type(C_PTR), intent(in), value :: param_name
        type(C_PTR), intent(in), value :: section

        character(:), allocatable :: param_name_f
        character(:), allocatable :: section_f
        character(len=CPL_STRING_MAX_LEN), dimension(:), allocatable, save:: string_array_param_f
        character, dimension(:,:), allocatable, target, save:: f


        call c_f_string(param_name, param_name_f)
        call c_f_string(section, section_f)

        call get_string_array_param(section_f, param_name_f, string_array_param_f)
        !allocate(f(1, 5))
        !f(1,1) = 'H'
        !f(1,2) = 'O'
        !f(1,3) = 'L'
        !f(1,4) = 'A'
        !f(1,5) = C_NULL_CHAR

        print*, "ENTRA"

        call f_c_string_array(string_array_param_f, string_array_param)

        array_len = size(string_array_param_f)
        print*, "SALE"
    end subroutine CPLC_get_string_array_param
#endif


    subroutine CPLC_init(calling_realm, returned_realm_comm) &
        bind (C, name="CPLC_init_Fort")
        use CPL, only: CPL_init
        implicit none
        
        integer(C_INT), value :: calling_realm
        integer(C_INT)        :: returned_realm_comm
        
        integer :: ierror

        call CPL_init(calling_realm, returned_realm_comm, ierror)

    end subroutine CPLC_init

    subroutine CPLC_finalize() &
        bind (C, name="CPLC_finalize")
        use CPL, only: CPL_finalize
        implicit none
               
        integer :: ierror

        call CPL_finalize(ierror)

    end subroutine CPLC_finalize

    subroutine CPLC_setup_cfd(icomm_grid, xyzL, xyz_orig, ncxyz) &
        bind (C, name="CPLC_setup_cfd_Fort")
        use CPL, only: CPL_setup_cfd
        implicit none

        ! Integers
        integer(C_INT), value :: icomm_grid
        type(C_PTR), value :: ncxyz ! (3)

        ! Reals
        type(C_PTR), value :: xyzL, xyz_orig ! (3)

        ! Fortran equivalent array pointers
        integer, dimension(:), pointer :: ncxyz_f
        real(kind(0.d0)), dimension(:), pointer :: xyzL_f, xyz_orig_f

        !n Integers
        call C_F_POINTER(ncxyz, ncxyz_f, [3])

        ! Reals
        call C_F_POINTER(xyzL, xyzL_f, [3])
        call C_F_POINTER(xyz_orig, xyz_orig_f, [3])

        call CPL_setup_cfd(icomm_grid, xyzL_f, xyz_orig_f, ncxyz_f)


    end subroutine CPLC_setup_cfd


    subroutine CPLC_cfd_init(icomm_grid, icoord, npxyz_cfd, xyzL, xyz_orig, &
                             ncxyz, ijkcmax, ijkcmin, iTmin, iTmax,  &
                             jTmin, jTmax, kTmin, kTmax, xgrid, ygrid, zgrid) &
        bind (C, name="CPLC_cfd_init")
        use CPL, only: coupler_cfd_init
        implicit none

        ! Integers
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
        type(C_PTR), value :: xyzL ! (3)
        type(C_PTR), value :: xyz_orig ! (3)
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
        real(kind(0.d0)), dimension(:), pointer :: xyz_orig_f
        real(kind(0.d0)), dimension(:,:,:), pointer :: xgrid_f, ygrid_f, zgrid_f

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
        call C_F_POINTER(xyz_orig, xyz_orig_f, [3])
        call C_F_POINTER(xgrid, xgrid_f, [ncx+1, ncy+1, ncz+1])
        call C_F_POINTER(ygrid, ygrid_f, [ncx+1, ncy+1, ncz+1])
        call C_F_POINTER(zgrid, zgrid_f, [ncx+1, ncy+1, ncz+1])

        ! Call the old routine
        call coupler_cfd_init(icomm_grid, icoord_f, npxyz_cfd_f,  &
                              xyzL_f, xyz_orig_f, ncxyz_f, ijkcmax_f, ijkcmin_f, &
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


    subroutine CPLC_setup_md(icomm_grid, xyzL, xyz_orig) &
        bind(C, name="CPLC_setup_md_Fort")
        use CPL, only: CPL_setup_md
        implicit none

        ! Integers
        integer(C_INT), value :: icomm_grid
       
        ! Reals
        type(C_PTR), value :: xyzL, xyz_orig! (3)

        ! Fortran equivalent arrays
        real(kind(0.d0)), dimension(:), pointer :: xyzL_f, xyz_orig_f

        ! C -> Fortran array pointer type conversions
        call C_F_POINTER(xyzL, xyzL_f, [3])
        call C_F_POINTER(xyz_orig, xyz_orig_f, [3])

        ! Call routine
        call CPL_setup_md(icomm_grid, xyzL_f, xyz_orig_f)
    
        
    end subroutine CPLC_setup_md


    subroutine CPLC_md_init(icomm_grid, icoord, & 
                            npxyz_md, globaldomain, xyz_orig) &
        bind(C, name="CPLC_md_init")
        use CPL, only: coupler_md_init
        implicit none
       
        ! Integers
        integer(C_INT), value :: icomm_grid
        type(C_PTR), value :: icoord ! (3, nprocs)
        type(C_PTR), value :: npxyz_md ! (3)
       
        ! Reals
        type(C_PTR), value :: globaldomain ! (3)
        type(C_PTR), value :: xyz_orig! (3)

        ! Fortran equivalent arrays
        integer, dimension(:), pointer :: npxyz_md_f
        integer, dimension(:,:), pointer :: icoord_f
        real(kind(0.d0)), dimension(:), pointer :: globaldomain_f
        real(kind(0.d0)), dimension(:), pointer :: xyz_orig_f 

        ! Other function internals
        integer :: nprocs

        ! C -> Fortran array pointer type conversions
        call C_F_POINTER(npxyz_md, npxyz_md_f, [3])
        call C_F_POINTER(globaldomain, globaldomain_f, [3])
        call C_F_POINTER(xyz_orig, xyz_orig_f, [3])

        ! Store convenience variables
        nprocs = npxyz_md_f(1)*npxyz_md_f(2)*npxyz_md_f(3)
   
        ! Remaining conversions
        call C_F_POINTER(icoord, icoord_f, [3, nprocs])

        ! Call routine
        call coupler_md_init(icomm_grid, icoord_f, npxyz_md_f, & 
                             globaldomain_f, xyz_orig_f)
    
    end subroutine CPLC_md_init


    subroutine CPLC_set_timing(initialstep, Nsteps, dt) &
        bind(C, name="CPLC_set_timing")
        use CPL, only: CPL_set_timing
        implicit none

        ! Integer
        integer(C_INT), value :: Nsteps, initialstep

        ! Real
        real(C_DOUBLE), value :: dt 

        call CPL_set_timing(initialstep, Nsteps, dt)

    end subroutine CPLC_set_timing


    subroutine CPLC_send(asend, asend_shape, limits, send_flag) &
        bind(C, name="CPLC_send")
        use CPL, only: CPL_send!, CPL_send_4d, CPL_send_3d
        implicit none

        ! Boolean
        logical(C_BOOL),intent(out) :: send_flag

        ! Inputs
        type(C_PTR), value :: asend
        type(C_PTR), value :: asend_shape
        type(C_PTR), value :: limits

        ! Fortran equivalent array pointers
        logical :: send_flag_f
        real(kind(0.d0)), dimension(:,:,:,:), pointer :: asend_f
        integer, dimension(:), pointer :: limits_f
        integer, dimension(:), pointer :: asend_shape_f

        ! Other useful variables internal to this subroutine
        integer :: s1, s2, s3, s4

        call C_F_POINTER(limits, limits_f, [6])
        limits_f = limits_f + 1

        call C_F_POINTER(asend_shape, asend_shape_f, [4])
        s1 = asend_shape_f(1) 
        s2 = asend_shape_f(2) 
        s3 = asend_shape_f(3) 
        s4 = asend_shape_f(4)
        call C_F_POINTER(asend, asend_f, [s1, s2, s3, s4])
        call C_F_BOOL(send_flag, send_flag_f)
        call CPL_send(asend_f, limits_f, send_flag_f)
        call F_C_BOOL(send_flag_f,send_flag)
        limits_f = limits_f - 1
    
    end subroutine CPLC_send


    subroutine CPLC_recv(arecv, arecv_shape, limits, recv_flag) &
        bind(C, name="CPLC_recv")
        use CPL, only: CPL_recv
        implicit none

        ! Boolean
        logical(C_BOOL),intent(out) :: recv_flag

        ! Inputs
        type(C_PTR), value :: arecv
        type(C_PTR), value :: limits
        type(C_PTR), value :: arecv_shape

        ! Fortran equivalent array pointers
        logical :: recv_flag_f
        real(kind(0.d0)), dimension(:,:,:,:), pointer :: arecv_f
        integer, dimension(:), pointer :: limits_f
        integer, dimension(:), pointer :: arecv_shape_f

        ! Other useful variables internal to this subroutine
        integer :: s1, s2, s3, s4

        call C_F_POINTER(limits, limits_f, [6])
        limits_f = limits_f + 1

        call C_F_POINTER(arecv_shape, arecv_shape_f, [4])
        s1 = arecv_shape_f(1) 
        s2 = arecv_shape_f(2) 
        s3 = arecv_shape_f(3) 
        s4 = arecv_shape_f(4) 
        call C_F_POINTER(arecv, arecv_f, [s1, s2, s3, s4])

        call C_F_BOOL(recv_flag, recv_flag_f)
        call CPL_recv(arecv_f, limits_f, recv_flag_f)
        call F_C_BOOL(recv_flag_f,recv_flag)
        limits_f = limits_f - 1

    end subroutine CPLC_recv

    subroutine CPLC_send_min(asend, asend_shape, send_flag) &
        bind(C, name="CPLC_send_min")
        use CPL, only: CPL_send!, CPL_send_4d, CPL_send_3d
        implicit none

        ! Boolean
        logical(C_BOOL),intent(out) :: send_flag

        ! Inputs
        type(C_PTR), value :: asend
        type(C_PTR), value :: asend_shape

        ! Fortran equivalent array pointers
        logical :: send_flag_f
        real(kind(0.d0)), dimension(:,:,:,:), pointer :: asend_f
        integer, dimension(:), pointer :: asend_shape_f

        ! Other useful variables internal to this subroutine
        integer :: s1, s2, s3, s4

        call C_F_POINTER(asend_shape, asend_shape_f, [4])
        s1 = asend_shape_f(1) 
        s2 = asend_shape_f(2) 
        s3 = asend_shape_f(3) 
        s4 = asend_shape_f(4)
        call C_F_POINTER(asend, asend_f, [s1, s2, s3, s4])
        call C_F_BOOL(send_flag, send_flag_f)
        call CPL_send(asend_f, send_flag_f)
        call F_C_BOOL(send_flag_f,send_flag)
    
    end subroutine CPLC_send_min

    subroutine CPLC_recv_min(arecv, arecv_shape, recv_flag) &
        bind(C, name="CPLC_recv_min")
        use CPL, only: CPL_recv
        implicit none

        ! Boolean
        logical(C_BOOL),intent(out) :: recv_flag

        ! Inputs
        type(C_PTR), value :: arecv
        type(C_PTR), value :: arecv_shape

        ! Fortran equivalent array pointers
        logical :: recv_flag_f
        real(kind(0.d0)), dimension(:,:,:,:), pointer :: arecv_f
        integer, dimension(:), pointer :: arecv_shape_f

        ! Other useful variables internal to this subroutine
        integer :: s1, s2, s3, s4

        call C_F_POINTER(arecv_shape, arecv_shape_f, [4])
        s1 = arecv_shape_f(1) 
        s2 = arecv_shape_f(2) 
        s3 = arecv_shape_f(3) 
        s4 = arecv_shape_f(4) 
        call C_F_POINTER(arecv, arecv_f, [s1, s2, s3, s4])

        call C_F_BOOL(recv_flag, recv_flag_f)
        call CPL_recv(arecv_f, recv_flag_f)
        call F_C_BOOL(recv_flag_f,recv_flag)

    end subroutine CPLC_recv_min


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

        limits_f = limits_f + 1

        ! Library call
        call CPL_scatter(scatterarray_f, sn, limits_f, recvarray_f)

        limits_f = limits_f - 1



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

        limits_f = limits_f + 1

        ! Library call
        call CPL_gather(gatherarray_f, rn, limits_f, recvarray_f)

        limits_f = limits_f - 1


    end subroutine CPLC_gather


    subroutine CPLC_swaphalos(A, A_shape) &
        bind(C, name="CPLC_swaphalos")
        use CPL, only: CPL_swaphalos
        implicit none

        ! Inputs
        type(C_PTR), value :: A
        type(C_PTR), value :: A_shape

        ! Fortran equivalent array pointers
        real(kind(0.d0)), dimension(:,:,:,:), pointer :: A_f
        integer, dimension(:), pointer :: A_shape_f

        ! Other useful variables internal to this subroutine
        integer :: s1, s2, s3, s4

        call C_F_POINTER(A_shape, A_shape_f, [4])
        s1 = A_shape_f(1) 
        s2 = A_shape_f(2) 
        s3 = A_shape_f(3) 
        s4 = A_shape_f(4)
        call C_F_POINTER(A, A_f, [s1, s2, s3, s4])
        call CPL_swaphalos(A_f)  

    end subroutine CPLC_swaphalos


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

        coord_f = coord_f + 1

        call CPL_proc_extents(coord_f, realm, extents_f)

        coord_f = coord_f - 1
        extents_f = extents_f - 1;

    end subroutine CPLC_proc_extents


    subroutine CPLC_my_proc_extents(extents) &
        bind(C, name="CPLC_my_proc_extents")
        use CPL, only: CPL_my_proc_extents
        
        ! Inputs
        type(C_PTR), value :: extents! (6)

        ! Fortran equivalent arrays
        integer, dimension(:), pointer :: extents_f

        call C_F_POINTER(extents, extents_f, [6])

        call CPL_my_proc_extents(extents_f)

        extents_f = extents_f - 1;

    end subroutine CPLC_my_proc_extents


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


        coord_f = coord_f + 1
        limits_f = limits_f + 1


        call CPL_proc_portion(coord_f, realm, limits_f, portion_f)

        portion_f = portion_f - 1;
        coord_f = coord_f - 1
        limits_f = limits_f - 1

    end subroutine CPLC_proc_portion


    subroutine CPLC_my_proc_portion(limits, portion) &
        bind(C, name="CPLC_my_proc_portion")
        use CPL, only: CPL_my_proc_portion
        
        ! Inputs
        type(C_PTR), value :: limits ! (6)
        type(C_PTR), value :: portion ! (6)

        ! Fortran equivalent arrays
        integer, dimension(:), pointer :: limits_f
        integer, dimension(:), pointer :: portion_f

        call C_F_POINTER(limits, limits_f, [6])
        call C_F_POINTER(portion, portion_f, [6])

        limits_f = limits_f + 1

        call CPL_my_proc_portion(limits_f, portion_f)

        portion_f = portion_f - 1;
        limits_f = limits_f - 1

    end subroutine CPLC_my_proc_portion

    logical(C_BOOL) function CPLC_is_proc_inside(region) &
        bind(C, name="CPLC_is_proc_inside")
        use CPL, only: CPL_is_proc_inside

        type(C_PTR), value :: region ! (6)
        integer, dimension(:), pointer :: region_f
        call C_F_POINTER(region, region_f, [6])

        region_f = region_f + 1

        CPLC_is_proc_inside = CPL_is_proc_inside(region_f)

        region_f = region_f - 1
    end function CPLC_is_proc_inside
        

    logical(C_BOOL) function CPLC_map_cfd2md_coord(coord_cfd, coord_md) &
        bind(C, name="CPLC_map_cfd2md_coord")
        use CPL, only: CPL_map_cfd2md_coord
    
        ! Input position
        type(C_PTR), value :: coord_cfd ! (3)
        type(C_PTR), value :: coord_md ! (3)
       
        ! Fortran equivalent array
        real(kind(0.d0)), dimension(:), pointer :: coord_cfd_f
        real(kind(0.d0)), dimension(:), pointer :: coord_md_f
       
        call C_F_POINTER(coord_cfd, coord_cfd_f, [3])
        call C_F_POINTER(coord_md, coord_md_f, [3])

        CPLC_map_cfd2md_coord = CPL_map_cfd2md_coord(coord_cfd_f, coord_md_f)

    end function CPLC_map_cfd2md_coord
    
    logical(C_BOOL) function CPLC_map_md2cfd_coord(coord_md, coord_cfd) &
        bind(C, name="CPLC_map_md2cfd_coord")
        use CPL, only: CPL_map_md2cfd_coord
    
        ! Input position
        type(C_PTR), value :: coord_cfd ! (3)
        type(C_PTR), value :: coord_md ! (3)
       
        ! Fortran equivalent array
        real(kind(0.d0)), dimension(:), pointer :: coord_cfd_f
        real(kind(0.d0)), dimension(:), pointer :: coord_md_f
       
        call C_F_POINTER(coord_cfd, coord_cfd_f, [3])
        call C_F_POINTER(coord_md, coord_md_f, [3])

        CPLC_map_md2cfd_coord = CPL_map_md2cfd_coord(coord_md_f, coord_cfd_f)

    end function CPLC_map_md2cfd_coord
 

    subroutine CPLC_map_cell2coord(i, j, k, coord_xyz) &
        bind(C, name="CPLC_map_cell2coord")
        use CPL, only: CPL_map_cell2coord
    
        ! Input position
        integer(C_INT), value :: i, j, k
        type(C_PTR), value :: coord_xyz ! (3)
       
        ! Fortran equivalent array
        real(kind(0.d0)), dimension(:), pointer :: coord_xyz_f
       
        call C_F_POINTER(coord_xyz, coord_xyz_f, [3])
        
        call CPL_map_cell2coord(i + 1, j + 1, k + 1, coord_xyz_f)
        
    end subroutine CPLC_map_cell2coord


    logical(C_BOOL) function CPLC_map_coord2cell(x, y, z, cell_ijk) &
        bind(C, name="CPLC_map_coord2cell")
        use CPL, only: CPL_map_coord2cell
    
        ! Input position
        real(C_DOUBLE), value :: x, y, z
        type(C_PTR), value :: cell_ijk ! (3)
       
        ! Fortran equivalent array
        integer, dimension(:), pointer :: cell_ijk_f
       
        call C_F_POINTER(cell_ijk, cell_ijk_f, [3])
       
        CPLC_map_coord2cell = CPL_map_coord2cell(x, y, z, cell_ijk_f)

        cell_ijk_f = cell_ijk_f - 1
        
    end function CPLC_map_coord2cell
    

    subroutine CPLC_get_no_cells(limits, no_cells) &
        bind(C, name="CPLC_get_no_cells")
        use CPL, only: CPL_get_no_cells
    
        ! Input position
        type(C_PTR), value :: limits ! (6)
        type(C_PTR), value :: no_cells ! (3)
       
        ! Fortran equivalent array
        integer, dimension(:), pointer :: limits_f
        integer, dimension(:), pointer :: no_cells_f
       
        call C_F_POINTER(limits, limits_f, [6])
        call C_F_POINTER(no_cells, no_cells_f, [3])

        limits_f = limits_f + 1

        call CPL_get_no_cells(limits_f, no_cells_f)

        limits_f = limits_f - 1

    end subroutine CPLC_get_no_cells

    logical(C_BOOL) function CPLC_map_glob2loc_cell(limits, glob_cell, loc_cell) &
        bind(C, name="CPLC_map_glob2loc_cell")
        use CPL, only: CPL_map_glob2loc_cell
    
        ! Input position
        type(C_PTR), value :: limits ! (6)
        type(C_PTR), value :: glob_cell, loc_cell ! (3)
       
        ! Fortran equivalent array
        integer, dimension(:), pointer :: limits_f
        integer, dimension(:), pointer :: glob_cell_f, loc_cell_f
       
        call C_F_POINTER(limits, limits_f, [6])
        call C_F_POINTER(glob_cell, glob_cell_f, [3])
        call C_F_POINTER(loc_cell, loc_cell_f, [3])

        limits_f = limits_f + 1
        glob_cell_f = glob_cell_f + 1

        CPLC_map_glob2loc_cell = CPL_map_glob2loc_cell(limits_f, glob_cell_f, loc_cell_f)

        limits_f = limits_f - 1
        glob_cell_f = glob_cell_f - 1
        loc_cell_f = loc_cell_f - 1

    end function CPLC_map_glob2loc_cell


    logical(C_BOOL) function CPLC_map_loc2glob_cell(limits, loc_cell, glob_cell) &
        bind(C, name="CPLC_map_loc2glob_cell")
        use CPL, only: CPL_map_loc2glob_cell
    
        ! Input position
        type(C_PTR), value :: limits ! (6)
        type(C_PTR), value :: glob_cell, loc_cell ! (3)
       
        ! Fortran equivalent array
        integer, dimension(:), pointer :: limits_f
        integer, dimension(:), pointer :: glob_cell_f, loc_cell_f
       
        call C_F_POINTER(limits, limits_f, [6])
        call C_F_POINTER(glob_cell, glob_cell_f, [3])
        call C_F_POINTER(loc_cell, loc_cell_f, [3])

        limits_f = limits_f + 1
        glob_cell_f = glob_cell_f + 1

        CPLC_map_loc2glob_cell = CPL_map_loc2glob_cell(limits_f, loc_cell_f, glob_cell_f)

        limits_f = limits_f - 1
        glob_cell_f = glob_cell_f - 1
        loc_cell_f = loc_cell_f - 1

    end function CPLC_map_loc2glob_cell
   

    subroutine CPLC_get_olap_limits(limits) &
        bind(C, name="CPLC_get_olap_limits")
        use CPL, only: CPL_get_olap_limits
    
        ! Input position
        type(C_PTR), value :: limits ! (6)
       
        ! Fortran equivalent array
        integer, dimension(:), pointer :: limits_f
       
        call C_F_POINTER(limits, limits_f, [6])


        call CPL_get_olap_limits(limits_f)

        limits_f = limits_f - 1

    end subroutine CPLC_get_olap_limits
   

    subroutine CPLC_get_cnst_limits(limits) &
        bind(C, name="CPLC_get_cnst_limits")
        use CPL, only: CPL_get_cnst_limits
    
        ! Input position
        type(C_PTR), value :: limits ! (6)
       
        ! Fortran equivalent array
        integer, dimension(:), pointer :: limits_f
       
        call C_F_POINTER(limits, limits_f, [6])

        call CPL_get_cnst_limits(limits_f)

        limits_f = limits_f - 1

    end subroutine CPLC_get_cnst_limits

    subroutine CPLC_get_bnry_limits(limits) &
        bind(C, name="CPLC_get_bnry_limits")
        use CPL, only: CPL_get_bnry_limits
    
        ! Input position
        type(C_PTR), value :: limits ! (6)
       
        ! Fortran equivalent array
        integer, dimension(:), pointer :: limits_f
       
        call C_F_POINTER(limits, limits_f, [6])

        call CPL_get_bnry_limits(limits_f)

        limits_f = limits_f - 1

    end subroutine CPLC_get_bnry_limits



    logical(C_BOOL) function CPLC_overlap() &
        bind(C, name="CPLC_overlap")
        use CPL, only: CPL_overlap
    
        call F_C_BOOL(CPL_overlap(), CPLC_overlap)

    end function CPLC_overlap


!    logical(C_BOOL) function CPLC_realm() &
!        bind(C, name="CPLC_realm")
!        use CPL, only: CPL_realm
!    
!        CPLC_overlap = CPL_realm()

!    end function CPLC_realm
    
   
    ! Setters:

    subroutine CPLC_set_output_mode(mode) &
        bind(C, name="CPLC_set_output_mode")
        use CPL, only: set_output_mode
        
        integer(C_INT), intent(in), value :: mode

        call set_output_mode(mode)
        
    end subroutine CPLC_set_output_mode


    ! Getters: integers

    integer(C_INT) function CPLC_nsteps_md() &
        bind(C, name="CPLC_nsteps_md")
        use CPL, only: nsteps_md
        implicit none

        CPLC_nsteps_md = nsteps_md
        
    end function CPLC_nsteps_md

    integer(C_INT) function CPLC_nsteps_cfd() &
        bind(C, name="CPLC_nsteps_cfd")
        use CPL, only: nsteps_cfd
        implicit none

        CPLC_nsteps_cfd = nsteps_cfd
        
    end function CPLC_nsteps_cfd

    integer(C_INT) function CPLC_nsteps_coupled() &
        bind(C, name="CPLC_nsteps_coupled")
        use CPL, only: nsteps_coupled
        implicit none

        CPLC_nsteps_coupled = nsteps_coupled
        
    end function CPLC_nsteps_coupled

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

        CPLC_icmin_olap = icmin_olap - 1
        
    end function CPLC_icmin_olap

    integer(C_INT) function CPLC_jcmin_olap() &
        bind(C, name="CPLC_jcmin_olap")
        use CPL, only: jcmin_olap
        implicit none

        CPLC_jcmin_olap = jcmin_olap - 1
        
    end function CPLC_jcmin_olap

    integer(C_INT) function CPLC_kcmin_olap() &
        bind(C, name="CPLC_kcmin_olap")
        use CPL, only: kcmin_olap
        implicit none

        CPLC_kcmin_olap = kcmin_olap - 1
        
    end function CPLC_kcmin_olap

    integer(C_INT) function CPLC_icmax_olap() &
        bind(C, name="CPLC_icmax_olap")
        use CPL, only: icmax_olap
        implicit none

        CPLC_icmax_olap = icmax_olap - 1
        
    end function CPLC_icmax_olap

    integer(C_INT) function CPLC_jcmax_olap() &
        bind(C, name="CPLC_jcmax_olap")
        use CPL, only: jcmax_olap
        implicit none

        CPLC_jcmax_olap = jcmax_olap - 1
        
    end function CPLC_jcmax_olap

    integer(C_INT) function CPLC_kcmax_olap() &
        bind(C, name="CPLC_kcmax_olap")
        use CPL, only: kcmax_olap
        implicit none

        CPLC_kcmax_olap = kcmax_olap - 1
        
    end function CPLC_kcmax_olap

    integer(C_INT) function CPLC_icmin_cnst() &
        bind(C, name="CPLC_icmin_cnst")
        use CPL, only: icmin_cnst
        implicit none

        CPLC_icmin_cnst = icmin_cnst - 1
        
    end function CPLC_icmin_cnst

    integer(C_INT) function CPLC_jcmin_cnst() &
        bind(C, name="CPLC_jcmin_cnst")
        use CPL, only: jcmin_cnst
        implicit none

        CPLC_jcmin_cnst = jcmin_cnst - 1
        
    end function CPLC_jcmin_cnst

    integer(C_INT) function CPLC_kcmin_cnst() &
        bind(C, name="CPLC_kcmin_cnst")
        use CPL, only: kcmin_cnst
        implicit none

        CPLC_kcmin_cnst = kcmin_cnst - 1
        
    end function CPLC_kcmin_cnst

    integer(C_INT) function CPLC_icmax_cnst() &
        bind(C, name="CPLC_icmax_cnst")
        use CPL, only: icmax_cnst
        implicit none

        CPLC_icmax_cnst = icmax_cnst - 1
        
    end function CPLC_icmax_cnst

    integer(C_INT) function CPLC_jcmax_cnst() &
        bind(C, name="CPLC_jcmax_cnst")
        use CPL, only: jcmax_cnst
        implicit none

        CPLC_jcmax_cnst = jcmax_cnst - 1
        
    end function CPLC_jcmax_cnst

    integer(C_INT) function CPLC_kcmax_cnst() &
        bind(C, name="CPLC_kcmax_cnst")
        use CPL, only: kcmax_cnst
        implicit none

        CPLC_kcmax_cnst = kcmax_cnst - 1
        
    end function CPLC_kcmax_cnst


    integer(C_INT) function CPLC_icmin_bnry() &
        bind(C, name="CPLC_icmin_bnry")
        use CPL, only: icmin_bnry
        implicit none

        CPLC_icmin_bnry = icmin_bnry - 1
        
    end function CPLC_icmin_bnry

    integer(C_INT) function CPLC_jcmin_bnry() &
        bind(C, name="CPLC_jcmin_bnry")
        use CPL, only: jcmin_bnry
        implicit none

        CPLC_jcmin_bnry = jcmin_bnry - 1
        
    end function CPLC_jcmin_bnry

    integer(C_INT) function CPLC_kcmin_bnry() &
        bind(C, name="CPLC_kcmin_bnry")
        use CPL, only: kcmin_bnry
        implicit none

        CPLC_kcmin_bnry = kcmin_bnry - 1
        
    end function CPLC_kcmin_bnry

    integer(C_INT) function CPLC_icmax_bnry() &
        bind(C, name="CPLC_icmax_bnry")
        use CPL, only: icmax_bnry
        implicit none

        CPLC_icmax_bnry = icmax_bnry - 1
        
    end function CPLC_icmax_bnry

    integer(C_INT) function CPLC_jcmax_bnry() &
        bind(C, name="CPLC_jcmax_bnry")
        use CPL, only: jcmax_bnry
        implicit none

        CPLC_jcmax_bnry = jcmax_bnry - 1
        
    end function CPLC_jcmax_bnry

    integer(C_INT) function CPLC_kcmax_bnry() &
        bind(C, name="CPLC_kcmax_bnry")
        use CPL, only: kcmax_bnry
        implicit none

        CPLC_kcmax_bnry = kcmax_bnry - 1
        
    end function CPLC_kcmax_bnry

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

    real(C_DOUBLE) function CPLC_x_orig_cfd() &
        bind(C, name="CPLC_x_orig_cfd")
        use CPL, only: x_orig_cfd
        implicit none

        CPLC_x_orig_cfd = x_orig_cfd
        
    end function CPLC_x_orig_cfd

    real(C_DOUBLE) function CPLC_y_orig_cfd() &
        bind(C, name="CPLC_y_orig_cfd")
        use CPL, only: y_orig_cfd
        implicit none

        CPLC_y_orig_cfd = y_orig_cfd
        
    end function CPLC_y_orig_cfd

    real(C_DOUBLE) function CPLC_z_orig_cfd() &
        bind(C, name="CPLC_z_orig_cfd")
        use CPL, only: z_orig_cfd
        implicit none

        CPLC_z_orig_cfd = z_orig_cfd
        
    end function CPLC_z_orig_cfd

    real(C_DOUBLE) function CPLC_x_orig_md() &
        bind(C, name="CPLC_x_orig_md")
        use CPL, only: x_orig_md
        implicit none

        CPLC_x_orig_md = x_orig_md
        
    end function CPLC_x_orig_md

    real(C_DOUBLE) function CPLC_y_orig_md() &
        bind(C, name="CPLC_y_orig_md")
        use CPL, only: y_orig_md
        implicit none

        CPLC_y_orig_md = y_orig_md
        
    end function CPLC_y_orig_md

    real(C_DOUBLE) function CPLC_z_orig_md() &
        bind(C, name="CPLC_z_orig_md")
        use CPL, only: z_orig_md
        implicit none

        CPLC_z_orig_md = z_orig_md
        
    end function CPLC_z_orig_md


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
