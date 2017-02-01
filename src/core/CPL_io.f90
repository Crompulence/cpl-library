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
! .. codeauthor:: Edward Smith Novemeber 2011 to present
! .. codeauthor:: Eduardo Ramos Fernandez 2015 to present
! .. codeauthor:: David Trevelyan September 2012 to December 2015
! .. codeauthor:: Lucian Anton, November 2011  
!
!=============================================================================

module io
    use json_module

    implicit none


    integer, protected :: fd_counter = 0 
    character(len=256), protected  :: param_fname
    type(json_file), private :: json
       

    public CPL_load_param_file

    interface get_file_param
        module procedure get_integer, get_logical, get_string, get_real,&
        get_real_array, get_integer_array, get_logical_array, get_string_array
    end interface get_file_param

    private get_integer, get_logical, get_string, get_real, get_real_array,&
            get_integer_array, get_logical_array, get_string_array, bcast_param_file


contains

function CPL_is_param_file_loaded() result(p)
    logical :: p

    if (fd_counter .gt. 0) p = .true.

endfunction CPL_is_param_file_loaded

subroutine CPL_load_param_file(fname)
    use mpi
    use coupler_module, only: error_abort

    implicit none

    integer :: myrank, ierr
    character(kind=json_CK, len=*), intent(in) :: fname


    fd_counter = fd_counter + 1
    param_fname = fname
    call json%initialize()

    call MPI_comm_rank(MPI_COMM_WORLD, myrank, ierr)

    if (myrank .eq. 0) then
        call json%load_file(filename = fname)
        if (json%failed()) then
            call json%print_error_message()
            call error_abort("CPL_load_param - Error loading param file.")
        endif
    endif

    call bcast_param_file()

endsubroutine CPL_load_param_file

subroutine CPL_close_param_file()
    implicit none

    call json%destroy()

endsubroutine CPL_close_param_file

subroutine bcast_param_file()
    use mpi
    
    implicit none

    character(kind=json_CK,len=:), allocatable :: json_string
    integer :: file_size, ierr, myrank

    call MPI_comm_rank(MPI_COMM_WORLD, myrank, ierr)
    
    if (myrank .eq. 0) then
        call json%print_to_string(json_string)
        if (json%failed()) stop 1
        file_size = len(json_string)
    endif

    call MPI_Bcast(file_size, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
    if (myrank .ne. 0) then 
        allocate(character(kind=json_CK,len=file_size)::json_string)
    endif

    call MPI_Bcast(json_string, file_size, MPI_CHAR, 0, MPI_COMM_WORLD, ierr)
    !print *, json_string
    call json%load_from_string(json_string)
    if (json%failed()) stop 1

endsubroutine bcast_param_file

subroutine get_real(section, var_name, ret)

    implicit none

    character(len=*), intent(in) :: section, var_name
    real(kind(0.d0)), intent(out) :: ret

    call json%get(section//'.'//var_name, ret)
    if (json%failed()) stop 1

endsubroutine get_real

subroutine get_real_array(section, var_name, ret)

    implicit none

    character(len=*), intent(in) :: section, var_name
    real(kind(0.d0)), dimension(:), allocatable, intent(out) :: ret

    call json%get(section//'.'//var_name, ret)
    if (json%failed()) stop 1

endsubroutine get_real_array

subroutine get_integer(section, var_name, ret)

    implicit none

    character(len=*), intent(in) :: section, var_name
    integer, intent(out) :: ret

    call json%get(section//'.'//var_name, ret)
    if (json%failed()) stop 1

endsubroutine get_integer

subroutine get_integer_array(section, var_name, ret)

    implicit none

    character(len=*), intent(in) :: section, var_name
    integer, dimension(:),allocatable,intent(out) :: ret

    call json%get(section//'.'//var_name, ret)
    if (json%failed()) stop 1

endsubroutine get_integer_array


subroutine get_string(section, var_name, ret)

    implicit none

    character(len=*), intent(in) :: section, var_name
    character(len=:), allocatable, intent(out) :: ret

    call json%get(section//'.'//var_name, ret)
    if (json%failed()) stop 1

endsubroutine get_string

subroutine get_string_array(section, var_name, ret)

    implicit none

    character(len=*), intent(in) :: section, var_name
    character(len=*),dimension(:),allocatable,intent(out) :: ret

    call json%get(section//'.'//var_name, ret)
    if (json%failed()) stop 1

endsubroutine get_string_array

subroutine get_logical(section, var_name, ret)

    implicit none

    character(len=*), intent(in) :: section, var_name
    logical, intent(out) :: ret

    call json%get(section//'.'//var_name, ret)
    if (json%failed()) stop 1

endsubroutine get_logical

subroutine get_logical_array(section, var_name, ret)

    implicit none

    character(len=*), intent(in) :: section, var_name
    logical, dimension(:), allocatable, intent(out) :: ret

    call json%get(section//'.'//var_name, ret)
    if (json%failed()) stop 1

endsubroutine get_logical_array


end module io
