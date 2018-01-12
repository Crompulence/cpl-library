!=============================================================================

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
    use commondefs

    implicit none


    integer, protected :: fd_counter = 0 
    character(len=CPL_STRING_MAX_LEN), protected  :: param_fname
    type(json_file), private :: json

    public CPL_load_param_file

    interface get_file_param
        module procedure get_int_param, get_boolean_param, get_string_param, get_real_param,&
        get_real_array_param, get_int_array_param 
    end interface get_file_param

    public get_int_param, get_boolean_param, get_string_param, get_real_param, get_real_array_param,&
           get_int_array_param 

    private bcast_param_file, read_line


contains

function read_line(aunit, InLine, trimmed) result(OK)
    integer, intent(IN) :: aunit
    character(LEN=:), allocatable, optional :: InLine
    logical, intent(in), optional :: trimmed
    integer, parameter :: line_buf_len= 1024*4
    character(LEN=line_buf_len) :: InS
    logical :: OK, set
    integer status, size

    OK = .false.
    set = .true.
    do
        read (aunit,'(a)',advance='NO',iostat=status, size=size) InS
        OK = .not. IS_IOSTAT_END(status)
        if (.not. OK) return
        if (present(InLine)) then
            if (set) then
                InLine = InS(1:size)
                set=.false.
            else
                InLine = InLine // InS(1:size)
            end if
        end if
        if (IS_IOSTAT_EOR(status)) exit
    end do
    if (present(trimmed) .and. present(InLine)) then
        if (trimmed) InLine = trim(adjustl(InLine))
    end if
end function read_line 

function CPL_is_param_file_loaded() result(p)
    logical :: p

    if (fd_counter .gt. 0) p = .true.

endfunction CPL_is_param_file_loaded

subroutine CPL_load_param_file(fname)
    use mpi
    use coupler_module, only: error_abort, get_new_fileunit,&
                              CPL_REALM_COMM, md_realm, realm

    implicit none

    integer :: myrank, ierr
    integer :: unitno_in, unitno_out
    character(kind=json_CK, len=*), intent(in) :: fname
    character(len=255) :: fname_tmp
    character(LEN=:), allocatable :: InLine

    fd_counter = fd_counter + 1
    param_fname = fname
    call json%initialize()

    call MPI_comm_rank(CPL_REALM_COMM, myrank, ierr)

    if (myrank .eq. 0) then
        unitno_in = get_new_fileunit() 
        open (unitno_in, file=fname, action="read")
        ! FILTER COMMENTS --> Generated 'config.tmp' comment-free
        if (realm .eq. md_realm) then
            fname_tmp = "config_md.tmp"
        else
            fname_tmp = "config_cfd.tmp"
        endif
        unitno_out = get_new_fileunit() 
        open (unitno_out, file=fname_tmp, action="write")
        do while(read_line(unitno_in, InLine))
            InLine = ADJUSTL(InLine)
            if (.not. InLine(1:1) .eq. '#') then
                write(unitno_out,*)  InLine
            end if
        end do 
        close(unitno_in)
        close(unitno_out)
        ! FILTER COMMENTS END
        call json%load_file(filename = fname_tmp)
        ! Delete 'config.tmp'
        open (unitno_out, file=fname_tmp, status='old')
        close(unitno_out, status='delete')
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
    use coupler_module, only: error_abort, CPL_REALM_COMM
    
    implicit none

    character(kind=json_CK,len=:), allocatable :: json_string
    integer :: file_size, ierr, myrank

    call MPI_comm_rank(CPL_REALM_COMM, myrank, ierr)
    
    if (myrank .eq. 0) then
        call json%print_to_string(json_string)
        if (json%failed()) then
            print *, "Stop print_to_string."
            stop 1
        endif
        file_size = len(json_string)
    endif

    call MPI_Bcast(file_size, 1, MPI_INTEGER, 0, CPL_REALM_COMM, ierr)
    if (myrank .ne. 0) then 
        allocate(character(kind=json_CK,len=file_size)::json_string)
    endif

    call MPI_Bcast(json_string, file_size, MPI_CHAR, 0, CPL_REALM_COMM, ierr)
    call json%load_from_string(json_string)
    if (json%failed()) then
        print *, "Stop load_from_string."
        stop 1
    endif

endsubroutine bcast_param_file

subroutine get_real_param(section, var_name, ret)
    use coupler_module, only: error_abort

    implicit none

    character(len=*), intent(in) :: section, var_name
    real(kind(0.d0)), intent(out) :: ret

    call json%get(section//'.'//var_name, ret)
    if (json%failed()) then
        print *, "Stop in get_real_param()."
        stop 1
    endif

endsubroutine get_real_param

subroutine get_real_array_param(section, var_name, ret)
    use coupler_module, only: error_abort

    implicit none

    character(len=*), intent(in) :: section, var_name
    real(kind(0.d0)), dimension(:), allocatable, intent(out) :: ret


    call json%get(section//'.'//var_name, ret)
    if (json%failed()) then
        print *, "Stop in get_real_array()."
        stop 1
    endif

endsubroutine get_real_array_param

subroutine get_int_param(section, var_name, ret)
    use coupler_module, only: error_abort

    implicit none

    character(len=*), intent(in) :: section, var_name
    integer, intent(out) :: ret

    call json%get(section//'.'//var_name, ret)
    if (json%failed()) then
        print *, "Stop in get_int()."
        stop 1
    endif

endsubroutine get_int_param

subroutine get_int_array_param(section, var_name, ret)
    use coupler_module, only: error_abort

    use coupler_module, only: error_abort
    implicit none

    character(len=*), intent(in) :: section, var_name
    integer, dimension(:),allocatable,intent(out) :: ret

    call json%get(section//'.'//var_name, ret)
    if (json%failed()) then
        call json%print_error_message()
        call error_abort("CPL_get_int_array_param - Error getting integer array param.")
    endif


endsubroutine get_int_array_param


subroutine get_string_param(section, var_name, ret)
    use coupler_module, only: error_abort

    implicit none

    character(len=*), intent(in) :: section, var_name
    character(len=:), allocatable, intent(out) :: ret

    call json%get(section//'.'//var_name, ret)
    if (json%failed())  then
        call json%print_error_message()
        call error_abort("CPL_get_string_param - Error getting string param.")
    endif

endsubroutine get_string_param

subroutine get_boolean_param(section, var_name, ret)
    use coupler_module, only: error_abort

    implicit none

    character(len=*), intent(in) :: section, var_name
    logical, intent(out) :: ret

    call json%get(section//'.'//var_name, ret)
    if (json%failed()) then
        call json%print_error_message()
        call error_abort("CPL_get_boolean_param - Error getting boolean param.")
    endif

endsubroutine get_boolean_param

! TODO NOTE: Those are not working at the moment
! subroutine get_string_array_param(section, var_name, ret)
!
!     implicit none
!
!     character(len=*), intent(in) :: section, var_name
!     character(len=*),dimension(:),allocatable,intent(out) :: ret
!
!     call json%get(section//'.'//var_name, ret)
!     if (json%failed()) stop 1
!
! endsubroutine get_string_array_param
!
! subroutine get_boolean_array_param(section, var_name, ret)
!
!     implicit none
!
!     character(len=*), intent(in) :: section, var_name
!     logical, dimension(:), allocatable, intent(out) :: ret
!
!     call json%get(section//'.'//var_name, ret)
!     if (json%failed()) stop 1
!
! endsubroutine get_boolean_array_param
!

end module io
