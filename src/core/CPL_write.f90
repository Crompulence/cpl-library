
module coupler_write
    use mpi

    integer :: ierr, nd=3

    !Generic interface so write arrays can be used with both integers and reals
    interface CPL_write_arrays
        module procedure rwrite_arrays_1, rwrite_arrays!, iwrite_arrays_1, iwrite_arrays
    end interface CPL_write_arrays
    private  rwrite_arrays_1, rwrite_arrays!, iwrite_arrays_1, iwrite_arrays
    
contains


!====================================================================
!           Array writing subroutines
!--------------------------------------------------------------------

!! --- Write integer arrays ---

!!   1D array wrapper
!subroutine iwrite_arrays_1(temp, nresults, outfile, outstep, icomm_grid, limits, nhb_)
!    implicit none

!    integer, intent(in)                     :: nresults, outstep, icomm_grid
!    integer, intent(in), dimension(6)       :: limits
!    integer, dimension(:,:,:),intent(in)    :: temp
!    character(*),intent(in)                 :: outfile

!    integer, intent(in), dimension(3), optional     :: nhb_

!    integer, dimension(:,:,:,:),allocatable :: some_array

!    allocate(some_array(size(temp,1),size(temp,2),size(temp,3),1))
!    some_array(:,:,:,1) = temp(:,:,:)
!    call iwrite_arrays(some_array, nresults, outfile, outstep, icomm_grid, limits, nhb_)
!    deallocate(some_array)

!end subroutine iwrite_arrays_1

!!   Main iwrite Routine
!subroutine iwrite_arrays(some_array, nresults, outfile, outstep, icomm_grid, limits, nhb_)
!    use coupler_module, only : iblock => iblock_realm, & 
!                               jblock => jblock_realm, & 
!                               kblock => kblock_realm, error_abort
!    use coupler, only : cpl_my_proc_portion, CPL_get_no_cells

!    implicit none

!    integer, intent(in)                     :: nresults, outstep, icomm_grid
!    integer, intent(in), dimension(6)       :: limits
!    integer, dimension(:,:,:,:),intent(in)  :: some_array!(nx,ny,nz,nresults)
!    character(*),intent(in)                 :: outfile

!    integer, intent(in), dimension(3), optional     :: nhb_

!    integer                                 :: n, fh, ixyz
!    integer                                 :: MEM_FLAG = 0
!    integer                                 :: FILE_FLAG = 0
!    integer                                 :: int_size,datatype, nproc
!    integer                                 :: status(mpi_status_size)
!    integer                                 :: filetype, memtype
!    integer (kind=MPI_offset_kind)          :: offset, global_cnt
!    integer, dimension(3)                   :: gsizes, lsizes, memsizes, nhb
!    integer, dimension(3)                   :: global_indices, local_indices
!    integer, dimension(6)                   :: portion
!    integer, allocatable,dimension(:,:)     :: proc_lsizes 
!    integer, allocatable,dimension(:,:,:)   :: OutBuffer
!    character(200)                          :: outfile_t

!    datatype = MPI_INTEGER
!    call MPI_TYPE_SIZE(datatype, int_size, ierr)

!    !--------- DEFINE LIMITS (FILE & LOCAL SUBARRAY) -------
!    !  Note:  MPI assumes here that numbering starts from zero
!    !  Since numbering starts from (one), subtract (one) from every index
!    !-------------------------------------------------------
!    !write a separate output files for each timestep
!    call get_Timestep_FileName(outstep-1, outfile, outfile_t)
!    global_cnt = 0; offset = 0

!    !Error if outstep less than one
!    if (outstep .le. 0) then
!        print*, "Error -- outstep ", outstep, " for filename ", & 
!                trim(outfile) , " results in ", trim(outfile_t)
!        call error_abort("Error in write_arrays -- requires outstep > 1 ")
!    endif

!    !Global "sizes", i.e. bins
!    call CPL_get_no_cells(limits, gsizes)
!    !Local "size"
!    call CPL_my_proc_portion(limits, portion)
!    call CPL_get_no_cells(portion, lsizes)
!    local_indices(:) = (/  0  , 0 , 0 /)
!    !Calculate global_indices
!    global_indices(:)= (/  0  , 0 , 0 /)
!    call MPI_comm_size(icomm_grid, nproc, ierr)
!    allocate(proc_lsizes(3,nproc))
!    call globalGather(lsizes, proc_lsizes, nd, icomm_grid)
!    do ixyz=1,nd
!        if (sum(proc_lsizes(ixyz,:)) .ne. gsizes(ixyz)) then 
!            call error_abort("Error in CPL_write iwrite_arrays -- not all cells have procs writing")
!        endif
!    enddo
!    global_indices(1) = sum(proc_lsizes(1,1:iblock-1))
!    global_indices(2) = sum(proc_lsizes(2,1:jblock-1))
!    global_indices(3) = sum(proc_lsizes(3,1:kblock-1))
!    deallocate(proc_lsizes)

!    !Allocate ouput buffer
!    allocate(OutBuffer(lsizes(1),lsizes(2),lsizes(3)))
!    memsizes = lsizes

!    !Open file to write
!    call MPI_file_open(icomm_grid, outfile_t, MPI_MODE_WRONLY+MPI_MODE_CREATE, &
!                                MPI_INFO_NULL, fh, ierr)

!    do n =1,nresults
!        !Copy to outbuffer
!        OutBuffer =  some_array(1+nhb(1):lsizes(1)+nhb(1), & 
!                                1+nhb(2):lsizes(2)+nhb(2), & 
!                                1+nhb(3):lsizes(3)+nhb(3),n)
!        !Update array datatype and reset fileview to correct section of array
!        CALL Create_commit_fileview(gsizes,lsizes,global_indices, & 
!                                    offset,datatype,FILE_FLAG,filetype,fh)
!        !Update local array datatype to ignore halo cells
!        CALL Create_commit_subarray(memsizes,lsizes,local_indices,datatype,MEM_FLAG,memtype)
!        !Write to file
!        CALL MPI_file_write_all(fh, OutBuffer, 1, memtype, status, ierr)

!        !Calculate global cnt offset
!        global_cnt = global_cnt + gsizes(1)*gsizes(2)*gsizes(3)
!        offset = global_cnt * int_size
!    
!    enddo

!    deallocate(OutBuffer)
!    CALL MPI_FILE_CLOSE(fh, ierr)

!    !==========================================================
!    !     FREE DATA TYPES
!    !----------------------------------------------------------
!    call MPI_TYPE_FREE(memtype,ierr) ; MEM_FLAG = 0
!    call MPI_TYPE_FREE(filetype,ierr); FILE_FLAG = 0

!end subroutine iwrite_arrays


! --- Double precision arrays ---

!   1D array wrapper
subroutine rwrite_arrays_1(temp, nresults, outfile, outstep, icomm_grid, limits, nhb_)
    implicit none

    integer, intent(in)                             :: nresults,outstep, icomm_grid
    integer, intent(in), dimension(6)               :: limits
    real(kind(0.d0)), dimension(:,:,:),intent(in)   :: temp
    character(*),intent(in)                         :: outfile

    integer, intent(in), dimension(3), optional     :: nhb_

    real(kind(0.d0)), dimension(:,:,:,:),allocatable    :: some_array

    allocate(some_array(size(temp,1),size(temp,2),size(temp,3),1))
    some_array(:,:,:,1) = temp(:,:,:)
    call rwrite_arrays(some_array, nresults, outfile, outstep, icomm_grid, limits, nhb_)
    deallocate(some_array)

end subroutine rwrite_arrays_1

!   Main rwrite Routine
subroutine rwrite_arrays(some_array, nresults, outfile, outstep, icomm_grid, limits, nhb_)
    use coupler_module, only : iblock => iblock_realm, & 
                               jblock => jblock_realm, & 
                               kblock => kblock_realm, VOID, error_abort
    use coupler, only : cpl_my_proc_portion, CPL_get_no_cells
    implicit none

    integer, intent(in)                             :: nresults, outstep, icomm_grid
    integer, intent(in), dimension(6)               :: limits
    real(kind(0.d0)), dimension(:,:,:,:),intent(in) :: some_array
    character(*),intent(in)                         :: outfile

    integer, intent(in), dimension(3), optional     :: nhb_

    integer                             :: n, fh
    integer                             :: MEM_FLAG = 0
    integer                             :: FILE_FLAG = 0
    integer                             :: datatype
    integer                             :: nproc
    integer                             :: status(mpi_status_size)
    integer                             :: filetype, memtype, dp_size_temp
    integer (kind=MPI_offset_kind)      :: offset, global_cnt, dp_size
    integer                             :: write_comm, CPL_write_comm
    integer, dimension(6)               :: portion
    integer, dimension(3)               :: gsizes, lsizes, memsizes, nhb
    integer, dimension(3)               :: global_indices, local_indices, blockmin
    integer, dimension(:,:),allocatable :: proc_lsizes 
    real(kind(0.d0)), allocatable,dimension(:,:,:)      :: OutBuffer
    character(200)                      :: outfile_t

    !Number of halos for each process (assumed 0)
    if (present(nhb_)) then
        nhb = nhb_
    else
        nhb = (/ 0, 0, 0 /)
    endif

    !Global "sizes", i.e. bins
    call CPL_get_no_cells(limits, gsizes)
    !Local "size"
    call CPL_my_proc_portion(limits, portion)
    if (any(portion .eq. VOID)) then
        lsizes = 0
        write_comm = 0
    else
        call CPL_get_no_cells(portion, lsizes)
        write_comm = 1
    endif
    !We need to split between writing and non-writing procs
    call MPI_comm_split(icomm_grid, write_comm, 0, CPL_write_comm, ierr)
    if (write_comm .eq. 0 ) then
        call MPI_Comm_free(CPL_write_comm, ierr)
        return
    endif


    datatype = MPI_DOUBLE_PRECISION
    call MPI_TYPE_SIZE(datatype, dp_size_temp, ierr)
    dp_size = dp_size_temp

    !--------- DEFINE LIMITS (FILE & LOCAL SUBARRAY) -------
    !  Note:  MPI assumes here that numbering starts from zero
    !  Since numbering starts from (one), subtract (one) from every index
    !-------------------------------------------------------
    ! write a separate output files for each timestep
    call get_Timestep_FileName(outstep-1,outfile,outfile_t)
    global_cnt = 0; offset = 0

    !Error if outstep less than one
    if (outstep .le. 0) then
        print*, "Error -- outstep ", outstep, " for filename ", & 
                trim(outfile) , " results in ", trim(outfile_t)
        call error_abort("Error in write_arrays -- requires outstep > 1 ")
    endif

    !Goes into MPI_TYPE_CREATE_SUBARRAY
    local_indices(:) = (/  0  , 0 , 0 /)                             
    !Calculate global_indices, goes into MPI_TYPE_CREATE_SUBARRAY
    global_indices(:)= (/  0  , 0 , 0 /)                             
    !Number of bins on each processor 
    call MPI_comm_size(CPL_write_comm, nproc, ierr)
    allocate(proc_lsizes(3,nproc))                                   
    call globalGather(lsizes, proc_lsizes, nd, CPL_write_comm)

    !We assume writing block is a contigous set of processors
    !Between minimum processor
    blockmin = (/ iblock, jblock, kblock /)
    call globalMinVect(blockmin, 3, CPL_write_comm)

    !Seems to be a set of offsets in x,y,z
    global_indices(1) = sum(proc_lsizes(1,1:iblock-blockmin(1)))
    global_indices(2) = sum(proc_lsizes(2,1:jblock-blockmin(2)))
    global_indices(3) = sum(proc_lsizes(3,1:kblock-blockmin(3)))
    deallocate(proc_lsizes)

    !Allocate ouput buffer
    allocate(OutBuffer(lsizes(1),lsizes(2),lsizes(3)))               !Buffer that's written
    memsizes = lsizes                                                !Not sure why redefining?

    !Open file to write
    call MPI_file_open(CPL_write_comm, outfile_t, MPI_MODE_WRONLY+MPI_MODE_CREATE, &
                                MPI_INFO_NULL, fh, ierr)

    do n =1,nresults
        !Copy to outbuffer
        OutBuffer =  some_array(n,1+nhb(1):lsizes(1)+nhb(1), & 
                                  1+nhb(2):lsizes(2)+nhb(2), & 
                                  1+nhb(3):lsizes(3)+nhb(3))
        !Update array datatype and reset fileview to correct section of array
        CALL Create_commit_fileview(gsizes,lsizes,global_indices,offset,datatype,FILE_FLAG,filetype,fh)
        !Update local array datatype to ignore halo cells
        CALL Create_commit_subarray(memsizes,lsizes,local_indices,datatype,MEM_FLAG,memtype)
        !Write to file
        CALL MPI_file_write_all(fh, OutBuffer, 1, memtype, status, ierr)

        !Calculate global cnt offset
        global_cnt = global_cnt + gsizes(1)*gsizes(2)*gsizes(3)
        offset = global_cnt * dp_size
    
    enddo

    deallocate(OutBuffer)
    CALL MPI_FILE_CLOSE(fh, ierr)

    !==========================================================
    !     FREE DATA TYPES
    !----------------------------------------------------------
    CALL MPI_BARRIER(CPL_write_comm,IERR)
    call MPI_TYPE_FREE(memtype,ierr) ; MEM_FLAG = 0
    call MPI_TYPE_FREE(filetype,ierr); FILE_FLAG = 0
    call MPI_Comm_free(CPL_write_comm, ierr)

end subroutine rwrite_arrays

subroutine Create_commit_fileview(gsizes,lsizes,global_indices, & 
                                  offset,datatype,FILE_FLAG,filetype,fh)
    implicit none

    integer,intent(inout)                       :: filetype
    integer,intent(in)                          :: datatype,fh
    integer (kind=MPI_offset_kind),intent(in)   :: offset
    integer, dimension(3),intent(in)            :: gsizes, lsizes, global_indices
    integer,intent(inout)                       :: FILE_FLAG

    if (FILE_FLAG.eq.1) then
        CALL MPI_TYPE_FREE(filetype,ierr)
        FILE_FLAG = 0
    end if
    CALL MPI_TYPE_CREATE_SUBARRAY(nd, gsizes, lsizes, global_indices, &
                                    MPI_ORDER_FORTRAN, datatype,  filetype, ierr)
    CALL MPI_TYPE_COMMIT(filetype, ierr)
    CALL MPI_FILE_SET_VIEW(fh, offset, datatype, filetype, &
                            'native', MPI_INFO_NULL, ierr)
    FILE_FLAG = 1
    
end subroutine Create_commit_fileview


subroutine Create_commit_subarray(memsizes,lsizes,local_indices,datatype,MEM_FLAG,memtype)
    implicit none

    integer,intent(inout)           :: memtype
    integer,intent(in)              :: datatype
    integer, dimension(3),intent(in):: memsizes, lsizes, local_indices
    integer,intent(inout)           :: MEM_FLAG


    if (MEM_FLAG.eq.1) then
        CALL MPI_TYPE_FREE(memtype,ierr)
        MEM_FLAG = 0
    end if
    CALL MPI_TYPE_CREATE_SUBARRAY(nd, memsizes, lsizes, local_indices, &
            MPI_ORDER_FORTRAN, datatype,  memtype, ierr)
    CALL MPI_TYPE_COMMIT(memtype, ierr)
    MEM_FLAG = 1
    
end subroutine Create_commit_subarray


!------------------------------------------------------------------------------
!Pure fortran subroutine to return an updated filename by appending
!the current timestep to that file
subroutine get_Timestep_FileName(timestep,basename,filename)
        implicit none

        integer,intent(in)             :: timestep
        character(*),intent(in)     :: basename
        character(*),intent(out)    :: filename

        if(timestep.le.9                                 ) &
        write(filename,'(a,a7,i1)') trim(basename),'.000000',timestep
        if(timestep.ge.10      .and. timestep.le.99     ) &
        write(filename,'(a,a6,i2)') trim(basename),'.00000' ,timestep
        if(timestep.ge.100     .and. timestep.le.999    ) &
        write(filename,'(a,a5,i3)') trim(basename),'.0000'  ,timestep
        if(timestep.ge.1000    .and. timestep.le.9999   ) &
        write(filename,'(a,a4,i4)') trim(basename),'.000'   ,timestep
        if(timestep.ge.10000   .and. timestep.le.99999  ) &
        write(filename,'(a,a3,i5)') trim(basename),'.00'    ,timestep
        if(timestep.ge.100000  .and. timestep.le.999999 ) &
        write(filename,'(a,a2,i6)') trim(basename),'.0'     ,timestep
        if(timestep.ge.1000000 .and. timestep.le.9999999) &
        write(filename,'(a,a1,i7)') trim(basename),'.'      ,timestep

        !Remove any surplus blanks
        filename = trim(filename)

end subroutine get_Timestep_FileName

subroutine globalGather(A, B, na, in_comm)
    implicit none

    integer, intent(in) :: na, in_comm

    integer, intent(in) :: A(na)
    integer, intent(out), dimension(:,:), allocatable :: B

    integer :: nproc

    call MPI_comm_size(in_comm, nproc, ierr)
    allocate(B(na, nproc))
    call MPI_Allgather (A, na, MPI_INTEGER, B, na, &
                        MPI_INTEGER,in_comm, ierr)

end subroutine globalGather

subroutine globalMinVect(A, na, in_comm)
    implicit none

    integer, intent(in) :: na, in_comm

    integer, intent(inout) :: A(na)
    integer, allocatable :: buf(:)

    allocate(buf(na))
    call MPI_AllReduce (A, buf, na, MPI_INTEGER, &
                        MPI_MIN, in_comm, ierr)
    A = buf
    deallocate(buf)

    return
end subroutine globalMinVect

end module coupler_write
