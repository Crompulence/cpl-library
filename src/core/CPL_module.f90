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
!! COUPLER MODULE:
!! A single coupler module for both codes - this contains the same information
!! on both md and cfd side
!!
!!  - Error handling
!!  - MPI communicators
!!  - Simulation realms
!!  - MPI processor IDs
!!  - Processor topologies
!!  - Processor cartesian coords
!!  - Global cell grid parameters
!!  - Processor cell ranges
!!  - Domain and cell dimensions
!!  - Positions of CFD grid lines
!!  - CFD to MD processor mapping
!!  - Simulation parameters
!!
!! The data is protected so only setup routines in this module can change it
!
!
!
!  SETUP SETUP SETUP SETUP SETUP SETUP SETUP SETUP SETUP SETUP SETUP SETUP
!  -----------------------------------------------------------------------
!
!! Setup routines which have access to coupler parameters
!!
!! - CPL_create_comm          (cfd+md)   splits MPI_COMM_WORLD, create inter -
!!                                       communicator between CFD and MD
!!
!! - CPL_create_map           (cfd+md)   creates correspondence maps between
!!                                       the CFD grid and MD domains
!!
!! - CPL_cfd_adjust_domain     (cfd)     adjust CFD tomain to an integer number
!!                                       FCC or similar MD initial layout
!
!! .. codeauthor:: David Trevelyan, Edward Smith
!! @see coupler
!=============================================================================

module coupler_module
    implicit none

    integer, parameter :: VOID=-666         !!VOID value for data initialisation
    integer, parameter :: cfd_realm = 1     !! CFD realm identifier
    integer, parameter :: md_realm  = 2     !! MD realm identifier
    character(len=*),parameter :: &
        realm_name(2) = (/ "CFD", "MD " /)  !! Used with realm identifier to get name

    !! Output options
    integer, parameter:: QUIET = 0
    integer, parameter:: NORMAL = 1

    !! Max
    integer, parameter :: maxgridsize = 2*1024**2 !! Maximum size to store global grid

    !! error codes
    integer, parameter :: & 
        COUPLER_ERROR_REALM  = 1        !! wrong realm value
    integer, parameter :: & 
        COUPLER_ERROR_ONE_REALM  = 2    !! one realm missing
    integer, parameter :: & 
        COUPLER_ERROR_INIT       = 3    !! initialisation error
    integer, parameter :: & 
        COUPLER_ERROR_INPUT_FILE = 4    !! wrong value in input file
    integer, parameter :: & 
        COUPLER_ERROR_READ_INPUT = 5    !! error in processing input file or data transfers
    integer, parameter :: & 
        COUPLER_ERROR_CONTINUUM_FORCE = 6   !!the region in which the continuum constrain force is apply spans over two MD domains
    integer, parameter :: & 
        COUPLER_ABORT_ON_REQUEST = 7    !! used in request_abort 
    integer, parameter :: & 
        COUPLER_ABORT_SEND_CFD   = 8    !! error in coupler_cfd_send
    integer, parameter :: & 
        COUPLER_ERROR_CART_COMM   = 9   !! Wrong comm value in CPL_Cart_coords
    integer, parameter :: & 
        COUPLER_ERROR_SETUP_INCOMPLETE   = 10   !! CPL_setup_md or CPL_setup_cfd

    !! Output mode flag
    integer :: output_mode = NORMAL

    !! MPI error flag
    integer     :: ierr 

    ! MPI Communicators
    integer,protected :: &
        CPL_WORLD_COMM !! Both CFD and MD realms;
    integer,protected :: &
        CPL_REALM_COMM !! INTRA communicators within MD/CFD realms;
    integer,protected :: &
        CPL_INTER_COMM !!  CFD/MD INTER communicator between realm comms;
    integer,protected :: &
        CPL_CART_COMM !!  Comm w/cartesian topology for each realm;
    integer,protected :: &
        CPL_OLAP_COMM !!  Local comm between only overlapping MD/CFD procs;
    integer,protected :: &
        CPL_REALM_OLAP_COMM !!  Local CFD/MD in overlapping region;
    integer,protected :: &
        CPL_GRAPH_COMM !!  Comm w/ graph topolgy between locally olapg procs;
    integer,protected :: &
        CPL_REALM_INTERSECTION_COMM !!  Intersecting MD/CFD procs in world;

    !! Simulation realms
    integer,protected :: &
        realm

    ! MPI processor IDs
    integer,protected :: &
        myid_world   !! Processor ID from 0 to nproc_world-1;
    integer,protected :: &
        rank_world   !! Processor rank from 1 to nproc_world;
    integer,protected :: &
        rootid_world !! Root processor in world;
    integer,protected :: &
        myid_realm   !! Processor ID from 0 to nproc_realm-1;
    integer,protected :: &
        rank_realm   !! Processor rank from 1 to nproc_realm;
    integer,protected :: &
        rootid_realm !! Root processor in each realm;
    integer,protected :: &
        myid_cart   !! Processor ID from 0 to nproc_cart-1;
    integer,protected :: &
        rank_cart   !! Processor rank from 1 to nproc_cart;
    integer,protected :: &
        rootid_cart !! Root processor in each cart topology;
    integer,protected :: &
        myid_olap   !! Processor ID from 0 to nproc_olap-1;
    integer,protected :: &
        rank_olap   !! Processor rank from 1 to nproc_olap;
    integer,protected :: &
        CFDid_olap  !! Root processor in overlap is the CFD processor;
    integer,protected :: &
        myid_graph  !! Processor ID from 0 to nproc_graph-1;
    integer,protected :: &
        rank_graph  !! Processor rank from 1 to nproc_graph;
    integer,protected :: &
        rank_intersect  !! Processor rank in intersection of overlaping proces;

    !! Get rank in CPL_world_COMM from rank in local COMM
    integer,protected, dimension(:), allocatable    :: &
        rank_world2rank_mdrealm,    &
        rank_world2rank_mdcart,     &
        rank_world2rank_cfdrealm,   &
        rank_world2rank_cfdcart,    &
        rank_world2rank_olap,       &
        rank_world2rank_graph,      &
        rank_world2rank_inter

    !! Get rank in local COMM from rank in CPL_world_COMM
    integer,protected, dimension(:), allocatable    :: &
         rank_mdrealm2rank_world,    &
          rank_mdcart2rank_world,    &
        rank_cfdrealm2rank_world,    &
         rank_cfdcart2rank_world,    &
            rank_olap2rank_world,    &
           rank_graph2rank_world,    &
           rank_inter2rank_world,    &
            rank_olap2rank_realm


    ! Processor topologies
    integer,protected :: &
        nproc_md     !! Total number of processor in md
    integer,protected :: &
        nproc_cfd    !! Total number of processor in cfd
    integer,protected :: &
        nproc_olap   !! Total number of processor in overlap region
    integer,protected :: &
        nproc_world  !! Total number of processor in world
    integer,protected :: &
        npx_md       !! Number of processor in x in the md
    integer,protected :: &
        npy_md       !! Number of processor in y in the md
    integer,protected :: &
        npz_md       !! Number of processor in z in the md
    integer,protected :: &
        npx_cfd      !! Number of processor in x in the cfd
    integer,protected :: &
        npy_cfd      !! Number of processor in y in the cfd
    integer,protected :: &
        npz_cfd      !! Number of processor in z in the cfd

    logical,protected, dimension(:), allocatable :: &
        olap_mask           !! Overlap mask specifying which processors overlap using world ranks
    integer,protected, dimension(:,:), allocatable :: &
        rank2coord_cfd,   &   !! Array containing coordinates for each cartesian rank 
        rank2coord_md
    integer,protected, dimension(:,:,:), allocatable :: &
        coord2rank_cfd,   &
        coord2rank_md

    !! Processor cartesian coords   
    integer,protected :: &
        iblock_realm,     &
        jblock_realm,     &
        kblock_realm

    ! Global cell grid parameters
    integer,protected :: &
        ncx,              &    ! Number of cells in domain
        ncy,              &
        ncz,              &
        icmin,            &    ! Domain cell extents
        icmax,            &
        jcmin,            &
        jcmax,            &
        kcmin,            &
        kcmax,            &
        icmin_olap,       &    ! Overlap region cell extents
        icmax_olap,       &
        jcmin_olap,       &
        jcmax_olap,       &
        kcmin_olap,       &
        kcmax_olap,       &
        ncx_olap,         &    ! Number of cells in overlap region
        ncy_olap,         &
        ncz_olap

    ! Constrained dynamics region flags and params
    integer, protected :: &
        constraint_algo,  &
        constraint_CVflag,&
        icmin_cnst,       &    ! Constrained dynamics region cell extents
        icmax_cnst,       &
        jcmin_cnst,       &
        jcmax_cnst,       &
        kcmin_cnst,       &
        kcmax_cnst

    ! Boundary region 
    integer, protected :: &
        boundary_algo,    &
        icmin_bnry,       &    
        icmax_bnry,       &
        jcmin_bnry,       &
        jcmax_bnry,       &
        kcmin_bnry,       &
        kcmax_bnry

    ! Coupling CFD boundary condition direction flags
    integer, protected :: &
        cpl_cfd_bc_x, &
        cpl_cfd_bc_y, &
        cpl_cfd_bc_z

    !Flag to check if setup has completed successfully
    integer :: CPL_setup_complete = 0
 
    
    ! Coupling constrained regions, average MD quantities 
    ! in spanwise direction (flags)
    integer, protected :: &
        cpl_md_bc_slice, &
        cpl_cfd_bc_slice ! (average MD values, not CFD)

    ! Constraint parameters 
    integer, parameter :: &
        constraint_off = 0,          &
        constraint_OT = 1,           &
        constraint_NCER = 2,         &
        constraint_Flekkoy = 3,      &
        constraint_CV = 4   


    !Sendtype flags
    integer, protected :: &
        sendtype_cfd_to_md,  &
        sendtype_md_to_cfd

    ! Processor cell ranges 
    integer,protected, dimension(:), allocatable :: &
        icPmin_md,        &
        icPmax_md,        &
        jcPmin_md,        &
        jcPmax_md,        &
        kcPmin_md,        &
        kcPmax_md,        &
        icPmin_cfd,       &
        icPmax_cfd,       &
        jcPmin_cfd,       &
        jcPmax_cfd,       &
        kcPmin_cfd,       &
        kcPmax_cfd
    
    ! Domain and cell dimensions
    real(kind(0.d0)),protected :: &
        xL_md,            &
        yL_md,            &
        zL_md,            &
        x_orig_md,        & ! Origin of the md domain
        y_orig_md,        &
        z_orig_md,        &
        xL_cfd,           &
        yL_cfd,           &
        zL_cfd,           &
        x_orig_cfd,       & ! Origin of the cfd domain
        y_orig_cfd,       &
        z_orig_cfd,       &
        xL_olap,          &
        yL_olap,          &
        zL_olap,          &
        xLl,              &
        yLl,              &
        zLl,              &
        dx,               &
        dy,               &
        dz,               &
        dymin,            &
        dymax


    ! Positions of CFD grid lines
    real(kind(0.d0)),protected, dimension(:,:,:), allocatable, target :: &
        xg,               &
        yg,               &
        zg

    ! CFD to MD processor mapping
    integer,protected, dimension(:,:), allocatable :: &
        cfd_icoord2olap_md_icoords, &
        cfd_jcoord2olap_md_jcoords, &
        cfd_kcoord2olap_md_kcoords

    ! Simulation parameters
    integer,protected :: &
        nsteps_md,        & !MD input steps
        nsteps_cfd,       & !CFD input steps
        nsteps_coupled,   & !Total number of steps for coupled simulation
        average_period=1, & ! average period for averages ( it must come from CFD !!!)
        save_period=10      ! save period (corresponts to tplot in CFD, revise please !!!)
    real(kind(0.d0)),protected :: &
        dt_md,            &
        dt_cfd,           &
        density_md,       &
        density_cfd
    integer,protected :: &
        timestep_ratio,        &
        md_cfd_match_cellsize, &
        testval
    logical,protected :: &
        staggered_averages(3) = (/.false.,.false.,.false./), &
        CPL_full_overlap

    ! Communications style
    integer, protected :: &
        comm_style
    integer, parameter :: &
        comm_style_send_recv = 0, &
        comm_style_gath_scat = 1

    !Interface for error handling functions
    interface error_abort
        module procedure error_abort_s, error_abort_si
    end interface error_abort

    private error_abort_si, error_abort_s

contains

!=============================================================================
!                    _____      _               
!                   /  ___|    | |              
!                   \ `--.  ___| |_ _   _ _ __  
!                    `--. \/ _ \ __| | | | '_ \
!                   /\__/ /  __/ |_| |_| | |_) |
!                   \____/ \___|\__|\__,_| .__/ 
!                                        | |    
!                                        |_|    
!=============================================================================






subroutine CPL_init(callingrealm, RETURNED_REALM_COMM, ierror)
! ----------------------------------------------------------------------------
!(cfd+md) Splits MPI_COMM_WORLD in both the CFD and MD code respectively 
!and create intercommunicator between CFD and MD
!
!**Remarks**
!
!Assumes MPI has been initialised `MPI_init` and communicator MPI_COMM_WORLD exists
!and contains all processors in both CFD and MD regions
!
!
!**Synopsis**
!
!.. code-block:: c
!
!  CPL_init(
!           callingrealm, 
!           RETURNED_REALM_COMM, 
!           ierror
!           )    
!
!**Inputs**
!
! - *callingrealm*
!
!   - Should identify calling processor as either CFD_REALM (integer with value 1) or MD_REALM (integer with value 2).
! 
!
!**Outputs**
!
! - RETURNED_REALM_COMM 
!
!   - Communicator based on callingrealm value local to CFD or MD processor and resulting from the split of MPI_COMM_WORLD
!
! - ierror
!
!   - Error flag
!
!**Example**
!
!.. literalinclude:: ../../examples/cpl_init/cfd_init.f90
!
!.. literalinclude:: ../../examples/cpl_init/md_init.f90
!
!
!**Errors**
!
!    COUPLER_ERROR_REALM  = 1       ! wrong realm value
!    COUPLER_ERROR_ONE_REALM = 2    ! one realm missing
!    COUPLER_ERROR_INIT = 3         ! initialisation error
!
! .. sectionauthor::Edward Smith, David Trevelyan, Eduardo Ramos Fernandez
! ------------------------------------
    use mpi
    implicit none

    integer, intent(in)   :: callingrealm ! CFD_REALM=1 or MD_REALM=2
    integer, intent(out)  :: RETURNED_REALM_COMM, ierror

    integer :: MPMD_mode
    logical :: MPI_initialised
    logical, save :: CPL_initialised=.false.

    call MPI_initialized(MPI_initialised, ierr)
    if (.not.MPI_initialised) then
        call error_abort("Error in CPL_init -- MPI not initialised") 
    endif

    if (.not.CPL_initialised) then
        !Set error to zero and copy realm
        ierror=0
        realm = callingrealm

        ! Test if we have an MPMD simulation with 
        ! CFD and a MD realm sharing a single MPI_COMM_WORLD
        ! or if two seperate codes with seperate MPI_COMM_WORLDs 
        ! exist and must be connected by opening a port.
        call test_realms(MPMD_mode)

        ! Create intercommunicator linking realms
        ! and intracommunicators in each of the realms
        call create_comm(MPMD_mode)

        !Return a duplicate to protect the internal CPL_REALM_COMM 
        call MPI_comm_dup(CPL_REALM_COMM, RETURNED_REALM_COMM, ierr)
        CPL_initialised = .true.

        !Print header if necessary
        if (output_mode .ne. QUIET) call print_cplheader
    else
        print*, "CPL_init has been called more than once. Returning same COMM"
        call MPI_comm_dup(CPL_REALM_COMM, RETURNED_REALM_COMM, ierr)    
        !Return an error flag but allow this
        ierror = COUPLER_ERROR_INIT
        !call error_abort("CPL_init is already initialised") 
    endif

contains

!-----------------------------------------------------------------------------
!   Print ascii header for CPL library
!-----------------------------------------------------------------------------

subroutine print_cplheader()
    use mpi
    implicit none

    integer :: ierr

    call MPI_Barrier(CPL_WORLD_COMM, ierr)
    if (myid_world .eq. rootid_world) then

        print*, "                                                                 "
        print*, "  ________/\\\\\\\\\__/\\\\\\\\\\\\\____/\\\_____________        "
        print*, "   _____/\\\////////__\/\\\/////////\\\_\/\\\_____________       "
        print*, "    ___/\\\/___________\/\\\_______\/\\\_\/\\\_____________      "
        print*, "     __/\\\_____________\/\\\\\\\\\\\\\/__\/\\\_____________     "
        print*, "      _\/\\\_____________\/\\\/////////____\/\\\_____________    "
        print*, "       _\//\\\____________\/\\\_____________\/\\\_____________   "
        print*, "        __\///\\\__________\/\\\_____________\/\\\_____________  "
        print*, "         ____\////\\\\\\\\\_\/\\\_____________\/\\\\\\\\\\\\\\\_ "
        print*, "          _______\/////////__\///______________\///////////////__"
        print*, "                                                                 "
        print*, "                     C P L  -  L I B R A R Y                     "
        print*, "                                                                 "

    endif
    call MPI_Barrier(CPL_WORLD_COMM, ierr)

end subroutine print_cplheader

!-----------------------------------------------------------------------------
!   Test if CFD and MD realms are assigned correctly
!-----------------------------------------------------------------------------

subroutine test_realms(MPMD_mode)
    implicit none

    integer,intent(out)  :: MPMD_mode
    integer              :: i, root, rank, nproc, ncfd, nmd
    integer, allocatable :: realm_list(:)


    !Get processor id in world comm
    call MPI_comm_rank(MPI_comm_world, rank, ierr)

    ! Allocate and gather array with realm (MD or CFD) of all
    ! processor in MPI_COMM_WORLD on the root processor
    root = 0
    if (rank .eq. root) then
        call MPI_comm_size(MPI_comm_world, nproc, ierr)
        allocate(realm_list(nproc))
    else
        allocate(realm_list(0))
    endif

    call MPI_gather(callingrealm, 1, MPI_INTEGER, realm_list, &
                    1, MPI_INTEGER, root, MPI_COMM_WORLD, ierr)

    !Check through array of processors on both realms
    !and return error if wrong values.
    if (rank .eq. root) then
        ncfd = 0; nmd = 0
        do i = 1, nproc
            if ( realm_list(i) .eq. cfd_realm ) then 
                ncfd = ncfd + 1
            else if ( realm_list(i) .eq. md_realm ) then
                nmd = nmd +1
            else
                ierror = COUPLER_ERROR_REALM
                print*, "Error in CPL_init --", realm_list(i), "is an unrecognised " // &
                           "callingrealm value Use ", cfd_realm, & 
                           " for CFD or ", md_realm, " for MD."
                call MPI_abort(MPI_COMM_WORLD, ierror, ierr)
            endif
        enddo

        ! Check that both CFD and MD codes exist in MPI_COMM_WORLD
        ! if not, then we need to open a port and link them
        if ( ncfd .eq. 0 ) then
            print*, "Only MD realm present in MPI_COMM_WORLD"
            MPMD_mode = 0
        elseif (nmd .eq. 0) then 
            print*, "Only CFD realm present in MPI_COMM_WORLD"
            MPMD_mode = 0
        else
            print*, "MPMD mode, CFD and MD both share MPI_COMM_WORLD"
            MPMD_mode = 1
        endif

    endif

    call MPI_BCAST(MPMD_mode, 1, MPI_INTEGER, root, MPI_COMM_WORLD, ierr)

end subroutine test_realms

!-----------------------------------------------------------------------------
! Create communicators for each realm and inter-communicator
!-----------------------------------------------------------------------------

subroutine create_comm(MPMD_mode)
    implicit none

    integer, intent(in) :: MPMD_mode

    integer ::  ibuf(2), jbuf(2), remote_leader, comm_size

    character(len=64) :: filename
    character(MPI_MAX_PORT_NAME) :: port
    integer :: port_connect, unitno, rsize
    logical :: portfileexists

    if (MPMD_mode .eq. 1) then
        ! Split MPI COMM WORLD ready to establish two communicators
        ! 1) A global intra-communicator in each realm for communication
        ! internally between CFD processes or between MD processes
        ! 2) An inter-communicator which allows communication between  
        ! the 'groups' of processors in MD and the group in the CFD
        call MPI_comm_dup(MPI_COMM_WORLD, CPL_WORLD_COMM, ierr)

        !Get processor id in world across both realms
        call MPI_comm_rank(CPL_WORLD_COMM, myid_world, ierr)
        rank_world = myid_world + 1; rootid_world = 0
        call MPI_comm_size(CPL_WORLD_COMM, nproc_world, ierr)

        !------------ create realm intra-communicators -----------------------
        ! Split CPL_WORLD_COMM into an intra-communicator for each realm 
        ! (used for any communication within each realm - e.g. broadcast from 
        !  an md process to all other md processes)
        call MPI_comm_split(CPL_WORLD_COMM, callingrealm, & 
                            myid_world, CPL_REALM_COMM, ierr)

        !Define processor ID
        call MPI_comm_rank(CPL_REALM_COMM, myid_realm, ierr)
        rank_realm = myid_realm + 1; rootid_realm = 0

        !------------ create realm inter-communicators -----------------------
        ! Create intercommunicator between the group of processor on each realm
        ! (used for any communication between realms - e.g. md group rank 2 sends
        ! to cfd group rank 5). inter-communication is by a single processor on 
        ! each group.

        ! Get the MPI_comm_world ranks that hold the largest ranks in cfd_comm and md_comm
        call MPI_comm_size(CPL_REALM_COMM, comm_size, ierr)
        ibuf(:) = -1
        jbuf(:) = -1
        if ( myid_realm .eq. comm_size - 1) then
            ibuf(realm) = myid_world
        endif

        call MPI_allreduce( ibuf ,jbuf, 2, MPI_INTEGER, MPI_MAX, &
                            CPL_WORLD_COMM, ierr)

        !Set this largest rank on each process to be the inter-communicators (WHY NOT 0??)
        select case (callingrealm)
        case (cfd_realm)
            remote_leader = jbuf(md_realm)
        case (md_realm)
            remote_leader = jbuf(cfd_realm)
        end select

        call MPI_intercomm_create(CPL_REALM_COMM, comm_size - 1, CPL_WORLD_COMM, &
                                  remote_leader, 1, CPL_INTER_COMM, ierr)

    elseif (MPMD_mode .eq. 0) then
        ! Here we need to open a port and wait for connection from 
        ! the other code which has been started as a seperate MPI
        ! instance with its own MPI_COMM_WORLD and exchange port
        ! information by writing to filename
        filename = "./port"

        !MPI_COMM_WORLD is the realm comm now
        call MPI_Comm_dup(MPI_COMM_WORLD, CPL_REALM_COMM, ierr)
        call MPI_Comm_size(CPL_REALM_COMM, comm_size, ierr)
        call MPI_Comm_rank(CPL_REALM_COMM, myid_realm, ierr)
        rank_realm = myid_realm + 1; rootid_realm = 0

        !Two seperate programs specified by the realm input argument
        if (callingrealm .eq. cfd_realm) then

            !Just root file needs to open a port
            if (myid_realm .eq. rootid_realm) then
                call MPI_Open_port(MPI_INFO_NULL, port, ierr)
                print*, "opened port: ", trim(port), & 
                        " Attempting to write to file:", trim(filename)

                !Write port file
                ierr = -1
                do while(ierr .ne. 0)
                    unitno = get_new_fileunit()
                    open(unitno, file=trim(filename), & 
                         action="write", iostat=ierr)
                    if (ierr .eq. 0) then
                        write(unitno,*) port
                        close(unitno)
                        print*, "Portname written to file ", filename
                    elseif (ierr .eq. 23) then
                        call sleep(3)
                    else
                        print*, "Error ", ierr , & 
                                "when attempting to write Port file: ", &
                                trim(filename)
                    endif
                enddo
            endif

            !All processors then attempt to establish connection
            call MPI_Comm_accept(port, MPI_INFO_NULL, 0, CPL_REALM_COMM, CPL_INTER_COMM, ierr)
            call MPI_Comm_remote_size(CPL_INTER_COMM, rsize, ierr)
            print*, "accepted connection on port to root of ", rsize, " procs."

        elseif (callingrealm .eq. md_realm) then
            ! Keep trying until successful connection to the port
            call MPI_Errhandler_set(CPL_REALM_COMM, MPI_ERRORS_RETURN, ierr)
            port_connect = MPI_ERR_PORT
            do while (port_connect .ne. MPI_SUCCESS)

                !Wait until a port file exists and is avaialble
                portfileexists = .false.
                do while(portfileexists .eqv. .false.)
                    inquire(file=trim(filename), exist=portfileexists)
                    if (portfileexists) then
                        unitno = get_new_fileunit()
                        open(unitno, file=trim(filename), & 
                             action="read", iostat=ierr)
                        if (ierr .eq. 0) then
                            read(unitno, *) port  
                        !Catch file already open error                       
                        elseif (ierr .eq. 23) then
                            call sleep(3)
                            cycle 
                        else
                            print*, "Error ", ierr , & 
                                   " when attempting to read Port file:", & 
                                   trim(filename)
                        endif
                    else
                        print*, "Cannot find port file: ", trim(filename), ". Waiting..."
                        call sleep(3)
                    endif
                enddo

                !Attempt to establish connection 
                call MPI_Comm_connect(port, MPI_INFO_NULL, 0, CPL_REALM_COMM, CPL_INTER_COMM, port_connect)
                if (port_connect .ne. MPI_SUCCESS) then
                    print*, "Error -- failed to connected to port ",  port 
                    call sleep(3)
                endif
            enddo

            !Delete file
            if (myid_realm .eq. rootid_realm) then
                close(unitno, status="delete")
            endif

            call MPI_Comm_remote_size(CPL_INTER_COMM, rsize, ierr)
            print*, "connection accepted to root of ", rsize , " procs." 
            call MPI_Errhandler_set(CPL_REALM_COMM, MPI_ERRORS_ARE_FATAL, ierr)
        else
            stop "Error -- Realm should be 1 or 2 "
        endif

        call MPI_Barrier(CPL_INTER_COMM, ierr)
        call MPI_Intercomm_merge(CPL_INTER_COMM, .true., CPL_WORLD_COMM, ierr)
        !Get processor id in world across both realms
        call MPI_comm_rank(CPL_WORLD_COMM, myid_world, ierr)
        rank_world = myid_world + 1; rootid_world = 0
        call MPI_comm_size(CPL_WORLD_COMM, nproc_world, ierr)
        print*, "Rank on realm ", realm, " is ", rank_realm, " of ", comm_size, & 
                "  and rank on intercomm is ", myid_world, " of ", nproc_world


    else

        call error_abort("Error in CPL_init -- MPMD_mode should be 0 or 1") 

    endif
    if (output_mode .ne. QUIET) then
        print*, 'Completed CPL communicator init for ', realm_name(realm), &
                ' , CPL_WORLD_COMM ID:', myid_world
    endif

end subroutine create_comm

end subroutine CPL_init


subroutine CPL_finalize(ierr)
    use mpi, only : MPI_COMM_NULL, MPI_Barrier, MPI_COMM_WORLD
    implicit none

    integer, intent(out) :: ierr

    !Comminicators setup by CPL_init()
    if (CPL_INTER_COMM .ne. MPI_COMM_NULL) call MPI_COMM_FREE(CPL_INTER_COMM, ierr)
    if (CPL_REALM_COMM .ne. MPI_COMM_NULL) call MPI_COMM_FREE(CPL_REALM_COMM, ierr)
    if (CPL_WORLD_COMM .ne. MPI_COMM_NULL) call MPI_COMM_FREE(CPL_WORLD_COMM, ierr)

    !Free communicators setup by CPL_setup
    if (CPL_setup_complete .eq. 1) then
        if (CPL_OLAP_COMM .ne. MPI_COMM_NULL) call MPI_COMM_FREE(CPL_OLAP_COMM, ierr)
        if (CPL_CART_COMM .ne. MPI_COMM_NULL) call MPI_COMM_FREE(CPL_CART_COMM, ierr)
        if (CPL_GRAPH_COMM .ne. MPI_COMM_NULL) call MPI_COMM_FREE(CPL_GRAPH_COMM, ierr)
    endif

    !Barrier over both CFD and MD realms
    call MPI_Barrier(MPI_COMM_WORLD, ierr)

end subroutine CPL_finalize

!=============================================================================
!! Read Coupler input file
!-----------------------------------------------------------------------------

subroutine read_coupler_input()
    implicit none

    integer :: infileid, readin
    logical :: found

    !Check all file ids until an unused one is found
    infileid = 100000
    do 
        inquire(unit=infileid,opened=found)
        if (.not.(found)) exit
        infileid = infileid + 1
    enddo

    !Open and read input file on all processes
    open(infileid,file='cpl/COUPLER.in',status="old",action="read", &
                  form="formatted")

    !Possible to specify full overlap so no need for other details
    call locate(infileid, 'FULL_OVERLAP', found)
    if (found) then
        read(infileid,*, IOSTAT=readin) CPL_full_overlap
    else
        CPL_full_overlap = .false.       
    endif

    if (CPL_full_overlap .eqv. .false.) then

        call locate(infileid, 'OVERLAP_EXTENTS', found)
        if (found) then
            read(infileid,*) icmin_olap
            read(infileid,*) icmax_olap
            read(infileid,*) jcmin_olap
            read(infileid,*) jcmax_olap
            read(infileid,*) kcmin_olap
            read(infileid,*) kcmax_olap
        else
            call error_abort("Ovelap extents unspecified in coupler input file.")
        end if

        call locate(infileid, 'CONSTRAINT_INFO', found)
        if (found) then
            read(infileid,*) constraint_algo
            if (constraint_algo .ne. 0) then
                read(infileid,*) constraint_CVflag
                read(infileid,*) icmin_cnst
                read(infileid,*) icmax_cnst
                read(infileid,*) jcmin_cnst
                read(infileid,*) jcmax_cnst
                read(infileid,*) kcmin_cnst
                read(infileid,*) kcmax_cnst
            endif
        else
            call error_abort("CONSTRAINT_INFO not specified in coupler input file")
        end if

        call locate(infileid, 'BOUNDARY_EXTENTS', found)
        if (found) then
            boundary_algo = 1
            read(infileid,*) icmin_bnry
            read(infileid,*) icmax_bnry
            read(infileid,*) jcmin_bnry
            read(infileid,*) jcmax_bnry
            read(infileid,*) kcmin_bnry
            read(infileid,*) kcmax_bnry
        else
            boundary_algo = 0
            icmin_bnry = VOID
            icmax_bnry = VOID
            jcmin_bnry = VOID
            jcmax_bnry = VOID
            kcmin_bnry = VOID
            kcmax_bnry = VOID
        end if

    else
        call locate(infileid, 'CONSTRAINT_INFO', found)
        if (found) then
            read(infileid,*) constraint_algo
        endif
    endif


    call locate(infileid, 'DENSITY_CFD', found)
    if (found) then
        read(infileid,*) density_cfd
    else
        call error_abort("Density not specified in coupler input file.")
    end if
    
    call locate(infileid, 'TIMESTEP_RATIO', found)
    if (found) then
        read(infileid,*) timestep_ratio !TODO name change
    else
        timestep_ratio = VOID
    end if
    
    call locate(infileid, 'MATCH_CELLSIZE', found)
    if (found) then
        read(infileid,*) md_cfd_match_cellsize
    else
        md_cfd_match_cellsize = 0
    end if

    call locate(infileid, 'CPL_CFD_BC_XYZ', found)
    if (found) then
        read(infileid,*) cpl_cfd_bc_x 
        read(infileid,*) cpl_cfd_bc_y 
        read(infileid,*) cpl_cfd_bc_z 
    else
        cpl_cfd_bc_x = 1
        cpl_cfd_bc_y = 0
        cpl_cfd_bc_z = 1
    end if

    call locate(infileid, 'CPL_CFD_BC_SLICE', found)
    if (found) then
        read(infileid,*) cpl_cfd_bc_slice
    else
        cpl_cfd_bc_slice = 0 
    end if

    call locate(infileid, 'CPL_MD_BC_SLICE', found)
    if (found) then
        read(infileid,*) cpl_md_bc_slice
    else
        cpl_md_bc_slice = 0 
    end if

    call locate(infileid, 'COMM_STYLE', found) 
    if (found) then
        read(infileid, *) comm_style
    else
        comm_style = comm_style_send_recv
    end if

    call locate(infileid, 'SENDTYPE_MD_TO_CFD', found) 
    if (found) then
        read(infileid, *) sendtype_md_to_cfd
    else
        sendtype_md_to_cfd = 1
    end if

    call locate(infileid, 'SENDTYPE_CFD_TO_MD', found) 
    if (found) then
        read(infileid, *) sendtype_cfd_to_md
    else
        sendtype_cfd_to_md = 1
    end if

    close(infileid,status="keep")

    if (myid_world .eq. rootid_world) then
        call CPL_write_header("./cpl/coupler_header")!todo better name/place
    end if

end subroutine read_coupler_input

!------------------------------------------------------------------------------
!                              CPL_write_header                               -
!------------------------------------------------------------------------------
!>
!! Writes header information to specified filename in the format
!! Variable description ; variable name ; variable
!!
!! - Synopsis
!!
!!  - CPL_write_header (header_filename)
!!
!! - Input
!!
!!  - header_filename
!!   - File name to write header to
!!
!! - Input/Output
!!  - NONE
!!
!! - Output
!!  - NONE
!! 
!! .. sectionauthor:: Edward Smith
!
! ----------------------------------------------------------------------------

subroutine CPL_write_header(header_filename)
    implicit none

    character(*),intent(in) :: header_filename

    logical             :: found
    integer             :: infileid
    Character(8)        :: the_date
    Character(10)       :: the_time

    !Check all file ids until an unused one is found
    infileid = 100000
    do 
        inquire(unit=infileid,opened=found)
        if (.not.(found)) exit
        infileid = infileid + 1
    enddo

    ! Write Simulation Parameter File contain all data required to completely recreate
    ! simulation and to be used for post processing
    call date_and_time(the_date, the_time)

    !Open and write input file on all processes
    !if (rank_world .eq. rootid_world+1) then

        open(infileid,file=trim(header_filename),action="write",form="formatted")

        write(infileid,*) 'Simulation run on Date;  sim_date ;', the_date
        write(infileid,*) 'Simulation start time ;  sim_start_time ;', the_time
        write(infileid,*) 'CFD density;  density_cfd ;', density_cfd
        write(infileid,*) 'choice of constraint algorithm;  constraint_algo ;', constraint_algo
        if (constraint_algo .ne. 0) then
            write(infileid,*) 'CV form of constraint;  constraint_CVflag ;', constraint_CVflag
            write(infileid,*) 'minimum x cell of constrained region;  icmin_cnst ;', icmin_cnst
            write(infileid,*) 'maximum x cell of constrained region;  icmax_cnst ;', icmax_cnst
            write(infileid,*) 'minimum y cell of constrained region;  jcmin_cnst ;', jcmin_cnst
            write(infileid,*) 'maximum y cell of constrained region;  jcmax_cnst ;', jcmax_cnst
            write(infileid,*) 'minimum z cell of constrained region;  kcmin_cnst ;', kcmin_cnst
            write(infileid,*) 'maximum z cell of constrained region;  kcmax_cnst ;', kcmax_cnst
        endif

        if (boundary_algo .ne. 0) then
            write(infileid,*) 'Boundary constrain specified;  boundary_algo ;', boundary_algo
            write(infileid,*) 'minimum x cell of boundary region;  icmin_bnry ;', icmin_bnry
            write(infileid,*) 'maximum x cell of boundary region;  icmax_bnry ;', icmax_bnry
            write(infileid,*) 'minimum y cell of boundary region;  jcmin_bnry ;', jcmin_bnry
            write(infileid,*) 'maximum y cell of boundary region;  jcmax_bnry ;', jcmax_bnry
            write(infileid,*) 'minimum z cell of boundary region;  kcmin_bnry ;', kcmin_bnry
            write(infileid,*) 'maximum z cell of boundary region;  kcmax_bnry ;', kcmax_bnry
        endif

        write(infileid,*) 'minimum x cell of overlap region;  icmin_olap ;', icmin_olap
        write(infileid,*) 'maximum x cell of overlap region;  icmax_olap ;', icmax_olap
        write(infileid,*) 'minimum y cell of overlap region;  jcmin_olap ;', jcmin_olap
        write(infileid,*) 'maximum y cell of overlap region;  jcmax_olap ;', jcmax_olap
        write(infileid,*) 'minimum z cell of overlap region;  kcmin_olap ;', kcmin_olap
        write(infileid,*) 'maximum z cell of overlap region;  kcmax_olap ;', kcmax_olap

        write(infileid,*) 'MD timesteps per CFD timestep;  timestep_ratio ;', timestep_ratio !TODO name change
        write(infileid,*) 'Enforce cellsize matching;  md_cfd_match_cellsize ;', md_cfd_match_cellsize

        close(infileid,status='keep')

    !endif

end subroutine  CPL_write_header


subroutine CPL_setup_cfd(icomm_grid, xyzL, xyz_orig, ncxyz)
! ----------------------------------------------------------------------------
!Initialisation routine for coupler module - Every variable is sent and stored
!to ensure both md and cfd region have an identical list of parameters
!
!**Remarks**
!
!Assumes CPL has been initialised `CPL_init` and communicator MD_REALM exists
!
!**Synopsis**
!
!.. code-block:: c
!
!  coupler_cfd_init(
!                  icomm_grid,
!                  xyzL,
!                  xyz_orig,
!                  ncxyz,
!                  )
!
!**Inputs**
!
! - *icomm_grid*
!
!   - Communicator based on CFD processor topology returned from a call to MPI_CART_CREATE.
! - *xyzL*
!
!   - CFD domain size.
! - *xyz_orig*
!
!   - CFD origin.
! - *ncxyz*
!
!   - Number of CFD cells in global domain.
!
! .. sectionauthor::Edward Smith, David Trevelyan, Eduardo Ramos Fernandez
! ------------------------------------
    use mpi
    implicit none           
    
    ! Params
    integer, intent(in)                           :: icomm_grid 
    integer, dimension(3), intent(in)             :: ncxyz
    real(kind(0.d0)), dimension(3),  intent(in)   :: xyzL, xyz_orig

    ! Vars
    integer, dimension(:), allocatable            :: iTmin, iTmax, jTmin,&
                                                     jTmax, kTmin, kTmax
    integer, dimension(:,:), allocatable          :: icoord
    real(kind(0.d0)), dimension(:,:,:), allocatable :: xgrid, ygrid, zgrid
    integer, dimension(3)                         :: ijkcmin, ijkcmax
    integer, dimension(3)                         :: npxyz_cfd, cart_coords
    logical, dimension(3)                         :: cart_periods
    real(kind(0.d0))                              :: dx, dy, dz
    integer                                       :: cart_nprocs, rank, x, y, z,&
                                                     ncxl, ncyl, nczl, i, j, k

    ! Ranges of cells 
    ijkcmin = (/1,1,1/)
    ijkcmax = ncxyz

    ! Get number of processors in each direction (periods and coords are not needed)
    call MPI_Cart_get(icomm_grid, 3, npxyz_cfd, cart_periods, cart_coords, ierr)
    call MPI_Comm_size(icomm_grid, cart_nprocs, ierr)
    
    allocate(icoord(3, cart_nprocs), stat=ierr)
    allocate(iTmin(npxyz_cfd(1)), stat=ierr)
    allocate(iTmax(npxyz_cfd(1)), stat=ierr)
    allocate(jTmin(npxyz_cfd(2)), stat=ierr)
    allocate(jTmax(npxyz_cfd(2)), stat=ierr)
    allocate(kTmin(npxyz_cfd(3)), stat=ierr)
    allocate(kTmax(npxyz_cfd(3)), stat=ierr)

    ! I am not sure if this should be checked for error here 
    ncxl = ncxyz(1) / npxyz_cfd(1)
    ncyl = ncxyz(2) / npxyz_cfd(2)
    nczl = ncxyz(3) / npxyz_cfd(3)
    
    ! Do the mapping from rank to cartesian coords
    do rank=0, cart_nprocs - 1
        call MPI_Cart_coords(icomm_grid, rank, 3, cart_coords, ierr)
        icoord(1:3, rank + 1) = cart_coords + 1
        x = cart_coords(1)
        y = cart_coords(2)
        z = cart_coords(3)
        iTmin(x+1) = x*ncxl + 1
        iTmax(x+1) = iTmin(x+1) + ncxl - 1
        jTmin(y+1) = y*ncyl + 1
        jTmax(y+1) = jTmin(y+1) + ncyl - 1
        kTmin(z+1) = z*nczl + 1
        kTmax(z+1) = kTmin(z+1) + nczl - 1
    enddo

    !Check if allocating global meshgrid is reasonable
    if (product(ncxyz) .lt. maxgridsize) then
        allocate(xgrid(ncxyz(1) + 1, ncxyz(2) + 1, ncxyz(3) + 1), stat=ierr)
        allocate(ygrid(ncxyz(1) + 1, ncxyz(2) + 1, ncxyz(3) + 1), stat=ierr)
        allocate(zgrid(ncxyz(1) + 1, ncxyz(2) + 1, ncxyz(3) + 1), stat=ierr)
        ! Construct cartesian grid
        dx =  xyzL(1) / ncxyz(1)
        dy =  xyzL(2) / ncxyz(2)
        dz =  xyzL(3) / ncxyz(3)
        do i=1, ncxyz(1) + 1
        do j=1, ncxyz(2) + 1
        do k=1, ncxyz(3) + 1
                xgrid(i, j, k) = (i-1) * dx + xyz_orig(1)
                ygrid(i, j, k) = (j-1) * dy + xyz_orig(2)
                zgrid(i, j, k) = (k-1) * dz + xyz_orig(3)
        enddo
        enddo
        enddo
    else
        !Other allocate 3 arrays for each dimension
        allocate(xgrid(ncxyz(1) + 1, 1, 1), stat=ierr)
        allocate(ygrid(1, ncxyz(2) + 1, 1), stat=ierr)
        allocate(zgrid(1, 1, ncxyz(3) + 1), stat=ierr)
        ! Construct cartesian grid
        dx =  xyzL(1) / ncxyz(1)
        dy =  xyzL(2) / ncxyz(2)
        dz =  xyzL(3) / ncxyz(3)
        do i=1, ncxyz(1) + 1
                xgrid(i, 1, 1) = (i-1) * dx + xyz_orig(1)
        enddo
        do j=1, ncxyz(2) + 1
                ygrid(1, j, 1) = (j-1) * dy + xyz_orig(2)
        enddo
        do k=1, ncxyz(3) + 1
                zgrid(1, 1, k) = (k-1) * dz + xyz_orig(3)
        enddo
    endif


    call coupler_cfd_init(icomm_grid, icoord, npxyz_cfd, xyzL, &
                          xyz_orig, ncxyz, ijkcmax, ijkcmin, iTmin, &
                          iTmax, jTmin, jTmax, kTmin, kTmax, & 
                          xgrid, ygrid, zgrid)
    deallocate(icoord)
    deallocate(iTmin)
    deallocate(iTmax)
    deallocate(jTmin)
    deallocate(jTmax)
    deallocate(kTmin)
    deallocate(kTmax)
    deallocate(xgrid)
    deallocate(ygrid)
    deallocate(zgrid)

    !Set flag to register setup is complete correctly
    CPL_setup_complete = 1
    call MPI_Barrier(MPI_COMM_WORLD, ierr)
 
end subroutine CPL_setup_cfd



!------------------------------------------------------------------------------
!                              coupler_cfd_init                               -
!------------------------------------------------------------------------------
!>
!! Initialisation routine for coupler module - Every variable is sent and stored
!! to ensure both md and cfd region have an identical list of parameters
!!
!! - Synopsis
!!
!!  - coupler_cfd_init(icomm_grid,icoord,npxyz_cfd,xyzL,ncxyz,
!!                     ijkcmax,ijkcmin,iTmin,iTmax,jTmin,
!!                     jTmax,kTmin,kTmax,xg,yg,zg)
!!
!! - Input
!!
!!  - icomm_grid
!!   - The MPI communicator setup by the MPI_CART_CREATE command in the 
!!     CFD region (integer)
!!  - icoord
!!   - The three coordinate for each rank in the domain (integer array nproc by 3)
!!  - npxyz_cfd
!!   - Number of processors in each cartesian dimension (integer array 3)
!!  - xyzL
!!   - Size of domain in each cartesian dimension (dp real array 3)
!!  - ncxyz
!!   - Global number of cells in each cartesian dimension (integer array 3)
!!  - ijkcmax
!!   - Global maximum cell in each cartesian dimension (integer array 3)
!!  - ijkcmin
!!   - Global minimum cell in each cartesian dimension (integer array 3)
!!  - iTmin
!!   - Local minimum cell for each rank (integer array no. procs in x)
!!  - iTmax
!!   - Local maximum cell for each rank (integer array no. procs in x)
!!  - jTmin
!!   - Local minimum cell for each rank (integer array no. procs in y)
!!  - jTmax
!!   - Local maximum cell for each rank (integer array no. procs in y)
!!  - kTmin
!!   - Local minimum cell for each rank (integer array no. procs in z)
!!  - kTmax
!!   - Local maximum cell for each rank (integer array no. procs in z)
!!  - xg
!!   - Array of cell vertices in the x direction 
!!     (no. cells in x + 1 by no. cells in y + 1, no. cells in z+1)
!!  - yg
!!   - Array of cell vertices in the y direction 
!!     (no. cells in x + 1 by no. cells in y + 1, no. cells in z+1)
!!  - zg
!!   - Array of cell vertices in the z direction 
!!     (no. cells in x + 1 by no. cells in y + 1, no. cells in z+1)
!!
!! - Input/Output
!!  - NONE
!!
!! - Output
!!  - NONE
!! 
!! .. sectionauthor:: Edward Smith
!
! ----------------------------------------------------------------------------


subroutine coupler_cfd_init(icomm_grid, icoord, npxyz_cfd, xyzL, xyz_orig, ncxyz, &
                            ijkcmax, ijkcmin, iTmin, iTmax, jTmin, & 
                            jTmax, kTmin, kTmax, xgrid, ygrid, zgrid)
    use mpi
    implicit none           

    integer,                        intent(in)    :: icomm_grid 
    integer,dimension(3),           intent(in)    :: ijkcmin,ijkcmax,npxyz_cfd,ncxyz
    integer,dimension(:),           intent(in)    :: iTmin,iTmax,jTmin,jTmax,kTmin,kTmax
    integer,dimension(:,:),         intent(in)    :: icoord
    real(kind(0.d0)),dimension(3),  intent(in)    :: xyzL, xyz_orig
    real(kind(0.d0)),dimension(:,:,:),intent(in)  :: xgrid,ygrid,zgrid


    integer                                         :: i,ib,jb,kb,pcoords(3),source,nproc
    integer,dimension(:),allocatable                :: buf,rank_world2rank_realm,rank_world2rank_cart
    real(kind=kind(0.d0))                           :: dxmin,dxmax,dzmin,dzmax
    real(kind=kind(0.d0)),dimension(:),allocatable  :: rbuf

    ! Read COUPLER.in input file
    call read_coupler_input()

    ! Duplicate grid communicator for coupler use
    call MPI_comm_dup(icomm_grid, CPL_CART_COMM, ierr)
    call MPI_comm_rank(CPL_CART_COMM, myid_cart, ierr) 
    rank_cart = myid_cart + 1; rootid_cart = 0
    !Send only from root processor
    if (myid_realm .eq. rootid_realm ) then
        source=MPI_ROOT
    else
        source=MPI_PROC_NULL
    endif

        ! ================ Exchange and store Data ==============================
    ! Data is stored to the coupler module with the same name in both realms
    ! Note - MPI Broadcast between intercommunicators is only supported by MPI-2

    ! ------------------------ Processor Topology ---------------------------
    ! Store & Send CFD number of processors
    npx_cfd = npxyz_cfd(1)
    npy_cfd = npxyz_cfd(2)
    npz_cfd = npxyz_cfd(3)
    nproc_cfd = npx_cfd * npy_cfd * npz_cfd
    call MPI_bcast(npxyz_cfd,3,MPI_INTEGER,source,CPL_INTER_COMM,ierr)  !Send

    ! Receive & Store MD number of processors
    allocate(buf(3))
    call MPI_bcast(buf, 3 , MPI_INTEGER, 0, CPL_INTER_COMM, ierr)  !Receive
    npx_md = buf(1)
    npy_md = buf(2)
    npz_md = buf(3)
    nproc_md = npx_md * npy_md * npz_md
    deallocate(buf)

    ! Store & Send CFD processor rank to coord
    allocate(rank2coord_cfd(3,nproc_cfd),stat=ierr); rank2coord_cfd = icoord
    
    iblock_realm=icoord(1,rank_realm) 
    jblock_realm=icoord(2,rank_realm)
    kblock_realm=icoord(3,rank_realm)
    allocate(buf(3*nproc_cfd)); buf = reshape(icoord, (/ 3*nproc_cfd /) )
    call MPI_bcast(buf,3*nproc_cfd,MPI_INTEGER,source,CPL_INTER_COMM,ierr)  !Send
    deallocate(buf)

    ! Receive & Store MD processor rank to coord
    allocate(buf(3*nproc_md))
    call MPI_bcast(buf,3*nproc_md ,MPI_INTEGER,  0   ,CPL_INTER_COMM,ierr)  !Receive
    allocate(rank2coord_md (3,nproc_md),stat=ierr); rank2coord_md = reshape(buf,(/ 3,nproc_md /))
    deallocate(buf)


    ! Setup CFD mapping from coordinate to rank
    ! Store & Send CFD mapping from coordinate to rank to MD
    allocate(coord2rank_cfd(npx_cfd,npy_cfd,npz_cfd))
    do ib = 1,npx_cfd
    do jb = 1,npy_cfd
    do kb = 1,npz_cfd
        pcoords = (/ ib, jb, kb /)-1
        call MPI_Cart_rank(CPL_CART_COMM,pcoords,i,ierr)
        coord2rank_cfd(ib,jb,kb) = i + 1
    enddo
    enddo
    enddo

    allocate(buf(nproc_cfd)); buf = reshape(coord2rank_cfd, (/ nproc_cfd /) )
    call MPI_bcast(coord2rank_cfd,nproc_cfd,MPI_INTEGER,source,CPL_INTER_COMM,ierr) !Send
    deallocate(buf)

    ! Receive & Store MD coordinate to rank mapping
    allocate(buf(nproc_md))
    call MPI_bcast(buf,nproc_md ,MPI_INTEGER,  0   ,CPL_INTER_COMM,ierr)    !Receive
    allocate(coord2rank_md (npx_md,npy_md,npz_md)) 
    coord2rank_md = reshape(buf,(/ npx_md,npy_md,npz_md /))
    deallocate(buf)

    ! Setup CFD mapping between realm & world rank 
    allocate(rank_cfdrealm2rank_world(nproc_cfd))
    allocate(rank_world2rank_realm(nproc_world))
    call CPL_rank_map(CPL_REALM_COMM,rank_realm,nproc, & 
                        rank_cfdrealm2rank_world,rank_world2rank_realm,ierr)

    !World to rank is the same on both realms
    allocate(rank_world2rank_cfdrealm(nproc_world))
    allocate(rank_world2rank_mdrealm(nproc_world))
    rank_world2rank_cfdrealm = rank_world2rank_realm
    rank_world2rank_mdrealm  = rank_world2rank_realm

    ! Send CFD mapping to MD
    call MPI_bcast(rank_cfdrealm2rank_world,nproc_cfd,MPI_integer,source,CPL_INTER_COMM,ierr)    !send

    ! Receive & Store MD mapping from realm to local rank from MD
    allocate(rank_mdrealm2rank_world(nproc_md))
    call MPI_bcast(rank_mdrealm2rank_world,nproc_md,MPI_integer,0,CPL_INTER_COMM,ierr)  !Receive

    ! Setup CFD mapping between cartesian topology & world rank 
    allocate(rank_cfdcart2rank_world(nproc_cfd))
    allocate(rank_world2rank_cart(nproc_world))
    call CPL_rank_map(CPL_CART_COMM,rank_cart,nproc, & 
                        rank_cfdcart2rank_world,rank_world2rank_cart,ierr)

    !World to rank is the same on both realms cart
    allocate(rank_world2rank_cfdcart(nproc_world))
    allocate(rank_world2rank_mdcart(nproc_world))
    rank_world2rank_cfdcart = rank_world2rank_cart
    rank_world2rank_mdcart  = rank_world2rank_cart

    ! Send CFD mapping to MD
    call MPI_bcast(rank_cfdcart2rank_world,nproc_cfd,MPI_integer,source,CPL_INTER_COMM,ierr)     !send

    ! Receive & Store MD mapping from cart to local rank from MD
    allocate(rank_mdcart2rank_world(nproc_md))
    call MPI_bcast(rank_mdcart2rank_world,nproc_md,MPI_integer,0,CPL_INTER_COMM,ierr)   !Receive

    ! ------------------ Send CFD grid extents ------------------------------

    ! Store & send CFD domain size
    xL_cfd = xyzL(1); yL_cfd = xyzL(2); zL_cfd = xyzL(3)
    call MPI_bcast(xyzL,3,MPI_double_precision,source,CPL_INTER_COMM,ierr)  !Send

    ! Store & send CFD domain origin coords
    x_orig_cfd = xyz_orig(1); y_orig_cfd = xyz_orig(2); z_orig_cfd = xyz_orig(3)
    call MPI_bcast(xyz_orig,3,MPI_double_precision,source,CPL_INTER_COMM,ierr)  !Send

    ! Receive & store MD domain size
    allocate(rbuf(3))
    call MPI_bcast(rbuf,3,MPI_double_precision,0,CPL_INTER_COMM,ierr)       !Receive
    xL_md = rbuf(1); yL_md = rbuf(2); zL_md = rbuf(3);
    deallocate(rbuf)

    ! Receive & store MD domain origin coords
    allocate(rbuf(3))
    call MPI_bcast(rbuf,3,MPI_double_precision,0,CPL_INTER_COMM,ierr)       !Receive
    x_orig_md = rbuf(1); y_orig_md = rbuf(2); z_orig_md = rbuf(3);
    deallocate(rbuf)

    ! Store & send CFD grid extents
    icmin = ijkcmin(1); jcmin = ijkcmin(2); kcmin = ijkcmin(3)
    icmax = ijkcmax(1); jcmax = ijkcmax(2); kcmax = ijkcmax(3)
    allocate(buf(6))
    buf = (/ icmin, icmax, jcmin, jcmax, kcmin, kcmax /)
    call MPI_bcast(buf,6,MPI_INTEGER,source,CPL_INTER_COMM,ierr) !Send
    deallocate(buf)

    !Set overlap, constraint and boundary to no. CFD cells if full overlap 
    if (CPL_full_overlap) then
        icmin_olap = icmin; icmax_olap = icmax
        jcmin_olap = jcmin; jcmax_olap = jcmax
        kcmin_olap = kcmin; kcmax_olap = kcmax
        
        icmin_cnst = icmin; icmax_cnst = icmax
        jcmin_cnst = jcmin; jcmax_cnst = jcmax 
        kcmin_cnst = kcmin; kcmax_cnst = kcmax

        icmin_bnry = icmin; icmax_bnry = icmax
        jcmin_bnry = jcmin; jcmax_bnry = jcmax 
        kcmin_bnry = kcmin; kcmax_bnry = kcmax
    endif

    ! Store & send global number of cells in CFD
    ncx = ncxyz(1); ncy = ncxyz(2); ncz = ncxyz(3)
    call MPI_bcast(ncxyz,3,MPI_INTEGER,source,CPL_INTER_COMM,ierr)              !Send

    ! Store & send array of global grid points
    allocate(xg(size(xgrid,1)+1,size(xgrid,2)+1,size(xgrid,3)+1),stat=ierr); xg = xgrid
    allocate(yg(size(xgrid,1)+1,size(xgrid,2)+1,size(xgrid,3)+1),stat=ierr); yg = ygrid
    allocate(zg(size(xgrid,1)+1,size(xgrid,2)+1,size(xgrid,3)+1),stat=ierr); zg = zgrid

    call MPI_bcast(xgrid,size(xgrid),MPI_double_precision,source,CPL_INTER_COMM,ierr) !Send
    call MPI_bcast(ygrid,size(ygrid),MPI_double_precision,source,CPL_INTER_COMM,ierr) !Send
    call MPI_bcast(zgrid,size(zgrid),MPI_double_precision,source,CPL_INTER_COMM,ierr) !Send

    ! Store & Send local (processor) CFD grid extents
    allocate(icPmin_cfd(npx_cfd),stat=ierr); icPmin_cfd(:) = iTmin(:)
    allocate(icPmax_cfd(npx_cfd),stat=ierr); icPmax_cfd(:) = iTmax(:)
    allocate(jcPmin_cfd(npy_cfd),stat=ierr); jcPmin_cfd(:) = jTmin(:)
    allocate(jcPmax_cfd(npy_cfd),stat=ierr); jcPmax_cfd(:) = jTmax(:)
    allocate(kcPmin_cfd(npz_cfd),stat=ierr); kcPmin_cfd(:) = kTmin(:)
    allocate(kcPmax_cfd(npz_cfd),stat=ierr); kcPmax_cfd(:) = kTmax(:)
    call MPI_bcast(icPmin_cfd,npx_cfd,MPI_INTEGER,source,CPL_INTER_COMM,ierr) !Send
    call MPI_bcast(icPmax_cfd,npx_cfd,MPI_INTEGER,source,CPL_INTER_COMM,ierr) !Send
    call MPI_bcast(jcPmin_cfd,npy_cfd,MPI_INTEGER,source,CPL_INTER_COMM,ierr) !Send
    call MPI_bcast(jcPmax_cfd,npy_cfd,MPI_INTEGER,source,CPL_INTER_COMM,ierr) !Send
    call MPI_bcast(kcPmin_cfd,npz_cfd,MPI_INTEGER,source,CPL_INTER_COMM,ierr) !Send
    call MPI_bcast(kcPmax_cfd,npz_cfd,MPI_INTEGER,source,CPL_INTER_COMM,ierr) !Send

    !Calculate the cell sizes dx,dy & dz
    dx = xg(2,1,1)-xg(1,1,1)
    dy = yg(1,2,1)-yg(1,1,1)
    dz = zg(1,1,2)-zg(1,1,1)

    !Calculate number of cells in overlap region
    ncx_olap = icmax_olap - icmin_olap + 1
    ncy_olap = jcmax_olap - jcmin_olap + 1
    ncz_olap = kcmax_olap - kcmin_olap + 1

    !Broadcast the overlap to CFD on intracommunicator
    call MPI_bcast(ncy_olap,1,MPI_INTEGER,rootid_realm,CPL_REALM_COMM,ierr)

    ! Establish mapping between MD and CFD
    call MPI_Barrier(CPL_WORLD_COMM, ierr)
    call CPL_create_map()

    !Check for grid strectching and terminate process if found
    call check_mesh()

contains

    subroutine check_mesh
        implicit none

        integer, dimension(3) :: sizex, sizey, sizez

        if (ncx*ncy*ncz .lt. maxgridsize) then
            sizex(1) = ncx + 1
            sizex(2) = ncy + 1
            sizex(3) = ncz + 1
            sizey(1) = ncx + 1
            sizey(2) = ncy + 1
            sizey(3) = ncz + 1
            sizez(1) = ncx + 1
            sizez(2) = ncy + 1
            sizez(3) = ncz + 1

            !Define cell sizes dx,dy & dz and check for grid stretching
            ! - - x - -
            dx = xg(2,1,1)-xg(1,1,1)
            dxmax = maxval(xg(2:ncx+1,2:ncy+1,2:ncz+1)-xg(1:ncx,1:ncy,1:ncz))
            dxmin = minval(xg(2:ncx+1,2:ncy+1,2:ncz+1)-xg(1:ncx,1:ncy,1:ncz))
            if (dxmax-dx.gt.0.00001d0 .or.dx-dxmin.gt.0.00001d0) then
                print'(3(a,f15.7))', 'Max dx = ', dxmax, ' dx = ', dx, ' Min dx = ',dxmin
                call error_abort("check_mesh error -  Grid stretching in x not supported")
            endif
            ! - - y - -
            dy = yg(1,2,1)-yg(1,1,1)
            dymax = maxval(yg(2:ncx+1,2:ncy+1,2:ncz+1)-yg(1:ncx,1:ncy,1:ncz))
            dymin = minval(yg(2:ncx+1,2:ncy+1,2:ncz+1)-yg(1:ncx,1:ncy,1:ncz))
            if (dymax-dy.gt.0.00001d0 .or. dy-dymin.gt.0.00001d0) then
                print'(3(a,f15.7))', 'Max dy = ', dymax, ' dy = ', dy, ' Min dy = ',dymin
                call error_abort("check_mesh error -  Grid stretching in y not supported")
            endif

            ! - - z - -
            dz = zg(1,1,2)-zg(1,1,1)
            dzmax = maxval(zg(2:ncx+1,2:ncy+1,2:ncz+1)-zg(1:ncx,1:ncy,1:ncz))
            dzmin = minval(zg(2:ncx+1,2:ncy+1,2:ncz+1)-zg(1:ncx,1:ncy,1:ncz))
            if (dzmax-dz.gt.0.00001d0 .or. dz-dzmin.gt.0.00001d0) then
                print'(3(a,f15.7))', 'Max dz = ', dzmax, ' dz = ', dz, ' Min dz = ',dzmin
                call error_abort("check_mesh error - Grid stretching in z not supported")
            endif

        else
            sizex = 1; sizey = 1; sizez = 1
            sizex(1) = ncx + 1
            sizey(2) = ncy + 1
            sizez(3) = ncz + 1
        endif
  
        ! Check grids are the right size
        if (size(xg,1) .ne. sizex(1) .or. & 
            size(xg,2) .ne. sizex(2) .or. &
            size(xg,3) .ne. sizex(3)) then
            call error_abort('check_mesh error - xg is the wrong size in cpl_cfd_init')
        end if
        if (size(yg,1) .ne. sizey(1) .or. & 
            size(yg,2) .ne. sizey(2) .or. &
            size(yg,3) .ne. sizey(3)) then
            call error_abort('check_mesh error - yg is the wrong size in cpl_cfd_init')
        end if
        if (size(zg,1) .ne. sizez(1) .or. & 
            size(zg,2) .ne. sizez(2) .or. &
            size(zg,3) .ne. sizez(3)) then
            call error_abort('check_mesh error - zg is the wrong size in cpl_cfd_init')
        end if

    end subroutine check_mesh

end subroutine coupler_cfd_init


subroutine CPL_setup_md(icomm_grid, xyzL, xyz_orig)
! ----------------------------------------------------------------------------
!Initialisation routine for coupler module - Every variable is sent and stored
!to ensure both md and cfd region have an identical list of parameters
!
!**Remarks**
!
!Assumes CPL has been initialised `CPL_init` and communicator MD_REALM exists
!
!**Synopsis**
!
!.. code-block:: c
!
!  coupler_md_init(
!                  icomm_grid,
!                  xyzL,
!                  xyz_orig,
!                  )
!
!**Inputs**
!
! - *icomm_grid*
!
!   - Communicator based on MD processor topology returned from a call to MPI_CART_CREATE.
! - *xyzL*
!
!   - MD domain size.
! - *xyz_orig*
!
!   - MD origin.
! - *density*
!
!   - Density of molecules in MD code.
!
! .. sectionauthor::Edward Smith, David Trevelyan, Eduardo Ramos Fernandez
! ------------------------------------
    use mpi
    implicit none

    !Params
    integer, intent(in)                             :: icomm_grid
    real(kind=kind(0.d0)), dimension(3), intent(in) :: xyzL, xyz_orig
    
    !Vars
    integer, dimension(3)                           :: npxyz_md, cart_coords
    logical, dimension(3)                           :: cart_periods
    integer, dimension(:,:), allocatable            :: icoord
    integer                                         :: cart_nprocs, rank

    ! Get number of processors in each direction (periods and coords are not needed)
    call MPI_Cart_get(icomm_grid, 3, npxyz_md, cart_periods, cart_coords, ierr)
    call MPI_Comm_size(icomm_grid, cart_nprocs, ierr)
    
    allocate(icoord(3, cart_nprocs), stat=ierr)

    ! Do the mapping from rank to cartesian coords
    do rank=0, cart_nprocs - 1
        call MPI_Cart_coords(icomm_grid, rank, 3, cart_coords, ierr)
        icoord(1:3, rank + 1) = cart_coords + 1
    enddo

    call coupler_md_init(icomm_grid, icoord, npxyz_md, xyzL, xyz_orig)
    deallocate(icoord)

    !Set flag to register setup is complete correctly
    CPL_setup_complete = 1
    call MPI_Barrier(MPI_COMM_WORLD, ierr)

end subroutine CPL_setup_md

!------------------------------------------------------------------------------
!                              coupler_md_init                               -
!------------------------------------------------------------------------------
!>
!! Initialisation routine for coupler module - Every variable is sent and stored
!! to ensure both md and cfd region have an identical list of parameters
!!
!! - Synopsis
!!
!!  - coupler_md_init(icomm_grid, icoord, npxyz_md, globaldomain)
!!
!! - Input
!!
!!  - icomm_grid
!!   - The MPI communicator setup by the MPI_CART_CREATE command in the 
!!     CFD region (integer)
!!  - icoord
!!   - The three coordinate for each rank in the domain (integer array nproc by 3)
!!  - npxyz_md
!!   - Number of processors in each cartesian dimension (integer array 3)
!!  - globaldomain
!!   - Size of domain in each cartesian dimension (dp real array 3)
!!  - density
!!   - Density of the CFD simulation (dp_real)
!!
!! - Input/Output
!!  - NONE
!!
!! - Output
!!  - NONE
!! 
!! .. sectionauthor:: Edward Smith
!
! ----------------------------------------------------------------------------

! ----------------------------------------------------------------------------
! Initialisation routine for coupler - Every variable is sent and stored
! to ensure both md and cfd region have an identical list of parameters

subroutine coupler_md_init(icomm_grid, icoord, npxyz_md, globaldomain, xyz_orig)
    use mpi
    implicit none

    integer, intent(in)                             :: icomm_grid
    integer, dimension(3), intent(in)               :: npxyz_md
    integer, dimension(:,:), intent(in)             :: icoord
    real(kind=kind(0.d0)),dimension(3),intent(in)   :: globaldomain, xyz_orig

    integer                                         :: i,ib,jb,kb,pcoords(3),source,nproc
    integer,dimension(:),allocatable                :: buf,rank_world2rank_realm,rank_world2rank_cart
    real(kind=kind(0.d0)),dimension(:),allocatable  :: rbuf

    ! Read COUPLER.in input file
    call read_coupler_input()

    ! Duplicate grid communicator for coupler use
    call MPI_comm_dup(icomm_grid, CPL_CART_COMM, ierr)
    call MPI_comm_rank(CPL_CART_COMM, myid_cart, ierr) 
    rank_cart = myid_cart + 1; rootid_cart = 0  
    !Send only from root processor
    if ( myid_realm .eq. rootid_realm ) then
        source=MPI_ROOT
    else
        source=MPI_PROC_NULL
    endif

    ! ================ Exchange and store Data ==============================
    ! Data is stored to the coupler module with the same name in both realms
    ! Note - MPI Broadcast between intercommunicators is only supported by MPI-2

    ! ------------------------ Processor Topology ---------------------------
    ! Receive & Store CFD number of processors
    allocate(buf(3))
    call MPI_bcast(  buf   ,3,MPI_INTEGER,  0   ,CPL_INTER_COMM,ierr) !Receive
    npx_cfd = buf(1); npy_cfd = buf(2); npz_cfd = buf(3)
    nproc_cfd = npx_cfd * npy_cfd * npz_cfd
    deallocate(buf)

    ! Store & Send MD number of processors
    npx_md = npxyz_md(1);   npy_md = npxyz_md(2);   npz_md = npxyz_md(3)    
    nproc_md = npx_md * npy_md * npz_md
    call MPI_bcast(npxyz_md, 3, MPI_INTEGER, source, CPL_INTER_COMM, ierr) !Send

    ! Receive & Store CFD processor rank to coord
    allocate(buf(3*nproc_cfd))
    call MPI_bcast(buf, 3*nproc_cfd, MPI_INTEGER,  0   ,CPL_INTER_COMM,ierr) !Receive
    allocate(rank2coord_cfd(3,nproc_cfd),stat=ierr); rank2coord_cfd = reshape(buf,(/ 3,nproc_cfd /))
    deallocate(buf)

    ! Store & Send MD processor rank to coord
    allocate(rank2coord_md(3,nproc_md),stat=ierr); rank2coord_md = icoord
    iblock_realm=icoord(1,rank_realm); jblock_realm=icoord(2,rank_realm); kblock_realm=icoord(3,rank_realm)
    allocate(buf(3*nproc_md)); buf = reshape(icoord,(/ 3*nproc_md /))
    call MPI_bcast(buf ,3*nproc_md,MPI_INTEGER,source,CPL_INTER_COMM,ierr) !Send
    deallocate(buf)

    ! Receive & Store CFD coordinate to rank mapping
    allocate(buf(nproc_cfd))
    call MPI_bcast(buf, nproc_cfd ,MPI_INTEGER, 0, CPL_INTER_COMM, ierr)   !Receive
    allocate(coord2rank_cfd (npx_cfd,npy_cfd,npz_cfd))
    coord2rank_cfd = reshape(buf,(/ npx_cfd,npy_cfd,npz_cfd /))
    deallocate(buf)

    ! Setup MD mapping from coordinate to rank, 
    ! Store & Send MD mapping from coordinate to rank to CFD
    allocate(coord2rank_md(npx_md,npy_md,npz_md))
    do ib = 1,npx_md
    do jb = 1,npy_md
    do kb = 1,npz_md
        pcoords = (/ ib, jb, kb /)-1
        call MPI_Cart_rank(CPL_CART_COMM,pcoords,i,ierr)
        coord2rank_md(ib,jb,kb) = i + 1
    enddo
    enddo
    enddo
    allocate(buf(nproc_md)); buf = reshape(coord2rank_md, (/ nproc_md /) )
    call MPI_bcast(coord2rank_md,nproc_md,MPI_INTEGER,source,CPL_INTER_COMM,ierr)   !Send
    deallocate(buf)

    ! Setup MD mapping between realm & world rank 
    allocate(rank_mdrealm2rank_world(nproc_md))
    allocate(rank_world2rank_realm(nproc_world))
    call CPL_rank_map(CPL_REALM_COMM,rank_realm,nproc, & 
                        rank_mdrealm2rank_world,rank_world2rank_realm,ierr)

    !World to rank is the same on both realms
    allocate(rank_world2rank_cfdrealm(nproc_world))
    allocate(rank_world2rank_mdrealm(nproc_world))
    rank_world2rank_cfdrealm = rank_world2rank_realm
    rank_world2rank_mdrealm  = rank_world2rank_realm

    ! Receive & Store CFD_mapping from realm to local rank
    allocate(rank_cfdrealm2rank_world(nproc_cfd))
    call MPI_bcast(rank_cfdrealm2rank_world,nproc_cfd,MPI_integer,0,CPL_INTER_COMM,ierr)    !Receive

    ! Send MD mapping from realm to local rank to CFD
    call MPI_bcast(rank_mdrealm2rank_world,nproc_md,MPI_integer,source,CPL_INTER_COMM,ierr)  !send

    ! Setup MD mapping between cartesian topology & world rank 
    allocate(rank_mdcart2rank_world(nproc_md))
    allocate(rank_world2rank_cart(nproc_world))
    call CPL_rank_map(CPL_CART_COMM,rank_cart,nproc, & 
                        rank_mdcart2rank_world,rank_world2rank_cart,ierr)

    !World to rank is the same on both realms cart
    allocate(rank_world2rank_cfdcart(nproc_world))
    allocate(rank_world2rank_mdcart(nproc_world))
    rank_world2rank_cfdcart = rank_world2rank_cart
    rank_world2rank_mdcart  = rank_world2rank_cart

    ! Receive & Store CFD_mapping from cart to local rank
    allocate(rank_cfdcart2rank_world(nproc_cfd))
    call MPI_bcast(rank_cfdcart2rank_world,nproc_cfd,MPI_integer,0,CPL_INTER_COMM,ierr) !Receive

    ! Send MD mapping from cart to local rank to CFD
    call MPI_bcast(rank_mdcart2rank_world,nproc_md,MPI_integer,source,CPL_INTER_COMM,ierr)   !send

    ! ------------------ Receive CFD grid extents ------------------------------

    ! Receive & store CFD domain size
    allocate(rbuf(3))
    call MPI_bcast(rbuf,3,MPI_double_precision,0,CPL_INTER_COMM,ierr)               !Receive
    xL_cfd = rbuf(1); yL_cfd = rbuf(2); zL_cfd = rbuf(3)
    deallocate(rbuf)

    ! Receive & store CFD domain origin coords
    allocate(rbuf(3))
    call MPI_bcast(rbuf,3,MPI_double_precision,0,CPL_INTER_COMM,ierr)               !Receive
    x_orig_cfd = rbuf(1); y_orig_cfd = rbuf(2); z_orig_cfd = rbuf(3)
    deallocate(rbuf)

    ! Store & send MD domain size
    xL_md = globaldomain(1); yL_md = globaldomain(2); zL_md = globaldomain(3) 
    call MPI_bcast(globaldomain,3,MPI_double_precision,source,CPL_INTER_COMM,ierr)  !Send

    ! Store & send MD domain origin coords
    x_orig_md = xyz_orig(1); y_orig_md = xyz_orig(2); z_orig_md = xyz_orig(3)
    call MPI_bcast(xyz_orig, 3, MPI_double_precision, source, CPL_INTER_COMM, ierr)  !Send

    ! Receive & Store global CFD grid extents
    allocate(buf(6))
    call MPI_bcast(buf, 6, MPI_INTEGER, 0, CPL_INTER_COMM,ierr) !Send
    icmin = buf(1); icmax = buf(2)
    jcmin = buf(3); jcmax = buf(4)
    kcmin = buf(5); kcmax = buf(6)
    deallocate(buf)

    !Set overlap, constraint and boundary to no. CFD cells if full overlap 
    if (CPL_full_overlap) then
        icmin_olap = icmin; icmax_olap = icmax
        jcmin_olap = jcmin; jcmax_olap = jcmax
        kcmin_olap = kcmin; kcmax_olap = kcmax
        
        icmin_cnst = icmin; icmax_cnst = icmax
        jcmin_cnst = jcmin; jcmax_cnst = jcmax 
        kcmin_cnst = kcmin; kcmax_cnst = kcmax

        icmin_bnry = icmin; icmax_bnry = icmax
        jcmin_bnry = jcmin; jcmax_bnry = jcmax 
        kcmin_bnry = kcmin; kcmax_bnry = kcmax
    endif

    ! Receive & Store array of global number of cells in CFD
    allocate(buf(3))
    call MPI_bcast(buf, 3, MPI_INTEGER, 0, CPL_INTER_COMM,ierr) !Receive
    ncx = buf(1); ncy = buf(2); ncz = buf(3)
    deallocate(buf)

    ! Receive & Store array of global grid points
    if (ncx*ncy*ncz .lt. maxgridsize) then
        allocate(xg(ncx+1,ncy+1,ncz+1),yg(ncx+1,ncy+1,ncz+1),zg(ncx+1,ncy+1,ncz+1))
    else
        allocate(xg(ncx+1,1,1),yg(1,ncy+1,1),zg(1,1,ncz+1))
    endif
    call MPI_bcast(xg,size(xg),MPI_double_precision,0,CPL_INTER_COMM,ierr) !Receive
    call MPI_bcast(yg,size(yg),MPI_double_precision,0,CPL_INTER_COMM,ierr) !Receive
    call MPI_bcast(zg,size(zg),MPI_double_precision,0,CPL_INTER_COMM,ierr) !Receive 

    ! Receive & Store local (processor) CFD grid extents
    allocate(icPmin_cfd(npx_cfd)); allocate(icPmax_cfd(npx_cfd));  
    allocate(jcPmin_cfd(npy_cfd)); allocate(jcPmax_cfd(npy_cfd));
    allocate(kcPmin_cfd(npz_cfd)); allocate(kcPmax_cfd(npz_cfd));
    call MPI_bcast(icPmin_cfd,npx_cfd,MPI_INTEGER,0,CPL_INTER_COMM,ierr) !Receive
    call MPI_bcast(icPmax_cfd,npx_cfd,MPI_INTEGER,0,CPL_INTER_COMM,ierr) !Receive
    call MPI_bcast(jcPmin_cfd,npy_cfd,MPI_INTEGER,0,CPL_INTER_COMM,ierr) !Receive
    call MPI_bcast(jcPmax_cfd,npy_cfd,MPI_INTEGER,0,CPL_INTER_COMM,ierr) !Receive
    call MPI_bcast(kcPmin_cfd,npz_cfd,MPI_INTEGER,0,CPL_INTER_COMM,ierr) !Receive
    call MPI_bcast(kcPmax_cfd,npz_cfd,MPI_INTEGER,0,CPL_INTER_COMM,ierr) !Receive

    !Calculate the cell sizes dx,dy & dz
    dx = xg(2,1,1)-xg(1,1,1)
    dy = yg(1,2,1)-yg(1,1,1) ! yL_cfd/ncy
    dz = zg(1,1,2)-zg(1,1,1)

    !Define number of cells in overlap region
    ncx_olap = icmax_olap - icmin_olap + 1
    ncy_olap = jcmax_olap - jcmin_olap + 1
    ncz_olap = kcmax_olap - kcmin_olap + 1

    ! Establish mapping between MD an CFD
    call MPI_Barrier(CPL_WORLD_COMM, ierr)
    call CPL_create_map()

end subroutine coupler_md_init

!-----------------------------------------------------------------------------

subroutine CPL_set_timing(initialstep, nsteps, dt)
! - *nsteps*
!
!   - Number of steps in MD simulation.
! - *initialstep*
!
!   - Initial steps in MD simulation.
! - *dt*
!
!   - Timestep in MD simulation.
    use mpi
    implicit none

    integer, intent(in) :: initialstep, nsteps
    real(kind=kind(0.d0)), intent(in) :: dt

    integer :: Nsteps_MDperCFD, source
    real(kind=kind(0.d0)) :: elapsedtime

    call error_abort("CPL_set_timing is depricated and should not be called")

    if ( myid_realm .eq. rootid_realm ) then
        source=MPI_ROOT
    else
        source=MPI_PROC_NULL
    endif

    ! ------------------ Timesteps and iterations ------------------------------
    if (realm .eq. md_realm) then
        ! Receive & store CFD nsteps and dt_cfd
        !print*, 'CFDbefore-nsteps: ', nsteps_cfd
        call MPI_bcast(nsteps_cfd, 1, MPI_integer, 0, CPL_INTER_COMM, ierr) !Receive
        !print*, 'CFDafter-nsteps: ', nsteps_cfd
        call MPI_bcast(dt_cfd, 1, MPI_double_precision, 0, CPL_INTER_COMM, ierr) !Receive
        !print*, 'CFD-dt: ', dt_cfd

        ! Store & send MD timestep to dt_md
        dt_MD = dt
        call MPI_bcast(dt_MD, 1, MPI_double_precision, source, CPL_INTER_COMM, ierr) !Send
        nsteps_MD = nsteps
        call MPI_bcast(nsteps_MD, 1, MPI_integer, source, CPL_INTER_COMM, ierr) !Send

    elseif (realm .eq. cfd_realm) then
    ! ------------------ Timesteps and iterations ------------------------------
        ! Store & send CFD nsteps and dt_cfd
        nsteps_cfd = nsteps
        !print*, 'CFD3() nsteps: ', nsteps_cfd

        call MPI_bcast(nsteps_cfd ,1,MPI_integer,source,CPL_INTER_COMM,ierr) !Send
        dt_cfd = dt
        call MPI_bcast(dt_cfd,1,MPI_double_precision,source,CPL_INTER_COMM,ierr) !Send

        ! Receive & store MD timestep dt_md
        call MPI_bcast(dt_md,1,MPI_double_precision,0,CPL_INTER_COMM,ierr)      !Receive
        call MPI_bcast(nsteps_md,1,MPI_integer,     0,CPL_INTER_COMM,ierr)      !Receive
    endif

    !Set number of MD timesteps per CFD using ratio of timestep or coupler value
    if (timestep_ratio .eq. VOID) then
        Nsteps_MDperCFD = int(dt_cfd/dt_MD)
    else 
        Nsteps_MDperCFD = timestep_ratio
    endif
    Nsteps_coupled = Nsteps_cfd

    !Set number of steps in MD simulation and final time elapsed
    Nsteps_md   = initialstep + Nsteps_cfd * Nsteps_MDperCFD
    elapsedtime = Nsteps_md * dt_MD

    if (rank_realm .eq. 1 .and. (output_mode .ne. QUIET)) then 
        print*, 'Nsteps in CFD is ', Nsteps_cfd
        print*, 'Nsteps in MD reset from ', Nsteps, ' to ', Nsteps_md
        print*, 'Total simulation time will be ', elapsedtime, '.'
    endif 

    call MPI_Barrier(MPI_COMM_WORLD, ierr)

    !Set corrected nsteps returned to MD
    !if (realm .eq. md_realm) Nsteps = Nsteps_md

end subroutine CPL_set_timing


!=============================================================================
!! Establish for all MD processors the mapping (if any) 
!! to coupled CFD processors
!-----------------------------------------------------------------------------

subroutine CPL_create_map()
    use mpi
    implicit none

    ! Check (as best as one can) that the inputs will work
    call check_config_feasibility()

    ! Get ranges of cells on each MD processor
    call get_md_cell_ranges()

    ! Get overlapping mapping for MD to CFD
    call get_overlap_blocks()

    ! Setup overlap communicators
    call prepare_overlap_comms()

    ! Setup graph topology
    call set_overlap_topology()

contains

subroutine check_config_feasibility()
    implicit none

    integer :: ival
    real(kind(0.d0)) :: rval, rtoler
    character(len=256) :: string
    logical :: error

    ! Check that CFD and MD domains agree in x and z directions
    rtoler = 1.d-4
    rval = 0.d0
    error = .false.
    rval = rval + abs(xL_md - xL_cfd)
    if (rval .gt. rtoler) error = .true.
    rval = rval + abs(zL_md - zL_cfd)
    if (rval .gt. rtoler) error = .true.
    if (error) then
        string = "CPL_create_map error - MD/CFD domain sizes do not match in both x and z "   // &
                 "directions. Aborting simulation. "
        print*, "xL_md = ",  xL_md
        print*, "xL_cfd = ", xL_cfd
        print*, "zL_md = ",  zL_md
        print*, "zL_cfd = ", zL_cfd
        call error_abort(string)
    end if


    !Check CFD and MD cell range
    if (npx_md .lt. npx_cfd) then
        call error_abort("CPL_create_map error - number of MD " // &
                         "processors must be greater than or equal" // &
                         " to CFD processors in x")
    endif
    if (npy_md .lt. npy_cfd) then
        call error_abort("CPL_create_map error - number of MD " // &
                         "processors must be greater than or equal" // &
                         " to CFD processors in y")
    endif
    if (npz_md .lt. npz_cfd) then
        call error_abort("CPL_create_map error - number of MD " // &
                         "processors must be greater than or equal" // &
                         " to CFD processors in z")
    endif

    ! Check there is only one overlap CFD proc in y
    ival = nint( dble(ncy) / dble(npy_cfd) )
    if (ncy_olap .gt. ival) then

        if (CPL_full_overlap .eqv. .false.) then
            string = "CPL_create_map error - This coupler will not work if there is more than one "// &
                     "CFD processor (y-coordinate) in the overlapping region. "       // &
                     "Aborting simulation."
            call error_abort(string)
        endif

    end if

    ! Check that MD processor size is an integer multiple of CFD cell size
    ! This test doesn't work if xL_xyz is slightly less than a multiple of dxyz
    ! We avoid this by adding the required tolerence to the mod and taking away after
    rtoler = 1.d-4
    rval = 0.d0
    error = .false.
    rval = abs(mod( xL_md+rtoler, dx )-rtoler)
    if (rval .gt. rtoler) error = .true.
    rval = abs(mod( yL_md+rtoler, dy )-rtoler)
    if (rval .gt. rtoler) error = .true.
    rval = abs(mod( zL_md+rtoler, dz )-rtoler)
    if (rval .gt. rtoler) error = .true.

    if (error) then
        print'(6(a,f10.5))', ' xL_md/dx = ',xL_md/dx, ' dx =', dx, & 
                             ' yL_md/dy = ',yL_md/dy, ' dy =', dy, &
                             ' zL_md/dz = ',zL_md/dz, ' dz =', dz
        string = "CPL_create_map error - MD region lengths must be an integer number of CFD " // &
                 "cell sizes (i.e. xL_md must be an integer multiple " // &
                 "of dx, etc. ), aborting simulation."
        call error_abort(string)

    end if 

    ! Check whether ncx,ncy,ncz are an integer multiple of npx_md, etc.
    ! - N.B. no need to check ncy/npy_md.
    ival = 0
    ival = ival + mod(ncx,npx_cfd)
    ival = ival + mod(ncy,npy_cfd)
    ival = ival + mod(ncz,npz_cfd)
    ival = ival + mod(ncx,npx_md)
    ival = ival + mod(ncz,npz_md)
    if (ival.ne.0) then 

        string = "CPL_create_map error - The number of cells in the cfd domain is not an "    // &
                 "integer multiple of the number of processors in "    // &
                 "the x and z directions. Aborting simulation." 
        call error_abort(string)

    end if

    ! Check that the MD region is large enough to cover overlap
    xL_olap = ncx_olap * dx
    yL_olap = ncy_olap * dy
    zL_olap = ncz_olap * dz
    ! Tolerance of half a cell width
    if (xL_md .lt. (xL_olap - dx/2.d0)) then
        print'(2(a,f15.6),a, i8,a,f15.6)',  "xL_md= ", xL_md, " xL_olap= ", xL_olap, " ncx_olap= ", ncx_olap , " dx= ",dx
        string = "CPL_create_map error - Overlap region is larger than the MD region. "       // &
                 "Aborting simulation."
        call error_abort(string)
    elseif (yL_md .lt. (yL_olap - dy/2.d0)) then
        print'(2(a,f15.6),a, i8,a,f15.6)',  "yL_md= ", yL_md, " yL_olap= ", yL_olap, " ncy_olap= ", ncy_olap , " dy= ",dy
        string = "CPL_create_map error - Overlap region is larger than the MD region. "       // &
                 "Aborting simulation."
        call error_abort(string)
    elseif (zL_md .lt. (zL_olap - dz/2.d0)) then
        print'(2(a,f15.6),a, i8,a,f15.6)',  "zL_md= ", zL_md, " zL_olap= ", zL_olap, " ncz_olap= ", ncz_olap , " dz= ",dz
        string = "CPL_create_map error - Overlap region is larger than the MD region. "       // &
                 "Aborting simulation."
        call error_abort(string)

    end if

    ! Check overlap lower limit are not greater than upper limits
    ival = 0
    if (icmin_olap.gt.icmax_olap) ival = ival + 1        
    if (jcmin_olap.gt.jcmax_olap) ival = ival + 1        
    if (kcmin_olap.gt.kcmax_olap) ival = ival + 1        
    if (ival.ne.0) then
        string = "CPL_create_map error - Overlap region lower limits are greater than upper"  // &
                 " limits for some directions. Aborting simulation."
        call error_abort(string)
    end if

    ! Check overlap cells are non-negative
    ival = 0
    if (icmin_olap.lt.0) ival = ival + 1        
    if (icmax_olap.lt.0) ival = ival + 1        
    if (jcmin_olap.lt.0) ival = ival + 1        
    if (jcmax_olap.lt.0) ival = ival + 1        
    if (kcmin_olap.lt.0) ival = ival + 1        
    if (kcmax_olap.lt.0) ival = ival + 1        
    if (ival.ne.0) then
        string = "CPL_create_map error - Overlap region limits contains a negative index"  // &
                 ". Aborting simulation."
        call error_abort(string)
    end if

    ! Check overlap cells are within CFD extents
    ival = 0
    if (icmin_olap.lt.icmin) ival = ival + 1        
    if (icmax_olap.gt.icmax) ival = ival + 1        
    if (jcmin_olap.lt.jcmin) ival = ival + 1        
    if (jcmax_olap.gt.jcmax) ival = ival + 1        
    if (kcmin_olap.lt.kcmin) ival = ival + 1        
    if (kcmax_olap.gt.kcmax) ival = ival + 1        
    if (ival.ne.0) then
        print'(a,6i10)', 'ijkcmin,ijkcmax = ' , icmin, icmax & 
                            , jcmin, jcmax &
                            , kcmin, kcmax
        print'(a,6i10)', 'olap extents    = ', icmin_olap, icmax_olap, jcmin_olap, &
                                   jcmax_olap, kcmin_olap, kcmax_olap
        string = "CPL_create_map error - Overlap region has been specified outside of the "  // &
                 "CFD region. Aborting simulation."
        call error_abort(string)
    end if
    
    ! Check MD/CFD ratios are integers in x and z
    if (mod(npx_md,npx_cfd) .ne. 0) then

        print'(a,i8,a,i8)', ' number of MD processors in x ', npx_md,     & 
                            ' number of CFD processors in x ', npx_cfd

        call error_abort("CPL_create_map error - get_overlap_blocks error - number of MD "    // & 
                         "processors must be an integer multiple "// &
                         "of number of CFD processors in x")

    elseif (mod(npy_md,npy_cfd) .ne. 0 .and. CPL_full_overlap) then

        print'(a,i8,a,i8)', ' number of MD processors in y ', npy_md,     & 
                            ' number of CFD processors in y ', npy_cfd

        call error_abort("CPL_create_map error - get_overlap_blocks error - number of MD "    // &
                         "processors must be an integer multiple "// &
                         "of number of CFD processors in y")


    elseif (mod(npz_md,npz_cfd) .ne. 0) then

        print'(a,i8,a,i8)', ' number of MD processors in z ', npz_md,     & 
                            ' number of CFD processors in z ', npz_cfd

        call error_abort("CPL_create_map error - get_overlap_blocks error - number of MD "    // &
                         "processors must be an integer multiple "// &
                         "of number of CFD processors in z")

    endif

end subroutine check_config_feasibility


!------------------------------------------------------------
!Calculate processor cell ranges of MD code on all processors
!------------------------------------------------------------------------------
!                      GET MD CELL RANGES                                     -
!------------------------------------------------------------------------------
!>
!! Store the minimum and maximum CFD cell coordinates that overlap each
!! MD processor.   
!!
!! - Synopsis
!!
!!  - get_md_cell_ranges()
!!
!! - Input
!!  - NONE 
!! - Input/Output
!!  - NONE 
!! - Output
!!  - NONE 
!! 
!! .. sectionauthor:: David Trevelyan 
subroutine get_md_cell_ranges()
    implicit none

    integer :: n
    integer :: olap_jmin_mdproc
    integer :: ncxl, ncyl, nczl
    integer :: ncy_mdonly, ncy_md, ncyP_md
    integer :: funit

    allocate(icPmin_md(npx_md)); icPmin_md = VOID
    allocate(jcPmin_md(npy_md)); jcPmin_md = VOID
    allocate(kcPmin_md(npz_md)); kcPmin_md = VOID
    allocate(icPmax_md(npx_md)); icPmax_md = VOID
    allocate(jcPmax_md(npy_md)); jcPmax_md = VOID
    allocate(kcPmax_md(npz_md)); kcPmax_md = VOID

    ! - - x - -
    ncxl = ceiling(dble(ncx)/dble(npx_md))
    do n=1,npx_md
        icPmax_md(n) = n * ncxl
        icPmin_md(n) = icPmax_md(n) - ncxl + 1
    end do  

    ! - - y - -
    if (CPL_full_overlap) then
        ncyl = ceiling(dble(ncy)/dble(npy_md))
        do n=1,npy_md
            jcPmax_md(n) = n * ncyl
            jcPmin_md(n) = jcPmax_md(n) - ncyl + 1
        end do
    else
        ncy_md   = nint(yL_md/dy)
        ncy_mdonly = ncy_md - ncy_olap
        ncyP_md = ncy_md / npy_md
        olap_jmin_mdproc = npy_md - ceiling(dble(ncy_olap)/dble(ncyP_md)) + 1
        do n = olap_jmin_mdproc,npy_md
            jcPmax_md(n) = n * ncyP_md - ncy_mdonly
            jcPmin_md(n) = jcPmax_md(n) - ncyP_md + 1
            if (jcPmin_md(n).le.0) jcPmin_md(n) = 1
        end do  
    endif

    ! - - z - -
    nczl = ceiling(dble(ncz)/dble(npz_md))
    do n=1,npz_md
        kcPmax_md(n) = n * nczl
        kcPmin_md(n) = kcPmax_md(n) - nczl + 1
    end do

    if (myid_world .eq. rootid_world) then

        funit = CPL_new_fileunit()
        open(funit, file="cpl/map_MD", action="write", status="replace")
        write(funit,*) ''
        write(funit,*) '==========================================='
        write(funit,*) '------------ M D   M A P ------------------'
        write(funit,*) '==========================================='
        write(funit,*) 'npx_md = ', npx_md
        write(funit,*) 'ncx    = ', ncx
        write(funit,*) 'ncxl   = ', ncxl
        write(funit,*) '-------------------------------------------'
        write(funit,*) '  icoord_md     icPmin_md     icPmax_md    '
        write(funit,*) '-------------------------------------------'
        do n=1,npx_md
            write(funit,'(1x,3i11)') n, icPmin_md(n), icPmax_md(n)
        end do  
        write(funit,*) '-------------------------------------------'
        write(funit,*) 'npy_md     = ', npy_md
        write(funit,*) 'ncy_md     = ', ncy_md
        write(funit,*) 'ncyP_md    = ', ncyP_md 
        write(funit,*) 'ncy_olap   = ', ncy_olap
        write(funit,*) 'ncy_mdonly = ', ncy_mdonly
        write(funit,*) 'olap_jmin_mdproc = ', olap_jmin_mdproc
        write(funit,*) 'dy         = ', dy
        write(funit,*) '-------------------------------------------'
        write(funit,*) '  jcoord_md     jcPmin_md       jcPmax_md  '
        write(funit,*) '-------------------------------------------'
        do n = 1,npy_md 
            write(funit,'(1x,3i11)') n, jcPmin_md(n), jcPmax_md(n)
        end do
        write(funit,*) '-------------------------------------------'
        write(funit,*) 'npz_md = ', npz_md
        write(funit,*) 'ncz    = ', ncz
        write(funit,*) 'nczl   = ', nczl
        write(funit,*) '-------------------------------------------'
        write(funit,*) '  kcoord_md     kcPmin_md       kcPmax_md  '
        write(funit,*) '-------------------------------------------'
        do n=1,npz_md
            write(funit,'(1x,3i11)') n, kcPmin_md(n), kcPmax_md(n)
        end do
        write(funit,*) '-------------------------------------------'
        close(funit, status="keep")

    endif

    !Sanity Check min not less than max
    if (icPmin_md(iblock_realm) .gt. icPmax_md(iblock_realm) .and. icPmax_md(iblock_realm) .ne. VOID) & 
        call error_abort("get_md_cell_ranges error - mapping failure imin greater than imax")
    do n = olap_jmin_mdproc,npy_md
        if (jcPmin_md(n) .gt. jcPmax_md(n) .and. jcPmax_md(n) .ne. VOID) &
            call error_abort("get_md_cell_ranges error - mapping failure jmin greater than jmax")
    enddo
    if (kcPmin_md(kblock_realm) .gt. kcPmax_md(kblock_realm) .and. kcPmax_md(kblock_realm) .ne. VOID) & 
        call error_abort("get_md_cell_ranges error - mapping failure kmin greater than kmax")

end subroutine get_md_cell_ranges

!------------------------------------------------------------
!Calculate processor overlap between CFD/MD on all processors

!------------------------------------------------------------------------------
!                          GET OVERLAP BLOCKS                                 -
!------------------------------------------------------------------------------
!>
!! Store MD processor coordinates that overlap each CFD processor coordinate. 
!!
!! - Synopsis
!!
!!  - get_overlap_blocks()
!!
!! - Input
!!  - NONE 
!! - Input/Output
!!  - NONE 
!! - Output
!!  - NONE 
!! 
!! .. sectionauthor:: David Trevelyan
subroutine get_overlap_blocks()
    implicit none

    integer             :: i,n,endproc,nolapsx,nolapsy,nolapsz
    integer :: funit

    real(kind(0.d0))    :: xLl_md, yLl_md, zLl_md, yLl_cfd

    xL_olap = ncx_olap * dx 
    yL_olap = ncy_olap * dy 
    zL_olap = ncz_olap * dz 

    xLl_md  = xL_md / npx_md
    yLl_md  = yL_md / npy_md
    zLl_md  = zL_md / npz_md

    if (realm .eq. md_realm) then
        xLl = xLl_md; yLl = yLl_md ; zLl = zLl_md 
    endif

    nolapsx = nint( dble( npx_md ) / dble( npx_cfd ) )
    if (CPL_full_overlap) then
        nolapsy = nint( dble( npy_md ) / dble( npy_cfd ) )
    else
        nolapsy = ceiling ( yL_olap / yLl_md ) 
    endif
    nolapsz = nint( dble( npz_md ) / dble( npz_cfd ) )

    !Get cartesian coordinate of overlapping md cells & cfd cells
    allocate(cfd_icoord2olap_md_icoords(npx_cfd,nolapsx)) 
    allocate(cfd_jcoord2olap_md_jcoords(npy_cfd,nolapsy)) 
    allocate(cfd_kcoord2olap_md_kcoords(npz_cfd,nolapsz)) 
    cfd_icoord2olap_md_icoords = VOID
    cfd_jcoord2olap_md_jcoords = VOID
    cfd_kcoord2olap_md_kcoords = VOID

    ! - - x - -
    do n = 1,npx_cfd
    do i = 1,nolapsx    
        cfd_icoord2olap_md_icoords(n,i) = (n-1)*nolapsx + i
    end do
    end do

    ! - - y - -
    if (CPL_full_overlap) then
        do n = 1,npy_cfd
        do i = 1,nolapsy    
            cfd_jcoord2olap_md_jcoords(n,i) = (n-1)*nolapsy + i
        end do
        end do
    else
        yLl_cfd = yL_cfd/npy_cfd
        endproc = ceiling(yL_olap/yLl_cfd)
        if (endproc .gt. npy_cfd) then
            print*, "get_overlap_blocks warning - in get_overlap_blocks"
            print*, "-- top processor in CFD greater than number of processors."
            print*, "This may be correct if some MD domain exists above CFD."
            endproc = npy_cfd
            nolapsy = 1
            print*, endproc, nolapsy
        endif
        do n = 1,endproc
        do i = 1,nolapsy
            cfd_jcoord2olap_md_jcoords(n,i) =   (n-1)*nolapsy + i &
                                              + (npy_md - nolapsy)
        end do
        end do
    endif


    ! - - z - -
    do n = 1,npz_cfd
    do i = 1,nolapsz    
        cfd_kcoord2olap_md_kcoords(n,i) = (n-1)*nolapsz + i
    end do
    end do


    if (myid_world .eq. rootid_world) then 
    
        funit = CPL_new_fileunit()
        open(funit, file="cpl/map_CFD", action="write", status="replace")
        write(funit,*) ''
        write(funit,*) '==========================================='
        write(funit,*) '------------ C F D   M A P ----------------'
        write(funit,*) '==========================================='
        write(funit,*) 'npx_cfd = ', npx_cfd
        write(funit,*) 'nolapsx = ', nolapsx
        write(funit,*) '-------------------------------------------'
        write(funit,*) '  icoord_cfd       olapmin     olapmax     ' 
        write(funit,*) '-------------------------------------------'
        do n=1,npx_cfd
            write(funit,'(1x,3i11)') n,               &
                  cfd_icoord2olap_md_icoords(n,1),         &
                  cfd_icoord2olap_md_icoords(n,nolapsx)
        end do  
        write(funit,*) '-------------------------------------------'

        write(funit,*) 'npy_cfd = ', npy_cfd
        write(funit,*) 'nolapsy = ', nolapsy
        write(funit,*) '-------------------------------------------'
        write(funit,*) '  jcoord_cfd       olapmin     olapmax     ' 
        write(funit,*) '-------------------------------------------'
        do n=1,npy_cfd
            write(funit,'(1x,3i11)') n,               &
                  cfd_jcoord2olap_md_jcoords(n,1),         &
                  cfd_jcoord2olap_md_jcoords(n,nolapsy)
        end do  
        write(funit,*) '-------------------------------------------'

        write(funit,*) 'npz_cfd = ', npz_cfd
        write(funit,*) 'nolapsz = ', nolapsz
        write(funit,*) '-------------------------------------------'
        write(funit,*) '  kcoord_cfd       olapmin     olapmax     ' 
        write(funit,*) '-------------------------------------------'
        do n=1,npz_cfd
            write(funit,'(1x,3i11)') n,               &
                  cfd_kcoord2olap_md_kcoords(n,1),         &
                  cfd_kcoord2olap_md_kcoords(n,nolapsz)
        end do  
        write(funit,*) '-------------------------------------------'
        close(funit, status="keep")

    endif

end subroutine get_overlap_blocks

!subroutine intersect_comm
!   use coupler_module
!   use mpi
!   implicit none

!   integer :: n
!   integer,dimension(3)   :: pcoords

    !Create communicator for all intersecting processors
!   if (realm .eq. cfd_realm) then
!       map%n = npx_md/npx_cfd
!       allocate(map%rank_list(map%n))
        !Get rank(s) of overlapping MD processor(s)
!       do n = 1,map%n
!           pcoords(1)=cfd_icoord2olap_md_icoords(rank2coord_cfd(1,rank_realm),n)
!           pcoords(2)=cfd_jcoord2olap_md_jcoords(rank2coord_cfd(2,rank_realm),n)
!           pcoords(3)=cfd_kcoord2olap_md_kcoords(rank2coord_cfd(3,rank_realm),n)
!           if (any(pcoords(:).eq.VOID)) then
!               map%n = 0; map%rank_list(:) = VOID
!           else
!               map%rank_list(n) = coord2rank_md(pcoords(1),pcoords(2),pcoords(3))
!           endif
!           write(250+rank_realm,'(2a,6i5)'), 'overlap',realm_name(realm),rank_realm,map%n,map%rank_list(n),pcoords
!       enddo
!   else if (realm .eq. md_realm) then
!       map%n = 1
!       allocate(map%rank_list(map%n))  
        !Get rank of overlapping CFD processor
!       pcoords(1) = rank2coord_md(1,rank_realm)*(dble(npx_cfd)/dble(npx_md))
!       pcoords(2) = npy_cfd !rank2coord_md(2,rank_realm)*(dble(npy_cfd)/dble(npy_md))
!       pcoords(3) = rank2coord_md(3,rank_realm)*(dble(npz_cfd)/dble(npz_md))
!!      map%rank_list(1) = coord2rank_cfd(pcoords(1),pcoords(2),pcoords(3))
!       write(300+rank_realm,'(2a,6i5)'), 'overlap',realm_name(realm),rank_realm,map%n,map%rank_list(1),pcoords
!   endif

!end subroutine intersect_comm

!=========================================================================

!------------------------------------------------------------------------------
!                    PREPARE OVERLAP COMMS                                    -
!------------------------------------------------------------------------------
!>
!! .. sectionauthor:: David Trevelyan 
!! Splits the world communicator into "overlap" communicators. Each overlap 
!! communicator consists of a CFD root processor and the MD processors which
!! lie on top of it in the domain. 
!!
!! - Synopsis
!!
!!  - prepare_overlap_comms()
!!
!! - Input
!!  - NONE
!! - Input/Output
!!  - NONE 
!! - Output
!!  - NONE 
!! 
subroutine prepare_overlap_comms()
    use mpi
    implicit none

    !General idea:
    !   1. loop over cfd cart ranks
    !   2. find cfd cart coords from cfd cart rank
    !   3. find overlapping md cart coords (from cfd_icoord2olap_md_icoords)
    !   4. find md cart rank from md cart coords (coord2rank_md) 
    !   5. find md world rank from md cart rank (rank_mdcart2rank_world)
    !   6. set group(md_world_rank) to cfd cart rank
    !   7. split world comm according to groups
    !      if group(world_rank) == 0, set olap_comm to null 

    integer :: i,j,k,ic,jc,kc,n
    integer :: trank_md, trank_cfd, trank_world, nolap, color
    integer, dimension(:), allocatable :: mdicoords, mdjcoords, mdkcoords
    integer, parameter :: olap_null = -666
    integer :: group(nproc_world)
    integer :: cfdcoord(3)
    integer :: tempsize
    !integer :: local_leader, remote_leader, comm_size
    !integer :: jbuf(2), ibuf(2)

    tempsize = size(cfd_icoord2olap_md_icoords,2)
    allocate(mdicoords(tempsize))
    tempsize = size(cfd_jcoord2olap_md_jcoords,2)
    allocate(mdjcoords(tempsize))
    tempsize = size(cfd_kcoord2olap_md_kcoords,2)
    allocate(mdkcoords(tempsize))

    allocate(olap_mask(nproc_world))

    !Set default values, must be done because coord2rank_md cannot
    !take "null" coordinates.
    group(:) = olap_null
    olap_mask(:) = .false.
    nolap = 0

    ! Every process loop over all cfd ranks
    do trank_cfd = 1,nproc_cfd

        ! Get cart coords of cfd rank
        cfdcoord(:)  = rank2coord_cfd(:,trank_cfd)

        ! Get md cart coords overlapping cfd proc
        mdicoords(:) = cfd_icoord2olap_md_icoords(cfdcoord(1),:)
        mdjcoords(:) = cfd_jcoord2olap_md_jcoords(cfdcoord(2),:)
        mdkcoords(:) = cfd_kcoord2olap_md_kcoords(cfdcoord(3),:)

        ! Set group and olap_mask for CFD processor if it overlaps
        if (any(mdicoords.ne.olap_null) .and. &
            any(mdjcoords.ne.olap_null) .and. &
            any(mdkcoords.ne.olap_null)) then

            trank_world = rank_cfdcart2rank_world(trank_cfd)
            olap_mask(trank_world) = .true.
            group    (trank_world) = trank_cfd

        end if

        ! Set group and olap_mask for MD processors
        do i = 1,size(mdicoords)
        do j = 1,size(mdjcoords)
        do k = 1,size(mdkcoords)

            ic = mdicoords(i)
            jc = mdjcoords(j)
            kc = mdkcoords(k)

            if (any((/ic,jc,kc/).eq.olap_null)) cycle

            trank_md = coord2rank_md(ic,jc,kc)
            trank_world = rank_mdcart2rank_world(trank_md)

            olap_mask(trank_world) = .true.
            group    (trank_world) = trank_cfd
            
        end do
        end do  
        end do


    end do
    
    call MPI_Barrier(CPL_REALM_COMM, ierr)

    ! Split world Comm into a set of comms for overlapping processors
    call MPI_comm_split(CPL_WORLD_COMM,group(rank_world),realm, &
                        CPL_OLAP_COMM,ierr)

    !Setup Overlap comm sizes and id
    if (realm.eq.cfd_realm) CFDid_olap = myid_olap
    call MPI_bcast(CFDid_olap,1,MPI_INTEGER,CFDid_olap,CPL_OLAP_COMM,ierr)

    ! USED ONLY FOR OUTPUT/TESTING??
    !if (myid_olap .eq. CFDid_olap) testval = group(rank_world)
    !call MPI_bcast(testval,1,MPI_INTEGER,CFDid_olap,CPL_OLAP_COMM,ierr)

    ! Set all non-overlapping processors to MPI_COMM_NULL
    if (olap_mask(rank_world).eqv..false.) then
        myid_olap = olap_null
        rank_olap = olap_null
        call MPI_COMM_FREE(CPL_OLAP_COMM, ierr)
    end if

    !Setup overlap map
    call CPL_rank_map(CPL_OLAP_COMM,rank_olap,nproc_olap, & 
                      rank_olap2rank_world,rank_world2rank_olap,ierr)
    myid_olap = rank_olap - 1

    deallocate(mdicoords)
    deallocate(mdjcoords)
    deallocate(mdkcoords)

    !Create a COMM for overlapping processes only
    if (olap_mask(rank_world).eqv..true.) then
        color = 1
    else
        color = 2
    endif
    call MPI_comm_split(CPL_WORLD_COMM, color,  myid_world, & 
                        CPL_REALM_INTERSECTION_COMM, ierr)

    call MPI_comm_rank(CPL_REALM_INTERSECTION_COMM, rank_intersect, ierr)
    rank_intersect = rank_intersect + 1

    !Create an intercomm for overlapping processes only
!    do n = 1,size(olap_mask,1)
!        if (olap_mask(n).eqv..true.) then
!            local_leader = n
!            exit
!        endif
!    enddo
!    ! Get the MPI_comm_world ranks that hold the largest ranks in cfd_comm and md_comm
!    call MPI_comm_size(CPL_REALM_COMM, comm_size, ierr)
!    ibuf(:) = -1
!    jbuf(:) = -1
!    if ( myid_realm .eq. comm_size - 1) then
!        ibuf(realm) = myid_world
!    endif

!    call MPI_allreduce( ibuf ,jbuf, 2, MPI_INTEGER, MPI_MAX, &
!                        CPL_WORLD_COMM, ierr)

!    !Set this largest rank on each process to be the inter-communicators (WHY NOT 0??)
!    select case (realm)
!    case (cfd_realm)
!        remote_leader = jbuf(md_realm)
!    case (md_realm)
!        remote_leader = jbuf(cfd_realm)
!    end select

!    print*, "intercomm create started", realm, local_leader, remote_leader,rank_world 
!    call MPI_intercomm_create(CPL_REALM_COMM, comm_size-1, CPL_OLAP_COMM, &
!                              remote_leader, 1, CPL_REALM_INTERSECTION_COMM, ierr)
!    print*, "intercomm create finished", realm, local_leader, remote_leader,rank_world 
!    
    !if (realm.eq.md_realm) call write_overlap_comms_md

end subroutine prepare_overlap_comms

!=========================================================================
!Setup topology graph of overlaps between CFD & MD processors

subroutine set_overlap_topology()
    use mpi
    implicit none

    integer                             :: i, n, nconnections
    integer, dimension(:),allocatable   :: index, edges
    logical                             :: reorder

    !Allow optimisations of ordering
    reorder = .true.

    !Get number of processors in communicating overlap region 
    if (olap_mask(rank_world).eqv..true.) then

        !CFD processor is root and has mapping to all MD processors
        allocate(index(nproc_olap))         ! Index for each processor
        allocate(edges(2*(nproc_olap)-1))   ! nproc_olap-1 for CFD and one for
                                            ! each of nproc_olap MD processors
        index = 0;  edges = 0

        !CFD processor has connections to nproc_olap MD processors
        nconnections = nproc_olap-1
        index(1) = nconnections
        do n = 1,nconnections
            edges(n) = n !olap_list(n+1) !CFD connected to all MD processors 1 to nconnections
        enddo

        !MD processor has a single connection to CFD
        nconnections = 1; i = 2
        do n = nproc_olap+1,2*(nproc_olap)-1
            index(i) = index(i-1) + nconnections !Each successive index incremented by one
            edges(n) = CFDid_olap   !Connected to CFD processor
            i = i + 1
        enddo

        !Create graph topology for overlap region
        call MPI_Graph_create(CPL_OLAP_COMM, nproc_olap, index, &
                              edges, reorder, CPL_GRAPH_COMM, ierr)
    else
        CPL_GRAPH_COMM = MPI_COMM_NULL
    endif

    ! Setup graph map
    call CPL_rank_map(CPL_GRAPH_COMM, rank_graph, nproc_olap, & 
                     rank_graph2rank_world, rank_world2rank_graph, ierr)
    myid_graph = rank_graph - 1

end subroutine set_overlap_topology


subroutine print_overlap_comms
    use mpi
    implicit none

    integer :: trank

    if (myid_world.eq.0) then
        write(7500+rank_realm,*) ''
        write(7500+rank_realm,*) '----------- OVERLAP COMMS INFO ------------'
        write(7500+rank_realm,*) '-------------------------------------------'
        write(7500+rank_realm,*) '        RANKS              BROADCAST TEST  '
        write(7500+rank_realm,*) '  world  realm  olap      testval( = group)'
        write(7500+rank_realm,*) '-------------------------------------------'
    end if
    
    do trank = 1,nproc_world
        if (rank_world.eq.trank) then
            write(7500+rank_realm,'(3i7,i16)') rank_world,rank_realm, &
                                rank_olap, testval  
        end if
    end do

    if (myid_world.eq.0) then
        write(7500+rank_realm,*) '-------- END OVERLAP COMMS INFO  ----------'
        write(7500+rank_realm,*) '==========================================='
    end if
    
end subroutine print_overlap_comms

end subroutine CPL_create_map

!-------------------------------------------------------------------
!                   CPL_rank_map                                   -
!-------------------------------------------------------------------

! Get COMM map for current communicator and relationship to 
! world rank used to link to others in the coupler hierachy

! - - - Synopsis - - -

! CPL_rank_map(COMM, rank, comm2world, world2comm, ierr)

! - - - Input Parameters - - -

!comm
!    communicator with cartesian structure (handle) 

! - - - Output Parameter - - -

!rank
!    rank of a process within group of comm (integer)
!    NOTE - fortran convention rank=1 to nproc  
!nproc
!    number of processes within group of comm (integer) 
!comm2world
!   Array of size nproc_world which for element at 
!   world_rank has local rank in COMM
!world2comm
!   Array of size nproc_COMM which for element at 
!   for local ranks in COMM has world rank 
!ierr
!    error flag


subroutine CPL_rank_map(COMM, rank, nproc, comm2world, world2comm, ierr)
    !use coupler_module, only : rank_world, nproc_world, CPL_WORLD_COMM, VOID
    use mpi
    implicit none

    integer, intent(in)                             :: COMM
    integer, intent(out)                            :: rank,nproc,ierr
    integer, dimension(:),allocatable,intent(out)   :: comm2world,world2comm

    allocate(world2comm( nproc_world))
    world2comm( nproc_world) = VOID

    if (COMM .ne. MPI_COMM_NULL) then

        !Mapping from comm rank to world rank
        call MPI_comm_rank(COMM,rank,ierr)
        rank = rank + 1
        call MPI_comm_size(COMM,nproc,ierr)
        allocate(comm2world(nproc))
        call MPI_allgather(rank_world,1,MPI_INTEGER, & 
                           comm2world,1,MPI_INTEGER,COMM,ierr)
    else
        rank = VOID
        allocate(comm2world(0))
    endif

    !Mapping from world rank to comm rank
    call MPI_allgather(rank      ,1,MPI_INTEGER, & 
                       world2comm,1,MPI_INTEGER,CPL_WORLD_COMM,ierr)

end subroutine CPL_rank_map

!---------------------------------------------------
! Locate file in input

subroutine locate(fileid, keyword, have_data)
    implicit none
    
    integer,intent(in)          :: fileid               ! File unit number
    character(len=*),intent(in) :: keyword              ! Input keyword 
    logical,intent(out)         :: have_data            ! Flag: input found

    character*(100)             :: linestring           ! First 100 chars
    integer                     :: keyword_length       ! Length of keyword
    integer                     :: io                   ! File status flag

    keyword_length = len(keyword)
    rewind(fileid)
    
    ! Loop until end of file or keyword found
    do
        ! Read first 100 characters of line
        read (fileid,'(a)',iostat=io) linestring

        ! If end of file is reached, exit
        if (io.ne.0) then 
            have_data = .false.
            exit
        end if
        
        ! If the first characters match keyword, exit
        if (linestring(1:keyword_length).eq.keyword) then
            have_data = .true.
            exit
        endif

    end do

end subroutine locate

!===========================================================================
!Error handling subroutines

subroutine error_abort_s(msg)
    use mpi
    use iso_fortran_env, only : error_unit
    implicit none

    character(len=*), intent(in), optional :: msg
   
    integer errcode,ierr
    
    ! In Unix systems 0 means exit with success
    errcode = 1

    if (present(msg)) then 
        write(error_unit ,*) msg
        flush(error_unit)
        call sleep(2)
    endif

    call MPI_Abort(CPL_WORLD_COMM,errcode,ierr)

end subroutine error_abort_s


subroutine error_abort_si(msg,i)
    use mpi
    use iso_fortran_env, only : error_unit
    implicit none

    character(len=*), intent(in) :: msg
    integer, intent(in) :: i

    integer errcode,ierr

    ! In Unix systems 0 means exit with success
    errcode = 1

    write(error_unit,*) msg,i
    flush(error_unit)
    call sleep(2)

    call MPI_Abort(CPL_WORLD_COMM,errcode,ierr)

end subroutine error_abort_si


subroutine messenger_lasterrorcheck
    use mpi
    use iso_fortran_env, only : error_unit
    implicit none

    integer resultlen
    character*12 err_buffer

    call MPI_Error_string(ierr,err_buffer,resultlen,ierr)
    print*, err_buffer

end subroutine messenger_lasterrorcheck


!--------------------------------------------------------------------------------------
! Prints formatted debug statements
subroutine printf(buf,dplaces_in)
    implicit none

    real(kind(0.d0)),dimension(:),intent(in):: buf
    integer, intent(in), optional           :: dplaces_in

    integer             :: n,dplaces, space
    real(kind(0.d0))    :: maxbuf,minbuf,order
    character*19        :: string
    character*42        :: buf_precision

    space = 2

    !Default number of decimal places if not supplied
    if (present(dplaces_in)) then
        if (dplaces_in .le. 9) then
            dplaces = dplaces_in
        else
            print*, 'Number of decimal places in printf if limited to 9'
            dplaces = 9 !Maximum
        endif
    else
        dplaces = 4
    endif

    !Find out required format to display maximum element in buffer
    maxbuf = maxval(buf); minbuf = minval(buf)
    maxbuf = max(maxbuf,10*abs(minbuf)) !10*Ensures extra space for minus sign
    order = 1.d0; n = 1
    do while (max(maxbuf,order) .ne. order)
        order = order*10.d0
        n = n + 1
    enddo
    if (maxbuf .lt. 0.d0 .and. maxbuf .gt. -1.d0) then
        n = n + 1 !For the case of -0.something
    endif
    if (n+dplaces+space .le. 9) then
        write(buf_precision,'(a,i1,a,i1)') 'f',n+dplaces+space,'.', dplaces
    else
        write(buf_precision,'(a,i2,a,i1)') 'f',n+dplaces+space,'.', dplaces
    endif

    ! Build up format specifier string based on size of passed array
    string='(a6,i3,   ' // trim(buf_precision) // ')'
    write(string(8:10),'(i3)') size(buf) 

    !Write formatted data 
    print(string),'printf',rank_world,buf

end subroutine printf


!--------------------------------------------------------------------------------------
!Write matrix in correct format

subroutine write_matrix_int(a,varname,fh)
    implicit none

    integer                 :: i,j,fh
    character(*)            :: varname
    integer, dimension(:,:) :: a

    write(fh,*) varname
    do i = lbound(a,1), ubound(a,1)
        write(fh,*) (a(i,j), j = lbound(a,2), ubound(a,2))
    end do

end subroutine write_matrix_int

subroutine write_matrix(a,varname,fh)
    implicit none

    integer                          :: i,j,fh
    character(*)                     :: varname
    real(kind(0.d0)), dimension(:,:) :: a

    write(fh,*) varname
    do i = lbound(a,1), ubound(a,1)
        write(fh,*) (a(i,j), j = lbound(a,2), ubound(a,2))
    end do

end subroutine write_matrix

!===========================================================================
! Subroutine that can be used to stop the code when reaching a given 
! point in coupler -- useful when coupling new codes
!---------------------------------------------------------------------------
!subroutine request_stop(tag)
!    use mpi
!    implicit none
!
!    character(len=*),intent(in) ::tag
!    integer myid, ierr
!
!    ! do nothing, get out quick 
!    if(.not. stop_request_activated ) return
!
!    if (tag /= stop_request_name) return
!
!    select case(stop_request_name)
!    case("create_comm","CREATE_COMM")
!        call mpi_comm_rank(CPL_REALM_COMM, myid,ierr)
!        write(0,*) 'stop as requested at ', trim(stop_request_name), ', realm',realm, 'rank', myid
!        call MPI_Finalize(ierr)
!        stop
!    case("create_map","CREATE_MAP")
!        call mpi_comm_rank(CPL_REALM_COMM, myid,ierr)
!        write(0,*) 'stop as requested at ', trim(stop_request_name), ', realm',realm, 'rank', myid
!        call MPI_Finalize(ierr)
!        stop    
!    case default
!        write(0,*) "WARNING: request abort activated, but the tag is unrecognized, check COUPLER.in"
!        write(0,*) "         accepted stop tags are: create_comm"
!    end select
!
!end subroutine request_stop

function CPL_new_fileunit() result (f)
    implicit none

    logical :: op
    integer :: f

    f = 1
    do 
        inquire(f,opened=op)
        if (op .eqv. .false.) exit
        f = f + 1
    enddo

end function


subroutine set_output_mode(mode)
    implicit none

    integer, intent(in) :: mode
    
    output_mode = mode
end subroutine

!--------------------------------------------------------------------------------------
! Return unused fileunit by checking all exisiting
function get_new_fileunit() result (f)
    implicit none

    logical    :: op
    integer    :: f

    f = 1
    do 
        inquire(f,opened=op)
        if (op .eqv. .false.) exit
        f = f + 1
    enddo

end function


end module coupler_module




!!=============================================================================
!!   Adjust CFD domain size to an integer number of lattice units used by  
!!   MD if sizes are given in sigma units
!!-----------------------------------------------------------------------------

!subroutine CPL_cfd_adjust_domain(xL, yL, zL, nx, ny, nz, density_output)
!    use mpi
!    !use coupler_module, only : density_cfd,CPL_REALM_COMM, rank_realm, ierr
!    implicit none

!    integer, optional, intent(inout)            :: nx, ny, nz
!    real(kind(0.d0)), optional, intent(inout)   :: xL,yL,zL
!    real(kind(0.d0)), optional, intent(inout)   :: density_output

!    ! Internal variables
!    integer                                     :: ierror, root
!    !character(1)                                :: direction
!    logical                                     :: changed

!    !Define root processes
!    root = 1

!    density_output = density_cfd

!    ! Check CFD domain and MD domain are compatible sizes to allow a
!    ! stable initial MD lattice structure - resize if possible or
!    ! stop code and demand a regeneration of grid if vary by more than 0.01
!    changed = .false.
!    if (present(xL)) then
!        call init_length(xL,resize=.true.,direction='x', &
!                         print_warning=changed)
!    endif

!    ! No need to adjust y because we can adjust DY in MD to
!    ! have an integer number of FCC units. ??????? What

!    if (present(zL)) then
!        call init_length(zL,resize=.true.,direction='z', &
!                         print_warning=changed)
!    endif

!    if ( changed ) then
!        print*, "CPL_cfd_adjust_domain error - Regenerate Grid with corrected sizes as above"
!        call MPI_Abort(CPL_WORLD_COMM,ierror,ierr)
!    endif

!    ! check id CFD cell sizes are larger than 2*sigma 
!    call test_cfd_cell_sizes

!contains

!!-----------------------------------------------------------------------------

!subroutine init_length(rout,resize,direction,print_warning)
!    !use coupler_module, only: dx,dy,dz,error_abort
!    implicit none
!            
!    real(kind=kind(0.d0)), intent(inout) :: rout
!    logical, intent(in)                  :: resize
!    character(*),intent(in)              :: direction
!    logical,intent(out)                  :: print_warning

!    real(kind(0.d0)) :: dxyz  ! dx, dy or dz
!    real(kind(0.d0)) :: rinit ! initial val of rout or rin for print

!    print_warning=.false.

!    select case (direction)
!    case('x','X')
!        dxyz = dx
!    case('y','Y')
!        dxyz = dy
!    case('z','Z')
!        dxyz = dz
!    case default
!        call error_abort('Wrong direction specified in init_length')
!    end select

!    if ( resize ) then

!        rinit = rout
!        rout = real(nint(rout/dxyz),kind(0.d0))*dxyz
!        print_warning = .true.
!        print*, direction, 'dxyz = ', dxyz 

!    endif

!    if (print_warning) then 

!        !if (rank_realm .eq. root) then 

!            write(*,'(3(a,/),3a,/,2(a,f20.10),/a,/,a)') &
!                    "*********************************************************************",    &
!                    "WARNING - this is a coupled run which resets CFD domain size         ",    &
!                    " to an integer number of MD initial cells:                           ",    &
!                    "   Domain resized in the ", direction, " direction                   ",    &
!                    " inital size =", rinit, " resized ", rout,                                 &
!                    "                                                                     ",    & 
!                    "*********************************************************************"   
!        !endif

!        !If resize is insignificant then return flag print_warning as false
!        if (abs(rinit-rout) .lt. 0.01) print_warning = .false.

!    end if
!    
!end subroutine init_length

!!-----------------------------------------------------------------------------

!subroutine test_cfd_cell_sizes
!    implicit none

!    if (rank_realm .eq. root) then
!        if (present(xL) .and. present(nx)) then
!            if (xL/nx < 2.0d0) then
!                write(0,*)" WARNING: CFD cell size in x direction is less that 2 * sigma. Does this make sense?" 
!                write(0,*)"          xL=",xL,"nx=",nx
!            endif
!        endif

!        if (present(yL) .and. present(ny)) then
!            if (yL/ny < 2.0d0) then
!                write(0,*)" WARNING: CFD cell size in y direction is less that 2 * sigma. Does this make sense?" 
!                write(0,*)"          yL=",yL,"nx=",ny
!            endif
!        endif

!        if (present(zL) .and. present(nz)) then
!            if (zL/nz < 2.0d0) then
!                write(0,*)" WARNING: CFD cell size in z direction is less that 2 * sigma. Does this make sense?" 
!                write(0,*)"          zL=",zL,"nx=",nz
!            endif
!        endif
!    end if
!        
!end subroutine test_cfd_cell_sizes
!            
!end subroutine CPL_cfd_adjust_domain













