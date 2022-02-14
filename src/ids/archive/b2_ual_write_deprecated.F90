!!-----------------------------------------------------------------------------
!! DOCUMENTATION (doxygen 1.8.8):
!>      @author
!>      Dejan Penko
!!
!>      @page b2uw_dep b2_ual_write_deprecated
!>      @section b2uw_dep_desc  Description
!!      @note   This script is the \b old b2_ual_write.F90 and was kept as an
!!              archive for possible future references and as an example for
!!              information purposes. The development of the old file continued
!!              under new codes named
!!              @ref b2uw_gsl "b2_ual_write_gsl",
!!              @ref b2uw_b2mod "b2_ual_write_b2mod" and the latest \b new
!!              @ref b2uw_b2.
!!              The development 'moved' to new codes due to different
!!              approach of writing data to IDS using IMAS GGD routines and
!!              b2mod scripts.
!!      @note   More on the b2_ual_writers is available in SOLPS-GUI
!!              documentation \b HOWTOs under section <b> 4.5 IMAS </b>.
!!
!!      b2_ual_write_deprecated.f90 code is used to generate
!!      b2_ual_write_deprecated.exe (main
!!      program), which is a post-processor for b2. It reads the plasma
!!      grid geometry (including Nodes and Cells grid subsets) and plasma
!!      state (electron density, electron temperature, ion temperature) and
!!      writes it to IDS database.
!!
!!      \b References:
!!          - @ref b2uw_dep_prog "b2_ual_write_deprecated file reference"
!!          - @ref b2uw_gsl "b2_ual_write_gsl"
!!          - @ref b2uw_b2mod "b2_ual_write_b2mod"
!!
!!-----------------------------------------------------------------------------

!!-----------------------------------------------------------------------------
!>      @section b2uw_dep_prog Program b2_ual_write_deprecated
!!      References:
!!          - @ref b2uw_dep "b2_ual_write_deprecated main page"
!!
!!      @subsection b2uw_dep_det  Details
!!      For more information see also routine b2cdca.
!!
!!      The complete program performs post-processing of the
!!      result of a b2 calculation.
!!      This program unit opens and closes the input/output units, and
!!      may perform some other system-dependent operations.
!!
!!      The input units are:
!!          - ninp(0): formatted, provides output control parameters.
!!          - ninp(1): un*formatted, provides the geometry.
!!          - ninp(2): un*formatted, provides the run parameters.
!!          - ninp(3): un*formatted, provides the plasma state.
!!          - ninp(4): unformatted, provides the detailed plasma state.
!!          - ninp(5): formatted, provides the run switches.
!!          - ninp(6): un*formatted, provides the atomic data.
!!
!!      The output units are:
!!          - nout(0): formatted, provides printed output.
!!
!!      @note   See routine b2cdca for the meaning of 'un*formatted'.
!!
!!      @subsection b2uw_pv  Parameters/variables
!!      @note   see also routine /b b2cdcv.F
!!
!!      @param  device - Device name of the IMAS IDS database
!!              (i. e. solps-iter, iter, aug)
!!      @param  edge_profiles - IDS designed to store data on edge plasma
!!              profiles  (includes the scrape-off layer and possibly part
!!              of the confined plasma)
!!      @param  idum
!!              <i> ( (0:9) integer array) </i>
!!      @param  idx - The returned identifier to be used in the subsequent
!!                    data access operation
!!      @param  io - iostat value
!!      @param  ne1 - Data field holding values in relation to quantity:
!!              Electron density
!!              <i> (real(kind=B2R8) array ) </i>
!!      @param  ninp - Specifies the input unit numbers
!!              <i> ( (0:6) integer array) </i>
!!      @param  nout - Specifies the output unit numbers
!!              <i> ( (0:2) integer array) </i>
!!      @param  ns - Specifies the number of atomic species in the calculation.
!!              The species are indexed by (0:ns-1). It will hold that 1.le.ns.
!!              <i> (integer) </i>
!!      @param  nx, ny - Specifies the number of interior cells along the
!!              first and the second coordinate, respectively.
!!              The total number of cells is (nx+2)*(ny+2), they are indexed
!!              by (-1:nx,-1:ny). It will hold that 0.le.nx and 0.le.ny
!!              <i> (integer) </i>
!!      @param  run - The run number of the database being created
!!      @param  shot - The shot number of the database being created
!!      @param  te1 - Data field holding values in relation to quantity:
!!              Electron temperature
!!              <i> (real(kind=B2R8) array ) </i>
!!      @param  ti1 - Data field holding values in relation to quantity:
!!              Ion temperature
!!              <i> (real(kind=B2R8) array ) </i>
!!      @param  treename - the name of the IMAS IDS database,
!!              (i.e. "edge_profiles" (mandatory) )
!!      @param  username - Creator/owner of the IMAS IDS database
!!      @param  version - Major version of the IMAS IDS database
!!
!!      @subsection b2uw_eind Error indicators
!!      In case an error condition is detected, a call is made to the
!!      routine \b xerrab. This causes an error message to be printed,
!!      after which the program halts.
!!
!!      @subsection b2uw_syx    Exceptional syntax explanation
!!      @code
!!          ! IGNORE    !! syntax used to ignore this module in list
!!                      !! dependency when compiling the code
!!      @endcode
program b2_ual_write_deprecated
    use b2mod_types , B2R8 => R8
    use b2mod_rates
    use b2mod_ual_io

    use ids_schemas     ! IGNORE
                        !! These are the Fortran type definitions for the
                        !! Physics Data Model
    use ids_routines    ! IGNORE
                        !! These are the Access Layer routines + management of
                        !! IDS structures

    implicit none

    !! Internal variables
    integer :: ninp(0:6)    !< Specifies the input unit numbers
    integer :: nout(0:2)    !< Specifies the output unit numbers
    integer :: nx   !< Specifies the number of interior cells along the \b first
        !< coordinate, respectively. The total number of cells is
        !< (nx+2)*(ny+2), they are indexed by (-1:nx,-1:ny). It will hold
        !< that 0.le.nx and 0.le.ny
    integer :: ny   !< Specifies the number of interior cells along the \b
        !< second coordinate, respectively. The total number of cells is
        !< (nx+2)*(ny+2), they are indexed by (-1:nx,-1:ny). It will hold
        !< that 0.le.nx and 0.le.ny
    integer :: ns   !< Specifies the number of atomic species in the
        !< calculation. The species are indexed by (0:ns-1). It will hold
        !< that 1.le.ns.
    integer :: idum(0:9)
    integer :: io   !< iostat value
    real(kind=B2R8), allocatable :: ne1(:)  !< Data field holding values in
        !< relation to quantity: Electron density
    real(kind=B2R8), allocatable :: te1(:)  !< Data field holding values in
        !< relation to quantity: Electron temperature
    real(kind=B2R8), allocatable :: ti1(:)  !< Data field holding values in
        !< relation to quantity: Ion temperature
    integer :: idx  !< The returned identifier to be used in the subsequent
        !< data access operation
    integer :: shot !< The shot number of the database being created
    integer :: run  !< The run number of the database being created
    character(len=24) :: treename !< The name of the IMAS IDS database
        !< (i.e. "edge_profiles" (mandatory) )
    character(len=24) :: username   !< Creator/owner of the IMAS IDS database
    character(len=24) :: device !< Device name of the IMAS IDS database
        !< (i. e. solps-iter, iter, aug)
    character(len=24) :: version    !< Major version of the IMAS IDS database
    type(ids_edge_profiles) :: edge_profiles !< IDS designed to store data on
        !< edge plasma profiles  (includes the scrape-off layer and possibly
        !< part of the confined plasma)
    !! Procedures
    external prgini, prgend, xerset, xertst, cfopen, cfruin
    !! Initialize input/output units
    data ninp/50, 51, 52, 53, 54, 55, 56/
    data nout/60, 61, 62/

    !! Read grid geometry data from b2fgmtry file
    call read_b2fgmtry()

    !! Read quantity data from b2fstate file
    call read_b2fstate( ninp, nout, nx, ny, ns, ne1, te1, ti1 )

    treename    = "ids"
    shot        = 16151
    run         = 1001
    username    = "penkod"
    device     = "solps-iter"
    version     = "3"

    !! Set the data to edge_profiles IDS
    call write_ids_edge_profiles( treename, shot, run, idx, username,   &
        &   device, version, ne1, te1, ti1 )

    ! call read_ids(treename, shot, run, idx, username, &
    !                                     & device, version)

contains

    !> Subroutine used to read \b b2fgmrty (file holding grid geometry data).
    subroutine read_b2fgmtry()
        !! Internal variables
        character(len=120) :: lblgm
        integer :: i, j, k

        !! Program start_up calls
        call prgini( 'b2_ual_write_deprecated' )

        !! Open files
        !! Input files
        ! call cfopen( ninp(0), 'b2_ual_write.dat', 'old', 'formatted' )
        call cfopen( ninp(1), 'b2fgmtry', 'old', 'un*formatted' )
        ! call cfopen( ninp(2), 'b2fparam', 'old', 'un*formatted' )
        call cfopen( ninp(3), 'b2fstate', 'old', 'un*formatted' )
        ! call cfopen( ninp(4), 'b2fplasma', 'old', 'unformatted' )
        ! call cfopen( ninp(5), 'b2mn.dat', 'old', 'formatted' )
        ! call cfopen( ninp(6), 'b2frates', 'old', 'formatted' )
        !!...output files
        call cfopen( nout(0), 'b2_ual_write_deprecated.prt', 'new', 'formatted' )
        call cfopen( nout(1), 'b2_ual_writ2.prt', 'new', 'formatted' )

        !!..mark output unit for error messages
        call xerset(0)
        !!..obtain version numbers
        call cfverr( ninp(1), b2fgmtry_version )
        ! call cfverr ( ninp(2), b2fparam_version )
        call cfverr( ninp(3), b2fstate_version )
        ! call cfverr ( ninp(6), b2frates_version )
        !!..obtain nx, ny, ns
        call cfruin( ninp(3), 3, idum, 'nx,ny,ns' )
        nx = idum(0)
        ny = idum(1)
        ns = idum(2)
        call xertst( 0.le.nx .and. 0.le.ny .and. 1.le.ns,   &
            &   'faulty input nx, ny, ns from plasma state file' )
        call cfruin ( ninp(1), 2, idum(0), 'nx,ny' )
        call xertst( idum(0).eq.nx .and. idum(1).eq.ny, &
            &   'faulty input nx, ny from geometry file' )

        call alloc_b2mod_geo( nx, ny )
        call alloc_b2mod_plasma( nx, ny, ns )
        call alloc_b2mod_indirect( nx, ny, nncutmax )
        call alloc_b2mod_sources( nx, ny, ns )

        call cfruch( ninp(1), 120, lblgm, 'label' )
        call cfruin( ninp(1), 1, idum, 'isymm' )
        call b2rugm( ninp(1), nx, ny, crx, cry, fpsi, ffbz, &
            &   bb, vol, hx, hy, qz, qc, qcb, gs, pbs, &
            &   wbbl, wbbr, wbbv, wbbc, cell_width, cell_height, gmap )

        call xertst(size(crx).eq.size(cry), &
            &   "The number of x and y coordinates is not the same." )

        corner: do k = 0, 3, 1
            xi: do j = -1, ny, 1
                yi: do i = -1, nx, 1
                    write(nout(0),*) i, j, k, crx(i, j, k)
                end do yi
            end do xi
        end do corner

        call prgend()

    end subroutine read_b2fgmtry

    !> Subroutine used to read \b b2fstate (file holding quantity data fields).
    subroutine read_b2fstate( ninp, nout, nx, ny, ns, ne, te, ti )
        integer :: ninp(0:6)    !< Specifies the input unit numbers
        integer :: nout(0:2)    !< Specifies the output unit numbers
        integer :: nx   !< Specifies the number of interior cells along the
                        !< first coordinate
        integer :: ny   !< Specifies the number of interior cells along the
                        !< second coordinate
        integer :: ns   !< Specifies the number of atomic species in the
                        !< calculation
        integer :: idum(0:9)   !! nscx, iscx(0:nscxmax-1), ismain
        real(kind=B2R8), intent(inout), allocatable :: ne(:)    !< Data field
            !< holding values in relation to quantity: Electron density
        real(kind=B2R8), intent(inout), allocatable :: te(:)    !< Data field
            !< holding values in relation to quantity: Electron temperature
        real(kind=B2R8), intent(inout), allocatable :: ti(:)    !< Data field
            !< holding values in relation to quantity: ion temperature
        !! Internal variables
        character :: lblgm*120
        character :: lblcp*120
        character :: lblmn*120
        character :: lblrc*120
        character :: exp*128
        character :: comment*128
        integer :: shot
        integer :: overwrite_shotnumber
        integer :: cf_sizes(22)
        logical :: timedep
        logical :: snapshot
        logical :: tallies
        logical :: movies
        !! Full description of most of the variables below is available
        !! in b2cdcv.F
        real (kind=B2R8), allocatable :: zamin(:)   !< Minimal atomic charge
                                                    !< range pf specie (is)
        real (kind=B2R8), allocatable :: zamax(:)   !< Maximum atomic charge
                                                    !< range pf specie (is)
        real (kind=B2R8), allocatable :: zn(:) !< Nuclear charge of specie (is)
        real (kind=B2R8), allocatable :: am(:) !< Atomic mass of specie (is)
        real (kind=B2R8), allocatable :: na(:) !< Ion density [particles/mˆ3].
        real (kind=B2R8), allocatable :: ua(:) !< Parallel velocity [m/s].
        real (kind=B2R8), allocatable :: uadia(:)   !< The total effective
            !< drift velocity of specie
        real (kind=B2R8), allocatable :: po(:)  !< Electric potential
        real (kind=B2R8), allocatable :: fna(:) !< Flux of atoms
        real (kind=B2R8), allocatable :: fhe(:) !< Electron heat flux
        real (kind=B2R8), allocatable :: fhi(:) !< All atom heat flux
        real (kind=B2R8), allocatable :: fch(:) !< Electric current
        real (kind=B2R8), allocatable :: fch_32(:)
        real (kind=B2R8), allocatable :: fch_52(:)
        real (kind=B2R8), allocatable :: kinrgy(:)
        real (kind=B2R8), allocatable :: fch_p(:) !< The product of the
            !< Parallel electric current and the poloidal magnetic field
            !< component.
        real(kind=B2R8) :: time
        !! Procedures
        external subini, subend, xertst, xerrab, cfruch, cfruin, cfrure
        !! Namelist
        namelist /b2md_namelist/ exp, shot, time, comment, timedep,     &
            &   snapshot, tallies, movies, overwrite_shotnumber

        !! Subprogram start-up calls
        call subini( 'b2mddr' )
        !! Test ninp, nout

        call xertst( 1 .le.ninp(0) .and. 1 .le.ninp(1) .and.                  &
            &   1 .le.ninp(2) .and. 1 .le.ninp(3) .and. 1 .le.ninp(5) .and.    &
            &   1 .le.ninp(6) .and. 1 .le.nout(0) .and. 1 .le.nout(1) .and.    &
            &   1 .le.nout(2), 'faulty argument ninp, nout' )
        !! Test dimensions

        call xertst( 0.le.nx .and. 0.le.ny, 'faulty argument nx, ny' )
        call xertst( 1.le.ns, 'faulty argument ns' )

        call xertst( ns .le.nsdecl, 'faulty parameter nsdecl' )

        !! Get array sizes from the b2fstate. We open the file again as a new
        !! unit to not disrupt the old unit.
        !! A routine should exist to do all that work but so far with no success
        !! of finding/properly use it.
        cf_sizes = read_additional_sizes(22)

        ! allocate( lblgm(cf_sizes(2) ) )
        allocate( zamin( cf_sizes(3) ) )
        allocate( zamax( cf_sizes(4) ) )
        allocate( zn( cf_sizes(5) ) )
        allocate( am( cf_sizes(6) ) )
        allocate( na( cf_sizes(7) ) )
        allocate( ne( cf_sizes(8) ) )
        allocate( ua( cf_sizes(9) ) )
        allocate( uadia( cf_sizes(10) ) )
        allocate( te( cf_sizes(11) ) )
        allocate( ti( cf_sizes(12) ) )
        allocate( po( cf_sizes(13) ) )
        allocate( fna( cf_sizes(14) ) )
        allocate( fhe( cf_sizes(15) ) )
        allocate( fhi( cf_sizes(16) ) )
        allocate( fch( cf_sizes(17) ) )
        allocate( fch_32( cf_sizes(18) ) )
        allocate( fch_52( cf_sizes(19) ) )
        !! The program gives an error saying the kinrgy is already allocated.
        !! Where, when?
        ! allocate( kinrgy( cf_sizes(20) ) )
        ! allocate( fch_p( cf_sizes(22) ) )

        ! call cfruch( ninp(3), cf_sizes(2), lblgm, 'label' )
        call cfruch( ninp(3), 120, lblgm, 'label' )
        call cfrure( ninp(3), cf_sizes(3), zamin, 'zamin' )
        call cfrure( ninp(3), cf_sizes(4), zamax, 'zamax' )
        call cfrure( ninp(3), cf_sizes(5), zn, 'zn' )
        call cfrure( ninp(3), cf_sizes(6), am, 'am' )
        call cfrure( ninp(3), cf_sizes(7), na, 'na' )
        call cfrure( ninp(3), cf_sizes(8), ne, 'ne' )
        call cfrure( ninp(3), cf_sizes(9), ua,'ua' )
        call cfrure( ninp(3), cf_sizes(10), uadia,'uadia' )
        call cfrure( ninp(3), cf_sizes(11), te, 'te' )
        call cfrure( ninp(3), cf_sizes(12), ti, 'ti' )
        call cfrure( ninp(3), cf_sizes(13), po, 'po' )
        call cfrure( ninp(3), cf_sizes(14), fna, 'fna' )
        call cfrure( ninp(3), cf_sizes(15), fhe, 'fhe' )
        call cfrure( ninp(3), cf_sizes(16), fhi, 'fhi' )
        call cfrure( ninp(3), cf_sizes(17), fch,'fch' )
        call cfrure( ninp(3), cf_sizes(18), fch_32, 'fch_32' )
        call cfrure( ninp(3), cf_sizes(19), fch_52, 'fch_52' )
        ! call cfrure( ninp(3), cf_sizes(20), kinrgy, 'kinrgy' )
        ! call cfrure( ninp(3), 1, time, 'time' )
        ! call cfrure( ninp(3), cf_sizes(22), fch_p, 'fch_p' )

    end subroutine read_b2fstate

    !> Get array sizes from the b2fstate.
    function read_additional_sizes( array_size )
        integer, intent(in) :: array_size   !< Size of the array
        character(len=256) :: rdline
        character(len=256) :: rdstring
        integer :: cfcount, stat
        integer :: cf_sizes(array_size)
        integer, dimension(:), allocatable :: read_additional_sizes

        allocate( read_additional_sizes( array_size ) )

        cfcount = 0
        open( 1, file='b2fstate', status="old" )
        do
            read(1,'(A)', iostat=io) rdline
            if( io/=0 ) exit
            if( rdline(1:3) .eq. "*cf") then
                cfcount = cfcount + 1
                rdstring = adjustl( rdline(20:28) )
                read( rdstring, *, iostat=stat )  cf_sizes( cfcount )
            end if
        end do
        close(1)
        read_additional_sizes = cf_sizes

    end function read_additional_sizes

    !> Subroutine used to put data to edge_profiles IDS
    subroutine write_ids_edge_profiles( treename, shot, run, idx, username, &
            &   device, version, ne, te, ti )
        !! Internal variables
        character(len=24) :: treename   !< The name of the IMAS IDS database
            !< (i.e. "edge_profiles" (mandatory) )
        character(len=24) :: username   !< Creator/owner of the IMAS IDS database
        character(len=24) :: device !< Device name of the IMAS IDS database
            !< (i. e. solps-iter, iter, aug)
        character(len=24) :: version    !< Major version of the IMAS IDS database
        character(len=255) :: grid_description  !< Grid description string
        character(len=132) :: gridSubset_name   !< Grid subset label
        integer :: shot !< The shot number of the database being created
        integer :: run  !< The run number of the database being created
        integer :: idx  !< The returned identifier to be used in the subsequent
            !< data access operation
        real (kind=B2R8), intent(in) :: ne(:)   !< Data field holding values in
            !< relation to quantity: Electron density
        real (kind=B2R8), intent(in) :: te(:)   !< Data field holding values in
            !< relation to quantity: Electron temperature
        real (kind=B2R8), intent(in) :: ti(:)   !< Data field holding values in
            !< relation to quantity: Ion temperature
        !! Internal variables
        integer :: i, j, k, icount, n
        integer :: num_nodes_all
        integer :: num_nodes
        integer :: num_gridSubsets
        integer :: gridSubset_index
        integer :: gridSubset_dim_index
        integer :: numCellsX
        integer :: numCellsY
        integer :: num_cells
        integer :: cellId
        integer :: num_ne_gridSubset
        integer :: num_ne_values
        integer :: num_te_gridSubset
        integer :: num_te_values
        integer :: num_ti_gridSubset
        integer :: num_ti_values
        integer :: num_ti_species
        integer :: ion_specie
        integer :: num_obj_0D
        integer :: obj_0D_id
        integer :: num_obj_2D
        integer :: obj_2D_id
        real(B2R8) :: time  !< Generic time
        type(ids_edge_profiles) :: edge_profiles    !< IDS designed to store
            !< data on edge plasma profiles  (includes the scrape-off layer and
            !< possibly part of the confined plasma)
        type(ids_edge_profiles_time_slice), pointer :: ggd  !< Type of IDS
            !< data structure, designed to store edge plasma quantities
            !< represented using the general grid description, for various
            !< time slices
        type(ids_generic_grid_dynamic), pointer :: grid !< Type of IDS
            !< data structure, designed for handling grid geometry data
        type(ids_generic_grid_dynamic_space), pointer :: space  !< Type of IDS
            !< data structure, designed for handling set of grid space data

        !! ===  SET UP IDS ===
        write(0,*) "IDS parameters"
        write(0,*) "shot: ", shot
        write(0,*) "run:", run
        write(0,*) "username: ", username
        write(0,*) "device: ", device
        write(0,*) "version: ", version

        !! We already went through the check if size(crx)==size(cry)
        !! so we can safely assume that num_nodes == size(crx)
        num_nodes_all   = size(crx) !! Number of all available coordinates

        write(0,*) "Setting data for edge_profiles IDS"

        !! Prepeare database
        !! In order to write to IDS database there are next steps that are
        !! absolutely neccessary to do, otherwise writing to IDS database
        !! will surely fail
        !! 1. Allocate profiles_1d
        allocate( edge_profiles%profiles_1d(1) )
        !! 2. Set homogeneous_time
        edge_profiles%ids_properties%homogeneous_time = 1
        !! 3. Allocate edge_profiles.time and Set desired value
        allocate( edge_profiles%time(1) )
        edge_profiles%time(1) = 1
        ! edge_profiles%time(1) = time

        !! Allocate ggd slice
        allocate( edge_profiles%ggd(1) )
        !! Set pointers to top IDS nodes/structures
        ggd => edge_profiles%ggd(1)
        grid => ggd%grid

        !! Allocate space
        allocate( grid%space(1) )
        !! Allocate objects_per_dimension substructure
        allocate( grid%space(1)%objects_per_dimension(3) )
        !! Set IDS grid description
        grid_description = "This IDS was written by b2_ual_write_deprecated"
        grid%identifier%description = grid_description

        !! === SET UP GRID ===
        !! The 2D structured grid is in our case composed out of one 2D
        !! structured space
        !! Set definition of the coordinate system of the space

        !! Allocate coordinates_type
        allocate( grid%space(1)%coordinates_type(2) )
        !! Set coordinate types
        grid%space(1)%coordinates_type(1) = 4 !! R
        grid%space(1)%coordinates_type(2) = 5 !! Z

        !! --- Set up grid space objects for Class 1 objects - points ---

        num_nodes = num_nodes_all
        num_obj_0D = num_nodes

        !! We are in 2D dimension space and the nodes are defined by coordinates
        !! in a way
        !! n1=[r1, z1], n2=[r2,z2], ..., nn=[rn, zn].

        allocate( grid%space(1)%objects_per_dimension(1)%object(num_obj_0D) )

        !! Set geometry (R,Z) coordinate of each node object
        do i = 1, num_obj_0D
            !! Set and fill list of indices for nodes
            allocate( grid%space(1)%objects_per_dimension(1)%object(i)%nodes(1) )
            grid%space(1)%objects_per_dimension(1)%object(i)%nodes(1) = i
            !! Set geometry (R,Z) coordinate of each node object
            allocate( grid%space(1)%objects_per_dimension(1)%object(i)%geometry(2) )
        enddo

        !! Get nodes coordinates
        icount = 0
        do k = 0, 3
            do j = -1, ny
                do i = -1, nx
                    icount = icount + 1
                    !! Fill geometry (R,Z) coordinate of each node object
                    grid%space(1)%objects_per_dimension(1)%     &
                        &   object(icount)%geometry(1) = crx(i,j,k)
                    grid%space(1)%objects_per_dimension(1)%     &
                        &   object(icount)%geometry(2) = cry(i,j,k)
                enddo
            enddo
        enddo

        !! --- (TODO) Set up connectivity array of grid space objects for   ---
        !! --- Class 2 objects - edges                                      ---
        !! TODO! Currently only placeholder edges are given

        !! Set (placeholder) list of indices for nodes defining each edge object
        allocate( grid%space(1)%objects_per_dimension(2)%object(1) )
        allocate( grid%space(1)%objects_per_dimension(2)%   &
            &   object(1)%nodes(1) )
        grid%space(1)%objects_per_dimension(2)%     &
            &   object(1)%nodes(1) = 0

        !! --- Set up connectivity array of grid space objects for      ---
        !! --- Class 3 objects - 2D cells                               ---

        !! Set list of indices for nodes defining each cell object
        !!
        !! We are are in 2D dimension and each cell is defined by 4 nodes
        !! Note: cell nodes must be ordered cyclically (here anti-clockwise)

        numCellsX = nx + 2
        numCellsY = ny + 2
        num_cells = numCellsX * numCellsY
        num_obj_2D = num_cells

        allocate( grid%space(1)%objects_per_dimension(3)%object(num_obj_2D) )

        do i = 1, num_obj_2D
            !! Set list of indices
            allocate( grid%space(1)%objects_per_dimension(3)%object(i)% &
                &   nodes(4) ) !! Each cell is formed by 4 nodes
        enddo
        cellId = 1
        do j = 1,numCellsY
            do i = 1, numCellsX
                !! Fill list of indices
                grid%space(1)%objects_per_dimension(3)%object(cellId)%  &
                    &   nodes(1) = cellId+0*numCellsX*numCellsY
                grid%space(1)%objects_per_dimension(3)%object(cellId)%  &
                    &   nodes(2) = cellId+1*numCellsX*numCellsY
                grid%space(1)%objects_per_dimension(3)%object(cellId)%  &
                    &   nodes(3) = cellId+3*numCellsX*numCellsY
                grid%space(1)%objects_per_dimension(3)%object(cellId)%  &
                    &   nodes(4) = cellId+2*numCellsX*numCellsY
                cellId = cellId + 1
            enddo
        enddo

        !! === SET UP GRID SUBSETS ===

        !! Two grid_subsets will be created, one for all nodes
        !! in the full domain, and one for
        !! all cells in the full domain

        !! Set up two grid subsets

        num_gridSubsets = 2         !! Number of grid subsets to write
                                    !! (Cells and Nodes)
        !! Allocate grid subsets
        allocate( grid%grid_subset(num_gridSubsets) )

        !! --- Set up grid_subset "Cells" for all cells in the domain ---
        !! (grid subset index: 1, objects forming the grid subset: cells, 2D)

        gridSubset_index = 1    !! Grid subset index of grid subset Cells
                                !! (Indexing of grid subsets follows the IDS
                                !! examples :
                                !! shot: 1, run:1, # device: iter, and
                                !! shot: 16151, run: 1000, device: aug
        gridSubset_name = "Cells"
        gridSubset_dim_index = 3 !! Grid subset Nodes consists of
                                 !! points -> 0D objects -> dimension index = 1
                                 !! (edges -> 1D objects -> dimension index = 2
                                 !! cells  -> 2D objects -> dimension index = 3)

        !! Set and fill grid subset name
        allocate( grid%grid_subset(gridSubset_index)%identifier%name(1) )
        grid%grid_subset(gridSubset_index)%identifier%name = gridSubset_name
        !! Fill subset base index ( this is special index, each subset has its
        !! unique index!)
        grid%grid_subset(gridSubset_index)%identifier%index = gridSubset_index

        !! Set subset element object reference to to appropriate space object
        !! node
        allocate( grid%grid_subset(gridSubset_index)%element(num_obj_2D) )
        do i = 1, num_obj_2D
            allocate( grid%grid_subset(gridSubset_index)%element(i)%object(1) )
            !! Fill subset space node index
            grid%grid_subset(gridSubset_index)%element(i)%object(1)%space = 1
            !! Fill subset object dimension i.e. object class
            grid%grid_subset(gridSubset_index)%element(i)%object(1)%dimension = &
                & gridSubset_dim_index
            !! Fill subgset object index
            grid%grid_subset(gridSubset_index)%element(i)%object(1)%index = i
        enddo

        !! --- Set up grid_subset "Nodes" for all nodes in the domain ---
        !! (grid subset index: 2, objects forming the grid subset: nodes/points)

        gridSubset_index = 2    !! Grid subset index of grid subset Nodes
                                !! (Indexing of grid subsets follows the IDS
                                !! examples :
                                !! shot: 1, run:1, # device: iter, and
                                !! shot: 16151, run: 1000, device: aug
        gridSubset_name = "Nodes"
        gridSubset_dim_index = 1    !! Grid subset Nodes consists of
                                    !! points -> 0D objects -> dimension index = 1
                                    !! (edges -> 1D objects -> dimension index = 2
                                    !! cells  -> 2D objects -> dimension index = 3)

        allocate( grid%grid_subset(gridSubset_index)%identifier%name(1) )
        grid%grid_subset(gridSubset_index)%identifier%name = gridSubset_name
        grid%grid_subset(gridSubset_index)%identifier%index = gridSubset_index

        allocate( grid%grid_subset(gridSubset_index)%element(num_obj_0D) )
        do i = 1, num_obj_0D
            allocate( grid%grid_subset(gridSubset_index)%element(i)%object(1) )
            grid%grid_subset(gridSubset_index)%element(i)%object(1)%space = 1
            grid%grid_subset(gridSubset_index)%element(i)%object(1)%    &
                &   dimension = gridSubset_dim_index
            grid%grid_subset(gridSubset_index)%element(i)%object(1)%index = i
        enddo

        !! === SET VALUES (ne, te, ti) for "Cells" grid subset ===
        !! --- Set ne (electron density) ---
        gridSubset_index = 1
        num_ne_gridSubset = 1
        num_ne_values = size(ne)
        allocate( ggd%electrons%density(num_ne_gridSubset) )
        ggd%electrons%density(num_ne_gridSubset)%grid_subset_index =    &
            &   gridSubset_index
        allocate( ggd%electrons%density(num_ne_gridSubset)%values(num_ne_values) )
        do n = 1, num_ne_values
            ggd%electrons%density(num_ne_gridSubset)%values(n) = ne(n)
        enddo

        !! --- Set te (electron temperature) ---
        num_te_gridSubset = 1
        num_te_values = size(te)
        allocate( ggd%electrons%temperature(num_te_gridSubset) )
        ggd%electrons%temperature(num_te_gridSubset)%grid_subset_index =    &
            &   gridSubset_index
        allocate( ggd%electrons%temperature(num_te_gridSubset)%values(num_te_values) )
        do n = 1, num_te_values
            !! convert to eV (1 J = 6.242e18 eV)
            ggd%electrons%temperature(num_te_gridSubset)%values(n) =    &
                &   te(n)*6.242e18
        enddo

        !! --- Set ti (ion temperature) ---
        num_ti_gridSubset   = 1
        num_ti_values       = size(ti)
        num_ti_species      = 1     !! Number of ion species, as in
                                    !! number of different ion charges.
        ion_specie = 1
        !! Ion specie is linked with the ion density of each ion charge,
        !! as ion temperature is taken as the same for all ion charges.
        allocate( ggd%ion(num_ti_species) )
        allocate( ggd%ion(ion_specie)%temperature(num_ti_gridSubset) )
        ggd%ion(ion_specie)%temperature(num_ti_gridSubset)%grid_subset_index =  &
            &   gridSubset_index
        allocate( ggd%ion(ion_specie)%temperature(num_ti_gridSubset)%   &
            &   values(num_ti_values) )
        do n = 1, num_ti_values
            !! convert to eV (1 J = 6.242e18 eV)
            ggd%ion(ion_specie)%temperature(num_ti_gridSubset)%values(n) =  &
                &   ti(n) * (6.242e18)
        enddo

        !! Set data to edge_profiles IDS
        write(0,*) "Writing to edge_profiles IDS"

        !! Create and modify new shot/run
        call imas_create_env(treename, shot, run, 0, 0, idx, username,  &
            &   device, version)

        !! Or open and modify existing shot/run (might work much faster than
        !! imas_create_env)
        ! call imas_open_env('treename', shot, run, idx, username, device, version)

        !! Put data to IDS
        call ids_put_slice( idx, "edge_profiles", edge_profiles )
        call ids_put( idx, "edge_profiles", edge_profiles )

        !! Close IDS
        call ids_deallocate( edge_profiles )
        call imas_close(idx)

        write(0,*) "IDS write finished"

    end subroutine write_ids_edge_profiles

    !> Example subroutine for reading data out from the edge_profiles IDS
    !! with Fortran90
    subroutine read_ids( treename, shot, run, idx, username, device, version )
        character(len=24) :: treename   !< The name of the IMAS IDS database
        character(len=24) :: username   !< Creator/owner of the IMAS IDS database
        character(len=24) :: device !< Device name of the IMAS IDS database
            !< (i. e. solps-iter, iter, aug)
        character(len=24) :: version    !< Major version of the IMAS IDS database
        integer :: shot !< The shot number of the database being created
        integer :: run  !< The run number of the database being created
        integer :: idx  !< The returned identifier to be used in the subsequent
        !! Internal variables
        integer :: gridSubset_index !< >Grid subset base index
        type(ids_edge_profiles) :: edge_profiles    !< IDS designed to store
            !< data in edge plasma profiles  (includes the scrape-off layer and
            !<  possibly part of the confined plasma)

        gridSubset_index = 2

        !! Open input datafile from local database
        write (0,*) "Started reading input IDS", idx, shot, run

        call imas_open_env('treename', shot, run, idx, username, device, version)
        call ids_get(idx, 'edge_profiles', edge_profiles)

        write(0,*) 'homogeneous_time = ',   &
            &   edge_profiles%ids_properties%homogeneous_time
        write(0,*) "Grid subset 2 name = ", &
            &   edge_profiles%ggd(1)%grid%grid_subset(gridSubset_index)%    &
            &   identifier%name
        write(0,*) "Grid subset 2 index = ",    &
            &   edge_profiles%ggd(1)%grid%grid_subset(gridSubset_index)%    &
            &   identifier%index
        ! write(0,*) "Time = ", edge_profiles%time(1)
        call ids_deallocate( edge_profiles )
        call imas_close(idx)
        write (0,*) "Finished reading input IDS"

    end subroutine read_ids

end program b2_ual_write_deprecated

!!!Local Variables:
!!! mode: f90
!!! End:
