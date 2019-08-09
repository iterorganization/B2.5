!!-----------------------------------------------------------------------------
!! DOCUMENTATION:
!>      @section b2uw_ualio_desc Description
!!      Module providing main routines for setting processed B2.5 data for
!!      ITM CPO edge or IMAS edge_profiles, edge_sources and edge_transport IDSs.
!!
!!      @detail The data comprises of grid geometry, grid subsets and data
!!              fields of various quantities.
!!
!!      @subsection b2uw_ualio_syx     Exceptional syntax explanation
!!      @code
!!          ! IGNORE    !! syntax used to ignore this module in list
!!                      !! dependency when compiling the code
!!      @endcode
!!
!!-----------------------------------------------------------------------------
module b2mod_ual_io

    !! B2 modules

    use b2mod_types
    use b2mod_b2cmpa
    use b2mod_geo
    use b2mod_diag
    use b2mod_rates
    use b2mod_plasma
    use b2mod_elements
    use b2mod_constants
    use b2mod_sources
    use b2mod_transport
    use b2mod_anomalous_transport
    use b2mod_boundary_namelist
    use b2mod_neutrals_namelist
    use b2mod_user_namelist
    use b2mod_indirect
    use b2mod_interp
    use b2mod_b2cmrc
    use b2mod_b2cmfs
    use b2mod_version
    use b2mod_grid_mapping
    use b2mod_b2plot &
    , only : textan, textmn, textin, nxtl, nxtr, jxa, jsep
#ifdef B25_EIRENE
    use eirmod_comusr &
    , only : natmi, nmoli, nioni, nmassa, nchara, nchrgi, nchari
#else
    use b2mod_b2plot &
    , only : natmi
#endif

    use logging

    !! UAL Access
#ifdef IMAS
    !! B2/CPO Mapping
    use b2mod_ual_io_data &
     & , only : b2_IMAS_Transform_Data_B2_To_IDS, &
     &          b2_IMAS_Transform_Data_B2_To_IDS_Vertex
    use b2mod_ual_io_grid &
     & , only : findGridSubsetByName, GridWriteData, &
     &          b2_IMAS_Fill_Grid_Desc
    use ids_routines       ! IGNORE
    use ids_schemas        ! IGNORE
#if IMAS_MINOR_VERSION > 14
    use ids_utility        ! IGNORE
#endif
    use ids_grid_common , &     ! IGNORE
        &   IDS_COORDTYPE_R => COORDTYPE_R,       &
        &   IDS_COORDTYPE_Z => COORDTYPE_Z,       &
        &   IDS_GRID_UNDEFINED => GRID_UNDEFINED
#else
#ifdef ITM_ENVIRONMENT_LOADED
    use euITM_schemas   ! IGNORE
    use euITM_routines  ! IGNORE
    use itm_grid_common ! IGNORE
#endif
#endif

  implicit none

#ifdef IMAS

  logical, parameter, private :: INCLUDE_GHOST_CELLS = .false.

contains

    !> Process B2.5 data and set it to IMAS IDS.
    !! @note    The \b B25_process_ids routine enables to store data for
    !!          specific time slice. By default it stores single default
    !!          time slice of time slice value 0.0.
    !!          \b num_time_slices_IN is required to beforehand allocate
    !!          required ggd(:) array of nodes structure and for additional
    !!          checks for correct use of the routine.
    !! @note    Time slice value is set as:
    !!          \b time_slice_value = \b time_step_IN * \b time_slice_ind_IN
    subroutine B25_process_ids( edge_profiles, edge_sources, edge_transport, &
            &   radiation, &
#if IMAS_MINOR_VERSION > 21
            &   summary, &
#endif
            &   time_IN, time_step_IN, time_slice_ind_IN, num_time_slices_IN )
#       include <git_version_B25.h>
        type (ids_edge_profiles) :: edge_profiles    !< IDS designed to
            !< store data on edge plasma profiles  (includes the scrape-off
            !< layer and possibly part of the confined plasma)
        type (ids_edge_sources) :: edge_sources !< IDS designed to store
            !< data on edge plasma sources. Energy terms correspond to the full
            !< kinetic energy equation (i.e. the energy flux takes into account
            !< the energy transported by the particle flux)
        type (ids_edge_transport) :: edge_transport !< IDS designed to store
            !< data on edge plasma transport. Energy terms correspond to the
            !< full kinetic energy equation (i.e. the energy flux takes into
            !< account the energy transported by the particle flux)
        type (ids_radiation) :: radiation !< IDS designed to store
            !< data on radiation emitted by the plasma species
#if IMAS_MINOR_VERSION > 21
        type (ids_summary) :: summary !< IDS designed to store
            !< run summary data
#endif
        real(IDS_real), intent(in), optional :: time_IN !< Time
        real(IDS_real), intent(in), optional :: time_step_IN !< Time step
        integer, intent(in), optional :: time_slice_ind_IN
            !< Time step index for the current time slice
        integer, intent(in), optional :: num_time_slices_IN
            !< Total number of time steps. It is required to beforehand allocate
            !< required ggd(:) array of nodes structure and for additional
            !< checks for correct use of the routine.

        !! Internal variables
        character(len=24) :: ion_label  !< Ion species label (e.g. D+1)
        character(len=24) :: mol_label  !< Molecule species label (e.g. D2)
        character(len=12) :: ion_charge !< Ion charge (e.g. '1', '2', etc.)
        character(len=24) :: source     !< Code source
        character(len=2)  :: plate_name(4) !< Divertor plate name
        integer :: ion_charge_int !< Ion charge (e.g. 1, 2, etc.)
        integer :: ion_label_tlen !< Length of the (trimmed) ion label
        integer :: ns   !< Total number of ion species
        integer :: nx   !< Specifies the number of interior cells
                        !< along the first coordinate
        integer :: ny   !< Specifies the number of interior cells
                        !< along the second coordinate
        integer :: n_process !< Number of radiation processes handled
        integer :: nelems    !< Number of elements present in a molecule or molecular ion
        integer :: is     !< Species index (iterator)
        integer :: i      !< Iterator
        integer :: j      !< Iterator
        integer :: k      !< Iterator
        integer :: ix     !< Iterator
        integer :: iy     !< Iterator
        integer :: iatm   !< Atom iterator
        integer :: iatm1  !< Hydrogenic atom index in molecule composition
        integer :: iatm2  !< Non-hydrogenic atom index in molecule composition
        integer :: istrai !< Stratum iterator
        integer :: is1    !< First ion of an isonuclear sequence
        integer :: is2    !< Last ion of an isonuclear sequence
        integer :: icnt   !< Boundary cell counter
        integer :: ntrgts !< Number of divertor targets
        integer :: o      !< Dummy integer
        integer :: p      !< Dummy integer
        integer :: ixpos(4), iypos(4) !< Target positions
        integer :: GeometryType !< Geometry identifier number
        integer :: iGsCoreBoundary  !< Variable to hold Core grid subset base
            !< index, later found by findGridSubsetByName() routine.
        integer :: iGsInnerMidplane !< Variable to hold Inner Midplane grid
            !< subset base index, later found by findGridSubsetByName() routine
        integer :: iGsOuterMidplane !< Variable to hold Outer Midplane grid
            !< subset base index, later found by findGridSubsetByName() routine
        integer :: iGsCore  !< Variable to hold Core grid
            !< subset base index, later found by findGridSubsetByName() routine
        integer :: iGsSOL   !< Variable to hold SOL grid
            !< subset base index, later found by findGridSubsetByName() routine
        integer :: iGsIDivertor     !< Variable to hold Inner Divertor grid
            !< subset base index, later found by findGridSubsetByName() routine
        integer :: iGsODivertor     !< Variable to hold Outer Divertor grid
            !< subset base index, later found by findGridSubsetByName() routine
        integer :: homogeneous_time !< Homogeneous time (0 or 1)
        integer :: num_time_slices  !< Total number of time slices.
        integer :: time_sind   !< Time slice index. Also General grid
            !< description slice identifier
        logical, parameter :: B2_WRITE_DATA = .true.
        real(IDS_real),   &
            &   dimension( -1:ubound( crx, 1 ), -1:ubound( crx, 2), 3, 3) :: e
        real(IDS_real) :: tmpFace( -1:ubound( na, 1), -1:ubound( na, 2 ), 0:1)
        real(IDS_real) :: tmpVx( -1:ubound( na, 1), -1:ubound( na, 2 ) )
        real(IDS_real) :: tmpCv( -1:ubound( na, 1), -1:ubound( na, 2 ) )
        real(IDS_real) :: zeff( -1:ubound( na, 1), -1:ubound( na, 2 ) )
        real(IDS_real) :: sna0( -1:ubound( na, 1), -1:ubound( na, 2 ), 0:1, &
                       &        -1:ubound( na, 3) )
        real(IDS_real) :: smo0( -1:ubound( na, 1), -1:ubound( na, 2 ), 0:3, &
                       &        -1:ubound( na, 3) )
        real(IDS_real) :: she0( -1:ubound( na, 1), -1:ubound( na, 2 ), 0:3)
        real(IDS_real) :: shi0( -1:ubound( na, 1), -1:ubound( na, 2 ), 0:3)
        real(IDS_real) :: time  !< Generic time
        real(IDS_real) :: time_step !< Time step
        real(IDS_real) :: time_slice_value   !< Time slice value
        real(IDS_real) :: b0, r0, b0r0, vtor, nisep, nasum, gsum, gmid, gbot, gtop
        type(B2GridMap) :: gmap !< Data structure holding an
            !< intermediate grid description to be transferred into a CPO or IDS
        type(ids_generic_grid_dynamic_grid_subset) :: gs_cell
        type(ids_generic_grid_dynamic_grid_subset) :: gs_face
        type(ids_generic_grid_dynamic_grid_subset) :: gs_bnd_core

        integer, save :: ismain = 1
        integer, save :: use_eirene = 0
        integer, save :: target_offset = 1
        integer, save :: nesepm_istra = -1
        real(IDS_real), save :: ndes = 0.0_IDS_real
        real(IDS_real), save :: ndes_sol = 0.0_IDS_real
        real(IDS_real), save :: nesepm_pfr = 0.0_IDS_real
        real(IDS_real), save :: nesepm_sol = 0.0_IDS_real
        real(IDS_real), save :: nepedm_sol = 0.0_IDS_real
        real(IDS_real), save :: volrec_sol = 0.0_IDS_real
        real(IDS_real), save :: private_flux_puff = 0.0_IDS_real
        real(IDS_real), save :: BoRiS
        character*8 date
        character*10 ctime
        character*5 zone
        integer tvalues(8)
        character*16 usrnam
        character*8 imas_version, ual_version
        logical match_found, streql
#ifdef USE_PXFGETENV
        integer lenval, ierror
#else
#ifdef NAGFOR
        integer lenval, ierror
#endif
#endif
        external usrnam, streql

        !! ===  SET UP IDS ===
        write(0,*) "Setting data for edge_profiles IDS"
        call ipgetr ('b2news_BoRiS', BoRiS)
        call ipgeti ('b2mndr_ismain', ismain)
        call ipgeti ('b2mndr_eirene', use_eirene)
        call ipgeti ('b2mwti_target_offset', target_offset)
        call ipgetr ('b2stbc_ndes', ndes)
        call ipgetr ('b2stbc_ndes_sol', ndes_sol)
        call ipgetr ('b2stbc_nesepm_pfr', nesepm_pfr)
        call ipgetr ('b2stbc_nesepm_sol', nesepm_sol)
        call ipgetr ('b2stbc_nepedm_sol', nepedm_sol)
        call ipgetr ('b2stbc_volrec_sol', volrec_sol)
        call ipgetr ('b2stbc_private_flux_puff', private_flux_puff)
        call ipgeti ('eirene_nesepm_istra', nesepm_istra)
        if(nesepm_istra.gt.0) then
         call xertst (nesepm_istra.le.nstrat,'faulty internal parameter nesepm_istra')
         call xertst (crcstra(nesepm_istra).eq.'C', &
           &  'Stratum nesepm_istra is not declared as a puff stratum!')
        else
         do istrai = 1, nstrat
           if(nesepm_istra.le.0.and.CRCSTRA(istrai).EQ.'C') then
             nesepm_istra=istrai
           endif
         enddo
        endif
        call date_and_time (date, ctime, zone, tvalues)
#ifdef NAGFOR
        call get_environment_variable('IMAS_VERSION',status=ierror, length=lenval)
        if (ierror.eq.0) call get_environment_variable('IMAS_VERSION',value=imas_version)
        call get_environment_variable('UAL_VERSION',status=ierror, length=lenval)
        if (ierror.eq.0) call get_environment_variable('UAL_VERSION',value=ual_version)
#else
#ifdef USE_PXFGETENV
        CALL PXFGETENV ('IMAS_VERSION', 0, imas_version, lenval, ierror)
        CALL PXFGETENV ('UAL_VERSION', 0, ual_version, lenval, ierror)
#else
        call getenv ('IMAS_VERSION', imas_version)
        call getenv ('UAL_VERSION', ual_version)
#endif
#endif

#ifdef B25_EIRENE
        source = "SOLPS-ITER"
#else
        source = "B2.5"
#endif
        ns = size( na, 3 )
        nx = ubound( na, 1 )
        ny = ubound( na, 2 )
        call b2xzef (nx, ny, ns, rz2, na, ne, zeff)

        !! Preparing database for writing
        !! Through practice it was disclosed that there are some mandatory
        !! steps to be done in order to assure for data to be successfully
        !! written to IDS. Without going through those steps errors and failed
        !! process of writing to IDS are to be expected.
        !! This can be done using setIDSFundamentals routine
        homogeneous_time = 1
        if ( present( time_IN ) ) then
            time = time_IN
        else
            time = 0.0_IDS_real
        end if

        !! Set default time step values
        time_sind = 1
        time_slice_value = 0.0_IDS_real
        time_step = IDS_REAL_INVALID
        num_time_slices = 1
        !! If present, set time step values
        if( present( time_step_IN ) ) time_step = time_step_IN
        if( present( time_slice_ind_IN ) ) then
            time_sind = time_slice_ind_IN
        !! Get time slice value
            time_slice_value = time_sind * time_step
        else
            time_slice_value = time
        end if
        if( present( num_time_slices_IN ) ) num_time_slices = num_time_slices_IN
        !! Check if num_time_slices >= time_sind
        call xertst( num_time_slices .ge. time_sind, &
            & "B25_process_ids: Time step index cannot be greater " // &
            & "than total number of time steps!" )

        !! Preparing edge_profiles IDS for writing
        !! In order to write to IDS database there are next steps that are
        !! mandatory to do, otherwise there is high change that writing to IDS
        !! database will fail
        !! 1. Set homogeneous_time to 0 or 1
        edge_profiles%ids_properties%homogeneous_time = homogeneous_time
        allocate( edge_profiles%ids_properties%comment(1) )
        edge_profiles%ids_properties%comment(1) = "Done by b2_ual_write"
        !! 2. Allocate edge_profiles.time and set it to desired values
        allocate( edge_profiles%time(num_time_slices) )
        do i = 1, num_time_slices
          edge_profiles%time(i) = time - (num_time_slices-i) * time_step
        end do

        !! Preparing edge_transport IDS for writing
        !! In order to write to IDS database there are next steps that are
        !! mandatory to do, otherwise there is high change that writing to IDS
        !! database will fail
        !! 1. Set homogeneous_time to 0 or 1
        edge_transport%ids_properties%homogeneous_time = homogeneous_time
        allocate( edge_transport%ids_properties%comment(1) )
        edge_transport%ids_properties%comment(1) = "Done by b2_ual_write"
        !! 2. Allocate edge_transport.time and set it to desired values
        allocate( edge_transport%time(num_time_slices) )
        do i = 1, num_time_slices
          edge_transport%time(i) = time - (num_time_slices-i) * time_step
        end do

        !! Preparing edge_sources IDS for writing
        !! In order to write to IDS database there are next steps that are
        !! mandatory to do, otherwise there is high change that writing to IDS
        !! database will fail
        !! 1. Set homogeneous_time to 0 or 1
        edge_sources%ids_properties%homogeneous_time = homogeneous_time
        allocate( edge_sources%ids_properties%comment(1) )
        edge_sources%ids_properties%comment(1) = "Done by b2_ual_write"
        !! 2. Allocate edge_sources.time and set it to desired values
        allocate( edge_sources%time(num_time_slices) )
        do i = 1, num_time_slices
          edge_sources%time(i) = time - (num_time_slices-i) * time_step
        end do

#if IMAS_MINOR_VERSION > 21
        !! Preparing summary IDS for writing
        !! In order to write to IDS database there are next steps that are
        !! mandatory to do, otherwise there is high change that writing to IDS
        !! database will fail
        !! 1. Set homogeneous_time to 0 or 1
        summary%ids_properties%homogeneous_time = homogeneous_time
        allocate( summary%ids_properties%comment(1) )
        summary%ids_properties%comment(1) = "Done by b2_ual_write"
        !! 2. Allocate summary.time and set it to desired values
        allocate( summary%time(num_time_slices) )
        do i = 1, num_time_slices
          summary%time(i) = time - (num_time_slices-i) * time_step
        end do
#endif

        !! Preparing radiation IDS for writing
        !! In order to write to IDS database there are next steps that are
        !! mandatory to do, otherwise there is high change that writing to IDS
        !! database will fail
        !! 1. Set homogeneous_time to 0 or 1
        radiation%ids_properties%homogeneous_time = homogeneous_time
        allocate( radiation%ids_properties%comment(1) )
        radiation%ids_properties%comment(1) = "Done by b2_ual_write_b2mod"
        !! 2. Allocate radiation.time and set it to desired values
        allocate( radiation%time(num_time_slices) )
        do i = 1, num_time_slices
          radiation%time(i) = time - (num_time_slices-i) * time_step
        end do
        !! 3. Allocate radiation.process
        !! Process 1: line and recombination radiation due to B2.5 species
        !! Process 2: bremsstrahlung recombination due to B2.5 species
        !! Process 3: line radiation due to Eirene neutrals (atoms and molecules)
        !! Process 4: line radiation dur to Eirene molecular ions
        if (use_eirene.ne.0) then
          n_process = 4
        else
          n_process = 2
        end if
        allocate( radiation%process(n_process) )
        do j = 1, n_process
          allocate( radiation%process(j)%identifier%name(1) )
          allocate( radiation%process(j)%identifier%description(1) )
        end do
        radiation%process(1)%identifier%name = 'line_radiation'
        radiation%process(1)%identifier%description = 'Line and rec. rad. from B2.5 species'
        radiation%process(2)%identifier%name = 'bremsstrahlung'
        radiation%process(2)%identifier%description = 'Bremsstrahlung from B2.5 species'
        radiation%process(1)%identifier%index = 2
        radiation%process(2)%identifier%index = 8
        if (use_eirene.ne.0) then
          radiation%process(3)%identifier%index = 1
          radiation%process(4)%identifier%index = 2
          radiation%process(3)%identifier%name = 'line_radiation'
          radiation%process(3)%identifier%description = 'Line radiation from Eirene neutrals'
          radiation%process(4)%identifier%name = 'line_radiation'
          radiation%process(4)%identifier%description = 'Line radiation from Eirene mol. ions'
        end if

        write(*,*) "Running b2CreateMap subroutine"
        !! Set up the B2<->IDS mappings
        call b2CreateMap( nx, ny, crx( -1:nx, -1:ny, : ),             &
            &   cry( -1:nx, -1:ny, : ), cflags, leftix, leftiy,       &
            &   rightix, rightiy, topix, topiy, bottomix,bottomiy,    &
            &   INCLUDE_GHOST_CELLS, gmap )
        mapInitialized = .true.

        !! To be done only on the run of the first time step. If already
        !! allocated structure is allocated again it deletes all previously
        !! set data!
        !! Check for edge_transport%model(1)%ggd and edge_sources%source(1)%ggd
        !! is not included as they contain beforehand model(:) / source(:)
        !! structures
        if( size( edge_profiles%ggd ) .ne. num_time_slices ) then
            !! Allocate ggd for number of different time steps
            time_sind = 1
            allocate( edge_profiles%ggd( num_time_slices ) )
#if IMAS_MINOR_VERSION > 14
            allocate( edge_profiles%grid_ggd( time_sind ) )
            allocate( edge_transport%grid_ggd( time_sind ) )
            allocate( edge_sources%grid_ggd( time_sind ) )
#if IMAS_MINOR_VERSION > 21
            allocate( radiation%grid_ggd( time_sind ) )
#endif
#endif
            allocate( edge_transport%model(1) )
            edge_transport%model(1)%identifier%index = 1
            allocate( edge_sources%source(1) )
            edge_sources%source(1)%identifier%index = 1

            allocate( edge_transport%model(1)%ggd( num_time_slices ) )
            allocate( edge_sources%source(1)%ggd( num_time_slices ) )
            allocate( edge_transport%model(1)%ggd( num_time_slices )% &
                &   electrons%energy%flux(1) )
        end if

        !! Allocate and init the IDS
        allocate( edge_profiles%code%name(1) )
        allocate( edge_transport%code%name(1) )
        allocate( edge_sources%code%name(1) )
        allocate( radiation%code%name(1) )
        edge_profiles%code%name = source
        edge_transport%code%name = source
        edge_sources%code%name = source
        radiation%code%name = source

        allocate( edge_profiles%code%version(1) )
        edge_profiles%code%version = newversion
        allocate( edge_transport%code%version(1) )
        edge_transport%code%version = newversion
        allocate( edge_sources%code%version(1) )
        edge_sources%code%version = newversion
        allocate( radiation%code%version(1) )
        radiation%code%version = newversion

        allocate( edge_profiles%code%commit(1) )
        edge_profiles%code%commit = git_version_B25
        allocate( edge_transport%code%commit(1) )
        edge_transport%code%commit = git_version_B25
        allocate( edge_sources%code%commit(1) )
        edge_sources%code%commit = git_version_B25
        allocate( radiation%code%commit(1) )
        radiation%code%commit = git_version_B25

        allocate( edge_profiles%code%repository(1) )
        edge_profiles%code%repository = "git.iter.org"
        allocate( edge_transport%code%repository(1) )
        edge_transport%code%repository = "git.iter.org"
        allocate( edge_sources%code%repository(1) )
        edge_sources%code%repository = "git.iter.org"
        allocate( radiation%code%repository(1) )
        radiation%code%repository = "git.iter.org"

        allocate( radiation%ids_properties%source(1) )
        radiation%ids_properties%source = b2frates_flag

        allocate( edge_profiles%ids_properties%provider(1) )
        edge_profiles%ids_properties%provider = usrnam()
        allocate( edge_transport%ids_properties%provider(1) )
        edge_transport%ids_properties%provider = usrnam()
        allocate( edge_sources%ids_properties%provider(1) )
        edge_sources%ids_properties%provider = usrnam()
        allocate( radiation%ids_properties%provider(1) )
        radiation%ids_properties%provider = usrnam()
#if IMAS_MINOR_VERSION > 21
        allocate( summary%code%name(1) )
        summary%code%name = source
        allocate( summary%code%version(1) )
        summary%code%version = newversion
        allocate( summary%code%commit(1) )
        summary%code%commit = git_version_B25
        allocate( summary%code%repository(1) )
        summary%code%repository = "git.iter.org"
        allocate( summary%ids_properties%provider(1) )
        summary%ids_properties%provider = usrnam()
#endif

        allocate( edge_profiles%ids_properties%creation_date(1) )
        edge_profiles%ids_properties%creation_date = &
                &   date//' '//ctime//' '//' '//zone
        allocate( edge_transport%ids_properties%creation_date(1) )
        edge_transport%ids_properties%creation_date = &
                &   date//' '//ctime//' '//' '//zone
        allocate( edge_sources%ids_properties%creation_date(1) )
        edge_sources%ids_properties%creation_date = &
                &   date//' '//ctime//' '//' '//zone
        allocate( radiation%ids_properties%creation_date(1) )
        radiation%ids_properties%creation_date = &
                &   date//' '//ctime//' '//' '//zone
#if IMAS_MINOR_VERSION > 21
        allocate( summary%ids_properties%creation_date(1) )
        summary%ids_properties%creation_date = &
                &   date//' '//ctime//' '//' '//zone
#endif

#if IMAS_MINOR_VERSION > 21
        allocate( edge_profiles%ids_properties%version_put%data_dictionary(1) )
        edge_profiles%ids_properties%version_put%data_dictionary = imas_version
        allocate( edge_transport%ids_properties%version_put%data_dictionary(1) )
        edge_transport%ids_properties%version_put%data_dictionary = imas_version
        allocate( edge_sources%ids_properties%version_put%data_dictionary(1) )
        edge_sources%ids_properties%version_put%data_dictionary = imas_version
        allocate( radiation%ids_properties%version_put%data_dictionary(1) )
        radiation%ids_properties%version_put%data_dictionary = imas_version
        allocate( summary%ids_properties%version_put%data_dictionary(1) )
        summary%ids_properties%version_put%data_dictionary = imas_version

        allocate( edge_profiles%ids_properties%version_put%access_layer(1) )
        edge_profiles%ids_properties%version_put%access_layer = ual_version
        allocate( edge_transport%ids_properties%version_put%access_layer(1) )
        edge_transport%ids_properties%version_put%access_layer = ual_version
        allocate( edge_sources%ids_properties%version_put%access_layer(1) )
        edge_sources%ids_properties%version_put%access_layer = ual_version
        allocate( radiation%ids_properties%version_put%access_layer(1) )
        radiation%ids_properties%version_put%access_layer = ual_version
        allocate( summary%ids_properties%version_put%access_layer(1) )
        summary%ids_properties%version_put%access_layer = ual_version

        allocate( edge_profiles%ids_properties%version_put%access_layer_language(1) )
        edge_profiles%ids_properties%version_put%access_layer_language = 'FORTRAN'
        allocate( edge_transport%ids_properties%version_put%access_layer_language(1) )
        edge_transport%ids_properties%version_put%access_layer_language = 'FORTRAN'
        allocate( edge_sources%ids_properties%version_put%access_layer_language(1) )
        edge_sources%ids_properties%version_put%access_layer_language = 'FORTRAN'
        allocate( radiation%ids_properties%version_put%access_layer_language(1) )
        radiation%ids_properties%version_put%access_layer_language = 'FORTRAN'
        allocate( summary%ids_properties%version_put%access_layer_language(1) )
        summary%ids_properties%version_put%access_layer_language = 'FORTRAN'

        i=index(git_version_B25,'-')
        allocate( summary%tag%name(1) )
        summary%tag%name = git_version_B25(1:i-1)
        r0 = 0.0_R8
        icnt = 0
        do ix = -1, nx
          if (on_closed_surface(ix,-1)) then
            icnt = icnt + 1
            if (isymm.eq.1.or.isymm.eq.2) then
              r0 = r0 + crx(ix,-1,0)
            else if (isymm.eq.3.or.isymm.eq.4) then
              r0 = r0 + cry(ix,-1,0)
            end if
          end if
        end do
        if (icnt.gt.0) then
          r0 = r0 / float(icnt)
          if (ffbz(jxa,-1,0).ne.0.0_R8) then
            b0 = ffbz(jxa,-1,0)/r0
          else if (isymm.eq.1 .or. isymm.eq.2) then
            b0r0 = bb(jxa,-1,2)*(crx(jxa,-1,0)+crx(jxa,-1,1)+ &
                              &  crx(jxa,-1,2)+crx(jxa,-1,3))/4.0
            b0 = b0r0 / r0
          else if (isymm.eq.3 .or. isymm.eq.4) then
            b0r0 = bb(jxa,-1,2)*(cry(jxa,-1,0)+cry(jxa,-1,1)+ &
                              &  cry(jxa,-1,2)+cry(jxa,-1,3))/4.0
            b0 = b0r0 / r0
          end if
        end if
        summary%global_quantities%r0%value = r0
        allocate( summary%global_quantities%r0%source(1) )
        summary%global_quantities%r0%source = source
        allocate( summary%global_quantities%b0%value( time_sind ) )
        summary%global_quantities%b0%value( time_sind ) = b0
        allocate( summary%global_quantities%b0%source(1) )
        summary%global_quantities%b0%source = source
#endif

        !! Write grid & grid subsets/subgrids
#if IMAS_MINOR_VERSION < 15
        call b2_IMAS_Fill_Grid_Desc( gmap,                                  &
            &   edge_profiles%ggd( time_sind )%grid,                        &
            &   nx, ny, crx(-1:nx, -1:ny, :), cry(-1:nx, -1:ny, : ),        &
            &   leftix, leftiy, rightix, rightiy, topix, topiy, bottomix,   &
            &   bottomiy, nnreg, topcut, region, cflags,                    &
            &   INCLUDE_GHOST_CELLS, vol, gs, qc )
        call b2_IMAS_Fill_Grid_Desc( gmap,                                  &
            &   edge_transport%model(1)%ggd( time_sind )%grid,              &
            &   nx, ny, crx(-1:nx, -1:ny, :), cry(-1:nx, -1:ny, : ),        &
            &   leftix, leftiy, rightix, rightiy, topix, topiy, bottomix,   &
            &   bottomiy, nnreg, topcut, region, cflags,                    &
            &   INCLUDE_GHOST_CELLS, vol, gs, qc )
        call b2_IMAS_Fill_Grid_Desc( gmap,                                  &
            &   edge_sources%source(1)%ggd( time_sind )%grid,               &
            &   nx, ny, crx(-1:nx, -1:ny, :), cry(-1:nx, -1:ny, : ),        &
            &   leftix, leftiy, rightix, rightiy, topix, topiy, bottomix,   &
            &   bottomiy, nnreg, topcut, region, cflags,                    &
            &   INCLUDE_GHOST_CELLS, vol, gs, qc )
#else
        call b2_IMAS_Fill_Grid_Desc( gmap,                                  &
            &   edge_profiles%grid_ggd( time_sind ),                        &
            &   nx, ny, crx(-1:nx, -1:ny, :), cry(-1:nx, -1:ny, : ),        &
            &   leftix, leftiy, rightix, rightiy, topix, topiy, bottomix,   &
            &   bottomiy, nnreg, topcut, region, cflags,                    &
            &   INCLUDE_GHOST_CELLS, vol, gs, qc )
        call b2_IMAS_Fill_Grid_Desc( gmap,                                  &
            &   edge_transport%grid_ggd( time_sind ),                       &
            &   nx, ny, crx(-1:nx, -1:ny, :), cry(-1:nx, -1:ny, : ),        &
            &   leftix, leftiy, rightix, rightiy, topix, topiy, bottomix,   &
            &   bottomiy, nnreg, topcut, region, cflags,                    &
            &   INCLUDE_GHOST_CELLS, vol, gs, qc )
        call b2_IMAS_Fill_Grid_Desc( gmap,                                  &
            &   edge_sources%grid_ggd( time_sind ),                         &
            &   nx, ny, crx(-1:nx, -1:ny, :), cry(-1:nx, -1:ny, : ),        &
            &   leftix, leftiy, rightix, rightiy, topix, topiy, bottomix,   &
            &   bottomiy, nnreg, topcut, region, cflags,                    &
            &   INCLUDE_GHOST_CELLS, vol, gs, qc )
#if IMAS_MINOR_VERSION > 21
        call b2_IMAS_Fill_Grid_Desc( gmap,                                  &
            &   radiation%grid_ggd( time_sind ),                            &
            &   nx, ny, crx(-1:nx, -1:ny, :), cry(-1:nx, -1:ny, : ),        &
            &   leftix, leftiy, rightix, rightiy, topix, topiy, bottomix,   &
            &   bottomiy, nnreg, topcut, region, cflags,                    &
            &   INCLUDE_GHOST_CELLS, vol, gs, qc )
#endif
#endif

        !! Allocate and set time slice value
#if IMAS_MINOR_VERSION > 14
        edge_profiles%grid_ggd( time_sind )%time = time_slice_value
        edge_transport%model(1)%ggd( time_sind )%time = time_slice_value
        edge_sources%grid_ggd( time_sind )%time = time_slice_value
#if IMAS_MINOR_VERSION > 21
        radiation%grid_ggd( time_sind )%time = time_slice_value
#endif
#endif
        edge_profiles%ggd( time_sind )%time = time_slice_value
        edge_transport%model(1)%ggd( time_sind )%time = time_slice_value
        edge_sources%source(1)%ggd( time_sind )%time = time_slice_value
#if IMAS_MINOR_VERSION > 21
        do j = 1, n_process
          allocate( radiation%process(j)%ggd( num_time_slices ) )
          radiation%process(j)%ggd( time_sind )%time = time_slice_value
        end do
#endif

        !! List of species
        allocate( edge_profiles%ggd( time_sind )%ion( ns ) )
#if IMAS_MINOR_VERSION > 21
        allocate( radiation%process(1)%ggd( time_sind )%ion( ns ) )
        allocate( radiation%process(2)%ggd( time_sind )%ion( ns ) )
#endif
        do is = 0, ns-1
            allocate( edge_profiles%ggd( time_sind )%ion( is + 1 )%label(1) )
            allocate( edge_profiles%ggd( time_sind )%ion( is + 1 )%state(1) )
            allocate( edge_profiles%ggd( time_sind )%ion( is + 1 )%state(1)%label(1) )
            allocate( edge_profiles%ggd( time_sind )%ion( is + 1 )%element(1) )

#if IMAS_MINOR_VERSION > 21
            allocate( radiation%process(1)%ggd( time_sind )%ion( is + 1 )%label(1) )
            allocate( radiation%process(1)%ggd( time_sind )%ion( is + 1 )%state(1) )
            allocate( radiation%process(1)%ggd( time_sind )%ion( is + 1 )%state(1)%label(1) )
            allocate( radiation%process(1)%ggd( time_sind )%ion( is + 1 )%element(1) )
            allocate( radiation%process(2)%ggd( time_sind )%ion( is + 1 )%label(1) )
            allocate( radiation%process(2)%ggd( time_sind )%ion( is + 1 )%state(1) )
            allocate( radiation%process(2)%ggd( time_sind )%ion( is + 1 )%state(1)%label(1) )
            allocate( radiation%process(2)%ggd( time_sind )%ion( is + 1 )%element(1) )
#endif
            ! Put label to ion(is + 1).state(1).label
            call species( is, edge_profiles%ggd( time_sind )%ion( is + 1 )% &
                &   state(1)%label, .false.)
#if IMAS_MINOR_VERSION > 21
            call species( is, radiation%process(1)%ggd( time_sind )%ion( is + 1 )% &
                &   state(1)%label, .false.)
            call species( is, radiation%process(2)%ggd( time_sind )%ion( is + 1 )% &
                &   state(1)%label, .false.)
#endif
            ! Set (previous) label
            ion_label = edge_profiles%ggd( time_sind )%ion( is + 1 )%state(1)%label(1)
            ! Trim label (remove whitespaces on the right side)
            ion_label_tlen = len_trim(ion_label)
            ! Set default values for variables marking the position of '+' and '0'
            p = 0
            o = 0
            ! Loop through characters in ion label string
            do i = 1, ion_label_tlen
                if (ion_label(i:i) .eq. '0') then
                    ! Set the variable, marking the position of '0' if present
                    o = i
                endif
                ! When '+' is found, remember the position (set new variable)
                if (ion_label(i:i) .eq. "+") then
                    ! Set the variable, marking the position of '+' if present
                    p = i
                endif
            enddo
            ! Find charge
            ! (and ion atom label and position of '+' if found - commented out)
            if (p > 0 .and. (p > 0 .or. o > 0)) then
                !! ion = ion_label(1:p-1)
                !! plus = ion_label(p:p)
                ion_charge = ion_label(p+1:)
            else if (o > 0) then
                !! ion = ion_label(1:o-1)
                ion_charge = ion_label(o:)
            endif
            ! Convert charge from string to integer
            read(ion_charge, *) ion_charge_int

            ! Put (complete) ion label identifying the species
            edge_profiles%ggd( time_sind )%ion( is + 1 )%label(1) = ion_label

            ! Put ion charge
            edge_profiles%ggd( time_sind )%ion( is + 1 )%z_ion = ion_charge_int

            ! Put mass of ion
            edge_profiles%ggd( time_sind )%ion( is + 1 )%element(1)%a =     &
                &   am( is )

            ! Put nuclear charge
            edge_profiles%ggd( time_sind )%ion( is + 1 )%element(1)%z_n =   &
                &   zn( is )

            ! Put number of atoms
            edge_profiles%ggd( time_sind )%ion( is + 1 )%element(1)%atoms_n = 1

            ! Put minimum Z of the charge state bundle
            ! (z_min = z_max = 0 for a neutral)
            edge_profiles%ggd( time_sind )%ion( is + 1 )%state(1)%z_min =   &
                &   zamin( is )

            ! Put maximum Z of the charge state bundle
            edge_profiles%ggd( time_sind )%ion( is + 1 )%state(1)%z_max =   &
                &   zamax( is )

#if IMAS_MINOR_VERSION > 21
            ! Put (complete) ion label identifying the species
            radiation%process(1)%ggd( time_sind )%ion( is + 1 )%label(1) = &
                &   ion_label
            radiation%process(2)%ggd( time_sind )%ion( is + 1 )%label(1) = &
                &   ion_label

            ! Put ion charge
            radiation%process(1)%ggd( time_sind )%ion( is + 1 )%z_ion = &
                &   ion_charge_int
            radiation%process(2)%ggd( time_sind )%ion( is + 1 )%z_ion = &
                &   ion_charge_int

            ! Put mass of ion
            radiation%process(1)%ggd( time_sind )%ion( is + 1 )%element(1)%a = &
                &   am( is )
            radiation%process(2)%ggd( time_sind )%ion( is + 1 )%element(1)%a = &
                &   am( is )

            ! Put nuclear charge
            radiation%process(1)%ggd( time_sind )%ion( is + 1 )%element(1)%z_n =   &
                &   zn( is )
            radiation%process(2)%ggd( time_sind )%ion( is + 1 )%element(1)%z_n =   &
                &   zn( is )

            ! Put number of atoms
            radiation%process(1)%ggd( time_sind )%ion( is + 1 )%element(1)%atoms_n = 1
            radiation%process(2)%ggd( time_sind )%ion( is + 1 )%element(1)%atoms_n = 1

            ! Put neutral index
            radiation%process(1)%ggd( time_sind )%ion( is + 1 )%neutral_index = &
                &   b2eatcr(is)
            radiation%process(2)%ggd( time_sind )%ion( is + 1 )%neutral_index = &
                &   b2eatcr(is)

            ! Put multiple states flag
            radiation%process(1)%ggd( time_sind )%ion( is + 1 )%multiple_states_flag = 0
            radiation%process(2)%ggd( time_sind )%ion( is + 1 )%multiple_states_flag = 0

            ! Put minimum Z of the charge state bundle
            ! (z_min = z_max = 0 for a neutral)
            radiation%process(1)%ggd( time_sind )%ion( is + 1 )%state(1)%z_min =   &
                &   zamin( is )
            radiation%process(2)%ggd( time_sind )%ion( is + 1 )%state(1)%z_min =   &
                &   zamin( is )

            ! Put maximum Z of the charge state bundle
            radiation%process(1)%ggd( time_sind )%ion( is + 1 )%state(1)%z_max =   &
                &   zamax( is )
            radiation%process(2)%ggd( time_sind )%ion( is + 1 )%state(1)%z_max =   &
                &   zamax( is )
#endif

        enddo

#if IMAS_MINOR_VERSION > 21
#ifdef B25_EIRENE
        if (use_eirene.ne.0) then
          allocate( radiation%process(3)%ggd( time_sind )%neutral( natmi + nmoli ) )

          !! List of Eirene atoms
          do is = 1, natmi
            allocate( radiation%process(3)%ggd( time_sind )%neutral( is )%element(1) )
            allocate( radiation%process(3)%ggd( time_sind )%neutral( is )%label(1) )
            allocate( radiation%process(3)%ggd( time_sind )%neutral( is )%state(1) )
            allocate( radiation%process(3)%ggd( time_sind )%neutral( is )%state(1)%label(1) )
            radiation%process(3)%ggd( time_sind )%neutral( is )%element(1)%a = nmassa( is )
            radiation%process(3)%ggd( time_sind )%neutral( is )%element(1)%z_n = nchara( is )
            radiation%process(3)%ggd( time_sind )%neutral( is )%element(1)%atoms_n = 1
            radiation%process(3)%ggd( time_sind )%neutral( is )%label(1) = &
                &    textan( is-1 )
            radiation%process(3)%ggd( time_sind )%neutral( is )%ion_index = eb2atcr( is ) + 1
            radiation%process(3)%ggd( time_sind )%neutral( is )%multiple_states_flag = 0
            radiation%process(3)%ggd( time_sind )%neutral( is )%state(1)%label(1) = &
                &    textan( is-1 )
          end do

          !! List of molecules
          do j = 1, nmoli
            is = natmi + j
            nelems = count ( mlcmp( 1:natmi, j ) > 0 )
            allocate( radiation%process(3)%ggd( time_sind )%neutral( is )%element( nelems ) )
            allocate( radiation%process(3)%ggd( time_sind )%neutral( is )%label(1) )
            allocate( radiation%process(3)%ggd( time_sind )%neutral( is )%state(1) )
            allocate( radiation%process(3)%ggd( time_sind )%neutral( is )%state(1)%label(1) )
            i = 0
            do k = 1, natmi
              if (mlcmp( k, j ) > 0 ) then
                i = i + 1
                radiation%process(3)%ggd( time_sind )%neutral( is )%element(i)%a = nmassa( k )
                radiation%process(3)%ggd( time_sind )%neutral( is )%element(i)%z_n = nchara( k )
                radiation%process(3)%ggd( time_sind )%neutral( is )%element(i)%atoms_n = mlcmp(k,j)
              end if
            end do
            radiation%process(3)%ggd( time_sind )%neutral( is )%label(1) = &
                &    textmn( j-1 )
            radiation%process(3)%ggd( time_sind )%neutral( is )%ion_index = &
                &    eb2atcr( lmolscl(j) ) + 1
            radiation%process(3)%ggd( time_sind )%neutral( is )%multiple_states_flag = 0
            radiation%process(3)%ggd( time_sind )%neutral( is )%state(1)%label(1) = &
                &    textmn( j-1 )
          end do

          !! List of molecular ions
          allocate( radiation%process(4)%ggd( time_sind )%ion( nioni ) )
          do is = 1, nioni
            allocate( radiation%process(4)%ggd( time_sind )%ion( is )%label(1) )
            allocate( radiation%process(4)%ggd( time_sind )%ion( is )%state(1) )
            allocate( radiation%process(4)%ggd( time_sind )%ion( is )%state(1)%label(1) )
            allocate( radiation%process(4)%ggd( time_sind )%ion( is )%element(1) )

            radiation%process(4)%ggd( time_sind )%ion( is )%state(1)%label(1) = textin( is-1 )
            ion_label = adjustl(textin( is-1 ))
            j = 1
            match_found = .false.
            do while (.not.match_found.and.j.le.nmoli)
              mol_label = textmn(j-1)
              mol_label = trim(adjustl(mol_label))//'+'
              if (streql(ion_label, mol_label)) then
                match_found = .true.
              end if
              if (.not.match_found) j = j+1
            end do
            nelems = count ( mlcmp( 1:natmi, j ) > 0 )
            allocate( radiation%process(2)%ggd( time_sind )%ion( is )%element( nelems ) )
            i = 0
            do k = 1, natmi
              if ( mlcmp( k, j ) > 0 ) then
                i = i + 1
                radiation%process(4)%ggd( time_sind )%ion( is )%element( i )%a = nmassa(k)
                radiation%process(4)%ggd( time_sind )%ion( is )%element( i )%z_n = nchara(k)
                radiation%process(4)%ggd( time_sind )%ion( is )%element( i )%atoms_n = mlcmp(k,j)
              end if
            end do
            radiation%process(4)%ggd( time_sind )%ion( is )%z_ion = nchrgi( is )
            radiation%process(4)%ggd( time_sind )%ion( is )%label = textin( is-1 )
            radiation%process(4)%ggd( time_sind )%ion( is )%neutral_index = lkindi( is )
            radiation%process(4)%ggd( time_sind )%ion( is )%multiple_states_flag = 0
            radiation%process(4)%ggd( time_sind )%ion( is )%state(1)%z_min = nchari( is )
            radiation%process(4)%ggd( time_sind )%ion( is )%state(1)%z_max = nchari( is )
          end do

        end if
#endif
#endif

        geometryType = geometryId(nnreg, isymm, periodic_bc, topcut)
#if IMAS_MINOR_VERSION > 21
        allocate( summary%configuration%value(1) )
        summary%configuration%value = geometryName( geometryType )
#endif

        !! Write plasma state
        if ( B2_WRITE_DATA ) then
            call logmsg( LOGDEBUG, &
            &   "b2mod_ual_io.B25_process_ids: writing plasma state" )

            !! Find grid subset base indices out of the available grid subset
            !! data stored in the IDS. That is done using IMAS GGD routine
            !! findGridSubsetByName().
#if IMAS_MINOR_VERSION < 15
            iGsCoreBoundary = findGridSubsetByName(         &
                &   edge_profiles%ggd( time_sind )%grid,    &
                &   "Core boundary" )
            iGsInnerMidplane = findGridSubsetByName( edge_profiles% &
                &   ggd( time_sind )%grid, "Inner Midplane" )
            iGsOuterMidplane = findGridSubsetByName( edge_profiles% &
                &   ggd( time_sind )%grid, "Outer Midplane" )
            iGsCore = findGridSubsetByName( edge_profiles%      &
                &   ggd( time_sind )%grid, "Core" )
            iGsSOL = findGridSubsetByName( edge_profiles%       &
                &   ggd( time_sind )%grid, "SOL" )
            iGsIDivertor = findGridSubsetByName( edge_profiles% &
                &   ggd( time_sind )%grid, "Inner divertor" )
            iGsODivertor = findGridSubsetByName( edge_profiles% &
                &   ggd( time_sind )%grid, "Outer divertor" )
#else
            iGsCoreBoundary = findGridSubsetByName(   &
                &   edge_profiles%grid_ggd( time_sind ), "Core boundary" )
            iGsInnerMidplane = findGridSubsetByName(  &
                &   edge_profiles%grid_ggd( time_sind ), "Inner Midplane" )
            iGsOuterMidplane = findGridSubsetByName(  &
                &   edge_profiles%grid_ggd( time_sind ), "Outer Midplane" )
            iGsCore = findGridSubsetByName(           &
                &   edge_profiles%grid_ggd( time_sind ), "Core" )
            iGsSOL = findGridSubsetByName(            &
                &   edge_profiles%grid_ggd( time_sind ), "SOL" )
            iGsIDivertor = findGridSubsetByName(      &
                &   edge_profiles%grid_ggd( time_sind ), "Inner divertor" )
            iGsODivertor = findGridSubsetByName(      &
                &   edge_profiles%grid_ggd( time_sind ), "Outer divertor" )
#endif
            !! TODO: The fluxes are currently in the edge_transport. They are
            !! supposed to be in the edge_profiles (24.10.2017)
            !! use electrons%particles or electrons%energy node?
            !! Currently using electrons%energy node

            !! edge_transport fundamentals
            !! TODO: create it as a subroutine
            allocate( edge_transport%model(1) )
            allocate( edge_transport%model(1)%ggd( num_time_slices ) )

            !! ne: Electron Density
            call write_quantity(                                            &
                &   val = edge_profiles%ggd( time_sind )%electrons%density, &
                &   fluxes = edge_transport%model(1)%ggd( time_sind )%      &
                &            electrons%particles%flux,                      &
                &   value = ne,                                             &
                &   flux = fne,                                             &
                &   time_sind = time_sind )
            call write_cell_scalar(                                         &
                &   scalar = edge_sources%source(1)%ggd( time_sind )%       &
                &            electrons%particles,                           &
                &   b2CellData = sne(:,:,0) + sne(:,:,1) * ne )

            !! ni (SOLPS 4.x) /
            !! na (SOLPS 5.x): Ion Density
            allocate( edge_transport%model(1)%ggd( time_sind )%ion( ns ) )
            allocate( edge_sources%source(1)%ggd( time_sind )%ion( ns ) )

            do is = 1, ns
                call write_quantity(                                          &
                    &   val = edge_profiles%ggd( time_sind )%ion(is)%density, &
                    &   fluxes = edge_transport%model(1)%ggd( time_sind )%    &
                    &   ion( is )%particles%flux,                             &
                    &   value = na(:,:, is - 1 ),                             &
                    &   flux = fna(:,:,:, is - 1 ),                           &
                    &   time_sind = time_sind )
                call write_cell_scalar( scalar = edge_sources%source(1)%      &
                    &   ggd( time_sind )%ion( is )%particles,                 &
                    &   b2CellData =                                          &
                    &   sna(:,:,0, is - 1 ) +                                 &
                    &   sna(:,:,1, is - 1 )*na(:,:, is - 1 ) )
            end do

!!$    !! ue: Parallel Electron Velocity
!!$    TODO: must be computed, refactor code from b2news into function
!!$    allocate(edge_profiles%fluid%ve%comps(1))
!!$    allocate(edge_profiles%fluid%ve%align(1))
!!$    allocate(edge_profiles%fluid%ve%alignid(1))
!!$    edge_profiles%fluid%ve%align(1) = VEC_ALIGN_PARALLEL
!!$    edge_profiles%fluid%ve%alignid(1) = VEC_ALIGN_PARALLEL_ID
!!$    call write_cell_scalar( edge_profiles%fluid%ve%comps(1)%value, ue(:,:) )
            !!$ call write_cell_vector_component(                           &
            !!$     &   vectorComponent = edge_profiles%ggd( time_sind )%   &
            !!$     &                     electrons%velocity,               &
            !!$     &   b2CellData = ue(:,:),                               &
            !!$     &   vectorID = VEC_ALIGN_PARALLEL_ID )

            !! ua: Parallel Ion Velocity
            do is = 1, ns
                allocate( edge_profiles%ggd( time_sind )%ion( is )%velocity(1) )

                call write_cell_vector_component(                           &
                    &   vectorComponent = edge_profiles%ggd( time_sind )%   &
                    &                     ion( is )%velocity,               &
                    &   b2CellData = ua(:,:, is - 1 ),                      &
                    &   vectorID = VEC_ALIGN_PARALLEL_ID )
            end do

            !! te: Electron Temperature
            call write_quantity(                                        &
                &   val = edge_profiles%ggd( time_sind )%electrons%     &
                &         temperature,                                  &
                &   fluxes = edge_transport%model(1)%ggd( time_sind )%  &
                &            electrons%energy%flux,                     &
                &   value = te/qe,                                      &
                &   flux = fhe,                                         &
                &   time_sind = time_sind )

            !! ti: Ion Temperature
            allocate( edge_profiles%ggd( time_sind )%ion(1)%temperature(1) )
            allocate( edge_transport%model(1)%ggd( time_sind )%ion( 1 ) )
            call write_quantity(                                            &
                &   val = edge_profiles%ggd( time_sind )%ion(1)%temperature,&
                &   fluxes = edge_transport%model(1)% ggd( time_sind )%     &
                &   ion(1)%energy%flux,                                     &
                &   value = ti/qe,                                          &
                &   flux = fhi,                                             &
                &   time_sind = time_sind )

            !! po: Electric Potential
            call write_cell_scalar( edge_profiles%ggd( time_sind )% &
                &   phi_potential, po )

            !! B (magnetic field vector)
            !! Compute unit basis vectors along the field directions
            call compute_Coordinate_Unit_Vectors(crx, cry, e(:,:,:,1), &
                &   e(:,:,:,2), e(:,:,:,3))

            !! Write the three unit basis vectors
            call write_cell_vector_component(                                 &
                &   vectorComponent = edge_profiles%ggd( time_sind )%e_field, &
                &   b2CellData = e(:,:,:,1),                                  &
                &   vectorID = VEC_ALIGN_POLOIDAL_ID )

            call write_cell_vector_component(                                 &
                &   vectorComponent = edge_profiles%ggd( time_sind )%e_field, &
                &   b2CellData = e(:,:,:,2),                                  &
                &   vectorID = VEC_ALIGN_RADIAL_ID )

            call write_cell_vector_component(                                 &
                &   vectorComponent = edge_profiles%ggd( time_sind )%e_field, &
                &   b2CellData = e(:,:,:,3),                                  &
                &   vectorID = VEC_ALIGN_TOROIDAL_ID )

            !! write the magnetic field vector in the b2 coordinate system
            call write_cell_vector_component(                                 &
                &   vectorComponent = edge_profiles%ggd( time_sind )%e_field, &
                &   b2CellData = bb(:,:,0:2),                                 &
                &   vectorID = "diamagnetic" )

#if IMAS_MINOR_VERSION > 21
            !! write the emissivity data
            !! Process 1. Line and recombination radiation from B2.5 ions
            do is = 0, ns-1
                call write_cell_scalar( scalar = radiation%process(1)%        &
                    &   ggd( time_sind )%ion( is + 1 )%emissivity,            &
                    &   b2CellData = rqrad(:,:,is) )
            end do
            !! Process 2. Bremsstrahlung from B2.5 ions
            do is = 0, ns-1
                call write_cell_scalar( scalar = radiation%process(2)%        &
                    &   ggd( time_sind )%ion( is + 1 )%emissivity,            &
                    &   b2CellData = rqbrm(:,:,is) )
            end do
#ifdef B25_EIRENE
            if (use_eirene.ne.0) then

              !! Process 3. Eirene neutrals (atoms and molecules)
              do is = 1, natmi
                do ix = -1, nx
                  do iy = -1, ny
                    tmpCv(ix,iy) = -1.0_R8*eneutrad(ix+1,iy+1,is,0)
                  end do
                end do
                call write_cell_scalar( scalar = radiation%process(3)%        &
                    &   ggd( time_sind )%neutral( is )%emissivity,            &
                    &   b2CellData = tmpCv(:,:) )
              end do
              do is = 1, nmoli
                do ix = -1, nx
                  do iy = -1, ny
                    tmpCv(ix,iy) = -1.0_R8*emolrad(ix+1,iy+1,is,0)
                  end do
                end do
                call write_cell_scalar( scalar = radiation%process(3)%        &
                    &   ggd( time_sind )%neutral( is + natmi )%emissivity,    &
                    &   b2CellData = tmpCv(:,:) )
              end do
              !! Process 4. Eirene molecular ions
              do is = 1, nioni
                do ix = -1, nx
                  do iy = -1, ny
                    tmpCv(ix,iy) = -1.0_R8*eionrad(ix+1,iy+1,is,0)
                  end do
                end do
                call write_cell_scalar( scalar = radiation%process(4)%        &
                    &   ggd( time_sind )%ion( is )%emissivity,                &
                    &   b2CellData = tmpCv(:,:) )
              end do
            end if
#endif
#endif
        end if

#if IMAS_MINOR_VERSION > 21
! Summary separatrix data
        if (maxval(abs(fpsi(-1:nx,-1:ny,0:3))).gt.0.0_R8) then
           allocate( summary%local%separatrix%position%psi( time_sind ) )
           summary%local%separatrix%position%psi( time_sind ) = fpsi(jxa,jsep,2)
        end if
        allocate( summary%local%separatrix%t_e%value( time_sind ) )
        summary%local%separatrix%t_e%value( time_sind ) = &
           &  0.5_R8 * (te(jxa,jsep)+ te(topix(jxa,jsep),topiy(jxa,jsep)))/ev
        allocate ( summary%local%separatrix%t_e%source(1) )
        summary%local%separatrix%t_e%source = source
        allocate( summary%local%separatrix%t_i_average%value( time_sind ) )
        summary%local%separatrix%t_i_average%value( time_sind ) = &
           &  0.5_R8 * (ti(jxa,jsep)+ ti(topix(jxa,jsep),topiy(jxa,jsep)))/ev
        allocate ( summary%local%separatrix%t_i_average%source(1) )
        summary%local%separatrix%t_i_average%source = source
        allocate( summary%local%separatrix%n_e%value( time_sind ) )
        summary%local%separatrix%n_e%value( time_sind ) = &
           &  0.5_R8 * (ne(jxa,jsep)+ ne(topix(jxa,jsep),topiy(jxa,jsep)))
        allocate ( summary%local%separatrix%n_e%source(1) )
        summary%local%separatrix%n_e%source = source
        do is = 1, nspecies
          is1 = eb2spcr(is)
          if (nint(zamax(is1)).eq.0) is1 = is1+1
          is2 = is1 + nfluids(is) - 1
          nisep = 0.0_R8
          nasum = 0.0_R8
          vtor = 0.0_R8
          do i = is1, is2
            nisep = nisep + &
              &  0.5_R8 * (na(jxa,jsep,i) + na(topix(jxa,jsep),topiy(jxa,jsep),i))
            nasum = nasum + na(jxa,jsep,i)
            vtor = vtor + ua(jxa,jsep,i)*abs(bb(jxa,jsep,2)/bb(jxa,jsep,3)) &
              &  * na(jxa,jsep,i)
          end do
          if (nasum.gt.0.0_R8) vtor = vtor / nasum
          select case (is_codes(eb2spcr(is)))
          case ('H')
            allocate( summary%local%separatrix%n_i%hydrogen%value( time_sind ))
            summary%local%separatrix%n_i%hydrogen%value( time_sind ) = nisep
            allocate ( summary%local%separatrix%n_i%hydrogen%source(1) )
            summary%local%separatrix%n_i%hydrogen%source = source
            allocate( summary%local%separatrix%velocity_tor%hydrogen%value( time_sind ))
            summary%local%separatrix%velocity_tor%hydrogen%value( time_sind ) = vtor
            allocate ( summary%local%separatrix%velocity_tor%hydrogen%source(1) )
            summary%local%separatrix%velocity_tor%hydrogen%source = source
          case ('D')
            allocate( summary%local%separatrix%n_i%deuterium%value( time_sind ))
            summary%local%separatrix%n_i%deuterium%value( time_sind ) = nisep
            allocate ( summary%local%separatrix%n_i%deuterium%source(1) )
            summary%local%separatrix%n_i%deuterium%source = source
            allocate( summary%local%separatrix%velocity_tor%deuterium%value( time_sind ))
            summary%local%separatrix%velocity_tor%deuterium%value( time_sind ) = vtor
            allocate ( summary%local%separatrix%velocity_tor%deuterium%source(1) )
            summary%local%separatrix%velocity_tor%deuterium%source = source
          case ('T')
            allocate( summary%local%separatrix%n_i%tritium%value( time_sind ))
            summary%local%separatrix%n_i%tritium%value( time_sind ) = nisep
            allocate ( summary%local%separatrix%n_i%tritium%source(1) )
            summary%local%separatrix%n_i%tritium%source = source
            allocate( summary%local%separatrix%velocity_tor%tritium%value( time_sind ))
            summary%local%separatrix%velocity_tor%tritium%value( time_sind ) = vtor
            allocate ( summary%local%separatrix%velocity_tor%tritium%source(1) )
            summary%local%separatrix%velocity_tor%tritium%source = source
          case ('He')
            if (nint(am(is)).eq.3) then
              allocate( summary%local%separatrix%n_i%helium_3%value( time_sind ))
              summary%local%separatrix%n_i%helium_3%value( time_sind ) = nisep
              allocate ( summary%local%separatrix%n_i%helium_3%source(1) )
              summary%local%separatrix%n_i%helium_3%source = source
              allocate( summary%local%separatrix%velocity_tor%helium_3%value( time_sind ))
              summary%local%separatrix%velocity_tor%helium_3%value( time_sind ) = vtor
              allocate ( summary%local%separatrix%velocity_tor%helium_3%source(1) )
              summary%local%separatrix%velocity_tor%helium_3%source = source
            else if (nint(am(is)).eq.4) then
              allocate( summary%local%separatrix%n_i%helium_4%value( time_sind ))
              summary%local%separatrix%n_i%helium_4%value( time_sind ) = nisep
              allocate ( summary%local%separatrix%n_i%helium_4%source(1) )
              summary%local%separatrix%n_i%helium_4%source = source
              allocate( summary%local%separatrix%velocity_tor%helium_4%value( time_sind ))
              summary%local%separatrix%velocity_tor%helium_4%value( time_sind ) = vtor
              allocate ( summary%local%separatrix%velocity_tor%helium_4%source(1) )
              summary%local%separatrix%velocity_tor%helium_4%source = source
            end if
          case ('Li')
            allocate( summary%local%separatrix%n_i%lithium%value( time_sind ))
            summary%local%separatrix%n_i%lithium%value( time_sind ) = nisep
            allocate ( summary%local%separatrix%n_i%lithium%source(1) )
            summary%local%separatrix%n_i%lithium%source = source
            allocate( summary%local%separatrix%velocity_tor%lithium%value( time_sind ))
            summary%local%separatrix%velocity_tor%lithium%value( time_sind ) = vtor
            allocate ( summary%local%separatrix%velocity_tor%lithium%source(1) )
            summary%local%separatrix%velocity_tor%lithium%source = source
          case ('Be')
            allocate( summary%local%separatrix%n_i%beryllium%value( time_sind ))
            summary%local%separatrix%n_i%beryllium%value( time_sind ) = nisep
            allocate ( summary%local%separatrix%n_i%beryllium%source(1) )
            summary%local%separatrix%n_i%beryllium%source = source
            allocate( summary%local%separatrix%velocity_tor%beryllium%value( time_sind ))
            summary%local%separatrix%velocity_tor%beryllium%value( time_sind ) = vtor
            allocate ( summary%local%separatrix%velocity_tor%beryllium%source(1) )
            summary%local%separatrix%velocity_tor%beryllium%source = source
          case ('C')
            allocate( summary%local%separatrix%n_i%carbon%value( time_sind ))
            summary%local%separatrix%n_i%carbon%value( time_sind ) = nisep
            allocate ( summary%local%separatrix%n_i%carbon%source(1) )
            summary%local%separatrix%n_i%carbon%source = source
            allocate( summary%local%separatrix%velocity_tor%carbon%value( time_sind ))
            summary%local%separatrix%velocity_tor%carbon%value( time_sind ) = vtor
            allocate ( summary%local%separatrix%velocity_tor%carbon%source(1) )
            summary%local%separatrix%velocity_tor%carbon%source = source
          case ('N')
            allocate( summary%local%separatrix%n_i%nitrogen%value( time_sind ))
            summary%local%separatrix%n_i%nitrogen%value( time_sind ) = nisep
            allocate ( summary%local%separatrix%n_i%nitrogen%source(1) )
            summary%local%separatrix%n_i%nitrogen%source = source
            allocate( summary%local%separatrix%velocity_tor%nitrogen%value( time_sind ))
            summary%local%separatrix%velocity_tor%nitrogen%value( time_sind ) = vtor
            allocate ( summary%local%separatrix%velocity_tor%nitrogen%source(1) )
            summary%local%separatrix%velocity_tor%nitrogen%source = source
          case ('O')
            allocate( summary%local%separatrix%n_i%oxygen%value( time_sind ))
            summary%local%separatrix%n_i%oxygen%value( time_sind ) = nisep
            allocate ( summary%local%separatrix%n_i%oxygen%source(1) )
            summary%local%separatrix%n_i%oxygen%source = source
            allocate( summary%local%separatrix%velocity_tor%oxygen%value( time_sind ))
            summary%local%separatrix%velocity_tor%oxygen%value( time_sind ) = vtor
            allocate ( summary%local%separatrix%velocity_tor%oxygen%source(1) )
            summary%local%separatrix%velocity_tor%oxygen%source = source
          case ('Ne')
            allocate( summary%local%separatrix%n_i%neon%value( time_sind ))
            summary%local%separatrix%n_i%neon%value( time_sind ) = nisep
            allocate ( summary%local%separatrix%n_i%neon%source(1) )
            summary%local%separatrix%n_i%neon%source = source
            allocate( summary%local%separatrix%velocity_tor%neon%value( time_sind ))
            summary%local%separatrix%velocity_tor%neon%value( time_sind ) = vtor
            allocate ( summary%local%separatrix%velocity_tor%neon%source(1) )
            summary%local%separatrix%velocity_tor%neon%source = source
          case ('Ar')
            allocate( summary%local%separatrix%n_i%argon%value( time_sind ))
            summary%local%separatrix%n_i%argon%value( time_sind ) = nisep
            allocate ( summary%local%separatrix%n_i%argon%source(1) )
            summary%local%separatrix%n_i%argon%source = source
            allocate( summary%local%separatrix%velocity_tor%argon%value( time_sind ))
            summary%local%separatrix%velocity_tor%argon%value( time_sind ) = vtor
            allocate ( summary%local%separatrix%velocity_tor%argon%source(1) )
            summary%local%separatrix%velocity_tor%argon%source = source
          case ('Xe')
            allocate( summary%local%separatrix%n_i%xenon%value( time_sind ))
            summary%local%separatrix%n_i%xenon%value( time_sind ) = nisep
            allocate ( summary%local%separatrix%n_i%xenon%source(1) )
            summary%local%separatrix%n_i%xenon%source = source
            allocate( summary%local%separatrix%velocity_tor%xenon%value( time_sind ))
            summary%local%separatrix%velocity_tor%xenon%value( time_sind ) = vtor
            allocate ( summary%local%separatrix%velocity_tor%xenon%source(1) )
            summary%local%separatrix%velocity_tor%xenon%source = source
          case ('W')
            allocate( summary%local%separatrix%n_i%tungsten%value( time_sind ))
            summary%local%separatrix%n_i%tungsten%value( time_sind ) = nisep
            allocate ( summary%local%separatrix%n_i%tungsten%source(1) )
            summary%local%separatrix%n_i%tungsten%source = source
            allocate( summary%local%separatrix%velocity_tor%tungsten%value( time_sind ))
            summary%local%separatrix%velocity_tor%tungsten%value( time_sind ) = vtor
            allocate ( summary%local%separatrix%velocity_tor%tungsten%source(1) )
            summary%local%separatrix%velocity_tor%tungsten%source = source
          end select
        end do
        allocate( summary%local%separatrix%n_i_total%value( time_sind ))
        summary%local%separatrix%n_i_total%value( time_sind ) = &
          & 0.5_R8 * (ni(jxa,jsep,1) + ni(topix(jxa,jsep),topiy(jxa,jsep),1))
        allocate ( summary%local%separatrix%n_i_total%source(1) )
        summary%local%separatrix%n_i_total%source = source
        allocate( summary%local%separatrix%zeff%value( time_sind ))
        summary%local%separatrix%zeff%value( time_sind ) = &
          & 0.5_R8 * (zeff(jxa,jsep) + zeff(topix(jxa,jsep),topiy(jxa,jsep)))
        allocate ( summary%local%separatrix%zeff%source(1) )
        summary%local%separatrix%zeff%source = source

! Summary divertor plate data
        if (nncut.eq.0) then
          if (geometryType.eq.GEOMETRY_LINEAR) then
            ntrgts=0
            if (use_eirene.ne.0) then
              do while (ltns(ntrgts+1).gt.0)
                ntrgts=ntrgts+1
              end do
            else if (boundary_namelist.ne.0) then
              do i=1,nbc
                if(bcchar(i).eq.'E'.or.bcchar(i).eq.'W') then
                  if(bcene(i).eq. 3.or.bcene(i).eq.12.or.bcene(i).eq.15) then
                    ntrgts=ntrgts+1
                    plate_name(ntrgts) = bcchar(i)
                    ixpos(ntrgts) = bcpos(i)
                    if(bcchar(i).eq.'W') then
                      ixpos(ntrgts) = ixpos(ntrgts)+target_offset
                    else if (bcchar(i).eq.'E') then
                      ixpos(ntrgts) = ixpos(ntrgts)-target_offset
                    end if
                    iypos(ntrgts) = jsep
                  end if
                end if
              end do
            end if
          else
            ntrgts = 0
          end if
        else
          ntrgts = 2*nncut
          plate_name(1) = 'LI'
          ixpos(1) = -1+target_offset
          iypos(1) = topcut(1)
          if (nncut.eq.2) then
            plate_name(2) = 'UI'
            ixpos(2) = nxtl-target_offset
            iypos(2) = topcut(2)
            plate_name(3) = 'UO'
            ixpos(3) = nxtr+target_offset
            iypos(3) = topcut(2)
          end if
          plate_name(ntrgts) = 'LO'
          ixpos(ntrgts) = nx-target_offset
          iypos(ntrgts) = topcut(1)
        end if
        if (ntrgts.gt.0) then
          allocate ( summary%local%divertor_plate( ntrgts ) )
          do i = 1, ntrgts
            allocate( summary%local%divertor_plate(i)%name%value( time_sind ) )
            summary%local%divertor_plate(i)%name%value( time_sind ) = plate_name(i)
            allocate( summary%local%divertor_plate(i)%name%source(1) )
            summary%local%divertor_plate(i)%name%source(1) = source
            allocate( summary%local%divertor_plate(i)%t_e%value( time_sind ) )
            summary%local%divertor_plate(i)%t_e%value( time_sind ) = &
              &  0.5_R8 * (te(ixpos(i),iypos(i))+  &
              &            te(topix(ixpos(i),iypos(i)),topiy(ixpos(i),iypos(i))))/ev
            allocate( summary%local%divertor_plate(i)%t_e%source(1) )
            summary%local%divertor_plate(i)%t_e%source(1) = source
            allocate( summary%local%divertor_plate(i)%t_i_average%value( time_sind ) )
            summary%local%divertor_plate(i)%t_i_average%value( time_sind ) = &
              &  0.5_R8 * (te(ixpos(i),iypos(i))+  &
              &            te(topix(ixpos(i),iypos(i)),topiy(ixpos(i),iypos(i))))/ev
            allocate( summary%local%divertor_plate(i)%t_i_average%source(1) )
            summary%local%divertor_plate(i)%t_i_average%source(1) = source
            allocate( summary%local%divertor_plate(i)%n_e%value( time_sind ) )
            summary%local%divertor_plate(i)%n_e%value( time_sind ) = &
              &  0.5_R8 * (ne(ixpos(i),iypos(i))+  &
              &            ne(topix(ixpos(i),iypos(i)),topiy(ixpos(i),iypos(i))))
            allocate( summary%local%divertor_plate(i)%n_e%source(1) )
            summary%local%divertor_plate(i)%n_e%source(1) = source
            do is = 1, nspecies
              is1 = eb2spcr(is)
              if (nint(zamax(is1)).eq.0) is1 = is1+1
              is2 = is1 + nfluids(is) - 1
              nisep = 0.0_R8
              do j = is1, is2
                nisep = nisep + &
                  &  0.5_R8 * (na(ixpos(i),iypos(i),j) + &
                  &            na(topix(ixpos(i),iypos(i)),topiy(ixpos(i),iypos(i)),j))
              end do
              select case (is_codes(eb2spcr(is)))
              case ('H')
                allocate( summary%local%divertor_plate(i)%n_i%hydrogen%value( time_sind ))
                summary%local%divertor_plate(i)%n_i%hydrogen%value( time_sind ) = nisep
                allocate ( summary%local%divertor_plate(i)%n_i%hydrogen%source(1) )
                summary%local%divertor_plate(i)%n_i%hydrogen%source = source
              case ('D')
                allocate( summary%local%divertor_plate(i)%n_i%deuterium%value( time_sind ))
                summary%local%divertor_plate(i)%n_i%deuterium%value( time_sind ) = nisep
                allocate ( summary%local%divertor_plate(i)%n_i%deuterium%source(1) )
                summary%local%divertor_plate(i)%n_i%deuterium%source = source
              case ('T')
                allocate( summary%local%divertor_plate(i)%n_i%tritium%value( time_sind ))
                summary%local%divertor_plate(i)%n_i%tritium%value( time_sind ) = nisep
                allocate ( summary%local%divertor_plate(i)%n_i%tritium%source(1) )
                summary%local%divertor_plate(i)%n_i%tritium%source = source
              case ('He')
                if (nint(am(is)).eq.3) then
                  allocate( summary%local%divertor_plate(i)%n_i%helium_3%value( time_sind ))
                  summary%local%divertor_plate(i)%n_i%helium_3%value( time_sind ) = nisep
                  allocate ( summary%local%divertor_plate(i)%n_i%helium_3%source(1) )
                  summary%local%divertor_plate(i)%n_i%helium_3%source = source
                else if (nint(am(is)).eq.4) then
                  allocate( summary%local%divertor_plate(i)%n_i%helium_4%value( time_sind ))
                  summary%local%divertor_plate(i)%n_i%helium_4%value( time_sind ) = nisep
                  allocate ( summary%local%divertor_plate(i)%n_i%helium_4%source(1) )
                  summary%local%divertor_plate(i)%n_i%helium_4%source = source
                end if
              case ('Li')
                allocate( summary%local%divertor_plate(i)%n_i%lithium%value( time_sind ))
                summary%local%divertor_plate(i)%n_i%lithium%value( time_sind ) = nisep
                allocate ( summary%local%divertor_plate(i)%n_i%lithium%source(1) )
                summary%local%divertor_plate(i)%n_i%lithium%source = source
              case ('Be')
                allocate( summary%local%divertor_plate(i)%n_i%beryllium%value( time_sind ))
                summary%local%divertor_plate(i)%n_i%beryllium%value( time_sind ) = nisep
                allocate ( summary%local%divertor_plate(i)%n_i%beryllium%source(1) )
                summary%local%divertor_plate(i)%n_i%beryllium%source = source
              case ('C')
                allocate( summary%local%divertor_plate(i)%n_i%carbon%value( time_sind ))
                summary%local%divertor_plate(i)%n_i%carbon%value( time_sind ) = nisep
                allocate ( summary%local%divertor_plate(i)%n_i%carbon%source(1) )
                summary%local%divertor_plate(i)%n_i%carbon%source = source
              case ('N')
                allocate( summary%local%divertor_plate(i)%n_i%nitrogen%value( time_sind ))
                summary%local%divertor_plate(i)%n_i%nitrogen%value( time_sind ) = nisep
                allocate ( summary%local%divertor_plate(i)%n_i%nitrogen%source(1) )
                summary%local%divertor_plate(i)%n_i%nitrogen%source = source
              case ('O')
                allocate( summary%local%divertor_plate(i)%n_i%oxygen%value( time_sind ))
                summary%local%divertor_plate(i)%n_i%oxygen%value( time_sind ) = nisep
                allocate ( summary%local%divertor_plate(i)%n_i%oxygen%source(1) )
                summary%local%divertor_plate(i)%n_i%oxygen%source = source
              case ('Ne')
                allocate( summary%local%divertor_plate(i)%n_i%neon%value( time_sind ))
                summary%local%divertor_plate(i)%n_i%neon%value( time_sind ) = nisep
                allocate ( summary%local%divertor_plate(i)%n_i%neon%source(1) )
                summary%local%divertor_plate(i)%n_i%neon%source = source
              case ('Ar')
                allocate( summary%local%divertor_plate(i)%n_i%argon%value( time_sind ))
                summary%local%divertor_plate(i)%n_i%argon%value( time_sind ) = nisep
                allocate ( summary%local%divertor_plate(i)%n_i%argon%source(1) )
                summary%local%divertor_plate(i)%n_i%argon%source = source
              case ('Xe')
                allocate( summary%local%divertor_plate(i)%n_i%xenon%value( time_sind ))
                summary%local%divertor_plate(i)%n_i%xenon%value( time_sind ) = nisep
                allocate ( summary%local%divertor_plate(i)%n_i%xenon%source(1) )
                summary%local%divertor_plate(i)%n_i%xenon%source = source
              case ('W')
                allocate( summary%local%divertor_plate(i)%n_i%tungsten%value( time_sind ))
                summary%local%divertor_plate(i)%n_i%tungsten%value( time_sind ) = nisep
                allocate ( summary%local%divertor_plate(i)%n_i%tungsten%source(1) )
                summary%local%divertor_plate(i)%n_i%tungsten%source = source
              end select
            end do
            allocate( summary%local%divertor_plate(i)%n_i_total%value( time_sind ))
            summary%local%divertor_plate(i)%n_i_total%value( time_sind ) = &
              & 0.5_R8 * (ni(ixpos(i),iypos(i),1) + &
              &           ni(topix(ixpos(i),iypos(i)),topiy(ixpos(i),iypos(i)),1))
            allocate ( summary%local%divertor_plate(i)%n_i_total%source(1) )
            summary%local%divertor_plate(i)%n_i_total%source = source
            allocate( summary%local%divertor_plate(i)%zeff%value( time_sind ))
            summary%local%divertor_plate(i)%zeff%value( time_sind ) = &
              & 0.5_R8 * (zeff(ixpos(i),iypos(i)) + &
              &           zeff(topix(ixpos(i),iypos(i)),topiy(ixpos(i),iypos(i))))
            allocate ( summary%local%divertor_plate(i)%zeff%source(1) )
            summary%local%divertor_plate(i)%zeff%source = source
          end do
        end if

        allocate( summary%boundary%type%value( time_sind ) )
        select case (GeometryType)
          case( GEOMETRY_LIMITER )
            summary%boundary%type%value( time_sind ) = 0
          case( GEOMETRY_SN )
            summary%boundary%type%value( time_sind ) = 1
            if ( isymm.eq.1 .or. isymm.eq.2 ) then
              if ( cry(leftcut(1),jsep,3).lt.0.0_R8) then
                summary%boundary%type%value( time_sind ) = 11
              else if ( cry(leftcut(1),jsep,3).gt.0.0_R8) then
                summary%boundary%type%value( time_sind ) = 12
              end if
            else if ( isymm.eq.3 .or. isymm.eq.4 ) then
              if ( crx(leftcut(1),jsep,3).lt.0.0_R8) then
                summary%boundary%type%value( time_sind ) = 11
              else if ( crx(leftcut(1),jsep,3).gt.0.0_R8) then
                summary%boundary%type%value( time_sind ) = 12
              end if
            end if
          case( GEOMETRY_CDN , GEOMETRY_DDN_BOTTOM , GEOMETRY_DDN_TOP )
            summary%boundary%type%value( time_sind ) = 13
        end select
        allocate( summary%boundary%type%source(1) )
        summary%boundary%type%source = source
        allocate( summary%boundary%strike_point_inner_r%value( time_sind ) )
        summary%boundary%strike_point_inner_r%value( time_sind ) = crx(-1,topcut(1),3)
        allocate( summary%boundary%strike_point_inner_r%source(1) )
        summary%boundary%strike_point_inner_r%source = source
        allocate( summary%boundary%strike_point_inner_z%value( time_sind ) )
        summary%boundary%strike_point_inner_z%value( time_sind ) = cry(-1,topcut(1),3)
        allocate( summary%boundary%strike_point_inner_z%source(1) )
        summary%boundary%strike_point_inner_z%source = source
        allocate( summary%boundary%strike_point_outer_r%value( time_sind ) )
        summary%boundary%strike_point_outer_r%value( time_sind ) = crx(nx,topcut(1),1)
        allocate( summary%boundary%strike_point_outer_r%source(1) )
        summary%boundary%strike_point_outer_r%source = source
        allocate( summary%boundary%strike_point_outer_z%value( time_sind ) )
        summary%boundary%strike_point_outer_z%value( time_sind ) = cry(nx,topcut(1),1)
        allocate( summary%boundary%strike_point_outer_z%source(1) )
        summary%boundary%strike_point_outer_z%source = source

        allocate( summary%fusion%power%value( time_sind ) )
        summary%fusion%power%value( time_sind ) = fusion_power
        allocate( summary%fusion%power%source(1) )
        summary%fusion%power%source = source

#ifdef B25_EIRENE
        call eirene_mc(nx, ny, ns, nxtl, nxtr, ismain, time_step, BoRiS, &
          &  sna0, smo0, she0, shi0, .false.)
        summary%gas_injection_rates%impurity_seeding%value = 0
        allocate( summary%gas_injection_rates%impurity_seeding%source(1) )
        summary%gas_injection_rates%impurity_seeding%source = source
        gsum = 0.0_R8
        gtop = 0.0_R8
        gmid = 0.0_R8
        gbot = 0.0_R8
        do istrai = 1, nstrai
          if (crcstra(istrai).eq.'C') then
            do iatm = 1, natmi
              gsum = gsum + tflux(istrai)*gpfc(iatm,istrai)*zn(eb2atcr(iatm))
              if (istrai.eq.nesepm_istra) then
                if (ndes.gt.0.0_R8 .or. nesepm_pfr.gt.0.0_R8 .or. &
                  & private_flux_puff.gt.0.0_R8) then
                  gbot = gbot + tflux(istrai)*gpfc(iatm,istrai)*zn(eb2atcr(iatm))
                else if (nesepm_sol.gt.0.0_R8 .or. volrec_sol.gt.0.0_R8 .or. &
                  & ndes_sol.gt.0.0_R8 .or. nepedm_sol.gt.0.0_R8) then
                  gmid = gmid + tflux(istrai)*gpfc(iatm,istrai)*zn(eb2atcr(iatm))
                end if
              end if
              if (gpfc(iatm,istrai).eq.1.0_R8) then
               select case (is_codes(eb2atcr(iatm)))
               case ('H')
                if ( associated( summary%gas_injection_rates%hydrogen%value ) ) then
                  summary%gas_injection_rates%hydrogen%value( time_sind ) = &
                    & summary%gas_injection_rates%hydrogen%value( time_sind ) + &
                    & tflux(istrai)*zn(eb2atcr(iatm))
                else
                  allocate( summary%gas_injection_rates%hydrogen%value( time_sind ))
                  summary%gas_injection_rates%hydrogen%value( time_sind ) = &
                    & tflux(istrai)*zn(eb2atcr(iatm))
                end if
                allocate ( summary%gas_injection_rates%hydrogen%source(1) )
                summary%gas_injection_rates%hydrogen%source = source
               case ('D')
                if ( associated( summary%gas_injection_rates%deuterium%value ) ) then
                  summary%gas_injection_rates%deuterium%value( time_sind ) = &
                    & summary%gas_injection_rates%deuterium%value( time_sind ) + &
                    & tflux(istrai)*zn(eb2atcr(iatm))
                else
                  allocate( summary%gas_injection_rates%deuterium%value( time_sind ))
                  summary%gas_injection_rates%deuterium%value( time_sind ) = &
                    & tflux(istrai)*zn(eb2atcr(iatm))
                end if
                allocate ( summary%gas_injection_rates%deuterium%source(1) )
                summary%gas_injection_rates%deuterium%source = source
               case ('T')
                if ( associated( summary%gas_injection_rates%tritium%value ) ) then
                  summary%gas_injection_rates%tritium%value( time_sind ) = &
                    & summary%gas_injection_rates%tritium%value( time_sind ) + &
                    & tflux(istrai)*zn(eb2atcr(iatm))
                else
                  allocate( summary%gas_injection_rates%tritium%value( time_sind ))
                  summary%gas_injection_rates%tritium%value( time_sind ) = &
                    & tflux(istrai)*zn(eb2atcr(iatm))
                end if
                allocate ( summary%gas_injection_rates%tritium%source(1) )
                summary%gas_injection_rates%tritium%source = source
               case ('He')
                if (nint(am(is)).eq.3) then
                  if ( associated( summary%gas_injection_rates%helium_3%value ) ) then
                    summary%gas_injection_rates%helium_3%value( time_sind ) = &
                      & summary%gas_injection_rates%helium_3%value( time_sind ) + &
                      & tflux(istrai)*zn(eb2atcr(iatm))
                  else
                    allocate( summary%gas_injection_rates%helium_3%value( time_sind ))
                    summary%gas_injection_rates%helium_3%value( time_sind ) = &
                    & tflux(istrai)*zn(eb2atcr(iatm))
                  end if
                  allocate ( summary%gas_injection_rates%helium_3%source(1) )
                  summary%gas_injection_rates%helium_3%source = source
                else if (nint(am(is)).eq.4) then
                  if ( associated( summary%gas_injection_rates%helium_4%value ) ) then
                    summary%gas_injection_rates%helium_4%value( time_sind ) = &
                      & summary%gas_injection_rates%helium_4%value( time_sind ) + &
                      & tflux(istrai)*zn(eb2atcr(iatm))
                  else
                    allocate( summary%gas_injection_rates%helium_4%value( time_sind ))
                    summary%gas_injection_rates%helium_4%value( time_sind ) = &
                    & tflux(istrai)*zn(eb2atcr(iatm))
                  end if
                  allocate ( summary%gas_injection_rates%helium_4%source(1) )
                  summary%gas_injection_rates%helium_4%source = source
                end if
               case ('Li')
                if ( associated( summary%gas_injection_rates%lithium%value ) ) then
                  summary%gas_injection_rates%lithium%value( time_sind ) = &
                    & summary%gas_injection_rates%lithium%value( time_sind ) + &
                    & tflux(istrai)*zn(eb2atcr(iatm))
                else
                  allocate( summary%gas_injection_rates%lithium%value( time_sind ))
                  summary%gas_injection_rates%lithium%value( time_sind ) = &
                    & tflux(istrai)*zn(eb2atcr(iatm))
                end if
                summary%gas_injection_rates%impurity_seeding%value = 1
                allocate ( summary%gas_injection_rates%lithium%source(1) )
                summary%gas_injection_rates%lithium%source = source
               case ('Be')
                if ( associated( summary%gas_injection_rates%beryllium%value ) ) then
                  summary%gas_injection_rates%beryllium%value( time_sind ) = &
                    & summary%gas_injection_rates%beryllium%value( time_sind ) + &
                    & tflux(istrai)*zn(eb2atcr(iatm))
                else
                  allocate( summary%gas_injection_rates%beryllium%value( time_sind ))
                  summary%gas_injection_rates%beryllium%value( time_sind ) = &
                    & tflux(istrai)*zn(eb2atcr(iatm))
                end if
                summary%gas_injection_rates%impurity_seeding%value = 1
                allocate ( summary%gas_injection_rates%beryllium%source(1) )
                summary%gas_injection_rates%beryllium%source = source
               case ('C')
                if ( associated( summary%gas_injection_rates%carbon%value ) ) then
                  summary%gas_injection_rates%carbon%value( time_sind ) = &
                    & summary%gas_injection_rates%carbon%value( time_sind ) + &
                    & tflux(istrai)*zn(eb2atcr(iatm))
                else
                  allocate( summary%gas_injection_rates%carbon%value( time_sind ))
                  summary%gas_injection_rates%carbon%value( time_sind ) = &
                    & tflux(istrai)*zn(eb2atcr(iatm))
                end if
                summary%gas_injection_rates%impurity_seeding%value = 1
                allocate ( summary%gas_injection_rates%carbon%source(1) )
                summary%gas_injection_rates%carbon%source = source
               case ('N')
                if ( associated( summary%gas_injection_rates%nitrogen%value ) ) then
                  summary%gas_injection_rates%nitrogen%value( time_sind ) = &
                    & summary%gas_injection_rates%nitrogen%value( time_sind ) + &
                    & tflux(istrai)*zn(eb2atcr(iatm))
                else
                  allocate( summary%gas_injection_rates%nitrogen%value( time_sind ))
                  summary%gas_injection_rates%nitrogen%value( time_sind ) = &
                    & tflux(istrai)*zn(eb2atcr(iatm))
                end if
                summary%gas_injection_rates%impurity_seeding%value = 1
                allocate ( summary%gas_injection_rates%nitrogen%source(1) )
                summary%gas_injection_rates%nitrogen%source = source
               case ('O')
                if ( associated( summary%gas_injection_rates%oxygen%value ) ) then
                  summary%gas_injection_rates%oxygen%value( time_sind ) = &
                    & summary%gas_injection_rates%oxygen%value( time_sind ) + &
                    & tflux(istrai)*zn(eb2atcr(iatm))
                else
                  allocate( summary%gas_injection_rates%oxygen%value( time_sind ))
                  summary%gas_injection_rates%oxygen%value( time_sind ) = &
                    & tflux(istrai)*zn(eb2atcr(iatm))
                end if
                summary%gas_injection_rates%impurity_seeding%value = 1
                allocate ( summary%gas_injection_rates%oxygen%source(1) )
                summary%gas_injection_rates%oxygen%source = source
               case ('Ne')
                if ( associated( summary%gas_injection_rates%neon%value ) ) then
                  summary%gas_injection_rates%neon%value( time_sind ) = &
                    & summary%gas_injection_rates%neon%value( time_sind ) + &
                    & tflux(istrai)*zn(eb2atcr(iatm))
                else
                  allocate( summary%gas_injection_rates%neon%value( time_sind ))
                  summary%gas_injection_rates%neon%value( time_sind ) = &
                    & tflux(istrai)*zn(eb2atcr(iatm))
                end if
                summary%gas_injection_rates%impurity_seeding%value = 1
                allocate ( summary%gas_injection_rates%neon%source(1) )
                summary%gas_injection_rates%neon%source = source
               case ('Ar')
                if ( associated( summary%gas_injection_rates%argon%value ) ) then
                  summary%gas_injection_rates%argon%value( time_sind ) = &
                    & summary%gas_injection_rates%argon%value( time_sind ) + &
                    & tflux(istrai)*zn(eb2atcr(iatm))
                else
                  allocate( summary%gas_injection_rates%argon%value( time_sind ))
                  summary%gas_injection_rates%argon%value( time_sind ) = &
                    & tflux(istrai)*zn(eb2atcr(iatm))
                end if
                summary%gas_injection_rates%impurity_seeding%value = 1
                allocate ( summary%gas_injection_rates%argon%source(1) )
                summary%gas_injection_rates%argon%source = source
               case ('Xe')
                if ( associated( summary%gas_injection_rates%xenon%value ) ) then
                  summary%gas_injection_rates%xenon%value( time_sind ) = &
                    & summary%gas_injection_rates%xenon%value( time_sind ) + &
                    & tflux(istrai)*zn(eb2atcr(iatm))
                else
                  allocate( summary%gas_injection_rates%xenon%value( time_sind ))
                  summary%gas_injection_rates%xenon%value( time_sind ) = &
                    & tflux(istrai)*zn(eb2atcr(iatm))
                end if
                summary%gas_injection_rates%impurity_seeding%value = 1
                allocate ( summary%gas_injection_rates%xenon%source(1) )
                summary%gas_injection_rates%xenon%source = source
               case ('Kr')
                if ( associated( summary%gas_injection_rates%krypton%value ) ) then
                  summary%gas_injection_rates%krypton%value( time_sind ) = &
                    & summary%gas_injection_rates%krypton%value( time_sind ) + &
                    & tflux(istrai)*zn(eb2atcr(iatm))
                else
                  allocate( summary%gas_injection_rates%krypton%value( time_sind ))
                  summary%gas_injection_rates%krypton%value( time_sind ) = &
                    & tflux(istrai)*zn(eb2atcr(iatm))
                end if
                summary%gas_injection_rates%impurity_seeding%value = 1
                allocate ( summary%gas_injection_rates%krypton%source(1) )
                summary%gas_injection_rates%krypton%source = source
               end select
              end if
            end do
            if (maxval(gpfc(1:natmi,istrai)).ne.1.0_R8) then
              do iatm = 1, natmi
                if (gpfc(iatm,istrai).gt.0.0_R8) then
                  if (nint(zn(eb2atcr(iatm))).eq.1) iatm1 = iatm
                  if (nint(zn(eb2atcr(iatm))).gt.1) iatm2 = iatm
                end if
              end do
              if (gpfc(iatm1,istrai)+gpfc(iatm2,istrai).le.0.9999_R8) then
                continue ! case not coded
              else
                if (nint(gpfc(iatm1,istrai)/gpfc(iatm2,istrai)).eq.4) then
                  if (nint(am(eb2atcr(iatm1))).eq.1 .and. &
                    & nint(am(eb2atcr(iatm2))).eq.12) then ! CH4
                     if ( associated( summary%gas_injection_rates%methane%value ) ) then
                       summary%gas_injection_rates%methane%value( time_sind ) = &
                        & summary%gas_injection_rates%methane%value( time_sind ) + &
                        & tflux(istrai)*10
                     else
                       allocate( summary%gas_injection_rates%methane%value( time_sind ))
                       summary%gas_injection_rates%methane%value( time_sind ) = &
                        & tflux(istrai)*10
                     end if
                     summary%gas_injection_rates%impurity_seeding%value = 1
                     allocate ( summary%gas_injection_rates%methane%source(1) )
                     summary%gas_injection_rates%methane%source = source
                  else if (nint(am(eb2atcr(iatm1))).eq.1 .and. &
                    &      nint(am(eb2atcr(iatm2))).eq.13) then ! 13CH4
                     if ( associated( summary%gas_injection_rates%methane_carbon_13%value ) ) then
                       summary%gas_injection_rates%methane_carbon_13%value( time_sind ) = &
                        & summary%gas_injection_rates%methane_carbon_13%value( time_sind ) + &
                        & tflux(istrai)*10
                     else
                       allocate( summary%gas_injection_rates%methane_carbon_13%value( time_sind ))
                       summary%gas_injection_rates%methane_carbon_13%value( time_sind ) = &
                        & tflux(istrai)*10
                     end if
                     summary%gas_injection_rates%impurity_seeding%value = 1
                     allocate ( summary%gas_injection_rates%methane_carbon_13%source(1) )
                     summary%gas_injection_rates%methane_carbon_13%source = source
                  else if (nint(am(eb2atcr(iatm1))).eq.2 .and. &
                    &      nint(am(eb2atcr(iatm2))).eq.12) then ! CD4
                     if ( associated( summary%gas_injection_rates%methane_deuterated%value ) ) then
                       summary%gas_injection_rates%methane_deuterated%value( time_sind ) = &
                        & summary%gas_injection_rates%methane_deuterated%value( time_sind ) + &
                        & tflux(istrai)*10
                     else
                       allocate( summary%gas_injection_rates%methane_deuterated%value( time_sind ))
                       summary%gas_injection_rates%methane_deuterated%value( time_sind ) = &
                        & tflux(istrai)*10
                     end if
                     summary%gas_injection_rates%impurity_seeding%value = 1
                     allocate ( summary%gas_injection_rates%methane_deuterated%source(1) )
                     summary%gas_injection_rates%methane_deuterated%source = source
                  else if (nint(am(eb2atcr(iatm1))).eq.1 .and. &
                    &      nint(am(eb2atcr(iatm2))).eq.28) then ! SiH4
                     if ( associated( summary%gas_injection_rates%silane%value ) ) then
                       summary%gas_injection_rates%silane%value( time_sind ) = &
                        & summary%gas_injection_rates%silane%value( time_sind ) + &
                        & tflux(istrai)*18
                     else
                       allocate( summary%gas_injection_rates%silane%value( time_sind ))
                       summary%gas_injection_rates%silane%value( time_sind ) = &
                        & tflux(istrai)*18
                     end if
                     summary%gas_injection_rates%impurity_seeding%value = 1
                     allocate ( summary%gas_injection_rates%silane%source(1) )
                     summary%gas_injection_rates%silane%source = source
                  end if
                else if (nint(gpfc(iatm1,istrai)/gpfc(iatm2,istrai)).eq.2) then
                  if (nint(am(eb2atcr(iatm1))).eq.1 .and. &
                    & nint(am(eb2atcr(iatm2))).eq.12) then ! C2H4
                     if ( associated( summary%gas_injection_rates%ethylene%value ) ) then
                       summary%gas_injection_rates%ethylene%value( time_sind ) = &
                        & summary%gas_injection_rates%ethylene%value( time_sind ) + &
                        & tflux(istrai)*16
                     else
                       allocate( summary%gas_injection_rates%ethylene%value( time_sind ))
                       summary%gas_injection_rates%ethylene%value( time_sind ) = &
                        & tflux(istrai)*16
                     end if
                     summary%gas_injection_rates%impurity_seeding%value = 1
                     allocate ( summary%gas_injection_rates%ethylene%source(1) )
                     summary%gas_injection_rates%ethylene%source = source
                  endif
                else if (nint(gpfc(iatm1,istrai)/gpfc(iatm2,istrai)).eq.3) then ! C2H6, C3H8, NH3, ND3
                  if (nint(am(eb2atcr(iatm1))).eq.1 .and. &
                    & nint(am(eb2atcr(iatm2))).eq.12) then
                    if (nint(gpfc(iatm2,istrai)*6).eq.2) then ! C2H6
                     if ( associated( summary%gas_injection_rates%ethane%value ) ) then
                       summary%gas_injection_rates%ethane%value( time_sind ) = &
                        & summary%gas_injection_rates%ethane%value( time_sind ) + &
                        & tflux(istrai)*18
                     else
                       allocate( summary%gas_injection_rates%ethane%value( time_sind ))
                       summary%gas_injection_rates%ethane%value( time_sind ) = &
                        & tflux(istrai)*18
                     end if
                     summary%gas_injection_rates%impurity_seeding%value = 1
                     allocate ( summary%gas_injection_rates%ethane%source(1) )
                     summary%gas_injection_rates%ethane%source = source
                    else if (nint(gpfc(iatm2,istrai)*8).eq.3) then ! C3H8
                     if ( associated( summary%gas_injection_rates%propane%value ) ) then
                       summary%gas_injection_rates%propane%value( time_sind ) = &
                        & summary%gas_injection_rates%propane%value( time_sind ) + &
                        & tflux(istrai)*26
                     else
                       allocate( summary%gas_injection_rates%propane%value( time_sind ))
                       summary%gas_injection_rates%propane%value( time_sind ) = &
                        & tflux(istrai)*26
                     end if
                     summary%gas_injection_rates%impurity_seeding%value = 1
                     allocate ( summary%gas_injection_rates%propane%source(1) )
                     summary%gas_injection_rates%propane%source = source
                    end if
                  else if (nint(am(eb2atcr(iatm1))).eq.1 .and. &
                    &      nint(am(eb2atcr(iatm2))).eq.14) then ! NH3
                     if ( associated( summary%gas_injection_rates%ammonia%value ) ) then
                       summary%gas_injection_rates%ammonia%value( time_sind ) = &
                        & summary%gas_injection_rates%ammonia%value( time_sind ) + &
                        & tflux(istrai)*10
                     else
                       allocate( summary%gas_injection_rates%ammonia%value( time_sind ))
                       summary%gas_injection_rates%ammonia%value( time_sind ) = &
                        & tflux(istrai)*10
                     end if
                     summary%gas_injection_rates%impurity_seeding%value = 1
                     allocate ( summary%gas_injection_rates%ammonia%source(1) )
                     summary%gas_injection_rates%ammonia%source = source
                  else if (nint(am(eb2atcr(iatm1))).eq.2 .and. &
                    &      nint(am(eb2atcr(iatm2))).eq.14) then ! ND3
                     if ( associated( summary%gas_injection_rates%ammonia_deuterated%value) ) then
                       summary%gas_injection_rates%ammonia_deuterated%value( time_sind ) = &
                        & summary%gas_injection_rates%ammonia_deuterated%value( time_sind ) + &
                        & tflux(istrai)*10
                     else
                       allocate( summary%gas_injection_rates%ammonia_deuterated%value( time_sind ))
                       summary%gas_injection_rates%ammonia_deuterated%value( time_sind ) = &
                        & tflux(istrai)*10
                     end if
                     summary%gas_injection_rates%impurity_seeding%value = 1
                     allocate ( summary%gas_injection_rates%ammonia_deuterated%source(1) )
                     summary%gas_injection_rates%ammonia_deuterated%source = source
                  end if
                end if
              endif
            end if
          end if
        end do
        gtop = gsum - gmid - gbot
        allocate( summary%gas_injection_rates%total%value( time_sind ) )
        summary%gas_injection_rates%total%value( time_sind ) = gsum
        allocate( summary%gas_injection_rates%total%source(1) )
        summary%gas_injection_rates%total%source = source
        allocate( summary%gas_injection_rates%midplane%value( time_sind ) )
        summary%gas_injection_rates%midplane%value( time_sind ) = gmid
        allocate( summary%gas_injection_rates%midplane%source(1) )
        summary%gas_injection_rates%midplane%source = source
        allocate( summary%gas_injection_rates%top%value( time_sind ) )
        summary%gas_injection_rates%top%value( time_sind ) = gtop
        allocate( summary%gas_injection_rates%top%source(1) )
        summary%gas_injection_rates%top%source = source
        allocate( summary%gas_injection_rates%bottom%value( time_sind ) )
        summary%gas_injection_rates%bottom%value( time_sind ) = gbot
        allocate( summary%gas_injection_rates%bottom%source(1) )
        summary%gas_injection_rates%bottom%source = source
#endif
#endif

        allocate( edge_profiles%code%output_flag( time_sind ) )
        edge_profiles%code%output_flag( time_sind ) = 0
        allocate( edge_transport%code%output_flag( time_sind ) )
        edge_transport%code%output_flag( time_sind ) = 0
        allocate( edge_sources%code%output_flag( time_sind ) )
        edge_sources%code%output_flag( time_sind ) = 0
        allocate( radiation%code%output_flag( time_sind ) )
        radiation%code%output_flag( time_sind ) = 0
#if IMAS_MINOR_VERSION > 21
        allocate( summary%code%output_flag( time_sind ) )
        summary%code%output_flag( time_sind ) = 0
#endif

        call logmsg( LOGDEBUG, "b2mod_ual_io.B25_process_ids: done" )

        contains

        !> Write scalar B2 cell quantity to 'ids_generic_grid_scalar'
        !! IMAS IDS data tree node.
        subroutine write_quantity( val, fluxes, value, flux, time_sind )
            use b2mod_interp
            type(ids_generic_grid_scalar), pointer, intent(inout) :: val(:)
                !< Type of IDS data structure, designed for scalar data handling
            ! real(IDS_real), pointer, intent(inout)  :: val(:)
            type (ids_generic_grid_scalar), pointer, intent(inout) :: fluxes(:)
                !< Type of IDS data structure, designed for scalar data handling
                !< (in this case scalars regarding fluxes)
            real(IDS_real), intent(in) :: value( -1:gmap%b2nx, -1:gmap%b2ny )
            real(IDS_real), intent(in) :: flux( -1:gmap%b2nx, -1:gmap%b2ny, 0:1 )
            integer, intent(in) :: time_sind    !< General grid description
                                                !< slice identifier
            real(IDS_real), dimension(:), pointer :: idsdata    !< Array for
                !< handing data field values
            integer :: ival

            !! Allocate data fields for 5+4 grid subsets
            allocate( val(9) )

            !! Write data for Cells grid subset
#if IMAS_MINOR_VERSION < 15
            idsdata => b2_IMAS_Transform_Data_B2_To_IDS( edge_profiles%   &
                &   ggd( time_sind )%grid, GRID_SUBSET_CELLS, gmap, value )
#else
            idsdata => b2_IMAS_Transform_Data_B2_To_IDS( edge_profiles%   &
                &   grid_ggd( time_sind ), GRID_SUBSET_CELLS, gmap, value )
#endif

            ival = 1
            call gridWriteData( val( ival ), GRID_SUBSET_CELLS, idsdata )
            deallocate( idsdata )

            !! Write data for Core Boundary grid subset
            ival = ival + 1
            tmpFace = 0.0_IDS_real
            call value_on_faces( nx, ny, vol, value, tmpFace)
#if IMAS_MINOR_VERSION < 15
            idsdata => b2_IMAS_Transform_Data_B2_To_IDS(                    &
                &   edge_profiles%ggd( time_sind )%grid, iGsCoreBoundary,   &
                &   gmap, tmpFace )
#else
            idsdata => b2_IMAS_Transform_Data_B2_To_IDS(                    &
                &   edge_profiles%grid_ggd( time_sind ), iGsCoreBoundary,   &
                &   gmap, tmpFace )
#endif
            call gridWriteData( val( ival ), iGsCoreBoundary, idsdata )
            deallocate( idsdata )

            !! Write data for Inner Midplane grid subset
            ival = ival + 1
            tmpVx = interpolateToVertices(  &
                &   gmap%b2nx, gmap%b2ny, VX_LOWERLEFT, value )
#if IMAS_MINOR_VERSION < 15
            idsdata => b2_IMAS_Transform_Data_B2_To_IDS_Vertex( &
                &   edge_profiles%ggd( time_sind )%grid,        &
                &   iGsInnerMidplane, gmap, tmpVx )
#else
            idsdata => b2_IMAS_Transform_Data_B2_To_IDS_Vertex( &
                &   edge_profiles%grid_ggd( time_sind ),        &
                &   iGsInnerMidplane, gmap, tmpVx )
#endif
            call gridWriteData( val( ival ), iGsInnerMidplane, idsdata )
            deallocate( idsdata )

            !! Write data for Outer Midplane grid subset
            ival = ival + 1
#if IMAS_MINOR_VERSION < 15
            idsdata => b2_IMAS_Transform_Data_B2_To_IDS_Vertex(             &
                &   edge_profiles%ggd( time_sind )%grid, iGsOuterMidplane,  &
                &   gmap, tmpVx )
#else
            idsdata => b2_IMAS_Transform_Data_B2_To_IDS_Vertex(             &
                &   edge_profiles%grid_ggd( time_sind ), iGsOuterMidplane,  &
                &   gmap, tmpVx )
#endif
            call gridWriteData( val( ival ), iGsOuterMidplane, idsdata )
            deallocate( idsdata )

            !! Write data for Nodes grid subset
            ival = ival + 1
#if IMAS_MINOR_VERSION < 15
            idsdata => b2_IMAS_Transform_Data_B2_To_IDS_Vertex(             &
                &   edge_profiles%ggd( time_sind )%grid, GRID_SUBSET_NODES, &
                &   gmap, tmpVx )
#else
            idsdata => b2_IMAS_Transform_Data_B2_To_IDS_Vertex(             &
                &   edge_profiles%grid_ggd( time_sind ), GRID_SUBSET_NODES, &
                &   gmap, tmpVx )
#endif
            call gridWriteData( val( ival ), GRID_SUBSET_NODES, idsdata )
            deallocate( idsdata )

            !! Write data on fluxes for y-aligned faces, x-aligned faces
            !! (those two are set as default) and Core boundary grid subset
            allocate( fluxes(2) )
            call write_face_vector( fluxes(1), flux, time_sind)
            call write_face_vector( fluxes(2), flux, time_sind, &
                &   gridSubsetId = iGsCoreBoundary )

            !! Write data for Core grid subset
            !! Note: Most of the grid subset is written OK, but at the
            !!       boundary with the seperatrix incorrect values seems
            !!       to be written (checked with ParaView ReadUALEdge plugin)!!
            !!      This is strange as the same procedure
            !!       works as it should for all other grid subsets (e.g. SOL,
            !!       Inner divertor and Outer divertor)
            ival = ival + 1
#if IMAS_MINOR_VERSION < 15
            idsdata => b2_IMAS_Transform_Data_B2_To_IDS( edge_profiles%   &
                &   ggd( time_sind )%grid, iGsCore, gmap, value )
#else
            idsdata => b2_IMAS_Transform_Data_B2_To_IDS( edge_profiles%   &
                &   grid_ggd( time_sind ), iGsCore, gmap, value )
#endif

            call gridWriteData( val( ival ), iGsCore, idsdata )
            deallocate( idsdata )

            !! Write data for SOL grid subset
            ival = ival + 1
#if IMAS_MINOR_VERSION < 15
            idsdata => b2_IMAS_Transform_Data_B2_To_IDS( edge_profiles%   &
                &   ggd( time_sind )%grid, iGsSOL, gmap, value )
#else
            idsdata => b2_IMAS_Transform_Data_B2_To_IDS( edge_profiles%   &
                &   grid_ggd( time_sind ), iGsSOL, gmap, value )
#endif

            call gridWriteData( val( ival ), iGsSOL, idsdata )
            deallocate( idsdata )

            !! Write data for Inner Divertor grid subset
            ival = ival + 1
#if IMAS_MINOR_VERSION < 15
            idsdata => b2_IMAS_Transform_Data_B2_To_IDS( edge_profiles%   &
                &   ggd( time_sind )%grid, iGsIDivertor, gmap, value )
#else
            idsdata => b2_IMAS_Transform_Data_B2_To_IDS( edge_profiles%   &
                &   grid_ggd( time_sind ), iGsIDivertor, gmap, value )
#endif
            call gridWriteData( val( ival ), iGsIDivertor, idsdata )
            deallocate( idsdata )

            !! Write data for Outer Divertor grid subset
            ival = ival + 1
#if IMAS_MINOR_VERSION < 15
            idsdata => b2_IMAS_Transform_Data_B2_To_IDS( edge_profiles%   &
                &   ggd( time_sind )%grid, iGsODivertor, gmap, value )
#else
            idsdata => b2_IMAS_Transform_Data_B2_To_IDS( edge_profiles%   &
                &   grid_ggd( time_sind ), iGsODivertor, gmap, value )
#endif

            call gridWriteData( val( ival ), iGsODivertor, idsdata )
            deallocate( idsdata )
        end subroutine write_quantity

        !> Write a scalar B2 cell quantity to ids_generic_grid_scalar
        subroutine write_cell_scalar( scalar, b2CellData )
            type(ids_generic_grid_scalar), intent(inout), pointer :: scalar(:)
                !< Type of IDS data structure, designed for scalar data handling
            real(IDS_real), intent(in) :: b2CellData(-1:gmap%b2nx, -1:gmap%b2ny)
            real(IDS_real), dimension(:), pointer :: idsdata    !< Array for
                !< handing data field values
            allocate( scalar(1) )

            !! TODO: add checks whether already allocated
#if IMAS_MINOR_VERSION < 15
            idsdata => b2_IMAS_Transform_Data_B2_To_IDS( edge_profiles% &
                &   ggd( time_sind )%grid, GRID_SUBSET_CELLS,           &
                &   gmap, b2CellData )
#else
            idsdata => b2_IMAS_Transform_Data_B2_To_IDS( edge_profiles% &
                &   grid_ggd( time_sind ), GRID_SUBSET_CELLS,           &
                &   gmap, b2CellData )
#endif
            call gridWriteData( scalar(1), GRID_SUBSET_CELLS, idsdata )
            deallocate(idsdata)
        end subroutine write_cell_scalar

        !> Write a vector component B2 cell quantity to ids_generic_grid_vector
        !! components
        !! @note Currently works only with parallel velocity data field
        !! @note Available IDS vector component data fields (vector IDs):
        !!          - VEC_ALIGN_RADIAL_ID ( "radial" ),
        !!          - "diamagnetic",
        !!          - VEC_ALIGN_PARALLEL_ID ( "parallel" ),
        !!          - VEC_ALIGN_POLOIDAL_ID ( "poloidal" ),
        !!          - VEC_ALIGN_TOROIDAL_ID ( "toroidal" )
        subroutine write_cell_vector_component( vectorComponent, b2CellData,    &
                &   vectorID )
            type(ids_generic_grid_vector_components), intent(inout),    &
                &   pointer :: vectorComponent(:) !< Type of IDS data structure,
                    !> designed for vector data handling
            real(IDS_real), intent(in) :: b2CellData(-1:gmap%b2nx, -1:gmap%b2ny)
            real(IDS_real), dimension(:), pointer :: idsdata    !< Array for
                !< handing data field values
            character(len=*), intent(in) :: vectorID    !< Vector ID (e.g.
                                                        !< VEC_ALIGN_RADIAL_ID)

            !! If required, allocate storage
            if ( .not. associated( vectorComponent ) ) then
                allocate( vectorComponent(1) )
            end if

            !! TODO: add checks whether already allocated
#if IMAS_MINOR_VERSION < 15
            idsdata => b2_IMAS_Transform_Data_B2_To_IDS( edge_profiles% &
                &   ggd( time_sind )%grid, GRID_SUBSET_CELLS,           &
                &   gmap, b2CellData )
#else
            idsdata => b2_IMAS_Transform_Data_B2_To_IDS( edge_profiles% &
                &   grid_ggd( time_sind ), GRID_SUBSET_CELLS,           &
                &   gmap, b2CellData )
#endif

            call B2grid_Write_Data_Vector_Components( vectorComponent(1),   &
                &   GRID_SUBSET_CELLS, vectorID, idsdata )
            deallocate(idsdata)
        end subroutine write_cell_vector_component

        !!$> TODO: add to GGD itself (ids_grid_data)!
        !> Write a scalar data field given as a scalar data representation to a
        !! generic grid vector component IDS data fields.
        !!
        !! @note    The routine will make sure the required storage is
        !!          allocated, and will deallocate and re-allocate fields as
        !!          necessary.
        !! @note Currently works only with parallel velocity data field
        !! @note Available IDS vector component data fields:
        !!          - VEC_ALIGN_RADIAL_ID ( "radial" ),
        !!          - "diamagnetic",
        !!          - VEC_ALIGN_PARALLEL_ID ( "parallel" ),
        !!          - VEC_ALIGN_POLOIDAL_ID ( "poloidal" ),
        !!          - VEC_ALIGN_TOROIDAL_ID ( "toroidal" )
        subroutine B2grid_Write_Data_Vector_Components( idsField_vcomp, &
                &   grid_subset_index, vectorID, data)
            type(ids_generic_grid_vector_components), intent(inout) ::  &
                &   idsField_vcomp
                !< Type of IDS data structure, designed for handling data
                !> regarding vector components (parallel, poloidal etc.)
            integer, intent(in) :: grid_subset_index    !< Base grid subset
                                                        !< index
            character(len=*), intent(in) :: vectorID    !< Vector ID (e.g. )
                                                        !< VEC_ALIGN_RADIAL_ID)
            real(ids_real), intent(in) :: data(:)   !< Data field to be written
                !< to IDS data structure leaf that corresponds to specified
                !< vector component

            !! set grid subset index
            idsField_vcomp%grid_subset_index = grid_subset_index

            select case( vectorID )
            case( VEC_ALIGN_RADIAL_ID )
                !! Writing radial quantity
                !! Make sure the data field is properly allocated
                if ( associated( idsField_vcomp%radial ) ) then
                    if ( .not. all( shape( idsField_vcomp%radial ) ==   &
                        &   shape(data) )) then
                        deallocate( idsField_vcomp%radial )
                    end if
                end if
                !! If required, allocate storage
                if ( .not. associated( idsField_vcomp%radial ) ) then
                    allocate(idsField_vcomp%radial( size(data, 1) ))
                end if
                !! copy radial data field
                idsField_vcomp%radial = data
            case( "diamagnetic" )
                !! Writing diamagnetic quantity
                !! Make sure the data field is properly allocated
                if ( associated( idsField_vcomp%diamagnetic ) ) then
                    if ( .not. all( shape( idsField_vcomp%diamagnetic) ==   &
                        &   shape(data) )) then
                        deallocate( idsField_vcomp%diamagnetic )
                    end if
                end if
                !! If required, allocate storage
                if ( .not. associated( idsField_vcomp%diamagnetic ) ) then
                    allocate( idsField_vcomp%diamagnetic( size(data, 1) ) )
                end if
                !! copy diamagnetic data field
                idsField_vcomp%diamagnetic = data
            case( VEC_ALIGN_PARALLEL_ID )
                !! Writing parallel quantity
                !! Make sure the data field is properly allocated
                if ( associated( idsField_vcomp%parallel ) ) then
                    if ( .not. all( shape( idsField_vcomp%parallel ) ==  &
                        &   shape(data) )) then
                        deallocate( idsField_vcomp%parallel )
                    end if
                end if
                !! If required, allocate storage
                if ( .not. associated( idsField_vcomp%parallel ) ) then
                    allocate(idsField_vcomp%parallel( size(data, 1) ))
                end if
                !! copy parallel data field
                idsField_vcomp%parallel = data
            case( VEC_ALIGN_POLOIDAL_ID )
                !! Writing poloidal quantity
                !! Make sure the data field is properly allocated
                if ( associated( idsField_vcomp%poloidal ) ) then
                    if ( .not. all( shape( idsField_vcomp%poloidal ) == &
                        &   shape(data) )) then
                        deallocate( idsField_vcomp%poloidal )
                    end if
                end if
                !! If required, allocate storage
                if ( .not. associated( idsField_vcomp%poloidal ) ) then
                    allocate( idsField_vcomp%poloidal( size(data, 1) ) )
                end if
                !! copy poloidal data field
                idsField_vcomp%poloidal = data
            case( VEC_ALIGN_TOROIDAL_ID )
                !! Writing toroidal quantity
                !! Make sure the data field is properly allocated
                if ( associated( idsField_vcomp%toroidal ) ) then
                    if ( .not. all( shape( idsField_vcomp%toroidal ) ==  &
                        &   shape(data) )) then
                        deallocate( idsField_vcomp%toroidal )
                    end if
                end if
                !! If required, allocate storage
                if ( .not. associated( idsField_vcomp%toroidal ) ) then
                    allocate(idsField_vcomp%toroidal( size(data, 1) ))
                end if
                !! copy toroidal data field
                idsField_vcomp%toroidal = data
            end select

        end subroutine B2grid_Write_Data_Vector_Components

#if 0
        !> Write a vector B2 cell quantity to a complexgrid_vector
        subroutine write_cell_vector( vector, align, alignid, vecdata )
            type(ids_generic_grid_vector), intent(inout) :: vector
            real(IDS_real), intent(in) :: vecdata(-1:gmap%b2nx, &
                &   -1:gmap%b2ny, 0:2)
            integer, intent(in) :: align(3)
            character(LEN=132), intent(in) :: alignid(3)
            real(IDS_real), dimension(:), pointer :: idsdata

            !! internal
            integer :: dim, i

            dim = size(vecdata, 3)

            !! ITM CPO versus IMAS IDS regarding the ITMs vector%comp,
            !! vector%align and vector%alignid:
            !! - ITM vector%comp:
            !!      Holds data on one of the vector components
            !!      ( parallel, poloidal, toroidal etc.). The %comp(:)
            !!      node can hold data for any of those components.
            !!      However the data inside that node must be properly
            !!      specified in order to provide necessary information
            !!      to which component this data relates to.
            !!      IDS does that differently. IDS has specially designed
            !!      nodes with node names being the same as names of the
            !!      components (for example
            !!      edge_profiles.ggd(:)%e_field(:)%parallel).
            !!      Each of those nodes hold data for its intended
            !!      component.
            !! - ITM vector%alignid:
            !!      Alignment information for vector components.
            !!      Describes vector component ID or label
            !!      ("parallel", "toroidal", etc.). In IDS this is not
            !!      needed as a node itself indicates to what vector
            !!      component the data relates to.
            !! - ITM vector%align:
            !!      Alignment information for vector components.
            !!      Holds vector component label (number tag). In IDS this
            !!      is probably not required as, same as for %alignid, a
            !!      node itself indicates to what vector component
            !!      the data relates to.

            !! Fill in vector component data
            do i = 1, dim
#ifdef GGD_OLD
                idsdata => b2_IMAS_Transform_Data_B2_To_IDS(        &
                    &   edge_profiles%ggd( time_sind )%grid,        &
                    &   GRID_SUBSET_CELLS, gmap, vecdata(:,:,i-1))
#else
                idsdata => b2_IMAS_Transform_Data_B2_To_IDS(        &
                    &   edge_profiles%grid_ggd( time_sind ),        &
                    &   GRID_SUBSET_CELLS, gmap, vecdata(:,:,i-1))
#endif
                call gridWriteData( vector, GRID_SUBSET_CELLS, idsdata )
                deallocate(idsdata)
            end do

        end subroutine write_cell_vector
#endif

        !> Write a vector B2 face quantity to a ids_generic_grid_vector
        !! @note    ITM CPO versus IMAS IDS regarding the ITMs vector%comp,
        !!          vector%align and vector%alignid:
        !!          - ITM vector%comp:
        !!               Holds data on one of the vector components
        !!               ( parallel, poloidal, toroidal etc.). The %comp(:)
        !!               node can hold data for any of those components.
        !!               However the data inside that node must be properly
        !!               specified in order to provide necessary information
        !!               to which component this data relates to.
        !!               IDS does that differently. IDS has specially designed
        !!               nodes with node names being the same as names of the
        !!               components (for example
        !!               edge_profiles.ggd(:)%e_field(:)%parallel).
        !!               Each of those nodes hold data for its intended
        !!               component.
        !!          - ITM vector%alignid:
        !!               Alignment information for vector components.
        !!               Describes vector component ID or label
        !!               ("parallel", "toroidal", etc.). In IDS this is not
        !!               needed as a node itself indicates to what vector
        !!               component the data relates to.
        !!          - ITM vector%align:
        !!               Alignment information for vector components.
        !!               Holds vector component label (number tag). In IDS this
        !!               is probably not required as, same as for %alignid, a
        !!               node itself indicates to what vector component
        !!               the data relates to.
        subroutine write_face_vector( vector, b2FaceData, time_sind,    &
                &   gridSubsetId )
            use ids_grid_data ! IGNORE
            type(ids_generic_grid_scalar), intent(inout) :: vector
                !< Type of IDS data structure, designed for scalar data handling
                !< (in this case 1D vector)
            real(IDS_real), intent(in) ::   &
                &   b2FaceData(-1:gmap%b2nx, -1:gmap%b2ny, 0:1)
            integer, intent(in), optional :: gridSubsetId   !< Base grid subset
                                                            !< index
            real(IDS_real), dimension(:), pointer :: idsdata    !< Dummy array
                !< for holding data field values
            integer, intent(in) :: time_sind    !< General grid description

            if ( .not. present(gridSubsetId) ) then
                !! Fill in vector component data
#if IMAS_MINOR_VERSION < 15
                idsdata => b2_IMAS_Transform_Data_B2_To_IDS(    &
                    &   edge_profiles%ggd( time_sind )%grid,    &
                    &   GRID_SUBSET_Y_ALIGNED_FACES, gmap, b2FaceData)
#else
                idsdata => b2_IMAS_Transform_Data_B2_To_IDS(    &
                    &   edge_profiles%grid_ggd( time_sind ),    &
                    &   GRID_SUBSET_Y_ALIGNED_FACES, gmap, b2FaceData)
#endif
                call gridWriteData( vector, GRID_SUBSET_Y_ALIGNED_FACES, idsdata )
                deallocate(idsdata)
#if IMAS_MINOR_VERSION < 15
                idsdata => b2_IMAS_Transform_Data_B2_To_IDS(    &
                    &   edge_profiles%ggd( time_sind )%grid,    &
                    &   GRID_SUBSET_X_ALIGNED_FACES, gmap, b2FaceData)
#else
                idsdata => b2_IMAS_Transform_Data_B2_To_IDS(    &
                    &   edge_profiles%grid_ggd( time_sind ),    &
                    &   GRID_SUBSET_X_ALIGNED_FACES, gmap, b2FaceData)
#endif
                call gridWriteData( vector, GRID_SUBSET_X_ALIGNED_FACES, idsdata )
                deallocate(idsdata)
            else
#if IMAS_MINOR_VERSION < 15
                idsdata => b2_IMAS_Transform_Data_B2_To_IDS(    &
                    &   edge_profiles%ggd( time_sind )%grid,    &
                    &   gridSubsetId, gmap, b2FaceData)
#else
                idsdata => b2_IMAS_Transform_Data_B2_To_IDS(    &
                    &   edge_profiles%grid_ggd( time_sind ),    &
                    &   gridSubsetId, gmap, b2FaceData)
#endif
                call gridWriteData( vector, gridSubsetId, idsdata )
                deallocate(idsdata)
            end if

        end subroutine write_face_vector
    end subroutine B25_process_ids

    !> From the B2 grid, compute the coordinate unit vectors
    !> (poloidal, radial, toroidal)
    subroutine compute_Coordinate_Unit_Vectors( crx, cry, e1, e2, e3 )
        real(IDS_real), intent(in), dimension(-1:,-1:,0:) :: crx    !< Horizontal
            !< coordinates of the four corners of the (ix, iy) cell
        real(IDS_real), intent(in), dimension(-1:,-1:,0:) :: cry    !< Vertical
            !< coordinates of the four corners of the (ix, iy) cell
        real(IDS_real), intent(out),    &
            &   dimension(-1:ubound(crx,1),-1:ubound(crx,2),3) :: e1
            !< First set of coordinates
        real(IDS_real), intent(out),    &
            &   dimension(-1:ubound(crx,1),-1:ubound(crx,2),3) :: e2
            !< Second set of coordinates
        real(IDS_real), intent(out),    &
            &   dimension(-1:ubound(crx,1),-1:ubound(crx,2),3) :: e3
            !< Third set of coordinates

        !! internal
        integer :: ix, iy, ixn, iyn, nx, ny
        real(IDS_real), dimension(0:1) :: cC, cN
        real(IDS_real) :: dir

        e1 = 0.0
        e2 = 0.0
        e3 = 0.0

        !! poloidal vectors
        nx = ubound(crx,1)
        ny = ubound(crx,2)
        do ix = -1, nx
            do iy = -1, ny

                cC = quadCentroid( &
                    & crx(ix, iy, 0), cry(ix, iy, 0), &
                    & crx(ix, iy, 1), cry(ix, iy, 1), &
                    & crx(ix, iy, 2), cry(ix, iy, 2), &
                    & crx(ix, iy, 3), cry(ix, iy, 3) )

                !! poloidal direction
                !! Try to find right neighbour
                dir = 1.0
                ixn = rightix( ix, iy )
                iyn = rightiy( ix, iy )

                if ( .not. isInDomain( nx, ny, ixn, iyn ) ) then
                    !! If not found, try to find left neighbour
                    !! ...and note to invert vector direction
                    dir = -1.0
                    ixn = leftix( ix, iy )
                    iyn = leftiy( ix, iy )
                end if
                if ( .not. isInDomain( nx, ny, ixn, iyn ) ) then
                    ! stop "compute_Coordinate_Unit_Vectors: "// &
                    ! & "not able to find poloidal neighbour for cell"
                    !! skip cell
                    cycle
                end if

                cN = quadCentroid( &
                    & crx(ixn, iyn, 0), cry(ixn, iyn, 0), &
                    & crx(ixn, iyn, 1), cry(ixn, iyn, 1), &
                    & crx(ixn, iyn, 2), cry(ixn, iyn, 2), &
                    & crx(ixn, iyn, 3), cry(ixn, iyn, 3) )

                !! compute vector from one centroid to the other
                e1(ix,iy,1) = cN(0) - cC(0)   !! R
                e1(ix,iy,2) = 0.0             !! phi
                e1(ix,iy,3) = cN(1) - cC(1)   !! Z

                e1(ix,iy,:) = e1(ix,iy,:) * dir  !! fix direction

                !! radial direction
                !! Try to find top neighbour
                dir = 1.0
                ixn = topix( ix, iy )
                iyn = topiy( ix, iy )

                if ( .not. isInDomain( nx, ny, ixn, iyn ) ) then
                    !! If not found, try to find bottom neighbour
                    !! ...and note to invert vector direction
                    dir = -1.0
                    ixn = bottomix( ix, iy )
                    iyn = bottomiy( ix, iy )
                end if
                if ( .not. isInDomain( nx, ny, ixn, iyn ) ) then
                    ! stop "compute_Coordinate_Unit_Vectors: "// &
                        ! &  "not able to find toroidal neighbour for cell"
                    !! skip cell
                    cycle
                end if

                cN = quadCentroid( &
                    & crx(ixn, iyn, 0), cry(ixn, iyn, 0), &
                    & crx(ixn, iyn, 1), cry(ixn, iyn, 1), &
                    & crx(ixn, iyn, 2), cry(ixn, iyn, 2), &
                    & crx(ixn, iyn, 3), cry(ixn, iyn, 3) )

                !! compute vector from one centroid to the other
                e2(ix,iy,1) = cN(0) - cC(0)   !! R
                e2(ix,iy,2) = 0.0             !! phi
                e2(ix,iy,3) = cN(1) - cC(1)   !! Z

                e2(ix,iy,:) = e2(ix,iy,:) * dir  !! fix direction

                !! toroidal direction
                e3(ix,iy,1) = 0.0   !! R
                e3(ix,iy,2) = 1.0   !! phi
                e3(ix,iy,3) = 0.0   !! Z

                !! make unit vectors
                e1(ix,iy,:) = unitVector(e1(ix,iy,:))
                e2(ix,iy,:) = unitVector(e2(ix,iy,:))
                e3(ix,iy,:) = unitVector(e3(ix,iy,:))

            end do
        end do

    end subroutine compute_Coordinate_Unit_Vectors

    !> Return unit vector along direction of given vector
    function unitVector(v) result(unitV)
        real(IDS_real), intent(in) :: v(:)  !< Vector
        real(IDS_real) :: unitV(size(v))    !< Unit vector

        unitV = v / sqrt( sum( v**2 ) )
    end function unitVector

#else
# ifdef ITM_ENVIRONMENT_LOADED

  logical, parameter, private :: INCLUDE_GHOST_CELLS = .false.

contains

  subroutine write_cpo(edgecpo)
    type (type_edge) :: edgecpo

    !! internal
    type(B2GridMap) :: gmap
    type(type_complexgrid_subgrid) :: sg_cell, sg_face, sg_bnd_core
    integer :: is, ns, nx, ny, i
    logical, parameter :: B2_WRITE_DATA = .true.
    real(ITM_R8), dimension(-1:ubound(crx,1),-1:ubound(crx,2),3,3) :: e
    integer :: iSgCore, iSgInnerMidplane, iSgOuterMidplane

    real(ITM_R8) :: tmpFace(-1:ubound(na, 1), -1:ubound(na, 2), 0:1)
    real(ITM_R8) :: tmpVx(-1:ubound(na, 1), -1:ubound(na, 2))


    !! allocate and init the cpo
    allocate(edgecpo%datainfo%dataprovider(1))
    edgecpo%datainfo%dataprovider="IPP"
    allocate(edgecpo%codeparam%codename(1))
    edgecpo%codeparam%codename(1)="B2.5"
    edgecpo%time= 0.0D0

    ns = size(na, 3)
    nx = ubound(na, 1)
    ny = ubound(na, 2)


!! species block
    allocate(edgecpo%species(ns))
    do is = 0, ns-1
       allocate(edgecpo%species(is+1)%label(1))
       call species(is, edgecpo%species(is+1)%label, .false.)
       edgecpo%species(is+1)%amn = am(is)
       edgecpo%species(is+1)%zn = zn(is)
       edgecpo%species(is+1)%zmin = zamin(is)
       edgecpo%species(is+1)%zmax = zamax(is)
    enddo


    !! set up the B2<->CPO mappings
    call b2ITMCreateMap( nx,ny,crx(-1:nx,-1:ny,: ),cry(-1:nx,-1:ny,:),&
        & cflags,leftix,leftiy,rightix,rightiy, &
        & topix,topiy,bottomix,bottomiy, INCLUDE_GHOST_CELLS, gmap )
    mapInitialized = .true.

    !! write grid & subgrids
    call b2ITMFillGridDescription( gmap, edgecpo%grid, &
        & nx,ny,crx(-1:nx,-1:ny,:),cry(-1:nx,-1:ny,:), &
        & leftix,leftiy,rightix,rightiy, &
        & topix,topiy,bottomix,bottomiy, &
        & nnreg, topcut, region, cflags, INCLUDE_GHOST_CELLS, vol, gs, qc )

    call xertst( geometryId( nnreg, periodic_bc, topcut ) == GEOMETRY_SN,   &
        &   "write_cpo: can only do single null" )

    !! Write plasma state

    if ( B2_WRITE_DATA ) then

        call logmsg( LOGDEBUG, "b2mod_ual_io.write_cpo: writing plasma state" )

        iSgCore = gridFindSubGridByName( edgecpo%grid, "Core boundary" )
        iSgInnerMidplane = gridFindSubGridByName( edgecpo%grid, "Inner midplane" )
        iSgOuterMidplane = gridFindSubGridByName( edgecpo%grid, "Outer midplane" )

        !! ne
        call write_quantity( edgecpo%fluid%ne%value, edgecpo%fluid%ne%flux, ne, fne )
        call write_cell_scalar( edgecpo%fluid%ne%source, &
            &   b2CellData = sne(:,:,0) + sne(:,:,1)*ne )

        !! na
        allocate(edgecpo%fluid%ni(ns))
        do is = 1, ns
            call write_quantity( edgecpo%fluid%ni(is)%value, &
                &   edgecpo%fluid%ni(is)%flux,               &
                &   value = na(:,:,is-1),                    &
                &   flux = fna(:,:,:,is-1) )
            call write_cell_scalar( edgecpo%fluid%ni(is)%source, &
                &   b2CellData = sna(:,:,0,is-1) +               &
                &                sna(:,:,1,is-1)*na(:,:,is-1) )
        end do

!!$    ! ue TODO: must be computed, refactor code from b2news into function
!!$    allocate(edgecpo%fluid%ve%comps(1))
!!$    allocate(edgecpo%fluid%ve%align(1))
!!$    allocate(edgecpo%fluid%ve%alignid(1))
!!$    edgecpo%fluid%ve%align(1) = VEC_ALIGN_PARALLEL
!!$    edgecpo%fluid%ve%alignid(1) = VEC_ALIGN_PARALLEL_ID
!!$    call write_cell_scalar( edgecpo%fluid%ve%comps(1)%value, b2CellData = ue(:,:) )

        !! ua
        allocate(edgecpo%fluid%vi(ns))
        do is = 1, ns
            allocate(edgecpo%fluid%vi(is)%comps(1))
            allocate(edgecpo%fluid%vi(is)%align(1))
            allocate(edgecpo%fluid%vi(is)%alignid(1))
            edgecpo%fluid%vi(is)%align(1) = VEC_ALIGN_PARALLEL
            edgecpo%fluid%vi(is)%alignid(1) = VEC_ALIGN_PARALLEL_ID

            call write_cell_scalar( edgecpo%fluid%vi(is)%comps(1)%value, &
                &   b2CellData = ua(:,:,is-1) )
        end do

        !! te
        call write_quantity( edgecpo%fluid%te%value, &
            &   edgecpo%fluid%te%flux,               &
            &   value = te/qe,                       &
            &   flux = fhe )

        !! ti
        allocate(edgecpo%fluid%ti(1))
        call write_quantity( edgecpo%fluid%ti(1)%value, &
            &   edgecpo%fluid%ti(1)%flux,               &
            &   value = ti/qe,                          &
            &   flux = fhi )

        !! po
        call write_cell_scalar( edgecpo%fluid%po%value, po )

        !! B (magnetic field vector)
        allocate(edgecpo%fluid%te_aniso%comps(4))

        !! Compute unit basis vectors along the field directions
        call compute_Coordinate_Unit_Vectors(crx, cry, e(:,:,:,1), &
            &                                          e(:,:,:,2), e(:,:,:,3))

        !! Write the three unit basis vectors
        do i = 1, 3
            allocate(edgecpo%fluid%te_aniso%comps(i)%flux(1))
            call write_cell_vector( edgecpo%fluid%te_aniso%comps(i)%flux(1), &
                & (/ VEC_ALIGN_DEFAULT, VEC_ALIGN_DEFAULT, VEC_ALIGN_DEFAULT /), &
                & (/ VEC_ALIGN_DEFAULT_ID, VEC_ALIGN_DEFAULT_ID, VEC_ALIGN_DEFAULT_ID /), &
                & e(:,:,:,i) )
        end do

        !! write the magnetic field vector in the b2 coordinate system
        allocate(edgecpo%fluid%te_aniso%comps(4)%flux(1))
        call write_cell_vector( edgecpo%fluid%te_aniso%comps(4)%flux(1), &
            & (/ VEC_ALIGN_POLOIDAL, VEC_ALIGN_RADIAL, VEC_ALIGN_TOROIDAL /), &
            & (/ VEC_ALIGN_POLOIDAL_ID, VEC_ALIGN_RADIAL_ID, VEC_ALIGN_TOROIDAL_ID /), &
            &  bb(:,:,0:2) )

    end if

    call logmsg( LOGDEBUG, "b2mod_ual_io.write_cpo: done" )

  contains

    !> Write a scalar B2 cell quantity to a complexgrid_scalar
    subroutine write_quantity( values, fluxes, value, flux )
      use b2mod_interp
      type(type_complexgrid_scalar), pointer, intent(inout) :: values(:)
      type(type_complexgrid_vector), pointer, intent(inout) :: fluxes(:)
      real(ITM_R8), intent(in) :: value(-1:gmap%b2nx, -1:gmap%b2ny)
      real(ITM_R8), intent(in) :: flux(-1:gmap%b2nx, -1:gmap%b2ny, 0:1)
      real(ITM_R8), dimension(:), pointer :: cpodata

      allocate(values(5))
      cpodata => b2ITMTransformDataB2ToCpo( edgecpo%grid, B2_SUBGRID_CELLS, gmap, value )
      call gridWriteData( values(1), B2_SUBGRID_CELLS, cpodata )
      deallocate(cpodata)
      tmpFace = 0.0_ITM_R8
      call value_on_faces(nx,ny,vol,value,tmpFace)
      cpodata => b2ITMTransformDataB2ToCpo( edgecpo%grid, iSgCore, gmap, tmpFace )
      call gridWriteData( values(2), iSgCore, cpodata )
      deallocate(cpodata)
      tmpVx = interpolateToVertices( gmap%b2nx, gmap%b2ny, VX_LOWERLEFT, value )
      cpodata => b2ITMTransformDataB2ToCpoVertex( edgecpo%grid, iSgInnerMidplane, gmap, tmpVx  )
      call gridWriteData( values(3), iSgInnerMidplane, cpodata )
      deallocate(cpodata)
      cpodata => b2ITMTransformDataB2ToCpoVertex( edgecpo%grid, iSgOuterMidplane, gmap, tmpVx )
      call gridWriteData( values(4), iSgOuterMidplane, cpodata )
      deallocate(cpodata)
      cpodata => b2ITMTransformDataB2ToCpoVertex( edgecpo%grid, B2_SUBGRID_NODES, gmap, tmpVx )
      call gridWriteData( values(5), B2_SUBGRID_NODES, cpodata )
      deallocate(cpodata)
      allocate( fluxes(2) )
      call write_face_vector( fluxes(1), flux )
      call write_face_vector( fluxes(2), flux, subgridInd = iSgCore )
    end subroutine write_quantity


    !> Write a scalar B2 cell quantity to a complexgrid_scalar
    subroutine write_cell_scalar(scalar, b2CellData)
      type(type_complexgrid_scalar), intent(inout), pointer :: scalar(:)
      real(ITM_R8), intent(in) :: b2CellData(-1:gmap%b2nx, -1:gmap%b2ny)
      real(ITM_R8), dimension(:), pointer :: cpodata

      !! TODO: add checks whether already allocated
      allocate(scalar(1))
      cpodata => b2ITMTransformDataB2ToCpo( edgecpo%grid, B2_SUBGRID_CELLS, gmap, b2CellData )
      call gridWriteData( scalar(1), B2_SUBGRID_CELLS, cpodata )
      deallocate(cpodata)
    end subroutine write_cell_scalar


    !> Write a vector B2 cell quantity to a complexgrid_vector
    subroutine write_cell_vector(vector, align, alignid, vecdata)
      type(type_complexgrid_vector), intent(inout) :: vector
      real(ITM_R8), intent(in) :: vecdata(-1:gmap%b2nx, -1:gmap%b2ny, 0:2)
      integer, intent(in) :: align(3)
      character(LEN=132), intent(in) :: alignid(3)
      real(ITM_R8), dimension(:), pointer :: cpodata

      !! internal
      integer :: dim, i

      dim = size(vecdata, 3)

      !! TODO: add checks whether already allocated
      allocate(vector%comp(dim))
      allocate(vector%align(dim))
      allocate(vector%alignid(dim))

      !! Fill in alignment information for vector components
      vector%align = align
      vector%alignid = alignid

      !! Fill in vector component data
      do i = 1, dim
         cpodata => b2ITMTransformDataB2ToCpo(edgecpo%grid, B2_SUBGRID_CELLS, gmap, vecdata(:,:,i-1))
         call gridWriteData( vector%comp(i), B2_SUBGRID_CELLS, cpodata )
         deallocate(cpodata)
      end do

    end subroutine write_cell_vector

    !> Write a vector B2 face quantity to a complexgrid_vector
    subroutine write_face_vector(vector, b2FaceData, subgridInd)
      type(type_complexgrid_vector), intent(inout) :: vector
      real(ITM_R8), intent(in) :: b2FaceData(-1:gmap%b2nx, -1:gmap%b2ny, 0:1)
      integer, intent(in), optional :: subgridInd
      real(ITM_R8), dimension(:), pointer :: cpodata

!!$      if ( .not. present(subgridInd) ) then
!!$          ! TODO: add checks whether already allocated
!!$          allocate(vector%comp(2))
!!$          allocate(vector%align(2))
!!$          allocate(vector%alignid(2))
!!$
!!$          ! Fill in alignment information for vector components
!!$          vector%align(1) = VEC_ALIGN_POLOIDAL
!!$          vector%alignid(1) = VEC_ALIGN_POLOIDAL_ID
!!$
!!$          vector%align(2) = VEC_ALIGN_RADIAL
!!$          vector%alignid(2) = VEC_ALIGN_RADIAL_ID
!!$
!!$          ! Fill in vector component data
!!$          cpodata => b2ITMTransformDataB2ToCpo(edgecpo%grid, B2_SUBGRID_FACES_Y, gmap, b2FaceData)
!!$          call gridWriteData( vector%comp(1), B2_SUBGRID_FACES_Y, cpodata )
!!$          deallocate(cpodata)
!!$          cpodata => b2ITMTransformDataB2ToCpo(edgecpo%grid, B2_SUBGRID_FACES_X, gmap, b2FaceData)
!!$          call gridWriteData( vector%comp(2), B2_SUBGRID_FACES_X, cpodata )
!!$          deallocate(cpodata)
!!$      else
!!$          allocate(vector%comp(1))
!!$          allocate(vector%align(1))
!!$          allocate(vector%alignid(1))
!!$
!!$          vector%align(1) = GRID_UNDEFINED
!!$          vector%alignid(1) = ""
!!$
!!$          cpodata => b2ITMTransformDataB2ToCpo(edgecpo%grid, subgridInd, gmap, b2FaceData)
!!$          call gridWriteData( vector%comp(1), subgridInd, cpodata )
!!$          deallocate(cpodata)
!!$      end if

    end subroutine write_face_vector

  end subroutine write_cpo

  !> From the B2 grid, compute the coordinate unit vectors (poloidal, radial. toroidal)
  subroutine compute_Coordinate_Unit_Vectors(crx, cry, e1, e2, e3)
    real(ITM_R8), intent(in), dimension(-1:,-1:,0:) :: crx, cry
    real(ITM_R8), intent(out), dimension(-1:ubound(crx,1),-1:ubound(crx,2),3) :: e1, e2, e3

    !! internal
    integer :: ix, iy, ixn, iyn, nx, ny
    real(ITM_R8), dimension(0:1) :: cC, cN
    real(ITM_R8) :: dir

    e1 = 0.0
    e2 = 0.0
    e3 = 0.0

    !! poloidal vectors

    nx = ubound(crx,1)
    ny = ubound(crx,2)
    do ix = -1, nx
        do iy = -1, ny

            cC = quadCentroid( &
                & crx(ix, iy, 0), cry(ix, iy, 0), &
                & crx(ix, iy, 1), cry(ix, iy, 1), &
                & crx(ix, iy, 2), cry(ix, iy, 2), &
                & crx(ix, iy, 3), cry(ix, iy, 3) )

            !! poloidal direction
            !! Try to find right neighbour
            dir = 1.0
            ixn = rightix( ix, iy )
            iyn = rightiy( ix, iy )

            if ( .not. isInDomain( nx, ny, ixn, iyn ) ) then
                !! If not found, try to find left neighbour
                !! ...and note to invert vector direction
                dir = -1.0
                ixn = leftix( ix, iy )
                iyn = leftiy( ix, iy )
            end if
            if ( .not. isInDomain( nx, ny, ixn, iyn ) ) then
                !stop "compute_Coordinate_Unit_Vectors: not able to find poloidal neighbour for cell"
                !! skip cell
                cycle
            end if

            cN = quadCentroid( &
                & crx(ixn, iyn, 0), cry(ixn, iyn, 0), &
                & crx(ixn, iyn, 1), cry(ixn, iyn, 1), &
                & crx(ixn, iyn, 2), cry(ixn, iyn, 2), &
                & crx(ixn, iyn, 3), cry(ixn, iyn, 3) )

            !! compute vector from one centroid to the other
            e1(ix,iy,1) = cN(0) - cC(0)   !! R
            e1(ix,iy,2) = 0.0             !! phi
            e1(ix,iy,3) = cN(1) - cC(1)   !! Z

            e1(ix,iy,:) = e1(ix,iy,:) * dir  !! fix direction


            !! radial direction
            !! Try to find top neighbour
            dir = 1.0
            ixn = topix( ix, iy )
            iyn = topiy( ix, iy )

            if ( .not. isInDomain( nx, ny, ixn, iyn ) ) then
                !! If not found, try to find bottom neighbour
                !! ...and note to invert vector direction
                dir = -1.0
                ixn = bottomix( ix, iy )
                iyn = bottomiy( ix, iy )
            end if
            if ( .not. isInDomain( nx, ny, ixn, iyn ) ) then
                !stop "compute_Coordinate_Unit_Vectors: not able to find toroidal neighbour for cell"
                !! skip cell
                cycle
            end if

            cN = quadCentroid( &
                & crx(ixn, iyn, 0), cry(ixn, iyn, 0), &
                & crx(ixn, iyn, 1), cry(ixn, iyn, 1), &
                & crx(ixn, iyn, 2), cry(ixn, iyn, 2), &
                & crx(ixn, iyn, 3), cry(ixn, iyn, 3) )

            !! compute vector from one centroid to the other
            e2(ix,iy,1) = cN(0) - cC(0)   !! R
            e2(ix,iy,2) = 0.0             !! phi
            e2(ix,iy,3) = cN(1) - cC(1)   !! Z

            e2(ix,iy,:) = e2(ix,iy,:) * dir  !! fix direction


            !! toroidal direction
            e3(ix,iy,1) = 0.0   !! R
            e3(ix,iy,2) = 1.0   !! phi
            e3(ix,iy,3) = 0.0   !! Z


            !! make unit vectors
            e1(ix,iy,:) = unitVector(e1(ix,iy,:))
            e2(ix,iy,:) = unitVector(e2(ix,iy,:))
            e3(ix,iy,:) = unitVector(e3(ix,iy,:))

        end do
    end do

  end subroutine compute_Coordinate_Unit_Vectors


  !> Return unit vector along direction of given vector
  function unitVector(v) result(unitV)
    real(ITM_R8), intent(in) :: v(:)
    real(ITM_R8) :: unitV(size(v))

    unitV = v / sqrt( sum( v**2 ) )
  end function unitVector

# endif
#endif

end module b2mod_ual_io

!!!Local Variables:
!!! mode: f90
!!! End:
