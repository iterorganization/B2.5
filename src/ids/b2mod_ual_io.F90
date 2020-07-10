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
    use b2mod_external
    use b2mod_interp
    use b2mod_b2cmrc
    use b2mod_b2cmfs
    use b2mod_version
    use b2mod_grid_mapping
    use b2mod_balance &
     & , only : eirene_mc_pael_sne_bal, eirene_mc_pmel_sne_bal, &
     &          eirene_mc_papl_sna_bal, eirene_mc_pmpl_sna_bal, &
     &          eirene_mc_paat_sna_bal, eirene_mc_pmat_sna_bal, &
     &          eirene_mc_piat_sna_bal, eirene_mc_paml_sna_bal, &
     &          eirene_mc_pmml_sna_bal, eirene_mc_piml_sna_bal, &
     &          eirene_mc_paio_sna_bal, eirene_mc_pmio_sna_bal, &
     &          eirene_mc_piio_sna_bal, &
     &          eirene_mc_mapl_smo_bal, eirene_mc_mmpl_smo_bal, &
     &          eirene_mc_eael_she_bal, eirene_mc_emel_she_bal, &
     &          eirene_mc_eapl_shi_bal, eirene_mc_empl_shi_bal, &
     &          b2stel_she_ion_bal, b2stel_she_rec_bal, b2npht_shei_bal, &
     &          b2stel_shi_ion_bal, b2stel_shi_rec_bal, &
     &          read_balance
    use b2mod_b2plot &
     & , only : nxtl, nxtr, jxa, jsep
#ifdef B25_EIRENE
    use eirmod_comusr &
    , only : natmi, nmoli, nioni, nmassa, nchara, nmassm, ncharm, nprt, nchrgi, nchari
#else
    use b2mod_b2plot &
    , only : natmi
#endif

    use logging

    !! UAL Access
    use b2mod_ual_io_grid &
     & , only : INCLUDE_GHOST_CELLS
#ifdef IMAS
#if IMAS_MINOR_VERSION > 11
    !! B2/CPO Mapping
    use b2mod_ual_io_data &
     & , only : b2_IMAS_Transform_Data_B2_To_IDS, &
     &          b2_IMAS_Transform_Data_B2_To_IDS_Vertex
    use b2mod_ual_io_grid &
     & , only : findGridSubsetByName, GridWriteData, &
     &          b2_IMAS_Fill_Grid_Desc
#endif
#if GGD_MINOR_VERSION < 9
    use b2mod_ual_io_grid &
     & , only : GRID_SUBSET_ACTIVE_SEPARATRIX, GRID_SUBSET_BETWEEN_SEPARATRICES, &
     &          GRID_SUBSET_MAIN_WALL, GRID_SUBSET_PFR_WALL, &
     &          GRID_SUBSET_OUTER_DIVERTOR_INACTIVE, GRID_SUBSET_INNER_DIVERTOR_INACTIVE, &
     &          GRID_SUBSET_SECOND_SEPARATRIX, &
     &          GRID_SUBSET_OUTER_BAFFLE_INACTIVE, GRID_SUBSET_INNER_BAFFLE_INACTIVE, &
     &          GRID_SUBSET_OUTER_PFR_WALL_INACTIVE, GRID_SUBSET_INNER_PFR_WALL_INACTIVE, &
     &          GRID_SUBSET_CORE_CUT_INACTIVE, GRID_SUBSET_PFR_CUT_INACTIVE, &
     &          GRID_SUBSET_OUTER_THROAT_INACTIVE, GRID_SUBSET_INNER_THROAT_INACTIVE, &
     &          GRID_SUBSET_OUTER_TARGET_INACTIVE, GRID_SUBSET_INNER_TARGET_INACTIVE
#endif
    use ids_routines       ! IGNORE
    use ids_schemas        ! IGNORE
#if IMAS_MINOR_VERSION > 14
    use ids_utility        ! IGNORE
#endif
#if IMAS_MINOR_VERSION > 11
    use ids_grid_common , &     ! IGNORE
        &   IDS_COORDTYPE_R => COORDTYPE_R,       &
        &   IDS_COORDTYPE_Z => COORDTYPE_Z,       &
        &   IDS_GRID_UNDEFINED => GRID_UNDEFINED
#endif
#else
#ifdef ITM_ENVIRONMENT_LOADED
    use euITM_schemas   ! IGNORE
    use euITM_routines  ! IGNORE
    use itm_grid_common ! IGNORE
#endif
#endif

  implicit none

#ifdef IMAS

#if IMAS_MINOR_VERSION < 9
  integer, parameter :: IDS_REAL = R8
  real(kind=R8), parameter :: IDS_REAL_INVALID = -9.0E40_R8
#endif

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
            &   radiation, description, &
#if IMAS_MINOR_VERSION > 21
            &   summary, &
#endif
#if IMAS_MINOR_VERSION > 25
            &   numerics, run_start_time_IN, run_end_time_IN, &
#endif
            &   time_IN, time_step_IN, shot, run, database, version, &
            &   time_slice_ind_IN, num_time_slices_IN )
#ifdef NO_OPT
!DIR$ NOOPTIMIZE
#endif
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
        type (ids_dataset_description) :: description !< IDS designed to store
            !< a description of the simulation
#if IMAS_MINOR_VERSION > 21
        type (ids_summary) :: summary !< IDS designed to store
            !< run summary data
#endif
#if IMAS_MINOR_VERSION > 25
        type (ids_numerics) :: numerics !< IDS designed to store
            !< run numerics data
        real(IDS_real), intent(in) :: run_start_time_IN, run_end_time_IN !< Run time bounds
#endif
        integer, intent(in) :: shot, run
        character(len=24), intent(in) :: database, version
        real(IDS_real), intent(in), optional :: time_IN !< Time
        real(IDS_real), intent(in), optional :: time_step_IN !< Time step
        integer, intent(in), optional :: time_slice_ind_IN
            !< Time step index for the current time slice
        integer, intent(in), optional :: num_time_slices_IN
            !< Total number of time steps. It is required to beforehand allocate
            !< required ggd(:) array of nodes structure and for additional
            !< checks for correct use of the routine.

        !! Internal variables
        character(len=120) :: AM_label  !< Description of A&M data source
        character(len=132) :: ion_label !< Ion species label (e.g. D+1)
        character(len=132) :: mol_label !< Molecule species label (e.g. D2)
        character(len=12) :: ion_charge !< Ion charge (e.g. '1', '2', etc.)
        character(len=24) :: source     !< Code source
        character(len=2)  :: plate_name(4) !< Divertor plate name
        integer :: ion_charge_int !< Ion charge (e.g. 1, 2, etc.)
        integer :: ion_label_tlen !< Length of the (trimmed) ion label
        integer :: ns    !< Total number of B2.5 species
        integer :: nsion !< Total number of IDS ion species
        integer :: nneut !< Total number of IDS neutral species
        integer :: nx    !< Specifies the number of interior cells
                         !< along the first coordinate
        integer :: ny    !< Specifies the number of interior cells
                         !< along the second coordinate
        integer :: n_process !< Number of radiation processes handled
        integer :: nelems    !< Number of elements present in a molecule or molecular ion
        integer :: is, js, ks !< Species indices (iterators)
        integer :: iss    !< State index
        integer :: ii, jj !< Iterators
        integer :: i      !< Iterator
        integer :: j      !< Iterator
        integer :: k      !< Iterator
        integer :: ix     !< Iterator
        integer :: iy     !< Iterator
        integer :: iatm   !< Atom iterator
        integer :: iatm1  !< Hydrogenic atom index in molecule composition
        integer :: iatm2  !< Non-hydrogenic atom index in molecule composition
        integer :: istrai !< Stratum iterator
        integer :: ntimes !< Number of previous timesteps in IDS
        integer :: is1    !< First ion of an isonuclear sequence
        integer :: is2    !< Last ion of an isonuclear sequence
        integer :: icnt   !< Boundary cell counter
        integer :: ib     !< Boundary condition index
        integer :: ntrgts !< Number of divertor targets
        integer :: o      !< Dummy integer
        integer :: p      !< Dummy integer
        integer, allocatable :: isstat(:) !< Mapping array
                                          !< from B2-Eirene species to IDS neutral states
        integer, allocatable :: imneut(:) !< Mapping array
                                          !< from Eirene molecules to IDS neutrals
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
        real(IDS_real) :: tmpFace( -1:ubound( na, 1), -1:ubound( na, 2), 0:1)
        real(IDS_real) :: tmpVx( -1:ubound( na, 1), -1:ubound( na, 2) )
        real(IDS_real) :: tmpCv( -1:ubound( na, 1), -1:ubound( na, 2) )
        real(IDS_real) :: pz( -1:ubound( na, 1), -1:ubound( na, 2) )
        real(IDS_real) :: pb( -1:ubound( na, 1), -1:ubound( na, 2) )
        real(IDS_real) :: pe( -1:ubound( na, 1), -1:ubound( na, 2) )
        real(IDS_real) :: ue( -1:ubound( na, 1), -1:ubound( na, 2) )
        real(IDS_real) :: roxa( -1:ubound( na, 1), -1:ubound( na, 2), 0:ubound( na, 3) )
        real(IDS_real) :: zeff( -1:ubound( na, 1), -1:ubound( na, 2) )
        real(IDS_real) :: sna0( -1:ubound( na, 1), -1:ubound( na, 2), 0:1, &
                       &        -1:ubound( na, 3) )
        real(IDS_real) :: smo0( -1:ubound( na, 1), -1:ubound( na, 2), 0:3, &
                       &        -1:ubound( na, 3) )
        real(IDS_real) :: she0( -1:ubound( na, 1), -1:ubound( na, 2), 0:3 )
        real(IDS_real) :: shi0( -1:ubound( na, 1), -1:ubound( na, 2), 0:3 )
        real(IDS_real) :: time  !< Generic time
        real(IDS_real) :: time_step !< Time step
        real(IDS_real) :: time_slice_value   !< Time slice value
        real(IDS_real) :: b0, r0, b0r0, b0r0_ref, nibnd, frac, &
            &             u, qetot, qitot, qemax, qimax, lambda, &
            &             vtor, nisep, nasum, gsum, gmid, gbot, gtop
        type(B2GridMap) :: gmap !< Data structure holding an
            !< intermediate grid description to be transferred into a CPO or IDS
        type(ids_generic_grid_dynamic_grid_subset) :: gs_cell
        type(ids_generic_grid_dynamic_grid_subset) :: gs_face
        type(ids_generic_grid_dynamic_grid_subset) :: gs_bnd_core

        integer, parameter :: nsources = 12
        integer, save :: style = 1
        integer, save :: ismain = 1
        integer, save :: ismain0 = 0
        integer, save :: ue_style = 2
        integer, save :: use_eirene = 0
        integer, save :: ids_from_43 = 0
        integer, save :: target_offset = 1
        integer, save :: nesepm_istra = -1
        integer, save :: balance_netcdf = 0
        real(IDS_real), save :: ndes = 0.0_IDS_real
        real(IDS_real), save :: ndes_sol = 0.0_IDS_real
        real(IDS_real), save :: nesepm_pfr = 0.0_IDS_real
        real(IDS_real), save :: nesepm_sol = 0.0_IDS_real
        real(IDS_real), save :: nepedm_sol = 0.0_IDS_real
        real(IDS_real), save :: volrec_sol = 0.0_IDS_real
        real(IDS_real), save :: private_flux_puff = 0.0_IDS_real
        real(IDS_real), save :: BoRiS = 0.0_IDS_real
        character*8 date
        character*10 ctime
        character*5 zone
        integer tvalues(8)
        character*16 usrnam
        character*8 imas_version, ual_version
        character*32 B25_git_version
        character*32 ADAS_git_version
        character*32 get_B25_hash
        character*32 get_ADAS_hash
        character*256 filename
        logical match_found, streql, exists
#ifdef B25_EIRENE
        logical, allocatable :: in_species(:)
        logical isadigit
#endif
#ifdef USE_PXFGETENV
        integer lenval, ierror
#else
#ifdef NAGFOR
        integer lenval, ierror
#endif
#endif
        external usrnam, streql, get_B25_hash, get_ADAS_hash

        !! ===  SET UP IDS ===
        write(0,*) "Setting data for edge_profiles IDS"
        call ipgetr ('b2news_BoRiS', BoRiS)
        call ipgeti ('b2mndt_style', style)
        call ipgeti ('b2mndr_ismain', ismain)
        call ipgeti ('b2sigp_style', ue_style)
        call ipgeti ('ids_from_43', ids_from_43)
        call ipgeti ('b2mndr_eirene', use_eirene)
        call ipgeti ('b2mwti_target_offset', target_offset)
        call ipgetr ('b2stbc_ndes', ndes)
        call ipgetr ('b2stbc_ndes_sol', ndes_sol)
        call ipgetr ('b2stbc_nesepm_pfr', nesepm_pfr)
        call ipgetr ('b2stbc_nesepm_sol', nesepm_sol)
        call ipgetr ('b2stbc_nepedm_sol', nepedm_sol)
        call ipgetr ('b2stbc_volrec_sol', volrec_sol)
        call ipgetr ('b2stbc_private_flux_puff', private_flux_puff)
        call ipgeti ('balance_netcdf', balance_netcdf)
        if (balance_netcdf.ne.0) then
          filename='balance.nc'
          call find_file(filename,exists)
          if (.not.exists) then
            write(*,*) 'Missing balance.nc file: data skipped!'
            balance_netcdf = 0
          end if
        end if
        match_found = .false.
        is = ismain
        do while (is.ge.0 .and. .not.match_found)
          if (is_neutral(is) .and. zn(is).eq.zn(ismain) .and. am(is).eq.am(ismain)) then
            ismain0 = is
            match_found = .true.
          end if
          is = is - 1
        end do
        if (.not.match_found.and.ismain.ne.1) ismain0 = ismain
        call ipgeti ('b2mwti_ismain0', ismain0)
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
        call get_environment_variable('IMAS_VERSION',status=ierror,length=lenval)
        if (ierror.eq.0) call get_environment_variable('IMAS_VERSION',value=imas_version)
        call get_environment_variable('UAL_VERSION',status=ierror,length=lenval)
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
        if (ids_from_43.eq.0) then
          source = "SOLPS-ITER"
        else
          source = "SOLPS4.3"
        end if
#else
        if (ids_from_43.eq.0) then
          source = "B2.5"
        else
          source = "B2"
        end if
#endif
        ns = size( na, 3 )
        nx = ubound( na, 1 )
        ny = ubound( na, 2 )
        do is = 0, ns - 1
           roxa(:,:,is) = am(is)*mp*na(:,:,is)
        end do
        call b2xzef (nx, ny, ns, rz2, na, ne, zeff)
        call b2xppz (nx, ny, ns, ne, na, te, ti, pz)

#ifdef B25_EIRENE
        call eirene_mc(nx, ny, ns, nxtl, nxtr, ismain, time_step, BoRiS, &
          &  sna0, smo0, she0, shi0, .false.)
#endif
        if (balance_netcdf.ne.0) call read_balance

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
        edge_profiles%ids_properties%comment(1) = label
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
        edge_transport%ids_properties%comment(1) = label
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
        edge_sources%ids_properties%comment(1) = label
        !! 2. Allocate edge_sources.time and set it to desired values
        allocate( edge_sources%time(num_time_slices) )
        do i = 1, num_time_slices
          edge_sources%time(i) = time - (num_time_slices-i) * time_step
        end do

        !! Preparing dataset_description IDS for writing
        !! In order to write to IDS database there are next steps that are
        !! mandatory to do, otherwise there is high change that writing to IDS
        !! database will fail
        !! 1. Set homogeneous_time to 0 or 1
        description%ids_properties%homogeneous_time = homogeneous_time
        allocate( description%ids_properties%comment(1) )
        description%ids_properties%comment(1) = label
        !! 2. Allocate description.time and set it to desired values
        allocate( description%time(num_time_slices) )
        do i = 1, num_time_slices
          description%time(i) = time - (num_time_slices-i) * time_step
        end do

#if IMAS_MINOR_VERSION > 21
        !! Preparing summary IDS for writing
        !! In order to write to IDS database there are next steps that are
        !! mandatory to do, otherwise there is high change that writing to IDS
        !! database will fail
        !! 1. Set homogeneous_time to 0 or 1
        summary%ids_properties%homogeneous_time = homogeneous_time
        allocate( summary%ids_properties%comment(1) )
        summary%ids_properties%comment(1) = label
        !! 2. Allocate summary.time and set it to desired values
        allocate( summary%time(num_time_slices) )
        do i = 1, num_time_slices
          summary%time(i) = time - (num_time_slices-i) * time_step
        end do
#endif

#if IMAS_MINOR_VERSION > 25
        !! Preparing numerics IDS for writing
        !! In order to write to IDS database there are next steps that are
        !! mandatory to do, otherwise there is high change that writing to IDS
        !! database will fail
        !! 1. Set homogeneous_time to 0 or 1
        numerics%ids_properties%homogeneous_time = homogeneous_time
        allocate( numerics%ids_properties%comment(1) )
        numerics%ids_properties%comment(1) = label
        !! 2. Allocate numerics.time and set it to desired values
        allocate( numerics%time(num_time_slices) )
        do i = 1, num_time_slices
          numerics%time(i) = time - (num_time_slices-i) * time_step
        end do
        allocate( numerics%time_start(num_time_slices) )
        numerics%time_start(num_time_slices) = run_start_time_IN
        allocate( numerics%time_step(num_time_slices) )
        numerics%time_step(num_time_slices) = time_step
        allocate( numerics%time_end(num_time_slices) )
        numerics%time_end(num_time_slices) = run_end_time_IN
#endif

        !! Preparing radiation IDS for writing
        !! In order to write to IDS database there are next steps that are
        !! mandatory to do, otherwise there is high change that writing to IDS
        !! database will fail
        !! 1. Set homogeneous_time to 0 or 1
        radiation%ids_properties%homogeneous_time = homogeneous_time
        allocate( radiation%ids_properties%comment(1) )
        radiation%ids_properties%comment(1) = label
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
        if ( associated( edge_profiles%ggd ) ) then
           ntimes = size( edge_profiles%ggd )
        else
           ntimes = 0
        end if
        if( ntimes .ne. num_time_slices ) then
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
            allocate( edge_transport%model(1)%identifier%name(1) )
            allocate( edge_transport%model(1)%identifier%description(1) )
            if (ids_from_43.eq.0) then
              if (style.eq.0) then
                edge_transport%model(1)%identifier%name(1) = "SOLPS5.0"
                edge_transport%model(1)%identifier%description(1) = "SOLPS5.0 physics model"
              else if (style.eq.1) then
                edge_transport%model(1)%identifier%name(1) = "SOLPS5.2"
                edge_transport%model(1)%identifier%description(1) = "SOLPS5.2 physics model"
              end if
              edge_transport%model(1)%flux_multiplier = 1.5_IDS_real + BoRiS
            else
              edge_transport%model(1)%identifier%name(1) = "SOLPS4.3"
              edge_transport%model(1)%identifier%description(1) = "SOLPS4.3 physics model"
              edge_transport%model(1)%flux_multiplier = 2.5_IDS_real
            end if
            allocate( edge_transport%model(1)%ggd( num_time_slices ) )
            allocate( edge_sources%source(nsources) )
            do is = 1, nsources
              allocate( edge_sources%source(is)%ggd( num_time_slices ) )
            end do
            !! Total sources
            edge_sources%source(1)%identifier%index = 1
            allocate( edge_sources%source(1)%identifier%name(1) )
            edge_sources%source(1)%identifier%name = "Total"
            allocate( edge_sources%source(1)%identifier%description(1) )
            edge_sources%source(1)%identifier%description = "Total source from "//trim(source)
            !! Background sources
            edge_sources%source(2)%identifier%index = 703
            allocate( edge_sources%source(2)%identifier%name(1) )
            edge_sources%source(2)%identifier%name = "Background"
            allocate( edge_sources%source(2)%identifier%description(1) )
            edge_sources%source(2)%identifier%description = "External sources from "//trim(source)
            !! Prescribed sources
            edge_sources%source(3)%identifier%index = 705
            allocate( edge_sources%source(3)%identifier%name(1) )
            edge_sources%source(3)%identifier%name = "Prescribed"
            allocate( edge_sources%source(3)%identifier%description(1) )
            edge_sources%source(3)%identifier%description = "Boundary conditions sources from "//trim(source)
            !! Time derivatives
            edge_sources%source(4)%identifier%index = 706
            allocate( edge_sources%source(4)%identifier%name(1) )
            edge_sources%source(4)%identifier%name = "Time derivative"
            allocate( edge_sources%source(4)%identifier%description(1) )
            edge_sources%source(4)%identifier%description = "Time derivative sources from "//trim(source)
            !! Atomic ionization
            edge_sources%source(5)%identifier%index = 707
            allocate( edge_sources%source(5)%identifier%name(1) )
            edge_sources%source(5)%identifier%name = "Atomic ionization"
            allocate( edge_sources%source(5)%identifier%description(1) )
            edge_sources%source(5)%identifier%description = "Atomic ionization sources from "//trim(source)
            !! Molecular ionization
            edge_sources%source(6)%identifier%index = 708
            allocate( edge_sources%source(6)%identifier%name(1) )
            edge_sources%source(6)%identifier%name = "Molecular ionization"
            allocate( edge_sources%source(6)%identifier%description(1) )
            edge_sources%source(6)%identifier%description = "Molecular ionization sources from "//trim(source)
            !! Ionization
            edge_sources%source(7)%identifier%index = 709
            allocate( edge_sources%source(7)%identifier%name(1) )
            edge_sources%source(7)%identifier%name = "Ionization"
            allocate( edge_sources%source(7)%identifier%description(1) )
            edge_sources%source(7)%identifier%description = "Ionization sources from "//trim(source)
            !! Recombination
            edge_sources%source(8)%identifier%index = 710
            allocate( edge_sources%source(8)%identifier%name(1) )
            edge_sources%source(8)%identifier%name = "Recombination"
            allocate( edge_sources%source(8)%identifier%description(1) )
            edge_sources%source(8)%identifier%description = "Recombination sources from "//trim(source)
            !! Charge exchange
            edge_sources%source(9)%identifier%index = 305
            allocate( edge_sources%source(9)%identifier%name(1) )
            edge_sources%source(9)%identifier%name = "Charge exchange"
            allocate( edge_sources%source(9)%identifier%description(1) )
            edge_sources%source(9)%identifier%description = "Charge exchange sources from "//trim(source)
            !! Collisional equipartition
            edge_sources%source(10)%identifier%index = 11
            allocate( edge_sources%source(10)%identifier%name(1) )
            edge_sources%source(10)%identifier%name = "Equipartition"
            allocate( edge_sources%source(10)%identifier%description(1) )
            edge_sources%source(10)%identifier%description = "Collisional equipartition sources from "//trim(source)
            !! Ohmic
            edge_sources%source(11)%identifier%index = 7
            allocate( edge_sources%source(11)%identifier%name(1) )
            edge_sources%source(11)%identifier%name = "Ohmic"
            allocate( edge_sources%source(11)%identifier%description(1) )
            edge_sources%source(11)%identifier%description = "Ohmic (Joule) sources from "//trim(source)
            !! Radiation
            edge_sources%source(12)%identifier%index = 200
            allocate( edge_sources%source(12)%identifier%name(1) )
            edge_sources%source(12)%identifier%name = "Radiation"
            allocate( edge_sources%source(12)%identifier%description(1) )
            edge_sources%source(12)%identifier%description = "Radiation sources from "//trim(source)
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

        B25_git_version = get_B25_hash()
        allocate( edge_profiles%code%commit(1) )
        edge_profiles%code%commit = B25_git_version
        allocate( edge_transport%code%commit(1) )
        edge_transport%code%commit = B25_git_version
        allocate( edge_sources%code%commit(1) )
        edge_sources%code%commit = B25_git_version
        allocate( radiation%code%commit(1) )
        if (streql(b2frates_flag,'adas')) then
          ADAS_git_version = get_ADAS_hash()
          radiation%code%commit = 'B25 : '//trim(B25_git_version)// &
                           &  ' + ADAS : '//trim(ADAS_git_version)
        else
          radiation%code%commit = B25_git_version
        endif

        allocate( edge_profiles%code%repository(1) )
        edge_profiles%code%repository = "git.iter.org"
        allocate( edge_transport%code%repository(1) )
        edge_transport%code%repository = "git.iter.org"
        allocate( edge_sources%code%repository(1) )
        edge_sources%code%repository = "git.iter.org"
        allocate( radiation%code%repository(1) )
        radiation%code%repository = "git.iter.org"

        allocate( radiation%ids_properties%source(1) )
        radiation%ids_properties%source = AM_label
        allocate( description%ids_properties%source(1) )
        description%ids_properties%source = source
#if IMAS_MINOR_VERSION > 25
        allocate( numerics%ids_properties%source(1) )
        numerics%ids_properties%source = source
#endif
#if IMAS_MINOR_VERSION > 14
        allocate( edge_profiles%ids_properties%provider(1) )
        edge_profiles%ids_properties%provider = usrnam()
        allocate( edge_transport%ids_properties%provider(1) )
        edge_transport%ids_properties%provider = usrnam()
        allocate( edge_sources%ids_properties%provider(1) )
        edge_sources%ids_properties%provider = usrnam()
        allocate( radiation%ids_properties%provider(1) )
        radiation%ids_properties%provider = usrnam()
        allocate( description%ids_properties%provider(1) )
        description%ids_properties%provider = usrnam()

#if IMAS_MINOR_VERSION > 21
        allocate( summary%code%name(1) )
        summary%code%name = source
        allocate( summary%code%version(1) )
        summary%code%version = newversion
        allocate( summary%code%commit(1) )
        summary%code%commit = get_B25_hash()
        allocate( summary%code%repository(1) )
        summary%code%repository = "git.iter.org"
        allocate( summary%ids_properties%provider(1) )
        summary%ids_properties%provider = usrnam()
#endif

#if IMAS_MINOR_VERSION > 25
        allocate( numerics%ids_properties%provider(1) )
        numerics%ids_properties%provider = usrnam()
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
        allocate( description%ids_properties%creation_date(1) )
        description%ids_properties%creation_date = &
                &   date//' '//ctime//' '//' '//zone
#if IMAS_MINOR_VERSION > 21
        allocate( summary%ids_properties%creation_date(1) )
        summary%ids_properties%creation_date = &
                &   date//' '//ctime//' '//' '//zone
#endif
#if IMAS_MINOR_VERSION > 25
        allocate( numerics%ids_properties%creation_date(1) )
        numerics%ids_properties%creation_date = &
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
        allocate( description%ids_properties%version_put%data_dictionary(1) )
        description%ids_properties%version_put%data_dictionary = imas_version
#if IMAS_MINOR_VERSION > 25
        allocate( numerics%ids_properties%version_put%data_dictionary(1) )
        numerics%ids_properties%version_put%data_dictionary = imas_version
#endif

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
        allocate( description%ids_properties%version_put%access_layer(1) )
        description%ids_properties%version_put%access_layer = ual_version
#if IMAS_MINOR_VERSION > 25
        allocate( numerics%ids_properties%version_put%access_layer(1) )
        numerics%ids_properties%version_put%access_layer = ual_version
#endif

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
        allocate( description%ids_properties%version_put%access_layer_language(1) )
        description%ids_properties%version_put%access_layer_language = 'FORTRAN'
#if IMAS_MINOR_VERSION > 25
        allocate( numerics%ids_properties%version_put%access_layer_language(1) )
        numerics%ids_properties%version_put%access_layer_language = 'FORTRAN'
#endif
#endif

        allocate( description%data_entry%user(1) )
        description%data_entry%user = usrnam()
        allocate( description%data_entry%machine(1) )
        description%data_entry%machine = database
        allocate( description%data_entry%pulse_type(1) )
        description%data_entry%pulse_type = "simulation"
        description%data_entry%pulse = shot
        description%data_entry%run = run
        allocate( description%imas_version(1) )
        description%imas_version = version
        allocate( description%dd_version(1) )
        description%dd_version = imas_version
        description%simulation%time_step = time_step_IN
        description%simulation%time_current = time_IN
        allocate( description%simulation%workflow(1) )
        description%simulation%workflow = source
#if IMAS_MINOR_VERSION > 25
        description%simulation%time_begin = run_start_time_IN
        description%simulation%time_end = run_end_time_IN
#endif

        i=index(B25_git_version,'-')
        allocate( summary%tag%name(1) )
        summary%tag%name = B25_git_version(1:i-1)
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
            if (isymm.eq.0) then
              b0r0 = ffbz(jxa,-1,0)
            else
              b0r0 = ffbz(jxa,-1,0)/(2.0_R8 * pi)
            end if
            b0 = b0r0/r0
          else if (isymm.eq.0) then
            b0 = bb(jxa,-1,2)
            b0r0 = b0*r0
          else if (isymm.eq.1 .or. isymm.eq.2) then
            b0r0 = bb(jxa,-1,2)*(crx(jxa,-1,0)+crx(jxa,-1,1)+ &
                              &  crx(jxa,-1,2)+crx(jxa,-1,3))/4.0
            b0 = b0r0 / r0
          else if (isymm.eq.3 .or. isymm.eq.4) then
            b0r0 = bb(jxa,-1,2)*(cry(jxa,-1,0)+cry(jxa,-1,1)+ &
                              &  cry(jxa,-1,2)+cry(jxa,-1,3))/4.0
            b0 = b0r0 / r0
          end if
        else
          b0 = bb(jxa,-1,2)
          if (isymm.eq.1 .or. isymm.eq.2) then
            b0r0 = bb(jxa,-1,2)*(crx(jxa,-1,0)+crx(jxa,-1,1)+ &
                              &  crx(jxa,-1,2)+crx(jxa,-1,3))/4.0
          else if (isymm.eq.3 .or. isymm.eq.4) then
            b0r0 = bb(jxa,-1,2)*(cry(jxa,-1,0)+cry(jxa,-1,1)+ &
                              &  cry(jxa,-1,2)+cry(jxa,-1,3))/4.0
          end if
        end if
        !> Careful: Sign convention for magnetic field in IDS
        !>          is OPPOSITE to that in SOLPS toroidal geometries
        if ( b0.ne.0.0_IDS_real ) then
          if (streql(database,'iter')) then
            b0r0_ref = 5.3_IDS_real * 6.2_IDS_real
            allocate( summary%global_quantities%ip%value( time_sind ) )
            allocate( edge_profiles%vacuum_toroidal_field%b0( time_sind ) )
            allocate( summary%global_quantities%b0%value( time_sind ) )
            allocate( summary%global_quantities%q_95%value( time_sind ) )
            i = nint(b0r0_ref/b0r0)
            select case (i)
            case (1)
              summary%global_quantities%ip%value( time_sind ) = -15.0e6_IDS_real
              edge_profiles%vacuum_toroidal_field%b0( time_sind ) = -5.3_IDS_real
              summary%global_quantities%b0%value( time_sind ) = -5.3_IDS_real
            case (2)
              summary%global_quantities%ip%value( time_sind ) =  -7.5e6_IDS_real
              edge_profiles%vacuum_toroidal_field%b0( time_sind ) = -2.65_IDS_real
              summary%global_quantities%b0%value( time_sind ) = -2.65_IDS_real
            case (3)
              summary%global_quantities%ip%value( time_sind ) =  -5.0e6_IDS_real
              edge_profiles%vacuum_toroidal_field%b0( time_sind ) = -1.8_IDS_real
              summary%global_quantities%b0%value( time_sind ) = -1.8_IDS_real
            case default
              summary%global_quantities%ip%value( time_sind ) = -15.0e6_IDS_real/nint(b0r0_ref/b0r0)
              edge_profiles%vacuum_toroidal_field%b0( time_sind ) = -b0r0 / 6.2_IDS_real
              summary%global_quantities%b0%value( time_sind ) = -b0r0 / 6.2_IDS_real
            end select
            allocate( summary%global_quantities%ip%source(1) )
            summary%global_quantities%ip%source = "ITER Baseline q95=3 equilibrium"
            summary%global_quantities%q_95%value( time_sind ) = 3.0_IDS_real
            allocate( summary%global_quantities%q_95%source(1) )
            summary%global_quantities%q_95%source = "ITER Baseline q95=3 equilibrium"
            edge_profiles%vacuum_toroidal_field%r0 = 6.2_IDS_real
            summary%global_quantities%r0%value = 6.2_IDS_real
          else
            edge_profiles%vacuum_toroidal_field%r0 = r0
            summary%global_quantities%r0%value = r0
            allocate( edge_profiles%vacuum_toroidal_field%b0( time_sind ) )
            allocate( summary%global_quantities%b0%value( time_sind ) )
            if (isymm.ne.0) then
              edge_profiles%vacuum_toroidal_field%b0( time_sind ) = -b0
              summary%global_quantities%b0%value( time_sind ) = -b0
            else
              edge_profiles%vacuum_toroidal_field%b0( time_sind ) = b0
              summary%global_quantities%b0%value( time_sind ) = b0
            end if
          end if
          allocate( summary%global_quantities%r0%source(1) )
          summary%global_quantities%r0%source = source
          allocate( summary%global_quantities%b0%source(1) )
          summary%global_quantities%b0%source = source
        end if
#endif

        !! Write grid & grid subsets/subgrids
#if IMAS_MINOR_VERSION > 11
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
        do is = 1, nsources
            call b2_IMAS_Fill_Grid_Desc( gmap,                                  &
                &   edge_sources%source(is)%ggd( time_sind )%grid,              &
                &   nx, ny, crx(-1:nx, -1:ny, :), cry(-1:nx, -1:ny, : ),        &
                &   leftix, leftiy, rightix, rightiy, topix, topiy, bottomix,   &
                &   bottomiy, nnreg, topcut, region, cflags,                    &
                &   INCLUDE_GHOST_CELLS, vol, gs, qc )
        end do
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
        do is = 1, nsources
          edge_sources%source(is)%ggd( time_sind )%time = time_slice_value
        end do
#if IMAS_MINOR_VERSION > 21
        do j = 1, n_process
          allocate( radiation%process(j)%ggd( num_time_slices ) )
          radiation%process(j)%ggd( time_sind )%time = time_slice_value
        end do
#endif

        !! List of species
        nsion = fluids_list(ns-1)
#ifdef B25_EIRENE
        if (use_eirene.ne.0) nsion = nsion + nioni
#endif
        allocate( edge_profiles%ggd( time_sind )%ion( nsion ) )
#if IMAS_MINOR_VERSION > 21
        allocate( radiation%process(1)%ggd( time_sind )%ion( nsion ) )
        allocate( radiation%process(2)%ggd( time_sind )%ion( nsion ) )
#endif
        do i = 1, nsources
            allocate( edge_sources%source(i)%ggd( time_sind )%ion( nsion ) )
        end do
        allocate( edge_transport%model(1)%ggd( time_sind )%ion( nsion ) )
        do is = 0, ns-1
            js = fluids_list(is)
            if (js.eq.0) cycle
            allocate( edge_profiles%ggd( time_sind )%ion( js )%label(1) )
            allocate( edge_profiles%ggd( time_sind )%ion( js )%state(1) )
            allocate( edge_profiles%ggd( time_sind )%ion( js )%state(1)%label(1) )
            allocate( edge_profiles%ggd( time_sind )%ion( js )%element(1) )
            do i = 1, nsources
                allocate( edge_sources%source(i)%ggd( time_sind )%ion( js )%label(1) )
                allocate( edge_sources%source(i)%ggd( time_sind )%ion( js )%state(1) )
                allocate( edge_sources%source(i)%ggd( time_sind )%ion( js )%state(1)%label(1) )
                allocate( edge_sources%source(i)%ggd( time_sind )%ion( js )%element(1) )
            end do
            allocate( edge_transport%model(1)%ggd( time_sind )%ion( js )%label(1) )
            allocate( edge_transport%model(1)%ggd( time_sind )%ion( js )%state(1) )
            allocate( edge_transport%model(1)%ggd( time_sind )%ion( js )%state(1)%label(1) )
            allocate( edge_transport%model(1)%ggd( time_sind )%ion( js )%element(1) )

#if IMAS_MINOR_VERSION > 21
            allocate( radiation%process(1)%ggd( time_sind )%ion( js )%label(1) )
            allocate( radiation%process(1)%ggd( time_sind )%ion( js )%state(1) )
            allocate( radiation%process(1)%ggd( time_sind )%ion( js )%state(1)%label(1) )
            allocate( radiation%process(1)%ggd( time_sind )%ion( js )%element(1) )
            allocate( radiation%process(2)%ggd( time_sind )%ion( js )%label(1) )
            allocate( radiation%process(2)%ggd( time_sind )%ion( js )%state(1) )
            allocate( radiation%process(2)%ggd( time_sind )%ion( js )%state(1)%label(1) )
            allocate( radiation%process(2)%ggd( time_sind )%ion( js )%element(1) )
#endif
            ! Put label to ion(is + 1).state(1).label
            call species( is, edge_profiles%ggd( time_sind )%ion( js )% &
                &   state(1)%label, .false.)
            do i = 1, nsources
                call species( is, edge_sources%source(i)%ggd( time_sind )%ion( js )% &
                    &   state(1)%label, .false.)
            end do
            call species( is, edge_transport%model(1)%ggd( time_sind )%ion( js )% &
                &   state(1)%label, .false.)
#if IMAS_MINOR_VERSION > 21
            call species( is, radiation%process(1)%ggd( time_sind )%ion( js )% &
                &   state(1)%label, .false.)
            call species( is, radiation%process(2)%ggd( time_sind )%ion( js )% &
                &   state(1)%label, .false.)
#endif
            ! Set (previous) label
            ion_label = edge_profiles%ggd( time_sind )%ion( js )%state(1)%label(1)
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
            edge_profiles%ggd( time_sind )%ion( js )%label = ion_label
            edge_transport%model(1)%ggd( time_sind )%ion( js )%label = ion_label

            ! Put ion charge
            edge_profiles%ggd( time_sind )%ion( js )%z_ion = ion_charge_int
            edge_transport%model(1)%ggd( time_sind )%ion( js )%z_ion = ion_charge_int

            ! Put mass of ion
            edge_profiles%ggd( time_sind )%ion( js )%element(1)%a = am( is )
            edge_transport%model(1)%ggd( time_sind )%ion( js )%element(1)%a = am( is )

            ! Put nuclear charge
            edge_profiles%ggd( time_sind )%ion( js )%element(1)%z_n = zn( is )
            edge_transport%model(1)%ggd( time_sind )%ion( js )%element(1)%z_n = zn( is )

            ! Put number of atoms
#if IMAS_MINOR_VERSION < 15
            edge_profiles%ggd( time_sind )%ion( js )%element(1)%multiplicity = 1.0_R8
            edge_transport%model(1)%ggd( time_sind )%ion( js )%element(1)%multiplicity = 1.0_R8
#else
            edge_profiles%ggd( time_sind )%ion( js )%element(1)%atoms_n = 1
            edge_transport%model(1)%ggd( time_sind )%ion( js )%element(1)%atoms_n = 1
#endif

            ! Put neutral index
            edge_profiles%ggd( time_sind )%ion( js )%neutral_index = b2eatcr(is)
            edge_transport%model(1)%ggd( time_sind )%ion( js )%neutral_index = b2eatcr(is)

            ! Put multiple states flag
            edge_profiles%ggd( time_sind )%ion( js )%multiple_states_flag = 0
            edge_transport%model(1)%ggd( time_sind )%ion( js )%multiple_states_flag = 0

            ! Put minimum Z of the charge state bundle
            ! (z_min = z_max = 0 for a neutral)
            edge_profiles%ggd( time_sind )%ion( js )%state(1)%z_min = zamin( is )
            edge_transport%model(1)%ggd( time_sind )%ion( js )%state(1)%z_min = zamin( is )

            ! Put maximum Z of the charge state bundle
            edge_profiles%ggd( time_sind )%ion( js )%state(1)%z_max = zamax( is )
            edge_transport%model(1)%ggd( time_sind )%ion( js )%state(1)%z_max = zamax( is )

            do i = 1, nsources
                ! Put (complete) ion label identifying the species
                edge_sources%source(i)%ggd( time_sind )%ion( js )%label = ion_label

                ! Put ion charge
                edge_sources%source(i)%ggd( time_sind )%ion( js )%z_ion = ion_charge_int

                ! Put mass of ion
                edge_sources%source(i)%ggd( time_sind )%ion( js )%element(1)%a = am( is )

                ! Put nuclear charge
                edge_sources%source(i)%ggd( time_sind )%ion( js )%element(1)%z_n = zn( is )

                ! Put number of atoms
#if IMAS_MINOR_VERSION < 15
                edge_sources%source(i)%ggd( time_sind )%ion( js )%element(1)%multiplicity = 1.0_R8
#else
                edge_sources%source(i)%ggd( time_sind )%ion( js )%element(1)%atoms_n = 1
#endif

                ! Put neutral index
                edge_sources%source(i)%ggd( time_sind )%ion( js )%neutral_index = b2eatcr(is)

                ! Put multiple states flag
                edge_sources%source(i)%ggd( time_sind )%ion( js )%multiple_states_flag = 0

                ! Put minimum Z of the charge state bundle
                ! (z_min = z_max = 0 for a neutral)
                edge_sources%source(i)%ggd( time_sind )%ion( js )%state(1)%z_min = zamin( is )

                ! Put maximum Z of the charge state bundle
                edge_sources%source(i)%ggd( time_sind )%ion( js )%state(1)%z_max = zamax( is )
            end do

#if IMAS_MINOR_VERSION > 21
            ! Put (complete) ion label identifying the species
            radiation%process(1)%ggd( time_sind )%ion( js )%label = ion_label
            radiation%process(2)%ggd( time_sind )%ion( js )%label = ion_label

            ! Put ion charge
            radiation%process(1)%ggd( time_sind )%ion( js )%z_ion = ion_charge_int
            radiation%process(2)%ggd( time_sind )%ion( js )%z_ion = ion_charge_int

            ! Put mass of ion
            radiation%process(1)%ggd( time_sind )%ion( js )%element(1)%a = am( is )
            radiation%process(2)%ggd( time_sind )%ion( js )%element(1)%a = am( is )

            ! Put nuclear charge
            radiation%process(1)%ggd( time_sind )%ion( js )%element(1)%z_n = zn( is )
            radiation%process(2)%ggd( time_sind )%ion( js )%element(1)%z_n = zn( is )

            ! Put number of atoms
            radiation%process(1)%ggd( time_sind )%ion( js )%element(1)%atoms_n = 1
            radiation%process(2)%ggd( time_sind )%ion( js )%element(1)%atoms_n = 1

            ! Put neutral index
            radiation%process(1)%ggd( time_sind )%ion( js )%neutral_index = b2eatcr(is)
            radiation%process(2)%ggd( time_sind )%ion( js )%neutral_index = b2eatcr(is)

            ! Put multiple states flag
            radiation%process(1)%ggd( time_sind )%ion( js )%multiple_states_flag = 0
            radiation%process(2)%ggd( time_sind )%ion( js )%multiple_states_flag = 0

            ! Put minimum Z of the charge state bundle
            ! (z_min = z_max = 0 for a neutral)
            radiation%process(1)%ggd( time_sind )%ion( js )%state(1)%z_min = zamin( is )
            radiation%process(2)%ggd( time_sind )%ion( js )%state(1)%z_min = zamin( is )

            ! Put maximum Z of the charge state bundle
            radiation%process(1)%ggd( time_sind )%ion( js )%state(1)%z_max = zamax( is )
            radiation%process(2)%ggd( time_sind )%ion( js )%state(1)%z_max = zamax( is )
#endif
        enddo

        if (use_eirene.ne.0) then
#ifdef B25_EIRENE
            do is = 1, nioni
               js = fluids_list(ns-1) + is
               allocate( edge_profiles%ggd( time_sind )%ion( js )%label(1) )
               allocate( edge_profiles%ggd( time_sind )%ion( js )%state(1) )
               allocate( edge_profiles%ggd( time_sind )%ion( js )%state(1)%label(1) )
               do i = 1, nsources
                   allocate( edge_sources%source(i)%ggd( time_sind )%ion( js )%label(1) )
                   allocate( edge_sources%source(i)%ggd( time_sind )%ion( js )%state(1) )
                   allocate( edge_sources%source(i)%ggd( time_sind )%ion( js )%state(1)%label(1) )
               end do
               allocate( edge_transport%model(1)%ggd( time_sind )%ion( js )%label(1) )
               allocate( edge_transport%model(1)%ggd( time_sind )%ion( js )%state(1) )
               allocate( edge_transport%model(1)%ggd( time_sind )%ion( js )%state(1)%label(1) )
               edge_profiles%ggd( time_sind )%ion( js )%state(1)%label = textin( is-1 )
               edge_transport%model(1)%ggd( time_sind )%ion( js )%state(1)%label = textin( is-1 )
               nelems = count ( micmp( 1:natmi, is ) > 0 )
               allocate( edge_profiles%ggd( time_sind )%ion( js )%element( nelems ) )
               do i = 1, nsources
                  allocate( edge_sources%source(i)%ggd( time_sind )%ion( js )%element( nelems ) )
               end do
               allocate( edge_transport%model(1)%ggd( time_sind )%ion( js )%element( nelems ) )
               i = 0
               do k = 1, natmi
                  if ( micmp( k, is ) > 0 ) then
                     i = i + 1
                     edge_profiles%ggd( time_sind )%ion( js )%element( i )%a = nmassa(k)
                     edge_profiles%ggd( time_sind )%ion( js )%element( i )%z_n = nchara(k)
#if IMAS_MINOR_VERSION < 15
                     edge_profiles%ggd( time_sind )%ion( js )%element( i )%multiplicity = micmp(k,is)
#else
                     edge_profiles%ggd( time_sind )%ion( js )%element( i )%atoms_n = micmp(k,is)
#endif
                     do ii = 1, nsources
                        edge_sources%source(ii)%ggd( time_sind )%ion( js )%element( i )%a = nmassa(k)
                        edge_sources%source(ii)%ggd( time_sind )%ion( js )%element( i )%z_n = nchara(k)
#if IMAS_MINOR_VERSION < 15
                        edge_sources%source(ii)%ggd( time_sind )%ion( js )%element( i )%multiplicity = &
                            &   micmp(k,is)
#else
                        edge_sources%source(ii)%ggd( time_sind )%ion( js )%element( i )%atoms_n = &
                            &   micmp(k,is)
#endif
                     end do
                     edge_transport%model(1)%ggd( time_sind )%ion( js )%element( i )%a = nmassa(k)
                     edge_transport%model(1)%ggd( time_sind )%ion( js )%element( i )%z_n = nchara(k)
#if IMAS_MINOR_VERSION < 15
                     edge_transport%model(1)%ggd( time_sind )%ion( js )%element( i )%multiplicity = micmp(k,is)
#else
                     edge_transport%model(1)%ggd( time_sind )%ion( js )%element( i )%atoms_n = micmp(k,is)
#endif
                  end if
               end do
               edge_profiles%ggd( time_sind )%ion( js )%z_ion = nchrgi( is )
               edge_profiles%ggd( time_sind )%ion( js )%label = textin( is-1 )
               edge_profiles%ggd( time_sind )%ion( js )%neutral_index = lkindi( is )
               edge_profiles%ggd( time_sind )%ion( js )%multiple_states_flag = 0
               edge_profiles%ggd( time_sind )%ion( js )%state(1)%z_min = nchari( is )
               edge_profiles%ggd( time_sind )%ion( js )%state(1)%z_max = nchari( is )
               do i = 1, nsources
                  edge_sources%source(i)%ggd( time_sind )%ion( js )%z_ion = nchrgi( is )
                  edge_sources%source(i)%ggd( time_sind )%ion( js )%label = textin( is-1 )
                  edge_sources%source(i)%ggd( time_sind )%ion( js )%neutral_index = lkindi( is )
                  edge_sources%source(i)%ggd( time_sind )%ion( js )%multiple_states_flag = 0
                  edge_sources%source(i)%ggd( time_sind )%ion( js )%state(1)%z_min = nchari( is )
                  edge_sources%source(i)%ggd( time_sind )%ion( js )%state(1)%z_max = nchari( is )
               end do
               edge_transport%model(1)%ggd( time_sind )%ion( js )%z_ion = nchrgi( is )
               edge_transport%model(1)%ggd( time_sind )%ion( js )%label = textin( is-1 )
               edge_transport%model(1)%ggd( time_sind )%ion( js )%neutral_index = lkindi( is )
               edge_transport%model(1)%ggd( time_sind )%ion( js )%multiple_states_flag = 0
               edge_transport%model(1)%ggd( time_sind )%ion( js )%state(1)%z_min = nchari( is )
               edge_transport%model(1)%ggd( time_sind )%ion( js )%state(1)%z_max = nchari( is )
            end do

            nneut = nspecies
            do j = 1, nmoli
               ks = 1
               do jj = 1, j-1
                  if ( nmassm(j).eq.nmassm(jj) .and. ncharm(j).eq.ncharm(jj) .and. &
                     & nprt(j).eq.nprt(jj) .and. lkindm(j).eq.lkindm(jj) ) then
                     ks = ks + 1
                  end if
               end do
               if (ks.eq.1) nneut = nneut + 1
            end do
            allocate( edge_profiles%ggd( time_sind )%neutral( nneut ) )
            do i = 1, nsources
               allocate( edge_sources%source(i)%ggd( time_sind )%neutral( nneut ) )
            end do
            allocate( edge_transport%model(1)%ggd( time_sind )%neutral( nneut ) )

            !! List of Eirene atoms
            allocate( in_species(nspecies) )
            allocate( isstat(natmi+nmoli) )
            allocate( imneut(nmoli) )
            in_species = .false.
            do is = 1, natmi
               js = latmscl(is)
               if ( in_species(js) ) then
                  isstat(is) = isstat(is-1) + 1
               else
                  in_species(js) = .true.
                  isstat(is) = 1
                  if (js.gt.1) &
                     &  allocate( edge_profiles%ggd( time_sind )%neutral( js-1 )%state( isstat(is-1) ) )
               end if
            end do
            ks = isstat(natmi)
            js = latmscl(natmi)
            allocate( edge_profiles%ggd( time_sind )%neutral( js )%state( ks ) )

            do js = 1, nspecies
               is = eb2spcr(js)
               allocate( edge_profiles%ggd( time_sind )%neutral( js )%element(1) )
               allocate( edge_profiles%ggd( time_sind )%neutral( js )%label(1) )
               edge_profiles%ggd( time_sind )%neutral( js )%element(1)%a = am( is )
               edge_profiles%ggd( time_sind )%neutral( js )%element(1)%z_n = zn( is )
#if IMAS_MINOR_VERSION < 15
               edge_profiles%ggd( time_sind )%neutral( js )%element(1)%multiplicity = 1.0_R8
#else
               edge_profiles%ggd( time_sind )%neutral( js )%element(1)%atoms_n = 1
#endif
               call species( is, edge_profiles%ggd( time_sind )%neutral( js )%label, &
                   &         .false. )
               edge_profiles%ggd( time_sind )%neutral( js )%ion_index = fluids_list( is+1 )
               do i = 1, nsources
                  allocate( edge_sources%source(i)%ggd( time_sind )%neutral( js )%element(1) )
                  allocate( edge_sources%source(i)%ggd( time_sind )%neutral( js )%label(1) )
                  edge_sources%source(i)%ggd( time_sind )%neutral( js )%element(1)%a = am( is )
                  edge_sources%source(i)%ggd( time_sind )%neutral( js )%element(1)%z_n = zn( is )
#if IMAS_MINOR_VERSION < 15
                  edge_sources%source(i)%ggd( time_sind )%neutral( js )%element(1)%multiplicity = 1.0_R8
#else
                  edge_sources%source(i)%ggd( time_sind )%neutral( js )%element(1)%atoms_n = 1
#endif
                  call species( is, edge_sources%source(i)%ggd( time_sind )%neutral( js )%label, &
                     &          .false. )
                  edge_sources%source(i)%ggd( time_sind )%neutral( js )%ion_index = fluids_list( is+1 )
               end do
               allocate( edge_transport%model(1)%ggd( time_sind )%neutral( js )%element(1) )
               allocate( edge_transport%model(1)%ggd( time_sind )%neutral( js )%label(1) )
               edge_transport%model(1)%ggd( time_sind )%neutral( js )%element(1)%a = am( is )
               edge_transport%model(1)%ggd( time_sind )%neutral( js )%element(1)%z_n = zn( is )
#if IMAS_MINOR_VERSION < 15
               edge_transport%model(1)%ggd( time_sind )%neutral( js )%element(1)%multiplicity = 1.0_R8
#else
               edge_transport%model(1)%ggd( time_sind )%neutral( js )%element(1)%atoms_n = 1
#endif
               call species( is, edge_transport%model(1)%ggd( time_sind )%neutral( js )%label, &
                   &         .false. )
               edge_transport%model(1)%ggd( time_sind )%neutral( js )%ion_index = fluids_list( is+1 )
               ks = size(edge_profiles%ggd( time_sind )%neutral( js )%state)
               do i = 1, nsources
                  allocate( edge_sources%source(i)%ggd( time_sind )%neutral( js )%state( ks ) )
               end do
               allocate( edge_transport%model(1)%ggd( time_sind )%neutral( js )%state( ks ) )
               do iss = 1, ks
                  iatm = b2eatcr(is) + iss - 1
                  allocate( edge_profiles%ggd( time_sind )%neutral( js )%state( iss )%label(1) )
                  allocate( edge_profiles%ggd( time_sind )%neutral( js )%state( iss )% &
                     &      neutral_type%name(1) )
                  allocate( edge_profiles%ggd( time_sind )%neutral( js )%state( iss )% &
                     &      neutral_type%description(1) )
                  edge_profiles%ggd( time_sind )%neutral( js )%state( iss )%neutral_type%name = &
                     &     "Kinetic"
                  edge_profiles%ggd( time_sind )%neutral( js )%state( iss )%neutral_type%index = -1
                  edge_profiles%ggd( time_sind )%neutral( js )%state( iss )%neutral_type%description = &
                     &     "Kinetic neutral atoms from Eirene"
                  if (ks.gt.1) then
                     edge_profiles%ggd( time_sind )%neutral( js )%multiple_states_flag = 1
                  else
                     edge_profiles%ggd( time_sind )%neutral( js )%multiple_states_flag = 0
                  end if
                  edge_profiles%ggd( time_sind )%neutral( js )%state( iss )%label = textan( iatm-1 )
                  do i = 1, nsources
                     allocate( edge_sources%source(i)%ggd( time_sind )%neutral( js )%state( iss )%label(1) )
                     allocate( edge_sources%source(i)%ggd( time_sind )%neutral( js )%state( iss )% &
                        &      neutral_type%name(1) )
                     allocate( edge_sources%source(i)%ggd( time_sind )%neutral( js )%state( iss )% &
                        &      neutral_type%description(1) )
                     edge_sources%source(i)%ggd( time_sind )%neutral( js )%state( iss )%neutral_type%name = &
                        &     "Kinetic"
                     edge_sources%source(i)%ggd( time_sind )%neutral( js )%state( iss )%neutral_type%index = -1
                     edge_sources%source(i)%ggd( time_sind )%neutral( js )%state( iss )%neutral_type%description = &
                        &     "Kinetic neutral atoms from Eirene"
                     if (ks.gt.1) then
                        edge_sources%source(i)%ggd( time_sind )%neutral( js )%multiple_states_flag = 1
                     else
                        edge_sources%source(i)%ggd( time_sind )%neutral( js )%multiple_states_flag = 0
                     end if
                     edge_sources%source(i)%ggd( time_sind )%neutral( js )%state( iss )%label = &
                     &   textan( iatm-1 )
                  end do
                  allocate( edge_transport%model(1)%ggd( time_sind )%neutral( js )%state( iss )%label(1) )
                  allocate( edge_transport%model(1)%ggd( time_sind )%neutral( js )%state( iss )% &
                     &      neutral_type%name(1) )
                  allocate( edge_transport%model(1)%ggd( time_sind )%neutral( js )%state( iss )% &
                     &      neutral_type%description(1) )
                  edge_transport%model(1)%ggd( time_sind )%neutral( js )%state( iss )%neutral_type%name = &
                     &     "Kinetic"
                  edge_transport%model(1)%ggd( time_sind )%neutral( js )%state( iss )%neutral_type%index = -1
                  edge_transport%model(1)%ggd( time_sind )%neutral( js )%state( iss )%neutral_type%description = &
                     &     "Kinetic neutral atoms from Eirene"
                  if (ks.gt.1) then
                     edge_transport%model(1)%ggd( time_sind )%neutral( js )%multiple_states_flag = 1
                  else
                     edge_transport%model(1)%ggd( time_sind )%neutral( js )%multiple_states_flag = 0
                  end if
                  edge_transport%model(1)%ggd( time_sind )%neutral( js )%state( iss )%label = textan( iatm-1 )
               end do
            end do

            !! List of molecules
            js = nspecies
            do j = 1, nmoli
               ks = 1
               do jj = 1, j-1
                  if ( nmassm(j).eq.nmassm(jj) .and. ncharm(j).eq.ncharm(jj) .and. &
                     & nprt(j).eq.nprt(jj) .and. lkindm(j).eq.lkindm(jj) ) then
                     ks = ks + 1
                  end if
               end do
               isstat(natmi+j) = ks
               if (ks.eq.1) then
                  if (j.gt.1) then
                     allocate( edge_profiles%ggd( time_sind )%neutral( js )%state( isstat(natmi+j-1) ) )
                     do i = 1, nsources
                        allocate( edge_sources%source(i)%ggd( time_sind )%neutral( js )%state( isstat(natmi+j-1) ) )
                     end do
                     allocate( edge_transport%model(1)%ggd( time_sind )%neutral( js )%state( isstat(natmi+j-1) ) )
                  end if
                  js = js + 1
               end if
               imneut(j) = js
            end do
            allocate( edge_profiles%ggd( time_sind )%neutral( js )%state( ks ) )
            do i = 1, nsources
               allocate( edge_sources%source(i)%ggd( time_sind )%neutral( js )%state( ks ) )
            end do
            allocate( edge_transport%model(1)%ggd( time_sind )%neutral( js )%state( ks ) )

            js = nspecies
            do j = 1, nmoli
               ks = isstat(natmi+j)
               if (ks.eq.1) js = js + 1
               nelems = count ( mlcmp( 1:natmi, j ) > 0 )
               allocate( edge_profiles%ggd( time_sind )%neutral( js )%element( nelems ) )
               allocate( edge_profiles%ggd( time_sind )%neutral( js )%label(1) )
               allocate( edge_profiles%ggd( time_sind )%neutral( js )%state( ks )%label(1) )
               do i = 1, nsources
                  allocate( edge_sources%source(i)%ggd( time_sind )%neutral( js )%element( nelems ) )
                  allocate( edge_sources%source(i)%ggd( time_sind )%neutral( js )%label(1) )
                  allocate( edge_sources%source(i)%ggd( time_sind )%neutral( js )%state( ks )%label(1) )
               end do
               allocate( edge_transport%model(1)%ggd( time_sind )%neutral( js )%element( nelems ) )
               allocate( edge_transport%model(1)%ggd( time_sind )%neutral( js )%label(1) )
               allocate( edge_transport%model(1)%ggd( time_sind )%neutral( js )%state( ks )%label(1) )
               i = 0
               do k = 1, natmi
                  if (mlcmp( k, j ) > 0 ) then
                      i = i + 1
                      edge_profiles%ggd( time_sind )%neutral( js )%element(i)%a = nmassa( k )
                      edge_profiles%ggd( time_sind )%neutral( js )%element(i)%z_n = &
                         &   nchara( k )
#if IMAS_MINOR_VERSION < 15
                      edge_profiles%ggd( time_sind )%neutral( js )%element(i)%multiplicity = &
                         &   mlcmp(k,j)
#else
                      edge_profiles%ggd( time_sind )%neutral( js )%element(i)%atoms_n = &
                         &   mlcmp(k,j)
#endif
                      do icnt = 1, nsources
                         edge_sources%source(icnt)%ggd( time_sind )%neutral( js )%element(i)%a = &
                            &   nmassa( k )
                         edge_sources%source(icnt)%ggd( time_sind )%neutral( js )%element(i)%z_n = &
                            &   nchara( k )
                         edge_sources%source(icnt)%ggd( time_sind )%neutral( js )%element(i)%atoms_n = &
                            &   mlcmp(k,j)
                      end do
                      edge_transport%model(1)%ggd( time_sind )%neutral( js )%element(i)%a = nmassa( k )
                      edge_transport%model(1)%ggd( time_sind )%neutral( js )%element(i)%z_n = &
                         &   nchara( k )
#if IMAS_MINOR_VERSION < 15
                      edge_transport%model(1)%ggd( time_sind )%neutral( js )%element(i)%multiplicity = &
                         &   mlcmp(k,j)
#else
                      edge_transport%model(1)%ggd( time_sind )%neutral( js )%element(i)%atoms_n = &
                         &   mlcmp(k,j)
#endif
                   end if
               end do
               edge_profiles%ggd( time_sind )%neutral( js )%label = textmn( j-1 )
               edge_profiles%ggd( time_sind )%neutral( js )%ion_index = &
                   &    eb2atcr( lmolscl(j) ) + 1
               edge_profiles%ggd( time_sind )%neutral( js )%multiple_states_flag = 0
               edge_profiles%ggd( time_sind )%neutral( js )%state( ks )%label = &
                   &    textmn( j-1 )
               allocate( edge_profiles%ggd( time_sind )%neutral( js )%state( ks )% &
                   &      neutral_type%name(1) )
               allocate( edge_profiles%ggd( time_sind )%neutral( js )%state( ks )% &
                   &      neutral_type%description(1) )
               edge_profiles%ggd( time_sind )%neutral( js )%state( ks )%neutral_type%name = &
                   &     "Kinetic"
               edge_profiles%ggd( time_sind )%neutral( js )%state( ks )%neutral_type%index = -1
               edge_profiles%ggd( time_sind )%neutral( js )%state( ks )%neutral_type%description = &
                   &     "Kinetic neutral molecules from Eirene"
               do i = 1, nsources
                  edge_sources%source(i)%ggd( time_sind )%neutral( js )%label = textmn( j-1 )
                  edge_sources%source(i)%ggd( time_sind )%neutral( js )%ion_index = &
                      &    eb2atcr( lmolscl(j) ) + 1
                  edge_sources%source(i)%ggd( time_sind )%neutral( js )%multiple_states_flag = 0
                  edge_sources%source(i)%ggd( time_sind )%neutral( js )%state( ks )%label = &
                      &    textmn( j-1 )
                  allocate( edge_sources%source(i)%ggd( time_sind )%neutral( js )%state( ks )% &
                      &      neutral_type%name(1) )
                  allocate( edge_sources%source(i)%ggd( time_sind )%neutral( js )%state( ks )% &
                      &      neutral_type%description(1) )
                  edge_sources%source(i)%ggd( time_sind )%neutral( js )%state( ks )%neutral_type%name = &
                      &     "Kinetic"
                  edge_sources%source(i)%ggd( time_sind )%neutral( js )%state( ks )%neutral_type%index = -1
                  edge_sources%source(i)%ggd( time_sind )%neutral( js )%state( ks )%neutral_type%description = &
                      &     "Kinetic neutral molecules from Eirene"
               end do
               edge_transport%model(1)%ggd( time_sind )%neutral( js )%label = textmn( j-1 )
               edge_transport%model(1)%ggd( time_sind )%neutral( js )%ion_index = &
                   &    eb2atcr( lmolscl(j) ) + 1
               edge_transport%model(1)%ggd( time_sind )%neutral( js )%multiple_states_flag = 0
               edge_transport%model(1)%ggd( time_sind )%neutral( js )%state( ks )%label = &
                   &    textmn( j-1 )
               allocate( edge_transport%model(1)%ggd( time_sind )%neutral( js )%state( ks )% &
                   &     neutral_type%name(1) )
               allocate( edge_transport%model(1)%ggd( time_sind )%neutral( js )%state( ks )% &
                   &     neutral_type%description(1) )
               edge_transport%model(1)%ggd( time_sind )%neutral( js )%state( ks )%neutral_type%name = &
                   &     "Kinetic"
               edge_transport%model(1)%ggd( time_sind )%neutral( js )%state( ks )%neutral_type%index = -1
               edge_transport%model(1)%ggd( time_sind )%neutral( js )%state( ks )%neutral_type%description = &
                   &     "Kinetic neutral molecules from Eirene"
            end do
#endif
        else
            allocate( edge_profiles%ggd( time_sind )%neutral( nspecies ) )
            do i = 1, nsources
               allocate( edge_sources%source(i)%ggd( time_sind )%neutral( nspecies ) )
            end do
            allocate( edge_transport%model(1)%ggd( time_sind )%neutral( nspecies ) )
            do js = 1, nspecies
               is = eb2spcr(js)
               allocate( edge_profiles%ggd( time_sind )%neutral( js )%element(1) )
               allocate( edge_profiles%ggd( time_sind )%neutral( js )%label(1) )
               allocate( edge_profiles%ggd( time_sind )%neutral( js )%state(1) )
               allocate( edge_profiles%ggd( time_sind )%neutral( js )%state(1)%label(1) )
               edge_profiles%ggd( time_sind )%neutral( js )%element(1)%a = am( is )
               edge_profiles%ggd( time_sind )%neutral( js )%element(1)%z_n = zn( is )
#if IMAS_MINOR_VERSION < 15
               edge_profiles%ggd( time_sind )%neutral( js )%element(1)%multiplicity = 1
#else
               edge_profiles%ggd( time_sind )%neutral( js )%element(1)%atoms_n = 1
#endif
               call species( is, edge_profiles%ggd( time_sind )%neutral( js )%label, &
                   &         .false. )
               edge_profiles%ggd( time_sind )%neutral( js )%ion_index = fluids_list( is+1 )
               edge_profiles%ggd( time_sind )%neutral( js )%multiple_states_flag = 0
               call species( is, edge_profiles%ggd( time_sind )%neutral( js )%state(1)%label, &
                   &         .false. )
               allocate( edge_profiles%ggd( time_sind )%neutral( js )%state(1)% &
                   &      neutral_type%name(1) )
               allocate( edge_profiles%ggd( time_sind )%neutral( js )%state(1)% &
                   &      neutral_type%description(1) )
               edge_profiles%ggd( time_sind )%neutral( js )%state(1)%neutral_type%name = &
                   &     "Thermal"
               edge_profiles%ggd( time_sind )%neutral( js )%state(1)%neutral_type%index = 2
               edge_profiles%ggd( time_sind )%neutral( js )%state(1)%neutral_type%description = &
                   &     "Fluid neutral species from B2.5"
               do i = 1, nsources
                  allocate( edge_sources%source(i)%ggd( time_sind )%neutral( js )%element(1) )
                  allocate( edge_sources%source(i)%ggd( time_sind )%neutral( js )%label(1) )
                  allocate( edge_sources%source(i)%ggd( time_sind )%neutral( js )%state(1) )
                  allocate( edge_sources%source(i)%ggd( time_sind )%neutral( js )%state(1)%label(1) )
                  edge_sources%source(i)%ggd( time_sind )%neutral( js )%element(1)%a = am( is )
                  edge_sources%source(i)%ggd( time_sind )%neutral( js )%element(1)%z_n = &
                     &   zn( is )
#if IMAS_MINOR_VERSION < 15
                  edge_sources%source(i)%ggd( time_sind )%neutral( js )%element(1)%multiplicity = 1.0_R8
#else
                  edge_sources%source(i)%ggd( time_sind )%neutral( js )%element(1)%atoms_n = 1
#endif
                  call species( is, edge_sources%source(i)%ggd( time_sind )%neutral( js )%label, &
                     &   .false. )
                  edge_sources%source(i)%ggd( time_sind )%neutral( js )%ion_index = &
                     &   fluids_list( is+1 )
                  edge_sources%source(i)%ggd( time_sind )%neutral( js )%multiple_states_flag = 0
                  call species( is, edge_sources%source(i)%ggd( time_sind )%neutral( js )%state(1)%label, &
                     &   .false. )
                  allocate( edge_sources%source(i)%ggd( time_sind )%neutral( js )%state(1)% &
                     &      neutral_type%name(1) )
                  allocate( edge_sources%source(i)%ggd( time_sind )%neutral( js )%state(1)% &
                     &      neutral_type%description(1) )
                  edge_sources%source(i)%ggd( time_sind )%neutral( js )%state(1)%neutral_type%name = &
                     &     "Thermal"
                  edge_sources%source(i)%ggd( time_sind )%neutral( js )%state(1)%neutral_type%index = 2
                  edge_sources%source(i)%ggd( time_sind )%neutral( js )%state(1)%neutral_type%description = &
                     &     "Fluid neutral species from B2.5"
               end do
               allocate( edge_transport%model(1)%ggd( time_sind )%neutral( js )%element(1) )
               allocate( edge_transport%model(1)%ggd( time_sind )%neutral( js )%label(1) )
               allocate( edge_transport%model(1)%ggd( time_sind )%neutral( js )%state(1) )
               allocate( edge_transport%model(1)%ggd( time_sind )%neutral( js )%state(1)%label(1) )
               edge_transport%model(1)%ggd( time_sind )%neutral( js )%element(1)%a = am( is )
               edge_transport%model(1)%ggd( time_sind )%neutral( js )%element(1)%z_n = zn( is )
#if IMAS_MINOR_VERSION < 15
               edge_transport%model(1)%ggd( time_sind )%neutral( js )%element(1)%multiplicity = 1
#else
               edge_transport%model(1)%ggd( time_sind )%neutral( js )%element(1)%atoms_n = 1
#endif
               call species( is, edge_transport%model(1)%ggd( time_sind )%neutral( js )%label, &
                   &         .false. )
               edge_transport%model(1)%ggd( time_sind )%neutral( js )%ion_index = fluids_list( is+1 )
               edge_transport%model(1)%ggd( time_sind )%neutral( js )%multiple_states_flag = 0
               call species( is, edge_transport%model(1)%ggd( time_sind )%neutral( js )%state(1)%label, &
                   &         .false. )
               allocate( edge_transport%model(1)%ggd( time_sind )%neutral( js )%state(1)% &
                   &      neutral_type%name(1) )
               allocate( edge_transport%model(1)%ggd( time_sind )%neutral( js )%state(1)% &
                   &      neutral_type%description(1) )
               edge_transport%model(1)%ggd( time_sind )%neutral( js )%state(1)%neutral_type%name = &
                   &     "Thermal"
               edge_transport%model(1)%ggd( time_sind )%neutral( js )%state(1)%neutral_type%index = 2
               edge_transport%model(1)%ggd( time_sind )%neutral( js )%state(1)%neutral_type%description = &
                   &     "Fluid neutral species from B2.5"
            end do
#if IMAS_MINOR_VERSION > 21
            allocate( radiation%process(1)%ggd( time_sind )%neutral( nspecies ) )
            do js = 1, nspecies
               is = eb2spcr(js)
               allocate( radiation%process(1)%ggd( time_sind )%neutral( js )%element(1) )
               allocate( radiation%process(1)%ggd( time_sind )%neutral( js )%label(1) )
               allocate( radiation%process(1)%ggd( time_sind )%neutral( js )%state(1) )
               allocate( radiation%process(1)%ggd( time_sind )%neutral( js )%state(1)%label(1) )
               radiation%process(1)%ggd( time_sind )%neutral( js )%element(1)%a = am(is)
               radiation%process(1)%ggd( time_sind )%neutral( js )%element(1)%z_n = zn(is)
               radiation%process(1)%ggd( time_sind )%neutral( js )%element(1)%atoms_n = 1
               call species( is, radiation%process(1)%ggd( time_sind )%neutral( js )%label, &
                   &        .false. )
               radiation%process(1)%ggd( time_sind )%neutral( js )%ion_index = fluids_list( is+1 )
               radiation%process(1)%ggd( time_sind )%neutral( js )%multiple_states_flag = 0
               call species( is, radiation%process(1)%ggd( time_sind )%neutral( js )%state(1)%label, &
                   &        .false. )
               allocate( radiation%process(1)%ggd( time_sind )%neutral( js )%state(1)% &
                   &     neutral_type%name(1) )
               allocate( radiation%process(1)%ggd( time_sind )%neutral( js )%state(1)% &
                   &     neutral_type%description(1) )
               radiation%process(1)%ggd( time_sind )%neutral( js )%state(1)%neutral_type%name = &
                   &     "Thermal"
               radiation%process(1)%ggd( time_sind )%neutral( js )%state(1)%neutral_type%index = 2
               radiation%process(1)%ggd( time_sind )%neutral( js )%state(1)%neutral_type%description = &
                   &     "Fluid neutral species from B2.5"
            end do
#endif
        end if

#if IMAS_MINOR_VERSION > 21
#ifdef B25_EIRENE
        if (use_eirene.ne.0) then
          allocate( radiation%process(3)%ggd( time_sind )%neutral( nneut ) )
          do is = 1, nneut
             ks = size( edge_profiles%ggd( time_sind )%neutral( is )%state )
             allocate( radiation%process(3)%ggd( time_sind )%neutral( is )%state( ks ) )
          end do

          !! List of Eirene atoms
          do is = 1, natmi
            js = latmscl(is)
            ks = isstat(is)
            allocate( radiation%process(3)%ggd( time_sind )%neutral( js )%element(1) )
            allocate( radiation%process(3)%ggd( time_sind )%neutral( js )%label(1) )
            allocate( radiation%process(3)%ggd( time_sind )%neutral( js )%state( ks )%label(1) )
            radiation%process(3)%ggd( time_sind )%neutral( js )%element(1)%a = nmassa( is )
            radiation%process(3)%ggd( time_sind )%neutral( js )%element(1)%z_n = nchara( is )
            radiation%process(3)%ggd( time_sind )%neutral( js )%element(1)%atoms_n = 1
            radiation%process(3)%ggd( time_sind )%neutral( js )%label = &
                &    species_list(js)
            radiation%process(3)%ggd( time_sind )%neutral( js )%ion_index = eb2atcr( is ) + 1
            if ( size(radiation%process(3)%ggd( time_sind )%neutral( js )%state).eq.1 ) then
              radiation%process(3)%ggd( time_sind )%neutral( js )%multiple_states_flag = 0
            else
              radiation%process(3)%ggd( time_sind )%neutral( js )%multiple_states_flag = 1
            end if
            radiation%process(3)%ggd( time_sind )%neutral( js )%state( ks )%label = &
                &    textan( is-1 )
            allocate( radiation%process(3)%ggd( time_sind )%neutral( js )%state( ks )% &
                &     neutral_type%name(1) )
            allocate( radiation%process(3)%ggd( time_sind )%neutral( js )%state( ks )% &
                &     neutral_type%description(1) )
            radiation%process(3)%ggd( time_sind )%neutral( js )%state( ks )%neutral_type%name = &
                &     "Kinetic"
            radiation%process(3)%ggd( time_sind )%neutral( js )%state( ks )%neutral_type%index = -1
            radiation%process(3)%ggd( time_sind )%neutral( js )%state( ks )%neutral_type%description = &
                &     "Kinetic neutral atoms from Eirene"
          end do

          !! List of molecules
          js = nspecies
          do j = 1, nmoli
            ks = isstat(natmi+j)
            if (ks.eq.1) js = js + 1
            nelems = count ( mlcmp( 1:natmi, j ) > 0 )
            allocate( radiation%process(3)%ggd( time_sind )%neutral( js )%element( nelems ) )
            allocate( radiation%process(3)%ggd( time_sind )%neutral( js )%label(1) )
            allocate( radiation%process(3)%ggd( time_sind )%neutral( js )%state( ks )%label(1) )
            i = 0
            do k = 1, natmi
              if (mlcmp( k, j ) > 0 ) then
                i = i + 1
                radiation%process(3)%ggd( time_sind )%neutral( js )%element(i)%a = nmassa( k )
                radiation%process(3)%ggd( time_sind )%neutral( js )%element(i)%z_n = nchara( k )
                radiation%process(3)%ggd( time_sind )%neutral( js )%element(i)%atoms_n = &
                    &   mlcmp(k,j)
              end if
            end do
            radiation%process(3)%ggd( time_sind )%neutral( js )%label = &
                &    textmn( j-1 )
            radiation%process(3)%ggd( time_sind )%neutral( js )%ion_index = &
                &    eb2atcr( lmolscl(j) ) + 1
            radiation%process(3)%ggd( time_sind )%neutral( js )%multiple_states_flag = 0
            radiation%process(3)%ggd( time_sind )%neutral( js )%state( ks )%label = &
                &    textmn( j-1 )
            allocate( radiation%process(3)%ggd( time_sind )%neutral( js )%state( ks )% &
                &     neutral_type%name(1) )
            allocate( radiation%process(3)%ggd( time_sind )%neutral( js )%state( ks )% &
                &     neutral_type%description(1) )
            radiation%process(3)%ggd( time_sind )%neutral( js )%state( ks )%neutral_type%name = &
                &     "Kinetic"
            radiation%process(3)%ggd( time_sind )%neutral( js )%state( ks )%neutral_type%index = -1
            radiation%process(3)%ggd( time_sind )%neutral( js )%state( ks )%neutral_type%description = &
                &     "Kinetic neutral molecules from Eirene"
          end do

          !! List of molecular ions
          allocate( radiation%process(4)%ggd( time_sind )%ion( nsion ) )
          do is = 1, nioni
            js = fluids_list(ns-1) + is
            allocate( radiation%process(4)%ggd( time_sind )%ion( js )%label(1) )
            allocate( radiation%process(4)%ggd( time_sind )%ion( js )%state(1) )
            allocate( radiation%process(4)%ggd( time_sind )%ion( js )%state(1)%label(1) )

            radiation%process(4)%ggd( time_sind )%ion( js )%state(1)%label = textin( is-1 )
            ion_label = adjustl(textin( is-1 ))
            nelems = count ( micmp( 1:natmi, is ) > 0 )
            allocate( radiation%process(4)%ggd( time_sind )%ion( js )%element( nelems ) )
            i = 0
            do k = 1, natmi
              if ( micmp( k, is ) > 0 ) then
                i = i + 1
                radiation%process(4)%ggd( time_sind )%ion( js )%element( i )%a = nmassa(k)
                radiation%process(4)%ggd( time_sind )%ion( js )%element( i )%z_n = nchara(k)
                radiation%process(4)%ggd( time_sind )%ion( js )%element( i )%atoms_n = micmp(k,is)
              end if
            end do
            radiation%process(4)%ggd( time_sind )%ion( js )%z_ion = nchrgi( is )
            radiation%process(4)%ggd( time_sind )%ion( js )%label = textin( is-1 )
            radiation%process(4)%ggd( time_sind )%ion( js )%neutral_index = lkindi( is )
            radiation%process(4)%ggd( time_sind )%ion( js )%multiple_states_flag = 0
            radiation%process(4)%ggd( time_sind )%ion( js )%state(1)%z_min = nchari( is )
            radiation%process(4)%ggd( time_sind )%ion( js )%state(1)%z_max = nchari( is )
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
#if IMAS_MINOR_VERSION > 11
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

            !! ne: Electron density
            call write_quantity(                                            &
                &   val = edge_profiles%ggd( time_sind )%electrons%density, &
                &   value = ne,                                             &
                &   time_sind = time_sind )
            !! fne: Electron particle flux
            do ix = -1, nx
              do iy = -1, ny
                if (fne(ix,iy,0).ne.0.0_R8) then
                  tmpFace(ix,iy,0) = fne(ix,iy,0)/gs(ix,iy,0)/qc(ix,iy)
                else
                  tmpFace(ix,iy,0) = 0.0_R8
                end if
                if (fne(ix,iy,1).ne.0.0_R8) then
                  tmpFace(ix,iy,1) = fne(ix,iy,1)/gs(ix,iy,1)
                else
                  tmpFace(ix,iy,1) = 0.0_R8
                end if
              end do
            end do
            call write_face_scalar(                                         &
                &   val = edge_transport%model(1)%ggd( time_sind )%         &
                &         electrons%particles%flux,                         &
                &   value = tmpFace,                                        &
                &   time_sind = time_sind )
            !! sne: Electron particle sources
            tmpCv(:,:) = ( sne(:,:,0) + sne(:,:,1) * ne(:,:) ) / vol(:,:)
            call write_cell_scalar(                                         &
                &   scalar = edge_sources%source(1)%ggd( time_sind )%       &
                &            electrons%particles,                           &
                &   b2CellData = tmpCv )
            tmpCv(:,:) = ext_sne(:,:) / vol(:,:)
            call write_cell_scalar(                                         &
                &   scalar = edge_sources%source(2)%ggd( time_sind )%       &
                &            electrons%particles,                           &
                &   b2CellData = tmpCv )
            tmpCv(:,:) = ( b2stbc_sne(:,:) + b2stbm_sne(:,:) ) / vol(:,:)
            call write_cell_scalar(                                         &
                &   scalar = edge_sources%source(3)%ggd( time_sind )%       &
                &            electrons%particles,                           &
                &   b2CellData = tmpCv )
            tmpCv(:,:) = ( snedt(:,:,0) + snedt(:,:,1) * ne(:,:) ) / vol(:,:)
            call write_cell_scalar(                                         &
                &   scalar = edge_sources%source(4)%ggd( time_sind )%       &
                &            electrons%particles,                           &
                &   b2CellData = tmpCv )
            if (use_eirene.ne.0 .and. balance_netcdf.ne.0) then
                tmpCv = 0.0_IDS_real
                do istrai = 1, size( eirene_mc_pael_sne_bal, 3)
                    tmpCv(:,:) = tmpCv(:,:) + eirene_mc_pael_sne_bal(:,:,istrai)
                end do
                tmpCv(:,:) = tmpCv(:,:) / vol(:,:)
                call write_cell_scalar(                                     &
                    &   scalar = edge_sources%source(5)%ggd( time_sind )%   &
                    &            electrons%particles,                       &
                    &   b2CellData = tmpCv )
                tmpCv = 0.0_IDS_real
                do istrai = 1, size( eirene_mc_pmel_sne_bal, 3)
                    tmpCv(:,:) = tmpCv(:,:) + eirene_mc_pmel_sne_bal(:,:,istrai)
                end do
                tmpCv(:,:) = tmpCv(:,:) / vol(:,:)
                call write_cell_scalar(                                     &
                    &   scalar = edge_sources%source(6)%ggd( time_sind )%   &
                    &            electrons%particles,                       &
                    &   b2CellData = tmpCv )
            else
                tmpCv(:,:) = b2stbr_sne(:,:) / vol(:,:)
                call write_cell_scalar(                                     &
                    &   scalar = edge_sources%source(5)%ggd( time_sind )%   &
                    &            electrons%particles,                       &
                    &   b2CellData = tmpCv )
            end if
            tmpCv = 0.0_IDS_real
            do is = 0, ns-1
              tmpCv(:,:) = tmpCv(:,:) + rsana(:,:,is)
            end do
            tmpCv(:,:) = tmpCv(:,:) / vol(:,:)
            call write_cell_scalar(                                         &
                &   scalar = edge_sources%source(7)%ggd( time_sind )%       &
                &            electrons%particles,                           &
                &   b2CellData = tmpCv )
            tmpCv = 0.0_IDS_real
            do is = 0, ns-1
              tmpCv(:,:) = tmpCv(:,:) + rrana(:,:,is)
            end do
            tmpCv(:,:) = tmpCv(:,:) / vol(:,:)
            call write_cell_scalar(                                         &
                &   scalar = edge_sources%source(8)%ggd( time_sind )%       &
                &            electrons%particles,                           &
                &   b2CellData = tmpCv )

            !! na: Ion density
            do is = 1, nsion
                if (is.le.fluids_list(ns-1)) then
                    js = ions_list(is)
                    call write_quantity(                                          &
                        &   val = edge_profiles%ggd( time_sind )%ion(is)%density, &
                        &   value = na(:,:,js),                                   &
                        &   time_sind = time_sind )
            !! fna: Ion particle flux
                    do ix = -1, nx
                      do iy = -1, ny
                        if (fna(ix,iy,0,js).ne.0.0_R8) then
                          tmpFace(ix,iy,0) = fna(ix,iy,0,js)/gs(ix,iy,0)/qc(ix,iy)
                        else
                          tmpFace(ix,iy,0) = 0.0_R8
                        end if
                        if (fna(ix,iy,1,js).ne.0.0_R8) then
                          tmpFace(ix,iy,1) = fna(ix,iy,1,js)/gs(ix,iy,1)
                        else
                          tmpFace(ix,iy,1) = 0.0_R8
                        end if
                      end do
                    end do
                    call write_face_scalar(                                   &
                        &   val = edge_transport%model(1)%ggd( time_sind )%   &
                        &         ion( is )%particles%flux,                   &
                        &   value = tmpFace,                                  &
                        &   time_sind = time_sind )
                    call write_face_scalar(                                   &
                        &   val = edge_transport%model(1)%ggd( time_sind )%   &
                        &         ion( is )%state(1)%particles%flux,          &
                        &   value = tmpFace,                                  &
                        &   time_sind = time_sind )
            !! cdna: Ion diffusivity
                    call write_face_scalar(                                   &
                        &   val = edge_transport%model(1)%ggd( time_sind )%   &
                        &         ion( is )%particles%d,                      &
                        &   value = cdna(:,:,:,js),                           &
                        &   time_sind = time_sind )
                    call write_face_scalar(                                   &
                        &   val = edge_transport%model(1)%ggd( time_sind )%   &
                        &         ion( is )%state(1)%particles%d,             &
                        &   value = cdna(:,:,:,js),                           &
                        &   time_sind = time_sind )
            !! cvla: Ion diffusivity
                    call write_face_scalar(                                   &
                        &   val = edge_transport%model(1)%ggd( time_sind )%   &
                        &         ion( is )%particles%v,                      &
                        &   value = cvla(:,:,:,js),                           &
                        &   time_sind = time_sind )
                    call write_face_scalar(                                   &
                        &   val = edge_transport%model(1)%ggd( time_sind )%   &
                        &         ion( is )%state(1)%particles%v,             &
                        &   value = cvla(:,:,:,js),                           &
                        &   time_sind = time_sind )
            !! fllim0fna: Ion flux limiter
                    call write_face_scalar(                                   &
                        &   val = edge_transport%model(1)%ggd( time_sind )%   &
                        &         ion( is )%particles%flux_limiter,           &
                        &   value = fllim0fna(:,:,:,js),                      &
                        &   time_sind = time_sind )
                    call write_face_scalar(                                   &
                        &   val = edge_transport%model(1)%ggd( time_sind )%   &
                        &         ion( is )%state(1)%particles%flux_limiter,  &
                        &   value = fllim0fna(:,:,:,js),                      &
                        &   time_sind = time_sind )
            !! sna: Ion particle sources
                    tmpCv(:,:) = ( sna(:,:,0, is - 1 ) +                      &
                        &          sna(:,:,1, is - 1 ) * na(:,:, is - 1 ) ) / vol(:,:)
                    call write_cell_scalar( scalar = edge_sources%source(1)%  &
                        &   ggd( time_sind )%ion( is )%particles,             &
                        &   b2CellData = tmpCv )
                    call write_cell_scalar( scalar = edge_sources%source(1)%  &
                        &   ggd( time_sind )%ion( is )%state(1)%particles,    &
                        &   b2CellData = tmpCv )
                    tmpCv(:,:) = ext_sna(:,:,js) / vol(:,:)
                    call write_cell_scalar(                                   &
                        &   scalar = edge_sources%source(2)%ggd( time_sind )% &
                        &            ion( is )%particles,                     &
                        &   b2CellData = tmpCv )
                    call write_cell_scalar(                                   &
                        &   scalar = edge_sources%source(2)%ggd( time_sind )% &
                        &            ion( is )%state(1)%particles,            &
                        &   b2CellData = tmpCv )
                    tmpCv(:,:) = ( b2stbc_sna(:,:,js) +                       &
                        &          b2stbm_sna(:,:,js) ) / vol(:,:)
                    call write_cell_scalar(                                   &
                        &   scalar = edge_sources%source(3)%ggd( time_sind )% &
                        &            ion( is )%particles,                     &
                        &   b2CellData = tmpCv )
                    call write_cell_scalar(                                   &
                        &   scalar = edge_sources%source(3)%ggd( time_sind )% &
                        &            ion( is )%state(1)%particles,            &
                        &   b2CellData = tmpCv )
                    tmpCv(:,:) = ( snadt(:,:,0,js) +                          &
                        &          snadt(:,:,1,js) * na(:,:,js) ) / vol(:,:)
                    call write_cell_scalar(                                   &
                        &   scalar = edge_sources%source(4)%ggd( time_sind )% &
                        &            ion( is )%particles,                     &
                        &   b2CellData = tmpCv )
                    call write_cell_scalar(                                   &
                        &   scalar = edge_sources%source(4)%ggd( time_sind )% &
                        &            ion( is )%state(1)%particles,            &
                        &   b2CellData = tmpCv )
                    if (use_eirene.ne.0 .and. balance_netcdf.ne.0) then
                        tmpCv = 0.0_IDS_real
                        do istrai = 1, size( eirene_mc_papl_sna_bal, 4)
                            tmpCv(:,:) = tmpCv(:,:)                           &
                               &       + eirene_mc_papl_sna_bal(:,:,js,istrai)
                        end do
                        tmpCv(:,:) = tmpCv(:,:) / vol(:,:)
                        call write_cell_scalar(                                     &
                            &   scalar = edge_sources%source(5)%ggd( time_sind )%   &
                            &            ion( is )%particles,                       &
                            &   b2CellData = tmpCv )
                        call write_cell_scalar(                                     &
                            &   scalar = edge_sources%source(5)%ggd( time_sind )%   &
                            &            ion( is )%state(1)%particles,              &
                            &   b2CellData = tmpCv )
                        tmpCv = 0.0_IDS_real
                        do istrai = 1, size( eirene_mc_pmpl_sna_bal, 4)
                            tmpCv(:,:) = tmpCv(:,:)                                 &
                               &       + eirene_mc_pmpl_sna_bal(:,:,js,istrai)
                        end do
                        tmpCv(:,:) = tmpCv(:,:) / vol(:,:)
                        call write_cell_scalar(                                     &
                            &   scalar = edge_sources%source(6)%ggd( time_sind )%   &
                            &            ion( is )%particles,                       &
                            &   b2CellData = tmpCv )
                        call write_cell_scalar(                                     &
                            &   scalar = edge_sources%source(6)%ggd( time_sind )%   &
                            &            ion( is )%state(1)%particles,              &
                            &   b2CellData = tmpCv )
                    else
                        tmpCv(:,:) = b2stbr_sna(:,:,js) / vol(:,:)
                        call write_cell_scalar(                                     &
                            &   scalar = edge_sources%source(5)%ggd( time_sind )%   &
                            &            ion( is )%particles,                       &
                            &   b2CellData = tmpCv )
                        call write_cell_scalar(                                     &
                            &   scalar = edge_sources%source(5)%ggd( time_sind )%   &
                            &            ion( is )%state(1)%particles,              &
                            &   b2CellData = tmpCv )
                    end if
                    tmpCv(:,:) = rsana(:,:,js) / vol(:,:)
                    call write_cell_scalar(                                       &
                        &   scalar = edge_sources%source(7)%ggd( time_sind )%     &
                        &            ion( is )%particles,                         &
                        &   b2CellData = tmpCv )
                    call write_cell_scalar(                                       &
                        &   scalar = edge_sources%source(7)%ggd( time_sind )%     &
                        &            ion( is )%state(1)%particles,                &
                        &   b2CellData = tmpCv )
                    tmpCv(:,:) = rrana(:,:,js) / vol(:,:)
                    call write_cell_scalar(                                       &
                        &   scalar = edge_sources%source(8)%ggd( time_sind )%     &
                        &            ion( is )%particles,                         &
                        &   b2CellData = tmpCv )
                    call write_cell_scalar(                                       &
                        &   scalar = edge_sources%source(8)%ggd( time_sind )%     &
                        &            ion( is )%state(1)%particles,                &
                        &   b2CellData = tmpCv )
                    tmpCv(:,:) = rcxna(:,:,js) / vol(:,:)
                    call write_cell_scalar(                                       &
                        &   scalar = edge_sources%source(9)%ggd( time_sind )%     &
                        &            ion( is )%particles,                         &
                        &   b2CellData = tmpCv )
                    call write_cell_scalar(                                       &
                        &   scalar = edge_sources%source(9)%ggd( time_sind )%     &
                        &            ion( is )%state(1)%particles,                &
                        &   b2CellData = tmpCv )
                else
                    js = is - fluids_list(ns-1)
                    tmpCv(-1:nx,-1:ny) = dib2(0:nx+1,0:ny+1,js,1)
                    call write_quantity(                                      &
                        &   val = edge_profiles%ggd( time_sind )%             &
                        &         ion(is)%state(1)%density,                   &
                        &   value = tmpCv,                                    &
                        &   time_sind = time_sind )
                    if (balance_netcdf.ne.0) then
                      tmpCv = 0.0_IDS_real
                      do istrai = 1, size(eirene_mc_paio_sna_bal,4)
                         tmpCv(:,:) = tmpCv(:,:)                              &
                           &       + eirene_mc_paio_sna_bal(:,:,js,istrai)    &
                           &       + eirene_mc_pmio_sna_bal(:,:,js,istrai)    &
                           &       + eirene_mc_piio_sna_bal(:,:,js,istrai)
                      end do
                      tmpCv(:,:) = tmpCv(:,:) / vol(:,:)
                      call write_cell_scalar(                                 &
                        &   scalar = edge_sources%source(1)%ggd( time_sind )% &
                        &            ion( is )%particles,                     &
                        &   b2CellData = tmpCv )
                      call write_cell_scalar(                                 &
                        &   scalar = edge_sources%source(1)%ggd( time_sind )% &
                        &            ion( is )%state(1)%particles,            &
                        &   b2CellData = tmpCv )
                    end if
                end if
            end do

            !! ue: Parallel Electron Velocity
            if (ue_style.eq.2) then
                call b2xpve (nx, ny, ns, qe, rza, rz2, alfx_c, sigx_c, po, &
                    &        ne, te, na, ua, ue_ext, ne_ext, ue)
            else if (ue_style.ne.2) then
              do iy=-1,ny
                do ix=-1,nx
                  if(leftix(ix,iy).ne.-2 .and. rightix(ix,iy).ne.nx+1) then
                    ue(ix,iy)=-(fch_p(ix,iy,0)/pbs(ix,iy,0)+             &
                      &         fch_p(rightix(ix,iy),rightiy(ix,iy),0)/  &
                      &         pbs(rightix(ix,iy),rightiy(ix,iy),0))/qe/2.0_IDS_real
                  elseif(leftix(ix,iy).eq.-2) then
                    ue(ix,iy)=-fch_p(rightix(ix,iy),rightiy(ix,iy),0)/   &
                      &        pbs(rightix(ix,iy),rightiy(ix,iy),0)/qe
                  elseif(rightix(ix,iy).eq.nx+1) then
                    ue(ix,iy)=-fch_p(ix,iy,0)/pbs(ix,iy,0)/qe
                  endif
                enddo
              enddo
              do is=0,ns-1
                ue(:,:)=ue(:,:)+rza(:,:,is)*ua(:,:,is)*na(:,:,is)
              enddo
              ue(:,:)=ue(:,:)+ue_ext(:,:)*ne_ext(:,:)
              ue(:,:)=ue(:,:)/ne(:,:)
            endif
            call write_cell_vector_component(                            &
                 &   vectorComponent = edge_profiles%ggd( time_sind )%   &
                 &                     electrons%velocity,               &
                 &   b2CellData = ue(:,:),                               &
                 &   vectorID = VEC_ALIGN_PARALLEL_ID )

            do is = 1, nsion
                if (is.le.fluids_list(ns-1)) then
                    js = ions_list(is)
                !! ua: Parallel ion velocity
                    call write_cell_vector_component(                           &
                        &   vectorComponent = edge_profiles%ggd( time_sind )%   &
                        &                     ion( is )%velocity,               &
                        &   b2CellData = ua(:,:,js),                            &
                        &   vectorID = VEC_ALIGN_PARALLEL_ID )
                    call write_cell_vector_component(                           &
                        &   vectorComponent = edge_profiles%ggd( time_sind )%   &
                        &                     ion( is )%state(1)%velocity,      &
                        &   b2CellData = ua(:,:,js),                            &
                        &   vectorID = VEC_ALIGN_PARALLEL_ID )
                !! wadia: Diamagnetic ion velocity
                    call write_cell_vector_component(                           &
                        &   vectorComponent = edge_profiles%ggd( time_sind )%   &
                        &                     ion( is )%state(1)%               &
                        &                     velocity_diamagnetic,             &
                        &   b2CellData = wadia(:,:,0,js),                       &
                        &   vectorID = VEC_ALIGN_POLOIDAL_ID )
                    call write_cell_vector_component(                           &
                        &   vectorComponent = edge_profiles%ggd( time_sind )%   &
                        &                     ion( is )%state(1)%               &
                        &                     velocity_diamagnetic,             &
                        &   b2CellData = wadia(:,:,1,js),                       &
                        &   vectorID = VEC_ALIGN_RADIAL_ID )
                !! vaecrb: ExB ion velocity
                    call write_cell_vector_component(                           &
                        &   vectorComponent = edge_profiles%ggd( time_sind )%   &
                        &                     ion( is )%state(1)%velocity_exb,  &
                        &   b2CellData = vaecrb(:,:,0,js),                      &
                        &   vectorID = VEC_ALIGN_POLOIDAL_ID )
                    call write_cell_vector_component(                           &
                        &   vectorComponent = edge_profiles%ggd( time_sind )%   &
                        &                     ion( is )%state(1)%velocity_exb,  &
                        &   b2CellData = vaecrb(:,:,1,js),                      &
                        &   vectorID = VEC_ALIGN_RADIAL_ID )
                !! cvsa: Ion diffusivity
                    call write_cell_vector_component(                           &
                        &   vectorComponent = edge_transport%model(1)%          &
                        &                     ggd( time_sind )%ion( is )%       &
                        &                     momentum%d,                       &
                        &   b2CellData = cvsa(:,:,0,js),                        &
                        &   vectorID = VEC_ALIGN_POLOIDAL_ID )
                    call write_cell_vector_component(                           &
                        &   vectorComponent = edge_transport%model(1)%          &
                        &                     ggd( time_sind )%ion( is )%       &
                        &                     momentum%d,                       &
                        &   b2CellData = cvsa(:,:,1,js),                        &
                        &   vectorID = VEC_ALIGN_RADIAL_ID )
                    call write_cell_vector_component(                           &
                        &   vectorComponent = edge_transport%model(1)%          &
                        &                     ggd( time_sind )%ion( is )%       &
                        &                     state(1)%momentum%d,              &
                        &   b2CellData = cvsa(:,:,0,js),                        &
                        &   vectorID = VEC_ALIGN_POLOIDAL_ID )
                    call write_cell_vector_component(                           &
                        &   vectorComponent = edge_transport%model(1)%          &
                        &                     ggd( time_sind )%ion( is )%       &
                        &                     state(1)%momentum%d,              &
                        &   b2CellData = cvsa(:,:,1,js),                        &
                        &   vectorID = VEC_ALIGN_RADIAL_ID )
                !! fmo: Ion momentum flux
                    do ix = -1, nx
                      do iy = -1, ny
                        if (fmo(ix,iy,0,js).ne.0.0_R8) then
                          tmpFace(ix,iy,0) = fmo(ix,iy,0,js)/gs(ix,iy,0)/qc(ix,iy)
                        else
                          tmpFace(ix,iy,0) = 0.0_R8
                        end if
                        if (fmo(ix,iy,1,js).ne.0.0_R8) then
                          tmpFace(ix,iy,1) = fmo(ix,iy,1,js)/gs(ix,iy,1)
                        else
                          tmpFace(ix,iy,1) = 0.0_R8
                        end if
                      end do
                    end do
                    call write_cell_vector_component(                           &
                        &   vectorComponent = edge_transport%model(1)%          &
                        &                     ggd( time_sind )%ion( is )%       &
                        &                     momentum%flux,                    &
                        &   b2CellData = tmpFace(:,:,0),                        &
                        &   vectorID = VEC_ALIGN_PARALLEL_ID )
                    call write_cell_vector_component(                           &
                        &   vectorComponent = edge_transport%model(1)%          &
                        &                     ggd( time_sind )%ion( is )%       &
                        &                     state(1)%momentum%flux,           &
                        &   b2CellData = tmpFace(:,:,0),                        &
                        &   vectorID = VEC_ALIGN_PARALLEL_ID )
                !! fllimvisc: Ion momentum transport flux limit
                    call write_cell_vector_component(                           &
                        &   vectorComponent = edge_transport%model(1)%          &
                        &                     ggd( time_sind )%ion( is )%       &
                        &                     momentum%flux_limiter,            &
                        &   b2CellData = fllimvisc(:,:,js),                     &
                        &   vectorID = VEC_ALIGN_PARALLEL_ID )
                    call write_cell_vector_component(                           &
                        &   vectorComponent = edge_transport%model(1)%          &
                        &                     ggd( time_sind )%ion( is )%       &
                        &                     state(1)%momentum%flux_limiter,   &
                        &   b2CellData = fllimvisc(:,:,js),                     &
                        &   vectorID = VEC_ALIGN_PARALLEL_ID )
                !! smo: Ion parallel momentum sources
                    tmpCv(:,:) = ( smo(:,:,0,js) +                              &
                        &          smo(:,:,1,js) * ua(:,:,js) +                 &
                        &          smo(:,:,2,js) * roxa(:,:,js) +               &
                        &          smo(:,:,3,js) * roxa(:,:,js) * ua(:,:,js) ) / vol(:,:)
                    call write_cell_vector_component(                               &
                        &   vectorComponent = edge_sources%source(1)%               &
                        &                     ggd( time_sind )%ion( is )%momentum,  &
                        &   b2CellData = tmpCv,                                     &
                        &   vectorID = VEC_ALIGN_PARALLEL_ID )
                    call write_cell_vector_component(                               &
                        &   vectorComponent = edge_sources%source(1)%               &
                        &                     ggd( time_sind )%ion( is )%           &
                        &                     state(1)%momentum,                    &
                        &   b2CellData = tmpCv,                                     &
                        &   vectorID = VEC_ALIGN_PARALLEL_ID )
                    tmpCv(:,:) = ext_smo(:,:,js) / vol(:,:)
                    call write_cell_vector_component(                               &
                        &   vectorComponent = edge_sources%source(2)%               &
                        &                     ggd( time_sind )%ion( is )%momentum,  &
                        &   b2CellData = tmpCv,                                     &
                        &   vectorID = VEC_ALIGN_PARALLEL_ID )
                    call write_cell_vector_component(                               &
                        &   vectorComponent = edge_sources%source(2)%               &
                        &                     ggd( time_sind )%ion( is )%           &
                        &                     state(1)%momentum,                    &
                        &   b2CellData = tmpCv,                                     &
                        &   vectorID = VEC_ALIGN_PARALLEL_ID )
                    tmpCv(:,:) = ( b2stbc_smo(:,:,js) +                             &
                        &          b2stbm_smo(:,:,js) ) / vol(:,:)
                    call write_cell_vector_component(                               &
                        &   vectorComponent = edge_sources%source(3)%               &
                        &                     ggd( time_sind )%ion( is )%momentum,  &
                        &   b2CellData = tmpCv,                                     &
                        &   vectorID = VEC_ALIGN_PARALLEL_ID )
                    call write_cell_vector_component(                               &
                        &   vectorComponent = edge_sources%source(3)%               &
                        &                     ggd( time_sind )%ion( is )%           &
                        &                     state(1)%momentum,                    &
                        &   b2CellData = tmpCv,                                     &
                        &   vectorID = VEC_ALIGN_PARALLEL_ID )
                    tmpCv(:,:) = ( smodt(:,:,0,js) +                                &
                        &          smodt(:,:,1,js) * ua(:,:,js) +                   &
                        &          smodt(:,:,2,js) * roxa(:,:,js) +                 &
                        &          smodt(:,:,3,js) * roxa(:,:,js) * ua(:,:,js) ) / vol(:,:)
                    call write_cell_vector_component(                               &
                        &   vectorComponent = edge_sources%source(4)%               &
                        &                     ggd( time_sind )%ion( is )%momentum,  &
                        &   b2CellData = tmpCv,                                     &
                        &   vectorID = VEC_ALIGN_PARALLEL_ID )
                    call write_cell_vector_component(                               &
                        &   vectorComponent = edge_sources%source(4)%               &
                        &                     ggd( time_sind )%ion( is )%           &
                        &                     state(1)%momentum,                    &
                        &   b2CellData = tmpCv,                                     &
                        &   vectorID = VEC_ALIGN_PARALLEL_ID )
                    if (use_eirene.ne.0 .and. balance_netcdf.ne.0) then
                        tmpCv = 0.0_IDS_real
                        do istrai = 1, size( eirene_mc_mapl_smo_bal, 4)
                            tmpCv(:,:) = tmpCv(:,:) + eirene_mc_mapl_smo_bal(:,:,js,istrai)
                        end do
                        tmpCv(:,:) = tmpCv(:,:) / vol(:,:)
                        call write_cell_vector_component(                           &
                            &   vectorComponent = edge_sources%source(5)%           &
                            &            ggd( time_sind )%ion( is )%momentum,       &
                            &   b2CellData = tmpCv,                                 &
                            &   vectorID = VEC_ALIGN_PARALLEL_ID )
                        call write_cell_vector_component(                           &
                            &   vectorComponent = edge_sources%source(5)%           &
                            &            ggd( time_sind )%ion( is )%                &
                            &            state(1)%momentum,                         &
                            &   b2CellData = tmpCv,                                 &
                            &   vectorID = VEC_ALIGN_PARALLEL_ID )
                        tmpCv = 0.0_IDS_real
                        do istrai = 1, size( eirene_mc_mmpl_smo_bal, 4)
                            tmpCv(:,:) = tmpCv(:,:) + eirene_mc_mmpl_smo_bal(:,:,js,istrai)
                        end do
                        tmpCv(:,:) = tmpCv(:,:) / vol(:,:)
                        call write_cell_vector_component(                           &
                            &   vectorComponent = edge_sources%source(6)%           &
                            &            ggd( time_sind )%ion( is )%momentum,       &
                            &   b2CellData = tmpCv,                                 &
                            &   vectorID = VEC_ALIGN_PARALLEL_ID )
                        call write_cell_vector_component(                           &
                            &   vectorComponent = edge_sources%source(6)%           &
                            &            ggd( time_sind )%ion( is )%                &
                            &            state(1)%momentum,                         &
                            &   b2CellData = tmpCv,                                 &
                            &   vectorID = VEC_ALIGN_PARALLEL_ID )
                    else
                        tmpCv(:,:) = b2stbr_smo(:,:,js) / vol(:,:)
                        call write_cell_vector_component(                           &
                            &   vectorComponent = edge_sources%source(5)%           &
                            &            ggd( time_sind )%ion( is )%momentum,       &
                            &   b2CellData = tmpCv,                                 &
                            &   vectorID = VEC_ALIGN_PARALLEL_ID )
                        call write_cell_vector_component(                           &
                            &   vectorComponent = edge_sources%source(5)%           &
                            &            ggd( time_sind )%ion( is )%                &
                            &            state(1)%momentum,                         &
                            &   b2CellData = tmpCv,                                 &
                            &   vectorID = VEC_ALIGN_PARALLEL_ID )
                    end if
                    tmpCv(:,:) = rsamo(:,:,js) / vol(:,:)
                    call write_cell_vector_component(                               &
                        &   vectorComponent = edge_sources%source(7)%               &
                        &                     ggd( time_sind )%ion( is )%momentum,  &
                        &   b2CellData = tmpCv,                                     &
                        &   vectorID = VEC_ALIGN_PARALLEL_ID )
                    call write_cell_vector_component(                               &
                        &   vectorComponent = edge_sources%source(7)%               &
                        &                     ggd( time_sind )%ion( is )%           &
                        &                     state(1)%momentum,                    &
                        &   b2CellData = tmpCv,                                     &
                        &   vectorID = VEC_ALIGN_PARALLEL_ID )
                    tmpCv(:,:) = rramo(:,:,js) / vol(:,:)
                    call write_cell_vector_component(                               &
                        &   vectorComponent = edge_sources%source(8)%               &
                        &                     ggd( time_sind )%ion( is )%momentum,  &
                        &   b2CellData = tmpCv,                                     &
                        &   vectorID = VEC_ALIGN_PARALLEL_ID )
                    call write_cell_vector_component(                               &
                        &   vectorComponent = edge_sources%source(8)%               &
                        &                     ggd( time_sind )%ion( is )%           &
                        &                     state(1)%momentum,                    &
                        &   b2CellData = tmpCv,                                     &
                        &   vectorID = VEC_ALIGN_PARALLEL_ID )
                    tmpCv(:,:) = rcxmo(:,:,js) / vol(:,:)
                    call write_cell_vector_component(                               &
                        &   vectorComponent = edge_sources%source(9)%               &
                        &                     ggd( time_sind )%ion( is )%momentum,  &
                        &   b2CellData = tmpCv,                                     &
                        &   vectorID = VEC_ALIGN_PARALLEL_ID )
                    call write_cell_vector_component(                               &
                        &   vectorComponent = edge_sources%source(9)%               &
                        &                     ggd( time_sind )%ion( is )%           &
                        &                     state(1)%momentum,                    &
                        &   b2CellData = tmpCv,                                     &
                        &   vectorID = VEC_ALIGN_PARALLEL_ID )
                else
                    js = is - fluids_list(ns-1)
                end if
            end do

            !! te: Electron Temperature
            tmpCv(:,:) = te(:,:)/qe
            call write_quantity(                                        &
                &   val = edge_profiles%ggd( time_sind )%electrons%     &
                &         temperature,                                  &
                &   value = tmpCv,                                      &
                &   time_sind = time_sind )
            call write_face_scalar(                                     &
                &   val = edge_transport%model(1)%ggd( time_sind )%     &
                &         electrons%energy%d,                           &
                &   value = chce,                                       &
                &   time_sind = time_sind )
            call write_face_scalar(                                     &
                &   val = edge_transport%model(1)%ggd( time_sind )%     &
                &         electrons%energy%v,                           &
                &   value = chve,                                       &
                &   time_sind = time_sind )
            do ix = -1, nx
              do iy = -1, ny
                if (fhe(ix,iy,0).ne.0.0_R8) then
                  tmpFace(ix,iy,0) = fhe(ix,iy,0)/gs(ix,iy,0)/qc(ix,iy)
                else
                  tmpFace(ix,iy,0) = 0.0_R8
                end if
                if (fhe(ix,iy,1).ne.0.0_R8) then
                  tmpFace(ix,iy,1) = fhe(ix,iy,1)/gs(ix,iy,1)
                else
                  tmpFace(ix,iy,1) = 0.0_R8
                end if
              end do
            end do
            call write_face_scalar(                                     &
                &   val = edge_transport%model(1)%ggd( time_sind )%     &
                &         electrons%energy%flux,                        &
                &   value = tmpFace,                                    &
                &   time_sind = time_sind )
            call write_cell_scalar(                                     &
                &   scalar = edge_transport%model(1)%ggd( time_sind )%  &
                &            electrons%energy%flux_limiter,             &
                &   b2CellData = fllime )
            tmpCv(:,:) = ( she(:,:,0) + she(:,:,1) * te(:,:) +          &
                &          she(:,:,2) * ne(:,:) +                       &
                &          she(:,:,3) * te(:,:) * ne(:,:) ) / vol(:,:)
            call write_cell_scalar(                                     &
                &   scalar = edge_sources%source(1)%ggd( time_sind )%   &
                &            electrons%energy,                          &
                &   b2CellData = tmpCv )
            tmpCv(:,:) = ext_she(:,:) / vol(:,:)
            call write_cell_scalar(                                     &
                &   scalar = edge_sources%source(2)%ggd( time_sind )%   &
                &            electrons%energy,                          &
                &   b2CellData = tmpCv )
            tmpCv(:,:) = ( b2stbc_she(:,:) + b2stbm_she(:,:) ) / vol(:,:)
            call write_cell_scalar(                                     &
                &   scalar = edge_sources%source(3)%ggd( time_sind )%   &
                &            electrons%energy,                          &
                &   b2CellData = tmpCv )
            tmpCv(:,:) = ( shedt(:,:,0) + shedt(:,:,1) * te(:,:) +      &
                &          shedt(:,:,2) * ne(:,:) +                     &
                &          shedt(:,:,3) * te(:,:) * ne(:,:) ) / vol(:,:)
            call write_cell_scalar(                                     &
                &   scalar = edge_sources%source(4)%ggd( time_sind )%   &
                &            electrons%energy,                          &
                &   b2CellData = tmpCv )
            if (use_eirene.ne.0 .and. balance_netcdf.ne.0) then
                tmpCv = 0.0_IDS_real
                do is = 1, size( eirene_mc_eael_she_bal, 3)
                    tmpCv(:,:) = tmpCv(:,:) + eirene_mc_eael_she_bal(:,:,is)
                end do
                tmpCv(:,:) = tmpCv(:,:) / vol(:,:)
                call write_cell_scalar(                                     &
                    &   scalar = edge_sources%source(5)%ggd( time_sind )%   &
                    &            electrons%energy,                          &
                    &   b2CellData = tmpCv )
                tmpCv = 0.0_IDS_real
                do is = 1, size( eirene_mc_emel_she_bal, 3)
                    tmpCv(:,:) = tmpCv(:,:) + eirene_mc_emel_she_bal(:,:,is)
                end do
                tmpCv(:,:) = tmpCv(:,:) / vol(:,:)
                call write_cell_scalar(                                     &
                    &   scalar = edge_sources%source(6)%ggd( time_sind )%   &
                    &            electrons%energy,                          &
                    &   b2CellData = tmpCv )
            else
                tmpCv(:,:) = b2stbr_she(:,:) / vol(:,:)
                call write_cell_scalar(                                     &
                    &   scalar = edge_sources%source(5)%ggd( time_sind )%   &
                    &            electrons%energy,                          &
                    &   b2CellData = tmpCv )
            end if
            if (balance_netcdf.ne.0) then
                tmpCv(:,:) = b2stel_she_ion_bal(:,:) / vol(:,:)
                call write_cell_scalar(                                     &
                    &   scalar = edge_sources%source(7)%ggd( time_sind )%   &
                    &            electrons%energy,                          &
                    &   b2CellData = tmpCv )
                tmpCv(:,:) = b2stel_she_rec_bal(:,:) / vol(:,:)
                call write_cell_scalar(                                     &
                    &   scalar = edge_sources%source(8)%ggd( time_sind )%   &
                    &            electrons%energy,                          &
                    &   b2CellData = tmpCv )
                tmpCv(:,:) = -b2npht_shei_bal(:,:) / vol(:,:)
                call write_cell_scalar(                                     &
                    &   scalar = edge_sources%source(10)%ggd( time_sind )%  &
                    &            electrons%energy,                          &
                    &   b2CellData = tmpCv )
            end if
            tmpCv(:,:) = b2sihs_joule(:,:) / vol(:,:)
            call write_cell_scalar(                                         &
                &   scalar = edge_sources%source(11)%ggd( time_sind )%      &
                &            electrons%energy,                              &
                &   b2CellData = tmpCv )
            tmpCv(:,:) = 0.0_IDS_real
            do is = 0, ns-1
              tmpCv(:,:) = tmpCv(:,:) + rqrad(:,:,is)
            end do
#ifdef B25_EIRENE
            do is = 1, natmi
              tmpCv(-1:nx,-1:ny) = tmpCv(-1:nx,-1:ny) + eneutrad(0:nx+1,0:ny+1,is,0)
            end do
            do is = 1, nmoli
              tmpCv(-1:nx,-1:ny) = tmpCv(-1:nx,-1:ny) + emolrad(0:nx+1,0:ny+1,is,0)
            end do
            do is = 1, nioni
              tmpCv(-1:nx,-1:ny) = tmpCv(-1:nx,-1:ny) + eionrad(0:nx+1,0:ny+1,is,0)
            end do
#endif
            tmpCv(:,:) = tmpCv(:,:) / vol(:,:)
            call write_cell_scalar(                                         &
                &   scalar = edge_sources%source(12)%ggd( time_sind )%      &
                &            electrons%energy,                              &
                &   b2CellData = tmpCv )

            !! pe: Electron pressure
            call b2xppe( nx, ny, ne, te, pe)
            call write_quantity(                                      &
                &   val = edge_profiles%ggd( time_sind )%electrons%   &
                &         pressure,                                   &
                &   value = pe,                                       &
                &   time_sind = time_sind )

            !! ti: (Common) Ion Temperature
            tmpCv(:,:) = ti(:,:)/qe
            call write_quantity(                                      &
                &   val = edge_profiles%ggd( time_sind )%t_i_average, &
                &   value = tmpCv,                                    &
                &   time_sind = time_sind )
            call write_face_scalar(                                   &
                &   val = edge_transport%model(1)%ggd( time_sind )%   &
                &         total_ion_energy%d,                         &
                &   value = chci,                                     &
                &   time_sind = time_sind )
            call write_face_scalar(                                   &
                 &   val = edge_transport%model(1)%ggd( time_sind )%  &
                 &         total_ion_energy%v,                        &
                 &   value = chvi,                                    &
                 &   time_sind = time_sind )
            !! fhi : Ion heat flux
            do ix = -1, nx
              do iy = -1, ny
                if (fhi(ix,iy,0).ne.0.0_R8) then
                  tmpFace(ix,iy,0) = fhi(ix,iy,0)/gs(ix,iy,0)/qc(ix,iy)
                else
                  tmpFace(ix,iy,0) = 0.0_R8
                end if
                if (fhi(ix,iy,1).ne.0.0_R8) then
                  tmpFace(ix,iy,1) = fhi(ix,iy,1)/gs(ix,iy,1)
                else
                  tmpFace(ix,iy,1) = 0.0_R8
                end if
              end do
            end do
            call write_face_scalar(                                     &
                &   val = edge_transport%model(1)%ggd( time_sind )%     &
                &         total_ion_energy%flux,                        &
                &   value = tmpFace,                                    &
                &   time_sind = time_sind )
            call write_cell_scalar(                                     &
                &   scalar = edge_transport%model(1)%ggd( time_sind )%  &
                &            total_ion_energy%flux_limiter,             &
                &   b2CellData = fllimi )
            !! Ion energy sources
            tmpCv(:,:) = ( shi(:,:,0) + shi(:,:,1) * ti(:,:) +          &
                &          shi(:,:,2) * ni(:,:,0) +                     &
                &          shi(:,:,3) * ni(:,:,0) * ti(:,:) ) / vol(:,:)
            call write_cell_scalar(                                     &
                &   scalar = edge_sources%source(1)%ggd( time_sind )%   &
                &            total_ion_energy,                          &
                &   b2CellData = tmpCv )
            tmpCv(:,:) = ext_shi(:,:) / vol(:,:)
            call write_cell_scalar(                                     &
                &   scalar = edge_sources%source(2)%ggd( time_sind )%   &
                &            total_ion_energy,                          &
                &   b2CellData = tmpCv )
            tmpCv(:,:) = ( b2stbc_shi(:,:) + b2stbm_shi(:,:) ) / vol(:,:)
            call write_cell_scalar(                                     &
                &   scalar = edge_sources%source(3)%ggd( time_sind )%   &
                &            total_ion_energy,                          &
                &   b2CellData = tmpCv )
            tmpCv(:,:) = ( shidt(:,:,0) + shidt(:,:,1) * ti(:,:) +      &
                &          shidt(:,:,2) * ni(:,:,0) +                   &
                &          shidt(:,:,3) * ni(:,:,0) * ti(:,:) ) / vol(:,:)
            call write_cell_scalar(                                     &
                &   scalar = edge_sources%source(4)%ggd( time_sind )%   &
                &            total_ion_energy,                          &
                &   b2CellData = tmpCv )
            if (use_eirene.ne.0 .and. balance_netcdf.ne.0) then
                tmpCv = 0.0_IDS_real
                do is = 1, size( eirene_mc_eapl_shi_bal, 3)
                    tmpCv(:,:) = tmpCv(:,:) + eirene_mc_eapl_shi_bal(:,:,is)
                end do
                tmpCv(:,:) = tmpCv(:,:) / vol(:,:)
                call write_cell_scalar(                                     &
                    &   scalar = edge_sources%source(5)%ggd( time_sind )%   &
                    &            total_ion_energy,                          &
                    &   b2CellData = tmpCv )
                tmpCv = 0.0_IDS_real
                do is = 1, size( eirene_mc_empl_shi_bal, 3)
                    tmpCv(:,:) = tmpCv(:,:) + eirene_mc_empl_shi_bal(:,:,is)
                end do
                tmpCv(:,:) = tmpCv(:,:) / vol(:,:)
                call write_cell_scalar(                                     &
                    &   scalar = edge_sources%source(6)%ggd( time_sind )%   &
                    &            total_ion_energy,                          &
                    &   b2CellData = tmpCv )
            else
                tmpCv(:,:) = b2stbr_shi(:,:) / vol(:,:)
                call write_cell_scalar(                                     &
                    &   scalar = edge_sources%source(5)%ggd( time_sind )%   &
                    &            total_ion_energy,                          &
                    &   b2CellData = tmpCv )
            end if
            if (balance_netcdf.ne.0) then
                tmpCv(:,:) = b2stel_shi_ion_bal(:,:) / vol(:,:)
                call write_cell_scalar(                                     &
                    &   scalar = edge_sources%source(7)%ggd( time_sind )%   &
                    &            total_ion_energy,                          &
                    &   b2CellData = tmpCv )
                tmpCv(:,:) = b2stel_shi_rec_bal(:,:) / vol(:,:)
                call write_cell_scalar(                                     &
                    &   scalar = edge_sources%source(8)%ggd( time_sind )%   &
                    &            total_ion_energy,                          &
                    &   b2CellData = tmpCv )
                tmpCv(:,:) = b2npht_shei_bal(:,:) / vol(:,:)
                call write_cell_scalar(                                     &
                    &   scalar = edge_sources%source(10)%ggd( time_sind )%  &
                    &            total_ion_energy,                          &
                    &   b2CellData = tmpCv )
            end if
            do is = 1, nsion
                if (is.le.fluids_list(ns-1)) then
                    js = ions_list(is)
            !! Ion energy sources resolved by species
                    tmpCv(:,:) = rsahi(:,:,js) / vol(:,:)
                    call write_cell_scalar(                                     &
                        &   scalar = edge_sources%source(7)%ggd( time_sind )%   &
                        &            ion( is )%energy,                          &
                        &   b2CellData = tmpCv )
                    call write_cell_scalar(                                     &
                        &   scalar = edge_sources%source(7)%ggd( time_sind )%   &
                        &            ion( is )%state(1)%energy,                 &
                        &   b2CellData = tmpCv )
                    tmpCv(:,:) = rrahi(:,:,js) / vol(:,:)
                    call write_cell_scalar(                                     &
                        &   scalar = edge_sources%source(8)%ggd( time_sind )%   &
                        &            ion( is )%energy,                          &
                        &   b2CellData = tmpCv )
                    call write_cell_scalar(                                     &
                        &   scalar = edge_sources%source(8)%ggd( time_sind )%   &
                        &            ion( is )%state(1)%energy,                 &
                        &   b2CellData = tmpCv )
                    tmpCv(:,:) = rcxhi(:,:,js) / vol(:,:)
                    call write_cell_scalar(                                     &
                        &   scalar = edge_sources%source(9)%ggd( time_sind )%   &
                        &            ion( is )%energy,                          &
                        &   b2CellData = tmpCv )
                    call write_cell_scalar(                                     &
                        &   scalar = edge_sources%source(9)%ggd( time_sind )%   &
                        &            ion( is )%state(1)%energy,                 &
                        &   b2CellData = tmpCv )
                else
                    js = is - fluids_list(ns-1)
#ifdef B25_EIRENE
            !! Test ion temperature
                    tmpCv(-1:nx,-1:ny) = tib2(0:nx+1,0:ny+1,js,1)/qe
                    call write_quantity(                                        &
                        &   val = edge_profiles%ggd( time_sind )%ion( is )%     &
                        &         temperature,                                  &
                        &   value = tmpCv,                                      &
                        &   time_sind = time_sind )
                    call write_quantity(                                        &
                        &   val = edge_profiles%ggd( time_sind )%ion( is )%     &
                        &         state(1)%temperature,                         &
                        &   value = tmpCv,                                      &
                        &   time_sind = time_sind )
#endif
                end if
            end do

            do is = 1, nsion
                if (is.le.fluids_list(ns-1)) then
                    js = ions_list(is)
                !! pb : Ion pressure
                    call b2xppb( nx, ny, rza(:,:,js), na(:,:,js), te, ti, pb)
                    call write_quantity(                                      &
                        &   val = edge_profiles%ggd( time_sind )%ion( is )%   &
                        &         pressure,                                   &
                        &   value = pb,                                       &
                        &   time_sind = time_sind )
                    call write_quantity(                                      &
                        &   val = edge_profiles%ggd( time_sind )%ion( is )%   &
                        &         state(1)%pressure,                          &
                        &   value = pb,                                       &
                        &   time_sind = time_sind )
                !! Kinetic energy density
                    tmpCv(:,:) = 0.5_IDS_real*am(js)*mp*(ua(:,:,js)**2)*na(:,:,js)
                    call write_quantity(                                      &
                        &   val = edge_profiles%ggd( time_sind )%ion( is )%   &
                        &         energy_density_kinetic,                     &
                        &   value = tmpCv,                                    &
                        &   time_sind = time_sind )
                    call write_quantity(                                      &
                        &   val = edge_profiles%ggd( time_sind )%ion( is )%   &
                        &         state(1)%energy_density_kinetic,            &
                        &   value = tmpCv,                                    &
                        &   time_sind = time_sind )
                !! Average charge
                    call write_quantity(                                      &
                        &   val = edge_profiles%ggd( time_sind )%ion( is )%   &
                        &         state(1)%z_average,                         &
                        &   value = rza(:,:,js),                              &
                        &   time_sind = time_sind )
                !! Average square charge
                    call write_quantity(                                      &
                        &   val = edge_profiles%ggd( time_sind )%ion( is )%   &
                        &         state(1)%z_square_average,                  &
                        &   value = rz2(:,:,js),                              &
                        &   time_sind = time_sind )
                !! Ionization potential
                    call write_quantity(                                      &
                        &   val = edge_profiles%ggd( time_sind )%ion( is )%   &
                        &         state(1)%ionisation_potential,              &
                        &   value = rpt(:,:,js),                              &
                        &   time_sind = time_sind )
                else
                    js = is - fluids_list(ns-1)
#ifdef B25_EIRENE
                !! Test ion pressure
                    tmpCv(-1:nx,-1:ny) = dib2(0:nx+1,0:ny+1,js,1)*tib2(0:nx+1,0:ny+1,js,1)
                    call write_quantity(                                      &
                        &   val = edge_profiles%ggd( time_sind )%ion( is )%   &
                        &         pressure,                                   &
                        &   value = tmpCv,                                    &
                        &   time_sind = time_sind )
                    call write_quantity(                                      &
                        &   val = edge_profiles%ggd( time_sind )%ion( is )%   &
                        &         state(1)%pressure,                          &
                        &   value = tmpCv,                                    &
                        &   time_sind = time_sind )
                !! Average charge
                    tmpCv(:,:) = nchrgi( js )
                    call write_quantity(                                      &
                        &   val = edge_profiles%ggd( time_sind )%ion( is )%   &
                        &         state(1)%z_average,                         &
                        &   value = tmpCv,                                    &
                        &   time_sind = time_sind )
                !! Average square charge
                    tmpCv(:,:) = nchrgi( js )**2
                    call write_quantity(                                      &
                        &   val = edge_profiles%ggd( time_sind )%ion( is )%   &
                        &         state(1)%z_square_average,                  &
                        &   value = tmpCv ,                                   &
                        &   time_sind = time_sind )
                !! Radiation source
                    tmpCv(-1:nx,-1:ny) = eionrad(0:nx+1,0:ny+1,js,0) / vol(-1:nx,-1:ny)
                    call write_cell_scalar(                                   &
                        &   scalar = edge_sources%source(12)%ggd( time_sind )%      &
                &            electrons%energy,                              &
                &   b2CellData = tmpCv )
#endif
               end if
            end do

            !! Ion energy sources per species
            do is = 1, nsion
                if (is.le.fluids_list(ns-1)) then
                    js = ions_list(is)
                    tmpCv(:,:) = rsahi(:,:,js) / vol(:,:)
                    call write_cell_scalar(                                     &
                        &   scalar = edge_sources%source(7)%ggd( time_sind )%   &
                        &            ion( is )%energy,                          &
                        &   b2CellData = tmpCv )
                    call write_cell_scalar(                                     &
                        &   scalar = edge_sources%source(7)%ggd( time_sind )%   &
                        &            ion( is )%state(1)%energy,                 &
                        &   b2CellData = tmpCv )
                    tmpCv(:,:) = rrahi(:,:,js) / vol(:,:)
                    call write_cell_scalar(                                     &
                        &   scalar = edge_sources%source(8)%ggd( time_sind )%   &
                        &            ion( is )%energy,                          &
                        &   b2CellData = tmpCv )
                    call write_cell_scalar(                                     &
                        &   scalar = edge_sources%source(8)%ggd( time_sind )%   &
                        &            ion( is )%state(1)%energy,                 &
                        &   b2CellData = tmpCv )
                    tmpCv(:,:) = rcxhi(:,:,js) / vol(:,:)
                    call write_cell_scalar(                                     &
                        &   scalar = edge_sources%source(9)%ggd( time_sind )%   &
                        &            ion( is )%energy,                          &
                        &   b2CellData = tmpCv )
                    call write_cell_scalar(                                     &
                        &   scalar = edge_sources%source(9)%ggd( time_sind )%   &
                        &            ion( is )%state(1)%energy,                 &
                        &   b2CellData = tmpCv )
                else
                    js = is - fluids_list(ns-1)
                    tmpCv(-1:nx,-1:ny) = tib2(0:nx+1,0:ny+1,js,1)/qe
                    call write_quantity(                                        &
                        &   val = edge_profiles%ggd( time_sind )%ion( is )%     &
                        &         temperature,                                  &
                        &   value = tmpCv,                                      &
                        &   time_sind = time_sind )
                    call write_quantity(                                        &
                        &   val = edge_profiles%ggd( time_sind )%ion( is )%     &
                        &         state(1)%temperature,                         &
                        &   value = tmpCv,                                      &
                        &   time_sind = time_sind )
                end if
            end do

            do is = 1, nsion
                if (is.le.fluids_list(ns-1)) then
                    js = ions_list(is)
                !! pb : Ion pressure
                    call b2xppb( nx, ny, rza(:,:,js), na(:,:,js), te, ti, pb)
                    call write_quantity(                                      &
                        &   val = edge_profiles%ggd( time_sind )%ion( is )%   &
                        &         pressure,                                   &
                        &   value = pb,                                       &
                        &   time_sind = time_sind )
                    call write_quantity(                                      &
                        &   val = edge_profiles%ggd( time_sind )%ion( is )%   &
                        &         state(1)%pressure,                          &
                        &   value = pb,                                       &
                        &   time_sind = time_sind )
                !! Kinetic energy density
                    tmpCv(:,:) = 0.5_IDS_real*am(js)*mp*(ua(:,:,js)**2)*na(:,:,js)
                    call write_quantity(                                      &
                        &   val = edge_profiles%ggd( time_sind )%ion( is )%   &
                        &         energy_density_kinetic,                     &
                        &   value = tmpCv,                                    &
                        &   time_sind = time_sind )
                    call write_quantity(                                      &
                        &   val = edge_profiles%ggd( time_sind )%ion( is )%   &
                        &         state(1)%energy_density_kinetic,            &
                        &   value = tmpCv,                                    &
                        &   time_sind = time_sind )
                !! Average charge
                    call write_quantity(                                      &
                        &   val = edge_profiles%ggd( time_sind )%ion( is )%   &
                        &         state(1)%z_average,                         &
                        &   value = rza(:,:,js),                              &
                        &   time_sind = time_sind )
                !! Average square charge
                    call write_quantity(                                      &
                        &   val = edge_profiles%ggd( time_sind )%ion( is )%   &
                        &         state(1)%z_square_average,                  &
                        &   value = rz2(:,:,js),                              &
                        &   time_sind = time_sind )
                !! Ionization potential
                    call write_quantity(                                      &
                        &   val = edge_profiles%ggd( time_sind )%ion( is )%   &
                        &         state(1)%ionisation_potential,              &
                        &   value = rpt(:,:,js),                              &
                        &   time_sind = time_sind )
                else
                    js = is - fluids_list(ns-1)
#ifdef B25_EIRENE
                !! Test ion pressure
                    tmpCv(-1:nx,-1:ny) = dib2(0:nx+1,0:ny+1,js,1)*tib2(0:nx+1,0:ny+1,js,1)
                    call write_quantity(                                      &
                        &   val = edge_profiles%ggd( time_sind )%ion( is )%   &
                        &         pressure,                                   &
                        &   value = tmpCv,                                    &
                        &   time_sind = time_sind )
                    call write_quantity(                                      &
                        &   val = edge_profiles%ggd( time_sind )%ion( is )%   &
                        &         state(1)%pressure,                          &
                        &   value = tmpCv,                                    &
                        &   time_sind = time_sind )
                !! Average charge
                    tmpCv(:,:) = nchrgi( js )
                    call write_quantity(                                      &
                        &   val = edge_profiles%ggd( time_sind )%ion( is )%   &
                        &         state(1)%z_average,                         &
                        &   value = tmpCv,                                    &
                        &   time_sind = time_sind )
                !! Average square charge
                    tmpCv(:,:) = nchrgi( js )**2
                    call write_quantity(                                      &
                        &   val = edge_profiles%ggd( time_sind )%ion( is )%   &
                        &         state(1)%z_square_average,                  &
                        &   value = tmpCv ,                                   &
                        &   time_sind = time_sind )
                !! Radiation source
                    tmpCv(-1:nx,-1:ny) = eionrad(0:nx+1,0:ny+1,js,0)
                    call write_cell_scalar(                                   &
                        &   scalar = edge_sources%source(12)%                 &
                        &            ggd( time_sind )%ion( is )%energy,       &
                        &   b2CellData = tmpCv )
                    call write_cell_scalar(                                   &
                        &   scalar = edge_sources%source(12)%                 &
                        &            ggd( time_sind )%ion( is )%              &
                        &            state(1)%energy,                         &
                        &   b2CellData = tmpCv )
#endif
                end if
            end do

            !! ni/ne: Total ion density over electron density
            tmpCv(:,:) = ni(:,:,1)/ne(:,:)
            call write_quantity(                                             &
                &   val = edge_profiles%ggd( time_sind )%n_i_total_over_n_e, &
                &   value = tmpCv,                                           &
                &   time_sind = time_sind )

            !! Zeff
            call write_quantity(                                             &
                &   val = edge_profiles%ggd( time_sind )%zeff,               &
                &   value = zeff,                                            &
                &   time_sind = time_sind )

            !! pz: Thermal plasma pressure (electrons+ions)
            call write_quantity(                                             &
                &   val = edge_profiles%ggd( time_sind )%pressure_thermal,   &
                &   value = pz,                                              &
                &   time_sind = time_sind )

            !! fchanml: Anomalous current
            call b2tanml (nx, ny, ns, vol, hx, hy, csig_an, po, fchanml)
            call write_cell_vector_component(                                &
                &   vectorComponent = edge_profiles%ggd( time_sind )%        &
                &                     j_anomalous,                           &
                &   b2CellData = fchanml(:,:,0),                             &
                &   vectorID = VEC_ALIGN_POLOIDAL_ID )
            call write_cell_vector_component(                                &
                &   vectorComponent = edge_profiles%ggd( time_sind )%        &
                &                     j_anomalous,                           &
                &   b2CellData = fchanml(:,:,1),                             &
                &   vectorID = VEC_ALIGN_RADIAL_ID )

            !! fchinert: Inertial current
            call write_cell_vector_component(                                &
                &   vectorComponent = edge_profiles%ggd( time_sind )%        &
                &                     j_inertial,                            &
                &   b2CellData = fchinert(:,:,0),                            &
                &   vectorID = VEC_ALIGN_POLOIDAL_ID )
            call write_cell_vector_component(                                &
                &   vectorComponent = edge_profiles%ggd( time_sind )%        &
                &                     j_inertial,                            &
                &   b2CellData = fchinert(:,:,1),                            &
                &   vectorID = VEC_ALIGN_RADIAL_ID )

            !! fchin: Ion-neutral friction current
            call write_cell_vector_component(                                &
                &   vectorComponent = edge_profiles%ggd( time_sind )%        &
                &                     j_ion_neutral_friction,                &
                &   b2CellData = fchin(:,:,0),                               &
                &   vectorID = VEC_ALIGN_POLOIDAL_ID )
            call write_cell_vector_component(                                &
                &   vectorComponent = edge_profiles%ggd( time_sind )%        &
                &                     j_ion_neutral_friction,                &
                &   b2CellData = fchin(:,:,1),                               &
                &   vectorID = VEC_ALIGN_RADIAL_ID )

            !! fchvispar: Parallel viscosity current
            call write_cell_vector_component(                                &
                &   vectorComponent = edge_profiles%ggd( time_sind )%        &
                &                     j_parallel_viscosity,                  &
                &   b2CellData = fchvispar(:,:,0),                           &
                &   vectorID = VEC_ALIGN_POLOIDAL_ID )
            call write_cell_vector_component(                                &
                &   vectorComponent = edge_profiles%ggd( time_sind )%        &
                &                     j_parallel_viscosity,                  &
                &   b2CellData = fchvispar(:,:,1),                           &
                &   vectorID = VEC_ALIGN_RADIAL_ID )

            !! fchvisper: Perpendicular viscosity current
            call write_cell_vector_component(                                &
                &   vectorComponent = edge_profiles%ggd( time_sind )%        &
                &                     j_perpendicular_viscosity,             &
                &   b2CellData = fchvisper(:,:,0),                           &
                &   vectorID = VEC_ALIGN_POLOIDAL_ID )
            call write_cell_vector_component(                                &
                &   vectorComponent = edge_profiles%ggd( time_sind )%        &
                &                     j_perpendicular_viscosity,             &
                &   b2CellData = fchvisper(:,:,1),                           &
                &   vectorID = VEC_ALIGN_RADIAL_ID )

            !! fchvisq: Heat viscosity current
            call write_cell_vector_component(                                &
                &   vectorComponent = edge_profiles%ggd( time_sind )%        &
                &                     j_heat_viscosity,                      &
                &   b2CellData = fchvisq(:,:,0),                             &
                &   vectorID = VEC_ALIGN_POLOIDAL_ID )
            call write_cell_vector_component(                                &
                &   vectorComponent = edge_profiles%ggd( time_sind )%        &
                &                     j_heat_viscosity,                      &
                &   b2CellData = fchvisq(:,:,1),                             &
                &   vectorID = VEC_ALIGN_RADIAL_ID )

            !! fchdia: Diamagnetic current
            call write_cell_vector_component(                                &
                &   vectorComponent = edge_profiles%ggd( time_sind )%        &
                &                     j_diamagnetic,                         &
                &   b2CellData = fchdia(:,:,0),                              &
                &   vectorID = VEC_ALIGN_POLOIDAL_ID )
            call write_cell_vector_component(                                &
                &   vectorComponent = edge_profiles%ggd( time_sind )%        &
                &                     j_diamagnetic,                         &
                &   b2CellData = fchdia(:,:,1),                              &
                &   vectorID = VEC_ALIGN_RADIAL_ID )

            if (use_eirene.ne.0) then
#ifdef B25_EIRENE
                !! Neutral pressure
                do is = 1, nspecies
                   tmpCv(:,:) = 0.0_IDS_real
                   tmpVx(:,:) = 0.0_IDS_real
                   tmpFace(:,:,:) = 0.0_IDS_real
                   do iss = 1, natmi
                      if (latmscl(iss).eq.is) then
                        tmpCv(-1:nx,-1:ny) = tmpCv(-1:nx,-1:ny) + &
                           &  dab2(0:nx+1,0:ny+1,iss,1)*tab2(0:nx+1,0:ny+1,iss,1)
                        tmpVx(-1:nx,-1:ny) = tmpVx(-1:nx,-1:ny) + &
                           &  dab2(0:nx+1,0:ny+1,iss,1)
                        tmpFace(-1:nx,-1:ny,0) = tmpFace(-1:nx,-1:ny,0) + &
                           &  pfluxa(0:nx+1,0:ny+1,iss,1)
                        tmpFace(-1:nx,-1:ny,1) = tmpFace(-1:nx,-1:ny,1) + &
                           &  rfluxa(0:nx+1,0:ny+1,iss,1)
                      end if
                   end do
                   call write_quantity(                                      &
                       &   val = edge_profiles%ggd( time_sind )%             &
                       &         neutral( is )%pressure,                     &
                       &   value = tmpCv,                                    &
                       &   time_sind = time_sind )
                   do ix = -1, nx
                      do iy = -1, ny
                        if (tmpVx(ix,iy).gt.0.0_IDS_real) then
                          tmpCv(ix,iy) = tmpCv(ix,iy)/tmpVx(ix,iy)/qe
                        else
                          tmpCv(ix,iy) = 1.0e-6_IDS_real
                        end if
                      end do
                   end do
                !! Neutral temperature
                   call write_quantity(                                      &
                       &   val = edge_profiles%ggd( time_sind )%             &
                       &         neutral( is )%temperature,                  &
                       &   value = tmpCv,                                    &
                       &   time_sind = time_sind )
                !! Neutral density
                   call write_quantity(                                      &
                      &   val = edge_profiles%ggd( time_sind )%              &
                      &         neutral( is )%density,                       &
                      &   value = tmpVx,                                     &
                      &   time_sind = time_sind )
                !! Neutral particle flux
                   call write_face_scalar(                                   &
                       &   val = edge_transport%model(1)%ggd( time_sind )%   &
                       &         neutral( is )%particles%flux,               &
                       &   value = tmpFace,                                  &
                       &   time_sind = time_sind )
                !! Neutral particle energy flux
                   tmpFace(:,:,:) = 0.0_IDS_real
                   do iss = 1, natmi
                      if (latmscl(iss).eq.is) then
                        tmpFace(-1:nx,-1:ny,0) = tmpFace(-1:nx,-1:ny,0) + &
                           &  pefluxa(0:nx+1,0:ny+1,iss,1)
                        tmpFace(-1:nx,-1:ny,1) = tmpFace(-1:nx,-1:ny,1) + &
                           &  refluxa(0:nx+1,0:ny+1,iss,1)
                      end if
                   end do
                   call write_face_scalar(                                   &
                       &   val = edge_transport%model(1)%ggd( time_sind )%   &
                       &         neutral( is )%energy%flux,                  &
                       &   value = tmpFace,                                  &
                       &   time_sind = time_sind )
                   if (balance_netcdf.ne.0) then
                     tmpCv = 0.0_IDS_real
                     do istrai = 1, size(eirene_mc_paat_sna_bal,4)
                       do iss = 1, natmi
                         if (latmscl(iss).eq.is) then
                           tmpCv(:,:) = tmpCv(:,:) &
                               &   + eirene_mc_paat_sna_bal(:,:,iss,istrai)  &
                               &   + eirene_mc_pmat_sna_bal(:,:,iss,istrai)  &
                               &   + eirene_mc_piat_sna_bal(:,:,iss,istrai)
                         end if
                       end do
                     end do
                     tmpCv(:,:) = tmpCv(:,:) / vol(:,:)
                     call write_cell_scalar(                                 &
                       &   scalar = edge_sources%source(1)%ggd( time_sind )% &
                       &            neutral( is )%particles,                 &
                       &   b2CellData = tmpCv )
                     tmpCv = 0.0_IDS_real
                     do istrai = 1, size(eirene_mc_paat_sna_bal,4)
                       do iss = 1, natmi
                         if (latmscl(iss).eq.is) then
                            tmpCv(:,:) = tmpCv(:,:)                          &
                               &       + eirene_mc_paat_sna_bal(:,:,iss,istrai)
                         end if
                       end do
                     end do
                     tmpCv(:,:) = tmpCv(:,:) / vol(:,:)
                     call write_cell_scalar(                                 &
                       &   scalar = edge_sources%source(5)%ggd( time_sind )% &
                       &            neutral( is )%particles,                 &
                       &   b2CellData = tmpCv )
                     tmpCv = 0.0_IDS_real
                     do istrai = 1, size(eirene_mc_pmat_sna_bal,4)
                       do iss = 1, natmi
                         if (latmscl(iss).eq.is) then
                            tmpCv(:,:) = tmpCv(:,:) &
                               &   + eirene_mc_pmat_sna_bal(:,:,iss,istrai)  &
                               &   + eirene_mc_piat_sna_bal(:,:,iss,istrai)
                         end if
                       end do
                     end do
                     tmpCv(:,:) = tmpCv(:,:) / vol(:,:)
                     call write_cell_scalar(                                 &
                       &   scalar = edge_sources%source(6)%ggd( time_sind )% &
                       &            neutral( is )%particles,                 &
                       &   b2CellData = tmpCv )
                   end if
                   tmpCv = 0.0_IDS_real
                   do iss = 1, natmi
                      if (latmscl(iss).eq.is) then
                         tmpCv(-1:nx,-1:ny) = tmpCv(-1:nx,-1:ny) + &
                            &                 eneutrad(0:nx+1,0:ny+1,iss,0)
                      end if
                   end do
                   tmpCv(:,:) = tmpCv(:,:) / vol(:,:)
                   call write_cell_scalar(                                   &
                       &   scalar = edge_sources%source(12)%ggd( time_sind )%&
                       &            neutral( is )%particles,                 &
                       &   b2CellData = tmpCv )
                 end do
                do is = 1, natmi
                   js = latmscl(is)
                   ks = isstat(is)
                   tmpCv(-1:nx,-1:ny) = tab2(0:nx+1,0:ny+1,is,1)/qe
                   call write_quantity(                                      &
                       &   val = edge_profiles%ggd( time_sind )%             &
                       &         neutral( js )%state( ks )%temperature,      &
                       &   value = tmpCv,                                    &
                       &   time_sind = time_sind )
                   tmpCv(-1:nx,-1:ny) = dab2(0:nx+1,0:ny+1,is,1)
                   call write_quantity(                                      &
                       &   val = edge_profiles%ggd( time_sind )%             &
                       &         neutral( js )%state( ks )%density,          &
                       &   value = tmpCv,                                    &
                       &   time_sind = time_sind )
                   tmpFace(-1:nx,-1:ny,0) = pfluxa(0:nx+1,0:ny+1,is,1)
                   tmpFace(-1:nx,-1:ny,1) = rfluxa(0:nx+1,0:ny+1,is,1)
                   call write_face_scalar(                                   &
                       &   val = edge_transport%model(1)%ggd( time_sind )%   &
                       &         neutral( js )%state( ks )%particles%flux,   &
                       &   value = tmpFace,                                  &
                       &   time_sind = time_sind )
                   tmpFace(-1:nx,-1:ny,0) = pefluxa(0:nx+1,0:ny+1,is,1)
                   tmpFace(-1:nx,-1:ny,1) = refluxa(0:nx+1,0:ny+1,is,1)
                   call write_face_scalar(                                   &
                       &   val = edge_transport%model(1)%ggd( time_sind )%   &
                       &         neutral( js )%state( ks )%energy%flux,      &
                       &   value = tmpFace,                                  &
                       &   time_sind = time_sind )
                   tmpCv(-1:nx,-1:ny) = dab2(0:nx+1,0:ny+1,is,1)*            &
                       &                tab2(0:nx+1,0:ny+1,is,1)
                   call write_quantity(                                      &
                       &   val = edge_profiles%ggd( time_sind )%             &
                       &         neutral( js )%state( ks )%pressure,         &
                       &   value = tmpCv,                                    &
                       &   time_sind = time_sind )
                   tmpCv = 0.0_IDS_real
                   call write_cell_vector_component(                         &
                       &   vectorComponent = edge_profiles%ggd( time_sind )% &
                       &                     neutral( js )%state( ks )%      &
                       &                     velocity_diamagnetic,           &
                       &   b2CellData = tmpCv,                               &
                       &   vectorID = VEC_ALIGN_POLOIDAL_ID )
                   call write_cell_vector_component(                         &
                       &   vectorComponent = edge_profiles%ggd( time_sind )% &
                       &                     neutral( js )%state( ks )%      &
                       &                     velocity_diamagnetic,           &
                       &   b2CellData = tmpCv,                               &
                       &   vectorID = VEC_ALIGN_RADIAL_ID )
                   call write_cell_vector_component(                         &
                       &   vectorComponent = edge_profiles%ggd( time_sind )% &
                       &                     neutral( js )%state( ks )%      &
                       &                     velocity_exb,                   &
                       &   b2CellData = tmpCv,                               &
                       &   vectorID = VEC_ALIGN_POLOIDAL_ID )
                   call write_cell_vector_component(                         &
                       &   vectorComponent = edge_profiles%ggd( time_sind )% &
                       &                     neutral( js )%state( ks )%      &
                       &                     velocity_exb,                   &
                       &   b2CellData = tmpCv,                               &
                       &   vectorID = VEC_ALIGN_RADIAL_ID )
                   if (balance_netcdf.ne.0) then
                     tmpCv = 0.0_IDS_real
                     do istrai = 1, size(eirene_mc_paat_sna_bal,4)
                       tmpCv(:,:) = tmpCv(:,:)                               &
                          &       + eirene_mc_paat_sna_bal(:,:,is,istrai)    &
                          &       + eirene_mc_pmat_sna_bal(:,:,is,istrai)    &
                          &       + eirene_mc_piat_sna_bal(:,:,is,istrai)
                     end do
                     tmpCv(:,:) = tmpCv(:,:) / vol(:,:)
                     call write_cell_scalar(                                 &
                       &   scalar = edge_sources%source(1)%ggd( time_sind )% &
                       &            neutral( js )%state( ks )%particles,     &
                       &   b2CellData = tmpCv )
                     tmpCv = 0.0_IDS_real
                     do istrai = 1, size(eirene_mc_paat_sna_bal,4)
                       tmpCv(:,:) = tmpCv(:,:)                               &
                          &       + eirene_mc_paat_sna_bal(:,:,is,istrai)
                     end do
                     tmpCv(:,:) = tmpCv(:,:) / vol(:,:)
                     call write_cell_scalar(                                 &
                       &   scalar = edge_sources%source(5)%ggd( time_sind )% &
                       &            neutral( js )%state( ks )%particles,     &
                       &   b2CellData = tmpCv )
                     tmpCv = 0.0_IDS_real
                     do istrai = 1, size(eirene_mc_pmat_sna_bal,4)
                       tmpCv(:,:) = tmpCv(:,:)                               &
                          &       + eirene_mc_pmat_sna_bal(:,:,is,istrai)    &
                          &       + eirene_mc_piat_sna_bal(:,:,is,istrai)
                     end do
                     tmpCv(:,:) = tmpCv(:,:) / vol(:,:)
                     call write_cell_scalar(                                 &
                       &   scalar = edge_sources%source(6)%ggd( time_sind )% &
                       &            neutral( js )%state( ks )%particles,     &
                       &   b2CellData = tmpCv )
                   end if
                   tmpCv(-1:nx,-1:ny) = eneutrad(0:nx+1,0:ny+1,is,0)/vol(-1:nx,-1:ny)
                   call write_cell_scalar(                                   &
                       &   scalar = edge_sources%source(12)%                 &
                       &            ggd( time_sind )%neutral( js )%          &
                       &            state( ks )%particles,                   &
                       &   b2CellData = tmpCv )
                end do

                !! Molecular quantities
                do js = nspecies+1, nneut
                   tmpCv(:,:) = 0.0_IDS_real
                   tmpVx(:,:) = 0.0_IDS_real
                   tmpFace(:,:,:) = 0.0_IDS_real
                   do is = 1, nmoli
                      if (imneut(is).eq.js) then
                        tmpCv(-1:nx,-1:ny) = tmpCv(-1:nx,-1:ny) + &
                           &  dmb2(0:nx+1,0:ny+1,is,1)*tmb2(0:nx+1,0:ny+1,is,1)
                        tmpVx(-1:nx,-1:ny) = tmpVx(-1:nx,-1:ny) + &
                           &  dmb2(0:nx+1,0:ny+1,is,1)
                        tmpFace(-1:nx,-1:ny,0) = tmpFace(-1:nx,-1:ny,0) + &
                           &  pfluxm(0:nx+1,0:ny+1,is,1)
                        tmpFace(-1:nx,-1:ny,1) = tmpFace(-1:nx,-1:ny,1) + &
                           &  rfluxm(0:nx+1,0:ny+1,is,1)
                      end if
                   end do
                 !! Molecular pressure
                   call write_quantity(                                      &
                       &   val = edge_profiles%ggd( time_sind )%             &
                       &         neutral( js )%pressure,                     &
                       &   value = tmpCv,                                    &
                       &   time_sind = time_sind )
                   do ix = -1, nx
                      do iy = -1, ny
                         if (tmpVx(ix,iy).gt.0.0_IDS_real) then
                             tmpCv(ix,iy) = tmpCv(ix,iy)/tmpVx(ix,iy)/qe
                         else
                             tmpCv(ix,iy) = 1.0e-6_IDS_real
                         end if
                      end do
                   end do
                 !! Molecular temperature
                   call write_quantity(                                      &
                       &   val = edge_profiles%ggd( time_sind )%             &
                       &         neutral( js )%temperature,                  &
                       &   value = tmpCv,                                    &
                       &   time_sind = time_sind )
                 !! Molecular density
                   call write_quantity(                                      &
                       &   val = edge_profiles%ggd( time_sind )%             &
                       &         neutral( js )%density,                      &
                       &   value = tmpVx,                                    &
                       &   time_sind = time_sind )
                 !! Molecular particular fluxes
                   call write_face_scalar(                                   &
                       &   val = edge_transport%model(1)%ggd( time_sind )%   &
                       &         neutral( js )%particles%flux,               &
                       &   value = tmpFace,                                  &
                       &   time_sind = time_sind )
                 !! Neutral particle energy flux
                   tmpFace(:,:,:) = 0.0_IDS_real
                   do is = 1, nmoli
                      if (imneut(is).eq.js) then
                        tmpFace(-1:nx,-1:ny,0) = tmpFace(-1:nx,-1:ny,0) + &
                           &  pefluxm(0:nx+1,0:ny+1,is,1)
                        tmpFace(-1:nx,-1:ny,1) = tmpFace(-1:nx,-1:ny,1) + &
                           &  refluxm(0:nx+1,0:ny+1,is,1)
                      end if
                   end do
                   call write_face_scalar(                                   &
                       &   val = edge_transport%model(1)%ggd( time_sind )%   &
                       &         neutral( js )%energy%flux,                  &
                       &   value = tmpFace,                                  &
                       &   time_sind = time_sind )
                   if (balance_netcdf.ne.0) then
                     tmpCv = 0.0_IDS_real
                     do istrai = 1, size(eirene_mc_paml_sna_bal,4)
                       do is = 1, nmoli
                         if (imneut(is).eq.js) then
                             tmpCv(:,:) = tmpCv(:,:)                         &
                          &       + eirene_mc_paml_sna_bal(:,:,is,istrai)    &
                          &       + eirene_mc_pmml_sna_bal(:,:,is,istrai)    &
                          &       + eirene_mc_piml_sna_bal(:,:,is,istrai)
                         end if
                       end do
                     end do
                     tmpCv(:,:) = tmpCv(:,:) / vol(:,:)
                     call write_cell_scalar(                                 &
                       &   scalar = edge_sources%source(1)%ggd( time_sind )% &
                       &            neutral( js )%particles,                 &
                       &   b2CellData = tmpCv )
                     tmpCv = 0.0_IDS_real
                     do istrai = 1, size(eirene_mc_paml_sna_bal,4)
                       do is = 1, nmoli
                         if (imneut(is).eq.js) then
                             tmpCv(:,:) = tmpCv(:,:)                         &
                                &       + eirene_mc_paml_sna_bal(:,:,is,istrai)
                         end if
                       end do
                     end do
                     tmpCv(:,:) = tmpCv(:,:) / vol(:,:)
                     call write_cell_scalar(                                 &
                       &   scalar = edge_sources%source(5)%ggd( time_sind )% &
                       &            neutral( js )%particles,                 &
                       &   b2CellData = tmpCv )
                     tmpCv = 0.0_IDS_real
                     do istrai = 1, size(eirene_mc_pmml_sna_bal,4)
                       do is = 1, nmoli
                         if (imneut(is).eq.js) then
                             tmpCv(:,:) = tmpCv(:,:)                         &
                           &      + eirene_mc_pmml_sna_bal(:,:,is,istrai)    &
                           &      + eirene_mc_piml_sna_bal(:,:,is,istrai)
                         end if
                       end do
                     end do
                     tmpCv(:,:) = tmpCv(:,:) / vol(:,:)
                     call write_cell_scalar(                                 &
                       &   scalar = edge_sources%source(6)%ggd( time_sind )% &
                       &            neutral( js )%particles,                 &
                       &   b2CellData = tmpCv )
                   end if
                   tmpCv = 0.0_IDS_real
                   do is = 1, nmoli
                      if (imneut(is).eq.js) then
                          tmpCv(-1:nx,-1:ny) = tmpCv(-1:nx,-1:ny)            &
                             &               + emolrad(0:nx+1,0:ny+1,is,0)
                      end if
                   end do
                   tmpCv(:,:) = tmpCv(:,:) / vol(:,:)
                   call write_cell_scalar(                                   &
                       &   scalar = edge_sources%source(12)%                 &
                       &            ggd( time_sind )%neutral( js )%particles,&
                       &   b2CellData = tmpCv )
                end do

                js = nspecies
                do is = 1, nmoli
                   ks = isstat(natmi+is)
                   if (ks.eq.1) js = js + 1
                   tmpCv(-1:nx,-1:ny) = tmb2(0:nx+1,0:ny+1,is,1)/qe
                   call write_quantity(                                      &
                       &   val = edge_profiles%ggd( time_sind )%             &
                       &         neutral( js )%state( ks )%temperature,      &
                       &   value = tmpCv,                                    &
                       &   time_sind = time_sind )
                   tmpCv(-1:nx,-1:ny) = dmb2(0:nx+1,0:ny+1,is,1)
                   call write_quantity(                                      &
                       &   val = edge_profiles%ggd( time_sind )%             &
                       &         neutral( js )%state( ks )%density,          &
                       &   value = tmpCv,                                    &
                       &   time_sind = time_sind )
                   tmpFace(-1:nx,-1:ny,0) = pfluxm(0:nx+1,0:ny+1,is,1)
                   tmpFace(-1:nx,-1:ny,1) = rfluxm(0:nx+1,0:ny+1,is,1)
                   call write_face_scalar(                                   &
                       &   val = edge_transport%model(1)%ggd( time_sind )%   &
                       &         neutral( js )%state( ks )%particles%flux,   &
                       &   value = tmpFace,                                  &
                       &   time_sind = time_sind )
                   tmpFace(-1:nx,-1:ny,0) = pefluxm(0:nx+1,0:ny+1,is,1)
                   tmpFace(-1:nx,-1:ny,1) = refluxm(0:nx+1,0:ny+1,is,1)
                   call write_face_scalar(                                   &
                       &   val = edge_transport%model(1)%ggd( time_sind )%   &
                       &         neutral( js )%state( ks )%energy%flux,      &
                       &   value = tmpFace,                                  &
                       &   time_sind = time_sind )
                   tmpCv(-1:nx,-1:ny) = dmb2(0:nx+1,0:ny+1,is,1)*            &
                       &                tmb2(0:nx+1,0:ny+1,is,1)
                   call write_quantity(                                      &
                       &   val = edge_profiles%ggd( time_sind )%             &
                       &         neutral( js )%state( ks )%pressure,         &
                       &   value = tmpCv,                                    &
                       &   time_sind = time_sind )
                   tmpCv = 0.0_IDS_real
                   call write_cell_vector_component(                         &
                       &   vectorComponent = edge_profiles%ggd( time_sind )% &
                       &                     neutral( js )%state( ks )%      &
                       &                     velocity_diamagnetic,           &
                       &   b2CellData = tmpCv,                               &
                       &   vectorID = VEC_ALIGN_POLOIDAL_ID )
                   call write_cell_vector_component(                         &
                       &   vectorComponent = edge_profiles%ggd( time_sind )% &
                       &                     neutral( js )%state( ks )%      &
                       &                     velocity_diamagnetic,           &
                       &   b2CellData = tmpCv,                               &
                       &   vectorID = VEC_ALIGN_RADIAL_ID )
                   call write_cell_vector_component(                         &
                       &   vectorComponent = edge_profiles%ggd( time_sind )% &
                       &                     neutral( js )%state( ks )%      &
                       &                     velocity_exb,                   &
                       &   b2CellData = tmpCv,                               &
                       &   vectorID = VEC_ALIGN_POLOIDAL_ID )
                   call write_cell_vector_component(                         &
                       &   vectorComponent = edge_profiles%ggd( time_sind )% &
                       &                     neutral( js )%state( ks )%      &
                       &                     velocity_exb,                   &
                       &   b2CellData = tmpCv,                               &
                       &   vectorID = VEC_ALIGN_RADIAL_ID )
                   if (balance_netcdf.ne.0) then
                     tmpCv = 0.0_IDS_real
                     do istrai = 1, size(eirene_mc_paml_sna_bal,4)
                       tmpCv(:,:) = tmpCv(:,:)                               &
                          &       + eirene_mc_paml_sna_bal(:,:,is,istrai)    &
                          &       + eirene_mc_pmml_sna_bal(:,:,is,istrai)    &
                          &       + eirene_mc_piml_sna_bal(:,:,is,istrai)
                     end do
                     tmpCv(:,:) = tmpCv(:,:) / vol(:,:)
                     call write_cell_scalar(                                 &
                       &   scalar = edge_sources%source(1)%ggd( time_sind )% &
                       &            neutral( js )%state( ks )%particles,     &
                       &   b2CellData = tmpCv )
                     tmpCv = 0.0_IDS_real
                     do istrai = 1, size(eirene_mc_paml_sna_bal,4)
                       tmpCv(:,:) = tmpCv(:,:)                               &
                          &       + eirene_mc_paml_sna_bal(:,:,is,istrai)
                     end do
                     tmpCv(:,:) = tmpCv(:,:) / vol(:,:)
                     call write_cell_scalar(                                 &
                       &   scalar = edge_sources%source(5)%ggd( time_sind )% &
                       &            neutral( js )%state( ks )%particles,     &
                       &   b2CellData = tmpCv )
                     tmpCv = 0.0_IDS_real
                     do istrai = 1, size(eirene_mc_pmml_sna_bal,4)
                       tmpCv(:,:) = tmpCv(:,:)                               &
                           &      + eirene_mc_pmml_sna_bal(:,:,is,istrai)    &
                           &      + eirene_mc_piml_sna_bal(:,:,is,istrai)
                     end do
                     tmpCv(:,:) = tmpCv(:,:) / vol(:,:)
                     call write_cell_scalar(                                 &
                       &   scalar = edge_sources%source(6)%ggd( time_sind )% &
                       &            neutral( js )%state( ks )%particles,     &
                       &   b2CellData = tmpCv )
                   end if
                   tmpCv(-1:nx,-1:ny) = emolrad(0:nx+1,0:ny+1,is,0)/vol(-1:nx,-1:ny)
                   call write_cell_scalar(                                   &
                       &   scalar = edge_sources%source(12)%                 &
                       &            ggd( time_sind )%neutral( js )%          &
                       &            state( ks )%particles,                   &
                       &   b2CellData = tmpCv )
                end do
#endif
            else
                do is = 1, nspecies
                    js = eb2spcr(is)
                !! na : Fluid neutral density
                    call write_quantity(                                      &
                        &   val = edge_profiles%ggd( time_sind )%             &
                        &         neutral( is )%density,                      &
                        &   value = na(:,:,js),                               &
                        &   time_sind = time_sind )
                    call write_quantity(                                      &
                        &   val = edge_profiles%ggd( time_sind )%             &
                        &         neutral( is )%state(1)%density,             &
                        &   value = na(:,:,js),                               &
                        &   time_sind = time_sind )
                !! ua: Parallel fluid neutral velocity
                    call write_cell_vector_component(                         &
                        &   vectorComponent = edge_profiles%ggd( time_sind )% &
                        &                     neutral( is )%velocity,         &
                        &   b2CellData = ua(:,:,js),                          &
                        &   vectorID = VEC_ALIGN_PARALLEL_ID )
                    call write_cell_vector_component(                         &
                        &   vectorComponent = edge_profiles%ggd( time_sind )% &
                        &                     neutral( is )%state(1)%velocity,&
                        &   b2CellData = ua(:,:,js),                          &
                        &   vectorID = VEC_ALIGN_PARALLEL_ID )
                !! wadia: Diamagnetic fluid neutral velocity
                    call write_cell_vector_component(                         &
                        &   vectorComponent = edge_profiles%ggd( time_sind )% &
                        &                     neutral( is )%state(1)%         &
                        &                     velocity_diamagnetic,           &
                        &   b2CellData = wadia(:,:,0,js),                     &
                        &   vectorID = VEC_ALIGN_POLOIDAL_ID )
                    call write_cell_vector_component(                         &
                        &   vectorComponent = edge_profiles%ggd( time_sind )% &
                        &                     neutral( is )%state(1)%         &
                        &                     velocity_diamagnetic,           &
                        &   b2CellData = wadia(:,:,1,js),                     &
                        &   vectorID = VEC_ALIGN_RADIAL_ID )
                !! vaecrb: ExB fluid neutral velocity
                    call write_cell_vector_component(                         &
                        &   vectorComponent = edge_profiles%ggd( time_sind )% &
                        &                     neutral( is )%state(1)%         &
                        &                     velocity_exb,                   &
                        &   b2CellData = vaecrb(:,:,0,js),                    &
                        &   vectorID = VEC_ALIGN_POLOIDAL_ID )
                    call write_cell_vector_component(                         &
                        &   vectorComponent = edge_profiles%ggd( time_sind )% &
                        &                     neutral( is )%state(1)%         &
                        &                     velocity_exb,                   &
                        &   b2CellData = vaecrb(:,:,1,js),                    &
                        &   vectorID = VEC_ALIGN_RADIAL_ID )
                !! fna: Fluid neutral particle flux
                    do ix = -1, nx
                      do iy = -1, ny
                        if (fna(ix,iy,0,js).ne.0.0_R8) then
                          tmpFace(ix,iy,0) = fna(ix,iy,0,js)/gs(ix,iy,0)/qc(ix,iy)
                        else
                          tmpFace(ix,iy,0) = 0.0_R8
                        end if
                        if (fna(ix,iy,1,js).ne.0.0_R8) then
                          tmpFace(ix,iy,1) = fna(ix,iy,1,js)/gs(ix,iy,1)
                        else
                          tmpFace(ix,iy,1) = 0.0_R8
                        end if
                      end do
                    end do
                    call write_face_scalar(                                   &
                        &   val = edge_transport%model(1)%ggd( time_sind )%   &
                        &         neutral( is )%particles%flux,               &
                        &   value = tmpFace,                                  &
                        &   time_sind = time_sind )
                    call write_face_scalar(                                   &
                        &   val = edge_transport%model(1)%ggd( time_sind )%   &
                        &         neutral( is )%state(1)%particles%flux,      &
                        &   value = tmpFace,                                  &
                        &   time_sind = time_sind )
                !! pb : Fluid neutral pressure
                    call b2xppb( nx, ny, rza(:,:,js), na(:,:,js), te, ti, pb)
                    call write_quantity(                                      &
                        &   val = edge_profiles%ggd( time_sind )%             &
                        &         neutral( is )%pressure,                     &
                        &   value = pb,                                       &
                        &   time_sind = time_sind )
                    call write_quantity(                                      &
                        &   val = edge_profiles%ggd( time_sind )%             &
                        &         neutral( is )%state(1)%pressure,            &
                        &   value = pb,                                       &
                        &   time_sind = time_sind )
                !! cdpa: Fluid neutral diffusivity
                    call write_face_scalar(                                   &
                        &   val = edge_transport%model(1)%ggd( time_sind )%   &
                        &         neutral( is )%particles%d,                  &
                        &   value = cdpa(:,:,:,js),                           &
                        &   time_sind = time_sind )
                    call write_face_scalar(                                   &
                        &   val = edge_transport%model(1)%ggd( time_sind )%   &
                        &         neutral( is )%state(1)%particles%d,         &
                        &   value = cdpa(:,:,:,js),                           &
                        &   time_sind = time_sind )
                !! Fluid neutral kinetic energy density
                    tmpCv(:,:) = 0.5_IDS_real*am(js)*mp*(ua(:,:,js)**2)*na(:,:,js)
                    call write_quantity(                                      &
                        &   val = edge_profiles%ggd( time_sind )%             &
                        &         neutral( is )%energy_density_kinetic,       &
                        &   value = tmpCv,                                    &
                        &   time_sind = time_sind )
                    call write_quantity(                                      &
                        &   val = edge_profiles%ggd( time_sind )%             &
                        &         neutral( is )%state(1)%                     &
                        &         energy_density_kinetic,                     &
                        &   value = tmpCv,                                    &
                        &   time_sind = time_sind )
                !! cvsa: Ion diffusivity
                    call write_cell_vector_component(                         &
                        &   vectorComponent = edge_transport%model(1)%        &
                        &                     ggd( time_sind )%neutral( is )% &
                        &                     momentum%d,                     &
                        &   b2CellData = cvsa(:,:,0,js),                      &
                        &   vectorID = VEC_ALIGN_POLOIDAL_ID )
                    call write_cell_vector_component(                         &
                        &   vectorComponent = edge_transport%model(1)%        &
                        &                     ggd( time_sind )%neutral( is )% &
                        &                     momentum%d,                     &
                        &   b2CellData = cvsa(:,:,1,js),                      &
                        &   vectorID = VEC_ALIGN_RADIAL_ID )
                    call write_cell_vector_component(                         &
                        &   vectorComponent = edge_transport%model(1)%        &
                        &                     ggd( time_sind )%neutral( is )% &
                        &                     state(1)%momentum%d,            &
                        &   b2CellData = cvsa(:,:,0,js),                      &
                        &   vectorID = VEC_ALIGN_POLOIDAL_ID )
                    call write_cell_vector_component(                         &
                        &   vectorComponent = edge_transport%model(1)%        &
                        &                     ggd( time_sind )%neutral( is )% &
                        &                     state(1)%momentum%d,            &
                        &   b2CellData = cvsa(:,:,1,js),                      &
                        &   vectorID = VEC_ALIGN_RADIAL_ID )
                !! fmo: Ion momentum flux
                    do ix = -1, nx
                      do iy = -1, ny
                        if (fmo(ix,iy,0,js).ne.0.0_R8) then
                          tmpFace(ix,iy,0) = fmo(ix,iy,0,js)/gs(ix,iy,0)/qc(ix,iy)
                        else
                          tmpFace(ix,iy,0) = 0.0_R8
                        end if
                        if (fmo(ix,iy,1,js).ne.0.0_R8) then
                          tmpFace(ix,iy,1) = fmo(ix,iy,1,js)/gs(ix,iy,1)
                        else
                          tmpFace(ix,iy,1) = 0.0_R8
                        end if
                      end do
                    end do
                    call write_cell_vector_component(                         &
                        &   vectorComponent = edge_transport%model(1)%        &
                        &                     ggd( time_sind )%neutral( is )% &
                        &                     momentum%flux,                  &
                        &   b2CellData = tmpFace(:,:,0),                      &
                        &   vectorID = VEC_ALIGN_PARALLEL_ID )
                    call write_cell_vector_component(                         &
                        &   vectorComponent = edge_transport%model(1)%        &
                        &                     ggd( time_sind )%neutral( is )% &
                        &                     state(1)%momentum%flux,         &
                        &   b2CellData = tmpFace(:,:,0),                      &
                        &   vectorID = VEC_ALIGN_PARALLEL_ID )
                !! fllim0fna: Fluid neutral flux limiter
                    call write_face_scalar(                                   &
                        &   val = edge_transport%model(1)%ggd( time_sind )%   &
                        &         neutral( is )%particles%flux_limiter,       &
                        &   value = fllim0fna(:,:,:,js),                      &
                        &   time_sind = time_sind )
                    call write_face_scalar(                                   &
                        &   val = edge_transport%model(1)%ggd( time_sind )%   &
                        &         neutral( is )%state(1)%                     &
                        &         particles%flux_limiter,                     &
                        &   value = fllim0fna(:,:,:,js),                      &
                        &   time_sind = time_sind )
                !! fllimvisc: Fluid neutral momentum transport flux limit
                    call write_cell_vector_component(                         &
                        &   vectorComponent = edge_transport%model(1)%        &
                        &                     ggd( time_sind )%neutral( is )% &
                        &                     momentum%flux_limiter,          &
                        &   b2CellData = fllimvisc(:,:,js),                   &
                        &   vectorID = VEC_ALIGN_PARALLEL_ID )
                    call write_cell_vector_component(                         &
                        &   vectorComponent = edge_transport%model(1)%        &
                        &                     ggd( time_sind )%neutral( is )% &
                        &                     state(1)%momentum%flux_limiter, &
                        &   b2CellData = fllimvisc(:,:,js),                   &
                        &   vectorID = VEC_ALIGN_PARALLEL_ID )
                !! sna: Fluid neutral particle sources
                    tmpCv(:,:) = ( sna(:,:,0,js) +                            &
                        &          sna(:,:,1,js) * na(:,:,js) ) / vol(:,:)
                    call write_cell_scalar( scalar = edge_sources%source(1)%  &
                        &   ggd( time_sind )%neutral( is )%particles,         &
                        &   b2CellData = tmpCv )
                    call write_cell_scalar( scalar = edge_sources%source(1)%  &
                        &   ggd( time_sind )%neutral( is )%state(1)%particles,&
                        &   b2CellData = tmpCv )
                    tmpCv(:,:) = ext_sna(:,:,js) / vol(:,:)
                    call write_cell_scalar(                                   &
                        &   scalar = edge_sources%source(2)%ggd( time_sind )% &
                        &            neutral( is )%particles,                 &
                        &   b2CellData = tmpCv )
                    call write_cell_scalar(                                   &
                        &   scalar = edge_sources%source(2)%ggd( time_sind )% &
                        &            neutral( is )%state(1)%particles,        &
                        &   b2CellData = tmpCv )
                    tmpCv(:,:) = ( b2stbc_sna(:,:,js) +                       &
                        &          b2stbm_sna(:,:,js) ) / vol(:,:)
                    call write_cell_scalar(                                   &
                        &   scalar = edge_sources%source(3)%ggd( time_sind )% &
                        &            neutral( is )%particles,                 &
                        &   b2CellData = tmpCv )
                    call write_cell_scalar(                                   &
                        &   scalar = edge_sources%source(3)%ggd( time_sind )% &
                        &            neutral( is )%state(1)%particles,        &
                        &   b2CellData = tmpCv )
                    tmpCv(:,:) = ( snadt(:,:,0,js) +                          &
                        &          snadt(:,:,1,js) * na(:,:,js) ) / vol(:,:)
                    call write_cell_scalar(                                   &
                        &   scalar = edge_sources%source(4)%ggd( time_sind )% &
                        &            neutral( is )%particles,                 &
                        &   b2CellData = tmpCv )
                    call write_cell_scalar(                                   &
                        &   scalar = edge_sources%source(4)%ggd( time_sind )% &
                        &            neutral( is )%state(1)%particles,        &
                        &   b2CellData = tmpCv )
                    tmpCv(:,:) = b2stbr_sna(:,:,js) / vol(:,:)
                    call write_cell_scalar(                                   &
                        &   scalar = edge_sources%source(5)%ggd( time_sind )% &
                        &            neutral( is )%particles,                 &
                        &   b2CellData = tmpCv )
                    call write_cell_scalar(                                   &
                        &   scalar = edge_sources%source(5)%ggd( time_sind )% &
                        &            neutral( is )%state(1)%particles,        &
                        &   b2CellData = tmpCv )
                    tmpCv(:,:) = rsana(:,:,js) / vol(:,:)
                    call write_cell_scalar(                                   &
                        &   scalar = edge_sources%source(7)%ggd( time_sind )% &
                        &            neutral( is )%particles,                 &
                        &   b2CellData = tmpCv )
                    call write_cell_scalar(                                   &
                        &   scalar = edge_sources%source(7)%ggd( time_sind )% &
                        &            neutral( is )%state(1)%particles,        &
                        &   b2CellData = tmpCv )
                    tmpCv(:,:) = rrana(:,:,js) / vol(:,:)
                    call write_cell_scalar(                                   &
                        &   scalar = edge_sources%source(8)%ggd( time_sind )% &
                        &            neutral( is )%particles,                 &
                        &   b2CellData = tmpCv )
                    call write_cell_scalar(                                   &
                        &   scalar = edge_sources%source(8)%ggd( time_sind )% &
                        &            neutral( is )%state(1)%particles,        &
                        &   b2CellData = tmpCv )
                    tmpCv(:,:) = rcxna(:,:,js) / vol(:,:)
                    call write_cell_scalar(                                   &
                        &   scalar = edge_sources%source(9)%ggd( time_sind )% &
                        &            neutral( is )%particles,                 &
                        &   b2CellData = tmpCv )
                    call write_cell_scalar(                                   &
                        &   scalar = edge_sources%source(9)%ggd( time_sind )% &
                        &            neutral( is )%state(1)%particles,        &
                        &   b2CellData = tmpCv )
                !! smo: Ion parallel momentum sources
                    tmpCv(:,:) = ( smo(:,:,0,js) +                            &
                        &          smo(:,:,1,js) * ua(:,:,js) +               &
                        &          smo(:,:,2,js) * roxa(:,:,js) +             &
                        &          smo(:,:,3,js) * roxa(:,:,js)               &
                        &                        * ua(:,:,js) ) / vol(:,:)
                    call write_cell_vector_component(                         &
                        &   vectorComponent = edge_sources%source(1)%         &
                        &                     ggd( time_sind )%               &
                        &                     neutral( is )%momentum,         &
                        &   b2CellData = tmpCv,                               &
                        &   vectorID = VEC_ALIGN_PARALLEL_ID )
                    call write_cell_vector_component(                         &
                        &   vectorComponent = edge_sources%source(1)%         &
                        &                     ggd( time_sind )%neutral( is )% &
                        &                     state(1)%momentum,              &
                        &   b2CellData = tmpCv,                               &
                        &   vectorID = VEC_ALIGN_PARALLEL_ID )
                    tmpCv(:,:) = ext_smo(:,:,js) / vol(:,:)
                    call write_cell_vector_component(                         &
                        &   vectorComponent = edge_sources%source(2)%         &
                        &                     ggd( time_sind )%neutral( is )% &
                        &                     momentum,                       &
                        &   b2CellData = tmpCv,                               &
                        &   vectorID = VEC_ALIGN_PARALLEL_ID )
                    call write_cell_vector_component(                         &
                        &   vectorComponent = edge_sources%source(2)%         &
                        &                     ggd( time_sind )%neutral( is )% &
                        &                     state(1)%momentum,              &
                        &   b2CellData = tmpCv,                               &
                        &   vectorID = VEC_ALIGN_PARALLEL_ID )
                    tmpCv(:,:) = ( b2stbc_smo(:,:,js) +                       &
                        &          b2stbm_smo(:,:,js) ) / vol(:,:)
                    call write_cell_vector_component(                         &
                        &   vectorComponent = edge_sources%source(3)%         &
                        &                     ggd( time_sind )%neutral( is )% &
                        &                     momentum,                       &
                        &   b2CellData = tmpCv,                               &
                        &   vectorID = VEC_ALIGN_PARALLEL_ID )
                    call write_cell_vector_component(                         &
                        &   vectorComponent = edge_sources%source(3)%         &
                        &                     ggd( time_sind )%neutral( is )% &
                        &                     state(1)%momentum,              &
                        &   b2CellData = tmpCv,                               &
                        &   vectorID = VEC_ALIGN_PARALLEL_ID )
                    tmpCv(:,:) = ( smodt(:,:,0,js) +                          &
                        &          smodt(:,:,1,js) * ua(:,:,js) +             &
                        &          smodt(:,:,2,js) * roxa(:,:,js) +           &
                        &          smodt(:,:,3,js) * roxa(:,:,js)             &
                        &                          * ua(:,:,js) ) / vol(:,:)
                    call write_cell_vector_component(                         &
                        &   vectorComponent = edge_sources%source(4)%         &
                        &                     ggd( time_sind )%neutral( is )% &
                        &                     momentum,                       &
                        &   b2CellData = tmpCv,                               &
                        &   vectorID = VEC_ALIGN_PARALLEL_ID )
                    call write_cell_vector_component(                         &
                        &   vectorComponent = edge_sources%source(4)%         &
                        &                     ggd( time_sind )%neutral( is )% &
                        &                     state(1)%momentum,              &
                        &   b2CellData = tmpCv,                               &
                        &   vectorID = VEC_ALIGN_PARALLEL_ID )
                    tmpCv(:,:) = b2stbr_smo(:,:,js) / vol(:,:)
                    call write_cell_vector_component(                         &
                        &   vectorComponent = edge_sources%source(5)%         &
                        &            ggd( time_sind )%neutral( is )%momentum, &
                        &   b2CellData = tmpCv,                               &
                        &   vectorID = VEC_ALIGN_PARALLEL_ID )
                    call write_cell_vector_component(                         &
                        &   vectorComponent = edge_sources%source(5)%         &
                        &            ggd( time_sind )%neutral( is )%          &
                        &            state(1)%momentum,                       &
                        &   b2CellData = tmpCv,                               &
                        &   vectorID = VEC_ALIGN_PARALLEL_ID )
                    tmpCv(:,:) = rsamo(:,:,js) / vol(:,:)
                    call write_cell_vector_component(                         &
                        &   vectorComponent = edge_sources%source(7)%         &
                        &                     ggd( time_sind )%neutral( is )% &
                        &                     momentum,                       &
                        &   b2CellData = tmpCv,                               &
                        &   vectorID = VEC_ALIGN_PARALLEL_ID )
                    call write_cell_vector_component(                         &
                        &   vectorComponent = edge_sources%source(7)%         &
                        &                     ggd( time_sind )%neutral( is )% &
                        &                     state(1)%momentum,              &
                        &   b2CellData = tmpCv,                               &
                        &   vectorID = VEC_ALIGN_PARALLEL_ID )
                    tmpCv(:,:) = rramo(:,:,js) / vol(:,:)
                    call write_cell_vector_component(                         &
                        &   vectorComponent = edge_sources%source(8)%         &
                        &                     ggd( time_sind )%neutral( is )% &
                        &                     momentum,                       &
                        &   b2CellData = tmpCv,                               &
                        &   vectorID = VEC_ALIGN_PARALLEL_ID )
                    call write_cell_vector_component(                         &
                        &   vectorComponent = edge_sources%source(8)%         &
                        &                     ggd( time_sind )%neutral( is )% &
                        &                     state(1)%momentum,              &
                        &   b2CellData = tmpCv,                               &
                        &   vectorID = VEC_ALIGN_PARALLEL_ID )
                    tmpCv(:,:) = rcxmo(:,:,js) / vol(:,:)
                    call write_cell_vector_component(                         &
                        &   vectorComponent = edge_sources%source(9)%         &
                        &                     ggd( time_sind )%neutral( is )% &
                        &                     momentum,                       &
                        &   b2CellData = tmpCv,                               &
                        &   vectorID = VEC_ALIGN_PARALLEL_ID )
                    call write_cell_vector_component(                         &
                        &   vectorComponent = edge_sources%source(9)%         &
                        &                     ggd( time_sind )%neutral( is )% &
                        &                     state(1)%momentum,              &
                        &   b2CellData = tmpCv,                               &
                        &   vectorID = VEC_ALIGN_PARALLEL_ID )
                end do
            end if

            !! po: Electric potential
            call write_quantity(                                        &
                &   val = edge_profiles%ggd( time_sind )%phi_potential, &
                &   value = po,                                         &
                &   time_sind = time_sind )
            !! Current sources
            tmpCv(:,:) = ( sch(:,:,0) + sch(:,:,1) * po(:,:) +          &
                &          sch(:,:,2) * ne(:,:) +                       &
                &          sch(:,:,3) * ne(:,:) * po(:,:) ) / vol(:,:)
            call write_cell_scalar(                                     &
                &   scalar = edge_sources%source(1)%ggd( time_sind )%   &
                &            current,                                   &
                &   b2CellData = tmpCv )
            tmpCv(:,:) = ext_sch(:,:) / vol(:,:)
            call write_cell_scalar(                                     &
                &   scalar = edge_sources%source(2)%ggd( time_sind )%   &
                &            current,                                   &
                &   b2CellData = tmpCv )
            tmpCv(:,:) = ( b2stbc_sch(:,:) + b2stbm_sch(:,:) ) / vol(:,:)
            call write_cell_scalar(                                     &
                &   scalar = edge_sources%source(3)%ggd( time_sind )%   &
                &            current,                                   &
                &   b2CellData = tmpCv )
            tmpCv(:,:) = ( schdt(:,:,0) + schdt(:,:,1) * po(:,:) +      &
                &          schdt(:,:,2) * ne(:,:) +                     &
                &          schdt(:,:,3) * ne(:,:) * po(:,:) ) / vol(:,:)
            call write_cell_scalar(                                     &
                &   scalar = edge_sources%source(4)%ggd( time_sind )%   &
                &            current,                                   &
                &   b2CellData = tmpCv )
            tmpCv(:,:) = b2stbr_sch(:,:) / vol(:,:)
            call write_cell_scalar(                                     &
                &   scalar = edge_sources%source(5)%ggd( time_sind )%   &
                &            current,                                   &
                &   b2CellData = tmpCv )
            !! csig : Electric conductivity
            call write_cell_vector_component(                           &
                &   vectorComponent = edge_transport%model(1)%          &
                &         ggd( time_sind )%conductivity,                &
                &   b2CellData = csig(:,:,0),                           &
                &   vectorID = VEC_ALIGN_POLOIDAL_ID )
            call write_cell_vector_component(                           &
                &   vectorComponent = edge_transport%model(1)%          &
                &         ggd( time_sind )%conductivity,                &
                &   b2CellData = csig(:,:,1),                           &
                &   vectorID = VEC_ALIGN_RADIAL_ID )

            !! B (magnetic field vector)
            !! Compute unit basis vectors along the field directions
            call compute_Coordinate_Unit_Vectors(crx, cry, e(:,:,:,1),  &
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

            !! write the magnetic field vector in the B2 coordinate system
            call write_cell_vector_component(                                 &
                &   vectorComponent = edge_profiles%ggd( time_sind )%e_field, &
                &   b2CellData = bb(:,:,0:2),                                 &
                &   vectorID = "diamagnetic" )

#if IMAS_MINOR_VERSION > 21
            !! write the emissivity data
            !! Process 1. Line and recombination radiation from B2.5 ions
            do is = 0, ns-1
                if (is_neutral(is) .and. use_eirene.eq.0) then
                   js = b2espcr(is)
                   tmpCv(:,:) = rqrad(:,:,is) / vol(:,:)
                   call write_cell_scalar( scalar = radiation%process(1)%     &
                       &   ggd( time_sind )%neutral( js )%emissivity,         &
                       &   b2CellData = tmpCv )
#if IMAS_MINOR_VERSION > 24
                   call write_cell_scalar( scalar = radiation%process(1)%     &
                       &   ggd( time_sind )%neutral( js )%                    &
                       &   state(1)%emissivity,                               &
                       &   b2CellData = tmpCv )
#endif
                else if (.not.is_neutral(is)) then
                   js = fluids_list(is)
                   tmpCv(:,:) = rqrad(:,:,is) / vol(:,:)
                   call write_cell_scalar( scalar = radiation%process(1)%     &
                       &   ggd( time_sind )%ion( js )%emissivity,             &
                       &   b2CellData = tmpCv )
#if IMAS_MINOR_VERSION > 24
                   call write_cell_scalar( scalar = radiation%process(1)%     &
                       &   ggd( time_sind )%ion( js )%state(1)%emissivity,    &
                       &   b2CellData = tmpCv )
#endif
                end if
            end do
            !! Process 2. Bremsstrahlung from B2.5 ions
            do is = 0, ns-1
                js = fluids_list(is)
                if (js.eq.0) cycle
                tmpCv(:,:) = rqbrm(:,:,is) / vol(:,:)
                call write_cell_scalar( scalar = radiation%process(2)%        &
                    &   ggd( time_sind )%ion( js )%emissivity,                &
                    &   b2CellData = tmpCv )
            end do
#ifdef B25_EIRENE
            if (use_eirene.ne.0) then

              !! Process 3. Eirene neutrals (atoms and molecules)
              do js = 1, nspecies
                tmpCv = 0.0_IDS_real
                do is = 1, natmi
                  if (latmscl(is).eq.js) then
                    do ix = -1, nx
                      do iy = -1, ny
                        tmpCv(ix,iy) = tmpCv(ix,iy)-eneutrad(ix+1,iy+1,is,0)
                      end do
                    end do
                  end if
                end do
                tmpCv(:,:) = tmpCv(:,:) / vol(:,:)
                call write_cell_scalar( scalar = radiation%process(3)%        &
                    &   ggd( time_sind )%neutral( js )%emissivity,            &
                    &   b2CellData = tmpCv )
              end do
              do is = 1, natmi
                js = latmscl(is)
                ks = isstat(is)
                tmpCv(-1:nx,-1:ny)=-eneutrad(0:nx+1,0:ny+1,is,0) / vol(-1:nx,-1:ny)
#if IMAS_MINOR_VERSION > 24
                call write_cell_scalar( scalar = radiation%process(3)%         &
                    &   ggd( time_sind )%neutral( js )%state( ks )%emissivity, &
                    &   b2CellData = tmpCv )
#endif
              end do
              do js = nspecies+1, nneut
                 tmpCv(:,:) = 0.0_IDS_real
                 do is = 1, nmoli
                    if (imneut(is).eq.js) then
                      do ix = -1, nx
                        do iy = -1, ny
                           tmpCv(ix,iy) = tmpCV(ix,iy)-emolrad(ix+1,iy+1,is,0)
                        end do
                      end do
                    end if
                 end do
                 tmpCv(:,:) = tmpCv(:,:) / vol(:,:)
                 call write_cell_scalar( scalar = radiation%process(3)%       &
                    &   ggd( time_sind )%neutral( js )%emissivity,            &
                    &   b2CellData = tmpCv )
              end do
              js = nspecies
              do is = 1, nmoli
                ks = isstat(natmi+is)
                if (ks.eq.1) js = js + 1
                do ix = -1, nx
                  do iy = -1, ny
                    tmpCv(ix,iy) = -emolrad(ix+1,iy+1,is,0) / vol(ix,iy)
                  end do
                end do
#if IMAS_MINOR_VERSION > 24
                call write_cell_scalar( scalar = radiation%process(3)%        &
                    &   ggd( time_sind )%neutral( js )%state( ks )%emissivity,&
                    &   b2CellData = tmpCv )
#endif
              end do
              !! Process 4. Eirene molecular ions
              do is = 1, nioni
                js = fluids_list(ns-1) + is
                do ix = -1, nx
                  do iy = -1, ny
                    tmpCv(ix,iy) = -eionrad(ix+1,iy+1,is,0) / vol(ix,iy)
                  end do
                end do
                call write_cell_scalar( scalar = radiation%process(4)%        &
                    &   ggd( time_sind )%ion( js )%emissivity,                &
                    &   b2CellData = tmpCv )
#if IMAS_MINOR_VERSION > 24
                call write_cell_scalar( scalar = radiation%process(4)%        &
                    &   ggd( time_sind )%ion( js )%state(1)%emissivity,       &
                    &   b2CellData = tmpCv )
#endif
              end do
            end if
#endif
#endif
#else
            call logmsg( LOGINFO, &
            &   "b2mod_ual_io.B25_process_ids: GGD not available, no plasma state writing" )
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
        allocate( summary%boundary%strike_point_inner_z%value( time_sind ) )
        allocate( summary%boundary%strike_point_outer_r%value( time_sind ) )
        allocate( summary%boundary%strike_point_outer_z%value( time_sind ) )
        if (LSN) then
          summary%boundary%strike_point_inner_r%value( time_sind ) = crx(-1,topcut(1),1)
          summary%boundary%strike_point_inner_z%value( time_sind ) = cry(-1,topcut(1),1)
          summary%boundary%strike_point_outer_r%value( time_sind ) = crx(nx,topcut(1),0)
          summary%boundary%strike_point_outer_z%value( time_sind ) = cry(nx,topcut(1),0)
        else
          summary%boundary%strike_point_inner_r%value( time_sind ) = crx(nx,topcut(1),0)
          summary%boundary%strike_point_inner_z%value( time_sind ) = cry(nx,topcut(1),0)
          summary%boundary%strike_point_outer_r%value( time_sind ) = crx(-1,topcut(1),1)
          summary%boundary%strike_point_outer_z%value( time_sind ) = cry(-1,topcut(1),1)
        endif
        allocate( summary%boundary%strike_point_outer_z%source(1) )
        summary%boundary%strike_point_outer_z%source = source
        allocate( summary%boundary%strike_point_inner_r%source(1) )
        summary%boundary%strike_point_inner_r%source = source
        allocate( summary%boundary%strike_point_inner_z%source(1) )
        summary%boundary%strike_point_inner_z%source = source
        allocate( summary%boundary%strike_point_outer_r%source(1) )
        summary%boundary%strike_point_outer_r%source = source

        allocate( summary%fusion%power%value( time_sind ) )
        summary%fusion%power%value( time_sind ) = fusion_power
        allocate( summary%fusion%power%source(1) )
        summary%fusion%power%source = source

#ifdef B25_EIRENE
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

        ib = 0
        do icnt = 1, nbc
          if (ib.ne.0) cycle
          if (bcchar(icnt).eq.'N') then
            if (bcstart(icnt).ne.-2 .and. &
             &  bcstart(icnt).le.jxa .and. jxa.le.bcend(icnt)) then
              ib = icnt
            else if (bcstart(icnt).eq.-2) then
              do ix = 1, bc_list_size(icnt)
                if (bc_list_x(ix,icnt).eq.jxa) ib = icnt
              end do
            end if
          end if
        end do
        if (ib.ne.0) then
          if (bcene(ib).eq.9.or.bcene(ib).eq.19) then
            allocate( summary%scrape_off_layer%t_e_decay_length%value( time_sind ) )
            summary%scrape_off_layer%t_e_decay_length%value( time_sind ) = enepar(ib,1)
            allocate( summary%scrape_off_layer%t_e_decay_length%source(1) )
            summary%scrape_off_layer%t_e_decay_length%source = source
          end if
          if (bceni(ib).eq.9.or.bceni(ib).eq.19) then
            allocate( summary%scrape_off_layer%t_i_average_decay_length%value( time_sind ) )
            summary%scrape_off_layer%t_i_average_decay_length%value( time_sind ) = enipar(ib,1)
            allocate( summary%scrape_off_layer%t_i_average_decay_length%source(1) )
            summary%scrape_off_layer%t_i_average_decay_length%source = source
          end if
          nibnd = IDS_REAL_INVALID
          match_found = .true.
          do is = 0, ns-1
            if (is_neutral(is).and.use_eirene.ne.0) cycle
            if (bccon(ib,is).eq.9) then
              if (nibnd.eq.IDS_REAL_INVALID) then
                nibnd = conpar(is,ib,1)
              else
                match_found = match_found.and.nibnd.eq.conpar(is,ib,1)
              end if
            end if
          end do
          if (match_found.and.nibnd.ne.IDS_REAL_INVALID) then
            allocate( summary%scrape_off_layer%n_e_decay_length%value( time_sind ) )
            summary%scrape_off_layer%n_e_decay_length%value( time_sind ) = nibnd
            allocate( summary%scrape_off_layer%n_e_decay_length%source(1) )
            summary%scrape_off_layer%n_e_decay_length%source = source
            allocate( summary%scrape_off_layer%n_i_total_decay_length%value( time_sind ) )
            summary%scrape_off_layer%n_i_total_decay_length%value( time_sind ) = nibnd
            allocate( summary%scrape_off_layer%n_i_total_decay_length%source(1) )
            summary%scrape_off_layer%n_i_total_decay_length%source = source
          end if
        end if
        u = 0.0_IDS_real
        do ix = 0, nx-1
          do iy = 0, ny-1
            if (geometryType.eq.GEOMETRY_LIMITER .or. &
             &  geometryType.eq.GEOMETRY_SN) then
              if (region(ix,iy,0).ne.2) cycle
            else if (geometryType.eq.GEOMETRY_STELLARATORISLAND) then
              if (region(ix,iy,0).ne.2 .and. region(ix,iy,0).ne.5) cycle
            else if (geometryType.eq.GEOMETRY_CDN .or. &
                  &  geometryType.eq.GEOMETRY_DDN_BOTTOM .or. &
                  &  geometryType.eq.GEOMETRY_DDN_TOP) then
              if (region(ix,iy,0).ne.2 .and. region(ix,iy,0).ne.6) cycle
            else
              cycle
            end if
            do is = 0, ns-1
              u = u + rqrad(ix,iy,is) + rqbrm(ix,iy,is)
            end do
#ifdef B25_EIRENE
            do is = 1, natmi
              u = u - eneutrad(ix+1,iy+1,is,0)
            end do
            do is = 1, nmoli
              u = u - emolrad(ix+1,iy+1,is,0)
            end do
            do is = 1, nioni
              u = u - eionrad(ix+1,iy+1,is,0)
            end do
#endif
          end do
        end do
        if (u.ne.0.0_IDS_real) then
          allocate( summary%scrape_off_layer%power_radiated%value( time_sind ) )
          summary%scrape_off_layer%power_radiated%value( time_sind ) = u
          allocate( summary%scrape_off_layer%power_radiated%source(1) )
          summary%scrape_off_layer%power_radiated%source(1) = source
        end if
        select case (geometryType)
        case (GEOMETRY_LIMITER)
          if (LSN) then
            ix = nx-1
          else
            ix = 0
          end if
        case (GEOMETRY_SN, GEOMETRY_CDN, GEOMETRY_DDN_BOTTOM)
          if (LSN) then
            ix = rightcut(1)-1
          else
            ix = leftcut(1)
          end if
        case (GEOMETRY_STELLARATORISLAND)
          ix = jxa
        case (GEOMETRY_DDN_TOP)
          ix = rightcut(2)
        case default
          ix = -2
        end select
        if (ix.ne.-2) then
          ii = ix
          jj = ix
          icnt = -1
          if (.not.LSN) icnt = 1
          allocate ( summary%scrape_off_layer%heat_flux_e_decay_length%value( time_sind ) )
          allocate ( summary%scrape_off_layer%heat_flux_i_decay_length%value( time_sind ) )
          qemax = 0.0_IDS_real
          qimax = 0.0_IDS_real
          do iy = jsep+1, ny-1
            qemax = qemax + u*fhe(ix,iy,0)
            qimax = qimax + u*fhi(ix,iy,0)
          end do
          do i = ix+icnt, jxa, icnt
            qetot = 0.0_IDS_real
            qitot = 0.0_IDS_real
            do iy = jsep+1, ny-1
              qetot = qetot + u*fhe(i,iy,0)
              qitot = qitot + u*fhi(i,iy,0)
            end do
            if (qetot*qemax.ge.0.0_IDS_real) then
              if (abs(qetot).gt.abs(qemax)) then
                ii = i
                qemax = qetot
              end if
            end if
            if (qitot*qimax.ge.0.0_IDS_real) then
              if (abs(qitot).gt.abs(qimax)) then
                jj = i
                qimax = qitot
              end if
            end if
          end do
          lambda = 0.0_IDS_real
          tmpFace(:,:,0) = abs(fhe(:,:,0))/gs(:,:,0)/qc(:,:)
          tmpFace(:,:,1) = abs(fhe(:,:,1))/gs(:,:,1)
          u = maxval(tmpFace(ii,jsep+1:ny,0))/2.0_IDS_real
          j = maxloc(tmpFace(ii,jsep+1:ny,0),dim=1)
          i = jsep+j
          do iy = jsep+j+1, ny
            if (i.gt.jsep+j) cycle
            if (tmpFace(ii,iy,0).lt.u) then
              i = iy
              frac = (tmpFace(ii,iy-1,0)-u)/(tmpFace(ii,iy-1,0)-tmpFace(ii,iy,0))
              lambda = lambda + frac*hy(jxa,iy-1)*qz(jxa,iy-1,1)
            else
              lambda = lambda + hy(jxa,iy-1)*qz(jxa,iy-1,1)
            end if
          end do
          lambda = lambda/log(2.0_IDS_real)
          summary%scrape_off_layer%heat_flux_e_decay_length%value( time_sind ) = lambda
          lambda = 0.0_IDS_real
          tmpFace(:,:,0) = abs(fhi(:,:,0))/gs(:,:,0)/qc(:,:)
          tmpFace(:,:,1) = abs(fhi(:,:,1))/gs(:,:,1)
          u = maxval(tmpFace(jj,jsep+1:ny,0))/2.0_IDS_real
          j = maxloc(tmpFace(jj,jsep+1:ny,0),dim=1)
          i = jsep+j
          do iy = jsep+j+1, ny
            if (i.gt.jsep+j) cycle
            if (tmpFace(jj,iy,0).lt.u) then
              i = iy
              frac = (tmpFace(jj,iy-1,0)-u)/(tmpFace(jj,iy-1,0)-tmpFace(jj,iy,0))
              lambda = lambda + frac*hy(jxa,iy-1)*qz(jxa,iy-1,1)
            else
              lambda = lambda + hy(jxa,iy-1)*qz(jxa,iy-1,1)
            end if
          end do
          lambda = lambda/log(2.0_IDS_real)
          summary%scrape_off_layer%heat_flux_i_decay_length%value( time_sind ) = lambda
          allocate ( summary%scrape_off_layer%heat_flux_e_decay_length%source(1) )
          allocate ( summary%scrape_off_layer%heat_flux_i_decay_length%source(1) )
          summary%scrape_off_layer%heat_flux_e_decay_length%source = source
          summary%scrape_off_layer%heat_flux_i_decay_length%source = source
        end if
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

        deallocate(isstat,imneut)
#ifdef B25_EIRENE
        deallocate(in_species)
#endif
        call logmsg( LOGDEBUG, "b2mod_ual_io.B25_process_ids: done" )

        contains

        integer function get_atom_number( compname )
            implicit none
            character*2 compname
            integer is, iatm

            get_atom_number = 0
            do iatm = 1, natmi
               if ( get_atom_number > 0 ) cycle
               is = eb2atcr(iatm)
               if (streql( is_codes( is ), compname ) ) get_atom_number = iatm
            end do
            return

        end function get_atom_number

#if IMAS_MINOR_VERSION > 11
        !> Write scalar B2 cell quantity to 'ids_generic_grid_scalar'
        !! IMAS IDS data tree node.
        subroutine write_quantity( val, value, time_sind )
            use b2mod_interp
            type(ids_generic_grid_scalar), pointer, intent(inout) :: val(:)
                !< Type of IDS data structure, designed for scalar data handling
            real(IDS_real), intent(in) :: value( -1:gmap%b2nx, -1:gmap%b2ny )
            integer, intent(in) :: time_sind    !< General grid description
                                                !< slice identifier
            real(IDS_real), dimension(:), pointer :: idsdata    !< Array for
                !< handing data field values
            real(IDS_real) :: weight( -1:gmap%b2nx, -1:gmap%b2ny, TO_SELF:TO_TOP )
            integer :: nSubsets  !< number of grid subsets to fill
            integer :: iSubset   !< Grid subset iterator
            integer :: iSubsetID !< Grid subset identifier index
            integer :: ggdID     !< Grid identifier index
            integer :: ndim      !< Grid subset dimension
            integer :: i         !< Iterator

#if IMAS_MINOR_VERSION < 15
            ggdId = edge_profiles%ggd(time_sind)%grid%identifier%index
            !! Assign 5+4 grid subsets
            nSubsets = 9
#else
            ggdId = edge_profiles%grid_ggd(time_sind)%identifier%index
            nSubsets = size(edge_profiles%grid_ggd(time_sind)%grid_subset)
#endif
            !! Allocate data fields for grid subsets
            allocate( val(nSubsets) )

            do iSubset = 1, nSubsets
#if IMAS_MINOR_VERSION < 15
               select case (iSubset)
                  case (1)
                    iSubsetID = GRID_SUBSET_CELLS
                    ndim = 3
                  case (2)
                    iSubsetID = iGsCoreBoundary
                    ndim = 2
                  case (3)
                    iSubsetID = iGsInnerMidplane
                    ndim = 1
                  case (4)
                    iSubsetID = iGsOuterMidplane
                    ndim = 1
                  case (5)
                    iSubsetID = GRID_SUBSET_NODES
                    ndim = 1
                  case (6)
                    iSubsetID = iGsCore
                    ndim = 3
                  case (7)
                    iSubsetID = iGsSOL
                    ndim = 3
                  case (8)
                    iSubsetID = iGsIDivertor
                    ndim = 3
                  case (9)
                    iSubsetID = iGsODivertor
                    ndim = 3
                  case default
                    iSubsetID = iSubset
                    ndim = IDS_INT_INVALID
               end select
#else
               ndim = edge_profiles%grid_ggd(time_sind)%grid_subset(iSubset)%dimension
               iSubsetID = edge_profiles%grid_ggd(time_sind)%grid_subset(iSubset)%identifier%index
#endif
               if (ndim.eq.IDS_INT_INVALID) then
                  select case (iSubsetID)
                  case( GRID_SUBSET_NODES, GRID_SUBSET_X_POINTS, &
                      & GRID_SUBSET_INNER_MIDPLANE, GRID_SUBSET_OUTER_MIDPLANE )
                     ndim = 1
                  case( GRID_SUBSET_FACES, &
                      & GRID_SUBSET_X_ALIGNED_FACES, GRID_SUBSET_Y_ALIGNED_FACES, &
                      & GRID_SUBSET_CORE_BOUNDARY, GRID_SUBSET_SEPARATRIX, &
                      & GRID_SUBSET_ACTIVE_SEPARATRIX, GRID_SUBSET_MAIN_CHAMBER_WALL, &
                      & GRID_SUBSET_OUTER_BAFFLE, GRID_SUBSET_INNER_BAFFLE, &
                      & GRID_SUBSET_OUTER_PFR_WALL, GRID_SUBSET_INNER_PFR_WALL, &
                      & GRID_SUBSET_MAIN_WALL, GRID_SUBSET_PFR_WALL, &
                      & GRID_SUBSET_SECOND_SEPARATRIX, &
                      & GRID_SUBSET_OUTER_BAFFLE_INACTIVE, &
                      & GRID_SUBSET_INNER_BAFFLE_INACTIVE, &
                      & GRID_SUBSET_OUTER_PFR_WALL_INACTIVE, &
                      & GRID_SUBSET_INNER_PFR_WALL_INACTIVE, &
                      & GRID_SUBSET_CORE_CUT, GRID_SUBSET_PFR_CUT, &
                      & GRID_SUBSET_OUTER_THROAT, GRID_SUBSET_INNER_THROAT, &
                      & GRID_SUBSET_OUTER_TARGET, GRID_SUBSET_INNER_TARGET, &
                      & GRID_SUBSET_CORE_CUT_INACTIVE, GRID_SUBSET_PFR_CUT_INACTIVE, &
                      & GRID_SUBSET_OUTER_THROAT_INACTIVE, &
                      & GRID_SUBSET_INNER_THROAT_INACTIVE, &
                      & GRID_SUBSET_OUTER_TARGET_INACTIVE, &
                      & GRID_SUBSET_INNER_TARGET_INACTIVE )
                     ndim = 2
                  case( GRID_SUBSET_CELLS, GRID_SUBSET_BETWEEN_SEPARATRICES, &
                      & GRID_SUBSET_CORE, GRID_SUBSET_SOL, &
                      & GRID_SUBSET_OUTER_DIVERTOR, GRID_SUBSET_INNER_DIVERTOR, &
                      & GRID_SUBSET_OUTER_DIVERTOR_INACTIVE, &
                      & GRID_SUBSET_INNER_DIVERTOR_INACTIVE )
                     ndim = 3
                  end select
               end if
               select case (ndim)
               case ( 1 ) !< Grid subset consists of nodes
                  tmpVx = interpolateToVertices(  &
                     &   gmap%b2nx, gmap%b2ny, VX_LOWERLEFT, value )
#if IMAS_MINOR_VERSION < 15
                  idsdata => b2_IMAS_Transform_Data_B2_To_IDS_Vertex(            &
                     &   edge_profiles%ggd( time_sind )%grid, iGsOuterMidplane,  &
                     &   gmap, tmpVx )
#else
                  idsdata => b2_IMAS_Transform_Data_B2_To_IDS_Vertex( &
                     &   edge_profiles%grid_ggd( time_sind ),         &
                     &   iSubset, gmap, tmpVx )
#endif
#if GGD_MINOR_VERSION > 8
                  call gridWriteData( val( iSubset ), ggdID, iSubsetID, idsdata )
#else
                  call gridWriteData( val( iSubset ), iSubsetID, idsdata )
#endif
                  deallocate( idsdata )
               case ( 2 ) !< Grid subset consists of faces
                  tmpFace = 0.0_IDS_real
                  do i = TO_SELF, TO_TOP
                     weight(:,:,i) = vol(:,:)
                  end do
                  call value_on_faces( nx, ny, weight, value, tmpFace)
#if IMAS_MINOR_VERSION < 15
                  idsdata => b2_IMAS_Transform_Data_B2_To_IDS(             &
                     &   edge_profiles%ggd( time_sind )%grid, iSubset,     &
                     &   gmap, tmpFace )
#else
                  idsdata => b2_IMAS_Transform_Data_B2_To_IDS(             &
                     &   edge_profiles%grid_ggd( time_sind ), iSubset,     &
                     &   gmap, tmpFace )
#endif
#if GGD_MINOR_VERSION > 8
                  call gridWriteData( val( iSubset ), ggdID, iSubsetID, idsdata )
#else
                  call gridWriteData( val( iSubset ), iSubsetID, idsdata )
#endif
                  deallocate( idsdata )
               case ( 3 ) !< Grid subset consists of cells
#if IMAS_MINOR_VERSION < 15
                  idsdata => b2_IMAS_Transform_Data_B2_To_IDS( edge_profiles%   &
                      &   ggd( time_sind )%grid, iSubset, gmap, value )
#else
                  idsdata => b2_IMAS_Transform_Data_B2_To_IDS( edge_profiles%   &
                      &   grid_ggd( time_sind ), iSubset, gmap, value )
#endif
#if GGD_MINOR_VERSION > 8
                  call gridWriteData( val( iSubset ), ggdID, iSubsetID, idsdata )
#else
                  call gridWriteData( val( iSubset ), iSubsetID, idsdata )
#endif
                  deallocate( idsdata )
               case default
                  call xerrab( 'Unknown grid subset '//int2str(iSubset)// &
                      &        ' dimension : '//int2str(ndim) )
               end select
            end do

        end subroutine write_quantity

        !> Write a scalar B2 face quantity to ids_generic_grid_scalar
        subroutine write_face_scalar( val, value, time_sind )
            type (ids_generic_grid_scalar), pointer, intent(inout) :: val(:)
                !< Type of IDS data structure, designed for scalar data handling
                !< (in this case scalars residing on grid faces)
            real(IDS_real), intent(in) :: value( -1:gmap%b2nx, -1:gmap%b2ny, 0:1 )
            integer, intent(in) :: time_sind    !< General grid description
                                                !< slice identifier
            integer :: nSubsets  !< number of grid subsets to fill
            integer :: iSubset   !< Grid subset iterator
            integer :: iSubsetID !< Grid subset identifier index
            integer :: ggdID     !< Grid identifier index
            integer :: ndim      !< Grid subset dimension

#if IMAS_MINOR_VERSION < 15
            ggdId = edge_profiles%ggd(time_sind)%grid%identifier%index
            nSubsets = 2
#else
            ggdId = edge_profiles%grid_ggd(time_sind)%identifier%index
            nSubsets = size(edge_profiles%grid_ggd(time_sind)%grid_subset)
#endif
            !! Allocate data fields for grid subsets
            allocate( val(nSubsets) )

            do iSubset = 1, nSubsets
#if IMAS_MINOR_VERSION < 15
               select case (iSubset)
                  case (1)
                    iSubsetID = GRID_SUBSET_X_ALIGNED_FACES
                    ndim = 2
                  case (2)
                    iSubsetID = GRID_SUBSET_Y_ALIGNED_FACES
                    ndim = 2
                  case default
                    iSubsetID = iSubset
                    ndim = IDS_INT_INVALID
               end select
#else
               ndim = edge_profiles%grid_ggd(time_sind)%grid_subset(iSubset)%dimension
               iSubsetID = edge_profiles%grid_ggd(time_sind)%grid_subset(iSubset)%identifier%index
#endif
               if (ndim.eq.IDS_INT_INVALID) then
                  select case (iSubsetID)
                  case( GRID_SUBSET_NODES, GRID_SUBSET_X_POINTS, &
                      & GRID_SUBSET_INNER_MIDPLANE, GRID_SUBSET_OUTER_MIDPLANE )
                     ndim = 1
                  case( GRID_SUBSET_FACES, &
                      & GRID_SUBSET_X_ALIGNED_FACES, GRID_SUBSET_Y_ALIGNED_FACES, &
                      & GRID_SUBSET_CORE_BOUNDARY, GRID_SUBSET_SEPARATRIX, &
                      & GRID_SUBSET_ACTIVE_SEPARATRIX, GRID_SUBSET_MAIN_CHAMBER_WALL, &
                      & GRID_SUBSET_OUTER_BAFFLE, GRID_SUBSET_INNER_BAFFLE, &
                      & GRID_SUBSET_OUTER_PFR_WALL, GRID_SUBSET_INNER_PFR_WALL, &
                      & GRID_SUBSET_MAIN_WALL, GRID_SUBSET_PFR_WALL, &
                      & GRID_SUBSET_SECOND_SEPARATRIX, &
                      & GRID_SUBSET_OUTER_BAFFLE_INACTIVE, &
                      & GRID_SUBSET_INNER_BAFFLE_INACTIVE, &
                      & GRID_SUBSET_OUTER_PFR_WALL_INACTIVE, &
                      & GRID_SUBSET_INNER_PFR_WALL_INACTIVE, &
                      & GRID_SUBSET_CORE_CUT, GRID_SUBSET_PFR_CUT, &
                      & GRID_SUBSET_OUTER_THROAT, GRID_SUBSET_INNER_THROAT, &
                      & GRID_SUBSET_OUTER_TARGET, GRID_SUBSET_INNER_TARGET, &
                      & GRID_SUBSET_CORE_CUT_INACTIVE, GRID_SUBSET_PFR_CUT_INACTIVE, &
                      & GRID_SUBSET_OUTER_THROAT_INACTIVE, &
                      & GRID_SUBSET_INNER_THROAT_INACTIVE, &
                      & GRID_SUBSET_OUTER_TARGET_INACTIVE, &
                      & GRID_SUBSET_INNER_TARGET_INACTIVE )
                     ndim = 2
                  case( GRID_SUBSET_CELLS, GRID_SUBSET_BETWEEN_SEPARATRICES, &
                      & GRID_SUBSET_CORE, GRID_SUBSET_SOL, &
                      & GRID_SUBSET_OUTER_DIVERTOR, GRID_SUBSET_INNER_DIVERTOR, &
                      & GRID_SUBSET_OUTER_DIVERTOR_INACTIVE, &
                      & GRID_SUBSET_INNER_DIVERTOR_INACTIVE )
                     ndim = 3
                  end select
               end if
               if (ndim.ne.2) cycle
               call write_face_vector( val( iSubset ), value, time_sind, ggdID, iSubsetID, iSubset )
            end do

        end subroutine write_face_scalar

        !> Write a scalar B2 cell quantity to ids_generic_grid_scalar
        subroutine write_cell_scalar( scalar, b2CellData )
            type(ids_generic_grid_scalar), intent(inout), pointer :: scalar(:)
                !< Type of IDS data structure, designed for scalar data handling
            real(IDS_real), intent(in) :: b2CellData(-1:gmap%b2nx, -1:gmap%b2ny)
            real(IDS_real), dimension(:), pointer :: idsdata    !< Array for
                !< handing data field values
            integer :: nSubsets  !< number of grid subsets to fill
            integer :: iSubset   !< Grid subset iterator
            integer :: ndim      !< Grid subset dimension
            integer :: iSubsetID !< Grid subset identifier index
            integer :: ggdID     !< Grid identifier index

#if IMAS_MINOR_VERSION < 15
            ggdId = edge_profiles%ggd(time_sind)%grid%identifier%index
            nSubsets = 1
#else
            ggdId = edge_profiles%grid_ggd(time_sind)%identifier%index
            nSubsets = size(edge_profiles%grid_ggd(time_sind)%grid_subset)
#endif
            !! Allocate data fields for grid subsets
            allocate( scalar(nSubsets) )

            do iSubset = 1, nSubsets
#if IMAS_MINOR_VERSION < 15
               ndim = 3
               iSubsetID = GRID_SUBSET_CELLS
#else
               ndim = edge_profiles%grid_ggd(time_sind)%grid_subset(iSubset)%dimension
               iSubsetID = edge_profiles%grid_ggd(time_sind)%grid_subset(iSubset)%identifier%index
#endif
               if (ndim.eq.IDS_INT_INVALID) then
                  select case (iSubsetID)
                  case( GRID_SUBSET_NODES, GRID_SUBSET_X_POINTS, &
                      & GRID_SUBSET_INNER_MIDPLANE, GRID_SUBSET_OUTER_MIDPLANE )
                     ndim = 1
                  case( GRID_SUBSET_FACES, &
                      & GRID_SUBSET_X_ALIGNED_FACES, GRID_SUBSET_Y_ALIGNED_FACES, &
                      & GRID_SUBSET_CORE_BOUNDARY, GRID_SUBSET_SEPARATRIX, &
                      & GRID_SUBSET_ACTIVE_SEPARATRIX, GRID_SUBSET_MAIN_CHAMBER_WALL, &
                      & GRID_SUBSET_OUTER_BAFFLE, GRID_SUBSET_INNER_BAFFLE, &
                      & GRID_SUBSET_OUTER_PFR_WALL, GRID_SUBSET_INNER_PFR_WALL, &
                      & GRID_SUBSET_MAIN_WALL, GRID_SUBSET_PFR_WALL, &
                      & GRID_SUBSET_SECOND_SEPARATRIX, &
                      & GRID_SUBSET_OUTER_BAFFLE_INACTIVE, &
                      & GRID_SUBSET_INNER_BAFFLE_INACTIVE, &
                      & GRID_SUBSET_OUTER_PFR_WALL_INACTIVE, &
                      & GRID_SUBSET_INNER_PFR_WALL_INACTIVE, &
                      & GRID_SUBSET_CORE_CUT, GRID_SUBSET_PFR_CUT, &
                      & GRID_SUBSET_OUTER_THROAT, GRID_SUBSET_INNER_THROAT, &
                      & GRID_SUBSET_OUTER_TARGET, GRID_SUBSET_INNER_TARGET, &
                      & GRID_SUBSET_CORE_CUT_INACTIVE, GRID_SUBSET_PFR_CUT_INACTIVE, &
                      & GRID_SUBSET_OUTER_THROAT_INACTIVE, &
                      & GRID_SUBSET_INNER_THROAT_INACTIVE, &
                      & GRID_SUBSET_OUTER_TARGET_INACTIVE, &
                      & GRID_SUBSET_INNER_TARGET_INACTIVE )
                     ndim = 2
                  case( GRID_SUBSET_CELLS, GRID_SUBSET_BETWEEN_SEPARATRICES, &
                      & GRID_SUBSET_CORE, GRID_SUBSET_SOL, &
                      & GRID_SUBSET_OUTER_DIVERTOR, GRID_SUBSET_INNER_DIVERTOR, &
                      & GRID_SUBSET_OUTER_DIVERTOR_INACTIVE, &
                      & GRID_SUBSET_INNER_DIVERTOR_INACTIVE )
                     ndim = 3
                  end select
               end if
               if (ndim.ne.3) cycle

            !! TODO: add checks whether already allocated
#if IMAS_MINOR_VERSION < 15
               idsdata => b2_IMAS_Transform_Data_B2_To_IDS( edge_profiles% &
                &   ggd( time_sind )%grid, iSubset,                        &
                &   gmap, b2CellData )
#else
               idsdata => b2_IMAS_Transform_Data_B2_To_IDS( edge_profiles% &
                &   grid_ggd( time_sind ), iSubset,                        &
                &   gmap, b2CellData )
#endif
#if GGD_MINOR_VERSION > 8
               call gridWriteData( scalar( iSubset ), ggdID, iSubsetID, idsdata )
#else
               call gridWriteData( scalar( iSubset ), iSubsetID, idsdata )
#endif
               deallocate(idsdata)
            end do

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
            integer :: nSubsets  !< number of grid subsets to fill
            integer :: iSubset   !< Grid subset iterator
            integer :: iSubsetID !< Grid subset identifier index
            integer :: ndim      !< Grid subset dimension
            integer :: ggdID     !< Grid identifier index

#if IMAS_MINOR_VERSION < 15
            ggdId = edge_profiles%ggd(time_sind)%grid%identifier%index
            nSubsets = 1
#else
            ggdId = edge_profiles%grid_ggd(time_sind)%identifier%index
            nSubsets = size(edge_profiles%grid_ggd(time_sind)%grid_subset)
#endif
            !! If required, allocate storage
            if ( .not. associated( vectorComponent ) ) then
                allocate( vectorComponent(nSubsets) )
            end if

            do iSubset = 1, nSubsets
#if IMAS_MINOR_VERSION < 15
               ndim = 3
               iSubsetID = GRID_SUBSET_CELLS
#else
               ndim = edge_profiles%grid_ggd(time_sind)%grid_subset(iSubset)%dimension
               iSubsetID = edge_profiles%grid_ggd(time_sind)%grid_subset(iSubset)%identifier%index
#endif
               if (ndim.eq.IDS_INT_INVALID) then
                  select case (iSubsetID)
                  case( GRID_SUBSET_NODES, GRID_SUBSET_X_POINTS, &
                      & GRID_SUBSET_INNER_MIDPLANE, GRID_SUBSET_OUTER_MIDPLANE )
                     ndim = 1
                  case( GRID_SUBSET_FACES, &
                      & GRID_SUBSET_X_ALIGNED_FACES, GRID_SUBSET_Y_ALIGNED_FACES, &
                      & GRID_SUBSET_CORE_BOUNDARY, GRID_SUBSET_SEPARATRIX, &
                      & GRID_SUBSET_ACTIVE_SEPARATRIX, GRID_SUBSET_MAIN_CHAMBER_WALL, &
                      & GRID_SUBSET_OUTER_BAFFLE, GRID_SUBSET_INNER_BAFFLE, &
                      & GRID_SUBSET_OUTER_PFR_WALL, GRID_SUBSET_INNER_PFR_WALL, &
                      & GRID_SUBSET_MAIN_WALL, GRID_SUBSET_PFR_WALL, &
                      & GRID_SUBSET_SECOND_SEPARATRIX, &
                      & GRID_SUBSET_OUTER_BAFFLE_INACTIVE, &
                      & GRID_SUBSET_INNER_BAFFLE_INACTIVE, &
                      & GRID_SUBSET_OUTER_PFR_WALL_INACTIVE, &
                      & GRID_SUBSET_INNER_PFR_WALL_INACTIVE, &
                      & GRID_SUBSET_CORE_CUT, GRID_SUBSET_PFR_CUT, &
                      & GRID_SUBSET_OUTER_THROAT, GRID_SUBSET_INNER_THROAT, &
                      & GRID_SUBSET_OUTER_TARGET, GRID_SUBSET_INNER_TARGET, &
                      & GRID_SUBSET_CORE_CUT_INACTIVE, GRID_SUBSET_PFR_CUT_INACTIVE, &
                      & GRID_SUBSET_OUTER_THROAT_INACTIVE, &
                      & GRID_SUBSET_INNER_THROAT_INACTIVE, &
                      & GRID_SUBSET_OUTER_TARGET_INACTIVE, &
                      & GRID_SUBSET_INNER_TARGET_INACTIVE )
                     ndim = 2
                  case( GRID_SUBSET_CELLS, GRID_SUBSET_BETWEEN_SEPARATRICES, &
                      & GRID_SUBSET_CORE, GRID_SUBSET_SOL, &
                      & GRID_SUBSET_OUTER_DIVERTOR, GRID_SUBSET_INNER_DIVERTOR, &
                      & GRID_SUBSET_OUTER_DIVERTOR_INACTIVE, &
                      & GRID_SUBSET_INNER_DIVERTOR_INACTIVE )
                     ndim = 3
                  end select
               end if
               if (ndim.ne.3) cycle
#if IMAS_MINOR_VERSION < 15
               idsdata => b2_IMAS_Transform_Data_B2_To_IDS( edge_profiles% &
                   &   ggd( time_sind )%grid, iSubset,                     &
                   &   gmap, b2CellData )
#else
               idsdata => b2_IMAS_Transform_Data_B2_To_IDS( edge_profiles% &
                   &   grid_ggd( time_sind ), iSubset,                     &
                   &   gmap, b2CellData )
#endif
               call B2grid_Write_Data_Vector_Components( vectorComponent(iSubset), &
                   &   ggdID, iSubsetID, vectorID, idsdata )
               deallocate(idsdata)
            end do

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
                &   grid_index, grid_subset_index, vectorID, data)
            type(ids_generic_grid_vector_components), intent(inout) ::  &
                &   idsField_vcomp
                !< Type of IDS data structure, designed for handling data
                !> regarding vector components (parallel, poloidal etc.)
            integer, intent(in) :: grid_index           !< Grid index
            integer, intent(in) :: grid_subset_index    !< Base grid subset
                                                        !< index
            character(len=*), intent(in) :: vectorID    !< Vector ID (e.g. )
                                                        !< VEC_ALIGN_RADIAL_ID)
            real(ids_real), intent(in) :: data(:)   !< Data field to be written
                !< to IDS data structure leaf that corresponds to specified
                !< vector component

            !! set grid index
            idsField_vcomp%grid_index = grid_index

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
#endif
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
            integer :: ggdId

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

            ggdId = edge_profiles%grid_ggd(time_sind)%identifier%index
            !! Fill in vector component data
            do i = 1, dim
#if GGD_MINOR_VERSION < 9
                idsdata => b2_IMAS_Transform_Data_B2_To_IDS(        &
                    &   edge_profiles%ggd( time_sind )%grid,        &
                    &   GRID_SUBSET_CELLS, gmap, vecdata(:,:,i-1))
#else
                idsdata => b2_IMAS_Transform_Data_B2_To_IDS(        &
                    &   edge_profiles%grid_ggd( time_sind ),        &
                    &   GRID_SUBSET_CELLS, gmap, vecdata(:,:,i-1))
#endif
#if GGD_MINOR_VERSION > 8
                call gridWriteData( vector, ggdId, GRID_SUBSET_CELLS, idsdata )
#else
                call gridWriteData( vector, GRID_SUBSET_CELLS, idsdata )
#endif
                deallocate(idsdata)
            end do

        end subroutine write_cell_vector
#endif
#if IMAS_MINOR_VERSION > 11
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
                &   gridID, gridSubsetID, gridSubsetInd )
            use ids_grid_data ! IGNORE
            type(ids_generic_grid_scalar), intent(inout) :: vector
                !< Type of IDS data structure, designed for scalar data handling
                !< (in this case 1D vector)
            real(IDS_real), intent(in) :: &
                &   b2FaceData(-1:gmap%b2nx, -1:gmap%b2ny, 0:1)
            integer, intent(in) :: gridID                    !< Grid identifier index
            integer, intent(in), optional :: gridSubsetID    !< Grid subset identifier index
            integer, intent(in), optional :: gridSubsetInd   !< Base grid subset index
            real(IDS_real), dimension(:), pointer :: idsdata !< Dummy array
                !< for holding data field values
            integer, intent(in) :: time_sind    !< General grid description

            if ( .not. present(gridSubsetInd) ) then
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
#if GGD_MINOR_VERSION > 8
                call gridWriteData( vector, gridId, GRID_SUBSET_Y_ALIGNED_FACES, idsdata )
#else
                call gridWriteData( vector, GRID_SUBSET_Y_ALIGNED_FACES, idsdata )
#endif
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
#if GGD_MINOR_VERSION > 8
                call gridWriteData( vector, gridId, GRID_SUBSET_X_ALIGNED_FACES, idsdata )
#else
                call gridWriteData( vector, GRID_SUBSET_X_ALIGNED_FACES, idsdata )
#endif
                deallocate(idsdata)
            else
#if IMAS_MINOR_VERSION < 15
                idsdata => b2_IMAS_Transform_Data_B2_To_IDS(    &
                    &   edge_profiles%ggd( time_sind )%grid,    &
                    &   gridSubsetInd, gmap, b2FaceData)
#else
                idsdata => b2_IMAS_Transform_Data_B2_To_IDS(    &
                    &   edge_profiles%grid_ggd( time_sind ),    &
                    &   gridSubsetInd, gmap, b2FaceData)
#endif
#if GGD_MINOR_VERSION > 8
                call gridWriteData( vector, gridId, gridSubsetID, idsdata )
#else
                call gridWriteData( vector, gridSubsetID, idsdata )
#endif
                deallocate(idsdata)
            end if

        end subroutine write_face_vector
#endif
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

        !! write the magnetic field vector in the B2 coordinate system
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
      real(ITM_R8) :: weight(-1:gmap%b2nx, -1:gmap%b2ny, TO_SELF:TO_TOP)
      integer i

      allocate(values(5))
      cpodata => b2ITMTransformDataB2ToCpo( edgecpo%grid, B2_SUBGRID_CELLS, gmap, value )
      call gridWriteData( values(1), B2_SUBGRID_CELLS, cpodata )
      deallocate(cpodata)
      tmpFace = 0.0_ITM_R8
      do i = TO_SELF, TO_TOP
        weight(:,:,i)=vol(:,:)
      end do
      call value_on_faces(nx,ny,weight,value,tmpFace)
      cpodata => b2ITMTransformDataB2ToCpo( edgecpo%grid, iSgCore, gmap, tmpFace )
      call gridWriteData( values(2), iSgCore, cpodata )
      deallocate(cpodata)
      tmpVx = interpolateToVertices( gmap%b2nx, gmap%b2ny, VX_LOWERLEFT, value )
      cpodata => b2ITMTransformDataB2ToCpoVertex( edgecpo%grid, iSgInnerMidplane, gmap, tmpVx )
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
