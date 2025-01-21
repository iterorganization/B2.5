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
    use b2us_map
    use b2us_geo
    use b2us_plasma
    use b2mod_elements
    use b2mod_constants
!    use b2mod_feedback
!    use b2mod_boundary_namelist
    use b2mod_neutrals_namelist
    use b2mod_user_namelist
    use b2mod_ipmain
    use b2mod_b2cmfs
    use b2mod_version
    use b2mod_switches
#ifdef IMAS
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
     &          b2stel_she_ion_bal, b2stel_shi_ion_bal, b2npht_shei_bal, &
     &          b2stel_she_rec_bal, b2stel_shi_rec_bal, &
     &          read_balance
    use b2mod_diag &
     & , only : nfluids, species_list, label
#endif
!    use b2mod_b2plot_wall_loading
#ifdef B25_EIRENE
!    use eirmod_ctrig
!    use eirmod_cestim
#ifdef IMAS
    use eirmod_cinit &
     & , only : fort_lc
    use eirmod_comusr &
     & , only : natmi, nmoli, nioni, nmassa, nchara, nmassm, ncharm, &
     &          nprt, nchrgi, nchari
    use b2mod_b2plot &
     & , only : triangle_vol, wklng, alloc_b2mod_b2plot_eirene
#endif
#elif defined(IMAS)
    use b2mod_b2plot &
     & , only : natmi
#endif
    use logging

#ifdef IMAS
    !! UAL Access
    use b2mod_ual_io_grid &
     & , only : INCLUDE_GHOST_CELLS, US_GRID_UNDEFINED,  &
     &          REGIONTYPE_EDGE,                         &
     &          GEOMETRY_LINEAR, GEOMETRY_CYLINDER,      &
     &          GEOMETRY_ANNULUS, GEOMETRY_LIMITER,      &
     &          GEOMETRY_SN, GEOMETRY_STELLARATORISLAND, &
     &          GEOMETRY_CDN,                            &
     &          GEOMETRY_DDN_BOTTOM, GEOMETRY_DDN_TOP,   &
     &          GEOMETRY_LFS_SNOWFLAKE_PLUS,             &
     &          GEOMETRY_LFS_SNOWFLAKE_MINUS,            &
     &          geometryId, geometryName, regionNumbers
#if ( IMAS_MINOR_VERSION > 11 || IMAS_MAJOR_VERSION > 3 ) && GGD_MAJOR_VERSION > 0
    !! B2/CPO Mapping
    use b2mod_ual_io_data &
     & , only : b2_IMAS_Transform_Data_B2_To_IDS_Cell,  &
     &          b2_IMAS_Transform_Data_B2_To_IDS_Face,  &
     &          b2_IMAS_Transform_Data_B2_To_IDS_Vertex
    use b2mod_ual_io_grid &
     & , only : b2_IMAS_Fill_Grid_Desc
    use ids_grid_subgrid  &     ! IGNORE
     & , only : findGridSubsetByName
    use ids_grid_structured &   ! IGNORE
     & , only : gridWriteData
    use ids_grid_common , &     ! IGNORE
        &   IDS_COORDTYPE_R => COORDTYPE_R,       &
        &   IDS_COORDTYPE_Z => COORDTYPE_Z,       &
        &   IDS_GRID_UNDEFINED => GRID_UNDEFINED
#endif
#if GGD_MAJOR_VERSION < 1
    use b2mod_ual_io_grid &
     & , only : VEC_ALIGN_RADIAL_ID,   &
     &          VEC_ALIGN_POLOIDAL_ID, &
     &          VEC_ALIGN_PARALLEL_ID, &
     &          VEC_ALIGN_TOROIDAL_ID, &
     &          VEC_ALIGN_DIAMAGNETIC_ID
#endif
#if ( GGD_MAJOR_VERSION < 2 && GGD_MINOR_VERSION < 11 && GGD_MICRO_VERSION < 3 )
    use b2mod_ual_io_grid &
     & , only : VEC_ALIGN_R_MAJOR_ID,   &
     &          VEC_ALIGN_Z_ID
#endif
#if ( GGD_MAJOR_VERSION < 2 && GGD_MINOR_VERSION < 13 )
    use b2mod_ual_io_grid &
     & , only : VEC_ALIGN_PHI_ID
#endif
#if GGD_MAJOR_VERSION > 0
#if ( GGD_MINOR_VERSION < 9 && GGD_MAJOR_VERSION < 2 )
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
#if ( GGD_MINOR_VERSION < 10 && GGD_MAJOR_VERSION == 1 )
    use b2mod_ual_io_grid &
     & , only : GRID_SUBSET_X_ALIGNED_EDGES, GRID_SUBSET_Y_ALIGNED_EDGES, &
     &          GRID_SUBSET_EDGES, GRID_SUBSET_VOLUMES
#endif
#if ( ( GGD_MINOR_VERSION < 10 || ( GGD_MINOR_VERSION == 10 && GGD_MICRO_VERSION < 2 ) ) && GGD_MAJOR_VERSION == 1 )
    use b2mod_ual_io_grid &
     & , only : GRID_SUBSET_MAGNETIC_AXIS, GRID_SUBSET_FULL_WALL
#endif
#if ( ( GGD_MINOR_VERSION < 10 || ( GGD_MINOR_VERSION == 10 && GGD_MICRO_VERSION < 3 ) ) && GGD_MAJOR_VERSION == 1 )
    use b2mod_ual_io_grid &
     & , only : GRID_SUBSET_OUTER_SF_LEG_ENTRANCE_1, &
     &          GRID_SUBSET_OUTER_SF_LEG_ENTRANCE_2,  &
     &          GRID_SUBSET_OUTER_SF_PFR_CONNECTION_1, &
     &          GRID_SUBSET_OUTER_SF_PFR_CONNECTION_2
#endif
#endif
    use ids_schemas &     ! IGNORE
     & , only : ids_string_length
#if ( IMAS_MINOR_VERSION > 8 || IMAS_MAJOR_VERSION > 3 )
    use ids_schemas &     ! IGNORE
     & , only : ids_real, ids_real_invalid
#endif
    use ids_schemas &     ! IGNORE
     & , only : ids_edge_profiles, ids_edge_sources, ids_edge_transport,     &
     &          ids_radiation, ids_dataset_description, ids_equilibrium,     &
     &          ids_ids_properties,                                          &
     &          ids_code, ids_signal_int_1d, ids_signal_flt_1d,              &
     &          ids_generic_grid_scalar, ids_generic_grid_vector_components, &
     &          ids_generic_grid_dynamic
#if ( IMAS_MINOR_VERSION > 14 || IMAS_MAJOR_VERSION > 3 )
    use ids_schemas &     ! IGNORE
     & , only : ids_generic_grid_AoS3_root
#endif
#if ( IMAS_MINOR_VERSION > 21 || IMAS_MAJOR_VERSION > 3 )
    use ids_schemas &     ! IGNORE
     & , only : ids_summary,                                                        &
     &          ids_summary_constant_flt_0d, ids_summary_constant_int_0d,           &
     &          ids_summary_dynamic_int_1d_root, ids_summary_dynamic_flt_1d_root,   &
     &          ids_summary_dynamic_flt_1d_root_parent_2, ids_summary_static_str_0d
#endif
#if ( IMAS_MINOR_VERSION > 25 && IMAS_MINOR_VERSION < 34 && IMAS_MAJOR_VERSION == 3 )
    use ids_schemas &     ! IGNORE
     & , only : ids_numerics
#endif
#if ( IMAS_MINOR_VERSION > 30 || IMAS_MAJOR_VERSION > 3 )
    use ids_schemas &     ! IGNORE
     & , only : ids_divertors
#if AL_MAJOR_VERSION < 5
    use imas_midplane_identifier &     ! IGNORE
     & , only : midplane_identifier
    use imas_neutrals_identifier &     ! IGNORE
     & , only : neutrals_identifier
    use imas_radiation_identifier &    ! IGNORE
     & , only : radiation_identifier
    use imas_edge_source_identifier &  ! IGNORE
     & , only : edge_source_identifier
#else
    use al_midplane_identifier &       ! IGNORE
     & , only : midplane_identifier
    use al_neutrals_identifier &       ! IGNORE
     & , only : neutrals_identifier
    use al_radiation_identifier &      ! IGNORE
     & , only : radiation_identifier
    use al_edge_source_identifier &    ! IGNORE
     & , only : edge_source_identifier
#endif
#endif
#if ( IMAS_MINOR_VERSION > 32 || IMAS_MAJOR_VERSION > 3 )
    use ids_utilities &   ! IGNORE
     & , only : ids_identifier_static
#endif
#if ( IMAS_MINOR_VERSION > 36 || IMAS_MAJOR_VERSION > 3 )
    use b2mod_math
    use ids_schemas &     ! IGNORE
     & , only : ids_summary_rz1d_dynamic
#endif
#if IMAS_MAJOR_VERSION > 3
    use ids_schemas &     ! IGNORE
     & , only : ids_code_constant
#endif
#if ( defined(AMNS) && ( IMAS_MINOR_VERSION > 29 || IMAS_MAJOR_VERSION > 3 ) )
    use amns_types  ! IGNORE
    use amns_module ! IGNORE
#endif
#elif defined(ITM_ENVIRONMENT_LOADED)
    use euITM_schemas   ! IGNORE
    use euITM_routines  ! IGNORE
    use itm_grid_common ! IGNORE
#endif

  public b25_process_ids, b25_av_ids
  integer, public :: num_time_slices  !< Total number of time slices.
  integer, public :: num_batch_slices !< Total number of batch slices.

  private

#ifdef IMAS

#if ( IMAS_MINOR_VERSION < 9 && IMAS_MAJOR_VERSION < 4 )
  integer, parameter :: IDS_REAL = R8
  real(kind=R8), parameter :: IDS_REAL_INVALID = -9.0E40_R8
#endif

  interface write_sourced_value
     module procedure write_sourced_value_root
     module procedure write_sourced_value_root_parent_2
  end interface write_sourced_value

  integer :: num_slices    !< Total number of time slices in IDS
  integer :: slice_index   !< Time slice index. Also General grid
            !< description slice identifier
  integer :: homogeneous_time !< Homogeneous time (0 or 1)
  integer, save :: nx    !< Specifies the number of interior cells
                         !< along the first coordinate
  integer, save :: ny    !< Specifies the number of interior cells
                         !< along the second coordinate
  integer, save :: midplane_id      !< Location of midplane:
                                    !< 1: Z equal to equilibrium O-point
                                    !< 2: Z at location of maximum major radius
                                    !< 3: Z at dR/dZ = 0 maximum R location
                                    !< 4: GGD grid subset defined by jxa value
  integer, save :: GeometryType !< Geometry identifier number
#if GGD_MAJOR_VERSION > 0
  integer, save :: iGsCoreBoundary  !< Variable to hold Core grid subset base
            !< index, later found by findGridSubsetByName() routine.
  integer, save :: iGsInnerMidplane !< Variable to hold Inner Midplane grid
            !< subset base index, later found by findGridSubsetByName() routine
  integer, save :: iGsOuterMidplane !< Variable to hold Outer Midplane grid
            !< subset base index, later found by findGridSubsetByName() routine
  integer, save :: iGsCore  !< Variable to hold Core grid
            !< subset base index, later found by findGridSubsetByName() routine
  integer, save :: iGsSOL   !< Variable to hold SOL grid
            !< subset base index, later found by findGridSubsetByName() routine
  integer, save :: iGsIDivertor     !< Variable to hold Inner Divertor grid
            !< subset base index, later found by findGridSubsetByName() routine
  integer, save :: iGsODivertor     !< Variable to hold Outer Divertor grid
            !< subset base index, later found by findGridSubsetByName() routine
#endif
  logical, parameter :: B2_WRITE_DATA = .true.
  real(IDS_real) :: time  !< Generic time
  real(IDS_real), save :: b0, r0, b0r0
  real(IDS_real), save, allocatable :: flux_expansion(:), &
      &  extension_r(:), extension_z(:), wetted_area(:)
  character(len=ids_string_length), save :: username  !< IDS user name
  character(len=ids_string_length), save :: source    !< Code source
  character(len=ids_string_length), save :: comment   !< IDS properties label
  character(len=ids_string_length), save :: create_date
  character(len=ids_string_length), save :: code_commit
  character(len=ids_string_length), save :: code_description
  character(len=ids_string_length), save :: configuration
  character(len=ids_string_length), save :: plate_name(4) !< Divertor plate name
  character*8, save :: imas_version, adas_version
  character*8, save :: date
  character*10, save :: ctime
  character*5, save :: zone
  character*32, save :: B25_git_version
  logical, save :: eq_found
#ifndef NO_GETENV
  integer lenval, ierror
#ifndef USE_PXFGETENV
  intrinsic get_environment_variable
#endif
#endif

contains

    subroutine IDS_init( mpg, geo )
    implicit none
    type( mapping ), intent(in) :: mpg
    type( geometry ), intent(in) :: geo
    integer tvalues(8)
    integer i, iCv, iFc
    real(IDS_real) :: r_min, r_max, z_min, z_max
    logical, save :: IDS_initialized = .false.
    character*16 usrnam
    character*32 get_B25_hash
    external usrnam
    external get_B25_hash

    if (IDS_initialized) return
    username = usrnam()
#ifdef NO_GETENV
    write(imas_version,'(i1,a1,i2,a1,i1)')  IMAS_MAJOR_VERSION,'.', &
                                      &     IMAS_MINOR_VERSION,'.', &
                                      &     IMAS_MICRO_VERSION
#elif defined(USE_PXFGETENV)
    CALL PXFGETENV ('IMAS_VERSION', 0, imas_version, lenval, ierror)
#else
    call get_environment_variable('IMAS_VERSION',status=ierror,length=lenval)
    if (ierror.eq.0) call get_environment_variable('IMAS_VERSION',value=imas_version)
#endif
    call date_and_time (date, ctime, zone, tvalues)
    create_date = date//' '//ctime//' '//' '//zone
    B25_git_version = get_B25_hash()
    code_commit = B25_git_version

    geometryType = geometryId( mpg, geo )
    configuration = geometryName( geometryType )
    select case ( geometryType )
    case ( GEOMETRY_LINEAR, GEOMETRY_CYLINDER )
      plate_name(1) = 'W'
      plate_name(2) = 'E'
    case ( GEOMETRY_LIMITER, GEOMETRY_SN, GEOMETRY_STELLARATORISLAND )
      if (geo%LSN) then
        plate_name(1) = 'W'
        plate_name(2) = 'E'
      else
        plate_name(1) = 'E'
        plate_name(2) = 'W'
      end if
    case ( GEOMETRY_CDN, GEOMETRY_DDN_TOP, GEOMETRY_DDN_BOTTOM )
      plate_name(1) = 'LI'
      plate_name(2) = 'UI'
      plate_name(3) = 'UO'
      plate_name(4) = 'LO'
    case ( GEOMETRY_LFS_SNOWFLAKE_PLUS, GEOMETRY_LFS_SNOWFLAKE_MINUS )
      plate_name(1) = 'E'
      plate_name(2) = 'W1'
      plate_name(3) = 'W2'
      plate_name(4) = 'W3'
    end select
    if (maxval(mpg%strDiv).gt.0) then
      allocate(wetted_area( maxval(mpg%strDiv) ) )
      allocate(extension_r( maxval(mpg%strDiv) ) )
      allocate(extension_z( maxval(mpg%strDiv) ) )
      allocate(flux_expansion( maxval(mpg%strDiv) ) )
    end if
    do i = 1, maxval(mpg%strDiv)
      r_min = huge(1.0_R8)
      r_max = 0.0_R8
      z_min = huge(1.0_R8)
      z_max =-huge(1.0_R8)
      wetted_area(i) = 0.0_IDS_real
      do iFc = 1, mpg%nFc
        if (mpg%fcReg(iFc).eq. &
          & regionNumbers(i,REGIONTYPE_EDGE,geometryType)) then
          if (mpg%fcCv(iFc,1).le.mpg%nCi) then
            iCv = mpg%fcCv(iFc,1)
          else
            iCv = mpg%fcCv(iFc,2)
          end if
          r_min = min( r_min, geo%vxX(mpg%fcVx(iFc,1)), &
               &              geo%vxX(mpg%fcVx(iFc,2)) )
          r_max = max( r_max, geo%vxX(mpg%fcVx(iFc,1)), &
               &              geo%vxX(mpg%fcVx(iFc,2)) )
          z_min = min( z_min, geo%vxY(mpg%fcVx(iFc,1)), &
               &              geo%vxY(mpg%fcVx(iFc,2)) )
          z_max = max( z_max, geo%vxY(mpg%fcVx(iFc,1)), &
               &              geo%vxY(mpg%fcVx(iFc,2)) )
          wetted_area(i) = wetted_area(i) + geo%fcS(iFc)
        end if
      end do
      extension_r(i) = r_max - r_min
      extension_z(i) = z_max - z_min
      if ( mpg%nnreg(0).eq.8 .and. (i.eq.1 .or. i.eq.2) ) then
        iFc = ifsepimp
      else
        iFc = ifsepomp
      end if
      if (mpg%nStr.gt.0) flux_expansion(i) = &
          & ( geo%fcBb(iFc,0)/geo%fcBb(iFc,3) ) / &
          & ( geo%vxBb(mpg%ivdiv(i),0)/geo%vxBb(mpg%ivdiv(i),3) )
    end do
    if (mpg%nStr.eq.0 .and. mpg%nTgc.ge.1 ) flux_expansion(1) = &
          & ( geo%fcBb(ifsepomp,0)/geo%fcBb(ifsepomp,3) ) / &
          & ( geo%fcBb(mpg%tgVx(1),0)/geo%fcBb(mpg%tgVx(1),3) )

    IDS_initialized = .true.
    return
    end subroutine IDS_init

    !> Process B2.5 data and set it to IMAS IDS.
    !! @note    The \b B25_process_ids routine enables to store data for
    !!          specific time slice. By default it stores single default
    !!          time slice of time slice value 0.0.
    !!          \b num_time_slices_IN is required to beforehand allocate
    !!          required ggd(:) array of nodes structure and for additional
    !!          checks for correct use of the routine.
    !! @note    Time slice value is set as:
    !!          \b time_slice_value = \b time_step_IN * \b time_slice_ind_IN
    subroutine B25_process_ids( &
            &   geo, mpg, state, state_ext, state_avg, switch, &
            &   edge_profiles, edge_sources, edge_transport, &
            &   radiation, description, equilibrium, &
#if ( IMAS_MINOR_VERSION > 21 || IMAS_MAJOR_VERSION > 3 )
            &   summary, &
#endif
#if ( IMAS_MINOR_VERSION > 25 && IMAS_MINOR_VERSION < 34 && IMAS_MAJOR_VERSION == 3 )
            &   numerics, &
#endif
#if ( IMAS_MINOR_VERSION > 21 || IMAS_MAJOR_VERSION > 3 )
            &   run_start_time_IN, run_end_time_IN, &
#endif
#if ( IMAS_MINOR_VERSION > 30 || IMAS_MAJOR_VERSION > 3 )
            &   divertors, &
#endif
            &   time_IN, time_step_IN, shot, run, database, version, &
            &   new_eq_ggd, &
            &   time_slice_ind_IN, num_time_slices_IN )
#ifdef NO_OPT
!DIR$ NOOPTIMIZE
#endif
        use b2mod_dimensions
        implicit none
        type (geometry), intent(in) :: geo
        type (mapping), intent(inout) :: mpg
        type (B2state), intent(inout) :: state
        type (B2stateExt), intent(inout) :: state_ext
        type (B2Average), intent(in) :: state_avg
        type (switches), intent(inout) :: switch
        type (ids_equilibrium) :: equilibrium !< IDS designed to
            !< store equilibrium data
        type (ids_edge_profiles) :: edge_profiles !< IDS designed to
            !< store data on edge plasma profiles (includes the scrape-off
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
#if ( IMAS_MINOR_VERSION > 21 || IMAS_MAJOR_VERSION > 3 )
        type (ids_summary) :: summary !< IDS designed to store
            !< run summary data
#endif
#if ( IMAS_MINOR_VERSION > 25 && IMAS_MINOR_VERSION < 34 && IMAS_MAJOR_VERSION == 3 )
        type (ids_numerics) :: numerics !< IDS designed to store
            !< run numerics data
#endif
#if ( IMAS_MINOR_VERSION > 25 || IMAS_MAJOR_VERSION > 3 )
        real(IDS_real), intent(in) :: run_start_time_IN, run_end_time_IN !< Run time bounds
#endif
#if ( IMAS_MINOR_VERSION > 30 || IMAS_MAJOR_VERSION > 3 )
        type (ids_divertors) :: divertors !< IDS designed to store
            !< data related to the divertor plates
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
        logical, intent(out) :: new_eq_ggd

        !! Internal variables
        character(len=13)  :: spclabel   !< Species label
        integer :: time_sind !< Time slice index
        integer :: ion_charge_int !< Ion charge (e.g. 1, 2, etc.)
        integer :: ns    !< Total number of B2.5 species
        integer :: nsion !< Total number of IDS ion species
        integer :: nneut !< Total number of IDS neutral species
        integer :: n_process !< Number of radiation processes handled
        integer :: is, js, ks !< Species indices (iterators)
        integer :: ii, jj !< Iterators
        integer :: i      !< Iterator
        integer :: j      !< Iterator
        integer :: k      !< Iterator
        integer :: iFc    !< Iterator on faces
        integer :: iCv    !< Iterator on control volumes
        integer :: iCv1, iCv2 !< Indices of the two CVs of both sides
                              !< of the OMP separatrix
        integer :: iVx    !< Iterator on vertices
        integer :: is1    !< First ion of an isonuclear sequence
        integer :: is2    !< Last ion of an isonuclear sequence
        integer :: icnt   !< Boundary cell counter
        integer :: ib     !< Boundary condition index
        integer :: isep(2) !< Array of separatrix regions
        integer, allocatable :: ionstt(:) !< Mapping array
                                          !< from B2-Eirene charged fluids to IDS ion states
        integer, allocatable :: istion(:) !< Number of IDS states for each ion
        integer, allocatable :: ispion(:,:) !< Mapping array
                                            !< from IDS ions and states to B2-Eirene ions
           !< ispion(i,j) contains the B2.5 species index for (ion i,state j) or
           !<                      the Eirene molecular ion index
#ifdef B25_EIRENE
        integer :: ind    !< Non-standard surface index in resolved list
        integer :: ias    !< Starting index for non-standard surface in resolved list
        integer :: iss    !< State index
        integer :: iatm   !< Atom iterator
        integer :: imol   !< Molecule iterator
        integer :: nelems !< Number of elements present in a molecule or molecular ion
        integer :: p      !< Dummy integer
        integer, allocatable :: isstat(:) !< Mapping array
                                          !< from Eirene atoms and molecules to IDS neutral states
        integer, allocatable :: imneut(:) !< Mapping array
                                          !< from Eirene molecules to IDS neutrals
        integer, allocatable :: imiion(:) !< Mapping array
                                          !< from Eirene molecular ions to IDS ion sequences
        integer :: ixx, iyy
#endif
        real(IDS_real) :: flxFace( mpg%nFc, 0:1 )
        real(IDS_real) :: totflux( mpg%nFc, 0:1 )
        real(IDS_real) :: tmpFace( mpg%nFc )
        real(IDS_real) :: tmpCv( mpg%nCv )
        real(IDS_real) :: totCv( mpg%nCv )
        real(IDS_real) :: tmpVx( mpg%nVx )
        real(IDS_real) :: pe( mpg%nCv )
        real(IDS_real) :: zeff( mpg%nCv )
        real(IDS_real) :: time_step !< Time step
        real(IDS_real) :: time_slice_value   !< Time slice value
        real(IDS_real) :: frac, u, v, psi_average,                     &
            &             vtor, nisep, nasum
        real(IDS_real), allocatable :: power_convected(:),             &
            &             power_conducted(:),                          &
            &             power_recombination_plasma(:),               &
            &             power_currents(:),                           &
            &             current_incident(:)
#if ( IMAS_MINOR_VERSION > 36 || IMAS_MAJOR_VERSION > 3 )
        integer :: iFs, iactive, inactive
        real(IDS_real) :: xi, yi, xo, yo
        logical :: found
#endif
#ifdef B25_EIRENE
#ifdef WG_TODO
        real(IDS_real), allocatable :: un0(:,:,:), um0(:,:,:)
#endif
#endif
 !< Type of IDS data structure, designed for handling grid geometry data
#if GGD_MAJOR_VERSION > 0
#if ( IMAS_MINOR_VERSION < 15 && IMAS_MAJOR_VERSION < 4 )
        type(ids_generic_grid_dynamic) :: edge_grid, transport_grid, &
            &  sources_grid
#else
        type(ids_generic_grid_aos3_root) :: edge_grid, transport_grid, &
            &  sources_grid
#if ( IMAS_MINOR_VERSION > 21 || IMAS_MAJOR_VERSION > 3 )
        type(ids_generic_grid_aos3_root) :: radiation_grid
#endif
#endif
        integer :: istrai !< Stratum iterator
#endif

#if ( ( IMAS_MINOR_VERSION > 38 || IMAS_MAJOR_VERSION > 3 ) && defined(B25_EIRENE) )
        integer, parameter :: nsources = 13
#else
        integer, parameter :: nsources = 12
#endif
        integer, save :: ncall = 0
        integer, save :: style = 1
        integer, save :: ismain = 1
        integer, save :: ismain0 = 0
        integer, save :: ids_from_43 = 0
        integer, save :: balance_netcdf = 0
        real(IDS_real), save :: dtim = 1.0_IDS_real
        character*132 radiation_commit
        character*256 filename
        character*5 hlp_frm
        logical match_found, streql, exists, wrong_flow
#ifdef B25_EIRENE
        character(len=132) :: mol_label !< Molecule species label (e.g. D2)
        character(len=132) :: ion_label !< Ion species label (e.g. D+1)
        logical, allocatable :: in_species(:)
#endif
        real (kind=R8) :: intvertex_s
        external b2xpne, b2xpni, b2xppb, b2xppe, b2xppz, b2xzef
        external b2sral, b2tral, b2trql, b2tanml
        external b2xpnn, b2tiner, b2tvspa
        external ipgetr, ipgeti, species, streql, xerrab, xertst
        external find_file, intvertex_s

        !! ===  SET UP IDS ===
        write(0,*) "Setting data for edge_profiles IDS"
        if (ncall.eq.0) then
          call ipgeti ('b2mndt_style', style)
          call ipgeti ('b2mndr_ismain', ismain)
          call ipgeti ('ids_from_43', ids_from_43)
          call ipgetr ('b2mndr_dtim', dtim)
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
            if (is_neutral(is) .and. nint(zn(is)).eq.nint(zn(ismain)) &
                             & .and. nint(am(is)).eq.nint(am(ismain))) then
              ismain0 = is
              match_found = .true.
            end if
            is = is - 1
          end do
          if (.not.match_found.and.ismain.ne.1) ismain0 = ismain
          call ipgeti ('b2mwti_ismain0', ismain0)
        end if
#ifdef B25_EIRENE
        if (ids_from_43.eq.0) then
          source = "SOLPS-ITER_WG"
        else
          source = "SOLPS4.3"
        end if
#else
        if (ids_from_43.eq.0) then
          source = "B2.5_WG"
        else
          source = "B2"
        end if
#endif

        call IDS_init( mpg, geo )
        ns = size( state%pl%na, 2 )
        do is = 0, ns-1
          call b2xppb( mpg%nCv, state%rt%rza(:,is),                 &
            &          state%pl%na(:,is), state%pl%te, state%pl%ti, &
            &          state%pl%tn, is, state%dv%pa(:,is) )
        end do
        call b2xpne( mpg%nCv, ns, state%rt%rza, state%pl%na,        &
            &        state_ext%ne, state%dv%ne)
        call b2xzef( mpg%nCv, ns, state%rt%rz2, state%pl%na,        &
            &        state%dv%ne, zeff, state_ext)
        call b2xppz( mpg%nCv, ns, state%dv%ne, state%pl%na,         &
            &        state%pl%te, state%pl%ti, state%dv%pz, state_ext)
        call b2tanml(mpg%nCv, mpg%nFc, ns, switch, geo, mpg,        &
            &        state%co%csig_an, state%pl%po, state%dv%fchanml)
        call b2tvspa(mpg%nCv, mpg%nFc, mpg%nVx, ns,                 &
            &        switch, geo, mpg, state%pl%ua,                 &
            &        state%co%vsaf_cl, state%dv%fac_vis, state%dv%fchvispar)
        call b2tiner(mpg%nCv, mpg%nFc, mpg%nVx, ns,                 &
            &        switch, geo, mpg, state%pl%na, state%pl%ua,    &
            &        state%dv%facdrift, state%dv%fchinert)

!   ..find nscx, iscx
!     (indices for neutral hydrogen species)
        state%rt%nscx = 0
        do is = 0, ns-1
          if (is_neutral(is).and.nint(zn(is)).eq.1) then
            call xertst (state%rt%nscx.lt.state%rt%nscxmax, &
                & 'too many neutral hydrogen species')
            state%rt%iscx(state%rt%nscx) = is
            state%rt%nscx = state%rt%nscx+1
          endif
        enddo
        do k = state%rt%nscx, state%rt%nscxmax-1
          state%rt%iscx(k) = -1
        enddo
        call b2xpni (mpg%nCv, ns, state%pl%na, state%dv%ni)
        call b2xpnn (mpg%nCv, ns, state%pl%na, state%dv%nn)
        call b2xpne (mpg%nCv, ns, state%rt%rz2, state%pl%na, &
            &        state_ext%ne2, state%dv%ne2)
!   ..compute flux limit coefficients
        call b2trql (mpg%nCv, mpg%nFc, ns, switch, geo, mpg, &
            &        state%pl, state%dv, state_ext,          &
            &        state%co%chvemx, state%co%chvimx)
        call b2tral (mpg%nCv, mpg%nFc, mpg%nVx, ns,          &
            &        state%rt%nscx, state%rt%nscxmax,        &
            &        state%rt%iscx, ismain,                  &
            &        switch, geo, mpg, state%pl, state%dv,   &
            &        state%rt, state_ext, state%co)
!   ..compute sources
        call b2sral ( mpg%nCv, mpg%nFc, mpg%nVx, ns,         &
     &   state%rt%nscx, state%rt%nscxmax, state%rt%iscx,     &
     &   ismain, ismain0, dtim,                              &
     &   switch, geo, mpg, state, state_ext, state_avg, wrong_flow, .false.)

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
        !! This routine only fills in one time slice at a time
        time_sind = 1
        slice_index = time_sind
        time_slice_value = time
        time_step = IDS_REAL_INVALID
        num_time_slices = 1
        num_slices = num_time_slices
        !! If present, set time step values
        if( present( time_step_IN ) ) time_step = time_step_IN
!        if( present( time_slice_ind_IN ) ) time_sind = time_slice_ind_IN
!        if( present( num_time_slices_IN ) ) num_time_slices = num_time_slices_IN
        !! Check if num_time_slices >= time_sind
        call xertst( num_time_slices .ge. time_sind, &
            & "B25_process_ids: Time step index cannot be greater " // &
            & "than total number of time steps!" )
        if( present( time_slice_ind_IN ) ) &
            & call xertst( time_slice_ind_IN .ge. 1, &
            & "faulty argument time_slice_ind_IN" )
        if( present( num_time_slices_IN ) ) &
            & call xertst( num_time_slices_IN .ge. 1, &
            & "faulty argument num_time_slices_IN" )

        !! Preparing IDSs for writing
        !! In order to write to IDS database there are next steps that are
        !! mandatory to do, otherwise there is high chance that writing to IDS
        !! database will fail
        comment = label
        !! 1. Set homogeneous_time to 0 or 1 and other properties
        call write_ids_properties( edge_profiles%ids_properties, &
          &  homogeneous_time )
        call write_ids_properties( edge_transport%ids_properties, &
          &  homogeneous_time )
        call write_ids_properties( edge_sources%ids_properties, &
          &  homogeneous_time )
        call write_ids_properties( radiation%ids_properties, &
          &  homogeneous_time )
        call write_ids_properties( description%ids_properties, 2 )
#if ( IMAS_MINOR_VERSION > 21 || IMAS_MAJOR_VERSION > 3 )
        call write_ids_properties( summary%ids_properties, &
          &  homogeneous_time )
#endif
#if ( IMAS_MINOR_VERSION > 25 && IMAS_MINOR_VERSION < 34 && IMAS_MAJOR_VERSION == 3 )
        call write_ids_properties( numerics%ids_properties, &
          &  homogeneous_time )
#endif
#if ( IMAS_MINOR_VERSION > 30 || IMAS_MAJOR_VERSION > 3 )
        call write_ids_properties( divertors%ids_properties, &
          &  homogeneous_time )
#endif

        !! 2. Set code and library data
#ifdef B25_EIRENE
        code_description = &
         & "Snapshot IDS from b2mod_ual_io routine (coupled SOLPS-ITER run)"
#else
        code_description = &
         & "Snapshot IDS from b2mod_ual_io routine (standalone B2.5 run)"
#endif
        radiation_commit = B25_git_version
        call write_ids_code( switch, edge_profiles%code, code_commit, code_description )
        call write_ids_code( switch, edge_transport%code, code_commit, code_description )
        call write_ids_code( switch, edge_sources%code, code_commit, code_description )
        call write_ids_code( switch, radiation%code, radiation_commit, code_description )
#if ( IMAS_MINOR_VERSION > 21 || IMAS_MAJOR_VERSION > 3 )
        call write_ids_code( switch, summary%code, code_commit, code_description )
#endif
#if ( IMAS_MINOR_VERSION > 30 || IMAS_MAJOR_VERSION > 3 )
        call write_ids_code( switch, divertors%code, code_commit, code_description )
#endif
#if IMAS_MAJOR_VERSION > 3
        call write_ids_code_constant( description%code, code_commit, code_description )
#endif
        allocate( edge_transport%model(1) )
        allocate( edge_transport%model(1)%identifier%name(1) )
        allocate( edge_transport%model(1)%identifier%description(1) )
        if (ids_from_43.eq.0) then
          if (style.eq.0) then
            edge_transport%model(1)%identifier%index = -2
            edge_transport%model(1)%identifier%name(1) = "SOLPS5.0"
            edge_transport%model(1)%identifier%description(1) = "SOLPS5.0 physics model"
          else if (style.ge.1) then
            edge_transport%model(1)%identifier%index = -3
            edge_transport%model(1)%identifier%name(1) = "SOLPS5.2"
            edge_transport%model(1)%identifier%description(1) = "SOLPS5.2 physics model"
          else if (style.eq.-1) then
            edge_transport%model(1)%identifier%index = -1
            edge_transport%model(1)%identifier%name(1) = "SOLPS4.3"
            edge_transport%model(1)%identifier%description(1) = "SOLPS4.3 physics model"
          end if
          edge_transport%model(1)%flux_multiplier = 1.5_IDS_real + switch%BoRiS
        else
          edge_transport%model(1)%identifier%index = -1
          edge_transport%model(1)%identifier%name(1) = "SOLPS4.3"
          edge_transport%model(1)%identifier%description(1) = "SOLPS4.3 physics model"
          edge_transport%model(1)%flux_multiplier = 2.5_IDS_real
        end if
#if ( IMAS_MINOR_VERSION > 29 || IMAS_MAJOR_VERSION > 3 )
        allocate( edge_transport%model(1)%code%name(1) )
        edge_transport%model(1)%code%name = source
#if ( IMAS_MINOR_VERSION > 38 || IMAS_MAJOR_VERSION > 3 )
        allocate( edge_transport%model(1)%code%description(1) )
        edge_transport%model(1)%code%description = &
            & "Snapshot IDS written by b2mod_ual_io routine"
#endif
        allocate( edge_transport%model(1)%code%version(1) )
        edge_transport%model(1)%code%version = newversion
        allocate( edge_transport%model(1)%code%commit(1) )
        edge_transport%model(1)%code%commit = B25_git_version
        allocate( edge_transport%model(1)%code%repository(1) )
        edge_transport%model(1)%code%repository = "ssh://git.iter.org/bnd/b2.5.git"
        call write_timed_integer( edge_transport%model(1)%code%output_flag, 0 )
#endif

        !! 3. Allocate IDS.time and set it to desired values
        allocate( edge_profiles%time(num_time_slices) )
        edge_profiles%time(time_sind) = time
        allocate( edge_transport%time(num_time_slices) )
        edge_transport%time(time_sind) = time
        allocate( edge_sources%time(num_time_slices) )
        edge_sources%time(time_sind) = time
#if IMAS_MAJOR_VERSION < 4
        allocate( description%time(num_time_slices) )
        description%time(time_sind) = time
#endif
#if ( IMAS_MINOR_VERSION > 21 || IMAS_MAJOR_VERSION > 3 )
        allocate( summary%time(num_time_slices) )
        summary%time(time_sind) = time
#endif
#if ( IMAS_MINOR_VERSION > 25 && IMAS_MINOR_VERSION < 34 && IMAS_MAJOR_VERSION == 3 )
        !! Allocate numerics.time and set it to desired values
        allocate( numerics%time(num_time_slices) )
        numerics%time(time_sind) = time
        allocate( numerics%time_start(num_time_slices) )
        numerics%time_start(time_sind) = run_start_time_IN
        allocate( numerics%time_step(num_time_slices) )
        numerics%time_step(time_sind) = time_step
        allocate( numerics%time_end(num_time_slices) )
        numerics%time_end(time_sind) = run_end_time_IN
#endif
#if ( IMAS_MINOR_VERSION > 30 || IMAS_MAJOR_VERSION > 3 )
        allocate( divertors%time(num_time_slices) )
        divertors%time(time_sind) = time
#endif
        allocate( radiation%time(num_time_slices) )
        radiation%time(time_sind) = time

        !! Allocate radiation.process
        !! Process 1: line and recombination radiation due to B2.5 species
        !! Process 2: bremsstrahlung recombination due to B2.5 species
        !! Process 3: line radiation due to Eirene neutrals (atoms and molecules)
        !! Process 4: line radiation dur to Eirene molecular ions
        if (switch%use_eirene.ne.0) then
          n_process = 4
        else
          n_process = 2
        end if
        allocate( radiation%process(n_process) )
        do j = 1, n_process
          allocate( radiation%process(j)%identifier%name(1) )
          allocate( radiation%process(j)%identifier%description(1) )
        end do
        radiation%process(1)%identifier%index = 2
        radiation%process(1)%identifier%name = 'line_radiation'
        radiation%process(1)%identifier%description = 'Line and rec. rad. from B2.5 species'
#if ( IMAS_MAJOR_VERSION > 3 || IMAS_MINOR_VERSION > 30 )
        radiation%process(2)%identifier%index = &
          &  radiation_identifier%bremsstrahlung
        radiation%process(2)%identifier%name = &
          &  radiation_identifier%name( radiation_identifier%bremsstrahlung )
        radiation%process(2)%identifier%description = &
          &  radiation_identifier%description( radiation_identifier%bremsstrahlung )
#else
        radiation%process(2)%identifier%index = 8
        radiation%process(2)%identifier%name = 'bremsstrahlung'
        radiation%process(2)%identifier%description = 'Bremsstrahlung from B2.5 species'
#endif
        if (switch%use_eirene.ne.0) then
          radiation%process(3)%identifier%index = 1
          radiation%process(4)%identifier%index = 2
          radiation%process(3)%identifier%name = 'line_radiation'
          radiation%process(3)%identifier%description = 'Line radiation from Eirene neutrals'
          radiation%process(4)%identifier%name = 'line_radiation'
          radiation%process(4)%identifier%description = 'Line radiation from Eirene mol. ions'
        end if

        !! Allocate ggd for number of different time steps
        allocate( edge_profiles%ggd( num_time_slices ) )
#if ( IMAS_MINOR_VERSION > 14 || IMAS_MAJOR_VERSION > 3 )
        allocate( edge_profiles%grid_ggd( num_time_slices ) )
        allocate( edge_transport%grid_ggd( num_time_slices ) )
        allocate( edge_sources%grid_ggd( num_time_slices ) )
#if ( IMAS_MINOR_VERSION > 21 || IMAS_MAJOR_VERSION > 3 )
        allocate( radiation%grid_ggd( num_time_slices ) )
#endif
#endif
        allocate( edge_transport%model(1)%ggd( num_time_slices ) )
        allocate( edge_sources%source(nsources) )
        do is = 1, nsources
          allocate( edge_sources%source(is)%ggd( num_time_slices ) )
          allocate( edge_sources%source(is)%identifier%name(1) )
          allocate( edge_sources%source(is)%identifier%description(1) )
        end do

#if ( IMAS_MAJOR_VERSION > 3 || IMAS_MINOR_VERSION > 30 )
        !! Total sources
        edge_sources%source(1)%identifier%index = edge_source_identifier%total
        edge_sources%source(1)%identifier%name = &
          &  edge_source_identifier%name( edge_source_identifier%total )
        edge_sources%source(1)%identifier%description = &
          &  edge_source_identifier%description( edge_source_identifier%total )
        !! Background sources
        edge_sources%source(2)%identifier%index = edge_source_identifier%background
        edge_sources%source(2)%identifier%name = &
          &  edge_source_identifier%name( edge_source_identifier%background )
        edge_sources%source(2)%identifier%description = &
          &  edge_source_identifier%description( edge_source_identifier%background )
        !! Prescribed sources
        edge_sources%source(3)%identifier%index = edge_source_identifier%prescribed
        edge_sources%source(3)%identifier%name = &
          &  edge_source_identifier%name( edge_source_identifier%prescribed )
        edge_sources%source(3)%identifier%description = &
          &  edge_source_identifier%description( edge_source_identifier%prescribed )
        !! Time derivatives
        edge_sources%source(4)%identifier%index = edge_source_identifier%time_derivative
        edge_sources%source(4)%identifier%name = &
          &  edge_source_identifier%name( edge_source_identifier%time_derivative )
        edge_sources%source(4)%identifier%description = &
          &  edge_source_identifier%description( edge_source_identifier%time_derivative )
        !! Atomic ionization
        edge_sources%source(5)%identifier%index = edge_source_identifier%atomic_ionization
        edge_sources%source(5)%identifier%name = &
          &  edge_source_identifier%name( edge_source_identifier%atomic_ionization )
        edge_sources%source(5)%identifier%description = &
          &  edge_source_identifier%description( edge_source_identifier%atomic_ionization )
        !! Molecular ionization
        edge_sources%source(6)%identifier%index = edge_source_identifier%molecular_ionization
        edge_sources%source(6)%identifier%name = &
          &  edge_source_identifier%name( edge_source_identifier%molecular_ionization )
        edge_sources%source(6)%identifier%description = &
          &  edge_source_identifier%description( edge_source_identifier%molecular_ionization )
        !! Ionization
        edge_sources%source(7)%identifier%index = edge_source_identifier%ionization
        edge_sources%source(7)%identifier%name = &
          &  edge_source_identifier%name( edge_source_identifier%ionization )
        edge_sources%source(7)%identifier%description = &
          &  edge_source_identifier%description( edge_source_identifier%ionization )
        !! Recombination
        edge_sources%source(8)%identifier%index = edge_source_identifier%recombination
        edge_sources%source(8)%identifier%name = &
          &  edge_source_identifier%name( edge_source_identifier%recombination )
        edge_sources%source(8)%identifier%description = &
          &  edge_source_identifier%description( edge_source_identifier%recombination )
        !! Charge exchange
        edge_sources%source(9)%identifier%index = edge_source_identifier%charge_exchange
        edge_sources%source(9)%identifier%name = &
          &  edge_source_identifier%name( edge_source_identifier%charge_exchange )
        edge_sources%source(9)%identifier%description = &
          &  edge_source_identifier%description( edge_source_identifier%charge_exchange )
        !! Collisional equipartition
        edge_sources%source(10)%identifier%index = edge_source_identifier%collisional_equipartition
        edge_sources%source(10)%identifier%name = &
          &  edge_source_identifier%name( edge_source_identifier%collisional_equipartition )
        edge_sources%source(10)%identifier%description = &
          &  edge_source_identifier%description( edge_source_identifier%collisional_equipartition )
        !! Ohmic
        edge_sources%source(11)%identifier%index = edge_source_identifier%ohmic
        edge_sources%source(11)%identifier%name = &
          &  edge_source_identifier%name( edge_source_identifier%ohmic )
        edge_sources%source(11)%identifier%description = &
          &  edge_source_identifier%description( edge_source_identifier%ohmic )
        !! Radiation
        edge_sources%source(12)%identifier%index = edge_source_identifier%radiation
        edge_sources%source(12)%identifier%name = &
          &  edge_source_identifier%name( edge_source_identifier%radiation )
        edge_sources%source(12)%identifier%description = &
          &  edge_source_identifier%description( edge_source_identifier%radiation )
#if ( ( IMAS_MINOR_VERSION > 38 || IMAS_MAJOR_VERSION > 3 ) && defined(B25_EIRENE) )
        !! Neutrals
        edge_sources%source(13)%identifier%index = edge_source_identifier%neutrals
        edge_sources%source(13)%identifier%name = &
          &  edge_source_identifier%name( edge_source_identifier%neutrals )
        edge_sources%source(13)%identifier%description = &
          &  edge_source_identifier%description( edge_source_identifier%neutrals )
#endif
#else
        !! Total sources
        edge_sources%source(1)%identifier%index = 1
        edge_sources%source(1)%identifier%name = "Total"
        edge_sources%source(1)%identifier%description = "Total source from "//trim(source)
        !! Background sources
        edge_sources%source(2)%identifier%index = 703
        edge_sources%source(2)%identifier%name = "Background"
        edge_sources%source(2)%identifier%description = "External sources from "//trim(source)
        !! Prescribed sources
        edge_sources%source(3)%identifier%index = 705
        edge_sources%source(3)%identifier%name = "Prescribed"
        edge_sources%source(3)%identifier%description = &
            & "Boundary conditions sources from "//trim(source)
        !! Time derivatives
        edge_sources%source(4)%identifier%index = 706
        edge_sources%source(4)%identifier%name = "Time derivative"
        edge_sources%source(4)%identifier%description = &
            & "Time derivative sources from "//trim(source)
        !! Atomic ionization
        edge_sources%source(5)%identifier%index = 707
        edge_sources%source(5)%identifier%name = "Atomic ionization"
        edge_sources%source(5)%identifier%description = &
            & "Atomic ionization sources from "//trim(source)
        !! Molecular ionization
        edge_sources%source(6)%identifier%index = 708
        edge_sources%source(6)%identifier%name = "Molecular ionization"
        edge_sources%source(6)%identifier%description = &
            & "Molecular ionization sources from "//trim(source)
        !! Ionization
        edge_sources%source(7)%identifier%index = 709
        edge_sources%source(7)%identifier%name = "Ionization"
        edge_sources%source(7)%identifier%description = &
            & "Ionization sources from "//trim(source)
        !! Recombination
        edge_sources%source(8)%identifier%index = 710
        edge_sources%source(8)%identifier%name = "Recombination"
        edge_sources%source(8)%identifier%description = &
            & "Recombination sources from "//trim(source)
        !! Charge exchange
        edge_sources%source(9)%identifier%index = 305
        edge_sources%source(9)%identifier%name = "Charge exchange"
        edge_sources%source(9)%identifier%description = &
            & "Charge exchange sources from "//trim(source)
        !! Collisional equipartition
        edge_sources%source(10)%identifier%index = 11
        edge_sources%source(10)%identifier%name = "Equipartition"
        edge_sources%source(10)%identifier%description = &
            & "Collisional equipartition sources from "//trim(source)
        !! Ohmic
        edge_sources%source(11)%identifier%index = 7
        edge_sources%source(11)%identifier%name = "Ohmic"
        edge_sources%source(11)%identifier%description = &
            & "Ohmic (Joule) sources from "//trim(source)
        !! Radiation
        edge_sources%source(12)%identifier%index = 200
        edge_sources%source(12)%identifier%name = "Radiation"
        edge_sources%source(12)%identifier%description = &
            & "Radiation sources from "//trim(source)
#if ( ( IMAS_MINOR_VERSION > 38 || IMAS_MAJOR_VERSION > 3 ) && defined(B25_EIRENE) )
        !! Neutrals
        edge_sources%source(13)%identifier%index = 701
        edge_sources%source(13)%identifier%name = "Neutrals"
        edge_sources%source(13)%identifier%description = &
            & "Total source due to plasma-neutral interactions from "//trim(source)
#endif
#endif

        call put_equilibrium_data ( mpg, geo, equilibrium, &
#if ( IMAS_MINOR_VERSION > 21 || IMAS_MAJOR_VERSION > 3 )
            &  summary, &
#endif
            &  edge_profiles, database, &
#if ( AL_MAJOR_VERSION > 4 && GGD_MAJOR_VERSION > 0 )
            &  time_sind, &
#endif
            &  time_slice_value, .true., new_eq_ggd )
        allocate( radiation%vacuum_toroidal_field%b0( num_time_slices ) )
        radiation%vacuum_toroidal_field%b0( time_sind ) = &
            &  edge_profiles%vacuum_toroidal_field%b0( time_sind )
        radiation%vacuum_toroidal_field%r0 = &
            &  edge_profiles%vacuum_toroidal_field%r0
#if ( IMAS_MINOR_VERSION > 21 && IMAS_MAJOR_VERSION == 3 )
        allocate( description%data_entry%user(1) )
        description%data_entry%user = username
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
#endif
#if ( IMAS_MINOR_VERSION > 21 || IMAS_MAJOR_VERSION > 3 )
        if ( present( time_step_IN ) ) &
          &  description%simulation%time_step = time_step_IN
        if ( present ( time_IN ) ) &
          &  description%simulation%time_current = time_IN
        allocate( description%simulation%workflow(1) )
        description%simulation%workflow = source
#if ( IMAS_MINOR_VERSION > 25 || IMAS_MAJOR_VERSION > 3 )
        description%simulation%time_begin = run_start_time_IN
        description%simulation%time_end = run_end_time_IN
#endif
#if IMAS_MAJOR_VERSION > 3
        description%type%index = 2
        allocate( description%type%name(1) )
        allocate( description%type%description(1) )
        description%type%name = "simulation"
        description%type%description = "Simulation results from "//trim(source)
        allocate( description%machine(1) )
        description%machine = database
        description%pulse = shot
#endif

        i=index(B25_git_version,'-')
        if (i.gt.0) then
          allocate( summary%tag%name(1) )
          write(hlp_frm,'(a,i2.2,a)') '(a',i-1,')'
          write(summary%tag%name,hlp_frm) B25_git_version(1:i-1)
        endif

#if ( IMAS_MINOR_VERSION > 32 || IMAS_MAJOR_VERSION > 3 )
        call write_ids_midplane( divertors%midplane, midplane_id )
        call write_ids_midplane( edge_profiles%midplane, midplane_id )
        call write_ids_midplane( edge_sources%midplane, midplane_id )
        call write_ids_midplane( edge_transport%midplane, midplane_id )
        call write_ids_midplane( summary%midplane, midplane_id )
#endif

        if (mpg%iFssep.ne.US_GRID_UNDEFINED) then
          u = 0.0_IDS_real
          do i = mpg%fsFcP(mpg%iFssep,1), &
               & mpg%fsFcP(mpg%iFssep,1) + mpg%fsFcP(mpg%iFssep,2) - 1
            iCv1 = mpg%fcCv(mpg%fsFc(i),1)
            iCv2 = mpg%fcCv(mpg%fsFc(i),2)
            if (mpg%cvOnClosedSurface(iCv1).neqv. &
              & mpg%cvOnClosedSurface(iCv2)) u = u + state%dv%fht(mpg%fsFc(i),1)
          end do
          if (mpg%iFssep2.ne.US_GRID_UNDEFINED) then
            do i = mpg%fsFcP(mpg%iFssep2,1), &
                 & mpg%fsFcP(mpg%iFssep2,1) + mpg%fsFcP(mpg%iFssep2,2) - 1
              iCv1 = mpg%fcCv(mpg%fsFc(i),1)
              iCv2 = mpg%fcCv(mpg%fsFc(i),2)
              if (mpg%cvOnClosedSurface(iCv1).neqv. &
                & mpg%cvOnClosedSurface(iCv2)) u = u + state%dv%fht(mpg%fsFc(i),1)
            end do
          end if
          if (u.ne.0.0_IDS_real) then
#if ( IMAS_MINOR_VERSION > 28 || IMAS_MAJOR_VERSION > 3 )
            call write_sourced_value( summary%global_quantities%power_loss, u )
#endif
          end if
        end if

        select case (GeometryType)
        case ( GEOMETRY_LINEAR, GEOMETRY_CYLINDER )
          u = 0.0_IDS_real
          frac = 0.0_IDS_real
          do iCv = 1, mpg%nCi
            if (mpg%iFssep.ne.US_GRID_UNDEFINED) then
              psi_average = 0.0_R8
              do i = mpg%cvVxP(iCv,1), mpg%cvVxP(iCv,1) + mpg%cvVxP(iCv,2) - 1
                iVx = mpg%cvVx(i)
                psi_average = psi_average + geo%vxFpsi(iVx)
              end do
              psi_average = psi_average / real(mpg%cvVxP(iCv,2))
              match_found = &
                & (geo%psi_increasing.and.psi_average.lt.geo%fsPsi(mpg%iFssep)).or. &
                & (.not.geo%psi_increasing.and.psi_average.gt.geo%fsPsi(mpg%iFssep))
            else
              match_found = .false.
            end if
            do is = 0, ns-1
              u = u + state%srw%rqrad(iCv,is) + state%srw%rqbrm(iCv,is)
              if (match_found) frac = frac + state%srw%rqrad(iCv,is) + state%srw%rqbrm(iCv,is)
            end do
#ifdef B25_EIRENE
            do is = 1, natmi
              u = u - eneutrad(iCv,is,0)
              if (match_found) frac = frac - eneutrad(iCv,is,0)
            end do
            do is = 1, nmoli
              u = u - emolrad(iCv,is,0)
              if (match_found) frac = frac - emolrad(iCv,is,0)
            end do
            do is = 1, nioni
              u = u - eionrad(iCv,is,0)
              if (match_found) frac = frac - eionrad(iCv,is,0)
            end do
#endif
          end do
          if (u.ne.0.0_IDS_real) then
            call write_sourced_value( summary%global_quantities%power_radiated, u )
          end if
#if ( IMAS_MINOR_VERSION > 30 || IMAS_MAJOR_VERSION > 3 )
          if (frac.ne.0.0_IDS_real) then
            call write_sourced_value( summary%global_quantities%power_radiated_inside_lcfs, frac )
            call write_sourced_value( summary%global_quantities%power_radiated_outside_lcfs, &
              &  u - frac )
          end if
#endif
        case ( GEOMETRY_LIMITER, GEOMETRY_SN, &
            &  GEOMETRY_STELLARATORISLAND, GEOMETRY_ANNULUS , &
            &  GEOMETRY_CDN, GEOMETRY_DDN_BOTTOM, GEOMETRY_DDN_TOP, &
            &  GEOMETRY_LFS_SNOWFLAKE_MINUS, GEOMETRY_LFS_SNOWFLAKE_PLUS)
          u = 0.0_IDS_real
          do iCv = 1, mpg%nCi
            if (mpg%cvOnClosedSurface(iCv).and.mpg%nnreg(0).ne.5) cycle
            if (mpg%nnreg(0).eq.5.and.mpg%cvReg(iCv).eq.1) cycle
            do is = 0, ns-1
              u = u + state%srw%rqrad(iCv,is) + state%srw%rqbrm(iCv,is)
            end do
#ifdef B25_EIRENE
            do is = 1, natmi
              u = u - eneutrad(iCv,is,0)
            end do
            do is = 1, nmoli
              u = u - emolrad(iCv,is,0)
            end do
            do is = 1, nioni
              u = u - eionrad(iCv,is,0)
            end do
#endif
          end do
#if ( IMAS_MINOR_VERSION > 30 || IMAS_MAJOR_VERSION > 3 )
          if (u.ne.0.0_IDS_real) then
            call write_sourced_value( summary%global_quantities%power_radiated_outside_lcfs, u )
          end if
#endif
        end select

#if ( IMAS_MINOR_VERSION > 36 || IMAS_MAJOR_VERSION > 3 )
        iactive = 0
        do i = 1, mpg%nXpt
          if (mpg%vxFs(mpg%Xpt(i)).eq.mpg%iFssep .or. &
           & (mpg%vxFs(mpg%Xpt(i)).eq.mpg%iFssep2.and.mpg%iFssep2.gt.0)) then
            if (iactive.eq.0) then
              iactive = i
            else if (geo%vxY(mpg%Xpt(iactive)).lt.geo%vxY(mpg%Xpt(i))) then
              iactive = i
            end if
          end if
        end do
        if (iactive.gt.0) &
            &  call write_sourced_rz( summary%boundary%x_point_main, &
            &   geo%vxX(mpg%Xpt(iactive)), geo%vxY(mpg%Xpt(iactive)) )
        select case (GeometryType)
        case ( GEOMETRY_CDN )
          call write_sourced_value( summary%boundary%distance_inner_outer_separatrices, 0.0_IDS_real )
        case ( GEOMETRY_DDN_BOTTOM, GEOMETRY_DDN_TOP, &
               GEOMETRY_LFS_SNOWFLAKE_MINUS )
          if (nomp.gt.0) then
            inactive = 0
            do i = 1, mpg%nXpt
              if (i.eq.iactive) cycle
              if (mpg%vxFs(mpg%Xpt(i)).eq.mpg%iFssep) cycle
              iFs = mpg%vxFs(mpg%Xpt(i))
              if (iFs.gt.0) then
                if (inactive.eq.0) then
                  inactive = i
                else if (iFs.lt.mpg%vxFs(mpg%Xpt(inactive))) then
                  inactive = i
                end if
              end if
            end do
            if (inactive.gt.0) then
              do i = mpg%cvFcP(omp(icsepomp),1), &
                   & mpg%cvFcP(omp(icsepomp),1) + mpg%cvFcP(omp(icsepomp),2) - 1
                iFs = mpg%fcFs(mpg%cvFc(i))
                if (iFs.eq.mpg%iFssep) then
                  xi = (geo%vxX(mpg%fcVx(mpg%cvFc(i),1)) + &
                     &  geo%vxX(mpg%fcVx(mpg%cvFc(i),2)))/2.0_R8
                  yi = (geo%vxY(mpg%fcVx(mpg%cvFc(i),1)) + &
                     &  geo%vxY(mpg%fcVx(mpg%cvFc(i),2)))/2.0_R8
                end if
              end do
              found = .false.
              do i = icsepomp, nomp
                if (found) cycle
                do j = mpg%cvFcP(omp(i),1), &
                     & mpg%cvFcP(omp(i),1) + mpg%cvFcP(omp(i),2) - 1
                  iFs = mpg%fcFs(mpg%cvFc(j))
                  if (iFs.eq.mpg%vxFs(mpg%Xpt(inactive))) then
                    found = .true.
                    xo = (geo%vxX(mpg%fcVx(mpg%cvFc(j),1)) + &
                       &  geo%vxX(mpg%fcVx(mpg%cvFc(j),2)))/2.0_R8
                    yo = (geo%vxY(mpg%fcVx(mpg%cvFc(j),1)) + &
                       &  geo%vxY(mpg%fcVx(mpg%cvFc(j),2)))/2.0_R8
                  end if
                end do
              end do
              if (found) then
                u = norm ( xo-xi, yo-yi )
                call write_sourced_value( &
                  &  summary%boundary%distance_inner_outer_separatrices, u )
              end if
            end if
          end if
        end select
#endif

#if ( IMAS_MINOR_VERSION > 30 || IMAS_MAJOR_VERSION > 3 )
        if ( mpg%nStr.gt.0 ) then
          allocate( power_currents( maxval(mpg%strDiv) ) )
          allocate( power_convected( maxval(mpg%strDiv) ) )
          allocate( power_conducted( maxval(mpg%strDiv) ) )
          allocate( power_recombination_plasma( maxval(mpg%strDiv) ) )
          allocate( current_incident( maxval(mpg%strDiv) ) )
          do i = 1, maxval(mpg%strDiv)
            power_conducted(i) = 0.0_IDS_real
            power_convected(i) = 0.0_IDS_real
            power_recombination_plasma(i) = 0.0_IDS_real
            power_currents(i) = 0.0_IDS_real
            current_incident(i) = 0.0_IDS_real
            do j = mpg%divFcP(i,1), mpg%divFcP(i,1) + mpg%divFcP(i,2) - 1
              iFc = mpg%divFc(j)
              u = 0.0_R8
              do is = 0, ns-1
                if (is_neutral(is).and.nint(zn(is)).eq.1) then
                  u = u + mpg%divFcOr(j)*state%pl%tn(iCv)* &
                    & (1.5_R8*(state%dv%fna_32(iFc,0,is) + &
                    &          state%dv%fna_32(iFc,1,is) ) + &
                    &  2.5_R8*(state%dv%fna_52(iFc,0,is) + &
                    &          state%dv%fna_52(iFc,1,is) ) )
                else
                  u = u + mpg%divFcOr(j)*  &
                    & (state%pl%ti(iCv) + &
                    &  state%pl%te(iCv)*state%rt%rza(iCv,is))* &
                    & (1.5_R8*(state%dv%fna_32(iFc,0,is) + &
                    &          state%dv%fna_32(iFc,1,is) ) + &
                    &  2.5_R8*(state%dv%fna_52(iFc,0,is) + &
                    &          state%dv%fna_52(iFc,1,is) ) )
                  power_recombination_plasma(i) = &
                    &  power_recombination_plasma(i) + mpg%divFcOr(j)* &
                    & (state%dv%fhp(iFc,0,is) + state%dv%fhp(iFc,1,is) )
                end if
              end do
              power_convected(i) = power_convected(i) + u
              power_conducted(i) = power_conducted(i) - u + mpg%divFcOr(j) *  &
                  & ( state%dv%fht(iFc,0) - state%dv%fhj(iFc,0) +  &
                  &   state%dv%fht(iFc,1) - state%dv%fhj(iFc,1) )
              power_currents(i) = power_currents(i) + mpg%divFcOr(j) * &
                  & ( state%dv%fhj(iFc,0) + state%dv%fhj(iFc,1) )
              current_incident(i) = current_incident(i) + mpg%divFcOr(j) * &
                  & ( state%dv%fch(iFc,0) + state%dv%fch(iFc,1) )
            end do
          end do
        end if
        select case ( GeometryType )
        case ( GEOMETRY_LINEAR, GEOMETRY_CYLINDER )
          if (mpg%nStr.gt.0) then
            allocate( divertors%divertor( mpg%nStr ) )
            do j = 1, mpg%nStr
              i = mpg%strDiv(j)
              allocate( divertors%divertor(i)%name(1) )
              allocate( divertors%divertor(i)%target(1) )
              allocate( divertors%divertor(i)%target(1)%name(1) )
              if (streql(plate_name(i),'W')) then
                divertors%divertor(i)%name = 'Western divertor'
                divertors%divertor(i)%target(1)%name = 'Western target'
              else if (streql(plate_name(i),'E')) then
                divertors%divertor(i)%name = 'Eastern divertor'
                divertors%divertor(i)%target(1)%name = 'Eastern target'
              else
                divertors%divertor(i)%name = 'Divertor '//int2str(i)
                divertors%divertor(i)%target(1)%name = 'Target '//int2str(i)
              end if
#if IMAS_MAJOR_VERSION > 3
              allocate( divertors%divertor(i)%description(1) )
              divertors%divertor(i)%description = plate_name(i)
              allocate( divertors%divertor(i)%target(1)%description(1) )
              divertors%divertor(i)%target(1)%description = plate_name(i)
#else
              allocate( divertors%divertor(i)%identifier(1) )
              divertors%divertor(i)%identifier = plate_name(i)
              allocate( divertors%divertor(i)%target(1)%identifier(1) )
              divertors%divertor(i)%target(1)%identifier = plate_name(i)
#endif
              divertors%divertor(i)%target(1)%extension_r = extension_r(i)
              divertors%divertor(i)%target(1)%extension_z = extension_z(i)
              call write_timed_value( &
                &  divertors%divertor(i)%target(1)%flux_expansion, &
                &  flux_expansion(i) )
              call write_timed_value( &
                &  divertors%divertor(i)%target(1)%wetted_area, &
                &  wetted_area(i) )
              call write_timed_value( &
                &  divertors%divertor(i)%wetted_area, &
                &  wetted_area(i) )
              call write_timed_value( &
                &  divertors%divertor(i)%target(1)%power_incident_fraction, &
                &  1.0_IDS_real )
              call write_timed_value( &
                &  divertors%divertor(i)%target(1)%power_conducted, &
                &  power_conducted(i) )
              call write_timed_value( &
                &  divertors%divertor(i)%power_conducted, &
                &  power_conducted(i) )
              call write_timed_value( &
                &  divertors%divertor(i)%target(1)%power_convected, &
                &  power_convected(i) )
              call write_timed_value( &
                &  divertors%divertor(i)%power_convected, &
                &  power_convected(i) )
              call write_timed_value( &
                &  divertors%divertor(i)%target(1)%power_recombination_plasma, &
                &  power_recombination_plasma(i) )
              call write_timed_value( &
                &  divertors%divertor(i)%power_recombination_plasma, &
                &  power_recombination_plasma(i) )
              call write_timed_value( &
                &  divertors%divertor(i)%target(1)%power_currents, &
                &  power_currents(i) )
              call write_timed_value( &
                &  divertors%divertor(i)%power_currents, &
                &  power_currents(i) )
#if ( IMAS_MINOR_VERSION > 32 || IMAS_MAJOR_VERSION > 3 )
              call write_timed_value( &
                &  divertors%divertor(i)%target(1)%current_incident, &
                &  current_incident(i) )
              call write_timed_value( &
                &  divertors%divertor(i)%current_incident, &
                &  current_incident(i) )
#endif
            end do
          end if
        case ( GEOMETRY_SN, GEOMETRY_STELLARATORISLAND )
          allocate( divertors%divertor(1) )
          allocate( divertors%divertor(1)%name(1) )
          if (geo%LSN) then
            divertors%divertor(1)%name = 'Lower divertor'
          else
            divertors%divertor(1)%name = 'Upper divertor'
          end if
          allocate( divertors%divertor(1)%target(2) )
          allocate( divertors%divertor(1)%target(1)%name(1) )
          allocate( divertors%divertor(1)%target(2)%name(1) )
#if IMAS_MAJOR_VERSION > 3
          allocate( divertors%divertor(1)%description(1) )
          allocate( divertors%divertor(1)%target(1)%description(1) )
          allocate( divertors%divertor(1)%target(2)%description(1) )
          if (geo%LSN) then
            divertors%divertor(1)%description = 'LSN'
            divertors%divertor(1)%target(1)%description = "ID"
            divertors%divertor(1)%target(2)%description = "OD"
          else
            divertors%divertor(1)%description = 'USN'
            divertors%divertor(1)%target(1)%description = "OD"
            divertors%divertor(1)%target(2)%description = "ID"
          end if
#else
          allocate( divertors%divertor(1)%identifier(1) )
          allocate( divertors%divertor(1)%target(1)%identifier(1) )
          allocate( divertors%divertor(1)%target(2)%identifier(1) )
          if (geo%LSN) then
            divertors%divertor(1)%identifier = 'LSN'
            divertors%divertor(1)%target(1)%identifier = "ID"
            divertors%divertor(1)%target(2)%identifier = "OD"
          else
            divertors%divertor(1)%identifier = 'USN'
            divertors%divertor(1)%target(1)%identifier = "OD"
            divertors%divertor(1)%target(2)%identifier = "ID"
          end if
#endif
          divertors%divertor(1)%target(1)%extension_r = extension_r(1)
          divertors%divertor(1)%target(1)%extension_z = extension_z(1)
          divertors%divertor(1)%target(2)%extension_r = extension_r(2)
          divertors%divertor(1)%target(2)%extension_z = extension_z(2)
          if (geo%LSN) then
            divertors%divertor(1)%target(1)%name = "Inner target"
            divertors%divertor(1)%target(2)%name = "Outer target"
          else
            divertors%divertor(1)%target(1)%name = "Outer target"
            divertors%divertor(1)%target(2)%name = "Inner target"
          end if
          do i = 1, maxval(mpg%strDiv)
            call write_timed_value( &
              &  divertors%divertor(1)%target(i)%flux_expansion, &
              &  flux_expansion(i) )
            call write_timed_value( &
              &  divertors%divertor(1)%target(i)%wetted_area, &
              &  wetted_area(i) )
            call write_timed_value( &
              &  divertors%divertor(1)%target(i)%power_conducted, &
              &  power_conducted(i) )
            call write_timed_value( &
              &  divertors%divertor(1)%target(i)%power_convected, &
              &  power_convected(i) )
            call write_timed_value( &
              &  divertors%divertor(1)%target(i)%power_recombination_plasma, &
              &  power_recombination_plasma(i) )
            call write_timed_value( &
              &  divertors%divertor(1)%target(i)%power_currents, &
              &  power_currents(i) )
#if ( IMAS_MINOR_VERSION > 32 || IMAS_MAJOR_VERSION > 3 )
            call write_timed_value( &
              &  divertors%divertor(1)%target(i)%current_incident, &
              &  current_incident(i) )
#endif
          end do
          call write_timed_value( &
            &  divertors%divertor(1)%wetted_area, &
            &  wetted_area(1)+wetted_area(2) )
          call write_timed_value( &
            &  divertors%divertor(1)%power_conducted, &
            &  power_conducted(1)+power_conducted(2) )
          call write_timed_value( &
            &  divertors%divertor(1)%power_convected, &
            &  power_convected(1)+power_convected(2) )
          call write_timed_value( &
            &  divertors%divertor(1)%power_recombination_plasma, &
            &  power_recombination_plasma(1)+power_recombination_plasma(2) )
          call write_timed_value( &
            &  divertors%divertor(1)%power_currents, &
            &  power_currents(1)+power_currents(2) )
          call write_timed_value( &
            &  divertors%divertor(1)%current_incident, &
            &  power_currents(1)+power_currents(2) )
#if ( IMAS_MINOR_VERSION > 32 || IMAS_MAJOR_VERSION > 3 )
          call write_timed_value( &
            &  divertors%divertor(1)%current_incident, &
            &  current_incident(1)+current_incident(2) )
#endif
        case ( GEOMETRY_CDN, GEOMETRY_DDN_BOTTOM, GEOMETRY_DDN_TOP, &
        &      GEOMETRY_LFS_SNOWFLAKE_MINUS, GEOMETRY_LFS_SNOWFLAKE_PLUS )
          allocate( divertors%divertor(2) )
          allocate( divertors%divertor(1)%name(1) )
          allocate( divertors%divertor(2)%name(1) )
          allocate( divertors%divertor(1)%target(2) )
          allocate( divertors%divertor(2)%target(2) )
          allocate( divertors%divertor(1)%target(1)%name(1) )
          allocate( divertors%divertor(1)%target(2)%name(1) )
          allocate( divertors%divertor(2)%target(1)%name(1) )
          allocate( divertors%divertor(2)%target(2)%name(1) )
#if IMAS_MAJOR_VERSION > 3
          allocate( divertors%divertor(1)%description(1) )
          allocate( divertors%divertor(2)%description(1) )
          allocate( divertors%divertor(1)%target(1)%description(1) )
          allocate( divertors%divertor(1)%target(2)%description(1) )
          allocate( divertors%divertor(2)%target(1)%description(1) )
          allocate( divertors%divertor(2)%target(2)%description(1) )
#else
          allocate( divertors%divertor(1)%identifier(1) )
          allocate( divertors%divertor(2)%identifier(1) )
          allocate( divertors%divertor(1)%target(1)%identifier(1) )
          allocate( divertors%divertor(1)%target(2)%identifier(1) )
          allocate( divertors%divertor(2)%target(1)%identifier(1) )
          allocate( divertors%divertor(2)%target(2)%identifier(1) )
#endif
          if (GeometryType == GEOMETRY_LFS_SNOWFLAKE_MINUS .or. &
          &   GeometryType == GEOMETRY_LFS_SNOWFLAKE_PLUS) then
            divertors%divertor(1)%name = 'Lower divertor'
            divertors%divertor(2)%name = 'Lower SF divertor'
            divertors%divertor(1)%target(1)%name = "Lower inner target"
            divertors%divertor(1)%target(2)%name = "Lower outer target"
            divertors%divertor(2)%target(1)%name = "Snowflake lower outer target"
            divertors%divertor(2)%target(2)%name = "Snowflake lower inner target"
#if IMAS_MAJOR_VERSION > 3
            divertors%divertor(1)%description = 'LD'
            divertors%divertor(2)%description = 'LSFD'
            divertors%divertor(1)%target(1)%description = "LID"
            divertors%divertor(1)%target(2)%description = "LOD"
            divertors%divertor(2)%target(1)%description = "LSFOD"
            divertors%divertor(2)%target(2)%description = "LSFID"
#else
            divertors%divertor(1)%identifier = 'LD'
            divertors%divertor(2)%identifier = 'LSFD'
            divertors%divertor(1)%target(1)%identifier = "LID"
            divertors%divertor(1)%target(2)%identifier = "LOD"
            divertors%divertor(2)%target(1)%identifier = "LSFOD"
            divertors%divertor(2)%target(2)%identifier = "LSFID"
#endif
          else
            divertors%divertor(1)%name = 'Lower divertor'
            divertors%divertor(2)%name = 'Upper divertor'
            divertors%divertor(1)%target(1)%name = "Lower inner target"
            divertors%divertor(1)%target(2)%name = "Lower outer target"
            divertors%divertor(2)%target(1)%name = "Upper inner target"
            divertors%divertor(2)%target(2)%name = "Upper outer target"
#if IMAS_MAJOR_VERSION > 3
            divertors%divertor(1)%description = 'LD'
            divertors%divertor(2)%description = 'UD'
            divertors%divertor(1)%target(1)%description = "LID"
            divertors%divertor(1)%target(2)%description = "LOD"
            divertors%divertor(2)%target(1)%description = "UID"
            divertors%divertor(2)%target(2)%description = "UOD"
#else
            divertors%divertor(1)%identifier = 'LD'
            divertors%divertor(2)%identifier = 'UD'
            divertors%divertor(1)%target(1)%identifier = "LID"
            divertors%divertor(1)%target(2)%identifier = "LOD"
            divertors%divertor(2)%target(1)%identifier = "UID"
            divertors%divertor(2)%target(2)%identifier = "UOD"
#endif
          endif
          do i = 1, maxval(mpg%strDiv)
            if (i.eq.1.or.i.eq.4) then
              k = 1
            else
              k = 2
            end if
            call write_timed_value( &
              &  divertors%divertor(k)%target(i)%flux_expansion, &
              &  flux_expansion(i) )
            call write_timed_value( &
              &  divertors%divertor(k)%target(i)%wetted_area, &
              &  wetted_area(i) )
            call write_timed_value( &
              &  divertors%divertor(k)%target(i)%power_conducted, &
              &  power_conducted(i) )
            call write_timed_value( &
              &  divertors%divertor(k)%target(i)%power_convected, &
              &  power_convected(i) )
            call write_timed_value( &
              &  divertors%divertor(k)%target(i)%power_recombination_plasma, &
              &  power_recombination_plasma(i) )
            call write_timed_value( &
              &  divertors%divertor(k)%target(i)%power_currents, &
              &  power_currents(i) )
#if ( IMAS_MINOR_VERSION > 32 || IMAS_MAJOR_VERSION > 3 )
            call write_timed_value( &
              &  divertors%divertor(k)%target(i)%current_incident, &
              &  current_incident(i) )
#endif
          end do
          call write_timed_value( &
            &  divertors%divertor(1)%wetted_area, &
            &  wetted_area(1)+wetted_area(4) )
          call write_timed_value( &
            &  divertors%divertor(1)%power_conducted, &
            &  power_conducted(1)+power_conducted(4) )
          call write_timed_value( &
            &  divertors%divertor(1)%power_convected, &
            &  power_convected(1)+power_convected(4) )
          call write_timed_value( &
            &  divertors%divertor(1)%power_recombination_plasma, &
            &  power_recombination_plasma(1)+power_recombination_plasma(4) )
          call write_timed_value( &
            &  divertors%divertor(1)%power_currents, &
            &  power_currents(1)+power_currents(4) )
          call write_timed_value( &
            &  divertors%divertor(1)%current_incident, &
            &  power_currents(1)+power_currents(4) )
#if ( IMAS_MINOR_VERSION > 32 || IMAS_MAJOR_VERSION > 3 )
          call write_timed_value( &
            &  divertors%divertor(1)%current_incident, &
            &  current_incident(1)+current_incident(4) )
#endif
          call write_timed_value( &
            &  divertors%divertor(2)%wetted_area, &
            &  wetted_area(2)+wetted_area(3) )
          call write_timed_value( &
            &  divertors%divertor(2)%power_conducted, &
            &  power_conducted(2)+power_conducted(3) )
          call write_timed_value( &
            &  divertors%divertor(2)%power_convected, &
            &  power_convected(2)+power_convected(3) )
          call write_timed_value( &
            &  divertors%divertor(2)%power_recombination_plasma, &
            &  power_recombination_plasma(2)+power_recombination_plasma(3) )
          call write_timed_value( &
            &  divertors%divertor(2)%power_currents, &
            &  power_currents(2)+power_currents(3) )
          call write_timed_value( &
            &  divertors%divertor(2)%current_incident, &
            &  power_currents(2)+power_currents(3) )
#if ( IMAS_MINOR_VERSION > 32 || IMAS_MAJOR_VERSION > 3 )
          call write_timed_value( &
            &  divertors%divertor(2)%current_incident, &
            &  current_incident(2)+current_incident(3) )
#endif
        end select
        deallocate( power_currents )
        deallocate( power_convected )
        deallocate( power_conducted )
        deallocate( power_recombination_plasma )
        deallocate( current_incident )
#endif
#endif

        !! Write grid & grid subsets/subgrids
#if ( IMAS_MINOR_VERSION > 11 || IMAS_MAJOR_VERSION > 3 ) && GGD_MAJOR_VERSION > 0
#if ( IMAS_MINOR_VERSION < 15 && IMAS_MAJOR_VERSION < 4 )
        call b2_IMAS_Fill_Grid_Desc( mpg, geo,                              &
            &   edge_profiles%ggd( time_sind )%grid )
        call b2_IMAS_Fill_Grid_Desc( mpg, geo,                              &
            &   edge_transport%model(1)%ggd( time_sind )%grid )
        do is = 1, nsources
            call b2_IMAS_Fill_Grid_Desc( mpg, geo,                          &
                &   edge_sources%source(is)%ggd( time_sind )%grid )
        end do
#else
        call b2_IMAS_Fill_Grid_Desc( mpg, geo,                              &
            &   edge_profiles%grid_ggd( time_sind ) )
#if AL_MAJOR_VERSION > 4
        allocate( edge_transport%grid_ggd( time_sind )%path(1) )
        edge_transport%grid_ggd( time_sind )%path = &
            &   "#edge_profiles/grid_ggd("//int2str(time_sind)//")"
        allocate( edge_sources%grid_ggd( time_sind )%path(1) )
        edge_sources%grid_ggd( time_sind )%path = &
            &   "#edge_profiles/grid_ggd("//int2str(time_sind)//")"
        allocate( radiation%grid_ggd( time_sind )%path(1) )
        radiation%grid_ggd( time_sind )%path = &
            &   "#edge_profiles/grid_ggd("//int2str(time_sind)//")"
#else
        call b2_IMAS_Fill_Grid_Desc( mpg, geo,                              &
            &   edge_transport%grid_ggd( time_sind ) )
        call b2_IMAS_Fill_Grid_Desc( mpg, geo,                              &
            &   edge_sources%grid_ggd( time_sind ) )
#if ( IMAS_MINOR_VERSION > 21 || IMAS_MAJOR_VERSION > 3 )
        call b2_IMAS_Fill_Grid_Desc( mpg, geo,                              &
            &   radiation%grid_ggd( time_sind ) )
#endif
#endif
#endif
#else
        write(0,*) 'Code was compiled without a GGD module'
        write(0,*) 'Most IDS output is disabled !'
#endif

        !! Allocate and set time slice value
#if ( IMAS_MINOR_VERSION > 14 || IMAS_MAJOR_VERSION > 3 )
        edge_profiles%grid_ggd( time_sind )%time = time_slice_value
        edge_transport%grid_ggd( time_sind )%time = time_slice_value
        edge_sources%grid_ggd( time_sind )%time = time_slice_value
#if ( IMAS_MINOR_VERSION > 21 || IMAS_MAJOR_VERSION > 3 )
        radiation%grid_ggd( time_sind )%time = time_slice_value
#endif
#endif
        edge_profiles%ggd( time_sind )%time = time_slice_value
        edge_transport%model(1)%ggd( time_sind )%time = time_slice_value
        do is = 1, nsources
          edge_sources%source(is)%ggd( time_sind )%time = time_slice_value
        end do
#if ( IMAS_MINOR_VERSION > 21 || IMAS_MAJOR_VERSION > 3 )
        do j = 1, n_process
          allocate( radiation%process(j)%ggd( num_time_slices ) )
          radiation%process(j)%ggd( time_sind )%time = time_slice_value
        end do
#endif

        !! List of species
        !! Careful here: ion in DD means isonuclear sequence !!
#ifdef B25_EIRENE
        allocate( ionstt(0:ns+nioni-1) )
        allocate( ispion(ns+nioni,max(ns,nioni)) )
#else
        allocate( ionstt(0:ns-1) )
        allocate( ispion(ns,ns) )
#endif
        ionstt = 0
        ispion = -1
        nsion = nspecies
        do is = 1, nspecies
          if (is.lt.nspecies) then
            ks = eb2spcr(is+1)-1
          else
            ks = ns-1
          end if
          js = 0
          do i = eb2spcr(is), ks
            if (is_neutral(i)) cycle
            js = js + 1
            ionstt(i) = js
            ispion(is,js) = i
          end do
        end do
#ifdef B25_EIRENE
        if (switch%use_eirene.ne.0) then
          allocate( imiion(nioni) )
          if (nioni.ge.1) then
            imiion = 0
            nsion = nsion + 1
            imiion(1) = nsion
            ispion(nspecies+1,1) = 1
            ionstt(ns) = 1
            do is = 2, nioni
              ks = 0
              do js = 1, is-1
                match_found = .false.
                do iatm = 1, natmi
                  match_found = match_found .and. &
                       &  micmp(iatm,js).eq.micmp(iatm,is)
                end do
                if (match_found) then
                  imiion(is)=imiion(js)
                  ks = js
                end if
              end do
              if (imiion(is).ne.0) then
                ionstt(ns-1+is) = ionstt(ns-1+ks)+1
                ispion(imiion(is),ionstt(ns-1+is)) = is
              else
                nsion = nsion + 1
                imiion(is) = nsion
                ionstt(ns-1+is) = 1
                ispion(imiion(is),1) = is
              end if
            end do
          end if
        end if
#endif
        allocate( istion( nsion ) )
        do is = 1, nspecies
          istion(is) = nfluids(is)
        end do
#ifdef B25_EIRENE
        if (switch%use_eirene.ne.0) then
          do is = nspecies+1, nsion
            do js = 1, nioni
              if (imiion(js).eq.is) istion(is) = ionstt(ns-1+js)
            end do
          end do
        end if
#endif
        allocate( edge_profiles%ggd( time_sind )%ion( nsion ) )
#if ( IMAS_MINOR_VERSION > 21 || IMAS_MAJOR_VERSION > 3 )
        allocate( radiation%process(1)%ggd( time_sind )%ion( nsion ) )
        allocate( radiation%process(2)%ggd( time_sind )%ion( nsion ) )
#endif
        do i = 1, nsources
          allocate( edge_sources%source(i)%ggd( time_sind )%ion( nsion ) )
        end do
        allocate( edge_transport%model(1)%ggd( time_sind )%ion( nsion ) )
        do js = 1, nspecies
          allocate( edge_profiles%ggd( time_sind )%ion( js )%state( nfluids(js) ) )
#if ( IMAS_MAJOR_VERSION < 4 && IMAS_MINOR_VERSION < 42 )
          allocate( edge_profiles%ggd( time_sind )%ion( js )%label(1) )
          do is = 1, nfluids(js)
            allocate( edge_profiles%ggd( time_sind )%ion( js )%state( is )%label(1) )
          end do
#else
          allocate( edge_profiles%ggd( time_sind )%ion( js )%name(1) )
          do is = 1, nfluids(js)
            allocate( edge_profiles%ggd( time_sind )%ion( js )%state( is )%name(1) )
          end do
#endif
          allocate( edge_profiles%ggd( time_sind )%ion( js )%element(1) )
          do i = 1, nsources
            allocate( edge_sources%source(i)%ggd( time_sind )%ion( js )%state( nfluids(js) ) )
#if ( IMAS_MAJOR_VERSION < 4 && IMAS_MINOR_VERSION < 42 )
            allocate( edge_sources%source(i)%ggd( time_sind )%ion( js )%label(1) )
            do is = 1, nfluids(js)
              allocate( edge_sources%source(i)%ggd( time_sind )%ion( js )%state( is )%label(1) )
            end do
#else
            allocate( edge_sources%source(i)%ggd( time_sind )%ion( js )%name(1) )
            do is = 1, nfluids(js)
              allocate( edge_sources%source(i)%ggd( time_sind )%ion( js )%state( is )%name(1) )
            end do
#endif
            allocate( edge_sources%source(i)%ggd( time_sind )%ion( js )%element(1) )
          end do
          allocate( edge_transport%model(1)%ggd( time_sind )%ion( js )%state( nfluids(js) ) )
#if ( IMAS_MAJOR_VERSION < 4 && IMAS_MINOR_VERSION < 42 )
          allocate( edge_transport%model(1)%ggd( time_sind )%ion( js )%label(1) )
          do is = 1, nfluids(js)
            allocate( edge_transport%model(1)%ggd( time_sind )%ion( js )%state( is )%label(1) )
          end do
#else
          allocate( edge_transport%model(1)%ggd( time_sind )%ion( js )%name(1) )
          do is = 1, nfluids(js)
            allocate( edge_transport%model(1)%ggd( time_sind )%ion( js )%state( is )%name(1) )
          end do
#endif
          allocate( edge_transport%model(1)%ggd( time_sind )%ion( js )%element(1) )

#if ( IMAS_MINOR_VERSION > 21 || IMAS_MAJOR_VERSION > 3 )
          allocate( radiation%process(1)%ggd( time_sind )%ion( js )%state( nfluids(js) ) )
          allocate( radiation%process(2)%ggd( time_sind )%ion( js )%state( nfluids(js) ) )
#if ( IMAS_MAJOR_VERSION < 4 && IMAS_MINOR_VERSION < 42 )
          allocate( radiation%process(1)%ggd( time_sind )%ion( js )%label(1) )
          allocate( radiation%process(2)%ggd( time_sind )%ion( js )%label(1) )
          do is = 1, nfluids(js)
            allocate( radiation%process(1)%ggd( time_sind )%ion( js )%state( is )%label(1) )
            allocate( radiation%process(2)%ggd( time_sind )%ion( js )%state( is )%label(1) )
          end do
#else
          allocate( radiation%process(1)%ggd( time_sind )%ion( js )%name(1) )
          allocate( radiation%process(2)%ggd( time_sind )%ion( js )%name(1) )
          do is = 1, nfluids(js)
            allocate( radiation%process(1)%ggd( time_sind )%ion( js )%state( is )%name(1) )
            allocate( radiation%process(2)%ggd( time_sind )%ion( js )%state( is )%name(1) )
          end do
#endif
          allocate( radiation%process(1)%ggd( time_sind )%ion( js )%element(1) )
          allocate( radiation%process(2)%ggd( time_sind )%ion( js )%element(1) )
#endif
          ! Put label to ion(js).state(is).label
          do is = 1, istion(js)
            ks = ispion(js,is)
            call species( ks, spclabel, .false.)
#if ( IMAS_MAJOR_VERSION < 4 && IMAS_MINOR_VERSION < 42 )
            edge_profiles%ggd( time_sind )%ion( js )%state( is )%label = spclabel
            do i = 1, nsources
              edge_sources%source(i)%ggd( time_sind )%ion( js )%state( is )%label = spclabel
            end do
            edge_transport%model(1)%ggd( time_sind )%ion( js )%state( is )%label = spclabel
#if ( IMAS_MINOR_VERSION > 21 || IMAS_MAJOR_VERSION > 3 )
            radiation%process(1)%ggd( time_sind )%ion( js )%state( is )%label = spclabel
            radiation%process(2)%ggd( time_sind )%ion( js )%state( is )%label = spclabel
#endif
#else
            edge_profiles%ggd( time_sind )%ion( js )%state( is )%name = spclabel
            do i = 1, nsources
              edge_sources%source(i)%ggd( time_sind )%ion( js )%state( is )%name = spclabel
            end do
            edge_transport%model(1)%ggd( time_sind )%ion( js )%state( is )%name = spclabel
#if ( IMAS_MINOR_VERSION > 21 || IMAS_MAJOR_VERSION > 3 )
            radiation%process(1)%ggd( time_sind )%ion( js )%state( is )%name = spclabel
            radiation%process(2)%ggd( time_sind )%ion( js )%state( is )%name = spclabel
#endif
#endif
          end do

          ! Put ion label identifying the species
#if ( IMAS_MAJOR_VERSION < 4 && IMAS_MINOR_VERSION < 42 )
          edge_profiles%ggd( time_sind )%ion( js )%label = species_list( js )
          edge_transport%model(1)%ggd( time_sind )%ion( js )%label = species_list( js )
#else
          edge_profiles%ggd( time_sind )%ion( js )%name = species_list( js )
          edge_transport%model(1)%ggd( time_sind )%ion( js )%name = species_list( js )
#endif

          ! Put ion charge if single ion in species
          is = ispion(js,1)
          if (istion(js).eq.1) then
            ion_charge_int = nint((zamin(is)+zamax(is))/2.0_R8)
            edge_profiles%ggd( time_sind )%ion( js )%z_ion = ion_charge_int
            edge_transport%model(1)%ggd( time_sind )%ion( js )%z_ion = ion_charge_int
          end if

          ! Put mass of species
          edge_profiles%ggd( time_sind )%ion( js )%element(1)%a = am( is )
          edge_transport%model(1)%ggd( time_sind )%ion( js )%element(1)%a = am( is )

          ! Put nuclear charge
#if IMAS_MAJOR_VERSION < 4
          edge_profiles%ggd( time_sind )%ion( js )%element(1)%z_n = zn( is )
          edge_transport%model(1)%ggd( time_sind )%ion( js )%element(1)%z_n = zn( is )
#else
          edge_profiles%ggd( time_sind )%ion( js )%element(1)%z_n = nint(zn(is))
          edge_transport%model(1)%ggd( time_sind )%ion( js )%element(1)%z_n = nint(zn(is))
#endif

          ! Put number of atoms
#if ( IMAS_MINOR_VERSION < 15 && IMAS_MAJOR_VERSION < 4 )
          edge_profiles%ggd( time_sind )%ion( js )%element(1)%multiplicity = 1.0_IDS_real
          edge_transport%model(1)%ggd( time_sind )%ion( js )%element(1)%multiplicity = 1.0_IDS_real
#else
          edge_profiles%ggd( time_sind )%ion( js )%element(1)%atoms_n = 1
          edge_transport%model(1)%ggd( time_sind )%ion( js )%element(1)%atoms_n = 1
#endif

          ! Put neutral index
          edge_profiles%ggd( time_sind )%ion( js )%neutral_index = b2eatcr(is)
          edge_transport%model(1)%ggd( time_sind )%ion( js )%neutral_index = b2eatcr(is)

          ! Put multiple states flag
          edge_profiles%ggd( time_sind )%ion( js )%multiple_states_flag = 1
          edge_transport%model(1)%ggd( time_sind )%ion( js )%multiple_states_flag = 1

          do is = 1, istion(js)
            ks = ispion(js,is)

            ! Put minimum Z of the charge state bundle
            ! (z_min = z_max = 0 for a neutral)
            edge_profiles%ggd( time_sind )%ion( js )%state( is )%z_min = zamin( ks )
            edge_transport%model(1)%ggd( time_sind )%ion( js )%state( is )%z_min = zamin( ks )

            ! Put maximum Z of the charge state bundle
            edge_profiles%ggd( time_sind )%ion( js )%state( is )%z_max = zamax( ks )
            edge_transport%model(1)%ggd( time_sind )%ion( js )%state( is )%z_max = zamax( ks )
          end do

          do i = 1, nsources
            ! Put ion label identifying the species
#if ( IMAS_MAJOR_VERSION < 4 && IMAS_MINOR_VERSION < 42 )
            edge_sources%source(i)%ggd( time_sind )%ion( js )%label = species_list( js )
#else
            edge_sources%source(i)%ggd( time_sind )%ion( js )%name = species_list( js )
#endif

            ! Put ion charge if single ion in species
            is = ispion(js,1)
            if (istion(js).eq.1) then
              ion_charge_int = nint((zamin(is)+zamax(is))/2.0_R8)
              edge_sources%source(i)%ggd( time_sind )%ion( js )%z_ion = ion_charge_int
            end if

            ! Put mass of ion
            edge_sources%source(i)%ggd( time_sind )%ion( js )%element(1)%a = am( is )

            ! Put nuclear charge
#if IMAS_MAJOR_VERSION < 4
            edge_sources%source(i)%ggd( time_sind )%ion( js )%element(1)%z_n = zn( is )
#else
            edge_sources%source(i)%ggd( time_sind )%ion( js )%element(1)%z_n = nint(zn(is))
#endif

            ! Put number of atoms
#if ( IMAS_MINOR_VERSION < 15 && IMAS_MAJOR_VERSION < 4 )
            edge_sources%source(i)%ggd( time_sind )%ion( js )%element(1)%multiplicity = 1.0_IDS_real
#else
            edge_sources%source(i)%ggd( time_sind )%ion( js )%element(1)%atoms_n = 1
#endif

            ! Put neutral index
            edge_sources%source(i)%ggd( time_sind )%ion( js )%neutral_index = b2eatcr(is)

            ! Put multiple states flag
            edge_sources%source(i)%ggd( time_sind )%ion( js )%multiple_states_flag = 1

            do is = 1, istion(js)
              ks = ispion(js,is)

              ! Put minimum Z of the charge state bundle
              ! (z_min = z_max = 0 for a neutral)
              edge_sources%source(i)%ggd( time_sind )%ion( js )%state( is )%z_min = &
                  &  zamin( ks )

              ! Put maximum Z of the charge state bundle
              edge_sources%source(i)%ggd( time_sind )%ion( js )%state( is )%z_max = &
                  &  zamax( ks )
            end do
          end do

#if ( IMAS_MINOR_VERSION > 21 || IMAS_MAJOR_VERSION > 3 )
          ! Put ion label identifying the species
#if ( IMAS_MAJOR_VERSION < 4 && IMAS_MINOR_VERSION < 42 )
          radiation%process(1)%ggd( time_sind )%ion( js )%label = species_list( js )
          radiation%process(2)%ggd( time_sind )%ion( js )%label = species_list( js )
#else
          radiation%process(1)%ggd( time_sind )%ion( js )%name = species_list( js )
          radiation%process(2)%ggd( time_sind )%ion( js )%name = species_list( js )
#endif

          ! Put ion charge if single ion in species
          is = ispion(js,1)
          if (istion(js).eq.1) then
            ion_charge_int = nint((zamin(is)+zamax(is))/2.0_R8)
            radiation%process(1)%ggd( time_sind )%ion( js )%z_ion = ion_charge_int
            radiation%process(2)%ggd( time_sind )%ion( js )%z_ion = ion_charge_int
          end if

          ! Put mass of ion
          radiation%process(1)%ggd( time_sind )%ion( js )%element(1)%a = am( is )
          radiation%process(2)%ggd( time_sind )%ion( js )%element(1)%a = am( is )

          ! Put nuclear charge
#if IMAS_MAJOR_VERSION < 4
          radiation%process(1)%ggd( time_sind )%ion( js )%element(1)%z_n = zn( is )
          radiation%process(2)%ggd( time_sind )%ion( js )%element(1)%z_n = zn( is )
#else
          radiation%process(1)%ggd( time_sind )%ion( js )%element(1)%z_n = nint(zn(is))
          radiation%process(2)%ggd( time_sind )%ion( js )%element(1)%z_n = nint(zn(is))
#endif

          ! Put number of atoms
          radiation%process(1)%ggd( time_sind )%ion( js )%element(1)%atoms_n = 1
          radiation%process(2)%ggd( time_sind )%ion( js )%element(1)%atoms_n = 1

          ! Put neutral index
          radiation%process(1)%ggd( time_sind )%ion( js )%neutral_index = b2eatcr(is)
          radiation%process(2)%ggd( time_sind )%ion( js )%neutral_index = b2eatcr(is)

          ! Put multiple states flag
          radiation%process(1)%ggd( time_sind )%ion( js )%multiple_states_flag = 1
          radiation%process(2)%ggd( time_sind )%ion( js )%multiple_states_flag = 1

          do is = 1, istion(js)
            ks = ispion(js,is)

            ! Put minimum Z of the charge state bundle
            ! (z_min = z_max = 0 for a neutral)
            radiation%process(1)%ggd( time_sind )%ion( js )%state( is )%z_min = zamin( ks )
            radiation%process(2)%ggd( time_sind )%ion( js )%state( is )%z_min = zamin( ks )

            ! Put maximum Z of the charge state bundle
            radiation%process(1)%ggd( time_sind )%ion( js )%state( is )%z_max = zamax( ks )
            radiation%process(2)%ggd( time_sind )%ion( js )%state( is )%z_max = zamax( ks )
          end do
#endif
        enddo

        if (switch%use_eirene.ne.0) then
#ifdef B25_EIRENE
          do js = nspecies+1, nsion
            allocate( edge_profiles%ggd( time_sind )%ion( js )%state( istion(js) ) )
#if ( IMAS_MAJOR_VERSION < 4 && IMAS_MINOR_VERSION < 42 )
            allocate( edge_profiles%ggd( time_sind )%ion( js )%label(1) )
            do ks = 1, istion(js)
              allocate( edge_profiles%ggd( time_sind )%ion( js )%state( ks )%label(1) )
            end do
#else
            allocate( edge_profiles%ggd( time_sind )%ion( js )%name(1) )
            do ks = 1, istion(js)
              allocate( edge_profiles%ggd( time_sind )%ion( js )%state( ks )%name(1) )
            end do
#endif
            do i = 1, nsources
              allocate( edge_sources%source(i)%ggd( time_sind )%ion( js )%state( istion(js) ) )
#if ( IMAS_MAJOR_VERSION < 4 && IMAS_MINOR_VERSION < 42 )
              allocate( edge_sources%source(i)%ggd( time_sind )%ion( js )%label(1) )
              do ks = 1, istion(js)
                allocate( edge_sources%source(i)%ggd( time_sind )%ion( js )%state( ks )%label(1) )
              end do
#else
              allocate( edge_sources%source(i)%ggd( time_sind )%ion( js )%name(1) )
              do ks = 1, istion(js)
                allocate( edge_sources%source(i)%ggd( time_sind )%ion( js )%state( ks )%name(1) )
              end do
#endif
            end do
            allocate( edge_transport%model(1)%ggd( time_sind )%ion( js )%state( istion(js) ) )
#if ( IMAS_MAJOR_VERSION < 4 && IMAS_MINOR_VERSION < 42 )
            allocate( edge_transport%model(1)%ggd( time_sind )%ion( js )%label(1) )
            do ks = 1, istion(js)
              allocate( edge_transport%model(1)%ggd( time_sind )%ion( js )%state( ks )%label(1) )
            end do
#else
            allocate( edge_transport%model(1)%ggd( time_sind )%ion( js )%name(1) )
            do ks = 1, istion(js)
              allocate( edge_transport%model(1)%ggd( time_sind )%ion( js )%state( ks )%name(1) )
            end do
#endif
            do ks = 1, istion(js)
              is = ispion(js,ks)
#if ( IMAS_MAJOR_VERSION < 4 && IMAS_MINOR_VERSION < 42 )
              edge_profiles%ggd( time_sind )%ion( js )%state( ks )%label = textin( is-1 )
              edge_transport%model(1)%ggd( time_sind )%ion( js )%state( ks )%label = &
                  &                                                        textin( is-1 )
              do i = 1, nsources
                edge_sources%source(i)%ggd( time_sind )%ion( js )%state( ks )%label = &
                    &                                                      textin( is-1 )
              end do
#else
              edge_profiles%ggd( time_sind )%ion( js )%state( ks )%name = textin( is-1 )
              edge_transport%model(1)%ggd( time_sind )%ion( js )%state( ks )%name = &
                  &                                                        textin( is-1 )
              do i = 1, nsources
                edge_sources%source(i)%ggd( time_sind )%ion( js )%state( ks )%name = &
                    &                                                      textin( is-1 )
              end do
#endif
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
#if IMAS_MAJOR_VERSION < 4
                  edge_profiles%ggd( time_sind )%ion( js )%element( i )%z_n = real(nchara(k),IDS_real)
#else
                  edge_profiles%ggd( time_sind )%ion( js )%element( i )%z_n = nchara(k)
#endif
#if ( IMAS_MINOR_VERSION < 15 && IMAS_MAJOR_VERSION < 4 )
                  edge_profiles%ggd( time_sind )%ion( js )%element( i )%multiplicity = micmp(k,is)
#else
                  edge_profiles%ggd( time_sind )%ion( js )%element( i )%atoms_n = micmp(k,is)
#endif
                  do ii = 1, nsources
                    edge_sources%source(ii)%ggd( time_sind )%ion( js )%element( i )%a = nmassa(k)
#if IMAS_MAJOR_VERSION < 4
                    edge_sources%source(ii)%ggd( time_sind )%ion( js )%element( i )%z_n = real(nchara(k),IDS_real)
#else
                    edge_sources%source(ii)%ggd( time_sind )%ion( js )%element( i )%z_n = nchara(k)
#endif
#if ( IMAS_MINOR_VERSION < 15 && IMAS_MAJOR_VERSION < 4 )
                    edge_sources%source(ii)%ggd( time_sind )%ion( js )%element( i )%multiplicity = &
                        &   micmp(k,is)
#else
                    edge_sources%source(ii)%ggd( time_sind )%ion( js )%element( i )%atoms_n = &
                        &   micmp(k,is)
#endif
                  end do
                  edge_transport%model(1)%ggd( time_sind )%ion( js )%element( i )%a = nmassa(k)
#if IMAS_MAJOR_VERSION < 4
                  edge_transport%model(1)%ggd( time_sind )%ion( js )%element( i )%z_n = real(nchara(k),IDS_real)
#else
                  edge_transport%model(1)%ggd( time_sind )%ion( js )%element( i )%z_n = nchara(k)
#endif
#if ( IMAS_MINOR_VERSION < 15 && IMAS_MAJOR_VERSION < 4 )
                  edge_transport%model(1)%ggd( time_sind )%ion( js )%element( i )%multiplicity = micmp(k,is)
#else
                  edge_transport%model(1)%ggd( time_sind )%ion( js )%element( i )%atoms_n = micmp(k,is)
#endif
                end if
              end do
            end do
            is = ispion(js,1)
            if (istion(js).eq.1) then
              edge_profiles%ggd( time_sind )%ion( js )%z_ion = nchrgi( is )
#if ( IMAS_MAJOR_VERSION < 4 && IMAS_MINOR_VERSION < 42 )
              edge_profiles%ggd( time_sind )%ion( js )%label = textin( is-1 )
#else
              edge_profiles%ggd( time_sind )%ion( js )%name = textin( is-1 )
#endif
              do i = 1, nsources
                edge_sources%source(i)%ggd( time_sind )%ion( js )%z_ion = nchrgi( is )
#if ( IMAS_MAJOR_VERSION < 4 && IMAS_MINOR_VERSION < 42 )
                edge_sources%source(i)%ggd( time_sind )%ion( js )%label = textin( is-1 )
#else
                edge_sources%source(i)%ggd( time_sind )%ion( js )%name = textin( is-1 )
#endif
              end do
              edge_transport%model(1)%ggd( time_sind )%ion( js )%z_ion = nchrgi( is )
#if ( IMAS_MAJOR_VERSION < 4 && IMAS_MINOR_VERSION < 42 )
              edge_transport%model(1)%ggd( time_sind )%ion( js )%label = textin( is-1 )
#else
              edge_transport%model(1)%ggd( time_sind )%ion( js )%name = textin( is-1 )
#endif
            else
              match_found = .false.
              do ks = 2, istion(js)
                match_found = match_found .and. nchrgi(is).eq.nchrgi(ispion(js,ks))
              end do
              if (match_found) then
                edge_profiles%ggd( time_sind )%ion( js )%z_ion = nchrgi( is )
                do i = 1, nsources
                  edge_sources%source(i)%ggd( time_sind )%ion( js )%z_ion = nchrgi( is )
                end do
                edge_transport%model(1)%ggd( time_sind )%ion( js )%z_ion = nchrgi( is )
              end if
              match_found = .false.
              p = index(textin(ispion(js,1)-1),'+')
              if (p.gt.1) then
                ion_label = textin(ispion(js,1)-1)(1:p-1)
              else
                ion_label = textin(ispion(js,1)-1)
              end if
              do ks = 2, istion(js)
                p = index(textin(ispion(js,ks)-1),'+')
                if (p.gt.1) then
                  match_found = match_found .and. &
                    &  streql(ion_label,textin(ispion(js,ks)-1)(1:p-1))
                else
                  match_found = match_found .and. &
                    &  streql(ion_label,textin(ispion(js,ks)-1))
                end if
                if (match_found) then
#if ( IMAS_MAJOR_VERSION < 4 && IMAS_MINOR_VERSION < 42 )
                  edge_profiles%ggd( time_sind )%ion( js )%label = ion_label
                  do i = 1, nsources
                    edge_sources%source(i)%ggd( time_sind )%ion( js )%label = ion_label
                  end do
                  edge_transport%model(1)%ggd( time_sind )%ion( js )%label = ion_label
#else
                  edge_profiles%ggd( time_sind )%ion( js )%name = ion_label
                  do i = 1, nsources
                    edge_sources%source(i)%ggd( time_sind )%ion( js )%name = ion_label
                  end do
                  edge_transport%model(1)%ggd( time_sind )%ion( js )%name = ion_label
#endif
                end if
              end do
            end if
            is = ispion(js,1)
            ion_label = adjustl(textin(is-1))
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
            if (match_found) then
              jj = nspecies + j
            else
              jj = lkindi( is )
            end if
            edge_profiles%ggd( time_sind )%ion( js )%neutral_index = jj
            edge_profiles%ggd( time_sind )%ion( js )%multiple_states_flag = 1
            do ks = 1, istion(js)
              is = ispion(js,ks)
              edge_profiles%ggd( time_sind )%ion( js )%state( ks )%z_min = nchrgi( is )
              edge_profiles%ggd( time_sind )%ion( js )%state( ks )%z_max = nchrgi( is )
            end do
            is = ispion(js,1)
            do i = 1, nsources
              edge_sources%source(i)%ggd( time_sind )%ion( js )%neutral_index = jj
              edge_sources%source(i)%ggd( time_sind )%ion( js )%multiple_states_flag = 1
              do ks = 1, istion(js)
                is = ispion(js,ks)
                edge_sources%source(i)%ggd( time_sind )%ion( js )%state( ks )%z_min = nchrgi( is )
                edge_sources%source(i)%ggd( time_sind )%ion( js )%state( ks )%z_max = nchrgi( is )
              end do
            end do
            is = ispion(js,1)
            edge_transport%model(1)%ggd( time_sind )%ion( js )%neutral_index = jj
            edge_transport%model(1)%ggd( time_sind )%ion( js )%multiple_states_flag = 1
            do ks = 1, istion(js)
              is = ispion(js,ks)
              edge_transport%model(1)%ggd( time_sind )%ion( js )%state( ks )%z_min = nchrgi( is )
              edge_transport%model(1)%ggd( time_sind )%ion( js )%state( ks )%z_max = nchrgi( is )
            end do
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
#if ( IMAS_MAJOR_VERSION < 4 && IMAS_MINOR_VERSION < 42 )
             allocate( edge_profiles%ggd( time_sind )%neutral( js )%label(1) )
             edge_profiles%ggd( time_sind )%neutral( js )%label = species_list( js )
#else
             allocate( edge_profiles%ggd( time_sind )%neutral( js )%name(1) )
             edge_profiles%ggd( time_sind )%neutral( js )%name = species_list( js )
#endif
             edge_profiles%ggd( time_sind )%neutral( js )%element(1)%a = am( is )
#if IMAS_MAJOR_VERSION < 4
             edge_profiles%ggd( time_sind )%neutral( js )%element(1)%z_n = zn( is )
#else
             edge_profiles%ggd( time_sind )%neutral( js )%element(1)%z_n = nint(zn(is))
#endif
#if ( IMAS_MINOR_VERSION < 15 && IMAS_MAJOR_VERSION < 4 )
             edge_profiles%ggd( time_sind )%neutral( js )%element(1)%multiplicity = 1.0_IDS_real
#else
             edge_profiles%ggd( time_sind )%neutral( js )%element(1)%atoms_n = 1
#endif
             edge_profiles%ggd( time_sind )%neutral( js )%ion_index = js
             do i = 1, nsources
                allocate( edge_sources%source(i)%ggd( time_sind )%neutral( js )%element(1) )
#if ( IMAS_MAJOR_VERSION < 4 && IMAS_MINOR_VERSION < 42 )
                allocate( edge_sources%source(i)%ggd( time_sind )%neutral( js )%label(1) )
                edge_sources%source(i)%ggd( time_sind )%neutral( js )%label = species_list( js )
#else
                allocate( edge_sources%source(i)%ggd( time_sind )%neutral( js )%name(1) )
                edge_sources%source(i)%ggd( time_sind )%neutral( js )%name = species_list( js )
#endif
                edge_sources%source(i)%ggd( time_sind )%neutral( js )%element(1)%a = am( is )
#if IMAS_MAJOR_VERSION < 4
                edge_sources%source(i)%ggd( time_sind )%neutral( js )%element(1)%z_n = zn( is )
#else
                edge_sources%source(i)%ggd( time_sind )%neutral( js )%element(1)%z_n = nint(zn(is))
#endif
#if ( IMAS_MINOR_VERSION < 15 && IMAS_MAJOR_VERSION < 4 )
                edge_sources%source(i)%ggd( time_sind )%neutral( js )%element(1)%multiplicity = 1.0_IDS_real
#else
                edge_sources%source(i)%ggd( time_sind )%neutral( js )%element(1)%atoms_n = 1
#endif
                edge_sources%source(i)%ggd( time_sind )%neutral( js )%ion_index = js
             end do
             allocate( edge_transport%model(1)%ggd( time_sind )%neutral( js )%element(1) )
#if ( IMAS_MAJOR_VERSION < 4 && IMAS_MINOR_VERSION < 42 )
             allocate( edge_transport%model(1)%ggd( time_sind )%neutral( js )%label(1) )
             edge_transport%model(1)%ggd( time_sind )%neutral( js )%label = species_list( js )
#else
             allocate( edge_transport%model(1)%ggd( time_sind )%neutral( js )%name(1) )
             edge_transport%model(1)%ggd( time_sind )%neutral( js )%name = species_list( js )
#endif
             edge_transport%model(1)%ggd( time_sind )%neutral( js )%element(1)%a = am( is )
#if IMAS_MAJOR_VERSION < 4
             edge_transport%model(1)%ggd( time_sind )%neutral( js )%element(1)%z_n = zn( is )
#else
             edge_transport%model(1)%ggd( time_sind )%neutral( js )%element(1)%z_n = nint(zn(is))
#endif
#if ( IMAS_MINOR_VERSION < 15 && IMAS_MAJOR_VERSION < 4 )
             edge_transport%model(1)%ggd( time_sind )%neutral( js )%element(1)%multiplicity = 1.0_IDS_real
#else
             edge_transport%model(1)%ggd( time_sind )%neutral( js )%element(1)%atoms_n = 1
#endif
             edge_transport%model(1)%ggd( time_sind )%neutral( js )%ion_index = js
             ks = size(edge_profiles%ggd( time_sind )%neutral( js )%state)
             do i = 1, nsources
                allocate( edge_sources%source(i)%ggd( time_sind )%neutral( js )%state( ks ) )
             end do
             allocate( edge_transport%model(1)%ggd( time_sind )%neutral( js )%state( ks ) )
             do iss = 1, ks
                iatm = b2eatcr(is) + iss - 1
#if ( IMAS_MAJOR_VERSION < 4 && IMAS_MINOR_VERSION < 42 )
                allocate( edge_profiles%ggd( time_sind )%neutral( js )%state( iss )%label(1) )
                edge_profiles%ggd( time_sind )%neutral( js )%state( iss )%label = textan( iatm-1 )
#else
                allocate( edge_profiles%ggd( time_sind )%neutral( js )%state( iss )%name(1) )
                edge_profiles%ggd( time_sind )%neutral( js )%state( iss )%name = textan( iatm-1 )
#endif
                allocate( edge_profiles%ggd( time_sind )%neutral( js )%state( iss )% &
                   &      neutral_type%name(1) )
                allocate( edge_profiles%ggd( time_sind )%neutral( js )%state( iss )% &
                   &      neutral_type%description(1) )
                edge_profiles%ggd( time_sind )%neutral( js )%state( iss )%neutral_type%name = &
                   &     "Kinetic"
                edge_profiles%ggd( time_sind )%neutral( js )%state( iss )%neutral_type%index = -1
                edge_profiles%ggd( time_sind )%neutral( js )%state( iss )%neutral_type%description = &
                   &     "Kinetic neutral atoms from Eirene"
                edge_profiles%ggd( time_sind )%neutral( js )%multiple_states_flag = 1
                do i = 1, nsources
#if ( IMAS_MAJOR_VERSION < 4 && IMAS_MINOR_VERSION < 42 )
                   allocate( edge_sources%source(i)%ggd( time_sind )%neutral( js )%state( iss )%label(1) )
                   edge_sources%source(i)%ggd( time_sind )%neutral( js )%state( iss )%label = &
                     &   textan( iatm-1 )
#else
                   allocate( edge_sources%source(i)%ggd( time_sind )%neutral( js )%state( iss )%name(1) )
                   edge_sources%source(i)%ggd( time_sind )%neutral( js )%state( iss )%name = &
                     &   textan( iatm-1 )
#endif
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
                end do
#if ( IMAS_MAJOR_VERSION < 4 && IMAS_MINOR_VERSION < 42 )
                allocate( edge_transport%model(1)%ggd( time_sind )%neutral( js )%state( iss )%label(1) )
                edge_transport%model(1)%ggd( time_sind )%neutral( js )%state( iss )%label = textan( iatm-1 )
#else
                allocate( edge_transport%model(1)%ggd( time_sind )%neutral( js )%state( iss )%name(1) )
                edge_transport%model(1)%ggd( time_sind )%neutral( js )%state( iss )%name = textan( iatm-1 )
#endif
                allocate( edge_transport%model(1)%ggd( time_sind )%neutral( js )%state( iss )% &
                   &      neutral_type%name(1) )
                allocate( edge_transport%model(1)%ggd( time_sind )%neutral( js )%state( iss )% &
                   &      neutral_type%description(1) )
                edge_transport%model(1)%ggd( time_sind )%neutral( js )%state( iss )%neutral_type%name = &
                   &     "Kinetic"
                edge_transport%model(1)%ggd( time_sind )%neutral( js )%state( iss )%neutral_type%index = -1
                edge_transport%model(1)%ggd( time_sind )%neutral( js )%state( iss )%neutral_type%description = &
                   &     "Kinetic neutral atoms from Eirene"
                edge_transport%model(1)%ggd( time_sind )%neutral( js )%multiple_states_flag = 1
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
#if ( IMAS_MAJOR_VERSION < 4 && IMAS_MINOR_VERSION < 42 )
               allocate( edge_profiles%ggd( time_sind )%neutral( js )%label(1) )
               allocate( edge_profiles%ggd( time_sind )%neutral( js )%state( ks )%label(1) )
#else
               allocate( edge_profiles%ggd( time_sind )%neutral( js )%name(1) )
               allocate( edge_profiles%ggd( time_sind )%neutral( js )%state( ks )%name(1) )
#endif
               do i = 1, nsources
                  allocate( edge_sources%source(i)%ggd( time_sind )%neutral( js )%element( nelems ) )
#if ( IMAS_MAJOR_VERSION < 4 && IMAS_MINOR_VERSION < 42 )
                  allocate( edge_sources%source(i)%ggd( time_sind )%neutral( js )%label(1) )
                  allocate( edge_sources%source(i)%ggd( time_sind )%neutral( js )%state( ks )%label(1) )
#else
                  allocate( edge_sources%source(i)%ggd( time_sind )%neutral( js )%name(1) )
                  allocate( edge_sources%source(i)%ggd( time_sind )%neutral( js )%state( ks )%name(1) )
#endif
               end do
               allocate( edge_transport%model(1)%ggd( time_sind )%neutral( js )%element( nelems ) )
#if ( IMAS_MAJOR_VERSION < 4 && IMAS_MINOR_VERSION < 42 )
               allocate( edge_transport%model(1)%ggd( time_sind )%neutral( js )%label(1) )
               allocate( edge_transport%model(1)%ggd( time_sind )%neutral( js )%state( ks )%label(1) )
#else
               allocate( edge_transport%model(1)%ggd( time_sind )%neutral( js )%name(1) )
               allocate( edge_transport%model(1)%ggd( time_sind )%neutral( js )%state( ks )%name(1) )
#endif
               i = 0
               do k = 1, natmi
                  if (mlcmp( k, j ) > 0 ) then
                      i = i + 1
                      edge_profiles%ggd( time_sind )%neutral( js )%element(i)%a = nmassa( k )
#if IMAS_MAJOR_VERSION < 4
                      edge_profiles%ggd( time_sind )%neutral( js )%element(i)%z_n = &
                         &   real(nchara(k),IDS_real)
#else
                      edge_profiles%ggd( time_sind )%neutral( js )%element(i)%z_n = &
                         &   nchara( k )
#endif
#if ( IMAS_MINOR_VERSION < 15 && IMAS_MAJOR_VERSION < 4 )
                      edge_profiles%ggd( time_sind )%neutral( js )%element(i)%multiplicity = &
                         &   mlcmp(k,j)
#else
                      edge_profiles%ggd( time_sind )%neutral( js )%element(i)%atoms_n = &
                         &   mlcmp(k,j)
#endif
                      do icnt = 1, nsources
                         edge_sources%source(icnt)%ggd( time_sind )%neutral( js )%element(i)%a = &
                            &   nmassa( k )
#if IMAS_MAJOR_VERSION < 4
                         edge_sources%source(icnt)%ggd( time_sind )%neutral( js )%element(i)%z_n = &
                            &   real(nchara(k),IDS_real)
#else
                         edge_sources%source(icnt)%ggd( time_sind )%neutral( js )%element(i)%z_n = &
                            &   nchara( k )
#endif
                         edge_sources%source(icnt)%ggd( time_sind )%neutral( js )%element(i)%atoms_n = &
                            &   mlcmp(k,j)
                      end do
                      edge_transport%model(1)%ggd( time_sind )%neutral( js )%element(i)%a = nmassa( k )
#if IMAS_MAJOR_VERSION < 4
                      edge_transport%model(1)%ggd( time_sind )%neutral( js )%element(i)%z_n = &
                         &   real(nchara(k),IDS_real)
#else
                      edge_transport%model(1)%ggd( time_sind )%neutral( js )%element(i)%z_n = &
                         &   nchara( k )
#endif
#if ( IMAS_MINOR_VERSION < 15 && IMAS_MAJOR_VERSION < 4 )
                      edge_transport%model(1)%ggd( time_sind )%neutral( js )%element(i)%multiplicity = &
                         &   mlcmp(k,j)
#else
                      edge_transport%model(1)%ggd( time_sind )%neutral( js )%element(i)%atoms_n = &
                         &   mlcmp(k,j)
#endif
                   end if
               end do
#if ( IMAS_MAJOR_VERSION < 4 && IMAS_MINOR_VERSION < 42 )
               edge_profiles%ggd( time_sind )%neutral( js )%label = textmn( j-1 )
#else
               edge_profiles%ggd( time_sind )%neutral( js )%name = textmn( j-1 )
#endif
               ion_label = trim(textmn( j-1 ))//'+'
               i = 1
               match_found = .false.
               do while (.not.match_found.and.i.le.nioni)
                 if (streql(ion_label,textin(i-1))) then
                   match_found = .true.
                 end if
                 if (.not.match_found) i = i+1
               end do
               if (match_found) then
                 k = nspecies+i
               else
                 k = latmscl(lmolscl(j))
               end if
               edge_profiles%ggd( time_sind )%neutral( js )%ion_index = k
               edge_profiles%ggd( time_sind )%neutral( js )%multiple_states_flag = 1
#if ( IMAS_MAJOR_VERSION < 4 && IMAS_MINOR_VERSION < 42 )
               edge_profiles%ggd( time_sind )%neutral( js )%state( ks )%label = &
                   &    textmn( j-1 )
#else
               edge_profiles%ggd( time_sind )%neutral( js )%state( ks )%name = &
                   &    textmn( j-1 )
#endif
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
#if ( IMAS_MAJOR_VERSION < 4 && IMAS_MINOR_VERSION < 42 )
                  edge_sources%source(i)%ggd( time_sind )%neutral( js )%label = textmn( j-1 )
                  edge_sources%source(i)%ggd( time_sind )%neutral( js )%state( ks )%label = &
                      &    textmn( j-1 )
#else
                  edge_sources%source(i)%ggd( time_sind )%neutral( js )%name = textmn( j-1 )
                  edge_sources%source(i)%ggd( time_sind )%neutral( js )%state( ks )%name = &
                      &    textmn( j-1 )
#endif
                  edge_sources%source(i)%ggd( time_sind )%neutral( js )%ion_index = k
                  edge_sources%source(i)%ggd( time_sind )%neutral( js )%multiple_states_flag = 1
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
#if ( IMAS_MAJOR_VERSION < 4 && IMAS_MINOR_VERSION < 42 )
               edge_transport%model(1)%ggd( time_sind )%neutral( js )%label = textmn( j-1 )
               edge_transport%model(1)%ggd( time_sind )%neutral( js )%state( ks )%label = &
                   &    textmn( j-1 )
#else
               edge_transport%model(1)%ggd( time_sind )%neutral( js )%name = textmn( j-1 )
               edge_transport%model(1)%ggd( time_sind )%neutral( js )%state( ks )%name = &
                   &    textmn( j-1 )
#endif
               edge_transport%model(1)%ggd( time_sind )%neutral( js )%ion_index = k
               edge_transport%model(1)%ggd( time_sind )%neutral( js )%multiple_states_flag = 1
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
            nneut = 0
            do is = 1, nspecies
              if (is_neutral(eb2spcr(is))) nneut=nneut+1
            end do
            allocate( edge_profiles%ggd( time_sind )%neutral( nneut ) )
            do i = 1, nsources
               allocate( edge_sources%source(i)%ggd( time_sind )%neutral( nneut ) )
            end do
            allocate( edge_transport%model(1)%ggd( time_sind )%neutral( nneut ) )
            j = 0
            do js = 1, nspecies
               is = eb2spcr(js)
               if (.not.is_neutral(is)) cycle
               j = j + 1
               call species( is, spclabel, .false. )
               allocate( edge_profiles%ggd( time_sind )%neutral( j )%element(1) )
               allocate( edge_profiles%ggd( time_sind )%neutral( j )%state(1) )
#if ( IMAS_MAJOR_VERSION < 4 && IMAS_MINOR_VERSION < 42 )
               allocate( edge_profiles%ggd( time_sind )%neutral( j )%label(1) )
               allocate( edge_profiles%ggd( time_sind )%neutral( j )%state(1)%label(1) )
               edge_profiles%ggd( time_sind )%neutral( j )%label = species_list( js )
               edge_profiles%ggd( time_sind )%neutral( j )%state(1)%label = spclabel
#else
               allocate( edge_profiles%ggd( time_sind )%neutral( j )%name(1) )
               allocate( edge_profiles%ggd( time_sind )%neutral( j )%state(1)%name(1) )
               edge_profiles%ggd( time_sind )%neutral( j )%name = species_list( js )
               edge_profiles%ggd( time_sind )%neutral( j )%state(1)%name = spclabel
#endif
               edge_profiles%ggd( time_sind )%neutral( j )%element(1)%a = am( is )
#if IMAS_MAJOR_VERSION < 4
               edge_profiles%ggd( time_sind )%neutral( j )%element(1)%z_n = zn( is )
#else
               edge_profiles%ggd( time_sind )%neutral( j )%element(1)%z_n = nint(zn(is))
#endif
#if ( IMAS_MINOR_VERSION < 15 && IMAS_MAJOR_VERSION < 4 )
               edge_profiles%ggd( time_sind )%neutral( j )%element(1)%multiplicity = 1
#else
               edge_profiles%ggd( time_sind )%neutral( j )%element(1)%atoms_n = 1
#endif
               edge_profiles%ggd( time_sind )%neutral( j )%ion_index = js
               edge_profiles%ggd( time_sind )%neutral( j )%multiple_states_flag = 1
               allocate( edge_profiles%ggd( time_sind )%neutral( j )%state(1)% &
                   &      neutral_type%name(1) )
               allocate( edge_profiles%ggd( time_sind )%neutral( j )%state(1)% &
                   &      neutral_type%description(1) )
#if ( IMAS_MINOR_VERSION > 30 || IMAS_MAJOR_VERSION > 3 )
               edge_profiles%ggd( time_sind )%neutral( j )%state(1)%neutral_type%index = &
                   &      neutrals_identifier%thermal
               edge_profiles%ggd( time_sind )%neutral( j )%state(1)%neutral_type%name = &
                   &      neutrals_identifier%name( neutrals_identifier%thermal )
               edge_profiles%ggd( time_sind )%neutral( j )%state(1)%neutral_type%description = &
                   &      neutrals_identifier%description( neutrals_identifier%thermal )
#else
               edge_profiles%ggd( time_sind )%neutral( j )%state(1)%neutral_type%index = 2
               edge_profiles%ggd( time_sind )%neutral( j )%state(1)%neutral_type%name = &
                   &     "Thermal"
               edge_profiles%ggd( time_sind )%neutral( j )%state(1)%neutral_type%description = &
                   &     "Fluid neutral species from B2.5"
#endif
               do i = 1, nsources
                  allocate( edge_sources%source(i)%ggd( time_sind )%neutral( j )%element(1) )
                  allocate( edge_sources%source(i)%ggd( time_sind )%neutral( j )%state(1) )
#if ( IMAS_MAJOR_VERSION < 4 && IMAS_MINOR_VERSION < 42 )
                  allocate( edge_sources%source(i)%ggd( time_sind )%neutral( j )%label(1) )
                  allocate( edge_sources%source(i)%ggd( time_sind )%neutral( j )%state(1)%label(1) )
                  edge_sources%source(i)%ggd( time_sind )%neutral( j )%label = species_list( js )
                  edge_sources%source(i)%ggd( time_sind )%neutral( j )%state(1)%label = spclabel
#else
                  allocate( edge_sources%source(i)%ggd( time_sind )%neutral( j )%name(1) )
                  allocate( edge_sources%source(i)%ggd( time_sind )%neutral( j )%state(1)%name(1) )
                  edge_sources%source(i)%ggd( time_sind )%neutral( j )%name = species_list( js )
                  edge_sources%source(i)%ggd( time_sind )%neutral( j )%state(1)%name = spclabel
#endif
                  edge_sources%source(i)%ggd( time_sind )%neutral( j )%element(1)%a = am( is )
#if IMAS_MAJOR_VERSION < 4
                  edge_sources%source(i)%ggd( time_sind )%neutral( j )%element(1)%z_n = &
                     &   zn( is )
#else
                  edge_sources%source(i)%ggd( time_sind )%neutral( j )%element(1)%z_n = &
                     &   nint(zn(is))
#endif
#if ( IMAS_MINOR_VERSION < 15 && IMAS_MAJOR_VERSION < 4 )
                  edge_sources%source(i)%ggd( time_sind )%neutral( j )%element(1)%multiplicity = 1.0_IDS_real
#else
                  edge_sources%source(i)%ggd( time_sind )%neutral( j )%element(1)%atoms_n = 1
#endif
                  edge_sources%source(i)%ggd( time_sind )%neutral( j )%ion_index = js
                  edge_sources%source(i)%ggd( time_sind )%neutral( j )%multiple_states_flag = 1
                  allocate( edge_sources%source(i)%ggd( time_sind )%neutral( j )%state(1)% &
                     &      neutral_type%name(1) )
                  allocate( edge_sources%source(i)%ggd( time_sind )%neutral( j )%state(1)% &
                     &      neutral_type%description(1) )
#if ( IMAS_MINOR_VERSION > 30 || IMAS_MAJOR_VERSION > 3 )
                  edge_sources%source(i)%ggd( time_sind )%neutral( j )%state(1)%neutral_type%index = &
                     &      neutrals_identifier%thermal
                  edge_sources%source(i)%ggd( time_sind )%neutral( j )%state(1)%neutral_type%name = &
                     &      neutrals_identifier%name( neutrals_identifier%thermal )
                  edge_sources%source(i)%ggd( time_sind )%neutral( j )%state(1)%neutral_type%description = &
                     &      neutrals_identifier%description( neutrals_identifier%thermal )
#else
                  edge_sources%source(i)%ggd( time_sind )%neutral( j )%state(1)%neutral_type%index = 2
                  edge_sources%source(i)%ggd( time_sind )%neutral( j )%state(1)%neutral_type%name = &
                     &     "Thermal"
                  edge_sources%source(i)%ggd( time_sind )%neutral( j )%state(1)%neutral_type%description = &
                     &     "Fluid neutral species from B2.5"
#endif
               end do
               allocate( edge_transport%model(1)%ggd( time_sind )%neutral( j )%element(1) )
               allocate( edge_transport%model(1)%ggd( time_sind )%neutral( j )%state(1) )
#if ( IMAS_MAJOR_VERSION < 4 && IMAS_MINOR_VERSION < 42 )
               allocate( edge_transport%model(1)%ggd( time_sind )%neutral( j )%label(1) )
               allocate( edge_transport%model(1)%ggd( time_sind )%neutral( j )%state(1)%label(1) )
               edge_transport%model(1)%ggd( time_sind )%neutral( j )%label = species_list( js )
               edge_transport%model(1)%ggd( time_sind )%neutral( j )%state(1)%label = spclabel
#else
               allocate( edge_transport%model(1)%ggd( time_sind )%neutral( j )%name(1) )
               allocate( edge_transport%model(1)%ggd( time_sind )%neutral( j )%state(1)%name(1) )
               edge_transport%model(1)%ggd( time_sind )%neutral( j )%name = species_list( js )
               edge_transport%model(1)%ggd( time_sind )%neutral( j )%state(1)%name = spclabel
#endif
               edge_transport%model(1)%ggd( time_sind )%neutral( j )%element(1)%a = am( is )
#if IMAS_MAJOR_VERSION < 4
               edge_transport%model(1)%ggd( time_sind )%neutral( j )%element(1)%z_n = zn( is )
#else
               edge_transport%model(1)%ggd( time_sind )%neutral( j )%element(1)%z_n = nint(zn(is))
#endif
#if ( IMAS_MINOR_VERSION < 15 && IMAS_MAJOR_VERSION < 4 )
               edge_transport%model(1)%ggd( time_sind )%neutral( j )%element(1)%multiplicity = 1
#else
               edge_transport%model(1)%ggd( time_sind )%neutral( j )%element(1)%atoms_n = 1
#endif
               edge_transport%model(1)%ggd( time_sind )%neutral( j )%ion_index = js
               edge_transport%model(1)%ggd( time_sind )%neutral( j )%multiple_states_flag = 1
               allocate( edge_transport%model(1)%ggd( time_sind )%neutral( j )%state(1)% &
                   &      neutral_type%name(1) )
               allocate( edge_transport%model(1)%ggd( time_sind )%neutral( j )%state(1)% &
                   &      neutral_type%description(1) )
#if ( IMAS_MINOR_VERSION > 30 || IMAS_MAJOR_VERSION > 3 )
               edge_transport%model(1)%ggd( time_sind )%neutral( j )%state(1)%neutral_type%index = &
                   &      neutrals_identifier%thermal
               edge_transport%model(1)%ggd( time_sind )%neutral( j )%state(1)%neutral_type%name = &
                   &      neutrals_identifier%name( neutrals_identifier%thermal )
               edge_transport%model(1)%ggd( time_sind )%neutral( j )%state(1)%neutral_type%description = &
                   &      neutrals_identifier%description( neutrals_identifier%thermal )
#else
               edge_transport%model(1)%ggd( time_sind )%neutral( j )%state(1)%neutral_type%index = 2
               edge_transport%model(1)%ggd( time_sind )%neutral( j )%state(1)%neutral_type%name = &
                   &     "Thermal"
               edge_transport%model(1)%ggd( time_sind )%neutral( j )%state(1)%neutral_type%description = &
                   &     "Fluid neutral species from B2.5"
#endif
            end do
#if ( IMAS_MINOR_VERSION > 21 || IMAS_MAJOR_VERSION > 3 )
            allocate( radiation%process(1)%ggd( time_sind )%neutral( nneut ) )
            j = 0
            do js = 1, nspecies
               is = eb2spcr(js)
               if (.not.is_neutral(is)) cycle
               j = j + 1
               call species( is, spclabel, .false. )
               allocate( radiation%process(1)%ggd( time_sind )%neutral( j )%element(1) )
               allocate( radiation%process(1)%ggd( time_sind )%neutral( j )%state(1) )
#if ( IMAS_MAJOR_VERSION < 4 && IMAS_MINOR_VERSION < 42 )
               allocate( radiation%process(1)%ggd( time_sind )%neutral( j )%label(1) )
               allocate( radiation%process(1)%ggd( time_sind )%neutral( j )%state(1)%label(1) )
               radiation%process(1)%ggd( time_sind )%neutral( j )%label = species_list(js)
               radiation%process(1)%ggd( time_sind )%neutral( j )%state(1)%label = spclabel
#else
               allocate( radiation%process(1)%ggd( time_sind )%neutral( j )%name(1) )
               allocate( radiation%process(1)%ggd( time_sind )%neutral( j )%state(1)%name(1) )
               radiation%process(1)%ggd( time_sind )%neutral( j )%name = species_list(js)
               radiation%process(1)%ggd( time_sind )%neutral( j )%state(1)%name = spclabel
#endif
               radiation%process(1)%ggd( time_sind )%neutral( j )%element(1)%a = am(is)
#if IMAS_MAJOR_VERSION < 4
               radiation%process(1)%ggd( time_sind )%neutral( j )%element(1)%z_n = zn(is)
#else
               radiation%process(1)%ggd( time_sind )%neutral( j )%element(1)%z_n = nint(zn(is))
#endif
               radiation%process(1)%ggd( time_sind )%neutral( j )%element(1)%atoms_n = 1
               radiation%process(1)%ggd( time_sind )%neutral( j )%ion_index = js
               radiation%process(1)%ggd( time_sind )%neutral( j )%multiple_states_flag = 1
               allocate( radiation%process(1)%ggd( time_sind )%neutral( j )%state(1)% &
                   &     neutral_type%name(1) )
               allocate( radiation%process(1)%ggd( time_sind )%neutral( j )%state(1)% &
                   &     neutral_type%description(1) )
#if ( IMAS_MINOR_VERSION > 30 || IMAS_MAJOR_VERSION > 3 )
               radiation%process(1)%ggd( time_sind )%neutral( j )%state(1)%neutral_type%index = &
                   &     neutrals_identifier%thermal
               radiation%process(1)%ggd( time_sind )%neutral( j )%state(1)%neutral_type%name = &
                   &     neutrals_identifier%name( neutrals_identifier%thermal )
               radiation%process(1)%ggd( time_sind )%neutral( j )%state(1)%neutral_type%description = &
                   &     neutrals_identifier%description( neutrals_identifier%thermal )
#else
               radiation%process(1)%ggd( time_sind )%neutral( j )%state(1)%neutral_type%index = 2
               radiation%process(1)%ggd( time_sind )%neutral( j )%state(1)%neutral_type%name = &
                   &     "Thermal"
               radiation%process(1)%ggd( time_sind )%neutral( j )%state(1)%neutral_type%description = &
                   &     "Fluid neutral species from B2.5"
#endif
            end do
#endif
        end if

#if ( defined(B25_EIRENE) && ( IMAS_MINOR_VERSION > 21 || IMAS_MAJOR_VERSION > 3 ) )
        if (switch%use_eirene.ne.0) then
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
#if ( IMAS_MAJOR_VERSION < 4 && IMAS_MINOR_VERSION < 42 )
            allocate( radiation%process(3)%ggd( time_sind )%neutral( js )%label(1) )
            allocate( radiation%process(3)%ggd( time_sind )%neutral( js )%state( ks )%label(1) )
            radiation%process(3)%ggd( time_sind )%neutral( js )%label = &
                &    species_list(js)
            radiation%process(3)%ggd( time_sind )%neutral( js )%state( ks )%label = &
                &    textan( is-1 )
#else
            allocate( radiation%process(3)%ggd( time_sind )%neutral( js )%name(1) )
            allocate( radiation%process(3)%ggd( time_sind )%neutral( js )%state( ks )%name(1) )
            radiation%process(3)%ggd( time_sind )%neutral( js )%name = &
                &    species_list(js)
            radiation%process(3)%ggd( time_sind )%neutral( js )%state( ks )%name = &
                &    textan( is-1 )
#endif
            radiation%process(3)%ggd( time_sind )%neutral( js )%element(1)%a = nmassa( is )
#if IMAS_MAJOR_VERSION < 4
            radiation%process(3)%ggd( time_sind )%neutral( js )%element(1)%z_n = real(nchara(is),IDS_real)
#else
            radiation%process(3)%ggd( time_sind )%neutral( js )%element(1)%z_n = nchara( is )
#endif
            radiation%process(3)%ggd( time_sind )%neutral( js )%element(1)%atoms_n = 1
            radiation%process(3)%ggd( time_sind )%neutral( js )%ion_index = js
            radiation%process(3)%ggd( time_sind )%neutral( js )%multiple_states_flag = 1
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
#if ( IMAS_MAJOR_VERSION < 4 && IMAS_MINOR_VERSION < 42 )
            allocate( radiation%process(3)%ggd( time_sind )%neutral( js )%label(1) )
            allocate( radiation%process(3)%ggd( time_sind )%neutral( js )%state( ks )%label(1) )
#else
            allocate( radiation%process(3)%ggd( time_sind )%neutral( js )%name(1) )
            allocate( radiation%process(3)%ggd( time_sind )%neutral( js )%state( ks )%name(1) )
#endif
            i = 0
            do k = 1, natmi
              if (mlcmp( k, j ) > 0 ) then
                i = i + 1
                radiation%process(3)%ggd( time_sind )%neutral( js )%element(i)%a = nmassa( k )
#if IMAS_MAJOR_VERSION < 4
                radiation%process(3)%ggd( time_sind )%neutral( js )%element(i)%z_n = real(nchara(k),IDS_real)
#else
                radiation%process(3)%ggd( time_sind )%neutral( js )%element(i)%z_n = nchara( k )
#endif
                radiation%process(3)%ggd( time_sind )%neutral( js )%element(i)%atoms_n = &
                    &   mlcmp(k,j)
              end if
            end do
            ion_label = trim(textmn( j-1 ))//'+'
            i = 1
            match_found = .false.
            do while (.not.match_found.and.i.le.nioni)
              if (streql(ion_label,textin(i-1))) then
                match_found = .true.
              end if
              if (.not.match_found) i = i+1
            end do
            if (match_found) then
              k = nspecies+i
            else
              k = latmscl(lmolscl(j))
            end if
#if ( IMAS_MAJOR_VERSION < 4 && IMAS_MINOR_VERSION < 42 )
            radiation%process(3)%ggd( time_sind )%neutral( js )%label = textmn( j-1 )
            radiation%process(3)%ggd( time_sind )%neutral( js )%state( ks )%label = &
                &    textmn( j-1 )
#else
            radiation%process(3)%ggd( time_sind )%neutral( js )%name = textmn( j-1 )
            radiation%process(3)%ggd( time_sind )%neutral( js )%state( ks )%name = &
                &    textmn( j-1 )
#endif
            radiation%process(3)%ggd( time_sind )%neutral( js )%ion_index = k
            radiation%process(3)%ggd( time_sind )%neutral( js )%multiple_states_flag = 0
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
          do js = nspecies+1, nsion
            allocate( radiation%process(4)%ggd( time_sind )%ion( js )%state( istion(js) ) )
#if ( IMAS_MAJOR_VERSION < 4 && IMAS_MINOR_VERSION < 42 )
            allocate( radiation%process(4)%ggd( time_sind )%ion( js )%label(1) )
            do ks = 1, istion(js)
              is = ispion(js,ks)
              allocate( radiation%process(4)%ggd( time_sind )%ion( js )%state( ks )%label(1) )
              radiation%process(4)%ggd( time_sind )%ion( js )%state( ks )%label = textin( is-1 )
            end do
#else
            allocate( radiation%process(4)%ggd( time_sind )%ion( js )%name(1) )
            do ks = 1, istion(js)
              is = ispion(js,ks)
              allocate( radiation%process(4)%ggd( time_sind )%ion( js )%state( ks )%name(1) )
              radiation%process(4)%ggd( time_sind )%ion( js )%state( ks )%name = textin( is-1 )
            end do
#endif
            is = ispion(js,1)
            nelems = count ( micmp( 1:natmi, is ) > 0 )
            allocate( radiation%process(4)%ggd( time_sind )%ion( js )%element( nelems ) )
            i = 0
            do k = 1, natmi
              if ( micmp( k, is ) > 0 ) then
                i = i + 1
                radiation%process(4)%ggd( time_sind )%ion( js )%element( i )%a = nmassa(k)
#if IMAS_MAJOR_VERSION < 4
                radiation%process(4)%ggd( time_sind )%ion( js )%element( i )%z_n = real(nchara(k),IDS_real)
#else
                radiation%process(4)%ggd( time_sind )%ion( js )%element( i )%z_n = nchara(k)
#endif
                radiation%process(4)%ggd( time_sind )%ion( js )%element( i )%atoms_n = micmp(k,is)
              end if
            end do
            if (istion(js).eq.1) then
              radiation%process(4)%ggd( time_sind )%ion( js )%z_ion = nchrgi( is )
#if ( IMAS_MAJOR_VERSION < 4 && IMAS_MINOR_VERSION < 42 )
              radiation%process(4)%ggd( time_sind )%ion( js )%label = textin( is-1 )
#else
              radiation%process(4)%ggd( time_sind )%ion( js )%name = textin( is-1 )
#endif
            else
              match_found = .false.
              do ks = 2, istion(js)
                match_found = match_found .and. nchrgi(is).eq.nchrgi(ispion(js,ks))
              end do
              if (match_found) then
                radiation%process(4)%ggd( time_sind )%ion( js )%z_ion = nchrgi( is )
              end if
              match_found = .false.
              p = index(textin(ispion(js,1)-1),'+')
              if (p.gt.1) then
                ion_label = textin(ispion(js,1)-1)(1:p-1)
              else
                ion_label = textin(ispion(js,1)-1)
              end if
              do ks = 2, istion(js)
                p = index(textin(ispion(js,ks)-1),'+')
                if (p.gt.1) then
                  match_found = match_found .and. &
                    &  streql(ion_label,textin(ispion(js,ks)-1)(1:p-1))
                else
                  match_found = match_found .and. &
                    &  streql(ion_label,textin(ispion(js,ks)-1))
                end if
                if (match_found) then
#if ( IMAS_MAJOR_VERSION < 4 && IMAS_MINOR_VERSION < 42 )
                  radiation%process(4)%ggd( time_sind )%ion( js )%label = ion_label
#else
                  radiation%process(4)%ggd( time_sind )%ion( js )%name = ion_label
#endif
                end if
              end do
            end if
            is = ispion(js,1)
            ion_label = adjustl(textin(is-1))
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
            if (match_found) then
              jj = nspecies + j
            else
              jj = lkindi( is )
            end if
            radiation%process(4)%ggd( time_sind )%ion( js )%neutral_index = jj
            radiation%process(4)%ggd( time_sind )%ion( js )%multiple_states_flag = 1
            do ks = 1, istion(js)
              is = ispion(js,ks)
              radiation%process(4)%ggd( time_sind )%ion( js )%state( ks )%z_min = nchrgi( is )
              radiation%process(4)%ggd( time_sind )%ion( js )%state( ks )%z_max = nchrgi( is )
            end do
          end do

        end if
#endif

#ifdef B25_EIRENE
#ifdef WG_TODO
!! Obtain the neutral velocities
!! Recall that P[XYZ]DEN[AM] are momentum densities in CGS units!
!! P[UV][XY] will only exist if fort.46 file was written with format 20170930 or later
        if (switch%use_eirene.ne.0 .and. allocated(pux)) then
          allocate(un0(mpg%nCv,0:2,natmi))
          allocate(um0(mpg%nCv,0:2,nmoli))
          un0 = 0.0_IDS_real
          um0 = 0.0_IDS_real
          DO IS = 1, NATMI
            totCv = 0.0_IDS_real
            do I = 1, NTRII
              IX = IXTRI(I)
              IY = IYTRI(I)
              IF (B2_CELL(IX,IY)) THEN
                IXX = ix_e2b(IX)
                IYY = IY-1
                UN0(IXX,IYY,0,IS) = UN0(IXX,IYY,0,IS) + &
                   &  TRIANGLE_VOL(I)*( &
                   &   VXDENA(IS,I)*PUX(I) + VYDENA(IS,I)*PUY(I) )
                UN0(IXX,IYY,1,IS) = UN0(IXX,IYY,1,IS) + &
                   &  TRIANGLE_VOL(I)*( &
                   &   VXDENA(IS,I)*PVX(I) + VYDENA(IS,I)*PVY(I) )
                UN0(IXX,IYY,2,IS) = UN0(IXX,IYY,2,IS) + &
                   &  TRIANGLE_VOL(I)*VZDENA(IS,I)
                totCv(IXX,IYY) = totCv(IXX,IYY) + &
                   &  TRIANGLE_VOL(I)*PDENA(IS,I)
              END IF
            end do
            do iy = -1, ny
              do ix = -1, nx
                if (totCv(ix,iy).gt.0.0_IDS_real) then
                  un0(ix,iy,:,is) = un0(ix,iy,:,is) / totCv(ix,iy) &
                     &  / (nmassa(is)*mp*1000.0_R8) / 100.0_R8
                end if
              end do
            end do
          END DO
          DO IS = 1, NMOLI
            totCv = 0.0_IDS_real
            do I = 1, NTRII
              IX = IXTRI(I)
              IY = IYTRI(I)
              IF (B2_CELL(IX,IY)) THEN
                IXX = ix_e2b(IX)
                IYY = IY-1
                UM0(IXX,IYY,0,IS) = UM0(IXX,IYY,0,IS) + &
                   &  TRIANGLE_VOL(I)*( &
                   &   VXDENM(IS,I)*PUX(I) + VYDENM(IS,I)*PUY(I) )
                UM0(IXX,IYY,1,IS) = UM0(IXX,IYY,1,IS) + &
                   &  TRIANGLE_VOL(I)*( &
                   &   VXDENM(IS,I)*PVX(I) + VYDENM(IS,I)*PVY(I) )
                UM0(IXX,IYY,2,IS) = UM0(IXX,IYY,2,IS) + &
                   &  TRIANGLE_VOL(I)*VZDENM(IS,I)
                totCv(IXX,IYY) = totCv(IXX,IYY) + &
                   &  TRIANGLE_VOL(I)*PDENM(IS,I)
              END IF
            end do
            do iy = -1, ny
              do ix = -1, nx
                if (totCv(ix,iy).gt.0.0_IDS_real) then
                  um0(ix,iy,:,is) = um0(ix,iy,:,is) / totCv(ix,iy) &
                     &  / (nmassm(is)*mp*1000.0_R8) / 100.0_R8
                end if
              end do
            end do
          END DO
        end if
#endif
#endif
#if ( IMAS_MINOR_VERSION > 21 || IMAS_MAJOR_VERSION > 3 )
        call write_sourced_string( summary%configuration, configuration )
#endif

        !! Write plasma state
        if ( B2_WRITE_DATA ) then
#if ( IMAS_MINOR_VERSION > 11 || IMAS_MAJOR_VERSION > 3 ) && GGD_MAJOR_VERSION > 0
            call logmsg( LOGDEBUG, &
            &   "b2mod_ual_io.B25_process_ids: writing plasma state" )

            !! Find grid subset base indices out of the available grid subset
            !! data stored in the IDS. That is done using IMAS GGD routine
            !! findGridSubsetByName().
#if ( IMAS_MINOR_VERSION < 15 && IMAS_MAJOR_VERSION < 4 )
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
            if (geo%LSN) then
              iGsIDivertor = findGridSubsetByName( edge_profiles% &
                &   ggd( time_sind )%grid, "Western divertor" )
              iGsODivertor = findGridSubsetByName( edge_profiles% &
                &   ggd( time_sind )%grid, "Eastern divertor" )
            else
              iGsIDivertor = findGridSubsetByName( edge_profiles% &
                &   ggd( time_sind )%grid, "Eastern divertor" )
              iGsODivertor = findGridSubsetByName( edge_profiles% &
                &   ggd( time_sind )%grid, "Western divertor" )
            end if
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
            if (geo%LSN) then
              iGsIDivertor = findGridSubsetByName(      &
                &   edge_profiles%grid_ggd( time_sind ), "Western divertor" )
              iGsODivertor = findGridSubsetByName(      &
                &   edge_profiles%grid_ggd( time_sind ), "Eastern divertor" )
            else
              iGsIDivertor = findGridSubsetByName(      &
                &   edge_profiles%grid_ggd( time_sind ), "Eastern divertor" )
              iGsODivertor = findGridSubsetByName(      &
                &   edge_profiles%grid_ggd( time_sind ), "Western divertor" )
            end if
#endif
#if ( IMAS_MINOR_VERSION < 15 && IMAS_MAJOR_VERSION < 4 )
            edge_grid = edge_profiles%ggd( time_sind )%grid
            transport_grid = edge_transport%ggd( time_sind )%grid
            sources_grid = edge_sources%ggd( time_sind )%grid
#else
#if AL_MAJOR_VERSION < 5
            edge_grid = edge_profiles%grid_ggd( time_sind )
            transport_grid = edge_transport%grid_ggd( time_sind )
            sources_grid = edge_sources%grid_ggd( time_sind )
#if ( IMAS_MINOR_VERSION > 21 || IMAS_MAJOR_VERSION > 3 )
            radiation_grid = radiation%grid_ggd( time_sind )
#endif
#else
            edge_grid = edge_profiles%grid_ggd( time_sind )
            transport_grid = edge_profiles%grid_ggd( time_sind )
            sources_grid = edge_profiles%grid_ggd( time_sind )
#if ( IMAS_MINOR_VERSION > 21 || IMAS_MAJOR_VERSION > 3 )
            radiation_grid = edge_profiles%grid_ggd( time_sind )
#endif
#endif
#endif
            !! ne: Electron density
            call write_IDS_quantity( edge_grid, mpg, geo,                   &
                &   val = edge_profiles%ggd( time_sind )%electrons%density, &
                &   value = state%dv%ne )
            !! fne: Electron particle flux
            call divide_by_area( mpg%nFc, geo, state%dv%fne, flxFace)
#if ( IMAS_MINOR_VERSION > 38 || IMAS_MAJOR_VERSION > 3 )
            call write_face_scalar( transport_grid, mpg,                    &
                &   val = edge_transport%model(1)%ggd( time_sind )%         &
                &         electrons%particles%flux_pol,                     &
                &   value = flxFace(:,0) )
            call write_face_scalar( transport_grid, mpg,                    &
                &   val = edge_transport%model(1)%ggd( time_sind )%         &
                &         electrons%particles%flux_radial,                  &
                &   value = flxFace(:,1) )
#else
            call split_face_flux( mpg, geo, flxFace, tmpFace )
            call write_face_scalar( transport_grid, mpg,                    &
                &   val = edge_transport%model(1)%ggd( time_sind )%         &
                &         electrons%particles%flux,                         &
                &   value = tmpFace )
#endif
            !! sne: Electron particle sources
            tmpCv(:) = ( state%sr%sne(:,0) + state%sr%sne(:,1) * state%dv%ne(:) ) / geo%cvVol(:)
            call write_cell_scalar( sources_grid, mpg,                      &
                &   scalar = edge_sources%source(1)%ggd( time_sind )%       &
                &            electrons%particles,                           &
                &   b2CellData = tmpCv )
            tmpCv(:) = state_ext%sne(:) / geo%cvVol(:)
            call write_cell_scalar( sources_grid, mpg,                      &
                &   scalar = edge_sources%source(2)%ggd( time_sind )%       &
                &            electrons%particles,                           &
                &   b2CellData = tmpCv )
            tmpCv(:) = ( state%srw%b2stbc_sne(:) + state%srw%b2stbm_sne(:) ) / geo%cvVol(:)
            call write_cell_scalar( sources_grid, mpg,                      &
                &   scalar = edge_sources%source(3)%ggd( time_sind )%       &
                &            electrons%particles,                           &
                &   b2CellData = tmpCv )
            tmpCv(:) = ( state%sr%snedt(:,0) +                              &
                &        state%sr%snedt(:,1) * state%dv%ne(:) ) / geo%cvVol(:)
            call write_cell_scalar( sources_grid, mpg,                      &
                &   scalar = edge_sources%source(4)%ggd( time_sind )%       &
                &            electrons%particles,                           &
                &   b2CellData = tmpCv )
            if (switch%use_eirene.ne.0 .and. balance_netcdf.ne.0) then
                tmpCv = 0.0_IDS_real
                do istrai = 1, size( eirene_mc_pael_sne_bal, 2)
                    tmpCv(:) = tmpCv(:) + eirene_mc_pael_sne_bal(:,istrai)
                end do
                tmpCv(:) = tmpCv(:) / geo%cvVol(:)
                call write_cell_scalar( sources_grid, mpg,                  &
                    &   scalar = edge_sources%source(5)%ggd( time_sind )%   &
                    &            electrons%particles,                       &
                    &   b2CellData = tmpCv )
                tmpCv = 0.0_IDS_real
                do istrai = 1, size( eirene_mc_pmel_sne_bal, 2)
                    tmpCv(:) = tmpCv(:) + eirene_mc_pmel_sne_bal(:,istrai)
                end do
                tmpCv(:) = tmpCv(:) / geo%cvVol(:)
                call write_cell_scalar( sources_grid, mpg,                  &
                    &   scalar = edge_sources%source(6)%ggd( time_sind )%   &
                    &            electrons%particles,                       &
                    &   b2CellData = tmpCv )
            else
                tmpCv(:) = state%srw%b2stbr_sne(:) / geo%cvVol(:)
                call write_cell_scalar( sources_grid, mpg,                  &
                    &   scalar = edge_sources%source(5)%ggd( time_sind )%   &
                    &            electrons%particles,                       &
                    &   b2CellData = tmpCv )
            end if
            tmpCv = 0.0_IDS_real
            do is = 0, ns-1
              tmpCv(:) = tmpCv(:) + state%srw%rsana(:,is)
            end do
            tmpCv(:) = tmpCv(:) / geo%cvVol(:)
            call write_cell_scalar( sources_grid, mpg,                      &
                &   scalar = edge_sources%source(7)%ggd( time_sind )%       &
                &            electrons%particles,                           &
                &   b2CellData = tmpCv )
            tmpCv = 0.0_IDS_real
            do is = 0, ns-1
              tmpCv(:) = tmpCv(:) + state%srw%rrana(:,is)
            end do
            tmpCv(:) = tmpCv(:) / geo%cvVol(:)
            call write_cell_scalar( sources_grid, mpg,                      &
                &   scalar = edge_sources%source(8)%ggd( time_sind )%       &
                &            electrons%particles,                           &
                &   b2CellData = tmpCv )
#if ( ( IMAS_MINOR_VERSION > 38 || IMAS_MAJOR_VERSION > 3 ) && defined(B25_EIRENE) )
            if (switch%use_eirene.ne.0) then
              tmpCv(:) = state%srw%sne0_eir_tot(:) / geo%cvVol(:)
              call write_cell_scalar( sources_grid, mpg,                    &
                &   scalar = edge_sources%source(13)%ggd( time_sind )%      &
                &            electrons%particles,                           &
                &   b2CellData = tmpCv )
            end if
#endif

            !! na: Ion density
            do is = 1, nsion
              if (is.le.nspecies) then
                tmpCv(:) = 0.0_IDS_real
                do js = 1, istion(is)
                  call write_IDS_quantity( edge_grid, mpg, geo,             &
                      &   val = edge_profiles%ggd( time_sind )%             &
                      &         ion( is )%state( js )%density,              &
                      &   value = state%pl%na(:,ispion(is,js)) )
                  tmpCv(:) = tmpCv(:) + state%pl%na(:,ispion(is,js))
                end do
                call write_IDS_quantity( edge_grid, mpg, geo,               &
                    &   val = edge_profiles%ggd( time_sind )%               &
                    &         ion(is)%density,                              &
                    &   value = tmpCv )
            !! fna: Ion particle flux
                totflux(:,0:1) = 0.0_IDS_real
                do js = 1, istion(is)
                  call divide_by_area( mpg%nFc, geo, state%dv%fna(:,:,ispion(is,js)), flxFace)
                  totflux(:,:) = totflux(:,:) + flxFace(:,:)
#if ( IMAS_MINOR_VERSION > 38 || IMAS_MAJOR_VERSION > 3 )
                  call write_face_scalar( transport_grid, mpg,               &
                      &   val = edge_transport%model(1)%ggd( time_sind )%    &
                      &         ion( is )%state( js )%particles%flux_pol,    &
                      &   value = flxFace(:,0) )
                  call write_face_scalar( transport_grid, mpg,               &
                      &   val = edge_transport%model(1)%ggd( time_sind )%    &
                      &         ion( is )%state( js )%particles%flux_radial, &
                      &   value = flxFace(:,1) )
#else
                  call split_face_flux( mpg, geo, flxFace, tmpFace )
                  call write_face_scalar( transport_grid, mpg,               &
                      &   val = edge_transport%model(1)%ggd( time_sind )%    &
                      &         ion( is )%state( js )%particles%flux,        &
                      &   value = tmpFace )
#endif
                end do
#if ( IMAS_MINOR_VERSION > 38 || IMAS_MAJOR_VERSION > 3 )
                call write_face_scalar( transport_grid, mpg,                 &
                    &   val = edge_transport%model(1)%ggd( time_sind )%      &
                    &         ion( is )%particles%flux_pol,                  &
                    &   value = totflux(:,0) )
                call write_face_scalar( transport_grid, mpg,                 &
                    &   val = edge_transport%model(1)%ggd( time_sind )%      &
                    &         ion( is )%particles%flux_radial,               &
                    &   value = totflux(:,1) )
#else
                call split_face_flux( mpg, geo, totflux, tmpFace )
                call write_face_scalar( transport_grid, mpg,                 &
                    &   val = edge_transport%model(1)%ggd( time_sind )%      &
                    &         ion( is )%particles%flux,                      &
                    &   value = tmpFace )
#endif
            !! cdna: Ion diffusivity
                do js = 1, istion(is)
#if ( IMAS_MINOR_VERSION > 38 || IMAS_MAJOR_VERSION > 3 )
                  call write_face_scalar( transport_grid, mpg,               &
                      &   val = edge_transport%model(1)%ggd( time_sind )%    &
                      &         ion( is )%state( js )%particles%d_pol,       &
                      &   value = state%co%cdna(:,0,ispion(is,js)) )
                  call write_face_scalar( transport_grid, mpg,               &
                      &   val = edge_transport%model(1)%ggd( time_sind )%    &
                      &         ion( is )%state( js )%particles%d_radial,    &
                      &   value = state%co%cdna(:,1,ispion(is,js)) )
            !! cvla: Ion diffusivity
                  call write_face_scalar( transport_grid, mpg,               &
                      &   val = edge_transport%model(1)%ggd( time_sind )%    &
                      &         ion( is )%state( js )%particles%v_pol,       &
                      &   value = state%co%cvla(:,0,ispion(is,js)) )
                  call write_face_scalar( transport_grid, mpg,               &
                      &   val = edge_transport%model(1)%ggd( time_sind )%    &
                      &         ion( is )%state( js )%particles%v_radial,    &
                      &   value = state%co%cvla(:,1,ispion(is,js)) )
#else
                  call split_face_flux( mpg, geo,                            &
                      &   state%co%cdna(:,:,ispion(is,js)), tmpFace )
                  call write_face_scalar( transport_grid, mpg,               &
                      &   val = edge_transport%model(1)%ggd( time_sind )%    &
                      &         ion( is )%state( js )%particles%d,           &
                      &   value = tmpFace )
                  call split_face_flux( mpg, geo,                            &
                      &   state%co%cvla(:,:,ispion(is,js)), tmpFace )
                  call write_face_scalar( transport_grid, mpg,               &
                      &   val = edge_transport%model(1)%ggd( time_sind )%    &
                      &         ion( is )%state( js )%particles%v,           &
                      &   value = tmpFace )
#endif
            !! fllim0fna: Ion flux limiter
                  call split_face_flux( mpg, geo,                            &
                      &   state%co%fllim0fna(:,:,ispion(is,js)), tmpFace )
                  call write_face_scalar( transport_grid, mpg,               &
                      &   val = edge_transport%model(1)%ggd( time_sind )%    &
                      &         ion( is )%state( js )%                       &
                      &         particles%flux_limiter,                      &
                      &   value = tmpFace )
                end do
            !! sna: Ion particle sources
                totCv(:) = 0.0_IDS_real
                do js = 1, istion(is)
                  tmpCv(:) = ( state%sr%sna(:,0,ispion(is,js)) +            &
                      &        state%sr%sna(:,1,ispion(is,js)) *            &
                      &        state%pl%na(:,ispion(is,js)) ) / geo%cvVol(:)
                  totCv(:) = totCv(:) + tmpCv(:)
                  call write_cell_scalar( sources_grid, mpg,                &
                      &   scalar = edge_sources%source(1)%                  &
                      &   ggd( time_sind )%ion( is )%state( js )%particles, &
                      &   b2CellData = tmpCv )
                end do
                call write_cell_scalar( sources_grid, mpg,                  &
                      &   scalar = edge_sources%source(1)%                  &
                      &   ggd( time_sind )%ion( is )%particles,             &
                      &   b2CellData = totCv )
                totCv(:) = 0.0_IDS_real
                do js = 1, istion(is)
                  tmpCv(:) = state_ext%sna(:,ispion(is,js)) / geo%cvVol(:)
                  totCv(:) = totCv(:) + tmpCv(:)
                  call write_cell_scalar( sources_grid, mpg,                &
                      &   scalar = edge_sources%source(2)%ggd( time_sind )% &
                      &            ion( is )%state( js )%particles,         &
                      &   b2CellData = tmpCv )
                end do
                call write_cell_scalar( sources_grid, mpg,                  &
                    &   scalar = edge_sources%source(2)%ggd( time_sind )%   &
                    &            ion( is )%particles,                       &
                    &   b2CellData = totCv )
                totCv(:) = 0.0_IDS_real
                do js = 1, istion(is)
                  tmpCv(:) = ( state%srw%b2stbc_sna(:,ispion(is,js)) +      &
                     &         state%srw%b2stbm_sna(:,ispion(is,js)) ) / geo%cvVol(:)
                  totCv(:) = totCv(:) + tmpCv(:)
                  call write_cell_scalar( sources_grid, mpg,                &
                      &   scalar = edge_sources%source(3)%ggd( time_sind )% &
                      &            ion( is )%state( js )%particles,         &
                      &   b2CellData = tmpCv )
                end do
                call write_cell_scalar( sources_grid, mpg,                  &
                    &   scalar = edge_sources%source(3)%ggd( time_sind )%   &
                    &            ion( is )%particles,                       &
                    &   b2CellData = totCv )
                totCv(:) = 0.0_IDS_real
                do js = 1, istion(is)
                  tmpCv(:) = ( state%sr%snadt(:,0,ispion(is,js)) +          &
                        &      state%sr%snadt(:,1,ispion(is,js)) *          &
                        &      state%pl%na(:,ispion(is,js)) ) / geo%cvVol(:)
                  totCv(:) = totCv(:) + tmpCv(:)
                  call write_cell_scalar( sources_grid, mpg,                &
                      &   scalar = edge_sources%source(4)%ggd( time_sind )% &
                      &            ion( is )%state( js )%particles,         &
                      &   b2CellData = tmpCv )
                end do
                call write_cell_scalar( sources_grid, mpg,                  &
                    &   scalar = edge_sources%source(4)%ggd( time_sind )%   &
                    &            ion( is )%particles,                       &
                    &   b2CellData = totCv )
                if (switch%use_eirene.ne.0 .and. balance_netcdf.ne.0) then
                  totCv(:) = 0.0_IDS_real
                  do js = 1, istion(is)
                    tmpCv = 0.0_IDS_real
                    do istrai = 1, size( eirene_mc_papl_sna_bal, 3)
                      tmpCv(:) = tmpCv(:)                                   &
                         &       + eirene_mc_papl_sna_bal(:,ispion(is,js),istrai)
                    end do
                    tmpCv(:) = tmpCv(:) / geo%cvVol(:)
                    totCv(:) = totCv(:) + tmpCv(:)
                    call write_cell_scalar( sources_grid, mpg,                &
                        &   scalar = edge_sources%source(5)%ggd( time_sind )% &
                        &            ion( is )%state( js )%particles,         &
                        &   b2CellData = tmpCv )
                  end do
                  call write_cell_scalar( sources_grid, mpg,                  &
                        &   scalar = edge_sources%source(5)%ggd( time_sind )% &
                        &            ion( is )%particles,                     &
                        &   b2CellData = totCv )
                  totCv = 0.0_IDS_real
                  do js = 1, istion(is)
                    tmpCv = 0.0_IDS_real
                    do istrai = 1, size( eirene_mc_pmpl_sna_bal, 3)
                      tmpCv(:) = tmpCv(:)                                     &
                         &       + eirene_mc_pmpl_sna_bal(:,ispion(is,js),istrai)
                    end do
                    tmpCv(:) = tmpCv(:) / geo%cvVol(:)
                    totCv(:) = totCv(:) + tmpCv(:)
                    call write_cell_scalar( sources_grid, mpg,                &
                        &   scalar = edge_sources%source(6)%ggd( time_sind )% &
                        &            ion( is )%state( js )%particles,         &
                        &   b2CellData = tmpCv )
                  end do
                  call write_cell_scalar( sources_grid, mpg,                  &
                      &   scalar = edge_sources%source(6)%ggd( time_sind )%   &
                      &            ion( is )%particles,                       &
                      &   b2CellData = totCv )
                else
                  totCv(:) = 0.0_IDS_real
                  do js = 1, istion(is)
                    tmpCv(:) = state%srw%b2stbr_sna(:,ispion(is,js)) / geo%cvVol(:)
                    totCv(:) = totCv(:) + tmpCv(:)
                    call write_cell_scalar( sources_grid, mpg,                &
                        &   scalar = edge_sources%source(5)%ggd( time_sind )% &
                        &            ion( is )%state( js )%particles,         &
                        &   b2CellData = tmpCv )
                  end do
                  call write_cell_scalar( sources_grid, mpg,                  &
                        &   scalar = edge_sources%source(5)%ggd( time_sind )% &
                        &            ion( is )%particles,                     &
                        &   b2CellData = totCv )
                end if
                do js = 1, istion(is)
                  tmpCv(:) = state%srw%rsana(:,ispion(is,js)) / geo%cvVol(:)
                  call write_cell_scalar( sources_grid, mpg,                  &
                      &   scalar = edge_sources%source(7)%ggd( time_sind )%   &
                      &            ion( is )%state( js )%particles,           &
                      &   b2CellData = tmpCv )
                  tmpCv(:) = state%srw%rrana(:,ispion(is,js)) / geo%cvVol(:)
                  call write_cell_scalar( sources_grid, mpg,                  &
                      &   scalar = edge_sources%source(8)%ggd( time_sind )%   &
                      &            ion( is )%state( js )%particles,           &
                      &   b2CellData = tmpCv )
                  tmpCv(:) = state%srw%rcxna(:,ispion(is,js)) / geo%cvVol(:)
                  call write_cell_scalar( sources_grid, mpg,                  &
                      &   scalar = edge_sources%source(9)%ggd( time_sind )%   &
                      &            ion( is )%state( js )%particles,           &
                      &   b2CellData = tmpCv )
                end do
#if ( ( IMAS_MINOR_VERSION > 38 || IMAS_MAJOR_VERSION > 3 ) && defined(B25_EIRENE) )
                if (switch%use_eirene.ne.0) then
                  do js = 1, istion(is)
                    tmpCv(:) = state%srw%sna0_eir_tot(:,ispion(is,js)) /      &
                      &   geo%cvVol(:)
                    call write_cell_scalar( sources_grid, mpg,                &
                      &   scalar = edge_sources%source(13)%ggd( time_sind )%  &
                      &            ion( is )%state( js )%particles,           &
                      &   b2CellData = tmpCv )
                  end do
                end if
#endif
              else
                totCv(:) = 0.0_IDS_real
                do js = 1, istion(is)
                  tmpCv(:) = dib2(:,ispion(is,js),1)
                  totCv(:) = totCv(:) + tmpCv(:)
                  call write_IDS_quantity( edge_grid, mpg, geo,             &
                      &   val = edge_profiles%ggd( time_sind )%             &
                      &         ion(is)%state( js )%density,                &
                      &   value = tmpCv )
                end do
                call write_IDS_quantity( edge_grid, mpg, geo,               &
                    &   val = edge_profiles%ggd( time_sind )%               &
                    &         ion(is)%density,                              &
                    &   value = totCv )
                if (balance_netcdf.ne.0) then
                  totCv = 0.0_IDS_real
                  do js = 1, istion(is)
                    tmpCv = 0.0_IDS_real
                    do istrai = 1, size(eirene_mc_paio_sna_bal,3)
                      tmpCv(:) = tmpCv(:)                                         &
                         &       + eirene_mc_paio_sna_bal(:,ispion(is,js),istrai) &
                         &       + eirene_mc_pmio_sna_bal(:,ispion(is,js),istrai) &
                         &       + eirene_mc_piio_sna_bal(:,ispion(is,js),istrai)
                    end do
                    tmpCv(:) = tmpCv(:) / geo%cvVol(:)
                    totCv(:) = totCv(:) + tmpCv(:)
                    call write_cell_scalar( sources_grid, mpg,                &
                        &   scalar = edge_sources%source(1)%ggd( time_sind )% &
                        &            ion( is )%state( js )%particles,         &
                        &   b2CellData = tmpCv )
                  end do
                  call write_cell_scalar( sources_grid, mpg,                  &
                      &   scalar = edge_sources%source(1)%ggd( time_sind )%   &
                      &            ion( is )%particles,                       &
                      &   b2CellData = totCv )
                end if
              end if
            end do

            !! ue: Parallel Electron Velocity
            call write_cell_vector_component( edge_grid, mpg,            &
                 &   vectorComponent = edge_profiles%ggd( time_sind )%   &
                 &                     electrons%velocity,               &
                 &   b2CellData = state%dv%ue(:),                        &
                 &   vectorID = VEC_ALIGN_PARALLEL_ID )

            do is = 1, nsion
              if (is.le.nspecies) then
                do js = 1, istion(is)
                !! ua: Parallel ion velocity
                  call write_cell_vector_component( edge_grid, mpg,         &
                      &   vectorComponent = edge_profiles%ggd( time_sind )% &
                      &                     ion( is )%state( js )%velocity, &
                      &   b2CellData = state%pl%ua(:,ispion(is,js)),        &
                      &   vectorID = VEC_ALIGN_PARALLEL_ID )
                !! wadia: Diamagnetic ion velocity
                  call write_face_vector_component( edge_grid, mpg,           &
                      &   vectorComponent = edge_profiles%ggd( time_sind )%   &
                      &         ion( is )%state( js )%velocity_diamagnetic,   &
                      &   b2FaceData = state%dv%wadia(:,0,ispion(is,js)),     &
                      &   vectorID = VEC_ALIGN_POLOIDAL_ID )
                  call write_face_vector_component( edge_grid, mpg,           &
                      &   vectorComponent = edge_profiles%ggd( time_sind )%   &
                      &         ion( is )%state( js )%velocity_diamagnetic,   &
                      &   b2FaceData = state%dv%wadia(:,1,ispion(is,js)),     &
                      &   vectorID = VEC_ALIGN_RADIAL_ID )
                !! vaecrb: ExB ion velocity
                  call write_face_vector_component( edge_grid, mpg,           &
                      &   vectorComponent = edge_profiles%ggd( time_sind )%   &
                      &         ion( is )%state( js )%velocity_exb,           &
                      &   b2FaceData = state%dv%vaecrb(:,0,ispion(is,js)),    &
                      &   vectorID = VEC_ALIGN_POLOIDAL_ID )
                  call write_face_vector_component( edge_grid, mpg,           &
                      &   vectorComponent = edge_profiles%ggd( time_sind )%   &
                      &         ion( is )%state( js )%velocity_exb,           &
                      &   b2FaceData = state%dv%vaecrb(:,1,ispion(is,js)),    &
                      &   vectorID = VEC_ALIGN_RADIAL_ID )
                !! cvsa: Ion diffusivity
                  call write_face_vector_component( transport_grid, mpg,      &
                      &   vectorComponent = edge_transport%model(1)%          &
                      &         ggd( time_sind )%ion( is )%state( js )%       &
                      &         momentum%d,                                   &
                      &   b2FaceData = state%co%cvsa(:,0,ispion(is,js)),      &
                      &   vectorID = VEC_ALIGN_POLOIDAL_ID )
                  call write_face_vector_component( transport_grid, mpg,      &
                      &   vectorComponent = edge_transport%model(1)%          &
                      &         ggd( time_sind )%ion( is )%state( js )%       &
                      &         momentum%d,                                   &
                      &   b2FaceData = state%co%cvsa(:,0,ispion(is,js)),      &
                      &   vectorID = VEC_ALIGN_RADIAL_ID )
                end do
                !! fmo: Ion momentum flux
                totflux(:,:) = 0.0_IDS_real
                do js = 1, istion(is)
                  call divide_by_area( mpg%nFc, geo, state%dv%fmo(1,0,ispion(is,js)), flxFace)
                  call write_face_vector_component( transport_grid, mpg,      &
                      &   vectorComponent = edge_transport%model(1)%          &
                      &                     ggd( time_sind )%ion( is )%       &
                      &                     state( js )%momentum%flux,        &
                      &   b2FaceData = flxFace,                               &
                      &   vectorID = VEC_ALIGN_PARALLEL_ID )
                  totflux(:,0) = totflux(:,0) + flxFace(:,0)
                  totflux(:,1) = totflux(:,1) + flxFace(:,1)
                end do
                call write_face_vector_component( transport_grid, mpg,        &
                    &   vectorComponent = edge_transport%model(1)%            &
                    &                     ggd( time_sind )%ion( is )%         &
                    &                     momentum%flux,                      &
                    &   b2FaceData = totflux,                                 &
                    &   vectorID = VEC_ALIGN_PARALLEL_ID )
                !! fllimvisc: Ion parallel momentum transport flux limit
                do js = 1, istion(is)
                  call write_face_vector_component( transport_grid, mpg,      &
                      &   vectorComponent = edge_transport%model(1)%          &
                      &                     ggd( time_sind )%ion( is )%       &
                      &                     state( js )%momentum%flux_limiter,&
                      &   b2FaceData = state%co%fllimvisc(:,ispion(is,js)),   &
                      &   vectorID = VEC_ALIGN_PARALLEL_ID )
                end do
                !! smo: Ion parallel momentum sources
                totCv(:) = 0.0_IDS_real
                do js = 1, istion(is)
                  do iCv = 1, mpg%nCv
                    tmpCv(iCv) = ( state%sr%smo(iCv,0,ispion(is,js)) +        &
                      &            state%sr%smo(iCv,1,ispion(is,js)) *        &
                      &            state%pl%ua(iCv,ispion(is,js)) +           &
                      &            state%sr%smo(iCv,2,ispion(is,js)) *        &
                      &                    roxa(iCv,ispion(is,js)) +          &
                      &            state%sr%smo(iCv,3,ispion(is,js)) *        &
                      &                    roxa(iCv,ispion(is,js)) *          &
                      &            state%pl%ua(iCv,ispion(is,js)) ) / geo%cvVol(iCv)
                  end do
                  totCv(:) = totCv(:) + tmpCv(:)
                  call write_cell_vector_component( sources_grid, mpg,  &
                      &   vectorComponent = edge_sources%source(1)%     &
                      &                     ggd( time_sind )%ion( is )% &
                      &                     state( js )%momentum,       &
                      &   b2CellData = tmpCv,                           &
                      &   vectorID = VEC_ALIGN_PARALLEL_ID )
                end do
                call write_cell_vector_component( sources_grid, mpg,           &
                    &   vectorComponent = edge_sources%source(1)%              &
                    &                     ggd( time_sind )%ion( is )%momentum, &
                    &   b2CellData = totCv,                                    &
                    &   vectorID = VEC_ALIGN_PARALLEL_ID )
                totCv(:) = 0.0_IDS_real
                do js = 1, istion(is)
                  tmpCv(:) = state_ext%smo(:,ispion(is,js)) / geo%cvVol(:)
                  totCv(:) = totCv(:) + tmpCv(:)
                  call write_cell_vector_component( sources_grid, mpg,  &
                      &   vectorComponent = edge_sources%source(2)%     &
                      &                     ggd( time_sind )%ion( is )% &
                      &                     state( js )%momentum,       &
                      &   b2CellData = tmpCv,                           &
                      &   vectorID = VEC_ALIGN_PARALLEL_ID )
                end do
                call write_cell_vector_component( sources_grid, mpg,    &
                      &   vectorComponent = edge_sources%source(2)%     &
                      &                     ggd( time_sind )%ion( is )% &
                      &                     momentum,                   &
                      &   b2CellData = totCv,                           &
                      &   vectorID = VEC_ALIGN_PARALLEL_ID )
                totCv(:) = 0.0_IDS_real
                do js = 1, istion(is)
                  tmpCv(:) = ( state%srw%b2stbc_smo(:,ispion(is,js)) +  &
                      &        state%srw%b2stbm_smo(:,ispion(is,js)) ) / geo%cvVol(:)
                  totCv(:) = totCv(:) + tmpCv(:)
                  call write_cell_vector_component( sources_grid, mpg,  &
                      &   vectorComponent = edge_sources%source(3)%     &
                      &                     ggd( time_sind )%ion( is )% &
                      &                     state( js )%momentum,       &
                      &   b2CellData = tmpCv,                           &
                      &   vectorID = VEC_ALIGN_PARALLEL_ID )
                end do
                call write_cell_vector_component( sources_grid, mpg,    &
                      &   vectorComponent = edge_sources%source(3)%     &
                      &                     ggd( time_sind )%ion( is )% &
                      &                     momentum,                   &
                      &   b2CellData = totCv,                           &
                      &   vectorID = VEC_ALIGN_PARALLEL_ID )
                totCv(:) = 0.0_IDS_real
                do js = 1, istion(is)
                  do iCv = 1, mpg%nCv
                    tmpCv(iCv) = ( state%sr%smodt(iCv,0,ispion(is,js)) + &
                      &            state%sr%smodt(iCv,1,ispion(is,js)) * &
                      &            state%pl%ua(iCv,ispion(is,js)) +      &
                      &            state%sr%smodt(iCv,2,ispion(is,js)) * &
                      &                     roxa(iCv,ispion(is,js)) +    &
                      &            state%sr%smodt(iCv,3,ispion(is,js)) * &
                      &                    roxa(iCv,ispion(is,js)) *     &
                      &            state%pl%ua(iCv,ispion(is,js)) ) / geo%cvVol(iCv)
                  end do
                  totCv(:) = totCv(:) + tmpCv(:)
                  call write_cell_vector_component( sources_grid, mpg,  &
                      &   vectorComponent = edge_sources%source(4)%     &
                      &                     ggd( time_sind )%ion( is )% &
                      &                     state( js )%momentum,       &
                      &   b2CellData = tmpCv,                           &
                      &   vectorID = VEC_ALIGN_PARALLEL_ID )
                end do
                call write_cell_vector_component( sources_grid, mpg,    &
                      &   vectorComponent = edge_sources%source(4)%     &
                      &                     ggd( time_sind )%ion( is )% &
                      &                     momentum,                   &
                      &   b2CellData = totCv,                           &
                      &   vectorID = VEC_ALIGN_PARALLEL_ID )
                if (switch%use_eirene.ne.0 .and. balance_netcdf.ne.0) then
                  totCv = 0.0_IDS_real
                  do js = 1, istion(is)
                    tmpCv = 0.0_IDS_real
                    do istrai = 1, size( eirene_mc_mapl_smo_bal, 3)
                      tmpCv(:) = tmpCv(:) + &
                         &  eirene_mc_mapl_smo_bal(:,ispion(is,js),istrai)
                    end do
                    tmpCv(:) = tmpCv(:) / geo%cvVol(:)
                    totCv(:) = totCv(:) + tmpCv(:)
                    call write_cell_vector_component( sources_grid, mpg,&
                        &   vectorComponent = edge_sources%source(5)%   &
                        &            ggd( time_sind )%ion( is )%        &
                        &            state( js )%momentum,              &
                        &   b2CellData = tmpCv,                         &
                        &   vectorID = VEC_ALIGN_PARALLEL_ID )
                  end do
                  call write_cell_vector_component( sources_grid, mpg,  &
                      &   vectorComponent = edge_sources%source(5)%     &
                      &            ggd( time_sind )%ion( is )%momentum, &
                      &   b2CellData = totCv,                           &
                      &   vectorID = VEC_ALIGN_PARALLEL_ID )
                  totCv = 0.0_IDS_real
                  do js = 1, istion(is)
                    tmpCv = 0.0_IDS_real
                    do istrai = 1, size( eirene_mc_mmpl_smo_bal, 3)
                      tmpCv(:) = tmpCv(:) + &
                         &  eirene_mc_mmpl_smo_bal(:,ispion(is,js),istrai)
                    end do
                    tmpCv(:) = tmpCv(:) / geo%cvVol(:)
                    totCv(:) = totCv(:) + tmpCv(:)
                    call write_cell_vector_component( sources_grid, mpg,&
                        &   vectorComponent = edge_sources%source(6)%   &
                        &            ggd( time_sind )%ion( is )%        &
                        &            state( js )%momentum,              &
                        &   b2CellData = tmpCv,                         &
                        &   vectorID = VEC_ALIGN_PARALLEL_ID )
                  end do
                  call write_cell_vector_component( sources_grid, mpg,    &
                        &   vectorComponent = edge_sources%source(6)%     &
                        &            ggd( time_sind )%ion( is )%momentum, &
                        &   b2CellData = totCv,                           &
                        &   vectorID = VEC_ALIGN_PARALLEL_ID )
                else
                  totCv = 0.0_IDS_real
                  do js = 1, istion(is)
                    tmpCv(:) = state%srw%b2stbr_smo(:,ispion(is,js)) / geo%cvVol(:)
                    totCv(:) = totCv(:) + tmpCv(:)
                    call write_cell_vector_component( sources_grid, mpg,&
                        &   vectorComponent = edge_sources%source(5)%   &
                        &            ggd( time_sind )%ion( is )%        &
                        &            state( js )%momentum,              &
                        &   b2CellData = tmpCv,                         &
                        &   vectorID = VEC_ALIGN_PARALLEL_ID )
                  end do
                  call write_cell_vector_component( sources_grid, mpg,    &
                        &   vectorComponent = edge_sources%source(5)%     &
                        &            ggd( time_sind )%ion( is )%momentum, &
                        &   b2CellData = totCv,                           &
                        &   vectorID = VEC_ALIGN_PARALLEL_ID )
                end if
                do js = 1, istion(is)
                  tmpCv(:) = state%srw%rsamo(:,ispion(is,js)) / geo%cvVol(:)
                  call write_cell_vector_component( sources_grid, mpg,  &
                      &   vectorComponent = edge_sources%source(7)%     &
                      &                     ggd( time_sind )%ion( is )% &
                      &                     state( js )%momentum,       &
                      &   b2CellData = tmpCv,                           &
                      &   vectorID = VEC_ALIGN_PARALLEL_ID )
                  tmpCv(:) = state%srw%rramo(:,ispion(is,js)) / geo%cvVol(:)
                  call write_cell_vector_component( sources_grid, mpg,  &
                      &   vectorComponent = edge_sources%source(8)%     &
                      &                     ggd( time_sind )%ion( is )% &
                      &                     state( js )%momentum,       &
                      &   b2CellData = tmpCv,                           &
                      &   vectorID = VEC_ALIGN_PARALLEL_ID )
                  tmpCv(:) = state%srw%rcxmo(:,ispion(is,js)) / geo%cvVol(:)
                  call write_cell_vector_component( sources_grid, mpg,  &
                      &   vectorComponent = edge_sources%source(9)%     &
                      &                     ggd( time_sind )%ion( is )% &
                      &                     state( js )%momentum,       &
                      &   b2CellData = tmpCv,                           &
                      &   vectorID = VEC_ALIGN_PARALLEL_ID )
                end do
#if ( ( IMAS_MINOR_VERSION > 38 || IMAS_MAJOR_VERSION > 3 ) && defined(B25_EIRENE) )
                if (switch%use_eirene.ne.0) then
                  do js = 1, istion(is)
                    tmpCv(:) = state%srw%smo0_eir_tot(:,ispion(is,js)) / &
                      &   geo%cvVol(:)
                    call write_cell_vector_component( sources_grid, mpg, &
                      &   vectorComponent = edge_sources%source(13)%     &
                      &                     ggd( time_sind )%ion( is )%  &
                      &                     state( js )%momentum,        &
                      &   b2CellData = tmpCv,                            &
                      &   vectorID = VEC_ALIGN_PARALLEL_ID )
                  end do
                end if
#endif
              end if
            end do

            !! te: Electron Temperature
            tmpCv(:) = state%pl%te(:)/qe
            call write_IDS_quantity( edge_grid, mpg, geo,               &
                &   val = edge_profiles%ggd( time_sind )%electrons%     &
                &         temperature,                                  &
                &   value = tmpCv )
            tmpCv(:) = state%co%hce0(:)/state%dv%ne(:)
            call write_IDS_quantity( transport_grid, mpg, geo,          &
                &   val = edge_transport%model(1)%ggd( time_sind )%     &
                &         electrons%energy%d,                           &
                &   value = tmpCv )
#if ( IMAS_MINOR_VERSION > 38 || IMAS_MAJOR_VERSION > 3 )
            call write_face_scalar( transport_grid, mpg,                &
                &   val = edge_transport%model(1)%ggd( time_sind )%     &
                &         electrons%energy%v_pol,                       &
                &   value = state%co%chve(:,0) )
            call write_face_scalar( transport_grid, mpg,                &
                &   val = edge_transport%model(1)%ggd( time_sind )%     &
                &         electrons%energy%v_radial,                    &
                &   value = state%co%chve(:,1) )
#else
            call split_face_flux( mpg, geo, state%co%chve, tmpFace )
            call write_face_scalar( transport_grid, mpg,                &
                &   val = edge_transport%model(1)%ggd( time_sind )%     &
                &         electrons%energy%v,                           &
                &   value = tmpFace )
#endif
            call divide_by_area( mpg%nFc, geo, state%dv%fhe, flxFace)
#if ( IMAS_MINOR_VERSION > 38 || IMAS_MAJOR_VERSION > 3 )
            call write_face_scalar( transport_grid, mpg,                &
                &   val = edge_transport%model(1)%ggd( time_sind )%     &
                &         electrons%energy%flux_pol,                    &
                &   value = flxFace(:,0) )
            call write_face_scalar( transport_grid, mpg,                &
                &   val = edge_transport%model(1)%ggd( time_sind )%     &
                &         electrons%energy%flux_radial,                 &
                &   value = flxFace(:,1) )
#else
            call split_face_flux( mpg, geo, flxFace, tmpFace )
            call write_face_scalar( transport_grid, mpg,                &
                &   val = edge_transport%model(1)%ggd( time_sind )%     &
                &         electrons%energy%flux,                        &
                &   value = tmpFace )
#endif
            call write_face_scalar( transport_grid, mpg,                &
                &   val = edge_transport%model(1)%ggd( time_sind )%     &
                &         electrons%energy%flux_limiter,                &
                &   value = state%dv%fllime )
            tmpCv(:) = ( state%sr%she(:,0) +                            &
                &        state%sr%she(:,1) * state%pl%te(:) +           &
                &        state%sr%she(:,2) * state%dv%ne(:) +           &
                &        state%sr%she(:,3) * state%pl%te(:) *           &
                &                            state%dv%ne(:) ) / geo%cvVol(:)
            call write_cell_scalar( sources_grid, mpg,                  &
                &   scalar = edge_sources%source(1)%ggd( time_sind )%   &
                &            electrons%energy,                          &
                &   b2CellData = tmpCv )
            tmpCv(:) = state_ext%she(:) / geo%cvVol(:)
            call write_cell_scalar( sources_grid, mpg,                  &
                &   scalar = edge_sources%source(2)%ggd( time_sind )%   &
                &            electrons%energy,                          &
                &   b2CellData = tmpCv )
            tmpCv(:) = ( state%srw%b2stbc_she(:) + state%srw%b2stbm_she(:) ) / geo%cvVol(:)
            call write_cell_scalar( sources_grid, mpg,                  &
                &   scalar = edge_sources%source(3)%ggd( time_sind )%   &
                &            electrons%energy,                          &
                &   b2CellData = tmpCv )
            tmpCv(:) = ( state%sr%shedt(:,0) +                          &
                &        state%sr%shedt(:,1) * state%pl%te(:) +         &
                &        state%sr%shedt(:,2) * state%dv%ne(:) +         &
                &        state%sr%shedt(:,3) * state%pl%te(:) *         &
                &                              state%dv%ne(:) ) / geo%cvVol(:)
            call write_cell_scalar( sources_grid, mpg,                  &
                &   scalar = edge_sources%source(4)%ggd( time_sind )%   &
                &            electrons%energy,                          &
                &   b2CellData = tmpCv )
            if (switch%use_eirene.ne.0 .and. balance_netcdf.ne.0) then
                tmpCv = 0.0_IDS_real
                do is = 1, size( eirene_mc_eael_she_bal, 2)
                    tmpCv(:) = tmpCv(:) + eirene_mc_eael_she_bal(:,is)
                end do
                tmpCv(:) = tmpCv(:) / geo%cvVol(:)
                call write_cell_scalar( sources_grid, mpg,                  &
                    &   scalar = edge_sources%source(5)%ggd( time_sind )%   &
                    &            electrons%energy,                          &
                    &   b2CellData = tmpCv )
                tmpCv = 0.0_IDS_real
                do is = 1, size( eirene_mc_emel_she_bal, 2)
                    tmpCv(:) = tmpCv(:) + eirene_mc_emel_she_bal(:,is)
                end do
                tmpCv(:) = tmpCv(:) / geo%cvVol(:)
                call write_cell_scalar( sources_grid, mpg,                  &
                    &   scalar = edge_sources%source(6)%ggd( time_sind )%   &
                    &            electrons%energy,                          &
                    &   b2CellData = tmpCv )
            else
                tmpCv(:) = state%srw%b2stbr_she(:) / geo%cvVol(:)
                call write_cell_scalar( sources_grid, mpg,                  &
                    &   scalar = edge_sources%source(5)%ggd( time_sind )%   &
                    &            electrons%energy,                          &
                    &   b2CellData = tmpCv )
            end if
            if (balance_netcdf.ne.0) then
                tmpCv(:) = b2stel_she_ion_bal(:) / geo%cvVol(:)
                call write_cell_scalar( sources_grid, mpg,                  &
                    &   scalar = edge_sources%source(7)%ggd( time_sind )%   &
                    &            electrons%energy,                          &
                    &   b2CellData = tmpCv )
                tmpCv(:) = b2stel_she_rec_bal(:) / geo%cvVol(:)
                call write_cell_scalar( sources_grid, mpg,                  &
                    &   scalar = edge_sources%source(8)%ggd( time_sind )%   &
                    &            electrons%energy,                          &
                    &   b2CellData = tmpCv )
                tmpCv(:) = -b2npht_shei_bal(:) / geo%cvVol(:)
                call write_cell_scalar( sources_grid, mpg,                  &
                    &   scalar = edge_sources%source(10)%ggd( time_sind )%  &
                    &            electrons%energy,                          &
                    &   b2CellData = tmpCv )
            end if
            tmpCv(:) = state%srw%b2sihs_joule(:) / geo%cvVol(:)
            call write_cell_scalar( sources_grid, mpg,                      &
                &   scalar = edge_sources%source(11)%ggd( time_sind )%      &
                &            electrons%energy,                              &
                &   b2CellData = tmpCv )
            tmpCv(:) = 0.0_IDS_real
            do is = 0, ns-1
              tmpCv(:) = tmpCv(:) + state%srw%rqrad(:,is)
            end do
#ifdef B25_EIRENE
            do is = 1, natmi
              tmpCv(:) = tmpCv(:) - eneutrad(:,is,0)
            end do
            do is = 1, nmoli
              tmpCv(:) = tmpCv(:) - emolrad(:,is,0)
            end do
            do is = 1, nioni
              tmpCv(:) = tmpCv(:) - eionrad(:,is,0)
            end do
#endif
            tmpCv(:) = tmpCv(:) / geo%cvVol(:)
            call write_cell_scalar( sources_grid, mpg,                      &
                &   scalar = edge_sources%source(12)%ggd( time_sind )%      &
                &            electrons%energy,                              &
                &   b2CellData = tmpCv )
#if ( ( IMAS_MINOR_VERSION > 38 || IMAS_MAJOR_VERSION > 3 ) && defined(B25_EIRENE) )
            if (switch%use_eirene.ne.0) then
              tmpCv(:) = state%srw%she0_eir_tot(:) / geo%cvVol(:)
              call write_cell_scalar( sources_grid, mpg,                    &
                &   scalar = edge_sources%source(13)%ggd( time_sind )%      &
                &            electrons%energy,                              &
                &   b2CellData = tmpCv )
            end if
#endif

            !! pe: Electron pressure
            call b2xppe( mpg%nCv, state%dv%ne, state%pl%te, pe)
            call write_IDS_quantity( edge_grid, mpg, geo,               &
                &   val = edge_profiles%ggd( time_sind )%electrons%     &
                &         pressure,                                     &
                &   value = pe )

            !! ti: (Common) Ion Temperature
            tmpCv(:) = state%pl%ti(:)/qe
            call write_IDS_quantity( edge_grid, mpg, geo,               &
                &   val = edge_profiles%ggd( time_sind )%t_i_average,   &
                &   value = tmpCv )
            if (switch%tn_style.eq.0) then
              tmpCv(:) = state%co%hci0(:)/state%dv%ni(:,0)
            else if (switch%tn_style.eq.1) then
              tmpCv(:) = state%co%hci0(:)/state%dv%ni(:,1)
            else if (switch%tn_style.eq.2) then
              tmpCv(:) = state%co%hci0(:)/(state%dv%ni(:,0)-state%dv%nn(:))
            end if
            call write_IDS_quantity( transport_grid, mpg, geo,          &
                &   val = edge_transport%model(1)%ggd( time_sind )%     &
                &         total_ion_energy%d,                           &
                &   value = tmpCv )
#if ( IMAS_MINOR_VERSION > 38 || IMAS_MAJOR_VERSION > 3 )
            call write_face_scalar( transport_grid, mpg,                &
                &   val = edge_transport%model(1)%ggd( time_sind )%     &
                &         total_ion_energy%v_pol,                       &
                &   value = state%co%chvi(:,0) )
            call write_face_scalar( transport_grid, mpg,                &
                &   val = edge_transport%model(1)%ggd( time_sind )%     &
                &         total_ion_energy%v_radial,                    &
                &   value = state%co%chvi(:,1) )
#else
            call split_face_flux( mpg, geo, state%co%chvi, tmpFace )
            call write_face_scalar( transport_grid, mpg,                &
                &   val = edge_transport%model(1)%ggd( time_sind )%     &
                &         total_ion_energy%v,                           &
                &   value = tmpFace )
#endif
            !! fhi : Ion heat flux
            call divide_by_area( mpg%nFc, geo, state%dv%fhi, flxFace)
#if ( IMAS_MINOR_VERSION > 38 || IMAS_MAJOR_VERSION > 3 )
            call write_face_scalar( transport_grid, mpg,                &
                &   val = edge_transport%model(1)%ggd( time_sind )%     &
                &         total_ion_energy%flux_pol,                    &
                &   value = flxFace(:,0) )
            call write_face_scalar( transport_grid, mpg,                &
                &   val = edge_transport%model(1)%ggd( time_sind )%     &
                &         total_ion_energy%flux_radial,                 &
                &   value = flxFace(:,1) )
#else
            call split_face_flux( mpg, geo, flxFace, tmpFace )
            call write_face_scalar( transport_grid, mpg,                &
                &   val = edge_transport%model(1)%ggd( time_sind )%     &
                &         total_ion_energy%flux,                        &
                &   value = tmpFace )
#endif
            call write_face_scalar( transport_grid, mpg,                &
                &   val = edge_transport%model(1)%ggd( time_sind )%     &
                &         total_ion_energy%flux_limiter,                &
                &   value = state%dv%fllimi )
            !! Ion energy sources
            tmpCv(:) = ( state%sr%shi(:,0) +                            &
                &        state%sr%shi(:,1) * state%pl%ti(:) +           &
                &        state%sr%shi(:,2) * state%dv%ni(:,0) +         &
                &        state%sr%shi(:,3) * state%dv%ni(:,0) *         &
                &        state%pl%ti(:) ) / geo%cvVol(:)
            call write_cell_scalar( sources_grid, mpg,                  &
                &   scalar = edge_sources%source(1)%ggd( time_sind )%   &
                &            total_ion_energy,                          &
                &   b2CellData = tmpCv )
            tmpCv(:) = state_ext%shi(:) / geo%cvVol(:)
            call write_cell_scalar( sources_grid, mpg,                  &
                &   scalar = edge_sources%source(2)%ggd( time_sind )%   &
                &            total_ion_energy,                          &
                &   b2CellData = tmpCv )
            tmpCv(:) = ( state%srw%b2stbc_shi(:) + state%srw%b2stbm_shi(:) ) / geo%cvVol(:)
            call write_cell_scalar( sources_grid, mpg,                  &
                &   scalar = edge_sources%source(3)%ggd( time_sind )%   &
                &            total_ion_energy,                          &
                &   b2CellData = tmpCv )
            tmpCv(:) = ( state%sr%shidt(:,0) +                          &
                &          state%sr%shidt(:,1) * state%pl%ti(:) +       &
                &          state%sr%shidt(:,2) * state%dv%ni(:,0) +     &
                &          state%sr%shidt(:,3) * state%dv%ni(:,0) *     &
                &          state%pl%ti(:) ) / geo%cvVol(:)
            call write_cell_scalar( sources_grid, mpg,                  &
                &   scalar = edge_sources%source(4)%ggd( time_sind )%   &
                &            total_ion_energy,                          &
                &   b2CellData = tmpCv )
            if (switch%use_eirene.ne.0 .and. balance_netcdf.ne.0) then
              tmpCv = 0.0_IDS_real
              do is = 1, size( eirene_mc_eapl_shi_bal, 2)
                tmpCv(:) = tmpCv(:) + eirene_mc_eapl_shi_bal(:,is)
              end do
              tmpCv(:) = tmpCv(:) / geo%cvVol(:)
              call write_cell_scalar( sources_grid, mpg,                  &
                  &   scalar = edge_sources%source(5)%ggd( time_sind )%   &
                  &            total_ion_energy,                          &
                  &   b2CellData = tmpCv )
              tmpCv = 0.0_IDS_real
              do is = 1, size( eirene_mc_empl_shi_bal, 2)
                tmpCv(:) = tmpCv(:) + eirene_mc_empl_shi_bal(:,is)
              end do
              tmpCv(:) = tmpCv(:) / geo%cvVol(:)
              call write_cell_scalar( sources_grid, mpg,                  &
                  &   scalar = edge_sources%source(6)%ggd( time_sind )%   &
                  &            total_ion_energy,                          &
                  &   b2CellData = tmpCv )
            else
              tmpCv(:) = state%srw%b2stbr_shi(:) / geo%cvVol(:)
              call write_cell_scalar( sources_grid, mpg,                  &
                  &   scalar = edge_sources%source(5)%ggd( time_sind )%   &
                  &            total_ion_energy,                          &
                  &   b2CellData = tmpCv )
            end if
            if (balance_netcdf.ne.0) then
              tmpCv(:) = b2stel_shi_ion_bal(:) / geo%cvVol(:)
              call write_cell_scalar( sources_grid, mpg,                  &
                  &   scalar = edge_sources%source(7)%ggd( time_sind )%   &
                  &            total_ion_energy,                          &
                  &   b2CellData = tmpCv )
              tmpCv(:) = b2stel_shi_rec_bal(:) / geo%cvVol(:)
              call write_cell_scalar( sources_grid, mpg,                  &
                  &   scalar = edge_sources%source(8)%ggd( time_sind )%   &
                  &            total_ion_energy,                          &
                  &   b2CellData = tmpCv )
              tmpCv(:) = b2npht_shei_bal(:) / geo%cvVol(:)
              call write_cell_scalar( sources_grid, mpg,                  &
                  &   scalar = edge_sources%source(10)%ggd( time_sind )%  &
                  &            total_ion_energy,                          &
                  &   b2CellData = tmpCv )
            end if
#if ( ( IMAS_MINOR_VERSION > 38 || IMAS_MAJOR_VERSION > 3 ) && defined(B25_EIRENE) )
            if (switch%use_eirene.ne.0) then
              tmpCv(:) = state%srw%shi0_eir_tot(:) / geo%cvVol(:)
              call write_cell_scalar( sources_grid, mpg,                  &
                  &   scalar = edge_sources%source(13)%ggd( time_sind )%  &
                  &            total_ion_energy,                          &
                  &   b2CellData = tmpCv )
            end if
#endif
            do is = 1, nsion
              if (is.le.nspecies) then
                totCv(:) = 0.0_IDS_real
                do js = 1, istion(is)
            !! Ion energy sources resolved by species
                  tmpCv(:) = state%srw%rsahi(:,ispion(is,js)) / geo%cvVol(:)
                  totCv(:) = totCv(:) + tmpCv(:)
                  call write_cell_scalar( sources_grid, mpg,                  &
                      &   scalar = edge_sources%source(7)%ggd( time_sind )%   &
                      &            ion( is )%state( js )%energy,              &
                      &   b2CellData = tmpCv )
                end do
                call write_cell_scalar( sources_grid, mpg,                    &
                      &   scalar = edge_sources%source(7)%ggd( time_sind )%   &
                      &            ion( is )%energy,                          &
                      &   b2CellData = totCv )
                totCv(:) = 0.0_IDS_real
                do js = 1, istion(is)
                  tmpCv(:) = state%srw%rrahi(:,ispion(is,js)) / geo%cvVol(:)
                  totCv(:) = totCv(:) + tmpCv(:)
                  call write_cell_scalar( sources_grid, mpg,                  &
                      &   scalar = edge_sources%source(8)%ggd( time_sind )%   &
                      &            ion( is )%state( js )%energy,              &
                      &   b2CellData = tmpCv )
                end do
                call write_cell_scalar( sources_grid, mpg,                    &
                    &   scalar = edge_sources%source(8)%ggd( time_sind )%     &
                    &            ion( is )%energy,                            &
                        &   b2CellData = totCv )
                totCv(:) = 0.0_IDS_real
                do js = 1, istion(is)
                  tmpCv(:) = state%srw%rcxhi(:,ispion(is,js)) / geo%cvVol(:)
                  totCv(:) = totCv(:) + tmpCv(:)
                  call write_cell_scalar( sources_grid, mpg,                  &
                      &   scalar = edge_sources%source(9)%ggd( time_sind )%   &
                      &            ion( is )%state( js )%energy,              &
                      &   b2CellData = tmpCv )
                end do
                call write_cell_scalar( sources_grid, mpg,                    &
                    &   scalar = edge_sources%source(9)%ggd( time_sind )%     &
                    &            ion( is )%energy,                            &
                    &   b2CellData = totCv )
              else
#ifdef B25_EIRENE
                do js = 1, istion(is)
            !! Test ion temperature
                  tmpCv(:) = tib2(:,ispion(is,js),1)/qe
                  call write_IDS_quantity( edge_grid, mpg, geo,               &
                      &   val = edge_profiles%ggd( time_sind )%ion( is )%     &
                      &         state( js )%temperature,                      &
                      &   value = tmpCv )
                end do
#endif
              end if
            end do

            do is = 1, nsion
              if (is.le.nspecies) then
                !! pb : Ion pressure
                totCv(:) = 0.0_IDS_real
                do js = 1, istion(is)
                  totCv(:) = totCv(:) + state%dv%pa(:,ispion(is,js))
                  call write_IDS_quantity( edge_grid, mpg, geo,             &
                      &   val = edge_profiles%ggd( time_sind )%ion( is )%   &
                      &         state( js )%pressure,                       &
                      &   value = state%dv%pa(:,ispion(is,js)) )
                end do
                call write_IDS_quantity( edge_grid, mpg, geo,               &
                    &   val = edge_profiles%ggd( time_sind )%ion( is )%     &
                    &         pressure,                                     &
                    &   value = totCv )
                !! Kinetic energy density
                totCv(:) = 0.0_IDS_real
                do js = 1, istion(is)
                  totCv(:) = totCv(:) + state%dv%kinrgy(:,ispion(is,js))
                  call write_IDS_quantity( edge_grid, mpg, geo,             &
                      &   val = edge_profiles%ggd( time_sind )%ion( is )%   &
                      &         state( js )%energy_density_kinetic,         &
                      &   value = state%dv%kinrgy(:,ispion(is,js)) )
                end do
                call write_IDS_quantity( edge_grid, mpg, geo,               &
                    &   val = edge_profiles%ggd( time_sind )%ion( is )%     &
                    &         energy_density_kinetic,                       &
                    &   value = totCv )
                do js = 1, istion(is)
                !! Average charge
                  call write_IDS_quantity( edge_grid, mpg, geo,             &
                      &   val = edge_profiles%ggd( time_sind )%ion( is )%   &
                      &         state( js )%z_average,                      &
                      &   value = state%rt%rza(:,ispion(is,js)) )
                !! Average square charge
                  call write_IDS_quantity( edge_grid, mpg, geo,             &
                      &   val = edge_profiles%ggd( time_sind )%ion( is )%   &
                      &         state( js )%z_square_average,               &
                      &   value = state%rt%rz2(:,ispion(is,js)) )
#if IMAS_MAJOR_VERSION > 3
                !! Ionization potential
                  call write_IDS_quantity( edge_grid, mpg, geo,             &
                      &   val = edge_profiles%ggd( time_sind )%ion( is )%   &
                      &         state( js )%ionization_potential,           &
                      &   value = state%rt%rpt(:,ispion(is,js)) )
#else
                !! Ionisation potential
                  call write_IDS_quantity( edge_grid, mpg, geo,             &
                      &   val = edge_profiles%ggd( time_sind )%ion( is )%   &
                      &         state( js )%ionisation_potential,           &
                      &   value = state%rt%rpt(:,ispion(is,js)) )
#endif
                end do
              else
#ifdef B25_EIRENE
                !! Test ion pressure
                totCv(:) = 0.0_IDS_real
                do js = 1, istion(is)
                  tmpCv(:) = dib2(:,ispion(is,js),1)                        &
                      &     *tib2(:,ispion(is,js),1)
                  totCv(:) = totCv(:) + tmpCv(:)
                  call write_IDS_quantity( edge_grid, mpg, geo,             &
                      &   val = edge_profiles%ggd( time_sind )%ion( is )%   &
                      &         state( js )%pressure,                       &
                      &   value = tmpCv )
                end do
                call write_IDS_quantity( edge_grid, mpg, geo,               &
                    &   val = edge_profiles%ggd( time_sind )%ion( is )%     &
                    &         pressure,                                     &
                    &   value = totCv )
                !! Average charge
                do js = 1, istion(is)
                  tmpCv(:) = nchrgi( ispion(is,js) )
                  call write_IDS_quantity( edge_grid, mpg, geo,             &
                      &   val = edge_profiles%ggd( time_sind )%ion( is )%   &
                      &         state( js )%z_average,                      &
                      &   value = tmpCv )
                !! Average square charge
                  tmpCv(:) = nchrgi( js )**2
                  call write_IDS_quantity( edge_grid, mpg, geo,             &
                      &   val = edge_profiles%ggd( time_sind )%ion( is )%   &
                      &         state( js )%z_square_average,               &
                      &   value = tmpCv )
                end do
#endif
              end if
            end do

            !! ni/ne: Total ion density over electron density
            tmpCv(:) = state%dv%ni(:,1)/state%dv%ne(:)
            call write_IDS_quantity( edge_grid, mpg, geo,                    &
                &   val = edge_profiles%ggd( time_sind )%n_i_total_over_n_e, &
                &   value = tmpCv )

            !! Zeff
            call write_IDS_quantity( edge_grid, mpg, geo,                    &
                &   val = edge_profiles%ggd( time_sind )%zeff,               &
                &   value = zeff )

            !! pz: Thermal plasma pressure (electrons+ions)
            call write_IDS_quantity( edge_grid, mpg, geo,                    &
                &   val = edge_profiles%ggd( time_sind )%pressure_thermal,   &
                &   value = state%dv%pz )

#if ( IMAS_MINOR_VERSION > 32 || IMAS_MAJOR_VERSION > 3 )
            !! fch: Total current
            call divide_by_area( mpg%nFc, geo, state%dv%fch, flxFace )
            call write_face_vector_component( edge_grid, mpg,                &
                &   vectorComponent = edge_profiles%ggd( time_sind )%        &
                &                     j_total,                               &
                &   b2FaceData = flxFace(:,0),                               &
                &   vectorID = VEC_ALIGN_POLOIDAL_ID )
            call write_face_vector_component( edge_grid, mpg,                &
                &   vectorComponent = edge_profiles%ggd( time_sind )%        &
                &                     j_total,                               &
                &   b2FaceData = flxFace(:,1),                               &
                &   vectorID = VEC_ALIGN_RADIAL_ID )

            !! fch_p: Parallel current
            call divide_by_area( mpg%nFc, geo, state%dv%fch_p, flxFace )
            call write_face_scalar( edge_grid, mpg,                          &
                &   val = edge_profiles%ggd( time_sind )%j_parallel,         &
                &   value = flxFace(:,0) )
#endif

            !! fchanml: Anomalous current
            call divide_by_area( mpg%nFc, geo, state%dv%fchanml, flxFace )
            call write_face_vector_component( edge_grid, mpg,                &
                &   vectorComponent = edge_profiles%ggd( time_sind )%        &
                &                     j_anomalous,                           &
                &   b2FaceData = flxFace(:,0),                               &
                &   vectorID = VEC_ALIGN_POLOIDAL_ID )
            call write_face_vector_component( edge_grid, mpg,                &
                &   vectorComponent = edge_profiles%ggd( time_sind )%        &
                &                     j_anomalous,                           &
                &   b2FaceData = flxFace(:,1),                               &
                &   vectorID = VEC_ALIGN_RADIAL_ID )

            !! fchinert: Inertial current
            call divide_by_area( mpg%nFc, geo, state%dv%fchinert, flxFace )
            call write_face_vector_component( edge_grid, mpg,                &
                &   vectorComponent = edge_profiles%ggd( time_sind )%        &
                &                     j_inertial,                            &
                &   b2FaceData = flxFace(:,0),                               &
                &   vectorID = VEC_ALIGN_POLOIDAL_ID )
            call write_face_vector_component( edge_grid, mpg,                &
                &   vectorComponent = edge_profiles%ggd( time_sind )%        &
                &                     j_inertial,                            &
                &   b2FaceData = flxFace(:,1),                               &
                &   vectorID = VEC_ALIGN_RADIAL_ID )

            !! fchin: Ion-neutral friction current
            call divide_by_area( mpg%nFc, geo, state%dv%fchin, flxFace )
            call write_face_vector_component( edge_grid, mpg,                &
                &   vectorComponent = edge_profiles%ggd( time_sind )%        &
                &                     j_ion_neutral_friction,                &
                &   b2FaceData = flxFace(:,0),                               &
                &   vectorID = VEC_ALIGN_POLOIDAL_ID )
            call write_face_vector_component( edge_grid, mpg,                &
                &   vectorComponent = edge_profiles%ggd( time_sind )%        &
                &                     j_ion_neutral_friction,                &
                &   b2FaceData = flxFace(:,1),                               &
                &   vectorID = VEC_ALIGN_RADIAL_ID )

            !! fchvispar: Parallel viscosity current
            call divide_by_area( mpg%nFc, geo, state%dv%fchvispar, flxFace )
            call write_face_vector_component( edge_grid, mpg,                &
                &   vectorComponent = edge_profiles%ggd( time_sind )%        &
                &                     j_parallel_viscosity,                  &
                &   b2FaceData = flxFace(:,0),                               &
                &   vectorID = VEC_ALIGN_POLOIDAL_ID )
            call write_face_vector_component( edge_grid, mpg,                &
                &   vectorComponent = edge_profiles%ggd( time_sind )%        &
                &                     j_parallel_viscosity,                  &
                &   b2FaceData = flxFace(:,1),                               &
                &   vectorID = VEC_ALIGN_RADIAL_ID )

            !! fchvisper: Perpendicular viscosity current
            call divide_by_area( mpg%nFc, geo, state%dv%fchvisper, flxFace )
            call write_face_vector_component( edge_grid, mpg,                &
                &   vectorComponent = edge_profiles%ggd( time_sind )%        &
                &                     j_perpendicular_viscosity,             &
                &   b2FaceData = flxFace(:,0),                               &
                &   vectorID = VEC_ALIGN_POLOIDAL_ID )
            call write_face_vector_component( edge_grid, mpg,                &
                &   vectorComponent = edge_profiles%ggd( time_sind )%        &
                &                     j_perpendicular_viscosity,             &
                &   b2FaceData = flxFace(:,1),                               &
                &   vectorID = VEC_ALIGN_RADIAL_ID )

            !! fchvisq: Heat viscosity current
            call divide_by_area( mpg%nFc, geo, state%dv%fchvisq, flxFace )
            call write_face_vector_component( edge_grid, mpg,                &
                &   vectorComponent = edge_profiles%ggd( time_sind )%        &
                &                     j_heat_viscosity,                      &
                &   b2FaceData = flxFace(:,0),                               &
                &   vectorID = VEC_ALIGN_POLOIDAL_ID )
            call write_face_vector_component( edge_grid, mpg,                &
                &   vectorComponent = edge_profiles%ggd( time_sind )%        &
                &                     j_heat_viscosity,                      &
                &   b2FaceData = flxFace(:,1),                               &
                &   vectorID = VEC_ALIGN_RADIAL_ID )

            !! fchdia: Diamagnetic current
            call divide_by_area( mpg%nFc, geo, state%dv%fchdia, flxFace )
            call write_face_vector_component( edge_grid, mpg,                &
                &   vectorComponent = edge_profiles%ggd( time_sind )%        &
                &                     j_diamagnetic,                         &
                &   b2FaceData = flxFace(:,0),                               &
                &   vectorID = VEC_ALIGN_POLOIDAL_ID )
            call write_face_vector_component( edge_grid, mpg,                &
                &   vectorComponent = edge_profiles%ggd( time_sind )%        &
                &                     j_diamagnetic,                         &
                &   b2FaceData = flxFace(:,1),                               &
                &   vectorID = VEC_ALIGN_RADIAL_ID )

            if (switch%use_eirene.ne.0) then
#ifdef B25_EIRENE
                do is = 1, nspecies
                   tmpCv(:) = 0.0_IDS_real
                   totCv(:) = 0.0_IDS_real
                   do iss = 1, natmi
                      if (latmscl(iss).eq.is) then
                        tmpCv(:) = tmpCv(:) + dab2(:,iss,1)*tab2(:,iss,1)
                        totCv(:) = totCv(:) + dab2(:,iss,1)
                      end if
                   end do
                !! Neutral pressure
                   call write_IDS_quantity( edge_grid, mpg, geo,             &
                       &   val = edge_profiles%ggd( time_sind )%             &
                       &         neutral( is )%pressure,                     &
                       &   value = tmpCv )
                !! Neutral density
                   call write_IDS_quantity( edge_grid, mpg, geo,             &
                       &   val = edge_profiles%ggd( time_sind )%             &
                       &         neutral( is )%density,                      &
                       &   value = totCv )
                !! Neutral radiation
                   tmpCv(:) = 0.0_IDS_real
                   do iss = 1, natmi
                      if (latmscl(iss).eq.is) then
                         tmpCv(:) = tmpCv(:) + eneutrad(:,iss,0)
                      end if
                   end do
                   tmpCv(:) = tmpCv(:) / geo%cvVol(:)
                   call write_cell_scalar( sources_grid, mpg,                &
                       &   scalar = edge_sources%source(12)%ggd( time_sind )%&
                       &            neutral( is )%energy,                    &
                       &   b2CellData = tmpCv )
                end do
                do is = 1, natmi
                   js = latmscl(is)
                   ks = isstat(is)
                   tmpCv(:) = tab2(:,is,1)/qe
                   call write_IDS_quantity( edge_grid, mpg, geo,             &
                       &   val = edge_profiles%ggd( time_sind )%             &
                       &         neutral( js )%state( ks )%temperature,      &
                       &   value = tmpCv )
                   call write_IDS_quantity( edge_grid, mpg, geo,             &
                       &   val = edge_profiles%ggd( time_sind )%             &
                       &         neutral( js )%state( ks )%density,          &
                       &   value = dab2(:,is,1) )
                   tmpCv(:) = dab2(:,is,1)*tab2(:,is,1)
                   call write_IDS_quantity( edge_grid, mpg, geo,             &
                       &   val = edge_profiles%ggd( time_sind )%             &
                       &         neutral( js )%state( ks )%pressure,         &
                       &   value = tmpCv )
                   flxFace = 0.0_IDS_real
                   call write_face_vector_component( edge_grid, mpg,         &
                       &   vectorComponent = edge_profiles%ggd( time_sind )% &
                       &                     neutral( js )%state( ks )%      &
                       &                     velocity_diamagnetic,           &
                       &   b2FaceData = flxFace(:,0),                        &
                       &   vectorID = VEC_ALIGN_POLOIDAL_ID )
                   call write_face_vector_component( edge_grid, mpg,         &
                       &   vectorComponent = edge_profiles%ggd( time_sind )% &
                       &                     neutral( js )%state( ks )%      &
                       &                     velocity_diamagnetic,           &
                       &   b2FaceData = flxFace(:,1),                        &
                       &   vectorID = VEC_ALIGN_RADIAL_ID )
                   call write_face_vector_component( edge_grid, mpg,         &
                       &   vectorComponent = edge_profiles%ggd( time_sind )% &
                       &                     neutral( js )%state( ks )%      &
                       &                     velocity_exb,                   &
                       &   b2FaceData = flxFace(:,0),                        &
                       &   vectorID = VEC_ALIGN_POLOIDAL_ID )
                   call write_face_vector_component( edge_grid, mpg,         &
                       &   vectorComponent = edge_profiles%ggd( time_sind )% &
                       &                     neutral( js )%state( ks )%      &
                       &                     velocity_exb,                   &
                       &   b2FaceData = flxFace(:,1),                        &
                       &   vectorID = VEC_ALIGN_RADIAL_ID )
                   if (balance_netcdf.ne.0) then
                     tmpCv = 0.0_IDS_real
                     do istrai = 1, size(eirene_mc_paat_sna_bal,3)
                       tmpCv(:) = tmpCv(:)                                 &
                          &       + eirene_mc_paat_sna_bal(:,is,istrai)    &
                          &       + eirene_mc_pmat_sna_bal(:,is,istrai)    &
                          &       + eirene_mc_piat_sna_bal(:,is,istrai)
                     end do
                     tmpCv(:) = tmpCv(:) / geo%cvVol(:)
                     call write_cell_scalar( sources_grid, mpg,              &
                       &   scalar = edge_sources%source(1)%ggd( time_sind )% &
                       &            neutral( js )%state( ks )%particles,     &
                       &   b2CellData = tmpCv )
                     tmpCv = 0.0_IDS_real
                     do istrai = 1, size(eirene_mc_paat_sna_bal,3)
                       tmpCv(:) = tmpCv(:)                                   &
                          &       + eirene_mc_paat_sna_bal(:,is,istrai)
                     end do
                     tmpCv(:) = tmpCv(:) / geo%cvVol(:)
                     call write_cell_scalar( sources_grid, mpg,              &
                       &   scalar = edge_sources%source(5)%ggd( time_sind )% &
                       &            neutral( js )%state( ks )%particles,     &
                       &   b2CellData = tmpCv )
                     tmpCv = 0.0_IDS_real
                     do istrai = 1, size(eirene_mc_pmat_sna_bal,3)
                       tmpCv(:) = tmpCv(:)                                   &
                          &       + eirene_mc_pmat_sna_bal(:,is,istrai)      &
                          &       + eirene_mc_piat_sna_bal(:,is,istrai)
                     end do
                     tmpCv(:) = tmpCv(:) / geo%cvVol(:)
                     call write_cell_scalar( sources_grid, mpg,              &
                       &   scalar = edge_sources%source(6)%ggd( time_sind )% &
                       &            neutral( js )%state( ks )%particles,     &
                       &   b2CellData = tmpCv )
                   end if
                   tmpCv(:) = eneutrad(:,is,0)/ geo%cvVol(:)
                   call write_cell_scalar( sources_grid, mpg,                &
                       &   scalar = edge_sources%source(12)%                 &
                       &            ggd( time_sind )%neutral( js )%          &
                       &            state( ks )%energy,                      &
                       &   b2CellData = tmpCv )
                end do

                !! Molecular quantities
                do js = nspecies+1, nneut
                   tmpCv(:) = 0.0_IDS_real
                   totCv(:) = 0.0_IDS_real
                   do is = 1, nmoli
                      if (imneut(is).eq.js) then
                        tmpCv(:) = tmpCv(:) + dmb2(:,is,1)*tmb2(:,is,1)
                        totCv(:) = totCv(:) + dmb2(:,is,1)
                      end if
                   end do
                 !! Molecular pressure
                   call write_IDS_quantity( edge_grid, mpg, geo,             &
                       &   val = edge_profiles%ggd( time_sind )%             &
                       &         neutral( js )%pressure,                     &
                       &   value = tmpCv )
                 !! Molecular density
                   call write_IDS_quantity( edge_grid, mpg, geo,             &
                       &   val = edge_profiles%ggd( time_sind )%             &
                       &         neutral( js )%density,                      &
                       &   value = totCv )
                end do

                js = nspecies
                do is = 1, nmoli
                   ks = isstat(natmi+is)
                   if (ks.eq.1) js = js + 1
                   tmpCv(:) = tmb2(:,is,1)/qe
                   call write_IDS_quantity( edge_grid, mpg, geo,             &
                       &   val = edge_profiles%ggd( time_sind )%             &
                       &         neutral( js )%state( ks )%temperature,      &
                       &   value = tmpCv )
                   call write_IDS_quantity( edge_grid, mpg, geo,             &
                       &   val = edge_profiles%ggd( time_sind )%             &
                       &         neutral( js )%state( ks )%density,          &
                       &   value = dmb2(:,is,1) )
                   tmpCv(:) = dmb2(:,is,1)*tmb2(:,is,1)
                   call write_IDS_quantity( edge_grid, mpg, geo,             &
                       &   val = edge_profiles%ggd( time_sind )%             &
                       &         neutral( js )%state( ks )%pressure,         &
                       &   value = tmpCv )
                   tmpFace = 0.0_IDS_real
                   call write_face_vector_component( edge_grid, mpg,         &
                       &   vectorComponent = edge_profiles%ggd( time_sind )% &
                       &                     neutral( js )%state( ks )%      &
                       &                     velocity_diamagnetic,           &
                       &   b2FaceData = tmpFace,                             &
                       &   vectorID = VEC_ALIGN_POLOIDAL_ID )
                   call write_face_vector_component( edge_grid, mpg,         &
                       &   vectorComponent = edge_profiles%ggd( time_sind )% &
                       &                     neutral( js )%state( ks )%      &
                       &                     velocity_diamagnetic,           &
                       &   b2FaceData = tmpFace,                             &
                       &   vectorID = VEC_ALIGN_RADIAL_ID )
                   call write_face_vector_component( edge_grid, mpg,         &
                       &   vectorComponent = edge_profiles%ggd( time_sind )% &
                       &                     neutral( js )%state( ks )%      &
                       &                     velocity_exb,                   &
                       &   b2FaceData = tmpFace,                             &
                       &   vectorID = VEC_ALIGN_POLOIDAL_ID )
                   call write_face_vector_component( edge_grid, mpg,         &
                       &   vectorComponent = edge_profiles%ggd( time_sind )% &
                       &                     neutral( js )%state( ks )%      &
                       &                     velocity_exb,                   &
                       &   b2FaceData = tmpFace,                             &
                       &   vectorID = VEC_ALIGN_RADIAL_ID )
                   if (balance_netcdf.ne.0) then
                     tmpCv = 0.0_IDS_real
                     do istrai = 1, size(eirene_mc_paml_sna_bal,3)
                       tmpCv(:) = tmpCv(:)                                   &
                          &       + eirene_mc_paml_sna_bal(:,is,istrai)      &
                          &       + eirene_mc_pmml_sna_bal(:,is,istrai)      &
                          &       + eirene_mc_piml_sna_bal(:,is,istrai)
                     end do
                     tmpCv(:) = tmpCv(:) / geo%cvVol(:)
                     call write_cell_scalar( sources_grid, mpg,              &
                       &   scalar = edge_sources%source(1)%ggd( time_sind )% &
                       &            neutral( js )%state( ks )%particles,     &
                       &   b2CellData = tmpCv )
                     tmpCv = 0.0_IDS_real
                     do istrai = 1, size(eirene_mc_paml_sna_bal,3)
                       tmpCv(:) = tmpCv(:)                                   &
                          &       + eirene_mc_paml_sna_bal(:,is,istrai)
                     end do
                     tmpCv(:) = tmpCv(:) / geo%cvVol(:)
                     call write_cell_scalar( sources_grid, mpg,              &
                       &   scalar = edge_sources%source(5)%ggd( time_sind )% &
                       &            neutral( js )%state( ks )%particles,     &
                       &   b2CellData = tmpCv )
                     tmpCv = 0.0_IDS_real
                     do istrai = 1, size(eirene_mc_pmml_sna_bal,3)
                       tmpCv(:) = tmpCv(:)                                   &
                           &      + eirene_mc_pmml_sna_bal(:,is,istrai)      &
                           &      + eirene_mc_piml_sna_bal(:,is,istrai)
                     end do
                     tmpCv(:) = tmpCv(:) / geo%cvVol(:)
                     call write_cell_scalar( sources_grid, mpg,              &
                       &   scalar = edge_sources%source(6)%ggd( time_sind )% &
                       &            neutral( js )%state( ks )%particles,     &
                       &   b2CellData = tmpCv )
                   end if
                end do
#endif
            else
                j = 0
                do is = 1, nspecies
                    js = eb2spcr(is)
                    if (.not.is_neutral(js)) cycle
                    j = j + 1
                !! na : Fluid neutral density
                    call write_IDS_quantity( edge_grid, mpg, geo,             &
                        &   val = edge_profiles%ggd( time_sind )%             &
                        &         neutral( j )%density,                       &
                        &   value = state%pl%na(:,js) )
                    call write_IDS_quantity( edge_grid, mpg, geo,             &
                        &   val = edge_profiles%ggd( time_sind )%             &
                        &         neutral( j )%state(1)%density,              &
                        &   value = state%pl%na(:,js) )
                !! ua: Parallel fluid neutral velocity
                    call write_cell_vector_component( edge_grid, mpg,         &
                        &   vectorComponent = edge_profiles%ggd( time_sind )% &
                        &                     neutral( j )%velocity,          &
                        &   b2CellData = state%pl%ua(:,js),                   &
                        &   vectorID = VEC_ALIGN_PARALLEL_ID )
                    call write_cell_vector_component( edge_grid, mpg,         &
                        &   vectorComponent = edge_profiles%ggd( time_sind )% &
                        &                     neutral( j )%state(1)%velocity, &
                        &   b2CellData = state%pl%ua(:,js),                   &
                        &   vectorID = VEC_ALIGN_PARALLEL_ID )
                !! wadia: Diamagnetic fluid neutral velocity
                    flxFace(:,0) = state%dv%wadia(:,0,js)
                    flxFace(:,1) = state%dv%wadia(:,1,js)
                    call write_face_vector_component( edge_grid, mpg,         &
                        &   vectorComponent = edge_profiles%ggd( time_sind )% &
                        &         neutral( j )%state(1)%velocity_diamagnetic, &
                        &   b2FaceData = flxFace(:,0),                        &
                        &   vectorID = VEC_ALIGN_POLOIDAL_ID )
                    call write_face_vector_component( edge_grid, mpg,         &
                        &   vectorComponent = edge_profiles%ggd( time_sind )% &
                        &         neutral( j )%state(1)%velocity_diamagnetic, &
                        &   b2FaceData = flxFace(:,1),                        &
                        &   vectorID = VEC_ALIGN_RADIAL_ID )
                !! vaecrb: ExB fluid neutral velocity
                    flxFace(:,0) = state%dv%vaecrb(:,0,js)
                    flxFace(:,1) = state%dv%vaecrb(:,1,js)
                    call write_face_vector_component( edge_grid, mpg,         &
                        &   vectorComponent = edge_profiles%ggd( time_sind )% &
                        &         neutral( j )%state(1)%velocity_exb,         &
                        &   b2FaceData = flxFace(:,0),                        &
                        &   vectorID = VEC_ALIGN_POLOIDAL_ID )
                    call write_face_vector_component( edge_grid, mpg,         &
                        &   vectorComponent = edge_profiles%ggd( time_sind )% &
                        &         neutral( j )%state(1)%velocity_exb,         &
                        &   b2FaceData = flxFace(:,1),                        &
                        &   vectorID = VEC_ALIGN_RADIAL_ID )
                !! fna: Fluid neutral particle flux
                    call divide_by_area( mpg%nFc, geo, state%dv%fna(:,:,js),  &
                        &   flxFace)
#if ( IMAS_MINOR_VERSION > 38 || IMAS_MAJOR_VERSION > 3 )
                    call write_face_scalar( transport_grid, mpg,              &
                        &   val = edge_transport%model(1)%ggd( time_sind )%   &
                        &         neutral( j )%particles%flux_pol,            &
                        &   value = flxFace(:,0) )
                    call write_face_scalar( transport_grid, mpg,              &
                        &   val = edge_transport%model(1)%ggd( time_sind )%   &
                        &         neutral( j )%particles%flux_radial,         &
                        &   value = flxFace(:,1) )
                    call write_face_scalar( transport_grid, mpg,              &
                        &   val = edge_transport%model(1)%ggd( time_sind )%   &
                        &         neutral( j )%state(1)%particles%flux_pol,   &
                        &   value = flxFace(:,0) )
                    call write_face_scalar( transport_grid, mpg,              &
                        &   val = edge_transport%model(1)%ggd( time_sind )%   &
                        &         neutral( j )%state(1)%particles%flux_radial,&
                        &   value = flxFace(:,1) )
#else
                    call split_face_flux( mpg, geo, flxFace, tmpFace )
                    call write_face_scalar( transport_grid, mpg,              &
                        &   val = edge_transport%model(1)%ggd( time_sind )%   &
                        &         neutral( j )%particles%flux,                &
                        &   value = tmpFace )
                    call write_face_scalar( transport_grid, mpg,              &
                        &   val = edge_transport%model(1)%ggd( time_sind )%   &
                        &         neutral( j )%state(1)%particles%flux,       &
                        &   value = tmpFace )
#endif
                !! pb : Fluid neutral pressure
                    call write_IDS_quantity( edge_grid, mpg, geo,             &
                        &   val = edge_profiles%ggd( time_sind )%             &
                        &         neutral( j )%pressure,                      &
                        &   value = state%dv%pa(:,js) )
                    call write_IDS_quantity( edge_grid, mpg, geo,             &
                        &   val = edge_profiles%ggd( time_sind )%             &
                        &         neutral( j )%state(1)%pressure,             &
                        &   value = state%dv%pa(:,js) )
                !! cdpa: Fluid neutral diffusivity
#if ( IMAS_MINOR_VERSION > 38 || IMAS_MAJOR_VERSION > 3 )
                    call write_face_scalar( transport_grid, mpg,              &
                        &   val = edge_transport%model(1)%ggd( time_sind )%   &
                        &         neutral( j )%particles%d_pol,               &
                        &   value = state%co%cdpa(:,0,js) )
                    call write_face_scalar( transport_grid, mpg,              &
                        &   val = edge_transport%model(1)%ggd( time_sind )%   &
                        &         neutral( j )%particles%d_pol,               &
                        &   value = state%co%cdpa(:,1,js) )
                    call write_face_scalar( transport_grid, mpg,              &
                        &   val = edge_transport%model(1)%ggd( time_sind )%   &
                        &         neutral( j )%state(1)%particles%d_pol,      &
                        &   value = state%co%cdpa(:,0,js) )
                    call write_face_scalar( transport_grid, mpg,              &
                        &   val = edge_transport%model(1)%ggd( time_sind )%   &
                        &         neutral( j )%state(1)%particles%d_radial,   &
                        &   value = state%co%cdpa(:,1,js) )
#else
                    call split_face_flux( mpg, geo, state%co%cdpa(:,:,js),    &
                        &   tmpFace )
                    call write_face_scalar( transport_grid, mpg,              &
                        &   val = edge_transport%model(1)%ggd( time_sind )%   &
                        &         neutral( j )%particles%d,                   &
                        &   value = tmpFace )
                    call write_face_scalar( transport_grid, mpg,              &
                        &   val = edge_transport%model(1)%ggd( time_sind )%   &
                        &         neutral( j )%state(1)%particles%d,          &
                        &   value = tmpFace )
#endif
                !! Fluid neutral kinetic energy density
                    call write_IDS_quantity( edge_grid, mpg, geo,             &
                        &   val = edge_profiles%ggd( time_sind )%             &
                        &         neutral( j )%energy_density_kinetic,        &
                        &   value = state%dv%kinrgy(:,js) )
                    call write_IDS_quantity( edge_grid, mpg, geo,             &
                        &   val = edge_profiles%ggd( time_sind )%             &
                        &         neutral( j )%state(1)%                      &
                        &         energy_density_kinetic,                     &
                        &   value =state%dv%kinrgy(:,js) )
                !! cvsa: Ion diffusivity
                    call write_face_vector_component( transport_grid, mpg,    &
                        &   vectorComponent = edge_transport%model(1)%        &
                        &                     ggd( time_sind )%neutral( j )%  &
                        &                     momentum%d,                     &
                        &   b2FaceData = state%co%cvsa(:,0,js),               &
                        &   vectorID = VEC_ALIGN_POLOIDAL_ID )
                    call write_face_vector_component( transport_grid, mpg,    &
                        &   vectorComponent = edge_transport%model(1)%        &
                        &                     ggd( time_sind )%neutral( j )%  &
                        &                     state(1)%momentum%d,            &
                        &   b2FaceData = state%co%cvsa(:,0,js),               &
                        &   vectorID = VEC_ALIGN_POLOIDAL_ID )
                    call write_face_vector_component( transport_grid, mpg,    &
                        &   vectorComponent = edge_transport%model(1)%        &
                        &                     ggd( time_sind )%neutral( j )%  &
                        &                     momentum%d,                     &
                        &   b2FaceData = state%co%cvsa(:,1,js),               &
                        &   vectorID = VEC_ALIGN_RADIAL_ID )
                    call write_face_vector_component( transport_grid, mpg,    &
                        &   vectorComponent = edge_transport%model(1)%        &
                        &                     ggd( time_sind )%neutral( j )%  &
                        &                     state(1)%momentum%d,            &
                        &   b2FaceData = state%co%cvsa(:,1,js),               &
                        &   vectorID = VEC_ALIGN_RADIAL_ID )
                !! fllim0fna: Fluid neutral flux limiter
                    call split_face_flux( mpg, geo,                           &
                        &   state%co%fllim0fna(:,:,js), tmpFace )
                    call write_face_scalar( transport_grid, mpg,              &
                        &   val = edge_transport%model(1)%ggd( time_sind )%   &
                        &         neutral( j )%particles%flux_limiter,        &
                        &   value = tmpFace )
                    call write_face_scalar( transport_grid, mpg,              &
                        &   val = edge_transport%model(1)%ggd( time_sind )%   &
                        &         neutral( j )%state(1)%                      &
                        &         particles%flux_limiter,                     &
                        &   value = tmpFace )
                !! fllimvisc: Fluid neutral momentum transport flux limit
                    call write_face_vector_component( transport_grid, mpg,    &
                        &   vectorComponent = edge_transport%model(1)%        &
                        &                     ggd( time_sind )%neutral( j )%  &
                        &                     momentum%flux_limiter,          &
                        &   b2FaceData = state%co%fllimvisc(:,js),            &
                        &   vectorID = VEC_ALIGN_PARALLEL_ID )
                    call write_face_vector_component( transport_grid, mpg,    &
                        &   vectorComponent = edge_transport%model(1)%        &
                        &                     ggd( time_sind )%neutral( j )%  &
                        &                     state(1)%momentum%flux_limiter, &
                        &   b2FaceData = state%co%fllimvisc(:,js),            &
                        &   vectorID = VEC_ALIGN_PARALLEL_ID )
                !! sna: Fluid neutral particle sources
                    tmpCv(:) = ( state%sr%sna(:,0,js) +                       &
                        &        state%sr%sna(:,1,js) * state%pl%na(:,js) ) / geo%cvVol(:)
                    call write_cell_scalar( sources_grid, mpg,                &
                        &   scalar = edge_sources%source(1)%                  &
                        &   ggd( time_sind )%neutral( j )%particles,          &
                        &   b2CellData = tmpCv )
                    call write_cell_scalar( sources_grid, mpg,                &
                        &   scalar = edge_sources%source(1)%                  &
                        &   ggd( time_sind )%neutral( j )%state(1)%particles, &
                        &   b2CellData = tmpCv )
                    tmpCv(:) = state_ext%sna(:,js) / geo%cvVol(:)
                    call write_cell_scalar( sources_grid, mpg,                &
                        &   scalar = edge_sources%source(2)%ggd( time_sind )% &
                        &            neutral( j )%particles,                  &
                        &   b2CellData = tmpCv )
                    call write_cell_scalar( sources_grid, mpg,                &
                        &   scalar = edge_sources%source(2)%ggd( time_sind )% &
                        &            neutral( j )%state(1)%particles,         &
                        &   b2CellData = tmpCv )
                    tmpCv(:) = ( state%srw%b2stbc_sna(:,js) +                 &
                        &        state%srw%b2stbm_sna(:,js) ) / geo%cvVol(:)
                    call write_cell_scalar( sources_grid, mpg,                &
                        &   scalar = edge_sources%source(3)%ggd( time_sind )% &
                        &            neutral( j )%particles,                  &
                        &   b2CellData = tmpCv )
                    call write_cell_scalar( sources_grid, mpg,                &
                        &   scalar = edge_sources%source(3)%ggd( time_sind )% &
                        &            neutral( j )%state(1)%particles,         &
                        &   b2CellData = tmpCv )
                    tmpCv(:) = ( state%sr%snadt(:,0,js) +                     &
                        &        state%sr%snadt(:,1,js) * state%pl%na(:,js) ) / geo%cvVol(:)
                    call write_cell_scalar( sources_grid, mpg,                &
                        &   scalar = edge_sources%source(4)%ggd( time_sind )% &
                        &            neutral( j )%particles,                  &
                        &   b2CellData = tmpCv )
                    call write_cell_scalar( sources_grid, mpg,                &
                        &   scalar = edge_sources%source(4)%ggd( time_sind )% &
                        &            neutral( j )%state(1)%particles,         &
                        &   b2CellData = tmpCv )
                    tmpCv(:) = state%srw%b2stbr_sna(:,js) / geo%cvVol(:)
                    call write_cell_scalar( sources_grid, mpg,                &
                        &   scalar = edge_sources%source(5)%ggd( time_sind )% &
                        &            neutral( j )%particles,                  &
                        &   b2CellData = tmpCv )
                    call write_cell_scalar( sources_grid, mpg,                &
                        &   scalar = edge_sources%source(5)%ggd( time_sind )% &
                        &            neutral( j )%state(1)%particles,         &
                        &   b2CellData = tmpCv )
                    tmpCv(:) = state%srw%rsana(:,js) / geo%cvVol(:)
                    call write_cell_scalar( sources_grid, mpg,                &
                        &   scalar = edge_sources%source(7)%ggd( time_sind )% &
                        &            neutral( j )%particles,                  &
                        &   b2CellData = tmpCv )
                    call write_cell_scalar( sources_grid, mpg,                &
                        &   scalar = edge_sources%source(7)%ggd( time_sind )% &
                        &            neutral( j )%state(1)%particles,         &
                        &   b2CellData = tmpCv )
                    tmpCv(:) = state%srw%rrana(:,js) / geo%cvVol(:)
                    call write_cell_scalar( sources_grid, mpg,                    &
                        &   scalar = edge_sources%source(8)%ggd( time_sind )% &
                        &            neutral( j )%particles,                  &
                        &   b2CellData = tmpCv )
                    call write_cell_scalar( sources_grid, mpg,                    &
                        &   scalar = edge_sources%source(8)%ggd( time_sind )% &
                        &            neutral( j )%state(1)%particles,         &
                        &   b2CellData = tmpCv )
                    tmpCv(:) = state%srw%rcxna(:,js) / geo%cvVol(:)
                    call write_cell_scalar( sources_grid, mpg,                &
                        &   scalar = edge_sources%source(9)%ggd( time_sind )% &
                        &            neutral( j )%particles,                  &
                        &   b2CellData = tmpCv )
                    call write_cell_scalar( sources_grid, mpg,                &
                        &   scalar = edge_sources%source(9)%ggd( time_sind )% &
                        &            neutral( j )%state(1)%particles,         &
                        &   b2CellData = tmpCv )
                !! smo: Neutral parallel momentum sources
                    do iCv = 1, mpg%nCv
                      tmpCv(iCv) = ( state%sr%smo(iCv,0,js) +                       &
                        &            state%sr%smo(iCv,1,js) * state%pl%ua(iCv,js) + &
                        &            state%sr%smo(iCv,2,js) * roxa(iCv,js) +        &
                        &            state%sr%smo(iCv,3,js) * roxa(iCv,js) *        &
                        &            state%pl%ua(iCv,js) ) / geo%cvVol(iCv)
                    end do
                    call write_cell_vector_component( sources_grid, mpg,      &
                        &   vectorComponent = edge_sources%source(1)%         &
                        &                     ggd( time_sind )%               &
                        &                     neutral( j )%momentum,          &
                        &   b2CellData = tmpCv,                               &
                        &   vectorID = VEC_ALIGN_PARALLEL_ID )
                    call write_cell_vector_component( sources_grid, mpg,      &
                        &   vectorComponent = edge_sources%source(1)%         &
                        &                     ggd( time_sind )%neutral( j )%  &
                        &                     state(1)%momentum,              &
                        &   b2CellData = tmpCv,                               &
                        &   vectorID = VEC_ALIGN_PARALLEL_ID )
                    tmpCv(:) = state_ext%smo(:,js) / geo%cvVol(:)
                    call write_cell_vector_component( sources_grid, mpg,      &
                        &   vectorComponent = edge_sources%source(2)%         &
                        &                     ggd( time_sind )%neutral( j )%  &
                        &                     momentum,                       &
                        &   b2CellData = tmpCv,                               &
                        &   vectorID = VEC_ALIGN_PARALLEL_ID )
                    call write_cell_vector_component( sources_grid, mpg,      &
                        &   vectorComponent = edge_sources%source(2)%         &
                        &                     ggd( time_sind )%neutral( j )%  &
                        &                     state(1)%momentum,              &
                        &   b2CellData = tmpCv,                               &
                        &   vectorID = VEC_ALIGN_PARALLEL_ID )
                    tmpCv(:) = ( state%srw%b2stbc_smo(:,js) +                 &
                        &        state%srw%b2stbm_smo(:,js) ) / geo%cvVol(:)
                    call write_cell_vector_component( sources_grid, mpg,      &
                        &   vectorComponent = edge_sources%source(3)%         &
                        &                     ggd( time_sind )%neutral( j )%  &
                        &                     momentum,                       &
                        &   b2CellData = tmpCv,                               &
                        &   vectorID = VEC_ALIGN_PARALLEL_ID )
                    call write_cell_vector_component( sources_grid, mpg,      &
                        &   vectorComponent = edge_sources%source(3)%         &
                        &                     ggd( time_sind )%neutral( j )%  &
                        &                     state(1)%momentum,              &
                        &   b2CellData = tmpCv,                               &
                        &   vectorID = VEC_ALIGN_PARALLEL_ID )
                    do iCv = 1, mpg%nCv
                      tmpCv(iCv) = ( state%sr%smodt(iCv,0,js) +                       &
                        &            state%sr%smodt(iCv,1,js) * state%pl%ua(iCv,js) + &
                        &            state%sr%smodt(iCv,2,js) * roxa(iCv,js) +        &
                        &            state%sr%smodt(iCv,3,js) * roxa(iCv,js) *        &
                        &            state%pl%ua(iCv,js) ) / geo%cvVol(iCv)
                    end do
                    call write_cell_vector_component( sources_grid, mpg,      &
                        &   vectorComponent = edge_sources%source(4)%         &
                        &                     ggd( time_sind )%neutral( j )%  &
                        &                     momentum,                       &
                        &   b2CellData = tmpCv,                               &
                        &   vectorID = VEC_ALIGN_PARALLEL_ID )
                    call write_cell_vector_component( sources_grid, mpg,      &
                        &   vectorComponent = edge_sources%source(4)%         &
                        &                     ggd( time_sind )%neutral( j )%  &
                        &                     state(1)%momentum,              &
                        &   b2CellData = tmpCv,                               &
                        &   vectorID = VEC_ALIGN_PARALLEL_ID )
                    tmpCv(:) = state%srw%b2stbr_smo(:,js) / geo%cvVol(:)
                    call write_cell_vector_component( sources_grid, mpg,      &
                        &   vectorComponent = edge_sources%source(5)%         &
                        &            ggd( time_sind )%neutral( j )%momentum,  &
                        &   b2CellData = tmpCv,                               &
                        &   vectorID = VEC_ALIGN_PARALLEL_ID )
                    call write_cell_vector_component( sources_grid, mpg,      &
                        &   vectorComponent = edge_sources%source(5)%         &
                        &            ggd( time_sind )%neutral( j )%           &
                        &            state(1)%momentum,                       &
                        &   b2CellData = tmpCv,                               &
                        &   vectorID = VEC_ALIGN_PARALLEL_ID )
                    tmpCv(:) = state%srw%rsamo(:,js) / geo%cvVol(:)
                    call write_cell_vector_component( sources_grid, mpg,      &
                        &   vectorComponent = edge_sources%source(7)%         &
                        &                     ggd( time_sind )%neutral( j )%  &
                        &                     momentum,                       &
                        &   b2CellData = tmpCv,                               &
                        &   vectorID = VEC_ALIGN_PARALLEL_ID )
                    call write_cell_vector_component( sources_grid, mpg,      &
                        &   vectorComponent = edge_sources%source(7)%         &
                        &                     ggd( time_sind )%neutral( j )%  &
                        &                     state(1)%momentum,              &
                        &   b2CellData = tmpCv,                               &
                        &   vectorID = VEC_ALIGN_PARALLEL_ID )
                    tmpCv(:) = state%srw%rramo(:,js) / geo%cvVol(:)
                    call write_cell_vector_component( sources_grid, mpg,      &
                        &   vectorComponent = edge_sources%source(8)%         &
                        &                     ggd( time_sind )%neutral( j )%  &
                        &                     momentum,                       &
                        &   b2CellData = tmpCv,                               &
                        &   vectorID = VEC_ALIGN_PARALLEL_ID )
                    call write_cell_vector_component( sources_grid, mpg,      &
                        &   vectorComponent = edge_sources%source(8)%         &
                        &                     ggd( time_sind )%neutral( j )%  &
                        &                     state(1)%momentum,              &
                        &   b2CellData = tmpCv,                               &
                        &   vectorID = VEC_ALIGN_PARALLEL_ID )
                    tmpCv(:) = state%srw%rcxmo(:,js) / geo%cvVol(:)
                    call write_cell_vector_component( sources_grid, mpg,      &
                        &   vectorComponent = edge_sources%source(9)%         &
                        &                     ggd( time_sind )%neutral( j )%  &
                        &                     momentum,                       &
                        &   b2CellData = tmpCv,                               &
                        &   vectorID = VEC_ALIGN_PARALLEL_ID )
                    call write_cell_vector_component( sources_grid, mpg,      &
                        &   vectorComponent = edge_sources%source(9)%         &
                        &                     ggd( time_sind )%neutral( j )%  &
                        &                     state(1)%momentum,              &
                        &   b2CellData = tmpCv,                               &
                        &   vectorID = VEC_ALIGN_PARALLEL_ID )
                !! tn: Fluid Neutral Temperature
                    if (nint(zn(js)).eq.1) then
                      tmpCv(:) = state%pl%tn(:)/qe
                    else
                      tmpCv(:) = state%pl%ti(:)/qe
                    end if
                    call write_IDS_quantity( edge_grid, mpg, geo,             &
                        &   val = edge_profiles%ggd( time_sind )%             &
                        &         neutral( j )%temperature,                   &
                        &   value = tmpCv )
                    call write_IDS_quantity( edge_grid, mpg, geo,             &
                        &   val = edge_profiles%ggd( time_sind )%             &
                        &         neutral( j )%state(1)%temperature,          &
                        &   value = tmpCv )
                    if (nint(zn(js)).eq.1) then
                      tmpCv(:) = state%co%hcn0(:)/state%dv%nn(:)
                      call write_IDS_quantity( transport_grid, mpg, geo,      &
                        &   val = edge_transport%model(1)%ggd( time_sind )%   &
                        &         neutral( j )%energy%d,                      &
                        &   value = tmpCv )
                      call write_IDS_quantity( transport_grid, mpg, geo,      &
                        &   val = edge_transport%model(1)%ggd( time_sind )%   &
                        &         neutral( j )%state(1)%energy%d,             &
                        &   value = tmpCv )
                !! fhn : Fluid neutral heat flux
                      call divide_by_area( mpg%nFc, geo, state%dv%fhn, flxFace)
#if ( IMAS_MINOR_VERSION > 38 || IMAS_MAJOR_VERSION > 3 )
                      call write_face_scalar( transport_grid, mpg,            &
                        &   val = edge_transport%model(1)%ggd( time_sind )%   &
                        &         neutral( j )%energy%flux_pol,               &
                        &   value = flxFace(:,0) )
                      call write_face_scalar( transport_grid, mpg,            &
                        &   val = edge_transport%model(1)%ggd( time_sind )%   &
                        &         neutral( j )%energy%flux_radial,            &
                        &   value = flxFace(:,1) )
                      call write_face_scalar( transport_grid, mpg,            &
                        &   val = edge_transport%model(1)%ggd( time_sind )%   &
                        &         neutral( j )%state(1)%energy%flux_pol,      &
                        &   value = flxFace(:,0) )
                      call write_face_scalar( transport_grid, mpg,            &
                        &   val = edge_transport%model(1)%ggd( time_sind )%   &
                        &         neutral( j )%state(1)%energy%flux_radial,   &
                        &   value = flxFace(:,1) )
#else
                      call split_face_flux( mpg, geo, flxFace, tmpFace )
                      call write_face_scalar( transport_grid, mpg,            &
                        &   val = edge_transport%model(1)%ggd( time_sind )%   &
                        &         neutral( j )%energy%flux,                   &
                        &   value = tmpFace )
                      call write_face_scalar( transport_grid, mpg,            &
                        &   val = edge_transport%model(1)%ggd( time_sind )%   &
                        &         neutral( j )%state(1)%energy%flux,          &
                        &   value = tmpFace )
#endif
                !! Fluid neutral energy sources
                      tmpCv(:) = ( state%sr%shn(:,0) +                        &
                        &          state%sr%shn(:,1) * state%pl%tn(:) +       &
                        &          state%sr%shn(:,2) * state%dv%nn(:) +       &
                        &          state%sr%shn(:,3) * state%dv%nn(:) *       &
                        &          state%pl%tn(:) ) / geo%cvVol(:)
                      call write_cell_scalar( sources_grid, mpg,              &
                        &   scalar = edge_sources%source(1)%ggd( time_sind )% &
                        &            neutral( j )%energy,                     &
                        &   b2CellData = tmpCv )
                      call write_cell_scalar( sources_grid, mpg,              &
                        &   scalar = edge_sources%source(1)%ggd( time_sind )% &
                        &            neutral( j )%state(1)%energy,            &
                        &   b2CellData = tmpCv )
                      tmpCv(:) = state%srw%b2stbc_shn(:) / geo%cvVol(:)
                      call write_cell_scalar( sources_grid, mpg,              &
                        &   scalar = edge_sources%source(3)%ggd( time_sind )% &
                        &            neutral( j )%energy,                     &
                        &   b2CellData = tmpCv )
                      call write_cell_scalar( sources_grid, mpg,              &
                        &   scalar = edge_sources%source(3)%ggd( time_sind )% &
                        &            neutral( j )%state(1)%energy,            &
                        &   b2CellData = tmpCv )
                      tmpCv(:) = ( state%sr%shndt(:,0) +                      &
                        &        state%sr%shndt(:,1) * state%pl%tn(:) +       &
                        &        state%sr%shndt(:,2) * state%dv%nn(:) +       &
                        &        state%sr%shndt(:,3) * state%dv%nn(:) *       &
                        &       state%pl%tn(:) ) / geo%cvVol(:)
                      call write_cell_scalar( sources_grid, mpg,              &
                        &   scalar = edge_sources%source(4)%ggd( time_sind )% &
                        &            neutral( j )%energy,                     &
                        &   b2CellData = tmpCv )
                      call write_cell_scalar( sources_grid, mpg,              &
                        &   scalar = edge_sources%source(4)%ggd( time_sind )% &
                        &            neutral( j )%state(1)%energy,            &
                        &   b2CellData = tmpCv )
                      tmpCv(:) = state%srw%b2stbr_shn(:) / geo%cvVol(:)
                      call write_cell_scalar( sources_grid, mpg,              &
                        &   scalar = edge_sources%source(5)%ggd( time_sind )% &
                        &            neutral( j )%energy,                     &
                        &   b2CellData = tmpCv )
                      call write_cell_scalar( sources_grid, mpg,              &
                        &   scalar = edge_sources%source(5)%ggd( time_sind )% &
                        &            neutral( j )%state(1)%energy,            &
                        &   b2CellData = tmpCv )
                    end if
                    tmpCv(:) = state%srw%rsahi(:,js) / geo%cvVol(:)
                    call write_cell_scalar( sources_grid, mpg,              &
                        &   scalar = edge_sources%source(7)%ggd( time_sind )% &
                        &            neutral( j )%energy,                     &
                        &   b2CellData = tmpCv )
                    call write_cell_scalar( sources_grid, mpg,              &
                        &   scalar = edge_sources%source(7)%ggd( time_sind )% &
                        &            neutral( j )%state(1)%energy,            &
                        &   b2CellData = tmpCv )
                    tmpCv(:) = state%srw%rrahi(:,js) / geo%cvVol(:)
                    call write_cell_scalar( sources_grid, mpg,                &
                        &   scalar = edge_sources%source(8)%ggd( time_sind )% &
                        &            neutral( j )%energy,                     &
                        &   b2CellData = tmpCv )
                    call write_cell_scalar( sources_grid, mpg,                &
                        &   scalar = edge_sources%source(8)%ggd( time_sind )% &
                        &            neutral( j )%state(1)%energy,            &
                        &   b2CellData = tmpCv )
                    tmpCv(:) = state%srw%rcxhi(:,js) / geo%cvVol(:)
                    call write_cell_scalar( sources_grid, mpg,                &
                        &   scalar = edge_sources%source(9)%ggd( time_sind )% &
                        &            neutral( j )%energy,                     &
                        &   b2CellData = tmpCv )
                    call write_cell_scalar( sources_grid, mpg,                &
                        &   scalar = edge_sources%source(9)%ggd( time_sind )% &
                        &            neutral( j )%state(1)%energy,            &
                        &   b2CellData = tmpCv )
                end do
            end if

            !! po: Electric potential
            call write_IDS_quantity( edge_grid, mpg, geo,               &
                &   val = edge_profiles%ggd( time_sind )%phi_potential, &
                &   value = state%pl%po )
            !! Current sources
            tmpCv(:) = ( state%sr%sch(:,0) +                            &
                &        state%sr%sch(:,1) * state%pl%po(:) +           &
                &        state%sr%sch(:,2) * state%dv%ne(:) +           &
                &        state%sr%sch(:,3) * state%dv%ne(:) *           &
                &        state%pl%po(:) ) / geo%cvVol(:)
            call write_cell_scalar( sources_grid, mpg,                  &
                &   scalar = edge_sources%source(1)%ggd( time_sind )%   &
                &            current,                                   &
                &   b2CellData = tmpCv )
            tmpCv(:) = state_ext%sch(:) / geo%cvVol(:)
            call write_cell_scalar( sources_grid, mpg,                  &
                &   scalar = edge_sources%source(2)%ggd( time_sind )%   &
                &            current,                                   &
                &   b2CellData = tmpCv )
            tmpCv(:) = ( state%srw%b2stbc_sch(:) + state%srw%b2stbm_sch(:) ) / geo%cvVol(:)
            call write_cell_scalar( sources_grid, mpg,                  &
                &   scalar = edge_sources%source(3)%ggd( time_sind )%   &
                &            current,                                   &
                &   b2CellData = tmpCv )
            tmpCv(:) = ( state%sr%schdt(:,0) +                          &
                &        state%sr%schdt(:,1) * state%pl%po(:) +         &
                &        state%sr%schdt(:,2) * state%dv%ne(:) +         &
                &        state%sr%schdt(:,3) * state%dv%ne(:) *         &
                &        state%pl%po(:) ) / geo%cvVol(:)
            call write_cell_scalar( sources_grid, mpg,                  &
                &   scalar = edge_sources%source(4)%ggd( time_sind )%   &
                &            current,                                   &
                &   b2CellData = tmpCv )
            tmpCv(:) = state%srw%b2stbr_sch(:) / geo%cvVol(:)
            call write_cell_scalar( sources_grid, mpg,                  &
                &   scalar = edge_sources%source(5)%ggd( time_sind )%   &
                &            current,                                   &
                &   b2CellData = tmpCv )
#if ( ( IMAS_MINOR_VERSION > 38 || IMAS_MAJOR_VERSION > 3 ) && defined(B25_EIRENE) )
            if (switch%use_eirene.ne.0) then
              tmpCv(:) = state%srw%sch0_eir_tot(:) / geo%cvVol(:)
              call write_cell_scalar( sources_grid, mpg,                &
                &   scalar = edge_sources%source(13)%ggd( time_sind )%  &
                &            current,                                   &
                &   b2CellData = tmpCv )
            end if
#endif
            !! csig : Electric conductivity
            call write_face_vector_component( transport_grid, mpg,      &
                &   vectorComponent = edge_transport%model(1)%          &
                &         ggd( time_sind )%conductivity,                &
                &   b2FaceData = state%co%csig(:,0),                    &
                &   vectorID = VEC_ALIGN_POLOIDAL_ID )
            call write_face_vector_component( transport_grid, mpg,      &
                &   vectorComponent = edge_transport%model(1)%          &
                &         ggd( time_sind )%conductivity,                &
                &   b2FaceData =state%co%csig(:,1),                     &
                &   vectorID = VEC_ALIGN_RADIAL_ID )

#if ( IMAS_MINOR_VERSION > 21 || IMAS_MAJOR_VERSION > 3 )
            !! write the emissivity data
            !! Process 1. Line and recombination radiation from B2.5 ions
            j = 0
            do is = 1, nspecies
              js = eb2spcr(is)
              if (is_neutral(js) .and. switch%use_eirene.eq.0) then
                j = j + 1
                tmpCv(:) = state%srw%rqrad(:,js) / geo%cvVol(:)
                call write_cell_scalar( radiation_grid, mpg,               &
                    &   scalar = radiation%process(1)%                     &
                    &   ggd( time_sind )%neutral( j )%emissivity,          &
                    &   b2CellData = tmpCv )
#if ( IMAS_MINOR_VERSION > 24 || IMAS_MAJOR_VERSION > 3 )
                call write_cell_scalar( radiation_grid, mpg,               &
                    &   scalar = radiation%process(1)%                     &
                    &   ggd( time_sind )%neutral( j )%                     &
                    &   state(1)%emissivity,                               &
                    &   b2CellData = tmpCv )
#endif
              end if
              totCv(:) = 0.0_IDS_real
              do js = 1, istion(is)
                tmpCv(:) = state%srw%rqrad(:,ispion(is,js)) / geo%cvVol(:)
                totCv(:) = totCv(:) + tmpCv(:)
#if ( IMAS_MINOR_VERSION > 24 || IMAS_MAJOR_VERSION > 3 )
                call write_cell_scalar( radiation_grid, mpg,               &
                    &   scalar = radiation%process(1)%                     &
                    &   ggd( time_sind )%ion( is )%state( js )%emissivity, &
                    &   b2CellData = tmpCv )
#endif
              end do
              call write_cell_scalar( radiation_grid, mpg,                 &
                    &   scalar = radiation%process(1)%                     &
                    &   ggd( time_sind )%ion( is )%emissivity,             &
                    &   b2CellData = totCv )
            end do
            !! Process 2. Bremsstrahlung from B2.5 ions
            totCv(:) = 0.0_IDS_real
            do is = 1, nspecies
              do js = 1, istion(is)
                tmpCv(:) = state%srw%rqbrm(:,ispion(is,js)) / geo%cvVol(:)
                totCv(:) = totCv(:) + tmpCv(:)
#if ( IMAS_MINOR_VERSION > 24 || IMAS_MAJOR_VERSION > 3 )
                call write_cell_scalar( radiation_grid, mpg,               &
                    &   scalar = radiation%process(2)%                     &
                    &   ggd( time_sind )%ion( is )%state( js )%emissivity, &
                    &   b2CellData = tmpCv )
#endif
              end do
              call write_cell_scalar( radiation_grid, mpg,                 &
                  &   scalar = radiation%process(2)%                       &
                  &   ggd( time_sind )%ion( is )%emissivity,               &
                  &   b2CellData = totCv )
            end do
#ifdef B25_EIRENE
            if (switch%use_eirene.ne.0) then
              !! Process 3. Eirene neutrals (atoms and molecules)
              do js = 1, nspecies
                tmpCv(:) = 0.0_IDS_real
                do is = 1, natmi
                  if (latmscl(is).eq.js) then
                    tmpCv(:) = tmpCv(:)-eneutrad(:,is,0)
                  end if
                end do
                tmpCv(:) = tmpCv(:) / geo%cvVol(:)
                call write_cell_scalar( radiation_grid, mpg,               &
                    &   scalar = radiation%process(3)%                     &
                    &   ggd( time_sind )%neutral( js )%emissivity,         &
                    &   b2CellData = tmpCv )
              end do
              do is = 1, natmi
                js = latmscl(is)
                ks = isstat(is)
                tmpCv(:)=-eneutrad(:,is,0) / geo%cvVol(:)
#if ( IMAS_MINOR_VERSION > 24 || IMAS_MAJOR_VERSION > 3 )
                call write_cell_scalar( radiation_grid, mpg,               &
                    &   scalar = radiation%process(3)%ggd( time_sind )%    &
                    &   neutral( js )%state( ks )%emissivity,              &
                    &   b2CellData = tmpCv )
#endif
              end do
              !! Process 4. Eirene molecular ions
            end if
#endif
#endif
#else
            call logmsg( LOGINFO, &
            &   "b2mod_ual_io.B25_process_ids: GGD not available, no plasma state writing" )
#endif
        end if

#if ( IMAS_MINOR_VERSION > 21 || IMAS_MAJOR_VERSION > 3 )
! Summary separatrix data
        allocate( summary%local%separatrix%position%psi( num_time_slices ) )
        summary%local%separatrix%position%psi( time_sind ) = geo%fsPsi(mpg%iFssep)
#if ( IMAS_MINOR_VERSION > 36 || IMAS_MAJOR_VERSION > 3 )
        allocate( summary%local%separatrix_average%position%psi( num_time_slices ) )
        summary%local%separatrix_average%position%psi( time_sind ) = geo%fsPsi(mpg%iFssep)
#endif
        iCv1 = mpg%fcCv(ifsepomp,1)
        iCv2 = mpg%fcCv(ifsepomp,2)
        call write_sourced_value( summary%local%separatrix%t_e, &
           &  0.5_R8 * (state%pl%te(iCv1) + state%pl%te(iCv2))/ev )
        call write_sourced_value( summary%local%separatrix%t_i_average, &
           &  0.5_R8 * (state%pl%ti(iCv1) + state%pl%ti(iCv2))/ev )
        call write_sourced_value( summary%local%separatrix%n_e, &
           &  0.5_R8 * (state%dv%ne(iCv1) + state%dv%ne(iCv2)) )
#if ( IMAS_MINOR_VERSION > 36 || IMAS_MAJOR_VERSION > 3 )
        tmpFace = 1.0_IDS_real
        u = separatrix_average( mpg, geo, state%pl%te, tmpFace )
        call write_sourced_value( summary%local%separatrix_average%t_e, u/ev )
        u = separatrix_average( mpg, geo, state%pl%ti, tmpFace )
        call write_sourced_value( summary%local%separatrix_average%t_i_average, u/ev )
        u = separatrix_average( mpg, geo, state%dv%ne, tmpFace )
        call write_sourced_value( summary%local%separatrix_average%n_e, u )
#endif
        do is = 1, nspecies
          is1 = eb2spcr(is)
          if (nint(zamax(is1)).eq.0) is1 = is1 + 1
          is2 = is1 + nfluids(is) - 1
          nisep = 0.0_R8
          nasum = 0.0_R8
          vtor = 0.0_R8
          do i = is1, is2
            nisep = nisep + &
              &  0.5_R8 * (state%pl%na(iCv1,i) + state%pl%na(iCv2,i))
            nasum = nasum + &
              &  0.5_R8 * (state%pl%na(iCv1,i) + state%pl%na(iCv2,i))
            vtor = vtor + state%pl%ua(omp(icsepomp-1),i)*                         &
              &    sign(geo%cvBb(omp(icsepomp-1),2),geo%cvEb(omp(icsepomp-1),2))/ &
              &    geo%cvBb(omp(icsepomp-1),3)*                                  &
              &    state%pl%na(omp(icsepomp-1),i)
          end do
          if (nasum.gt.0.0_R8) vtor = vtor / nasum
          totCv = 0.0_IDS_real
          tmpCv = 0.0_IDS_real
          do i = is1, is2
            tmpCv(:) = tmpCv(:) + state%pl%na(:,i)
            totCv(:) = totCv(:) + state%pl%ua(:,i)*state%pl%na(:,i)* &
                &                 sign(geo%cvBb(:,2),geo%cvEb(:,2))/ &
                &                 geo%cvBb(:,3)
          end do
          if (nasum.gt.0.0_R8) totCv(:) = totCv(:)/tmpCv(:)
          u = separatrix_average( mpg, geo, tmpCv, tmpFace )
          v = separatrix_average( mpg, geo, totCv, tmpFace )
          select case (is_codes(eb2spcr(is)))
          case ('H')
            call write_sourced_value( summary%local%separatrix%n_i%hydrogen, nisep )
#if ( IMAS_MINOR_VERSION > 41 || IMAS_MAJOR_VERSION > 3 )
            call write_sourced_value( summary%local%separatrix%velocity_phi%hydrogen, -vtor )
#else
            call write_sourced_value( summary%local%separatrix%velocity_tor%hydrogen, -vtor )
#endif
#if ( IMAS_MINOR_VERSION > 36 || IMAS_MAJOR_VERSION > 3 )
            call write_sourced_value( summary%local%separatrix_average%n_i%hydrogen, u )
#if ( IMAS_MINOR_VERSION > 41 || IMAS_MAJOR_VERSION > 3 )
            call write_sourced_value( summary%local%separatrix_average%velocity_phi%hydrogen, -v )
#else
            call write_sourced_value( summary%local%separatrix_average%velocity_tor%hydrogen, -v )
#endif
#endif
          case ('D')
            call write_sourced_value( summary%local%separatrix%n_i%deuterium, nisep )
#if ( IMAS_MINOR_VERSION > 41 || IMAS_MAJOR_VERSION > 3 )
            call write_sourced_value( summary%local%separatrix%velocity_phi%deuterium, -vtor )
#else
            call write_sourced_value( summary%local%separatrix%velocity_tor%deuterium, -vtor )
#endif
#if ( IMAS_MINOR_VERSION > 36 || IMAS_MAJOR_VERSION > 3 )
            call write_sourced_value( summary%local%separatrix_average%n_i%deuterium, u )
#if ( IMAS_MINOR_VERSION > 41 || IMAS_MAJOR_VERSION > 3 )
            call write_sourced_value( summary%local%separatrix_average%velocity_phi%deuterium, -v )
#else
            call write_sourced_value( summary%local%separatrix_average%velocity_tor%deuterium, -v )
#endif
#endif
          case ('T')
            call write_sourced_value( summary%local%separatrix%n_i%tritium, nisep )
#if ( IMAS_MINOR_VERSION > 41 || IMAS_MAJOR_VERSION > 3 )
            call write_sourced_value( summary%local%separatrix%velocity_phi%tritium, -vtor )
#else
            call write_sourced_value( summary%local%separatrix%velocity_tor%tritium, -vtor )
#endif
#if ( IMAS_MINOR_VERSION > 36 || IMAS_MAJOR_VERSION > 3 )
            call write_sourced_value( summary%local%separatrix_average%n_i%tritium, u )
#if ( IMAS_MINOR_VERSION > 41 || IMAS_MAJOR_VERSION > 3 )
            call write_sourced_value( summary%local%separatrix_average%velocity_phi%tritium, -v )
#else
            call write_sourced_value( summary%local%separatrix_average%velocity_tor%tritium, -v )
#endif
#endif
          case ('He')
            if (nint(am(eb2spcr(is))).eq.3) then
              call write_sourced_value( summary%local%separatrix%n_i%helium_3, nisep )
#if ( IMAS_MINOR_VERSION > 41 || IMAS_MAJOR_VERSION > 3 )
              call write_sourced_value( summary%local%separatrix%velocity_phi%helium_3, -vtor )
#else
              call write_sourced_value( summary%local%separatrix%velocity_tor%helium_3, -vtor )
#endif
#if ( IMAS_MINOR_VERSION > 36 || IMAS_MAJOR_VERSION > 3 )
              call write_sourced_value( summary%local%separatrix_average%n_i%helium_3, u )
#if ( IMAS_MINOR_VERSION > 41 || IMAS_MAJOR_VERSION > 3 )
              call write_sourced_value( summary%local%separatrix_average%velocity_phi%helium_3, -v )
#else
              call write_sourced_value( summary%local%separatrix_average%velocity_tor%helium_3, -v )
#endif
#endif
            else if (nint(am(eb2spcr(is))).eq.4) then
              call write_sourced_value( summary%local%separatrix%n_i%helium_4, nisep )
#if ( IMAS_MINOR_VERSION > 41 || IMAS_MAJOR_VERSION > 3 )
              call write_sourced_value( summary%local%separatrix%velocity_phi%helium_4, -vtor )
#else
              call write_sourced_value( summary%local%separatrix%velocity_tor%helium_4, -vtor )
#endif
#if ( IMAS_MINOR_VERSION > 36 || IMAS_MAJOR_VERSION > 3 )
            call write_sourced_value( summary%local%separatrix_average%n_i%helium_4, u )
#if ( IMAS_MINOR_VERSION > 41 || IMAS_MAJOR_VERSION > 3 )
            call write_sourced_value( summary%local%separatrix_average%velocity_phi%helium_4, -v )
#else
            call write_sourced_value( summary%local%separatrix_average%velocity_tor%helium_4, -v )
#endif
#endif
            end if
          case ('Li')
            call write_sourced_value( summary%local%separatrix%n_i%lithium, nisep )
#if ( IMAS_MINOR_VERSION > 41 || IMAS_MAJOR_VERSION > 3 )
            call write_sourced_value( summary%local%separatrix%velocity_phi%lithium, -vtor )
#else
            call write_sourced_value( summary%local%separatrix%velocity_tor%lithium, -vtor )
#endif
#if ( IMAS_MINOR_VERSION > 36 || IMAS_MAJOR_VERSION > 3 )
            call write_sourced_value( summary%local%separatrix_average%n_i%lithium, u )
#if ( IMAS_MINOR_VERSION > 41 || IMAS_MAJOR_VERSION > 3 )
            call write_sourced_value( summary%local%separatrix_average%velocity_phi%lithium, -v )
#else
            call write_sourced_value( summary%local%separatrix_average%velocity_tor%lithium, -v )
#endif
#endif
          case ('Be')
            call write_sourced_value( summary%local%separatrix%n_i%beryllium, nisep )
#if ( IMAS_MINOR_VERSION > 41 || IMAS_MAJOR_VERSION > 3 )
            call write_sourced_value( summary%local%separatrix%velocity_phi%beryllium, -vtor )
#else
            call write_sourced_value( summary%local%separatrix%velocity_tor%beryllium, -vtor )
#endif
#if ( IMAS_MINOR_VERSION > 36 || IMAS_MAJOR_VERSION > 3 )
            call write_sourced_value( summary%local%separatrix_average%n_i%beryllium, u )
#if ( IMAS_MINOR_VERSION > 41 || IMAS_MAJOR_VERSION > 3 )
            call write_sourced_value( summary%local%separatrix_average%velocity_phi%beryllium, -v )
#else
            call write_sourced_value( summary%local%separatrix_average%velocity_tor%beryllium, -v )
#endif
#endif
          case ('C')
            call write_sourced_value( summary%local%separatrix%n_i%carbon, nisep )
#if ( IMAS_MINOR_VERSION > 41 || IMAS_MAJOR_VERSION > 3 )
            call write_sourced_value( summary%local%separatrix%velocity_phi%carbon, -vtor )
#else
            call write_sourced_value( summary%local%separatrix%velocity_tor%carbon, -vtor )
#endif
#if ( IMAS_MINOR_VERSION > 36 || IMAS_MAJOR_VERSION > 3 )
            call write_sourced_value( summary%local%separatrix_average%n_i%carbon, u )
#if ( IMAS_MINOR_VERSION > 41 || IMAS_MAJOR_VERSION > 3 )
            call write_sourced_value( summary%local%separatrix_average%velocity_phi%carbon, -v )
#else
            call write_sourced_value( summary%local%separatrix_average%velocity_tor%carbon, -v )
#endif
#endif
          case ('N')
            call write_sourced_value( summary%local%separatrix%n_i%nitrogen, nisep )
#if ( IMAS_MINOR_VERSION > 41 || IMAS_MAJOR_VERSION > 3 )
            call write_sourced_value( summary%local%separatrix%velocity_phi%nitrogen, -vtor )
#else
            call write_sourced_value( summary%local%separatrix%velocity_tor%nitrogen, -vtor )
#endif
#if ( IMAS_MINOR_VERSION > 36 || IMAS_MAJOR_VERSION > 3 )
            call write_sourced_value( summary%local%separatrix_average%n_i%nitrogen, u )
#if ( IMAS_MINOR_VERSION > 41 || IMAS_MAJOR_VERSION > 3 )
            call write_sourced_value( summary%local%separatrix_average%velocity_phi%nitrogen, -v )
#else
            call write_sourced_value( summary%local%separatrix_average%velocity_tor%nitrogen, -v )
#endif
#endif
          case ('O')
            call write_sourced_value( summary%local%separatrix%n_i%oxygen, nisep )
#if ( IMAS_MINOR_VERSION > 41 || IMAS_MAJOR_VERSION > 3 )
            call write_sourced_value( summary%local%separatrix%velocity_phi%oxygen, -vtor )
#else
            call write_sourced_value( summary%local%separatrix%velocity_tor%oxygen, -vtor )
#endif
#if ( IMAS_MINOR_VERSION > 36 || IMAS_MAJOR_VERSION > 3 )
            call write_sourced_value( summary%local%separatrix_average%n_i%oxygen, u )
#if ( IMAS_MINOR_VERSION > 41 || IMAS_MAJOR_VERSION > 3 )
            call write_sourced_value( summary%local%separatrix_average%velocity_phi%oxygen, -v )
#else
            call write_sourced_value( summary%local%separatrix_average%velocity_tor%oxygen, -v )
#endif
#endif
          case ('Ne')
            call write_sourced_value( summary%local%separatrix%n_i%neon, nisep )
#if ( IMAS_MINOR_VERSION > 41 || IMAS_MAJOR_VERSION > 3 )
            call write_sourced_value( summary%local%separatrix%velocity_phi%neon, -vtor )
#else
            call write_sourced_value( summary%local%separatrix%velocity_tor%neon, -vtor )
#endif
#if ( IMAS_MINOR_VERSION > 36 || IMAS_MAJOR_VERSION > 3 )
            call write_sourced_value( summary%local%separatrix_average%n_i%neon, u )
#if ( IMAS_MINOR_VERSION > 41 || IMAS_MAJOR_VERSION > 3 )
            call write_sourced_value( summary%local%separatrix_average%velocity_phi%neon, -v )
#else
            call write_sourced_value( summary%local%separatrix_average%velocity_tor%neon, -v )
#endif
#endif
          case ('Ar')
            call write_sourced_value( summary%local%separatrix%n_i%argon, nisep )
#if ( IMAS_MINOR_VERSION > 41 || IMAS_MAJOR_VERSION > 3 )
            call write_sourced_value( summary%local%separatrix%velocity_phi%argon, -vtor )
#else
            call write_sourced_value( summary%local%separatrix%velocity_tor%argon, -vtor )
#endif
#if ( IMAS_MINOR_VERSION > 36 || IMAS_MAJOR_VERSION > 3 )
            call write_sourced_value( summary%local%separatrix_average%n_i%argon, u )
#if ( IMAS_MINOR_VERSION > 41 || IMAS_MAJOR_VERSION > 3 )
            call write_sourced_value( summary%local%separatrix_average%velocity_phi%argon, -v )
#else
            call write_sourced_value( summary%local%separatrix_average%velocity_tor%argon, -v )
#endif
#endif
#if ( IMAS_MINOR_VERSION > 30 || IMAS_MAJOR_VERSION > 3 )
          case ('Fe')
            call write_sourced_value( summary%local%separatrix%n_i%iron, nisep )
#if ( IMAS_MINOR_VERSION > 41 || IMAS_MAJOR_VERSION > 3 )
            call write_sourced_value( summary%local%separatrix%velocity_phi%iron, -vtor )
#else
            call write_sourced_value( summary%local%separatrix%velocity_tor%iron, -vtor )
#endif
#if ( IMAS_MINOR_VERSION > 36 || IMAS_MAJOR_VERSION > 3 )
            call write_sourced_value( summary%local%separatrix_average%n_i%iron, u )
#if ( IMAS_MINOR_VERSION > 41 || IMAS_MAJOR_VERSION > 3 )
            call write_sourced_value( summary%local%separatrix_average%velocity_phi%iron, -v )
#else
            call write_sourced_value( summary%local%separatrix_average%velocity_tor%iron, -v )
#endif
#endif
          case ('Kr')
            call write_sourced_value( summary%local%separatrix%n_i%krypton, nisep )
#if ( IMAS_MINOR_VERSION > 41 || IMAS_MAJOR_VERSION > 3 )
            call write_sourced_value( summary%local%separatrix%velocity_phi%krypton, -vtor )
#else
            call write_sourced_value( summary%local%separatrix%velocity_tor%krypton, -vtor )
#endif
#if ( IMAS_MINOR_VERSION > 36 || IMAS_MAJOR_VERSION > 3 )
            call write_sourced_value( summary%local%separatrix_average%n_i%krypton, u )
#if ( IMAS_MINOR_VERSION > 41 || IMAS_MAJOR_VERSION > 3 )
            call write_sourced_value( summary%local%separatrix_average%velocity_phi%krypton, -v )
#else
            call write_sourced_value( summary%local%separatrix_average%velocity_tor%krypton, -v )
#endif
#endif
#endif
          case ('Xe')
            call write_sourced_value( summary%local%separatrix%n_i%xenon, nisep )
#if ( IMAS_MINOR_VERSION > 41 || IMAS_MAJOR_VERSION > 3 )
            call write_sourced_value( summary%local%separatrix%velocity_phi%xenon, -vtor )
#else
            call write_sourced_value( summary%local%separatrix%velocity_tor%xenon, -vtor )
#endif
#if ( IMAS_MINOR_VERSION > 36 || IMAS_MAJOR_VERSION > 3 )
            call write_sourced_value( summary%local%separatrix_average%n_i%xenon, u )
#if ( IMAS_MINOR_VERSION > 41 || IMAS_MAJOR_VERSION > 3 )
            call write_sourced_value( summary%local%separatrix_average%velocity_phi%xenon, -v )
#else
            call write_sourced_value( summary%local%separatrix_average%velocity_tor%xenon, -v )
#endif
#endif
          case ('W')
            call write_sourced_value( summary%local%separatrix%n_i%tungsten, nisep )
#if ( IMAS_MINOR_VERSION > 41 || IMAS_MAJOR_VERSION > 3 )
            call write_sourced_value( summary%local%separatrix%velocity_phi%tungsten, -vtor )
#else
            call write_sourced_value( summary%local%separatrix%velocity_tor%tungsten, -vtor )
#endif
#if ( IMAS_MINOR_VERSION > 36 || IMAS_MAJOR_VERSION > 3 )
            call write_sourced_value( summary%local%separatrix_average%n_i%tungsten, u )
#if ( IMAS_MINOR_VERSION > 41 || IMAS_MAJOR_VERSION > 3 )
            call write_sourced_value( summary%local%separatrix_average%velocity_phi%tungsten, -v )
#else
            call write_sourced_value( summary%local%separatrix_average%velocity_tor%tungsten, -v )
#endif
#endif
          end select
        end do
        call write_sourced_value( summary%local%separatrix%n_i_total, &
          & 0.5_R8 * (state%dv%ni(iCv1,1) + state%dv%ni(iCv2,1)) )
        u = separatrix_average( mpg, geo, state%dv%ni(:,1), tmpFace )
#if ( IMAS_MINOR_VERSION > 36 || IMAS_MAJOR_VERSION > 3 )
        call write_sourced_value( summary%local%separatrix_average%n_i_total, u )
#endif
        call write_sourced_value( summary%local%separatrix%zeff, &
          & 0.5_R8 * (zeff(iCv1) + zeff(iCv2)) )
        u = separatrix_average( mpg, geo, zeff, tmpFace )
#if ( IMAS_MINOR_VERSION > 36 || IMAS_MAJOR_VERSION > 3 )
        call write_sourced_value( summary%local%separatrix_average%zeff, u )
#endif

! Data at limiter tangency point
        if (geometryType.eq.GEOMETRY_LIMITER) then
          u = intvertex_s( mpg%tgVx(1), mpg%nCv, mpg%nVx, mpg, geo%vxVol, state%pl%te )/ev
          call write_sourced_value( summary%local%limiter%t_e, u )
          u = intvertex_s( mpg%tgVx(1), mpg%nCv, mpg%nVx, mpg, geo%vxVol, state%pl%ti )/ev
          call write_sourced_value( summary%local%limiter%t_i_average, u )
          u = intvertex_s( mpg%tgVx(1), mpg%nCv, mpg%nVx, mpg, geo%vxVol, state%dv%ne )
          call write_sourced_value( summary%local%limiter%n_e, u )
          do is = 1, nspecies
            is1 = eb2spcr(is)
            if (nint(zamax(is1)).eq.0) is1 = is1 + 1
            is2 = is1 + nfluids(is) - 1
            nisep = 0.0_R8
            do i = is1, is2
              u = intvertex_s( mpg%tgVx(1), mpg%nCv, mpg%nVx, mpg, geo%vxVol, state%pl%na(1,i) )
              nisep = nisep + u
            end do
            select case (is_codes(eb2spcr(is)))
            case ('H')
              call write_sourced_value( summary%local%limiter%n_i%hydrogen, nisep )
            case ('D')
              call write_sourced_value( summary%local%limiter%n_i%deuterium, nisep )
            case ('T')
              call write_sourced_value( summary%local%limiter%n_i%tritium, nisep )
            case ('He')
              if (nint(am(eb2spcr(is))).eq.3) then
                call write_sourced_value( summary%local%limiter%n_i%helium_3, nisep )
              else if (nint(am(eb2spcr(is))).eq.4) then
                call write_sourced_value( summary%local%limiter%n_i%helium_4, nisep )
              end if
            case ('Li')
              call write_sourced_value( summary%local%limiter%n_i%lithium, nisep )
            case ('Be')
              call write_sourced_value( summary%local%limiter%n_i%beryllium, nisep )
            case ('C')
              call write_sourced_value( summary%local%limiter%n_i%carbon, nisep )
            case ('N')
              call write_sourced_value( summary%local%limiter%n_i%nitrogen, nisep )
            case ('O')
              call write_sourced_value( summary%local%limiter%n_i%oxygen, nisep )
            case ('Ne')
              call write_sourced_value( summary%local%limiter%n_i%neon, nisep )
            case ('Ar')
              call write_sourced_value( summary%local%limiter%n_i%argon, nisep )
#if ( IMAS_MINOR_VERSION > 30 || IMAS_MAJOR_VERSION > 3 )
            case ('Fe')
              call write_sourced_value( summary%local%limiter%n_i%iron, nisep )
            case ('Kr')
              call write_sourced_value( summary%local%limiter%n_i%krypton, nisep )
#endif
            case ('Xe')
              call write_sourced_value( summary%local%limiter%n_i%xenon, nisep )
            case ('W')
              call write_sourced_value( summary%local%limiter%n_i%tungsten, nisep )
            end select
          end do
          u = intvertex_s( mpg%tgVx(1), mpg%nCv, mpg%nVx, mpg, geo%vxVol, state%dv%ni(1,1) )
          call write_sourced_value( summary%local%limiter%n_i_total, u )
          u = intvertex_s( mpg%tgVx(1), mpg%nCv, mpg%nVx, mpg, geo%vxVol, zeff )
          call write_sourced_value( summary%local%limiter%zeff, u )
          call write_sourced_value( summary%local%limiter%flux_expansion, flux_expansion(mpg%strDiv(1)) )
        end if

! Summary divertor plate data
        if (maxval(mpg%strDiv).gt.0) then
#if ( IMAS_MINOR_VERSION > 34 || IMAS_MAJOR_VERSION > 3 )
          allocate ( summary%local%divertor_target( maxval(mpg%strDiv) ) )
#else
          allocate ( summary%local%divertor_plate( maxval(mpg%strDiv) ) )
#endif
          do i = 1, maxval(mpg%strDiv)
#if ( IMAS_MINOR_VERSION > 34 || IMAS_MAJOR_VERSION > 3 )
            call write_sourced_string( summary%local%divertor_target(i)%name, plate_name(i) )
            u = intvertex_s( mpg%ivdiv(i), mpg%nCv, mpg%nVx, mpg, geo%vxVol, state%pl%te )/ev
            call write_sourced_value( summary%local%divertor_target(i)%t_e, u )
            u = intvertex_s( mpg%ivdiv(i), mpg%nCv, mpg%nVx, mpg, geo%vxVol, state%pl%ti )/ev
            call write_sourced_value( summary%local%divertor_target(i)%t_i_average, &
              &  u )
            u = intvertex_s( mpg%ivdiv(i), mpg%nCv, mpg%nVx, mpg, geo%vxVol, state%dv%ne )
            call write_sourced_value( summary%local%divertor_target(i)%n_e, u )
#else
            call write_sourced_string( summary%local%divertor_plate(i)%name, plate_name(i) )
            u = intvertex_s( mpg%ivdiv(i), mpg%nCv, mpg%nVx, mpg, geo%vxVol, state%pl%te )/ev
            call write_sourced_value( summary%local%divertor_plate(i)%t_e, u )
            u = intvertex_s( mpg%ivdiv(i), mpg%nCv, mpg%nVx, mpg, geo%vxVol, state%pl%ti )/ev
            call write_sourced_value( summary%local%divertor_plate(i)%t_i_average, &
              &   )
            u = intvertex_s( mpg%ivdiv(i), mpg%nCv, mpg%nVx, mpg, geo%vxVol, state%dv%ne )
            call write_sourced_value( summary%local%divertor_plate(i)%n_e, u )
#endif
            do is = 1, nspecies
              is1 = eb2spcr(is)
              if (nint(zamax(is1)).eq.0) is1 = is1+1
              is2 = is1 + nfluids(is) - 1
              nisep = 0.0_R8
              do j = is1, is2
                u = intvertex_s( mpg%ivdiv(i), mpg%nCv, mpg%nVx, mpg, geo%vxVol, state%pl%na(1,j) )
                nisep = nisep + u
              end do
              select case (is_codes(eb2spcr(is)))
              case ('H')
#if ( IMAS_MINOR_VERSION > 34 || IMAS_MAJOR_VERSION > 3 )
                call write_sourced_value( summary%local%divertor_target(i)%n_i%hydrogen, nisep )
#else
                call write_sourced_value( summary%local%divertor_plate(i)%n_i%hydrogen, nisep )
#endif
              case ('D')
#if ( IMAS_MINOR_VERSION > 34 || IMAS_MAJOR_VERSION > 3 )
                call write_sourced_value( summary%local%divertor_target(i)%n_i%deuterium, nisep )
#else
                call write_sourced_value( summary%local%divertor_plate(i)%n_i%deuterium, nisep )
#endif
              case ('T')
#if ( IMAS_MINOR_VERSION > 34 || IMAS_MAJOR_VERSION > 3 )
                call write_sourced_value( summary%local%divertor_target(i)%n_i%tritium, nisep )
#else
                call write_sourced_value( summary%local%divertor_plate(i)%n_i%tritium, nisep )
#endif
              case ('He')
                if (nint(am(eb2spcr(is))).eq.3) then
#if ( IMAS_MINOR_VERSION > 34 || IMAS_MAJOR_VERSION > 3 )
                  call write_sourced_value( summary%local%divertor_target(i)%n_i%helium_3, nisep )
#else
                  call write_sourced_value( summary%local%divertor_plate(i)%n_i%helium_3, nisep )
#endif
                else if (nint(am(eb2spcr(is))).eq.4) then
#if ( IMAS_MINOR_VERSION > 34 || IMAS_MAJOR_VERSION > 3 )
                  call write_sourced_value( summary%local%divertor_target(i)%n_i%helium_4, nisep )
#else
                  call write_sourced_value( summary%local%divertor_plate(i)%n_i%helium_4, nisep )
#endif
                end if
              case ('Li')
#if ( IMAS_MINOR_VERSION > 34 || IMAS_MAJOR_VERSION > 3 )
                call write_sourced_value( summary%local%divertor_target(i)%n_i%lithium, nisep )
#else
                call write_sourced_value( summary%local%divertor_plate(i)%n_i%lithium, nisep )
#endif
              case ('Be')
#if ( IMAS_MINOR_VERSION > 34 || IMAS_MAJOR_VERSION > 3 )
                call write_sourced_value( summary%local%divertor_target(i)%n_i%beryllium, nisep )
#else
                call write_sourced_value( summary%local%divertor_plate(i)%n_i%beryllium, nisep )
#endif
              case ('C')
#if ( IMAS_MINOR_VERSION > 34 || IMAS_MAJOR_VERSION > 3 )
                call write_sourced_value( summary%local%divertor_target(i)%n_i%carbon, nisep )
#else
                call write_sourced_value( summary%local%divertor_plate(i)%n_i%carbon, nisep )
#endif
              case ('N')
#if ( IMAS_MINOR_VERSION > 34 || IMAS_MAJOR_VERSION > 3 )
                call write_sourced_value( summary%local%divertor_target(i)%n_i%nitrogen, nisep )
#else
                call write_sourced_value( summary%local%divertor_plate(i)%n_i%nitrogen, nisep )
#endif
              case ('O')
#if ( IMAS_MINOR_VERSION > 34 || IMAS_MAJOR_VERSION > 3 )
                call write_sourced_value( summary%local%divertor_target(i)%n_i%oxygen, nisep )
#else
                call write_sourced_value( summary%local%divertor_plate(i)%n_i%oxygen, nisep )
#endif
              case ('Ne')
#if ( IMAS_MINOR_VERSION > 34 || IMAS_MAJOR_VERSION > 3 )
                call write_sourced_value( summary%local%divertor_target(i)%n_i%neon, nisep )
#else
                call write_sourced_value( summary%local%divertor_plate(i)%n_i%neon, nisep )
#endif
              case ('Ar')
#if ( IMAS_MINOR_VERSION > 34 || IMAS_MAJOR_VERSION > 3 )
                call write_sourced_value( summary%local%divertor_target(i)%n_i%argon, nisep )
#else
                call write_sourced_value( summary%local%divertor_plate(i)%n_i%argon, nisep )
#endif
#if ( IMAS_MINOR_VERSION > 30 || IMAS_MAJOR_VERSION > 3 )
              case ('Fe')
#if ( IMAS_MINOR_VERSION > 34 || IMAS_MAJOR_VERSION > 3 )
                call write_sourced_value( summary%local%divertor_target(i)%n_i%iron, nisep )
#else
                call write_sourced_value( summary%local%divertor_plate(i)%n_i%iron, nisep )
#endif
              case ('Kr')
#if ( IMAS_MINOR_VERSION > 34 || IMAS_MAJOR_VERSION > 3 )
                call write_sourced_value( summary%local%divertor_target(i)%n_i%krypton, nisep )
#else
                call write_sourced_value( summary%local%divertor_plate(i)%n_i%krypton, nisep )
#endif
#endif
              case ('Xe')
#if ( IMAS_MINOR_VERSION > 34 || IMAS_MAJOR_VERSION > 3 )
                call write_sourced_value( summary%local%divertor_target(i)%n_i%xenon, nisep )
#else
                call write_sourced_value( summary%local%divertor_plate(i)%n_i%xenon, nisep )
#endif
              case ('W')
#if ( IMAS_MINOR_VERSION > 34 || IMAS_MAJOR_VERSION > 3 )
                call write_sourced_value( summary%local%divertor_target(i)%n_i%tungsten, nisep )
#else
                call write_sourced_value( summary%local%divertor_plate(i)%n_i%tungsten, nisep )
#endif
              end select
            end do
#if ( IMAS_MINOR_VERSION > 34 || IMAS_MAJOR_VERSION > 3 )
            u = intvertex_s( mpg%ivdiv(i), mpg%nCv, mpg%nVx, mpg, geo%vxVol, state%dv%ni(1,1) )
            call write_sourced_value( summary%local%divertor_target(i)%n_i_total, u )
            u = intvertex_s( mpg%ivdiv(i), mpg%nCv, mpg%nVx, mpg, geo%vxVol, zeff )
            call write_sourced_value( summary%local%divertor_target(i)%zeff, u )
            call write_sourced_value( summary%local%divertor_target(i)%flux_expansion, &
              & flux_expansion(i) )
#else
            u = intvertex_s( mpg%ivdiv(i), mpg%nCv, mpg%nVx, mpg, geo%vxVol, state%dv%ni(1,1) )
            call write_sourced_value( summary%local%divertor_plate(i)%n_i_total, u )
            u = intvertex_s( mpg%ivdiv(i), mpg%nCv, mpg%nVx, mpg, geo%vxVol, zeff )
            call write_sourced_value( summary%local%divertor_plate(i)%zeff, u )
            call write_sourced_value( summary%local%divertor_plate(i)%flux_expansion, &
              & flux_expansion(i) )
#endif
          end do
        end if

        call fill_summary_data( geo, summary )
        u = 0.0_IDS_real
        do iCv = 1, mpg%nCi
          if (geometryType.eq.GEOMETRY_LIMITER .or. &
           &  geometryType.eq.GEOMETRY_SN .or. &
           &  geometryType.eq.GEOMETRY_LFS_SNOWFLAKE_MINUS .or. &
           &  geometryType.eq.GEOMETRY_LFS_SNOWFLAKE_PLUS) then
            if (mpg%cvReg(iCv).ne.2) cycle
          else if (geometryType.eq.GEOMETRY_STELLARATORISLAND) then
            if (mpg%cvReg(iCv).ne.2 .and. mpg%cvReg(iCv).ne.5) cycle
          else if (geometryType.eq.GEOMETRY_CDN .or. &
               &  geometryType.eq.GEOMETRY_DDN_BOTTOM .or. &
               &  geometryType.eq.GEOMETRY_DDN_TOP) then
            if (mpg%cvReg(iCv).ne.2 .and. mpg%cvReg(iCv).ne.6) cycle
          else
            cycle
          end if
          do is = 0, ns-1
            u = u + state%srw%rqrad(iCv,is) + state%srw%rqbrm(iCv,is)
          end do
#ifdef B25_EIRENE
          do is = 1, natmi
            u = u - eneutrad(iCv,is,0)
          end do
#endif
        end do
        if (u.ne.0.0_IDS_real) then
          call write_sourced_value( summary%scrape_off_layer%power_radiated, u )
        end if
#endif

        deallocate(ionstt,istion,ispion)
#ifdef B25_EIRENE
        if (switch%use_eirene.ne.0) then
          deallocate(isstat,imneut,imiion)
          deallocate(in_species)
        end if
#endif
        call logmsg( LOGDEBUG, "b2mod_ual_io.B25_process_ids: done" )
        ncall = ncall + 1

        contains

        function roxa(iCv, is)
        implicit none
        integer, intent(in) :: iCv, is
        real(kind=R8) :: roxa

        roxa = am(is)*mp*state%pl%na(iCv,is)
        return
        end function roxa

#if ( IMAS_MINOR_VERSION > 29 || IMAS_MAJOR_VERSION > 3 )
        subroutine write_timed_integer( ival, ivalue )
            type(ids_signal_int_1d), intent(inout) :: ival
                !< Type of IDS data structure, designed for integer data handling
            integer, intent(in) :: ivalue

            allocate( ival%data( num_slices ) )
            ival%data( slice_index ) = ivalue
            allocate( ival%time( num_slices ) )
            ival%time( slice_index ) = time_slice_value

            return

        end subroutine write_timed_integer
#endif

#if ( IMAS_MINOR_VERSION > 30 || IMAS_MAJOR_VERSION > 3 )
        subroutine write_timed_value( val, value )
            type(ids_signal_flt_1d), intent(inout) :: val
                !< Type of IDS data structure, designed for scalar data handling
            real(IDS_real), intent(in) :: value

            allocate( val%data( num_slices ) )
            val%data( slice_index ) = value
            allocate( val%time( num_slices ) )
            val%time( slice_index ) = time_slice_value

            return

        end subroutine write_timed_value
#endif

    end subroutine B25_process_ids

    !> Process averaged B2.5 data and set it to IMAS IDS.
    !! @note    The \b B25_av_ids routine enables to store data for
    !!          specific batch average. By default it stores single default
    !!          time slice of time slice value 0.0.
    !!          \b num_batch_slices_IN is required to beforehand allocate
    !!          required ggd(:) array of nodes structure and for additional
    !!          checks for correct use of the routine.
    !! @note    Time slice value is set as:
    !!          \b time_slice_value = \b time_step_IN * \b time_slice_ind_IN
    subroutine B25_av_ids( do_description, &
            &   geo, mpg, state, state_ext, state_avg, switch, &
            &   batch_profiles, batch_sources, &
            &   description, equilibrium, &
#if ( IMAS_MINOR_VERSION > 21 || IMAS_MAJOR_VERSION > 3 )
            &   summary, &
#endif
            &   time_IN, shot, run, database, version, &
            &   new_eq_ggd, &
            &   batch_ind_IN, num_batch_slices_IN )
#ifdef NO_OPT
!DIR$ NOOPTIMIZE
#endif
        use b2mod_mwti
        implicit none
        type (geometry), intent(in) :: geo
        type (mapping), intent(inout) :: mpg
        type (B2state), intent(inout) :: state
        type (B2stateExt), intent(inout) :: state_ext
        type (B2Average), intent(in) :: state_avg
        type (switches), intent(inout) :: switch
        logical, intent(in) :: do_description
        type (ids_equilibrium) :: equilibrium !< IDS designed to
            !< store equilibrium data
        type (ids_edge_profiles) :: batch_profiles !< IDS designed to
            !< store data on edge plasma profiles (includes the scrape-off
            !< layer and possibly part of the confined plasma)
        type (ids_edge_sources) :: batch_sources !< IDS designed to store
            !< data on edge plasma sources. Energy terms correspond to the full
            !< kinetic energy equation (i.e. the energy flux takes into account
            !< the energy transported by the particle flux)
        type (ids_dataset_description) :: description !< IDS designed to store
            !< a description of the simulation
#if ( IMAS_MINOR_VERSION > 21 || IMAS_MAJOR_VERSION > 3 )
        type (ids_summary) :: summary !< IDS designed to store
            !< run summary data
#endif
        integer, intent(in) :: shot, run
        character(len=24), intent(in) :: database, version
        real(IDS_real), intent(in), optional :: time_IN !< Time
        integer, intent(in), optional :: batch_ind_IN
            !< Batch index for the current time slice
        integer, intent(in), optional :: num_batch_slices_IN
            !< Total number of batches. It is required to beforehand allocate
            !< required ggd(:) array of nodes structure and for additional
            !< checks for correct use of the routine.
        logical, intent(out) :: new_eq_ggd

        !! Internal variables
        integer :: i, is, js, ks, ion_charge_int, nc
        integer :: batch_index
        real(IDS_real) :: batch_slice_value   !< Time slice value
        character(len=13) :: spclabel         !< Species label
        character(len=5) :: hlp_frm
 !< Type of IDS data structure, designed for handling grid geometry data
#if GGD_MAJOR_VERSION > 0
#if ( IMAS_MINOR_VERSION < 15 && IMAS_MAJOR_VERSION < 4 )
        type(ids_generic_grid_dynamic) :: batch_grid, sources_grid
#else
        type(ids_generic_grid_aos3_root) :: batch_grid, sources_grid
#endif
        real(IDS_real) :: tmpCv( mpg%nCv )
        real(IDS_real) :: totCv( mpg%nCv )
#endif
        !! Procedures
        external species, xertst

        !! ===  SET UP IDS ===
        write(0,*) "Setting data for batch_profiles IDS"
        source = "SOLPS-ITER (batch averaged)"

        !! Preparing database for writing
        call IDS_init( mpg, geo )
        homogeneous_time = 1
        if ( present( time_IN ) ) then
            time = time_IN
        else
            time = 0.0_IDS_real
        end if

        !! Set default time step values
        !! This routine only fills in one time slice at a time
        batch_index = 1
        slice_index = batch_index
        batch_slice_value = time
        num_batch_slices = 1
        num_slices = num_batch_slices
        call xertst( num_batch_slices .ge. batch_index, &
            & "B25_av_ids: Batch index cannot be greater " // &
            & "than total number of batches!" )
        if( present( batch_ind_IN ) ) &
            & call xertst( batch_ind_IN .ge. 1, &
            & "faulty argument batch_ind_IN" )
        if( present( num_batch_slices_IN ) ) &
            & call xertst( num_batch_slices_IN .ge. 1, &
            & "faulty argument num_batch_slices_IN" )

        !! Preparing IDSs for writing
        comment = trim(label)//' (batch averaged)'
        !! 1. Set homogeneous_time to 0 or 1 and other properties
        call write_ids_properties( batch_profiles%ids_properties, &
          &  homogeneous_time )
        call write_ids_properties( batch_sources%ids_properties, &
          &  homogeneous_time )
        if ( do_description ) then
          call write_ids_properties( description%ids_properties, 2 )
#if ( IMAS_MINOR_VERSION > 21 || IMAS_MAJOR_VERSION > 3 )
          call write_ids_properties( summary%ids_properties, &
            &  homogeneous_time )
#endif
        end if

        !! 2. Set code and library data
        code_description = "Batch-averaged IDS from b2mod_ual_io routine"
        call write_ids_code( switch, batch_profiles%code, code_commit, code_description )
        call write_ids_code( switch, batch_sources%code, code_commit, code_description )
#if ( IMAS_MINOR_VERSION > 21 || IMAS_MAJOR_VERSION > 3 )
        if (do_description) &
          &  call write_ids_code( switch, summary%code, code_commit, code_description )
#endif
#if IMAS_MAJOR_VERSION > 3
        if (do_description) &
          &  call write_ids_code_constant( description%code, code_commit, code_description )
#endif

        !! 3. Allocate IDS.time and set it to desired values
        allocate( batch_profiles%time(num_batch_slices) )
        batch_profiles%time(batch_index) = time
        allocate( batch_sources%time(num_batch_slices) )
        batch_sources%time(batch_index) = time
        if (do_description) then
#if IMAS_MAJOR_VERSION < 4
          allocate( description%time(num_batch_slices) )
          description%time(batch_index) = time
#endif
#if ( IMAS_MINOR_VERSION > 21 || IMAS_MAJOR_VERSION > 3 )
          allocate( summary%time(num_batch_slices) )
          summary%time(batch_index) = time
#endif
        end if

        !! Allocate ggd for number of different time steps
        allocate( batch_profiles%ggd( num_batch_slices ) )
#if ( IMAS_MINOR_VERSION > 14 || IMAS_MAJOR_VERSION > 3 )
        allocate( batch_profiles%grid_ggd( num_batch_slices ) )
        allocate( batch_sources%grid_ggd( num_batch_slices ) )
#endif
        allocate (batch_sources%source(1) )
        allocate (batch_sources%source(1)%ggd( num_batch_slices ) )
#ifdef B25_EIRENE
        allocate( batch_sources%source(1)%identifier%name(1) )
        allocate( batch_sources%source(1)%identifier%description(1) )
#if ( IMAS_MINOR_VERSION > 38 || IMAS_MAJOR_VERSION > 3 )
        !! Neutrals
        batch_sources%source(1)%identifier%index = edge_source_identifier%neutrals
        batch_sources%source(1)%identifier%name = &
          &  edge_source_identifier%name( edge_source_identifier%neutrals )
        batch_sources%source(1)%identifier%description = &
          &  edge_source_identifier%description( edge_source_identifier%neutrals )
#else
        !! Total sources due to Eirene species
        batch_sources%source(1)%identifier%index = 0
        batch_sources%source(1)%identifier%name = "Eirene"
        batch_sources%source(1)%identifier%description = &
             &  "Total source due to Eirene species "//trim(source)
#endif
#endif

        call put_equilibrium_data ( mpg, geo, equilibrium, &
#if ( IMAS_MINOR_VERSION > 21 || IMAS_MAJOR_VERSION > 3 )
            &  summary, &
#endif
            &  batch_profiles, database, &
#if ( AL_MAJOR_VERSION > 4 && GGD_MAJOR_VERSION > 0 )
            &  batch_index, &
#endif
            &  time, do_description, new_eq_ggd )
#if ( IMAS_MINOR_VERSION > 21 || IMAS_MAJOR_VERSION > 3 )
        if (do_description) then
#if IMAS_MAJOR_VERSION == 3
          allocate( description%data_entry%user(1) )
          description%data_entry%user = username
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
#elif IMAS_MAJOR_VERSION > 3
          description%type%index = 2
          allocate( description%type%name(1) )
          allocate( description%type%description(1) )
          description%type%name = "simulation"
          description%type%description = "Simulation results from "//trim(source)
          allocate( description%machine(1) )
          description%machine = database
          description%pulse = shot
#endif
          if ( present( time_IN ) ) &
            &  description%simulation%time_current = time_IN
          allocate( description%simulation%workflow(1) )
          description%simulation%workflow = source

          i=index(B25_git_version,'-')
          if (i.gt.0) then
            allocate( summary%tag%name(1) )
            write(hlp_frm,'(a,i2.2,a)') '(a',i-1,')'
            write(summary%tag%name,hlp_frm) B25_git_version(1:i-1)
          end if
        end if
#endif

#if ( IMAS_MINOR_VERSION > 32 || IMAS_MAJOR_VERSION > 3 )
        call write_ids_midplane( batch_profiles%midplane, midplane_id )
        call write_ids_midplane( batch_sources%midplane, midplane_id )
        if (do_description) &
          & call write_ids_midplane( summary%midplane, midplane_id )
#endif

        !! Write grid & grid subsets/subgrids
#if ( IMAS_MINOR_VERSION > 11 || IMAS_MAJOR_VERSION > 3 ) && GGD_MAJOR_VERSION > 0
#if ( IMAS_MINOR_VERSION < 15 && IMAS_MAJOR_VERSION < 4 )
        call b2_IMAS_Fill_Grid_Desc( mpg, geo, &
            &   batch_profiles%ggd( batch_index )%grid )
        call b2_IMAS_Fill_Grid_Desc( mpg, geo, &
            &   batch_sources%ggd( batch_index )%grid )
#else
        call b2_IMAS_Fill_Grid_Desc( mpg, geo, &
            &   batch_profiles%grid_ggd( batch_index ) )
#if AL_MAJOR_VERSION > 4
        allocate( batch_sources%grid_ggd( batch_index )%path(1) )
        batch_sources%grid_ggd( batch_index )%path = &
            &   "#batch_profiles/grid_ggd("//int2str(batch_index)//")"
#else
        call b2_IMAS_Fill_Grid_Desc( mpg, geo, &
            &   batch_sources%grid_ggd( batch_index ) )
#endif
#endif
#endif
#if ( IMAS_MINOR_VERSION > 21 || IMAS_MAJOR_VERSION > 3 )
        if (do_description) &
          & call write_sourced_string( summary%configuration, configuration )
#endif

        !! Allocate and set time slice value
#if ( IMAS_MINOR_VERSION > 14 || IMAS_MAJOR_VERSION > 3 )
        batch_profiles%grid_ggd( batch_index )%time = batch_slice_value
        batch_sources%grid_ggd( batch_index )%time = batch_slice_value
#endif
        batch_profiles%ggd( batch_index )%time = batch_slice_value
        batch_sources%source(1)%ggd( batch_index )%time = batch_slice_value

        !! List of species
        !! Careful here: ion in DD means isonuclear sequence !!
        allocate( batch_profiles%ggd( batch_index )%ion( nspecies ) )
        allocate( batch_profiles%ggd( batch_index )%neutral( nspecies ) )
        allocate( batch_sources%source(1)%ggd( batch_index )%ion( nspecies ) )
        allocate( batch_sources%source(1)%ggd( batch_index )%neutral( nspecies ) )
        ks = 0
        do js = 1, nspecies
          if (is_neutral(ks)) ks = ks + 1  ! Skip the neutral species
          allocate( batch_profiles%ggd( batch_index )%ion( js )%state( nfluids(js) ) )
          allocate( batch_profiles%ggd( batch_index )%ion( js )%element(1) )
          allocate( batch_sources%source(1)%ggd( batch_index )%ion( js )%state( nfluids(js) ) )
          allocate( batch_sources%source(1)%ggd( batch_index )%ion( js )%element(1) )
#if ( IMAS_MAJOR_VERSION < 4 && IMAS_MINOR_VERSION < 42 )
          allocate( batch_profiles%ggd( batch_index )%ion( js )%label(1) )
          allocate( batch_sources%source(1)%ggd( batch_index )%ion( js )%label(1) )
          do is = 1, nfluids(js)
            allocate( batch_profiles%ggd( batch_index )%ion( js )%state( is )%label(1) )
            allocate( batch_sources%source(1)%ggd( batch_index )%ion( js )%state( is )%label(1) )
          end do
#else
          allocate( batch_profiles%ggd( batch_index )%ion( js )%name(1) )
          allocate( batch_sources%source(1)%ggd( batch_index )%ion( js )%name(1) )
          do is = 1, nfluids(js)
            allocate( batch_profiles%ggd( batch_index )%ion( js )%state( is )%name(1) )
            allocate( batch_sources%source(1)%ggd( batch_index )%ion( js )%state( is )%name(1) )
          end do
#endif
          do is = 1, nfluids(js)
            call species( ks, spclabel, .false.)
#if ( IMAS_MAJOR_VERSION < 4 && IMAS_MINOR_VERSION < 42 )
            batch_profiles%ggd( batch_index )%ion( js )%state( is )%label = spclabel
            batch_sources%source(1)%ggd( batch_index )%ion( js )%state( is )%label = spclabel
#else
            batch_profiles%ggd( batch_index )%ion( js )%state( is )%name = spclabel
            batch_sources%source(1)%ggd( batch_index )%ion( js )%state( is )%name = spclabel
#endif
            ! Put minimum Z of the charge state bundle
            ! (z_min = z_max = 0 for a neutral)
            batch_profiles%ggd( batch_index )%ion( js )%state( is )%z_min = zamin( ks )
            batch_sources%source(1)%ggd( batch_index )%ion( js )%state( is )%z_min = &
                  &  zamin( ks )
            ! Put maximum Z of the charge state bundle
            batch_profiles%ggd( batch_index )%ion( js )%state( is )%z_max = zamax( ks )
            batch_sources%source(1)%ggd( batch_index )%ion( js )%state( is )%z_max = &
                  &  zamax( ks )
            if (is.lt.nfluids(js)) ks = ks + 1
          end do

          ! Put ion label identifying the species
#if ( IMAS_MAJOR_VERSION < 4 && IMAS_MINOR_VERSION < 42 )
          batch_profiles%ggd( batch_index )%ion( js )%label = species_list( js )
          batch_sources%source(1)%ggd( batch_index )%ion( js )%label = species_list( js )
#else
          batch_profiles%ggd( batch_index )%ion( js )%name = species_list( js )
          batch_sources%source(1)%ggd( batch_index )%ion( js )%name = species_list( js )
#endif
          ! Put ion charge if single ion in species
          if (nfluids(js).eq.1) then
            ion_charge_int = nint((zamin(ks)+zamax(ks))/2.0_R8)
            batch_profiles%ggd( batch_index )%ion( js )%z_ion = ion_charge_int
            batch_sources%source(1)%ggd( batch_index )%ion( js )%z_ion = ion_charge_int
          end if
          ! Put mass of species
          batch_profiles%ggd( batch_index )%ion( js )%element(1)%a = am( ks )
          batch_sources%source(1)%ggd( batch_index )%ion( js )%element(1)%a = am( ks )
          ! Put nuclear charge
#if IMAS_MAJOR_VERSION < 4
          batch_profiles%ggd( batch_index )%ion( js )%element(1)%z_n = zn( ks )
          batch_sources%source(1)%ggd( batch_index )%ion( js )%element(1)%z_n = zn( ks )
#else
          batch_profiles%ggd( batch_index )%ion( js )%element(1)%z_n = nint(zn(ks))
          batch_sources%source(1)%ggd( batch_index )%ion( js )%element(1)%z_n = nint(zn(ks))
#endif
          ! Put number of atoms
#if ( IMAS_MINOR_VERSION < 15 && IMAS_MAJOR_VERSION < 4 )
          batch_profiles%ggd( batch_index )%ion( js )%element(1)%multiplicity = 1.0_IDS_real
          batch_sources%source(1)%ggd( batch_index )%ion( js )%element(1)%multiplicity = 1.0_IDS_real
#else
          batch_profiles%ggd( batch_index )%ion( js )%element(1)%atoms_n = 1
          batch_sources%source(1)%ggd( batch_index )%ion( js )%element(1)%atoms_n = 1
#endif
          ! Put neutral index
          batch_profiles%ggd( batch_index )%ion( js )%neutral_index = js
          batch_sources%source(1)%ggd( batch_index )%ion( js )%neutral_index = js
          ! Put multiple states flag
          batch_profiles%ggd( batch_index )%ion( js )%multiple_states_flag = 1
          batch_sources%source(1)%ggd( batch_index )%ion( js )%multiple_states_flag = 1

          !! List of neutrals
          allocate( batch_profiles%ggd( batch_index )%neutral( js )%element(1) )
          batch_profiles%ggd( batch_index )%neutral( js )%element(1)%a = am( ks )
#if IMAS_MAJOR_VERSION < 4
          batch_profiles%ggd( batch_index )%neutral( js )%element(1)%z_n = zn( ks )
#else
          batch_profiles%ggd( batch_index )%neutral( js )%element(1)%z_n = nint(zn(ks))
#endif
#if ( IMAS_MINOR_VERSION < 15 && IMAS_MAJOR_VERSION < 4 )
          batch_profiles%ggd( batch_index )%neutral( js )%element(1)%multiplicity = 1.0_IDS_real
#else
          batch_profiles%ggd( batch_index )%neutral( js )%element(1)%atoms_n = 1
#endif
#if ( IMAS_MAJOR_VERSION < 4 && IMAS_MINOR_VERSION < 42 )
          allocate( batch_profiles%ggd( batch_index )%neutral( js )%label(1) )
          allocate( batch_sources%source(1)%ggd( batch_index )%neutral( js )%label(1) )
          batch_profiles%ggd( batch_index )%neutral( js )%label = species_list( js )
          batch_sources%source(1)%ggd( batch_index )%neutral( js )%label = species_list( js )
#else
          allocate( batch_profiles%ggd( batch_index )%neutral( js )%name(1) )
          allocate( batch_sources%source(1)%ggd( batch_index )%neutral( js )%name(1) )
          batch_profiles%ggd( batch_index )%neutral( js )%name = species_list( js )
          batch_sources%source(1)%ggd( batch_index )%neutral( js )%name = species_list( js )
#endif
          batch_profiles%ggd( batch_index )%neutral( js )%ion_index = js
          batch_profiles%ggd( batch_index )%neutral( js )%multiple_states_flag = 0
          allocate( batch_sources%source(1)%ggd( batch_index )%neutral( js )%element(1) )
          batch_sources%source(1)%ggd( batch_index )%neutral( js )%element(1)%a = am( ks )
#if IMAS_MAJOR_VERSION < 4
          batch_sources%source(1)%ggd( batch_index )%neutral( js )%element(1)%z_n = zn( ks )
#else
          batch_sources%source(1)%ggd( batch_index )%neutral( js )%element(1)%z_n = nint(zn(ks))
#endif
#if ( IMAS_MINOR_VERSION < 15 && IMAS_MAJOR_VERSION < 4 )
          batch_sources%source(1)%ggd( batch_index )%neutral( js )%element(1)%multiplicity = 1.0_IDS_real
#else
          batch_sources%source(1)%ggd( batch_index )%neutral( js )%element(1)%atoms_n = 1
#endif
          batch_sources%source(1)%ggd( batch_index )%neutral( js )%ion_index = js
          batch_sources%source(1)%ggd( batch_index )%neutral( js )%multiple_states_flag = 0

          ks = ks + 1
        end do

        !! Write plasma state
        if ( B2_WRITE_DATA ) then
#if ( IMAS_MINOR_VERSION > 11 || IMAS_MAJOR_VERSION > 3 ) && GGD_MAJOR_VERSION > 0
          call logmsg( LOGDEBUG, &
            &   "b2mod_ual_io.B25_av_ids: writing averaged plasma state" )
            !! Find grid subset base indices out of the available grid subset
            !! data stored in the IDS. That is done using IMAS GGD routine
            !! findGridSubsetByName().
#if ( IMAS_MINOR_VERSION < 15 && IMAS_MAJOR_VERSION < 4 )
          iGsCoreBoundary = findGridSubsetByName( batch_profiles%  &
             &   ggd( batch_index )%grid, "Core boundary" )
          iGsInnerMidplane = findGridSubsetByName( batch_profiles% &
             &   ggd( batch_index )%grid, "Inner Midplane" )
          iGsOuterMidplane = findGridSubsetByName( batch_profiles% &
             &   ggd( batch_index )%grid, "Outer Midplane" )
          iGsCore = findGridSubsetByName( batch_profiles%      &
             &   ggd( batch_index )%grid, "Core" )
          iGsSOL = findGridSubsetByName( batch_profiles%       &
             &   ggd( batch_index )%grid, "SOL" )
          if (geo%LSN) then
            iGsIDivertor = findGridSubsetByName( batch_profiles% &
             &   ggd( batch_index )%grid, "Western divertor" )
            iGsODivertor = findGridSubsetByName( batch_profiles% &
             &   ggd( batch_index )%grid, "Eastern divertor" )
          else
            iGsIDivertor = findGridSubsetByName( batch_profiles% &
             &   ggd( batch_index )%grid, "Eastern divertor" )
            iGsODivertor = findGridSubsetByName( batch_profiles% &
             &   ggd( batch_index )%grid, "Western divertor" )
          end if
#else
          iGsCoreBoundary = findGridSubsetByName(   &
             &   batch_profiles%grid_ggd( batch_index ), "Core boundary" )
          iGsInnerMidplane = findGridSubsetByName(  &
             &   batch_profiles%grid_ggd( batch_index ), "Inner Midplane" )
          iGsOuterMidplane = findGridSubsetByName(  &
             &   batch_profiles%grid_ggd( batch_index ), "Outer Midplane" )
          iGsCore = findGridSubsetByName(           &
             &   batch_profiles%grid_ggd( batch_index ), "Core" )
          iGsSOL = findGridSubsetByName(            &
             &   batch_profiles%grid_ggd( batch_index ), "SOL" )
          if (geo%LSN) then
            iGsIDivertor = findGridSubsetByName(      &
             &   batch_profiles%grid_ggd( batch_index ), "Western divertor" )
            iGsODivertor = findGridSubsetByName(      &
             &   batch_profiles%grid_ggd( batch_index ), "Eastern divertor" )
          else
            iGsIDivertor = findGridSubsetByName(      &
             &   batch_profiles%grid_ggd( batch_index ), "Eastern divertor" )
            iGsODivertor = findGridSubsetByName(      &
             &   batch_profiles%grid_ggd( batch_index ), "Western divertor" )
          end if
#endif

#if ( IMAS_MINOR_VERSION < 15 && IMAS_MAJOR_VERSION < 4 )
          batch_grid = batch_profiles%ggd( batch_index )%grid
          sources_grid = batch_sources%ggd( batch_index )%grid
#else
#if AL_MAJOR_VERSION < 5
          batch_grid = batch_profiles%grid_ggd( batch_index )
          sources_grid = batch_sources%grid_ggd( batch_index )
#else
          batch_grid = batch_profiles%grid_ggd( batch_index )
          sources_grid = batch_profiles%grid_ggd( batch_index )
#endif
#endif
          !! na: Ion density
          ks = 0
          do is = 1, nspecies
            if (is_neutral(ks)) ks = ks + 1
            tmpCv(:) = 0.0_IDS_real
            do js = 1, nfluids(is)
              call write_IDS_quantity( batch_grid, mpg, geo,            &
                  &   val = batch_profiles%ggd( batch_index )%          &
                  &         ion( is )%state( js )%density,              &
                  &   value = state_avg%na_mean(:,ks) )
              tmpCv(:) = tmpCv(:) + state_avg%na_mean(:,ks)
              ks = ks + 1
            end do
            call write_IDS_quantity( batch_grid, mpg, geo,              &
                &   val = batch_profiles%ggd( batch_index )%            &
                &         ion(is)%density,                              &
                &   value = tmpCv )
          end do
          !! sna: Ion particle sources
          ks = 0
          do is = 1, nspecies
            if (is_neutral(ks)) ks = ks + 1
            totCv(:) = 0.0_IDS_real
            do js = 1, nfluids(is)
              tmpCv(:) = state_avg%sna_mean(:,ks) / geo%cvVol(:)
              totCv(:) = totCv(:) + tmpCv(:)
              call write_cell_scalar( sources_grid, mpg,                  &
                  &   scalar = batch_sources%source(1)%                   &
                  &   ggd( batch_index )%ion( is )%state( js )%particles, &
                  &   b2CellData = tmpCv )
              ks = ks + 1
            end do
            call write_cell_scalar( sources_grid, mpg,                    &
                &   scalar = batch_sources%source(1)%                     &
                &   ggd( batch_index )%ion( is )%particles,               &
                &   b2CellData = totCv )
          end do
          !! ua: Parallel ion velocity
          ks = 0
          do is = 1, nspecies
            if (is_neutral(ks)) ks = ks + 1
            do js = 1, nfluids(is)
              call write_cell_vector_component( batch_grid, mpg,          &
                  &   vectorComponent = batch_profiles%                   &
                  &   ggd( batch_index )%ion( is )%state( js )%velocity,  &
                  &   b2CellData = state_avg%ua_mean(:,ks),               &
                  &   vectorID = VEC_ALIGN_PARALLEL_ID )
              ks = ks + 1
            end do
          end do
          !! smo: Ion parallel momentum sources
          ks = 0
          do is = 1, nspecies
            if (is_neutral(ks)) ks = ks + 1
            totCv(:) = 0.0_IDS_real
            do js = 1, nfluids(is)
              tmpCv(:) = state_avg%smo_mean(:,ks) / geo%cvVol(:)
              totCv(:) = totCv(:) + tmpCv(:)
              call write_cell_vector_component( sources_grid, mpg,        &
                  &   vectorComponent = batch_sources%source(1)%          &
                  &                     ggd( batch_index )%ion( is )%     &
                  &                     state( js )%momentum,             &
                  &   b2CellData = tmpCv,                                 &
                  &   vectorID = VEC_ALIGN_PARALLEL_ID )
              ks = ks + 1
            end do
            call write_cell_vector_component( sources_grid, mpg,          &
                &   vectorComponent = batch_sources%source(1)%            &
                &                     ggd( batch_index )%ion( is )%       &
                &                     momentum,                           &
                &   b2CellData = totCv,                                   &
                &   vectorID = VEC_ALIGN_PARALLEL_ID )
          end do
          !! te: Electron Temperature
          tmpCv(:) = state_avg%te_mean(:)/qe
          call write_IDS_quantity( batch_grid, mpg, geo,                 &
              &   val = batch_profiles%ggd( batch_index )%electrons%     &
              &         temperature,                                     &
              &   value = tmpCv )
          tmpCv(:) = state_avg%she_mean(:) / geo%cvVol(:)
          call write_cell_scalar( sources_grid, mpg,                     &
              &   scalar = batch_sources%source(1)%ggd( batch_index )%   &
              &            electrons%energy,                             &
              &   b2CellData = tmpCv )
          !! ti: (Common) Ion Temperature
          tmpCv(:) = state_avg%ti_mean(:)/qe
          call write_IDS_quantity( batch_grid, mpg, geo,                 &
              &   val = batch_profiles%ggd( batch_index )%t_i_average,   &
              &   value = tmpCv )
          !! Ion energy sources
          tmpCv(:) = state_avg%shi_mean(:) / geo%cvVol(:)
          call write_cell_scalar( sources_grid, mpg,                     &
              &   scalar = batch_sources%source(1)%ggd( batch_index )%   &
              &            total_ion_energy,                             &
              &   b2CellData = tmpCv )

          js = 0
          do is = 1, nspecies
            !! sna: Neutral particle sources
            tmpCv(:) = state_avg%sna_mean(:,js) / geo%cvVol(:)
            call write_cell_scalar( sources_grid, mpg,                   &
                &   scalar = batch_sources%source(1)%                    &
                &   ggd( batch_index )%neutral( is )%particles,          &
                &   b2CellData = tmpCv )
            !! smo: Neutral parallel momentum sources
            tmpCv(:) = state_avg%smo_mean(:,js) / geo%cvVol(:)
            call write_cell_vector_component( sources_grid, mpg,         &
                &   vectorComponent = batch_sources%source(1)%           &
                &                     ggd( batch_index )%                &
                &                     neutral( is )%momentum,            &
                &   b2CellData = tmpCv,                                  &
                &   vectorID = VEC_ALIGN_PARALLEL_ID )
            js = js + nfluids(is) + 1
          end do

          !! po: Electric potential
          call write_IDS_quantity( batch_grid, mpg, geo,                 &
              &   val = batch_profiles%ggd( batch_index )%phi_potential, &
              &   value = state_avg%po_mean )
#else
          call logmsg( LOGINFO, &
              &   "b2mod_ual_io.B25_av_ids: GGD not available, no averaged plasma state writing" )
#endif
        end if

        nc = max(mpg%nncut,1)
#if ( IMAS_MINOR_VERSION > 21 || IMAS_MAJOR_VERSION > 3 )
        if (do_description) then
! Summary separatrix data
          allocate( summary%local%separatrix%position%psi( num_batch_slices ) )
          summary%local%separatrix%position%psi( batch_index ) = geo%fsPsi(mpg%iFssep)
#if ( IMAS_MINOR_VERSION > 36 || IMAS_MAJOR_VERSION > 3 )
          allocate( summary%local%separatrix_average%position%psi( num_batch_slices ) )
          summary%local%separatrix_average%position%psi( batch_index ) = geo%fsPsi(mpg%iFssep)
#endif
          call write_errored_value( summary%local%separatrix%t_e,         &
              &  tesepm_av(nc), tesepm_std(nc) )
          call write_errored_value( summary%local%separatrix%t_i_average, &
              &  tisepm_av(nc), tisepm_std(nc) )
          call write_errored_value( summary%local%separatrix%n_e,         &
              &  nesepm_av(nc), nesepm_std(nc) )
! Summary divertor plate data
          if (maxval(mpg%strDiv).gt.0) then
#if ( IMAS_MINOR_VERSION > 34 || IMAS_MAJOR_VERSION > 3 )
            allocate ( summary%local%divertor_target( maxval(mpg%strDiv) ) )
#else
            allocate ( summary%local%divertor_plate( maxval(mpg%strDiv) ) )
#endif
            do i = 1, maxval(mpg%strDiv)
#if ( IMAS_MINOR_VERSION > 34 || IMAS_MAJOR_VERSION > 3 )
              call write_sourced_string( summary%local%divertor_target(i)%name, plate_name(i) )
#else
              call write_sourced_string( summary%local%divertor_plate(i)%name, plate_name(i) )
#endif
            end do
#if ( IMAS_MINOR_VERSION > 34 || IMAS_MAJOR_VERSION > 3 )
            call write_errored_value( summary%local%divertor_target(1)%t_e, &
              &  tesepi_av(1), tesepi_std(1) )
            call write_errored_value( summary%local%divertor_target(1)%t_i_average, &
              &  tisepi_av(1), tisepi_std(1) )
            call write_errored_value( summary%local%divertor_target(1)%n_e, &
              &  nesepi_av(1), nesepi_std(1) )
            if (mpg%nncut.eq.2) then
              call write_errored_value( summary%local%divertor_target(2)%t_e, &
                &  tesepi_av(2), tesepi_std(2) )
              call write_errored_value( summary%local%divertor_target(2)%t_i_average, &
                &  tisepi_av(2), tisepi_std(2) )
              call write_errored_value( summary%local%divertor_target(2)%n_e, &
                &  nesepi_av(2), nesepi_std(2) )
              call write_errored_value( summary%local%divertor_target(3)%t_e, &
                &  tesepa_av(1), tesepa_std(1) )
              call write_errored_value( summary%local%divertor_target(3)%t_i_average, &
                &  tisepa_av(1), tisepa_std(1) )
              call write_errored_value( summary%local%divertor_target(3)%n_e, &
                &  nesepa_av(1), nesepa_std(1) )
            end if
            call write_errored_value( summary%local%divertor_target(maxval(mpg%strDiv))%t_e, &
              &  tesepa_av(nc), tesepa_std(nc) )
            call write_errored_value( summary%local%divertor_target(maxval(mpg%strDiv))%t_i_average, &
              &  tisepa_av(nc), tisepa_std(nc) )
            call write_errored_value( summary%local%divertor_target(maxval(mpg%strDiv))%n_e, &
              &  nesepa_av(nc), nesepi_std(nc) )
#else
            call write_errored_value( summary%local%divertor_plate(1)%t_e, &
              &  tesepi_av(1), tesepi_std(1) )
            call write_errored_value( summary%local%divertor_plate(1)%t_i_average, &
              &  tisepi_av(1), tisepi_std(1) )
            call write_errored_value( summary%local%divertor_plate(1)%n_e, &
              &  nesepi_av(1), nesepi_std(1) )
            if (mpg%nncut.eq.2) then
              call write_errored_value( summary%local%divertor_plate(2)%t_e, &
                &  tesepi_av(2), tesepi_std(2) )
              call write_errored_value( summary%local%divertor_plate(2)%t_i_average, &
                &  tisepi_av(2), tisepi_std(2) )
              call write_errored_value( summary%local%divertor_plate(2)%n_e, &
                &  nesepi_av(2), nesepi_std(2) )
              call write_errored_value( summary%local%divertor_plate(3)%t_e, &
                &  tesepa_av(1), tesepa_std(1) )
              call write_errored_value( summary%local%divertor_plate(3)%t_i_average, &
                &  tisepa_av(1), tisepa_std(1) )
              call write_errored_value( summary%local%divertor_plate(3)%n_e, &
                &  nesepa_av(1), nesepa_std(1) )
            end if
            call write_errored_value( summary%local%divertor_plate(maxval(mpg%strDiv))%t_e, &
              &  tesepa_av(nc), tesepa_std(nc) )
            call write_errored_value( summary%local%divertor_plate(maxval(mpg%strDiv))%t_i_average, &
              &  tisepa_av(nc), tisepa_std(nc) )
            call write_errored_value( summary%local%divertor_plate(maxval(mpg%strDiv))%n_e, &
              &  nesepa_av(nc), nesepi_std(nc) )
#endif
          end if
        end if

        call fill_summary_data( geo, summary )
#endif
        call logmsg( LOGDEBUG, "b2mod_ual_io.B25_av_ids: done" )

        return
    end subroutine B25_av_ids

    function separatrix_average( mpg, geo, field, weight )
    ! This function is devoted to obtain the weighted average along the active separatrix
    ! of a plasma field quantity
    ! The average is made using face-centered quantities on the cell faces forming the separatrix
    ! The weighting automatically includes the areas of the cell faces
    implicit none
    type( mapping ), intent(in) :: mpg
    type( geometry ), intent(in) :: geo
    real(kind=IDS_real) :: separatrix_average
    real(kind=IDS_real), intent(in) :: field(mpg%nCv), weight(mpg%nFc)
    integer :: i, iFc, iCv1, iCv2
    real(kind=IDS_real) :: sum, area_sum

    separatrix_average = IDS_REAL_INVALID
    sum = 0.0_IDS_real
    area_sum = 0.0_IDS_real
    do i = mpg%fsFcP(mpg%iFssep,1), &
       &   mpg%fsFcP(mpg%iFssep,1) + mpg%fsFcP(mpg%iFssep,2) - 1
      iFc = mpg%fsFc(i)
      iCv1 = mpg%fcCv(iFc,1)
      iCv2 = mpg%fcCv(iFc,2)
      if ( mpg%cvReg(iCv1).eq.1 .or. mpg%cvReg(iCv2).eq.1 .or. &
        & (mpg%cvReg(iCv1).eq.5 .and. mpg%nnreg(0).eq.8) .or.  &
        & (mpg%cvReg(iCv2).eq.5 .and. mpg%nnreg(0).eq.8) ) then
        sum = sum + geo%fcS(iFc) * weight(iFc) * &
            & ( field(iCv1) + field(iCv2) ) / 2.0_IDS_real
        area_sum = area_sum + geo%fcS(iFc) * weight(iFc)
      end if
    end do
    if (mpg%iFssep2.gt.0) then
      do i = mpg%fsFcP(mpg%iFssep2,1), &
         &   mpg%fsFcP(mpg%iFssep2,1) + mpg%fsFcP(mpg%iFssep2,2) - 1
        iFc = mpg%fsFc(i)
        iCv1 = mpg%fcCv(iFc,1)
        iCv2 = mpg%fcCv(iFc,2)
        if ( mpg%cvReg(iCv1).eq.1 .or. mpg%cvReg(iCv2).eq.1 .or. &
          & (mpg%cvReg(iCv1).eq.5 .and. mpg%nnreg(0).eq.8) .or.  &
          & (mpg%cvReg(iCv2).eq.5 .and. mpg%nnreg(0).eq.8) ) then
          sum = sum + geo%fcS(iFc) * weight(iFc) * &
              & ( field(iCv1) + field(iCv2) ) / 2.0_IDS_real
          area_sum = area_sum + geo%fcS(iFc) * weight(iFc)
        end if
      end do
    end if
    if (area_sum.ne.0.0_IDS_real) separatrix_average = sum / area_sum

    return
    end function separatrix_average

    subroutine write_ids_properties( properties, homo )
    implicit none
    type(ids_ids_properties), intent(inout) :: properties
                !< Type of IDS data structure, designed for IDS properties
    integer, intent(in) :: homo

    properties%homogeneous_time = homo
    allocate( properties%comment(1) )
    properties%comment = comment
#if ( IMAS_MINOR_VERSION > 33 || IMAS_MAJOR_VERSION > 3 )
    allocate( properties%provenance%node(1) )
#if ( IMAS_MINOR_VERSION > 41 || IMAS_MAJOR_VERSION > 3 )
    allocate( properties%provenance%node(1)%reference(1) )
    allocate( properties%provenance%node(1)%reference(1)%name(1) )
    properties%provenance%node(1)%reference(1)%name(1) = source
#else
    allocate( properties%provenance%node(1)%sources(1) )
    properties%provenance%node(1)%sources(1) = source
#endif
#else
    allocate( properties%source(1) )
    properties%source = source
#endif
    allocate( properties%creation_date(1) )
    properties%creation_date = create_date
#if ( IMAS_MINOR_VERSION > 14 || IMAS_MAJOR_VERSION > 3 )
    allocate( properties%provider(1) )
    properties%provider = username
#endif
    return

    end subroutine write_ids_properties

    subroutine write_ids_code( switch, code, commit, description )
    implicit none
    type (switches), intent(in) :: switch
    type(ids_code), intent(inout) :: code
                !< Type of IDS data structure, designed for code data handling
    character(len=ids_string_length), intent(in) :: commit
    character(len=ids_string_length), intent(in) :: description
#if ( IMAS_MINOR_VERSION > 29 || IMAS_MAJOR_VERSION > 3 )
    integer :: nlibs !< Number of declared libraries in IDS description
    character*8 ggd_version, mscl_version
    character*32 SOLPS_git_version
    character*32 get_SOLPS_hash
    character(len=ids_string_length) :: repository
#ifdef B25_EIRENE
    integer p
    character*8 eirene_version
    character*32 Eirene_git_version
    character*32 get_Eir_hash
#endif
#ifdef AMNS
    type (amns_handle_type) :: amns
    type (amns_query_type) :: query
    type (amns_answer_type) :: answer
    type (amns_error_type) :: amns_status
#endif
#endif
    logical streql
    external streql

    allocate( code%name(1) )
    code%name = source
#if ( IMAS_MINOR_VERSION > 38 || IMAS_MAJOR_VERSION > 3 )
    allocate( code%description(1) )
    code%description = description
#endif
    allocate( code%version(1) )
    code%version = newversion
    allocate( code%commit(1) )
    code%commit = commit
    allocate( code%repository(1) )
    code%repository(1) = "ssh://git.iter.org/bnd/b2.5.git"
    allocate( code%output_flag( num_slices ) )
    code%output_flag( slice_index ) = 0

#if ( IMAS_MINOR_VERSION > 29 || IMAS_MAJOR_VERSION > 3 )
    nlibs = 1
#ifdef B25_EIRENE
    if (switch%use_eirene.ne.0) then
      nlibs = nlibs + 1
      Eirene_git_version = get_Eir_hash()
      p = index(Eirene_git_version,'-')
      if (p.eq.0) then
        eirene_version = trim(Eirene_git_version)
      else if (p.gt.1) then
        eirene_version = Eirene_git_version(1:p-1)
      else
        eirene_version = ''
      end if
    end if
#endif
    SOLPS_git_version = get_SOLPS_hash()
    if (.not.streql(SOLPS_git_version,'0.0.0-0-g0000000')) &
      & nlibs = nlibs + 1

    mscl_version='0.0.0'
#ifdef NO_GETENV
    write(ggd_version,'(i1,a1,i2,a1,i1)') GGD_MAJOR_VERSION,'.', &
                                        & GGD_MINOR_VERSION,'.', &
                                        & GGD_MICRO_VERSION
#elif defined(USE_PXFGETENV)
    CALL PXFGETENV ('GGD_VERSION', 0, ggd_version, lenval, ierror)
    CALL PXFGETENV ('EBVERSIONMSCL', 0, mscl_version, lenval, ierror)
#else
    call get_environment_variable('GGD_VERSION', &
        &  status=ierror,length=lenval)
    if (ierror.eq.0) call get_environment_variable('GGD_VERSION', &
        &  value=ggd_version)
    call get_environment_variable('EBVERSIONMSCL', &
        &  status=ierror,length=lenval)
    if (ierror.eq.0) call get_environment_variable('EBVERSIONMSCL', &
        &  value=mscl_version)
#endif
    if (.not.streql(mscl_version,'0.0.0')) nlibs = nlibs + 1

    allocate( code%library( nlibs ) )

    nlibs = 1
    allocate( code%library( nlibs )%name(1) )
    code%library( nlibs )%name = 'GGD'
    allocate( code%library( nlibs )%version(1) )
    code%library( nlibs )%version = ggd_version
    allocate( code%library( nlibs )%repository(1) )
    repository = "ssh://git.iter.org/imex/ggd.git"
    code%library( nlibs )%repository = repository
#ifdef B25_EIRENE
    if (switch%use_eirene.ne.0) then
      nlibs = nlibs + 1
      allocate( code%library( nlibs )%name(1) )
      code%library( nlibs )%name = 'EIRENE'
      allocate( code%library( nlibs )%version(1) )
      code%library( nlibs )%version = eirene_version
      allocate( code%library( nlibs )%commit(1) )
      code%library( nlibs )%commit = Eirene_git_version
      allocate( code%library( nlibs )%repository(1) )
      repository = "ssh://git.iter.org/bnd/eirene.git"
      code%library( nlibs )%repository = repository
    end if
#endif
    if (.not.streql(SOLPS_git_version,'0.0.0-0-g0000000')) then
      nlibs = nlibs + 1
      allocate( code%library( nlibs )%name(1) )
      code%library( nlibs )%name = 'SOLPS-ITER'
      allocate( code%library( nlibs )%version(1) )
      code%library( nlibs )%version = newversion
      allocate( code%library( nlibs )%commit(1) )
      code%library( nlibs )%commit = SOLPS_git_version
      allocate( code%library( nlibs )%repository(1) )
      repository = "ssh://git.iter.org/bnd/solps-iter.git"
      code%library( nlibs )%repository = repository
    end if
    if (.not.streql(mscl_version,'0.0.0')) then
      nlibs = nlibs + 1
      allocate( code%library( nlibs )%name(1) )
      code%library( nlibs )%name = 'MSCL'
      allocate( code%library( nlibs )%version(1) )
      code%library( nlibs )%version = mscl_version
      allocate( code%library( nlibs )%commit(1) )
      code%library( nlibs )%commit = SOLPS_git_version
      allocate( code%library( nlibs )%repository(1) )
      repository = "ssh://git.iter.org/lib/mscl.git"
      code%library( nlibs )%repository = repository
    end if
#endif
    return
    end subroutine write_ids_code

#if IMAS_MAJOR_VERSION > 3
    subroutine write_ids_code_constant( code, commit, description )
    implicit none
    type(ids_code_constant), intent(inout) :: code
                !< Type of IDS data structure, designed for code data handling
    character(len=ids_string_length), intent(in) :: commit
    character(len=ids_string_length), intent(in) :: description

    allocate( code%name(1) )
    code%name = source
    allocate( code%description(1) )
    code%description = description
    allocate( code%version(1) )
    code%version = newversion
    allocate( code%commit(1) )
    code%commit = commit
    allocate( code%repository(1) )
    code%repository(1) = "ssh://git.iter.org/bnd/b2.5.git"

    return
    end subroutine write_ids_code_constant
#endif

    subroutine put_equilibrium_data ( mpg, geo, equilibrium, &
#if ( IMAS_MINOR_VERSION > 21 || IMAS_MAJOR_VERSION > 3 )
       &  summary, &
#endif
       &  edgeprof, database, &
#if ( AL_MAJOR_VERSION > 4 && GGD_MAJOR_VERSION > 0 )
       &  time_sind, &
#endif
       &  time_slice_value, do_summary_data, new_eq_ggd )
#if ( IMAS_MINOR_VERSION > 14 || IMAS_MAJOR_VERSION > 3 ) && GGD_MAJOR_VERSION > 0
    use b2mod_ual_io_grid &
       & , only: GGD_copy_AoS3Root_to_Dynamic
#endif
    implicit none
    type (mapping) :: mpg !< The grid mapping
    type (geometry) :: geo !< The grid geometry
    type (ids_equilibrium) :: equilibrium !< IDS designed to store
            !< equilibrium data
#if ( IMAS_MINOR_VERSION > 21 || IMAS_MAJOR_VERSION > 3 )
    type (ids_summary) :: summary !< IDS designed to store
            !< run summary data
#endif
    type (ids_edge_profiles) :: edgeprof !< IDS designed to store
            !< edge profiles data
    character(len=24), intent(in) :: database
#if ( AL_MAJOR_VERSION > 4 && GGD_MAJOR_VERSION > 0 )
    integer, intent(in) :: time_sind     !< Corresponding time slice index
                                         !< in edge_profiles IDS
#endif
    real(IDS_real), intent(in) :: time_slice_value   !< Time slice value
    logical, intent(in) :: do_summary_data
    logical, intent(out) :: new_eq_ggd
#if GGD_MAJOR_VERSION > 0
#if ( IMAS_MINOR_VERSION < 15 && IMAS_MAJOR_VERSION < 4 )
    type(ids_generic_grid_dynamic) :: eq_grid !< Type of IDS
        !< data structure, designed for handling equilibrium grid geometry data
#else
    type(ids_generic_grid_aos3_root) :: eq_grid !< Type of IDS
        !< data structure, designed for handling equilibrium grid geometry data
#endif
    real(IDS_real) :: tmpVx( mpg%nVx )
    real(IDS_real) :: tmpFace( mpg%nFc )
    real(IDS_real) :: tmpCv( mpg%nCv )
#endif
    integer :: i, ic, iCv, iCv1, iCv2, iFc, iFt, iVx, inode, icrmax
    integer :: idum(0:3)
    integer, save :: ncall = 0
    real(IDS_real) :: parg(0:99)
    real(IDS_real), save :: pit_rescale = 1.0_IDS_real
    real(IDS_real) :: b0r0_ref, z_eq, z_min, z_max, r_max, rFc
    character*8 id
    character*80 cnamip, cvalip
    character*132 eq_source
    character*256 filename
    character*500 line, ligne
    logical exists
    logical is_comment, streql
    external is_comment, streql
    external b2agx0, find_file, ipgetr, strip_spaces

    eq_found = .false.
    new_eq_ggd = .false.
    eq_source = "ITER Baseline q95=3 equilibrium"
    if ( associated( equilibrium%time_slice ) ) then
      if ( size( equilibrium%time_slice ).ge.slice_index ) then
        eq_found = .true.
#if ( IMAS_MAJOR_VERSION > 3 || IMAS_MINOR_VERSION > 33 )
        if ( associated( equilibrium%ids_properties%provenance%node ) ) then
#if ( IMAS_MAJOR_VERSION > 3 || IMAS_MINOR_VERSION > 41 )
          if ( associated( equilibrium%ids_properties%provenance%node(1)%reference ) ) &
            & eq_source = equilibrium%ids_properties%provenance%node(1)%reference(1)%name(1)
#else
          if ( associated( equilibrium%ids_properties%provenance%node(1)%sources ) ) &
            & eq_source = equilibrium%ids_properties%provenance%node(1)%sources(1)
#endif
#if ( IMAS_MAJOR_VERSION == 3 && IMAS_MINOR_VERSION > 33 )
        else if ( associated( equilibrium%ids_properties%source ) ) then
          eq_source = equilibrium%ids_properties%source(1)
#endif
        end if
#else
        if ( associated( equilibrium%ids_properties%source ) ) &
           & eq_source = equilibrium%ids_properties%source(1)
#endif
        r0 = equilibrium%vacuum_toroidal_field%r0
        b0 = equilibrium%vacuum_toroidal_field%b0( slice_index )
        b0r0 = b0 * r0
      end if
    end if
    if (.not.eq_found) then
      r0 = 0.0_R8
      i = 1
      do while ( i.le.mpg%nCi .and. .not.mpg%cvOnClosedSurface(i) )
        i = i + 1
      end do
      if ( i.le.mpg%nCi ) then
        iFt = mpg%cvFt(i)
        do ic = mpg%ftCvP(iFt,1), mpg%ftCvP(iFt,1)+mpg%ftCvP(iFt,2)-1
          iCv = mpg%ftCv(ic)
          if (isymm.eq.1.or.isymm.eq.2) then
            r0 = r0 + geo%cvX(iCv)
          else if (isymm.eq.3.or.isymm.eq.4) then
            r0 = r0 + geo%cvY(iCv)
          end if
        end do
        r0 = r0 / float(mpg%ftCvP(iFt,2))
      end if
      if (r0.gt.0.0_R8) then
        if ( geo%vxFfbz(1).ne.0.0_R8) then
          if (isymm.eq.0) then
            b0r0 = geo%vxFfbz(1)
          else
            b0r0 = geo%vxFfbz(1)/(2.0_R8 * pi)
          end if
          b0 = b0r0 / r0
        else if (isymm.eq.0) then
          b0 = sign(geo%cvBb(1,2),geo%cvEb(1,2))
          b0r0 = b0*r0
        else if (isymm.eq.1 .or. isymm.eq.2) then
          b0r0 = sign(geo%cvBb(1,2),geo%cvEb(1,2))*geo%cvX(1)
          b0 = b0r0 / r0
        else if (isymm.eq.3 .or. isymm.eq.4) then
          b0r0 = sign(geo%cvBb(1,2),geo%cvEb(1,2))*geo%cvY(1)
          b0 = b0r0 / r0
        end if
      else
        b0 = sign(geo%cvBb(1,2),geo%cvEb(1,2))
        if (isymm.eq.1 .or. isymm.eq.2) then
          b0r0 = sign(geo%cvBb(1,2),geo%cvEb(1,2))*geo%cvX(1)
        else if (isymm.eq.3 .or. isymm.eq.4) then
          b0r0 = sign(geo%cvBb(1,2),geo%cvEb(1,2))*geo%cvY(1)
        end if
        if (b0.ne.0.0_R8) then
          r0 = b0r0 / b0
        else
          r0 = IDS_REAL_INVALID
        end if
      end if
    end if

    if (ncall.eq.0) then
      call ipgetr ('b2agfs_pit_rescale', pit_rescale)
      if (pit_rescale.eq.1.0_R8) then
        filename='b2ag.dat'
        call find_file(filename,exists)
        if (exists) then
          open(99,file=filename)
          call b2agx0 (99, idum(0), idum(1), idum(2))
          read (99,'(a8)',err=2) id
          read (99,*,err=2) parg
    1     continue
          read (99,'(a)',end=2,err=2) line
          if (.not.is_comment(line)) then
            ligne = line
            call strip_spaces(ligne)
            if (ligne(1:1).eq.'''') then
              read (line,*) cnamip, cvalip
              call ipsetc (cnamip, cvalip)
            endif
          endif
          goto 1
    2     continue
          close(99)
          call ipgetr ('b2agfs_pit_rescale', pit_rescale)
        end if
      end if
    end if

    !> Careful: Sign convention for magnetic field in IDS
    !>          is OPPOSITE to that in SOLPS toroidal geometries
    if ( b0.ne.0.0_IDS_real ) then
      if (streql(database,'ITER').and..not.eq_found) then
        b0r0_ref = 5.3_IDS_real * 6.2_IDS_real
        allocate( edgeprof%vacuum_toroidal_field%b0( num_slices ) )
        if ( abs(pit_rescale).eq.1.0_IDS_real ) then
          i = nint(b0r0_ref/b0r0)
          select case (i)
          case (1)
#if ( IMAS_MINOR_VERSION > 21 || IMAS_MAJOR_VERSION > 3 )
            if (do_summary_data) then
              call write_sourced_value( summary%global_quantities%ip, &
                  & -15.0e6_IDS_real )
              call write_sourced_value( summary%global_quantities%b0, &
                  & -5.3_IDS_real )
            endif
#endif
            edgeprof%vacuum_toroidal_field%b0( slice_index ) = -5.3_IDS_real
          case (2)
#if ( IMAS_MINOR_VERSION > 21 || IMAS_MAJOR_VERSION > 3 )
            if (do_summary_data) then
              call write_sourced_value( summary%global_quantities%ip, &
                  & -7.5e6_IDS_real )
              call write_sourced_value( summary%global_quantities%b0, &
                  & -2.65_IDS_real )
            end if
#endif
            edgeprof%vacuum_toroidal_field%b0( slice_index ) = -2.65_IDS_real
          case (3)
#if ( IMAS_MINOR_VERSION > 21 || IMAS_MAJOR_VERSION > 3 )
            if (do_summary_data) then
              call write_sourced_value( summary%global_quantities%ip, &
                  & -5.0e6_IDS_real )
              call write_sourced_value( summary%global_quantities%b0, &
                  & -1.8_IDS_real )
            endif
#endif
            edgeprof%vacuum_toroidal_field%b0( slice_index ) = -1.8_IDS_real
          case default
#if ( IMAS_MINOR_VERSION > 21 || IMAS_MAJOR_VERSION > 3 )
            if (do_summary_data) then
              call write_sourced_value( summary%global_quantities%ip, &
                  & -15.0e6_IDS_real/nint(b0r0_ref/b0r0) )
              call write_sourced_value( summary%global_quantities%b0, &
                  & -b0r0 / 6.2e6_IDS_real )
            endif
#endif
            edgeprof%vacuum_toroidal_field%b0( slice_index ) = -b0r0 / 6.2_IDS_real
          end select
#if ( IMAS_MINOR_VERSION > 21 || IMAS_MAJOR_VERSION > 3 )
          if (do_summary_data) then
            summary%global_quantities%b0%source = &
                  & "ITER Baseline q95=3 equilibrium"
            summary%global_quantities%ip%source = &
                  & "ITER Baseline q95=3 equilibrium"
            call write_sourced_value( summary%global_quantities%q_95, &
                  & 3.0_IDS_real )
            summary%global_quantities%q_95%source = &
                  & "ITER Baseline q95=3 equilibrium"
          endif
#endif
        else
          edgeprof%vacuum_toroidal_field%b0( slice_index ) = -b0r0 / 6.2_IDS_real
#if ( IMAS_MINOR_VERSION > 21 || IMAS_MAJOR_VERSION > 3 )
          if (do_summary_data) then
            call write_sourced_value( summary%global_quantities%b0, &
                  & -b0r0 / 6.2_IDS_real )
            write(eq_source, '(a,1pe12.5,a)' ) &
                  & "ITER Baseline q95=3 equilibrium"// &
                  &  " (current rescaled by ",abs(pit_rescale),")"
            summary%global_quantities%b0%source = eq_source
            call write_sourced_value( summary%global_quantities%ip, &
                  &  -15.0e6_IDS_real*abs(pit_rescale) )
            summary%global_quantities%ip%source = eq_source
            call write_sourced_value( summary%global_quantities%q_95, &
                  &   3.0_IDS_real/abs(pit_rescale) )
            summary%global_quantities%q_95%source = eq_source
          endif
#endif
        end if
#if ( IMAS_MINOR_VERSION > 21 || IMAS_MAJOR_VERSION > 3 )
        if (do_summary_data) then
          call write_sourced_constant( summary%global_quantities%r0, &
                  & 6.2_IDS_real )
          summary%global_quantities%r0%source = &
                  & "ITER Baseline q95=3 equilibrium"
        end if
#endif
        edgeprof%vacuum_toroidal_field%r0 = 6.2_IDS_real
      else
        allocate( edgeprof%vacuum_toroidal_field%b0( num_slices ) )
        if (eq_found) then
          edgeprof%vacuum_toroidal_field%b0( slice_index ) = b0
          edgeprof%vacuum_toroidal_field%r0 = r0
#if ( IMAS_MINOR_VERSION > 21 || IMAS_MAJOR_VERSION > 3 )
          if (do_summary_data) then
            call write_sourced_value( summary%global_quantities%b0, b0 )
            summary%global_quantities%b0%source = eq_source
            if ( equilibrium%time_slice( slice_index )%global_quantities%  &
               & ip .ne. IDS_REAL_INVALID ) then
              call write_sourced_value( summary%global_quantities%ip,    &
                &  equilibrium%time_slice( slice_index )%global_quantities%ip )
              summary%global_quantities%ip%source = eq_source
            end if
            if ( equilibrium%time_slice( slice_index )%global_quantities%  &
               & q_95 .ne. IDS_REAL_INVALID ) then
              call write_sourced_value( summary%global_quantities%q_95,  &
                &  equilibrium%time_slice( slice_index )%global_quantities%q_95 )
              summary%global_quantities%q_95%source = eq_source
            else if (streql(eq_source,"ITER Baseline q95=3 equilibrium")) then
              call write_sourced_value( summary%global_quantities%q_95,  &
                &  3.0_IDS_real )
              summary%global_quantities%q_95%source = eq_source
            end if
            call write_sourced_constant( summary%global_quantities%r0, r0 )
            summary%global_quantities%r0%source = eq_source
          end if
#endif
#if ( IMAS_MINOR_VERSION > 11 || IMAS_MAJOR_VERSION > 3 ) && GGD_MAJOR_VERSION > 0
          new_eq_ggd = .not.associated( equilibrium%grids_ggd )
          if ( .not.new_eq_ggd ) new_eq_ggd = &
            &  .not.associated( equilibrium%grids_ggd( slice_index )%grid )
#else
          new_eq_ggd = .false.
#endif
          if ( new_eq_ggd ) then
#if ( IMAS_MINOR_VERSION > 11 || IMAS_MAJOR_VERSION > 3 ) && GGD_MAJOR_VERSION > 0
            if (.not.associated( equilibrium%grids_ggd ) ) &
              &  allocate( equilibrium%grids_ggd( num_time_slices ) )
            allocate( equilibrium%grids_ggd( slice_index )%grid(1) )
            call b2_IMAS_Fill_Grid_Desc( mpg, geo, eq_grid )
#if ( IMAS_MINOR_VERSION > 14 || IMAS_MAJOR_VERSION > 3 )
#if AL_MAJOR_VERSION > 4
            allocate( equilibrium%grids_ggd( slice_index )%grid(1)%path(1) )
            equilibrium%grids_ggd( slice_index )%grid(1)%path = &
              &   "#edge_profiles/grid_ggd("//int2str(time_sind)//")"
#else
            call GGD_copy_AoS3Root_to_Dynamic( eq_grid, &
              &   equilibrium%grids_ggd( slice_index )%grid(1) )
#endif
#else
            equilibrium%grids_ggd( slice_index )%grid(1) = eq_grid
#endif
            equilibrium%grids_ggd( slice_index )%time = time_slice_value
#endif
#if ( IMAS_MINOR_VERSION > 33 || IMAS_MAJOR_VERSION > 3 )
            if (.not.associated( equilibrium%ids_properties%provenance%node ) ) then
              inode = 0
            else
              inode = size( equilibrium%ids_properties%provenance%node )
            endif
            allocate( equilibrium%ids_properties%provenance%node(inode + 1) )
            allocate( &
               & equilibrium%ids_properties%provenance%node(inode+1)%path(1) )
            equilibrium%ids_properties%provenance%node(inode+1)%path =        &
               & "grids_ggd"
#if ( IMAS_MINOR_VERSION > 41 || IMAS_MAJOR_VERSION > 3 )
            allocate( &
               & equilibrium%ids_properties%provenance%node(inode+1)%reference(1) )
            allocate( &
               & equilibrium%ids_properties%provenance%node(inode+1)%reference(1)%name(1) )
            equilibrium%ids_properties%provenance%node(inode+1)%reference(1)%name(1) =  &
               &  source
#else
            allocate( &
               & equilibrium%ids_properties%provenance%node(inode+1)%sources(1) )
            equilibrium%ids_properties%provenance%node(inode+1)%sources(1) =  &
               &  source
#endif
#endif
          end if
          if (.not.associated( equilibrium%time_slice )) then
            allocate( equilibrium%time_slice( num_time_slices ) )
          end if
          if ( equilibrium%time_slice( slice_index )%time.eq.IDS_REAL_INVALID ) &
            &  equilibrium%time_slice( slice_index )%time = time_slice_value
          if (.not.associated( equilibrium%time_slice( slice_index )%ggd ) )  &
            & then
            allocate( equilibrium%time_slice( slice_index )%ggd(1) )
          end if
#if ( IMAS_MINOR_VERSION > 11 || IMAS_MAJOR_VERSION > 3 ) && GGD_MAJOR_VERSION > 0
          if (.not.associated(                                                &
            &  equilibrium%time_slice( slice_index )%ggd(1)%r ) ) then
            call write_vertex_scalar( eq_grid, mpg,                           &
                &   scalar = equilibrium%time_slice( slice_index )%ggd(1)%r,  &
                &   b2VertexData = geo%vxX )
            tmpFace(:) = ( geo%vxX(mpg%fcVx(:,1)) + geo%vxX(mpg%fcVx(:,2)) )/2.0_IDS_real
            call write_face_scalar( eq_grid, mpg,                             &
                &   val = equilibrium%time_slice( slice_index )%ggd(1)%r,     &
                &   value = tmpFace )
            call write_cell_scalar( eq_grid, mpg,                             &
                &   scalar = equilibrium%time_slice( slice_index )%ggd(1)%r,  &
                &   b2CellData = geo%cvX )
          end if
          if (.not.associated(                                                &
            &  equilibrium%time_slice( slice_index )%ggd(1)%z ) ) then
            call write_vertex_scalar( eq_grid, mpg,                           &
                &   scalar = equilibrium%time_slice( slice_index )%ggd(1)%z,  &
                &   b2VertexData = geo%vxY )
            tmpFace(:) = ( geo%vxY(mpg%fcVx(:,1)) + geo%vxY(mpg%fcVx(:,2)) )/2.0_IDS_real
            call write_face_scalar( eq_grid, mpg,                             &
                &   val = equilibrium%time_slice( slice_index )%ggd(1)%z,     &
                &   value = tmpFace )
            call write_cell_scalar( eq_grid, mpg,                             &
                &   scalar = equilibrium%time_slice( slice_index )%ggd(1)%z,  &
                &   b2CellData = geo%cvY )
          end if
          if (maxval(abs(geo%vxFpsi)).ne.0.0_IDS_real .and. .not.associated(        &
            &  equilibrium%time_slice( slice_index )%ggd(1)%psi ) ) then
            call write_vertex_scalar( eq_grid, mpg,                           &
                &   scalar = equilibrium%time_slice( slice_index )%ggd(1)%    &
                &            psi,                                             &
                &   b2VertexData = geo%vxFpsi )
            tmpFace(:) = ( geo%vxFpsi(mpg%fcVx(:,1)) + geo%vxFpsi(mpg%fcVx(:,2)) )/2.0_IDS_real
            call write_face_scalar( eq_grid, mpg,                             &
                &   val = equilibrium%time_slice( slice_index )%ggd(1)%psi,   &
                &   value = tmpFace )
#ifdef WG_TODO
            tmpCv(:,:) = (fpsi(:,:,0) + fpsi(:,:,1) +                         &
                &         fpsi(:,:,2) + fpsi(:,:,3) )/4.0_IDS_real
            call write_cell_scalar( eq_grid, mpg,                             &
                &   scalar = equilibrium%time_slice( slice_index )%ggd(1)%    &
                &            psi,   &
                &   b2CellData = tmpCv )
#endif
          end if
          if (maxval(abs(geo%vxFfbz)).ne.0.0_IDS_real .and. .not.associated(        &
            &  equilibrium%time_slice( slice_index )%ggd(1)%phi ) ) then
            tmpFace(:) = ( geo%vxFfbz(mpg%fcVx(:,1)) + geo%vxFfbz(mpg%fcVx(:,2)) )/2.0_IDS_real
            tmpCv = 0.0_IDS_real
            do iCv = 1, mpg%nCv
              do i = mpg%cvVxP(iCv,1), mpg%cvVxP(iCv,1)+mpg%cvVxP(iCv,2)-1
                iVx = mpg%cvVx(iCv)
                tmpCv(iCv) = tmpCv(iCv) + geo%vxFFbz(iVx)
              end do
              tmpCv(iCv) = tmpCv(iCv) / float( mpg%cvVxP(iCv,2) )
            end do
            call write_vertex_scalar( eq_grid, mpg,                           &
                &   scalar = equilibrium%time_slice( slice_index )%ggd(1)%    &
                &            phi,                                             &
                &   b2VertexData = geo%vxFfbz )
            call write_face_scalar( eq_grid, mpg,                             &
                &   val = equilibrium%time_slice( slice_index )%ggd(1)%phi,   &
                &   value = tmpFace )
            call write_cell_scalar( eq_grid, mpg,                             &
                &   scalar = equilibrium%time_slice( slice_index )%ggd(1)%    &
                &            phi,                                             &
                &   b2CellData = tmpCv )
          end if
          if (.not.associated(                                                &
            &  equilibrium%time_slice( slice_index )%ggd(1)%b_field_r ) ) then
            tmpVx(:) = geo%vxBb(:,0)*geo%vxEb(:,0)
            tmpFace(:) = geo%fcBb(:,0)*geo%fcEb(:,0)
            tmpCv(:) = geo%cvBb(:,0)*geo%cvEb(:,0)
            call write_vertex_scalar( eq_grid, mpg,                           &
                &   scalar = equilibrium%time_slice( slice_index )%ggd(1)%    &
                &         b_field_r,                                          &
                &   b2VertexData = tmpVx )
            call write_face_scalar( eq_grid, mpg,                             &
                &   val = equilibrium%time_slice( slice_index )%ggd(1)%       &
                &         b_field_r,                                          &
                &   value = tmpFace )
            call write_cell_scalar( eq_grid, mpg,                             &
                &   scalar = equilibrium%time_slice( slice_index )%ggd(1)%    &
                &         b_field_r,                                          &
                &   b2CellData = tmpCv )
          end if
          if (.not.associated(                                                &
            &  equilibrium%time_slice( slice_index )%ggd(1)%b_field_z ) ) then
            tmpVx(:) = geo%vxBb(:,0)*geo%vxEb(:,1)
            tmpFace(:) = geo%fcBb(:,0)*geo%fcEb(:,1)
            tmpCv(:) = geo%cvBb(:,0)*geo%cvEb(:,1)
            call write_vertex_scalar( eq_grid, mpg,                           &
                &   scalar = equilibrium%time_slice( slice_index )%ggd(1)%    &
                &         b_field_z,                                          &
                &   b2VertexData = tmpVx )
            call write_face_scalar( eq_grid, mpg,                             &
                &   val = equilibrium%time_slice( slice_index )%ggd(1)%       &
                &         b_field_z,                                          &
                &   value = tmpFace )
            call write_cell_scalar( eq_grid, mpg,                             &
                &   scalar = equilibrium%time_slice( slice_index )%ggd(1)%    &
                &         b_field_z,                                          &
                &   b2CellData = tmpCv )
          end if
#if ( IMAS_MINOR_VERSION > 41 || IMAS_MAJOR_VERSION > 3 )
          if (.not.associated(                                                &
            &  equilibrium%time_slice( slice_index )%ggd(1)%b_field_phi ) ) then
            tmpVx(:) = sign(geo%vxBb(:,2),-geo%vxEb(:,2))
            tmpFace(:) = sign(geo%fcBb(:,2),-geo%fcEb(:,2))
            tmpCv(:) = sign(geo%cvBb(:,2),-geo%cvEb(:,2))
            call write_vertex_scalar( eq_grid, mpg,                           &
                &   scalar = equilibrium%time_slice( slice_index )%ggd(1)%    &
                &         b_field_phi,                                        &
                &   b2VertexData = tmpVx )
            call write_face_scalar( eq_grid, mpg,                             &
                &   val = equilibrium%time_slice( slice_index )%ggd(1)%       &
                &         b_field_phi,                                        &
                &   value = tmpFace )
            call write_cell_scalar( eq_grid, mpg,                             &
                &   scalar = equilibrium%time_slice( slice_index )%ggd(1)%    &
                &         b_field_phi,                                        &
                &   b2CellData = tmpCv )
          end if
#else
          if (.not.associated(                                                &
            &  equilibrium%time_slice( slice_index )%ggd(1)%b_field_tor ) ) then
            tmpVx(:) = sign(geo%vxBb(:,2),-geo%vxEb(:,2))
            tmpFace(:) = sign(geo%fcBb(:,2),-geo%fcEb(:,2))
            tmpCv(:) = sign(geo%cvBb(:,2),-geo%cvEb(:,2))
            call write_vertex_scalar( eq_grid, mpg,                           &
                &   scalar = equilibrium%time_slice( slice_index )%ggd(1)%    &
                &         b_field_tor,                                        &
                &   b2VertexData = tmpVx )
            call write_face_scalar( eq_grid, mpg,                             &
                &   val = equilibrium%time_slice( slice_index )%ggd(1)%       &
                &         b_field_tor,                                        &
                &   value = tmpFace )
            call write_cell_scalar( eq_grid, mpg,                             &
                &   scalar = equilibrium%time_slice( slice_index )%ggd(1)%    &
                &         b_field_tor,                                        &
                &   b2CellData = tmpCv )
          end if
#endif
#endif
          if ( equilibrium%time( slice_index ).eq.0.0_IDS_real ) then
            equilibrium%time( slice_index ) = time_slice_value
          end if
        else if (isymm.ne.0) then
#if ( IMAS_MINOR_VERSION > 21 || IMAS_MAJOR_VERSION > 3 )
          if (do_summary_data) then
            call write_sourced_value( summary%global_quantities%b0, -b0 )
            call write_sourced_constant( summary%global_quantities%r0, r0 )
          end if
#endif
          edgeprof%vacuum_toroidal_field%b0( slice_index ) = -b0
          edgeprof%vacuum_toroidal_field%r0 = r0
        else
#if ( IMAS_MINOR_VERSION > 21 || IMAS_MAJOR_VERSION > 3 )
          if (do_summary_data) then
            call write_sourced_value( summary%global_quantities%b0, b0 )
            call write_sourced_constant( summary%global_quantities%r0, r0 )
          end if
#endif
          edgeprof%vacuum_toroidal_field%b0( slice_index ) = b0
          edgeprof%vacuum_toroidal_field%r0 = r0
        end if
      end if
    end if

    if (GeometryType .eq. GEOMETRY_LINEAR) then
      midplane_id = 4
    else
      if ( eq_found ) then
        z_eq = equilibrium%time_slice( slice_index )%global_quantities%  &
             &   magnetic_axis%z
      else
        z_eq = IDS_REAL_INVALID
      end if
      iVx = mpg%cvVx(mpg%cvVxP(omp(icsepomp-1),1))
      z_min = geo%vxY(iVx)
      z_max = geo%vxY(iVx)
      do i = mpg%cvVxP(omp(icsepomp-1),1) + 1, &
     &       mpg%cvVxP(omp(icsepomp-1),1) + mpg%cvVxP(omp(icsepomp-1),2) - 1
        iVx = mpg%cvVx(i)
        z_min = min(z_min,geo%vxY(iVx))
        z_max = max(z_max,geo%vxY(iVx))
      end do
      iFc = mpg%fsFc(mpg%fsFcP(mpg%iFssep,1))
      r_max = ( geo%vxX(mpg%fcVx(iFc,1)) + geo%vxX(mpg%fcVx(iFc,2)) ) / 2.0_R8
      do i = mpg%fsFcP(mpg%iFssep,1) + 1, &
     &       mpg%fsFcP(mpg%iFssep,1) + mpg%fsFcP(mpg%iFssep,2) - 1
        iFc = mpg%fsFc(i)
        rFc = ( geo%vxX(mpg%fcVx(iFc,1)) + geo%vxX(mpg%fcVx(iFc,2)) ) / 2.0_R8
        iCv1 = mpg%fcCv(iFc,1)
        iCv2 = mpg%fcCv(iFc,2)
        iCv = 0
        if ( mpg%cvReg(iCv1).eq.1 .or. &
     &     ( mpg%cvReg(iCv1).eq.5 .and. mpg%nnreg(0).eq.8 ) ) iCv = iCv1
        if ( mpg%cvReg(iCv2).eq.1 .or. &
     &     ( mpg%cvReg(iCv2).eq.5 .and. mpg%nnreg(0).eq.8 ) ) iCv = iCv2
        if ( iCv.eq.0 ) cycle ! skip the divertor legs
        if ( rFc .gt. r_max) then
          r_max = rFc
          icrmax = iCv
        end if
      end do
      if (mpg%iFssep2.gt.0) then
        iFc = mpg%fsFc(mpg%fsFcP(mpg%iFssep2,1))
        do i = mpg%fsFcP(mpg%iFssep2,1), &
     &         mpg%fsFcP(mpg%iFssep2,1) + mpg%fsFcP(mpg%iFssep2,2) - 1
          iFc = mpg%fsFc(i)
          rFc = ( geo%vxX(mpg%fcVx(iFc,1)) + geo%vxX(mpg%fcVx(iFc,2)) ) / 2.0_R8
          iCv1 = mpg%fcCv(iFc,1)
          iCv2 = mpg%fcCv(iFc,2)
          iCv = 0
          if ( mpg%cvReg(iCv1).eq.1 .or. &
     &       ( mpg%cvReg(iCv1).eq.5 .and. mpg%nnreg(0).eq.8 ) ) iCv = iCv1
          if ( mpg%cvReg(iCv2).eq.1 .or. &
     &       ( mpg%cvReg(iCv2).eq.5 .and. mpg%nnreg(0).eq.8 ) ) iCv = iCv2
          if ( iCv.eq.0 ) cycle ! skip the divertor legs
          if ( rFc .gt. r_max) then
            r_max = rFc
            icrmax = iCv
          end if
        end do
      end if
      if ( z_eq.ne.IDS_REAL_INVALID .and. &
         & z_min.le.z_eq .and. z_max.ge.z_eq ) then
        midplane_id = 1
      else if ( omp(icsepomp-1) .eq. icrmax ) then
        midplane_id = 2
      else if ( z_min*z_max .le. 0.0_R8 ) then
        midplane_id = 3
      else
        midplane_id = 4
      end if
    end if

    ncall = ncall + 1
    return
    end subroutine put_equilibrium_data

#if ( IMAS_MINOR_VERSION > 21 || IMAS_MAJOR_VERSION > 3 )
    subroutine fill_summary_data( geo, summary )
    implicit none
    type (geometry), intent(in) :: geo
    type (ids_summary), intent(inout) :: summary

    select case (GeometryType)
    case( GEOMETRY_LIMITER )
      call write_sourced_integer( summary%boundary%type, 0 )
    case( GEOMETRY_SN )
      if ( isymm.ne.0 ) then
        if ( geo%LSN ) then
          call write_sourced_integer( summary%boundary%type, 11 )
        else
          call write_sourced_integer( summary%boundary%type, 12 )
        end if
      else
        call write_sourced_integer( summary%boundary%type, 1 )
      end if
    case( GEOMETRY_CDN , GEOMETRY_DDN_BOTTOM , GEOMETRY_DDN_TOP )
      call write_sourced_integer( summary%boundary%type, 13 )
    case( GEOMETRY_LFS_SNOWFLAKE_MINUS , &
        & GEOMETRY_LFS_SNOWFLAKE_PLUS )
      call write_sourced_integer( summary%boundary%type, 14 )
    end select

    call write_sourced_value( summary%fusion%power, fusion_power*1.0e6_IDS_real )

    return
    end subroutine fill_summary_data
#endif

#if ( IMAS_MINOR_VERSION > 32 || IMAS_MAJOR_VERSION > 3 )
    subroutine write_ids_midplane( midplane, midplane_id )
    implicit none
    type(ids_identifier_static) :: midplane
    integer, intent(in) :: midplane_id

    midplane%index = midplane_id
    allocate( midplane%name(1) )
    allocate( midplane%description(1) )
#if ( IMAS_MINOR_VERSION > 30 || IMAS_MAJOR_VERSION > 3 )
    midplane%name = midplane_identifier%name( midplane_id )
    midplane%description = midplane_identifier%description( midplane_id )
#else
    select case (midplane_id)
    case (1)
      midplane%name = 'magnetic_axis'
      midplane%description = 'Height of equilibrium O-point'
    case (2)
      midplane%name = 'dr_dz_zero_sep'
      midplane%description = 'Maximum radius location along separatrix'
    case (3)
      midplane%name = 'z_zero'
      midplane%description = 'Z = 0 plane'
    case (4)
      midplane%name = 'ggd_subset'
      midplane%description = &
         &  'Location specified by GGD outer midplane grid subset'
    end select
#endif
    return

    end subroutine write_ids_midplane
#endif

#if ( IMAS_MINOR_VERSION > 11 || IMAS_MAJOR_VERSION > 3 ) && GGD_MAJOR_VERSION > 0
    !> Write scalar B2 cell quantity to 'ids_generic_grid_scalar'
    !! IMAS IDS data tree node.
    subroutine write_IDS_quantity( basegrid, mpg, geo, val, value )
    use b2mod_interp
    implicit none
#if ( IMAS_MINOR_VERSION < 15 && IMAS_MAJOR_VERSION < 4 )
    type(ids_generic_grid_dynamic), intent(in) :: basegrid !< Type of IDS
        !< data structure, designed for handling grid geometry data
#else
    type(ids_generic_grid_aos3_root), intent(in) :: basegrid !< Type of IDS
        !< data structure, designed for handling grid geometry data
#endif
    type(mapping), intent(in) :: mpg
    type(geometry), intent(in) :: geo
    type(ids_generic_grid_scalar), pointer, intent(inout) :: val(:)
        !< Type of IDS data structure, designed for scalar data handling
    real(IDS_real), intent(in) :: value( mpg%nCv )
    real(IDS_real), dimension(:), pointer :: idsdata    !< Array for
        !< handling data field values
    real(IDS_real) :: tmpFace( mpg%nFc, 0:1)
    real(IDS_real) :: tmpVx( mpg%nVx )
    integer :: nSubsets  !< number of grid subsets to fill
    integer :: iSubset   !< Grid subset iterator
    integer :: iSubsetID !< Grid subset identifier index
    integer :: ggdID     !< Grid identifier index
    integer :: ndim      !< Grid subset dimension
    integer :: i         !< Iterator
    external xerrab

    ggdId = basegrid%identifier%index
#if ( IMAS_MINOR_VERSION < 15 && IMAS_MAJOR_VERSION < 4 )
    !! Assign 5+4 grid subsets
    nSubsets = 9
#else
    nSubsets = size(basegrid%grid_subset)
    if (nSubsets.eq.0) return
#endif
    !! Interpolate data to vertices
    call intvertex( mpg%nCv, mpg%nVx, mpg, geo%vxVol, value, tmpVx )

    !! Interpolate data to cell faces, using a volume weighting
    tmpFace = 0.0_IDS_real
    call intface( mpg%nCv, mpg%nFc, mpg%fcCv, geo%fcVol, value, tmpFace)

    !! Allocate data fields for grid subsets
    if (.not.associated( val ) ) then
      allocate( val(nSubsets) )
    end if

    do iSubset = 1, nSubsets
#if ( IMAS_MINOR_VERSION < 15 && IMAS_MAJOR_VERSION < 4 )
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
      ndim = basegrid%grid_subset(iSubset)%dimension
      iSubsetID = basegrid%grid_subset(iSubset)%identifier%index
#endif
      if (ndim.eq.IDS_INT_INVALID) then
        select case (iSubsetID)
        case( GRID_SUBSET_NODES, GRID_SUBSET_X_POINTS, &
            & GRID_SUBSET_MAGNETIC_AXIS,               &
            & GRID_SUBSET_INNER_MIDPLANE_SEPARATRIX,   &
            & GRID_SUBSET_OUTER_MIDPLANE_SEPARATRIX,   &
            & GRID_SUBSET_INNER_STRIKEPOINT,           &
            & GRID_SUBSET_OUTER_STRIKEPOINT,           &
            & GRID_SUBSET_INNER_STRIKEPOINT_INACTIVE,  &
            & GRID_SUBSET_OUTER_STRIKEPOINT_INACTIVE )
          ndim = 1
        case( GRID_SUBSET_EDGES, &
            & GRID_SUBSET_X_ALIGNED_EDGES, GRID_SUBSET_Y_ALIGNED_EDGES, &
            & GRID_SUBSET_CORE_BOUNDARY, GRID_SUBSET_SEPARATRIX, &
            & GRID_SUBSET_ACTIVE_SEPARATRIX, GRID_SUBSET_MAIN_CHAMBER_WALL, &
            & GRID_SUBSET_OUTER_BAFFLE, GRID_SUBSET_INNER_BAFFLE, &
            & GRID_SUBSET_OUTER_PFR_WALL, GRID_SUBSET_INNER_PFR_WALL, &
            & GRID_SUBSET_MAIN_WALL, GRID_SUBSET_PFR_WALL, &
            & GRID_SUBSET_FULL_WALL, &
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
            & GRID_SUBSET_INNER_TARGET_INACTIVE, &
            & GRID_SUBSET_OUTER_SF_LEG_ENTRANCE_1, &
            & GRID_SUBSET_OUTER_SF_LEG_ENTRANCE_2, &
            & GRID_SUBSET_OUTER_SF_PFR_CONNECTION_1, &
            & GRID_SUBSET_OUTER_SF_PFR_CONNECTION_2)
          ndim = 2
        case( GRID_SUBSET_CELLS, GRID_SUBSET_BETWEEN_SEPARATRICES, &
            & GRID_SUBSET_CORE, GRID_SUBSET_SOL, &
            & GRID_SUBSET_INNER_MIDPLANE, GRID_SUBSET_OUTER_MIDPLANE, &
            & GRID_SUBSET_OUTER_DIVERTOR, GRID_SUBSET_INNER_DIVERTOR, &
            & GRID_SUBSET_OUTER_DIVERTOR_INACTIVE, &
            & GRID_SUBSET_INNER_DIVERTOR_INACTIVE )
          ndim = 3
        case( GRID_SUBSET_VOLUMES )
          ndim = 4
        end select
      end if
      select case (ndim)
      case ( 1 ) !< Grid subset consists of nodes
        idsdata => b2_IMAS_Transform_Data_B2_To_IDS_Vertex(        &
                     &   basegrid, iSubset, mpg, tmpVx )
        if ( size( idsdata ) > 0 ) then
#if ( GGD_MINOR_VERSION > 8 || GGD_MAJOR_VERSION > 1 )
          call gridWriteData( val( iSubset ), ggdID, iSubsetID, idsdata )
#else
          val(iSubset)%grid_index = ggdId
          call gridWriteData( val( iSubset ), iSubsetID, idsdata )
#endif
        else
          val(iSubset)%grid_index = ggdId
          val(iSubset)%grid_subset_index = iSubsetID
        end if
        deallocate( idsdata )
      case ( 2 ) !< Grid subset consists of faces
        idsdata => b2_IMAS_Transform_Data_B2_To_IDS_Face(          &
                     &   basegrid, iSubset, mpg, tmpFace )
        if ( size( idsdata ) > 0 ) then
#if ( GGD_MINOR_VERSION > 8 || GGD_MAJOR_VERSION > 1 )
          call gridWriteData( val( iSubset ), ggdID, iSubsetID, idsdata )
#else
          val(iSubset)%grid_index = ggdId
          call gridWriteData( val( iSubset ), iSubsetID, idsdata )
#endif
        else
          val(iSubset)%grid_index = ggdId
          val(iSubset)%grid_subset_index = iSubsetID
        end if
        deallocate( idsdata )
      case ( 3 ) !< Grid subset consists of cells
        idsdata => b2_IMAS_Transform_Data_B2_To_IDS_Cell(          &
                      &  basegrid, iSubset, mpg, value )
        if ( size( idsdata ) > 0 ) then
#if ( GGD_MINOR_VERSION > 8 || GGD_MAJOR_VERSION > 1 )
          call gridWriteData( val( iSubset ), ggdID, iSubsetID, idsdata )
#else
          val(iSubset)%grid_index = ggdId
          call gridWriteData( val( iSubset ), iSubsetID, idsdata )
#endif
        else
          val(iSubset)%grid_index = ggdId
          val(iSubset)%grid_subset_index = iSubsetID
        end if
        deallocate( idsdata )
      case ( 4 ) !< Grid subset consists of volumes
        call xerrab( 'Case for volumetric grid subsets not yet coded !')
      case default
        call xerrab( 'Unknown grid subset '//int2str(iSubset)// &
            &        ' dimension : '//int2str(ndim) )
      end select
    end do

    return
    end subroutine write_IDS_quantity

    !> Write a scalar B2 cell quantity to ids_generic_grid_scalar
    subroutine write_cell_scalar( basegrid, mpg, scalar, b2CellData )
    implicit none
#if ( IMAS_MINOR_VERSION < 15 && IMAS_MAJOR_VERSION < 4 )
    type(ids_generic_grid_dynamic), intent(in) :: basegrid !< Type of IDS
        !< data structure, designed for handling grid geometry data
#else
    type(ids_generic_grid_aos3_root), intent(in) :: basegrid !< Type of IDS
        !< data structure, designed for handling grid geometry data
#endif
    type(mapping), intent(in) :: mpg
    type(ids_generic_grid_scalar), intent(inout), pointer :: scalar(:)
        !< Type of IDS data structure, designed for scalar data handling
    real(IDS_real), intent(in) :: b2CellData( mpg%nCv )
    real(IDS_real), dimension(:), pointer :: idsdata    !< Array for
        !< handling data field values
    integer :: nSubsets  !< number of grid subsets to fill
    integer :: iSubset   !< Grid subset iterator
    integer :: ndim      !< Grid subset dimension
    integer :: iSubsetID !< Grid subset identifier index
    integer :: ggdID     !< Grid identifier index

    ggdId = basegrid%identifier%index
#if ( IMAS_MINOR_VERSION < 15 && IMAS_MAJOR_VERSION < 4 )
    nSubsets = 1
#else
    nSubsets = size(basegrid%grid_subset)
    if (nSubsets.eq.0) return
#endif
    !! Allocate data fields for grid subsets
    if (.not.associated( scalar ) ) then
      allocate( scalar(nSubsets) )
    end if

    do iSubset = 1, nSubsets
#if ( IMAS_MINOR_VERSION < 15 && IMAS_MAJOR_VERSION < 4 )
      ndim = 3
      iSubsetID = GRID_SUBSET_CELLS
#else
      ndim = basegrid%grid_subset(iSubset)%dimension
      iSubsetID = basegrid%grid_subset(iSubset)%identifier%index
#endif
      if (ndim.eq.IDS_INT_INVALID) then
        select case (iSubsetID)
        case( GRID_SUBSET_NODES, GRID_SUBSET_X_POINTS, &
            & GRID_SUBSET_MAGNETIC_AXIS,               &
            & GRID_SUBSET_INNER_MIDPLANE_SEPARATRIX,   &
            & GRID_SUBSET_OUTER_MIDPLANE_SEPARATRIX,   &
            & GRID_SUBSET_INNER_STRIKEPOINT,           &
            & GRID_SUBSET_OUTER_STRIKEPOINT,           &
            & GRID_SUBSET_INNER_STRIKEPOINT_INACTIVE,  &
            & GRID_SUBSET_OUTER_STRIKEPOINT_INACTIVE )
          ndim = 1
        case( GRID_SUBSET_EDGES, &
            & GRID_SUBSET_X_ALIGNED_EDGES, GRID_SUBSET_Y_ALIGNED_EDGES, &
            & GRID_SUBSET_CORE_BOUNDARY, GRID_SUBSET_SEPARATRIX, &
            & GRID_SUBSET_ACTIVE_SEPARATRIX, GRID_SUBSET_MAIN_CHAMBER_WALL, &
            & GRID_SUBSET_OUTER_BAFFLE, GRID_SUBSET_INNER_BAFFLE, &
            & GRID_SUBSET_OUTER_PFR_WALL, GRID_SUBSET_INNER_PFR_WALL, &
            & GRID_SUBSET_MAIN_WALL, GRID_SUBSET_PFR_WALL, &
            & GRID_SUBSET_FULL_WALL, &
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
            & GRID_SUBSET_INNER_TARGET_INACTIVE, &
            & GRID_SUBSET_OUTER_SF_LEG_ENTRANCE_1, &
            & GRID_SUBSET_OUTER_SF_LEG_ENTRANCE_2, &
            & GRID_SUBSET_OUTER_SF_PFR_CONNECTION_1, &
            & GRID_SUBSET_OUTER_SF_PFR_CONNECTION_2)
          ndim = 2
        case( GRID_SUBSET_CELLS, GRID_SUBSET_BETWEEN_SEPARATRICES, &
            & GRID_SUBSET_CORE, GRID_SUBSET_SOL, &
            & GRID_SUBSET_INNER_MIDPLANE, GRID_SUBSET_OUTER_MIDPLANE, &
            & GRID_SUBSET_OUTER_DIVERTOR, GRID_SUBSET_INNER_DIVERTOR, &
            & GRID_SUBSET_OUTER_DIVERTOR_INACTIVE, &
            & GRID_SUBSET_INNER_DIVERTOR_INACTIVE )
          ndim = 3
        case( GRID_SUBSET_VOLUMES )
          ndim = 4
        end select
      end if
      if (ndim.ne.3) cycle

      idsdata => b2_IMAS_Transform_Data_B2_To_IDS_Cell(  &
          &   basegrid, iSubset, mpg, b2CellData )
      if ( size( idsdata ) > 0 ) then
#if ( GGD_MINOR_VERSION > 8 || GGD_MAJOR_VERSION > 1 )
        call gridWriteData( scalar( iSubset ), ggdID, iSubsetID, idsdata )
#else
        scalar(iSubset)%grid_index = ggdId
        call gridWriteData( scalar( iSubset ), iSubsetID, idsdata )
#endif
      else
        scalar(iSubset)%grid_index = ggdId
        scalar(iSubset)%grid_subset_index = iSubsetID
      end if
      deallocate(idsdata)
    end do

    return
    end subroutine write_cell_scalar

    !> Write a vector component B2 cell quantity to ids_generic_grid_vector
    !! components
    !! @note Available IDS vector component data fields (vector IDs):
    !!          - VEC_ALIGN_RADIAL_ID ( "Radial" ),
    !!          - VEC_ALIGN_DIAMAGNETIC_ID ( "Diamagnetic" ),
    !!          - VEC_ALIGN_PARALLEL_ID ( "Parallel" ),
    !!          - VEC_ALIGN_POLOIDAL_ID ( "Poloidal" ),
    !!          - VEC_ALIGN_TOROIDAL_ID ( "Toroidal" ),
    !!          - VEC_ALIGN_PHI_ID ( "Phi" ),
    !!          - VEC_ALIGN_R_MAJOR_ID ( "R" ),
    !!          - VEC_ALIGN_Z_ID ( "Z" )
    subroutine write_cell_vector_component( basegrid, mpg, &
       &  vectorComponent, b2CellData, vectorID )
    implicit none
#if ( IMAS_MINOR_VERSION < 15 && IMAS_MAJOR_VERSION < 4 )
    type(ids_generic_grid_dynamic), intent(in) :: basegrid !< Type of IDS
        !< data structure, designed for handling grid geometry data
#else
    type(ids_generic_grid_aos3_root), intent(in) :: basegrid !< Type of IDS
        !< data structure, designed for handling grid geometry data
#endif
    type(mapping), intent(in) :: mpg
    type(ids_generic_grid_vector_components), intent(inout),    &
              &   pointer :: vectorComponent(:) !< Type of IDS data structure,
                    !> designed for vector data handling
    real(IDS_real), intent(in) :: b2CellData( mpg%nCv )
    real(IDS_real), dimension(:), pointer :: idsdata    !< Array for
                    !< handling data field values
    character(len=*), intent(in) :: vectorID    !< Vector ID (e.g.
                                                !< VEC_ALIGN_RADIAL_ID)
    integer :: nSubsets  !< number of grid subsets to fill
    integer :: iSubset   !< Grid subset iterator
    integer :: iSubsetID !< Grid subset identifier index
    integer :: ndim      !< Grid subset dimension
    integer :: ggdID     !< Grid identifier index

    ggdId = basegrid%identifier%index
#if ( IMAS_MINOR_VERSION < 15 && IMAS_MAJOR_VERSION < 4 )
    nSubsets = 1
#else
    nSubsets = size(basegrid%grid_subset)
    if (nSubsets.eq.0) return
#endif
    !! If required, allocate storage
    if ( .not. associated( vectorComponent ) ) then
      allocate( vectorComponent(nSubsets) )
    end if

    do iSubset = 1, nSubsets
#if ( IMAS_MINOR_VERSION < 15 && IMAS_MAJOR_VERSION < 4 )
      ndim = 3
      iSubsetID = GRID_SUBSET_CELLS
#else
      ndim = basegrid%grid_subset(iSubset)%dimension
      iSubsetID = basegrid%grid_subset(iSubset)%identifier%index
#endif
      if (ndim.eq.IDS_INT_INVALID) then
        select case (iSubsetID)
        case( GRID_SUBSET_NODES, GRID_SUBSET_X_POINTS, &
            & GRID_SUBSET_MAGNETIC_AXIS,               &
            & GRID_SUBSET_INNER_MIDPLANE_SEPARATRIX,   &
            & GRID_SUBSET_OUTER_MIDPLANE_SEPARATRIX,   &
            & GRID_SUBSET_INNER_STRIKEPOINT,           &
            & GRID_SUBSET_OUTER_STRIKEPOINT,           &
            & GRID_SUBSET_INNER_STRIKEPOINT_INACTIVE,  &
            & GRID_SUBSET_OUTER_STRIKEPOINT_INACTIVE )
          ndim = 1
        case( GRID_SUBSET_EDGES, &
            & GRID_SUBSET_X_ALIGNED_EDGES, GRID_SUBSET_Y_ALIGNED_EDGES, &
            & GRID_SUBSET_CORE_BOUNDARY, GRID_SUBSET_SEPARATRIX, &
            & GRID_SUBSET_ACTIVE_SEPARATRIX, GRID_SUBSET_MAIN_CHAMBER_WALL, &
            & GRID_SUBSET_OUTER_BAFFLE, GRID_SUBSET_INNER_BAFFLE, &
            & GRID_SUBSET_OUTER_PFR_WALL, GRID_SUBSET_INNER_PFR_WALL, &
            & GRID_SUBSET_MAIN_WALL, GRID_SUBSET_PFR_WALL, &
            & GRID_SUBSET_FULL_WALL, &
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
            & GRID_SUBSET_INNER_TARGET_INACTIVE, &
            & GRID_SUBSET_OUTER_SF_LEG_ENTRANCE_1, &
            & GRID_SUBSET_OUTER_SF_LEG_ENTRANCE_2, &
            & GRID_SUBSET_OUTER_SF_PFR_CONNECTION_1, &
            & GRID_SUBSET_OUTER_SF_PFR_CONNECTION_2 )
          ndim = 2
        case( GRID_SUBSET_CELLS, GRID_SUBSET_BETWEEN_SEPARATRICES, &
            & GRID_SUBSET_CORE, GRID_SUBSET_SOL, &
            & GRID_SUBSET_INNER_MIDPLANE, GRID_SUBSET_OUTER_MIDPLANE, &
            & GRID_SUBSET_OUTER_DIVERTOR, GRID_SUBSET_INNER_DIVERTOR, &
            & GRID_SUBSET_OUTER_DIVERTOR_INACTIVE, &
            & GRID_SUBSET_INNER_DIVERTOR_INACTIVE )
          ndim = 3
        case( GRID_SUBSET_VOLUMES )
          ndim = 4
        end select
      end if
      if (ndim.ne.3) cycle
      idsdata => b2_IMAS_Transform_Data_B2_To_IDS_Cell(       &
                   &   basegrid, iSubset, mpg, b2CellData )
      if ( size(idsdata) > 0 ) then
        call B2grid_Write_Data_Vector_Components( vectorComponent(iSubset), &
            &   ggdID, iSubsetID, vectorID, idsdata )
      else
        vectorComponent(iSubset)%grid_index = ggdId
        vectorComponent(iSubset)%grid_subset_index = iSubsetID
      end if
      deallocate(idsdata)
    end do

    return
    end subroutine write_cell_vector_component

    !> Write a scalar B2 face quantity to ids_generic_grid_scalar
    subroutine write_face_scalar( basegrid, mpg, val, value )
    implicit none
#if ( IMAS_MINOR_VERSION < 15 && IMAS_MAJOR_VERSION < 4 )
    type (ids_generic_grid_dynamic), intent(in) :: basegrid !< Type of IDS
        !< data structure, designed for handling grid geometry data
#else
    type (ids_generic_grid_aos3_root), intent(in) :: basegrid !< Type of IDS
        !< data structure, designed for handling grid geometry data
#endif
    type (mapping), intent(in) :: mpg
    type (ids_generic_grid_scalar), pointer, intent(inout) :: val(:)
        !< Type of IDS data structure, designed for scalar data handling
        !< (in this case scalars residing on grid faces)
    real(IDS_real), intent(in) :: value( mpg%nFc )
    integer :: nSubsets  !< number of grid subsets to fill
    integer :: iSubset   !< Grid subset iterator
    integer :: iSubsetID !< Grid subset identifier index
    integer :: ggdID     !< Grid identifier index
    integer :: ndim      !< Grid subset dimension

    ggdId = basegrid%identifier%index
#if ( IMAS_MINOR_VERSION < 15 && IMAS_MAJOR_VERSION < 4 )
    nSubsets = 2
#else
    nSubsets = size(basegrid%grid_subset)
    if (nSubsets.eq.0) return
#endif
    !! Allocate data fields for grid subsets
    if (.not.associated( val ) ) then
      allocate( val(nSubsets) )
    end if

    do iSubset = 1, nSubsets
#if ( IMAS_MINOR_VERSION < 15 && IMAS_MAJOR_VERSION < 4 )
      select case (iSubset)
      case (1)
        iSubsetID = GRID_SUBSET_X_ALIGNED_EDGES
        ndim = 2
      case (2)
        iSubsetID = GRID_SUBSET_Y_ALIGNED_EDGES
        ndim = 2
      case default
        iSubsetID = iSubset
        ndim = IDS_INT_INVALID
       end select
#else
       ndim = basegrid%grid_subset(iSubset)%dimension
       iSubsetID = basegrid%grid_subset(iSubset)%identifier%index
#endif
       if (ndim.eq.IDS_INT_INVALID) then
         select case (iSubsetID)
         case( GRID_SUBSET_NODES, GRID_SUBSET_X_POINTS, &
             & GRID_SUBSET_MAGNETIC_AXIS,               &
             & GRID_SUBSET_INNER_MIDPLANE_SEPARATRIX,   &
             & GRID_SUBSET_OUTER_MIDPLANE_SEPARATRIX,   &
             & GRID_SUBSET_INNER_STRIKEPOINT,           &
             & GRID_SUBSET_OUTER_STRIKEPOINT,           &
             & GRID_SUBSET_INNER_STRIKEPOINT_INACTIVE,  &
             & GRID_SUBSET_OUTER_STRIKEPOINT_INACTIVE )
           ndim = 1
         case( GRID_SUBSET_EDGES, &
             & GRID_SUBSET_X_ALIGNED_EDGES, GRID_SUBSET_Y_ALIGNED_EDGES, &
             & GRID_SUBSET_CORE_BOUNDARY, GRID_SUBSET_SEPARATRIX, &
             & GRID_SUBSET_ACTIVE_SEPARATRIX, GRID_SUBSET_MAIN_CHAMBER_WALL, &
             & GRID_SUBSET_OUTER_BAFFLE, GRID_SUBSET_INNER_BAFFLE, &
             & GRID_SUBSET_OUTER_PFR_WALL, GRID_SUBSET_INNER_PFR_WALL, &
             & GRID_SUBSET_MAIN_WALL, GRID_SUBSET_PFR_WALL, &
             & GRID_SUBSET_FULL_WALL, &
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
             & GRID_SUBSET_INNER_TARGET_INACTIVE, &
             & GRID_SUBSET_OUTER_SF_LEG_ENTRANCE_1, &
             & GRID_SUBSET_OUTER_SF_LEG_ENTRANCE_2, &
             & GRID_SUBSET_OUTER_SF_PFR_CONNECTION_1, &
             & GRID_SUBSET_OUTER_SF_PFR_CONNECTION_2)
           ndim = 2
         case( GRID_SUBSET_CELLS, GRID_SUBSET_BETWEEN_SEPARATRICES, &
             & GRID_SUBSET_CORE, GRID_SUBSET_SOL, &
             & GRID_SUBSET_INNER_MIDPLANE, GRID_SUBSET_OUTER_MIDPLANE, &
             & GRID_SUBSET_OUTER_DIVERTOR, GRID_SUBSET_INNER_DIVERTOR, &
             & GRID_SUBSET_OUTER_DIVERTOR_INACTIVE, &
             & GRID_SUBSET_INNER_DIVERTOR_INACTIVE )
           ndim = 3
         case( GRID_SUBSET_VOLUMES )
           ndim = 4
         end select
       end if
       if (ndim.ne.2) cycle
       call write_face_vector( basegrid, mpg, val( iSubset ), value, &
           &    ggdID, iSubsetID, iSubset )
    end do

    return
    end subroutine write_face_scalar

    !> Write a vector component B2 face quantity to ids_generic_grid_vector
    !! components
    !! @note Available IDS vector component data fields (vector IDs):
    !!          - VEC_ALIGN_RADIAL_ID ( "Radial" ),
    !!          - VEC_ALIGN_DIAMAGNETIC_ID ( "Diamagnetic" ),
    !!          - VEC_ALIGN_PARALLEL_ID ( "Parallel" ),
    !!          - VEC_ALIGN_POLOIDAL_ID ( "Poloidal" ),
    !!          - VEC_ALIGN_TOROIDAL_ID ( "Toroidal" ),
    !!          - VEC_ALIGN_PHI_ID ( "Phi" ),
    !!          - VEC_ALIGN_R_MAJOR_ID ( "R" ),
    !!          - VEC_ALIGN_Z_ID ( "Z" )
    subroutine write_face_vector_component( basegrid, mpg, &
       &   vectorComponent, b2FaceData, vectorID )
    implicit none
#if ( IMAS_MINOR_VERSION < 15 && IMAS_MAJOR_VERSION < 4 )
    type(ids_generic_grid_dynamic), intent(in) :: basegrid !< Type of IDS
        !< data structure, designed for handling grid geometry data
#else
    type(ids_generic_grid_aos3_root), intent(in) :: basegrid !< Type of IDS
        !< data structure, designed for handling grid geometry data
#endif
    type(mapping), intent(in) :: mpg
    type(ids_generic_grid_vector_components), intent(inout),    &
       &   pointer :: vectorComponent(:) !< Type of IDS data structure,
                                         !> designed for vector data handling
    real(IDS_real), intent(in) :: b2FaceData( mpg%nFc )
    real(IDS_real), dimension(:), pointer :: idsdata    !< Array for
                                         !< handling data field values
    character(len=*), intent(in) :: vectorID    !< Vector ID (e.g.
                                                !< VEC_ALIGN_RADIAL_ID)
    integer :: nSubsets  !< number of grid subsets to fill
    integer :: iSubset   !< Grid subset iterator
    integer :: iSubsetID !< Grid subset identifier index
    integer :: ndim      !< Grid subset dimension
    integer :: ggdID     !< Grid identifier index

    ggdId = basegrid%identifier%index
#if ( IMAS_MINOR_VERSION < 15 && IMAS_MAJOR_VERSION < 4 )
    nSubsets = 1
#else
    nSubsets = size(basegrid%grid_subset)
    if (nSubsets.eq.0) return
#endif
    !! If required, allocate storage
    if ( .not. associated( vectorComponent ) ) then
      allocate( vectorComponent(nSubsets) )
    end if

    do iSubset = 1, nSubsets
#if ( IMAS_MINOR_VERSION < 15 && IMAS_MAJOR_VERSION < 4 )
      ndim = 2
      iSubsetID = GRID_SUBSET_FACES
#else
      ndim = basegrid%grid_subset(iSubset)%dimension
      iSubsetID = basegrid%grid_subset(iSubset)%identifier%index
#endif
      if (ndim.eq.IDS_INT_INVALID) then
        select case (iSubsetID)
        case( GRID_SUBSET_NODES, GRID_SUBSET_X_POINTS, &
            & GRID_SUBSET_MAGNETIC_AXIS,               &
            & GRID_SUBSET_INNER_MIDPLANE_SEPARATRIX,   &
            & GRID_SUBSET_OUTER_MIDPLANE_SEPARATRIX,   &
            & GRID_SUBSET_INNER_STRIKEPOINT,           &
            & GRID_SUBSET_OUTER_STRIKEPOINT,           &
            & GRID_SUBSET_INNER_STRIKEPOINT_INACTIVE,  &
            & GRID_SUBSET_OUTER_STRIKEPOINT_INACTIVE )
          ndim = 1
        case( GRID_SUBSET_EDGES, &
            & GRID_SUBSET_X_ALIGNED_EDGES, GRID_SUBSET_Y_ALIGNED_EDGES, &
            & GRID_SUBSET_CORE_BOUNDARY, GRID_SUBSET_SEPARATRIX, &
            & GRID_SUBSET_ACTIVE_SEPARATRIX, GRID_SUBSET_MAIN_CHAMBER_WALL, &
            & GRID_SUBSET_OUTER_BAFFLE, GRID_SUBSET_INNER_BAFFLE, &
            & GRID_SUBSET_OUTER_PFR_WALL, GRID_SUBSET_INNER_PFR_WALL, &
            & GRID_SUBSET_MAIN_WALL, GRID_SUBSET_PFR_WALL, &
            & GRID_SUBSET_FULL_WALL, &
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
            & GRID_SUBSET_INNER_TARGET_INACTIVE, &
            & GRID_SUBSET_OUTER_SF_LEG_ENTRANCE_1, &
            & GRID_SUBSET_OUTER_SF_LEG_ENTRANCE_2, &
            & GRID_SUBSET_OUTER_SF_PFR_CONNECTION_1, &
            & GRID_SUBSET_OUTER_SF_PFR_CONNECTION_2 )
          ndim = 2
        case( GRID_SUBSET_CELLS, GRID_SUBSET_BETWEEN_SEPARATRICES, &
            & GRID_SUBSET_CORE, GRID_SUBSET_SOL, &
            & GRID_SUBSET_INNER_MIDPLANE, GRID_SUBSET_OUTER_MIDPLANE, &
            & GRID_SUBSET_OUTER_DIVERTOR, GRID_SUBSET_INNER_DIVERTOR, &
            & GRID_SUBSET_OUTER_DIVERTOR_INACTIVE, &
            & GRID_SUBSET_INNER_DIVERTOR_INACTIVE )
          ndim = 3
        case( GRID_SUBSET_VOLUMES )
          ndim = 4
        end select
      end if
      if (ndim.ne.2) cycle
      idsdata => b2_IMAS_Transform_Data_B2_To_IDS_Face( &
                   &   basegrid, iSubset, mpg, b2FaceData )
      if ( size(idsdata) > 0 ) then
        call B2grid_Write_Data_Vector_Components( vectorComponent(iSubset), &
                   &   ggdID, iSubsetID, vectorID, idsdata )
      else
        vectorComponent(iSubset)%grid_index = ggdId
        vectorComponent(iSubset)%grid_subset_index = iSubsetID
      end if
      deallocate(idsdata)
    end do

    return
    end subroutine write_face_vector_component

    !> Write a scalar B2 vertex quantity to ids_generic_grid_scalar
    subroutine write_vertex_scalar( basegrid, mpg, scalar, b2VertexData )
    implicit none
#if ( IMAS_MINOR_VERSION < 15 && IMAS_MAJOR_VERSION < 4 )
    type(ids_generic_grid_dynamic), intent(in) :: basegrid !< Type of IDS
        !< data structure, designed for handling grid geometry data
#else
    type(ids_generic_grid_aos3_root), intent(in) :: basegrid !< Type of IDS
        !< data structure, designed for handling grid geometry data
#endif
    type(mapping), intent(in) :: mpg
    type(ids_generic_grid_scalar), intent(inout), pointer :: scalar(:)
        !< Type of IDS data structure, designed for scalar data handling
    real(IDS_real), intent(in) :: b2VertexData( mpg%nVx )
    real(IDS_real), dimension(:), pointer :: idsdata    !< Array for
        !< handling data field values
    integer :: nSubsets  !< number of grid subsets to fill
    integer :: iSubset   !< Grid subset iterator
    integer :: ndim      !< Grid subset dimension
    integer :: iSubsetID !< Grid subset identifier index
    integer :: ggdID     !< Grid identifier index

    ggdId = basegrid%identifier%index
#if ( IMAS_MINOR_VERSION < 15 && IMAS_MAJOR_VERSION < 4 )
    nSubsets = 1
#else
    nSubsets = size(basegrid%grid_subset)
    if (nSubsets.eq.0) return
#endif
    !! Allocate data fields for grid subsets
    if (.not.associated( scalar ) ) then
      allocate( scalar(nSubsets) )
    end if

    do iSubset = 1, nSubsets
#if ( IMAS_MINOR_VERSION < 15 && IMAS_MAJOR_VERSION < 4 )
      ndim = 3
      iSubsetID = GRID_SUBSET_CELLS
#else
      ndim = basegrid%grid_subset(iSubset)%dimension
      iSubsetID = basegrid%grid_subset(iSubset)%identifier%index
#endif
      if (ndim.eq.IDS_INT_INVALID) then
        select case (iSubsetID)
        case( GRID_SUBSET_NODES, GRID_SUBSET_X_POINTS, &
            & GRID_SUBSET_MAGNETIC_AXIS,               &
            & GRID_SUBSET_INNER_MIDPLANE_SEPARATRIX,   &
            & GRID_SUBSET_OUTER_MIDPLANE_SEPARATRIX,   &
            & GRID_SUBSET_INNER_STRIKEPOINT,           &
            & GRID_SUBSET_OUTER_STRIKEPOINT,           &
            & GRID_SUBSET_INNER_STRIKEPOINT_INACTIVE,  &
            & GRID_SUBSET_OUTER_STRIKEPOINT_INACTIVE )
          ndim = 1
        case( GRID_SUBSET_EDGES, &
            & GRID_SUBSET_X_ALIGNED_EDGES, GRID_SUBSET_Y_ALIGNED_EDGES, &
            & GRID_SUBSET_CORE_BOUNDARY, GRID_SUBSET_SEPARATRIX, &
            & GRID_SUBSET_ACTIVE_SEPARATRIX, GRID_SUBSET_MAIN_CHAMBER_WALL, &
            & GRID_SUBSET_OUTER_BAFFLE, GRID_SUBSET_INNER_BAFFLE, &
            & GRID_SUBSET_OUTER_PFR_WALL, GRID_SUBSET_INNER_PFR_WALL, &
            & GRID_SUBSET_MAIN_WALL, GRID_SUBSET_PFR_WALL, &
            & GRID_SUBSET_FULL_WALL, &
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
            & GRID_SUBSET_INNER_TARGET_INACTIVE, &
            & GRID_SUBSET_OUTER_SF_LEG_ENTRANCE_1, &
            & GRID_SUBSET_OUTER_SF_LEG_ENTRANCE_2, &
            & GRID_SUBSET_OUTER_SF_PFR_CONNECTION_1, &
            & GRID_SUBSET_OUTER_SF_PFR_CONNECTION_2)
          ndim = 2
        case( GRID_SUBSET_CELLS, GRID_SUBSET_BETWEEN_SEPARATRICES, &
            & GRID_SUBSET_CORE, GRID_SUBSET_SOL, &
            & GRID_SUBSET_INNER_MIDPLANE, GRID_SUBSET_OUTER_MIDPLANE, &
            & GRID_SUBSET_OUTER_DIVERTOR, GRID_SUBSET_INNER_DIVERTOR, &
            & GRID_SUBSET_OUTER_DIVERTOR_INACTIVE, &
            & GRID_SUBSET_INNER_DIVERTOR_INACTIVE )
          ndim = 3
        case( GRID_SUBSET_VOLUMES )
          ndim = 4
        end select
      end if
      if (ndim.ne.1) cycle

      idsdata => b2_IMAS_Transform_Data_B2_To_IDS_Vertex(  &
          &   basegrid, iSubset, mpg, b2VertexData )
      if ( size( idsdata ) > 0 ) then
#if ( GGD_MINOR_VERSION > 8 || GGD_MAJOR_VERSION > 1 )
        call gridWriteData( scalar( iSubset ), ggdID, iSubsetID, idsdata )
#else
        scalar(iSubset)%grid_index = ggdId
        call gridWriteData( scalar( iSubset ), iSubsetID, idsdata )
#endif
      else
        scalar(iSubset)%grid_index = ggdId
        scalar(iSubset)%grid_subset_index = iSubsetID
      end if
      deallocate(idsdata)
    end do

    return
    end subroutine write_vertex_scalar

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
    subroutine write_face_vector( basegrid, mpg, vector, b2FaceData, &
       &   gridID, gridSubsetID, gridSubsetInd )
    implicit none
#if ( IMAS_MINOR_VERSION < 15 && IMAS_MAJOR_VERSION < 4 )
    type(ids_generic_grid_dynamic), intent(in) :: basegrid !< Type of IDS
        !< data structure, designed for handling grid geometry data
#else
    type(ids_generic_grid_aos3_root), intent(in) :: basegrid !< Type of IDS
        !< data structure, designed for handling grid geometry data
#endif
    type(mapping), intent(in) :: mpg
    type(ids_generic_grid_scalar), intent(inout) :: vector
        !< Type of IDS data structure, designed for scalar data handling
        !< (in this case 1D vector)
    real(IDS_real), intent(in) :: b2FaceData( mpg%nFc, 0:1 )
    integer, intent(in) :: gridID                    !< Grid identifier index
    integer, intent(in), optional :: gridSubsetID    !< Grid subset identifier index
    integer, intent(in), optional :: gridSubsetInd   !< Base grid subset index
    real(IDS_real), dimension(:), pointer :: idsdata !< Dummy array
        !< for holding data field values

    if ( .not. present(gridSubsetInd) ) then
      !! Fill in vector component data
      idsdata => b2_IMAS_Transform_Data_B2_To_IDS_Face(    &
               &   basegrid, GRID_SUBSET_Y_ALIGNED_EDGES, mpg, b2FaceData)
      if ( size( idsdata ) > 0 ) then
#if ( GGD_MINOR_VERSION > 8 || GGD_MAJOR_VERSION > 1 )
        call gridWriteData( vector, gridId, GRID_SUBSET_Y_ALIGNED_EDGES, idsdata )
#else
        vector%grid_index = gridId
        call gridWriteData( vector, GRID_SUBSET_Y_ALIGNED_EDGES, idsdata )
#endif
      else
        vector%grid_index = gridId
        vector%grid_subset_index = GRID_SUBSET_Y_ALIGNED_EDGES
      end if
      deallocate(idsdata)
      idsdata => b2_IMAS_Transform_Data_B2_To_IDS_Face(    &
               &   basegrid, GRID_SUBSET_X_ALIGNED_EDGES, mpg, b2FaceData)
      if ( size( idsdata ) > 0 ) then
#if ( GGD_MINOR_VERSION > 8 || GGD_MAJOR_VERSION > 1 )
        call gridWriteData( vector, gridId, GRID_SUBSET_X_ALIGNED_EDGES, idsdata )
#else
        vector%grid_index = gridId
        call gridWriteData( vector, GRID_SUBSET_X_ALIGNED_EDGES, idsdata )
#endif
      else
        vector%grid_index = gridId
        vector%grid_subset_index = GRID_SUBSET_X_ALIGNED_EDGES
      end if
      deallocate(idsdata)
    else
      idsdata => b2_IMAS_Transform_Data_B2_To_IDS_Face(    &
               &   basegrid, gridSubsetInd, mpg, b2FaceData)
      if ( size( idsdata ) > 0 ) then
#if ( GGD_MINOR_VERSION > 8 || GGD_MAJOR_VERSION > 1 )
        call gridWriteData( vector, gridId, gridSubsetID, idsdata )
#else
        vector%grid_index = gridId
        call gridWriteData( vector, gridSubsetID, idsdata )
#endif
      else
        vector%grid_index = gridId
        vector%grid_subset_index = gridSubsetID
      end if
      deallocate(idsdata)
    end if

    return
    end subroutine write_face_vector
#endif

    subroutine write_sourced_constant( val, value )
    implicit none
    type(ids_summary_constant_flt_0d) :: val
        !< Type of IDS data structure, designed for sourced real constant data handling
    real(IDS_real), intent(in) :: value

    val%value = value
    allocate( val%source(1) )
    val%source = source

    return
    end subroutine write_sourced_constant

    subroutine write_sourced_int_constant( ival, ivalue )
    implicit none
    type(ids_summary_constant_int_0d) :: ival
        !< Type of IDS data structure, designed for sourced integer constant data handling
    integer, intent(in) :: ivalue

    ival%value = ivalue
    allocate( ival%source(1) )
    ival%source = source

    return
    end subroutine write_sourced_int_constant

    subroutine write_sourced_integer( ival, ivalue )
    implicit none
    type(ids_summary_dynamic_int_1d_root) :: ival
        !< Type of IDS data structure, designed for sourced integer data handling
    integer, intent(in) :: ivalue

    allocate( ival%value( num_slices ) )
    ival%value( slice_index ) = ivalue
    allocate( ival%source(1) )
    ival%source = source

    return
    end subroutine write_sourced_integer

    subroutine write_sourced_string( val, string )
    implicit none
    type(ids_summary_static_str_0d) :: val
        !< Type of IDS data structure, designed for sourced string data handling
    character(len=ids_string_length), intent(in) :: string

    allocate( val%value(1) )
    val%value = string
    allocate( val%source(1) )
    val%source = source

    return
    end subroutine write_sourced_string

    subroutine add_sourced_value( val, value )
    implicit none
    type(ids_summary_dynamic_flt_1d_root_parent_2) :: val
        !< Type of IDS data structure, designed for sourced float data handling
    real(IDS_real), intent(in) :: value

    if ( associated( val%value ) ) then
      val%value( slice_index ) = val%value( slice_index ) + value
    else
      call write_sourced_value( val, value )
    end if

    return
    end subroutine add_sourced_value

#if ( IMAS_MINOR_VERSION > 11 || IMAS_MAJOR_VERSION > 3 ) && GGD_MAJOR_VERSION > 0
    !!$> TODO: add to GGD itself (ids_grid_data)!
    !> Write a scalar data field given as a scalar data representation to a
    !! generic grid vector component IDS data fields.
    !!
    !! @note    The routine will make sure the required storage is
    !!          allocated, and will deallocate and re-allocate fields as
    !!          necessary.
    !! @note Available IDS vector component data fields:
    !!          - VEC_ALIGN_RADIAL_ID ( "Radial" ),
    !!          - VEC_ALIGN_DIAMAGNETIC_ID ( "Diamagnetic" ),
    !!          - VEC_ALIGN_PARALLEL_ID ( "Parallel" ),
    !!          - VEC_ALIGN_POLOIDAL_ID ( "Poloidal" ),
    !!          - VEC_ALIGN_TOROIDAL_ID ( "Toroidal" ),
    !!          - VEC_ALIGN_PHI_ID ( "Phi" ),
    !!          - VEC_ALIGN_R_MAJOR_ID ( "R" ),
    !!          - VEC_ALIGN_Z_ID ( "Z" )
    subroutine B2grid_Write_Data_Vector_Components( idsField_vcomp, &
             &   grid_index, grid_subset_index, vectorID, data)
    implicit none
    type(ids_generic_grid_vector_components), intent(inout) ::  &
        &   idsField_vcomp
        !< Type of IDS data structure, designed for handling data
        !< regarding vector components (parallel, poloidal etc.)
    integer, intent(in) :: grid_index           !< Grid index
    integer, intent(in) :: grid_subset_index    !< Base grid subset
                                                !< index
    character(len=*), intent(in) :: vectorID    !< Vector ID (e.g. )
                                                !< VEC_ALIGN_RADIAL_ID)
    real(IDS_real), intent(in) :: data(:)   !< Data field to be written
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
    case( VEC_ALIGN_DIAMAGNETIC_ID )
      !! Writing diamagnetic quantity
      !! Make sure the data field is properly allocated
      if ( associated( idsField_vcomp%diamagnetic ) ) then
        if ( .not. all( shape( idsField_vcomp%diamagnetic ) ==   &
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
    case( VEC_ALIGN_TOROIDAL_ID, VEC_ALIGN_PHI_ID )
      !! Writing toroidal quantity
      !! Make sure the data field is properly allocated
#if ( IMAS_MINOR_VERSION > 41 || IMAS_MAJOR_VERSION > 3 )
      if ( associated( idsField_vcomp%phi ) ) then
        if ( .not. all( shape( idsField_vcomp%phi ) ==  &
                    &   shape(data) )) then
          deallocate( idsField_vcomp%phi )
        end if
      end if
      !! If required, allocate storage
      if ( .not. associated( idsField_vcomp%phi ) ) then
        allocate(idsField_vcomp%phi( size(data, 1) ))
      end if
      !! copy toroidal data field
      idsField_vcomp%phi = data
#else
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
#endif
#if ( IMAS_MINOR_VERSION > 37 || ( IMAS_MINOR_VERSION == 37 && IMAS_MICRO_VERSION > 0 ) || IMAS_MAJOR_VERSION > 3 )
    case( VEC_ALIGN_R_MAJOR_ID )
      !! Writing major radius aligned quantity
      !! Make sure the data field is properly allocated
      if ( associated( idsField_vcomp%r ) ) then
        if ( .not. all( shape( idsField_vcomp%r ) ==  &
                    &   shape(data) )) then
          deallocate( idsField_vcomp%r )
        end if
      end if
      !! If required, allocate storage
      if ( .not. associated( idsField_vcomp%r ) ) then
        allocate(idsField_vcomp%r( size(data, 1) ))
      end if
      !! copy major radius aligned data field
      idsField_vcomp%r = data
    case( VEC_ALIGN_Z_ID )
      !! Writing vertical quantity
      !! Make sure the data field is properly allocated
      if ( associated( idsField_vcomp%z ) ) then
        if ( .not. all( shape( idsField_vcomp%z ) ==  &
                    &   shape(data) )) then
          deallocate( idsField_vcomp%z )
        end if
      end if
      !! If required, allocate storage
      if ( .not. associated( idsField_vcomp%z ) ) then
        allocate(idsField_vcomp%z( size(data, 1) ))
      end if
      !! copy vertical data field
      idsField_vcomp%z = data
#endif
    end select

    return
    end subroutine B2grid_Write_Data_Vector_Components
#endif

#if GGD_MAJOR_VERSION > 0
    !> Return unit vector along direction of given vector
    function unitVector(v) result(unitV)
    implicit none
    real(IDS_real), intent(in) :: v(:)  !< Vector
    real(IDS_real) :: unitV(size(v))    !< Unit vector

    unitV = v / sqrt( sum( v**2 ) )
    return
    end function unitVector
#endif

    subroutine write_sourced_value_root( val, value )
    implicit none
    type(ids_summary_dynamic_flt_1d_root) :: val
        !< Type of IDS data structure, designed for sourced float data handling
    real(IDS_real), intent(in) :: value

    allocate( val%value( num_slices ) )
    val%value( slice_index ) = value
    allocate( val%source(1) )
    val%source = source

    return
    end subroutine write_sourced_value_root

#if ( IMAS_MINOR_VERSION > 36 || IMAS_MAJOR_VERSION > 3 )
    subroutine write_sourced_rz( val, rvalue, zvalue )
    implicit none
    type(ids_summary_rz1d_dynamic) :: val
        !< Type of IDS data structure, designed for sourced float data handling
    real(IDS_real), intent(in) :: rvalue, zvalue

    allocate( val%r( num_slices ), val%z( num_slices ) )
    val%r( slice_index ) = rvalue
    val%z( slice_index ) = zvalue
    allocate( val%source(1) )
    val%source = source

    return
    end subroutine write_sourced_rz
#endif

    subroutine write_errored_value( val, value, error )
    implicit none
    type(ids_summary_dynamic_flt_1d_root) :: val
        !< Type of IDS data structure, designed for sourced float data handling
    real(IDS_real), intent(in) :: value
    real(IDS_real), intent(in) :: error

    allocate( val%value( num_slices ) )
    val%value( slice_index ) = value
    allocate( val%value_error_upper( num_slices ) )
    val%value_error_upper( slice_index ) = error
    allocate( val%source(1) )
    val%source = source

    return
    end subroutine write_errored_value

    subroutine write_sourced_value_root_parent_2( val, value )
        type(ids_summary_dynamic_flt_1d_root_parent_2) :: val
            !< Type of IDS data structure, designed for sourced float data handling
        real(ids_real), intent(in) :: value

        allocate( val%value( num_slices ) )
        val%value( slice_index ) = value
        allocate( val%source(1) )
        val%source = source

        return

    end subroutine write_sourced_value_root_parent_2

#elif defined(ITM_ENVIRONMENT_LOADED)

  logical, parameter, private :: INCLUDE_GHOST_CELLS = .false.

contains

  subroutine write_cpo(edgecpo)
    type (type_edge) :: edgecpo

    !! internal
    type(B2GridMap), save :: CPOmap
    logical, save :: CPOmapInitialized = .false.
    type(type_complexgrid_subgrid) :: sg_cell, sg_face, sg_bnd_core
    integer :: is, ns, nx, ny, i
    logical, parameter :: B2_WRITE_DATA = .true.
    real(ITM_R8), dimension(-1:ubound(crx,1),-1:ubound(crx,2),3,3) :: e
    integer :: iSgCore, iSgInnerMidplane, iSgOuterMidplane

    real(ITM_R8) :: tmpFace(-1:ubound(na, 1), -1:ubound(na, 2), 0:1)
    real(ITM_R8) :: tmpVx(-1:ubound(na, 1), -1:ubound(na, 2))
    character(len=13) :: spclabel

    !! allocate and init the CPO
    allocate(edgecpo%datainfo%dataprovider(1))
    edgecpo%datainfo%dataprovider="ITER"
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
       call species(is, spclabel, .false.)
       edgecpo%species(is+1)%label = spclabel
       edgecpo%species(is+1)%amn = am(is)
       edgecpo%species(is+1)%zn = zn(is)
       edgecpo%species(is+1)%zmin = zamin(is)
       edgecpo%species(is+1)%zmax = zamax(is)
    enddo

    !! set up the B2<->CPO mappings
    if (.not.CPOmapInitialized) &
      & call b2ITMCreateMap( nx,ny,crx(-1:nx,-1:ny,: ),cry(-1:nx,-1:ny,:),&
        & cflags,leftix,leftiy,rightix,rightiy, &
        & topix,topiy,bottomix,bottomiy, INCLUDE_GHOST_CELLS, CPOmap )
    CPOmapInitialized = .true.

    !! write grid & subgrids
    call b2ITMFillGridDescription( CPOmap, edgecpo%grid, &
        & nx,ny,crx(-1:nx,-1:ny,:),cry(-1:nx,-1:ny,:), &
        & leftix,leftiy,rightix,rightiy, &
        & topix,topiy,bottomix,bottomiy, &
        & nnreg, topcut, region, cflags, INCLUDE_GHOST_CELLS, vol, gs, qc )

    call xertst( geometryId( mpg, geo ) == GEOMETRY_SN,   &
        &   "write_cpo: can only do single null" )

    !! Write plasma state

    if ( B2_WRITE_DATA ) then

        call logmsg( LOGDEBUG, "b2mod_ual_io.write_cpo: writing plasma state" )

        iSgCore = gridFindSubGridByName( edgecpo%grid, "Core boundary" )
        iSgInnerMidplane = gridFindSubGridByName( edgecpo%grid, "Inner midplane" )
        iSgOuterMidplane = gridFindSubGridByName( edgecpo%grid, "Outer midplane" )

        !! ne
        call write_CPO_quantity( edgecpo%fluid%ne%value, edgecpo%fluid%ne%flux, ne, fne )
        call write_CPO_cell_scalar( edgecpo%fluid%ne%source, &
            &   b2CellData = sne(:,:,0) + sne(:,:,1)*ne )

        !! na
        allocate(edgecpo%fluid%ni(ns))
        do is = 1, ns
            call write_CPO_quantity( edgecpo%fluid%ni(is)%value, &
                &   edgecpo%fluid%ni(is)%flux,               &
                &   value = na(:,:,is-1),                    &
                &   flux = fna(:,:,:,is-1) )
            call write_CPO_cell_scalar( edgecpo%fluid%ni(is)%source, &
                &   b2CellData = sna(:,:,0,is-1) +               &
                &                sna(:,:,1,is-1)*na(:,:,is-1) )
        end do

        !! ue
        allocate(edgecpo%fluid%ve)
        allocate(edgecpo%fluid%ve%comps(1))
        allocate(edgecpo%fluid%ve%align(1))
        allocate(edgecpo%fluid%ve%alignid(1))
        edgecpo%fluid%ve%align(1) = VEC_ALIGN_PARALLEL
        edgecpo%fluid%ve%alignid(1) = VEC_ALIGN_PARALLEL_ID

        call write_CPO_cell_scalar( edgecpo%fluid%ve%comps(1)%value, &
            &   b2CellData = ue(:,:) )

        !! ua
        allocate(edgecpo%fluid%vi(ns))
        do is = 1, ns
            allocate(edgecpo%fluid%vi(is)%comps(1))
            allocate(edgecpo%fluid%vi(is)%align(1))
            allocate(edgecpo%fluid%vi(is)%alignid(1))
            edgecpo%fluid%vi(is)%align(1) = VEC_ALIGN_PARALLEL
            edgecpo%fluid%vi(is)%alignid(1) = VEC_ALIGN_PARALLEL_ID

            call write_CPO_cell_scalar( edgecpo%fluid%vi(is)%comps(1)%value, &
                &   b2CellData = ua(:,:,is-1) )
        end do

        !! te
        call write_CPO_quantity( edgecpo%fluid%te%value, &
            &   edgecpo%fluid%te%flux,                   &
            &   value = te/qe,                           &
            &   flux = fhe )

        !! ti
        allocate(edgecpo%fluid%ti(1))
        call write_CPO_quantity( edgecpo%fluid%ti(1)%value, &
            &   edgecpo%fluid%ti(1)%flux,                   &
            &   value = ti/qe,                              &
            &   flux = fhi )

        !! po
        call write_CPO_cell_scalar( edgecpo%fluid%po%value, po )

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
    subroutine write_CPO_quantity( values, fluxes, value, flux )
      use b2mod_interp
      type(type_complexgrid_scalar), pointer, intent(inout) :: values(:)
      type(type_complexgrid_vector), pointer, intent(inout) :: fluxes(:)
      real(ITM_R8), intent(in) :: value(-1:CPOmap%b2nx, -1:CPOmap%b2ny)
      real(ITM_R8), intent(in) :: flux(-1:CPOmap%b2nx, -1:CPOmap%b2ny, 0:1)
      real(ITM_R8), dimension(:), pointer :: cpodata
      real(ITM_R8) :: weight(-1:CPOmap%b2nx, -1:CPOmap%b2ny, TO_SELF:TO_TOP)
      integer i

      allocate(values(5))
      cpodata => b2ITMTransformDataB2ToCPO( edgecpo%grid, B2_SUBGRID_CELLS, &
               &  CPOmap, value )
      if ( size(cpodata) > 0 ) &
         & call gridWriteData( values(1), B2_SUBGRID_CELLS, cpodata )
      deallocate(cpodata)
      tmpFace = 0.0_ITM_R8
      do i = TO_SELF, TO_TOP
        weight(:,:,i)=vol(:,:)
      end do
      call value_on_faces(nx,ny,weight,value,tmpFace)
      cpodata => b2ITMTransformDataB2ToCPO( edgecpo%grid, &
               & iSgCore, CPOmap, tmpFace )
      if ( size(cpodata) > 0 ) &
         & call gridWriteData( values(2), iSgCore, cpodata )
      deallocate(cpodata)
      tmpVx = interpolateToVertices( CPOmap%b2nx, CPOmap%b2ny, &
               & VX_LOWERLEFT, value )
      cpodata => b2ITMTransformDataB2ToCPOVertex( edgecpo%grid, &
               & iSgInnerMidplane, CPOmap, tmpVx )
      if ( size(cpodata) > 0 ) &
         & call gridWriteData( values(3), iSgInnerMidplane, cpodata )
      deallocate(cpodata)
      cpodata => b2ITMTransformDataB2ToCPOVertex( edgecpo%grid, &
               & iSgOuterMidplane, CPOmap, tmpVx )
      if ( size(cpodata) > 0 ) &
         & call gridWriteData( values(4), iSgOuterMidplane, cpodata )
      deallocate(cpodata)
      cpodata => b2ITMTransformDataB2ToCPOVertex( edgecpo%grid, &
               & B2_SUBGRID_NODES, CPOmap, tmpVx )
      if ( size(cpodata) > 0 ) &
         & call gridWriteData( values(5), B2_SUBGRID_NODES, cpodata )
      deallocate(cpodata)
      allocate( fluxes(2) )
      call write_CPO_face_vector( fluxes(1), flux )
      call write_CPO_face_vector( fluxes(2), flux, subgridInd = iSgCore )
    end subroutine write_CPO_quantity


    !> Write a scalar B2 cell quantity to a complexgrid_scalar
    subroutine write_CPO_cell_scalar(scalar, b2CellData)
      type(type_complexgrid_scalar), intent(inout), pointer :: scalar(:)
      real(ITM_R8), intent(in) :: b2CellData(-1:CPOmap%b2nx, -1:CPOmap%b2ny)
      real(ITM_R8), dimension(:), pointer :: cpodata

      !! TODO: add checks whether already allocated
      allocate(scalar(1))
      cpodata => b2ITMTransformDataB2ToCPO( edgecpo%grid, &
               & B2_SUBGRID_CELLS, CPOmap, b2CellData )
      if ( size(cpodata) > 0 ) &
         & call gridWriteData( scalar(1), B2_SUBGRID_CELLS, cpodata )
      deallocate(cpodata)
    end subroutine write_CPO_cell_scalar


    !> Write a vector B2 cell quantity to a complexgrid_vector
    subroutine write_cell_vector(vector, align, alignid, vecdata)
      type(type_complexgrid_vector), intent(inout) :: vector
      real(ITM_R8), intent(in) :: vecdata(-1:CPOmap%b2nx, -1:CPOmap%b2ny, 0:2)
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
         cpodata => b2ITMTransformDataB2ToCPO(edgecpo%grid, &
                  & B2_SUBGRID_CELLS, CPOmap, vecdata(:,:,i-1))
         if ( size(cpodata) > 0 ) &
            & call gridWriteData( vector%comp(i), B2_SUBGRID_CELLS, cpodata )
         deallocate(cpodata)
      end do

    end subroutine write_cell_vector

    !> Write a vector B2 face quantity to a complexgrid_vector
    subroutine write_CPO_face_vector(vector, b2FaceData, subgridInd)
      type(type_complexgrid_vector), intent(inout) :: vector
      real(ITM_R8), intent(in) :: b2FaceData(-1:CPOmap%b2nx, -1:CPOmap%b2ny, 0:1)
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
!!$          cpodata => b2ITMTransformDataB2ToCPO(edgecpo%grid, &
!!$                   & B2_SUBGRID_EDGES_Y, CPOmap, b2FaceData)
!!$          if ( size(cpodata) > 0 ) &
!!$             & call gridWriteData( vector%comp(1), B2_SUBGRID_EDGES_Y, cpodata )
!!$          deallocate(cpodata)
!!$          cpodata => b2ITMTransformDataB2ToCPO(edgecpo%grid, &
!!$                   & B2_SUBGRID_EDGES_X, CPOmap, b2FaceData)
!!$          if ( size(cpodata) > 0 ) &
!!$             & call gridWriteData( vector%comp(2), B2_SUBGRID_EDGES_X, cpodata )
!!$          deallocate(cpodata)
!!$      else
!!$          allocate(vector%comp(1))
!!$          allocate(vector%align(1))
!!$          allocate(vector%alignid(1))
!!$
!!$          vector%align(1) = GRID_UNDEFINED
!!$          vector%alignid(1) = ""
!!$
!!$          cpodata => b2ITMTransformDataB2ToCPO(edgecpo%grid, &
!!$                   & subgridInd, CPOmap, b2FaceData)
!!$          if ( size(cpodata) > 0 ) &
!!$             & call gridWriteData( vector%comp(1), subgridInd, cpodata )
!!$          deallocate(cpodata)
!!$      end if

    end subroutine write_CPO_face_vector

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
                !! call xerrab ( "compute_Coordinate_Unit_Vectors: "// &
                !! & "not able to find poloidal neighbour for cell" )
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
                !! call xerrab ( "compute_Coordinate_Unit_Vectors: "// &
                !! & "not able to find toroidal neighbour for cell" )
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

#endif

end module b2mod_ual_io

!!!Local Variables:
!!! mode: f90
!!! End:
