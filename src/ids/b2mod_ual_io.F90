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
    use b2mod_b2cmrc
    use b2mod_geo
    use b2mod_work
    use b2mod_diag
    use b2mod_mwti
    use b2mod_rates
    use b2mod_plasma
    use b2mod_elements
    use b2mod_constants
    use b2mod_sources
    use b2mod_running_average
    use b2mod_feedback
    use b2mod_transport
    use b2mod_transport_nspecies
    use b2mod_anomalous_transport
    use b2mod_boundary_namelist
    use b2mod_neutrals_namelist
    use b2mod_user_namelist
    use b2mod_indirect
    use b2mod_external
    use b2mod_interp
    use b2mod_ipmain
    use b2mod_b2cmfs
    use b2mod_b2cmpb
    use b2mod_version
    use b2mod_grid_mapping
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
    use b2mod_b2plot &
     & , only : nxtl, nxtr, jxi, jxa, jsep
#endif
    use b2mod_b2plot_wall_loading
#ifdef B25_EIRENE
    use eirmod_ctrig
    use eirmod_cestim
    use eirmod_wneutrals
#ifdef IMAS
    use eirmod_comusr &
     & , only : lkindi, lkindm
    use eirmod_cinit &
     & , only : fort_lc
    use eirmod_comusr &
     & , only : natmi, nmoli, nioni, nmassa, nchara, nmassm, ncharm, &
     &          nprt, nchrgi, nchari
    use b2mod_b2plot &
     & , only : triangle_vol, ix_e2b, wklng, alloc_b2mod_b2plot_eirene
#endif
#elif defined(IMAS)
    use b2mod_b2plot &
     & , only : natmi
#endif
    use logging

#ifdef IMAS
    !! UAL Access
    use b2mod_ual_io_grid &
     & , only : INCLUDE_GHOST_CELLS
#if ( IMAS_MINOR_VERSION > 11 || IMAS_MAJOR_VERSION > 3 ) && GGD_MAJOR_VERSION > 0
    !! B2/CPO Mapping
    use b2mod_ual_io_data &
     & , only : b2_IMAS_Transform_Data_B2_To_IDS,       &
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
     &          ids_radiation, ids_equilibrium, ids_ids_properties,          &
     &          ids_code, ids_signal_int_1d, ids_signal_flt_1d,              &
     &          ids_generic_grid_scalar, ids_generic_grid_vector_components, &
     &          ids_generic_grid_dynamic,                                    &
     &          ids_plasma_composition_neutral_element,                      &
     &          ids_identifier_dynamic_aos3
#if ( IMAS_MAJOR_VERSION < 4 || ( IMAS_MAJOR_VERSION == 4 && IMAS_MINOR_VERSION < 1 ) )
    use ids_schemas &     ! IGNORE
     & , only : ids_dataset_description
#endif
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
#if ( IMAS_MINOR_VERSION > 29 || IMAS_MAJOR_VERSION > 3 )
    use ids_schemas &     ! IGNORE
     & , only : ids_code_with_timebase
#endif
#if ( IMAS_MINOR_VERSION > 30 || IMAS_MAJOR_VERSION > 3 )
    use ids_schemas &     ! IGNORE
     & , only : ids_divertors
#if ( IMAS_MAJOR_VERSION > 4 || ( IMAS_MAJOR_VERSION == 4 && IMAS_MINOR_VERSION > 0 ) )
    use al_midplane_identifier &       ! IGNORE
     & , only : set_midplane_identifier => set_identifier, &
     &          get_midplane_name => get_name
    use al_neutrals_identifier &       ! IGNORE
     & , only : set_neutral_type_identifier => set_identifier
    use al_radiation_identifier &      ! IGNORE
     & , only : set_radiation_identifier => set_identifier
    use al_edge_source_identifier &    ! IGNORE
     & , only : set_edge_source_identifier => set_identifier
    use al_plasma_source_identifier &  ! IGNORE
     & , only : set_plasma_source_identifier => set_identifier
#elif AL_MAJOR_VERSION > 4
    use al_midplane_identifier &       ! IGNORE
     & , only : midplane_identifier
    use al_neutrals_identifier &       ! IGNORE
     & , only : neutrals_identifier
    use al_radiation_identifier &      ! IGNORE
     & , only : radiation_identifier
    use al_edge_source_identifier &    ! IGNORE
     & , only : edge_source_identifier
#else
    use imas_midplane_identifier &     ! IGNORE
     & , only : midplane_identifier
    use imas_neutrals_identifier &     ! IGNORE
     & , only : neutrals_identifier
    use imas_radiation_identifier &    ! IGNORE
     & , only : radiation_identifier
    use imas_edge_source_identifier &  ! IGNORE
     & , only : edge_source_identifier
#endif
#endif
#if ( IMAS_MINOR_VERSION > 32 || IMAS_MAJOR_VERSION > 3 )
    use ids_utilities &   ! IGNORE
     & , only : ids_identifier, ids_identifier_static
#endif
#if ( IMAS_MINOR_VERSION > 36 || IMAS_MAJOR_VERSION > 3 )
    use ids_schemas &     ! IGNORE
     & , only : ids_summary_rz1d_dynamic
#endif
#if IMAS_MAJOR_VERSION > 3
    use ids_schemas &     ! IGNORE
     & , only : ids_code_constant, &
     &          ids_plasma_profiles, ids_plasma_sources, ids_plasma_transport
#if IMAS_MAJOR_VERSION > 4 || IMAS_MINOR_VERSION > 0
    use ids_schemas &     ! IGNORE
     & , only : ids_plasma_profiles_time_slice, &
     &          ids_plasma_sources_source_ggd, &
     &          ids_plasma_transport_model_ggd, &
     &          ids_summary_constant_flt_0d_2
#endif
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
  integer, save :: nneut !< Total number of IDS neutral species
  integer, save :: nsion !< Total number of IDS ion species
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
  integer, allocatable, save :: ionstt(:) !< Mapping array
                                          !< from B2-Eirene charged fluids to IDS ion states
  integer, allocatable, save :: istion(:) !< Number of IDS states for each ion
  integer, allocatable, save :: ispion(:,:) !< Mapping array
                                            !< from IDS ions and states to B2-Eirene ions
           !< ispion(i,j) contains the B2.5 species index for (ion i,state j) or
           !<                      the Eirene molecular ion index
#ifdef B25_EIRENE
  integer, allocatable, save :: isstat(:) !< Mapping array
                             !< from Eirene atoms and molecules to IDS neutral states
  integer, allocatable, save :: imneut(:) !< Mapping array
                             !< from Eirene molecules to IDS neutrals
  integer, allocatable, save :: imiion(:) !< Mapping array
                             !< from Eirene molecular ions to IDS ion sequences
#endif
  logical, parameter :: B2_WRITE_DATA = .true.
  integer, save :: ixpos(4), ifpos(4), iypos(4) !< Target positions
  integer, save :: idir(4), iysep(4), ixmid(4), ixmax(4)
  integer, save :: target_offset = 1
  integer, save :: nesepm_istra = -1
  integer, save :: pfrregno1 = 0
  integer, save :: pfrregno2 = 2
  integer, save :: use_eirene = 0
  real(IDS_real), save :: ndes = 0.0_IDS_real
  real(IDS_real), save :: ndes_sol = 0.0_IDS_real
  real(IDS_real), save :: nesepm_pfr = 0.0_IDS_real
  real(IDS_real), save :: nesepm_sol = 0.0_IDS_real
  real(IDS_real), save :: nepedm_sol = 0.0_IDS_real
  real(IDS_real), save :: volrec_sol = 0.0_IDS_real
  real(IDS_real), save :: private_flux_puff = 0.0_IDS_real
  real(IDS_real), save :: time  !< Generic time
  real(IDS_real), save :: time_slice_value   !< Time slice value
  real(IDS_real), save :: b0, r0, b0r0, z_eq
  real(IDS_real), save :: flux_expansion(4), extension_r(4), extension_z(4), &
      &                   wetted_area(4)
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
  character*32, save :: ADAS_git_version
  logical, save :: IDSmapInitialized = .false.
  logical, save :: eq_found
#ifndef NO_GETENV
  integer lenval, ierror
#ifndef USE_PXFGETENV
  intrinsic get_environment_variable
#endif
#endif
  type(B2GridMap), save :: IDSmap

contains

    subroutine IDS_init
    implicit none
    integer tvalues(8)
    integer i, istrai, p
    real(IDS_real) :: r_min, r_max, z_min, z_max
    logical, save :: IDS_initialized = .false.
    character*16 usrnam
    character*32 get_B25_hash
    character*32 get_ADAS_hash
    logical target_east, target_west
    logical streql
    external ipgeti, ipgetr, streql, usrnam, xertst
    external get_B25_hash, get_ADAS_hash

    if (IDS_initialized) return
    call ipgeti ('b2mndr_eirene', use_eirene)
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
    if (streql(b2frates_flag,'adas')) then
      ADAS_git_version = get_ADAS_hash()
      p = index(ADAS_git_version,'-')
      if (p.eq.0) then
        adas_version = trim(ADAS_git_version)
      else if (p.gt.1) then
        adas_version = ADAS_git_version(1:p-1)
      else
        adas_version = ''
      end if
    end if

    write(*,*) "Running b2CreateMap subroutine"
    !! Set up the B2<->IDS mappings
    nx = ubound( na, 1 )
    ny = ubound( na, 2 )
    geometryType = geometryId(nnreg, isymm, periodic_bc, topcut)
    configuration = geometryName( geometryType )
    if (.not.IDSmapInitialized)                                       &
       &  call b2CreateMap( nx, ny, crx( -1:nx, -1:ny, : ),           &
            &   cry( -1:nx, -1:ny, : ), cflags, leftix, leftiy,       &
            &   rightix, rightiy, topix, topiy, bottomix,bottomiy,    &
            &   INCLUDE_GHOST_CELLS, IDSmap, .false. )
    IDSmapInitialized = .true.

! Determine divertor plate generic information
    call ipgeti ('b2stbc_pfrregno1', pfrregno1)
    call ipgeti ('b2stbc_pfrregno2', pfrregno2)
    call ipgetr ('b2stbc_ndes', ndes)
    call ipgetr ('b2stbc_ndes_sol', ndes_sol)
    call ipgetr ('b2stbc_nesepm_pfr', nesepm_pfr)
    call ipgetr ('b2stbc_nesepm_sol', nesepm_sol)
    call ipgetr ('b2stbc_nepedm_sol', nepedm_sol)
    call ipgetr ('b2stbc_volrec_sol', volrec_sol)
    call ipgeti ('b2mwti_target_offset', target_offset)
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
    ntrgts=0
    if (nncut.eq.0) then
      if (geometryType.eq.GEOMETRY_LINEAR .or. &
        & geometryType.eq.GEOMETRY_CYLINDER) then
        if (boundary_namelist.ne.0) then
          target_east = .false.
          target_west = .false.
          do i=1,nbc
            if(bcchar(i).eq.'E'.and.bcpos(i).eq.-1) then
              target_east = target_east.or. &
                 &  bcene(i).eq.3.or.bcene(i).eq.12.or.bcene(i).eq.15
            end if
            if(bcchar(i).eq.'W'.and.bcpos(i).eq.nx) then
              target_west = target_west.or. &
                 &  bcene(i).eq.3.or.bcene(i).eq.12.or.bcene(i).eq.15
            end if
          end do
          if(target_west) then
            ntrgts=1
            plate_name(ntrgts) = bcchar(i)
            ixpos(ntrgts) = bcpos(i)
            itrg(ntrgts) = 1
            ixpos(ntrgts) = ixpos(ntrgts)+target_offset
            ifpos(ntrgts) = ixpos(ntrgts)+1
            iypos(ntrgts) = jsep
            idir(ntrgts) = -1
            ixmid(ntrgts) = jxa
            iysep(ntrgts) = jsep
            flux_expansion(ntrgts) =                                            &
                & ( wbbv(topix(jxa,jsep),topiy(jxa,jsep),0)/                    &
                &   wbbv(topix(jxa,jsep),topiy(jxa,jsep),3) )/                  &
                & ( wbbc(rightix(topix(bcpos(i),jsep),topiy(bcpos(i),jsep)),    &
                &        rightiy(topix(bcpos(i),jsep),topiy(bcpos(i),jsep)),0)/ &
                &   wbbc(rightix(topix(bcpos(i),jsep),topiy(bcpos(i),jsep)),    &
                &        rightiy(topix(bcpos(i),jsep),topiy(bcpos(i),jsep)),3) )
            r_max = max(maxval(crx(bc_list_x(1:bc_list_size(i),i),              &
                &                  bc_list_y(1:bc_list_size(i),i),1)),          &
                &       maxval(crx(bc_list_x(1:bc_list_size(i),i),              &
                &                  bc_list_y(1:bc_list_size(i),i),3)))
            r_min = min(minval(crx(bc_list_x(1:bc_list_size(i),i),              &
                &                  bc_list_y(1:bc_list_size(i),i),1)),          &
                &       minval(crx(bc_list_x(1:bc_list_size(i),i),              &
                &                  bc_list_y(1:bc_list_size(i),i),3)))
            z_max = max(maxval(cry(bc_list_x(1:bc_list_size(i),i),              &
                &                  bc_list_y(1:bc_list_size(i),i),1)),          &
                &       maxval(cry(bc_list_x(1:bc_list_size(i),i),              &
                &                  bc_list_y(1:bc_list_size(i),i),3)))
            z_min = min(minval(cry(bc_list_x(1:bc_list_size(i),i),              &
                &                  bc_list_y(1:bc_list_size(i),i),1)),          &
                &       minval(cry(bc_list_x(1:bc_list_size(i),i),              &
                &                  bc_list_y(1:bc_list_size(i),i),3)))
            extension_r(ntrgts) = r_max - r_min
            extension_z(ntrgts) = z_max - z_min
          end if
          if (target_east) then
            ntrgts = ntrgts+1
            itrg(ntrgts) = 2
            ixpos(ntrgts) = ixpos(ntrgts)-target_offset
            ifpos(ntrgts) = ixpos(ntrgts)
            iypos(ntrgts) = jsep
            idir(ntrgts) = 1
            ixmid(ntrgts) = jxa
            iysep(ntrgts) = jsep
            flux_expansion(ntrgts) =                                         &
                & ( wbbv(topix(jxa,jsep),topiy(jxa,jsep),0)/                 &
                &   wbbv(topix(jxa,jsep),topiy(jxa,jsep),3) )/               &
                & ( wbbc(topix(bcpos(i),jsep),topiy(bcpos(i),jsep),0)/       &
                &   wbbc(topix(bcpos(i),jsep),topiy(bcpos(i),jsep),3) )
            r_max = max(maxval(crx(bc_list_x(1:bc_list_size(i),i),           &
                &                  bc_list_y(1:bc_list_size(i),i),0)),       &
                &       maxval(crx(bc_list_x(1:bc_list_size(i),i),           &
                &                  bc_list_y(1:bc_list_size(i),i),2)))
            r_min = min(minval(crx(bc_list_x(1:bc_list_size(i),i),           &
                &                  bc_list_y(1:bc_list_size(i),i),0)),       &
                &       minval(crx(bc_list_x(1:bc_list_size(i),i),           &
                &                  bc_list_y(1:bc_list_size(i),i),2)))
            z_max = max(maxval(cry(bc_list_x(1:bc_list_size(i),i),           &
                &                  bc_list_y(1:bc_list_size(i),i),0)),       &
                &       maxval(cry(bc_list_x(1:bc_list_size(i),i),           &
                &                  bc_list_y(1:bc_list_size(i),i),2)))
            z_min = min(minval(cry(bc_list_x(1:bc_list_size(i),i),           &
                &                  bc_list_y(1:bc_list_size(i),i),0)),       &
                &       minval(cry(bc_list_x(1:bc_list_size(i),i),           &
                &                  bc_list_y(1:bc_list_size(i),i),2)))
            extension_r(ntrgts) = r_max - r_min
            extension_z(ntrgts) = z_max - z_min
          end if
        end if
      else if (geometryType.eq.GEOMETRY_LIMITER) then
        ntrgts = 2
        itrg(1) = 1
        ixpos(1) = -1+target_offset
        ifpos(1) = 0
        iypos(1) = topcut(1)
        idir(1) = -1
        iysep(1) = topcut(1)
        itrg(2) = 2
        ixpos(2) = nx-target_offset
        ifpos(2) = nx
        iypos(2) = topcut(1)
        idir(2) = 1
        iysep(2) = topcut(1)
        ixmid(1) = jxa
        ixmid(2) = jxa
        flux_expansion(1) = ( wbbv(topix(jxa,jsep),topiy(jxa,jsep),0)/   &
            &                 wbbv(topix(jxa,jsep),topiy(jxa,jsep),3) )/ &
            &               ( wbbc(topix(0,jsep),topiy(0,jsep),0)/       &
            &                 wbbc(topix(0,jsep),topiy(0,jsep),3) )
        flux_expansion(2) = ( wbbv(topix(jxa,jsep),topiy(jxa,jsep),0)/   &
            &                 wbbv(topix(jxa,jsep),topiy(jxa,jsep),3) )/ &
            &               ( wbbc(topix(nx,jsep),topiy(nx,jsep),0)/     &
            &                 wbbc(topix(nx,jsep),topiy(nx,jsep),3) )
      end if
    else
      ntrgts = 2*nncut
      plate_name(1) = 'LI'
      itrg(1) = 1
      ixpos(1) = -1+target_offset
      ifpos(1) = 0
      iypos(1) = topcut(1)
      idir(1) = -1
      iysep(1) = topcut(1)
      r_max = max(maxval(crx(0,0:ny-1,0)),maxval(crx(0,0:ny-1,2)))
      r_min = min(minval(crx(0,0:ny-1,0)),minval(crx(0,0:ny-1,2)))
      z_max = max(maxval(cry(0,0:ny-1,0)),maxval(cry(0,0:ny-1,2)))
      z_min = min(minval(cry(0,0:ny-1,0)),minval(cry(0,0:ny-1,2)))
      extension_r(1) = r_max - r_min
      extension_z(1) = z_max - z_min
      if (nncut.eq.2) then
        itrg(2) = 2
        plate_name(2) = 'UI'
        ixpos(2) = nxtl-target_offset
        ifpos(2) = nxtl
        iypos(2) = topcut(2)
        idir(2) = 1
        iysep(2) = topcut(2)
        r_max = max(maxval(crx(nxtl,0:ny-1,0)),maxval(crx(nxtl,0:ny-1,2)))
        r_min = min(minval(crx(nxtl,0:ny-1,0)),minval(crx(nxtl,0:ny-1,2)))
        z_max = max(maxval(cry(nxtl,0:ny-1,0)),maxval(cry(nxtl,0:ny-1,2)))
        z_min = min(minval(cry(nxtl,0:ny-1,0)),minval(cry(nxtl,0:ny-1,2)))
        extension_r(2) = r_max - r_min
        extension_z(2) = z_max - z_min
        itrg(3) = 3
        plate_name(3) = 'UO'
        ixpos(3) = nxtr+target_offset
        ifpos(3) = nxtr+1
        iypos(3) = topcut(2)
        idir(3) = -1
        iysep(3) = topcut(2)
        r_max = max(maxval(crx(nxtr,0:ny-1,1)),maxval(crx(nxtr,0:ny-1,3)))
        r_min = min(minval(crx(nxtr,0:ny-1,1)),minval(crx(nxtr,0:ny-1,3)))
        z_max = max(maxval(cry(nxtr,0:ny-1,1)),maxval(cry(nxtr,0:ny-1,3)))
        z_min = min(minval(cry(nxtr,0:ny-1,1)),minval(cry(nxtr,0:ny-1,3)))
        extension_r(3) = r_max - r_min
        extension_z(3) = z_max - z_min
      end if
      itrg(ntrgts) = ntrgts
      plate_name(ntrgts) = 'LO'
      ixpos(ntrgts) = nx-target_offset
      ifpos(ntrgts) = nx
      iypos(ntrgts) = topcut(1)
      idir(ntrgts) = 1
      iysep(ntrgts) = topcut(1)
      r_max = max(maxval(crx(nx,0:ny-1,0)),maxval(crx(nx,0:ny-1,2)))
      r_min = min(minval(crx(nx,0:ny-1,0)),minval(crx(nx,0:ny-1,2)))
      z_max = max(maxval(cry(nx,0:ny-1,0)),maxval(cry(nx,0:ny-1,2)))
      z_min = min(minval(cry(nx,0:ny-1,0)),minval(cry(nx,0:ny-1,2)))
      extension_r(ntrgts) = r_max - r_min
      extension_z(ntrgts) = z_max - z_min
      if (nncut.eq.1) then
        ixmid(1) = jxa
        flux_expansion(1) = ( wbbv(topix(jxa,jsep),topiy(jxa,jsep),0)/   &
            &                 wbbv(topix(jxa,jsep),topiy(jxa,jsep),3) )/ &
            &               ( wbbc(topix(0,jsep),topiy(0,jsep),0)/       &
            &                 wbbc(topix(0,jsep),topiy(0,jsep),3) )
        ixmid(2) = jxa
        flux_expansion(2) = ( wbbv(topix(jxa,jsep),topiy(jxa,jsep),0)/   &
            &                 wbbv(topix(jxa,jsep),topiy(jxa,jsep),3) )/ &
            &               ( wbbc(topix(nx,jsep),topiy(nx,jsep),0)/     &
            &                 wbbc(topix(nx,jsep),topiy(nx,jsep),3) )
      elseif (nnreg(0).eq.7) then
        ixmid(1) = jxa
        flux_expansion(1) = ( wbbv(topix(jxa,jsep),topiy(jxa,jsep),0)/   &
            &                 wbbv(topix(jxa,jsep),topiy(jxa,jsep),3) )/ &
            &               ( wbbc(topix(0,jsep),topiy(0,jsep),0)/       &
            &                 wbbc(topix(0,jsep),topiy(0,jsep),3) )
        ixmid(2) = jxa
        flux_expansion(2) = ( wbbv(topix(jxa,jsep),topiy(jxa,jsep),0)/   &
            &                 wbbv(topix(jxa,jsep),topiy(jxa,jsep),3) )/ &
            &               ( wbbc(topix(nxtl,jsep),topiy(nxtl,jsep),0)/ &
            &                 wbbc(topix(nxtl,jsep),topiy(nxtl,jsep),3) )
        ixmid(3) = jxa
        flux_expansion(3) = ( wbbv(topix(jxa,jsep),topiy(jxa,jsep),0)/   &
            &                 wbbv(topix(jxa,jsep),topiy(jxa,jsep),3) )/ &
            &               ( wbbc(topix(nxtr,jsep),topiy(nxtr,jsep),0)/ &
            &                 wbbc(topix(nxtr,jsep),topiy(nxtr,jsep),3) )
        ixmid(4) = jxa
        flux_expansion(4) = ( wbbv(topix(jxa,jsep),topiy(jxa,jsep),0)/   &
            &                 wbbv(topix(jxa,jsep),topiy(jxa,jsep),3) )/ &
            &               ( wbbc(topix(nx,jsep),topiy(nx,jsep),0)/     &
            &                 wbbc(topix(nx,jsep),topiy(nx,jsep),3) )
      else
        if (topcut(1).lt.topcut(2)) then
          ixmid(1) = jxa
          flux_expansion(1) = &
              &  ( wbbv(topix(jxa,topcut(1)),topiy(jxa,topcut(1)),0)/    &
              &    wbbv(topix(jxa,topcut(1)),topiy(jxa,topcut(1)),3) )/  &
              &  ( wbbc(topix(0,topcut(1)),topiy(0,topcut(1)),0)/        &
              &    wbbc(topix(0,topcut(1)),topiy(0,topcut(1)),3) )
        else
          ixmid(1) = jxi
          flux_expansion(1) = &
              & ( wbbv(topix(jxi,topcut(1)),topiy(jxi,topcut(1)),0)/    &
              &   wbbv(topix(jxi,topcut(1)),topiy(jxi,topcut(1)),3) )/  &
              & ( wbbc(topix(0,topcut(1)),topiy(0,topcut(1)),0)/        &
              &   wbbc(topix(0,topcut(1)),topiy(0,topcut(1)),3) )
        end if
        if (topcut(2).lt.topcut(1)) then
          ixmid(2) = jxa
          flux_expansion(2) = &
              & ( wbbv(topix(jxa,topcut(2)),topiy(jxa,topcut(2)),0)/    &
              &   wbbv(topix(jxa,topcut(2)),topiy(jxa,topcut(2)),3) )/  &
              & ( wbbc(topix(nxtl,topcut(2)),topiy(nxtl,topcut(2)),0)/  &
              &   wbbc(topix(nxtl,topcut(2)),topiy(nxtl,topcut(2)),3) )
        else
          ixmid(2) = jxi
          flux_expansion(2) = &
              & ( wbbv(topix(jxi,topcut(2)),topiy(jxi,topcut(2)),0)/    &
              &   wbbv(topix(jxi,topcut(2)),topiy(jxi,topcut(2)),3) )/  &
              & ( wbbc(topix(nxtl,topcut(2)),topiy(nxtl,topcut(2)),0)/  &
              &   wbbc(topix(nxtl,topcut(2)),topiy(nxtl,topcut(2)),3) )
        end if
        ixmid(3) = jxa
        flux_expansion(3) = &
            & ( wbbv(topix(jxa,topcut(2)),topiy(jxa,topcut(2)),0)/    &
            &   wbbv(topix(jxa,topcut(2)),topiy(jxa,topcut(2)),3) )/  &
            & ( wbbc(rightix(topix(nxtr,topcut(2)),                   &
            &                topiy(nxtr,topcut(2))),                  &
            &        rightiy(topix(nxtr,topcut(2)),                   &
            &                topiy(nxtr,topcut(2))),0)/               &
            &   wbbc(rightix(topix(nxtr,topcut(2)),                   &
            &                topiy(nxtr,topcut(2))),                  &
            &        rightiy(topix(nxtr,topcut(2)),                   &
            &                topiy(nxtr,topcut(2))),3) )
        ixmid(4) = jxa
        flux_expansion(4) = &
            & ( wbbv(topix(jxa,topcut(1)),topiy(jxa,topcut(1)),0)/    &
            &   wbbv(topix(jxa,topcut(1)),topiy(jxa,topcut(1)),3) )/  &
            & ( wbbc(topix(nx,topcut(1)),topiy(nx,topcut(1)),0)/      &
            &   wbbc(topix(nx,topcut(1)),topiy(nx,topcut(1)),3) )
      endif
    end if

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
            &   edge_profiles, edge_sources, edge_transport, &
#if IMAS_MAJOR_VERSION > 3
            &   plasma_profiles, plasma_sources, plasma_transport, &
#endif
            &   radiation, &
#if ( IMAS_MAJOR_VERSION < 4 || ( IMAS_MAJOR_VERSION == 4 && IMAS_MINOR_VERSION < 1 ) )
            &   description, &
#endif
            &   equilibrium, &
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
#if IMAS_MAJOR_VERSION > 3
            &   time_IN, time_step_IN, shot, database, &
#else
            &   time_IN, time_step_IN, shot, run, database, version, &
#endif
            &   new_eq_ggd, &
            &   time_slice_ind_IN, num_time_slices_IN )
#ifdef NO_OPT
!DIR$ NOOPTIMIZE
#endif
        use b2mod_dimensions
        implicit none
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
#if IMAS_MAJOR_VERSION > 3
        type (ids_plasma_profiles) :: plasma_profiles !< IDS designed to
            !< store data on plasma profiles
        type (ids_plasma_sources) :: plasma_sources !< IDS designed to store
            !< data on plasma sources. Energy terms correspond to the full
            !< kinetic energy equation (i.e. the energy flux takes into account
            !< the energy transported by the particle flux)
        type (ids_plasma_transport) :: plasma_transport !< IDS designed to store
            !< data on plasma transport. Energy terms correspond to the
            !< full kinetic energy equation (i.e. the energy flux takes into
            !< account the energy transported by the particle flux)
#endif
        type (ids_radiation) :: radiation !< IDS designed to store
            !< data on radiation emitted by the plasma species
#if ( IMAS_MAJOR_VERSION < 4 || ( IMAS_MAJOR_VERSION == 4 && IMAS_MINOR_VERSION < 1 ) )
        type (ids_dataset_description) :: description !< IDS designed to store
            !< a description of the simulation
#endif
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
#if IMAS_MAJOR_VERSION < 4
        integer, intent(in) :: run
        character(len=24), intent(in) :: version
#endif
        integer, intent(in) :: shot
        character(len=24), intent(in) :: database
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
        integer :: n_process !< Number of radiation processes handled
        integer :: is, js, ks !< Species indices (iterators)
        integer :: ii, jj !< Iterators
        integer :: i      !< Iterator
        integer :: j      !< Iterator
        integer :: k      !< Iterator
        integer :: ix     !< Iterator
        integer :: iy     !< Iterator
        integer :: is1    !< First ion of an isonuclear sequence
        integer :: is2    !< Last ion of an isonuclear sequence
        integer :: icnt   !< Boundary cell counter
        integer :: ib     !< Boundary condition index
        integer :: isep(2) !< Array of separatrix regions
#ifdef B25_EIRENE
        integer :: ind    !< Non-standard surface index in resolved list
        integer :: ias    !< Starting index for non-standard surface in resolved list
        integer :: iss    !< State index
        integer :: iatm   !< Atom iterator
        integer :: imol   !< Molecule iterator
        integer :: nelems !< Number of elements present in a molecule or molecular ion
        integer :: p      !< Dummy integer
        integer :: ixx, iyy
#endif
        integer :: nscx, iscx(0:nscxmax-1)
#if ( IMAS_MINOR_VERSION > 38 || IMAS_MAJOR_VERSION > 3 )
        real(IDS_real) :: flxFace( -1:ubound( na, 1), -1:ubound( na, 2), 0:1, 0:1 )
        real(IDS_real) :: totflux( -1:ubound( na, 1), -1:ubound( na, 2), 0:1, 0:1 )
#endif
        real(IDS_real) :: hlpFace( -1:ubound( na, 1), -1:ubound( na, 2), 0:1 )
        real(IDS_real) :: tmpFace( -1:ubound( na, 1), -1:ubound( na, 2), 0:1 )
        real(IDS_real) :: totFace( -1:ubound( na, 1), -1:ubound( na, 2), 0:1 )
        real(IDS_real) :: tmpVx( -1:ubound( na, 1), -1:ubound( na, 2) )
        real(IDS_real) :: tmpCv( -1:ubound( na, 1), -1:ubound( na, 2) )
        real(IDS_real) :: totCv( -1:ubound( na, 1), -1:ubound( na, 2) )
        real(IDS_real) :: lnlam( -1:ubound( na, 1), -1:ubound( na, 2) )
        real(IDS_real) :: pz( -1:ubound( na, 1), -1:ubound( na, 2) )
        real(IDS_real) :: zeff( -1:ubound( na, 1), -1:ubound( na, 2) )
        real(IDS_real) :: time_step !< Time step
        real(IDS_real) :: nibnd, frac, u, v,                                 &
            &             qtot, qetot, qitot, qmax, qemax, qimax, lambda,    &
            &             vtor, nisep, nasum
        real(IDS_real) :: power_convected(4),                                &
            &             power_conducted(4), power_neutrals(4),             &
            &             power_incident(4), power_flux_peak(4),             &
            &             power_recomb_neutrals(4), power_radiated(4),       &
            &             recycled_flux(4)
#if ( IMAS_MAJOR_VERSION > 4 || ( IMAS_MAJOR_VERSION == 4 && IMAS_MINOR_VERSION > 0 ) )
        real(IDS_real) :: nesum
#endif
        real(IDS_real), allocatable :: wrdtrg(:,:,:)
#ifdef B25_EIRENE
        real(IDS_real), allocatable :: un0(:,:,:,:), um0(:,:,:,:)
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
#if ( IMAS_MAJOR_VERSION > 4 || ( IMAS_MAJOR_VERSION == 4 && IMAS_MINOR_VERSION > 0 ) )
        type(ids_plasma_profiles_time_slice) :: profiles_ggd
        type(ids_plasma_sources_source_ggd), allocatable :: sources_ggd(:)
        type(ids_plasma_transport_model_ggd), allocatable :: transport_ggd(:)
#endif
#endif
        integer :: istrai !< Stratum iterator
        real(IDS_real),   &
            &   dimension( -1:ubound( crx, 1 ), -1:ubound( crx, 2), 3, 3) :: e
        real(IDS_real) :: pb( -1:ubound( na, 1), -1:ubound( na, 2) )
        real(IDS_real) :: pe( -1:ubound( na, 1), -1:ubound( na, 2) )
        real(IDS_real) :: ue( -1:ubound( na, 1), -1:ubound( na, 2) )
#endif
#if ( ( IMAS_MINOR_VERSION > 38 || IMAS_MAJOR_VERSION > 3 ) && defined(B25_EIRENE) )
        integer, parameter :: nsources = 13
#else
        integer, parameter :: nsources = 12
#endif
        integer, save :: ncall = 0
        integer, save :: ismain = 1
        integer, save :: ismain0 = 0
        integer, save :: ue_style = 2
        integer, save :: ids_from_43 = 0
        integer, save :: balance_netcdf = 0
        integer, save :: drift_style
        real(IDS_real), save :: dtim = 1.0_IDS_real
        real(IDS_real), save :: neutral_sources_rescale = 1.0_IDS_real
        real(IDS_real), save :: BoRiS = 0.0_IDS_real
        real(IDS_real), save :: flux_multiplier
        character*132 radiation_commit
        character*256 filename
        character*5 hlp_frm
        logical match_found, streql, exists, wrong_flow
#ifdef B25_EIRENE
        character(len=132) :: mol_label !< Molecule species label (e.g. D2)
        character(len=132) :: ion_label !< Ion species label (e.g. D+1)
        logical, allocatable :: in_species(:)
#endif
        external b2xpne, b2xpni, b2xppb, b2xppe, b2xppz, b2xpve, b2xzef
        external b2sral, b2spcx, b2tral, b2tanml, b2ptrdl
        external ipgetr, ipgeti, species, streql, xerrab, xertst
        external find_file

        !! ===  SET UP IDS ===
        write(0,*) "Setting data for edge_profiles IDS"
        if (ncall.eq.0) then
          call ipgetr ('b2news_BoRiS', BoRiS)
          call ipgeti ('b2mndr_ismain', ismain)
          call ipgeti ('b2sigp_style', ue_style)
          call ipgetr ('b2mndr_dtim', dtim)
          call ipgeti ('ids_from_43', ids_from_43)
          call ipgetr ('b2mndr_rescale_neutrals_sources', &
              &                 neutral_sources_rescale)
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
          if (redef_gmtry.eq.0) then
            drift_style = 1
          else
            drift_style = 2
          end if
          call ipgeti ('b2tfnb_drift_style', drift_style)
        end if
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

        call IDS_init
        ns = size( na, 3 )
        call alloc_b2mod_user(nx, ns, nlim, nmol, ntns)
        call b2xzef (nx, ny, ns, rz2, na, ne, zeff)
        call b2xppz (nx, ny, ns, ne, na, te, ti, pz)

!   ..find nscx, iscx
!     (indices for neutral hydrogen species)
        nscx = 0
        do is = 0, ns-1
          if (is_neutral(is).and.nint(zn(is)).eq.1) then
            call xertst (nscx.lt.nscxmax, 'too many neutral hydrogen species')
            iscx(nscx) = is
            nscx = nscx+1
          endif
        enddo
        do k = nscx, nscxmax-1
          iscx(k) = -1
        enddo
        call b2xpni (nx, ny, ns, na, ni)
        call b2xpne (nx, ny, ns, rza, na, ne_ext, ne)
        call b2xpne (nx, ny, ns, rz2, na, ne2_ext, ne2)
!   ..compute flux limit coefficients
        call b2tral (nx, ny, ns, nscx, nscxmax, iscx, ismain,                     &
            &        bb, conn, vol, gs, hx, hy, hz, qz, qc, qs,                   &
            &        pbs, crx, cry, bzb, lnlam,                                   &
            &        fch, na, ua, te, ti, po, ne, ni, ne2, chvemx, chvimx,        &
            &        cdna, cdpa, cddi, cdde, cvla, cvsa, chce, chve, chci,        & !som 03.11.21
            &        chvi, csig, csigin, calf, cthe, cthi,                        &
            &        cdnahz, cdpahz, cvlahz, cvmahz, cvsahz, cvsa_cl, cvsahz_cl,  &
            &        cvsahz_drho, cvsa_drho,                                      & !som 18.08.21
            &        chci_al_ast, chci_a, cvsahz_hAdp_albe, cvsahz_hBdp_al,       & !som 08.10.21
            &        cvsa_hAdp_albe, cvsa_hBdp_al,                                & !som 20.10.21
            &        cvsahz_uAdp_albe, cvsahz_uBdp_al,                            & !som 22.12.23
            &        cvsa_uAdp_albe, cvsa_uBdp_al,                                & !som 26.04.24
            &        fllim0fhi, fllimvisc, csig_cl, calf_cl,                      &
            &        csig_stoch, chce_stoch)
!  ..compute log-log charge exchange rate coefficients
        do k = 0, nscx-1
           call b2spcx (nx, ny, ns, ev, am(iscx(k)), ti, ne, rlcx(-1,-1,0,0,k))
        enddo
!   ..compute sources
        call b2sral (nx, ny, ns,                                                 &
            &        nscx, nscxmax, 0, ns, iscx, ismain, ismain0,                &
            &        dtim, BoRiS, facdrift, fac_ExB, fac_vis,                    &
            &        vol, hx, hy, hy1, hz, qz, qc, qs, gs, pbs, bb, lnlam,       &
            &        na, ua,                                                     &
            &        uadia, vedia, vadia, wadia, veecrb, vaecrb, ve, wedia,      &
            &        te, ti, po, ne, ni, kinrgy, floe_noc, floi_noc,             &
            &        fna, fna_32, fna_52, fni_32, fni_52, fne_32, fne_52,        &
            &        fna_mdf, fhe_mdf, fhi_mdf, fna_fcor, fna_nodrift, fna_he,   &
            &        fhe, fhi, fhm, fht, fnaPSch, fhePSch, fhiPSch, fch,         &
            &        fchanml_a, fchinert_a, fchvispar_a, fchvisper_a, fchvisq_a, & !srv 08.09.21
            &        fchin_a, fchdia, fchin, fch_p, fchvispar, fchvisper,        &
            &        fchvisq, fchinert, fchanml,                                 &
            &        fna_eir, fne_eir, fhe_eir, fhi_eir,                         &
            &        cdna, cdpa, cvsa_cl, cvsa_drho,                             & !som 18.08.21
            &        chci_al_ast, chci_a, cvsa_hAdp_albe, cvsa_hBdp_al,          & !som 20.10.21
            &        cvsa_uAdp_albe, cvsa_uBdp_al,                               & !som 26.04.24
            &        cvla, chce, chve, chci, chvi, calf, csigin,                 &
            &        rlsa, rlra, rlqa, rlcx, rlrd, rlbr,                         &
            &        rlza, rlz2, rlpt, rlpi,                                     &
            &        rza, rz2, rpt, rpi, sna, smo, smq, she, shi, sch, sne,      &
            &        wrong_flow, .false.)
#ifdef B25_EIRENE
        if (use_eirene.ne.0) then
          filename=fort_lc//'46'
          call find_file(filename,exists)
          if(exists.and.use_eirene.ne.0) then
            if (.not.allocated(xtrian)) then ! Eirene has not been called
              open(unit=46,file=filename)
              call ntread
            else if (.not.allocated(triangle_vol)) then ! this is the first pass
              call alloc_b2mod_b2plot_eirene(natmi,nmoli,nioni,ntrii,wklng)
              call compute_triangle_area
              call compute_triangle_vol
            end if
          end if
        end if
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
#if IMAS_MAJOR_VERSION > 3
        call write_ids_properties( plasma_profiles%ids_properties, &
          &  homogeneous_time )
        call write_ids_properties( plasma_transport%ids_properties, &
          &  homogeneous_time )
        call write_ids_properties( plasma_sources%ids_properties, &
          &  homogeneous_time )
#endif
        call write_ids_properties( radiation%ids_properties, &
          &  homogeneous_time )
#if ( IMAS_MAJOR_VERSION < 4 || ( IMAS_MAJOR_VERSION == 4 && IMAS_MINOR_VERSION < 1 ) )
        call write_ids_properties( description%ids_properties, 2 )
#endif
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
        if (streql(b2frates_flag,'adas')) then
          radiation_commit = 'B25 : '//trim(B25_git_version)// &
                      &  ' + ADAS : '//trim(ADAS_git_version)
        else
          radiation_commit = B25_git_version
        endif
        call write_ids_code( edge_profiles%code, code_commit, code_description )
        call write_ids_code( edge_transport%code, code_commit, code_description )
        call write_ids_code( edge_sources%code, code_commit, code_description )
#if IMAS_MAJOR_VERSION > 3
        call write_ids_code( plasma_profiles%code, code_commit, code_description )
        call write_ids_code( plasma_transport%code, code_commit, code_description )
        call write_ids_code( plasma_sources%code, code_commit, code_description )
#endif
        call write_ids_code( radiation%code, radiation_commit, code_description )
#if ( IMAS_MINOR_VERSION > 21 || IMAS_MAJOR_VERSION > 3 )
        call write_ids_code( summary%code, code_commit, code_description )
#endif
#if ( IMAS_MINOR_VERSION > 30 || IMAS_MAJOR_VERSION > 3 )
        call write_ids_code( divertors%code, code_commit, code_description )
#endif
#if ( IMAS_MAJOR_VERSION == 4 && IMAS_MINOR_VERSION == 0 )
        call write_ids_code_constant( description%code, code_commit, code_description )
#endif
        allocate( edge_transport%model(1) )
        call write_model_identifier( edge_transport%model(1)%identifier )
        if (ids_from_43.eq.0) then
          flux_multiplier = 1.5_IDS_real + BoRiS
        else
          flux_multiplier = 2.5_IDS_real
        end if
        edge_transport%model(1)%flux_multiplier = flux_multiplier
        code_commit = B25_git_version
        code_description = "Snapshot IDS written by b2mod_ual_io routine"
#if ( IMAS_MINOR_VERSION > 29 || IMAS_MAJOR_VERSION > 3 )
        call write_ids_code_timed( edge_transport%model(1)%code, &
            & code_commit, code_description )
#endif
#if IMAS_MAJOR_VERSION > 3
        allocate( plasma_transport%model(1) )
        call write_model_identifier( plasma_transport%model(1)%identifier )
        plasma_transport%model(1)%flux_multiplier = flux_multiplier
        call write_ids_code_timed( plasma_transport%model(1)%code, &
            & code_commit, code_description )
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
#else
        allocate( plasma_profiles%time(num_time_slices) )
        plasma_profiles%time(time_sind) = time
        allocate( plasma_transport%time(num_time_slices) )
        plasma_transport%time(time_sind) = time
        allocate( plasma_sources%time(num_time_slices) )
        plasma_sources%time(time_sind) = time
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
        if (use_eirene.ne.0) then
          n_process = 4
        else
          n_process = 2
        end if
        allocate( radiation%process(n_process) )
        do j = 1, n_process
#if ( IMAS_MAJOR_VERSION > 4 || ( IMAS_MAJOR_VERSION == 4 && IMAS_MINOR_VERSION > 0 ) )
          if (j.eq.2) cycle ! Now taken care of by set_identifier API
#endif
          allocate( radiation%process(j)%identifier%name(1) )
          allocate( radiation%process(j)%identifier%description(1) )
        end do
        radiation%process(1)%identifier%index = 2
        radiation%process(1)%identifier%name = 'line_radiation'
        radiation%process(1)%identifier%description = 'Line and rec. rad. from B2.5 species'
#if ( IMAS_MAJOR_VERSION > 4 || ( IMAS_MAJOR_VERSION == 4 && IMAS_MINOR_VERSION > 0 ) )
        call set_radiation_identifier( radiation%process(2)%identifier, 'bremsstrahlung' )
#elif ( IMAS_MAJOR_VERSION > 3 || IMAS_MINOR_VERSION > 30 )
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
        if (use_eirene.ne.0) then
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
#if IMAS_MAJOR_VERSION > 3
        allocate( plasma_profiles%ggd( num_time_slices ) )
        allocate( plasma_profiles%grid_ggd( num_time_slices ) )
        allocate( plasma_transport%grid_ggd( num_time_slices ) )
        allocate( plasma_sources%grid_ggd( num_time_slices ) )
#endif
#if ( IMAS_MINOR_VERSION > 21 || IMAS_MAJOR_VERSION > 3 )
        allocate( radiation%grid_ggd( num_time_slices ) )
#endif
#endif
        allocate( edge_transport%model(1)%ggd( num_time_slices ) )
#if IMAS_MAJOR_VERSION > 3
        allocate( plasma_transport%model(1)%ggd( num_time_slices ) )
#endif
        allocate( edge_sources%source(nsources) )
        do is = 1, nsources
          allocate( edge_sources%source(is)%ggd( num_time_slices ) )
#if ( IMAS_MAJOR_VERSION == 3 && IMAS_MINOR_VERSION < 31 )
          allocate( edge_sources%source(is)%identifier%name(1) )
          allocate( edge_sources%source(is)%identifier%description(1) )
#endif
        end do
#if ( IMAS_MAJOR_VERSION > 4 || ( IMAS_MAJOR_VERSION == 4 && IMAS_MINOR_VERSION > 0 ) )
        !! Total sources
        call set_edge_source_identifier( edge_sources%source(1)%identifier, "Total" )
        !! Background sources
        call set_edge_source_identifier( edge_sources%source(2)%identifier, "Background" )
        !! Prescribed sources
        call set_edge_source_identifier( edge_sources%source(3)%identifier, "Prescribed" )
        !! Time derivative
        call set_edge_source_identifier( edge_sources%source(4)%identifier, "Time derivative" )
        !! Atomic ionization
        call set_edge_source_identifier( edge_sources%source(5)%identifier, "Atomic ionization" )
        !! Molecular ionization
        call set_edge_source_identifier( edge_sources%source(6)%identifier, "Molecular ionization" )
        !! Ionization
        call set_edge_source_identifier( edge_sources%source(7)%identifier, "Ionization" )
        !! Recombination
        call set_edge_source_identifier( edge_sources%source(8)%identifier, "Recombination" )
        !! Charge exchange
        call set_edge_source_identifier( edge_sources%source(9)%identifier, "Charge exchange" )
        !! Collisional equipartition
        call set_edge_source_identifier( edge_sources%source(10)%identifier, "Equipartition" )
        !! Ohmic
        call set_edge_source_identifier( edge_sources%source(11)%identifier, "Ohmic" )
        !! Radiation
        call set_edge_source_identifier( edge_sources%source(12)%identifier, "Radiation" )
#ifdef B25_EIRENE
        !! Neutrals
        call set_edge_source_identifier( edge_sources%source(13)%identifier, "Neutrals" )
#endif
#elif ( IMAS_MAJOR_VERSION > 3 || IMAS_MINOR_VERSION > 30 )
        !! Total sources
        call write_source_identifier( &
          &  edge_sources%source(1)%identifier, edge_source_identifier%total )
        !! Background sources
        call write_source_identifier( &
          &  edge_sources%source(2)%identifier, edge_source_identifier%background )
        !! Prescribed sources
        call write_source_identifier( &
          &  edge_sources%source(3)%identifier, edge_source_identifier%prescribed )
        !! Time derivative
        call write_source_identifier( &
          &  edge_sources%source(4)%identifier, edge_source_identifier%time_derivative )
        !! Atomic ionization
        call write_source_identifier( &
          &  edge_sources%source(5)%identifier, edge_source_identifier%atomic_ionization )
        !! Molecular ionization
        call write_source_identifier( &
          &  edge_sources%source(6)%identifier, edge_source_identifier%molecular_ionization )
        !! Ionization
        call write_source_identifier( &
          &  edge_sources%source(7)%identifier, edge_source_identifier%ionization )
        !! Recombination
        call write_source_identifier( &
          &  edge_sources%source(8)%identifier, edge_source_identifier%recombination )
        !! Charge exchange
        call write_source_identifier( &
          &  edge_sources%source(9)%identifier, edge_source_identifier%charge_exchange )
        !! Collisional equipartition
        call write_source_identifier( &
          &  edge_sources%source(10)%identifier, edge_source_identifier%collisional_equipartition )
        !! Ohmic
        call write_source_identifier( &
          &  edge_sources%source(11)%identifier, edge_source_identifier%ohmic )
        !! Radiation
        call write_source_identifier( &
          &  edge_sources%source(12)%identifier, edge_source_identifier%radiation )
#if ( ( IMAS_MINOR_VERSION > 38 || IMAS_MAJOR_VERSION > 3 ) && defined(B25_EIRENE) )
        !! Neutrals
        call write_source_identifier( &
          &  edge_sources%source(13)%identifier, edge_source_identifier%neutrals )
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
#if IMAS_MAJOR_VERSION > 3
        allocate( plasma_sources%source(nsources) )
        do is = 1, nsources
          allocate( plasma_sources%source(is)%ggd( num_time_slices ) )
        end do
#if ( IMAS_MAJOR_VERSION > 4 || ( IMAS_MAJOR_VERSION == 4 && IMAS_MINOR_VERSION > 0 ) )
        !! Total sources
        call set_plasma_source_identifier( plasma_sources%source(1)%identifier, "Total" )
        !! Background sources
        call set_plasma_source_identifier( plasma_sources%source(2)%identifier, "Background" )
        !! Prescribed sources
        call set_plasma_source_identifier( plasma_sources%source(3)%identifier, "Prescribed" )
        !! Time derivative
        call set_plasma_source_identifier( plasma_sources%source(4)%identifier, "Time derivative" )
        !! Atomic ionization
        call set_plasma_source_identifier( plasma_sources%source(5)%identifier, "Atomic ionization" )
        !! Molecular ionization
        call set_plasma_source_identifier( plasma_sources%source(6)%identifier, "Molecular ionization" )
        !! Ionization
        call set_plasma_source_identifier( plasma_sources%source(7)%identifier, "Ionization" )
        !! Recombination
        call set_plasma_source_identifier( plasma_sources%source(8)%identifier, "Recombination" )
        !! Charge exchange
        call set_plasma_source_identifier( plasma_sources%source(9)%identifier, "Charge exchange" )
        !! Collisional equipartition
        call set_plasma_source_identifier( plasma_sources%source(10)%identifier, "Equipartition" )
        !! Ohmic
        call set_plasma_source_identifier( plasma_sources%source(11)%identifier, "Ohmic" )
        !! Radiation
        call set_plasma_source_identifier( plasma_sources%source(12)%identifier, "Radiation" )
#ifdef B25_EIRENE
        !! Neutrals
        call set_plasma_source_identifier( plasma_sources%source(13)%identifier, "Neutrals" )
#endif
#else
        !! Total sources
        call write_source_identifier( &
          &  plasma_sources%source(1)%identifier, edge_source_identifier%total )
        !! Background sources
        call write_source_identifier( &
          &  plasma_sources%source(2)%identifier, edge_source_identifier%background )
        !! Prescribed sources
        call write_source_identifier( &
          &  plasma_sources%source(3)%identifier, edge_source_identifier%prescribed )
        !! Time derivative
        call write_source_identifier( &
          &  plasma_sources%source(4)%identifier, edge_source_identifier%time_derivative )
        !! Atomic ionization
        call write_source_identifier( &
          &  plasma_sources%source(5)%identifier, edge_source_identifier%atomic_ionization )
        !! Molecular ionization
        call write_source_identifier( &
          &  plasma_sources%source(6)%identifier, edge_source_identifier%molecular_ionization )
        !! Ionization
        call write_source_identifier( &
          &  plasma_sources%source(7)%identifier, edge_source_identifier%ionization )
        !! Recombination
        call write_source_identifier( &
          &  plasma_sources%source(8)%identifier, edge_source_identifier%recombination )
        !! Charge exchange
        call write_source_identifier( &
          &  plasma_sources%source(9)%identifier, edge_source_identifier%charge_exchange )
        !! Collisional equipartition
        call write_source_identifier( &
          &  plasma_sources%source(10)%identifier, edge_source_identifier%collisional_equipartition )
        !! Ohmic
        call write_source_identifier( &
          &  plasma_sources%source(11)%identifier, edge_source_identifier%ohmic )
        !! Radiation
        call write_source_identifier( &
          &  plasma_sources%source(12)%identifier, edge_source_identifier%radiation )
#ifdef B25_EIRENE
        !! Neutrals
        call write_source_identifier( &
          &  plasma_sources%source(13)%identifier, edge_source_identifier%neutrals )
#endif
#endif
#endif

        call put_equilibrium_data ( equilibrium, &
#if ( IMAS_MINOR_VERSION > 21 || IMAS_MAJOR_VERSION > 3 )
            &  summary, &
#endif
            &  edge_profiles, database, &
#if ( AL_MAJOR_VERSION > 4 && GGD_MAJOR_VERSION > 0 )
            &  time_sind, &
#endif
            &  .true., new_eq_ggd )
#if IMAS_MAJOR_VERSION > 3
        allocate( plasma_profiles%vacuum_toroidal_field%b0( num_time_slices ) )
        plasma_profiles%vacuum_toroidal_field%b0( time_sind ) = &
            &  edge_profiles%vacuum_toroidal_field%b0( time_sind )
        plasma_profiles%vacuum_toroidal_field%r0 = &
            &  edge_profiles%vacuum_toroidal_field%r0
#endif
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
#if ( IMAS_MAJOR_VERSION < 4 || ( IMAS_MAJOR_VERSION == 4 && IMAS_MINOR_VERSION < 1 ) )
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
#else
        if ( present( time_step_IN ) ) &
          &  summary%simulation%time_step = time_step_IN
        if ( present ( time_IN ) ) &
          &  summary%simulation%time_current = time_IN
        allocate( summary%simulation%workflow(1) )
        summary%simulation%workflow = source
        summary%simulation%time_begin = run_start_time_IN
        summary%simulation%time_end = run_end_time_IN
        summary%type%index = 2
        allocate( summary%type%name(1) )
        allocate( summary%type%description(1) )
        summary%type%name = "simulation"
        summary%type%description = "Simulation results from "//trim(source)
        allocate( summary%machine(1) )
        summary%machine = database
        summary%pulse = shot
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
#if IMAS_MAJOR_VERSION > 3
        call write_ids_midplane( plasma_profiles%midplane, midplane_id )
        call write_ids_midplane( plasma_sources%midplane, midplane_id )
        call write_ids_midplane( plasma_transport%midplane, midplane_id )
#endif
        call write_ids_midplane( summary%midplane, midplane_id )
#endif

#if ( IMAS_MAJOR_VERSION > 4 || ( IMAS_MAJOR_VERSION == 4 && IMAS_MINOR_VERSION > 0 ) )
        nesum = 0.0_IDS_real
        do iy = -1, ny
          do ix = -1, nx
            if (.not.on_closed_surface(ix,iy)) cycle
            nesum = nesum + ne(ix,iy)*vol(ix,iy)
          end do
        end do
        do i = 1, nspecies
          frac = 0.0_IDS_real
          do is = eb2spcr(i), eb2spcr(i)+nfluids(i)+1
            if (is.ge.ns) cycle
            if (is_neutral(is)) cycle
            if (.not.(is.eq.eb2spcr(i).or.lnext(eb2spcr(i),is))) cycle
            do iy = -1, ny
              do ix = -1, nx
                if (.not.on_closed_surface(ix,iy)) cycle
                frac = frac + na(ix,iy,is)*vol(ix,iy)
              end do
            end do
          end do
          if (nesum.gt.0.0_IDS_real) frac = frac / nesum
          select case (is_codes(eb2spcr(is)))
          case ('H')
            call write_sourced_constant_2( summary%composition%hydrogen, frac )
          case ('D')
            call write_sourced_constant_2( summary%composition%deuterium, frac )
          case ('T')
            call write_sourced_constant_2( summary%composition%tritium, frac )
          case ('DT')
            call write_sourced_constant_2( summary%composition%deuterium_tritium, frac )
          case ('He')
            if (nint(am(eb2spcr(is))).eq.3) then
              call write_sourced_constant_2( summary%composition%helium_3, frac )
            else if (nint(am(eb2spcr(is))).eq.4) then
              call write_sourced_constant_2( summary%composition%helium_4, frac )
            end if
          case ('Li')
            call write_sourced_constant_2( summary%composition%lithium, frac )
          case ('Be')
            call write_sourced_constant_2( summary%composition%beryllium, frac )
          case ('B')
            call write_sourced_constant_2( summary%composition%boron, frac )
          case ('C')
            call write_sourced_constant_2( summary%composition%carbon, frac )
          case ('N')
#if ( IMAS_MAJOR_VERSION == 4 && IMAS_MINOR_VERSION == 1 )
            call write_sourced_constant( summary%composition%nitrogen, frac )
#else
            call write_sourced_constant_2( summary%composition%nitrogen, frac )
#endif
          case ('O')
            call write_sourced_constant_2( summary%composition%oxygen, frac )
          case ('Ne')
            call write_sourced_constant_2( summary%composition%neon, frac )
          case ('Ar')
            call write_sourced_constant_2( summary%composition%argon, frac )
          case ('Fe')
            call write_sourced_constant_2( summary%composition%iron, frac )
          case ('Xe')
            call write_sourced_constant_2( summary%composition%xenon, frac )
          case ('W')
            call write_sourced_constant_2( summary%composition%tungsten, frac )
          case ('Kr')
            call write_sourced_constant_2( summary%composition%krypton, frac )
          end select
        end do
#endif

        select case (GeometryType)
        case ( GEOMETRY_CYLINDER, GEOMETRY_LIMITER, GEOMETRY_ANNULUS )
          icnt = 1
          isep(1) = 2
        case ( GEOMETRY_SN, GEOMETRY_STELLARATORISLAND,  &
        &      GEOMETRY_LFS_SNOWFLAKE_MINUS, GEOMETRY_LFS_SNOWFLAKE_PLUS )
          icnt = 1
          isep(1) = 4
        case ( GEOMETRY_CDN, GEOMETRY_DDN_BOTTOM, GEOMETRY_DDN_TOP )
          icnt = 2
          isep(1) = 4
          isep(2) = 11
        case default
          icnt = 0
        end select
        u = 0.0_IDS_real
        do i = 1, icnt
          do ix = -1, nx
            do iy = -1, ny
              if (region(ix,iy,2).eq.isep(i)) u = u + fht(ix,iy,1,1)
            end do
          end do
        end do
        if (u.ne.0.0_IDS_real) then
#if ( IMAS_MINOR_VERSION > 28 || IMAS_MAJOR_VERSION > 3 )
          call write_sourced_value( summary%global_quantities%power_loss, u )
#endif
        end if

        select case (GeometryType)
        case ( GEOMETRY_LINEAR, GEOMETRY_CYLINDER )
          u = 0.0_IDS_real
          frac = 0.0_IDS_real
          do ix = 0, nx-1
            do iy = 0, ny-1
              match_found = jsep.gt.-1 .and. jsep.le.ny
              match_found = match_found .and. iy.le.jsep
              do is = 0, ns-1
                u = u + rqrad(ix,iy,is) + rqbrm(ix,iy,is)
                if (match_found) frac = frac + rqrad(ix,iy,is) + rqbrm(ix,iy,is)
              end do
#ifdef B25_EIRENE
              do is = 1, natmi
                u = u - eneutrad(ix+1,iy+1,is,0)
                if (match_found) frac = frac - eneutrad(ix+1,iy+1,is,0)
              end do
              do is = 1, nmoli
                u = u - emolrad(ix+1,iy+1,is,0)
                if (match_found) frac = frac - emolrad(ix+1,iy+1,is,0)
              end do
              do is = 1, nioni
                u = u - eionrad(ix+1,iy+1,is,0)
                if (match_found) frac = frac - eionrad(ix+1,iy+1,is,0)
              end do
#endif
            end do
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
          do ix = 0, nx-1
            do iy = 0, ny-1
              if (on_closed_surface(ix,iy) .and. iy.le.jsep) cycle
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
#if ( IMAS_MINOR_VERSION > 30 || IMAS_MAJOR_VERSION > 3 )
          if (u.ne.0.0_IDS_real) then
            call write_sourced_value( summary%global_quantities%power_radiated_outside_lcfs, u )
          end if
#endif
        end select

#if ( IMAS_MINOR_VERSION > 36 || IMAS_MAJOR_VERSION > 3 )
        select case (GeometryType)
        case ( GEOMETRY_SN, GEOMETRY_CDN, GEOMETRY_DDN_BOTTOM, &
            &  GEOMETRY_LFS_SNOWFLAKE_MINUS, GEOMETRY_LFS_SNOWFLAKE_PLUS)
          call write_sourced_rz( summary%boundary%x_point_main, &
            &   crx(leftcut(1),topcut(1),0), cry(leftcut(1),topcut(1),0) )
        case ( GEOMETRY_DDN_TOP )
          call write_sourced_rz( summary%boundary%x_point_main, &
            &   crx(leftcut(2),topcut(2),0), cry(leftcut(2),topcut(2),0) )
        end select
        select case (GeometryType)
        case ( GEOMETRY_CDN )
          call write_sourced_value( summary%boundary%distance_inner_outer_separatrices, 0.0_IDS_real )
        case ( GEOMETRY_DDN_BOTTOM, GEOMETRY_DDN_TOP, &
             & GEOMETRY_LFS_SNOWFLAKE_MINUS )
          u = norm( (crx(nmdpl,topcut(2),0) + crx(nmdpl,topcut(2),1))/2.0_R8 - &
            &       (crx(nmdpl,topcut(1),0) + crx(nmdpl,topcut(1),1))/2.0_R8, &
            &       (cry(nmdpl,topcut(2),0) + cry(nmdpl,topcut(2),1))/2.0_R8 - &
            &       (cry(nmdpl,topcut(1),0) + cry(nmdpl,topcut(1),1))/2.0_R8 )
          call write_sourced_value( summary%boundary%distance_inner_outer_separatrices, u )
        end select
#endif

        totFace(-1:nx,-1:ny,0)=abs(fht(-1:nx,-1:ny,0,0))
        totFace(-1:nx,-1:ny,1)=abs(fht(-1:nx,-1:ny,1,1))
        call divide_by_poloidal_areas(nx,ny,totFace,tmpFace)
        call divide_by_contact_areas(nx,ny,totFace,hlpFace)
        call alloc_b2plot_wall_loading(nlim,nsgmx)
        allocate(wrdtrg(0:ny-1,ntrgsx,0:DEF_NATM))
        wrdtrg = 0.0_IDS_real
        call b2ptrdl(wrdtrg)
        do i = 1, 2*max(1,nncut)
          if (itrg(i).eq.0) cycle
          qmax = 1.0_IDS_real
          ixmax(itrg(i)) = ixmid(itrg(i))
          do ix = ixmid(itrg(i)), ixpos(itrg(i)), idir(itrg(i))
            qtot = 0.0_IDS_real
            do iy = iysep(itrg(i))+1, ny
              qtot = qtot + idir(itrg(i))*(fht(ix,iy,0,0)+fht(ix,iy,0,1))
            end do
            if (abs(qtot).gt.abs(qmax).and.qtot*qmax.gt.0.0_IDS_real) then
              qmax = qtot
              ixmax(itrg(i)) = ix
            end if
          end do
          wetted_area(itrg(i)) = 0.0_IDS_real
          u = maxval(tmpFace(ixmax(itrg(i)),iysep(itrg(i))+1:ny,0)) &
            &  *exp(-1.0_IDS_real)
          j = maxloc(tmpFace(ixmax(itrg(i)),iysep(itrg(i))+1:ny,0),dim=1)
          k = iysep(itrg(i))+j
          do iy = iysep(itrg(i))+1, ny
            if (k.gt.iysep(itrg(i))+j) cycle
            if (tmpFace(ixmax(itrg(i)),iy,0).lt.u) then
              k = max(k,iy)
              frac = (tmpFace(ixmax(itrg(i)),iy-1,0)-u)/ &
                   & (tmpFace(ixmax(itrg(i)),iy-1,0)-    &
                   &  tmpFace(ixmax(itrg(i)),iy,0))
              wetted_area(itrg(i)) = wetted_area(itrg(i)) + &
                   &  frac*gs(ifpos(itrg(i)),iy,0)
            else
              wetted_area(itrg(i)) = wetted_area(itrg(i)) + &
                   &  gs(ifpos(itrg(i)),iy,0)
            end if
          end do
          recycled_flux(itrg(i)) = 0.0_IDS_real
          power_neutrals(itrg(i)) = 0.0_IDS_real
          power_incident(itrg(i)) = 0.0_IDS_real
          power_radiated(itrg(i)) = 0.0_IDS_real
          power_conducted(itrg(i)) = 0.0_IDS_real
          power_convected(itrg(i)) = 0.0_IDS_real
          power_recomb_neutrals(itrg(i)) = 0.0_IDS_real
          do iy = 0, ny-1
            u = 0.0_R8
            do is = 0, ns-1
              u = u + idir(itrg(i))* &
                 & (ti(ixpos(itrg(i)),iy) + &
                 &  te(ixpos(itrg(i)),iy)*rza(ixpos(itrg(i)),iy,is))* &
                 & (1.5_R8*(fna_32(ifpos(itrg(i)),iy,0,0,is)+ &
                 &          fna_32(ifpos(itrg(i)),iy,0,1,is)) + &
                 &  2.5_R8*(fna_52(ifpos(itrg(i)),iy,0,0,is)+ &
                 &          fna_52(ifpos(itrg(i)),iy,0,1,is)))
              if (is_neutral(is)) &
                &  power_neutrals(itrg(i)) = power_neutrals(itrg(i)) + &
                &    idir(itrg(i))*(ti(ixpos(itrg(i)),iy)* &
                &   (1.5_R8*(fna_32(ifpos(itrg(i)),iy,0,0,is)+ &
                &            fna_32(ifpos(itrg(i)),iy,0,1,is)) + &
                &    2.5_R8*(fna_52(ifpos(itrg(i)),iy,0,0,is)+ &
                &            fna_52(ifpos(itrg(i)),iy,0,1,is))) + &
                &              (fhm(ifpos(itrg(i)),iy,0,0,is) + &
                &               fhm(ifpos(itrg(i)),iy,0,1,is)))
              match_found = .true.
              do j = 1, nstrai
                if (streql(crcstra(j),'E').or.streql(crcstra(j),'W')) then
                  if (rcpos(j).eq. &
                     &  (ixpos(itrg(i))+target_offset*idir(itrg(i))).and. &
                     &  (rcstart(j).le.iy .and. rcend(j).ge.iy)) then
                    recycled_flux(itrg(i)) = recycled_flux(itrg(i)) + &
                      &  neutral_sources_rescale*recyc(is,j)* &
                      &  max(0.0_R8,idir(itrg(i))* &
                      &  (fna(ifpos(itrg(i)),iy,0,0,is)+ &
                      &   fna(ifpos(itrg(i)),iy,0,1,is)))*zn(is)
                  end if
                end if
              end do
            end do
            power_radiated(itrg(i)) = power_radiated(itrg(i)) + &
                 &  wrdtrg(iy,itrg(i),0)
            power_convected(itrg(i)) = power_convected(itrg(i)) + u
            power_conducted(itrg(i)) = power_conducted(itrg(i))   &
                 & - u + idir(itrg(i))* &
                 & (fht(ifpos(itrg(i)),iy,0,0)-fhj(ifpos(itrg(i)),iy,0,0)+ &
                 &  fht(ifpos(itrg(i)),iy,0,1)-fhj(ifpos(itrg(i)),iy,0,1))
            power_incident(itrg(i)) = power_incident(itrg(i)) + &
                 &  wrdtrg(iy,itrg(i),0) + &
                 &  idir(itrg(i))*(fht(ifpos(itrg(i)),iy,0,0) + &
                 &                 fht(ifpos(itrg(i)),iy,0,1))
          end do
#ifdef B25_EIRENE
          do js = 1, nsts
            if (eirdiag_nds_typ(js).ne.2) cycle
            if (eirdiag_nds_srf(js).ne.ifpos(itrg(i))) cycle
            do iy = 0, ny-1
              if (eirdiag_nds_start(js).gt.(iy+1)) cycle
              if (eirdiag_nds_end(js).lt.(iy+1)) cycle
              ind = eirdiag_nds_ind(js)
              ias = ind+1-(eirdiag_nds_start(js)-1)
              power_neutrals(itrg(i)) = power_neutrals(itrg(i)) + &
                   &  ewldt_res(iy+ias)
              power_incident(itrg(i)) = power_incident(itrg(i)) + &
                   &  ewldt_res(iy+ias)
              do imol = 1, nmoli
                power_recomb_neutrals(itrg(i)) = &
                   &  power_recomb_neutrals(itrg(i)) + ewldmr_res(imol,iy+ias)
              end do
              hlpFace(ifpos(itrg(i)),iy,0) = hlpFace(ifpos(itrg(i)),iy,0) + &
                &  ewldt_res(iy+ias)/gs(ifpos(itrg(i)),iy,0)
            end do
            do iatm = 1, natmi
              recycled_flux(itrg(i)) = recycled_flux(itrg(i)) + &
                &  wldpa(nlim+js,iatm,0)*zn(eb2atcr(iatm))
            end do
            do imol = 1, nmoli
              do iatm = 1, natmi
                recycled_flux(itrg(i)) = recycled_flux(itrg(i)) + &
                  &  wldpm(nlim+js,imol,0)*mlcmp(iatm,imol)*zn(eb2atcr(iatm))
              end do
            end do
          end do
#endif
          power_flux_peak(itrg(i)) = maxval(hlpFace(ifpos(itrg(i)),0:ny-1,0))
        end do
#if ( IMAS_MINOR_VERSION > 30 || IMAS_MAJOR_VERSION > 3 )
        select case ( GeometryType )
        case ( GEOMETRY_LINEAR, GEOMETRY_CYLINDER )
          if (ntrgts.gt.0) then
            allocate( divertors%divertor(ntrgts) )
            do i = 1, ntrgts
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
!! FIXME: Should represent the full extent of the physical divertor
              divertors%divertor(i)%target(1)%extension_r = extension_r(i)
              divertors%divertor(i)%target(1)%extension_z = extension_z(i)
            end do
            call write_timed_value( &
              &  divertors%divertor(1)%target(1)%power_flux_peak, &
              &  power_flux_peak(1) )
            call write_timed_value( &
              &  divertors%divertor(1)%target(1)%flux_expansion, &
              &  flux_expansion(1) )
            call write_timed_value( &
              &  divertors%divertor(1)%target(1)%wetted_area, &
              &  wetted_area(1) )
            call write_timed_value( &
              &  divertors%divertor(1)%wetted_area, &
              &  wetted_area(1) )
            call write_timed_value( &
              &  divertors%divertor(1)%target(1)%power_incident_fraction, &
              &  1.0_IDS_real )
            call write_timed_value( &
              &  divertors%divertor(1)%target(1)%power_incident, &
              &  power_incident(1) )
            call write_timed_value( &
              &  divertors%divertor(1)%power_incident, &
              &  power_incident(1) )
            call write_timed_value( &
              &  divertors%divertor(1)%target(1)%power_conducted, &
              &  power_conducted(1) )
            call write_timed_value( &
              &  divertors%divertor(1)%power_conducted, &
              &  power_conducted(1) )
            call write_timed_value( &
              &  divertors%divertor(1)%target(1)%power_convected, &
              &  power_convected(1) )
            call write_timed_value( &
              &  divertors%divertor(1)%power_convected, &
              &  power_convected(1) )
            call write_timed_value( &
              &  divertors%divertor(1)%target(1)%power_radiated, &
              &  power_radiated(1) )
            call write_timed_value( &
              &  divertors%divertor(1)%power_radiated, &
              &  power_radiated(1) )
            call write_timed_value( &
              &  divertors%divertor(1)%target(1)%power_neutrals, &
              &  power_neutrals(1) )
            call write_timed_value( &
              &  divertors%divertor(1)%power_neutrals, &
              &  power_neutrals(1) )
            u = idir(1)*sum(fhp(ifpos(1),:,0,:,:))
            call write_timed_value( &
              &  divertors%divertor(1)%target(1)%power_recombination_plasma, u )
            call write_timed_value( &
              &  divertors%divertor(1)%power_recombination_plasma, u )
            call write_timed_value( &
              &  divertors%divertor(1)%target(1)%power_recombination_neutrals, &
              &  power_recomb_neutrals(1) )
            call write_timed_value( &
              &  divertors%divertor(1)%power_recombination_neutrals, &
              &  power_recomb_neutrals(1) )
            u = idir(1)*sum(fhj(ifpos(1),:,0,:))
            call write_timed_value( &
              &  divertors%divertor(1)%target(1)%power_currents, u )
            call write_timed_value( &
              &  divertors%divertor(1)%power_currents, u )
#if ( IMAS_MINOR_VERSION > 32 || IMAS_MAJOR_VERSION > 3 )
            u = idir(1)*sum(fch(ifpos(1),:,0,:))
            call write_timed_value( &
              &  divertors%divertor(1)%target(1)%current_incident, u )
            call write_timed_value( &
              &  divertors%divertor(1)%current_incident, u )
#endif
            call write_timed_value( &
              &  divertors%divertor(1)%particle_flux_recycled_total, &
              &  recycled_flux(1) )
            if (ntrgts.eq.2) then
              call write_timed_value( &
                &  divertors%divertor(2)%target(1)%power_flux_peak, &
                &  power_flux_peak(2) )
              call write_timed_value( &
                &  divertors%divertor(2)%target(1)%flux_expansion, &
                &  flux_expansion(2) )
              call write_timed_value( &
                &  divertors%divertor(2)%target(1)%wetted_area, &
                &  wetted_area(2) )
              call write_timed_value( &
                &  divertors%divertor(2)%wetted_area, &
                &  wetted_area(2) )
              call write_timed_value( &
                &  divertors%divertor(2)%target(1)%power_incident_fraction, &
                &  1.0_IDS_real )
              call write_timed_value( &
                &  divertors%divertor(2)%target(1)%power_incident, &
                &  power_incident(2) )
              call write_timed_value( &
                &  divertors%divertor(2)%power_incident, &
                &  power_incident(2) )
              call write_timed_value( &
                &  divertors%divertor(2)%target(1)%power_conducted, &
                &  power_conducted(2) )
              call write_timed_value( &
                &  divertors%divertor(2)%power_conducted, &
                &  power_conducted(2) )
              call write_timed_value( &
                &  divertors%divertor(2)%target(1)%power_convected, &
                &  power_convected(2) )
              call write_timed_value( &
                &  divertors%divertor(2)%power_convected, &
                &  power_convected(2) )
              call write_timed_value( &
                &  divertors%divertor(2)%target(1)%power_radiated, &
                &  power_radiated(2) )
              call write_timed_value( &
                &  divertors%divertor(2)%power_radiated, &
                &  power_radiated(2) )
              call write_timed_value( &
                &  divertors%divertor(2)%target(1)%power_neutrals, &
                &  power_neutrals(2) )
              call write_timed_value( &
                &  divertors%divertor(2)%power_neutrals, &
                &  power_neutrals(2) )
              u = idir(2)*sum(fhp(ifpos(2),:,0,:,:))
              call write_timed_value( &
                &  divertors%divertor(2)%target(1)%power_recombination_plasma, u )
              call write_timed_value( &
                &  divertors%divertor(2)%power_recombination_plasma, u )
              call write_timed_value( &
                &  divertors%divertor(2)%target(1)%power_recombination_neutrals, &
                &  power_recomb_neutrals(2) )
              call write_timed_value( &
                &  divertors%divertor(2)%power_recombination_neutrals, &
                &  power_recomb_neutrals(2) )
              u = idir(2)*sum(fhj(ifpos(2),:,0,:))
              call write_timed_value( &
                &  divertors%divertor(2)%target(1)%power_currents, u )
              call write_timed_value( &
                &  divertors%divertor(2)%power_currents, u )
#if ( IMAS_MINOR_VERSION > 32 || IMAS_MAJOR_VERSION > 3 )
              u = idir(2)*sum(fch(ifpos(2),:,0,:))
              call write_timed_value( &
                &  divertors%divertor(2)%target(1)%current_incident, u )
              call write_timed_value( &
                &  divertors%divertor(2)%current_incident, u )
#endif
              call write_timed_value( &
                &  divertors%divertor(2)%particle_flux_recycled_total, &
                &  recycled_flux(2) )
            end if
          end if
        case ( GEOMETRY_SN, GEOMETRY_STELLARATORISLAND )
          allocate( divertors%divertor(1) )
          allocate( divertors%divertor(1)%name(1) )
          if (LSN) then
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
          if (LSN) then
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
          if (LSN) then
            divertors%divertor(1)%identifier = 'LSN'
            divertors%divertor(1)%target(1)%identifier = "ID"
            divertors%divertor(1)%target(2)%identifier = "OD"
          else
            divertors%divertor(1)%identifier = 'USN'
            divertors%divertor(1)%target(1)%identifier = "OD"
            divertors%divertor(1)%target(2)%identifier = "ID"
          end if
#endif
!! FIXME: Should represent the full extent of the physical divertor
          divertors%divertor(1)%target(1)%extension_r = extension_r(1)
          divertors%divertor(1)%target(1)%extension_z = extension_z(1)
          divertors%divertor(1)%target(2)%extension_r = extension_r(2)
          divertors%divertor(1)%target(2)%extension_z = extension_z(2)
          if (LSN) then
            divertors%divertor(1)%target(1)%name = "Inner target"
            divertors%divertor(1)%target(2)%name = "Outer target"
          else
            divertors%divertor(1)%target(1)%name = "Outer target"
            divertors%divertor(1)%target(2)%name = "Inner target"
          end if
          call write_timed_value( &
            &  divertors%divertor(1)%target(1)%power_flux_peak, &
            &  power_flux_peak(1) )
          call write_timed_value( &
            &  divertors%divertor(1)%target(2)%power_flux_peak, &
            &  power_flux_peak(2) )
          call write_timed_value( &
            &  divertors%divertor(1)%target(1)%flux_expansion, &
            &  flux_expansion(1) )
          call write_timed_value( &
            &  divertors%divertor(1)%target(2)%flux_expansion, &
            &  flux_expansion(2) )
          call write_timed_value( &
            &  divertors%divertor(1)%target(1)%wetted_area, &
            &  wetted_area(1) )
          call write_timed_value( &
            &  divertors%divertor(1)%target(2)%wetted_area, &
            &  wetted_area(2) )
          call write_timed_value( &
            &  divertors%divertor(1)%wetted_area, &
            &  wetted_area(1)+wetted_area(2) )
          u = power_incident(1) / (power_incident(1) + power_incident(2))
          call write_timed_value( &
            &  divertors%divertor(1)%target(1)%power_incident_fraction, u )
          u = power_incident(2) / (power_incident(1) + power_incident(2))
          call write_timed_value( &
            &  divertors%divertor(1)%target(2)%power_incident_fraction, u )
          call write_timed_value( &
            &  divertors%divertor(1)%target(1)%power_incident, &
            &  power_incident(1) )
          call write_timed_value( &
            &  divertors%divertor(1)%target(2)%power_incident, &
            &  power_incident(2) )
          call write_timed_value( &
            &  divertors%divertor(1)%power_incident, &
            &  power_incident(1)+power_incident(2) )
          call write_timed_value( &
            &  divertors%divertor(1)%target(1)%power_conducted, &
            &  power_conducted(1) )
          call write_timed_value( &
            &  divertors%divertor(1)%target(2)%power_conducted, &
            &  power_conducted(2) )
          call write_timed_value( &
            &  divertors%divertor(1)%power_conducted, &
            &  power_conducted(1)+power_conducted(2) )
          call write_timed_value( &
            &  divertors%divertor(1)%target(1)%power_convected, &
            &  power_convected(1) )
          call write_timed_value( &
            &  divertors%divertor(1)%target(2)%power_convected, &
            &  power_convected(2) )
          call write_timed_value( &
            &  divertors%divertor(1)%power_convected, &
            &  power_convected(1)+power_convected(2) )
          call write_timed_value( &
            &  divertors%divertor(1)%target(1)%power_radiated, &
            &  power_radiated(1) )
          call write_timed_value( &
            &  divertors%divertor(1)%target(2)%power_radiated, &
            &  power_radiated(2) )
          call write_timed_value( &
            &  divertors%divertor(1)%power_radiated, &
            &  power_radiated(1)+power_radiated(2) )
          call write_timed_value( &
            &  divertors%divertor(1)%target(1)%power_neutrals, &
            &  power_neutrals(1) )
          call write_timed_value( &
            &  divertors%divertor(1)%target(2)%power_neutrals, &
            &  power_neutrals(2) )
          call write_timed_value( &
            &  divertors%divertor(1)%power_neutrals, &
            &  power_neutrals(1)+power_neutrals(2) )
          u = idir(1)*sum(fhp(ifpos(1),:,0,:,:))
          call write_timed_value( &
            &  divertors%divertor(1)%target(1)%power_recombination_plasma, u )
          v = idir(2)*sum(fhp(ifpos(2),:,0,:,:))
          call write_timed_value( &
            &  divertors%divertor(1)%target(2)%power_recombination_plasma, v )
          call write_timed_value( &
            &  divertors%divertor(1)%power_recombination_plasma, u+v )
          call write_timed_value( &
            &  divertors%divertor(1)%target(1)%power_recombination_neutrals, &
            &  power_recomb_neutrals(1) )
          call write_timed_value( &
            &  divertors%divertor(1)%target(2)%power_recombination_neutrals, &
            &  power_recomb_neutrals(2) )
          call write_timed_value( &
            &  divertors%divertor(1)%power_recombination_neutrals, &
            &  power_recomb_neutrals(1)+power_recomb_neutrals(2) )
          u = idir(1)*sum(fhj(ifpos(1),:,0,:))
          call write_timed_value( &
            &  divertors%divertor(1)%target(1)%power_currents, u )
          v = idir(2)*sum(fhj(ifpos(2),:,0,:))
          call write_timed_value( &
            &  divertors%divertor(1)%target(2)%power_currents, v )
          call write_timed_value( &
            &  divertors%divertor(1)%power_currents, u+v )
#if ( IMAS_MINOR_VERSION > 32 || IMAS_MAJOR_VERSION > 3 )
          u = idir(1)*sum(fch(ifpos(1),:,0,:))
          call write_timed_value( &
            &  divertors%divertor(1)%target(1)%current_incident, u )
          v = idir(2)*sum(fch(ifpos(2),:,0,:))
          call write_timed_value( &
            &  divertors%divertor(1)%target(2)%current_incident, v )
          call write_timed_value( &
            &  divertors%divertor(1)%current_incident, u+v )
#endif
          call write_timed_value( &
            &  divertors%divertor(1)%particle_flux_recycled_total, &
            &  recycled_flux(1) )
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
!! FIXME: Should represent the full extent of the physical divertor
          divertors%divertor(1)%target(1)%extension_r = extension_r(1)
          divertors%divertor(1)%target(1)%extension_z = extension_z(1)
          divertors%divertor(1)%target(2)%extension_r = extension_r(4)
          divertors%divertor(1)%target(2)%extension_z = extension_z(4)
          divertors%divertor(2)%target(1)%extension_r = extension_r(2)
          divertors%divertor(2)%target(1)%extension_z = extension_z(2)
          divertors%divertor(2)%target(2)%extension_r = extension_r(3)
          divertors%divertor(2)%target(2)%extension_z = extension_z(3)
          call write_timed_value( &
            &  divertors%divertor(1)%target(1)%power_flux_peak, &
            &  power_flux_peak(1) )
          call write_timed_value( &
            &  divertors%divertor(1)%target(2)%power_flux_peak, &
            &  power_flux_peak(4) )
          call write_timed_value( &
            &  divertors%divertor(1)%target(1)%flux_expansion, &
            &  flux_expansion(1) )
          call write_timed_value( &
            &  divertors%divertor(1)%target(2)%flux_expansion, &
            &  flux_expansion(4) )
          call write_timed_value( &
            &  divertors%divertor(1)%target(1)%wetted_area, &
            &  wetted_area(1) )
          call write_timed_value( &
            &  divertors%divertor(1)%target(2)%wetted_area, &
            &  wetted_area(4) )
          call write_timed_value( &
            &  divertors%divertor(1)%wetted_area, &
            &  wetted_area(1)+wetted_area(4) )
          u = power_incident(1) / (power_incident(1) + power_incident(4))
          call write_timed_value( &
            &  divertors%divertor(1)%target(1)%power_incident_fraction, u )
          u = power_incident(4) / (power_incident(1) + power_incident(4))
          call write_timed_value( &
            &  divertors%divertor(1)%target(2)%power_incident_fraction, u )
          call write_timed_value( &
            &  divertors%divertor(1)%target(1)%power_incident, &
            &  power_incident(1) )
          call write_timed_value( &
            &  divertors%divertor(1)%target(2)%power_incident, &
            &  power_incident(4) )
          call write_timed_value( &
            &  divertors%divertor(1)%power_incident, &
            &  power_incident(1)+power_incident(4) )
          call write_timed_value( &
            &  divertors%divertor(1)%target(1)%power_conducted, &
            &  power_conducted(1) )
          call write_timed_value( &
            &  divertors%divertor(1)%target(2)%power_conducted, &
            &  power_conducted(4) )
          call write_timed_value( &
            &  divertors%divertor(1)%power_conducted, &
            &  power_conducted(1)+power_conducted(4) )
          call write_timed_value( &
            &  divertors%divertor(1)%target(1)%power_convected, &
            &  power_convected(1) )
          call write_timed_value( &
            &  divertors%divertor(1)%target(2)%power_convected, &
            &  power_convected(4) )
          call write_timed_value( &
            &  divertors%divertor(1)%power_convected, &
            &  power_convected(1)+power_convected(4) )
          call write_timed_value( &
            &  divertors%divertor(1)%target(1)%power_radiated, &
            &  power_radiated(1) )
          call write_timed_value( &
            &  divertors%divertor(1)%target(2)%power_radiated, &
            &  power_radiated(4) )
          call write_timed_value( &
            &  divertors%divertor(1)%power_radiated, &
            &  power_radiated(1)+power_radiated(4) )
          call write_timed_value( &
            &  divertors%divertor(1)%target(1)%power_neutrals, &
            &  power_neutrals(1) )
          call write_timed_value( &
            &  divertors%divertor(1)%target(2)%power_neutrals, &
            &  power_neutrals(4) )
          call write_timed_value( &
            &  divertors%divertor(1)%power_neutrals, &
            &  power_neutrals(1)+power_neutrals(4) )
          u = idir(1)*sum(fhp(ifpos(1),:,0,:,:))
          call write_timed_value( &
            &  divertors%divertor(1)%target(1)%power_recombination_plasma, u )
          v = idir(4)*sum(fhp(ifpos(4),:,0,:,:))
          call write_timed_value( &
            &  divertors%divertor(1)%target(2)%power_recombination_plasma, v )
          call write_timed_value( &
            &  divertors%divertor(1)%power_recombination_plasma, u+v )
          call write_timed_value( &
            &  divertors%divertor(1)%target(1)%power_recombination_neutrals, &
            &  power_recomb_neutrals(1) )
          call write_timed_value( &
            &  divertors%divertor(1)%target(2)%power_recombination_neutrals, &
            &  power_recomb_neutrals(4) )
          call write_timed_value( &
            &  divertors%divertor(1)%power_recombination_neutrals, &
            &  power_recomb_neutrals(1)+power_recomb_neutrals(4) )
          u = idir(1)*sum(fhj(ifpos(1),:,0,:))
          call write_timed_value( &
            &  divertors%divertor(1)%target(1)%power_currents, u )
          v = idir(4)*sum(fhj(ifpos(4),:,0,:))
          call write_timed_value( &
            &  divertors%divertor(1)%target(2)%power_currents, v )
          call write_timed_value( &
            &  divertors%divertor(1)%power_currents, u+v )
#if ( IMAS_MINOR_VERSION > 32 || IMAS_MAJOR_VERSION > 3 )
          u = idir(1)*sum(fch(ifpos(1),:,0,:))
          call write_timed_value( &
            &  divertors%divertor(1)%target(1)%current_incident, u )
          v = idir(4)*sum(fch(ifpos(4),:,0,:))
          call write_timed_value( &
            &  divertors%divertor(1)%target(2)%current_incident, v )
          call write_timed_value( &
            &  divertors%divertor(1)%current_incident, u+v )
#endif
          call write_timed_value( &
            &  divertors%divertor(1)%particle_flux_recycled_total, &
            &  recycled_flux(1)+recycled_flux(4) )
          call write_timed_value( &
            &  divertors%divertor(2)%target(1)%power_flux_peak, &
            &  power_flux_peak(2) )
          call write_timed_value( &
            &  divertors%divertor(2)%target(2)%power_flux_peak, &
            &  power_flux_peak(3) )
          call write_timed_value( &
            &  divertors%divertor(2)%target(1)%flux_expansion, &
            &  flux_expansion(2) )
          call write_timed_value( &
            &  divertors%divertor(2)%target(2)%flux_expansion, &
            &  flux_expansion(3) )
          call write_timed_value( &
            &  divertors%divertor(2)%target(1)%wetted_area, &
            &  wetted_area(2) )
          call write_timed_value( &
            &  divertors%divertor(2)%target(2)%wetted_area, &
            &  wetted_area(3) )
          call write_timed_value( &
            &  divertors%divertor(2)%wetted_area, &
            &  wetted_area(2)+wetted_area(3) )
          u = power_incident(2) / (power_incident(2) + power_incident(3))
          call write_timed_value( &
            &  divertors%divertor(2)%target(1)%power_incident_fraction, u )
          u = power_incident(3) / (power_incident(2) + power_incident(3))
          call write_timed_value( &
            &  divertors%divertor(2)%target(2)%power_incident_fraction, u )
          call write_timed_value( &
            &  divertors%divertor(2)%target(1)%power_incident, &
            &  power_incident(2) )
          call write_timed_value( &
            &  divertors%divertor(2)%target(2)%power_incident, &
            &  power_incident(3) )
          call write_timed_value( &
            &  divertors%divertor(2)%power_incident, &
            &  power_incident(2)+power_incident(3) )
          call write_timed_value( &
            &  divertors%divertor(2)%target(1)%power_conducted, &
            &  power_conducted(2) )
          call write_timed_value( &
            &  divertors%divertor(2)%target(2)%power_conducted, &
            &  power_conducted(3) )
          call write_timed_value( &
            &  divertors%divertor(2)%power_conducted, &
            &  power_conducted(2)+power_conducted(3) )
          call write_timed_value( &
            &  divertors%divertor(2)%target(1)%power_convected, &
            &  power_convected(2) )
          call write_timed_value( &
            &  divertors%divertor(2)%target(2)%power_convected, &
            &  power_convected(3) )
          call write_timed_value( &
            &  divertors%divertor(2)%power_convected, &
            &  power_convected(2)+power_convected(3) )
          call write_timed_value( &
            &  divertors%divertor(2)%target(1)%power_radiated, &
            &  power_radiated(2) )
          call write_timed_value( &
            &  divertors%divertor(2)%target(2)%power_radiated, &
            &  power_radiated(3) )
          call write_timed_value( &
            &  divertors%divertor(2)%power_radiated, &
            &  power_radiated(2)+power_radiated(3) )
          call write_timed_value( &
            &  divertors%divertor(2)%target(1)%power_neutrals, &
            &  power_neutrals(2) )
          call write_timed_value( &
            &  divertors%divertor(2)%target(2)%power_neutrals, &
            &  power_neutrals(3) )
          call write_timed_value( &
            &  divertors%divertor(2)%power_neutrals, &
            &  power_neutrals(2)+power_neutrals(3) )
          u = idir(2)*sum(fhp(ifpos(2),:,0,:,:))
          call write_timed_value( &
            &  divertors%divertor(2)%target(1)%power_recombination_plasma, u )
          v = idir(3)*sum(fhp(ifpos(3),:,0,:,:))
          call write_timed_value( &
            &  divertors%divertor(2)%target(2)%power_recombination_plasma, v )
          call write_timed_value( &
            &  divertors%divertor(2)%power_recombination_plasma, u+v )
          call write_timed_value( &
            &  divertors%divertor(2)%target(1)%power_recombination_neutrals, &
            &  power_recomb_neutrals(2) )
          call write_timed_value( &
            &  divertors%divertor(2)%target(2)%power_recombination_neutrals, &
            &  power_recomb_neutrals(3) )
          call write_timed_value( &
            &  divertors%divertor(2)%power_recombination_neutrals, &
            &  power_recomb_neutrals(2)+power_recomb_neutrals(3) )
          u = idir(2)*sum(fhj(ifpos(2),:,0,:))
          call write_timed_value( &
            &  divertors%divertor(2)%target(1)%power_currents, u )
          v = idir(3)*sum(fhj(ifpos(3),:,0,:))
          call write_timed_value( &
            &  divertors%divertor(2)%target(2)%power_currents, v )
          call write_timed_value( &
            &  divertors%divertor(2)%power_currents, u+v )
#if ( IMAS_MINOR_VERSION > 32 || IMAS_MAJOR_VERSION > 3 )
          u = idir(2)*sum(fch(ifpos(2),:,0,:))
          call write_timed_value( &
            &  divertors%divertor(2)%target(1)%current_incident, u )
          v = idir(3)*sum(fch(ifpos(3),:,0,:))
          call write_timed_value( &
            &  divertors%divertor(2)%target(2)%current_incident, v )
          call write_timed_value( &
            &  divertors%divertor(2)%current_incident, u+v )
#endif
          call write_timed_value( &
            &  divertors%divertor(2)%particle_flux_recycled_total, &
            &  recycled_flux(2)+recycled_flux(3) )
        end select
#endif
        deallocate(wrdtrg)
        call dealloc_b2plot_wall_loading
#endif

        !! Write grid & grid subsets/subgrids
#if ( IMAS_MINOR_VERSION > 11 || IMAS_MAJOR_VERSION > 3 ) && GGD_MAJOR_VERSION > 0
#if ( IMAS_MINOR_VERSION < 15 && IMAS_MAJOR_VERSION < 4 )
        call b2_IMAS_Fill_Grid_Desc( IDSmap,                                &
            &   edge_profiles%ggd( time_sind )%grid,                        &
            &   nx, ny, crx(-1:nx, -1:ny, :), cry(-1:nx, -1:ny, : ),        &
            &   leftix, leftiy, rightix, rightiy, topix, topiy, bottomix,   &
            &   bottomiy, nnreg, topcut, region, cflags,                    &
            &   INCLUDE_GHOST_CELLS, midplane_id, z_eq, vol, gs, qc )
        call b2_IMAS_Fill_Grid_Desc( IDSmap,                                &
            &   edge_transport%model(1)%ggd( time_sind )%grid,              &
            &   nx, ny, crx(-1:nx, -1:ny, :), cry(-1:nx, -1:ny, : ),        &
            &   leftix, leftiy, rightix, rightiy, topix, topiy, bottomix,   &
            &   bottomiy, nnreg, topcut, region, cflags,                    &
            &   INCLUDE_GHOST_CELLS, midplane_id, z_eq, vol, gs, qc )
        do is = 1, nsources
            call b2_IMAS_Fill_Grid_Desc( IDSmap,                              &
                &   edge_sources%source(is)%ggd( time_sind )%grid,            &
                &   nx, ny, crx(-1:nx, -1:ny, :), cry(-1:nx, -1:ny, : ),      &
                &   leftix, leftiy, rightix, rightiy, topix, topiy, bottomix, &
                &   bottomiy, nnreg, topcut, region, cflags,                  &
                &   INCLUDE_GHOST_CELLS, midplane_id, z_eq, vol, gs, qc )
        end do
#else
        call b2_IMAS_Fill_Grid_Desc( IDSmap,                                &
            &   edge_profiles%grid_ggd( time_sind ),                        &
            &   nx, ny, crx(-1:nx, -1:ny, :), cry(-1:nx, -1:ny, : ),        &
            &   leftix, leftiy, rightix, rightiy, topix, topiy, bottomix,   &
            &   bottomiy, nnreg, topcut, region, cflags,                    &
            &   INCLUDE_GHOST_CELLS, midplane_id, z_eq, vol, gs, qc )
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
#if IMAS_MAJOR_VERSION > 3
        allocate( plasma_profiles%grid_ggd( time_sind )%path(1) )
        plasma_profiles%grid_ggd( time_sind )%path = &
            &   "#edge_profiles/grid_ggd("//int2str(time_sind)//")"
        allocate( plasma_transport%grid_ggd( time_sind )%path(1) )
        plasma_transport%grid_ggd( time_sind )%path = &
            &   "#edge_profiles/grid_ggd("//int2str(time_sind)//")"
        allocate( plasma_sources%grid_ggd( time_sind )%path(1) )
        plasma_sources%grid_ggd( time_sind )%path = &
            &   "#edge_profiles/grid_ggd("//int2str(time_sind)//")"
#endif
#else
        call b2_IMAS_Fill_Grid_Desc( IDSmap,                                &
            &   edge_transport%grid_ggd( time_sind ),                       &
            &   nx, ny, crx(-1:nx, -1:ny, :), cry(-1:nx, -1:ny, : ),        &
            &   leftix, leftiy, rightix, rightiy, topix, topiy, bottomix,   &
            &   bottomiy, nnreg, topcut, region, cflags,                    &
            &   INCLUDE_GHOST_CELLS, midplane_id, z_eq, vol, gs, qc )
        call b2_IMAS_Fill_Grid_Desc( IDSmap,                                &
            &   edge_sources%grid_ggd( time_sind ),                         &
            &   nx, ny, crx(-1:nx, -1:ny, :), cry(-1:nx, -1:ny, : ),        &
            &   leftix, leftiy, rightix, rightiy, topix, topiy, bottomix,   &
            &   bottomiy, nnreg, topcut, region, cflags,                    &
            &   INCLUDE_GHOST_CELLS, midplane_id, z_eq, vol, gs, qc )
#if ( IMAS_MINOR_VERSION > 21 || IMAS_MAJOR_VERSION > 3 )
        call b2_IMAS_Fill_Grid_Desc( IDSmap,                                &
            &   radiation%grid_ggd( time_sind ),                            &
            &   nx, ny, crx(-1:nx, -1:ny, :), cry(-1:nx, -1:ny, : ),        &
            &   leftix, leftiy, rightix, rightiy, topix, topiy, bottomix,   &
            &   bottomiy, nnreg, topcut, region, cflags,                    &
            &   INCLUDE_GHOST_CELLS, midplane_id, z_eq, vol, gs, qc )
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
        edge_profiles%ggd( time_sind )%time = time_slice_value
        edge_transport%grid_ggd( time_sind )%time = time_slice_value
        edge_sources%grid_ggd( time_sind )%time = time_slice_value
#if IMAS_MAJOR_VERSION > 3
        plasma_profiles%grid_ggd( time_sind )%time = time_slice_value
        plasma_transport%grid_ggd( time_sind )%time = time_slice_value
        plasma_sources%grid_ggd( time_sind )%time = time_slice_value
#endif
#if ( IMAS_MINOR_VERSION > 21 || IMAS_MAJOR_VERSION > 3 )
        radiation%grid_ggd( time_sind )%time = time_slice_value
#endif
#endif

#if ( IMAS_MAJOR_VERSION < 4 || ( IMAS_MAJOR_VERSION == 4 && IMAS_MINOR_VERSION < 1 ) )
        edge_transport%model(1)%ggd( time_sind )%time = time_slice_value
        do is = 1, nsources
          edge_sources%source(is)%ggd( time_sind )%time = time_slice_value
        end do
#if IMAS_MAJOR_VERSION > 3
        plasma_profiles%ggd( time_sind )%time = time_slice_value
#if ( IMAS_MAJOR_VERSION == 4 && IMAS_MINOR_VERSION < 1 )
        plasma_transport%model(1)%ggd( time_sind )%time = time_slice_value
        do is = 1, nsources
          plasma_sources%source(is)%ggd( time_sind )%time = time_slice_value
        end do
#endif
#endif
#else
        allocate( sources_ggd( nsources ) )
        allocate( transport_ggd(1) )
        profiles_ggd%time = time_slice_value
        transport_ggd(1)%time = time_slice_value
        do is = 1, nsources
          sources_ggd(is)%time = time_slice_value
        end do
#endif
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
        if (use_eirene.ne.0) then
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
        if (use_eirene.ne.0) then
          do is = nspecies+1, nsion
            do js = 1, nioni
              if (imiion(js).eq.is) istion(is) = ionstt(ns-1+js)
            end do
          end do
        end if
#endif
#if ( IMAS_MAJOR_VERSION < 4 || ( IMAS_MAJOR_VERSION == 4 && IMAS_MINOR_VERSION < 1 ) )
        allocate( edge_profiles%ggd( time_sind )%ion( nsion ) )
        allocate( edge_transport%model(1)%ggd( time_sind )%ion( nsion ) )
        do i = 1, nsources
          allocate( edge_sources%source(i)%ggd( time_sind )%ion( nsion ) )
        end do
#else
        allocate( profiles_ggd%ion( nsion ) )
        allocate( transport_ggd(1)%ion( nsion ) )
        do i = 1, nsources
          allocate( sources_ggd(i)%ion( nsion ) )
        end do
#endif
#if ( IMAS_MINOR_VERSION > 21 || IMAS_MAJOR_VERSION > 3 )
        allocate( radiation%process(1)%ggd( time_sind )%ion( nsion ) )
        allocate( radiation%process(2)%ggd( time_sind )%ion( nsion ) )
#endif
        do js = 1, nspecies
#if ( IMAS_MAJOR_VERSION < 4 || ( IMAS_MAJOR_VERSION == 4 && IMAS_MINOR_VERSION < 1 ) )
          allocate( edge_profiles%ggd( time_sind )%ion( js )%element(1) )
          allocate( edge_profiles%ggd( time_sind )%ion( js )%state( nfluids(js) ) )
          allocate( edge_transport%model(1)%ggd( time_sind )%ion( js )%element(1) )
          allocate( edge_transport%model(1)%ggd( time_sind )%ion( js )%state( nfluids(js) ) )
#if ( IMAS_MAJOR_VERSION < 4 && IMAS_MINOR_VERSION < 42 )
          allocate( edge_profiles%ggd( time_sind )%ion( js )%label(1) )
          allocate( edge_transport%model(1)%ggd( time_sind )%ion( js )%label(1) )
          do is = 1, nfluids(js)
            allocate( edge_profiles%ggd( time_sind )%ion( js )%state( is )%label(1) )
            allocate( edge_transport%model(1)%ggd( time_sind )%ion( js )%state( is )%label(1) )
          end do
#else
          allocate( edge_profiles%ggd( time_sind )%ion( js )%name(1) )
          allocate( edge_transport%model(1)%ggd( time_sind )%ion( js )%name(1) )
          do is = 1, nfluids(js)
            allocate( edge_profiles%ggd( time_sind )%ion( js )%state( is )%name(1) )
            allocate( edge_transport%model(1)%ggd( time_sind )%ion( js )%state( is )%name(1) )
          end do
#endif
          do i = 1, nsources
            allocate( edge_sources%source(i)%ggd( time_sind )%ion( js )%state( nfluids(js) ) )
            allocate( edge_sources%source(i)%ggd( time_sind )%ion( js )%element(1) )
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
          end do
#else
          allocate( profiles_ggd%ion( js )%element(1) )
          allocate( profiles_ggd%ion( js )%state( nfluids(js) ) )
          allocate( profiles_ggd%ion( js )%name(1) )
          allocate( transport_ggd(1)%ion( js )%element(1) )
          allocate( transport_ggd(1)%ion( js )%state( nfluids(js) ) )
          allocate( transport_ggd(1)%ion( js )%name(1) )
          do is = 1, nfluids(js)
            allocate( profiles_ggd%ion( js )%state( is )%name(1) )
            allocate( transport_ggd(1)%ion( js )%state( is )%name(1) )
          end do
          do i = 1, nsources
            allocate( sources_ggd(i)%ion( js )%element(1) )
            allocate( sources_ggd(i)%ion( js )%state( nfluids(js) ) )
            allocate( sources_ggd(i)%ion( js )%name(1) )
            do is = 1, nfluids(js)
              allocate( sources_ggd(i)%ion( js )%state( is )%name(1) )
            end do
          end do
#endif

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
            call shrink_label(spclabel)
#if ( IMAS_MAJOR_VERSION < 4 || ( IMAS_MAJOR_VERSION == 4 && IMAS_MINOR_VERSION < 1 ) )
#if ( IMAS_MAJOR_VERSION < 4 && IMAS_MINOR_VERSION < 42 )
            edge_profiles%ggd( time_sind )%ion( js )%state( is )%label = spclabel
            edge_transport%model(1)%ggd( time_sind )%ion( js )%state( is )%label = spclabel
            do i = 1, nsources
              edge_sources%source(i)%ggd( time_sind )%ion( js )%state( is )%label = spclabel
            end do
#if ( IMAS_MINOR_VERSION > 21 || IMAS_MAJOR_VERSION > 3 )
            radiation%process(1)%ggd( time_sind )%ion( js )%state( is )%label = spclabel
            radiation%process(2)%ggd( time_sind )%ion( js )%state( is )%label = spclabel
#endif
#else
            edge_profiles%ggd( time_sind )%ion( js )%state( is )%name = spclabel
            edge_transport%model(1)%ggd( time_sind )%ion( js )%state( is )%name = spclabel
            do i = 1, nsources
              edge_sources%source(i)%ggd( time_sind )%ion( js )%state( is )%name = spclabel
            end do
            radiation%process(1)%ggd( time_sind )%ion( js )%state( is )%name = spclabel
            radiation%process(2)%ggd( time_sind )%ion( js )%state( is )%name = spclabel
#endif
#else
            profiles_ggd%ion( js )%state( is )%name = spclabel
            transport_ggd(1)%ion( js )%state( is )%name = spclabel
            do i = 1, nsources
              sources_ggd(i)%ion( js )%state( is )%name = spclabel
            end do
            radiation%process(1)%ggd( time_sind )%ion( js )%state( is )%name = spclabel
            radiation%process(2)%ggd( time_sind )%ion( js )%state( is )%name = spclabel
#endif
          end do

          ! Put ion label identifying the species
#if ( IMAS_MAJOR_VERSION < 4 || ( IMAS_MAJOR_VERSION == 4 && IMAS_MINOR_VERSION < 1 ) )
#if ( IMAS_MAJOR_VERSION < 4 && IMAS_MINOR_VERSION < 42 )
          edge_profiles%ggd( time_sind )%ion( js )%label = species_list( js )
          edge_transport%model(1)%ggd( time_sind )%ion( js )%label = species_list( js )
#else
          edge_profiles%ggd( time_sind )%ion( js )%name = species_list( js )
          edge_transport%model(1)%ggd( time_sind )%ion( js )%name = species_list( js )
#endif
#else
          profiles_ggd%ion( js )%name = species_list( js )
          transport_ggd(1)%ion( js )%name = species_list( js )
#endif

          ! Put ion charge if single ion in species
          is = ispion(js,1)
          if (istion(js).eq.1) then
            ion_charge_int = nint((zamin(is)+zamax(is))/2.0_R8)
#if ( IMAS_MAJOR_VERSION < 4 || ( IMAS_MAJOR_VERSION == 4 && IMAS_MINOR_VERSION < 1 ) )
            edge_profiles%ggd( time_sind )%ion( js )%z_ion = ion_charge_int
            edge_transport%model(1)%ggd( time_sind )%ion( js )%z_ion = ion_charge_int
#else
            profiles_ggd%ion( js )%z_ion = ion_charge_int
            transport_ggd(1)%ion( js )%z_ion = ion_charge_int
#endif
          end if

#if ( IMAS_MAJOR_VERSION < 4 || ( IMAS_MAJOR_VERSION == 4 && IMAS_MINOR_VERSION < 1 ) )
          ! Put element data
          call fill_neutral_element( is, js, &
            &  edge_profiles%ggd( time_sind )%ion( js )%element(1) )
          call fill_neutral_element( is, js, &
            &  edge_transport%model(1)%ggd( time_sind )%ion( js )%element(1) )

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
#else
          ! Put element data
          call fill_neutral_element( is, js, profiles_ggd%ion( js )%element(1) )
          call fill_neutral_element( is, js, transport_ggd(1)%ion( js )%element(1) )

          ! Put neutral index
          profiles_ggd%ion( js )%neutral_index = b2eatcr(is)
          transport_ggd(1)%ion( js )%neutral_index = b2eatcr(is)

          ! Put multiple states flag
          profiles_ggd%ion( js )%multiple_states_flag = 1
          transport_ggd(1)%ion( js )%multiple_states_flag = 1

          do is = 1, istion(js)
            ks = ispion(js,is)

            ! Put minimum Z of the charge state bundle
            ! (z_min = z_max = 0 for a neutral)
            profiles_ggd%ion( js )%state( is )%z_min = zamin( ks )
            transport_ggd(1)%ion( js )%state( is )%z_min = zamin( ks )

            ! Put maximum Z of the charge state bundle
            profiles_ggd%ion( js )%state( is )%z_max = zamax( ks )
            transport_ggd(1)%ion( js )%state( is )%z_max = zamax( ks )
          end do
#endif

          do i = 1, nsources
#if ( IMAS_MAJOR_VERSION < 4 || ( IMAS_MAJOR_VERSION == 4 && IMAS_MINOR_VERSION < 1 ) )
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

            ! Put element data
            call fill_neutral_element( is, js, &
              &  edge_sources%source(i)%ggd( time_sind )%ion( js )%element(1) )

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
#else
            ! Put ion label identifying the species
            sources_ggd(i)%ion( js )%name = species_list( js )

            ! Put ion charge if single ion in species
            is = ispion(js,1)
            if (istion(js).eq.1) then
              ion_charge_int = nint((zamin(is)+zamax(is))/2.0_R8)
              sources_ggd(i)%ion( js )%z_ion = ion_charge_int
            end if

            ! Put element data
            call fill_neutral_element( is, js, sources_ggd(i)%ion( js )%element(1) )

            ! Put neutral index
            sources_ggd(i)%ion( js )%neutral_index = b2eatcr(is)

            ! Put multiple states flag
            sources_ggd(i)%ion( js )%multiple_states_flag = 1

            do is = 1, istion(js)
              ks = ispion(js,is)

              ! Put minimum Z of the charge state bundle
              ! (z_min = z_max = 0 for a neutral)
              sources_ggd(i)%ion( js )%state( is )%z_min = zamin( ks )

              ! Put maximum Z of the charge state bundle
              sources_ggd(i)%ion( js )%state( is )%z_max = zamax( ks )
            end do
#endif
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

          ! Put element data
          call fill_neutral_element( is, js, &
            &  radiation%process(1)%ggd( time_sind )%ion( js )%element(1) )
          call fill_neutral_element( is, js, &
            &  radiation%process(2)%ggd( time_sind )%ion( js )%element(1) )

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

        if (use_eirene.ne.0) then
#ifdef B25_EIRENE
          do js = nspecies+1, nsion
#if ( IMAS_MAJOR_VERSION < 4 || ( IMAS_MAJOR_VERSION == 4 && IMAS_MINOR_VERSION < 1 ) )
            allocate( edge_profiles%ggd( time_sind )%ion( js )%state( istion(js) ) )
            allocate( edge_transport%model(1)%ggd( time_sind )%ion( js )%state( istion(js) ) )
#if ( IMAS_MAJOR_VERSION < 4 && IMAS_MINOR_VERSION < 42 )
            allocate( edge_profiles%ggd( time_sind )%ion( js )%label(1) )
            allocate( edge_transport%model(1)%ggd( time_sind )%ion( js )%label(1) )
            do ks = 1, istion(js)
              allocate( edge_profiles%ggd( time_sind )%ion( js )%state( ks )%label(1) )
              allocate( edge_transport%model(1)%ggd( time_sind )%ion( js )%state( ks )%label(1) )
            end do
#else
            allocate( edge_profiles%ggd( time_sind )%ion( js )%name(1) )
            allocate( edge_transport%model(1)%ggd( time_sind )%ion( js )%name(1) )
            do ks = 1, istion(js)
              allocate( edge_profiles%ggd( time_sind )%ion( js )%state( ks )%name(1) )
              allocate( edge_transport%model(1)%ggd( time_sind )%ion( js )%state( ks )%name(1) )
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
#else
            allocate( profiles_ggd%ion( js )%state( istion(js) ) )
            allocate( profiles_ggd%ion( js )%name(1) )
            allocate( transport_ggd(1)%ion( js )%state( istion(js) ) )
            allocate( transport_ggd(1)%ion( js )%name(1) )
            do ks = 1, istion(js)
              allocate( profiles_ggd%ion( js )%state( ks )%name(1) )
              allocate( transport_ggd(1)%ion( js )%state( ks )%name(1) )
            end do
            do i = 1, nsources
              allocate( sources_ggd(i)%ion( js )%state( istion(js) ) )
              allocate( sources_ggd(i)%ion( js )%name(1) )
              do ks = 1, istion(js)
                allocate( sources_ggd(i)%ion( js )%state( ks )%name(1) )
              end do
            end do
#endif
            do ks = 1, istion(js)
              is = ispion(js,ks)
              spclabel = trim(textin(is-1))
              call shrink_label(spclabel)
#if ( IMAS_MAJOR_VERSION < 4 || ( IMAS_MAJOR_VERSION == 4 && IMAS_MINOR_VERSION < 1 ) )
#if ( IMAS_MAJOR_VERSION < 4 && IMAS_MINOR_VERSION < 42 )
              edge_profiles%ggd( time_sind )%ion( js )%state( ks )%label = spclabel
              edge_transport%model(1)%ggd( time_sind )%ion( js )%state( ks )%label = &
                  &                                                        spclabel
              do i = 1, nsources
                edge_sources%source(i)%ggd( time_sind )%ion( js )%state( ks )%label = &
                    &                                                      spclabel
              end do
#else
              edge_profiles%ggd( time_sind )%ion( js )%state( ks )%name = spclabel
              edge_transport%model(1)%ggd( time_sind )%ion( js )%state( ks )%name = &
                  &                                                       spclabel
              do i = 1, nsources
                edge_sources%source(i)%ggd( time_sind )%ion( js )%state( ks )%name = &
                    &                                                     spclabel
              end do
#endif
              nelems = count ( micmp( 1:natmi, is ) > 0 )
              allocate( edge_profiles%ggd( time_sind )%ion( js )%element( nelems ) )
              allocate( edge_transport%model(1)%ggd( time_sind )%ion( js )%element( nelems ) )
              call fill_mol_ion_elements( nelems, is, &
                &  edge_profiles%ggd( time_sind )%ion( js )%element(:) )
              call fill_mol_ion_elements( nelems, is, &
                &  edge_transport%model(1)%ggd( time_sind )%ion( js )%element(:) )
              do i = 1, nsources
                allocate( edge_sources%source(i)%ggd( time_sind )%ion( js )%element( nelems ) )
                call fill_mol_ion_elements( nelems, is, &
                  &  edge_sources%source(i)%ggd( time_sind )%ion( js )%element(:) )
              end do
#else
              profiles_ggd%ion( js )%state( ks )%name = spclabel
              transport_ggd(1)%ion( js )%state( ks )%name = spclabel
              nelems = count ( micmp( 1:natmi, is ) > 0 )
              allocate( profiles_ggd%ion( js )%element( nelems ) )
              allocate( transport_ggd(1)%ion( js )%element( nelems ) )
              call fill_mol_ion_elements( nelems, is, &
                &  profiles_ggd%ion( js )%element(:) )
              call fill_mol_ion_elements( nelems, is, &
                &  transport_ggd(1)%ion( js )%element(:) )
              do i = 1, nsources
                sources_ggd(i)%ion( js )%state( ks )%name = spclabel
                allocate( sources_ggd(i)%ion( js )%element( nelems ) )
                call fill_mol_ion_elements( nelems, is, &
                  &  sources_ggd(i)%ion( js )%element(:) )
              end do
#endif
            end do
            is = ispion(js,1)
            spclabel = trim(textin(is-1))
            call shrink_label(spclabel)
            if (istion(js).eq.1) then
#if ( IMAS_MAJOR_VERSION < 4 || ( IMAS_MAJOR_VERSION == 4 && IMAS_MINOR_VERSION < 1 ) )
              edge_profiles%ggd( time_sind )%ion( js )%z_ion = nchrgi( is )
              edge_transport%model(1)%ggd( time_sind )%ion( js )%z_ion = nchrgi( is )
#if ( IMAS_MAJOR_VERSION < 4 && IMAS_MINOR_VERSION < 42 )
              edge_profiles%ggd( time_sind )%ion( js )%label = spclabel
              edge_transport%model(1)%ggd( time_sind )%ion( js )%label = spclabel
#else
              edge_profiles%ggd( time_sind )%ion( js )%name = spclabel
              edge_transport%model(1)%ggd( time_sind )%ion( js )%name = spclabel
#endif
              do i = 1, nsources
                edge_sources%source(i)%ggd( time_sind )%ion( js )%z_ion = nchrgi( is )
#if ( IMAS_MAJOR_VERSION < 4 && IMAS_MINOR_VERSION < 42 )
                edge_sources%source(i)%ggd( time_sind )%ion( js )%label = spclabel
#else
                edge_sources%source(i)%ggd( time_sind )%ion( js )%name = spclabel
#endif
              end do
#else
              profiles_ggd%ion( js )%z_ion = nchrgi( is )
              profiles_ggd%ion( js )%name = spclabel
              transport_ggd(1)%ion( js )%z_ion = nchrgi( is )
              transport_ggd(1)%ion( js )%name = spclabel
              do i = 1, nsources
                sources_ggd(i)%ion( js )%z_ion = nchrgi( is )
                sources_ggd(i)%ion( js )%name = spclabel
              end do
#endif
            else
              match_found = .false.
              do ks = 2, istion(js)
                match_found = match_found .and. nchrgi(is).eq.nchrgi(ispion(js,ks))
              end do
              if (match_found) then
#if ( IMAS_MAJOR_VERSION < 4 || ( IMAS_MAJOR_VERSION == 4 && IMAS_MINOR_VERSION < 1 ) )
                edge_profiles%ggd( time_sind )%ion( js )%z_ion = nchrgi( is )
                edge_transport%model(1)%ggd( time_sind )%ion( js )%z_ion = nchrgi( is )
                do i = 1, nsources
                  edge_sources%source(i)%ggd( time_sind )%ion( js )%z_ion = nchrgi( is )
                end do
#else
                profiles_ggd%ion( js )%z_ion = nchrgi( is )
                transport_ggd(1)%ion( js )%z_ion = nchrgi( is )
                do i = 1, nsources
                  sources_ggd(i)%ion( js )%z_ion = nchrgi( is )
                end do
#endif
              end if
              match_found = .false.
              p = index(textin(ispion(js,1)-1),'+')
              if (p.gt.1) then
                ion_label = textin(ispion(js,1)-1)(1:p-1)
              else
                ion_label = textin(ispion(js,1)-1)
              end if
              call shrink_label(ion_label)
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
#if ( IMAS_MAJOR_VERSION < 4 || ( IMAS_MAJOR_VERSION == 4 && IMAS_MINOR_VERSION < 1 ) )
#if ( IMAS_MAJOR_VERSION < 4 && IMAS_MINOR_VERSION < 42 )
                  edge_profiles%ggd( time_sind )%ion( js )%label = ion_label
                  edge_transport%model(1)%ggd( time_sind )%ion( js )%label = ion_label
                  do i = 1, nsources
                    edge_sources%source(i)%ggd( time_sind )%ion( js )%label = ion_label
                  end do
#else
                  edge_profiles%ggd( time_sind )%ion( js )%name = ion_label
                  edge_transport%model(1)%ggd( time_sind )%ion( js )%name = ion_label
                  do i = 1, nsources
                    edge_sources%source(i)%ggd( time_sind )%ion( js )%name = ion_label
                  end do
#endif
#else
                  profiles_ggd%ion( js )%name = ion_label
                  transport_ggd(1)%ion( js )%name = ion_label
                  do i = 1, nsources
                    sources_ggd(i)%ion( js )%name = ion_label
                  end do
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
#if ( IMAS_MAJOR_VERSION < 4 || ( IMAS_MAJOR_VERSION == 4 && IMAS_MINOR_VERSION < 1 ) )
            edge_profiles%ggd( time_sind )%ion( js )%neutral_index = jj
            edge_profiles%ggd( time_sind )%ion( js )%multiple_states_flag = 1
            edge_transport%model(1)%ggd( time_sind )%ion( js )%neutral_index = jj
            edge_transport%model(1)%ggd( time_sind )%ion( js )%multiple_states_flag = 1
            do ks = 1, istion(js)
              is = ispion(js,ks)
              edge_profiles%ggd( time_sind )%ion( js )%state( ks )%z_min = nchrgi( is )
              edge_profiles%ggd( time_sind )%ion( js )%state( ks )%z_max = nchrgi( is )
              edge_transport%model(1)%ggd( time_sind )%ion( js )%state( ks )%z_min = nchrgi( is )
              edge_transport%model(1)%ggd( time_sind )%ion( js )%state( ks )%z_max = nchrgi( is )
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
#else
            profiles_ggd%ion( js )%neutral_index = jj
            profiles_ggd%ion( js )%multiple_states_flag = 1
            transport_ggd(1)%ion( js )%neutral_index = jj
            transport_ggd(1)%ion( js )%multiple_states_flag = 1
            do ks = 1, istion(js)
              is = ispion(js,ks)
              profiles_ggd%ion( js )%state( ks )%z_min = nchrgi( is )
              profiles_ggd%ion( js )%state( ks )%z_max = nchrgi( is )
              transport_ggd(1)%ion( js )%state( ks )%z_min = nchrgi( is )
              transport_ggd(1)%ion( js )%state( ks )%z_max = nchrgi( is )
            end do
            do i = 1, nsources
              sources_ggd(i)%ion( js )%neutral_index = jj
              sources_ggd(i)%ion( js )%multiple_states_flag = 1
              do ks = 1, istion(js)
                is = ispion(js,ks)
                sources_ggd(i)%ion( js )%state( ks )%z_min = nchrgi( is )
                sources_ggd(i)%ion( js )%state( ks )%z_max = nchrgi( is )
              end do
            end do
#endif
          end do

          nneut = nspecies
          do j = 1, nmoli
             ks = 1
             do jj = 1, j-1
                if ( nmassm(j).eq.nmassm(jj) .and. ncharm(j).eq.ncharm(jj) .and. &
                  &  nprt(j).eq.nprt(jj) .and. lkindm(j).eq.lkindm(jj) ) then
                   ks = ks + 1
                end if
             end do
             if (ks.eq.1) nneut = nneut + 1
          end do
#if ( IMAS_MAJOR_VERSION < 4 || ( IMAS_MAJOR_VERSION == 4 && IMAS_MINOR_VERSION < 1 ) )
          allocate( edge_profiles%ggd( time_sind )%neutral( nneut ) )
          allocate( edge_transport%model(1)%ggd( time_sind )%neutral( nneut ) )
          do i = 1, nsources
             allocate( edge_sources%source(i)%ggd( time_sind )%neutral( nneut ) )
          end do
#else
          allocate( profiles_ggd%neutral( nneut ) )
          allocate( transport_ggd(1)%neutral( nneut ) )
          do i = 1, nsources
            allocate( sources_ggd(i)%neutral( nneut ) )
          end do
#endif
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
                if (js.gt.1) then
#if ( IMAS_MAJOR_VERSION < 4 || ( IMAS_MAJOR_VERSION == 4 && IMAS_MINOR_VERSION < 1 ) )
                  allocate( edge_profiles%ggd( time_sind )%neutral( js-1 )%state( isstat(is-1) ) )
                  allocate( edge_transport%model(1)%ggd( time_sind )%neutral( js-1 )%state( isstat(is-1) ) )
                  do i = 1, nsources
                    allocate( edge_sources%source(i)%ggd( time_sind )%neutral( js-1 )%state( isstat(is-1) ) )
                  end do
#else
                  allocate( profiles_ggd%neutral( js-1 )%state( isstat(is-1) ) )
                  allocate( transport_ggd(1)%neutral( js-1 )%state( isstat(is-1) ) )
                  do i = 1, nsources
                    allocate( sources_ggd(i)%neutral( js-1 )%state( isstat(is-1) ) )
                  end do
#endif
                end if
             end if
          end do
          ks = isstat(natmi)
          js = latmscl(natmi)
#if ( IMAS_MAJOR_VERSION < 4 || ( IMAS_MAJOR_VERSION == 4 && IMAS_MINOR_VERSION < 1 ) )
          allocate( edge_profiles%ggd( time_sind )%neutral( js )%state( ks ) )
          allocate( edge_transport%model(1)%ggd( time_sind )%neutral( js )%state( ks ) )
          do i = 1, nsources
            allocate( edge_sources%source(i)%ggd( time_sind )%neutral( js )%state( ks ) )
          end do
#else
          allocate( profiles_ggd%neutral( js )%state( ks ) )
          allocate( transport_ggd(1)%neutral( js )%state( ks ) )
          do i = 1, nsources
            allocate( sources_ggd(i)%neutral( js )%state( ks ) )
          end do
#endif

          do js = 1, nspecies
             is = eb2spcr(js)
#if ( IMAS_MAJOR_VERSION < 4 || ( IMAS_MAJOR_VERSION == 4 && IMAS_MINOR_VERSION < 1 ) )
             allocate( edge_profiles%ggd( time_sind )%neutral( js )%element(1) )
             allocate( edge_transport%model(1)%ggd( time_sind )%neutral( js )%element(1) )
             call fill_neutral_element( is, js, &
               &  edge_profiles%ggd( time_sind )%neutral( js )%element(1) )
             call fill_neutral_element( is, js, &
               &  edge_transport%model(1)%ggd( time_sind )%neutral( js )%element(1) )
#if ( IMAS_MAJOR_VERSION < 4 && IMAS_MINOR_VERSION < 42 )
             allocate( edge_profiles%ggd( time_sind )%neutral( js )%label(1) )
             allocate( edge_transport%model(1)%ggd( time_sind )%neutral( js )%label(1) )
             edge_profiles%ggd( time_sind )%neutral( js )%label = species_list( js )
             edge_transport%model(1)%ggd( time_sind )%neutral( js )%label = species_list( js )
#else
             allocate( edge_profiles%ggd( time_sind )%neutral( js )%name(1) )
             allocate( edge_transport%model(1)%ggd( time_sind )%neutral( js )%name(1) )
             edge_profiles%ggd( time_sind )%neutral( js )%name = species_list( js )
             edge_transport%model(1)%ggd( time_sind )%neutral( js )%name = species_list( js )
#endif
             edge_profiles%ggd( time_sind )%neutral( js )%ion_index = js
             edge_transport%model(1)%ggd( time_sind )%neutral( js )%ion_index = js
             do i = 1, nsources
                allocate( edge_sources%source(i)%ggd( time_sind )%neutral( js )%element(1) )
                call fill_neutral_element( is, js, &
                  &  edge_sources%source(i)%ggd( time_sind )%neutral( js )%element(1) )
#if ( IMAS_MAJOR_VERSION < 4 && IMAS_MINOR_VERSION < 42 )
                allocate( edge_sources%source(i)%ggd( time_sind )%neutral( js )%label(1) )
                edge_sources%source(i)%ggd( time_sind )%neutral( js )%label = species_list( js )
#else
                allocate( edge_sources%source(i)%ggd( time_sind )%neutral( js )%name(1) )
                edge_sources%source(i)%ggd( time_sind )%neutral( js )%name = species_list( js )
#endif
                edge_sources%source(i)%ggd( time_sind )%neutral( js )%ion_index = js
             end do
#else
             allocate( profiles_ggd%neutral( js )%element(1) )
             allocate( profiles_ggd%neutral( js )%name(1) )
             allocate( transport_ggd(1)%neutral( js )%element(1) )
             allocate( transport_ggd(1)%neutral( js )%name(1) )
             call fill_neutral_element( is, js, &
               &  profiles_ggd%neutral( js )%element(1) )
             call fill_neutral_element( is, js, &
               &  transport_ggd(1)%neutral( js )%element(1) )
             profiles_ggd%neutral( js )%name = species_list( js )
             profiles_ggd%neutral( js )%ion_index = js
             transport_ggd(1)%neutral( js )%name = species_list( js )
             transport_ggd(1)%neutral( js )%ion_index = js
             do i = 1, nsources
                allocate( sources_ggd(i)%neutral( js )%element(1) )
                allocate( sources_ggd(i)%neutral( js )%name(1) )
                call fill_neutral_element( is, js, &
                 &  sources_ggd(i)%neutral( js )%element(1) )
                sources_ggd(i)%neutral( js )%name = species_list( js )
                sources_ggd(i)%neutral( js )%ion_index = js
             end do
#endif
             do iss = 1, ks
                iatm = b2eatcr(is) + iss - 1
                spclabel = trim(textan(iatm-1))
                call shrink_label(spclabel)
#if ( IMAS_MAJOR_VERSION < 4 || ( IMAS_MAJOR_VERSION == 4 && IMAS_MINOR_VERSION < 1 ) )
#if ( IMAS_MAJOR_VERSION < 4 && IMAS_MINOR_VERSION < 42 )
                allocate( edge_profiles%ggd( time_sind )%neutral( js )%state( iss )%label(1) )
                allocate( edge_transport%model(1)%ggd( time_sind )%neutral( js )%state( iss )%label(1) )
                edge_profiles%ggd( time_sind )%neutral( js )%state( iss )%label = spclabel
                edge_transport%model(1)%ggd( time_sind )%neutral( js )%state( iss )%label = spclabel
#else
                allocate( edge_profiles%ggd( time_sind )%neutral( js )%state( iss )%name(1) )
                allocate( edge_transport%model(1)%ggd( time_sind )%neutral( js )%state( iss )%name(1) )
                edge_profiles%ggd( time_sind )%neutral( js )%state( iss )%name = spclabel
                edge_transport%model(1)%ggd( time_sind )%neutral( js )%state( iss )%name = spclabel
#endif
                call fill_atom_neutral_type( iatm, &
                  &  edge_profiles%ggd( time_sind )%neutral( js )%state( iss )%neutral_type )
                call fill_atom_neutral_type( iatm, &
                  &  edge_transport%model(1)%ggd( time_sind )%neutral( js )%state( iss )%neutral_type )
                edge_profiles%ggd( time_sind )%neutral( js )%multiple_states_flag = 1
                edge_transport%model(1)%ggd( time_sind )%neutral( js )%multiple_states_flag = 1
                do i = 1, nsources
#if ( IMAS_MAJOR_VERSION < 4 && IMAS_MINOR_VERSION < 42 )
                   allocate( edge_sources%source(i)%ggd( time_sind )%neutral( js )%state( iss )%label(1) )
                   edge_sources%source(i)%ggd( time_sind )%neutral( js )%state( iss )%label = spclabel
#else
                   allocate( edge_sources%source(i)%ggd( time_sind )%neutral( js )%state( iss )%name(1) )
                   edge_sources%source(i)%ggd( time_sind )%neutral( js )%state( iss )%name = spclabel
#endif
                   call fill_atom_neutral_type( iatm, &
                  &  edge_sources%source(i)%ggd( time_sind )%neutral( js )%state( iss )%neutral_type )
                   edge_sources%source(i)%ggd( time_sind )%neutral( js )%multiple_states_flag = 1
                end do
#else
                allocate( profiles_ggd%neutral( js )%state( iss )%name(1) )
                allocate( transport_ggd(1)%neutral( js )%state( iss )%name(1) )
                profiles_ggd%neutral( js )%state( iss )%name = spclabel
                transport_ggd(1)%neutral( js )%state( iss )%name = spclabel
                call fill_atom_neutral_type( iatm, &
                  &  profiles_ggd%neutral( js )%state( iss )%neutral_type )
                call fill_atom_neutral_type( iatm, &
                  &  transport_ggd(1)%neutral( js )%state( iss )%neutral_type )
                profiles_ggd%neutral( js )%multiple_states_flag = 1
                transport_ggd(1)%neutral( js )%multiple_states_flag = 1
                do i = 1, nsources
                   allocate( sources_ggd(i)%neutral( js )%state( iss )%name(1) )
                   sources_ggd(i)%neutral( js )%state( iss )%name = spclabel
                   call fill_atom_neutral_type( iatm, &
                  &  sources_ggd(i)%neutral( js )%state( iss )%neutral_type )
                   sources_ggd(i)%neutral( js )%multiple_states_flag = 1
                end do
#endif
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
#if ( IMAS_MAJOR_VERSION < 4 || ( IMAS_MAJOR_VERSION == 4 && IMAS_MINOR_VERSION < 1 ) )
                     allocate( edge_profiles%ggd( time_sind )%neutral( js )%state( isstat(natmi+j-1) ) )
                     allocate( edge_transport%model(1)%ggd( time_sind )%neutral( js )%state( isstat(natmi+j-1) ) )
                     do i = 1, nsources
                        allocate( edge_sources%source(i)%ggd( time_sind )%neutral( js )%state( isstat(natmi+j-1) ) )
                     end do
#else
                     allocate( profiles_ggd%neutral( js )%state( isstat(natmi+j-1) ) )
                     allocate( transport_ggd(1)%neutral( js )%state( isstat(natmi+j-1) ) )
                     do i = 1, nsources
                        allocate( sources_ggd(i)%neutral( js )%state( isstat(natmi+j-1) ) )
                     end do
#endif
                  end if
                  js = js + 1
               end if
               imneut(j) = js
            end do
#if ( IMAS_MAJOR_VERSION < 4 || ( IMAS_MAJOR_VERSION == 4 && IMAS_MINOR_VERSION < 1 ) )
            allocate( edge_profiles%ggd( time_sind )%neutral( js )%state( ks ) )
            allocate( edge_transport%model(1)%ggd( time_sind )%neutral( js )%state( ks ) )
            do i = 1, nsources
               allocate( edge_sources%source(i)%ggd( time_sind )%neutral( js )%state( ks ) )
            end do
#else
            allocate( profiles_ggd%neutral( js )%state( ks ) )
            allocate( transport_ggd(1)%neutral( js )%state( ks ) )
            do i = 1, nsources
               allocate( sources_ggd(i)%neutral( js )%state( ks ) )
            end do
#endif

            js = nspecies
            do j = 1, nmoli
               ks = isstat(natmi+j)
               if (ks.eq.1) js = js + 1
               spclabel = trim(textmn(j-1))
               call shrink_label(spclabel)
               nelems = count ( mlcmp( 1:natmi, j ) > 0 )
#if ( IMAS_MAJOR_VERSION < 4 || ( IMAS_MAJOR_VERSION == 4 && IMAS_MINOR_VERSION < 1 ) )
               allocate( edge_profiles%ggd( time_sind )%neutral( js )%element( nelems ) )
               allocate( edge_transport%model(1)%ggd( time_sind )%neutral( js )%element( nelems ) )
               call fill_molecule_elements( nelems, j, &
                 &  edge_profiles%ggd( time_sind )%neutral( js )%element(:) )
               call fill_molecule_elements( nelems, j, &
                 &  edge_transport%model(1)%ggd( time_sind )%neutral( js )%element(:) )
#if ( IMAS_MAJOR_VERSION < 4 && IMAS_MINOR_VERSION < 42 )
               allocate( edge_profiles%ggd( time_sind )%neutral( js )%label(1) )
               allocate( edge_transport%model(1)%ggd( time_sind )%neutral( js )%label(1) )
               edge_profiles%ggd( time_sind )%neutral( js )%label = spclabel
               edge_transport%model(1)%ggd( time_sind )%neutral( js )%label = spclabel
               allocate( edge_profiles%ggd( time_sind )%neutral( js )%state( ks )%label(1) )
               allocate( edge_transport%model(1)%ggd( time_sind )%neutral( js )%state( ks )%label(1) )
               edge_profiles%ggd( time_sind )%neutral( js )%state( ks )%label = spclabel
               edge_transport%model(1)%ggd( time_sind )%neutral( js )%state( ks )%label = spclabel
#else
               allocate( edge_profiles%ggd( time_sind )%neutral( js )%name(1) )
               allocate( edge_transport%model(1)%ggd( time_sind )%neutral( js )%name(1) )
               edge_profiles%ggd( time_sind )%neutral( js )%name = spclabel
               edge_transport%model(1)%ggd( time_sind )%neutral( js )%name = spclabel
               allocate( edge_profiles%ggd( time_sind )%neutral( js )%state( ks )%name(1) )
               allocate( edge_transport%model(1)%ggd( time_sind )%neutral( js )%state( ks )%name(1) )
               edge_profiles%ggd( time_sind )%neutral( js )%state( ks )%name = spclabel
               edge_transport%model(1)%ggd( time_sind )%neutral( js )%state( ks )%name = spclabel
#endif
               call fill_molecule_type( &
                 &  edge_profiles%ggd( time_sind )%neutral( js )%state( ks )%neutral_type )
               call fill_molecule_type( &
                 &  edge_transport%model(1)%ggd( time_sind )%neutral( js )%state( ks )%neutral_type )
               do i = 1, nsources
                  allocate( edge_sources%source(i)%ggd( time_sind )%neutral( js )%element( nelems ) )
                  call fill_molecule_elements( nelems, j, &
                    &  edge_sources%source(i)%ggd( time_sind )%neutral( js )%element(:) )
#if ( IMAS_MAJOR_VERSION < 4 && IMAS_MINOR_VERSION < 42 )
                  allocate( edge_sources%source(i)%ggd( time_sind )%neutral( js )%label(1) )
                  allocate( edge_sources%source(i)%ggd( time_sind )%neutral( js )%state( ks )%label(1) )
                  edge_sources%source(i)%ggd( time_sind )%neutral( js )%label = spclabel
                  edge_sources%source(i)%ggd( time_sind )%neutral( js )%state( ks )%label = spclabel
#else
                  allocate( edge_sources%source(i)%ggd( time_sind )%neutral( js )%name(1) )
                  allocate( edge_sources%source(i)%ggd( time_sind )%neutral( js )%state( ks )%name(1) )
                  edge_sources%source(i)%ggd( time_sind )%neutral( js )%name = spclabel
                  edge_sources%source(i)%ggd( time_sind )%neutral( js )%state( ks )%name = spclabel
#endif
                  call fill_molecule_type( &
                    &  edge_sources%source(i)%ggd( time_sind )%neutral( js )%state( ks )%neutral_type )
               end do
#else
               allocate( profiles_ggd%neutral( js )%element( nelems ) )
               allocate( profiles_ggd%neutral( js )%name(1) )
               allocate( profiles_ggd%neutral( js )%state( ks )%name(1) )
               allocate( transport_ggd(1)%neutral( js )%element( nelems ) )
               allocate( transport_ggd(1)%neutral( js )%name(1) )
               allocate( transport_ggd(1)%neutral( js )%state( ks )%name(1) )
               call fill_molecule_elements( nelems, j, &
                 &  profiles_ggd%neutral( js )%element(:) )
               call fill_molecule_elements( nelems, j, &
                 &  transport_ggd(1)%neutral( js )%element(:) )
               profiles_ggd%neutral( js )%name = spclabel
               profiles_ggd%neutral( js )%state( ks )%name = spclabel
               call fill_molecule_type( &
                 &  profiles_ggd%neutral( js )%state( ks )%neutral_type )
               call fill_molecule_type( &
                 &  transport_ggd(1)%neutral( js )%state( ks )%neutral_type )
               do i = 1, nsources
                  allocate( sources_ggd(i)%neutral( js )%element( nelems ) )
                  allocate( sources_ggd(i)%neutral( js )%name(1) )
                  allocate( sources_ggd(i)%neutral( js )%state( ks )%name(1) )
                  call fill_molecule_elements( nelems, j, &
                    &  sources_ggd(i)%neutral( js )%element(:) )
                  sources_ggd(i)%neutral( js )%name = spclabel
                  sources_ggd(i)%neutral( js )%state( ks )%name = spclabel
                  call fill_molecule_type( &
                    &  sources_ggd(i)%neutral( js )%state( ks )%neutral_type )
               end do
#endif
               ion_label = trim(spclabel)//'+'
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
#if ( IMAS_MAJOR_VERSION < 4 || ( IMAS_MAJOR_VERSION == 4 && IMAS_MINOR_VERSION < 1 ) )
               edge_profiles%ggd( time_sind )%neutral( js )%ion_index = k
               edge_profiles%ggd( time_sind )%neutral( js )%multiple_states_flag = 1
               edge_transport%model(1)%ggd( time_sind )%neutral( js )%ion_index = k
               edge_transport%model(1)%ggd( time_sind )%neutral( js )%multiple_states_flag = 1
               do i = 1, nsources
                  edge_sources%source(i)%ggd( time_sind )%neutral( js )%ion_index = k
                  edge_sources%source(i)%ggd( time_sind )%neutral( js )%multiple_states_flag = 1
               end do
#else
               profiles_ggd%neutral( js )%ion_index = k
               profiles_ggd%neutral( js )%multiple_states_flag = 1
               transport_ggd(1)%neutral( js )%ion_index = k
               transport_ggd(1)%neutral( js )%multiple_states_flag = 1
               do i = 1, nsources
                  sources_ggd(i)%neutral( js )%ion_index = k
                  sources_ggd(i)%neutral( js )%multiple_states_flag = 1
               end do
#endif
            end do
#endif
        else
            nneut = 0
            do is = 1, nspecies
              if (is_neutral(eb2spcr(is))) nneut=nneut+1
            end do
#if ( IMAS_MAJOR_VERSION < 4 || ( IMAS_MAJOR_VERSION == 4 && IMAS_MINOR_VERSION < 1 ) )
            allocate( edge_profiles%ggd( time_sind )%neutral( nneut ) )
            allocate( edge_transport%model(1)%ggd( time_sind )%neutral( nneut ) )
            do i = 1, nsources
               allocate( edge_sources%source(i)%ggd( time_sind )%neutral( nneut ) )
            end do
#else
            allocate( profiles_ggd%neutral( nneut ) )
            allocate( transport_ggd(1)%neutral( nneut ) )
            do i = 1, nsources
               allocate( sources_ggd(i)%neutral( nneut ) )
            end do
#endif
#if ( IMAS_MINOR_VERSION > 21 || IMAS_MAJOR_VERSION > 3 )
            allocate( radiation%process(1)%ggd( time_sind )%neutral( nneut ) )
#endif
            j = 0
            do js = 1, nspecies
               is = eb2spcr(js)
               if (.not.is_neutral(is)) cycle
               j = j + 1
               call species( is, spclabel, .false. )
               call shrink_label(spclabel)
#if ( IMAS_MAJOR_VERSION < 4 || ( IMAS_MAJOR_VERSION == 4 && IMAS_MINOR_VERSION < 1 ) )
               allocate( edge_profiles%ggd( time_sind )%neutral( j )%element(1) )
               allocate( edge_transport%model(1)%ggd( time_sind )%neutral( j )%element(1) )
               call fill_neutral_element( is, js, &
                 &  edge_profiles%ggd( time_sind )%neutral( j )%element(1) )
               call fill_neutral_element( is, js, &
                 &  edge_transport%model(1)%ggd( time_sind )%neutral( j )%element(1) )
#if ( IMAS_MAJOR_VERSION < 4 && IMAS_MINOR_VERSION < 42 )
               allocate( edge_profiles%ggd( time_sind )%neutral( j )%label(1) )
               allocate( edge_transport%model(1)%ggd( time_sind )%neutral( j )%label(1) )
               edge_profiles%ggd( time_sind )%neutral( j )%label = species_list( js )
               edge_transport%model(1)%ggd( time_sind )%neutral( j )%label = species_list( js )
#else
               allocate( edge_profiles%ggd( time_sind )%neutral( j )%name(1) )
               allocate( edge_transport%model(1)%ggd( time_sind )%neutral( j )%name(1) )
               edge_profiles%ggd( time_sind )%neutral( j )%name = species_list( js )
               edge_transport%model(1)%ggd( time_sind )%neutral( j )%name = species_list( js )
#endif
               allocate( edge_profiles%ggd( time_sind )%neutral( j )%state(1) )
               allocate( edge_transport%model(1)%ggd( time_sind )%neutral( j )%state(1) )
#if ( IMAS_MAJOR_VERSION < 4 && IMAS_MINOR_VERSION < 42 )
               allocate( edge_profiles%ggd( time_sind )%neutral( j )%state(1)%label(1) )
               allocate( edge_transport%model(1)%ggd( time_sind )%neutral( j )%state(1)%label(1) )
               edge_profiles%ggd( time_sind )%neutral( j )%state(1)%label = spclabel
               edge_transport%model(1)%ggd( time_sind )%neutral( j )%state(1)%label = spclabel
#else
               allocate( edge_profiles%ggd( time_sind )%neutral( j )%state(1)%name(1) )
               allocate( edge_transport%model(1)%ggd( time_sind )%neutral( j )%state(1)%name(1) )
               edge_profiles%ggd( time_sind )%neutral( j )%state(1)%name = spclabel
               edge_transport%model(1)%ggd( time_sind )%neutral( j )%state(1)%name = spclabel
#endif
               call fill_atom_neutral_type( js, &
                 &  edge_profiles%ggd( time_sind )%neutral( j )%state( 1 )%neutral_type )
               call fill_atom_neutral_type( js, &
                 &  edge_transport%model(1)%ggd( time_sind )%neutral( j )%state(1)%neutral_type )
               edge_profiles%ggd( time_sind )%neutral( j )%ion_index = js
               edge_profiles%ggd( time_sind )%neutral( j )%multiple_states_flag = 1
               edge_transport%model(1)%ggd( time_sind )%neutral( j )%ion_index = js
               edge_transport%model(1)%ggd( time_sind )%neutral( j )%multiple_states_flag = 1
               do i = 1, nsources
                  allocate( edge_sources%source(i)%ggd( time_sind )%neutral( j )%element(1) )
                  call fill_neutral_element( is, js, &
                    &  edge_sources%source(i)%ggd( time_sind )%neutral( j )%element(1) )
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
                  call fill_atom_neutral_type( js, &
                    &  edge_sources%source(i)%ggd( time_sind )%neutral( j )%state(1)%neutral_type )
                  edge_sources%source(i)%ggd( time_sind )%neutral( j )%ion_index = js
                  edge_sources%source(i)%ggd( time_sind )%neutral( j )%multiple_states_flag = 1
               end do
#else
               allocate( profiles_ggd%neutral( j )%element(1) )
               allocate( profiles_ggd%neutral( j )%name(1) )
               allocate( transport_ggd(1)%neutral( j )%element(1) )
               allocate( transport_ggd(1)%neutral( j )%name(1) )
               call fill_neutral_element( is, js, &
                 &  profiles_ggd%neutral( j )%element(1) )
               call fill_neutral_element( is, js, &
                 &  transport_ggd(1)%neutral( j )%element(1) )
               profiles_ggd%neutral( j )%name = species_list( js )
               transport_ggd(1)%neutral( j )%name = species_list( js )
               allocate( profiles_ggd%neutral( j )%state(1) )
               allocate( profiles_ggd%neutral( j )%state(1)%name(1) )
               allocate( transport_ggd(1)%neutral( j )%state(1) )
               allocate( transport_ggd(1)%neutral( j )%state(1)%name(1) )
               profiles_ggd%neutral( j )%state(1)%name = spclabel
               transport_ggd(1)%neutral( j )%state(1)%name = spclabel
               call fill_atom_neutral_type( js, &
                 &  profiles_ggd%neutral( j )%state( 1 )%neutral_type )
               call fill_atom_neutral_type( js, &
                 &  transport_ggd(1)%neutral( j )%state(1)%neutral_type )
               profiles_ggd%neutral( j )%ion_index = js
               profiles_ggd%neutral( j )%multiple_states_flag = 1
               transport_ggd(1)%neutral( j )%ion_index = js
               transport_ggd(1)%neutral( j )%multiple_states_flag = 1
               do i = 1, nsources
                  allocate( sources_ggd(i)%neutral( j )%element(1) )
                  allocate( sources_ggd(i)%neutral( j )%name(1) )
                  call fill_neutral_element( is, js, &
                 &  sources_ggd(i)%neutral( j )%element(1) )
                  sources_ggd(i)%neutral( j )%name = species_list( js )
                  allocate( sources_ggd(i)%neutral( j )%state(1) )
                  allocate( sources_ggd(i)%neutral( j )%state(1) )
                  sources_ggd(i)%neutral( j )%state(1)%name = spclabel
                  call fill_atom_neutral_type( js, &
                    &  sources_ggd(i)%neutral( j )%state(1)%neutral_type )
                  sources_ggd(i)%neutral( j )%ion_index = js
                  sources_ggd(i)%neutral( j )%multiple_states_flag = 1
               end do
#endif
#if ( IMAS_MINOR_VERSION > 21 || IMAS_MAJOR_VERSION > 3 )
               allocate( radiation%process(1)%ggd( time_sind )%neutral( j )%element(1) )
               call fill_neutral_element( is, js, &
                 &  radiation%process(1)%ggd( time_sind )%neutral( j )%element(1) )
               allocate( radiation%process(1)%ggd( time_sind )%neutral( j )%state(1) )
#if ( IMAS_MAJOR_VERSION < 4 && IMAS_MINOR_VERSION < 42 )
               allocate( radiation%process(1)%ggd( time_sind )%neutral( j )%label(1) )
               allocate( radiation%process(1)%ggd( time_sind )%neutral( j )%state(1)%label(1) )
               radiation%process(1)%ggd( time_sind )%neutral( j )%label = species_list( js )
               radiation%process(1)%ggd( time_sind )%neutral( j )%state(1)%label = spclabel
#else
               allocate( radiation%process(1)%ggd( time_sind )%neutral( j )%name(1) )
               allocate( radiation%process(1)%ggd( time_sind )%neutral( j )%state(1)%name(1) )
               radiation%process(1)%ggd( time_sind )%neutral( j )%name = species_list( js )
               radiation%process(1)%ggd( time_sind )%neutral( j )%state(1)%name = spclabel
#endif
               radiation%process(1)%ggd( time_sind )%neutral( j )%ion_index = js
               radiation%process(1)%ggd( time_sind )%neutral( j )%multiple_states_flag = 1
               call fill_atom_neutral_type( js, &
                 &  radiation%process(1)%ggd( time_sind )%neutral( j )%state(1)%neutral_type )
#endif
            end do
        end if

#if IMAS_MAJOR_VERSION > 3
        allocate( plasma_profiles%ggd( time_sind )%ion( nsion ) )
        do i = 1, nsources
          allocate( plasma_sources%source(i)%ggd( time_sind )%ion( nsion ) )
          allocate( plasma_sources%source(i)%ggd( time_sind )%neutral( nneut ) )
        end do
        allocate( plasma_transport%model(1)%ggd( time_sind )%ion( nsion ) )
        allocate( plasma_transport%model(1)%ggd( time_sind )%neutral( nneut ) )
#endif

#if ( defined(B25_EIRENE) && ( IMAS_MINOR_VERSION > 21 || IMAS_MAJOR_VERSION > 3 ) )
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
            spclabel = trim(textan(is-1))
            call shrink_label(spclabel)
            allocate( radiation%process(3)%ggd( time_sind )%neutral( js )%element(1) )
            call fill_neutral_element( eb2atcr(is), js, &
              &  radiation%process(3)%ggd( time_sind )%neutral( js )%element(1) )
#if ( IMAS_MAJOR_VERSION < 4 && IMAS_MINOR_VERSION < 42 )
            allocate( radiation%process(3)%ggd( time_sind )%neutral( js )%label(1) )
            allocate( radiation%process(3)%ggd( time_sind )%neutral( js )%state( ks )%label(1) )
            radiation%process(3)%ggd( time_sind )%neutral( js )%label = species_list( js )
            radiation%process(3)%ggd( time_sind )%neutral( js )%state( ks )%label = spclabel
#else
            allocate( radiation%process(3)%ggd( time_sind )%neutral( js )%name(1) )
            allocate( radiation%process(3)%ggd( time_sind )%neutral( js )%state( ks )%name(1) )
            radiation%process(3)%ggd( time_sind )%neutral( js )%name = species_list( js )
            radiation%process(3)%ggd( time_sind )%neutral( js )%state( ks )%name = spclabel
#endif
            call fill_atom_neutral_type( is, &
              &  radiation%process(3)%ggd( time_sind )%neutral( js )%state( ks )%neutral_type )
            radiation%process(3)%ggd( time_sind )%neutral( js )%ion_index = js
            radiation%process(3)%ggd( time_sind )%neutral( js )%multiple_states_flag = 1
          end do

          !! List of molecules
          js = nspecies
          do j = 1, nmoli
            ks = isstat(natmi+j)
            if (ks.eq.1) js = js + 1
            spclabel = trim(textmn(j-1))
            call shrink_label(spclabel)
            nelems = count ( mlcmp( 1:natmi, j ) > 0 )
            allocate( radiation%process(3)%ggd( time_sind )%neutral( js )%element( nelems ) )
            call fill_molecule_elements( nelems, j, &
              &  radiation%process(3)%ggd( time_sind )%neutral( js )%element(:) )
#if ( IMAS_MAJOR_VERSION < 4 && IMAS_MINOR_VERSION < 42 )
            allocate( radiation%process(3)%ggd( time_sind )%neutral( js )%label(1) )
            allocate( radiation%process(3)%ggd( time_sind )%neutral( js )%state( ks )%label(1) )
            radiation%process(3)%ggd( time_sind )%neutral( js )%label = spclabel
            radiation%process(3)%ggd( time_sind )%neutral( js )%state( ks )%label = spclabel
#else
            allocate( radiation%process(3)%ggd( time_sind )%neutral( js )%name(1) )
            allocate( radiation%process(3)%ggd( time_sind )%neutral( js )%state( ks )%name(1) )
            radiation%process(3)%ggd( time_sind )%neutral( js )%name = spclabel
            radiation%process(3)%ggd( time_sind )%neutral( js )%state( ks )%name = spclabel
#endif
            call fill_molecule_type( &
              &  radiation%process(3)%ggd( time_sind )%neutral( js )%state( ks )%neutral_type )
            ion_label = trim(spclabel)//'+'
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
            radiation%process(3)%ggd( time_sind )%neutral( js )%ion_index = k
            radiation%process(3)%ggd( time_sind )%neutral( js )%multiple_states_flag = 0
          end do

          !! List of molecular ions
          allocate( radiation%process(4)%ggd( time_sind )%ion( nsion ) )
          do js = nspecies+1, nsion
            allocate( radiation%process(4)%ggd( time_sind )%ion( js )%state( istion(js) ) )
#if ( IMAS_MAJOR_VERSION < 4 && IMAS_MINOR_VERSION < 42 )
            allocate( radiation%process(4)%ggd( time_sind )%ion( js )%label(1) )
            do ks = 1, istion(js)
              is = ispion(js,ks)
              spclabel = trim(textin(is-1))
              call shrink_label(spclabel)
              allocate( radiation%process(4)%ggd( time_sind )%ion( js )%state( ks )%label(1) )
              radiation%process(4)%ggd( time_sind )%ion( js )%state( ks )%label = spclabel
            end do
#else
            allocate( radiation%process(4)%ggd( time_sind )%ion( js )%name(1) )
            do ks = 1, istion(js)
              is = ispion(js,ks)
              spclabel = trim(textin(is-1))
              call shrink_label(spclabel)
              allocate( radiation%process(4)%ggd( time_sind )%ion( js )%state( ks )%name(1) )
              radiation%process(4)%ggd( time_sind )%ion( js )%state( ks )%name = spclabel
            end do
#endif
            is = ispion(js,1)
            spclabel = trim(textin(is-1))
            call shrink_label(spclabel)
            nelems = count ( micmp( 1:natmi, is ) > 0 )
            allocate( radiation%process(4)%ggd( time_sind )%ion( js )%element( nelems ) )
            call fill_mol_ion_elements( nelems, is, &
              &  radiation%process(4)%ggd( time_sind )%ion( js )%element(:) )
            if (istion(js).eq.1) then
              radiation%process(4)%ggd( time_sind )%ion( js )%z_ion = nchrgi( is )
#if ( IMAS_MAJOR_VERSION < 4 && IMAS_MINOR_VERSION < 42 )
              radiation%process(4)%ggd( time_sind )%ion( js )%label = spclabel
#else
              radiation%process(4)%ggd( time_sind )%ion( js )%name = spclabel
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
              call shrink_label(ion_label)
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
!! Obtain the neutral velocities
!! Recall that P[XYZ]DEN[AM] are momentum densities in CGS units!
!! P[UV][XY] will only exist if fort.46 file was written with format 20170930 or later
        if (use_eirene.ne.0 .and. allocated(pux)) then
          allocate(un0(-1:nx,-1:ny,0:2,natmi))
          allocate(um0(-1:nx,-1:ny,0:2,nmoli))
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
            if (LSN) then
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
            if (LSN) then
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
#if ( IMAS_MAJOR_VERSION < 4 || ( IMAS_MAJOR_VERSION == 4 && IMAS_MINOR_VERSION < 1 ) )
            !! ne: Electron density
            call write_quantity( edge_grid,                                 &
                &   val = edge_profiles%ggd( time_sind )%electrons%density, &
                &   value = ne )
            !! fne: Electron particle flux
#if ( IMAS_MINOR_VERSION > 38 || IMAS_MAJOR_VERSION > 3 )
            call divide_by_contact_areas(nx,ny,fne(:,:,:,0),flxFace(:,:,:,0))
            call divide_by_contact_areas(nx,ny,fne(:,:,:,1),flxFace(:,:,:,1))
            call write_face_scalar( transport_grid,                         &
                &   val = edge_transport%model(1)%ggd( time_sind )%         &
                &         electrons%particles%flux_pol,                     &
                &   value = flxFace(:,:,:,0) )
            call write_face_scalar( transport_grid,                         &
                &   val = edge_transport%model(1)%ggd( time_sind )%         &
                &         electrons%particles%flux_radial,                  &
                &   value = flxFace(:,:,:,1) )
#else
            hlpFace(:,:,0) = fne(:,:,0,0)
            hlpFace(:,:,1) = fne(:,:,1,1)
            call divide_by_contact_areas(nx,ny,hlpFace,totFace)
            call write_face_scalar( transport_grid,                         &
                &   val = edge_transport%model(1)%ggd( time_sind )%         &
                &         electrons%particles%flux,                         &
                &   value = totFace )
#endif
            !! sne: Electron particle sources
            tmpCv(:,:) = ( sne(:,:,0) + sne(:,:,1) * ne(:,:) ) / vol(:,:)
            call write_cell_scalar( sources_grid,                           &
                &   scalar = edge_sources%source(1)%ggd( time_sind )%       &
                &            electrons%particles,                           &
                &   b2CellData = tmpCv )
            tmpCv(:,:) = ext_sne(:,:) / vol(:,:)
            call write_cell_scalar( sources_grid,                           &
                &   scalar = edge_sources%source(2)%ggd( time_sind )%       &
                &            electrons%particles,                           &
                &   b2CellData = tmpCv )
            tmpCv(:,:) = ( b2stbc_sne(:,:) + b2stbm_sne(:,:) ) / vol(:,:)
            call write_cell_scalar( sources_grid,                           &
                &   scalar = edge_sources%source(3)%ggd( time_sind )%       &
                &            electrons%particles,                           &
                &   b2CellData = tmpCv )
            tmpCv(:,:) = ( snedt(:,:,0) + snedt(:,:,1) * ne(:,:) ) / vol(:,:)
            call write_cell_scalar( sources_grid,                           &
                &   scalar = edge_sources%source(4)%ggd( time_sind )%       &
                &            electrons%particles,                           &
                &   b2CellData = tmpCv )
            if (use_eirene.ne.0 .and. balance_netcdf.ne.0) then
                tmpCv = 0.0_IDS_real
                do istrai = 1, size( eirene_mc_pael_sne_bal, 3)
                    tmpCv(:,:) = tmpCv(:,:) + eirene_mc_pael_sne_bal(:,:,istrai)
                end do
                tmpCv(:,:) = tmpCv(:,:) / vol(:,:)
                call write_cell_scalar( sources_grid,                       &
                    &   scalar = edge_sources%source(5)%ggd( time_sind )%   &
                    &            electrons%particles,                       &
                    &   b2CellData = tmpCv )
                tmpCv = 0.0_IDS_real
                do istrai = 1, size( eirene_mc_pmel_sne_bal, 3)
                    tmpCv(:,:) = tmpCv(:,:) + eirene_mc_pmel_sne_bal(:,:,istrai)
                end do
                tmpCv(:,:) = tmpCv(:,:) / vol(:,:)
                call write_cell_scalar( sources_grid,                       &
                    &   scalar = edge_sources%source(6)%ggd( time_sind )%   &
                    &            electrons%particles,                       &
                    &   b2CellData = tmpCv )
            else
                tmpCv(:,:) = b2stbr_sne(:,:) / vol(:,:)
                call write_cell_scalar( sources_grid,                       &
                    &   scalar = edge_sources%source(5)%ggd( time_sind )%   &
                    &            electrons%particles,                       &
                    &   b2CellData = tmpCv )
            end if
            tmpCv = 0.0_IDS_real
            do is = 0, ns-1
              tmpCv(:,:) = tmpCv(:,:) + rsana(:,:,is)
            end do
            tmpCv(:,:) = tmpCv(:,:) / vol(:,:)
            call write_cell_scalar( sources_grid,                           &
                &   scalar = edge_sources%source(7)%ggd( time_sind )%       &
                &            electrons%particles,                           &
                &   b2CellData = tmpCv )
            tmpCv = 0.0_IDS_real
            do is = 0, ns-1
              tmpCv(:,:) = tmpCv(:,:) + rrana(:,:,is)
            end do
            tmpCv(:,:) = tmpCv(:,:) / vol(:,:)
            call write_cell_scalar( sources_grid,                           &
                &   scalar = edge_sources%source(8)%ggd( time_sind )%       &
                &            electrons%particles,                           &
                &   b2CellData = tmpCv )
#if ( ( IMAS_MINOR_VERSION > 38 || IMAS_MAJOR_VERSION > 3 ) && defined(B25_EIRENE) )
            if (use_eirene.ne.0) then
              tmpCv(:,:) = sne0_eir_tot(:,:) / vol(:,:)
              call write_cell_scalar( sources_grid,                         &
                &   scalar = edge_sources%source(13)%ggd( time_sind )%      &
                &            electrons%particles,                           &
                &   b2CellData = tmpCv )
            end if
#endif
#else
            !! ne: Electron density
            call write_quantity( edge_grid,                                 &
                &   val = profiles_ggd%electrons%density,                   &
                &   value = ne )
            !! fne: Electron particle flux
            call divide_by_contact_areas(nx,ny,fne(:,:,:,0),flxFace(:,:,:,0))
            call divide_by_contact_areas(nx,ny,fne(:,:,:,1),flxFace(:,:,:,1))
            call write_face_scalar( transport_grid,                         &
                &   val = transport_ggd(1)%electrons%particles%flux_pol,    &
                &   value = flxFace(:,:,:,0) )
            call write_face_scalar( transport_grid,                         &
                &   val = transport_ggd(1)%electrons%particles%flux_radial, &
                &   value = flxFace(:,:,:,1) )
            !! sne: Electron particle sources
            tmpCv(:,:) = ( sne(:,:,0) + sne(:,:,1) * ne(:,:) ) / vol(:,:)
            call write_cell_scalar( sources_grid,                           &
                &   scalar = sources_ggd(1)%electrons%particles,            &
                &   b2CellData = tmpCv )
            tmpCv(:,:) = ext_sne(:,:) / vol(:,:)
            call write_cell_scalar( sources_grid,                           &
                &   scalar = sources_ggd(2)%electrons%particles,            &
                &   b2CellData = tmpCv )
            tmpCv(:,:) = ( b2stbc_sne(:,:) + b2stbm_sne(:,:) ) / vol(:,:)
            call write_cell_scalar( sources_grid,                           &
                &   scalar = sources_ggd(3)%electrons%particles,            &
                &   b2CellData = tmpCv )
            tmpCv(:,:) = ( snedt(:,:,0) + snedt(:,:,1) * ne(:,:) ) / vol(:,:)
            call write_cell_scalar( sources_grid,                           &
                &   scalar = sources_ggd(4)%electrons%particles,            &
                &   b2CellData = tmpCv )
            if (use_eirene.ne.0 .and. balance_netcdf.ne.0) then
                tmpCv = 0.0_IDS_real
                do istrai = 1, size( eirene_mc_pael_sne_bal, 3)
                    tmpCv(:,:) = tmpCv(:,:) + eirene_mc_pael_sne_bal(:,:,istrai)
                end do
                tmpCv(:,:) = tmpCv(:,:) / vol(:,:)
                call write_cell_scalar( sources_grid,                       &
                    &   scalar = sources_ggd(5)%electrons%particles,        &
                    &   b2CellData = tmpCv )
                tmpCv = 0.0_IDS_real
                do istrai = 1, size( eirene_mc_pmel_sne_bal, 3)
                    tmpCv(:,:) = tmpCv(:,:) + eirene_mc_pmel_sne_bal(:,:,istrai)
                end do
                tmpCv(:,:) = tmpCv(:,:) / vol(:,:)
                call write_cell_scalar( sources_grid,                       &
                    &   scalar = sources_ggd(6)%electrons%particles,        &
                    &   b2CellData = tmpCv )
            else
                tmpCv(:,:) = b2stbr_sne(:,:) / vol(:,:)
                call write_cell_scalar( sources_grid,                       &
                    &   scalar = sources_ggd(5)%electrons%particles,        &
                    &   b2CellData = tmpCv )
            end if
            tmpCv = 0.0_IDS_real
            do is = 0, ns-1
              tmpCv(:,:) = tmpCv(:,:) + rsana(:,:,is)
            end do
            tmpCv(:,:) = tmpCv(:,:) / vol(:,:)
            call write_cell_scalar( sources_grid,                           &
                &   scalar = sources_ggd(7)%electrons%particles,            &
                &   b2CellData = tmpCv )
            tmpCv = 0.0_IDS_real
            do is = 0, ns-1
              tmpCv(:,:) = tmpCv(:,:) + rrana(:,:,is)
            end do
            tmpCv(:,:) = tmpCv(:,:) / vol(:,:)
            call write_cell_scalar( sources_grid,                           &
                &   scalar = sources_ggd(8)%electrons%particles,            &
                &   b2CellData = tmpCv )
#ifdef B25_EIRENE
            if (use_eirene.ne.0) then
              tmpCv(:,:) = sne0_eir_tot(:,:) / vol(:,:)
              call write_cell_scalar( sources_grid,                         &
                &   scalar = sources_ggd(13)%electrons%particles,           &
                &   b2CellData = tmpCv )
            end if
#endif
#endif

            !! na: Ion density
            do is = 1, nsion
              if (is.le.nspecies) then
#if ( IMAS_MAJOR_VERSION < 4 || ( IMAS_MAJOR_VERSION == 4 && IMAS_MINOR_VERSION < 1 ) )
                tmpCv(:,:) = 0.0_IDS_real
                do js = 1, istion(is)
                  call write_quantity( edge_grid,                           &
                      &   val = edge_profiles%ggd( time_sind )%             &
                      &         ion( is )%state( js )%density,              &
                      &   value = na(:,:,ispion(is,js)) )
                  tmpCv(:,:) = tmpCv(:,:) + na(:,:,ispion(is,js))
                end do
                call write_quantity( edge_grid,                             &
                    &   val = edge_profiles%ggd( time_sind )%               &
                    &         ion(is)%density,                              &
                    &   value = tmpCv )
            !! fna: Ion particle flux
#if ( IMAS_MINOR_VERSION > 38 || IMAS_MAJOR_VERSION > 3 )
                totflux = 0.0_IDS_real
                do js = 1, istion(is)
                  flxFace = fna(:,:,:,:,ispion(is,js))
                  call divide_by_contact_areas(nx,ny,                       &
                      &   fna(:,:,:,0,ispion(is,js)),flxFace(:,:,:,0))
                  call divide_by_contact_areas(nx,ny,                       &
                      &   fna(:,:,:,1,ispion(is,js)),flxFace(:,:,:,1))
                  totflux = totflux + flxFace
                  call write_face_scalar( transport_grid,                   &
                      &   val = edge_transport%model(1)%ggd( time_sind )%   &
                      &         ion( is )%state( js )%particles%flux_pol,   &
                      &   value = flxFace(:,:,:,0) )
                  call write_face_scalar( transport_grid,                   &
                      &   val = edge_transport%model(1)%ggd( time_sind )%   &
                      &         ion( is )%state( js )%particles%flux_radial,&
                      &   value = flxFace(:,:,:,1) )
                end do
                call write_face_scalar( transport_grid,                     &
                    &   val = edge_transport%model(1)%ggd( time_sind )%     &
                    &         ion( is )%particles%flux_pol,                 &
                    &   value = totflux(:,:,:,0) )
                call write_face_scalar( transport_grid,                     &
                    &   val = edge_transport%model(1)%ggd( time_sind )%     &
                    &         ion( is )%particles%flux_radial,              &
                    &   value = totflux(:,:,:,1) )
#else
                totFace = 0.0_IDS_real
                do js = 1, istion(is)
                  hlpFace(-1:nx,-1:ny,0:1) =                                &
                     &   fna(-1:nx,-1:ny,0:1,0,ispion(is,js)) +             &
                     &   fna(-1:nx,-1:ny,0:1,1,ispion(is,js))
                  call divide_by_contact_areas(nx,ny,hlpFace,tmpFace)
                  totFace(:,:,0) = totFace(:,:,0) + tmpFace(:,:,0)
                  totFace(:,:,1) = totFace(:,:,1) + tmpFace(:,:,1)
                  call write_face_scalar( transport_grid,                   &
                      &   val = edge_transport%model(1)%ggd( time_sind )%   &
                      &         ion( is )%state( js )%particles%flux,       &
                      &   value = tmpFace )
                end do
                call write_face_scalar( transport_grid,                     &
                    &   val = edge_transport%model(1)%ggd( time_sind )%     &
                    &         ion( is )%particles%flux,                     &
                    &   value = totFace )
#endif
            !! cdna: Ion diffusivity
                do js = 1, istion(is)
#if ( IMAS_MINOR_VERSION > 38 || IMAS_MAJOR_VERSION > 3 )
                  call write_face_scalar( transport_grid,                   &
                      &   val = edge_transport%model(1)%ggd( time_sind )%   &
                      &         ion( is )%state( js )%particles%d_pol,      &
                      &   value = cdna(:,:,:,0,ispion(is,js)) )
                  call write_face_scalar( transport_grid,                   &
                      &   val = edge_transport%model(1)%ggd( time_sind )%   &
                      &         ion( is )%state( js )%particles%d_radial,   &
                      &   value = cdna(:,:,:,1,ispion(is,js)) )
            !! cvla: Ion diffusivity
                  call write_face_scalar( transport_grid,                   &
                      &   val = edge_transport%model(1)%ggd( time_sind )%   &
                      &         ion( is )%state( js )%particles%v_pol,      &
                      &   value = cvla(:,:,:,0,ispion(is,js)) )
                  call write_face_scalar( transport_grid,                   &
                      &   val = edge_transport%model(1)%ggd( time_sind )%   &
                      &         ion( is )%state( js )%particles%v_radial,   &
                      &   value = cvla(:,:,:,1,ispion(is,js)) )
#else
                  hlpFace(:,:,0) = cdna(:,:,0,0,ispion(is,js))
                  hlpFace(:,:,1) = cdna(:,:,1,1,ispion(is,js))
                  call write_face_scalar( transport_grid,                   &
                      &   val = edge_transport%model(1)%ggd( time_sind )%   &
                      &         ion( is )%state( js )%particles%d,          &
                      &   value = hlpFace )
            !! cvla: Ion diffusivity
                  hlpFace(:,:,0) = cvla(:,:,0,0,ispion(is,js))
                  hlpFace(:,:,1) = cvla(:,:,1,1,ispion(is,js))
                  call write_face_scalar( transport_grid,                   &
                      &   val = edge_transport%model(1)%ggd( time_sind )%   &
                      &         ion( is )%state( js )%particles%v,          &
                      &   value = hlpFace )
#endif
            !! fllim0fna: Ion flux limiter
                  hlpFace(:,:,0) = fllim0fna(:,:,0,0,ispion(is,js))
                  hlpFace(:,:,1) = fllim0fna(:,:,1,1,ispion(is,js))
                  call write_face_scalar( transport_grid,                   &
                      &   val = edge_transport%model(1)%ggd( time_sind )%   &
                      &         ion( is )%state( js )%                      &
                      &         particles%flux_limiter,                     &
                      &   value = hlpFace )
                end do
            !! sna: Ion particle sources
                totCv(:,:) = 0.0_IDS_real
                do js = 1, istion(is)
                  tmpCv(:,:) = ( sna(:,:,0, ispion(is,js) ) +               &
                      &          sna(:,:,1, ispion(is,js) ) * na(:,:, ispion(is,js) ) ) / vol(:,:)
                  totCv(:,:) = totCv(:,:) + tmpCv(:,:)
                  call write_cell_scalar( sources_grid,                     &
                      &   scalar = edge_sources%source(1)%                  &
                      &   ggd( time_sind )%ion( is )%state( js )%particles, &
                      &   b2CellData = tmpCv )
                end do
                call write_cell_scalar( sources_grid,                       &
                      &   scalar = edge_sources%source(1)%                  &
                      &   ggd( time_sind )%ion( is )%particles,             &
                      &   b2CellData = totCv )
                totCv(:,:) = 0.0_IDS_real
                do js = 1, istion(is)
                  tmpCv(:,:) = ext_sna(:,:,ispion(is,js)) / vol(:,:)
                  totCv(:,:) = totCv(:,:) + tmpCv(:,:)
                  call write_cell_scalar( sources_grid,                     &
                      &   scalar = edge_sources%source(2)%ggd( time_sind )% &
                      &            ion( is )%state( js )%particles,         &
                      &   b2CellData = tmpCv )
                end do
                call write_cell_scalar( sources_grid,                       &
                    &   scalar = edge_sources%source(2)%ggd( time_sind )%   &
                    &            ion( is )%particles,                       &
                    &   b2CellData = totCv )
                totCv(:,:) = 0.0_IDS_real
                do js = 1, istion(is)
                  tmpCv(:,:) = ( b2stbc_sna(:,:,ispion(is,js)) +            &
                        &        b2stbm_sna(:,:,ispion(is,js)) ) / vol(:,:)
                  totCv(:,:) = totCv(:,:) + tmpCv(:,:)
                  call write_cell_scalar( sources_grid,                     &
                      &   scalar = edge_sources%source(3)%ggd( time_sind )% &
                      &            ion( is )%state( js )%particles,         &
                      &   b2CellData = tmpCv )
                end do
                call write_cell_scalar( sources_grid,                       &
                    &   scalar = edge_sources%source(3)%ggd( time_sind )%   &
                    &            ion( is )%particles,                       &
                    &   b2CellData = totCv )
                totCv(:,:) = 0.0_IDS_real
                do js = 1, istion(is)
                  tmpCv(:,:) = ( snadt(:,:,0,ispion(is,js)) +               &
                        &        snadt(:,:,1,ispion(is,js)) * na(:,:,ispion(is,js)) ) / vol(:,:)
                  totCv(:,:) = totCv(:,:) + tmpCv(:,:)
                  call write_cell_scalar( sources_grid,                     &
                      &   scalar = edge_sources%source(4)%ggd( time_sind )% &
                      &            ion( is )%state( js )%particles,         &
                      &   b2CellData = tmpCv )
                end do
                call write_cell_scalar( sources_grid,                       &
                    &   scalar = edge_sources%source(4)%ggd( time_sind )%   &
                    &            ion( is )%particles,                       &
                    &   b2CellData = totCv )
                if (use_eirene.ne.0 .and. balance_netcdf.ne.0) then
                  totCv(:,:) = 0.0_IDS_real
                  do js = 1, istion(is)
                    tmpCv = 0.0_IDS_real
                    do istrai = 1, size( eirene_mc_papl_sna_bal, 4)
                      tmpCv(:,:) = tmpCv(:,:)                               &
                         &       + eirene_mc_papl_sna_bal(:,:,ispion(is,js),istrai)
                    end do
                    tmpCv(:,:) = tmpCv(:,:) / vol(:,:)
                    totCv(:,:) = totCv(:,:) + tmpCv(:,:)
                    call write_cell_scalar( sources_grid,                     &
                        &   scalar = edge_sources%source(5)%ggd( time_sind )% &
                        &            ion( is )%state( js )%particles,         &
                        &   b2CellData = tmpCv )
                  end do
                  call write_cell_scalar( sources_grid,                       &
                        &   scalar = edge_sources%source(5)%ggd( time_sind )% &
                        &            ion( is )%particles,                     &
                        &   b2CellData = totCv )
                  totCv = 0.0_IDS_real
                  do js = 1, istion(is)
                    tmpCv = 0.0_IDS_real
                    do istrai = 1, size( eirene_mc_pmpl_sna_bal, 4)
                      tmpCv(:,:) = tmpCv(:,:)                                 &
                         &       + eirene_mc_pmpl_sna_bal(:,:,ispion(is,js),istrai)
                    end do
                    tmpCv(:,:) = tmpCv(:,:) / vol(:,:)
                    totCv(:,:) = totCv(:,:) + tmpCv(:,:)
                    call write_cell_scalar( sources_grid,                     &
                        &   scalar = edge_sources%source(6)%ggd( time_sind )% &
                        &            ion( is )%state( js )%particles,         &
                        &   b2CellData = tmpCv )
                  end do
                  call write_cell_scalar( sources_grid,                       &
                      &   scalar = edge_sources%source(6)%ggd( time_sind )%   &
                      &            ion( is )%particles,                       &
                      &   b2CellData = totCv )
                else
                  totCv(:,:) = 0.0_IDS_real
                  do js = 1, istion(is)
                    tmpCv(:,:) = b2stbr_sna(:,:,ispion(is,js)) / vol(:,:)
                    totCv(:,:) = totCv(:,:) + tmpCv(:,:)
                    call write_cell_scalar( sources_grid,                     &
                        &   scalar = edge_sources%source(5)%ggd( time_sind )% &
                        &            ion( is )%state( js )%particles,         &
                        &   b2CellData = tmpCv )
                  end do
                  call write_cell_scalar( sources_grid,                       &
                        &   scalar = edge_sources%source(5)%ggd( time_sind )% &
                        &            ion( is )%particles,                     &
                        &   b2CellData = totCv )
                end if
                do js = 1, istion(is)
                  tmpCv(:,:) = rsana(:,:,ispion(is,js)) / vol(:,:)
                  call write_cell_scalar( sources_grid,                       &
                      &   scalar = edge_sources%source(7)%ggd( time_sind )%   &
                      &            ion( is )%state( js )%particles,           &
                      &   b2CellData = tmpCv )
                  tmpCv(:,:) = rrana(:,:,ispion(is,js)) / vol(:,:)
                  call write_cell_scalar( sources_grid,                       &
                      &   scalar = edge_sources%source(8)%ggd( time_sind )%   &
                      &            ion( is )%state( js )%particles,           &
                      &   b2CellData = tmpCv )
                  tmpCv(:,:) = rcxna(:,:,ispion(is,js)) / vol(:,:)
                  call write_cell_scalar( sources_grid,                       &
                      &   scalar = edge_sources%source(9)%ggd( time_sind )%   &
                      &            ion( is )%state( js )%particles,           &
                      &   b2CellData = tmpCv )
                end do
#if ( ( IMAS_MINOR_VERSION > 38 || IMAS_MAJOR_VERSION > 3 ) && defined(B25_EIRENE) )
                if (use_eirene.ne.0) then
                  do js = 1, istion(is)
                    tmpCv(:,:) = sna0_eir_tot(:,:,ispion(is,js)) / vol(:,:)
                    call write_cell_scalar( sources_grid,                     &
                      &   scalar = edge_sources%source(13)%ggd( time_sind )%  &
                      &            ion( is )%state( js )%particles,           &
                      &   b2CellData = tmpCv )
                  end do
                end if
#endif
#else
                tmpCv(:,:) = 0.0_IDS_real
                do js = 1, istion(is)
                  call write_quantity( edge_grid,                           &
                      &   val = profiles_ggd%ion( is )%state( js )%density, &
                      &   value = na(:,:,ispion(is,js)) )
                  tmpCv(:,:) = tmpCv(:,:) + na(:,:,ispion(is,js))
                end do
                call write_quantity( edge_grid,                             &
                    &   val = profiles_ggd%ion(is)%density,                 &
                    &   value = tmpCv )
            !! fna: Ion particle flux
                totflux = 0.0_IDS_real
                do js = 1, istion(is)
                  flxFace = fna(:,:,:,:,ispion(is,js))
                  call divide_by_contact_areas(nx,ny,                       &
                      &   fna(:,:,:,0,ispion(is,js)),flxFace(:,:,:,0))
                  call divide_by_contact_areas(nx,ny,                       &
                      &   fna(:,:,:,1,ispion(is,js)),flxFace(:,:,:,1))
                  totflux = totflux + flxFace
                  call write_face_scalar( transport_grid,                   &
                      &   val = transport_ggd(1)%                           &
                      &         ion( is )%state( js )%particles%flux_pol,   &
                      &   value = flxFace(:,:,:,0) )
                  call write_face_scalar( transport_grid,                   &
                      &   val = transport_ggd(1)%                           &
                      &         ion( is )%state( js )%particles%flux_radial,&
                      &   value = flxFace(:,:,:,1) )
                end do
                call write_face_scalar( transport_grid,                     &
                    &   val = transport_ggd(1)%                             &
                    &         ion( is )%particles%flux_pol,                 &
                    &   value = totflux(:,:,:,0) )
                call write_face_scalar( transport_grid,                     &
                    &   val = transport_ggd(1)%                             &
                    &         ion( is )%particles%flux_radial,              &
                    &   value = totflux(:,:,:,1) )
            !! cdna: Ion diffusivity
                do js = 1, istion(is)
                  call write_face_scalar( transport_grid,                   &
                      &   val = transport_ggd(1)%                           &
                      &         ion( is )%state( js )%particles%d_pol,      &
                      &   value = cdna(:,:,:,0,ispion(is,js)) )
                  call write_face_scalar( transport_grid,                   &
                      &   val = transport_ggd(1)%                           &
                      &         ion( is )%state( js )%particles%d_radial,   &
                      &   value = cdna(:,:,:,1,ispion(is,js)) )
            !! cvla: Ion diffusivity
                  call write_face_scalar( transport_grid,                   &
                      &   val = transport_ggd(1)%                           &
                      &         ion( is )%state( js )%particles%v_pol,      &
                      &   value = cvla(:,:,:,0,ispion(is,js)) )
                  call write_face_scalar( transport_grid,                   &
                      &   val = transport_ggd(1)%                           &
                      &         ion( is )%state( js )%particles%v_radial,   &
                      &   value = cvla(:,:,:,1,ispion(is,js)) )
            !! fllim0fna: Ion flux limiter
                  hlpFace(:,:,0) = fllim0fna(:,:,0,0,ispion(is,js))
                  hlpFace(:,:,1) = fllim0fna(:,:,1,1,ispion(is,js))
                  call write_face_scalar( transport_grid,                   &
                      &   val = transport_ggd(1)%ion( is )%state( js )%     &
                      &         particles%flux_limiter,                     &
                      &   value = hlpFace )
                end do
            !! sna: Ion particle sources
                totCv(:,:) = 0.0_IDS_real
                do js = 1, istion(is)
                  tmpCv(:,:) = ( sna(:,:,0, ispion(is,js) ) +               &
                      &          sna(:,:,1, ispion(is,js) ) * na(:,:, ispion(is,js) ) ) / vol(:,:)
                  totCv(:,:) = totCv(:,:) + tmpCv(:,:)
                  call write_cell_scalar( sources_grid,                     &
                      &   scalar = sources_ggd(1)%                          &
                      &            ion( is )%state( js )%particles,         &
                      &   b2CellData = tmpCv )
                end do
                call write_cell_scalar( sources_grid,                       &
                      &   scalar = sources_ggd(1)%ion( is )%particles,      &
                      &   b2CellData = totCv )
                totCv(:,:) = 0.0_IDS_real
                do js = 1, istion(is)
                  tmpCv(:,:) = ext_sna(:,:,ispion(is,js)) / vol(:,:)
                  totCv(:,:) = totCv(:,:) + tmpCv(:,:)
                  call write_cell_scalar( sources_grid,                     &
                      &   scalar = sources_ggd(2)%                          &
                      &            ion( is )%state( js )%particles,         &
                      &   b2CellData = tmpCv )
                end do
                call write_cell_scalar( sources_grid,                       &
                    &   scalar = sources_ggd(2)%ion( is )%particles,        &
                    &   b2CellData = totCv )
                totCv(:,:) = 0.0_IDS_real
                do js = 1, istion(is)
                  tmpCv(:,:) = ( b2stbc_sna(:,:,ispion(is,js)) +            &
                        &        b2stbm_sna(:,:,ispion(is,js)) ) / vol(:,:)
                  totCv(:,:) = totCv(:,:) + tmpCv(:,:)
                  call write_cell_scalar( sources_grid,                     &
                      &   scalar = sources_ggd(3)%                          &
                      &            ion( is )%state( js )%particles,         &
                      &   b2CellData = tmpCv )
                end do
                call write_cell_scalar( sources_grid,                       &
                    &   scalar = sources_ggd(3)%ion( is )%particles,        &
                    &   b2CellData = totCv )
                totCv(:,:) = 0.0_IDS_real
                do js = 1, istion(is)
                  tmpCv(:,:) = ( snadt(:,:,0,ispion(is,js)) +               &
                        &        snadt(:,:,1,ispion(is,js)) * na(:,:,ispion(is,js)) ) / vol(:,:)
                  totCv(:,:) = totCv(:,:) + tmpCv(:,:)
                  call write_cell_scalar( sources_grid,                     &
                      &   scalar = sources_ggd(4)%                          &
                      &            ion( is )%state( js )%particles,         &
                      &   b2CellData = tmpCv )
                end do
                call write_cell_scalar( sources_grid,                       &
                    &   scalar = sources_ggd(4)%ion( is )%particles,        &
                    &   b2CellData = totCv )
                if (use_eirene.ne.0 .and. balance_netcdf.ne.0) then
                  totCv(:,:) = 0.0_IDS_real
                  do js = 1, istion(is)
                    tmpCv = 0.0_IDS_real
                    do istrai = 1, size( eirene_mc_papl_sna_bal, 4)
                      tmpCv(:,:) = tmpCv(:,:)                               &
                         &       + eirene_mc_papl_sna_bal(:,:,ispion(is,js),istrai)
                    end do
                    tmpCv(:,:) = tmpCv(:,:) / vol(:,:)
                    totCv(:,:) = totCv(:,:) + tmpCv(:,:)
                    call write_cell_scalar( sources_grid,                   &
                        &   scalar = sources_ggd(5)%                        &
                        &            ion( is )%state( js )%particles,       &
                        &   b2CellData = tmpCv )
                  end do
                  call write_cell_scalar( sources_grid,                     &
                        &   scalar = sources_ggd(5)%ion( is )%particles,    &
                        &   b2CellData = totCv )
                  totCv = 0.0_IDS_real
                  do js = 1, istion(is)
                    tmpCv = 0.0_IDS_real
                    do istrai = 1, size( eirene_mc_pmpl_sna_bal, 4)
                      tmpCv(:,:) = tmpCv(:,:)                               &
                         &       + eirene_mc_pmpl_sna_bal(:,:,ispion(is,js),istrai)
                    end do
                    tmpCv(:,:) = tmpCv(:,:) / vol(:,:)
                    totCv(:,:) = totCv(:,:) + tmpCv(:,:)
                    call write_cell_scalar( sources_grid,                   &
                        &   scalar = sources_ggd(6)%                        &
                        &            ion( is )%state( js )%particles,       &
                        &   b2CellData = tmpCv )
                  end do
                  call write_cell_scalar( sources_grid,                     &
                      &   scalar = sources_ggd(6)%ion( is )%particles,      &
                      &   b2CellData = totCv )
                else
                  totCv(:,:) = 0.0_IDS_real
                  do js = 1, istion(is)
                    tmpCv(:,:) = b2stbr_sna(:,:,ispion(is,js)) / vol(:,:)
                    totCv(:,:) = totCv(:,:) + tmpCv(:,:)
                    call write_cell_scalar( sources_grid,                   &
                        &   scalar = sources_ggd(5)%                        &
                        &            ion( is )%state( js )%particles,       &
                        &   b2CellData = tmpCv )
                  end do
                  call write_cell_scalar( sources_grid,                     &
                        &   scalar = sources_ggd(5)%ion( is )%particles,    &
                        &   b2CellData = totCv )
                end if
                do js = 1, istion(is)
                  tmpCv(:,:) = rsana(:,:,ispion(is,js)) / vol(:,:)
                  call write_cell_scalar( sources_grid,                     &
                      &   scalar = sources_ggd(7)%                          &
                      &            ion( is )%state( js )%particles,         &
                      &   b2CellData = tmpCv )
                  tmpCv(:,:) = rrana(:,:,ispion(is,js)) / vol(:,:)
                  call write_cell_scalar( sources_grid,                     &
                      &   scalar = sources_ggd(8)%                          &
                      &            ion( is )%state( js )%particles,         &
                      &   b2CellData = tmpCv )
                  tmpCv(:,:) = rcxna(:,:,ispion(is,js)) / vol(:,:)
                  call write_cell_scalar( sources_grid,                     &
                      &   scalar = sources_ggd(9)%                          &
                      &            ion( is )%state( js )%particles,         &
                      &   b2CellData = tmpCv )
                end do
#ifdef B25_EIRENE
                if (use_eirene.ne.0) then
                  do js = 1, istion(is)
                    tmpCv(:,:) = sna0_eir_tot(:,:,ispion(is,js)) / vol(:,:)
                    call write_cell_scalar( sources_grid,                   &
                      &   scalar = sources_ggd(13)%                         &
                      &            ion( is )%state( js )%particles,         &
                      &   b2CellData = tmpCv )
                  end do
                end if
#endif
#endif
              else
                totCv(:,:) = 0.0_IDS_real
                do js = 1, istion(is)
                  tmpCv(-1:nx,-1:ny) = dib2(0:nx+1,0:ny+1,ispion(is,js),1)
                  totCv(-1:nx,-1:ny) = totCv(-1:nx,-1:ny) + tmpCv(-1:nx,-1:ny)
#if ( IMAS_MAJOR_VERSION < 4 || ( IMAS_MAJOR_VERSION == 4 && IMAS_MINOR_VERSION < 1 ) )
                  call write_quantity( edge_grid,                           &
                      &   val = edge_profiles%ggd( time_sind )%             &
                      &         ion(is)%state( js )%density,                &
                      &   value = tmpCv )
#else
                  call write_quantity( edge_grid,                           &
                      &   val = profiles_ggd%ion(is)%state( js )%density,   &
                      &   value = tmpCv )
#endif
                end do
#if ( IMAS_MAJOR_VERSION < 4 || ( IMAS_MAJOR_VERSION == 4 && IMAS_MINOR_VERSION < 1 ) )
                call write_quantity( edge_grid,                             &
                    &   val = edge_profiles%ggd( time_sind )%               &
                    &         ion(is)%density,                              &
                    &   value = totCv )
#else
                call write_quantity( edge_grid,                             &
                    &   val = profiles_ggd%ion(is)%density,                 &
                    &   value = totCv )
#endif
                if (balance_netcdf.ne.0) then
                  totCv = 0.0_IDS_real
                  do js = 1, istion(is)
                    tmpCv = 0.0_IDS_real
                    do istrai = 1, size(eirene_mc_paio_sna_bal,4)
                      tmpCv(:,:) = tmpCv(:,:)                               &
                         &       + eirene_mc_paio_sna_bal(:,:,ispion(is,js),istrai)    &
                         &       + eirene_mc_pmio_sna_bal(:,:,ispion(is,js),istrai)    &
                         &       + eirene_mc_piio_sna_bal(:,:,ispion(is,js),istrai)
                    end do
                    tmpCv(:,:) = tmpCv(:,:) / vol(:,:)
                    totCv(:,:) = totCv(:,:) + tmpCv(:,:)
#if ( IMAS_MAJOR_VERSION < 4 || ( IMAS_MAJOR_VERSION == 4 && IMAS_MINOR_VERSION < 1 ) )
                    call write_cell_scalar( sources_grid,                     &
                        &   scalar = edge_sources%source(1)%ggd( time_sind )% &
                        &            ion( is )%state( js )%particles,         &
                        &   b2CellData = tmpCv )
#else
                    call write_cell_scalar( sources_grid,                     &
                        &   scalar = sources_ggd(1)%                          &
                        &            ion( is )%state( js )%particles,         &
                        &   b2CellData = tmpCv )
#endif
                  end do
#if ( IMAS_MAJOR_VERSION < 4 || ( IMAS_MAJOR_VERSION == 4 && IMAS_MINOR_VERSION < 1 ) )
                  call write_cell_scalar( sources_grid,                       &
                      &   scalar = edge_sources%source(1)%ggd( time_sind )%   &
                      &            ion( is )%particles,                       &
                      &   b2CellData = totCv )
#else
                  call write_cell_scalar( sources_grid,                       &
                      &   scalar = sources_ggd(1)%ion( is )%particles,        &
                      &   b2CellData = totCv )
#endif
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
                    ue(ix,iy)=-(fch_p(ix,iy,0,0)/pbs(ix,iy,0)+             &
                      &         fch_p(rightix(ix,iy),rightiy(ix,iy),0,0)/  &
                      &         pbs(rightix(ix,iy),rightiy(ix,iy),0))/qe/2.0_IDS_real
                  elseif(leftix(ix,iy).eq.-2) then
                    ue(ix,iy)=-fch_p(rightix(ix,iy),rightiy(ix,iy),0,0)/   &
                      &        pbs(rightix(ix,iy),rightiy(ix,iy),0)/qe
                  elseif(rightix(ix,iy).eq.nx+1) then
                    ue(ix,iy)=-fch_p(ix,iy,0,0)/pbs(ix,iy,0)/qe
                  endif
                enddo
              enddo
              do is=0,ns-1
                ue(:,:)=ue(:,:)+rza(:,:,is)*ua(:,:,is)*na(:,:,is)
              enddo
              ue(:,:)=ue(:,:)+ue_ext(:,:)*ne_ext(:,:)
              ue(:,:)=ue(:,:)/ne(:,:)
            endif
#if ( IMAS_MAJOR_VERSION < 4 || ( IMAS_MAJOR_VERSION == 4 && IMAS_MINOR_VERSION < 1 ) )
            call write_cell_vector_component( edge_grid,                 &
                 &   vectorComponent = edge_profiles%ggd( time_sind )%   &
                 &                     electrons%velocity,               &
                 &   b2CellData = ue(:,:),                               &
                 &   vectorID = VEC_ALIGN_PARALLEL_ID )
#else
            call write_cell_vector_component( edge_grid,                 &
                 &   vectorComponent = profiles_ggd%electrons%velocity,  &
                 &   b2CellData = ue(:,:),                               &
                 &   vectorID = VEC_ALIGN_PARALLEL_ID )
#endif

            do is = 1, nsion
              if (is.le.nspecies) then
                do js = 1, istion(is)
#if ( IMAS_MAJOR_VERSION < 4 || ( IMAS_MAJOR_VERSION == 4 && IMAS_MINOR_VERSION < 1 ) )
                !! ua: Parallel ion velocity
                  call write_cell_vector_component( edge_grid,              &
                      &   vectorComponent = edge_profiles%ggd( time_sind )% &
                      &                     ion( is )%state( js )%velocity, &
                      &   b2CellData = ua(:,:,ispion(is,js)),               &
                      &   vectorID = VEC_ALIGN_PARALLEL_ID )
                  if (drift_style.eq.0) then
                !! wadia: Diamagnetic ion velocity
                    call write_cell_vector_component( edge_grid,            &
                      &   vectorComponent = edge_profiles%ggd( time_sind )% &
                      &                     ion( is )%state( js )%          &
                      &                     velocity_diamagnetic,           &
                      &   b2CellData = wadia(:,:,0,0,ispion(is,js)),        &
                      &   vectorID = VEC_ALIGN_POLOIDAL_ID )
                    call write_cell_vector_component( edge_grid,            &
                      &   vectorComponent = edge_profiles%ggd( time_sind )% &
                      &                     ion( is )%state( js )%          &
                      &                     velocity_diamagnetic,           &
                      &   b2CellData = wadia(:,:,1,1,ispion(is,js)),        &
                      &   vectorID = VEC_ALIGN_RADIAL_ID )
                !! vaecrb: ExB ion velocity
                    call write_cell_vector_component( edge_grid,            &
                      &   vectorComponent = edge_profiles%ggd( time_sind )% &
                      &                     ion( is )%state( js )%          &
                      &                     velocity_exb,                   &
                      &   b2CellData = vaecrb(:,:,0,0,ispion(is,js)),       &
                      &   vectorID = VEC_ALIGN_POLOIDAL_ID )
                    call write_cell_vector_component( edge_grid,            &
                      &   vectorComponent = edge_profiles%ggd( time_sind )% &
                      &                     ion( is )%state( js )%          &
                      &                     velocity_exb,                   &
                      &   b2CellData = vaecrb(:,:,1,1,ispion(is,js)),       &
                      &   vectorID = VEC_ALIGN_RADIAL_ID )
                  else
                !! wadia: Diamagnetic ion velocity
                    tmpFace(:,:,0) = wadia(:,:,0,0,ispion(is,js))
                    tmpFace(:,:,1) = wadia(:,:,1,0,ispion(is,js))
                    call write_face_vector_component( edge_grid,              &
                      &   vectorComponent = edge_profiles%ggd( time_sind )%   &
                      &         ion( is )%state( js )%velocity_diamagnetic,   &
                      &   b2FaceData = tmpFace,                               &
                      &   vectorID = VEC_ALIGN_POLOIDAL_ID )
                    tmpFace(:,:,0) = wadia(:,:,0,1,ispion(is,js))
                    tmpFace(:,:,1) = wadia(:,:,1,1,ispion(is,js))
                    call write_face_vector_component( edge_grid,              &
                      &   vectorComponent = edge_profiles%ggd( time_sind )%   &
                      &         ion( is )%state( js )%velocity_diamagnetic,   &
                      &   b2FaceData = tmpFace,                               &
                      &   vectorID = VEC_ALIGN_RADIAL_ID )
                !! vaecrb: ExB ion velocity
                    tmpFace(:,:,0) = vaecrb(:,:,0,0,ispion(is,js))
                    tmpFace(:,:,1) = vaecrb(:,:,1,0,ispion(is,js))
                    call write_face_vector_component( edge_grid,              &
                      &   vectorComponent = edge_profiles%ggd( time_sind )%   &
                      &         ion( is )%state( js )%velocity_exb,           &
                      &   b2FaceData = tmpFace,                               &
                      &   vectorID = VEC_ALIGN_POLOIDAL_ID )
                    tmpFace(:,:,0) = vaecrb(:,:,0,1,ispion(is,js))
                    tmpFace(:,:,1) = vaecrb(:,:,1,1,ispion(is,js))
                    call write_face_vector_component( edge_grid,              &
                      &   vectorComponent = edge_profiles%ggd( time_sind )%   &
                      &         ion( is )%state( js )%velocity_exb,           &
                      &   b2FaceData = tmpFace,                               &
                      &   vectorID = VEC_ALIGN_RADIAL_ID )
                  end if
                !! cvsa: Ion diffusivity
                  tmpFace(:,:,0) = cvsa(:,:,0,0,ispion(is,js))
                  tmpFace(:,:,1) = cvsa(:,:,1,1,ispion(is,js))
                  call write_face_vector_component( transport_grid,           &
                      &   vectorComponent = edge_transport%model(1)%          &
                      &         ggd( time_sind )%ion( is )%state( js )%       &
                      &         momentum%d,                                   &
                      &   b2FaceData = tmpFace,                               &
                      &   vectorID = VEC_ALIGN_PARALLEL_ID )
#else
                !! ua: Parallel ion velocity
                  call write_cell_vector_component( edge_grid,              &
                      &   vectorComponent = profiles_ggd%                   &
                      &                     ion( is )%state( js )%velocity, &
                      &   b2CellData = ua(:,:,ispion(is,js)),               &
                      &   vectorID = VEC_ALIGN_PARALLEL_ID )
                  if (drift_style.eq.0) then
                !! wadia: Diamagnetic ion velocity
                    call write_cell_vector_component( edge_grid,                &
                      &   vectorComponent = profiles_ggd%ion( is )%state( js )% &
                      &                     velocity_diamagnetic,               &
                      &   b2CellData = wadia(:,:,0,0,ispion(is,js)),            &
                      &   vectorID = VEC_ALIGN_POLOIDAL_ID )
                    call write_cell_vector_component( edge_grid,                &
                      &   vectorComponent = profiles_ggd%ion( is )%state( js )% &
                      &                     velocity_diamagnetic,               &
                      &   b2CellData = wadia(:,:,1,1,ispion(is,js)),            &
                      &   vectorID = VEC_ALIGN_RADIAL_ID )
                !! vaecrb: ExB ion velocity
                    call write_cell_vector_component( edge_grid,                &
                      &   vectorComponent = profiles_ggd%ion( is )%state( js )% &
                      &                     velocity_exb,                       &
                      &   b2CellData = vaecrb(:,:,0,0,ispion(is,js)),           &
                      &   vectorID = VEC_ALIGN_POLOIDAL_ID )
                    call write_cell_vector_component( edge_grid,                &
                      &   vectorComponent = profiles_ggd%ion( is )%state( js )% &
                      &                     velocity_exb,                       &
                      &   b2CellData = vaecrb(:,:,1,1,ispion(is,js)),           &
                      &   vectorID = VEC_ALIGN_RADIAL_ID )
                  else
                !! wadia: Diamagnetic ion velocity
                    tmpFace(:,:,0) = wadia(:,:,0,0,ispion(is,js))
                    tmpFace(:,:,1) = wadia(:,:,1,0,ispion(is,js))
                    call write_face_vector_component( edge_grid,              &
                      &   vectorComponent = profiles_ggd%                     &
                      &         ion( is )%state( js )%velocity_diamagnetic,   &
                      &   b2FaceData = tmpFace,                               &
                      &   vectorID = VEC_ALIGN_POLOIDAL_ID )
                    tmpFace(:,:,0) = wadia(:,:,0,1,ispion(is,js))
                    tmpFace(:,:,1) = wadia(:,:,1,1,ispion(is,js))
                    call write_face_vector_component( edge_grid,              &
                      &   vectorComponent = profiles_ggd%                     &
                      &         ion( is )%state( js )%velocity_diamagnetic,   &
                      &   b2FaceData = tmpFace,                               &
                      &   vectorID = VEC_ALIGN_RADIAL_ID )
                !! vaecrb: ExB ion velocity
                    tmpFace(:,:,0) = vaecrb(:,:,0,0,ispion(is,js))
                    tmpFace(:,:,1) = vaecrb(:,:,1,0,ispion(is,js))
                    call write_face_vector_component( edge_grid,              &
                      &   vectorComponent = profiles_ggd%                     &
                      &         ion( is )%state( js )%velocity_exb,           &
                      &   b2FaceData = tmpFace,                               &
                      &   vectorID = VEC_ALIGN_POLOIDAL_ID )
                    tmpFace(:,:,0) = vaecrb(:,:,0,1,ispion(is,js))
                    tmpFace(:,:,1) = vaecrb(:,:,1,1,ispion(is,js))
                    call write_face_vector_component( edge_grid,              &
                      &   vectorComponent = profiles_ggd%                     &
                      &         ion( is )%state( js )%velocity_exb,           &
                      &   b2FaceData = tmpFace,                               &
                      &   vectorID = VEC_ALIGN_RADIAL_ID )
                  end if
                !! cvsa: Ion diffusivity
                  tmpFace(:,:,0) = cvsa(:,:,0,0,ispion(is,js))
                  tmpFace(:,:,1) = cvsa(:,:,1,1,ispion(is,js))
                  call write_face_vector_component( transport_grid,           &
                      &   vectorComponent = transport_ggd(1)%                 &
                      &         ion( is )%state( js )%momentum%d,             &
                      &   b2FaceData = tmpFace,                               &
                      &   vectorID = VEC_ALIGN_PARALLEL_ID )
#endif
                end do
#if ( IMAS_MAJOR_VERSION < 4 || ( IMAS_MAJOR_VERSION == 4 && IMAS_MINOR_VERSION < 1 ) )
                !! fmo: Ion momentum flux
                totFace(:,:,:) = 0.0_IDS_real
                do js = 1, istion(is)
                  hlpFace(:,:,:) = fmo(:,:,:,0,ispion(is,js)) +               &
                     &             fmo(:,:,:,1,ispion(is,js))
                  call divide_by_contact_areas(nx,ny,hlpFace,tmpFace)
                  call write_face_vector_component( transport_grid,           &
                      &   vectorComponent = edge_transport%model(1)%          &
                      &                     ggd( time_sind )%ion( is )%       &
                      &                     state( js )%momentum%flux,        &
                      &   b2FaceData = tmpFace,                               &
                      &   vectorID = VEC_ALIGN_PARALLEL_ID )
                  totFace(:,:,0) = totFace(:,:,0) + tmpFace(:,:,0)
                  totFace(:,:,1) = totFace(:,:,1) + tmpFace(:,:,1)
                end do
                call write_face_vector_component( transport_grid,             &
                    &   vectorComponent = edge_transport%model(1)%            &
                    &                     ggd( time_sind )%ion( is )%         &
                    &                     momentum%flux,                      &
                    &   b2FaceData = totFace,                                 &
                    &   vectorID = VEC_ALIGN_PARALLEL_ID )
                !! fllimvisc: Ion parallel momentum transport flux limit
                do js = 1, istion(is)
                  tmpFace(:,:,0) = fllimvisc(:,:,ispion(is,js))
                  tmpFace(:,:,1) = 1.0_IDS_real
                  call write_face_vector_component( transport_grid,           &
                      &   vectorComponent = edge_transport%model(1)%          &
                      &                     ggd( time_sind )%ion( is )%       &
                      &                     state( js )%momentum%flux_limiter,&
                      &   b2FaceData = tmpFace,                               &
                      &   vectorID = VEC_ALIGN_PARALLEL_ID )
                end do
                !! smo: Ion parallel momentum sources
                totCv(:,:) = 0.0_IDS_real
                do js = 1, istion(is)
                  do iy = -1, ny
                    do ix = -1, nx
                      tmpCv(ix,iy) = ( smo(ix,iy,0,ispion(is,js)) + &
                      &                smo(ix,iy,1,ispion(is,js)) * &
                      &                   ua(ix,iy,ispion(is,js)) + &
                      &                smo(ix,iy,2,ispion(is,js)) * &
                      &                 roxa(ix,iy,ispion(is,js)) + &
                      &                smo(ix,iy,3,ispion(is,js)) * &
                      &                 roxa(ix,iy,ispion(is,js)) * &
                      &                   ua(ix,iy,ispion(is,js)) ) / vol(ix,iy)
                    end do
                  end do
                  totCv(:,:) = totCv(:,:) + tmpCv(:,:)
                  call write_cell_vector_component( sources_grid,       &
                      &   vectorComponent = edge_sources%source(1)%     &
                      &                     ggd( time_sind )%ion( is )% &
                      &                     state( js )%momentum,       &
                      &   b2CellData = tmpCv,                           &
                      &   vectorID = VEC_ALIGN_PARALLEL_ID )
                end do
                call write_cell_vector_component( sources_grid,                &
                    &   vectorComponent = edge_sources%source(1)%              &
                    &                     ggd( time_sind )%ion( is )%momentum, &
                    &   b2CellData = totCv,                                    &
                    &   vectorID = VEC_ALIGN_PARALLEL_ID )
                totCv(:,:) = 0.0_IDS_real
                do js = 1, istion(is)
                  tmpCv(:,:) = ext_smo(:,:,ispion(is,js)) / vol(:,:)
                  totCv(:,:) = totCv(:,:) + tmpCv(:,:)
                  call write_cell_vector_component( sources_grid,       &
                      &   vectorComponent = edge_sources%source(2)%     &
                      &                     ggd( time_sind )%ion( is )% &
                      &                     state( js )%momentum,       &
                      &   b2CellData = tmpCv,                           &
                      &   vectorID = VEC_ALIGN_PARALLEL_ID )
                end do
                call write_cell_vector_component( sources_grid,         &
                      &   vectorComponent = edge_sources%source(2)%     &
                      &                     ggd( time_sind )%ion( is )% &
                      &                     momentum,                   &
                      &   b2CellData = totCv,                           &
                      &   vectorID = VEC_ALIGN_PARALLEL_ID )
                totCv(:,:) = 0.0_IDS_real
                do js = 1, istion(is)
                  tmpCv(:,:) = ( b2stbc_smo(:,:,ispion(is,js)) +        &
                        &        b2stbm_smo(:,:,ispion(is,js)) ) / vol(:,:)
                  totCv(:,:) = totCv(:,:) + tmpCv(:,:)
                  call write_cell_vector_component( sources_grid,       &
                      &   vectorComponent = edge_sources%source(3)%     &
                      &                     ggd( time_sind )%ion( is )% &
                      &                     state( js )%momentum,       &
                      &   b2CellData = tmpCv,                           &
                      &   vectorID = VEC_ALIGN_PARALLEL_ID )
                end do
                call write_cell_vector_component( sources_grid,         &
                      &   vectorComponent = edge_sources%source(3)%     &
                      &                     ggd( time_sind )%ion( is )% &
                      &                     momentum,                   &
                      &   b2CellData = totCv,                           &
                      &   vectorID = VEC_ALIGN_PARALLEL_ID )
                totCv(:,:) = 0.0_IDS_real
                do js = 1, istion(is)
                  do iy = -1, ny
                    do ix = -1, nx
                      tmpCv(ix,iy) = ( smodt(ix,iy,0,ispion(is,js)) + &
                      &                smodt(ix,iy,1,ispion(is,js)) * &
                      &                     ua(ix,iy,ispion(is,js)) + &
                      &                smodt(ix,iy,2,ispion(is,js)) * &
                      &                   roxa(ix,iy,ispion(is,js)) + &
                      &                smodt(ix,iy,3,ispion(is,js)) * &
                      &                   roxa(ix,iy,ispion(is,js)) * &
                      &                     ua(ix,iy,ispion(is,js)) ) / vol(ix,iy)
                    end do
                  end do
                  totCv(:,:) = totCv(:,:) + tmpCv(:,:)
                  call write_cell_vector_component( sources_grid,       &
                      &   vectorComponent = edge_sources%source(4)%     &
                      &                     ggd( time_sind )%ion( is )% &
                      &                     state( js )%momentum,       &
                      &   b2CellData = tmpCv,                           &
                      &   vectorID = VEC_ALIGN_PARALLEL_ID )
                end do
                call write_cell_vector_component( sources_grid,         &
                      &   vectorComponent = edge_sources%source(4)%     &
                      &                     ggd( time_sind )%ion( is )% &
                      &                     momentum,                   &
                      &   b2CellData = totCv,                           &
                      &   vectorID = VEC_ALIGN_PARALLEL_ID )
                if (use_eirene.ne.0 .and. balance_netcdf.ne.0) then
                  totCv = 0.0_IDS_real
                  do js = 1, istion(is)
                    tmpCv = 0.0_IDS_real
                    do istrai = 1, size( eirene_mc_mapl_smo_bal, 4)
                      tmpCv(:,:) = tmpCv(:,:) + &
                         &  eirene_mc_mapl_smo_bal(:,:,ispion(is,js),istrai)
                    end do
                    tmpCv(:,:) = tmpCv(:,:) / vol(:,:)
                    totCv(:,:) = totCv(:,:) + tmpCv(:,:)
                    call write_cell_vector_component( sources_grid,     &
                        &   vectorComponent = edge_sources%source(5)%   &
                        &            ggd( time_sind )%ion( is )%        &
                        &            state( js )%momentum,              &
                        &   b2CellData = tmpCv,                         &
                        &   vectorID = VEC_ALIGN_PARALLEL_ID )
                  end do
                  call write_cell_vector_component( sources_grid,       &
                      &   vectorComponent = edge_sources%source(5)%     &
                      &            ggd( time_sind )%ion( is )%momentum, &
                      &   b2CellData = totCv,                           &
                      &   vectorID = VEC_ALIGN_PARALLEL_ID )
                  totCv = 0.0_IDS_real
                  do js = 1, istion(is)
                    tmpCv = 0.0_IDS_real
                    do istrai = 1, size( eirene_mc_mmpl_smo_bal, 4)
                      tmpCv(:,:) = tmpCv(:,:) + &
                         &  eirene_mc_mmpl_smo_bal(:,:,ispion(is,js),istrai)
                    end do
                    tmpCv(:,:) = tmpCv(:,:) / vol(:,:)
                    totCv(:,:) = totCv(:,:) + tmpCv(:,:)
                    call write_cell_vector_component( sources_grid,     &
                        &   vectorComponent = edge_sources%source(6)%   &
                        &            ggd( time_sind )%ion( is )%        &
                        &            state( js )%momentum,              &
                        &   b2CellData = tmpCv,                         &
                        &   vectorID = VEC_ALIGN_PARALLEL_ID )
                  end do
                  call write_cell_vector_component( sources_grid,         &
                        &   vectorComponent = edge_sources%source(6)%     &
                        &            ggd( time_sind )%ion( is )%momentum, &
                        &   b2CellData = totCv,                           &
                        &   vectorID = VEC_ALIGN_PARALLEL_ID )
                else
                  totCv = 0.0_IDS_real
                  do js = 1, istion(is)
                    tmpCv(:,:) = b2stbr_smo(:,:,ispion(is,js)) / vol(:,:)
                    totCv(:,:) = totCv(:,:) + tmpCv(:,:)
                    call write_cell_vector_component( sources_grid,     &
                        &   vectorComponent = edge_sources%source(5)%   &
                        &            ggd( time_sind )%ion( is )%        &
                        &            state( js )%momentum,              &
                        &   b2CellData = tmpCv,                         &
                        &   vectorID = VEC_ALIGN_PARALLEL_ID )
                  end do
                  call write_cell_vector_component( sources_grid,         &
                        &   vectorComponent = edge_sources%source(5)%     &
                        &            ggd( time_sind )%ion( is )%momentum, &
                        &   b2CellData = totCv,                           &
                        &   vectorID = VEC_ALIGN_PARALLEL_ID )
                end if
                do js = 1, istion(is)
                  tmpCv(:,:) = rsamo(:,:,ispion(is,js)) / vol(:,:)
                  call write_cell_vector_component( sources_grid,       &
                      &   vectorComponent = edge_sources%source(7)%     &
                      &                     ggd( time_sind )%ion( is )% &
                      &                     state( js )%momentum,       &
                      &   b2CellData = tmpCv,                           &
                      &   vectorID = VEC_ALIGN_PARALLEL_ID )
                  tmpCv(:,:) = rramo(:,:,ispion(is,js)) / vol(:,:)
                  call write_cell_vector_component( sources_grid,       &
                      &   vectorComponent = edge_sources%source(8)%     &
                      &                     ggd( time_sind )%ion( is )% &
                      &                     state( js )%momentum,       &
                      &   b2CellData = tmpCv,                           &
                      &   vectorID = VEC_ALIGN_PARALLEL_ID )
                  tmpCv(:,:) = rcxmo(:,:,ispion(is,js)) / vol(:,:)
                  call write_cell_vector_component( sources_grid,       &
                      &   vectorComponent = edge_sources%source(9)%     &
                      &                     ggd( time_sind )%ion( is )% &
                      &                     state( js )%momentum,       &
                      &   b2CellData = tmpCv,                           &
                      &   vectorID = VEC_ALIGN_PARALLEL_ID )
                end do
#if ( ( IMAS_MINOR_VERSION > 38 || IMAS_MAJOR_VERSION > 3 ) && defined(B25_EIRENE) )
                if (use_eirene.ne.0) then
                  do js = 1, istion(is)
                    tmpCv(:,:) = smo0_eir_tot(:,:,ispion(is,js)) / vol(:,:)
                    call write_cell_vector_component( sources_grid,     &
                      &   vectorComponent = edge_sources%source(13)%    &
                      &                     ggd( time_sind )%ion( is )% &
                      &                     state( js )%momentum,       &
                      &   b2CellData = tmpCv,                           &
                      &   vectorID = VEC_ALIGN_PARALLEL_ID )
                    tmpCv(:,:) = smr0_eir_tot(:,:,ispion(is,js)) / vol(:,:)
                    call write_cell_vector_component( sources_grid,     &
                      &   vectorComponent = edge_sources%source(13)%    &
                      &                     ggd( time_sind )%ion( is )% &
                      &                     state( js )%momentum,       &
                      &   b2CellData = tmpCv,                           &
                      &   vectorID = VEC_ALIGN_RADIAL_ID )
                    tmpCv(:,:) = smd0_eir_tot(:,:,ispion(is,js)) / vol(:,:)
                    call write_cell_vector_component( sources_grid,     &
                      &   vectorComponent = edge_sources%source(13)%    &
                      &                     ggd( time_sind )%ion( is )% &
                      &                     state( js )%momentum,       &
                      &   b2CellData = tmpCv,                           &
                      &   vectorID = VEC_ALIGN_DIAMAGNETIC_ID )
                  end do
                end if
#endif
#else
                !! fmo: Ion momentum flux
                totFace(:,:,:) = 0.0_IDS_real
                do js = 1, istion(is)
                  hlpFace(:,:,:) = fmo(:,:,:,0,ispion(is,js)) +               &
                     &             fmo(:,:,:,1,ispion(is,js))
                  call divide_by_contact_areas(nx,ny,hlpFace,tmpFace)
                  call write_face_vector_component( transport_grid,           &
                      &   vectorComponent = transport_ggd(1)%ion( is )%       &
                      &                     state( js )%momentum%flux,        &
                      &   b2FaceData = tmpFace,                               &
                      &   vectorID = VEC_ALIGN_PARALLEL_ID )
                  totFace(:,:,0) = totFace(:,:,0) + tmpFace(:,:,0)
                  totFace(:,:,1) = totFace(:,:,1) + tmpFace(:,:,1)
                end do
                call write_face_vector_component( transport_grid,             &
                    &   vectorComponent = transport_ggd(1)%ion( is )%         &
                    &                     momentum%flux,                      &
                    &   b2FaceData = totFace,                                 &
                    &   vectorID = VEC_ALIGN_PARALLEL_ID )
                !! fllimvisc: Ion parallel momentum transport flux limit
                do js = 1, istion(is)
                  tmpFace(:,:,0) = fllimvisc(:,:,ispion(is,js))
                  tmpFace(:,:,1) = 1.0_IDS_real
                  call write_face_vector_component( transport_grid,           &
                      &   vectorComponent = transport_ggd(1)%ion( is )%       &
                      &                     state( js )%momentum%flux_limiter,&
                      &   b2FaceData = tmpFace,                               &
                      &   vectorID = VEC_ALIGN_PARALLEL_ID )
                end do
                !! smo: Ion parallel momentum sources
                totCv(:,:) = 0.0_IDS_real
                do js = 1, istion(is)
                  do iy = -1, ny
                    do ix = -1, nx
                      tmpCv(ix,iy) = ( smo(ix,iy,0,ispion(is,js)) + &
                      &                smo(ix,iy,1,ispion(is,js)) * &
                      &                   ua(ix,iy,ispion(is,js)) + &
                      &                smo(ix,iy,2,ispion(is,js)) * &
                      &                 roxa(ix,iy,ispion(is,js)) + &
                      &                smo(ix,iy,3,ispion(is,js)) * &
                      &                 roxa(ix,iy,ispion(is,js)) * &
                      &                   ua(ix,iy,ispion(is,js)) ) / vol(ix,iy)
                    end do
                  end do
                  totCv(:,:) = totCv(:,:) + tmpCv(:,:)
                  call write_cell_vector_component( sources_grid,            &
                      &   vectorComponent = sources_ggd(1)%ion( is )%        &
                      &                     state( js )%momentum,            &
                      &   b2CellData = tmpCv,                                &
                      &   vectorID = VEC_ALIGN_PARALLEL_ID )
                end do
                call write_cell_vector_component( sources_grid,              &
                    &   vectorComponent = sources_ggd(1)%ion( is )%momentum, &
                    &   b2CellData = totCv,                                  &
                    &   vectorID = VEC_ALIGN_PARALLEL_ID )
                totCv(:,:) = 0.0_IDS_real
                do js = 1, istion(is)
                  tmpCv(:,:) = ext_smo(:,:,ispion(is,js)) / vol(:,:)
                  totCv(:,:) = totCv(:,:) + tmpCv(:,:)
                  call write_cell_vector_component( sources_grid,            &
                      &   vectorComponent = sources_ggd(2)%ion( is )%        &
                      &                     state( js )%momentum,            &
                      &   b2CellData = tmpCv,                                &
                      &   vectorID = VEC_ALIGN_PARALLEL_ID )
                end do
                call write_cell_vector_component( sources_grid,              &
                    &   vectorComponent = sources_ggd(2)%ion( is )%momentum, &
                    &   b2CellData = totCv,                                  &
                    &   vectorID = VEC_ALIGN_PARALLEL_ID )
                totCv(:,:) = 0.0_IDS_real
                do js = 1, istion(is)
                  tmpCv(:,:) = ( b2stbc_smo(:,:,ispion(is,js)) +             &
                        &        b2stbm_smo(:,:,ispion(is,js)) ) / vol(:,:)
                  totCv(:,:) = totCv(:,:) + tmpCv(:,:)
                  call write_cell_vector_component( sources_grid,            &
                      &   vectorComponent = sources_ggd(3)%ion( is )%        &
                      &                     state( js )%momentum,            &
                      &   b2CellData = tmpCv,                                &
                      &   vectorID = VEC_ALIGN_PARALLEL_ID )
                end do
                call write_cell_vector_component( sources_grid,              &
                    &   vectorComponent = sources_ggd(3)%ion( is )%momentum, &
                    &   b2CellData = totCv,                                  &
                    &   vectorID = VEC_ALIGN_PARALLEL_ID )
                totCv(:,:) = 0.0_IDS_real
                do js = 1, istion(is)
                  do iy = -1, ny
                    do ix = -1, nx
                      tmpCv(ix,iy) = ( smodt(ix,iy,0,ispion(is,js)) + &
                      &                smodt(ix,iy,1,ispion(is,js)) * &
                      &                     ua(ix,iy,ispion(is,js)) + &
                      &                smodt(ix,iy,2,ispion(is,js)) * &
                      &                   roxa(ix,iy,ispion(is,js)) + &
                      &                smodt(ix,iy,3,ispion(is,js)) * &
                      &                   roxa(ix,iy,ispion(is,js)) * &
                      &                     ua(ix,iy,ispion(is,js)) ) / vol(ix,iy)
                    end do
                  end do
                  totCv(:,:) = totCv(:,:) + tmpCv(:,:)
                  call write_cell_vector_component( sources_grid,       &
                      &   vectorComponent = sources_ggd(4)%ion( is )%   &
                      &                     state( js )%momentum,       &
                      &   b2CellData = tmpCv,                           &
                      &   vectorID = VEC_ALIGN_PARALLEL_ID )
                end do
                call write_cell_vector_component( sources_grid,         &
                      &   vectorComponent = sources_ggd(4)%ion( is )%   &
                      &                     momentum,                   &
                      &   b2CellData = totCv,                           &
                      &   vectorID = VEC_ALIGN_PARALLEL_ID )
                if (use_eirene.ne.0 .and. balance_netcdf.ne.0) then
                  totCv = 0.0_IDS_real
                  do js = 1, istion(is)
                    tmpCv = 0.0_IDS_real
                    do istrai = 1, size( eirene_mc_mapl_smo_bal, 4)
                      tmpCv(:,:) = tmpCv(:,:) + &
                         &  eirene_mc_mapl_smo_bal(:,:,ispion(is,js),istrai)
                    end do
                    tmpCv(:,:) = tmpCv(:,:) / vol(:,:)
                    totCv(:,:) = totCv(:,:) + tmpCv(:,:)
                    call write_cell_vector_component( sources_grid,     &
                        &   vectorComponent = sources_ggd(5)%ion( is )% &
                        &            state( js )%momentum,              &
                        &   b2CellData = tmpCv,                         &
                        &   vectorID = VEC_ALIGN_PARALLEL_ID )
                  end do
                  call write_cell_vector_component( sources_grid,       &
                      &   vectorComponent = sources_ggd(5)%ion( is )%   &
                      &                     momentum,                   &
                      &   b2CellData = totCv,                           &
                      &   vectorID = VEC_ALIGN_PARALLEL_ID )
                  totCv = 0.0_IDS_real
                  do js = 1, istion(is)
                    tmpCv = 0.0_IDS_real
                    do istrai = 1, size( eirene_mc_mmpl_smo_bal, 4)
                      tmpCv(:,:) = tmpCv(:,:) + &
                         &  eirene_mc_mmpl_smo_bal(:,:,ispion(is,js),istrai)
                    end do
                    tmpCv(:,:) = tmpCv(:,:) / vol(:,:)
                    totCv(:,:) = totCv(:,:) + tmpCv(:,:)
                    call write_cell_vector_component( sources_grid,     &
                        &   vectorComponent = sources_ggd(6)%ion( is )% &
                        &            state( js )%momentum,              &
                        &   b2CellData = tmpCv,                         &
                        &   vectorID = VEC_ALIGN_PARALLEL_ID )
                  end do
                  call write_cell_vector_component( sources_grid,       &
                        &   vectorComponent = sources_ggd(6)%ion( is )% &
                        &                     momentum,                 &
                        &   b2CellData = totCv,                         &
                        &   vectorID = VEC_ALIGN_PARALLEL_ID )
                else
                  totCv = 0.0_IDS_real
                  do js = 1, istion(is)
                    tmpCv(:,:) = b2stbr_smo(:,:,ispion(is,js)) / vol(:,:)
                    totCv(:,:) = totCv(:,:) + tmpCv(:,:)
                    call write_cell_vector_component( sources_grid,     &
                        &   vectorComponent = sources_ggd(5)%ion( is )% &
                        &            state( js )%momentum,              &
                        &   b2CellData = tmpCv,                         &
                        &   vectorID = VEC_ALIGN_PARALLEL_ID )
                  end do
                  call write_cell_vector_component( sources_grid,       &
                        &   vectorComponent = sources_ggd(5)%ion( is )% &
                        &                     momentum,                 &
                        &   b2CellData = totCv,                         &
                        &   vectorID = VEC_ALIGN_PARALLEL_ID )
                end if
                do js = 1, istion(is)
                  tmpCv(:,:) = rsamo(:,:,ispion(is,js)) / vol(:,:)
                  call write_cell_vector_component( sources_grid,       &
                      &   vectorComponent = sources_ggd(7)%ion( is )%   &
                      &                     state( js )%momentum,       &
                      &   b2CellData = tmpCv,                           &
                      &   vectorID = VEC_ALIGN_PARALLEL_ID )
                  tmpCv(:,:) = rramo(:,:,ispion(is,js)) / vol(:,:)
                  call write_cell_vector_component( sources_grid,       &
                      &   vectorComponent = sources_ggd(8)%ion( is )%   &
                      &                     state( js )%momentum,       &
                      &   b2CellData = tmpCv,                           &
                      &   vectorID = VEC_ALIGN_PARALLEL_ID )
                  tmpCv(:,:) = rcxmo(:,:,ispion(is,js)) / vol(:,:)
                  call write_cell_vector_component( sources_grid,       &
                      &   vectorComponent = sources_ggd(9)%ion( is )%   &
                      &                     state( js )%momentum,       &
                      &   b2CellData = tmpCv,                           &
                      &   vectorID = VEC_ALIGN_PARALLEL_ID )
                end do
#ifdef B25_EIRENE
                if (use_eirene.ne.0) then
                  do js = 1, istion(is)
                    tmpCv(:,:) = smo0_eir_tot(:,:,ispion(is,js)) / vol(:,:)
                    call write_cell_vector_component( sources_grid,     &
                      &   vectorComponent = sources_ggd(13)%ion( is )%  &
                      &                     state( js )%momentum,       &
                      &   b2CellData = tmpCv,                           &
                      &   vectorID = VEC_ALIGN_PARALLEL_ID )
                    tmpCv(:,:) = smr0_eir_tot(:,:,ispion(is,js)) / vol(:,:)
                    call write_cell_vector_component( sources_grid,     &
                      &   vectorComponent = sources_ggd(13)%ion( is )%  &
                      &                     state( js )%momentum,       &
                      &   b2CellData = tmpCv,                           &
                      &   vectorID = VEC_ALIGN_RADIAL_ID )
                    tmpCv(:,:) = smd0_eir_tot(:,:,ispion(is,js)) / vol(:,:)
                    call write_cell_vector_component( sources_grid,     &
                      &   vectorComponent = sources_ggd(13)%ion( is )%  &
                      &                     state( js )%momentum,       &
                      &   b2CellData = tmpCv,                           &
                      &   vectorID = VEC_ALIGN_DIAMAGNETIC_ID )
                  end do
                end if
#endif
#endif
              end if
            end do

#if ( IMAS_MAJOR_VERSION < 4 || ( IMAS_MAJOR_VERSION == 4 && IMAS_MINOR_VERSION < 1 ) )
            !! te: Electron Temperature
            tmpCv(:,:) = te(:,:)/qe
            call write_quantity( edge_grid,                             &
                &   val = edge_profiles%ggd( time_sind )%electrons%     &
                &         temperature,                                  &
                &   value = tmpCv )
            tmpCv(:,:) = hce0(:,:)/ne(:,:)
            call write_quantity( transport_grid,                        &
                &   val = edge_transport%model(1)%ggd( time_sind )%     &
                &         electrons%energy%d,                           &
                &   value = tmpCv )
#if ( IMAS_MINOR_VERSION > 38 || IMAS_MAJOR_VERSION > 3 )
            call write_face_scalar( transport_grid,                     &
                &   val = edge_transport%model(1)%ggd( time_sind )%     &
                &         electrons%energy%v_pol,                       &
                &   value = chve(:,:,:,0) )
            call write_face_scalar( transport_grid,                     &
                &   val = edge_transport%model(1)%ggd( time_sind )%     &
                &         electrons%energy%v_radial,                    &
                &   value = chve(:,:,:,1) )
            call divide_by_contact_areas(nx,ny,fhe(:,:,:,0),flxFace(:,:,:,0))
            call divide_by_contact_areas(nx,ny,fhe(:,:,:,1),flxFace(:,:,:,1))
            call write_face_scalar( transport_grid,                     &
                &   val = edge_transport%model(1)%ggd( time_sind )%     &
                &         electrons%energy%flux_pol,                    &
                &   value = flxFace(:,:,:,0) )
            call write_face_scalar( transport_grid,                     &
                &   val = edge_transport%model(1)%ggd( time_sind )%     &
                &         electrons%energy%flux_radial,                 &
                &   value = flxFace(:,:,:,1) )
#else
            tmpFace(:,:,0) = chve(:,:,0,0)
            tmpFace(:,:,1) = chve(:,:,1,1)
            call write_face_scalar( transport_grid,                     &
                &   val = edge_transport%model(1)%ggd( time_sind )%     &
                &         electrons%energy%v,                           &
                &   value = tmpFace )
            tmpFace(:,:,:) = fhe(:,:,:,0) + fhe(:,:,:,1)
            call divide_by_contact_areas(nx,ny,tmpFace,totFace)
            call write_face_scalar( transport_grid,                     &
                &   val = edge_transport%model(1)%ggd( time_sind )%     &
                &         electrons%energy%flux,                        &
                &   value = totFace )
#endif
            call write_cell_scalar( transport_grid,                     &
                &   scalar = edge_transport%model(1)%ggd( time_sind )%  &
                &            electrons%energy%flux_limiter,             &
                &   b2CellData = fllime )
            tmpCv(:,:) = ( she(:,:,0) + she(:,:,1) * te(:,:) +          &
                &          she(:,:,2) * ne(:,:) +                       &
                &          she(:,:,3) * te(:,:) * ne(:,:) ) / vol(:,:)
            call write_cell_scalar( sources_grid,                       &
                &   scalar = edge_sources%source(1)%ggd( time_sind )%   &
                &            electrons%energy,                          &
                &   b2CellData = tmpCv )
            tmpCv(:,:) = ext_she(:,:) / vol(:,:)
            call write_cell_scalar( sources_grid,                       &
                &   scalar = edge_sources%source(2)%ggd( time_sind )%   &
                &            electrons%energy,                          &
                &   b2CellData = tmpCv )
            tmpCv(:,:) = ( b2stbc_she(:,:) + b2stbm_she(:,:) ) / vol(:,:)
            call write_cell_scalar( sources_grid,                       &
                &   scalar = edge_sources%source(3)%ggd( time_sind )%   &
                &            electrons%energy,                          &
                &   b2CellData = tmpCv )
            tmpCv(:,:) = ( shedt(:,:,0) + shedt(:,:,1) * te(:,:) +      &
                &          shedt(:,:,2) * ne(:,:) +                     &
                &          shedt(:,:,3) * te(:,:) * ne(:,:) ) / vol(:,:)
            call write_cell_scalar( sources_grid,                       &
                &   scalar = edge_sources%source(4)%ggd( time_sind )%   &
                &            electrons%energy,                          &
                &   b2CellData = tmpCv )
            if (use_eirene.ne.0 .and. balance_netcdf.ne.0) then
                tmpCv = 0.0_IDS_real
                do is = 1, size( eirene_mc_eael_she_bal, 3)
                    tmpCv(:,:) = tmpCv(:,:) + eirene_mc_eael_she_bal(:,:,is)
                end do
                tmpCv(:,:) = tmpCv(:,:) / vol(:,:)
                call write_cell_scalar( sources_grid,                       &
                    &   scalar = edge_sources%source(5)%ggd( time_sind )%   &
                    &            electrons%energy,                          &
                    &   b2CellData = tmpCv )
                tmpCv = 0.0_IDS_real
                do is = 1, size( eirene_mc_emel_she_bal, 3)
                    tmpCv(:,:) = tmpCv(:,:) + eirene_mc_emel_she_bal(:,:,is)
                end do
                tmpCv(:,:) = tmpCv(:,:) / vol(:,:)
                call write_cell_scalar( sources_grid,                       &
                    &   scalar = edge_sources%source(6)%ggd( time_sind )%   &
                    &            electrons%energy,                          &
                    &   b2CellData = tmpCv )
            else
                tmpCv(:,:) = b2stbr_she(:,:) / vol(:,:)
                call write_cell_scalar( sources_grid,                       &
                    &   scalar = edge_sources%source(5)%ggd( time_sind )%   &
                    &            electrons%energy,                          &
                    &   b2CellData = tmpCv )
            end if
            if (balance_netcdf.ne.0) then
                tmpCv(:,:) = b2stel_she_ion_bal(:,:) / vol(:,:)
                call write_cell_scalar( sources_grid,                       &
                    &   scalar = edge_sources%source(7)%ggd( time_sind )%   &
                    &            electrons%energy,                          &
                    &   b2CellData = tmpCv )
                tmpCv(:,:) = b2stel_she_rec_bal(:,:) / vol(:,:)
                call write_cell_scalar( sources_grid,                       &
                    &   scalar = edge_sources%source(8)%ggd( time_sind )%   &
                    &            electrons%energy,                          &
                    &   b2CellData = tmpCv )
                tmpCv(:,:) = -b2npht_shei_bal(:,:) / vol(:,:)
                call write_cell_scalar( sources_grid,                       &
                    &   scalar = edge_sources%source(10)%ggd( time_sind )%  &
                    &            electrons%energy,                          &
                    &   b2CellData = tmpCv )
            end if
            tmpCv(:,:) = b2sihs_joule(:,:) / vol(:,:)
            call write_cell_scalar( sources_grid,                           &
                &   scalar = edge_sources%source(11)%ggd( time_sind )%      &
                &            electrons%energy,                              &
                &   b2CellData = tmpCv )
            tmpCv(:,:) = 0.0_IDS_real
            do is = 0, ns-1
              tmpCv(:,:) = tmpCv(:,:) + rqrad(:,:,is)
            end do
#ifdef B25_EIRENE
            do is = 1, natmi
              tmpCv(-1:nx,-1:ny) = tmpCv(-1:nx,-1:ny) - eneutrad(0:nx+1,0:ny+1,is,0)
            end do
            do is = 1, nmoli
              tmpCv(-1:nx,-1:ny) = tmpCv(-1:nx,-1:ny) - emolrad(0:nx+1,0:ny+1,is,0)
            end do
            do is = 1, nioni
              tmpCv(-1:nx,-1:ny) = tmpCv(-1:nx,-1:ny) - eionrad(0:nx+1,0:ny+1,is,0)
            end do
#endif
            tmpCv(:,:) = tmpCv(:,:) / vol(:,:)
            call write_cell_scalar( sources_grid,                           &
                &   scalar = edge_sources%source(12)%ggd( time_sind )%      &
                &            electrons%energy,                              &
                &   b2CellData = tmpCv )
#if ( ( IMAS_MINOR_VERSION > 38 || IMAS_MAJOR_VERSION > 3 ) && defined(B25_EIRENE) )
            if (use_eirene.ne.0) then
              tmpCv(:,:) = she0_eir_tot(:,:) / vol(:,:)
              call write_cell_scalar( sources_grid,                         &
                &   scalar = edge_sources%source(13)%ggd( time_sind )%      &
                &            electrons%energy,                              &
                &   b2CellData = tmpCv )
            end if
#endif
#else
            !! te: Electron Temperature
            tmpCv(:,:) = te(:,:)/qe
            call write_quantity( edge_grid,                             &
                &   val = profiles_ggd%electrons%temperature,           &
                &   value = tmpCv )
            tmpCv(:,:) = hce0(:,:)/ne(:,:)
            call write_quantity( transport_grid,                        &
                &   val = transport_ggd(1)%electrons%energy%d,          &
                &   value = tmpCv )
            call write_face_scalar( transport_grid,                     &
                &   val = transport_ggd(1)%electrons%energy%v_pol,      &
                &   value = chve(:,:,:,0) )
            call write_face_scalar( transport_grid,                     &
                &   val = transport_ggd(1)%electrons%energy%v_radial,   &
                &   value = chve(:,:,:,1) )
            call divide_by_contact_areas(nx,ny,fhe(:,:,:,0),flxFace(:,:,:,0))
            call divide_by_contact_areas(nx,ny,fhe(:,:,:,1),flxFace(:,:,:,1))
            call write_face_scalar( transport_grid,                     &
                &   val = transport_ggd(1)%electrons%energy%flux_pol,   &
                &   value = flxFace(:,:,:,0) )
            call write_face_scalar( transport_grid,                     &
                &   val = transport_ggd(1)%electrons%energy%flux_radial, &
                &   value = flxFace(:,:,:,1) )
            call write_cell_scalar( transport_grid,                     &
                &   scalar = transport_ggd(1)%                          &
                &            electrons%energy%flux_limiter,             &
                &   b2CellData = fllime )
            tmpCv(:,:) = ( she(:,:,0) + she(:,:,1) * te(:,:) +          &
                &          she(:,:,2) * ne(:,:) +                       &
                &          she(:,:,3) * te(:,:) * ne(:,:) ) / vol(:,:)
            call write_cell_scalar( sources_grid,                       &
                &   scalar = sources_ggd(1)%electrons%energy,           &
                &   b2CellData = tmpCv )
            tmpCv(:,:) = ext_she(:,:) / vol(:,:)
            call write_cell_scalar( sources_grid,                       &
                &   scalar = sources_ggd(2)%electrons%energy,           &
                &   b2CellData = tmpCv )
            tmpCv(:,:) = ( b2stbc_she(:,:) + b2stbm_she(:,:) ) / vol(:,:)
            call write_cell_scalar( sources_grid,                       &
                &   scalar = sources_ggd(3)%electrons%energy,           &
                &   b2CellData = tmpCv )
            tmpCv(:,:) = ( shedt(:,:,0) + shedt(:,:,1) * te(:,:) +      &
                &          shedt(:,:,2) * ne(:,:) +                     &
                &          shedt(:,:,3) * te(:,:) * ne(:,:) ) / vol(:,:)
            call write_cell_scalar( sources_grid,                       &
                &   scalar = sources_ggd(4)%electrons%energy,           &
                &   b2CellData = tmpCv )
            if (use_eirene.ne.0 .and. balance_netcdf.ne.0) then
                tmpCv = 0.0_IDS_real
                do is = 1, size( eirene_mc_eael_she_bal, 3)
                    tmpCv(:,:) = tmpCv(:,:) + eirene_mc_eael_she_bal(:,:,is)
                end do
                tmpCv(:,:) = tmpCv(:,:) / vol(:,:)
                call write_cell_scalar( sources_grid,                   &
                    &   scalar = sources_ggd(5)%electrons%energy,       &
                    &   b2CellData = tmpCv )
                tmpCv = 0.0_IDS_real
                do is = 1, size( eirene_mc_emel_she_bal, 3)
                    tmpCv(:,:) = tmpCv(:,:) + eirene_mc_emel_she_bal(:,:,is)
                end do
                tmpCv(:,:) = tmpCv(:,:) / vol(:,:)
                call write_cell_scalar( sources_grid,                   &
                    &   scalar = sources_ggd(6)%electrons%energy,       &
                    &   b2CellData = tmpCv )
            else
                tmpCv(:,:) = b2stbr_she(:,:) / vol(:,:)
                call write_cell_scalar( sources_grid,                   &
                    &   scalar = sources_ggd(5)%electrons%energy,       &
                    &   b2CellData = tmpCv )
            end if
            if (balance_netcdf.ne.0) then
                tmpCv(:,:) = b2stel_she_ion_bal(:,:) / vol(:,:)
                call write_cell_scalar( sources_grid,                   &
                    &   scalar = sources_ggd(7)%electrons%energy,       &
                    &   b2CellData = tmpCv )
                tmpCv(:,:) = b2stel_she_rec_bal(:,:) / vol(:,:)
                call write_cell_scalar( sources_grid,                   &
                    &   scalar = sources_ggd(8)%electrons%energy,       &
                    &   b2CellData = tmpCv )
                tmpCv(:,:) = -b2npht_shei_bal(:,:) / vol(:,:)
                call write_cell_scalar( sources_grid,                   &
                    &   scalar = sources_ggd(10)%electrons%energy,      &
                    &   b2CellData = tmpCv )
            end if
            tmpCv(:,:) = b2sihs_joule(:,:) / vol(:,:)
            call write_cell_scalar( sources_grid,                       &
                &   scalar = sources_ggd(11)%electrons%energy,          &
                &   b2CellData = tmpCv )
            tmpCv(:,:) = 0.0_IDS_real
            do is = 0, ns-1
              tmpCv(:,:) = tmpCv(:,:) + rqrad(:,:,is)
            end do
#ifdef B25_EIRENE
            do is = 1, natmi
              tmpCv(-1:nx,-1:ny) = tmpCv(-1:nx,-1:ny) - eneutrad(0:nx+1,0:ny+1,is,0)
            end do
            do is = 1, nmoli
              tmpCv(-1:nx,-1:ny) = tmpCv(-1:nx,-1:ny) - emolrad(0:nx+1,0:ny+1,is,0)
            end do
            do is = 1, nioni
              tmpCv(-1:nx,-1:ny) = tmpCv(-1:nx,-1:ny) - eionrad(0:nx+1,0:ny+1,is,0)
            end do
#endif
            tmpCv(:,:) = tmpCv(:,:) / vol(:,:)
            call write_cell_scalar( sources_grid,                       &
                &   scalar = sources_ggd(12)%electrons%energy,          &
                &   b2CellData = tmpCv )
#ifdef B25_EIRENE
            if (use_eirene.ne.0) then
              tmpCv(:,:) = she0_eir_tot(:,:) / vol(:,:)
              call write_cell_scalar( sources_grid,                     &
                &   scalar = sources_ggd(13)%electrons%energy,          &
                &   b2CellData = tmpCv )
            end if
#endif
#endif
#if ( IMAS_MAJOR_VERSION < 4 || ( IMAS_MAJOR_VERSION == 4 && IMAS_MINOR_VERSION < 1 ) )
            !! pe: Electron pressure
            call b2xppe( nx, ny, ne, te, pe)
            call write_quantity( edge_grid,                           &
                &   val = edge_profiles%ggd( time_sind )%electrons%   &
                &         pressure,                                   &
                &   value = pe )

            !! ti: (Common) Ion Temperature
            tmpCv(:,:) = ti(:,:)/qe
            call write_quantity( edge_grid,                           &
                &   val = edge_profiles%ggd( time_sind )%t_i_average, &
                &   value = tmpCv )
            tmpCv(:,:) = hci0(:,:)/ni(:,:,0)
            call write_quantity( transport_grid,                      &
                &   val = edge_transport%model(1)%ggd( time_sind )%   &
                &         total_ion_energy%d,                         &
                &   value = tmpCv )
#if ( IMAS_MINOR_VERSION > 38 || IMAS_MAJOR_VERSION > 3 )
            call write_face_scalar( transport_grid,                   &
                 &   val = edge_transport%model(1)%ggd( time_sind )%  &
                 &         total_ion_energy%v_pol,                    &
                 &   value = chvi(:,:,:,0) )
            call write_face_scalar( transport_grid,                   &
                 &   val = edge_transport%model(1)%ggd( time_sind )%  &
                 &         total_ion_energy%v_radial,                 &
                 &   value = chvi(:,:,:,1) )
            !! fhi : Ion heat flux
            call divide_by_contact_areas(nx,ny,fhi(:,:,:,0),flxFace(:,:,:,0))
            call divide_by_contact_areas(nx,ny,fhi(:,:,:,1),flxFace(:,:,:,1))
            call write_face_scalar( transport_grid,                     &
                &   val = edge_transport%model(1)%ggd( time_sind )%     &
                &         total_ion_energy%flux_pol,                    &
                &   value = flxFace(:,:,:,0) )
            call write_face_scalar( transport_grid,                     &
                &   val = edge_transport%model(1)%ggd( time_sind )%     &
                &         total_ion_energy%flux_radial,                 &
                &   value = flxFace(:,:,:,1) )
#else
            tmpFace(:,:,0) = chvi(:,:,0,0)
            tmpFace(:,:,1) = chvi(:,:,1,1)
            call write_face_scalar( transport_grid,                   &
                 &   val = edge_transport%model(1)%ggd( time_sind )%  &
                 &         total_ion_energy%v,                        &
                 &   value = tmpFace )
            !! fhi : Ion heat flux
            tmpFace(:,:,:) = fhi(:,:,:,0) + fhi(:,:,:,1)
            call divide_by_contact_areas(nx,ny,tmpFace,totFace)
            call write_face_scalar( transport_grid,                     &
                &   val = edge_transport%model(1)%ggd( time_sind )%     &
                &         total_ion_energy%flux,                        &
                &   value = totFace )
#endif
            call write_cell_scalar( transport_grid,                     &
                &   scalar = edge_transport%model(1)%ggd( time_sind )%  &
                &            total_ion_energy%flux_limiter,             &
                &   b2CellData = fllimi )
            !! Ion energy sources
            tmpCv(:,:) = ( shi(:,:,0) + shi(:,:,1) * ti(:,:) +          &
                &          shi(:,:,2) * ni(:,:,0) +                     &
                &          shi(:,:,3) * ni(:,:,0) * ti(:,:) ) / vol(:,:)
            call write_cell_scalar( sources_grid,                       &
                &   scalar = edge_sources%source(1)%ggd( time_sind )%   &
                &            total_ion_energy,                          &
                &   b2CellData = tmpCv )
            tmpCv(:,:) = ext_shi(:,:) / vol(:,:)
            call write_cell_scalar( sources_grid,                       &
                &   scalar = edge_sources%source(2)%ggd( time_sind )%   &
                &            total_ion_energy,                          &
                &   b2CellData = tmpCv )
            tmpCv(:,:) = ( b2stbc_shi(:,:) + b2stbm_shi(:,:) ) / vol(:,:)
            call write_cell_scalar( sources_grid,                       &
                &   scalar = edge_sources%source(3)%ggd( time_sind )%   &
                &            total_ion_energy,                          &
                &   b2CellData = tmpCv )
            tmpCv(:,:) = ( shidt(:,:,0) + shidt(:,:,1) * ti(:,:) +      &
                &          shidt(:,:,2) * ni(:,:,0) +                   &
                &          shidt(:,:,3) * ni(:,:,0) * ti(:,:) ) / vol(:,:)
            call write_cell_scalar( sources_grid,                       &
                &   scalar = edge_sources%source(4)%ggd( time_sind )%   &
                &            total_ion_energy,                          &
                &   b2CellData = tmpCv )
            if (use_eirene.ne.0 .and. balance_netcdf.ne.0) then
              tmpCv = 0.0_IDS_real
              do is = 1, size( eirene_mc_eapl_shi_bal, 3)
                tmpCv(:,:) = tmpCv(:,:) + eirene_mc_eapl_shi_bal(:,:,is)
              end do
              tmpCv(:,:) = tmpCv(:,:) / vol(:,:)
              call write_cell_scalar( sources_grid,                       &
                  &   scalar = edge_sources%source(5)%ggd( time_sind )%   &
                  &            total_ion_energy,                          &
                  &   b2CellData = tmpCv )
              tmpCv = 0.0_IDS_real
              do is = 1, size( eirene_mc_empl_shi_bal, 3)
                tmpCv(:,:) = tmpCv(:,:) + eirene_mc_empl_shi_bal(:,:,is)
              end do
              tmpCv(:,:) = tmpCv(:,:) / vol(:,:)
              call write_cell_scalar( sources_grid,                       &
                  &   scalar = edge_sources%source(6)%ggd( time_sind )%   &
                  &            total_ion_energy,                          &
                  &   b2CellData = tmpCv )
            else
              tmpCv(:,:) = b2stbr_shi(:,:) / vol(:,:)
              call write_cell_scalar( sources_grid,                       &
                  &   scalar = edge_sources%source(5)%ggd( time_sind )%   &
                  &            total_ion_energy,                          &
                  &   b2CellData = tmpCv )
            end if
            if (balance_netcdf.ne.0) then
              tmpCv(:,:) = b2stel_shi_ion_bal(:,:) / vol(:,:)
              call write_cell_scalar( sources_grid,                       &
                  &   scalar = edge_sources%source(7)%ggd( time_sind )%   &
                  &            total_ion_energy,                          &
                  &   b2CellData = tmpCv )
              tmpCv(:,:) = b2stel_shi_rec_bal(:,:) / vol(:,:)
              call write_cell_scalar( sources_grid,                       &
                  &   scalar = edge_sources%source(8)%ggd( time_sind )%   &
                  &            total_ion_energy,                          &
                  &   b2CellData = tmpCv )
              tmpCv(:,:) = b2npht_shei_bal(:,:) / vol(:,:)
              call write_cell_scalar( sources_grid,                       &
                  &   scalar = edge_sources%source(10)%ggd( time_sind )%  &
                  &            total_ion_energy,                          &
                  &   b2CellData = tmpCv )
            end if
#if ( ( IMAS_MINOR_VERSION > 38 || IMAS_MAJOR_VERSION > 3 ) && defined(B25_EIRENE) )
            if (use_eirene.ne.0) then
              tmpCv(:,:) = shi0_eir_tot(:,:) / vol(:,:)
              call write_cell_scalar( sources_grid,                       &
                  &   scalar = edge_sources%source(13)%ggd( time_sind )%  &
                  &            total_ion_energy,                          &
                  &   b2CellData = tmpCv )
            end if
#endif
#else
            !! pe: Electron pressure
            call b2xppe( nx, ny, ne, te, pe)
            call write_quantity( edge_grid,                           &
                &   val = profiles_ggd%electrons%pressure,            &
                &   value = pe )

            !! ti: (Common) Ion Temperature
            tmpCv(:,:) = ti(:,:)/qe
            call write_quantity( edge_grid,                           &
                &   val = profiles_ggd%t_i_average,                   &
                &   value = tmpCv )
            tmpCv(:,:) = hci0(:,:)/ni(:,:,0)
            call write_quantity( transport_grid,                      &
                &   val = transport_ggd(1)%total_ion_energy%d,        &
                &   value = tmpCv )
            call write_face_scalar( transport_grid,                   &
                 &   val = transport_ggd(1)%total_ion_energy%v_pol,   &
                 &   value = chvi(:,:,:,0) )
            call write_face_scalar( transport_grid,                   &
                 &   val = transport_ggd(1)%total_ion_energy%v_radial, &
                 &   value = chvi(:,:,:,1) )
            !! fhi : Ion heat flux
            call divide_by_contact_areas(nx,ny,fhi(:,:,:,0),flxFace(:,:,:,0))
            call divide_by_contact_areas(nx,ny,fhi(:,:,:,1),flxFace(:,:,:,1))
            call write_face_scalar( transport_grid,                   &
                &   val = transport_ggd(1)%total_ion_energy%flux_pol, &
                &   value = flxFace(:,:,:,0) )
            call write_face_scalar( transport_grid,                   &
                &   val = transport_ggd(1)%total_ion_energy%flux_radial, &
                &   value = flxFace(:,:,:,1) )
            call write_cell_scalar( transport_grid,                   &
                &   scalar = transport_ggd(1)%total_ion_energy%       &
                &            flux_limiter,                            &
                &   b2CellData = fllimi )
            !! Ion energy sources
            tmpCv(:,:) = ( shi(:,:,0) + shi(:,:,1) * ti(:,:) +        &
                &          shi(:,:,2) * ni(:,:,0) +                   &
                &          shi(:,:,3) * ni(:,:,0) * ti(:,:) ) / vol(:,:)
            call write_cell_scalar( sources_grid,                     &
                &   scalar = sources_ggd(1)%total_ion_energy,         &
                &   b2CellData = tmpCv )
            tmpCv(:,:) = ext_shi(:,:) / vol(:,:)
            call write_cell_scalar( sources_grid,                     &
                &   scalar = sources_ggd(2)%total_ion_energy,         &
                &   b2CellData = tmpCv )
            tmpCv(:,:) = ( b2stbc_shi(:,:) + b2stbm_shi(:,:) ) / vol(:,:)
            call write_cell_scalar( sources_grid,                     &
                &   scalar = sources_ggd(3)%total_ion_energy,         &
                &   b2CellData = tmpCv )
            tmpCv(:,:) = ( shidt(:,:,0) + shidt(:,:,1) * ti(:,:) +    &
                &          shidt(:,:,2) * ni(:,:,0) +                 &
                &          shidt(:,:,3) * ni(:,:,0) * ti(:,:) ) / vol(:,:)
            call write_cell_scalar( sources_grid,                     &
                &   scalar = sources_ggd(4)%total_ion_energy,         &
                &   b2CellData = tmpCv )
            if (use_eirene.ne.0 .and. balance_netcdf.ne.0) then
              tmpCv = 0.0_IDS_real
              do is = 1, size( eirene_mc_eapl_shi_bal, 3)
                tmpCv(:,:) = tmpCv(:,:) + eirene_mc_eapl_shi_bal(:,:,is)
              end do
              tmpCv(:,:) = tmpCv(:,:) / vol(:,:)
              call write_cell_scalar( sources_grid,                   &
                  &   scalar = sources_ggd(5)%total_ion_energy,       &
                  &   b2CellData = tmpCv )
              tmpCv = 0.0_IDS_real
              do is = 1, size( eirene_mc_empl_shi_bal, 3)
                tmpCv(:,:) = tmpCv(:,:) + eirene_mc_empl_shi_bal(:,:,is)
              end do
              tmpCv(:,:) = tmpCv(:,:) / vol(:,:)
              call write_cell_scalar( sources_grid,                   &
                  &   scalar = sources_ggd(6)%total_ion_energy,       &
                  &   b2CellData = tmpCv )
            else
              tmpCv(:,:) = b2stbr_shi(:,:) / vol(:,:)
              call write_cell_scalar( sources_grid,                   &
                  &   scalar = sources_ggd(5)%total_ion_energy,       &
                  &   b2CellData = tmpCv )
            end if
            if (balance_netcdf.ne.0) then
              tmpCv(:,:) = b2stel_shi_ion_bal(:,:) / vol(:,:)
              call write_cell_scalar( sources_grid,                   &
                  &   scalar = sources_ggd(7)%total_ion_energy,       &
                  &   b2CellData = tmpCv )
              tmpCv(:,:) = b2stel_shi_rec_bal(:,:) / vol(:,:)
              call write_cell_scalar( sources_grid,                   &
                  &   scalar = sources_ggd(8)%total_ion_energy,       &
                  &   b2CellData = tmpCv )
              tmpCv(:,:) = b2npht_shei_bal(:,:) / vol(:,:)
              call write_cell_scalar( sources_grid,                   &
                  &   scalar = sources_ggd(10)%total_ion_energy,      &
                  &   b2CellData = tmpCv )
            end if
#ifdef B25_EIRENE
            if (use_eirene.ne.0) then
              tmpCv(:,:) = shi0_eir_tot(:,:) / vol(:,:)
              call write_cell_scalar( sources_grid,                   &
                  &   scalar = sources_ggd(13)%total_ion_energy,      &
                  &   b2CellData = tmpCv )
            end if
#endif
#endif
            do is = 1, nsion
              if (is.le.nspecies) then
#if ( IMAS_MAJOR_VERSION < 4 || ( IMAS_MAJOR_VERSION == 4 && IMAS_MINOR_VERSION < 1 ) )
                totCv(:,:) = 0.0_IDS_real
                do js = 1, istion(is)
            !! Ion energy sources resolved by species
                  tmpCv(:,:) = rsahi(:,:,ispion(is,js)) / vol(:,:)
                  totCv(:,:) = totCv(:,:) + tmpCv(:,:)
                  call write_cell_scalar( sources_grid,                       &
                      &   scalar = edge_sources%source(7)%ggd( time_sind )%   &
                      &            ion( is )%state( js )%energy,              &
                      &   b2CellData = tmpCv )
                end do
                call write_cell_scalar( sources_grid,                         &
                      &   scalar = edge_sources%source(7)%ggd( time_sind )%   &
                      &            ion( is )%energy,                          &
                      &   b2CellData = totCv )
                totCv(:,:) = 0.0_IDS_real
                do js = 1, istion(is)
                  tmpCv(:,:) = rrahi(:,:,ispion(is,js)) / vol(:,:)
                  totCv(:,:) = totCv(:,:) + tmpCv(:,:)
                  call write_cell_scalar( sources_grid,                       &
                      &   scalar = edge_sources%source(8)%ggd( time_sind )%   &
                      &            ion( is )%state( js )%energy,              &
                      &   b2CellData = tmpCv )
                end do
                call write_cell_scalar( sources_grid,                         &
                    &   scalar = edge_sources%source(8)%ggd( time_sind )%     &
                    &            ion( is )%energy,                            &
                        &   b2CellData = totCv )
                totCv(:,:) = 0.0_IDS_real
                do js = 1, istion(is)
                  tmpCv(:,:) = rcxhi(:,:,ispion(is,js)) / vol(:,:)
                  totCv(:,:) = totCv(:,:) + tmpCv(:,:)
                  call write_cell_scalar( sources_grid,                       &
                      &   scalar = edge_sources%source(9)%ggd( time_sind )%   &
                      &            ion( is )%state( js )%energy,              &
                      &   b2CellData = tmpCv )
                end do
                call write_cell_scalar( sources_grid,                         &
                    &   scalar = edge_sources%source(9)%ggd( time_sind )%     &
                    &            ion( is )%energy,                            &
                    &   b2CellData = totCv )
#else
                totCv(:,:) = 0.0_IDS_real
                do js = 1, istion(is)
            !! Ion energy sources resolved by species
                  tmpCv(:,:) = rsahi(:,:,ispion(is,js)) / vol(:,:)
                  totCv(:,:) = totCv(:,:) + tmpCv(:,:)
                  call write_cell_scalar( sources_grid,               &
                      &   scalar = sources_ggd(7)%                    &
                      &            ion( is )%state( js )%energy,      &
                      &   b2CellData = tmpCv )
                end do
                call write_cell_scalar( sources_grid,                 &
                      &   scalar = sources_ggd(7)%ion( is )%energy,   &
                      &   b2CellData = totCv )
                totCv(:,:) = 0.0_IDS_real
                do js = 1, istion(is)
                  tmpCv(:,:) = rrahi(:,:,ispion(is,js)) / vol(:,:)
                  totCv(:,:) = totCv(:,:) + tmpCv(:,:)
                  call write_cell_scalar( sources_grid,               &
                      &   scalar = sources_ggd(8)%                    &
                      &            ion( is )%state( js )%energy,      &
                      &   b2CellData = tmpCv )
                end do
                call write_cell_scalar( sources_grid,                 &
                    &   scalar = sources_ggd(8)%ion( is )%energy,     &
                        &   b2CellData = totCv )
                totCv(:,:) = 0.0_IDS_real
                do js = 1, istion(is)
                  tmpCv(:,:) = rcxhi(:,:,ispion(is,js)) / vol(:,:)
                  totCv(:,:) = totCv(:,:) + tmpCv(:,:)
                  call write_cell_scalar( sources_grid,               &
                      &   scalar = sources_ggd(9)%                    &
                      &            ion( is )%state( js )%energy,      &
                      &   b2CellData = tmpCv )
                end do
                call write_cell_scalar( sources_grid,                 &
                    &   scalar = sources_ggd(9)%ion( is )%energy,     &
                    &   b2CellData = totCv )
#endif
              else
#ifdef B25_EIRENE
                do js = 1, istion(is)
            !! Test ion temperature
                  tmpCv(-1:nx,-1:ny) = tib2(0:nx+1,0:ny+1,ispion(is,js),1)/qe
#if ( IMAS_MAJOR_VERSION < 4 || ( IMAS_MAJOR_VERSION == 4 && IMAS_MINOR_VERSION < 1 ) )
                  call write_quantity( edge_grid,                               &
                      &   val = edge_profiles%ggd( time_sind )%ion( is )%       &
                      &         state( js )%temperature,                        &
                      &   value = tmpCv )
#else
                  call write_quantity( edge_grid,                               &
                      &   val = profiles_ggd%ion( is )%state( js )%temperature, &
                      &   value = tmpCv )
#endif
                end do
#endif
              end if
            end do

            do is = 1, nsion
              if (is.le.nspecies) then
#if ( IMAS_MAJOR_VERSION < 4 || ( IMAS_MAJOR_VERSION == 4 && IMAS_MINOR_VERSION < 1 ) )
                !! pb : Ion pressure
                totCv(:,:) = 0.0_IDS_real
                do js = 1, istion(is)
                  call b2xppb( nx, ny, rza(:,:,ispion(is,js)),              &
                      &                 na(:,:,ispion(is,js)), te, ti, pb)
                  totCv(:,:) = totCv(:,:) + pb(:,:)
                  call write_quantity( edge_grid,                           &
                      &   val = edge_profiles%ggd( time_sind )%ion( is )%   &
                      &         state( js )%pressure,                       &
                      &   value = pb )
                end do
                call write_quantity( edge_grid,                             &
                    &   val = edge_profiles%ggd( time_sind )%ion( is )%     &
                    &         pressure,                                     &
                    &   value = totCv )
                !! Kinetic energy density
                totCv(:,:) = 0.0_IDS_real
                do js = 1, istion(is)
                  tmpCv(:,:) = 0.5_IDS_real*am(ispion(is,js))*mp*           &
                      &       (ua(:,:,ispion(is,js))**2)*na(:,:,ispion(is,js))
                  totCv(:,:) = totCv(:,:) + tmpCv(:,:)
                  call write_quantity( edge_grid,                           &
                      &   val = edge_profiles%ggd( time_sind )%ion( is )%   &
                      &         state( js )%energy_density_kinetic,         &
                      &   value = tmpCv )
                end do
                call write_quantity( edge_grid,                             &
                    &   val = edge_profiles%ggd( time_sind )%ion( is )%     &
                    &         energy_density_kinetic,                       &
                    &   value = totCv )
                do js = 1, istion(is)
                !! Average charge
                  call write_quantity( edge_grid,                           &
                      &   val = edge_profiles%ggd( time_sind )%ion( is )%   &
                      &         state( js )%z_average,                      &
                      &   value = rza(:,:,ispion(is,js)) )
                !! Average square charge
                  call write_quantity( edge_grid,                           &
                      &   val = edge_profiles%ggd( time_sind )%ion( is )%   &
                      &         state( js )%z_square_average,               &
                      &   value = rz2(:,:,ispion(is,js)) )
#if IMAS_MAJOR_VERSION > 3
                !! Ionization potential
                  call write_quantity( edge_grid,                           &
                      &   val = edge_profiles%ggd( time_sind )%ion( is )%   &
                      &         state( js )%ionization_potential,           &
                      &   value = rpt(:,:,ispion(is,js)) )
#else
                !! Ionisation potential
                  call write_quantity( edge_grid,                           &
                      &   val = edge_profiles%ggd( time_sind )%ion( is )%   &
                      &         state( js )%ionisation_potential,           &
                      &   value = rpt(:,:,ispion(is,js)) )
#endif
#else
                !! pb : Ion pressure
                totCv(:,:) = 0.0_IDS_real
                do js = 1, istion(is)
                  call b2xppb( nx, ny, rza(:,:,ispion(is,js)),                &
                      &                 na(:,:,ispion(is,js)), te, ti, pb)
                  totCv(:,:) = totCv(:,:) + pb(:,:)
                  call write_quantity( edge_grid,                             &
                      &   val = profiles_ggd%ion( is )%state( js )%pressure,  &
                      &   value = pb )
                end do
                call write_quantity( edge_grid,                               &
                    &   val = profiles_ggd%ion( is )%pressure,                &
                    &   value = totCv )
                !! Kinetic energy density
                totCv(:,:) = 0.0_IDS_real
                do js = 1, istion(is)
                  tmpCv(:,:) = 0.5_IDS_real*am(ispion(is,js))*mp*             &
                      &       (ua(:,:,ispion(is,js))**2)*na(:,:,ispion(is,js))
                  totCv(:,:) = totCv(:,:) + tmpCv(:,:)
                  call write_quantity( edge_grid,                             &
                      &   val = profiles_ggd%ion( is )%                       &
                      &         state( js )%energy_density_kinetic,           &
                      &   value = tmpCv )
                end do
                call write_quantity( edge_grid,                               &
                    &   val = profiles_ggd%ion( is )%energy_density_kinetic,  &
                    &   value = totCv )
                do js = 1, istion(is)
                !! Average charge
                  call write_quantity( edge_grid,                             &
                      &   val = profiles_ggd%ion( is )%state( js )%z_average, &
                      &   value = rza(:,:,ispion(is,js)) )
                !! Average square charge
                  call write_quantity( edge_grid,                             &
                      &   val = profiles_ggd%ion( is )%state( js )%           &
                      &         z_square_average,                             &
                      &   value = rz2(:,:,ispion(is,js)) )
                !! Ionization potential
#if ( IMAS_MAJOR_VERSION == 4 && IMAS_MINOR_VERSION == 1 )
                  call write_quantity( edge_grid,                             &
                      &   val = profiles_ggd%ion( is )%                       &
                      &         state( js )%ionisation_potential,             &
                      &   value = rpt(:,:,ispion(is,js)) )
#else
                  call write_quantity( edge_grid,                             &
                      &   val = profiles_ggd%ion( is )%                       &
                      &         state( js )%ionization_potential,             &
                      &   value = rpt(:,:,ispion(is,js)) )
#endif
#endif
                end do
#ifdef B25_EIRENE
              else
#if ( IMAS_MAJOR_VERSION < 4 || ( IMAS_MAJOR_VERSION == 4 && IMAS_MINOR_VERSION < 1 ) )
                !! Test ion pressure
                totCv(:,:) = 0.0_IDS_real
                do js = 1, istion(is)
                  tmpCv(-1:nx,-1:ny) = dib2(0:nx+1,0:ny+1,ispion(is,js),1)  &
                      &               *tib2(0:nx+1,0:ny+1,ispion(is,js),1)
                  totCv(-1:nx,-1:ny) = totCv(-1:nx,-1:ny) + tmpCv(-1:nx,-1:ny)
                  call write_quantity( edge_grid,                           &
                      &   val = edge_profiles%ggd( time_sind )%ion( is )%   &
                      &         state( js )%pressure,                       &
                      &   value = tmpCv )
                end do
                call write_quantity( edge_grid,                             &
                    &   val = edge_profiles%ggd( time_sind )%ion( is )%     &
                    &         pressure,                                     &
                    &   value = totCv )
                !! Average charge
                do js = 1, istion(is)
                  tmpCv(:,:) = nchrgi( ispion(is,js) )
                  call write_quantity( edge_grid,                           &
                      &   val = edge_profiles%ggd( time_sind )%ion( is )%   &
                      &         state( js )%z_average,                      &
                      &   value = tmpCv )
                !! Average square charge
                  tmpCv(:,:) = nchrgi( js )**2
                  call write_quantity( edge_grid,                           &
                      &   val = edge_profiles%ggd( time_sind )%ion( is )%   &
                      &         state( js )%z_square_average,               &
                      &   value = tmpCv )
                end do
                !! Radiation source
                totCv = 0.0_IDS_real
                do js = 1, istion(is)
                  tmpCv(-1:nx,-1:ny) = eionrad(0:nx+1,0:ny+1,ispion(is,js),0) / vol(-1:nx,-1:ny)
                  totCv(-1:nx,-1:ny) = totCv(-1:nx,-1:ny) + tmpCv(-1:nx,-1:ny)
                  call write_cell_scalar( sources_grid,                     &
                      &   scalar = edge_sources%source(12)%                 &
                      &            ggd( time_sind )%ion( is )%              &
                      &            state( js )%energy,                      &
                      &   b2CellData = tmpCv )
                end do
                call write_cell_scalar( sources_grid,                       &
                    &   scalar = edge_sources%source(12)%                   &
                    &            ggd( time_sind )%ion( is )%energy,         &
                    &   b2CellData = totCv )
#else
                !! Test ion pressure
                totCv(:,:) = 0.0_IDS_real
                do js = 1, istion(is)
                  tmpCv(-1:nx,-1:ny) = dib2(0:nx+1,0:ny+1,ispion(is,js),1)    &
                      &               *tib2(0:nx+1,0:ny+1,ispion(is,js),1)
                  totCv(-1:nx,-1:ny) = totCv(-1:nx,-1:ny) + tmpCv(-1:nx,-1:ny)
                  call write_quantity( edge_grid,                             &
                      &   val = profiles_ggd%ion( is )%state( js )%pressure,  &
                      &   value = tmpCv )
                end do
                call write_quantity( edge_grid,                               &
                    &   val = profiles_ggd%ion( is )%pressure,                &
                    &   value = totCv )
                !! Average charge
                do js = 1, istion(is)
                  tmpCv(:,:) = nchrgi( ispion(is,js) )
                  call write_quantity( edge_grid,                             &
                      &   val = profiles_ggd%ion( is )%state( js )%z_average, &
                      &   value = tmpCv )
                !! Average square charge
                  tmpCv(:,:) = nchrgi( js )**2
                  call write_quantity( edge_grid,                             &
                      &   val = profiles_ggd%ion( is )%state( js )%           &
                      &         z_square_average,                             &
                      &   value = tmpCv )
                end do
                !! Radiation source
                totCv = 0.0_IDS_real
                do js = 1, istion(is)
                  tmpCv(-1:nx,-1:ny) = eionrad(0:nx+1,0:ny+1,ispion(is,js),0) / vol(-1:nx,-1:ny)
                  totCv(-1:nx,-1:ny) = totCv(-1:nx,-1:ny) + tmpCv(-1:nx,-1:ny)
                  call write_cell_scalar( sources_grid,                      &
                      &   scalar = sources_ggd(12)%ion( is )%                &
                      &            state( js )%energy,                       &
                      &   b2CellData = tmpCv )
                end do
                call write_cell_scalar( sources_grid,                        &
                    &   scalar = sources_ggd(12)%ion( is )%energy,           &
                    &   b2CellData = totCv )
#endif
#endif
              end if
            end do

#if ( IMAS_MAJOR_VERSION < 4 || ( IMAS_MAJOR_VERSION == 4 && IMAS_MINOR_VERSION < 1 ) )
            !! ni/ne: Total ion density over electron density
            tmpCv(:,:) = ni(:,:,1)/ne(:,:)
            call write_quantity( edge_grid,                                  &
                &   val = edge_profiles%ggd( time_sind )%n_i_total_over_n_e, &
                &   value = tmpCv )

            !! Zeff
            call write_quantity( edge_grid,                                  &
                &   val = edge_profiles%ggd( time_sind )%zeff,               &
                &   value = zeff )

            !! pz: Thermal plasma pressure (electrons+ions)
            call write_quantity( edge_grid,                                  &
                &   val = edge_profiles%ggd( time_sind )%pressure_thermal,   &
                &   value = pz )

#if ( IMAS_MINOR_VERSION > 32 || IMAS_MAJOR_VERSION > 3 )
            !! fch: Total current
            call divide_by_contact_areas(nx,ny,fch,tmpFace)
            totFace(:,:,0) = tmpFace(:,:,0)
            totFace(:,:,1) = IDS_REAL_INVALID
            call write_face_vector_component( edge_grid,                     &
                &   vectorComponent = edge_profiles%ggd( time_sind )%        &
                &                     j_total,                               &
                &   b2FaceData = totFace,                                    &
                &   vectorID = VEC_ALIGN_POLOIDAL_ID )
            totFace(:,:,0) = IDS_REAL_INVALID
            totFace(:,:,1) = tmpFace(:,:,1)
            call write_face_vector_component( edge_grid,                     &
                &   vectorComponent = edge_profiles%ggd( time_sind )%        &
                &                     j_total,                               &
                &   b2FaceData = totFace,                                    &
                &   vectorID = VEC_ALIGN_RADIAL_ID )

            !! fch_p: Parallel current
            hlpFace(:,:,0) = fch_p(:,:,0,0)
            hlpFace(:,:,1) = fch_p(:,:,1,1)
            call divide_by_poloidal_areas(nx,ny,hlpFace,tmpFace)
            call write_face_scalar( edge_grid,                               &
                &   val = edge_profiles%ggd( time_sind )%j_parallel,         &
                &   value = tmpFace )
#endif

            !! fchanml: Anomalous current
            call b2tanml (nx, ny, ns, ismain, vol, hx, hy, qs,               & !som 02.11.21
                &         csig_an, po, ne, na, fchanml_a, fchanml)
            hlpFace(:,:,0) = fchanml(:,:,0,0)
            hlpFace(:,:,1) = fchanml(:,:,1,0)
            call divide_by_contact_areas(nx,ny,hlpFace,totFace)
            call write_face_vector_component( edge_grid,                     &
                &   vectorComponent = edge_profiles%ggd( time_sind )%        &
                &                     j_anomalous,                           &
                &   b2FaceData = totFace,                                    &
                &   vectorID = VEC_ALIGN_POLOIDAL_ID )
            hlpFace(:,:,0) = fchanml(:,:,0,1)
            hlpFace(:,:,1) = fchanml(:,:,1,1)
            call divide_by_contact_areas(nx,ny,hlpFace,totFace)
            call write_face_vector_component( edge_grid,                     &
                &   vectorComponent = edge_profiles%ggd( time_sind )%        &
                &                     j_anomalous,                           &
                &   b2FaceData = totFace,                                    &
                &   vectorID = VEC_ALIGN_RADIAL_ID )

            !! fchinert: Inertial current
            hlpFace(:,:,0) = fchinert(:,:,0,0)
            hlpFace(:,:,1) = fchinert(:,:,1,0)
            call divide_by_contact_areas(nx,ny,hlpFace,totFace)
            call write_face_vector_component( edge_grid,                     &
                &   vectorComponent = edge_profiles%ggd( time_sind )%        &
                &                     j_inertial,                            &
                &   b2FaceData = totFace,                                    &
                &   vectorID = VEC_ALIGN_POLOIDAL_ID )
            hlpFace(:,:,0) = fchinert(:,:,0,1)
            hlpFace(:,:,1) = fchinert(:,:,1,1)
            call divide_by_contact_areas(nx,ny,hlpFace,totFace)
            call write_face_vector_component( edge_grid,                     &
                &   vectorComponent = edge_profiles%ggd( time_sind )%        &
                &                     j_inertial,                            &
                &   b2FaceData = totFace,                                    &
                &   vectorID = VEC_ALIGN_RADIAL_ID )

            !! fchin: Ion-neutral friction current
            hlpFace(:,:,0) = fchin(:,:,0,0)
            hlpFace(:,:,1) = fchin(:,:,1,0)
            call divide_by_contact_areas(nx,ny,hlpFace,totFace)
            call write_face_vector_component( edge_grid,                     &
                &   vectorComponent = edge_profiles%ggd( time_sind )%        &
                &                     j_ion_neutral_friction,                &
                &   b2FaceData = totFace,                                    &
                &   vectorID = VEC_ALIGN_POLOIDAL_ID )
            hlpFace(:,:,0) = fchin(:,:,0,1)
            hlpFace(:,:,1) = fchin(:,:,1,1)
            call divide_by_contact_areas(nx,ny,hlpFace,totFace)
            call write_face_vector_component( edge_grid,                     &
                &   vectorComponent = edge_profiles%ggd( time_sind )%        &
                &                     j_ion_neutral_friction,                &
                &   b2FaceData = totFace,                                    &
                &   vectorID = VEC_ALIGN_RADIAL_ID )

            !! fchvispar: Parallel viscosity current
            hlpFace(:,:,0) = fchvispar(:,:,0,0)
            hlpFace(:,:,1) = fchvispar(:,:,1,0)
            call divide_by_contact_areas(nx,ny,hlpFace,totFace)
            call write_face_vector_component( edge_grid,                     &
                &   vectorComponent = edge_profiles%ggd( time_sind )%        &
                &                     j_parallel_viscosity,                  &
                &   b2FaceData = totFace,                                    &
                &   vectorID = VEC_ALIGN_POLOIDAL_ID )
            hlpFace(:,:,0) = fchvispar(:,:,0,1)
            hlpFace(:,:,1) = fchvispar(:,:,1,1)
            call divide_by_contact_areas(nx,ny,hlpFace,totFace)
            call write_face_vector_component( edge_grid,                     &
                &   vectorComponent = edge_profiles%ggd( time_sind )%        &
                &                     j_parallel_viscosity,                  &
                &   b2FaceData = totFace,                                    &
                &   vectorID = VEC_ALIGN_RADIAL_ID )

            !! fchvisper: Perpendicular viscosity current
            hlpFace(:,:,0) = fchvisper(:,:,0,0)
            hlpFace(:,:,1) = fchvisper(:,:,1,0)
            call divide_by_contact_areas(nx,ny,hlpFace,totFace)
            call write_face_vector_component( edge_grid,                     &
                &   vectorComponent = edge_profiles%ggd( time_sind )%        &
                &                     j_perpendicular_viscosity,             &
                &   b2FaceData = totFace,                                    &
                &   vectorID = VEC_ALIGN_POLOIDAL_ID )
            hlpFace(:,:,0) = fchvisper(:,:,0,1)
            hlpFace(:,:,1) = fchvisper(:,:,1,1)
            call divide_by_contact_areas(nx,ny,hlpFace,totFace)
            call write_face_vector_component( edge_grid,                     &
                &   vectorComponent = edge_profiles%ggd( time_sind )%        &
                &                     j_perpendicular_viscosity,             &
                &   b2FaceData = totFace,                                    &
                &   vectorID = VEC_ALIGN_RADIAL_ID )

            !! fchvisq: Heat viscosity current
            hlpFace(:,:,0) = fchvisq(:,:,0,0)
            hlpFace(:,:,1) = fchvisq(:,:,1,0)
            call divide_by_contact_areas(nx,ny,hlpFace,totFace)
            call write_face_vector_component( edge_grid,                     &
                &   vectorComponent = edge_profiles%ggd( time_sind )%        &
                &                     j_heat_viscosity,                      &
                &   b2FaceData = totFace,                                    &
                &   vectorID = VEC_ALIGN_POLOIDAL_ID )
            hlpFace(:,:,0) = fchvisq(:,:,0,1)
            hlpFace(:,:,1) = fchvisq(:,:,1,1)
            call divide_by_contact_areas(nx,ny,hlpFace,totFace)
            call write_face_vector_component( edge_grid,                     &
                &   vectorComponent = edge_profiles%ggd( time_sind )%        &
                &                     j_heat_viscosity,                      &
                &   b2FaceData = totFace,                                    &
                &   vectorID = VEC_ALIGN_RADIAL_ID )

            !! fchdia: Diamagnetic current
            hlpFace(:,:,0) = fchdia(:,:,0,0)
            hlpFace(:,:,1) = fchdia(:,:,1,0)
            call divide_by_contact_areas(nx,ny,hlpFace,totFace)
            call write_face_vector_component( edge_grid,                     &
                &   vectorComponent = edge_profiles%ggd( time_sind )%        &
                &                     j_diamagnetic,                         &
                &   b2FaceData = totFace,                                    &
                &   vectorID = VEC_ALIGN_POLOIDAL_ID )
            hlpFace(:,:,0) = fchdia(:,:,1,0)
            hlpFace(:,:,1) = fchdia(:,:,1,1)
            call divide_by_contact_areas(nx,ny,hlpFace,totFace)
            call write_face_vector_component( edge_grid,                     &
                &   vectorComponent = edge_profiles%ggd( time_sind )%        &
                &                     j_diamagnetic,                         &
                &   b2FaceData = totFace,                                    &
                &   vectorID = VEC_ALIGN_RADIAL_ID )
#else
            !! ni/ne: Total ion density over electron density
            tmpCv(:,:) = ni(:,:,1)/ne(:,:)
            call write_quantity( edge_grid,                                  &
                &   val = profiles_ggd%n_i_total_over_n_e,                   &
                &   value = tmpCv )

            !! Zeff
            call write_quantity( edge_grid,                                  &
                &   val = profiles_ggd%zeff,                                 &
                &   value = zeff )

            !! pz: Thermal plasma pressure (electrons+ions)
            call write_quantity( edge_grid,                                  &
                &   val = profiles_ggd%pressure_thermal,                     &
                &   value = pz )

            !! fch: Total current
            call divide_by_contact_areas(nx,ny,fch,tmpFace)
            totFace(:,:,0) = tmpFace(:,:,0)
            totFace(:,:,1) = IDS_REAL_INVALID
            call write_face_vector_component( edge_grid,                     &
                &   vectorComponent = profiles_ggd%j_total,                  &
                &   b2FaceData = totFace,                                    &
                &   vectorID = VEC_ALIGN_POLOIDAL_ID )
            totFace(:,:,0) = IDS_REAL_INVALID
            totFace(:,:,1) = tmpFace(:,:,1)
            call write_face_vector_component( edge_grid,                     &
                &   vectorComponent = profiles_ggd%j_total,                  &
                &   b2FaceData = totFace,                                    &
                &   vectorID = VEC_ALIGN_RADIAL_ID )

            !! fch_p: Parallel current
            hlpFace(:,:,0) = fch_p(:,:,0,0)
            hlpFace(:,:,1) = fch_p(:,:,1,1)
            call divide_by_poloidal_areas(nx,ny,hlpFace,tmpFace)
            call write_face_scalar( edge_grid,                               &
                &   val = profiles_ggd%j_parallel,                           &
                &   value = tmpFace )

            !! fchanml: Anomalous current
            call b2tanml (nx, ny, ns, ismain, vol, hx, hy, qs,               & !som 02.11.21
                &         csig_an, po, ne, na, fchanml_a, fchanml)
            hlpFace(:,:,0) = fchanml(:,:,0,0)
            hlpFace(:,:,1) = fchanml(:,:,1,0)
            call divide_by_contact_areas(nx,ny,hlpFace,totFace)
            call write_face_vector_component( edge_grid,                     &
                &   vectorComponent = profiles_ggd%j_anomalous,              &
                &   b2FaceData = totFace,                                    &
                &   vectorID = VEC_ALIGN_POLOIDAL_ID )
            hlpFace(:,:,0) = fchanml(:,:,0,1)
            hlpFace(:,:,1) = fchanml(:,:,1,1)
            call divide_by_contact_areas(nx,ny,hlpFace,totFace)
            call write_face_vector_component( edge_grid,                     &
                &   vectorComponent = profiles_ggd%j_anomalous,              &
                &   b2FaceData = totFace,                                    &
                &   vectorID = VEC_ALIGN_RADIAL_ID )

            !! fchinert: Inertial current
            hlpFace(:,:,0) = fchinert(:,:,0,0)
            hlpFace(:,:,1) = fchinert(:,:,1,0)
            call divide_by_contact_areas(nx,ny,hlpFace,totFace)
            call write_face_vector_component( edge_grid,                     &
                &   vectorComponent = profiles_ggd%j_inertial,               &
                &   b2FaceData = totFace,                                    &
                &   vectorID = VEC_ALIGN_POLOIDAL_ID )
            hlpFace(:,:,0) = fchinert(:,:,0,1)
            hlpFace(:,:,1) = fchinert(:,:,1,1)
            call divide_by_contact_areas(nx,ny,hlpFace,totFace)
            call write_face_vector_component( edge_grid,                     &
                &   vectorComponent = profiles_ggd%j_inertial,               &
                &   b2FaceData = totFace,                                    &
                &   vectorID = VEC_ALIGN_RADIAL_ID )

            !! fchin: Ion-neutral friction current
            hlpFace(:,:,0) = fchin(:,:,0,0)
            hlpFace(:,:,1) = fchin(:,:,1,0)
            call divide_by_contact_areas(nx,ny,hlpFace,totFace)
            call write_face_vector_component( edge_grid,                     &
                &   vectorComponent = profiles_ggd%j_ion_neutral_friction,   &
                &   b2FaceData = totFace,                                    &
                &   vectorID = VEC_ALIGN_POLOIDAL_ID )
            hlpFace(:,:,0) = fchin(:,:,0,1)
            hlpFace(:,:,1) = fchin(:,:,1,1)
            call divide_by_contact_areas(nx,ny,hlpFace,totFace)
            call write_face_vector_component( edge_grid,                     &
                &   vectorComponent = profiles_ggd%j_ion_neutral_friction,   &
                &   b2FaceData = totFace,                                    &
                &   vectorID = VEC_ALIGN_RADIAL_ID )

            !! fchvispar: Parallel viscosity current
            hlpFace(:,:,0) = fchvispar(:,:,0,0)
            hlpFace(:,:,1) = fchvispar(:,:,1,0)
            call divide_by_contact_areas(nx,ny,hlpFace,totFace)
            call write_face_vector_component( edge_grid,                     &
                &   vectorComponent = profiles_ggd%j_parallel_viscosity,     &
                &   b2FaceData = totFace,                                    &
                &   vectorID = VEC_ALIGN_POLOIDAL_ID )
            hlpFace(:,:,0) = fchvispar(:,:,0,1)
            hlpFace(:,:,1) = fchvispar(:,:,1,1)
            call divide_by_contact_areas(nx,ny,hlpFace,totFace)
            call write_face_vector_component( edge_grid,                     &
                &   vectorComponent = profiles_ggd%j_parallel_viscosity,     &
                &   b2FaceData = totFace,                                    &
                &   vectorID = VEC_ALIGN_RADIAL_ID )

            !! fchvisper: Perpendicular viscosity current
            hlpFace(:,:,0) = fchvisper(:,:,0,0)
            hlpFace(:,:,1) = fchvisper(:,:,1,0)
            call divide_by_contact_areas(nx,ny,hlpFace,totFace)
            call write_face_vector_component( edge_grid,                     &
                &   vectorComponent = profiles_ggd%                          &
                &                     j_perpendicular_viscosity,             &
                &   b2FaceData = totFace,                                    &
                &   vectorID = VEC_ALIGN_POLOIDAL_ID )
            hlpFace(:,:,0) = fchvisper(:,:,0,1)
            hlpFace(:,:,1) = fchvisper(:,:,1,1)
            call divide_by_contact_areas(nx,ny,hlpFace,totFace)
            call write_face_vector_component( edge_grid,                     &
                &   vectorComponent = profiles_ggd%                          &
                &                     j_perpendicular_viscosity,             &
                &   b2FaceData = totFace,                                    &
                &   vectorID = VEC_ALIGN_RADIAL_ID )

            !! fchvisq: Heat viscosity current
            hlpFace(:,:,0) = fchvisq(:,:,0,0)
            hlpFace(:,:,1) = fchvisq(:,:,1,0)
            call divide_by_contact_areas(nx,ny,hlpFace,totFace)
            call write_face_vector_component( edge_grid,                     &
                &   vectorComponent = profiles_ggd%j_heat_viscosity,         &
                &   b2FaceData = totFace,                                    &
                &   vectorID = VEC_ALIGN_POLOIDAL_ID )
            hlpFace(:,:,0) = fchvisq(:,:,0,1)
            hlpFace(:,:,1) = fchvisq(:,:,1,1)
            call divide_by_contact_areas(nx,ny,hlpFace,totFace)
            call write_face_vector_component( edge_grid,                     &
                &   vectorComponent = profiles_ggd%j_heat_viscosity,         &
                &   b2FaceData = totFace,                                    &
                &   vectorID = VEC_ALIGN_RADIAL_ID )

            !! fchdia: Diamagnetic current
            hlpFace(:,:,0) = fchdia(:,:,0,0)
            hlpFace(:,:,1) = fchdia(:,:,1,0)
            call divide_by_contact_areas(nx,ny,hlpFace,totFace)
            call write_face_vector_component( edge_grid,                     &
                &   vectorComponent = profiles_ggd%j_diamagnetic,            &
                &   b2FaceData = totFace,                                    &
                &   vectorID = VEC_ALIGN_POLOIDAL_ID )
            hlpFace(:,:,0) = fchdia(:,:,1,0)
            hlpFace(:,:,1) = fchdia(:,:,1,1)
            call divide_by_contact_areas(nx,ny,hlpFace,totFace)
            call write_face_vector_component( edge_grid,                     &
                &   vectorComponent = profiles_ggd%j_diamagnetic,            &
                &   b2FaceData = totFace,                                    &
                &   vectorID = VEC_ALIGN_RADIAL_ID )
#endif

            if (use_eirene.ne.0) then
#ifdef B25_EIRENE
                do is = 1, nspecies
                   tmpCv(:,:) = 0.0_IDS_real
                   totCv(:,:) = 0.0_IDS_real
                   tmpFace(:,:,:) = 0.0_IDS_real
                   do iss = 1, natmi
                      if (latmscl(iss).eq.is) then
                        tmpCv(-1:nx,-1:ny) = tmpCv(-1:nx,-1:ny) +            &
                           &  dab2(0:nx+1,0:ny+1,iss,1)*tab2(0:nx+1,0:ny+1,iss,1)
                        totCv(-1:nx,-1:ny) = totCv(-1:nx,-1:ny) +            &
                           &  dab2(0:nx+1,0:ny+1,iss,1)
                        tmpFace(-1:nx,-1:ny,0) = tmpFace(-1:nx,-1:ny,0) +    &
                           &  pfluxa(0:nx+1,0:ny+1,iss,1)
                        tmpFace(-1:nx,-1:ny,1) = tmpFace(-1:nx,-1:ny,1) +    &
                           &  rfluxa(0:nx+1,0:ny+1,iss,1)
                      end if
                   end do
#if ( IMAS_MAJOR_VERSION < 4 || ( IMAS_MAJOR_VERSION == 4 && IMAS_MINOR_VERSION < 1 ) )
                !! Neutral pressure
                   call write_quantity( edge_grid,                           &
                       &   val = edge_profiles%ggd( time_sind )%             &
                       &         neutral( is )%pressure,                     &
                       &   value = tmpCv )
                !! Neutral density
                   call write_quantity( edge_grid,                           &
                       &   val = edge_profiles%ggd( time_sind )%             &
                       &         neutral( is )%density,                      &
                       &   value = totCv )
                !! Neutral particle flux
                   call write_face_scalar( transport_grid,                   &
                       &   val = edge_transport%model(1)%ggd( time_sind )%   &
                       &         neutral( is )%particles%flux,               &
                       &   value = tmpFace )
#else
                !! Neutral pressure
                   call write_quantity( edge_grid,                           &
                       &   val = profiles_ggd%neutral( is )%pressure,        &
                       &   value = tmpCv )
                !! Neutral density
                   call write_quantity( edge_grid,                           &
                       &   val = profiles_ggd%neutral( is )%density,         &
                       &   value = totCv )
                !! Neutral particle flux
                   call write_face_scalar( transport_grid,                   &
                       &   val = transport_ggd(1)%                           &
                       &         neutral( is )%particles%flux,               &
                       &   value = tmpFace )
#endif
                !! Neutral velocity (poloidal projection)
                   if (allocated(un0)) then
                     tmpCv(:,:) = 0.0_IDS_real
                     do iss = 1, natmi
                       if (latmscl(iss).eq.is) then
                         tmpCv(-1:nx,-1:ny) = tmpCv(-1:nx,-1:ny) +           &
                           &  un0(-1:nx,-1:ny,0,iss)*dab2(0:nx+1,0:ny+1,iss,1)
                       end if
                     end do
                     do iy = -1, ny
                       do ix = -1, nx
                         if (totCv(ix,iy).gt.0.0_IDS_real)                   &
                           &  tmpCv(ix,iy) = tmpCv(ix,iy) / totCv(ix,iy)
                       end do
                     end do
#if ( IMAS_MAJOR_VERSION < 4 || ( IMAS_MAJOR_VERSION == 4 && IMAS_MINOR_VERSION < 1 ) )
                     call write_cell_vector_component( edge_grid,            &
                       &   vectorComponent = edge_profiles%ggd( time_sind )% &
                       &                     neutral( is )%velocity,         &
                       &   b2CellData = tmpCv(:,:),                          &
                       &   vectorID = VEC_ALIGN_POLOIDAL_ID )
#else
                     call write_cell_vector_component( edge_grid,            &
                       &   vectorComponent = profiles_ggd%                   &
                       &                     neutral( is )%velocity,         &
                       &   b2CellData = tmpCv(:,:),                          &
                       &   vectorID = VEC_ALIGN_POLOIDAL_ID )
#endif
                !! Neutral velocity (radial projection)
                     tmpCv(:,:) = 0.0_IDS_real
                     do iss = 1, natmi
                       if (latmscl(iss).eq.is) then
                         tmpCv(-1:nx,-1:ny) = tmpCv(-1:nx,-1:ny) +           &
                           &  un0(-1:nx,-1:ny,1,iss)*dab2(0:nx+1,0:ny+1,iss,1)
                       end if
                     end do
                     do iy = -1, ny
                       do ix = -1, nx
                         if (totCv(ix,iy).gt.0.0_IDS_real)                   &
                           &  tmpCv(ix,iy) = tmpCv(ix,iy) / totCv(ix,iy)
                       end do
                     end do
#if ( IMAS_MAJOR_VERSION < 4 || ( IMAS_MAJOR_VERSION == 4 && IMAS_MINOR_VERSION < 1 ) )
                     call write_cell_vector_component( edge_grid,            &
                       &   vectorComponent = edge_profiles%ggd( time_sind )% &
                       &                     neutral( is )%velocity,         &
                       &   b2CellData = tmpCv(:,:),                          &
                       &   vectorID = VEC_ALIGN_RADIAL_ID )
#else
                     call write_cell_vector_component( edge_grid,            &
                       &   vectorComponent = profiles_ggd%                   &
                       &                     neutral( is )%velocity,         &
                       &   b2CellData = tmpCv(:,:),                          &
                       &   vectorID = VEC_ALIGN_RADIAL_ID )
#endif
                !! Neutral velocity (toroidal projection)
                     tmpCv(:,:) = 0.0_IDS_real
                     do iss = 1, natmi
                       if (latmscl(iss).eq.is) then
                         tmpCv(-1:nx,-1:ny) = tmpCv(-1:nx,-1:ny) +           &
                           &  un0(-1:nx,-1:ny,2,iss)*dab2(0:nx+1,0:ny+1,iss,1)
                       end if
                     end do
                     do iy = -1, ny
                       do ix = -1, nx
                         if (totCv(ix,iy).gt.0.0_IDS_real)                   &
                           &  tmpCv(ix,iy) = tmpCv(ix,iy) / totCv(ix,iy)
                       end do
                     end do
#if ( GGD_MAJOR_VERSION < 2 && GGD_MINOR_VERSION < 13 )
                     call write_cell_vector_component( edge_grid,            &
                       &   vectorComponent = edge_profiles%ggd( time_sind )% &
                       &                     neutral( is )%velocity,         &
                       &   b2CellData = tmpCv(:,:),                          &
                       &   vectorID = VEC_ALIGN_TOROIDAL_ID )
#else
#if ( IMAS_MAJOR_VERSION < 4 || ( IMAS_MAJOR_VERSION == 4 && IMAS_MINOR_VERSION < 1 ) )
                     call write_cell_vector_component( edge_grid,            &
                       &   vectorComponent = edge_profiles%ggd( time_sind )% &
                       &                     neutral( is )%velocity,         &
                       &   b2CellData = tmpCv(:,:),                          &
                       &   vectorID = VEC_ALIGN_PHI_ID )
#else
                     call write_cell_vector_component( edge_grid,            &
                       &   vectorComponent = profiles_ggd%                   &
                       &                     neutral( is )%velocity,         &
                       &   b2CellData = tmpCv(:,:),                          &
                       &   vectorID = VEC_ALIGN_PHI_ID )
#endif
#endif
                   end if
#if ( IMAS_MAJOR_VERSION < 4 || ( IMAS_MAJOR_VERSION == 4 && IMAS_MINOR_VERSION < 1 ) )
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
                   call write_face_scalar( transport_grid,                   &
                       &   val = edge_transport%model(1)%ggd( time_sind )%   &
                       &         neutral( is )%energy%flux,                  &
                       &   value = tmpFace )
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
                     call write_cell_scalar( sources_grid,                   &
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
                     call write_cell_scalar( sources_grid,                   &
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
                     call write_cell_scalar( sources_grid,                   &
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
                   call write_cell_scalar( sources_grid,                     &
                       &   scalar = edge_sources%source(12)%ggd( time_sind )%&
                       &            neutral( is )%energy,                    &
                       &   b2CellData = tmpCv )
#else
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
                   call write_face_scalar( transport_grid,                   &
                       &   val = transport_ggd(1)%neutral( is )%energy%flux, &
                       &   value = tmpFace )
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
                     call write_cell_scalar( sources_grid,                   &
                       &   scalar = sources_ggd(1)%neutral( is )%particles,  &
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
                     call write_cell_scalar( sources_grid,                   &
                       &   scalar = sources_ggd(5)%neutral( is )%particles,  &
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
                     call write_cell_scalar( sources_grid,                   &
                       &   scalar = sources_ggd(6)%neutral( is )%particles,  &
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
                   call write_cell_scalar( sources_grid,                     &
                       &   scalar = sources_ggd(12)%neutral( is )%particles, &
                       &   b2CellData = tmpCv )
#endif
                end do
                do is = 1, natmi
                   js = latmscl(is)
                   ks = isstat(is)
#if ( IMAS_MAJOR_VERSION < 4 || ( IMAS_MAJOR_VERSION == 4 && IMAS_MINOR_VERSION < 1 ) )
                   tmpCv(-1:nx,-1:ny) = tab2(0:nx+1,0:ny+1,is,1)/qe
                   call write_quantity( edge_grid,                           &
                       &   val = edge_profiles%ggd( time_sind )%             &
                       &         neutral( js )%state( ks )%temperature,      &
                       &   value = tmpCv )
                   tmpCv(-1:nx,-1:ny) = dab2(0:nx+1,0:ny+1,is,1)
                   call write_quantity( edge_grid,                           &
                       &   val = edge_profiles%ggd( time_sind )%             &
                       &         neutral( js )%state( ks )%density,          &
                       &   value = tmpCv )
                   tmpFace(-1:nx,-1:ny,0) = pfluxa(0:nx+1,0:ny+1,is,1)
                   tmpFace(-1:nx,-1:ny,1) = rfluxa(0:nx+1,0:ny+1,is,1)
                   call write_face_scalar( transport_grid,                   &
                       &   val = edge_transport%model(1)%ggd( time_sind )%   &
                       &         neutral( js )%state( ks )%particles%flux,   &
                       &   value = tmpFace )
                   tmpFace(-1:nx,-1:ny,0) = pefluxa(0:nx+1,0:ny+1,is,1)
                   tmpFace(-1:nx,-1:ny,1) = refluxa(0:nx+1,0:ny+1,is,1)
                   call write_face_scalar( transport_grid,                   &
                       &   val = edge_transport%model(1)%ggd( time_sind )%   &
                       &         neutral( js )%state( ks )%energy%flux,      &
                       &   value = tmpFace )
                   tmpCv(-1:nx,-1:ny) = dab2(0:nx+1,0:ny+1,is,1)*            &
                       &                tab2(0:nx+1,0:ny+1,is,1)
                   call write_quantity( edge_grid,                           &
                       &   val = edge_profiles%ggd( time_sind )%             &
                       &         neutral( js )%state( ks )%pressure,         &
                       &   value = tmpCv )
                   if (allocated(un0)) then
                     call write_cell_vector_component( edge_grid,            &
                       &   vectorComponent = edge_profiles%ggd( time_sind )% &
                       &         neutral( js )%state( ks )%velocity,         &
                       &   b2CellData = un0(:,:,0,is),                       &
                       &   vectorID = VEC_ALIGN_POLOIDAL_ID )
                     call write_cell_vector_component( edge_grid,            &
                       &   vectorComponent = edge_profiles%ggd( time_sind )% &
                       &         neutral( js )%state( ks )%velocity,         &
                       &   b2CellData = un0(:,:,1,is),                       &
                       &   vectorID = VEC_ALIGN_RADIAL_ID )
#if ( GGD_MAJOR_VERSION < 2 && GGD_MINOR_VERSION < 13 )
                     call write_cell_vector_component( edge_grid,            &
                       &   vectorComponent = edge_profiles%ggd( time_sind )% &
                       &         neutral( js )%state( ks )%velocity,         &
                       &   b2CellData = un0(:,:,2,is),                       &
                       &   vectorID = VEC_ALIGN_TOROIDAL_ID )
#else
                     call write_cell_vector_component( edge_grid,            &
                       &   vectorComponent = edge_profiles%ggd( time_sind )% &
                       &         neutral( js )%state( ks )%velocity,         &
                       &   b2CellData = un0(:,:,2,is),                       &
                       &   vectorID = VEC_ALIGN_PHI_ID )
#endif
                   end if
                   if (drift_style.eq.0) then
                     tmpCv = 0.0_IDS_real
                     call write_cell_vector_component( edge_grid,            &
                       &   vectorComponent = edge_profiles%ggd( time_sind )% &
                       &                     neutral( js )%state( ks )%      &
                       &                     velocity_diamagnetic,           &
                       &   b2CellData = tmpCv,                               &
                       &   vectorID = VEC_ALIGN_POLOIDAL_ID )
                     call write_cell_vector_component( edge_grid,            &
                       &   vectorComponent = edge_profiles%ggd( time_sind )% &
                       &                     neutral( js )%state( ks )%      &
                       &                     velocity_diamagnetic,           &
                       &   b2CellData = tmpCv,                               &
                       &   vectorID = VEC_ALIGN_RADIAL_ID )
                     call write_cell_vector_component( edge_grid,            &
                       &   vectorComponent = edge_profiles%ggd( time_sind )% &
                       &                     neutral( js )%state( ks )%      &
                       &                     velocity_exb,                   &
                       &   b2CellData = tmpCv,                               &
                       &   vectorID = VEC_ALIGN_POLOIDAL_ID )
                     call write_cell_vector_component( edge_grid,            &
                       &   vectorComponent = edge_profiles%ggd( time_sind )% &
                       &                     neutral( js )%state( ks )%      &
                       &                     velocity_exb,                   &
                       &   b2CellData = tmpCv,                               &
                       &   vectorID = VEC_ALIGN_RADIAL_ID )
                   else
                     tmpFace = 0.0_IDS_real
                     call write_face_vector_component( edge_grid,            &
                       &   vectorComponent = edge_profiles%ggd( time_sind )% &
                       &                     neutral( js )%state( ks )%      &
                       &                     velocity_diamagnetic,           &
                       &   b2FaceData = tmpFace,                             &
                       &   vectorID = VEC_ALIGN_POLOIDAL_ID )
                     call write_face_vector_component( edge_grid,            &
                       &   vectorComponent = edge_profiles%ggd( time_sind )% &
                       &                     neutral( js )%state( ks )%      &
                       &                     velocity_diamagnetic,           &
                       &   b2FaceData = tmpFace,                             &
                       &   vectorID = VEC_ALIGN_RADIAL_ID )
                     call write_face_vector_component( edge_grid,            &
                       &   vectorComponent = edge_profiles%ggd( time_sind )% &
                       &                     neutral( js )%state( ks )%      &
                       &                     velocity_exb,                   &
                       &   b2FaceData = tmpFace,                             &
                       &   vectorID = VEC_ALIGN_POLOIDAL_ID )
                     call write_face_vector_component( edge_grid,            &
                       &   vectorComponent = edge_profiles%ggd( time_sind )% &
                       &                     neutral( js )%state( ks )%      &
                       &                     velocity_exb,                   &
                       &   b2FaceData = tmpFace,                             &
                       &   vectorID = VEC_ALIGN_RADIAL_ID )
                   end if
#else
                   tmpCv(-1:nx,-1:ny) = tab2(0:nx+1,0:ny+1,is,1)/qe
                   call write_quantity( edge_grid,                           &
                       &   val = profiles_ggd%                               &
                       &         neutral( js )%state( ks )%temperature,      &
                       &   value = tmpCv )
                   tmpCv(-1:nx,-1:ny) = dab2(0:nx+1,0:ny+1,is,1)
                   call write_quantity( edge_grid,                           &
                       &   val = profiles_ggd%                               &
                       &         neutral( js )%state( ks )%density,          &
                       &   value = tmpCv )
                   tmpFace(-1:nx,-1:ny,0) = pfluxa(0:nx+1,0:ny+1,is,1)
                   tmpFace(-1:nx,-1:ny,1) = rfluxa(0:nx+1,0:ny+1,is,1)
                   call write_face_scalar( transport_grid,                   &
                       &   val = transport_ggd(1)%                           &
                       &         neutral( js )%state( ks )%particles%flux,   &
                       &   value = tmpFace )
                   tmpFace(-1:nx,-1:ny,0) = pefluxa(0:nx+1,0:ny+1,is,1)
                   tmpFace(-1:nx,-1:ny,1) = refluxa(0:nx+1,0:ny+1,is,1)
                   call write_face_scalar( transport_grid,                   &
                       &   val = transport_ggd(1)%                           &
                       &         neutral( js )%state( ks )%energy%flux,      &
                       &   value = tmpFace )
                   tmpCv(-1:nx,-1:ny) = dab2(0:nx+1,0:ny+1,is,1)*            &
                       &                tab2(0:nx+1,0:ny+1,is,1)
                   call write_quantity( edge_grid,                           &
                       &   val = profiles_ggd%                               &
                       &         neutral( js )%state( ks )%pressure,         &
                       &   value = tmpCv )
                   if (allocated(un0)) then
                     call write_cell_vector_component( edge_grid,            &
                       &   vectorComponent = profiles_ggd%                   &
                       &         neutral( js )%state( ks )%velocity,         &
                       &   b2CellData = un0(:,:,0,is),                       &
                       &   vectorID = VEC_ALIGN_POLOIDAL_ID )
                     call write_cell_vector_component( edge_grid,            &
                       &   vectorComponent = profiles_ggd%                   &
                       &         neutral( js )%state( ks )%velocity,         &
                       &   b2CellData = un0(:,:,1,is),                       &
                       &   vectorID = VEC_ALIGN_RADIAL_ID )
                     call write_cell_vector_component( edge_grid,            &
                       &   vectorComponent = profiles_ggd%                   &
                       &         neutral( js )%state( ks )%velocity,         &
                       &   b2CellData = un0(:,:,2,is),                       &
                       &   vectorID = VEC_ALIGN_PHI_ID )
                   end if
                   if (drift_style.eq.0) then
                     tmpCv = 0.0_IDS_real
                     call write_cell_vector_component( edge_grid,            &
                       &   vectorComponent = profiles_ggd%                   &
                       &                     neutral( js )%state( ks )%      &
                       &                     velocity_diamagnetic,           &
                       &   b2CellData = tmpCv,                               &
                       &   vectorID = VEC_ALIGN_POLOIDAL_ID )
                     call write_cell_vector_component( edge_grid,            &
                       &   vectorComponent = profiles_ggd%                   &
                       &                     neutral( js )%state( ks )%      &
                       &                     velocity_diamagnetic,           &
                       &   b2CellData = tmpCv,                               &
                       &   vectorID = VEC_ALIGN_RADIAL_ID )
                     call write_cell_vector_component( edge_grid,            &
                       &   vectorComponent = profiles_ggd%                   &
                       &                     neutral( js )%state( ks )%      &
                       &                     velocity_exb,                   &
                       &   b2CellData = tmpCv,                               &
                       &   vectorID = VEC_ALIGN_POLOIDAL_ID )
                     call write_cell_vector_component( edge_grid,            &
                       &   vectorComponent = profiles_ggd%                   &
                       &                     neutral( js )%state( ks )%      &
                       &                     velocity_exb,                   &
                       &   b2CellData = tmpCv,                               &
                       &   vectorID = VEC_ALIGN_RADIAL_ID )
                   else
                     tmpFace = 0.0_IDS_real
                     call write_face_vector_component( edge_grid,            &
                       &   vectorComponent = profiles_ggd%                   &
                       &                     neutral( js )%state( ks )%      &
                       &                     velocity_diamagnetic,           &
                       &   b2FaceData = tmpFace,                             &
                       &   vectorID = VEC_ALIGN_POLOIDAL_ID )
                     call write_face_vector_component( edge_grid,            &
                       &   vectorComponent = profiles_ggd%                   &
                       &                     neutral( js )%state( ks )%      &
                       &                     velocity_diamagnetic,           &
                       &   b2FaceData = tmpFace,                             &
                       &   vectorID = VEC_ALIGN_RADIAL_ID )
                     call write_face_vector_component( edge_grid,            &
                       &   vectorComponent = profiles_ggd%                   &
                       &                     neutral( js )%state( ks )%      &
                       &                     velocity_exb,                   &
                       &   b2FaceData = tmpFace,                             &
                       &   vectorID = VEC_ALIGN_POLOIDAL_ID )
                     call write_face_vector_component( edge_grid,            &
                       &   vectorComponent = profiles_ggd%                   &
                       &                     neutral( js )%state( ks )%      &
                       &                     velocity_exb,                   &
                       &   b2FaceData = tmpFace,                             &
                       &   vectorID = VEC_ALIGN_RADIAL_ID )
                   end if
#endif
#if ( IMAS_MAJOR_VERSION < 4 || ( IMAS_MAJOR_VERSION == 4 && IMAS_MINOR_VERSION < 1 ) )
                   if (balance_netcdf.ne.0) then
                     tmpCv = 0.0_IDS_real
                     do istrai = 1, size(eirene_mc_paat_sna_bal,4)
                       tmpCv(:,:) = tmpCv(:,:)                               &
                          &       + eirene_mc_paat_sna_bal(:,:,is,istrai)    &
                          &       + eirene_mc_pmat_sna_bal(:,:,is,istrai)    &
                          &       + eirene_mc_piat_sna_bal(:,:,is,istrai)
                     end do
                     tmpCv(:,:) = tmpCv(:,:) / vol(:,:)
                     call write_cell_scalar( sources_grid,                   &
                       &   scalar = edge_sources%source(1)%ggd( time_sind )% &
                       &            neutral( js )%state( ks )%particles,     &
                       &   b2CellData = tmpCv )
                     tmpCv = 0.0_IDS_real
                     do istrai = 1, size(eirene_mc_paat_sna_bal,4)
                       tmpCv(:,:) = tmpCv(:,:)                               &
                          &       + eirene_mc_paat_sna_bal(:,:,is,istrai)
                     end do
                     tmpCv(:,:) = tmpCv(:,:) / vol(:,:)
                     call write_cell_scalar( sources_grid,                   &
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
                     call write_cell_scalar( sources_grid,                   &
                       &   scalar = edge_sources%source(6)%ggd( time_sind )% &
                       &            neutral( js )%state( ks )%particles,     &
                       &   b2CellData = tmpCv )
                   end if
                   tmpCv(-1:nx,-1:ny) = eneutrad(0:nx+1,0:ny+1,is,0)/        &
                       &                vol(-1:nx,-1:ny)
                   call write_cell_scalar( sources_grid,                     &
                       &   scalar = edge_sources%source(12)%                 &
                       &            ggd( time_sind )%neutral( js )%          &
                       &            state( ks )%energy,                      &
                       &   b2CellData = tmpCv )
#else
                   if (balance_netcdf.ne.0) then
                     tmpCv = 0.0_IDS_real
                     do istrai = 1, size(eirene_mc_paat_sna_bal,4)
                       tmpCv(:,:) = tmpCv(:,:)                               &
                          &       + eirene_mc_paat_sna_bal(:,:,is,istrai)    &
                          &       + eirene_mc_pmat_sna_bal(:,:,is,istrai)    &
                          &       + eirene_mc_piat_sna_bal(:,:,is,istrai)
                     end do
                     tmpCv(:,:) = tmpCv(:,:) / vol(:,:)
                     call write_cell_scalar( sources_grid,                   &
                       &   scalar = sources_ggd(1)%                          &
                       &            neutral( js )%state( ks )%particles,     &
                       &   b2CellData = tmpCv )
                     tmpCv = 0.0_IDS_real
                     do istrai = 1, size(eirene_mc_paat_sna_bal,4)
                       tmpCv(:,:) = tmpCv(:,:)                               &
                          &       + eirene_mc_paat_sna_bal(:,:,is,istrai)
                     end do
                     tmpCv(:,:) = tmpCv(:,:) / vol(:,:)
                     call write_cell_scalar( sources_grid,                   &
                       &   scalar = sources_ggd(5)%                          &
                       &            neutral( js )%state( ks )%particles,     &
                       &   b2CellData = tmpCv )
                     tmpCv = 0.0_IDS_real
                     do istrai = 1, size(eirene_mc_pmat_sna_bal,4)
                       tmpCv(:,:) = tmpCv(:,:)                               &
                          &       + eirene_mc_pmat_sna_bal(:,:,is,istrai)    &
                          &       + eirene_mc_piat_sna_bal(:,:,is,istrai)
                     end do
                     tmpCv(:,:) = tmpCv(:,:) / vol(:,:)
                     call write_cell_scalar( sources_grid,                   &
                       &   scalar = sources_ggd(6)%                          &
                       &            neutral( js )%state( ks )%particles,     &
                       &   b2CellData = tmpCv )
                   end if
                   tmpCv(-1:nx,-1:ny) = eneutrad(0:nx+1,0:ny+1,is,0)/        &
                       &                vol(-1:nx,-1:ny)
                   call write_cell_scalar( sources_grid,                     &
                       &   scalar = sources_ggd(12)%                         &
                       &            neutral( js )%state( ks )%energy,        &
                       &   b2CellData = tmpCv )
#endif
                end do

                !! Molecular quantities
                do js = nspecies+1, nneut
                   tmpCv(:,:) = 0.0_IDS_real
                   totCv(:,:) = 0.0_IDS_real
                   tmpFace(:,:,:) = 0.0_IDS_real
                   do is = 1, nmoli
                      if (imneut(is).eq.js) then
                        tmpCv(-1:nx,-1:ny) = tmpCv(-1:nx,-1:ny) +            &
                           &  dmb2(0:nx+1,0:ny+1,is,1)*tmb2(0:nx+1,0:ny+1,is,1)
                        totCv(-1:nx,-1:ny) = totCv(-1:nx,-1:ny) +            &
                           &  dmb2(0:nx+1,0:ny+1,is,1)
                        tmpFace(-1:nx,-1:ny,0) = tmpFace(-1:nx,-1:ny,0) +    &
                           &  pfluxm(0:nx+1,0:ny+1,is,1)
                        tmpFace(-1:nx,-1:ny,1) = tmpFace(-1:nx,-1:ny,1) +    &
                           &  rfluxm(0:nx+1,0:ny+1,is,1)
                      end if
                   end do
#if ( IMAS_MAJOR_VERSION < 4 || ( IMAS_MAJOR_VERSION == 4 && IMAS_MINOR_VERSION < 1 ) )
                 !! Molecular pressure
                   call write_quantity( edge_grid,                           &
                       &   val = edge_profiles%ggd( time_sind )%             &
                       &         neutral( js )%pressure,                     &
                       &   value = tmpCv )
                 !! Molecular density
                   call write_quantity( edge_grid,                           &
                       &   val = edge_profiles%ggd( time_sind )%             &
                       &         neutral( js )%density,                      &
                       &   value = totCv )
                 !! Molecular particular flux
                   call write_face_scalar( transport_grid,                   &
                       &   val = edge_transport%model(1)%ggd( time_sind )%   &
                       &         neutral( js )%particles%flux,               &
                       &   value = tmpFace )
#else
                 !! Molecular pressure
                   call write_quantity( edge_grid,                           &
                       &   val = profiles_ggd%neutral( js )%pressure,        &
                       &   value = tmpCv )
                 !! Molecular density
                   call write_quantity( edge_grid,                           &
                       &   val = profiles_ggd%neutral( js )%density,         &
                       &   value = totCv )
                 !! Molecular particular flux
                   call write_face_scalar( transport_grid,                   &
                       &   val = transport_ggd(1)%                           &
                       &         neutral( js )%particles%flux,               &
                       &   value = tmpFace )
#endif
                !! ua: Molecular velocity (poloidal projection)
                   if (allocated(um0)) then
                     tmpCv(:,:) = 0.0_IDS_real
                     do is = 1, nmoli
                       if (imneut(is).eq.js) then
                         tmpCv(-1:nx,-1:ny) = tmpCv(-1:nx,-1:ny) +           &
                           &  um0(-1:nx,-1:ny,0,is)*dmb2(0:nx+1,0:ny+1,is,1)
                       end if
                     end do
                     do iy = -1, ny
                       do ix = -1, nx
                         if (totCv(ix,iy).gt.0.0_IDS_real)                   &
                           &  tmpCv(ix,iy) = tmpCv(ix,iy) / totCv(ix,iy)
                       end do
                     end do
#if ( IMAS_MAJOR_VERSION < 4 || ( IMAS_MAJOR_VERSION == 4 && IMAS_MINOR_VERSION < 1 ) )
                     call write_cell_vector_component( edge_grid,            &
                       &   vectorComponent = edge_profiles%ggd( time_sind )% &
                       &                     neutral( js )%velocity,         &
                       &   b2CellData = tmpCv(:,:),                          &
                       &   vectorID = VEC_ALIGN_POLOIDAL_ID )
#else
                     call write_cell_vector_component( edge_grid,            &
                       &   vectorComponent = profiles_ggd%                   &
                       &                     neutral( js )%velocity,         &
                       &   b2CellData = tmpCv(:,:),                          &
                       &   vectorID = VEC_ALIGN_POLOIDAL_ID )
#endif
                !! Molecular velocity (radial projection)
                     tmpCv(:,:) = 0.0_IDS_real
                     do is = 1, nmoli
                       if (imneut(is).eq.js) then
                         tmpCv(-1:nx,-1:ny) = tmpCv(-1:nx,-1:ny) +           &
                           &  um0(-1:nx,-1:ny,1,is)*dmb2(0:nx+1,0:ny+1,is,1)
                       end if
                     end do
                     do iy = -1, ny
                       do ix = -1, nx
                         if (totCv(ix,iy).gt.0.0_IDS_real)                   &
                           &  tmpCv(ix,iy) = tmpCv(ix,iy) / totCv(ix,iy)
                       end do
                     end do
#if ( IMAS_MAJOR_VERSION < 4 || ( IMAS_MAJOR_VERSION == 4 && IMAS_MINOR_VERSION < 1 ) )
                     call write_cell_vector_component( edge_grid,            &
                       &   vectorComponent = edge_profiles%ggd( time_sind )% &
                       &                     neutral( js )%velocity,         &
                       &   b2CellData = tmpCv(:,:),                          &
                       &   vectorID = VEC_ALIGN_RADIAL_ID )
#else
                     call write_cell_vector_component( edge_grid,            &
                       &   vectorComponent = profiles_ggd%                   &
                       &                     neutral( js )%velocity,         &
                       &   b2CellData = tmpCv(:,:),                          &
                       &   vectorID = VEC_ALIGN_RADIAL_ID )
#endif
                !! Molecular velocity (toroidal projection)
                     tmpCv(:,:) = 0.0_IDS_real
                     do is = 1, nmoli
                       if (imneut(is).eq.js) then
                         tmpCv(-1:nx,-1:ny) = tmpCv(-1:nx,-1:ny) +           &
                           &  um0(-1:nx,-1:ny,2,is)*dmb2(0:nx+1,0:ny+1,is,1)
                       end if
                     end do
                     do iy = -1, ny
                       do ix = -1, nx
                         if (totCv(ix,iy).gt.0.0_IDS_real)                   &
                           &  tmpCv(ix,iy) = tmpCv(ix,iy) / totCv(ix,iy)
                       end do
                     end do
#if ( IMAS_MAJOR_VERSION < 4 || ( IMAS_MAJOR_VERSION == 4 && IMAS_MINOR_VERSION < 1 ) )
#if ( GGD_MAJOR_VERSION < 2 && GGD_MINOR_VERSION < 13 )
                     call write_cell_vector_component( edge_grid,            &
                       &   vectorComponent = edge_profiles%ggd( time_sind )% &
                       &                     neutral( js )%velocity,         &
                       &   b2CellData = tmpCv(:,:),                          &
                       &   vectorID = VEC_ALIGN_TOROIDAL_ID )
#else
                     call write_cell_vector_component( edge_grid,            &
                       &   vectorComponent = edge_profiles%ggd( time_sind )% &
                       &                     neutral( js )%velocity,         &
                       &   b2CellData = tmpCv(:,:),                          &
                       &   vectorID = VEC_ALIGN_PHI_ID )
#endif
#else
                     call write_cell_vector_component( edge_grid,            &
                       &   vectorComponent = profiles_ggd%                   &
                       &                     neutral( js )%velocity,         &
                       &   b2CellData = tmpCv(:,:),                          &
                       &   vectorID = VEC_ALIGN_PHI_ID )
#endif
                   end if
#if ( IMAS_MAJOR_VERSION < 4 || ( IMAS_MAJOR_VERSION == 4 && IMAS_MINOR_VERSION < 1 ) )
                 !! Molecular particle energy flux
                   tmpFace(:,:,:) = 0.0_IDS_real
                   do is = 1, nmoli
                      if (imneut(is).eq.js) then
                        tmpFace(-1:nx,-1:ny,0) = tmpFace(-1:nx,-1:ny,0) + &
                           &  pefluxm(0:nx+1,0:ny+1,is,1)
                        tmpFace(-1:nx,-1:ny,1) = tmpFace(-1:nx,-1:ny,1) + &
                           &  refluxm(0:nx+1,0:ny+1,is,1)
                      end if
                   end do
                   call write_face_scalar( transport_grid,                   &
                       &   val = edge_transport%model(1)%ggd( time_sind )%   &
                       &         neutral( js )%energy%flux,                  &
                       &   value = tmpFace )
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
                     call write_cell_scalar( sources_grid,                   &
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
                     call write_cell_scalar( sources_grid,                   &
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
                     call write_cell_scalar( sources_grid,                   &
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
                   call write_cell_scalar( sources_grid,                     &
                       &   scalar = edge_sources%source(12)%                 &
                       &            ggd( time_sind )%neutral( js )%energy,   &
                       &   b2CellData = tmpCv )
#else
                 !! Molecular particle energy flux
                   tmpFace(:,:,:) = 0.0_IDS_real
                   do is = 1, nmoli
                      if (imneut(is).eq.js) then
                        tmpFace(-1:nx,-1:ny,0) = tmpFace(-1:nx,-1:ny,0) + &
                           &  pefluxm(0:nx+1,0:ny+1,is,1)
                        tmpFace(-1:nx,-1:ny,1) = tmpFace(-1:nx,-1:ny,1) + &
                           &  refluxm(0:nx+1,0:ny+1,is,1)
                      end if
                   end do
                   call write_face_scalar( transport_grid,                   &
                       &   val = transport_ggd(1)%neutral( js )%energy%flux, &
                       &   value = tmpFace )
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
                     call write_cell_scalar( sources_grid,                   &
                       &   scalar = sources_ggd(1)%neutral( js )%particles,  &
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
                     call write_cell_scalar( sources_grid,                   &
                       &   scalar = sources_ggd(5)%neutral( js )%particles,  &
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
                     call write_cell_scalar( sources_grid,                   &
                       &   scalar = sources_ggd(6)%neutral( js )%particles,  &
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
                   call write_cell_scalar( sources_grid,                     &
                       &   scalar = sources_ggd(12)%neutral( js )%energy,    &
                       &   b2CellData = tmpCv )
#endif
                end do

                js = nspecies
                do is = 1, nmoli
                   ks = isstat(natmi+is)
                   if (ks.eq.1) js = js + 1
#if ( IMAS_MAJOR_VERSION < 4 || ( IMAS_MAJOR_VERSION == 4 && IMAS_MINOR_VERSION < 1 ) )
                   tmpCv(-1:nx,-1:ny) = tmb2(0:nx+1,0:ny+1,is,1)/qe
                   call write_quantity( edge_grid,                           &
                       &   val = edge_profiles%ggd( time_sind )%             &
                       &         neutral( js )%state( ks )%temperature,      &
                       &   value = tmpCv )
                   tmpCv(-1:nx,-1:ny) = dmb2(0:nx+1,0:ny+1,is,1)
                   call write_quantity( edge_grid,                           &
                       &   val = edge_profiles%ggd( time_sind )%             &
                       &         neutral( js )%state( ks )%density,          &
                       &   value = tmpCv )
                   tmpFace(-1:nx,-1:ny,0) = pfluxm(0:nx+1,0:ny+1,is,1)
                   tmpFace(-1:nx,-1:ny,1) = rfluxm(0:nx+1,0:ny+1,is,1)
                   call write_face_scalar( transport_grid,                   &
                       &   val = edge_transport%model(1)%ggd( time_sind )%   &
                       &         neutral( js )%state( ks )%particles%flux,   &
                       &   value = tmpFace )
                   tmpFace(-1:nx,-1:ny,0) = pefluxm(0:nx+1,0:ny+1,is,1)
                   tmpFace(-1:nx,-1:ny,1) = refluxm(0:nx+1,0:ny+1,is,1)
                   call write_face_scalar( transport_grid,                   &
                       &   val = edge_transport%model(1)%ggd( time_sind )%   &
                       &         neutral( js )%state( ks )%energy%flux,      &
                       &   value = tmpFace )
                   tmpCv(-1:nx,-1:ny) = dmb2(0:nx+1,0:ny+1,is,1)*            &
                       &                tmb2(0:nx+1,0:ny+1,is,1)
                   call write_quantity( edge_grid,                           &
                       &   val = edge_profiles%ggd( time_sind )%             &
                       &         neutral( js )%state( ks )%pressure,         &
                       &   value = tmpCv )
                   if (allocated(um0)) then
                     call write_cell_vector_component( edge_grid,            &
                       &   vectorComponent = edge_profiles%ggd( time_sind )% &
                       &         neutral( js )%state( ks )%velocity,         &
                       &   b2CellData = um0(:,:,0,is),                       &
                       &   vectorID = VEC_ALIGN_POLOIDAL_ID )
                     call write_cell_vector_component( edge_grid,            &
                       &   vectorComponent = edge_profiles%ggd( time_sind )% &
                       &         neutral( js )%state( ks )%velocity,         &
                       &   b2CellData = um0(:,:,1,is),                       &
                       &   vectorID = VEC_ALIGN_RADIAL_ID )
#if ( GGD_MAJOR_VERSION < 2 && GGD_MINOR_VERSION < 13 )
                     call write_cell_vector_component( edge_grid,            &
                       &   vectorComponent = edge_profiles%ggd( time_sind )% &
                       &         neutral( js )%state( ks )%velocity,         &
                       &   b2CellData = um0(:,:,2,is),                       &
                       &   vectorID = VEC_ALIGN_TOROIDAL_ID )
#else
                     call write_cell_vector_component( edge_grid,            &
                       &   vectorComponent = edge_profiles%ggd( time_sind )% &
                       &         neutral( js )%state( ks )%velocity,         &
                       &   b2CellData = um0(:,:,2,is),                       &
                       &   vectorID = VEC_ALIGN_PHI_ID )
#endif
                   end if
                   if (drift_style.eq.0) then
                     tmpCv = 0.0_IDS_real
                     call write_cell_vector_component( edge_grid,            &
                       &   vectorComponent = edge_profiles%ggd( time_sind )% &
                       &                     neutral( js )%state( ks )%      &
                       &                     velocity_diamagnetic,           &
                       &   b2CellData = tmpCv,                               &
                       &   vectorID = VEC_ALIGN_POLOIDAL_ID )
                     call write_cell_vector_component( edge_grid,            &
                       &   vectorComponent = edge_profiles%ggd( time_sind )% &
                       &                     neutral( js )%state( ks )%      &
                       &                     velocity_diamagnetic,           &
                       &   b2CellData = tmpCv,                               &
                       &   vectorID = VEC_ALIGN_RADIAL_ID )
                     call write_cell_vector_component( edge_grid,            &
                       &   vectorComponent = edge_profiles%ggd( time_sind )% &
                       &                     neutral( js )%state( ks )%      &
                       &                     velocity_exb,                   &
                       &   b2CellData = tmpCv,                               &
                       &   vectorID = VEC_ALIGN_POLOIDAL_ID )
                     call write_cell_vector_component( edge_grid,            &
                       &   vectorComponent = edge_profiles%ggd( time_sind )% &
                       &                     neutral( js )%state( ks )%      &
                       &                     velocity_exb,                   &
                       &   b2CellData = tmpCv,                               &
                       &   vectorID = VEC_ALIGN_RADIAL_ID )
                   else
                     tmpFace = 0.0_IDS_real
                     call write_face_vector_component( edge_grid,            &
                       &   vectorComponent = edge_profiles%ggd( time_sind )% &
                       &                     neutral( js )%state( ks )%      &
                       &                     velocity_diamagnetic,           &
                       &   b2FaceData = tmpFace,                             &
                       &   vectorID = VEC_ALIGN_POLOIDAL_ID )
                     call write_face_vector_component( edge_grid,            &
                       &   vectorComponent = edge_profiles%ggd( time_sind )% &
                       &                     neutral( js )%state( ks )%      &
                       &                     velocity_diamagnetic,           &
                       &   b2FaceData = tmpFace,                             &
                       &   vectorID = VEC_ALIGN_RADIAL_ID )
                     call write_face_vector_component( edge_grid,            &
                       &   vectorComponent = edge_profiles%ggd( time_sind )% &
                       &                     neutral( js )%state( ks )%      &
                       &                     velocity_exb,                   &
                       &   b2FaceData = tmpFace,                             &
                       &   vectorID = VEC_ALIGN_POLOIDAL_ID )
                     call write_face_vector_component( edge_grid,            &
                       &   vectorComponent = edge_profiles%ggd( time_sind )% &
                       &                     neutral( js )%state( ks )%      &
                       &                     velocity_exb,                   &
                       &   b2FaceData = tmpFace,                             &
                       &   vectorID = VEC_ALIGN_RADIAL_ID )
                   end if
#else
                   tmpCv(-1:nx,-1:ny) = tmb2(0:nx+1,0:ny+1,is,1)/qe
                   call write_quantity( edge_grid,                           &
                       &   val = profiles_ggd%                               &
                       &         neutral( js )%state( ks )%temperature,      &
                       &   value = tmpCv )
                   tmpCv(-1:nx,-1:ny) = dmb2(0:nx+1,0:ny+1,is,1)
                   call write_quantity( edge_grid,                           &
                       &   val = profiles_ggd%                               &
                       &         neutral( js )%state( ks )%density,          &
                       &   value = tmpCv )
                   tmpFace(-1:nx,-1:ny,0) = pfluxm(0:nx+1,0:ny+1,is,1)
                   tmpFace(-1:nx,-1:ny,1) = rfluxm(0:nx+1,0:ny+1,is,1)
                   call write_face_scalar( transport_grid,                   &
                       &   val = transport_ggd(1)%                           &
                       &         neutral( js )%state( ks )%particles%flux,   &
                       &   value = tmpFace )
                   tmpFace(-1:nx,-1:ny,0) = pefluxm(0:nx+1,0:ny+1,is,1)
                   tmpFace(-1:nx,-1:ny,1) = refluxm(0:nx+1,0:ny+1,is,1)
                   call write_face_scalar( transport_grid,                   &
                       &   val = transport_ggd(1)%                           &
                       &         neutral( js )%state( ks )%energy%flux,      &
                       &   value = tmpFace )
                   tmpCv(-1:nx,-1:ny) = dmb2(0:nx+1,0:ny+1,is,1)*            &
                       &                tmb2(0:nx+1,0:ny+1,is,1)
                   call write_quantity( edge_grid,                           &
                       &   val = profiles_ggd%                               &
                       &         neutral( js )%state( ks )%pressure,         &
                       &   value = tmpCv )
                   if (allocated(um0)) then
                     call write_cell_vector_component( edge_grid,            &
                       &   vectorComponent = profiles_ggd%                   &
                       &         neutral( js )%state( ks )%velocity,         &
                       &   b2CellData = um0(:,:,0,is),                       &
                       &   vectorID = VEC_ALIGN_POLOIDAL_ID )
                     call write_cell_vector_component( edge_grid,            &
                       &   vectorComponent = profiles_ggd%                   &
                       &         neutral( js )%state( ks )%velocity,         &
                       &   b2CellData = um0(:,:,1,is),                       &
                       &   vectorID = VEC_ALIGN_RADIAL_ID )
                     call write_cell_vector_component( edge_grid,            &
                       &   vectorComponent = profiles_ggd%                   &
                       &         neutral( js )%state( ks )%velocity,         &
                       &   b2CellData = um0(:,:,2,is),                       &
                       &   vectorID = VEC_ALIGN_PHI_ID )
                   end if
                   if (drift_style.eq.0) then
                     tmpCv = 0.0_IDS_real
                     call write_cell_vector_component( edge_grid,            &
                       &   vectorComponent = profiles_ggd%                   &
                       &                     neutral( js )%state( ks )%      &
                       &                     velocity_diamagnetic,           &
                       &   b2CellData = tmpCv,                               &
                       &   vectorID = VEC_ALIGN_POLOIDAL_ID )
                     call write_cell_vector_component( edge_grid,            &
                       &   vectorComponent = profiles_ggd%                   &
                       &                     neutral( js )%state( ks )%      &
                       &                     velocity_diamagnetic,           &
                       &   b2CellData = tmpCv,                               &
                       &   vectorID = VEC_ALIGN_RADIAL_ID )
                     call write_cell_vector_component( edge_grid,            &
                       &   vectorComponent = profiles_ggd%                   &
                       &                     neutral( js )%state( ks )%      &
                       &                     velocity_exb,                   &
                       &   b2CellData = tmpCv,                               &
                       &   vectorID = VEC_ALIGN_POLOIDAL_ID )
                     call write_cell_vector_component( edge_grid,            &
                       &   vectorComponent = profiles_ggd%                   &
                       &                     neutral( js )%state( ks )%      &
                       &                     velocity_exb,                   &
                       &   b2CellData = tmpCv,                               &
                       &   vectorID = VEC_ALIGN_RADIAL_ID )
                   else
                     tmpFace = 0.0_IDS_real
                     call write_face_vector_component( edge_grid,            &
                       &   vectorComponent = profiles_ggd%                   &
                       &                     neutral( js )%state( ks )%      &
                       &                     velocity_diamagnetic,           &
                       &   b2FaceData = tmpFace,                             &
                       &   vectorID = VEC_ALIGN_POLOIDAL_ID )
                     call write_face_vector_component( edge_grid,            &
                       &   vectorComponent = profiles_ggd%                   &
                       &                     neutral( js )%state( ks )%      &
                       &                     velocity_diamagnetic,           &
                       &   b2FaceData = tmpFace,                             &
                       &   vectorID = VEC_ALIGN_RADIAL_ID )
                     call write_face_vector_component( edge_grid,            &
                       &   vectorComponent = profiles_ggd%                   &
                       &                     neutral( js )%state( ks )%      &
                       &                     velocity_exb,                   &
                       &   b2FaceData = tmpFace,                             &
                       &   vectorID = VEC_ALIGN_POLOIDAL_ID )
                     call write_face_vector_component( edge_grid,            &
                       &   vectorComponent = profiles_ggd%                   &
                       &                     neutral( js )%state( ks )%      &
                       &                     velocity_exb,                   &
                       &   b2FaceData = tmpFace,                             &
                       &   vectorID = VEC_ALIGN_RADIAL_ID )
                   end if
#endif
#if ( IMAS_MAJOR_VERSION < 4 || ( IMAS_MAJOR_VERSION == 4 && IMAS_MINOR_VERSION < 1 ) )
                   if (balance_netcdf.ne.0) then
                     tmpCv = 0.0_IDS_real
                     do istrai = 1, size(eirene_mc_paml_sna_bal,4)
                       tmpCv(:,:) = tmpCv(:,:)                               &
                          &       + eirene_mc_paml_sna_bal(:,:,is,istrai)    &
                          &       + eirene_mc_pmml_sna_bal(:,:,is,istrai)    &
                          &       + eirene_mc_piml_sna_bal(:,:,is,istrai)
                     end do
                     tmpCv(:,:) = tmpCv(:,:) / vol(:,:)
                     call write_cell_scalar( sources_grid,                   &
                       &   scalar = edge_sources%source(1)%ggd( time_sind )% &
                       &            neutral( js )%state( ks )%particles,     &
                       &   b2CellData = tmpCv )
                     tmpCv = 0.0_IDS_real
                     do istrai = 1, size(eirene_mc_paml_sna_bal,4)
                       tmpCv(:,:) = tmpCv(:,:)                               &
                          &       + eirene_mc_paml_sna_bal(:,:,is,istrai)
                     end do
                     tmpCv(:,:) = tmpCv(:,:) / vol(:,:)
                     call write_cell_scalar( sources_grid,                   &
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
                     call write_cell_scalar( sources_grid,                   &
                       &   scalar = edge_sources%source(6)%ggd( time_sind )% &
                       &            neutral( js )%state( ks )%particles,     &
                       &   b2CellData = tmpCv )
                   end if
                   tmpCv(-1:nx,-1:ny) = emolrad(0:nx+1,0:ny+1,is,0)/vol(-1:nx,-1:ny)
                   call write_cell_scalar( sources_grid,                     &
                       &   scalar = edge_sources%source(12)%                 &
                       &            ggd( time_sind )%neutral( js )%          &
                       &            state( ks )%energy,                      &
                       &   b2CellData = tmpCv )
#else
                   if (balance_netcdf.ne.0) then
                     tmpCv = 0.0_IDS_real
                     do istrai = 1, size(eirene_mc_paml_sna_bal,4)
                       tmpCv(:,:) = tmpCv(:,:)                               &
                          &       + eirene_mc_paml_sna_bal(:,:,is,istrai)    &
                          &       + eirene_mc_pmml_sna_bal(:,:,is,istrai)    &
                          &       + eirene_mc_piml_sna_bal(:,:,is,istrai)
                     end do
                     tmpCv(:,:) = tmpCv(:,:) / vol(:,:)
                     call write_cell_scalar( sources_grid,                   &
                       &   scalar = sources_ggd(1)%                          &
                       &            neutral( js )%state( ks )%particles,     &
                       &   b2CellData = tmpCv )
                     tmpCv = 0.0_IDS_real
                     do istrai = 1, size(eirene_mc_paml_sna_bal,4)
                       tmpCv(:,:) = tmpCv(:,:)                               &
                          &       + eirene_mc_paml_sna_bal(:,:,is,istrai)
                     end do
                     tmpCv(:,:) = tmpCv(:,:) / vol(:,:)
                     call write_cell_scalar( sources_grid,                   &
                       &   scalar = sources_ggd(5)%                          &
                       &            neutral( js )%state( ks )%particles,     &
                       &   b2CellData = tmpCv )
                     tmpCv = 0.0_IDS_real
                     do istrai = 1, size(eirene_mc_pmml_sna_bal,4)
                       tmpCv(:,:) = tmpCv(:,:)                               &
                           &      + eirene_mc_pmml_sna_bal(:,:,is,istrai)    &
                           &      + eirene_mc_piml_sna_bal(:,:,is,istrai)
                     end do
                     tmpCv(:,:) = tmpCv(:,:) / vol(:,:)
                     call write_cell_scalar( sources_grid,                   &
                       &   scalar = sources_ggd(6)%                          &
                       &            neutral( js )%state( ks )%particles,     &
                       &   b2CellData = tmpCv )
                   end if
                   tmpCv(-1:nx,-1:ny) = emolrad(0:nx+1,0:ny+1,is,0)/vol(-1:nx,-1:ny)
                   call write_cell_scalar( sources_grid,                     &
                       &   scalar = sources_ggd(12)%neutral( js )%           &
                       &            state( ks )%energy,                      &
                       &   b2CellData = tmpCv )
#endif
                end do
#endif
            else
                j = 0
                do is = 1, nspecies
                    js = eb2spcr(is)
                    if (.not.is_neutral(js)) cycle
                    j = j + 1
#if ( IMAS_MAJOR_VERSION < 4 || ( IMAS_MAJOR_VERSION == 4 && IMAS_MINOR_VERSION < 1 ) )
                !! na : Fluid neutral density
                    call write_quantity( edge_grid,                           &
                        &   val = edge_profiles%ggd( time_sind )%             &
                        &         neutral( j )%density,                       &
                        &   value = na(:,:,js) )
                    call write_quantity( edge_grid,                           &
                        &   val = edge_profiles%ggd( time_sind )%             &
                        &         neutral( j )%state(1)%density,              &
                        &   value = na(:,:,js) )
                !! ua: Parallel fluid neutral velocity
                    call write_cell_vector_component( edge_grid,              &
                        &   vectorComponent = edge_profiles%ggd( time_sind )% &
                        &                     neutral( j )%velocity,          &
                        &   b2CellData = ua(:,:,js),                          &
                        &   vectorID = VEC_ALIGN_PARALLEL_ID )
                    call write_cell_vector_component( edge_grid,              &
                        &   vectorComponent = edge_profiles%ggd( time_sind )% &
                        &                     neutral( j )%state(1)%velocity, &
                        &   b2CellData = ua(:,:,js),                          &
                        &   vectorID = VEC_ALIGN_PARALLEL_ID )
                    if (drift_style.eq.0) then
                !! wadia: Diamagnetic fluid neutral velocity
                      call write_cell_vector_component( edge_grid,            &
                        &   vectorComponent = edge_profiles%ggd( time_sind )% &
                        &                     neutral( j )%state(1)%          &
                        &                     velocity_diamagnetic,           &
                        &   b2CellData = wadia(:,:,0,0,js),                   &
                        &   vectorID = VEC_ALIGN_POLOIDAL_ID )
                      call write_cell_vector_component( edge_grid,            &
                        &   vectorComponent = edge_profiles%ggd( time_sind )% &
                        &                     neutral( j )%state(1)%          &
                        &                     velocity_diamagnetic,           &
                        &   b2CellData = wadia(:,:,1,1,js),                   &
                        &   vectorID = VEC_ALIGN_RADIAL_ID )
                !! vaecrb: ExB fluid neutral velocity
                      call write_cell_vector_component( edge_grid,            &
                        &   vectorComponent = edge_profiles%ggd( time_sind )% &
                        &                     neutral( j )%state(1)%          &
                        &                     velocity_exb,                   &
                        &   b2CellData = vaecrb(:,:,0,0,js),                  &
                        &   vectorID = VEC_ALIGN_POLOIDAL_ID )
                      call write_cell_vector_component( edge_grid,            &
                        &   vectorComponent = edge_profiles%ggd( time_sind )% &
                        &                     neutral( j )%state(1)%          &
                        &                     velocity_exb,                   &
                        &   b2CellData = vaecrb(:,:,1,1,js),                  &
                        &   vectorID = VEC_ALIGN_RADIAL_ID )
                    else
                !! wadia: Diamagnetic fluid neutral velocity
                      tmpFace(:,:,0) = wadia(:,:,0,0,js)
                      tmpFace(:,:,1) = wadia(:,:,1,0,js)
                      call write_face_vector_component( edge_grid,            &
                        &   vectorComponent = edge_profiles%ggd( time_sind )% &
                        &         neutral( j )%state(1)%velocity_diamagnetic, &
                        &   b2FaceData = tmpFace,                             &
                        &   vectorID = VEC_ALIGN_POLOIDAL_ID )
                      tmpFace(:,:,0) = wadia(:,:,0,1,js)
                      tmpFace(:,:,1) = wadia(:,:,1,1,js)
                      call write_face_vector_component( edge_grid,            &
                        &   vectorComponent = edge_profiles%ggd( time_sind )% &
                        &         neutral( j )%state(1)%velocity_diamagnetic, &
                        &   b2FaceData = tmpFace,                             &
                        &   vectorID = VEC_ALIGN_RADIAL_ID )
                !! vaecrb: ExB fluid neutral velocity
                      tmpFace(:,:,0) = vaecrb(:,:,0,0,js)
                      tmpFace(:,:,1) = vaecrb(:,:,1,0,js)
                      call write_face_vector_component( edge_grid,            &
                        &   vectorComponent = edge_profiles%ggd( time_sind )% &
                        &         neutral( j )%state(1)%velocity_exb,         &
                        &   b2FaceData = tmpFace,                             &
                        &   vectorID = VEC_ALIGN_POLOIDAL_ID )
                      tmpFace(:,:,0) = vaecrb(:,:,0,1,js)
                      tmpFace(:,:,1) = vaecrb(:,:,1,1,js)
                      call write_face_vector_component( edge_grid,            &
                        &   vectorComponent = edge_profiles%ggd( time_sind )% &
                        &         neutral( j )%state(1)%velocity_exb,         &
                        &   b2FaceData = tmpFace,                             &
                        &   vectorID = VEC_ALIGN_RADIAL_ID )
                    end if
                !! fna: Fluid neutral particle flux
#if ( IMAS_MINOR_VERSION > 38 || IMAS_MAJOR_VERSION > 3 )
                    call divide_by_contact_areas(nx,ny,fna(:,:,:,0,js),       &
                        &   flxFace(:,:,:,0))
                    call divide_by_contact_areas(nx,ny,fna(:,:,:,1,js),       &
                        &   flxFace(:,:,:,1))
                    call write_face_scalar( transport_grid,                   &
                        &   val = edge_transport%model(1)%ggd( time_sind )%   &
                        &         neutral( j )%particles%flux_pol,            &
                        &   value = flxFace(:,:,:,0) )
                    call write_face_scalar( transport_grid,                   &
                        &   val = edge_transport%model(1)%ggd( time_sind )%   &
                        &         neutral( j )%particles%flux_radial,         &
                        &   value = flxFace(:,:,:,1) )
                    call write_face_scalar( transport_grid,                   &
                        &   val = edge_transport%model(1)%ggd( time_sind )%   &
                        &         neutral( j )%state(1)%particles%flux_pol,   &
                        &   value = flxFace(:,:,:,0) )
                    call write_face_scalar( transport_grid,                   &
                        &   val = edge_transport%model(1)%ggd( time_sind )%   &
                        &         neutral( j )%state(1)%particles%flux_radial,&
                        &   value = flxFace(:,:,:,1) )
#else
                    tmpFace(:,:,0) = fna(:,:,0,0,js) + fna(:,:,0,1,js)
                    tmpFace(:,:,1) = fna(:,:,1,0,js) + fna(:,:,1,1,js)
                    call divide_by_contact_areas(nx,ny,tmpFace,totFace)
                    call write_face_scalar( transport_grid,                   &
                        &   val = edge_transport%model(1)%ggd( time_sind )%   &
                        &         neutral( j )%particles%flux,                &
                        &   value = totFace )
                    call write_face_scalar( transport_grid,                   &
                        &   val = edge_transport%model(1)%ggd( time_sind )%   &
                        &         neutral( j )%state(1)%particles%flux,       &
                        &   value = totFace )
#endif
                !! pb : Fluid neutral pressure
                    call b2xppb( nx, ny, rza(:,:,js), na(:,:,js), te, ti, pb)
                    call write_quantity( edge_grid,                           &
                        &   val = edge_profiles%ggd( time_sind )%             &
                        &         neutral( j )%pressure,                      &
                        &   value = pb )
                    call write_quantity( edge_grid,                           &
                        &   val = edge_profiles%ggd( time_sind )%             &
                        &         neutral( j )%state(1)%pressure,             &
                        &   value = pb )
                !! cdpa: Fluid neutral diffusivity
#if ( IMAS_MINOR_VERSION > 38 || IMAS_MAJOR_VERSION > 3 )
                    call write_face_scalar( transport_grid,                   &
                        &   val = edge_transport%model(1)%ggd( time_sind )%   &
                        &         neutral( j )%particles%d_pol,               &
                        &   value = cdpa(:,:,:,0,js) )
                    call write_face_scalar( transport_grid,                   &
                        &   val = edge_transport%model(1)%ggd( time_sind )%   &
                        &         neutral( j )%particles%d_radial,            &
                        &   value = cdpa(:,:,:,1,js) )
#else
                    tmpFace(:,:,0) = cdpa(:,:,0,0,js)
                    tmpFace(:,:,1) = cdpa(:,:,1,1,js)
                    call write_face_scalar( transport_grid,                   &
                        &   val = edge_transport%model(1)%ggd( time_sind )%   &
                        &         neutral( j )%particles%d,                   &
                        &   value = tmpFace )
                    call write_face_scalar( transport_grid,                   &
                        &   val = edge_transport%model(1)%ggd( time_sind )%   &
                        &         neutral( j )%state(1)%particles%d,          &
                        &   value = tmpFace )
#endif
                !! Fluid neutral kinetic energy density
                    tmpCv(:,:) = 0.5_IDS_real*am(js)*mp*(ua(:,:,js)**2)*na(:,:,js)
                    call write_quantity( edge_grid,                           &
                        &   val = edge_profiles%ggd( time_sind )%             &
                        &         neutral( j )%energy_density_kinetic,        &
                        &   value = tmpCv )
                    call write_quantity( edge_grid,                           &
                        &   val = edge_profiles%ggd( time_sind )%             &
                        &         neutral( j )%state(1)%                      &
                        &         energy_density_kinetic,                     &
                        &   value = tmpCv )
#else
                !! na : Fluid neutral density
                    call write_quantity( edge_grid,                           &
                        &   val = profiles_ggd%neutral( j )%density,          &
                        &   value = na(:,:,js) )
                    call write_quantity( edge_grid,                           &
                        &   val = profiles_ggd%neutral( j )%state(1)%density, &
                        &   value = na(:,:,js) )
                !! ua: Parallel fluid neutral velocity
                    call write_cell_vector_component( edge_grid,              &
                        &   vectorComponent = profiles_ggd%                   &
                        &                     neutral( j )%velocity,          &
                        &   b2CellData = ua(:,:,js),                          &
                        &   vectorID = VEC_ALIGN_PARALLEL_ID )
                    call write_cell_vector_component( edge_grid,              &
                        &   vectorComponent = profiles_ggd%                   &
                        &                     neutral( j )%state(1)%velocity, &
                        &   b2CellData = ua(:,:,js),                          &
                        &   vectorID = VEC_ALIGN_PARALLEL_ID )
                    if (drift_style.eq.0) then
                !! wadia: Diamagnetic fluid neutral velocity
                      call write_cell_vector_component( edge_grid,            &
                        &   vectorComponent = profiles_ggd%                   &
                        &                     neutral( j )%state(1)%          &
                        &                     velocity_diamagnetic,           &
                        &   b2CellData = wadia(:,:,0,0,js),                   &
                        &   vectorID = VEC_ALIGN_POLOIDAL_ID )
                      call write_cell_vector_component( edge_grid,            &
                        &   vectorComponent = profiles_ggd%                   &
                        &                     neutral( j )%state(1)%          &
                        &                     velocity_diamagnetic,           &
                        &   b2CellData = wadia(:,:,1,1,js),                   &
                        &   vectorID = VEC_ALIGN_RADIAL_ID )
                !! vaecrb: ExB fluid neutral velocity
                      call write_cell_vector_component( edge_grid,            &
                        &   vectorComponent = profiles_ggd%                   &
                        &                     neutral( j )%state(1)%          &
                        &                     velocity_exb,                   &
                        &   b2CellData = vaecrb(:,:,0,0,js),                  &
                        &   vectorID = VEC_ALIGN_POLOIDAL_ID )
                      call write_cell_vector_component( edge_grid,            &
                        &   vectorComponent = profiles_ggd%                   &
                        &                     neutral( j )%state(1)%          &
                        &                     velocity_exb,                   &
                        &   b2CellData = vaecrb(:,:,1,1,js),                  &
                        &   vectorID = VEC_ALIGN_RADIAL_ID )
                    else
                !! wadia: Diamagnetic fluid neutral velocity
                      tmpFace(:,:,0) = wadia(:,:,0,0,js)
                      tmpFace(:,:,1) = wadia(:,:,1,0,js)
                      call write_face_vector_component( edge_grid,            &
                        &   vectorComponent = profiles_ggd%                   &
                        &         neutral( j )%state(1)%velocity_diamagnetic, &
                        &   b2FaceData = tmpFace,                             &
                        &   vectorID = VEC_ALIGN_POLOIDAL_ID )
                      tmpFace(:,:,0) = wadia(:,:,0,1,js)
                      tmpFace(:,:,1) = wadia(:,:,1,1,js)
                      call write_face_vector_component( edge_grid,            &
                        &   vectorComponent = profiles_ggd%                   &
                        &         neutral( j )%state(1)%velocity_diamagnetic, &
                        &   b2FaceData = tmpFace,                             &
                        &   vectorID = VEC_ALIGN_RADIAL_ID )
                !! vaecrb: ExB fluid neutral velocity
                      tmpFace(:,:,0) = vaecrb(:,:,0,0,js)
                      tmpFace(:,:,1) = vaecrb(:,:,1,0,js)
                      call write_face_vector_component( edge_grid,            &
                        &   vectorComponent = profiles_ggd%                   &
                        &         neutral( j )%state(1)%velocity_exb,         &
                        &   b2FaceData = tmpFace,                             &
                        &   vectorID = VEC_ALIGN_POLOIDAL_ID )
                      tmpFace(:,:,0) = vaecrb(:,:,0,1,js)
                      tmpFace(:,:,1) = vaecrb(:,:,1,1,js)
                      call write_face_vector_component( edge_grid,            &
                        &   vectorComponent = profiles_ggd%                   &
                        &         neutral( j )%state(1)%velocity_exb,         &
                        &   b2FaceData = tmpFace,                             &
                        &   vectorID = VEC_ALIGN_RADIAL_ID )
                    end if
                !! fna: Fluid neutral particle flux
                    call divide_by_contact_areas(nx,ny,fna(:,:,:,0,js),       &
                        &   flxFace(:,:,:,0))
                    call divide_by_contact_areas(nx,ny,fna(:,:,:,1,js),       &
                        &   flxFace(:,:,:,1))
                    call write_face_scalar( transport_grid,                   &
                        &   val = transport_ggd(1)%                           &
                        &         neutral( j )%particles%flux_pol,            &
                        &   value = flxFace(:,:,:,0) )
                    call write_face_scalar( transport_grid,                   &
                        &   val = transport_ggd(1)%                           &
                        &         neutral( j )%particles%flux_radial,         &
                        &   value = flxFace(:,:,:,1) )
                    call write_face_scalar( transport_grid,                   &
                        &   val = transport_ggd(1)%                           &
                        &         neutral( j )%state(1)%particles%flux_pol,   &
                        &   value = flxFace(:,:,:,0) )
                    call write_face_scalar( transport_grid,                   &
                        &   val = transport_ggd(1)%                           &
                        &         neutral( j )%state(1)%particles%flux_radial,&
                        &   value = flxFace(:,:,:,1) )
                !! pb : Fluid neutral pressure
                    call b2xppb( nx, ny, rza(:,:,js), na(:,:,js), te, ti, pb)
                    call write_quantity( edge_grid,                           &
                        &   val = profiles_ggd%neutral( j )%pressure,         &
                        &   value = pb )
                    call write_quantity( edge_grid,                           &
                        &   val = profiles_ggd%                               &
                        &         neutral( j )%state(1)%pressure,             &
                        &   value = pb )
                !! cdpa: Fluid neutral diffusivity
                    call write_face_scalar( transport_grid,                   &
                        &   val = transport_ggd(1)%neutral( j )%particles%d_pol, &
                        &   value = cdpa(:,:,:,0,js) )
                    call write_face_scalar( transport_grid,                   &
                        &   val = transport_ggd(1)%neutral( j )%state(1)%     &
                        &         particles%d_pol,                            &
                        &   value = cdpa(:,:,:,0,js) )
                    call write_face_scalar( transport_grid,                   &
                        &   val = transport_ggd(1)%neutral( j )%particles%d_radial, &
                        &   value = cdpa(:,:,:,1,js) )
                        &   val = transport_ggd(1)%neutral( j )%state(1)%     &
                        &         particles%d_radial,                         &
                        &   value = cdpa(:,:,:,1,js) )
                !! Fluid neutral kinetic energy density
                    tmpCv(:,:) = 0.5_IDS_real*am(js)*mp*(ua(:,:,js)**2)*na(:,:,js)
                    call write_quantity( edge_grid,                           &
                        &   val = profiles_ggd%                               &
                        &         neutral( j )%energy_density_kinetic,        &
                        &   value = tmpCv )
                    call write_quantity( edge_grid,                           &
                        &   val = profiles_ggd%neutral( j )%state(1)%         &
                        &         energy_density_kinetic,                     &
                        &   value = tmpCv )
#endif
#if ( IMAS_MAJOR_VERSION < 4 || ( IMAS_MAJOR_VERSION == 4 && IMAS_MINOR_VERSION < 1 ) )
                !! cvsa: Ion diffusivity
                    tmpFace(:,:,0) = cvsa(:,:,0,0,js)
                    tmpFace(:,:,1) = cvsa(:,:,1,1,js)
                    call write_face_vector_component( transport_grid,         &
                        &   vectorComponent = edge_transport%model(1)%        &
                        &                     ggd( time_sind )%neutral( j )%  &
                        &                     momentum%d,                     &
                        &   b2FaceData = tmpFace,                             &
                        &   vectorID = VEC_ALIGN_PARALLEL_ID )
                    call write_face_vector_component( transport_grid,         &
                        &   vectorComponent = edge_transport%model(1)%        &
                        &                     ggd( time_sind )%neutral( j )%  &
                        &                     state(1)%momentum%d,            &
                        &   b2FaceData = tmpFace,                             &
                        &   vectorID = VEC_ALIGN_PARALLEL_ID )
                !! fmo: Ion momentum flux
                    hlpFace(:,:,:) = fmo(:,:,:,0,js) + fmo(:,:,:,1,js)
                    call divide_by_contact_areas(nx,ny,hlpFace,tmpFace)
                    call write_face_vector_component( transport_grid,         &
                        &   vectorComponent = edge_transport%model(1)%        &
                        &                     ggd( time_sind )%neutral( j )%  &
                        &                     momentum%flux,                  &
                        &   b2FaceData = tmpFace,                             &
                        &   vectorID = VEC_ALIGN_PARALLEL_ID )
                    call write_face_vector_component( transport_grid,         &
                        &   vectorComponent = edge_transport%model(1)%        &
                        &                     ggd( time_sind )%neutral( j )%  &
                        &                     state(1)%momentum%flux,         &
                        &   b2FaceData = tmpFace,                             &
                        &   vectorID = VEC_ALIGN_PARALLEL_ID )
                !! fllim0fna: Fluid neutral flux limiter
                    tmpFace(:,:,0) = fllim0fna(:,:,0,0,js)
                    tmpFace(:,:,1) = fllim0fna(:,:,1,1,js)
                    call write_face_scalar( transport_grid,                   &
                        &   val = edge_transport%model(1)%ggd( time_sind )%   &
                        &         neutral( j )%particles%flux_limiter,        &
                        &   value = tmpFace )
                    call write_face_scalar( transport_grid,                   &
                        &   val = edge_transport%model(1)%ggd( time_sind )%   &
                        &         neutral( j )%state(1)%                      &
                        &         particles%flux_limiter,                     &
                        &   value = tmpFace )
                !! fllimvisc: Fluid neutral momentum transport flux limit
                    tmpFace(:,:,0) = fllimvisc(:,:,js)
                    tmpFace(:,:,1) = 1.0_IDS_real
                    call write_face_vector_component( transport_grid,         &
                        &   vectorComponent = edge_transport%model(1)%        &
                        &                     ggd( time_sind )%neutral( j )%  &
                        &                     momentum%flux_limiter,          &
                        &   b2FaceData = fllimvisc(:,:,js),                   &
                        &   vectorID = VEC_ALIGN_PARALLEL_ID )
                    call write_face_vector_component( transport_grid,         &
                        &   vectorComponent = edge_transport%model(1)%        &
                        &                     ggd( time_sind )%neutral( j )%  &
                        &                     state(1)%momentum%flux_limiter, &
                        &   b2FaceData = fllimvisc(:,:,js),                   &
                        &   vectorID = VEC_ALIGN_PARALLEL_ID )
                !! sna: Fluid neutral particle sources
                    tmpCv(:,:) = ( sna(:,:,0,js) +                            &
                        &          sna(:,:,1,js) * na(:,:,js) ) / vol(:,:)
                    call write_cell_scalar( sources_grid,                     &
                        &   scalar = edge_sources%source(1)%                  &
                        &   ggd( time_sind )%neutral( j )%particles,          &
                        &   b2CellData = tmpCv )
                    call write_cell_scalar( sources_grid,                     &
                        &   scalar = edge_sources%source(1)%                  &
                        &   ggd( time_sind )%neutral( j )%state(1)%particles, &
                        &   b2CellData = tmpCv )
                    tmpCv(:,:) = ext_sna(:,:,js) / vol(:,:)
                    call write_cell_scalar( sources_grid,                     &
                        &   scalar = edge_sources%source(2)%ggd( time_sind )% &
                        &            neutral( j )%particles,                  &
                        &   b2CellData = tmpCv )
                    call write_cell_scalar( sources_grid,                     &
                        &   scalar = edge_sources%source(2)%ggd( time_sind )% &
                        &            neutral( j )%state(1)%particles,         &
                        &   b2CellData = tmpCv )
                    tmpCv(:,:) = ( b2stbc_sna(:,:,js) +                       &
                        &          b2stbm_sna(:,:,js) ) / vol(:,:)
                    call write_cell_scalar( sources_grid,                     &
                        &   scalar = edge_sources%source(3)%ggd( time_sind )% &
                        &            neutral( j )%particles,                  &
                        &   b2CellData = tmpCv )
                    call write_cell_scalar( sources_grid,                     &
                        &   scalar = edge_sources%source(3)%ggd( time_sind )% &
                        &            neutral( j )%state(1)%particles,         &
                        &   b2CellData = tmpCv )
                    tmpCv(:,:) = ( snadt(:,:,0,js) +                          &
                        &          snadt(:,:,1,js) * na(:,:,js) ) / vol(:,:)
                    call write_cell_scalar( sources_grid,                     &
                        &   scalar = edge_sources%source(4)%ggd( time_sind )% &
                        &            neutral( j )%particles,                  &
                        &   b2CellData = tmpCv )
                    call write_cell_scalar( sources_grid,                     &
                        &   scalar = edge_sources%source(4)%ggd( time_sind )% &
                        &            neutral( j )%state(1)%particles,         &
                        &   b2CellData = tmpCv )
                    tmpCv(:,:) = b2stbr_sna(:,:,js) / vol(:,:)
                    call write_cell_scalar( sources_grid,                     &
                        &   scalar = edge_sources%source(5)%ggd( time_sind )% &
                        &            neutral( j )%particles,                  &
                        &   b2CellData = tmpCv )
                    call write_cell_scalar( sources_grid,                     &
                        &   scalar = edge_sources%source(5)%ggd( time_sind )% &
                        &            neutral( j )%state(1)%particles,         &
                        &   b2CellData = tmpCv )
                    tmpCv(:,:) = rsana(:,:,js) / vol(:,:)
                    call write_cell_scalar( sources_grid,                     &
                        &   scalar = edge_sources%source(7)%ggd( time_sind )% &
                        &            neutral( j )%particles,                  &
                        &   b2CellData = tmpCv )
                    call write_cell_scalar( sources_grid,                     &
                        &   scalar = edge_sources%source(7)%ggd( time_sind )% &
                        &            neutral( j )%state(1)%particles,         &
                        &   b2CellData = tmpCv )
                    tmpCv(:,:) = rrana(:,:,js) / vol(:,:)
                    call write_cell_scalar( sources_grid,                     &
                        &   scalar = edge_sources%source(8)%ggd( time_sind )% &
                        &            neutral( j )%particles,                  &
                        &   b2CellData = tmpCv )
                    call write_cell_scalar( sources_grid,                     &
                        &   scalar = edge_sources%source(8)%ggd( time_sind )% &
                        &            neutral( j )%state(1)%particles,         &
                        &   b2CellData = tmpCv )
                    tmpCv(:,:) = rcxna(:,:,js) / vol(:,:)
                    call write_cell_scalar( sources_grid,                     &
                        &   scalar = edge_sources%source(9)%ggd( time_sind )% &
                        &            neutral( j )%particles,                  &
                        &   b2CellData = tmpCv )
                    call write_cell_scalar( sources_grid,                     &
                        &   scalar = edge_sources%source(9)%ggd( time_sind )% &
                        &            neutral( j )%state(1)%particles,         &
                        &   b2CellData = tmpCv )
                !! smo: Neutral parallel momentum sources
                    do iy = -1, ny
                      do ix = -1, nx
                        tmpCv(ix,iy) = ( smo(ix,iy,0,js) +                  &
                        &                smo(ix,iy,1,js) * ua(ix,iy,js) +   &
                        &                smo(ix,iy,2,js) * roxa(ix,iy,js) + &
                        &                smo(ix,iy,3,js) * roxa(ix,iy,js)   &
                        &                                * ua(ix,iy,js) ) / vol(ix,iy)
                      end do
                    end do
                    call write_cell_vector_component( sources_grid,           &
                        &   vectorComponent = edge_sources%source(1)%         &
                        &                     ggd( time_sind )%               &
                        &                     neutral( j )%momentum,          &
                        &   b2CellData = tmpCv,                               &
                        &   vectorID = VEC_ALIGN_PARALLEL_ID )
                    call write_cell_vector_component( sources_grid,           &
                        &   vectorComponent = edge_sources%source(1)%         &
                        &                     ggd( time_sind )%neutral( j )%  &
                        &                     state(1)%momentum,              &
                        &   b2CellData = tmpCv,                               &
                        &   vectorID = VEC_ALIGN_PARALLEL_ID )
                    tmpCv(:,:) = ext_smo(:,:,js) / vol(:,:)
                    call write_cell_vector_component( sources_grid,           &
                        &   vectorComponent = edge_sources%source(2)%         &
                        &                     ggd( time_sind )%neutral( j )%  &
                        &                     momentum,                       &
                        &   b2CellData = tmpCv,                               &
                        &   vectorID = VEC_ALIGN_PARALLEL_ID )
                    call write_cell_vector_component( sources_grid,           &
                        &   vectorComponent = edge_sources%source(2)%         &
                        &                     ggd( time_sind )%neutral( j )%  &
                        &                     state(1)%momentum,              &
                        &   b2CellData = tmpCv,                               &
                        &   vectorID = VEC_ALIGN_PARALLEL_ID )
                    tmpCv(:,:) = ( b2stbc_smo(:,:,js) +                       &
                        &          b2stbm_smo(:,:,js) ) / vol(:,:)
                    call write_cell_vector_component( sources_grid,           &
                        &   vectorComponent = edge_sources%source(3)%         &
                        &                     ggd( time_sind )%neutral( j )%  &
                        &                     momentum,                       &
                        &   b2CellData = tmpCv,                               &
                        &   vectorID = VEC_ALIGN_PARALLEL_ID )
                    call write_cell_vector_component( sources_grid,           &
                        &   vectorComponent = edge_sources%source(3)%         &
                        &                     ggd( time_sind )%neutral( j )%  &
                        &                     state(1)%momentum,              &
                        &   b2CellData = tmpCv,                               &
                        &   vectorID = VEC_ALIGN_PARALLEL_ID )
                    do iy = -1, ny
                      do ix = -1, nx
                        tmpCv(ix,iy) = ( smodt(ix,iy,0,js) +                  &
                        &                smodt(ix,iy,1,js) * ua(ix,iy,js) +   &
                        &                smodt(ix,iy,2,js) * roxa(ix,iy,js) + &
                        &                smodt(ix,iy,3,js) * roxa(ix,iy,js)   &
                        &                                  * ua(ix,iy,js) ) / vol(ix,iy)
                      end do
                    end do
                    call write_cell_vector_component( sources_grid,           &
                        &   vectorComponent = edge_sources%source(4)%         &
                        &                     ggd( time_sind )%neutral( j )%  &
                        &                     momentum,                       &
                        &   b2CellData = tmpCv,                               &
                        &   vectorID = VEC_ALIGN_PARALLEL_ID )
                    call write_cell_vector_component( sources_grid,           &
                        &   vectorComponent = edge_sources%source(4)%         &
                        &                     ggd( time_sind )%neutral( j )%  &
                        &                     state(1)%momentum,              &
                        &   b2CellData = tmpCv,                               &
                        &   vectorID = VEC_ALIGN_PARALLEL_ID )
                    tmpCv(:,:) = b2stbr_smo(:,:,js) / vol(:,:)
                    call write_cell_vector_component( sources_grid,           &
                        &   vectorComponent = edge_sources%source(5)%         &
                        &            ggd( time_sind )%neutral( j )%momentum,  &
                        &   b2CellData = tmpCv,                               &
                        &   vectorID = VEC_ALIGN_PARALLEL_ID )
                    call write_cell_vector_component( sources_grid,           &
                        &   vectorComponent = edge_sources%source(5)%         &
                        &            ggd( time_sind )%neutral( j )%           &
                        &            state(1)%momentum,                       &
                        &   b2CellData = tmpCv,                               &
                        &   vectorID = VEC_ALIGN_PARALLEL_ID )
                    tmpCv(:,:) = rsamo(:,:,js) / vol(:,:)
                    call write_cell_vector_component( sources_grid,           &
                        &   vectorComponent = edge_sources%source(7)%         &
                        &                     ggd( time_sind )%neutral( j )%  &
                        &                     momentum,                       &
                        &   b2CellData = tmpCv,                               &
                        &   vectorID = VEC_ALIGN_PARALLEL_ID )
                    call write_cell_vector_component( sources_grid,           &
                        &   vectorComponent = edge_sources%source(7)%         &
                        &                     ggd( time_sind )%neutral( j )%  &
                        &                     state(1)%momentum,              &
                        &   b2CellData = tmpCv,                               &
                        &   vectorID = VEC_ALIGN_PARALLEL_ID )
                    tmpCv(:,:) = rramo(:,:,js) / vol(:,:)
                    call write_cell_vector_component( sources_grid,           &
                        &   vectorComponent = edge_sources%source(8)%         &
                        &                     ggd( time_sind )%neutral( j )%  &
                        &                     momentum,                       &
                        &   b2CellData = tmpCv,                               &
                        &   vectorID = VEC_ALIGN_PARALLEL_ID )
                    call write_cell_vector_component( sources_grid,           &
                        &   vectorComponent = edge_sources%source(8)%         &
                        &                     ggd( time_sind )%neutral( j )%  &
                        &                     state(1)%momentum,              &
                        &   b2CellData = tmpCv,                               &
                        &   vectorID = VEC_ALIGN_PARALLEL_ID )
                    tmpCv(:,:) = rcxmo(:,:,js) / vol(:,:)
                    call write_cell_vector_component( sources_grid,           &
                        &   vectorComponent = edge_sources%source(9)%         &
                        &                     ggd( time_sind )%neutral( j )%  &
                        &                     momentum,                       &
                        &   b2CellData = tmpCv,                               &
                        &   vectorID = VEC_ALIGN_PARALLEL_ID )
                    call write_cell_vector_component( sources_grid,           &
                        &   vectorComponent = edge_sources%source(9)%         &
                        &                     ggd( time_sind )%neutral( j )%  &
                        &                     state(1)%momentum,              &
                        &   b2CellData = tmpCv,                               &
                        &   vectorID = VEC_ALIGN_PARALLEL_ID )
#else
                !! cvsa: Ion diffusivity
                    tmpFace(:,:,0) = cvsa(:,:,0,0,js)
                    tmpFace(:,:,1) = cvsa(:,:,1,1,js)
                    call write_face_vector_component( transport_grid,         &
                        &   vectorComponent = transport_ggd(1)%neutral( j )%  &
                        &                     momentum%d,                     &
                        &   b2FaceData = tmpFace,                             &
                        &   vectorID = VEC_ALIGN_PARALLEL_ID )
                    call write_face_vector_component( transport_grid,         &
                        &   vectorComponent = transport_ggd(1)%neutral( j )%  &
                        &                     state(1)%momentum%d,            &
                        &   b2FaceData = tmpFace,                             &
                        &   vectorID = VEC_ALIGN_PARALLEL_ID )
                !! fmo: Ion momentum flux
                    hlpFace(:,:,:) = fmo(:,:,:,0,js) + fmo(:,:,:,1,js)
                    call divide_by_contact_areas(nx,ny,hlpFace,tmpFace)
                    call write_face_vector_component( transport_grid,         &
                        &   vectorComponent = transport_ggd(1)%neutral( j )%  &
                        &                     momentum%flux,                  &
                        &   b2FaceData = tmpFace,                             &
                        &   vectorID = VEC_ALIGN_PARALLEL_ID )
                    call write_face_vector_component( transport_grid,         &
                        &   vectorComponent = transport_ggd(1)%neutral( j )%  &
                        &                     state(1)%momentum%flux,         &
                        &   b2FaceData = tmpFace,                             &
                        &   vectorID = VEC_ALIGN_PARALLEL_ID )
                !! fllim0fna: Fluid neutral flux limiter
                    tmpFace(:,:,0) = fllim0fna(:,:,0,0,js)
                    tmpFace(:,:,1) = fllim0fna(:,:,1,1,js)
                    call write_face_scalar( transport_grid,                   &
                        &   val = transport_ggd(1)%                           &
                        &         neutral( j )%particles%flux_limiter,        &
                        &   value = tmpFace )
                    call write_face_scalar( transport_grid,                   &
                        &   val = transport_ggd(1)%neutral( j )%state(1)%     &
                        &         particles%flux_limiter,                     &
                        &   value = tmpFace )
                !! fllimvisc: Fluid neutral momentum transport flux limit
                    tmpFace(:,:,0) = fllimvisc(:,:,js)
                    tmpFace(:,:,1) = 1.0_IDS_real
                    call write_face_vector_component( transport_grid,         &
                        &   vectorComponent = transport_ggd(1)%neutral( j )%  &
                        &                     momentum%flux_limiter,          &
                        &   b2FaceData = fllimvisc(:,:,js),                   &
                        &   vectorID = VEC_ALIGN_PARALLEL_ID )
                    call write_face_vector_component( transport_grid,         &
                        &   vectorComponent = transport_ggd(1)%neutral( j )%  &
                        &                     state(1)%momentum%flux_limiter, &
                        &   b2FaceData = fllimvisc(:,:,js),                   &
                        &   vectorID = VEC_ALIGN_PARALLEL_ID )
                !! sna: Fluid neutral particle sources
                    tmpCv(:,:) = ( sna(:,:,0,js) +                            &
                        &          sna(:,:,1,js) * na(:,:,js) ) / vol(:,:)
                    call write_cell_scalar( sources_grid,                     &
                        &   scalar = sources_ggd(1)%neutral( j )%particles,   &
                        &   b2CellData = tmpCv )
                    call write_cell_scalar( sources_grid,                     &
                        &   scalar = sources_ggd(1)%neutral( j )%state(1)%    &
                        &            particles, &
                        &   b2CellData = tmpCv )
                    tmpCv(:,:) = ext_sna(:,:,js) / vol(:,:)
                    call write_cell_scalar( sources_grid,                     &
                        &   scalar = sources_ggd(2)%neutral( j )%particles,   &
                        &   b2CellData = tmpCv )
                    call write_cell_scalar( sources_grid,                     &
                        &   scalar = sources_ggd(2)%                          &
                        &            neutral( j )%state(1)%particles,         &
                        &   b2CellData = tmpCv )
                    tmpCv(:,:) = ( b2stbc_sna(:,:,js) +                       &
                        &          b2stbm_sna(:,:,js) ) / vol(:,:)
                    call write_cell_scalar( sources_grid,                     &
                        &   scalar = sources_ggd(3)%neutral( j )%particles,   &
                        &   b2CellData = tmpCv )
                    call write_cell_scalar( sources_grid,                     &
                        &   scalar = sources_ggd(3)%                          &
                        &            neutral( j )%state(1)%particles,         &
                        &   b2CellData = tmpCv )
                    tmpCv(:,:) = ( snadt(:,:,0,js) +                          &
                        &          snadt(:,:,1,js) * na(:,:,js) ) / vol(:,:)
                    call write_cell_scalar( sources_grid,                     &
                        &   scalar = sources_ggd(4)%neutral( j )%particles,   &
                        &   b2CellData = tmpCv )
                    call write_cell_scalar( sources_grid,                     &
                        &   scalar = sources_ggd(4)%                          &
                        &            neutral( j )%state(1)%particles,         &
                        &   b2CellData = tmpCv )
                    tmpCv(:,:) = b2stbr_sna(:,:,js) / vol(:,:)
                    call write_cell_scalar( sources_grid,                     &
                        &   scalar = sources_ggd(5)%neutral( j )%particles,   &
                        &   b2CellData = tmpCv )
                    call write_cell_scalar( sources_grid,                     &
                        &   scalar = sources_ggd(5)%                          &
                        &            neutral( j )%state(1)%particles,         &
                        &   b2CellData = tmpCv )
                    tmpCv(:,:) = rsana(:,:,js) / vol(:,:)
                    call write_cell_scalar( sources_grid,                     &
                        &   scalar = sources_ggd(7)%neutral( j )%particles,   &
                        &   b2CellData = tmpCv )
                    call write_cell_scalar( sources_grid,                     &
                        &   scalar = sources_ggd(7)%                          &
                        &            neutral( j )%state(1)%particles,         &
                        &   b2CellData = tmpCv )
                    tmpCv(:,:) = rrana(:,:,js) / vol(:,:)
                    call write_cell_scalar( sources_grid,                     &
                        &   scalar = sources_ggd(8)%neutral( j )%particles,   &
                        &   b2CellData = tmpCv )
                    call write_cell_scalar( sources_grid,                     &
                        &   scalar = sources_ggd(8)%                          &
                        &            neutral( j )%state(1)%particles,         &
                        &   b2CellData = tmpCv )
                    tmpCv(:,:) = rcxna(:,:,js) / vol(:,:)
                    call write_cell_scalar( sources_grid,                     &
                        &   scalar = sources_ggd(9)%neutral( j )%particles,   &
                        &   b2CellData = tmpCv )
                    call write_cell_scalar( sources_grid,                     &
                        &   scalar = sources_ggd(9)%                          &
                        &            neutral( j )%state(1)%particles,         &
                        &   b2CellData = tmpCv )
                !! smo: Neutral parallel momentum sources
                    do iy = -1, ny
                      do ix = -1, nx
                        tmpCv(ix,iy) = ( smo(ix,iy,0,js) +                  &
                        &                smo(ix,iy,1,js) * ua(ix,iy,js) +   &
                        &                smo(ix,iy,2,js) * roxa(ix,iy,js) + &
                        &                smo(ix,iy,3,js) * roxa(ix,iy,js)   &
                        &                                * ua(ix,iy,js) ) / vol(ix,iy)
                      end do
                    end do
                    call write_cell_vector_component( sources_grid,           &
                        &   vectorComponent = sources_ggd(1)%                 &
                        &                     neutral( j )%momentum,          &
                        &   b2CellData = tmpCv,                               &
                        &   vectorID = VEC_ALIGN_PARALLEL_ID )
                    call write_cell_vector_component( sources_grid,           &
                        &   vectorComponent = sources_ggd(1)%neutral( j )%    &
                        &                     state(1)%momentum,              &
                        &   b2CellData = tmpCv,                               &
                        &   vectorID = VEC_ALIGN_PARALLEL_ID )
                    tmpCv(:,:) = ext_smo(:,:,js) / vol(:,:)
                    call write_cell_vector_component( sources_grid,           &
                        &   vectorComponent = sources_ggd(2)%neutral( j )%    &
                        &                     momentum,                       &
                        &   b2CellData = tmpCv,                               &
                        &   vectorID = VEC_ALIGN_PARALLEL_ID )
                    call write_cell_vector_component( sources_grid,           &
                        &   vectorComponent = sources_ggd(2)%neutral( j )%    &
                        &                     state(1)%momentum,              &
                        &   b2CellData = tmpCv,                               &
                        &   vectorID = VEC_ALIGN_PARALLEL_ID )
                    tmpCv(:,:) = ( b2stbc_smo(:,:,js) +                       &
                        &          b2stbm_smo(:,:,js) ) / vol(:,:)
                    call write_cell_vector_component( sources_grid,           &
                        &   vectorComponent = sources_ggd(3)%neutral( j )%    &
                        &                     momentum,                       &
                        &   b2CellData = tmpCv,                               &
                        &   vectorID = VEC_ALIGN_PARALLEL_ID )
                    call write_cell_vector_component( sources_grid,           &
                        &   vectorComponent = sources_ggd(3)%neutral( j )%    &
                        &                     state(1)%momentum,              &
                        &   b2CellData = tmpCv,                               &
                        &   vectorID = VEC_ALIGN_PARALLEL_ID )
                    do iy = -1, ny
                      do ix = -1, nx
                        tmpCv(ix,iy) = ( smodt(ix,iy,0,js) +                  &
                        &                smodt(ix,iy,1,js) * ua(ix,iy,js) +   &
                        &                smodt(ix,iy,2,js) * roxa(ix,iy,js) + &
                        &                smodt(ix,iy,3,js) * roxa(ix,iy,js)   &
                        &                                  * ua(ix,iy,js) ) / vol(ix,iy)
                      end do
                    end do
                    call write_cell_vector_component( sources_grid,           &
                        &   vectorComponent = sources_ggd(4)%neutral( j )%    &
                        &                     momentum,                       &
                        &   b2CellData = tmpCv,                               &
                        &   vectorID = VEC_ALIGN_PARALLEL_ID )
                    call write_cell_vector_component( sources_grid,           &
                        &   vectorComponent = sources_ggd(4)%neutral( j )%    &
                        &                     state(1)%momentum,              &
                        &   b2CellData = tmpCv,                               &
                        &   vectorID = VEC_ALIGN_PARALLEL_ID )
                    tmpCv(:,:) = b2stbr_smo(:,:,js) / vol(:,:)
                    call write_cell_vector_component( sources_grid,           &
                        &   vectorComponent = sources_ggd(5)%neutral( j )%    &
                        &                     momentum,                       &
                        &   b2CellData = tmpCv,                               &
                        &   vectorID = VEC_ALIGN_PARALLEL_ID )
                    call write_cell_vector_component( sources_grid,           &
                        &   vectorComponent = sources_ggd(5)%neutral( j )%    &
                        &                     state(1)%momentum,              &
                        &   b2CellData = tmpCv,                               &
                        &   vectorID = VEC_ALIGN_PARALLEL_ID )
                    tmpCv(:,:) = rsamo(:,:,js) / vol(:,:)
                    call write_cell_vector_component( sources_grid,           &
                        &   vectorComponent = sources_ggd(7)%neutral( j )%    &
                        &                     momentum,                       &
                        &   b2CellData = tmpCv,                               &
                        &   vectorID = VEC_ALIGN_PARALLEL_ID )
                    call write_cell_vector_component( sources_grid,           &
                        &   vectorComponent = sources_ggd(7)%neutral( j )%    &
                        &                     state(1)%momentum,              &
                        &   b2CellData = tmpCv,                               &
                        &   vectorID = VEC_ALIGN_PARALLEL_ID )
                    tmpCv(:,:) = rramo(:,:,js) / vol(:,:)
                    call write_cell_vector_component( sources_grid,           &
                        &   vectorComponent = sources_ggd(8)%neutral( j )%    &
                        &                     momentum,                       &
                        &   b2CellData = tmpCv,                               &
                        &   vectorID = VEC_ALIGN_PARALLEL_ID )
                    call write_cell_vector_component( sources_grid,           &
                        &   vectorComponent = sources_ggd(8)%neutral( j )%    &
                        &                     state(1)%momentum,              &
                        &   b2CellData = tmpCv,                               &
                        &   vectorID = VEC_ALIGN_PARALLEL_ID )
                    tmpCv(:,:) = rcxmo(:,:,js) / vol(:,:)
                    call write_cell_vector_component( sources_grid,           &
                        &   vectorComponent = sources_ggd(9)%neutral( j )%    &
                        &                     momentum,                       &
                        &   b2CellData = tmpCv,                               &
                        &   vectorID = VEC_ALIGN_PARALLEL_ID )
                    call write_cell_vector_component( sources_grid,           &
                        &   vectorComponent = sources_ggd(9)%neutral( j )%    &
                        &                     state(1)%momentum,              &
                        &   b2CellData = tmpCv,                               &
                        &   vectorID = VEC_ALIGN_PARALLEL_ID )
#endif
                end do
            end if

#if ( IMAS_MAJOR_VERSION < 4 || ( IMAS_MAJOR_VERSION == 4 && IMAS_MINOR_VERSION < 1 ) )
            !! po: Electric potential
            call write_quantity( edge_grid,                             &
                &   val = edge_profiles%ggd( time_sind )%phi_potential, &
                &   value = po )
            !! Current sources
            tmpCv(:,:) = ( sch(:,:,0) + sch(:,:,1) * po(:,:) +          &
                &          sch(:,:,2) * ne(:,:) +                       &
                &          sch(:,:,3) * ne(:,:) * po(:,:) ) / vol(:,:)
            call write_cell_scalar( sources_grid,                       &
                &   scalar = edge_sources%source(1)%ggd( time_sind )%   &
                &            current,                                   &
                &   b2CellData = tmpCv )
            tmpCv(:,:) = ext_sch(:,:) / vol(:,:)
            call write_cell_scalar( sources_grid,                       &
                &   scalar = edge_sources%source(2)%ggd( time_sind )%   &
                &            current,                                   &
                &   b2CellData = tmpCv )
            tmpCv(:,:) = ( b2stbc_sch(:,:) + b2stbm_sch(:,:) ) / vol(:,:)
            call write_cell_scalar( sources_grid,                       &
                &   scalar = edge_sources%source(3)%ggd( time_sind )%   &
                &            current,                                   &
                &   b2CellData = tmpCv )
            tmpCv(:,:) = ( schdt(:,:,0) + schdt(:,:,1) * po(:,:) +      &
                &          schdt(:,:,2) * ne(:,:) +                     &
                &          schdt(:,:,3) * ne(:,:) * po(:,:) ) / vol(:,:)
            call write_cell_scalar( sources_grid,                       &
                &   scalar = edge_sources%source(4)%ggd( time_sind )%   &
                &            current,                                   &
                &   b2CellData = tmpCv )
            tmpCv(:,:) = b2stbr_sch(:,:) / vol(:,:)
            call write_cell_scalar( sources_grid,                       &
                &   scalar = edge_sources%source(5)%ggd( time_sind )%   &
                &            current,                                   &
                &   b2CellData = tmpCv )
#if ( ( IMAS_MINOR_VERSION > 38 || IMAS_MAJOR_VERSION > 3 ) && defined(B25_EIRENE) )
            if (use_eirene.ne.0) then
              tmpCv(:,:) = sch0_eir_tot(:,:) / vol(:,:)
              call write_cell_scalar( sources_grid,                     &
                &   scalar = edge_sources%source(13)%ggd( time_sind )%  &
                &            current,                                   &
                &   b2CellData = tmpCv )
            end if
#endif
            !! csig : Electric conductivity
            tmpFace(:,:,0) = csig(:,:,0,0)
            tmpFace(:,:,1) = csig(:,:,1,0)
            call write_face_vector_component( transport_grid,           &
                &   vectorComponent = edge_transport%model(1)%          &
                &         ggd( time_sind )%conductivity,                &
                &   b2FaceData = tmpFace,                               &
                &   vectorID = VEC_ALIGN_POLOIDAL_ID )
            tmpFace(:,:,0) = csig(:,:,0,1)
            tmpFace(:,:,1) = csig(:,:,1,1)
            call write_face_vector_component( transport_grid,           &
                &   vectorComponent = edge_transport%model(1)%          &
                &         ggd( time_sind )%conductivity,                &
                &   b2FaceData = tmpFace,                               &
                &   vectorID = VEC_ALIGN_RADIAL_ID )
#else
            !! po: Electric potential
            call write_quantity( edge_grid,                             &
                &   val = profiles_ggd%phi_potential,                   &
                &   value = po )
            !! Current sources
            tmpCv(:,:) = ( sch(:,:,0) + sch(:,:,1) * po(:,:) +          &
                &          sch(:,:,2) * ne(:,:) +                       &
                &          sch(:,:,3) * ne(:,:) * po(:,:) ) / vol(:,:)
            call write_cell_scalar( sources_grid,                       &
                &   scalar = sources_ggd(1)%current,                    &
                &   b2CellData = tmpCv )
            tmpCv(:,:) = ext_sch(:,:) / vol(:,:)
            call write_cell_scalar( sources_grid,                       &
                &   scalar = sources_ggd(2)%current,                    &
                &   b2CellData = tmpCv )
            tmpCv(:,:) = ( b2stbc_sch(:,:) + b2stbm_sch(:,:) ) / vol(:,:)
            call write_cell_scalar( sources_grid,                       &
                &   scalar = sources_ggd(3)%current,                    &
                &   b2CellData = tmpCv )
            tmpCv(:,:) = ( schdt(:,:,0) + schdt(:,:,1) * po(:,:) +      &
                &          schdt(:,:,2) * ne(:,:) +                     &
                &          schdt(:,:,3) * ne(:,:) * po(:,:) ) / vol(:,:)
            call write_cell_scalar( sources_grid,                       &
                &   scalar = sources_ggd(4)%current,                    &
                &   b2CellData = tmpCv )
            tmpCv(:,:) = b2stbr_sch(:,:) / vol(:,:)
            call write_cell_scalar( sources_grid,                       &
                &   scalar = sources_ggd(5)%current,                    &
                &   b2CellData = tmpCv )
#ifdef B25_EIRENE
            if (use_eirene.ne.0) then
              tmpCv(:,:) = sch0_eir_tot(:,:) / vol(:,:)
              call write_cell_scalar( sources_grid,                     &
                &   scalar = sources_ggd(13)%current,                   &
                &   b2CellData = tmpCv )
            end if
#endif
            !! csig : Electric conductivity
            tmpFace(:,:,0) = csig(:,:,0,0)
            tmpFace(:,:,1) = csig(:,:,1,0)
            call write_face_vector_component( transport_grid,           &
                &   vectorComponent = transport_ggd(1)%conductivity,    &
                &   b2FaceData = tmpFace,                               &
                &   vectorID = VEC_ALIGN_POLOIDAL_ID )
            tmpFace(:,:,0) = csig(:,:,0,1)
            tmpFace(:,:,1) = csig(:,:,1,1)
            call write_face_vector_component( transport_grid,           &
                &   vectorComponent = transport_ggd(1)%conductivity,    &
                &   b2FaceData = tmpFace,                               &
                &   vectorID = VEC_ALIGN_RADIAL_ID )
#endif

            !! B (magnetic field vector)
            !! Compute unit basis vectors along the field directions
            call compute_Coordinate_Unit_Vectors(crx, cry, e(:,:,:,1),  &
                &   e(:,:,:,2), e(:,:,:,3))

#if ( IMAS_MINOR_VERSION > 21 || IMAS_MAJOR_VERSION > 3 )
            !! write the emissivity data
            !! Process 1. Line and recombination radiation from B2.5 ions
            j = 0
            do is = 1, nspecies
              js = eb2spcr(is)
              if (is_neutral(js) .and. use_eirene.eq.0) then
                j = j + 1
                tmpCv(:,:) = rqrad(:,:,js) / vol(:,:)
                call write_cell_scalar( radiation_grid,                    &
                    &   scalar = radiation%process(1)%                     &
                    &   ggd( time_sind )%neutral( j )%emissivity,          &
                    &   b2CellData = tmpCv )
#if ( IMAS_MINOR_VERSION > 24 || IMAS_MAJOR_VERSION > 3 )
                call write_cell_scalar( radiation_grid,                    &
                    &   scalar = radiation%process(1)%                     &
                    &   ggd( time_sind )%neutral( j )%                     &
                    &   state(1)%emissivity,                               &
                    &   b2CellData = tmpCv )
#endif
              end if
              totCv(:,:) = 0.0_IDS_real
              do js = 1, istion(is)
                tmpCv(:,:) = rqrad(:,:,ispion(is,js)) / vol(:,:)
                totCv(:,:) = totCv(:,:) + tmpCv(:,:)
#if ( IMAS_MINOR_VERSION > 24 || IMAS_MAJOR_VERSION > 3 )
                call write_cell_scalar( radiation_grid,                    &
                    &   scalar = radiation%process(1)%                     &
                    &   ggd( time_sind )%ion( is )%state( js )%emissivity, &
                    &   b2CellData = tmpCv )
#endif
              end do
              call write_cell_scalar( radiation_grid,                      &
                    &   scalar = radiation%process(1)%                     &
                    &   ggd( time_sind )%ion( is )%emissivity,             &
                    &   b2CellData = totCv )
            end do
            !! Process 2. Bremsstrahlung from B2.5 ions
            totCv(:,:) = 0.0_IDS_real
            do is = 1, nspecies
              do js = 1, istion(is)
                tmpCv(:,:) = rqbrm(:,:,ispion(is,js)) / vol(:,:)
                totCv(:,:) = totCv(:,:) + tmpCv(:,:)
#if ( IMAS_MINOR_VERSION > 24 || IMAS_MAJOR_VERSION > 3 )
                call write_cell_scalar( radiation_grid,                    &
                    &   scalar = radiation%process(2)%                     &
                    &   ggd( time_sind )%ion( is )%state( js )%emissivity, &
                    &   b2CellData = tmpCv )
#endif
              end do
              call write_cell_scalar( radiation_grid,                      &
                  &   scalar = radiation%process(2)%                       &
                  &   ggd( time_sind )%ion( is )%emissivity,               &
                  &   b2CellData = totCv )
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
                call write_cell_scalar( radiation_grid,                    &
                    &   scalar = radiation%process(3)%                     &
                    &   ggd( time_sind )%neutral( js )%emissivity,         &
                    &   b2CellData = tmpCv )
              end do
              do is = 1, natmi
                js = latmscl(is)
                ks = isstat(is)
                tmpCv(-1:nx,-1:ny)=-eneutrad(0:nx+1,0:ny+1,is,0) / vol(-1:nx,-1:ny)
#if ( IMAS_MINOR_VERSION > 24 || IMAS_MAJOR_VERSION > 3 )
                call write_cell_scalar( radiation_grid,                    &
                    &   scalar = radiation%process(3)%ggd( time_sind )%    &
                    &   neutral( js )%state( ks )%emissivity,              &
                    &   b2CellData = tmpCv )
#endif
              end do
              do js = nspecies+1, nneut
                 tmpCv(:,:) = 0.0_IDS_real
                 do is = 1, nmoli
                    if (imneut(is).eq.js) then
                      do ix = -1, nx
                        do iy = -1, ny
                           tmpCv(ix,iy) = tmpCv(ix,iy)-emolrad(ix+1,iy+1,is,0)
                        end do
                      end do
                    end if
                 end do
                 tmpCv(:,:) = tmpCv(:,:) / vol(:,:)
                 call write_cell_scalar( radiation_grid,                      &
                    &   scalar = radiation%process(3)%                        &
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
#if ( IMAS_MINOR_VERSION > 24 || IMAS_MAJOR_VERSION > 3 )
                call write_cell_scalar( radiation_grid,                      &
                    &   scalar = radiation%process(3)%ggd( time_sind )%      &
                    &   neutral( js )%state( ks )%emissivity,                &
                    &   b2CellData = tmpCv )
#endif
              end do
              !! Process 4. Eirene molecular ions
              do is = nspecies+1, nsion
                totCv(:,:) = 0.0_R8
                do js = 1, istion(is)
                  do ix = -1, nx
                    do iy = -1, ny
                      tmpCv(ix,iy) = -eionrad(ix+1,iy+1,ispion(is,js),0) / vol(ix,iy)
                    end do
                  end do
                  totCv(:,:) = totCv(:,:) + tmpCv(:,:)
#if ( IMAS_MINOR_VERSION > 24 || IMAS_MAJOR_VERSION > 3 )
                  call write_cell_scalar( radiation_grid,                     &
                      &   scalar = radiation%process(4)%                      &
                      &   ggd( time_sind )%ion( is )%state( js )%emissivity,  &
                      &   b2CellData = tmpCv )
#endif
                end do
                call write_cell_scalar( radiation_grid,                       &
                    &   scalar = radiation%process(4)%                        &
                    &   ggd( time_sind )%ion( is )%emissivity,                &
                    &   b2CellData = totCv )
              end do
            end if
#endif
#endif
#else
            call logmsg( LOGINFO, &
            &   "b2mod_ual_io.B25_process_ids: GGD not available, no plasma state writing" )
#endif
        end if

! Copy in profiles, sources and transport data
#if ( IMAS_MAJOR_VERSION > 4 || ( IMAS_MAJOR_VERSION == 4 && IMAS_MINOR_VERSION > 0 ) )
        edge_profiles%ggd( time_sind ) = profiles_ggd
        plasma_profiles%ggd( time_sind ) = profiles_ggd
        do i = 1, nsources
          edge_sources%source(i)%ggd( time_sind ) = sources_ggd(i)
          plasma_sources%source(i)%ggd( time_sind ) = sources_ggd(i)
        end do
        edge_transport%model(1)%ggd( time_sind ) = transport_ggd(1)
        plasma_transport%model(1)%ggd( time_sind ) = transport_ggd(1)
#endif

#if ( IMAS_MINOR_VERSION > 21 || IMAS_MAJOR_VERSION > 3 )
! Summary separatrix data
        if (maxval(abs(fpsi(-1:nx,-1:ny,0:3))).gt.0.0_R8) then
           allocate( summary%local%separatrix%position%psi( num_time_slices ) )
           summary%local%separatrix%position%psi( time_sind ) = fpsi(jxa,jsep,2)
#if ( IMAS_MINOR_VERSION > 36 || IMAS_MAJOR_VERSION > 3 )
           allocate( summary%local%separatrix_average%position%psi( num_time_slices ) )
           summary%local%separatrix_average%position%psi( time_sind ) = fpsi(jxa,jsep,2)
#endif
        end if
        call write_sourced_value( summary%local%separatrix%t_e, &
           &  0.5_R8 * (te(jxa,jsep)+ te(topix(jxa,jsep),topiy(jxa,jsep)))/ev )
        call write_sourced_value( summary%local%separatrix%t_i_average, &
           &  0.5_R8 * (ti(jxa,jsep)+ ti(topix(jxa,jsep),topiy(jxa,jsep)))/ev )
        call write_sourced_value( summary%local%separatrix%n_e, &
           &  0.5_R8 * (ne(jxa,jsep)+ ne(topix(jxa,jsep),topiy(jxa,jsep))) )
#if ( IMAS_MINOR_VERSION > 36 || IMAS_MAJOR_VERSION > 3 )
        tmpCv = 1.0_IDS_real
        u = separatrix_average( te, tmpCv )
        call write_sourced_value( summary%local%separatrix_average%t_e, u/ev )
        tmpCv = 1.0_IDS_real
        u = separatrix_average( ti, tmpCv )
        call write_sourced_value( summary%local%separatrix_average%t_i_average, u/ev )
        tmpCv = 1.0_IDS_real
        u = separatrix_average( ne, tmpCv )
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
              &  0.5_R8 * (na(jxa,jsep,i) + na(topix(jxa,jsep),topiy(jxa,jsep),i))
            nasum = nasum + na(jxa,jsep,i)
            vtor = vtor + ua(jxa,jsep,i)*abs(bb(jxa,jsep,2)/bb(jxa,jsep,3)) &
              &  * na(jxa,jsep,i)
          end do
          if (nasum.gt.0.0_R8) vtor = vtor / nasum
          totCv = 0.0_IDS_real
          tmpVx = 0.0_IDS_real
          tmpCv = 1.0_IDS_real
          do i = is1, is2
            tmpVx(:,:) = tmpVx(:,:) + na(:,:,i)
            totCv(:,:) = totCv(:,:) + ua(:,:,1)*na(:,:,i)*abs(bb(:,:,2)/bb(:,:,3))
          end do
          if (nasum.gt.0.0_R8) totCv(:,:) = totCv(:,:)/tmpVx(:,:)
          u = separatrix_average( tmpVx, tmpCv )
          v = separatrix_average( totCv, tmpCv )
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
#if ( IMAS_MAJOR_VERSION > 4 || ( IMAS_MAJOR_VERSION == 4 && IMAS_MINOR_VERSION > 0 ) )
          case ('DT')
            call write_sourced_value( summary%local%separatrix%n_i%deuterium_tritium, nisep )
            call write_sourced_value( summary%local%separatrix%velocity_phi%deuterium_tritium, -vtor )
            call write_sourced_value( summary%local%separatrix_average%n_i%deuterium_tritium, u )
            call write_sourced_value( summary%local%separatrix_average%velocity_phi%deuterium_tritium, -v )
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
          & 0.5_R8 * (ni(jxa,jsep,1) + ni(topix(jxa,jsep),topiy(jxa,jsep),1)) )
        u = separatrix_average( ni(:,:,1), tmpCv )
#if ( IMAS_MINOR_VERSION > 36 || IMAS_MAJOR_VERSION > 3 )
        call write_sourced_value( summary%local%separatrix_average%n_i_total, u )
#endif
        u = separatrix_average( zeff, tmpCv )
        call write_sourced_value( summary%local%separatrix%zeff, &
          & 0.5_R8 * (zeff(jxa,jsep) + zeff(topix(jxa,jsep),topiy(jxa,jsep))) )
#if ( IMAS_MINOR_VERSION > 36 || IMAS_MAJOR_VERSION > 3 )
        call write_sourced_value( summary%local%separatrix_average%zeff, u )
#endif

! Data at limiter tangency point
        if (geometryType.eq.GEOMETRY_LIMITER) then
          call write_sourced_value( summary%local%limiter%t_e, &
           & 0.25_R8 * ( te(0,jsep) + te(0,jsep+1) + &
           &             te(nx-1,jsep) + te(nx-1,jsep+1) )/ev )
          call write_sourced_value( summary%local%limiter%t_i_average, &
           & 0.25_R8 * ( ti(0,jsep) + ti(0,jsep+1) + &
           &             ti(nx-1,jsep) + ti(nx-1,jsep+1) )/ev )
          call write_sourced_value( summary%local%limiter%n_e, &
           & 0.25_R8 * ( ne(0,jsep) + ne(0,jsep+1) + &
           &             ne(nx-1,jsep) + ne(nx-1,jsep+1) ) )
          do is = 1, nspecies
            is1 = eb2spcr(is)
            if (nint(zamax(is1)).eq.0) is1 = is1 + 1
            is2 = is1 + nfluids(is) - 1
            nisep = 0.0_R8
            do i = is1, is2
              nisep = nisep + &
                &  0.25_R8 * ( na(0,jsep,i) + na(0,jsep+1,i) + &
                &              na(nx-1,jsep,i) + na(nx-1,jsep+1,i) )
            end do
            select case (is_codes(eb2spcr(is)))
            case ('H')
              call write_sourced_value( summary%local%limiter%n_i%hydrogen, nisep )
            case ('D')
              call write_sourced_value( summary%local%limiter%n_i%deuterium, nisep )
            case ('T')
              call write_sourced_value( summary%local%limiter%n_i%tritium, nisep )
#if ( IMAS_MAJOR_VERSION > 4 || ( IMAS_MAJOR_VERSION == 4 && IMAS_MINOR_VERSION > 0 ) )
            case ('DT')
              call write_sourced_value( summary%local%limiter%n_i%deuterium_tritium, nisep )
#endif
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
          call write_sourced_value( summary%local%limiter%n_i_total, &
           & 0.25_R8 * ( ni(0,jsep,1) + ni(0,jsep+1,1) + &
           &             ni(nx-1,jsep,1) + ni(nx-1,jsep+1,1) ) )
          call write_sourced_value( summary%local%limiter%zeff, &
           & 0.25_R8 * ( zeff(0,jsep) + zeff(0,jsep+1) + &
           &             zeff(nx-1,jsep) + zeff(nx-1,jsep+1) ) )
          call write_sourced_value( summary%local%limiter%flux_expansion, flux_expansion(itrg(1)) )
#if ( IMAS_MINOR_VERSION > 31 || IMAS_MAJOR_VERSION > 3 )
          u = max( power_flux_peak(itrg(1)), power_flux_peak(itrg(2)) )
          call write_sourced_value( summary%local%limiter%power_flux_peak, u )
#endif
        end if

! Summary divertor plate data
        if (ntrgts.gt.0) then
#if ( IMAS_MINOR_VERSION > 34 || IMAS_MAJOR_VERSION > 3 )
          allocate ( summary%local%divertor_target( ntrgts ) )
#else
          allocate ( summary%local%divertor_plate( ntrgts ) )
#endif
          do i = 1, ntrgts
#if ( IMAS_MINOR_VERSION > 34 || IMAS_MAJOR_VERSION > 3 )
            call write_sourced_string( summary%local%divertor_target(i)%name, plate_name(i) )
            call write_sourced_value( summary%local%divertor_target(i)%t_e, &
              &  0.5_R8 * (te(ixpos(itrg(i)),iypos(itrg(i)))+  &
              &            te(topix(ixpos(itrg(i)),iypos(itrg(i))), &
              &               topiy(ixpos(itrg(i)),iypos(itrg(i)))))/ev )
            call write_sourced_value( summary%local%divertor_target(i)%t_i_average, &
              &  0.5_R8 * (ti(ixpos(itrg(i)),iypos(itrg(i)))+  &
              &            ti(topix(ixpos(itrg(i)),iypos(itrg(i))), &
              &               topiy(ixpos(itrg(i)),iypos(itrg(i)))))/ev )
            call write_sourced_value( summary%local%divertor_target(i)%n_e, &
              &  0.5_R8 * (ne(ixpos(itrg(i)),iypos(itrg(i)))+  &
              &            ne(topix(ixpos(itrg(i)),iypos(itrg(i))), &
              &               topiy(ixpos(itrg(i)),iypos(itrg(i))))) )
#else
            call write_sourced_string( summary%local%divertor_plate(i)%name, plate_name(i) )
            call write_sourced_value( summary%local%divertor_plate(i)%t_e, &
              &  0.5_R8 * (te(ixpos(itrg(i)),iypos(itrg(i)))+  &
              &            te(topix(ixpos(itrg(i)),iypos(itrg(i))), &
              &               topiy(ixpos(itrg(i)),iypos(itrg(i)))))/ev )
            call write_sourced_value( summary%local%divertor_plate(i)%t_i_average, &
              &  0.5_R8 * (ti(ixpos(itrg(i)),iypos(itrg(i)))+  &
              &            ti(topix(ixpos(itrg(i)),iypos(itrg(i))), &
              &               topiy(ixpos(itrg(i)),iypos(itrg(i)))))/ev )
            call write_sourced_value( summary%local%divertor_plate(i)%n_e, &
              &  0.5_R8 * (ne(ixpos(itrg(i)),iypos(itrg(i)))+  &
              &            ne(topix(ixpos(itrg(i)),iypos(itrg(i))), &
              &               topiy(ixpos(itrg(i)),iypos(itrg(i))))) )
#endif
            do is = 1, nspecies
              is1 = eb2spcr(is)
              if (nint(zamax(is1)).eq.0) is1 = is1+1
              is2 = is1 + nfluids(is) - 1
              nisep = 0.0_R8
              do j = is1, is2
                nisep = nisep + &
                  &  0.5_R8 * (na(ixpos(itrg(i)),iypos(itrg(i)),j) + &
                  &            na(topix(ixpos(itrg(i)),iypos(itrg(i))), &
                  &               topiy(ixpos(itrg(i)),iypos(itrg(i))),j))
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
#if ( IMAS_MAJOR_VERSION > 4 || ( IMAS_MAJOR_VERSION == 4 && IMAS_MINOR_VERSION > 0 ) )
              case ('DT')
                call write_sourced_value( summary%local%divertor_target(i)%n_i%deuterium_tritium, nisep )
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
            call write_sourced_value( summary%local%divertor_target(i)%n_i_total, &
              & 0.5_R8 * (ni(ixpos(itrg(i)),iypos(itrg(i)),1) + &
              &           ni(topix(ixpos(itrg(i)),iypos(itrg(i))), &
              &              topiy(ixpos(itrg(i)),iypos(itrg(i))),1)) )
            call write_sourced_value( summary%local%divertor_target(i)%zeff, &
              & 0.5_R8 * (zeff(ixpos(itrg(i)),iypos(itrg(i))) + &
              &           zeff(topix(ixpos(itrg(i)),iypos(itrg(i))), &
              &                topiy(ixpos(itrg(i)),iypos(itrg(i))))) )
            call write_sourced_value( summary%local%divertor_target(i)%flux_expansion, &
              & flux_expansion(itrg(i)) )
            call write_sourced_value( summary%local%divertor_target(i)%power_flux_peak, &
              & power_flux_peak(itrg(i)) )
#else
            call write_sourced_value( summary%local%divertor_plate(i)%n_i_total, &
              & 0.5_R8 * (ni(ixpos(itrg(i)),iypos(itrg(i)),1) + &
              &           ni(topix(ixpos(itrg(i)),iypos(itrg(i))), &
              &              topiy(ixpos(itrg(i)),iypos(itrg(i))),1)) )
            call write_sourced_value( summary%local%divertor_plate(i)%zeff, &
              & 0.5_R8 * (zeff(ixpos(itrg(i)),iypos(itrg(i))) + &
              &           zeff(topix(ixpos(itrg(i)),iypos(itrg(i))), &
              &                topiy(ixpos(itrg(i)),iypos(itrg(i))))) )
            call write_sourced_value( summary%local%divertor_plate(i)%flux_expansion, &
              & flux_expansion(itrg(i)) )
#if ( IMAS_MINOR_VERSION > 31 || IMAS_MAJOR_VERSION > 3 )
            call write_sourced_value( summary%local%divertor_plate(i)%power_flux_peak, &
              & power_flux_peak(itrg(i)) )
#endif
#endif
          end do
        end if

        call fill_summary_data( summary )
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
            call write_sourced_value( summary%scrape_off_layer%t_e_decay_length, enepar(ib,1) )
          end if
          if (bceni(ib).eq.9.or.bceni(ib).eq.19) then
            call write_sourced_value( summary%scrape_off_layer%t_i_average_decay_length, &
               & enipar(ib,1) )
          end if
          nibnd = IDS_REAL_INVALID
          match_found = .true.
          do is = 0, ns-1
            if (is_neutral(is).and.use_eirene.ne.0) cycle
            if (bccon(is,ib).eq.9) then
              if (nibnd.eq.IDS_REAL_INVALID) then
                nibnd = conpar(is,ib,1)
              else
                match_found = match_found.and.nibnd.eq.conpar(is,ib,1)
              end if
            end if
          end do
          if (match_found.and.nibnd.ne.IDS_REAL_INVALID) then
            call write_sourced_value( summary%scrape_off_layer%n_e_decay_length, nibnd )
            call write_sourced_value( summary%scrape_off_layer%n_i_total_decay_length, nibnd )
          end if
        end if
        u = 0.0_IDS_real
        do ix = 0, nx-1
          do iy = 0, ny-1
            if (geometryType.eq.GEOMETRY_LIMITER .or. &
             &  geometryType.eq.GEOMETRY_SN .or. &
             &  geometryType.eq.GEOMETRY_LFS_SNOWFLAKE_MINUS .or. &
             &  geometryType.eq.GEOMETRY_LFS_SNOWFLAKE_PLUS) then
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
          call write_sourced_value( summary%scrape_off_layer%power_radiated, u )
        end if
        select case (geometryType)
        case (GEOMETRY_LINEAR, GEOMETRY_CYLINDER)
          if (ntrgts.ge.1) then
            ix = jxa
          else
            ix = -2
          end if
        case (GEOMETRY_LIMITER)
          if (LSN) then
            ix = nx-1
          else
            ix = 0
          end if
        case (GEOMETRY_SN, GEOMETRY_CDN, GEOMETRY_DDN_BOTTOM, &
        &     GEOMETRY_LFS_SNOWFLAKE_MINUS, GEOMETRY_LFS_SNOWFLAKE_PLUS)
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
          qemax = 0.0_IDS_real
          qimax = 0.0_IDS_real
          do iy = jsep+1, ny-1
            qemax = qemax + (-icnt)*fhe(ix,iy,0,0)
            qimax = qimax + (-icnt)*fhi(ix,iy,0,0)
          end do
          do i = ix+icnt, jxa, icnt
            qetot = 0.0_IDS_real
            qitot = 0.0_IDS_real
            do iy = jsep+1, ny-1
              qetot = qetot + (-icnt)*fhe(i,iy,0,0)
              qitot = qitot + (-icnt)*fhi(i,iy,0,0)
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
          totFace(:,:,0)=abs(fhe(:,:,0,0))
          totFace(:,:,1)=abs(fhe(:,:,1,1))
          call divide_by_poloidal_areas(nx,ny,totFace,tmpFace)
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
          call write_sourced_value( summary%scrape_off_layer%heat_flux_e_decay_length, lambda )
          lambda = 0.0_IDS_real
          totFace(:,:,0)=abs(fhi(:,:,0,0))
          totFace(:,:,1)=abs(fhi(:,:,1,1))
          call divide_by_poloidal_areas(nx,ny,totFace,tmpFace)
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
          call write_sourced_value( summary%scrape_off_layer%heat_flux_i_decay_length, lambda )
        end if
#endif

        deallocate(ionstt,istion,ispion)
#ifdef B25_EIRENE
        if (use_eirene.ne.0) then
          deallocate(isstat,imneut,imiion)
          deallocate(in_species)
          if (allocated(un0)) deallocate(un0,um0)
        end if
#endif
        call logmsg( LOGDEBUG, "b2mod_ual_io.B25_process_ids: done" )
        ncall = ncall + 1

        contains

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
            &   batch_profiles, batch_sources, &
#if IMAS_MAJOR_VERSION > 3
            &   batch_plasma_profiles, batch_plasma_sources, &
#endif
#if ( IMAS_MAJOR_VERSION < 4 || ( IMAS_MAJOR_VERSION == 4 && IMAS_MINOR_VERSION < 1 ) )
            &   description, &
#endif
            &   equilibrium, &
#if ( IMAS_MINOR_VERSION > 21 || IMAS_MAJOR_VERSION > 3 )
            &   summary, &
#if IMAS_MAJOR_VERSION == 3
            &   time_IN, shot, run, database, version, &
#else
            &   time_IN, shot, database, &
#endif
#endif
            &   new_eq_ggd, &
            &   batch_ind_IN, num_batch_slices_IN )
#ifdef NO_OPT
!DIR$ NOOPTIMIZE
#endif
        implicit none
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
#if IMAS_MAJOR_VERSION > 3
        type (ids_plasma_profiles) :: batch_plasma_profiles !< IDS designed to
            !< store data on plasma profiles
        type (ids_plasma_sources) :: batch_plasma_sources !< IDS designed to
            !< store data on plasma sources. Energy terms correspond to the full
            !< kinetic energy equation (i.e. the energy flux takes into account
            !< the energy transported by the particle flux)
#endif
#if ( IMAS_MAJOR_VERSION < 4 || ( IMAS_MAJOR_VERSION == 4 && IMAS_MINOR_VERSION < 1 ) )
        type (ids_dataset_description) :: description !< IDS designed to store
            !< a description of the simulation
#endif
#if ( IMAS_MINOR_VERSION > 21 || IMAS_MAJOR_VERSION > 3 )
        type (ids_summary) :: summary !< IDS designed to store
            !< run summary data
#if IMAS_MAJOR_VERSION == 3
        integer, intent(in) :: run
        character(len=24), intent(in) :: version
#endif
#endif
        integer, intent(in) :: shot
        character(len=24), intent(in) :: database
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
        integer :: ns    !< Total number of B2.5 species
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
#if ( IMAS_MAJOR_VERSION > 4 || ( IMAS_MAJOR_VERSION == 4 && IMAS_MINOR_VERSION > 0 ) )
        type(ids_plasma_profiles_time_slice) :: profiles_ggd
        type(ids_plasma_sources_source_ggd), allocatable :: sources_ggd(:)
#endif
        real(IDS_real) :: tmpCv( -1:ubound( na, 1), -1:ubound( na, 2) )
        real(IDS_real) :: totCv( -1:ubound( na, 1), -1:ubound( na, 2) )
#endif
#if ( IMAS_MAJOR_VERSION > 4 || ( IMAS_MAJOR_VERSION == 4 && IMAS_MINOR_VERSION > 0 ) )
        integer :: ix, iy
        real(IDS_real) :: u, frac
#endif
        !! Procedures
        external species, xertst

        !! ===  SET UP IDS ===
        write(0,*) "Setting data for batch_profiles IDS"
        source = "SOLPS-ITER (batch averaged)"

        !! Preparing database for writing
        call IDS_init
        ns = size( na, 3 )
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
#if IMAS_MAJOR_VERSION > 3
        call write_ids_properties( batch_plasma_profiles%ids_properties, &
          &  homogeneous_time )
        call write_ids_properties( batch_plasma_sources%ids_properties, &
          &  homogeneous_time )
#endif
        if ( do_description ) then
#if ( IMAS_MAJOR_VERSION < 4 || ( IMAS_MAJOR_VERSION == 4 && IMAS_MINOR_VERSION < 1 ) )
          call write_ids_properties( description%ids_properties, 2 )
#endif
#if ( IMAS_MINOR_VERSION > 21 || IMAS_MAJOR_VERSION > 3 )
          call write_ids_properties( summary%ids_properties, &
            &  homogeneous_time )
#endif
        end if

        !! 2. Set code and library data
        code_description = "Batch-averaged IDS from b2mod_ual_io routine"
        call write_ids_code( batch_profiles%code, code_commit, code_description )
        call write_ids_code( batch_sources%code, code_commit, code_description )
#if IMAS_MAJOR_VERSION > 3
        call write_ids_code( batch_plasma_profiles%code, code_commit, code_description )
        call write_ids_code( batch_plasma_sources%code, code_commit, code_description )
#endif
#if ( IMAS_MINOR_VERSION > 21 || IMAS_MAJOR_VERSION > 3 )
        if (do_description) &
          &  call write_ids_code( summary%code, code_commit, code_description )
#endif
#if ( IMAS_MAJOR_VERSION == 4 && IMAS_MINOR_VERSION == 0 )
        if (do_description) &
          &  call write_ids_code_constant( description%code, code_commit, code_description )
#endif

        !! 3. Allocate IDS.time and set it to desired values
        allocate( batch_profiles%time(num_batch_slices) )
        batch_profiles%time(batch_index) = time
        allocate( batch_sources%time(num_batch_slices) )
        batch_sources%time(batch_index) = time
#if IMAS_MAJOR_VERSION > 3
        allocate( batch_plasma_profiles%time(num_batch_slices) )
        batch_plasma_profiles%time(batch_index) = time
        allocate( batch_plasma_sources%time(num_batch_slices) )
        batch_plasma_sources%time(batch_index) = time
#endif
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
        allocate( batch_sources%source(1) )
        allocate( batch_sources%source(1)%ggd( num_batch_slices ) )
#if IMAS_MAJOR_VERSION > 3
        allocate( batch_plasma_profiles%ggd( num_batch_slices ) )
        allocate( batch_plasma_profiles%grid_ggd( num_batch_slices ) )
        allocate( batch_plasma_sources%grid_ggd( num_batch_slices ) )
        allocate( batch_plasma_sources%source(1) )
        allocate( batch_plasma_sources%source(1)%ggd( num_batch_slices ) )
#if ( IMAS_MAJOR_VERSION > 4 || ( IMAS_MAJOR_VERSION == 4 && IMAS_MINOR_VERSION > 0 ) )
        allocate( sources_ggd(1) )
#endif
#endif
#ifdef B25_EIRENE
#if ( IMAS_MAJOR_VERSION > 4 || ( IMAS_MAJOR_VERSION == 4 && IMAS_MINOR_VERSION > 0 ) )
        !! Neutrals
        call set_edge_source_identifier( batch_sources%source(1)%identifier, "Neutrals" )
#elif ( IMAS_MINOR_VERSION > 38 || IMAS_MAJOR_VERSION > 3 )
        !! Neutrals
        call write_source_identifier( &
          &  batch_sources%source(1)%identifier, edge_source_identifier%neutrals )
#if IMAS_MAJOR_VERSION > 3
        call write_source_identifier( &
          &  batch_plasma_sources%source(1)%identifier, edge_source_identifier%neutrals )
#endif
#else
        !! Total sources due to Eirene species
        allocate( batch_sources%source(1)%identifier%name(1) )
        allocate( batch_sources%source(1)%identifier%description(1) )
        batch_sources%source(1)%identifier%index = 0
        batch_sources%source(1)%identifier%name = "Eirene"
        batch_sources%source(1)%identifier%description = &
             &  "Total source due to Eirene species "//trim(source)
#endif
#endif

        call put_equilibrium_data ( equilibrium, &
#if ( IMAS_MINOR_VERSION > 21 || IMAS_MAJOR_VERSION > 3 )
            &  summary, &
#endif
            &  batch_profiles, database, &
#if ( AL_MAJOR_VERSION > 4 && GGD_MAJOR_VERSION > 0 )
            &  batch_index, &
#endif
            &  do_description, new_eq_ggd )
#if IMAS_MAJOR_VERSION > 3
        allocate( batch_plasma_profiles%vacuum_toroidal_field%b0( num_time_slices ) )
        batch_plasma_profiles%vacuum_toroidal_field%b0( batch_index ) = &
            &  batch_profiles%vacuum_toroidal_field%b0( batch_index )
        batch_plasma_profiles%vacuum_toroidal_field%r0 = &
            &  batch_profiles%vacuum_toroidal_field%r0
#endif
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
#elif ( IMAS_MAJOR_VERSION == 4 && IMAS_MINOR_VERSION == 0 )
          description%type%index = 2
          allocate( description%type%name(1) )
          allocate( description%type%description(1) )
          description%type%name = "simulation"
          description%type%description = "Simulation results from "//trim(source)
          allocate( description%machine(1) )
          description%machine = database
          description%pulse = shot
#elif ( IMAS_MAJOR_VERSION > 4 || ( IMAS_MAJOR_VERSION == 4 && IMAS_MINOR_VERSION > 0 ) )
          summary%type%index = 2
          allocate( summary%type%name(1) )
          allocate( summary%type%description(1) )
          summary%type%name = "simulation"
          summary%type%description = "Simulation results from "//trim(source)
          allocate( summary%machine(1) )
          summary%machine = database
          summary%pulse = shot
#endif
#if ( IMAS_MAJOR_VERSION < 4 || ( IMAS_MAJOR_VERSION == 4 && IMAS_MINOR_VERSION < 1 ) )
          if ( present( time_IN ) ) &
            &  description%simulation%time_current = time_IN
          allocate( description%simulation%workflow(1) )
          description%simulation%workflow = source
#else
          allocate( summary%simulation%workflow(1) )
          summary%simulation%workflow = source
#endif

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
#if IMAS_MAJOR_VERSION > 3
        call write_ids_midplane( batch_plasma_profiles%midplane, midplane_id )
        call write_ids_midplane( batch_plasma_sources%midplane, midplane_id )
#endif
        if (do_description) &
          & call write_ids_midplane( summary%midplane, midplane_id )
#endif

        !! Write grid & grid subsets/subgrids
#if ( IMAS_MINOR_VERSION > 11 || IMAS_MAJOR_VERSION > 3 ) && GGD_MAJOR_VERSION > 0
#if ( IMAS_MINOR_VERSION < 15 && IMAS_MAJOR_VERSION < 4 )
        call b2_IMAS_Fill_Grid_Desc( IDSmap,                                &
            &   batch_profiles%ggd( batch_index )%grid,                     &
            &   nx, ny, crx(-1:nx, -1:ny, :), cry(-1:nx, -1:ny, : ),        &
            &   leftix, leftiy, rightix, rightiy, topix, topiy, bottomix,   &
            &   bottomiy, nnreg, topcut, region, cflags,                    &
            &   INCLUDE_GHOST_CELLS, midplane_id, z_eq, vol, gs, qc )
        call b2_IMAS_Fill_Grid_Desc( IDSmap,                                &
            &   batch_sources%source(1)%ggd( batch_index )%grid,            &
            &   nx, ny, crx(-1:nx, -1:ny, :), cry(-1:nx, -1:ny, : ),        &
            &   leftix, leftiy, rightix, rightiy, topix, topiy, bottomix,   &
            &   bottomiy, nnreg, topcut, region, cflags,                    &
            &   INCLUDE_GHOST_CELLS, midplane_id, z_eq, vol, gs, qc )
#else
        call b2_IMAS_Fill_Grid_Desc( IDSmap,                                &
            &   batch_profiles%grid_ggd( batch_index ),                     &
            &   nx, ny, crx(-1:nx, -1:ny, :), cry(-1:nx, -1:ny, : ),        &
            &   leftix, leftiy, rightix, rightiy, topix, topiy, bottomix,   &
            &   bottomiy, nnreg, topcut, region, cflags,                    &
            &   INCLUDE_GHOST_CELLS, midplane_id, z_eq, vol, gs, qc )
#if AL_MAJOR_VERSION > 4
        allocate( batch_sources%grid_ggd( batch_index )%path(1) )
        batch_sources%grid_ggd( batch_index )%path = &
            &   "#edge_profiles(1)/grid_ggd("//int2str(batch_index)//")"
#if IMAS_MAJOR_VERSION > 3
        allocate( batch_plasma_profiles%grid_ggd( batch_index )%path(1) )
        batch_plasma_profiles%grid_ggd( batch_index )%path = &
            &   "#edge_profiles(1)/grid_ggd("//int2str(batch_index)//")"
        allocate( batch_plasma_sources%grid_ggd( batch_index )%path(1) )
        batch_plasma_sources%grid_ggd( batch_index )%path = &
            &   "#edge_profiles(1)/grid_ggd("//int2str(batch_index)//")"
#endif
#else
        call b2_IMAS_Fill_Grid_Desc( IDSmap,                                &
            &   batch_sources%grid_ggd( batch_index ),                      &
            &   nx, ny, crx(-1:nx, -1:ny, :), cry(-1:nx, -1:ny, : ),        &
            &   leftix, leftiy, rightix, rightiy, topix, topiy, bottomix,   &
            &   bottomiy, nnreg, topcut, region, cflags,                    &
            &   INCLUDE_GHOST_CELLS, midplane_id, z_eq, vol, gs, qc )
#endif
#endif
#else
        if (do_description) then
          write(0,*) 'Code was compiled without a GGD module'
          write(0,*) 'Most IDS output is disabled !'
        end if
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
#if ( IMAS_MAJOR_VERSION < 4 || ( IMAS_MAJOR_VERSION == 4 && IMAS_MINOR_VERSION < 1 ) )
        batch_profiles%ggd( batch_index )%time = batch_slice_value
        batch_sources%source(1)%ggd( batch_index )%time = batch_slice_value
#else
        profiles_ggd%time = batch_slice_value
        sources_ggd(1)%time = batch_slice_value
#endif
#if IMAS_MAJOR_VERSION > 3
        batch_plasma_profiles%grid_ggd( batch_index )%time = batch_slice_value
#if ( IMAS_MAJOR_VERSION == 4 && IMAS_MINOR_VERSION < 1 )
        batch_plasma_profiles%ggd( batch_index )%time = batch_slice_value
        batch_plasma_sources%grid_ggd( batch_index )%time = batch_slice_value
        batch_plasma_sources%source(1)%ggd( batch_index )%time = batch_slice_value
#endif
#endif

        !! List of species
        !! Careful here: ion in DD means isonuclear sequence !!
#if ( IMAS_MAJOR_VERSION < 4 || ( IMAS_MAJOR_VERSION == 4 && IMAS_MINOR_VERSION < 1 ) )
        allocate( batch_profiles%ggd( batch_index )%ion( nspecies ) )
        allocate( batch_profiles%ggd( batch_index )%neutral( nspecies ) )
        allocate( batch_sources%source(1)%ggd( batch_index )%ion( nspecies ) )
        allocate( batch_sources%source(1)%ggd( batch_index )%neutral( nspecies ) )
#else
        allocate( profiles_ggd%ion( nspecies ) )
        allocate( profiles_ggd%neutral( nspecies ) )
        allocate( sources_ggd(1)%ion( nspecies ) )
        allocate( sources_ggd(1)%neutral( nspecies ) )
#endif
        ks = 0
        do js = 1, nspecies
          if (is_neutral(ks)) ks = ks + 1  ! Skip the neutral species
#if ( IMAS_MAJOR_VERSION < 4 || ( IMAS_MAJOR_VERSION == 4 && IMAS_MINOR_VERSION < 1 ) )
          allocate( batch_profiles%ggd( batch_index )%ion( js )%state( nfluids(js) ) )
          allocate( batch_profiles%ggd( batch_index )%ion( js )%element(1) )
          allocate( batch_sources%source(1)%ggd( batch_index )%ion( js )%state( nfluids(js) ) )
          allocate( batch_sources%source(1)%ggd( batch_index )%ion( js )%element(1) )
#else
          allocate( profiles_ggd%ion( js )%state( nfluids(js) ) )
          allocate( profiles_ggd%ion( js )%element(1) )
          allocate( sources_ggd(1)%ion( js )%state( nfluids(js) ) )
          allocate( sources_ggd(1)%ion( js )%element(1) )
#endif
#if ( IMAS_MAJOR_VERSION < 4 && IMAS_MINOR_VERSION < 42 )
          allocate( batch_profiles%ggd( batch_index )%ion( js )%label(1) )
          allocate( batch_sources%source(1)%ggd( batch_index )%ion( js )%label(1) )
          do is = 1, nfluids(js)
            allocate( batch_profiles%ggd( batch_index )%ion( js )%state( is )%label(1) )
            allocate( batch_sources%source(1)%ggd( batch_index )%ion( js )%state( is )%label(1) )
          end do
#elif ( IMAS_MAJOR_VERSION < 4 || ( IMAS_MAJOR_VERSION == 4 && IMAS_MINOR_VERSION < 1 ) )
          allocate( batch_profiles%ggd( batch_index )%ion( js )%name(1) )
          allocate( batch_sources%source(1)%ggd( batch_index )%ion( js )%name(1) )
          do is = 1, nfluids(js)
            allocate( batch_profiles%ggd( batch_index )%ion( js )%state( is )%name(1) )
            allocate( batch_sources%source(1)%ggd( batch_index )%ion( js )%state( is )%name(1) )
          end do
#else
          allocate( profiles_ggd%ion( js )%name(1) )
          allocate( sources_ggd(1)%ion( js )%name(1) )
          do is = 1, nfluids(js)
            allocate( profiles_ggd%ion( js )%state( is )%name(1) )
            allocate( sources_ggd(1)%ion( js )%state( is )%name(1) )
          end do
#endif
          do is = 1, nfluids(js)
            call species( ks, spclabel, .false.)
            call shrink_label(spclabel)
#if ( IMAS_MAJOR_VERSION < 4 && IMAS_MINOR_VERSION < 42 )
            batch_profiles%ggd( batch_index )%ion( js )%state( is )%label = spclabel
            batch_sources%source(1)%ggd( batch_index )%ion( js )%state( is )%label = spclabel
#elif ( IMAS_MAJOR_VERSION < 4 || ( IMAS_MAJOR_VERSION == 4 && IMAS_MINOR_VERSION < 1 ) )
            batch_profiles%ggd( batch_index )%ion( js )%state( is )%name = spclabel
            batch_sources%source(1)%ggd( batch_index )%ion( js )%state( is )%name = spclabel
#else
            profiles_ggd%ion( js )%state( is )%name = spclabel
            sources_ggd(1)%ion( js )%state( is )%name = spclabel
#endif
            ! Put minimum Z of the charge state bundle
            ! (z_min = z_max = 0 for a neutral)
#if ( IMAS_MAJOR_VERSION < 4 || ( IMAS_MAJOR_VERSION == 4 && IMAS_MINOR_VERSION < 1 ) )
            batch_profiles%ggd( batch_index )%ion( js )%state( is )%z_min = zamin( ks )
            batch_sources%source(1)%ggd( batch_index )%ion( js )%state( is )%z_min = &
                  &  zamin( ks )
#else
            profiles_ggd%ion( js )%state( is )%z_min = zamin( ks )
            sources_ggd(1)%ion( js )%state( is )%z_min = zamin( ks )
#endif
            ! Put maximum Z of the charge state bundle
#if ( IMAS_MAJOR_VERSION < 4 || ( IMAS_MAJOR_VERSION == 4 && IMAS_MINOR_VERSION < 1 ) )
            batch_profiles%ggd( batch_index )%ion( js )%state( is )%z_max = zamax( ks )
            batch_sources%source(1)%ggd( batch_index )%ion( js )%state( is )%z_max = &
                  &  zamax( ks )
#else
            profiles_ggd%ion( js )%state( is )%z_max = zamax( ks )
            sources_ggd(1)%ion( js )%state( is )%z_max = zamax( ks )
#endif
            if (is.lt.nfluids(js)) ks = ks + 1
          end do

          ! Put ion label identifying the species
#if ( IMAS_MAJOR_VERSION < 4 && IMAS_MINOR_VERSION < 42 )
          batch_profiles%ggd( batch_index )%ion( js )%label = species_list( js )
          batch_sources%source(1)%ggd( batch_index )%ion( js )%label = species_list( js )
#elif ( IMAS_MAJOR_VERSION < 4 || ( IMAS_MAJOR_VERSION == 4 && IMAS_MINOR_VERSION < 1 ) )
          batch_profiles%ggd( batch_index )%ion( js )%name = species_list( js )
          batch_sources%source(1)%ggd( batch_index )%ion( js )%name = species_list( js )
#else
          profiles_ggd%ion( js )%name = species_list( js )
          sources_ggd(1)%ion( js )%name = species_list( js )
#endif
          ! Put ion charge if single ion in species
          if (nfluids(js).eq.1) then
            ion_charge_int = nint((zamin(ks)+zamax(ks))/2.0_R8)
#if ( IMAS_MAJOR_VERSION < 4 || ( IMAS_MAJOR_VERSION == 4 && IMAS_MINOR_VERSION < 1 ) )
            batch_profiles%ggd( batch_index )%ion( js )%z_ion = ion_charge_int
            batch_sources%source(1)%ggd( batch_index )%ion( js )%z_ion = ion_charge_int
#else
            profiles_ggd%ion( js )%z_ion = ion_charge_int
            sources_ggd(1)%ion( js )%z_ion = ion_charge_int
#endif
          end if

          ! Put element data
#if ( IMAS_MAJOR_VERSION < 4 || ( IMAS_MAJOR_VERSION == 4 && IMAS_MINOR_VERSION < 1 ) )
          call fill_neutral_element( ks, js, &
            &  batch_profiles%ggd( batch_index )%ion( js )%element(1) )
          call fill_neutral_element( ks, js, &
            &  batch_sources%source(1)%ggd( batch_index )%ion( js )%element(1) )
#else
          call fill_neutral_element( ks, js, &
            &  profiles_ggd%ion( js )%element(1) )
          call fill_neutral_element( ks, js, &
            &  sources_ggd(1)%ion( js )%element(1) )
#endif

          ! Put neutral index
#if ( IMAS_MAJOR_VERSION < 4 || ( IMAS_MAJOR_VERSION == 4 && IMAS_MINOR_VERSION < 1 ) )
          batch_profiles%ggd( batch_index )%ion( js )%neutral_index = js
          batch_sources%source(1)%ggd( batch_index )%ion( js )%neutral_index = js
#else
          profiles_ggd%ion( js )%neutral_index = js
          sources_ggd(1)%ion( js )%neutral_index = js
#endif
          ! Put multiple states flag
#if ( IMAS_MAJOR_VERSION < 4 || ( IMAS_MAJOR_VERSION == 4 && IMAS_MINOR_VERSION < 1 ) )
          batch_profiles%ggd( batch_index )%ion( js )%multiple_states_flag = 1
          batch_sources%source(1)%ggd( batch_index )%ion( js )%multiple_states_flag = 1
#else
          profiles_ggd%ion( js )%multiple_states_flag = 1
          sources_ggd(1)%ion( js )%multiple_states_flag = 1
#endif

          !! List of neutrals
#if ( IMAS_MAJOR_VERSION < 4 || ( IMAS_MAJOR_VERSION == 4 && IMAS_MINOR_VERSION < 1 ) )
          allocate( batch_profiles%ggd( batch_index )%neutral( js )%element(1) )
          allocate( batch_sources%source(1)%ggd( batch_index )%neutral( js )%element(1) )
          call fill_neutral_element( ks, js, &
            &  batch_profiles%ggd( batch_index )%neutral( js )%element(1) )
          call fill_neutral_element( ks, js, &
            &  batch_sources%source(1)%ggd( batch_index )%neutral( js )%element(1) )
          batch_profiles%ggd( batch_index )%neutral( js )%ion_index = js
          batch_profiles%ggd( batch_index )%neutral( js )%multiple_states_flag = 0
          batch_sources%source(1)%ggd( batch_index )%neutral( js )%ion_index = js
          batch_sources%source(1)%ggd( batch_index )%neutral( js )%multiple_states_flag = 0
#else
          allocate( profiles_ggd%neutral( js )%element(1) )
          call fill_neutral_element( ks, js, &
            &  profiles_ggd%neutral( js )%element(1) )
          profiles_ggd%neutral( js )%ion_index = js
          profiles_ggd%neutral( js )%multiple_states_flag = 0
          allocate( sources_ggd(1)%neutral( js )%element(1) )
          call fill_neutral_element( ks, js, &
            &  sources_ggd(1)%neutral( js )%element(1) )
          sources_ggd(1)%neutral( js )%ion_index = js
          sources_ggd(1)%neutral( js )%multiple_states_flag = 0
#endif
#if ( IMAS_MAJOR_VERSION < 4 && IMAS_MINOR_VERSION < 42 )
          allocate( batch_profiles%ggd( batch_index )%neutral( js )%label(1) )
          allocate( batch_sources%source(1)%ggd( batch_index )%neutral( js )%label(1) )
          batch_profiles%ggd( batch_index )%neutral( js )%label = species_list( js )
          batch_sources%source(1)%ggd( batch_index )%neutral( js )%label = species_list( js )
#elif ( IMAS_MAJOR_VERSION < 4 || ( IMAS_MAJOR_VERSION == 4 && IMAS_MINOR_VERSION < 1 ) )
          allocate( batch_profiles%ggd( batch_index )%neutral( js )%name(1) )
          batch_profiles%ggd( batch_index )%neutral( js )%name = species_list( js )
          allocate( batch_sources%source(1)%ggd( batch_index )%neutral( js )%name(1) )
          batch_sources%source(1)%ggd( batch_index )%neutral( js )%name = species_list( js )
#else
          allocate( profiles_ggd%neutral( js )%name(1) )
          profiles_ggd%neutral( js )%name = species_list( js )
          allocate( sources_ggd(1)%neutral( js )%name(1) )
          sources_ggd(1)%neutral( js )%name = species_list( js )
#endif

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
          if (LSN) then
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
          if (LSN) then
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
          !! sne: Electron particle sources
          tmpCv(:,:) = sne_mean(:,:) / vol(:,:)
#if ( IMAS_MAJOR_VERSION < 4 || ( IMAS_MAJOR_VERSION == 4 && IMAS_MINOR_VERSION < 1 ) )
          call write_cell_scalar( sources_grid,                         &
              &   scalar = batch_sources%source(1)%ggd( batch_index )%  &
              &            electrons%particles,                         &
              &   b2CellData = tmpCv )
#else
          call write_cell_scalar( sources_grid,                         &
              &   scalar = sources_ggd(1)%electrons%particles,          &
              &   b2CellData = tmpCv )
#endif

#if ( IMAS_MAJOR_VERSION < 4 || ( IMAS_MAJOR_VERSION == 4 && IMAS_MINOR_VERSION < 1 ) )
          !! na: Ion density
          ks = 0
          do is = 1, nspecies
            if (is_neutral(ks)) ks = ks + 1
            tmpCv(:,:) = 0.0_IDS_real
            do js = 1, nfluids(is)
              call write_quantity( batch_grid,                          &
                  &   val = batch_profiles%ggd( batch_index )%          &
                  &         ion( is )%state( js )%density,              &
                  &   value = na_mean(:,:,ks) )
              tmpCv(:,:) = tmpCv(:,:) + na_mean(:,:,ks)
              ks = ks + 1
            end do
            call write_quantity( batch_grid,                            &
                &   val = batch_profiles%ggd( batch_index )%            &
                &         ion(is)%density,                              &
                &   value = tmpCv )
          end do
          !! sna: Ion particle sources
          ks = 0
          do is = 1, nspecies
            if (is_neutral(ks)) ks = ks + 1
            totCv(:,:) = 0.0_IDS_real
            do js = 1, nfluids(is)
              tmpCv(:,:) = sna_mean(:,:,ks) / vol(:,:)
              totCv(:,:) = totCv(:,:) + tmpCv(:,:)
              call write_cell_scalar( sources_grid,                       &
                  &   scalar = batch_sources%source(1)%                   &
                  &   ggd( batch_index )%ion( is )%state( js )%particles, &
                  &   b2CellData = tmpCv )
              ks = ks + 1
            end do
            call write_cell_scalar( sources_grid,                         &
                &   scalar = batch_sources%source(1)%                     &
                &   ggd( batch_index )%ion( is )%particles,               &
                &   b2CellData = totCv )
          end do
          !! ua: Parallel ion velocity
          ks = 0
          do is = 1, nspecies
            if (is_neutral(ks)) ks = ks + 1
            do js = 1, nfluids(is)
              call write_cell_vector_component( batch_grid,               &
                  &   vectorComponent = batch_profiles%                   &
                  &   ggd( batch_index )%ion( is )%state( js )%velocity,  &
                  &   b2CellData = ua_mean(:,:,ks),                       &
                  &   vectorID = VEC_ALIGN_PARALLEL_ID )
              ks = ks + 1
            end do
          end do
          !! smo: Ion parallel momentum sources
          ks = 0
          do is = 1, nspecies
            if (is_neutral(ks)) ks = ks + 1
            totCv(:,:) = 0.0_IDS_real
            do js = 1, nfluids(is)
              tmpCv(:,:) = smo_mean(:,:,ks) / vol(:,:)
              totCv(:,:) = totCv(:,:) + tmpCv(:,:)
              call write_cell_vector_component( sources_grid,             &
                  &   vectorComponent = batch_sources%source(1)%          &
                  &                     ggd( batch_index )%ion( is )%     &
                  &                     state( js )%momentum,             &
                  &   b2CellData = tmpCv,                                 &
                  &   vectorID = VEC_ALIGN_PARALLEL_ID )
              ks = ks + 1
            end do
            call write_cell_vector_component( sources_grid,               &
                &   vectorComponent = batch_sources%source(1)%            &
                &                     ggd( batch_index )%ion( is )%       &
                &                     momentum,                           &
                &   b2CellData = totCv,                                   &
                &   vectorID = VEC_ALIGN_PARALLEL_ID )
          end do
#else
          !! na: Ion density
          ks = 0
          do is = 1, nspecies
            if (is_neutral(ks)) ks = ks + 1
            tmpCv(:,:) = 0.0_IDS_real
            do js = 1, nfluids(is)
              call write_quantity( batch_grid,                          &
                  &   val = profiles_ggd%ion( is )%state( js )%density, &
                  &   value = na_mean(:,:,ks) )
              tmpCv(:,:) = tmpCv(:,:) + na_mean(:,:,ks)
              ks = ks + 1
            end do
            call write_quantity( batch_grid,                            &
                &   val = profiles_ggd%ion(is)%density,                 &
                &   value = tmpCv )
          end do
          !! sna: Ion particle sources
          ks = 0
          do is = 1, nspecies
            if (is_neutral(ks)) ks = ks + 1
            totCv(:,:) = 0.0_IDS_real
            do js = 1, nfluids(is)
              tmpCv(:,:) = sna_mean(:,:,ks) / vol(:,:)
              totCv(:,:) = totCv(:,:) + tmpCv(:,:)
              call write_cell_scalar( sources_grid,                       &
                  &   scalar = sources_ggd(1)%ion( is )%state( js )%      &
                  &            particles,                                 &
                  &   b2CellData = tmpCv )
              ks = ks + 1
            end do
            call write_cell_scalar( sources_grid,                         &
                &   scalar = sources_ggd(1)%ion( is )%particles,          &
                &   b2CellData = totCv )
          end do
          !! ua: Parallel ion velocity
          ks = 0
          do is = 1, nspecies
            if (is_neutral(ks)) ks = ks + 1
            do js = 1, nfluids(is)
              call write_cell_vector_component( batch_grid,               &
                  &   vectorComponent = profiles_ggd%                     &
                  &   ion( is )%state( js )%velocity,                     &
                  &   b2CellData = ua_mean(:,:,ks),                       &
                  &   vectorID = VEC_ALIGN_PARALLEL_ID )
              ks = ks + 1
            end do
          end do
          !! smo: Ion parallel momentum sources
          ks = 0
          do is = 1, nspecies
            if (is_neutral(ks)) ks = ks + 1
            totCv(:,:) = 0.0_IDS_real
            do js = 1, nfluids(is)
              tmpCv(:,:) = smo_mean(:,:,ks) / vol(:,:)
              totCv(:,:) = totCv(:,:) + tmpCv(:,:)
              call write_cell_vector_component( sources_grid,             &
                  &   vectorComponent = sources_ggd(1)%ion( is )%         &
                  &                     state( js )%momentum,             &
                  &   b2CellData = tmpCv,                                 &
                  &   vectorID = VEC_ALIGN_PARALLEL_ID )
              ks = ks + 1
            end do
            call write_cell_vector_component( sources_grid,               &
                &   vectorComponent = sources_ggd(1)%ion( is )%momentum,  &
                &   b2CellData = totCv,                                   &
                &   vectorID = VEC_ALIGN_PARALLEL_ID )
          end do
#endif
#if ( IMAS_MAJOR_VERSION < 4 || ( IMAS_MAJOR_VERSION == 4 && IMAS_MINOR_VERSION < 1 ) )
          !! te: Electron Temperature
          tmpCv(:,:) = te_mean(:,:)/qe
          call write_quantity( batch_grid,                               &
              &   val = batch_profiles%ggd( batch_index )%electrons%     &
              &         temperature,                                     &
              &   value = tmpCv )
          tmpCv(:,:) = she_mean(:,:) / vol(:,:)
          call write_cell_scalar( sources_grid,                          &
              &   scalar = batch_sources%source(1)%ggd( batch_index )%   &
              &            electrons%energy,                             &
              &   b2CellData = tmpCv )
          !! ti: (Common) Ion Temperature
          tmpCv(:,:) = ti_mean(:,:)/qe
          call write_quantity( batch_grid,                               &
              &   val = batch_profiles%ggd( batch_index )%t_i_average,   &
              &   value = tmpCv )
#else
          !! te: Electron Temperature
          tmpCv(:,:) = te_mean(:,:)/qe
          call write_quantity( batch_grid,                               &
              &   val = profiles_ggd%electrons%temperature,              &
              &   value = tmpCv )
          tmpCv(:,:) = she_mean(:,:) / vol(:,:)
          call write_cell_scalar( sources_grid,                          &
              &   scalar = sources_ggd(1)%electrons%energy,              &
              &   b2CellData = tmpCv )
          !! ti: (Common) Ion Temperature
          tmpCv(:,:) = ti_mean(:,:)/qe
          call write_quantity( batch_grid,                               &
              &   val = profiles_ggd%t_i_average,                        &
              &   value = tmpCv )
#endif
#if ( IMAS_MAJOR_VERSION < 4 || ( IMAS_MAJOR_VERSION == 4 && IMAS_MINOR_VERSION < 1 ) )
          !! Ion energy sources
          tmpCv(:,:) = shi_mean(:,:) / vol(:,:)
          call write_cell_scalar( sources_grid,                          &
              &   scalar = batch_sources%source(1)%ggd( batch_index )%   &
              &            total_ion_energy,                             &
              &   b2CellData = tmpCv )

          js = 0
          do is = 1, nspecies
            !! sna: Neutral particle sources
            tmpCv(:,:) = sna_mean(:,:,js) / vol(:,:)
            call write_cell_scalar( sources_grid,                        &
                &   scalar = batch_sources%source(1)%                    &
                &   ggd( batch_index )%neutral( is )%particles,          &
                &   b2CellData = tmpCv )
            !! smo: Neutral parallel momentum sources
            tmpCv(:,:) = smo_mean(:,:,js) / vol(:,:)
            call write_cell_vector_component( sources_grid,              &
                &   vectorComponent = batch_sources%source(1)%           &
                &                     ggd( batch_index )%                &
                &                     neutral( is )%momentum,            &
                &   b2CellData = tmpCv,                                  &
                &   vectorID = VEC_ALIGN_PARALLEL_ID )
            js = js + nfluids(is) + 1
          end do
#else
          !! Ion energy sources
          tmpCv(:,:) = shi_mean(:,:) / vol(:,:)
          call write_cell_scalar( sources_grid,                          &
              &   scalar = sources_ggd(1)%total_ion_energy,              &
              &   b2CellData = tmpCv )

          js = 0
          do is = 1, nspecies
            !! sna: Neutral particle sources
            tmpCv(:,:) = sna_mean(:,:,js) / vol(:,:)
            call write_cell_scalar( sources_grid,                        &
                &   scalar = sources_ggd(1)%neutral( is )%particles,     &
                &   b2CellData = tmpCv )
            !! smo: Neutral parallel momentum sources
            tmpCv(:,:) = smo_mean(:,:,js) / vol(:,:)
            call write_cell_vector_component( sources_grid,              &
                &   vectorComponent = sources_ggd(1)%                    &
                &                     neutral( is )%momentum,            &
                &   b2CellData = tmpCv,                                  &
                &   vectorID = VEC_ALIGN_PARALLEL_ID )
            js = js + nfluids(is) + 1
          end do
#endif

#if ( IMAS_MAJOR_VERSION < 4 || ( IMAS_MAJOR_VERSION == 4 && IMAS_MINOR_VERSION < 1 ) )
          !! po: Electric potential
          call write_quantity( batch_grid,                               &
              &   val = batch_profiles%ggd( batch_index )%phi_potential, &
              &   value = po_mean )
          !! sch: Current sources
          tmpCv(:,:) = sch_mean(:,:) / vol(:,:)
          call write_cell_scalar( sources_grid,                          &
              &   scalar = batch_sources%source(1)%ggd( batch_index)%    &
              &            current,                                      &
              &   b2CellData = tmpCv )
#else
          !! po: Electric potential
          call write_quantity( batch_grid,                               &
              &   val = profiles_ggd%phi_potential,                      &
              &   value = po_mean )
          !! sch: Current sources
          tmpCv(:,:) = sch_mean(:,:) / vol(:,:)
          call write_cell_scalar( sources_grid,                          &
              &   scalar = sources_ggd(1)%current,                       &
              &   b2CellData = tmpCv )
#endif
#else
          call logmsg( LOGINFO, &
              &   "b2mod_ual_io.B25_av_ids: GGD not available, no averaged plasma state writing" )
#endif
        end if

! Copy in source data
#if ( IMAS_MAJOR_VERSION > 4 || ( IMAS_MAJOR_VERSION == 4 && IMAS_MINOR_VERSION > 0 ) )
        batch_profiles%ggd( batch_index ) = profiles_ggd
        batch_plasma_profiles%ggd( batch_index ) = profiles_ggd
        batch_sources%source(1)%ggd( batch_index ) = sources_ggd(1)
        batch_plasma_sources%source(1)%ggd( batch_index ) = sources_ggd(1)
#endif

        nc = max(nncut,1)
#if ( IMAS_MINOR_VERSION > 21 || IMAS_MAJOR_VERSION > 3 )
        if (do_description) then
! Summary plasma composition
#if ( IMAS_MAJOR_VERSION > 4 || ( IMAS_MAJOR_VERSION == 4 && IMAS_MINOR_VERSION > 0 ) )
          u = 0.0_IDS_real
          do is = 0, ns-1
            do iy = -1, ny
              do ix = -1, nx
                if (.not.on_closed_surface(ix,iy)) cycle
                u = u + na_mean(ix,iy,is)*rza(ix,iy,is)*vol(ix,iy)
              end do
            end do
          end do
          do i = 1, nspecies
            frac = 0.0_IDS_real
            do is = eb2spcr(i), eb2spcr(i)+nfluids(i)+1
              if (is.ge.ns) cycle
              if (is_neutral(is)) cycle
              if (.not.(is.eq.eb2spcr(i).or.lnext(eb2spcr(i),is))) cycle
              do iy = -1, ny
                do ix = -1, nx
                  if (.not.on_closed_surface(ix,iy)) cycle
                  frac = frac + na_mean(ix,iy,is)*vol(ix,iy)
                end do
              end do
            end do
            if (u.gt.0.0_IDS_real) frac = frac / u
            select case (is_codes(eb2spcr(is)))
            case ('H')
              call write_sourced_constant_2( summary%composition%hydrogen, frac )
            case ('D')
              call write_sourced_constant_2( summary%composition%deuterium, frac )
            case ('T')
              call write_sourced_constant_2( summary%composition%tritium, frac )
            case ('DT')
              call write_sourced_constant_2( summary%composition%deuterium_tritium, frac )
            case ('He')
              if (nint(am(eb2spcr(is))).eq.3) then
                call write_sourced_constant_2( summary%composition%helium_3, frac )
              else if (nint(am(eb2spcr(is))).eq.4) then
                call write_sourced_constant_2( summary%composition%helium_4, frac )
              end if
            case ('Li')
              call write_sourced_constant_2( summary%composition%lithium, frac )
            case ('Be')
              call write_sourced_constant_2( summary%composition%beryllium, frac )
            case ('B')
              call write_sourced_constant_2( summary%composition%boron, frac )
            case ('C')
              call write_sourced_constant_2( summary%composition%carbon, frac )
            case ('N')
#if ( IMAS_MAJOR_VERSION == 4 && IMAS_MINOR_VERSION == 1 )
              call write_sourced_constant( summary%composition%nitrogen, frac )
#else
              call write_sourced_constant_2( summary%composition%nitrogen, frac )
#endif
            case ('O')
              call write_sourced_constant_2( summary%composition%oxygen, frac )
            case ('Ne')
              call write_sourced_constant_2( summary%composition%neon, frac )
            case ('Ar')
              call write_sourced_constant_2( summary%composition%argon, frac )
            case ('Fe')
              call write_sourced_constant_2( summary%composition%iron, frac )
            case ('Xe')
              call write_sourced_constant_2( summary%composition%xenon, frac )
            case ('W')
              call write_sourced_constant_2( summary%composition%tungsten, frac )
            case ('Kr')
              call write_sourced_constant_2( summary%composition%krypton, frac )
            end select
          end do
#endif
! Summary separatrix data
          if (maxval(abs(fpsi(-1:nx,-1:ny,0:3))).gt.0.0_R8) then
            allocate( summary%local%separatrix%position%psi( num_batch_slices ) )
            summary%local%separatrix%position%psi( batch_index ) = fpsi(jxa,jsep,2)
#if ( IMAS_MINOR_VERSION > 36 || IMAS_MAJOR_VERSION > 3 )
            allocate( summary%local%separatrix_average%position%psi( num_batch_slices ) )
            summary%local%separatrix_average%position%psi( batch_index ) = fpsi(jxa,jsep,2)
#endif
          end if
          call write_errored_value( summary%local%separatrix%t_e,         &
              &  tesepm_av(nc), tesepm_std(nc) )
          call write_errored_value( summary%local%separatrix%t_i_average, &
              &  tisepm_av(nc), tisepm_std(nc) )
          call write_errored_value( summary%local%separatrix%n_e,         &
              &  nesepm_av(nc), nesepm_std(nc) )
! Summary divertor plate data
          if (ntrgts.gt.0) then
#if ( IMAS_MINOR_VERSION > 34 || IMAS_MAJOR_VERSION > 3 )
            allocate ( summary%local%divertor_target( ntrgts ) )
#else
            allocate ( summary%local%divertor_plate( ntrgts ) )
#endif
            do i = 1, ntrgts
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
            if (nncut.eq.2) then
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
            call write_errored_value( summary%local%divertor_target(ntrgts)%t_e, &
              &  tesepa_av(nc), tesepa_std(nc) )
            call write_errored_value( summary%local%divertor_target(ntrgts)%t_i_average, &
              &  tisepa_av(nc), tisepa_std(nc) )
            call write_errored_value( summary%local%divertor_target(ntrgts)%n_e, &
              &  nesepa_av(nc), nesepi_std(nc) )
#else
            call write_errored_value( summary%local%divertor_plate(1)%t_e, &
              &  tesepi_av(1), tesepi_std(1) )
            call write_errored_value( summary%local%divertor_plate(1)%t_i_average, &
              &  tisepi_av(1), tisepi_std(1) )
            call write_errored_value( summary%local%divertor_plate(1)%n_e, &
              &  nesepi_av(1), nesepi_std(1) )
            if (nncut.eq.2) then
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
            call write_errored_value( summary%local%divertor_plate(ntrgts)%t_e, &
              &  tesepa_av(nc), tesepa_std(nc) )
            call write_errored_value( summary%local%divertor_plate(ntrgts)%t_i_average, &
              &  tisepa_av(nc), tisepa_std(nc) )
            call write_errored_value( summary%local%divertor_plate(ntrgts)%n_e, &
              &  nesepa_av(nc), nesepi_std(nc) )
#endif
          end if
        end if

        call fill_summary_data( summary )
#endif
        call logmsg( LOGDEBUG, "b2mod_ual_io.B25_av_ids: done" )

        return
    end subroutine B25_av_ids

    function separatrix_average( field, weight )
    ! This function is devoted to obtain the weighted average along the active separatrix
    ! of a plasma field quantity
    ! The average is made using face-centered quantities on the cell faces forming the separatrix
    ! The weighting automatically includes the areas of the cell faces
    implicit none
    real(kind=IDS_real) :: separatrix_average
    real(kind=IDS_real), intent(in) :: field(-1:nx,-1:ny), weight(-1:nx,-1:ny)
    real(kind=IDS_real) :: sum, area_sum
    integer ix

    separatrix_average = IDS_REAL_INVALID
    sum = 0.0_IDS_real
    area_sum = 0.0_IDS_real
    do ix = -1, nx
      if (mod(region(ix,jsep,0),4).eq.1) then
        sum = sum + gs(topix(ix,jsep),topiy(ix,jsep),1) * &
            & ( field(ix,jsep)*weight(ix,jsep) +          &
            &   field(topix(ix,jsep),topiy(ix,jsep))*     &
            &  weight(topix(ix,jsep),topiy(ix,jsep)) )/2.0_IDS_real
        area_sum = area_sum + gs(topix(ix,jsep),topiy(ix,jsep),1)
      end if
    end do
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

    subroutine write_ids_code( code, commit, description )
    implicit none
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
    if (streql(b2frates_flag,'adas')) nlibs = nlibs + 1
#ifdef AMNS
    if (streql(b2frates_flag,'amns')) nlibs = nlibs + 1
#endif
#ifdef B25_EIRENE
    if (use_eirene.ne.0) then
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
    if (streql(b2frates_flag,'adas')) then
      nlibs = nlibs + 1
      allocate( code%library( nlibs )%name(1) )
      code%library( nlibs )%name = 'ADAS'
      allocate( code%library( nlibs )%version(1) )
      code%library( nlibs )%version = adas_version
      allocate( code%library( nlibs )%commit(1) )
      code%library( nlibs )%commit = ADAS_git_version
      allocate( code%library( nlibs )%repository(1) )
      repository = "ssh://git.iter.org/imex/amns-adas.git"
      code%library( nlibs )%repository = repository
    end if
#ifdef AMNS
    if (streql(b2frates_flag,'amns')) then
      nlibs = nlibs + 1
      call IMAS_AMNS_SETUP(amns)
      allocate( code%library( nlibs )%name(1) )
      code%library( nlibs )%name = 'AMNS'
      query%string = 'code_version'
      call IMAS_AMNS_QUERY(amns,query,answer,amns_status)
      if (.not.amns_status%flag) then
        allocate( code%library( nlibs )%version(1) )
        code%library( nlibs )%version = answer%string
      end if
      query%string = 'code_commit'
      call IMAS_AMNS_QUERY(amns,query,answer,amns_status)
      if (.not.amns_status%flag) then
        allocate( code%library( nlibs )%commit(1) )
        code%library( nlibs )%commit = answer%string
      end if
      query%string = 'code_repository'
      call IMAS_AMNS_QUERY(amns,query,answer,amns_status)
      if (.not.amns_status%flag) then
        allocate( code%library( nlibs )%repository(1) )
        code%library( nlibs )%repository = answer%string
      end if
      call IMAS_AMNS_FINISH(amns)
    end if
#endif
#ifdef B25_EIRENE
    if (use_eirene.ne.0) then
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

#if ( IMAS_MINOR_VERSION > 29 || IMAS_MAJOR_VERSION > 3 )
    subroutine write_ids_code_timed( code, commit, description )
    implicit none
    type(ids_code_with_timebase), intent(inout) :: code
                !< Type of IDS data structure, designed for code data handling
    character(len=ids_string_length), intent(in) :: commit
    character(len=ids_string_length), intent(in) :: description

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
    call write_timed_integer( code%output_flag, 0 )

    return
    end subroutine write_ids_code_timed

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

    subroutine write_model_identifier( model_id )
    implicit none
    type(ids_identifier) :: model_id
    integer, save :: ncall = 0
    integer, save :: style = 1
    integer, save :: ids_from_43 = 0

    if (ncall.eq.0) then
      call ipgeti ('ids_from_43', ids_from_43)
      call ipgeti ('b2mndt_style', style)
    endif

    allocate( model_id%name(1) )
    allocate( model_id%description(1) )
    if (ids_from_43.eq.0) then
      if (style.eq.0) then
        model_id%index = -2
        model_id%name(1) = "SOLPS5.0"
        model_id%description(1) = "SOLPS5.0 physics model"
      else if (style.ge.1) then
        model_id%index = -3
        model_id%name(1) = "SOLPS5.2"
        model_id%description(1) = "SOLPS5.2 physics model"
      else if (style.eq.-1) then
        model_id%index = -1
        model_id%name(1) = "SOLPS4.3"
        model_id%description(1) = "SOLPS4.3 physics model"
      end if
    else
      model_id%index = -1
      model_id%name(1) = "SOLPS4.3"
      model_id%description(1) = "SOLPS4.3 physics model"
    end if

    ncall = ncall + 1
    return
    end subroutine write_model_identifier

#if ( ( IMAS_MAJOR_VERSION == 3 && IMAS_MINOR_VERSION > 30 ) || ( IMAS_MAJOR_VERSION == 4 && IMAS_MINOR_VERSION == 0 ) )
    subroutine write_source_identifier( source_id, id_index )
    implicit none
    type(ids_identifier) :: source_id
    integer :: id_index

    allocate( source_id%name(1) )
    allocate( source_id%description(1) )
    source_id%index = id_index
    source_id%name = edge_source_identifier%name( id_index )
    source_id%description = edge_source_identifier%description( id_index )

    return
    end subroutine write_source_identifier
#endif

    subroutine put_equilibrium_data ( equilibrium, &
#if ( IMAS_MINOR_VERSION > 21 || IMAS_MAJOR_VERSION > 3 )
       &  summary, &
#endif
       &  edgeprof, database, &
#if ( AL_MAJOR_VERSION > 4 && GGD_MAJOR_VERSION > 0 )
       &  time_sind, &
#endif
       &  do_summary_data, new_eq_ggd )
#if ( IMAS_MINOR_VERSION > 14 || IMAS_MAJOR_VERSION > 3 ) && GGD_MAJOR_VERSION > 0
    use b2mod_ual_io_grid &
       & , only: GGD_copy_AoS3Root_to_Dynamic
#endif
    implicit none
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
    integer :: iy
    real(IDS_real) :: er_Vx( -1:ubound( na, 1), -1:ubound( na, 2) )
    real(IDS_real) :: er_Fc( -1:ubound( na, 1), -1:ubound( na, 2), 0:1)
    real(IDS_real) :: er_Cv( -1:ubound( na, 1), -1:ubound( na, 2) )
    real(IDS_real) :: ez_Vx( -1:ubound( na, 1), -1:ubound( na, 2) )
    real(IDS_real) :: ez_Fc( -1:ubound( na, 1), -1:ubound( na, 2), 0:1)
    real(IDS_real) :: ez_Cv( -1:ubound( na, 1), -1:ubound( na, 2) )
    real(IDS_real) :: tmpVx( -1:ubound( na, 1), -1:ubound( na, 2) )
    real(IDS_real) :: tmpFace( -1:ubound( na, 1), -1:ubound( na, 2), 0:1)
    real(IDS_real) :: tmpCv( -1:ubound( na, 1), -1:ubound( na, 2) )
#endif
    integer :: i, ix, icnt, inode
    integer :: idum(0:3)
    integer, save :: ncall = 0
    real(IDS_real) :: parg(0:99)
    real(IDS_real), save :: pit_rescale = 1.0_IDS_real
    real(IDS_real) :: b0r0_ref, zmid
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
                            &  crx(jxa,-1,2)+crx(jxa,-1,3))/4.0_R8
          b0 = b0r0 / r0
        else if (isymm.eq.3 .or. isymm.eq.4) then
          b0r0 = bb(jxa,-1,2)*(cry(jxa,-1,0)+cry(jxa,-1,1)+ &
                            &  cry(jxa,-1,2)+cry(jxa,-1,3))/4.0_R8
          b0 = b0r0 / r0
        end if
      else
        b0 = bb(jxa,-1,2)
        if (isymm.eq.1 .or. isymm.eq.2) then
          b0r0 = bb(jxa,-1,2)*(crx(jxa,-1,0)+crx(jxa,-1,1)+ &
                            &  crx(jxa,-1,2)+crx(jxa,-1,3))/4.0_R8
        else if (isymm.eq.3 .or. isymm.eq.4) then
          b0r0 = bb(jxa,-1,2)*(cry(jxa,-1,0)+cry(jxa,-1,1)+ &
                            &  cry(jxa,-1,2)+cry(jxa,-1,3))/4.0_R8
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
          call b2agx0 (99, idum(0), idum(1), idum(2), idum(3))
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
            call b2_IMAS_Fill_Grid_Desc( IDSmap, eq_grid,                     &
              &   nx, ny, crx(-1:nx, -1:ny, :), cry(-1:nx, -1:ny, : ),        &
              &   leftix, leftiy, rightix, rightiy, topix, topiy, bottomix,   &
              &   bottomiy, nnreg, topcut, region, cflags,                    &
              &   INCLUDE_GHOST_CELLS, midplane_id, z_eq, vol, gs, qc )
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
            do iy = -1, ny
             do ix = -1, nx
              tmpVx(ix,iy) = crx(ix,iy,0)
              tmpFace(ix,iy,0) = (crx(ix,iy,0) + crx(ix,iy,1))/2.0_IDS_real
              tmpFace(ix,iy,1) = (crx(ix,iy,0) + crx(ix,iy,2))/2.0_IDS_real
              tmpCv(ix,iy) = cr(ix,iy)
             end do
            end do
            call write_vertex_scalar( eq_grid,                                &
                &   scalar = equilibrium%time_slice( slice_index )%ggd(1)%r,  &
                &   b2VertexData = tmpVx )
            call write_face_scalar( eq_grid,                                  &
                &   val = equilibrium%time_slice( slice_index )%ggd(1)%r,     &
                &   value = tmpFace )
            call write_cell_scalar( eq_grid,                                  &
                &   scalar = equilibrium%time_slice( slice_index )%ggd(1)%r,  &
                &   b2CellData = tmpCv )
          end if
          if (.not.associated(                                                &
            &  equilibrium%time_slice( slice_index )%ggd(1)%z ) ) then
            do iy = -1, ny
             do ix = -1, nx
              tmpVx(ix,iy) = cry(ix,iy,0)
              tmpFace(ix,iy,0) = (cry(ix,iy,0) + cry(ix,iy,1))/2.0_IDS_real
              tmpFace(ix,iy,1) = (cry(ix,iy,0) + cry(ix,iy,2))/2.0_IDS_real
              tmpCv(ix,iy) = cz(ix,iy)
             end do
            end do
            call write_vertex_scalar( eq_grid,                                &
                &   scalar = equilibrium%time_slice( slice_index )%ggd(1)%z,  &
                &   b2VertexData = tmpVx )
            call write_face_scalar( eq_grid,                                  &
                &   val = equilibrium%time_slice( slice_index )%ggd(1)%z,     &
                &   value = tmpFace )
            call write_cell_scalar( eq_grid,                                  &
                &   scalar = equilibrium%time_slice( slice_index )%ggd(1)%z,  &
                &   b2CellData = tmpCv )
          end if
          if (maxval(abs(fpsi)).ne.0.0_IDS_real .and. .not.associated(        &
            &  equilibrium%time_slice( slice_index )%ggd(1)%psi ) ) then
            tmpVx(:,:) = fpsi(:,:,0)
            tmpFace(:,:,0) = (fpsi(:,:,0) + fpsi(:,:,1))/2.0_IDS_real
            tmpFace(:,:,1) = (fpsi(:,:,0) + fpsi(:,:,2))/2.0_IDS_real
            tmpCv(:,:) = (fpsi(:,:,0) + fpsi(:,:,1) +                         &
                &         fpsi(:,:,2) + fpsi(:,:,3) )/4.0_IDS_real
            call write_vertex_scalar( eq_grid,                                &
                &   scalar = equilibrium%time_slice( slice_index )%ggd(1)%    &
                &            psi,                                             &
                &   b2VertexData = tmpVx )
            call write_face_scalar( eq_grid,                                  &
                &   val = equilibrium%time_slice( slice_index )%ggd(1)%psi,   &
                &   value = tmpFace )
            call write_cell_scalar( eq_grid,                                  &
                &   scalar = equilibrium%time_slice( slice_index )%ggd(1)%    &
                &            psi,   &
                &   b2CellData = tmpCv )
          end if
          if (maxval(abs(ffbz)).ne.0.0_IDS_real .and. .not.associated(        &
            &  equilibrium%time_slice( slice_index )%ggd(1)%phi ) ) then
            tmpVx(:,:) = ffbz(:,:,0)
            tmpFace(:,:,0) = (ffbz(:,:,0) + ffbz(:,:,1))/2.0_IDS_real
            tmpFace(:,:,1) = (ffbz(:,:,0) + ffbz(:,:,2))/2.0_IDS_real
            tmpCv(:,:) = (ffbz(:,:,0) + ffbz(:,:,1) +                         &
                &         ffbz(:,:,2) + ffbz(:,:,3) )/4.0_IDS_real
            call write_vertex_scalar( eq_grid,                                &
                &   scalar = equilibrium%time_slice( slice_index )%ggd(1)%    &
                &            phi,                                             &
                &   b2VertexData = tmpVx )
            call write_face_scalar( eq_grid,                                  &
                &   val = equilibrium%time_slice( slice_index )%ggd(1)%phi,   &
                &   value = tmpFace )
            call write_cell_scalar( eq_grid,                                  &
                &   scalar = equilibrium%time_slice( slice_index )%ggd(1)%    &
                &            phi,                                             &
                &   b2CellData = tmpCv )
          end if
          if (.not.associated(                                                &
            &  equilibrium%time_slice( slice_index )%ggd(1)%b_field_r ) ) then
            call compute_er( nx, ny, er_Vx, er_Fc, er_Cv )
            tmpVx(-1:nx,-1:ny) = wbbc(-1:nx,-1:ny,0)*er_Vx(-1:nx,-1:ny)
            tmpFace(:,:,0) = wbbl(:,:,0)*er_Fc(:,:,0)
            tmpFace(:,:,1) = wbbv(:,:,0)*er_Fc(:,:,1)
            tmpCv(:,:) = bb(:,:,0)*er_Cv(:,:)
            call write_vertex_scalar( eq_grid,                                &
                &   scalar = equilibrium%time_slice( slice_index )%ggd(1)%    &
                &         b_field_r,                                          &
                &   b2VertexData = tmpVx )
            call write_face_scalar( eq_grid,                                  &
                &   val = equilibrium%time_slice( slice_index )%ggd(1)%       &
                &         b_field_r,                                          &
                &   value = tmpFace )
            call write_cell_scalar( eq_grid,                                  &
                &   scalar = equilibrium%time_slice( slice_index )%ggd(1)%    &
                &         b_field_r,                                          &
                &   b2CellData = tmpCv )
          end if
          if (.not.associated(                                                &
            &  equilibrium%time_slice( slice_index )%ggd(1)%b_field_z ) ) then
            call compute_ez( nx, ny, ez_Vx, ez_Fc, ez_Cv )
            tmpVx(-1:nx,-1:ny) = wbbc(-1:nx,-1:ny,0)*ez_Vx(-1:nx,-1:ny)
            tmpFace(:,:,0) = wbbl(:,:,0)*ez_Fc(:,:,0)
            tmpFace(:,:,1) = wbbv(:,:,0)*ez_Fc(:,:,1)
            tmpCv(:,:) = bb(:,:,0)*ez_Cv(:,:)
            call write_vertex_scalar( eq_grid,                                &
                &   scalar = equilibrium%time_slice( slice_index )%ggd(1)%    &
                &         b_field_z,                                          &
                &   b2VertexData = tmpVx )
            call write_face_scalar( eq_grid,                                  &
                &   val = equilibrium%time_slice( slice_index )%ggd(1)%       &
                &         b_field_z,                                          &
                &   value = tmpFace )
            call write_cell_scalar( eq_grid,                                  &
                &   scalar = equilibrium%time_slice( slice_index )%ggd(1)%    &
                &         b_field_z,                                          &
                &   b2CellData = tmpCv )
          end if
#if ( IMAS_MINOR_VERSION > 41 || IMAS_MAJOR_VERSION > 3 )
          if (.not.associated(                                                &
            &  equilibrium%time_slice( slice_index )%ggd(1)%b_field_phi ) ) then
            tmpVx(-1:nx,-1:ny) = -wbbc(-1:nx,-1:ny,2)
            tmpFace(:,:,0) = -wbbl(:,:,2)
            tmpFace(:,:,1) = -wbbv(:,:,2)
            tmpCv(:,:) = -bb(:,:,2)
            call write_vertex_scalar( eq_grid,                                &
                &   scalar = equilibrium%time_slice( slice_index )%ggd(1)%    &
                &         b_field_phi,                                        &
                &   b2VertexData = tmpVx )
            call write_face_scalar( eq_grid,                                  &
                &   val = equilibrium%time_slice( slice_index )%ggd(1)%       &
                &         b_field_phi,                                        &
                &   value = tmpFace )
            call write_cell_scalar( eq_grid,                                  &
                &   scalar = equilibrium%time_slice( slice_index )%ggd(1)%    &
                &         b_field_phi,                                        &
                &   b2CellData = tmpCv )
          end if
#else
          if (.not.associated(                                                &
            &  equilibrium%time_slice( slice_index )%ggd(1)%b_field_tor ) ) then
            tmpVx(-1:nx,-1:ny) = -wbbc(-1:nx,-1:ny,2)
            tmpFace(:,:,0) = -wbbl(:,:,2)
            tmpFace(:,:,1) = -wbbv(:,:,2)
            tmpCv(:,:) = -bb(:,:,2)
            call write_vertex_scalar( eq_grid,                                &
                &   scalar = equilibrium%time_slice( slice_index )%ggd(1)%    &
                &         b_field_tor,                                        &
                &   b2VertexData = tmpVx )
            call write_face_scalar( eq_grid,                                  &
                &   val = equilibrium%time_slice( slice_index )%ggd(1)%       &
                &         b_field_tor,                                        &
                &   value = tmpFace )
            call write_cell_scalar( eq_grid,                                  &
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
      if ( z_eq.ne.IDS_REAL_INVALID .and. &
         & (cry(jxa,jsep,2)-z_eq)*(cry(jxa,jsep,3)-z_eq).le.0.0_R8 ) then
        midplane_id = 1
      else if ( jxa .eq. nmdpl ) then
        zmid = 0.5_R8*(cry(jxa,jsep,2)+cry(jxa,jsep,3))
        if ( jxi.eq.-2 ) then
          midplane_id = 2
        else if ( (cry(jxi,jsep,2)-zmid)*(cry(jxi,jsep,3)-zmid).le.0.0_R8) then
          midplane_id = 2
        else
          midplane_id = 4
        end if
      else if ( cry(jxa,jsep,2)*cry(jxa,jsep,3).le.0.0_R8 ) then
        midplane_id = 3
      else
        midplane_id = 4
      end if
    end if

    ncall = ncall + 1
    return
    end subroutine put_equilibrium_data

    subroutine fill_neutral_element( is, js, neutral_element )
    implicit none
    integer, intent(in) :: is !< B2.5 species index
    integer, intent(in) :: js !< B2.5 isonuclear species index
    type(ids_plasma_composition_neutral_element) :: neutral_element

    neutral_element%a = am(is)
#if IMAS_MAJOR_VERSION < 4
    neutral_element%z_n = zn(is)
#else
    neutral_element%z_n = nint(zn(is))
#endif
#if ( IMAS_MINOR_VERSION < 15 && IMAS_MAJOR_VERSION < 4 )
    neutral_element%multiplicity = 1.0_IDS_real
#else
    neutral_element%atoms_n = 1
#endif

    return
    end subroutine fill_neutral_element

    subroutine fill_atom_neutral_type( iatm, neutral_type )
    implicit none
    integer, intent(in) :: iatm !< Eirene atom species index or
                                !< B2.5 species index
    type(ids_identifier_dynamic_aos3) :: neutral_type
    character(len=13)  :: spclabel   !< Species label

    if (use_eirene.ne.0) then
      neutral_type%index = -1
      allocate( neutral_type%name(1) )
      neutral_type%name = "Kinetic"
      allocate( neutral_type%description(1) )
      neutral_type%description = "Kinetic neutral atoms from Eirene"
    else
#if ( IMAS_MAJOR_VERSION > 4 || ( IMAS_MAJOR_VERSION == 4 && IMAS_MINOR_VERSION > 0 ) )
      call set_neutral_type_identifier( neutral_type, "Thermal")
#elif ( IMAS_MINOR_VERSION > 30 || IMAS_MAJOR_VERSION > 3 )
      allocate( neutral_type%name(1) )
      allocate( neutral_type%description(1) )
      neutral_type%index = neutrals_identifier%thermal
      neutral_type%name = neutrals_identifier%name( neutrals_identifier%thermal )
      neutral_type%description = neutrals_identifier%description( neutrals_identifier%thermal )
#else
      allocate( neutral_type%name(1) )
      allocate( neutral_type%description(1) )
      neutral_type%index = 2
      neutral_type%name = "Thermal"
      neutral_type%description = "Fluid neutral species from B2.5"
#endif
    end if

    return
    end subroutine fill_atom_neutral_type

#ifdef B25_EIRENE
    subroutine fill_molecule_elements( nelems, imol, neutral_element )
    implicit none
    integer, intent(in) :: nelems !< Number of elements in molecule
    integer, intent(in) :: imol !< Eirene molecular species index
    type(ids_plasma_composition_neutral_element) :: neutral_element(nelems)
    integer :: i, k

    i = 0
    do k = 1, natmi
      if (mlcmp( k, imol ) > 0 ) then
        i = i + 1
        neutral_element(i)%a = nmassa( k )
#if IMAS_MAJOR_VERSION < 4
        neutral_element(i)%z_n = real(nchara(k),IDS_real)
#else
        neutral_element(i)%z_n = nchara(k)
#endif
#if ( IMAS_MINOR_VERSION < 15 && IMAS_MAJOR_VERSION < 4 )
        neutral_element(i)%multiplicity = mlcmp(k,imol)
#else
        neutral_element(i)%atoms_n = mlcmp(k,imol)
#endif
      end if
    end do

    return
    end subroutine fill_molecule_elements

    subroutine fill_molecule_type( neutral_type )
    implicit none
    type(ids_identifier_dynamic_aos3) :: neutral_type

    neutral_type%index = -1
    allocate( neutral_type%name(1) )
    neutral_type%name = "Kinetic"
    allocate( neutral_type%description(1) )
    neutral_type%description = "Kinetic neutral molecules from Eirene"

    return
    end subroutine fill_molecule_type

    subroutine fill_mol_ion_elements( nelems, iion, neutral_element )
    implicit none
    integer, intent(in) :: nelems !< Number of elements in molecule
    integer, intent(in) :: iion !< Eirene molecular ion species index
    type(ids_plasma_composition_neutral_element) :: neutral_element(nelems)
    integer :: i, k

    i = 0
    do k = 1, natmi
      if (micmp( k, iion ) > 0 ) then
        i = i + 1
        neutral_element(i)%a = nmassa( k )
#if IMAS_MAJOR_VERSION < 4
        neutral_element(i)%z_n = real(nchara(k),IDS_real)
#else
        neutral_element(i)%z_n = nchara(k)
#endif
#if ( IMAS_MINOR_VERSION < 15 && IMAS_MAJOR_VERSION < 4 )
        neutral_element(i)%multiplicity = micmp(k,iion)
#else
        neutral_element(i)%atoms_n = micmp(k,iion)
#endif
      end if
    end do

    return
    end subroutine fill_mol_ion_elements
#endif

#if ( IMAS_MINOR_VERSION > 21 || IMAS_MAJOR_VERSION > 3 )
    subroutine fill_summary_data( summary )
    implicit none
    type (ids_summary), intent(inout) :: summary
    integer :: i, ib, ireg, is, ix, ngpn
    integer :: istrai
    integer :: iatm
    integer :: iatm1  !< Hydrogenic atom index in molecule composition
    integer :: iatm2  !< Non-hydrogenic atom index in molecule composition
    real(IDS_real) :: gpff, gsum, gmid, gbot, gtop, area
    logical at_top, at_bot, at_mid

    select case (GeometryType)
    case( GEOMETRY_LIMITER )
      call write_sourced_integer( summary%boundary%type, 0 )
    case( GEOMETRY_SN )
      if ( isymm.eq.1 .or. isymm.eq.2 ) then
        if ( cry(leftcut(1),jsep,3).lt.0.0_R8 ) then
          call write_sourced_integer( summary%boundary%type, 11 )
        else if ( cry(leftcut(1),jsep,3).gt.0.0_R8 ) then
          call write_sourced_integer( summary%boundary%type, 12 )
        end if
      else if ( isymm.eq.3 .or. isymm.eq.4 ) then
        if ( crx(leftcut(1),jsep,3).lt.0.0_R8 ) then
          call write_sourced_integer( summary%boundary%type, 11 )
        else if ( crx(leftcut(1),jsep,3).gt.0.0_R8 ) then
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
    if (LSN) then
      call write_sourced_value( summary%boundary%strike_point_inner_r, crx(-1,topcut(1),1) )
      call write_sourced_value( summary%boundary%strike_point_inner_z, cry(-1,topcut(1),1) )
      call write_sourced_value( summary%boundary%strike_point_outer_r, crx(nx,topcut(1),0) )
      call write_sourced_value( summary%boundary%strike_point_outer_z, cry(nx,topcut(1),0) )
    else
      call write_sourced_value( summary%boundary%strike_point_inner_r, crx(nx,topcut(1),0) )
      call write_sourced_value( summary%boundary%strike_point_inner_z, cry(nx,topcut(1),0) )
      call write_sourced_value( summary%boundary%strike_point_outer_r, crx(-1,topcut(1),1) )
      call write_sourced_value( summary%boundary%strike_point_outer_z, cry(-1,topcut(1),1) )
    endif

    call write_sourced_value( summary%fusion%power, fusion_power*1.0e6_IDS_real/5.0_IDS_real )

    call write_sourced_int_constant( summary%gas_injection_rates%impurity_seeding, 0 )
    gsum = 0.0_IDS_real
    gtop = 0.0_IDS_real
    gmid = 0.0_IDS_real
    gbot = 0.0_IDS_real
    do i = 1, nspecies
      is = eb2spcr(i)
      if (is.lt.0) cycle
      if (.not.is_neutral(is)) cycle
      do ib = 1, nbc
        ireg = region(bc_list_x(1,ib),bc_list_y(1,ib),0)
        if (ireg.eq.0) cycle
        at_top = .false.
        at_bot = .false.
        at_mid = .false.
        select case (GeometryType)
        case (GEOMETRY_LINEAR)
          at_mid = bcchar(ib).eq.'N'.or.bcchar(ib).eq.'W'.or.bcchar(ib).eq.'E'
          at_bot = bcchar(ib).eq.'S'.and.LSN
          at_top = bcchar(ib).eq.'S'.and..not.LSN
        case (GEOMETRY_SN)
          at_mid = bcchar(ib).eq.'N'
          at_bot = bcchar(ib).eq.'S'.and.LSN.and.(ireg.eq.3.or.ireg.eq.4)
          at_top = bcchar(ib).eq.'S'.and..not.LSN.and.(ireg.eq.3.or.ireg.eq.4)
        case (GEOMETRY_LFS_SNOWFLAKE_MINUS, GEOMETRY_LFS_SNOWFLAKE_PLUS)
          at_mid = bcchar(ib).eq.'N'
          at_bot = bcchar(ib).eq.'S'.and.(ireg.eq.3 .or. ireg.eq.4 .or. &
        &                                 ireg.eq.5 .or. ireg.eq.6 .or. ireg.eq.7)
        case (GEOMETRY_CDN, GEOMETRY_DDN_BOTTOM, GEOMETRY_DDN_TOP)
          at_mid = bcchar(ib).eq.'N'
          at_bot = bcchar(ib).eq.'S'.and.(ireg.eq.3.or.ireg.eq.8)
          at_top = bcchar(ib).eq.'S'.and.(ireg.eq.4.or.ireg.eq.7)
        case (GEOMETRY_CYLINDER, GEOMETRY_LIMITER, GEOMETRY_ANNULUS)
          at_mid = bcchar(ib).eq.'N'
        case (GEOMETRY_STELLARATORISLAND)
          at_mid = bcchar(ib).eq.'S'.and.(ireg.eq.3.or.ireg.eq.4)
        end select
        if (.not.(at_top.or.at_bot.or.at_mid)) cycle
        gpff = 0.0_R8
        select case (bccon(is,ib))
        case(5)
          area = 0.0_R8
          do ix = 1, bc_list_size(ib)
            if (bcchar(ib).eq.'S') then
              area = area + gs(bc_list_x(ix,ib),bc_list_y(ix,ib),1)
            else if (bcchar(ib).eq.'N') then
              area = area + gs(topix(bc_list_x(ix,ib),bc_list_y(ix,ib)), &
                             & topiy(bc_list_x(ix,ib),bc_list_y(ix,ib)),1)
            else if (bcchar(ib).eq.'W') then
              area = area + gs(bc_list_x(ix,ib),bc_list_y(ix,ib),0)
            else if (bcchar(ib).eq.'E') then
              area = area + gs(rightix(bc_list_x(ix,ib),bc_list_y(ix,ib)), &
                             & rightiy(bc_list_x(ix,ib),bc_list_y(ix,ib)),0)
            end if
          end do
          gpff = area*conpar(is,ib,1)
        case(6, 8, 13, 16, 18, 19, 22, 23, 26, 27)
          gpff = conpar(is,ib,1)
        case(11, 12)
          if (bcchar(ib).eq.'S') then
            ix = cbirso
            do while (ix.lt.cbirso+cbnrso-1.and.cbrbrk(ix)*nx.lt. &
              & maxval(bc_list_x(1:bc_list_size(ib),ib)+0.5_R8))
              ix = ix + 1
            end do
          else if (bcchar(ib).eq.'N') then
            ix = cbirno
            do while (ix.lt.cbirno+cbnrno-1.and.cbrbrk(ix)*nx.lt. &
              & maxval(bc_list_x(1:bc_list_size(ib),ib)+0.5_R8))
              ix = ix + 1
            end do
          else
            ix = -1
          end if
          if (ix.ge.0.and.use_eirene.eq.0) gpff = cbsna(0,is,ix)
        end select
        if (gpff.eq.0.0_R8) cycle
        gsum = gsum + gpff*zn(is)
        if (at_top) then
          gtop = gtop + gpff*zn(is)
        else if (at_bot) then
          gbot = gbot + gpff*zn(is)
        else if (at_mid) then
          gmid = gmid + gpff*zn(is)
        end if
        select case (is_codes(is))
        case ('H')
          call add_sourced_value( summary%gas_injection_rates%hydrogen, gpff )
        case ('D')
          call add_sourced_value( summary%gas_injection_rates%deuterium, gpff )
        case ('T')
          call add_sourced_value( summary%gas_injection_rates%tritium, gpff )
        case ('DT')
          call add_sourced_value( summary%gas_injection_rates%deuterium, gpff/2.0_R8 )
          call add_sourced_value( summary%gas_injection_rates%tritium, gpff/2.0_R8 )
        case ('He')
          if (nint(am(is)).eq.3) then
            call add_sourced_value( summary%gas_injection_rates%helium_3, gpff*zn(is) )
          else if (nint(am(is)).eq.4) then
            call add_sourced_value( summary%gas_injection_rates%helium_4, gpff*zn(is) )
          end if
        case ('Li')
          call add_sourced_value( summary%gas_injection_rates%lithium, gpff*zn(is) )
          summary%gas_injection_rates%impurity_seeding%value = 1
        case ('Be')
          call add_sourced_value( summary%gas_injection_rates%beryllium, gpff*zn(is) )
          summary%gas_injection_rates%impurity_seeding%value = 1
        case ('C')
          call add_sourced_value( summary%gas_injection_rates%carbon, gpff*zn(is) )
          summary%gas_injection_rates%impurity_seeding%value = 1
        case ('N')
          call add_sourced_value( summary%gas_injection_rates%nitrogen, gpff*zn(is) )
          summary%gas_injection_rates%impurity_seeding%value = 1
        case ('O')
          call add_sourced_value( summary%gas_injection_rates%oxygen, gpff*zn(is) )
          summary%gas_injection_rates%impurity_seeding%value = 1
        case ('Ne')
          call add_sourced_value( summary%gas_injection_rates%neon, gpff*zn(is) )
          summary%gas_injection_rates%impurity_seeding%value = 1
        case ('Ar')
          call add_sourced_value( summary%gas_injection_rates%argon, gpff*zn(is) )
          summary%gas_injection_rates%impurity_seeding%value = 1
        case ('Xe')
          call add_sourced_value( summary%gas_injection_rates%xenon, gpff*zn(is) )
          summary%gas_injection_rates%impurity_seeding%value = 1
        case ('Kr')
          call add_sourced_value( summary%gas_injection_rates%krypton, gpff*zn(is) )
          summary%gas_injection_rates%impurity_seeding%value = 1
        end select
      end do
    end do
    do istrai = 1, nstrai
      if (crcstra(istrai).eq.'C') then
        if (use_eirene.ne.0) then
          ngpn=natmi
        else
          ngpn=nspecies
        end if
        do iatm = 1, ngpn
          gpff = tflux(istrai)*gpfc(iatm,istrai)*zn(eb2atcr(iatm))
          if (gpff.eq.0.0_R8) cycle
          gsum = gsum + gpff
          if (istrai.eq.nesepm_istra) then
            if (ndes.gt.0.0_R8 .or. nesepm_pfr.gt.0.0_R8 .or. &
              & private_flux_puff.gt.0.0_R8) then
              if (.not.LSN) then
                gtop = gtop + gpff
              else
                if (pfrregno1.eq.0 .and. &
                 & (pfrregno2.eq.pfrregno1 .or. &
                 & (pfrregno2.eq.2 .and. (GeometryType.eq.GEOMETRY_SN .or. &
                 &  GeometryType.eq.GEOMETRY_LFS_SNOWFLAKE_MINUS .or. &
                 &  GeometryType.eq.GEOMETRY_LFS_SNOWFLAKE_PLUS)) .or. &
                 & (pfrregno2.eq.5 .and.(GeometryType.eq.GEOMETRY_CDN .or. &
                 &  GeometryType.eq.GEOMETRY_DDN_BOTTOM .or. &
                 &  GeometryType.eq.GEOMETRY_DDN_TOP)))) then
                  gbot = gbot + gpff
                else if (pfrregno1.eq.2 .and. pfrregno2.eq.pfrregno1 .and. &
                 & GeometryType.eq.GEOMETRY_SN) then
                  gbot = gbot + gpff
                else if ((pfrregno1.eq.2 .or. pfrregno1.eq.3).and. &
                 & (GeometryType.eq.GEOMETRY_CDN .or. &
                 &  GeometryType.eq.GEOMETRY_DDN_TOP .or. &
                 &  GeometryType.eq.GEOMETRY_DDN_BOTTOM)) then
                  gtop = gtop + gpff
                end if
              end if
            else if (nesepm_sol.gt.0.0_R8 .or. volrec_sol.gt.0.0_R8 .or. &
              & ndes_sol.gt.0.0_R8 .or. nepedm_sol.gt.0.0_R8) then
              gmid = gmid + gpff
            else if (feedback_strata(latmscl(iatm)-1).eq.istrai) then
              ib = na_feedback_ib(latmscl(iatm)-1)
              ireg = region(bc_list_x(1,ib),bc_list_y(1,ib),0)
              at_top = .false.
              at_bot = .false.
              at_mid = .false.
              select case (GeometryType)
              case (GEOMETRY_LINEAR)
                at_mid = bcchar(ib).eq.'N'.or.bcchar(ib).eq.'W'.or.bcchar(ib).eq.'E'
                at_bot = bcchar(ib).eq.'S'.and.LSN
                at_top = bcchar(ib).eq.'S'.and..not.LSN
              case (GEOMETRY_SN, GEOMETRY_LFS_SNOWFLAKE_MINUS, GEOMETRY_LFS_SNOWFLAKE_PLUS)
                at_mid = bcchar(ib).eq.'N'
                at_bot = bcchar(ib).eq.'S'.and.LSN.and.(ireg.eq.3.or.ireg.eq.4.or. &
        &                                               ireg.eq.5.or.ireg.eq.6.or.ireg.eq.7)
                at_top = bcchar(ib).eq.'S'.and..not.LSN.and.(ireg.eq.3.or.  &
        &                                  ireg.eq.4.or.ireg.eq.5.or.ireg.eq.6.or.ireg.eq.7)
              case (GEOMETRY_CDN, GEOMETRY_DDN_BOTTOM, GEOMETRY_DDN_TOP)
                at_mid = bcchar(ib).eq.'N'
                at_bot = bcchar(ib).eq.'S'.and.(ireg.eq.3.or.ireg.eq.8)
                at_top = bcchar(ib).eq.'S'.and.(ireg.eq.4.or.ireg.eq.7)
              case (GEOMETRY_CYLINDER, GEOMETRY_LIMITER, GEOMETRY_ANNULUS)
                at_mid = bcchar(ib).eq.'N'
              case (GEOMETRY_STELLARATORISLAND)
                at_mid = bcchar(ib).eq.'S'.and.(ireg.eq.3.or.ireg.eq.4)
              end select
              if (at_top) then
                gtop = gtop + gpff
              else if (at_bot) then
                gbot = gbot + gpff
              else if (at_mid) then
                gmid = gmid + gpff
              end if
            else !! FIXME : assign to midplane until we know where the stratum is located
              gmid = gmid + gpff
            end if
          else !! FIXME : assign to midplane until we know where the stratum is located
            gmid = gmid + gpff
          end if
          if (gpfc(iatm,istrai).eq.sum(gpfc(:,istrai))) then
            select case (is_codes(eb2atcr(iatm)))
            case ('H')
              call add_sourced_value( summary%gas_injection_rates%hydrogen, &
                  & tflux(istrai)*zn(eb2atcr(iatm)) )
            case ('D')
              call add_sourced_value( summary%gas_injection_rates%deuterium, &
                  & tflux(istrai)*zn(eb2atcr(iatm)) )
            case ('T')
              call add_sourced_value( summary%gas_injection_rates%tritium, &
                  & tflux(istrai)*zn(eb2atcr(iatm)) )
            case ('DT')
              call add_sourced_value( summary%gas_injection_rates%deuterium, &
                  & tflux(istrai)*zn(eb2atcr(iatm))/2.0_R8 )
              call add_sourced_value( summary%gas_injection_rates%tritium, &
                  & tflux(istrai)*zn(eb2atcr(iatm))/2.0_R8 )
            case ('He')
              if (nint(am(is)).eq.3) then
                call add_sourced_value( summary%gas_injection_rates%helium_3, &
                    & tflux(istrai)*zn(eb2atcr(iatm)) )
              else if (nint(am(is)).eq.4) then
                call add_sourced_value( summary%gas_injection_rates%helium_4, &
                    & tflux(istrai)*zn(eb2atcr(iatm)) )
              end if
            case ('Li')
              call add_sourced_value( summary%gas_injection_rates%lithium, &
                  & tflux(istrai)*zn(eb2atcr(iatm)) )
              summary%gas_injection_rates%impurity_seeding%value = 1
            case ('Be')
              call add_sourced_value( summary%gas_injection_rates%beryllium, &
                  & tflux(istrai)*zn(eb2atcr(iatm)) )
              summary%gas_injection_rates%impurity_seeding%value = 1
            case ('C')
              call add_sourced_value( summary%gas_injection_rates%carbon, &
                  & tflux(istrai)*zn(eb2atcr(iatm)) )
              summary%gas_injection_rates%impurity_seeding%value = 1
            case ('N')
              call add_sourced_value( summary%gas_injection_rates%nitrogen, &
                  & tflux(istrai)*zn(eb2atcr(iatm)) )
              summary%gas_injection_rates%impurity_seeding%value = 1
            case ('O')
              call add_sourced_value( summary%gas_injection_rates%oxygen, &
                  & tflux(istrai)*zn(eb2atcr(iatm)) )
              summary%gas_injection_rates%impurity_seeding%value = 1
            case ('Ne')
              call add_sourced_value( summary%gas_injection_rates%neon, &
                  & tflux(istrai)*zn(eb2atcr(iatm)) )
              summary%gas_injection_rates%impurity_seeding%value = 1
            case ('Ar')
              call add_sourced_value( summary%gas_injection_rates%argon, &
                  & tflux(istrai)*zn(eb2atcr(iatm)) )
              summary%gas_injection_rates%impurity_seeding%value = 1
            case ('Xe')
              call add_sourced_value( summary%gas_injection_rates%xenon, &
                  & tflux(istrai)*zn(eb2atcr(iatm)) )
              summary%gas_injection_rates%impurity_seeding%value = 1
            case ('Kr')
              call add_sourced_value( summary%gas_injection_rates%krypton, &
                  & tflux(istrai)*zn(eb2atcr(iatm)) )
              summary%gas_injection_rates%impurity_seeding%value = 1
            end select
          end if
        end do
        if (count(gpfc(:,istrai).gt.0.0_R8).gt.1) then
          iatm1 = 0
          iatm2 = 0
          do iatm = 1, ngpn
            if (gpfc(iatm,istrai).gt.0.0_R8) then
              if (nint(zn(eb2atcr(iatm))).eq.1) iatm1 = iatm
              if (nint(zn(eb2atcr(iatm))).gt.1) iatm2 = iatm
            end if
          end do
          if (iatm1.eq.0.or.iatm2.eq.0) then
            continue ! molecule not identified
          else if (gpfc(iatm1,istrai)+gpfc(iatm2,istrai).le.0.9999_R8) then
            continue ! case not coded
          else
            if (nint(gpfc(iatm1,istrai)/gpfc(iatm2,istrai)).eq.4) then
              if (nint(am(eb2atcr(iatm1))).eq.1 .and. &
                & nint(am(eb2atcr(iatm2))).eq.12) then ! CH4
                call add_sourced_value( summary%gas_injection_rates%methane, &
                    & tflux(istrai)*10 )
                summary%gas_injection_rates%impurity_seeding%value = 1
              else if (nint(am(eb2atcr(iatm1))).eq.1 .and. &
                &      nint(am(eb2atcr(iatm2))).eq.13) then ! 13CH4
                call add_sourced_value( summary%gas_injection_rates%methane_carbon_13, &
                    & tflux(istrai)*10 )
                summary%gas_injection_rates%impurity_seeding%value = 1
              else if (nint(am(eb2atcr(iatm1))).eq.2 .and. &
                &      nint(am(eb2atcr(iatm2))).eq.12) then ! CD4
                call add_sourced_value( summary%gas_injection_rates%methane_deuterated, &
                    & tflux(istrai)*10 )
                summary%gas_injection_rates%impurity_seeding%value = 1
              else if (nint(am(eb2atcr(iatm1))).eq.1 .and. &
                &      nint(am(eb2atcr(iatm2))).eq.28) then ! SiH4
                call add_sourced_value( summary%gas_injection_rates%silane, &
                    & tflux(istrai)*18 )
                summary%gas_injection_rates%impurity_seeding%value = 1
              end if
            else if (nint(gpfc(iatm1,istrai)/gpfc(iatm2,istrai)).eq.2) then
              if (nint(am(eb2atcr(iatm1))).eq.1 .and. &
                & nint(am(eb2atcr(iatm2))).eq.12) then ! C2H4
                call add_sourced_value( summary%gas_injection_rates%ethylene, &
                    & tflux(istrai)*16 )
                summary%gas_injection_rates%impurity_seeding%value = 1
              endif
            else if (nint(gpfc(iatm1,istrai)/gpfc(iatm2,istrai)).eq.3) then ! C2H6, C3H8, NH3, ND3
              if (nint(am(eb2atcr(iatm1))).eq.1 .and. &
                & nint(am(eb2atcr(iatm2))).eq.12) then
                if (nint(gpfc(iatm2,istrai)*6).eq.2) then ! C2H6
                  call add_sourced_value( summary%gas_injection_rates%ethane, &
                      & tflux(istrai)*18 )
                  summary%gas_injection_rates%impurity_seeding%value = 1
                else if (nint(gpfc(iatm2,istrai)*8).eq.3) then ! C3H8
                  call add_sourced_value( summary%gas_injection_rates%propane, &
                      & tflux(istrai)*26 )
                  summary%gas_injection_rates%impurity_seeding%value = 1
                end if
              else if (nint(am(eb2atcr(iatm1))).eq.1 .and. &
                &      nint(am(eb2atcr(iatm2))).eq.14) then ! NH3
                call add_sourced_value( summary%gas_injection_rates%ammonia, &
                    & tflux(istrai)*10 )
                summary%gas_injection_rates%impurity_seeding%value = 1
              else if (nint(am(eb2atcr(iatm1))).eq.2 .and. &
                &      nint(am(eb2atcr(iatm2))).eq.14) then ! ND3
                call add_sourced_value( summary%gas_injection_rates%ammonia_deuterated, &
                    & tflux(istrai)*10 )
                summary%gas_injection_rates%impurity_seeding%value = 1
              end if
            end if
          end if
        end if
      end if
    end do
    call write_sourced_value( summary%gas_injection_rates%total, gsum )
    call write_sourced_value( summary%gas_injection_rates%midplane, gmid )
    call write_sourced_value( summary%gas_injection_rates%top, gtop )
    call write_sourced_value( summary%gas_injection_rates%bottom, gbot )

    return
    end subroutine fill_summary_data
#endif

#if ( IMAS_MINOR_VERSION > 32 || IMAS_MAJOR_VERSION > 3 )
    subroutine write_ids_midplane( midplane, midplane_id )
    implicit none
    type(ids_identifier_static) :: midplane
    integer, intent(in) :: midplane_id

#if ( IMAS_MAJOR_VERSION > 4 || ( IMAS_MAJOR_VERSION == 4 && IMAS_MINOR_VERSION > 0 ) )
    call set_midplane_identifier( midplane, get_midplane_name(midplane_id) )
#else
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
#endif
    return

    end subroutine write_ids_midplane
#endif

#if ( IMAS_MINOR_VERSION > 11 || IMAS_MAJOR_VERSION > 3 ) && GGD_MAJOR_VERSION > 0
    !> Write scalar B2 cell quantity to 'ids_generic_grid_scalar'
    !! IMAS IDS data tree node.
    subroutine write_quantity( basegrid, val, value )
    use b2mod_interp
    implicit none
#if ( IMAS_MINOR_VERSION < 15 && IMAS_MAJOR_VERSION < 4 )
    type(ids_generic_grid_dynamic), intent(in) :: basegrid !< Type of IDS
        !< data structure, designed for handling grid geometry data
#else
    type(ids_generic_grid_aos3_root), intent(in) :: basegrid !< Type of IDS
        !< data structure, designed for handling grid geometry data
#endif
    type(ids_generic_grid_scalar), pointer, intent(inout) :: val(:)
        !< Type of IDS data structure, designed for scalar data handling
    real(IDS_real), intent(in) :: value( -1:IDSmap%b2nx, -1:IDSmap%b2ny )
    real(IDS_real), dimension(:), pointer :: idsdata    !< Array for
        !< handling data field values
    real(IDS_real) :: weight( -1:IDSmap%b2nx, -1:IDSmap%b2ny, TO_SELF:TO_TOP )
    real(IDS_real) :: tmpFace( -1:ubound( na, 1), -1:ubound( na, 2), 0:1)
    real(IDS_real) :: tmpVx( -1:ubound( na, 1), -1:ubound( na, 2) )
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
    tmpVx = interpolateToVertices(  &
          &   IDSmap%b2nx, IDSmap%b2ny, VX_LOWERLEFT, value )

    !! Interpolate data to cell faces, using a volume weighting
    tmpFace = 0.0_IDS_real
    do i = TO_SELF, TO_TOP
       weight(:,:,i) = vol(:,:)
    end do
    call value_on_faces( nx, ny, weight, value, tmpFace)

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
            & GRID_SUBSET_INNER_MIDPLANE, GRID_SUBSET_OUTER_MIDPLANE, &
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
                     &   basegrid, iSubset, IDSmap, tmpVx )
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
                     &   basegrid, iSubset, IDSmap, tmpFace )
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
        idsdata => b2_IMAS_Transform_Data_B2_To_IDS(               &
                      &  basegrid, iSubset, IDSmap, value )
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
    end subroutine write_quantity

    !> Write a scalar B2 cell quantity to ids_generic_grid_scalar
    subroutine write_cell_scalar( basegrid, scalar, b2CellData )
    implicit none
#if ( IMAS_MINOR_VERSION < 15 && IMAS_MAJOR_VERSION < 4 )
    type(ids_generic_grid_dynamic), intent(in) :: basegrid !< Type of IDS
        !< data structure, designed for handling grid geometry data
#else
    type(ids_generic_grid_aos3_root), intent(in) :: basegrid !< Type of IDS
        !< data structure, designed for handling grid geometry data
#endif
    type(ids_generic_grid_scalar), intent(inout), pointer :: scalar(:)
        !< Type of IDS data structure, designed for scalar data handling
    real(IDS_real), intent(in) :: b2CellData(-1:IDSmap%b2nx, -1:IDSmap%b2ny)
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
            & GRID_SUBSET_INNER_MIDPLANE, GRID_SUBSET_OUTER_MIDPLANE, &
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
            & GRID_SUBSET_OUTER_DIVERTOR, GRID_SUBSET_INNER_DIVERTOR, &
            & GRID_SUBSET_OUTER_DIVERTOR_INACTIVE, &
            & GRID_SUBSET_INNER_DIVERTOR_INACTIVE )
          ndim = 3
        case( GRID_SUBSET_VOLUMES )
          ndim = 4
        end select
      end if
      if (ndim.ne.3) cycle

      idsdata => b2_IMAS_Transform_Data_B2_To_IDS(  &
          &   basegrid, iSubset, IDSmap, b2CellData )
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
    subroutine write_cell_vector_component( basegrid, &
       &  vectorComponent, b2CellData, vectorID )
    implicit none
#if ( IMAS_MINOR_VERSION < 15 && IMAS_MAJOR_VERSION < 4 )
    type(ids_generic_grid_dynamic), intent(in) :: basegrid !< Type of IDS
        !< data structure, designed for handling grid geometry data
#else
    type(ids_generic_grid_aos3_root), intent(in) :: basegrid !< Type of IDS
        !< data structure, designed for handling grid geometry data
#endif
    type(ids_generic_grid_vector_components), intent(inout),    &
              &   pointer :: vectorComponent(:) !< Type of IDS data structure,
                    !> designed for vector data handling
    real(IDS_real), intent(in) :: b2CellData(-1:IDSmap%b2nx, -1:IDSmap%b2ny)
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
            & GRID_SUBSET_INNER_MIDPLANE, GRID_SUBSET_OUTER_MIDPLANE, &
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
            & GRID_SUBSET_OUTER_DIVERTOR, GRID_SUBSET_INNER_DIVERTOR, &
            & GRID_SUBSET_OUTER_DIVERTOR_INACTIVE, &
            & GRID_SUBSET_INNER_DIVERTOR_INACTIVE )
          ndim = 3
        case( GRID_SUBSET_VOLUMES )
          ndim = 4
        end select
      end if
      if (ndim.ne.3) cycle
      idsdata => b2_IMAS_Transform_Data_B2_To_IDS(           &
                   &   basegrid, iSubset, IDSmap, b2CellData )
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
    subroutine write_face_scalar( basegrid, val, value )
    implicit none
#if ( IMAS_MINOR_VERSION < 15 && IMAS_MAJOR_VERSION < 4 )
    type (ids_generic_grid_dynamic), intent(in) :: basegrid !< Type of IDS
        !< data structure, designed for handling grid geometry data
#else
    type (ids_generic_grid_aos3_root), intent(in) :: basegrid !< Type of IDS
        !< data structure, designed for handling grid geometry data
#endif
    type (ids_generic_grid_scalar), pointer, intent(inout) :: val(:)
        !< Type of IDS data structure, designed for scalar data handling
        !< (in this case scalars residing on grid faces)
    real(IDS_real), intent(in) :: value( -1:IDSmap%b2nx, -1:IDSmap%b2ny, 0:1 )
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
             & GRID_SUBSET_INNER_MIDPLANE, GRID_SUBSET_OUTER_MIDPLANE, &
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
             & GRID_SUBSET_OUTER_DIVERTOR, GRID_SUBSET_INNER_DIVERTOR, &
             & GRID_SUBSET_OUTER_DIVERTOR_INACTIVE, &
             & GRID_SUBSET_INNER_DIVERTOR_INACTIVE )
           ndim = 3
         case( GRID_SUBSET_VOLUMES )
           ndim = 4
         end select
       end if
       if (ndim.ne.2) cycle
       call write_face_vector( basegrid, val( iSubset ), value, &
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
    subroutine write_face_vector_component( basegrid, &
       &   vectorComponent, b2FaceData, vectorID )
    implicit none
#if ( IMAS_MINOR_VERSION < 15 && IMAS_MAJOR_VERSION < 4 )
    type(ids_generic_grid_dynamic), intent(in) :: basegrid !< Type of IDS
        !< data structure, designed for handling grid geometry data
#else
    type(ids_generic_grid_aos3_root), intent(in) :: basegrid !< Type of IDS
        !< data structure, designed for handling grid geometry data
#endif
    type(ids_generic_grid_vector_components), intent(inout),    &
       &   pointer :: vectorComponent(:) !< Type of IDS data structure,
                                         !> designed for vector data handling
    real(IDS_real), intent(in) ::  &
       &   b2FaceData(-1:IDSmap%b2nx, -1:IDSmap%b2ny, 0:1)
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
            & GRID_SUBSET_INNER_MIDPLANE, GRID_SUBSET_OUTER_MIDPLANE, &
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
            & GRID_SUBSET_OUTER_DIVERTOR, GRID_SUBSET_INNER_DIVERTOR, &
            & GRID_SUBSET_OUTER_DIVERTOR_INACTIVE, &
            & GRID_SUBSET_INNER_DIVERTOR_INACTIVE )
          ndim = 3
        case( GRID_SUBSET_VOLUMES )
          ndim = 4
        end select
      end if
      if (ndim.ne.2) cycle
      idsdata => b2_IMAS_Transform_Data_B2_To_IDS( &
                   &   basegrid, iSubset, IDSmap, b2FaceData )
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
    subroutine write_vertex_scalar( basegrid, scalar, b2VertexData )
    implicit none
#if ( IMAS_MINOR_VERSION < 15 && IMAS_MAJOR_VERSION < 4 )
    type(ids_generic_grid_dynamic), intent(in) :: basegrid !< Type of IDS
        !< data structure, designed for handling grid geometry data
#else
    type(ids_generic_grid_aos3_root), intent(in) :: basegrid !< Type of IDS
        !< data structure, designed for handling grid geometry data
#endif
    type(ids_generic_grid_scalar), intent(inout), pointer :: scalar(:)
        !< Type of IDS data structure, designed for scalar data handling
    real(IDS_real), intent(in) :: b2VertexData(-1:IDSmap%b2nx, -1:IDSmap%b2ny)
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
            & GRID_SUBSET_INNER_MIDPLANE, GRID_SUBSET_OUTER_MIDPLANE, &
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
          &   basegrid, iSubset, IDSmap, b2VertexData )
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
    subroutine write_face_vector( basegrid, vector, b2FaceData, &
       &   gridID, gridSubsetID, gridSubsetInd )
    implicit none
#if ( IMAS_MINOR_VERSION < 15 && IMAS_MAJOR_VERSION < 4 )
    type(ids_generic_grid_dynamic), intent(in) :: basegrid !< Type of IDS
        !< data structure, designed for handling grid geometry data
#else
    type(ids_generic_grid_aos3_root), intent(in) :: basegrid !< Type of IDS
        !< data structure, designed for handling grid geometry data
#endif
    type(ids_generic_grid_scalar), intent(inout) :: vector
        !< Type of IDS data structure, designed for scalar data handling
        !< (in this case 1D vector)
    real(IDS_real), intent(in) :: &
        &   b2FaceData(-1:IDSmap%b2nx, -1:IDSmap%b2ny, 0:1)
    integer, intent(in) :: gridID                    !< Grid identifier index
    integer, intent(in), optional :: gridSubsetID    !< Grid subset identifier index
    integer, intent(in), optional :: gridSubsetInd   !< Base grid subset index
    real(IDS_real), dimension(:), pointer :: idsdata !< Dummy array
        !< for holding data field values

    if ( .not. present(gridSubsetInd) ) then
      !! Fill in vector component data
      idsdata => b2_IMAS_Transform_Data_B2_To_IDS(    &
               &   basegrid, GRID_SUBSET_Y_ALIGNED_EDGES, IDSmap, b2FaceData)
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
      idsdata => b2_IMAS_Transform_Data_B2_To_IDS(    &
               &   basegrid, GRID_SUBSET_X_ALIGNED_EDGES, IDSmap, b2FaceData)
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
      idsdata => b2_IMAS_Transform_Data_B2_To_IDS(    &
               &   basegrid, gridSubsetInd, IDSmap, b2FaceData)
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

    subroutine write_sourced_constant_2( val, value )
    implicit none
    type(ids_summary_constant_flt_0d_2) :: val
        !< Type of IDS data structure, designed for sourced real constant data handling
    real(IDS_real), intent(in) :: value

    val%value = value
    allocate( val%source(1) )
    val%source = source

    return
    end subroutine write_sourced_constant_2

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
    !> From the B2 grid, compute the coordinate unit vectors
    !> (poloidal, radial, toroidal)
    subroutine compute_Coordinate_Unit_Vectors( crx, cry, e1, e2, e3 )
    implicit none
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

    e1 = 0.0_IDS_real
    e2 = 0.0_IDS_real
    e3 = 0.0_IDS_real

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
        dir = 1.0_IDS_real
        ixn = rightix( ix, iy )
        iyn = rightiy( ix, iy )

        if ( .not. isInDomain( nx, ny, ixn, iyn ) ) then
          !! If not found, try to find left neighbour
          !! ...and note to invert vector direction
          dir = -1.0_IDS_real
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
        e1(ix,iy,2) = 0.0_IDS_real    !! phi
        e1(ix,iy,3) = cN(1) - cC(1)   !! Z

        e1(ix,iy,:) = e1(ix,iy,:) * dir  !! fix direction

        !! radial direction
        !! Try to find top neighbour
        dir = 1.0_IDS_real
        ixn = topix( ix, iy )
        iyn = topiy( ix, iy )

        if ( .not. isInDomain( nx, ny, ixn, iyn ) ) then
          !! If not found, try to find bottom neighbour
          !! ...and note to invert vector direction
          dir = -1.0_IDS_real
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
        e2(ix,iy,2) = 0.0_IDS_real    !! phi
        e2(ix,iy,3) = cN(1) - cC(1)   !! Z

        e2(ix,iy,:) = e2(ix,iy,:) * dir  !! fix direction

        !! toroidal direction
        e3(ix,iy,1) = 0.0_IDS_real   !! R
        e3(ix,iy,2) = 1.0_IDS_real   !! phi
        e3(ix,iy,3) = 0.0_IDS_real   !! Z

        !! make unit vectors
        e1(ix,iy,:) = unitVector(e1(ix,iy,:))
        e2(ix,iy,:) = unitVector(e2(ix,iy,:))
        e3(ix,iy,:) = unitVector(e3(ix,iy,:))

      end do
    end do

    return
    end subroutine compute_Coordinate_Unit_Vectors

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
       call shrink_label(spclabel)
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

        !! ue
        allocate(edgecpo%fluid%ve)
        allocate(edgecpo%fluid%ve%comps(1))
        allocate(edgecpo%fluid%ve%align(1))
        allocate(edgecpo%fluid%ve%alignid(1))
        edgecpo%fluid%ve%align(1) = VEC_ALIGN_PARALLEL
        edgecpo%fluid%ve%alignid(1) = VEC_ALIGN_PARALLEL_ID

        call write_cell_scalar( edgecpo%fluid%ve%comps(1)%value, &
            &   b2CellData = ue(:,:) )

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
      call write_face_vector( fluxes(1), flux )
      call write_face_vector( fluxes(2), flux, subgridInd = iSgCore )
    end subroutine write_quantity


    !> Write a scalar B2 cell quantity to a complexgrid_scalar
    subroutine write_cell_scalar(scalar, b2CellData)
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
    end subroutine write_cell_scalar


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
    subroutine write_face_vector(vector, b2FaceData, subgridInd)
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
