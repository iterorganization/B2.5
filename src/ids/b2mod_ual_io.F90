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
    use b2mod_rates
    use b2mod_plasma
    use b2mod_elements
    use b2mod_constants
    use b2mod_sources
    use b2mod_feedback
    use b2mod_transport
    use b2mod_anomalous_transport
    use b2mod_boundary_namelist
    use b2mod_neutrals_namelist
    use b2mod_user_namelist
    use b2mod_indirect
    use b2mod_external
    use b2mod_interp
    use b2mod_ipmain
    use b2mod_b2cmrc
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
    use eirmod_cinit &
     & , only : fort_lc
    use eirmod_comusr &
     & , only : natmi, nmoli, nioni, nmassa, nchara, nmassm, ncharm, &
     &          nprt, nchrgi, nchari
    use b2mod_b2plot &
     & , only : triangle_vol, ix_e2b
#else
#ifdef IMAS
    use b2mod_b2plot &
     & , only : natmi
#endif
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
     & , only : b2_IMAS_Fill_Grid_Desc
    use ids_grid_subgrid  &     ! IGNORE
     & , only : findGridSubsetByName
    use ids_grid_structured &   ! IGNORE
     & , only : GridWriteData
    use ids_grid_common , &     ! IGNORE
        &   IDS_COORDTYPE_R => COORDTYPE_R,       &
        &   IDS_COORDTYPE_Z => COORDTYPE_Z,       &
        &   IDS_GRID_UNDEFINED => GRID_UNDEFINED
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
#if GGD_MINOR_VERSION < 10
    use b2mod_ual_io_grid &
     & , only : GRID_SUBSET_X_ALIGNED_EDGES, GRID_SUBSET_Y_ALIGNED_EDGES, &
     &          GRID_SUBSET_EDGES
#endif
#if IMAS_MINOR_VERSION > 8
    use ids_schemas &     ! IGNORE
     & , only : ids_real, ids_real_invalid
#endif
    use ids_schemas &     ! IGNORE
     & , only : ids_edge_profiles, ids_edge_sources, ids_edge_transport,    &
     &          ids_radiation, ids_dataset_description, ids_equilibrium,    &
     &          ids_ids_properties, &
     &          ids_code, ids_signal_int_1d, ids_signal_flt_1d,             &
     &          ids_generic_grid_scalar, ids_generic_grid_vector_components
#if IMAS_MINOR_VERSION > 21
    use ids_schemas &     ! IGNORE
     & , only : ids_summary,                                                        &
     &          ids_summary_constant_flt_0d, ids_summary_constant_int_0d,           &
     &          ids_summary_dynamic_int_1d_root, ids_summary_dynamic_flt_1d_root,   &
     &          ids_summary_dynamic_flt_1d_root_parent_2, ids_summary_static_str_0d
#endif
#if ( IMAS_MINOR_VERSION > 25 && IMAS_MINOR_VERSION < 34 )
    use ids_schemas &     ! IGNORE
     & , only : ids_numerics
#endif
#if IMAS_MINOR_VERSION > 30
    use ids_schemas &     ! IGNORE
     & , only : ids_divertors
#endif
#if IMAS_MINOR_VERSION > 32
    use ids_utilities &   ! IGNORE
     & , only : ids_identifier_static
#endif
#if IMAS_MINOR_VERSION > 29
#ifdef AMNS
    use amns_types  ! IGNORE
    use amns_module ! IGNORE
#endif
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

  interface write_sourced_value
     module procedure write_sourced_value_root
     module procedure write_sourced_value_root_parent_2
  end interface write_sourced_value

  integer :: num_time_slices  !< Total number of time slices.
  integer :: time_sind   !< Time slice index. Also General grid
            !< description slice identifier
  character(len=132) :: source    !< Code source

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
            &   radiation, description, equilibrium, &
#if IMAS_MINOR_VERSION > 21
            &   summary, &
#endif
#if ( IMAS_MINOR_VERSION > 25 && IMAS_MINOR_VERSION < 34 )
            &   numerics, run_start_time_IN, run_end_time_IN, &
#endif
#if IMAS_MINOR_VERSION > 30
            &   divertors, &
#endif
            &   time_IN, time_step_IN, shot, run, database, version, &
            &   time_slice_ind_IN, num_time_slices_IN )
#ifdef NO_OPT
!DIR$ NOOPTIMIZE
#endif
#include <DIMENSIONS.F>
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
#if IMAS_MINOR_VERSION > 21
        type (ids_summary) :: summary !< IDS designed to store
            !< run summary data
#endif
#if ( IMAS_MINOR_VERSION > 25 && IMAS_MINOR_VERSION < 34 )
        type (ids_numerics) :: numerics !< IDS designed to store
            !< run numerics data
        real(IDS_real), intent(in) :: run_start_time_IN, run_end_time_IN !< Run time bounds
#endif
#if IMAS_MINOR_VERSION > 30
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

        !! Internal variables
        character(len=132) :: comment   !< IDS properties label
        character(len=132) :: username  !< IDS user name
        character(len=132) :: plate_name(4) !< Divertor plate name
        character(len=13)  :: spclabel   !< Species label
        integer :: ion_charge_int !< Ion charge (e.g. 1, 2, etc.)
        integer :: ns    !< Total number of B2.5 species
        integer :: nsion !< Total number of IDS ion species
        integer :: nneut !< Total number of IDS neutral species
        integer :: nx    !< Specifies the number of interior cells
                         !< along the first coordinate
        integer :: ny    !< Specifies the number of interior cells
                         !< along the second coordinate
        integer :: n_process !< Number of radiation processes handled
        integer :: is, js, ks !< Species indices (iterators)
        integer :: ii, jj !< Iterators
        integer :: i      !< Iterator
        integer :: j      !< Iterator
        integer :: k      !< Iterator
        integer :: ix     !< Iterator
        integer :: iy     !< Iterator
        integer :: istrai !< Stratum iterator
        integer :: ireg   !< Region index
        integer :: ntimes !< Number of previous timesteps in IDS
        integer :: is1    !< First ion of an isonuclear sequence
        integer :: is2    !< Last ion of an isonuclear sequence
        integer :: icnt   !< Boundary cell counter
        integer :: ib     !< Boundary condition index
        integer :: ntrgts !< Number of divertor targets
        integer :: p      !< Dummy integer
        integer :: isep(2) !< Array of separatrix regions
        integer :: itrg(4) !< Array of target indices
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
        integer :: iatm1  !< Hydrogenic atom index in molecule composition
        integer :: iatm2  !< Non-hydrogenic atom index in molecule composition
        integer :: nelems !< Number of elements present in a molecule or molecular ion
        integer, allocatable :: isstat(:) !< Mapping array
                                          !< from Eirene atoms and molecules to IDS neutral states
        integer, allocatable :: imneut(:) !< Mapping array
                                          !< from Eirene molecules to IDS neutrals
        integer, allocatable :: imiion(:) !< Mapping array
                                          !< from Eirene molecular ions to IDS ion sequences
        integer :: ixx, iyy
#endif
        integer :: nscx, iscx(0:nscxmax-1)
        integer :: ixpos(4), ifpos(4), iypos(4) !< Target positions
        integer :: idir(4), iysep(4), ixmid(4), ixmax(4)
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
        integer :: midplane_id      !< Location of midplane:
                                    !< 1: Z equal to equilibrium O-point
                                    !< 2: Z at location of maximum major radius
                                    !< 3: Z at dR/dZ = 0 maximum R location
                                    !< 4: GGD grid subset defined by jxa value
        logical, parameter :: B2_WRITE_DATA = .true.
        real(IDS_real),   &
            &   dimension( -1:ubound( crx, 1 ), -1:ubound( crx, 2), 3, 3) :: e
        real(IDS_real) :: flxFace( -1:ubound( na, 1), -1:ubound( na, 2), 0:1)
        real(IDS_real) :: tmpFace( -1:ubound( na, 1), -1:ubound( na, 2), 0:1)
        real(IDS_real) :: totFace( -1:ubound( na, 1), -1:ubound( na, 2), 0:1)
        real(IDS_real) :: tmpVx( -1:ubound( na, 1), -1:ubound( na, 2) )
        real(IDS_real) :: tmpCv( -1:ubound( na, 1), -1:ubound( na, 2) )
        real(IDS_real) :: totCv( -1:ubound( na, 1), -1:ubound( na, 2) )
        real(IDS_real) :: lnlam( -1:ubound( na, 1), -1:ubound( na, 2) )
        real(IDS_real) :: pz( -1:ubound( na, 1), -1:ubound( na, 2) )
        real(IDS_real) :: pb( -1:ubound( na, 1), -1:ubound( na, 2) )
        real(IDS_real) :: pe( -1:ubound( na, 1), -1:ubound( na, 2) )
        real(IDS_real) :: ue( -1:ubound( na, 1), -1:ubound( na, 2) )
        real(IDS_real) :: zeff( -1:ubound( na, 1), -1:ubound( na, 2) )
        real(IDS_real) :: time  !< Generic time
        real(IDS_real) :: time_step !< Time step
        real(IDS_real) :: time_slice_value   !< Time slice value
        real(IDS_real) :: b0, r0, b0r0, b0r0_ref, nibnd, frac, u, v,         &
            &             qtot, qetot, qitot, qmax, qemax, qimax, lambda,    &
            &             vtor, nisep, nasum, area
        real(IDS_real) :: gpff, gsum, gmid, gbot, gtop
        real(IDS_real) :: r_min, r_max, z_min, z_max, z_eq
        real(IDS_real) :: flux_expansion(4), extension_r(4), extension_z(4), &
            &             wetted_area(4), power_convected(4),                &
            &             power_conducted(4), power_neutrals(4),             &
            &             power_incident(4), power_flux_peak(4),             &
            &             power_recomb_neutrals(4), power_radiated(4),       &
            &             recycled_flux(4)
        real(IDS_real), allocatable :: wrdtrg(:,:,:)
#ifdef B25_EIRENE
        real(IDS_real), allocatable :: un0(:,:,:,:), um0(:,:,:,:)
#endif

        type(B2GridMap) :: gmap !< Data structure holding an
            !< intermediate grid description to be transferred into a CPO or IDS

        integer, parameter :: nsources = 12
        integer, save :: style = 1
        integer, save :: ismain = 1
        integer, save :: ismain0 = 0
        integer, save :: ue_style = 2
        integer, save :: use_eirene = 0
        integer, save :: pfrregno1 = 0
        integer, save :: pfrregno2 = 2
        integer, save :: ids_from_43 = 0
        integer, save :: target_offset = 1
        integer, save :: nesepm_istra = -1
        integer, save :: balance_netcdf = 0
        integer, save :: drift_style
        real(IDS_real), save :: dtim = 1.0_IDS_real
        real(IDS_real), save :: ndes = 0.0_IDS_real
        real(IDS_real), save :: ndes_sol = 0.0_IDS_real
        real(IDS_real), save :: nesepm_pfr = 0.0_IDS_real
        real(IDS_real), save :: nesepm_sol = 0.0_IDS_real
        real(IDS_real), save :: nepedm_sol = 0.0_IDS_real
        real(IDS_real), save :: volrec_sol = 0.0_IDS_real
        real(IDS_real), save :: pit_rescale = 1.0_IDS_real
        real(IDS_real), save :: private_flux_puff = 0.0_IDS_real
        real(IDS_real), save :: neutral_sources_rescale = 1.0_IDS_real
        real(IDS_real), save :: BoRiS = 0.0_IDS_real
        integer :: idum(0:3)
        real(IDS_real) :: parg(0:99)
        character*8 date
        character*10 ctime
        character*5 zone
        character*132 create_date
        integer tvalues(8)
        character*16 usrnam
        character*8 imas_version, ual_version, adas_version
        character*32 B25_git_version
        character*32 ADAS_git_version
        character*32 get_B25_hash
        character*32 get_ADAS_hash
        character*8 id
        character*80 cnamip, cvalip
        character*132 code_commit, radiation_commit, eq_source
        character*256 filename
        character*500 line, ligne
        logical match_found, streql, exists, wrong_flow, eq_found
        logical at_top, at_bot, at_mid, target_east, target_west
#ifdef B25_EIRENE
        character(len=132) :: mol_label !< Molecule species label (e.g. D2)
        character(len=132) :: ion_label !< Ion species label (e.g. D+1)
        logical, allocatable :: in_species(:)
#endif
#ifdef USE_PXFGETENV
        integer lenval, ierror
#else
#ifdef NAGFOR
        integer lenval, ierror
#endif
#endif
        logical is_comment
        external is_comment
        external b2xpne, b2xpni, b2xppb, b2xppe, b2xppz, b2xpve, b2xzef
        external b2agx0, b2sral, b2spcx, b2tral, b2tanml, b2ptrdl
        external ipgetr, ipgeti, species, usrnam, streql, xerrab, xertst
        external find_file, strip_spaces, get_B25_hash, get_ADAS_hash

        !! ===  SET UP IDS ===
        write(0,*) "Setting data for edge_profiles IDS"
        call ipgetr ('b2news_BoRiS', BoRiS)
        call ipgeti ('b2mndt_style', style)
        call ipgeti ('b2mndr_ismain', ismain)
        call ipgeti ('b2sigp_style', ue_style)
        call ipgeti ('ids_from_43', ids_from_43)
        call ipgeti ('b2mndr_eirene', use_eirene)
        call ipgeti ('b2mwti_target_offset', target_offset)
        call ipgetr ('b2mndr_dtim', dtim)
        call ipgetr ('b2stbc_ndes', ndes)
        call ipgetr ('b2stbc_ndes_sol', ndes_sol)
        call ipgeti ('b2stbc_pfrregno1', pfrregno1)
        call ipgeti ('b2stbc_pfrregno2', pfrregno2)
        call ipgetr ('b2stbc_nesepm_pfr', nesepm_pfr)
        call ipgetr ('b2stbc_nesepm_sol', nesepm_sol)
        call ipgetr ('b2stbc_nepedm_sol', nepedm_sol)
        call ipgetr ('b2stbc_volrec_sol', volrec_sol)
        call ipgetr ('b2stbc_private_flux_puff', private_flux_puff)
        call ipgetr ('b2mndr_rescale_neutrals_sources', neutral_sources_rescale)
        call ipgeti ('balance_netcdf', balance_netcdf)
        call ipgetr ('b2agfs_pit_rescale', pit_rescale)
        if (pit_rescale.eq.1.0_R8) then
          filename='b2ag.dat'
          call find_file(filename,exists)
          if (exists) then
            open(99,file=filename)
            call b2agx0 (99, idum(0), idum(1), idum(2), idum(3))
            read (99,'(a8)',err=2) id
            read (99,*,err=2) parg
    1       continue
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
    2       continue
            close(99)
            call ipgetr ('b2agfs_pit_rescale', pit_rescale)
          end if
        end if
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
        if (redef_gmtry.eq.0) then
         drift_style = 1
        else
         drift_style = 2
        end if
        call ipgeti ('b2tfnb_drift_style', drift_style)
        call date_and_time (date, ctime, zone, tvalues)
        username = usrnam()
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
        call alloc_b2mod_user(nx, ns, nlim, nmol, ntns)
        call b2xzef (nx, ny, ns, rz2, na, ne, zeff)
        call b2xppz (nx, ny, ns, ne, na, te, ti, pz)
        geometryType = geometryId(nnreg, isymm, periodic_bc, topcut)

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
            &        bb, conn, vol, gs, hx, hy, hz, qz, qc,                       &
            &        pbs, crx, cry, bzb, lnlam,                                   &
            &        fch, na, ua, te, ti, po, ne, ni, ne2, chvemx, chvimx,        &
            &        cdna, cdpa, cddi, cvla, cvsa, chce, chve, chci,              &
            &        chvi, csig, csigin, calf, cthe, cthi,                        &
            &        cdnahz, cdpahz, cvlahz, cvmahz, cvsahz, cvsa_cl, cvsahz_cl,  &
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
            &        vol, hx, hy, hz, qz, qc, gs, pbs, bb, lnlam,                &
            &        na, ua,                                                     &
            &        uadia, vedia, vadia, wadia, veecrb, vaecrb, ve, wedia,      &
            &        te, ti, po, ne, ni, kinrgy, floe_noc, floi_noc,             &
            &        fna, fna_32, fna_52, fni_32, fni_52, fne_32, fne_52,        &
            &        fna_mdf, fhe_mdf, fhi_mdf, fna_fcor, fna_nodrift, fna_he,   &
            &        fhe, fhi, fhm, fht, fnaPSch, fhePSch, fhiPSch, fch,         &
            &        fchanml_a, fchinert_a, fchvispar_a, fchvisper_a, fchvisq_a, &      !srv 08.09.21
            &        fchdia, fchin, fch_p, fchvispar, fchvisper, fchvisq,        &
            &        fchinert, fchanml, fna_eir, fne_eir, fhe_eir, fhi_eir,      &
            &        cdna, cdpa, cvsa_cl, cvla, chce, chve, chci, chvi, calf,    &
            &        rlsa, rlra, rlqa, rlcx, rlrd, rlbr,                         &
            &        rlza, rlz2, rlpt, rlpi,                                     &
            &        rza, rz2, rpt, rpi, sna, smo, smq, she, shi, sch, sne,      &
            &        wrong_flow, .false.)
#ifdef B25_EIRENE
        if (use_eirene.ne.0) then
          filename=fort_lc//'46'
          call find_file(filename,exists)
          if(exists.and.use_eirene.ne.0) then
            open(unit=46,file=filename)
            call ntread
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
        time_slice_value = time
        time_step = IDS_REAL_INVALID
        num_time_slices = 1
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
        create_date = date//' '//ctime//' '//' '//zone
        !! 1. Set homogeneous_time to 0 or 1 and other properties
        call write_ids_properties( edge_profiles%ids_properties, &
          &  homogeneous_time, comment, source, create_date )
        call write_ids_properties( edge_transport%ids_properties, &
          &  homogeneous_time, comment, source, create_date )
        call write_ids_properties( edge_sources%ids_properties, &
          &  homogeneous_time, comment, source, create_date )
        call write_ids_properties( radiation%ids_properties, &
          &  homogeneous_time, comment, AM_label, create_date )
        call write_ids_properties( description%ids_properties, &
          &  homogeneous_time, comment, source, create_date )
#if IMAS_MINOR_VERSION > 21
        call write_ids_properties( summary%ids_properties, &
          &  homogeneous_time, comment, source, create_date )
#endif
#if ( IMAS_MINOR_VERSION > 25 && IMAS_MINOR_VERSION < 34 )
        call write_ids_properties( numerics%ids_properties, &
          &  homogeneous_time, comment, source, create_date )
#endif
#if IMAS_MINOR_VERSION > 30
        call write_ids_properties( divertors%ids_properties, &
          &  homogeneous_time, comment, source, create_date )
#endif

        !! 2. Set code and library data
        B25_git_version = get_B25_hash()
        code_commit = B25_git_version
        if (streql(b2frates_flag,'adas')) then
          ADAS_git_version = get_ADAS_hash()
          radiation_commit = 'B25 : '//trim(B25_git_version)// &
                      &  ' + ADAS : '//trim(ADAS_git_version)
          p = index(ADAS_git_version,'-')
          if (p.eq.0) then
            adas_version = trim(ADAS_git_version)
          else if (p.gt.1) then
            adas_version = ADAS_git_version(1:p-1)
          else
            adas_version = ''
          end if
        else
          radiation_commit = B25_git_version
        endif
        call write_ids_code( edge_profiles%code, code_commit )
        call write_ids_code( edge_transport%code, code_commit )
        call write_ids_code( edge_sources%code, code_commit )
        call write_ids_code( radiation%code, radiation_commit )
#if IMAS_MINOR_VERSION > 21
        call write_ids_code( summary%code, code_commit )
#endif
#if IMAS_MINOR_VERSION > 30
        call write_ids_code( divertors%code, code_commit )
#endif
        allocate( edge_transport%model(1) )
        edge_transport%model(1)%identifier%index = 1
        allocate( edge_transport%model(1)%identifier%name(1) )
        allocate( edge_transport%model(1)%identifier%description(1) )
        if (ids_from_43.eq.0) then
          if (style.eq.0) then
            edge_transport%model(1)%identifier%name(1) = "SOLPS5.0"
            edge_transport%model(1)%identifier%description(1) = "SOLPS5.0 physics model"
          else if (style.ge.1) then
            edge_transport%model(1)%identifier%name(1) = "SOLPS5.2"
            edge_transport%model(1)%identifier%description(1) = "SOLPS5.2 physics model"
          else if (style.eq.-1) then
            edge_transport%model(1)%identifier%name(1) = "SOLPS4.3"
            edge_transport%model(1)%identifier%description(1) = "SOLPS4.3 physics model"
          end if
          edge_transport%model(1)%flux_multiplier = 1.5_IDS_real + BoRiS
        else
          edge_transport%model(1)%identifier%name(1) = "SOLPS4.3"
          edge_transport%model(1)%identifier%description(1) = "SOLPS4.3 physics model"
          edge_transport%model(1)%flux_multiplier = 2.5_IDS_real
        end if
#if IMAS_MINOR_VERSION > 29
        allocate( edge_transport%model(1)%code%name(1) )
        edge_transport%model(1)%code%name = source
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
        allocate( description%time(num_time_slices) )
        description%time(time_sind) = time
#if IMAS_MINOR_VERSION > 21
        allocate( summary%time(num_time_slices) )
        summary%time(time_sind) = time
#endif
#if ( IMAS_MINOR_VERSION > 25 && IMAS_MINOR_VERSION < 34 )
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
#if IMAS_MINOR_VERSION > 30
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
            allocate( edge_profiles%ggd( num_time_slices ) )
#if IMAS_MINOR_VERSION > 14
            allocate( edge_profiles%grid_ggd( num_time_slices ) )
            allocate( edge_transport%grid_ggd( num_time_slices ) )
            allocate( edge_sources%grid_ggd( num_time_slices ) )
#if IMAS_MINOR_VERSION > 21
            allocate( radiation%grid_ggd( num_time_slices ) )
#endif
#endif
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

#if IMAS_MINOR_VERSION > 21
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
        description%simulation%time_step = time_step_IN
        description%simulation%time_current = time_IN
        allocate( description%simulation%workflow(1) )
        description%simulation%workflow = source
#if ( IMAS_MINOR_VERSION > 25 && IMAS_MINOR_VERSION < 34 )
        description%simulation%time_begin = run_start_time_IN
        description%simulation%time_end = run_end_time_IN
#endif

        i=index(B25_git_version,'-')
        allocate( summary%tag%name(1) )
        summary%tag%name = B25_git_version(1:i-1)
        eq_found = .false.
        eq_source = "ITER Baseline q95=3 equilibrium"
        if ( associated( equilibrium%time_slice ) ) then
          if ( size( equilibrium%time_slice ).ge.time_sind ) then
            eq_found = .true.
            if ( associated( equilibrium%ids_properties%source ) ) &
               & eq_source = equilibrium%ids_properties%source(1)
            r0 = equilibrium%vacuum_toroidal_field%r0
            b0 = equilibrium%vacuum_toroidal_field%b0( time_sind )
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
          end if
        end if
        !> Careful: Sign convention for magnetic field in IDS
        !>          is OPPOSITE to that in SOLPS toroidal geometries
        if ( b0.ne.0.0_IDS_real ) then
          if (streql(database,'ITER').and..not.eq_found) then
            b0r0_ref = 5.3_IDS_real * 6.2_IDS_real
            allocate( edge_profiles%vacuum_toroidal_field%b0( num_time_slices ) )
            if ( pit_rescale.eq.1.0_IDS_real ) then
              i = nint(b0r0_ref/b0r0)
              select case (i)
              case (1)
                call write_sourced_value( summary%global_quantities%ip, -15.0e6_IDS_real )
                call write_sourced_value( summary%global_quantities%b0, -5.3_IDS_real )
                edge_profiles%vacuum_toroidal_field%b0( time_sind ) = -5.3_IDS_real
              case (2)
                call write_sourced_value( summary%global_quantities%ip, -7.5e6_IDS_real )
                call write_sourced_value( summary%global_quantities%b0, -2.65_IDS_real )
                edge_profiles%vacuum_toroidal_field%b0( time_sind ) = -2.65_IDS_real
              case (3)
                call write_sourced_value( summary%global_quantities%ip, -5.0e6_IDS_real )
                call write_sourced_value( summary%global_quantities%b0, -1.8_IDS_real )
                edge_profiles%vacuum_toroidal_field%b0( time_sind ) = -1.8_IDS_real
              case default
                call write_sourced_value( summary%global_quantities%ip, &
                  &  -15.0e6_IDS_real/nint(b0r0_ref/b0r0) )
                call write_sourced_value( summary%global_quantities%b0, &
                  &  -b0r0 / 6.2e6_IDS_real )
                edge_profiles%vacuum_toroidal_field%b0( time_sind ) = -b0r0 / 6.2_IDS_real
              end select
              summary%global_quantities%ip%source = "ITER Baseline q95=3 equilibrium"
              call write_sourced_value( summary%global_quantities%q_95, 3.0_IDS_real )
              summary%global_quantities%q_95%source = "ITER Baseline q95=3 equilibrium"
            else
              call write_sourced_value( summary%global_quantities%b0, -b0 )
              edge_profiles%vacuum_toroidal_field%b0( time_sind ) = -b0
              write(eq_source, '(a,1pe12.5,a)' ) "ITER Baseline q95=3 equilibrium"// &
                  &  " (current rescaled by ",pit_rescale,")"
              call write_sourced_value( summary%global_quantities%ip, &
                  &  -15.0e6_IDS_real*pit_rescale )
              summary%global_quantities%ip%source = eq_source
              call write_sourced_value( summary%global_quantities%q_95, &
                  &   3.0_IDS_real/pit_rescale )
              summary%global_quantities%q_95%source = eq_source
            end if
            call write_sourced_constant( summary%global_quantities%r0, 6.2_IDS_real )
            edge_profiles%vacuum_toroidal_field%r0 = 6.2_IDS_real
          else
            call write_sourced_constant( summary%global_quantities%r0, r0 )
            edge_profiles%vacuum_toroidal_field%r0 = r0
            allocate( edge_profiles%vacuum_toroidal_field%b0( num_time_slices ) )
            if (eq_found) then
              call write_sourced_value( summary%global_quantities%b0, b0 )
              edge_profiles%vacuum_toroidal_field%b0( time_sind ) = b0
              if ( equilibrium%time_slice( time_sind )%global_quantities%  &
                 & ip .ne. IDS_REAL_INVALID ) then
                call write_sourced_value( summary%global_quantities%ip,    &
                  &  equilibrium%time_slice( time_sind )%global_quantities%ip )
                summary%global_quantities%ip%source = eq_source
              end if
              if ( equilibrium%time_slice( time_sind )%global_quantities%  &
                 & q_95 .ne. IDS_REAL_INVALID ) then
                call write_sourced_value( summary%global_quantities%q_95,  &
                  &  equilibrium%time_slice( time_sind )%global_quantities%q_95 )
                summary%global_quantities%q_95%source = eq_source
              else if (streql(eq_source,"ITER Baseline q95=3 equilibrium")) then
                call write_sourced_value( summary%global_quantities%q_95,  &
                  &  3.0_IDS_real )
                summary%global_quantities%q_95%source = eq_source
              end if
            else if (isymm.ne.0) then
              call write_sourced_value( summary%global_quantities%b0, -b0 )
              edge_profiles%vacuum_toroidal_field%b0( time_sind ) = -b0
            else
              call write_sourced_value( summary%global_quantities%b0, -b0 )
              edge_profiles%vacuum_toroidal_field%b0( time_sind ) = b0
            end if
          end if
        end if

        if (GeometryType .eq. GEOMETRY_LINEAR) then
          midplane_id = 4
        else
          if ( eq_found ) then
            z_eq = equilibrium%time_slice( time_sind )%global_quantities%  &
               &   magnetic_axis%z
          else
            z_eq = IDS_REAL_INVALID
          end if
          if ( z_eq.ne.IDS_REAL_INVALID .and. &
             & (cry(jxa,jsep,2)-z_eq)*(cry(jxa,jsep,3)-z_eq).lt.0.0_R8 ) then
            midplane_id = 1
          else if ( jxa .eq. nmdpl ) then
            midplane_id = 2
          else if ( cry(jxa,jsep,2)*cry(jxa,jsep,3).lt.0.0_R8 ) then
            midplane_id = 3
          else
            midplane_id = 4
          end if
        end if
#if IMAS_MINOR_VERSION > 32
        call write_ids_midplane( divertors%midplane, midplane_id )
        call write_ids_midplane( edge_profiles%midplane, midplane_id )
        call write_ids_midplane( edge_sources%midplane, midplane_id )
        call write_ids_midplane( edge_transport%midplane, midplane_id )
        call write_ids_midplane( summary%midplane, midplane_id )
#endif

        select case (GeometryType)
        case ( GEOMETRY_CYLINDER, GEOMETRY_LIMITER, GEOMETRY_ANNULUS )
          icnt = 1
          isep(1) = 2
        case ( GEOMETRY_SN, GEOMETRY_STELLARATORISLAND )
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
              if (region(ix,iy,2).eq.isep(i)) u = u + fht(ix,iy,1)
            end do
          end do
        end do
        if (u.ne.0.0_IDS_real) then
#if IMAS_MINOR_VERSION > 28
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
#if IMAS_MINOR_VERSION > 30
          if (frac.ne.0.0_IDS_real) then
            call write_sourced_value( summary%global_quantities%power_radiated_inside_lcfs, frac )
            call write_sourced_value( summary%global_quantities%power_radiated_outside_lcfs, &
              &  u - frac )
          end if
#endif
        case ( GEOMETRY_LIMITER, GEOMETRY_SN, &
            &  GEOMETRY_STELLARATORISLAND, GEOMETRY_ANNULUS , &
            &  GEOMETRY_CDN, GEOMETRY_DDN_BOTTOM, GEOMETRY_DDN_TOP )
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
#if IMAS_MINOR_VERSION > 30
          if (u.ne.0.0_IDS_real) then
            call write_sourced_value( summary%global_quantities%power_radiated_outside_lcfs, u )
          end if
#endif
        end select

! Determine divertor plate generic information
        if (nncut.eq.0) then
          if (geometryType.eq.GEOMETRY_LINEAR .or. &
            & geometryType.eq.GEOMETRY_CYLINDER) then
            ntrgts=0
            if (boundary_namelist.ne.0) then
              target_east = .false.
              target_west = .false.
              do i=1,nbc
                if(bcchar(i).eq.'E'.and.bcpos(i).eq.-1) then
                  target_east = bcene(i).eq. 3.or. &
                        &       bcene(i).eq.12.or.bcene(i).eq.15
                end if
                if(bcchar(i).eq.'W'.and.bcpos(i).eq.nx) then
                  target_west = bcene(i).eq. 3.or. &
                        &       bcene(i).eq.12.or.bcene(i).eq.15
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
                flux_expansion(ntrgts) =                                         &
                 & ( wbbv(topix(jxa,jsep),topiy(jxa,jsep),0)/                    &
                 &   wbbv(topix(jxa,jsep),topiy(jxa,jsep),3) )/                  &
                 & ( wbbc(rightix(topix(bcpos(i),jsep),topiy(bcpos(i),jsep)),    &
                 &        rightiy(topix(bcpos(i),jsep),topiy(bcpos(i),jsep)),0)/ &
                 &   wbbc(rightix(topix(bcpos(i),jsep),topiy(bcpos(i),jsep)),    &
                 &        rightiy(topix(bcpos(i),jsep),topiy(bcpos(i),jsep)),3) )
                r_max = max(maxval(crx(bc_list_x(1:bc_list_size(i),i),           &
                 &                     bc_list_y(1:bc_list_size(i),i),1)),       &
                 &          maxval(crx(bc_list_x(1:bc_list_size(i),i),           &
                 &                     bc_list_y(1:bc_list_size(i),i),3)))
                r_min = min(minval(crx(bc_list_x(1:bc_list_size(i),i),           &
                 &                     bc_list_y(1:bc_list_size(i),i),1)),       &
                 &          minval(crx(bc_list_x(1:bc_list_size(i),i),           &
                 &                     bc_list_y(1:bc_list_size(i),i),3)))
                z_max = max(maxval(cry(bc_list_x(1:bc_list_size(i),i),           &
                 &                     bc_list_y(1:bc_list_size(i),i),1)),       &
                 &          maxval(cry(bc_list_x(1:bc_list_size(i),i),           &
                 &                     bc_list_y(1:bc_list_size(i),i),3)))
                z_min = min(minval(cry(bc_list_x(1:bc_list_size(i),i),           &
                 &                     bc_list_y(1:bc_list_size(i),i),1)),       &
                 &          minval(cry(bc_list_x(1:bc_list_size(i),i),           &
                 &                     bc_list_y(1:bc_list_size(i),i),3)))
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
                 & ( wbbv(topix(jxa,jsep),topiy(jxa,jsep),0)/                    &
                 &   wbbv(topix(jxa,jsep),topiy(jxa,jsep),3) )/                  &
                 & ( wbbc(topix(bcpos(i),jsep),topiy(bcpos(i),jsep),0)/          &
                 &   wbbc(topix(bcpos(i),jsep),topiy(bcpos(i),jsep),3) )
                r_max = max(maxval(crx(bc_list_x(1:bc_list_size(i),i),           &
                 &                     bc_list_y(1:bc_list_size(i),i),0)),       &
                 &          maxval(crx(bc_list_x(1:bc_list_size(i),i),           &
                 &                     bc_list_y(1:bc_list_size(i),i),2)))
                r_min = min(minval(crx(bc_list_x(1:bc_list_size(i),i),           &
                 &                     bc_list_y(1:bc_list_size(i),i),0)),       &
                 &          minval(crx(bc_list_x(1:bc_list_size(i),i),           &
                 &                     bc_list_y(1:bc_list_size(i),i),2)))
                z_max = max(maxval(cry(bc_list_x(1:bc_list_size(i),i),           &
                 &                     bc_list_y(1:bc_list_size(i),i),0)),       &
                 &          maxval(cry(bc_list_x(1:bc_list_size(i),i),           &
                 &                     bc_list_y(1:bc_list_size(i),i),2)))
                z_min = min(minval(cry(bc_list_x(1:bc_list_size(i),i),           &
                 &                     bc_list_y(1:bc_list_size(i),i),0)),       &
                 &          minval(cry(bc_list_x(1:bc_list_size(i),i),           &
                 &                     bc_list_y(1:bc_list_size(i),i),2)))
                extension_r(ntrgts) = r_max - r_min
                extension_z(ntrgts) = z_max - z_min
              end if
            end if
          else
            ntrgts = 0
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
              &                   wbbv(topix(jxa,jsep),topiy(jxa,jsep),3) )/ &
              &                 ( wbbc(topix(nx,jsep),topiy(nx,jsep),0)/     &
              &                   wbbc(topix(nx,jsep),topiy(nx,jsep),3) )
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
        totFace=abs(fht)
        call divide_by_poloidal_areas(nx,ny,totFace,tmpFace)
        call divide_by_contact_areas(nx,ny,totFace,flxFace)
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
              qtot = qtot + idir(itrg(i))*fht(ix,iy,0)
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
                 & (1.5_R8*fna_32(ifpos(itrg(i)),iy,0,is) + &
                 &  2.5_R8*fna_52(ifpos(itrg(i)),iy,0,is))
              if (is_neutral(is)) &
                &  power_neutrals(itrg(i)) = power_neutrals(itrg(i)) + &
                &    idir(itrg(i))*(ti(ixpos(itrg(i)),iy)* &
                &   (1.5_R8*fna_32(ifpos(itrg(i)),iy,0,is) +  &
                &    2.5_R8*fna_52(ifpos(itrg(i)),iy,0,is)) + &
                &              fhm(ifpos(itrg(i)),iy,0,is))
              match_found = .true.
              do j = 1, nstrai
                if (streql(crcstra(j),'E').or.streql(crcstra(j),'W')) then
                  if (rcpos(j).eq. &
                     &  (ixpos(itrg(i))+target_offset*idir(itrg(i))).and. &
                     &  (rcstart(j).le.iy .and. rcend(j).ge.iy)) then
                    recycled_flux(itrg(i)) = recycled_flux(itrg(i)) + &
                      &  neutral_sources_rescale*recyc(is,j)* &
                      &  max(0.0_R8,idir(itrg(i))* &
                      &  fna(ifpos(itrg(i)),iy,0,is))*zn(is)
                  end if
                end if
              end do
            end do
            power_radiated(itrg(i)) = power_radiated(itrg(i)) + &
                 &  wrdtrg(iy,itrg(i),0)
            power_convected(itrg(i)) = power_convected(itrg(i)) + u
            power_conducted(itrg(i)) = power_conducted(itrg(i))   &
                 & - u + idir(itrg(i))* &
                 & (fht(ifpos(itrg(i)),iy,0)-fhj(ifpos(itrg(i)),iy,0))
            power_incident(itrg(i)) = power_incident(itrg(i)) + &
                 &  wrdtrg(iy,itrg(i),0) + &
                 &  idir(itrg(i))*fht(ifpos(itrg(i)),iy,0)
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
              flxFace(ifpos(itrg(i)),iy,0) = flxFace(ifpos(itrg(i)),iy,0) + &
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
          power_flux_peak(itrg(i)) = maxval(flxFace(ifpos(itrg(i)),0:ny-1,0))
        end do
#if IMAS_MINOR_VERSION > 30
        select case ( GeometryType )
        case ( GEOMETRY_LINEAR, GEOMETRY_CYLINDER )
          if (ntrgts.gt.0) then
            allocate( divertors%divertor(ntrgts) )
            do i = 1, ntrgts
              allocate( divertors%divertor(i)%name(1) )
              allocate( divertors%divertor(i)%identifier(1) )
              allocate( divertors%divertor(i)%target(1) )
              allocate( divertors%divertor(i)%target(1)%name(1) )
              allocate( divertors%divertor(i)%target(1)%identifier(1) )
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
              divertors%divertor(i)%identifier = plate_name(i)
              divertors%divertor(i)%target(1)%identifier = plate_name(i)
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
            u = idir(1)*sum(fhp(ifpos(1),:,0,:))
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
            u = idir(1)*sum(fhj(ifpos(1),:,0))
            call write_timed_value( &
              &  divertors%divertor(1)%target(1)%power_currents, u )
            call write_timed_value( &
              &  divertors%divertor(1)%power_currents, u )
#if IMAS_MINOR_VERSION > 32
            u = idir(1)*sum(fch(ifpos(1),:,0))
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
                & divertors%divertor(2)%target(1)%power_flux_peak, &
                & power_flux_peak(2) )
              call write_timed_value( &
                & divertors%divertor(2)%target(1)%flux_expansion, &
                & flux_expansion(2) )
              call write_timed_value( &
                & divertors%divertor(2)%target(1)%wetted_area, &
                & wetted_area(2) )
              call write_timed_value( &
                & divertors%divertor(2)%wetted_area, &
                & wetted_area(2) )
              call write_timed_value( &
                & divertors%divertor(2)%target(1)%power_incident_fraction, &
                & 1.0_IDS_real )
              call write_timed_value( &
                & divertors%divertor(2)%target(1)%power_incident, &
                & power_incident(2) )
              call write_timed_value( &
                & divertors%divertor(2)%power_incident, &
                & power_incident(2) )
              call write_timed_value( &
                & divertors%divertor(2)%target(1)%power_conducted, &
                & power_conducted(2) )
              call write_timed_value( &
                & divertors%divertor(2)%power_conducted, &
                & power_conducted(2) )
              call write_timed_value( &
                & divertors%divertor(2)%target(1)%power_convected, &
                & power_convected(2) )
              call write_timed_value( &
                & divertors%divertor(2)%power_convected, &
                & power_convected(2) )
              call write_timed_value( &
                & divertors%divertor(2)%target(1)%power_radiated, &
                & power_radiated(2) )
              call write_timed_value( &
                & divertors%divertor(2)%power_radiated, &
                & power_radiated(2) )
              call write_timed_value( &
                & divertors%divertor(2)%target(1)%power_neutrals, &
                & power_neutrals(2) )
              call write_timed_value( &
                & divertors%divertor(2)%power_neutrals, &
                & power_neutrals(2) )
              u = idir(2)*sum(fhp(ifpos(2),:,0,:))
              call write_timed_value( &
                & divertors%divertor(2)%target(1)%power_recombination_plasma, u )
              call write_timed_value( &
                & divertors%divertor(2)%power_recombination_plasma, u )
              call write_timed_value( &
                & divertors%divertor(2)%target(1)%power_recombination_neutrals, &
                & power_recomb_neutrals(2) )
              call write_timed_value( &
                & divertors%divertor(2)%power_recombination_neutrals, &
                & power_recomb_neutrals(2) )
              u = idir(2)*sum(fhj(ifpos(2),:,0))
              call write_timed_value( &
                & divertors%divertor(2)%target(1)%power_currents, u )
              call write_timed_value( &
                & divertors%divertor(2)%power_currents, u )
#if IMAS_MINOR_VERSION > 32
              u = idir(2)*sum(fch(ifpos(2),:,0))
              call write_timed_value( &
                & divertors%divertor(2)%target(1)%current_incident, u )
              call write_timed_value( &
                & divertors%divertor(2)%current_incident, u )
#endif
              call write_timed_value( &
                & divertors%divertor(2)%particle_flux_recycled_total, &
                & recycled_flux(2) )
            end if
          end if
        case ( GEOMETRY_SN, GEOMETRY_STELLARATORISLAND )
          allocate( divertors%divertor(1) )
          allocate( divertors%divertor(1)%name(1) )
          allocate( divertors%divertor(1)%identifier(1) )
          if (LSN) then
            divertors%divertor(1)%name = 'Lower divertor'
            divertors%divertor(1)%identifier = 'LSN'
          else
            divertors%divertor(1)%name = 'Upper divertor'
            divertors%divertor(1)%identifier = 'USN'
          end if
          allocate( divertors%divertor(1)%target(2) )
          allocate( divertors%divertor(1)%target(1)%name(1) )
          allocate( divertors%divertor(1)%target(1)%identifier(1) )
          allocate( divertors%divertor(1)%target(2)%name(1) )
          allocate( divertors%divertor(1)%target(2)%identifier(1) )
!! FIXME: Should represent the full extent of the physical divertor
          divertors%divertor(1)%target(1)%extension_r = extension_r(1)
          divertors%divertor(1)%target(1)%extension_z = extension_z(1)
          divertors%divertor(1)%target(2)%extension_r = extension_r(2)
          divertors%divertor(1)%target(2)%extension_z = extension_z(2)
          if (LSN) then
            divertors%divertor(1)%target(1)%name = "Inner target"
            divertors%divertor(1)%target(1)%identifier = "ID"
            divertors%divertor(1)%target(2)%name = "Outer target"
            divertors%divertor(1)%target(2)%identifier = "OD"
          else
            divertors%divertor(1)%target(1)%name = "Outer target"
            divertors%divertor(1)%target(1)%identifier = "OD"
            divertors%divertor(1)%target(2)%name = "Inner target"
            divertors%divertor(1)%target(2)%identifier = "ID"
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
          u = idir(1)*sum(fhp(ifpos(1),:,0,:))
          call write_timed_value( &
            &  divertors%divertor(1)%target(1)%power_recombination_plasma, u )
          v = idir(2)*sum(fhp(ifpos(2),:,0,:))
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
          u = idir(1)*sum(fhj(ifpos(1),:,0))
          call write_timed_value( &
            &  divertors%divertor(1)%target(1)%power_currents, u )
          v = idir(2)*sum(fhj(ifpos(2),:,0))
          call write_timed_value( &
            &  divertors%divertor(1)%target(2)%power_currents, v )
          call write_timed_value( &
            &  divertors%divertor(1)%power_currents, u+v )
#if IMAS_MINOR_VERSION > 32
          u = idir(1)*sum(fch(ifpos(1),:,0))
          call write_timed_value( &
            &  divertors%divertor(1)%target(1)%current_incident, u )
          v = idir(2)*sum(fch(ifpos(2),:,0))
          call write_timed_value( &
            &  divertors%divertor(1)%target(2)%current_incident, v )
          call write_timed_value( &
            &  divertors%divertor(1)%current_incident, u+v )
#endif
          call write_timed_value( &
            &  divertors%divertor(1)%particle_flux_recycled_total, &
            &  recycled_flux(1) )
        case ( GEOMETRY_CDN, GEOMETRY_DDN_BOTTOM, GEOMETRY_DDN_TOP )
          allocate( divertors%divertor(2) )
          allocate( divertors%divertor(1)%name(1) )
          allocate( divertors%divertor(1)%identifier(1) )
          allocate( divertors%divertor(2)%name(1) )
          allocate( divertors%divertor(2)%identifier(1) )
          divertors%divertor(1)%name = 'Lower divertor'
          divertors%divertor(1)%identifier = 'LD'
          divertors%divertor(2)%name = 'Upper divertor'
          divertors%divertor(2)%identifier = 'UD'
          allocate( divertors%divertor(1)%target(2) )
          allocate( divertors%divertor(2)%target(2) )
          allocate( divertors%divertor(1)%target(1)%name(1) )
          allocate( divertors%divertor(1)%target(1)%identifier(1) )
          allocate( divertors%divertor(1)%target(2)%name(1) )
          allocate( divertors%divertor(1)%target(2)%identifier(1) )
          allocate( divertors%divertor(2)%target(1)%name(1) )
          allocate( divertors%divertor(2)%target(1)%identifier(1) )
          allocate( divertors%divertor(2)%target(2)%name(1) )
          allocate( divertors%divertor(2)%target(2)%identifier(1) )
!! FIXME: Should represent the full extent of the physical divertor
          divertors%divertor(1)%target(1)%extension_r = extension_r(1)
          divertors%divertor(1)%target(1)%extension_z = extension_z(1)
          divertors%divertor(1)%target(2)%extension_r = extension_r(4)
          divertors%divertor(1)%target(2)%extension_z = extension_z(4)
          divertors%divertor(2)%target(1)%extension_r = extension_r(2)
          divertors%divertor(2)%target(1)%extension_z = extension_z(2)
          divertors%divertor(2)%target(2)%extension_r = extension_r(3)
          divertors%divertor(2)%target(2)%extension_z = extension_z(3)
          divertors%divertor(1)%target(1)%name = "Lower inner target"
          divertors%divertor(1)%target(1)%identifier = "LID"
          divertors%divertor(1)%target(2)%name = "Lower outer target"
          divertors%divertor(1)%target(2)%identifier = "LOD"
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
          u = idir(1)*sum(fhp(ifpos(1),:,0,:))
          call write_timed_value( &
            &  divertors%divertor(1)%target(1)%power_recombination_plasma, u )
          v = idir(4)*sum(fhp(ifpos(4),:,0,:))
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
          u = idir(1)*sum(fhj(ifpos(1),:,0))
          call write_timed_value( &
            &  divertors%divertor(1)%target(1)%power_currents, u )
          v = idir(4)*sum(fhj(ifpos(4),:,0))
          call write_timed_value( &
            &  divertors%divertor(1)%target(2)%power_currents, v )
          call write_timed_value( &
            &  divertors%divertor(1)%power_currents, u+v )
#if IMAS_MINOR_VERSION > 32
          u = idir(1)*sum(fch(ifpos(1),:,0))
          call write_timed_value( &
            &  divertors%divertor(1)%target(1)%current_incident, u )
          v = idir(4)*sum(fch(ifpos(4),:,0))
          call write_timed_value( &
            &  divertors%divertor(1)%target(2)%current_incident, v )
          call write_timed_value( &
            &  divertors%divertor(1)%current_incident, u+v )
#endif
          call write_timed_value( &
            &  divertors%divertor(1)%particle_flux_recycled_total, &
            &  recycled_flux(1)+recycled_flux(4) )
          divertors%divertor(2)%target(1)%name = "Upper inner target"
          divertors%divertor(2)%target(1)%identifier = "UID"
          divertors%divertor(2)%target(2)%name = "Upper outer target"
          divertors%divertor(2)%target(2)%identifier = "UOD"
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
          u = idir(2)*sum(fhp(ifpos(2),:,0,:))
          call write_timed_value( &
            &  divertors%divertor(2)%target(1)%power_recombination_plasma, u )
          v = idir(3)*sum(fhp(ifpos(3),:,0,:))
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
          u = idir(2)*sum(fhj(ifpos(2),:,0))
          call write_timed_value( &
            &  divertors%divertor(2)%target(1)%power_currents, u )
          v = idir(3)*sum(fhj(ifpos(3),:,0))
          call write_timed_value( &
            &  divertors%divertor(2)%target(2)%power_currents, v )
          call write_timed_value( &
            &  divertors%divertor(2)%power_currents, u+v )
#if IMAS_MINOR_VERSION > 32
          u = idir(2)*sum(fch(ifpos(2),:,0))
          call write_timed_value( &
            &  divertors%divertor(2)%target(1)%current_incident, u )
          v = idir(3)*sum(fch(ifpos(3),:,0))
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
            call b2_IMAS_Fill_Grid_Desc( gmap,                                &
                &   edge_sources%source(is)%ggd( time_sind )%grid,            &
                &   nx, ny, crx(-1:nx, -1:ny, :), cry(-1:nx, -1:ny, : ),      &
                &   leftix, leftiy, rightix, rightiy, topix, topiy, bottomix, &
                &   bottomiy, nnreg, topcut, region, cflags,                  &
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
        allocate( edge_profiles%ggd( time_sind )%ion( nsion ) )
#if IMAS_MINOR_VERSION > 21
        allocate( radiation%process(1)%ggd( time_sind )%ion( nsion ) )
        allocate( radiation%process(2)%ggd( time_sind )%ion( nsion ) )
#endif
        do i = 1, nsources
          allocate( edge_sources%source(i)%ggd( time_sind )%ion( nsion ) )
        end do
        allocate( edge_transport%model(1)%ggd( time_sind )%ion( nsion ) )
        do js = 1, nspecies
          allocate( edge_profiles%ggd( time_sind )%ion( js )%label(1) )
          allocate( edge_profiles%ggd( time_sind )%ion( js )%state( nfluids(js) ) )
          do is = 1, nfluids(js)
            allocate( edge_profiles%ggd( time_sind )%ion( js )%state( is )%label(1) )
          end do
          allocate( edge_profiles%ggd( time_sind )%ion( js )%element(1) )
          do i = 1, nsources
            allocate( edge_sources%source(i)%ggd( time_sind )%ion( js )%label(1) )
            allocate( edge_sources%source(i)%ggd( time_sind )%ion( js )%state( nfluids(js) ) )
            do is = 1, nfluids(js)
              allocate( edge_sources%source(i)%ggd( time_sind )%ion( js )%state( is )%label(1) )
            end do
            allocate( edge_sources%source(i)%ggd( time_sind )%ion( js )%element(1) )
          end do
          allocate( edge_transport%model(1)%ggd( time_sind )%ion( js )%label(1) )
          allocate( edge_transport%model(1)%ggd( time_sind )%ion( js )%state( nfluids(js) ) )
          do is = 1, nfluids(js)
            allocate( edge_transport%model(1)%ggd( time_sind )%ion( js )%state( is )%label(1) )
          end do
          allocate( edge_transport%model(1)%ggd( time_sind )%ion( js )%element(1) )

#if IMAS_MINOR_VERSION > 21
          allocate( radiation%process(1)%ggd( time_sind )%ion( js )%label(1) )
          allocate( radiation%process(1)%ggd( time_sind )%ion( js )%state( nfluids(js) ) )
          do is = 1, nfluids(js)
            allocate( radiation%process(1)%ggd( time_sind )%ion( js )%state( is )%label(1) )
          end do
          allocate( radiation%process(1)%ggd( time_sind )%ion( js )%element(1) )
          allocate( radiation%process(2)%ggd( time_sind )%ion( js )%label(1) )
          allocate( radiation%process(2)%ggd( time_sind )%ion( js )%state( nfluids(js) ) )
          do is = 1, nfluids(js)
            allocate( radiation%process(2)%ggd( time_sind )%ion( js )%state( is )%label(1) )
          end do
          allocate( radiation%process(2)%ggd( time_sind )%ion( js )%element(1) )
#endif
          ! Put label to ion(js).state(is).label
          do is = 1, istion(js)
            ks = ispion(js,is)
            call species( ks, spclabel, .false.)
            edge_profiles%ggd( time_sind )%ion( js )%state( is )%label = spclabel
            do i = 1, nsources
              call species( ks, spclabel, .false.)
              edge_sources%source(i)%ggd( time_sind )%ion( js )%state( is )%label = spclabel
            end do
            call species( ks, spclabel, .false.)
            edge_transport%model(1)%ggd( time_sind )%ion( js )%state( is )%label = spclabel
#if IMAS_MINOR_VERSION > 21
            call species( ks, spclabel, .false.)
            radiation%process(1)%ggd( time_sind )%ion( js )%state( is )%label = spclabel
            call species( ks, spclabel, .false.)
            radiation%process(2)%ggd( time_sind )%ion( js )%state( is )%label = spclabel
#endif
          end do

          ! Put ion label identifying the species
          edge_profiles%ggd( time_sind )%ion( js )%label = species_list( js )
          edge_transport%model(1)%ggd( time_sind )%ion( js )%label = species_list( js )

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
          edge_profiles%ggd( time_sind )%ion( js )%element(1)%z_n = zn( is )
          edge_transport%model(1)%ggd( time_sind )%ion( js )%element(1)%z_n = zn( is )

          ! Put number of atoms
#if IMAS_MINOR_VERSION < 15
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
            edge_sources%source(i)%ggd( time_sind )%ion( js )%label = species_list( js )

            ! Put ion charge if single ion in species
            is = ispion(js,1)
            if (istion(js).eq.1) then
              ion_charge_int = nint((zamin(is)+zamax(is))/2.0_R8)
              edge_sources%source(i)%ggd( time_sind )%ion( js )%z_ion = ion_charge_int
            end if

            ! Put mass of ion
            edge_sources%source(i)%ggd( time_sind )%ion( js )%element(1)%a = am( is )

            ! Put nuclear charge
            edge_sources%source(i)%ggd( time_sind )%ion( js )%element(1)%z_n = zn( is )

            ! Put number of atoms
#if IMAS_MINOR_VERSION < 15
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

#if IMAS_MINOR_VERSION > 21
          ! Put ion label identifying the species
          radiation%process(1)%ggd( time_sind )%ion( js )%label = species_list( js )
          radiation%process(2)%ggd( time_sind )%ion( js )%label = species_list( js )

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
          radiation%process(1)%ggd( time_sind )%ion( js )%element(1)%z_n = zn( is )
          radiation%process(2)%ggd( time_sind )%ion( js )%element(1)%z_n = zn( is )

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

        if (use_eirene.ne.0) then
#ifdef B25_EIRENE
          do js = nspecies+1, nsion
            allocate( edge_profiles%ggd( time_sind )%ion( js )%label(1) )
            allocate( edge_profiles%ggd( time_sind )%ion( js )%state( istion(js) ) )
            do ks = 1, istion(js)
              allocate( edge_profiles%ggd( time_sind )%ion( js )%state( ks )%label(1) )
            end do
            do i = 1, nsources
              allocate( edge_sources%source(i)%ggd( time_sind )%ion( js )%label(1) )
              allocate( edge_sources%source(i)%ggd( time_sind )%ion( js )%state( istion(js) ) )
              do ks = 1, istion(js)
                allocate( edge_sources%source(i)%ggd( time_sind )%ion( js )%state( ks )%label(1) )
              end do
            end do
            allocate( edge_transport%model(1)%ggd( time_sind )%ion( js )%label(1) )
            allocate( edge_transport%model(1)%ggd( time_sind )%ion( js )%state( istion(js) ) )
            do ks = 1, istion(js)
              allocate( edge_transport%model(1)%ggd( time_sind )%ion( js )%state( ks )%label(1) )
            end do
            do ks = 1, istion(js)
              is = ispion(js,ks)
              edge_profiles%ggd( time_sind )%ion( js )%state( ks )%label = textin( is-1 )
              edge_transport%model(1)%ggd( time_sind )%ion( js )%state( ks )%label = &
                  &                                                        textin( is-1 )
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
            end do
            is = ispion(js,1)
            if (istion(js).eq.1) then
              edge_profiles%ggd( time_sind )%ion( js )%z_ion = nchrgi( is )
              edge_profiles%ggd( time_sind )%ion( js )%label = textin( is-1 )
              do i = 1, nsources
                edge_sources%source(i)%ggd( time_sind )%ion( js )%z_ion = nchrgi( is )
                edge_sources%source(i)%ggd( time_sind )%ion( js )%label = textin( is-1 )
              end do
              edge_transport%model(1)%ggd( time_sind )%ion( js )%z_ion = nchrgi( is )
              edge_transport%model(1)%ggd( time_sind )%ion( js )%label = textin( is-1 )
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
                  edge_profiles%ggd( time_sind )%ion( js )%label = ion_label
                  do i = 1, nsources
                    edge_sources%source(i)%ggd( time_sind )%ion( js )%label = ion_label
                  end do
                  edge_transport%model(1)%ggd( time_sind )%ion( js )%label = ion_label
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
               allocate( edge_profiles%ggd( time_sind )%neutral( js )%label(1) )
               edge_profiles%ggd( time_sind )%neutral( js )%element(1)%a = am( is )
               edge_profiles%ggd( time_sind )%neutral( js )%element(1)%z_n = zn( is )
#if IMAS_MINOR_VERSION < 15
               edge_profiles%ggd( time_sind )%neutral( js )%element(1)%multiplicity = 1.0_IDS_real
#else
               edge_profiles%ggd( time_sind )%neutral( js )%element(1)%atoms_n = 1
#endif
               edge_profiles%ggd( time_sind )%neutral( js )%label = species_list( js )
               edge_profiles%ggd( time_sind )%neutral( js )%ion_index = js
               do i = 1, nsources
                  allocate( edge_sources%source(i)%ggd( time_sind )%neutral( js )%element(1) )
                  allocate( edge_sources%source(i)%ggd( time_sind )%neutral( js )%label(1) )
                  edge_sources%source(i)%ggd( time_sind )%neutral( js )%element(1)%a = am( is )
                  edge_sources%source(i)%ggd( time_sind )%neutral( js )%element(1)%z_n = zn( is )
#if IMAS_MINOR_VERSION < 15
                  edge_sources%source(i)%ggd( time_sind )%neutral( js )%element(1)%multiplicity = 1.0_IDS_real
#else
                  edge_sources%source(i)%ggd( time_sind )%neutral( js )%element(1)%atoms_n = 1
#endif
                  edge_sources%source(i)%ggd( time_sind )%neutral( js )%label = species_list( js )
                  edge_sources%source(i)%ggd( time_sind )%neutral( js )%ion_index = js
               end do
               allocate( edge_transport%model(1)%ggd( time_sind )%neutral( js )%element(1) )
               allocate( edge_transport%model(1)%ggd( time_sind )%neutral( js )%label(1) )
               edge_transport%model(1)%ggd( time_sind )%neutral( js )%element(1)%a = am( is )
               edge_transport%model(1)%ggd( time_sind )%neutral( js )%element(1)%z_n = zn( is )
#if IMAS_MINOR_VERSION < 15
               edge_transport%model(1)%ggd( time_sind )%neutral( js )%element(1)%multiplicity = 1.0_IDS_real
#else
               edge_transport%model(1)%ggd( time_sind )%neutral( js )%element(1)%atoms_n = 1
#endif
               edge_transport%model(1)%ggd( time_sind )%neutral( js )%label = species_list( js )
               edge_transport%model(1)%ggd( time_sind )%neutral( js )%ion_index = js
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
                  edge_profiles%ggd( time_sind )%neutral( js )%multiple_states_flag = 1
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
                  edge_transport%model(1)%ggd( time_sind )%neutral( js )%multiple_states_flag = 1
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
                  edge_sources%source(i)%ggd( time_sind )%neutral( js )%ion_index = k
                  edge_sources%source(i)%ggd( time_sind )%neutral( js )%multiple_states_flag = 1
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
               edge_transport%model(1)%ggd( time_sind )%neutral( js )%ion_index = k
               edge_transport%model(1)%ggd( time_sind )%neutral( js )%multiple_states_flag = 1
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
               allocate( edge_profiles%ggd( time_sind )%neutral( j )%element(1) )
               allocate( edge_profiles%ggd( time_sind )%neutral( j )%label(1) )
               allocate( edge_profiles%ggd( time_sind )%neutral( j )%state(1) )
               allocate( edge_profiles%ggd( time_sind )%neutral( j )%state(1)%label(1) )
               edge_profiles%ggd( time_sind )%neutral( j )%element(1)%a = am( is )
               edge_profiles%ggd( time_sind )%neutral( j )%element(1)%z_n = zn( is )
#if IMAS_MINOR_VERSION < 15
               edge_profiles%ggd( time_sind )%neutral( j )%element(1)%multiplicity = 1
#else
               edge_profiles%ggd( time_sind )%neutral( j )%element(1)%atoms_n = 1
#endif
               edge_profiles%ggd( time_sind )%neutral( j )%label = species_list( js )
               edge_profiles%ggd( time_sind )%neutral( j )%ion_index = js
               edge_profiles%ggd( time_sind )%neutral( j )%multiple_states_flag = 1
               call species( is, spclabel, .false. )
               edge_profiles%ggd( time_sind )%neutral( j )%state(1)%label = spclabel
               allocate( edge_profiles%ggd( time_sind )%neutral( j )%state(1)% &
                   &      neutral_type%name(1) )
               allocate( edge_profiles%ggd( time_sind )%neutral( j )%state(1)% &
                   &      neutral_type%description(1) )
               edge_profiles%ggd( time_sind )%neutral( j )%state(1)%neutral_type%name = &
                   &     "Thermal"
               edge_profiles%ggd( time_sind )%neutral( j )%state(1)%neutral_type%index = 2
               edge_profiles%ggd( time_sind )%neutral( j )%state(1)%neutral_type%description = &
                   &     "Fluid neutral species from B2.5"
               do i = 1, nsources
                  allocate( edge_sources%source(i)%ggd( time_sind )%neutral( j )%element(1) )
                  allocate( edge_sources%source(i)%ggd( time_sind )%neutral( j )%label(1) )
                  allocate( edge_sources%source(i)%ggd( time_sind )%neutral( j )%state(1) )
                  allocate( edge_sources%source(i)%ggd( time_sind )%neutral( j )%state(1)%label(1) )
                  edge_sources%source(i)%ggd( time_sind )%neutral( j )%element(1)%a = am( is )
                  edge_sources%source(i)%ggd( time_sind )%neutral( j )%element(1)%z_n = &
                     &   zn( is )
#if IMAS_MINOR_VERSION < 15
                  edge_sources%source(i)%ggd( time_sind )%neutral( j )%element(1)%multiplicity = 1.0_IDS_real
#else
                  edge_sources%source(i)%ggd( time_sind )%neutral( j )%element(1)%atoms_n = 1
#endif
                  edge_sources%source(i)%ggd( time_sind )%neutral( j )%label = species_list( js )
                  edge_sources%source(i)%ggd( time_sind )%neutral( j )%ion_index = js
                  edge_sources%source(i)%ggd( time_sind )%neutral( j )%multiple_states_flag = 1
                  call species( is, spclabel, .false. )
                  edge_sources%source(i)%ggd( time_sind )%neutral( j )%state(1)%label = spclabel
                  allocate( edge_sources%source(i)%ggd( time_sind )%neutral( j )%state(1)% &
                     &      neutral_type%name(1) )
                  allocate( edge_sources%source(i)%ggd( time_sind )%neutral( j )%state(1)% &
                     &      neutral_type%description(1) )
                  edge_sources%source(i)%ggd( time_sind )%neutral( j )%state(1)%neutral_type%name = &
                     &     "Thermal"
                  edge_sources%source(i)%ggd( time_sind )%neutral( j )%state(1)%neutral_type%index = 2
                  edge_sources%source(i)%ggd( time_sind )%neutral( j )%state(1)%neutral_type%description = &
                     &     "Fluid neutral species from B2.5"
               end do
               allocate( edge_transport%model(1)%ggd( time_sind )%neutral( j )%element(1) )
               allocate( edge_transport%model(1)%ggd( time_sind )%neutral( j )%label(1) )
               allocate( edge_transport%model(1)%ggd( time_sind )%neutral( j )%state(1) )
               allocate( edge_transport%model(1)%ggd( time_sind )%neutral( j )%state(1)%label(1) )
               edge_transport%model(1)%ggd( time_sind )%neutral( j )%element(1)%a = am( is )
               edge_transport%model(1)%ggd( time_sind )%neutral( j )%element(1)%z_n = zn( is )
#if IMAS_MINOR_VERSION < 15
               edge_transport%model(1)%ggd( time_sind )%neutral( j )%element(1)%multiplicity = 1
#else
               edge_transport%model(1)%ggd( time_sind )%neutral( j )%element(1)%atoms_n = 1
#endif
               edge_transport%model(1)%ggd( time_sind )%neutral( j )%label = species_list( js )
               edge_transport%model(1)%ggd( time_sind )%neutral( j )%ion_index = js
               edge_transport%model(1)%ggd( time_sind )%neutral( j )%multiple_states_flag = 1
               call species( is, spclabel, .false. )
               edge_transport%model(1)%ggd( time_sind )%neutral( j )%state(1)%label = spclabel
               allocate( edge_transport%model(1)%ggd( time_sind )%neutral( j )%state(1)% &
                   &      neutral_type%name(1) )
               allocate( edge_transport%model(1)%ggd( time_sind )%neutral( j )%state(1)% &
                   &      neutral_type%description(1) )
               edge_transport%model(1)%ggd( time_sind )%neutral( j )%state(1)%neutral_type%name = &
                   &     "Thermal"
               edge_transport%model(1)%ggd( time_sind )%neutral( j )%state(1)%neutral_type%index = 2
               edge_transport%model(1)%ggd( time_sind )%neutral( j )%state(1)%neutral_type%description = &
                   &     "Fluid neutral species from B2.5"
            end do
#if IMAS_MINOR_VERSION > 21
            allocate( radiation%process(1)%ggd( time_sind )%neutral( nneut ) )
            j = 0
            do js = 1, nspecies
               is = eb2spcr(js)
               if (.not.is_neutral(is)) cycle
               j = j + 1
               allocate( radiation%process(1)%ggd( time_sind )%neutral( j )%element(1) )
               allocate( radiation%process(1)%ggd( time_sind )%neutral( j )%label(1) )
               allocate( radiation%process(1)%ggd( time_sind )%neutral( j )%state(1) )
               allocate( radiation%process(1)%ggd( time_sind )%neutral( j )%state(1)%label(1) )
               radiation%process(1)%ggd( time_sind )%neutral( j )%element(1)%a = am(is)
               radiation%process(1)%ggd( time_sind )%neutral( j )%element(1)%z_n = zn(is)
               radiation%process(1)%ggd( time_sind )%neutral( j )%element(1)%atoms_n = 1
               radiation%process(1)%ggd( time_sind )%neutral( j )%label = species_list(js)
               radiation%process(1)%ggd( time_sind )%neutral( j )%ion_index = js
               radiation%process(1)%ggd( time_sind )%neutral( j )%multiple_states_flag = 1
               call species( is, spclabel, .false. )
               radiation%process(1)%ggd( time_sind )%neutral( j )%state(1)%label = spclabel
               allocate( radiation%process(1)%ggd( time_sind )%neutral( j )%state(1)% &
                   &     neutral_type%name(1) )
               allocate( radiation%process(1)%ggd( time_sind )%neutral( j )%state(1)% &
                   &     neutral_type%description(1) )
               radiation%process(1)%ggd( time_sind )%neutral( j )%state(1)%neutral_type%name = &
                   &     "Thermal"
               radiation%process(1)%ggd( time_sind )%neutral( j )%state(1)%neutral_type%index = 2
               radiation%process(1)%ggd( time_sind )%neutral( j )%state(1)%neutral_type%description = &
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
            radiation%process(3)%ggd( time_sind )%neutral( js )%ion_index = js
            radiation%process(3)%ggd( time_sind )%neutral( js )%multiple_states_flag = 1
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
            radiation%process(3)%ggd( time_sind )%neutral( js )%label = textmn( j-1 )
            radiation%process(3)%ggd( time_sind )%neutral( js )%ion_index = k
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
          do js = nspecies+1, nsion
            allocate( radiation%process(4)%ggd( time_sind )%ion( js )%label(1) )
            allocate( radiation%process(4)%ggd( time_sind )%ion( js )%state( istion(js) ) )
            do ks = 1, istion(js)
              is = ispion(js,ks)
              allocate( radiation%process(4)%ggd( time_sind )%ion( js )%state( ks )%label(1) )
              radiation%process(4)%ggd( time_sind )%ion( js )%state( ks )%label = textin( is-1 )
            end do
            is = ispion(js,1)
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
            if (istion(js).eq.1) then
              radiation%process(4)%ggd( time_sind )%ion( js )%z_ion = nchrgi( is )
              radiation%process(4)%ggd( time_sind )%ion( js )%label = textin( is-1 )
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
                  radiation%process(4)%ggd( time_sind )%ion( js )%label = ion_label
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
            call divide_by_contact_areas(nx,ny,fne,totFace)
            call write_face_scalar(                                         &
                &   val = edge_transport%model(1)%ggd( time_sind )%         &
                &         electrons%particles%flux,                         &
                &   value = totFace,                                        &
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
              if (is.le.nspecies) then
                tmpCv(:,:) = 0.0_IDS_real
                do js = 1, istion(is)
                  call write_quantity(                                      &
                      &   val = edge_profiles%ggd( time_sind )%             &
                      &         ion( is )%state( js )%density,              &
                      &   value = na(:,:,ispion(is,js)),                    &
                      &   time_sind = time_sind )
                  tmpCv(:,:) = tmpCv(:,:) + na(:,:,ispion(is,js))
                end do
                call write_quantity(                                          &
                    &   val = edge_profiles%ggd( time_sind )%ion(is)%density, &
                    &   value = tmpCv,                                        &
                    &   time_sind = time_sind )
            !! fna: Ion particle flux
                totFace(:,:,0:1) = 0.0_IDS_real
                do js = 1, istion(is)
                  call divide_by_contact_areas                              &
                      &  (nx,ny,fna(-1,-1,0,ispion(is,js)),tmpFace)
                  totFace(:,:,0) = totFace(:,:,0) + tmpFace(:,:,0)
                  totFace(:,:,1) = totFace(:,:,1) + tmpFace(:,:,1)
                  call write_face_scalar(                                   &
                      &   val = edge_transport%model(1)%ggd( time_sind )%   &
                      &         ion( is )%state( js )%particles%flux,       &
                      &   value = tmpFace,                                  &
                      &   time_sind = time_sind )
                end do
                call write_face_scalar(                                   &
                    &   val = edge_transport%model(1)%ggd( time_sind )%   &
                    &         ion( is )%particles%flux,                   &
                    &   value = totFace,                                  &
                    &   time_sind = time_sind )
            !! cdna: Ion diffusivity
                do js = 1, istion(is)
                  call write_face_scalar(                                   &
                      &   val = edge_transport%model(1)%ggd( time_sind )%   &
                      &         ion( is )%state( js )%particles%d,          &
                      &   value = cdna(:,:,:,ispion(is,js)),                &
                      &   time_sind = time_sind )
            !! cvla: Ion diffusivity
                  call write_face_scalar(                                   &
                      &   val = edge_transport%model(1)%ggd( time_sind )%   &
                      &         ion( is )%state( js )%particles%v,          &
                      &   value = cvla(:,:,:,ispion(is,js)),                &
                      &   time_sind = time_sind )
            !! fllim0fna: Ion flux limiter
                  call write_face_scalar(                                   &
                      &   val = edge_transport%model(1)%ggd( time_sind )%   &
                      &         ion( is )%state( js )%                      &
                      &         particles%flux_limiter,                     &
                      &   value = fllim0fna(:,:,:,ispion(is,js)),           &
                      &   time_sind = time_sind )
                end do
            !! sna: Ion particle sources
                totCv(:,:) = 0.0_IDS_real
                do js = 1, istion(is)
                  tmpCv(:,:) = ( sna(:,:,0, ispion(is,js) ) +               &
                      &          sna(:,:,1, ispion(is,js) ) * na(:,:, ispion(is,js) ) ) / vol(:,:)
                  totCv(:,:) = totCv(:,:) + tmpCv(:,:)
                  call write_cell_scalar( scalar = edge_sources%source(1)%  &
                      &   ggd( time_sind )%ion( is )%state( js )%particles, &
                      &   b2CellData = tmpCv )
                end do
                call write_cell_scalar( scalar = edge_sources%source(1)%    &
                    &   ggd( time_sind )%ion( is )%particles,               &
                    &   b2CellData = totCv )
                totCv(:,:) = 0.0_IDS_real
                do js = 1, istion(is)
                  tmpCv(:,:) = ext_sna(:,:,ispion(is,js)) / vol(:,:)
                  totCv(:,:) = totCv(:,:) + tmpCv(:,:)
                  call write_cell_scalar(                                   &
                      &   scalar = edge_sources%source(2)%ggd( time_sind )% &
                      &            ion( is )%state( js )%particles,         &
                      &   b2CellData = tmpCv )
                end do
                call write_cell_scalar(                                     &
                    &   scalar = edge_sources%source(2)%ggd( time_sind )%   &
                    &            ion( is )%particles,                       &
                    &   b2CellData = totCv )
                totCv(:,:) = 0.0_IDS_real
                do js = 1, istion(is)
                  tmpCv(:,:) = ( b2stbc_sna(:,:,ispion(is,js)) +            &
                        &        b2stbm_sna(:,:,ispion(is,js)) ) / vol(:,:)
                  totCv(:,:) = totCv(:,:) + tmpCv(:,:)
                  call write_cell_scalar(                                   &
                      &   scalar = edge_sources%source(3)%ggd( time_sind )% &
                      &            ion( is )%state( js )%particles,         &
                      &   b2CellData = tmpCv )
                end do
                call write_cell_scalar(                                     &
                    &   scalar = edge_sources%source(3)%ggd( time_sind )%   &
                    &            ion( is )%particles,                       &
                    &   b2CellData = totCv )
                totCv(:,:) = 0.0_IDS_real
                do js = 1, istion(is)
                  tmpCv(:,:) = ( snadt(:,:,0,ispion(is,js)) +                &
                        &        snadt(:,:,1,ispion(is,js)) * na(:,:,ispion(is,js)) ) / vol(:,:)
                  totCv(:,:) = totCv(:,:) + tmpCv(:,:)
                  call write_cell_scalar(                                   &
                      &   scalar = edge_sources%source(4)%ggd( time_sind )% &
                      &            ion( is )%state( js )%particles,         &
                      &   b2CellData = tmpCv )
                end do
                call write_cell_scalar(                                     &
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
                    call write_cell_scalar(                                     &
                        &   scalar = edge_sources%source(5)%ggd( time_sind )%   &
                        &            ion( is )%state( js )%particles,           &
                        &   b2CellData = tmpCv )
                  end do
                  call write_cell_scalar(                                       &
                        &   scalar = edge_sources%source(5)%ggd( time_sind )%   &
                        &            ion( is )%particles,                       &
                        &   b2CellData = totCv )
                  totCv = 0.0_IDS_real
                  do js = 1, istion(is)
                    tmpCv = 0.0_IDS_real
                    do istrai = 1, size( eirene_mc_pmpl_sna_bal, 4)
                      tmpCv(:,:) = tmpCv(:,:)                                   &
                         &       + eirene_mc_pmpl_sna_bal(:,:,ispion(is,js),istrai)
                    end do
                    tmpCv(:,:) = tmpCv(:,:) / vol(:,:)
                    totCv(:,:) = totCv(:,:) + tmpCv(:,:)
                    call write_cell_scalar(                                     &
                        &   scalar = edge_sources%source(6)%ggd( time_sind )%   &
                        &            ion( is )%state( js )%particles,           &
                        &   b2CellData = tmpCv )
                  end do
                  call write_cell_scalar(                                       &
                      &   scalar = edge_sources%source(6)%ggd( time_sind )%     &
                      &            ion( is )%particles,                         &
                      &   b2CellData = totCv )
                else
                  totCv(:,:) = 0.0_IDS_real
                  do js = 1, istion(is)
                    tmpCv(:,:) = b2stbr_sna(:,:,ispion(is,js)) / vol(:,:)
                    totCv(:,:) = totCv(:,:) + tmpCv(:,:)
                    call write_cell_scalar(                                     &
                        &   scalar = edge_sources%source(5)%ggd( time_sind )%   &
                        &            ion( is )%state( js )%particles,           &
                        &   b2CellData = tmpCv )
                  end do
                  call write_cell_scalar(                                       &
                        &   scalar = edge_sources%source(5)%ggd( time_sind )%   &
                        &            ion( is )%particles,                       &
                        &   b2CellData = totCv )
                end if
                do js = 1, istion(is)
                  tmpCv(:,:) = rsana(:,:,ispion(is,js)) / vol(:,:)
                  call write_cell_scalar(                                       &
                      &   scalar = edge_sources%source(7)%ggd( time_sind )%     &
                      &            ion( is )%state( js )%particles,             &
                      &   b2CellData = tmpCv )
                  tmpCv(:,:) = rrana(:,:,ispion(is,js)) / vol(:,:)
                  call write_cell_scalar(                                       &
                      &   scalar = edge_sources%source(8)%ggd( time_sind )%     &
                      &            ion( is )%state( js )%particles,             &
                      &   b2CellData = tmpCv )
                  tmpCv(:,:) = rcxna(:,:,ispion(is,js)) / vol(:,:)
                  call write_cell_scalar(                                       &
                      &   scalar = edge_sources%source(9)%ggd( time_sind )%     &
                      &            ion( is )%state( js )%particles,             &
                      &   b2CellData = tmpCv )
                end do
              else
                totCv(:,:) = 0.0_IDS_real
                do js = 1, istion(is)
                  tmpCv(-1:nx,-1:ny) = dib2(0:nx+1,0:ny+1,ispion(is,js),1)
                  totCv(-1:nx,-1:ny) = totCv(-1:nx,-1:ny) + tmpCv(-1:nx,-1:ny)
                  call write_quantity(                                      &
                      &   val = edge_profiles%ggd( time_sind )%             &
                      &         ion(is)%state( js )%density,                &
                      &   value = tmpCv,                                    &
                      &   time_sind = time_sind )
                end do
                call write_quantity(                                        &
                    &   val = edge_profiles%ggd( time_sind )%               &
                    &         ion(is)%density,                              &
                    &   value = totCv,                                      &
                    &   time_sind = time_sind )
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
                    call write_cell_scalar(                                   &
                        &   scalar = edge_sources%source(1)%ggd( time_sind )% &
                        &            ion( is )%state( js )%particles,         &
                        &   b2CellData = tmpCv )
                  end do
                  call write_cell_scalar(                                     &
                      &   scalar = edge_sources%source(1)%ggd( time_sind )%   &
                      &            ion( is )%particles,                       &
                      &   b2CellData = totCv )
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
              if (is.le.nspecies) then
                do js = 1, istion(is)
                !! ua: Parallel ion velocity
                  call write_cell_vector_component(                           &
                      &   vectorComponent = edge_profiles%ggd( time_sind )%   &
                      &                     ion( is )%state( js )%velocity,   &
                      &   b2CellData = ua(:,:,ispion(is,js)),                 &
                      &   vectorID = VEC_ALIGN_PARALLEL_ID )
                  if (drift_style.eq.0) then
                !! wadia: Diamagnetic ion velocity
                    call write_cell_vector_component(                         &
                      &   vectorComponent = edge_profiles%ggd( time_sind )%   &
                      &                     ion( is )%state( js )%            &
                      &                     velocity_diamagnetic,             &
                      &   b2CellData = wadia(:,:,0,ispion(is,js)),            &
                      &   vectorID = VEC_ALIGN_POLOIDAL_ID )
                    call write_cell_vector_component(                         &
                      &   vectorComponent = edge_profiles%ggd( time_sind )%   &
                      &                     ion( is )%state( js )%            &
                      &                     velocity_diamagnetic,             &
                      &   b2CellData = wadia(:,:,1,ispion(is,js)),            &
                      &   vectorID = VEC_ALIGN_RADIAL_ID )
                !! vaecrb: ExB ion velocity
                    call write_cell_vector_component(                         &
                      &   vectorComponent = edge_profiles%ggd( time_sind )%   &
                      &                     ion( is )%state( js )%            &
                      &                     velocity_exb,                     &
                      &   b2CellData = vaecrb(:,:,0,ispion(is,js)),           &
                      &   vectorID = VEC_ALIGN_POLOIDAL_ID )
                    call write_cell_vector_component(                         &
                      &   vectorComponent = edge_profiles%ggd( time_sind )%   &
                      &                     ion( is )%state( js )%            &
                      &                     velocity_exb,                     &
                      &   b2CellData = vaecrb(:,:,1,ispion(is,js)),           &
                      &   vectorID = VEC_ALIGN_RADIAL_ID )
                  else
                !! wadia: Diamagnetic ion velocity
                    tmpFace(:,:,0) = wadia(:,:,0,ispion(is,js))
                    tmpFace(:,:,1) = IDS_REAL_INVALID
                    call write_face_vector_component(                         &
                      &   vectorComponent = edge_profiles%ggd( time_sind )%   &
                      &         ion( is )%state( js )%velocity_diamagnetic,   &
                      &   b2FaceData = tmpFace,                               &
                      &   vectorID = VEC_ALIGN_POLOIDAL_ID )
                    tmpFace(:,:,0) = IDS_REAL_INVALID
                    tmpFace(:,:,1) = wadia(:,:,1,ispion(is,js))
                    call write_face_vector_component(                         &
                      &   vectorComponent = edge_profiles%ggd( time_sind )%   &
                      &         ion( is )%state( js )%velocity_diamagnetic,   &
                      &   b2FaceData = tmpFace,                               &
                      &   vectorID = VEC_ALIGN_RADIAL_ID )
                !! vaecrb: ExB ion velocity
                    tmpFace(:,:,0) = vaecrb(:,:,0,ispion(is,js))
                    tmpFace(:,:,1) = IDS_REAL_INVALID
                    call write_face_vector_component(                         &
                      &   vectorComponent = edge_profiles%ggd( time_sind )%   &
                      &         ion( is )%state( js )%velocity_exb,           &
                      &   b2FaceData = tmpFace,                               &
                      &   vectorID = VEC_ALIGN_POLOIDAL_ID )
                    tmpFace(:,:,0) = IDS_REAL_INVALID
                    tmpFace(:,:,1) = vaecrb(:,:,1,ispion(is,js))
                    call write_face_vector_component(                         &
                      &   vectorComponent = edge_profiles%ggd( time_sind )%   &
                      &         ion( is )%state( js )%velocity_exb,           &
                      &   b2FaceData = tmpFace,                               &
                      &   vectorID = VEC_ALIGN_RADIAL_ID )
                  end if
                !! cvsa: Ion diffusivity
                  call write_face_vector_component(                           &
                      &   vectorComponent = edge_transport%model(1)%          &
                      &         ggd( time_sind )%ion( is )%state( js )%       &
                      &         momentum%d,                                   &
                      &   b2FaceData = cvsa(:,:,:,ispion(is,js)),             &
                      &   vectorID = VEC_ALIGN_PARALLEL_ID )
                end do
                !! fmo: Ion momentum flux
                totFace(:,:,:) = 0.0_IDS_real
                do js = 1, istion(is)
                  call divide_by_contact_areas                                &
                      &  (nx,ny,fmo(-1,-1,0,ispion(is,js)),tmpFace)
                  call write_face_vector_component(                           &
                      &   vectorComponent = edge_transport%model(1)%          &
                      &                     ggd( time_sind )%ion( is )%       &
                      &                     state( js )%momentum%flux,        &
                      &   b2FaceData = tmpFace,                               &
                      &   vectorID = VEC_ALIGN_PARALLEL_ID )
                  totFace(:,:,0) = totFace(:,:,0) + tmpFace(:,:,0)
                  totFace(:,:,1) = totFace(:,:,1) + tmpFace(:,:,1)
                end do
                call write_face_vector_component(                             &
                    &   vectorComponent = edge_transport%model(1)%            &
                    &                     ggd( time_sind )%ion( is )%         &
                    &                     momentum%flux,                      &
                    &   b2FaceData = totFace,                                 &
                    &   vectorID = VEC_ALIGN_PARALLEL_ID )
                !! fllimvisc: Ion parallel momentum transport flux limit
                do js = 1, istion(is)
                  tmpFace(:,:,0) = fllimvisc(:,:,ispion(is,js))
                  tmpFace(:,:,1) = 1.0_IDS_real
                  call write_face_vector_component(                           &
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
                  call write_cell_vector_component(                            &
                      &   vectorComponent = edge_sources%source(1)%            &
                      &                     ggd( time_sind )%ion( is )%        &
                      &                     state( js )%momentum,              &
                      &   b2CellData = tmpCv,                                  &
                      &   vectorID = VEC_ALIGN_PARALLEL_ID )
                end do
                call write_cell_vector_component(                              &
                    &   vectorComponent = edge_sources%source(1)%              &
                    &                     ggd( time_sind )%ion( is )%momentum, &
                    &   b2CellData = totCv,                                    &
                    &   vectorID = VEC_ALIGN_PARALLEL_ID )
                totCv(:,:) = 0.0_IDS_real
                do js = 1, istion(is)
                  tmpCv(:,:) = ext_smo(:,:,ispion(is,js)) / vol(:,:)
                  totCv(:,:) = totCv(:,:) + tmpCv(:,:)
                  call write_cell_vector_component(                            &
                      &   vectorComponent = edge_sources%source(2)%            &
                      &                     ggd( time_sind )%ion( is )%        &
                      &                     state( js )%momentum,              &
                      &   b2CellData = tmpCv,                                  &
                      &   vectorID = VEC_ALIGN_PARALLEL_ID )
                end do
                call write_cell_vector_component(                              &
                      &   vectorComponent = edge_sources%source(2)%            &
                      &                     ggd( time_sind )%ion( is )%        &
                      &                     momentum,                          &
                      &   b2CellData = totCv,                                  &
                      &   vectorID = VEC_ALIGN_PARALLEL_ID )
                totCv(:,:) = 0.0_IDS_real
                do js = 1, istion(is)
                  tmpCv(:,:) = ( b2stbc_smo(:,:,ispion(is,js)) +               &
                        &        b2stbm_smo(:,:,ispion(is,js)) ) / vol(:,:)
                  totCv(:,:) = totCv(:,:) + tmpCv(:,:)
                  call write_cell_vector_component(                            &
                      &   vectorComponent = edge_sources%source(3)%            &
                      &                     ggd( time_sind )%ion( is )%        &
                      &                     state( js )%momentum,              &
                      &   b2CellData = tmpCv,                                  &
                      &   vectorID = VEC_ALIGN_PARALLEL_ID )
                end do
                call write_cell_vector_component(                              &
                      &   vectorComponent = edge_sources%source(3)%            &
                      &                     ggd( time_sind )%ion( is )%        &
                      &                     momentum,                          &
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
                  call write_cell_vector_component(                            &
                      &   vectorComponent = edge_sources%source(4)%            &
                      &                     ggd( time_sind )%ion( is )%        &
                      &                     state( js )%momentum,              &
                      &   b2CellData = tmpCv,                                  &
                      &   vectorID = VEC_ALIGN_PARALLEL_ID )
                end do
                call write_cell_vector_component(                              &
                      &   vectorComponent = edge_sources%source(4)%            &
                      &                     ggd( time_sind )%ion( is )%        &
                      &                     momentum,                          &
                      &   b2CellData = totCv,                                  &
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
                    call write_cell_vector_component(                          &
                        &   vectorComponent = edge_sources%source(5)%          &
                        &            ggd( time_sind )%ion( is )%               &
                        &            state( js )%momentum,                     &
                        &   b2CellData = tmpCv,                                &
                        &   vectorID = VEC_ALIGN_PARALLEL_ID )
                  end do
                  call write_cell_vector_component(                            &
                      &   vectorComponent = edge_sources%source(5)%            &
                      &            ggd( time_sind )%ion( is )%momentum,        &
                      &   b2CellData = totCv,                                  &
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
                    call write_cell_vector_component(                          &
                        &   vectorComponent = edge_sources%source(6)%          &
                        &            ggd( time_sind )%ion( is )%               &
                        &            state( js )%momentum,                     &
                        &   b2CellData = tmpCv,                                &
                        &   vectorID = VEC_ALIGN_PARALLEL_ID )
                  end do
                  call write_cell_vector_component(                            &
                        &   vectorComponent = edge_sources%source(6)%          &
                        &            ggd( time_sind )%ion( is )%momentum,      &
                        &   b2CellData = totCv,                                &
                        &   vectorID = VEC_ALIGN_PARALLEL_ID )
                else
                  totCv = 0.0_IDS_real
                  do js = 1, istion(is)
                    tmpCv(:,:) = b2stbr_smo(:,:,ispion(is,js)) / vol(:,:)
                    totCv(:,:) = totCv(:,:) + tmpCv(:,:)
                    call write_cell_vector_component(                          &
                        &   vectorComponent = edge_sources%source(5)%          &
                        &            ggd( time_sind )%ion( is )%               &
                        &            state( js )%momentum,                     &
                        &   b2CellData = tmpCv,                                &
                        &   vectorID = VEC_ALIGN_PARALLEL_ID )
                  end do
                  call write_cell_vector_component(                            &
                        &   vectorComponent = edge_sources%source(5)%          &
                        &            ggd( time_sind )%ion( is )%momentum,      &
                        &   b2CellData = totCv,                                &
                        &   vectorID = VEC_ALIGN_PARALLEL_ID )
                end if
                do js = 1, istion(is)
                  tmpCv(:,:) = rsamo(:,:,ispion(is,js)) / vol(:,:)
                  call write_cell_vector_component(                            &
                      &   vectorComponent = edge_sources%source(7)%            &
                      &                     ggd( time_sind )%ion( is )%        &
                      &                     state( js )%momentum,              &
                      &   b2CellData = tmpCv,                                  &
                      &   vectorID = VEC_ALIGN_PARALLEL_ID )
                  tmpCv(:,:) = rramo(:,:,ispion(is,js)) / vol(:,:)
                  call write_cell_vector_component(                            &
                      &   vectorComponent = edge_sources%source(8)%            &
                      &                     ggd( time_sind )%ion( is )%        &
                      &                     state( js )%momentum,              &
                      &   b2CellData = tmpCv,                                  &
                      &   vectorID = VEC_ALIGN_PARALLEL_ID )
                  tmpCv(:,:) = rcxmo(:,:,ispion(is,js)) / vol(:,:)
                  call write_cell_vector_component(                            &
                      &   vectorComponent = edge_sources%source(9)%            &
                      &                     ggd( time_sind )%ion( is )%        &
                      &                     state( js )%momentum,              &
                      &   b2CellData = tmpCv,                                  &
                      &   vectorID = VEC_ALIGN_PARALLEL_ID )
                end do
              end if
            end do

            !! te: Electron Temperature
            tmpCv(:,:) = te(:,:)/qe
            call write_quantity(                                        &
                &   val = edge_profiles%ggd( time_sind )%electrons%     &
                &         temperature,                                  &
                &   value = tmpCv,                                      &
                &   time_sind = time_sind )
            tmpCv(:,:) = hce0(:,:)/ne(:,:)
            call write_quantity(                                        &
                &   val = edge_transport%model(1)%ggd( time_sind )%     &
                &         electrons%energy%d,                           &
                &   value = tmpCv,                                      &
                &   time_sind = time_sind )
            call write_face_scalar(                                     &
                &   val = edge_transport%model(1)%ggd( time_sind )%     &
                &         electrons%energy%v,                           &
                &   value = chve,                                       &
                &   time_sind = time_sind )
            call divide_by_contact_areas(nx,ny,fhe,totFace)
            call write_face_scalar(                                     &
                &   val = edge_transport%model(1)%ggd( time_sind )%     &
                &         electrons%energy%flux,                        &
                &   value = totFace,                                    &
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
            tmpCv(:,:) = hci0(:,:)/ni(:,:,0)
            call write_quantity(                                      &
                &   val = edge_transport%model(1)%ggd( time_sind )%   &
                &         total_ion_energy%d,                         &
                &   value = tmpCv,                                    &
                &   time_sind = time_sind )
            call write_face_scalar(                                   &
                 &   val = edge_transport%model(1)%ggd( time_sind )%  &
                 &         total_ion_energy%v,                        &
                 &   value = chvi,                                    &
                 &   time_sind = time_sind )
            !! fhi : Ion heat flux
            call divide_by_contact_areas(nx,ny,fhi,totFace)
            call write_face_scalar(                                     &
                &   val = edge_transport%model(1)%ggd( time_sind )%     &
                &         total_ion_energy%flux,                        &
                &   value = totFace,                                    &
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
              if (is.le.nspecies) then
                totCv(:,:) = 0.0_IDS_real
                do js = 1, istion(is)
            !! Ion energy sources resolved by species
                  tmpCv(:,:) = rsahi(:,:,ispion(is,js)) / vol(:,:)
                  totCv(:,:) = totCv(:,:) + tmpCv(:,:)
                  call write_cell_scalar(                                     &
                      &   scalar = edge_sources%source(7)%ggd( time_sind )%   &
                      &            ion( is )%state( js )%energy,              &
                      &   b2CellData = tmpCv )
                end do
                call write_cell_scalar(                                       &
                      &   scalar = edge_sources%source(7)%ggd( time_sind )%   &
                      &            ion( is )%energy,                          &
                      &   b2CellData = totCv )
                totCv(:,:) = 0.0_IDS_real
                do js = 1, istion(is)
                  tmpCv(:,:) = rrahi(:,:,ispion(is,js)) / vol(:,:)
                  totCv(:,:) = totCv(:,:) + tmpCv(:,:)
                  call write_cell_scalar(                                     &
                      &   scalar = edge_sources%source(8)%ggd( time_sind )%   &
                      &            ion( is )%state( js )%energy,              &
                      &   b2CellData = tmpCv )
                end do
                call write_cell_scalar(                                       &
                    &   scalar = edge_sources%source(8)%ggd( time_sind )%     &
                    &            ion( is )%energy,                            &
                        &   b2CellData = totCv )
                totCv(:,:) = 0.0_IDS_real
                do js = 1, istion(is)
                  tmpCv(:,:) = rcxhi(:,:,ispion(is,js)) / vol(:,:)
                  totCv(:,:) = totCv(:,:) + tmpCv(:,:)
                  call write_cell_scalar(                                     &
                      &   scalar = edge_sources%source(9)%ggd( time_sind )%   &
                      &            ion( is )%state( js )%energy,              &
                      &   b2CellData = tmpCv )
                end do
                call write_cell_scalar(                                       &
                    &   scalar = edge_sources%source(9)%ggd( time_sind )%     &
                    &            ion( is )%energy,                            &
                    &   b2CellData = totCv )
              else
#ifdef B25_EIRENE
                do js = 1, istion(is)
            !! Test ion temperature
                  tmpCv(-1:nx,-1:ny) = tib2(0:nx+1,0:ny+1,ispion(is,js),1)/qe
                  call write_quantity(                                        &
                      &   val = edge_profiles%ggd( time_sind )%ion( is )%     &
                      &         state( js )%temperature,                      &
                      &   value = tmpCv,                                      &
                      &   time_sind = time_sind )
                end do
#endif
              end if
            end do

            do is = 1, nsion
              if (is.le.nspecies) then
                !! pb : Ion pressure
                totCv(:,:) = 0.0_IDS_real
                do js = 1, istion(is)
                  call b2xppb( nx, ny, rza(:,:,ispion(is,js)),              &
                      &                 na(:,:,ispion(is,js)), te, ti, pb)
                  totCv(:,:) = totCv(:,:) + pb(:,:)
                  call write_quantity(                                      &
                      &   val = edge_profiles%ggd( time_sind )%ion( is )%   &
                      &         state( js )%pressure,                       &
                      &   value = pb,                                       &
                      &   time_sind = time_sind )
                end do
                call write_quantity(                                        &
                    &   val = edge_profiles%ggd( time_sind )%ion( is )%     &
                    &         pressure,                                     &
                    &   value = totCv,                                      &
                    &   time_sind = time_sind )
                !! Kinetic energy density
                totCv(:,:) = 0.0_IDS_real
                do js = 1, istion(is)
                  tmpCv(:,:) = 0.5_IDS_real*am(ispion(is,js))*mp*           &
                      &       (ua(:,:,ispion(is,js))**2)*na(:,:,ispion(is,js))
                  totCv(:,:) = totCv(:,:) + tmpCv(:,:)
                  call write_quantity(                                      &
                      &   val = edge_profiles%ggd( time_sind )%ion( is )%   &
                      &         state( js )%energy_density_kinetic,         &
                      &   value = tmpCv,                                    &
                      &   time_sind = time_sind )
                end do
                call write_quantity(                                        &
                    &   val = edge_profiles%ggd( time_sind )%ion( is )%     &
                    &         energy_density_kinetic,                       &
                    &   value = totCv,                                      &
                    &   time_sind = time_sind )
                do js = 1, istion(is)
                !! Average charge
                  call write_quantity(                                      &
                      &   val = edge_profiles%ggd( time_sind )%ion( is )%   &
                      &         state( js )%z_average,                      &
                      &   value = rza(:,:,ispion(is,js)),                   &
                      &   time_sind = time_sind )
                !! Average square charge
                  call write_quantity(                                      &
                      &   val = edge_profiles%ggd( time_sind )%ion( is )%   &
                      &         state( js )%z_square_average,               &
                      &   value = rz2(:,:,ispion(is,js)),                   &
                      &   time_sind = time_sind )
                !! Ionization potential
                  call write_quantity(                                      &
                      &   val = edge_profiles%ggd( time_sind )%ion( is )%   &
                      &         state( js )%ionisation_potential,           &
                      &   value = rpt(:,:,ispion(is,js)),                   &
                      &   time_sind = time_sind )
                end do
              else
#ifdef B25_EIRENE
                !! Test ion pressure
                totCv(:,:) = 0.0_IDS_real
                do js = 1, istion(is)
                  tmpCv(-1:nx,-1:ny) = dib2(0:nx+1,0:ny+1,ispion(is,js),1)  &
                      &               *tib2(0:nx+1,0:ny+1,ispion(is,js),1)
                  totCv(-1:nx,-1:ny) = totCv(-1:nx,-1:ny) + tmpCv(-1:nx,-1:ny)
                  call write_quantity(                                      &
                      &   val = edge_profiles%ggd( time_sind )%ion( is )%   &
                      &         state( js )%pressure,                       &
                      &   value = tmpCv,                                    &
                      &   time_sind = time_sind )
                end do
                call write_quantity(                                        &
                    &   val = edge_profiles%ggd( time_sind )%ion( is )%     &
                    &         pressure,                                     &
                    &   value = totCv,                                      &
                    &   time_sind = time_sind )
                !! Average charge
                do js = 1, istion(is)
                  tmpCv(:,:) = nchrgi( ispion(is,js) )
                  call write_quantity(                                      &
                      &   val = edge_profiles%ggd( time_sind )%ion( is )%   &
                      &         state( js )%z_average,                      &
                      &   value = tmpCv,                                    &
                      &   time_sind = time_sind )
                !! Average square charge
                  tmpCv(:,:) = nchrgi( js )**2
                  call write_quantity(                                      &
                      &   val = edge_profiles%ggd( time_sind )%ion( is )%   &
                      &         state( js )%z_square_average,               &
                      &   value = tmpCv,                                    &
                      &   time_sind = time_sind )
                end do
                !! Radiation source
                totCv = 0.0_IDS_real
                do js = 1, istion(is)
                  tmpCv(-1:nx,-1:ny) = eionrad(0:nx+1,0:ny+1,ispion(is,js),0) / vol(-1:nx,-1:ny)
                  totCv(-1:nx,-1:ny) = totCv(-1:nx,-1:ny) + tmpCv(-1:nx,-1:ny)
                  call write_cell_scalar(                                   &
                      &   scalar = edge_sources%source(12)%                 &
                      &            ggd( time_sind )%ion( is )%              &
                      &            state( js )%energy,                      &
                      &   b2CellData = tmpCv )
                end do
                call write_cell_scalar(                                     &
                    &   scalar = edge_sources%source(12)%                   &
                    &            ggd( time_sind )%ion( is )%energy,         &
                    &   b2CellData = totCv )
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

#if IMAS_MINOR_VERSION > 32
            !! fch: Total current
            call divide_by_contact_areas(nx,ny,fch,tmpFace)
            totFace(:,:,0) = tmpFace(:,:,0)
            totFace(:,:,1) = IDS_REAL_INVALID
            call write_face_vector_component(                                &
                &   vectorComponent = edge_profiles%ggd( time_sind )%        &
                &                     j_total,                               &
                &   b2FaceData = totFace,                                    &
                &   vectorID = VEC_ALIGN_POLOIDAL_ID )
            totFace(:,:,0) = IDS_REAL_INVALID
            totFace(:,:,1) = tmpFace(:,:,1)
            call write_face_vector_component(                                &
                &   vectorComponent = edge_profiles%ggd( time_sind )%        &
                &                     j_total,                               &
                &   b2FaceData = totFace,                                    &
                &   vectorID = VEC_ALIGN_RADIAL_ID )

            !! fch_p: Parallel current
            call divide_by_poloidal_areas(nx,ny,fch_p,tmpFace)
            call write_face_scalar(                                          &
                &   val = edge_profiles%ggd( time_sind )%j_parallel,         &
                &   value = tmpFace,                                         &
                &   time_sind = time_sind )
#endif

            !! fchanml: Anomalous current
            call b2tanml (nx, ny, ns, csig_an, po, fchanml_a, fchanml)
            call divide_by_contact_areas(nx,ny,fchanml,tmpFace)
            totFace(:,:,0) = tmpFace(:,:,0)
            totFace(:,:,1) = IDS_REAL_INVALID
            call write_face_vector_component(                                &
                &   vectorComponent = edge_profiles%ggd( time_sind )%        &
                &                     j_anomalous,                           &
                &   b2FaceData = totFace,                                    &
                &   vectorID = VEC_ALIGN_POLOIDAL_ID )
            totFace(:,:,0) = IDS_REAL_INVALID
            totFace(:,:,1) = tmpFace(:,:,1)
            call write_face_vector_component(                                &
                &   vectorComponent = edge_profiles%ggd( time_sind )%        &
                &                     j_anomalous,                           &
                &   b2FaceData = totFace,                                    &
                &   vectorID = VEC_ALIGN_RADIAL_ID )

            !! fchinert: Inertial current
            call divide_by_contact_areas(nx,ny,fchinert,tmpFace)
            totFace(:,:,0) = tmpFace(:,:,0)
            totFace(:,:,1) = IDS_REAL_INVALID
            call write_face_vector_component(                                &
                &   vectorComponent = edge_profiles%ggd( time_sind )%        &
                &                     j_inertial,                            &
                &   b2FaceData = totFace,                                    &
                &   vectorID = VEC_ALIGN_POLOIDAL_ID )
            totFace(:,:,0) = IDS_REAL_INVALID
            totFace(:,:,1) = tmpFace(:,:,1)
            call write_face_vector_component(                                &
                &   vectorComponent = edge_profiles%ggd( time_sind )%        &
                &                     j_inertial,                            &
                &   b2FaceData = totFace,                                    &
                &   vectorID = VEC_ALIGN_RADIAL_ID )

            !! fchin: Ion-neutral friction current
            call divide_by_contact_areas(nx,ny,fchin,tmpFace)
            totFace(:,:,0) = tmpFace(:,:,0)
            totFace(:,:,1) = IDS_REAL_INVALID
            call write_face_vector_component(                                &
                &   vectorComponent = edge_profiles%ggd( time_sind )%        &
                &                     j_ion_neutral_friction,                &
                &   b2FaceData = totFace,                                    &
                &   vectorID = VEC_ALIGN_POLOIDAL_ID )
            totFace(:,:,0) = IDS_REAL_INVALID
            totFace(:,:,1) = tmpFace(:,:,1)
            call write_face_vector_component(                                &
                &   vectorComponent = edge_profiles%ggd( time_sind )%        &
                &                     j_ion_neutral_friction,                &
                &   b2FaceData = totFace,                                    &
                &   vectorID = VEC_ALIGN_RADIAL_ID )

            !! fchvispar: Parallel viscosity current
            call divide_by_contact_areas(nx,ny,fchvispar,tmpFace)
            totFace(:,:,0) = tmpFace(:,:,0)
            totFace(:,:,1) = IDS_REAL_INVALID
            call write_face_vector_component(                                &
                &   vectorComponent = edge_profiles%ggd( time_sind )%        &
                &                     j_parallel_viscosity,                  &
                &   b2FaceData = totFace,                                    &
                &   vectorID = VEC_ALIGN_POLOIDAL_ID )
            totFace(:,:,0) = IDS_REAL_INVALID
            totFace(:,:,1) = tmpFace(:,:,1)
            call write_face_vector_component(                                &
                &   vectorComponent = edge_profiles%ggd( time_sind )%        &
                &                     j_parallel_viscosity,                  &
                &   b2FaceData = totFace,                                    &
                &   vectorID = VEC_ALIGN_RADIAL_ID )

            !! fchvisper: Perpendicular viscosity current
            call divide_by_contact_areas(nx,ny,fchvisper,tmpFace)
            totFace(:,:,0) = tmpFace(:,:,0)
            totFace(:,:,1) = IDS_REAL_INVALID
            call write_face_vector_component(                                &
                &   vectorComponent = edge_profiles%ggd( time_sind )%        &
                &                     j_perpendicular_viscosity,             &
                &   b2FaceData = totFace,                                    &
                &   vectorID = VEC_ALIGN_POLOIDAL_ID )
            totFace(:,:,0) = IDS_REAL_INVALID
            totFace(:,:,1) = tmpFace(:,:,1)
            call write_face_vector_component(                                &
                &   vectorComponent = edge_profiles%ggd( time_sind )%        &
                &                     j_perpendicular_viscosity,             &
                &   b2FaceData = totFace,                                    &
                &   vectorID = VEC_ALIGN_RADIAL_ID )

            !! fchvisq: Heat viscosity current
            call divide_by_contact_areas(nx,ny,fchvisq,tmpFace)
            totFace(:,:,0) = tmpFace(:,:,0)
            totFace(:,:,1) = IDS_REAL_INVALID
            call write_face_vector_component(                                &
                &   vectorComponent = edge_profiles%ggd( time_sind )%        &
                &                     j_heat_viscosity,                      &
                &   b2FaceData = totFace,                                    &
                &   vectorID = VEC_ALIGN_POLOIDAL_ID )
            totFace(:,:,0) = IDS_REAL_INVALID
            totFace(:,:,1) = tmpFace(:,:,1)
            call write_face_vector_component(                                &
                &   vectorComponent = edge_profiles%ggd( time_sind )%        &
                &                     j_heat_viscosity,                      &
                &   b2FaceData = totFace,                                    &
                &   vectorID = VEC_ALIGN_RADIAL_ID )

            !! fchdia: Diamagnetic current
            call divide_by_contact_areas(nx,ny,fchdia,tmpFace)
            totFace(:,:,0) = tmpFace(:,:,0)
            totFace(:,:,1) = IDS_REAL_INVALID
            call write_face_vector_component(                                &
                &   vectorComponent = edge_profiles%ggd( time_sind )%        &
                &                     j_diamagnetic,                         &
                &   b2FaceData = totFace,                                    &
                &   vectorID = VEC_ALIGN_POLOIDAL_ID )
            totFace(:,:,0) = IDS_REAL_INVALID
            totFace(:,:,1) = tmpFace(:,:,1)
            call write_face_vector_component(                                &
                &   vectorComponent = edge_profiles%ggd( time_sind )%        &
                &                     j_diamagnetic,                         &
                &   b2FaceData = totFace,                                    &
                &   vectorID = VEC_ALIGN_RADIAL_ID )

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
                !! Neutral pressure
                   call write_quantity(                                      &
                       &   val = edge_profiles%ggd( time_sind )%             &
                       &         neutral( is )%pressure,                     &
                       &   value = tmpCv,                                    &
                       &   time_sind = time_sind )
                !! Neutral density
                   call write_quantity(                                      &
                       &   val = edge_profiles%ggd( time_sind )%             &
                       &         neutral( is )%density,                      &
                       &   value = totCv,                                    &
                       &   time_sind = time_sind )
                !! Neutral particle flux
                   call write_face_scalar(                                   &
                       &   val = edge_transport%model(1)%ggd( time_sind )%   &
                       &         neutral( is )%particles%flux,               &
                       &   value = tmpFace,                                  &
                       &   time_sind = time_sind )
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
                     call write_cell_vector_component(                       &
                       &   vectorComponent = edge_profiles%ggd( time_sind )% &
                       &                     neutral( is )%velocity,         &
                       &   b2CellData = tmpCv(:,:),                          &
                       &   vectorID = VEC_ALIGN_POLOIDAL_ID )
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
                     call write_cell_vector_component(                       &
                       &   vectorComponent = edge_profiles%ggd( time_sind )% &
                       &                     neutral( is )%velocity,         &
                       &   b2CellData = tmpCv(:,:),                          &
                       &   vectorID = VEC_ALIGN_RADIAL_ID )
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
                     call write_cell_vector_component(                       &
                       &   vectorComponent = edge_profiles%ggd( time_sind )% &
                       &                     neutral( is )%velocity,         &
                       &   b2CellData = tmpCv(:,:),                          &
                       &   vectorID = VEC_ALIGN_TOROIDAL_ID )
                   end if
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
                   if (allocated(un0)) then
                     call write_cell_vector_component(                       &
                       &   vectorComponent = edge_profiles%ggd( time_sind )% &
                       &         neutral( js )%state( ks )%velocity,         &
                       &   b2CellData = un0(:,:,0,is),                       &
                       &   vectorID = VEC_ALIGN_POLOIDAL_ID )
                     call write_cell_vector_component(                       &
                       &   vectorComponent = edge_profiles%ggd( time_sind )% &
                       &         neutral( js )%state( ks )%velocity,         &
                       &   b2CellData = un0(:,:,1,is),                       &
                       &   vectorID = VEC_ALIGN_RADIAL_ID )
                     call write_cell_vector_component(                       &
                       &   vectorComponent = edge_profiles%ggd( time_sind )% &
                       &         neutral( js )%state( ks )%velocity,         &
                       &   b2CellData = un0(:,:,2,is),                       &
                       &   vectorID = VEC_ALIGN_TOROIDAL_ID )
                   end if
                   if (drift_style.eq.0) then
                     tmpCv = 0.0_IDS_real
                     call write_cell_vector_component(                       &
                       &   vectorComponent = edge_profiles%ggd( time_sind )% &
                       &                     neutral( js )%state( ks )%      &
                       &                     velocity_diamagnetic,           &
                       &   b2CellData = tmpCv,                               &
                       &   vectorID = VEC_ALIGN_POLOIDAL_ID )
                     call write_cell_vector_component(                       &
                       &   vectorComponent = edge_profiles%ggd( time_sind )% &
                       &                     neutral( js )%state( ks )%      &
                       &                     velocity_diamagnetic,           &
                       &   b2CellData = tmpCv,                               &
                       &   vectorID = VEC_ALIGN_RADIAL_ID )
                     call write_cell_vector_component(                       &
                       &   vectorComponent = edge_profiles%ggd( time_sind )% &
                       &                     neutral( js )%state( ks )%      &
                       &                     velocity_exb,                   &
                       &   b2CellData = tmpCv,                               &
                       &   vectorID = VEC_ALIGN_POLOIDAL_ID )
                     call write_cell_vector_component(                       &
                       &   vectorComponent = edge_profiles%ggd( time_sind )% &
                       &                     neutral( js )%state( ks )%      &
                       &                     velocity_exb,                   &
                       &   b2CellData = tmpCv,                               &
                       &   vectorID = VEC_ALIGN_RADIAL_ID )
                   else
                     tmpFace = 0.0_IDS_real
                     call write_face_vector_component(                       &
                       &   vectorComponent = edge_profiles%ggd( time_sind )% &
                       &                     neutral( js )%state( ks )%      &
                       &                     velocity_diamagnetic,           &
                       &   b2FaceData = tmpFace,                             &
                       &   vectorID = VEC_ALIGN_POLOIDAL_ID )
                     call write_face_vector_component(                       &
                       &   vectorComponent = edge_profiles%ggd( time_sind )% &
                       &                     neutral( js )%state( ks )%      &
                       &                     velocity_diamagnetic,           &
                       &   b2FaceData = tmpFace,                             &
                       &   vectorID = VEC_ALIGN_RADIAL_ID )
                     call write_face_vector_component(                       &
                       &   vectorComponent = edge_profiles%ggd( time_sind )% &
                       &                     neutral( js )%state( ks )%      &
                       &                     velocity_exb,                   &
                       &   b2FaceData = tmpFace,                             &
                       &   vectorID = VEC_ALIGN_POLOIDAL_ID )
                     call write_face_vector_component(                       &
                       &   vectorComponent = edge_profiles%ggd( time_sind )% &
                       &                     neutral( js )%state( ks )%      &
                       &                     velocity_exb,                   &
                       &   b2FaceData = tmpFace,                             &
                       &   vectorID = VEC_ALIGN_RADIAL_ID )
                   end if
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
                 !! Molecular pressure
                   call write_quantity(                                      &
                       &   val = edge_profiles%ggd( time_sind )%             &
                       &         neutral( js )%pressure,                     &
                       &   value = tmpCv,                                    &
                       &   time_sind = time_sind )
                 !! Molecular density
                   call write_quantity(                                      &
                       &   val = edge_profiles%ggd( time_sind )%             &
                       &         neutral( js )%density,                      &
                       &   value = totCv,                                    &
                       &   time_sind = time_sind )
                 !! Molecular particular flux
                   call write_face_scalar(                                   &
                       &   val = edge_transport%model(1)%ggd( time_sind )%   &
                       &         neutral( js )%particles%flux,               &
                       &   value = tmpFace,                                  &
                       &   time_sind = time_sind )
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
                     call write_cell_vector_component(                       &
                       &   vectorComponent = edge_profiles%ggd( time_sind )% &
                       &                     neutral( js )%velocity,         &
                       &   b2CellData = tmpCv(:,:),                          &
                       &   vectorID = VEC_ALIGN_POLOIDAL_ID )
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
                     call write_cell_vector_component(                       &
                       &   vectorComponent = edge_profiles%ggd( time_sind )% &
                       &                     neutral( js )%velocity,         &
                       &   b2CellData = tmpCv(:,:),                          &
                       &   vectorID = VEC_ALIGN_RADIAL_ID )
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
                     call write_cell_vector_component(                       &
                       &   vectorComponent = edge_profiles%ggd( time_sind )% &
                       &                     neutral( js )%velocity,         &
                       &   b2CellData = tmpCv(:,:),                          &
                       &   vectorID = VEC_ALIGN_TOROIDAL_ID )
                   end if
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
                   if (allocated(um0)) then
                     call write_cell_vector_component(                       &
                       &   vectorComponent = edge_profiles%ggd( time_sind )% &
                       &         neutral( js )%state( ks )%velocity,         &
                       &   b2CellData = um0(:,:,0,is),                       &
                       &   vectorID = VEC_ALIGN_POLOIDAL_ID )
                     call write_cell_vector_component(                       &
                       &   vectorComponent = edge_profiles%ggd( time_sind )% &
                       &         neutral( js )%state( ks )%velocity,         &
                       &   b2CellData = um0(:,:,1,is),                       &
                       &   vectorID = VEC_ALIGN_RADIAL_ID )
                     call write_cell_vector_component(                       &
                       &   vectorComponent = edge_profiles%ggd( time_sind )% &
                       &         neutral( js )%state( ks )%velocity,         &
                       &   b2CellData = um0(:,:,2,is),                       &
                       &   vectorID = VEC_ALIGN_TOROIDAL_ID )
                   end if
                   if (drift_style.eq.0) then
                     tmpCv = 0.0_IDS_real
                     call write_cell_vector_component(                       &
                       &   vectorComponent = edge_profiles%ggd( time_sind )% &
                       &                     neutral( js )%state( ks )%      &
                       &                     velocity_diamagnetic,           &
                       &   b2CellData = tmpCv,                               &
                       &   vectorID = VEC_ALIGN_POLOIDAL_ID )
                     call write_cell_vector_component(                       &
                       &   vectorComponent = edge_profiles%ggd( time_sind )% &
                       &                     neutral( js )%state( ks )%      &
                       &                     velocity_diamagnetic,           &
                       &   b2CellData = tmpCv,                               &
                       &   vectorID = VEC_ALIGN_RADIAL_ID )
                     call write_cell_vector_component(                       &
                       &   vectorComponent = edge_profiles%ggd( time_sind )% &
                       &                     neutral( js )%state( ks )%      &
                       &                     velocity_exb,                   &
                       &   b2CellData = tmpCv,                               &
                       &   vectorID = VEC_ALIGN_POLOIDAL_ID )
                     call write_cell_vector_component(                       &
                       &   vectorComponent = edge_profiles%ggd( time_sind )% &
                       &                     neutral( js )%state( ks )%      &
                       &                     velocity_exb,                   &
                       &   b2CellData = tmpCv,                               &
                       &   vectorID = VEC_ALIGN_RADIAL_ID )
                   else
                     tmpFace = 0.0_IDS_real
                     call write_face_vector_component(                       &
                       &   vectorComponent = edge_profiles%ggd( time_sind )% &
                       &                     neutral( js )%state( ks )%      &
                       &                     velocity_diamagnetic,           &
                       &   b2FaceData = tmpFace,                             &
                       &   vectorID = VEC_ALIGN_POLOIDAL_ID )
                     call write_face_vector_component(                       &
                       &   vectorComponent = edge_profiles%ggd( time_sind )% &
                       &                     neutral( js )%state( ks )%      &
                       &                     velocity_diamagnetic,           &
                       &   b2FaceData = tmpFace,                             &
                       &   vectorID = VEC_ALIGN_RADIAL_ID )
                     call write_face_vector_component(                       &
                       &   vectorComponent = edge_profiles%ggd( time_sind )% &
                       &                     neutral( js )%state( ks )%      &
                       &                     velocity_exb,                   &
                       &   b2FaceData = tmpFace,                             &
                       &   vectorID = VEC_ALIGN_POLOIDAL_ID )
                     call write_face_vector_component(                       &
                       &   vectorComponent = edge_profiles%ggd( time_sind )% &
                       &                     neutral( js )%state( ks )%      &
                       &                     velocity_exb,                   &
                       &   b2FaceData = tmpFace,                             &
                       &   vectorID = VEC_ALIGN_RADIAL_ID )
                   end if
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
                j = 0
                do is = 1, nspecies
                    js = eb2spcr(is)
                    if (.not.is_neutral(js)) cycle
                    j = j + 1
                !! na : Fluid neutral density
                    call write_quantity(                                      &
                        &   val = edge_profiles%ggd( time_sind )%             &
                        &         neutral( j )%density,                       &
                        &   value = na(:,:,js),                               &
                        &   time_sind = time_sind )
                    call write_quantity(                                      &
                        &   val = edge_profiles%ggd( time_sind )%             &
                        &         neutral( j )%state(1)%density,              &
                        &   value = na(:,:,js),                               &
                        &   time_sind = time_sind )
                !! ua: Parallel fluid neutral velocity
                    call write_cell_vector_component(                         &
                        &   vectorComponent = edge_profiles%ggd( time_sind )% &
                        &                     neutral( j )%velocity,          &
                        &   b2CellData = ua(:,:,js),                          &
                        &   vectorID = VEC_ALIGN_PARALLEL_ID )
                    call write_cell_vector_component(                         &
                        &   vectorComponent = edge_profiles%ggd( time_sind )% &
                        &                     neutral( j )%state(1)%velocity, &
                        &   b2CellData = ua(:,:,js),                          &
                        &   vectorID = VEC_ALIGN_PARALLEL_ID )
                    if (drift_style.eq.0) then
                !! wadia: Diamagnetic fluid neutral velocity
                      call write_cell_vector_component(                       &
                        &   vectorComponent = edge_profiles%ggd( time_sind )% &
                        &                     neutral( j )%state(1)%          &
                        &                     velocity_diamagnetic,           &
                        &   b2CellData = wadia(:,:,0,js),                     &
                        &   vectorID = VEC_ALIGN_POLOIDAL_ID )
                      call write_cell_vector_component(                       &
                        &   vectorComponent = edge_profiles%ggd( time_sind )% &
                        &                     neutral( j )%state(1)%          &
                        &                     velocity_diamagnetic,           &
                        &   b2CellData = wadia(:,:,1,js),                     &
                        &   vectorID = VEC_ALIGN_RADIAL_ID )
                !! vaecrb: ExB fluid neutral velocity
                      call write_cell_vector_component(                       &
                        &   vectorComponent = edge_profiles%ggd( time_sind )% &
                        &                     neutral( j )%state(1)%          &
                        &                     velocity_exb,                   &
                        &   b2CellData = vaecrb(:,:,0,js),                    &
                        &   vectorID = VEC_ALIGN_POLOIDAL_ID )
                      call write_cell_vector_component(                       &
                        &   vectorComponent = edge_profiles%ggd( time_sind )% &
                        &                     neutral( j )%state(1)%          &
                        &                     velocity_exb,                   &
                        &   b2CellData = vaecrb(:,:,1,js),                    &
                        &   vectorID = VEC_ALIGN_RADIAL_ID )
                    else
                !! wadia: Diamagnetic fluid neutral velocity
                      tmpFace(:,:,0) = wadia(:,:,0,js)
                      tmpFace(:,:,1) = IDS_REAL_INVALID
                      call write_face_vector_component(                       &
                        &   vectorComponent = edge_profiles%ggd( time_sind )% &
                        &         neutral( j )%state(1)%velocity_diamagnetic, &
                        &   b2FaceData = tmpFace,                             &
                        &   vectorID = VEC_ALIGN_POLOIDAL_ID )
                      tmpFace(:,:,0) = IDS_REAL_INVALID
                      tmpFace(:,:,1) = wadia(:,:,1,js)
                      call write_face_vector_component(                       &
                        &   vectorComponent = edge_profiles%ggd( time_sind )% &
                        &         neutral( j )%state(1)%velocity_diamagnetic, &
                        &   b2FaceData = tmpFace,                             &
                        &   vectorID = VEC_ALIGN_RADIAL_ID )
                !! vaecrb: ExB fluid neutral velocity
                      tmpFace(:,:,0) = vaecrb(:,:,0,js)
                      tmpFace(:,:,1) = IDS_REAL_INVALID
                      call write_face_vector_component(                       &
                        &   vectorComponent = edge_profiles%ggd( time_sind )% &
                        &         neutral( j )%state(1)%velocity_exb,         &
                        &   b2FaceData = tmpFace,                             &
                        &   vectorID = VEC_ALIGN_POLOIDAL_ID )
                      tmpFace(:,:,0) = IDS_REAL_INVALID
                      tmpFace(:,:,1) = vaecrb(:,:,1,js)
                      call write_face_vector_component(                       &
                        &   vectorComponent = edge_profiles%ggd( time_sind )% &
                        &         neutral( j )%state(1)%velocity_exb,         &
                        &   b2FaceData = tmpFace,                             &
                        &   vectorID = VEC_ALIGN_RADIAL_ID )
                    end if
                !! fna: Fluid neutral particle flux
                    call divide_by_contact_areas(nx,ny,fna(-1,-1,0,js),totFace)
                    call write_face_scalar(                                   &
                        &   val = edge_transport%model(1)%ggd( time_sind )%   &
                        &         neutral( j )%particles%flux,                &
                        &   value = totFace,                                  &
                        &   time_sind = time_sind )
                    call write_face_scalar(                                   &
                        &   val = edge_transport%model(1)%ggd( time_sind )%   &
                        &         neutral( j )%state(1)%particles%flux,       &
                        &   value = totFace,                                  &
                        &   time_sind = time_sind )
                !! pb : Fluid neutral pressure
                    call b2xppb( nx, ny, rza(:,:,js), na(:,:,js), te, ti, pb)
                    call write_quantity(                                      &
                        &   val = edge_profiles%ggd( time_sind )%             &
                        &         neutral( j )%pressure,                      &
                        &   value = pb,                                       &
                        &   time_sind = time_sind )
                    call write_quantity(                                      &
                        &   val = edge_profiles%ggd( time_sind )%             &
                        &         neutral( j )%state(1)%pressure,             &
                        &   value = pb,                                       &
                        &   time_sind = time_sind )
                !! cdpa: Fluid neutral diffusivity
                    call write_face_scalar(                                   &
                        &   val = edge_transport%model(1)%ggd( time_sind )%   &
                        &         neutral( j )%particles%d,                   &
                        &   value = cdpa(:,:,:,js),                           &
                        &   time_sind = time_sind )
                    call write_face_scalar(                                   &
                        &   val = edge_transport%model(1)%ggd( time_sind )%   &
                        &         neutral( j )%state(1)%particles%d,          &
                        &   value = cdpa(:,:,:,js),                           &
                        &   time_sind = time_sind )
                !! Fluid neutral kinetic energy density
                    tmpCv(:,:) = 0.5_IDS_real*am(js)*mp*(ua(:,:,js)**2)*na(:,:,js)
                    call write_quantity(                                      &
                        &   val = edge_profiles%ggd( time_sind )%             &
                        &         neutral( j )%energy_density_kinetic,        &
                        &   value = tmpCv,                                    &
                        &   time_sind = time_sind )
                    call write_quantity(                                      &
                        &   val = edge_profiles%ggd( time_sind )%             &
                        &         neutral( j )%state(1)%                      &
                        &         energy_density_kinetic,                     &
                        &   value = tmpCv,                                    &
                        &   time_sind = time_sind )
                !! cvsa: Ion diffusivity
                    call write_face_vector_component(                         &
                        &   vectorComponent = edge_transport%model(1)%        &
                        &                     ggd( time_sind )%neutral( j )%  &
                        &                     momentum%d,                     &
                        &   b2FaceData = cvsa(:,:,:,js),                      &
                        &   vectorID = VEC_ALIGN_PARALLEL_ID )
                    call write_face_vector_component(                         &
                        &   vectorComponent = edge_transport%model(1)%        &
                        &                     ggd( time_sind )%neutral( j )%  &
                        &                     state(1)%momentum%d,            &
                        &   b2FaceData = cvsa(:,:,:,js),                      &
                        &   vectorID = VEC_ALIGN_PARALLEL_ID )
                !! fmo: Ion momentum flux
                    call divide_by_contact_areas(nx,ny,fmo(-1,-1,0,js),tmpFace)
                    call write_face_vector_component(                         &
                        &   vectorComponent = edge_transport%model(1)%        &
                        &                     ggd( time_sind )%neutral( j )%  &
                        &                     momentum%flux,                  &
                        &   b2FaceData = tmpFace,                             &
                        &   vectorID = VEC_ALIGN_PARALLEL_ID )
                    call write_face_vector_component(                         &
                        &   vectorComponent = edge_transport%model(1)%        &
                        &                     ggd( time_sind )%neutral( j )%  &
                        &                     state(1)%momentum%flux,         &
                        &   b2FaceData = tmpFace,                             &
                        &   vectorID = VEC_ALIGN_PARALLEL_ID )
                !! fllim0fna: Fluid neutral flux limiter
                    call write_face_scalar(                                   &
                        &   val = edge_transport%model(1)%ggd( time_sind )%   &
                        &         neutral( j )%particles%flux_limiter,        &
                        &   value = fllim0fna(:,:,:,js),                      &
                        &   time_sind = time_sind )
                    call write_face_scalar(                                   &
                        &   val = edge_transport%model(1)%ggd( time_sind )%   &
                        &         neutral( j )%state(1)%                      &
                        &         particles%flux_limiter,                     &
                        &   value = fllim0fna(:,:,:,js),                      &
                        &   time_sind = time_sind )
                !! fllimvisc: Fluid neutral momentum transport flux limit
                    tmpFace(:,:,0) = fllimvisc(:,:,js)
                    tmpFace(:,:,1) = 1.0_IDS_real
                    call write_face_vector_component(                         &
                        &   vectorComponent = edge_transport%model(1)%        &
                        &                     ggd( time_sind )%neutral( j )%  &
                        &                     momentum%flux_limiter,          &
                        &   b2FaceData = fllimvisc(:,:,js),                   &
                        &   vectorID = VEC_ALIGN_PARALLEL_ID )
                    call write_face_vector_component(                         &
                        &   vectorComponent = edge_transport%model(1)%        &
                        &                     ggd( time_sind )%neutral( j )%  &
                        &                     state(1)%momentum%flux_limiter, &
                        &   b2FaceData = fllimvisc(:,:,js),                   &
                        &   vectorID = VEC_ALIGN_PARALLEL_ID )
                !! sna: Fluid neutral particle sources
                    tmpCv(:,:) = ( sna(:,:,0,js) +                            &
                        &          sna(:,:,1,js) * na(:,:,js) ) / vol(:,:)
                    call write_cell_scalar( scalar = edge_sources%source(1)%  &
                        &   ggd( time_sind )%neutral( j )%particles,          &
                        &   b2CellData = tmpCv )
                    call write_cell_scalar( scalar = edge_sources%source(1)%  &
                        &   ggd( time_sind )%neutral( j )%state(1)%particles, &
                        &   b2CellData = tmpCv )
                    tmpCv(:,:) = ext_sna(:,:,js) / vol(:,:)
                    call write_cell_scalar(                                   &
                        &   scalar = edge_sources%source(2)%ggd( time_sind )% &
                        &            neutral( j )%particles,                  &
                        &   b2CellData = tmpCv )
                    call write_cell_scalar(                                   &
                        &   scalar = edge_sources%source(2)%ggd( time_sind )% &
                        &            neutral( j )%state(1)%particles,         &
                        &   b2CellData = tmpCv )
                    tmpCv(:,:) = ( b2stbc_sna(:,:,js) +                       &
                        &          b2stbm_sna(:,:,js) ) / vol(:,:)
                    call write_cell_scalar(                                   &
                        &   scalar = edge_sources%source(3)%ggd( time_sind )% &
                        &            neutral( j )%particles,                  &
                        &   b2CellData = tmpCv )
                    call write_cell_scalar(                                   &
                        &   scalar = edge_sources%source(3)%ggd( time_sind )% &
                        &            neutral( j )%state(1)%particles,         &
                        &   b2CellData = tmpCv )
                    tmpCv(:,:) = ( snadt(:,:,0,js) +                          &
                        &          snadt(:,:,1,js) * na(:,:,js) ) / vol(:,:)
                    call write_cell_scalar(                                   &
                        &   scalar = edge_sources%source(4)%ggd( time_sind )% &
                        &            neutral( j )%particles,                  &
                        &   b2CellData = tmpCv )
                    call write_cell_scalar(                                   &
                        &   scalar = edge_sources%source(4)%ggd( time_sind )% &
                        &            neutral( j )%state(1)%particles,         &
                        &   b2CellData = tmpCv )
                    tmpCv(:,:) = b2stbr_sna(:,:,js) / vol(:,:)
                    call write_cell_scalar(                                   &
                        &   scalar = edge_sources%source(5)%ggd( time_sind )% &
                        &            neutral( j )%particles,                  &
                        &   b2CellData = tmpCv )
                    call write_cell_scalar(                                   &
                        &   scalar = edge_sources%source(5)%ggd( time_sind )% &
                        &            neutral( j )%state(1)%particles,         &
                        &   b2CellData = tmpCv )
                    tmpCv(:,:) = rsana(:,:,js) / vol(:,:)
                    call write_cell_scalar(                                   &
                        &   scalar = edge_sources%source(7)%ggd( time_sind )% &
                        &            neutral( j )%particles,                  &
                        &   b2CellData = tmpCv )
                    call write_cell_scalar(                                   &
                        &   scalar = edge_sources%source(7)%ggd( time_sind )% &
                        &            neutral( j )%state(1)%particles,         &
                        &   b2CellData = tmpCv )
                    tmpCv(:,:) = rrana(:,:,js) / vol(:,:)
                    call write_cell_scalar(                                   &
                        &   scalar = edge_sources%source(8)%ggd( time_sind )% &
                        &            neutral( j )%particles,                  &
                        &   b2CellData = tmpCv )
                    call write_cell_scalar(                                   &
                        &   scalar = edge_sources%source(8)%ggd( time_sind )% &
                        &            neutral( j )%state(1)%particles,         &
                        &   b2CellData = tmpCv )
                    tmpCv(:,:) = rcxna(:,:,js) / vol(:,:)
                    call write_cell_scalar(                                   &
                        &   scalar = edge_sources%source(9)%ggd( time_sind )% &
                        &            neutral( j )%particles,                  &
                        &   b2CellData = tmpCv )
                    call write_cell_scalar(                                   &
                        &   scalar = edge_sources%source(9)%ggd( time_sind )% &
                        &            neutral( j )%state(1)%particles,         &
                        &   b2CellData = tmpCv )
                !! smo: Ion parallel momentum sources
                    do iy = -1, ny
                      do ix = -1, nx
                        tmpCv(ix,iy) = ( smo(ix,iy,0,js) +                  &
                        &                smo(ix,iy,1,js) * ua(ix,iy,js) +   &
                        &                smo(ix,iy,2,js) * roxa(ix,iy,js) + &
                        &                smo(ix,iy,3,js) * roxa(ix,iy,js)   &
                        &                                * ua(ix,iy,js) ) / vol(ix,iy)
                      end do
                    end do
                    call write_cell_vector_component(                         &
                        &   vectorComponent = edge_sources%source(1)%         &
                        &                     ggd( time_sind )%               &
                        &                     neutral( j )%momentum,          &
                        &   b2CellData = tmpCv,                               &
                        &   vectorID = VEC_ALIGN_PARALLEL_ID )
                    call write_cell_vector_component(                         &
                        &   vectorComponent = edge_sources%source(1)%         &
                        &                     ggd( time_sind )%neutral( j )%  &
                        &                     state(1)%momentum,              &
                        &   b2CellData = tmpCv,                               &
                        &   vectorID = VEC_ALIGN_PARALLEL_ID )
                    tmpCv(:,:) = ext_smo(:,:,js) / vol(:,:)
                    call write_cell_vector_component(                         &
                        &   vectorComponent = edge_sources%source(2)%         &
                        &                     ggd( time_sind )%neutral( j )%  &
                        &                     momentum,                       &
                        &   b2CellData = tmpCv,                               &
                        &   vectorID = VEC_ALIGN_PARALLEL_ID )
                    call write_cell_vector_component(                         &
                        &   vectorComponent = edge_sources%source(2)%         &
                        &                     ggd( time_sind )%neutral( j )%  &
                        &                     state(1)%momentum,              &
                        &   b2CellData = tmpCv,                               &
                        &   vectorID = VEC_ALIGN_PARALLEL_ID )
                    tmpCv(:,:) = ( b2stbc_smo(:,:,js) +                       &
                        &          b2stbm_smo(:,:,js) ) / vol(:,:)
                    call write_cell_vector_component(                         &
                        &   vectorComponent = edge_sources%source(3)%         &
                        &                     ggd( time_sind )%neutral( j )%  &
                        &                     momentum,                       &
                        &   b2CellData = tmpCv,                               &
                        &   vectorID = VEC_ALIGN_PARALLEL_ID )
                    call write_cell_vector_component(                         &
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
                    call write_cell_vector_component(                         &
                        &   vectorComponent = edge_sources%source(4)%         &
                        &                     ggd( time_sind )%neutral( j )%  &
                        &                     momentum,                       &
                        &   b2CellData = tmpCv,                               &
                        &   vectorID = VEC_ALIGN_PARALLEL_ID )
                    call write_cell_vector_component(                         &
                        &   vectorComponent = edge_sources%source(4)%         &
                        &                     ggd( time_sind )%neutral( j )%  &
                        &                     state(1)%momentum,              &
                        &   b2CellData = tmpCv,                               &
                        &   vectorID = VEC_ALIGN_PARALLEL_ID )
                    tmpCv(:,:) = b2stbr_smo(:,:,js) / vol(:,:)
                    call write_cell_vector_component(                         &
                        &   vectorComponent = edge_sources%source(5)%         &
                        &            ggd( time_sind )%neutral( j )%momentum,  &
                        &   b2CellData = tmpCv,                               &
                        &   vectorID = VEC_ALIGN_PARALLEL_ID )
                    call write_cell_vector_component(                         &
                        &   vectorComponent = edge_sources%source(5)%         &
                        &            ggd( time_sind )%neutral( j )%           &
                        &            state(1)%momentum,                       &
                        &   b2CellData = tmpCv,                               &
                        &   vectorID = VEC_ALIGN_PARALLEL_ID )
                    tmpCv(:,:) = rsamo(:,:,js) / vol(:,:)
                    call write_cell_vector_component(                         &
                        &   vectorComponent = edge_sources%source(7)%         &
                        &                     ggd( time_sind )%neutral( j )%  &
                        &                     momentum,                       &
                        &   b2CellData = tmpCv,                               &
                        &   vectorID = VEC_ALIGN_PARALLEL_ID )
                    call write_cell_vector_component(                         &
                        &   vectorComponent = edge_sources%source(7)%         &
                        &                     ggd( time_sind )%neutral( j )%  &
                        &                     state(1)%momentum,              &
                        &   b2CellData = tmpCv,                               &
                        &   vectorID = VEC_ALIGN_PARALLEL_ID )
                    tmpCv(:,:) = rramo(:,:,js) / vol(:,:)
                    call write_cell_vector_component(                         &
                        &   vectorComponent = edge_sources%source(8)%         &
                        &                     ggd( time_sind )%neutral( j )%  &
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
                        &                     ggd( time_sind )%neutral( j )%  &
                        &                     momentum,                       &
                        &   b2CellData = tmpCv,                               &
                        &   vectorID = VEC_ALIGN_PARALLEL_ID )
                    call write_cell_vector_component(                         &
                        &   vectorComponent = edge_sources%source(9)%         &
                        &                     ggd( time_sind )%neutral( j )%  &
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
            tmpFace(:,:,0) = csig(:,:,0)
            tmpFace(:,:,1) = IDS_REAL_INVALID
            call write_face_vector_component(                           &
                &   vectorComponent = edge_transport%model(1)%          &
                &         ggd( time_sind )%conductivity,                &
                &   b2FaceData = tmpFace,                               &
                &   vectorID = VEC_ALIGN_POLOIDAL_ID )
            tmpFace(:,:,0) = IDS_REAL_INVALID
            tmpFace(:,:,1) = csig(:,:,1)
            call write_face_vector_component(                           &
                &   vectorComponent = edge_transport%model(1)%          &
                &         ggd( time_sind )%conductivity,                &
                &   b2FaceData = tmpFace,                               &
                &   vectorID = VEC_ALIGN_RADIAL_ID )

            !! B (magnetic field vector)
            !! Compute unit basis vectors along the field directions
            call compute_Coordinate_Unit_Vectors(crx, cry, e(:,:,:,1),  &
                &   e(:,:,:,2), e(:,:,:,3))

#if IMAS_MINOR_VERSION > 21
            !! write the emissivity data
            !! Process 1. Line and recombination radiation from B2.5 ions
            j = 0
            do is = 1, nspecies
              js = eb2spcr(is)
              if (is_neutral(js) .and. use_eirene.eq.0) then
                j = j + 1
                tmpCv(:,:) = rqrad(:,:,js) / vol(:,:)
                call write_cell_scalar( scalar = radiation%process(1)%     &
                    &   ggd( time_sind )%neutral( j )%emissivity,          &
                    &   b2CellData = tmpCv )
#if IMAS_MINOR_VERSION > 24
                call write_cell_scalar( scalar = radiation%process(1)%     &
                    &   ggd( time_sind )%neutral( j )%                     &
                    &   state(1)%emissivity,                               &
                    &   b2CellData = tmpCv )
#endif
              end if
              totCv(:,:) = 0.0_IDS_real
              do js = 1, istion(is)
                tmpCv(:,:) = rqrad(:,:,ispion(is,js)) / vol(:,:)
                totCv(:,:) = totCv(:,:) + tmpCv(:,:)
#if IMAS_MINOR_VERSION > 24
                call write_cell_scalar( scalar = radiation%process(1)%     &
                    &   ggd( time_sind )%ion( is )%state( js )%emissivity, &
                    &   b2CellData = tmpCv )
#endif
              end do
              call write_cell_scalar( scalar = radiation%process(1)%       &
                    &   ggd( time_sind )%ion( is )%emissivity,             &
                    &   b2CellData = totCv )
            end do
            !! Process 2. Bremsstrahlung from B2.5 ions
            totCv(:,:) = 0.0_IDS_real
            do is = 1, nspecies
              do js = 1, istion(is)
                tmpCv(:,:) = rqbrm(:,:,ispion(is,js)) / vol(:,:)
                totCv(:,:) = totCv(:,:) + tmpCv(:,:)
#if IMAS_MINOR_VERSION > 24
                call write_cell_scalar( scalar = radiation%process(2)%     &
                    &   ggd( time_sind )%ion( is )%state( js )%emissivity, &
                    &   b2CellData = tmpCv )
#endif
              end do
              call write_cell_scalar( scalar = radiation%process(2)%       &
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
              do is = nspecies+1, nsion
                totCv(:,:) = 0.0_R8
                do js = 1, istion(is)
                  do ix = -1, nx
                    do iy = -1, ny
                      tmpCv(ix,iy) = -eionrad(ix+1,iy+1,ispion(is,js),0) / vol(ix,iy)
                    end do
                  end do
                  totCv(:,:) = totCv(:,:) + tmpCv(:,:)
#if IMAS_MINOR_VERSION > 24
                  call write_cell_scalar( scalar = radiation%process(4)%      &
                      &   ggd( time_sind )%ion( is )%state( js )%emissivity,  &
                      &   b2CellData = tmpCv )
#endif
                end do
                call write_cell_scalar( scalar = radiation%process(4)%        &
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

#if IMAS_MINOR_VERSION > 21
! Summary separatrix data
        if (maxval(abs(fpsi(-1:nx,-1:ny,0:3))).gt.0.0_R8) then
           allocate( summary%local%separatrix%position%psi( num_time_slices ) )
           summary%local%separatrix%position%psi( time_sind ) = fpsi(jxa,jsep,2)
#if IMAS_MINOR_VERSION > 34
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
#if IMAS_MINOR_VERSION > 34
        tmpCv = ne(:,:)
        u = separatrix_average( te, tmpCv )
        call write_sourced_value( summary%local%separatrix_average%t_e, u )
        tmpCv = ni(:,:,1)
        u = separatrix_average( ti, tmpCv )
        call write_sourced_value( summary%local%separatrix_average%t_i_average, u )
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
            call write_sourced_value( summary%local%separatrix%velocity_tor%hydrogen, vtor )
#if IMAS_MINOR_VERSION > 34
            call write_sourced_value( summary%local%separatrix_average%n_i%hydrogen, u )
            call write_sourced_value( summary%local%separatrix_average%velocity_tor%hydrogen, v )
#endif
          case ('D')
            call write_sourced_value( summary%local%separatrix%n_i%deuterium, nisep )
            call write_sourced_value( summary%local%separatrix%velocity_tor%deuterium, vtor )
#if IMAS_MINOR_VERSION > 34
            call write_sourced_value( summary%local%separatrix_average%n_i%deuterium, u )
            call write_sourced_value( summary%local%separatrix_average%velocity_tor%deuterium, v )
#endif
          case ('T')
            call write_sourced_value( summary%local%separatrix%n_i%tritium, nisep )
            call write_sourced_value( summary%local%separatrix%velocity_tor%tritium, vtor )
#if IMAS_MINOR_VERSION > 34
            call write_sourced_value( summary%local%separatrix_average%n_i%tritium, u )
            call write_sourced_value( summary%local%separatrix_average%velocity_tor%tritium, v )
#endif
          case ('He')
            if (nint(am(eb2spcr(is))).eq.3) then
              call write_sourced_value( summary%local%separatrix%n_i%helium_3, nisep )
              call write_sourced_value( summary%local%separatrix%velocity_tor%helium_3, vtor )
#if IMAS_MINOR_VERSION > 34
              call write_sourced_value( summary%local%separatrix_average%n_i%helium_3, u )
              call write_sourced_value( summary%local%separatrix_average%velocity_tor%helium_3, v )
#endif
            else if (nint(am(eb2spcr(is))).eq.4) then
              call write_sourced_value( summary%local%separatrix%n_i%helium_4, nisep )
              call write_sourced_value( summary%local%separatrix%velocity_tor%helium_4, vtor )
#if IMAS_MINOR_VERSION > 34
            call write_sourced_value( summary%local%separatrix_average%n_i%helium_4, u )
            call write_sourced_value( summary%local%separatrix_average%velocity_tor%helium_4, v )
#endif
            end if
          case ('Li')
            call write_sourced_value( summary%local%separatrix%n_i%lithium, nisep )
            call write_sourced_value( summary%local%separatrix%velocity_tor%lithium, vtor )
#if IMAS_MINOR_VERSION > 34
            call write_sourced_value( summary%local%separatrix_average%n_i%lithium, u )
            call write_sourced_value( summary%local%separatrix_average%velocity_tor%lithium, v )
#endif
          case ('Be')
            call write_sourced_value( summary%local%separatrix%n_i%beryllium, nisep )
            call write_sourced_value( summary%local%separatrix%velocity_tor%beryllium, vtor )
#if IMAS_MINOR_VERSION > 34
            call write_sourced_value( summary%local%separatrix_average%n_i%beryllium, u )
            call write_sourced_value( summary%local%separatrix_average%velocity_tor%beryllium, v )
#endif
          case ('C')
            call write_sourced_value( summary%local%separatrix%n_i%carbon, nisep )
            call write_sourced_value( summary%local%separatrix%velocity_tor%carbon, vtor )
#if IMAS_MINOR_VERSION > 34
            call write_sourced_value( summary%local%separatrix_average%n_i%carbon, u )
            call write_sourced_value( summary%local%separatrix_average%velocity_tor%carbon, v )
#endif
          case ('N')
            call write_sourced_value( summary%local%separatrix%n_i%nitrogen, nisep )
            call write_sourced_value( summary%local%separatrix%velocity_tor%nitrogen, vtor )
#if IMAS_MINOR_VERSION > 34
            call write_sourced_value( summary%local%separatrix_average%n_i%nitrogen, u )
            call write_sourced_value( summary%local%separatrix_average%velocity_tor%nitrogen, v )
#endif
          case ('O')
            call write_sourced_value( summary%local%separatrix%n_i%oxygen, nisep )
            call write_sourced_value( summary%local%separatrix%velocity_tor%oxygen, vtor )
#if IMAS_MINOR_VERSION > 34
            call write_sourced_value( summary%local%separatrix_average%n_i%oxygen, u )
            call write_sourced_value( summary%local%separatrix_average%velocity_tor%oxygen, v )
#endif
          case ('Ne')
            call write_sourced_value( summary%local%separatrix%n_i%neon, nisep )
            call write_sourced_value( summary%local%separatrix%velocity_tor%neon, vtor )
#if IMAS_MINOR_VERSION > 34
            call write_sourced_value( summary%local%separatrix_average%n_i%neon, u )
            call write_sourced_value( summary%local%separatrix_average%velocity_tor%neon, v )
#endif
          case ('Ar')
            call write_sourced_value( summary%local%separatrix%n_i%argon, nisep )
            call write_sourced_value( summary%local%separatrix%velocity_tor%argon, vtor )
#if IMAS_MINOR_VERSION > 34
            call write_sourced_value( summary%local%separatrix_average%n_i%argon, u )
            call write_sourced_value( summary%local%separatrix_average%velocity_tor%argon, v )
#endif
#if IMAS_MINOR_VERSION > 30
          case ('Fe')
            call write_sourced_value( summary%local%separatrix%n_i%iron, nisep )
            call write_sourced_value( summary%local%separatrix%velocity_tor%iron, vtor )
#if IMAS_MINOR_VERSION > 34
            call write_sourced_value( summary%local%separatrix_average%n_i%iron, u )
            call write_sourced_value( summary%local%separatrix_average%velocity_tor%iron, v )
#endif
          case ('Kr')
            call write_sourced_value( summary%local%separatrix%n_i%krypton, nisep )
            call write_sourced_value( summary%local%separatrix%velocity_tor%krypton, vtor )
#if IMAS_MINOR_VERSION > 34
            call write_sourced_value( summary%local%separatrix_average%n_i%krypton, u )
            call write_sourced_value( summary%local%separatrix_average%velocity_tor%krypton, v )
#endif
#endif
          case ('Xe')
            call write_sourced_value( summary%local%separatrix%n_i%xenon, nisep )
            call write_sourced_value( summary%local%separatrix%velocity_tor%xenon, vtor )
#if IMAS_MINOR_VERSION > 34
            call write_sourced_value( summary%local%separatrix_average%n_i%xenon, u )
            call write_sourced_value( summary%local%separatrix_average%velocity_tor%xenon, v )
#endif
          case ('W')
            call write_sourced_value( summary%local%separatrix%n_i%tungsten, nisep )
            call write_sourced_value( summary%local%separatrix%velocity_tor%tungsten, vtor )
#if IMAS_MINOR_VERSION > 34
            call write_sourced_value( summary%local%separatrix_average%n_i%tungsten, u )
            call write_sourced_value( summary%local%separatrix_average%velocity_tor%tungsten, v )
#endif
          end select
        end do
        call write_sourced_value( summary%local%separatrix%n_i_total, &
          & 0.5_R8 * (ni(jxa,jsep,1) + ni(topix(jxa,jsep),topiy(jxa,jsep),1)) )
        u = separatrix_average( ni(:,:,1), tmpCv )
#if IMAS_MINOR_VERSION > 34
        call write_sourced_value( summary%local%separatrix_average%n_i_total, u )
#endif
        u = separatrix_average( zeff, tmpCv )
        call write_sourced_value( summary%local%separatrix%zeff, &
          & 0.5_R8 * (zeff(jxa,jsep) + zeff(topix(jxa,jsep),topiy(jxa,jsep))) )
#if IMAS_MINOR_VERSION > 34
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
#if IMAS_MINOR_VERSION > 30
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
#if IMAS_MINOR_VERSION > 31
          u = max( power_flux_peak(itrg(1)), power_flux_peak(itrg(2)) )
          call write_sourced_value( summary%local%limiter%power_flux_peak, u )
#endif
        end if

! Summary divertor plate data
        if (ntrgts.gt.0) then
#if IMAS_MINOR_VERSION > 34
          allocate ( summary%local%divertor_target( ntrgts ) )
#else
          allocate ( summary%local%divertor_plate( ntrgts ) )
#endif
          do i = 1, ntrgts
#if IMAS_MINOR_VERSION > 34
            call write_sourced_string( summary%local%divertor_target(i)%name, plate_name(i) )
            call write_sourced_value( summary%local%divertor_target(i)%t_e, &
              &  0.5_R8 * (te(ixpos(itrg(i)),iypos(itrg(i)))+  &
              &            te(topix(ixpos(itrg(i)),iypos(itrg(i))), &
              &               topiy(ixpos(itrg(i)),iypos(itrg(i)))))/ev )
            call write_sourced_value( summary%local%divertor_target(i)%t_i_average, &
              &  0.5_R8 * (te(ixpos(itrg(i)),iypos(itrg(i)))+  &
              &            te(topix(ixpos(itrg(i)),iypos(itrg(i))), &
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
              &  0.5_R8 * (te(ixpos(itrg(i)),iypos(itrg(i)))+  &
              &            te(topix(ixpos(itrg(i)),iypos(itrg(i))), &
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
#if IMAS_MINOR_VERSION > 34
                call write_sourced_value( summary%local%divertor_target(i)%n_i%hydrogen, nisep )
#else
                call write_sourced_value( summary%local%divertor_plate(i)%n_i%hydrogen, nisep )
#endif
              case ('D')
#if IMAS_MINOR_VERSION > 34
                call write_sourced_value( summary%local%divertor_target(i)%n_i%deuterium, nisep )
#else
                call write_sourced_value( summary%local%divertor_plate(i)%n_i%deuterium, nisep )
#endif
              case ('T')
#if IMAS_MINOR_VERSION > 34
                call write_sourced_value( summary%local%divertor_target(i)%n_i%tritium, nisep )
#else
                call write_sourced_value( summary%local%divertor_plate(i)%n_i%tritium, nisep )
#endif
              case ('He')
                if (nint(am(eb2spcr(is))).eq.3) then
#if IMAS_MINOR_VERSION > 34
                  call write_sourced_value( summary%local%divertor_target(i)%n_i%helium_3, nisep )
#else
                  call write_sourced_value( summary%local%divertor_plate(i)%n_i%helium_3, nisep )
#endif
                else if (nint(am(eb2spcr(is))).eq.4) then
#if IMAS_MINOR_VERSION > 34
                  call write_sourced_value( summary%local%divertor_target(i)%n_i%helium_4, nisep )
#else
                  call write_sourced_value( summary%local%divertor_plate(i)%n_i%helium_4, nisep )
#endif
                end if
              case ('Li')
#if IMAS_MINOR_VERSION > 34
                call write_sourced_value( summary%local%divertor_target(i)%n_i%lithium, nisep )
#else
                call write_sourced_value( summary%local%divertor_plate(i)%n_i%lithium, nisep )
#endif
              case ('Be')
#if IMAS_MINOR_VERSION > 34
                call write_sourced_value( summary%local%divertor_target(i)%n_i%beryllium, nisep )
#else
                call write_sourced_value( summary%local%divertor_plate(i)%n_i%beryllium, nisep )
#endif
              case ('C')
#if IMAS_MINOR_VERSION > 34
                call write_sourced_value( summary%local%divertor_target(i)%n_i%carbon, nisep )
#else
                call write_sourced_value( summary%local%divertor_plate(i)%n_i%carbon, nisep )
#endif
              case ('N')
#if IMAS_MINOR_VERSION > 34
                call write_sourced_value( summary%local%divertor_target(i)%n_i%nitrogen, nisep )
#else
                call write_sourced_value( summary%local%divertor_plate(i)%n_i%nitrogen, nisep )
#endif
              case ('O')
#if IMAS_MINOR_VERSION > 34
                call write_sourced_value( summary%local%divertor_target(i)%n_i%oxygen, nisep )
#else
                call write_sourced_value( summary%local%divertor_plate(i)%n_i%oxygen, nisep )
#endif
              case ('Ne')
#if IMAS_MINOR_VERSION > 34
                call write_sourced_value( summary%local%divertor_target(i)%n_i%neon, nisep )
#else
                call write_sourced_value( summary%local%divertor_plate(i)%n_i%neon, nisep )
#endif
              case ('Ar')
#if IMAS_MINOR_VERSION > 34
                call write_sourced_value( summary%local%divertor_target(i)%n_i%argon, nisep )
#else
                call write_sourced_value( summary%local%divertor_plate(i)%n_i%argon, nisep )
#endif
#if IMAS_MINOR_VERSION > 30
              case ('Fe')
#if IMAS_MINOR_VERSION > 34
                call write_sourced_value( summary%local%divertor_target(i)%n_i%iron, nisep )
#else
                call write_sourced_value( summary%local%divertor_plate(i)%n_i%iron, nisep )
#endif
              case ('Kr')
#if IMAS_MINOR_VERSION > 34
                call write_sourced_value( summary%local%divertor_target(i)%n_i%krypton, nisep )
#else
                call write_sourced_value( summary%local%divertor_plate(i)%n_i%krypton, nisep )
#endif
#endif
              case ('Xe')
#if IMAS_MINOR_VERSION > 34
                call write_sourced_value( summary%local%divertor_target(i)%n_i%xenon, nisep )
#else
                call write_sourced_value( summary%local%divertor_plate(i)%n_i%xenon, nisep )
#endif
              case ('W')
#if IMAS_MINOR_VERSION > 34
                call write_sourced_value( summary%local%divertor_target(i)%n_i%tungsten, nisep )
#else
                call write_sourced_value( summary%local%divertor_plate(i)%n_i%tungsten, nisep )
#endif
              end select
            end do
#if IMAS_MINOR_VERSION > 34
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
#if IMAS_MINOR_VERSION > 31
            call write_sourced_value( summary%local%divertor_plate(i)%power_flux_peak, &
              & power_flux_peak(itrg(i)) )
#endif
#endif
          end do
        end if

        select case (GeometryType)
          case( GEOMETRY_LIMITER )
            call write_sourced_integer( summary%boundary%type, 0 )
          case( GEOMETRY_SN )
            if ( isymm.eq.1 .or. isymm.eq.2 ) then
              if ( cry(leftcut(1),jsep,3).lt.0.0_R8) then
                call write_sourced_integer( summary%boundary%type, 11 )
              else if ( cry(leftcut(1),jsep,3).gt.0.0_R8) then
                call write_sourced_integer( summary%boundary%type, 12 )
              end if
            else if ( isymm.eq.3 .or. isymm.eq.4 ) then
              if ( crx(leftcut(1),jsep,3).lt.0.0_R8) then
                call write_sourced_integer( summary%boundary%type, 11 )
              else if ( crx(leftcut(1),jsep,3).gt.0.0_R8) then
                call write_sourced_integer( summary%boundary%type, 12 )
              end if
            else
              call write_sourced_integer( summary%boundary%type, 1 )
            end if
          case( GEOMETRY_CDN , GEOMETRY_DDN_BOTTOM , GEOMETRY_DDN_TOP )
            call write_sourced_integer( summary%boundary%type, 13 )
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

        call write_sourced_value( summary%fusion%power, fusion_power*1.0e6_IDS_real )

        call write_sourced_int_constant( summary%gas_injection_rates%impurity_seeding, 0 )
        gsum = 0.0_R8
        gtop = 0.0_R8
        gmid = 0.0_R8
        gbot = 0.0_R8
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
#ifdef B25_EIRENE
        do istrai = 1, nstrai
          if (crcstra(istrai).eq.'C') then
            do iatm = 1, natmi
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
                     & (pfrregno2.eq.2 .and. GeometryType.eq.GEOMETRY_SN) .or. &
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
                  case (GEOMETRY_SN)
                    at_mid = bcchar(ib).eq.'N'
                    at_bot = bcchar(ib).eq.'S'.and.LSN.and.(ireg.eq.3.or.ireg.eq.4)
                    at_top = bcchar(ib).eq.'S'.and..not.LSN.and.(ireg.eq.3.or.ireg.eq.4)
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
              if (gpfc(iatm,istrai).eq.1.0_R8) then
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
              endif
            end if
          end if
        end do
#endif
        call write_sourced_value( summary%gas_injection_rates%total, gsum )
        call write_sourced_value( summary%gas_injection_rates%midplane, gmid )
        call write_sourced_value( summary%gas_injection_rates%top, gtop )
        call write_sourced_value( summary%gas_injection_rates%bottom, gbot )

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
          qemax = 0.0_IDS_real
          qimax = 0.0_IDS_real
          do iy = jsep+1, ny-1
            qemax = qemax + (-icnt)*fhe(ix,iy,0)
            qimax = qimax + (-icnt)*fhi(ix,iy,0)
          end do
          do i = ix+icnt, jxa, icnt
            qetot = 0.0_IDS_real
            qitot = 0.0_IDS_real
            do iy = jsep+1, ny-1
              qetot = qetot + (-icnt)*fhe(i,iy,0)
              qitot = qitot + (-icnt)*fhi(i,iy,0)
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
          totFace=abs(fhe)
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
          totFace=abs(fhi)
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

        contains

        integer function get_atom_number( compname )
            implicit none
            character*2 compname
            integer is, iatm
            logical streql
            external streql

            get_atom_number = 0
            do iatm = 1, natmi
               if ( get_atom_number > 0 ) cycle
               is = eb2atcr(iatm)
               if (streql( is_codes( is ), compname ) ) get_atom_number = iatm
            end do
            return

        end function get_atom_number

        function separatrix_average( field, weight )
        ! This function is devoted to obtain the weighted average along the active separatrix
        ! of a plasma field quantity
        ! The average is made using face-centered quantities on the cell faces forming the separatrix
        ! The weighting automatically includes the areas of the cell faces
        implicit none
        real(kind=IDS_real) :: separatrix_average
        real(kind=IDS_real), intent(in) :: field(nx,ny), weight(nx,ny)
        real(kind=IDS_real) :: sum, area_sum

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

        subroutine write_ids_properties( properties, homo, &
          & comment, source, create_date)
            type(ids_ids_properties), intent(inout) :: properties
                !< Type of IDS data structure, designed for IDS properties
            integer, intent(in) :: homo
            character(len=ids_string_length), intent(in) :: comment
            character(len=ids_string_length), intent(in) :: source
            character(len=ids_string_length), intent(in) :: create_date

            properties%homogeneous_time = homo
            allocate( properties%comment(1) )
            properties%comment = comment
#if ( IMAS_MINOR_VERSION > 33 && 0 )
            allocate( properties%provenance%node(1) )
            allocate( properties%provenance%node(1)%sources(1) )
            properties%provenance%node(1)%sources(1) = source
#else
            allocate( properties%source(1) )
            properties%source = source
#endif
            allocate( properties%creation_date(1) )
            properties%creation_date = create_date
#if IMAS_MINOR_VERSION > 14
            allocate( properties%provider(1) )
            properties%provider = username
#endif
#if IMAS_MINOR_VERSION > 21
            allocate( properties%version_put%data_dictionary(1) )
            properties%version_put%data_dictionary = imas_version
            allocate( properties%version_put%access_layer(1) )
            properties%version_put%access_layer = ual_version
            allocate( properties%version_put%access_layer_language(1) )
            properties%version_put%access_layer_language = 'FORTRAN'
#endif
            return

        end subroutine write_ids_properties

        subroutine write_ids_code( code, commit )
            type(ids_code), intent(inout) :: code
                !< Type of IDS data structure, designed for code data handling
            character(len=ids_string_length), intent(in) :: commit
#if IMAS_MINOR_VERSION > 29
            integer :: nlibs !< Number of declared libraries in IDS description
            character*8 ggd_version, mscl_version
            character*32 SOLPS_git_version
            character*32 get_SOLPS_hash
            character(len=ids_string_length) :: repository
#ifdef B25_EIRENE
            character*8 eirene_version
            character*31 Eirene_git_version
            character*31 get_Eir_hash
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
            allocate( code%version(1) )
            code%version = newversion
            allocate( code%commit(1) )
            code%commit = commit
            allocate( code%repository(1) )
            code%repository(1) = "ssh://git.iter.org/bnd/b2.5.git"
            allocate( code%output_flag( num_time_slices ) )
            code%output_flag( time_sind ) = 0

#if IMAS_MINOR_VERSION > 29
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
#ifdef NAGFOR
            call get_environment_variable('GGD_VERSION', &
              &  status=ierror,length=lenval)
            if (ierror.eq.0) call get_environment_variable('GGD_VERSION', &
              &  value=ggd_version)
            call get_environment_variable('EBVERSIONMSCL', &
              &  status=ierror,length=lenval)
            if (ierror.eq.0) call get_environment_variable('EBVERSIONMSCL', &
              &  value=mscl_version)
#else
#ifdef USE_PXFGETENV
            CALL PXFGETENV ('GGD_VERSION', 0, ggd_version, lenval, ierror)
            CALL PXFGETENV ('EBVERSIONMSCL', 0, mscl_version, lenval, ierror)
#else
            call getenv ('GGD_VERSION', ggd_version)
            call getenv ('EBVERSIONMSCL', mscl_version)
#endif
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

#if IMAS_MINOR_VERSION > 32
        subroutine write_ids_midplane( midplane, midplane_id )
            type(ids_identifier_static) :: midplane
            integer, intent(in) :: midplane_id

            midplane%index = midplane_id
            allocate( midplane%name(1) )
            allocate( midplane%description(1) )
            select case (midplane_id)
            case (1)
              midplane%name = 'magnetic_axis'
              midplane%description = &
                 &  'Height of equilibrium O-point'
            case (2)
              midplane%name = 'dr_dz_zero_sep'
              midplane%description = &
                 &  'Maximum radius location along separatrix'
            case (3)
              midplane%name = 'z_zero'
              midplane%description = &
                 &  'Z = 0 plane'
            case (4)
              midplane%name = 'ggd_subset'
              midplane%description = &
                 &  'Location specified by GGD outer midplane grid subset'
            end select
            return

        end subroutine write_ids_midplane
#endif

        subroutine write_sourced_constant( val, value )
            type(ids_summary_constant_flt_0d) :: val
                !< Type of IDS data structure, designed for sourced real constant data handling
            real(ids_real), intent(in) :: value

            val%value = value
            allocate( val%source(1) )
            val%source = source

            return

        end subroutine write_sourced_constant

        subroutine write_sourced_int_constant( ival, ivalue )
            type(ids_summary_constant_int_0d) :: ival
                !< Type of IDS data structure, designed for sourced integer constant data handling
            integer, intent(in) :: ivalue

            ival%value = ivalue
            allocate( ival%source(1) )
            ival%source = source

            return

        end subroutine write_sourced_int_constant

        subroutine write_sourced_integer( ival, ivalue )
            type(ids_summary_dynamic_int_1d_root) :: ival
                !< Type of IDS data structure, designed for sourced integer data handling
            integer, intent(in) :: ivalue

            allocate( ival%value( num_time_slices ) )
            ival%value( time_sind ) = ivalue
            allocate( ival%source(1) )
            ival%source = source

            return

        end subroutine write_sourced_integer

        subroutine write_sourced_string( val, string )
            type(ids_summary_static_str_0d) :: val
                !< Type of IDS data structure, designed for sourced string data handling
            character(len=ids_string_length), intent(in) :: string

            allocate( val%value(1) )
            val%value = string
            allocate( val%source(1) )
            val%source = source

            return

        end subroutine write_sourced_string

#if IMAS_MINOR_VERSION > 29
        subroutine write_timed_integer( ival, ivalue )
            type(ids_signal_int_1d), intent(inout) :: ival
                !< Type of IDS data structure, designed for integer data handling
            integer, intent(in) :: ivalue

            allocate( ival%data( num_time_slices ) )
            ival%data( time_sind ) = ivalue
            allocate( ival%time( num_time_slices ) )
            ival%time( time_sind ) = time_slice_value

            return

        end subroutine write_timed_integer
#endif

#if IMAS_MINOR_VERSION > 30
        subroutine write_timed_value( val, value )
            type(ids_signal_flt_1d), intent(inout) :: val
                !< Type of IDS data structure, designed for scalar data handling
            real(IDS_real), intent(in) :: value

            allocate( val%data( num_time_slices ) )
            val%data( time_sind ) = value
            allocate( val%time( num_time_slices ) )
            val%time( time_sind ) = time_slice_value

            return

        end subroutine write_timed_value
#endif

        subroutine add_sourced_value( val, value )
            type(ids_summary_dynamic_flt_1d_root_parent_2) :: val
                !< Type of IDS data structure, designed for sourced float data handling
            real(ids_real), intent(in) :: value

            if ( associated( val%value ) ) then
              val%value( time_sind ) = val%value( time_sind ) + value
            else
              call write_sourced_value( val, value )
            end if

            return

        end subroutine add_sourced_value

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
                !< handling data field values
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
            !! Interpolate data to vertices
            tmpVx = interpolateToVertices(  &
                  &   gmap%b2nx, gmap%b2ny, VX_LOWERLEFT, value )

            !! Interpolate data to cell faces, using a volume weighting
            tmpFace = 0.0_IDS_real
            do i = TO_SELF, TO_TOP
                weight(:,:,i) = vol(:,:)
            end do
            call value_on_faces( nx, ny, weight, value, tmpFace)

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
                  case( GRID_SUBSET_EDGES, &
                      & GRID_SUBSET_X_ALIGNED_EDGES, GRID_SUBSET_Y_ALIGNED_EDGES, &
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
#if IMAS_MINOR_VERSION < 15
                  idsdata => b2_IMAS_Transform_Data_B2_To_IDS_Vertex( &
                     &   edge_profiles%ggd( time_sind )%grid,         &
                     &   iGsOuterMidplane, gmap, tmpVx )
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
                  idsdata => b2_IMAS_Transform_Data_B2_To_IDS(             &
                      &  edge_profiles% ggd( time_sind )%grid, iSubset,    &
                      &  gmap, value )
#else
                  idsdata => b2_IMAS_Transform_Data_B2_To_IDS(             &
                      &  edge_profiles%grid_ggd( time_sind ), iSubset,     &
                      &  gmap, value )
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
               ndim = edge_profiles%grid_ggd(time_sind)%grid_subset(iSubset)%dimension
               iSubsetID = edge_profiles%grid_ggd(time_sind)%grid_subset(iSubset)%identifier%index
#endif
               if (ndim.eq.IDS_INT_INVALID) then
                  select case (iSubsetID)
                  case( GRID_SUBSET_NODES, GRID_SUBSET_X_POINTS, &
                      & GRID_SUBSET_INNER_MIDPLANE, GRID_SUBSET_OUTER_MIDPLANE )
                     ndim = 1
                  case( GRID_SUBSET_EDGES, &
                      & GRID_SUBSET_X_ALIGNED_EDGES, GRID_SUBSET_Y_ALIGNED_EDGES, &
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
                !< handling data field values
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
                  case( GRID_SUBSET_EDGES, &
                      & GRID_SUBSET_X_ALIGNED_EDGES, GRID_SUBSET_Y_ALIGNED_EDGES, &
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
        !! @note Available IDS vector component data fields (vector IDs):
        !!          - VEC_ALIGN_RADIAL_ID ( "radial" ),
        !!          - "diamagnetic",
        !!          - VEC_ALIGN_PARALLEL_ID ( "parallel" ),
        !!          - VEC_ALIGN_POLOIDAL_ID ( "poloidal" ),
        !!          - VEC_ALIGN_TOROIDAL_ID ( "toroidal" )
        subroutine write_cell_vector_component( vectorComponent, b2CellData,   &
                &   vectorID )
            type(ids_generic_grid_vector_components), intent(inout),    &
                &   pointer :: vectorComponent(:) !< Type of IDS data structure,
                    !> designed for vector data handling
            real(IDS_real), intent(in) :: b2CellData(-1:gmap%b2nx, -1:gmap%b2ny)
            real(IDS_real), dimension(:), pointer :: idsdata    !< Array for
                !< handling data field values
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
                  case( GRID_SUBSET_EDGES, &
                      & GRID_SUBSET_X_ALIGNED_EDGES, GRID_SUBSET_Y_ALIGNED_EDGES, &
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

        !> Write a vector component B2 face quantity to ids_generic_grid_vector
        !! components
        !! @note Available IDS vector component data fields (vector IDs):
        !!          - VEC_ALIGN_RADIAL_ID ( "radial" ),
        !!          - "diamagnetic",
        !!          - VEC_ALIGN_PARALLEL_ID ( "parallel" ),
        !!          - VEC_ALIGN_POLOIDAL_ID ( "poloidal" ),
        !!          - VEC_ALIGN_TOROIDAL_ID ( "toroidal" )
        subroutine write_face_vector_component( vectorComponent, b2FaceData,   &
                &   vectorID )
            type(ids_generic_grid_vector_components), intent(inout),    &
                &   pointer :: vectorComponent(:) !< Type of IDS data structure,
                    !> designed for vector data handling
            real(IDS_real), intent(in) ::  &
                &   b2FaceData(-1:gmap%b2nx, -1:gmap%b2ny, 0:1)
            real(IDS_real), dimension(:), pointer :: idsdata    !< Array for
                !< handling data field values
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
               ndim = 2
               iSubsetID = GRID_SUBSET_FACES
#else
               ndim = edge_profiles%grid_ggd(time_sind)%grid_subset(iSubset)%dimension
               iSubsetID = edge_profiles%grid_ggd(time_sind)%grid_subset(iSubset)%identifier%index
#endif
               if (ndim.eq.IDS_INT_INVALID) then
                  select case (iSubsetID)
                  case( GRID_SUBSET_NODES, GRID_SUBSET_X_POINTS, &
                      & GRID_SUBSET_INNER_MIDPLANE, GRID_SUBSET_OUTER_MIDPLANE )
                     ndim = 1
                  case( GRID_SUBSET_EDGES, &
                      & GRID_SUBSET_X_ALIGNED_EDGES, GRID_SUBSET_Y_ALIGNED_EDGES, &
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
#if IMAS_MINOR_VERSION < 15
               idsdata => b2_IMAS_Transform_Data_B2_To_IDS( edge_profiles% &
                   &   ggd( time_sind )%grid, iSubset,                     &
                   &   gmap, b2FaceData )
#else
               idsdata => b2_IMAS_Transform_Data_B2_To_IDS( edge_profiles% &
                   &   grid_ggd( time_sind ), iSubset,                     &
                   &   gmap, b2FaceData )
#endif
               call B2grid_Write_Data_Vector_Components( vectorComponent(iSubset), &
                   &   ggdID, iSubsetID, vectorID, idsdata )
               deallocate(idsdata)
            end do

        end subroutine write_face_vector_component

        !!$> TODO: add to GGD itself (ids_grid_data)!
        !> Write a scalar data field given as a scalar data representation to a
        !! generic grid vector component IDS data fields.
        !!
        !! @note    The routine will make sure the required storage is
        !!          allocated, and will deallocate and re-allocate fields as
        !!          necessary.
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
            integer, intent(in) :: time_sind    !< Time slice index

            if ( .not. present(gridSubsetInd) ) then
                !! Fill in vector component data
#if IMAS_MINOR_VERSION < 15
                idsdata => b2_IMAS_Transform_Data_B2_To_IDS(    &
                    &   edge_profiles%ggd( time_sind )%grid,    &
                    &   GRID_SUBSET_Y_ALIGNED_EDGES, gmap, b2FaceData)
#else
                idsdata => b2_IMAS_Transform_Data_B2_To_IDS(    &
                    &   edge_profiles%grid_ggd( time_sind ),    &
                    &   GRID_SUBSET_Y_ALIGNED_EDGES, gmap, b2FaceData)
#endif
#if GGD_MINOR_VERSION > 8
                call gridWriteData( vector, gridId, GRID_SUBSET_Y_ALIGNED_EDGES, idsdata )
#else
                call gridWriteData( vector, GRID_SUBSET_Y_ALIGNED_EDGES, idsdata )
#endif
                deallocate(idsdata)
#if IMAS_MINOR_VERSION < 15
                idsdata => b2_IMAS_Transform_Data_B2_To_IDS(    &
                    &   edge_profiles%ggd( time_sind )%grid,    &
                    &   GRID_SUBSET_X_ALIGNED_EDGES, gmap, b2FaceData)
#else
                idsdata => b2_IMAS_Transform_Data_B2_To_IDS(    &
                    &   edge_profiles%grid_ggd( time_sind ),    &
                    &   GRID_SUBSET_X_ALIGNED_EDGES, gmap, b2FaceData)
#endif
#if GGD_MINOR_VERSION > 8
                call gridWriteData( vector, gridId, GRID_SUBSET_X_ALIGNED_EDGES, idsdata )
#else
                call gridWriteData( vector, GRID_SUBSET_X_ALIGNED_EDGES, idsdata )
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

    subroutine write_sourced_value_root( val, value )
        type(ids_summary_dynamic_flt_1d_root) :: val
            !< Type of IDS data structure, designed for sourced float data handling
        real(ids_real), intent(in) :: value

        allocate( val%value( num_time_slices ) )
        val%value( time_sind ) = value
        allocate( val%source(1) )
        val%source = source

        return

    end subroutine write_sourced_value_root

    subroutine write_sourced_value_root_parent_2( val, value )
        type(ids_summary_dynamic_flt_1d_root_parent_2) :: val
            !< Type of IDS data structure, designed for sourced float data handling
        real(ids_real), intent(in) :: value

        allocate( val%value( num_time_slices ) )
        val%value( time_sind ) = value
        allocate( val%source(1) )
        val%source = source

        return

    end subroutine write_sourced_value_root_parent_2

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
      real(ITM_R8), intent(in) :: value(-1:gmap%b2nx, -1:gmap%b2ny)
      real(ITM_R8), intent(in) :: flux(-1:gmap%b2nx, -1:gmap%b2ny, 0:1)
      real(ITM_R8), dimension(:), pointer :: cpodata
      real(ITM_R8) :: weight(-1:gmap%b2nx, -1:gmap%b2ny, TO_SELF:TO_TOP)
      integer i

      allocate(values(5))
      cpodata => b2ITMTransformDataB2ToCPO( edgecpo%grid, B2_SUBGRID_CELLS, gmap, value )
      call gridWriteData( values(1), B2_SUBGRID_CELLS, cpodata )
      deallocate(cpodata)
      tmpFace = 0.0_ITM_R8
      do i = TO_SELF, TO_TOP
        weight(:,:,i)=vol(:,:)
      end do
      call value_on_faces(nx,ny,weight,value,tmpFace)
      cpodata => b2ITMTransformDataB2ToCPO( edgecpo%grid, iSgCore, gmap, tmpFace )
      call gridWriteData( values(2), iSgCore, cpodata )
      deallocate(cpodata)
      tmpVx = interpolateToVertices( gmap%b2nx, gmap%b2ny, VX_LOWERLEFT, value )
      cpodata => b2ITMTransformDataB2ToCPOVertex( edgecpo%grid, iSgInnerMidplane, gmap, tmpVx )
      call gridWriteData( values(3), iSgInnerMidplane, cpodata )
      deallocate(cpodata)
      cpodata => b2ITMTransformDataB2ToCPOVertex( edgecpo%grid, iSgOuterMidplane, gmap, tmpVx )
      call gridWriteData( values(4), iSgOuterMidplane, cpodata )
      deallocate(cpodata)
      cpodata => b2ITMTransformDataB2ToCPOVertex( edgecpo%grid, B2_SUBGRID_NODES, gmap, tmpVx )
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
      cpodata => b2ITMTransformDataB2ToCPO( edgecpo%grid, B2_SUBGRID_CELLS, gmap, b2CellData )
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
         cpodata => b2ITMTransformDataB2ToCPO(edgecpo%grid, B2_SUBGRID_CELLS, gmap, vecdata(:,:,i-1))
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
!!$          cpodata => b2ITMTransformDataB2ToCPO(edgecpo%grid, B2_SUBGRID_EDGES_Y, gmap, b2FaceData)
!!$          call gridWriteData( vector%comp(1), B2_SUBGRID_EDGES_Y, cpodata )
!!$          deallocate(cpodata)
!!$          cpodata => b2ITMTransformDataB2ToCPO(edgecpo%grid, B2_SUBGRID_EDGES_X, gmap, b2FaceData)
!!$          call gridWriteData( vector%comp(2), B2_SUBGRID_EDGES_X, cpodata )
!!$          deallocate(cpodata)
!!$      else
!!$          allocate(vector%comp(1))
!!$          allocate(vector%align(1))
!!$          allocate(vector%alignid(1))
!!$
!!$          vector%align(1) = GRID_UNDEFINED
!!$          vector%alignid(1) = ""
!!$
!!$          cpodata => b2ITMTransformDataB2ToCPO(edgecpo%grid, subgridInd, gmap, b2FaceData)
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
