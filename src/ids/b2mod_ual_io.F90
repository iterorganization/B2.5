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
    use b2mod_constants
    !use b2mod_residuals
    use b2mod_sources
    use b2mod_transport
    use b2mod_anomalous_transport
    use b2mod_boundary_namelist
    use b2mod_neutrals_namelist
    use b2mod_user_namelist
    use b2mod_indirect
    use b2mod_external
    use b2mod_interp
    use b2mod_ipmain
    use b2mod_b2cmfs
    use b2mod_version
    use b2mod_grid_mapping
#ifdef IMAS
    use b2mod_balance &
     & , only : read_balance
    use b2mod_b2plot &
     & , only : nxtl, nxtr, jxi, jxa, jsep
#endif
    use b2mod_b2plot_wall_loading
#ifdef B25_EIRENE
    use eirmod_ctrig
    use eirmod_cestim
    use eirmod_wneutrals
#ifdef IMAS
    use eirmod_cinit &
     & , only : fort_lc
    use eirmod_comusr &
     & , only : natmi, nmoli, nioni, nmassa, nchara, nmassm, ncharm, &
     &          nprt, nchrgi, nchari
    use b2mod_b2plot &
     & , only : triangle_vol, ix_e2b, wklng, alloc_b2mod_b2plot_eirene
#endif
#else
#ifdef IMAS
    use b2mod_b2plot &
     & , only : natmi
#endif
#endif
    use logging

#ifdef IMAS
    !! UAL Access
    use b2mod_ual_io_grid &
     & , only : INCLUDE_GHOST_CELLS
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

  public b25_process_ids
  integer, public :: num_time_slices  !< Total number of time slices.
  integer, public :: num_batch_slices !< Total number of batch slices.

  private

#ifdef IMAS

#if IMAS_MINOR_VERSION < 9
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
  logical, parameter :: B2_WRITE_DATA = .true.
  integer, save :: ixpos(4), ifpos(4), iypos(4) !< Target positions
  integer, save :: idir(4), iysep(4), ixmid(4)
  integer, save :: target_offset = 1
  integer, save :: nesepm_istra = -1
  integer, save :: pfrregno1 = 0
  integer, save :: pfrregno2 = 2
  integer, save :: use_eirene = 0
  integer, save :: boundary_namelist = 0
  real(IDS_real), save :: ndes = 0.0_IDS_real
  real(IDS_real), save :: ndes_sol = 0.0_IDS_real
  real(IDS_real), save :: nesepm_pfr = 0.0_IDS_real
  real(IDS_real), save :: nesepm_sol = 0.0_IDS_real
  real(IDS_real), save :: nepedm_sol = 0.0_IDS_real
  real(IDS_real), save :: volrec_sol = 0.0_IDS_real
  real(IDS_real), save :: private_flux_puff = 0.0_IDS_real
  real(IDS_real), save :: b0, r0, b0r0
  real(IDS_real), save :: flux_expansion(4), extension_r(4), extension_z(4)
  character(len=ids_string_length), save :: username  !< IDS user name
  character(len=ids_string_length), save :: source    !< Code source
  character(len=ids_string_length), save :: comment   !< IDS properties label
  character(len=ids_string_length), save :: create_date
  character(len=ids_string_length), save :: code_commit
  character(len=ids_string_length), save :: configuration
  character(len=ids_string_length), save :: plate_name(4) !< Divertor plate name
  character*8, save :: imas_version, ual_version, adas_version
  character*8, save :: date
  character*10, save :: ctime
  character*5, save :: zone
  character*32, save :: B25_git_version
  character*32, save :: ADAS_git_version
  logical, save :: IDSmapInitialized = .false.
  logical, save :: eq_found
#ifdef USE_PXFGETENV
  integer lenval, ierror
#else
#ifdef NAGFOR
  integer lenval, ierror
#endif
#endif
  type(B2GridMap), save :: IDSmap

contains

    subroutine IDS_init
    implicit none
    integer tvalues(8)
    integer i, istrai, p
    integer itrg(4)
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
    call ipgeti ('b2stbc_boundary_namelist', boundary_namelist)
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
        implicit none
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
        integer :: time_sind !< Time slice index
        integer :: ns   !< Total number of B2.5 species
        integer :: nx   !< Specifies the number of interior cells
                        !< along the first coordinate
        integer :: ny   !< Specifies the number of interior cells
                        !< along the second coordinate
        integer :: n_process !< Number of radiation processes handled
        integer :: is   !< Species index (iterator)
        integer :: i    !< Iterator
        integer :: j    !< Iterator
        integer :: k      !< Iterator
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
        logical, parameter :: B2_WRITE_DATA = .true.
        integer :: nscx, iscx(0:nscxmax-1)
        real(IDS_real),   &
            &   dimension( -1:ubound( crx, 1 ), -1:ubound( crx, 2), 3, 3) :: e
        real(IDS_real) :: tmpCv( -1:ubound( na, 1), -1:ubound( na, 2) )
        real(IDS_real) :: lnlam( -1:ubound( na, 1), -1:ubound( na, 2) )
        real(IDS_real) :: zeff( -1:ubound( na, 1), -1:ubound( na, 2) )
        real(IDS_real) :: time_step !< Time step
        real(IDS_real) :: time_slice_value   !< Time slice value
        real(IDS_real) :: time  !< Generic time
        type(B2GridMap) :: IDSmap !< Data structure holding an
            !< intermediate grid description to be transferred into a CPO or IDS

        integer, parameter :: nsources = 12
        integer, save :: ncall = 0
        integer, save :: style = 1
        integer, save :: ismain = 1
        integer, save :: ismain0 = 0
        integer, save :: ue_style = 2
        integer, save :: ids_from_43 = 0
        integer, save :: balance_netcdf = 0
        integer, save :: drift_style
        real(IDS_real), save :: dtim = 1.0_IDS_real
        real(IDS_real), save :: neutral_sources_rescale = 1.0_IDS_real
        real(IDS_real), save :: BoRiS = 0.0_IDS_real
        character*132 radiation_commit
        character*256 filename
        logical match_found, streql, exists, wrong_flow
        external streql

        !! ===  SET UP IDS ===
        write(0,*) "Setting data for edge_profiles IDS"
        if (ncall.eq.0) then
          call ipgetr ('b2news_BoRiS', BoRiS)
          call ipgeti ('b2mndt_style', style)
          call ipgeti ('b2mndr_ismain', ismain)
          call ipgeti ('b2sigp_style', ue_style)
          call ipgeti ('ids_from_43', ids_from_43)
          call ipgetr ('b2mndr_dtim', dtim)
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
            if (is_neutral(is) .and. zn(is).eq.zn(ismain) &
                             & .and. am(is).eq.am(ismain)) then
              ismain0 = is
              match_found = .true.
            end if
            is = is - 1
          end do
          if (.not.match_found.and.ismain.ne.1) ismain0 = ismain
          call ipgeti ('b2mwti_ismain0', ismain0)
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
            &        fch_p, na, ua, te, ti, po, ne, ni, ne2, chvemx, chvimx,      &
            &        cdna, cdpa, cddi, cvla, cvsa, chce, chve, chci,              &
            &        chvi, csig, csigin, calf, cthe, cthi,                        &
            &        cdnahz, cdpahz, cvlahz, cvmahz, cvsahz, cvsa_cl, cvsahz_cl,  &
            &        fllim0fhi, fllimvisc, csig_cl, calf_cl)
!  ..compute log-log charge exchange rate coefficients
        do k = 0, nscx-1
           call b2spcx (nx, ny, ns, ev, am(iscx(k)), ti, ne, rlcx(-1,-1,0,0,k))
        enddo
!   ..compute sources
        call b2sral (nx, ny, ns, nxtl, nxtr,                                     &
            &        nscx, nscxmax, 0, ns, iscx, ismain, ismain0,                &
            &        dtim, BoRiS, facdrift, fac_ExB, fac_vis,                    &
            &        vol, crx, hx, hy, hz, qz, qc, qs, gs, pbs, pbshz, bb,       &
            &        na, ua,                                                     &
            &        uadia, vedia, vadia, wadia, veecrb, vaecrb, ve, wedia,      &
            &        te, ti, po, ne, ni, kinrgy, floe_noc, floi_noc,             &
            &        fna, fna_32, fna_52, fni_32, fni_52, fne_32, fne_52,        &
            &        fna_mdf, fhe_mdf, fhi_mdf, fna_fcor, fna_nodrift, fna_he,   &
            &        fhe, fhi, fnaPSch, fhePSch, fhiPSch, fch,                   &
            &        fchdia, fchin, fch_p, fchvispar, fchvisper, fchvisq,        &
            &        fchinert, fchanml, fna_eir, fne_eir, fhe_eir, fhi_eir,      &
            &        cdna, cdpa, cvsa_cl,                                        &
            &        cvla, chce, chve, chci, chvi, calf,                         &
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
              close(46)
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
        call write_ids_properties( radiation%ids_properties, &
          &  homogeneous_time )
        call write_ids_properties( description%ids_properties, &
          &  homogeneous_time )
#if IMAS_MINOR_VERSION > 21
        call write_ids_properties( summary%ids_properties, &
          &  homogeneous_time )
#endif
#if ( IMAS_MINOR_VERSION > 25 && IMAS_MINOR_VERSION < 34 )
        call write_ids_properties( numerics%ids_properties, &
          &  homogeneous_time )
#endif
#if IMAS_MINOR_VERSION > 30
        call write_ids_properties( divertors%ids_properties, &
          &  homogeneous_time )
#endif

        !! 2. Set code and library data
        if (streql(b2frates_flag,'adas')) then
          radiation_commit = 'B25 : '//trim(B25_git_version)// &
                      &  ' + ADAS : '//trim(ADAS_git_version)
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
        edge_sources%source(3)%identifier%description = &
            & "Boundary conditions sources from "//trim(source)
        !! Time derivatives
        edge_sources%source(4)%identifier%index = 706
        allocate( edge_sources%source(4)%identifier%name(1) )
        edge_sources%source(4)%identifier%name = "Time derivative"
        allocate( edge_sources%source(4)%identifier%description(1) )
        edge_sources%source(4)%identifier%description = &
            & "Time derivative sources from "//trim(source)
        !! Atomic ionization
        edge_sources%source(5)%identifier%index = 707
        allocate( edge_sources%source(5)%identifier%name(1) )
        edge_sources%source(5)%identifier%name = "Atomic ionization"
        allocate( edge_sources%source(5)%identifier%description(1) )
        edge_sources%source(5)%identifier%description = &
            & "Atomic ionization sources from "//trim(source)
        !! Molecular ionization
        edge_sources%source(6)%identifier%index = 708
        allocate( edge_sources%source(6)%identifier%name(1) )
        edge_sources%source(6)%identifier%name = "Molecular ionization"
        allocate( edge_sources%source(6)%identifier%description(1) )
        edge_sources%source(6)%identifier%description = &
            & "Molecular ionization sources from "//trim(source)
        !! Ionization
        edge_sources%source(7)%identifier%index = 709
        allocate( edge_sources%source(7)%identifier%name(1) )
        edge_sources%source(7)%identifier%name = "Ionization"
        allocate( edge_sources%source(7)%identifier%description(1) )
        edge_sources%source(7)%identifier%description = &
            & "Ionization sources from "//trim(source)
        !! Recombination
        edge_sources%source(8)%identifier%index = 710
        allocate( edge_sources%source(8)%identifier%name(1) )
        edge_sources%source(8)%identifier%name = "Recombination"
        allocate( edge_sources%source(8)%identifier%description(1) )
        edge_sources%source(8)%identifier%description = &
            & "Recombination sources from "//trim(source)
        !! Charge exchange
        edge_sources%source(9)%identifier%index = 305
        allocate( edge_sources%source(9)%identifier%name(1) )
        edge_sources%source(9)%identifier%name = "Charge exchange"
        allocate( edge_sources%source(9)%identifier%description(1) )
        edge_sources%source(9)%identifier%description = &
            & "Charge exchange sources from "//trim(source)
        !! Collisional equipartition
        edge_sources%source(10)%identifier%index = 11
        allocate( edge_sources%source(10)%identifier%name(1) )
        edge_sources%source(10)%identifier%name = "Equipartition"
        allocate( edge_sources%source(10)%identifier%description(1) )
        edge_sources%source(10)%identifier%description = &
            & "Collisional equipartition sources from "//trim(source)
        !! Ohmic
        edge_sources%source(11)%identifier%index = 7
        allocate( edge_sources%source(11)%identifier%name(1) )
        edge_sources%source(11)%identifier%name = "Ohmic"
        allocate( edge_sources%source(11)%identifier%description(1) )
        edge_sources%source(11)%identifier%description = &
            & "Ohmic (Joule) sources from "//trim(source)
        !! Radiation
        edge_sources%source(12)%identifier%index = 200
        allocate( edge_sources%source(12)%identifier%name(1) )
        edge_sources%source(12)%identifier%name = "Radiation"
        allocate( edge_sources%source(12)%identifier%description(1) )
        edge_sources%source(12)%identifier%description = &
            & "Radiation sources from "//trim(source)

        call put_equilibrium_data ( equilibrium, &
#if IMAS_MINOR_VERSION > 21
            &  summary, &
#endif
            &  edge_profiles, database, .true. )
        allocate( radiation%vacuum_toroidal_field%b0( num_time_slices ) )
        radiation%vacuum_toroidal_field%b0( time_sind ) = &
            &  edge_profiles%vacuum_toroidal_field%b0( time_sind )
        radiation%vacuum_toroidal_field%r0 = &
            &  edge_profiles%vacuum_toroidal_field%r0
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
        if ( present( time_step_IN ) ) &
          &  description%simulation%time_step = time_step_IN
        if ( present ( time_IN ) ) &
          &  description%simulation%time_current = time_IN
        allocate( description%simulation%workflow(1) )
        description%simulation%workflow = source
#if ( IMAS_MINOR_VERSION > 25 && IMAS_MINOR_VERSION < 34 )
        description%simulation%time_begin = run_start_time_IN
        description%simulation%time_end = run_end_time_IN
#endif

        i=index(B25_git_version,'-')
        allocate( summary%tag%name(1) )
        summary%tag%name = B25_git_version(1:i-1)

#if IMAS_MINOR_VERSION > 32
        call write_ids_midplane( divertors%midplane, midplane_id )
        call write_ids_midplane( edge_profiles%midplane, midplane_id )
        call write_ids_midplane( edge_sources%midplane, midplane_id )
        call write_ids_midplane( edge_transport%midplane, midplane_id )
        call write_ids_midplane( summary%midplane, midplane_id )
#endif

        ns = size( na, 3 )
        nx = ubound( na, 1 )
        ny = ubound( na, 2 )

        !! List of species
        allocate( edge_profiles%ggd( time_sind )%ion( ns ) )
        do is = 0, ns-1
            allocate( edge_profiles%ggd( time_sind )%ion( is + 1 )%state(1) )
            allocate( edge_profiles%ggd( time_sind )%ion( is + 1 )%element(1) )

            call species( is, edge_profiles%ggd( time_sind )%ion( is + 1 )% &
                &   state(1)%label, .false.)
            edge_profiles%ggd( time_sind )%ion( is + 1 )%element(1)%a = am( is )
            edge_profiles%ggd( time_sind )%ion( is + 1 )%element(1)%z_n =   &
                &   zn( is )
            edge_profiles%ggd( time_sind )%ion( is + 1 )%state(1)%z_min =   &
                &   zamin( is )
            edge_profiles%ggd( time_sind )%ion( is + 1 )%state(1)%z_max =   &
                &   zamax( is )
        enddo

        write(*,*) "Running b2CreateMap subroutine"
        !! Set up the B2<->IDS mappings
        call b2CreateMap( nx, ny, crx( -1:nx, -1:ny, : ),             &
            &   cry( -1:nx, -1:ny, : ), cflags, leftix, leftiy,       &
            &   rightix, rightiy, topix, topiy, bottomix,bottomiy,    &
            &   INCLUDE_GHOST_CELLS, IDSmap, .false. )
        mapInitialized = .true.

        !! Write grid & grid subsets/subgrids
#if IMAS_MINOR_VERSION > 11
#if IMAS_MINOR_VERSION < 15
        call b2_IMAS_Fill_Grid_Desc( IDSmap,                                &
            &   edge_profiles%ggd( time_sind )%grid,                        &
            &   nx, ny, crx(-1:nx, -1:ny, :), cry(-1:nx, -1:ny, : ),        &
            &   leftix, leftiy, rightix, rightiy, topix, topiy, bottomix,   &
            &   bottomiy, nnreg, topcut, region, cflags,                    &
            &   INCLUDE_GHOST_CELLS, vol, gs, qc )
        call b2_IMAS_Fill_Grid_Desc( IDSmap,                                &
            &   edge_transport%model(1)%ggd( time_sind )%grid,              &
            &   nx, ny, crx(-1:nx, -1:ny, :), cry(-1:nx, -1:ny, : ),        &
            &   leftix, leftiy, rightix, rightiy, topix, topiy, bottomix,   &
            &   bottomiy, nnreg, topcut, region, cflags,                    &
            &   INCLUDE_GHOST_CELLS, vol, gs, qc )
        do is = 1, nsources
            call b2_IMAS_Fill_Grid_Desc( IDSmap,                              &
                &   edge_sources%source(is)%ggd( time_sind )%grid,            &
                &   nx, ny, crx(-1:nx, -1:ny, :), cry(-1:nx, -1:ny, : ),      &
                &   leftix, leftiy, rightix, rightiy, topix, topiy, bottomix, &
                &   bottomiy, nnreg, topcut, region, cflags,                  &
                &   INCLUDE_GHOST_CELLS, vol, gs, qc )
        end do
#else
        call b2_IMAS_Fill_Grid_Desc( IDSmap,                                &
            &   edge_profiles%grid_ggd( time_sind ),                        &
            &   nx, ny, crx(-1:nx, -1:ny, :), cry(-1:nx, -1:ny, : ),        &
            &   leftix, leftiy, rightix, rightiy, topix, topiy, bottomix,   &
            &   bottomiy, nnreg, topcut, region, cflags,                    &
            &   INCLUDE_GHOST_CELLS, vol, gs, qc )
        call b2_IMAS_Fill_Grid_Desc( IDSmap,                                &
            &   edge_transport%grid_ggd( time_sind ),                       &
            &   nx, ny, crx(-1:nx, -1:ny, :), cry(-1:nx, -1:ny, : ),        &
            &   leftix, leftiy, rightix, rightiy, topix, topiy, bottomix,   &
            &   bottomiy, nnreg, topcut, region, cflags,                    &
            &   INCLUDE_GHOST_CELLS, vol, gs, qc )
        call b2_IMAS_Fill_Grid_Desc( IDSmap,                                &
            &   edge_sources%grid_ggd( time_sind ),                         &
            &   nx, ny, crx(-1:nx, -1:ny, :), cry(-1:nx, -1:ny, : ),        &
            &   leftix, leftiy, rightix, rightiy, topix, topiy, bottomix,   &
            &   bottomiy, nnreg, topcut, region, cflags,                    &
            &   INCLUDE_GHOST_CELLS, vol, gs, qc )
#if IMAS_MINOR_VERSION > 21
        call b2_IMAS_Fill_Grid_Desc( IDSmap,                                &
            &   radiation%grid_ggd( time_sind ),                            &
            &   nx, ny, crx(-1:nx, -1:ny, :), cry(-1:nx, -1:ny, : ),        &
            &   leftix, leftiy, rightix, rightiy, topix, topiy, bottomix,   &
            &   bottomiy, nnreg, topcut, region, cflags,                    &
            &   INCLUDE_GHOST_CELLS, vol, gs, qc )
#endif
#endif
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

            !! ne: Electron Density
            call write_quantity( edge_profiles,                             &
                &   val = edge_profiles%ggd( time_sind )%electrons%density, &
                &   value = ne )
            !! sne: Electron particle sources
            tmpCv(:,:) = ( sne(:,:,0) + sne(:,:,1) * ne(:,:) ) / vol(:,:)
            call write_cell_scalar( edge_profiles,                          &
                &   scalar = edge_sources%source(1)%ggd( time_sind )%       &
                &            electrons%particles,                           &
                &   b2CellData = tmpCv )

            !! ni (SOLPS 4.x) /
            !! na (SOLPS 5.x): Ion Density
            allocate( edge_profiles%ggd( time_sind )%ion( ns ) )
            allocate( edge_transport%model(1)%ggd( time_sind )%ion( ns ) )
            allocate( edge_sources%source(1)%ggd( time_sind )%ion( ns ) )

            do is = 1, ns
                call write_quantity( edge_profiles,                             &
                    &   val = edge_profiles%ggd( time_sind )%ion(is)%density,   &
                    &   value = na(:,:, is - 1 ) )
               tmpCv(:,:) = ( sna(:,:,0,is - 1 ) + &
                    &         sna(:,:,1,is - 1 ) * na(:,:, is - 1) ) / vol(:,:)
               call write_cell_scalar( edge_profiles,                           &
                    &   scalar = edge_sources%source(1)%ggd( time_sind )%       &
                    &   ion( is )%particles,                                    &
                    &   b2CellData = tmpCv )
            end do


            !! ua: Parallel Ion Velocity
            do is = 1, ns
                call write_cell_vector_component( edge_profiles,            &
                    &   vectorComponent = edge_profiles%ggd( time_sind )%   &
                    &                     ion( is )%velocity,               &
                    &   b2CellData = ua(:,:, is - 1 ),                      &
                    &   vectorID = VEC_ALIGN_PARALLEL_ID )
            end do

            !! te: Electron Temperature
            tmpCv(:,:) = te(:,:)/qe
            call write_quantity( edge_profiles,                         &
                &   val = edge_profiles%ggd( time_sind )%electrons%     &
                &         temperature,                                  &
                &   value = tmpCv )

            !! ti: Ion Temperature
            tmpCv(:,:) = ti(:,:)/qe
            call write_quantity( edge_profiles,                         &
                &   val = edge_profiles%ggd( time_sind )%t_i_average,   &
                &   value = tmpCv )

            !! po: Electric Potential
            call write_quantity( edge_profiles,                         &
                &   val = edge_profiles%ggd( time_sind )%phi_potential, &
                &   value = po )

            !! B (magnetic field vector)
            !! Compute unit basis vectors along the field directions
            call compute_Coordinate_Unit_Vectors(crx, cry, e(:,:,:,1), &
                &   e(:,:,:,2), e(:,:,:,3))

        end if

        call logmsg( LOGDEBUG, "b2mod_ual_io.B25_process_ids: done" )

        contains

#endif

#if IMAS_MINOR_VERSION > 29
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
#endif
#if 0
        !> Write a vector B2 cell quantity to a complexgrid_vector
        subroutine write_cell_vector( vector, align, alignid, vecdata )
            type(ids_generic_grid_vector), intent(inout) :: vector
            real(IDS_real), intent(in) :: vecdata(-1:IDSmap%b2nx, &
                &   -1:IDSmap%b2ny, 0:2)
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

            ggdId = edge_profiles%grid_ggd(slice_index)%identifier%index
            !! Fill in vector component data
            do i = 1, dim
#if GGD_MINOR_VERSION < 9
                idsdata => b2_IMAS_Transform_Data_B2_To_IDS(        &
                    &   edge_profiles%ggd( slice_index )%grid,      &
                    &   GRID_SUBSET_CELLS, IDSmap, vecdata(:,:,i-1))
#else
                idsdata => b2_IMAS_Transform_Data_B2_To_IDS(        &
                    &   edge_profiles%grid_ggd( slice_index ),      &
                    &   GRID_SUBSET_CELLS, IDSmap, vecdata(:,:,i-1))
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
    end subroutine B25_process_ids

    subroutine write_ids_properties( properties, homo )
    implicit none
    type(ids_ids_properties), intent(inout) :: properties
                !< Type of IDS data structure, designed for IDS properties
    integer, intent(in) :: homo

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
    implicit none
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
    integer p
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
    allocate( code%output_flag( num_slices ) )
    code%output_flag( slice_index ) = 0

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

    subroutine put_equilibrium_data ( equilibrium, &
#if IMAS_MINOR_VERSION > 21
       &  summary, &
#endif
       &  edgeprof, database, do_summary_data )
    implicit none
    type (ids_equilibrium) :: equilibrium !< IDS designed to store
            !< equilibrium data
#if IMAS_MINOR_VERSION > 21
    type (ids_summary) :: summary !< IDS designed to store
            !< run summary data
#endif
    type (ids_edge_profiles) :: edgeprof !< IDS designed to store
            !< edge profiles data
    character(len=24), intent(in) :: database
    logical, intent(in) :: do_summary_data
    integer :: i, ix, icnt
    integer :: idum(0:3)
    integer, save :: ncall = 0
    real(IDS_real) :: parg(0:99)
    real(IDS_real), save :: pit_rescale = 1.0_IDS_real
    real(IDS_real) :: b0r0_ref, z_eq
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
    eq_source = "ITER Baseline q95=3 equilibrium"
    if ( associated( equilibrium%time_slice ) ) then
      if ( size( equilibrium%time_slice ).ge.slice_index ) then
        eq_found = .true.
        if ( associated( equilibrium%ids_properties%source ) ) &
           & eq_source = equilibrium%ids_properties%source(1)
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
#if IMAS_MINOR_VERSION > 21
            if (do_summary_data) then
              call write_sourced_value( summary%global_quantities%ip, &
                  & -15.0e6_IDS_real )
              call write_sourced_value( summary%global_quantities%b0, &
                  & -5.3_IDS_real )
            endif
#endif
            edgeprof%vacuum_toroidal_field%b0( slice_index ) = -5.3_IDS_real
          case (2)
#if IMAS_MINOR_VERSION > 21
            if (do_summary_data) then
              call write_sourced_value( summary%global_quantities%ip, &
                  & -7.5e6_IDS_real )
              call write_sourced_value( summary%global_quantities%b0, &
                  & -2.65_IDS_real )
            end if
#endif
            edgeprof%vacuum_toroidal_field%b0( slice_index ) = -2.65_IDS_real
          case (3)
#if IMAS_MINOR_VERSION > 21
            if (do_summary_data) then
              call write_sourced_value( summary%global_quantities%ip, &
                  & -5.0e6_IDS_real )
              call write_sourced_value( summary%global_quantities%b0, &
                  & -1.8_IDS_real )
            endif
#endif
            edgeprof%vacuum_toroidal_field%b0( slice_index ) = -1.8_IDS_real
          case default
#if IMAS_MINOR_VERSION > 21
            if (do_summary_data) then
              call write_sourced_value( summary%global_quantities%ip, &
                  & -15.0e6_IDS_real/nint(b0r0_ref/b0r0) )
              call write_sourced_value( summary%global_quantities%b0, &
                  & -b0r0 / 6.2e6_IDS_real )
            endif
#endif
            edgeprof%vacuum_toroidal_field%b0( slice_index ) = -b0r0 / 6.2_IDS_real
          end select
#if IMAS_MINOR_VERSION > 21
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
#if IMAS_MINOR_VERSION > 21
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
#if IMAS_MINOR_VERSION > 21
        if (do_summary_data) then
          call write_sourced_constant( summary%global_quantities%r0, &
                  & 6.2_IDS_real )
          summary%global_quantities%r0%source = &
                  & "ITER Baseline q95=3 equilibrium"
        end if
#endif
        edgeprof%vacuum_toroidal_field%r0 = 6.2_IDS_real
      else
        if (eq_found) then
          allocate( edgeprof%vacuum_toroidal_field%b0( num_slices ) )
          edgeprof%vacuum_toroidal_field%b0( slice_index ) = b0
          edgeprof%vacuum_toroidal_field%r0 = r0
#if IMAS_MINOR_VERSION > 21
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
        else if (isymm.ne.0) then
#if IMAS_MINOR_VERSION > 21
          if (do_summary_data) &
            & call write_sourced_value( summary%global_quantities%b0, -b0 )
#endif
          edgeprof%vacuum_toroidal_field%b0( slice_index ) = -b0
        else
#if IMAS_MINOR_VERSION > 21
          if (do_summary_data) &
            & call write_sourced_value( summary%global_quantities%b0, -b0 )
#endif
          edgeprof%vacuum_toroidal_field%b0( slice_index ) = b0
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

    ncall = ncall + 1
    return
    end subroutine put_equilibrium_data

#if IMAS_MINOR_VERSION > 32
    subroutine write_ids_midplane( midplane, midplane_id )
    implicit none
    type(ids_identifier_static) :: midplane
    integer, intent(in) :: midplane_id

    midplane%index = midplane_id
    allocate( midplane%name(1) )
    allocate( midplane%description(1) )
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
    return

    end subroutine write_ids_midplane
#endif

#if IMAS_MINOR_VERSION > 11
    !> Write scalar B2 cell quantity to 'ids_generic_grid_scalar'
    !! IMAS IDS data tree node.
    subroutine write_quantity( edgeprof, val, value )
    use b2mod_interp
    implicit none
    type(ids_edge_profiles), intent(inout) :: edgeprof
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

#if IMAS_MINOR_VERSION < 15
    ggdId = edgeprof%ggd(slice_index)%grid%identifier%index
    !! Assign 5+4 grid subsets
    nSubsets = 9
#else
    ggdId = edgeprof%grid_ggd(slice_index)%identifier%index
    nSubsets = size(edgeprof%grid_ggd(slice_index)%grid_subset)
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
      ndim = edgeprof%grid_ggd(slice_index)%grid_subset(iSubset)%dimension
      iSubsetID = edgeprof%grid_ggd(slice_index)%grid_subset(iSubset)%identifier%index
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
        idsdata => b2_IMAS_Transform_Data_B2_To_IDS_Vertex(        &
                     &   edgeprof%ggd( slice_index )%grid,         &
                     &   iGsOuterMidplane, IDSmap, tmpVx )
#else
        idsdata => b2_IMAS_Transform_Data_B2_To_IDS_Vertex(        &
                     &   edgeprof%grid_ggd( slice_index ),         &
                     &   iSubset, IDSmap, tmpVx )
#endif
#if GGD_MINOR_VERSION > 8
        call gridWriteData( val( iSubset ), ggdID, iSubsetID, idsdata )
#else
        call gridWriteData( val( iSubset ), iSubsetID, idsdata )
#endif
        deallocate( idsdata )
      case ( 2 ) !< Grid subset consists of faces
#if IMAS_MINOR_VERSION < 15
        idsdata => b2_IMAS_Transform_Data_B2_To_IDS(                &
                     &   edgeprof%ggd( slice_index )%grid, iSubset, &
                     &   IDSmap, tmpFace )
#else
        idsdata => b2_IMAS_Transform_Data_B2_To_IDS(                &
                     &   edgeprof%grid_ggd( slice_index ), iSubset, &
                     &   IDSmap, tmpFace )
#endif
#if GGD_MINOR_VERSION > 8
        call gridWriteData( val( iSubset ), ggdID, iSubsetID, idsdata )
#else
        call gridWriteData( val( iSubset ), iSubsetID, idsdata )
#endif
        deallocate( idsdata )
      case ( 3 ) !< Grid subset consists of cells
#if IMAS_MINOR_VERSION < 15
        idsdata => b2_IMAS_Transform_Data_B2_To_IDS(                 &
                      &  edgeprof% ggd( slice_index )%grid, iSubset, &
                      &  IDSmap, value )
#else
        idsdata => b2_IMAS_Transform_Data_B2_To_IDS(                 &
                      &  edgeprof%grid_ggd( slice_index ), iSubset,  &
                      &  IDSmap, value )
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

    return
    end subroutine write_quantity

    !> Write a scalar B2 cell quantity to ids_generic_grid_scalar
    subroutine write_cell_scalar( edgeprof, scalar, b2CellData )
    implicit none
    type (ids_edge_profiles), intent(in) :: edgeprof !< IDS designed to store
        !< edge profiles data
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

#if IMAS_MINOR_VERSION < 15
    ggdId = edgeprof%ggd(slice_index)%grid%identifier%index
    nSubsets = 1
#else
    ggdId = edgeprof%grid_ggd(slice_index)%identifier%index
    nSubsets = size(edgeprof%grid_ggd(slice_index)%grid_subset)
#endif
    !! Allocate data fields for grid subsets
    allocate( scalar(nSubsets) )

    do iSubset = 1, nSubsets
#if IMAS_MINOR_VERSION < 15
       ndim = 3
       iSubsetID = GRID_SUBSET_CELLS
#else
       ndim = edgeprof%grid_ggd(slice_index)%grid_subset(iSubset)%dimension
       iSubsetID = edgeprof%grid_ggd(slice_index)%grid_subset(iSubset)%identifier%index
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
       idsdata => b2_IMAS_Transform_Data_B2_To_IDS( edgeprof% &
          &   ggd( slice_index )%grid, iSubset, IDSmap, b2CellData )
#else
       idsdata => b2_IMAS_Transform_Data_B2_To_IDS( edgeprof% &
          &   grid_ggd( slice_index ), iSubset, IDSmap, b2CellData )
#endif
#if GGD_MINOR_VERSION > 8
       call gridWriteData( scalar( iSubset ), ggdID, iSubsetID, idsdata )
#else
       call gridWriteData( scalar( iSubset ), iSubsetID, idsdata )
#endif
       deallocate(idsdata)
    end do

    return
    end subroutine write_cell_scalar

    !> Write a vector component B2 cell quantity to ids_generic_grid_vector
    !! components
    !! @note Available IDS vector component data fields (vector IDs):
    !!          - VEC_ALIGN_RADIAL_ID ( "radial" ),
    !!          - "diamagnetic",
    !!          - VEC_ALIGN_PARALLEL_ID ( "parallel" ),
    !!          - VEC_ALIGN_POLOIDAL_ID ( "poloidal" ),
    !!          - VEC_ALIGN_TOROIDAL_ID ( "toroidal" )
    subroutine write_cell_vector_component( edgeprof, &
       &  vectorComponent, b2CellData, vectorID )
    implicit none
    type(ids_edge_profiles), intent(inout) :: edgeprof
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

#if IMAS_MINOR_VERSION < 15
    ggdId = edgeprof%ggd(slice_index)%grid%identifier%index
    nSubsets = 1
#else
    ggdId = edgeprof%grid_ggd(slice_index)%identifier%index
    nSubsets = size(edgeprof%grid_ggd(slice_index)%grid_subset)
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
      ndim = edgeprof%grid_ggd(slice_index)%grid_subset(iSubset)%dimension
      iSubsetID = edgeprof%grid_ggd(slice_index)%grid_subset(iSubset)%identifier%index
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
      idsdata => b2_IMAS_Transform_Data_B2_To_IDS( edgeprof% &
                   &   ggd( slice_index )%grid, iSubset,     &
                   &   IDSmap, b2CellData )
#else
      idsdata => b2_IMAS_Transform_Data_B2_To_IDS( edgeprof% &
                   &   grid_ggd( slice_index ), iSubset,     &
                   &   IDSmap, b2CellData )
#endif
      call B2grid_Write_Data_Vector_Components( vectorComponent(iSubset), &
          &   ggdID, iSubsetID, vectorID, idsdata )
      deallocate(idsdata)
    end do

    return
    end subroutine write_cell_vector_component


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
    subroutine write_face_vector( edgeprof, vector, b2FaceData, &
       &   gridID, gridSubsetID, gridSubsetInd )
    implicit none
    type(ids_edge_profiles), intent(inout) :: edgeprof
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
#if IMAS_MINOR_VERSION < 15
      idsdata => b2_IMAS_Transform_Data_B2_To_IDS(    &
               &   edgeprof%ggd( slice_index )%grid,    &
               &   GRID_SUBSET_Y_ALIGNED_EDGES, IDSmap, b2FaceData)
#else
      idsdata => b2_IMAS_Transform_Data_B2_To_IDS(    &
               &   edgeprof%grid_ggd( slice_index ),    &
               &   GRID_SUBSET_Y_ALIGNED_EDGES, IDSmap, b2FaceData)
#endif
#if GGD_MINOR_VERSION > 8
      call gridWriteData( vector, gridId, GRID_SUBSET_Y_ALIGNED_EDGES, idsdata )
#else
      call gridWriteData( vector, GRID_SUBSET_Y_ALIGNED_EDGES, idsdata )
#endif
      deallocate(idsdata)
#if IMAS_MINOR_VERSION < 15
      idsdata => b2_IMAS_Transform_Data_B2_To_IDS(    &
               &   edgeprof%ggd( slice_index )%grid,    &
               &   GRID_SUBSET_X_ALIGNED_EDGES, IDSmap, b2FaceData)
#else
      idsdata => b2_IMAS_Transform_Data_B2_To_IDS(    &
               &   edgeprof%grid_ggd( slice_index ),    &
               &   GRID_SUBSET_X_ALIGNED_EDGES, IDSmap, b2FaceData)
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
               &   edgeprof%ggd( slice_index )%grid,  &
               &   gridSubsetInd, IDSmap, b2FaceData)
#else
      idsdata => b2_IMAS_Transform_Data_B2_To_IDS(    &
               &   edgeprof%grid_ggd( slice_index ),  &
               &   gridSubsetInd, IDSmap, b2FaceData)
#endif
#if GGD_MINOR_VERSION > 8
      call gridWriteData( vector, gridId, gridSubsetID, idsdata )
#else
      call gridWriteData( vector, gridSubsetID, idsdata )
#endif
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

    return
    end subroutine B2grid_Write_Data_Vector_Components

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
                    ! stop "compute_Coordinate_Unit_Vectors: not able to find&
                    ! &   poloidal neighbour for cell"
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
                    ! stop "compute_Coordinate_Unit_Vectors: not able to find&
                        ! &    toroidal neighbour for cell"
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

    end subroutine compute_Coordinate_Unit_Vectors

    !> Return unit vector along direction of given vector
    function unitVector(v) result(unitV)
    implicit none
    real(IDS_real), intent(in) :: v(:)  !< Vector
    real(IDS_real) :: unitV(size(v))    !< Unit vector

        unitV = v / sqrt( sum( v**2 ) )
    end function unitVector

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

#else
# ifdef ITM_ENVIRONMENT_LOADED

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


    !! allocate and init the cpo
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
       call species(is, edgecpo%species(is+1)%label, .false.)
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
        call write_cell_scalar( edgecpo%fluid%ne%source, sne(:,:,0) + sne(:,:,1)*ne )

        !! na
        allocate(edgecpo%fluid%ni(ns))
        do is = 1, ns
            call write_quantity( edgecpo%fluid%ni(is)%value, edgecpo%fluid%ni(is)%flux, na(:,:,is-1), fna(:,:,:,is-1) )
            call write_cell_scalar( edgecpo%fluid%ni(is)%source, sna(:,:,0,is-1) + sna(:,:,1,is-1)*na(:,:,is-1) )
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

            call write_cell_scalar( edgecpo%fluid%vi(is)%comps(1)%value, ua(:,:,is-1) )
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
        call compute_Coordinate_Unit_Vectors(crx, cry, e(:,:,:,1), e(:,:,:,2), e(:,:,:,3))

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
      real(ITM_R8), intent(in) :: value(-1:CPOmap%b2nx, -1:CPOmap%b2ny)
      real(ITM_R8), intent(in) :: flux(-1:CPOmap%b2nx, -1:CPOmap%b2ny, 0:1)
      real(ITM_R8), dimension(:), pointer :: cpodata
      real(ITM_R8) :: weight(-1:CPOmap%b2nx, -1:CPOmap%b2ny, TO_SELF:TO_TOP)
      integer i

      allocate(values(5))
      cpodata => b2ITMTransformDataB2ToCpo( edgecpo%grid, B2_SUBGRID_CELLS, CPOmap, value )
      call gridWriteData( values(1), B2_SUBGRID_CELLS, cpodata )
      deallocate(cpodata)
      tmpFace = 0.0_ITM_R8
      do i = TO_SELF, TO_TOP
        weight(:,:,i)=vol(:,:)
      end do
      call value_on_faces(nx,ny,weight,value,tmpFace)
      cpodata => b2ITMTransformDataB2ToCpo( edgecpo%grid, iSgCore, CPOmap, tmpFace )
      call gridWriteData( values(2), iSgCore, cpodata )
      deallocate(cpodata)
      tmpVx = interpolateToVertices( CPOmap%b2nx, CPOmap%b2ny, VX_LOWERLEFT, value )
      cpodata => b2ITMTransformDataB2ToCpoVertex( edgecpo%grid, iSgInnerMidplane, CPOmap, tmpVx  )
      call gridWriteData( values(3), iSgInnerMidplane, cpodata )
      deallocate(cpodata)
      cpodata => b2ITMTransformDataB2ToCpoVertex( edgecpo%grid, iSgOuterMidplane, CPOmap, tmpVx )
      call gridWriteData( values(4), iSgOuterMidplane, cpodata )
      deallocate(cpodata)
      cpodata => b2ITMTransformDataB2ToCpoVertex( edgecpo%grid, B2_SUBGRID_NODES, CPOmap, tmpVx )
      call gridWriteData( values(5), B2_SUBGRID_NODES, cpodata )
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
      cpodata => b2ITMTransformDataB2ToCpo( edgecpo%grid, B2_SUBGRID_CELLS, CPOmap, b2CellData )
      call gridWriteData( scalar(1), B2_SUBGRID_CELLS, cpodata )
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
         cpodata => b2ITMTransformDataB2ToCpo(edgecpo%grid, B2_SUBGRID_CELLS, CPOmap, vecdata(:,:,i-1))
         call gridWriteData( vector%comp(i), B2_SUBGRID_CELLS, cpodata )
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
!!$          cpodata => b2ITMTransformDataB2ToCpo(edgecpo%grid, B2_SUBGRID_EDGES_Y, CPOmap, b2FaceData)
!!$          call gridWriteData( vector%comp(1), B2_SUBGRID_EDGES_Y, cpodata )
!!$          deallocate(cpodata)
!!$          cpodata => b2ITMTransformDataB2ToCpo(edgecpo%grid, B2_SUBGRID_EDGES_X, CPOmap, b2FaceData)
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
!!$          cpodata => b2ITMTransformDataB2ToCpo(edgecpo%grid, subgridInd, CPOmap, b2FaceData)
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
