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
    use b2mod_plasma
    use b2mod_constants
    use b2mod_sources
    use b2mod_transport
    use b2mod_anomalous_transport
    use b2mod_indirect
    use b2mod_interp
    use b2mod_b2cmrc
    use b2mod_version
    use b2mod_grid_mapping

    use logging

    !! UAL Access
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
    use ids_routines       ! IGNORE
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
            &   time_IN, time_step_IN, &
            &   time_slice_ind_IN, num_time_slices_IN )
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
        real(IDS_real), intent(in), optional :: time_IN !< Time 
        real(IDS_real), intent(in), optional :: time_step_IN !< Time step
        integer, intent(in), optional :: time_slice_ind_IN
            !< Time step index for the current time slice
        integer, intent(in), optional :: num_time_slices_IN
            !< Total number of time steps. It is required to beforehand allocate
            !< required ggd(:) array of nodes structure and for additional
            !< checks for correct use of the routine.

        !! Internal variables
        character(len=132) :: ion_label !< Ion species label (e.g. D+1)
        character(len=12) :: ion_charge !< Ion charge (e.g. '1', '2', etc.)
        integer :: ion_charge_int !< Ion charge (e.g. 1, 2, etc.)
        integer :: ion_label_tlen !< Length of the (trimmed) ion label
        integer :: ns     !< Total number of ion species
        integer :: nx     !< Specifies the number of interior cells
                          !< along the first coordinate
        integer :: ny     !< Specifies the number of interior cells
                          !< along the second coordinate
        integer :: is     !< Species index (iterator)
        integer :: i      !< Iterator
        integer :: ntimes !< Number of previous timesteps in IDS
        integer :: o      !< Dummy integer
        integer :: p      !< Dummy integer
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
        real(IDS_real) :: time  !< Generic time
        real(IDS_real) :: time_step !< Time step
        real(IDS_real) :: time_slice_value   !< Time slice value
        type(B2GridMap) :: gmap !< Data structure holding an
            !< intermediate grid description to be transferred into a CPO or IDS

        integer, save :: use_eirene = 0
        character*8 date
        character*10 ctime
        character*5 zone
        integer tvalues(8)
        character*16 usrnam
        character*8 imas_version, ual_version
        character*32 B25_git_version
        character*32 get_B25_hash
        logical streql
#ifdef USE_PXFGETENV
        integer lenval, ierror
#else
#ifdef NAGFOR
        integer lenval, ierror
#endif
#endif
        external usrnam, streql, get_B25_hash

        !! ===  SET UP IDS ===
        write(0,*) "Setting data for edge_profiles IDS"
        call ipgeti('b2mndr_eirene', use_eirene)
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

        ns = size( na, 3 )
        nx = ubound( na, 1 )
        ny = ubound( na, 2 )

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
#endif
            allocate( edge_transport%model(1) )
            edge_transport%model(1)%identifier%index = 1
            allocate( edge_sources%source(1) )
            edge_sources%source(1)%identifier%index = 1
            allocate( edge_transport%model(1)%ggd( num_time_slices ) )
            allocate( edge_sources%source(1)%ggd( num_time_slices ) )
            allocate( edge_transport%model(1)%ggd( num_time_slices )% &
                &   electrons%energy%flux(1) )

            !! Allocate and init the IDS
            allocate( edge_profiles%code%name(1) )
            allocate( edge_transport%code%name(1) )
            allocate( edge_sources%code%name(1) )
# ifdef B25_EIRENE
            edge_profiles%code%name = "SOLPS-ITER"
            edge_transport%code%name = "SOLPS-ITER"
            edge_sources%code%name = "SOLPS-ITER"
# else
            edge_profiles%code%name = "B2.5"
            edge_transport%code%name = "B2.5"
            edge_sources%code%name = "B2.5"
# endif
            allocate( edge_profiles%code%version(1) )
            edge_profiles%code%version = newversion
            allocate( edge_transport%code%version(1) )
            edge_transport%code%version = newversion
            allocate( edge_sources%code%version(1) )
            edge_sources%code%version = newversion

            B25_git_version = get_B25_hash()
            allocate( edge_profiles%code%commit(1) )
            edge_profiles%code%commit = B25_git_version
            allocate( edge_transport%code%commit(1) )
            edge_transport%code%commit = B25_git_version
            allocate( edge_sources%code%commit(1) )
            edge_sources%code%commit = B25_git_version

            allocate( edge_profiles%code%repository(1) )
            edge_profiles%code%repository = "git.iter.org"
            allocate( edge_transport%code%repository(1) )
            edge_transport%code%repository = "git.iter.org"
            allocate( edge_sources%code%repository(1) )
            edge_sources%code%repository = "git.iter.org"

#if IMAS_MINOR_VERSION > 14
            allocate( edge_profiles%ids_properties%provider(1) )
            edge_profiles%ids_properties%provider = usrnam()
            allocate( edge_transport%ids_properties%provider(1) )
            edge_transport%ids_properties%provider = usrnam()
            allocate( edge_sources%ids_properties%provider(1) )
            edge_sources%ids_properties%provider = usrnam()

            allocate( edge_profiles%ids_properties%creation_date(1) )
            edge_profiles%ids_properties%creation_date = &
                &   date//' '//ctime//' '//' '//zone
            allocate( edge_transport%ids_properties%creation_date(1) )
            edge_transport%ids_properties%creation_date = &
                &   date//' '//ctime//' '//' '//zone
            allocate( edge_sources%ids_properties%creation_date(1) )
            edge_sources%ids_properties%creation_date = &
                &   date//' '//ctime//' '//' '//zone
#endif
#if IMAS_MINOR_VERSION > 21
            allocate( edge_profiles%ids_properties%version_put%data_dictionary(1) )
            edge_profiles%ids_properties%version_put%data_dictionary = imas_version
            allocate( edge_transport%ids_properties%version_put%data_dictionary(1) )
            edge_transport%ids_properties%version_put%data_dictionary = imas_version
            allocate( edge_sources%ids_properties%version_put%data_dictionary(1) )
            edge_sources%ids_properties%version_put%data_dictionary = imas_version

            allocate( edge_profiles%ids_properties%version_put%access_layer(1) )
            edge_profiles%ids_properties%version_put%access_layer = ual_version
            allocate( edge_transport%ids_properties%version_put%access_layer(1) )
            edge_transport%ids_properties%version_put%access_layer = ual_version
            allocate( edge_sources%ids_properties%version_put%access_layer(1) )
            edge_sources%ids_properties%version_put%access_layer = ual_version

            allocate( edge_profiles%ids_properties%version_put%access_layer_language(1) )
            edge_profiles%ids_properties%version_put%access_layer_language = 'FORTRAN'
            allocate( edge_transport%ids_properties%version_put%access_layer_language(1) )
            edge_transport%ids_properties%version_put%access_layer_language = 'FORTRAN'
            allocate( edge_sources%ids_properties%version_put%access_layer_language(1) )
            edge_sources%ids_properties%version_put%access_layer_language = 'FORTRAN'
#endif
        end if

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
#endif
#endif

        !! Allocate and set time slice value
#if IMAS_MINOR_VERSION > 14
        edge_profiles%grid_ggd( time_sind )%time = time_slice_value
        edge_transport%model(1)%ggd( time_sind )%time = time_slice_value
        edge_sources%grid_ggd( time_sind )%time = time_slice_value
#endif
        edge_profiles%ggd( time_sind )%time = time_slice_value
        edge_transport%model(1)%ggd( time_sind )%time = time_slice_value
        edge_sources%source(1)%ggd( time_sind )%time = time_slice_value

        !! List of species
        allocate( edge_profiles%ggd( time_sind )%ion( ns ) )
        do is = 0, ns-1
            allocate( edge_profiles%ggd( time_sind )%ion( is + 1 )%label(1) )
            allocate( edge_profiles%ggd( time_sind )%ion( is + 1 )%state(1) )
            allocate( edge_profiles%ggd( time_sind )%ion( is + 1 )%state(1)%label(1) )
            allocate( edge_profiles%ggd( time_sind )%ion( is + 1 )%element(1) )

            ! Put label to ion(is + 1).state(1).label
            call species( is, edge_profiles%ggd( time_sind )%ion( is + 1 )% &
                &   state(1)%label, .false.)
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
            edge_profiles%ggd( time_sind )%ion( is + 1 )%element(1)%a = am( is )

            ! Put nuclear charge
            edge_profiles%ggd( time_sind )%ion( is + 1 )%element(1)%z_n = zn( is )

            ! Put number of atoms
#if IMAS_MINOR_VERSION < 15
            edge_profiles%ggd( time_sind )%ion( is + 1 )%element(1)%multiplicity = 1.0_R8
#else
            edge_profiles%ggd( time_sind )%ion( is + 1 )%element(1)%atoms_n = 1
#endif

            ! Put minimum Z of the charge state bundle
            ! (z_min = z_max = 0 for a neutral)
            edge_profiles%ggd( time_sind )%ion( is + 1 )%state(1)%z_min = zamin( is )

            ! Put maximum Z of the charge state bundle
            edge_profiles%ggd( time_sind )%ion( is + 1 )%state(1)%z_max = zamax( is )

        enddo

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
            !! TODO: The fluxes are currently in the edge_transport. They are
            !! supposed to be in the edge_profiles (24.10.2017)
            !! use electrons%particles or electrons%energy node?
            !! Currently using electrons%energy node

            !! edge_transport fundamentals
            !! TODO: create it as a subroutine
            allocate( edge_transport%model(1) )
            allocate( edge_transport%model(1)%ggd( num_time_slices ) )

            !! ne: Electron density
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

            !! ti: (Common) Ion Temperature
            allocate( edge_profiles%ggd( time_sind )%ion(1)%temperature(1) )
            allocate( edge_transport%model(1)%ggd( time_sind )%ion( 1 ) )
            call write_quantity(                                            &
                &   val = edge_profiles%ggd( time_sind )%ion(1)%temperature,&
                &   fluxes = edge_transport%model(1)% ggd( time_sind )%     &
                &   ion(1)%energy%flux,                                     &
                &   value = ti/qe,                                          &
                &   flux = fhi,                                             &
                &   time_sind = time_sind )

            !! po: Electric potential
            call write_cell_scalar( edge_profiles%ggd( time_sind )% &
                &   phi_potential, po )

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

            !! write the magnetic field vector in the b2 coordinate system
            call write_cell_vector_component(                                 &
                &   vectorComponent = edge_profiles%ggd( time_sind )%e_field, &
                &   b2CellData = bb(:,:,0:2),                                 &
                &   vectorID = "diamagnetic" )
#else
            call logmsg( LOGINFO, &
            &   "b2mod_ual_io.B25_process_ids: GGD not available, no plasma state writing" )
#endif
        end if

        allocate( edge_profiles%code%output_flag( time_sind ) )
        edge_profiles%code%output_flag( time_sind ) = 0
        allocate( edge_transport%code%output_flag( time_sind ) )
        edge_transport%code%output_flag( time_sind ) = 0
        allocate( edge_sources%code%output_flag( time_sind ) )
        edge_sources%code%output_flag( time_sind ) = 0

        call logmsg( LOGDEBUG, "b2mod_ual_io.B25_process_ids: done" )

        contains

#if IMAS_MINOR_VERSION > 11
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
            integer :: ggdID     !< Grid identifier index
            integer :: ival

#if IMAS_MINOR_VERSION < 15
            ggdId = edge_profiles%ggd(time_sind)%grid%identifier%index
#else
            ggdId = edge_profiles%grid_ggd(time_sind)%identifier%index
#endif
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
#if GGD_MINOR_VERSION > 8
            call gridWriteData( val( ival ), ggdID, GRID_SUBSET_CELLS, idsdata )
#else
            call gridWriteData( val( ival ), GRID_SUBSET_CELLS, idsdata )
#endif
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
#if GGD_MINOR_VERSION > 8
            call gridWriteData( val( ival ), ggdID, iGsCoreBoundary, idsdata )
#else
            call gridWriteData( val( ival ), iGsCoreBoundary, idsdata )
#endif
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
#if GGD_MINOR_VERSION > 8
            call gridWriteData( val( ival ), ggdID, iGsInnerMidplane, idsdata )
#else
            call gridWriteData( val( ival ), iGsInnerMidplane, idsdata )
#endif
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
#if GGD_MINOR_VERSION > 8
            call gridWriteData( val( ival ), ggdID, iGsOuterMidplane, idsdata )
#else
            call gridWriteData( val( ival ), iGsOuterMidplane, idsdata )
#endif
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
#if GGD_MINOR_VERSION > 8
            call gridWriteData( val( ival ), ggdID, GRID_SUBSET_NODES, idsdata )
#else
            call gridWriteData( val( ival ), GRID_SUBSET_NODES, idsdata )
#endif
            deallocate( idsdata )

            !! Write data on fluxes for y-aligned faces, x-aligned faces
            !! (those two are set as default) and Core boundary grid subset
            allocate( fluxes(2) )
            call write_face_vector( fluxes(1), flux, time_sind, ggdID)
            call write_face_vector( fluxes(2), flux, time_sind, ggdID, &
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

#if GGD_MINOR_VERSION > 8
            call gridWriteData( val( ival ), ggdID, iGsCore, idsdata )
#else
            call gridWriteData( val( ival ), iGsCore, idsdata )
#endif
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

#if GGD_MINOR_VERSION > 8
            call gridWriteData( val( ival ), ggdID, iGsSOL, idsdata )
#else
            call gridWriteData( val( ival ), iGsSOL, idsdata )
#endif
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
#if GGD_MINOR_VERSION > 8
            call gridWriteData( val( ival ), ggdID, iGsIDivertor, idsdata )
#else
            call gridWriteData( val( ival ), iGsIDivertor, idsdata )
#endif
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

#if GGD_MINOR_VERSION > 8
            call gridWriteData( val( ival ), ggdID, iGsODivertor, idsdata )
#else
            call gridWriteData( val( ival ), iGsODivertor, idsdata )
#endif
            deallocate( idsdata )
        end subroutine write_quantity

        !> Write a scalar B2 cell quantity to ids_generic_grid_scalar
        subroutine write_cell_scalar( scalar, b2CellData )
            type(ids_generic_grid_scalar), intent(inout), pointer :: scalar(:)
                !< Type of IDS data structure, designed for scalar data handling
            real(IDS_real), intent(in) :: b2CellData(-1:gmap%b2nx, -1:gmap%b2ny)
            real(IDS_real), dimension(:), pointer :: idsdata    !< Array for
                !< handing data field values
            integer :: ggdID     !< Grid identifier index

#if IMAS_MINOR_VERSION < 15
            ggdId = edge_profiles%ggd(time_sind)%grid%identifier%index
#else
            ggdId = edge_profiles%grid_ggd(time_sind)%identifier%index
#endif
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
#if GGD_MINOR_VERSION > 8
            call gridWriteData( scalar(1), ggdID, GRID_SUBSET_CELLS, idsdata )
#else
            call gridWriteData( scalar(1), GRID_SUBSET_CELLS, idsdata )
#endif
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
            integer :: ggdID     !< Grid identifier index
            character(len=*), intent(in) :: vectorID    !< Vector ID (e.g.
                                                        !< VEC_ALIGN_RADIAL_ID)

#if IMAS_MINOR_VERSION < 15
            ggdId = edge_profiles%ggd(time_sind)%grid%identifier%index
#else
            ggdId = edge_profiles%grid_ggd(time_sind)%identifier%index
#endif
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
                &   ggdID, GRID_SUBSET_CELLS, vectorID, idsdata )
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
                &   gridID, gridSubsetId )
            use ids_grid_data ! IGNORE
            type(ids_generic_grid_scalar), intent(inout) :: vector
                !< Type of IDS data structure, designed for scalar data handling
                !< (in this case 1D vector)
            real(IDS_real), intent(in) ::   &
                &   b2FaceData(-1:gmap%b2nx, -1:gmap%b2ny, 0:1)
            integer, intent(in) :: gridID                    !< Grid identifier index
            integer, intent(in), optional :: gridSubsetId    !< Base grid subset index
            real(IDS_real), dimension(:), pointer :: idsdata !< Dummy array
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
                    &   gridSubsetId, gmap, b2FaceData)
#else
                idsdata => b2_IMAS_Transform_Data_B2_To_IDS(    &
                    &   edge_profiles%grid_ggd( time_sind ),    &
                    &   gridSubsetId, gmap, b2FaceData)
#endif
#if GGD_MINOR_VERSION > 8
                call gridWriteData( vector, gridId, gridSubsetId, idsdata )
#else
                call gridWriteData( vector, gridSubsetId, idsdata )
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
