!!-----------------------------------------------------------------------------
!! DOCUMENTATION (doxygen 1.8.8):
!>      @author
!>      Dejan Penko
!!
!>      @page b2uw_b2 b2_ual_write
!>      @section b2uw_b2_desc   Description
!!      b2_ual_write code is used to generate b2_ual_write.exe
!!      (main program), which is intended to be used within SOLPS-GUI.
!!      The code reads the plasma grid
!!      geometry (full geometry descriptions of all available grid subsets)
!!      and plasma state (electron density/temperature, ion temperature,
!!      velocity etc.). The code then writes the obtained data to IDS database
!!      with the use of b2mod scripts that utilize IMAS GGD Grid Service
!!      Library routines.
!!
!!      @note   More on the b2_ual_writers is available in SOLPS-GUI
!!              documentation \b HOWTOs under section <b> 4.5 IMAS </b>.
!!      @note   More information on this b2_ual_write is available in SOLPS-GUI
!!              documentation \b HOWTOs under section <b> 4.6 Put IDS and Get
!!              IDS functions </b>.
!!
!!      Variables inherited from b2mod_driver module
!!      character(len=24) :: treename   !< The name of the IMAS IDS database
!!        !< (i.e. "edge_profiles" (mandatory) )
!!      character(len=24) :: username   !< Creator/owner of the IMAS IDS database
!!      character(len=24) :: database   !< IMAS IDS database name
!!        !< (i. e. solps-iter, ITER, aug)
!!      character(len=24) :: version    !< Major version of the IMAS IDS database
!!      integer :: idx    !< The returned identifier to be used in the subsequent
!!        !< data access operation
!!      integer :: shot   !< The shot number of the database being created
!!      integer :: run    !< The run number of the database being created
!!      integer :: status !< Returned status for UAL commands
!!      type(ids_edge_profiles) :: edge_profiles !< IDS designed to store data on
!!        !< edge plasma profiles (includes the scrape-off layer and possibly
!!        !< part of the confined plasma)
!!      type (ids_edge_profiles) :: old_edge_profiles
!!      type (ids_edge_sources) :: edge_sources !< IDS designed to store
!!        !< data on edge plasma sources. Energy terms correspond to the full
!!        !< kinetic energy equation (i.e. the energy flux takes into account
!!        !< the energy transported by the particle flux)
!!      type (ids_edge_transport) :: edge_transport !< IDS designed to store
!!        !< data on edge plasma transport. Energy terms correspond to the
!!        !< full kinetic energy equation (i.e. the energy flux takes into
!!        !< account the energy transported by the particle flux)
!!      type (ids_radiation) :: radiation !< IDS designed to store
!!        !< data on radiation emitted by the plasma species
!!      type (ids_dataset_description) :: description !< IDS designed to store
!!        !< a description of the simulation
!!      type (ids_dataset_description) :: old_description
!!      type (ids_summary) :: summary !< IDS designed to store
!!        !< run summary data
!!      type (ids_numerics) :: numerics !< IDS designed to store
!!        !< run numerics data
!!      type (ids_divertors) :: divertors !< IDS designed to store
!!        !< divertor data
!!      integer num_time_slices, time_slice_index
!!      real(IDS_real) :: old_start_time, old_end_time, ids_end_time
!!      logical continued
!!
!!-----------------------------------------------------------------------------

program b2_ual_write

    use b2mod_geo
    use b2mod_main
    use b2mod_driver
    use b2mod_plasma
    use b2mod_sources
    use b2mod_switches
    use b2mod_transport
    use b2mod_grid_mapping
    use b2us_geo
    use b2us_map
    use b2us_plasma
    use b2mod_ual    &
     & , only : put_ids_edge, b25_process_ids, &
     &          ids_edge_profiles, ids_edge_sources, ids_edge_transport, &
     &          ids_radiation, ids_dataset_description, ids_equilibrium
    use b2mod_ual_io
#if IMAS_MINOR_VERSION > 21
    use b2mod_ual    &
     & , only : ids_summary
#endif
#if IMAS_MINOR_VERSION > 30
    use b2mod_ual    &
     & , only : ids_divertors
#endif
#ifdef B25_EIRENE
    use eirmod_comusr
    use eirmod_extrab25
#endif
    use b2mod_ipmain
    implicit none
#ifndef NO_GETENV
    character(len=24) :: device_env
    integer lenval, ierror
#ifndef USE_PXFGETENV
    intrinsic get_environment_variable
#endif
#endif
    logical streql
    external ipgeti, streql

    !! Local variables
    type (geometry) :: geo
    type (mapping) :: mpg
    type (B2state) :: state
    type (B2StateExt) :: state_ext
    type (switches) :: switch
    character(len=24) :: shot_string
    character(len=24) :: run_string
    character(len=24) :: argName
    integer nx, ny
    integer narg, cptArg
    character*16 usrnam
    external usrnam

    !! Set default value for IMAS major version and IDS treename
    status = 0
    write(version,'(i1)') IMAS_MAJOR_VERSION
    treename = 'ids'
    write (*,*) 'Starting b2mn init'
    call b2mn_init (switch, geo, mpg, state, state_ext)
    ! call b2mn_step(0)
#ifdef B25_EIRENE
    CALL EIRENE_ALLOC_COMUSR(1)
    call eirene_extrab25_eirpbls_init(1,natm,nmol,nion,npls,nstrat,npls)
#endif
    ! read plasma state
    call cfopen(56,'b2fplasma','old','unformatted')
    call cfverr(56, b2fplasma_version)
    call read_b2mod_geo(nx, ny, 56)
    call read_b2mod_plasma(nx, ny, ns, 56)
    call read_b2mod_residuals(56)
    call read_b2mod_sources(56)
    call read_b2mod_transport(nx, ny, ns, 56)

    call ipgeti('b2mndr_shot_number', shot )
    call ipgeti('b2mndr_run_number', run )
    username = usrnam()
    call ipgetc('b2mndr_user', username )
    database = 'solps-iter'
#ifdef NO_GETENV
    device_env = ' '
#else
#ifdef USE_PXFGETENV
    CALL PXFGETENV ('DEVICE', 0, device_env, lenval, ierror)
    CALL PXFGETENV ('IMAS_VERSION', 0, imas_version, lenval, ierror)
#else
    call get_environment_variable ('DEVICE', status=ierror, length=lenval)
    if (ierror.eq.0) call get_environment_variable('DEVICE', value=device_env)
    call get_environment_variable ('IMAS_VERSION', status=ierror, length=lenval)
    if (ierror.eq.0) call get_environment_variable('IMAS_VERSION', value=imas_version)
#endif
    if (.not.streql(device_env,' ')) database = device_env
    if (streql(database,'iter')) database = 'ITER'
#endif
    call ipgetc('b2mndr_device', database )
    call ipgetc('b2mndr_database', database )
    ! Check for optional command line arguments
    ! which will supersede input from b2mn.dat if present
    narg = command_argument_count()
    do cptArg = 1, narg
      call get_command_argument( cptArg, argName )
      select case( adjustl( argName ) )
        case("--shot","-s")
          call get_command_argument( cptArg + 1, shot_string )
          !! Transform dummy string variable to integer
          read( shot_string, *) shot
        case("--run","-r")
          call get_command_argument( cptArg + 1, run_string )
          !! Transform dummy string variable to integer
          read( run_string, *) run
        case("--username","-u")
          call get_command_argument( cptArg + 1, username )
        case("--database","--device","-d")
          call get_command_argument( cptArg + 1, database )
        case("--version","-v")
          call get_command_argument( cptArg + 1, version )
      end select
    end do

    call xertst( 0.lt.shot.and.shot.le.214748, 'Invalid shot number')
    call xertst( 0.le.run.and.run.le.99999, 'Invalid run number')
    call xertst( .not.streql(username,' '), 'User name not defined !')
    call xertst( .not.streql(database,' '), 'Database not defined !')

    write(*,'(a,i8,a,i8,4a)') 'Shot: ', shot, ' Run: ', run, &
        & ' User: ', trim(username), ' Database: ', trim(database)

    !! Process B2.5 data and set it to IMAS IDS
    write(*,*) "START B25_process_ids"
    write (0,*) "Checking if IMAS data-entry already exists : ", trim(database), shot, run
    call imas_open_env(treename, shot, run, idx, username, database, version, status)
    if (status.ne.0) then
      if (database.eq.'ITER') then
        write(*,*) "Did not find ITER database IDS file."
        write(*,*) "Checking if old ''iter'' case exists."
        call imas_open_env(treename, shot, run, idx, username, 'iter', version, status)
        if (status.eq.0) then
          database = 'iter'
          write(*,*) "Old database case found."
          write(*,*) "Will be rewritten in new location."
        end if
      else if (database.eq.'iter') then
        call imas_open_env(treename, shot, run, idx, username, 'ITER', version, status)
        database = 'ITER'
      end if
    end if
    !! If this is a time continuation run, append the new data to the IDS
    if ( status.eq.0 .and. idx.ne.0 ) then
      write (0,*) "Reading old IMAS data-entry: ", trim(database), shot, run
      call ids_get( idx, "equilibrium", equilibrium, status)
      call ids_get( idx, "edge_profiles", old_edge_profiles, status)
      if ( status.ne.0 ) then
        write (0,*) 'Error opening old edge_profiles IDS ! Will create a new one.'
        idx = 0
        continued = .false.
      else
        num_time_slices = size(old_edge_profiles%time)
        if (num_time_slices.gt.0) then
          ids_end_time = old_edge_profiles%time(num_time_slices)
        else
          ids_end_time = IDS_REAL_INVALID
        end if
        call ids_deallocate( old_edge_profiles )
        old_start_time = 0.0_IDS_real
        old_end_time = IDS_REAL_INVALID
        old_imas_version = 'x.xx.x'
        call ids_get( idx, "dataset_description", old_description, status)
        if ( status.ne.0 ) then
          write (0,*) 'Error opening old dataset_description IDS !'
        else if (associated(old_description%dd_version)) then
#if ( IMAS_MINOR_VERSION > 25 && IMAS_MINOR_VERSION < 34 )
          old_start_time = old_description%simulation%time_begin
          old_end_time = old_description%simulation%time_end
#endif
          old_imas_version = old_description%dd_version(1)
          call ids_deallocate( old_description )
        end if
        continued = run_start_time.eq.IDS_REAL_INVALID .and. &
           &       (ids_end_time.lt.tim .and. ids_end_time.ne.IDS_REAL_INVALID)
        continued = continued .or. &
           &       (run_start_time.ge.ids_end_time .and. &
           &       (ids_end_time.lt.tim .and. ids_end_time.ne.IDS_REAL_INVALID))
        if (continued.or.database.eq.'iter') then
          if (.not.streql(old_imas_version,imas_version).or.database.eq.'iter') then
            if (.not.streql(old_imas_version,imas_version)) then
              write(*,*) &
               & 'Old IMAS data-entry was written using Data Dictionary version '// &
               &  trim(old_imas_version)//'.'
              write(*,*) &
               & 'Recreating using Data Dictionary version '// &
               &  trim(imas_version)//'.'
            end if
            if (database.eq.'iter') &
              &  write(*,*) 'IDS file will be moved to ITER database.'
            call close_ual(idx)
            idx = 0
!xpb Copy the IDS to a temporary location with the new DD and then bring it back
            tmp_run = run
            if (database.ne.'iter') tmp_run = run + 1000
#if IMAS_MINOR_VERSION > 31
            write(systemarg,'(a,i7,a,i4,a,i7,a,i4,a,a,a,a)') &
               & 'idscp --setDatasetVersion'//                 &
               &       ' -si ',shot,' -ri ',run,               &
               &       ' -so ',shot,' -ro ',tmp_run,           &
               &       ' -d ',trim(database),' -u ',trim(username)
#else
            write(systemarg,'(a,i7,a,i4,a,i7,a,i4,a,a,a,a)') &
               & 'idscp -si ',shot,' -ri ',run,                &
               &      ' -so ',shot,' -ro ',tmp_run,            &
               &      ' -d ',trim(database),' -u ',trim(username)
#endif
            if (database.eq.'iter') systemarg = trim(systemarg)//' -do ITER'
#ifdef NAGFOR
            call system(systemarg, status, ierror)
#else
            call system(systemarg)
#endif
            if (database.ne.'iter') then
#if IMAS_MINOR_VERSION > 31
              write(systemarg,'(a,i7,a,i4,a,i7,a,i4,a,a,a,a)') &
               & 'idscp --setDatasetVersion'//                 &
               &       ' -si ',shot,' -ri ',tmp_run,           &
               &       ' -so ',shot,' -ro ',run,               &
               &       ' -d ',trim(database),' -u ',trim(username)
#else
              write(systemarg,'(a,i7,a,i4,a,i7,a,i4,a,a,a,a)') &
               & 'idscp -si ',shot,' -ri ',tmp_run,            &
               &      ' -so ',shot,' -ro ',run,                &
               &      ' -d ',trim(database),' -u ',trim(username)
#endif
#ifdef NAGFOR
              call system(systemarg, status, ierror)
#else
              call system(systemarg)
#endif
            end if
            if (database.eq.'iter') database = 'ITER'
            call imas_open_env(treename, shot, run, idx, &
             &                 username, database, version, status)
          end if
          if (continued) then
            write (0,*) "Appending a new time slice at t = ", tim, " s."
            num_time_slices = num_time_slices + 1
          end if
          time_slice_index = num_time_slices
          call B25_process_ids( geo, mpg, state, switch, &
             &  edge_profiles, edge_sources, edge_transport, &
             &  radiation, description, equilibrium, &
#if IMAS_MINOR_VERSION > 21
             &  summary, &
#endif
#if ( IMAS_MINOR_VERSION > 25 && IMAS_MINOR_VERSION < 34 )
             &  numerics, old_start_time, run_end_time, &
#endif
#if IMAS_MINOR_VERSION > 30
             &  divertors, &
#endif
             &  tim, dtim, shot, run, database, version, new_eq_ggd, &
             &  time_slice_index, num_time_slices )
        else
          write (0,*) "Not a time continuation, IDS will be overwritten !"
          idx = 0
        end if
      end if
    else
      write (0,*) "No previous IMAS data-entry found, a new one will be created"
      idx = 0
      if (database.eq.'iter') database = 'ITER'
    end if
    if ( status.ne.0 .or. idx.eq.0 ) then
      call B25_process_ids( geo, mpg, state, switch, &
         &  edge_profiles, edge_sources, edge_transport, &
         &  radiation, description, equilibrium, &
#if IMAS_MINOR_VERSION > 21
         &  summary, &
#endif
#if ( IMAS_MINOR_VERSION > 25 && IMAS_MINOR_VERSION < 34 )
         &  numerics, run_start_time, run_end_time, &
#endif
#if IMAS_MINOR_VERSION > 30
         &  divertors, &
#endif
         &  tim, dtim, shot, run, database, version, new_eq_ggd )
    end if

    !! Create/Write the set data to IDSs
    write(*,*) "START put_ids_edge"
    call put_ids_edge( edge_profiles, edge_sources, edge_transport, &
        &   radiation, description, equilibrium, &
#if IMAS_MINOR_VERSION > 21
        &   summary, &
#endif
#if ( IMAS_MINOR_VERSION > 25 && IMAS_MINOR_VERSION < 34 )
        &   numerics, &
#endif
#if IMAS_MINOR_VERSION > 30
        &   divertors, &
#endif
        &   treename, shot, run, idx, username, database, version, &
        &   new_eq_ggd )
    call dealloc_ids_edge( edge_profiles, edge_sources, edge_transport, &
#if ( IMAS_MINOR_VERSION > 25 && IMAS_MINOR_VERSION < 34 )
        &   numerics, &
#endif
#if IMAS_MINOR_VERSION > 30
        &   divertors, &
#endif
        &   radiation )
    call dealloc_batch_edge( batch_profiles, batch_sources, &
#if IMAS_MINOR_VERSION > 21
        &   summary, &
#endif
        &   description )
    call close_ual(idx)
    idx = 0

end program b2_ual_write

!!!Local Variables:
!!! mode: f90
!!! End:
