!!-----------------------------------------------------------------------------
!! DOCUMENTATION (doxygen 1.8.8):
!>      @author
!>      Dejan Penko
!!
!>      @page b2uw_b2 b2_ual_rewrite
!>      @section b2uw_b2_desc   Description
!!      b2_ual_rewrite code is used to generate b2_ual_rewrite.exe
!!      (main program), which is intended to be used within SOLPS-GUI.
!!      The code reads the plasma grid
!!      geometry (full geometry descriptions of all available grid subsets)
!!      and plasma state (electron density/temperature, ion temperature,
!!      velocity etc.). The code then writes the obtained data to IDS database
!!      with the use of b2mod scripts that utilize IMAS GGD Grid Service
!!      Library routines. It erases older IDS versions and rewrites the
!!      data in its current corrected form.
!!
!!      @note   More on b2_ual_rewrite is available in SOLPS-GUI
!!              documentation \b HOWTOs under section <b> 4.5 IMAS </b>.
!!      @note   More information on this b2_ual_rewrite is available in
!!              SOLPS-GUI documentation \b HOWTOs under section
!!              <b> 4.6 Put IDS and Get IDS functions </b>.
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
!!
!!-----------------------------------------------------------------------------

program b2_ual_rewrite

    use b2mod_main
    use b2mod_driver
    use b2mod_grid_mapping
    use ids_routines &  ! IGNORE
     & , only : imas_create_env
    use ids_schemas &   ! IGNORE
     & , only : ids_edge_profiles, ids_edge_sources, ids_edge_transport, &
     &          ids_radiation, ids_dataset_description, ids_equilibrium
    use b2mod_ual &
     & , only : new_ids_edge, delete_ids_edge
    use b2mod_ual_io &
     & , only : b25_process_ids
#if IMAS_MINOR_VERSION > 21
    use ids_schemas &   ! IGNORE
     & , only : ids_summary
#endif
#if ( IMAS_MINOR_VERSION > 25 && IMAS_MINOR_VERSION < 34 )
    use ids_schemas &   ! IGNORE
     & , only : ids_numerics
#endif
#if IMAS_MINOR_VERSION > 30
    use ids_schemas &   ! IGNORE
     & , only : ids_divertors
#endif
#ifdef B25_EIRENE
    use eirmod_comusr
    use eirmod_extrab25
#endif
    use b2mod_ipmain
    implicit none
#ifdef USE_PXFGETENV
    integer lenval, ierror
#else
#ifdef NAGFOR
    integer lenval, ierror
#endif
#endif
#ifndef NO_GETENV
    character(len=24) :: device_env
#endif
    logical streql
    external ipgeti, streql

    !! Local variables
    character(len=24) :: shot_string
    character(len=24) :: run_string
    character(len=24) :: new_run_string
    character(len=24) :: argName
    integer narg, cptArg, new_run
    character*16 usrnam
    logical same_run_number
    data new_run / 0 /
    external usrnam

    !! Set default value for IMAS major version and IDS treename
    status = 0
    write(version,'(i1)') IMAS_MAJOR_VERSION
    treename = 'ids'
    same_run_number = .true.
    write (*,*) 'Starting b2mn init'
    call b2mn_init
    ! call b2mn_step(0)
#ifdef B25_EIRENE
    CALL EIRENE_ALLOC_COMUSR(1)
    call eirene_extrab25_eirpbls_init(nmol,nion,npls)
#endif
    ! read plasma state
    call cfopen(56,'b2fplasma','old','unformatted')
    call cfverr(56, b2fplasma_version)
    call read_b2mod_geo(nx, ny, 56)
    call read_b2mod_plasma(nx, ny, ns, 56)
    call read_b2mod_residuals(56)
    call read_b2mod_sources(56)
    call read_b2mod_transport(56)

    call ipgeti('b2mndr_shot_number', shot )
    if (shot.gt.0) then
      write(shot_string,'(i8)') shot
      call strip_spaces(shot_string)
    end if
    call ipgeti('b2mndr_run_number', run )
    if (run.gt.0) then
      write(run_string,'(i5)') run
      call strip_spaces(run_string)
      new_run = run
      write(new_run_string,'(i5)') new_run
      call strip_spaces(new_run_string)
    end if
    username = usrnam()
    call ipgetc('b2mndr_user', username )
    database = 'solps-iter'
#ifndef NO_GETENV
    device_env = ' '
#ifdef NAGFOR
    call get_environment_variable('DEVICE', status=ierror, length=lenval)
    if (ierror.eq.0) call get_environment_variable('DEVICE', value=device_env)
    call get_environment_variable('IMAS_VERSION', status=ierror, length=lenval)
    if (ierror.eq.0) call get_environment_variable('IMAS_VERSION', value=imas_version)
#else
#ifdef USE_PXFGETENV
    CALL PXFGETENV ('DEVICE', 0, device_env, lenval, ierror)
    CALL PXFGETENV ('IMAS_VERSION', 0, imas_version, lenval, ierror)
#else
    call getenv ('DEVICE', device_env)
    call getenv ('IMAS_VERSION', imas_version)
#endif
#endif
    if (.not.streql(device_env,' ')) database = device_env
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
          if (same_run_number) then
            new_run_string = run_string
            new_run = run
          end if
        case("--newrun","-n")
          call get_command_argument( cptArg + 1, new_run_string )
          !! Transform dummy string variable to integer
          read( new_run_string, *) new_run
          same_run_number = new_run.eq.run
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
    call xertst( 0.le.new_run.and.new_run.le.99999, 'Invalid new run number')
    call xertst( new_run.ge.run, 'New run number must be larger than old one!')
    call xertst( .not.streql(username,' '), 'User name not defined !')
    call xertst( .not.streql(database,' '), 'Database not defined !')

    write(*,'(a,i6,a,i5,4a)') 'Shot: ', shot, ' Run: ', run, &
        & ' User: ', trim(username), ' Database: ', trim(database)
    if (.not.same_run_number) then
      write(*,'(a,i5)') ' will be rewritten with run number ',new_run
      if (database.eq.'iter') write(*,'(a)') ' in ITER database'
    else if (database.eq.'iter') then
      write(*,'(a)') ' will be rewritten in ITER database'
    end if

    !! Process B2.5 data and set it to IMAS IDS
    write(*,*) "START B25_process_ids"
    write (0,*) "Checking if IMAS data-entry already exists : ", trim(database), shot, run
    call imas_open_env(treename, shot, run, idx, &
      &                username, database, version, status)
    if ( status.eq.0 .and. idx.ne.0 ) then
      write (0,*) "Reading old IMAS data-entry: ", trim(database), shot, run
      call ids_get( idx, "equilibrium", equilibrium, status)
      if(status.ne.0) write(0,*) 'Error opening equilibrium IDS !'
      call ids_get( idx, "dataset_description", old_description, status)
      old_imas_version = 'x.xx.x'
      if ( status.ne.0 ) then
        write (0,*) 'Error opening old dataset_description IDS !'
      else if (associated(old_description%dd_version)) then
        old_imas_version = old_description%dd_version(1)
      end if
      call ids_deallocate( old_description )
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
! Copy the IDS to a temporary location with the new DD and then bring it back
        tmp_run = run
        if (new_run.eq.run .and. database.ne.'iter') then
          tmp_run = run + 1000
#if IMAS_MINOR_VERSION > 31
          write(systemarg,'(a,i7,a,i4,a,i7,a,i4,a,a,a,a)') &
     &     'idscp --setDatasetVersion'//                   &
     &        ' -si ',shot,' -ri ',run,                    &
     &        ' -so ',shot,' -ro ',tmp_run,                &
     &        ' -d ',trim(database),' -u ',trim(username)
#else
          write(systemarg,'(a,i7,a,i4,a,i7,a,i4,a,a,a,a)') &
     &     'idscp -si ',shot,' -ri ',run,                  &
     &          ' -so ',shot,' -ro ',tmp_run,              &
     &          ' -d ',trim(database),' -u ',trim(username)
#endif
          if (database.eq.'iter') systemarg = trim(systemarg)//' -do ITER'
#ifdef NAGFOR
          call system(systemarg, status, ierror)
#else
          call system(systemarg)
#endif
          if (database.eq.'iter') database = 'ITER'
        end if
#if IMAS_MINOR_VERSION > 31
        write(systemarg,'(a,i7,a,i4,a,i7,a,i4,a,a,a,a)') &
     &   'idscp --setDatasetVersion'//                   &
     &        ' -si ',shot,' -ri ',tmp_run,              &
     &        ' -so ',shot,' -ro ',new_run,              &
     &        ' -d ',trim(database),' -u ',trim(username)
#else
        write(systemarg,'(a,i7,a,i4,a,i7,a,i4,a,a,a,a)') &
     &   'idscp -si ',shot,' -ri ',tmp_run,              &
     &        ' -so ',shot,' -ro ',new_run,              &
     &        ' -d ',trim(database),' -u ',trim(username)
#endif
        if (database.eq.'iter') systemarg = trim(systemarg)//' -do ITER'
#ifdef NAGFOR
        call system(systemarg, status, ierror)
#else
        call system(systemarg)
#endif
        if (database.eq.'iter') database = 'ITER'
        call imas_open_env(treename, shot, new_run, idx, &
          &                username, database, version, status )
        call xertst( status.eq.0, 'Error recreating IDS with new DD version !')
      else if (.not.same_run_number) then
        call close_ual(idx)
        idx = 0
        call imas_open_env(treename, shot, new_run, idx, &
          &                username, database, version, status)
        if ( status.ne.0 .or. idx.eq.0 .or. database.eq.'iter') then ! New run IDS must be created
          idx = 0
          if (database.eq.'iter') database = 'ITER'
          call imas_create_env(treename, shot, new_run, 0, 0, idx, &
            &                  username, database, version, status )
          call xertst( status.eq.0, 'Error creating new IDS file !')
        end if
      end if
    else
      write (0,*) "No previous IMAS data-entry found, a new one will be created"
      idx = 0
      if (database.eq.'iter') database = 'ITER'
      call imas_create_env(treename, shot, run, 0, 0, idx, &
        &                  username, database, version, status )
      call xertst( status.eq.0, 'Error creating new IDS !')
    end if
    !! Create/Write the set data to IDSs
    call B25_process_ids( edge_profiles, edge_sources, edge_transport, &
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
      &  tim, dteff, shot, new_run, database, version, new_eq_ggd )

    write(*,*) "START new_ids_edge"
    call new_ids_edge( edge_profiles, edge_sources, edge_transport, &
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
        &   idx, new_eq_ggd )
    systemarg = 'create_db_entry -u '//trim(username)//' -d '//trim(database) &
        &  //' -s '//trim(shot_string)//' -r '//trim(new_run_string)
    write(*,*) trim(systemarg)
#ifdef NAGFOR
    call system(systemarg, status, ierror)
#else
    call system(systemarg)
#endif
    if (.not.same_run_number) then
! Add superceding information to .yaml file
      systemarg = 'IDS_yaml_replace '//trim(shot_string)//' '// &
        &  trim(run_string)//' '//trim(new_run_string)
      write(*,*) trim(systemarg)
#ifdef NAGFOR
      call system(systemarg, status, ierror)
#else
      call system(systemarg)
#endif
    end if
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

end program b2_ual_rewrite

!!!Local Variables:
!!! mode: f90
!!! End:
