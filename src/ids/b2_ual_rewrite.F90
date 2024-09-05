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
!!      character(len=256) :: path      !< Directory path of the IMAS data entry
!!      integer :: shot   !< The pulse (previously shot) number of the database
!!        !< being created
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

    use b2mod_main &
     & , only : ns, &
     &          b2mn_init
    use b2mod_driver &
     & , only : idx, ids_path, dtim, &
     &          shot, run, username, database, version, &
     &          old_description, equilibrium, &
     &          old_imas_version, imas_version, new_eq_ggd, &
     &          description, &
     &          edge_profiles, edge_sources, edge_transport, radiation, &
     &          batch_profiles, batch_sources
    use b2mod_time
    use b2mod_version
    use b2mod_switches
    use b2mod_grid_mapping
    use b2us_io
    use b2us_geo
    use b2us_map
    use b2us_plasma
    use b2mod_solpstop
    use b2mod_numerics_namelist
    use ids_routines &  ! IGNORE
     & , only : imas_create_env
    use ids_schemas &   ! IGNORE
     & , only : ids_edge_profiles, ids_edge_sources, ids_edge_transport, &
     &          ids_radiation, ids_dataset_description, ids_equilibrium
    use b2mod_ual &
     & , only : new_ids_edge, delete_ids_edge, &
     &          dealloc_ids_edge, dealloc_batch_edge, close_ual
    use b2mod_ual_io &
     & , only : b25_process_ids
#if ( IMAS_MINOR_VERSION > 21 || IMAS_MAJOR_VERSION > 3 )
    use ids_schemas &   ! IGNORE
     & , only : ids_summary
    use b2mod_driver &
     & , only : summary
#endif
#if ( IMAS_MINOR_VERSION > 25 && IMAS_MINOR_VERSION < 34 && IMAS_MAJOR_VERSION == 3 )
    use ids_schemas &   ! IGNORE
     & , only : ids_numerics
    use b2mod_driver &
     & , only : numerics
#endif
#if ( IMAS_MINOR_VERSION > 30 || IMAS_MAJOR_VERSION > 3 )
    use ids_schemas &   ! IGNORE
     & , only : ids_divertors
    use b2mod_driver &
     & , only : divertors
#endif
    use ids_routines &  ! IGNORE
     & , only : ids_get, ids_deallocate
#if AL_MAJOR_VERSION > 4
    use ids_routines &  ! IGNORE
     & , only : imas_open, al_build_uri_from_legacy_parameters, &
     &          OPEN_PULSE, FORCE_CREATE_PULSE, STRMAXLEN, MDSPLUS_BACKEND
#else
    use ids_routines &  ! IGNORE
     & , only : imas_open_env
    use b2mod_driver &
     & , only : treename
#endif
#ifdef B25_EIRENE
    use eirmod_parmmod
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
    intrinsic index
    external ipgeti, streql

    !! Local variables
    type (geometry) :: geo
    type (mapping) :: mpg
    type (B2state) :: state
    type (B2StateExt) :: state_ext
    type (B2Average) :: state_avg

    type (switches) :: switch
    character(len=24) :: shot_string
    character(len=24) :: run_string
    character(len=24) :: new_run_string
    character(len=24) :: argName
    integer narg, cptArg, l, m, new_run, status, tmp_run
    integer idum(0:2)
    character*16 usrnam, run_user
    character*256 path, new_path, imasdir, home_dir, systemarg
    logical same_run_number, absolute_path, not_default
#if AL_MAJOR_VERSION > 4
    character(len=STRMAXLEN) :: uri, uri_source, uri_dest
    character(len=:), allocatable :: message
#endif
    data new_run / 0 /
    integer len_of_digits
    external usrnam, len_of_digits

    !! Set default value for IMAS major version and IDS treename
    status = 0
    run_user = usrnam()
    call ipgetc('b2mndr_user', run_user )
    call xertst( .not.streql(run_user,' '), 'User name not defined !')
    not_default = .not.streql(run_user, usrnam())
    username = run_user
    home_dir = '/home/'//trim(run_user)
    database = 'solps-iter'
    write(version,'(i1)') IMAS_MAJOR_VERSION
#ifndef NO_GETENV
    device_env = ' '
#ifdef USE_PXFGETENV
    CALL PXFGETENV ('HOME', 0, home_dir, lenval, ierror)
    CALL PXFGETENV ('DEVICE', 0, device_env, lenval, ierror)
#else
    call get_environment_variable('HOME', status=ierror, length=lenval)
    if (ierror.eq.0) call get_environment_variable('HOME', value=home_dir)
    call get_environment_variable('DEVICE', status=ierror, length=lenval)
    if (ierror.eq.0) call get_environment_variable('DEVICE', value=device_env)
#endif
    if (.not.streql(device_env,' ')) database = device_env
    if (streql(database,'iter')) database = 'ITER'
#endif
    call ipgetc('b2mndr_device', database )
    call ipgetc('b2mndr_database', database )
    call xertst( .not.streql(database,' '), 'Database not defined !')
    not_default = not_default.or. &
       &         (.not.streql(database, device_env).and. &
       &          .not.streql(database,'solps-iter'))
#if AL_MAJOR_VERSION > 4
    imasdir = trim(home_dir)//'/public/imasdb/'//trim(database)//'/'//trim(version)
#else
    imasdir = trim(home_dir)//'/public/imasdb/'//trim(database)//'/'//trim(version)//'/0'
#endif
#ifndef NO_GETENV
#ifdef USE_PXFGETENV
    CALL PXFGETENV ('IMASDIR', 0, imasdir, lenval, ierror)
#else
    call get_environment_variable('IMASDIR', status=ierror, length=lenval)
    if (ierror.eq.0) call get_environment_variable('IMASDIR',value=imasdir)
#endif
#endif
    ids_path = ' '
    same_run_number = .true.
    write (*,*) 'Starting b2mn init'
    call b2mn_init (switch, geo, mpg, state, state_ext, state_avg)
    ! call b2mn_step(0)
#ifdef B25_EIRENE
    CALL EIRENE_ALLOC_COMUSR(1)
    call eirene_extrab25_eirpbls_init(nmol,nion,npls)
#endif
    ! read plasma state
    call cfopen(56,'b2fplasma','old','unformatted')
    call cfverr(56, b2fplasma_version)
    ! obtain parameters from b2fplasma file
    call cfruin (56,3,idum,'nCv,nFc,ns')
    call xertst (idum(0).eq.mpg%nCv.and.idum(1).eq.mpg%nFc.and. &
     &           idum(2).eq.state%pl%ns, &
     &          'faulty input nCv, nFc, ns from b2fplasma file')
    call read_b2fplasma(56, mpg%nCv, mpg%nFc, ns, state)

    call ipgeti('b2mndr_pulse_number', shot )
    if (shot.eq.0) call ipgeti('b2mndr_shot_number', shot )
    if (shot.gt.0) then
      write(shot_string,'(i8)') shot
      call strip_spaces(shot_string)
    end if
    call ipgeti('b2mndr_run_number', run )
    if (run.ge.0) then
      write(run_string,'(i5)') run
      call strip_spaces(run_string)
      new_run = run
      write(new_run_string,'(i5)') new_run
      call strip_spaces(new_run_string)
    end if
    path = ' '
    absolute_path = .false.
#if AL_MAJOR_VERSION > 4
    call ipgetc('b2mndr_ids_path', path )
    not_default = not_default.or..not.streql(path,' ')
#endif
    call strip_spaces(path)
    if (.not.streql(path,' ')) absolute_path = path(1:1).eq.'/'
#ifdef NO_GETENV
    write(imas_version,'(i1,a1,i2,a1,i1)') IMAS_MAJOR_VERSION,'.', &
                                         & IMAS_MINOR_VERSION,'.', &
                                         & IMAS_MICRO_VERSION
#else
#ifdef USE_PXFGETENV
    CALL PXFGETENV ('IMAS_VERSION', 0, imas_version, lenval, ierror)
#else
    call get_environment_variable('IMAS_VERSION', status=ierror, length=lenval)
    if (ierror.eq.0) call get_environment_variable('IMAS_VERSION', value=imas_version)
#endif
#endif
    ! Check for optional command line arguments
    ! which will supersede input from b2mn.dat if present
    narg = command_argument_count()
    do cptArg = 1, narg
      call get_command_argument( cptArg, argName )
      select case( adjustl( argName ) )
#if AL_MAJOR_VERSION > 4
        case("--path","-p")
          call get_command_argument( cptArg + 1, path )
          !! Parse HOME and SOLPSTOP and remove IMASDIR prefix if present
          l=index(path,'$HOME')
          if (l.gt.0) then
            ids_path = trim(home_dir)//trim(path(l+5:256))
            path = ids_path
          end if
          l=index(path,'$SOLPSTOP')
          if (l.gt.0) then
            ids_path = trim(solpstop())//trim(path(l+9:256))
            path = ids_path
          end if
          l=index(path,'$HOMEDIR')
          if (l.gt.0) then
            ids_path = trim(imasdir)//trim(path(l+8:256))
            path = ids_path
          end if
          l=index(path,trim(imasdir))
          if (l.eq.0) then
            m=index(path,'/')
          else
            m=index(path(l+len_trim(imasdir):256),'/')
          end if
          absolute_path = l.eq.1.or.(l.eq.0.and.m.eq.1)
          if (absolute_path) then
            ids_path = path
          else if (l.eq.0) then
            ids_path = trim(imasdir)//'/'//trim(path)
          end if
#endif
        case("--pulse","--shot","-s")
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

#if AL_MAJOR_VERSION > 4
    call xertst( 0.le.shot, 'Invalid pulse number')
    call xertst( 0.le.run, 'Invalid run number')
    call xertst( 0.le.new_run, 'Invalid new run number')
#else
    call xertst( 0.lt.shot.and.shot.le.214748, 'Invalid shot number')
    call xertst( 0.le.run.and.run.le.99999, 'Invalid run number')
    call xertst( 0.le.new_run.and.new_run.le.99999, 'Invalid new run number')
#endif
    call xertst( .not.streql(username,' '), 'User name not defined !')
    call xertst( .not.streql(database,' '), 'Database not defined !')
    call xertst( new_run.ge.run, 'New run number must be larger than old one!')
    if (index(imasdir,trim(username)).eq.0) then
      l=index(imasdir,trim(run_user))
      m=index(imasdir(l+len_trim(run_user):256),'/')
      write(imasdir,'(a)') imasdir(1:l)//trim(username)//trim(imasdir(m+l:256))
    end if
    if (index(imasdir,'imasdb/'//trim(database)).eq.0) then
      l=index(imasdir,'imasdb/')
      m=index(imasdir(l+7:256),'/')
      write(imasdir,'(a)') imasdir(1:l+6)//trim(database)//trim(imasdir(m+l+6:256))
    end if
    if (.not.streql(version,int2str(IMAS_MAJOR_VERSION))) then
      l=len_trim(version)
      m=len_trim(imasdir)
#if AL_MAJOR_VERSION > 4
      if (.not.streql(imasdir(m:m),version)) &
        & write(imasdir,'(a)') imasdir(1:m-1)//trim(version)
#else
      if (.not.streql(imasdir(m-2:m-2),version)) &
        & write(imasdir,'(a)') imasdir(1:m-3)//trim(version)//'/'//int2str(run/10000)
    else if (run.ge.10000) then
      m=len_trim(imasdir)
      write(imasdir,'(a)') imasdir(1:m-1)//int2str(run/10000)
#endif
    end if
    if (streql(path,' ')) then
#if AL_MAJOR_VERSION > 4
      call xertst( 0.lt.shot, 'Invalid pulse number !')
      ids_path = trim(imasdir)//'/'//int2str(shot)//'/'//int2str(run)
      new_path = trim(imasdir)//'/'//int2str(shot)//'/'//int2str(new_run)
#else
      ids_path = imasdir
      new_path = imasdir(1:len_trim(imasdir)-1)//int2str(new_run/10000)
#endif
    else if (.not.absolute_path) then
      ids_path = trim(imasdir)//'/'//trim(path)
      if (new_run.gt.0.and.index(path,int2str(run)).gt.0) then
        l=len_of_digits(run)
        new_path = trim(imasdir)//'/'//path(1:len_trim(path)-l)//int2str(new_run)
      else
        new_path = ids_path
      end if
    else
      ids_path = path
      if (new_run.gt.0.and.index(path,int2str(run)).gt.0) then
        l=len_of_digits(run)
        new_path = path(1:len_trim(path)-l)//int2str(new_run)
      else
        new_path = ids_path
      end if
    end if

#if AL_MAJOR_VERSION > 4
    write(0,'(6a)') 'Path: ', trim(ids_path), &
        & ' User: ', trim(username), ' Database: ', trim(database)
#else
    write(0,'(a,i6,a,i5,4a)') 'Shot: ', shot, ' Run: ', run, &
        & ' User: ', trim(username), ' Database: ', trim(database)
#endif
    if (.not.same_run_number.and.new_run.gt.0) then
      write(0,'(a,i5)') ' will be rewritten with run number ',new_run
      if (database.eq.'iter') write(*,'(a)') ' in ITER database'
    else if (database.eq.'iter') then
      write(0,'(a)') ' will be rewritten in ITER database'
    end if

    !! Process B2.5 data and set it to IMAS IDS
    write(*,*) "START B25_process_ids"
#if AL_MAJOR_VERSION > 4
    write(0,'(2a)') "Checking if IMAS data entry already exists : ", &
      &  trim(ids_path)
    uri = 'imas:mdsplus?path='//trim(ids_path)
    call imas_open( uri, OPEN_PULSE, idx, status, message )
#else
    write(0,'(2a,2i8)') &
      & "Checking if IMAS data entry already exists : ", &
      &  trim(database), shot, run
    call imas_open_env(treename, shot, run, idx, &
      &                username, database, version, status)
#endif
    if ( status.eq.0 .and. idx.ne.0 ) then
      write(0,*) "Reading old IMAS data entry: ", trim(database), shot, run
      call ids_get( idx, "equilibrium", equilibrium, status )
      if(status.ne.0) write(0,*) 'Error opening equilibrium IDS !'
      call ids_get( idx, "dataset_description", old_description, status )
      old_imas_version = 'x.xx.x'
      if ( status.ne.0 ) then
        write(0,*) 'Error opening old dataset_description IDS !'
      else if (associated(old_description%ids_properties% &
                      &   version_put%data_dictionary)) then
        old_imas_version = old_description%ids_properties% &
                      &   version_put%data_dictionary(1)
      end if
      call ids_deallocate( old_description )
      if (.not.streql(old_imas_version,imas_version).or.database.eq.'iter') then
        if (.not.streql(old_imas_version,imas_version)) then
          write(0,*) &
            & 'Old IMAS data entry was written using Data Dictionary version '// &
            &  trim(old_imas_version)//'.'
          write(0,*) &
            & 'Recreating using Data Dictionary version '// &
            &  trim(imas_version)//'.'
        end if
        if (database.eq.'iter'.or.index(ids_path,'imasdb/iter').gt.0) &
          &  write(0,*) 'IDS file will be moved to ITER database.'
        call close_ual(idx)
! Copy the IDS to a temporary location with the new DD and then bring it back
        tmp_run = run
        if (new_run.eq.run .and. database.ne.'iter' .and. &
     &      index(new_path,'imasdb/iter').eq.0) then
          tmp_run = run + 1000
#if AL_MAJOR_VERSION > 4
          call al_build_uri_from_legacy_parameters              &
     &     ( MDSPLUS_BACKEND, shot, run, trim(username),        &
     &       trim(database), int2str(IMAS_MAJOR_VERSION), '',   &
     &       uri_source, status )
          if (database.eq.'iter'.or.index(ids_path,'imasdb/iter').gt.0) then
            call al_build_uri_from_legacy_parameters            &
     &       ( MDSPLUS_BACKEND, shot, tmp_run, trim(username),  &
     &         'ITER', int2str(IMAS_MAJOR_VERSION), '',         &
     &         uri_dest, status )
          else
            call al_build_uri_from_legacy_parameters            &
     &       ( MDSPLUS_BACKEND, shot, tmp_run, trim(username),  &
     &         trim(database), int2str(IMAS_MAJOR_VERSION), '', &
     &         uri_dest, status )
          end if
          write(systemarg,'(a,a,a,a)')                     &
     &      'idscp --setDatasetVersion'//                  &
     &         ' -s ',trim(uri_source),                    &
     &         ' -d ',trim(uri_dest)
#else
#if ( IMAS_MINOR_VERSION > 31 || IMAS_MAJOR_VERSION > 3 )
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
          if (database.eq.'iter'.or.index(new_path,'imasdb/iter').gt.0) &
     &     systemarg = trim(systemarg)//' -do ITER'
#endif
#ifdef NAGFOR
          call system(systemarg, status, ierror)
#else
          call system(systemarg)
#endif
          if (database.eq.'iter') database = 'ITER'
          if (index(new_path,'imasdb/iter').gt.0) then
            l=index(new_path,'imasdb/iter')
            write(new_path(l+7:l+10),'(a4)') 'ITER'
          end if
        end if
#if AL_MAJOR_VERSION > 4
        call al_build_uri_from_legacy_parameters              &
     &   ( MDSPLUS_BACKEND, shot, tmp_run, trim(username),    &
     &     trim(database), int2str(IMAS_MAJOR_VERSION), '',   &
     &     uri_source, status )
        if (database.eq.'iter'.or.index(new_path,'imasdb/iter').gt.0) then
          call al_build_uri_from_legacy_parameters            &
     &     ( MDSPLUS_BACKEND, shot, new_run, trim(username),  &
     &       'ITER', int2str(IMAS_MAJOR_VERSION), '',         &
     &       uri_dest, status )
        else
          call al_build_uri_from_legacy_parameters            &
     &     ( MDSPLUS_BACKEND, shot, new_run, trim(username),  &
     &       trim(database), int2str(IMAS_MAJOR_VERSION), '', &
     &       uri_dest, status )
        end if
        write(systemarg,'(a,a,a,a)')                     &
     &    'idscp --setDatasetVersion'//                  &
     &       ' -s ',trim(uri_source),                    &
     &       ' -d ',trim(uri_dest)
#else
#if ( IMAS_MINOR_VERSION > 31 || IMAS_MAJOR_VERSION > 3 )
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
        if (database.eq.'iter'.or.index(new_path,'imasdb/iter').gt.0) &
     &   systemarg = trim(systemarg)//' -do ITER'
#endif
#ifdef NAGFOR
        call system(systemarg, status, ierror)
#else
        call system(systemarg)
#endif
        if (database.eq.'iter') database = 'ITER'
        if (index(new_path,'imasdb/iter').gt.0) then
          l=index(new_path,'imasdb/iter')
          write(new_path(l+7:l+10),'(a4)') 'ITER'
        end if
#if AL_MAJOR_VERSION > 4
        uri = 'imas:mdsplus?path='//trim(new_path)
        call imas_open( uri, FORCE_CREATE_PULSE, idx, status, message )
        call xertst ( status.eq.0, trim(message) )
#else
        call imas_open_env(treename, shot, new_run, idx, &
          &                username, database, version, status )
        call xertst( status.eq.0, 'Error recreating IDS with new DD version !')
#endif
      else if (.not.same_run_number) then
        call close_ual(idx)
#if AL_MAJOR_VERSION > 4
        uri = 'imas:mdsplus?path='//trim(new_path)
        call imas_open( uri, OPEN_PULSE, idx, status, message )
#else
        call imas_open_env(treename, shot, new_run, idx, &
          &                username, database, version, status)
#endif
        if ( status.ne.0 .or. idx.eq.0 .or. database.eq.'iter' .or. &
          &  index(new_path,'imasdb/iter').gt.0) then ! New run IDS must be created
          idx = 0
          if (database.eq.'iter') database = 'ITER'
          if (index(new_path,'imasdb/iter').gt.0) then
            l=index(new_path,'imasdb/iter')
            write(new_path(l+7:l+10),'(a4)') 'ITER'
          end if
#if AL_MAJOR_VERSION > 4
          uri = 'imas:mdsplus?path='//trim(new_path)
          call imas_open( uri, FORCE_CREATE_PULSE, idx, status, message )
          call xertst( status.eq.0, trim(message) )
#else
          call imas_create_env(treename, shot, new_run, 0, 0, idx, &
            &                  username, database, version, status )
          call xertst( status.eq.0, 'Error creating new IDS file !')
#endif
        end if
      end if
    else
      write(0,*) "No previous IMAS data entry found, a new one will be created"
      idx = 0
      if (database.eq.'iter') database = 'ITER'
      if (index(new_path,'imasdb/iter').gt.0) then
        l=index(new_path,'imasdb/iter')
        write(new_path(l+7:l+10),'(a4)') 'ITER'
      end if
#if AL_MAJOR_VERSION > 4
      uri = 'imas:mdsplus?path='//trim(new_path)
      call imas_open( uri, FORCE_CREATE_PULSE, idx, status, message )
      call xertst( status.eq.0, trim(message) )
#else
      call imas_create_env(treename, shot, run, 0, 0, idx, &
        &                  username, database, version, status )
      call xertst( status.eq.0, 'Error creating new IDS !')
#endif
    end if
    !! Create/Write the set data to IDSs
    call B25_process_ids( geo, mpg, state, state_ext, state_avg, switch, &
      &  edge_profiles, edge_sources, edge_transport, &
      &  radiation, description, equilibrium, &
#if ( IMAS_MINOR_VERSION > 21 || IMAS_MAJOR_VERSION > 3 )
      &  summary, &
#endif
#if ( IMAS_MINOR_VERSION > 25 && IMAS_MINOR_VERSION < 34 && IMAS_MAJOR_VERSION == 3 )
      &  numerics, &
#endif
#if ( IMAS_MINOR_VERSION > 25 || IMAS_MAJOR_VERSION > 3 )
     &   run_start_time, run_end_time, &
#endif
#if ( IMAS_MINOR_VERSION > 30 || IMAS_MAJOR_VERSION > 3 )
      &  divertors, &
#endif
      &  tim, dtim, shot, new_run, database, version, new_eq_ggd )

    write(*,*) "START new_ids_edge"
    call new_ids_edge( edge_profiles, edge_sources, edge_transport, &
        &   radiation, description, equilibrium, &
#if ( IMAS_MINOR_VERSION > 21 || IMAS_MAJOR_VERSION > 3 )
        &   summary, &
#endif
#if ( IMAS_MINOR_VERSION > 25 && IMAS_MINOR_VERSION < 34 && IMAS_MAJOR_VERSION == 3 )
        &   numerics, &
#endif
#if ( IMAS_MINOR_VERSION > 30 || IMAS_MAJOR_VERSION > 3 )
        &   divertors, &
#endif
#if IMAS_MAJOR_VERSION > 3
        &   uri, &
#endif
        &   idx, new_eq_ggd )
    systemarg = 'create_db_entry -u '//trim(username)//' -d '//trim(database) &
        &  //' -s '//trim(shot_string)//' -r '//trim(new_run_string)
    write(0,*) trim(systemarg)
#ifdef NAGFOR
    call system(systemarg, status, ierror)
#else
    call system(systemarg)
#endif
    if (.not.same_run_number) then
! Add superceding information to .yaml file
      systemarg = 'IDS_yaml_replace '//trim(shot_string)//' '// &
        &  trim(run_string)//' '//trim(new_run_string)
      write(0,*) trim(systemarg)
#ifdef NAGFOR
      call system(systemarg, status, ierror)
#else
      call system(systemarg)
#endif
    end if
    call dealloc_ids_edge( edge_profiles, edge_sources, edge_transport, &
#if ( IMAS_MINOR_VERSION > 25 && IMAS_MINOR_VERSION < 34 && IMAS_MAJOR_VERSION == 3 )
        &   numerics, &
#endif
#if ( IMAS_MINOR_VERSION > 30 || IMAS_MAJOR_VERSION > 3 )
        &   divertors, &
#endif
        &   radiation )
    call dealloc_batch_edge( batch_profiles, batch_sources, &
#if ( IMAS_MINOR_VERSION > 21 || IMAS_MAJOR_VERSION > 3 )
        &   summary, &
#endif
        &   description )
    call close_ual(idx)

end program b2_ual_rewrite

!!!Local Variables:
!!! mode: f90
!!! End:
