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
!!      @note   More on the b2_ual writers is available in SOLPS-GUI
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
!!      integer num_time_slices, time_slice_index
!!      real(IDS_real) :: old_start_time, old_end_time, ids_end_time
!!      logical continued
!!
!!-----------------------------------------------------------------------------

program b2_ual_write

    use b2mod_main &
     & , only : b2mn_init
    use b2mod_driver &
     & , only : idx, ids_path, dtim, continued, &
     &          shot, run, username, database, version, &
     &          old_description, old_edge_profiles, equilibrium, &
     &          old_imas_version, old_start_time, old_end_time, &
     &          description, imas_version, ids_end_time, new_eq_ggd, &
     &          edge_profiles, edge_sources, edge_transport, radiation, &
     &          batch_profiles, batch_sources
    use b2mod_time
    use b2mod_switches
    use b2mod_version
    use b2mod_solpstop
    use b2mod_grid_mapping
    use b2mod_numerics_namelist
    use b2us_io
    use b2us_geo
    use b2us_map
    use b2us_plasma
    use b2mod_ual    &
     & , only : put_ids_edge, dealloc_ids_edge, dealloc_batch_edge, &
     &          b25_process_ids, close_ual, &
     &          ids_edge_profiles, ids_edge_sources, ids_edge_transport, &
     &          ids_radiation, ids_dataset_description, ids_equilibrium
    use b2mod_ual_io
#if ( IMAS_MINOR_VERSION > 21 || IMAS_MAJOR_VERSION > 3 )
    use b2mod_ual    &
     & , only : ids_summary
    use b2mod_driver &
     & , only : summary
#endif
#if ( IMAS_MINOR_VERSION > 30 || IMAS_MAJOR_VERSION > 3 )
    use b2mod_ual    &
     & , only : ids_divertors
    use b2mod_driver &
     & , only : divertors
#endif
#if AL_MAJOR_VERSION > 4
    use ids_routines &  ! IGNORE
     & , only : imas_open, al_build_uri_from_legacy_parameters, &
     &          OPEN_PULSE, STRMAXLEN, MDSPLUS_BACKEND
#else
    use ids_routines &  ! IGNORE
     & , only : imas_open_env
    use b2mod_driver &
     & , only : treename
#endif
    use ids_routines &  ! IGNORE
     & , only : ids_get, ids_deallocate
    use ids_schemas  &  ! IGNORE
     & , only : IDS_real, IDS_REAL_INVALID
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
    character(len=24) :: argName
    integer narg, cptArg, l, m, status, tmp_run, time_slice_index
#if AL_MAJOR_VERSION > 4
    character*256 olddir
    character(len=STRMAXLEN) :: uri, uri_source, uri_dest
    character(len=:), allocatable :: message
#endif
    integer idum(0:2)
    character*16 usrnam, run_user
    character*256 imasdir, home_dir, path, systemarg
    logical absolute_path, not_default
    external usrnam

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
    not_default = not_default.or. &
       &         (.not.streql(database, device_env).and. &
       &          .not.streql(database,'solps-iter'))
#if AL_MAJOR_VERSION > 4
    imasdir = trim(home_dir)//'/public/imasdb/'//trim(database)//'/'//trim(version)
#else
    treename = 'ids'
    imasdir = trim(home_dir)//'/public/imasdb/'//trim(database)//'/'//trim(version)//'/0'
#endif
#ifdef NO_GETENV
    write(imas_version,'(i1,a1,i2,a1,i1)') IMAS_MAJOR_VERSION,'.', &
                                         & IMAS_MINOR_VERSION,'.', &
                                         & IMAS_MICRO_VERSION
#else
#ifdef USE_PXFGETENV
    CALL PXFGETENV ('IMASDIR', 0, imasdir, lenval, ierror)
    CALL PXFGETENV ('IMAS_VERSION', 0, imas_version, lenval, ierror)
#else
    call get_environment_variable('IMASDIR', status=ierror, length=lenval)
    if (ierror.eq.0) call get_environment_variable('IMASDIR',value=imasdir)
    call get_environment_variable('IMAS_VERSION', status=ierror, length=lenval)
    if (ierror.eq.0) call get_environment_variable('IMAS_VERSION', value=imas_version)
#endif
#endif
    write(*,*) 'Starting b2mn init'
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
    call read_b2fplasma(56, mpg%nCv, mpg%nFc, state%pl%ns, state)

    call ipgeti('b2mndr_pulse_number', shot )
    if (shot.eq.0) call ipgeti('b2mndr_shot_number', shot )
    call ipgeti('b2mndr_run_number', run )
    path = ' '
    absolute_path = .false.
#if AL_MAJOR_VERSION > 4
    call ipgetc('b2mndr_ids_path', path )
    not_default = not_default.or..not.streql(path,' ')
#endif
    call strip_spaces(path)
    if (.not.streql(path,' ')) absolute_path = path(1:1).eq.'/'

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
          l=index(path,'$IMASDIR')
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
#else
    call xertst( 0.lt.shot.and.shot.le.214748, 'Invalid shot number')
    call xertst( 0.le.run.and.run.le.99999, 'Invalid run number')
    call xertst( .not.streql(username,' '), 'User name not defined !')
#endif
    call xertst( .not.streql(database,' '), 'Database not defined !')
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
#else
      ids_path = imasdir
#endif
    else if (.not.absolute_path) then
      ids_path = trim(imasdir)//'/'//trim(path)
    else
      ids_path = path
    end if

#if AL_MAJOR_VERSION > 4
    write(0,'(6a)') 'Path: ', trim(ids_path), &
        & ' User: ', trim(username), ' Database: ', trim(database)
#else
    write(0,'(a,i8,a,i8,4a)') 'Shot: ', shot, ' Run: ', run, &
        & ' User: ', trim(username), ' Database: ', trim(database)
#endif

    !! Process B2.5 data and set it to IMAS IDS
    write(*,*) "START B25_process_ids"
#if AL_MAJOR_VERSION > 4
    uri = 'imas:mdsplus?path='//trim(ids_path)
    write(0,'(2a)') "Checking if IMAS data entry already exists : ", &
      &  trim(imasdir)//'/'//trim(ids_path)
    call imas_open( uri, OPEN_PULSE, idx, status, message )
#else
    write(0,'(2a,2i8)') "Checking if IMAS data entry already exists : ", &
      &  trim(database), shot, run
    call imas_open_env(treename, shot, run, idx, username, database, version, status)
#endif
    if (status.ne.0) then
      if (database.eq.'ITER'.or.index(ids_path,'imasdb/ITER').gt.0) then
        write(0,*) "Did not find ITER database IDS file."
        write(0,*) "Checking if old ''iter'' case exists."
#if AL_MAJOR_VERSION > 4
        l=index(ids_path,'imasdb/ITER')
        write(olddir,'(a)') ids_path(1:l-1)//'imasdb/iter'//ids_path(l+11:256)
        uri = 'imas:mdsplus?path='//trim(olddir)
        call imas_open( uri, OPEN_PULSE, idx, status, message )
#else
        call imas_open_env(treename, shot, run, idx, username, 'iter', version, status)
#endif
        if (status.eq.0) then
          database = 'iter'
          write(0,*) "Old database case found."
          write(0,*) "Will be rewritten in new location."
        end if
      else if (database.eq.'iter'.or.index(ids_path,'imasdb/iter').gt.0) then
#if AL_MAJOR_VERSION > 4
        l=index(ids_path,'imasdb/iter')
        write(olddir,'(a)') ids_path(1:l-1)//'imasdb/ITER'//ids_path(l+11:256)
        uri = 'imas:mdsplus?path='//trim(olddir)
        call imas_open( uri, OPEN_PULSE, idx, status, message )
#else
        call imas_open_env(treename, shot, run, idx, username, 'ITER', version, status)
#endif
        database = 'ITER'
      end if
    end if
    !! If this is a time continuation run, append the new data to the IDS
    if ( status.eq.0 .and. idx.ne.0 ) then
      write(0,*) "Reading old IMAS data entry: ", trim(database), shot, run
      call ids_get( idx, "equilibrium", equilibrium, status)
      call ids_get( idx, "edge_profiles", old_edge_profiles, status)
      if ( status.ne.0 ) then
        write(0,*) 'Error opening old edge_profiles IDS ! Will create a new one.'
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
          write(0,*) 'Error opening old dataset_description IDS !'
        else if (associated(old_description%ids_properties% &
                        &   version_put%data_dictionary)) then
#if ( IMAS_MINOR_VERSION > 25 || IMAS_MAJOR_VERSION > 3 )
          old_start_time = old_description%simulation%time_begin
          old_end_time = old_description%simulation%time_end
#endif
          old_imas_version = old_description%ids_properties% &
                          &  version_put%data_dictionary(1)
          call ids_deallocate( old_description )
        end if
        continued = run_start_time.eq.IDS_REAL_INVALID .and. &
           &       (ids_end_time.lt.tim .and. ids_end_time.ne.IDS_REAL_INVALID)
        continued = continued .or. &
           &       (run_start_time.ge.ids_end_time .and. &
           &       (ids_end_time.lt.tim .and. ids_end_time.ne.IDS_REAL_INVALID))
        if (continued.or.database.eq.'iter'.or.index(ids_path,'imasdb/iter').gt.0) then
          if (.not.streql(old_imas_version,imas_version).or.database.eq.'iter'.or. &
            & index(ids_path,'imasdb/iter').gt.0) then
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
!xpb Copy the IDS to a temporary location with the new DD and then bring it back
            tmp_run = run
            if (database.ne.'iter'.and.index(ids_path,'imasdb/iter').eq.0) &
              & tmp_run = run + 1000
#if AL_MAJOR_VERSION > 4
            call al_build_uri_from_legacy_parameters                 &
              & ( MDSPLUS_BACKEND, shot, run, trim(username),        &
              &   trim(database), int2str(IMAS_MAJOR_VERSION), '',   &
              &   uri_source, status )
            if (database.eq.'iter'.or.index(ids_path,'imasdb/iter').gt.0) then
              call al_build_uri_from_legacy_parameters               &
                & ( MDSPLUS_BACKEND, shot, tmp_run, trim(username),  &
                &   'ITER', int2str(IMAS_MAJOR_VERSION), '',         &
                &   uri_dest, status )
            else
              call al_build_uri_from_legacy_parameters               &
                & ( MDSPLUS_BACKEND, shot, tmp_run, trim(username),  &
                &   trim(database), int2str(IMAS_MAJOR_VERSION), '', &
                &   uri_dest, status )
            end if
#if IMAS_MAJOR_VERSION < 4
            write(systemarg,'(a,a,a,a,a)')                    &
              & 'idscp --set-dataset-version'//               &
              &      ' -s ',trim(uri_source),                 &
              &      ' -d ',trim(uri_dest),' --dd-update'
#else
            write(systemarg,'(a,a,a,a,a)')                    &
              & 'idscp -s ',trim(uri_source),                 &
              &      ' -d ',trim(uri_dest),' --dd-update'
#endif
#else
#if ( IMAS_MINOR_VERSION > 31 || IMAS_MAJOR_VERSION > 3 )
            write(systemarg,'(a,i7,a,i4,a,i7,a,i4,a,a,a,a)')  &
              & 'idscp --setDatasetVersion'//                 &
              &      ' -si ',shot,' -ri ',run,                &
              &      ' -so ',shot,' -ro ',tmp_run,            &
              &      ' -d ',trim(database),' -u ',trim(username)
#else
            write(systemarg,'(a,i7,a,i4,a,i7,a,i4,a,a,a,a)')  &
              & 'idscp -si ',shot,' -ri ',run,                &
              &      ' -so ',shot,' -ro ',tmp_run,            &
              &      ' -d ',trim(database),' -u ',trim(username)
#endif
            if (database.eq.'iter'.or.index(ids_path,'imasdb/iter').gt.0) &
              & systemarg = trim(systemarg)//' -do ITER'
#endif
#ifdef NAGFOR
            call system(systemarg, status, ierror)
#else
            call system(systemarg)
#endif
            if (database.ne.'iter'.and.index(ids_path,'imasdb/iter').eq.0) then
#if AL_MAJOR_VERSION > 4
              call al_build_uri_from_legacy_parameters               &
                & ( MDSPLUS_BACKEND, shot, tmp_run, trim(username),  &
                &   trim(database), int2str(IMAS_MAJOR_VERSION), '', &
                &   uri_source, status )
              call al_build_uri_from_legacy_parameters               &
                & ( MDSPLUS_BACKEND, shot, run, trim(username),      &
                &   trim(database), int2str(IMAS_MAJOR_VERSION), '', &
                &   uri_dest, status )
              write(systemarg,'(a,a,a,a)')                     &
                & 'idscp -s ',trim(uri_source),                &
                &      ' -d ',trim(uri_dest)
#else
#if ( IMAS_MINOR_VERSION > 31 || IMAS_MAJOR_VERSION > 3 )
              write(systemarg,'(a,i7,a,i4,a,i7,a,i4,a,a,a,a)') &
                & 'idscp --setDatasetVersion'//                &
                &      ' -si ',shot,' -ri ',tmp_run,           &
                &      ' -so ',shot,' -ro ',run,               &
                &      ' -d ',trim(database),' -u ',trim(username)
#else
              write(systemarg,'(a,i7,a,i4,a,i7,a,i4,a,a,a,a)') &
                & 'idscp -si ',shot,' -ri ',tmp_run,           &
                &      ' -so ',shot,' -ro ',run,               &
                &      ' -d ',trim(database),' -u ',trim(username)
#endif
#endif
#ifdef NAGFOR
              call system(systemarg, status, ierror)
#else
              call system(systemarg)
#endif
            end if
            if (database.eq.'iter') database = 'ITER'
            if (index(ids_path,'imasdb/iter').gt.0) then
              l=index(ids_path,'imasdb/iter')
              write(ids_path(l+7:l+10),'(a4)') 'ITER'
            end if
#if AL_MAJOR_VERSION > 4
            uri = 'imas:mdsplus?path='//trim(ids_path)
            call imas_open( uri, OPEN_PULSE, idx, status, message )
#else
            call imas_open_env(treename, shot, run, idx, &
             &                 username, database, version, status)
#endif
          end if
          if (continued) then
            write(0,*) "Appending a new time slice at t = ", tim, " s."
            num_time_slices = num_time_slices + 1
          end if
          time_slice_index = num_time_slices
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
            &   old_start_time, run_end_time, &
#endif
#if ( IMAS_MINOR_VERSION > 30 || IMAS_MAJOR_VERSION > 3 )
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
      write(0,*) "No previous IMAS data entry found, a new one will be created"
      idx = 0
      if (database.eq.'iter') database = 'ITER'
    end if
    if ( status.ne.0 .or. idx.eq.0 ) then
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
         &  run_start_time, run_end_time, &
#endif
#if ( IMAS_MINOR_VERSION > 30 || IMAS_MAJOR_VERSION > 3 )
         &  divertors, &
#endif
         &  tim, dtim, shot, run, database, version, new_eq_ggd )
    end if

    !! Create/Write the set data to IDSs
    write(*,*) "START put_ids_edge"
    call put_ids_edge( edge_profiles, edge_sources, edge_transport, &
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
#if AL_MAJOR_VERSION > 4
        &   ids_path, &
#else
        &   treename, shot, run, username, database, version, &
#endif
        &   idx, new_eq_ggd )
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

end program b2_ual_write

!!!Local Variables:
!!! mode: f90
!!! End:
