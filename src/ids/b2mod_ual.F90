!!-----------------------------------------------------------------------------
!! DOCUMENTATION:
!>      @section b2mod_ual_desc  Description
!!      Module providing Basic UAL routines shared by the entire SOLPS
!!      application chain.
!!
!!      @note   For IMAS: IMAS GGD (GSL) library provides similar and greater
!!              number of routines that are being used for the same purposes.
!!
!!-----------------------------------------------------------------------------
module b2mod_ual

    use b2mod_types
#ifdef IMAS
    use ids_routines &  ! IGNORE
     & , only : ids_deallocate, ids_put, ids_get, ids_delete, ids_put_slice, &
     &          CLOSE_PULSE
#if AL_MAJOR_VERSION > 4
    use ids_routines &  ! IGNORE
     & , only : imas_open, imas_close, al_build_uri_from_legacy_parameters, &
     &          HDF5_BACKEND, MDSPLUS_BACKEND, &
     &          FORCE_CREATE_PULSE, OPEN_PULSE, STRMAXLEN
    use ids_schemas &     ! IGNORE
     & , only : ids_string_length
#elif AL_MAJOR_VERSION == 4
    use ids_routines &  ! IGNORE
     & , only : imas_open_env, imas_create_env, &
     &          ual_begin_pulse_action, ual_open_pulse, ual_close_pulse
# if AL_MINOR_VERSION > 8
    use ids_routines &  ! IGNORE
     & , only : HDF5_BACKEND, FORCE_CREATE_PULSE, OPEN_PULSE
# endif
#else
    use ids_routines &  ! IGNORE
     & , only : imas_open_env, imas_create_env, &
     &          imas_open
#endif
    use ids_schemas &   ! IGNORE
     & , only : ids_edge_profiles, ids_edge_sources, ids_edge_transport, &
     &          ids_radiation, ids_equilibrium
#if ( IMAS_MAJOR_VERSION < 4 || ( IMAS_MAJOR_VERSION == 4 && IMAS_MINOR_VERSION < 1 ) )
    use ids_schemas &   ! IGNORE
     & , only : ids_dataset_description
#endif
#if ( IMAS_MINOR_VERSION > 21 || IMAS_MAJOR_VERSION > 3 )
    use ids_schemas &   ! IGNORE
     & , only : ids_summary
#endif
#if ( IMAS_MINOR_VERSION > 25 && IMAS_MINOR_VERSION < 34 && IMAS_MAJOR_VERSION == 3 )
    use ids_schemas &   ! IGNORE
     & , only : ids_numerics
#endif
#if ( IMAS_MINOR_VERSION > 30 || IMAS_MAJOR_VERSION > 3 )
    use ids_schemas &   ! IGNORE
     & , only : ids_divertors
#endif
#if IMAS_MAJOR_VERSION > 3
    use ids_schemas &   ! IGNORE
     & , only : ids_plasma_profiles, ids_plasma_sources, ids_plasma_transport
#endif
    use b2mod_ual_io &
     & , only : b25_process_ids
#elif defined(ITM_ENVIRONMENT_LOADED)
    use euITM_schemas   ! IGNORE
    use euITM_routines  ! IGNORE
#endif

  implicit none

  private

  public open_ual, close_ual
#ifdef IMAS
  public put_ids_edge, new_ids_edge, delete_ids_edge
  public dealloc_ids_edge, dealloc_batch_edge
  public put_batch_edge, new_batch_edge
  public b25_process_ids, read_ids
  public ids_edge_profiles, ids_edge_sources, ids_edge_transport, &
    &    ids_radiation, ids_equilibrium
#if ( IMAS_MAJOR_VERSION < 4 || ( IMAS_MAJOR_VERSION == 4 && IMAS_MINOR_VERSION < 1 ) )
  public ids_dataset_description
#endif
#if ( IMAS_MINOR_VERSION > 21 || IMAS_MAJOR_VERSION > 3 )
  public ids_summary
#endif
#if ( IMAS_MINOR_VERSION > 25 && IMAS_MINOR_VERSION < 34 && IMAS_MAJOR_VERSION == 3 )
  public ids_numerics
#endif
#if ( IMAS_MINOR_VERSION > 30 || IMAS_MAJOR_VERSION > 3 )
  public ids_divertors
#endif
#endif

contains

#ifdef IMAS
    !> Subroutine used to put data to edge_profiles, edge_sources and
    !! edge_transport IDSs.
    subroutine put_ids_edge( &
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
#if ( IMAS_MINOR_VERSION > 30 || IMAS_MAJOR_VERSION > 3 )
            &   divertors, &
#endif
#if AL_MAJOR_VERSION > 4
            &   ids_path, &
#else
            &   treename, shot, run, username, database, version, &
#endif
            &   idx, new_eq_ggd )
        type (ids_edge_profiles), intent(inout) :: edge_profiles   !< IDS
            !< designed to store data on edge plasma profiles (includes the
            !< scrape-off layer and possibly part of the confined plasma)
        type (ids_edge_sources), intent(inout) :: edge_sources     !< IDS
            !< designed to store data on edge plasma sources. Energy terms
            !< correspond to the full kinetic energy equation (i.e. the energy
            !< flux takes into account the energy transported by the particle
            !< flux)
        type (ids_edge_transport), intent(inout) :: edge_transport !< IDS
            !< designed to store data on edge plasma transport. Energy terms
            !< correspond to the full kinetic energy equation (i.e. the energy
            !< flux takes into account the energy transported by the particle
            !< flux)
#if IMAS_MAJOR_VERSION > 3
        type (ids_plasma_profiles), intent(inout) :: plasma_profiles   !< IDS
            !< designed to store data on plasma profiles
        type (ids_plasma_sources), intent(inout) :: plasma_sources     !< IDS
            !< designed to store data on plasma sources. Energy terms
            !< correspond to the full kinetic energy equation (i.e. the energy
            !< flux takes into account the energy transported by the particle
            !< flux)
        type (ids_plasma_transport), intent(inout) :: plasma_transport !< IDS
            !< designed to store data on plasma transport. Energy terms
            !< correspond to the full kinetic energy equation (i.e. the energy
            !< flux takes into account the energy transported by the particle
            !< flux)
#endif
        type (ids_radiation), intent(inout) :: radiation !< IDS
            !< designed to store data about plasma radiation
#if ( IMAS_MAJOR_VERSION < 4 || ( IMAS_MAJOR_VERSION == 4 && IMAS_MINOR_VERSION < 1 ) )
        type (ids_dataset_description), intent(inout) :: description !< IDS
            !<  designed to store a description of the simulation
#endif
        type (ids_equilibrium), intent(inout) :: equilibrium !< IDS
            !< designed to store a description of the equilibrium
#if ( IMAS_MINOR_VERSION > 21 || IMAS_MAJOR_VERSION > 3 )
        type (ids_summary), intent(inout) :: summary !< IDS
            !< designed to store run summary data
#endif
#if ( IMAS_MINOR_VERSION > 25 && IMAS_MINOR_VERSION < 34 && IMAS_MAJOR_VERSION == 3 )
        type (ids_numerics), intent(inout) :: numerics !< IDS designed to store
            !< run numerics data
#endif
#if ( IMAS_MINOR_VERSION > 30 || IMAS_MAJOR_VERSION > 3 )
        type (ids_divertors), intent(inout) :: divertors !< IDS
            !< designed to store run data related to the divertor plates
#endif
#if AL_MAJOR_VERSION > 4
        character(len=256), intent(in) :: ids_path  !< The path to the IMAS data entry
#else
        character(len=24), intent(in) :: treename   !< The name of the IMAS IDS database
        integer, intent(in) :: shot   !< The shot number of the database being created
        integer, intent(in) :: run    !< The run number of the database being created
        character(len=24), intent(in) :: username   !< Creator/owner of the IMAS IDS
            !< database
        character(len=24), intent(in) :: database   !< IMAS database name
            !< (i. e. solps-iter, ITER, aug)
        character(len=24), intent(in) :: version    !< Major version of the IMAS IDS
#endif
        integer, intent(inout) :: idx !< The returned identifier to be used in the
            !< subsequent data access operation
        logical, intent(in) :: new_eq_ggd
            !< database
        integer :: status
        character*256 :: ids_list
#if AL_MAJOR_VERSION > 4
        character(len=:), allocatable :: message
        character(len=STRMAXLEN) :: uri
        logical, save :: first_pass = .true.
#endif

            !< procedures
        external xertst, xerrab

        ids_list = "edge_profiles, edge_sources, edge_transport"
#if IMAS_MAJOR_VERSION > 3
        ids_list = trim(ids_list)// &
          & ", plasma_profiles, plasma_sources, plasma_transport"
#endif
#if ( IMAS_MINOR_VERSION > 21 || IMAS_MAJOR_VERSION > 3 )
        ids_list = trim(ids_list)//", summary"
#endif
#if ( IMAS_MINOR_VERSION > 25 && IMAS_MINOR_VERSION < 34 && IMAS_MAJOR_VERSION == 3 )
        ids_list = trim(ids_list)//", numerics"
#endif
#if ( IMAS_MINOR_VERSION > 30 || IMAS_MAJOR_VERSION > 3 )
        ids_list = trim(ids_list)//", divertors"
#endif
#if ( IMAS_MAJOR_VERSION < 4 || ( IMAS_MAJOR_VERSION == 4 && IMAS_MINOR_VERSION < 1 ) )
#if AL_MAJOR_VERSION > 4
        if (first_pass) ids_list = trim(ids_list)//", dataset_description"
#else
        ids_list = trim(ids_list)//", dataset_description"
#endif
#endif
        ids_list = trim(ids_list)//", and radiation"
        !! Set data to edge_profiles IDS
        write(*,'(1x,a)') "Writing "//trim(ids_list)//" IDS"

        !! Create and modify new shot/run
        if ( idx.eq.0 ) then
#if AL_MAJOR_VERSION > 4
          uri = 'imas:mdsplus?path='//trim(ids_path)
#if ( IMAS_MAJOR_VERSION < 4 || ( IMAS_MAJOR_VERSION == 4 && IMAS_MINOR_VERSION < 1 ) )
#if IMAS_MAJOR_VERSION > 3
          allocate( description%uri(1) )
          description%uri = trim(uri)
#endif
#else
          allocate( summary%identifier(1) )
          summary%identifier = trim(uri)
#endif
          call imas_open( uri, FORCE_CREATE_PULSE, idx, status, message )
#else
          call imas_create_env( treename, shot, run, 0, 0, idx, username, &
             & database, version, status )
#endif
          if (status.ne.0) then
            write(0,*) 'Opening IMAS database failed !'
#if AL_MAJOR_VERSION > 4
            write(0,*) 'Make sure the IDS path directory exists.'
            write(0,*) 'IDS path requested is : '//trim(ids_path)
            call xerrab( trim(message) )
#else
            write(0,*) 'Make sure it exists or create it with the command:'
            write(0,*) 'imasdb '//trim(database)
            call xerrab( 'Error opening IMAS database !')
#endif
          endif

        !! Put data to IDS
          write(*,*) 'Putting edge_profiles IDS'
          call ids_put( idx, "edge_profiles", edge_profiles, status )
          call xertst( status.eq.0, 'Error putting edge_profiles IDS !')
          write(*,*) 'Putting egde_sources IDS'
          call ids_put( idx, "edge_sources", edge_sources, status )
          call xertst( status.eq.0, 'Error putting edge_sources IDS !')
          write(*,*) 'Putting edge_transport IDS'
          call ids_put( idx, "edge_transport", edge_transport, status )
          call xertst( status.eq.0, 'Error putting edge_transport IDS !')
#if IMAS_MAJOR_VERSION > 3
          write(*,*) 'Putting plasma_profiles IDS'
          call ids_put( idx, "plasma_profiles", plasma_profiles, status )
          call xertst( status.eq.0, 'Error putting plasma_profiles IDS !')
          write(*,*) 'Putting plasma_sources IDS'
          call ids_put( idx, "plasma_sources", plasma_sources, status )
          call xertst( status.eq.0, 'Error putting plasma_sources IDS !')
          write(*,*) 'Putting plasma_transport IDS'
          call ids_put( idx, "plasma_transport", plasma_transport, status )
          call xertst( status.eq.0, 'Error putting plasma_transport IDS !')
#endif
          write(*,*) 'Putting radiation IDS'
          call ids_put( idx, "radiation", radiation, status )
          call xertst( status.eq.0, 'Error putting radiation IDS !')
#if ( IMAS_MAJOR_VERSION < 4 || ( IMAS_MAJOR_VERSION == 4 && IMAS_MINOR_VERSION < 1 ) )
          write(*,*) 'Putting dataset_description IDS'
          call ids_put( idx, "dataset_description", description, status )
          call xertst( status.eq.0, 'Error putting dataset_description IDS !')
#endif
#if ( IMAS_MINOR_VERSION > 21 || IMAS_MAJOR_VERSION > 3 )
          write(*,*) 'Putting summary IDS'
          call ids_put( idx, "summary", summary, status )
          call xertst( status.eq.0, 'Error putting summary IDS !')
#endif
#if ( IMAS_MINOR_VERSION > 25 && IMAS_MINOR_VERSION < 34 && IMAS_MAJOR_VERSION == 3 )
          write(*,*) 'Putting numerics IDS'
          call ids_put( idx, "numerics", numerics, status )
          call xertst( status.eq.0, 'Error putting numerics IDS !')
#endif
#if ( IMAS_MINOR_VERSION > 30 || IMAS_MAJOR_VERSION > 3 )
          write(*,*) 'Putting divertors IDS'
          call ids_put( idx, "divertors", divertors, status )
          call xertst( status.eq.0, 'Error putting divertors IDS !')
#endif
        else
        !! Or open and modify existing shot/run
        !! (might work much faster than imas_create_env)

          if (new_eq_ggd) then
            write(*,'(1x,a)') "Adding GGD data to equilibrium IDS"
            call ids_put_slice( idx, "equilibrium", equilibrium, status )
            call xertst( status.eq.0, 'Error putting slice in equilibrium IDS !')
          end if

        !! Put data to IDS
          write(*,*) 'Putting edge_profiles IDS slice'
          call ids_put_slice( idx, "edge_profiles", edge_profiles, status )
          call xertst( status.eq.0, 'Error putting slice in edge_profiles IDS !')
          write(*,*) 'Putting edge_sources IDS slice'
          call ids_put_slice( idx, "edge_sources", edge_sources, status )
          call xertst( status.eq.0, 'Error putting slice in edge_sources IDS !')
          write(*,*) 'Putting edge_transport IDS slice'
          call ids_put_slice( idx, "edge_transport", edge_transport, status )
          call xertst( status.eq.0, 'Error putting slice in edge_transport IDS !')
#if IMAS_MAJOR_VERSION > 3
          write(*,*) 'Putting plasma_profiles IDS slice'
          call ids_put_slice( idx, "plasma_profiles", plasma_profiles, status )
          call xertst( status.eq.0, 'Error putting slice in plasma_profiles IDS !')
          write(*,*) 'Putting plasma_sources IDS slice'
          call ids_put_slice( idx, "plasma_sources", plasma_sources, status )
          call xertst( status.eq.0, 'Error putting slice in plasma_sources IDS !')
          write(*,*) 'Putting plasma_transport IDS slice'
          call ids_put_slice( idx, "plasma_transport", plasma_transport, status )
          call xertst( status.eq.0, 'Error putting slice in plasma_transport IDS !')
#endif
          write(*,*) 'Putting radiation IDS slice'
          call ids_put_slice( idx, "radiation", radiation, status )
          call xertst( status.eq.0, 'Error putting slice in radiation IDS !')
#if ( IMAS_MAJOR_VERSION < 4 || ( IMAS_MAJOR_VERSION == 4 && IMAS_MINOR_VERSION < 1 ) )
#if AL_MAJOR_VERSION > 4
          if (first_pass) then
            write(*,*) 'Putting dataset_description IDS'
            call ids_put( idx, "dataset_description", description, status )
            call xertst( status.eq.0, 'Error putting dataset_description IDS !')
          endif
#else
          write(*,*) 'Putting dataset_description IDS slice'
          call ids_put_slice( idx, "dataset_description", description, status )
          call xertst( status.eq.0, 'Error putting slice in dataset_description IDS !')
#endif
#endif
#if ( IMAS_MINOR_VERSION > 21 || IMAS_MAJOR_VERSION > 3 )
          write(*,*) 'Putting summary IDS slice'
          call ids_put_slice( idx, "summary", summary, status )
          call xertst( status.eq.0, 'Error putting slice in summary IDS !')
#endif
#if ( IMAS_MINOR_VERSION > 25 && IMAS_MINOR_VERSION < 34 && IMAS_MAJOR_VERSION == 3 )
          write(*,*) 'Putting numerics IDS slice'
          call ids_put_slice( idx, "numerics", numerics, status )
          call xertst( status.eq.0, 'Error putting slice in numerics IDS !')
#endif
#if ( IMAS_MINOR_VERSION > 30 || IMAS_MAJOR_VERSION > 3 )
          write(*,*) 'Putting divertors IDS slice'
          call ids_put_slice( idx, "divertors", divertors, status )
          call xertst( status.eq.0, 'Error putting slice in divertors IDS !')
#endif
        end if

#if AL_MAJOR_VERSION > 4
       first_pass = .false.
#endif
        write(*,*) "IDS write finished"
        return

    end subroutine put_ids_edge

    subroutine dealloc_ids_edge( &
            &   edge_profiles, edge_sources, edge_transport, &
#if IMAS_MAJOR_VERSION > 3
            &   plasma_profiles, plasma_sources, plasma_transport, &
#endif
#if ( IMAS_MINOR_VERSION > 25 && IMAS_MINOR_VERSION < 34 && IMAS_MAJOR_VERSION == 3 )
            &   numerics, &
#endif
#if ( IMAS_MINOR_VERSION > 30 || IMAS_MAJOR_VERSION > 3 )
            &   divertors, &
#endif
            &   radiation )
        implicit none
        type (ids_edge_profiles), intent(inout) :: edge_profiles   !< IDS
            !< designed to store data on edge plasma profiles (includes the
            !< scrape-off layer and possibly part of the confined plasma)
        type (ids_edge_sources), intent(inout) :: edge_sources     !< IDS
            !< designed to store data on edge plasma sources. Energy terms
            !< correspond to the full kinetic energy equation (i.e. the energy
            !< flux takes into account the energy transported by the particle
            !< flux)
        type (ids_edge_transport), intent(inout) :: edge_transport !< IDS
            !< designed to store data on edge plasma transport. Energy terms
            !< correspond to the full kinetic energy equation (i.e. the energy
            !< flux takes into account the energy transported by the particle
            !< flux)
#if IMAS_MAJOR_VERSION > 3
        type (ids_plasma_profiles), intent(inout) :: plasma_profiles   !< IDS
            !< designed to store data on plasma profiles
        type (ids_plasma_sources), intent(inout) :: plasma_sources     !< IDS
            !< designed to store data on plasma sources. Energy terms
            !< correspond to the full kinetic energy equation (i.e. the energy
            !< flux takes into account the energy transported by the particle
            !< flux)
        type (ids_plasma_transport), intent(inout) :: plasma_transport !< IDS
            !< designed to store data on plasma transport. Energy terms
            !< correspond to the full kinetic energy equation (i.e. the energy
            !< flux takes into account the energy transported by the particle
            !< flux)
#endif
        type (ids_radiation), intent(inout) :: radiation !< IDS
            !< designed to store data about plasma radiation
#if ( IMAS_MINOR_VERSION > 25 && IMAS_MINOR_VERSION < 34 && IMAS_MAJOR_VERSION == 3 )
        type (ids_numerics), intent(inout) :: numerics !< IDS designed to store
            !< run numerics data
#endif
#if ( IMAS_MINOR_VERSION > 30 || IMAS_MAJOR_VERSION > 3 )
        type (ids_divertors), intent(inout) :: divertors !< IDS
            !< designed to store run data related to the divertor plates
#endif
        if (.not.associated( edge_profiles%ids_properties%comment )) return
        call ids_deallocate( edge_profiles )
        call ids_deallocate( edge_sources )
        call ids_deallocate( edge_transport )
#if IMAS_MAJOR_VERSION > 3
        call ids_deallocate( plasma_profiles )
        call ids_deallocate( plasma_sources )
        call ids_deallocate( plasma_transport )
#endif
        call ids_deallocate( radiation )
#if ( IMAS_MINOR_VERSION > 25 && IMAS_MINOR_VERSION < 34 && IMAs_MAJOR_VERSION == 3 )
        call ids_deallocate( numerics )
#endif
#if ( IMAS_MINOR_VERSION > 30 || IMAS_MAJOR_VERSION > 3 )
        call ids_deallocate( divertors )
#endif
         return
    end subroutine dealloc_ids_edge

    subroutine put_batch_edge( &
#if AL_MAJOR_VERSION > 4
            &   ids_path, &
#else
            &   treename, shot, run, username, database, version, &
#endif
            &   idx, batch_profiles, batch_sources, equilibrium, &
#if ( IMAS_MINOR_VERSION > 21 || IMAS_MAJOR_VERSION > 3 )
            &   summary, &
#endif
#if ( IMAS_MAJOR_VERSION < 4 || ( IMAS_MAJOR_VERSION == 4 && IMAS_MINOR_VERSION < 1 ) )
            &   description, &
#endif
            &   do_summary, new_eq_ggd )
        type (ids_edge_profiles), intent(inout) :: batch_profiles   !< IDS
            !< designed to store data on edge plasma profiles (includes the
            !< scrape-off layer and possibly part of the confined plasma)
        type (ids_edge_sources), intent(inout) :: batch_sources     !< IDS
            !< designed to store data on edge plasma sources. Energy terms
            !< correspond to the full kinetic energy equation (i.e. the energy
            !< flux takes into account the energy transported by the particle
            !< flux)
#if ( IMAS_MAJOR_VERSION < 4 || ( IMAS_MAJOR_VERSION == 4 && IMAS_MINOR_VERSION < 1 ) )
        type (ids_dataset_description), intent(inout) :: description
            !< IDS designed to store a description of the simulation
#endif
        type (ids_equilibrium), intent(inout) :: equilibrium
            !< IDS designed to store a description of the equilibrium
#if ( IMAS_MINOR_VERSION > 21 || IMAS_MAJOR_VERSION > 3 )
        type (ids_summary), intent(inout) :: summary !< IDS
            !< designed to store run summary data
#endif
#if AL_MAJOR_VERSION > 4
        character(len=256), intent(in) :: ids_path  !< The path to the IMAS data entry
#else
        character(len=24), intent(in) :: treename   !< The name of the IMAS IDS database
        integer, intent(in) :: shot   !< The shot number of the database being created
        integer, intent(in) :: run    !< The run number of the database being created
        character(len=24), intent(in) :: username   !< Creator/owner of the IMAS IDS
            !< database
        character(len=24), intent(in) :: database   !< IMAS database name
            !< (i. e. solps-iter, ITER, aug)
        character(len=24), intent(in) :: version    !< Major version of the IMAS IDS
#endif
        integer, intent(inout) :: idx !< The returned identifier to be used in the
            !< subsequent data access operation
        logical, intent(in) :: do_summary, new_eq_ggd
            !< database
        integer :: status
        character*256 :: ids_list
#if AL_MAJOR_VERSION > 4
        character(len=:), allocatable :: message
        character(len=STRMAXLEN) :: uri
        logical, save :: first_pass = .true.
#endif

            !< procedures
        external xertst, xerrab

        !! Set data to edge_profiles IDS
        ids_list = "batch_profiles and batch_sources"
        if (do_summary) then
          ids_list = "batch_profiles, batch_sources"
#if ( IMAS_MINOR_VERSION > 21 || IMAS_MAJOR_VERSION > 3 )
          ids_list = trim(ids_list)//", summary"
#endif
#if ( IMAS_MAJOR_VERSION < 4 || ( IMAS_MAJOR_VERSION == 4 && IMAS_MINOR_VERSION < 1 ) )
#if AL_MAJOR_VERSION > 4
          if (first_pass) ids_list = trim(ids_list)//", and dataset_description"
#else
          ids_list = trim(ids_list)//", and dataset_description"
#endif
#endif
        end if
        write(*,'(1x,a)') "Writing "//trim(ids_list)//" IDS"

        !! Create and modify new shot/run
        if ( idx.eq.0 ) then
#if AL_MAJOR_VERSION > 4
          uri = 'imas:mdsplus?path='//trim(ids_path)
#if IMAS_MAJOR_VERSION > 3
          if (do_summary) then
#if ( IMAS_MAJOR_VERSION < 4 || ( IMAS_MAJOR_VERSION == 4 && IMAS_MINOR_VERSION < 1 ) )
            allocate( description%uri(1) )
            description%uri = trim(uri)
#else
            allocate( summary%identifier(1) )
            summary%identifier = trim(uri)
#endif
          end if
#endif
          call imas_open( uri, FORCE_CREATE_PULSE, idx, status, message )
#else
          call imas_create_env( treename, shot, run, 0, 0, idx, username, &
             & database, version, status )
#endif
          if (status.ne.0) then
            write(0,*) 'Opening IMAS database failed !'
#if AL_MAJOR_VERSION > 4
            write(0,*) 'Make sure the IDS path directory exists.'
            write(0,*) 'IDS path requested is : '//trim(ids_path)
            call xerrab( trim(message) )
#else
            write(0,*) 'Make sure it exists or create it with the command:'
            write(0,*) 'imasdb '//trim(database)
            call xerrab( 'Error opening IMAS database !')
#endif
          endif

        !! Put data to IDS
          call ids_put( idx, "edge_profiles/1", batch_profiles, status )
          call xertst( status.eq.0, 'Error putting batch_profiles IDS !')
          call ids_put( idx, "edge_sources/1", batch_sources, status )
          call xertst( status.eq.0, 'Error putting batch_sources IDS !')
          if (do_summary) then
#if ( IMAS_MAJOR_VERSION < 4 || ( IMAS_MAJOR_VERSION == 4 && IMAS_MINOR_VERSION < 1 ) )
            call ids_put( idx, "dataset_description", description, status )
            call xertst( status.eq.0, 'Error putting dataset_description IDS !')
#endif
#if ( IMAS_MINOR_VERSION > 21 || IMAS_MAJOR_VERSION > 3 )
            call ids_put( idx, "summary", summary, status )
            call xertst( status.eq.0, 'Error putting summary IDS !')
#endif
          end if
        else
        !! Or open and modify existing shot/run
        !! (might work much faster than imas_create_env)

          if (new_eq_ggd) then
            write(*,'(1x,a)') "Adding GGD data to equilibrium IDS"
            call ids_put_slice( idx, "equilibrium", equilibrium, status )
            call xertst ( status.eq.0, 'Error putting slice in equilibrium IDS !')
          end if

        !! Put data to IDS
          call ids_put_slice( idx, "edge_profiles/1", batch_profiles, status )
          call xertst( status.eq.0, 'Error putting slice in batch_profiles IDS !')
          call ids_put_slice( idx, "edge_sources/1", batch_sources, status )
          call xertst( status.eq.0, 'Error putting slice in batch_sources IDS !')
          if (do_summary) then
#if ( IMAS_MAJOR_VERSION < 4 || ( IMAS_MAJOR_VERSION == 4 && IMAS_MINOR_VERSION < 1 ) )
#if AL_MAJOR_VERSION > 4
            if (first_pass) then
              call ids_put( idx, "dataset_description", description, status )
              call xertst( status.eq.0, 'Error putting dataset_description IDS !')
            end if
#else
            call ids_put_slice( idx, "dataset_description", description, status )
            call xertst( status.eq.0, 'Error putting slice in dataset_description IDS !')
#endif
#endif
#if ( IMAS_MINOR_VERSION > 21 || IMAS_MAJOR_VERSION > 3 )
            call ids_put_slice( idx, "summary", summary, status )
            call xertst( status.eq.0, 'Error putting slice in summary IDS !')
#endif
          end if
        end if

#if AL_MAJOR_VERSION > 4
       first_pass = .false.
#endif
        write(*,*) "IDS write finished for batch averages"
        return

    end subroutine put_batch_edge

    subroutine dealloc_batch_edge( &
            &   batch_profiles, batch_sources, &
#if IMAS_MAJOR_VERSION > 3
            &   batch_plasma_profiles, batch_plasma_sources, &
#endif
#if ( IMAS_MINOR_VERSION > 21 || IMAS_MAJOR_VERSION > 3 )
            &   summary, &
#endif
#if ( IMAS_MAJOR_VERSION < 4 || ( IMAS_MAJOR_VERSION == 4 && IMAS_MINOR_VERSION < 1 ) )
            &   description &
#endif
            &   )
        implicit none
        type (ids_edge_profiles), intent(inout) :: batch_profiles   !< IDS
            !< designed to store data on edge plasma profiles (includes the
            !< scrape-off layer and possibly part of the confined plasma)
        type (ids_edge_sources), intent(inout) :: batch_sources     !< IDS
            !< designed to store data on edge plasma sources. Energy terms
            !< correspond to the full kinetic energy equation (i.e. the energy
            !< flux takes into account the energy transported by the particle
            !< flux)
#if IMAS_MAJOR_VERSION > 3
        type (ids_plasma_profiles), intent(inout) :: batch_plasma_profiles
            !< IDS designed to store data on plasma profiles
        type (ids_plasma_sources), intent(inout) :: batch_plasma_sources
            !< IDS designed to store data on plasma sources. Energy terms
            !< correspond to the full kinetic energy equation (i.e. the energy
            !< flux takes into account the energy transported by the particle
            !< flux)
#endif
#if ( IMAS_MAJOR_VERSION < 4 || ( IMAS_MAJOR_VERSION == 4 && IMAS_MINOR_VERSION < 1 ) )
        type (ids_dataset_description), intent(inout) :: description
            !< IDS designed to store a description of the simulation
#endif
#if ( IMAS_MINOR_VERSION > 21 || IMAS_MAJOR_VERSION > 3 )
        type (ids_summary), intent(inout) :: summary !< IDS
            !< designed to store run summary data
#endif
        if (associated( batch_profiles%ids_properties%comment ) ) then
          call ids_deallocate( batch_profiles )
          call ids_deallocate( batch_sources )
#if IMAS_MAJOR_VERSION > 3
          call ids_deallocate( batch_plasma_profiles )
          call ids_deallocate( batch_plasma_sources )
#endif
        end if
#if ( IMAS_MAJOR_VERSION < 4 || ( IMAS_MAJOR_VERSION == 4 && IMAS_MINOR_VERSION < 1 ) )
        if (associated( description%ids_properties%comment ) ) then
          call ids_deallocate( description )
#else
        if (associated( summary%ids_properties%comment ) ) then
#endif
#if ( IMAS_MINOR_VERSION > 21 || IMAS_MAJOR_VERSION > 3 )
          call ids_deallocate( summary )
#endif
        end if
        return
    end subroutine dealloc_batch_edge

    !> Subroutine used to delete data from edge_profiles, edge_sources and
    !! edge_transport IDSs.
    subroutine delete_ids_edge( &
            &   edge_profiles, edge_sources, edge_transport, &
#if IMAS_MAJOR_VERSION > 3
            &   plasma_profiles, plasma_sources, plasma_transport, &
#endif
            &   radiation, &
#if ( IMAS_MAJOR_VERSION < 4 || ( IMAS_MAJOR_VERSION == 4 && IMAS_MINOR_VERSION < 1 ) )
            &   description, &
#endif
#if ( IMAS_MINOR_VERSION > 21 || IMAS_MAJOR_VERSION > 3 )
            &   summary, &
#endif
#if ( IMAS_MINOR_VERSION > 25 && IMAS_MINOR_VERSION < 34 && IMAS_MAJOR_VERSION == 3 )
            &   numerics, &
#endif
#if ( IMAS_MINOR_VERSION > 30 || IMAS_MAJOR_VERSION > 3 )
            &   divertors, &
#endif
            &   idx )
        type (ids_edge_profiles), intent(inout) :: edge_profiles   !< IDS
            !< designed to store data on edge plasma profiles (includes the
            !< scrape-off layer and possibly part of the confined plasma)
        type (ids_edge_sources), intent(inout) :: edge_sources     !< IDS
            !< designed to store data on edge plasma sources. Energy terms
            !< correspond to the full kinetic energy equation (i.e. the energy
            !< flux takes into account the energy transported by the particle
            !< flux)
        type (ids_edge_transport), intent(inout) :: edge_transport !< IDS
            !< designed to store data on edge plasma transport. Energy terms
            !< correspond to the full kinetic energy equation (i.e. the energy
            !< flux takes into account the energy transported by the particle
            !< flux)
#if IMAS_MAJOR_VERSION > 3
        type (ids_plasma_profiles), intent(inout) :: plasma_profiles   !< IDS
            !< designed to store data on plasma profiles
        type (ids_plasma_sources), intent(inout) :: plasma_sources     !< IDS
            !< designed to store data on plasma sources. Energy terms
            !< correspond to the full kinetic energy equation (i.e. the energy
            !< flux takes into account the energy transported by the particle
            !< flux)
        type (ids_plasma_transport), intent(inout) :: plasma_transport !< IDS
            !< designed to store data on plasma transport. Energy terms
            !< correspond to the full kinetic energy equation (i.e. the energy
            !< flux takes into account the energy transported by the particle
            !< flux)
#endif
        type (ids_radiation), intent(inout) :: radiation !< IDS
            !< designed to store data about plasma radiation
#if ( IMAS_MAJOR_VERSION < 4 || ( IMAS_MAJOR_VERSION == 4 && IMAS_MINOR_VERSION < 1 ) )
        type (ids_dataset_description) :: description !< IDS designed to store
            !< a description of the simulation
#endif
#if ( IMAS_MINOR_VERSION > 21 || IMAS_MAJOR_VERSION > 3 )
        type (ids_summary), intent(inout) :: summary !< IDS
            !< designed to store run summary data
#endif
#if ( IMAS_MINOR_VERSION > 25 && IMAS_MINOR_VERSION < 34 && IMAS_MAJOR_VERSION == 3 )
        type (ids_numerics), intent(inout) :: numerics !< IDS designed to store
            !< run numerics data
#endif
#if ( IMAS_MINOR_VERSION > 30 || IMAS_MAJOR_VERSION > 3 )
        type (ids_divertors), intent(inout) :: divertors !< IDS
            !< designed to store data related to divertor plates
#endif
        integer, intent(inout) :: idx !< The returned identifier to be used in the
            !< subsequent data access operation

        if ( idx.ne.0 ) then

        !! Delete data from IDS
          call ids_delete( idx, "edge_profiles", edge_profiles)
          call ids_delete( idx, "edge_sources", edge_sources)
          call ids_delete( idx, "edge_transport", edge_transport)
#if IMAS_MAJOR_VERSION > 3
          call ids_delete( idx, "plasma_profiles", plasma_profiles)
          call ids_delete( idx, "plasma_sources", plasma_sources)
          call ids_delete( idx, "plasma_transport", plasma_transport)
#endif
          call ids_delete( idx, "radiation", radiation)
#if ( IMAS_MAJOR_VERSION < 4 || ( IMAS_MAJOR_VERSION == 4 && IMAS_MINOR_VERSION < 1 ) )
          call ids_delete( idx, "dataset_description", description)
#endif
#if ( IMAS_MINOR_VERSION > 21 || IMAS_MAJOR_VERSION > 3 )
          call ids_delete( idx, "summary", summary)
#endif
#if ( IMAS_MINOR_VERSION > 25 && IMAS_MINOR_VERSION < 34 && IMAS_MAJOR_VERSION == 3 )
          call ids_delete( idx, "numerics", numerics)
#endif
#if ( IMAS_MINOR_VERSION > 30 || IMAS_MAJOR_VERSION > 3 )
          call ids_delete( idx, "divertors", divertors)
#endif
          write(*,*) "IDS delete finished"
        end if
        return

    end subroutine delete_ids_edge

    !> Subroutine used to rewrite data to edge_profiles, edge_sources and
    !! edge_transport IDSs.
    subroutine new_ids_edge( &
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
#if ( IMAS_MINOR_VERSION > 30 || IMAS_MAJOR_VERSION > 3 )
            &   divertors, &
#endif
#if IMAS_MAJOR_VERSION > 3
            &   uri, &
#endif
            &   idx, new_eq_ggd )
        type (ids_edge_profiles), intent(inout) :: edge_profiles   !< IDS
            !< designed to store data on edge plasma profiles (includes the
            !< scrape-off layer and possibly part of the confined plasma)
        type (ids_edge_sources), intent(inout) :: edge_sources     !< IDS
            !< designed to store data on edge plasma sources. Energy terms
            !< correspond to the full kinetic energy equation (i.e. the energy
            !< flux takes into account the energy transported by the particle
            !< flux)
        type (ids_edge_transport), intent(inout) :: edge_transport !< IDS
            !< designed to store data on edge plasma transport. Energy terms
            !< correspond to the full kinetic energy equation (i.e. the energy
            !< flux takes into account the energy transported by the particle
            !< flux)
#if IMAS_MAJOR_VERSION > 3
        type (ids_plasma_profiles), intent(inout) :: plasma_profiles   !< IDS
            !< designed to store data on plasma profiles
        type (ids_plasma_sources), intent(inout) :: plasma_sources     !< IDS
            !< designed to store data on plasma sources. Energy terms
            !< correspond to the full kinetic energy equation (i.e. the energy
            !< flux takes into account the energy transported by the particle
            !< flux)
        type (ids_plasma_transport), intent(inout) :: plasma_transport !< IDS
            !< designed to store data on edge transport. Energy terms
            !< correspond to the full kinetic energy equation (i.e. the energy
            !< flux takes into account the energy transported by the particle
            !< flux)
#endif
        type (ids_radiation), intent(inout) :: radiation !< IDS
            !< designed to store data about plasma radiation
#if ( IMAS_MAJOR_VERSION < 4 || ( IMAS_MAJOR_VERSION == 4 && IMAS_MINOR_VERSION < 1 ) )
        type (ids_dataset_description) :: description !< IDS designed to store
            !< a description of the simulation
#endif
        type (ids_equilibrium) :: equilibrium !< IDS designed to store
            !< a description of the equilibrium
#if ( IMAS_MINOR_VERSION > 21 || IMAS_MAJOR_VERSION > 3 )
        type (ids_summary), intent(inout) :: summary !< IDS
            !< designed to store run summary data
#endif
#if ( IMAS_MINOR_VERSION > 25 && IMAS_MINOR_VERSION < 34 && IMAS_MAJOR_VERSION == 3 )
        type (ids_numerics), intent(inout) :: numerics !< IDS designed to store
            !< run numerics data
#endif
#if ( IMAS_MINOR_VERSION > 30 || IMAS_MAJOR_VERSION > 3 )
        type (ids_divertors), intent(inout) :: divertors !< IDS
            !< designed to store data related to the divertor plates
#endif
#if IMAS_MAJOR_VERSION > 3
        character(len=STRMAXLEN), intent(in) :: uri
#endif
        integer, intent(inout) :: idx !< The returned identifier to be used in the
            !< subsequent data access operation
        logical, intent(in) :: new_eq_ggd
        integer :: status

            !< procedures
        external xertst

        !! Set data to edge_profiles IDS
        write(*,'(1x,a)') "Writing edge_profiles, edge_sources, edge_transport, "// &
#if IMAS_MAJOR_VERSION > 3
          &  "plasma_profiles, plasma_sources, plasma_transport, "// &
#endif
#if ( IMAS_MINOR_VERSION > 21 || IMAS_MAJOR_VERSION > 3 )
          &  "summary, "// &
#endif
#if ( IMAS_MINOR_VERSION > 25 && IMAS_MINOR_VERSION < 34 && IMAS_MAJOR_VERSION == 3 )
          &  "numerics, "// &
#endif
#if ( IMAS_MINOR_VERSION > 30 || IMAS_MAJOR_VERSION > 3 )
          &  "divertors, "// &
#endif
#if ( IMAS_MAJOR_VERSION < 4 || ( IMAS_MAJOR_VERSION == 4 && IMAS_MINOR_VERSION < 1 ) )
          &  "dataset_description, "// &
#endif
          &  "and radiation IDS"

#if ( IMAS_MAJOR_VERSION == 4 && IMAS_MINOR_VERSION < 1 )
        allocate( description%uri(1) )
        description%uri = trim(uri)
#elif ( IMAS_MAJOR_VERSION > 4 || ( IMAS_MAJOR_VERSION == 4 && IMAS_MINOR_VERSION > 0 ) )
        allocate( summary%identifier(1) )
        summary%identifier = trim(uri)
#endif
        if (new_eq_ggd) then
          write(*,'(1x,a)') "Adding GGD data to equilibrium IDS"
          call ids_put( idx, "equilibrium", equilibrium, status )
          call xertst( status.eq.0, 'Error putting equilibrium IDS !')
        end if

        !! Put data to IDS
        call ids_put( idx, "edge_profiles", edge_profiles, status )
        call xertst( status.eq.0, 'Error putting edge_profiles IDS !')
        call ids_put( idx, "edge_sources", edge_sources, status )
        call xertst( status.eq.0, 'Error putting edge_sources IDS !')
        call ids_put( idx, "edge_transport", edge_transport, status )
        call xertst( status.eq.0, 'Error putting edge_transport IDS !')
#if IMAS_MAJOR_VERSION > 3
        call ids_put( idx, "plasma_profiles", plasma_profiles, status )
        call xertst( status.eq.0, 'Error putting plasma_profiles IDS !')
        call ids_put( idx, "plasma_sources", plasma_sources, status )
        call xertst( status.eq.0, 'Error putting plasma_sources IDS !')
        call ids_put( idx, "plasma_transport", plasma_transport, status )
        call xertst( status.eq.0, 'Error putting plasma_transport IDS !')
#endif
        call ids_put( idx, "radiation", radiation, status )
        call xertst( status.eq.0, 'Error putting radiation IDS !')
#if ( IMAS_MAJOR_VERSION < 4 || ( IMAS_MAJOR_VERSION == 4 && IMAS_MINOR_VERSION < 1 ) )
        call ids_put( idx, "dataset_description", description, status )
        call xertst( status.eq.0, 'Error putting dataset_description IDS !')
#endif
#if ( IMAS_MINOR_VERSION > 21 || IMAS_MAJOR_VERSION > 3 )
        call ids_put( idx, "summary", summary, status )
        call xertst( status.eq.0, 'Error putting summary IDS !')
#endif
#if ( IMAS_MINOR_VERSION > 25 && IMAS_MINOR_VERSION < 34 && IMAS_MAJOR_VERSION == 3 )
        call ids_put( idx, "numerics", numerics, status )
        call xertst( status.eq.0, 'Error putting numerics IDS !')
#endif
#if ( IMAS_MINOR_VERSION > 30 || IMAS_MAJOR_VERSION > 3 )
        call ids_put( idx, "divertors", divertors, status )
        call xertst( status.eq.0, 'Error putting divertors IDS !')
#endif

        write(*,*) "IDS rewrite finished"
        return

    end subroutine new_ids_edge

    subroutine new_batch_edge( idx, &
#if IMAS_MAJOR_VERSION > 3
            &   uri, &
#endif
            &   batch_profiles, batch_sources, &
#if IMAS_MAJOR_VERSION > 3
            &   batch_plasma_profiles, batch_plasma_sources, &
#endif
            &   equilibrium, &
#if ( IMAS_MINOR_VERSION > 21 || IMAS_MAJOR_VERSION > 3 )
            &   summary, &
#endif
#if ( IMAS_MAJOR_VERSION < 4 || ( IMAS_MAJOR_VERSION == 4 && IMAS_MINOR_VERSION < 1 ) )
            &   description, &
#endif
            &   do_summary, new_eq_ggd )
        type (ids_edge_profiles), intent(inout) :: batch_profiles   !< IDS
            !< designed to store data on edge plasma profiles (includes the
            !< scrape-off layer and possibly part of the confined plasma)
        type (ids_edge_sources), intent(inout) :: batch_sources     !< IDS
            !< designed to store data on edge plasma sources. Energy terms
            !< correspond to the full kinetic energy equation (i.e. the energy
            !< flux takes into account the energy transported by the particle
            !< flux)
#if IMAS_MAJOR_VERSION > 3
        type (ids_plasma_profiles), intent(inout) :: batch_plasma_profiles
            !< IDS designed to store data on plasma profiles
        type (ids_plasma_sources), intent(inout) :: batch_plasma_sources
            !< IDS designed to store data on plasma sources. Energy terms
            !< correspond to the full kinetic energy equation (i.e. the energy
            !< flux takes into account the energy transported by the particle
            !< flux)
#endif
#if ( IMAS_MAJOR_VERSION < 4 || ( IMAS_MAJOR_VERSION == 4 && IMAS_MINOR_VERSION < 1 ) )
        type (ids_dataset_description) :: description !< IDS
            !< designed to store a description of the simulation
#endif
        type (ids_equilibrium) :: equilibrium !< IDS
            !< designed to store a description of the equilibrium
#if ( IMAS_MINOR_VERSION > 21 || IMAS_MAJOR_VERSION > 3 )
        type (ids_summary), intent(inout) :: summary !< IDS
            !< designed to store run summary data
#endif
        integer, intent(inout) :: idx !< The returned identifier to be used in the
#if IMAS_MAJOR_VERSION > 3
        character(len=STRMAXLEN), intent(in) :: uri
#endif
        logical, intent(in) :: do_summary, new_eq_ggd
            !< subsequent data access operation
        integer :: status

            !< procedures
        external xertst

        !! Set data to edge_profiles IDS
        if (do_summary) then
          write(*,'(1x,a)') "Writing batch_profiles, batch_sources, "// &
#if IMAS_MAJOR_VERSION > 3
            &  "batch_plasma_profiles, batch_plasma_sources, "// &
#endif
#if ( IMAS_MAJOR_VERSION < 4 || ( IMAS_MAJOR_VERSION == 4 && IMAS_MINOR_VERSION < 1 ) )
#if ( IMAS_MINOR_VERSION > 21 || IMAS_MAJOR_VERSION > 3 )
            &  "summary, "// &
#endif
            &  "and dataset_description IDS"
#else
            &  "and summary IDS"
#endif
        else
          write(*,'(1x,a)') "Writing batch_profiles and batch_sources IDS "
        end if

        if (new_eq_ggd) then
          write(*,'(1x,a)') "Adding GGD data to equilibrium IDS"
          call ids_put( idx, "equilibrium", equilibrium, status )
          call xertst( status.eq.0, 'Error putting equilibrium IDS !')
        end if

        !! Put data to IDS
        call ids_put( idx, "edge_profiles/1", batch_profiles, status )
        call xertst( status.eq.0, 'Error putting batch_profiles IDS !')
        call ids_put( idx, "edge_sources/1", batch_sources, status )
        call xertst( status.eq.0, 'Error putting batch_sources IDS !')
#if IMAS_MAJOR_VERSION > 3
        call ids_put( idx, "plasma_profiles/1", batch_plasma_profiles, status )
        call xertst( status.eq.0, 'Error putting batch_plasma_profiles IDS !')
        call ids_put( idx, "plasma_sources/1", batch_plasma_sources, status )
        call xertst( status.eq.0, 'Error putting batch_plasma_sources IDS !')
#endif
        if (do_summary) then
#if ( IMAS_MAJOR_VERSION == 4 && IMAS_MINOR_VERSION == 0 )
          allocate( description%uri(1) )
          description%uri = trim(uri)
#elif ( IMAS_MAJOR_VERSION > 4 || ( IMAS_MAJOR_VERSION == 4 && IMAS_MINOR_VERSION > 0 ) )
          allocate( summary%identifier(1) )
          summary%identifier = trim(uri)
#endif
#if ( IMAS_MAJOR_VERSION < 4 || ( IMAS_MAJOR_VERSION == 4 && IMAS_MINOR_VERSION < 1 ) )
          call ids_put( idx, "dataset_description", description, status )
          call xertst( status.eq.0, 'Error putting dataset_description IDS !')
#endif
#if ( IMAS_MINOR_VERSION > 21 || IMAS_MAJOR_VERSION > 3 )
          call ids_put( idx, "summary", summary, status )
          call xertst( status.eq.0, 'Error putting summary IDS !')
#endif
        end if

        write(*,*) "IDS rewrite finished for averaged solution"
        return

    end subroutine new_batch_edge

    !> Example subroutine for reading edge_profiles IDS
    !! with Fortran90
    subroutine read_ids( idx, &
#if AL_MAJOR_VERSION > 4
         & ids_path )
#else
         & treename, shot, run, username, database, version )
#endif
        use ids_routines &  ! IGNORE
         & , only : imas_close
        implicit none
        integer, intent(out) :: idx !< The returned identifier to be used in the subsequent
#if AL_MAJOR_VERSION > 4
        character(len=256), intent(in) :: ids_path  !< The path to the IMAS data entry
#else
        character(len=24), intent(in) :: treename   !< The name of the IMAS IDS database
        integer, intent(in) :: shot !< The shot number of the database being created
        integer, intent(in) :: run  !< The run number of the database being created
        character(len=24), intent(in) :: username   !< Creator/owner of the IMAS IDS database
        character(len=24), intent(in) :: database   !< IMAS IDS database name
            !< (i. e. solps-iter, ITER, aug)
        character(len=24), intent(in) :: version    !< Major version of the IMAS IDS database
#endif
        !! Internal variables
#if AL_MAJOR_VERSION > 4
        character(len=:), allocatable :: message
        character(len=STRMAXLEN) :: uri
#endif
        integer :: gridSubset_index !< >Grid subset base index
        type(ids_edge_profiles) :: edge_profiles    !< IDS designed to store
            !< data in edge plasma profiles (includes the scrape-off layer and
            !<  possibly part of the confined plasma)
        integer :: status

        gridSubset_index = 3

        !! Open input datafile from local database
#if AL_MAJOR_VERSION > 4
        uri = 'imas:mdsplus?path='//trim(ids_path)
        write(0,*) "Started reading input IMAS data entry", trim(uri)
        call imas_open( uri, OPEN_PULSE, idx, status, message )
        call xertst ( status.eq.0, trim(message) )
#else
        write(0,*) "Started reading input IMAS data entry", idx, shot, run
        call imas_open_env(treename, shot, run, idx, username, &
            &   database, version, status )
        call xertst ( status.eq.0, 'Error opening IMAS database !')
#endif
        call ids_get(idx, "edge_profiles", edge_profiles, status)
        call xertst ( status.eq.0, 'Error opening edge_profiles IDS !')

        write(0,*) "homogeneous_time = ",   &
            &   edge_profiles%ids_properties%homogeneous_time
#if ( IMAS_MINOR_VERSION < 15 && IMAS_MAJOR_VERSION < 4 )
        write(0,*) "Grid subset 3 name = ", edge_profiles%ggd(1)%grid%  &
            &   grid_subset(gridSubset_index)%identifier%name
        write(0,*) "Grid subset 3 index = ", edge_profiles%ggd(1)%grid% &
            &   grid_subset(gridSubset_index)%identifier%index
#else
        write(0,*) "Grid subset 3 name = ", edge_profiles%grid_ggd(1)%  &
            &   grid_subset(gridSubset_index)%identifier%name
        write(0,*) "Grid subset 3 index = ", edge_profiles%grid_ggd(1)% &
            &   grid_subset(gridSubset_index)%identifier%index
#endif
        ! write(0,*) "Time = ", edge_profiles%time(1)
        call ids_deallocate( edge_profiles )
        call imas_close( idx, status )
        call xertst ( status.eq.0, 'Error closing IMAS database !')
        write(0,*) "Finished reading input IMAS data entry"

    end subroutine read_ids

#endif

    !> Routine to open UAL database.
    !! @note For IMAS IDS is recommended use of IMAS GGD library routine
    !! "exampleOpenIDS"
    subroutine open_ual( idx, shot, run, time, user, database, dataversion,  &
        &   doCreate, useHdf5, nmlFile )
        integer, intent(out) :: idx !< The returned identifier to be used in the
                                    !< subsequent data access operation
        integer, intent(in), optional :: shot   !< The pulse (previously shot) number of the
                                                !< database being created
        integer, intent(in), optional :: run    !< The run number of the
                                                !< database being created
        real(R8), intent(out), optional :: time !< Generic time.
                                                !! Time is special: it is not
                                                !! used here, but can be read
                                                !! from the namelist and
                                                !! returned
        character(*), intent(in), optional :: user  !< Creator/owner of the
                                                    !< database
        character(*), intent(in), optional :: database !< Database name
                                                       !< (i. e. solps-iter, ITER, aug)
        character(*), intent(in), optional :: dataversion   !< Major version of
                                                            !< the database
        logical, intent(in), optional :: doCreate
        logical, intent(in), optional :: useHdf5
        character(*), intent(in), optional :: nmlFile

        !! Internal variables

        character(*), parameter :: NAMELIST_FILE = "ual.namelist"
        integer, parameter :: NAMELIST_UNIT = 979

        integer :: lRefshot = 0, lRefrun = 0
#ifdef IMAS
        integer :: lStatus = 0
        character(32) :: lTreename = "ids"
# if ( AL_MAJOR_VERSION < 4 || ( AL_MAJOR_VERSION == 4 && AL_MINOR_VERSION < 9 ) )
        character(13) :: hlp_frm
        character(80) :: message
        integer len_of_digits
        external len_of_digits
# elif AL_MAJOR_VERSION > 4
        character(len=:), allocatable :: message
        character(len=STRMAXLEN) :: uri
# endif
#elif defined(ITM_ENVIRONMENT_LOADED)
        character(32) :: lTreename = "euitm"
#else
        character(32) :: lTreename = "none"
#endif

        integer :: lShot = 1, lRun = 0
        real(R8) :: lTime = 0.0_R8
        character(32) :: lUser="unspecified", lTokamak="unspecified",   &
            &   lDataversion="unspecified"
        logical :: lDoCreate = .false., lUseHdf5 = .false.

        logical :: namelistExists, openEnv = .false.

        external xerrab, xertst

        namelist /ual_namelist/ lTreename, lShot, lRun, lTime, lRefshot,    &
            &   lRefrun, lUser, lTokamak, lDataversion, openEnv, lDoCreate, &
            &   lUseHdf5

        if( present( shot ) ) lShot = shot
        if( present( run ) ) lRun = run
        if( present( user ) ) lUser = user
        if( present( database ) ) lTokamak = database
        if( present( dataversion ) ) lDataversion = dataversion
        if( present( doCreate ) ) lDoCreate = doCreate
        if( present( useHdf5 ) ) lUseHdf5 = useHdf5

        if( present( user ) ) openEnv = .true.

        !! If file exists, read namelist from configuration file
        !! If not, write out namelist
        if( present( nmlFile ) ) then
            inquire( file=nmlFile, exist=namelistExists )
        else
            inquire( file=NAMELIST_FILE, exist=namelistExists )
        end if
        if( namelistExists ) then
            if( present( nmlFile ) ) then
                open( unit=NAMELIST_UNIT, file=nmlFile )
            else
                open( unit=NAMELIST_UNIT, file=NAMELIST_FILE )
            end if
            read( NAMELIST_UNIT, nml=ual_namelist )
            close( unit=NAMELIST_UNIT )
        else
            if( present( nmlFile ) ) then
                open( unit=NAMELIST_UNIT, file=nmlFile, status="new", &
                    &   action="write" )
            else
                open( unit=NAMELIST_UNIT, file=NAMELIST_FILE, status="new", &
                    &   action="write" )
            end if
            write( NAMELIST_UNIT, nml=ual_namelist )
            close( unit=NAMELIST_UNIT )
        end if

        !! establish UAL access
        if( lDoCreate ) then
#ifdef IMAS
            if( lUseHdf5 ) then
# if AL_MAJOR_VERSION > 4
                call al_build_uri_from_legacy_parameters &
                  & ( HDF5_BACKEND, lShot, lRun, lUser, lTokamak, lDataversion, &
                  &   '', uri, lStatus )
                call imas_open( uri, FORCE_CREATE_PULSE, idx, lStatus, message )
                call xertst ( lStatus.eq.0, trim(message) )
# else
#  if ( AL_MAJOR_VERSION == 4 && AL_MINOR_VERSION > 8 )
                call ual_begin_pulse_action( HDF5_BACKEND, lShot, lRun, lUser, &
                   &   lTokamak, lDataversion, idx )
                call ual_open_pulse( idx, FORCE_CREATE_PULSE, '', lStatus )
                call xertst ( lStatus.eq.0, 'Error opening IMAS database !')
#  else
                write(hlp_frm,'(a,i1,a)') &
                   &  '(a,i1,a,i',len_of_digits(AL_MINOR_VERSION),',a)'
                write(message,hlp_frm) &
                   &  'HDF5 backend not supported with AL v', &
                   &   AL_MAJOR_VERSION,'.',AL_MINOR_VERSION,'!'
                call xerrab ( trim(message) )
#  endif
# endif
            else
                if( openEnv ) then
# if AL_MAJOR_VERSION > 4
                    call al_build_uri_from_legacy_parameters &
                        & ( MDSPLUS_BACKEND, lShot, lRun, lUser, lTokamak, lDataversion, &
                        &   '', uri, lStatus )
                    call imas_open( uri, FORCE_CREATE_PULSE, idx, lStatus, message )
                    call xertst ( lStatus.eq.0, trim(message) )
# else
                    call imas_create_env( lTreename, lShot, lRun, lRefshot, &
                        &   lRefrun, idx, lUser, lTokamak, lDataversion,    &
                        &   lStatus)
                    call xertst ( lStatus.eq.0, 'Error opening IMAS database !')
# endif
                else
# if AL_MAJOR_VERSION < 4
                    call imas_create( lTreename, lShot, lRun, lRefshot, &
                        &   lRefrun, idx)
# else
                    call xerrab ('Must define username!')
# endif
                end if
            end if
        else
            if( lUseHdf5 ) then
# if AL_MAJOR_VERSION > 4
                call al_build_uri_from_legacy_parameters &
                  & ( HDF5_BACKEND, lShot, lRun, lUser, lTokamak, lDataversion, &
                  &   '', uri, lStatus )
                call imas_open ( uri, OPEN_PULSE, idx, lStatus, message )
                call xertst ( lStatus.eq.0, trim(message) )
# else
#  if ( AL_MAJOR_VERSION == 4 && AL_MINOR_VERSION > 8 )
                call ual_begin_pulse_action( HDF5_BACKEND, lShot, lRun, lUser, &
                        &    lTokamak, lDataversion, idx )
                call ual_open_pulse( idx, OPEN_PULSE, '', lStatus )
                call xertst ( lStatus.eq.0, 'Error opening IMAS data entry !')
#  else
                write(hlp_frm,'(a,i1,a)') &
                   &  '(a,i1,a,i',len_of_digits(AL_MINOR_VERSION),',a)'
                write(message,hlp_frm) &
                   &  'HDF5 backend not supported with AL v', &
                   &   AL_MAJOR_VERSION,'.',AL_MINOR_VERSION,'!'
                call xerrab ( trim(message) )
#  endif
# endif
            else
                if( openEnv ) then
# if AL_MAJOR_VERSION > 4
                    call al_build_uri_from_legacy_parameters &
                        & ( MDSPLUS_BACKEND, lShot, lRun, lUser, lTokamak, lDataversion, &
                        &   '', uri, lStatus )
                    call imas_open( uri, OPEN_PULSE, idx, lStatus, message )
                    call xertst ( lStatus.eq.0, trim(message) )
# else
                    call imas_open_env(lTreename, lShot, lRun, idx, lUser, &
                        &   lTokamak, lDataversion, lStatus)
                    call xertst ( lStatus.eq.0, 'Error opening IMAS data entry !')
# endif
                else
# if AL_MAJOR_VERSION < 4
                    call imas_open(lTreename, lShot, lRun, idx)
# else
                    call xerrab ('Must define username!')
# endif
                end if
            end if
#elif defined(ITM_ENVIRONMENT_LOADED)
            if( lUseHdf5 ) then
                call euITM_create_hdf5(lTreename, lShot, lRun, lRefshot, &
                        &   lRefrun, idx)
            else
                if( openEnv ) then
                    call euITM_create_env(lTreename, lShot, lRun, lRefshot, &
                        &   lRefrun, idx, lUser, lTokamak, lDataversion)
                else
                    call euITM_create(lTreename, lShot, lRun, lRefshot, &
                        &   lRefrun, idx)
                end if
            end if
        else
            if( lUseHdf5 ) then
                call euITM_open_hdf5(lTreename, lShot, lRun, idx)
            else
                if( openEnv) then
                    call euITM_open_env(lTreename, lShot, lRun, idx,  lUser,  &
                        &   lTokamak, lDataversion)
                else
                    call euITM_open(lTreename, lShot, lRun, lRefshot, &
                        &   lRefrun, idx)
                end if
            end if
#else
            idx = 0
#endif
        end if

        !! Return time if requested
        if( present( time ) ) time = lTime

    end subroutine open_ual

    !> Close UAL.
    !! @note  For IMAS edge_profiles IDS "examplePutIDS" IMAS GGD library routine
    !! can be used instead (that routine also writes the set data to IDS and then
    !! closes the IDS)
    subroutine close_ual(idx)
        integer, intent(inout) :: idx  !< The returned identifier to be used in the
                                       !< subsequent data access operation
#ifdef IMAS
        integer :: status
        external xertst

        !! Close IDS
# if AL_MAJOR_VERSION > 4
        call imas_close( idx, status )
# else
        call ual_close_pulse( idx, CLOSE_PULSE, '', status )
# endif
        call xertst( status.eq.0, 'Error closing IMAS database !' )
#elif defined(ITM_ENVIRONMENT_LOADED)
        call euITM_close( idx )
#endif
        idx = 0

    end subroutine close_ual

end module b2mod_ual

!!!Local Variables:
!!! mode: f90
!!! End:
