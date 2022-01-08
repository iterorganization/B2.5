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
     & ,only: imas_open_env, imas_create_env, imas_close, &
     &        ids_deallocate, ids_get, ids_put, ids_delete, ids_put_slice
    use ids_schemas &   ! IGNORE
     & ,only: ids_edge_profiles, ids_edge_sources, ids_edge_transport, &
     &        ids_radiation, ids_dataset_description, ids_equilibrium
    use b2mod_ual_io &
     & ,only: b25_process_ids
#if IMAS_MINOR_VERSION > 21
    use ids_schemas &   ! IGNORE
     & ,only: ids_summary
#endif
#if ( IMAS_MINOR_VERSION > 25 && IMAS_MINOR_VERSION < 34 )
    use ids_schemas &   ! IGNORE
     & ,only: ids_numerics
#endif
#if IMAS_MINOR_VERSION > 30
    use ids_schemas &   ! IGNORE
     & ,only: ids_divertors
#endif
#else
# ifdef ITM_ENVIRONMENT_LOADED
    use euITM_schemas   ! IGNORE
    use euITM_routines  ! IGNORE
# endif
#endif

  implicit none

  private

  public open_ual, close_ual
#ifdef IMAS
  public put_ids_edge, new_ids_edge, delete_ids_edge
  public dealloc_ids_edge, dealloc_batch_edge
  public put_batch_edge, new_batch_edge
  public b25_process_ids
  public ids_edge_profiles, ids_edge_sources, ids_edge_transport, &
    &    ids_radiation, ids_dataset_description, ids_equilibrium
#if IMAS_MINOR_VERSION > 21
  public ids_summary
#endif
#if ( IMAS_MINOR_VERSION > 25 && IMAS_MINOR_VERSION < 34 )
  public ids_numerics
#endif
#if IMAS_MINOR_VERSION > 30
  public ids_divertors
#endif
#endif


contains

#ifdef IMAS
    !> Subroutine used to put data to edge_profiles, edge_sources and
    !! edge_transport IDSs.
    subroutine put_ids_edge( edge_profiles, edge_sources, edge_transport, &
            &   radiation, description, &
#if IMAS_MINOR_VERSION > 21
            &   summary, &
#endif
#if ( IMAS_MINOR_VERSION > 25 && IMAS_MINOR_VERSION < 34 )
            &   numerics, &
#endif
#if IMAS_MINOR_VERSION > 30
            &   divertors, &
#endif
            &   treename, shot, run, idx, username, database, version )
        type (ids_edge_profiles), intent(inout) :: edge_profiles   !< IDS
            !< designed to store data on edge plasma profiles (includes the
            !< scrape-off layer and possibly part of the confined plasma)
        type (ids_edge_sources), intent(inout) :: edge_sources     !< IDS
            !< designed to store data on edge plasma sources. Energy terms
            !< correspond to the full kinetic energy equation (i.e. the energy
            !< flux takes into account the energy transported by the particle
            !< flux)
        type (ids_edge_transport), intent(inout) :: edge_transport !< IDS
            !< designed to store  data on edge plasma transport. Energy terms
            !< correspond to the full kinetic energy equation (i.e. the energy
            !< flux takes into account the energy transported by the particle
            !< flux)
        type (ids_radiation), intent(inout) :: radiation !< IDS
            !< designed to store data about plasma radiation
        type (ids_dataset_description), intent(inout) :: description !< IDS
            !<  designed to store a description of the simulation
#if IMAS_MINOR_VERSION > 21
        type (ids_summary), intent(inout) :: summary !< IDS
            !< designed to store run summary data
#endif
#if ( IMAS_MINOR_VERSION > 25 && IMAS_MINOR_VERSION < 34 )
        type (ids_numerics), intent(inout) :: numerics !< IDS designed to store
            !< run numerics data
#endif
#if IMAS_MINOR_VERSION > 30
        type (ids_divertors), intent(inout) :: divertors !< IDS
            !< designed to store run data related to the divertor plates
#endif
        character(len=24), intent(in) :: treename   !< The name of the IMAS IDS database
            !< (i.e. "edge_profiles" (mandatory) )
        integer, intent(in) :: shot   !< The shot number of the database being created
        integer, intent(in) :: run    !< The run number of the database being created
        integer, intent(inout) :: idx !< The returned identifier to be used in the
            !< subsequent data access operation
        character(len=24), intent(in) :: username   !< Creator/owner of the IMAS IDS
            !< database
        character(len=24), intent(in) :: database   !< IMAS database name
            !< (i. e. solps-iter, ITER, aug)
        character(len=24), intent(in) :: version    !< Major version of the IMAS IDS
            !< database
        integer :: status

            !< procedures
        external xertst, xerrab

        !! Set data to edge_profiles IDS
        write(*,'(1x,a)') "Writing edge_profiles, edge_sources, edge_transport, "// &
#if IMAS_MINOR_VERSION > 21
          &  "summary, "// &
#endif
#if ( IMAS_MINOR_VERSION > 25 && IMAS_MINOR_VERSION < 34 )
          &  "numerics, "// &
#endif
#if IMAS_MINOR_VERSION > 30
          &  "divertors, "// &
#endif
          &  "dataset_description, and radiation IDS"

        !! Create and modify new shot/run
        if ( idx.eq.0 ) then
          call imas_create_env( treename, shot, run, 0, 0, idx, username, &
             & database, version, status )
          call xertst( status.eq.0, 'Error opening IMAS database !')

        !! Put data to IDS
          call ids_put( idx, "edge_profiles", edge_profiles, status )
          call xertst( status.eq.0, 'Error putting edge_profiles IDS !')
          call ids_put( idx, "edge_sources", edge_sources, status )
          call xertst( status.eq.0, 'Error putting edge_sources IDS !')
          call ids_put( idx, "edge_transport", edge_transport, status )
          call xertst( status.eq.0, 'Error putting edge_transport IDS !')
          call ids_put( idx, "radiation", radiation, status )
          call xertst( status.eq.0, 'Error putting radiation IDS !')
          call ids_put( idx, "dataset_description", description, status )
          call xertst( status.eq.0, 'Error putting dataset_description IDS !')
#if IMAS_MINOR_VERSION > 21
          call ids_put( idx, "summary", summary, status )
          call xertst( status.eq.0, 'Error putting summary IDS !')
#endif
#if ( IMAS_MINOR_VERSION > 25 && IMAS_MINOR_VERSION < 34 )
          call ids_put( idx, "numerics", numerics, status )
          call xertst( status.eq.0, 'Error putting numerics IDS !')
#endif
#if IMAS_MINOR_VERSION > 30
          call ids_put( idx, "divertors", divertors, status )
          call xertst( status.eq.0, 'Error putting divertors IDS !')
#endif
        else
        !! Or open and modify existing shot/run
        !! (might work much faster than imas_create_env)
        ! call imas_open_env(treename, shot, run, idx, username, &
        !  database, version, status )

        !! Put data to IDS
          call ids_put_slice( idx, "edge_profiles", edge_profiles, status )
          call xertst( status.eq.0, 'Error putting slice in edge_profiles IDS !')
          call ids_put_slice( idx, "edge_sources", edge_sources, status )
          call xertst( status.eq.0, 'Error putting slice in edge_sources IDS !')
          call ids_put_slice( idx, "edge_transport", edge_transport, status )
          call xertst( status.eq.0, 'Error putting slice in edge_transport IDS !')
          call ids_put_slice( idx, "radiation", radiation, status )
          call xertst( status.eq.0, 'Error putting slice in radiation IDS !')
          call ids_put_slice( idx, "dataset_description", description, status )
          call xertst( status.eq.0, 'Error putting slice in dataset_description IDS !')
#if IMAS_MINOR_VERSION > 21
          call ids_put_slice( idx, "summary", summary, status )
          call xertst( status.eq.0, 'Error putting slice in summary IDS !')
#endif
#if ( IMAS_MINOR_VERSION > 25 && IMAS_MINOR_VERSION < 34 )
          call ids_put_slice( idx, "numerics", numerics, status )
          call xertst( status.eq.0, 'Error putting slice in numerics IDS !')
#endif
#if IMAS_MINOR_VERSION > 30
          call ids_put_slice( idx, "divertors", divertors, status )
          call xertst( status.eq.0, 'Error putting slice in divertors IDS !')
#endif
        end if

        write(*,*) "IDS write finished"
        return

    end subroutine put_ids_edge

    subroutine dealloc_ids_edge( edge_profiles, edge_sources, edge_transport, &
#if ( IMAS_MINOR_VERSION > 25 && IMAS_MINOR_VERSION < 34 )
            &   numerics, &
#endif
#if IMAS_MINOR_VERSION > 30
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
            !< designed to store  data on edge plasma transport. Energy terms
            !< correspond to the full kinetic energy equation (i.e. the energy
            !< flux takes into account the energy transported by the particle
            !< flux)
        type (ids_radiation), intent(inout) :: radiation !< IDS
            !< designed to store data about plasma radiation
#if ( IMAS_MINOR_VERSION > 25 && IMAS_MINOR_VERSION < 34 )
        type (ids_numerics), intent(inout) :: numerics !< IDS designed to store
            !< run numerics data
#endif
#if IMAS_MINOR_VERSION > 30
        type (ids_divertors), intent(inout) :: divertors !< IDS
            !< designed to store run data related to the divertor plates
#endif
        if (.not.associated( edge_profiles%ids_properties%comment )) return
        call ids_deallocate( edge_profiles )
        call ids_deallocate( edge_sources )
        call ids_deallocate( edge_transport )
        call ids_deallocate( radiation )
#if ( IMAS_MINOR_VERSION > 25 && IMAS_MINOR_VERSION < 34 )
        call ids_deallocate( numerics )
#endif
#if IMAS_MINOR_VERSION > 30
        call ids_deallocate( divertors )
#endif
         return
    end subroutine dealloc_ids_edge

    subroutine put_batch_edge( &
            &   treename, shot, run, idx, username, database, version, &
            &   batch_profiles, batch_sources, &
#if IMAS_MINOR_VERSION > 21
            &   summary, &
#endif
            &   description, do_summary )
        type (ids_edge_profiles), intent(inout) :: batch_profiles   !< IDS
            !< designed to store data on edge plasma profiles (includes the
            !< scrape-off layer and possibly part of the confined plasma)
        type (ids_edge_sources), intent(inout) :: batch_sources     !< IDS
            !< designed to store data on edge plasma sources. Energy terms
            !< correspond to the full kinetic energy equation (i.e. the energy
            !< flux takes into account the energy transported by the particle
            !< flux)
        type (ids_dataset_description), intent(inout) :: description
            !< IDS designed to store a description of the simulation
#if IMAS_MINOR_VERSION > 21
        type (ids_summary), intent(inout) :: summary !< IDS
            !< designed to store run summary data
#endif
        character(len=24), intent(in) :: treename   !< The name of the IMAS IDS database
            !< (i.e. "edge_profiles" (mandatory) )
        integer, intent(in) :: shot   !< The shot number of the database being created
        integer, intent(in) :: run    !< The run number of the database being created
        integer, intent(inout) :: idx !< The returned identifier to be used in the
            !< subsequent data access operation
        character(len=24), intent(in) :: username   !< Creator/owner of the IMAS IDS
            !< database
        character(len=24), intent(in) :: database   !< IMAS database name
            !< (i. e. solps-iter, ITER, aug)
        character(len=24), intent(in) :: version    !< Major version of the IMAS IDS
        logical, intent(in) :: do_summary
            !< database
        integer :: status

            !< procedures
        external xertst, xerrab

        !! Set data to edge_profiles IDS
        if (do_summary) then
          write(*,'(1x,a)') "Writing batch_profiles, batch_sources, "// &
#if IMAS_MINOR_VERSION > 21
            &  "summary, "// &
#endif
            &  "and dataset_description IDS"
        else
          write(*,'(1x,a)') "Writing batch_profiles and batch_sources IDS "
        end if

        !! Create and modify new shot/run
        if ( idx.eq.0 ) then
          call imas_create_env( treename, shot, run, 0, 0, idx, username, &
             & database, version, status )
          call xertst( status.eq.0, 'Error opening IMAS database !')

        !! Put data to IDS
          call ids_put( idx, "edge_profiles/1", batch_profiles, status )
          call xertst( status.eq.0, 'Error putting batch_profiles IDS !')
          call ids_put( idx, "edge_sources/1", batch_sources, status )
          call xertst( status.eq.0, 'Error putting batch_sources IDS !')
          if (do_summary) then
            call ids_put( idx, "dataset_description", description, status )
            call xertst( status.eq.0, 'Error putting dataset_description IDS !')
#if IMAS_MINOR_VERSION > 21
            call ids_put( idx, "summary", summary, status )
            call xertst( status.eq.0, 'Error putting summary IDS !')
#endif
          end if
        else
        !! Or open and modify existing shot/run
        !! (might work much faster than imas_create_env)

        !! Put data to IDS
          call ids_put_slice( idx, "edge_profiles/1", batch_profiles, status )
          call xertst( status.eq.0, 'Error putting slice in batch_profiles IDS !')
          call ids_put_slice( idx, "edge_sources/1", batch_sources, status )
          call xertst( status.eq.0, 'Error putting slice in batch_sources IDS !')
          if (do_summary) then
            call ids_put_slice( idx, "dataset_description", description, status )
            call xertst( status.eq.0, 'Error putting slice in dataset_description IDS !')
#if IMAS_MINOR_VERSION > 21
            call ids_put_slice( idx, "summary", summary, status )
            call xertst( status.eq.0, 'Error putting slice in summary IDS !')
#endif
          end if
        end if

        write(*,*) "IDS write finished for batch averages"
        return

    end subroutine put_batch_edge

    subroutine dealloc_batch_edge( batch_profiles, batch_sources, &
#if IMAS_MINOR_VERSION > 21
            &   summary, &
#endif
            &   description )
        implicit none
        type (ids_edge_profiles), intent(inout) :: batch_profiles   !< IDS
            !< designed to store data on edge plasma profiles (includes the
            !< scrape-off layer and possibly part of the confined plasma)
        type (ids_edge_sources), intent(inout) :: batch_sources     !< IDS
            !< designed to store data on edge plasma sources. Energy terms
            !< correspond to the full kinetic energy equation (i.e. the energy
            !< flux takes into account the energy transported by the particle
            !< flux)
        type (ids_dataset_description), intent(inout) :: description
            !< IDS designed to store a description of the simulation
#if IMAS_MINOR_VERSION > 21
        type (ids_summary), intent(inout) :: summary !< IDS
            !< designed to store run summary data
#endif
        if (associated( batch_profiles%ids_properties%comment ) ) then
          call ids_deallocate( batch_profiles )
          call ids_deallocate( batch_sources )
        end if
        if (associated( description%ids_properties%comment ) ) then
          call ids_deallocate( description )
#if IMAS_MINOR_VERSION > 21
          call ids_deallocate( summary )
#endif
        end if
        return
    end subroutine dealloc_batch_edge

    !> Subroutine used to delete data from edge_profiles, edge_sources and
    !! edge_transport IDSs.
    subroutine delete_ids_edge( edge_profiles, edge_sources, edge_transport, &
            &   radiation, description, &
#if IMAS_MINOR_VERSION > 21
            &   summary, &
#endif
#if ( IMAS_MINOR_VERSION > 25 && IMAS_MINOR_VERSION < 34 )
            &   numerics, &
#endif
#if IMAS_MINOR_VERSION > 30
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
            !< designed to store  data on edge plasma transport. Energy terms
            !< correspond to the full kinetic energy equation (i.e. the energy
            !< flux takes into account the energy transported by the particle
            !< flux)
        type (ids_radiation), intent(inout) :: radiation !< IDS
            !< designed to store data about plasma radiation
        type (ids_dataset_description) :: description !< IDS designed to store
            !< a description of the simulation
#if IMAS_MINOR_VERSION > 21
        type (ids_summary), intent(inout) :: summary !< IDS
            !< designed to store run summary data
#endif
#if ( IMAS_MINOR_VERSION > 25 && IMAS_MINOR_VERSION < 34 )
        type (ids_numerics), intent(inout) :: numerics !< IDS designed to store
            !< run numerics data
#endif
#if IMAS_MINOR_VERSION > 30
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
          call ids_delete( idx, "radiation", radiation)
          call ids_delete( idx, "dataset_description", description)
#if IMAS_MINOR_VERSION > 21
          call ids_delete( idx, "summary", summary)
#endif
#if ( IMAS_MINOR_VERSION > 25 && IMAS_MINOR_VERSION < 34 )
          call ids_delete( idx, "numerics", numerics)
#endif
#if IMAS_MINOR_VERSION > 30
          call ids_delete( idx, "divertors", divertors)
#endif
          write(*,*) "IDS delete finished"
        end if
        return

    end subroutine delete_ids_edge

    !> Subroutine used to rewrite data to edge_profiles, edge_sources and
    !! edge_transport IDSs.
    subroutine new_ids_edge( edge_profiles, edge_sources, edge_transport, &
            &   radiation, description, &
#if IMAS_MINOR_VERSION > 21
            &   summary, &
#endif
#if ( IMAS_MINOR_VERSION > 25 && IMAS_MINOR_VERSION < 34 )
            &   numerics, &
#endif
#if IMAS_MINOR_VERSION > 30
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
            !< designed to store  data on edge plasma transport. Energy terms
            !< correspond to the full kinetic energy equation (i.e. the energy
            !< flux takes into account the energy transported by the particle
            !< flux)
        type (ids_radiation), intent(inout) :: radiation !< IDS
            !< designed to store data about plasma radiation
        type (ids_dataset_description) :: description !< IDS designed to store
            !< a description of the simulation
#if IMAS_MINOR_VERSION > 21
        type (ids_summary), intent(inout) :: summary !< IDS
            !< designed to store run summary data
#endif
#if ( IMAS_MINOR_VERSION > 25 && IMAS_MINOR_VERSION < 34 )
        type (ids_numerics), intent(inout) :: numerics !< IDS designed to store
            !< run numerics data
#endif
#if IMAS_MINOR_VERSION > 30
        type (ids_divertors), intent(inout) :: divertors !< IDS
            !< designed to store data related to the divertor plates
#endif
        integer, intent(inout) :: idx !< The returned identifier to be used in the
            !< subsequent data access operation
        integer :: status

            !< procedures
        external xertst

        !! Set data to edge_profiles IDS
        write(*,'(1x,a)') "Writing edge_profiles, edge_sources, edge_transport, "// &
#if IMAS_MINOR_VERSION > 21
          &  "summary, "// &
#endif
#if ( IMAS_MINOR_VERSION > 25 && IMAS_MINOR_VERSION < 34 )
          &  "numerics, "// &
#endif
#if IMAS_MINOR_VERSION > 30
          &  "divertors, "// &
#endif
          &  "dataset_description, and radiation IDS"

        !! Put data to IDS
        call ids_put( idx, "edge_profiles", edge_profiles, status )
        call xertst( status.eq.0, 'Error putting edge_profiles IDS !')
        call ids_put( idx, "edge_sources", edge_sources, status )
        call xertst( status.eq.0, 'Error putting edge_sources IDS !')
        call ids_put( idx, "edge_transport", edge_transport, status )
        call xertst( status.eq.0, 'Error putting edge_transport IDS !')
        call ids_put( idx, "radiation", radiation, status )
        call xertst( status.eq.0, 'Error putting radiation IDS !')
        call ids_put( idx, "dataset_description", description, status )
        call xertst( status.eq.0, 'Error putting dataset_description IDS !')
#if IMAS_MINOR_VERSION > 21
        call ids_put( idx, "summary", summary, status )
        call xertst( status.eq.0, 'Error putting summary IDS !')
#endif
#if ( IMAS_MINOR_VERSION > 25 && IMAS_MINOR_VERSION < 34 )
        call ids_put( idx, "numerics", numerics, status )
        call xertst( status.eq.0, 'Error putting numerics IDS !')
#endif
#if IMAS_MINOR_VERSION > 30
        call ids_put( idx, "divertors", divertors, status )
        call xertst( status.eq.0, 'Error putting divertors IDS !')
#endif

        write(*,*) "IDS rewrite finished"
        return

    end subroutine new_ids_edge

    subroutine new_batch_edge( idx, batch_profiles, batch_sources, &
#if IMAS_MINOR_VERSION > 21
            &   summary, &
#endif
            &   description, do_summary )
        type (ids_edge_profiles), intent(inout) :: batch_profiles   !< IDS
            !< designed to store data on edge plasma profiles (includes the
            !< scrape-off layer and possibly part of the confined plasma)
        type (ids_edge_sources), intent(inout) :: batch_sources     !< IDS
            !< designed to store data on edge plasma sources. Energy terms
            !< correspond to the full kinetic energy equation (i.e. the energy
            !< flux takes into account the energy transported by the particle
            !< flux)
        type (ids_dataset_description) :: description !< IDS
            !< designed to store a description of the simulation
#if IMAS_MINOR_VERSION > 21
        type (ids_summary), intent(inout) :: summary !< IDS
            !< designed to store run summary data
#endif
        integer, intent(inout) :: idx !< The returned identifier to be used in the
        logical, intent(in) :: do_summary
            !< subsequent data access operation
        integer :: status

            !< procedures
        external xertst

        !! Set data to edge_profiles IDS
        if (do_summary) then
          write(*,'(1x,a)') "Writing batch_profiles, batch_sources, "// &
#if IMAS_MINOR_VERSION > 21
            &  "summary, "// &
#endif
            &  "and dataset_description IDS"
        else
          write(*,'(1x,a)') "Writing batch_profiles and batch_sources IDS "
        end if

        !! Put data to IDS
        call ids_put( idx, "edge_profiles/1", batch_profiles, status )
        call xertst( status.eq.0, 'Error putting batch_profiles IDS !')
        call ids_put( idx, "edge_sources/1", batch_sources, status )
        call xertst( status.eq.0, 'Error putting batch_sources IDS !')
        if (do_summary) then
          call ids_put( idx, "dataset_description", description, status )
          call xertst( status.eq.0, 'Error putting dataset_description IDS !')
#if IMAS_MINOR_VERSION > 21
          call ids_put( idx, "summary", summary, status )
          call xertst( status.eq.0, 'Error putting summary IDS !')
#endif
        end if

        write(*,*) "IDS rewrite finished for averaged solution"
        return

    end subroutine new_batch_edge

#endif

    !> Routine to open UAL database.
    !! @note For IMAS IDS is recommended use of IMAS GGD library routine
    !! "exampleOpenIDS"
    subroutine open_ual( idx, shot, run, time, user, database, dataversion,  &
        &   doCreate, useHdf5, nmlFile )
        integer, intent(out) :: idx !< The returned identifier to be used in the
                                    !< subsequent data access operation
        integer, intent(in), optional :: shot   !< The shot number of the
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
#else
# ifdef ITM_ENVIRONMENT_LOADED
        character(32) :: lTreename = "euitm"
# else
        character(32) :: lTreename = "none"
# endif
#endif

        integer :: lShot = 1, lRun = 0
        real(R8) :: lTime = 0.0_R8
        character(32) :: luser="unspecified", lTokamak="unspecified",   &
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
        if( lDoCreate) then
#ifdef IMAS
            if( lUseHdf5) then
# if UAL_MAJOR_VERSION < 4
                call imas_create_hdf5(lTreename, lShot, lRun, lRefshot, &
                        &   lRefrun, idx)
# else
                call xerrab ('HDF5 IMAS format not supported with UAL v4!')
# endif
            else
                if( openEnv) then
                    call imas_create_env(lTreename, lShot, lRun, lRefshot,  &
                        &   lRefrun, idx, lUser, lTokamak, lDataversion,    &
                        &   lStatus)
                    call xertst ( lStatus.eq.0, 'Error opening IMAS database !')
                else
# if UAL_MAJOR_VERSION < 4
                    call imas_create(lTreename, lShot, lRun, lRefshot, &
                        &   lRefrun, idx)
# else
                    call xerrab ('Must define username!')
# endif
                end if
            end if
        else
            if( lUseHdf5) then
# if UAL_MAJOR_VERSION < 4
                call imas_open_hdf5(lTreename, lShot, lRun, idx)
# else
                call xerrab ('HDF5 IMAS format not supported with UAL v4!')
# endif
            else
                if( openEnv) then
                    call imas_open_env(lTreename, lShot, lRun, idx, lUser, &
                        &   lTokamak, lDataversion, lStatus)
                    call xertst ( lStatus.eq.0, 'Error opening IMAS database !')
                else
# if UAL_MAJOR_VERSION < 4
                    call imas_open(lTreename, lShot, lRun, lRefshot, &
                        &   lRefrun, idx)
# else
                    call xerrab ('Must define username!')
# endif
                end if
            end if
#else
# ifdef ITM_ENVIRONMENT_LOADED
            if( lUseHdf5) then
                call euITM_create_hdf5(lTreename, lShot, lRun, lRefshot, &
                        &   lRefrun, idx)
            else
                if( openEnv) then
                    call euITM_create_env(lTreename, lShot, lRun, lRefshot, &
                        &   lRefrun, idx, lUser, lTokamak, lDataversion)
                else
                    call euITM_create(lTreename, lShot, lRun, lRefshot, &
                        &   lRefrun, idx)
                end if
            end if
        else
            if( lUseHdf5) then
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
# else
            idx = 0
# endif
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
        integer, intent(in) :: idx  !< The returned identifier to be used in the
                                    !< subsequent data access operation
#ifdef IMAS
        integer :: status
        external xertst

        !! Close IDS
        call imas_close( idx, status )
        call xertst( status.eq.0, 'Error closing IMAS database !' )
#else
# ifdef ITM_ENVIRONMENT_LOADED
        call euITM_close( idx )
# endif
#endif
    end subroutine close_ual

end module b2mod_ual

!!!Local Variables:
!!! mode: f90
!!! End:
