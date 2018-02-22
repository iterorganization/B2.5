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
!!      geometry ( full geometry descriptions of all available grid subsets )
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
!!-----------------------------------------------------------------------------

program b2_ual_write

    use b2mod_main
    use b2mod_ual
    use b2mod_grid_mapping
    use b2mod_ual_io
    use ids_schemas     ! IGNORE
                        !! These are the Fortran type definitions for the
                        !! Physics Data Model
    use ids_routines    ! IGNORE
                        !! These are the Access Layer routines + management of
                        !! IDS structures
    use ids_grid_common &       ! IGNORE
        & , IDS_COORDTYPE_R => COORDTYPE_R    &
        & , IDS_COORDTYPE_Z => COORDTYPE_Z
        ! &   GRID_UNDEFINED  => B2_GRID_UNDEFINED
    use ids_string              ! IGNORE
    use ids_grid_subgrid        ! IGNORE
    use ids_grid_objectlist     ! IGNORE
    use ids_grid_examples       ! IGNORE
    use ids_grid_unstructured   ! IGNORE
    use ids_grid_structured     ! IGNORE

#ifdef USE_PXFGETENV
    integer lenval, ierror
#endif
    implicit none

    external ipgeti, ipgetc

    !! Local variables
    character(len=24) :: treename   !< The name of the IMAS IDS database
        !< (i.e. "edge_profiles" (mandatory) )
    character(len=24) :: username   !< Creator/owner of the IMAS IDS database
    character(len=24) :: device     !< Device name of the IMAS IDS database
        !< (i. e. solps-iter, iter, aug)
    character(len=24) :: version    !< Major version of the IMAS IDS database
    integer :: idx  !< The returned identifier to be used in the subsequent
        !< data access operation
    integer :: shot !< The shot number of the database being created
    integer :: run  !< The run number of the database being created
    type(ids_edge_profiles) :: edge_profiles    !< IDS designed to store data on
        !< edge plasma profiles  (includes the scrape-off layer and possibly
        !< part of the confined plasma)
    type (ids_edge_sources) :: edge_sources !< IDS designed to store
        !< data on edge plasma sources. Energy terms correspond to the full
        !< kinetic energy equation (i.e. the energy flux takes into account
        !< the energy transported by the particle flux)
    type (ids_edge_transport) :: edge_transport !< IDS designed to store
        !< data on edge plasma transport. Energy terms correspond to the
        !< full kinetic energy equation (i.e. the energy flux takes into
        !< account the energy transported by the particle flux)
    character*256 systemarg

    !! Set default value for IMAS major version and IDS treename
    version = '3'
    treename = 'ids'
    write (*,*) 'Starting b2mn init'
    call b2mn_init
    ! call b2mn_step(0)

    call ipgeti( 'b2mndr_shot_number', shot )
    call xertst( 0.lt.shot.and.shot.le.214748, 'Invalid shot number')
    call ipgeti( 'b2mndr_run_number', run )
    call xertst( 0.le.run.and.run.le.9999, 'Invalid run number')
#ifdef NO_GETENV
    username=' '
#else
#ifdef USE_PXFGETENV
    CALL PXFGETENV ('USER', 0, username, lenval, ierror)
#else
    call getenv ('USER', username)
#endif
#endif
    call ipgetc( 'b2mndr_user', username )
    device = 'solps-iter'
#ifndef NO_GETENV
#ifdef USE_PXFGETENV
    CALL PXFGETENV ('DEVICE', 0, device, lenval, ierror)
#else
    call getenv ('DEVICE', device)
#endif
#endif
    call ipgetc( 'b2mndr_device', device )
    systemarg='imasdb '//trim(device)
#ifdef IMAS
    call system(systemarg)
#endif

    write(*,'(a,i8,a,i8,4a)') 'Shot: ', shot, ' Run: ', run, &
        & ' User: ', trim(username), ' Device: ', trim(device)

    !! Process B2.5 data and set it to IMAS IDS
    write(*,*) "START B25_process_ids"
    call B25_process_ids( edge_profiles, edge_sources, edge_transport )

    !! Create Write the set data to IDSs
    write(*,*) "START put_ids_edge"
    call put_ids_edge( edge_profiles, edge_sources, edge_transport, treename,   &
        &   shot, run, idx, username, device, version )

end program b2_ual_write

!!!Local Variables:
!!! mode: f90
!!! End:
