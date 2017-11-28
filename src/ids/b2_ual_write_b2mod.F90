!!-----------------------------------------------------------------------------
!! DOCUMENTATION (doxygen 1.8.8):
!>      @author
!>      Dejan Penko
!!
!>      @section desc   Description
!!      b2_ual_write_b2mod.f90 script is used to generate b2_ual_write_b2mod.exe
!!      (main program), which is a post-processor for b2.
!!      The script currently reads the plasma grid
!!      geometry ( full geometry descriptions of all available grid subsets )
!!      and plasma state ( electron density ) and writes it to IDS database
!!      with the use of b2mod scripts that utilize IMAS GGD Grid Service
!!      Library routines.
!!
!!      @subsection det   Details
!!      For more information see also routine b2cdca.
!!
!!      The complete program performs post-processing of the
!!      result of a b2 calculation.
!!      This program unit opens and closes the input/output units, and
!!      may perform some other system-dependent operations.
!!
!!      The input units are:
!!          - ninp(0): formatted; provides output control parameters.
!!          - ninp(1): un*formatted; provides the geometry.
!!          - ninp(2): un*formatted; provides the run parameters.
!!          - ninp(3): un*formatted; provides the plasma state.
!!          - ninp(4): unformatted; provides the detailed plasma state.
!!          - ninp(5): formatted; provides the run switches.
!!          - ninp(6): un*formatted; provides the atomic data.
!!
!!      The output units are:
!!          - nout(0): formatted; provides printed output.
!!
!!      @note   See routine b2cdca for the meaning of "un*formatted".
!!
!!      @subsection pv    Parameters/variables
!!      @note   see also routine b2cdcv
!!
!!      @param  device - Device name of the IMAS IDS database
!!              (i. e. solps-iter, iter, aug)
!!      @param  edge_profiles - IDS designed to store data on edge plasma
!!              profiles  (includes the scrape-off layer and possibly part
!!              of the confined plasma)
!!      @param  edge_sources - IDS designed to store data on edge plasma
!!              sources. Energy terms correspond to the full kinetic energy
!!              equation (i.e. the energy flux takes into account the energy
!!              transported by the particle flux)
!!      @param  edge_transport - IDS designed to store data on edge plasma
!!              transport. Energy terms correspond to the full kinetic energy
!!              equation (i.e. the energy flux takes into account the energy
!!              transported by the particle flux)
!!      @param  idx - The returned identifier to be used in the subsequent
!!                    data access operation
!!      @param  run - The run number of the database being created
!!      @param  shot - The shot number of the database being created
!!      @param  treename - the name of the IMAS IDS database,
!!              (i.e. "edge_profiles" (mandatory) )
!!      @param  username - Creator/owner of the IMAS IDS database
!!      @param  version - Major version of the IMAS IDS database
!!
!!      @subsection eind  Error indicators
!!      In case an error condition is detected, a call is made to the
!!      routine \b xerrab. This causes an error message to be printed,
!!      after which the program halts.
!!
!!      @subsection syx     Exceptional syntax explanation
!!      @code
!!          ! IGNORE    !! syntax used to ignore this module in list
!!                      !! dependency when compiling the code
!!      @endcode
!!
!!-----------------------------------------------------------------------------

program b2_ual_write_b2mod

    use b2mod_main
    use b2mod_grid_mapping
    use b2mod_ual_io
    use ids_schemas     ! IGNORE
                        !! These are the Fortran type definitions for the
                        !! Physics Data Model
    use ids_routines    ! IGNORE
                        !! These are the Access Layer routines + management of
                        !! IDS structures
    use ids_assert      ! IGNORE
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

    implicit none

    !! Local variables
    character(len=24) :: treename, username, device, version
    integer :: idx, i
    integer :: shot, run
    type(ids_edge_profiles) :: edge_profiles
    type(ids_edge_sources) :: edge_sources
    type(ids_edge_transport) :: edge_transport

    !! Check if supposed new file already exists and delete it
    call checkFileAndDelete( "b2fparam" )
    call checkFileAndDelete( "b2mn.prt" )
    call checkFileAndDelete( "b2fstate" )
    call checkFileAndDelete( "b2fmovie" )
    call checkFileAndDelete( "b2ftrace" )
    call checkFileAndDelete( "b2ftrack" )

    !! Run main b2 routine to process and read the b2 data
    write(0,*) "Running b2mn_init"
    call b2mn_init
    write(0,*) "b2mn_init completed"

    write(0,*) "Running b2mn_step(0)"
    call b2mn_step(0)
    write(0,*) "b2mn_step(0) completed"

    ! write(0,*) " Running b2mn_fin"
    ! call b2mn_fin
    ! write(0,*) "b2mn_fin completed"

    treename = "ids"
    shot = 100
    run = 7
    username = "penkod"
    device = "solps-iter"
    version = "3"

    !! b2mod routine write_ids
    write(*,*) "START write_ids"
    call write_ids( edge_profiles, edge_sources, edge_transport )

    write(*,*) "START put_ids_edge"
    call put_ids_edge( edge_profiles, edge_sources, edge_transport, treename,   &
        &   shot, run, idx, username, device, version )

    ! call read_ids(treename, shot, run, idx, username, &
    !                                     & device, version )

contains

    !> Subroutine intended to check if supposed new file already exists. If
    !! the file exists it deletes it.
    !! @param   filename - Name of the file to be checked
    subroutine checkFileAndDelete( fileName )
        !! Internal variables
        character(len=*), intent(in) :: fileName
        !! Local variables
        integer :: stat
        logical :: file_exists

        inquire( file=fileName, exist=file_exists )
        if ( file_exists ) then
            write(*,*) "Deleting old ", fileName
            open(unit=1234, iostat=stat, file=fileName, status='old')
            if (stat == 0) close(1234, status='delete')
        endif

    end subroutine

    !> Subroutine used to put data to edge_profiles, edge_sources and
    !! edge_transport IDSs.
    !!      @param  edge_profiles - IDS designed to store data on edge plasma
    !!              profiles  (includes the scrape-off layer and possibly part
    !!              of the confined plasma)
    !!      @param  edge_sources - IDS designed to store data on edge plasma
    !!              sources. Energy terms correspond to the full kinetic energy
    !!              equation (i.e. the energy flux takes into account the energy
    !!              transported by the particle flux)
    !!      @param  edge_transport - IDS designed to store data on edge plasma
    !!              transport. Energy terms correspond to the full kinetic energy
    !!              equation (i.e. the energy flux takes into account the energy
    !!              transported by the particle flux)
    !!      @param  treename - the name of the IMAS IDS database,
    !!              (i.e. "edge_profiles" (mandatory) )
    !!      @param  shot - The shot number of the database being created
    !!      @param  run - The run number of the database being created
    !!      @param  idx - The returned identifier to be used in the subsequent
    !!              data access operation
    !!      @param  username - Creator/owner of the IMAS IDS database
    !!      @param  device - Device name of the IMAS IDS database
    !!              (i. e. solps-iter, iter, aug)
    !!      @param  version - Major version of the IMAS IDS database
    subroutine put_ids_edge( edge_profiles, edge_sources, edge_transport,   &
            &   treename, shot, run, idx, username, device, version )
        !! Internal variables
        integer :: shot, run, idx
        type(ids_edge_profiles), intent(inout) :: edge_profiles
        type(ids_edge_sources), intent(inout) :: edge_sources
        type(ids_edge_transport), intent(inout) :: edge_transport
        character(len=24) :: treename, username, device, version

        !! Set data to edge_profiles IDS
        write(0,*) "Writing to edge_profiles, edge_sources and edge_transport IDS"

        !! Create and modify new shot/run
        call imas_create_env( treename, shot, run, 0, 0, idx, username, &
            device, version )

        !! Or open and modify existing shot/run (might work much faster than
        !! imas_create_env)
        ! call imas_open_env('treename', shot, run, idx, username, device, version )

        !! Put data to IDS
        ! call ids_put_slice( idx, "edge_profiles", edge_profiles )
        ! call ids_put_slice( idx, "edge_transport", edge_sources )
        ! call ids_put_slice( idx, "edge_transport", edge_transport )
        call ids_put( idx, "edge_profiles", edge_profiles )
        call ids_put( idx, "edge_sources", edge_sources )
        call ids_put( idx, "edge_transport", edge_transport )

        !! Close IDS
        call ids_deallocate( edge_profiles )
        call ids_deallocate( edge_sources )
        call ids_deallocate( edge_transport )
        call imas_close( idx )

        write(0,*) "IDS write finished"

    end subroutine put_ids_edge

    !> Example subroutine for reading edge_profiles IDS
    !! with Fortran90
    !!      @param  treename - the name of the IMAS IDS database,
    !!              (i.e. "edge_profiles" (mandatory) )
    !!      @param  shot - The shot number of the database being created
    !!      @param  run - The run number of the database being created
    !!      @param  idx - The returned identifier to be used in the subsequent
    !!                    data access operation
    !!      @param  username - Creator/owner of the IMAS IDS database
    !!      @param  device - Device name of the IMAS IDS database
    !!              (i. e. solps-iter, iter, aug)
    !!      @param  version - Major version of the IMAS IDS database
    subroutine read_ids( treename, shot, run, idx, username, device, version )
        !! Internal variables
        character(len=24)       ::  treename, username, device, version
        integer                 ::  shot, run, idx
        integer                 ::  gridSubset_index
        !! Local variables
        type(ids_edge_profiles) ::  edge_profiles

        gridSubset_index = 3

        !! Open input datafile from local database
        write (0,*) "Started reading input IDS", idx, shot, run

        call imas_open_env('treename', shot, run, idx, username, device, version )
        call ids_get(idx, "edge_profiles", edge_profiles)

        write(0,*) "homogeneous_time = ",   &
            &   edge_profiles%ids_properties%homogeneous_time
        write(0,*) "Grid subset 3 name = ", edge_profiles%ggd(1)%grid%  &
            &   grid_subset(gridSubset_index)%identifier%name
        write(0,*) "Grid subset 3 index = ", edge_profiles%ggd(1)%grid% &
            &   grid_subset(gridSubset_index)%identifier%index
        ! write(0,*) "Time = ", edge_profiles%time(1)
        call ids_deallocate( edge_profiles )
        call imas_close( idx )
        write (0,*) "Finished reading input IDS"

    end subroutine read_ids

end program b2_ual_write_b2mod

!!!Local Variables:
!!! mode: f90
!!! End:
