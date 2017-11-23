!!  Legend:
!!     !> ................ Documentation comment (file description, function
!!                         description etc.). Also intended for doxygen
!!                         generated documentation
!!     !> @note .......... Documentation notes, intended for doxygen generated
!!                         documentation
!!     !!  ............... variables description, additional (helpful)
!!                         information etc.
!!     ! IGNORE    ....... Used to ignore this module in list dependency when
!!                         building
!!     !   ............... Commented part of code
!!-----------------------------------------------------------------------------
!! DOCUMENTATION:
!> 1. purpose
!>
!>      b2_ual_write_b2mod.f90 script is used to generate b2_ual_write_b2mod.exe
!>      (main program), which is a post-processor for b2.
!>      The script currently reads the plasma grid
!>      geometry ( full geometry descriptions of all available grid subsets )
!>      and plasma state ( electron density ) and writes it to IDS database
!>      with the use of b2mod scripts that utilize IMAS GGD Grid Service
!>      Library routines.
!>
!>
!> 2. specification
!>
!>      Main program.
!>
!>
!> 3. description (see also routine b2cdca)
!>
!>      The complete program performs post-processing of the
!>      result of a b2 calculation.
!>      This program unit opens and closes the input/output units, and
!>      may perform some other system-dependent operations.
!>
!>      The input units are:
!>      ninp(0): formatted; provides output control parameters.
!>      ninp(1): un*formatted; provides the geometry.
!>      ninp(2): un*formatted; provides the run parameters.
!>      ninp(3): un*formatted; provides the plasma state.
!>      ninp(4): unformatted; provides the detailed plasma state.
!>      ninp(5): formatted; provides the run switches.
!>      ninp(6): un*formatted; provides the atomic data.
!>
!>      The output units are:
!>      nout(0): formatted; provides printed output.
!>
!>      (See routine b2cdca for the meaning of "un*formatted".)
!>
!>
!> 4. parameters (see also routine b2cdcv)
!>
!>      None.
!>
!>
!> 5. error indicators
!>
!>      In case an error condition is detected, a call is made to the
!>      routine xerrab. This causes an error message to be printed,
!>      after which the program halts.
!>
!!-----------------------------------------------------------------------------

program b2_ual_write_b2mod

    use b2mod_main
    use b2mod_grid_mapping
    use b2mod_ual_io

    use ids_schemas     ! IGNORE
                        !> These are the Fortran type definitions for the
                        !> Physics Data Model
    use ids_routines    ! IGNORE
                        !> These are the Access Layer routines + management of
                        !> IDS structures
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

    !!--------------------------------------------------------------------------

    !!.declarations

    !!..common blocks

    !!..local variables
    integer ::  idx, i
    integer ::  shot, run
    character(len=24)        ::  treename, username, device, version
    !> character(len=255)    ::  imas_connect_url
    type(ids_edge_profiles)  ::  edge_profiles
    type(ids_edge_sources)   ::  edge_sources
    type(ids_edge_transport) ::  edge_transport

    !!--------------------------------------------------------------------------
    !!.documentation-internal
    !!
    !!      The following common blocks have their outermost declaration in
    !!      this routine; they need not be preserved between calls.
    !!
    !!      Description of some local variables:
    !!
    !!      ninp - (0:6) integer array.
    !!      ninp specifies the input unit numbers.
    !!
    !!      nout - (0:2) integer array.
    !!      nout specifies the output unit numbers.
    !!
    !!      nx, ny - integer.
    !!      nx and ny specify the number of interior cells along the first
    !!      and the second coordinate, respectively. The total number of
    !!      cells is (nx+2)*(ny+2); they are indexed by (-1:nx,-1:ny).
    !!      It will hold that 0.le.nx and 0.le.ny.
    !!
    !!      ns - integer.
    !!      ns specifies the number of atomic species in the calculation.
    !!      The species are indexed by (0:ns-1).
    !!      It will hold that 1.le.ns.
    !!
    !!--------------------------------------------------------------------------
    !!.computation

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

    treename    = "ids"
    ! shot        = 16151
    ! run         = 1001
    shot    = 100
    run     = 7
    username    = "penkod"
    device     = "solps-iter"
    version     = "3"

    !! b2mod routine write_ids
    write(*,*) "START write_ids"
    call write_ids( edge_profiles, edge_sources, edge_transport )

    write(*,*) "START put_ids_edge"
    call put_ids_edge( edge_profiles, edge_sources, edge_transport, treename,   &
        &   shot, run, idx, username, device, version )

    ! call read_ids(treename, shot, run, idx, username, &
    !                                     & device, version )

contains

    !> Subroutine intended to check if supposed new file already exists
    !> and then delete it
    subroutine checkFileAndDelete( fileName )
        character(len=*), intent(in) :: fileName
        logical :: file_exists
        integer :: stat

        inquire( file=fileName, exist=file_exists )
        if ( file_exists ) then
            write(*,*) "Deleting old ", fileName
            open(unit=1234, iostat=stat, file=fileName, status='old')
            if (stat == 0) close(1234, status='delete')
        endif

    end subroutine

    !> Subroutine used to put data to edge_profiles IDS
    subroutine put_ids_edge( edge_profiles, edge_sources, edge_transport,   &
            &   treename, shot, run, idx, username, device, version )
        integer :: shot, run, idx
        type(ids_edge_profiles), intent(inout) :: edge_profiles
        type(ids_edge_sources), intent(inout) :: edge_sources
        type(ids_edge_transport), intent(inout) :: edge_transport
        character(len=24) :: treename, username, device, version

        !! Set data to edge_profiles IDS
        write(0,*) "Writing to edge_profiles, edge_sources and edge_transport IDS"

        !! Create and modify new shot/run
        call imas_create_env(treename, shot, run, 0, 0, idx, username, &
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
    subroutine read_ids( treename, shot, run, idx, username, device, version )
        !> Example for reading IDS database in Fortran90
        integer                 ::  shot, run, idx
        integer                 ::  gridSubset_index
        character(len=24)       ::  treename, username, device, version
        type(ids_edge_profiles) ::  edge_profiles
        character(len=255)      ::  imas_connect_url

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
