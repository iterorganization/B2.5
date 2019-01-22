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
    use b2mod_ual_io &
     & , only : b25_process_ids, ids_put, ids_deallocate, &
     &          ids_edge_profiles, ids_edge_sources, ids_edge_transport
#else
# ifdef ITM
    use euITM_schemas  ! IGNORE
    use euITM_routines ! IGNORE
# endif
#endif

  implicit none

  private

  public open_ual, close_ual
#ifdef IMAS
  public put_ids_edge
  public b25_process_ids
  public ids_edge_profiles, ids_edge_sources, ids_edge_transport
#endif


contains

#ifdef IMAS
    !> Subroutine used to put data to edge_profiles, edge_sources and
    !! edge_transport IDSs.
    subroutine put_ids_edge( edge_profiles, edge_sources, edge_transport,   &
            &   treename, shot, run, idx, username, device, version )
        type(ids_edge_profiles), intent(inout)  :: edge_profiles    !< IDS
            !< designed to store data on edge plasma profiles  (includes the
            !< scrape-off layer and possibly part of the confined plasma)
        type (ids_edge_sources), intent(inout)  :: edge_sources     !< IDS
            !< designed to store data on edge plasma sources. Energy terms
            !< correspond to the full kinetic energy equation (i.e. the energy
            !< flux takes into account the energy transported by the particle
            !< flux)
        type (ids_edge_transport), intent(inout)  :: edge_transport !< IDS
            !< designed to store  data on edge plasma transport. Energy terms
            !< correspond to the full kinetic energy equation (i.e. the energy
            !< flux takes into account the energy transported by the particle
            !< flux)
        character(len=24), intent(in) :: treename   !< The name of the IMAS IDS database
            !< (i.e. "edge_profiles" (mandatory) )
        integer, intent(in) :: shot !< The shot number of the database being created
        integer, intent(in) :: run  !< The run number of the database being created
        integer, intent(in) :: idx  !< The returned identifier to be used in the subsequent
            !< data access operation
        character(len=24), intent(in) :: username   !< Creator/owner of the IMAS IDS
            !< database
        character(len=24), intent(in) :: device     !< Device name of the IMAS IDS database
            !< (i. e. solps-iter, iter, aug)
        character(len=24), intent(in) :: version    !< Major version of the IMAS IDS
            !< database

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
#endif

    !> Routine to open UAL database.
    !! @note For IMAS IDS is recommended use of IMAS GGD library routine
    !! "exampleOpenIDS"
    subroutine open_ual( idx, shot, run, time, user, tokamak, dataversion,  &
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
        character(*), intent(in), optional :: tokamak !< Device name of the
            !< database (i. e. solps-iter, iter, aug)
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
        character(32) :: lTreename = "ids"
#else
# ifdef ITM
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

        namelist /ual_namelist/ lTreename, lShot, lRun, lTime, lRefshot,    &
            &   lRefrun, lUser, lTokamak, lDataversion, openEnv, lDoCreate, &
            &   lUseHdf5

        if( present( shot ) ) lShot = shot
        if( present( run ) ) lRun = run
        if( present( user ) ) lUser = user
        if( present( tokamak ) ) lTokamak = tokamak
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
                call imas_create_hdf5(lTreename, lShot, lRun, lRefshot, &
                        &   lRefrun, idx)
            else
                if( openEnv) then
                    call imas_create_env(lTreename, lShot, lRun, lRefshot,  &
                        &   lRefrun, idx, lUser, lTokamak, lDataversion)
                else
                    call imas_create(lTreename, lShot, lRun, lRefshot, &
                        &   lRefrun, idx)
                end if
            end if
        else
            if( lUseHdf5) then
                call imas_open_hdf5(lTreename, lShot, lRun, idx)
            else
                if( openEnv) then
                    call imas_open_env(lTreename, lShot, lRun, idx, lUser, &
                        &   lTokamak, lDataversion)
                else
                    call imas_open(lTreename, lShot, lRun, lRefshot, &
                        &   lRefrun, idx)
                end if
            end if
#else
#ifdef ITM
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
#else
            idx = 0
#endif
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
        call imas_close(idx)
#else
#ifdef ITM
        call euITM_close(idx)
#endif
#endif
    end subroutine close_ual

end module b2mod_ual

!!!Local Variables:
!!! mode: f90
!!! End:
