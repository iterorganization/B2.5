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

module b2mod_ual

  !! Basic UAL routines shared by the entire SOLPS application chain

  use b2mod_types
#ifdef IMAS
  use ids_schemas  ! IGNORE
  use ids_routines ! IGNORE
#else
# ifdef ITM
  use euITM_schemas  ! IGNORE
  use euITM_routines ! IGNORE
# endif
#endif

  implicit none

  private

  public open_ual, close_ual


contains

  !> Routine to open UAL database.
  !> @note For IMAS IDS is recommended use of IMAS GGD library routine
  !> "exampleOpenIDS"
  subroutine open_ual(idx, shot, run, time, user, tokamak, dataversion, &
        &   doCreate, useHdf5, nmlFile)
    integer, intent(out) :: idx
    integer, intent(in), optional :: shot, run
    real(R8), intent(out), optional :: time     !> Time is special: it is not
                                                !> used here, but can be read
                                                !> from the namelist and
                                                !> returned
    character(*), intent(in), optional :: user
    character(*), intent(in), optional :: tokamak
    character(*), intent(in), optional :: dataversion
    logical, intent(in), optional :: doCreate
    logical, intent(in), optional :: useHdf5
    character(*), intent(in), optional :: nmlFile

    !! internal

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
        &   lRefrun, lUser, lTokamak, lDataversion, openEnv, lDoCreate, lUseHdf5

    if (present(shot)) lShot = shot
    if (present(run)) lRun = run
    if (present(user)) lUser = user
    if (present(tokamak)) lTokamak = tokamak
    if (present(dataversion)) lDataversion = dataversion
    if (present(doCreate)) lDoCreate = doCreate
    if (present(useHdf5)) lUseHdf5 = useHdf5

    if (present(user)) openEnv = .true.

    !! If file exists, read namelist from configuration file
    !! If not, write out namelist
    if (present(nmlFile)) then
        inquire(file=nmlFile, exist=namelistExists)
    else
        inquire(file=NAMELIST_FILE, exist=namelistExists)
    end if
    if (namelistExists) then
        if (present(nmlFile)) then
            open(unit=NAMELIST_UNIT, file=nmlFile)
        else
            open(unit=NAMELIST_UNIT, file=NAMELIST_FILE)
        end if
        read (NAMELIST_UNIT, nml=ual_namelist)
        close(unit=NAMELIST_UNIT)
    else
        if (present(nmlFile)) then
            open(unit=NAMELIST_UNIT, file=nmlFile, status="new", action="write")
        else
            open(unit=NAMELIST_UNIT, file=NAMELIST_FILE, status="new",  &
                &   action="write")
        end if
        write (NAMELIST_UNIT, nml=ual_namelist)
        close(unit=NAMELIST_UNIT)
    end if

    !! establish UAL access
    if (lDoCreate) then
#ifdef IMAS
        if (lUseHdf5) then
            call imas_create_hdf5(lTreename, lShot, lRun, lRefshot, lRefrun, idx)
        else
            if (openEnv) then
                call imas_create_env(lTreename, lShot, lRun, lRefshot,  &
                    &   lRefrun, idx, lUser, lTokamak, lDataversion)
            else
                call imas_create(lTreename, lShot, lRun, lRefshot, lRefrun, idx)
            end if
        end if
    else
        if (lUseHdf5) then
            call imas_open_hdf5(lTreename, lShot, lRun, idx)
        else
            if (openEnv) then
                call imas_open_env(lTreename, lShot, lRun, idx,  lUser, &
                    &   lTokamak, lDataversion)
            else
                call imas_open(lTreename, lShot, lRun, lRefshot, lRefrun, idx)
            end if
        end if
#else
#ifdef ITM
        if (lUseHdf5) then
            call euITM_create_hdf5(lTreename, lShot, lRun, lRefshot, lRefrun, idx)
        else
            if (openEnv) then
                call euITM_create_env(lTreename, lShot, lRun, lRefshot, &
                    &   lRefrun, idx, lUser, lTokamak, lDataversion)
            else
                call euITM_create(lTreename, lShot, lRun, lRefshot, lRefrun, idx)
            end if
        end if
    else
        if (lUseHdf5) then
            call euITM_open_hdf5(lTreename, lShot, lRun, idx)
        else
            if (openEnv) then
                call euITM_open_env(lTreename, lShot, lRun, idx,  lUser,    &
                    &   lTokamak, lDataversion)
            else
                call euITM_open(lTreename, lShot, lRun, lRefshot, lRefrun, idx)
            end if
        end if
#else
        idx = 0
#endif
#endif
    end if

    !! Return time if requested
    if (present(time)) time = lTime

  end subroutine open_ual

  !> Close UAL.
  !> @note  For IMAS edge_profiles IDS "examplePutIDS" IMAS GGD library routine
  !> can be used instead (that routine also writes the set data to IDS and then
  !> closes the IDS)
  subroutine close_ual(idx)
    integer, intent(in) :: idx

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
