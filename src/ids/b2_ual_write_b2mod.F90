!!-----------------------------------------------------------------------------
!! DOCUMENTATION (doxygen 1.8.8):
!>      @author
!>      Dejan Penko
!!
!>      @page b2uw_b2mod b2_ual_write_b2mod
!>      @section b2uw_b2mod_desc   Description
!!      b2_ual_write_b2mod code is used to generate b2_ual_write_b2mod.exe
!!      (main program), which is a post-processor for B2.
!!      The code reads the plasma grid
!!      geometry (full geometry descriptions of all available grid subsets)
!!      and plasma state (electron density/temperature, ion temperature,
!!      velocity etc.). The code then writes the obtained data to IDS database
!!      with the use of b2mod scripts that utilize IMAS GGD Grid Service
!!      Library routines.
!!
!!      @note   More on the b2_ual writers is available in SOLPS-GUI
!!              documentation \b HOWTOs under section <b> 4.5 IMAS </b>.
!!
!!      @note   A short video tutorial on the use of the B2.5 writer is
!!              available <a href="https://youtu.be/5IuADXPAgkQ">here</a>.
!!
!!      @section b2uw_b2mod_cpo2ids Mapping CPO -> IDS
!!      Most B2.5 routines were originally developed to work solely with ITM
!!      CPOs. Those same routines were modified or created anew, providing
!!      necessary tools for working with IMAS IDSs.
!!      Since CPO and IDS data structure is not the same a lot of proper
!!      adjustments had to be made, mostly in modules
!!      @ref b2uw_ualio_grid_desc "b2mod_ual_io_grid" and
!!      @ref b2uw_ualio_desc      "b2mod_ual_io".
!!
!!
!!      @subsection b2uw_b2mod_cpo2ids_constants Constants, classes etc.
!!      Below are listed CPO constants, classes etc. and corresponding IDS
!!      ones that were used in IMAS IDS B2.5 routines.
!!
!!      CPO                                     | IDS
!!      -------------------------------         | ------------------
!!      CLASS_NODE = (/ 0, 0 /)                 | IDS_CLASS_NODE = 1
!!      CLASS_RZ_EDGE = (/ 1, 0 /)              | IDS_CLASS_RZ_EDGE = 2
!!      CLASS_PHI_EDGE = (/ 0,10 /)             | IDS_CLASS_PHI_EDGE = 2
!!      CLASS_POLOIDALRADIAL_EDGE = (/ 1, 1 /)  | IDS_CLASS_POLOIDALRADIAL_EDGE = 2
!!      CLASS_TOROIDAL_EDGE = (/ 2, 0 /)        | IDS_CLASS_TOROIDAL_EDGE = 2
!!      CLASS_CELL = (/ 2, 1 /)                 | IDS_CLASS_CELL = 3
!!
!!      <b> Grid subset IDs </b>:
!!
!!      B2.5 ITM routines use grid subset IDs (B2_SUBGRID_UNSPECIFIED,
!!      B2_SUBGRID_NODES, B2_SUBGRID_CELLS, etc.) defined in
!!      @ref b2uw_ualio_grid_desc "b2mod_ual_io_grid.F90", while B2.5 IDS
!!      uses grid subset IDs defined in IMAS GGD (ids_grid_common.f90).
!!
!!      @subsection b2uw_b2mod_cpo2ids_nodes Data tree nodes
!!      Below are listed CPO nodes and corresponding IDS nodes to which the
!!      data was written instead.
!!
!!      CPO edge.grid. ...              | IDS edge_profiles.ggd(:).grid. ...
!!      ------------------------------- | -------------------------------------
!!      spaces(:).coordtype             | space(:).coordinates_type
!!      spaces(:).objects               | space(:).objects_per_dimension(:).object
!!      spaces(:).objects(:).geo        | space(:).objects_per_dimension(:).object(:).geometry
!!      spaces(:).objects(:).boundary   | space(:).objects_per_dimension(:).object(:).boundary
!!      spaces(:).objects(:).neighbour  | space(:).objects_per_dimension(:).object(:).boundary(:).neighbours
!!      spaces(:).xpoints               | No node for data on x-points was found
!!      subgrids                        | grid_subset
!!
!!      CPO edge. ...           | IDS
!!      ----------------------- | ---------------------------------------------
!!      fluid.ne.value          | edge_profiles.ggd(:).electrons.density
!!      fluid.ne.flux           | edge_transport.model(:).ggd(:).electrons.particles.flux
!!      fluid.ne.source         | edge_sources.source(:).ggd(:).electrons.particles
!!      fluid.ni.value          | edge_profiles.ggd(:).ion(:).density
!!      fluid.ni.flux           | edge_transport.model(:).ggd(:).ion(:).particles.flux
!!      fluid.ni.source         | edge_sources.source(:).ggd(:).ion(:).particles
!!      fluid.vi(:).comps(1)    | edge_profiles.ggd(:).ion(:).velocity(:).radial
!!      fluid.vi(:).comps(2)    | edge_profiles.ggd(:).ion(:).velocity(:).poloidal
!!      fluid.vi(:).comps(3)    | edge_profiles.ggd(:).ion(:).velocity(:).toroidal
!!      fluid.vi(:).align(:)    | Probably not required as the leaf label refers to vector component itself
!!      fluid.vi(:).alignid     | Refers to the label of the node velocity(:) leaf
!!      fluid.te(:).value       | edge_profiles.ggd(:).electrons.temperature
!!      fluid.te(:).flux        | edge_transport.model(:).ggd(:).electrons.energy.flux
!!      fluid.ti(:).value       | edge_profiles.ggd(:).ion(:).temperature
!!      fluid.ti(:).flux        | edge_transport.model(:).ggd(:).ion(:).energy.flux
!!      fluid.po.value          | edge_profiles.ggd(:).phi_potential
!!      fluid.te_aniso.comps(1) | edge_profiles.ggd(:).e_field.poloidal
!!      fluid.te_aniso.comps(2) | edge_profiles.ggd(:).e_field.radial
!!      fluid.te_aniso.comps(3) | edge_profiles.ggd(:).e_field.toroidal
!!      fluid.te_aniso.comps(4) | edge_profiles.ggd(:).e_field.diamagnetic
!!
!!      @note   In the future, IDS data structure nodes that correspond to
!!              flux data fields are to be moved from edge_transport IDS to
!!              edge_profiles IDS.
!!
!!      @section b2uw_b2mod_env Compiling and setting the environment
!!      In order compile the B2.5 writer code use the commands below:
!!
!!      @verbatim
!!          cd $HOME/solps-iter
!!          tcsh
!!          source setup.csh
!!          cd modules/B2.5
!!          make ids
!!      @endverbatim
!!
!!      @note   At the time of writing this manual IMAS module
!!              imas/3.13.0/ual/3.6.3 was used. This IMAS module provides
!!              also GGD support as it includes IMAS GGD library routines
!!              (Fortran90).
!!
!!      @note   b2_ual_write_deprecated and b2_ual_write_gsl are OUTDATED codes,
!!              but were left in the repository for documentation purposes and
!!              as extra examples.
!!
!!      @subsection b2uw_b2mod_run Running the code:
!!      The examples are available on ITER portal:
!!      <a href="https://portal.iter.org/departments/POP/CM/IMAS/Forms/AllItems.aspx?RootFolder=%2Fdepartments%2FPOP%2FCM%2FIMAS%2FSOLPS-ITER%2FExamples">B2.5 examples link</a>
!!
!!      In terminal navigate to directory containing the case required data files
!!      (b2fgmtry, b2fstate etc.) and run the following command:
!!
!!      @verbatim
!!          $SOLPSTOP/modules/B2.5/builds/standalone.$HOST_NAME.$COMPILER/b2_ual_write_b2mod.exe
!!          --pulse <pulse> --run <run> --username <username> --database <database>
!!          --version <version> --step
!!      @endverbatim
!!
!!      The arguments marked with < ... > are the parameters of the IDS database
!!      where the data is to be stored:
!!          - \b path:     The path to the IMAS data entry being created
!!          - \b pulse:    The pulse (previously shot) number of the database being created
!!          - \b run:      The run number of the database being created
!!          - \b username: Creator/owner of the IMAS IDS database
!!          - \b database: IMAS IDS database name (i. e. solps-iter, ITER, aug)
!!          - \b version:  Major version of the IMAS IDS database
!!
!!      Example of the command:
!!      @verbatim
!!          $SOLPSTOP/modules/B2.5/builds/standalone.$HOST_NAME.$COMPILER/b2_ual_write_b2mod.exe
!!          --pulse 1512 --run 6 --username penkod --database solps-iter --version 3 --step
!!      @endverbatim
!!
!!      \b References:
!!          - @ref b2uw_b2mod_prog "b2_ual_write_b2mod file reference"
!!          - @ref b2uw_ualio_grid_desc "b2mod_ual_io_grid module "
!!          - @ref b2uw_ualio_desc      "b2mod_ual_io module".
!!
!!-----------------------------------------------------------------------------

!!-----------------------------------------------------------------------------
!>      @section b2uw_b2mod_prog Program b2_ual_write_b2mod
!!      References:
!!          - @ref b2uw_b2mod "b2_ual_write_b2mod main page"
!!
!!      @subsection b2uw_b2mod_det   Details
!!      For more information see also routine \b b2cdca.
!!
!!      The complete program performs post-processing of the
!!      result of a B2 calculation.
!!      This program unit opens and closes the input/output units, and
!!      may perform some other system-dependent operations.
!!
!!      The input units are:
!!          - ninp(0): formatted, provides output control parameters.
!!          - ninp(1): un*formatted, provides the geometry.
!!          - ninp(2): un*formatted, provides the run parameters.
!!          - ninp(3): un*formatted, provides the plasma state.
!!          - ninp(4): unformatted, provides the detailed plasma state.
!!          - ninp(5): formatted, provides the run switches.
!!          - ninp(6): un*formatted, provides the atomic data.
!!
!!      The output units are:
!!          - nout(0): formatted, provides printed output.
!!
!!      @note   See routine b2cdca for the meaning of "un*formatted".
!!
!!      @subsection b2uw_b2mod_pv    Parameters/variables
!!      @note   see also routine \b b2cdcv
!!
!!      @param  database - IMAS IDS database name (i. e. solps-iter, ITER, aug)
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
!!              data access operation
!!      @param  run - The run number of the database being created
!!      @param  shot - The pulse (previously shot) number of the database
!!              being created
!!      @param  treename - The name of the IMAS IDS database
!!              (local $DEVICE environment variable by default if defined,
!!               otherwise "solps-iter" if argument is not provided)
!!      @param  username - Creator/owner of the IMAS IDS database
!!      @param  version - Major version of the IMAS IDS database
!!
!!      @subsection b2uw_b2mod_eind  Error indicators
!!      In case an error condition is detected, a call is made to the
!!      routine \b xerrab. This causes an error message to be printed,
!!      after which the program halts.
!!
!!      @subsection b2uw_b2mod_syx     Exceptional syntax explanation
!!      @code
!!          ! IGNORE    !! syntax used to ignore this module in list
!!                      !! dependency when compiling the code
!!      @endcode
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
!!      integer :: shot   !< The pulse (previsouly shot) number of the database
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

program b2_ual_write_b2mod

    use b2mod_main &
     & , only : b2mn_init, b2mn_step, b2mn_fin
    use b2mod_driver &
     & , only : idx, ids_path, dtim, imas_version, continued, &
     &          shot, run, username, database, version, &
     &          description, old_imas_version, &
     &          old_description, old_edge_profiles, equilibrium, wall, &
     &          ids_end_time, old_start_time, old_end_time, new_eq_ggd, &
     &          edge_profiles, edge_sources, edge_transport, radiation, &
     &          batch_profiles, batch_sources
    use b2mod_time
    use b2mod_ipmain
    use b2mod_solpstop
    use b2mod_switches
    use b2us_geo
    use b2us_map
    use b2us_data
    use b2us_plasma
    use b2mod_grid_mapping
    use ids_routines &  ! IGNORE
     & , only : imas_create_env
    use ids_schemas &   ! IGNORE
     & , only : ids_edge_profiles, ids_edge_sources, ids_edge_transport, &
     &          ids_radiation, ids_equilibrium, ids_wall
    use b2mod_ual &
     & , only : new_ids_edge, put_ids_edge, close_ual, &
     &          dealloc_ids_edge, dealloc_batch_edge
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
#if ( IMAS_MAJOR_VERSION < 4 || ( IMAS_MAJOR_VERSION == 4 && IMAS_MINOR_VERSION < 1 ) )
    use ids_schemas &   ! IGNORE
     & , only : ids_dataset_description
    use b2mod_driver &
     & , only : old_description, description
#else
    use b2mod_driver &
     & , only : old_summary
#endif
#if ( IMAS_MINOR_VERSION > 11 && IMAS_MINOR_VERSION < 15 && IMAS_MAJOR_VERSION == 3 )
    use ids_grid_examples       ! IGNORE
#endif
#if AL_MAJOR_VERSION > 4
    use ids_routines &  ! IGNORE
     & , only : imas_open, al_build_uri_from_legacy_parameters, &
     &          OPEN_PULSE, STRMAXLEN, MDSPLUS_BACKEND
#else
    use ids_routines &  ! IGNORE
     & , only : imas_open_env
#endif
    use ids_routines &  ! IGNORE
     & , only : ids_deallocate, ids_get
    use ids_schemas  &  ! IGNORE
     & , only : IDS_real, IDS_REAL_INVALID
    use b2mod_ad &
        & , only : nncf
#ifdef B25_EIRENE
    use eirmod_parmmod
    use eirmod_comusr
    use eirmod_extrab25
#endif

    implicit none
#ifndef NO_GETENV
    character(len=24) :: device_env
    integer lenval, ierror
#ifndef USE_PXFGETENV
    intrinsic get_environment_variable
#endif
#endif

    !! Local variables
    integer :: narg     !< Total Number of input arguments (shot, run, etc.)
    integer :: cptArg
    integer :: l, m, num_time_slices, tmp_run, time_slice_index, status

    !! Dummy variables
    character(len=16) :: run_user
    character(len=24) :: shot_string
    character(len=24) :: run_string
    character(len=24) :: argName
    character(len=24) :: treename
    character(len=256) :: systemarg, imasdir, home_dir, path
    real(kind=r8) :: J(nncf)
    logical :: new_run, absolute_path
#if AL_MAJOR_VERSION > 4
    character*256 olddir
    character(len=STRMAXLEN) :: uri, uri_source, uri_dest
    character(len=:), allocatable :: message
#endif

    !! Procedures
    character*16 usrnam
    logical streql
    external usrnam, streql

    !! Check if supposed new file already exists and delete it
    call checkFileAndDelete( "b2fparam" )
    call checkFileAndDelete( "b2mn.prt" )
    call checkFileAndDelete( "b2fstate" )
    call checkFileAndDelete( "b2fmovie" )
    call checkFileAndDelete( "b2ftrace" )
    call checkFileAndDelete( "b2ftrack" )
    call checkFileAndDelete( "b2fplasma" )

    !! Run main B2 routine to process and read the B2 data
    write(0,*) "Running b2mn_init"
    call b2mn_init
    write(0,*) "b2mn_init completed"

    !! Set default value for IMAS major version
    new_run = .false.
    status = 0
    write(version,'(i1)') IMAS_MAJOR_VERSION
    treename = 'ids'
    run_user = usrnam()
    username = run_user
    call ipgetc('b2mndr_user', username )
    call xertst( .not.streql(username,' '), 'User name not defined !')
    home_dir = '/home/'//trim(run_user)
    database = 'solps-iter'
    path = ' '
#ifdef NO_GETENV
    write(imas_version,'(i1,a1,i2,a1,i1)') IMAS_MAJOR_VERSION,'.', &
                                         & IMAS_MINOR_VERSION,'.', &
                                         & IMAS_MICRO_VERSION
#else
    device_env = ' '
#ifdef USE_PXFGETENV
    CALL PXFGETENV ('HOME', 0, home_dir, lenval, ierror)
    CALL PXFGETENV ('DEVICE', 0, device_env, lenval, ierror)
    CALL PXFGETENV ('IMAS_VERSION', 0, imas_version, lenval, ierror)
#else
    call get_environment_variable('HOME', status=ierror, length=lenval)
    if (ierror.eq.0) call get_environment_variable('HOME', value=home_dir)
    call get_environment_variable('DEVICE', status=ierror, length=lenval)
    if (ierror.eq.0) call get_environment_variable('DEVICE', value=device_env)
    call get_environment_variable('IMAS_VERSION', status=ierror, length=lenval)
    if (ierror.eq.0) call get_environment_variable('IMAS_VERSION', value=imas_version)
#endif
    if (.not.streql(device_env,' ')) database = device_env
    if (streql(database,'iter')) database = 'ITER'
#endif
#ifndef NO_GETENV
#ifdef USE_PXFGETENV
    CALL PXFGETENV ('IMASDIR', 0, imasdir, lenval, ierror)
#else
    call get_environment_variable('IMASDIR', status=ierror, length=lenval)
    if (ierror.eq.0) call get_environment_variable('IMASDIR', value=imasdir)
#endif
    call ipgetc('b2mndr_device', database )
    call ipgetc('b2mndr_database', database )
#endif
#if AL_MAJOR_VERSION > 4
    imasdir = trim(home_dir)//'/public/imasdb/'//trim(database)//'/'//trim(version)
#else
    imasdir = trim(home_dir)//'/public/imasdb/'//trim(database)//'/'//trim(version)//'/0'
#endif

    !! Check if arguments are found
    narg = command_argument_count()

    if( narg > 0 ) then
        !! Loop across arguments and allocate proper values to variables
        do cptArg = 1, narg
            call get_command_argument( cptArg, argName )
            select case( adjustl( argName ) )
#if AL_MAJOR_VERSION > 4
                case("--path","-p")
                    call get_command_argument( cptArg + 1, path )
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
                case("--device","--database","-d")
                    call get_command_argument( cptArg + 1, database )
                case("--version","-v")
                    call get_command_argument( cptArg + 1, version )
                case("--step","-S")
                    new_run = .true.
            end select
        end do
    end if
    if (narg.lt.2) then
        write(0,*) "ERROR! In order to run b2_ual_write_b2mod "// &
            & "input IMAS data entry "// &
#if AL_MAJOR_VERSION > 4
            & "user, database and version and either "// &
            & "pulse and run or path variables must "// &
            & "be defined. Example (terminal): "
        write(0,*) "$SOLPSTOP/modules/B2.5/builds/standalone.ITER.ifort64/"// &
            & "b2_ual_write_b2mod.exe --pulse 1512 --run 6 "// &
            & "--username penkod --database solps-iter --version 3 --step"
#else
            & "user, database, version, shot and run variables must "// &
            & "be defined. Example (terminal): "
        write(0,*) "$SOLPSTOP/modules/B2.5/builds/standalone.ITER.ifort64/"// &
            & "b2_ual_write_b2mod.exe --shot 1512 --run 6 "// &
            & "--username penkod --database solps-iter --version 3 --step"
#endif
        call exit(0)
    end if

#if AL_MAJOR_VERSION > 4
    call xertst( 0.lt.shot.or..not.streql(path,' '), 'Invalid pulse number')
    call xertst( 0.le.run, 'Invalid run number')
#else
    call xertst( 0.lt.shot.and.shot.le.214748, 'Invalid shot number')
    call xertst( 0.le.run.and.run.le.99999, 'Invalid run number')
#endif
    call xertst( .not.streql(username,' '), 'User name not defined !')
    call xertst( .not.streql(database,' '), 'Database not defined !')

    !! Parse HOME and SOLPSTOP and remove IMASDIR prefix from path if present
    !! and correct from input if necessary
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
      path = ids_path
    end if
    if (index(imasdir,trim(username)).eq.0) then
      l=index(imasdir,trim(run_user))
      m=index(imasdir(l+len_trim(run_user):256),'/')
      write(imasdir,'(a)') &
        & imasdir(1:l-1)//trim(username)//trim(imasdir(m+l+len_trim(run_user)-1:256))
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
      ids_path = trim(imasdir)//'/'//int2str(shot)//'/'//int2str(run)
#else
      ids_path = imasdir
#endif
    else if (.not.absolute_path) then
      ids_path = trim(imasdir)//'/'//trim(path)
    else
      ids_path = path
    end if

    !! If step was defined then run the b2mn_step routine
    if( new_run ) then
        write(0,*) "Running b2mn_step()"
        call b2mn_step( J )
        write(0,*) "b2mn_step() completed"
#ifdef B25_EIRENE
    else
      CALL EIRENE_ALLOC_COMUSR(1)
      call eirene_extrab25_eirpbls_init(nmol,nion,npls)
      call ntread
#endif
    end if

    !! Process B2.5 data and set it to IMAS IDS
    write(*,*) "START B25_process_ids"
#if AL_MAJOR_VERSION > 4
    uri = 'imas:mdsplus?path='//trim(ids_path)
    write(0,'(2a)') "Checking if IMAS data entry already exists : ",trim(ids_path)
    call imas_open( uri, OPEN_PULSE, idx, status, message )
#else
    write(0,'(2a,2i8)') "Checking if IDS already exists : ", trim(database), shot, run
    call imas_create_env( treename, shot, run, 0, 0, idx, username, &
       &     database, version, status )
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
        call imas_create_env( treename, shot, run, 0, 0, idx, username, &
       &     'iter', version, status )
#endif
        if (status.eq.0) then
          database = 'iter'
          write(0,*) "Old database case found."
          write(0,*) "Will be rewritten in new location."
        end if
      else if (database.eq.'iter'.or.index(ids_path,'imasdb/iter').gt.0) then
#if AL_MAJOR_VERSION > 4
        l=index(ids_path,'imasdb/ITER')
        write(olddir,'(a)') ids_path(1:l-1)//'imasdb/iter'//ids_path(l+11:256)
        uri = 'imas:mdsplus?path='//trim(olddir)
        call imas_open( uri, OPEN_PULSE, idx, status, message )
#else
        call imas_create_env( treename, shot, run, 0, 0, idx, username, &
       &     'ITER', version, status )
#endif
        database = 'ITER'
      end if
    end if
    !! If this is a time continuation run, append the new data to the IDS
    if ( status.eq.0 .and. idx.ne.0 ) then
      write (0,*) "Reading old IDS ", trim(database), shot, run
      call ids_get( idx, "equilibrium", equilibrium, status)
      call ids_get( idx, "edge_profiles", old_edge_profiles, status)
      if ( status.ne.0 ) then
        write (0,*) 'Error opening old edge_profiles IDS ! Will create a new one.'
        idx = 0
        continued = .false.
      else
        num_time_slices = size(edge_profiles%time)
        if (num_time_slices.gt.0) then
          ids_end_time = edge_profiles%time(num_time_slices)
        else
          ids_end_time = IDS_REAL_INVALID
        end if
        call ids_deallocate( old_edge_profiles )
        old_start_time = 0.0_IDS_real
        old_end_time = IDS_REAL_INVALID
        old_imas_version = 'x.xx.x'
#if ( IMAS_MAJOR_VERSION > 4 || ( IMAS_MAJOR_VERSION == 4 && IMAS_MINOR_VERSION > 1 ) )
        call ids_get( idx, "summary", old_summary, status )
        if ( status.ne.0 ) then
          write (0,*) 'Error opening old summary IDS !'
        else if (associated(old_summary%ids_properties% &
                        &   version_put%data_dictionary)) then
          old_start_time = old_summary%simulation%time_begin
          old_end_time = old_summary%simulation%time_end
          old_imas_version = old_summary%ids_properties% &
                          &  version_put%data_dictionary(1)
          call ids_deallocate( old_summary )
#else
        call ids_get( idx, "dataset_description", old_description, status )
        if ( status.ne.0 ) then
          write (0,*) 'Error opening old dataset_description IDS !'
        else if (associated(old_description%ids_properties% &
                        &   version_put%data_dictionary)) then
#if ( IMAS_MINOR_VERSION > 25 || IMAS_MAJOR_VERSION > 3 )
          old_start_time = old_description%simulation%time_begin
          old_end_time = old_description%simulation%time_end
#endif
          old_imas_version = old_description%ids_properties% &
                          &  version_put%data_dictionary(1)
          call ids_deallocate( old_description )
#endif
        else if ( streql(old_imas_version,'x.xx.x') ) then
          call xerrab ('Old IMAS data entry is incomplete !')
        end if
        continued = run_start_time.eq.IDS_REAL_INVALID .and. &
           &       (ids_end_time.lt.tim .and. ids_end_time.ne.IDS_REAL_INVALID)
        continued = continued .or. &
           &        run_start_time.ge.ids_end_time
        if (continued.or.database.eq.'iter'.or.index(ids_path,'imasdb/iter').gt.0) then
          if (.not.streql(old_imas_version,imas_version).or.database.eq.'iter'.or. &
            & index(ids_path,'imasdb/iter').gt.0) then
            if (.not.streql(old_imas_version,imas_version)) then
              write(*,*) &
               & 'Old IDS was written using IMAS version '// &
               &  trim(old_imas_version)//'.'
              write(*,*) &
               & 'Recreating using IMAS version '// &
               &  trim(imas_version)//'.'
            end if
            if (database.eq.'iter'.or.index(ids_path,'imasdb/iter').gt.0) &
               &  write(*,*) 'IDS file will be moved to ITER database.'
            call close_ual(idx)
!xpb Copy the IDS to a temporary location with the new DD and then bring it back
            tmp_run = run
            if (database.ne.'iter'.and.index(ids_path,'imasdb/iter').eq.0) then
              tmp_run = run + 1000
#if AL_MAJOR_VERSION > 4
              call al_build_uri_from_legacy_parameters               &
                & ( MDSPLUS_BACKEND, shot, run, trim(username),      &
                &   trim(database), int2str(IMAS_MAJOR_VERSION), '', &
                &   uri_source, status )
              call al_build_uri_from_legacy_parameters               &
                & ( MDSPLUS_BACKEND, shot, tmp_run, trim(username),  &
                &   trim(database), int2str(IMAS_MAJOR_VERSION), '', &
                &   uri_source, status )
#if IMAS_MAJOR_VERSION < 4
              write(systemarg,'(a,a,a,a,a)')                   &
                & 'idscp --set-dataset-version'//              &
                &      ' -s ',trim(uri_source),                &
                &      ' -d ',trim(uri_dest),' --dd-update'
#else
              write(systemarg,'(a,a,a,a,a)')                   &
                & 'idscp -s ',trim(uri_source),                &
                &      ' -d ',trim(uri_dest),' --dd-update'
#endif
#else
#if ( IMAS_MINOR_VERSION > 31 || IMAS_MAJOR_VERSION > 3 )
              write(systemarg,'(a,i7,a,i4,a,i7,a,i4,a,a,a,a)') &
                & 'idscp --setDatasetVersion'//                &
                &      ' -si ',shot,' -ri ',run,               &
                &      ' -so ',shot,' -ro ',tmp_run,           &
                &      ' -d ',trim(database),' -u ',trim(username)
#else
              write(systemarg,'(a,i7,a,i4,a,i7,a,i4,a,a,a,a)') &
                & 'idscp -si ',shot,' -ri ',run,               &
                &      ' -so ',shot,' -ro ',tmp_run,           &
                &      ' -d ',trim(database),' -u ',trim(username)
#endif
#endif
#ifdef NAGFOR
              call system(systemarg, status, ierror)
#else
              call system(systemarg)
#endif
            end if
#if AL_MAJOR_VERSION > 4
            call al_build_uri_from_legacy_parameters                 &
              & ( MDSPLUS_BACKEND, shot, tmp_run, trim(username),    &
              &   trim(database), int2str(IMAS_MAJOR_VERSION), '',   &
              &   uri_source, status )
            if (database.eq.'iter'.or.index(ids_path,'imasdb/iter').gt.0) then
              call al_build_uri_from_legacy_parameters               &
                & ( MDSPLUS_BACKEND, shot, run, trim(username),      &
                &   'ITER', int2str(IMAS_MAJOR_VERSION), '',         &
                &   uri_source, status )
            else
              call al_build_uri_from_legacy_parameters               &
                & ( MDSPLUS_BACKEND, shot, run, trim(username),      &
                &   trim(database), int2str(IMAS_MAJOR_VERSION), '', &
                &   uri_source, status )
            end if
            write(systemarg,'(a,a,a,a)')                     &
              & 'idscp -s ',trim(uri_source),                &
              &      ' -d ',trim(uri_dest)
#else
#if ( IMAS_MINOR_VERSION > 31 || IMAS_MAJOR_VERSION > 3 )
            write(systemarg,'(a,i7,a,i4,a,i7,a,i4,a,a,a,a)') &
             & 'idscp --setDatasetVersion'//                 &
             &       ' -si ',shot,' -ri ',tmp_run,           &
             &       ' -so ',shot,' -ro ',run,               &
             &       ' -d ',trim(database),' -u ',trim(username)
#else
            write(systemarg,'(a,i7,a,i4,a,i7,a,i4,a,a,a,a)') &
             & 'idscp -si ',shot,' -ri ',tmp_run,            &
             &      ' -so ',shot,' -ro ',run,                &
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
          call B25_process_ids( &
             &  edge_profiles, edge_sources, edge_transport, &
             &  radiation, wall, &
#if ( IMAS_MAJOR_VERSION < 4 || ( IMAS_MAJOR_VERSION == 4 && IMAS_MINOR_VERSION < 1 ) )
             &  description, &
#endif
             &  equilibrium, &
#if ( IMAS_MINOR_VERSION > 21 || IMAS_MAJOR_VERSION > 3 )
             &  summary, &
#endif
#if ( IMAS_MINOR_VERSION > 25 && IMAS_MINOR_VERSION < 34 && IMAS_MAJOR_VERSION == 3 )
             &  numerics,
#endif
#if ( IMAS_MINOR_VERSION > 25 || IMAS_MAJOR_VERSION > 3 )
             &  old_start_time, run_end_time, &
#endif
#if ( IMAS_MINOR_VERSION > 30 || IMAS_MAJOR_VERSION > 3 )
             &  divertors, &
#endif
#if IMAS_MAJOR_VERSION > 3
             &  tim, dtim, shot, database, new_eq_ggd, &
#else
             &  tim, dtim, shot, run, database, version, new_eq_ggd, &
#endif
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
      if (index(ids_path,'imasdb/iter').gt.0) then
        l=index(ids_path,'imasdb/iter')
        write(ids_path(l+7:l+10),'(a4)') 'ITER'
      end if
    end if
    if ( status.ne.0 .or. idx.eq.0 ) then
      call B25_process_ids( &
         &  edge_profiles, edge_sources, edge_transport, &
         &  radiation, wall, &
#if ( IMAS_MAJOR_VERSION < 4 || ( IMAS_MAJOR_VERSION == 4 && IMAS_MINOR_VERSION < 1 ) )
         &  description, &
#endif
         &  equilibrium, &
#if ( IMAS_MINOR_VERSION > 21 || IMAS_MAJOR_VERSION > 3 )
         &  summary, &
#endif
#if ( IMAS_MINOR_VERSION > 25 && IMAS_MINOR_VERSION < 34 && IMAS_MAJOR_VERSION == 3 )
         &  numerics,
#endif
#if ( IMAS_MINOR_VERSION > 25 || IMAS_MAJOR_VERSION > 3 )
         &  run_start_time, run_end_time, &
#endif
#if ( IMAS_MINOR_VERSION > 30 || IMAS_MAJOR_VERSION > 3 )
         &  divertors, &
#endif
#if IMAS_MAJOR_VERSION > 3
         &  tim, dtim, shot, database, new_eq_ggd )
#else
         &  tim, dtim, shot, run, database, version, new_eq_ggd )
#endif
    end if

    !! Create Write the set data to IDSs
    write(*,*) "START put_ids_edge"
    call put_ids_edge( &
        &   edge_profiles, edge_sources, edge_transport, &
        &   radiation, wall, &
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
    call dealloc_ids_edge( edge_profiles, edge_sources, edge_transport, &
#if ( IMAS_MINOR_VERSION > 25 && IMAS_MINOR_VERSION < 34 && IMAS_MAJOR_VERSION == 3 )
        &   numerics, &
#endif
#if ( IMAS_MINOR_VERSION > 30 || IMAS_MAJOR_VERSION > 3 )
        &   divertors, &
#endif
        &   radiation, wall )
    call dealloc_batch_edge( batch_profiles, batch_sources, &
#if ( IMAS_MINOR_VERSION > 21 || IMAS_MAJOR_VERSION > 3 )
        &   summary, &
#endif
#if ( IMAS_MAJOR_VERSION < 4 || ( IMAS_MAJOR_VERSION == 4 && IMAS_MINOR_VERSION < 1 ) )
        &   description &
#endif
        &   )
    call close_ual(idx)

    write(0,*) " Running b2mn_fin"
    call b2mn_fin
    write(0,*) "b2mn_fin completed"

contains

    !> Subroutine intended to check if supposed new file already exists. If
    !! the file exists it deletes it.
    subroutine checkFileAndDelete( fileName )
        character(len=*), intent(in) :: fileName    !< Name of the file to be
                                                    !< checked
        !! Internal variables
        integer :: stat
        logical :: file_exists

        inquire( file=fileName, exist=file_exists )
        if ( file_exists ) then
            write(*,*) "Deleting old ", fileName
            open(unit=1234, iostat=stat, file=fileName, status='old')
            if (stat == 0) close(1234, status='delete')
        endif

    end subroutine

end program b2_ual_write_b2mod

!!!Local Variables:
!!! mode: f90
!!! End:
