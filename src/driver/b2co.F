*-----------------------------------------------------------------------
*.specification

      program b2co
      use b2mod_types
      use b2mod_version
      use b2mod_time
!      use b2mod_diag
!      use b2mod_b2plot
!      use b2mod_external
!      use b2mod_elements
      use b2mod_constants
      use b2mod_user_namelist
      use b2mod_boundary_namelist
      use b2mod_neutrals_namelist
      use b2mod_neutr_src_scaling
      use b2us_geo
      use b2us_map
      use b2us_plasma
      use b2us_io
#ifdef NAGF95
      use f90_unix           ! IGNORE
#endif
      use b2mod_b2cmfs
      use b2mod_b2cmpa
      use b2mod_file
      use b2mod_xerset
      use b2mod_ipmain
#ifdef B25_EIRENE
      use eirmod_parmmod
     .   , ncv_eir => ncv
      use eirmod_clogau
      use eirmod_ctrcei
      use eirmod_braeir
      use eirmod_cinit
     , , only: eirene_init_cinit, fort_lc
#ifdef WG_TODO
      use eirmod_comprt
     , , only: iunout
#endif
      use eirmod_eirbra
     & , smo_eir => smo
#else
      use b2mod_braeir
      use b2mod_eirbra
     & , smo_eir => smo
#endif
      implicit none
#include "DIMENSIONS.F"

*-----------------------------------------------------------------------
*.documentation
*
*  1. purpose
*
*     B2CO reads in a formatted plasmastate file, and writes out a
*     a b2fstate file.
*
*
*-----------------------------------------------------------------------
*.declarations

*   ..common blocks

*   ..main data structures
      type(mapping), save    :: mpg
      type(geometry), save   :: geo
      type(B2State), save    :: state
!      type(B2StateExt), save :: state_ext
!      type(switches), save   :: switch

*   ..local variables
#ifndef NO_GETENV
      integer lenval, ierror
#ifndef USE_PXFGETENV
      intrinsic get_environment_variable
#endif
#endif
      integer ninp(0:5), nout(0:7), nCv, nFc, ns, ns0, idum(0:9), io
#ifndef NAGF95
#ifndef NAGFOR
      integer iargc
#endif
#endif
#ifdef WG_TODO
      integer iflag
      integer :: b2time = 1, ntim_batch = 500
      integer :: ismain = 1, b2wall_netcdf = 0, mvnum = 0, tally = 1
      integer :: use_eirene = 0, external_species = 0
      integer :: boundary_namelist, neutrals_namelist = 0
      real (kind=R8) :: BoRiS = 0.0_R8
      real (kind=R8) :: dtim = 1.0_R8
#endif
      character savefile*256, arg*256, lblgm*120, lblmn*120, lblcp*120
#ifdef WG_TODO
      character filename*256
      character b2fstate_only*256
#ifdef B25_EIRENE
      character extension*4
#endif
      logical file_exists, todo
      integer lnblnk
      external lnblnk
#endif
      character cnamip*80, cvalip*80, line*500
      logical file_open
*   ..procedures
      logical is_comment, streql
      external is_comment, streql
      external prgini, prgend, xertst, cfopen, cfruin
      external b2rucp, b2rups, b2wuzd, b2wups, b2xvps, ipgeti, ipgetr
*   ..initialise input/output unit numbers
      data ninp/50,51,52,53,54,55/
      data nout/60,61,62,63,64,65,66,67/
      data savefile/'plasmastate.xxxx'/

*-----------------------------------------------------------------------
*.computation

*   ..program start-up calls
      call prgini ('b2co')
*   ..open files
      call cfopen (ninp(0),'b2mn.dat','old','formatted')
      call cfopen (ninp(1),'b2fgmtry','old','un*formatted')
      call cfopen (ninp(2),'b2fpardf','old','un*formatted')
      call cfopen (ninp(3),'b2fstati','old','un*formatted')
#ifndef NAGFOR
      call xertst (iargc().ge.1,
     &  'program should be supplied with a state number')
#ifdef USE_PXFGETENV
      call pxfgetarg(1,arg,lenval,ierror)
#else
      call getarg(1,arg)
#endif
#else
      call xertst (command_argument_count().ge.1,
     &  'program should be supplied with a state number')
      call get_command_argument(1,arg,lenval,ierror)
#endif
      savefile=arg
      call cfopen (ninp(4), trim(savefile), 'old','unformatted')
      call cfopen (nout(0),'b2co.prt','new','formatted')
      call cfopen (nout(1),'b2co.files','new','formatted')
      call cfopen (nout(2),'b2fstate','new','un*formatted')
c      call cfopen (nout(7),'b2fplasma','new','unformatted')
*   ..mark output unit for error messages
      call xerset (0)
*   ..obtain version numbers                                           !xpb
      call cfverr (ninp(1), b2fgmtry_version)                          !xpb
      call cfverr (ninp(2), b2fparam_version)                          !xpb
      call cfverr (ninp(3), b2fstati_version)                          !xpb
*   ..obtain mapping/geometry dimensions from b2fgmtry file
      call cfruin (ninp(1),7,idum,'nCi,nCg,nCv,nFc,nVx,nFs,nFt')
      mpg%nCi = idum(0) 
      mpg%nCg = idum(1) 
      mpg%nCv = idum(2)
      mpg%nFc = idum(3)
      mpg%nVx = idum(4)
      mpg%nFs = idum(5)
      mpg%nFt = idum(6)
      call cfruin (ninp(1),5,idum,'nCmxVx,nCmxFc,nVmxCv,nVmxFc,nCmxNv')
      mpg%nCmxVx = idum(0)
      mpg%nCmxFc = idum(1)
      mpg%nVmxCv = idum(2)
      mpg%nVmxFc = idum(3)
      mpg%nCmxNv = idum(4)
      call cfruin (ninp(1),1,idum,'isClassicalGrid')
      mpg%isClassicalGrid = idum(0)
      call cfruin (ninp(1),3,idum,'nx,ny,nncut')
      mpg%nx    = idum(0)
      mpg%ny    = idum(1)
      mpg%nncut = idum(2)
      call xertst (0.le.mpg%nCv.and.0.le.mpg%nFc,
     &  'faulty input nCv, nFc from geometry file')
*   ..obtain ns from parameters file
      call cfruin (ninp(2), 1, idum, 'ns')
      ns = idum(0)
      call xertst (1.le.ns, 'faulty input ns from parameters file')
*   ..obtain nCv, nFc, ns from b2fstate file
      call cfruin (ninp(3), 3, idum, 'nCv,nFc,ns')
      nCv = idum(0)
      nFc = idum(1)
      ns0 = idum(2)
      call xertst (mpg%nCv.eq.nCv.and.mpg%nFc.eq.nFc,
     &  'faulty input nCv, nFc, ns from b2fstate file')
      call xertst (1.le.ns0.and.ns0.le.ns,
     &  'faulty input ns from plasma state file')
*   ..call principal driver routine
      call alloc_mapping (mpg)
      call alloc_geometry (geo, mpg%nCv, mpg%nFc, mpg%nVx, mpg%nFs, 
     .  mpg%nFt, mpg%nVmxCv)
      call createB2State (nCv, nFc, ns, state)
!      call createB2StateExt ( mpg%nCv, mpg%nFc, ns, state_ext )
*   ..obtain lblmn
      call b2xxid ('b2mn', lblmn(1:60))
      call b2rflb (ninp(0), lblmn(61:120))
      label = lblmn
*   ..read physics common
      call cfruch (ninp(2), 120, lblcp, 'label')
      call b2rucp (ninp(2), b2fparam_version, ns)
*   ..read revised parameters in physics common
*     (ninp(0) need supply only parameters that are to be changed from
*     the values read in as default physics common.)
      call b2rfcp (ninp(0), ns)
*   ..test parameters
      call b2xvcp (ns)
*   ..read and echo code internal parameters
    1 continue
       read (ninp(0),'(a)',end=2,err=91) line
       if (.not.is_comment(line)) then
         read (line,*) cnamip, cvalip
         call ipsetc (cnamip, cvalip)
       end if
       goto 1
    2 continue
*   ..read header for geometry
      call cfruch (ninp(1), 120, lblgm, 'label')
      call cfruin (ninp(1), 1, idum, 'isymm')
      isymm = idum(0)
*     (isymm is now stored in the b2cmfs common)               ! xpb 10.08.00
*   ..read geometry and mapping
      call read_b2fgmtry (ninp(1), mpg, geo)
!WG_TODO      call get_nxt(nx, nxtl, nxtr)
*   ..read initial state
      call cfverr (ninp(4), plasmastate_version)                        !xpb
      call read_b2fstate (ninp(4), nCv, nFc, ns, state)
*   ..obtain external species data if any
!WG_TODO      call ipgeti ('b2mndr_external_species', external_species)
!WG_TODO      if (external_species.ne.0) then
!WG_TODO        call read_b2mod_external(nx,ny,ns,ismain)
!WG_TODO      else
!WG_TODO        call alloc_b2mod_external_dummy(nx,ny)
!WG_TODO      endif
*   ..test plasma state and fluxes
      call b2xvps (nCv, nFc, ns, state % pl, state % dv)
*   ..save final state
      newversion = plasmastate_version
      write(lblmn,*) 'Recovered from ',trim(savefile)
      idum(0) = nCv
      idum(1) = nFc
      idum(2) = ns
      call cfverw (nout(2), plasmastate_version)                        !xpb
      call cfwuin (nout(2), 3, idum, 'nCv,nFc,ns')
      call cfwuch (nout(2), 120, lblmn, 'label')
      call b2wuzd (nout(2), plasmastate_version, ns,
     &  zamin, zamax, zn, am)
      call write_b2fstate (nout(2), mpg%nCv, mpg%nFc, ns, state)
*   ..close graphics files
*     (none)
*   ..write plasma state
c      call cfverw(nout(7), plasmastate_version)
c      idum(0) = nCv
c      idum(1) = nFc
c      idum(2) = ns
c      call cfwuin (nout(7), 3, idum, 'nCv,nFc,ns')
c      call write_b2fplasma (nout(7), nCv, nFc, ns, state)
      write(nout(0),*) 'Converted ',trim(savefile),' to b2fstate'
*   ..truncate tracing files
#ifdef WG_TODO
      call ipgeti ('ank_tracing', ank_tracing)
      call ipgeti ('b2mndr_ismain', ismain)
      call ipgetr ('b2mndr_dtim', dtim)
      call ipgetr ('b2news_BoRiS', BoRiS)
      call ipgeti ('b2mndr_eirene', use_eirene)
      call ipgeti ('b2stbr_sput_dst', sput_dst)
      call ipgeti ('b2stbc_boundary_namelist', boundary_namelist)
      call ipgeti ('b2stbr_neutrals_namelist', neutrals_namelist)
#ifdef B25_EIRENE
      filename='input.dat'
      call find_file(filename,file_exists)
      if (file_exists.and.use_eirene.ne.0) then
        open(unit=45,file=filename)
        read(45,'(A72)',END=6999) line
        rewind(45)
        CALL EIRENE_FIND_PARAM(45)
        close(45)
        CALL EIRENE_SET_PARMMOD(1)
        write(*,*) 'Called FIND_PARAM'
        write(*,*) '   NDX = ',NDX
        write(*,*) '   NDY = ',NDY
        write(*,*) ' NSTRA = ',NSTRA
        write(*,*) '   NFL = ',NFL
        write(*,*) '  NATM = ',NATM
        write(*,*) '  NMOL = ',NMOL
        write(*,*) '  NION = ',NION
        write(*,*) 'NLIMPS = ',NLIMPS
      endif
      goto 7000
 6999 continue
      file_exists = .false.
 7000 continue
      if (.not.file_exists.or.use_eirene.eq.0) then
        ndx=DEF_NXD+DEF_NCUT*5
        ndy=DEF_NYD
        nfl=DEF_NFL
        natm=DEF_NATM
        nmol=DEF_NMOL
        nion=DEF_NION
        npls=DEF_NPLS
        nlim=DEF_NLIM
        nsts=DEF_NSTS
        nstra=DEF_NSTRA
        nsrfs=DEF_NSRFS
        N1ST=DEF_NYD+1
        N2ND=DEF_NXD+MAX(1,DEF_NCUT/2*3)+
     .       DEF_ISOEXTRA*MAX(0,DEF_NCUT/2-1)
        N3RD=DEF_NPHID
        NREAC=7
        IUNOUT=6
        CALL EIRENE_SET_PARMMOD(1)
        CALL EIRENE_ALLOC_CLOGAU
        CALL EIRENE_INIT_CINIT
      endif
#endif
      natmi = natm
      nmoli = nmol
      nioni = nion
      NFLAI = NFL
#ifdef B25_EIRENE
      CALL EIRENE_ALLOC_BRAEIR(NDX,NDY,NFL,ifoff)
      CALL EIRENE_ALLOC_EIRBRA(NDX,NDY,NFL,NSTRA,ifoff)
      CALL EIRENE_ALLOC_CTRCEI
      CALL EIRENE_SETCON
      CALL EIRENE_EXTRAB25_WNEUINIT
#else
      call alloc_b2mod_braeir(nx,ny,ns,nstra)
!PB      call alloc_b2mod_eirbra(nx,ny,ns,natm,npls,nstra)
      call alloc_b2mod_eirbra(nx*ny,ns,natm,npls,nstra)
#endif
      call alloc_b2mod_elements(ns,zn,am)
      call alloc_b2mod_b2plot(nx,ny,ns,ns_ext,wklng,
     1 natm,nmol,nion,nlimps)
      call alloc_b2mod_boundary(nx, ny, ns, nbcd, nnisod)
      if (boundary_namelist.eq.1)
     > call read_b2mod_boundary_namelist(nx, ny, ns, .false.)
      call alloc_b2mod_neutrals
     . (nx, ny, ns, natm, nmol, nion, nstraid, nsts, sput_dst)
      if (neutrals_namelist.eq.1 .or. use_eirene.ne.0)
     > call read_neutrals_namelist(nx, ny, ns, .false.)
!!WG_TODO      call alloc_b2mod_diag(nx,ny,ns,nnreg(0),nstraid,natm)
      write_nml_user = .false.
      call alloc_b2mod_wall_init(nx,ny,nncut)
#ifdef WG_TODO
      call init_wall(nx,ny)
#endif
      first_time_step = .true.  ! IYS 21.06.2022
      if (use_eirene.eq.0) then
        nnatmi = 0
        nnmoli = 0
      else
#ifdef B25_EIRENE
        extension=savefile(lnblnk(savefile)-3:lnblnk(savefile))
        filename=fort_lc//'44.'//extension
        call find_file(filename,file_exists)
#else
        file_exists = .false.
#endif
        if (file_exists) then
          open(unit=44,file=filename)
          call ngread(ismain,iflag)
          nnatmi = natmi
          nnmoli = nmoli
        else
          nnatmi = 0
          nnmoli = 0
        end if
      end if
      if (ank_tracing.ne.0) then
        call b2trci(ns,nxtl,nxtr,geo,mpg,switch)
        call b2trcc(nout(1),nxtl,nxtr,geo,mpg,switch)
        call b2filei(0,ns,nnreg(0))
        call b2filec(0,ns,nnreg(0))
      end if
*   ..rewrite NetCDF files
#ifdef NO_GETENV
      todo=.true.
#else
#ifdef USE_PXFGETENV
      CALL PXFGETENV ('B2FSTATE_ONLY', 0, b2fstate_only, lenval, ierror)
      if (ierror.ne.0) b2fstate_only='no'
#else
      call get_environment_variable('B2FSTATE_ONLY',
     . status=ierror, length=lenval)
      if (ierror.eq.0) then
        call get_environment_variable('B2FSTATE_ONLY',
     .                           value=b2fstate_only)
      elseif (ierror.eq.1) then
        b2fstate_only='no'
      endif
#endif
      todo=.not.streql(b2fstate_only,'yes')
#endif
      call ipgeti ('b2mndr_b2time', b2time)
      call ipgeti ('b2mndt_av_ntim_batch', ntim_batch)
      if (todo.and.b2time.gt.0) then
        filename='b2time.nc'
        call find_file(filename,file_exists)
        if(file_exists) then
          write(*,*) 'Truncating ', trim(filename)
          call b2co_time_cdf(nx, ny, ns, tim)
          write(*,*) 'Returned from b2co_time_cdf'
        else
          write(*,*) 'Skipped calling b2co_time_cdf because '//
     .               'no b2time.nc file found'
        end if
      end if
      if (todo.and.ntim_batch.gt.0) then
        filename='b2batch.nc'
        call find_file(filename,file_exists)
        if(file_exists) then
          write(*,*) 'Truncating ', trim(filename)
          call b2co_time_cdf(nx, ny, ns, tim)
          write(*,*) 'Returned from b2co_time_cdf'
        else
          write(*,*) 'Skipped calling b2co_time_cdf because '//
     .               'no b2batch.nc file found'
        end if
      end if
      call ipgeti ('b2stbr_b2wall_netcdf', b2wall_netcdf)
      if (todo.and.b2wall_netcdf.gt.0) then
        filename='b2wall.nc'
        call find_file(filename,file_exists)
        if(file_exists) then
          write(*,*) 'Truncating ', trim(filename)
          call b2co_wall_cdf(nx, ny, ns, tim)
          write(*,*) 'Returned from b2co_wall_cdf'
        else
          write(*,*) 'Skipped calling b2co_wall_cdf because '//
     .               'no b2wall.nc file found'
        end if
      end if
      call ipgeti ('b2mndr_tally', tally)
      if (todo.and.tally.gt.0) then
        filename='b2tallies.nc'
        call find_file(filename,file_exists)
        if(file_exists) then
          write(*,*) 'Truncating ', trim(filename)
          call b2co_tally_cdf(nx, ny, ns, tim)
          write(*,*) 'Returned from b2co_tally_cdf'
        else
          write(*,*) 'Skipped calling b2co_tally_cdf because '//
     .               'no b2tallies.nc file found'
        end if
      end if
      call ipgeti ('b2mndr_mvnum', mvnum)
      if (todo.and.mvnum.gt.0) then
        filename='b2movies.nc'
        call find_file(filename,file_exists)
        if(file_exists) then
          write(*,*) 'Truncating ', trim(filename)
          call b2co_movies_cdf(nx, ny, nx, nnatmi, nnmoli, tim)
          write(*,*) 'Returned from b2co_movies_cdf'
        else
          write(*,*) 'Skipped calling b2co_movies_cdf because '//
     .               'no b2movies.nc file found'
        end if
      end if
#endif
*   ..close other files
      do io=0,size(ninp)-1
         inquire(unit=ninp(io),opened=file_open)
         if (file_open) close(ninp(io))
      enddo
      do io=0,size(nout)-1
         inquire(unit=nout(io),opened=file_open)
         if (file_open) close(nout(io))
      enddo
* ..end of computation
      call destroyB2State (state)
!      call destroyB2StateExt (state_ext)
      call dealloc_geometry (geo)
      call dealloc_mapping (mpg)
      call prgend ()
      stop 'b2co'

*-----------------------------------------------------------------------
*.error conditions

   91 call xerrab ('error trying to read internal parameters')

*-----------------------------------------------------------------------
*.end b2co

      end program b2co

