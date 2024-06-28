      subroutine b2mds(mpg,gm,state,state_ext,nsd_local,switch)
      use b2mod_types
      use b2mod_constants
      use b2mod_b2plot
     * , only : textpl, textwl, textsp, textbl,
     *          textcomponents, text0123,
     *          natmi, nmoli, nioni, nlimi, wklng,
     *          solpstop, alloc_b2mod_b2plot, dealloc_b2mod_b2plot
      use b2mod_wall
      use b2mod_work
      use b2mod_eirene_globals
      use b2mod_user_namelist
      use b2mod_boundary_namelist
      use b2mod_neutrals_namelist
      use b2mod_transport_namelist
      use b2mod_elements
      use b2mod_sputter
#ifdef NAGF95
      use f90_unix           ! IGNORE
#endif
#ifdef MDSPLUS
      use b2mod_mdsplus
#endif
      use b2mod_version
      use b2us_map
      use b2us_geo
      use b2us_plasma
      use b2us_prep
      use b2mod_switches
#ifdef B25_EIRENE
      use eirmod_parmmod
     * , only : ndx, ndy, nstra, nfl, natm, nmol, nion, npls,
     *          nsrfs, n1st, n2nd, n3rd, nreac, 
     *          eirene_set_parmmod
      use eirmod_cinit
     * , only : eirene_init_cinit, fort_lc
      use eirmod_clogau
      use eirmod_ctrcei
      use eirmod_comprt
      use eirmod_braeir
      use eirmod_eirbra
     * , only : eirene_alloc_eirbra, eirene_dealloc_eirbra
      use eirmod_infcop
      use eirmod_wneutrals
#else
      use b2mod_braeir
      use b2mod_eirbra
     & , smo_eir => smo
#endif
      use b2mod_b2cmpa
      use b2mod_b2cmfs
      use b2mod_external
      use b2mod_subsys
      use b2mod_dimensions
      implicit none

#ifdef MDSPLUS
#     include "mdslib.inc"
#endif
      integer nsd_local
      character*24 hlp_frm
      integer iflag
      character*256 filename
#ifdef B25_EIRENE
      character*72 zeile
#endif
      type(mapping) :: mpg
      type(geometry) :: gm
      type(switches) :: switch
      type(B2State) :: state
      type(B2StateExt) :: state_ext
c
c  version : 27.07.96 20:37
c
      integer ispecies(0:nsd_local-1)

      integer is,iss,isp,nfla,ns
c
      logical ex
c
      integer il,imi,iFc
      real (kind=R8) ::  pr(mpg%nCv)
#ifdef MDSPLUS
      character*512 textpl_cat, textwl_cat
      character*12 :: textdrift(0:1)
      integer :: status
      character :: mds_signal_name*(64), tdi*(172)
      real (kind=R8), allocatable ::
     *   vesl(:,:)
      integer i, j, iCv, icnt, nsw, iVx1, iVx2
#ifdef WG_TODO
      integer itmp(-1:nxd,-1:nyd)
#endif
      real (kind=R8) ::
     *   rdummy, rtmp, temp(mpg%nCv),
     5   tmp(mpg%nCv,0:nsd_local-1),
     5   crf(mpg%nFc), czf(mpg%nFc),
     7   flux_x(mpg%nFc,0:nsd_local-1),
     7   flux_y(mpg%nFc,0:nsd_local-1),
#ifdef WG_TODO
     8   dsrad(-1:nxd,-1:nyd),dspol(-1:nxd,-1:nyd),dspar(-1:nxd,-1:nyd),
#endif
     9   fhesep,fhisep,fnesep,fnisep,
     *   pot(0:nsd_local-1),poti(0:nsd_local-1),
     7   za(0:nsd_local-1)
#endif
c
      integer nCv, nFc, nVx, iCv1, iCv2
      integer ismain, refl_model, sput_phys_model, use_eirene,
     . sput_phys_col, boundary_namelist, neutrals_namelist
      data ismain/1/, refl_model/1/, sput_phys_col/3/, use_eirene/0/,
     3     sput_phys_model/1/, boundary_namelist/0/,
     4     neutrals_namelist/0/
      save ismain, refl_model, sput_phys_model,
     2     sput_phys_col, boundary_namelist, neutrals_namelist,
     3     use_eirene

#ifndef NO_GETENV
      integer lenval, ierror
#ifndef USE_PXFGETENV
      intrinsic get_environment_variable
#endif
#endif
*!!   (We avoid to declare, as external or intrinsic, the subroutine
*     getenv; it is external on some systems (SunOS, Unicos, ...) and
*     intrinsic on others (AIX, HP-UX, ...).  If getenv is declared
*     external on AIX or HP-UX it returns some non-printing garbage,
*     without error message (AIX/XLF Version 2.2, HP-UX Release 9.0))
*     ------------------------------------------------------------------
      external ipgeti, xertst

c Function to check the correct parameters!
c
      call subini('b2mds')
      nCv = mpg%nCv
      nFc = mpg%nFc
      nVx = mpg%nVx
      call ipgeti ('b2mndr_eirene', use_eirene)
#ifdef B25_EIRENE
      CALL EIRENE_DEFAULTS_USR
      filename='input.dat'
      call find_file(filename,ex)
      if (ex.and.use_eirene.ne.0) then
        open(unit=45,file=filename)
        read(45,'(A72)',END=6999) ZEILE
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
      ex = .false.
 7000 continue
      if (.not.ex.or.use_eirene.eq.0) then
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
#else
      call xertst(use_eirene.eq.0,
     &    'Compile code with B25_EIRENE option for coupled runs')
#endif
      NFLAI = NFL
#ifdef B25_EIRENE
      CALL EIRENE_ALLOC_BRAEIR(nCv,nFc,NFL)
      CALL EIRENE_ALLOC_EIRBRA(nCv,NFL,NATM,NMOL,NION,NSTRA)
      CALL EIRENE_ALLOC_CTRCEI(1)
      CALL EIRENE_ALLOC_CTRCEI(2)
      CALL EIRENE_SETCON
      CALL EIRENE_WNEUTRALS_INIT
#else
      call alloc_b2mod_braeir(nCv,nFc,nsd_local,nstra)
      call alloc_b2mod_eirbra(nCv,nsd_local,natm,nmol,nion,npls,nstra)
#endif
      call alloc_b2mod_b2plot(nCv,nFc,nVx,nsd_local,ns_ext,wklng,
     1 natm,nmol,nion,nlimps)
      solpstop=' '
#ifndef NO_GETENV
#ifdef USE_PXFGETENV
      CALL PXFGETENV ('SOLPSTOP', 0, solpstop, lenval, ierror)
#else
      call get_environment_variable('SOLPSTOP',
     .   status=ierror, length=lenval)
      if (ierror.eq.0) call get_environment_variable('SOLPSTOP',
     .   value=solpstop)
#endif
#endif
#ifdef MDSPLUS
      do iFc = 1, nFc
        iVx1 = mpg%fcVx(iFc,1)
        iVx2 = mpg%fcVx(iFc,2)
        crf(iFc) = (gm%vxX(iVx1)+gm%vxX(iVx2))/2.0_R8
        czf(iFc) = (gm%vxY(iVx1)+gm%vxY(iVx2))/2.0_R8
      end do
#endif
c
      ns=max(nsd_local,1)
      write(hlp_frm,'(a6,i3,a3)') '(1x,a,',ns,'i3)'
      write(*,hlp_frm) 'zn: ',(nint(zn(is)),is=0,ns-1)
      do is=0,ns-1
        call species(is,textpl(is),.false.)
        call species(is,textwl(is),.false.)
      end do
      iss=0
      il=0
      do is=0,ns-1
        if(nint(zn(is)) .ne. il) then
          il=nint(zn(is))
          if(il .eq. 0) then
            textsp(iss)=' '//'n'
          else if(il .eq. 1) then
            imi=nint(am(is))
            if(imi .eq. 1) textsp(iss)=' '//'H'
            if(imi .eq. 2) textsp(iss)=' '//'D'
            if(imi .eq. 3) textsp(iss)=' '//'T'
          else
            textsp(iss)=' '//elements(il)
          end if
          iss=iss+1
        end if
      end do
      call ipgeti ('b2mndr_ismain', ismain)
      call xertst (0.le.ismain.and.ismain.lt.ns,
     &         'faulty b2mndr_ismain')
      call xertst (.not.is_neutral(ismain),
     &         'ismain cannot refer to a neutral species')
      textbl(0)=' '
      textcomponents(0)='x (poloidal)'
      textcomponents(1)='y'
      textcomponents(2)='z (toroidal)'
      textcomponents(3)='(total)'
      text0123(0)='0'
      text0123(1)='1'
      text0123(2)='2'
      text0123(3)='3'

cxpb  This call is moved up from 'f31' so that we use the original cell-faced
cxpb  quantities before their mapping to cell centers
      call b2xpni (nCv, ns, state%pl%na, state%dv%ni)
      call b2xppr (nCv, ns,
     .  state%dv%ne, state%pl%na, state%pl%te, state%pl%ti, state%pl%tn,
     .  pr, state_ext)
      nfla=0
      do is=0,ns-1
        if(.not.is_neutral(is)) nfla=nfla+1
      enddo
      if(nfla.gt.NFL) then
        write(*,*) 'copy_background not called because nfla.gt.NFL',
     1   nfla,NFL
      else
        call copy_background(nCv, nFc, ns, switch, mpg, gm,
     &                       state%pl, state%dv, state_ext)
      endif
c
      call initialise_potentials(ns)
      isp=0
      ispecies(0)=isp
      do is=1,ns-1
        if(nint(zn(is)).ne.nint(zn(is-1)) .or.
     &     nint(am(is)).ne.nint(am(is-1))) then
          isp=isp+1
        endif
        ispecies(is)=isp
      enddo
      do is=0,ns-1
        write(*,'(i3,1x,i3,1x,f5.1,1x,f5.1,1x,f5.1,4f11.3)')
     1   is,ispecies(is),zamin(is),zamax(is),zn(is),
     2   potimin(is),potimax(is),potmin(is),potmax(is)
      end do
      call alloc_b2mod_elements(ns, zn, am)
c calculate additional fluxes (fhj, fhm, fhp, fht)
      call b2xpen (nCv, nFc, ns, mpg,
     1  state%dv%fna, state%dv%fna_53, state%dv%fne_53, state%dv%fch,
     1  state%dv%fhe, state%dv%fhi, state%dv%fhn, state%rt%rpt,
     2  state%pl%te, state%pl%ti, state%pl%tn, state%pl%po,
     2  switch%BoRiS,
     2  state%dv%fhm, state%dv%fhp, state%dv%fhj, state%dv%fht)
      state%dv%fht=state%dv%fht+state_ext%fhi
      do is=0,state_ext%ns-1
        do iFc=1,nFc
          iCv1 = mpg%fcCv(iFc,1)
          iCv2 = mpg%fcCv(iFc,2)
          state%dv%fht(iFc,0)=state%dv%fht(iFc,0)+
     .      0.5_R8*state_ext%fa(iFc,0,is)*
     .     (1.0_R8-switch%BoRis)*(0.5_R8*state_ext%am(is)*mp*
     .     (state_ext%ua(iCv1,is)**2+state_ext%ua(iCv2,is)**2)+
     .     (state_ext%ta(iCv1,is)   +state_ext%ta(iCv2,is)))
          state%dv%fht(iFc,0)=state%dv%fht(iFc,0)+
     .      0.5_R8*state_ext%fa(iFc,0,is)*
     .     (state_ext%pt(iCv1,is)+state_ext%pt(iCv2,is))*ev
          state%dv%fht(iFc,1)=state%dv%fht(iFc,1)+
     .      0.5_R8*state_ext%fa(iFc,1,is)*
     .     (1.0_R8-switch%BoRis)*(0.5_R8*state_ext%am(is)*mp*
     .     (state_ext%ua(iCv1,is)**2+state_ext%ua(iCv2,is)**2)+
     .     (state_ext%ta(iCv1,is)   +state_ext%ta(iCv2,is)))
          state%dv%fht(iFc,1)=state%dv%fht(iFc,1)+
     .      0.5_R8*state_ext%fa(iFc,1,is)*
     .     (state_ext%pt(iCv1,is)+state_ext%pt(iCv2,is))*ev
        enddo
      enddo
c
cxpb  Bring in wall data and compute oft-needed wall quantities
      call ipgeti ('b2stbr_refl_model', refl_model)
      call ipgeti ('b2stbr_plate_model', plate_model)
      if (plate_model.eq.0) then
        write(*,*)
     .   'Using O-D time-independent plate temperature model'
      else if (plate_model.eq.1) then
        write(*,*) 'Using 1-D time-dependent plate temperature model'
      else if (plate_model.eq.2) then
        write(*,*) 'Using 2-D time-dependent plate temperature model'
      else
        call xerrab ('Invalid value for b2stbr_plate_model')
      endif
      call ipgeti ('b2stbr_plate_option', plate_option)
      call xertst (1.le.plate_option .and.
     &             3.ge.plate_option, 'faulty b2stbr_plate_option')
      if (plate_option.eq.2) plate_option=3 ! we do not want to recompute the values!
      call ipgetr ('b2stbr_plate_temp', plate_temp)
      call xertst (0.0_R8.lt.plate_temp, 'faulty b2stbr_plate_temp')
      call ipgetr ('b2stbr_plate_thick', plate_thick)
      call xertst (0.0_R8.le.plate_thick, 'faulty b2stbr_plate_thick')
      call ipgeti ('b2stbr_sput_dst', sput_dst)
      call ipgeti ('b2stbr_sput_dst2', sput_dst2)
      call ipgeti ('b2stbr_sput_dst3', sput_dst3)
      call xertst (sput_dst.lt.ns, 'faulty b2stbr_sput_dst')
      call xertst (sput_dst2.lt.ns, 'faulty b2stbr_sput_dst2')
      call xertst (sput_dst3.lt.ns, 'faulty b2stbr_sput_dst3')
      call ipgeti ('b2stbr_sput_phys_col', sput_phys_col)
c the TRIM sputtering files have 9 columns corresponding to
c 0 15 30 45 55 65 75 80 85 degree impact
      call xertst (0.le.sput_phys_col.and.sput_phys_col.le.9,
     1 'faulty b2stbr_sput_phys_col')
      call ipgeti ('b2stbr_sput_phys_model', sput_phys_model)
      call alloc_b2mod_boundary(1, 1, ns, nbcd, nnisod)
      call ipgeti ('b2stbc_boundary_namelist', boundary_namelist)
      if (boundary_namelist.eq.1)
     > call read_b2mod_boundary_namelist(ns, mpg, .true.)
      call xertst (nbc.le.nbcd, 'increase DEF_NBC in DIMENSIONS.F')
      call xertst (nniso.le.nnisod, 'increase DEF_NISO in DIMENSIONS.F')
      call ipgeti ('b2stbr_neutrals_namelist', neutrals_namelist)
      call alloc_b2mod_neutrals
     > (ns, natm, nmol, nion, nstraid, nsts, sput_dst)
      if (neutrals_namelist.eq.1 .or. use_eirene.ne.0)
     > call read_neutrals_namelist(ns, mpg, switch, .false.)
      call xertst (nstrat.le.nstraid,
     . 'Increase DEF_NSTRA in DIMENSIONS.F')
#ifdef WG_TODO
cxpb  Bring in wall data and compute oft-needed wall quantities
      write(*,*) 'Calling alloc_b2mod_sputter'
      call alloc_b2mod_sputter(nx,ny,ns)
      call alloc_b2mod_layer(nx,ny,nncut,nlim,nsts)
      call alloc_b2mod_wall_init(nx,ny,nncut)
      call alloc_b2mod_wall(ns)
      call init_wall(nx,ny)
#endif
      call read_b2mod_wall_namelist(99,'    ',ex)
      call xertst (ndepth_nml.gt.1.or.plate_model.eq.0,
     .  'Inconsistent values of ndepth_nml and plate_model')
#ifdef WG_TODO
      call get_wall_data(nx,ny,ns,switch%BoRiS,ex)
#endif
#ifdef B25_EIRENE
      filename=fort_lc//'44'
      call find_file(filename,ex)
#else
      ex = .false.
#endif
      if (.not.ex.or.use_eirene.eq.0) then
        natmi=0
        textan(0)=' ZERO'
        nmoli=0
        textmn(0)=' ZERO'
        nioni=0
        textin(0)=' ZERO'
        call eirene_mc_init(nCv,nFc)
        nlimi=0
      else
        open(unit=44,file=filename)
        call ngread(mpg,gm,ns,switch,ismain,iflag)
        if (natmi.gt.0)
     .   textwl(ns:ns+natmi-1)=textan(0:natmi-1)
        if (nmoli.gt.0)
     .   textwl(ns+natmi:ns+natmi+nmoli-1)=textmn(0:nmoli-1)
        if (nioni.gt.0)
     .   textwl(ns+natmi+nmoli:ns+natmi+nmoli+nioni-1)=textin(0:nioni-1)
      endif

!XXX
      write(*,*) 'b2mds ...'

#ifdef MDSPLUS
      textdrift(0)='Diamagnetic'
      textdrift(1)='Radial'
      nsw=ns+natmi+nmoli+nioni

      call mds_put_int('\TOP.SNAPSHOT.DIMENSIONS:NCV',nCv)
      call mds_put_int('\TOP.SNAPSHOT.DIMENSIONS:NFC',nFc)
      call mds_put_int('\TOP.SNAPSHOT.DIMENSIONS:NVX',nVx)
      call mds_put_int('\TOP.SNAPSHOT.DIMENSIONS:NS',ns)
      call mds_put_int('\TOP.SNAPSHOT.DIMENSIONS:NATM',natmi)
      call mds_put_int('\TOP.SNAPSHOT.DIMENSIONS:NMOL',nmoli)
      call mds_put_int('\TOP.SNAPSHOT.DIMENSIONS:NION',nioni)
      call mds_put_int('\TOP.SNAPSHOT.DIMENSIONS:NSW',nsw)
      call mds_put_int('\TOP.SNAPSHOT.DIMENSIONS:NWALL',nbnd)
#ifdef WG_TODO
      call mds_put_int('\TOP.SNAPSHOT.DIMENSIONS:TARG1',1)
      call mds_put_int('\TOP.SNAPSHOT.DIMENSIONS:TARG2',-1)
      call mds_put_int('\TOP.SNAPSHOT.DIMENSIONS:TARG3',-1)
      call mds_put_int('\TOP.SNAPSHOT.DIMENSIONS:TARG4',nx)
#endif
      call mds_put_int('\TOP.SNAPSHOT.DIMENSIONS:IMP',ifsepimp)
      call mds_put_int('\TOP.SNAPSHOT.DIMENSIONS:OMP',ifsepomp)
      call mds_put_int('\TOP.SNAPSHOT.DIMENSIONS:SEP',mpg%iFssep)

      tdi='build_with_units($,"m")'

      mds_signal_name='\TOP.SNAPSHOT.GRID:CR'
      write(*,*) 'mdsput2 of ',trim(mds_signal_name),' with ',trim(tdi)
#if ( defined(GFORTRAN) && COMPILER_MAJOR_VERSION > 9 ) || defined(G95)
      status=mdsput2(
     * trim(mds_signal_name)//char(0),trim(tdi)//char(0),
     * descr2_i1(idtype_double,nCv,0),gm%cvX,0)
#else
      status=mdsput2(
     * trim(mds_signal_name)//char(0),trim(tdi)//char(0),
     * descr2(idtype_double,nCv,0),gm%cvX,0)
#endif
      write(*,*) 'mdsput2'
      call mds_err_msg(status)

      mds_signal_name='\TOP.SNAPSHOT.GRID:CZ'
      write(*,*) 'mdsput2 of ',trim(mds_signal_name),' with ',trim(tdi)
#if ( defined(GFORTRAN) && COMPILER_MAJOR_VERSION > 9 ) || defined(G95)
      status=mdsput2(
     * trim(mds_signal_name)//char(0),
     * trim(tdi)//char(0),
     * descr2_i1(idtype_double,nCv,0),gm%cvY,0)
#else
      status=mdsput2(
     * trim(mds_signal_name)//char(0),
     * trim(tdi)//char(0),
     * descr2(idtype_double,nCv,0),gm%cvY,0)
#endif
      write(*,*) 'mdsput2'
      call mds_err_msg(status)

      mds_signal_name='\TOP.SNAPSHOT.GRID:CR_F'
      write(*,*) 'mdsput2 of ',trim(mds_signal_name),' with ',trim(tdi)
#if ( defined(GFORTRAN) && COMPILER_MAJOR_VERSION > 9 ) || defined(G95)
      status=mdsput2(
     * trim(mds_signal_name)//char(0),trim(tdi)//char(0),
     * descr2_i1(idtype_double,nFc,0),crf,0)
#else
      status=mdsput2(
     * trim(mds_signal_name)//char(0),trim(tdi)//char(0),
     * descr2(idtype_double,nFc,0),crf,0)
#endif
      write(*,*) 'mdsput2'
      call mds_err_msg(status)

      mds_signal_name='\TOP.SNAPSHOT.GRID:CZ_F'
      write(*,*) 'mdsput2 of ',trim(mds_signal_name),' with ',trim(tdi)
#if ( defined(GFORTRAN) && COMPILER_MAJOR_VERSION > 9 ) || defined(G95)
      status=mdsput2(
     * trim(mds_signal_name)//char(0),
     * trim(tdi)//char(0),
     * descr2_i1(idtype_double,nFc,0),czf,0)
#else
      status=mdsput2(
     * trim(mds_signal_name)//char(0),
     * trim(tdi)//char(0),
     * descr2(idtype_double,nFc,0),czf,0)
#endif
      write(*,*) 'mdsput2'
      call mds_err_msg(status)

      mds_signal_name='\TOP.SNAPSHOT.GRID:R'
      write(*,*) 'mdsput2 of ',trim(mds_signal_name),' with ',trim(tdi)
#if ( defined(GFORTRAN) && COMPILER_MAJOR_VERSION > 9 ) || defined(G95)
      status=mdsput2(
     * trim(mds_signal_name)//char(0),
     * trim(tdi)//char(0),
     * descr2_i1(idtype_double,nVx,0),gm%vxX,
     * 0)
#else
      status=mdsput2(
     * trim(mds_signal_name)//char(0),
     * trim(tdi)//char(0),
     * descr2(idtype_double,nVx,0),gm%vxX,
     * 0)
#endif
      write(*,*) 'mdsput2'
      call mds_err_msg(status)

      mds_signal_name='\TOP.SNAPSHOT.GRID:Z'
      write(*,*) 'mdsput2 of ',trim(mds_signal_name),' with ',trim(tdi)
#if ( defined(GFORTRAN) && COMPILER_MAJOR_VERSION > 9 ) || defined(G95)
      status=mdsput2(
     * trim(mds_signal_name)//char(0),
     * trim(tdi)//char(0),
     * descr2_i1(idtype_double,nVx,0),gm%vxY,
     * 0)
#else
      status=mdsput2(
     * trim(mds_signal_name)//char(0),
     * trim(tdi)//char(0),
     * descr2(idtype_double,nVx,0),gm%vxY,
     * 0)
#endif
      write(*,*) 'mdsput2'
      call mds_err_msg(status)

#ifdef WG_TODO
      do ix=-1,nx
        dsrad(ix,jsep+1)=hy1(ix,jsep+1)/2.0_R8
        do iy=jsep+2,ny
          dsrad(ix,iy)=dsrad(ix,iy-1)+(hy1(ix,iy-1)+hy1(ix,iy))/2.0_R8
        enddo
        dsrad(ix,jsep)=-hy1(ix,jsep)/2.0_R8
        do iy=jsep-1,-1,-1
          dsrad(ix,iy)=dsrad(ix,iy+1)-(hy1(ix,iy+1)+hy1(ix,iy))/2.0_R8
        enddo
      enddo
      write(*,*) 'dsrad :',(dsrad(jxa,iy),iy=-1,ny)
      do iy=-1,ny
        dspol(-1,iy)=0.0_R8
        dspar(-1,iy)=0.0_R8
        do ix=0,nx
          dspol(ix,iy)=dspol(ix-1,iy)+(hx(ix-1,iy)+hx(ix,iy))/2.0_R8
          dspar(ix,iy)=dspar(ix-1,iy)+
     1     (hx(ix-1,iy)*abs(bb(ix-1,iy,3)/bb(ix-1,iy,0))+
     2      hx(ix,iy)*abs(bb(ix,iy,3)/bb(ix,iy,0)))/2.0_R8
        enddo
      enddo
      write(*,*) 'dspol :',(dspol(ix,jsep),ix=-1,nx)
      write(*,*) 'dspar :',(dspar(ix,jsep),ix=-1,nx)
      mds_signal_name='\TOP.SNAPSHOT.GRID:DSRAD'
      write(*,*) 'mdsput2 of ',trim(mds_signal_name),' with ',trim(tdi)
#if ( defined(GFORTRAN) && COMPILER_MAJOR_VERSION > 9 ) || defined(G95)
      status=mdsput2(
     * trim(mds_signal_name)//char(0),trim(tdi)//char(0),
     * descr2_i2(idtype_double,nx+2,ny+2,0),dsrad,0)
#else
      status=mdsput2(
     * trim(mds_signal_name)//char(0),trim(tdi)//char(0),
     * descr2(idtype_double,nx+2,ny+2,0),dsrad,0)
#endif
      write(*,*) 'mdsput2'
      call mds_err_msg(status)
      mds_signal_name='\TOP.SNAPSHOT.GRID:DSPOL'
      write(*,*) 'mdsput2 of ',trim(mds_signal_name),' with ',trim(tdi)
#if ( defined(GFORTRAN) && COMPILER_MAJOR_VERSION > 9 ) || defined(G95)
      status=mdsput2(
     * trim(mds_signal_name)//char(0),trim(tdi)//char(0),
     * descr2_i2(idtype_double,nx+2,ny+2,0),dspol,0)
#else
      status=mdsput2(
     * trim(mds_signal_name)//char(0),trim(tdi)//char(0),
     * descr2(idtype_double,nx+2,ny+2,0),dspol,0)
#endif
      write(*,*) 'mdsput2'
      call mds_err_msg(status)
      mds_signal_name='\TOP.SNAPSHOT.GRID:DSPAR'
      write(*,*) 'mdsput2 of ',trim(mds_signal_name),' with ',trim(tdi)
#if ( defined(GFORTRAN) && COMPILER_MAJOR_VERSION > 9 ) || defined(G95)
      status=mdsput2(
     * trim(mds_signal_name)//char(0),trim(tdi)//char(0),
     * descr2_i2(idtype_double,nx+2,ny+2,0),dspar,0)
#else
      status=mdsput2(
     * trim(mds_signal_name)//char(0),trim(tdi)//char(0),
     * descr2(idtype_double,nx+2,ny+2,0),dspar,0)
#endif
      write(*,*) 'mdsput2'
      call mds_err_msg(status)
#endif
c
      dbsolpsversion=b2fplasma_version
      call mds_put_string('\TOP.IDENT:SOLPSVERSION',dbsolpsversion)

c change dbversion below when MDS+ output format/data is amended/augmented
c 1: fht uses fhe and fhi instead of fhet and fhit
c 2: added fna{xy}_{35}2
c 3: fixed definition of \top.ident.impsep:kyi0 and \top.ident.ompsep:kyi0
c 4: added poti and pot
c 5: changed hcib and hci0 to include the radial neutral convective flux limits
c 6: added sxperp and qz
c 7: added new fields to identify b2f* file versions and svn version
c 8: changed from hy to hy1 in dsrad calculation
c 9: changed definition of za, pot and poti
c    added zamin, zamax, rza, rz2, rpt and rpi
c 10: added uadia, fch{xy}_{35}2, fch_p and kinrgy for recovery purposes
c 11: added fllimi & fllime
c 12: added ceqp
c 13: correctly save tmhacore, tmhasol & tmhadiv
c 14: includes contributions from external species in ne, pr, fht
c 15: correction of sign of electric energy flux
c 16: add \TOP.TIMEDEP.INPUT:DS \TOP.SNAPSHOT:FHJX \TOP.SNAPSHOT:FHJY
c     \TOP.TALLIES:fhjxreg \TOP.TALLIES:fhjyreg
c 17: corrected cr and cz locations for inner and outer separatrix midplanes
c 18: introduced new quantities from SOLPS_ITER
c     ENEUTRAD now with third dimension for atomic species
c 19: added GIT version
c 20: output na0, ua0, ne0, te0, ti0, po0
c 21: output smaf, smag, smav, smfr, smo, smpr, smpt, smq, smth, resco, reshe, reshi, resmo, resmt, respo
c 22: added various quantities for recovery of advanced SOLPS-ITER features
c     fna_mdf, fhe_mdf, fhi_mdf, fna_fcor, fna_nodrift, fna_he,
c     fnaPSch, fhePSch, fhiPSch, fna_eir, fne_eir, fhe_eir, fhi_eir,
c     fni_32, fni_52, fne_32, fne_52,
c     fchdia, fchin, fchvispar, fchvisper, fchvisq, fchinert,
c     vaecrb, vadia, wadia, veecrb, vedia, floe_noc, floi_noc
c 23: fht definition changed to 5/2 formulation
c 25: debugged indexing error for eneutrad
c 26: using total field for computing dspar
c 27: Added EMOLRAD and EIONRAD arrays
c 29: Same indexing error was also present in emolrad and eionrad
c 30: changed quantities to unstructured grid
c 32: Corrected units for srcml and edissml

      dbversion=32
      call mds_put_int('\TOP.IDENT:VERSION',dbversion)

      textpl_cat=textpl(0)
      do is=1,ns-1
        textpl_cat=trim(textpl_cat) // ' ' // textpl(is)
      enddo
      call mds_put_string_ns('\TOP.SNAPSHOT:TEXTPL',ns,textpl)
      call mds_put_string('\TOP.IDENT:SPECIES',textpl_cat)
      dbspecies=trim(textpl_cat)
c
      textwl_cat=textwl(0)
      do is=1,ns-1
        textwl_cat=trim(textwl_cat) // ' ' // textwl(is)
      enddo
      if (natmi.gt.0) then
        do is=ns,ns+natmi-1
          textwl_cat=trim(textwl_cat) // ' ' // textwl(is)
        enddo
      end if
      if (nmoli.gt.0) then
        do is=ns+natmi,ns+natmi+nmoli-1
          textwl_cat=trim(textwl_cat) // ' ' // textwl(is)
        enddo
      end if
      if (nioni.gt.0) then
        do is=ns+natmi+nmoli,ns+natmi+nmoli+nioni-1
          textwl_cat=trim(textwl_cat) // ' ' // textwl(is)
        enddo
      end if
      call mds_put_string_ns('\TOP.SNAPSHOT:TEXTWL',
     1 ns+natmi+nmoli+nioni,textwl)
      call mds_put_string('\TOP.IDENT:WALL_SPECIES',textwl_cat)
c
      tdi='$'
      call mds_put_double_ns('\TOP.SNAPSHOT.GRID:ZAMIN','',
     1 ns,zamin(0:ns-1))
      call mds_put_double_ns('\TOP.SNAPSHOT.GRID:ZAMAX','',
     1 ns,zamax(0:ns-1))
      za(0:ns-1)=(zamin(0:ns-1)+zamax(0:ns-1))/2.0_R8
      call mds_put_double_ns('\TOP.SNAPSHOT.GRID:ZA','',ns,za(0:ns-1))
      call mds_put_double_ns('\TOP.SNAPSHOT.GRID:ZN','',ns,zn(0:ns-1))
      call mds_put_double_ns('\TOP.SNAPSHOT.GRID:AM','',ns,am(0:ns-1))
      poti(0:ns-1)=(potimin(0:ns-1)+potimax(0:ns-1))/2.0_R8
      call mds_put_double_ns('\TOP.SNAPSHOT.GRID:POTI','eV',ns,poti)
      pot(0:ns-1)=(potmin(0:ns-1)+potmax(0:ns-1))/2.0_R8
      call mds_put_double_ns('\TOP.SNAPSHOT.GRID:POT','eV',ns,pot)
      if(natmi.gt.0)
     1 call mds_put_string_ns('\TOP.SNAPSHOT:TEXTAN',natmi,textan)
      if(nmoli.gt.0)
     1 call mds_put_string_ns('\TOP.SNAPSHOT:TEXTMN',nmoli,textmn)
      if(nioni.gt.0)
     1 call mds_put_string_ns('\TOP.SNAPSHOT:TEXTIN',nioni,textin)
      call mds_put_string_ns('\TOP.SNAPSHOT:TEXTCOMP',4,textcomponents)
      call mds_put_string_ns('\TOP.SNAPSHOT:TEXTDRIFT',2,textdrift)
      call mds_put_string_ns('\TOP.SNAPSHOT:TEXTLIN',4,text0123)
c
      filename='mesh.extra'
      call find_file(filename,ex)
      if(ex) then
        open(99,file=filename,status='old')
        icnt=0
   10   continue
        read(99,*,end=20) (rdummy,j=1,4)
        icnt=icnt+1
        goto 10
   20   continue
        write(*,*) 'Vessel size = ',icnt
        if (icnt.gt.0) then
          allocate(vesl(4,icnt))
          close(99)
          open(99,file='mesh.extra',status='old')
          do i=1,icnt
            read(99,*) (vesl(j,i),j=1,4)
          enddo
          close(99)
          tdi='build_with_units($,"m")'
          write(*,'(a,a,a,a)')
     1     'mdsput2 of ','\TOP.SNAPSHOT.GRID:VESSEL',' with ',trim(tdi)
#if ( defined(GFORTRAN) && COMPILER_MAJOR_VERSION > 9 ) || defined(G95)
          status=mdsput2(
     *     '\TOP.SNAPSHOT.GRID:VESSEL'//char(0),trim(tdi)//char(0),
     *     descr2_i2(idtype_double,4,icnt,0),vesl,
     *     0)
#else
          status=mdsput2(
     *     '\TOP.SNAPSHOT.GRID:VESSEL'//char(0),trim(tdi)//char(0),
     *     descr2(idtype_double,4,icnt,0),vesl,
     *     0)
#endif
          call mds_err_msg(status)
          deallocate(vesl)
        else
          write(*,*) 'Empty mesh.extra file found !'
        endif
      else
        write(*,*) 'No mesh.extra file found !'
      endif

#ifdef WG_TODO
      itmp=leftix+1
      call mds_put_int_rz(
     1 '\TOP.SNAPSHOT.GRID:LEFTIX',nx+2,ny+2,itmp)
      itmp=leftiy+1
      call mds_put_int_rz(
     1 '\TOP.SNAPSHOT.GRID:LEFTIY',nx+2,ny+2,itmp)
      itmp=rightix+1
      call mds_put_int_rz(
     1 '\TOP.SNAPSHOT.GRID:RIGHTIX',nx+2,ny+2,itmp)
      itmp=rightiy+1
      call mds_put_int_rz(
     1 '\TOP.SNAPSHOT.GRID:RIGHTIY',nx+2,ny+2,itmp)
      itmp=bottomix+1
      call mds_put_int_rz(
     1 '\TOP.SNAPSHOT.GRID:BOTTOMIX',nx+2,ny+2,itmp)
      itmp=bottomiy+1
      call mds_put_int_rz(
     1 '\TOP.SNAPSHOT.GRID:BOTTOMIY',nx+2,ny+2,itmp)
      itmp=topix+1
      call mds_put_int_rz(
     1 '\TOP.SNAPSHOT.GRID:TOPIX',nx+2,ny+2,itmp)
      itmp=topiy+1
      call mds_put_int_rz(
     1 '\TOP.SNAPSHOT.GRID:TOPIY',nx+2,ny+2,itmp)
#endif
      call mds_put_int_us(
     1 '\TOP.SNAPSHOT.GRID:REGVOL',nCv,mpg%cvReg)
      call mds_put_int_us(
     1 '\TOP.SNAPSHOT.GRID:REGFL',nFc,mpg%fcReg)
#ifdef WG_TODO
      call mds_put_int_rz(
     1 '\TOP.SNAPSHOT.GRID:XYMAP',nx+2,ny+2,xymap)
      call mds_put_int_ns(
     1 '\TOP.SNAPSHOT.GRID:IMAPX',nbnd,imapx)
      call mds_put_int_ns(
     1 '\TOP.SNAPSHOT.GRID:IMAPY',nbnd,imapy)
#endif

      temp=state%pl%te/qe
      call mds_put_double_us('\TOP.SNAPSHOT:TE','eV',nCv,temp,'')
      temp=state%psnl%te/qe
      call mds_put_double_us('\TOP.SNAPSHOT:TE0','eV',nCv,temp,'')
      temp=state%pl%ti/qe
      call mds_put_double_us('\TOP.SNAPSHOT:TI','eV',nCv,temp,'')
      temp=state%psnl%ti/qe
      call mds_put_double_us('\TOP.SNAPSHOT:TI0','eV',nCv,temp,'')
      call mds_put_double_us('\TOP.SNAPSHOT:NE','/m^3',
     &     nCv,state%dv%ne,'')
      call mds_put_double_us('\TOP.SNAPSHOT:NE0','/m^3',
     &     nCv,state%psnl%ne,'')
      call mds_put_double_us('\TOP.SNAPSHOT:PR','Pa',nCv,pr,'')
      call mds_put_double_us('\TOP.SNAPSHOT:PO','V',nCv,state%pl%po,'')
      call mds_put_double_us('\TOP.SNAPSHOT:PO0','V',
     &     nCv,state%psnl%po,'')
      call mds_put_double_us('\TOP.SNAPSHOT:VOL','m^3',nCv,gm%cvVol,'')
      call mds_put_double_us('\TOP.SNAPSHOT:HX','m',nCv,gm%cvHx,'')
      call mds_put_double_us('\TOP.SNAPSHOT:fllime','',
     1     nFc,state%dv%fllime,'')
      call mds_put_double_us('\TOP.SNAPSHOT:fllimi','',
     1     nFc,state%dv%fllimi,'')
      call mds_put_double_us('\TOP.SNAPSHOT:ceqp','',
     1     nCv,state%co%ceqp,'')
      call mds_put_double_usns('\TOP.SNAPSHOT:RZA','',
     1     nCv,ns,state%rt%rza,'TEXTPL','')
      call mds_put_double_usns('\TOP.SNAPSHOT:RZ2','',
     1     nCv,ns,state%rt%rz2,'TEXTPL','')
      call mds_put_double_usns('\TOP.SNAPSHOT:RPT','eV',
     1     nCv,ns,state%rt%rpt,'TEXTPL','')
      call mds_put_double_usns('\TOP.SNAPSHOT:RPI','eV',
     1     nCv,ns,state%rt%rpi,'TEXTPL','')
      call mds_put_double_usns('\TOP.SNAPSHOT:NA','/m^3',
     1     nCv,ns,state%pl%na,'TEXTPL','')
      call mds_put_double_usns('\TOP.SNAPSHOT:NA0','/m^3',
     1     nCv,ns,state%psnl%na,'TEXTPL','')
      call mds_put_double_usns('\TOP.SNAPSHOT:UA','m/s',
     1     nCv,ns,state%pl%ua,'TEXTPL','')
      call mds_put_double_usns('\TOP.SNAPSHOT:UA0','m/s',
     1     nCv,ns,state%psnl%ua,'TEXTPL','')
      call mds_put_double_usdns('\TOP.SNAPSHOT:UADIA','m/s',
     1     nFc,2,ns,state%dv%uadia,'TEXTDRIFT','TEXTPL','')
      call mds_put_double_usdns('\TOP.SNAPSHOT:VADIA','m/s',
     1     nFc,2,ns,state%dv%vadia,'TEXTDRIFT','TEXTPL','')
      call mds_put_double_usdns('\TOP.SNAPSHOT:WADIA','m/s',
     1     nFc,2,ns,state%dv%wadia,'TEXTDRIFT','TEXTPL','')
      call mds_put_double_usdns('\TOP.SNAPSHOT:VAECRB','m/s',
     1     nFc,2,ns,state%dv%vaecrb,'TEXTDRIFT','TEXTPL','')
      call mds_put_double_usns('\TOP.SNAPSHOT:VEDIA','m/s',
     1     nFc,2,state%dv%vedia,'TEXTDRIFT','')
      call mds_put_double_usns('\TOP.SNAPSHOT:VEECRB','m/s',
     1     nFc,2,state%dv%veecrb,'TEXTDRIFT','')
      call mds_put_double_usns('\TOP.SNAPSHOT:KINRGY','J',
     1     nCv,ns,state%dv%kinrgy,'TEXTPL','')
      call mds_put_double_usns('\TOP.SNAPSHOT:B','T',
     1     nCv,4,gm%cvBb,'TEXTCOMP','')
      call mds_put_double_usns('\TOP.SNAPSHOT:QZ','',
     1     nCv,2,gm%cvQgam,'TEXTCOMP','')
      call mds_put_double_us('\TOP.SNAPSHOT:BRHE','W',
     1     nCv,state%srw%b2stbr_she,'')
      call mds_put_double_us('\TOP.SNAPSHOT:BRHI','W',
     1     nCv,state%srw%b2stbr_shi,'')
      call mds_put_double_usns('\TOP.SNAPSHOT:BRNA','particles/s',
     1     nCv,ns,state%srw%b2stbr_sna,'TEXTPL','')
      call mds_put_double_usns('\TOP.SNAPSHOT:BRMO','N',
     1     nCv,ns,state%srw%b2stbr_smo,'TEXTPL','')

      call mds_put_double_usns('\TOP.SNAPSHOT:RQAHE','W',
     1     nCv,ns,state%srw%rqahe,'TEXTPL','')
      call mds_put_double_usns('\TOP.SNAPSHOT:RQRAD','W',
     1     nCv,ns,state%srw%rqrad,'TEXTPL','')
      call mds_put_double_usns('\TOP.SNAPSHOT:RQBRM','W',
     1     nCv,ns,state%srw%rqbrm,'TEXTPL','')
      call mds_put_double_usns('\TOP.SNAPSHOT:RSAHI','W',
     1     nCv,ns,state%srw%rsahi,'TEXTPL','')
      call mds_put_double_usns('\TOP.SNAPSHOT:RSANA','/s',
     1     nCv,ns,state%srw%rsana,'TEXTPL','')
      call mds_put_double_usns('\TOP.SNAPSHOT:RSAMO','m/s^2',
     1     nCv,ns,state%srw%rsamo,'TEXTPL','')
      call mds_put_double_usns('\TOP.SNAPSHOT:RRAHI','W',
     1     nCv,ns,state%srw%rrahi,'TEXTPL','')
      call mds_put_double_usns('\TOP.SNAPSHOT:RRANA','/s',
     1     nCv,ns,state%srw%rrana,'TEXTPL','')
      call mds_put_double_usns('\TOP.SNAPSHOT:RRAMO','m/s^2',
     1     nCv,ns,state%srw%rramo,'TEXTPL','')
      call mds_put_double_usns('\TOP.SNAPSHOT:RCXHI','W',
     1     nCv,ns,state%srw%rcxhi,'TEXTPL','')
      call mds_put_double_usns('\TOP.SNAPSHOT:RCXNA','/s',
     1     nCv,ns,state%srw%rcxna,'TEXTPL','')
      call mds_put_double_usns('\TOP.SNAPSHOT:RCXMO','m/s^2',
     1     nCv,ns,state%srw%rcxmo,'TEXTPL','')

      call mds_put_double_usns('\TOP.SNAPSHOT:D','m^2.s^-1',
     1     nCv,ns,state%co%dna0,'TEXTPL','')
      do is = 0, ns-1
        if (is_neutral(is)) then
          tmp(:,is)=state%co%dpa0(:,is)*state%pl%tn(:)
        else
          tmp(:,is)=state%co%dpa0(:,is)*
     *             (state%rt%rza(:,is)*state%pl%te(:)+state%pl%ti(:))
        end if
      end do
      call mds_put_double_usns('\TOP.SNAPSHOT:DP','m^2.s^-1',
     1     nCv,ns,tmp,'TEXTPL','')
      temp=state%co%hce0/state%dv%ne
      call mds_put_double_us('\TOP.SNAPSHOT:KYE','m^2.s^-1',
     1     nCv,temp,'')
      forall (iCv=1:nCv,is=0:ns-1)
     1 tmp(iCv,is)=state%co%hcib(iCv,is)/state%pl%na(iCv,is)
      call mds_put_double_usns('\TOP.SNAPSHOT:KYI','m^2.s^-1',
     1     nCv,ns,tmp,'TEXTPL','')
      temp=state%co%hci0(:)/state%dv%ni(:,0)
      call mds_put_double_us('\TOP.SNAPSHOT:KYI0','m^2.s^-1',
     1     nCv,temp,'')
      tmp(:,:)=state%co%vla0(:,0,:)
      call mds_put_double_usns('\TOP.SNAPSHOT:VLAX','m.s^-1',
     1     nCv,ns,tmp,'TEXTPL','')
      tmp(:,:)=state%co%vla0(:,1,:)
      call mds_put_double_usns('\TOP.SNAPSHOT:VLAY','m.s^-1',
     1     nCv,ns,tmp,'TEXTPL','')
      forall (iCv=1:nCv,is=0:ns-1)
     1 tmp(iCv,is)=state%co%vsa0(iCv,is)/(mp*am(is)*state%pl%na(iCv,is))
      call mds_put_double_usns('\TOP.SNAPSHOT:VISC','m.s^-1',
     1     nCv,ns,tmp,'TEXTPL','')
      temp=state%co%sig0/(qe*state%dv%ne)
      call mds_put_double_us('\TOP.SNAPSHOT:SIG','mho.m^-1',
     1     nCv,temp,'')
      temp=state%co%alf0(:)/(state%dv%ne(:)*sqrt(qe/state%pl%te(:)))
      call mds_put_double_us('\TOP.SNAPSHOT:ALF','mho.m^-1',
     1     nCv,temp,'')

      flux_x(:,:)=state%dv%fna(:,0,:)
      call mds_put_double_usns('\TOP.SNAPSHOT:FNAX','/s',
     1     nFc,ns,flux_x,'TEXTPL','_F')
      flux_y(:,:)=state%dv%fna(:,1,:)
      call mds_put_double_usns('\TOP.SNAPSHOT:FNAY','/s',
     1     nFc,ns,flux_y,'TEXTPL','_F')
      flux_x(:,:)=state%dv%fna_32(:,0,:)
      call mds_put_double_usns('\TOP.SNAPSHOT:FNAX_32','/s',
     1     nFc,ns,flux_x,'TEXTPL','_F')
      flux_y(:,:)=state%dv%fna_32(:,1,:)
      call mds_put_double_usns('\TOP.SNAPSHOT:FNAY_32','/s',
     1     nFc,ns,flux_y,'TEXTPL','_F')
      flux_x(:,:)=state%dv%fna_52(:,0,:)
      call mds_put_double_usns('\TOP.SNAPSHOT:FNAX_52','/s',
     1     nFc,ns,flux_x,'TEXTPL','_F')
      flux_y(:,:)=state%dv%fna_52(:,1,:)
      call mds_put_double_usns('\TOP.SNAPSHOT:FNAY_52','/s',
     1     nFc,ns,flux_y,'TEXTPL','_F')
      flux_x(:,:)=state%dv%fna_mdf(:,0,:)
      call mds_put_double_usns('\TOP.SNAPSHOT:FNAX_MDF','/s',
     1     nFc,ns,flux_x,'TEXTPL','_F')
      flux_y(:,:)=state%dv%fna_mdf(:,1,:)
      call mds_put_double_usns('\TOP.SNAPSHOT:FNAY_MDF','/s',
     1     nFc,ns,flux_y,'TEXTPL','_F')
      flux_x(:,:)=state%dv%fna_fcor(:,0,:)
      call mds_put_double_usns('\TOP.SNAPSHOT:FNAX_FCOR','/s',
     1     nFc,ns,flux_x,'TEXTPL','_F')
      flux_y(:,:)=state%dv%fna_fcor(:,1,:)
      call mds_put_double_usns('\TOP.SNAPSHOT:FNAY_FCOR','/s',
     1     nFc,ns,flux_y,'TEXTPL','_F')
      flux_x(:,:)=state%dv%fna_nodrift(:,0,:)
      call mds_put_double_usns('\TOP.SNAPSHOT:FNAX_NODRIFT','/s',
     1     nFc,ns,flux_x,'TEXTPL','_F')
      flux_y(:,:)=state%dv%fna_nodrift(:,1,:)
      call mds_put_double_usns('\TOP.SNAPSHOT:FNAY_NODRIFT','/s',
     1     nFc,ns,flux_y,'TEXTPL','_F')
      flux_x(:,:)=state%dv%fna_he(:,0,:)
      call mds_put_double_usns('\TOP.SNAPSHOT:FNAX_HE','/s',
     1     nFc,ns,flux_x,'TEXTPL','_F')
      flux_y(:,:)=state%dv%fna_he(:,1,:)
      call mds_put_double_usns('\TOP.SNAPSHOT:FNAY_HE','/s',
     1     nFc,ns,flux_y,'TEXTPL','_F')
      flux_x(:,:)=state%dv%fnaPSch(:,0,:)
      call mds_put_double_usns('\TOP.SNAPSHOT:FNAPSCHX','/s',
     1     nFc,ns,flux_x,'TEXTPL','_F')
      flux_y(:,:)=state%dv%fnaPSch(:,1,:)
      call mds_put_double_usns('\TOP.SNAPSHOT:FNAPSCHY','/s',
     1     nFc,ns,flux_y,'TEXTPL','_F')
      flux_x(:,:)=state%dv%fna_eir(:,0,:)
      call mds_put_double_usns('\TOP.SNAPSHOT:FNAX_EIR','/s',
     1     nFc,ns,flux_x,'TEXTPL','_F')
      flux_y(:,:)=state%dv%fna_eir(:,1,:)
      call mds_put_double_usns('\TOP.SNAPSHOT:FNAY_EIR','/s',
     1     nFc,ns,flux_y,'TEXTPL','_F')
      flux_x(:,:)=state%dv%fmo(:,0,:)
      call mds_put_double_usns('\TOP.SNAPSHOT:FMOX','N',
     1     nFc,ns,flux_x,'TEXTPL','_F')
      flux_y(:,:)=state%dv%fmo(:,1,:)
      call mds_put_double_usns('\TOP.SNAPSHOT:FMOY','N',
     1     nFc,ns,flux_y,'TEXTPL','_F')
      flux_x(:,:)=state%dv%fhm(:,0,:)
      call mds_put_double_usns('\TOP.SNAPSHOT:FHMX','W',
     1     nFc,ns,flux_x,'TEXTPL','_F')
      flux_y(:,:)=state%dv%fhm(:,1,:)
      call mds_put_double_usns('\TOP.SNAPSHOT:FHMY','W',
     1     nFc,ns,flux_y,'TEXTPL','_F')
      flux_x(:,:)=state%dv%fhp(:,0,:)
      call mds_put_double_usns('\TOP.SNAPSHOT:FHPX','W',
     1     nFc,ns,flux_x,'TEXTPL','_F')
      flux_y(:,:)=state%dv%fhp(:,1,:)
      call mds_put_double_usns('\TOP.SNAPSHOT:FHPY','W',
     1     nFc,ns,flux_y,'TEXTPL','_F')
      flux_x(:,0)=state%dv%fne_eir(:,0)
      call mds_put_double_us('\TOP.SNAPSHOT:FNEX_EIR','/s',
     1     nFc,flux_x,'_F')
      flux_y(:,0)=state%dv%fne_eir(:,1)
      call mds_put_double_us('\TOP.SNAPSHOT:FNEY_EIR','/s',
     1     nFc,flux_y,'_F')
      flux_x(:,0)=state%dv%fne_32(:,0)
      call mds_put_double_us('\TOP.SNAPSHOT:FNEX_32','/s',
     1     nFc,flux_x,'_F')
      flux_y(:,0)=state%dv%fne_32(:,1)
      call mds_put_double_us('\TOP.SNAPSHOT:FNEY_32','/s',
     1     nFc,flux_y,'_F')
      flux_x(:,0)=state%dv%fne_52(:,0)
      call mds_put_double_us('\TOP.SNAPSHOT:FNEX_52','/s',
     1     nFc,flux_x,'_F')
      flux_y(:,0)=state%dv%fne_52(:,1)
      call mds_put_double_us('\TOP.SNAPSHOT:FNEY_52','/s',
     1     nFc,flux_y,'_F')
      flux_x(:,0)=state%dv%fni_32(:,0)
      call mds_put_double_us('\TOP.SNAPSHOT:FNIX_32','/s',
     1     nFc,flux_x,'_F')
      flux_y(:,0)=state%dv%fni_32(:,1)
      call mds_put_double_us('\TOP.SNAPSHOT:FNIY_32','/s',
     1     nFc,flux_y,'_F')
      flux_x(:,0)=state%dv%fni_52(:,0)
      call mds_put_double_us('\TOP.SNAPSHOT:FNIX_52','/s',
     1     nFc,flux_x,'_F')
      flux_y(:,0)=state%dv%fni_52(:,1)
      call mds_put_double_us('\TOP.SNAPSHOT:FNIY_52','/s',
     1     nFc,flux_y,'_F')
      flux_x(:,0)=state%dv%fhe(:,0)
      call mds_put_double_us('\TOP.SNAPSHOT:FHEX','W',
     1     nFc,flux_x,'_F')
      flux_y(:,0)=state%dv%fhe(:,1)
      call mds_put_double_us('\TOP.SNAPSHOT:FHEY','W',
     1     nFc,flux_y,'_F')
      flux_x(:,0)=state%dv%fhe_mdf(:,0)
      call mds_put_double_us('\TOP.SNAPSHOT:FHEX_MDF','W',
     1     nFc,flux_x,'_F')
      flux_y(:,0)=state%dv%fhe_mdf(:,1)
      call mds_put_double_us('\TOP.SNAPSHOT:FHEY_MDF','W',
     1     nFc,flux_y,'_F')
      flux_x(:,0)=state%dv%fhePSch(:,0)
      call mds_put_double_us('\TOP.SNAPSHOT:FHEPSCHX','W',
     1     nFc,flux_x,'_F')
      flux_y(:,0)=state%dv%fhePSch(:,1)
      call mds_put_double_us('\TOP.SNAPSHOT:FHEPSCHY','W',
     1     nFc,flux_y,'_F')
      flux_x(:,0)=state%dv%fhe_eir(:,0)
      call mds_put_double_us('\TOP.SNAPSHOT:FHEX_EIR','W',
     1     nFc,flux_x,'_F')
      flux_y(:,0)=state%dv%fhe_eir(:,1)
      call mds_put_double_us('\TOP.SNAPSHOT:FHEY_EIR','W',
     1     nFc,flux_y,'_F')
      flux_x(:,0)=state%dv%fhi(:,0)
      call mds_put_double_us('\TOP.SNAPSHOT:FHIX','W',
     1     nFc,flux_x,'_F')
      flux_y(:,0)=state%dv%fhi(:,1)
      call mds_put_double_us('\TOP.SNAPSHOT:FHIY','W',
     1     nFc,flux_y,'_F')
      flux_x(:,0)=state%dv%fhi_mdf(:,0)
      call mds_put_double_us('\TOP.SNAPSHOT:FHIX_MDF','W',
     1     nFc,flux_x,'_F')
      flux_y(:,0)=state%dv%fhi_mdf(:,1)
      call mds_put_double_us('\TOP.SNAPSHOT:FHIY_MDF','W',
     1     nFc,flux_y,'_F')
      flux_x(:,0)=state%dv%fhiPSch(:,0)
      call mds_put_double_us('\TOP.SNAPSHOT:FHIPSCHX','W',
     1     nFc,flux_x,'_F')
      flux_y(:,0)=state%dv%fhiPSch(:,1)
      call mds_put_double_us('\TOP.SNAPSHOT:FHIPSCHY','W',
     1     nFc,flux_y,'_F')
      flux_x(:,0)=state%dv%fhi_eir(:,0)
      call mds_put_double_us('\TOP.SNAPSHOT:FHIX_EIR','W',
     1     nFc,flux_x,'_F')
      flux_y(:,0)=state%dv%fhi_eir(:,1)
      call mds_put_double_us('\TOP.SNAPSHOT:FHIY_EIR','W',
     1     nFc,flux_y,'_F')
      flux_x(:,0)=state%dv%fch(:,0)
      call mds_put_double_us('\TOP.SNAPSHOT:FCHX','A',
     1     nFc,flux_x,'_F')
      flux_y(:,0)=state%dv%fch(:,1)
      call mds_put_double_us('\TOP.SNAPSHOT:FCHY','A',
     1     nFc,flux_y,'_F')
      flux_x(:,0)=state%dv%fch_32(:,0)
      call mds_put_double_us('\TOP.SNAPSHOT:FCHX_32','A',
     1     nFc,flux_x,'_F')
      flux_y(:,0)=state%dv%fch_32(:,1)
      call mds_put_double_us('\TOP.SNAPSHOT:FCHY_32','A',
     1     nFc,flux_y,'_F')
      flux_x(:,0)=state%dv%fch_52(:,0)
      call mds_put_double_us('\TOP.SNAPSHOT:FCHX_52','A',
     1     nFc,flux_x,'_F')
      flux_y(:,0)=state%dv%fch_52(:,1)
      call mds_put_double_us('\TOP.SNAPSHOT:FCHY_52','A',
     1     nFc,flux_y,'_F')
      flux_x(:,0)=state%dv%fch_p(:,0)
      call mds_put_double_us('\TOP.SNAPSHOT:FCH_P','A',
     1     nFc,flux_x,'_F')
      flux_x(:,0)=state%dv%fchdia(:,0)
      call mds_put_double_us('\TOP.SNAPSHOT:FCHDIAX','A',
     1     nFc,flux_x,'_F')
      flux_y(:,0)=state%dv%fchdia(:,1)
      call mds_put_double_us('\TOP.SNAPSHOT:FCHDIAY','A',
     1     nFc,flux_y,'_F')
      flux_x(:,0)=state%dv%fchin(:,0)
      call mds_put_double_us('\TOP.SNAPSHOT:FCHINX','A',
     1     nFc,flux_x,'_F')
      flux_y(:,0)=state%dv%fchin(:,1)
      call mds_put_double_us('\TOP.SNAPSHOT:FCHINY','A',
     1     nFc,flux_y,'_F')
      flux_x(:,0)=state%dv%fchvispar(:,0)
      call mds_put_double_us('\TOP.SNAPSHOT:FCHVISPARX','A',
     1     nFc,flux_x,'_F')
      flux_y(:,0)=state%dv%fchvispar(:,1)
      call mds_put_double_us('\TOP.SNAPSHOT:FCHVISPARY','A',
     1     nFc,flux_y,'_F')
      flux_x(:,0)=state%dv%fchvisper(:,0)
      call mds_put_double_us('\TOP.SNAPSHOT:FCHVISPERX','A',
     1     nFc,flux_x,'_F')
      flux_y(:,0)=state%dv%fchvisper(:,1)
      call mds_put_double_us('\TOP.SNAPSHOT:FCHVISPERY','A',
     1     nFc,flux_y,'_F')
      flux_x(:,0)=state%dv%fchvisq(:,0)
      call mds_put_double_us('\TOP.SNAPSHOT:FCHVISQX','A',
     1     nFc,flux_x,'_F')
      flux_y(:,0)=state%dv%fchvisq(:,1)
      call mds_put_double_us('\TOP.SNAPSHOT:FCHVISQY','A',
     1     nFc,flux_y,'_F')
      flux_x(:,0)=state%dv%fchinert(:,0)
      call mds_put_double_us('\TOP.SNAPSHOT:FCHINERTX','A',
     1     nFc,flux_x,'_F')
      flux_y(:,0)=state%dv%fchinert(:,1)
      call mds_put_double_us('\TOP.SNAPSHOT:FCHINERTY','A',
     1     nFc,flux_y,'_F')
      flux_x(:,0)=state%dv%fht(:,0)
      call mds_put_double_us('\TOP.SNAPSHOT:FHTX','W',
     1     nFc,flux_x,'_F')
      flux_y(:,0)=state%dv%fht(:,1)
      call mds_put_double_us('\TOP.SNAPSHOT:FHTY','W',
     1     nFc,flux_y,'_F')
      flux_x(:,0)=state%dv%fhj(:,0)
      call mds_put_double_us('\TOP.SNAPSHOT:FHJX','W',
     1     nFc,flux_x,'_F')
      flux_y(:,0)=state%dv%fhj(:,1)
      call mds_put_double_us('\TOP.SNAPSHOT:FHJY','W',
     1     nFc,flux_y,'_F')
      flux_x(:,0)=state%dv%floe_noc(:,0)
      call mds_put_double_us('\TOP.SNAPSHOT:FLOEX_NOC','/s',
     1     nFc,flux_x,'_F')
      flux_y(:,0)=state%dv%floe_noc(:,1)
      call mds_put_double_us('\TOP.SNAPSHOT:FLOEY_NOC','/s',
     1     nFc,flux_y,'_F')
      flux_x(:,0)=state%dv%floi_noc(:,0)
      call mds_put_double_us('\TOP.SNAPSHOT:FLOIX_NOC','/s',
     1     nFc,flux_x,'_F')
      flux_y(:,0)=state%dv%floi_noc(:,1)
      call mds_put_double_us('\TOP.SNAPSHOT:FLOIY_NOC','/s',
     1     nFc,flux_y,'_F')

#ifdef WG_TODO
      flux_x(0:nx,-1:ny,0)=gs(0:nx,-1:ny,0)
      call mds_put_double_rz('\TOP.SNAPSHOT:SX','m^2',
     1 nx+1,ny+2,flux_x,'_X')
      do iy=-1,ny
        do ix=0,nx
          flux_x(ix,iy,0)=gs(ix,iy,0)*
     1     (vol(ix,iy)*qz(leftix(ix,iy),leftiy(ix,iy),1)+
     2     vol(leftix(ix,iy),leftiy(ix,iy))*qz(ix,iy,1))/
     3     (vol(leftix(ix,iy),leftiy(ix,iy))+vol(ix,iy))
        enddo
      enddo
      call mds_put_double_rz('\TOP.SNAPSHOT:SXPERP','m^2',
     1 nx+1,ny+2,flux_x,'_X')
      flux_y(-1:nx,0:ny,0)=gs(-1:nx,0:ny,1)
      call mds_put_double_rz('\TOP.SNAPSHOT:SY','m^2',
     1 nx+2,ny+1,flux_y,'_Y')

      tdi='build_with_units($,"K")'
      mds_signal_name='\TOP.SNAPSHOT:WALL_TEMP'
#if ( defined(GFORTRAN) && COMPILER_MAJOR_VERSION > 9 ) || defined(G95)
      status=mdsput2(
     * trim(mds_signal_name)//char(0),trim(tdi)//char(0),
     * descr2_i2(idtype_double,nbnd,ndepth_nml,0),target_temp,0)
#else
      status=mdsput2(
     * trim(mds_signal_name)//char(0),trim(tdi)//char(0),
     * descr2(idtype_double,nbnd,ndepth_nml,0),target_temp,0)
#endif
      write(*,*) 'mdsput2 of ',trim(mds_signal_name),' with ',trim(tdi)
      call mds_err_msg(status)

      call mds_put_double_rzdns('\TOP.SNAPSHOT:SMAF','N',
     1 nx+2,ny+2,4,ns,b2npmo_smaf,'TEXTLIN','TEXTPL','')
      call mds_put_double_rzdns('\TOP.SNAPSHOT:SMAG','N',
     1 nx+2,ny+2,4,ns,b2npmo_smag,'TEXTLIN','TEXTPL','')
      call mds_put_double_rzdns('\TOP.SNAPSHOT:SMAV','N',
     1 nx+2,ny+2,4,ns,b2npmo_smav,'TEXTLIN','TEXTPL','')
      call mds_put_double_rzns('\TOP.SNAPSHOT:SMFR','N',
     1 nx+2,ny+1,ns,smfr,'TEXTPL','')
#endif
      call mds_put_double_usdns('\TOP.SNAPSHOT:SMO','N',
     1     nCv,4,ns,state%sr%smo,'TEXTLIN','TEXTPL','')
#ifdef WG_TODO
      call mds_put_double_rzns('\TOP.SNAPSHOT:SMPR','N',
     1 nx+2,ny+1,ns,smpr,'TEXTPL','')
      call mds_put_double_rzns('\TOP.SNAPSHOT:SMPT','N',
     1 nx+2,ny+1,ns,smpt,'TEXTPL','')
#endif
      call mds_put_double_usdns('\TOP.SNAPSHOT:SMQ','N',
     1     nCv,4,ns,state%sr%smq,'TEXTLIN','TEXTPL','')
#ifdef WG_TODO
      tmp=smpr+smpt
      call mds_put_double_rzns('\TOP.SNAPSHOT:SMTH','N',
     1 nx+2,ny+1,ns,tmp,'TEXTPL','')
#endif
      call mds_put_double_usns('\TOP.SNAPSHOT:RESCO','#/s',
     1     nCv,ns,state%dv%resco,'TEXTPL','')
      call mds_put_double_us('\TOP.SNAPSHOT:RESHE','W',
     1     nCv,state%dv%reshe,'')
      call mds_put_double_us('\TOP.SNAPSHOT:RESHI','W',
     1     nCv,state%dv%reshi,'')
      call mds_put_double_usns('\TOP.SNAPSHOT:RESMO','N',
     1     nCv,ns,state%dv%resmo,'TEXTPL','')
      call mds_put_double_us('\TOP.SNAPSHOT:RESMT','N',
     1     nCv,state%dv%resmt,'')
      call mds_put_double_us('\TOP.SNAPSHOT:RESPO','V',
     1     nCv,state%dv%respo,'')

#ifdef B25_EIRENE
      if(natmi.gt.0) then
        forall (iCv=1:nCv,is=0:natmi-1) tmp(iCv,is)=dab2(iCv,is+1,1)
        call mds_put_double_usns('\TOP.SNAPSHOT:DAB2','/m^3',
     1       nCv,natmi,tmp,'TEXTAN','')
        forall (iCv=1:nCv,is=0:natmi-1) tmp(iCv,is)=tab2(iCv,is+1,1)/ev
        call mds_put_double_usns('\TOP.SNAPSHOT:TAB2','eV',
     1       nCv,natmi,tmp,'TEXTAN','')
        forall (iCv=1:nCv,is=0:natmi-1) tmp(iCv,is)=rfluxa(iCv,is+1,1)
        call mds_put_double_usns('\TOP.SNAPSHOT:RFLA','s^-1.m^-2',
     1       nCv,natmi,tmp,'TEXTAN','')
        forall (iCv=1:nCv,is=0:natmi-1) tmp(iCv,is)=refluxa(iCv,is+1,1)
        call mds_put_double_usns('\TOP.SNAPSHOT:REFA','W.m^-2',
     1       nCv,natmi,tmp,'TEXTAN','')
        forall (iCv=1:nCv,is=0:natmi-1) tmp(iCv,is)=pfluxa(iCv,is+1,1)
        call mds_put_double_usns('\TOP.SNAPSHOT:PFLA','s^-1.m^-2',
     1       nCv,natmi,tmp,'TEXTAN','')
        forall (iCv=1:nCv,is=0:natmi-1) tmp(iCv,is)=pefluxa(iCv,is+1,1)
        call mds_put_double_usns('\TOP.SNAPSHOT:PEFA','W.m^-2',
     1       nCv,natmi,tmp,'TEXTAN','')
        temp(:)=emiss(:,1,1)
        call mds_put_double_us('\TOP.SNAPSHOT:EMISS','ph.m^-3.s^-1',
     1       nCv,temp,'')
        forall (iCv=1:nCv,is=0:natmi-1) tmp(iCv,is)=eneutrad(iCv,is+1,0)
        call mds_put_double_usns('\TOP.SNAPSHOT:ENEUTRAD','W',
     1       nCv,natmi,tmp,'TEXTAN','')
      endif
      if(nmoli.gt.0) then
        forall (iCv=1:nCv,is=0:nmoli-1) tmp(iCv,is)=dmb2(iCv,is+1,1)
        call mds_put_double_usns('\TOP.SNAPSHOT:DMB2','/m^3',
     1       nCv,nmoli,tmp,'TEXTMN','')
        forall (iCv=1:nCv,is=0:nmoli-1) tmp(iCv,is)=tmb2(iCv,is+1,1)/ev
        call mds_put_double_usns('\TOP.SNAPSHOT:TMB2','eV',
     1       nCv,nmoli,tmp,'TEXTMN','')
        forall (iCv=1:nCv,is=0:nmoli-1) tmp(iCv,is)=rfluxm(iCv,is+1,1)
        call mds_put_double_usns('\TOP.SNAPSHOT:RFLM','s^-1.m^-2',
     1       nCv,nmoli,tmp,'TEXTMN','')
        forall (iCv=1:nCv,is=0:nmoli-1) tmp(iCv,is)=refluxm(iCV,is+1,1)
        call mds_put_double_usns('\TOP.SNAPSHOT:REFM','W.m^-2',
     1       nCv,nmoli,tmp,'TEXTMN','')
        forall (iCv=1:nCv,is=0:nmoli-1) tmp(iCv,is)=pfluxm(iCv,is+1,1)
        call mds_put_double_usns('\TOP.SNAPSHOT:PFLM','s^-1.m^-2',
     1       nCv,nmoli,tmp,'TEXTMN','')
        forall (iCv=1:nCv,is=0:nmoli-1) tmp(iCv,is)=pefluxm(iCv,is+1,1)
        call mds_put_double_usns('\TOP.SNAPSHOT:PEFM','W.m^-2',
     1       nCv,nmoli,tmp,'TEXTMN','')
        temp(:)=emissmol(:,1,1)
        call mds_put_double_us('\TOP.SNAPSHOT:EMISSMOL','ph.m^-3.s^-1',
     1       nCv,temp,'')
        if (jvft44.lt.20240627) then
          forall (iCv=1:nCv,is=0:nmoli-1) tmp(iCv,is)=srcml(iCv,is+1,1)
     1     *1.0e6_R8*geo%cvVol(iCv)/qe
        else
          forall (iCv=1:nCv,is=0:nmoli-1) tmp(iCv,is)=srcml(iCv,is+1,1)
        endif
        call mds_put_double_usns('\TOP.SNAPSHOT:SRCML','s^-1',
     1       nCv,nmoli,tmp,'TEXTMN','')
        if (jvft44.lt.20240627) then
          forall (iCv=1:nCv,is=0:nmoli-1) tmp(iCv,is)=edissml(iCv,is+1,1)
     1     *1.0e6_R8*geo%cvVol(iCv)
        else
          forall (iCv=1:nCv,is=0:nmoli-1) tmp(iCv,is)=edissml(iCv,is+1,1)
        endif
        call mds_put_double_usns('\TOP.SNAPSHOT:EDISSML','W',
     1       nCv,nmoli,tmp,'TEXTMN','')
        forall (iCv=1:nCv,is=0:nmoli-1) tmp(iCv,is)=emolrad(iCv,is+1,0)
        call mds_put_double_usns('\TOP.SNAPSHOT:EMOLRAD','W',
     1       nCv,nmoli,tmp,'TEXTMN','')
      endif
      if(nioni.gt.0) then
        forall (iCv=1:nCv,is=0:nioni-1) tmp(iCv,is)=dib2(iCv,is+1,1)
        call mds_put_double_usns('\TOP.SNAPSHOT:DIB2','/m^3',
     1       nCv,nioni,tmp,'TEXTIN','')
        forall (iCv=1:nCv,is=0:nioni-1) tmp(iCv,is)=tib2(iCv,is+1,1)/ev
        call mds_put_double_usns('\TOP.SNAPSHOT:TIB2','eV',
     1       nCv,nioni,tmp,'TEXTIN','')
        forall (iCv=1:nCv,is=0:nioni-1) tmp(iCv,is)=eionrad(iCv,is+1,0)
        call mds_put_double_usns('\TOP.SNAPSHOT:EIONRAD','W',
     1       nCv,nioni,tmp,'TEXTIN','')
      endif
#endif

      if (nimp.gt.0) then
        iCv1 = mpg%fcCv(ifsepimp,1)
        iCv2 = mpg%fcCv(ifsepimp,2) 
        rtmp=(state%pl%te(iCv1)/qe*gm%cvVol(iCv2)+
     2        state%pl%te(iCv2)/qe*gm%cvVol(iCv1))/
     3       (gm%cvVol(iCv1)+gm%cvVol(iCv2))
        call mds_put_double('\top.ident.impsep:te','eV',rtmp)
        dbimpsep_te=rtmp
        rtmp=(state%pl%ti(iCv1)/qe*gm%cvVol(iCv2)+
     2        state%pl%ti(iCv2)/qe*gm%cvVol(iCv1))/
     3       (gm%cvVol(iCv1)+gm%cvVol(iCv2))
        call mds_put_double('\top.ident.impsep:ti','eV',rtmp)
        dbimpsep_ti=rtmp
        rtmp=(state%dv%ne(iCv1)*gm%cvVol(iCv2)+
     2        state%dv%ne(iCv2)*gm%cvVol(iCv1))/
     3       (gm%cvVol(iCv1)+gm%cvVol(iCv2))
        call mds_put_double('\top.ident.impsep:ne','m^-3',rtmp)
        dbimpsep_ne=rtmp
        rtmp=(state%co%dna0(iCv1,ismain)*gm%cvVol(iCv2)+
     2        state%co%dna0(iCv2,ismain)*gm%cvVol(iCv1))/
     3       (gm%cvVol(iCv1)+gm%cvVol(iCv2))
        call mds_put_double('\top.ident.impsep:d','m^2/s',rtmp)
        dbimpsep_d=rtmp
        rtmp=(state%co%hce0(iCv1)/state%dv%ne(iCv1)*gm%cvVol(iCv2)+
     2        state%co%hce0(iCv2)/state%dv%ne(iCv2)*gm%cvVol(iCv1))/
     3       (gm%cvVol(iCv1)+gm%cvVol(iCv2))
        call mds_put_double('\top.ident.impsep:kye','m^2/s',rtmp)
        dbimpsep_kye=rtmp
        rtmp=(state%co%hcib(iCv1,ismain)/state%pl%na(iCv1,ismain)*
     1        gm%cvVol(iCv2)+
     2        state%co%hcib(iCv2,ismain)/state%pl%na(iCv2,ismain)*
     2        gm%cvVol(iCv1))/
     3       (gm%cvVol(iCv1)+gm%cvVol(iCv2))
        call mds_put_double('\top.ident.impsep:kyi','m^2/s',rtmp)
        dbimpsep_kyi=rtmp
        rtmp=(state%co%hci0(iCv1)/state%dv%ni(iCv1,0)*gm%cvVol(iCv2)+
     2        state%co%hci0(iCv2)/state%dv%ni(iCv2,0)*gm%cvVol(iCv1))/
     3       (gm%cvVol(iCv1)+gm%cvVol(iCv2))
        call mds_put_double('\top.ident.impsep:kyi0','m^2/s',rtmp)
        dbimpsep_kyi0=rtmp
        rtmp=crf(ifsepimp)
        call mds_put_double('\top.ident.impsep:cr','m',rtmp)
        dbimpsep_cr=rtmp
        rtmp=czf(ifsepimp)
        call mds_put_double('\top.ident.impsep:cz','m',rtmp)
        dbimpsep_cz=rtmp
      end if

      if (nomp.gt.0) then
        iCv1 = mpg%fcCv(ifsepomp,1)
        iCv2 = mpg%fcCv(ifsepomp,2) 
        rtmp=(state%pl%te(iCv1)/qe*gm%cvVol(iCv2)+
     2        state%pl%te(iCv2)/qe*gm%cvVol(iCv1))/
     3       (gm%cvVol(iCv1)+gm%cvVol(iCv2))
        call mds_put_double('\top.ident.ompsep:te','eV',rtmp)
        dbompsep_te=rtmp
        rtmp=(state%pl%ti(iCv1)/qe*gm%cvVol(iCv2)+
     2        state%pl%ti(iCv2)/qe*gm%cvVol(iCv1))/
     3       (gm%cvVol(iCv1)+gm%cvVol(iCv2))
        call mds_put_double('\top.ident.ompsep:ti','eV',rtmp)
        dbompsep_ti=rtmp
        rtmp=(state%dv%ne(iCv1)*gm%cvVol(iCv2)+
     2        state%dv%ne(iCv2)*gm%cvVol(iCv1))/
     3       (gm%cvVol(iCv1)+gm%cvVol(iCv2))
        call mds_put_double('\top.ident.ompsep:ne','eV',rtmp)
        dbompsep_ne=rtmp
        rtmp=(state%co%dna0(iCv1,ismain)*gm%cvVol(iCv2)+
     2        state%co%dna0(iCv2,ismain)*gm%cvVol(iCv1))/
     3       (gm%cvVol(iCv1)+gm%cvVol(iCv2))
        call mds_put_double('\top.ident.ompsep:d','m^2/s',rtmp)
        dbompsep_d=rtmp
        rtmp=(state%co%hce0(iCv1)/state%dv%ne(iCv1)*gm%cvVol(iCv2)+
     2        state%co%hce0(iCv2)/state%dv%ne(iCv2)*gm%cvVol(iCv1))/
     3       (gm%cvVol(iCv1)+gm%cvVol(iCv2))
        call mds_put_double('\top.ident.ompsep:kye','m^2/s',rtmp)
        dbompsep_kye=rtmp
        rtmp=(state%co%hcib(iCv1,ismain)/state%pl%na(iCv1,ismain)*
     1        gm%cvVol(iCv2)+
     2        state%co%hcib(iCv2,ismain)/state%pl%na(iCv2,ismain)*
     2        gm%cvVol(iCv1))/
     3       (gm%cvVol(iCv1)+gm%cvVol(iCv2))
        call mds_put_double('\top.ident.ompsep:kyi','m^2/s',rtmp)
        dbompsep_kyi=rtmp
        rtmp=(state%co%hci0(iCv1)/state%dv%ni(iCv1,0)*gm%cvVol(iCv2)+
     2        state%co%hci0(iCv2)/state%dv%ni(iCv2,0)*gm%cvVol(iCv1))/
     3       (gm%cvVol(iCv1)+gm%cvVol(iCv2))
        call mds_put_double('\top.ident.ompsep:kyi0','m^2/s',rtmp)
        dbompsep_kyi0=rtmp
        rtmp=crf(ifsepomp)
        call mds_put_double('\top.ident.ompsep:cr','m',rtmp)
        dbompsep_cr=rtmp
        rtmp=czf(ifsepomp)
        call mds_put_double('\top.ident.ompsep:cz','m',rtmp)
        dbompsep_cz=rtmp
      end if

      fhesep=0.0_R8
      fhisep=0.0_R8
      fnesep=0.0_R8
      fnisep=0.0_R8
      if (mpg%iFssep.ne.0) then
        do i = mpg%fsFcP(mpg%iFssep,1),
     1         mpg%fsFcP(mpg%iFssep,1) + mpg%fsFcP(mpg%iFssep,2) - 1
          iCv1 = mpg%fcCv(mpg%fsFc(i),1)
          iCv2 = mpg%fcCv(mpg%fsFc(i),2)
          if(mpg%cvOnClosedSurface(iCv1).neqv.
     1       mpg%cvOnClosedSurface(iCv2)) then
            fhesep=fhesep+state%dv%fhe(mpg%fsFc(i),1)
            fhisep=fhisep+state%dv%fhi(mpg%fsFc(i),1)
            fnesep=fnesep+state%dv%fne(mpg%fsFc(i),1)
            fnisep=fnisep+state%dv%fni(mpg%fsFc(i),1)
          endif
        enddo
      endif
      if (mpg%iFssep2.ne.0) then
        do i = mpg%fsFcP(mpg%iFssep2,1),
     1         mpg%fsFcP(mpg%iFssep2,1) + mpg%fsFcP(mpg%iFssep2,2) - 1
          iCv1 = mpg%fcCv(mpg%fsFc(i),1)
          iCv2 = mpg%fcCv(mpg%fsFc(i),2)
          if(mpg%cvOnClosedSurface(iCv1).neqv.
     1       mpg%cvOnClosedSurface(iCv2)) then
            fhesep=fhesep+state%dv%fhe(mpg%fsFc(i),1)
            fhisep=fhisep+state%dv%fhi(mpg%fsFc(i),1)
            fnesep=fnesep+state%dv%fne(mpg%fsFc(i),1)
            fnisep=fnisep+state%dv%fni(mpg%fsFc(i),1)
          endif
        enddo
      endif
      call mds_put_double('\top.ident.sep:fhe','W',fhesep)
      dbsep_fhe=fhesep
      call mds_put_double('\top.ident.sep:fhi','W',fhisep)
      dbsep_fhi=fhisep
      call mds_put_double('\top.ident.sep:fne','#/s',fnesep)
      dbsep_fne=fnesep
      call mds_put_double('\top.ident.sep:fni','#/s',fnisep)
      dbsep_fni=fnisep

      write(*,*) '... b2mds'

#else
      call xerrab('b2mds: ' //
     1 'MDSPLUS routine called when not compiled with MDSPLUS')
#endif
      call dealloc_b2mod_wall
      call dealloc_b2mod_layer
      call dealloc_b2mod_b2plot
      call dealloc_sputter_data
      call dealloc_b2mod_sputter
      call dealloc_b2mod_elements
      call dealloc_b2mod_boundary
      call dealloc_b2mod_neutrals
      call dealloc_b2mod_wall_init
#ifdef B25_EIRENE
      CALL EIRENE_DEALLOC_CTRCEI
      CALL EIRENE_DEALLOC_BRAEIR
      CALL EIRENE_DEALLOC_EIRBRA
      CALL EIRENE_DEALLOC_INFCOP
      CALL EIRENE_EXTRAB25_CLEANUP
      CALL EIRENE_WNEUTRALS_DEALLOC
#else
      call dealloc_b2mod_braeir
      call dealloc_b2mod_eirbra
      call dealloc_b2mod_eirdiag
#endif
      call subend()
      return
      end subroutine b2mds
