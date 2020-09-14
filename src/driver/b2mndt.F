*-----------------------------------------------------------------------
*.specification

      subroutine b2mndt (nout, nCv, nFc, nVx, ns, nxtl, nxtr,
     . ismain, ismain0, nscx, nscxmax, iscx, itim, dtim, ntim, switch, 
     . geo, mpg, st, st_ext, ierr)
      use b2mod_types
!      use b2mod_anomalous_transport
!      use b2mod_plasma
!      use b2mod_rates
!      use b2mod_residuals
!      use b2mod_sources
!      use b2mod_transport
      use b2mod_constants
      use b2mod_diag
     &   , only : nnatmi, alloc_b2mod_diag
!      use b2mod_work
      use b2mod_time
      use b2mod_b2cmpa
!      use b2mod_neoclassical
      use b2mod_user_namelist
      use b2mod_boundary_namelist
      use b2mod_numerics_namelist
      use b2mod_equation_sources                                        !srv 26.02.18
      use b2mod_boundary_sources                                        !srv 26.02.18
      use b2mod_switches
      use b2us_geo
      use b2us_map
      use b2us_plasma
      use b2mod_running_average
      use b2mod_batch_average
     & , only : batch_av_all
      implicit none
*   ..input arguments (unchanged on exit)
      integer nout(0:10), nCv, nFc, nVx, ns, nxtl, nxtr, ismain,
     *  ismain0, nscx, nscxmax, iscx(0:nscxmax-1), itim, ntim
      integer, save :: ntim_batch = 500
      integer, save :: ibatch_av_all = 0
      real (kind=R8) ::
     *  dtim
      type (geometry), intent (in) :: geo
      type (mapping), intent (in) :: mpg
      type (switches), intent (in) :: switch
*   ..input/output arguments
      type (B2State), intent (inout) :: st
      type (B2StateExt), intent (inout) :: st_ext
*   ..output arguments (unspecified on entry)
      integer ierr
*   ..common blocks

*-----------------------------------------------------------------------
*.documentation
*
*  1. purpose
*
*     B2MNDT performs one implicit timestep for the b2 system of
*     equations.
*
*
*  3. description (see also routine b2cdca)
*
*     For a detailed description of the equations and the method of
*     solution, see the references cited in the main program unit.
*
*     This routine performs one timestep for the b2 equations.
*     An outline of the procedure follows.
*     ..test input arguments
*     ..compute auxiliary quantities
*     ..prepare source computation
*     ..do while ()
*      ..compute log-log linearised rate coefficients
*      ..do while ()
*       ..compute source linearisation
*       ..do while ()
*        ..perform one inner iteration
*        ..re-compute auxiliary quantities
*        ..produce monitoring output
*       ..enddo
*      ..enddo
*     ..enddo
*
*
*-----------------------------------------------------------------------
*.declarations

*   ..local variables
      integer ncall, nstg(0:2), moqtlv, moitlv, itcnt, density_control,  !xpb
     *  ramp_slow, ier0, no_solve, iCv, boundary_namelist,
     *  nsmin, nsmax
      integer, save :: ntim_step_out = 1                                 !srv 31.10.07
      integer, save :: use_b2srst = 1                                    !srv 13.07.12
      integer, save :: ion_vlct_restrict = 0, tn_style = 0                             !iyv 06.02.14
      integer, save :: boundary_sources=0, equation_sources=0            !srv 26.02.18
      real (kind=R8) ::
     *  rxf, BoRiS
      real (kind=R8) ::
     *  nstg_areshe, nstg_areshi, nstg_aresco
      real (kind=R8), save :: facdrift_start, facExB_start, facvis_start
      save ncall, nstg, moqtlv, moitlv, rxf, itcnt, nsmin, nsmax,        !xpb
     *  density_control, ramp_slow, BoRiS, no_solve,
     *  nstg_areshe, nstg_areshi, nstg_aresco, boundary_namelist
      logical domoqt, domoit, wrong_flow, ramp_call                      !xpb
      logical test_residual, quit_residual, main_call
      logical solving(4)
      save test_residual
      integer istg(0:2), is, k
      real (kind=R8) ::
     *  acorpa(0:ns-1), acorua(0:ns-1), acorut,
     &  acorte, acorti, acorpo, nasum, na0sum, 
     &  pz(nCv), rz(nCv), cssb(nCv,0:1)         !iyv 06.02.14
! csc UNCOMMENT WITH AD
!  REAL(kind=r8) :: acorpad(0:ns-1), acoruad(0:ns-1), acorutd, acorted, &
!& acortid, acorpod
!WG_TODO: add some sources locally here, until a module/data-type for WG
!WG_TODO: is available
      real (kind=R8) ::                                                  !xpb
     &  snasm(nCv,0:1,0:ns-1),                                    !xpb
     &  smosm(nCv,0:3,0:ns-1),                                    !xpb
     &  shesm(nCv,0:3), shism(nCv,0:3)!WG_TODO,                    !xpb
!WG_TODO     &  snadt(-1:nx,-1:ny,0:1,0:ns-1),                                    !xpb
!WG_TODO     &  smodt(-1:nx,-1:ny,0:3,0:ns-1),                                    !xpb
!WG_TODO     &  shedt(-1:nx,-1:ny,0:3), shidt(-1:nx,-1:ny,0:3),                    !xpb
!WG_TODO     &  schdt(-1:nx,-1:ny,0:3),                                           !xpb
!WG_TODO     &  snedt(-1:nx,-1:ny,0:1)
*   ..procedures
      real (kind=R8) ::
     *  smin, smax, b2sasum, samax
      external subini, subend, xertst, ipgeti, ipgetr, smin,
     &  smax, b2sasum, samax
      external b2xvsg, b2xvps, b2xpne, b2xpni,
     &  b2spel, b2spcx, b2sral, b2srdt, b2srst, b2news, b2mxzr,
     &  b2mxar, b2mxac, b2mwqt, b2mwit, b2tral, b2srsm
*   ..initialisation
      data ncall/0/, nstg/1,1,1/, moqtlv/3/, moitlv/-1/, rxf/0.5_R8/
      data itcnt/0/, density_control/0/, ramp_slow/0/, BoRiS/0.0_R8/
      data no_solve/0/, boundary_namelist/0/
      data nstg_areshe/0.0_R8/, nstg_areshi/0.0_R8/, nstg_aresco/0.0_R8/
*
*-----------------------------------------------------------------------
*.computation

* ..preliminaries
*   ..subprogram start-up calls
      call subini ('b2mndt')
*   ..set internal parameters on first call
      if (ncall.eq.0) then
       nsmin=0
       nsmax=ns
       if (ntim.ne.0) then
         call ipgeti ('b2mndt_nstg0', nstg(0))
         call ipgeti ('b2mndt_nstg1', nstg(1))
         call ipgeti ('b2mndt_nstg2', nstg(2))
         call xertst (1.le.nstg(0).and.1.le.nstg(1).and.1.le.nstg(2),
     &    'faulty internal parameter nstg(0:2)')
       else
         nstg(0) = 1
         nstg(1) = 1
         nstg(2) = 0
       end if
       call ipgeti ('b2mndt_moqtlv', moqtlv)
       call xertst (-1.le.moqtlv.and.moqtlv.le.3,
     &   'faulty internal parameter moqtlv')
       call ipgeti ('b2mndt_moitlv', moitlv)
       call xertst (-1.le.moitlv.and.moitlv.le.3,
     &   'faulty internal parameter moitlv')
       call ipgetr ('b2mndt_rxf', rxf)
       call xertst (0.0_R8.le.rxf.and.rxf.le.1.0_R8,
     &   'faulty internal parameter rxf')
       call ipgeti ('b2mndt_density_control',density_control)
       call ipgeti ('b2news_ramp_slow', ramp_slow)
       call ipgetr ('b2news_BoRiS', BoRiS)
       call ipgeti ('b2news_no_solve', no_solve)
       call xertst (-15.le.no_solve.and.no_solve.le.1,
     &   'faulty input no_solve')
       if (no_solve.eq.1) nstg(0:2)=1
       call ipgetr ('b2mndt_nstg_areshe', nstg_areshe)
       call xertst (nstg_areshe.ge.0.0_R8,
     &      'faulty internal parameter nstg_areshe')
       call ipgetr ('b2mndt_nstg_areshi', nstg_areshi)
       call xertst (nstg_areshi.ge.0.0_R8,
     &      'faulty internal parameter nstg_areshi')
       call ipgetr ('b2mndt_nstg_aresco', nstg_aresco)
       call xertst (nstg_aresco.ge.0.0_R8,
     &      'faulty internal parameter nstg_aresco')
       test_residual=nstg_areshe.gt.0.0_R8.or.nstg_areshi.gt.0.0_R8.or.
     &               nstg_aresco.gt.0.0_R8
       call ipgeti ('b2npmo_ion_vlct_restrict', ion_vlct_restrict)      !iyv 06.02.14
       if (switch%use_eirene) then
         call alloc_b2mod_diag(nCv,nFc,ns,mpg%nnreg(0),nstrat,nnatmi)
       endif
       call ipgeti ('b2mndt_ntim_step_out', ntim_step_out)              !srv 31.01.07 {
       call xertst (1.le.ntim_step_out,
     &   'faulty internal parameter ntim_step_out')                     !srv 31.01.07 }
       call ipgeti ('b2mndt_use_b2srst', use_b2srst)                    !srv 13.07.12
       call ipgeti ('b2mndr_equation_sources', equation_sources)        !srv 26.02.18 {
       if(equation_sources.ne.0) call read_equation_sources(nCv, ns)
       call ipgeti ('b2mndr_boundary_sources', boundary_sources)
       if(boundary_sources.ne.0) call read_boundary_sources(nCv, ns)    !srv 26.02.18 }
       call ipgeti ('b2stbc_boundary_namelist', boundary_namelist)
       call ipgeti ('b2news_nsmin', nsmin)
       if(nsmin.gt.ns) then
         nsmin=ns
       elseif(nsmin.lt.0) then
         nsmin=0
       endif
       call ipgeti ('b2news_nsmax', nsmax)
       if(nsmax.gt.ns) then
         nsmax=ns
       elseif(nsmax.lt.0) then
         nsmax=0
       endif
       call ipgeti ('b2mndt_av_ntim_batch', ntim_batch)
       call ipgeti ('b2mndt_av_batch_all', ibatch_av_all)
       call ipgeti ('b2mn_tn_style', tn_style)
       if (boundary_namelist.ge.1) then
!WG_TODO         call initialize_b2mod_boundary(ns,mpg,tim)
       else
         call xerrab('b2mndr -- specify BCs using boundary_namelist')
!WG_TODO         call alloc_b2mod_boundary(nx,ny,ns,nbcd,nnisod)
       endif
      else
       if (boundary_namelist.ge.1) then
!WG_TODO         call catchup_b2mod_boundary(ns,mpg,tim)
       endif
      endif
*   ..test nCv, nFc, ns, ismain
      call xertst (0.le.nCv.and.0.le.nFc, 'faulty argument nCv, nFc')
      call xertst (1.le.ns, 'faulty argument ns')
      call xertst (0.le.ismain.and.ismain.lt.ns,
     . 'invalid main plasma species index ismain')
      call xertst (.not.is_neutral(ismain),
     . 'invalid main plasma species ismain; must not be neutral')
      call xertst (0.le.ismain0.and.ismain0.lt.ns,
     . 'invalid main neutral species index ismain0')
      call xertst (is_neutral(ismain0).or.ismain.eq.ismain0,
     . 'invalid main neutral species ismain0; must be neutral')
*   ..test itim, dtim
      call xertst (0.le.itim.and.0.0_R8.lt.dtim,
     &  'faulty argument itim, dtim')
*   ..test sign of bb
      call b2xvsg (nCv, geo%cvBb(1,3), 1, 'bb3', '.gt.')
*   ..test sign of vol, hx, hy
      call b2xvsg (nCv, geo%cvVol, 1, 'vol', '.gt.')
*   ..test old state
      call b2xvsg (nCv*ns, st%psnc%na, 1, 'na0', '.gt.')
      call b2xvsg (nCv, st%psnc%te, 1, 'te0', '.gt.')
      call b2xvsg (nCv, st%psnc%ti, 1, 'ti0', '.gt.')
*   ..test present state
      call b2xvps (nCv, nFc, ns, st%pl, st%dv)
*   ..test facdrift and fac_ExB (not on first call)                      !xpb
      if (ncall.gt.0) then                                               !xpb
        call xertst (0.0_R8.le.minval(st%dv%facdrift) .and.
     &    1.0_R8.ge.maxval(st%dv%facdrift), 'faulty argument facdrift')
        call xertst (0.0_R8.le.minval(st%dv%fac_ExB) .and.
     &    1.0_R8.ge.maxval(st%dv%fac_ExB), 'faulty argument fac_ExB')
      else
c the main initialisation is in b2news --- unfortunately a few routines
c that need an initial value get called before the first call to b2news
c therefore initialise here
cwdk    For unstructured solver, start with simple reading of start-values
        call ipgetr ('b2news_facExB_start',   facExB_start)
        call ipgetr ('b2news_facdrift_start', facdrift_start)
        call ipgetr ('b2news_facvis_start',   facvis_start)
        st%dv%fac_ExB  = facExB_start
        st%dv%facdrift = facdrift_start
        st%dv%fac_vis  = facvis_start

c checks for use of new energy equations
       if (tn_style.ne.0) then
         call xertst(ns.eq.2,
     &        'new energy equations only for single species')
       endif
#ifdef WG_TODO
        call b2news_facdrift_fac_ExB_init(nx,ny,ramp_slow)
#endif
      endif
*   ..initialise error parameter
      ierr = 0
      quit_residual = .false.

* ..initialisation
*   ..initialise sources to 0
      call sfill (nCv*2*ns, 0.0_R8, snasm, 1)                        !srv 11.09.09 {
      call sfill (nCv*4*ns, 0.0_R8, smosm, 1)
      call sfill (nCv*4, 0.0_R8, shesm, 1)
      call sfill (nCv*4, 0.0_R8, shism, 1)
!WG_TODO      call sfill (n2*2*ns, 0.0_R8, snadt, 1)
!WG_TODO      call sfill (n2*4*ns, 0.0_R8, smodt, 1)
!WG_TODO      call sfill (n2*4, 0.0_R8, shedt, 1)
!WG_TODO      call sfill (n2*4, 0.0_R8, shidt, 1)
!WG_TODO      call sfill (n2*4, 0.0_R8, schdt, 1)
!WG_TODO      call sfill (n2*2, 0.0_R8, snedt, 1)
*   ..compute ni0
      call b2xpni (nCv, ns, st%psnc%na, st%psnc%ni)
*   ..compute nasum, na0sum
      nasum = b2sasum(nCv*ns, st%pl%na, 1)
      na0sum = b2sasum(nCv*ns, st%psnc%na, 1)
*   ..rescale densities
      if (density_control.ne.0) then
        st%dv%ne = st%dv%ne*na0sum/nasum
        st%pl%na = st%pl%na*na0sum/nasum
        st%dv%fna = st%dv%fna*na0sum/nasum
        st%dv%kinrgy = st%dv%kinrgy*na0sum/nasum
        write(*,*) 'Densities rescaled by :',na0sum/nasum
      endif
*   ..compute initial ni
      call b2xpni (nCv, ns, st%pl%na, st%dv%ni)
* ..main loop
*   ..outer iteration
      istg(0) = 0
    1 if (istg(0).lt.nstg(0) .and. .not.quit_residual) then
*    ..compute log-log electron rate coefficients
        call b2spel (nCv, ns, ev, st%pl%te, st%dv%ne, st%rt)
*    ..compute log-log charge exchange rate coefficients
        do k = 0, nscx-1
         call b2spcx (nCv, ns, ev, am(iscx(k)), st%pl%ti, st%dv%ne,
     &     st%rt%rlcx(1,0,0,k))
        enddo
ctmp.dpc
        if (nscx.ge.1) then
         write(*,*) 'DPC: b2mndt: rlcx'
         do k = 0, nscx-1
           do is=0, ns-1
             write(*,'(a,2i4,1p,2g14.7)') 'is0,is,rlcx_0,rlcx_1 ',
     1            k,is,samax(nCv,st%rt%rlcx(1,0,is,k),1),
     2                 samax(nCv,st%rt%rlcx(1,1,is,k),1)
           enddo
         enddo
        else
         write(*,*) 'DPC: b2mndt: no CX reactions'
        endif
ctmp.dpc
*    ..intermediate iteration
        istg(1) = 0
    2   if (istg(1).lt.nstg(1) .and. .not.quit_residual) then
*   ..calculate neoclassical transport only in intermediate loop
!WG_TODO          newneo=.true.

*   ..one-off calculation of transport coefficients
          call b2xpni (nCv, ns, st%pl%na, st%dv%ni)
          call b2xpne (nCv, ns, st%rt%rza, st%pl%na, st_ext%ne, 
     .      st%dv%ne)
          call b2xpne (nCv, ns, st%rt%rz2, st%pl%na, st_ext%ne2, 
     .      st%dv%ne2)
*   ..compute flux limit coefficients                                    !srv 11.10.02 03.06.03 {
          call b2trql (nCv, nFc, ns, switch, geo, mpg, st%pl, st%dv, 
     .      st_ext, st%co%chvemx, st%co%chvimx)
          call b2tral (nCv, nFc, nVx, ns, nscx, nscxmax, iscx, ismain,   !srv 23.11.10
     1      switch, geo, mpg, st%pl, st%dv, st%rt, st_ext, st%co)
!WG_TODO         newneo=.true.
           call b2tanml (nCv, nFc, ns, switch, geo, mpg, st%co%csig_an,
     &                    st%pl%po, st%dv%fchanml)
          call b2tvspa (nCv, nFc, nVx, ns, switch, geo, mpg, st%pl%ua,
     &       st%co%vsaf_cl, st%dv%fac_vis, st%dv%fchvispar)
          call b2tiner (nCv, nFc, nVx, ns, switch, geo, mpg, st%pl%na,
     &       st%pl%ua, st%dv%facdrift, st%dv%fchinert)
#ifdef WG_TODO
          src_ni=0.0_R8
          src_ne=0.0_R8
          src_mo=0.0_R8
          src_ee=0.0_R8
          src_ei=0.0_R8
          src_ms=0.0_R8
          src_es=0.0_R8
          cor_ni=0.0_R8
          cor_ne=0.0_R8
          cor_mo=0.0_R8
          cor_ee=0.0_R8
          cor_ei=0.0_R8
#endif
*     ..compute source linearisation
          main_call = ntim.ne.0
          call b2sral (nCv, nFc, nVx, ns, nxtl, nxtr, nscx, 
     &      nscxmax, nsmin, nsmax, iscx, ismain, ismain0, dtim, BoRiS,
     &      switch, geo, mpg, st%pl, st%dv, st%co, st%rt, st%rtw, 
     &      st_ext, st%sr, st%srw, wrong_flow, main_call)
          if(wrong_flow) write(*,*)
     &         'b2mndt: wrong_flow returned from b2sral'
*     ..smooth out sources
!WG_TODO          call b2srsm (nx, ny, ns,
!WG_TODO     &     dtim, vol, na, ua, te, ti, ne, ni,
!WG_TODO     &     sna, smo, she, shi, sne, .true.)                              !srv 11.09.09 {
*     ..include contributions due to timestep
          call b2srdt (nCv, ns, dtim, BoRiS, switch, geo, mpg,
     &    st%psnc%na, st%psnc%ua, st%psnc%te, st%psnc%ti, 
     &    st%psnc%tn,
     &    st%psnc%ne, st%psnc%ni, st%psnc%kinrgy,
     &    st%pl%na, st%pl%ua, st%pl%te, st%pl%ti, st%dv%ne, 
     &    st%dv%ni, st%dv%kinrgy,  st%sr, .false.)
*   ..add to previous contributions
!WG_TODO This is already done in b2srdt
!WG_TODO Still need *dt-arrays for output?
!WG_TODO          call b2saxpy (n2*2*ns, 1.0_R8, snadt, 1, sna, 1)
!WG_TODO          call b2saxpy (n2*4*ns, 1.0_R8, smodt, 1, smo, 1)
!WG_TODO          call b2saxpy (n2*2, 1.0_R8, snedt, 1, sne, 1)
!WG_TODO          call b2saxpy (n2*4, 1.0_R8, shedt, 1, she, 1)
!WG_TODO          call b2saxpy (n2*4, 1.0_R8, shidt, 1, shi, 1)
!WG_TODO          call b2saxpy (n2*4, 1.0_R8, schdt, 1, sch, 1)                  !srv 11.09.09 }
          if(equation_sources.ne.0) then                                 !srv 26.02.18 {
	    do iCv = 1,nCv
              st%sr%sna(iCv,0:1,0:ns-1) = st%sr%sna(iCv,0:1,0:ns-1)+
     &          geo%cvVol(iCv)*art_sna(iCv,0:1,0:ns-1) 
              st%sr%smo(iCv,0:3,0:ns-1) = st%sr%smo(iCv,0:3,0:ns-1)+
     &          geo%cvVol(iCv)*geo%cvHz(iCv)*art_smo(iCv,0:3,0:ns-1) 
              st%sr%she(iCv,0:3) = st%sr%she(iCv,0:3)+
     &          geo%cvVol(iCv)*art_she(iCv,0:3) 
              st%sr%shi(iCv,0:3) = st%sr%shi(iCv,0:3)+
     &          geo%cvVol(iCv)*art_shi(iCv,0:3) 
              st%sr%sch(iCv,0:3) = st%sr%sch(iCv,0:3)+
     &          geo%cvVol(iCv)*art_sch(iCv,0:3) 
            enddo
          endif                                                          !srv 26.02.18 }
*     ..stabilise the source coefficients
          if(use_b2srst.eq.1)
     &     call b2srst (nCv, ns, switch, mpg, st%pl%na, st%pl%ua,    !srv 01.07.08
     &      st%pl%te, st%pl%ti, st%pl%tn, st%pl%po, st%dv%ne, st%dv%ni,
     &      st%sr)
*       (smq is inherently stable.)
*     ..inner iteration
          istg(2) = 0
    3     if (istg(2).lt.nstg(2) .and. .not.quit_residual) then
*      ..basic iterative step
#ifdef WG_TODO
           ramp_call = ramp_slow.eq.0.or.(istg(0)+istg(1)+istg(2)).eq.0
#endif
           call b2news_(nCv, nFc, nVx, ns, nxtl, nxtr, nscx, iscx,
     &       nscxmax, ismain, ismain0, dtim, rxf, switch, geo, mpg,
     &       st%pl, st%dv, st%co, st%rt, st%rtw, st%sr, st%srw, 
     &       st_ext, st%psnc, ier0)
           flush(6)
           if(ier0.eq.10) then
             ierr=ier0
             goto 99999
           endif
           call xertst (ier0.eq.0, 'error return from b2news')
*!!      (should develop other error handling)
           ierr = ier0
*      ..re-compute ne, ne2, ni
           call b2xpne (nCv, ns, st%rt%rza, st%pl%na, st_ext%ne, 
     .       st%dv%ne)
           call b2xpne (nCv, ns, st%rt%rz2, st%pl%na, st_ext%ne2, 
     .       st%dv%ne2)
           call b2xpni (nCv, ns, st%pl%na, st%dv%ni)
*      ..zero out residuals on special regions
*        (this in order to avoid their contribution to the norms)
           call b2mxzr (nCv, ns, mpg, st%dv%resco, st%dv%resmo, 
     &      st%dv%resmt, st%dv%reshe, st%dv%reshi, st%dv%respo)
*      ..compute norms of the residuals
           call b2mxar (nCv, ns, nnreg(0), BoRiS, switch, geo, mpg, 
     &      st%pl, st%dv, st%diag)
           quit_residual=test_residual .and.
     &      ((st%diag%areshe*nstg_areshe).le.(nstg_areshe*nstg_areshe)
     &       .or.nstg_areshe.eq.0.0_R8).and.
     &      ((st%diag%areshi*nstg_areshi).le.(nstg_areshi*nstg_areshi)
     &       .or.nstg_areshi.eq.0.0_R8).and.
     &      ((st%diag%aresco(ismain)*nstg_aresco)
     &       .le.(nstg_aresco*nstg_aresco)
     &       .or.nstg_aresco.eq.0.0_R8)
*      ..compute norms of the corrections
           call b2mxac (nCv, ns, st%dv, st%diag)
! csc UNCOMMENT WITH AD
!!*      ..compute norms of the differentiated residuals
!          CALL B2MXAR_DIFF(ncv, ns, nnreg(0), boris, geo, mpg, std%pl, &
!&                      std%dv, std%diag)
!!*      ..compute norms of the differentiated corrections
!          CALL B2MXAC_DIFF(ncv, ns, std%dv, std%diag)
*      ..decide what output
           domoqt = 0.le.moqtlv.and.
     &      (1.le.moqtlv.or.istg(0).eq.0).and.
     &      (2.le.moqtlv.or.istg(1).eq.0).and.
     &      (3.le.moqtlv.or.istg(2).eq.0)
           domoit = 0.le.moitlv.and.
     &      (1.le.moitlv.or.istg(0).eq.0).and.
     &      (2.le.moitlv.or.istg(1).eq.0).and.
     &      (3.le.moitlv.or.istg(2).eq.0)
*      ..maybe write quick trace
           if (domoqt) then
            if (mod(ncall,ntim_step_out) .eq. 0) then                    !srv 31.01.07
             call b2mwqt (nout(4), nCv, ns, itim, itcnt, BoRiS, 
     &        switch, geo, st%pl, st%dv, st%diag)
             flush(nout(4))
! csc UNCOMMENT WITH AD
!              CALL B2MWQT_DIFF(nout(10), ncv, ns, itim, itcnt, boris, &
!&                          geo, std%pl, std%dv, std%diag)
!              FLUSH(nout(10)) 
            endif
           endif
#ifdef WG_TODO
*      ..maybe detailed monitor
           if (domoit) then
             call b2mwit (nout(5), nx, ny, ns, itim, itcnt, resco,
     &        resmo, resmt, reshe, reshi, respo,
     &        corpa, corua, corut, corte, corti,
     &        corpo)
             flush(nout(5))
           endif
*      ..write the tracing data
           if (ank_tracing.ne.0) then
             if (ncall.eq.0.and.istg(2).eq.0) then
               call alloc_b2mod_user(nx,ns,nlim,nmol,ntns)
               call b2filei(1,ns,nnreg)
             end if
             if (mod(itim,ank_tracing).eq.1.or.ank_tracing.eq.1)
     >        call b2filew(istg(2),ns,nnreg)
           endif
#endif
*      ..update iteration counter
           itcnt = itcnt+1
*      ..end inner iteration
           istg(2) = istg(2)+1
           goto 3
         else if (ntim.eq.0) then ! compute derived quantities
#ifdef WG_TODO
           if(ion_vlct_restrict .ne. 0 ) then                            !iyv 06.02.14
             call b2xprz (nx, ny, ns, mp, am, na, rz)
             call b2xppz (nx, ny, ns, ne, na, te, ti, pz)
             wrk0(:,:)=sqrt(pz(:,:)/rz(:,:))*abs(bb(:,:,0)/bb(:,:,3))
             call intfaceh(nx,ny,leftix,leftiy,vol,wrk0,cssb(-1,-1,0))
             call intfacev(nx,ny,bottomix,bottomiy,vol,wrk0,
     &                     cssb(-1,-1,1))
           else
              cssb = 0.0_R8
           endif                                                         !iyv 06.02.14
           do is=0,ns-1
*      ..compute partial pressure
             call b2xppb (nx, ny, rza(-1,-1,is), na(-1,-1,is),
     &                    te, ti, pb)                                    !srv 01.10.99
*      ..compute particle fluxes
             call b2tfnb (nx, ny, is, ismain,
     &          guard_flows, facdrift, fac_ExB, am(is)*mp, vol, cssb,    !srv 07.06.00 09.01.01 08.09.01 19.05.03 04.02.05
     &          rza(-1,-1,is), pbs, pbshz, na(-1,-1,is),                 !srv 09.01.01 20.09.06
     &          pb, ua(-1,-1,is), uadia(-1,-1,0,0,is), 
     &          vaecrb(-1,-1,0,0,is), vadia(-1,-1,0,0,is),
     &          wadia(-1,-1,0,0,is), ne, te, ti, bb,
     &          hx, hy, hz, qz, gs, qc, qs, po,                          !srv 19.05.03 08.10.04 04.02.05
     &          cvla(-1,-1,0,0,is), cdna(-1,-1,0,0,is), 
     &          cdpa(-1,-1,0,0,is),
     &          cddi(-1,-1,0,0,is), cvlahz(-1,-1,0,0,is),
     &          cvmahz(-1,-1,0,0,is), cdnahz(-1,-1,0,0,is),              !srv 09.01.01 15.12.05 13.10.06
     &          cdpahz(-1,-1,0,0,is),
     &          fchvispar, fchvisper, fchinert, fchanml, fchvisq,        !srv 20.09.06
     &          fna(-1,-1,0,0,is), fna_nodrift(-1,-1,0,0,is),
     &          fna_fcor(-1,-1,0,0,is), fna_he(-1,-1,0,0,is),
     &          fnaPSch(-1,-1,0,0,is), fna_mdf(-1,-1,0,0,is),
     &          flob, conb, fna_32(-1,-1,0,0,is), fna_52(-1,-1,0,0,is),
     &          fna_32nd(-1,-1,0,0,is), fna_52nd(-1,-1,0,0,is),
     &          kinrgy(-1,-1,is), fllim0fna(-1,-1,0,0,is),
     &          .false.)                                                 !djm Jan2017
           end do
*      ..compute electron drift velocities
        call b2tfed (nCv, nFc, nVx, switch, geo, mpg, dv%facdrift,
     &       dv%fac_ExB, dv%ne, pl%te, pl%po, dv%vedia, dv%veecrb)
*      ..compute other fluxes
           call b2xpen (nx, ny, ns,
     1      fna, fne, fch, fhe, fhi, kinrgy, rpt, te, ti, po, BoRiS,
     2      fhm, fhp, fhj, fht)
*      ..compute ue, electron parallel velocity
           ue=0.0_R8
           do iy=-1,ny
            do ix=-1,nx
             if(leftix(ix,iy).ne.-2 .and. rightix(ix,iy).ne.nx+1 .and.
     &         bottomiy(ix,iy).ne.-2 .and. topiy(ix,iy).ne.ny+1 .and.
     &         region(ix,iy,0).ne.0) then
              ue(ix,iy)=-(
     1           fch_p(ix,iy,0,0)/pbs(ix,iy,0)+                          !srv 10.06.08
     1           fch_p(rightix(ix,iy),rightiy(ix,iy),0,0)/               !srv 10.06.08
     2             pbs(rightix(ix,iy),rightiy(ix,iy),0))/qe/2.0_R8
             elseif(leftix(ix,iy).eq.-2 .and.
     &         bottomiy(ix,iy).ne.-2 .and. topiy(ix,iy).ne.ny+1 .and.
     &         region(ix,iy,0).ne.0) then
              ue(ix,iy)=-fch_p(rightix(ix,iy),rightiy(ix,iy),0,0)/       !srv 10.06.08
     2           pbs(rightix(ix,iy),rightiy(ix,iy),0)/qe
             elseif(rightix(ix,iy).eq.nx+1 .and.
     &         bottomiy(ix,iy).ne.-2 .and. topiy(ix,iy).ne.ny+1 .and.
     &         region(ix,iy,0).ne.0) then
              ue(ix,iy)=-fch_p(ix,iy,0,0)/pbs(ix,iy,0)/qe                !srv 10.06.08
             endif
            end do
           end do
           do is=0,ns-1
            ue(:,:)=ue(:,:)+rza(:,:,is)*ua(:,:,is)*na(:,:,is)
           enddo
           ue(:,:)=ue(:,:)+ue_ext(:,:)*ne_ext(:,:)
           ue(:,:)=ue(:,:)/ne(:,:)
*   ..compute ceqp
           call b2xpnm (nx, ny, ns, rz2, am, na, wrk0)
           call b2treq (nx, ny, vol, te, ti, ne, ni, wrk0, lnlam, ceqp) !srv 01.07.09 20.09.11
*      ..compute energy sources
           solving = .false.
           call b2npht (nx, ny, ns, nnreg(0), ncall, ismain, solving,  !srv 05.07.17
     &      solveee, solveei, solveen, solvepo, solvemo,
     &      rxf, BoRiS, fac_ExB, hx, hy, qz, qs, vol, bb,
     &      pbs, lnlam, floe, cone, floi, coni, ceqp, pccm, she, shi,  !srv 01.07.09
     &      fne, fch, cvsa, fchdia, uadia, wadia, vadia, vaecrb,       !srv 09.07.08 07.07.09
     &      na, ua, ue, te, ti, po, ne, ni, fhe, fhi,                  !srv 09.07.08
     &      reshe, reshi, wrk0, corte, corti, ier0)
#endif
         endif
*     ..end intermediate iteration
         istg(1) = istg(1)+1
         goto 2
       endif
*    ..end outer iteration
       istg(0) = istg(0)+1
       goto 1
      endif


      if (ibatch_av_all > 0) call batch_av_all(itim,ntim_batch,nCv,ns)

* ..return
99999 continue
      ncall = ncall+1
      call subend ()
      return

*-----------------------------------------------------------------------
*.end b2mndt

      end subroutine b2mndt
