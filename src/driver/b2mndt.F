*-----------------------------------------------------------------------
*.specification

      subroutine b2mndt (nout, nx, ny, ns, ismain, ismain0, 
     . nscx, iscx, itim, dtim, ierr)
      use b2mod_types
      use b2mod_geo
      use b2mod_anomalous_transport
      use b2mod_plasma
      use b2mod_rates
      use b2mod_residuals
      use b2mod_sources
      use b2mod_transport
      use b2mod_indirect
      use b2mod_constants
      use b2mod_diag
      use b2mod_work
      use b2mod_b2cmpa
      use b2mod_external
      use b2mod_neoclassical
      implicit none
*   ..input arguments (unchanged on exit)
      integer nout(0:5), nx, ny, ns, ismain, ismain0,  
     *  nscx, iscx(0:nscxmax-1), itim
      real (kind=R8) ::
     *  dtim
*   ..input/output arguments
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
     *  ramp_slow, ier0, no_solve, ix, iy
      integer, save :: style = 0, guard_flows = 2                        !srv 18.06.08
      integer, save :: anomalous_core_only = 0                           !srv 19.04.07
      integer, save :: ntim_step_out = 1                                 !srv 31.10.07
      integer, save :: use_b2srst = 1                                    !srv 13.07.12                                 
      real (kind=R8) ::
     *  rxf, BoRiS
      real (kind=R8) ::
     *  nstg_areshe, nstg_areshi, nstg_aresco
      save ncall, nstg, moqtlv, moitlv, rxf, itcnt,                      !xpb
     *  density_control, ramp_slow, BoRiS, no_solve,
     *  nstg_areshe, nstg_areshi, nstg_aresco
      logical domoqt, domoit, wrong_flow, ramp_call !, b2fstati_addexists  !xpb !srv 15.07.08
      logical test_residual, quit_residual
      save test_residual
      integer n2, istg(0:2), is, k
      real (kind=R8) ::
     *  acorpa(0:ns-1), acorua(0:ns-1), acorut,
     &  acorte, acorti, acorpo, nasum, na0sum,
     &  vsal(-1:nx,-1:ny,0:1,0:ns-1), pz(-1:nx,-1:ny),lnlam(-1:nx,-1:ny) !srv 11.10.02
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
      data no_solve/0/
      data nstg_areshe/0.0_R8/, nstg_areshi/0.0_R8/, nstg_aresco/0.0_R8/
*
*-----------------------------------------------------------------------
*.computation

* ..preliminaries
*   ..subprogram start-up calls
      call subini ('b2mndt')
*   ..set internal parameters on first call
      if (ncall.eq.0) then
       call ipgeti ('b2mndt_style', style)                            !srv 18.06.08
       call ipgeti ('b2news_guard_flows', guard_flows)                !srv 18.06.08
       call ipgeti ('b2tfnb_anomalous_core_only', anomalous_core_only) !srv 19.04.07
       call ipgeti ('b2mndt_nstg0', nstg(0))
       call ipgeti ('b2mndt_nstg1', nstg(1))
       call ipgeti ('b2mndt_nstg2', nstg(2))
       call xertst (1.le.nstg(0).and.1.le.nstg(1).and.1.le.nstg(2),
     &   'faulty internal parameter nstg(0:2)')
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
       call alloc_b2mod_diag(nx,ny,ns,nnreg(0),nstrat,natm)
       call ipgeti ('b2mndt_ntim_step_out', ntim_step_out)            !srv 31.01.07 {
       call xertst (1.le.ntim_step_out, 
     &   'faulty internal parameter ntim_step_out')                   !srv 31.01.07 }
       call ipgeti ('b2mndt_use_b2srst', use_b2srst)            !srv 13.07.12
      endif
*   ..set n2
      n2 = (nx+2)*(ny+2)
*   ..test nx, ny, ns, ismain
      call xertst (0.le.nx.and.0.le.ny, 'faulty argument nx, ny')
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
      call b2xvsg (n2, bb(-1,-1,3), 1, 'bb3', '.gt.')
*   ..test sign of vol, hx, hy
      call b2xvsg (n2, vol, 1, 'vol', '.gt.')
      call b2xvsg (n2, hx, 1, 'hx', '.gt.')
      call b2xvsg (n2, hy, 1, 'hy', '.gt.')
*   ..test range of qz
      call xertst (-1.0_R8.lt.smin(n2,qz(-1,-1,0),1).and.
     &             smax(n2,qz(-1,-1,0),1).lt.1.0_R8, 
     &             'faulty argument range qz(,,0)')
      call xertst (0.0_R8.lt.smin(n2,qz(-1,-1,1),1).and.
     &             smax(n2,qz(-1,-1,1),1).le.1.0_R8, 
     &             'faulty argument range qz(,,1)')
*   ..test old state
      call b2xvsg (n2*ns, na0, 1, 'na0', '.gt.')
      call b2xvsg (n2, te0, 1, 'te', '.gt.')
      call b2xvsg (n2, ti0, 1, 'ti', '.gt.')
*   ..test present state
      call b2xvps (nx, ny, ns, 
     &  ne, na, ua, te, ti, po, fna, fhe, fhi, fch)
      call b2xvsg (n2*ns, kinrgy, 1, 'kinrgy', '.ge.')
*   ..test facdrift and fac_ExB (not on first call)                      !xpb
      if (ncall.gt.0) then                                               !xpb
        call xertst (0.0_R8.le.minval(facdrift) .and.
     &       1.0_R8.ge.maxval(facdrift), 'faulty argument facdrift')
        call xertst (0.0_R8.le.minval(fac_ExB) .and.
     &       1.0_R8.ge.maxval(fac_ExB), 'faulty argument fac_ExB')
      else
c the main initialisation is in b2news --- unfortunately a few routines
c that need an initial value get called before the first call to b2news
c therefore initialise here
         call b2news_facdrift_fac_ExB_init(nx,ny,ramp_slow)
      endif
*   ..initialise error parameter
      ierr = 0
      quit_residual = .false.

* ..initialisation
*   ..initialise sources to 0
      call sfill (n2*2*ns, 0.0_R8, snasm, 1)                        !srv 11.09.09 {
      call sfill (n2*4*ns, 0.0_R8, smosm, 1)
      call sfill (n2*4, 0.0_R8, shesm, 1)
      call sfill (n2*4, 0.0_R8, shism, 1)
      call sfill (n2*2*ns, 0.0_R8, snadt, 1)
      call sfill (n2*4*ns, 0.0_R8, smodt, 1)
      call sfill (n2*4, 0.0_R8, shedt, 1)
      call sfill (n2*4, 0.0_R8, shidt, 1)
      call sfill (n2*4, 0.0_R8, schdt, 1)
      call sfill (n2*2, 0.0_R8, snedt, 1)
*   ..compute ni0
      call b2xpni (nx, ny, ns, na0, ni0)
*   ..compute nasum, na0sum
      nasum = b2sasum(n2*ns, na, 1)
      na0sum = b2sasum(n2*ns, na0, 1)
*   ..rescale densities
      if (density_control.ne.0) then
        ne = ne*na0sum/nasum
        na = na*na0sum/nasum
        fna = fna*na0sum/nasum
        kinrgy = kinrgy*na0sum/nasum
        write(*,*) 'Densities rescaled by :',na0sum/nasum
      endif
*   ..compute initial ni
      call b2xpni (nx, ny, ns, na, ni)
* ..main loop
*   ..outer iteration
      istg(0) = 0
    1 if (istg(0).lt.nstg(0) .and. .not.quit_residual) then
*    ..compute log-log electron rate coefficients
        call b2spel (nx, ny, ns, ev, te, ne, rlsa, rlra,
     &   rlqa, rlrd, rlbr, rlza, rlz2, rlpt, rlpi)
*    ..compute log-log charge-exchange rate coefficients
        do k = 0, nscx-1
         call b2spcx (nx, ny, ns, ev, am(iscx(k)), ti, ne,
     &     rlcx(-1,-1,0,0,k))
        enddo
ctmp.dpc
       if (nscx.ge.1) then
         write(*,*) 'DPC: b2mndt: rlcx'
         do k = 0, nscx-1
           do is=0, ns-1
             write(*,'(a,2i4,1p,2g14.7)') 'is0,is,rlcx_0,rlcx_1 ',
     1            k,is,samax(n2,rlcx(-1,-1,0,is,k),1),
     2            samax(n2,rlcx(-1,-1,1,is,k),1)
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
          newneo=.true.

          if(itcnt.eq.0) then
*   ..one-off calculation of transport coefficients
            call b2xpni (nx, ny, ns, na, ni)
            call b2xpne (nx, ny, ns, rza, na, ne_ext, ne)
            call b2xpne (nx, ny, ns, rz2, na, ne2_ext, ne2)
*   ..compute flux limit coefficients                                  !srv 11.10.02 03.06.03 {
            call b2trql (nx, ny, ns, bb, gs, vol, hx, hy, qz, qc, pbs, 
     &        na, ua, te, ti, po, ne, ni, ne2, chvemx, chvimx)         !srv 11.10.02 03.06.03 }
            call b2tral (nx, ny, ns, nscx, nscxmax, iscx, ismain,      !srv 23.11.10 
     1       bb, conn, vol, gs, hx, hy, hz, qz, qc, 
     1       pbs, crx, cry, bzb, lnlam,                                !srv 13.10.06 01.07.09
     1       fch, na, ua, te, ti, po, ne, ni, ne2, chvemx, chvimx,     !srv 03.06.03
!     1     fch_p, na, ua, te, ti, po, ne, ni, ne2, chvemx, chvimx,    !srv 03.06.03
     2       cdna, cdpa, cddi, cvla, cvsa, chce, chve, chci,           !srv 13.10.06
     3       chvi, csig, csigin, calf, cthe, cthi, 
     4       cdnahz, cdpahz, cvlahz, cvmahz, cvsahz, cvsa_cl, cvsahz_cl, !srv 09.01.01 15.12.05
     5       vsal, fllim0fhi, fllimvisc)
            newneo=.true.
            if(style.ne.0) then                                          !srv 18.06.08 {
!             INQUIRE (FILE = 'b2fstati_add', EXIST = b2fstati_addexists) !srv 15.07.08
!             IF (.NOT. b2fstati_addexists) THEN
!              WRITE (*,'(/2A/)') '  Cannot find file ', 'b2fstati_add'
!	        ! if there is no the file b2fstati_add then it is possible
!	        ! a jump of residuals when continuing of the run  
              if((maxval(facdrift).ne.0.0_R8 .or. 
     &            maxval(fac_ExB).ne.0.0_R8)) then

               call b2tvspa (nx, ny, ns, vol, hx, hy, qz, bb, ua, cvsa,   !srv 20.09.06  {
     &                        facdrift, fchvispar)
               call b2tiner (nx, ny, ns, vol,
     &                       gs, hx, hy, bb, na, ua, facdrift, fchinert)
               call b2tanml (nx, ny, ns, vol, hx, hy, csig_an, po,
     &                        fchanml)                                     !srv 20.09.06  }
               if(anomalous_core_only.ne.0) then                          !srv 19.04.07  {
                 do iy = 0, ny                                             !lk 19.09.07
                   do ix = -1, nx
                     if(mod(region(ix,bottomiy(ix,iy),0),4).ne.1) then     !lk 19.09.07
                       fchanml(ix,iy,1) = 0.0e0_R8
                     endif 
                   enddo
                 enddo
               endif                                                      !srv 19.04.07  }
              else
               fchvispar = 0.0_R8
               fchinert  = 0.0_R8
               fchanml   = 0.0_R8
              endif
              do is=0,ns-1                                               !srv 01.10.99 25.02.00
*   ..         compute partial pressure
               call b2xppb (nx, ny, rza(-1,-1,is), na(-1,-1,is), 
     &                      te, ti, pb)                                  !srv 01.10.99
*   ..         compute particle flux
               call b2tfnb (nx, ny, is, ismain, 
     &          guard_flows, facdrift, fac_ExB,  !srv 07.06.00 09.01.01 08.09.01 19.05.03 04.02.05
     &          am(is)*mp, vol, rza(-1,-1,is), pbs, pbshz, na(-1,-1,is), !srv 09.01.01 20.09.06
     &          pb, ua(-1,-1,is), uadia(-1,-1,0,is), vaecrb(-1,-1,0,is),   
     &          vadia(-1,-1,0,is), wadia(-1,-1,0,is), ne, te, ti, bb, 
     &          hx, hy, hz, qz, gs, qc, po,                              !srv 19.05.03 08.10.04 04.02.05
     &          cvla(-1,-1,0,is), cdna(-1,-1,0,is), cdpa(-1,-1,0,is),
     &          cddi(-1,-1,0,is), cvlahz(-1,-1,0,is), 
     &          cvmahz(-1,-1,0,is), cdnahz(-1,-1,0,is),                  !srv 09.01.01 15.12.05 13.10.06
     &          cdpahz(-1,-1,0,is), fchvispar, fchinert, fchanml,        !srv 20.09.06
     &          fna(-1,-1,0,is), fna_nodrift(-1,-1,0,is), 
     &          fna_fcor(-1,-1,0,is), fna_he(-1,-1,0,is), 
     &          fnaPSch(-1,-1,0,is), fna_mdf(-1,-1,0,is), 
     &          flob, conb, fna_32(-1,-1,0,is), fna_52(-1,-1,0,is), 
     &          fna_32nd(-1,-1,0,is), fna_52nd(-1,-1,0,is), 
     &          kinrgy(-1,-1,is), fllim0fna(-1,-1,0,is))
               enddo                                                     !srv 01.10.99 25.02.00 13.10.06
*    ..get current fch_p                                                !srv 11.10.02 {
              call b2tcpa (nx, ny, te, po, ne, csig, calf, fch_p)
              call b2xppz (nx, ny, ns, ne, na, te, ti, pz)

              call b2xpfe (nx, ny, ns, ns_ext, qe, 
     &         rza, za_ext, fna_he, fa_ext, fch_p, fne)
           
              call b2tfhe_ (nx, ny, ns, ismain, ismain0, facdrift,        
     1         vol, gs, hx, hy, qz, qc, bb, 
     1         te, po, ne, ne2, lnlam, pz, ti,                          !srv 01.07.09
     1         ni(-1,-1,1), na, fna, vadia, ua, fne, csig, csigin, chce, 
     &         chve, chci, calf, chvemx, cvsa, cdna, cdpa, cddi, fch,
     &         fhe, fhe_mdf, fhePSch, conc, avisper, bvisper, cvisper,  !srv 13.10.06
     &         floe, cone, fchdia, fchin, fch_p, fchvispar, fchvisper,
     &         fchvisq, fchinert, fchanml)                          !srv 24.03.07
              floe_noc = floe                                
              call b2xpfi (nx, ny, ns, fna_he, fni)                      !srv 28.10.02 {
csw 30apr2014 added fna_he, 17sep2014 added ua,cvsa
              call b2tfhi_ (nx, ny, ns, facdrift, vol, hx, hy, qz, bb,   
     &         ti, ni, na, fni, fna_he,ua,chci,cddi,chvi,chvimx,cvsa,
     &         vadia,                                               !srv 13.10.06
     &         fchvispar, fchvisq, fchinert, fchanml,
     &         fhi, fhi_mdf, fhiPSch, floi, coni)
              floi_noc = floi                                       !srv 24.10.02
!	       ELSE
!              open (70, file='b2fstati_add',form='unformatted')
!	          read(70,end=10) vaecrb, vadia, floe_noc, floi_noc, fna, 
!     &          fna_fcor, fna_nodrift, fna_he, fhe, fhi, fch,
!     &          fchdia, fchin, fch_p, fchvispar, fchvisper, fchvisq,
!     &          fchinert, fnaPSch, fna_mdf, fhePSch, fhiPSch,               !srv 13.10.06
!     &          fhe_mdf, fhi_mdf                                          !srv 29.10.09 {
!	          goto 20
!   10	          write(*,*) ' File b2fstati_add is old. '
!   	          write(*,*) ' fhe_mdf = fhe and fhi_mdf = fhi'
!                fhe_mdf = fhe
!	          fhi_mdf = fhi
!   20	        close(70)                                                   !srv 29.10.09 }
!	       ENDIF
            else ! if(style)
             vaecrb=0.0_R8; vadia=0.0_R8; floe_noc=0.0_R8;                !srv 19.11.09 {
             floi_noc=0.0_R8; fna_fcor=0.0_R8; fna_nodrift=0.0_R8;
             fna_he=0.0_R8; fchdia=0.0_R8; fchin=0.0_R8; fch_p=0.0_R8;
             fchvispar=0.0_R8; fchvisper=0.0_R8; fchvisq=0.0_R8;
             fchinert=0.0_R8; fnaPSch=0.0_R8; fhePSch=0.0_R8; 
             fhiPSch=0.0_R8; fna_mdf=0.0_R8; fhe_mdf=0.0_R8;
             fhi_mdf=0.0_R8;                                               !srv 19.11.09 }
            endif                                                        !srv 18.06.08 }
          endif
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
*     ..compute source linearisation
         call b2sral (nx, ny, ns, nscx, nscxmax, iscx, ismain, ismain0,
     &    dtim, BoRiS, facdrift, fac_ExB,
     &    vol, crx, hx, hy, hz, qz, qc, gs, pbs, pbshz, bb,                !srv 09.01.01
     &    na, ua, uadia, vedia, vadia, wadia, veecrb, vaecrb, ve, wedia,
     &    te, ti, po, ne, ni, kinrgy, floe_noc, floi_noc,                  !srv 28.10.02
     &    fna, fna_32, fna_52, fni_32, fni_52, fne_32, fne_52, 
     &    fna_mdf, fhe_mdf, fhi_mdf, fna_fcor, fna_nodrift, fna_he,        !srv 13.10.06 25.12.06 23.10.09
     &    fhe, fhi, fnaPSch, fhePSch, fhiPSch, fch, 
     &    fchdia, fchin, fch_p, fchvispar, fchvisper, fchvisq, fchinert,   !srv 21.10.02
     &    fchanml, fna_eir, fne_eir, fhe_eir, fhi_eir, 
     &    cdna, cvsa, cvla, chce, chve, chci, chvi, calf,   !xpb
     &    rlsa, rlra, rlqa, rlcx, rlrd, rlbr, rlza, rlz2, rlpt, rlpi,
     &    rza, rz2, rpt, rpi, sna, smo, smq, she, shi, sch, sne,
     &    wrong_flow, .true.)
          if(wrong_flow) write(*,*) 
     &         'b2mndt: wrong_flow returned from b2sral'
*     ..smooth out sources
csw 17oct2012 incorrect call to b2srsm !
c          call b2srsm (nx, ny, ns, dtim, vol, na, ua, te, ti, ne, ni,
c     &     snasm, smosm, shesm, shism, snesm, .true.)                       !srv 11.09.09 {
c*   ..add to previous contributions
c          call b2saxpy (n2*2*ns, 1.0_R8, snasm, 1, sna, 1)
c          call b2saxpy (n2*4*ns, 1.0_R8, smosm, 1, smo, 1)
c          call b2saxpy (n2*4, 1.0_R8, shesm, 1, she, 1)
c          call b2saxpy (n2*4, 1.0_R8, shism, 1, shi, 1)                   !srv 11.09.09 }
csw 17oct2012 this is the corrected call
          call b2srsm (nx, ny, ns, 
     &     dtim, vol, na, ua, te, ti, ne, ni,
     &     sna, smo, she, shi, sne, .true.)                       !srv 11.09.09 {

*     ..include contributions due to timestep
          call b2srdt (nx, ny, ns, guard_flows, dtim, BoRiS, hz,         !srv 09.01.01 01.07.08
     &     vol, na0, ua0, te0, ti0, ne0, ni0, kinrgy0,
     &     na, ua, te, ti, ne, ni, kinrgy,                          !srv 16.07.10
     &     snadt, smodt, shedt, shidt, schdt, snedt, .false.)                       !srv 11.09.09
*   ..add to previous contributions
          call b2saxpy (n2*2*ns, 1.0_R8, snadt, 1, sna, 1)
          call b2saxpy (n2*4*ns, 1.0_R8, smodt, 1, smo, 1)
          call b2saxpy (n2*4, 1.0_R8, snedt, 1, sne, 1)
          call b2saxpy (n2*4, 1.0_R8, shedt, 1, she, 1)
          call b2saxpy (n2*4, 1.0_R8, shidt, 1, shi, 1)
          call b2saxpy (n2*4, 1.0_R8, schdt, 1, sch, 1)                   !srv 11.09.09 }
*     ..stabilise the source coefficients
          if(use_b2srst.eq.1)
     &    call b2srst (nx, ny, ns, guard_flows, na, ua, te, ti, po, ne,  !srv 01.07.08
     &     ni, sna, smo, she, shi, sch, sne)
*       (smq is inherently stable.)
*     ..inner iteration
         istg(2) = 0
    3    if (istg(2).lt.nstg(2) .and. .not.quit_residual) then
*      ..basic iterative step
           if(style.eq.0) then                                          !srv 01.07.08
            ramp_call = ramp_slow.eq.0.or.(istg(0)+istg(1)+istg(2)).eq.0
            call b2news (nx, ny, ns, nscx, iscx, ismain, ismain0,
     &       natm, nmol, nion, nlim, nsts, itim,
     &       dtim, rxf, ier0, wrong_flow, ramp_call)                    !xpb
           else                                                         !srv 01.07.08 {
            call b2news_(nx, ny, ns, nscx, iscx, ismain, ismain0,
     &       natm, nmol, nion, nlim, nsts, itim,
     &       dtim, rxf,ier0)                                            !srv 16.07.10 23.11.10
	   endif                                                        !srv 01.07.08 )
           if(ier0.eq.10) then
             ierr=ier0
             goto 99999
           endif   
           call xertst (ier0.eq.0, 'error return from b2news')
*!!      (should develop other error handling)
           ierr = ier0
*      ..re-compute ne, ne2, ni
           call b2xpne (nx, ny, ns, rza, na, ne_ext, ne)
           call b2xpne (nx, ny, ns, rz2, na, ne2_ext, ne2)
           call b2xpni (nx, ny, ns, na, ni)
*      ..zero out residuals on special regions
*        (this in order to avoid their contribution to the norms)
           call b2mxzr (nx, ny, ns, resco, resmo, resmt,
     &      reshe, reshi, respo)
*      ..compute norms of the residuals
           call b2mxar (nx, ny, ns, nnreg(0), BoRiS, 
     &      vol, ne, na, te, ti, 
     &      resco, resmo, resmt, reshe, reshi, respo,
     &      rescoreg, resmoreg, reshereg, reshireg,
     &      aresco, aresmo, aresmt, areshe, areshi, arespo)
           quit_residual=test_residual .and.
     &      ((areshe*nstg_areshe).le.(nstg_areshe*nstg_areshe).or.
     &       nstg_areshe.eq.0.0_R8).and.
     &      ((areshi*nstg_areshi).le.(nstg_areshi*nstg_areshi).or.
     &       nstg_areshi.eq.0.0_R8).and.
     &      ((aresco(ismain)*nstg_aresco).le.(nstg_aresco*nstg_aresco)
     &       .or.nstg_aresco.eq.0.0_R8)
*      ..compute norms of the corrections
           call b2mxac (nx, ny, ns, 
     &      corpa, corua, corut, corte, corti, corpo,
     &      acorpa, acorua, acorut, acorte, acorti, acorpo)
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
            if (mod(ncall,ntim_step_out) .eq. 0)                         !srv 31.01.07
     &      call b2mwqt (nout(4), nx, ny, ns, itim, itcnt, vol,
     &       na, ua, te, ti, po, ne, ni, fna, fhe, fhi, fch,
     &       aresco, aresmo, aresmt, areshe, areshi, arespo,
     &       acorpa, acorua, acorut, acorte, acorti, acorpo)
           endif
*      ..maybe detailed monitor
           if (domoit) then
             call b2mwit (nout(5), nx, ny, ns, itim, itcnt, resco,
     &        resmo, resmt, reshe, reshi, respo,
     &        corpa, corua, corut, corte, corti,
     &        corpo)
           endif
*      ..write the tracing data
           if (ank_tracing.ne.0) then
             if (ncall.eq.0) call b2filei(ns)
             if (mod(itim,ank_tracing).eq.1.or.ank_tracing.eq.1) 
     >        call b2filew(istg(2),ns) 
           endif
*      ..update iteration counter
           itcnt = itcnt+1
*      ..end inner iteration
           istg(2) = istg(2)+1
           goto 3
         endif
*     ..end intermediate iteration
         istg(1) = istg(1)+1
         goto 2
       endif
*    ..end outer iteration
       istg(0) = istg(0)+1
       goto 1
      endif

* ..return
99999 continue
      ncall = ncall+1
      call subend ()
      return

*-----------------------------------------------------------------------
*.end b2mndt

      end subroutine b2mndt
