*-----------------------------------------------------------------------
*.specification

      subroutine b2tfch_ (nCv, nFc, nVx, ns, ismain, ismain0, geo, mpg, !srv 29.01.20 {
     &  pl, dv, co, rt)                                                 
      use b2mod_types
      use b2mod_boundary_namelist
      use b2mod_constants
      use b2mod_b2cmfs
      use b2mod_b2cmpa
      use b2us_geo
      use b2us_map
      use b2us_plasma
      implicit none
*   ..input arguments (unchanged on exit)
      integer nCv, nFc, nVx, ns, ismain, ismain0
      type (geometry), intent (in) :: geo
      type (mapping), intent (in) :: mpg
      type(B2Plasma), intent (in) :: pl 
      type(B2Rates), intent (in) :: rt
*   ..output arguments (unspecified on entry)
*   ..input/output arguments (unspecified on entry)
      type(B2Coeff), intent (inout) :: co
      type(B2Derivatives), intent (inout) :: dv
*   ..local variables
*   ..common blocks
*-----------------------------------------------------------------------
*.documentation
*
*  1. purpose
*
*     B2TFCH_ computes the currents and the coefficients that specify the
*     linearized response to changes in electric potential.
*
*  8. further comments
*
*-----------------------------------------------------------------------
*.declarations

*   ..local variables
      integer ncall, is, is0
      integer, save :: iout = 0, iout_b2wdat = 0,                        !srv 16.09.02
     &                 iout_b2npp7 = 0, iout_b2nppo = 0
      integer, save :: pot_eq = 1                                        !srv 26.07.06
      integer, save :: no_current= 0                                     !srv 27.02.13
csw 10oct2012
      integer, save :: b2_upwind = 0

      character(3) :: charns, charns0
      real (kind=R8), save :: drift_hyb = 1.0_R8
      real (kind=R8) :: csig_stoch(nFc), fchstoch(nFc)                   !lk 31.01.08 {
      real (kind=R8), save :: fch_stochastic = 0.0e0_R8                  !lk 31.01.08 }
      real (kind=R8) ::
     &    ehx(nFc), dia_cur, vis_par, vis_q, vis_per                     !srv 17.05.99 11.07.99
      real (kind=R8) ::
     *    fch_ion_neutral, fch_anomalous, fch_inert                      !srv 30.09.02 07.07.99
      real (kind=R8) :: facdriftm, fac_ExBm                              !srv 23.12.19
      real (kind=R8) ::
     *    wrk(nFc,0:1), dumm0(nFc,0:1), dumm1(nFc,0:1)

      save ncall, fch_ion_neutral, fch_anomalous, dia_cur, fch_inert,    !srv 07.07.99
     &     vis_par, vis_per, vis_q                                       !srv 17.05.99 11.07.99 22.11.99
*   ..procedures
      intrinsic max, abs
      external subini, subend, xertst
      external b2xvsg, b2xvff, b2xvfx, b2ttia                            !srv 11.07.99
*   ..initialization
      data ncall/0/, fch_ion_neutral/0.0e0_R8/, fch_inert/1.0e0_R8/,     !srv 07.07.99 22.03.02
     *   fch_anomalous/1.0e0_R8/, dia_cur/1.0e0_R8/, vis_par/0.0e0_R8/,  !srv 22.03.02
     *   vis_q/1.0e0_R8/, vis_per/0.0e0_R8/                              !srv 11.07.99 22.11.99 22.03.02

*-----------------------------------------------------------------------
*.computation

* ..preliminaries
*   ..subprogram start-up calls
      call subini ('b2tfch_')
*   ..set internal parameters on first call 
      if (ncall.eq.0) then 
       call ipgetr ('b2tfhe_neutral', fch_ion_neutral) 
       call ipgetr ('b2tfhe_anomalous', fch_anomalous) 
       if (fch_anomalous.eq.1.0_R8)
     >  call ipgetr ('b2tanml_anomalous', fch_anomalous) ! for backward compatibility
       call ipgetr ('b2tfhe_dia_cur', dia_cur) 
       call ipgetr ('b2tfhe_vis_par', vis_par) 
       call ipgetr ('b2tfhe_vis_per', vis_per)                           !srv 22.11.99
       call ipgetr ('b2tfhe_vis_q', vis_q)                               !srv 11.07.99
       call ipgetr ('b2tiner_inert', fch_inert)                          !srv 07.07.99
       call ipgetr ('b2tfhe_stochastic', fch_stochastic)                 !lk 31.01.08
       call xertst(fch_ion_neutral .ne. 0.0e0_R8 .or. 
     *             fch_anomalous .ne. 0.0e0_R8 .or.
     *             vis_q .ne. 0.0_R8 .or. fch_inert .ne. 0.0_R8 .or.
     *             vis_per .ne. 0.0_R8 .or.
     *             fch_stochastic .ne. 0.0_R8 .or.
     *             dia_cur .ne. 0.0e0_R8 .or. vis_par .ne. 0.0_R8, 
     *             'No radial current terms!') 
       call ipgeti ('b2tfch_iout', iout)                                 !srv 16.09.02
       call ipgeti ('b2news_poteq', pot_eq)                              !srv 26.07.06
       call ipgeti ('b2npp7_iout', iout_b2npp7)
       call ipgeti ('b2nppo_iout', iout_b2nppo)
       call ipgeti ('b2wdat_iout', iout_b2wdat)                          !iys 17.04.17
       call ipgeti ('b2tfhe_no_current',no_current)
       call xertst (pot_eq.eq.1.or.no_current.eq.1,
     &  'If b2news_poteq.ne.1, then b2tfhe_no_current must be 1!')
       call ipgeti ('b2tfhe_upwind', b2_upwind)
       call ipgetr ('b2tfhe_drift_hyb',drift_hyb)
       endif 
*   ..test nCv, nFc
      call xertst (0.le.nCv.and.0.le.nFc, 'faulty argument nCv, nFc')
*   ..extensive tests on first few calls
      if (ncall.lt.3) then
*    ..test sign of ne, te, csig, chce
       call b2xvsg (nCv, dv%ne, 1, 'ne', '.gt.')
       call b2xvsg (nCv, pl%te, 1, 'te', '.gt.')
       call b2xvsg (nFc, co%csig(:,0)*geo%fcQalf(:,0),1,'csig0','.ge.')
       call b2xvsg (nFc, co%csig(:,1)*geo%fcQalf(:,1),1,'csig1','.ge.')
       call b2xvsg (nFc, co%chce(:,0)*geo%fcQalf(:,0),1,'chce0','.ge.')
       call b2xvsg (nFc, co%chce(:,1)*geo%fcQalf(:,1),1,'chce1','.ge.')
      endif

      facdriftm = maxval(dv%facdrift)
      fac_ExBm  = maxval(dv%fac_ExB)                                     !srv 23.12.19
      dv%fchin     = 0.0_R8
      dv%fchdia    = 0.0_R8
      dv%fchvispar = 0.0_R8
      dv%fchvisper = 0.0_R8
      dv%fchvisq   = 0.0_R8
      dv%fchinert  = 0.0_R8
      dv%fch_P     = 0.0_R8
!WG_TODO      avisper = 0.0_R8
!WG_TODO      bvisper = 0.0_R8
!WG_TODO      cvisper = 0.0_R8
      fchstoch     = 0.0_R8
      csig_stoch   = 0.0_R8
*   ..calculate the contribution of ion-neutral current                  !srv 11.02.04 {
      call b2tinnt (nCv, nFc, nVx, ns, geo, mpg, rt, dv%facdrift,
     &       dv%fac_ExB, ismain, ismain0, pl%te, pl%po, pl%ti, pl%na,
     &       dv%fna, co%csigin, dv%fchin)                                !srv 11.05.09
*   ..compute x,y-components of divergent part of diamagnetic current
      call b2tdia(nCv, nFc, nVx, ns, geo, mpg, dv%pz, dv%fna, dv%uadia,
     &            dv%facdrift, dv%fchdia)
*   .. calculate x,y-components of viscosity parallel current
      call b2tvspa (nCv, nFc, nVx, ns, geo, mpg, pl%ua, co%vsaf_cl, 
     &  dv%fac_vis, dv%fchvispar)
#ifdef WG_TODO
*   ..calculate the contribution of viscosity perpendicular current      !srv 22.11.99
      call b2tvspr(nx, ny, ns, fac_vis, vol, hx, hy, qz, bb,             !srv 10.10.17
     &             po, ti, na, avisper, bvisper, cvisper, fchvisper)
*   ..calculate the parallel viscosity current connected with the heat flux
      call b2tvsq (nx, ny, ns, fac_vis, vol, hx, hy, qz, qs, bb, ne2,    !srv 10.10.17
     &             lnlam, te, ti, ni, chci, fchvisq)
#endif
*  .. calculate inertia current                                          !srv 07.07.99
      call b2tiner (nCv, nFc, nVx, ns, geo, mpg, pl%na, pl%ua,
     &      dv%facdrift, dv%fchinert)
*   ..calculate the x,y-components of anomalous current
      call b2tanml (nCv, nFc, ns, geo, mpg, co%csig_an,pl%po,dv%fchanml)
#ifdef WG_TODO
*   ..calculate the y-component of electron current associated with
*     stochastization of magnetic field
      call b2tstch (nx, ny, vol, hy, qz, bb, po, te, ne, csig_stoch,     !srv 07.07.14
     &              fchstoch)                                            !srv 25.01.02 {
#endif
      call b2tcpa (nCv, nFc, nVx, geo, mpg, pl%te, pl%po, dv%ne,
     *    co%csig, co%calf, ehx, dv%fch_p)                               !srv 13.01.17

*   ..compute currents
      dv%fch = dv%fch_p + dv%fchin + dv%fchdia + dv%fchvispar + 
     &  dv%fchvisq + dv%fchinert +dv%fchanml 
      dv%fch(:,1) = dv%fch(:,1) + dv%fchvisper(:,1) + fchstoch

*   ..compute conc
      dumm0 = 0.0e0_R8
      wrk   = co%csig
      if(fch_ion_neutral.ne.0.0e0_R8 .and.
     &     (facdriftm.gt.0.0_R8 .or. fac_ExBm.gt.0.0_R8)) then           !srv 23.12.19
        do is = 0, ns-1
          if (.not.is_neutral(is).and.nint(zn(is)).eq.1) then
            do is0 = 0, ns-1
              if (is_neutral(is0).and.nint(zn(is0)).eq.1) then
                wrk(:,0) = wrk(:,0) + co%csigin(:,0,is,is0)*
     &            (geo%fcBb(:,2)/geo%fcBb(:,3))**2
                wrk(:,1) = wrk(:,1) + co%csigin(:,1,is,is0)
              endif
            enddo
          endif
        enddo
      endif
      if(fch_stochastic.ne.0.0e0_R8 ) then                               !lk 31.01.08 {
        wrk(:,1) = wrk(:,1) + csig_stoch
      endif                                                              !lk 31.01.08 }
      call calccoef(nCv,nFc,nVx,0,geo,mpg,dumm0,wrk,dumm1,dv%conc)

      if(iout.ne.0) then                                                 !srv 18.06.08 {
       call my_out_us(70,nFc,1,dv%fch(1,0),'b2tfch__fchx')
       call my_out_us(70,nFc,1,dv%fch(1,1),'b2tfch__fchy')
       call my_out_us(70,nFc,1,dv%fchin(1,0),'b2tfch__fchinx')
       call my_out_us(70,nFc,1,dv%fchin(1,1),'b2tfch__fchiny')
       call my_out_us(70,nFc,1,dv%fchvispar(1,0),'b2tfch__fchvisparx')
       call my_out_us(70,nFc,1,dv%fchvispar(1,1),'b2tfch__fchvispary')
       call my_out_us(70,nFc,1,dv%fchvisq(1,0),'b2tfch__fchvisqx')
       call my_out_us(70,nFc,1,dv%fchvisq(1,1),'b2tfch__fchvisqy')
       call my_out_us(70,nFc,1,dv%fchvisper(1,0),'b2tfch__fchvisperx')
       call my_out_us(70,nFc,1,dv%fchvisper(1,1),'b2tfch__fchvispery')
       call my_out_us(70,nFc,1,dv%fchinert(1,0),'b2tfch__fchinertx')
       call my_out_us(70,nFc,1,dv%fchinert(1,1),'b2tfch__fchinerty')
       call my_out_us(70,nFc,1,dv%fchdia(1,0),'b2tfch__fchdiax')
       call my_out_us(70,nFc,1,dv%fchdia(1,1),'b2tfch__fchdiay')
       call my_out_us(70,nFc,1,dv%fchanml(1,0),'b2tfch__fchanmlx')
       call my_out_us(70,nFc,1,dv%fchanml(1,1),'b2tfch__fchanmly')
       call my_out_us(70,nFc,1,fchstoch,'b2tfch__fchstochy')             !lk 31.01.08
       call my_out_us(70,nFc,1,dv%conc(1,0),'b2tfch__concx')
       call my_out_us(70,nFc,1,dv%conc(1,1),'b2tfch__concy')
       call my_out_us(70,nFc,1,dv%fch_p(1,0),'b2tfch__fch_px')
       call my_out_us(70,nFc,1,dv%fch_p(1,1),'b2tfch__fch_py')
       call my_out_us(70,nFc,1,co%csig_an(1,0),'b2tfch__csig_anx')
       call my_out_us(70,nFc,1,co%csig_an(1,1),'b2tfch__csig_any')
       call my_out_us(70,nFc,1,co%csig(1,0),'b2tfch__csigx')
       call my_out_us(70,nFc,1,co%csig(1,1),'b2tfch__csigy')
       call my_out_us(70,nFc,1,co%calf(1,0),'b2tfch__calfx')
       call my_out_us(70,nFc,1,co%calf(1,1),'b2tfch__calfy')
       call my_out_us(70,nCv,0,pl%te,'b2tfch__te')
       call my_out_us(70,nCv,0,dv%ne,'b2tfch__ne')
       call my_out_us(70,nCv,0,pl%po,'b2tfch__po')
       do is = 0, ns-1
        write (charns,'(i3.3)') is
        call my_out_us(70,nCv,0,pl%na(1,is),'b2tfch__na'//charns)
        do is0=0,ns-1
         write (charns0,'(i3.3)') is0
         call my_out_us(70,nFc,1,co%csigin(1,0,is,is0),'b2tfch__csiginx'
     &        //charns//'_'//charns0)
         call my_out_us(70,nFc,1,co%csigin(1,1,is,is0),'b2tfch__csiginy'
     &        //charns//'_'//charns0)
        enddo
       enddo
      endif                                                              !srv 18.06.08 }
      
#ifdef WG_TODO
      if(iout_b2npp7.ne.0 .or. iout_b2nppo.ne.0
     &                          .or. iout_b2wdat.eq.4) then
!       call my_out(70,nx,ny,fch(-1,-1,0,0),'b2tfch__fchx')
!       call my_out(70,nx,ny,fch(-1,-1,1,1),'b2tfch__fchy')
       call my_out(70,nx,ny,fchin(-1,-1,0,0),'b2tfch__fchinx')
       call my_out(70,nx,ny,fchin(-1,-1,1,1),'b2tfch__fchiny')
       call my_out(70,nx,ny,fchvispar(-1,-1,0,0),'b2tfch__fchvisparx')
       call my_out(70,nx,ny,fchvispar(-1,-1,1,1),'b2tfch__fchvispary')
!       call my_out(70,nx,ny,fchvisq(-1,-1,0,0),'b2tfch__fchvisqx')
       call my_out(70,nx,ny,fchvisq(-1,-1,1,1),'b2tfch__fchvisqy')
!       call my_out(70,nx,ny,fchvisper(-1,-1,0,0),'b2tfch__fchvisperx')
       call my_out(70,nx,ny,fchvisper(-1,-1,1,1),'b2tfch__fchvispery')
       call my_out(70,nx,ny,fchinert(-1,-1,0,0),'b2tfch__fchinertx')
       call my_out(70,nx,ny,fchinert(-1,-1,1,1),'b2tfch__fchinerty')
       call my_out(70,nx,ny,fchdia(-1,-1,0,0),'b2tfch__fchdiax')
       call my_out(70,nx,ny,fchdia(-1,-1,1,1),'b2tfch__fchdiay')
       call my_out(70,nx,ny,fchanml(-1,-1,0,0),'b2tfch__fchanmlx')
       call my_out(70,nx,ny,fchanml(-1,-1,1,1),'b2tfch__fchanmly')
       call my_out(70,nx,ny,fch_p(-1,-1,0,0),'b2tfch__fch_px')
!      call my_out(70,nx,ny,fch_p(-1,-1,1,1),'b2tfch__fch_py')
       call my_out(70,nx,ny,csig_an(-1,-1,0,0),'b2tfch__csig_anx')
       call my_out(70,nx,ny,csig_an(-1,-1,1,1),'b2tfch__csig_any')
       call my_out(70,nx,ny,fchstoch,'b2tfch__fchstochy')                !lk 31.01.08
      endif                                                              !srv 11.09.09 }
#endif
* ..return
      ncall = ncall+1
      call subend ()
      return

*-----------------------------------------------------------------------
*.end b2tfch_

      end                                                                !srv 29.01.20 }
