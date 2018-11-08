*-----------------------------------------------------------------------
*.specification

      subroutine b2tfhi (nx, ny, ns, guard_flows, facdrift, BoRiS,
     1  vol, hx, hy, gs, qz, qc, bb, ti, ni, na,
     2  fna, fni_32, fni_52, kinrgy,
     3  chci, cddi, chvi, chvimx, cvsa, vadia, vaecrb,
     &  fchvispar, fchvisper, fchvisq, fchinert, fchanml,            !lk 30.05.11
     &  fhi, fhi_mdf, fhiPSch, fhit, floi, coni, fllimi,             !srv 13.10.06
     &  balance)                                                     !djm Jan2017
      use b2mod_types
      use b2mod_boundary_namelist
      use b2mod_indirect
      use b2mod_constants
      use b2mod_b2cmpa
      use b2mod_balance                                              !djm Jan2017
     & , only : fhi_32, fhi_52, fhi_cond, fhi_dia, fhi_ecrb,
     &          fhi_strange, fhi_pschused, fhi_kevis, balance_netcdf,
     &          fhi_inert, fhi_vispar, fhi_visper, fhi_visq, fhi_anml
      implicit none
*   ..input arguments (unchanged on exit)
      integer nx, ny, ns, guard_flows
      real (kind=R8) ::
     *  ti(-1:nx,-1:ny), ni(-1:nx,-1:ny,0:1),
     &  fna(-1:nx,-1:ny,0:1,0:ns-1), cddi(-1:nx,-1:ny,0:1,0:ns-1),    !srv 13.10.06
     1  fni_32(-1:nx,-1:ny,0:1), fni_52(-1:nx,-1:ny,0:1),
     &  fchvispar(-1:nx,-1:ny,0:1), fchvisper(-1:nx,-1:ny,0:1),       !lk 30.05.11
     &  fchinert(-1:nx,-1:ny,0:1), fchanml(-1:nx,-1:ny,0:1),          !lk 30.05.11
     &  fchvisq(-1:nx,-1:ny,0:1),
     &  vol(-1:nx,-1:ny), chci(-1:nx,-1:ny,0:1), chvi(-1:nx,-1:ny,0:1),
     &  chvimx(-1:nx,-1:ny), cvsa(-1:nx,-1:ny,0:1,0:nx-1),
     &  hx(-1:nx,-1:ny), hy(-1:nx,-1:ny),
     &  gs(-1:nx,-1:ny,0:2), qc(-1:nx,-1:ny), 
     1  qz(-1:nx,-1:ny,0:1), bb(-1:nx,-1:ny,0:3), 
     .  na(-1:nx,-1:ny,0:ns-1), 
     2  vadia(-1:nx,-1:ny,0:1,0:ns-1), vaecrb(-1:nx,-1:ny,0:1,0:ns-1),
     3  facdrift(-1:nx,-1:ny), kinrgy(-1:nx,-1:ny,0:ns-1), BoRiS
      logical balance                                                !djm Jan2017
*   ..output arguments (unspecified on entry)
      real (kind=R8) ::
     *  fhi(-1:nx,-1:ny,0:1), fhit(-1:nx,-1:ny,0:1), 
     &  fhi_mdf(-1:nx,-1:ny,0:1), fhiPSch(-1:nx,-1:ny,0:1),          !srv 13.10.06
     1  floi(-1:nx,-1:ny,0:1), coni(-1:nx,-1:ny,0:1,0:2),            !srv 02.01.07
     2  fllimi(-1:nx,-1:ny)
*   ..local variables
      real (kind=R8) :: upfac, hybfac, scur !djm Jan2017
*   ..common blocks
*-----------------------------------------------------------------------
*.documentation
*
*  1. purpose
*
*     B2TFHI computes the all atom heat flux between neighboring cells,
*     and the coefficients that specify its linearisation.
*
*
*  3. description (see also routine b2cdca)
*
*     This routine computes the atom heat flux, fhi, and the
*     coefficients floi and coni that linearise the response of the
*     atom heat flux to changes in the atom temperature.
*
*-----------------------------------------------------------------------
*.declarations

*   ..local variables
      integer ncall, n2, ix, iy, is, ib, conduction_only, no_hybr, 
     &        ixl, iyl, ixb, iyb
      integer, save :: iout = 0                                      !srv 16.09.02
      integer, save :: lim_flux = 0, mdf_fhi = 0                     !srv 04.02.03 13.10.06
      integer, save :: drift_style = 1
      integer, save :: conductive_limit = 1
      real (kind=R8), save :: PSch=1.0e0_R8                          !srv 17.01.07
      real (kind=R8), save :: xcur = 0.0_R8                          !xpb 27.09.17
      real (kind=R8), save :: ycur = 1.0_R8                          !srv 07.07.08
      real (kind=R8), save :: vis_per = 1.0_R8, vis_q = 1.0_R8       !srv 02.10.17
      real (kind=R8) ::
     *  flo, dif, floi0, flomx, coni0, t0, t1, dti, dti0,
     1  cum_flo_vadia, cum_flo_vaecrb, prefactor, cutlo,
     &   floi0_mdf, coni0_mdf, t0_mdf,dti0_mdf, t1_mdf               !srv 13.10.06
      real (kind=R8) :: w0(-1:nx,-1:ny,0:ns-1),
     1  w1(-1:nx,-1:ny,0:ns-1), w2(-1:nx,-1:ny,0:ns-1)
      parameter (cutlo=1.0e-60_R8)
      save ncall, conduction_only, no_hybr
*   ..procedures
      intrinsic min, max, dim, abs, sqrt, sign
      real (kind=R8) ::
     *  hybr, upwind, hy1
      external subini, subend, xertst, ipgeti
      external b2xvsg, b2xvff
      hybr(flo,dif) = max(dif,abs(flo/2),sqrt(dif**2+flo**2/4)) !!!
      upwind(flo,t0,t1) = max(flo,0.0_R8)*t0+min(flo,0.0_R8)*t1
      hy1(ix,iy) = hy(ix,iy)*qz(ix,iy,1)
*   ..initialisation
      data ncall/0/, conduction_only/0/, no_hybr/0/

*-----------------------------------------------------------------------
*.computation

* ..preliminaries
*   ..subprogram start-up calls
      call subini ('b2tfhi')
*   ..set n2
      n2 = (nx+2)*(ny+2)
*   ..test nx, ny
      call xertst (0.le.nx.and.0.le.ny, 'faulty argument nx, ny')
*   ..test facdrift
      call xertst (0.0_R8.le.minval(facdrift) .and.
     &       1.0_R8.ge.maxval(facdrift), 'faulty argument facdrift')
*   ..get input parameters
      if (ncall.eq.0) then
        call ipgeti ('b2tfhe_conduction_only', conduction_only)
        call ipgeti ('b2tfnb_drift_style', drift_style)
        call ipgeti ('b2tfhi_no_hybr', no_hybr)
        call ipgeti ('b2tfhi_mdf', mdf_fhi)                          !srv 13.10.06
        call ipgetr ('b2tfhi_PSch', PSch)                            !srv 17.01.07
        call ipgetr ('b2tfnb_xcur', xcur)                            !xpb 27.09.17
        call ipgetr ('b2tfnb_ycur', ycur)                            !srv 17.07.08
        call ipgetr ('b2tfnb_vis_per', vis_per)                      !srv 02.10.17
        call ipgetr ('b2tfnb_vis_q', vis_q)                          !srv 02.10.17
        call ipgeti ('b2tfhi_iout', iout)                            !srv 16.09.02
        call ipgeti ('b2tfhi_lim_flux', lim_flux)
        call ipgeti ('b2trcl_conductive_limit', conductive_limit) 
      endif
*   ..extensive tests on first few calls
      if (ncall.lt.3) then
*    ..test edge values of fni
       call b2xvff (nx, ny, fni_32, 'fni_32')
       call b2xvff (nx, ny, fni_52, 'fni_52')
*    ..test sign of na, ni, ti, chci
       call b2xvsg (n2*ns, na, 1, 'na', '.gt.')
       call b2xvsg (n2*2, ni, 1, 'ni', '.gt.')
       call b2xvsg (n2, ti, 1, 'ti', '.gt.')
       call b2xvsg (n2*2, chci, 1, 'chci', '.ge.')
*    ..test edge values of chci, chvi
       call b2xvff (nx, ny, chci, 'chci')
       call b2xvff (nx, ny, chvi, 'chvi')
      endif

      call b2tipsch (nx, ny, ns, facdrift, vol, hx, hy, qz,
     1  na, ti, cddi, coni, fhiPSch)
      
      fhi_mdf = 0.0e0_R8                                              !srv 13.10.06

* ..compute fluxes and linearisation
      do is = 0, ns-1
        w0(:,:,is) = kinrgy(:,:,is)/(am(is)*mp)
      end do
*   ..compute fhix, floix, conix
      do iy = -1, ny
       do ix = -1, nx
        if(leftix(ix,iy).eq.-2) then
         floi(ix,iy,0) = 0.0_R8
         coni(ix,iy,0,0) = 0.0_R8                                    !srv 02.01.07
         fhi_mdf(ix,iy,0) = 0.0_R8
         fhi(ix,iy,0) = 0.0_R8
         fhit(ix,iy,0) = 0.0_R8
         fllimi(ix,iy) = 1.0_R8
         w1(ix,iy,:) = 0.0_R8
         if (balance.and.(balance_netcdf.ne.0)) then
           fhi_32(ix,iy,0) = 0.0_R8                               !djm Jan2017
           fhi_52(ix,iy,0) = 0.0_R8
           fhi_cond(ix,iy,0) = 0.0_R8
           fhi_dia(ix,iy,0) = 0.0_R8
           fhi_ecrb(ix,iy,0) = 0.0_R8
           fhi_strange(ix,iy,0) = 0.0_R8
           fhi_pschused(ix,iy,0) = 0.0_R8
           fhi_inert(ix,iy,0) = 0.0_R8
           fhi_vispar(ix,iy,0) = 0.0_R8
           fhi_visper(ix,iy,0) = 0.0_R8
           fhi_visq(ix,iy,0) = 0.0_R8
           fhi_anml(ix,iy,0) = 0.0_R8
         endif
        else
         dti = ti(ix,iy)-ti(leftix(ix,iy),leftiy(ix,iy))
         w1(ix,iy,:) = 0.0_R8
cxpb Changed the solution to (3/2+BoRiS) for benchmarking !!!!!!!!!!!!!!!!!!!
         if (conduction_only.eq.1) then
           floi0 = 0.0_R8
           cum_flo_vadia=0.0_R8
           cum_flo_vaecrb=0.0_R8
         else
           if (BoRiS.eq.1.0_R8) then   ! Add the kinetic and viscous terms
             do is = 0, ns-1
!wdk use of cvsa in coni0 ensures consistency w.r.t. flux limits and discretization
               floi0 = am(is)*mp*fna(ix,iy,0,is)
               coni0 = 4.0e0_R8/3.0e0_R8*cvsa(ix,iy,0,is)
               t1 = dim ( coni0,
     .          abs(floi0) *
     .          max(0.0_R8,
     .           min(w0(ix,iy,is),w0(leftix(ix,iy),leftiy(ix,iy),is)))/
     .          max(cutlo,
     .              w0(ix,iy,is)+w0(leftix(ix,iy),leftiy(ix,iy),is)) )
               w1(ix,iy,is) =
     .          upwind(floi0,
     .                 w0(leftix(ix,iy),leftiy(ix,iy),is),w0(ix,iy,is))
     .           - t1*(w0(ix,iy,is)-w0(leftix(ix,iy),leftiy(ix,iy),is))
             enddo
           endif
           scur =
     &      1.5e0_R8/qe*(fchinert(ix,iy,0)+fchvispar(ix,iy,0)
     &                  +fchanml(ix,iy,0)+fchvisper(ix,iy,0)
     &                  +fchvisq(ix,iy,0))*xcur
           floi0 = (1.5_R8+BoRiS)*fni_32(ix,iy,0)+
     1      2.5_R8*fni_52(ix,iy,0)+
     &      chvi(ix,iy,0)*(ni(leftix(ix,iy),leftiy(ix,iy),0)+
     1                     ni(ix,iy,0))/2.0_R8+scur
           cum_flo_vadia=0.0_R8
           cum_flo_vaecrb=0.0_R8
           do is=0,ns-1
             prefactor=
     1        (vol(leftix(ix,iy),leftiy(ix,iy))+vol(ix,iy))     !srv 17.05.99
     &        /(hx(leftix(ix,iy),leftiy(ix,iy))+hx(ix,iy))      !srv 17.05.99
     &        *(na(leftix(ix,iy),leftiy(ix,iy),is)+na(ix,iy,is))/2.0_R8
             if (drift_style.eq.0) then
               cum_flo_vadia=cum_flo_vadia + 0.5_R8*prefactor*
     1           (vadia(ix,iy,0,is) + 
     2            vadia(leftix(ix,iy),leftiy(ix,iy),0,is))
               cum_flo_vaecrb=cum_flo_vaecrb + 0.5_R8*prefactor*
     1           (vaecrb(ix,iy,0,is) + 
     2            vaecrb(leftix(ix,iy),leftiy(ix,iy),0,is))
             else
               cum_flo_vadia=cum_flo_vadia+
     1           prefactor*vadia(ix,iy,0,is)
               cum_flo_vaecrb=cum_flo_vaecrb+
     1           prefactor*vaecrb(ix,iy,0,is)
             endif
           enddo
         endif
         floi0_mdf = floi0+1.5_R8*cum_flo_vaecrb                !srv 18.08.08
         floi0 = floi0+2.5_R8*cum_flo_vadia+1.5_R8*cum_flo_vaecrb
         flomx = chvimx(ix,iy)*(ni(leftix(ix,iy),leftiy(ix,iy),0)
     &                         +ni(ix,iy,0))/2.0_R8
         if (no_hybr.eq.0) then
           coni0 = hybr(floi0,chci(ix,iy,0))-abs(floi0/2.0_R8)
           coni0_mdf = hybr(floi0_mdf,chci(ix,iy,0))-
     &                  abs(floi0_mdf/2.0_R8)                   !srv 18.08.08
         else
           coni0 = chci(ix,iy,0)
           coni0_mdf = chci(ix,iy,0)
         endif
         if (flomx.gt.0.0_R8 .and.
     &       lim_flux.eq.1 .and. conductive_limit.eq.0) then
           t0 = coni0*abs(dti)/
     &      (flomx*max(ti(leftix(ix,iy),leftiy(ix,iy)),ti(ix,iy)))
           t0_mdf = coni0_mdf*abs(dti)/                             !srv 18.08.08
     &      (flomx*max(ti(leftix(ix,iy),leftiy(ix,iy)),ti(ix,iy)))
         else
           t0 = 0.0_R8
           t0_mdf = 0.0_R8                                      !srv 18.08.08
         endif
         dti0 = dti/(1.0_R8+t0)
         t1 = t0/(1.0_R8+t0)
         dti0_mdf = dti/(1.0_R8+t0_mdf)             !srv 18.08.08 {
         t1_mdf = t0_mdf/(1.0_R8+t0_mdf)
         if(mdf_fhi.eq.0) then
          fllimi(ix,iy) = 1.0_R8+t0
         else
          fllimi(ix,iy) = 1.0_R8+t0_mdf
         endif                                      !srv 18.08.08 }
         if (no_hybr.eq.0) then
           fhi(ix,iy,0) = floi0*(ti(leftix(ix,iy),leftiy(ix,iy))
     &                          +ti(ix,iy))/2.0_R8-
     &                  abs(floi0/2.0_R8)*dti-coni0*dti0
           fhi_mdf(ix,iy,0) =
     .                    floi0_mdf*(ti(leftix(ix,iy),leftiy(ix,iy))
     &                          +ti(ix,iy))/2.0_R8-
     &                  abs(floi0_mdf/2.0_R8)*dti-coni0_mdf*dti0_mdf
     &                       - fhiPSch(ix,iy,0)                      !srv 13.10.06 }
           if (balance.and.
     &         (balance_netcdf.ne.0).and.
     &         (mdf_fhi.eq.0)) then !djm Jan2017
             !djm Jan2017 share the hybrid part according to relative contributions
             if (floi0.ge.0.0_R8) then
               hybfac=1.0_R8
             else
               hybfac=-1.0_R8
             endif             
             fhi_32(ix,iy,0) = (1.5_R8+BoRiS)*fni_32(ix,iy,0)*
     &         (ti(leftix(ix,iy),leftiy(ix,iy))+ti(ix,iy))/2.0_R8-
     &         hybfac*(1.5_R8+BoRiS)*fni_32(ix,iy,0)/2.0_R8*dti
             fhi_52(ix,iy,0) = 2.5_R8*fni_52(ix,iy,0)*
     &         (ti(leftix(ix,iy),leftiy(ix,iy))+ti(ix,iy))/2.0_R8-
     &         hybfac*2.5_R8*fni_52(ix,iy,0)/2.0_R8*dti
             fhi_strange(ix,iy,0) = chvi(ix,iy,0)*
     &         (ni(leftix(ix,iy),leftiy(ix,iy),0)+ni(ix,iy,0))/2.0_R8*
     &         (ti(leftix(ix,iy),leftiy(ix,iy))+ti(ix,iy))/2.0_R8-
     &         hybfac*chvi(ix,iy,0)*
     &         (ni(leftix(ix,iy),leftiy(ix,iy),0)+ni(ix,iy,0))/2.0_R8/
     &         2.0_R8*dti
             fhi_ecrb(ix,iy,0) = 1.5_R8*cum_flo_vaecrb*
     &         (ti(leftix(ix,iy),leftiy(ix,iy))+ti(ix,iy))/2.0_R8-
     &         hybfac*(1.5_R8*cum_flo_vaecrb)/2.0_R8*dti
             fhi_dia(ix,iy,0) = 2.5_R8*cum_flo_vadia*
     &         (ti(leftix(ix,iy),leftiy(ix,iy))+ti(ix,iy))/2.0_R8-
     &         hybfac*2.5_R8*cum_flo_vadia/2.0_R8*dti
             fhi_cond(ix,iy,0) = -coni0*dti0
             fhi_pschused(ix,iy,0) = 0.0_R8
             fhi_inert(ix,iy,0) = 1.5e0_R8/qe*fchinert(ix,iy,0)*xcur*
     &         (ti(leftix(ix,iy),leftiy(ix,iy))+ti(ix,iy))/2.0_R8-
     &         hybfac*1.5e0_R8/qe*fchinert(ix,iy,0)*xcur/2.0_R8*dti
             fhi_vispar(ix,iy,0) = 1.5e0_R8/qe*fchvispar(ix,iy,0)*xcur*
     &         (ti(leftix(ix,iy),leftiy(ix,iy))+ti(ix,iy))/2.0_R8-
     &         hybfac*1.5e0_R8/qe*fchvispar(ix,iy,0)*xcur/2.0_R8*dti
             fhi_visper(ix,iy,0) = 1.5e0_R8/qe*fchvisper(ix,iy,0)*xcur*
     &         (ti(leftix(ix,iy),leftiy(ix,iy))+ti(ix,iy))/2.0_R8-
     &         hybfac*1.5e0_R8/qe*fchvisper(ix,iy,0)*xcur/2.0_R8*dti
             fhi_visq(ix,iy,0) = 1.5e0_R8/qe*fchvisq(ix,iy,0)*xcur*
     &         (ti(leftix(ix,iy),leftiy(ix,iy))+ti(ix,iy))/2.0_R8-
     &         hybfac*1.5e0_R8/qe*fchvisq(ix,iy,0)*xcur/2.0_R8*dti
             fhi_anml(ix,iy,0) = 1.5e0_R8/qe*fchanml(ix,iy,0)*xcur*
     &         (ti(leftix(ix,iy),leftiy(ix,iy))+ti(ix,iy))/2.0_R8-
     &         hybfac*1.5e0_R8/qe*fchanml(ix,iy,0)*xcur/2.0_R8*dti
           endif
           if (balance.and.
     &         (balance_netcdf.ne.0).and.
     &         (mdf_fhi.ne.0)) then !djm Jan2017
             if (floi0_mdf.ge.0.0_R8) then
               hybfac=1.0_R8
             else
               hybfac=-1.0_R8
             endif             
             fhi_32(ix,iy,0) = (1.5_R8+BoRiS)*fni_32(ix,iy,0)*
     &         (ti(leftix(ix,iy),leftiy(ix,iy))+ti(ix,iy))/2.0_R8-
     &         hybfac*(1.5_R8+BoRiS)*fni_32(ix,iy,0)/2.0_R8*dti
             fhi_52(ix,iy,0) = 2.5_R8*fni_52(ix,iy,0)*
     &         (ti(leftix(ix,iy),leftiy(ix,iy))+ti(ix,iy))/2.0_R8-
     &         hybfac*2.5_R8*fni_52(ix,iy,0)/2.0_R8*dti
             fhi_strange(ix,iy,0) = chvi(ix,iy,0)*
     &         (ni(leftix(ix,iy),leftiy(ix,iy),0)+ni(ix,iy,0))/2.0_R8*
     &         (ti(leftix(ix,iy),leftiy(ix,iy))+ti(ix,iy))/2.0_R8-
     &         hybfac*chvi(ix,iy,0)*
     &         (ni(leftix(ix,iy),leftiy(ix,iy),0)+ni(ix,iy,0))/2.0_R8/
     &         2.0_R8*dti
             fhi_ecrb(ix,iy,0) = 1.5_R8*cum_flo_vaecrb*
     &         (ti(leftix(ix,iy),leftiy(ix,iy))+ti(ix,iy))/2.0_R8-
     &         hybfac*(1.5_R8*cum_flo_vaecrb)/2.0_R8*dti
             fhi_dia(ix,iy,0) = 0.0_R8
             fhi_cond(ix,iy,0) = -coni0_mdf*dti0_mdf
             fhi_pschused(ix,iy,0) = -fhiPSch(ix,iy,0)
             fhi_inert(ix,iy,0) = 1.5e0_R8/qe*fchinert(ix,iy,0)*xcur*
     &         (ti(leftix(ix,iy),leftiy(ix,iy))+ti(ix,iy))/2.0_R8-
     &         hybfac*1.5e0_R8/qe*fchinert(ix,iy,0)*xcur/2.0_R8*dti
             fhi_vispar(ix,iy,0) = 1.5e0_R8/qe*fchvispar(ix,iy,0)*xcur*
     &         (ti(leftix(ix,iy),leftiy(ix,iy))+ti(ix,iy))/2.0_R8-
     &         hybfac*1.5e0_R8/qe*fchvispar(ix,iy,0)*xcur/2.0_R8*dti
             fhi_visper(ix,iy,0) = 1.5e0_R8/qe*fchvisper(ix,iy,0)*xcur*
     &         (ti(leftix(ix,iy),leftiy(ix,iy))+ti(ix,iy))/2.0_R8-
     &         hybfac*1.5e0_R8/qe*fchvisper(ix,iy,0)*xcur/2.0_R8*dti
             fhi_visq(ix,iy,0) = 1.5e0_R8/qe*fchvisq(ix,iy,0)*xcur*
     &         (ti(leftix(ix,iy),leftiy(ix,iy))+ti(ix,iy))/2.0_R8-
     &         hybfac*1.5e0_R8/qe*fchvisq(ix,iy,0)*xcur/2.0_R8*dti
             fhi_anml(ix,iy,0) = 1.5e0_R8/qe*fchanml(ix,iy,0)*xcur*
     &         (ti(leftix(ix,iy),leftiy(ix,iy))+ti(ix,iy))/2.0_R8-
     &         hybfac*1.5e0_R8/qe*fchanml(ix,iy,0)*xcur/2.0_R8*dti
           endif
         else
           t1 = dim ( coni0/fllimi(ix,iy),
     .         abs(floi0) *
     .         max(0.0_R8, 
     .             min(ti(ix,iy),ti(leftix(ix,iy),leftiy(ix,iy))))/
     .         max(cutlo, ti(ix,iy)+ti(leftix(ix,iy),leftiy(ix,iy))) )
           fhi(ix,iy,0) =
     .         upwind(floi0, ti(leftix(ix,iy),leftiy(ix,iy)), ti(ix,iy))
     .         - t1 * (ti(ix,iy)-ti(leftix(ix,iy),leftiy(ix,iy)))
           fhi_mdf(ix,iy,0) =
     .         upwind(floi0_mdf,
     .                 ti(leftix(ix,iy),leftiy(ix,iy)), ti(ix,iy))
     .         - t1 * (ti(ix,iy)-ti(leftix(ix,iy),leftiy(ix,iy)))
           if (balance.and.
     &         (balance_netcdf.ne.0).and.
     &         (mdf_fhi.eq.0)) then !djm Jan2017
             if(floi0.ge.0.0_R8) then
               upfac=1.0_R8
             else
               upfac=0.0_R8
             endif
             fhi_32(ix,iy,0) = (1.5_R8+BoRiS)*fni_32(ix,iy,0)*
     &                         (upfac*ti(leftix(ix,iy),leftiy(ix,iy))+
     &                          (1.0_R8-upfac)*ti(ix,iy))
             fhi_52(ix,iy,0) = 2.5_R8*fni_52(ix,iy,0)*
     &                         (upfac*ti(leftix(ix,iy),leftiy(ix,iy))+
     &                          (1.0_R8-upfac)*ti(ix,iy))
             fhi_strange(ix,iy,0) = chvi(ix,iy,0)*
     &         (ni(leftix(ix,iy),leftiy(ix,iy),0)+ni(ix,iy,0))/2.0_R8*
     &         (upfac*ti(leftix(ix,iy),leftiy(ix,iy))+
     &          (1.0_R8-upfac)*ti(ix,iy))
             fhi_ecrb(ix,iy,0) = 1.5_R8*cum_flo_vaecrb*
     &                           (upfac*ti(leftix(ix,iy),leftiy(ix,iy))+
     &                            (1.0_R8-upfac)*ti(ix,iy))
             fhi_dia(ix,iy,0) = 2.5_R8*cum_flo_vadia*
     &                          (upfac*ti(leftix(ix,iy),leftiy(ix,iy))+
     &                           (1.0_R8-upfac)*ti(ix,iy))
             fhi_cond(ix,iy,0) = -coni0*dti0
             fhi_pschused(ix,iy,0) = 0.0_R8
             fhi_inert(ix,iy,0) = 1.5e0_R8/qe*fchinert(ix,iy,0)*xcur*
     &                         (upfac*ti(leftix(ix,iy),leftiy(ix,iy))+
     &                          (1.0_R8-upfac)*ti(ix,iy))
             fhi_vispar(ix,iy,0) = 1.5e0_R8/qe*fchvispar(ix,iy,0)*xcur*
     &                         (upfac*ti(leftix(ix,iy),leftiy(ix,iy))+
     &                          (1.0_R8-upfac)*ti(ix,iy))
             fhi_visper(ix,iy,0) = 1.5e0_R8/qe*fchvisper(ix,iy,0)*xcur*
     &                         (upfac*ti(leftix(ix,iy),leftiy(ix,iy))+
     &                          (1.0_R8-upfac)*ti(ix,iy))
             fhi_visq(ix,iy,0) = 1.5e0_R8/qe*fchvisq(ix,iy,0)*xcur*
     &                         (upfac*ti(leftix(ix,iy),leftiy(ix,iy))+
     &                          (1.0_R8-upfac)*ti(ix,iy))
             fhi_anml(ix,iy,0) = 1.5e0_R8/qe*fchanml(ix,iy,0)*xcur*
     &                         (upfac*ti(leftix(ix,iy),leftiy(ix,iy))+
     &                          (1.0_R8-upfac)*ti(ix,iy))
           endif
           if (balance.and.
     &         (balance_netcdf.ne.0).and.
     &         (mdf_fhi.ne.0)) then !djm Jan2017
             if(floi0_mdf.ge.0.0_R8) then
               upfac=1.0_R8
             else
               upfac=0.0_R8
             endif
             fhi_32(ix,iy,0) = (1.5_R8+BoRiS)*fni_32(ix,iy,0)*
     &                         (upfac*ti(leftix(ix,iy),leftiy(ix,iy))+
     &                          (1.0_R8-upfac)*ti(ix,iy))
             fhi_52(ix,iy,0) = 2.5_R8*fni_52(ix,iy,0)*
     &                         (upfac*ti(leftix(ix,iy),leftiy(ix,iy))+
     &                          (1.0_R8-upfac)*ti(ix,iy))
             fhi_strange(ix,iy,0) = chvi(ix,iy,0)*
     &         (ni(leftix(ix,iy),leftiy(ix,iy),0)+ni(ix,iy,0))/2.0_R8*
     &         (upfac*ti(leftix(ix,iy),leftiy(ix,iy))+
     &          (1.0_R8-upfac)*ti(ix,iy))
             fhi_ecrb(ix,iy,0) = 1.5_R8*cum_flo_vaecrb*
     &                           (upfac*ti(leftix(ix,iy),leftiy(ix,iy))+
     &                            (1.0_R8-upfac)*ti(ix,iy))
             fhi_dia(ix,iy,0) = 0.0_R8
             fhi_cond(ix,iy,0) = -coni0_mdf*dti0_mdf
             fhi_pschused(ix,iy,0) = -fhiPSch(ix,iy,0)
             fhi_inert(ix,iy,0) = 1.5e0_R8/qe*fchinert(ix,iy,0)*xcur*
     &                         (upfac*ti(leftix(ix,iy),leftiy(ix,iy))+
     &                          (1.0_R8-upfac)*ti(ix,iy))
             fhi_vispar(ix,iy,0) = 1.5e0_R8/qe*fchvispar(ix,iy,0)*xcur*
     &                         (upfac*ti(leftix(ix,iy),leftiy(ix,iy))+
     &                          (1.0_R8-upfac)*ti(ix,iy))
             fhi_visper(ix,iy,0) = 1.5e0_R8/qe*fchvisper(ix,iy,0)*xcur*
     &                         (upfac*ti(leftix(ix,iy),leftiy(ix,iy))+
     &                          (1.0_R8-upfac)*ti(ix,iy))
             fhi_visq(ix,iy,0) = 1.5e0_R8/qe*fchvisq(ix,iy,0)*xcur*
     &                         (upfac*ti(leftix(ix,iy),leftiy(ix,iy))+
     &                          (1.0_R8-upfac)*ti(ix,iy))
             fhi_anml(ix,iy,0) = 1.5e0_R8/qe*fchanml(ix,iy,0)*xcur*
     &                         (upfac*ti(leftix(ix,iy),leftiy(ix,iy))+
     &                          (1.0_R8-upfac)*ti(ix,iy))
           endif
         endif
         do is = 0, ns-1
           fhi(ix,iy,0) = fhi(ix,iy,0) + w1(ix,iy,is)
           fhi_mdf(ix,iy,0) = fhi_mdf(ix,iy,0) + w1(ix,iy,is)
         end do
         if (balance_netcdf.ne.0) then
           fhi_kevis(ix,iy,0) = 0.0_R8
           do is = 0, ns-1
             fhi_kevis(ix,iy,0) = fhi_kevis(ix,iy,0)+w1(ix,iy,is)
           enddo
         endif
         fhit(ix,iy,0) = fhi(ix,iy,0) +
     1    ((1.0_R8-BoRiS)*fni_32(ix,iy,0)+cum_flo_vaecrb)*
     &    (ti(leftix(ix,iy),leftiy(ix,iy))+ti(ix,iy))/2.0_R8
         if(mdf_fhi.eq.0) then
          floi(ix,iy,0) = floi0-
     &    sign(flomx,ti(ix,iy)-ti(leftix(ix,iy),leftiy(ix,iy)))*
     &    min(0.0_R8,t1)**2
         else
          floi(ix,iy,0) = floi0_mdf-
     &    sign(flomx,ti(ix,iy)-ti(leftix(ix,iy),leftiy(ix,iy)))*
     &    min(0.0_R8,t1_mdf)**2
         endif
         if (no_hybr.eq.0) then
          if(mdf_fhi.eq.0) then
           coni(ix,iy,0,0) = abs(floi0/2.0_R8)+coni0/(1.0_R8+t0)     !srv 02.01.07
          else
           coni(ix,iy,0,0) = abs(floi0_mdf/2.0_R8)+
     &                       coni0_mdf/(1.0_R8+t0_mdf)     !srv 02.01.07
          endif
         else
           coni(ix,iy,0,0) = coni0/(1.0_R8+t0)                       !srv 02.01.07
         endif
        endif
       enddo
      enddo

*   ..compute fhiy, floiy, coniy
      do iy = -1, ny
       do ix = -1, nx
        if(bottomiy(ix,iy).eq.-2) then
         floi(ix,iy,1) = 0.0_R8
         coni(ix,iy,1,0) = 0.0_R8                                    !srv 02.01.07
         fhi_mdf(ix,iy,1) = 0.0e0_R8                                 !srv 13.10.06
         fhi(ix,iy,1) = 0.0_R8
         fhit(ix,iy,1) = 0.0_R8
         w2(ix,iy,:) = 0.0_R8
         if (balance.and.(balance_netcdf.ne.0)) then
           fhi_32(ix,iy,1) = 0.0_R8                               !djm Jan2017
           fhi_52(ix,iy,1) = 0.0_R8
           fhi_cond(ix,iy,1) = 0.0_R8
           fhi_dia(ix,iy,1) = 0.0_R8
           fhi_ecrb(ix,iy,1) = 0.0_R8
           fhi_strange(ix,iy,1) = 0.0_R8
           fhi_pschused(ix,iy,1) = 0.0_R8
           fhi_inert(ix,iy,1) = 0.0_R8
           fhi_vispar(ix,iy,1) = 0.0_R8
           fhi_visper(ix,iy,1) = 0.0_R8
           fhi_visq(ix,iy,1) = 0.0_R8
           fhi_anml(ix,iy,1) = 0.0_R8
         endif
        else
         w2(ix,iy,:) = 0.0_R8
         if (conduction_only.eq.1) then
           floi0 = 0.0_R8
           floi0_mdf = 0.0_R8
           cum_flo_vadia=0.0_R8
           cum_flo_vaecrb=0.0_R8
         else
           if (BoRiS.eq.1.0_R8) then   ! Add the kinetic and viscous terms
             do is = 0, ns-1
               floi0 = am(is)*mp*fna(ix,iy,1,is)
               coni0 = cvsa(ix,iy,1,is)
               t1 = dim ( coni0,
     .                   abs(floi0) *
     .          max(0.0_R8, 
     .             min(w0(ix,iy,is),
     .                 w0(bottomix(ix,iy),bottomiy(ix,iy),is)))/
     .          max(cutlo, 
     .           w0(ix,iy,is)+w0(bottomix(ix,iy),bottomiy(ix,iy),is)))
               w2(ix,iy,is) =
     .          upwind(floi0,
     .           w0(bottomix(ix,iy),bottomiy(ix,iy),is), w0(ix,iy,is))
     .           - t1 * (w0(ix,iy,is)-
     .                   w0(bottomix(ix,iy),bottomiy(ix,iy),is))
             enddo
           endif
           scur =
     &      1.5e0_R8/qe*(fchinert(ix,iy,1)+fchvispar(ix,iy,1)
     &                  +fchanml(ix,iy,1)+vis_per*fchvisper(ix,iy,1)      !srv 02.10.17
     &                  +vis_q*fchvisq(ix,iy,1))*ycur                     !srv 02.10.17
           floi0 = (1.5_R8+BoRiS)*fni_32(ix,iy,1)+
     1      2.5_R8*fni_52(ix,iy,1)+
     2      chvi(ix,iy,1)*(ni(bottomix(ix,iy),bottomiy(ix,iy),0)
     &                    +ni(ix,iy,0))/2.0_R8+scur
           cum_flo_vadia=0.0_R8
           cum_flo_vaecrb=0.0_R8
           do is=0,ns-1
             prefactor=
     &         (vol(bottomix(ix,iy),bottomiy(ix,iy))+vol(ix,iy))   !srv 17.05.99
     &        /(hy1(bottomix(ix,iy),bottomiy(ix,iy))+hy1(ix,iy))   !srv 17.05.99
     &        *(na(bottomix(ix,iy),bottomiy(ix,iy),is)             !srv 17.05.99
     &         +na(ix,iy,is))/2.0_R8                               !srv 17.05.99
             if (drift_style.eq.0) then
               cum_flo_vadia=cum_flo_vadia+
     1          0.5_R8*prefactor*(vadia(ix,iy,1,is) + 
     2          vadia(bottomix(ix,iy),bottomiy(ix,iy),1,is))
               cum_flo_vaecrb=cum_flo_vaecrb+
     1          0.5_R8*prefactor*(vaecrb(ix,iy,1,is) + 
     2          vaecrb(bottomix(ix,iy),bottomiy(ix,iy),1,is))
             else
               cum_flo_vadia=cum_flo_vadia+
     1          prefactor*vadia(ix,iy,1,is)
               cum_flo_vaecrb=cum_flo_vaecrb+
     1          prefactor*vaecrb(ix,iy,1,is)
             endif
           enddo
         endif
         if(mdf_fhi.eq.0) then
          floi(ix,iy,1) = floi0+
     1     2.5_R8*cum_flo_vadia+1.5_R8*cum_flo_vaecrb
         else
          floi(ix,iy,1) = floi0+1.5_R8*cum_flo_vaecrb
         endif
         if (no_hybr.eq.0) then
           coni(ix,iy,1,0) = hybr(floi(ix,iy,1),chci(ix,iy,1))    !srv 02.01.07
           fhi(ix,iy,1) = floi(ix,iy,1)*
     1      (ti(bottomix(ix,iy),bottomiy(ix,iy))
     &                                +ti(ix,iy))/2.0_R8-
     &       coni(ix,iy,1,0)*                                     !srv 02.01.07
     1      (ti(ix,iy)-ti(bottomix(ix,iy),bottomiy(ix,iy)))
           fhi_mdf(ix,iy,1) = floi(ix,iy,1)*
     1      (ti(bottomix(ix,iy),bottomiy(ix,iy))
     &                                +ti(ix,iy))/2.0_R8-
     &       coni(ix,iy,1,0)*                                     !srv 02.01.07
     1      (ti(ix,iy)-ti(bottomix(ix,iy),bottomiy(ix,iy)))
     &      -fhiPSch(ix,iy,1)                                     !srv 13.10.06 }
           if (balance.and.(balance_netcdf.ne.0)) then !djm Jan2017
             fhi_32(ix,iy,1) = (1.5_R8+BoRiS)*fni_32(ix,iy,1)*
     &      (ti(bottomix(ix,iy),bottomiy(ix,iy))+ti(ix,iy))/2.0_R8
             fhi_52(ix,iy,1) = 2.5_R8*fni_52(ix,iy,1)*
     &      (ti(bottomix(ix,iy),bottomiy(ix,iy))+ti(ix,iy))/2.0_R8
             fhi_strange(ix,iy,1) = chvi(ix,iy,1)*
     &      (ni(bottomix(ix,iy),bottomiy(ix,iy),0)+ni(ix,iy,0))/2.0_R8*
     &      (ti(bottomix(ix,iy),bottomiy(ix,iy))+ti(ix,iy))/2.0_R8
             if (mdf_fhi.eq.0) then
               fhi_dia(ix,iy,1) = 2.5_R8*cum_flo_vadia*
     &         (ti(bottomix(ix,iy),bottomiy(ix,iy))+ti(ix,iy))/2.0_R8
               fhi_pschused(ix,iy,1) = 0.0_R8
             else
               fhi_dia(ix,iy,1) = 0.0_R8
               fhi_pschused(ix,iy,1) = -fhiPSch(ix,iy,1)
             endif
             fhi_ecrb(ix,iy,1) = 1.5_R8*cum_flo_vaecrb*
     &        (ti(bottomix(ix,iy),bottomiy(ix,iy))+ti(ix,iy))/2.0_R8
             fhi_cond(ix,iy,1) = -coni(ix,iy,1,0)*
     &        (ti(ix,iy)-ti(bottomix(ix,iy),bottomiy(ix,iy)))
             fhi_inert(ix,iy,1) = 1.5e0_R8/qe*fchinert(ix,iy,1)*ycur*
     &        (ti(bottomix(ix,iy),bottomiy(ix,iy))+ti(ix,iy))/2.0_R8
             fhi_vispar(ix,iy,1) = 1.5e0_R8/qe*fchvispar(ix,iy,1)*ycur*
     &        (ti(bottomix(ix,iy),bottomiy(ix,iy))+ti(ix,iy))/2.0_R8
             fhi_visper(ix,iy,1) =
     &         1.5e0_R8/qe*vis_per*fchvisper(ix,iy,1)*ycur*         !srv 02.10.17
     &        (ti(bottomix(ix,iy),bottomiy(ix,iy))+ti(ix,iy))/2.0_R8
             fhi_visq(ix,iy,1) =
     &         1.5e0_R8/qe*vis_q*fchvisq(ix,iy,1)*ycur*             !srv 02.10.17
     &        (ti(bottomix(ix,iy),bottomiy(ix,iy))+ti(ix,iy))/2.0_R8
             fhi_anml(ix,iy,1) = 1.5e0_R8/qe*fchanml(ix,iy,1)*ycur*
     &        (ti(bottomix(ix,iy),bottomiy(ix,iy))+ti(ix,iy))/2.0_R8
           endif
         else
           coni(ix,iy,1,0) = chci(ix,iy,1)                        !srv 02.01.07
           t1 = dim ( coni(ix,iy,1,0),                            !srv 02.01.07
     .         abs(floi(ix,iy,1)) *
     .         max(0.0_R8,
     .             min(ti(ix,iy),ti(bottomix(ix,iy),bottomiy(ix,iy))))/
     .         max(cutlo, 
     .                 ti(ix,iy)+ti(bottomix(ix,iy),bottomiy(ix,iy))))
           fhi(ix,iy,1) =
     .         upwind(floi(ix,iy,1),
     .                ti(bottomix(ix,iy),bottomiy(ix,iy)), ti(ix,iy))
     .         - t1 * (ti(ix,iy)-ti(bottomix(ix,iy),bottomiy(ix,iy)))
            fhi_mdf(ix,iy,1) =
     .         upwind(floi(ix,iy,1),
     .                ti(bottomix(ix,iy),bottomiy(ix,iy)), ti(ix,iy))
     .         - t1 * (ti(ix,iy)-ti(bottomix(ix,iy),bottomiy(ix,iy)))
        endif
         do is = 0, ns-1
           fhi(ix,iy,1) = fhi(ix,iy,1) + w2(ix,iy,is)
           fhi_mdf(ix,iy,1) = fhi_mdf(ix,iy,1) + w2(ix,iy,is)
         end do
         fhit(ix,iy,1) = fhi(ix,iy,1)+
     1    ((1.0_R8-BoRiS)*fni_32(ix,iy,1)+cum_flo_vaecrb)*
     &    (ti(bottomix(ix,iy),bottomiy(ix,iy))+ti(ix,iy))/2.0_R8
        endif
       enddo
      enddo

      if(guard_flows.eq.0) then                                    !xpb
cxpb turn off flows between boundary cells
        do iy = -1, ny
          do ix = -1, nx
            if(leftix(ix,iy).eq.-2 .or. rightix(ix,iy).eq.nx+1) then
              fhi(ix,iy,1)=0.0_R8
              if (balance.and.(balance_netcdf.ne.0)) then !djm Jan2017
                fhi_32(ix,iy,1) = 0.0_R8
                fhi_52(ix,iy,1) = 0.0_R8
                fhi_cond(ix,iy,1) = 0.0_R8
                fhi_dia(ix,iy,1) = 0.0_R8
                fhi_ecrb(ix,iy,1) = 0.0_R8
                fhi_strange(ix,iy,1) = 0.0_R8
                fhi_pschused(ix,iy,1) = 0.0_R8
                fhi_inert(ix,iy,1) = 0.0_R8
                fhi_vispar(ix,iy,1) = 0.0_R8
                fhi_visper(ix,iy,1) = 0.0_R8
                fhi_visq(ix,iy,1) = 0.0_R8
                fhi_anml(ix,iy,1) = 0.0_R8
             endif
              floi(ix,iy,1)=0.0_R8
              fhit(ix,iy,1)=0.0_R8
              coni(ix,iy,1,:)=0.0_R8
            endif
          enddo
        enddo

cxpb turn off poloidal flows in core boundary
        do ix = -1, nx
          if(on_closed_surface(ix,-1)) then
            fhi(ix,-1,0) = 0.0_R8
            if (balance.and.(balance_netcdf.ne.0)) then !djm Jan2017
              fhi_32(ix,-1,0) = 0.0_R8
              fhi_52(ix,-1,0) = 0.0_R8
              fhi_cond(ix,-1,0) = 0.0_R8
              fhi_dia(ix,-1,0) = 0.0_R8
              fhi_ecrb(ix,-1,0) = 0.0_R8
              fhi_strange(ix,-1,0) = 0.0_R8
              fhi_pschused(ix,-1,0) = 0.0_R8
              fhi_inert(ix,-1,0) = 0.0_R8
              fhi_vispar(ix,-1,0) = 0.0_R8
              fhi_visper(ix,-1,0) = 0.0_R8
              fhi_visq(ix,-1,0) = 0.0_R8
              fhi_anml(ix,-1,0) = 0.0_R8
            endif
            floi(ix,-1,0) = 0.0_R8
            fhit(ix,-1,0) = 0.0_R8
            coni(ix,-1,0,:) = 0.0_R8
            fllimi(ix,-1) = 1.0_R8
          endif
        enddo
      endif                                                        !xpb

*  INTERIOR BOUNDARY AT DEAD REGIONS
      DO IB=1,NNISO
        DO IX=NXISO1(IB), NXISO2(IB)
          DO IY=NYISO1(IB), NYISO2(IB)
            fhi(ix,iy,0)=0.0_R8
            if (balance.and.(balance_netcdf.ne.0)) then !djm Jan2017
              fhi_32(ix,iy,0) = 0.0_R8
              fhi_52(ix,iy,0) = 0.0_R8
              fhi_cond(ix,iy,0) = 0.0_R8
              fhi_dia(ix,iy,0) = 0.0_R8
              fhi_ecrb(ix,iy,0) = 0.0_R8
              fhi_strange(ix,iy,0) = 0.0_R8
              fhi_pschused(ix,iy,0) = 0.0_R8
              fhi_inert(ix,iy,0) = 0.0_R8
              fhi_vispar(ix,iy,0) = 0.0_R8
              fhi_visper(ix,iy,0) = 0.0_R8
              fhi_visq(ix,iy,0) = 0.0_R8
              fhi_anml(ix,iy,0) = 0.0_R8
              fhi_32(ix,iy,1) = 0.0_R8
              fhi_52(ix,iy,1) = 0.0_R8
              fhi_cond(ix,iy,1) = 0.0_R8
              fhi_dia(ix,iy,1) = 0.0_R8
              fhi_ecrb(ix,iy,1) = 0.0_R8
              fhi_strange(ix,iy,1) = 0.0_R8
              fhi_pschused(ix,iy,1) = 0.0_R8
              fhi_inert(ix,iy,1) = 0.0_R8
              fhi_vispar(ix,iy,1) = 0.0_R8
              fhi_visper(ix,iy,1) = 0.0_R8
              fhi_visq(ix,iy,1) = 0.0_R8
              fhi_anml(ix,iy,1) = 0.0_R8
            endif
            floi(ix,iy,0)=0.0_R8
            fhit(ix,iy,0)=0.0_R8
            coni(ix,iy,0,:)=0.0_R8
            fllimi(ix,iy)=1.0_R8
            fhi(ix,iy,1)=0.0_R8
            floi(ix,iy,1)=0.0_R8
            fhit(ix,iy,1)=0.0_R8
            coni(ix,iy,1,:)=0.0_R8
          end do
        end do
      end do

      if(iout.ne.0) then                                             !srv 18.06.08 {
        call my_out(70,nx,ny,floi(-1,-1,0),'b2tfhi_floix')
        call my_out(70,nx,ny,floi(-1,-1,1),'b2tfhi_floiy')
        call my_out(70,nx,ny,coni(-1,-1,0,0),'b2tfhi_conix0')
        call my_out(70,nx,ny,coni(-1,-1,0,1),'b2tfhi_conix1')
        call my_out(70,nx,ny,coni(-1,-1,0,2),'b2tfhi_conix2')
        call my_out(70,nx,ny,coni(-1,-1,1,0),'b2tfhi_coniy0')
        call my_out(70,nx,ny,coni(-1,-1,1,1),'b2tfhi_coniy1')
        call my_out(70,nx,ny,coni(-1,-1,1,2),'b2tfhi_coniy2')
        call my_out(70,nx,ny,fhi(-1,-1,0),'b2tfhi_fhix')
        call my_out(70,nx,ny,fhi(-1,-1,1),'b2tfhi_fhiy')
        call my_out(70,nx,ny,fhi_mdf(-1,-1,0),'b2tfhi_fhi_mdfx')
        call my_out(70,nx,ny,fhi_mdf(-1,-1,1),'b2tfhi_fhi_mdfy')
        call my_out(70,nx,ny,fhiPSch(-1,-1,0),'b2tfhi_fhiPSchx')
        call my_out(70,nx,ny,fhiPSch(-1,-1,1),'b2tfhi_fhiPSchy')
      endif                                                          !srv 18.06.08 }

* ..return
      ncall = ncall+1
      call subend ()
      return

*-----------------------------------------------------------------------
*.scribbles

*!!   Reconsider factor 2.5_R8 versus 1.5_R8 in transport and linearisation.
cxpb  Factor changed to 1.5_R8.
cxpb  This version of the routine uses a (1.5+BoRiS) factor for the
cxpb  poloidal convective heat flux instead of the usual 3/2 or 5/2 
cxpb  so as to solve the same equation as in BoRiS for benchmark purposes. 
cxpb  The default is BoRiS = 0.0.
cxpb  For comparison, set BoRiS = 1.0.
cxpb  Routine rewritten to explicitly separate 3/2 and 5/2 terms.
cxpb  The contributions from fni_32 and fni_52 are NO DRIFTS
cxpb  contributions, i.e. diamagnetic and ExB terms are treated
cxpb  separately.
cxpb  The BoRiS treatment must be completed with addition of the 
cxpb  viscous and kinetic terms, as per the SOLPS4.0 treatment.

*-----------------------------------------------------------------------
*.end b2tfhi

      end subroutine b2tfhi
