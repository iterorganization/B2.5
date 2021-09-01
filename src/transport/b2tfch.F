*-----------------------------------------------------------------------
*.specification

      subroutine b2tfch (nx, ny, ns, ismain, ismain0,
     1  guard_flows, facdrift, fac_ExB, fac_vis, BoRiS,
     1  vol, gs, hx, hy, qz, qc, bb, te, po, ne, ne2, lnlam, pz,         !xpb
     1  ti, ni, na, fna, uadia, vedia, wadia, vaecrb,
     2  veecrb, ua, fne_32, fne_52, csig, csigin, csig_stoch,            !xpb !srv 24.08.21
     &  chce, chve, chci,
     &  calf, chvemx, cvsa_cl, fch, conc, fchstoch, fchdia,              !srv 17.12.13 16.10.17
     &  fch_p, fchvispar, fchvisq, fchin, fchinert, fch_32, fch_52,      !srv 25.11.08
     1  avisper, bvisper, cvisper, fchvisper, csig_cl, calf_cl)
      use b2mod_types
      use b2mod_math
      use b2mod_rates
      use b2mod_boundary_namelist
      use b2mod_indirect
      use b2mod_constants
      use b2mod_anomalous_transport                                      !xpb
      use b2mod_geo_corner                                               !xpb
      use b2mod_b2cmfs
      use b2mod_b2cmpa
      use b2mod_subsys
      implicit none
*   ..input arguments (unchanged on exit)
      integer nx, ny, ns, ismain, ismain0, guard_flows                   !xpb
      real (kind=R8) ::
     * te(-1:nx,-1:ny), po(-1:nx,-1:ny), ni(-1:nx,-1:ny),
     & ne(-1:nx,-1:ny), ne2(-1:nx,-1:ny), lnlam(-1:nx,-1:ny),
     & pz(-1:nx,-1:ny), vol(-1:nx,-1:ny),
     & fna(-1:nx,-1:ny,0:1,0:ns-1),
     & gs(-1:nx,-1:ny,0:2), qz(-1:nx,-1:ny,0:1), qc(-1:nx,-1:ny),        !xpb
     & uadia(-1:nx,-1:ny,0:1,0:ns-1), 
     & ua(-1:nx,-1:ny,0:ns-1),
     & wadia(-1:nx,-1:ny,0:1,0:ns-1),
     & vaecrb(-1:nx,-1:ny,0:1,0:ns-1),                                   !xpb
     & fne_32(-1:nx,-1:ny,0:1), fne_52(-1:nx,-1:ny,0:1),
     1 csig(-1:nx,-1:ny,0:1), chci(-1:nx,-1:ny,0:1),
     & csigin(-1:nx,-1:ny,0:1,0:ns-1,0:ns-1),                             !xpb
     & chce(-1:nx,-1:ny,0:1), chve(-1:nx,-1:ny,0:1),
     & calf(-1:nx,-1:ny,0:1), chvemx(-1:nx,-1:ny),
     & cvsa_cl(-1:nx,-1:ny,0:1,0:ns-1), na(-1:nx,-1:ny,0:ns-1),          !srv 16.10.17
     1 hx(-1:nx,-1:ny), hy(-1:nx,-1:ny),
     1 bb(-1:nx,-1:ny,0:3), ti(-1:nx,-1:ny), 
     2 facdrift(-1:nx,-1:ny), BoRiS,
     3 fac_ExB(-1:nx,-1:ny), fac_vis(-1:nx,-1:ny),
     4 csig_cl(-1:nx,-1:ny,0:1), calf_cl(-1:nx,-1:ny,0:1)
      real (kind=R8) :: csig_stoch(-1:nx,-1:ny,0:1)                     !lk 31.01.08 } 24.08.21
*   ..output arguments (unspecified on entry)
      real (kind=R8) ::
     *  fch(-1:nx,-1:ny,0:1), 
     &  conc(-1:nx,-1:ny,0:1),
     *  vedia(-1:nx,-1:ny,0:1), veecrb(-1:nx,-1:ny,0:1),
     2  fchdia(-1:nx,-1:ny,0:1), 
     3  fch_p(-1:nx,-1:ny,0:1), fchvispar(-1:nx,-1:ny,0:1),              !xpb srv 10.06.08
     4  fch_32(-1:nx,-1:ny,0:1), fch_52(-1:nx,-1:ny,0:1),
     5  avisper(-1:nx,-1:ny,0:1), bvisper(-1:nx,-1:ny,0:1),
     6  cvisper(-1:nx,-1:ny,0:1), fchvisper(-1:nx,-1:ny,0:1),
     7  fchvisq(-1:nx,-1:ny,0:1), fchin(-1:nx,-1:ny,0:1),                !srv 25.11.08
     &  fchstoch(-1:nx,-1:ny,0:1), fchinert(-1:nx,-1:ny,0:1)             !srv 24.08.21
*   ..common blocks
*-----------------------------------------------------------------------
*.documentation
*
*  1. purpose
*
*     B2TFCH computes the electric current between neighboring cells,
*     and the coefficients that specify the linearised response to
*     changes in electric potential.
*
*-----------------------------------------------------------------------
*.declarations

*   ..local variables
      integer ncall, n2, ix, iy, ib, is, is0, 
     *   no_current, conduction_only, no_hybr, pot_eq
      integer, save :: iout = 0                                         !srv 18.06.08
      integer, save :: lim_flux = 0
      integer, save :: conductive_limit = 1
      integer, save :: style = 2                                        !srv 13.01.17
      real (kind=R8), save :: fch_stochastic = 0.0e0_R8                 !lk 31.01.08 {
      real (kind=R8), save :: fch_inert = 1.0e0_R8                      !lk 31.01.08 {
      real (kind=R8), save :: dia_cur = 1.0e0_R8                        !srv 17.05.99
      real (kind=R8), save :: vis_q = 1.0e0_R8                          !srv 11.07.99
      real (kind=R8), save :: prl_cur = 1.0_R8, phm2 = 1.0_R8           !srv 15.10.18
      real (kind=R8) ::
     *   ehx, ehy, flomx,
     &   tex, ehx_po, ehy_po,
     1   floe0, floe0_nd, cone0, cone0_nd, t0, t0_nd, tex0, tex0_nd     !srv 15.10.18
      real (kind=R8) :: facdriftm, fac_ExBm                             !srv 23.12.19
      real (kind=R8) ::                                                 !xpb
     * fch_ion_neutral, fch_anomalous, vis_par, vis_per, sig_stoch,     !xpb
     * ehx_(-1:nx,-1:ny), ehy_(-1:nx,-1:ny), zeff(-1:nx,-1:ny)          !srv 13.01.17
      save ncall, fch_ion_neutral, fch_anomalous, vis_par, pot_eq,
     *   no_current, conduction_only, no_hybr, vis_per, sig_stoch
*   ..procedures
      external xertst, ipgetr, ipgeti
      external b2xvsg, b2xvff, b2xvfx, intcor
*   ..initialisation
      data ncall/0/, fch_ion_neutral/0.0_R8/, fch_anomalous/1.0_R8/,    !xpb
     * vis_par/1.0_R8/, no_current/0/, conduction_only/0/,
     1 vis_per/1.0_R8/, no_hybr/0/, pot_eq/1/, sig_stoch/1.0_R8/

*-----------------------------------------------------------------------
*.computation

* ..preliminaries
*   ..subprogram start-up calls
      call subini ('b2tfch')
      if (ncall.eq.0) then                                              !xpb
       call ipgetr ('b2tfhe_neutral', fch_ion_neutral)                  !xpb
       call ipgetr ('b2tfhe_anomalous', fch_anomalous)                  !xpb
       call ipgetr ('b2tfhe_stochastic', fch_stochastic)                !lk 31.01.08
       call ipgetr ('b2tstch_sig', sig_stoch)
       call ipgeti ('b2tfhe_lim_flux', lim_flux)
       call ipgeti ('b2trcl_conductive_limit', conductive_limit) 
       if (fch_anomalous.eq.1.0_R8)
     >  call ipgetr('b2tanml_anomalous', fch_anomalous) ! for backward compatibility
       call ipgetr ('b2tfhe_vis_q', vis_q)                              !srv 11.07.99
       call ipgetr ('b2tfhe_vis_par', vis_par)                          !xpb
       call ipgetr ('b2tfhe_vis_per', vis_per)                          !xpb
       call ipgetr ('b2tfhe_dia_cur', dia_cur)                          !srv
       call ipgetr ('b2tiner_inert', fch_inert)                         !srv 07.07.99
       call xertst (fch_anomalous.ne.0.0_R8
     &         .or.(fch_stochastic.ne.0.0_R8.and.sig_stoch.ne.0.0_R8)
     &         .or.(fch_inert.ne.0.0_R8.and.maxval(facdrift).ne.0.0_R8)
     &         .or.(vis_q.ne.0.0_R8.and.maxval(fac_vis).ne.0.0_R8)
     &         .or.(vis_par.ne.0.0_R8.and.maxval(fac_vis).ne.0.0_R8)
     &         .or.(fch_ion_neutral.ne.0.0_R8.and.
     &                                (maxval(facdrift).ne.0.0_R8.or.
     &                                 maxval(fac_ExB).ne.0.0_R8))
     &         .or.(maxval(facdrift).ne.0.0_R8.and.dia_cur.ne.0.0_R8)
     &         .or.(vis_per.ne.0.0_R8.and.maxval(fac_vis).ne.0.0_R8),
     &         'No radial current terms activated!')                    !xpb
       call xertst (0.0_R8.le.dia_cur .and. 1.0_R8.ge.dia_cur,          !xpb
     &  'faulty argument dia_cur')                                      !xpb
       call xertst (0.0_R8.le.vis_q .and. 1.0_R8.ge.vis_q,              !xpb
     &  'faulty argument vis_q')                                        !xpb
       call xertst (0.0_R8.le.vis_par .and. 1.0_R8.ge.vis_par,          !xpb
     &  'faulty argument vis_par')                                      !xpb
       call xertst (0.0_R8.le.vis_per .and. 1.0_R8.ge.vis_per,          !xpb
     &  'faulty argument vis_per')                                      !xpb
       call xertst (0.0_R8.le.fch_ion_neutral .and.                     !xpb
     &  1.0_R8.ge.fch_ion_neutral, 'faulty argument fch_ion_neutral')   !xpb
       call ipgeti ('b2tfhe_no_current', no_current)                    !xpb
       call ipgeti ('b2news_poteq', pot_eq)                             !srv 27.02.13
       call xertst (pot_eq.eq.1.or.no_current.eq.1,
     &  'If b2news_poteq.ne.1, then b2tfhe_no_current must be 1!')
       call ipgeti ('b2tfhe_conduction_only', conduction_only)          !xpb
       call ipgeti ('b2tfhe_no_hybr', no_hybr)                          !xpb
       call ipgeti ('b2tfch_iout', iout)                                !srv 18.06.08
       call ipgeti ('b2sigp_style', style)                              !srv 13.01.17
       call ipgetr ('b2sifr_phm2', phm2)                                !srv 15.10.18
       call ipgetr ('b2tfhe_prl_cur', prl_cur)                          !srv 15.10.18
      endif                                                             !xpb
*   ..set n2
      n2 = (nx+2)*(ny+2)
*   ..test nx, ny
      call xertst (0.le.nx.and.0.le.ny, 'faulty argument nx, ny')
*   ..test facdrift
      call xertst (0.0_R8.le.minval(facdrift) .and.
     &     1.0_R8.ge.maxval(facdrift), 'faulty argument facdrift')
*   ..extensive tests on first few calls
      if (ncall.lt.3) then
*    ..test edge values of fne
       call b2xvff (nx, ny, fne_32, 'fne_32')
       call b2xvff (nx, ny, fne_52, 'fne_52')
*    ..test sign of ne, te, csig, csig_an, calf_an, csigin, chce
       call b2xvsg (n2, ne, 1, 'ne', '.gt.')
       call b2xvsg (n2, te, 1, 'te', '.gt.')
       call b2xvsg (n2*2, csig, 1, 'csig', '.ge.')
       call b2xvsg (n2*2, csig_an, 1, 'csig_an', '.ge.')
       call b2xvsg (n2*2, calf_an, 1, 'calf_an', '.ge.')
       call b2xvsg (n2*2*ns*ns, csigin, 1, 'csigin', '.ge.')
       call b2xvsg (n2*2, chce, 1, 'chce', '.ge.')
*    ..test edge values of csig, calf, csig_an, calf_an, csigin, chce, chve, chvemx
       call b2xvff (nx, ny, csig, 'csig')
       call b2xvff (nx, ny, calf, 'calf')
       call b2xvff (nx, ny, csig_an, 'csig_an')
       call b2xvff (nx, ny, calf_an, 'calf_an')
       do is = 0, ns-1
         do is0 = 0, ns-1
           call b2xvff (nx, ny, csigin(-1,-1,0,is,is0), 'csigin')
         enddo
       enddo
       call b2xvff (nx, ny, chce, 'chce')
       call b2xvff (nx, ny, chve, 'chve')
       call b2xvfx (nx, ny, chvemx, 'chvemx')
      endif

c Drift terms
      facdriftm = maxval(facdrift)                                       !srv 23.12.19
      fac_ExBm  = maxval(fac_ExB)                                        !srv 23.12.19
*   ..calculate the contribution of ion-neutral current                  !srv 11.02.04 {
      call b2tinnt (nx, ny, ns, facdrift, fac_ExB, ismain, ismain0,
     &              vol, gs, hx, hy, qz, qc, bb,
     &              te, po, ti, na, fna, csigin, fchin)                  !srv 11.05.09
*   .. compute x,y-components of divergent part of diamagnetic current
      call b2tdia(nx, ny, ns, vol, gs, bb,
     & pz, fna, uadia, facdrift, fchdia)
*   .. calculate x,y-components of viscosity parallel current
      call b2tvspa(nx, ny, ns, bb, ua,
     & cvsa_cl, fac_vis, fchvispar)                                      !srv 20.09.06 16.10.17
*   .. calculate the contribution of viscosity perpendicular current     !srv 22.11.99
      call b2tvspr(nx, ny, ns, fac_vis, vol,                             !srv 10.10.17
     & po, ti, na, avisper, bvisper, cvisper, fchvisper)
*   .. calculate the parallel viscosity current connected with the heat flux
      call b2tvsq (nx, ny, ns, fac_vis, ! vol, hx, hy, qz, bb,           !srv 10.10.17
     & ne2, lnlam, te, ti, ni, chci, fchvisq)
*  ..  calculate inertia current                                         !srv 07.07.99
      call b2tiner (nx, ny, ns, vol, gs, bb, na, ua,                     !srv 22.11.12
     & facdrift, fchinert)
*  .. calculate the y-component of electron current associated with
*     stochastization of magnetic field
      call b2tstch (nx, ny, vol, po, te, ne, csig_stoch, fchstoch)       !srv 07.07.14
*   ..compute parallel current
      call b2tcpa (nx, ny, ns, te, po, ne, na, ua, csig, calf,           !srv 29.01.20
     & ehx_, zeff, fch_p)
* ..compute fluxes and linearisation
*  ..do not forget equivalent code in B2TFHE !!!!!!!!!!!!!!!!!!!!!
*   ..compute fchx, concx
      do iy = -1, ny
        do ix = -1, nx
          if(leftix(ix,iy).eq.-2) then
            conc(ix,iy,0)   = 0.0_R8
            fch(ix,iy,0)    = 0.0_R8
            fch_32(ix,iy,0) = 0.0_R8
            fch_52(ix,iy,0) = 0.0_R8
          else
!xpb        ehx = -(po(ix,iy)-po(leftix(ix,iy),leftiy(ix,iy)))+
!xpb     &         (ne(ix,iy)*te(ix,iy)
!xpb     &         -ne(leftix(ix,iy),leftiy(ix,iy))*
!xpb     &         te(leftix(ix,iy),leftiy(ix,iy)))*
!xpb     &         (2/qe)/(ne(leftix(ix,iy),leftiy(ix,iy))+ne(ix,iy))
            if (no_current.eq.1) then
              ehx = 0.0_R8
              ehx_po = 0.0_R8
            else
              ehx = ehx_(ix,iy)                                          !srv 13.01.17
              ehx_po = -(po(ix,iy)-po(leftix(ix,iy),leftiy(ix,iy)))      !xpb
            endif
            tex = te(ix,iy)-te(leftix(ix,iy),leftiy(ix,iy))
cxpb Changed the solution to (3/2+BoRiS) for benchmarking !!!!!!!!!!!!!!!!!!!
            if (conduction_only.eq.1) then
              floe0_nd = 0.0_R8
            else
              floe0_nd = (1.5_R8+BoRiS)*fne_32(ix,iy,0)+
     1         2.5_R8*fne_52(ix,iy,0)+
     &         chve(ix,iy,0)*
     1         (ne(leftix(ix,iy),leftiy(ix,iy))+ne(ix,iy))/2.0_R8
            endif
            floe0 = floe0_nd
            if (no_hybr.eq.0) then
              cone0    = hybr(floe0   ,chce(ix,iy,0))-abs(floe0/2)
              cone0_nd = hybr(floe0_nd,chce(ix,iy,0))-abs(floe0_nd/2)
            else
              cone0    = chce(ix,iy,0)
              cone0_nd = chce(ix,iy,0)
            endif
cwdk 10.03.16 unified treatment of chvemx in b2tlhe.F. Test for b2tqce_model no longer required
            flomx = chvemx(ix,iy)*(ne(leftix(ix,iy),leftiy(ix,iy))+
     &        ne(ix,iy))/2.0_R8
            if (flomx.gt.0.0_R8 .and. 
     &        lim_flux.eq.1 .and. conductive_limit.eq.0) then
              t0    = cone0   *abs(tex)/(flomx*max(
     &         te(leftix(ix,iy),leftiy(ix,iy)),te(ix,iy)))
              t0_nd = cone0_nd*abs(tex)/(flomx*max(
     &         te(leftix(ix,iy),leftiy(ix,iy)),te(ix,iy)))
            else
              t0    = 0.0_R8
              t0_nd = 0.0_R8
            endif
            tex0    = tex/(1.0_R8+t0)
            tex0_nd = tex/(1.0_R8+t0_nd)
            if (style.eq.0) then                                         !srv 15.10.18 {
              fch_p(ix,iy,0) =
     &         ((csig(ix,iy,0)-csig_an(ix,iy,0))*ehx-
     &         phm2*calf(ix,iy,0)*tex0)*prl_cur
              fch(ix,iy,0) = ((csig(ix,iy,0)-csig_an(ix,iy,0))*ehx-
     &         phm2*calf(ix,iy,0)*tex0)*prl_cur
     &       +fch_anomalous*(csig_an(ix,iy,0)*ehx_po-
     &                  phm2*calf_an(ix,iy,0)*tex0)*prl_cur
     &       +fchdia(ix,iy,0)+fchin(ix,iy,0)+fchvispar(ix,iy,0)          !xpb
     &       +fchvisper(ix,iy,0)
            else                                                         !srv 15.10.18 }
              fch(ix,iy,0) = fch_p(ix,iy,0)                              !srv 13.01.17
     &       +fchdia(ix,iy,0)+fchin(ix,iy,0)+fchvispar(ix,iy,0)          !xpb
     &       +fchvisper(ix,iy,0)
            endif                                                        !srv 15.10.18
            fch_32(ix,iy,0) = csig_cl(ix,iy,0)*ehx
     1       +fchin(ix,iy,0)+fchvispar(ix,iy,0)
     &       +fchvisper(ix,iy,0)
            if (style.eq.0) then
              fch_52(ix,iy,0) =
     &           fch_anomalous*(csig_an(ix,iy,0)*ehx_po
     &                         -calf_an(ix,iy,0)*tex0_nd)
     &                         -calf(ix,iy,0)*tex0_nd
            else
              fch_52(ix,iy,0) = csig_an(ix,iy,0)*ehx_po
     &                         -calf(ix,iy,0)*tex0_nd
            end if
            conc(ix,iy,0) = csig(ix,iy,0)
            if(fch_ion_neutral.ne.0.0_R8 .and.                           !xpb
     &        (facdriftm.gt.0.0_R8 .or. fac_ExBm.gt.0.0_R8)) then        !srv 23.12.19
              do is = 0, ns-1
                if (.not.is_neutral(is) .and. nint(zn(is)).eq.1) then
                  do is0 = 0, ns-1
                    if (is_neutral(is0) .and. nint(zn(is0)).eq.1)
     >               conc(ix,iy,0) =
     &                conc(ix,iy,0) + csigin(ix,iy,0,is,is0)
                  enddo
                endif
              enddo
            endif                                                        !xpb
          endif
        enddo
      enddo
*   ..compute fchy, concy
      call b2xehy (nx, ny, te, po, ne, ehy_)                             !srv 13.01.17
      do iy = -1, ny
        do ix = -1, nx
          if (bottomiy(ix,iy).eq.-2) then
            conc(ix,iy,1)   = 0.0_R8
            fch(ix,iy,1)    = 0.0_R8
            fch_32(ix,iy,1) = 0.0_R8
            fch_52(ix,iy,1) = 0.0_R8
          else
            if (no_current.eq.1) then
              ehy = 0.0_R8
              ehy_po = 0.0_R8
            else
              ehy = ehy_(ix,iy)                                          !srv 13.01.17
              ehy_po = -(po(ix,iy)-po(bottomix(ix,iy),bottomiy(ix,iy)))  !xpb
            endif
            fch(ix,iy,1)= csig_cl(ix,iy,1)*ehy
     .       +csig_an(ix,iy,1)*ehy_po
     .       -calf(ix,iy,1)*
     .       (te(ix,iy)-te(bottomix(ix,iy),bottomiy(ix,iy)))             !xpb
     .       +fchdia(ix,iy,1)+fchin(ix,iy,1)+fchvispar(ix,iy,1)          !xpb
     &       +fchvisper(ix,iy,1)
     8       +fchstoch(ix,iy,1)                                          !lk 31.01.08 24.08.21
            fch_32(ix,iy,1)= csig_cl(ix,iy,1)*ehy
     .       +fchin(ix,iy,1)+fchvispar(ix,iy,1)                          !xpb
     &       +fchvisper(ix,iy,1)
            fch_52(ix,iy,1)=csig_an(ix,iy,1)*ehy_po-calf(ix,iy,1)*
     .       (te(ix,iy)-te(bottomix(ix,iy),bottomiy(ix,iy))) 
            conc(ix,iy,1) = csig(ix,iy,1)
            if(fch_ion_neutral.ne.0.0_R8 .and.                           !xpb
     &        (facdriftm.gt.0.0_R8 .or. fac_ExBm.gt.0.0_R8)) then        !srv 23.12.19
              do is = 0, ns-1
                if (.not.is_neutral(is)) then
                  do is0 = 0, ns-1
                    if (is_neutral(is0)) conc(ix,iy,1) =
     .               conc(ix,iy,1) + csigin(ix,iy,1,is,is0)              !xpb
                  enddo
                endif
              enddo
            endif
            if(fch_stochastic.ne.0.0e0_R8 ) then                         !lk 31.01.08 {
              conc(ix,iy,1) = conc(ix,iy,1) + csig_stoch(ix,iy,1)
            endif                                                        !lk 31.01.08 }
          endif                                                          !xpb
        enddo
      enddo

      if(guard_flows.eq.0) then                                          !xpb
cxpb turn off flows between boundary cells
        do iy = -1, ny
          do ix = -1, nx
            if(leftix(ix,iy).eq.-2 .or. rightix(ix,iy).eq.nx+1) then
              conc(ix,iy,1)=0.0_R8
              fch(ix,iy,1)=0.0_R8
              fch_32(ix,iy,1)=0.0_R8
              fch_52(ix,iy,1)=0.0_R8
              fchdia(ix,iy,1)=0.0_R8
              fchin(ix,iy,1)=0.0_R8
              fchvispar(ix,iy,1)=0.0_R8
              fchvisper(ix,iy,1)=0.0_R8
              fchvisq(ix,iy,1)=0.0_R8
            endif
          enddo        
        enddo

cxpb turn off poloidal flows in core boundary
        if (boundary_namelist.eq.1) then
          do ib = 1, nbc
            if(core_boundary(ib)) then
              do ix = 1, bc_list_size(ib)
                conc(bc_list_x(ix,ib),bc_list_y(ix,ib),0) = 0.0_R8
                fch(bc_list_x(ix,ib),bc_list_y(ix,ib),0) = 0.0_R8
                fch_32(bc_list_x(ix,ib),bc_list_y(ix,ib),0) = 0.0_R8
                fch_52(bc_list_x(ix,ib),bc_list_y(ix,ib),0) = 0.0_R8
                fch_p(bc_list_x(ix,ib),bc_list_y(ix,ib),0) = 0.0_R8      !srv 10.06.08
                fchdia(bc_list_x(ix,ib),bc_list_y(ix,ib),0) = 0.0_R8
                fchin(bc_list_x(ix,ib),bc_list_y(ix,ib),0) = 0.0_R8
                fchvispar(bc_list_x(ix,ib),bc_list_y(ix,ib),0) = 0.0_R8
                fchvisper(bc_list_x(ix,ib),bc_list_y(ix,ib),0) = 0.0_R8
                fchvisq(bc_list_x(ix,ib),bc_list_y(ix,ib),0) = 0.0_R8
              enddo
            endif
          enddo
        else
          do ix = -1, nx
            if (on_closed_surface(ix,-1)) then
              conc(ix,-1,0) = 0.0_R8
              fch(ix,-1,0) = 0.0_R8
              fch_32(ix,-1,0) = 0.0_R8
              fch_52(ix,-1,0) = 0.0_R8
              fch_p(ix,-1,0) = 0.0_R8
              fchdia(ix,-1,0) = 0.0_R8
              fchin(ix,-1,0) = 0.0_R8
              fchvispar(ix,-1,0) = 0.0_R8
              fchvisper(ix,-1,0) = 0.0_R8
              fchvisq(ix,-1,0) = 0.0_R8
            endif
          enddo
        endif
      endif                                                              !xpb

cxpb turn off current terms if requested
      if (no_current.eq.1) then
        conc=0.0_R8
        fch=0.0_R8
        fch_p=0.0_R8
        fch_32=0.0_R8
        fch_52=0.0_R8
        fchdia=0.0_R8
        fchin=0.0_R8
        fchvispar=0.0_R8
        fchvisper=0.0_R8
        fchvisq=0.0_R8
      endif

*  INTERIOR BOUNDARY AT DEAD REGIONS
      DO IB=1,NNISO
        DO IX=NXISO1(IB), NXISO2(IB)
          DO IY=NYISO1(IB), NYISO2(IB)
            conc(ix,iy,0)=0.0_R8
            fch(ix,iy,0)=0.0_R8
            fch_32(ix,iy,0)=0.0_R8
            fch_52(ix,iy,0)=0.0_R8
            fchdia(ix,iy,0)=0.0_R8
            fch_p(ix,iy,:) =0.0_R8
            fchvispar(ix,iy,0)=0.0_R8
            fchvisper(ix,iy,0)=0.0_R8
            fchvisq(ix,iy,0)=0.0_R8
            fchin(ix,iy,0)=0.0_R8
            conc(ix,iy,1)=0.0_R8
            fch(ix,iy,1)=0.0_R8
            fch_32(ix,iy,1)=0.0_R8
            fch_52(ix,iy,1)=0.0_R8
            fchdia(ix,iy,1)=0.0_R8
            fchin(ix,iy,1)=0.0_R8
            fchvispar(ix,iy,1)=0.0_R8
            fchvisper(ix,iy,1)=0.0_R8
            fchvisq(ix,iy,1)=0.0_R8
          end do
        end do
      end do

      if(iout.ne.0) then                                             !srv 18.06.08 {
       call my_out(70,nx,ny,fch(-1,-1,0),'b2tfch_fchx')
       call my_out(70,nx,ny,fch(-1,-1,1),'b2tfch_fchy')
       call my_out(70,nx,ny,fchin(-1,-1,0),'b2tfch_chinx')
       call my_out(70,nx,ny,fchin(-1,-1,1),'b2tfch_chiny')
       call my_out(70,nx,ny,fchvispar(-1,-1,0),'b2tfch_chvisparx')
       call my_out(70,nx,ny,fchvispar(-1,-1,1),'b2tfch_chvispary')
!       call my_out(70,nx,ny,fchvisq(-1,-1,0),'b2tfch_chvisqx')
!       call my_out(70,nx,ny,fchvisq(-1,-1,1),'b2tfch_chvisqy')
       call my_out(70,nx,ny,fchvisper(-1,-1,0),'b2tfch_chvisperx')
       call my_out(70,nx,ny,fchvisper(-1,-1,1),'b2tfch_chvispery')
!       call my_out(70,nx,ny,fchinert(-1,-1,0),'b2tfch_chinertx')
!       call my_out(70,nx,ny,fchinert(-1,-1,1),'b2tfch_chinerty')
       call my_out(70,nx,ny,fchdia(-1,-1,0),'b2tfch_chdiax')
       call my_out(70,nx,ny,fchdia(-1,-1,1),'b2tfch_chdiay')
!       call my_out(70,nx,ny,fchanml(-1,-1,0),'b2tfch_chanmlx')
!       call my_out(70,nx,ny,fchanml(-1,-1,1),'b2tfch_chanmly')
       call my_out(70,nx,ny,conc(-1,-1,0),'b2tfch_concx')
       call my_out(70,nx,ny,conc(-1,-1,1),'b2tfch_concy')
       call my_out(70,nx,ny,fch_p(-1,-1,0),'b2tfch_ch_px')
       call my_out(70,nx,ny,fch_p(-1,-1,1),'b2tfch_ch_py')
       call my_out(70,nx,ny,avisper(-1,-1,0),'b2tfch_avisperx')
       call my_out(70,nx,ny,avisper(-1,-1,1),'b2tfch_avispery')
       call my_out(70,nx,ny,bvisper(-1,-1,0),'b2tfch_bvisperx')
       call my_out(70,nx,ny,bvisper(-1,-1,1),'b2tfch_bvispery')
       call my_out(70,nx,ny,cvisper(-1,-1,0),'b2tfch_cvisperx')
       call my_out(70,nx,ny,cvisper(-1,-1,1),'b2tfch_cvispery')
       call my_out(70,nx,ny,fchstoch(-1,-1,1),'fchstochy')           !lk 31.01.08 24.08.21
      endif                                                          !srv 18.06.08 }

* ..return
      ncall = ncall+1
      call subend ()
      return

*-----------------------------------------------------------------------
*.scribbles

cxpb  This version of the routine uses a (1.5+BoRiS) factor for the
cxpb  poloidal convective heat flux instead of the usual 3/2 or 5/2 
cxpb  so as to solve the same equation as in BoRiS for benchmark purposes. 
cxpb  The default is BoRiS = 0.0.
cxpb  For comparison, set BoRiS = 1.0.
cxpb  Routine rewritten to explicitly separate 3/2 and 5/2 terms.
cxpb  The contributions from fne_32 and fne_52 are NO DRIFTS
cxpb  contributions, i.e. diamagnetic and ExB terms are treated
cxpb  separately. By the same token, the computed fch_32 and fch_52
cxpb  terms are also NO DRIFTS currents, i.e. without contribution
cxpb  from diamagnetic or ExB terms.

*-----------------------------------------------------------------------
*.end b2tfch

      end subroutine b2tfch
