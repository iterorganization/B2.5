*-----------------------------------------------------------------------
*.specification

      subroutine b2tfch_ (nx, ny, ns, nscx, ismain, ismain0, iscx,       !srv 09.01.20 {
     & facdrift, fac_ExB, fac_vis,
     1 vol, gs, hx, hy, hy1, hz, qs, qz, qc, bb,
     1 te, po, ne, ne2, lnlam, pz, ti, ni, na, fna, uadia, ua,
     1 csig, csigin, csig_stoch,                                         !srv 24.08.21
     & chci, calf, cvsa_cl, cvsa_drho,                                   !som 18.08.21
     & chci_al_ast, chci_a, cvsa_hAdp_albe, cvsa_hBdp_al,                !som 20.10.21
     & cvsa_uAdp_albe, cvsa_uBdp_al,                                     !som 26.04.24
     & conc, avisper_a, bvisper_a, cvisper_a,
     & fchanml_a, fchinert_a, fchvispar_a, fchvisper_a, fchvisq_a,       !srv 08.09.21
     & fchin_a, fch, fch_p, fchdia, fchin,
     & fchvispar, fchvisper, fchvisq, fchinert, fchanml, fchstoch)
      use b2mod_types
      use b2mod_boundary_namelist
      use b2mod_indirect
      use b2mod_constants
      use b2mod_anomalous_transport
      use b2mod_geo_corner
      use b2mod_rates
      use b2mod_b2cmfs
      use b2mod_b2cmpa
      use b2mod_subsys
      implicit none
*   ..input arguments (unchanged on exit)
      integer nx, ny, ns, nscx, ismain, ismain0, iscx(0:nscxmax-1)
      real (kind=R8) ::
     *   te(-1:nx,-1:ny), po(-1:nx,-1:ny), ne(-1:nx,-1:ny),
     &   vol(-1:nx,-1:ny), fna(-1:nx,-1:ny,0:1,0:1,0:ns-1),
     &   gs(-1:nx,-1:ny,0:2), ne2(-1:nx,-1:ny), pz(-1:nx,-1:ny),         !srv 11.07.99 18.01.02
     &   uadia(-1:nx,-1:ny,0:1,0:1,0:ns-1), ua(-1:nx,-1:ny,0:ns-1),      !srv 17.05.99 25.01.02
     &   csig(-1:nx,-1:ny,0:1,0:1),
     &   csigin(-1:nx,-1:ny,0:1,0:1,0:ns-1,0:ns-1),
     &   chci(-1:nx,-1:ny,0:1,0:1), calf(-1:nx,-1:ny,0:1,0:1),           !srv 11.07.99 24.06.08
     &   cvsa_cl(-1:nx,-1:ny,0:1,0:1,0:ns-1), na(-1:nx,-1:ny,0:ns-1),    !srv 16.10.17
     &   cvsa_drho(-1:nx,-1:ny,0:1,0:1,0:ns-1),                          !som 18.08.21
     &   chci_al_ast(-1:nx,-1:ny,0:1,0:1,0:nspecies-1),                  !som 20.10.21
     &   chci_a(-1:nx,-1:ny,0:1,0:1,0:ns-1),                             !som 20.10.21
     &   cvsa_hAdp_albe(-1:nx,-1:ny,0:1,0:1,0:nspecies-1,0:nspecies-1),  !som 20.10.21
     &   cvsa_hBdp_al(-1:nx,-1:ny,0:1,0:1,0:nspecies-1),                 !som 20.10.21
     &   cvsa_uAdp_albe(-1:nx,-1:ny,0:1,0:1,0:nspecies-1,0:nspecies-1),  !som 26.04.24
     &   cvsa_uBdp_al(-1:nx,-1:ny,0:1,0:1,0:nspecies-1),                 !som 26.04.24
     &   lnlam(-1:nx,-1:ny), qs(-1:nx,-1:ny,0:1),                        !srv 13.10.06 01.07.09
     1   hx(-1:nx,-1:ny), hy(-1:nx,-1:ny), hy1(-1:nx,-1:ny),
     1   hz(-1:nx,-1:ny), qz(-1:nx,-1:ny,0:1), qc(-1:nx,-1:ny,0:1),
     1   bb(-1:nx,-1:ny,0:3), ti(-1:nx,-1:ny),
     2   ni(-1:nx,-1:ny), facdrift(-1:nx,-1:ny), fac_ExB(-1:nx,-1:ny),
     2   fac_vis(-1:nx,-1:ny)
      real (kind=R8) :: csig_stoch(-1:nx,-1:ny,0:1,0:1)                  !lk 31.01.08 { 24.08.21
*   ..output arguments (unspecified on entry)
      real (kind=R8) ::
     &   fchanml_a(-1:nx,-1:ny,0:1,0:1,0:ns-1),                          !srv 08.09.21 {
     &   fchinert_a(-1:nx,-1:ny,0:1,0:1,0:ns-1),
     &   fchin_a(-1:nx,-1:ny,0:1,0:1,0:ns-1),
     &   fchvispar_a(-1:nx,-1:ny,0:1,0:1,0:ns-1),
     &   fchvisper_a(-1:nx,-1:ny,0:1,0:1,0:ns-1),
     &   fchvisq_a(-1:nx,-1:ny,0:1,0:1,0:ns-1),                          !srv 08.09.21 }
     *   fch(-1:nx,-1:ny,0:1,0:1), fch_p(-1:nx,-1:ny,0:1,0:1),
     &   fchstoch(-1:nx,-1:ny,0:1,0:1), fchvispar(-1:nx,-1:ny,0:1,0:1),  !srv 24.08.21
     &   fchvisper(-1:nx,-1:ny,0:1,0:1), fchin(-1:nx,-1:ny,0:1,0:1),     !srv 22.11.99 30.09.02
     &   conc(-1:nx,-1:ny,0:1,0:1), fchvisq(-1:nx,-1:ny,0:1,0:1),        !srv 11.07.99
     2   fchdia(-1:nx,-1:ny,0:1,0:1), fchinert(-1:nx,-1:ny,0:1,0:1),     !srv 07.07.99
     3   avisper_a(-1:nx,-1:ny,0:1,0:ns-1),
     &   bvisper_a(-1:nx,-1:ny,0:1,0:ns-1),                              !srv 27.01.00
     &   cvisper_a(-1:nx,-1:ny,0:1,0:ns-1), fchanml(-1:nx,-1:ny,0:1,0:1) !srv 24.03.07
*   ..local variables
*   ..common blocks
*-----------------------------------------------------------------------
*.documentation
*
*  1. purpose
*
*     B2TFCH_ computes the electric current between neighboring cells,
*     and the coefficients that specify the linearized response to
*     changes in electric potential.
*
*-----------------------------------------------------------------------
*.declarations

*   ..local variables
      integer ncall, n2, ix, iy, ib, is, is0
      integer, save :: iout = 0, iout_b2wdat = 0,
     &                 iout_b2npp7 = 0, iout_b2nppo = 0
      integer, save :: pot_eq = 1
      integer, save :: no_current = 0
      integer, save :: use_9pt_stencil = 1

      character(3) :: chns, chns0
      real (kind=R8), save :: fch_stochastic = 0.0e0_R8                  !lk 31.01.08 }
      real (kind=R8) ::
     &  dia_cur, vis_par, vis_q,                                         !srv 17.05.99 11.07.99
     &  vis_per, sig_stoch                                               !srv 22.11.99 13.10.06 17.06.09
      real (kind=R8) ::
     *    fch_ion_neutral, fch_anomalous, fch_inert,                     !srv 30.09.02 07.07.99
     &    ehx_(-1:nx,-1:ny,0:1), zeff(-1:nx,-1:ny)                       !srv 13.01.17
      real (kind=R8) :: facdriftm, fac_ExBm, t0                          !srv 23.12.19

      save ncall, fch_ion_neutral, fch_anomalous, dia_cur, fch_inert,    !srv 07.07.99
     &     vis_par, vis_per, vis_q, sig_stoch                            !srv 17.05.99 11.07.99 22.11.99
*   ..procedures
      external xertst
      external b2xvsg, b2xvff, b2xvfx, b2ttia                            !srv 11.07.99
*   ..initialization
      data ncall/0/, fch_ion_neutral/0.0_R8/, fch_inert/1.0_R8/          !srv 07.07.99 22.03.02
     *        , fch_anomalous/1.0_R8/                                    !srv 22.03.02
     *        , dia_cur/1.0_R8/, vis_par/1.0_R8/,                        !srv 17.05.99 22.03.02
     *          vis_q/1.0_R8/, vis_per/0.0_R8/, sig_stoch/1.0_R8/        !srv 11.07.99 22.11.99 22.03.02

*-----------------------------------------------------------------------
*.computation

* ..preliminaries
*   ..subprogram start-up calls
      call subini ('b2tfch_')
*   ..test input
      call xertst (0.le.nx.and.0.le.ny, 'faulty argument nx, ny')
      call xertst (1.le.ns, 'faulty argument ns')
      call xertst (0.le.ismain.and.ismain.lt.ns,
     . 'invalid main plasma species index ismain')
      call xertst (0.le.ismain0.and.ismain0.lt.ns.and.
     & (is_neutral(ismain0).or.ismain0.eq.ismain),
     . 'invalid main neutral species index ismain0')
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
       call ipgetr ('b2tstch_sig', sig_stoch)
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
     &         'No radial current terms activated!')                     !xpb
       call ipgeti ('b2tfch_iout', iout)                                 !srv 18.06.08
       call ipgeti ('b2npp7_iout', iout_b2npp7)
       call ipgeti ('b2nppo_iout', iout_b2nppo)
       call ipgeti ('b2wdat_iout', iout_b2wdat)                          !iys 17.04.17
       call ipgeti ('b2tfhe_no_current', no_current)
       call xertst (pot_eq.eq.1.or.no_current.eq.1,
     &  'If b2news_poteq.ne.1, then b2tfhe_no_current must be 1!')
       call ipgeti ('b2mndr_use_9pt_stencil', use_9pt_stencil)
       call xertst (.not.(use_9pt_stencil.ne.0.and.vis_per.ne.0.0_R8),
     &  'Cannot request 9-pt stencil and perp. viscosity together!')
      endif
*   ..set n2
      n2 = (nx+2)*(ny+2)
*   ..test nx, ny
      call xertst (0.le.nx.and.0.le.ny, 'faulty argument nx, ny')
*   ..extensive tests on first few calls
      if (ncall.lt.3) then
*    ..test sign of ne, te, csig, chce
       call b2xvsg (n2, ne, 1, 'ne', '.gt.')
       call b2xvsg (n2, te, 1, 'te', '.gt.')
       call b2xvsg (n2*2, csig, 1, 'csig', '.ge.')
*    ..test edge values of csig, calf
       call b2xvff (nx, ny, csig, 'csig')
       call b2xvff (nx, ny, calf, 'calf')
      endif

      facdriftm = maxval(facdrift)
      fac_ExBm  = maxval(fac_ExB)                                        !srv 23.12.19
*   ..calculate the contribution of ion-neutral current                  !srv 11.02.04 {
      call b2tinnt (nx, ny, ns, nscx, iscx,
     &              facdrift, fac_ExB,
     &              vol, gs, hx, hy, hy1, hz, qz, qc, qs, bb,
     &              te, po, ti, na, fna, csigin, fchin_a, fchin)         !srv 11.05.09, IYS 17.10.2023
*   .. compute x,y-components of divergent part of diamagnetic current
      call b2tdia(nx, ny, ns, vol, gs, hx, hy, hy1, qz, qs, bb,
     &            pz, fna, uadia, facdrift, fchdia)
*   .. calculate x,y-components of viscosity parallel current
      call b2tvspa (nx, ny, ns, ismain, vol, hx, hy, hy1, qs, bb,        !srv 19.02.21
     &              na, ua, ti, cvsa_cl, cvsa_drho, fac_vis,             !som 18.08.21
     &              cvsa_uAdp_albe, cvsa_uBdp_al,                        !som 26.04.24
     &              fchvispar_a, fchvispar)
*   .. calculate the contribution of viscosity perpendicular current     !srv 22.11.99
      call b2tvspr(nx, ny, ns, ismain, fac_vis,                          !srv 10.10.17 19.02.21
     &             vol, hy1, po, ti, na,
     &             avisper_a, bvisper_a, cvisper_a,
     &             fchvisper_a, fchvisper)
*   .. calculate the parallel viscosity current connected with the heat flux
      call b2tvsq (nx, ny, ns, ismain, fac_vis,                          !srv 10.10.17 19.02.21
     &             ne2, lnlam, te, ti, ni, chci, na,                     !som 21.10.21
     &             chci_al_ast, chci_a, cvsa_hAdp_albe, cvsa_hBdp_al,    !som 20.10.21
     &             fchvisq_a, fchvisq)
*  .. calculate inertia current                                          !srv 07.07.99
      call b2tiner (nx, ny, ns, ismain, vol, gs, hx, hy, qs, bb,         !srv 22.11.12
     &              na, ua, facdrift, fchinert_a, fchinert)
*  .. calculate the x,y-components of anomalous current
      call b2tanml (nx, ny, ns, ismain, vol, hx, hy, qs,                 !som 02.11.21
     &              csig_an, po, ne, na, fchanml_a, fchanml)             !srv 08.09.21 }
*  .. calculate the y-component of electron current associated with
*     stochastization of magnetic field
      call b2tstch (nx, ny, vol, hx, hy, hy1, qc, qs, po, te, ne,        !srv 07.07.14
     &              csig_stoch, fchstoch)
*  .. compute parallel current
      call b2tcpa (nx, ny, ns, te, po, ne, na, ua, csig, calf,
     &  ehx_, zeff, fch_p)
*  .. compute total current without fchstoch, which will be added below
      fch = fch_p + fchin + fchdia + fchvisper + fchvispar +
     &      fchvisq + fchinert + fchanml

*   ..compute concx
      do iy = -1, ny
       do ix = -1, nx                                                    !srv 26.11.02 {
        if(leftix(ix,iy).eq.-2)then                                      !srv 26.11.02
         conc(ix,iy,0,0:1)    = 0.0_R8
        else                                                             !srv 26.11.02 }
         conc(ix,iy,0,0) = csig(ix,iy,0,0)
         conc(ix,iy,0,1) = csig(ix,iy,0,1)
         if(fch_ion_neutral.ne.0.0e0_R8 .and.
     &     (facdriftm.gt.0.0_R8 .or. fac_ExBm.gt.0.0_R8)) then           !srv 23.12.19
          do is = 0, ns-1                                                !srv 24.06.08 {
!           if (.not.is_neutral(is)) then
            if (.not.is_neutral(is).and.nint(zn(is)).eq.1) then          !srv 20.10.17
             do is0 = 0, ns-1
              if (is_neutral(is0).and.nint(zn(is0)).eq.1) then           !srv 20.10.17
               conc(ix,iy,0,0) = conc(ix,iy,0,0) +
     &            csigin(ix,iy,0,0,is,is0) !*facdrift                    !srv 07.02.03
     &          *(wbbl(ix,iy,2)/wbbl(ix,iy,3))**2                        !srv 15.02.11
!     1        *((bb(ix,iy,2)/bb(ix,iy,3))**2+                           !srv 25.02.03 {
!     2        (bb(leftix(ix,iy),leftiy(ix,iy),2)/
!     3        bb(leftix(ix,iy),leftiy(ix,iy),3))**2)*0.5e0_R8           !srv 25.02.03 }
              conc(ix,iy,0,1) = conc(ix,iy,0,1)+csigin(ix,iy,0,1,is,is0)
              endif                                                      !srv 20.10.17
             enddo
            endif
          enddo                                                          !srv 24.06.08 }
         endif                                                           !srv
         if(fch_stochastic.ne.0.0_R8 ) then                              !lk 31.01.08 {
          conc(ix,iy,0,1) = conc(ix,iy,0,1) + csig_stoch(ix,iy,0,1)
          fch(ix,iy,0,1) = fch(ix,iy,0,1) + fchstoch(ix,iy,0,1)
         endif                                                           !lk 31.01.08 }
        endif                                                            !srv 26.11.02
       enddo
      enddo

*   ..compute concy
      do ix = -1, nx
       conc(ix,-1,1,0:1) = 0.0_R8
      enddo
      do iy = 0, ny
       do ix = -1, nx
         conc(ix,iy,1,1) = csig(ix,iy,1,1)                               !srv
         if (use_9pt_stencil.ne.0) then
cwdk for hybrid 9-point scheme, need to have only "effective diffusion" in coni
c    (for multiplication with full radial difference), and floi will be used in
c    an upwind scheme
c    implementation ensures backwards compatibility for perfectly orthogonal grids
           t0 = -(hy(ix,iy) + hy(bottomix(ix,iy),bottomiy(ix,iy)))/
     &           (hx(ix,iy) + hx(bottomix(ix,iy),bottomiy(ix,iy)))*
     &           qs(ix,iy,1)
           conc(ix,iy,1,0) = csig(ix,iy,1,1)*t0
         endif
         if(fch_ion_neutral.ne.0.0_R8 .and.
     &     (facdriftm.gt.0.0_R8 .or. fac_ExBm.gt.0.0_R8)) then           !srv 23.12.19
         do is = 0, ns-1                                                 !srv 24.06.08 {
!          if (.not.is_neutral(is)) then
          if (.not.is_neutral(is).and.nint(zn(is)).eq.1) then            !srv 20.10.17
           do is0 = 0, ns-1
            if (is_neutral(is0).and.nint(zn(is0)).eq.1) then             !srv 20.10.17
            conc(ix,iy,1,1) = conc(ix,iy,1,1) + csigin(ix,iy,1,1,is,is0) !*facdrift !srv 27.03.99 07.02.03 31.01.08
            if (use_9pt_stencil.ne.0) conc(ix,iy,1,0) = conc(ix,iy,1,0)
     &        +csigin(ix,iy,1,1,is,is0)*t0
            endif                                                        !srv 20.10.17
           enddo
          endif
         enddo                                                           !srv 24.06.08 }
        endif                                                            !srv
        if(fch_stochastic.ne.0.0_R8 ) then                               !lk 31.01.08 {
          conc(ix,iy,1,1) = conc(ix,iy,1,1) + csig_stoch(ix,iy,1,1)
          if (use_9pt_stencil.ne.0)
     >     conc(ix,iy,1,0) = conc(ix,iy,1,0) + csig_stoch(ix,iy,1,1)*t0
          fch(ix,iy,1,1) = fch(ix,iy,1,1) + fchstoch(ix,iy,1,1)
        endif                                                            !lk 31.01.08 }
       enddo
      enddo

*  INTERIOR BOUNDARY AT DEAD REGIONS
      DO IB=1,NNISO
        DO IX=NXISO1(IB), NXISO2(IB)
          DO IY=NYISO1(IB), NYISO2(IB)
            fch(ix,iy,0,0:1)  = 0.0_R8
            conc(ix,iy,0,0:1) = 0.0_R8
          enddo
        enddo
        DO IX=NXISO1(IB), NXISO2(IB)
          DO IY=NYISO1(IB), NYISO2(IB)
            fch(ix,iy,1,0:1)  = 0.0_R8
            conc(ix,iy,1,0:1) = 0.0_R8
          end do
        enddo
      enddo

      if(iout.ne.0) then                                                 !srv 18.06.08 {
       call my_out(70,nx,ny,fch(-1,-1,0,0),'b2tfch__fchx')
       call my_out(70,nx,ny,fch(-1,-1,1,1),'b2tfch__fchy')
       call my_out(70,nx,ny,fchin(-1,-1,0,0),'b2tfch__fchinx')
       call my_out(70,nx,ny,fchin(-1,-1,1,1),'b2tfch__fchiny')
       call my_out(70,nx,ny,fchvispar(-1,-1,0,0),'b2tfch__fchvisparx')
       call my_out(70,nx,ny,fchvispar(-1,-1,1,1),'b2tfch__fchvispary')
       call my_out(70,nx,ny,fchvisq(-1,-1,0,0),'b2tfch__fchvisqx')
       call my_out(70,nx,ny,fchvisq(-1,-1,1,1),'b2tfch__fchvisqy')
       call my_out(70,nx,ny,fchvisper(-1,-1,0,0),'b2tfch__fchvisperx')
       call my_out(70,nx,ny,fchvisper(-1,-1,1,1),'b2tfch__fchvispery')
       call my_out(70,nx,ny,fchinert(-1,-1,0,0),'b2tfch__fchinertx')
       call my_out(70,nx,ny,fchinert(-1,-1,1,1),'b2tfch__fchinerty')
       call my_out(70,nx,ny,fchdia(-1,-1,0,0),'b2tfch__fchdiax')
       call my_out(70,nx,ny,fchdia(-1,-1,1,1),'b2tfch__fchdiay')
       call my_out(70,nx,ny,fchanml(-1,-1,0,0),'b2tfch__fchanmlx')
       call my_out(70,nx,ny,fchanml(-1,-1,1,1),'b2tfch__fchanmly')
       call my_out(70,nx,ny,fchstoch(-1,-1,1,1),'b2tfch__fchstochy')     !lk 31.01.08 24.08.21
       call my_out(70,nx,ny,conc(-1,-1,0,0),'b2tfch__concx')
       call my_out(70,nx,ny,conc(-1,-1,1,1),'b2tfch__concy')
       call my_out(70,nx,ny,fch_p(-1,-1,0,0),'b2tfch__fch_px')
       call my_out(70,nx,ny,fch_p(-1,-1,1,1),'b2tfch__fch_py')
       call my_out(70,nx,ny,csig_an(-1,-1,0,0),'b2tfch__csig_anx')
       call my_out(70,nx,ny,csig_an(-1,-1,1,1),'b2tfch__csig_any')
       call my_out(70,nx,ny,csig(-1,-1,0,0),'b2tfch__csigx')
       call my_out(70,nx,ny,csig(-1,-1,1,1),'b2tfch__csigy')
       call my_out(70,nx,ny,calf(-1,-1,0,0),'b2tfch__calfx')
       call my_out(70,nx,ny,calf(-1,-1,1,1),'b2tfch__calfy')
       call my_out(70,nx,ny,te,'b2tfch__te')
       call my_out(70,nx,ny,ne,'b2tfch__ne')
       call my_out(70,nx,ny,po,'b2tfch__po')
       do is = 0, ns-1
        write (chns,'(i3.3)') is
        call my_out(70,nx,ny,avisper_a(-1,-1,0,is),
     &  'b2tfch__avisper_ax'//chns)
        call my_out(70,nx,ny,avisper_a(-1,-1,1,is),
     &  'b2tfch__avisper_ay'//chns)
        call my_out(70,nx,ny,bvisper_a(-1,-1,0,is),
     &  'b2tfch__bvisper_ax'//chns)
        call my_out(70,nx,ny,bvisper_a(-1,-1,1,is),
     &  'b2tfch__bvisper_ay'//chns)
        call my_out(70,nx,ny,cvisper_a(-1,-1,0,is),
     &  'b2tfch__cvisper_ax'//chns)
        call my_out(70,nx,ny,cvisper_a(-1,-1,1,is),
     &  'b2tfch__cvisper_ay'//chns)
        call my_out(70,nx,ny,na(-1,-1,is),'b2tfch__na'//chns)
        do is0 = 0, ns-1
         write (chns0,'(i3.3)') is0
         call my_out(70,nx,ny,csigin(-1,-1,0,0,is,is0),'b2tfch__csiginx'
     &        //chns//'_'//chns0)
         call my_out(70,nx,ny,csigin(-1,-1,1,1,is,is0),'b2tfch__csiginy'
     &        //chns//'_'//chns0)
        enddo
       enddo
       call my_out(70,nx,ny,ehx_,'b2tfch__ehx_')                         !srv 13.01.17
      endif                                                              !srv 18.06.08 }

      if(iout_b2npp7.ne.0 .or. iout_b2nppo.ne.0) then
!       call my_out(70,nx,ny,fch(-1,-1,0,0),'b2tfch__fchx')
!       call my_out(70,nx,ny,fch(-1,-1,1,1),'b2tfch__fchy')
       call my_out(70,nx,ny,fchin(-1,-1,0,0),'b2tfch__for_po_fchinx')
       call my_out(70,nx,ny,fchin(-1,-1,1,1),'b2tfch__for_po_fchiny')
       call my_out(70,nx,ny,fchvispar(-1,-1,0,0),
     &   'b2tfch__for_po_fchvisparx')
       call my_out(70,nx,ny,fchvispar(-1,-1,1,1),
     &   'b2tfch__for_po_fchvispary')
!       call my_out(70,nx,ny,fchvisq(-1,-1,0,0),
!     &   'b2tfch__for_po_fchvisqx')
       call my_out(70,nx,ny,fchvisq(-1,-1,1,1),
     &   'b2tfch__for_po_fchvisqy')
!       call my_out(70,nx,ny,fchvisper(-1,-1,0,0),
!     &    'b2tfch__for_po_fchvisperx')
       call my_out(70,nx,ny,fchvisper(-1,-1,1,1),
     &    'b2tfch__for_po_fchvispery')
       call my_out(70,nx,ny,fchinert(-1,-1,0,0),
     &    'b2tfch__for_po_fchinertx')
       call my_out(70,nx,ny,fchinert(-1,-1,1,1),
     &    'b2tfch__for_po_fchinerty')
       call my_out(70,nx,ny,fchdia(-1,-1,0,0),'b2tfch__for_po_fchdiax')
       call my_out(70,nx,ny,fchdia(-1,-1,1,1),'b2tfch__for_po_fchdiay')
       call my_out(70,nx,ny,fchanml(-1,-1,0,0),
     &    'b2tfch__for_po_fchanmlx')
       call my_out(70,nx,ny,fchanml(-1,-1,1,1),
     &    'b2tfch__for_po_fchanmly')
       call my_out(70,nx,ny,fch_p(-1,-1,0,0),'b2tfch__for_po_fch_px')
!      call my_out(70,nx,ny,fch_p(-1,-1,1,1),'b2tfch__for_po_fch_py')
       call my_out(70,nx,ny,csig_an(-1,-1,0,0),
     &    'b2tfch__for_po_csig_anx')
       call my_out(70,nx,ny,csig_an(-1,-1,1,1),
     &    'b2tfch__for_po_csig_any')
       call my_out(70,nx,ny,fchstoch(-1,-1,1,1),
     &    'b2tfch__for_po_fchstochy')                                    !lk 31.01.08 24.08.21
      endif                                                              !srv 11.09.09 }

      if(iout_b2wdat.eq.4) then
!       call my_out(70,nx,ny,fch(-1,-1,0,0),'b2tfch__fchx')
!       call my_out(70,nx,ny,fch(-1,-1,1,1),'b2tfch__fchy')
       call my_out(70,nx,ny,fchin(-1,-1,0,0),'b2tfch__fchinx')
       call my_out(70,nx,ny,fchin(-1,-1,1,1),'b2tfch__fchiny')
       call my_out(70,nx,ny,fchvispar(-1,-1,0,0),'b2tfch__fchvisparx')
       call my_out(70,nx,ny,fchvispar(-1,-1,1,1),'b2tfch__fchvispary')
       call my_out(70,nx,ny,fchvisq(-1,-1,0,0),'b2tfch__fchvisqx')
       call my_out(70,nx,ny,fchvisq(-1,-1,1,1),'b2tfch__fchvisqy')
       call my_out(70,nx,ny,fchvisper(-1,-1,0,0),'b2tfch__fchvisperx')
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
       call my_out(70,nx,ny,fchstoch(-1,-1,1,1),'b2tfch__fchstochy')     !lk 31.01.08 24.08.21
      endif                                                              !srv 11.09.09 }

      if(iout.ne.0 .or. (iout_b2npp7.ne.0 .or. iout_b2nppo.ne.0) .or.    !srv 08.09.21 {
     &   iout_b2wdat.eq.4 ) then
       do is = 0, ns-1
        write (chns,'(i3.3)') is
        call my_out(70,nx,ny,fchanml_a(-1,-1,0,0,is),
     &                       'b2tfch__fchanml_ax'//chns)
        call my_out(70,nx,ny,fchanml_a(-1,-1,1,1,is),
     &                       'b2tfch__fchanml_ay'//chns)
        call my_out(70,nx,ny,fchinert_a(-1,-1,0,0,is),
     &                      'b2tfch__fchinert_ax'//chns)
        call my_out(70,nx,ny,fchinert_a(-1,-1,1,1,is),
     &                      'b2tfch__fchinert_ay'//chns)
        call my_out(70,nx,ny,fchvispar_a(-1,-1,0,0,is),
     &                      'b2tfch__fchvispar_ax'//chns)
        call my_out(70,nx,ny,fchvispar_a(-1,-1,1,1,is),
     &                      'b2tfch__fchvispar_ay'//chns)
        call my_out(70,nx,ny,fchvisper_a(-1,-1,0,0,is),
     &                      'b2tfch__fchvisper_ax'//chns)
        call my_out(70,nx,ny,fchvisper_a(-1,-1,1,1,is),
     &                      'b2tfch__fchvisper_ay'//chns)
        call my_out(70,nx,ny,fchvisq_a(-1,-1,0,0,is),
     &                      'b2tfch__fchvisq_ax'//chns)
        call my_out(70,nx,ny,fchvisq_a(-1,-1,1,1,is),
     &                      'b2tfch__fchvisq_ay'//chns)
        call my_out(70,nx,ny,fchin_a(-1,-1,0,0,is),
     &                      'b2tfch__fchin_ax'//chns)                    !IYS 17.10.2023
        call my_out(70,nx,ny,fchin_a(-1,-1,1,1,is),
     &                      'b2tfch__fchin_ay'//chns)                    !IYS 17.10.2023
        enddo
       endif                                                             !srv 08.09.21 }

* ..return
      ncall = ncall+1
      call subend ()
      return

*-----------------------------------------------------------------------
*.end b2tfch_

      end subroutine b2tfch_                                             !srv 09.01.20 }
