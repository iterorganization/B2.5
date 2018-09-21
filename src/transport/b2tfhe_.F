*-----------------------------------------------------------------------
*.specification

      subroutine b2tfhe_ (nCv, nFc, nVx, ns, ismain, ismain0, geo, mpg, 
     &  pl, dv, co, rt, balance)                                                          !djm Jan2017
      use b2mod_types
      use b2mod_boundary_namelist
      use b2mod_constants
!      use b2mod_anomalous_transport
!      use b2mod_geo_corner
      use b2mod_b2cmfs
      use b2mod_b2cmpa
      use b2us_geo
      use b2us_map
      use b2us_plasma
!WG_TODO      use b2mod_balance                                                  !djm Jan2017
!WG_TODO     & , only : fhe_32, fhe_52, fhe_thermj, fhe_cond, fhe_dia, fhe_ecrb,
!WG_TODO     &          fhe_strange, fhe_pschused, balance_netcdf
!WG_TODO      use b2mod_transport_fun                                            !srv 13.01.17
      implicit none
*   ..input arguments (unchanged on exit)
      integer nCv, nFc, nVx, ns, ismain, ismain0
      type (geometry), intent (in) :: geo
      type (mapping), intent (in) :: mpg
      type(B2Plasma), intent (in) :: pl 
      type(B2Rates), intent (in) :: rt
      real (kind=R8) :: cutlo
      parameter (cutlo = 1.0e-60_R8)
      logical balance                                                    !djm Jan2017
*   ..output arguments (unspecified on entry)
*   ..input/output arguments (unspecified on entry)
      type(B2Coeff), intent (inout) :: co
      type(B2Derivatives), intent (inout) :: dv
*   ..local variables
      real (kind=R8) :: upfac, hybfac !djm Jan2017
*   ..common blocks
*-----------------------------------------------------------------------
*.documentation
*
*  1. purpose
*
*     B2TFHE_ computes the electric current and the electron heat flux
*     between neighboring cells, and the coefficients that specify the
*     linearized response to changes in electric potential and electron
*     temperature.
*
*
*
*  8. further comments
*
*     In the evaluation of the thermo-electric term in the electron
*     heat flux, the expression harm(te(),te()) cannot be replaced
*     by the arithmetic average (this is an empirical result).
*
*-----------------------------------------------------------------------
*.declarations

*   ..local variables
      integer ncall, iCv, iFc, ib, is, ixl, iyl, ixb, iyb, is0
      integer, save :: iout = 0, iout_b2wdat = 0,                        !srv 16.09.02
     &                 iout_b2nph9 = 0, iout_b2npht = 0,
     &                 iout_b2npp7 = 0, iout_b2nppo = 0
      integer, save :: lim_flux = 0, mdf_fhe = 0                         !srv 04.02.03 13.10.06
      integer, save :: no_Ptncr_x_he = 0, no_Ptncr_y_he = 0              !srv 06.07.06
      integer, save :: pot_eq = 1                                        !srv 26.07.06
      integer, save :: conductive_limit = 1                              !srv 16.06.08 17.01.12
      integer, save :: style = 2                                         !srv 13.01.17
      integer, save :: no_current= 0                                     !srv 27.02.13
      integer, save :: meth = 2
csw 10oct2012
      integer, save :: b2_upwind = 0
      integer, save :: b2tfhe_hybr2 = 0, flo53 = 0

      character(3) :: charns, charns0
      real (kind=R8), save :: alfTeEh = 0.0_R8                           !srv 28.04.08 30.05.08
      real (kind=R8), save :: fch_pTe = 1.0_R8                           !srv 30.12.09
!WG_TODO      real (kind=R8) :: csig_stoch(-1:nx,-1:ny),                         !lk 31.01.08 {
!WG_TODO     &                  fchstoch(-1:nx,-1:ny)
      real (kind=R8), save :: fch_stochastic = 0.0e0_R8                  !lk 31.01.08 }
      real (kind=R8), save :: PSch=1.0e0_R8                              !srv 17.01.07
      real (kind=R8) ::
     *   ehx(nFc), ehy(nFc), floe0(nFc,0:1), flomx, alfe0(nFc),
     &  dia_cur, vis_par, vis_q, vdia_par,                    !srv 17.05.99 11.07.99
     &  vis_per, floe0_mdf(nFc,0:1)                          !srv 22.11.99 13.10.06 17.06.09
      real (kind=R8) ::
     *    fch_ion_neutral, fch_anomalous, fch_inert,                     !srv 30.09.02 07.07.99
     *    fledia(nFc,0:1), c071(nFc)!,                                   !srv 22.11.99 25.01.02
 !WG_TODO    &    zeff(-1:nx,-1:ny)                           !srv 13.01.17
      real (kind=R8) ::
     *    wrkv(nVx), tef(nFc), nef(nFc), tefh(nFc),
     *    dte(nFc,0:1), dpo(nFc,0:1), dpe(nFc,0:1), wght(nFc,2),
     *    fhe0_mdf(nFc,0:1), dumm0(nFc,0:1), dumm1(nFc,0:1)
      real (kind=R8) :: facdriftm
      real (kind=R8), allocatable ::
     *    w1(:,:), w2(:,:), w3(:,:), w4(:,:), w5(:,:), w6(:,:)           !srv 16.09.02

      save ncall, fch_ion_neutral, fch_anomalous, dia_cur, fch_inert,    !srv 07.07.99
     &     vis_par, vis_per, vis_q, vdia_par                             !srv 17.05.99 11.07.99 22.11.99
*   ..procedures
      intrinsic max, abs
      external subini, subend, xertst
      external b2xvsg, b2xvff, b2xvfx, b2ttia                            !srv 11.07.99
*   ..initialization
      data ncall/0/, fch_ion_neutral/0.0e0_R8/, fch_inert/1.0e0_R8/      !srv 07.07.99 22.03.02
     *        , fch_anomalous/1.0e0_R8/                                  !srv 22.03.02
     *        , dia_cur/1.0e0_R8/, vis_par/0.0e0_R8/, vdia_par/1.0e0_R8/ !srv 17.05.99 22.03.02
     *        , vis_q/1.0e0_R8/, vis_per/0.0e0_R8/                       !srv 11.07.99 22.11.99 22.03.02

*-----------------------------------------------------------------------
*.computation

* ..preliminaries
*   ..subprogram start-up calls
      call subini ('b2tfhe_')
*   ..set internal parameters on first call                              !srv
      if (ncall.eq.0) then                                               !srv
       call ipgetr ('b2tfhe_neutral', fch_ion_neutral)                   !srv
       call ipgetr ('b2tfhe_anomalous', fch_anomalous)                   !srv
       if (fch_anomalous.eq.1.0_R8)
     >  call ipgetr ('b2tanml_anomalous', fch_anomalous) ! for backward compatibility
       call ipgetr ('b2tfhe_dia_cur', dia_cur)                           !srv
       call ipgetr ('b2tfhe_vis_par', vis_par)                           !srv
       call ipgetr ('b2tfhe_vis_per', vis_per)                           !srv 22.11.99
       call ipgetr ('b2tfhe_vis_q', vis_q)                               !srv 11.07.99
       call ipgetr ('b2tfhe_vdia_par', vdia_par)                         !srv 17.05.99
       call ipgetr ('b2tiner_inert', fch_inert)                          !srv 07.07.99
       call ipgeti ('b2tfhe_lim_flux', lim_flux)                         !srv 04.02.03
       call ipgetr ('b2tfhe_stochastic', fch_stochastic)                 !lk 31.01.08
       call xertst(fch_ion_neutral .ne. 0.0e0_R8 .or.                    !srv
     *             fch_anomalous .ne. 0.0e0_R8 .or.
     *             vis_q .ne. 0.0_R8 .or. fch_inert .ne. 0.0_R8 .or.
     *             vis_per .ne. 0.0_R8 .or.
     *             fch_stochastic .ne. 0.0_R8 .or.
     *             dia_cur .ne. 0.0e0_R8 .or. vis_par .ne. 0.0_R8,       !srv
     *             'No radial current terms!')                           !srv
       call ipgeti ('b2tfhe_iout', iout)                                 !srv 16.09.02
       call ipgeti ('b2tfhe_no_Ptncr_x', no_Ptncr_x_he)                  !srv 06.07.06
       call ipgeti ('b2tfhe_no_Ptncr_y', no_Ptncr_y_he)                  !srv 06.07.06
       call ipgeti ('b2news_poteq', pot_eq)                              !srv 26.07.06
       call ipgeti ('b2tfhe_mdf', mdf_fhe)                               !srv 13.10.06
       call ipgetr ('b2tfhe_PSch', PSch)                                 !srv 17.01.07
       call ipgetr ('b2tfhe_alfTeEh', alfTeEh)                           !srv 28.04.08
       call ipgetr ('b2tfhe_fch_pTe', fch_pTe)                           !srv 30.12.09
       if(fch_pTe .eq. 1.0_R8) alfTeEh = 0.0_R8                          !srv 16.04.10
       call xertst (.not.(0.0_R8.ne.alfTeEh.and.0.0_R8.ne.fch_pTe),
     &        ' Do not set 0.0_R8.ne.alfTeEh.and.0.0_R8.ne.fch_pTe')
       call ipgeti ('b2nph9_iout', iout_b2nph9)
       call ipgeti ('b2npht_iout', iout_b2npht)
       call ipgeti ('b2npp7_iout', iout_b2npp7)
       call ipgeti ('b2nppo_iout', iout_b2nppo)
       call ipgeti ('b2wdat_iout', iout_b2wdat)                          !iys 17.04.17
       call ipgeti ('b2trcl_conductive_limit', conductive_limit)         !srv 16.06.08 17.01.12
       call ipgeti ('b2tfhe_no_current',no_current)
       call xertst (pot_eq.eq.1.or.no_current.eq.1,
     &  'If b2news_poteq.ne.1, then b2tfhe_no_current must be 1!')
       call ipgeti ('b2tfhe_upwind', b2_upwind)
       call ipgeti ('b2tfhe_hybr2', b2tfhe_hybr2)
       call ipgeti ('b2tfhe_flo53', flo53)
       call ipgeti ('b2sigp_style', style)                               !srv 13.01.17
       endif                                                             !srv
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
      dv%fchin     = 0.0_R8
      dv%fchdia    = 0.0_R8
      dv%fchvispar = 0.0_R8
      dv%fchvisper = 0.0_R8
      dv%fchvisq   = 0.0_R8
      dv%fchinert  = 0.0_R8
      dv%fhePSch   = 0.0_R8
!WG_TODO      avisper = 0.0_R8
!WG_TODO      bvisper = 0.0_R8
!WG_TODO      cvisper = 0.0_R8
!WG_TODO      fchstoch = 0.0_R8
      fledia = 0.0_R8 
#ifdef WG_TODO
*   ..calculate the contribution of ion-neutral current                  !srv 11.02.04 {
      call b2tinnt (nx, ny, ns, facdrift, fac_ExB, ismain, ismain0,
     &              vol, gs, hx, hy, qz, qc, qs, bb,
     &              te, po, ti, na, fna, csigin, fchin)                  !srv 11.05.09
c Drift terms
*   .. compute x,y-components of divergent part of diamagnetic current
      call b2tdia(nx, ny, ns, vol, gs, hx, hy, qz, qs, bb,
     &            pz, fna, uadia, facdrift, fchdia)
*   .. calculate x,y-components of viscosity parallel current
      call b2tvspa (nx, ny, ns, vol, hx, hy, qz, qs, bb, ua,
     &              cvsa_cl, fac_vis, fchvispar)                         !srv 20.09.06 16.10.17
*   .. calculate the contribution of viscosity perpendicular current     !srv 22.11.99
      call b2tvspr(nx, ny, ns, fac_vis, vol, hx, hy, qz, bb,             !srv 10.10.17
     &             po, ti, na, avisper, bvisper, cvisper, fchvisper)
*   .. calculate the parallel viscosity current connected with the heat flux
      call b2tvsq (nx, ny, ns, fac_vis, vol, hx, hy, qz, qs, bb, ne2,    !srv 10.10.17
     &             lnlam, te, ti, ni, chci, fchvisq)
*  ..  calculate inertia current                                         !srv 07.07.99
      call b2tiner (nx, ny, ns, vol, gs, hx, hy, qs, bb, na, ua,         !srv 22.11.12
     &              facdrift, fchinert)
#endif
*  .. calculate the x,y-components of anomalous current
      call b2tanml (nCv, nFc, ns, geo, mpg, co%csig_an,pl%po,dv%fchanml)
#ifdef WG_TODO
*  .. calculate the y-component of electron current associated with
*     stochastization of magnetic field
      call b2tstch (nx, ny, vol, hy, qz, bb, po, te, ne, csig_stoch,     !srv 07.07.14
     &              fchstoch)                                                 !srv 25.01.02 {
      if(facdriftm.ne.0.0e0_R8) then                                     !srv 18.01.07
       do is = 0,ns-1                                                    !srv 12.05.11
        if(nint(zn(is)).eq.1.and.nint(zamax(is)).eq.1) then              !srv 12.05.11
!         is = ismain   ! 1    ! for is = 1 (za = 1)                     !srv 12.05.11
* ..compute flediax
         do iy = -1, ny
          do ix = 0, nx
           if(leftix(ix,iy).ne.-2) then                                  !srv 26.11.02
            fledia(ix,iy,0,0) =
     &        (vol(leftix(ix,iy),leftiy(ix,iy))+vol(ix,iy))              !srv 25.01.02 was 5/2
     &       /(hx(leftix(ix,iy),leftiy(ix,iy))+hx(ix,iy))                !srv 17.05.99 11.03.00
     &       *vedia(ix,iy,0,0)*vdia_par                                  !srv 17.05.99 11.03.00 21.04.00 25.01.02 18.01.07
!     1       *(bb(ix,iy,2)/bb(ix,iy,3)+                                  !srv 13.09.02 {
!     2       bb(leftix(ix,iy),leftiy(ix,iy),2)/
!     3       bb(leftix(ix,iy),leftiy(ix,iy),3))/2.0e0_R8                 !srv 13.09.02 }
            fledia(ix,iy,0,1) =
     &       -(vol(leftix(ix,iy),leftiy(ix,iy))+vol(ix,iy))
     &       /(hx(leftix(ix,iy),leftiy(ix,iy))+hx(ix,iy))
     &       *vedia(ix,iy,0,1)*vdia_par*qs(ix,iy,0)/qc(ix,iy,0)
!     1       *(bb(ix,iy,2)/bb(ix,iy,3)+                                  !srv 13.09.02 {
!     2       bb(leftix(ix,iy),leftiy(ix,iy),2)/
!     3       bb(leftix(ix,iy),leftiy(ix,iy),3))/2.0e0_R8                 !srv 13.09.02 }
           endif                                                         !srv 26.11.02
          enddo
         enddo
* ..compute flediay
         do iy = 0, ny
          do ix = -1, nx
           fledia(ix,iy,1,1) =
     &       (vol(bottomix(ix,iy),bottomiy(ix,iy))+vol(ix,iy))           !srv 17.05.99 11.03.00 25.01.02 was 5/2
     &      /(hy1(bottomix(ix,iy),bottomiy(ix,iy))+hy1(ix,iy))           !srv 17.05.99 11.03.00
     &      *vedia(ix,iy,1,1)*vdia_par                                   !srv 17.05.99 11.03.00 18.01.07
!     1      *(bb(ix,iy,2)/bb(ix,iy,3)+                                   !srv 13.09.02 {
!     2      bb(bottomix(ix,iy),bottomiy(ix,iy),2)/
!     3      bb(bottomix(ix,iy),bottomiy(ix,iy),3))/2.0e0_R8              !srv 13.09.02 }
          enddo
         enddo                                                           !srv 25.01.02 }
        endif
       enddo
      endif                                                              !srv 18.01.07

      if(iout_b2nph9.ne.0 .or. iout_b2npht.ne.0
     &                         .or. iout_b2wdat.eq.4) then               !srv 16.09.02 {
       allocate(
     *    w1(-1:nx,-1:ny), w2(-1:nx,-1:ny), w3(-1:nx,-1:ny),
     *    w4(-1:nx,-1:ny), w5(-1:nx,-1:ny), w6(-1:nx,-1:ny)
     *         )
       w1=0.0e0_R8; w2=0.0e0_R8; w3=0.0e0_R8; w4=0.0e0_R8; w5=0.0e0_R8;
       w6=0.0e0_R8;
      endif                                                              !srv 16.09.02 }

      call b2tepsch(nx, ny, ns, hx, hy, qs, facdrift, cddi, te, ne, 
     .  cone, fhePSch)                                                   !srv 16.10.17

      fhe_mdf = 0.0e0_R8                                                 !srv 13.10.06
#endif
* ..compute parallel current
      call b2xehx (nCv, nFc, nVx, geo, mpg, pl%te, pl%po, dv%ne, ehx)    !srv 13.01.17 {
      call b2xehy (nCv, nFc, nVx, geo, mpg, pl%te, pl%po, dv%ne, ehy)    !srv 13.01.17
      call b2tcpa (nCv, nFc, nVx, geo, mpg, pl%te, pl%po, dv%ne, 
     &             co%csig, co%calf, ehx, dv%fch_p)
      if(style.eq.2) then
        call xerrab ('b2sigp_style 2 not yet supported')
#ifdef WG_TODO
        call b2xzef (nx, ny, ns, rz2, na, ne, zeff)
        call b2xpic (nx, ny, ns, qe, rza, bb, na, ua, zeff, fch_pi_c)
        call b2xpif (nx, ny, hx, hy, vol, fch_pi_c, fch_pi_f)
        do iy = -1, ny
          do ix = -1, nx
            c071(ix,iy) = fce2n(zeff(ix,iy))
          end do
        end do
#endif
      else
        dv%fch_pi_f = 0.0_R8                                           !srv 04.04.17
        c071 = 0.71e0_R8                                               !srv 13.01.17
      endif
      dv%fch_p = dv%fch_p - dv%fch_pi_f                                !srv 13.01.17 }

      if(pot_eq.ne.1 .and. no_current.ne.0) then                       !srv 27.02.13
        dv%fch_p = 0.0e0_R8
        ehx      = 0.0_R8                                              !srv 27.02.13
        ehy      = 0.0_R8                                              !srv 27.02.13
      endif                                                            !srv 30.12.09 }


*   ..interpolate ne and te to cell faces
      wght = 1.0_R8
      call intface(nCv,nFc,mpg%fcCv,wght,dv%ne,nef)
      call intface(nCv,nFc,mpg%fcCv,wght,pl%te,tef)

*   ..compute differences of te, po, and pe
      call diff(nCv, nFc, nVx, 0, geo, mpg, pl%te, wrkv, dte)
      call diff(nCv, nFc, nVx, 0, geo, mpg, pl%po, wrkv, dpo)
      call diff(nCv, nFc, nVx, 0, geo, mpg, pl%te*dv%ne, wrkv, dpe)

* ..compute fluxes and linearization
       floe0(:,0) = 1.5e0_R8*dv%fne(:,0) +
     &    (co%chve(:,0) + 2.5e0_R8*fledia(:,0))*nef                    !srv 13.09.02 18.01.07 11.01.13  *signmf
     &    -c071*dv%fch_p(:,0)/qe*fch_pTe                                  !srv 30.12.09 13.01.17
       floe0(:,1) = 1.5e0_R8*dv%fne(:,1) +                           
     &    (co%chve(:,1) + 2.5e0_R8*fledia(:,1))*nef !WG_TODO-                  !srv 11.05.99
!WG_TODO      &      5.0e0_R8/3.0e0_R8*fchstoch/qe)

       floe0_mdf(:,0) = 1.5e0_R8*dv%fne(:,0) +                          !srv 13.06.06 {
     &    co%chve(:,0)*nef - c071*dv%fch_p(:,0)/qe*fch_pTe                     !srv 30.12.09 13.01.17
       floe0_mdf(:,1) = 1.5e0_R8*dv%fne(:,1) +
     &    co%chve(:,1)*nef !WG_TODO-                  !srv 11.05.99
!WG_TODO      &      5.0e0_R8/3.0e0_R8*fchstoch/qe)

*   ..apply discretization scheme
      call calcflow(nCv,nFc,nVx,meth,geo,mpg,pl%te,floe0_mdf,
     &  co%chce,fhe0_mdf,dumm0,dumm1)

*   ..compute fluxes
      alfe0 = co%calf(:,0)
      tefh  = tef !WG_TODO: harmonic average to faces

      dv%fhe(:,0) = fhe0_mdf(:,0) +  (floe0(:,0) - floe0_mdf(:,0))*tef + 
     &  alfTeEh*alfe0*tefh*ehx
      dv%fhe(:,1) = fhe0_mdf(:,1) +  (floe0(:,1) - floe0_mdf(:,1))*tef + 
     &  alfTeEh*co%calf(:,1)*tefh*ehy

      dv%fhe_mdf(:,0)  = fhe0_mdf(:,0) - dv%fhePSch(:,0)  + 
     &  alfTeEh*alfe0*tefh*ehx
      dv%fhe_mdf(:,1)  = fhe0_mdf(:,1) - dv%fhePSch(:,1)  + 
     &  alfTeEh*co%calf(:,1)*tefh*ehy

*   ..compute floe, cone
      call calccoef(nCv,nFc,nVx,meth,geo,mpg,floe0_mdf,co%chce,
     &  dv%floe,dv%cone)

*   ..compute currents
       dv%fch = dv%fch_p + dv%fchin + dv%fchdia + dv%fchvispar + 
     &   dv%fchvisq + dv%fchinert +dv%fchanml 
       dv%fch(:,1) = dv%fch(:,1) + dv%fchvisper(:,1)!WG_TODO + fchstoch

*   ..compute conc
!WG_TODO: add contribution from csigin and cstoch first!
      dumm0 = 0.0e0_R8
      call calccoef(nCv,nFc,nVx,0,geo,mpg,dumm0,co%csig,dumm1,dv%conc)

#ifdef WG_TODO
         if(fch_ion_neutral.ne.0.0e0_R8 ) then                           !srv
          do is = 0, ns-1                                                !srv 24.06.08 {
!           if (.not.is_neutral(is)) then
            if (.not.is_neutral(is).and.nint(zn(is)).eq.1) then          !srv 20.10.17
             do is0 = 0, ns-1
              if (is_neutral(is0).and.nint(zn(is0)).eq.1) then           !srv 20.10.17
               conc(ix,iy,0,0) = conc(ix,iy,0,0) + 
     &            csigin(ix,iy,0,0,is,is0) !*facdrift                    !srv 07.02.03
     &          *(wbbl(ix,iy,2)/wbbl(ix,iy,3))**2                        !srv 15.02.11
              conc(ix,iy,0,1) = conc(ix,iy,0,1)+csigin(ix,iy,0,1,is,is0)
!     1        *((bb(ix,iy,2)/bb(ix,iy,3))**2+                           !srv 25.02.03 {
!     2        (bb(leftix(ix,iy),leftiy(ix,iy),2)/
!     3        bb(leftix(ix,iy),leftiy(ix,iy),3))**2)*0.5e0_R8           !srv 25.02.03 }
              endif                                                      !srv 20.10.17
             enddo
            endif
          enddo                                                          !srv 24.06.08 }
         endif                                                           !srv
        endif                                                            !srv 26.11.02
       enddo
      enddo



        if(fch_ion_neutral.ne.0.0e0_R8 ) then                            !srv
         do is = 0, ns-1                                                 !srv 24.06.08 {
!          if (.not.is_neutral(is)) then
          if (.not.is_neutral(is).and.nint(zn(is)).eq.1) then            !srv 20.10.17
           do is0 = 0, ns-1
            if (is_neutral(is0).and.nint(zn(is0)).eq.1) then             !srv 20.10.17
            conc(ix,iy,1,1) = conc(ix,iy,1,1) + csigin(ix,iy,1,1,is,is0) !*facdrift !srv 27.03.99 07.02.03 31.01.08
            conc(ix,iy,1,0) = conc(ix,iy,1,0) + 
     &         csigin(ix,iy,1,1,is,is0)*t0
            endif                                                        !srv 20.10.17
           enddo
          endif
         enddo                                                           !srv 24.06.08 }
        endif                                                            !srv
        if(fch_stochastic.ne.0.0e0_R8 ) then                             !lk 31.01.08 {
          conc(ix,iy,1,1) = conc(ix,iy,1,1) + csig_stoch(ix,iy)
          conc(ix,iy,1,0) = conc(ix,iy,1,0) + csig_stoch(ix,iy)*t0
        endif                                                            !lk 31.01.08 }

       enddo
      enddo


#endif


#ifdef WG_TODO
            if (balance.and.(balance_netcdf.ne.0)) then !djm Jan2017
              fhe_32(ix,iy,0,0) = 1.5e0_R8*fne(ix,iy,0,0)*
     &               (te(leftix(ix,iy),leftiy(ix,iy))+te(ix,iy))/2.0_R8
              fhe_52(ix,iy,0,0) = 0.0_R8
              fhe_strange(ix,iy,0,0) = chve(ix,iy,0,0)*
     &         (ne(leftix(ix,iy),leftiy(ix,iy))+ne(ix,iy))*0.5e0_R8*
     &         (te(leftix(ix,iy),leftiy(ix,iy))+te(ix,iy))/2.0_R8
              if (mdf_fhe.eq.0) then
                fhe_dia(ix,iy,0,0) = 2.5e0_R8*fledia(ix,iy,0,0)*
     &           (ne(leftix(ix,iy),leftiy(ix,iy))+ne(ix,iy))*0.5e0_R8*
     &           (te(leftix(ix,iy),leftiy(ix,iy))+te(ix,iy))/2.0_R8
                fhe_cond(ix,iy,0,0) = -cone0*tex0
                fhe_pschused(ix,iy,0,0) = 0.0_R8
              else
                fhe_dia(ix,iy,0,0) = 0.0_R8
                fhe_cond(ix,iy,0,0) = -cone0_mdf*tex0
                fhe_pschused(ix,iy,0,0) = -fhePSch(ix,iy,0,0)
              endif
              fhe_thermj(ix,iy,0,0) = -c071(ix,iy)*fch_p(ix,iy,0,0)/qe*
     &          fch_pTe*(te(leftix(ix,iy),leftiy(ix,iy))+te(ix,iy))
     &          /2.0_R8 + alfTeEh*alfe0*
     &          harm(te(leftix(ix,iy),leftiy(ix,iy)),te(ix,iy))*ehx
              fhe_ecrb(ix,iy,0,0) = 0.0_R8
            endif
          else
            fhe(ix,iy,0,0) =
     &         upwind(floe0,te(leftix(ix,iy),leftiy(ix,iy)) , te(ix,iy))
     &         -cone0*tex0 + alfTeEh*
     &         alfe0*harm(te(leftix(ix,iy),leftiy(ix,iy)),te(ix,iy))*ehx
            fhe_mdf(ix,iy,0,0) =                                         !srv 13.06.06 {
     &         upwind(floe0_mdf,
     &                        te(leftix(ix,iy),leftiy(ix,iy)),te(ix,iy))
     &         -cone0_mdf*tex0 + alfTeEh*
     &         alfe0*harm(te(leftix(ix,iy),leftiy(ix,iy)),te(ix,iy))*ehx
     &        -fhePSch(ix,iy,0,0)                                        !srv 13.06.06 }
            if (balance.and.
     &          (balance_netcdf.ne.0).and.
     &          (mdf_fhe.eq.0)) then !djm Jan2017
              if(floe0.ge.0.0_R8) then
                upfac=1.0_R8
              else
                upfac=0.0_R8
              endif
              fhe_32(ix,iy,0,0) = 1.5e0_R8*fne(ix,iy,0,0)*
     &                          (upfac*te(leftix(ix,iy),leftiy(ix,iy))+
     &                           (1.0_R8-upfac)*te(ix,iy))
              fhe_52(ix,iy,0,0) = 0.0_R8
              fhe_strange(ix,iy,0,0) = chve(ix,iy,0,0)*
     &           (ne(leftix(ix,iy),leftiy(ix,iy))+ne(ix,iy))*0.5e0_R8*
     &           (upfac*te(leftix(ix,iy),leftiy(ix,iy))+
     &            (1.0_R8-upfac)*te(ix,iy))
              fhe_dia(ix,iy,0,0) = 2.5e0_R8*fledia(ix,iy,0,0)*
     &           (ne(leftix(ix,iy),leftiy(ix,iy))+ne(ix,iy))*0.5e0_R8*
     &           (upfac*te(leftix(ix,iy),leftiy(ix,iy))+
     &            (1.0_R8-upfac)*te(ix,iy))
              fhe_thermj(ix,iy,0,0) =
     &           -c071(ix,iy)*fch_p(ix,iy,0,0)/qe*fch_pTe*
     &           (upfac*te(leftix(ix,iy),leftiy(ix,iy))+
     &            (1.0_R8-upfac)*te(ix,iy))+
     &           alfTeEh*alfe0*
     &           harm(te(leftix(ix,iy),leftiy(ix,iy)),te(ix,iy))*ehx
              fhe_ecrb(ix,iy,0,0) = 0.0_R8
              fhe_cond(ix,iy,0,0) = -cone0*tex0
              fhe_pschused(ix,iy,0,0) = 0.0_R8
            endif
            if (balance.and.
     &          (balance_netcdf.ne.0).and.
     &          (mdf_fhe.ne.0)) then !djm Jan2017
              if(floe0_mdf.ge.0.0_R8) then
                upfac=1.0_R8
              else
                upfac=0.0_R8
              endif
              fhe_32(ix,iy,0,0) = 1.5e0_R8*fne(ix,iy,0,0)*
     &                          (upfac*te(leftix(ix,iy),leftiy(ix,iy))+
     &                           (1.0_R8-upfac)*te(ix,iy))
              fhe_52(ix,iy,0,0) = 0.0_R8
              fhe_strange(ix,iy,0,0) = chve(ix,iy,0,0)*
     &           (ne(leftix(ix,iy),leftiy(ix,iy))+ne(ix,iy))*0.5e0_R8*
     &           (upfac*te(leftix(ix,iy),leftiy(ix,iy))+
     &            (1.0_R8-upfac)*te(ix,iy))
              fhe_dia(ix,iy,0,0) = 0.0_R8
              fhe_thermj(ix,iy,0,0) =
     &           -c071(ix,iy)*fch_p(ix,iy,0,0)/qe*fch_pTe*
     &           (upfac*te(leftix(ix,iy),leftiy(ix,iy))+
     &            (1.0_R8-upfac)*te(ix,iy))+
     &           alfTeEh*alfe0*
     &           harm(te(leftix(ix,iy),leftiy(ix,iy)),te(ix,iy))*ehx
              fhe_ecrb(ix,iy,0,0) = 0.0_R8
              fhe_cond(ix,iy,0,0) = -cone0_mdf*tex0
              fhe_pschused(ix,iy,0,0) = -fhePSch(ix,iy,0,0)
            endif
          endif


          if(iout_b2nph9.ne.0 .or. iout_b2npht.ne.0
     &                             .or. iout_b2wdat.eq.4) then           !srv 16.09.02 {
           if(mdf_fhe.eq.0) then                                         !srv 13.10.06 {
            w1(ix,iy)=1.5e0_R8*fne(ix,iy,0,0)
     &           *(te(leftix(ix,iy),leftiy(ix,iy))+te(ix,iy))/2.0_R8
            w2(ix,iy)=(chve(ix,iy,0,0)+2.5e0_R8*fledia(ix,iy,0,0))*      !srv 13.09.02 11.01.13 *signmf
     &      (ne(leftix(ix,iy),leftiy(ix,iy))+ne(ix,iy))*0.5e0_R8         !srv 13.09.02
     &          *(te(leftix(ix,iy),leftiy(ix,iy))+te(ix,iy))/2.0_R8
            w3(ix,iy)= -cone0*tex0
            if(alfTeEh.eq.1.0_R8)
     &       w4(ix,iy)=
     *        alfe0*harm(te(leftix(ix,iy),leftiy(ix,iy)),te(ix,iy))*ehx
            if(fch_pTe.eq.1.0_R8)
     &       w4(ix,iy)=-c071(ix,iy)*fch_p(ix,iy,0,0)/qe*                 !srv 13.01.17
     *        (te(leftix(ix,iy),leftiy(ix,iy))+te(ix,iy))/2.0_R8
            w6(ix,iy) = cone0
           else
            w1(ix,iy)=1.5e0_R8*fne(ix,iy,0,0)
     &           *(te(leftix(ix,iy),leftiy(ix,iy))+te(ix,iy))/2.0_R8
            w2(ix,iy)= chve(ix,iy,0,0)*                                  !srv 13.09.02
     &       (ne(leftix(ix,iy),leftiy(ix,iy))+ne(ix,iy))*0.5e0_R8        !srv 13.09.02
     &       *(te(leftix(ix,iy),leftiy(ix,iy))+te(ix,iy))/2.0_R8
            w3(ix,iy)= -cone0_mdf*tex0
            w4(ix,iy)=
     *       alfe0*harm(te(leftix(ix,iy),leftiy(ix,iy)),te(ix,iy))*ehx
            w6(ix,iy) = cone0_mdf
           endif
           w5(ix,iy) = 2.5e0_R8*fledia(ix,iy,0,0)*                       !srv 13.09.02 11.01.13 *signmf
     &       (ne(leftix(ix,iy),leftiy(ix,iy))+ne(ix,iy))*0.5e0_R8        !srv 13.09.02
     &          *(te(leftix(ix,iy),leftiy(ix,iy))+te(ix,iy))/2.0_R8
          endif                                                          !srv 16.09.02 }

#endif 


#ifdef WG_TODO
        if(lim_flux.eq.1 .and. conductive_limit.eq.0) then             !srv 04.02.03 17.01.12

          if(no_Ptncr_x_he.eq.0) then
           cone0 = hybr(floe0,chce(ix,iy,0,0))-abs(floe0/2.0_R8)
          else
           cone0 = chce(ix,iy,0,0)-abs(floe0/2.0_R8)
          endif

          if (chce(ix,iy,0,0).eq.0) then
           alfe0 = 0.0_R8
          else
           alfe0 = calf(ix,iy,0,0)*cone0/chce(ix,iy,0,0)
          endif

cwdk 10.03.16 unified treatment of chvemx in b2tlhe.F. Test for b2tqce_model no longer required
          flomx = chvemx(ix,iy)*(ne(leftix(ix,iy),leftiy(ix,iy))+
     &      ne(ix,iy))/2
          if(flomx.gt.0.0_R8) then
            t0 = cone0*abs(tex)/(flomx*max(
     &        te(leftix(ix,iy),leftiy(ix,iy)),te(ix,iy)))
          else
            t0=0.0_R8
          endif

          tex0 = tex/(1.d0+t0)

csw 13mar2013 SOLPS4 modified f. flux limited hybrid style (overwrite cone0)
          if(b2tfhe_hybr2 .ne. 0) then

csw 30apr2014 flo53 ONLY FOR TESTING (SOLPS4)
           if(flo53 .ne. 0) then
            floe0=floe0*5.d0/3.d0
           endif

           cone0=dim(chce(ix,iy,0,0)/(1.d0+t0), abs(floe0)*
     .      max(0.0d0,min(te(ix,iy),te(leftix(ix,iy),leftiy(ix,iy))))/
     .      max(cutlo,te(ix,iy)+te(leftix(ix,iy), leftiy(ix,iy))))
          endif


csw 27nov2012 add check for upwind (SOLPS4 style)
          if(b2_upwind .ne. 0) then
             fhe(ix,iy,0,0) =
     &        upwind(floe0,te(leftix(ix,iy),leftiy(ix,iy)), te(ix,iy))
     &         -cone0*(te(ix,iy)-te(leftix(ix,iy),leftiy(ix,iy)))
     &         + alfTeEh*
     &         alfe0*harm(te(leftix(ix,iy),leftiy(ix,iy)),te(ix,iy))*ehx
             if (balance.and.
     &           (balance_netcdf.ne.0).and.
     &           (mdf_fhe.eq.0)) then !djm Jan2017
               if(floe0.ge.0.0_R8) then
                 upfac=1.0_R8
               else
                 upfac=0.0_R8
               endif
               fhe_32(ix,iy,0,0) = 1.5e0_R8*fne(ix,iy,0,0)*
     &                           (upfac*te(leftix(ix,iy),leftiy(ix,iy))+
     &                            (1.0_R8-upfac)*te(ix,iy))
               fhe_52(ix,iy,0,0) = 0.0_R8
               fhe_strange(ix,iy,0,0) = chve(ix,iy,0,0)*
     &             (ne(leftix(ix,iy),leftiy(ix,iy))+ne(ix,iy))*0.5e0_R8*
     &             (upfac*te(leftix(ix,iy),leftiy(ix,iy))+
     &              (1.0_R8-upfac)*te(ix,iy))
               fhe_dia(ix,iy,0,0) = 2.5e0_R8*fledia(ix,iy,0,0)*
     &             (ne(leftix(ix,iy),leftiy(ix,iy))+ne(ix,iy))*0.5e0_R8*
     &             (upfac*te(leftix(ix,iy),leftiy(ix,iy))+
     &              (1.0_R8-upfac)*te(ix,iy))
               fhe_thermj(ix,iy,0,0) =
     &              -c071(ix,iy)*fch_p(ix,iy,0,0)/qe*fch_pTe*
     &              (upfac*te(leftix(ix,iy),leftiy(ix,iy))+
     &               (1.0_R8-upfac)*te(ix,iy))+
     &              alfTeEh*alfe0*
     &              harm(te(leftix(ix,iy),leftiy(ix,iy)),te(ix,iy))*ehx
               fhe_ecrb(ix,iy,0,0) = 0.0_R8
               fhe_cond(ix,iy,0,0) = -cone0*
     &                      (te(ix,iy)-te(leftix(ix,iy),leftiy(ix,iy)))
               fhe_pschused(ix,iy,0,0) = 0.0_R8
             endif
          else
             fhe(ix,iy,0,0) =
     &         floe0*(te(leftix(ix,iy),leftiy(ix,iy))+te(ix,iy))/2.0_R8
     &         -cone0*tex0 + alfTeEh*
     &         alfe0*harm(te(leftix(ix,iy),leftiy(ix,iy)),te(ix,iy))*ehx
             if (balance.and.
     &           (balance_netcdf.ne.0).and.
     &           (mdf_fhe.eq.0)) then !djm Jan2017
               fhe_32(ix,iy,0,0) = 1.5e0_R8*fne(ix,iy,0,0)*
     &           (te(leftix(ix,iy),leftiy(ix,iy))+te(ix,iy))/2.0_R8
               fhe_52(ix,iy,0,0) = 0.0_R8
               fhe_strange(ix,iy,0,0) = chve(ix,iy,0,0)*
     &           (ne(leftix(ix,iy),leftiy(ix,iy))+ne(ix,iy))*0.5e0_R8*
     &           (te(leftix(ix,iy),leftiy(ix,iy))+te(ix,iy))/2.0_R8
               fhe_dia(ix,iy,0,0) = 2.5e0_R8*fledia(ix,iy,0,0)*
     &           (ne(leftix(ix,iy),leftiy(ix,iy))+ne(ix,iy))*0.5e0_R8*
     &           (te(leftix(ix,iy),leftiy(ix,iy))+te(ix,iy))/2.0_R8
               fhe_thermj(ix,iy,0,0) = -c071(ix,iy)*fch_p(ix,iy,0,0)/qe*
     &            fch_pTe*
     &            (te(leftix(ix,iy),leftiy(ix,iy))+te(ix,iy))/2.0_R8+
     &            alfTeEh*alfe0*
     &            harm(te(leftix(ix,iy),leftiy(ix,iy)),te(ix,iy))*ehx
               fhe_ecrb(ix,iy,0,0) = 0.0_R8
               fhe_cond(ix,iy,0,0) = -cone0*tex0
               fhe_pschused(ix,iy,0,0) = 0.0_R8
             endif
          endif

          if(mdf_fhe.eq.0) then                                         !srv 17.06.09 {

csw 30apr2014
            if(b2tfhe_hybr2 .ne.0) then
              floe(ix,iy,0,0) = floe0
              cone(ix,iy,0,0,0) = cone0
            else
              floe(ix,iy,0,0) = floe0
              cone(ix,iy,0,0,0) = abs(floe0/2.0_R8)+cone0/(1.0_R8+t0)   !srv 02.01.07
            endif

          else

csw 17feb2012
            if(no_Ptncr_x_he.eq.0) then
              cone0_mdf = hybr(floe0_mdf,chce(ix,iy,0,0))-
     &                    abs(floe0_mdf/2.0_R8)
            else
              cone0_mdf = chce(ix,iy,0,0)-abs(floe0_mdf/2.0_R8)
            endif

csw 13mar2013 SOLPS4 modified f. flux limited hybrid style (overwrite cone0_mdf)
            if(b2tfhe_hybr2 .ne. 0) then
              cone0_mdf=dim(chce(ix,iy,0,0), abs(floe0_mdf)*
     .         max(0._R8,
     .         min(te(ix,iy),te(leftix(ix,iy),leftiy(ix,iy))))/
     .         max(1.e-30_R8,te(ix,iy)+te(leftix(ix,iy),leftiy(ix,iy))))
            endif

            if (chce(ix,iy,0,0).eq.0) then
              alfe0_mdf = 0.0_R8
            else
              alfe0_mdf = calf(ix,iy,0,0)*cone0_mdf/chce(ix,iy,0,0)
            endif

            t0_mdf = cone0_mdf*abs(tex)/(flomx*max(
     &         te(leftix(ix,iy),leftiy(ix,iy)),te(ix,iy)))
            tex0_mdf = tex/(1.0_R8+t0_mdf)

csw 27nov2012 added check for upwind (SOLPS4 style)
            if(b2_upwind .ne. 0) then
              fhe_mdf(ix,iy,0,0) =
     &         upwind(floe0_mdf,
     &                te(leftix(ix,iy),leftiy(ix,iy)),te(ix,iy))
     &         -cone0_mdf*tex0_mdf + alfTeEh*alfe0_mdf*
     &          harm(te(leftix(ix,iy),leftiy(ix,iy)),te(ix,iy))*ehx
              if (balance.and.
     &           (balance_netcdf.ne.0).and.
     &           (mdf_fhe.ne.0)) then !djm Jan2017
                if(floe0_mdf.ge.0.0_R8) then
                  upfac=1.0_R8
                else
                  upfac=0.0_R8
                endif
                fhe_32(ix,iy,0,0) = 1.5e0_R8*fne(ix,iy,0,0)*
     &                          (upfac*te(leftix(ix,iy),leftiy(ix,iy))+
     &                          (1.0_R8-upfac)*te(ix,iy))
                fhe_52(ix,iy,0,0) = 0.0_R8
                fhe_strange(ix,iy,0,0) = chve(ix,iy,0,0)*
     &             (ne(leftix(ix,iy),leftiy(ix,iy))+ne(ix,iy))*0.5e0_R8*
     &             (upfac*te(leftix(ix,iy),leftiy(ix,iy))+
     &              (1.0_R8-upfac)*te(ix,iy))
                fhe_dia(ix,iy,0,0) = 0.0_R8
                fhe_thermj(ix,iy,0,0) =
     &             -c071(ix,iy)*fch_p(ix,iy,0,0)/qe*fch_pTe*
     &             (upfac*te(leftix(ix,iy),leftiy(ix,iy))+
     &              (1.0_R8-upfac)*te(ix,iy))+
     &             alfTeEh*alfe0_mdf*
     &             harm(te(leftix(ix,iy),leftiy(ix,iy)),te(ix,iy))*ehx
                fhe_ecrb(ix,iy,0,0) = 0.0_R8
                fhe_cond(ix,iy,0,0) = -cone0_mdf*tex0_mdf
                fhe_pschused(ix,iy,0,0) = 0.0_R8
              endif
            else
              fhe_mdf(ix,iy,0,0) = floe0_mdf*
     &         (te(leftix(ix,iy),leftiy(ix,iy))+te(ix,iy))/2.0_R8
     &         -abs(floe0_mdf/2.0_R8)*tex
     &         -cone0_mdf*tex0_mdf + alfTeEh*alfe0_mdf*
     &          harm(te(leftix(ix,iy),leftiy(ix,iy)),te(ix,iy))*ehx
              if (balance.and.
     &           (balance_netcdf.ne.0).and.
     &           (mdf_fhe.ne.0)) then !djm Jan2017
               !djm Jan2017 share the hybrid part according to relative contributions
                if (floe0_mdf.ge.0.0_R8) then
                  hybfac=1.0_R8
                else
                  hybfac=-1.0_R8
                endif
                fhe_32(ix,iy,0,0) = 1.5e0_R8*fne(ix,iy,0,0)*
     &           (te(leftix(ix,iy),leftiy(ix,iy))+te(ix,iy))/2.0_R8-
     &            hybfac*1.5e0_R8*fne(ix,iy,0,0)/2.0_R8*tex
                fhe_52(ix,iy,0,0) = 0.0_R8
                fhe_strange(ix,iy,0,0) = chve(ix,iy,0,0)*
     &           (ne(leftix(ix,iy),leftiy(ix,iy))+ne(ix,iy))*0.5_R8*
     &           (te(leftix(ix,iy),leftiy(ix,iy))+te(ix,iy))/2.0_R8-
     &           hybfac*chve(ix,iy,0,0)*
     &           (ne(leftix(ix,iy),leftiy(ix,iy))+ne(ix,iy))*0.5_R8/
     &           2.0_R8*tex
                fhe_dia(ix,iy,0,0) = 0.0_R8
                fhe_thermj(ix,iy,0,0) =
     &           -c071(ix,iy)*fch_p(ix,iy,0,0)/qe*fch_pTe*
     &           (te(leftix(ix,iy),leftiy(ix,iy))+te(ix,iy))/2.0_R8+
     &           hybfac*c071(ix,iy)*fch_p(ix,iy,0,0)/qe*
     &           fch_pTe/2.0_R8*tex+
     &           alfTeEh*alfe0_mdf*
     &           harm(te(leftix(ix,iy),leftiy(ix,iy)),te(ix,iy))*ehx
                fhe_ecrb(ix,iy,0,0) = 0.0_R8
                fhe_cond(ix,iy,0,0) = -cone0_mdf*tex0_mdf
                fhe_pschused(ix,iy,0,0) = 0.0_R8
              endif
            endif
            floe(ix,iy,0,0) = floe0_mdf
            cone(ix,iy,0,0,0) = cone0_mdf/(1.0_R8+t0_mdf)                !srv 02.01.07
          endif                                                          !srv 17.06.09 }
          if(pot_eq.ne.1 .and. no_current.ne.0) then                     !srv 27.02.13
            fch_p = 0.0e0_R8
            ehx = 0.0_R8                                                 !srv 27.02.13
          endif                                                          !srv 30.12.09 }

        else !lim_flux 
#endif 






#ifdef WG_TODO
      if(iout_b2nph9.ne.0 .or. iout_b2npht.ne.0
     &                      .or. iout_b2wdat.eq.4) then                  !srv 16.09.02 {
       call my_out(70,nx,ny,w1,'b2tfhe__qe_32x')
!       call my_out(70,nx,ny,w2,'b2tfhe__qe_strx')
       call my_out(70,nx,ny,w3,'b2tfhe__qe_ke_gTx')
       call my_out(70,nx,ny,w4,'b2tfhe__qe_alphaTehx')
!       w5=w1+w2+w3+w4
!       call my_out(70,nx,ny,w5,'b2tfhe__qex5')
       call my_out(70,nx,ny,w5,'b2tfhe__qediax')
       call my_out(70,nx,ny,w6,'b2tfhe__chcex_Pat')

       w1=0.0e0_R8; w2=0.0e0_R8; w3=0.0e0_R8; w4=0.0e0_R8; w5=0.0e0_R8;
       w6=0.0e0_R8;

      endif                                                              !srv 16.09.02 }
#endif
#ifdef WG_TODO

        if(iout_b2nph9.ne.0 .or. iout_b2npht.ne.0
     &                        .or. iout_b2wdat.eq.4) then                !srv 16.09.02 {
         if(mdf_fhe.eq.0) then                                           !srv 13.10.06 {
          w1(ix,iy)=1.5e0_R8*fne(ix,iy,1,1)*
     1      (te(bottomix(ix,iy),bottomiy(ix,iy))+te(ix,iy))/2.0_R8
          w2(ix,iy)=2.5e0_R8*fledia(ix,iy,1,1)*
     &      (ne(bottomix(ix,iy),bottomiy(ix,iy))+ne(ix,iy))*0.5e0_R8
     1     *(te(bottomix(ix,iy),bottomiy(ix,iy))+te(ix,iy))/2.0_R8
          w3(ix,iy)=-cone0*
     1      (te(ix,iy)-te(bottomix(ix,iy),bottomiy(ix,iy)))
          w4(ix,iy)=calf(ix,iy,1,1)*
     &       harm(te(bottomix(ix,iy),bottomiy(ix,iy)),te(ix,iy))*ehy
          w6(ix,iy) = cone0
         else
          w1(ix,iy)=1.5e0_R8*fne(ix,iy,1,1)*
     1      (te(bottomix(ix,iy),bottomiy(ix,iy))+te(ix,iy))/2.0_R8
          w2(ix,iy)=0.0e0_R8
          w3(ix,iy)=-cone0_mdf*
     1      (te(ix,iy)-te(bottomix(ix,iy),bottomiy(ix,iy)))
          w4(ix,iy)=calf(ix,iy,1,1)*
     &       harm(te(bottomix(ix,iy),bottomiy(ix,iy)),te(ix,iy))*ehy
          w6(ix,iy) = cone0_mdf
         endif
         w5(ix,iy)=2.5e0_R8*fledia(ix,iy,1,1)*
     &      (ne(bottomix(ix,iy),bottomiy(ix,iy))+ne(ix,iy))*0.5e0_R8
     1     *(te(bottomix(ix,iy),bottomiy(ix,iy))+te(ix,iy))/2.0_R8
        endif                                                            !srv 16.09.02 }
#endif

#ifdef WG_TODO

         if (balance.and.(balance_netcdf.ne.0)) then !djm Jan2017
           fhe_32(ix,iy,1,1) = 1.5e0_R8*fne(ix,iy,1,1)*
     &         (te(bottomix(ix,iy),bottomiy(ix,iy))+te(ix,iy))/2.0_R8
           fhe_52(ix,iy,1,1) = 0.0_R8
           fhe_strange(ix,iy,1,1) = chve(ix,iy,1,1)*
     &         (ne(bottomix(ix,iy),bottomiy(ix,iy))+ne(ix,iy))*0.5e0_R8*
     &         (te(bottomix(ix,iy),bottomiy(ix,iy))+te(ix,iy))/2.0_R8
           if (mdf_fhe.eq.0) then
             fhe_cond(ix,iy,1,1) = -cone0*
     &                (te(ix,iy)-te(bottomix(ix,iy),bottomiy(ix,iy)))
             fhe_dia(ix,iy,1,1) = 2.5e0_R8*fledia(ix,iy,1,1)*
     &         (ne(bottomix(ix,iy),bottomiy(ix,iy))+ne(ix,iy))*0.5e0_R8*
     &         (te(bottomix(ix,iy),bottomiy(ix,iy))+te(ix,iy))/2.0_R8
          else
             fhe_cond(ix,iy,1,1) = -cone0_mdf*
     &                (te(ix,iy)-te(bottomix(ix,iy),bottomiy(ix,iy)))
             fhe_dia(ix,iy,1,1) = 0.0_R8
           endif
           fhe_pschused(ix,iy,1,1) = -fhePSch(ix,iy,1,1)
           fhe_thermj(ix,iy,1,1) = alfTeEh*calf(ix,iy,1,1)*
     &         harm(te(bottomix(ix,iy),bottomiy(ix,iy)),te(ix,iy))*ehy
           fhe_ecrb(ix,iy,1,1) = 0.0_R8
         endif
        endif
#endif


#ifdef WG_TODO

      if(iout_b2nph9.ne.0 .or. iout_b2npht.ne.0
     &                                .or. iout_b2wdat.eq.4) then   !srv 16.09.02 {
       call my_out(70,nx,ny,w1,'b2tfhe__qe_32y')        !srv 16.09.02 {
!       call my_out(70,nx,ny,w2,'b2tfhe__qe_stry')
       call my_out(70,nx,ny,w3,'b2tfhe__qe_ke_gTy')
       call my_out(70,nx,ny,w4,'b2tfhe__qe_alphaTehy')
!       w5=w1+w2+w3+w4
!       call my_out(70,nx,ny,w5,'b2tfhe__qey5')
       call my_out(70,nx,ny,w5,'b2tfhe__qediay')
       call my_out(70,nx,ny,w6,'b2tfhe__chcey_Pat')
      endif                                             !srv 16.09.02 }

      if(iout_b2nph9.ne.0 .or. iout_b2npht.ne.0
     &                                .or. iout_b2wdat.eq.4) then   !srv 16.09.02 {
       deallocate( w1, w2, w3, w4, w5, w6 )
      endif                                             !srv 16.09.02 }
#endif


      if(iout.ne.0) then                                !srv 18.06.08 {
       call my_out_us(70,nFc,1,dv%fch(1,0),'b2tfhe__fchx')
       call my_out_us(70,nFc,1,dv%fch(1,1),'b2tfhe__fchy')
       call my_out_us(70,nFc,1,dv%fchin(1,0),'b2tfhe__fchinx')
       call my_out_us(70,nFc,1,dv%fchin(1,1),'b2tfhe__fchiny')
       call my_out_us(70,nFc,1,dv%fchvispar(1,0),'b2tfhe__fchvisparx')
       call my_out_us(70,nFc,1,dv%fchvispar(1,1),'b2tfhe__fchvispary')
       call my_out_us(70,nFc,1,dv%fchvisq(1,0),'b2tfhe__fchvisqx')
       call my_out_us(70,nFc,1,dv%fchvisq(1,1),'b2tfhe__fchvisqy')
       call my_out_us(70,nFc,1,dv%fchvisper(1,0),'b2tfhe__fchvisperx')
       call my_out_us(70,nFc,1,dv%fchvisper(1,1),'b2tfhe__fchvispery')
       call my_out_us(70,nFc,1,dv%fchinert(1,0),'b2tfhe__fchinertx')
       call my_out_us(70,nFc,1,dv%fchinert(1,1),'b2tfhe__fchinerty')
       call my_out_us(70,nFc,1,dv%fchdia(1,0),'b2tfhe__fchdiax')
       call my_out_us(70,nFc,1,dv%fchdia(1,1),'b2tfhe__fchdiay')
       call my_out_us(70,nFc,1,dv%fchanml(1,0),'b2tfhe__fchanmlx')
       call my_out_us(70,nFc,1,dv%fchanml(1,1),'b2tfhe__fchanmly')
!WG_TODO       call my_out_us(70,nFc,1,fchstoch,'b2tfhe__fchstochy')        !lk 31.01.08
       call my_out_us(70,nFc,1,dv%fne(1,0),'b2tfhe__fnex')
       call my_out_us(70,nFc,1,dv%fne(1,1),'b2tfhe__fney')
       call my_out_us(70,nFc,1,dv%conc(1,0),'b2tfhe__concx')
       call my_out_us(70,nFc,1,dv%conc(1,1),'b2tfhe__concy')
       call my_out_us(70,nFc,1,dv%floe(1,0),'b2tfhe__floex')
       call my_out_us(70,nFc,1,dv%floe(1,1),'b2tfhe__floey')
       call my_out_us(70,nFc,1,dv%cone(1,0,0),'b2tfhe__conex0')
       call my_out_us(70,nFc,1,dv%cone(1,0,1),'b2tfhe__conex1')
       call my_out_us(70,nFc,1,dv%cone(1,0,2),'b2tfhe__conex2')
       call my_out_us(70,nFc,1,dv%cone(1,1,0),'b2tfhe__coney0')
       call my_out_us(70,nFc,1,dv%cone(1,1,1),'b2tfhe__coney1')
       call my_out_us(70,nFc,1,dv%cone(1,1,2),'b2tfhe__coney2')
!WG_TODO       call my_out_us(70,nFc,1,avisper(1),'b2tfhe__avisperx')
!WG_TODO       call my_out_us(70,nFc,1,avisper(1),'b2tfhe__avispery')
!WG_TODO      call my_out_us(70,nFc,1,bvisper(1),'b2tfhe__bvisperx')
!WG_TODO       call my_out_us(70,nFc,1,bvisper(1),'b2tfhe__bvispery')
!WG_TODO       call my_out_us(70,nFc,1,cvisper(1),'b2tfhe__cvisperx')
!WG_TODO       call my_out_us(70,nFc,1,cvisper(1),'b2tfhe__cvispery')
       call my_out_us(70,nFc,1,dv%fch_p(1,0),'b2tfhe__fch_px')
       call my_out_us(70,nFc,1,dv%fch_p(1,1),'b2tfhe__fch_py')
       call my_out_us(70,nFc,1,co%csig_an(1,0),'b2tfhe__csig_anx')
       call my_out_us(70,nFc,1,co%csig_an(1,1),'b2tfhe__csig_any')
       call my_out_us(70,nFc,1,co%chve(1,0),'b2tfhe__chvex')
       call my_out_us(70,nFc,1,co%chve(1,1),'b2tfhe__chvey')
       call my_out_us(70,nFc,1,co%csig(1,0),'b2tfhe__csigx')
       call my_out_us(70,nFc,1,co%csig(1,1),'b2tfhe__csigy')
       call my_out_us(70,nFc,1,co%calf(1,0),'b2tfhe__calfx')
       call my_out_us(70,nFc,1,co%calf(1,1),'b2tfhe__calfy')
       call my_out_us(70,nFc,1,co%chci(1,0),'b2tfhe__chcix')
       call my_out_us(70,nFc,1,co%chci(1,1),'b2tfhe__chciy')
       call my_out_us(70,nFc,1,co%chce(1,0),'b2tfhe__chcex')
       call my_out_us(70,nFc,1,co%chce(1,1),'b2tfhe__chcey')
       call my_out_us(70,nCv,0,pl%te,'b2tfhe__te')
       call my_out_us(70,nCv,0,dv%ne,'b2tfhe__ne')
       call my_out_us(70,nCv,0,pl%po,'b2tfhe__po')
       do is = 0, ns-1
        write (charns,'(i3.3)') is
        call my_out_us(70,nCv,0,pl%na(1,is),'b2tfhe__na'//charns)
        do is0=0,ns-1
         write (charns0,'(i3.3)') is0
         call my_out_us(70,nFc,1,co%csigin(1,0,is,is0),'b2tfhe__csiginx'
     &        //charns//'_'//charns0)
         call my_out_us(70,nFc,1,co%csigin(1,1,is,is0),'b2tfhe__csiginy'
     &        //charns//'_'//charns0)
        enddo
       enddo
       call my_out_us(70,nFc,1,dv%fhe(1,0),'b2tfhe__fhex')
       call my_out_us(70,nFc,1,dv%fhe(1,1),'b2tfhe__fhey')
       call my_out_us(70,nFc,1,dv%fhe_mdf(1,0),'b2tfhe__fhe_mdfx')
       call my_out_us(70,nFc,1,dv%fhe_mdf(1,1),'b2tfhe__fhe_mdfy')
       call my_out_us(70,nFc,1,dv%fhePSch(1,0),'b2tfhe__fhePSchx')
       call my_out_us(70,nFc,1,dv%fhePSch(1,1),'b2tfhe__fhePSchy')
       call my_out_us(70,nFc,1,ehx,'b2tfhe__ehx_')            !srv 13.01.17
       call my_out_us(70,nFc,1,ehy,'b2tfhe__ehy_')            !srv 13.01.17
      endif                                                          !srv 18.06.08 }
#ifdef WG_TODO
      if(iout_b2nph9.ne.0 .or. iout_b2npht.ne.0
     &                          .or. iout_b2wdat.eq.4) then          !srv 11.09.09 {
       call my_out(70,nx,ny,fhe(-1,-1,0,0),'b2tfhe__fhe_no_mdfx')
       call my_out(70,nx,ny,fhe(-1,-1,1,1),'b2tfhe__fhe_no_mdfy')
!       call my_out(70,nx,ny,fhe_mdf(-1,-1,0,0),'b2tfhe__fhe_mdfx')
!       call my_out(70,nx,ny,fhe_mdf(-1,-1,1,1),'b2tfhe__fhe_mdfy')
       call my_out(70,nx,ny,fhePSch(-1,-1,0,0),'b2tfhe__fhePSchx')
       call my_out(70,nx,ny,fhePSch(-1,-1,1,1),'b2tfhe__fhePSchy')
      endif

      if(iout_b2npp7.ne.0 .or. iout_b2nppo.ne.0
     &                          .or. iout_b2wdat.eq.4) then
!       call my_out(70,nx,ny,fch(-1,-1,0,0),'b2tfhe__fchx')
!       call my_out(70,nx,ny,fch(-1,-1,1,1),'b2tfhe__fchy')
       call my_out(70,nx,ny,fchin(-1,-1,0,0),'b2tfhe__fchinx')
       call my_out(70,nx,ny,fchin(-1,-1,1,1),'b2tfhe__fchiny')
       call my_out(70,nx,ny,fchvispar(-1,-1,0,0),'b2tfhe__fchvisparx')
       call my_out(70,nx,ny,fchvispar(-1,-1,1,1),'b2tfhe__fchvispary')
!       call my_out(70,nx,ny,fchvisq(-1,-1,0,0),'b2tfhe__fchvisqx')
       call my_out(70,nx,ny,fchvisq(-1,-1,1,1),'b2tfhe__fchvisqy')
!       call my_out(70,nx,ny,fchvisper(-1,-1,0,0),'b2tfhe__fchvisperx')
       call my_out(70,nx,ny,fchvisper(-1,-1,1,1),'b2tfhe__fchvispery')
       call my_out(70,nx,ny,fchinert(-1,-1,0,0),'b2tfhe__fchinertx')
       call my_out(70,nx,ny,fchinert(-1,-1,1,1),'b2tfhe__fchinerty')
       call my_out(70,nx,ny,fchdia(-1,-1,0,0),'b2tfhe__fchdiax')
       call my_out(70,nx,ny,fchdia(-1,-1,1,1),'b2tfhe__fchdiay')
       call my_out(70,nx,ny,fchanml(-1,-1,0,0),'b2tfhe__fchanmlx')
       call my_out(70,nx,ny,fchanml(-1,-1,1,1),'b2tfhe__fchanmly')
       call my_out(70,nx,ny,fch_p(-1,-1,0,0),'b2tfhe__fch_px')
!      call my_out(70,nx,ny,fch_p(-1,-1,1,1),'b2tfhe__fch_py')
       call my_out(70,nx,ny,csig_an(-1,-1,0,0),'b2tfhe__csig_anx')
       call my_out(70,nx,ny,csig_an(-1,-1,1,1),'b2tfhe__csig_any')
       call my_out(70,nx,ny,fchstoch,'b2tfhe__fchstochy')            !lk 31.01.08
      endif                                                          !srv 11.09.09 }
#endif
* ..return
      ncall = ncall+1
      call subend ()
      return

*-----------------------------------------------------------------------
*.end b2tfhe_

      end
