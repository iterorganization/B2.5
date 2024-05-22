*-----------------------------------------------------------------------
*.specification

      subroutine b2sihs_ (nx, ny, ns, BoRiS, ismain, fac_ExB,
     *  hx, hy, qc, qs, vol, gs, bb, pbs, lnlam, na, ua, ue,          !srv 07.06.00 02.06.99 07.05.08 05.07.17
     *  te, ti, po, ne, ni, fne, fch, fchstoch, cvsa, cvsa_cl,        !som 16.08.21 !srv 24.08.21
     *  cvsa_drho,                                                    !som 16.08.21
     *  cvsa_uAdp_albe, cvsa_uBdp_al,                                 !som 08.05.24
     *  she0, shi0)                                                   !srv 28.06.03 07.07.09
      use b2mod_types
      use b2mod_rates
      use b2mod_sources
      use b2mod_tallies
      use b2mod_indirect
      use b2mod_geo_corner
      use b2mod_anomalous_transport
      use b2mod_constants
      use b2mod_b2cmfs
      use b2mod_b2cmpa
      use b2mod_b2zhco
     & , only: zhd_corr, zh_vis_type, zhcscorr                        !som 21.03.21  21.12.23
      use b2mod_zhfrtf
     & , only: avm_u, nal,                                            !som 21.03.21 16.08.21
     &         b2mod_zhfrtf_avm_u, b2mod_zhfrtf_prep                  !som 18.08.21
      use b2mod_b2mrel
     & , only : s2nucl, nucl2s, znnucl_int
      use b2mod_geo
     & , only : hy1
      use b2mod_balance !djm Jan2017
     & , only : b2sihs_divue0to3, b2sihs_exbe0to3,
     &          b2sihs_joule0to3, b2sihs_divua0to3,
     &          b2sihs_exba0to3, b2sihs_visa0to3, b2sihs_fraa0to3,
     &          balance_netcdf
      use b2mod_transport_fun                                         !srv 05.07.17
      use b2mod_external                                              !srv 05.07.17
      use b2mod_subsys
      implicit none
*   ..input arguments (unchanged on exit)
      integer nx, ny, ns, ismain                                      !srv 05.07.17
      real (kind=R8) :: BoRiS,
     *  fac_ExB(-1:nx,-1:ny), hx(-1:nx,-1:ny), hy(-1:nx,-1:ny),       !srv 07.06.00 02.06.99
     *  vol(-1:nx,-1:ny), gs(-1:nx,-1:ny,0:2),
     &  na(-1:nx,-1:ny,0:ns-1),
     &  ua(-1:nx,-1:ny,0:ns-1), ue(-1:nx,-1:ny),                      !srv 05.03.11
     &  te(-1:nx,-1:ny), ti(-1:nx,-1:ny),
     &  po(-1:nx,-1:ny), ne(-1:nx,-1:ny), ni(-1:nx,-1:ny,0:1),
     &  qc(-1:nx,-1:ny,0:1), qs(-1:nx,-1:ny,0:1),
     &  lnlam(-1:nx,-1:ny)                                            !srv 07.05.08 01.07.09
      real (kind=R8) ::
     *  fne(-1:nx,-1:ny,0:1,0:1), fchstoch(-1:nx,-1:ny,0:1,0:1),      !srv 24.08.21
     &  fch(-1:nx,-1:ny,0:1,0:1), cvsa(-1:nx,-1:ny,0:1,0:1,0:ns-1),
     &  cvsa_cl(-1:nx,-1:ny,0:1,0:1,0:ns-1),
     &  cvsa_drho(-1:nx,-1:ny,0:1,0:1,0:ns-1),                        !som 16.08.21
     &  cvsa_uAdp_albe(-1:nx,-1:ny,0:1,0:1,0:nspecies-1,0:nspecies-1), !som 08.05.24
     &  cvsa_uBdp_al(-1:nx,-1:ny,0:1,0:1,0:nspecies-1),               !som 08.05.24
     2  bb(-1:nx,-1:ny,0:3), pbs(-1:nx,-1:ny,0:1)                     !srv 28.06.03
*   ..output arguments (unspecified on entry)
      real (kind=R8) ::
     *   she0(-1:nx,-1:ny,0:3), shi0(-1:nx,-1:ny,0:3)
*-----------------------------------------------------------------------
*.documentation
*
*  1. purpose
*
*     B2SIHS computes the source coefficients for a set of terms in the
*     heat equations. These terms specify the contributions of the
*     divergence of the flow field, of viscous heating, and of the
*     electron-atom and atom-atom friction.
*
*-----------------------------------------------------------------------
*.declarations

*   ..global named constants

*   ..local variables
      integer ncall, style, n2, ix, iy, is, isb, k, inucl, iz
      integer, save :: iout = 0, iout_b2nph9 = 0, iout_b2npht = 0,    !srv 11.09.09
     &                 iout_b2wdat = 0
      integer, save :: icase = 0
      integer, save :: istyle_Joule_heating = 1                       !srv 17.01.14
      integer, save :: b2sigp_style = 2                               !srv 05.07.17
      integer, save :: use_eirene = 0
      integer, save :: style_int_vel = 0                              !srv 01.04.19
      integer, save :: impr_form_tf, impr_form_fr                     !som 13.11.20
      integer, save :: Zhdanov_closure = 0, Zhdanov_test = 0,         !som 21.03.21
     &                 Zhdanov_tf_fr = 1                              !som 01.07.21
      integer, save :: Zhdanov_visu = 1                               !som 16.08.21
      real (kind=R8) :: mm(0:ns-1), am1(0:ns-1)                       !som 13.11.20
      real (kind=R8), save :: mass_multiplicator                      !som 13.11.20
      real (kind=R8) :: ft                                            !som 16.08.21
      real (kind=R8) ::
     *  vti32(-1:nx,-1:ny)
      real (kind=R8) ::
     *  rf0, rf1, rf2, rf3, rf4, phm0, phm1, phm2, phm3,
     &  phm4, phm5, phm6, wrk0(-1:nx,-1:ny),                          !srv 02.06.99
     &  phm7, phm8, wrk1(-1:nx,-1:ny),                                !srv 28.08.02
     &  dpo(-1:nx,-1:ny,0:1,0:1), dOnedBsq(-1:nx,-1:ny,0:1,0:1),      !srv 03.04.07
     &  wght(-1:nx,-1:ny)
      real (kind=R8), save :: fch_stochastic = 0.0e0_R8               !lk 31.01.08 {
      real (kind=R8) ::
     *  shidu(-1:nx,-1:ny,0:3), shidd(-1:nx,-1:ny,0:3),               !srv 11.09.09 {
     *  shiva(-1:nx,-1:ny,0:3), shivc(-1:nx,-1:ny,0:3),
     *  shifr(-1:nx,-1:ny,0:3), shedu(-1:nx,-1:ny,0:3),
     *  shedd(-1:nx,-1:ny,0:3), shefr(-1:nx,-1:ny,0:3),
     *  shest(-1:nx,-1:ny,0:3), duar(-1:nx,-1:ny,0:1),                !srv 17.12.13
     *  wrks(-1:nx,-1:ny), ehx(-1:nx,-1:ny,0:1), ehy(-1:nx,-1:ny,0:1) !srv 13.01.17
      real (kind=R8) ::
     3  zeff(-1:nx,-1:ny), zmff(-1:nx,-1:ny), ce1(-1:nx,-1:ny),       !srv 13.01.17 05.07.17 {
     4  ce2(-1:nx,-1:ny), cimp1(-1:nx,-1:ny,0:ns-1),
     *  cimp2(-1:nx,-1:ny,0:ns-1),                                    !som 13.11.20
     5  zetap(-1:nx,-1:ny), zetae(-1:nx,-1:ny),
     6  gti(-1:nx,-1:ny), gte(-1:nx,-1:ny), ne2(-1:nx,-1:ny),
     7  smbch(-1:nx,-1:ny,0:3), smbtf(-1:nx,-1:ny,0:3),
     *  smbfria(-1:nx,-1:ny,0:3), smbfrea(-1:nx,-1:ny,0:3),
     *  smbtfia(-1:nx,-1:ny,0:3), smbtfea(-1:nx,-1:ny,0:3),           !srv 13.01.17 05.07.17 }
     8  zmc(-1:nx,-1:ny), z11(-1:nx,-1:ny),                           !som 13.11.20
     9  z12(-1:nx,-1:ny), z22(-1:nx,-1:ny),
     1  zf1(-1:nx,-1:ny), zf2(-1:nx,-1:ny)
      character chk*1                                                 !srv 11.09.09 } !srv 05.07.17
      save ncall, style, rf0, rf1, rf2, rf3, rf4, phm0, phm1, phm2,
     &  phm3, phm4, phm5, phm6, phm7, phm8                            !srv 02.06.99 27.12.01 28.08.02
      real (kind=R8) ::
     *  t0, t1, ctaup, coef, czetap, czetae                           !srv 05.07.17
*   ..procedures
      intrinsic min, max, sqrt
      external xertst, ipgeti, ipgetr, sfill
      external b2xvff, b2xvfx, b2xvfv, b2xvsg
*   ..initialization
      data ncall/0/, style/0/
      data rf0/1.0e0_R8/, rf1/1.0e0_R8/, rf2/1.0e0_R8/, rf3/1.0e0_R8/,
     &   rf4/1.0e0_R8/, phm0/1.0e0_R8/, phm1/1.0e0_R8/, phm2/1.0e0_R8/,
     &  phm3/1.0e0_R8/, phm4/1.0e0_R8/, phm5/1.0e0_R8/, phm6/1.0e0_R8/,!srv 02.06.99
     &  phm7/0.0e0_R8/, phm8/1.0e0_R8/                                 !srv 27.12.01 28.08.02
*     Note 0: style determines the form of the strange electron-atom
*     energy transfer term.
*     Note 1: rf0-5 are numerical parameters, see the code. All must
*     be non-negative and of order unity; larger values are meant to
*     provide more stabilization. Some or all may be set to zero and
*     removed eventually. I have had some experience to indicate that
*     rf0, rf2, rf4 may have a useful function; however, later changes
*     in the code may have made all these parameters superfluous.

*-----------------------------------------------------------------------
*.computation

* ..preliminaries
*   ..subprogram start-up calls
      call subini ('b2sihs_')
*   ..set internal parameters on first call
      if (ncall.eq.0) then
       call ipgeti ('b2sihs_style', style)
       call ipgetr ('b2sihs_rf0', rf0)
       call ipgetr ('b2sihs_rf1', rf1)
       call ipgetr ('b2sihs_rf2', rf2)
       call ipgetr ('b2sihs_rf3', rf3)
       call ipgetr ('b2sihs_rf4', rf4)
       call ipgetr ('b2sihs_phm0', phm0)
       call ipgetr ('b2sihs_phm1', phm1)
       call ipgetr ('b2sihs_phm2', phm2)
       call ipgetr ('b2sihs_phm3', phm3)
       call ipgetr ('b2sihs_phm4', phm4)
       call ipgetr ('b2sihs_phm5', phm5)                               !srv 02.06.99
       call ipgetr ('b2sihs_phm6', phm6)                               !srv 02.06.99
       call ipgetr ('b2sihs_phm7', phm7)                               !srv 27.12.01
       call ipgetr ('b2sihs_phm8', phm8)                               !srv 28.08.02
       call xertst (0.le.style, 'faulty internal parameter style')
       call xertst (0.le.min(rf0,rf1,rf2,rf3,rf4),
     &   'faulty internal parameter rf0/1/2/3/4')
       call xertst (0.le.min(phm0,phm1,phm2,phm3,phm4,phm5,phm6),      !srv 02.06.99
     &   'faulty internal parameter phm0/1/2/3/4/5/6')                 !srv 02.06.99
       call ipgeti ('b2sihs_iout', iout)                               !srv 10.09.03
       call ipgeti ('b2nph9_iout', iout_b2nph9)                        !srv 11.09.09
       call ipgeti ('b2npht_iout', iout_b2npht)                        !srv 11.09.09
       call ipgeti ('b2wdat_iout', iout_b2wdat)                        !iys 17.04.17
       call ipgeti ('b2tlnl_ii', icase)
       call xertst (0.eq.icase .or. 1.eq.icase,
     &               'faulty input b2tlnl_ii')
       call ipgetr ('b2tfhe_stochastic', fch_stochastic)               !lk 31.01.08
       call ipgeti ('b2sihs_istyle_Joule_heating',
     &                      istyle_Joule_heating)                      !srv 17.01.14
       call ipgeti ('b2sigp_style', b2sigp_style)                      !srv 05.07.17
       call ipgeti ('b2mndr_eirene', use_eirene)
       call ipgeti ('b2sihs_style_int_vel', style_int_vel)             !srv 01.04.19
       impr_form_tf = 0                                                !som 13.11.20
       impr_form_fr = 0
       mass_multiplicator = 1.0_R8
       call ipgeti ('b2npmo_impr_form_tf', impr_form_tf)
       call ipgeti ('b2npmo_impr_form_fr', impr_form_fr)
       call ipgetr ('b2npmo_mass_multiplicator', mass_multiplicator)
       call ipgeti ('b2tral_Zhdanov_closure', Zhdanov_closure)         !som 21.03.21
       call ipgeti ('b2sifr_Zhdanov_test', Zhdanov_test)               !som 21.03.21
       call ipgeti ('b2sifr_Zhdanov_tf_fr', Zhdanov_tf_fr)             !som 01.07.21
       call ipgeti ('b2tral_Zhdanov_visu', Zhdanov_visu)               !som 18.08.21
      endif
*   ..set n2
      n2 = (nx+2)*(ny+2)
*   ..test nx, ny, ns
      call xertst (0.le.nx.and.0.le.ny, 'faulty argument nx, ny')
      call xertst (1.le.ns, 'faulty argument ns')
*   ..extensive tests on first few calls
      if (ncall.lt.3) then
*    ..test sign of vol
       call b2xvsg (n2, vol, 1, 'vol', '.gt.')
!*    ..test edge values of ua                                     !srv 15.04.05 }
!       do is = 0, ns-1
!        call b2xvfx (nx, ny, ua(-1,-1,is), 'ua')
!       enddo                                                      !srv 15.04.05 }
*    ..test sign of na, te, ti, ne, ni
       call b2xvsg (n2*ns, na, 1, 'na', '.gt.')
       call b2xvsg (n2, te, 1, 'te', '.gt.')
       call b2xvsg (n2, ti, 1, 'ti', '.gt.')
       call b2xvsg (n2, ne, 1, 'ne', '.gt.')
       call b2xvsg (n2, ni, 1, 'ni', '.gt.')
*    ..test edge values of fna, fne, fch
!       do is = 0, ns-1                                            !srv 28.06.03 {
!        call b2xvff (nx, ny, fna(-1,-1,0,0,is), 'fna')
!       enddo                                                      !srv 28.06.03 }
       call b2xvff (nx, ny, fne, 'fne')
       call b2xvff (nx, ny, fch, 'fch')
*    ..test edge values of cvsa
       do is = 0, ns-1
        call b2xvff (nx, ny, cvsa(-1,-1,0,0,is), 'cvsa')           !srv 01.07.05
       enddo
      endif
*   ..compute ctaup
*     Following Balescu; Braginskii has ctaup different by a factor
*     sqrt(2). The mathematical expression is:
*       t0 = (3/(4*sqrt(2*pi)))*(4*pi*eps0/qe**2)**2*ev*sqrt(ev)
*       ctaup = t0*sqrt(mp)/lnlam
*     The numerical value in SI units is: ctaup = 1.4743d13/lnlam.
*     The following calculation of t0 is designed to avoid overflow
*     even for IEEE single precision arithmetic.
      t1 = (4*pi*eps0/qe)/qe
      t0 = (3/(4*sqrt(2*pi)))*t1*(t1*ev*sqrt(ev))
      ctaup = t0*sqrt(mp) ! /lnlam it is taken below     !srv 01.07.09
*   ..compute vti32
      do iy = -1, ny
       do ix = -1, nx
        vti32(ix,iy) = vol(ix,iy)/((ti(ix,iy)/ev)*sqrt(ti(ix,iy)/ev))
       enddo
      enddo

* ..compute heat source terms
*   ..initialize sources to zero
      call sfill (n2*4, 0.0e0_R8, she0, 1)
      call sfill (n2*4, 0.0e0_R8, shi0, 1)
      call sfill (n2*4, 0.0e0_R8, shedu, 1)                         !srv 11.09.09 {
      call sfill (n2*4, 0.0e0_R8, shedd, 1)
      call sfill (n2*4, 0.0e0_R8, shefr, 1)
      call sfill (n2*4, 0.0e0_R8, shidu, 1)
      call sfill (n2*4, 0.0e0_R8, shidd, 1)
      call sfill (n2*4, 0.0e0_R8, shifr, 1)
      call sfill (n2*4, 0.0e0_R8, shiva, 1)
      call sfill (n2*4, 0.0e0_R8, shivc, 1)                         !srv 11.09.09 }
      call sfill (n2*4, 0.0e0_R8, shest, 1)                         !srv 17.12.13
*   ..compute the Coulomb logarithm
      call b2tlnl (nx, ny, icase, te, ti, ne, lnlam)                !srv 01.07.09  20.09.11
*   ..contribution from stochastic electron flow
      if(fch_stochastic.ne.0.0e0_R8) then                           !lk 31.01.08 { !srv 17.12.13
       do iy = 0, ny-1
        do ix = 0, nx-1
         if(leftix(ix,iy).ne.-2 .and. rightix(ix,iy).ne.nx+1) then
          shest(ix,iy,0) = shest(ix,iy,0)-
     &    (ne(topix(ix,iy),topiy(ix,iy))*te(topix(ix,iy),topiy(ix,iy))
     &    -ne(bottomix(ix,iy),bottomiy(ix,iy))
     &    *te(bottomix(ix,iy),bottomiy(ix,iy)))
     &    /hy(ix,iy)/qe/ne(ix,iy)*vol(ix,iy)
     &    *(fchstoch(ix,iy,1,1)*(hy1(ix,iy)                         !srv 24.08.21
     &    +hy1(bottomix(ix,iy),bottomiy(ix,iy)))
     &    /(vol(ix,iy) + vol(bottomix(ix,iy),bottomiy(ix,iy)))
     &    +fchstoch(topix(ix,iy),topiy(ix,iy),1,1)*(hy1(ix,iy)      !srv 24.08.21
     &    +hy1(topix(ix,iy),topiy(ix,iy)))
     &    /(vol(ix,iy) + vol(topix(ix,iy),topiy(ix,iy))))/4.0e0_R8-
     &    (ne(rightix(ix,iy),rightiy(ix,iy))
     &    *te(rightix(ix,iy),rightiy(ix,iy))
     &    -ne(leftix(ix,iy),leftiy(ix,iy))
     &    *te(leftix(ix,iy),leftiy(ix,iy)))
     &    /hy(ix,iy)/qe/ne(ix,iy)*vol(ix,iy)
     &    *(fchstoch(ix,iy,0,1)*(hy1(ix,iy)                         !srv 24.08.21
     &    +hy1(leftix(ix,iy),leftiy(ix,iy)))
     &    *(hy(ix,iy)+hy(leftix(ix,iy),leftiy(ix,iy)))/
     &     (hx(ix,iy)+hx(leftix(ix,iy),leftiy(ix,iy)))*
     &      qc(ix,iy,0)
     &    /(vol(ix,iy) + vol(leftix(ix,iy),leftiy(ix,iy)))
     &    +fchstoch(rightix(ix,iy),rightiy(ix,iy),0,1)*(hy1(ix,iy)  !srv 24.08.21
     &    +hy1(rightix(ix,iy),rightiy(ix,iy)))
     &    *(hy(ix,iy)+hy(rightix(ix,iy),rightiy(ix,iy)))/
     &     (hx(ix,iy)+hx(rightix(ix,iy),rightiy(ix,iy)))*
     &      qc(rightix(ix,iy),rightiy(ix,iy),0)
     &    /(vol(ix,iy) + vol(rightix(ix,iy),rightiy(ix,iy))))/4.0e0_R8
         endif
        enddo
       enddo
       call b2saxpy (n2*4, 1.0e0_R8, shest, 1, she0, 1)             !srv 17.12.13 24.08.21
      endif                                                         !lk 31.01.08 }
*   ..contribution from divergence ue
      do iy = 0, ny-1                                               !srv 05.03.11 {
        do ix = 0, nx
          if(leftix(ix,iy).eq.-2) then
            wrk0(ix,iy) = 0.0e0_R8
          else
            if(style_int_vel.eq.0) then                             !srv 01.04.19 {
              wrk0(ix,iy)=(hx(ix,iy)*ue(leftix(ix,iy),leftiy(ix,iy))+
     1                   hx(leftix(ix,iy),leftiy(ix,iy))*ue(ix,iy))/
     2        (hx(leftix(ix,iy),leftiy(ix,iy))+hx(ix,iy))*pbs(ix,iy,0)
            else
              wrk0(ix,iy)=(vol(ix,iy)*ue(leftix(ix,iy),leftiy(ix,iy))+
     1                   vol(leftix(ix,iy),leftiy(ix,iy))*ue(ix,iy))/
     2        (vol(leftix(ix,iy),leftiy(ix,iy))+vol(ix,iy))*pbs(ix,iy,0)
            endif                                                     !srv 01.04.19 }
          endif
        enddo
      enddo                                                           !srv 05.03.11 }
      b2divue=0.0_R8
      b2sihs_divue=0.0_R8
      do iy = 0, ny-1
       do ix = 0, nx-1
        if(leftix(ix,iy).ne.-2 .and. rightix(ix,iy).ne.nx+1) then     !srv 26.11.02
         t0 = -phm0*(1.0_R8-BoRiS)*
     &   (wrk0(rightix(ix,iy),rightiy(ix,iy))-wrk0(ix,iy))            !srv 31.05.99 27.06.03 05.03.11
!     &              +(uesy(topix(ix,iy),topiy(ix,iy))-uesy(ix,iy)))  !srv 31.05.99 27.06.03
         shedu(ix,iy,0) = shedu(ix,iy,0)+
     &    (max(0.0e0_R8,t0)+rf0*abs(t0))*ne(ix,iy)*te(ix,iy)
         shedu(ix,iy,3) = shedu(ix,iy,3)+min(0.0e0_R8,t0)-rf0*abs(t0)
         b2divue(region(ix,iy,0))=b2divue(region(ix,iy,0))+
     1     t0*ne(ix,iy)*te(ix,iy)
         b2sihs_divue(ix,iy)=t0*ne(ix,iy)*te(ix,iy)
        endif                                                         !srv 26.11.02
       enddo
      enddo
      if (balance_netcdf.ne.0) then !djm Jan2017
        b2sihs_divue0to3 = shedu
      endif
      call b2saxpy (n2*4, 1.0e0_R8, shedu, 1, she0, 1)                !srv 11.09.09
*   ..contribution from divergence ua
      b2divua=0.0_R8
      b2sihs_divua=0.0_R8
      do is = 0, ns-1
       do iy = 0, ny-1                                                !srv 21.07.05 {
        do ix = 0, nx
          if(leftix(ix,iy).eq.-2) then
            wrk0(ix,iy) = 0.0e0_R8
          else
            if(style_int_vel.eq.0) then                               !srv 01.04.19 {
             wrk0(ix,iy)=(hx(ix,iy)*ua(leftix(ix,iy),leftiy(ix,iy),is)+
     1                hx(leftix(ix,iy),leftiy(ix,iy))*ua(ix,iy,is))/
     2        (hx(leftix(ix,iy),leftiy(ix,iy))+hx(ix,iy))*pbs(ix,iy,0)
            else
             wrk0(ix,iy)=(vol(ix,iy)*ua(leftix(ix,iy),leftiy(ix,iy),is)+
     1                vol(leftix(ix,iy),leftiy(ix,iy))*ua(ix,iy,is))/
     2       (vol(leftix(ix,iy),leftiy(ix,iy))+vol(ix,iy))*pbs(ix,iy,0)
            endif                                                     !srv 01.04.19 }
          endif                                                       !srv 21.07.05 }
!         if(leftix(ix,iy).eq.-2) then
!           wrk0(ix,iy) = 0.0e0_R8
!         else
!           wrk0(ix,iy) = (vol(ix,iy)*bb(ix,iy,0)/bb(ix,iy,3)*ua(ix,iy,is)-
!     &                   vol(leftix(ix,iy),leftiy(ix,iy))*
!     &                    bb(leftix(ix,iy),leftiy(ix,iy),0)/
!     &                    bb(leftix(ix,iy),leftiy(ix,iy),3)*
!     &                    ua(leftix(ix,iy),leftiy(ix,iy),is))/
!     &             (hx(ix,iy)+hx(leftix(ix,iy),leftiy(ix,iy)))*2.0e0_R8
!         endif
        enddo
       enddo
       do iy = 0, ny-1
        do ix = 0, nx-1
         if(leftix(ix,iy).ne.-2 .and. rightix(ix,iy).ne.nx+1) then    !srv 26.11.02
          t0 = -phm1*(na(ix,iy,is)/ni(ix,iy,0))*(1.0_R8-BoRiS)*       !srv 21.07.05 {
     &            ( wrk0(rightix(ix,iy),rightiy(ix,iy))-wrk0(ix,iy))  !srv 21.07.05 }
!          t0 = -phm1*(na(ix,iy,is)/ni(ix,iy,0))*(1.0_R8-BoRiS)*(     !srv 31.05.99
!     &        pbs(rightix(ix,iy),rightiy(ix,iy),0)*                  !srv 28.06.03 {
!     &        ua(rightix(ix,iy),rightiy(ix,iy),is)-
!     &        pbs(ix,iy,0)*ua(ix,iy,is)
!     &                                          )                    !srv 28.06.03 }
!     &      (uasx(rightix(ix,iy),rightiy(ix,iy),is)-uasx(ix,iy,is))) !srv 31.05.99
!     &      (uasx(rightix(ix,iy),rightiy(ix,iy),is)-uasx(ix,iy,is))  !srv 31.05.99 27.06.03
!     &     +(uasy(topix(ix,iy),topiy(ix,iy),is)-uasy(ix,iy,is)))     !srv 31.05.99 27.06.03
          shidu(ix,iy,0) = shidu(ix,iy,0)+
     &      (max(0.0e0_R8,t0)+rf1*abs(t0))*ni(ix,iy,0)*ti(ix,iy)
          shidu(ix,iy,3) = shidu(ix,iy,3)+min(0.0e0_R8,t0)-rf1*abs(t0)
          b2divua(region(ix,iy,0))=b2divua(region(ix,iy,0))
     1         +t0*ni(ix,iy,0)*ti(ix,iy)
          b2sihs_divua(ix,iy)=b2sihs_divua(ix,iy)
     1         +t0*ni(ix,iy,0)*ti(ix,iy)
         endif                                                        !srv 26.11.02
        enddo
       enddo
      enddo
      if (balance_netcdf.ne.0) then !djm Jan2017
        b2sihs_divua0to3 = shidu
      endif
      call b2saxpy (n2*4, 1.0e0_R8, shidu, 1, shi0, 1)                !srv 11.09.09
*   ..contribution from divergence of the electrical drift            !srv 02.06.99
      wght = 1.0_R8
      call grad(nx, ny, hx, hy, wght, qs, po, dpo)
      call grad(nx, ny, hx, hy, wght, qs, OnedBsq, dOnedBsq)
!      call grad_p(nx, ny, 1, hx, wght, OnedBsq, OnedBsqc,
!     .  dOnedBsq(-1,-1,0,0))
!      call grad_r(nx, ny, 1, hx, hy, wght, qs, OnedBsq, OnedBsqc,
!     .  dOnedBsq(-1,-1,0,1))

      b2exbe=0.0_R8
      b2sihs_exbe=0.0_R8
      do iy = 0, ny-1
       do ix = 0, nx-1
        if(leftix(ix,iy).ne.-2 .and. rightix(ix,iy).ne.nx+1) then     !srv 26.11.02
!         wrk0(ix,iy) = bb(ix,iy,3)*vol(ix,iy)/hx(ix,iy)/hy(ix,iy)*(
!     &      (po(topix(ix,iy),topiy(ix,iy))-
!     &       po(bottomix(ix,iy),bottomiy(ix,iy)))
!     &     *(OnedBsq(rightix(ix,iy),rightiy(ix,iy))-
!     &       OnedBsq(leftix(ix,iy),leftiy(ix,iy)))
!     &     -(po(rightix(ix,iy),rightiy(ix,iy))-
!     &       po(leftix(ix,iy),leftiy(ix,iy)))
!     &     *(OnedBsq(topix(ix,iy),topiy(ix,iy))-
!     &       OnedBsq(bottomix(ix,iy),bottomiy(ix,iy)))
!     &                                                )*0.25e0_R8*signmf !srv 21.04.00
         wrk0(ix,iy) = bb(ix,iy,2)*vol(ix,iy)*                        !srv 11.01.13
     &   ((dpo(ix,iy,1,1)+dpo(topix(ix,iy),topiy(ix,iy),1,1))*
     &   (dOnedBsq(ix,iy,0,0)+
     &    dOnedBsq(rightix(ix,iy),rightiy(ix,iy),0,0)) -
     &    (dpo(ix,iy,0,0)+dpo(rightix(ix,iy),rightiy(ix,iy),0,0))*
     &   (dOnedBsq(ix,iy,1,1)+dOnedBsq(topix(ix,iy),topiy(ix,iy),1,1)))
     &    *0.25e0_R8 ! *signmf                                        !srv 11.01.13 03.04.07 }
         t0 = phm5*fac_ExB(ix,iy)*wrk0(ix,iy)                         !srv 07.06.00 02.06.99
         shedd(ix,iy,0) = shedd(ix,iy,0)+
     &     (max(0.0e0_R8,t0)+rf0*abs(t0))*ne(ix,iy)*te(ix,iy)
         shedd(ix,iy,3) = shedd(ix,iy,3)+min(0.0e0_R8,t0)-rf0*abs(t0)
         b2exbe(region(ix,iy,0))=b2exbe(region(ix,iy,0))+
     1       t0*ne(ix,iy)*te(ix,iy)
         b2sihs_exbe(ix,iy)=t0*ne(ix,iy)*te(ix,iy)
        endif                                                         !srv 26.11.02
       enddo
      enddo
      if (balance_netcdf.ne.0) then !djm Jan2017
        b2sihs_exbe0to3 = shedd
      endif
      call b2saxpy (n2*4, 1.0e0_R8, shedd, 1, she0, 1)                !srv 11.09.09

*   ..contribution from divergence of the electrical drift            !srv 02.06.99                                               !srv 02.06.99
      b2exba=0.0_R8
      b2sihs_exba=0.0_R8
      do is = 0, ns-1
       if(.not.is_neutral(is)) then
        do iy = 0, ny-1
         do ix = 0, nx-1
          if(leftix(ix,iy).ne.-2 .and. rightix(ix,iy).ne.nx+1) then   !srv 26.11.02
           t0=phm6*fac_ExB(ix,iy)*
     &      (na(ix,iy,is)/ni(ix,iy,0))*wrk0(ix,iy)                    !srv 07.06.00 02.06.99
           shidd(ix,iy,0) = shidd(ix,iy,0)+
     &      (max(0.0e0_R8,t0)+rf1*abs(t0))*ni(ix,iy,0)*ti(ix,iy)
!          shi0(ix,iy,0) = shi0(ix,iy,0)+
!     &      (max(0.0e0_R8,t0)+rf1*abs(t0))*ni(ix,iy,0)*ti(ix,iy)
           shidd(ix,iy,3) = shidd(ix,iy,3)+min(0.0e0_R8,t0)-rf1*abs(t0)
!          shi0(ix,iy,3) = shi0(ix,iy,3)+min(0.0e0_R8,t0)-rf1*abs(t0)
           b2exba(region(ix,iy,0))=b2exba(region(ix,iy,0))+
     1          t0*ni(ix,iy,0)*ti(ix,iy)
           b2sihs_exba(ix,iy)=b2sihs_exba(ix,iy)+
     1          t0*ni(ix,iy,0)*ti(ix,iy)
          endif                                                       !srv 26.11.02
         enddo
        enddo
       endif
      enddo
      if (balance_netcdf.ne.0) then !djm Jan2017
        b2sihs_exba0to3 = shidd
      endif
      call b2saxpy (n2*4, 1.0e0_R8, shidd, 1, shi0, 1)                !srv 11.09.09
*   ..contribution from viscous heating
      if ((Zhdanov_closure.eq.1) .and. (Zhdanov_test.eq.0).and.
     &  Zhdanov_visu.eq.1) then                                       !som 18.08.21
*   ..re-calculate mass-average velocity
        call b2mod_zhfrtf_avm_u(nx, ny, ns, na, ua)
*   ..re-calculate isonuclear density
        call b2mod_zhfrtf_prep(nx,ny,ns,na,ua)                        !som 09.05.24
      endif
      b2visa=0.0_R8
      b2sihs_visa=0.0_R8
      ft = 4.0e0_R8/3.0e0_R8                                          !srv 16.10.17
      do is = 0, ns-1
       if(.not.is_neutral(is).or.use_eirene.eq.0) then                !srv 29.08.02
cwdk compute radial differences of ua at cell faces
        call diff_r(nx,ny,0,hx,hy,vol,qs,ua(-1,-1,is),wrk0,duar)
        wrk0 = 0.0e0_R8
        if ((Zhdanov_closure.eq.0) .or. (Zhdanov_test.eq.1).or.
     &  Zhdanov_visu.eq.0 .or. is_neutral(is)) then                   !som 16.08.21
         do iy = -1, ny                                               !srv 28.08.02 21.07.05 16.10.17 {
          do ix = -1, nx
           if(leftix(ix,iy).eq.-2) then
             wrk0(ix,iy)=0.0e0_R8
           else
             wrk0(ix,iy) = cvsa(ix,iy,0,0,is)*                        !srv 16.10.17
     &         (ua(ix,iy,is)-ua(leftix(ix,iy),leftiy(ix,iy),is))**2
           endif
          enddo
         enddo
        else
         if (zh_vis_type.eq.0) then                                   !som 13.01.22
          do iy = -1, ny                                              !srv 28.08.02 21.07.05 16.10.17 {
           do ix = -1, nx
            if(leftix(ix,iy).eq.-2) then
             wrk0(ix,iy)=0.0e0_R8
            else
             wrk0(ix,iy) = cvsa(ix,iy,0,0,is)*                        !srv 16.10.17
     &         (ua(ix,iy,is)-ua(leftix(ix,iy),leftiy(ix,iy),is))**2
            endif
           enddo
          enddo
         else if (zh_vis_type.eq.1) then
          do iy = -1, ny
           do ix = -1, nx
            if(leftix(ix,iy).eq.-2) then
             wrk0(ix,iy)=0.0e0_R8
            else
             wrk0(ix,iy) = (cvsa(ix,iy,0,0,is)-
     &                ft*cvsa_cl(ix,iy,0,0,is))*
     &        (ua(ix,iy,is)-ua(leftix(ix,iy),leftiy(ix,iy),is))**2+
     &         ft*cvsa_drho(ix,iy,0,0,is)*
     &        (avm_u(ix,iy)-avm_u(leftix(ix,iy),leftiy(ix,iy)))**2    !som 16.08.21
            endif
           enddo
          enddo
         else if (zh_vis_type.eq.2) then                              !som 08.05.24
          do iy = -1, ny
           do ix = -1, nx
            if(leftix(ix,iy).eq.-2) then
             wrk0(ix,iy)=0.0e0_R8
            else
             wrk0(ix,iy) = cvsa(ix,iy,0,0,is)*
     &         (ua(ix,iy,is)-ua(leftix(ix,iy),leftiy(ix,iy),is))**2

             do inucl=0,nspecies-1                                     !som 09.05.24
              if(inucl.ne.s2nucl(is)-1) then
               t0=0.0e0_R8
               do iz= 1, znnucl_int(inucl)
                t0 = t0 + 0.5e0_R8*(na(ix,iy,nucl2s(inucl,iz))*
     &          ti(ix,iy)+
     &          na(leftix(ix,iy),leftiy(ix,iy),nucl2s(inucl,iz))*
     &          ti(leftix(ix,iy),leftiy(ix,iy)))*
     &          (ua(ix,iy,nucl2s(inucl,iz))
     &          -ua(leftix(ix,iy),leftiy(ix,iy),nucl2s(inucl,iz)))
               enddo

               wrk0(ix,iy) = wrk0(ix,iy) +
     &         (na(ix,iy,is)/nal(ix,iy,s2nucl(is)-1)+
     &         na(leftix(ix,iy),leftiy(ix,iy),is)/
     &         nal(leftix(ix,iy),leftiy(ix,iy),s2nucl(is)-1))*
     &         0.5e0_R8*cvsa_uAdp_albe(ix,iy,0,0,s2nucl(is)-1,inucl)*
     &         t0*
     &         (ua(ix,iy,is)-ua(leftix(ix,iy),leftiy(ix,iy),is))
              endif
             enddo

             t0=0.0e0_R8                                               !som 09.05.24
             do iz= 1, znnucl_int(s2nucl(is)-1)
              if(iz.ne.nint(rza(ix,iy,is))) then
               t0 = t0 + 0.5e0_R8*(na(ix,iy,nucl2s(s2nucl(is)-1,iz))*
     &         ti(ix,iy)+
     &         na(leftix(ix,iy),leftiy(ix,iy),
     &         nucl2s(s2nucl(is)-1,iz))*
     &         ti(leftix(ix,iy),leftiy(ix,iy)))*
     &         (ua(ix,iy,nucl2s(s2nucl(is)-1,iz)) -
     &         ua(leftix(ix,iy),leftiy(ix,iy),
     &         nucl2s(s2nucl(is)-1,iz)))
              endif
             enddo

             wrk0(ix,iy) = wrk0(ix,iy) +
     &        (na(ix,iy,is)/nal(ix,iy,s2nucl(is)-1)+
     &        na(leftix(ix,iy),leftiy(ix,iy),is)/
     &        nal(leftix(ix,iy),leftiy(ix,iy),s2nucl(is)-1))*
     &        0.5e0_R8*
     &        (cvsa_uAdp_albe(ix,iy,0,0,s2nucl(is)-1,s2nucl(is)-1)
     &        -zhcscorr*cvsa_uBdp_al(ix,iy,0,0,s2nucl(is)-1))*
     &        t0*
     &        (ua(ix,iy,is)-ua(leftix(ix,iy),leftiy(ix,iy),is))

            endif
           enddo
          enddo
         endif
        endif
        do iy = -1, ny
         do ix = -1, nx
          if(bottomiy(ix,iy).eq.-2) then
            wrk1(ix,iy) = 0.0e0_R8
          else
            wrk1(ix,iy) = cvsa(ix,iy,1,1,is)*duar(ix,iy,1)**2         !srv 16.10.17
          endif
         enddo
        enddo                                                         !srv 28.08.02 21.07.05 16.10.17 }
        do iy = 0, ny-1
         do ix = 0, nx-1
          if(leftix(ix,iy).ne.-2 .and. rightix(ix,iy).ne.nx+1) then   !srv 26.11.02
           t0 = phm2*(1.0_R8-BoRiS)*
     &       (wrk0(ix,iy)+wrk0(rightix(ix,iy),rightiy(ix,iy)))        !srv 21.07.05 }
     &           *0.5e0_R8
           t1 = phm8*(1.0_R8-BoRiS)*!                                 vsa0(ix,iy,is)*vol(ix,iy)*
     &       (wrk1(ix,iy)+wrk1(topix(ix,iy),topiy(ix,iy)))
     &           *0.5e0_R8                                            !srv 21.07.05 16.10.17 }
           shivc(ix,iy,0) = shivc(ix,iy,0)+(1.0e0_R8+rf2)*t0
           shivc(ix,iy,3) = shivc(ix,iy,3)-
     &                        rf2*t0/(ni(ix,iy,0)*ti(ix,iy))
           shiva(ix,iy,0) = shiva(ix,iy,0)+(1.0e0_R8+rf2)*t1
           shiva(ix,iy,3) = shiva(ix,iy,3)-
     &                        rf2*t1/(ni(ix,iy,0)*ti(ix,iy))
           b2visa(region(ix,iy,0))=b2visa(region(ix,iy,0))+t0+t1
           b2sihs_visa(ix,iy)=b2sihs_visa(ix,iy)+t0+t1
          endif                                                       !srv 26.11.02
         enddo
        enddo
       endif                                                          !srv 29.08.02
      enddo
      if (balance_netcdf.ne.0) then !djm Jan2017
        b2sihs_visa0to3 = shivc+shiva
      endif
      call b2saxpy (n2*4, 1.0e0_R8, shivc, 1, shi0, 1)                !srv 11.09.09
      call b2saxpy (n2*4, 1.0e0_R8, shiva, 1, shi0, 1)                !srv 11.09.09

*   ..contribution from electron-ion and ion-ion friction in case of arbitrary impurity amount
*     in case of small impurity amount see below for b2sigp_style.ne.2
      if(b2sigp_style.eq.2) then                                      !srv 05.07.17
       call b2xpne (nx, ny, ns, rz2, na, ne2_ext, ne2)
       call b2xzef (nx, ny, ns, rz2, na, ne, zeff)                    !som 13.11.20 {

       do iy = -1, ny
         do ix = -1, nx
           ce1(ix,iy)=fce1(zeff(ix,iy))
           ce2(ix,iy)=fce2(zeff(ix,iy))
         enddo
       enddo

       if (Zhdanov_closure.eq.0 .or. Zhdanov_test.eq.1 .or.
     &     Zhdanov_tf_fr.eq.0) then                                   !som 01.07.21
         mm(0:ns-1) = mass_multiplicator
         mm(ismain) = 1.0_R8
         am1(0:ns-1)= mm(0:ns-1)*am(0:ns-1)
         call b2xzmf (nx, ny, ns, ismain, rz2, na, zmff)
         call b2xzmc (nx, ny, ns, ismain, rz2, na, am1, zmc)
         call b2xzf1 (nx, ny, ns, ismain, rz2, na, am1, zf1)
         call b2xzf2 (nx, ny, ns, ismain, rz2, na, am1, zf2)
         call b2xz11 (nx, ny, ns, ismain, rz2, na, am1, z11)
         call b2xz12 (nx, ny, ns, ismain, rz2, na, am1, z12)
         call b2xz22 (nx, ny, ns, ismain, rz2, na, am1, z22)

         if (impr_form_fr.eq.0) then
           do iy = -1, ny
             do ix = -1, nx
               cimp1(ix,iy,0:ns-1)=fce1(zmff(ix,iy))
             enddo
           enddo
         else
           do iy = -1, ny
             do ix = -1, nx
               do is = 0, ns-1
                 cimp1(ix,iy,is)=fcim1(ismain,is,ns,am1,z11(ix,iy),
     &            z12(ix,iy),z22(ix,iy),zf1(ix,iy),zf2(ix,iy))
               enddo
             enddo
           enddo
         endif

         if (impr_form_tf.eq.0) then
           do iy = -1, ny
             do ix = -1, nx
               cimp2(ix,iy,0:ns-1)=fce2(zmff(ix,iy))
             enddo
           enddo
         else
           do iy = -1, ny
             do ix = -1, nx
               do is = 0, ns-1
                 cimp2(ix,iy,is)=fcim2(ismain,is,ns,am1,zmc(ix,iy),
     &            z11(ix,iy),z12(ix,iy),z22(ix,iy))
               enddo
             enddo
           enddo
         endif
       endif                                                          !som 13.11.20 }

*   ..compute the Coulomb logarithm
       call b2tlnl (nx, ny, icase, te, ti, ne, lnlam)
*   ..compute zetap and zetae
       t1 = (4.0_R8*pi*eps0/qe)/qe
       t0 = (3.0_R8/(4.0_R8*sqrt(2.0_R8*pi)))*t1*(t1*ev*sqrt(ev))
       czetap = t0*sqrt(mp)
       czetae = t0*sqrt(me)
       zetap = czetap*(ti/ev)*sqrt(ti/ev)/lnlam
       zetae = czetae*(te/ev)*sqrt(te/ev)/lnlam
*   ..compute gti,gte
       if (Zhdanov_closure.eq.0 .or. Zhdanov_test.eq.1 .or.
     &     Zhdanov_tf_fr.eq.0) then                                    !som 01.07.21
        do iy = -1, ny
         do ix = -1, nx
          if(leftix(ix,iy).eq.-2) then
           gti(ix,iy) = 0.0_R8
          else
           gti(ix,iy) = (ti(ix,iy)-ti(leftix(ix,iy),leftiy(ix,iy)))/
     &       (hx(ix,iy)+hx(leftix(ix,iy),leftiy(ix,iy)))*2.0e0_R8
          endif
         enddo
        enddo
       endif
       do iy = -1, ny
        do ix = -1, nx
         if(leftix(ix,iy).eq.-2) then
          gte(ix,iy) = 0.0_R8
         else
          gte(ix,iy) = (te(ix,iy)-te(leftix(ix,iy),leftiy(ix,iy)))/
     &      (hx(ix,iy)+hx(leftix(ix,iy),leftiy(ix,iy)))*2.0e0_R8
         endif
        enddo
       enddo
       shifr = 0.0_R8
       shefr = 0.0_R8
       do is = 0,ns-1
         if(.not. is_neutral(is)) then
           call b2sifrtf (nx, ny, ns, is, ismain, bb, na, ua,
     &      ue, te, ti, ni, ne, ne2,
     &      zetap, zetae, gti, gte, ce1, ce2, cimp1, cimp2,           !srv 05.07.17
     &      zeff, smbfrea, smbfria, smbtfea, smbtfia,                 !srv 05.07.17
     &      smbch, smbtf)
           wrks =
     &      (smbfrea(:,:,0)+
     &       smbfrea(:,:,1)*ua(:,:,is)+
     &       smbfrea(:,:,2)*na(:,:,is)*mp*am(is)+
     &       smbfrea(:,:,3)*ua(:,:,is)*na(:,:,is)*mp*am(is))+
     &      (smbtfea(:,:,0)+
     &       smbtfea(:,:,1)*ua(:,:,is)+
     &       smbtfea(:,:,2)*na(:,:,is)*mp*am(is)+
     &       smbtfea(:,:,3)*ua(:,:,is)*na(:,:,is)*mp*am(is))
           shefr(:,:,0) = shefr(:,:,0) + wrks
           shifr(:,:,0) = shifr(:,:,0) + ua(:,:,is)*(wrks+
     &      (smbfria(:,:,0)+
     &       smbfria(:,:,1)*ua(:,:,is)+
     &       smbfria(:,:,2)*na(:,:,is)*mp*am(is)+
     &       smbfria(:,:,3)*ua(:,:,is)*na(:,:,is)*mp*am(is))+
     &      (smbtfia(:,:,0)+
     &       smbtfia(:,:,1)*ua(:,:,is)+
     &       smbtfia(:,:,2)*na(:,:,is)*mp*am(is)+
     &       smbtfia(:,:,3)*ua(:,:,is)*na(:,:,is)*mp*am(is))
     &                                               )
         endif
       enddo
*   ..add contribution due to Zhdanov multi-species approach
       if(Zhdanov_closure.eq.1.and.Zhdanov_test.eq.0.and.
     &    Zhdanov_tf_fr.eq.1) then                                    !som 01.07.21
         shifr(:,:,0) = shifr(:,:,0)-zhd_corr*avm_u*shefr(:,:,0)
         shefr(:,:,0) = shefr(:,:,0)-zhd_corr*(avm_u/ue)*shefr(:,:,0)
       endif
       shefr(:,:,0) = ue*shefr(:,:,0)

       b2joule=0.0_R8
       b2sihs_joule=0.0_R8
       b2fraa=0.0_R8
       b2sihs_fraa=0.0_R8
       do iy = 0,ny-1
         do ix =0,nx-1
           if(leftix(ix,iy).ne.-2 .and. rightix(ix,iy).ne.nx+1 .and.
     &      bottomiy(ix,iy).ne.-2 .and. topiy(ix,iy).ne.ny+1 .and.
     &      region(ix,iy,0).ne.0) then
             t0 = phm3*vol(ix,iy)*shefr(ix,iy,0)
!            shefr(ix,iy,0) = t0
             shefr(ix,iy,0) = max(0.0e0_R8,t0)+rf3*abs(t0)
             shefr(ix,iy,3) =
     &        (min(0.0e0_R8,t0)-rf3*abs(t0))/(ne(ix,iy)*te(ix,iy))
             b2joule(region(ix,iy,0))=b2joule(region(ix,iy,0))+t0
             b2sihs_joule(ix,iy)=t0
             t0 = - phm4*(1.0_R8-BoRiS)*vol(ix,iy)*shifr(ix,iy,0)
!            shifr(ix,iy,0) = t0
!             shifr(ix,iy,0) = (2.5e0_R8+rf4)*t0
!             shifr(ix,iy,3) = -(1.5e0_R8+rf4)*t0/(ni(ix,iy,0)*ti(ix,iy))
             shifr(ix,iy,0) = max(0.0e0_R8,t0)+rf4*abs(t0)
             shifr(ix,iy,3) =
     &        (min(0.0e0_R8,t0)-rf4*abs(t0))/(ni(ix,iy,0)*ti(ix,iy))
             b2fraa(region(ix,iy,0))=b2fraa(region(ix,iy,0))+t0
             b2sihs_fraa(ix,iy)=t0
           endif
         enddo
       enddo
      endif                                                           !srv 05.07.17

*   ..contribution from electron-ion friction
*     (a.k.a. Joule heating) in case of small impurity amount;
*     in case of arbitrary impurity amount see above for b2sigp_style.eq.2
      if(b2sigp_style.ne.2) then                                      !srv 05.07.17
       call b2xehx (nx, ny, te, po, ne, ehx)
       if (istyle_Joule_heating.eq.0 .or. istyle_Joule_heating.eq.3)
     &   call b2xehy (nx, ny, te, po, ne, ehy)                        !srv 13.01.17
       b2joule=0.0_R8
       b2sihs_joule=0.0_R8
       do iy = 0, ny-1
        do ix = 0, nx-1
         if(leftix(ix,iy).ne.-2 .and. rightix(ix,iy).ne.nx+1 .and.
     &     bottomiy(ix,iy).ne.-2 .and. topiy(ix,iy).ne.ny+1 .and.
     &     region(ix,iy,0).ne.0) then
          if (istyle_Joule_heating.le.1) then
           t0 = fch(ix,iy,0,0)*ehx(ix,iy,0)/
     &       (vol(leftix(ix,iy),leftiy(ix,iy))+vol(ix,iy))+
     &        fch(rightix(ix,iy),rightiy(ix,iy),0,0)*
     &        ehx(rightix(ix,iy),rightiy(ix,iy),0)/
     &       (vol(ix,iy)+vol(rightix(ix,iy),rightiy(ix,iy)))
          else
           t0 = fch(ix,iy,0,0)*ehx(ix,iy,0)/(gs(ix,iy,0)*qc(ix,iy,0)*
     &      (hx(leftix(ix,iy),leftiy(ix,iy))+hx(ix,iy))) +
     &      fch(rightix(ix,iy),rightiy(ix,iy),0,0)*
     &      ehx(rightix(ix,iy),rightiy(ix,iy),0)/(
     &      gs(rightix(ix,iy),rightiy(ix,iy),0)*
     &      qc(rightix(ix,iy),rightiy(ix,iy),0)*
     &      (hx(ix,iy)+hx(rightix(ix,iy),rightiy(ix,iy))))
          endif
          if (istyle_Joule_heating.eq.0) t0 = t0 +
     &      fch(ix,iy,1,1)*ehy(ix,iy,1)/
     &      (vol(bottomix(ix,iy),bottomiy(ix,iy))+vol(ix,iy))+
     &       fch(topix(ix,iy),topiy(ix,iy),1,1)*
     &       ehy(topix(ix,iy),topiy(ix,iy),1)/
     &      (vol(ix,iy)+vol(topix(ix,iy),topiy(ix,iy)))
          if (istyle_Joule_heating.eq.2) t0 = t0 +
     &      fch(ix,iy,1,1)*ehy(ix,iy,1)/(gs(ix,iy,1)*
     &     (hy1(bottomix(ix,iy),bottomiy(ix,iy))+hy1(ix,iy))) +
     &      fch(topix(ix,iy),topiy(ix,iy),1,1)*
     &      ehy(topix(ix,iy),topiy(ix,iy),1)/(
     &      gs(topix(ix,iy),topiy(ix,iy),1)*
     &     (hy1(ix,iy)+hy1(topix(ix,iy),topiy(ix,iy))))
          t0 = phm3*vol(ix,iy)*t0
          shefr(ix,iy,0) = shefr(ix,iy,0)+max(0.0e0_R8,t0)+rf3*abs(t0)
          shefr(ix,iy,3) = shefr(ix,iy,3)+
     &    (min(0.0e0_R8,t0)-rf3*abs(t0))/(ne(ix,iy)*te(ix,iy))
          b2joule(region(ix,iy,0))=b2joule(region(ix,iy,0))+t0
          b2sihs_joule(ix,iy)=t0
         endif                                                        !srv 26.11.02
        enddo
       enddo
      endif                                                           !srv 05.07.17
      if (balance_netcdf.ne.0) then !djm Jan2017
        b2sihs_joule0to3 = shefr
      endif
      call b2saxpy (n2*4, 1.0e0_R8, shefr, 1, she0, 1)                !srv 11.09.09
*   ..contribution from atom-atom friction
      if(b2sigp_style.ne.2) then                                      !srv 05.07.17
       b2fraa=0.0_R8
       b2sihs_fraa=0.0_R8
       do isb = 0, ns-1
        do is = isb+1, ns-1
*     ..consider distinct cases
         if (.not.is_neutral(is).and..not.is_neutral(isb)) then
*      ..contribution from friction between charged species in case of small impurity amount;
*       in case of arbitrary impurity amount see above for b2sigp_style.eq.2
          coef = phm4/ctaup*mp*sqrt(am(isb)*am(is)/(am(isb)+am(is)))  !srv 01.07.09 11.01.12
          do iy = 0, ny-1
           do ix = 0, nx-1
            if(leftix(ix,iy).ne.-2 .and. rightix(ix,iy).ne.nx+1) then !srv 26.11.02
             t0 = rz2(ix,iy,isb)*rz2(ix,iy,is)*
     &          (coef*lnlam(ix,iy)/4.0_R8)*(1.0_R8-BoRiS)*            !srv 11.01.12
     &           na(ix,iy,is)*na(ix,iy,isb)*vti32(ix,iy)*
     &           ((ua(ix,iy,is)-ua(ix,iy,isb))**2)                    !srv 22.07.05
             shifr(ix,iy,0) = shifr(ix,iy,0)+(2.5e0_R8+rf4)*t0
             shifr(ix,iy,3) = shifr(ix,iy,3)-
     &        (1.5e0_R8+rf4)*t0/(ni(ix,iy,0)*ti(ix,iy))
             b2fraa(region(ix,iy,0))=b2fraa(region(ix,iy,0))+t0
             b2sihs_fraa(ix,iy)=b2sihs_fraa(ix,iy)+t0
            endif                                                     !srv 26.11.02
           enddo
          enddo
         else
*      ..contribution from friction involving neutral species
*        (We ignore such a term, except for charge exchange friction
*        between H0 and H1. That term is calculated via the atomic
*        physics routines.)
         endif
        enddo
       enddo
      endif                                                           !srv 05.07.17
      if (balance_netcdf.ne.0) then !djm Jan2017
        b2sihs_fraa0to3 = shifr
      endif
      call b2saxpy (n2*4, 1.0e0_R8, shifr, 1, shi0, 1)                !srv 11.09.09
      b2str=0.0_R8
      b2sihs_str=0.0_R8
      b2wrong1=0.0_R8
* ..return
      if(iout.eq.1) then
        do k=0,3
          write (chk,'(i1)') k
          call my_out(70,nx,ny,shedu(-1,-1,k),'b2sihs__shedu'//chk)
          call my_out(70,nx,ny,shedd(-1,-1,k),'b2sihs__shedd'//chk)
          call my_out(70,nx,ny,shefr(-1,-1,k),'b2sihs__shefr'//chk)
!         call my_out(70,nx,ny,she0(-1,-1,k),'b2sihs__she'//chk)

          call my_out(70,nx,ny,shidu(-1,-1,k),'b2sihs__shidu'//chk)
          call my_out(70,nx,ny,shidd(-1,-1,k),'b2sihs__shidd'//chk)
          call my_out(70,nx,ny,shifr(-1,-1,k),'b2sihs__shifr'//chk)
          call my_out(70,nx,ny,shiva(-1,-1,k),'b2sihs__shiva'//chk)
          call my_out(70,nx,ny,shivc(-1,-1,k),'b2sihs__shivc'//chk)
!         call my_out(70,nx,ny,shi0(-1,-1,k),'b2sihs__shi'//chk)
        enddo
      endif
      if (iout_b2nph9.eq.1 .or. iout_b2npht.eq.1 .or. iout_b2wdat.eq.4)
     &                                            then
        wrks=she0(:,:,0)+she0(:,:,1)*te+
     &       she0(:,:,2)*ne + she0(:,:,3)*te*ne
        call my_out(70,nx,ny,wrks,'b2sihs__she')
        wrks=shest(:,:,0)+shest(:,:,1)*te+
     &       shest(:,:,2)*ne + shest(:,:,3)*te*ne
        call my_out(70,nx,ny,shest,'b2sihs__shest')        !srv 17.12.13
        wrks=shedu(:,:,0)+shedu(:,:,1)*te+
     &       shedu(:,:,2)*ne + shedu(:,:,3)*te*ne
        call my_out(70,nx,ny,wrks,'b2sihs__shedu')
        wrks=shedd(:,:,0)+shedd(:,:,1)*te+
     &       shedd(:,:,2)*ne + shedd(:,:,3)*te*ne
        call my_out(70,nx,ny,wrks,'b2sihs__shedd')
        wrks=shefr(:,:,0)+shefr(:,:,1)*te+
     &  shefr(:,:,2)*ne + shefr(:,:,3)*te*ne
        call my_out(70,nx,ny,wrks,'b2sihs__shefr')

        wrks=shi0(:,:,0)+shi0(:,:,1)*ti+
     &       shi0(:,:,2)*ni(:,:,0)+shi0(:,:,3)*ti*ni(:,:,0)
        call my_out(70,nx,ny,wrks,'b2sihs__shi')
        wrks=shidu(:,:,0)+shidu(:,:,1)*ti+
     &       shidu(:,:,2)*ni(:,:,0)+shidu(:,:,3)*ti*ni(:,:,0)
        call my_out(70,nx,ny,wrks,'b2sihs__shidu')
        wrks=shidd(:,:,0)+shidd(:,:,1)*ti+
     &       shidd(:,:,2)*ni(:,:,0)+shidd(:,:,3)*ti*ni(:,:,0)
        call my_out(70,nx,ny,wrks,'b2sihs__shidd')
        wrks=shifr(:,:,0)+shifr(:,:,1)*ti+
     &       shifr(:,:,2)*ni(:,:,0)+shifr(:,:,3)*ti*ni(:,:,0)
        call my_out(70,nx,ny,wrks,'b2sihs__shifr')
        wrks=shiva(:,:,0)+shiva(:,:,1)*ti+
     &       shiva(:,:,2)*ni(:,:,0)+shiva(:,:,3)*ti*ni(:,:,0)
        call my_out(70,nx,ny,wrks,'b2sihs__shiva')
        wrks=shivc(:,:,0)+shivc(:,:,1)*ti+
     &       shivc(:,:,2)*ni(:,:,0)+shivc(:,:,3)*ti*ni(:,:,0)
        call my_out(70,nx,ny,wrks,'b2sihs__shivc')
      endif
      ncall = ncall+1
      call subend ()
      return

*-----------------------------------------------------------------------
*.end b2sihs_

      end subroutine b2sihs_
