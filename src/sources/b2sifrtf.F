*-----------------------------------------------------------------------
*.specification

      subroutine b2sifrtf (nx, ny, ns, isb, ismain, bb, na,             !srv 09.01.01 18.10.02 13.01.17 05.07.17
     &   ua, ue, te, ti, ni, ne, ne2, zetap, zetae, gti, gte, ce1, ce2, !srv 01.07.09 11.09.09 05.07.17
     &   cimp1, cimp2, zeff, smbfrea, smbfria, smbtfea, smbtfia,        !srv 05.07.17
     &   smbch, smbtf)
      use b2mod_types
      use b2mod_rates
      use b2mod_indirect
      use b2mod_constants
      use b2mod_b2cmpa
      use b2mod_b2cmpt
      use b2mod_zhfrtf
     & , only: b2mod_zhfrtf_tf, b2mod_zhfrtf_fr
      use b2mod_b2zhco
     & , only: zhd_corr, zhd_corrT                                      !som 05.03.21
      use b2mod_plasma
     , , only : roxa
      use b2mod_transport
     , , only : f_luc_sg, alfx_c, sigx_c                                !srv 08.10.14
      use b2mod_external
      use b2mod_geo
     , , only : hz, vol                                                 !srv 05.07.17
      use b2mod_balance                                                 !djm Jun2019
     & , only : b2sifr_smofrea0to3, b2sifr_smofria0to3,
     &          b2sifr_smotfea0to3, b2sifr_smotfia0to3, balance_netcdf
      use b2mod_openmp
      use b2mod_subsys
      use b2mod_frtf_nccorr
     , , only : corr_fria, corr_tfia, alpha_HS_style, g_HS_style
#ifdef WITH_FTIMINGS
      use timings
#endif
#ifndef NO_IEEE
      use ieee_arithmetic    ! IGNORE
#endif
      implicit none
*   ..input arguments (unchanged on exit)
      integer nx, ny, ns, isb, ismain
      real (kind=R8) ::
     &  bb(-1:nx,-1:ny,0:3),
     &  na(-1:nx,-1:ny,0:ns-1), ua(-1:nx,-1:ny,0:ns-1),
     &  ue(-1:nx,-1:ny),
     &  ti(-1:nx,-1:ny), ni(-1:nx,-1:ny,0:1), ne(-1:nx,-1:ny),
     &  ne2(-1:nx,-1:ny), gti(-1:nx,-1:ny), gte(-1:nx,-1:ny),
     &  te(-1:nx,-1:ny), zetap(-1:nx,-1:ny), zetae(-1:nx,-1:ny),        !srv 13.01.17
     &  ce1(-1:nx,-1:ny), ce2(-1:nx,-1:ny),                             !srv 13.01.17
     &  cimp1(-1:nx,-1:ny,0:ns-1), cimp2(-1:nx,-1:ny,0:ns-1),           !som 13.11.20
     &  zeff(-1:nx,-1:ny)
*   ..output arguments (unspecified on entry)
      real (kind=R8) ::
     *  smbfria(-1:nx,-1:ny,0:3), smbfrea(-1:nx,-1:ny,0:3),
     *  smbtfia(-1:nx,-1:ny,0:3), smbtfea(-1:nx,-1:ny,0:3)              !srv 13.01.17 05.07.17
      real (kind=R8) ::
     *  smbch(-1:nx,-1:ny,0:3), smbtf(-1:nx,-1:ny,0:3)
*   ..common blocks
*-----------------------------------------------------------------------
*.documentation
*
*  1. purpose
*
*     B2SIFRTF computes a linearized representation of the friction
*     and thermal forces.
*
*
*  3. description (see also routine b2cdca)
*
*     This routine computes a linear representation of the force
*     exerted on one species by all others. This includes the ordinary
*     friction force proportional to the velocity difference, as well
*     as electron-atom and atom-atom thermal forces. It is an important
*     property that these forces add up to 0 when summed over all
*     species. This property must be maintained numerically.
*
*-----------------------------------------------------------------------
*.declarations

*   ..local variables
      integer ncall, styl0, n2, ix, iy, is
      integer, save :: iout = 0, iout_b2npmo = 0, iout_b2wdat = 0
      integer, save :: Zhdanov_closure = 0, Zhdanov_test = 0,           !som 23.02.21
     &  Zhdanov_tf_fr = 1                                               !som 01.07.21
      integer, save :: Zhdanov_nc = 0
      integer, save :: icase = 0                                        !som 05.03.21
#ifdef _OPENMP
!$OMP THREADPRIVATE(iout, iout_b2npmo, iout_b2wdat)
#endif
      character chns*3, chis*3
      real (kind=R8) ::
     *  phm0, phm1, phm2, phm3, coef, t0, calfab, fllim_al_ab, jabmax
      real (kind=R8) ::
     1  delti, delte
      real (kind=R8) ::
     &  KabVp(-1:nx,-1:ny), Kabtf(-1:nx,-1:ny), Kbatf(-1:nx,-1:ny),     !srv 13.01.17 {
     &  ka(-1:nx,-1:ny), kb(-1:nx,-1:ny),
     &  wrk(-1:nx,-1:ny)
! Local arrays smb*al have been introduced specially.
! The usage of smb*a instead of smb*al should be checked.
      real (kind=R8) ::
     *  smbfrial(-1:nx,-1:ny,0:3), smbfreal(-1:nx,-1:ny,0:3),
     *  smbtfial(-1:nx,-1:ny,0:3), smbtfeal(-1:nx,-1:ny,0:3),           !srv 13.01.17 }
     *  smbfrial_test(-1:nx,-1:ny,0:3),smbtfial_test(-1:nx,-1:ny,0:3),  !som 05.03.21
     *  smbtfial_nofl(-1:nx,-1:ny,0:3),                                 !som 06.03.21
     *  smbtfial_test_nofl(-1:nx,-1:ny,0:3)                             !som 06.03.21
      real (kind=R8) :: temp
      real (kind=R8) :: rz2_temp(1:ns*(nx+2)*(ny+2)),
     &                  na_temp(1:ns*(nx+2)*(ny+2)),
     &                  am_sqrt(0:ns-1, 0:ns-1)
      integer :: ctr, i, j
#if ( defined(WITH_FTIMINGS) || ( defined(_OPENMP) && !defined(NO_OPENMP_B2SIFRTF) ) )
      logical :: make_parallel
#endif
*   ..procedures
      intrinsic sqrt
      external xertst, ipgeti, ipgetr, sfill
      external b2xvff, b2xvfx, b2xvsg
*   ..initialization
      save ncall, styl0, phm0, phm1, phm2, phm3
#ifdef _OPENMP
!$OMP THREADPRIVATE(ncall, phm1, phm2)
#endif
      data ncall/0/, styl0/0/, phm0/1.0e0_R8/,
     &  phm1/1.0e0_R8/, phm2/1.0e0_R8/, phm3/1.0e0_R8/

*-----------------------------------------------------------------------
*.computation

* ..preliminaries
*   ..subprogram start-up calls
      call subini ('b2sifrtf')
*   ..set internal parameters on first call
      if (ncall.eq.0) then
       call ipgeti ('b2sifr_styl0', styl0)
       call ipgetr ('b2sifr_phm0', phm0)
       call ipgetr ('b2sifr_phm1', phm1)
       call ipgetr ('b2sifr_phm2', phm2)
       call ipgetr ('b2sifr_phm3', phm3)
       call xertst (0.le.styl0.and.styl0.le.2,
     &  'faulty internal parameter styl0')
       call xertst (0.0_R8.le.phm0.and.0.0_R8.le.phm1.and.
     &              0.0_R8.le.phm2.and.0.0_R8.le.phm3,
     &  'faulty internal parameter phm0, phm1, phm2, phm3')
       call ipgeti ('b2sifrtf_iout', iout)
       if ((iout.ne.0) .and. in_parallel()) then
#ifdef _OPENMP
!$OMP MASTER
#endif
          write(*,*) 'B2SIFRTF OpenMP warning: no file output in ',
     &               'parallel mode'
#ifdef _OPENMP
!$OMP END MASTER
#endif
          iout = 0
        end if
       call ipgeti ('b2wdat_iout', iout_b2wdat)                         !iys 17.04.17
       if ((iout_b2wdat==4) .and. in_parallel()) then
#ifdef _OPENMP
!$OMP MASTER
#endif
          write(*,*) 'B2SIFRTF OpenMP warning: no file output in ',
     &               'parallel mode (b2wdat)'
#ifdef _OPENMP
!$OMP END MASTER
#endif
          iout_b2wdat = 0
        end if
       call ipgeti ('b2npmo_iout', iout_b2npmo)
       call ipgeti ('b2tral_Zhdanov_closure', Zhdanov_closure)          !som 23.02.21
       call ipgeti ('b2sifr_Zhdanov_test', Zhdanov_test)                !som 12.03.21
       call ipgeti ('b2tral_Zhdanov_nc', Zhdanov_nc)
       call ipgeti ('b2sifr_Zhdanov_tf_fr', Zhdanov_tf_fr)              !som 01.07.21
       call ipgeti ('b2tral_zhd_corr', zhd_corr)                        !som 05.03.21
       call ipgeti ('b2tral_zhd_corrT', zhd_corrT)                      !som 03.09.21
       call ipgeti ('b2tlnl_ii', icase)                                 !som 05.03.21
       call xertst (0.eq.icase .or. 1.eq.icase,
     &              'faulty input b2tlnl_ii')
      endif
#if ( defined(WITH_FTIMINGS) || ( defined(_OPENMP) && !defined(NO_OPENMP_B2SIFRTF) ) )
      make_parallel = .not. in_parallel()
#endif
*   ..set n2
      n2 = (nx+2)*(ny+2)
*   ..test nx, ny
      call xertst (0.le.nx.and.0.le.ny, 'faulty argument nx, ny')
*   ..extensive tests on first few calls
      if (ncall.lt.3) then
*    ..test sign of vol
       call b2xvsg (n2, vol, 1, 'vol', '.gt.')
*    ..test sign of na, ni, ti, ne, ne2
       call b2xvsg (n2*ns, na, 1, 'na', '.gt.')
       call b2xvsg (n2*2, ni, 1, 'ni', '.gt.')
       call b2xvsg (n2, te, 1, 'te', '.gt.')
       call b2xvsg (n2, ti, 1, 'ti', '.gt.')
       call b2xvsg (n2, ne, 1, 'ne', '.gt.')
       call b2xvsg (n2, ne2, 1, 'ne2', '.gt.')
      endif

*   ..initialize to 0
      call sfill (n2*4, 0.0e0_R8, smbch, 1)                             !srv 11.09.09
      call sfill (n2*4, 0.0e0_R8, smbtf, 1)                             !srv 11.09.09
      call sfill (n2*4, 0.0e0_R8, smbfrial, 1)                          !srv 13.01.17 { !srv 05.07.17
      call sfill (n2*4, 0.0e0_R8, smbfreal, 1)
      call sfill (n2*4, 0.0e0_R8, smbtfial, 1)
      call sfill (n2*4, 0.0e0_R8, smbtfial_nofl, 1)
      call sfill (n2*4, 0.0e0_R8, smbtfeal, 1)                          !srv 13.01.17 }

*   ..make 3rd dimension contiguous by copying to 1-d array (Gaurav Saxena: BSC ACH)
      ctr = 1
      do iy = -1, ny
        do ix = -1, nx
#if ( defined(_OPENMP) && !defined(LEGACYCOMP) && !defined(NAGFOR) )
!$OMP SIMD
#endif
          do is = 0, ns-1
            rz2_temp(ctr) = rz2(ix,iy,is)
            na_temp(ctr) = na(ix,iy,is)
            ctr = ctr + 1
          end do
        end do
      end do
      do j = 0, ns-1
        do i = 0, ns-1
          am_sqrt(i,j) = sqrt(am(j)*am(i)/(am(j)+am(i)))
        end do
      end do

* ..compute friction
*   ..loop over all other species
#ifdef WITH_FTIMINGS
      if (make_parallel) call timer%start("parallel")
#endif
*     ..contribution from friction between ions
      if (Zhdanov_closure.eq.1.and.Zhdanov_test.eq.0.and.
     &    Zhdanov_tf_fr.eq.1) then                                      !som 01.07.2021
*     ..actual use of the Zhdanov closure
        do iy = 0, ny-1
          do ix = 0, nx-1
           if(leftix(ix,iy).ne.-2 .and. rightix(ix,iy).ne.nx+1) then
            temp=rz2(ix,iy,isb)*na(ix,iy,isb)*mp/zetap(ix,iy)
            call b2mod_zhfrtf_fr(ix,iy,isb,temp,
     &       smbfrial(ix,iy,0),smbfrial(ix,iy,1),
     &       corr_fria(ix,iy,isb))
           endif
          enddo
        enddo
      elseif (Zhdanov_closure.eq.1.and.(Zhdanov_test.eq.1)
     &  .and.Zhdanov_nc.eq.1) then
        do iy = 0, ny-1
          do ix = 0, nx-1
           if(leftix(ix,iy).ne.-2 .and. rightix(ix,iy).ne.nx+1) then
            temp=rz2(ix,iy,isb)*na(ix,iy,isb)*mp/zetap(ix,iy)
            call b2mod_zhfrtf_fr(ix,iy,isb,temp,
     &       smbfrial(ix,iy,0),smbfrial(ix,iy,1),
     &       corr_fria(ix,iy,isb))
           endif
          enddo
        enddo
      elseif (Zhdanov_closure.eq.1.and.(Zhdanov_test.eq.1.or.
     &        Zhdanov_tf_fr.eq.0)) then                                 !som 01.07.2021
*     ..test run of the Zhdanov closure
        smbfrial_test=0
        do iy = 0, ny-1
          do ix = 0, nx-1
           if(leftix(ix,iy).ne.-2 .and. rightix(ix,iy).ne.nx+1) then
            temp=rz2(ix,iy,isb)*na(ix,iy,isb)*mp/zetap(ix,iy)
            call b2mod_zhfrtf_fr(ix,iy,isb,temp,
     &       smbfrial_test(ix,iy,0),smbfrial_test(ix,iy,1),
     &       corr_fria(ix,iy,isb))
           endif
          enddo
        enddo
#if ( defined(_OPENMP) && !defined(NO_OPENMP_B2SIFRTF) )
!$OMP PARALLEL DO IF (make_parallel)
!$OMP& DEFAULT(NONE)
!$OMP& PRIVATE(KabVp, is, ix, iy, coef, t0)
!$OMP& SHARED(ns, isb, is_neutral, nx, ny, leftix, rightix)
!$OMP& SHARED(rz2, am, zetap, na, ua)
!$OMP& SHARED(ismain, cimp1)
!$OMP& FIRSTPRIVATE(styl0, phm0)
!$OMP& REDUCTION(+: smbfrial)
#endif
        do is = 0, ns-1
         if(isb.ne.is .and.
     &    .not.is_neutral(isb).and. .not.is_neutral(is)) then
          do iy = 0, ny-1
           do ix = 0, nx-1
            KabVp(ix,iy)=fkabvp(ix,iy,isb,is)
            if(leftix(ix,iy).ne.-2 .and. rightix(ix,iy).ne.nx+1) then   !srv 26.11.02
             coef = phm0*rz2(ix,iy,isb)*rz2(ix,iy,is)*
     &           mp*sqrt(am(isb)*am(is)/(am(isb)+am(is)))/
     &           zetap(ix,iy)*KabVp(ix,iy)
             t0 = na(ix,iy,isb)*na(ix,iy,is)   !   *vol(ix,iy)*hz(ix,iy) ! it will be multiplied in b2npmo !srv 05.07.17
             if (styl0.eq.0) then
               smbfrial(ix,iy,0) = smbfrial(ix,iy,0)+
     &          coef*t0*ua(ix,iy,is)
               smbfrial(ix,iy,1) = smbfrial(ix,iy,1)+(-coef*t0)
             else if (styl0.eq.1) then
               smbfrial(ix,iy,0) = smbfrial(ix,iy,0)+
     &          coef*t0*ua(ix,iy,is)
               smbfrial(ix,iy,3) = smbfrial(ix,iy,3)+(-coef*t0/
     &                       roxa(ix,iy,isb))
             else if (styl0.eq.2) then
               smbfrial(ix,iy,2) = smbfrial(ix,iy,2)+
     &          coef*t0*ua(ix,iy,is)/roxa(ix,iy,isb)
               smbfrial(ix,iy,3) = smbfrial(ix,iy,3)+(-coef*t0/
     &                       roxa(ix,iy,isb))
             endif
            endif
           enddo
          enddo
         endif
        enddo
      else
*     ..standard friction force calculation
#if ( defined(_OPENMP) && !defined(NO_OPENMP_B2SIFRTF) )
!$OMP PARALLEL DO IF (make_parallel)
!$OMP& DEFAULT(NONE)
!$OMP& PRIVATE(KabVp, is, ix, iy, coef, t0)
!$OMP& SHARED(ns, isb, is_neutral, nx, ny, leftix, rightix)
!$OMP& SHARED(rz2, am, zetap, na, ua)
!$OMP& SHARED(ismain, cimp1)
!$OMP& FIRSTPRIVATE(styl0, phm0)
!$OMP& REDUCTION(+: smbfrial)
#endif
        do is = 0, ns-1
         if(isb.ne.is .and.
     &    .not.is_neutral(isb).and. .not.is_neutral(is)) then
*       ..contribution from friction between ions
          do iy = 0, ny-1
           do ix = 0, nx-1
            KabVp(ix,iy)=fkabvp(ix,iy,isb,is)
            if(leftix(ix,iy).ne.-2 .and. rightix(ix,iy).ne.nx+1) then   !srv 26.11.02
             coef = phm0*rz2(ix,iy,isb)*rz2(ix,iy,is)*
     &           mp*sqrt(am(isb)*am(is)/(am(isb)+am(is)))/
     &           zetap(ix,iy)*KabVp(ix,iy)
             t0 = na(ix,iy,isb)*na(ix,iy,is)   !   *vol(ix,iy)*hz(ix,iy) ! it will be multiplied in b2npmo !srv 05.07.17
             if (styl0.eq.0) then
               smbfrial(ix,iy,0) = smbfrial(ix,iy,0)+
     &          coef*t0*ua(ix,iy,is)
               smbfrial(ix,iy,1) = smbfrial(ix,iy,1)+(-coef*t0)
             else if (styl0.eq.1) then
               smbfrial(ix,iy,0) = smbfrial(ix,iy,0)+
     &          coef*t0*ua(ix,iy,is)
               smbfrial(ix,iy,3) = smbfrial(ix,iy,3)+(-coef*t0/
     &                       roxa(ix,iy,isb))
             else if (styl0.eq.2) then
               smbfrial(ix,iy,2) = smbfrial(ix,iy,2)+
     &          coef*t0*ua(ix,iy,is)/roxa(ix,iy,isb)
               smbfrial(ix,iy,3) = smbfrial(ix,iy,3)+(-coef*t0/
     &                       roxa(ix,iy,isb))
             endif
            endif
           enddo
          enddo
         endif
        enddo
      endif
      wrk(-1:nx,-1:ny) = smbfrial(-1:nx,-1:ny,0)
      call nancheck(nx, ny, wrk, 'smbfrial(,,0)')
#ifdef WITH_FTIMINGS
      if (make_parallel) call timer%stop("parallel")
#endif
!*   ..compute ue, electron parallel velocity
!      call b2xpve (nx, ny, ns, qe, rza, rz2, alfx_c, sigx_c, po,       !srv 05.07.17
!     &             ne, te, na, ua, ue_ext, ne_ext, ue)
*     ..contribution from friction between between electron and ions
      if(.not.is_neutral(isb)) then
       do iy = 0, ny-1
        do ix = 0, nx-1
         if(leftix(ix,iy).ne.-2 .and. rightix(ix,iy).ne.nx+1) then
          coef = phm1*(bb(ix,iy,0)/bb(ix,iy,3))**2*rz2(ix,iy,isb)*qe**2
     &           /(sigx_c(ix,iy)*zeff(ix,iy))
          t0 = na(ix,iy,isb)*ne(ix,iy)   !   *vol(ix,iy)*hz(ix,iy) ! it will be multiplied in b2npmo !srv 05.07.17
          smbfreal(ix,iy,0) = smbfreal(ix,iy,0)+coef*t0*ue(ix,iy)
          smbfreal(ix,iy,1) = smbfreal(ix,iy,1)-coef*t0
         endif
        enddo
       enddo
      endif

      smbch = smbfreal + smbfrial

* ..compute thermal force

*   ..contribution from electron thermal force
      if(.not.is_neutral(isb)) then
       do iy = 0, ny-1
        do ix = 0, nx-1
         if(leftix(ix,iy).ne.-2 .and. rightix(ix,iy).ne.nx+1) then
          coef = -phm2*alfx_c(ix,iy)*ce1(ix,iy)*rz2(ix,iy,isb)*(me/qe)/
     &           zetae(ix,iy)/(bb(ix,iy,0)/bb(ix,iy,3))/f_luc_sg(ix,iy)   !srv f_luc_sg should be taken in cell center ?
          delte=(gte(ix,iy)+gte(rightix(ix,iy),rightiy(ix,iy)))*0.5e0_R8
          t0 = na(ix,iy,isb)*delte    !   *vol(ix,iy)*hz(ix,iy) ! it will be multiplied in b2npmo !srv 05.07.17
          smbtfeal(ix,iy,0) = smbtfeal(ix,iy,0)+coef*t0
         endif
        enddo
       enddo
      endif
*   ..contribution from ions thermal force
#ifdef WITH_FTIMINGS
      if (make_parallel) call timer%start("parallel")
#endif
      if(Zhdanov_closure.eq.1.and.Zhdanov_test.eq.0.and.
     &   Zhdanov_tf_fr.eq.1) then                                       !som 01.07.2021
*     ..actual use of the Zhdanov closure
          do iy = 0, ny-1
           do ix = 0, nx-1
            if(leftix(ix,iy).ne.-2 .and. rightix(ix,iy).ne.nx+1) then
             call b2mod_zhfrtf_tf(ix,iy,isb,na(ix,iy,isb),
     &        smbtfial(ix,iy,0),smbtfial_nofl(ix,iy,0),
     &        corr_tfia(ix,iy,isb))
            endif
           enddo
          enddo
      elseif (Zhdanov_closure.eq.1.and.(Zhdanov_test.eq.1.or.
     &        Zhdanov_tf_fr.eq.0)) then                                 !som 01.07.2021
*     ..test run of the Zhdanov closure
        smbtfial_test=0
        smbtfial_test_nofl=0
        do iy = 0, ny-1
         do ix = 0, nx-1
           if(leftix(ix,iy).ne.-2 .and. rightix(ix,iy).ne.nx+1) then
             call b2mod_zhfrtf_tf(ix,iy,isb,na(ix,iy,isb),
     &        smbtfial_test(ix,iy,0),smbtfial_test_nofl(ix,iy,0),
     &        corr_tfia(ix,iy,isb))
           endif
         enddo
        enddo
        if (Zhdanov_nc.eq.1) then
          call b2scopy (n2, smbtfial_test(-1,-1,0), 1,
     &    smbtfial(:,:,0), 1)
          call b2scopy (n2, smbtfial_test_nofl(-1,-1,0), 1,
     &    smbtfial_nofl(:,:,0), 1)
        else
        do iy = 0, ny-1
         do ix = 0, nx-1
           ka(ix,iy) = fka(ix,iy,isb)
         enddo
        enddo
#if ( defined(_OPENMP) && !defined(NO_OPENMP_B2SIFRTF) )
!$OMP PARALLEL DO IF (make_parallel)
!$OMP& DEFAULT(NONE)
!$OMP& PRIVATE(KabVp, Kbatf, Kabtf, kb, is, ix, iy, coef, t0)
!$OMP& PRIVATE(calfab, delti, jabmax, fllim_al_ab)
!$OMP& SHARED(ns, isb, is_neutral, nx, ny, leftix, rightix)
!$OMP& SHARED(rz2, am, zetap, na, ua, bb, gti, cflim, ti, rightiy)
!$OMP& SHARED(on_closed_surface, ka)
!$OMP& SHARED(rz2_temp, na_temp, am_sqrt)
!$OMP& SHARED(ismain, cimp1, cimp2)
!$OMP& FIRSTPRIVATE(phm3)
!$OMP& REDUCTION(+: smbtfial, smbtfial_nofl)
#endif
        do is = 0, ns-1
         if(isb.ne.is .and.
     &      .not.is_neutral(isb) .and. .not.is_neutral(is)) then
          do iy = 0, ny-1
           do ix = 0, nx-1
            if(leftix(ix,iy).ne.-2 .and. rightix(ix,iy).ne.nx+1) then
             Kabtf(ix,iy) = fkabtf(ix,iy,isb,is)
             Kbatf(ix,iy) = fkabtf(ix,iy,is,isb)
c            kb(ix,iy)    = fka(ix,iy,is)
c            kb(ix,iy)    = fka_new(ix,iy,is,rz2_temp,na_temp)
             kb(ix,iy)    =
     &                 fka_new_opt(ix,iy,is,rz2_temp,na_temp,am_sqrt)
             KabVp(ix,iy) = fkabvp(ix,iy,isb,is)
c        .. calculations could be simplified
c        .. but we did not do it in order to have the same formulae as manual
             calfab = zetap(ix,iy)/sqrt(mp)*
     &            qe*(na(ix,iy,isb)+na(ix,iy,is))*
     &            sqrt(am(isb)*am(is)/(am(isb)+am(is)))*
     &           (Kabtf(ix,iy)/kb(ix,iy)-Kbatf(ix,iy)/ka(ix,iy))*
     &            bb(ix,iy,0)/bb(ix,iy,3)
             coef = phm3*rz2(ix,iy,isb)*rz2(ix,iy,is)*calfab*
     &            mp/(zetap(ix,iy)*qe)
             delti=
     &         (gti(ix,iy)+gti(rightix(ix,iy),rightiy(ix,iy)))*0.5e0_R8
c   .... apply flux limit
             if(cflim(4).ne.0.0_R8.and.
     &         .not.on_closed_surface(ix,iy)) then
               jabmax = qe*KabVp(ix,iy)*sqrt(ti(ix,iy)/mp)*
     &          (na(ix,iy,isb)+na(ix,iy,is))
               t0 = abs(calfab*delti)/(cflim(4)*jabmax)
               fllim_al_ab = 1.0_R8/(1.0e0_R8+t0)
             else
               fllim_al_ab = 1.0_R8
             endif

             t0 = na(ix,iy,isb)*na(ix,iy,is)/    !   *vol(ix,iy)*hz(ix,iy) ! it will be multiplied in b2npmo !srv 05.07.17
     &                              ((na(ix,iy,isb)+na(ix,iy,is)))*delti
             smbtfial(ix,iy,0) = smbtfial(ix,iy,0) + fllim_al_ab*coef*t0
             smbtfial_nofl(ix,iy,0) = smbtfial_nofl(ix,iy,0) + coef*t0
            endif
           enddo
          enddo
         endif
        enddo
        endif
      else
*     ..standard thermal force calculation
        do iy = 0, ny-1
         do ix = 0, nx-1
           ka(ix,iy) = fka(ix,iy,isb)
         enddo
        enddo
#if ( defined(_OPENMP) && !defined(NO_OPENMP_B2SIFRTF) )
!$OMP PARALLEL DO IF (make_parallel)
!$OMP& DEFAULT(NONE)
!$OMP& PRIVATE(KabVp, Kbatf, Kabtf, kb, is, ix, iy, coef, t0)
!$OMP& PRIVATE(calfab, delti, jabmax, fllim_al_ab)
!$OMP& SHARED(ns, isb, is_neutral, nx, ny, leftix, rightix)
!$OMP& SHARED(rz2, am, zetap, na, ua, bb, gti, cflim, ti, rightiy)
!$OMP& SHARED(ismain, cimp1, cimp2)
!$OMP& SHARED(on_closed_surface, ka)
!$OMP& SHARED(rz2_temp, na_temp, am_sqrt)
!$OMP& FIRSTPRIVATE(phm3)
!$OMP& REDUCTION(+: smbtfial, smbtfial_nofl)
#endif
        do is = 0, ns-1
         if(isb.ne.is .and.
     &      .not.is_neutral(isb) .and. .not.is_neutral(is)) then

          do iy = 0, ny-1
           do ix = 0, nx-1
            if(leftix(ix,iy).ne.-2 .and. rightix(ix,iy).ne.nx+1) then
             Kabtf(ix,iy) = fkabtf(ix,iy,isb,is)
             Kbatf(ix,iy) = fkabtf(ix,iy,is,isb)
c            kb(ix,iy)    = fka(ix,iy,is)
c            kb(ix,iy)    = fka_new(ix,iy,is,rz2_temp,na_temp)
             kb(ix,iy)    =
     &                 fka_new_opt(ix,iy,is,rz2_temp,na_temp,am_sqrt)
             KabVp(ix,iy) = fkabvp(ix,iy,isb,is)
c        .. calculations could be simplified
c        .. but we did not do it in order to have the same formulae as manual
             calfab = zetap(ix,iy)/sqrt(mp)*
     &            qe*(na(ix,iy,isb)+na(ix,iy,is))*
     &            sqrt(am(isb)*am(is)/(am(isb)+am(is)))*
     &           (Kabtf(ix,iy)/kb(ix,iy)-Kbatf(ix,iy)/ka(ix,iy))*
     &            bb(ix,iy,0)/bb(ix,iy,3)
             coef = phm3*rz2(ix,iy,isb)*rz2(ix,iy,is)*calfab*
     &            mp/(zetap(ix,iy)*qe)
             delti=
     &         (gti(ix,iy)+gti(rightix(ix,iy),rightiy(ix,iy)))*0.5e0_R8
c   .... apply flux limit
             if(cflim(4).ne.0.0_R8.and.
     &         .not.on_closed_surface(ix,iy)) then
               jabmax = qe*KabVp(ix,iy)*sqrt(ti(ix,iy)/mp)*
     &          (na(ix,iy,isb)+na(ix,iy,is))
               t0 = abs(calfab*delti)/(cflim(4)*jabmax)
               fllim_al_ab = 1.0_R8/(1.0e0_R8+t0)
             else
               fllim_al_ab = 1.0_R8
             endif

             t0 = na(ix,iy,isb)*na(ix,iy,is)/    !   *vol(ix,iy)*hz(ix,iy) ! it will be multiplied in b2npmo !srv 05.07.17
     &                              ((na(ix,iy,isb)+na(ix,iy,is)))*delti
             smbtfial(ix,iy,0) = smbtfial(ix,iy,0) + fllim_al_ab*coef*t0
             smbtfial_nofl(ix,iy,0) = smbtfial_nofl(ix,iy,0) + coef*t0
            endif
           enddo
          enddo
         endif
        enddo
      endif
      wrk(-1:nx,-1:ny) = smbtfial(-1:nx,-1:ny,0)
      call nancheck(nx, ny, wrk, 'smbtfial(,,0)')
#ifdef WITH_FTIMINGS
      if (make_parallel) call timer%stop("parallel")
#endif
      smbtf = smbtfeal + smbtfial

      smbfria = smbfrial
      smbfrea = smbfreal
      smbtfia = smbtfial
      smbtfea = smbtfeal

      if (balance_netcdf.ne.0) then !djm Jun2019
        do iy=-1,ny
          do ix =-1,nx
            b2sifr_smofria0to3(ix,iy,0:3,isb)=vol(ix,iy)*hz(ix,iy)*
     &                                        smbfria(ix,iy,0:3)
            b2sifr_smofrea0to3(ix,iy,0:3,isb)=vol(ix,iy)*hz(ix,iy)*
     &                                        smbfrea(ix,iy,0:3)
            b2sifr_smotfia0to3(ix,iy,0:3,isb)=vol(ix,iy)*hz(ix,iy)*
     &                                        smbtfia(ix,iy,0:3)
            b2sifr_smotfea0to3(ix,iy,0:3,isb)=vol(ix,iy)*hz(ix,iy)*
     &                                        smbtfea(ix,iy,0:3)
          enddo
        enddo
      endif

      if (iout.ne.0) then !  .or. iout_b2wdat.eq.4) then                !srv 05.07.17
        write (chns,'(i3.3)') isb

        wrk = vol*hz*(smbfria(:,:,0)+smbfria(:,:,1)*ua(:,:,isb)+        !srv 05.07.17 {
     &        smbfria(:,:,2)*na(:,:,isb)*mp*am(isb)+
     &        smbfria(:,:,3)*ua(:,:,isb)*na(:,:,isb)*mp*am(isb))
        call my_out(70,nx,ny,wrk,'b2sifrtf_smofria'//chns)

        wrk = vol*hz*(smbfreal(:,:,0)+smbfreal(:,:,1)*ua(:,:,isb))
        call my_out(70,nx,ny,wrk,'b2sifrtf_smofrea'//chns)
        wrk = vol*hz*smbtfial(:,:,0)
        call my_out(70,nx,ny,wrk,'b2sifrtf_smotfia'//chns)
        wrk = vol*hz*smbtfeal(:,:,0)
        call my_out(70,nx,ny,wrk,'b2sifrtf_smotfea'//chns)              !srv 05.07.17 }
        call my_out(70,nx,ny,ue(:,:),'b2sifrtf_ue')
        if (Zhdanov_nc.eq.1) then
          call my_out(70,nx,ny,corr_fria(:,:,isb),'b2sifrtf_cf'//chns)
          call my_out(70,nx,ny,corr_tfia(:,:,isb),'b2sifrtf_ct'//chns)
          call my_out(70,nx,ny,alpha_HS_style(:,:,isb),
     &  'b2sifrtf_alpha'//chns)
          call my_out(70,nx,ny,g_HS_style(:,:),
     &  'b2sifrtf_g')
        endif

        wrk = vol*hz*(smbfrial_test(:,:,0)+
     &        smbfrial_test(:,:,1)*ua(:,:,isb)+
     &        smbfrial_test(:,:,2)*na(:,:,isb)*mp*am(isb)+
     &        smbfrial_test(:,:,3)*ua(:,:,isb)*na(:,:,isb)*mp*am(isb))
        call my_out(70,nx,ny,wrk,'b2sifrtf_smofria_test'//chns)
        wrk = vol*hz*smbtfial_test(:,:,0)
        call my_out(70,nx,ny,wrk,'b2sifrtf_smotfia_test'//chns)
        wrk = vol*hz*smbtfial_nofl(:,:,0)
        call my_out(70,nx,ny,wrk,'b2sifrtf_smotfia_nofl'//chns)
        wrk = vol*hz*smbtfial_test_nofl(:,:,0)
        call my_out(70,nx,ny,wrk,'b2sifrtf_smotfia_test_nofl'//chns)
      endif

      if (iout.gt.1) then                                               !srv 05.07.17
        do is = 0,ns-1
          write (chis,'(i3.3)') is
          do iy = -1, ny
            do ix = -1, nx
              wrk(ix,iy) = fkabvp(ix, iy, isb, is)
            end do
          end do
          call my_out(70,nx,ny,wrk,'b2sifrtf_kabvp_'//chns//'_'//chis)
          do iy = -1, ny
            do ix = -1, nx
              wrk(ix,iy) = fkabtf(ix, iy, isb, is)
            end do
          end do
          call my_out(70,nx,ny,wrk,'b2sifrtf_kabtf_'//chns//'_'//chis)
          do iy = -1, ny
            do ix = -1, nx
              wrk(ix,iy) = fkabtf(ix, iy, is, isb)
            end do
          end do
          call my_out(70,nx,ny,wrk,'b2sifrtf_kbatf_'//chns//'_'//chis)
          do iy = -1, ny
            do ix = -1, nx
              wrk(ix,iy) = fka(ix, iy, isb)
            end do
          end do
          call my_out(70,nx,ny,wrk,'b2sifrtf_ka_'//chns)
          do iy = -1, ny
            do ix = -1, nx
              wrk(ix,iy) = fka(ix, iy, is)
            end do
          end do
          call my_out(70,nx,ny,wrk,'b2sifrtf_kb_'//chis)
        enddo

        call my_out(70,nx,ny,f_luc_sg,'b2sifrtf_f_luc_sg')
        call my_out(70,nx,ny,ce1,'b2sifrtf_ce1')
        call my_out(70,nx,ny,ce2,'b2sifrtf_ce2')
        call my_out(70,nx,ny,cimp1(:,:,isb),'b2sifrtf_cimp1'//chns)     !som 13.11.20
        call my_out(70,nx,ny,cimp2(:,:,isb),'b2sifrtf_cimp2'//chns)     !som 13.11.20
        call my_out(70,nx,ny,alfx_c,'b2sifrtf_alfx_c')
        call my_out(70,nx,ny,sigx_c,'b2sifrtf_sigx_c')
        call my_out(70,nx,ny,zetae,'b2sifrtf_zetae')

        call my_out(70,nx,ny,ua(:,:,isb),'b2sifrtf_ua'//chns)
      endif


* ..return
      ncall = ncall+1
      call subend ()
      return

      contains

      function fkabvp(ix, iy, a, b)
      implicit none
      integer, intent(in) :: ix, iy, a, b
      real (kind=R8) :: fkabvp
      if((a.eq.ismain) .and. (b.ne.a) .and.
     &  (.not.is_neutral(a)).and.(.not.is_neutral(b))) then
        fkabvp=cimp1(ix,iy,b)                                 !som 13.11.20
      elseif((b.eq.ismain) .and. (a.ne.b) .and.
     &  (.not.is_neutral(a)).and.(.not.is_neutral(b))) then
        fkabvp=cimp1(ix,iy,a)                                 !som 13.11.20
      elseif((a.ne.b .and. a.ne.ismain) .and. (b.ne.ismain) .and.
     &  (.not.is_neutral(a)).and.(.not.is_neutral(b))) then
        fkabvp=1.0_R8
      elseif((a.eq.b) .and.
     &  (.not.is_neutral(a)).and.(.not.is_neutral(b))) then
        fkabvp=0.0_R8
      else
        fkabvp=0.0_R8
      endif
      return
      end function fkabvp

      function fkabtf(ix, iy, a, b)
      implicit none
      integer, intent(in) :: ix, iy, a, b
      real (kind=R8) :: fkabtf
      if((b.eq.ismain) .and. (b.ne.a) .and.
     &  (.not.is_neutral(a)).and.(.not.is_neutral(b))) then
        fkabtf=cimp2(ix,iy,a)                                 !som 13.11.20
      elseif((a.ne.ismain) .and. (b.ne.ismain) .and.
     &  (.not.is_neutral(a)).and.(.not.is_neutral(b))) then
        fkabtf=0.0_R8
      elseif((a.eq.ismain)) then
        fkabtf=0.0_R8
      else
        fkabtf=0.0_R8
      endif
      return
      end function fkabtf

      function fka(ix, iy, a)
      implicit none
      integer, intent(in) :: ix, iy, a
      real (kind=R8) :: fka
      integer r
      fka = 0.0_R8
#if ( defined(_OPENMP) && !defined(LEGACYCOMP) && !defined(NAGFOR) )
!$OMP SIMD REDUCTION(+:fka)
#endif
      do r = 0, ns-1
        fka = fka + rz2(ix,iy,r)*na(ix,iy,r)*sqrt(mp)*
     &          sqrt(am(a)*am(r)/(am(a)+am(r)))
      enddo
      fka = fka*rz2(ix,iy,a)
      return
      end function fka

      function fka_new(ix, iy, a, rz2_temp, na_temp)
      implicit none
      integer, intent(in) :: ix, iy, a
      real(kind=R8) :: fka_new
      real(kind=R8) :: rz2_temp(1:ns*(nx+2)*(ny+2))
      real(kind=R8) :: na_temp(1:ns*(nx+2)*(ny+2))
      integer r, ind
      fka_new = 0.0_R8
      ind = (ix+1)*ns+(iy+1)*ns*(nx+2)
#if ( defined(_OPENMP) && !defined(LEGACYCOMP) && !defined(NAGFOR) )
!$OMP SIMD REDUCTION(+:fka_new)
#endif
      do r = 0, ns-1
        fka_new = fka_new + rz2_temp(ind+r+1)*na_temp(ind+r+1)*
     &   sqrt(mp)*sqrt(am(a)*am(r)/(am(a)+am(r)))
      end do
      fka_new = fka_new*rz2_temp(ind+a+1)
      return
      end function fka_new

      function fka_new_opt(ix, iy, a, rz2_temp, na_temp, am_sqrt)
      implicit none
      integer, intent(in) :: ix, iy, a
      real(kind=R8) :: fka_new_opt
      real(kind=R8) :: rz2_temp(1:ns*(nx+2)*(ny+2))
      real(kind=R8) :: na_temp(1:ns*(nx+2)*(ny+2))
      real(kind=R8) :: am_sqrt(0:ns-1,0:ns-1)
      integer r, ind
      fka_new_opt = 0.0_R8
      ind = (ix+1)*ns+(iy+1)*ns*(nx+2)
#if ( defined(_OPENMP) && !defined(LEGACYCOMP) && !defined(NAGFOR) )
!$OMP SIMD REDUCTION(+:fka_new_opt)
#endif
      do r = 0, ns-1
        fka_new_opt = fka_new_opt + rz2_temp(ind+r+1)*na_temp(ind+r+1)*
     &   sqrt(mp)*am_sqrt(r,a)
      end do
      fka_new_opt = fka_new_opt*rz2_temp(ind+a+1)
      return
      end function fka_new_opt

*-----------------------------------------------------------------------
*.end b2sifrtf

      end subroutine b2sifrtf                                           !srv 01.07.05 }
