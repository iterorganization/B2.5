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
      use b2mod_plasma
     , , only : roxa
      use b2mod_transport
     , , only : f_luc_sg, alfx_c, sigx_c                                !srv 08.10.14
      use b2mod_external
      use b2mod_geo
     , , only : hz, vol                                                 !srv 05.07.17
      use b2mod_openmp
      use b2mod_subsys
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
     &  cimp1(-1:nx,-1:ny), cimp2(-1:nx,-1:ny),
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
#ifdef NO_OPENMP_B2NPMO
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
     *  smbtfial(-1:nx,-1:ny,0:3), smbtfeal(-1:nx,-1:ny,0:3)            !srv 13.01.17 }
#ifdef WITH_FTIMINGS
      logical :: do_parallel
#endif
*   ..procedures
      intrinsic sqrt
      external xertst, ipgeti, ipgetr, sfill
      external b2xvff, b2xvfx, b2xvsg
#ifdef NO_IEEE
      logical, external :: my_isnan
#endif
*   ..initialization
      save ncall, styl0, phm0, phm1, phm2, phm3
#ifdef NO_OPENMP_B2NPMO
      !$OMP THREADPRIVATE(ncall, styl0, phm0, phm1, phm2, phm3)
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
          !$OMP MASTER
          write(*,*) 'B2SIFRTF OpenMP warning: no file output in ',
     &               'parallel mode'
          !$OMP END MASTER
          iout = 0
        end if
       call ipgeti ('b2wdat_iout', iout_b2wdat)                         !iys 17.04.17
       if ((iout_b2wdat==4) .and. in_parallel()) then
          !$OMP MASTER
          write(*,*) 'B2SIFRTF OpenMP warning: no file output in ',
     &               'parallel mode (b2wdat)'
          !$OMP END MASTER
          iout_b2wdat = 0
        end if
       call ipgeti ('b2npmo_iout', iout_b2npmo)
      endif
#ifdef WITH_FTIMINGS
      do_parallel = .not. in_parallel()
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
      call sfill (n2*4, 0.0e0_R8, smbtfeal, 1)                          !srv 13.01.17 }
* ..compute friction
*   ..loop over all other species
#ifdef WITH_FTIMINGS
      if (do_parallel) call timer%start("parallel")
#endif
#ifndef NO_OPENMP_B2SIFRTF
!$OMP PARALLEL DO IF (.not. OMP_IN_PARALLEL())
!$OMP& DEFAULT(NONE)
!$OMP& PRIVATE(KabVp, is, ix, iy, coef, t0)
!$OMP& SHARED(ns, isb, is_neutral, nx, ny, leftix, rightix)
!$OMP& SHARED(rz2, am, zetap, na, ua)
#ifndef NO_OPENMP_B2NPMO
!$OMP& PRIVATE(styl0, phm0)
#else
!$OMP& COPYIN(styl0, phm0)
#endif
!$OMP& REDUCTION(+: smbfrial)
#endif
      do is = 0, ns-1
       if(isb.ne.is .and. 
     &    .not.is_neutral(isb).and. .not.is_neutral(is)) then
        KabVp=fkabvp(isb,is)
*     ..contribution from friction between ions
        do iy = 0, ny-1
         do ix = 0, nx-1
          if(leftix(ix,iy).ne.-2 .and. rightix(ix,iy).ne.nx+1) then     !srv 26.11.02
           coef = phm0*rz2(ix,iy,isb)*rz2(ix,iy,is)*
     &           mp*sqrt(am(isb)*am(is)/(am(isb)+am(is)))/
     &           zetap(ix,iy)*KabVp(ix,iy)
           t0 = na(ix,iy,isb)*na(ix,iy,is)   !   *vol(ix,iy)*hz(ix,iy) ! it will be multiplied in b2npmo !srv 05.07.17
             !$omp critical (ionfriction)
           if (styl0.eq.0) then
             smbfrial(ix,iy,0) = smbfrial(ix,iy,0)+coef*t0*ua(ix,iy,is)
             smbfrial(ix,iy,1) = smbfrial(ix,iy,1)-coef*t0
           else if (styl0.eq.1) then
             smbfrial(ix,iy,0) = smbfrial(ix,iy,0)+coef*t0*ua(ix,iy,is)
             smbfrial(ix,iy,3) = smbfrial(ix,iy,3)-
     &                           coef*t0/roxa(ix,iy,isb)
           else if (styl0.eq.2) then
             smbfrial(ix,iy,2) = smbfrial(ix,iy,2)+coef*t0*ua(ix,iy,is)/
     &                       roxa(ix,iy,isb)
             smbfrial(ix,iy,3) = smbfrial(ix,iy,3)-coef*t0/
     &                       roxa(ix,iy,isb)
           endif
             !$omp end critical (ionfriction)
          endif
         enddo
        enddo
       endif
      enddo
#ifndef NO_OPENMP_B2SIFRTF
!$OMP END PARALLEL DO
#endif
      call nancheck(nx, ny, smbfrial(:,:,0), 'smbfrial(,,0)')
#ifdef WITH_FTIMINGS
      if (do_parallel) call timer%stop("parallel")
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
      if (do_parallel) call timer%start("parallel")
#endif
#ifndef NO_OPENMP_B2SIFRTF
!$OMP PARALLEL DO IF (.not.OMP_IN_PARALLEL())
!$OMP& DEFAULT(NONE)
!$OMP& PRIVATE(KabVp, Kbatf, Kabtf, ka, kb, is, ix, iy, coef, t0)
!$OMP& PRIVATE(calfab, delti, jabmax, fllim_al_ab)
!$OMP& SHARED(ns, isb, is_neutral, nx, ny, leftix, rightix)
!$OMP& SHARED(rz2, am, zetap, na, ua, bb, gti, cflim, ti, rightiy)
!$OMP& SHARED(on_closed_surface)
#ifdef NO_OPENMP_B2NPMO
!$OMP& PRIVATE(phm3)
#else
!$OMP& SHARED(phm3)
#endif
!$OMP& REDUCTION(+: smbtfial)
#endif
      do is = 0, ns-1
       if(isb.ne.is .and.
     &    .not.is_neutral(isb) .and. .not.is_neutral(is)) then
        Kabtf = fkabtf(isb,is)
        Kbatf = fkabtf(is,isb)
        ka    = fka(isb)
        kb    = fka(is)
        KabVp = fkabvp(isb,is)

        do iy = 0, ny-1
         do ix = 0, nx-1
          if(leftix(ix,iy).ne.-2 .and. rightix(ix,iy).ne.nx+1) then
c      .. calculations could be simplified
c      .. but we did not do it in order to have the same formulae as manual
           calfab = zetap(ix,iy)/sqrt(mp)*
     &          qe*(na(ix,iy,isb)+na(ix,iy,is))*
     &          sqrt(am(isb)*am(is)/(am(isb)+am(is)))*
     &         (Kabtf(ix,iy)/kb(ix,iy)-Kbatf(ix,iy)/ka(ix,iy))*
     &          bb(ix,iy,0)/bb(ix,iy,3)
           coef = phm3*rz2(ix,iy,isb)*rz2(ix,iy,is)*calfab*
     &          mp/(zetap(ix,iy)*qe)
           delti=
     &       (gti(ix,iy)+gti(rightix(ix,iy),rightiy(ix,iy)))*0.5e0_R8
c   .. apply flux limit
           if(cflim(4).ne.0.0_R8.and..not.on_closed_surface(ix,iy)) then
             jabmax = qe*KabVp(ix,iy)*sqrt(ti(ix,iy)/mp)*
     &        (na(ix,iy,isb)+na(ix,iy,is))
             t0 = abs(calfab*delti)/(cflim(4)*jabmax)
             fllim_al_ab = 1.0_R8/(1.0e0_R8+t0)
           else
             fllim_al_ab = 1.0_R8
           endif

           t0 = na(ix,iy,isb)*na(ix,iy,is)/    !   *vol(ix,iy)*hz(ix,iy) ! it will be multiplied in b2npmo !srv 05.07.17
     &                              ((na(ix,iy,isb)+na(ix,iy,is)))*delti
           smbtfial(ix,iy,0) = smbtfial(ix,iy,0) + fllim_al_ab*coef*t0
          endif
         enddo
        enddo
       endif
      enddo
#ifndef NO_OPENMP_B2SIFRTF
!$OMP END PARALLEL DO
#endif
      call nancheck(nx, ny, smbtfial(:,:,0), 'smbtfial(,,0)')
#ifdef WITH_FTIMINGS
      if (do_parallel) call timer%stop("parallel")
#endif
      smbtf = smbtfeal + smbtfial

      smbfria = smbfrial
      smbfrea = smbfreal
      smbtfia = smbtfial
      smbtfea = smbtfeal

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
      endif

      if (iout.gt.1) then                                               !srv 05.07.17
        do is = 0,ns-1
          write (chis,'(i3.3)') is
          wrk = fkabvp(isb, is)
          call my_out(70,nx,ny,wrk,'b2sifrtf_kabvp_'//chns//'_'//chis)
          wrk = fkabtf(isb, is)
          call my_out(70,nx,ny,wrk,'b2sifrtf_kabtf_'//chns//'_'//chis)
          wrk = fkabtf(is, isb)
          call my_out(70,nx,ny,wrk,'b2sifrtf_kbatf_'//chns//'_'//chis)
          wrk = fka(isb)
          call my_out(70,nx,ny,wrk,'b2sifrtf_ka_'//chns)
          wrk = fka(is)
          call my_out(70,nx,ny,wrk,'b2sifrtf_kb_'//chis)
        enddo

        call my_out(70,nx,ny,f_luc_sg,'b2sifrtf_f_luc_sg')
        call my_out(70,nx,ny,ce1,'b2sifrtf_ce1')
        call my_out(70,nx,ny,ce2,'b2sifrtf_ce2')
        call my_out(70,nx,ny,cimp1,'b2sifrtf_cimp1')
        call my_out(70,nx,ny,cimp2,'b2sifrtf_cimp2')
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

      function fkabvp(a, b)
      implicit none
      integer a, b
      real (kind=R8) :: fkabvp(-1:nx,-1:ny)
      if((a.eq.ismain) .and. (b.ne.a) .and.
     &  (.not.is_neutral(a)).and.(.not.is_neutral(b))) then
        fkabvp=cimp1
      elseif((b.eq.ismain) .and. (a.ne.b) .and.
     &  (.not.is_neutral(a)).and.(.not.is_neutral(b))) then
        fkabvp=cimp1
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

      function fkabtf(a, b)
      implicit none
      integer a, b
      real (kind=R8) :: fkabtf(-1:nx,-1:ny)
      if((b.eq.ismain) .and. (b.ne.a) .and.
     &  (.not.is_neutral(a)).and.(.not.is_neutral(b))) then
        fkabtf=cimp2
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

      function fka(a)
      implicit none
      integer a
      real (kind=R8) :: fka(-1:nx,-1:ny)
      integer r
      fka = 0.0_R8
      do r = 0, ns-1
        fka = fka + rz2(:,:,r)*na(:,:,r)*sqrt(mp)*
     &          sqrt(am(a)*am(r)/(am(a)+am(r)))
      enddo
      fka = fka*rz2(:,:,a)
      return
      end function fka

*-----------------------------------------------------------------------
*.end b2sifrtf

      end subroutine b2sifrtf                                           !srv 01.07.05 }
