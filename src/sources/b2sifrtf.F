*-----------------------------------------------------------------------
*.specification

      subroutine b2sifrtf (nCv, nFc, nVx, ns, isb, ismain, switch, geo,
     &   mpg, na, ua, ue, te, ti, ni, ne, ne2, rza, rz2, po,
     &   zetap, zetae, gti, gte, ce1, ce2, cimp1, cimp2, zeff, 
     &   f_luc_sg, alfx_c, sigx_c, st_ext,
     &   smbfrea, smbfria, smbtfea, smbtfia, smbch, smbtf)
      use b2mod_types
      use b2mod_constants
      use b2mod_b2cmpa
      use b2mod_b2cmpt
      use b2mod_switches
      use b2us_geo
      use b2us_map
      use b2us_plasma
! csc The following are not necessary for computation but are needed
!     for adjoint AD to avoid side-effect variables
      use b2mod_ad
     & , only : ncall_b2sifrtf
      use b2mod_subsys
      implicit none
*   ..input arguments (unchanged on exit)
      integer nCv, nFc, nVx, ns, isb, ismain
      type (switches), intent (in) :: switch
      type (geometry), intent (in) :: geo
      type (mapping), intent (in) :: mpg
      type (B2StateExt), intent (in) :: st_ext
      real (kind=R8) ::
     &  na(nCv,0:ns-1), ua(nCv,0:ns-1), ue(nCv),
     &  te(nCv), ti(nCv), ni(nCv,0:1), ne(nCv), ne2(nCv),
     &  rza(nCv,0:ns-1), rz2(nCv,0:ns-1), po(nCv),
     &  zetap(nCv), zetae(nCv), gti(nCv), gte(nCv),
     &  ce1(nCv), ce2(nCv), cimp1(nCv), cimp2(nCv), zeff(nCv),          !srv 13.01.17
     &  f_luc_sg(nFc), alfx_c(nCv), sigx_c(nCv)
*   ..output arguments (unspecified on entry)
      real (kind=R8) ::
     *  smbfria(nCv,0:3), smbfrea(nCv,0:3),
     *  smbtfia(nCv,0:3), smbtfea(nCv,0:3)                              !srv 13.01.17 05.07.17
      real (kind=R8) ::
     *  smbch(nCv,0:3), smbtf(nCv,0:3)
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
      integer iCv, is
      integer, save :: iout = 0, iout_b2npmo = 0, iout_b2wdat =0
      character chns*3, chis*3
      real (kind=R8) ::
     *  coef, t0, calfab, fllim_al_ab, jabmax
      real (kind=R8) ::
     &  KabVp(nCv), Kabtf(nCv), Kbatf(nCv), ka(nCv), kb(nCv), wrk(nCv)  !srv 13.01.17 {
! Local arrays smb*al have been introduced specially.
! The usage of smb*a instaed of smb*al should be checked.
      real (kind=R8) ::
     *  smbfrial(nCv,0:3), smbfreal(nCv,0:3),
     *  smbtfial(nCv,0:3), smbtfeal(nCv,0:3)                            !srv 13.01.17 }
*   ..procedures
      intrinsic sqrt
      external xertst, ipgeti, ipgetr, sfill
      external b2xvff, b2xvfx, b2xvsg

*-----------------------------------------------------------------------
*.computation

* ..preliminaries
*   ..subprogram start-up calls
      call subini ('b2sifrtf')
*   ..set internal parameters on first call
      if (ncall_b2sifrtf.eq.0) then
       call xertst (0.le.switch%b2sifr_styl0.and.
     &              switch%b2sifr_styl0.le.2,
     &  'faulty internal parameter styl0')
       call xertst (0.0_R8.le.switch%b2sifr_phm0.and.
     &              0.0_R8.le.switch%b2sifr_phm1.and.
     &              0.0_R8.le.switch%b2sifr_phm2.and.
     &              0.0_R8.le.switch%b2sifr_phm3,
     &  'faulty internal parameter phm0, phm1, phm2, phm3')
       call ipgeti ('b2sifrtf_iout', iout)
       call ipgeti ('b2wdat_iout', iout_b2wdat)                         !iys 17.04.17
       call ipgeti ('b2npmo_iout', iout_b2npmo)
      endif
*   ..test nCv, nFc
      call xertst (0.le.nCv.and.0.le.nFc, 'faulty argument nCv, nFc')
*   ..extensive tests on first few calls
      if (ncall_b2sifrtf.lt.3) then
*    ..test sign of vol
       call b2xvsg (nCv, geo%cvVol, 1, 'vol', '.gt.')
*    ..test sign of na, ni, ti, ne, ne2
       call b2xvsg (nCv*ns, na, 1, 'na', '.gt.')
       call b2xvsg (nCv*2, ni, 1, 'ni', '.gt.')
       call b2xvsg (nCv, ti, 1, 'ti', '.gt.')
       call b2xvsg (nCv, ne, 1, 'ne', '.gt.')
       call b2xvsg (nCv, ne2, 1, 'ne2', '.gt.')
      endif

*   ..initialize to 0
      call sfill (nCv*4, 0.0e0_R8, smbch, 1)                            !srv 11.09.09
      call sfill (nCv*4, 0.0e0_R8, smbtf, 1)                            !srv 11.09.09
      call sfill (nCv*4, 0.0e0_R8, smbfrial, 1)                         !srv 13.01.17 { !srv 05.07.17
      call sfill (nCv*4, 0.0e0_R8, smbfreal, 1)
      call sfill (nCv*4, 0.0e0_R8, smbtfial, 1)
      call sfill (nCv*4, 0.0e0_R8, smbtfeal, 1)                         !srv 13.01.17 }
* ..compute friction
*   ..loop over all other species
      do is = 0, ns-1
       if(isb.ne.is .and.
     &    .not.is_neutral(isb).and. .not.is_neutral(is)) then
*     ..contribution from friction between ions
        do iCv = 1, mpg%nCi
         KabVp(iCv)=fkabvp(iCv,isb,is)
         coef = switch%b2sifr_phm0*rz2(iCv,isb)*rz2(iCv,is)*
     &          mp*sqrt(am(isb)*am(is)/(am(isb)+am(is)))/
     &          zetap(iCv)*KabVp(iCv)
         t0 = na(iCv,isb)*na(iCv,is)
         if (switch%b2sifr_styl0.eq.0) then
          smbfrial(iCv,0) = smbfrial(iCv,0)+coef*t0*ua(iCv,is)
          smbfrial(iCv,1) = smbfrial(iCv,1)-coef*t0
         else if (switch%b2sifr_styl0.eq.1) then
          smbfrial(iCv,0) = smbfrial(iCv,0)+coef*t0*ua(iCv,is)
          smbfrial(iCv,3) = smbfrial(iCv,3)-
     &                      coef*t0/roxa(iCv,isb)
         else if (switch%b2sifr_styl0.eq.2) then
          smbfrial(iCv,2) = smbfrial(iCv,2)+coef*t0*ua(iCv,is)/
     &                      roxa(iCv,isb)
          smbfrial(iCv,3) = smbfrial(iCv,3)-coef*t0/
     &                      roxa(iCv,isb)
         endif
        enddo
       endif
      enddo

!*   ..compute ue, electron parallel velocity
!      call b2xpve (nCv, nFc, ns, switch, geo, mpg, qe, rza, rz2,
!     & alfx_c, sigx_c, fch_pi_c, po, ne, te, na, ua, st_ext, ue)
*     ..contribution from friction between between electron and ions
      if(.not.is_neutral(isb)) then
       do iCv = 1, mpg%nCi
        coef = switch%b2sifr_phm1*(geo%cvBb(iCv,0)/geo%cvBb(iCv,3))**2*
     &         rz2(iCv,isb)*qe**2/(sigx_c(iCv)*zeff(iCv))
        t0 = na(iCv,isb)*ne(iCv)
        smbfreal(iCv,0) = smbfreal(iCv,0)+coef*t0*ue(iCv)
        smbfreal(iCv,1) = smbfreal(iCv,1)-coef*t0
       enddo
      endif

      smbch = smbfreal + smbfrial

* ..compute thermal force

*   ..contribution from electron thermal force
      if(.not.is_neutral(isb)) then
       !wdk interpolate f_luc_sg to cell center
       call intcell(nFc,nCv,mpg,mpg%intcellP,f_luc_sg,wrk)
       do iCv = 1, mpg%nCi
        coef = -switch%b2sifr_phm2*alfx_c(iCv)*ce1(iCv)*rz2(iCv,isb)*
     &    (me/qe)/zetae(iCv)/(geo%cvBb(iCv,0)/geo%cvBb(iCv,3))/wrk(iCv)
        t0 = na(iCv,isb)*gte(iCv)
        smbtfeal(iCv,0) = smbtfeal(iCv,0)+coef*t0
       enddo
      endif
*   ..contribution from ions thermal force
      do is = 0, ns-1
       if(isb.ne.is .and.
     &    .not.is_neutral(isb) .and. .not.is_neutral(is)) then
        do iCv = 1, mpg%nCi
          Kabtf(iCv) = fkabtf(iCv,isb,is)
          Kbatf(iCv) = fkabtf(iCv,is,isb)
          ka(iCv)    = fka(iCv,isb)
          kb(iCv)    = fka(iCv,is)
          KabVp(iCv) = fkabvp(iCv,isb,is)
c      .. calculations could be simplified
c      .. but we did not do it in order to have the same formulae as manual
          calfab = zetap(iCv)/sqrt(mp)*
     &         qe*(na(iCv,isb)+na(iCv,is))*
     &         sqrt(am(isb)*am(is)/(am(isb)+am(is)))*
     &        (Kabtf(iCv)/kb(iCv)-Kbatf(iCv)/ka(iCv))*
     &         geo%cvBb(iCv,0)/geo%cvBb(iCv,3)
          coef = switch%b2sifr_phm3*rz2(iCv,isb)*rz2(iCv,is)*calfab*
     &         mp/(zetap(iCv)*qe)
c      .. apply flux limit
          if(cflim(4).ne.0.0_R8.and..not.mpg%cvOnClosedSurface(iCv))then
            jabmax = qe*KabVp(iCv)*sqrt(ti(iCv)/mp)*
     &       (na(iCv,isb)+na(iCv,is))
            t0 = abs(calfab*gti(iCv))/(cflim(4)*jabmax)
            fllim_al_ab = 1.0_R8/(1.0e0_R8+t0)
          else
            fllim_al_ab = 1.0_R8
          endif

          t0 = na(iCv,isb)*na(iCv,is)/
     &                            ((na(iCv,isb)+na(iCv,is)))*gti(iCv)
          smbtfial(iCv,0) = smbtfial(iCv,0) + fllim_al_ab*coef*t0
        enddo
       endif
      enddo

      smbtf = smbtfeal + smbtfial

      smbfria = smbfrial
      smbfrea = smbfreal
      smbtfia = smbtfial
      smbtfea = smbtfeal

      if (iout.ne.0) then !  . or. iout_b2wdat.eq.4) then              !srv 05.07.17
       write (chns,'(i3.3)') isb

       wrk = geo%cvVol*geo%cvHz*
     &       (smbfria(:,0)+smbfria(:,1)*ua(:,isb)+        !srv 05.07.17 {
     &        smbfria(:,2)*na(:,isb)*mp*am(isb)+
     &        smbfria(:,3)*ua(:,isb)*na(:,isb)*mp*am(isb))
       call my_out_us(70,nCv,0,wrk,'b2sifrtf_smofria'//chns)
       wrk = geo%cvVol*geo%cvHz*(smbfreal(:,0)+smbfreal(:,1)*ua(:,isb))
       call my_out_us(70,nCv,0,wrk,'b2sifrtf_smofrea'//chns)
       wrk = geo%cvVol*geo%cvHz*smbtfial(:,0)
       call my_out_us(70,nCv,0,wrk,'b2sifrtf_smotfia'//chns)
       wrk = geo%cvVol*geo%cvHz*smbtfeal(:,0)
       call my_out_us(70,nCv,0,wrk,'b2sifrtf_smotfea'//chns)            !srv 05.07.17 }
       call my_out_us(70,nCv,0,ue,'b2sifrtf_ue')
      endif

      if (iout.gt.1) then                                               !srv 05.07.17
        do is = 0,ns-1
          write (chis,'(i3.3)') is
          do iCv = 1, nCv
            wrk(iCv) = fkabvp(iCv, isb, is)
          end do
          call my_out_us(70,nCv,0,wrk,
     &              'b2sifrtf_kabvp_'//chns//'_'//chis)
          do iCv = 1, nCv
            wrk(iCv) = fkabtf(iCv, isb, is)
          end do
          call my_out_us(70,nCv,0,wrk,
     &              'b2sifrtf_kabtf_'//chns//'_'//chis)
          do iCv = 1, nCv
            wrk(iCv) = fkabtf(iCv, is, isb)
          end do
          call my_out_us(70,nCv,0,wrk,
     &              'b2sifrtf_kbatf_'//chns//'_'//chis)
          do iCv = 1, nCv
            wrk(iCv) = fka(iCv, isb)
          end do
          call my_out_us(70,nCv,0,wrk,'b2sifrtf_ka_'//chns)
          do iCv = 1, nCv
            wrk(iCv) = fka(iCv, is)
          end do
          call my_out_us(70,nCv,0,wrk,'b2sifrtf_kb_'//chis)
        enddo

        call intcell(nFc,nCv,mpg,mpg%intcellP,f_luc_sg,wrk)
        call my_out_us(70,nCv,0,wrk,'b2sifrtf_f_luc_sg')
        call my_out_us(70,nCv,0,ce1,'b2sifrtf_ce1')
        call my_out_us(70,nCv,0,ce2,'b2sifrtf_ce2')
        call my_out_us(70,nCv,0,cimp1,'b2sifrtf_cimp1')
        call my_out_us(70,nCv,0,cimp2,'b2sifrtf_cimp2')
        call my_out_us(70,nCv,0,alfx_c,'b2sifrtf_alfx_c')
        call my_out_us(70,nCv,0,sigx_c,'b2sifrtf_sigx_c')
        call my_out_us(70,nCv,0,zetae,'b2sifrtf_zetae')

        call my_out_us(70,nCv,0,ua(:,isb),'b2sifrtf_ua'//chns)
      endif


* ..return
      ncall_b2sifrtf = ncall_b2sifrtf+1
      call subend ()
      return

      contains

      function fkabvp(iCv, a, b)
      implicit none
      integer, intent(in) :: iCv, a, b
      real (kind=R8) :: fkabvp
      if((a.eq.ismain) .and. (b.ne.a) .and.
     &  (.not.is_neutral(a)).and.(.not.is_neutral(b))) then
        fkabvp(iCv)=cimp1(iCv)
      elseif((b.eq.ismain) .and. (a.ne.b) .and.
     &  (.not.is_neutral(a)).and.(.not.is_neutral(b))) then
        fkabvp(iCv)=cimp1(iCv)
      elseif((a.ne.b .and. a.ne.ismain) .and. (b.ne.ismain) .and.
     &  (.not.is_neutral(a)).and.(.not.is_neutral(b))) then
        fkabvp(iCv)=1.0_R8
      elseif((a.eq.b) .and.
     &  (.not.is_neutral(a)).and.(.not.is_neutral(b))) then
        fkabvp(iCv)=0.0_R8
      else
        fkabvp(iCv)=0.0_R8
      endif
      return
      end function fkabvp

      function fkabtf(iCv, a, b)
      implicit none
      integer, intent(in) :: iCv, a, b
      real (kind=R8) :: fkabtf
      if((b.eq.ismain) .and. (b.ne.a) .and.
     &  (.not.is_neutral(a)).and.(.not.is_neutral(b))) then
        fkabtf(iCv)=cimp2(iCv)
      elseif((a.ne.ismain) .and. (b.ne.ismain) .and.
     &  (.not.is_neutral(a)).and.(.not.is_neutral(b))) then
        fkabtf(iCv)=0.0_R8
      elseif((a.eq.ismain)) then
        fkabtf(iCv)=0.0_R8
      else
        fkabtf(iCv)=0.0_R8
      endif
      return
      end function fkabtf

      function fka(iCv, a)
      implicit none
      integer, intent(in) :: iCv, a
      real (kind=R8) :: fka
      integer r
      fka(iCv) = 0.0_R8
      do r = 0, ns-1
        fka(iCv) = fka(iCv) + rz2(iCv,r)*na(iCv,r)*sqrt(mp)*
     &          sqrt(am(a)*am(r)/(am(a)+am(r)))
      enddo
      fka(iCv) = fka(iCv)*rz2(iCv,a)
      return
      end function fka

      function roxa(iCv, is)
      implicit none
      integer, intent(in) :: iCv, is
      real (kind=R8) :: roxa

      roxa = am(is)*mp*na(iCv,is)
      return
      end function roxa

*-----------------------------------------------------------------------
*.end b2sifrtf

      end subroutine b2sifrtf                                           !srv 01.07.05 }
