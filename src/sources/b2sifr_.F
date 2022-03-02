*-----------------------------------------------------------------------
*.specification

      subroutine b2sifr_ (nCv, nFc, nVx, ns, isb, switch, geo, mpg,
     &   na, ua, te, ti, ni, ne, ne2, rza, rz2, ehxp, lnlam, cthe, cthi,
     &   st_ext, smbch, smbtf, smprb, smptb, smfrb, vti32, gti, gte)
      use b2mod_types
!WG_RM      use b2mod_rates
!WG_RM      use b2mod_indirect
      use b2mod_constants
!WG_RM      use b2mod_external
      use b2mod_b2cmpa
!WG_TODO      use b2mod_balance                                               !djm Jan2017
!WG_TODO     & , only : b2sifr_smoch0to3, b2sifr_smotf0to3, balance_netcdf
      use b2mod_switches
      use b2us_geo
      use b2us_map
      use b2us_plasma
! csc The following are not necessary for computation but are needed
!     for adjoint AD to avoid side-effect variables
      use b2mod_ad
     & , only : ncall_b2sifr, ncall_b2tlnl
      implicit none
*   ..input arguments (unchanged on exit)
      integer nCv, nFc, nVx, ns, isb
      type (switches), intent (in) :: switch
      type (geometry), intent (in) :: geo
      type (mapping), intent (in) :: mpg
      type (B2StateExt), intent (in) :: st_ext
      real (kind=R8) ::
     &  na(nCv,0:ns-1), ua(nCv,0:ns-1), te(nCv), ti(nCv), ni(nCv,0:1),
     &  ne(nCv), ne2(nCv), rza(nCv,0:ns-1), rz2(nCv,0:ns-1), ehxp(nCv),
     &  lnlam(nCv), cthe(nCv,0:ns-1), cthi(nCv,0:ns-1)
*   ..output arguments (unspecified on entry)
      real (kind=R8) ::
     *   smbch(nCv,0:3), smbtf(nCv,0:3),                                !srv 11.09.09
     &   smprb(nCv), smptb(nCv), smfrb(nCv)
*   ..workspace arguments (unspecified on entry and on exit)
      real (kind=R8) ::
     *   vti32(nCv), gti(nCv), gte(nCv)
*   ..common blocks
*-----------------------------------------------------------------------
*.documentation
*
*  1. purpose
*
*     B2SIFR computes a linearised representation of the friction
*     force between one species and all others.
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

*   ..global named constants
      
*   ..local variables
      integer iCv, iFc, is, k                                           !srv 17.06.02
      integer, save :: iout = 0                                         !srv 17.06.02
      integer, save :: icase = 0
      character chns*3, chk*1                                           !srv 17.06.02
      real (kind=R8) ::
     *   ctaup, coef, t0, t1,
     *   wrkv(nVx)
      real (kind=R8) ::
     1 delti, delte, zlamii, zlamee, dtiapr, dteapr
*   ..procedures
      intrinsic sqrt
      real (kind=R8) ::
     *   roxa
      external subini, subend, xertst, xerrab, ipgeti, ipgetr, sfill
      external b2xvff, b2xvfx, b2xvsg
*   ..initialisation
      roxa(iCv,is) = am(is)*mp*na(iCv,is)                               !srv 01.07.05
*-----------------------------------------------------------------------
*.computation

* ..preliminaries
*   ..subprogram start-up calls
      call subini ('b2sifr_')
*   ..set internal parameters on first call
      if (ncall_b2sifr.eq.0) then
csw 24sep2013
       call xertst (
     &   0.le.switch%b2sifr_styl0.and.switch%b2sifr_styl0.lt.9, 
     &   'faulty internal parameter styl0')
       call xertst (
     &  0.0_R8.le.switch%b2sifr_phm0.and.
     &  0.0_R8.le.switch%b2sifr_phm1.and.
     &  0.0_R8.le.switch%b2sifr_phm2.and.
     &  0.0_R8.le.switch%b2sifr_phm3,
     &  'faulty internal parameter phm0, phm1, phm2, phm3')
       call xertst (0.0_R8.le.switch%b2sifr_limthii.and.
     &              0.0_R8.le.switch%b2sifr_limthee,
     &  'faulty internal parameter limthii, limthee')
       call ipgeti ('b2sifr_iout', iout)                             !srv 17.06.02
       call ipgeti ('b2tlnl_ii', icase)
       call xertst (0.eq.icase .or. 1.eq.icase, 
     &               'faulty input b2tlnl_ii')
      endif
*   ..test nCv, nFc
      call xertst (0.le.nCv.and.0.le.nFc, 'faulty argument nCv, nFc')
*   ..extensive tests on first few calls
      if (ncall_b2sifr.lt.3) then
*    ..test sign of vol
       call b2xvsg (nCv, geo%cvVol, 1, 'vol', '.gt.')
*    ..test sign of na, ni, ti, ne, ne2
       call b2xvsg (nCv*ns, na, 1, 'na', '.gt.')
       call b2xvsg (nCv*2, ni, 1, 'ni', '.gt.')
       call b2xvsg (nCv, ti, 1, 'ti', '.gt.')
       call b2xvsg (nCv, ne, 1, 'ne', '.gt.')
       call b2xvsg (nCv, ne2, 1, 'ne2', '.gt.')
      endif
*   ..compute ctaup
*     Following Balescu; Braginskii has ctaup different by a factor
*     sqrt(2). The mathematical expression is:
*       t0 = (3/(4*sqrt(2*pi)))*(4*pi*eps0/qe**2)**2*ev*sqrt(ev)
*       ctaup = t0*sqrt(mp)/lnlam
*     The numerical value in SI units is: ctaup = 1.4743e13/lnlam.
*     The following calculation of t0 is designed to avoid overflow
*     even for IEEE single precision arithmetic.
      t1 = (4*pi*eps0/qe)/qe
      t0 = (3/(4*sqrt(2*pi)))*t1*(t1*ev*sqrt(ev))
      ctaup = t0*sqrt(mp) ! /lnlam                !srv 01.07.09 it is taken below
*   ..compute the Coulomb logarithm
      call b2tlnl (nCv, switch, icase, te, ti, ne, lnlam)
*   ..compute vti32
      vti32 = geo%cvVol/((ti/ev)*sqrt(ti/ev))*geo%cvHz
*   ..compute gti,gte
      call gradc_p(nCv, nFc, nVx, 0, geo, mpg, ti, wrkv, gti)
      call gradc_p(nCv, nFc, nVx, 0, geo, mpg, te, wrkv, gte)
* ..compute friction
*   ..initialise to 0
      smbch = 0.0_R8
      smbtf = 0.0_R8
      smprb = 0.0_R8
      smptb = 0.0_R8
      smfrb = 0.0_R8
*   ..loop over all other species
      do is = 0, ns-1
*    ..consider distinct cases
       if (is.ne.isb.and.
     &    (.not.is_neutral(is)).and.(.not.is_neutral(isb))) then
*     ..contribution from friction between charged species
        do iCv = 1, mpg%nCi
          coef = switch%b2sifr_phm0*(lnlam(iCv)/ctaup)*rz2(iCv,isb)*    !srv 01.07.09
     &           rz2(iCv,is)*mp*sqrt(am(isb)*am(is)/(am(isb)+am(is)))
          t0 = na(iCv,isb)*na(iCv,is)*vti32(iCv)
          if (switch%b2sifr_styl0.eq.0) then
            smbch(iCv,0) = smbch(iCv,0)+coef*t0*ua(iCv,is)
            smbch(iCv,1) = smbch(iCv,1)-coef*t0
          else if (switch%b2sifr_styl0.eq.1) then
            smbch(iCv,0) = smbch(iCv,0)+coef*t0*ua(iCv,is)
            smbch(iCv,3) = smbch(iCv,3)-coef*t0/roxa(iCv,isb)
          else if (switch%b2sifr_styl0.eq.2) then
            smbch(iCv,2) = smbch(iCv,2)+coef*t0*ua(iCv,is)/
     &                     roxa(iCv,isb)
            smbch(iCv,3) = smbch(iCv,3)-coef*t0/
     &                     roxa(iCv,isb)
          endif
          smfrb(iCv) = smfrb(iCv) + 
     &      coef*t0*(ua(iCv,is)-ua(iCv,isb))
        enddo
       else
*     ..contribution from friction involving neutral species
*       (We ignore such a term, except for charge exchange friction
*       between H0 and H1. That term is calculated via the atomic
*       physics routines.)
       endif
      enddo
*   ..compute thermal force
      if (.not.is_neutral(isb)) then
       do iCv = 1, mpg%nCi
         t0 = na(iCv,isb)*
     &      (rza(iCv,isb)/(ne(iCv)-st_ext%ne(iCv))-
     &       rz2(iCv,isb)/(ne2(iCv)-st_ext%ne2(iCv)))
         delte=gte(iCv)
         if(switch%b2sifr_limthee.gt.0.0_R8) then
           zlamee=1.5e16_R8*(te(iCv)/ev)**2/ne(iCv)
           dteapr=switch%b2sifr_limthee/zlamee*te(iCv)/geo%cvBb(iCv,0)*
     &      geo%cvBb(iCv,3) !srv 12.04.10 11.01.13
           delte=sign(min(abs(delte),dteapr),delte)
         endif
         delti=gti(iCv)
         if(switch%b2sifr_limthii.gt.0.0_R8) then
           zlamii=1.5e16_R8*(ti(iCv)/ev)**2/ne2(iCv)
           dtiapr=switch%b2sifr_limthii/zlamii*ti(iCv)/geo%cvBb(iCv,0)*
     &      geo%cvBb(iCv,3) !srv 12.04.10 11.01.13
           delti=sign(min(abs(delti),dtiapr),delti)
         endif
c        if(cthi(iCv,isb).ne.2.65e0_R8)
c     1    write(*,*) iCv,isb,cthi(iCv,isb)
         smbtf(iCv,0) = smbtf(iCv,0)+t0*
     &    (switch%b2sifr_phm1*ehxp(iCv)-
     &    (switch%b2sifr_phm2*cthe(iCv,isb)*delte*ne(iCv)
     &    +switch%b2sifr_phm3*cthi(iCv,isb)*delti*
     &    (ne(iCv)-st_ext%ne(iCv)))*
     &    geo%cvVol(iCv)*geo%cvBb(iCv,0)/geo%cvBb(iCv,3))*geo%cvHz(iCv)  !srv 01.07.05
         smprb(iCv) = smprb(iCv) + 
     &        t0*switch%b2sifr_phm1*ehxp(iCv)*geo%cvHz(iCv) 
         smptb(iCv) = smptb(iCv) - t0*geo%cvHz(iCv)*
     &        (switch%b2sifr_phm2*cthe(iCv,isb)*delte*ne(iCv)
     &        +switch%b2sifr_phm3*cthi(iCv,isb)*delti*
     &        (ne(iCv)-st_ext%ne(iCv)))*
     1        geo%cvVol(iCv)*geo%cvBb(iCv,0)/geo%cvBb(iCv,3)
       enddo
      endif
      if (iout.ne.0) then                                           !srv 17.06.02 {
         write (chns,'(i3.3)') isb    
         call my_out_us(70,nCv,0,ehxp,'b2sifr_ehxp'//chns)
         call my_out_us(70,nCv,0,cthe(1,isb),'b2sifr_cthe'//chns)
         call my_out_us(70,nCv,0,cthi(1,isb),'b2sifr_cthi'//chns)
         do k=0,3
          write (chk,'(i1)') k  
          call my_out_us(70,nCv,0,smbtf(1,k),
     *                'b2sifr_smotf'//chk//chns)                    !srv 07.02.12
          call my_out_us(70,nCv,0,smbch(1,k),
     *                'b2sifr_smoch'//chk//chns)                       !srv 07.02.12
         enddo
      endif                                                            !srv 17.06.02 }

cdjm Jan2017 Store the linearised sources for balance
!WG_TODO      if (balance_netcdf.ne.0) then
!WG_TODO        b2sifr_smoch0to3(1:Cv,0:3,isb)=smbch
!WG_TODO        b2sifr_smotf0to3(1:nCv,0:3,isb)=smbtf
!WG_TODO      endif

* ..return
      ncall_b2sifr = ncall_b2sifr+1
      call subend ()
      return

*-----------------------------------------------------------------------
*.end b2sifr_

      end                                                             !srv 01.07.05 }
