*-----------------------------------------------------------------------
*.specification

      subroutine b2stel (nCv, nFc, ns, ismain,
     &  switch, geo, mpg, pl, dv, rt, rtw, srw)
      use b2mod_types
      use b2mod_diag
      use b2mod_tallies
      use b2mod_constants
      use b2mod_eirene_globals
      use b2mod_b2cmpa
      use b2mod_b2cmrc
#ifndef FIXED_POINT
      use b2mod_balance !djm Jan2017
     & , only : b2stel_sna_ion0to1, b2stel_sna_rec0to1,
     &          b2stel_smq_ion0to3, b2stel_smq_rec0to3,
     &          b2stel_she0to3,
     &          b2stel_she_ion0to3, b2stel_she_rec0to3,
     &          b2stel_shi_ion0to3, b2stel_shi_rec0to3,
     &          balance_netcdf
#endif
#ifdef ASTRA
      use b2mod_b2_to_astra
     , , only : svol, se, seh, sih, b2_to_astra_src  !srv 28.07.11 11.06.13
#endif
      use b2mod_switches
      use b2us_geo
      use b2us_map
      use b2us_plasma
! csc The following are not necessary for computation but are needed
!     for adjoint AD to avoid side-effect variables
      use b2mod_ad
     & , only : ncall_b2stel
      use b2mod_subsys
      use b2mod_dimensions
      implicit none
*   ..input arguments (unchanged on exit)
      integer nCv, nFc, ns, ismain
      type (switches), intent (in) :: switch
      type (geometry), intent (in) :: geo
      type (mapping), intent (in) :: mpg
      type (B2Plasma), intent (in) :: pl
      type (B2Derivatives), intent (inout) :: dv
      type (B2Rates), intent (in) :: rt
      type (B2RatesWork), intent (in) :: rtw
*   ..output arguments (unspecified on entry)
      type (B2SourceWork), intent (inout) :: srw
*-----------------------------------------------------------------------
*.documentation
*
*  1. purpose
*
*     B2STEL provides linearised expressions for the sources of
*     density, parallel momentum, electron heat and atom heat, due to
*     the electron-atom collision processes of ionisation,
*     recombination, and radiation.
*
*-----------------------------------------------------------------------
*.declarations

*   ..local variables
      integer iCv, is!, k
      integer, save :: nimp_tmp, jimp_tmp(DEF_NATM)
      character chns*3!, chk*1                                          !srv 30.06.08
      real (kind=R8) ::
     *  rf0, t0, t1, t2, tkin, t1i, t1n
      real (kind=R8) :: wrk0(nCv)
      real (kind=R8) ::
     *  sna0_ion(nCv,0:1,0:ns-1), sna0_rec(nCv,0:1,0:ns-1),
     *  smq0_ion(nCv,0:3,0:ns-1), smq0_rec(nCv,0:3,0:ns-1),
     &  she0_ion(nCv,0:3),        she0_rec(nCv,0:3),
     &  shi0_ion(nCv,0:3),        shi0_rec(nCv,0:3),
     &  shn0_ion(nCv,0:3),        shn0_rec(nCv,0:3)
*   ..procedures
      real (kind=R8) :: damax
      external xertst, ipgeti, ipgetr, damax, sfill
      external b2xvsg
*   ..initialisation

*-----------------------------------------------------------------------
*.computation

* ..preliminaries
*   ..subprogram start-up calls
      call subini ('b2stel')
*   ..set internal parameters on first call
*   ..identify impurity isonuclear sequences
      nimp_tmp=0
      jimp_tmp=0
      do is = 1, nspecies !{
        if(nfluids(is).gt.1) then !{
          nimp_tmp=nimp_tmp+1
          jimp_tmp(is)=nimp_tmp
        end if !}
      end do !}
*   ..test nCv, ns
      call xertst (0.lt.nCv, 'faulty argument nCv')
      call xertst (1.le.ns, 'faulty argument ns')
      call xertst (0.le.ismain.and.ismain.lt.ns,
     * 'faulty argument ismain')
*   ..extensive tests on first few calls
      if (ncall_b2stel.lt.3) then
*    ..test sign of vol
       call b2xvsg (nCv, geo%cvBb(1,3), 1, 'vol', '.gt.')
       call b2xvsg (nCv, geo%cvVol, 1, 'vol', '.gt.')
       call b2xvsg (nFc*2, geo%fcVol, 1, 'vol', '.gt.')
*    ..test sign of na, ni, ne, te, ti
       call b2xvsg (nCv*ns, pl%na, 1, 'na', '.gt.')
       call b2xvsg (nCv*2, dv%ni, 1, 'ni', '.gt.')
       call b2xvsg (nCv, dv%ne, 1, 'ne', '.gt.')
       call b2xvsg (nCv, pl%te, 1, 'te', '.gt.')
       call b2xvsg (nCv, pl%ti, 1, 'ti', '.gt.')
       call b2xvsg (nCv, pl%tn, 1, 'tn', '.gt.')
*   ..test facdrift, fac_ExB
      call xertst (0.0_R8.le.minval(dv%facdrift) .and.
     &     1.0_R8.ge.maxval(dv%facdrift), 'faulty argument facdrift')
      call xertst (0.0_R8.le.minval(dv%fac_ExB) .and.
     &     1.0_R8.ge.maxval(dv%fac_ExB), 'faulty argument fac_ExB')
*    ..test rate coefficients
       do is = 0, ns-2
        if (.not.lnext(is,is+1)) then
         call xertst (damax(nCv,rtw%rsa(1,is),1).eq.0.0_R8,
     &     'faulty argument rsa: nonzero for unrelated species')
         call xertst (damax(nCv,rtw%rra(1,is+1),1).eq.0.0_R8,
     &     'faulty argument rra: nonzero for unrelated species')
        endif
       enddo
      endif

* ..compute source terms
*   ..initialise sources to 0
      call sfill (nCv*2*ns, 0.0_R8, sna0_ion, 1)                        !srv 11.09.09
      call sfill (nCv*2*ns, 0.0_R8, sna0_rec, 1)
      call sfill (nCv*4*ns, 0.0_R8, smq0_ion, 1)
      call sfill (nCv*4*ns, 0.0_R8, smq0_rec, 1)
      call sfill (nCv*4, 0.0_R8, srw%she0, 1)
      call sfill (nCv*4, 0.0_R8, she0_ion, 1)
      call sfill (nCv*4, 0.0_R8, she0_rec, 1)
      call sfill (nCv*4, 0.0_R8, shi0_ion, 1)
      call sfill (nCv*4, 0.0_R8, shi0_rec, 1)
      call sfill (nCv*4, 0.0_R8, shn0_ion, 1)
      call sfill (nCv*4, 0.0_R8, shn0_rec, 1)
*   ..compute sources due to ionisation
      rsanareg=0.0_R8
      rsahireg=0.0_R8
      rsamoreg=0.0_R8
      srw%rsana=0.0_R8
      srw%rsahi=0.0_R8
      srw%rsamo=0.0_R8
      do is = 0, ns-2
       if (lnext(is,is+1)) then
*     ..particle source and heat source
        do iCv = 1, mpg%nCi ! internal cells only
          if(mpg%cvReg(iCv).ne.0) then
*       ..compute damping coefficient
           rf0 = switch%b2stel_rxm0*max(0.0_R8,rt%rlsa(iCv,1,is))
*       ..compute rate for is->is+1
           t0 = switch%b2stel_phm0*
     &                rtw%rsa(iCv,is)*geo%cvVol(iCv)*dv%ne(iCv)         !srv 05.06.18
*       ..compute particle source
           sna0_ion(iCv,0,is+1) = sna0_ion(iCv,0,is+1)+                 !srv 11.09.09 {
     &      (1.0_R8+rf0)*t0*pl%na(iCv,is)
           sna0_ion(iCv,1,is+1) = sna0_ion(iCv,1,is+1)-
     &      rf0*t0*pl%na(iCv,is)/pl%na(iCv,is+1)
           sna0_ion(iCv,0,is) = sna0_ion(iCv,0,is)+
     &      rf0*t0*pl%na(iCv,is)
           sna0_ion(iCv,1,is) = sna0_ion(iCv,1,is)-
     &      (1.0_R8+rf0)*t0  !srv 11.09.09 }
           srw%rsana(iCv,is) = t0*pl%na(iCv,is)
           rsanareg(mpg%cvReg(iCv),is)=rsanareg(mpg%cvReg(iCv),is)+
     1      t0*pl%na(iCv,is)
*       ..compute electron energy sink
            t2 = t0*rt%rpi(iCv,is)*ev
            she0_ion(iCv,0) = she0_ion(iCv,0)+(1.0_R8+rf0)*t2
            she0_ion(iCv,3) = she0_ion(iCv,3)-
     &                          rf0*t2/(dv%ne(iCv)*pl%te(iCv))
*       ..compute atom heat source
           t1 = t0*pl%na(iCv,is)*(1.0_R8-switch%BoRiS)*
     &        ((am(is)+am(is+1))/2.0_R8)*(mp/2.0_R8)*
     &        ((pl%ua(iCv,is)-pl%ua(iCv,is+1))**2)
           if (switch%tn_style.eq.0.or.nint(zn(is)).ne.1
     &                             .or..not.is_neutral(is)) then        !! combined ion-neutral, and/or non-hydrogenic (default)
             shi0_ion(iCv,0) = shi0_ion(iCv,0)+(1.0_R8+rf0)*t1
             shi0_ion(iCv,3) = shi0_ion(iCv,3)-
     &                         rf0*t1/(dv%ni(iCv,0)*pl%ti(iCv))
             rsahireg(mpg%cvReg(iCv),is)=rsahireg(mpg%cvReg(iCv),is)+t1
           elseif (switch%tn_style.eq.1) then !! pure ion energy eq.
             !! case of hydrogenic atoms
             t1i = t1 + t0*pl%na(iCv,is)*1.5_R8*pl%tn(iCv)
             shi0_ion(iCv,0) = shi0_ion(iCv,0)+(1.0_R8+rf0)*t1i
             shi0_ion(iCv,3) = shi0_ion(iCv,3)-
     &                         rf0*t1i/(dv%ni(iCv,0)*pl%ti(iCv))
             rsahireg(mpg%cvReg(iCv),is)=rsahireg(mpg%cvReg(iCv),is)+t1i
           else !! separate energy eq. for hydrogenic atoms
             !! case of hydrogenic atoms
             t1i = t1 + t0*pl%na(iCv,is)*1.5_R8*pl%tn(iCv)
             t1n = t0*pl%na(iCv,is)*(-1.5_R8*pl%tn(iCv))
             shi0_ion(iCv,0) = shi0_ion(iCv,0)+(1.0_R8+rf0)*t1i
             shi0_ion(iCv,3) = shi0_ion(iCv,3)-
     &                         rf0*t1i/(dv%ni(iCv,0)*pl%ti(iCv))
             shn0_ion(iCv,0) = shn0_ion(iCv,0)+(1.0_R8+rf0)*t1n
             shn0_ion(iCv,3) = shn0_ion(iCv,3)-
     &                         rf0*t1n/(dv%nn(iCv)*pl%tn(iCv))
             rsahireg(mpg%cvReg(iCv),is)=rsahireg(mpg%cvReg(iCv),is)+t1i
           endif
           srw%rsahi(iCv,is) = t1
          endif
        enddo
*     ..parallel momentum source
        do iCv = 1, mpg%nCi ! internal cells only
          if(mpg%cvReg(iCv).ne.0) then
*       ..compute damping coefficient
             rf0 = switch%b2stel_rg0 +
     &             switch%b2stel_rxm1*max(0.0_R8,rt%rlsa(iCv,1,is))
*       ..compute rate for is->is+1
             t0 = switch%b2stel_phm0*
     &                 rtw%rsa(iCv,is)*geo%cvVol(iCv)*dv%ne(iCv)*       !srv 05.06.18
     &           pl%na(iCv,is)*geo%cvHz(iCv)                            !srv 09.01.01
*       ..compute parallel momentum source
             t1 = ((am(is)+am(is+1))/2.0_R8)*mp*t0
             if (switch%b2stel_styl0.eq.0) then
                smq0_ion(iCv,0,is+1) = smq0_ion(iCv,0,is+1)+
     &           t1*(pl%ua(iCv,is)+rf0*pl%ua(iCv,is+1))
                smq0_ion(iCv,1,is+1) = smq0_ion(iCv,1,is+1)-
     &           rf0*t1
                smq0_ion(iCv,0,is) = smq0_ion(iCv,0,is)+
     &           rf0*t1*pl%ua(iCv,is)
                smq0_ion(iCv,1,is) = smq0_ion(iCv,1,is)-
     &           (1.0_R8+rf0)*t1
             else if (switch%b2stel_styl0.eq.1) then
*       ..compute parallel momentum source
               smq0_ion(iCv,2,is+1) = smq0_ion(iCv,2,is+1)+
     &          t1*(pl%ua(iCv,is)+rf0*pl%ua(iCv,is+1))/roxa(iCv,is+1)
               smq0_ion(iCv,3,is+1) = smq0_ion(iCv,3,is+1)-
     &          rf0*t1/roxa(iCv,is+1)
               smq0_ion(iCv,2,is) = smq0_ion(iCv,2,is)+
     &          rf0*t1*pl%ua(iCv,is)/roxa(iCv,is)
               smq0_ion(iCv,3,is) = smq0_ion(iCv,3,is)-
     &          (1.0_R8+rf0)*t1/roxa(iCv,is)
             else if (switch%b2stel_styl0.eq.2) then
               smq0_ion(iCv,0,is+1) = smq0_ion(iCv,0,is+1)+
     &          t1*(pl%ua(iCv,is)+rf0*pl%ua(iCv,is+1))
               smq0_ion(iCv,3,is+1) = smq0_ion(iCv,3,is+1)-
     &          rf0*t1/roxa(iCv,is+1)
               smq0_ion(iCv,0,is) = smq0_ion(iCv,0,is)+
     &          rf0*t1*pl%ua(iCv,is)
               smq0_ion(iCv,3,is) = smq0_ion(iCv,3,is)-
     &          (1.0_R8+rf0)*t1/roxa(iCv,is)
             else if (switch%b2stel_styl0.eq.3) then
               smq0_ion(iCv,2,is+1) = smq0_ion(iCv,2,is+1)+
     &          t1*(pl%ua(iCv,is)+rf0*pl%ua(iCv,is+1))/roxa(iCv,is+1)
               smq0_ion(iCv,1,is+1) = smq0_ion(iCv,1,is+1)-
     &          rf0*t1
               smq0_ion(iCv,2,is) = smq0_ion(iCv,2,is)+
     &          rf0*t1*pl%ua(iCv,is)/roxa(iCv,is)
               smq0_ion(iCv,1,is) = smq0_ion(iCv,1,is)-
     &         (1.0_R8+rf0)*t1
             endif
             srw%rsamo(iCv,is) = srw%rsamo(iCv,is)+t1*pl%ua(iCv,is)
             rsamoreg(mpg%cvReg(iCv),is) =
     &        rsamoreg(mpg%cvReg(iCv),is)+t1*pl%ua(iCv,is)
             srw%rsamo(iCv,is+1) = srw%rsamo(iCv,is+1)-t1*pl%ua(iCv,is)
             rsamoreg(mpg%cvReg(iCv),is+1) =
     &        rsamoreg(mpg%cvReg(iCv),is+1)-t1*pl%ua(iCv,is)
          endif
        enddo
       endif
      enddo

*   ..compute sources due to recombination
      rranareg=0.0_R8
      rrahireg=0.0_R8
      rramoreg=0.0_R8
      rqahereg=0.0_R8
      srw%rrana=0.0_R8
      srw%rrahi=0.0_R8
      srw%rramo=0.0_R8
      srw%rqahe=0.0_R8

      do is = 1, ns-1
       if (lnext(is-1,is)) then
*     ..particle source and heat source
        do iCv = 1, mpg%nCi ! internal cells only
         if(mpg%cvReg(iCv).ne.0) then
*       ..compute damping coefficient
           rf0 = switch%b2stel_rxm2*max(0.0_R8,-rt%rlra(iCv,1,is))
*       ..compute rate for is->is-1
           t0 = switch%b2stel_phm0*
     &                    rtw%rra(iCv,is)*geo%cvVol(iCv)*dv%ne(iCv)     !srv 05.06.18
*       ..compute particle source
           sna0_rec(iCv,0,is) = sna0_rec(iCv,0,is)+                     !srv 11.09.09 {
     &       rf0*t0*pl%na(iCv,is)
           sna0_rec(iCv,1,is) = sna0_rec(iCv,1,is)-(1.0_R8+rf0)*t0
           sna0_rec(iCv,0,is-1) = sna0_rec(iCv,0,is-1)+
     &      (1.0_R8+rf0)*t0*pl%na(iCv,is)
           sna0_rec(iCv,1,is-1) = sna0_rec(iCv,1,is-1)-                 !srv 11.09.09 {
     &      rf0*t0*pl%na(iCv,is)/pl%na(iCv,is-1)
           srw%rrana(iCv,is) = t0*pl%na(iCv,is)
           rranareg(mpg%cvReg(iCv),is)=rranareg(mpg%cvReg(iCv),is)+
     1      t0*pl%na(iCv,is)
*       ..compute electron heat source
           t2 = t0*pl%na(iCv,is)*rt%rpi(iCv,is)*ev
           she0_rec(iCv,0) = she0_rec(iCv,0)+(1.0_R8+rf0)*t2
           she0_rec(iCv,0) = she0_rec(iCv,3)-
     &                         rf0*t2/(dv%ne(iCv)*pl%te(iCv))
*       ..compute atom heat source
           if (switch%tn_style.eq.0.or.nint(zn(is)).ne.1
     &                             .or..not.is_neutral(is-1)) then      !! combined ion-neutral, and/or non-hydrogenic (default)
             t1 = t0*pl%na(iCv,is)*(1.0_R8-switch%BoRiS)*
     &         ((am(is)+am(is-1))/2.0_R8)*(mp/2.0_R8)*
     &         ((pl%ua(iCv,is)-pl%ua(iCv,is-1))**2)
             shi0_rec(iCv,0) = shi0_rec(iCv,0)+(1.0_R8+rf0)*t1
             shi0_rec(iCv,3) = shi0_rec(iCv,3)-
     *                         rf0*t1/(dv%ni(iCv,0)*pl%ti(iCv))
           elseif (switch%tn_style.eq.1) then !! pure ion
             !! case of hydrogenic atoms
             t1i = t0*pl%na(iCv,is)*(-1.5_R8*pl%ti(iCv))*
     *         (1.0_R8-switch%BoRiS)
             shi0_rec(iCv,0) = shi0_rec(iCv,0)+(1.0_R8+rf0)*t1i
             shi0_rec(iCv,3) = shi0_rec(iCv,3)-
     *                         rf0*t1i/(dv%ni(iCv,0)*pl%ti(iCv))
           else !! separate energy eq. for hydrogenic atoms
             !! case of hydrogenic atoms
             tkin = (am(is)+am(is-1))/2.0_R8*(mp/2.0_R8)*
     &              ((pl%ua(iCv,is)-pl%ua(iCv,is-1))**2)
             t1i = t0*pl%na(iCv,is)*(-1.5_R8*pl%ti(iCv))
             t1n = t0*pl%na(iCv,is)*(+1.5_R8*pl%ti(iCv)+tkin)
             shi0_rec(iCv,0) = shi0_rec(iCv,0)+(1.0_R8+rf0)*t1i
             shi0_rec(iCv,3) = shi0_rec(iCv,3)-
     *                         rf0*t1i/(dv%ni(iCv,0)*pl%ti(iCv))
             shn0_rec(iCv,0) = shn0_rec(iCv,0)+(1.0_R8+rf0)*t1n
             shn0_rec(icv,3) = shn0_rec(iCv,3)-
     *                         rf0*t1n/((dv%ni(iCv,0)-dv%ni(iCv,1))*
     *                         pl%tn(iCv))
           endif
           srw%rrahi(iCv,is) = t1
           rrahireg(mpg%cvReg(iCv),is) = rrahireg(mpg%cvReg(iCv),is)+t1
*       ..compute electron heat source for is->is-1, if not already included from ADPAK
           if(fix_recomb.ne.0.and.nint(zn(is)).le.zmax_recomb) then
c the default treatment is to have the entire ionisation potential
c difference lost, and the only additional cooling rate correction
c to be the Bremsstrahlung
c this treatment changes this --- the ionisation potential difference
c is given back as energy to the electrons and additional radiative
c losses should be included (in ADAS terms, the PRB file should be
c accounted for)
             srw%she0(iCv,0) = srw%she0(iCv,0)+(1.0_R8+rf0)*
     &        t0*rt%rpi(iCv,is)*ev*pl%na(iCv,is)
             srw%she0(iCv,3) = srw%she0(iCv,3)-rf0*
     &        t0*rt%rpi(iCv,is)*ev*pl%na(iCv,is)/(dv%ne(iCv)*pl%te(iCv))
             srw%rqahe(iCv,is) = srw%rqahe(iCv,is) -
     &        t0*rt%rpi(iCv,is)*ev*pl%na(iCv,is)
             rqahereg(mpg%cvReg(iCv),is) =
     &        rqahereg(mpg%cvReg(iCv),is) -
     &        t0*rt%rpi(iCv,is)*ev*pl%na(iCv,is)
           endif
         endif
        enddo
*     ..parallel momentum source
        do iCv = 1, mpg%nCi ! internal cells only
         if(mpg%cvReg(iCv).ne.0) then
*       ..compute damping coefficient
           rf0 = switch%b2stel_rg1 +
     &           switch%b2stel_rxm3*max(0.0_R8,-rt%rlra(iCv,1,is))
*       ..compute rate for is->is-1
           t0 = switch%b2stel_phm0*
     &              rtw%rra(iCv,is)*geo%cvVol(iCv)*dv%ne(iCv)*      !srv 05.06.18
     &         pl%na(iCv,is)*geo%cvHz(iCv)                          !srv 09.01.01
*       ..compute parallel momentum source
           t1 = ((am(is)+am(is-1))/2.0_R8)*mp*t0
           smq0_rec(iCv,0,is) = smq0_rec(iCv,0,is)+
     &      rf0*t1*pl%ua(iCv,is)
           smq0_rec(iCv,1,is) = smq0_rec(iCv,1,is)-(1.0_R8+rf0)*t1
           smq0_rec(iCv,0,is-1) = smq0_rec(iCv,0,is-1)+
     &      t1*(pl%ua(iCv,is)+rf0*pl%ua(iCv,is-1))
           smq0_rec(iCv,1,is-1) = smq0_rec(iCv,1,is-1)-rf0*t1
           srw%rramo(iCv,is) = srw%rramo(iCv,is)+t1*pl%ua(iCv,is)
           rramoreg(mpg%cvReg(iCv),is) =
     &      rramoreg(mpg%cvReg(iCv),is)+t1*pl%ua(iCv,is)
           srw%rramo(iCv,is-1) = srw%rramo(iCv,is-1)+t1*pl%ua(iCv,is)
           rramoreg(mpg%cvReg(iCv),is-1) =
     &      rramoreg(mpg%cvReg(iCv),is-1)+t1*pl%ua(iCv,is)
         endif
        enddo
       endif
      enddo

*   ..compute source terms for electron heat loss and radiation
      b2rad=0.0_R8
      b2radreg=0.0_R8
      b2bremreg=0.0_R8
      rqradreg=0.0_R8
      rqbrmreg=0.0_R8
      srw%rqrad=0.0_R8
      b2brem=0.0_R8
      rad_imp=0.0_R8
      do is = 0, ns-1
        do iCv = 1, mpg%nCi ! internal cells only
         if(mpg%cvReg(iCv).ne.0) then
*      ..compute damping coefficient
          rf0 = switch%b2stel_rxm4*max(0.0_R8,rt%rlqa(iCv,1,is))
*      ..compute electron cooling rate for is->any
*      ..includes line radiation, bremsstrahlung and ionisation contributions
*      ..for adpak rates, also contains recombination and CX terms
          t0 = switch%b2stel_phm0*(rtw%rqa(iCv,is)+rtw%rqr(iCv,is))
     &         *geo%cvVol(iCv)*dv%ne(iCv)*pl%na(iCv,is)    !srv 05.06.18
*      ..compute electron heat sink
          srw%she0(iCv,0) = srw%she0(iCv,0)+rf0*t0
          srw%she0(iCv,3) = srw%she0(iCv,3)-
     &     (1.0_R8+rf0)*t0/(dv%ne(iCv)*pl%te(iCv))
#ifndef FIXED_POINT
          if (balance_netcdf.ne.0) then
            b2stel_she0to3(iCv,0,is) = rf0*t0
            b2stel_she0to3(iCv,3,is) =
     &           -(1.0_R8+rf0)*t0/(dv%ne(iCv)*pl%te(iCv))
          endif
#endif
          srw%rqahe(iCv,is) = srw%rqahe(iCv,is) + t0
          rqahereg(mpg%cvReg(iCv),is)=rqahereg(mpg%cvReg(iCv),is)+t0
*      ..compute line radiation rate for is->any
          t0 = switch%b2stel_phm0*
     &         rtw%rrd(iCv,is)*geo%cvVol(iCv)*dv%ne(iCv)*
     &                                             pl%na(iCv,is)    !srv 05.06.18
*      ..compute line radiation
          srw%rqrad(iCv,is) = t0
          rqradreg(mpg%cvReg(iCv),is)=rqradreg(mpg%cvReg(iCv),is)+t0
          b2radreg(mpg%cvReg(iCv))=b2radreg(mpg%cvReg(iCv))-t0
          b2rad=b2rad-t0
          if (.not.is_neutral(is) .and. jimp_tmp(b2espcr(is)).ne.0) then
            rad_imp(jimp_tmp(b2espcr(is)))=
     .       rad_imp(jimp_tmp(b2espcr(is)))-t0
          endif
*      ..compute bremsstrahlung radiation rate for is->any
          t0 = switch%b2stel_phm0*
     &         rtw%rbr(iCv,is)*geo%cvVol(iCv)*dv%ne(iCv)*
     &                                             pl%na(iCv,is)    !srv 05.06.18
*      ..compute bremsstrahlung radiation
          srw%rqbrm(iCv,is) = t0
          rqbrmreg(mpg%cvReg(iCv),is)=rqbrmreg(mpg%cvReg(iCv),is)+t0
          b2bremreg(mpg%cvReg(iCv))=b2bremreg(mpg%cvReg(iCv))-t0
          b2brem=b2brem-t0
         endif
        enddo
      enddo

*  ..copy to srw
      srw%sna0 = sna0_ion + sna0_rec
      srw%smq0 = smq0_ion + smq0_rec
      srw%shi0 = shi0_ion + shi0_rec
      srw%shn0 = shn0_ion + shn0_rec

#ifdef ASTRA
      do is = 0, ns-1                                                 !srv 28.07.11 {
        if(is.eq.ismain) then
          wrk0=sna0_ion(:,0,is)+sna0_ion(:,1,is)*pl%na(:,is)
          call b2_to_astra_src(mpg, nx, ny, geo%cvVol, wrk0, svol, se)
          wrk0=sna0_rec(:,:,0,is)+sna0_rec(:,1,is)*pl%na(:,is)
          call b2_to_astra_src(mpg, nx, ny, geo%cvVol, wrk0, svol, se)
        endif
      enddo                                                           !srv 28.07.11 }

      wrk0=srw%she0(:,0)+
     &     srw%she0(:,1)*pl%te(:)+srw%she0(:,2)*dv%ne+
     &     srw%she0(:,3)*dv%ne*pl%te
      call b2_to_astra_src(mpg, nx, ny, geo%cvVol, wrk0, svol, seh)              !srv 11.06.13

      wrk0=shi0_ion(:,0)+
     &     shi0_ion(:,1)*pl%ti+shi0_ion(:,2)*dv%ni(:,0)+
     &     shi0_ion(:,3)*dv%ni(:,0)*pl%ti
      call b2_to_astra_src(mpg, nx, ny, geo%cvVol, wrk0, svol, sih)              !srv 11.06.13

      wrk0=shi0_rec(:,0)+
     &     shi0_rec(:,1)*pl%ti+shi0_rec(:,2)*dv%ni(:,0)+
     &     shi0_rec(:,3)*dv%ni(:,0)*pl%ti
      call b2_to_astra_src(mpg, nx, ny, geo%cvVol, wrk0, svol, sih)              !srv 11.06.13
#endif

      if (switch%b2stel_iout.eq.1 .or. switch%b2npmo_iout.eq.1 .or.
     &    switch%iout_b2wdat .eq. 4) then                              !srv 30.06.08 11.09.09 {
        do is = 0, ns-1
          write (chns,'(i3.3)') is
!         do k=0,3
!          write (chk,'(i1)') k
!          call my_out(70,nx,ny,smq0_ion(-1,-1,k,is),
!     &                       'b2stel_smq_ion'//chk//chns)
!          call my_out(70,nx,ny,smq0_rec(-1,-1,k,is),
!     &                       'b2stel_smq_rec'//chk//chns)
!         enddo
          wrk0=smq0_ion(:,0,is)+
     &         smq0_ion(:,1,is)*pl%ua(:,is)+
     &         smq0_ion(:,2,is)*pl%na(:,is)*mp*am(is)+
     &         smq0_ion(:,3,is)*pl%ua(:,is)*pl%na(:,is)*mp*am(is)
          call my_out_us(70,nCv,0,wrk0,'b2stel_smq_ion'//chns)
          wrk0=smq0_rec(:,0,is)+
     &         smq0_rec(:,1,is)*pl%ua(:,is)+
     &         smq0_rec(:,2,is)*pl%na(:,is)*mp*am(is)+
     &         smq0_rec(:,3,is)*pl%ua(:,is)*pl%na(:,is)*mp*am(is)
          call my_out_us(70,nCv,0,wrk0,'b2stel_smq_rec'//chns)
        enddo
      endif
      if (switch%b2stel_iout.eq.1 .or. switch%b2npco_iout.eq.1
     &                           .or. switch%iout_b2wdat .eq. 4) then
        do is = 0, ns-1
          write (chns,'(i3.3)') is
!         do k=0,1
!          write (chk,'(i1)') k
!          call my_out(70,nx,ny,sna0_ion(-1,-1,k,is),
!     &                       'b2stel_sna_ion'//chk//chns)
!          call my_out(70,nx,ny,sna0_rec(-1,-1,k,is),
!     &                       'b2stel_sna_rec'//chk//chns)
!         enddo
          wrk0=sna0_ion(:,0,is)+sna0_ion(:,1,is)*pl%na(:,is)
          call my_out_us(70,nCv,0,wrk0,'b2stel_sna_ion'//chns)
          wrk0=sna0_rec(:,0,is)+sna0_rec(:,1,is)*pl%na(:,is)
          call my_out_us(70,nCv,0,wrk0,'b2stel_sna_rec'//chns)

        enddo
      endif
      wrk0=srw%she0(:,0)+srw%she0(:,1)*pl%te+
     &     srw%she0(:,2)*dv%ne + srw%she0(:,3)*pl%te*dv%ne
      if (switch%b2stel_iout.eq.1 .or. switch%b2npht_iout.eq.1
     &                            .or. switch%iout_b2wdat .eq. 4) then
!        do k=0,3
!         write (chk,'(i1)') k
!         call my_out(70,nx,ny,shi0_ion(-1,-1,k),'b2stel_shi_ion'//chk)
!         call my_out(70,nx,ny,shi0_rec(-1,-1,k),'b2stel_shi_rec'//chk)
!        enddo
!      call my_out(70,nx,ny,srw%she0(-1,-1,k),'b2stel_she')
        call my_out_us(70,nCv,0,wrk0,'b2stel_she_rad')
        wrk0=shi0_ion(:,0)+shi0_ion(:,1)*pl%ti+
     &       shi0_ion(:,2)*dv%ni(:,0)+
     &       shi0_ion(:,3)*pl%ti*dv%ni(:,0)
        call my_out_us(70,nCv,0,wrk0,'b2stel_shi_ion')                  !srv 01.11.10
        wrk0=shi0_rec(:,0)+shi0_rec(:,1)*pl%ti+
     &       shi0_rec(:,2)*dv%ni(:,0)+
     &       shi0_rec(:,3)*pl%ti*dv%ni(:,0)
        call my_out_us(70,nCv,0,wrk0,'b2stel_shi_rec')                  !srv 01.11.10

        do is=0,ns-1                                                 !srv 03.02.12 {
          write(chns,'(i3.3)') is
          call my_out_us(70,nCv,0,srw%rqrad(1,is),'b2stel_rqrad'//chns)
          call my_out_us(70,nCv,0,srw%rqbrm(1,is),'b2stel_rqbrm'//chns)
        enddo                                                        !srv 03.02.12 }
      endif                                                          !srv 30.06.08 11.09.09 }

#ifndef FIXED_POINT
cdjm Jan2017 Keep linearised sources for balance
      if (balance_netcdf.ne.0) then
        b2stel_sna_ion0to1 = sna0_ion
        b2stel_sna_rec0to1 = sna0_rec
        b2stel_smq_ion0to3 = smq0_ion
        b2stel_smq_rec0to3 = smq0_rec
        b2stel_she_ion0to3 = she0_ion
        b2stel_she_rec0to3 = she0_rec
        b2stel_shi_ion0to3 = shi0_ion
        b2stel_shi_rec0to3 = shi0_rec
      endif
#endif

* ..return
      ncall_b2stel = ncall_b2stel+1
      call subend ()
      return

*-----------------------------------------------------------------------
*.end b2stel

      contains

      function roxa(iCv,is)
      implicit none
      real (kind=R8) :: roxa
      integer, intent(in) :: iCv, is

      roxa = am(is)*mp*pl%na(iCv,is)

      return
      end function roxa

      end subroutine b2stel
