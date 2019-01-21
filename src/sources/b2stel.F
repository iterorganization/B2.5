*-----------------------------------------------------------------------
*.specification

      subroutine b2stel (nCv, nFc, nVx, ns, ismain, BoRiS,
     &  geo, mpg, pl, dv, rt, rtw, srw)
      use b2mod_types
!      use b2mod_diag
      use b2mod_tallies
      use b2mod_constants
!      use b2mod_sources
      use b2mod_eirene_globals
      use b2mod_b2cmpa
      use b2mod_b2cmrc
!WG_TODO      use b2mod_balance !djm Jan2017
!WG_TODO     & , only : b2stel_sna_ion0to1, b2stel_sna_rec0to1,
!WG_TODO     &          b2stel_smq_ion0to3, b2stel_smq_rec0to3,
!WG_TODO     &          b2stel_she0to3,
!WG_TODO     &          b2stel_shi_ion0to3, b2stel_shi_rec0to3,
!WG_TODO     &          balance_netcdf
#ifdef ASTRA
      use b2mod_b2_to_astra
     , , only : svol, se, seh, sih, b2_to_astra_src  !srv 28.07.11 11.06.13
#endif
      use b2us_geo
      use b2us_map
      use b2us_plasma
      implicit none
*   ..input arguments (unchanged on exit)
      integer nCv, nFc, nVx, ns, ismain
      real (kind=R8) ::
     *  BoRiS
      type (geometry), intent (in) :: geo
      type (mapping), intent (in) :: mpg
      type (B2Plasma), intent (in) :: pl
      type (B2Derivatives), intent (inout) :: dv
      type (B2Rates), intent (in) :: rt
      type (B2RatesWork), intent (in) :: rtw
*   ..output arguments (unspecified on entry)
      type (B2SourceWork), intent (inout) :: srw
*   ..common blocks
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

#     include <DIMENSIONS.F>
*   ..local variables
      integer ncall, styl0, iCv, is, js, fix_recomb_energy, k
!WG_TODO      integer, save :: nimp_tmp, jimp_tmp(DEF_NATM)
      integer, save :: iout = 0, iout_b2wdat = 0,
     &                 iout_b2npc9 = 0, iout_b2npco = 0,
     &                 iout_b2nph9 = 0, iout_b2npht = 0,
     &                 iout_b2npmo =0
      character chns*3, chk*1                                           !srv 30.06.08
      real (kind=R8) ::
     *  rg0, rg1, rxm0, rxm1, rxm2, rxm3, rxm4, phm0
      save ncall, styl0, rg0, rg1, rxm0, rxm1, rxm2, rxm3, rxm4, phm0,
     *  fix_recomb_energy
      real (kind=R8) ::
     *  rf0, t0, t1
      real (kind=R8) :: wrk0(nCv)
      real (kind=R8) ::
     *  sna0_ion(nCv,0:1,0:ns-1), sna0_rec(nCv,0:1,0:ns-1),
     *  smq0_ion(nCv,0:3,0:ns-1), smq0_rec(nCv,0:3,0:ns-1),
     &  shi0_ion(nCv,0:3),        shi0_rec(nCv,0:3)
*   ..procedures
      intrinsic abs
      logical lnext
      real (kind=R8) ::
     *   samax, roxa
      external subini, subend, xertst, ipgeti, ipgetr, samax, sfill
      external b2xvsg
      lnext(is,js) = zamax(is).lt.zamin(js).and.zn(is).eq.zn(js).and.
     &  abs(am(is)-am(js)).lt.1.0e-3_R8
      roxa(iCv,is) = am(is)*mp*pl%na(iCv,is)
*   ..initialisation
      data ncall/0/, styl0/0/, rg0/1.0_R8/,
     &  rg1/1.0_R8/, rxm0/0.0_R8/,
     &  rxm1/0.0_R8/, rxm2/0.0_R8/,
     &  rxm3/0.0_R8/, rxm4/0.0_R8/, phm0/1.0_R8/,           !srv 05.06.18
     &  fix_recomb_energy/0/

*-----------------------------------------------------------------------
*.computation

* ..preliminaries
*   ..subprogram start-up calls
      call subini ('b2stel')
*   ..set internal parameters on first call
      if (ncall.eq.0) then
       call ipgetr ('b2stel_rg0', rg0)
       call ipgetr ('b2stel_rg1', rg1)
       call xertst (0.0_R8.le.rg0.and.0.0_R8.le.rg1,
     &   'faulty internal parameter rg0/1')
       call ipgetr ('b2stel_rxm0', rxm0)
       call ipgetr ('b2stel_rxm1', rxm1)
       call ipgetr ('b2stel_rxm2', rxm2)
       call ipgetr ('b2stel_rxm3', rxm3)
       call ipgetr ('b2stel_rxm4', rxm4)
       call ipgetr ('b2stel_phm0', phm0)
       call xertst (0.0_R8.le.rxm0.and.
     &   0.0_R8.le.rxm1.and.0.0_R8.le.rxm2.and.
     &   0.0_R8.le.rxm3.and.0.0_R8.le.rxm4,
     &   'faulty internal parameter rxm0/1/2/3/4')
       call xertst (0.0_R8.le.phm0, 'faulty parameter phm0')
       call ipgeti ('b2stel_styl0', styl0)
       call xertst (0.le.styl0.and.styl0.lt.4,
     &   'faulty internal parameter styl0')
       call ipgeti ('b2stel_fix_recomb_energy', fix_recomb_energy)
       call xertst (.not.(adpak_used.and.fix_recomb_energy.ne.0),
     &   'fix_recomb_energy should not be used with ADPAK rates!')
       if(fix_recomb_energy.ne.0) then
         write(*,*)
     &    'WARNING: fix_recomb_energy <> 0 not fully tested yet!'
       endif
*   ..identify impurity homonuclear sequences
!WG_TODO       nimp_tmp=0
!WG_TODO       jimp_tmp=0
!WG_TODO       do is = 1, nspecies !{
!WG_TODO         if(nfluids(is).gt.1) then !{
!WG_TODO           nimp_tmp=nimp_tmp+1
!WG_TODO           jimp_tmp(is)=nimp_tmp
!WG_TODO         end if !}
!WG_TODO       end do !}
       call ipgeti ('b2stel_iout', iout)                         !srv 30.06.08
       call ipgeti ('b2npc9_iout', iout_b2npc9)                  !srv 11.09.09 {
       call ipgeti ('b2npco_iout', iout_b2npco)
       call ipgeti ('b2npmo_iout', iout_b2npmo)
       call ipgeti ('b2nph9_iout', iout_b2nph9)
       call ipgeti ('b2npht_iout', iout_b2npht)                  !srv 11.09.09 }
       call ipgeti ('b2wdat_iout', iout_b2wdat)                  !iys 17.04.17
      endif
*   ..test nCv, ns
      call xertst (0.le.nCv, 'faulty argument nx, ny')
      call xertst (1.le.ns, 'faulty argument ns')
*   ..extensive tests on first few calls
      if (ncall.lt.3) then
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
*   ..test facdrift, fac_ExB
      call xertst (0.0_R8.le.minval(dv%facdrift) .and.
     &     1.0_R8.ge.maxval(dv%facdrift), 'faulty argument facdrift')
      call xertst (0.0_R8.le.minval(dv%fac_ExB) .and.
     &     1.0_R8.ge.maxval(dv%fac_ExB), 'faulty argument fac_ExB')
*    ..test rate coefficients
       do is = 0, ns-2
        if (.not.lnext(is,is+1)) then
         call xertst (samax(nCv,rtw%rsa(1,is),1).eq.0.0_R8,
     &     'faulty argument rsa: nonzero for unrelated species')
         call xertst (samax(nCv,rtw%rra(1,is+1),1).eq.0.0_R8,
     &     'faulty argument rra: nonzero for unrelated species')
        endif
       enddo
      endif

* ..compute source terms
*   ..initialise sources to 0
      call sfill (nCv*2*ns, 0.0_R8, sna0_ion, 1)                      !srv 11.09.09
      call sfill (nCv*2*ns, 0.0_R8, sna0_rec, 1)
      call sfill (nCv*4*ns, 0.0_R8, smq0_ion, 1)
      call sfill (nCv*4*ns, 0.0_R8, smq0_rec, 1)
      call sfill (nCv*4, 0.0_R8, srw%she0, 1)
      call sfill (nCv*4, 0.0_R8, shi0_ion, 1)
      call sfill (nCv*4, 0.0_R8, shi0_rec, 1)
*   ..compute sources due to ionisation
      rsanareg=0.0_R8
      rsahireg=0.0_R8
      rsamoreg=0.0_R8
!WG_TODO      rsana=0.0_R8
!WG_TODO      rsahi=0.0_R8
!WG_TODO      rsamo=0.0_R8
      do is = 0, ns-2
       if (lnext(is,is+1)) then
*     ..particle source and heat source
        do iCv = 1, mpg%nCi ! internal cells only
          if(mpg%cvReg(iCv).ne.0) then
*       ..compute damping coefficient
           rf0 = rxm0*max(0.0_R8,rt%rlsa(iCv,1,is))
*       ..compute rate for is->is+1
           t0 = phm0*rtw%rsa(iCv,is)*geo%cvVol(iCv)*dv%ne(iCv)      !srv 05.06.18
*       ..compute particle source
           sna0_ion(iCv,0,is+1) = sna0_ion(iCv,0,is+1)+             !srv 11.09.09 {
     &      (1.0_R8+rf0)*t0*pl%na(iCv,is)
           sna0_ion(iCv,1,is+1) = sna0_ion(iCv,1,is+1)-
     &      rf0*t0*pl%na(iCv,is)/pl%na(iCv,is+1)
           sna0_ion(iCv,0,is) = sna0_ion(iCv,0,is)+
     &      rf0*t0*pl%na(iCv,is)
           sna0_ion(iCv,1,is) = sna0_ion(iCv,1,is)-
     &      (1.0_R8+rf0)*t0  !srv 11.09.09 }
!WG_TODO           rsana(iCv,is) = t0*pl%na(iCv,is)
           rsanareg(mpg%cvReg(iCv),is)=rsanareg(mpg%cvReg(iCv),is)+
     1      t0*pl%na(iCv,is)
*       ..compute atom heat source
           t1 = t0*pl%na(iCv,is)*(1.0_R8-BoRiS)*
     &       ((am(is)+am(is+1))/2.0_R8)*(mp/2.0_R8)*
     &       ((pl%ua(iCv,is)-pl%ua(iCv,is+1))**2)
           shi0_ion(iCv,0) = shi0_ion(iCv,0)+(1.0_R8+rf0)*t1
           shi0_ion(iCv,3) = shi0_ion(iCv,3)-
     &                         rf0*t1/(dv%ni(iCv,0)*pl%ti(iCv))
!WG_TODO           rsahi(iCv,is) = t1
           rsahireg(mpg%cvReg(iCv),is)=rsahireg(mpg%cvReg(iCv),is)+t1
          endif
        enddo
*     ..parallel momentum source
        do iCv = 1, mpg%nCi ! internal cells only
          if(mpg%cvReg(iCv).ne.0) then
*       ..compute damping coefficient
             rf0 = rg0+rxm1*max(0.0_R8,rt%rlsa(iCv,1,is))
*       ..compute rate for is->is+1
             t0 = phm0*rtw%rsa(iCv,is)*geo%cvVol(iCv)*dv%ne(iCv)*   !srv 05.06.18
     &           pl%na(iCv,is)*geo%cvHz(iCv)                        !srv 09.01.01
*       ..compute parallel momentum source
             t1 = ((am(is)+am(is+1))/2.0_R8)*mp*t0
             if (styl0.eq.0) then
                smq0_ion(iCv,0,is+1) = smq0_ion(iCv,0,is+1)+
     &           t1*(pl%ua(iCv,is)+rf0*pl%ua(iCv,is+1))
                smq0_ion(iCv,1,is+1) = smq0_ion(iCv,1,is+1)-
     &           rf0*t1
                smq0_ion(iCv,0,is) = smq0_ion(iCv,0,is)+
     &           rf0*t1*pl%ua(iCv,is)
                smq0_ion(iCv,1,is) = smq0_ion(iCv,1,is)-
     &           (1.0_R8+rf0)*t1
             else if (styl0.eq.1) then
*       ..compute parallel momentum source
               smq0_ion(iCv,2,is+1) = smq0_ion(iCv,2,is+1)+
     &          t1*(pl%ua(iCv,is)+rf0*pl%ua(iCv,is+1))/roxa(iCv,is+1)
               smq0_ion(iCv,3,is+1) = smq0_ion(iCv,3,is+1)-
     &          rf0*t1/roxa(iCv,is+1)
               smq0_ion(iCv,2,is) = smq0_ion(iCv,2,is)+
     &          rf0*t1*pl%ua(iCv,is)/roxa(iCv,is)
               smq0_ion(iCv,3,is) = smq0_ion(iCv,3,is)-
     &          (1.0_R8+rf0)*t1/roxa(iCv,is)
             else if (styl0.eq.2) then
               smq0_ion(iCv,0,is+1) = smq0_ion(iCv,0,is+1)+
     &          t1*(pl%ua(iCv,is)+rf0*pl%ua(iCv,is+1))
               smq0_ion(iCv,3,is+1) = smq0_ion(iCv,3,is+1)-
     &          rf0*t1/roxa(iCv,is+1)
               smq0_ion(iCv,0,is) = smq0_ion(iCv,0,is)+
     &          rf0*t1*pl%ua(iCv,is)
               smq0_ion(iCv,3,is) = smq0_ion(iCv,3,is)-
     &          (1.0_R8+rf0)*t1/roxa(iCv,is)
             else if (styl0.eq.3) then
               smq0_ion(iCv,2,is+1) = smq0_ion(iCv,2,is+1)+
     &          t1*(pl%ua(iCv,is)+rf0*pl%ua(iCv,is+1))/roxa(iCv,is+1)
               smq0_ion(iCv,1,is+1) = smq0_ion(iCv,1,is+1)-
     &          rf0*t1
               smq0_ion(iCv,2,is) = smq0_ion(iCv,2,is)+
     &          rf0*t1*pl%ua(iCv,is)/roxa(iCv,is)
               smq0_ion(iCv,1,is) = smq0_ion(iCv,1,is)-
     &         (1.0_R8+rf0)*t1
             endif
!WG_TODO             rsamo(iCv,is) = rsamo(iCv,is)+t1*pl%ua(iCv,is)
             rsamoreg(mpg%cvReg(iCv),is) =
     &        rsamoreg(mpg%cvReg(iCv),is)+t1*pl%ua(iCv,is)
!WG_TODO             rsamo(iCv,is+1) = rsamo(iCv,is+1)-t1*pl%ua(iCv,is)
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
!WG_TODO      rrana=0.0_R8
!WG_TODO      rrahi=0.0_R8
!WG_TODO      rramo=0.0_R8
!WG_TODO      rqahe=0.0_R8
      do is = 1, ns-1
       if (lnext(is-1,is)) then
*     ..particle source and heat source
        do iCv = 1, mpg%nCi ! internal cells only
         if(mpg%cvReg(iCv).ne.0) then
*       ..compute damping coefficient
           rf0 = rxm2*max(0.0_R8,-rt%rlra(iCv,1,is))
*       ..compute rate for is->is-1
           t0 = phm0*rtw%rra(iCv,is)*geo%cvVol(iCv)*dv%ne(iCv)    !srv 05.06.18
*       ..compute particle source
           sna0_rec(iCv,0,is) = sna0_rec(iCv,0,is)+               !srv 11.09.09 {
     &       rf0*t0*pl%na(iCv,is)
           sna0_rec(iCv,1,is) = sna0_rec(iCv,1,is)-(1.0_R8+rf0)*t0
           sna0_rec(iCv,0,is-1) = sna0_rec(iCv,0,is-1)+
     &      (1.0_R8+rf0)*t0*pl%na(iCv,is)
           sna0_rec(iCv,1,is-1) = sna0_rec(iCv,1,is-1)-           !srv 11.09.09 {
     &      rf0*t0*pl%na(iCv,is)/pl%na(iCv,is-1)
!WG_TODO           rrana(iCv,is) = t0*pl%na(iCv,is)
           rranareg(mpg%cvReg(iCv),is)=rranareg(mpg%cvReg(iCv),is)+
     1      t0*pl%na(iCv,is)
*       ..compute atom heat source
           t1 = t0*pl%na(iCv,is)*(1.0_R8-BoRiS)*
     &       ((am(is)+am(is-1))/2.0_R8)*(mp/2.0_R8)*
     &       ((pl%ua(iCv,is)-pl%ua(iCv,is-1))**2)
           shi0_rec(iCv,0) = shi0_rec(iCv,0)+(1.0_R8+rf0)*t1
           shi0_rec(iCv,3) = shi0_rec(iCv,3)-
     *                         rf0*t1/(dv%ni(iCv,0)*pl%ti(iCv))
!WG_TODO           rrahi(iCv,is) = t1
           rrahireg(mpg%cvReg(iCv),is)=rrahireg(mpg%cvReg(iCv),is)+t1
*       ..compute electron heat source for is->is-1, if not already included from ADPAK
           if(fix_recomb_energy.ne.0) then
c the default treatment is to have the entire ionization potential
c difference lost, and the only additional cooling rate correction
c to be the Bremsstrahlung
c this treatment changes this --- the ionization potential difference
c is given back as energy to the electrons and additional radiative
c losses should be included (in ADAS terms, the PRB file should be
c accounted for)
             srw%she0(iCv,0)=srw%she0(iCv,0)+(1.0_R8+rf0)*
     &        t0*rt%rpi(iCv,is)*ev*pl%na(iCv,is)
             srw%she0(iCv,3)=srw%she0(iCv,3)-rf0*
     &        t0*rt%rpi(iCv,is)*ev*pl%na(iCv,is)/(dv%ne(iCv)*pl%te(iCv))
!WG_TODO             rqahe(iCv,is) = rqahe(iCv,is) -
!WG_TODO     &        t0*rt%rpi(iCv,is)*ev*pl%na(iCv,is)
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
           rf0 = rg1+rxm3*max(0.0_R8,-rt%rlra(iCv,1,is))
*       ..compute rate for is->is-1
           t0 =phm0*rtw%rra(iCv,is)*geo%cvVol(iCv)*dv%ne(iCv)*      !srv 05.06.18
     &         pl%na(iCv,is)*geo%cvHz(iCv)                          !srv 09.01.01
*       ..compute parallel momentum source
           t1 = ((am(is)+am(is-1))/2.0_R8)*mp*t0
           smq0_rec(iCv,0,is) = smq0_rec(iCv,0,is)+
     &      rf0*t1*pl%ua(iCv,is)
           smq0_rec(iCv,1,is) = smq0_rec(iCv,1,is)-(1.0_R8+rf0)*t1
           smq0_rec(iCv,0,is-1) = smq0_rec(iCv,0,is-1)+
     &      t1*(pl%ua(iCv,is)+rf0*pl%ua(iCv,is-1))
           smq0_rec(iCv,1,is-1) = smq0_rec(iCv,1,is-1)-rf0*t1
!WG_TODO           rramo(iCv,is) = rramo(iCv,is)+t1*pl%ua(iCv,is)
           rramoreg(mpg%cvReg(iCv),is) =
     &      rramoreg(mpg%cvReg(iCv),is)+t1*pl%ua(iCv,is)
!WG_TODO           rramo(iCv,is-1) = rramo(iCv,is-1)+t1*pl%ua(iCv,is)
           rramoreg(mpg%cvReg(iCv),is-1) =
     &      rramoreg(mpg%cvReg(iCv),is-1)+t1*pl%ua(iCv,is)
         endif
        enddo
       endif
      enddo
*   ..compute source terms for electron heat loss and radiation
!WG_TODO      b2rad=0.0_R8
      b2radreg=0.0_R8
      b2bremreg=0.0_R8
      rqradreg=0.0_R8
      rqbrmreg=0.0_R8
!WG_TODO      rqrad=0.0_R8
!WG_TODO      b2brem=0.0_R8
!WG_TODO      rad_imp=0.0_R8
      do is = 0, ns-1
        do iCv = 1, mpg%nCi ! internal cells only
         if(mpg%cvReg(iCv).ne.0) then
*      ..compute damping coefficient
          rf0 = rxm4*max(0.0_R8,rt%rlqa(iCv,1,is))
*      ..compute electron cooling rate for is->any
*      ..includes line radiation, bremsstrahlung and ionisation contributions
*      ..if adpak_used, also contains recombination and CX terms
          t0 = phm0*rtw%rqa(iCv,is)*geo%cvVol(iCv)*dv%ne(iCv)*
     &                                             pl%na(iCv,is)    !srv 05.06.18
*      ..compute electron heat sink
          srw%she0(iCv,0) = srw%she0(iCv,0)+rf0*t0
          srw%she0(iCv,3) = srw%she0(iCv,3)-
     &     (1.0_R8+rf0)*t0/(dv%ne(iCv)*pl%te(iCv))
!WG_TODO          rqahe(iCv,is) = rqahe(iCv,is) + t0
          rqahereg(mpg%cvReg(iCv),is)=rqahereg(mpg%cvReg(iCv),is)+t0
*      ..compute line radiation rate for is->any
          t0 = phm0*rtw%rrd(iCv,is)*geo%cvVol(iCv)*dv%ne(iCv)*
     &                                             pl%na(iCv,is)    !srv 05.06.18
*      ..compute line radiation
!WG_TODO          rqrad(iCv,is) = t0
          rqradreg(mpg%cvReg(iCv),is)=rqradreg(mpg%cvReg(iCv),is)+t0
          b2radreg(mpg%cvReg(iCv))=b2radreg(mpg%cvReg(iCv))-t0
!WG_TODO          b2rad=b2rad-t0
!WG_TODO          if (.not.is_neutral(is) .and. jimp_tmp(b2espcr(is)).ne.0) then
!WG_TODO            rad_imp(jimp_tmp(b2espcr(is)))=
!WG_TODO     .       rad_imp(jimp_tmp(b2espcr(is)))-t0
!WG_TODO          endif
*      ..compute bremsstrahlung radiation rate for is->any
          t0 = phm0*rtw%rbr(iCv,is)*geo%cvVol(iCv)*dv%ne(iCv)*
     &                                             pl%na(iCv,is)    !srv 05.06.18
*      ..compute bremsstrahlung radiation
!WG_TODO          rqbrm(iCv,is) = t0
          rqbrmreg(mpg%cvReg(iCv),is)=rqbrmreg(mpg%cvReg(iCv),is)+t0
          b2bremreg(mpg%cvReg(iCv))=b2bremreg(mpg%cvReg(iCv))-t0
!WG_TODO          b2brem=b2brem-t0
         endif
        enddo
      enddo

*  ..copy to srw
      srw%sna0 = sna0_ion + sna0_rec
      srw%smq0 = smq0_ion + smq0_rec
      srw%shi0 = shi0_ion + shi0_rec

#ifdef ASTRA
      do is = 0, ns-1                                                 !srv 28.07.11 {
        if(is.eq.ismain) then
          wrk0=sna0_ion(:,:,0,is)+sna0_ion(:,:,1,is)*pl%na(:,:,is)
          call b2_to_astra_src(nx, ny, vol, wrk0, svol, se)
          wrk0=sna0_rec(:,:,0,is)+sna0_rec(:,:,1,is)*pl%na(:,:,is)
          call b2_to_astra_src(nx, ny, vol, wrk0, svol, se)
        endif
      enddo                                                           !srv 28.07.11 }

      wrk0=srw%she0(:,0)+
     &     srw%she0(:,1)*te(:)+srw%she0(:,2)*ne+
     &     srw%she0(:,3)*ne*te
      call b2_to_astra_src(nx, ny, vol, wrk0, svol, seh)              !srv 11.06.13

      wrk0=shi0_ion(:,0)+
     &     shi0_ion(:,1)*ti+shi0_ion(:,2)*ni(:,0)+
     &     shi0_ion(:,3)*ni(:,0)*ti
      call b2_to_astra_src(nx, ny, vol, wrk0, svol, sih)              !srv 11.06.13

      wrk0=shi0_rec(:,0)+
     &     shi0_rec(:,1)*ti+shi0_rec(:,2)*ni(:,0)+
     &     shi0_rec(:,3)*ni(:,0)*ti
      call b2_to_astra_src(nx, ny, vol, wrk0, svol, sih)              !srv 11.06.13
#endif

      if (iout.eq.1 .or. iout_b2npmo.eq.1 .or.
     &    iout_b2wdat .eq. 4) then                                    !srv 30.06.08 11.09.09 {
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
      if (iout.eq.1 .or. (iout_b2npc9.eq.1 .or. iout_b2npco.eq.1)
     &                 .or. iout_b2wdat .eq. 4) then
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
      if (iout.eq.1 .or. (iout_b2nph9.eq.1 .or. iout_b2npht.eq.1)
     &                 .or. iout_b2wdat .eq. 4) then
!        do k=0,3
!         write (chk,'(i1)') k
!         call my_out(70,nx,ny,shi0_ion(-1,-1,k),'b2stel_shi_ion'//chk)
!         call my_out(70,nx,ny,shi0_rec(-1,-1,k),'b2stel_shi_rec'//chk)
!        enddo
!      call my_out(70,nx,ny,srw%she0(-1,-1,k),'b2stel_she')
        wrk0=srw%she0(:,0)+srw%she0(:,1)*pl%te+
     &       srw%she0(:,2)*dv%ne + srw%she0(:,3)*pl%te*dv%ne
        call my_out_us(70,nCv,0,wrk0,'b2stel_she_rad')
        wrk0=shi0_ion(:,0)+shi0_ion(:,1)*pl%ti+
     &       shi0_ion(:,2)*dv%ni(:,0)+
     &       shi0_ion(:,3)*pl%ti*dv%ni(:,0)
        call my_out_us(70,nCv,0,wrk0,'b2stel_shi_ion')                  !srv 01.11.10
        wrk0=shi0_rec(:,0)+shi0_rec(:,1)*pl%ti+
     &       shi0_rec(:,2)*dv%ni(:,0)+
     &       shi0_rec(:,3)*pl%ti*dv%ni(:,0)
        call my_out_us(70,nCv,0,wrk0,'b2stel_shi_rec')                  !srv 01.11.10

!WG_TODO        do is=0,ns-1                                                 !srv 03.02.12 {
!WG_TODO          write(chns,'(i3.3)') is
!WG_TODO          call my_out_us(70,nCv,0,rqrad(1,is),'b2stel_rqrad'//chns)
!WG_TODO          call my_out_us(70,nCv,0,rqbrm(1,is),'b2stel_rqbrm'//chns)
!WG_TODO        enddo                                                        !srv 03.02.12 }
      endif                                                          !srv 30.06.08 11.09.09 }

cdjm Jan2017 Keep linearised sources for balance
#ifdef WG_TODO
      if (balance_netcdf.ne.0) then
        b2stel_sna_ion0to1 = sna0_ion
        b2stel_sna_rec0to1 = sna0_rec
        b2stel_smq_ion0to3 = smq0_ion
        b2stel_smq_rec0to3 = smq0_rec
        b2stel_she0to3=srw%she0
        b2stel_shi_ion0to3=shi0_ion
        b2stel_shi_rec0to3=shi0_rec
      endif
#endif

* ..return
      ncall = ncall+1
      call subend ()
      return

*-----------------------------------------------------------------------
*.scribbles

*!!   Consider to calculate sna0, she0, shi0 over (-1:nx,-1:ny)
*     instead of over (0:nx-1,0:ny-1).
*xpb  Went back to using just (0:nx-1,0:ny-1). Otherwise sources (in
*xpb  particular recombination sources) in guard cells are not balanced
*xpb  by anything and this can lead to particle accumulation and
*xpb  temperature collapse.

*-----------------------------------------------------------------------
*.end b2stel

      end subroutine b2stel
