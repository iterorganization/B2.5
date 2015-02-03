*-----------------------------------------------------------------------
*.specification

      subroutine b2stel (nx, ny, ns, ismain,                             !srv 28.07.11
     & BoRiS, bb, vol, hx, hy, hz, na, ua,                               !srv 09.01.01
     & wadia, vadia, vaecrb, te, ti, po, ne, ni, facdrift, fac_ExB,      !srv 07.07.09
     & rsa, rra, rqa, rrd, rbr, rza, rz2, rpt, rpi, rlsa, rlra, rlqa, 
     & sna0_ion, sna0_rec, smq0_ion, smq0_rec, she0, shi0_ion, shi0_rec) !srv 11.09.09
      use b2mod_types
      use b2mod_diag
      use b2mod_tallies
      use b2mod_indirect
      use b2mod_constants
      use b2mod_sources
      use b2mod_eirene_globals
      use b2mod_b2cmpa
      use b2mod_b2cmrc
#ifdef ASTRA
      use b2mod_b2_to_astra
     , , only : svol, se, seh, sih, b2_to_astra_src  !srv 28.07.11 11.06.13
#endif	                       
      implicit none
*   ..input arguments (unchanged on exit)
      integer nx, ny, ns, ismain
      real (kind=R8) :: 
     *  BoRiS, facdrift(-1:nx,-1:ny), fac_ExB(-1:nx,-1:ny), 
     &  vol(-1:nx,-1:ny), na(-1:nx,-1:ny,0:ns-1), hz(-1:nx,-1:ny),    !srv 09.01.01
     &  bb(-1:nx,-1:ny,0:3), hx(-1:nx,-1:ny), hy(-1:nx,-1:ny),
     &  ua(-1:nx,-1:ny,0:ns-1), te(-1:nx,-1:ny), ti(-1:nx,-1:ny),
     &  ne(-1:nx,-1:ny), ni(-1:nx,-1:ny,0:1), po(-1:nx,-1:ny),
     &  wadia(-1:nx,-1:ny,0:1,0:ns-1), vaecrb(-1:nx,-1:ny,0:1,0:ns-1),
     &  vadia(-1:nx,-1:ny,0:1,0:ns-1)                                 !srv 07.07.09
      real (kind=R8) ::
     *  rsa(-1:nx,-1:ny,0:ns-1), rra(-1:nx,-1:ny,0:ns-1),
     &  rqa(-1:nx,-1:ny,0:ns-1), rrd(-1:nx,-1:ny,0:ns-1), 
     &  rbr(-1:nx,-1:ny,0:ns-1), 
     &  rza(-1:nx,-1:ny,0:ns-1), rz2(-1:nx,-1:ny,0:ns-1), 
     &  rpt(-1:nx,-1:ny,0:ns-1), rpi(-1:nx,-1:ny,0:ns-1),
     &  rlsa(-1:nx,-1:ny,0:1,0:ns-1),
     &  rlra(-1:nx,-1:ny,0:1,0:ns-1), rlqa(-1:nx,-1:ny,0:1,0:ns-1)
*   ..output arguments (unspecified on entry)
      real (kind=R8) ::
     *  sna0_ion(-1:nx,-1:ny,0:1,0:ns-1),                            !srv 11.09.09
     &  sna0_rec(-1:nx,-1:ny,0:1,0:ns-1),
     *  smq0_ion(-1:nx,-1:ny,0:3,0:ns-1),
     *  smq0_rec(-1:nx,-1:ny,0:3,0:ns-1),
     &  she0(-1:nx,-1:ny,0:3),
     &  shi0_ion(-1:nx,-1:ny,0:3), shi0_rec(-1:nx,-1:ny,0:3)
*   ..common blocks
*-----------------------------------------------------------------------
*.documentation
*
*  1. purpose
*
*     B2STEL provides linearised expressions for the sources of
*     density, momentum, electron heat and atom heat, due to the
*     electron-atom collision processes of ionisation, recombination,
*     and radiation.
*
*-----------------------------------------------------------------------
*.declarations

*   ..local variables
      integer ncall, styl0, n2, ix, iy, is, js, fix_recomb_energy, k
      integer, save :: nimp_tmp, jimp_tmp(natm)
	integer, save :: iout = 0,
     &                 iout_b2npc9 = 0, iout_b2npco = 0,
     &                 iout_b2nph9 = 0, iout_b2npht = 0,
     &                 iout_b2npmo =0
	character chns*2, chk*1                                           !srv 30.06.08
      real (kind=R8) ::
     *  rg0, rg1, rxm0, rxm1, rxm2, rxm3, rxm4, phm0
      save ncall, styl0, rg0, rg1, rxm0, rxm1, rxm2, rxm3, rxm4, phm0,
     *  fix_recomb_energy
      real (kind=R8) ::
     *  rf0, t0, t1
*   ..procedures
      intrinsic abs, max
      logical lnext
      real (kind=R8) ::
     *   samax, roxa
      external subini, subend, xertst, ipgeti, ipgetr, samax, sfill
      external b2xvsg
      lnext(is,js) = zamax(is).lt.zamin(js).and.zn(is).eq.zn(js).and.
     &  abs(am(is)-am(js)).lt.1.0e-3_R8
      roxa(ix,iy,is) = am(is)*mp*na(ix,iy,is)
*   ..initialisation
      data ncall/0/, styl0/0/, rg0/1.0_R8/, 
     &  rg1/1.0_R8/, rxm0/0.0_R8/,
     &  rxm1/0.0_R8/, rxm2/0.0_R8/,
     &  rxm3/0.0_R8/, rxm4/0.0_R8/, phm0/0.0_R8/,
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
       nimp_tmp=0
       jimp_tmp=0
       do is = 1, nspecies !{
         if(nfluids(is).gt.1) then !{
           nimp_tmp=nimp_tmp+1
           jimp_tmp(is)=nimp_tmp
         end if !}
       end do !}
       call ipgeti ('b2stel_iout', iout)                         !srv 30.06.08
       call ipgeti ('b2npc9_iout', iout_b2npc9)                  !srv 11.09.09 {
       call ipgeti ('b2npco_iout', iout_b2npco)
       call ipgeti ('b2npmo_iout', iout_b2npmo)
       call ipgeti ('b2nph9_iout', iout_b2nph9)
       call ipgeti ('b2npht_iout', iout_b2npht)                  !srv 11.09.09 }
      endif
*   ..set n2
      n2 = (nx+2)*(ny+2)
*   ..test nx, ny, ns
      call xertst (0.le.nx.and.0.le.ny, 'faulty argument nx, ny')
      call xertst (1.le.ns, 'faulty argument ns')
*   ..extensive tests on first few calls
      if (ncall.lt.3) then
*    ..test sign of bb
       call b2xvsg (n2, bb(-1,-1,3), 1, 'bb3', '.gt.')
*    ..test vol, hx, hy, qz
       call b2xvsg (n2, vol, 1, 'vol', '.gt.')
       call b2xvsg (n2, hx, 1, 'hx', '.gt.')
       call b2xvsg (n2, hy, 1, 'hy', '.gt.')
*    ..test sign of na, te, ti, ne, ni
       call b2xvsg (n2*ns, na, 1, 'na', '.gt.')
       call b2xvsg (n2, te, 1, 'te', '.gt.')
       call b2xvsg (n2, ti, 1, 'ti', '.gt.')
       call b2xvsg (n2, ne, 1, 'ne', '.gt.')
       call b2xvsg (n2*2, ni, 1, 'ni', '.gt.')
*   ..test facdrift, fac_ExB
      call xertst (0.0_R8.le.minval(facdrift) .and.
     &     1.0_R8.ge.maxval(facdrift), 'faulty argument facdrift')
      call xertst (0.0_R8.le.minval(fac_ExB) .and.
     &     1.0_R8.ge.maxval(fac_ExB), 'faulty argument fac_ExB')
*    ..test rate coefficients
       do is = 0, ns-2
        if (.not.lnext(is,is+1)) then
         call xertst (samax(n2,rsa(-1,-1,is),1).eq.0.0_R8,
     &     'faulty argument rsa: non-zero for unrelated species')
         call xertst (samax(n2,rra(-1,-1,is+1),1).eq.0.0_R8,
     &     'faulty argument rra: non-zero for unrelated species')
        endif
       enddo
      endif

* ..compute source terms
*   ..initialise sources to 0
      call sfill (n2*2*ns, 0.0_R8, sna0_ion, 1)                      !srv 11.09.09
      call sfill (n2*2*ns, 0.0_R8, sna0_rec, 1)
      call sfill (n2*4*ns, 0.0_R8, smq0_ion, 1)
      call sfill (n2*4*ns, 0.0_R8, smq0_rec, 1)
      call sfill (n2*4, 0.0_R8, she0, 1)
      call sfill (n2*4, 0.0_R8, shi0_ion, 1)
      call sfill (n2*4, 0.0_R8, shi0_rec, 1)
*   ..compute sources due to ionisation
      rsanareg=0.0_R8
      rsahireg=0.0_R8
      rsamoreg=0.0_R8
      rsana=0.0_R8
      rsahi=0.0_R8
      rsamo=0.0_R8
      do is = 0, ns-2
       if (lnext(is,is+1)) then
*     ..particle source and heat source
        do iy = -1, ny
         do ix = -1, nx
          if(rightix(ix,iy).ne.nx+1 .and. leftix(ix,iy).ne.-2 .and.
     &       topiy(ix,iy).ne.ny+1 .and. bottomiy(ix,iy).ne.-2 .and.
     &       region(ix,iy,0).ne.0) then
*       ..compute damping coefficient
           rf0 = rxm0*max(0.0_R8,rlsa(ix,iy,1,is))
*       ..compute rate for is->is+1
           t0 = rsa(ix,iy,is)*vol(ix,iy)*ne(ix,iy)
*       ..compute particle source
           sna0_ion(ix,iy,0,is+1) = sna0_ion(ix,iy,0,is+1)+             !srv 11.09.09 {
     &      (1.0_R8+rf0)*t0*na(ix,iy,is)
           sna0_ion(ix,iy,1,is+1) = sna0_ion(ix,iy,1,is+1)-
     &      rf0*t0*na(ix,iy,is)/na(ix,iy,is+1)
           sna0_ion(ix,iy,0,is) = sna0_ion(ix,iy,0,is)+
     &      rf0*t0*na(ix,iy,is)
           sna0_ion(ix,iy,1,is) = sna0_ion(ix,iy,1,is)-(1.0_R8+rf0)*t0  !srv 11.09.09 }
           rsana(ix,iy,is) = t0*na(ix,iy,is)
           rsanareg(region(ix,iy,0),is)=rsanareg(region(ix,iy,0),is)+
     1      t0*na(ix,iy,is)
*       ..compute atom heat source
           t1 = t0*na(ix,iy,is)*(1.0_R8-BoRiS)*
     &       ((am(is)+am(is+1))/2.0_R8)*(mp/2.0_R8)*
     &       ((ua(ix,iy,is)-ua(ix,iy,is+1))**2)
           shi0_ion(ix,iy,0) = shi0_ion(ix,iy,0)+(1.0_R8+rf0)*t1
           shi0_ion(ix,iy,3) = shi0_ion(ix,iy,3)-
     &                         rf0*t1/(ni(ix,iy,0)*ti(ix,iy))
           rsahi(ix,iy,is) = t1
           rsahireg(region(ix,iy,0),is)=rsahireg(region(ix,iy,0),is)+t1
          endif
         enddo
        enddo
*     ..parallel momentum source
        do iy = -1, ny
          do ix = -1, nx
            if(rightix(ix,iy).ne.nx+1 .and. leftix(ix,iy).ne.-2 .and.
     &         topiy(ix,iy).ne.ny+1 .and. bottomiy(ix,iy).ne.-2 .and.
     &         region(ix,iy,0).ne.0) then
*       ..compute damping coefficient
              rf0 = rg0+rxm1*max(0.0_R8,rlsa(ix,iy,1,is))
*       ..compute rate for is->is+1
              t0 = rsa(ix,iy,is)*vol(ix,iy)*ne(ix,iy)*na(ix,iy,is)
     &            *hz(ix,iy)                                         !srv 09.01.01
*       ..compute momentum source
              t1 = ((am(is)+am(is+1))/2.0_R8)*mp*t0
              if (styl0.eq.0) then
                smq0_ion(ix,iy,0,is+1) = smq0_ion(ix,iy,0,is+1)+
     &           t1*(ua(ix,iy,is)+rf0*ua(ix,iy,is+1))
                smq0_ion(ix,iy,1,is+1) = smq0_ion(ix,iy,1,is+1)-
     &           rf0*t1
                smq0_ion(ix,iy,0,is) = smq0_ion(ix,iy,0,is)+
     &           rf0*t1*ua(ix,iy,is)
                smq0_ion(ix,iy,1,is) = smq0_ion(ix,iy,1,is)-
     &           (1.0_R8+rf0)*t1
              else if (styl0.eq.1) then
*       ..compute momentum source
               smq0_ion(ix,iy,2,is+1) = smq0_ion(ix,iy,2,is+1)+
     &          t1*(ua(ix,iy,is)+rf0*ua(ix,iy,is+1))/roxa(ix,iy,is+1)
               smq0_ion(ix,iy,3,is+1) = smq0_ion(ix,iy,3,is+1)-
     &          rf0*t1/roxa(ix,iy,is+1)
               smq0_ion(ix,iy,2,is) = smq0_ion(ix,iy,2,is)+
     &          rf0*t1*ua(ix,iy,is)/roxa(ix,iy,is)
               smq0_ion(ix,iy,3,is) = smq0_ion(ix,iy,3,is)-
     &          (1.0_R8+rf0)*t1/roxa(ix,iy,is)
             else if (styl0.eq.2) then
               smq0_ion(ix,iy,0,is+1) = smq0_ion(ix,iy,0,is+1)+
     &          t1*(ua(ix,iy,is)+rf0*ua(ix,iy,is+1))
               smq0_ion(ix,iy,3,is+1) = smq0_ion(ix,iy,3,is+1)-
     &          rf0*t1/roxa(ix,iy,is+1)
               smq0_ion(ix,iy,0,is) = smq0_ion(ix,iy,0,is)+
     &          rf0*t1*ua(ix,iy,is)
               smq0_ion(ix,iy,3,is) = smq0_ion(ix,iy,3,is)-
     &          (1.0_R8+rf0)*t1/roxa(ix,iy,is)
             else if (styl0.eq.3) then
               smq0_ion(ix,iy,2,is+1) = smq0_ion(ix,iy,2,is+1)+
     &          t1*(ua(ix,iy,is)+rf0*ua(ix,iy,is+1))/roxa(ix,iy,is+1)
               smq0_ion(ix,iy,1,is+1) = smq0_ion(ix,iy,1,is+1)-
     &          rf0*t1
               smq0_ion(ix,iy,2,is) = smq0_ion(ix,iy,2,is)+
     &          rf0*t1*ua(ix,iy,is)/roxa(ix,iy,is)
               smq0_ion(ix,iy,1,is) = smq0_ion(ix,iy,1,is)-
     &         (1.0_R8+rf0)*t1
             endif
             rsamo(ix,iy,is) = rsamo(ix,iy,is)+t1*ua(ix,iy,is)
             rsamoreg(region(ix,iy,0),is) = 
     &        rsamoreg(region(ix,iy,0),is)+t1*ua(ix,iy,is)
             rsamo(ix,iy,is+1) = rsamo(ix,iy,is+1)-t1*ua(ix,iy,is)
             rsamoreg(region(ix,iy,0),is+1) = 
     &        rsamoreg(region(ix,iy,0),is+1)-t1*ua(ix,iy,is)
           endif
          enddo
        enddo
       endif
      enddo
*   ..compute sources due to recombination
      rranareg=0.0_R8
      rrahireg=0.0_R8
      rramoreg=0.0_R8
      rqahereg=0.0_R8
      rrana=0.0_R8
      rrahi=0.0_R8
      rramo=0.0_R8
      rqahe=0.0_R8
      do is = 1, ns-1
       if (lnext(is-1,is)) then
*     ..particle source and heat source
        do iy = -1, ny
         do ix = -1, nx
          if(rightix(ix,iy).ne.nx+1 .and. leftix(ix,iy).ne.-2 .and.
     &       topiy(ix,iy).ne.ny+1 .and. bottomiy(ix,iy).ne.-2 .and.
     &       region(ix,iy,0).ne.0) then
*       ..compute damping coefficient
           rf0 = rxm2*max(0.0_R8,-rlra(ix,iy,1,is))
*       ..compute rate for is->is-1
           t0 = rra(ix,iy,is)*vol(ix,iy)*ne(ix,iy)
*       ..compute particle source
           sna0_rec(ix,iy,0,is) = sna0_rec(ix,iy,0,is)+               !srv 11.09.09 {
     &       rf0*t0*na(ix,iy,is)
           sna0_rec(ix,iy,1,is) = sna0_rec(ix,iy,1,is)-(1.0_R8+rf0)*t0
           sna0_rec(ix,iy,0,is-1) = sna0_rec(ix,iy,0,is-1)+
     &      (1.0_R8+rf0)*t0*na(ix,iy,is)
           sna0_rec(ix,iy,1,is-1) = sna0_rec(ix,iy,1,is-1)-           !srv 11.09.09 {
     &      rf0*t0*na(ix,iy,is)/na(ix,iy,is-1)
           rrana(ix,iy,is) = t0*na(ix,iy,is)
           rranareg(region(ix,iy,0),is)=rranareg(region(ix,iy,0),is)+
     1      t0*na(ix,iy,is)
*       ..compute atom heat source
           t1 = t0*na(ix,iy,is)*(1.0_R8-BoRiS)*
     &       ((am(is)+am(is-1))/2.0_R8)*(mp/2.0_R8)*
     &       ((ua(ix,iy,is)-ua(ix,iy,is-1))**2)
           shi0_rec(ix,iy,0) = shi0_rec(ix,iy,0)+(1.0_R8+rf0)*t1
           shi0_rec(ix,iy,3) = shi0_rec(ix,iy,3)-
     *                         rf0*t1/(ni(ix,iy,0)*ti(ix,iy))
           rrahi(ix,iy,is) = t1
           rrahireg(region(ix,iy,0),is)=rrahireg(region(ix,iy,0),is)+t1
*       ..compute electron heat source for is->is-1, if not already included from ADPAK
           if(fix_recomb_energy.ne.0) then 
c the default treatment is to have the entire ionization potential 
c difference lost, and the only additional cooling rate correction 
c to be the Bremsstrahlung 
c this treatment changes this --- the ionization potential difference 
c is given back as energy to the electrons and additional radiative 
c losses should be included (in ADAS terms, the PRB file should be 
c accounted for) 
             she0(ix,iy,0)=she0(ix,iy,0)+(1.0_R8+rf0)*
     &        t0*rpi(ix,iy,is)*ev*na(ix,iy,is) 
             she0(ix,iy,3)=she0(ix,iy,0)-rf0*
     &        t0*rpi(ix,iy,is)*ev*na(ix,iy,is)/(ne(ix,iy)*te(ix,iy))
             rqahe(ix,iy,is) = rqahe(ix,iy,is) - 
     &        t0*rpi(ix,iy,is)*ev*na(ix,iy,is)
             rqahereg(region(ix,iy,0),is) = 
     &        rqahereg(region(ix,iy,0),is) -
     &        t0*rpi(ix,iy,is)*ev*na(ix,iy,is)
           endif 
          endif
         enddo
        enddo
*     ..parallel momentum source
        do iy = -1, ny
         do ix = -1, nx
          if(rightix(ix,iy).ne.nx+1 .and. leftix(ix,iy).ne.-2 .and.
     &       topiy(ix,iy).ne.ny+1 .and. bottomiy(ix,iy).ne.-2 .and.
     &       region(ix,iy,0).ne.0) then
*       ..compute damping coefficient
           rf0 = rg1+rxm3*max(0.0_R8,-rlra(ix,iy,1,is))
*       ..compute rate for is->is-1
           t0 =rra(ix,iy,is)*vol(ix,iy)*ne(ix,iy)*na(ix,iy,is)
     &         *hz(ix,iy)                                            !srv 09.01.01
*       ..compute momentum source
           t1 = ((am(is)+am(is-1))/2.0_R8)*mp*t0
           smq0_rec(ix,iy,0,is) = smq0_rec(ix,iy,0,is)+
     &      rf0*t1*ua(ix,iy,is)
           smq0_rec(ix,iy,1,is) = smq0_rec(ix,iy,1,is)-(1.0_R8+rf0)*t1
           smq0_rec(ix,iy,0,is-1) = smq0_rec(ix,iy,0,is-1)+
     &      t1*(ua(ix,iy,is)+rf0*ua(ix,iy,is-1))
           smq0_rec(ix,iy,1,is-1) = smq0_rec(ix,iy,1,is-1)-rf0*t1
           rramo(ix,iy,is) = rramo(ix,iy,is)+t1*ua(ix,iy,is)
           rramoreg(region(ix,iy,0),is) = 
     &      rramoreg(region(ix,iy,0),is)+t1*ua(ix,iy,is)
           rramo(ix,iy,is-1) = rramo(ix,iy,is-1)+t1*ua(ix,iy,is)
           rramoreg(region(ix,iy,0),is-1) = 
     &      rramoreg(region(ix,iy,0),is-1)+t1*ua(ix,iy,is)
          endif
         enddo
        enddo
       endif
      enddo
*   ..compute source terms for electron heat loss and radiation
      b2rad=0.0_R8
      b2radreg=0.0_R8
      b2bremreg=0.0_R8
      rqradreg=0.0_R8
      rqbrmreg=0.0_R8
      rqrad=0.0_R8
      b2brem=0.0_R8
      rad_imp=0.0_R8
      do is = 0, ns-1
       do iy = -1, ny
        do ix = -1, nx
         if(rightix(ix,iy).ne.nx+1 .and. leftix(ix,iy).ne.-2 .and.
     &      topiy(ix,iy).ne.ny+1 .and. bottomiy(ix,iy).ne.-2 .and.
     &      region(ix,iy,0).ne.0) then
*      ..compute damping coefficient
          rf0 = rxm4*max(0.0_R8,rlqa(ix,iy,1,is))
*      ..compute electron cooling rate for is->any
*      ..includes line radiation, bremsstrahlung and ionisation contributions
*      ..if adpak_used, also contains recombination and CX terms
          t0 = rqa(ix,iy,is)*vol(ix,iy)*ne(ix,iy)*na(ix,iy,is)
*      ..compute electron heat sink
          she0(ix,iy,0) = she0(ix,iy,0)+rf0*t0
          she0(ix,iy,3) = she0(ix,iy,3)-
     &     (1.0_R8+rf0)*t0/(ne(ix,iy)*te(ix,iy))
          rqahe(ix,iy,is) = rqahe(ix,iy,is) + t0
          rqahereg(region(ix,iy,0),is)=rqahereg(region(ix,iy,0),is)+t0
*      ..compute line radiation rate for is->any
          t0 = rrd(ix,iy,is)*vol(ix,iy)*ne(ix,iy)*na(ix,iy,is)
*      ..compute line radiation
          rqrad(ix,iy,is) = t0
          rqradreg(region(ix,iy,0),is)=rqradreg(region(ix,iy,0),is)+t0
          b2radreg(region(ix,iy,0))=b2radreg(region(ix,iy,0))-t0
          b2rad=b2rad-t0
          if (.not.is_neutral(is) .and. jimp_tmp(b2espcr(is)).ne.0) then
            rad_imp(jimp_tmp(b2espcr(is)))=
     .       rad_imp(jimp_tmp(b2espcr(is)))-t0
          endif
*      ..compute bremsstrahlung radiation rate for is->any
          t0 = rbr(ix,iy,is)*vol(ix,iy)*ne(ix,iy)*na(ix,iy,is)
*      ..compute bremsstrahlung radiation
          rqbrm(ix,iy,is) = t0
          rqbrmreg(region(ix,iy,0),is)=rqbrmreg(region(ix,iy,0),is)+t0
          b2bremreg(region(ix,iy,0))=b2bremreg(region(ix,iy,0))-t0
          b2brem=b2brem-t0
         endif
        enddo
       enddo
      enddo

#ifdef ASTRA

      do is = 0, ns-1                                                 !srv 28.07.11 {
	 if(is.eq.ismain) then
        call b2_to_astra_src(nx, ny, vol, sna0_ion(:,:,0,is)+
     &   sna0_ion(:,:,1,is)*na(:,:,is), svol, se)
        call b2_to_astra_src(nx, ny, vol, sna0_rec(:,:,0,is)+
     &   sna0_rec(:,:,1,is)*na(:,:,is), svol, se)
	 endif
	enddo                                                           !srv 28.07.11 }

      call b2_to_astra_src(nx, ny, vol, she0(:,:,0)+
     &   she0(:,:,1)*te(:,:)+she0(:,:,2)*ne+
     &   she0(:,:,3)*ne*te, svol, seh)                            !srv 11.06.13

      call b2_to_astra_src(nx, ny, vol, shi0_ion(:,:,0)+
     &   shi0_ion(:,:,1)*ti+shi0_ion(:,:,2)*ni(:,:,0)+
     &   shi0_ion(:,:,3)*ni(:,:,0)*ti, svol, sih)                     !srv 11.06.13

      call b2_to_astra_src(nx, ny, vol, shi0_rec(:,:,0)+
     &   shi0_rec(:,:,1)*ti+shi0_rec(:,:,2)*ni(:,:,0)+
     &   shi0_rec(:,:,3)*ni(:,:,0)*ti, svol, sih)                     !srv 11.06.13
#endif	                       

	if (iout.eq.1 .or. iout_b2npmo.eq.1) then                       !srv 30.06.08 11.09.09 {
       do is = 0, ns-1
	  write (chns,'(i2.2)') is
!	  do k=0,3
!	   write (chk,'(i1)') k
!	   call my_out(70,nx,ny,smq0_ion(-1,-1,k,is),
!     &                       'b2stel_smq_ion'//chk//chns)
!	   call my_out(70,nx,ny,smq0_rec(-1,-1,k,is),
!     &                       'b2stel_smq_rec'//chk//chns)
!	  enddo
	 call my_out(70,nx,ny,smq0_ion(:,:,0,is)+
     &        smq0_ion(:,:,1,is)*ua(:,:,is)+
     &        smq0_ion(:,:,2,is)*na(:,:,is)*mp*am(is)+
     &        smq0_ion(:,:,3,is)*ua(:,:,is)*na(:,:,is)*mp*am(is),
     &        'b2stel_smq_ion'//chns)
	 call my_out(70,nx,ny,smq0_rec(:,:,0,is)+
     &        smq0_rec(:,:,1,is)*ua(:,:,is)+
     &        smq0_rec(:,:,2,is)*na(:,:,is)*mp*am(is)+
     &        smq0_rec(:,:,3,is)*ua(:,:,is)*na(:,:,is)*mp*am(is),
     &        'b2stel_smq_rec'//chns)
	 enddo
	endif
	if (iout.eq.1 .or. (iout_b2npc9.eq.1 .or. iout_b2npco.eq.1)) then
       do is = 0, ns-1
	  write (chns,'(i2.2)') is
!	  do k=0,1
!	   write (chk,'(i1)') k
!	   call my_out(70,nx,ny,sna0_ion(-1,-1,k,is),
!     &                       'b2stel_sna_ion'//chk//chns)
!	   call my_out(70,nx,ny,sna0_rec(-1,-1,k,is),
!     &                       'b2stel_sna_rec'//chk//chns)
!	  enddo
	  call my_out(70,nx,ny,sna0_ion(:,:,0,is)+
     &             sna0_ion(:,:,1,is)*na(:,:,is),'b2stel_sna_ion'//chns)
	  call my_out(70,nx,ny,sna0_rec(:,:,0,is)+
     &             sna0_rec(:,:,1,is)*na(:,:,is),'b2stel_sna_rec'//chns)

	 enddo
	endif
	if (iout.eq.1 .or. (iout_b2nph9.eq.1 .or. iout_b2npht.eq.1)) then
!	 do k=0,3
!	  write (chk,'(i1)') k
!	  call my_out(70,nx,ny,shi0_ion(-1,-1,k),'b2stel_shi_ion'//chk)
!	  call my_out(70,nx,ny,shi0_rec(-1,-1,k),'b2stel_shi_rec'//chk)
!	 enddo
!      call my_out(70,nx,ny,she0(-1,-1,k),'b2stel_she')
	 call my_out(70,nx,ny,she0(:,:,0)+she0(:,:,1)*te+
     &   she0(:,:,2)*ne + she0(:,:,3)*te*ne,'b2stel_she_rad')
	 call my_out(70,nx,ny,shi0_ion(:,:,0)+shi0_ion(:,:,1)*ti+
     &   shi0_ion(:,:,2)*ni(:,:,0)+shi0_ion(:,:,3)*ti*ni(:,:,0),
     &   'b2stel_shi_ion')                                            !srv 01.11.10
	 call my_out(70,nx,ny,shi0_rec(:,:,0)+shi0_rec(:,:,1)*ti+
     &   shi0_rec(:,:,2)*ni(:,:,0)+shi0_rec(:,:,3)*ti*ni(:,:,0),
     &   'b2stel_shi_rec')                                            !srv 01.11.10

       do is=0,ns-1                                                !srv 03.02.12 {
        write(chns,'(i2.2)') is
	  call my_out(70,nx,ny,rqrad(-1,-1,is),'b2stel_rqrad'//chns)
	  call my_out(70,nx,ny,rqbrm(-1,-1,is),'b2stel_rqbrm'//chns)
       enddo                                                       !srv 03.02.12 }
	endif                                                           !srv 30.06.08 11.09.09 }

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
