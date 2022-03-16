*-----------------------------------------------------------------------
*.specification

      subroutine b2mwqt (nput, nx, ny, ns, itim, itcnt, BoRiS, vol,
     &  na, ua, te, ti, po, ne, ni, fna, fhe, fhi, fch,
     &  aresco, aresmo, aresmt, areshe, areshi, arespo,
     &  acorpa, acorua, acorut, acorte, acorti, acorpo)
      use b2mod_types
      use b2mod_constants
      use b2mod_b2cmpa
      use b2mod_subsys
      implicit none
*   ..input arguments (unchanged on exit)
      integer nput, nx, ny, ns, itim, itcnt
      real (kind=R8) :: BoRiS,
     *  na(-1:nx,-1:ny,0:ns-1), ua(-1:nx,-1:ny,0:ns-1),
     &  te(-1:nx,-1:ny), ti(-1:nx,-1:ny), po(-1:nx,-1:ny),
     &  ne(-1:nx,-1:ny), ni(-1:nx,-1:ny,0:1), vol(-1:nx,-1:ny)
      real (kind=R8) ::
     *  fna(-1:nx,-1:ny,0:1,0:1,0:ns-1), fhe(-1:nx,-1:ny,0:1,0:1),
     &  fhi(-1:nx,-1:ny,0:1,0:1), fch(-1:nx,-1:ny,0:1,0:1)
      real (kind=R8) ::
     *  aresco(0:ns-1), aresmo(0:ns-1), aresmt, areshe, areshi,
     &  arespo, acorpa(0:ns-1), acorua(0:ns-1), acorut, acorte,
     &  acorti, acorpo
*   ..common blocks
*-----------------------------------------------------------------------
*.documentation
*
*  1. purpose
*
*     B2MWQT writes brief data concerning the numerical evolution.
*
*
*  3. description (see also routine b2cdca)
*
*     This routine writes to unit (nput) the 'quick trace' evolution
*     data. Unit (nput) must be connected to an un*formatted file.
*
*
*  5. parameters (see also routine b2cdcv)
*
*     nput - integer, input.
*     nput specifies a unit number for un*formatted output.
*
*     ns - integer, input.
*     ns specifies the number of atomic species in the calculation.
*     The species are indexed by (0:ns-1).
*     It will hold that 1.le.ns.
*
*     itim - integer, input.
*     itim specifies the timestep counter. It is initialised to 0 and
*     incremented by 1 following each successful call to b2mndt.
*
*     itcnt - integer, input.
*     itcnt counts the number of inner iterations. It is initialised to
*     0 and incremented by 1 following each call to b2news.
*
*     na, ua, te, ti, po, ne, ni - real*8 array, input.
*     See the calling routine for description.
*
*     fna, fhe, fhi, fch - real*8 array, input.
*     See the calling routine for description.
*
*     aresco - (0:ns-1) real*8 array, input.
*     For (is) in (0:ns-1), aresco(is) specifies the scaled norm of the
*     residual of the continuity equation for species (is).
*
*     aresmo - (0:ns-1) real*8 array, input.
*     For (is) in (0:ns-1), aresmo(is) specifies the scaled norm of the
*     residual of the parallel momentum equation for species (is).
*
*     aresmt - real, input.
*     aresmt specifies the scaled norm of the residual of the total
*     parallel momentum equation.
*
*     areshe - real, input.
*     areshe specifies the scaled norm of the residual of the electron
*     heat balance equation.
*
*     areshi - real, input.
*     areshi specifies the scaled norm of the residual of the all atom
*     heat balance equation.
*
*     arespo - real, input.
*     arespo specifies the scaled norm of the residual of the electric
*     potential equation.
*
*     acorpa - (0:ns-1) real*8 array, input.
*     For (is) in (0:ns-1), acorpa(is) specifies the norm of the last
*     computed correction for the continuity equation for species (is).
*
*     acorua - (0:ns-1) real*8 array, input.
*     For (is) in (0:ns-1), acorua(is) specifies the norm of the last
*     computed correction for the parallel momentum balance equation for
*     species (is).
*
*     acorut - real, input.
*     acorut specifies the norm of the last computed correction for the
*     parallel momentum balance equation for all species simultaneously.
*
*     acorte - real, input.
*     acorte specifies the norm of the last computed correction for the
*     electron heat balance equation.
*
*     acorti - real, input.
*     acorti specifies the norm of the last computed correction for the
*     all atom heat balance equation.
*
*     acorpo - real, input.
*     acorpo specifies the norm of the last computed correction for the
*     electric potential equation.
*
*     (parameters in common:)
*
*     /b2cmpa/ - input.
*     /b2cmpa/ specifies basic physics parameters.
*     See the calling routine for description.
*
*-----------------------------------------------------------------------
*.declarations

*   ..local variables
      integer ncall, nfldmx, style, nfld
      parameter (nfldmx=1000)
      save ncall, style
      integer k, ix, iy, is, iix, iiy, il, ix0(-1:5), iy0(-1:5)
      real (kind=R8) ::
     *   z(0:nfldmx-1), t0, t1, tv
*   ..procedures
      intrinsic min, max
      real (kind=R8) ::
     *   smin, b2ssum
      external xertst, ipgeti, smin, b2ssum, cfwure
*   ..initialisation
      data ncall/0/, style/1/

*-----------------------------------------------------------------------
*.documentation-internal
*
*     Description of some local variables:
*
*     style - integer.
*     style specifies the amount of output produced on each call.
*     Legitimate values are in the range 0..4. style.eq.0 implies
*     no output, and the amount of output is an ascending function
*     of the value of style. See the code for precise information.

*-----------------------------------------------------------------------
*.computation

* ..preliminaries
*   ..subprogram start-up calls
      call subini ('b2mwqt')
*   ..set internal parameters on first call
      if (ncall.eq.0) then
       call ipgeti ('b2mwqt_style', style)
       call xertst (0.le.style, 'faulty internal parameter style')
      endif
*   ..test nput
      call xertst (1.le.nput, 'faulty argument nput')
*   ..test ns
      call xertst (1.le.ns, 'faulty argument ns')
*   ..test itim, itcnt
      call xertst (0.le.itim.and.0.le.itcnt,
     &  'faulty argument itim, itcnt')
*   ..test residuals
      call xertst (0.0_R8.le.smin(ns,aresco,1).and.
     &  0.0_R8.le.areshe.and.0.0_R8.le.areshi.and.0.0_R8.le.arespo,
     &  'faulty argument aresco/ee/ei/po')
*   ..test corrections
      call xertst (0.0_R8.le.smin(ns,acorpa,1).and.
     &  0.0_R8.le.acorte.and.0.0_R8.le.acorti.and.0.0_R8.le.acorpo,
     &  'faulty argument acorpa/te/ti/po')
*   ..set and test nfld
      nfld = 0
      if (1.le.style) then
       nfld = nfld+2+2*(2*ns+4)
      endif
      if (2.le.style) then
       nfld = nfld+36*(ns+3)
      endif
      if (3.le.style) then
       nfld = nfld+40*(ns+3)
      endif
      if (4.le.style) then
       nfld = nfld+72*(ns+3)
      endif
      call xertst (nfld.le.nfldmx, 'faulty parameter nfldmx')
*   ..define region boundaries
      ix0(-1) = -1
      do iix = 0, 4
       ix0(iix) = (iix*nx)/4
      enddo
      ix0(5) = ix0(4)+1
      iy0(-1) = -1
      do iiy = 0, 4
       iy0(iiy) = (iiy*ny)/4
      enddo
      iy0(5) = iy0(4)+1

* ..get data for output
      k = 0
      if (1.le.style) then
*    ..iteration counters
       z(k) = itim
       k = k+1
       z(k) = itcnt
       k = k+1
*    ..norms of residuals
       do is = 0, ns-1
        z(k) = aresco(is)
        k = k+1
       enddo
       do is = 0, ns-1
        z(k) = aresmo(is)
        k = k+1
       enddo
       z(k) = aresmt
       k = k+1
       z(k) = areshe
       k = k+1
       z(k) = areshi
       k = k+1
       z(k) = arespo
       k = k+1
*    ..norms of corrections
       do is = 0, ns-1
        z(k) = acorpa(is)
        k = k+1
       enddo
       do is = 0, ns-1
        z(k) = acorua(is)
        k = k+1
       enddo
       z(k) = acorut
       k = k+1
       z(k) = acorte/ev
       k = k+1
       z(k) = acorti/ev
       k = k+1
       z(k) = acorpo
       k = k+1
      endif
      if (2.le.style) then
*    ..total content on regions
       do is = 0, ns-1
        do iiy = -1, 4
         do iix = -1, 4
          t0 = 0.0_R8
          tv = 0.0_R8
          do iy = iy0(iiy), iy0(iiy+1)-1
           do ix = ix0(iix), ix0(iix+1)-1
            t0 = t0+na(ix,iy,is)*vol(ix,iy)
            tv = tv+vol(ix,iy)
           enddo
          enddo
          z(k) = t0      ! /tv
          k = k+1
         enddo
        enddo
       enddo
       do iiy = -1, 4
        do iix = -1, 4
         t0 = 0.0_R8
         tv = 0.0_R8
         do iy = iy0(iiy), iy0(iiy+1)-1
          do ix = ix0(iix), ix0(iix+1)-1
           t0 = t0+(1.5_R8+BoRiS)*ne(ix,iy)*te(ix,iy)*vol(ix,iy)
           tv = tv+vol(ix,iy)
          enddo
         enddo
         z(k) = t0      ! /tv
         k = k+1
        enddo
       enddo
       do iiy = -1, 4
        do iix = -1, 4
         t0 = 0.0_R8
         tv = 0.0_R8
         do iy = iy0(iiy), iy0(iiy+1)-1
          do ix = ix0(iix), ix0(iix+1)-1
           t0 = t0+(1.5_R8+BoRiS)*ni(ix,iy,0)*ti(ix,iy)*vol(ix,iy)
           tv = tv+vol(ix,iy)
          enddo
         enddo
         z(k) = t0      ! /tv
         k = k+1
        enddo
       enddo
       do iiy = -1, 4
        do iix = -1, 4
         t0 = 0.0_R8
         do iy = iy0(iiy), iy0(iiy+1)-1
          do ix = ix0(iix), ix0(iix+1)-1
           t0 = t0+qe*ne(ix,iy)*po(ix,iy)
          enddo
         enddo
         z(k) = t0
         k = k+1
        enddo
       enddo
      endif
      if (3.le.style) then
*    ..integrated fluxes
       do is = 0, ns-1
        do iiy = 0, 3
         do iix = 0, 4
          ix = ix0(iix)
          iy = iy0(iiy)
          il = iy0(iiy+1)-iy0(iiy)
          z(k) = b2ssum(il,fna(ix,iy,0,0,is),nx+2)
          k = k+1
         enddo
        enddo
        do iiy = 0, 4
         do iix = 0, 3
          ix = ix0(iix)
          iy = iy0(iiy)
          il = ix0(iix+1)-ix0(iix)
          z(k) = b2ssum(il,fna(ix,iy,1,1,is),1)
          k = k+1
         enddo
        enddo
       enddo
       do iiy = 0, 3
        do iix = 0, 4
         ix = ix0(iix)
         iy = iy0(iiy)
         il = iy0(iiy+1)-iy0(iiy)
         z(k) = b2ssum(il,fhe(ix,iy,0,0),nx+2)
         k = k+1
        enddo
       enddo
       do iiy = 0, 4
        do iix = 0, 3
         ix = ix0(iix)
         iy = iy0(iiy)
         il = ix0(iix+1)-ix0(iix)
         z(k) = b2ssum(il,fhe(ix,iy,1,1),1)
         k = k+1
        enddo
       enddo
       do iiy = 0, 3
        do iix = 0, 4
         ix = ix0(iix)
         iy = iy0(iiy)
         il = iy0(iiy+1)-iy0(iiy)
         z(k) = b2ssum(il,fhi(ix,iy,0,0),nx+2)
         k = k+1
        enddo
       enddo
       do iiy = 0, 4
        do iix = 0, 3
         ix = ix0(iix)
         iy = iy0(iiy)
         il = ix0(iix+1)-ix0(iix)
         z(k) = b2ssum(il,fhi(ix,iy,1,1),1)
         k = k+1
        enddo
       enddo
       do iiy = 0, 3
        do iix = 0, 4
         ix = ix0(iix)
         iy = iy0(iiy)
         il = iy0(iiy+1)-iy0(iiy)
         z(k) = b2ssum(il,fch(ix,iy,0,0),nx+2)
         k = k+1
        enddo
       enddo
       do iiy = 0, 4
        do iix = 0, 3
         ix = ix0(iix)
         iy = iy0(iiy)
         il = ix0(iix+1)-ix0(iix)
         z(k) = b2ssum(il,fch(ix,iy,1,1),1)
         k = k+1
        enddo
       enddo
      endif
      if (4.le.style) then
*    ..min and max on regions
       do is = 0, ns-1
        do iiy = -1, 4
         do iix = -1, 4
          t0 = na(ix0(iix),iy0(iiy),is)
          t1 = na(ix0(iix),iy0(iiy),is)
          do iy = iy0(iiy), iy0(iiy+1)-1
           do ix = ix0(iix), ix0(iix+1)-1
            t0 = min(t0,na(ix,iy,is))
            t1 = max(t1,na(ix,iy,is))
           enddo
          enddo
          z(k) = t0
          k = k+1
          z(k) = t1
          k = k+1
         enddo
        enddo
       enddo
       do iiy = -1, 4
        do iix = -1, 4
         t0 = te(ix0(iix),iy0(iiy))
         t1 = te(ix0(iix),iy0(iiy))
         do iy = iy0(iiy), iy0(iiy+1)-1
          do ix = ix0(iix), ix0(iix+1)-1
           t0 = min(t0,te(ix,iy))
           t1 = max(t1,te(ix,iy))
          enddo
         enddo
         z(k) = t0/ev
         k = k+1
         z(k) = t1/ev
         k = k+1
        enddo
       enddo
       do iiy = -1, 4
        do iix = -1, 4
         t0 = ti(ix0(iix),iy0(iiy))
         t1 = ti(ix0(iix),iy0(iiy))
         do iy = iy0(iiy), iy0(iiy+1)-1
          do ix = ix0(iix), ix0(iix+1)-1
           t0 = min(t0,ti(ix,iy))
           t1 = max(t1,ti(ix,iy))
          enddo
         enddo
         z(k) = t0/ev
         k = k+1
         z(k) = t1/ev
         k = k+1
        enddo
       enddo
       do iiy = -1, 4
        do iix = -1, 4
         t0 = po(ix0(iix),iy0(iiy))
         t1 = po(ix0(iix),iy0(iiy))
         do iy = iy0(iiy), iy0(iiy+1)-1
          do ix = ix0(iix), ix0(iix+1)-1
           t0 = min(t0,po(ix,iy))
           t1 = max(t1,po(ix,iy))
          enddo
         enddo
         z(k) = t0
         k = k+1
         z(k) = t1
         k = k+1
        enddo
       enddo
      endif
*   ..test k
      call xertst (k.eq.nfld, 'faulty parameter nfld')

* ..produce output
      call cfwure (nput, nfld, z, 'data')

* ..return
      ncall = ncall+1
      call subend ()
      return

*-----------------------------------------------------------------------
*.end b2mwqt

      end subroutine b2mwqt
