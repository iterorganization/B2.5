*-----------------------------------------------------------------------
*.specification

      subroutine b2mwqt (nput, nCv, ns, itim, itcnt, BoRiS, geo, pl, dv,
     &  diag)
      use b2mod_types
      use b2mod_constants
      use b2mod_b2cmpa
      use b2us_geo
      use b2us_plasma
      implicit none
*   ..input arguments (unchanged on exit)
      integer nput, nCv, ns, itim, itcnt
      real (kind=R8) :: BoRiS
      type (geometry), intent (in) :: geo
      type (B2Plasma), intent (in) :: pl
      type (B2Derivatives), intent (in) :: dv
      type (B2Diagnostic), intent (in) :: diag
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
*     areshn - real, input.
*     areshn specifies the scaled norm of the residual of the fluid neutral
*     heat balance equation. Only calculated if a separate neutral energy 
*     equation is solved. In this case areshi contains only the ions.
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
      integer ncall, nfldmx, style, nfld, add_areshn
      parameter (nfldmx=1000)
      save ncall, style, add_areshn
      integer k, icv, is, iicv, il, ix0(-1:5), iy0(-1:5)
      real (kind=R8) ::
     *   z(0:nfldmx-1), t0, t1, tv
*   ..procedures
      intrinsic min, max
      real (kind=R8) ::
     *   smin, b2ssum
      external subini, subend, xertst, ipgeti, smin, b2ssum, cfwure
*   ..initialisation
      data ncall/0/, style/1/, add_areshn/0/

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
       call ipgeti ('b2mwqt_add_areshn', add_areshn)
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
      call xertst (0.0_R8.le.smin(ns,diag%aresco,1).and.
     &  0.0_R8.le.diag%areshe.and.0.0_R8.le.diag%areshi.and.
     &  0.0_R8.le.diag%arespo,'faulty argument aresco/ee/ei/po')
*   ..test corrections
      call xertst (0.0_R8.le.smin(ns,diag%acorpa,1).and.
     &  0.0_R8.le.diag%acorte.and.0.0_R8.le.diag%acorti.and.
     &  0.0_R8.le.diag%acorpo,'faulty argument acorpa/te/ti/po')
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
      if (add_areshn.eq.1) then
        nfld = nfld + 1
      endif
      call xertst (nfld.le.nfldmx, 'faulty parameter nfldmx')
*   ..define region boundaries
#ifdef WG_TODO
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
#endif
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
        z(k) = diag%aresco(is)
        k = k+1
       enddo
       do is = 0, ns-1
        z(k) = diag%aresmo(is)
        k = k+1
       enddo
       z(k) = diag%aresmt
       k = k+1
       z(k) = diag%areshe
       k = k+1
       z(k) = diag%areshi
       k = k+1
       if (add_areshn.eq.1) then
         z(k) = diag%areshn
         k = k+1
       endif
       z(k) = diag%arespo
       k = k+1
*    ..norms of corrections
       do is = 0, ns-1
        z(k) = diag%acorpa(is)
        k = k+1
       enddo
       do is = 0, ns-1
        z(k) = diag%acorua(is)
        k = k+1
       enddo
       z(k) = diag%acorut
       k = k+1
       z(k) = diag%acorte/ev
       k = k+1
       z(k) = diag%acorti/ev
       k = k+1
       z(k) = diag%acorpo
       k = k+1
      endif
      if (2.le.style) then
*    ..total content on regions
       call xerrab ('WG - to be done!')
#ifdef WG_TODO
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
#endif
      endif
      if (3.le.style) then
*    ..integrated fluxes
       call xerrab ('WG - to be done!')
#ifdef WG_TODO
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
#endif
      endif
      if (4.le.style) then
*    ..min and max on regions
       call xerrab ('WG - to be done!')
#ifdef WG_TODO
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
#endif
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
