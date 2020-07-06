*-----------------------------------------------------------------------
*.specification

      subroutine b2mwqt_diffv (nput, nCv, ns, itim, itcnt, BoRiS, geo, 
     &  pl, dv, diag)
      use b2mod_types
      use b2mod_constants
      use b2mod_b2cmpa_diffv
      use b2us_geo_diffv
      use b2us_plasma_diffv
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
* This is a modified version of b2mwqt used in the B2.5 code differentiated
* with TAPENADE in order to writes brief data concerning the numerical 
* evolution of the differentiated equations. The main difference wrt b2mwqt
* is that here no checks on signs of quantities is done
*
*-----------------------------------------------------------------------
*.declarations

*   ..local variables
      integer ncall, nfldmx, style, nfld
      parameter (nfldmx=1000)
      save ncall, style
      integer k, icv, is, iicv, il, ix0(-1:5), iy0(-1:5)
      real (kind=R8) ::
     *   z(0:nfldmx-1), t0, t1, tv
*   ..procedures
      intrinsic min, max
      real (kind=R8) ::
     *   smin, b2ssum
      external subini, subend, xertst, ipgeti, smin, b2ssum, cfwure
*   ..initialisation
      data ncall/0/, style/1/
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

      end subroutine b2mwqt_diffv
