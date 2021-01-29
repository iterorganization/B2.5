*-----------------------------------------------------------------------
*.specification

      subroutine b2mwqt_diff (nput, nCv, ns, itim, itcnt, BoRiS,
     &  switch, geo, pl, dv, diag)
      use b2mod_types
      use b2mod_constants
      use b2mod_b2cmpa_diff
      use b2mod_switches_diff
      use b2us_geo_diff
      use b2us_plasma_diff
      implicit none
*   ..input arguments (unchanged on exit)
      integer nput, nCv, ns, itim, itcnt
      real (kind=R8) :: BoRiS
      type (switches), intent (in) :: switch
      type (geometry), intent (in) :: geo
      type (B2Plasma_diff), intent (in) :: pl
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
      integer nfldmx, nfld
      parameter (nfldmx=1000)
      integer k, icv, is, iicv, il, ix0(-1:5), iy0(-1:5)
      real (kind=R8) ::
     *   z(0:nfldmx-1), t0, t1, tv
*   ..procedures
      intrinsic min, max
      real (kind=R8) ::
     *   smin, b2ssum
      external subini, subend, xertst, ipgeti, smin, b2ssum, cfwure
*.computation

* ..preliminaries
*   ..subprogram start-up calls
      call subini ('b2mwqt')
*   ..test nput
      call xertst (1.le.nput, 'faulty argument nput')
*   ..test ns
      call xertst (1.le.ns, 'faulty argument ns')
*   ..test itim, itcnt
      call xertst (0.le.itim.and.0.le.itcnt,
     &  'faulty argument itim, itcnt')
*   ..set and test nfld
      nfld = 0
      if (1.le.switch%b2mwqt_style) then
       nfld = nfld+2+2*(2*ns+6)
      endif
      if (2.le.switch%b2mwqt_style) then
       nfld = nfld+36*(ns+3)
      endif
      if (3.le.switch%b2mwqt_style) then
       nfld = nfld+40*(ns+3)
      endif
      if (4.le.switch%b2mwqt_style) then
       nfld = nfld+72*(ns+3)
      endif
      call xertst (nfld.le.nfldmx, 'faulty parameter nfldmx')
*   ..define region boundaries
* ..get data for output
      k = 0
      if (1.le.switch%b2mwqt_style) then
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
       z(k) = diag%areshn
       k = k+1
       z(k) = diag%arespo
       k = k+1
       z(k) = diag%areskt
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
       z(k) = diag%acortn/ev
       k = k+1
       z(k) = diag%acorpo
       k = k+1
       z(k) = diag%acorkt/ev
       k = k+1
      endif
      if (2.le.switch%b2mwqt_style) then
*    ..total content on regions
       call xerrab ('WG - to be done!')
      endif
      if (3.le.switch%b2mwqt_style) then
*    ..integrated fluxes
       call xerrab ('WG - to be done!')
      endif
      if (4.le.switch%b2mwqt_style) then
*    ..min and max on regions
       call xerrab ('WG - to be done!')
      endif
*   ..test k
      call xertst (k.eq.nfld, 'faulty parameter nfld')

* ..produce output
      call cfwure (nput, nfld, z, 'data')

* ..return
      call subend ()
      return

*-----------------------------------------------------------------------
*.end b2mwqt

      end subroutine b2mwqt_diff
