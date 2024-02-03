      subroutine b2txsx (nCv, nFc, geo, mpg, fcVol, fcS, fun, funsx)
      use b2mod_types
      use b2us_geo
      use b2us_map
      use b2mod_subsys
      implicit none
      integer nCv, nFc
      type (geometry), intent (in) :: geo
      type (mapping), intent (in) :: mpg
      real (kind=R8) ::
     &  fcVol(nFc,2), fcS(nFc), fun(nCv), funsx(nFc)
*     ------------------------------------------------------------------
*     B2TXSX computes the integrals over cell x-faces of a function
*     that is provided on the cells.
*!!   (Reconsider specification and use of this routine for the case
*     of a non-orthogonal grid.)
*     ------------------------------------------------------------------
      integer iFc
*     ------------------------------------------------------------------
      call subini ('b2txsx')
*   ..test nCv, nFc
      call xertst (0.le.nCv .and. 0.le.nFc, 'faulty argument nCv,nFc')
*   ..compute funsx
      do iFc = 1, nFc
        funsx(iFc) = fcS(iFc)*geo%fcQalf(iFc,0)*
     &   (fcVol(iFc,2)*fun(mpg%fcCv(iFc,1)) + 
     &    fcVol(iFc,1)*fun(mpg%fcCv(iFc,2)))/
     &   (fcVol(iFc,1) + fcVol(iFc,2))
      enddo
*   ..return
      call subend ()
      return
*     ------------------------------------------------------------------
      end subroutine b2txsx
