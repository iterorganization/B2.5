      subroutine b2txsy (nCv, nFc, geo, mpg, fcVol, fcS, fun, funsy)
      use b2mod_types
      use b2us_geo
      use b2us_map
      implicit none
      integer nCv, nFc
      type (geometry), intent (in) :: geo
      type (mapping), intent (in) :: mpg
      real (kind=R8) ::
     &  fcVol(nFc,2), fcS(nFc), fun(nCv), funsy(nFc)
*     ------------------------------------------------------------------
*     B2TXSY computes the integrals over cell y-faces of a function
*     that is provided on the cells.
*!!   (Re-consider specification and use of this routine for the case
*     of a non-orthogonal grid.)
*     ------------------------------------------------------------------
      integer iFc
      external subini, subend, b2xxgs
*     ------------------------------------------------------------------
      call subini ('b2txsy')
*   ..compute funsy
      do iFc = 1, nFc
        funsy(iFc) = fcS(iFc)*geo%fcQalf(iFc,1)*
     &   (fcVol(iFc,2)*fun(mpg%fcCv(iFc,1)) + 
     &    fcVol(iFc,1)*fun(mpg%fcCv(iFc,2)))/
     &   (fcVol(iFc,1) + fcVol(iFc,2))
      enddo
*   ..reset null regions
!      call b2xxgs (nx, ny, 0.0_R8, funsy, 2)
*   ..return
      call subend ()
      return
*     ------------------------------------------------------------------
      end
