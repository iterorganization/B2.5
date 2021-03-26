*-----------------------------------------------------------------------
*.specification

      subroutine diff(nCv, nFc, nVx, mode, geo, mpg, fun, funv, dfun)
      use b2mod_types
      use b2us_geo
      use b2us_map
      implicit none
*   ..input arguments
      integer nCv, nFc, nVx, mode
      type (geometry), intent (in) :: geo
      type (mapping), intent (in) :: mpg
      real (kind=R8) ::
     *   fun(nCv), funv(nVx)
*   ..output arguments
      real (kind=R8) ::
     *   dfun(nFc,0:1)
*-----------------------------------------------------------------------
*.documentation
*
*     DIFF computes differences of a cell centered quantity fun on
*     faces of a cell.
*     dfun(:,0) will hold poloidal differences
*     dfun(:,1) will hold radial differences
*
*     The routine is a simple interface to DIFF_P and DIFF_R, that will
*     compute the actual poloidal and radial differences
*
*-----------------------------------------------------------------------
*.declarations

*   ..procedures
      external diff_p, diff_r

*-----------------------------------------------------------------------
*.computation      

*   ..test mode
      call xertst (0.le.mode.and.mode.le.1,
     .  'diff--faulty argument mode')

*   ..interpolate fun to cell vertices, using volume weighted averaging
      if (mode.eq.0) then
        call intvertex(nCv,nVx,mpg,geo%vxVol,fun,funv)
      endif

*   ..compute poloidal and radial differences
      call diff_p(nCv, nFc, nVx, 1, geo, mpg, fun, funv, dfun(:,0))
      call diff_r(nCv, nFc, nVx, 1, geo, mpg, fun, funv, dfun(:,1))

*   ..return
      return

*-----------------------------------------------------------------------
*.end diff

      end subroutine diff
