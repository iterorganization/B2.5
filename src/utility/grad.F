*-----------------------------------------------------------------------
*.specification

      subroutine grad(nCv, nFc, nVx, mode, geo, mpg, fun, funv, gfun)
      use b2mod_types
      use b2us_geo
      use b2us_map
      implicit none
*   ..input arguments
      integer, intent (in) :: nCv, nFc, nVx, mode
      type (geometry), intent (in) :: geo
      type (mapping), intent (in) :: mpg
      real (kind=R8) ::
     *   fun(nCv), funv(nVx)
*   ..output arguments
      real (kind=R8) ::
     *   gfun(nFc,0:1)
*-----------------------------------------------------------------------
*.documentation
*
*     GRAD computes gradients of a cell centered quantity fun on
*     faces of a cell.
*     gfun(:,0) will hold poloidal gradients
*     gfun(:,1) will hold radial gradients
*
*     The routine is a simple interface to GRAD_P and GRAD_R, that will
*     compute the actual poloidal and radial gradients
*
*-----------------------------------------------------------------------
*.declarations

*   ..procedures
      external grad_p, grad_r

*-----------------------------------------------------------------------
*.computation      

*   ..test mode
      call xertst (0.le.mode.and.mode.le.1,
     .  'grad--faulty argument mode')

*   ..interpolate fun to cell vertices, using volume weighted averaging
      if (mode.eq.0) then
        call intvertex(nCv,nVx,mpg,geo%vxVol,fun,funv)
      endif

*   ..compute poloidal and radial gradients
      call grad_p(nCv, nFc, nVx, 1, geo, mpg, fun, funv, gfun(1,0))
      call grad_r(nCv, nFc, nVx, 1, geo, mpg, fun, funv, gfun(1,1))

*   ..return
      return

*-----------------------------------------------------------------------
*.end grad

      end subroutine grad
