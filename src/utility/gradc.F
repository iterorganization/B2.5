*-----------------------------------------------------------------------
*.specification

      subroutine gradc(nCv, nFc, nVx, mode, geo, mpg, fun, funv, gfun)
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
     *   gfun(nCv,0:1)
*-----------------------------------------------------------------------
*.documentation
*
*     GRADC computes gradients of a cell centered quantity fun in cell
*     centers.
*     dfun(:,:,0) will hold the gradient in the poloidal direction
*     dfun(:,:,1) will hold the gradient in the radial direction
*
*     The routine is a simple interface to GRADC_P and GRADC_R, that 
*     compute the actual poloidal and radial gradients
*
*-----------------------------------------------------------------------
*.declarations

*   ..local variables

*   ..procedures
      external gradc_p, gradc_r

*-----------------------------------------------------------------------
*.computation      

*   ..test mode
      call xertst (0.le.mode.and.mode.le.1,
     .  'gradc--faulty argument mode')

*   ..interpolate fun to cell vertices, using volume weighted averaging
      if (mode.eq.0) then
        call intvertex(nCv,nVx,mpg,geo%vxVol,fun,funv)
      endif

*   ..compute poloidal and radial gradients
      call gradc_p(nCv, nFc, nVx, 1, geo, mpg, fun, funv, gfun(:,0))
      call gradc_r(nCv, nFc, nVx, 1, geo, mpg, fun, funv, gfun(:,1))

*   ..return
      return

*-----------------------------------------------------------------------
*.end gradc

      end subroutine gradc
