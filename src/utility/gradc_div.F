*-----------------------------------------------------------------------
*.specification

      subroutine gradc_div(nCv, nFc, nVx, mode, geo, mpg, fun, funf,
     .  gfun)
      use b2mod_types
      use b2us_geo
      use b2us_map
      implicit none
*   ..input arguments
      integer nCv, nFc, nVx, mode
      type (geometry), intent (in) :: geo
      type (mapping), intent (in) :: mpg
      real (kind=R8) ::
     *   fun(nCv), funf(nFc)
*   ..output arguments
      real (kind=R8) ::
     *   gfun(nCv,0:1)
*-----------------------------------------------------------------------
*.documentation
*
*     GRADC_DIV computes gradients of a cell centered quantity fun in
*     cell centers.  Computation is based on the divergence theorem.
*     dfun(:,:,0) will hold the gradient in the poloidal direction
*     dfun(:,:,1) will hold the gradient in the radial direction
*
*     The routine is a simple interface to GRADC_DIV_P and GRADC_DIV_R, 
*     that compute the actual poloidal and radial gradients.
*
*-----------------------------------------------------------------------
*.declarations

*   ..local variables

*   ..procedures
      external gradc_div_p, gradc_div_r

*-----------------------------------------------------------------------
*.computation      

*   ..test mode
      call xertst (0.le.mode.and.mode.le.1,
     .  'gradc_div--faulty argument mode')

*   ..interpolate fun to cell faces, using volume weighted averaging
      if (mode.eq.0) then
        call intface(nCv,nFc,mpg%fcCv,geo%fcVol,fun,funf)
      endif

*   ..compute poloidal and radial gradients
      call gradc_div_p(nCv, nFc, nVx, 1, geo, mpg, fun, funf, gfun(:,0))
      call gradc_div_r(nCv, nFc, nVx, 1, geo, mpg, fun, funf, gfun(:,1))

*   ..return
      return

*-----------------------------------------------------------------------
*.end gradc_div

      end subroutine gradc_div
