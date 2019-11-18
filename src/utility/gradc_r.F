*-----------------------------------------------------------------------
*.specification

      subroutine gradc_r(nCv, nFc, nVx, mode, geo, mpg, fun, funv,
     &  gfunr)
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
     *   gfunr(nCv)
*-----------------------------------------------------------------------
*.documentation
*
*     GRADC_R computes gradients of a cell centered quantity fun in the 
*     radial direction, in cell centers. The routine first computes
*     the gradients on faces, and then interpolates to cell centers.
*
*-----------------------------------------------------------------------
*.declarations

*   ..local variables
      real(kind=R8) :: gfunrf(nFc)

*-----------------------------------------------------------------------
*.computation      

*   ..test mode
      call xertst (0.le.mode.and.mode.le.1,
     .  'gradc_r--faulty argument mode')

*   ..interpolate fun to cell vertices, using volume weighted averaging
      if (mode.eq.0) then
        call intvertex(nCv,nVx,mpg,geo%vxVol,fun,funv)
      endif

*   ..compute gradients on faces
      call grad_r(nCv, nFc, nVx, 1, geo, mpg, fun, funv, gfunrf)

*   ..interpolate to centers
      call intcell(nFc, nCv, mpg, mpg%intcellR, gfunrf,gfunr)

*   ..return
      return

*-----------------------------------------------------------------------
*.end gradc_r

      end subroutine gradc_r
