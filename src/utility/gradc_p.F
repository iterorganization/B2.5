*-----------------------------------------------------------------------
*.specification

      subroutine gradc_p(nCv, nFc, nVx, mode, geo, mpg, fun, funv,
     &  gfunp)
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
     *   gfunp(nCv)
*-----------------------------------------------------------------------
*.documentation
*
*     GRADC_P computes gradients of a cell centered quantity fun in the 
*     poloidal direction, in cell centers. The routine first computes
*     the gradients on faces, and then interpolates to cell centers.
*
*-----------------------------------------------------------------------
*.declarations

*   ..local variables
      real(kind=R8) :: gfunpf(nFc)

*-----------------------------------------------------------------------
*.computation      

*   ..test mode
      call xertst (0.le.mode.and.mode.le.1,
     .  'gradc_p--faulty argument mode')

*   ..interpolate fun to cell vertices, using volume weighted averaging
      if (mode.eq.0) then
        call intvertex(nCv,nVx,mpg,geo%vxVol,fun,funv)
      endif

*   ..compute gradients on faces
      call grad_p(nCv, nFc, nVx, 1, geo, mpg, fun, funv, gfunpf)

*   ..interpolate to centers
      call intcell(nFc, nCv, mpg, mpg%intcellP, gfunpf,gfunp)

*   ..return
      return

*-----------------------------------------------------------------------
*.end gradc_p

      end subroutine gradc_p
