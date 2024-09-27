*-----------------------------------------------------------------------
*.specification

      subroutine grad_r(nCv, nFc, nVx, mode, geo, mpg, fun, funv, gfunr)
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
     *   gfunr(nFc)
*-----------------------------------------------------------------------
*.documentation
*
*     GRAD_R computes gradients of a cell centered quantity fun in the 
*     radial direction, on faces of a cell.
*
*     To compute values of the field in cell vertices, the routine relies 
*     calls to intvertex.
*
*     If mode.eq.0, the routine will interpolate fun to vertices in funv
*     If mode.eq.1, the routine assumess funv has been pre-computed
*
*-----------------------------------------------------------------------
*.declarations

*   ..local variables
      integer iFc

*   ..procedures
      external intvertex

*-----------------------------------------------------------------------
*.computation      

*   ..test mode
      call xertst (0.le.mode.and.mode.le.1,
     .  'grad_r--faulty argument mode')

*   ..interpolate fun to cell vertices, using volume weighted averaging
      if (mode.eq.0) then
        call intvertex(nCv,nVx,mpg,geo%vxVol,fun,funv)
      endif

*   ..compute poloidal gradients at faces
      do iFc = 1, nFc
       gfunr(iFc) = 
     &  (fun(mpg%fcCv(iFc,2))  - fun(mpg%fcCv(iFc,1)))*
     &   geo%fcQalf(iFc,1)/(geo%fcQgam(iFc,0)*
     &  (geo%fcHc(iFc,1) + geo%fcHc(iFc,2))) +
     &  (funv(mpg%fcVx(iFc,2)) - funv(mpg%fcVx(iFc,1)))*
     &   geo%fcQbet(iFc,0)/(geo%fcQgam(iFc,0)*geo%fcHt(iFc))
      enddo

*   ..return
      return

*-----------------------------------------------------------------------
*.end grad_r

      end subroutine grad_r
