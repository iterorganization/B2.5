*-----------------------------------------------------------------------
*.specification
! $AD NOCHECKPOINT
      subroutine grad_p(nCv, nFc, nVx, mode, geo, mpg, fun, funv, gfunp)
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
     *   gfunp(nFc)
      integer cvs(2)
      real (kind=R8) :: radcoef
*-----------------------------------------------------------------------
*.documentation
*
*     GRAD_P computes gradients of a cell centered quantity fun in the 
*     poloidal direction, on faces of a cell.
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
     .  'grad_p--faulty argument mode')

*   ..interpolate fun to cell vertices, using volume weighted averaging
      if (mode.eq.0) then
        call intvertex(nCv,nVx,mpg,geo%vxVol,fun,funv)
      endif

*   ..compute poloidal gradients at faces
      do iFc = 1,nFc
       cvs = mpg%fcCv(iFc,:)
       if (mpg%cvFt(cvs(1)).eq.mpg%cvFt(cvs(2))) then
         radcoef = 0.0_R8
       else
         radcoef = 1.0_R8
       endif
       gfunp(iFc) = 
     &  (fun(mpg%fcCv(iFc,2))  - fun(mpg%fcCv(iFc,1)))*
     &   geo%fcQalf(iFc,0)/(geo%fcQgam(iFc,0)*
     &  (geo%fcHc(iFc,1) + geo%fcHc(iFc,2))) +
     &  (funv(mpg%fcVx(iFc,2)) - funv(mpg%fcVx(iFc,1)))*
     &   geo%fcQbet(iFc,1)/(geo%fcQgam(iFc,0)*geo%fcHt(iFc))*radcoef
      enddo

*   ..return
      return

*-----------------------------------------------------------------------
*.end grad_p

      end subroutine grad_p
