*-----------------------------------------------------------------------
*.specification

      subroutine diff_p(nCv, nFc, nVx, mode, geo, mpg, fun, funv, dfunp)
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
     *   dfunp(nFc)
*-----------------------------------------------------------------------
*.documentation
*
*     DIFF_P computes differences of a cell centered quantity fun in the 
*     poloidal direction, on faces of a cell. 
*
*     To compute values of the field in cell vertices, the routine relies 
*     calls to intvertex.
*
*     If mode.eq.0, the routine will interpolate fun to vertices (funv)
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
     .  'diff_p--faulty argument mode')

*   ..interpolate fun to cell vertices, using volume weighted averaging
      if (mode.eq.0) then
        call intvertex(nCv,nVx,mpg,geo%vxVol,fun,funv)
      endif

*   ..compute poloidal difference at faces
      do iFc = 1,nFc
       dfunp(iFc) = 
     &  (fun(mpg%fcCv(iFc,2))  - fun(mpg%fcCv(iFc,1)))*
     &   geo%fcQalf(iFc,0) +
     &  (funv(mpg%fcVx(iFc,2)) - funv(mpg%fcVx(iFc,1)))*
     &  (geo%fcHc(iFc,1) + geo%fcHc(iFc,2))/geo%fcHt(iFc)*
     &   geo%fcQbet(iFc,1)
      enddo

*   ..return
      return

*-----------------------------------------------------------------------
*.end diff_p

      end subroutine diff_p
