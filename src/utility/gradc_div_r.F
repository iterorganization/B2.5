*-----------------------------------------------------------------------
*.specification

      subroutine gradc_div_r(nCv, nFc, nVx, mode, geo, mpg, fun, funf,
     &  gfunr)
      use b2mod_types
      use b2us_geo
      use b2us_map
      implicit none
*   ..input arguments
      integer, intent (in) :: nCv, nFc, nVx, mode
      type (geometry), intent (in) :: geo
      type (mapping), intent (in) :: mpg
      real (kind=R8) ::
     *   fun(nCv), funf(nFc)
*   ..output arguments
      real (kind=R8) ::
     *   gfunr(nCv)
*-----------------------------------------------------------------------
*.documentation
*
*     GRADC_DIV_R computes gradients of a cell centered quantity fun in 
*     the radial direction, in cell centers. Computation is based on
*     the divergence theorem.
*     For guard cells, the value of the gradient is set to 0.
*
*-----------------------------------------------------------------------
*.declarations

*   ..local variables
      real(kind=R8) :: wrkc(nCv), wrkf(nFc)

*-----------------------------------------------------------------------
*.computation      

*   ..test mode
      call xertst (0.le.mode.and.mode.le.1,
     .  'gradc_div_r--faulty argument mode')

*   ..interpolate fun to cell faces, using volume weighted averaging
      if (mode.eq.0) then
        call intface(nCv,nFc,mpg%fcCv,geo%fcVol,fun,funf)
      endif

*   ..compute gradient
      wrkf = funf*geo%fcS*geo%fcQalf(:,1)
      call div0(nCv, nFc, mpg, wrkf, gfunr)
      wrkf = geo%fcS*geo%fcQalf(:,0)
      call div0(nCv, nFc, mpg, wrkf,wrkc)
      gfunr = (gfunr - fun*wrkc)/geo%cvVol

*   ..set zero in guard cells (no reliable information from div0 there!)
      gfunr(mpg%nCi+1:mpg%nCv) = 0.0_R8

*   ..return
      return

*-----------------------------------------------------------------------
*.end gradc_div_r

      end subroutine gradc_div_r
