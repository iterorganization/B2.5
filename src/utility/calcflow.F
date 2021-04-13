
*-----------------------------------------------------------------------
*.specification
! $AD NOCHECKPOINT
      subroutine calcflow(nCv,nFc,nVx,meth,geo,mpg,fun,flo,con,
     *  flow,flowf,flowc)
      use b2mod_types
      use b2us_geo
      use b2us_map
      implicit none
*   ..input arguments
      integer nCv, nFc, nVx, meth
      type (geometry), intent (in) :: geo
      type (mapping), intent (in) :: mpg
      real (kind=R8) ::
     *   fun(nCv), flo(nFc,0:1), con(nFc,0:1)
*   ..output arguments
      real (kind=R8) ::
     *   flow(nFc,0:1), flowf(nFc,0:1), flowc(nFc,0:1)
*-----------------------------------------------------------------------
*.documentation
*
*     CALCFLOW discretizes the convective-conductive flow of quantity
*     fun, on cell faces.
*
*     Input:
*     flo(nFc,0) contains the convective coefficient, poloidal direction
*     flo(nFc,1) contains the convective coefficient, radial direction
*     con(nFc,0) contains the conductive coefficient, poloidal direction
*     con(nFc,1) contains the conductive coefficient, radial direction
*     
*     Output:
*     flow(nFc,0) : the poloidal flux on cell faces
*     flow(nFc,1) : the radial flux on cell faces
* 
*     flowf: convective part of the flow
*     flowc: conductive part of the flow
*
*-----------------------------------------------------------------------
*.declarations

*   ..local variables
      integer iFc
      real (kind=R8) ::
     *  dfun(nFc,0:1), funv(nVx), fl, df, t0, t1
      real (kind=R8), parameter :: eps = 1.0e-60_R8

*   ..procedures
      real (kind=R8) ::
     *   hybr, upwind, cvmgp
      external cvmgp
      intrinsic abs, max, min, mod, sqrt, sign
      hybr(fl,df) = max(abs(df),abs(fl/2.),sqrt(df**2+fl**2/4.)) !!!
      upwind(fl,t0,t1) = max(fl,0.0_R8)*t0+min(fl,0.0_R8)*t1

*-----------------------------------------------------------------------
*.computation      

      if (meth .eq. 0) then ! CENTRAL scheme

*    ..differences of fun at cell faces
       call diff(nCv, nFc, nVx, 0, geo, mpg, fun, funv, dfun)

*    ..convective piece
       do iFc = 1,nFc
*    ..poloidal
        flowf(iFc,0) =
     &    0.5_R8*flo(iFc,0)*(fun(mpg%fcCv(iFc,1))+fun(mpg%fcCv(iFc,2)))
*    ..radial
        flowf(iFc,1) =
     &    0.5_R8*flo(iFc,1)*(fun(mpg%fcCv(iFc,1))+fun(mpg%fcCv(iFc,2)))
       enddo

*    ..add conductive piece
       flow = flowf - con*dfun

      elseif (meth .eq. 1) then ! UPWIND scheme

*    ..differences of fun at cell faces
       call diff(nCv, nFc, nVx, 0, geo, mpg, fun, funv, dfun)

*    ..convective piece
       do iFc = 1,nFc
*    ..poloidal
        flowf(iFc,0) =
     &    upwind(flo(iFc,0),fun(mpg%fcCv(iFc,1)),fun(mpg%fcCv(iFc,2))) 
*    ..radial
        flowf(iFc,1) = 
     &    upwind(flo(iFc,1),fun(mpg%fcCv(iFc,1)),fun(mpg%fcCv(iFc,2)))
       enddo

*    ..add conductive piece
       flow = flowf - con*dfun

      elseif (meth .eq. 2) then ! SOLPS5 hybrid scheme

*    ..differences of fun at cell faces
       call diff(nCv, nFc, nVx, 0, geo, mpg, fun, funv, dfun)

       do iFc = 1,nFc

*    ..poloidal flow
        t0 = (hybr(flo(iFc,0),con(iFc,0)*geo%fcQalf(iFc,0)) - 
     &    0.5_R8*abs(flo(iFc,0)))/(geo%fcQalf(iFc,0) + eps)!! 

        flowf(iFc,0) =
     &    upwind(flo(iFc,0),fun(mpg%fcCv(iFc,1)),fun(mpg%fcCv(iFc,2))) 
        flow(iFc,0) = 
     &    flowf(iFc,0) - t0*dfun(iFc,0)

*    ..radial flow
        t0 = (hybr(flo(iFc,1),con(iFc,1)*geo%fcQalf(iFc,1)) - 
     &    0.5_R8*abs(flo(iFc,1)))/(geo%fcQalf(iFc,1) + eps)  
        flowf(iFc,1) = 
     &    upwind(flo(iFc,1),fun(mpg%fcCv(iFc,1)),fun(mpg%fcCv(iFc,2)))
        flow(iFc,1) = 
     &    flowf(iFc,1) - t0*dfun(iFc,1)

       enddo

      elseif (meth .eq. 3) then ! SOLPS4 continuity scheme

        call xerrab('calcflow -- meth 3 to be implemented')      

      else
        call xerrab('calcflow -- wrong value of meth')
      endif

*   ..conductive part of the flow
      flowc = flow - flowf

*   ..return
      return

*-----------------------------------------------------------------------
*.end calcflow

      end subroutine calcflow



