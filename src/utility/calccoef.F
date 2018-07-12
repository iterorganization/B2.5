*-----------------------------------------------------------------------
*.specification

      subroutine calccoef(nCv,nFc,nVx,meth,geo,mpg,flo,con,flo0,con0)
      use b2mod_types
      use b2us_geo
      use b2us_map
      implicit none
*   ..input arguments
      integer nCv, nFc, nVx, meth
      type (geometry), intent (in) :: geo
      type (mapping), intent (in) :: mpg
      real (kind=R8) ::
     *   flo(nFc,0:1), con(nFc,0:1)
*   ..output arguments
      real (kind=R8) ::
     *   flo0(nFc,0:1), con0(nFc,0:1)
*-----------------------------------------------------------------------
*.documentation
*
*     CALCCOEF computes the coefficients flo0 and con0 for the 
*     linearization of the total (poloidal + radial) convective-conductive
*     flow across a face in the standardized form
*     
*     flow = 0.5*flo0(,0)*(phi_c1 + phi_c2) - con0(,0)*(phi_c2 - phi_c1)
*          + 0.5*flo0(,1)*(phi_v1 + phi_v2) - con0(,1)*(phi_v2 - phi_v1)
*
*     Input:
*     flo(nFc,0) contains the convective coefficient, poloidal direction
*     flo(nFc,1) contains the convective coefficient, radial direction
*     con(nFc,0) contains the conductive coefficient, poloidal direction
*     con(nFc,1) contains the conductive coefficient, radial direction
*     
*
*-----------------------------------------------------------------------
*.declarations

*   ..local variables
      integer iFc
      real (kind=R8) ::
     *  fl, df, t0, t1, t00, t01, t10, t11
      real (kind=R8), save :: eps = 1.0e-60_R8

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

       do iFc = 1,nFc
*    ..geometric factors for cell and vertex contributions
*      to differences of cell faces (consistent with diff_p, diff_r)
        t00 = geo%fcQalf(iFc,0)
        t01 = (geo%fcHc(iFc,1) + geo%fcHc(iFc,2))/geo%fcHt(iFc)*
     &         geo%fcQbet(iFc,1)
        t10 = geo%fcQalf(iFc,1) 
        t11 = (geo%fcHc(iFc,1) + geo%fcHc(iFc,2))/geo%fcHt(iFc)*
     &         geo%fcQbet(iFc,0)
*    ..convective contributions to flob0, conb0
        flo0(iFc,0) = flo(iFc,0) + flo(iFc,1)
        flo0(iFc,1) = 0.0_R8
!        flo0(iFc,0) = 0.5_R8*(flo(iFc,0) + flo(iFc,1))
!        flo0(iFc,1) = 0.5_R8*(flo(iFc,0) + flo(iFc,1))
*    ..conductive contributions to flob0, conb0
        con0(iFc,0) = con(iFc,0)*t00 + con(iFc,1)*t10
        con0(iFc,1) = con(iFc,0)*t01 + con(iFc,1)*t11
       enddo

      elseif (meth .eq. 1) then ! UPWIND scheme

        call xerrab('calccoef -- meth 1 to be implemented')      

      elseif (meth .eq. 2) then ! SOLPS5 hybrid scheme

       do iFc = 1,nFc
*    ..effective diffusion coefficients
        t0 = sign(hybr(flo(iFc,0),con(iFc,0)*geo%fcQalf(iFc,0)) - 
     &    0.5_R8*abs(flo(iFc,0)), con(iFc,0))/(geo%fcQalf(iFc,0) + eps)  
        t1 = sign(hybr(flo(iFc,1),con(iFc,1)*geo%fcQalf(iFc,1)) - 
     &    0.5_R8*abs(flo(iFc,1)), con(iFc,1))/(geo%fcQalf(iFc,1) + eps)  
*    ..compute geometric factors for cell and vertex contributions
*      to differences of cell faces (consistent with diff_p, diff_r)
        t00 = geo%fcQalf(iFc,0)
        t01 = (geo%fcHc(iFc,1) + geo%fcHc(iFc,2))/geo%fcHt(iFc)*
     &         geo%fcQbet(iFc,1)
        t10 = geo%fcQalf(iFc,1) 
        t11 = (geo%fcHc(iFc,1) + geo%fcHc(iFc,2))/geo%fcHt(iFc)*
     &         geo%fcQbet(iFc,0)
*    ..convective contributions to flob0, conb0
        flo0(iFc,0) = flo(iFc,0) + flo(iFc,1)
        flo0(iFc,1) = 0.0_R8
!        flo0(iFc,0) = 0.5_R8*(flo(iFc,0) + flo(iFc,1))
!        flo0(iFc,1) = 0.5_R8*(flo(iFc,0) + flo(iFc,1))
*    ..conductive contributions to flob0, conb0
        con0(iFc,0) = (t0 + 0.5_R8*abs(flo(iFc,0)))*t00 + 
     &                (t1 + 0.5_R8*abs(flo(iFc,1)))*t10
        con0(iFc,1) = t0*t01 + t1*t11
       enddo

      elseif (meth .eq. 3) then ! SOLPS4 continuity scheme

        call xerrab('calccoef -- meth 3 to be implemented')      

      else
        call xerrab('calccoef -- wrong value of meth')
      endif


*   ..return
      return

*-----------------------------------------------------------------------
*.end calccoef

      end subroutine calccoef

