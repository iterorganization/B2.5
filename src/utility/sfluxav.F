      function sfluxav(quant,mpg,geo,iFt,nlambda,bmax)

*---------------------------------------------------------------------
*.documentation
*
*   1. purpose
*
*      SFLUXAV computes a specific flux tube average of a quantity
*
*   2. parameters
*
*      quant - array(1:mpg%nCv) quantity to be averaged
*      iFt   - label of the flux tube
*
*   3. description
*
*      This function computes and returns a specific flux tube average
*      of the quantity quant on the flux tube with label iFt.
*      It assumes that the flux tube is closed, i.e., inside the core
*      or the the pedestal region.
*----------------------------------------------------------------------
*.declarations
      use b2mod_types
      use b2us_map
      use b2us_geo
      implicit none
      type(mapping), intent(in) :: mpg
      type(geometry), intent(in) :: geo
      real (kind=R8) precision, sfluxav
      parameter (precision=1e-8_R8)
*   ..input arguments (unchanged on exit)
      real (kind=R8) ::
     *   quant(mpg%nCv), nlambda, bmax, quantn(mpg%nCv)
*   ..local variables
      real (kind=R8) :: 
     *   qsum, wsum
      integer i, iFt, iCv
*   ..procedures
      external xertst

*.flux surface averaging
      call xertst (1.le.iFt.and.iFt.le.mpg%nFt, 'faulty argument iFt')
      qsum=0.0_R8
      wsum=0.0_R8
*   ..working variable quantn
      quantn=quant
      do i=1,mpg%ftCvP(iFt,2)
       iCv=mpg%ftCv(mpg%ftCvP(iFt,1)+i-1)
       if (mpg%cvOnclosedsurface(iCv)) then
        quantn(iCv)=sqrt(1._R8-nlambda*geo%cvBb(iCv,3)/bmax)
#ifdef WG_TODO
        if (1.0_R8-nlambda*bb(topix(ix,iy),topiy(ix,iy),3)/bmax
     &   .lt.0.0_R8)  then
          quantn(topix(ix,iy),topiy(ix,iy))=0.0_R8
         else
          quantn(topix(ix,iy),topiy(ix,iy))=real(sqrt(1.0_R8-nlambda*
     &           bb(topix(ix,iy),topiy(ix,iy),3)/bmax),R8)
        endif
        if (1.0_R8-nlambda*bb(bottomix(ix,iy),bottomiy(ix,iy),3)/bmax
     &   .lt.0.0_R8) then
          quantn(bottomix(ix,iy),bottomiy(ix,iy))=0.0_R8
         else
          quantn(bottomix(ix,iy),bottomiy(ix,iy))=real(
     &          sqrt(1.0_R8-nlambda*bb(bottomix(ix,iy),
     &                                 bottomiy(ix,iy),3)/bmax),R8)
        endif
#endif
        qsum=qsum+quantn(iCv)*geo%cvHx(iCv)/geo%cvBb(iCv,0)
        wsum=wsum+geo%cvHx(iCv)/geo%cvBb(iCv,0)
       endif
      enddo
      if(abs(wsum).lt.precision) then  ! psauter
        sfluxav=0.0_R8
      else
        sfluxav=qsum/wsum
      endif

*   ..return
      return
      end function sfluxav
