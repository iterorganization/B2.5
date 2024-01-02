      function tfluxav(quant,mpg,geo,iFt)

*---------------------------------------------------------------------
*.documentation
*
*   1. purpose
*
*      TFLUXAV computes the flux tube average of a quantity
*
*   2. parameters
*
*      quant - array(1:mpg%nCv) quantity to be averaged
*      iFt   - label of the flux tube
*
*   3. description
*
*      This function computes and returns the flux tube average
*      of the quantity quant on the flux tube with the label iFt.
*      It allows for a special treatment if Ift is the first flux
*      tube outside the separatrix, otherwise it assumes the flux
*      tube is closed.
*----------------------------------------------------------------------
*.declarations
      use b2mod_types
      use b2us_map
      use b2us_geo
      implicit none
      type(mapping), intent(in) :: mpg
      type(geometry), intent(in) :: geo
      real (kind=R8) precision, tfluxav
      parameter (precision=1e-8_R8)
*   ..input arguments (unchanged on exit)
      real (kind=R8) ::
     *   quant(mpg%nCv)
*   ..local variables
      real (kind=R8) :: 
     *   qsum, wsum
      integer i, j, iFt, iCv, iCn, iFc
      logical on_separatrix
*   ..procedures
      external xertst

*.flux surface averaging
      call xertst (1.le.iFt.and.iFt.le.mpg%nFt, 'faulty argument iFt')
      qsum=0.0_R8
      wsum=0.0_R8
      do i=1,mpg%ftCvP(iFt,2)
       iCv=mpg%ftCv(mpg%ftCvP(iFt,1)+i-1)
       on_separatrix=.true.
       if(.not.mpg%cvOnClosedSurface(iCv)) then
         on_separatrix=.false.
         do j=1,mpg%cvFcP(iCv,2)
           if (on_separatrix) cycle
           iFc=mpg%cvFc(mpg%cvFcP(iCv,1)+j-1)
           if(mpg%fcCv(iFc,1).eq.iCv) then
             iCn=mpg%fcCv(iFc,2)
           else
             iCn=mpg%fcCv(iFc,1)
           endif
           on_separatrix=mpg%cvOnClosedSurface(iCn)
         end do
       endif
       if (.not.on_separatrix) cycle
       qsum=qsum+quant(iCv)*geo%cvHx(iCv)/geo%cvBb(iCv,0)
       wsum=wsum+geo%cvHx(iCv)/geo%cvBb(iCv,0)
      enddo
      if(abs(wsum).lt.precision) then  ! psauter
        tfluxav=0.0_R8
      else
        tfluxav=qsum/wsum
      endif

*   ..return
      return
      end function tfluxav
