      function fluxav(quant,mpg,geo,iFt)

*---------------------------------------------------------------------
*.documentation
*
*   1. purpose
*
*      FLUXAV computes the flux tube average of a quantity
*
*   2. parameters
*
*      quant - array(1:mpg%nCv) quantity to be averaged
*      iFt   - label of the flux tube
*
*   3. description
*
*      This function computes and returns the flux tube average
*      of the quantity quant on the flux tube with index iFt.
*      It assumes that the flux tube is closed, i.e., inside the core
*      or the the pedestal region.
*      The integration is made using the trapeze rule.
*----------------------------------------------------------------------
*.declarations
      use b2mod_types
      use b2us_map
      use b2us_geo
      implicit none
      type(mapping), intent(in) :: mpg
      type(geometry), intent(in) :: geo

      real (kind=R8) :: precision, fluxav
      parameter (precision=1e-8_R8)
*   ..input arguments (unchanged on exit)
      integer iFt
      real (kind=R8) ::
     *   quant(mpg%nCv)
*   ..local variables
      integer i, iCv, iCm, iCp
      real (kind=R8) :: 
     *   qm, qp, qsum, wsum
*   ..procedures
      external xertst

*.flux tube averaging
      call xertst (1.le.iFt.and.iFt.le.mpg%nFt, 'faulty argument iFt')
      qsum=0.0_R8
      wsum=0.0_R8
      do i=1,mpg%ftCvP(iFt,2)
        iCv=mpg%ftCv(mpg%ftCvP(iFt,1)+i-1)
        if (mpg%cvOnClosedSurface(iCv)) then
          if (i.gt.1) then
            iCm=mpg%ftCv(mpg%ftCvP(iFt,1)+i-2)
          else
            iCm=mpg%ftCv(mpg%ftCvP(iFt,1)+mpg%ftCvP(iFt,2)-1)
          end if
          if (i.lt.mpg%ftCvP(iFt,2)) then
            iCp=mpg%ftCv(mpg%ftCvP(iFt,1)+i)
          else
            iCp=mpg%ftCv(mpg%ftCvP(iFt,1))
          end if
          qm=(quant(iCm)*geo%cvHx(iCv)/geo%cvBb(iCv,0)+
     .        quant(iCv)*geo%cvHx(iCm)/geo%cvBb(iCm,0))/
     .       (geo%cvHx(iCm)/geo%cvBb(iCm,0)+
     .        geo%cvHx(iCv)/geo%cvBb(iCv,0))
          qp=(quant(iCp)*geo%cvHx(iCv)/geo%cvBb(iCv,0)+
     .        quant(iCv)*geo%cvHx(iCp)/geo%cvBb(iCp,0))/
     .       (geo%cvHx(iCv)/geo%cvBb(iCv,0)+
     .        geo%cvHx(iCp)/geo%cvBb(iCp,0))
          qsum=qsum+0.5_R8*geo%cvHx(iCv)/geo%cvBb(iCv,0)*
     .       (0.5_R8*(qm+quant(iCv))+0.5_R8*(quant(iCv)+qp))
          wsum=wsum+geo%cvHx(iCv)/geo%cvBb(iCv,0)
        endif
      enddo
      if(abs(wsum).lt.precision) then  ! psauter
        fluxav=0.0_R8
      else
        fluxav=qsum/wsum
      endif

*   ..return
      return
      end function fluxav
