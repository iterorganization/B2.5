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
*      The integration is made using the trapeze rule.
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
     *   qp, qm, qsum, wsum
      integer i, iFt, iCv, iCp, iCm
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
          quantn(iCm)=sqrt(1.0_R8-nlambda*geo%cvBb(iCm,3)/bmax)
          quantn(iCv)=sqrt(1.0_R8-nlambda*geo%cvBb(iCv,3)/bmax)
          quantn(iCp)=sqrt(1.0_R8-nlambda*geo%cvBb(iCp,3)/bmax)
          qm=(quantn(iCm)*geo%cvHx(iCv)/geo%cvBb(iCv,0)+
     .        quantn(iCv)*geo%cvHx(iCm)/geo%cvBb(iCm,0))/
     .       (geo%cvHx(iCm)/geo%cvBb(iCm,0)+
     .        geo%cvHx(iCv)/geo%cvBb(iCv,0))
          qp=(quantn(iCp)*geo%cvHx(iCv)/geo%cvBb(iCv,0)+
     .        quantn(iCv)*geo%cvHx(iCp)/geo%cvBb(iCp,0))/
     .       (geo%cvHx(iCv)/geo%cvBb(iCv,0)+
     .        geo%cvHx(iCp)/geo%cvBb(iCp,0))
          qsum=qsum+0.5_R8*geo%cvHx(iCv)/geo%cvBb(iCv,0)*
     .       (0.5_R8*(qm+quant(iCv))+0.5_R8*(quant(iCv)+qp))
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
