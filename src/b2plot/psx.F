c
      subroutine psx()
      use b2mod_types
      use b2us_data
      use b2mod_b2plot
      use b2mod_b2plot_values
      implicit none
      integer iCv,iFc,iFt1,iFt2,i,is
      real(kind=R8) :: sx
      intrinsic min,max
      character*256 stmp

      values(9) = huge(1.0_R8)
      values(10)=-huge(1.0_R8)
      do is=0,ndim(2,iwrk)-1
        do iCv=1,mpg%nCv
          sx=0.0_R8
          do i=1,mpg%cvFcP(iCv,2)
            iFc=mpg%cvFc(mpg%cvFcP(iCv,1)+i-1)
            iFt1=mpg%cvFt(mpg%fcCv(iFc,1))
            iFt2=mpg%cvFt(mpg%fcCv(iFc,2))
            if (iFt1.eq.iFt2.and.iFt1.ne.0.and.iFt2.ne.0)
     >       sx=sx+geo%fcHt(iFc)
          end do
          if (iCv.le.mpg%nCi) sx=0.5_R8*sx
          work(iCv,is,iwrk)=work(iCv,is,iwrk)/sx/geo%cvHz(iCv)
          values(9)=min(values(9),work(iCv,is,iwrk))
          values(10)=max(values(10),work(iCv,is,iwrk))
        enddo
      enddo
      stmp=trim(graphlabel(iwrk))//'/(m^{2})'
      call customise_label(stmp,graphlabel(iwrk))
      return
      end subroutine psx
c
      subroutine psxperp()
      use b2mod_types
      use b2us_data
      use b2mod_b2plot
      use b2mod_b2plot_values
      implicit none
      integer iCv,is
      intrinsic min,max
      character*256 stmp

      values(9) = huge(1.0_R8)
      values(10)=-huge(1.0_R8)
      do is=0,ndim(2,iwrk)-1
        do iCv=1,mpg%nCv
          work(iCv,is,iwrk)=work(iCv,is,iwrk)/
     .     geo%cvHy(iCv)/geo%cvHz(iCv)
          values(9)=min(values(9),work(iCv,is,iwrk))
          values(10)=max(values(10),work(iCv,is,iwrk))
        enddo
      enddo
      stmp=trim(graphlabel(iwrk))//'/(m^{2})'
      call customise_label(stmp,graphlabel(iwrk))
      return
      end subroutine psxperp
