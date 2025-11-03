c
      subroutine psy()
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
     .     geo%cvHx(iCv)/geo%cvHz(iCv)
          values(9)=min(values(9),work(iCv,is,iwrk))
          values(10)=max(values(10),work(iCv,is,iwrk))
        enddo
      enddo
      stmp=trim(graphlabel(iwrk))//'/(m^{2})'
      call customise_label(stmp,graphlabel(iwrk))
      return
      end subroutine psy
