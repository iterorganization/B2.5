c
      subroutine pvol()
      use b2mod_types
      use b2us_data
      use b2mod_b2plot
      use b2mod_b2plot_values
      implicit none
      integer is,lnblnk
      character*256 stmp
      external lnblnk

      do is=0,ndim(2,iwrk)-1
        work(:,is,iwrk)=work(:,is,iwrk)/geo%cvVol(:)
      enddo
      values(9)=minval(work(:,0:ndim(2,iwrk)-1,iwrk))
      values(10)=maxval(work(:,0:ndim(2,iwrk)-1,iwrk))
      stmp=graphlabel(iwrk)(1:lnblnk(graphlabel(iwrk)))//'/(m^{3})'
      call customise_label(stmp,graphlabel(iwrk))
      return
      end subroutine pvol
