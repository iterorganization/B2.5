c
      subroutine pswa()
      use b2mod_types
      use b2mod_wall
      use b2mod_b2plot
      use b2mod_b2plot_values
      implicit none
      integer iw,is,lnblnk
      character*256 stmp
      intrinsic min,max
      external lnblnk

      values(14)=wall_stack(1,0,iwwrk)/wall_area(1)
      values(15)=wall_stack(1,0,iwwrk)/wall_area(1)
      do is=0,nwdim(2,iwwrk)-1
        do iw=1,nbnd
          wall_stack(iw,is,iwwrk)=wall_stack(iw,is,iwwrk)/wall_area(iw)
          values(14)=min(values(14),wall_stack(iw,is,iwwrk))
          values(15)=max(values(15),wall_stack(iw,is,iwwrk))
        enddo
      enddo
      stmp=wall_label(iwwrk)(1:lnblnk(wall_label(iwwrk)))//'/(m^{2})'
      call customise_label(stmp,wall_label(iwwrk))
      return
      end subroutine pswa
