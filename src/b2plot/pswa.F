c
      subroutine pswa()
      use b2mod_types
      use b2mod_wall
      use b2mod_b2plot
      use b2mod_b2plot_values
      implicit none
      integer iw,is,lnblnk
      intrinsic min,max
      external lnblnk

      values(9)=wall_stack(1,0,iwwrk)/wall_area(1)
      values(10)=wall_stack(1,0,iwwrk)/wall_area(1)
      do is=0,nwdim(2,iwwrk)-1
        do iw=1,nbnd
          wall_stack(iw,is,iwwrk)=wall_stack(iw,is,iwwrk)/wall_area(iw)
          values(9)=min(values(9),wall_stack(iw,is,iwwrk))
          values(10)=max(values(10),wall_stack(iw,is,iwwrk))
        enddo
      enddo
      wall_label(iwwrk)=wall_label(iwwrk)(1:lnblnk(wall_label(iwwrk)))//
     1 '/(m^2)'
      return
      end
