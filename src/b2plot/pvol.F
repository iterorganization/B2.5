c
      subroutine pvol()
      use b2mod_types
      use b2mod_geo
      use b2mod_b2plot
      use b2mod_b2plot_values
      implicit none
      integer ix,iy,is,lnblnk
      character*256 stmp
      external lnblnk

      values(9)=work(-1,-1,0,iwrk)/vol(-1,-1)
      values(10)=work(-1,-1,0,iwrk)/vol(-1,-1)
      do is=0,ndim(3,iwrk)-1
        do iy=-1,ny
          do ix=-1,nx
            work(ix,iy,is,iwrk)=work(ix,iy,is,iwrk)/vol(ix,iy)
            values(9)=min(values(9),work(ix,iy,is,iwrk))
            values(10)=max(values(10),work(ix,iy,is,iwrk))
          enddo
        enddo
      enddo
      stmp=graphlabel(iwrk)(1:lnblnk(graphlabel(iwrk)))//'/(m^{3})'
      call customise_label(stmp,graphlabel(iwrk))
      return
      end
