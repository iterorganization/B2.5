c
      subroutine psy()
      use b2mod_types
      use b2mod_geo
      use b2mod_b2plot
      use b2mod_b2plot_values
      implicit none
#ifdef WG_TODO
      integer ix,iy,is
      intrinsic min,max
#endif
      integer lnblnk
      character*256 stmp
      external lnblnk

#ifdef WG_TODO
      values(9)=work(-1,-1,0,iwrk)/gs(-1,-1,1)
      values(10)=work(-1,-1,0,iwrk)/gs(-1,-1,1)
      do is=0,ndim(2,iwrk)-1
        do iy=-1,ny
          do ix=-1,nx
            work(ix,iy,is,iwrk)=work(ix,iy,is,iwrk)/gs(ix,iy,1)
            values(9)=min(values(9),work(ix,iy,is,iwrk))
            values(10)=max(values(10),work(ix,iy,is,iwrk))
          enddo
        enddo
      enddo
#endif
      stmp=graphlabel(iwrk)(1:lnblnk(graphlabel(iwrk)))//'/(m^{2})'
      call customise_label(stmp,graphlabel(iwrk))
      return
      end subroutine psy
