c
      subroutine psx()
      use b2mod_types
      use b2mod_geo
      use b2mod_b2plot
      use b2mod_b2plot_values
      implicit none
      integer ix,iy,is,lnblnk
      character*256 stmp
      intrinsic min,max
      external lnblnk

      values(9)=work(-1,-1,0,iwrk)/gs(-1,-1,0)
      values(10)=work(-1,-1,0,iwrk)/gs(-1,-1,0)
      do is=0,ndim(3,iwrk)-1
        do iy=-1,ny
          do ix=-1,nx
            work(ix,iy,is,iwrk)=work(ix,iy,is,iwrk)/gs(ix,iy,0)
            values(9)=min(values(9),work(ix,iy,is,iwrk))
            values(10)=max(values(10),work(ix,iy,is,iwrk))
          enddo
        enddo
      enddo
      stmp=graphlabel(iwrk)(1:lnblnk(graphlabel(iwrk)))//'/(m^{2})'
      call customise_label(stmp,graphlabel(iwrk))
      return
      end
c
      subroutine psxperp()
      use b2mod_types
      use b2mod_geo
      use b2mod_b2plot
      use b2mod_b2plot_values
      implicit none
      integer ix,iy,is,lnblnk
      character*256 stmp
      intrinsic min,max
      external lnblnk

      values(9)=work(-1,-1,0,iwrk)/gs(-1,-1,0)/qc(-1,-1)
      values(10)=work(-1,-1,0,iwrk)/gs(-1,-1,0)/qc(-1,-1)
      do is=0,ndim(3,iwrk)-1
        do iy=-1,ny
          do ix=-1,nx
	    work(ix,iy,is,iwrk)=work(ix,iy,is,iwrk)/
     1	     gs(ix,iy,0)/qc(ix,iy)
            values(9)=min(values(9),work(ix,iy,is,iwrk))
            values(10)=max(values(10),work(ix,iy,is,iwrk))
          enddo
        enddo
      enddo
      stmp=graphlabel(iwrk)(1:lnblnk(graphlabel(iwrk)))//'/(m^{2})'
      call customise_label(stmp,graphlabel(iwrk))
      return
      end
