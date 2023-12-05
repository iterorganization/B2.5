c
      subroutine sumy()
      use b2mod_types
      use b2mod_b2plot
      use b2mod_b2plot_values
      implicit none
#ifdef WG_TODO
      integer ix,iy,is
c
c apply operator
c
      do is=0,ndim(2,iwrk)-1
        do ix=-1,nx
          do iy=0,ny
            work(ix,-1,is,iwrk)=work(ix,-1,is,iwrk)+work(ix,iy,is,iwrk)
          enddo
        enddo
      enddo
      values(9)=work(-1,-1,0,iwrk)
      values(10)=work(-1,-1,0,iwrk)
      do is=0,ndim(2,iwrk)-1
        do ix=-1,nx
          do iy=0,ny
            work(ix,iy,is,iwrk)=work(ix,-1,is,iwrk)
          enddo
          values(9)=min(values(9),work(ix,-1,is,iwrk))
          values(10)=max(values(10),work(ix,-1,is,iwrk))
        enddo
      enddo
#endif

      return
      end subroutine sumy
c
      subroutine inty()
      use b2mod_types
      use b2mod_b2plot
      use b2mod_b2plot_values
      implicit none
#ifdef WG_TODO
      integer ix,iy,is
c
c apply operator
c
      do is=0,ndim(2,iwrk)-1
        do ix=-1,nx
          do iy=1,ny
            work(ix,iy,is,iwrk)=work(ix,iy-1,is,iwrk)+
     1       work(ix,iy,is,iwrk)
          enddo
        enddo
      enddo
      values(9)=work(-1,ny,0,iwrk)
      values(10)=work(-1,ny,0,iwrk)
      do is=0,ndim(2,iwrk)-1
        do ix=-1,nx
          values(9)=min(values(9),work(ix,ny,is,iwrk))
          values(10)=max(values(10),work(ix,ny,is,iwrk))
        enddo
      enddo
#endif

      return
      end subroutine inty
