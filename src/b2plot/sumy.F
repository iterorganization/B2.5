c
      subroutine sumy()
      use b2mod_types
      use b2mod_b2plot
      use b2mod_b2plot_values
      implicit none
      integer ix,iy,is
      if(ndim(1,iwrk) .lt. 0 .or. ndim(2,iwrk) .lt. 0) then
        write(*,*) 'To apply this operator, a two dimensional array ',
     1   'must be specified'
        return
      endif
c
c apply operator
c
      do is=0,ndim(3,iwrk)-1
        do ix=-1,nx
          do iy=0,ny
            work(ix,-1,is,iwrk)=work(ix,-1,is,iwrk)+work(ix,iy,is,iwrk)
          enddo
        enddo
      enddo
      values(9)=work(-1,-1,0,iwrk)
      values(10)=work(-1,-1,0,iwrk)
      do is=0,ndim(3,iwrk)-1
        do ix=-1,nx
          do iy=0,ny
            work(ix,iy,is,iwrk)=work(ix,-1,is,iwrk)
          enddo
          values(9)=min(values(9),work(ix,-1,is,iwrk))
          values(10)=max(values(10),work(ix,-1,is,iwrk))
        enddo
      enddo

      return
      end subroutine sumy
c
      subroutine inty()
      use b2mod_types
      use b2mod_b2plot
      use b2mod_b2plot_values
      implicit none
      integer ix,iy,is
      if(ndim(1,iwrk) .lt. 0 .or. ndim(2,iwrk) .lt. 0) then
        write(*,*) 'To apply this operator, a two dimensional array ',
     1   'must be specified'
        return
      endif
c
c apply operator
c
      do is=0,ndim(3,iwrk)-1
        do ix=-1,nx
          do iy=0,ny
            work(ix,iy,is,iwrk)=work(ix,iy-1,is,iwrk)+
     1       work(ix,iy,is,iwrk)
          enddo
        enddo
      enddo
      values(9)=work(-1,ny,0,iwrk)
      values(10)=work(-1,ny,0,iwrk)
      do is=0,ndim(3,iwrk)-1
        do ix=-1,nx
          values(9)=min(values(9),work(ix,ny,is,iwrk))
          values(10)=max(values(10),work(ix,ny,is,iwrk))
        enddo
      enddo

      return
      end subroutine inty
