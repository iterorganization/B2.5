c
      subroutine sumx()
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
        do iy=-1,ny
          do ix=0,nx
            work(-1,iy,is,iwrk)=work(-1,iy,is,iwrk)+work(ix,iy,is,iwrk)
          enddo
        enddo
      enddo
      values(9)=work(-1,-1,0,iwrk)
      values(10)=work(-1,-1,0,iwrk)
      do is=0,ndim(3,iwrk)-1
        do iy=-1,ny
          do ix=0,nx
            work(ix,iy,is,iwrk)=work(-1,iy,is,iwrk)
          enddo
          values(9)=min(values(9),work(-1,iy,is,iwrk))
          values(10)=max(values(10),work(-1,iy,is,iwrk))
        enddo
      enddo

      return
      end
c
      subroutine sumw()
      use b2mod_types
      use b2mod_b2plot
      use b2mod_b2plot_values
      implicit none
      integer i,is
      if(nwdim(1,iwwrk) .lt. 0) then
        write(*,*) 'To apply this operator, an array ',
     1   'must be specified'
        return
      endif
c
c apply operator
c
      do is=0,nwdim(2,iwwrk)-1
        do i= 2,nwdim(1,iwwrk)
          wall_stack(1,is,iwwrk)=
     .     wall_stack(1,is,iwwrk)+wall_stack(i,is,iwwrk)
        enddo
      enddo
      values(9)=wall_stack(1,0,iwwrk)
      values(10)=wall_stack(1,0,iwwrk)
      do is=0,nwdim(2,iwwrk)-1
        do i=2,nwdim(1,iwwrk)
          wall_stack(i,is,iwwrk)=wall_stack(1,is,iwwrk)
        enddo
        values(9)=min(values(9),wall_stack(1,is,iwwrk))
        values(10)=max(values(10),wall_stack(1,is,iwwrk))
      enddo

      return
      end
c
      subroutine intx()
      use b2mod_types
      use b2mod_b2plot
      use b2mod_indirect
      use b2mod_b2plot_values
      implicit none
      integer ix,ixl,ixr,iy,is
      intrinsic min, max
      if(ndim(1,iwrk) .le. 1 .or. ndim(2,iwrk) .le. 1) then
        write(*,*) 'To apply this operator, a two dimensional array ',
     1   'must be specified'
        return
      endif
c
c apply operator
c
      do is=0,ndim(3,iwrk)-1
        do iy=-1,ny
          if(nnreg(0).le.5 .or. iy.gt.jsep) then
            do ix = 0, nx
              ixl=leftix(ix,iy)
              ixr=rightix(ix,iy)
              if(ixl.eq.ix-1)
     1         work(ix,iy,is,iwrk)=work(ix,iy,is,iwrk)+
     1                           work(ixl,iy,is,iwrk)
              if(ixr.lt.ix) exit
            enddo
            ixl=ixr
            do ix = nx-1, ixl, -1
              ixr=rightix(ix,iy)
              if(ixr.eq.0 .or. ixr.eq.-1) exit
              if(ixr.eq.ix+1)
     1         work(ix,iy,is,iwrk)=work(ix,iy,is,iwrk)+
     1                           work(ixr,iy,is,iwrk)
              if(leftix(ix,iy).lt.ix-1) exit
            enddo
          else
            ix=0
            do while (leftix(ix,iy).eq.ix-1 .and. ix.lt.nx)
              work(ix,iy,is,iwrk)=work(ix,iy,is,iwrk)+
     1         work(leftix(ix,iy),leftiy(ix,iy),is,iwrk)
              ix=ix+1
            enddo
            ix=ix+1     
            do while (leftix(ix,iy).eq.ix-1 .and. ix.lt.nx)
              work(ix,iy,is,iwrk)=work(ix,iy,is,iwrk)+
     1         work(leftix(ix,iy),leftiy(ix,iy),is,iwrk)
              ix=ix+1
            enddo
            ixr=ix
            do while (rightix(ixr,iy).ne.nx+1)
              ixr=ixr+1
            enddo
            ix=ixr-1
            do while (rightix(ix,iy).eq.ix+1 .and. ix.gt.-1)
              work(ix,iy,is,iwrk)=work(ix,iy,is,iwrk)+
     1         work(rightix(ix,iy),rightiy(ix,iy),is,iwrk)
              ix=ix-1
            enddo
            ixl=ixr
            do while (leftix(ixl,iy).ne.-2)
              ixl=ixl+1
            enddo
            ix=ixl+1   
            do while (leftix(ix,iy).eq.ix-1 .and. ix.lt.nx)
              work(ix,iy,is,iwrk)=work(ix,iy,is,iwrk)+
     1         work(leftix(ix,iy),leftiy(ix,iy),is,iwrk)
              ix=ix+1
            enddo
            ix=ix+1     
            do while (leftix(ix,iy).eq.ix-1 .and. ix.lt.nx)
              work(ix,iy,is,iwrk)=work(ix,iy,is,iwrk)+
     1         work(leftix(ix,iy),leftiy(ix,iy),is,iwrk)
              ix=ix+1
            enddo
            ix=nx
            do while (rightix(ix,iy).eq.ix+1 .and. ix.gt.-1)
              work(ix,iy,is,iwrk)=work(ix,iy,is,iwrk)+
     1         work(rightix(ix,iy),rightiy(ix,iy),is,iwrk)
              ix=ix-1
            enddo
          endif
        enddo
      enddo
      values(9)=work(-1,-1,0,iwrk)
      values(10)=work(-1,-1,0,iwrk)
      do is=0,ndim(3,iwrk)-1
        do iy=-1,ny
          do ix=-1,nx
            values(9)=min(values(9),work(ix,iy,is,iwrk))
            values(10)=max(values(10),work(ix,iy,is,iwrk))
          enddo
        enddo
      enddo

      return
      end
c
      subroutine intw()
      use b2mod_types
      use b2mod_wall
      use b2mod_b2plot
      use b2mod_b2plot_values
      implicit none
      integer iw,iseg,is
      intrinsic min, max
      if(nwdim(1,iwwrk) .le. 1) then
        write(*,*) 'To apply this operator, an array ',
     1   'must be specified'
        return
      endif
c
c apply operator
c
      do is=0,nwdim(2,iwwrk)-1
        do iseg = 1, nsegments
          do iw=wall_begin(iseg)+1,wall_end(iseg)
            wall_stack(iw,is,iwwrk)=
     .       wall_stack(iw-1,is,iwwrk)+
     .       wall_stack(iw,is,iwwrk)
          enddo
        enddo
      enddo
      values(9)=wall_stack(1,0,iwwrk)
      values(10)=wall_stack(1,0,iwwrk)
      do is=0,nwdim(2,iwwrk)-1
        do iw=1,nbnd
          values(9)=min(values(9),wall_stack(iw,is,iwwrk))
          values(10)=max(values(10),wall_stack(iw,is,iwwrk))
        enddo
      enddo

      return
      end
c
      subroutine intwstar()
c same as intw but no account of segments
      use b2mod_types
      use b2mod_wall
      use b2mod_b2plot
      use b2mod_b2plot_values
      implicit none
      integer iw,is
      intrinsic min, max
      if(nwdim(1,iwwrk) .le. 1) then
        write(*,*) 'To apply this operator, an array ',
     1   'must be specified'
        return
      endif
c
c apply operator
c
      do is=0,nwdim(2,iwwrk)-1
        do iw = 2, nwdim(1,iwwrk)
          wall_stack(iw,is,iwwrk)=
     1     wall_stack(iw-1,is,iwwrk)+wall_stack(iw,is,iwwrk)
        enddo
      enddo
      values(9)=wall_stack(1,0,iwwrk)
      values(10)=wall_stack(1,0,iwwrk)
      do is=0,nwdim(2,iwwrk)-1
        do iw=1,nbnd
          values(9)=min(values(9),wall_stack(iw,is,iwwrk))
          values(10)=max(values(10),wall_stack(iw,is,iwwrk))
        enddo
      enddo

      return
      end
