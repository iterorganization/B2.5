      subroutine hbout(a1,nx,ny,ns,y1,y1a,y2,y2a,z1,z1a,z2,z2a,
     &     z3,z3a,z4,z4a)
      use b2mod_types
      integer nx, ny, ns, ix, iy, is
      real (kind=R8) ::
     *   z1(-1:nx,-1:ny,0:ns-1), y1(-1:nx,-1:ny), y2(-1:nx,-1:ny),
     &  z2(-1:nx,-1:ny,0:ns-1), z3(-1:nx,-1:ny,0:ns-1),
     &  z4(-1:nx,-1:ny,0:ns-1),
     &  ma, mi, x(-1:nx)
      character*(*) :: y1a,y2a,z1a,z2a,z3a,z4a
      character*12 :: a1
      open(99,file=trim(a1))

      write (99,*) y1a,'start'
      ma=0.
      mi=0.
      do iy=-1,ny
         do ix=-1,nx
           x(ix)=y1(ix,iy)
           ma=MAX(ma,x(ix))
           mi=MIN(mi,x(ix))
         enddo
         write (99,*) iy,ix,(x(ix),ix=-1,nx)
      enddo
      write (99,*) y1a,'ende:',mi,ma

      write (99,*) y2a,'start'
      ma=0.
      mi=0.
      do iy=-1,ny
         do ix=-1,nx
           x(ix)=y2(ix,iy)
           ma=MAX(ma,x(ix))
           mi=MIN(mi,x(ix))
         enddo
         write (99,*) iy,ix,(x(ix),ix=-1,nx)
      enddo
      write (99,*) y2a,'ende:',mi,ma

      write (99,*) z1a,'start'
      ma=0.
      mi=0.
      do is=0,ns-1
      do iy=-1,ny
         do ix=-1,nx
           x(ix)=z1(ix,iy,is)
           ma=MAX(ma,x(ix))
           mi=MIN(mi,x(ix))
         enddo
         write (99,*) is,iy,ix,(x(ix),ix=-1,nx)
      enddo
      enddo
      write (99,*) z1a,'ende:',mi,ma

      write (99,*) z2a,'start'
      ma=0.
      mi=0.
      do is=0,ns-1
      do iy=-1,ny
         do ix=-1,nx
           x(ix)=z2(ix,iy,is)
           ma=MAX(ma,x(ix))
           mi=MIN(mi,x(ix))
         enddo
         write (99,*) is,iy,ix,(x(ix),ix=-1,nx)
      enddo
      enddo
      write (99,*) z2a,'ende:',mi,ma

      write (99,*) z3a,'start'
      ma=0.
      mi=0.
      do is=0,ns-1
      do iy=-1,ny
         do ix=-1,nx
           x(ix)=z3(ix,iy,is)
           ma=MAX(ma,x(ix))
           mi=MIN(mi,x(ix))
         enddo
         write (99,*) is,iy,ix,(x(ix),ix=-1,nx)
      enddo
      enddo
      write (99,*) z3a,'ende:',mi,ma

      write (99,*) z4a,'start'
      ma=0.
      mi=0.
      do is=0,ns-1
      do iy=-1,ny
         do ix=-1,nx
           x(ix)=z4(ix,iy,is)
           ma=MAX(ma,x(ix))
           mi=MIN(mi,x(ix))
         enddo
         write (99,*) is,iy,ix,(x(ix),ix=-1,nx)
      enddo
      enddo
      write (99,*) z4a,'ende:',mi,ma

      close (99)

      return
      end subroutine hbout
