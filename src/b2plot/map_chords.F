c
      subroutine map_chords(filename)

      use b2mod_types
      use b2mod_geo
      use b2mod_b2plot
      use b2mod_indirect
      use b2mod_b2plot_values
      implicit none

      character*(*) filename

      real (kind=R8) ::
     *   mach(-1:nx,-1:ny)
      character*80 dataname,ccol
      logical open_file,found
      real (kind=R8) ::
     *   rcol,r,ro,
     1   x,x1,x2,xo,y,y1,y2,yo,z,z1,z2,zo
      integer npoints,i,ix,iy
      integer lnblnk
      external lnblnk

      if(.not.open_file(90,filename,
     1 'data.local/chords|data/chords/')) return
      read(90,*,err=30,end=30) dataname,npoints
      npoints=abs(npoints)

   10 read(90,*,err=30,end=40) x1,y1,z1,x2,y2,z2,rcol
      write(ccol,'(1p,g12.5)') rcol
      x=x1
      y=y1
      z=z1
      r=sqrt(x**2+y**2)
      do ix=-1,nx
        do iy=-1,ny
          mach(ix,iy)=0.0_R8
        enddo
      enddo
      call init_page(
     1 dataname(1:lnblnk(dataname))//' '//ccol(1:lnblnk(ccol)),
     1 values(5),values(6),values(7),values(8),1,.true.,.true.,.true.)
      do i=1,npoints
        xo=x
        yo=y
        zo=z
        ro=r
        x=x1*real(npoints-i)/real(npoints)+x2*real(i)/real(npoints)
        y=y1*real(npoints-i)/real(npoints)+y2*real(i)/real(npoints)
        z=z1*real(npoints-i)/real(npoints)+z2*real(i)/real(npoints)
        r=sqrt(x**2+y**2)
        call chord_find(nx,ny,crx,cry,
     1     (x+xo)/2.0_R8,(y+yo)/2.0_R8,(z+zo)/2.0_R8,
     2     found,ix,iy)
        if(found) then
          mach(ix,iy)=mach(ix,iy)+1.0_R8/real(npoints)
        endif
        call lined_d(ro,zo,r,z)
      enddo
      call plotpatch(nx,ny,0,crx,cry,region,mach,
     1 0.0_R8,1.0_R8,.true.,linlog,1)
      call data_copy(0)
      call additional_elements(nx,ny,.true.)
      call end_page
      goto 10
   30 write(*,*) 'Error in reading ',filename(1:lnblnk(filename))
      return
   40 continue
      close(90)
      return
      end
