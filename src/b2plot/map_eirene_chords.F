c
      subroutine map_eirene_chords(filename)

      use b2mod_types
      use b2mod_geo
      use b2mod_b2plot
      use b2mod_indirect
      use b2mod_b2plot_values
#ifdef B25_EIRENE
      use eirmod_ctrig
#endif
      implicit none

      character*(*) filename

#ifdef B25_EIRENE
      real (kind=R8) :: wrkarr(ntrii)
      character*80 dataname,ccol
      logical open_file,found
      real (kind=R8) ::
     *   rcol,r,ro,
     1   x,x1,x2,xo,y,y1,y2,yo,z,z1,z2,zo
      integer npoints,i,it
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
      wrkarr(1:ntrii)=0.0_R8
      call init_page(
     1 dataname(1:lnblnk(dataname))//' '//ccol(1:lnblnk(ccol)),
     1 values(5),values(6),values(7),values(8),1,.true.,.true.,.true.)
      do i=1,npoints
        xo=x
        yo=y
        zo=z
        ro=r
        x=x1*real(npoints-i,R8)/real(npoints,R8)+
     .    x2*real(i,R8)/real(npoints,R8)
        y=y1*real(npoints-i,R8)/real(npoints,R8)+
     .    y2*real(i,R8)/real(npoints,R8)
        z=z1*real(npoints-i,R8)/real(npoints,R8)+
     .    z2*real(i,R8)/real(npoints,R8)
        r=sqrt(x**2+y**2)
        call chord_tria_find(ntrii,nrknot,necke,xtrian,ytrian,
     1     (x+xo)/2.0_R8,(y+yo)/2.0_R8,(z+zo)/2.0_R8,
     2     found,it)
        if(found) then
          wrkarr(it)=wrkarr(it)+1.0_R8/real(npoints)
        endif
        call lined_d(ro,zo,r,z)
      enddo
      call plottripatch(nx,ny,ntrii,nrknot,0,
     2 xtrian,ytrian,necke,wrkarr,
     1 0.0_R8,1.0_R8,linlog)
      call data_copy(0)
      call additional_elements(nx,ny,.true.)
      call end_page
      goto 10
   30 write(*,*) 'Error in reading ',filename(1:lnblnk(filename))
      return
   40 continue
      close(90)
#endif
      return
      end subroutine map_eirene_chords
