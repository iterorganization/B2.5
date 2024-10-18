c
      subroutine plot_chords(filename)
      use b2mod_types
      use b2mod_b2plot
      use b2mod_b2plot_values
      use b2mod_chord_shift
      implicit none
      character*(*) filename
      character*80 dataname
      logical open_file,lcol
      real (kind=R8) ::
     1  rcol,x,y,x1,y1,z1,x2,y2,z2,r,z,ro,zo
      integer icol,npoints,npointsl,i

      if(.not.open_file(90,filename,
     1 'data.local/chords|data/chords')) return
      read(90,*,err=30,end=30) dataname,npoints
      lcol=npoints .lt. 0
      npoints=abs(npoints)

      call init_page(trim(dataname),
     1 values(5),values(6),values(7),values(8),1,.true.,.true.,.true.)
      call additional_elements(nx,ny,.true.)

      if(lwidth .ne. 1.0_R8) call gslwsc_d(lwidth)

   10 read(90,*,err=30,end=40) x1,y1,z1,x2,y2,z2,rcol
      npointsl=npoints
      if(
     1 ((x1 .eq. x2) .and. (x1 .eq. 0.0_R8)) .or. 
     2 ((y1 .eq. y2) .and. (y1 .eq. 0.0_R8))) then
        do while(npointsl .gt. 10)
          npointsl=npointsl/2
        enddo
      endif
      x=x1
      y=y1
      z=z1+zshift
      r=sqrt(x**2+y**2)+rshift
      if(lcol) then
        call lastd
        if(color_set) then
          icol=nint(rcol)
          call gsplci(icol)
        else
          call gslwsc_d(rcol)
        endif
      endif
      do i=1,npointsl
        ro=r
        zo=z
        x=x1*real(npointsl-i,R8)/real(npointsl,R8)+
     .    x2*real(i,R8)/real(npointsl,R8)
        y=y1*real(npointsl-i,R8)/real(npointsl,R8)+
     .    y2*real(i,R8)/real(npointsl,R8)
        z=z1*real(npointsl-i,R8)/real(npointsl,R8)+
     .    z2*real(i,R8)/real(npointsl,R8)
        r=sqrt(x**2+y**2)
        r=r+rshift
        z=z+zshift
        call lined_d(ro,zo,r,z)
      enddo
      if(lcol) then
        call lastd
        if(color_set) then
          call gsplci(1)
        else
          call gslwsc_d(1.0_R8)
        endif
      endif
      goto 10
   30 write(*,*) 'Error in reading ',trim(filename)
      return
   40 continue
      close(90)
      if(lwidth .ne. 1.0_R8) call gslwsc_d(1.0_R8)
      call end_page
      return
      end subroutine plot_chords
