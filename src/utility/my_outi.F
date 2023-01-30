      subroutine my_outi(nwrite, nx, ny, iarr, first_name) !srv 19.12.07
      integer nwrite, nx, ny
      integer iarr(-1:nx,-1:ny)
      character *(*) first_name

      character *64 filename
      character(23) hlp_frm                               !srv 11.03.13
      integer ix, iy
      integer, save :: ncall = 0
      integer, save :: iappend = 0
      character(7), save :: folder                        !srv 26.09.17
      logical file_opened

      if (ncall.eq.0) then
        call ipgeti('b2wdat_append', iappend)
        ncall = 1
*    ..try to create a file in optional output folder
        inquire(unit=nwrite,opened=file_opened)
        if (file_opened) then
         write(*,fmt="(a23,i3,a20)"),
     & 'Error in my_outi: unit ',nwrite,' is already opened'
         return
        endif
        open(nwrite,file='output/__1__.dat',err=10)
*    ..folder exists: future output will be done in folder
        folder='output/'
        write(nwrite,*) 'Test file to check on output folder'
        close(nwrite)
        go to 20
*    ..folder does not exist: future output will be done in current folder
   10   folder=''
   20   continue
      end if
      filename = trim(folder)//trim(first_name)//'.dat'   !srv 26.09.17
      inquire(unit=nwrite,opened=file_opened)
      if (file_opened) then
         write(*,fmt="(a23,i3,a20)"),
     & 'Error in my_outi: unit ',nwrite,' is already opened'
        return
      endif
      if (iappend.eq.0) then
        open (nwrite,file=filename)
      else
        open (nwrite,file=filename,position='append')
      end if
      write(hlp_frm,"(a4,i4,a7)") '(5x,',nx+2,'( i4 ))'
      write (nwrite,hlp_frm) (ix,ix=-1,nx)
      write (nwrite,"(a4,i4,a7)")
      write (hlp_frm,"(a10,i4,a3)") '(1x,i3,1x,',nx+2,'i4)'
      do iy=ny,-1,-1
        write(nwrite,hlp_frm) iy,(iarr(ix,iy),ix=-1,nx)
      enddo
      close (nwrite)

      return
      end subroutine my_outi
