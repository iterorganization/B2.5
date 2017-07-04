      subroutine my_outi(nwrite, nx, ny, iarr, first_name) !srv 19.12.07
      integer nwrite, nx, ny
      integer iarr(-1:nx,-1:ny)
      character *(*) first_name

      character *64 filename
      character(23) hlp_frm                               !srv 11.03.13
      integer ix, iy
      integer, save :: ncall = 0
      integer, save :: iappend = 0

      if (ncall.eq.0) then
        call ipgeti('b2wdat_append', iappend)
        ncall = 1
      end if
      filename = trim(first_name)//'.dat'
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
      end
