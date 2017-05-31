      subroutine my_out(nwrite, nx, ny, fun, first_name)  !srv 19.12.07
      use b2mod_types
      integer nwrite, nx, ny
      real (kind=R8) :: fun(-1:nx,-1:ny)
      character *(*) first_name

      character *64 filename                              !srv 05.10.09
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
      write(hlp_frm,"(a4,i4,a11)") '(5x,',nx+2,'(6x,i3,5x))'
      write (nwrite,hlp_frm) (ix,ix=-1,nx)
      write(hlp_frm,"(a10,i4,a9)") '(1x,i3,1x,',nx+2,'(es14.6))'
      do iy=ny,-1,-1
        write(nwrite,hlp_frm) iy, (fun(ix,iy), ix = -1, nx)
      enddo
      close (nwrite)

      return
      end
