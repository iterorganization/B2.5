      subroutine my_out(nwrite, nx, ny, fun, first_name)  !srv 19.12.07
      use b2mod_types
      integer nwrite, nx, ny
      real (kind=R8) :: fun(-1:nx,-1:ny)
      character *(*) first_name

      character *64 filename                              !srv 05.10.09
      character(24) hlp_frm                               !srv 11.03.13
      integer ix, iy
      integer, save :: ncall = 0
      integer, save :: iappend = 0
#ifdef DBG
      integer, save :: ndigits = 15
#else
      integer, save :: ndigits = 6
#endif
      character(7), save :: folder                        !srv 26.09.17

      if (ncall.eq.0) then
        call ipgeti('b2wdat_append', iappend)
        call ipgeti('my_out_digits', ndigits)
        call xertst(0.le.ndigits, 'invalid ndigits value')
        ncall = 1
*    ..try to create a file in optional output folder
        open(70,file='output/__1__.dat',err=10)
*    ..folder exists: future output will be done in folder
        folder='output/'
        write(70,*) 'Test file to check on output folder'
        close(70)
        go to 20
*    ..folder does not exist: future output will be done in current folder
   10   folder=''
   20   continue
      end if
      filename = trim(folder)//trim(first_name)//'.dat'   !srv 26.09.17
      if (iappend.eq.0) then
        open (nwrite,file=filename)
      else
        open (nwrite,file=filename,position='append')
      end if
      write(hlp_frm,"(a4,i4,a1,i2,a9)")
     . '(5x,',nx+2,'(',ndigits,'x,i3,5x))'
      write (nwrite,hlp_frm) (ix,ix=-1,nx)
      write(hlp_frm,"(a10,i4,a3,i2.2,a1,i2.2,a2)")
     . '(1x,i3,1x,',nx+2,'(es',ndigits+8,'.',ndigits,'))'
      do iy=ny,-1,-1
        write(nwrite,hlp_frm) iy, (fun(ix,iy), ix = -1, nx)
      enddo
      close (nwrite)

      return
      end
