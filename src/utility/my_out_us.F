      subroutine my_out_us(nwrite, nCv, fun, first_name)
cwdk  This routine is intended to write cell-centered fields on an
c     unstructured grid similar to the functionality of the my_out-routine
c     for "original" B2.5 grids.
c     File names will be appended with "*_us.dat".
      use b2mod_types
      integer nwrite, nCv
      real (kind=R8) :: fun(1:nCv)
      character *(*) first_name

      character *64 filename
      character(23) hlp_frm
      integer ic
      integer, save :: ncall = 0
      integer, save :: iappend = 0
      logical ex
      character(7), save :: folder

      if (ncall.eq.0) then
        call ipgeti('b2wdat_append', iappend)
        ncall = 1
*    ..try to create a file in optional output folder
        open(70,file='output/__1__us.dat',err=10)
*    ..folder exists: future output will be done in folder
        folder='output/'
        write(70,*) 'Test file to check on output folder'
        close(70)
        go to 20
*    ..folder does not exist: future output will be done in current folder
 10     folder=''
 20     continue
      end if
      filename = trim(folder)//trim(first_name)//'_us.dat'
      if (iappend.eq.0) then
        open (nwrite,file=filename)
      else
        open (nwrite,file=filename,position='append')
      end if
      write(hlp_frm,"(a4,i4,a11)") '(5x,',nCv,'(6x,i3,5x))'
      write (nwrite,hlp_frm) (ic,ic=1,nCv)
      write (hlp_frm,"(a10,i4,a9)") '(5x,',nCv,'(es14.6))'
      write (nwrite,hlp_frm) (fun(ic), ic = 1, nCv)
      close (nwrite)

      return
      end
