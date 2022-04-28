      subroutine get_profile(filename, npoints, positions, values)  
      use b2mod_types
      implicit none
      logical file_ok
      character*260 filename
      integer :: npoints, iy
c      real (kind=R8), allocatable :: positions(:), values(:)
      real (kind=R8) :: positions(0:npoints-1), values(0:npoints-1)
      external xerrab
      inquire(file=filename,exist=file_ok)
      if(file_ok) then
        open(99,file=filename)
        read(99,*,err=100) iy
        do iy = 0, npoints-1
          read(99,*,err=101) positions(iy), values(iy)
        end do
        close(99)
      else
        write (*,*) 'get_profile: filename =', filename
        call xerrab('file did not exist')         
      end if
      return
 100  write (*,*) 'get_profile: filename =', filename
      call XERRAB('get_profile, FILE READING ERROR')
 101  call XERRAB('get_profile, reached EOF')

      end subroutine 

      subroutine read_first_line(filename, sz)
*open and read first line, an integer value indicating the number of lines to follow
      use b2mod_types
      implicit none
      integer sz
      logical file_ok
      character*260 filename
      external xerrab

      inquire(file=filename,exist=file_ok)
      if (file_ok) then
        open(99,file=filename)
        read(99,*,err=100) sz
        close(99)
*check input
        if (sz.le.3) then
          call XERRAB ('read_first_line, more than 3 points
     .          needed for profile')
        end if
      else
        write (*,*) 'read_first_line: filename =', filename
        call xerrab('file did not exist')
      end if
      return
 100  call XERRAB('READ_FIRST_LINE, FILE READING ERROR')
      
      end subroutine read_first_line   
