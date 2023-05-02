      integer function lnblnk(string)
      implicit none
      integer i
      intrinsic len
c
c returns the position of the last non-blank character in "string"
c
      character*(*) string
      do i=len(string),1,-1
        if(string(i:i).ne.' ') then
          lnblnk=i
          return
        endif
      enddo
      lnblnk=1
      return
      end function lnblnk
c
