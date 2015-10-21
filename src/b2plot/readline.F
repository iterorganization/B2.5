c
      subroutine readline(textline,iun,ipos)
      use b2mod_types
      implicit none
      integer iun, ipos
      character*500 textline
      logical is_comment
      external is_comment
   10 read(iun,'(a500)',end=100) textline
      if(is_comment(textline)) goto 10
      ipos=1
      return
  100 textline='eof'
      ipos=1
      return
      end
