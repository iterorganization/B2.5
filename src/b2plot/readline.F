c
      subroutine readline(textline,iun,ipos)
      use b2mod_types
      character*500 textline
   10 read(iun,'(a500)',end=100) textline
      if(is_comment(text_line)) goto 10
      ipos=1
      return
  100 textline='eof'
      ipos=1
      return
      end
