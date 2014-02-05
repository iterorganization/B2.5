c
      subroutine skipblanks(textline,iun,ipos)
      use b2mod_types
      integer iun, ipos
      character*500 textline
      character*1 c
   10 if(ipos .gt. 500) call readline(textline,iun,ipos)
      c=textline(ipos:ipos)
      if(c .eq. ' ') then
        ipos=ipos+1
        goto 10
      endif
      return
      end
