c
      function isnegative(textline,iun,ipos)
      use b2mod_types
      integer iun, ipos
      character*500 textline
      character*1 c
      logical isnegative, isplusorminus
      c=textline(ipos:ipos)
      isnegative=.false.
      if(isplusorminus(c)) then
        if(c .eq. '-') isnegative=.true.
        ipos=ipos+1
        if(ipos .gt. 500) call readline(textline,iun,ipos)
      endif
      return
      end function isnegative
