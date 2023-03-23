c
      logical function isplusorminus(c)
      use b2mod_types
      character*1 c
      isplusorminus=c .eq. '+' .or. c .eq. '-'
      return
      end function isplusorminus
