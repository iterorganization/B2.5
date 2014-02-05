c
      logical function isadigit(c)
      use b2mod_types
      character*1 c
      isadigit=c .ge. '0' .and. c .le. '9'
      return
      end
