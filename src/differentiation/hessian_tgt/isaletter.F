c
      logical function isaletter(c)
      use b2mod_types
      character*1 c
      isaletter=(c.ge.'A' .and. c.le.'Z') .or.
     &          (c.ge.'a' .and. c.le.'z')
      return
      end function isaletter
