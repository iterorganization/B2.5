c
      subroutine ipolate(onflag)
      use b2mod_types
      use b2mod_graphics
      implicit none
      logical onflag

      if(onflag) then
        intpol=.true.
      else
        intpol=.false.
      endif
      return
      end subroutine ipolate
