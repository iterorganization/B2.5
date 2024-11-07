c
      subroutine separatrixon(onflag)
      use b2mod_types
      use b2mod_graphics
      implicit none
      logical onflag

      if(onflag) then
        sptrix=.not.sptrix
        write (*,*) 'separatrix state changed'
      else
        sptrix=.false.
        write (*,*) 'separatrix off'
      endif
      return
      end subroutine separatrixon
