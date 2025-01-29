c
      subroutine clvl(onflag)
      use b2mod_types
      use b2mod_graphics
      implicit none
      logical onflag
      if(onflag) then
        contours=.true.
        write (*,*) 'contour level lines on'
      else
        contours=.false.
        write (*,*) 'contour level lines off'
      endif
      return
      end subroutine clvl
