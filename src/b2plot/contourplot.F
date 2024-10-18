c
      subroutine contourplot(onflag)
      use b2mod_types
      use b2mod_graphics
      implicit none
      logical onflag

      if(onflag) then
        contplot=.true.
        write (*,*) 'contourplot on'
      else
        contplot=.false.
        write (*,*) 'contourplot off'
      endif
      return
      end subroutine contourplot
