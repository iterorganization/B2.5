c
      subroutine lblalongcl(along)
      use b2mod_types
      use b2mod_graphics
      implicit none
      logical along

      if(along) then
        lalong=.true.
        write (*,*) 'Labels: along Contour Lines'
      else
        lalong=.false.
        write (*,*) 'Labels: normal'
      endif
      return
      end subroutine lblalongcl
