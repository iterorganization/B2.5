c
      subroutine fly(ix)
      use b2mod_types
      use b2mod_b2plot
      implicit none
      integer ix
      integer iadd
c
c apply operator
c
      iadd=0
      if (linlog.le.2) iadd=2
      linlog=linlog+iadd
      call fy(ix)
      linlog=linlog-iadd
      return
      end subroutine fly
