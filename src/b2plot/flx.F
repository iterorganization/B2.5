c
      subroutine flx(iy)
      use b2mod_types
      use b2mod_b2plot
      implicit none
      integer iy
      integer iadd
c
c apply operator
c
      iadd=0
      if (linlog.le.2) iadd=2
      linlog=linlog+iadd
      call fx(iy)
      linlog=linlog-iadd
      return
      end subroutine flx
