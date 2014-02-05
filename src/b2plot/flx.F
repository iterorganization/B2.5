c
      subroutine flx(iy)
      use b2mod_types
      use b2mod_b2plot
      implicit none
      integer iy
      integer iadd
      if(ndim(1,iwrk) .lt. 0 .or. ndim(2,iwrk) .lt. 0) then
        write(*,*) 'To apply this operator, a two dimensional array ',
     1   'must be specified'
        return
      endif
c
c apply operator
c
      iadd=0
      if (linlog.le.2) iadd=2
      linlog=linlog+iadd
      call fx(iy)
      linlog=linlog-iadd
      return
      end
