c
      subroutine avey()
      use b2mod_types
      use b2mod_b2plot
      implicit none
      if(ndim(1,iwrk) .lt. 0 .or. ndim(2,iwrk) .lt. 0) then
        write(*,*) 'To apply this operator, a two dimensional array ',
     1   'must be specified'
        return
      endif
c
c apply operator
c
c      ndim=1
      return
      end
