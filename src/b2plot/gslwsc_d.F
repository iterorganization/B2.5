
      subroutine gslwsc_d(dw)
      use b2mod_types
      real (kind=R8) :: dw
      real (kind=R4) :: w
      w=real(dw,R4)
#ifndef PLPLOT
      call lastd
      call gslwsc(w)
#endif
      return
      end
