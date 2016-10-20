      function scal10 (y)
      use b2mod_types
      implicit none
      real (kind=R8) ::
     *   scal10
      real (kind=R8) ::
     *   y
*     ------------------------------------------------------------------
*     SCAL10 is meant to compute a suitable scale for graphical output.
*     scal10(y) returns the least in magnitude real*8 value r that
*     satisfies abs(y).le.abs(r) and that has a decimal representation
*     of the form  r = 0  or  r = sgn(y)*i0*10.0**i1 ,
*     where i0 is 2 or 5 or 10, and i1 is an integer.
*     ------------------------------------------------------------------
      integer i0, i1
      real (kind=R8) ::
     *   t0
      intrinsic log10, abs, nint, sign
*     ------------------------------------------------------------------
      if (y.eq.0.0_R8) then
       scal10 = 0.0_R8
      else
       i1 = nint(log10(abs(y))-0.5_R8)
       t0 = abs(y)*10.0_R8**(-i1)
*      (1.le.t0.le.10, to machine accuracy)
       if (t0.le.2.0_R8) then
        i0 = 2
       else if (t0.le.5.0_R8) then
        i0 = 5
       else
        i0 = 10
       endif
       scal10 = sign(1.0_R8,y)*i0*10.0_R8**i1
      endif
      return
*     ------------------------------------------------------------------
      end
