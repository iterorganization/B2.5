      subroutine erf_b (x, xb, resb)
      use b2mod_types
      use b2mod_constants
      implicit none
      real (kind=R8) :: x, xb, resb
 
*     ------------------------------------------------------------------
*      Adjoint differentiation of intrinsic function erf which is  
*      needed for differentiation with TAPENADE       
*      Example: yy = f(xx) = erf(xx)
*      xb = xb + yb*df(xx) = xb + yb*2/sqrt(pi)*exp(-xx^2)
*     ------------------------------------------------------------------
      intrinsic sqrt, exp, erf

      xb = xb + resb*2.0_R8/sqrt(pi)*exp(-(x**2))
 
      return
*     ------------------------------------------------------------------
      end subroutine erf_b

