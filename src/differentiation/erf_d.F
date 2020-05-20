      function erf_d (x, xd, res)
      use b2mod_types
      use b2mod_constants
      implicit none
      real (kind=R8) :: erf_d, x, xd, res
 
*     ------------------------------------------------------------------
*      Forward differentiation of intrinsic function erf which is  
*      needed for differentiation with TAPENADE       
*     ------------------------------------------------------------------
      intrinsic sqrt, exp, erf

      res = erf(x)

      erf_d = xd*2.0_R8/sqrt(pi)*exp(-(x**2))
 
      return
*     ------------------------------------------------------------------
      end

