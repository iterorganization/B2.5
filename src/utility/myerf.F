      function myerf (x)
      use b2mod_types
      implicit none
      real (kind=R8) :: myerf, x
 
*     ------------------------------------------------------------------
*      B2.5 representation of intrinsic function erf which is 
*      needed for differentiation with TAPENADE        
*     ------------------------------------------------------------------
      intrinsic erf

      myerf = erf(x)    
 
      return
*     ------------------------------------------------------------------
      end

