! $AD NOCHECKPOINT
      function upwind (fl,t0,t1)
      use b2mod_types
      implicit none
      real (kind=R8) :: upwind, fl, t0, t1
      intrinsic min, max

      upwind = max(fl,0.0_R8)*t0+min(fl,0.0_R8)*t1

      return
*     ------------------------------------------------------------------
      end
