! $AD NOCHECKPOINT
      function harm (t0,t1)
      use b2mod_types
      implicit none
      real (kind=R8) :: harm, t0, t1
      intrinsic min, max

      harm = 2.0_R8*t0*t1/(t0+t1)

      return
*     ------------------------------------------------------------------
      end
