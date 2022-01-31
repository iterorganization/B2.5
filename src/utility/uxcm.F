! $AD NOCHECKPOINT
      function uxcm (t0, t1)
      use b2mod_types
      implicit none
      real (kind=R8) :: uxcm, t0, t1
      intrinsic sign, max, abs

      uxcm = sign(max(abs(t0),abs(t1)),t0+t1)

      return
*     ------------------------------------------------------------------
      end
