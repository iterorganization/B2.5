! $AD NOCHECKPOINT
      function med (t0, t1, t2)
      use b2mod_types
      implicit none
      real (kind=R8) :: t0, t1, t2, med
      intrinsic min, max

      med = min(max(t0,t1),max(t1,t2),max(t2,t0))
 
      return
*     ------------------------------------------------------------------
      end
