! $AD NOCHECKPOINT
      function trimg (aa,bb,cc)
      use b2mod_types
      implicit none
      real (kind=R8) :: trimg, aa, bb, cc
      intrinsic min, max

      trimg = max(aa,min(bb,cc))
 
      return
*     ------------------------------------------------------------------
      end
