! $AD NOCHECKPOINT
      function hybr (fl,df)
      use b2mod_types
      implicit none
      real (kind=R8) :: hybr, fl, df
      intrinsic max, abs, sqrt

      hybr = max(abs(df),abs(fl/2.0_R8),sqrt(df**2+fl**2/4.0_R8))

      return
*     ------------------------------------------------------------------
      end
