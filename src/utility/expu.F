! $AD NOCHECKPOINT
      function expu (arg)
      use b2mod_types
      implicit none
      real (kind=R8) :: expu, arg, cutll, cutlo, MACHSFR
      intrinsic dim, exp, max, log
      external machsfr

      cutll = 0.0_R8
      cutlo = 0.0_R8

      cutll = log(machsfr())
      cutlo = max(machsfr(),exp(cutll))

      expu = dim(exp(max(arg,cutll)),cutlo)
 
      return
*     ------------------------------------------------------------------
      end
