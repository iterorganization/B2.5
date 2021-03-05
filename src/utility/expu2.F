      function expu2 (arg)
      use b2mod_types
      implicit none
      real (kind=R8) :: expu2, arg, cutll, cutlo, MACHSFR
      intrinsic dim, exp, max, log, min
      external machsfr

      cutll = 0.0_R8
      cutlo = 0.0_R8

      cutll = log(machsfr())
      cutlo = max(machsfr(),exp(cutll))

      expu2 = dim(exp(min(20.0_R8,max(arg,cutll))),cutlo)

      return
*     ------------------------------------------------------------------
      end
