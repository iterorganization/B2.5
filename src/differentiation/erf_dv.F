      subroutine erf_dv (x, xd, res, resd, nbdirs)
      use b2mod_types
      use b2mod_constants
      use diffsizes
      implicit none
      integer nbdirs
      real (kind=R8) :: x, res
      real (kind=R8) :: xd(nbdirsmax), resd(nbdirsmax)
*     ------------------------------------------------------------------
*      Forward vector differentiation of intrinsic function erf which is  
*      needed for differentiation with TAPENADE       
*     ------------------------------------------------------------------
      integer nd
      intrinsic sqrt, exp, erf

      res = erf(x)
      do nd=1,nbdirs
        resd(nd) = xd(nd)*2.0_R8/sqrt(pi)*exp(-(x**2))
      end do

      return
*     ------------------------------------------------------------------
      end subroutine erf_dv

