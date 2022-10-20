      subroutine b2xpmc (nCv, ub, pr, ro, mcb)
      use b2mod_types
      implicit none
      integer nCv
      real (kind=R8) ::
     *   ub(nCv), pr(nCv), ro(nCv), mcb(nCv)
*     ------------------------------------------------------------------
*     B2XPMC computes the mach number of the parallel flow.
*     ------------------------------------------------------------------
      intrinsic sqrt
*     ------------------------------------------------------------------
*$$$  call subini ('b2xpmc')
      mcb = ub/sqrt(pr/ro)
*$$$  call subend ()
      return
*     ------------------------------------------------------------------
      end subroutine b2xpmc
