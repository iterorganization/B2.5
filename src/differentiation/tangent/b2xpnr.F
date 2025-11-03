      subroutine b2xpnr (nCv, ns, na, nirm)
      use b2mod_types
      use b2mod_b2cmpa_diffv
      implicit none
      integer nCv, ns
      real (kind=R8) ::
     *   na(nCv,0:ns-1), nirm(nCv)
*     ------------------------------------------------------------------
*     B2XPNR computes the atom density weighted by 1/sqrt(am), nirm:
*       nirm() = (sum is :: na(,is)/sqrt(am(is))) .
*     ------------------------------------------------------------------
      integer is
      intrinsic sqrt
*     ------------------------------------------------------------------
*$$$  call subini ('b2xpnr')
      nirm = 0.0_R8
      do is = 0, ns-1
       if (.not.is_neutral(is)) then
        nirm = nirm + na(:,is)/sqrt(am(is))
       endif
      enddo
*$$$  call subend ()
      return
*     ------------------------------------------------------------------
      end subroutine b2xpnr
