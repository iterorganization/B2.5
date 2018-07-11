      subroutine b2xpnr (nCv, ns, am, na, nirm)
      use b2mod_types
      implicit none
      integer nCv, ns
      real (kind=R8) ::
     *   am(0:ns-1), na(nCv,0:ns-1), nirm(nCv)
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
        nirm = nirm + na(:,is)/sqrt(am(is))
      enddo
*$$$  call subend ()
      return
*     ------------------------------------------------------------------
      end
