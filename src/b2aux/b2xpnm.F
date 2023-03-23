      subroutine b2xpnm (nCv, ns, rz2, am, na, ne2m)
      use b2mod_types
      implicit none
      integer nCv, ns
      real (kind=R8) ::
     *  rz2(nCv,0:ns-1), am(0:ns-1), na(nCv,0:ns-1),
     &  ne2m(nCv)
*     ------------------------------------------------------------------
*     B2XPNM computes the atom density weighted by rz2/am, ne2m:
*       ne2m() = (sum is :: (rz2(,is)/am(is))*na(,is)) .
*     ------------------------------------------------------------------
      integer is
      external sfill
*     ------------------------------------------------------------------
*$$$  call subini ('b2xpnm')
      call sfill (nCv, 0.0_R8, ne2m, 1)
      do is = 0, ns-1
       ne2m = ne2m + (rz2(:,is)/am(is)) * na(:,is)
      enddo
*$$$  call subend ()
      return
*     ------------------------------------------------------------------
      end subroutine b2xpnm
