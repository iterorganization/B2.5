      subroutine b2xpnm (nx, ny, ns, rz2, am, na, ne2m)
      use b2mod_types
      implicit none
      integer nx, ny, ns
      real (kind=R8) ::
     *  rz2(-1:nx,-1:ny,0:ns-1), am(0:ns-1), na(-1:nx,-1:ny,0:ns-1),
     &  ne2m(-1:nx,-1:ny)
*     ------------------------------------------------------------------
*     B2XPNM computes the atom density weighted by rz2/am, ne2m:
*       ne2m(,) = (sum is :: (rz2(,,is)/am(is))*na(,,is)) .
*     ------------------------------------------------------------------
      integer n2, is
      external sfill
*     ------------------------------------------------------------------
*$$$  call subini ('b2xpnm')
      n2 = (nx+2)*(ny+2)
      call sfill (n2, 0.0_R8, ne2m, 1)
      do is = 0, ns-1
       ne2m(-1:nx,-1:ny) = ne2m(-1:nx,-1:ny) +
     .  (rz2(-1:nx,-1:ny,is)/am(is)) * na(-1:nx,-1:ny,is)
      enddo
*$$$  call subend ()
      return
*     ------------------------------------------------------------------
      end
