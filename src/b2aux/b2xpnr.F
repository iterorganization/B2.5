      subroutine b2xpnr (nx, ny, ns, am, na, nirm)
      use b2mod_types
      implicit none
      integer nx, ny, ns
      real (kind=R8) ::
     *   am(0:ns-1), na(-1:nx,-1:ny,0:ns-1), nirm(-1:nx,-1:ny)
*     ------------------------------------------------------------------
*     B2XPNR computes the atom density weighted by 1/sqrt(am), nirm:
*       nirm(,) = (sum is :: na(,,is)/sqrt(am(is))) .
*     ------------------------------------------------------------------
      integer n2, is
      real (kind=R8) ::
     *   t0
      intrinsic sqrt
      external sfill, saxpy
*     ------------------------------------------------------------------
*$$$  call subini ('b2xpnr')
      n2 = (nx+2)*(ny+2)
      call sfill (n2, 0.0_R8, nirm, 1)
      do is = 0, ns-1
       t0 = 1.0_R8/sqrt(am(is))
       call saxpy (n2, t0, na(-1,-1,is), 1, nirm, 1)
      enddo
*$$$  call subend ()
      return
*     ------------------------------------------------------------------
      end
