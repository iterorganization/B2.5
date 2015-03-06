      subroutine b2xpni (nx, ny, ns, na, ni)
      use b2mod_types
      use b2mod_constants
      use b2mod_b2cmpa
      implicit none
      integer nx, ny, ns
      real (kind=R8) ::
     *   na(-1:nx,-1:ny,0:ns-1), ni(-1:nx,-1:ny,0:1)
*   ..common blocks
*     ------------------------------------------------------------------
*     B2XPNI computes the total atom density, ni0:
*       ni(,,0) = (sum is :: na(,,is)) .
*     and the total ion density, ni1:
*       ni(,,1) = (sum is :: na(,,is) if is_neutral(is).eqv..false.) .
*     ------------------------------------------------------------------
      integer n2, is
      external sfill, b2saxpy
*     ------------------------------------------------------------------
*$$$  call subini ('b2xpni')
      n2 = (nx+2)*(ny+2)
      call sfill (n2*2, 0.0_R8, ni, 1)
      do is = 0, ns-1
       call b2saxpy (n2, 1.0_R8, na(-1,-1,is), 1, ni, 1)
       if(.not.is_neutral(is))
     1 call b2saxpy (n2, 1.0_R8, na(-1,-1,is), 1, ni(-1,-1,1), 1)
      enddo
*$$$  call subend ()
      return
*     ------------------------------------------------------------------
      end
