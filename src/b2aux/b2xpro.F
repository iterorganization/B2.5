      subroutine b2xpro (nx, ny, ns, mp, am, na, ro)
      use b2mod_types
      use b2mod_external
      implicit none
      integer nx, ny, ns
      real (kind=R8) ::
     *   mp, am(0:ns-1), na(-1:nx,-1:ny,0:ns-1), ro(-1:nx,-1:ny)
*     ------------------------------------------------------------------
*     B2XPRO computes the mass density, ro:
*       ro(,) = (sum is :: am(is)*mp*na(,,is)) .
*     ------------------------------------------------------------------
      integer n2, is
      external sfill, b2saxpy
*     ------------------------------------------------------------------
*$$$  call subini ('b2xpro')
      n2 = (nx+2)*(ny+2)
      call sfill (n2, 0.0_R8, ro, 1)
      do is = 0, ns-1
       call b2saxpy (n2, mp*am(is), na(-1,-1,is), 1, ro, 1)
      enddo
      do is = 0, ns_ext-1
       call b2saxpy (n2, mp*am_ext(is), na_ext(-1,-1,is), 1, ro, 1)
      enddo
*$$$  call subend ()
      return
*     ------------------------------------------------------------------
      end subroutine b2xpro
