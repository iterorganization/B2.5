      subroutine b2xprz (nx, ny, ns, mp, am, na, rz)
      use b2mod_types
      use b2mod_b2cmpa
     * , only : is_neutral
      use b2mod_external
      implicit none
      integer nx, ny, ns
      real (kind=R8) ::
     *   mp, am(0:ns-1), na(-1:nx,-1:ny,0:ns-1), rz(-1:nx,-1:ny)
*     ------------------------------------------------------------------
*     B2XPRZ computes the mass density of the charged particles, rz:
*       rz(,) = (sum is :: am(is)*mp*na(,,is) for za(is).ne.0) .
*     ------------------------------------------------------------------
      integer n2, is
      external sfill, b2saxpy
*     ------------------------------------------------------------------
*$$$  call subini ('b2xprz')
      n2 = (nx+2)*(ny+2)
      call sfill (n2, 0.0_R8, rz, 1)
      do is = 0, ns-1
       if (.not.is_neutral(is)) then
        call b2saxpy (n2, mp*am(is), na(-1,-1,is), 1, rz, 1)
       endif
      enddo
      do is = 0, ns_ext-1
       if (.not.is_ext_neutral(is)) then
        call b2saxpy (n2, mp*am_ext(is), na_ext(-1,-1,is), 1, rz, 1)
       endif
      enddo
*$$$  call subend ()
      return
*     ------------------------------------------------------------------
      end subroutine b2xprz
