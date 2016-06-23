      subroutine b2xppe (nx, ny, ne, te, pe)
      use b2mod_types
      implicit none
      integer nx, ny
      real (kind=R8) ::
     *  ne(-1:nx,-1:ny), 
     &  te(-1:nx,-1:ny), pe(-1:nx,-1:ny)
*     ------------------------------------------------------------------
*     B2XPPE computes the electron pressure, pe = ne*te.
*     ------------------------------------------------------------------
      integer n2
      external sfill
*     ------------------------------------------------------------------
*$$$  call subini ('b2xppe')
      n2 = (nx+2)*(ny+2)
      call sfill (n2, 0.0_R8, pe, 1)
      pe = ne*te
*$$$  call subend ()
      return
*     ------------------------------------------------------------------
      end subroutine b2xppe
