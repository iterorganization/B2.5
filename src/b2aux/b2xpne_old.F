      subroutine b2xpne_old (nx, ny, ns, rza, na, ne_ext, ne)
      use b2mod_types
      implicit none
      integer nx, ny, ns
      real (kind=R8) ::
     *  rza(-1:nx,-1:ny,0:ns-1), na(-1:nx,-1:ny,0:ns-1), 
     *  ne_ext(-1:nx,-1:ny), ne(-1:nx,-1:ny)
*     ------------------------------------------------------------------
*     B2XPNE computes the electron density, ne:
*       ne(,) = (sum is :: rza(,,is)*na(,,is)) .
*     A contribution from external charged species not followed by B2
*     can be included via use of the ne_ext array.
*     ------------------------------------------------------------------
      integer n2, is
      external sfill
*     ------------------------------------------------------------------
*$$$  call subini ('b2xpne')
      n2 = (nx+2)*(ny+2)
      call sfill (n2, 0.0_R8, ne, 1)
      do is = 0, ns-1
       ne(-1:nx,-1:ny) = ne(-1:nx,-1:ny) +
     .  rza(-1:nx,-1:ny,is)*na(-1:nx,-1:ny,is)
      enddo
      ne = ne + ne_ext
*$$$  call subend ()
      return
*     ------------------------------------------------------------------
      end subroutine b2xpne_old
