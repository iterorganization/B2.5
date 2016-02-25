      subroutine b2xpfi (nx, ny, ns, fna, fni)
      use b2mod_types
      use b2mod_constants
      use b2mod_b2cmpa
      implicit none
      integer nx, ny, ns
      real (kind=R8) ::
     *   fna(-1:nx,-1:ny,0:1,0:1,0:ns-1), fni(-1:nx,-1:ny,0:1,0:1)
*     ------------------------------------------------------------------
*     B2XPFI computes the total atom particle flux, fni:
*       fni(,,) = (sum is :: fna(,,,is)) .
*     ------------------------------------------------------------------
      integer n2, is
      external sfill, b2saxpy
*     ------------------------------------------------------------------
*$$$  call subini ('b2xpfi')
      n2 = (nx+2)*(ny+2)
      call sfill (n2*4, 0.0_R8, fni, 1)
      do is = 0, ns-1
       call b2saxpy (n2*4, 1.0_R8, fna(-1,-1,0,0,is), 1, fni, 1)
      enddo
*$$$  call subend ()
      return
*     ------------------------------------------------------------------
      end
