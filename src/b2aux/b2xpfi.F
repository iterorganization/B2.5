      subroutine b2xpfi (nFc, ns, fna, fni)
      use b2mod_types
      use b2mod_constants
      use b2mod_b2cmpa
      implicit none
      integer nFc, ns
      real (kind=R8) ::
     *   fna(nFc,0:1,0:ns-1), fni(nFc,0:1)
*     ------------------------------------------------------------------
*     B2XPFI computes the total atom particle flux, fni:
*       fni(,) = (sum is :: fna(,,is)) .
*     ------------------------------------------------------------------
      integer is
      external sfill, b2saxpy
*     ------------------------------------------------------------------
*$$$  call subini ('b2xpfi')
      call sfill (nFc*2, 0.0_R8, fni, 1)
      do is = 0, ns-1
       call b2saxpy (nFc*2, 1.0_R8, fna(1,0,is), 1, fni, 1)
      enddo
*$$$  call subend ()
      return
*     ------------------------------------------------------------------
      end subroutine b2xpfi
