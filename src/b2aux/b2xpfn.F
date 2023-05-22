      subroutine b2xpfn (nFc, ns, fna, fnn)
      use b2mod_types
      use b2mod_constants
      use b2mod_b2cmpa
      implicit none
      integer nFc, ns
      real (kind=R8) ::
     *   fna(nFc,0:1,0:ns-1), fnn(nFc,0:1)
*     ------------------------------------------------------------------
*     B2XPFN computes the total particle flux of hydrogenic atoms, fnn:
*       fnn(,,) = (sum is, is_neutral(is).and.zn(is).eq.1 :: fna(,,is)).
*     ------------------------------------------------------------------
      integer is
      external sfill, b2saxpy
*     ------------------------------------------------------------------
*$$$  call subini ('b2xpfn')
      call sfill (nFc*2, 0.0_R8, fnn, 1)
      do is = 0, ns-1
        if (is_neutral(is).and.zn(is).eq.1.0_R8)
     .    call b2saxpy (nFc*2, 1.0_R8, fna(1,0,is), 1, fnn, 1)
      enddo
*$$$  call subend ()
      return
*     ------------------------------------------------------------------
      end subroutine b2xpfn
