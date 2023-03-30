      subroutine b2xpfn (nFc, ns, ismain0, fna, fnn)
      use b2mod_types
      use b2mod_constants
      use b2mod_b2cmpa
      implicit none
      integer nFc, ns, ismain0
      real (kind=R8) ::
     *   fna(nFc,0:1,0:ns-1), fnn(nFc,0:1)
*     ------------------------------------------------------------------
*     B2XPFN computes the total particle flux of main hydrogenic atoms, fnn
*     ------------------------------------------------------------------
      integer is
      external sfill, b2saxpy
*     ------------------------------------------------------------------
*$$$  call subini ('b2xpfn')
      call sfill (nFc*2, 0.0_R8, fnn, 1)
      if (is_neutral(ismain0).and.nint(zn(ismain0)).eq.1)
     > call b2saxpy (nFc*2, 1.0_R8, fna(1,0,ismain0), 1, fnn, 1)
*$$$  call subend ()
      return
*     ------------------------------------------------------------------
      end subroutine b2xpfn
