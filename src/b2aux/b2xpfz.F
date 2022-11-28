      subroutine b2xpfz (nFc, ns, fna, fni)
      use b2mod_types
      use b2mod_constants
      use b2mod_b2cmpa
      implicit none
      integer nFc, ns
      real (kind=R8) ::
     *   fna(nFc,0:1,0:ns-1), fni(nFc,0:1)
*     ------------------------------------------------------------------
*     B2XPFZ computes the total ion particle flux, fni:
*       fni(,) = (sum is, .not.is_neutral(is) :: fna(,,is)) .
*     ------------------------------------------------------------------
      integer is
      external sfill, b2saxpy
*     ------------------------------------------------------------------
*$$$  call subini ('b2xpfz')
      call sfill (nFc*2, 0.0_R8, fni, 1)
      do is = 0, ns-1
        if (.not.is_neutral(is))
     .    call b2saxpy (nFc*2, 1.0_R8, fna(1,0,is), 1, fni, 1)
      enddo
*$$$  call subend ()
      return
*     ------------------------------------------------------------------
      end subroutine b2xpfz
