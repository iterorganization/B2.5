      subroutine ipseti (cnam, ival)
      use b2mod_types
      implicit none
      character cnam*(*)
      integer ival
*     ------------------------------------------------------------------
*     IPSETI is an interface to ipsetc; it saves an integer value in
*     the table of internal parameters.
*     ------------------------------------------------------------------
      character cval0*32
      external ipsetc
*     ------------------------------------------------------------------
      write (cval0,'(i32)') ival
      call ipsetc (cnam, cval0)
      return
*     ------------------------------------------------------------------
      end subroutine ipseti
