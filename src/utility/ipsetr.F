      subroutine ipsetr (cnam, rval)
      use b2mod_types
      implicit none
      character cnam*(*)
      real (kind=R8) ::
     *   rval
*     ------------------------------------------------------------------
*     IPSETR is an interface to ipsetc; it saves a real*8 value in the
*     table of internal parameters.
*     ------------------------------------------------------------------
      character cval0*32
      external ipsetc
*     ------------------------------------------------------------------
      write (cval0,'(1pe32.15)') rval
      call ipsetc (cnam, cval0)
      return
*     ------------------------------------------------------------------
      end
