      subroutine ipgetr (cnam, rval)
      use b2mod_types
      implicit none
      character cnam*(*)
      real (kind=R8) ::
     *   rval
*     ------------------------------------------------------------------
*     IPGETR is an interface to ipgetc; it retrieves a real*8 value
*     rval from the table of internal parameters, or returns rval
*     unchanged if the name is not found or if the corresponding entry
*     is a blank string. The internal parameter may have been set by an
*     earlier call to ipsetr or to ipsetc.
*     ------------------------------------------------------------------
      character cnam0*32, cval0*32
      external xerrab, ipgetc
*     ------------------------------------------------------------------
      cnam0 = cnam
      cval0 = ' '
      call ipgetc (cnam, cval0)
      if (cval0.ne.' ') then
       read (cval0,'(e32.0)',err=91) rval
      endif
      return
*     ------------------------------------------------------------------
   91 call xerrab ('ipgetr--invalid stored value for: '//cnam0)
*     ------------------------------------------------------------------
      end subroutine ipgetr
