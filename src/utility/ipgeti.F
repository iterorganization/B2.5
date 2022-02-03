      subroutine ipgeti (cnam, ival)
      use b2mod_types
      use b2mod_ipmain
      implicit none
      character cnam*(*)
      integer ival
*     ------------------------------------------------------------------
*     IPGETI is an interface to ipgetc; it retrieves an integer value
*     ival from the table of internal parameters, or returns ival
*     unchanged if the name is not found or if the corresponding entry
*     is a blank string. The internal parameter may have been set by an
*     earlier call to ipseti or to ipsetc.
*     ------------------------------------------------------------------
      character cnam0*32, cval0*32
      external xerrab
*     ------------------------------------------------------------------
      cnam0 = cnam
      cval0 = ' '
      call ipgetc (cnam, cval0)
      if (cval0.ne.' ') then
       read (cval0,'(i32)',err=91) ival
      endif
      return
*     ------------------------------------------------------------------
   91 call xerrab ('ipgeti--invalid stored value for: '//cnam0)
*     ------------------------------------------------------------------
      end subroutine ipgeti
