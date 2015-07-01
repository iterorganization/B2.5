      subroutine b2xvsg (n, fun, inc, vrname, chsg)
      use b2mod_types
      implicit none
      integer n, inc
      real (kind=R8) ::
     *   fun(0:*)
      character vrname*(*), chsg*(*)
*     -----------------------------------------------------------------
*     B2XVSG verifies the sign of some vector quantity.
*     The vector fun is represented by (n,fun,inc), using BLAS storage
*     conventions. The sign of each element of fun must be as specified
*     by chsg, which may be '.lt.', '.le.', '.gt.' or '.ge.'. If this
*     is violated the code will halt via a call to xertst/xerrab.
*     vrname must be the name of the program variable in the calling
*     routine that is associated with fun.
*     -----------------------------------------------------------------
!      integer, save :: ncall = 0, no_test_sign = 1                       !srv 26.10.08  11.01.13
      character vrnam0*8
      real (kind=R8) ::
     *   smin, smax
      external xertst, xerrab, smin, smax
*     -----------------------------------------------------------------
*$$$  call subini ('b2xvsg')
*   ..set internal parameters on first call
!      if (ncall.eq.0) then
!       call ipgeti ('b2xvsg_no_test_sign', no_test_sign)                 !srv 26.10.08
!	 write(*,*) ' ncall =', ncall, ' no_test_sign =',no_test_sign
!      endif
      call xertst (0.lt.n.and.inc.ne.0,
     &  'b2xvsg--faulty argument n, inc')
      vrnam0 = vrname
!	if(no_test_sign.eq.0) then                                         !srv 26.10.08
      if (chsg.eq.'.lt.') then
        call xertst (smax(n,fun,inc).lt.0,
     &   'b2xvsg--faulty argument sign '//vrnam0)
      else if (chsg.eq.'.le.') then
        call xertst (smax(n,fun,inc).le.0,
     &   'b2xvsg--faulty argument sign '//vrnam0)
      else if (chsg.eq.'.gt.') then
        call xertst (smin(n,fun,inc).gt.0,
     &   'b2xvsg--faulty argument sign '//vrnam0)
      else if (chsg.eq.'.ge.') then
        call xertst (smin(n,fun,inc).ge.0,
     &   'b2xvsg--faulty argument sign '//vrnam0)
      else
        call xerrab ('b2xvsg--faulty argument chsg')
      endif
!	endif
!	ncall = ncall + 1
*$$$  call subend ()
      return
*     -----------------------------------------------------------------
      end
