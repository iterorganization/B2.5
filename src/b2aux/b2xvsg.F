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
      end subroutine b2xvsg

      subroutine b2xvsg1(nx1, ny1, fun, inc, vrname, chsg)
      use b2mod_types
      use b2mod_geo2
      use b2mod_connectivity
      implicit none
      integer nx1, ny1, inc
      real (kind=R8) :: fun(-1:nx1,-1:ny1)
      character vrname*(*), chsg*4
*     -----------------------------------------------------------------
*     B2XVSG verifies the sign of some vector quantity over the B2 
*     coarse grid.
*     The vector fun is parsed to keep only the values over valid cells
*     and mapped to a one-dimensional list fun1d of size n.
*     The sign of each element of fun must be as specified
*     by chsg, which may be '.lt.', '.le.', '.gt.' or '.ge.'. If this
*     is violated the code will halt via a call to xertst/xerrab.
*     vrname must be the name of the program variable in the calling
*     routine that is associated with fun.
*     -----------------------------------------------------------------
      integer n
      integer i, ix, iy
      real (kind=R8), allocatable :: fun1d(:)
      external xertst
*     -----------------------------------------------------------------
*$$$  call subini ('b2xvsg1')
      call xertst (1.lt.nx1.and.1.lt.ny1.and.inc.ne.0,
     &  'b2xvsg1--faulty argument nx1, ny1, inc')
      n = (nx1+2)*(ny1+2)
      allocate(fun1d(0:n-1))
      i = 0
      do ix = -1, nx1
        do iy = -1, ny1
          if (isUnusedCell(cflags1(ix, iy, CELLFLAG_TYPE))) cycle
          fun1d(i) = fun(ix,iy)
          i = i+1
        end do
      end do
      call xertst(0.lt.i .and. i.le.n, 
     . 'b2xvsg1--Wrong number of used cells!')
      n = i-1
      call b2xvs1 (n, fun1d, inc, vrname, chsg) 
      deallocate(fun1d)
*$$$  call subend ()
      return
*     -----------------------------------------------------------------
      end subroutine b2xvsg1

      subroutine b2xvs1 (n, fun, inc, vrname, chsg)
      use b2mod_types
      implicit none
      integer n, inc
      real (kind=R8) ::
     *   fun(0:n-1)
      character vrname*(*), chsg*4
*     -----------------------------------------------------------------
*     B2XVS1 verifies the sign of some one-dimensional vector quantity.
*     The vector fun is represented by (n,fun,inc), using BLAS storage
*     conventions. The sign of each element of fun must be as specified
*     by chsg, which may be '.lt.', '.le.', '.gt.' or '.ge.'. If this
*     is violated the code will halt via a call to xertst/xerrab.
*     vrname must be the name of the program variable in the calling
*     routine that is associated with fun.
*     -----------------------------------------------------------------
      character*12 vrnam0
      real (kind=R8) ::
     *   smin, smax
      external xertst, xerrab, smin, smax
*     -----------------------------------------------------------------
*$$$  call subini ('b2xvs1')
      call xertst (0.lt.n.and.inc.ne.0,
     &  'b2xvs1--faulty argument n, inc')
      vrnam0 = trim(vrname)
      if (chsg.eq.'.lt.') then
       call xertst (smax(n,fun,inc).lt.0.0_R8,
     &   'b2xvs1--faulty argument sign '//vrnam0)
      else if (chsg.eq.'.le.') then
       call xertst (smax(n,fun,inc).le.0.0_R8,
     &   'b2xvs1--faulty argument sign '//vrnam0)
      else if (chsg.eq.'.gt.') then
       call xertst (smin(n,fun,inc).gt.0.0_R8,
     &   'b2xvs1--faulty argument sign '//vrnam0)
      else if (chsg.eq.'.ge.') then
       call xertst (smin(n,fun,inc).ge.0.0_R8,
     &   'b2xvs1--faulty argument sign '//vrnam0)
      else
       call xerrab ('b2xvs1--faulty argument chsg')
      endif
*$$$  call subend ()
      return
*     -----------------------------------------------------------------
      end subroutine b2xvs1
