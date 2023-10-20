      function b2ssum (n, sx, incx)
      use b2mod_types
      use b2mod_math
      implicit none
      real (kind=R8) ::
     *   b2ssum
      integer n, incx
      real (kind=R8) ::
     *   sx(0:n*incx-1)
*     ------------------------------------------------------------------
*     B2SSUM (BLAS-like routine) sums the elements of a vector.
*     b2ssum(n,sx,incx) returns the sum of the elements in the d.p. real
*     vector described by (n,sx,incx). BLAS storage conventions.
*     Return with b2ssum.eq.0 if n.lt.1.
*     Undefined if incx.eq.0.
*     ------------------------------------------------------------------
      integer i, ix
      real (kind=R8) ::
     *   t
      external xerrab
*     ------------------------------------------------------------------
      if (incx.eq.0) then
       call xerrab ('b2ssum--error, incx.eq.0')
      else if (incx.eq.1) then
       t = 0.0_R8
       do i = 0, n-1
        t = t+sx(i)
       enddo
       b2ssum = t
      else
       ix = ibase(n,incx)
       t = 0.0_R8
       do i = 0, n-1
        t = t+sx(ix)
        ix = ix+incx
       enddo
       b2ssum = t
      endif
      return
*     ------------------------------------------------------------------
      end function b2ssum
      function b2sasum (n, sx, incx)
      use b2mod_types
      use b2mod_math
      implicit none
      real (kind=R8) ::
     *   b2sasum
      integer n, incx
      real (kind=R8) ::
     *   sx(0:*)
*     ------------------------------------------------------------------
*     B2SASUM (BLAS routine) computes the L1-norm of a vector.
*     b2sasum(n,sx,incx) returns the sum of absolute values of the
*     elements in the d.p. real vector described by (n,sx,incx).
*     BLAS storage conventions. Return with b2sasum.eq.0 if n.lt.1.
*     Undefined if incx.eq.0.
*     ------------------------------------------------------------------
      integer i, ix
      real (kind=R8) ::
     *   t
      intrinsic abs
      external xerrab
*     ------------------------------------------------------------------
      if (incx.eq.0) then
       call xerrab ('b2sasum--error, incx.eq.0')
      else if (incx.eq.1) then
       t = 0.0_R8
       do i = 0, n-1
        t = t+abs(sx(i))
       enddo
       b2sasum = t
      else
       ix = ibase(n,incx)
       t = 0.0_R8
       do i = 0, n-1
        t = t+abs(sx(ix))
        ix = ix+incx
       enddo
       b2sasum = t
      endif
      return
*     ------------------------------------------------------------------
      end function b2sasum
      function b2s2sum (n1, n2, sx, inc1x, inc2x)
      use b2mod_types
      use b2mod_math
      implicit none
      real (kind=R8) ::
     *   b2s2sum
      integer n1, n2, inc1x, inc2x
      real (kind=R8) ::
     *   sx(*)
*     ------------------------------------------------------------------
*     COMPUTE SUM OF VALUES IN THE D.P. REAL MATRIX DESCRIBED BY
*     (N1, N2, SX, INC1X, INC2X). BLAS**2 STORAGE CONVENTIONS.
*     RETURN WITH B2S2SUM = 0 WHEN N1 .LT. 1 OR N2 .LT. 1.
*     ------------------------------------------------------------------
      integer i, i2
      real (kind=R8) ::
     *   t, b2ssum
      external b2ssum
*     ------------------------------------------------------------------
      if (n1 .lt. 1 .or. n2 .lt. 1) then
         b2s2sum = 0.0_R8
         return
*        - - - - - - -
      else if (n1*inc1x .eq. inc2x) then
         b2s2sum = b2ssum (n1*n2, sx, inc1x)
         return
*        - - - - - - -
      else
         t = 0.0_R8
         i2 = ibase1(n2, inc2x)
         do 1, i = 1, n2
            t = t + b2ssum (n1, sx(i2), inc1x)
            i2 = i2 + inc2x
    1    continue
         b2s2sum = t
         return
*        - - - - - - -
      end if
*     ---------------------------------------------------BJB-83/04/06---
      end function b2s2sum
      subroutine b2sscal (n, sa, sx, incx)
      use b2mod_types
      use b2mod_math
      implicit none
      integer n, incx
      real (kind=R8) ::
     *   sa, sx(0:*)
*     ------------------------------------------------------------------
*     B2SSCAL (BLAS routine) performs the vector operation x:=sa*x.
*     x is an d.p. real vector, described by (n,sx,incx). BLAS storage
*     conventions. Return without effects if sa.eq.1.
*     Undefined if incx.eq.0.
*     ------------------------------------------------------------------
      integer i, ix
      external xerrab
*     ------------------------------------------------------------------
      if (incx.eq.0) then
       call xerrab ('b2sscal--error, incx.eq.0')
      else if (sa.eq.1.0_R8) then
      else if (incx.eq.1) then
       do i = 0, n-1
        sx(i) = sa*sx(i)
       enddo
      else
       ix = ibase(n,incx)
       do i = 0, n-1
        sx(ix) = sa*sx(ix)
        ix = ix+incx
       enddo
      endif
      return
*     ------------------------------------------------------------------
      end subroutine b2sscal
      subroutine b2scopy (n, sx, incx, sy, incy)
      use b2mod_types
      use b2mod_math
      implicit none
      integer n, incx, incy
      real (kind=R8) ::
     *   sx(0:*), sy(0:*)
*     ------------------------------------------------------------------
*     B2SCOPY (BLAS routine) performs the vector operation y:=x.
*     x and y are d.p. real vectors described by (n,sx,incx) and
*     (n,sy,incy). BLAS storage conventions.
*     Undefined if incy.eq.0.
*     ------------------------------------------------------------------
      integer i, ix, iy
      external xerrab
*     ------------------------------------------------------------------
      if (incy.eq.0) then
       call xerrab ('b2scopy--error, incy.eq.0')
      else if (incx.eq.1.and.incy.eq.1) then
       do i = 0, n-1
        sy(i) = sx(i)
       enddo
      else
       ix = ibase(n,incx)
       iy = ibase(n,incy)
       do i = 0, n-1
        sy(iy) = sx(ix)
        ix = ix+incx
        iy = iy+incy
       enddo
      endif
      return
*     ------------------------------------------------------------------
      end subroutine b2scopy
      subroutine b2saxpy (n, sa, sx, incx, sy, incy)
      use b2mod_types
      use b2mod_math
      use b2mod_openmp
#ifdef WITH_FTIMINGS
      use timings
#endif
      implicit none
      integer n, incx, incy
      real (kind=R8) ::
     *   sa, sx(0:*), sy(0:*)
#if ( !defined(NO_OPENMP_MYBLAS) || defined(WITH_FTIMINGS) )
      logical do_parallel
#endif
*     ------------------------------------------------------------------
*     B2SAXPY (BLAS routine) performs the vector operation y:=y+sa*x.
*     x and y are d.p. real vectors, described by (n,sx,incx) and
*     (n,sy,incy). BLAS storage conventions.
*     Return without effects if sa.eq.0.
*     Undefined if incy.eq.0.
*     ------------------------------------------------------------------
      integer i, ix, iy
      external xerrab
*     ------------------------------------------------------------------
#if ( !defined(NO_OPENMP_MYBLAS) || defined(WITH_FTIMINGS) )
      do_parallel =  (.not. in_parallel()) .and. (n > 7000)
#endif
#ifdef WITH_FTIMINGS
      call timer%start("b2saxpy")
      if (do_parallel) call timer%start("parallel")
#endif
      if (incy.eq.0) then
       call xerrab ('b2saxpy--error, incy.eq.0')
      else if (sa.eq.0.0_R8) then
      else if (incx.eq.1.and.incy.eq.1) then
#ifndef NO_OPENMP_MYBLAS
!$OMP PARALLEL DO IF (do_parallel)
!$OMP& SCHEDULE(static)
!$OMP& DEFAULT(none)
!$OMP& PRIVATE(i)
!$OMP& SHARED(sx,sy,sa,n)
#endif
       do i = 0, n-1
        sy(i) = sy(i)+sa*sx(i)
       enddo
      else
       ix = ibase(n,incx)
       iy = ibase(n,incy)
       do i = 0, n-1
        sy(iy) = sy(iy)+sa*sx(ix)
        ix = ix+incx
        iy = iy+incy
       enddo
      endif
#ifdef WITH_FTIMINGS
      if(do_parallel) call timer%stop("parallel")
      call timer%stop("b2saxpy")
#endif
      return
*     ------------------------------------------------------------------
      end subroutine b2saxpy
      function b2s2min (n1, n2, sx, inc1x, inc2x)
      use b2mod_types
      use b2mod_math
      implicit none
      real (kind=R8) ::
     *   b2s2min
      integer n1, n2, inc1x, inc2x
      real (kind=R8) ::
     *   sx(*)
*     ------------------------------------------------------------------
*     COMPUTE MINIMUM VALUE IN THE D.P. REAL MATRIX DESCRIBED BY
*     (N1, N2, SX, INC1X, INC2X). BLAS**2 STORAGE CONVENTIONS.
*     N1 AND N2 MUST BE POSITIVE.
*     ------------------------------------------------------------------
      integer i, i2
      real (kind=R8) ::
     *   t, smin
      intrinsic min
      external smin, xerrab
*     ------------------------------------------------------------------
      if (n1 .lt. 1 .or. n2 .lt. 1) then
         call xerrab ('b2s2min: invalid n1/n2 values')
*        - - - - - - -
      else if (n1*inc1x .eq. inc2x) then
         b2s2min = smin (n1*n2, sx, inc1x)
         return
*        - - - - - - -
      else
         i2 = ibase1(n2, inc2x)
         t = smin (n1, sx(i2), inc1x)
         i2 = i2 + inc2x
         do 1, i = 2, n2
            t = min (t, smin (n1, sx(i2), inc1x))
            i2 = i2 + inc2x
    1    continue
         b2s2min = t
         return
*        - - - - - - -
      end if
*     ---------------------------------------------------BJB-83/04/06---
      end function b2s2min
      function b2sign (a, b)
      use b2mod_types
      implicit none
      real (kind=R8) ::
     *   b2sign
      real (kind=R8) ::
     *   a, b
*     ------------------------------------------------------------------
*     MODIFIED SIGN FUNCTION THAT RETURNS ZERO IF B is ZERO.
*     ------------------------------------------------------------------
      intrinsic sign
*     ------------------------------------------------------------------
      if (b.eq.0.0_R8) then
        b2sign = 0.0_R8
      else
        b2sign = sign(a,b)
      end if
*     -------------------------------------------------XPB-2013/09/16---
      end function b2sign
