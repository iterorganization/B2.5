      function smax (n, sx, incx)
      use b2mod_types
      use b2mod_math_diffv
      implicit none
      real (kind=R8) ::
     *   smax
      integer n, incx
      real (kind=R8) ::
     *   sx(0:*)
*     ------------------------------------------------------------------
*     SMAX (BLAS-like routine) computes the maximum value in a vector.
*     smax(n,sx,incx) returns the maximum value in the s.p. real*8 vector
*     described by (n,sx,incx). BLAS storage conventions.
*     n must be positive.
*     ------------------------------------------------------------------
      integer i, ix
      real (kind=R8) ::
     *   t
      intrinsic max
      external xerrab
*     ------------------------------------------------------------------
      if (n.lt.1.or.incx.eq.0) then
       call xerrab ('smax--error, n.lt.1.or.incx.eq.0')
      else if (incx.eq.1) then
       t = sx(0)
       do i = 1, n-1
        t = max(t,sx(i))
       enddo
       smax = t
      else
       ix = ibase(n,incx)
       t = sx(ix)
       ix = ix+incx
       do i = 1, n-1
        t = max(t,sx(ix))
        ix = ix+incx
       enddo
       smax = t
      endif
      return
*     ------------------------------------------------------------------
      end function smax
