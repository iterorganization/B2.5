      function damax (n, sx, incx)
      use b2mod_types
      use b2mod_math_diff
      implicit none
      real (kind=R8) ::
     *   damax
      integer n, incx
      real (kind=R8) ::
     *   sx(0:*)
*     ------------------------------------------------------------------
*     DAMAX (BLAS-like routine) computes the maximum absolute value in
*     a vector.
*     damax(n,sx,incx) returns the maximum absolute value in the
*     real*8 vector described by (n,sx,incx). BLAS storage conventions.
*     n must be positive.
*     ------------------------------------------------------------------
      integer i, ix
      real (kind=R8) ::
     *   t
      intrinsic max, abs
*     ------------------------------------------------------------------
      if (n.lt.1.or.incx.eq.0) then
       call xerrab ('damax--error, n.lt.1.or.incx.eq.0')
      else if (incx.eq.1) then
       t = abs(sx(0))
       do i = 1, n-1
        t = max(t,abs(sx(i)))
       enddo
       damax = t
      else
       ix = ibase(n,incx)
       t = abs(sx(ix))
       ix = ix+incx
       do i = 1, n-1
        t = max(t,abs(sx(ix)))
        ix = ix+incx
       enddo
       damax = t
      endif
      return
*     ------------------------------------------------------------------
      end function damax
