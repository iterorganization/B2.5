      function samax (n, sx, incx)
      use b2mod_types
      implicit none
      real (kind=R8) ::
     *   samax
      integer n, incx
      real (kind=R8) ::
     *   sx(0:*)
*     ------------------------------------------------------------------
*     SAMAX (BLAS-like routine) computes the maximum absolute value in
*     a vector.
*     samax(n,sx,incx) returns the maximum absolute value in the s.p.
*     real*8 vector described by (n,sx,incx). BLAS storage conventions.
*     n must be positive.
*     ------------------------------------------------------------------
      integer i, ix
      real (kind=R8) ::
     *   t
      intrinsic max, abs
      integer ibase
      ibase(n,incx) = max(0,(1-n)*incx)
*     ------------------------------------------------------------------
      if (n.lt.1.or.incx.eq.0) then
       call xerrab ('samax--error, n.lt.1.or.incx.eq.0')
      else if (incx.eq.1) then
       t = abs(sx(0))
       do i = 1, n-1
        t = max(t,abs(sx(i)))
       enddo
       samax = t
      else
       ix = ibase(n,incx)
       t = abs(sx(ix))
       ix = ix+incx
       do i = 1, n-1
        t = max(t,abs(sx(ix)))
        ix = ix+incx
       enddo
       samax = t
      endif
      return
*     ------------------------------------------------------------------
      end
