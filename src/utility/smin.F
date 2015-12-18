      function smin (n, sx, incx)
      use b2mod_types
      implicit none
      real (kind=R8) ::
     *   smin
      integer n, incx
      real (kind=R8) ::
     *   sx(0:*)
*     ------------------------------------------------------------------
*     SMIN (BLAS-like routine) computes the minimum value in a vector.
*     smin(n,sx,incx) returns the minimum value in the s.p. real*8 vector
*     described by (n,sx,incx). BLAS storage conventions.
*     n must be positive.
*     ------------------------------------------------------------------
      integer i, ix
      real (kind=R8) ::
     *   t
      intrinsic min, max
      integer ibase
      ibase(n,incx) = max(0,(1-n)*incx)
*     ------------------------------------------------------------------
      if (n.lt.1.or.incx.eq.0) then
       call xerrab ('smin--error, n.lt.1.or.incx.eq.0')
      else if (incx.eq.1) then
       t = sx(0)
       do i = 1, n-1
        t = min(t,sx(i))
       enddo
       smin = t
      else
       ix = ibase(n,incx)
       t = sx(ix)
       ix = ix+incx
       do i = 1, n-1
        t = min(t,sx(ix))
        ix = ix+incx
       enddo
       smin = t
      endif
      return
*     ------------------------------------------------------------------
      end
