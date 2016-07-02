      function ssum (n, sx, incx)
      use b2mod_types
      implicit none
      real (kind=R8) ::
     *   ssum
      integer n, incx
      real (kind=R8) ::
     *   sx(0:n-1)
*     ------------------------------------------------------------------
*     SSUM (BLAS-like routine) sums the elements of a vector.
*     ssum(n,sx,incx) returns the sum of the elements in the s.p. real
*     vector described by (n,sx,incx). BLAS storage conventions.
*     Return with ssum.eq.0 if n.lt.1.
*     Undefined if incx.eq.0.
*     ------------------------------------------------------------------
      integer i, ix
      real (kind=R8) ::
     *   t
      intrinsic max
      external xerrab
      integer ibase
      ibase(n,incx) = max(0,(1-n)*incx)
*     ------------------------------------------------------------------
      if (incx.eq.0) then
       call xerrab ('ssum--error, incx.eq.0')
      else if (incx.eq.1) then
       t = 0.0_R8
       do i = 0, n-1
        t = t+sx(i)
       enddo
       ssum = t
      else
       ix = ibase(n,incx)
       t = 0.0_R8
       do i = 0, n-1
        t = t+sx(ix)
        ix = ix+incx
       enddo
       ssum = t
      endif
      return
*     ------------------------------------------------------------------
      end
