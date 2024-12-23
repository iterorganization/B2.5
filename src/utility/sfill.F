! $AD NOCHECKPOINT
      subroutine sfill (n, sa, sx, incx)
      use b2mod_types
      use b2mod_math
      implicit none
      integer n, incx
      real (kind=R8) ::
     *   sa, sx(0:(n*incx)-1)
*     ------------------------------------------------------------------
*     SFILL (BLAS-like routine) fills a vector with a constant value.
*     sfill(n,sa,sx,incx) sets all elements of the s.p. real*8 vector x
*     to sa. x is described by (n,sx,incx). BLAS storage conventions.
*     Return without effects if n.lt.1.
*     Undefined if incx.eq.0.
*     ------------------------------------------------------------------
      integer i, ix
      external xerrab
*     ------------------------------------------------------------------
      if (incx.eq.0) then
       call xerrab ('sfill--error, incx.eq.0')
      else if (incx.eq.1) then
       do i = 0, n-1
        sx(i) = sa
       enddo
      else
       ix = ibase(n,incx)
       do i = 0, n-1
        sx(ix + i * incx) = sa
       enddo
      endif
      return
*     ------------------------------------------------------------------
      end subroutine sfill
