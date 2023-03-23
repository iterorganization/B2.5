      subroutine sbias (n, sa, sx, incx)
      use b2mod_types
      use b2mod_math
      implicit none
      integer n, incx
      real (kind=R8) ::
     *   sa, sx(0:(n-1)*incx)
*     ------------------------------------------------------------------
*     SBIAS (BLAS-like routine) adds a constant value to a vector.
*     sbias(n,sa,sx,incx) replaces the s.p. real*8 vector x by sa+x.
*     x is described by (n,sx,incx). BLAS storage conventions.
*     Return without effects if n.lt.1.
*     Undefined if incx.eq.0.
*     ------------------------------------------------------------------
      integer i, ix
*     ------------------------------------------------------------------
      if (incx.eq.0) then
       call xerrab ('sbias--error, incx.eq.0')
      else if (incx.eq.1) then
       do i = 0, n-1
        sx(i) = sa+sx(i)
       enddo
      else
       ix = ibase(n,incx)
       do i = 0, n-1
        sx(ix) = sa+sx(ix)
        ix = ix+incx
       enddo
      endif
      return
*     ------------------------------------------------------------------
      end subroutine sbias
