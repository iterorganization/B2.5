      subroutine sfill (n, sa, sx, incx)
      use b2mod_types
      use b2mod_openmp
      use b2mod_math
#ifdef WITH_FTIMINGS
      use timings
#endif
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
#if ( !defined(NO_OPENMP_SFILL) || defined(WITH_FTIMINGS) )
      logical make_parallel
#endif
      external xerrab
*     ------------------------------------------------------------------
#if ( !defined(NO_OPENMP_SFILL) || defined(WITH_FTIMINGS) )
      make_parallel = (.not. in_parallel())
#endif
#ifdef WITH_FTIMINGS
      call timer%start("sfill")
      if (make_parallel) call timer%start("parallel")
#endif
      if (incx.eq.0) then
       call xerrab ('sfill--error, incx.eq.0')
      else if (incx.eq.1) then
#ifndef NO_OPENMP_SFILL
!$OMP PARALLEL DO IF (make_parallel)
!$OMP& SCHEDULE(static)
!$OMP& DEFAULT(none)
!$OMP& PRIVATE(i)
!$OMP& SHARED(sx,sa,n)
#endif
       do i = 0, n-1
        sx(i) = sa
       enddo
      else
       ix = ibase(n,incx)
#ifndef NO_OPENMP_SFILL
!$OMP PARALLEL DO IF (make_parallel)
!$OMP& SCHEDULE(static)
!$OMP& DEFAULT(none)
!$OMP& PRIVATE(i)
!$OMP& SHARED(sx,sa,n,incx,ix)
#endif
       do i = 0, n-1
        sx(ix + i * incx) = sa
       enddo
      endif
#ifdef WITH_FTIMINGS
      if (make_parallel) call timer%stop("parallel")
      call timer%stop("sfill")
#endif
      return
*     ------------------------------------------------------------------
      end subroutine sfill
