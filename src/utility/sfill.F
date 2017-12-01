      subroutine sfill (n, sa, sx, incx)
      use b2mod_types
      use b2mod_openmp
#ifdef WITH_FTIMINGS
      use timings
#endif
      implicit none
      integer n, incx
      real (kind=R8) ::
     *   sa, sx(0:*)
*     ------------------------------------------------------------------
*     SFILL (BLAS-like routine) fills a vector with a constant value.
*     sfill(n,sa,sx,incx) sets all elements of the s.p. real*8 vector x
*     to sa. x is described by (n,sx,incx). BLAS storage conventions.
*     Return without effects if n.lt.1.
*     Undefined if incx.eq.0.
*     ------------------------------------------------------------------
      integer i, ix
      intrinsic max
      logical do_parallel
      integer ibase
      ibase(n,incx) = max(0,(1-n)*incx)
*     ------------------------------------------------------------------

#ifndef NO_OPENMP_SFILL
      do_parallel =
     &  (.not. in_parallel())
     &  .and. (n > 16384)
#endif

#ifdef WITH_FTIMINGS
      call timer%start("sfill")
      if (do_parallel) call timer%start("parallel")
#endif
      if (incx.eq.0) then
       call xerrab ('sfill--error, incx.eq.0')
      else if (incx.eq.1) then
#ifndef NO_OPENMP_SFILL
!$OMP PARALLEL DO IF (do_parallel)
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
!$OMP PARALLEL DO IF (do_parallel)
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
      if (do_parallel) call timer%stop("parallel")
      call timer%stop("sfill")
#endif
      return
*     ------------------------------------------------------------------
      end
