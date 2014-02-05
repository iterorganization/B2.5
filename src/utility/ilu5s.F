      subroutine ilu5s (m, n, a, b, c, d, e, r, x)
      use b2mod_types
      implicit none
      integer m, n
      real (kind=R8) ::
     *   a(0:m*n-1), b(0:m*n-1), c(0:m*n-1), d(0:m*n-1),
     &  e(0:m*n-1), r(0:m*n-1), x(0:m*n-1)
*     ------------------------------------------------------------------
*     ILU5S solves an equation of the form L*U*x=r, where L and U are
*     the factors in an ILU factorisation of a 2d-5point matrix.
*     ------------------------------------------------------------------
      integer i
      external xertst
*     ------------------------------------------------------------------
      call xertst (1.le.m.and.1.le.n, 'ilu5s--faulty argument m, n')
      x(0) = r(0)/c(0)
      do i = 1, m-1
       x(i) = (r(i)-b(i)*x(i-1))/c(i)
      enddo
      do i = m, m*n-1
       x(i) = (r(i)-b(i)*x(i-1)-a(i)*x(i-m))/c(i)
      enddo
      x(m*n-1) = x(m*n-1)
      do i = m*n-2, m*(n-1), -1
       x(i) = x(i)-d(i)*x(i+1)
      enddo
      do i = m*(n-1)-1, 0, -1
       x(i) = x(i)-d(i)*x(i+1)-e(i)*x(i+m)
      enddo
      return
*     ------------------------------------------------------------------
      end
