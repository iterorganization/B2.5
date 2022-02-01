      subroutine ilu5g (m, n, a, b, c, d, e, r, x)
      use b2mod_types
      implicit none
      integer m, n
      real (kind=R8) ::
     *   a(0:m*n-1), b(0:m*n-1), c(0:m*n-1), d(0:m*n-1),
     &  e(0:m*n-1), r(0:m*n-1), x(0:m*n-1)
*     ------------------------------------------------------------------
*     ILU5G computes the ILU factorisation of a 2d 5-point matrix and
*     applies the factors for one iteration of an update equation.
*     ------------------------------------------------------------------
      integer i
      external xertst
*     ------------------------------------------------------------------
      call xertst (1.le.m.and.1.le.n, 'ilu5g--faulty argument m, n')
*     a(0) = a(0)
*     b(0) = b(0)
*     c(0) = c(0)
      d(0) = d(0)/c(0)
      e(0) = e(0)/c(0)
      x(0) = r(0)/c(0)
      do i = 1, m-1
*      a(i) = a(i)
*      b(i) = b(i)
       c(i) = c(i)-b(i)*d(i-1)
       d(i) = d(i)/c(i)
       e(i) = e(i)/c(i)
       x(i) = (r(i)-b(i)*x(i-1))/c(i)
      enddo
      do i = m, m*n-1
*      a(i) = a(i)
*      b(i) = b(i)
       c(i) = c(i)-a(i)*e(i-m)-b(i)*d(i-1)
       d(i) = d(i)/c(i)
       e(i) = e(i)/c(i)
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
      end subroutine ilu5g
