      subroutine ilu5f (m, n, a, b, c, d, e)
      use b2mod_types
      implicit none
      integer m, n
      real (kind=R8) ::
     *   a(0:m*n-1), b(0:m*n-1), c(0:m*n-1), d(0:m*n-1),
     &  e(0:m*n-1)
*     ------------------------------------------------------------------
*     ILU5F computes the ILU factorisation of a 2d-5point matrix.
*     ------------------------------------------------------------------
      integer i
      external xertst
*     ------------------------------------------------------------------
      call xertst (1.le.m.and.1.le.n, 'ilu5f--faulty argument m, n')
*     a(0) = a(0)
*     b(0) = b(0)
*     c(0) = c(0)
      d(0) = d(0)/c(0)
      e(0) = e(0)/c(0)
      do i = 1, m-1
*      a(i) = a(i)
*      b(i) = b(i)
       c(i) = c(i)-b(i)*d(i-1)
       d(i) = d(i)/c(i)
       e(i) = e(i)/c(i)
      enddo
      do i = m, m*n-1
*      a(i) = a(i)
*      b(i) = b(i)
       c(i) = c(i)-a(i)*e(i-m)-b(i)*d(i-1)
       d(i) = d(i)/c(i)
       e(i) = e(i)/c(i)
      enddo
      return
*     ------------------------------------------------------------------
      end
