      subroutine sip5g (m, n, alf, a, b, c, d, e, r, x)
      use b2mod_types
      implicit none
      integer m, n
      real (kind=R8) ::
     *   alf, a(0:m*n-1), b(0:m*n-1), c(0:m*n-1), d(0:m*n-1),
     &  e(0:m*n-1), r(0:m*n-1), x(0:m*n-1)
*     ------------------------------------------------------------------
*     SIP5G computes the SIP-ILU factorisation of a 2d 5-point matrix,
*     and applies the factors for one iteration of an update equation.
*     ------------------------------------------------------------------
*     /1/ Herbert L. Stone: "Iterative Solution of Implicit Approxima-
*     tions of Multidimensional Partial Differential Equations",
*     SIAM J. Numer. Anal. 5 (1968) 530-558.
*     /2/ D. A. H. Jacobs: "SI205: A Subroutine for Iterating Using a
*     Five Diagonal Matrix to Solve a System of Algebraic Equations",
*     Central Electricity Generating Board (UK) Report RD/L/P8/79,
*     January 1980.
*     /3/ Numerical Algorithms Group: NAG Fortran Library Routine
*     D03UAF.
*     ------------------------------------------------------------------
      integer i
      external xertst, ilu5g
*     ------------------------------------------------------------------
      call xertst (1.le.m.and.1.le.n, 'sip5g--faulty argument m, n')
      if (alf.eq.0) then
*      (faster code for common special case)
       call ilu5g (m, n, a, b, c, d, e, r, x)
      else if (m.lt.3.or.n.lt.3) then
*      (standard ilu is an exact factorisation in this case)
       call ilu5g (m, n, a, b, c, d, e, r, x)
      else
       a(0) = a(0)
       b(0) = b(0)
       c(0) = c(0)
       d(0) = d(0)/c(0)
       e(0) = e(0)/c(0)
       x(0) = r(0)/c(0)
       do i = 1, m-1
        a(i) = a(i)
        b(i) = b(i)/(1+alf*e(i-1))
        c(i) = c(i)-b(i)*d(i-1)+alf*b(i)*e(i-1)
        d(i) = d(i)/c(i)
        e(i) = (e(i)-alf*b(i)*e(i-1))/c(i)
        x(i) = (r(i)-b(i)*x(i-1))/c(i)
       enddo
       do i = m, m*n-1
        a(i) = a(i)/(1+alf*d(i-m))
        b(i) = b(i)/(1+alf*e(i-1))
        c(i) = c(i)-a(i)*e(i-m)-b(i)*d(i-1)+
     &    alf*(a(i)*d(i-m)+b(i)*e(i-1))
        d(i) = (d(i)-alf*a(i)*d(i-m))/c(i)
        e(i) = (e(i)-alf*b(i)*e(i-1))/c(i)
        x(i) = (r(i)-b(i)*x(i-1)-a(i)*x(i-m))/c(i)
       enddo
       x(m*n-1) = x(m*n-1)
       do i = m*n-2, m*(n-1), -1
        x(i) = x(i)-d(i)*x(i+1)
       enddo
       do i = m*(n-1)-1, 0, -1
        x(i) = x(i)-d(i)*x(i+1)-e(i)*x(i+m)
       enddo
      endif
      return
*     ------------------------------------------------------------------
      end subroutine sip5g
