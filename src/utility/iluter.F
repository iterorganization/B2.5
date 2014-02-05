*-----------------------------------------------------------------------
*.specification

      subroutine iluter (m, n, aa5, rhs, iters, acpar, itpar,
     &  sol, nrmc, bb5, wrk0, wrk1, wrk2)
      use b2mod_types
      implicit none
*   ..input arguments (unchanged on exit)
      integer m, n, iters, itpar
      real (kind=R8) ::
     *   acpar, aa5(0:m-1,0:n-1,-2:2), rhs(0:m-1,0:n-1)
*   ..output arguments (unspecified on entry)
      real (kind=R8) ::
     *   sol(0:m-1,0:n-1), nrmc(0:iters-1)
*   ..workspace arguments (unspecified on entry and on exit)
      real (kind=R8) ::
     *   bb5(0:m-1,0:n-1,-2:2), wrk0(0:m-1,0:n-1),
     &  wrk1(0:m-1,0:n-1), wrk2(0:m-1,0:n-1)

*-----------------------------------------------------------------------
*.documentation
*
*  1. purpose
*
*     ILUTER performs a number of iterations of the Stone strongly
*     implicit procedure for a system of equations of five-point form.
*
*
*  3. description
*
*     This routine employs the Strongly Implicit Procedure (SIP-ILU)
*     to obtain an approximate solution to a five-point system of
*     equations. The system is
*       aa5(i,j,-2)*sol(i,j-1)+aa5(i,j,-1)*sol(i-1,j)+
*         aa5(i,j,0)*sol(i,j)+aa5(i,j,1)*sol(i+1,j)+
*         aa5(i,j,2)*sol(i,j+1) = rhs(i,j)
*     for (i,j) in (0:m-1,0:n-1); any term involving indices outside
*     (0:m-1,0:n-1) is to be understood as 0.
*
*     The initial guess for the solution is identically zero, so when
*     this routine is embedded in another iterative procedure, rhs will
*     represent a residual, and sol will represent a correction.
*
*
*  4. references
*
*     /1/ Herbert L. Stone: "Iterative Solution of Implicit Approxima-
*     tions of Multidimensional Partial Differential Equations",
*     SIAM J. Numer. Anal. 5 (1968) 530-558.
*     /2/ D. A. H. Jacobs: "SI205: A Subroutine for Iterating Using a
*     Five Diagonal Matrix to Solve a System of Algebraic Equations",
*     Central Electricity Generating Board (UK) Report RD/L/P8/79,
*     January 1980.
*     /3/ Numerical Algorithms Group: NAG Fortran Library Routine
*     D03UAF.
*
*
*  5. parameters
*
*     m, n - integer, input.
*     m and n specify the number of mesh points along the first and the
*     second coordinate direction, respectively. The mesh points are
*     indexed by (0:m-1,0:n-1). It will hold that 2.le.m and 2.le.n.
*
*     aa5 - (0:m-1,0:n-1,-2:2) real*8 array, input.
*     aa5 specifies the coefficients of the five-point matrix as
*     described in Section 3. Any coefficient that would multiply an
*     entry outside the mesh must have the value 0.
*
*     rhs - (0:m-1,0:n-1) real*8 array, input.
*     rhs specifies the right hand side of the equation.
*
*     iters - integer, input.
*     iters specifies the number of iterations.
*
*     acpar - real, input.
*     acpar specifies the iteration parameter as suggested by Jacobs
*     and employed in the NAG routine D03UAF.
*
*     itpar - integer, input.
*     itpar specifies the initial value of the iteration parameter.
*
*     sol - (0:m-1,0:n-1) real*8 array, output.
*     sol specifies the calculated solution at the end of the
*     iterations.
*
*     nrmc - (0:iters-1) real*8 array, output.
*     For 0.le.it.lt.iters, nrmc(it) specifies the L1-norm of the
*     correction at iteration number (it).
*
*     bb5 - (0:m-1,0:n-1,-2:2) real*8 array, workspace.
*     bb5 will be used to store a copy of aa5 (perhaps with a different
*     ordering of the entries) that will be destroyed in the ILU
*     factorisation process.
*
*     wrk0, wrk1, wrk2 - (0:m-1,0:n-1) real*8 array, workspace.
*
*
*  8. further comments
*
*     Consider to use one-dimensional indexing for the mesh; i.e.,
*     0:m-1,0:n-1 becomes 0:m*n-1.
*
*
*-----------------------------------------------------------------------
*.declarations

*   ..local scalars
      integer iter, itp, i, j
      real (kind=R8) ::
     *   rxf
      parameter (rxf=1.0_R8)
      real (kind=R8) ::
     *   alp(0:8), alf
      save alp
*   ..procedures
      intrinsic max, mod
      real (kind=R8) ::
     *   smin, samax
      external xerrab, smin, samax, sfill, scopy, saxpy, sip5g
      data alp/0.0_R8,0.625_R8,0.250_R8,
     1 0.875_R8,0.500_R8,0.125_R8,0.750_R8,
     1 0.375_R8,1.0_R8/

*-----------------------------------------------------------------------
*.computation

* ..preliminaries
*   ..test m, n, iters, acpar, itpar
      if (m.lt.2.or.n.lt.2) then
       call xerrab ('iluter--faulty argument m, n')
      else if (iters.lt.0) then
       call xerrab ('iluter--faulty argument iters')
      else if (acpar.le.0) then
       call xerrab ('iluter--faulty argument acpar')
      else if (itpar.lt.0) then
       call xerrab ('iluter--faulty argument itpar')
      endif
*   ..test edge values of side-diagonals
      if (0.lt.samax(m,aa5(0,0,-2),1).or.
     &  0.lt.samax(m,aa5(0,n-1,2),1)) then
       call xerrab ('iluter--faulty argument aa5(,,-2), aa5(,,2)')
      else if (0.lt.samax(n,aa5(0,0,-1),m).or.
     &   0.lt.samax(n,aa5(m-1,0,1),m)) then
       call xerrab ('iluter--faulty argument aa5(,,-1), aa5(,,1)')
      endif
*   ..test sign of main diagonal
      if (smin(m*n,aa5(0,0,0),1).lt.0) then
       call xerrab ('iluter--faulty argument sign aa5(,,0)')
      endif

* ..perform the iterations
*   ..initialise sol
      call sfill (m*n, 0.0_R8, sol, 1)
*   ..iterate
      do iter = 0, iters-1
*    ..compute the residual
       call scopy (m*n, rhs, 1, wrk0, 1)
       if (0.lt.iter) then
        do j = 0, n-1
         wrk0(0,j) = wrk0(0,j)-aa5(0,j,0)*sol(0,j)
         do i = 1, m-1
          wrk0(i,j) = wrk0(i,j)-aa5(i,j,0)*sol(i,j)-
     &      aa5(i,j,-1)*sol(i-1,j)
          wrk0(i-1,j) = wrk0(i-1,j)-aa5(i-1,j,1)*sol(i,j)
         enddo
        enddo
        do j = 1, n-1
         do i = 0, m-1
          wrk0(i,j) = wrk0(i,j)-aa5(i,j,-2)*sol(i,j-1)
          wrk0(i,j-1) = wrk0(i,j-1)-aa5(i,j-1,2)*sol(i,j)
         enddo
        enddo
       endif
*    ..compute iteration matrix and rhs
       if (mod(itpar+iter,2).eq.0) then
        call scopy (m*n, wrk0, 1, wrk1, 1)
        call scopy (5*m*n, aa5, 1, bb5, 1)
       else
        do j = 0, n-1
         call scopy (m, wrk0(0,j), 1, wrk1(0,n-1-j), 1)
         call scopy (m, aa5(0,j,-2), 1, bb5(0,n-1-j,2), 1)
         call scopy (m, aa5(0,j,-1), 1, bb5(0,n-1-j,-1), 1)
         call scopy (m, aa5(0,j,0), 1, bb5(0,n-1-j,0), 1)
         call scopy (m, aa5(0,j,1), 1, bb5(0,n-1-j,1), 1)
         call scopy (m, aa5(0,j,2), 1, bb5(0,n-1-j,-2), 1)
        enddo
       endif
       call stabeq (m, n, 1.0_R8, bb5)
*    ..compute sip-ilu solution
       itp = mod(itpar+iter,18)/2
       alf = max(0.0_R8,1.0_R8-(2*acpar/(m**2+n**2))**alp(itp))
       call sip5g (m, n, alf, bb5(0,0,-2), bb5(0,0,-1), bb5(0,0,0),
     &   bb5(0,0,1), bb5(0,0,2), wrk1, wrk2)
*    ..compute norm
       nrmc(iter) = samax(m*n,wrk2,1)
*    ..update sol
       if (mod(itpar+iter,2).eq.0) then
        call saxpy (m*n, rxf, wrk2, 1, sol, 1)
       else
        do j = 0, n-1
         call saxpy (m, rxf, wrk2(0,j), 1, sol(0,n-1-j), 1)
        enddo
       endif
*   ..end 'do iter'
      enddo

* ..return
      return

*-----------------------------------------------------------------------
*.end iluter

      end subroutine iluter
