*-----------------------------------------------------------------------
*.specification

      subroutine slv5pt (m, n, aa5, rhs, iters, sol, nrmc, k, aafac,
     &  wrk0, iwrk)
      use b2mod_types
      implicit none
*   ..input arguments (unchanged on exit)
      integer m, n, iters, k
      real (kind=R8) ::
     *   aa5(0:m-1,0:n-1,-2:2), rhs(0:m-1,0:n-1)
*   ..output arguments (unspecified on entry)
      real (kind=R8) ::
     *   sol(0:m-1,0:n-1), nrmc(0:iters-1)
*   ..workspace arguments (unspecified on entry and on exit)
      real (kind=R8) ::
     *   aafac(-2*k:k,0:m*n-1), wrk0(0:m*n-1)
      integer iwrk(0:m*n-1)

*-----------------------------------------------------------------------
*.documentation
*
*  1. purpose
*
*     SLV5PT solves a five-point system of equations by a direct
*     method.
*
*
*  3. description
*
*     This routine employs a direct band matrix solver to solve a
*     five-point system of equations. The system is
*       aa5(i,j,-2)*sol(i,j-1)+aa5(i,j,-1)*sol(i-1,j)+
*         aa5(i,j,0)*sol(i,j)+aa5(i,j,1)*sol(i+1,j)+
*         aa5(i,j,2)*sol(i,j+1) = rhs(i,j)
*     for (i,j) in (0:m-1,0:n-1); any term involving indices outside
*     (0:m-1,0:n-1) is to be understood as 0.
*
*     The matrix aa5 is expanded into a dense band matrix aafac, which
*     is decomposed into lower and upper factors. The factored matrix
*     is employed for one or more iterations to obtain the solution.
*     The number of iterations is specified by the parameter iters and
*     should normally be set to one.
*
*
*  5. parameters
*
*     m, n - integer, input.
*     m and n specify the number of mesh points along the first and
*     the second coordinate direction, respectively. The mesh points
*     are indexed by (0:m-1,0:n-1). It will hold that 2.le.m and
*     2.le.n.
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
*     iters specifies the number of iterations. This parameter should
*     normally be set to 1.
*
*     sol - (0:m-1,0:n-1) real*8 array, output.
*     sol specifies the solution.
*
*     nrmc - (0:iters-1) real*8 array, output.
*     For 0.le.it.lt.iters, nrmc(it) specifies the L1-norm of the
*     correction at iteration number (it).
*     
*     k - integer, input.
*     k is used in the dimension declaration of aafac.  It is required
*     that min(m,n).le.k.
*
*     aafac - (-2*k:k,0:m*n-1) real*8 array, workspace.
*     aafac is used to store aa5 in dense band storage mode, with
*     min(m,n) additional lower bands that are needed as space for the
*     L*U factorisation because of possible pivoting.
*
*     wrk0 - (0:m*n-1) real*8 array, workspace.
*     wrk0 is used to store the residual during the iteration.
*
*     iwrk - (0:m*n-1) integer array, workspace.
*     iwrk is used to store pivoting data.
*
*
*-----------------------------------------------------------------------
*.declarations

*   ..local scalars
      integer it, i, j, ib, ip
      integer info
      real (kind=R8) ::
     *   rxf
      parameter (rxf=1.0_R8)
*   ..procedures
      integer indx0, indx1
      real (kind=R8) ::
     *   smin, samax
      external xertst, smin, samax, sfill, sgbtrf, sgbtrs
      indx0(i,j) = i+m*j
      indx1(i,j) = n*i+j

*-----------------------------------------------------------------------
*.computation

* ..preliminaries
*   ..subprogram start-up calls
      call subini ('slv5pt')
*   ..test m, n, iters
      if (m.lt.2.or.n.lt.2) then
       call xerrab ('slv5pt--faulty argument m, n')
      else if (iters.le.0) then
       call xerrab ('slv5pt--faulty argument iters')
      endif
*   ..test edge values of side-diagonals
      if (0.lt.samax(m,aa5(0,0,-2),1).or.
     &  0.lt.samax(m,aa5(0,n-1,2),1)) then
       call xerrab ('slv5pt--faulty argument aa5(,,-2), aa5(,,2)')
      else if (0.lt.samax(n,aa5(0,0,-1),m).or.
     &   0.lt.samax(n,aa5(m-1,0,1),m)) then
       call xerrab ('slv5pt--faulty argument aa5(,,-1), aa5(,,1)')
      endif
*   ..test sign of main diagonal
      if (smin(m*n,aa5(0,0,0),1).lt.0) then
       call xerrab ('slv5pt--faulty argument sign aa5(,,0)')
      endif

* ..compute iteration matrix
*   ..branch on ordering
      if (m.le.n) then
*    ..use natural ordering
       do ip = 0, m*n-1
        do ib = -2*m, m
         aafac(ib,ip) = 0
        enddo
       enddo
       do j = 1, n-1
        do i = 0, m-1
         aafac(-m,indx0(i,j)) = aa5(i,j-1,2)
        enddo
       enddo
       do j = 0, n-1
        do i = 1, m-1
         aafac(-1,indx0(i,j)) = aa5(i-1,j,1)
        enddo
       enddo
       do j = 0, n-1
        do i = 0, m-1
         aafac(0,indx0(i,j)) = aa5(i,j,0)
        enddo
       enddo
       do j = 0, n-1
        do i = 0, m-2
         aafac(1,indx0(i,j)) = aa5(i+1,j,-1)
        enddo
       enddo
       do j = 0, n-2
        do i = 0, m-1
         aafac(m,indx0(i,j)) = aa5(i,j+1,-2)
        enddo
       enddo
*    ..decompose aafac
*      (LINPACK call has been commented out; we now use LAPACK)
*$$$   call sgbco (aafac(-2*m,0), 3*k+1, m*n, m, m, iwrk, rcond, wrk0)
*$$$   if (rcond.eq.0) then
*$$$    call xerrab ('slv5pt--error, singular matrix')
*$$$   endif
       call sgbtrf (m*n, m*n, m, m, aafac(-2*m,0), 3*k+1, iwrk, info)
       call xertst (info.eq.0, 'slv5pt--error return from sgbtrf')
      else
*    ..use ordering interchanging the two coordinates
       do ip = 0, m*n-1
        do ib = -2*n, n
         aafac(ib,ip) = 0
        enddo
       enddo
       do j = 1, n-1
        do i = 0, m-1
         aafac(-1,indx1(i,j)) = aa5(i,j-1,2)
        enddo
       enddo
       do j = 0, n-1
        do i = 1, m-1
         aafac(-n,indx1(i,j)) = aa5(i-1,j,1)
        enddo
       enddo
       do j = 0, n-1
        do i = 0, m-1
         aafac(0,indx1(i,j)) = aa5(i,j,0)
        enddo
       enddo
       do j = 0, n-1
        do i = 0, m-2
         aafac(n,indx1(i,j)) = aa5(i+1,j,-1)
        enddo
       enddo
       do j = 0, n-2
        do i = 0, m-1
         aafac(1,indx1(i,j)) = aa5(i,j+1,-2)
        enddo
       enddo
*    ..decompose aafac
*      (LINPACK call has been commented out; we now use LAPACK)
*$$$   call sgbco (aafac(-2*n,0), 3*k+1, m*n, n, n, iwrk, rcond, wrk0)
*$$$   if (rcond.eq.0) then
*$$$    call xerrab ('slv5pt--error, singular matrix')
*$$$   endif
       call sgbtrf (m*n, m*n, n, n, aafac(-2*n,0), 3*k+1, iwrk, info)
       call xertst (info.eq.0, 'slv5pt--error return from sgbtrf')
      endif

* ..perform the iterations
*   ..initialise sol
      call sfill (m*n, 0.0_R8, sol, 1)
*   ..iterate
      do it = 0, iters-1
*    ..compute the residual
       if (m.le.n) then
*     ..natural order
        do j = 0, n-1
         do i = 0, m-1
          wrk0(indx0(i,j)) = rhs(i,j)
         enddo
        enddo
        if (0.lt.it) then
         do j = 0, n-1
          wrk0(indx0(0,j)) = wrk0(indx0(0,j))-aa5(0,j,0)*sol(0,j)
          do i = 1, m-1
           wrk0(indx0(i,j)) = wrk0(indx0(i,j))-aa5(i,j,0)*sol(i,j)-
     &       aa5(i,j,-1)*sol(i-1,j)
           wrk0(indx0(i-1,j)) = wrk0(indx0(i-1,j))-aa5(i-1,j,1)*sol(i,j)
          enddo
         enddo
         do j = 1, n-1
          do i = 0, m-1
           wrk0(indx0(i,j)) = wrk0(indx0(i,j))-aa5(i,j,-2)*sol(i,j-1)
           wrk0(indx0(i,j-1)) = wrk0(indx0(i,j-1))-aa5(i,j-1,2)*sol(i,j)
          enddo
         enddo
        endif
       else
*     ..inverted order
        do j = 0, n-1
         do i = 0, m-1
          wrk0(indx1(i,j)) = rhs(i,j)
         enddo
        enddo
        if (0.lt.it) then
         do j = 0, n-1
          wrk0(indx1(0,j)) = wrk0(indx1(0,j))-aa5(0,j,0)*sol(0,j)
          do i = 1, m-1
           wrk0(indx1(i,j)) = wrk0(indx1(i,j))-aa5(i,j,0)*sol(i,j)-
     &       aa5(i,j,-1)*sol(i-1,j)
           wrk0(indx1(i-1,j)) = wrk0(indx1(i-1,j))-aa5(i-1,j,1)*sol(i,j)
          enddo
         enddo
         do j = 1, n-1
          do i = 0, m-1
           wrk0(indx1(i,j)) = wrk0(indx1(i,j))-aa5(i,j,-2)*sol(i,j-1)
           wrk0(indx1(i,j-1)) = wrk0(indx1(i,j-1))-aa5(i,j-1,2)*sol(i,j)
          enddo
         enddo
        endif
       endif
*    ..solve the equation
       if (m.le.n) then
*     ..natural order
*       (LINPACK call has been commented out)
*$$$    call sgbsl (aafac(-2*m,0), 3*k+1, m*n, m, m, iwrk, wrk0, 0)
        call sgbtrs ('N', m*n, m, m, 1, aafac(-2*m,0), 3*k+1, iwrk,
     &    wrk0, m*n, info)
        call xertst (info.eq.0, 'slv5pt--error return from sgbtrs')
       else
*     ..inverted order
*       (LINPACK call has been commented out)
*$$$    call sgbsl (aafac(-2*n,0), 3*k+1, m*n, n, n, iwrk, wrk0, 0)
        call sgbtrs ('N', m*n, n, n, 1, aafac(-2*n,0), 3*k+1, iwrk,
     &    wrk0, m*n, info)
        call xertst (info.eq.0, 'slv5pt--error return from sgbtrs')
       endif
*    ..compute norm
       nrmc(it) = samax(m*n,wrk0,1)
*    ..update sol
       if (m.le.n) then
        do j = 0, n-1
         do i = 0, m-1
          sol(i,j) = sol(i,j)+rxf*wrk0(indx0(i,j))
         enddo
        enddo
       else
        do j = 0, n-1
         do i = 0, m-1
          sol(i,j) = sol(i,j)+rxf*wrk0(indx1(i,j))
         enddo
        enddo
       endif
      enddo

* ..return
      call subend ()
      return

*-----------------------------------------------------------------------
*.end slv5pt

      end subroutine slv5pt
