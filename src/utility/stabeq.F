*-----------------------------------------------------------------------
*.begin stabeq

      subroutine stabeq (m, n, rf, aa5)
      use b2mod_types
      implicit none
*   ..input arguments (unchanged on exit)
      integer m, n
      real (kind=R8) ::
     *   rf
*   ..input/output arguments
      real (kind=R8) ::
     *   aa5(0:m-1,0:n-1,-2:2)

*-----------------------------------------------------------------------
*.documentation
*
*  1. purpose
*
*     STABEQ adjusts the diagonal of a five-point matrix to satisfy
*     a property related to diagonal dominance.
*
*
*  3. description
*
*     This routine modifies the diagonal of the five-point matrix aa5
*     to ensure that the expression
*       aa5(i,j,0)+rf*(aa5(i,j,-2)+aa5(i,j,-1)+aa5(i,j,1)+aa5(i,j,2))
*     is everywhere non-negative.
*     (If rf.eq.1.0 this amounts to aa5 being diagonally dominant by
*     rows.)
*
*
*  5. parameters
*
*     m, n - integer, input.
*     m and n specify the size of the rectangular mesh on which the
*     five-point operator is defined. This mesh has m points along the
*     first coordinate direction and n points along the second
*     coordinate direction.
*     m and n are unchanged on exit.
*
*     rf - real, input.
*     rf specifies the coefficient described in Section 3. It will hold
*     that 0.le.rf.
*     rf is unchanged on exit.
*
*     aa5(0:m-1,0:n-1,-2:2) real*8 array, input/output.
*     For (i,j) in (0:m-1,0:n-1), aa5(i,j,-2:2) specify the five
*     coefficients for the equation at the (i,j) location.
*     On entry the central coefficient, aa5(,,0), must be non-negative.
*     On exit, for (i,j) in (0:m-1,0:n-1), the central coefficient
*     aa5(i,j,0) will have been replaced by the stabilised value
*       max(aa5(i,j,0),
*         -rf*(aa5(i,j,-2)+aa5(i,j,-1)+aa5(i,j,1)+aa5(i,j,2)))
*     The other entries in aa5 are unchanged on exit.
*
*-----------------------------------------------------------------------
*.declarations

*   ..local variables
      integer i, j
*   ..procedures
      intrinsic max
      real (kind=R8) ::
     *   smin
      external xerrab, smin

*-----------------------------------------------------------------------
*.computation

*   ..test m, n
      if (m.lt.1.or.n.lt.1) then
       call xerrab ('stabeq--faulty argument m, n')
      endif
*   ..test rf
      if (rf.lt.0) then
       call xerrab ('stabeq--faulty argument sign rf')
      endif
*   ..test aa5(,,0)
      if (smin(m*n,aa5(0,0,0),1).lt.0) then
       call xerrab ('stabeq--faulty argument sign aa5(,,0)')
      endif
*   ..adjust aa5(,,0)
      do j = 0, n-1
       do i = 0, m-1
        aa5(i,j,0) = max(aa5(i,j,0),
     &    -rf*(aa5(i,j,-2)+aa5(i,j,-1)+aa5(i,j,1)+aa5(i,j,2)))
       enddo
      enddo

* ..return
      return

*-----------------------------------------------------------------------
*.end stabeq

      end subroutine stabeq
