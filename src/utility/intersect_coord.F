      subroutine intersect_coord(mpg, p1, q1, p2, q2, XI, YI)
c     Subroutine to calculate the intersection of two line segments
c     Segment 1: P1 = [X1, Y1], P2 = [X2, Y2]
c     Segment 2: Q1 = [X3, Y3], Q2 = [X4, Y4]
c     Inputs:
c     P1(2), P2(2): Endpoints of the first segment
c     Q1(2), Q2(2): Endpoints of the second segment
c     Outputs:
c     XI, YI: Intersection points if they exist
c     iflag: 0 -> No intersection
c            1 -> Segments intersect at a single point
c            2 -> Segments overlap (collinear)
      use b2mod_types
      use b2us_map
      use b2mod_subsys
      type (mapping), intent(in) :: mpg
      real (kind = R8), intent(in) :: p1(2), q1(2), p2(2), q2(2)
      real (kind = R8), intent(out) :: XI, YI
      integer :: iflag
c     Internal variables
      real (kind = R8) :: denom, num1, num2, t1, t2
      real (kind = R8) :: x1, y1, x2, y2, x3, y3, x4, y4

      call subini('intersect_coord')

c Extract coordinates for clarity
      x1 = p1(1)
      y1 = p1(2)
      x2 = p2(1)
      y2 = p2(2)
      x3 = q1(1)
      y3 = q1(2)
      x4 = p2(1)
      y4 = q2(2)

c Calculate denominator
      denom = (x1-x2)*(y3-y4)-(y1-y2)*(x3-x4)

      if(abs(denom).lt.EPSILON(1.0_R8)) then
c Segments are parallel or collinear
         if((x1-x3)*(y3-y4)-(y1-y3)*(x3-x4).eq.0) then
c Collinear segments
            iflag = 2
         else
c Parallel but not collinear
            iflag = 0
         endif
         return
      endif

c Compute numerators for the parametric equations
      num1 = (x1-x3)*(y3-y4)-(y1-y3)*(x3-x4)
      num2 = (x1-x3)*(y1-y2)-(y1-y3)*(x1-x2)

c Compute parameters t1 and t2
      t1 = num1 / denom
      t2 = num2 / denom

c Check if the intersection point lies within both segments
      if(t1.ge.0.and.t1.le.1.and.t2.ge.0.and.t2.le.1) then
c Intersection exists
         XI = x1+t1*(x2-x1)
         YI = y1+t1*(y2-y1)
         iflag = 1
      else
c Intersection point is outside the segments
         iflag = 0
      endif

*   ..return
      call subend ()
      return

*-----------------------------------------------------------------------
*.end intersect_coord

      end subroutine intersect_coord
