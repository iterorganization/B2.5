      function intersects(p1,q1,p2,q2)
      use b2mod_types
      implicit none
      logical :: intersects
      real(kind=R8) :: p1(2),q1(2),p2(2),q2(2)
***
*     This function checks if segments p1-q1 and p2-q2 intersects each other
***
      integer :: o1, o2, o3, o4, orient
      logical :: onSegment
      external orient, onSegment
cc    Find the four orientations needed for general and special cases 
      o1 = orient(p1, q1, p2)
      o2 = orient(p1, q1, q2)
      o3 = orient(p2, q2, p1)
      o4 = orient(p2, q2, q1)
      
      intersects = .true. 
cc    general case (p1, q1, p2) and (p1, q1, q2) have different orientations and
cc                 (p2, q2, p1) and (p2, q2, q1) have different orientations. 
      if (o1 .ne. o2 .and. o3 .ne. o4) return 
  
cc    special cases
cc    p1, q1 and p2 are colinear and p2 lies on segment p1q1 
      if (o1 .eq. 0 .and. onSegment(p1, p2, q1)) return 
  
cc    p1, q1 and q2 are colinear and q2 lies on segment p1q1 
      if (o2 .eq. 0 .and. onSegment(p1, q2, q1)) return
  
cc    p2, q2 and p1 are colinear and p1 lies on segment p2q2 
      if (o3 .eq. 0 .and. onSegment(p2, p1, q2)) return 
  
cc    p2, q2 and q1 are colinear and q1 lies on segment p2q2 
      if (o4 .eq. 0 .and. onSegment(p2, q1, q2)) return 
      intersects = .false.
      return 
      end function intersects

      function orient(p, q, r)
      use b2mod_types
      implicit none
      integer :: orient
      real(kind=R8) :: p(2), q(2), r(2), val

      val = (q(2) - p(2)) * (r(1) - q(1)) - 
     .      (q(1) - p(1)) * (r(2) - q(2))
 
      if (val .gt. 0.0_R8) then
        orient = 1 !clockwise
        return
      elseif (val .le. 0.0_R8) then
        orient = 2 !counterclockwise
        return
      else
        orient = 0 !colinear
        return
      endif
      end function orient

      function onSegment(p, q, r)
      use b2mod_types
      implicit none
      logical :: onSegment
      real(kind=R8) :: p(2), q(2), r(2)

      onSegment = .false.
      if ( (q(1) .le. max(p(1), r(1))) .and. 
     &     (q(1) .ge. min(p(1), r(1))) .and.
     &     (q(2) .le. max(p(2), r(2))) .and. 
     &     (q(2) .ge. min(p(2), r(2))) )
     .   onSegment = .true.
      return
      end function onSegment
