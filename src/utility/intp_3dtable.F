

      function intp_3dtable(tab,r1,r2,r3,dr1,dr2,dr3,
     &     nr1,nr2,nr3,p1,p2,p3,ir1,
     &     ir2,ir3,ip1,ip2,ip3) result(res)
c  linear interpolation in 3 parameter table tab
c  linear extrapolation at all 8 parameter boundaries

c input:
c  table: tab(r1,r2,r3)
c  r1(0:nr1):first variable
c  r2(0:nr2):second variable
c  r3(0:nr3):third variable
c  dr1(0:nr1-1),dr2(0:nr2-1); dr3(0:nr3-1); inverse of table spacing
c  parameters:  p1, p2, p3

c output:
c  res,  interpolated value at (p1,p2)
c  ir1,ir2,ir3:  Lowest indices of interpolated table values
c  ip1 = 0:  parameter p1 is in valid range of table
c  ip1 =-1:  parameter p1 is below left boundary of tabulated range
c  ip1 =+1:  parameter p1 is above right boundary of tabulated range

c  ip2    :  same as ip1, for parameter p2.
c  ip3    :  same as ip1, for parameter p3.

 
      use B2MOD_TYPES
 
      implicit none
      integer,intent(in):: nr1,nr2,nr3
      real(kind=R8), intent(in) :: tab(0:nr1,0:nr2,0:nr3),
     &   r1(0:nr1),r2(0:nr2),r3(0:nr3),
     &   dr1(0:nr1-1),dr2(0:nr2-1),dr3(0:nr3-1)
      real(kind=R8), intent(in) :: p1,p2,p3
      integer, intent(out) :: ir1,ir2,ir3,ip1,ip2,ip3
      integer:: ne1,ne2,ne3
      real(kind=R8) :: res, rx, ry, rz, ff(2,2,2)

c  Arrays shifted from zero-based to 1-based indexing in subroutines!
      ne1 = nr1+1
      ne2 = nr2+1
      ne3 = nr3+1
c  first parameter, find cell ir1, ir1+1 
      if (p1 < r1(0)) then
        ir1 = 0
        ip1 = -1
      else if (p1 > r1(nr1)) then
        ir1 = nr1-1
        ip1 = 1
      else
        call binsearch_2 (r1, ne1, p1, ir1)
        ir1 = ir1-1
        ip1= 0
      end if
c  second parameter, find cell ir2, ir2+1  
      if (p2 < r2(0)) then
        ir2 = 0
        ip2 = -1
      else if (p1 > r1(nr1)) then
        ir2 = nr2-1
        ip2 = 1
      else
        call binsearch_2 (r2, ne2, p2, ir2)
        ir2 = ir2-1
        ip2 = 0
      end if
c  third parameter, find cell ir3, ir3+1  
      if (p3 < r3(0)) then
        ir3 = 0
        ip3 = -1
      else if (p1 > r1(nr1)) then
        ir3 = nr3-1
        ip3 = 1
      else
        call binsearch_2 (r3, ne3, p3, ir3)
        ir3 = ir3-1
        ip3 = 0
      end if

c  cell indices ir1, ir2 and ir3 are found. 
c  the parameter (p1,p2) is in cell [ir1, ir1+1, ir2, ir2+1, ir3, ir3+1].

      rx = (r1(ir1+1) - p1) * dr1(ir1)
      ry = (r2(ir2+1) - p2) * dr2(ir2)
      rz = (r3(ir3+1) - p3) * dr3(ir3)
c  fill ff: intermediate 2 by 2 by 2 matrix of tabulated values at cell vertices
      ff = tab(ir1:ir1+1,ir2:ir2+1,ir3:ir3+1)

      call trilinear_int (ff, rx, ry, rz, res)
 
      return
      end function intp_3dtable

      function intp_3dtable0(tab,r1,r2,r3,dr1,dr2,dr3,
     &     nr1,nr2,nr3,p1,p2,p3,ir1,
     &     ir2,ir3) result(res)
c  interpolates a table for which the closest indices in the table are already known from a prior application of intp_3dtable
c  trilinear interpolation in 3 parameter table tab

c input:
c  ff = tab(ir1:ir1+1,ir2:ir2+1,ir3:ir3+1)
c  r1(0:nr1):first variable
c  r2(0:nr2):second variable
c  r3(0:nr3):third variable
c  dr1(0:nr1-1),dr2(0:nr2-1),dr3(0:nr3-1); inverse of table spacing
c  parameters:  p1, p2, p3
c  ir1,ir2,ir3:  Lowest indices of table values to be interpolated

c output:
c  res,  interpolated value at (p1,p2,p3)


      use B2MOD_TYPES
 
      implicit none
      integer,intent(in):: nr1,nr2,nr3,ir1,ir2,ir3
      real(kind=R8), intent(in) :: tab(0:nr1,0:nr2,0:nr3),
     &   r1(0:nr1),r2(0:nr2),r3(0:nr3),
     &   dr1(0:nr1-1),dr2(0:nr2-1),dr3(0:nr3-1)
      real(kind=R8), intent(in) :: p1,p2,p3
      real(kind=R8) :: res, rx, ry,rz,ff(2,2,2)

c  cell indices ir1, ir2 and ir3 are know from a previous call to intp_3dtable are 
c  the parameter (p1,p2,p3) is in cell [ir1, ir1+1, ir2, ir2+1, ir3, ir3+1].


      rx = (r1(ir1+1) - p1) * dr1(ir1)
      ry = (r2(ir2+1) - p2) * dr2(ir2)
      rz = (r3(ir3+1) - p3) * dr3(ir3)
c  fill ff: intermediate 2 by 2 matrix of tabulated values at cell vertices
      ff = tab(ir1:ir1+1,ir2:ir2+1,ir3:ir3+1)

      call trilinear_int (ff, rx, ry, rz, res)
 
      return
      end function intp_3dtable0
