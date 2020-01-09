!mb  20/07/2017 New routine intp_2dtable for bilinear interpolation of tables based on eirene intp_adas routine

!pb  070109  ff(2,2) introduced to avoid array temporaries

cdr  oct.15  comments. Parameters ip1,ip2 added, to indicate extrapolation

      function intp_2dtable(tab,r1,r2,dr1,dr2,nr1,nr2,p1,p2,ir1,
     &     ir2,ip1,ip2) result(res)
c  bilinear interpolation in 2 parameter table tab
c  linear extrapolation at all 4 parameter boundaries

c input:
c  table: tab(r1,r2)
c  r1(0:nr1):first variable
c  r2(0:nr2):second variable
c  dr1(0:nr1-1),dr2(0:nr2-1); inverse of table spacing
c  parameters:  p1, p2

c output:
c  res,  interpolated value at (p1,p2)
c  ir1,ir2:  Lowest indices of interpolated table values
c  ip1 = 0:  parameter p1 is in valid range of table
c  ip1 =-1:  parameter p1 is below left boundary of tabulated range
c  ip1 =+1:  parameter p1 is above right boundary of tabulated range

c  ip2    :  same as ip1, for parameter p2.

 
      use B2MOD_TYPES
 
      implicit none
      integer,intent(in):: nr1,nr2
      real(kind=R8), intent(in) :: tab(0:nr1,0:nr2),r1(0:nr1),r2(0:nr2),
     &   dr1(0:nr1-1),dr2(0:nr2-1)
      real(kind=R8), intent(in) :: p1,p2
      integer, intent(out) :: ir1,ir2,ip1,ip2 
      integer:: ne1,ne2
      real(kind=R8) :: res, rx, ry, ff(2,2)

c  Arrays shifted from zero-based to 1-based indexing in subroutines!
      ne1 = nr1+1
      ne2 = nr2+1
c  first parameter, find cell ir1, ir1+1 
      if (p1 <= r1(0)) then
        ir1 = 0
        ip1 = -1
      else if (p1 >= r1(nr1)) then
        ir1 = nr1-1
        ip1 = 1
      else
        call binsearch_2 (r1, ne1, p1, ir1)
        ir1 = ir1-1
        ip1= 0
      end if
c  second parameter, find cell ir2, ir2+1  
      if (p2 <= r2(0)) then
        ir2 = 0
        ip2 = -1
      else if (p2 >= r2(nr2)) then
        ir2 = nr2-1
        ip2 = 1
      else
        call binsearch_2 (r2, ne2, p2, ir2)
        ir2 = ir2-1
        ip2 = 0
      end if

c  both cell indices ir1 and ir2 are found. 
c  the parameter (p1,p2) is either in cell [ir1, ir1+1, ir2, ir2+1],
c  or outside the parameter range of this table.

 
      rx = (r1(ir1+1) - p1) * dr1(ir1)
      ry = (r2(ir2+1) - p2) * dr2(ir2)
c  fill ff: intermediate 2 by 2 matrix of tabulated values at cell vertices
      ff = tab(ir1:ir1+1,ir2:ir2+1)

      call bilinear_int (ff, rx, ry, res)
 
      return
      end function intp_2dtable

      function intp_2dtable0(tab,r1,r2,dr1,dr2,nr1,nr2,p1,p2,ir1,
     &     ir2,ip1,ip2) result(res)
c  interpolates a table for which the closest indices in the table are already known from a prior application of intp_2dtable
c  bilinear interpolation in 2 parameter table tab
c  linear extrapolation at all 4 parameter boundaries

c input:
c  ff = tab(ir1:ir1+1,ir2:ir2+1)
c  r1(0:nr1):first variable
c  r2(0:nr2):second variable
c  dr1(0:nr1-1),dr2(0:nr2-1); inverse of table spacing
c  parameters:  p1, p2
c  ir1,ir2:  Lowest indices of table values to be interpolated

c output:
c  res,  interpolated value at (p1,p2)
c  ip1 = 0:  parameter p1 is in valid range of table
c  ip1 =-1:  parameter p1 is below left boundary of tabulated range
c  ip1 =+1:  parameter p1 is above right boundary of tabulated range

c  ip2    :  same as ip1, for parameter p2.

 
      use B2MOD_TYPES
 
      implicit none
      integer,intent(in):: nr1,nr2,ir1,ir2
      real(kind=R8), intent(in) :: tab(0:nr1,0:nr2),r1(0:nr1),r2(0:nr2),
     &   dr1(0:nr1-1),dr2(0:nr2-1)
      real(kind=R8), intent(in) :: p1,p2
      integer, intent(in) :: ip1,ip2 
      real(kind=R8) :: res, rx, ry,ff(2,2)

c  both cell indices ir1 and ir2 are found. 
c  the parameter (p1,p2) is either in cell [ir1, ir1+1, ir2, ir2+1],
c  or outside the parameter range of this table.

 
      rx = (r1(ir1+1) - p1) * dr1(ir1)
      ry = (r2(ir2+1) - p2) * dr2(ir2)
c  fill ff: intermediate 2 by 2 matrix of tabulated values at cell vertices
      ff = tab(ir1:ir1+1,ir2:ir2+1)

      call bilinear_int (ff, rx, ry, res)
 
      return
      end function intp_2dtable0
