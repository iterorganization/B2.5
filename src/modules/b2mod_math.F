      module b2mod_math
      use b2mod_types
      implicit none

      real(kind=R8), save :: cutll, cutlo, cutl4
      real(kind=R4), save :: small_R4_constant
      logical, save :: b2mod_math_initialised = .false.

*-----------------------------------------------------------------------
*.documentation-internal
*
*     Let (p0x,p0y) specify a vector p0, and (p1x,p1y) specify
*     a vector p1. Then norm(p0x,p0y) is the standard norm of p0,
*     and cosang(p0x,p0y,p1x,p1y) is cos(t), where t is the complementary
*     angle between p0 and p1.
*     t=pi/2 when p0 and p1 are parallel, t=0 or t=pi when p0 and p1
*     are orthogonal, and t=-pi/2 when p0 and p1 are anti-parallel.
*
*-----------------------------------------------------------------------

      contains

      function addup(t0,t1)
      implicit none
      real(kind=R8) :: addup, t0, t1
      addup = t0+t1
      return
      end function addup

      function cosang(p0x,p0y,p1x,p1y)
      implicit none
      real(kind=R8) :: cosang, p0x, p0y, p1x, p1y
      cosang = trim1((p0x*p1y-p0y*p1x)/(norm(p0x,p0y)*norm(p1x,p1y)))
      return
      end function cosang

      function expu(t0)
*     (expu is the exponential function with underflow set to 0)
      implicit none
      real(kind=R8) :: expu, t0
      if (.not.b2mod_math_initialised) call b2mod_math_init
      expu = dim(exp(max(t0,cutll)),cutlo)
      return
      end function expu

      function expu1(t0)
*     (expu1 is the exponential function with underflow set to a tiny
*     positive value)
      implicit none
      real(kind=R8) :: expu1, t0
      if (.not.b2mod_math_initialised) call b2mod_math_init
      expu1 = exp(max(t0,cutll))
      return
      end function expu1

      function expu2(t0)
*     (expu2 is the exponential function with underflow set to 0 and
*     large values limited at a somewhat arbitrary level.)
      implicit none
      real(kind=R8) :: expu2, t0
      if (.not.b2mod_math_initialised) call b2mod_math_init
      expu2 = dim(exp(min(20.0_R8,max(t0,cutll))),cutlo)
      return
      end function expu2

      function expu4(t0)
*     (expu4 is the exponential function with underflow set to a tiny
*     single-precision positive value)
      implicit none
      real(kind=R8) :: expu4, t0
      if (.not.b2mod_math_initialised) call b2mod_math_init
      expu4 = exp(max(t0,cutl4))
      return
      end function expu4

      integer function ibase(n,incx)
      implicit none
      integer n, incx
      intrinsic max
      ibase = max(0,(1-n)*incx)
      return
      end function ibase

      integer function ibase1(n1,inc1x)
      implicit none
      integer n1, inc1x
      intrinsic max
      ibase1 = max(1,1+(1-n1)*inc1x)
      return
      end function ibase1

! $AD NOCHECKPOINT
      function med(t0,t1,t2)
      implicit none
      real(kind=R8) :: med, t0, t1, t2
      med = min(max(t0,t1),max(t1,t2),max(t2,t0))
      return
      end function med

      function norm(p0x,p0y)
      implicit none
      real(kind=R8) :: norm, p0x, p0y
      norm = sqrt(p0x**2+p0y**2)
      return
      end function norm

#ifdef LEGACYCOMP
      function norm2(n,vector)
      use b2mod_types
      implicit none
      real (kind=R8) :: norm2
      integer, intent(in) :: n
      real (kind=R8), intent(in) :: vector(n)
      integer i

      norm2 = 0.0_R8
      do i = 1, n
        norm2 = norm2 + vector(i)*vector(i)
      end do
      norm2 = sqrt(norm2)
      return
      end function norm2
#endif

      function sinang(p0x,p0y,p1x,p1y)
      implicit none
      real(kind=R8) :: sinang, p0x, p0y, p1x, p1y
      sinang = trim1((p0x*p1x+p0y*p1y)/(norm(p0x,p0y)*norm(p1x,p1y)))
      return
      end function sinang

      function trim1(p0x)
      implicit none
      real(kind=R8) :: trim1, p0x
      trim1 = max(-1.0_R8,min(1.0_R8,p0x))
      return
      end function trim1

*     (trimg will be used to trim the logarithmic gradients to avoid
*     some pathetic behaviour near the underflow limit.)
! $AD NOCHECKPOINT
      function trimg(t0,t1,tt)
      implicit none
      real(kind=R8) :: trimg, t0, t1, tt
      trimg = max(t0,min(t1,tt))
      return
      end function trimg

      function trunc(t0,edg)
      implicit none
      real(kind=R8) :: trunc, t0, edg
      trunc = max(-edg,min(1+edg,t0))
      return
      end function trunc

      logical function ltst(t0,t1)
      implicit none
      real(kind=R8) :: t0, t1
      ltst = abs(t0-t1).le.1.0e-6_R8*(abs(t0)+abs(t1))
      return
      end function ltst

  !> Area of triangle given by vertex coordinates
      function triangleArea( x0, y0, x1, y1, x2, y2 )
      implicit none
      real (kind=R8) :: triangleArea
      real (kind=R8) :: x0, y0, x1, y1, x2, y2
      intrinsic abs
      triangleArea = 0.5_R8 * abs( (x1-x0)*(y2-y0) - (y1-y0)*(x2-x0) )
      return
      end function triangleArea

  !> Coordinates of centroid of a triangle given by vertex coordinates
      function triangleCentroid( x0, y0, x1, y1, x2, y2 )
      implicit none
      real (kind=R8) :: triangleCentroid(0:1)
      real (kind=R8) :: x0, y0, x1, y1, x2, y2
      triangleCentroid(0) = ( x0 + x1 + x2 ) / 3.0_R8
      triangleCentroid(1) = ( y0 + y1 + y2 ) / 3.0_R8
      return
      end function triangleCentroid

  !> Area of quadrangle given by vertex coordinates
  !> The vertices must be in sequence! (=> structured)
      function quadArea( x0, y0, x1, y1, x2, y2, x3, y3 )
      implicit none
      real (kind=R8) :: quadArea
      real (kind=R8) :: x0, y0, x1, y1, x2, y2, x3, y3

      quadArea = TriangleArea(x0, y0, x1, y1, x2, y2) +
     .           TriangleArea(x2, y2, x3, y3, x0, y0)

      end function quadArea

  !> Coordinates of centroid of a quadrangle given by vertex coordinates
  !> Corners are numbered in the B2 convention 2-3
  !>                                           0-1
      function quadCentroid( x0, y0, x1, y1, x2, y2, x3, y3 )
      implicit none
      real (kind=R8) :: quadCentroid(0:1)
      real (kind=R8) :: x0, y0, x1, y1, x2, y2, x3, y3

      quadCentroid = ( triangleCentroid(x0,y0,x1,y1,x3,y3) 
     &   * triangleArea(x0,y0,x1,y1,x3,y3) 
     &   + triangleCentroid(x0,y0,x2,y2,x3,y3) 
     &   * triangleArea(x0,y0,x2,y2,x3,y3) ) 
     &   / ( triangleArea(x0,y0,x1,y1,x3,y3) 
     &   + triangleArea(x0,y0,x2,y2,x3,y3) )
      return
      end function quadCentroid

      !> Area of a convex pentagonal cell
      !> The vertices must be in sequence

      function pentArea(x0, y0, x1, y1, x2, y2, x3, y3, x4, y4)
      implicit none
      real (kind=R8) :: pentArea
      real (kind=R8) :: x0, y0, x1, y1, x2, y2, x3, y3, x4, y4

      pentArea = triangleArea(x0, y0, x1, y1, x2, y2) + 
     .         triangleArea(x2, y2, x3, y3, x4, y4) +
     .         triangleArea(x2, y2, x0, y0, x4, y4)
 
      return
      end function pentArea

      !> Centroid of convex pentagonal cell
      !> The vertices must be in sequence      
      function pentCentroid( x0,y0,x1,y1,x2,y2,x3,y3,x4,y4)
      implicit none
      real (kind=R8) :: pentCentroid(0:1)
      real (kind=R8) :: x0, y0, x1, y1, x2, y2, x3, y3, x4, y4

      pentCentroid = (triangleCentroid(x0, y0, x1, y1, x2, y2)
     &  * triangleArea(x0, y0, x1, y1, x2, y2) 
     &  + triangleCentroid(x2, y2, x3, y3, x4, y4) 
     &  * triangleArea(x2, y2, x3, y3, x4, y4)  
     &  + triangleCentroid(x2, y2, x0, y0, x4, y4)
     &  * triangleArea(x2, y2, x0, y0, x4, y4) )
     &  / (triangleArea(x0, y0, x1, y1, x2, y2)  
     &  + triangleArea(x2, y2, x3, y3, x4, y4)
     &  + triangleArea(x2, y2, x0, y0, x4, y4))

      return
      end function pentCentroid

  !> Average of a quantity given on the vertices of
  !> an arbitrary quadrangular cell
  !> Corners are numbered in the B2 convention 2-3
  !>                                           0-1
      function averageVertex( f0, x0, y0, f1, x1, y1, 
     &                        f2, x2, y2, f3, x3, y3 )
      implicit none
      real (kind=R8) :: averageVertex
      real (kind=R8) :: f0, f1, f2, f3
      real (kind=R8) :: x0, y0, x1, y1, x2, y2, x3, y3
      real (kind=R8) :: a0, a1, a2, a3
  
      a0 = TriangleArea( x2, y2, x0, y0, x1, y1 )
      a1 = TriangleArea( x0, y0, x1, y1, x3, y3 )
      a2 = TriangleArea( x0, y0, x2, y2, x3, y3 )
      a3 = TriangleArea( x2, y2, x3, y3, x1, y1 )
 
      averageVertex = (f0*a0 + f1*a1 + f2*a2 + f3*a3) /  
     &                (a0 + a1 + a2 + a3)
 
      return
      end function averageVertex

      subroutine b2mod_math_init
      implicit none
      real(kind=R8) :: machsfr
      real(kind=R4) :: machsfr4
      external machsfr, machsfr4
*   ..set cutll, cutlo, cutl4, small_R4_constant
      cutll = log(machsfr())
      cutlo = max(machsfr(),exp(cutll))
      cutl4 = real(log(machsfr4()),R8)
      small_R4_constant = max(machsfr4(),real(exp(cutl4),R4))
      b2mod_math_initialised = .true.
      return
      end subroutine b2mod_math_init

      end module b2mod_math
