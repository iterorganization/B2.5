      module b2mod_math
      use b2mod_types
      implicit none

      real(kind=R8), save :: cutll, cutlo
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
      expu4 = exp(max(t0,cutll))+small_R4_constant
      return
      end function expu4

      function harm(t1,t2)
      implicit none
      real(kind=R8) :: harm, t1, t2
      harm = 2.0_R8*t1*t2/(t1+t2)
      return
      end function harm

      function hybr(flo,dif)
      implicit none
      real(kind=R8) :: hybr, flo, dif
      hybr = max(dif,abs(flo/2.0_R8),sqrt(dif**2+flo**2/4.0_R8))
      return
      end function hybr

      integer function ibase(n,incx)
      implicit none
      integer n, incx
      intrinsic max
      ibase = max(0,(1-n)*incx)
      return
      end function ibase

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

csw 2013 uwpind scheme (for modified hybrid, SOLPS4 style)
      function upwind(flo,t1,t2)
      implicit none
      real(kind=R8) :: upwind, flo, t1, t2
      upwind = max(flo,0.0_R8)*t1+min(flo,0.0_R8)*t2
      return
      end function upwind

      function uxcm(t0,t1)
      implicit none
      real(kind=R8) :: uxcm, t0, t1
      uxcm = sign(max(abs(t0),abs(t1)),t0+t1)
      return
      end function uxcm

      subroutine b2mod_math_init
      implicit none
      real(kind=R8) :: machsfr
      external machsfr
*   ..set cutll, cutlo
      cutll = log(machsfr())
      cutlo = max(machsfr(),exp(cutll))
      small_R4_constant = tiny(small_R4_constant)
      b2mod_math_initialised = .true.
      return
      end subroutine b2mod_math_init

      end module b2mod_math
