      module b2mod_math
      use b2mod_types
      implicit none

      real(kind=R8), save :: cutll, cutlo
      real(kind=R4), save :: small_R4_constant
      logical, save :: b2mod_math_initialised = .false.

      contains

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

      integer function ibase(n,incx)
      implicit none
      integer n, incx
      intrinsic max
      ibase = max(0,(1-n)*incx)
      return
      end function ibase

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
