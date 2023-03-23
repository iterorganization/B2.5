      function macheps ()
      use b2mod_types
      implicit none
      real (kind=R8) ::
     *   macheps
*     ------------------------------------------------------------------
*     MACHEPS returns the machine precision.
*!!   This routine may contain system dependencies.
*     ------------------------------------------------------------------
      logical lshort
      parameter (lshort=1.0_R8+1.0e-12_R8.eq.1.0_R8)
      integer diver
      real (kind=R8) ::
     *   mcheps
      equivalence (mcheps,diver)
      real (kind=R8) ::
     *   x02ajf
      external x02ajf
      data diver/872415232/
*     ------------------------------------------------------------------
*     This code is tricky. The parameter lshort is intended to indicate
*     whether real*8 variables are four-byte or eight-byte. In the case
*     of four-byte reals it is assumed that we have IEEE arithmetic,
*     and the real*8 value mcheps is used (mcheps was taken from code for
*     the subroutine r1mach in the core section on netlib, April 1992).
*     In the case of eight-byte reals the NAG library routine x02ajf is
*     employed.
*     ------------------------------------------------------------------
      if (lshort) then
       macheps = mcheps
      else
       macheps = x02ajf()
      endif
      return
*     ------------------------------------------------------------------
      end function macheps
