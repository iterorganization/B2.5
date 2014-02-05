      function machsfr ()
      use b2mod_types
      implicit none
      real (kind=R8) ::
     *   machsfr
*     ------------------------------------------------------------------
*     MACHSFR returns the safe numerical range (a tiny positive machine
*     value e0 for which -e0, 1/e0, -1/e0 are all computable).
*!!   This routine may contain system dependencies.
*     ------------------------------------------------------------------
      logical lshort
      parameter (lshort=1.0_R8+1.0e-12_R8.eq.1.0_R8)
      integer small
      real (kind=R8) ::
     *   mchsfr
      equivalence (mchsfr,small)
      real (kind=R8) ::
     *   x02amf
      external x02amf
      data small/8388614/
*$$$  data small/8388608/ ! this value is given in r1mach, see below.
*     ------------------------------------------------------------------
*     This code is tricky. The parameter lshort is intended to indicate
*     whether real*8 variables are four-byte or eight-byte. In the case
*     of four-byte reals it is assumed that we have IEEE arithmetic,
*     and the real*8 value mchsfr is used (mchsfr was taken from code for
*     the subroutine r1mach in the core section on netlib, April 1992).
*     In the case of eight-byte reals the NAG library routine x02amf is
*     employed.
*     It was found that if mchsfr is set exactly as in r1mach, i.e.
*     using small=8388608, then the expression exp(log(mchsfr)) raises
*     an underflow exception on the Sparc (f77 with libraries SC1.0).
*     Using a slightly larger value, small=8388614, avoids that problem.
*     The decimal equivalent value is about 1.1755e-38 in either case.
*     ------------------------------------------------------------------
      if (lshort) then
       machsfr = mchsfr
      else
       machsfr = x02amf()
      endif
      return
*     ------------------------------------------------------------------
      end
