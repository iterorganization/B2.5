*-----------------------------------------------------------------------
*.specification

      subroutine rattai (nt, ndum, tailep, rtt, rtorg, rtmod)
      use b2mod_types
      implicit none
*   ..input arguments (unchanged on exit)
      integer nt, ndum
      real (kind=R8) ::
     *   tailep, rtt(0:nt), rtorg(0:nt,0:ndum-1)
*   ..output arguments (unspecified on entry)
      real (kind=R8) ::
     *   rtmod(0:nt,0:ndum-1)

*-----------------------------------------------------------------------
*.documentation
*
*  1. purpose
*
*     RATTAI computes a smoothed set of rate coefficients to represent
*     the effect of a power-law tail on the distribution function.
*
*
*  3. description
*
*     This routine ...
*
*
*  5. parameters
*
*     nt - integer, input.
*     nt specifies the upper dimension bound for the array rtt.
*     It will hold that 0.le.nt.
*
*     ndum - integer, input.
*     ndum specifies the number of elements in the (insignificant)
*     second dimension of rtorg and rtmod.
*     It will hold that 1.le.ndum.
*
*     tailep - real, input.
*     tailep specifies a smoothing parameter for the temperature
*     dependence of some rate coefficient; the smoothing models the
*     effect of a power-law tail on the distribution function, with
*     exponent equal to -1/tailep.
*     It is required that 0.le.tailep. tailep.eq.0 implies no smoothing.
*
*     rtt - (0:nt) real*8 array, input.
*     rtt(0:nt) specifies the sequence of values of the temperature
*     variable, expressed in ev, for which the rate coefficients are
*     evaluated. rtt(0:nt) must be positive.
*
*     rtorg - (0:nt,0:ndum-1) real*8 array, input.
*     For 0.le.k.lt.ndum, rtorg(0:nt,k) specifies a set of original
*     rate coefficients, representing a function of temperature.
*
*     rtmod - (0:nt,0:ndum-1) real*8 array, output.
*     For 0.le.k.lt.ndum, rtt(0:nt,k) specifies the smoothed rate
*     coefficients corresponding to rtorg(0:nt,k).
*
*-----------------------------------------------------------------------
*.declarations

*   ..local variables
      integer i, j, k
      real (kind=R8) ::
     *   t0, t1, ds
*   ..procedures
      intrinsic log

*-----------------------------------------------------------------------
*.computation

* ..preliminaries
*   ..test nt, ndum, tailep
      if (.not.(0.le.nt)) then
       stop 'rattai--faulty argument nt'
      else if (.not.(1.le.ndum)) then
       stop 'rattai--faulty argument ndum'
      else if (.not.(0.0_R8.le.tailep)) then
       stop 'rattai--faulty argument tailep'
      endif
*   ..test rtt, rtorg
      do i = 0, nt
       if (.not.0.0_R8.lt.rtt(i)) then
        stop 'rattai--faulty argument array rtt'
       endif
      enddo
      do k = 0, ndum-1
       do i = 0, nt
        if (.not.0.0_R8.le.rtorg(i,k)) then
         stop 'rattai--faulty argument array rtorg'
        endif
       enddo
      enddo

* ..compute smoothed rate coefficients
*   ..loop over table
      do i = 0, nt
       do k = 0, ndum-1
        rtmod(i,k) = rtorg(i,k)
       enddo
       if (0.0_R8.lt.tailep) then
        do j = i, nt-1
         t0 = (rtt(j)/rtt(i))**(-1/tailep)
         t1 = (rtt(j+1)/rtt(i))**(-1/tailep)
         ds = log(rtt(j+1)/rtt(j))
         do k = 0, ndum-1
          rtmod(i,k) = rtmod(i,k)+(ds/2.0_R8)*
     &     ((rtorg(j,k)-rtorg(i,k))*t0+
     &      (rtorg(j+1,k)-rtorg(i,k))*t1)
         enddo
        enddo
       endif
      enddo

* ..return
      return

*-----------------------------------------------------------------------
*.end rattai

      end subroutine rattai
