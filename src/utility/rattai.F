*-----------------------------------------------------------------------
*.specification

      subroutine rattai (nt, ndum, tailep, rtt, rtorg, rtmod)
      use b2mod_types
      use b2mod_math
      use b2mod_subsys
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
     *   t0, t1, t2, rtm, rt0, rt1, rtt0, rtt1, ds
*   ..procedures
      intrinsic log
      external xertst

*-----------------------------------------------------------------------
*.computation

* ..preliminaries
      call subini ('rattai')
*   ..set cutll, cutlo
      if (.not.b2mod_math_initialised) call b2mod_math_init
*   ..test nt, ndum, tailep
      call xertst (0.le.nt, 'rattai--faulty argument nt')
      call xertst (1.le.ndum, 'rattai--faulty argument ndum')
      call xertst (0.0_R8.le.tailep, 'rattai--faulty argument tailep')
*   ..test rtt, rtorg
      do i = 0, nt
       call xertst (0.0_R8.lt.rtt(i),
     .  'rattai--faulty argument array rtt')
      enddo
      do k = 0, ndum-1
       do i = 0, nt
        call xertst (0.0_R8.le.rtorg(i,k),
     .   'rattai--faulty argument array rtorg')
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
          rt0 = rtorg(j,k)-rtorg(i,k)
          rt1 = rtorg(j+1,k) - rtorg(i,k)
          if (rt0.gt.0.0_R8) then
            rtt0 = log(rt0)+log(t0)
          else
            rtt0 = cutll
          endif
          if (rt1.gt.0.0_R8) then
            rtt1 = log(rt1)+log(t1)
          else
            rtt1 = cutll
          endif
          if (rtt0.gt.cutll) then
            rtt0 = expu(rtt0)
          else
            rtt0 = 0.0_R8
          end if
          if (rtt1.gt.cutll) then
            rtt1 = expu(rtt1)
          else
            rtt1 = 0.0_R8
          end if
          if (rtt0+rtt1.gt.cutlo) then
            t2 = log(ds/2.0_R8)+log(rtt0+rtt1)
          else
            t2 = cutll
          end if
          if (t2.gt.cutll) then
            rtm = exp(t2)
          else
            rtm = 0.0_R8
          endif
          rtmod(i,k) = rtmod(i,k)+rtm
         enddo
        enddo
       endif
      enddo

* ..return
      call subend ()
      return

*-----------------------------------------------------------------------
*.end rattai

      end subroutine rattai
