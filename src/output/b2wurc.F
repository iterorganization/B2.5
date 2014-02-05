*-----------------------------------------------------------------------
*.specification

      subroutine b2wurc (nput, version)
      use b2mod_types
      use b2mod_b2cmrc
      implicit none
*   ..input arguments (unchanged on exit)
      integer nput
      character*10 version
*   ..output arguments (unspecified on entry)
*     (none)
*   ..common blocks
*-----------------------------------------------------------------------
*.documentation
*
*  1. purpose
*
*     B2WURC writes a table of atomic rate coefficients to an
*     un*formatted file.
*
*
*  3. description (see also routine b2cdca)
*
*     This routine writes to unit (nput) a table of atomic rate
*     coefficients for the processes of ionisation, recombination,
*     electron heat loss, charge exchange and radiation with 
*     neutral hydrogen.
*     Unit nput must be connected to an un*formatted file.
*
*
*  5. parameters (see also routine b2cdcv)
*
*     nput - integer, input.
*     nput specifies a unit number for un*formatted output.
*
*     (parameters in common:)
*
*     /b2cmrc/ - input.
*     /b2cmrc/ contains a table of atomic rate coefficients.
*     See the calling routine for description.
*
*
*  8. further comments
*
*     See routine b2cdca for explanation of 'un*formatted file'.
*
*-----------------------------------------------------------------------
*.declarations

*   ..local variables
      logical lt0
      integer k, k1, idum(0:9)
      real (kind=R8) ::
     *   cutlo, t0
      real (kind=R8), allocatable :: rtza(:), rtza2(:)
*   ..procedures
      intrinsic log
      logical lnext
      real (kind=R8) ::
     *   machsfr, smax
      external subini, subend, xertst, machsfr, cfverw, cfwure, cfwuin
      external smax
      lnext(k,k1) = rtzmax(k).lt.rtzmin(k1).and.rtzn(k).eq.rtzn(k1)

*-----------------------------------------------------------------------
*.computation

* ..preliminaries
*   ..subprogram start-up calls
      call subini ('b2wurc')
*   ..set cutlo
      cutlo = machsfr()
*   ..test nput
      call xertst (1.le.nput, 'faulty argument nput')
*   ..test /b2cmrc/ dimensions
      call xertst (1.le.rtns,
     &  'faulty /b2cmrc/ rtns')
*   ..test /b2cmrc/ rtzmin, rtzmax, rtzn
      do k = 0, rtns-1
       call xertst (0.0_R8.le.rtzmin(k).and.0.0_R8.lt.rtzn(k).and.
     &   rtzmin(k).le.rtzmax(k).and.rtzmax(k).le.rtzn(k),
     &   'faulty /b2cmrc/ rtzmin, rtzmax, rtzn')
      enddo

* ..test rate coefficients
*   ..test ionisation rate coefficients
      do k = 0, rtns-1
       if (k.eq.rtns-1) then
        lt0 = .true.
       else
        lt0 = .not.lnext(k,k+1)
       endif
*      (the meaning of lt0 is 'species k cannot ionise to k+1')
       if (lt0) then
*       (ionisation rate coefficient must be 0)
        t0 = smax((rtnt+1)*(rtnn+1),rtlsa(0,0,k),1)
        call xertst (t0.lt.log(cutlo), 'faulty entry in rtlsa table')
       endif
      enddo
*   ..test recombination rate coefficients
      do k = 0, rtns-1
       if (k.eq.0) then
        lt0 = .true.
       else
        lt0 = .not.lnext(k-1,k)
       endif
*      (the meaning of lt0 is 'species k can not recombine to k-1')
       if (lt0) then
*       (recombination rate coefficient must be 0)
        t0 = smax((rtnt+1)*(rtnn+1),rtlra(0,0,k),1)
        call xertst (t0.lt.log(cutlo), 'faulty entry in rtlra table')
       endif
      enddo
*   ..test charge exchange rate coefficients
      do k = 0, rtns-1
       if (k.eq.0) then
        lt0 = .true.
       else
        lt0 = .not.lnext(k-1,k)
       endif
*      (the meaning of lt0 is 'species k can not recombine to k-1')
       if (lt0) then
*       (charge exchange rate coefficient must be 0)
        t0 = smax((rtnt+1)*(rtnn+1),rtlcx(0,0,k),1)
        call xertst (t0.lt.log(cutlo), 'faulty entry in rtlcx table')
       endif
      enddo

* ..save the table of rate coefficients
*   ..save dimensions
      idum(0) = rtnt
      idum(1) = rtnn
      idum(2) = rtns
      call cfwuin (nput, 3, idum, 'rtnt,rtnn,rtns')
*   ..save species data
      if (version.lt.'01.001.024') then
        allocate(rtza(0:rtns-1))
        allocate(rtza2(0:rtns-1))
        rtza(0:rtns-1) = (rtzmin(0:rtns-1)+rtzmax(0:rtns-1))/2.0_R8
        rtza2(0:rtns-1) = 
     &                  ((rtzmin(0:rtns-1)+rtzmax(0:rtns-1))/2.0_R8)**2
        call cfwure (nput, rtns, rtza, 'rtza')
        call cfwure (nput, rtns, rtzn, 'rtzn')
        call cfwure (nput, rtns, rtza2, 'rtza2')
        deallocate(rtza, rtza2)
      else
        call cfwure (nput, rtns, rtzmin, 'rtzmin')
        call cfwure (nput, rtns, rtzmax, 'rtzmax')
        call cfwure (nput, rtns, rtzn, 'rtzn')
      endif
*   ..save abscissae
      call cfwure (nput, rtnt+1, rtt, 'rtt')
      call cfwure (nput, rtnn+1, rtn, 'rtn')
      call cfwure (nput, rtnt+1, rtlt, 'rtlt')
      call cfwure (nput, rtnn+1, rtln, 'rtln')
*   ..save rate coefficients
      call cfwure (nput, (rtnt+1)*(rtnn+1)*rtns, rtlsa, 'rtlsa')
      call cfwure (nput, (rtnt+1)*(rtnn+1)*rtns, rtlra, 'rtlra')
      call cfwure (nput, (rtnt+1)*(rtnn+1)*rtns, rtlqa, 'rtlqa')
      call cfwure (nput, (rtnt+1)*(rtnn+1)*rtns, rtlcx, 'rtlcx')
      if (version.ge.'01.001.016') then
        call cfwure (nput, (rtnt+1)*(rtnn+1)*rtns, rtlrd, 'rtlrd')
      end if
      if (version.ge.'01.001.019') then
        call cfwure (nput, (rtnt+1)*(rtnn+1)*rtns, rtlbr, 'rtlbr')
      end if
      if (version.ge.'01.001.024') then
        call cfwure (nput, (rtnt+1)*(rtnn+1)*rtns, rtlza, 'rtlza')
        call cfwure (nput, (rtnt+1)*(rtnn+1)*rtns, rtlz2, 'rtlz2')
        call cfwure (nput, (rtnt+1)*(rtnn+1)*rtns, rtlpt, 'rtlpt')
        call cfwure (nput, (rtnt+1)*(rtnn+1)*rtns, rtlpi, 'rtlpi')
      end if

* ..return
      call subend ()
      return

*-----------------------------------------------------------------------
*.end b2wurc

      end subroutine b2wurc
