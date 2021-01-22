*-----------------------------------------------------------------------
*.specification

      subroutine b2sqcx (nCv, ns, ev, switch, amh0, th0, thn0, rlcx0, 
     &                  rcx)
      use b2mod_switches
      use b2mod_types
      implicit none
*   ..input arguments (unchanged on exit)
      type (switches), intent (in) :: switch
      integer nCv, ns
      real (kind=R8) ::
     *   ev, amh0, th0(nCv), thn0(nCv), rlcx0(nCv,0:1,0:ns-1)
*   ..output arguments (unspecified on entry)
      real (kind=R8) ::
     *   rcx(nCv,0:ns-1)

*-----------------------------------------------------------------------
*.documentation
*
*  1. purpose
*
*     B2SQCX computes on all cells of the mesh the rate coefficients
*     for charge exchange with one neutral hydrogen species.
*
*
*  3. description (see also routine b2cdca)
*
*     This routine employs the log-log linear representation provided
*     in the array rlcx0 to evaluate on all cells of the mesh the rate
*     coefficients for charge exchange with one neutral hydrogen
*     species, characterised by atomic mass amh0.
*
*
*  8. further comments
*
*     In case a rate coefficient vanishes, the representation of the
*     log of the rate coefficient will provide a large negative value,
*     for which the exponential function should underflow. The present
*     code tests for this condition and explicitly sets those rate
*     coefficients to 0.
*
*-----------------------------------------------------------------------
*.declarations

*   ..local variables
      integer ncall, n2, iCv, is
      real (kind=R8) :: t_av
      real (kind=R8) ::
     *   wrk0(nCv), expu
      save ncall
*   ..procedures
      intrinsic log, exp, max, dim
      external subini, subend, xertst, expu
      external b2xvsg
*   ..initialisation
      data ncall/0/

*-----------------------------------------------------------------------
*.computation

* ..preliminaries
*   ..subprogram start-up calls
      call subini ('b2sqcx')
*   ..test nCv, ns
      call xertst (0.le.nCv, 'faulty argument nCv')
      call xertst (1.le.ns, 'faulty argument ns')
*   ..test ev, amh0
      call xertst (0.0_R8.lt.ev, 'faulty argument ev')
      call xertst (1.0_R8.le.amh0.and.amh0.le.3.0_R8,
     &  'faulty argument range amh0')
*   ..extensive tests on first few calls
      if (ncall.lt.3) then
*    ..test th0
       call b2xvsg (nCv, th0, 1, 'th0', '.gt.')
      endif
*   ..compute wrk0(,) = log(th0(,)/(amh0*ev))
      do iCv = 1, nCv
        t_av = 0.5_R8*(th0(iCv)+thn0(iCv))
        wrk0(iCv) = log(t_av/(amh0*ev))
      enddo

* ..main computation
*   ..compute rate coefficients
      do is = 0, ns-1
       if (switch%b2sqcx_styl0.eq.0.or.is.lt.2) then
        do iCv = 1, nCv
*      ..compute charge exchange rate coefficient
          rcx(iCv,is) = switch%b2sqcx_phm0*
     &     expu(rlcx0(iCv,0,is)+rlcx0(iCv,1,is)*wrk0(iCv))
        enddo
       else
*      ..compute charge exchange rate coefficient
          rcx(:,is) = 0.0_R8
       endif
      enddo

* ..return
      ncall = ncall+1
      call subend ()
      return

*-----------------------------------------------------------------------
*.end b2sqcx

      end subroutine b2sqcx
