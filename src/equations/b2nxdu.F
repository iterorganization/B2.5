*-----------------------------------------------------------------------
*.specification

      subroutine b2nxdu (nCv, ns, na, ua, smq, smqdu, wrk0, wrk1)
      use b2mod_types
      use b2mod_constants
      use b2mod_b2cmpa
! csc The following are not necessary for computation but are needed
!     for adjoint AD to avoid side-effect variables
      use b2mod_ad
     & , only : ncall_b2nxdu
      use b2mod_subsys
      implicit none
*   ..input arguments (unchanged on exit)
      integer nCv, ns
      real (kind=R8) ::
     *   na(nCv,0:ns-1), ua(nCv,0:ns-1),
     &  smq(nCv,0:3,0:ns-1)
*   ..output arguments (unspecified on entry)
      real (kind=R8) ::
     *   smqdu(nCv)
*   ..workspace arguments (unspecified on entry and on exit)
      real (kind=R8) ::
     *   wrk0(nCv), wrk1(nCv)
*   ..common blocks
*-----------------------------------------------------------------------
*.documentation
*
*  1. purpose
*
*     B2NXDU computes the artificial velocity shift that is
*     required to make the total numerical friction force cancel.
*
*-----------------------------------------------------------------------
*.declarations

*   ..local variables
      integer iCv, is
*   ..procedures
      external xertst, sfill
      external b2xvsg, b2xvfx
*   ..initialisation

*-----------------------------------------------------------------------
*.computation

* ..preliminaries
*   ..subprogram start-up calls
      call subini ('b2nxdu')
*   ..test nCv
      call xertst (0.lt.nCv, 'faulty argument nCv')
*   ..extensive tests on first few calls
      if (ncall_b2nxdu.lt.3) then
*    ..test sign of smq
       do is = 0, ns-1
        call b2xvsg (nCv, smq(1,1,is), 1, 'smq1', '.le.')
        call b2xvsg (nCv, smq(1,3,is), 1, 'smq3', '.le.')
       enddo
      endif

* ..compute smqdu
*   ..compute auxiliaries wrk0 and wrk1
      call sfill (nCv, 0.0_R8, wrk0, 1)
      call sfill (nCv, 0.0_R8, wrk1, 1)
      do is = 0, ns-1
       wrk0 = wrk0 + smq(:,1,is) + smq(:,3,is)*am(is)*mp*na(:,is)
       wrk1 = wrk1 + smq(:,0,is) +
     &               smq(:,1,is)*ua(:,is) +
     &               smq(:,2,is)*am(is)*mp*na(:,is) +
     &               smq(:,3,is)*am(is)*mp*na(:,is)*ua(:,is)
      enddo
*   ..compute smqdu
      do iCv = 1, nCv
        if (wrk0(iCv).eq.0.0_R8) then
         smqdu(iCv) = 0.0_R8
        else
         smqdu(iCv) = wrk1(iCv)/wrk0(iCv)
        endif
      enddo

* ..return
      ncall_b2nxdu = ncall_b2nxdu+1
      call subend ()
      return

*-----------------------------------------------------------------------
*.end b2nxdu

      end subroutine b2nxdu
