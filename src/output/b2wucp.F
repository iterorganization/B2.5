*-----------------------------------------------------------------------
*.specification

      subroutine b2wucp (nput, version, ns)
      use b2mod_types
      use b2mod_constants
      use b2mod_b2cmpa
      use b2mod_b2cmpb
      use b2mod_b2cmpt
      use b2mod_subsys
      implicit none
*   ..input arguments (unchanged on exit)
      integer nput, ns
      character*10 version
*   ..output arguments (unspecified on entry)
*     (none)
*   ..common blocks
*-----------------------------------------------------------------------
*.documentation
*
*  1. purpose
*
*     B2WUCP writes to an un*formatted file the parameters in the
*     physics common blocks.
*
*
*  3. description (see also routine b2cdca)
*
*     This routine writes to unit (nput) the parameters in the
*     physics common blocks: /b2cmpa/, /b2cmpb/ and /b2cmpt/.
*     Unit nput must be connected to an un*formatted file.
*
*
*  5. parameters (see also routine b2cdcv)
*
*     nput - integer, input.
*     nput specifies a unit number for un*formatted output.
*
*     ns - integer, input.
*     ns specifies the number of atomic species in the calculation.
*     The species are indexed by (0:ns-1).
*     It will hold that 1.le.ns.
*
*     (parameters in common:)
*
*     /b2cmpa/ - input.
*     /b2cmpa/ specifies basic physics parameters.
*     See the calling routine for description.
*
*     /b2cmpb/ - input.
*     /b2cmpb/ specifies physics parameters for boundary conditions.
*     See the calling routine for description.
*
*     /b2cmpt/ - input.
*     /b2cmpt/ specifies physics parameters for transport.
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
      integer ireg, idum(0:9)
*   ..procedures
      external xertst, cfwuin, cfwure
      external b2xvcp, b2wuzd

*-----------------------------------------------------------------------
*.computation

* ..preliminaries
*   ..subprogram start-up calls
      call subini ('b2wucp')
*   ..test nput
      call xertst (1.le.nput, 'faulty argument nput')
*   ..test ns
      call xertst (1.le.ns, 'faulty argument ns')
*   ..test /b2cmpa/, /b2cmpb/, /b2cmpt/
      call b2xvcp (ns)

* ..save common blocks
*   ..save /b2cmpa/
      call b2wuzd (nput, version, ns, zamin, zamax, zn, am)
*   ..save /b2cmpb/
      if (cbnrso.gt.0.and.cbnrno.gt.0.and.
     &    cbnrwe.gt.0.and.cbnrea.gt.0) then
       idum(0) = cbirso
       idum(1) = cbirno
       idum(2) = cbirwe
       idum(3) = cbirea
       idum(4) = cbnrso
       idum(5) = cbnrno
       idum(6) = cbnrwe
       idum(7) = cbnrea
       call cfwuin (nput, 8, idum, 'cbir,cbnr')
       call cfwure (nput, cbnrso+cbnrno+cbnrwe+cbnrea, cbrbrk, 'cbrbrk')
       do ireg = 0, cbnrso+cbnrno+cbnrwe+cbnrea-1
        call cfwure (nput, 8*ns, cbsna(0,0,ireg), 'cbsna')
        call cfwure (nput, 8*ns, cbsmo(0,0,ireg), 'cbsmo')
        call cfwure (nput, 8*ns, cbshi(0,0,ireg), 'cbshi')
        call cfwure (nput, 8, cbshe(0,ireg), 'cbshe')
        call cfwure (nput, 8, cbsch(0,ireg), 'cbsch')
        call cfwure (nput, 8*ns, cbrec(0,0,ireg), 'cbrec')
        call cfwure (nput, 8, cbmsa(0,ireg), 'cbmsa')
        call cfwure (nput, 8*ns, cbmsb(0,0,ireg), 'cbmsb')
       enddo
      endif
*   ..save /b2cmpt/
      call cfwure (nput, 8*ns, cfdf0, 'cfdf0')
      call cfwure (nput, 8*ns, cfdna, 'cfdna')
      call cfwure (nput, 8*ns, cfdpa, 'cfdpa')
      call cfwure (nput, 8*ns, cfvla, 'cfvla')
      call cfwure (nput, 8*ns, cfvsa, 'cfvsa')
      call cfwure (nput, 8*ns, cfhci, 'cfhci')
      call cfwure (nput, 8, cfhce, 'cfhce')
      call cfwure (nput, 8, cfsig, 'cfsig')
      call cfwure (nput, 8, cfalf, 'cfalf')
      call cfwure (nput, 8, cflim, 'cflim')

* ..return
      call subend ()
      return

*-----------------------------------------------------------------------
*.end b2wucp

      end subroutine b2wucp
