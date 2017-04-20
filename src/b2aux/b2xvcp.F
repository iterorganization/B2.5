*-----------------------------------------------------------------------
*.specification

      subroutine b2xvcp (nx, ns)
      use b2mod_types
      use b2mod_indirect
      use b2mod_constants
      use b2mod_b2cmpa
      use b2mod_b2cmpb
      use b2mod_b2cmpt
      implicit none
*   ..input arguments (unchanged on exit)
      integer nx, ns
*   ..common blocks

*-----------------------------------------------------------------------
*.documentation
*
*  1. purpose
*
*     B2XVCP tests the validity of parameters in the physics common
*     blocks.
*
*
*  3. description (see also routine b2cdca)
*
*     This routine tests the validity of the parameters in the common
*     blocks /b2cmpa/, /b2cmpb/, and /b2cmpt/.
*
*
*  5. parameters (see also routine b2cdcv)
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
*-----------------------------------------------------------------------
*.declarations

*   ..local variables
      integer i, ic, is, ireg
      logical found
*   ..procedures
      intrinsic nint
      external subini, subend, xertst
      external b2xvsg

*-----------------------------------------------------------------------
*.computation

*   ..subprogram start-up calls
      call subini ('b2xvcp')
*   ..test ns
      call xertst (1.le.ns, 'faulty argument ns')
*   ..test nsdecl
      call xertst (ns.le.nsdecl, 'faulty parameter nsdecl')
*   ..test parameters in /b2cmpa/
      call b2xvsg (ns, zn, 1, 'zn', '.gt.')
      call b2xvsg (ns, am, 1, 'am', '.gt.')
      call b2xvsg (ns, zamin, 1, 'zamin', '.ge.')
      call b2xvsg (ns, zamax, 1, 'zamax', '.ge.')
      do is = 0, ns-1
       call xertst (zn(is).eq.anint(zn(is)),
     &   'faulty argument zn; not integral')
       call xertst (zamin(is).eq.anint(zamin(is)),
     &   'faulty argument zamin; not integral')
       call xertst (zamax(is).eq.anint(zamax(is)),
     &   'faulty argument zamax; not integral')
       call xertst (0.0_R8.le.zamin(is).and.zamin(is).le.zamax(is),
     &   'faulty species range zamin, zamax')
       call xertst (zamax(is).le.zn(is),
     &   'mismatch in species range zamax, zn')
       if (zamin(is).eq.0.0_R8) then
         call xertst (zamax(is).eq.0.0_R8,
     &    'no bundling of neutral species allowed')
       endif
      enddo
*   ..test region boundaries in /b2cmpb/
      call xertst (cbirso.eq.0.and.cbirso+cbnrso.eq.cbirno.and.
     &  cbirno+cbnrno.eq.cbirwe.and.cbirwe+cbnrwe.eq.cbirea.and.
     &  cbirea+cbnrea.le.cbnrmx, 'faulty argument cbir, cbnr')
      do ireg = cbirso, cbirso+cbnrso-2
       call xertst (cbrbrk(ireg).le.cbrbrk(ireg+1),
     &   'faulty argument cbrbrk, south, not monotone')
      enddo
      call xertst (1.0_R8.lt.cbrbrk(cbirso+cbnrso-1),
     &  'faulty argument cbrbrk, south, endpoint')
      do ireg = cbirno, cbirno+cbnrno-2
       call xertst (cbrbrk(ireg).le.cbrbrk(ireg+1),
     &   'faulty argument cbrbrk, north, not monotone')
      enddo
      call xertst (1.0_R8.lt.cbrbrk(cbirno+cbnrno-1),
     &  'faulty argument cbrbrk, north, endpoint')
      do ireg = cbirwe, cbirwe+cbnrwe-2
       call xertst (cbrbrk(ireg).le.cbrbrk(ireg+1),
     &   'faulty argument cbrbrk, west, not monotone')
      enddo
      call xertst (1.0_R8.lt.cbrbrk(cbirwe+cbnrwe-1),
     &  'faulty argument cbrbrk, west, endpoint')
      do ireg = cbirea, cbirea+cbnrea-2
       call xertst (cbrbrk(ireg).le.cbrbrk(ireg+1),
     &   'faulty argument cbrbrk, east, not monotone')
      enddo
      call xertst (1.0_R8.lt.cbrbrk(cbirea+cbnrea-1),
     &  'faulty argument cbrbrk, east, endpoint')
*   ..test boundary cross-consistency in /b2cmpb/
      do ic = 1, nncut
        found = .false.
        ireg = 0
        do while (.not.found .and. ireg.lt.cbnrso)
          found = nint(cbrbrk(ireg)*nx).eq.leftcut(ic)
          ireg = ireg + 1
        end do
        call xertst (found, 
     .   'Mismatch between leftcut and cbrbrk (South)!')
        found = .false.
        ireg = 0
        do while (.not.found .and. ireg.lt.cbnrso)
          found = nint(cbrbrk(ireg)*nx).eq.rightcut(ic)
          ireg = ireg + 1
        end do
        call xertst (found, 
     .   'Mismatch between rightcut and cbrbrk (South)!')
      end do
*   ..test parameters in /b2cmpb/
      do ireg = 0, cbnrso+cbnrno+cbnrwe+cbnrea-1
*    ..test cbsna, cbsmo, cbshi
       do is = 0, ns-1
        call xertst (0.0_R8.le.cbsna(0,is,ireg).and.
     &    cbsna(1,is,ireg).le.0.0_R8.and.cbsna(2,is,ireg).le.0.0_R8.and.
     &    cbsna(3,is,ireg).le.0.0_R8.and.cbsna(4,is,ireg).le.0.0_R8.and.
     &    cbsna(5,is,ireg).le.0.0_R8.and.0.0_R8.le.cbsna(6,is,ireg).and.
     &    cbsna(7,is,ireg).le.0.0_R8, 'faulty argument sign cbsna')
        call xertst (cbsmo(3,is,ireg).le.0.0_R8.and.
     &    cbsmo(4,is,ireg).le.0.0_R8.and.cbsmo(6,is,ireg).le.0.0_R8,
     &    'faulty argument sign cbsmo')
        call xertst (0.0_R8.le.cbshi(0,is,ireg).and.
     &    cbshi(1,is,ireg).le.0.0_R8.and.cbshi(2,is,ireg).le.0.0_R8.and.
     &    cbshi(3,is,ireg).le.0.0_R8.and.cbshi(4,is,ireg).le.0.0_R8.and.
     &    cbshi(5,is,ireg).le.0.0_R8.and.0.0_R8.le.cbshi(6,is,ireg).and.
     &    cbshi(7,is,ireg).le.0.0_R8, 'faulty argument sign cbshi')
       enddo
*    ..test cbshe
       call xertst (0.0_R8.le.cbshe(0,ireg).and.
     &   cbshe(1,ireg).le.0.0_R8.and.cbshe(2,ireg).le.0.0_R8.and.
     &   cbshe(3,ireg).le.0.0_R8.and.cbshe(4,ireg).le.0.0_R8.and.
     &   cbshe(5,ireg).le.0.0_R8.and.0.0_R8.le.cbshe(6,ireg).and.
     &   cbshe(7,ireg).le.0.0_R8, 'faulty argument sign cbshe')
*    ..test cbsch
       call xertst (cbsch(1,ireg).le.0.0_R8.and.
     &   cbsch(3,ireg).le.0.0_R8.and.cbsch(5,ireg).le.0.0_R8.and.
     &   cbsch(6,ireg).ge.0.0_R8.and.cbsch(7,ireg).ge.0.0_R8,
     &   'faulty argument sign cbsch')
*    ..test cbrec, cbmsa, cbmsb
*      (there are no tests for cbrec, cbmsa, cbmsb)
      enddo
*   ..test parameters in /b2cmpt/
      do is = 0, ns-1
       call b2xvsg (3, cfdf0(0,is), 1, 'cfdf0', '.ge.')
       call b2xvsg (3, cfdna(0,is), 1, 'cfdna', '.ge.')
       call b2xvsg (3, cfdpa(0,is), 1, 'cfdpa', '.ge.')
       call b2xvsg (3, cfvsa(0,is), 1, 'cfvsa', '.ge.')
       call b2xvsg (3, cfhci(0,is), 1, 'cfhci', '.ge.')
      enddo
      call b2xvsg (3, cfhce(0), 1, 'cfhce', '.ge.')
      call b2xvsg (3, cfsig(0), 1, 'cfsig', '.ge.')
      do i = 0, 3
       call xertst (cfalf(i)**2.le.cfhce(i)*cfsig(i),
     &   'faulty argument, cfhce*cfsig.lt.cfalf**2')
      enddo
      call b2xvsg (3, cflim(0), 1, 'cflim', '.ge.')

* ..return
      call subend ()
      return

*-----------------------------------------------------------------------
*.end b2xvcp

      end
