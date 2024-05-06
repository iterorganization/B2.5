*-----------------------------------------------------------------------
*.specification

      program b2ah
      use b2mod_types
      use b2mod_version
      use b2mod_xerset
      implicit none

*-----------------------------------------------------------------------
*.documentation
*
*  1. purpose
*
*     B2AH (main program) prepares default physics parameters for B2.
*     This program calls the driver routine b2ahdr for the actual
*     computation.
*
*
*  2. specification
*
*     Main program.
*
*
*  3. description (see also routine b2cdca)
*
*     The complete program b2ah prepares a file of default physics
*     parameters for a B2 calculation.
*     This program unit opens and closes the input/output units, and
*     may perform some other system-dependent operations. It calls
*     routine b2ahdr to control the actual computation.
*
*     The input units are:
*     ninp(0): formatted; provides physics parameters.
*
*     The output units are:
*     nout(0): formatted; provides summary printed output.
*     nout(1): un*formatted; provides default physics parameters.
*
*     (See routine b2cdca for the meaning of 'un*formatted'.)
*
*
*  5. parameters (see also routine b2cdcv)
*
*     None.
*
*
*  6. error indicators
*
*     In case an error condition is detected, a call is made to the
*     routine xerrab. This causes an error message to be printed,
*     after which the program halts.
*
*-----------------------------------------------------------------------
*.declarations

*   ..common blocks

*   ..local variables
      integer ninp(0:0), nout(0:1), ns, io
      logical file_open
*   ..procedures
      external prgini, prgend, cfopen, cfruin, cfruch
      external b2ahx0, b2ahdr
*   ..initialise input/output units
      data ninp/50/
      data nout/60,61/

*-----------------------------------------------------------------------
*.documentation-internal
*
*     The following common blocks have their outermost declaration in
*     this routine; they need not be preserved between calls.
*
*     Description of some local variables:
*
*     ninp - (0:0) integer array.
*     ninp specifies the input unit numbers.
*
*     nout - (0:1) integer array.
*     nout specifies the output unit numbers.
*
*     ns - integer.
*     ns specifies the number of atomic species in the calculation.
*     The species are indexed by (0:ns-1).
*     It will hold that 1.le.ns.
*
*-----------------------------------------------------------------------
*.computation

*   ..program start-up calls
      call prgini ('b2ah')
*   ..open files
      call cfopen (ninp(0),'b2ah.dat','old','formatted')
      call cfopen (nout(0),'b2ah.prt','new','formatted')
      call cfopen (nout(1),'b2fpardf','new','un*formatted')
*   ..mark output unit for error messages
      call xerset (0)
*   ..obtain ns
      call b2ahx0 (ninp(0), ns)
*   ..call principal driver routine
      call b2ahdr (ninp, nout, ns)
*   ..close graphics files
*     (none)
*   ..close other files
      do io=0,size(ninp)-1
         inquire(unit=ninp(io),opened=file_open)
         if (file_open) close(ninp(io))
      enddo
      do io=0,size(nout)-1
         inquire(unit=nout(io),opened=file_open)
         if (file_open) close(nout(io))
      enddo

* ..end of computation
      call prgend ()
      stop 'b2ah'

*-----------------------------------------------------------------------
*.end b2ah

      end program b2ah
