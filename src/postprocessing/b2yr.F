*-----------------------------------------------------------------------
*.specification

      program b2yr
      use b2mod_types
      use b2mod_geo
      use b2mod_plasma
      use b2mod_rates
      use b2mod_residuals
      use b2mod_sources
      use b2mod_transport
      use b2mod_work
      use b2mod_xerset
      implicit none

*-----------------------------------------------------------------------
*.documentation
*
*  1. purpose
*
*     B2YR (main program) is a postprocessor for b2. It reads and
*     displays a table of atomic rate coefficients. This program calls
*     the driver routine b2yrdr for the actual computation.
*
*
*  2. specification
*
*     Main program.
*
*
*  3. description (see also routine b2cdca)
*
*     This program displays a table of atomic rate coefficients that
*     has been prepared for the b2 code.
*
*     This program unit opens and closes the input/output units, and
*     may perform some other system-dependent operations. It calls
*     routine b2yrdr to control the actual computation.
*
*     The input units are:
*     ninp(0): formatted; provides output control.
*     ninp(1): un*formatted; provides the table of rate coefficients.
*     ninp(2): un*formatted; provides default physics parameters.
*     ninp(3): un*formatted; provides the geometry.
*
*     The output units are:
*     nout(0): formatted; provides printed summary of rate coefs.
*     nout(1): graphical: provides graphical summary of rate coefs.
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
      integer ninp(0:3), nout(0:1), io
      logical file_open
*   ..procedures
      external prgini, prgend, cfopen, gnopen, gopks, gclks, gclwk
      external b2yrdr
*   ..initialise input/output unit numbers
      data ninp/50,51,52,53/
      data nout/60,61/

*-----------------------------------------------------------------------
*.documentation-internal
*
*     The following common blocks have their outermost declaration in
*     this routine; they need not be preserved between calls.
*
*     Description of some local variables:
*
*     ninp - (0:3) integer array.
*     ninp specifies the input unit numbers.
*
*     nout - (0:1) integer array.
*     nout specifies the output unit numbers.
*
*-----------------------------------------------------------------------
*.computation

*   ..program start-up calls
      call prgini ('b2yr')
*   ..open gks
      call gopks (0, 0)
*   ..open files
      call cfopen (ninp(0),'b2yr.dat','old','formatted')
      call cfopen (ninp(1),'b2frates','old','un*formatted')
      call cfopen (ninp(2),'b2fpardf','old','un*formatted')
      call cfopen (ninp(3),'b2fgmtry','old','un*formatted')
      call cfopen (nout(0),'b2yr.prt','new','formatted')
      call gnopen (nout(1),'b2yr.plt')
*   ..mark output unit for error messages
      call xerset (0)
*   ..display table
      call b2yrdr (ninp, nout)
*   ..close graphics files
      call gclwk (nout(1))
*   ..close gks
      call gclks ()
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
      stop 'b2yr'

*-----------------------------------------------------------------------
*.end b2yr

      end program b2yr
