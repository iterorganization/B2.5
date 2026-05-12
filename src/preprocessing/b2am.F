*-----------------------------------------------------------------------
*.specification

      program b2am
      use b2mod_types
      use b2mod_xerset
      implicit none

*-----------------------------------------------------------------------
*.documentation
*
*  1. purpose
*
*     B2AM (main program) prepares a table of molecular rate coefficients
*     for B2. This program calls the driver routine b2amdr for the
*     actual computation.
*
*
*  2. specification
*
*     Main program.
*
*
*  3. description (see also routine b2cdca)
*
*     B2.5 is composed of several independent programs that communicate
*     through external files. b2am is one of these programs; it
*     computes and writes to an external file a table of molecular rate
*     coefficients. 
*
*     This program unit opens and closes the input/output units, and
*     may perform some other system-dependent operations. It calls
*     routine b2amdr to control the actual computation.
*
*     The input units are:
*     ninp(0): formatted; provides species parameters.
cc     ninp(1): formatted; provides input data for STRAHL calculation.
cc     ninp(2): formatted; provides input data from Weisheit code.
*     ninp(3): un*formatted; provides default physics parameters.
*     ninp(4): un*formatted; provides the geometry.
*
*     The output units are:
*     nout(0): formatted; provides printed summary of rate coefs.
*     nout(1): un*formatted; provides the table of rate coefficients.
*
*     (See routine b2cdca for the meaning of 'un*formatted'.)
*
*
*  4. references
*
*     !! TODO
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
*
*  7. accuracy
*
*     !! TODO
*
*
*  8. further comments
*
*-----------------------------------------------------------------------
*.declarations

*   ..common blocks

*   ..local variables
      integer ninp, nout(0:1)
      integer io
      logical file_open
*   ..procedures
      external prgini, prgend, cfopen
      external b2amdr
*   ..initialise input/output unit numbers
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
*     ninp - integer .
*     ninp specifies the input unit number.
*
*     nout - (0:1) integer array.
*     nout specifies the output unit numbers.
*
*-----------------------------------------------------------------------
*.computation

*   ..program start-up calls
      call prgini ('b2am')
*   ..open files
      call cfopen (ninp,'b2am.dat','old','formatted')
      call cfopen (nout(0),'b2am.prt','new','formatted')
      call cfopen (nout(1),'b2frates_mol','new','un*formatted')
*   ..mark output unit for error messages
      call xerset (0)
*   ..compute, display and save table
      call b2amdr (ninp, nout)
*   ..close files
      inquire(unit=ninp,opened=file_open)
      if (file_open) close (ninp)
      do io = 0, size(nout)-1
        inquire(unit=nout(io),opened=file_open)
        if (file_open) close(nout(io))
      enddo

* ..end of computation
      call prgend ()
      stop 'b2am'

*-----------------------------------------------------------------------
*.end b2am

      end program b2am
