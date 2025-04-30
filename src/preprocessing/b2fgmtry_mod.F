*-----------------------------------------------------------------------
*.specification

      program b2fgmtry_mod
      use b2mod_types
      use b2mod_xerset
      implicit none

*-----------------------------------------------------------------------
*.documentation
*
*  1. purpose
*
*     B2FGMTRY_MOD (main program) modifies a geometry saved in b2fgmtry
*     and creates a new file b2fgmtry_mod for B2.
*
*
*  2. specification
*
*     Main program.
*
*
*  3. description (see also routine b2cdca)
*
*     ... to be written...
*
*     The input units are:
*     ninp(0): formatted; provides geometric parameters.
*     ninp(1): un*formatted; old specification of the geometry
*
*     The output units are:
*     nout(0): formatted; provides summary printed output.
*     nout(1): un*formatted; new specification of the geometry.
*
*     (See routine b2cdca for the meaning of 'un*formatted'.)
*
*
*  4. references
*
*     (Refs. /1/ and /2/ are out of date; !! to be replaced.)
*     /1/ B. J. Braams, Computational Studies in Tokamak Equilibrium
*     and Transport, thesis, State University at Utrecht, June 1986.
*     /2/ B. J. Braams, A Multi-Fluid Code for Simulation of the Edge
*     Plasma in Tokamaks, NET Report Nr 68, January 1987.
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
*     The computation performed in this program depends on the
*     application, and it is difficult to provide a useful general
*     statement concerning the accuracy of the calculation.
*     Numerical differentiation is employed to compute the metric
*     parameters of the grid and to compute the magnetic field. So as
*     to retain sufficient accuracy when performing this numerical
*     differentiation, cells should not be too thin. Be careful about
*     the boundary cells in particular.
*
*
*  8. further comments
*
*     adapted by s.wiesen@fz-juelich.de February 2014
*
*-----------------------------------------------------------------------
*.declarations

*   ..common blocks

*   ..local variables
      integer ninp(0:1), nout(0:1), io
      logical file_open
*   ..procedures
      external prgini, prgend, cfopen
*   ..initialise input/output units
      data ninp/50,51/
      data nout/60,61/

*-----------------------------------------------------------------------
*.documentation-internal
*
*     The following common blocks have their outermost declaration in
*     this routine; they need not be preserved between calls.
*
*     Description of some local variables:
*
*     ninp - (0:1) integer array.
*     ninp specifies the input unit numbers.
*
*     nout - (0:1) integer array.
*     nout specifies the output unit numbers.
*
*-----------------------------------------------------------------------
*.computation

*   ..program start-up calls
      call prgini ('b2fgmtry_mod')
*   ..open files
      call cfopen (ninp(0),'b2fgmtry_mod.dat','old','formatted')
      call cfopen (ninp(1),'b2fgmtry','old','un*formatted')
      call cfopen (nout(0),'b2fgmtry_mod.prt','new','formatted')
      call cfopen (nout(1),'b2fgmtry_mod','new','un*formatted')
*   ..mark output unit for error messages
      call xerset (0)
*   ..call principal driver routine
      call b2fgmtry_moddr (ninp, nout)
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
      stop 'b2fgmtry_mod'

*-----------------------------------------------------------------------
*.end b2fgmtry_mod

      end program b2fgmtry_mod
