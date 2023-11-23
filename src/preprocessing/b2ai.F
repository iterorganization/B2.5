*-----------------------------------------------------------------------
*.specification

      program b2ai
      use b2mod_types
      use b2us_geo
      use b2us_map
      use b2us_plasma
      use b2mod_xerset
      implicit none

*-----------------------------------------------------------------------
*.documentation
*
*  1. purpose
*
*     B2AI (main program) prepares an initial plasma state for B2.
*     This program calls the driver routine b2aidr for the actual
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
*     The complete program b2ai prepares an initial state for a B2
*     calculation.
*     This program unit opens and closes the input/output units, and
*     may perform some other system-dependent operations. It calls
*     routine b2aidr to control the actual computation.
*
*     The input units are:
*     ninp(0): formatted; provides parameters describing the state.
*     ninp(1): formatted; provides parameters describing the geometry.
*     ninp(2): un*formatted; provides the geometry.
*
*     The output units are:
*     nout(0): formatted; provides summary printed output.
*     nout(1): un*formatted; provides the initial state.
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
*  8. further comments
*
*-----------------------------------------------------------------------
*.declarations

*   ..common blocks

*   ..local variables
      integer ninp(0:2), nout(0:1), nCv, nFc, ns, io                    !srv 22.12.09
      logical file_open
      type(mapping)  :: mpg
      type(geometry) :: geo
      type(B2State)  :: st
*   ..procedures
      external prgini, prgend, cfopen
      external b2aix0, b2aidr
*   ..initialise input/output units
      data ninp/50,51,52/                                               !srv 22.12.09
      data nout/60,61/

*-----------------------------------------------------------------------
*.documentation-internal
*
*     The following common blocks have their outermost declaration in
*     this routine; they need not be preserved between calls.
*
*     Description of some local variables:
*
*     ninp - (0:2) integer array.
*     ninp specifies the input unit numbers.
*
*     nout - (0:1) integer array.
*     nout specifies the output unit numbers.
*
*     nCv, nFc - integer.
*     nCv and nFc specify the number of control volumes and cell faces
*     in the unstructured mesh.
*     It will hold that 0.le.nCv and 0.le.nFc.
*
*     ns - integer.
*     ns specifies the number of atomic species in the calculation.
*     The species are indexed by (0:ns-1).
*     It will hold that 1.le.ns.
*
*-----------------------------------------------------------------------
*.computation

*   ..program start-up calls
      call prgini ('b2ai')
*   ..open files
      call cfopen (ninp(0),'b2ai.dat','old','formatted')
      call cfopen (ninp(1),'b2fgmtry','old','un*formatted')
      call cfopen (nout(0),'b2ai.prt','new','formatted')
      call cfopen (nout(1),'b2fstati','new','un*formatted')
*   ..mark output unit for error messages
      call xerset (0)
*   ..obtain nCv, nFc, ns
      call b2aix0 (ninp, nCv, nFc, ns, mpg)
*   ..allocate structures
      call alloc_mapping(mpg)
      call alloc_geometry(geo,mpg%nCv,mpg%nFc,mpg%nVx,mpg%nFs,
     .  mpg%nFt,mpg%nVmxCv)
      call createB2State (nCv, nFc, ns, st)
*   ..call principal driver routine
      call b2aidr (ninp, nout, nCv, nFc, ns, geo, mpg, st)
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
      call destroyB2State (st)
      call dealloc_geometry (geo)
      call dealloc_mapping (mpg)
      call prgend ()
      stop 'b2ai'

*-----------------------------------------------------------------------
*.end b2ai

      end program b2ai
