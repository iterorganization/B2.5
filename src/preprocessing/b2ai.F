*-----------------------------------------------------------------------
*.specification

      program b2ai
      use b2mod_types
      use b2mod_rates
      use b2mod_geo
      use b2mod_work
      use b2mod_plasma
      use b2mod_indirect
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
      integer ninp(0:3), nout(0:1), nx, ny, ns, io                   !srv 22.12.09
      logical file_open
*   ..procedures
      external prgini, prgend, cfopen
      external b2aix0, b2aidr
*   ..initialise input/output units
      data ninp/50,51,52,53/                                         !srv 22.12.09
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
*     nx, ny - integer.
*     nx and ny specify the number of interior cells along the first
*     and the second coordinate, respectively. The total number of
*     cells is (nx+2)*(ny+2); they are indexed by (-1:nx,-1:ny).
*     It will hold that 0.le.nx and 0.le.ny.
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
      call cfopen (ninp(0),'b2ag.dat','old','formatted')
      call cfopen (ninp(1),'b2ai.dat','old','formatted')
      call cfopen (ninp(2),'b2fgmtry_st','old','un*formatted')
      call cfopen (nout(0),'b2ai.prt','new','formatted')
      call cfopen (nout(1),'b2fstati_st','new','un*formatted')
*   ..mark output unit for error messages
      call xerset (0)
*   ..obtain nx, ny, ns
      call b2aix0 (ninp, nx, ny, ns)
*   ..call principal driver routine
      call alloc_b2mod_geo(nx,ny)
      call alloc_b2mod_work(nx,ny,ns)
      call alloc_b2mod_rates(nx,ny,ns)
      call alloc_b2mod_plasma(nx,ny,ns)
      call alloc_b2mod_indirect(nx,ny,nncutmax)
      call b2aidr (ninp, nout, nx, ny, ns)
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
      call dealloc_b2mod_geo
      call dealloc_b2mod_work
      call dealloc_b2mod_rates
      call dealloc_b2mod_plasma
      call dealloc_b2mod_indirect
      call prgend ()
      stop 'b2ai'

*-----------------------------------------------------------------------
*.end b2ai

      end program b2ai
