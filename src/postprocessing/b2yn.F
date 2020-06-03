*-----------------------------------------------------------------------
*.specification

      program b2yn
      use b2mod_types
      use b2mod_rates
      use b2mod_version
      use b2mod_geo
      use b2mod_plasma
      use b2mod_rates
      use b2mod_residuals
      use b2mod_sources
      use b2mod_transport
      use b2mod_work
      use b2mod_ppout
      implicit none

*-----------------------------------------------------------------------
*.documentation
*
*  1. purpose
*
*     INCLUDE 'doc/b2yn.purp'
*
*
*  2. specification
*
*     Main program.
*
*
*  3. description (see also routine b2cdca)
*
*     The complete program b2yn displays the progress of the inner
*     B2 iterations.
*     This program unit opens and closes the input/output units, and
*     may perform some other system-dependent operations. It calls
*     routine b2yndr to control the actual computation.
*
*     The input units are:
*     ninp(0): formatted; provides output control parameters.
*     ninp(1): un*formatted; provides the corrections and residuals.
*     ninp(2): un*formatted; provides atomic rate data.
*     ninp(3): un*formatted; provides the plasma state.
*     ninp(4): formatted; provides run switches.
*
*     The output units are:
*     nout(0): formatted; provides printed output.
*     nout(1): graphical; provides graphical output.
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
      integer ninp(0:4), nout(0:1), nx, ny, ns, idum(0:9), io
      logical file_open
*   ..procedures
      external prgini, prgfin, xerset, xertst, cfopen, gnopen, cfruin,
     &  gopks, gclks, gclwk
      external b2yndr
*   ..initialise input/output units
      data ninp/50,51,52,53,54/
      data nout/60,61/

*-----------------------------------------------------------------------
*.documentation-internal
*
*     The following common blocks have their outermost declaration in
*     this routine; they need not be preserved between calls.
*
*     Description of some local variables:
*
*     ninp - (0:4) integer array.
*     ninp specifies the input unit numbers.
*
*     nout - (0:1) integer array.
*     nout specifies the output unit numbers.
*
*     nx, ny - integers.
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
      call prgini ('b2yn')
*   ..open gks
      call gopks (0, 0)
*   ..open files
      call cfopen (ninp(0),'b2yn.dat','old','formatted')
      call cfopen (ninp(1),'b2ftrack','old','un*formatted')
      call cfopen (ninp(2),'b2frates','old','un*formatted')
      call cfopen (ninp(3),'b2fstate','old','un*formatted')
      call cfopen (ninp(4),'b2mn.dat','old','formatted')
      call cfopen (nout(0),'b2yn.prt','new','formatted')
      call gnopen (nout(1),'b2yn.plt')
*   ..mark output unit for error messages
      call xerset (0)
*   ..obtain version numbers    					 !xpb
      call cfverr (ninp(1), b2ftrack_version)				 !xpb
      call cfverr (ninp(2), b2frates_version)				 !xpb
      call cfverr (ninp(3), b2fstate_version)				 !xpb
*   ..obtain nx, ny, ns
      call cfruin (ninp(1), 3, idum, 'nx,ny,ns')
      nx = idum(0)
      ny = idum(1)
      ns = idum(2)
      call xertst (0.le.nx.and.0.le.ny.and.1.le.ns,
     &  'faulty input nx, ny, ns')
*   ..obtain nx, ny, ns from plasma state file
      call cfruin (ninp(3), 3, idum, 'nx,ny,ns')
      call xertst (
     & nx.eq.idum(0) .and. ny.eq.idum(1) .and. ns.eq.idum(2),
     & 'mistmatched inputs nx, ny, ns')
*   ..call principal driver routine
      call alloc_b2mod_geo(nx,ny)
      call alloc_b2mod_work(nx,ny,ns)
      call alloc_b2mod_rates(nx,ny,ns)
      call alloc_b2mod_plasma(nx,ny,ns)
      call alloc_b2mod_ppout(nx,ny,ns)
      call b2yndr (ninp, nout, nx, ny, ns)
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
      call dealloc_b2mod_geo
      call dealloc_b2mod_work
      call dealloc_b2mod_rates
      call dealloc_b2mod_ppout
      call dealloc_b2mod_plasma
      call prgfin ()
      stop 'b2yn'

*-----------------------------------------------------------------------
*.end b2yn

      end program b2yn
