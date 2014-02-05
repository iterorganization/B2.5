*-----------------------------------------------------------------------
*.specification

      program b2yp
      use b2mod_types
      use b2mod_rates
      use b2mod_version
      use b2mod_geo
      use b2mod_plasma
      use b2mod_indirect
      use b2mod_ppout
      implicit none

*-----------------------------------------------------------------------
*.documentation
*
*  1. purpose
*
*     B2YP (main program) is a post-processor for b2. It displays the
*     plasma state. This program calls the driver routine b2ypdr for
*     the actual computation.
*
*
*  2. specification
*
*     Main program.
*
*
*  3. description (see also routine b2cdca)
*
*     The complete program b2yp performs post-processing of the
*     result of a b2 calculation.
*     This program unit opens and closes the input/output units, and
*     may perform some other system-dependent operations. It calls
*     routine b2ypdr to control the actual computation.
*
*     The input units are:
*     ninp(0): formatted; provides output control parameters.
*     ninp(1): un*formatted; provides the geometry.
*     ninp(2): un*formatted; provides the physics parameters.
*     ninp(3): un*formatted; provides the plasma state.
*     ninp(4): un*formatted; provides atomic rate data.
*
*     The output units are:
*     nout(0): formatted; provides printed output.
*     nout(1): graphical; provides graphical output.
*     nout(2): formatted; may provide special purpose output.
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
      integer ninp(0:4), nout(0:2), nx, ny, ns, idum(0:9)
*   ..procedures
      external prgini, prgend, xerset, xertst, cfopen, gnopen, cfruin,
     &  gopks, gclks, gclwk
      external b2ypdr
*   ..initialise input/output units
      data ninp/50,51,52,53,54/
      data nout/60,61,62/

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
*     nout - (0:2) integer array.
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
      call prgini ('b2yp')
*   ..open gks
      call gopks (0, 0)
*   ..open files
      call cfopen (ninp(0),'b2yp.dat','old','formatted')
      call cfopen (ninp(1),'b2fgmtry','old','un*formatted')
      call cfopen (ninp(2),'b2fparam','old','un*formatted')
      call cfopen (ninp(3),'b2fstate','old','un*formatted')
      call cfopen (ninp(4),'b2frates','old','un*formatted')
      call cfopen (nout(0),'b2yp.prt','new','formatted')
      call gnopen (nout(1),'b2yp.plt')
      call cfopen (nout(2),'b2specp','new','formatted')
*   ..mark output unit for error messages
      call xerset (0)
*   ..obtain version numbers						 !xpb
      call cfverr (ninp(1), b2fgmtry_version)				 !xpb
      call cfverr (ninp(2), b2fparam_version)				 !xpb
      call cfverr (ninp(3), b2fstate_version)				 !xpb
      call cfverr (ninp(4), b2frates_version)				 !xpb
*   ..obtain nx, ny, ns
      call cfruin (ninp(3), 3, idum, 'nx,ny,ns')
      nx = idum(0)
      ny = idum(1)
      ns = idum(2)
      call xertst (0.le.nx.and.0.le.ny.and.1.le.ns,
     &  'faulty input nx, ny, ns from plasma state file')
      call cfruin (ninp(1), 2, idum(0), 'nx,ny')
      call cfruin (ninp(2), 1, idum(2), 'ns')
      call xertst (idum(0).eq.nx.and.idum(1).eq.ny,
     &  'faulty input nx, ny from geometry file')
      call xertst (idum(2).eq.ns,
     &  'faulty input ns from physics parameters file')
*   ..call principal driver routine
      call alloc_b2mod_geo(nx,ny)
      call alloc_b2mod_rates(nx,ny,ns)
      call alloc_b2mod_plasma(nx,ny,ns)
      call alloc_b2mod_indirect(nx,ny,nncutmax)                          !xpb
      call alloc_b2mod_ppout(nx,ny,ns)
      call b2ypdr (ninp, nout, nx, ny, ns)
*   ..close graphics files
      call gclwk (nout(1))
*   ..close gks
      call gclks ()

* ..end of computation
      call prgend ()
      stop 'b2yp'

*-----------------------------------------------------------------------
*.end b2yp

      end program b2yp
