*-----------------------------------------------------------------------
*.specification

      program b2yt
      use b2mod_types
      use b2mod_geo
      use b2mod_geo2
      use b2mod_work
      use b2mod_rates
      use b2mod_plasma
      use b2mod_b2plot
      use b2mod_version
      use b2mod_indirect
      use b2mod_xerset
      implicit none
#include "DIMENSIONS.F"

*-----------------------------------------------------------------------
*.documentation
*
*  1. purpose
*
*     B2YT (main program) is a postprocessor for B2. It transfers
*     the plasma state to a mesh of different size and/or a different
*     species bundling. It also attempts
*     to re-write a new set of input files consistent with the new
*     mesh size and species bundle distribution. This program calls
*     the driver routine b2ytdr for the actual computation.
*
*
*  2. specification
*
*     Main program.
*
*
*  3. description (see also routine b2cdca)
*
*     The complete program b2yt transfers the result of a B2
*     calculation to a mesh of different size.
*     This program unit opens and closes the input/output units, and
*     may perform some other system-dependent operations. It calls
*     routine b2ytdr to control the actual computation.
*
*     The input units are:
*     ninp(0): formatted; specifies the output mesh size.
*     ninp(1): un*formatted; provides the plasma state.
*     ninp(2): formatted; provides the run switches.
*     ninp(3): un*formatted; provides the (new) geometry.
*     ninp(4): un*formatted; provides the (new) atomic rates.
*     ninp(5): formatted; provides the transport switches.
*     ninp(6): formatted; provides the initial state.
*     ninp(7): formatted; provides geometric parameters.
*     ninp(8): un*formatted; provides default physics parameters.
*     For non-conformal cases, we additionally have:
*     70: un*formatted; provides the (old) geometry.
*
*     The output units are:
*     nout(0): formatted; provides printed output.
*     nout(1): un*formatted; provides the transferred plasma state.
*     nout(2): formatted; provides modified run switches.
*     nout(3): formatted; provides modified transport parameters.
*     nout(4): formatted; provides modified initial state.
*     nout(5): formatted; provides list of modified files.
*     nout(6): formatted; provides modified (inverted field) geometry.
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
      integer ninp(0:8), nout(0:6), nx1, ny1, idum(0:9)
      integer io
      logical file_open
*   ..procedures
      external prgini, prgfin, xertst, cfopen, cfruin
      external b2ytdr, b2ytx0, ipgeti
*   ..initialise input/output units
      data ninp/50,51,52,53,54,55,56,57,58/
      data nout/60,61,62,63,64,65,66/

*-----------------------------------------------------------------------
*.documentation-internal
*
*     The following common blocks have their outermost declaration in
*     this routine; they need not be preserved between calls.
*
*     Description of some local variables:
*
*     ninp - (0:8) integer array.
*     ninp specifies the input unit numbers.
*
*     nout - (0:5) integer array.
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
*     nx1, ny1 - integers.
*     nx1 and ny1 specify the mesh dimensions for the output state.
*     It will hold that 0.le.nx1 and 0.le.ny1.
*
*-----------------------------------------------------------------------
*.computation

*   ..program start-up calls
      call prgini ('b2yt')
*   ..open files
      call cfopen (ninp(0),'b2yt.dat','old','formatted')
      call cfopen (ninp(1),'b2fstate','old','un*formatted')
      call cfopen (ninp(2),'b2mn.dat','old','formatted')
      call cfopen (ninp(3),'b2fgmtry','old','un*formatted')
      call cfopen (ninp(4),'b2frates','old','un*formatted')
      call cfopen (ninp(5),'b2ah.dat','old','formatted')
      call cfopen (ninp(6),'b2ai.dat','old','formatted')
      call cfopen (ninp(7),'b2ag.dat','old','formatted')
      call cfopen (ninp(8),'b2fpardf','old','un*formatted')
      call cfopen (nout(0),'b2yt.prt','new','formatted')
      call cfopen (nout(1),'b2fstatt','new','un*formatted')
      call cfopen (nout(2),'b2mn.dat.2','new','formatted')
      call cfopen (nout(5),'b2yt.files','new','formatted')
*   ..mark output unit for error messages
      call xerset (0)
*   ..obtain version numbers                                             !xpb
      call cfverr (ninp(1), b2fstate_version)                            !xpb
      call cfverr (ninp(3), b2fgmtry_version)                            !xpb
      call cfverr (ninp(4), b2frates_version)                            !xpb
      call cfverr (ninp(8), b2fparam_version)                            !xpb
*   ..obtain nx, ny, ns from old state file
      call cfruin (ninp(1), 3, idum, 'nx,ny,ns')
      nx = idum(0)
      ny = idum(1)
      ns = idum(2)
      call xertst (0.le.nx.and.0.le.ny.and.1.le.ns,
     &  'faulty input nx, ny, ns from plasma state file')
*   ..obtain nx1, ny1
      call b2ytx0 (ninp(0), nx1, ny1)
      call xertst (0.le.nx1.and.0.le.ny1,
     &  'faulty input nx1, ny1')
*   ..obtain nx, ny from geometry file
      call cfruin (ninp(3), 2, idum, 'nx,ny')
      call xertst (nx1.eq.idum(0), 
     & 'inconsistent nx definition between b2fgmtry and b2yt.dat!')
      call xertst (ny1.eq.idum(1),
     & 'inconsistent ny definition between b2fgmtry and b2yt.dat!')
*   ..read and check ns from parameters file
      call cfruin (ninp(8), 1, idum, 'ns')
*   ..populate the first elements of b2yt.files
      write(nout(5),'(a)') 'b2fstatt'
      write(nout(5),'(a)') 'b2mn.dat.2'
*   ..call principal driver routine
      call alloc_b2mod_plasma(nx,ny,ns)
      call alloc_b2mod_geo(nx,ny)
      call alloc_b2mod_geo2(nx1,ny1,DEF_NCUT)
      call alloc_b2mod_work(nx,ny,ns)
      call alloc_b2mod_rates(nx,ny,ns)
      call alloc_b2mod_indirect(nx,ny,nncutmax)
      call b2ytdr (ninp, nout, nx, ny, nx1, ny1)
*   ..close files
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
      call dealloc_b2mod_geo2
      call dealloc_b2mod_work
      call dealloc_b2mod_rates
      call dealloc_b2mod_plasma
      call dealloc_b2mod_indirect
      call prgfin ()
      stop 'b2yt'

*-----------------------------------------------------------------------
*.end b2yt

      end program b2yt
