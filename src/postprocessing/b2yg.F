*-----------------------------------------------------------------------
*.specification

      program b2yg
      use b2mod_types
      use b2mod_geo
      use b2mod_indirect
      use b2mod_xerset
      implicit none

*-----------------------------------------------------------------------
*.documentation
*
*  1. purpose
*
*     B2YG (main program) is a postprocessor for b2. It reads and
*     displays the geometry and magnetic field. This program calls the
*     driver routine b2ygdr for the actual computation.
*
*
*  2. specification
*
*     Main program.
*
*
*  3. description (see also routine b2cdca)
*
*     The complete program b2yg displays a b2 geometry.
*     This program unit opens and closes the input/output units, and
*     may perform some other system-dependent operations. It calls
*     routine b2ygdr to control the actual computation.
*
*     The input units are:
*     ninp(0): formatted; provides output control.
*     ninp(1): un*formatted; full specification of the geometry.
*
*     The output units are:
*     nout(0): formatted; provides summary printed output.
*     nout(1): graphical; displays the geometry and the flux functions.
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
      integer ninp(0:1), nout(0:1), nx, ny, idum(0:9), io
      logical file_open
*   ..procedures
      external prgini, prgend, xertst, cfopen, gnopen, cfruin,
     &  gopks, gclks, gclwk
      external b2ygdr
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
*     nx, ny - integer.
*     nx and ny specify the number of interior cells along the first
*     and the second coordinate, respectively. The total number of
*     cells is (nx+2)*(ny+2); they are indexed by (-1:nx,-1:ny).
*     It will hold that 0.le.nx and 0.le.ny.
*
*-----------------------------------------------------------------------
*.computation

*   ..program start-up calls
      call prgini ('b2yg')
*   ..open gks
      call gopks (0, 0)
*   ..open files
      call cfopen (ninp(0),'b2yg.dat','old','formatted')
      call cfopen (ninp(1),'b2fgmtry','old','un*formatted')
      call cfopen (nout(0),'b2yg.prt','new','formatted')
      call gnopen (nout(1),'b2yg.plt')
*   ..mark output unit for error messages
      call xerset (0)
*   ..obtain version number						 !xpb
      call cfverr (ninp(1), b2fgmtry_version)				 !xpb
*   ..obtain nx, ny
      call cfruin (ninp(1), 2, idum, 'nx,ny')
      nx = idum(0)
      ny = idum(1)
      call xertst (0.le.nx.and.0.le.ny, 'faulty input nx, ny')
*   ..call principal driver routine
      call alloc_b2mod_geo(nx,ny)
      call alloc_b2mod_indirect(nx,ny,nncutmax)                 	 !xpb
      call b2ygdr (ninp, nout, nx, ny)
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
      call dealloc_b2mod_indirect
      call prgend ()
      stop 'b2yg'

*-----------------------------------------------------------------------
*.end b2yg

      end program b2yg
