*-----------------------------------------------------------------------
*.specification

      program b2yi
      use b2mod_types
      use b2mod_rates
      use b2mod_version
      use b2mod_geo
      use b2mod_work
      use b2mod_plasma
      use b2mod_indirect
      use b2mod_ppout
      use b2mod_xerset
      implicit none

*-----------------------------------------------------------------------
*.documentation
*
*  1. purpose
*
*     B2YI (main program) is a postprocessor for B2. It reads and
*     displays the plasma state. This program calls the driver routine
*     b2yidr for the actual computation.
*
*
*  2. specification
*
*     Main program.
*
*
*  3. description (see also routine b2cdca)
*
*     The complete program b2yi performs postprocessing of the
*     result of a B2 calculation.
*     This program unit opens and closes the input/output units, and
*     may perform some other system-dependent operations. It calls
*     routine b2yidr to control the actual computation.
*
*     The input units are:
*     ninp(0): formatted; provides output control parameters.
*     ninp(1): un*formatted; provides the plasma state.
*     ninp(2): un*formatted; provides atomic rate data.
*     ninp(3): un*formatted; provides the geometry.
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
      logical file_open, ex
*   ..procedures
      external prgini, prgfin, xertst, cfopen, gnopen, cfruin,
     &  gopks, gclks, gclwk
      external b2yidr
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
      call prgini ('b2yi')
*   ..open gks
      call gopks (0, 0)
*   ..open files
      call cfopen (ninp(0),'b2yi.dat','old','formatted')
      inquire(file='b2fstate',exist=ex)
      if (ex) then
        call cfopen (ninp(1),'b2fstate','old','un*formatted')
      else
        call cfopen (ninp(1),'b2fstati','old','un*formatted')
      endif
      call cfopen (ninp(2),'b2frates','old','un*formatted')
      call cfopen (ninp(3),'b2fgmtry','old','un*formatted')
      inquire(file='b2mn.dat',exist=ex)
      if (ex) call cfopen (ninp(4),'b2mn.dat','old','formatted')
      call cfopen (nout(0),'b2yi.prt','new','formatted')
      call gnopen (nout(1),'b2yi.plt')
*   ..mark output unit for error messages
      call xerset (0)
*   ..obtain version numbers                                             !xpb
      call cfverr (ninp(1), b2fstate_version)                            !xpb
      call cfverr (ninp(2), b2frates_version)
      call cfverr (ninp(3), b2fgmtry_version)                            !xpb
*   ..obtain nx, ny from geometry file
      call cfruin (ninp(3), 2, idum, 'nx,ny')
      nx = idum(0)
      ny = idum(1)
      call xertst (0.le.nx.and.0.le.ny,
     &  'faulty input nx, ny from geometry file')
*   ..obtain nx, ny, ns from plasma state file
      call cfruin (ninp(1), 3, idum, 'nx,ny,ns')
      ns = idum(2)
      call xertst (idum(0).eq.nx.and.idum(1).eq.ny.and.1.le.ns,
     &  'faulty input nx, ny, ns from plasma state file')
*   ..call principal driver routine
      call alloc_b2mod_geo(nx,ny)
      call alloc_b2mod_work(nx,ny,ns)
      call alloc_b2mod_rates(nx,ny,ns)
      call alloc_b2mod_plasma(nx,ny,ns)
      call alloc_b2mod_indirect(nx,ny,nncutmax)                          !xpb
      call alloc_b2mod_ppout(nx,ny,ns)
      call b2yidr (ninp, nout, nx, ny, ns)
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
      call dealloc_b2mod_indirect
      call prgfin ()
      stop 'b2yi'

*-----------------------------------------------------------------------
*.end b2yi

      end program b2yi
