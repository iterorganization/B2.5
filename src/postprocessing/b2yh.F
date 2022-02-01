*-----------------------------------------------------------------------
*.specification

      program b2yh
      use b2mod_types
      use b2mod_geo
      use b2mod_plasma
      use b2mod_rates
      use b2mod_residuals
      use b2mod_sources
      use b2mod_transport
      use b2mod_indirect
      use b2mod_xerset
      implicit none

*-----------------------------------------------------------------------
*.documentation
*
*  1. purpose
*
*     B2YH (main program) is a postprocessor for B2. It reads and
*     displays the table of physics parameters. This program calls the
*     driver routine b2yhdr for the actual computation.
*
*
*  2. specification
*
*     Main program.
*
*
*  3. description (see also routine b2cdca)
*
*     The complete program b2yh performs postprocessing of the
*     result of a B2 calculation. It displays in a convenient format
*     the table of physics parameters.
*     This program unit opens and closes the input/output units, and
*     may perform some other system-dependent operations. It calls
*     routine b2yhdr to control the actual computation.
*
*     The input units are:
*     ninp(0): formatted; provides output control.
*     ninp(1): un*formatted; provides the physics parameters.
*     ninp(2): un*formatted; provides the geometry.
*
*     The output units are:
*     nout(0): formatted; provides printed output.
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
      integer ninp(0:2), nout(0:0), nx, ny, ns, idum(0:9), io
      character*120 lblgm
      logical ex
      logical file_open
*   ..procedures
      external prgini, prgend, xertst, cfopen, cfruin, cfruch
      external b2yhdr, b2rugm
*   ..initialise input/output units
      data ninp/50,51,52/
      data nout/60/

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
*     nout - (0:0) integer array.
*     nout specifies the output unit numbers.
*
*     ns - integer.
*     ns specifies the number of atomic species in the calculation.
*     The species are indexed by (0:ns-1).
*     It will hold that 1.le.ns.
*
*-----------------------------------------------------------------------
*.computation

*   ..program start-up calls
      call prgini ('b2yh')
*   ..open files
      call cfopen (ninp(0),'b2yh.dat','old','formatted')
      inquire(file='b2fparam',exist=ex)
      if (ex) then
        call cfopen (ninp(1),'b2fparam','old','un*formatted')
      else
        call cfopen (ninp(1),'b2fpardf','old','un*formatted')
      endif
      call cfopen (ninp(2),'b2fgmtry','old','un*formatted')
      call cfopen (nout(0),'b2yh.prt','new','formatted')
*   ..mark output unit for error messages
      call xerset (0)
*   ..obtain version numbers                                             !xpb
      call cfverr (ninp(1), b2fparam_version)                            !xpb
      call cfverr (ninp(2), b2fgmtry_version)                            !xpb
*   ..obtain nx, ny from geometry file
      call cfruin (ninp(2), 2, idum, 'nx,ny')
      nx = idum(0)
      ny = idum(1)
      call xertst (0.le.nx.and.0.le.ny,
     &  'faulty input nx, ny from geometry file')
*   ..obtain ns
      call cfruin (ninp(1), 1, idum, 'ns')
      ns = idum(0)
      call xertst (1.le.ns,
     &  'faulty input ns from physics parameters file')
*   ..read geometry and metric
      call alloc_b2mod_geo(nx,ny)
      call alloc_b2mod_indirect(nx,ny,nncutmax)
      call cfruch (ninp(2), 120, lblgm, 'label')
      call cfruin (ninp(2), 1, idum, 'isymm')
      call b2rugm (ninp(2), nx, ny, crx, cry, fpsi, ffbz,
     &  bb, vol, hx, hy, hz, qz, qc, qcb, gs, pbs,
     &  hyl, hxb, hzl, hzb, wbbl, wbbr, wbbv, wbbc,
     &  cell_width, cell_height, gmap)
*   ..call principal driver routine
      call b2yhdr (ninp, nout, nx, ns)
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
      call dealloc_b2mod_indirect
      call prgend ()
      stop 'b2yh'

*-----------------------------------------------------------------------
*.end b2yh

      end program b2yh
