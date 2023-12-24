*-----------------------------------------------------------------------
*.specification

      program b2xd
      use b2mod_types
      use b2mod_version
      use b2mod_xerset
      use b2us_geo
      use b2us_map
      use b2us_io
      use b2us_plasma
      implicit none

*-----------------------------------------------------------------------
*.documentation
*
*  1. purpose
*
*     B2XD (main program) is a postprocessor for B2. It reads and
*     displays the plasma state. This program calls the driver routine
*     b2xddr for the actual computation.
*
*
*  2. specification
*
*     Main program.
*
*
*  3. description (see also routine b2cdca)
*
*     The complete program b2xd performs postprocessing of the
*     result of a B2 calculation.
*     This program unit opens and closes the input/output units, and
*     may perform some other system-dependent operations. It calls
*     routine b2xddr to control the actual computation.
*
*     The input units are:
*     ninp(0): formatted; provides output control parameters.
*     ninp(1): un*formatted; provides the plasma state.
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

*   ..main data structures
      type(mapping), save  :: mpg
      type(geometry), save :: geo
      type(B2State), save  :: state
*   ..local variables
      integer ninp(0:4), nout(0:2), io, nCv, nFc, ns, idum(0:9)
      logical file_open
*   ..procedures
      external prgini, prgend, xertst, cfopen, cfruin
      external b2xddr
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
      call prgini ('b2xd')
*   ..open files
c      call cfopen (ninp(0),'b2yp.dat','old','formatted')
      call cfopen (ninp(1),'b2fgmtry','old','un*formatted')
      call cfopen (ninp(2),'b2fparam','old','un*formatted')
      call cfopen (ninp(3),'b2fstate','old','un*formatted')
      call cfopen (ninp(4),'b2fplasma','old','unformatted')
*   ..mark output unit for error messages
      call xerset (0)
*   ..obtain version numbers                                             !xpb
      call cfverr (ninp(1), b2fgmtry_version)                            !xpb
      call cfverr (ninp(2), b2fparam_version)                            !xpb
      call cfverr (ninp(3), b2fstate_version)                            !xpb
      call cfverr (ninp(4), b2fplasma_version)
*   ..obtain nCv, nFc, ns from plasma state file
      call cfruin (ninp(3), 3, idum, 'nCv,nFc,ns')
      nCv = idum(0)
      nFc = idum(1)
      ns  = idum(2)
      call xertst (0.le.nCv.and.0.le.nFc.and.1.le.ns,
     &  'faulty input nx, ny, ns from plasma state file')
*   ..obtain mapping/geometry dimensions from geometry file
      call cfruin (ninp(1),7,idum,'nCi,nCg,nCv,nFc,nVx,nFs,nFt')
      mpg%nCi = idum(0)
      mpg%nCg = idum(1)
      mpg%nCv = idum(2)
      mpg%nFc = idum(3)
      mpg%nVx = idum(4)
      mpg%nFs = idum(5)
      mpg%nFt = idum(6)
      call cfruin (ninp(1),5,idum,'nCmxVx,nCmxFc,nVmxCv,nVmxFc,nCmxNv')
      mpg%nCmxVx = idum(0)
      mpg%nCmxFc = idum(1)
      mpg%nVmxCv = idum(2)
      mpg%nVmxFc = idum(3)
      mpg%nCmxNv = idum(4)
      call cfruin (ninp(1),1,idum,'isClassicalGrid')
      mpg%isClassicalGrid = idum(0)
      call cfruin (ninp(1),3,idum,'nx,ny,nncut')
      mpg%nx    = idum(0)
      mpg%ny    = idum(1)
      mpg%nncut = idum(2)
      call xertst (mpg%nCv.eq.nCv.and.mpg%nFc.eq.nFc,
     &  'faulty input nCv, nFc from geometry file')
*   ..obtain parameters from b2fparam file
      call cfruin (ninp(2),1,idum,'ns')
      call xertst (idum(0).eq.ns,
     &  'faulty input ns from physics parameters file')
*   ..obtain parameters from b2fplasma file
      call cfruin (ninp(4),3,idum,'nCv,nFc,ns')
      call xertst (idum(0).eq.nCv.and.idum(1).eq.nFc.and.idum(2).eq.ns,
     &  'faulty input nCv, nFc, ns from b2fplasma file')
*   ..call principal driver routine
      call alloc_mapping (mpg)
      call alloc_geometry (geo, mpg%nCv, mpg%nFc, mpg%nVx, mpg%nFs,
     .  mpg%nFt, mpg%nVmxCv)
      call createB2State (nCv, nFc, ns, state)
*   ..read plasma state
      call read_b2fplasma(ninp(4), nCv, nFc, ns, state)
*   .. call converter
      call b2xddr (ninp, nout, ns, mpg, geo, state)
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
      call destroyB2State (state)
      call dealloc_geometry (geo)
      call dealloc_mapping (mpg)
      call prgend ()
      stop 'b2xd'

*-----------------------------------------------------------------------
*.end b2xd

      end program b2xd
