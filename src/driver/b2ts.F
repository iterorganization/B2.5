*-----------------------------------------------------------------------
*.specification

      program b2ts
      use b2mod_types
      use b2mod_version
      use b2mod_geo
      use b2mod_plasma
      use b2mod_rates
      use b2mod_residuals
      use b2mod_sources
      use b2mod_transport
      use b2mod_anomalous_transport
      use b2mod_work
      use b2mod_indirect
      use b2mod_neutrals_namelist
      use b2mod_neutr_src_scaling
      implicit none

*-----------------------------------------------------------------------
*.documentation
*
*  1. purpose
*
*     B2TS (main program) is a postprocessor for b2.
*
*
*  2. specification
*
*     Main program.
*
*
*  3. description (see also routine b2cdca)
*
*     The complete program b2ts performs postprocessing of the
*     result of a b2 calculation.
*     This program unit opens and closes the input units, and
*     may perform some other system-dependent operations. It calls
*     routine b2tsdr to control the actual computation.
*
*     The input units are:
*     ninp(0): formatted; provides output control parameters.
*     ninp(1): un*formatted; provides the plasma state.
*     ninp(2): un*formatted; provides the physics parameters.
*     ninp(3): un*formatted; provides the plasma state.
*     ninp(4): un*formatted;
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
      integer ninp(0:4), nx, ny, ns, idum(0:9), io
      logical file_open
*   ..procedures
      external prgini, prgend, xerset, xertst, cfopen, cfruin
!      external b2tsdr
*   ..initialise input units
      data ninp/50,51,52,53,54/

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
      call prgini ('b2ts')
*   ..open files
c      call cfopen (ninp(0),'b2ts.dat','old','formatted')
      call cfopen (ninp(1),'b2fgmtry','old','un*formatted')
      call cfopen (ninp(2),'b2fparam','old','un*formatted')
      call cfopen (ninp(3),'b2fstate','old','un*formatted')
      call cfopen (ninp(4),'b2fplasma','old','unformatted')
*   ..mark output unit for error messages
      call xerset (0)
*   ..obtain version numbers						 !xpb
      call cfverr (ninp(1), b2fgmtry_version)				 !xpb
      call cfverr (ninp(2), b2fparam_version)				 !xpb
      call cfverr (ninp(3), b2fstate_version)				 !xpb
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
      call alloc_b2mod_plasma(nx,ny,ns)
      call alloc_b2mod_rates(nx,ny,ns)
      call alloc_b2mod_residuals(nx,ny,ns)
      call alloc_b2mod_sources(nx,ny,ns)
      call alloc_b2mod_transport(nx,ny,ns)
      call alloc_b2mod_anomalous_transport(nx,ny,ns)
      call alloc_b2mod_work(nx,ny,ns)
      call alloc_b2mod_indirect(nx,ny,nncutmax)
*   ..read plasma state
      call cfverr(ninp(4),b2fplasma_version)
      call read_b2mod_geo(nx,ny,ninp(4))
      call read_b2mod_plasma(nx,ny,ns,ninp(4))
      call read_b2mod_residuals(ninp(4))
      call read_b2mod_sources(ninp(4))
      call read_b2mod_transport(ninp(4))
      call read_b2mod_anomalous_transport(ninp(4))
      call read_b2mod_neutr_src_scaling(ninp(4),ns,nstrai)
      write(*,*) crx(5*nx/8-2,0,0:3)
      write(*,*) cry(5*nx/8-2,0,0:3)
      write(*,*) crx(5*nx/8-1,0,0:3)
      write(*,*) cry(5*nx/8-1,0,0:3)
      write(*,*) crx(5*nx/8,0,0:3)
      write(*,*) cry(5*nx/8,0,0:3)
      write(*,*) crx(5*nx/8+1,0,0:3)
      write(*,*) cry(5*nx/8+1,0,0:3)
*   ..close files
      do io=0,size(ninp)-1
         inquire(unit=ninp(io),opened=file_open)
         if (file_open) close(ninp(io))
      enddo

* ..end of computation
      call dealloc_b2mod_geo
      call dealloc_b2mod_work
      call dealloc_b2mod_rates
      call dealloc_b2mod_plasma
      call dealloc_b2mod_sources
      call dealloc_b2mod_indirect
      call dealloc_b2mod_residuals
      call dealloc_b2mod_transport
      call dealloc_anomalous_transport
      call dealloc_b2mod_neutr_src_scaling
      call prgend ()
      stop 'b2ts'

*-----------------------------------------------------------------------
*.end b2ts

      end program b2ts
