*-----------------------------------------------------------------------
*.specification

      program b2ym
      use b2mod_types
      use b2mod_version
      use b2mod_geo
      use b2mod_plasma
      use b2mod_rates
      use b2mod_residuals
      use b2mod_sources
      use b2mod_transport
      use b2mod_work
      implicit none
#ifdef USE_MPI
      include 'mpif.h'
#endif

*-----------------------------------------------------------------------
*.documentation
*
*  1. purpose
*
*     INCLUDE 'doc/b2ym.purp'
*
*
*  2. specification
*
*     Main program.
*
*
*  3. description (see also routine b2cdca)
*
*     The complete program b2ym is meant to produce a movie to display
*     the time-evolution of a b2 calculation. Presently only a sequence
*     of contour plots is produced.
*     This program unit opens and closes the input/output units, and
*     may perform some other system-dependent operations. It calls
*     routine b2ymdr to control the actual computation.
*
*     The input units are:
*     ninp(0): formatted; provides output control parameters.
*     ninp(1): un*formatted; provides the plasma state evolution.
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
      integer ninp(0:1), nout(0:1), nx, ny, ns, idum(0:9)
#ifdef USE_MPI
      integer rank_mpi,size_mpi,ierr_mpi
#endif
*   ..procedures
      external prgini, prgend, xerset, xertst, cfopen, gnopen, cfruin,
     &  gopks, gclks, gclwk
      external b2ymdr
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
#ifdef USE_MPI
      call mpi_init(ierr_mpi)
      call mpi_comm_rank(MPI_COMM_WORLD,rank_mpi,ierr_mpi)
      call mpi_comm_size(MPI_COMM_WORLD,size_mpi,ierr_mpi)
      if(rank_mpi .ne. 0) then
        call mpi_barrier(MPI_COMM_WORLD,ierr_mpi)
        call mpi_finalize(ierr_mpi)
        stop 'b2ym_mpi_child'
      endif
#endif
      call prgini ('b2ym')
*   ..open gks
      call gopks (0, 0)
*   ..open files
      call cfopen (ninp(0),'b2ym.dat','old','formatted')
      call cfopen (ninp(1),'b2fmovie','old','un*formatted')
      call cfopen (nout(0),'b2ym.prt','new','formatted')
      call gnopen (nout(1),'b2ym.plt')
*   ..mark output unit for error messages
      call xerset (0)
*   ..obtain version number						 !xpb
      call cfverr (ninp(1), b2fmovie_version)				 !xpb
*   ..obtain nx, ny, ns
      call cfruin (ninp(1), 3, idum, 'nx,ny,ns')
      nx = idum(0)
      ny = idum(1)
      ns = idum(2)
      call xertst (0.le.nx.and.0.le.ny.and.1.le.ns,
     &  'faulty input nx, ny, ns')
*   ..call principal driver routine
      call alloc_b2mod_geo(nx,ny)
      call alloc_b2mod_plasma(nx,ny,ns)
      call b2ymdr (ninp, nout, nx, ny, ns)
*   ..close graphics files
      call gclwk (nout(1))
*   ..close gks
      call gclks ()

* ..end of computation
      call dealloc_b2mod_geo
      call dealloc_b2mod_plasma
      call prgend ()
#ifdef USE_MPI
      call mpi_barrier(MPI_COMM_WORLD,ierr_mpi)
      call mpi_finalize(ierr_mpi)
#endif
      stop 'b2ym'

*-----------------------------------------------------------------------
*.end b2ym

      end program b2ym
