*-----------------------------------------------------------------------
*.specification

      program b2mnastra
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
      use b2mod_ma28
      use b2mod_indirect
      use b2mod_geo_corner
      use b2mod_neutrals_namelist
      use b2mod_neutr_src_scaling
      use b2mod_ppout
      use b2mod_astra_to_b2
      use b2mod_b2_to_astra
      use b2mod_driver
      use b2mod_xerset
#ifdef USE_MPI
#if MPI_MOD > 0
      use mpi      ! IGNORE
#endif
#endif
      implicit none
#ifdef USE_MPI
#if MPI_MOD < 1
      include 'mpif.h'
#endif
#endif

*-----------------------------------------------------------------------
*.documentation
*
*  1. purpose
*
*     B2MNASTRA is the main program for the B2 code. The code solves a
*     system of fluid equations that models a multi-species plasma
*     in a two-dimensional geometry. This program calls the driver
*     routines b2mndr_[0-2] for the actual computation.
*
*
*  2. specification
*
*     Main program.
*
*
*  3. description (see also routine b2cdca)
*
*     The B2 code /1,2/ solves a system of fluid equations that models
*     a multi-species plasma in a two-dimensional geometry. Its main
*     application is simulation of the edge plasma in tokamaks.
*
*     B2 is composed of several independent programs that communicate
*     through external files. The main sequence involves the programs
*     b2ag (prepares a geometry and magnetic field), b2ah (prepares
*     a table of default physics parameters), b2ai (prepares an
*     initial plasma state), b2mn (principal calculation), and b2yp
*     (postprocessing).
*     Additional programs may prepare tables of rate coefficients for
*     atomic processes, tables of plasma-wall interaction processes,
*     external source terms, and other physics, or may perform further
*     postprocessing.
*
*     The complete program b2mn is designed to solve the equations.
*     This program unit opens and closes the input/output units, and
*     may perform some other system-dependent operations. It calls
*     routines b2mndr_[0-2] to control the actual computation.
*
*     The input units are:
*     ninp(0): formatted; provides physics and numerics parameters.
*     ninp(1): un*formatted; provides the geometry.
*     ninp(2): un*formatted; provides default physics parameters.
*     ninp(3): un*formatted; provides the initial state.
*     ninp(4): un*formatted; provides atomic rate data.
*     ninp(5): un*formatted; provides source term data.
*
*     The output units are:
*     nout(0): formatted; provides summary print output.
*     nout(1): un*formatted; provides actual physics parameters.
*     nout(2): un*formatted; provides the final state.
*     nout(3): un*formatted; may provide data for movie.
*     nout(4): un*formatted; provides data on the numerical evolution.
*     nout(5): un*formatted; may provide detailed convergence data.
*     nout(6): un*formatted; save file for restart
*     nout(7): un*formatted; full plasma description with derived quantities
*     nout(8): un*formatted; ASTRA coupling data file
*     nout(9): formatted; postprocessor output. [JET]
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
*  7. accuracy
*
*     As a measure of convergence we employ the scaled L1-norms of the
*     residuals of the individual equations (continuity for each
*     species, momentum balance for each species, electron heat
*     balance, atom heat balance, and charge balance). The scaling
*     provides numbers that have the dimension of (1/time). A reliable
*     estimate of the accuracy of a calculation can be obtained only
*     by continuing the calculation to better convergence.
*
*
*  8. further comments
*
*-----------------------------------------------------------------------
*.declarations

*   ..common blocks

*   ..local variables
#ifdef JET
      integer ninp(0:6), nout(0:9), nx, ny, ns, ns0, idum(0:9), io   !srv 02.05.07
#else
      integer ninp(0:6), nout(0:8), nx, ny, ns, ns0, idum(0:9), io   !srv 02.05.07
#endif
#ifdef USE_MPI
      integer rank_mpi,size_mpi,ierr_mpi
#endif
      real (kind=R4) :: cpu_s, cpu_f
      logical file_open
*   ..procedures
      external prgini, prgend, cfopen, cfruin
*   ..initialise input/output unit numbers
      data ninp/50,51,52,53,54,55,56/                               !srv 02.05.07
#ifdef JET
      data nout/60,61,62,63,64,65,66,67,68,69/
#else
      data nout/60,61,62,63,64,65,66,67,68/
#endif
*-----------------------------------------------------------------------
*.computation
crfs !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
crfs put here the other main programs b2ag (prepares geometry and
crfs  magnetic field), b2ah (prepares tables of default physics
crfs  parameters), b2ai (prepares an initial plasma state),
crfs !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

*   ..program start-up calls
#ifdef USE_MPI
      call mpi_init(ierr_mpi)
      call mpi_comm_rank(MPI_COMM_WORLD,rank_mpi,ierr_mpi)
      call mpi_comm_size(MPI_COMM_WORLD,size_mpi,ierr_mpi)
      if(rank_mpi .ne. 0) then
        call mpi_barrier(MPI_COMM_WORLD,ierr_mpi)
        call mpi_finalize(ierr_mpi)
        stop 'b2mnastra_mpi_child'
      endif
#endif
      call prgini ('b2mnastra')
#ifdef _OPENMP
      cpu_s = real(omp_get_wtime(),R4)
#else
      call cpu_time(cpu_s)
#endif
*   ..open files
      call cfopen (ninp(0),'b2mn.dat','old','formatted')
      call cfopen (ninp(1),'b2fgmtry','old','un*formatted')
      call cfopen (ninp(2),'b2fpardf','old','un*formatted')
      call cfopen (ninp(3),'b2fstati','old','un*formatted')
      call cfopen (ninp(4),'b2frates','old','un*formatted')
*$$$  call cfopen (ninp(5),'b2fsourc','old','un*formatted')
*     (ninp(5) is reserved for interface to a neutral gas code)
      call cfopen (ninp(6),'astra_to_b2.dat','old','un*formatted')
      call cfopen (nout(0),'b2mn.prt','replace','formatted')       !srv 02.05.09 {
      call cfopen (nout(1),'b2fparam','replace','un*formatted')
      call cfopen (nout(2),'b2fstate','replace','un*formatted')
      call cfopen (nout(3),'b2fmovie','replace','un*formatted')
      call cfopen (nout(4),'b2ftrace','replace','un*formatted')
      call cfopen (nout(5),'b2ftrack','replace','un*formatted')
*   ..mark output unit for error messages
      call xerset (0)

*   ..obtain version numbers                                             !xpb
      call cfverr (ninp(1), b2fgmtry_version)                            !xpb
      call cfverr (ninp(2), b2fparam_version)                            !xpb
      call cfverr (ninp(3), b2fstati_version)                            !xpb
*   ..obtain nx, ny from geometry file
      call cfruin (ninp(1), 2, idum, 'nx,ny')
      nx = idum(0)
      ny = idum(1)
      call xertst (0.le.nx.and.0.le.ny,
     &  'faulty input nx, ny from geometry file')
*   ..obtain ns from parameters file
      call cfruin (ninp(2), 1, idum, 'ns')
      ns = idum(0)
      call xertst (1.le.ns, 'faulty input ns from parameters file')
*   ..obtain nx, ny, ns0 from plasma state file
      call cfruin (ninp(3), 3, idum, 'nx,ny,ns')
      ns0 = idum(2)
      call xertst (idum(0).eq.nx.and.idum(1).eq.ny,
     &  'faulty input nx, ny from plasma state file')
      call xertst (1.le.ns0.and.ns0.le.ns,
     &  'faulty input ns from plasma state file')
*   ..call principal driver routine
      call alloc_b2mod_geo(nx,ny)
      call alloc_b2mod_plasma(nx,ny,ns)
      call alloc_b2mod_rates(nx,ny,ns)
      call alloc_b2mod_residuals(nx,ny,ns)
      call alloc_b2mod_sources(nx,ny,ns)
      call alloc_b2mod_transport(nx,ny,ns)
      call alloc_b2mod_anomalous_transport(nx,ny,ns)
      call alloc_b2mod_work(nx,ny,ns)
      call alloc_b2mod_indirect(nx,ny,nncutmax)                          !xpb
      call alloc_b2mod_ppout(nx,ny,ns)
      call alloc_b2mod_astra_to_b2(nx, ny)                      !srv 28.07.11 {
      call alloc_b2mod_b2_to_astra(nx, ny, ns)
      call read_astra_to_b2(ny, ninp(6))
      open(70,file='check_astra_to_b2.prt')
      call write_to_check_astra_to_b2(70)
      close(70)                                                 !srv 28.07.11 }
*   ..call driver routines
!WG_TODO      call b2mndr_0 (ninp, nout, nx, ny, ns, ns0)
!WG_TODO      call b2mndr_1 (nout, nx, ny, ns, 0)
!WG_TODO      call b2mndr_2 (nout, nx, ny, ns)
*   ..write plasma state
      call cfopen (nout(7),'b2fplasma','replace','unformatted')
      call cfverw(nout(7), newversion)
      call write_b2mod_geo(nout(7))
      call write_b2mod_plasma(nout(7))
      call write_b2mod_residuals(nout(7))
      call write_b2mod_sources(nout(7))
      call write_b2mod_transport(nout(7))
      call write_b2mod_anomalous_transport(nout(7))
      call write_b2mod_neutr_src_scaling(nout(7),nstrai)
      close(nout(7))

#ifdef JET
*   ..set and write ppout
      call cfopen (nout(9),'tran','new','unformatted')
      call set_b2mod_ppout
      call write_b2mod_ppout(nout(9))
      close (nout(9))
#endif
*   ..close files
      do io=0,size(ninp)-1
         inquire(unit=ninp(io),opened=file_open)
         if (file_open) close(ninp(io))
      enddo
      do io=0,size(nout)-1
         inquire(unit=nout(io),opened=file_open)
         if (file_open) close(nout(io))
      enddo

#ifdef _OPENMP
      cpu_f = real(omp_get_wtime(),R4)
#else
      call cpu_time(cpu_f)
#endif
      write(*,"(/' b2mnastra takes ',f12.1,' sec'/)") cpu_f - cpu_s

* ..end of computation
      call dealloc_b2mod_geo
      call dealloc_b2mod_work
      call dealloc_b2mod_ppout
      call dealloc_b2mod_rates
      call dealloc_b2mod_plasma
      call dealloc_b2mod_sources
      call dealloc_b2mod_indirect
      call dealloc_b2mod_residuals
      call dealloc_b2mod_transport
      call dealloc_b2mod_astra_to_b2
      call dealloc_b2mod_b2_to_astra
      call dealloc_anomalous_transport
      call prgend ()
#ifdef USE_MPI
      call mpi_barrier(MPI_COMM_WORLD,ierr_mpi)
      call mpi_finalize(ierr_mpi)
#endif
      stop 'b2mnastra'

*-----------------------------------------------------------------------
*.end b2mnastra

      end program b2mnastra

