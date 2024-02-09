*-----------------------------------------------------------------------
*.specification

      program b2mn_opt
      use b2mod_types
      use b2mod_version
      use b2mod_ma28
      use b2mod_ppout
      use b2mod_par_opt
      use b2mod_exp_file_names
      use b2mod_driver
      use b2us_geo
      use b2us_map
      use b2us_plasma
      use b2us_prep
      use b2us_io
      use b2mod_switches
      use b2mod_xerset
#if ( defined(USE_MPI) && MPI_MOD > 0 )
      use mpi      ! IGNORE
#endif
      implicit none
#if ( defined(USE_MPI) && MPI_MOD < 1 )
      include 'mpif.h'
#endif

*-----------------------------------------------------------------------
*.documentation
*
*  1. purpose
*
*     B2MN_OPT is the main program for the B2 code. The code solves a
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
*     ninp(6): un*formatted; provides the initial differentiated state.
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
*     nout(9): un*formatted; provides the final differentiated state.
*     nout(10): un*formatted; provides data on the numerical evolution of differentiated quantities.
*     nout(11): formatted; postprocessor output. [JET]
*
*     (See routine b2cdca for the meaning of 'un*formatted'.)
*
*
*  4. references
*
*     (Refs. /1/ and /2/ are out of date; !! to be replaced.)
*     /1/ B. J. Braams, "Computational Studies in Tokamak Equilibrium
*     and Transport", thesis, State University at Utrecht, June 1986.
*     /2/ B. J. Braams, "A Multi-Fluid Code for Simulation of the Edge
*     Plasma in Tokamaks", NET Report Nr 68, January 1987.
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

*   ..main data structures
      type(mapping), save    :: mpg
      type(geometry), save   :: geo
      type(B2State), save    :: state
      type(B2StateExt), save :: state_ext
      type(B2Average), save  :: state_avg
      type(switches), save   :: switch
*   ..local variables
      integer ix, iy, ic
#ifdef JET
      integer ninp(0:6), nout(0:11), nCv, nFc, ns, ns0, idum(0:9), io   !srv 02.05.07
#else
      integer ninp(0:6), nout(0:10), nCv, nFc, ns, ns0, idum(0:9), io   !srv 02.05.07
#endif
#ifdef USE_MPI
      integer rank_mpi,size_mpi,ierr_mpi
#endif
      logical file_open
#ifdef _OPENMP
      real(kind=R8) :: cpu_s, cpu_f
#else
      real(kind=R4) :: cpu_s, cpu_f
#endif
*   ..procedures
      external prgini, prgend, cfopen, cfruin
*   ..initialise input/output unit numbers
      data ninp/50,51,52,53,54,55,56/                                !srv 02.05.07
#ifdef JET
      data nout/60,61,62,63,64,65,66,67,68,69,70,71/
#else
      data nout/60,61,62,63,64,65,66,67,68,69,70/
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
        stop 'b2mn_opt_mpi_child'
      endif
#endif
      call prgini ('b2mn_opt')
#ifdef _OPENMP
      cpu_s = omp_get_wtime()
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
      call cfopen (ninp(6),'b2fparop.dat','old','unformatted')   !srv 02.05.07
      call cfopen (nout(0),'b2mn.prt','replace','formatted')     !srv 02.05.09 {
      call cfopen (nout(1),'b2fparam','replace','un*formatted')
      call cfopen (nout(2),'b2fstate','replace','un*formatted')
      call cfopen (nout(3),'b2fmovie','replace','un*formatted')
      call cfopen (nout(4),'b2ftrace','replace','un*formatted')
      call cfopen (nout(5),'b2ftrack','replace','un*formatted')
*   ..mark output unit for error messages
      call xerset (0)

*   ..obtain the optimization parameters
      open(70,file='b2op.dat')
      read(70,'(i3)') npar_opt
      allocate(par_opt0(1:npar_opt), par_opt(1:npar_opt),
     &         par_opt_phys(1:npar_opt))
      close(70)
      read(ninp(6)) par_opt_phys                               !srv 02.05.09

*   ..obtain version numbers						 !xpb
      call cfverr (ninp(1), b2fgmtry_version)				 !xpb
      call cfverr (ninp(2), b2fparam_version)				 !xpb
      call cfverr (ninp(3), b2fstati_version)				 !xpb
*   ..obtain nCv, nFc, ns from b2fstate file
      call cfruin (ninp(3), 3, idum, 'nCv,nFc,ns')
      nCv = idum(0)
      nFc = idum(1)
      ns  = idum(2)
      call xertst (0.le.nCv.and.0.le.nFc.and.1.le.ns,
     &  'faulty input nCv, nFc, ns from plasma state file')
*   ..obtain mapping/geometry dimensions from b2fgmtry file
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
*   ..obtain ns from b2fparam file
      call cfruin (ninp(2),1,idum,'ns')
      call xertst (idum(0).eq.ns,
     &  'faulty input ns from physics parameters file')
*   ..call principal driver routine
      call alloc_mapping (mpg)
      call alloc_geometry (geo, mpg%nCv, mpg%nFc, mpg%nVx, mpg%nFs,
     .  mpg%nFt, mpg%nVmxCv)
      call createB2State (nCv, nFc, ns, state)
      call alloc_b2mod_ppout(DEF_NXD, DEF_NYD, ns)
*   ..call driver routines
      call b2mndr_0 (ninp, nout, ns, ns0, switch, geo, mpg,
     &  state, state_ext, state_avg)
      call b2mndr_1 (nout, ns, switch, geo, mpg,
     &  state, state_ext, state_avg, 0)
      call b2mndr_2 (nout, ns, switch, geo, mpg,
     &  state, state_ext, state_avg)
*   ..write plasma state
      call cfopen (nout(7),'b2fplasma','replace','unformatted')  !srv 02.05.09 }
      call cfverw(nout(7), newversion)
      idum(0) = nCv
      idum(1) = nFc
      idum(2) = ns
      call cfwuin (nout(7), 3, idum, 'nCv,nFc,ns')
      call write_b2fplasma (nout(7), nCv, nFc, ns, state)
      close(nout(7))
*   ..write data for optimization
#ifdef WG_TODO
      call write_results(nx, ny)
#endif
#ifdef JET
*   ..set and write ppout
      call cfopen (nout(11),'tran','new','unformatted')
      call set_b2mod_ppout
      call write_b2mod_ppout(nout(11))
      close (nout(11))
#endif
*   ..close files
      do io=0,size(ninp)-1
         inquire(unit=ninp(io),opened=file_open)
         if (file_open) close(ninp(io))
      enddo
      do io=0,5
         inquire(unit=nout(io),opened=file_open)
         if (file_open) close(nout(io))
      enddo

#ifdef _OPENMP
      cpu_f = omp_get_wtime()
#else
      call cpu_time(cpu_f)
#endif
      write(*,"(/' b2mn_opt takes ',f12.1,' sec'/)") cpu_f - cpu_s
      deallocate(par_opt0, par_opt, par_opt_phys)

* ..end of computation
      call dealloc_b2mod_ppout
      call destroyB2State (state)
      call destroyB2StateExt (state_ext)
      call dealloc_geometry (geo)
      call dealloc_mapping (mpg)
      call prgend ()
#ifdef USE_MPI
      call mpi_barrier(MPI_COMM_WORLD,ierr_mpi)
      call mpi_finalize(ierr_mpi)
#endif
      stop 'b2mn_opt'

*-----------------------------------------------------------------------
*.end b2mn_opt

      end program b2mn_opt

