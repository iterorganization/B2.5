!> This is a restructured version of b2mn
!>

      module b2mod_main

      use b2mod_types
      use b2mod_driver
      use b2mod_version
      use b2mod_neutrals_namelist
      use b2mod_neutr_src_scaling
      use b2mod_ppout
!WG_TODO      use b2mod_b2_to_astra
!WG_TODO      use b2mod_astra_to_b2
      use b2us_geo
      use b2us_map
      use b2us_plasma
      use b2us_io
      use b2mod_switches
      use b2mod_recycle
      use b2mod_xerset
#ifdef B25_EIRENE
#ifdef USE_MPI
      use eirmod_parmmod
     , , only : ifoff
      use eirmod_cpes
     , , only : my_pe
#endif
#endif
#if MPI_MOD > 0
#ifdef USE_MPI
      use mpi  ! IGNORE
#else
#ifdef USE_IMPGYRO
      use mpi  ! IGNORE
#endif
#endif
#endif
      implicit none
#if MPI_MOD < 1
#ifdef USE_MPI
      include 'mpif.h'
#else
#ifdef USE_IMPGYRO
      include 'mpif.h'
#endif
#endif
#endif

*-----------------------------------------------------------------------
*.documentation
*
*  1. purpose
*
*     B2MOD_MAIN contains the main program for the B2 code. The code solves a
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
*     species, parallel momentum balance for each species, electron heat
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
      integer ninp(0:6), nout(0:11), nCv, nFc, nVx, ns, ns0
#else
      integer ninp(0:6), nout(0:10), nCv, nFc, nVx, ns, ns0
#endif
#ifdef USE_MPI
      integer rank_mpi,size_mpi,ierr_mpi
#ifdef B25_EIRENE
      real (kind=R8) :: dumarr(100)
      logical ltrigger
#endif
#else
#ifdef USE_IMPGYRO
      integer rank_mpi,size_mpi,ierr_mpi
#endif
#endif
      character*32 B25_git_version
      logical file_open
*   ..procedures
      character*32 get_B25_hash
      external get_B25_hash
      external prgini, prgend, cfopen, cfruin
*   ..initialise input/output unit numbers
      data ninp/50,51,52,53,54,55,56/
#ifdef JET
      data nout/60,61,62,63,64,65,66,67,68,69,70,71/
#else
      data nout/60,61,62,63,64,65,66,67,68,69,70/
#endif

      contains

      subroutine b2mn_init (switch, geo, mpg, state, state_ext)
      implicit none
      type(switches), intent (inout)   :: switch
      type(mapping), intent (inout)    :: mpg
      type(geometry), intent (inout)   :: geo
      type(B2State), intent (inout)    :: state
      type(B2StateExt), intent (inout) :: state_ext
      integer idum(0:9)

*   ..program start-up calls
#ifdef USE_MPI
csw 28feb2011 mpi
      call mpi_init(ierr_mpi)
      call mpi_comm_rank(MPI_COMM_WORLD,rank_mpi,ierr_mpi)
      call mpi_comm_size(MPI_COMM_WORLD,size_mpi,ierr_mpi)
      if(rank_mpi .ne. 0) then
#ifdef B25_EIRENE
        my_pe=rank_mpi
        iunout=7+ifoff
        call eirene_find_param(0)
        do
          call mpi_bcast(ltrigger,1,MPI_LOGICAL,0,
     .                   MPI_COMM_WORLD,ierr_mpi)
          if(ltrigger) then
            call eirene_eirsrt(.false.,.false.,
     .       0.0_R8,dumarr,0.0_R8,0.0_R8,0.0_R8,0.0_R8,0.0_R8)
          else
            exit
          endif
        enddo
csw
#endif
        call mpi_barrier(MPI_COMM_WORLD,ierr_mpi)
        call mpi_finalize(ierr_mpi)
        stop 'b2mn_mpi_child'
      endif
#else
#ifdef USE_IMPGYRO
      write(*,*) 'b2ig: running B2.5-EIRENE'
      call mpi_init(ierr_mpi)
      call mpi_comm_rank(MPI_COMM_WORLD,rank_mpi,ierr_mpi)
      call mpi_comm_size(MPI_COMM_WORLD,size_mpi,ierr_mpi)
      write(*,*) 'b2ig: b2mn size_mpi =', size_mpi
#endif
#endif

      call prgini ('b2mn')
*   ..open files
      call cfopen (ninp(0),'b2mn.dat','old','formatted')
      call cfopen (ninp(1),'b2fgmtry','old','un*formatted')
      call cfopen (ninp(2),'b2fpardf','old','un*formatted')
      call cfopen (ninp(3),'b2fstati','old','un*formatted')
      call cfopen (ninp(4),'b2frates','old','un*formatted')
*$$$  call cfopen (ninp(5),'b2fsourc','old','un*formatted')
*     (ninp(5) is reserved for interface to a neutral gas code)
      call cfopen (nout(0),'b2mn.prt','new','formatted')
      call cfopen (nout(1),'b2fparam','new','un*formatted')
      call cfopen (nout(2),'b2fstate','new','un*formatted')
      call cfopen (nout(3),'b2fmovie','new','un*formatted')
      call cfopen (nout(4),'b2ftrace','new','un*formatted')
      call cfopen (nout(5),'b2ftrack','new','un*formatted')
*   ..mark output unit for error messages
      call xerset (0)
*   ..obtain version numbers                                             !xpb
      B25_git_version = get_B25_hash()
      call cfverw (nout(0), newversion)
      write(6,'(4a)')
     &  'Running with VERSION ',newversion,' ', trim(B25_git_version)
      call cfverr (ninp(1), b2fgmtry_version)                            !xpb
      call cfverr (ninp(2), b2fparam_version)                            !xpb
      call cfverr (ninp(3), b2fstati_version)                            !xpb
      if (b2fgmtry_version .lt. '03.002.000'  .or.
     &    b2fstati_version .lt. '03.002.000') then
        call xerrab('Using old input format. Run b2us first.')
      endif
*   ..obtain mapping/geometry dimensions from file
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
*   ..obtain ns from parameters file
      call cfruin (ninp(2), 1, idum, 'ns')
      ns = idum(0)
      call xertst (1.le.ns, 'faulty input ns from parameters file')
*   ..obtain plasma state from file
      call cfruin (ninp(3), 3, idum, 'nCv,nFc,ns')
      ns0 = idum(2)
      call xertst (idum(0).eq.mpg%nCv.and.idum(1).eq.mpg%nFc,
     &  'faulty input nCv, nFc from plasma state file')
      call xertst (1.le.ns0.and.ns0.le.ns,
     &  'faulty input ns from plasma state file')
*   ..obtain external state
      call alloc_mapping (mpg)
      call alloc_geometry (geo, mpg%nCv, mpg%nFc, mpg%nVx, mpg%nFs, 
     .  mpg%nFt, mpg%nVmxCv)
      call createB2State (mpg%nCv, mpg%nFc, ns, state)
      call createB2StateExt ( mpg%nCv, mpg%nFc, ns, state_ext )
*   ..call principal driver routine
      call alloc_b2mod_ppout(1,1,ns)
      call alloc_b2mod_recycle(mpg%nFc,ns)   !nh 30.03.21
!WG_TODO      call alloc_b2mod_astra_to_b2(nx, ny)
!WG_TODO      call alloc_b2mod_b2_to_astra(nx, ny, ns)

      call b2mndr_0 (ninp, nout, ns, ns0, switch, geo, mpg, 
     &               state, state_ext)

      return
      end subroutine b2mn_init

      subroutine b2mn_step(switch, geo, mpg, state, state_ext, J)
      use b2mod_types
      use b2mod_b2cmpt
      use b2mod_b2cmrc
      use b2mod_boundary_namelist
      use b2mod_neutrals_namelist
      use b2mod_user_namelist
     , , only : sigma, nsigma
      use b2mod_input_profile
     , , only : tdata
      use b2mod_ad
     , , only: nncf
      use b2mod_par_opt
     & , only : par_opt_phys, npar_opt
      implicit none
      type(switches), intent (inout)   :: switch
      type(mapping), intent (inout)    :: mpg
      type(geometry), intent (inout)   :: geo
      type(B2State), intent (inout)    :: state
      type(B2StateExt), intent (inout) :: state_ext
      real(kind=R8) :: J(nncf)

      call b2mndr_1 (nout, ns, switch, geo, mpg, state, state_ext, J)

      return
      end subroutine b2mn_step

      subroutine b2mn_fin (switch, geo, mpg, state, state_ext)
      implicit none
      type(switches), intent (inout)   :: switch
      type(mapping), intent (inout)    :: mpg
      type(geometry), intent (inout)   :: geo
      type(B2State), intent (inout)    :: state
      type(B2StateExt), intent (inout) :: state_ext
      integer io

      call b2mndr_2 (nout, ns, switch, geo, mpg, state, state_ext)
#ifdef JET
*   ..set and write ppout
      call cfopen (nout(10),'tran','new','unformatted')
      call set_b2mod_ppout
      call write_b2mod_ppout(nout(10))
      close (nout(10))
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

* ..end of computation
      call destroyB2State (state)
      call destroyB2StateExt (state_ext)
      call dealloc_geometry (geo)
      call dealloc_mapping (mpg)
      call dealloc_b2mod_ppout
!WG_TODO      call dealloc_b2mod_astra_to_b2
!WG_TODO      call dealloc_b2mod_b2_to_astra
      call dealloc_b2mod_neutr_src_scaling
      call dealloc_b2mod_recycle   !nh 30.03.21
      call prgend ()

#ifdef USE_MPI
#ifdef B25_EIRENE
      ltrigger=.false.
      call mpi_bcast(ltrigger,1,MPI_LOGICAL,0,
     .                   MPI_COMM_WORLD,ierr_mpi)
#endif
      call mpi_barrier(MPI_COMM_WORLD,ierr_mpi)
      call mpi_finalize(ierr_mpi)
#else
#ifdef USE_IMPGYRO
      call mpi_barrier(MPI_COMM_WORLD,ierr_mpi)
      call mpi_finalize(ierr_mpi)
#endif
#endif
      return
      end subroutine b2mn_fin

      end module b2mod_main

