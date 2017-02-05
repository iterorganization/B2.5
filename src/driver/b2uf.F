*-----------------------------------------------------------------------
*.specification

      program b2uf
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
      use b2mod_neutrals_namelist
      use b2mod_neutr_src_scaling
      use b2mod_b2cmpa
      implicit none
#ifdef USE_MPI
      include 'mpif.h'
#endif

*-----------------------------------------------------------------------
*.documentation
*
*  1. purpose
*
*     B2UF reads in an unformatted b2fplasma file, and writes out a
*     a formatted b2fplasma file.
*
*
*-----------------------------------------------------------------------
*.declarations

*   ..common blocks

*   ..local variables
      integer ninp(0:5), nout(0:7), nx, ny, ns, idum(0:9), io
#ifdef USE_MPI
      integer rank_mpi,size_mpi,ierr_mpi
#endif
      logical file_open
*   ..procedures
      external prgini, prgend, xerset, xertst, cfopen, cfruin
*   ..initialise input/output unit numbers
      data ninp/50,51,52,53,54,55/
      data nout/60,61,62,63,64,65,66,67/

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
        stop 'b2uf_mpi_child'
      endif
#endif
      call prgini ('b2uf')
*   ..open files
      call cfopen (ninp(1),'b2fgmtry','old','un*formatted')
      call cfopen (ninp(2),'b2fparam','old','un*formatted')
      call cfopen (ninp(3),'b2fstate','old','un*formatted')
      call cfopen (ninp(4),'b2fplasma','old','unformatted')
      call cfopen (nout(0),'b2uf.prt','new','formatted')
      call cfopen (nout(7),'b2fplasmf','new','formatted')
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
      call alloc_b2mod_indirect(nx,ny,nncutmax)                 	 !xpb
*   ..read plasma state
      call cfverr(ninp(4), b2fplasma_version)
      call read_b2mod_geo(nx,ny,ninp(4))
      call read_b2mod_plasma(nx,ny,ns,ninp(4))
      call read_b2mod_residuals(ninp(4))
      call read_b2mod_sources(ninp(4))
      call read_b2mod_transport(ninp(4))
      call read_b2mod_anomalous_transport(ninp(4))
      call read_b2mod_neutr_src_scaling(ninp(4),ns,nstrai)
      call read_b2mod_eirene_sources(nx,ny,ns,nstrai,ninp(4))
*   ..write plasma state
      newversion=b2fplasma_version
      call cfverw(nout(7), newversion)
      call write_b2mod_geo(nout(7))
      call write_b2mod_plasma(nout(7))
      call write_b2mod_residuals(nout(7))
      call write_b2mod_sources(nout(7))
      call write_b2mod_transport(nout(7))
      call write_b2mod_anomalous_transport(nout(7))
      call write_b2mod_neutr_src_scaling(nout(7),nstrai)
      call write_b2mod_eirene_sources(nout(7))
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
      call dealloc_b2mod_work
      call dealloc_b2mod_rates
      call dealloc_b2mod_plasma
      call dealloc_b2mod_sources
      call dealloc_b2mod_indirect
      call dealloc_b2mod_residuals
      call dealloc_b2mod_transport
      call dealloc_anomalous_transport
      call dealloc_b2mod_neutr_src_scaling
      call dealloc_b2mod_eirene_sources
      call prgend ()
#ifdef USE_MPI
      call mpi_barrier(MPI_COMM_WORLD,ierr_mpi)
      call mpi_finalize(ierr_mpi)
#endif
      stop 'b2uf'

*-----------------------------------------------------------------------
*.end b2uf

      end program b2uf
