*-----------------------------------------------------------------------
*.specification

      program b2co
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
      use b2mod_constants
      use b2mod_neutrals_namelist
      use b2mod_neutr_src_scaling
#ifdef NAGF95
      use f90_unix           ! IGNORE
#endif
      use b2mod_b2cmfs
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
*     B2CO reads in a formatted plasmastate file, and writes out a
*     a b2fstate file.
*
*
*-----------------------------------------------------------------------
*.declarations

*   ..common blocks

*   ..local variables
#ifdef USE_PXFGETENV
      integer lenval, ierror
#endif
      integer ninp(0:5), nout(0:7), nx, ny, ns, ns0, idum(0:9), io
#ifdef USE_MPI
      integer rank_mpi,size_mpi,ierr_mpi
#endif
#ifndef NAGF95
      integer iargc
#endif
      character savefile*16, arg*256, lblgm*120, lblmn*120, lblcp*120
      logical file_open
*   ..procedures
      external prgini, prgend, xerset, xertst, cfopen, cfruin
      external b2rucp, b2rups, b2wuzd, b2wups, b2xvps
*   ..initialise input/output unit numbers
      data ninp/50,51,52,53,54,55/
      data nout/60,61,62,63,64,65,66,67/
      data savefile/'plasmastate.xxxx'/

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
        stop 'b2co_mpi_child'
      endif
#endif
      call prgini ('b2co')
*   ..open files
      call cfopen (ninp(1),'b2fgmtry','old','un*formatted')
      call cfopen (ninp(2),'b2fpardf','old','un*formatted')
      call cfopen (ninp(3),'b2fstati','old','un*formatted')
      call xertst (iargc().ge.1,
     &  'program should be supplied with a state number')
#ifdef USE_PXFGETENV
      call pxfgetarg(1,arg,lenval,ierror)
#else
      call getarg(1,arg)
#endif
      savefile=arg
      call cfopen (ninp(4), savefile, 'old','unformatted')
      call cfopen (nout(0),'b2co.prt','new','formatted')
      call cfopen (nout(2),'b2fstate','new','un*formatted')
c      call cfopen (nout(7),'b2fplasma','new','unformatted')
*   ..mark output unit for error messages
      call xerset (0)
*   ..obtain version numbers                                    !xpb
      call cfverr (ninp(1), b2fgmtry_version)                   !xpb
      call cfverr (ninp(2), b2fparam_version)                   !xpb
      call cfverr (ninp(3), b2fstati_version)                   !xpb
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
      call alloc_b2mod_indirect(nx,ny,nncutmax)                         !xpb
*   ..read header for geometry
      call cfruch (ninp(1), 120, lblgm, 'label')
      call cfruin (ninp(1), 1, idum, 'isymm')
      isymm = idum(0)
*     (isymm is now stored in the b2cmfs common)               ! xpb 10.08.00
*   ..read geometry and metric
      call b2rugm (ninp(1), nx, ny, crx, cry, fpsi,
     &  ffbz, bb, vol, hx, hy, qz, qc, gs, pbs, wbbl, wbbr, wbbv)
*   ..read physics common
      call cfruch (ninp(2), 120, lblcp, 'label')
      call b2rucp (ninp(2), b2fparam_version, nx, ns)
*   ..read initial state
      call cfverr (ninp(4), plasmastate_version)                        !xpb
      call b2rups (ninp(4), nx, ny, ns, ne, na, ua, uadia, te, ti,      !xpb
     &  po, fna, fhe, fhi, fch, fch_32, fch_52, fch_p, kinrgy,
     &  fna_mdf, fhe_mdf, fhi_mdf, fna_fcor, fna_nodrift, fna_he,
     &  fnaPSch, fhePSch, fhiPSch, fna_eir, fne_eir, fhe_eir, fhi_eir,
     &  fna_32, fna_52, fni_32, fni_52, fne_32, fne_52,
     &  fchdia, fchin, fchvispar, fchvisper, fchvisq, fchinert,
     &  vaecrb, vadia, wadia, floe_noc, floi_noc,
     1  plasmastate_version)
*   ..test plasma state and fluxes
      call b2xvps (nx, ny, ns,
     &  ne, na, ua, te, ti, po, fna, fhe, fhi, fch)
*   ..save final state
      newversion = plasmastate_version
      write(lblmn,*) 'Recovered from ',savefile
      idum(0) = nx
      idum(1) = ny
      idum(2) = ns
      call cfverw (nout(2), plasmastate_version)                        !xpb
      call cfwuin (nout(2), 3, idum, 'nx,ny,ns')
      call cfwuch (nout(2), 120, lblmn, 'label')
      call b2wuzd (nout(2), plasmastate_version, ns,
     &  zamin, zamax, zn, am)
      call b2wups (nout(2), nx, ny, ns, ne, na, ua, uadia, te,          !xpb
     &  ti, po, fna, fhe, fhi, fch, fch_32, fch_52, fch_p, kinrgy,
     &  fna_mdf, fhe_mdf, fhi_mdf, fna_fcor, fna_nodrift, fna_he,
     &  fnaPSch, fhePSch, fhiPSch, fna_eir, fne_eir, fhe_eir, fhi_eir,
     &  fna_32, fna_52, fni_32, fni_52, fne_32, fne_52,
     &  fchdia, fchin, fchvispar, fchvisper, fchvisq, fchinert,
     &  vaecrb, vadia, wadia, floe_noc, floi_noc)
*   ..close graphics files
*     (none)
*   ..write plasma state
c      call cfverw(nout(7), plasmastate_version)
c      call write_b2mod_geo(nout(7))
c      call write_b2mod_plasma(nout(7))
c      call write_b2mod_residuals(nout(7))
c      call write_b2mod_sources(nout(7))
c      call write_b2mod_transport(nout(7))
c      call write_b2mod_anomalous_transport(nout(7))
c      call write_b2mod_neutr_src_scaling(nout(7),nstra)
      write(nout(0),*) 'Converted ',savefile,' to b2fstate'
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
      call dealloc_b2mod_plasma
      call dealloc_b2mod_sources
      call dealloc_b2mod_indirect
      call dealloc_b2mod_residuals
      call dealloc_b2mod_transport
      call dealloc_anomalous_transport
      call prgend ()
#ifdef USE_MPI
      call mpi_barrier(MPI_COMM_WORLD,ierr_mpi)
      call mpi_finalize(ierr_mpi)
#endif
      stop 'b2co'

*-----------------------------------------------------------------------
*.end b2co

      end program b2co

