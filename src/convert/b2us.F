*-----------------------------------------------------------------------
*.specification

      program b2us
      use b2mod_types
      use b2mod_version
      use b2mod_time
      use b2mod_geo
      use b2mod_plasma
      use b2mod_rates
      use b2mod_residuals
      use b2mod_sources
      use b2mod_sputter
      use b2mod_transport
      use b2mod_anomalous_transport
      use b2mod_work
      use b2mod_diag
      use b2mod_b2plot
      use b2mod_indirect
      use b2mod_external
      use b2mod_elements
      use b2mod_constants
      use b2mod_user_namelist
      use b2mod_boundary_namelist
      use b2mod_neutrals_namelist
      use b2mod_neutr_src_scaling

      use b2us_map
      use b2us_geo
      use b2us_plasma
      use b2us_prep
#ifdef NAGF95
      use f90_unix           ! IGNORE
#endif
      use b2mod_b2cmfs
      use b2mod_b2cmpa
#ifdef B25_EIRENE
      use eirmod_parmmod
     ,   , ncv_eir => ncv
      use eirmod_clogau
      use eirmod_ctrcei
      use eirmod_braeir
      use eirmod_cinit
     , , only: eirene_init_cinit
      use eirmod_comprt
     , , only: iunout
      use eirmod_eirbra
     & , smo_eir => smo
#else
      use b2mod_braeir
      use b2mod_eirbra
     & , smo_eir => smo
#endif
      implicit none
#ifdef USE_MPI
      include 'mpif.h'
#endif
#include "DIMENSIONS.F"

*-----------------------------------------------------------------------
*.documentation
*
*  1. purpose
*
*     B2US reads in structured geometry and plasmastate files, 
*     converts the data to unstructured format and writes 
*     unstructured geometry and plasmastate files
*
*
*-----------------------------------------------------------------------
*.declarations

*   ..common blocks

*   ..local variables
#ifdef USE_PXFGETENV
      integer lenval, ierror
#endif
      integer ninp(0:5), nout(0:7), ns0, idum(0:9), io
      integer ix, iy
      integer iflag
#ifdef USE_MPI
      integer rank_mpi,size_mpi,ierr_mpi
#endif
#ifndef NAGF95
      integer iargc
#endif
      integer :: b2time = 1, ntim_batch = 500
      integer :: ismain = 1, b2wall_netcdf = 0, mvnum = 0, tally = 1
      integer :: use_eirene = 0, external_species = 0
      integer :: boundary_namelist=0, neutrals_namelist = 0
      real (kind=R8) :: BoRiS = 0.0_R8
      real (kind=R8) :: dtim = 1.0_R8
      real (kind=R8), save :: fac_hz = 0.0e0_R8                      !srv 09.01.01 09.06.08
      real (kind=R8) :: hzconst                                      !srv 30.01.10
      real (kind=R8) ::
     *   zpdum(0:nsdecl-1,4), t0, t1
      character savefile*256, arg*256, lblgm*120, lblmn*120, lblcp*120,
     *          lblps*120
      character filename*256
      character extension*4
      character cnamip*80, cvalip*80, line*500
      character bcp_filename*260
#ifdef B25_EIRENE
      character*72 zeile
#endif
      logical file_open, file_exists, ex, 
     .        linq_33, linq_34, linq_35, lcheck, file_exist

      type(mapping), save :: mpg
      type(geometry), save :: gm
      type(B2State), save :: state
      type(B2StateExt), save :: state_ext
      integer :: ntm, is, io_33, io_34, io_35, ico, iel, iel2, natm_loc
*   ..procedures
      integer lnblnk
      logical is_comment
      external lnblnk, is_comment
      external prgini, prgend, xerset, xertst, cfopen, cfruin
      external b2rucp, b2rups, b2wuzd, b2wups, b2xvps, ipgeti, ipgetr
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
        stop 'b2us_mpi_child'
      endif
#endif
      call prgini ('b2us')
*   ..open files
      call cfopen (ninp(0),'b2mn.dat','old','formatted')
      call cfopen (ninp(1),'b2fgmtry','old','un*formatted')
      call cfopen (ninp(2),'b2fpardf','old','un*formatted')
      call cfopen (ninp(3),'b2fstati','old','un*formatted')
      call cfopen (nout(0),'b2us.prt','new','formatted')
      call cfopen (nout(1),'b2fgmtry_us','new','formatted')
      call cfopen (nout(2),'b2fstati_us','new','un*formatted')
*   ..mark output unit for error messages
      call xerset (0)
*   ..obtain version numbers                                           !xpb
      call cfverr (ninp(1), b2fgmtry_version)                          !xpb
      call cfverr (ninp(2), b2fparam_version)                          !xpb
      call cfverr (ninp(3), b2fstati_version)                          !xpb
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
!WG_TODO      call alloc_b2mod_sputter(nx,ny,ns)
      call alloc_b2mod_transport(nx,ny,ns)
      call alloc_b2mod_anomalous_transport(nx,ny,ns)
      call alloc_b2mod_work(nx,ny,ns)
      call alloc_b2mod_indirect(nx,ny,nncutmax)                         !xpb
*   ..obtain lblmn
      call b2xxid ('b2mn', lblmn(1:60))
      call b2rflb (ninp(0), lblmn(61:120))
      label = lblmn
*   ..read physics common
      call cfruch (ninp(2), 120, lblcp, 'label')
      call b2rucp (ninp(2), b2fparam_version, ns)
*   ..read revised parameters in physics common
*     (ninp(0) need supply only parameters that are to be changed from
*     the values read in as default physics common.)
      call b2rfcp (ninp(0), ns)
*   ..test parameters
      call b2xvcp (ns)
*   ..read and echo code internal parameters
    1 continue
       read (ninp(0),'(a)',end=2,err=91) line
       if (.not.is_comment(line)) then
         read (line,*) cnamip, cvalip
         call ipsetc (cnamip, cvalip)
       end if
       goto 1
    2 continue
*   ..read header for geometry
      call cfruch (ninp(1), 120, lblgm, 'label')
      call cfruin (ninp(1), 1, idum, 'isymm')
      isymm = idum(0)
*     (isymm is now stored in the b2cmfs common)               ! xpb 10.08.00
*   ..read geometry and metric
      call b2rugm (ninp(1), nx, ny, isymm, crx, cry, fpsi, ffbz,
     &  bb, vol, hx, hy, hc, ht, qz, qc, qs, qac, qas, gs, pbs, wbbl,
     &  wbbr, wbbv, wbbc, ebc, cell_width, cell_height, gmap)
      call get_nxt(nx, nxtl, nxtr)

      call ipgeti ('b2mndr_eirene', use_eirene)	
      call ipgetr ('b2mndr_hz', fac_hz)                       !srv 09.01.01 09.06.08 {
      call xertst (0.0.le.fac_hz .and. 1.0.ge.fac_hz,
     &  'faulty parameter b2mndr_hz')                          !srv 09.01.01 09.06.08 }
*   ..compute hz                                                         !srv 09.01.01  09.06.08 {
!WG_TODO probably should always write correct hz here,
!WG_TODO and overwrite with value 1.0 in b2mn if b2mndr_hz = 0.0
c      hz = 1.0e0_R8*(1.0e0_R8-fac_hz)+fac_hz*vol/hx/hy
!      hz = 1.0e0_R8*(1.0e0_R8-fac_hz)+fac_hz*vol/gs(:,:,2)              !srv 09.01.01
      hzconst = 2.0_R8*pi*6.0_R8   ! for ITER for testing
      hzconst = 2.0_R8*pi*1.6_R8   ! for ASDEX-UPGRADE for testing
      hzconst = 1.0_R8                                                   !srv 30.01.10
      hz = hzconst*(1.0e0_R8-fac_hz)+fac_hz*vol/gs(:,:,2)                !srv 16.06.09
*   ..read initial state
      call cfruch (ninp(3), 120, lblps, 'label')
      call b2ruzd (ninp(3), b2fstati_version, ns0, 
     & zpdum(0,1), zpdum(0,2), zpdum(0,3), zpdum(0,4), .false.)
      call b2rups (ninp(3), nx, ny, ns, ne, na, ua, uadia, te, ti,      !xpb
     &  tn, po, kt, fna, fhe, fhi, fch, fch_32, fch_52, fch_p, kinrgy,
     &  fna_mdf, fhe_mdf, fhi_mdf, fna_fcor, fna_nodrift, fna_he,
     &  fnaPSch, fhePSch, fhiPSch, fna_eir, fne_eir, fhe_eir, fhi_eir,
     &  fna_32, fna_52, fni_32, fni_52, fne_32, fne_52,
     &  fchdia, fchin, fchvispar, fchvisper, fchvisq, fchinert,
     &  vaecrb, vadia, wadia, veecrb, vedia, floe_noc, floi_noc,
     1  b2fstati_version)
*   ..obtain external species data if any
      call ipgeti ('b2mndr_external_species', external_species)
      if (external_species.ne.0) then
        call read_b2mod_external(nx,ny,ns,ismain)
      else
        call alloc_b2mod_external_dummy(nx,ny)
      endif

*   conversion to new, unstructured data-types
      call calc_mapping (nx, ny, mpg)

      write (*,*) 'Mapping created'
      
      call prep_geometry (nx, ny,
     .        mpg%nCv, mpg%nCi, mpg%nFc, mpg%nVx, mpg%nFs, mpg%nFt,
     .        mpg%nVmxCv, gm, mpg)
      write (*,*) 'Geometry created'
      

      call prep_plasma (nx, ny, ns,
     .                    mpg%nCv, mpg%nCi, mpg%nFc, mpg, state)
      write (*,*) 'State created'

      call prep_plasma_ext(nx, ny, ns,
     .                    mpg%nCv, mpg%nCi, mpg%nFc, state_ext)
      write (*,*) 'State_ext created'

*   ..test plasma state and fluxes
!      call b2xvps (nx, ny, ns,
!     &  ne, na, ua, te, ti, po, fna, fhe, fhi, fch)
      call b2xvps (mpg%nCv, mpg%nFc, ns, state%pl, state%dv)

*   ..write headers for output files
      call cfverw (nout(1), newversion)
      idum(0) = mpg%nCi
      idum(1) = mpg%nCg
      idum(2) = mpg%nCv
      idum(3) = mpg%nFc
      idum(4) = mpg%nVx
      idum(5) = mpg%nFs
      idum(6) = mpg%nFt
      call cfwuin (nout(1),7,idum,'nCi,nCg,nCv,nFc,nVx,nFs,nFt')
      idum(0) = mpg%nCmxVx
      idum(1) = mpg%nCmxFc
      idum(2) = mpg%nVmxCv
      idum(3) = mpg%nVmxFc
      idum(4) = mpg%nCmxNv
      call cfwuin (nout(1),5,idum,'nCmxVx,nCmxFc,nVmxCv,nVmxFc,nCmxNv')
      idum(0) = mpg%isClassicalGrid
      call cfwuin (nout(1),1,idum,'isClassicalGrid')
      if (mpg%isClassicalGrid == 1) then
        idum(0) = mpg%nx
        idum(1) = mpg%ny
        idum(2) = mpg%nncut
      else
        ! dummy values for unstructured grid
        idum(0) = -1
        idum(1) = -1
        idum(2) = -1
      endif
      call cfwuin (nout(1), 3, idum, 'nx,ny,nncut')
      call cfwuch (nout(1), 120, lblgm, 'label')
      idum(0) = isymm
      call cfwuin (nout(1), 1, idum, 'isymm')
      idum(0) = mpg%nCv
      idum(1) = mpg%nFc
      idum(2) = ns
      call cfverw (nout(2), newversion)
      call cfwuin (nout(2), 3, idum, 'nCv,nFc,ns')
      call cfwuch (nout(2), 120, lblps, 'label')
      call b2wuzd (nout(2), newversion, ns, zamin, zamax, zn, am)

*   ..write state and geometry
      call write_map_geo(nout(1), mpg, gm)
      call write_plasma_state (nout(2), ns, mpg%nCv, mpg%nFc, state)

      call ipgeti ('b2stbc_boundary_namelist', boundary_namelist)
      call ipgeti ('b2stbr_neutrals_namelist', neutrals_namelist)
      call alloc_b2mod_elements(ns,zn,am)

      if (boundary_namelist.eq.1) then
        bcp_filename=repeat(' ',260)
        bcp_filename(1:22) = 'b2.boundary.parameters'
        call get_version(bcp_filename,boundary_version)
        if (boundary_version.lt.'03.002.000') then
          call alloc_b2mod_boundary(nx,ny,ns,nbcd,nnisod)
          call read_b2mod_boundary_namelist(nx, ny, ns, .false.)
          call conv_bcparms (nx, ny, ns, mpg, gm)
          call write_bcparms(ns,mpg)
        else
          call alloc_b2mod_boundary(mpg%nCi/2,mpg%nCi/2,ns,nbcd,nnisod)
          call read_b2mod_boundary_namelist(ns, mpg, .false.)
        end if
        call bc_to_struct(mpg)
      end if

#ifdef B25_EIRENE
      filename='input.dat'
      call find_file(filename,file_exists)
      if (file_exists.and.use_eirene.ne.0) then
        open(unit=45,file=filename)
        read(45,'(A72)',END=6999) line
        rewind(45)
        CALL EIRENE_FIND_PARAM(45)
        close(45)
        CALL EIRENE_SET_PARMMOD(1)
        write(*,*) 'Called FIND_PARAM'
        write(*,*) ' NSTRA = ',NSTRA
        write(*,*) '   NFL = ',NFL
        write(*,*) '  NATM = ',NATM
        write(*,*) '  NMOL = ',NMOL
        write(*,*) '  NION = ',NION
        write(*,*) 'NLIMPS = ',NLIMPS
        write(*,*) '  NLIM = ',NLIM
        write(*,*) '  NSTS = ',NSTS
      endif
      goto 7000
 6999 continue
      file_exist=.false.
 7000 continue
      if (.not.file_exists) then
        nstra=DEF_NSTRA
        nfl=DEF_NFL
        natm=DEF_NATM
        nmol=DEF_NMOL
        nion=DEF_NION
        npls=DEF_NPLS
        nlim=DEF_NLIM
        nsts=DEF_NSTS
        nlimps=nlim+nsts
!        call xerrab (' b2us -- no Eirene input available')
      endif
#endif 

      if (neutrals_namelist.eq.1 .or. use_eirene.ne.0) then
        bcp_filename=repeat(' ',260)
        bcp_filename(1:22) = 'b2.neutrals.parameters'
        call get_version(bcp_filename,neutrals_version)
        ntm = 0
        do is = 0, ns-1
          if (is_neutral(is)) ntm = ntm+1
        enddo
        natm_loc = max(natm, ntm)
        if (neutrals_version.lt.'03.002.000') then
          call alloc_b2mod_neutrals
     .     (nx, ny, ns, natm_loc, nmol, nion, nstraid, nsts, sput_dst)
          call read_neutrals_namelist(nx, ny, ns, .false.)
          call conv_rcparms (nx, ny, ns, mpg, gm)
          call write_rcparms(ns,mpg)
        else
          call alloc_b2mod_neutrals (mpg%nCi/2,mpg%nCi/2, ns, 
     .      natm_loc, nmol, nion, nstraid, nsts, sput_dst)
          call read_neutrals_namelist(ns, mpg, .false.)
        end if
        call rc_to_struct(mpg)
      end if
#ifdef B25_EIRENE
      open (unit=33)
      read (33,*,iostat=io_33) ico
      linq_33 = (io_33 == 0)
      close (33)
      open (unit=34)
      read (34,*,iostat=io_34) iel
      linq_34 = (io_34 == 0)
      close (34)
      open (unit=35)
      read (35,*,iostat=io_35) iel2
      linq_35 = (io_35 == 0)
      close (35)
      lcheck = (ico > 0) .and. (iel > 0) .and. (iel2 > 0) .and.
     .         (iel == iel2)
      if (linq_33 .and. linq_34 .and. linq_35 .and. lcheck) then
        write (*,*) ' files fort.33, fort.34 and fort.35 found '
        write (*,*) ' convert files to represent unstructured mesh '
        call conv_triangles(nout(0), mpg, gm, nx, ny)
      else
        write (*,*) ' no files fort.33, fort.34 and fort.35 available'
        write (*,*) ' create files fort.33, fort.34 and fort.35 '
        call write_triangles(nout(0), mpg, gm, nx, ny)
      end if
#else
!      if (use_eirene.ne.0) then
        call write_triangles(nout(0),mpg, gm, nx, ny)
!      endif
#endif
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
!WG_TODO      call dealloc_b2mod_sputter
      call dealloc_b2mod_indirect
      call dealloc_b2mod_residuals
      call dealloc_b2mod_transport
      call dealloc_b2mod_wall_init
      call dealloc_anomalous_transport
      call dealloc_b2mod_boundary
      call dealloc_b2mod_neutrals
      call dealloc_b2mod_neutrals_save
      call dealloc_geometry(gm)
      call dealloc_mapping(mpg)
      call destroyB2State(state)
      call prgend ()
#ifdef USE_MPI
      call mpi_barrier(MPI_COMM_WORLD,ierr_mpi)
      call mpi_finalize(ierr_mpi)
#endif
      stop 'b2us'

*-----------------------------------------------------------------------
*.error conditions

   91 call xerrab ('error trying to read internal parameters')

*-----------------------------------------------------------------------
*.end b2us

      end program b2us

