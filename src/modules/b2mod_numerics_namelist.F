      module b2mod_numerics_namelist
      use b2mod_types
      use b2mod_indirect
      implicit none

#ifdef ALLOCATE_AND_NAMELIST
      real (kind=R8), allocatable, save ::
     *   dtco(:,:), dtco_SOL(:),
     *   dtmo(:,:), dtmo_SOL(:),
     *   dtee(:), dtee_SOL,
     *   dtei(:), dtei_SOL
     *   , corr_core_dn(:)                                 ! IYS 19.12.2017
      logical, allocatable, save ::
     *   solveco(:,:),
     *   solvemo(:,:),
     *   solvemt(:),
     *   solvepo(:),
     *   solveee(:),
     *   solveei(:),
     *   solveet(:),
     *   last_solve_5(:),
     *   last_solve_7(:),                                  !srv 02.04.18
     *   last_solve_9(:),
     *   last_solve_11(:),                                 !srv 02.04.18
     *   add_te_corr_to_po(:)  ! IYS 21.09.2018
#else
#     include "DIMENSIONS.F"
      integer NSDMAX
      parameter (NSDMAX=DEF_NSD)
      real (kind=R8), save :: 
     &    dtco(0:nsdmax-1,0:NNREGMAX), dtco_SOL(0:nsdmax-1)
      real (kind=R8), save :: 
     &    dtmo(0:nsdmax-1,0:NNREGMAX), dtmo_SOL(0:nsdmax-1)
      real (kind=R8), save :: dtee(0:NNREGMAX), dtee_SOL
      real (kind=R8), save :: dtei(0:NNREGMAX), dtei_SOL
      real (kind=R8), save :: corr_core_dn(0:nsdmax-1)     ! IYS 19.12.2017
      logical, save :: solveco(0:nsdmax-1,0:NNREGMAX)
      logical, save :: solvemo(0:nsdmax-1,0:NNREGMAX)
      logical, save :: solvemt(0:NNREGMAX)
      logical, save :: solvepo(0:NNREGMAX)
      logical, save :: solveee(0:NNREGMAX)
      logical, save :: solveei(0:NNREGMAX)
      logical, save :: solveet(0:NNREGMAX)
      logical, save :: last_solve_5(0:NNREGMAX)
      logical, save :: last_solve_7(0:NNREGMAX)            !srv 02.04.18
      logical, save :: last_solve_9(0:NNREGMAX)
      logical, save :: last_solve_11(0:NNREGMAX)           !srv 02.04.18
      logical, save :: add_te_corr_to_po(0:NNREGMAX)       ! IYS 21.09.2018
#endif
      !$OMP THREADPRIVATE(last_solve_5, last_solve_7)
      !$OMP THREADPRIVATE(last_solve_9, last_solve_11)
      real (kind=R8), save :: corr_core_dt                 ! IYS 19.12.2017
      real (kind=R8), save :: dteff
      real (kind=R8), save :: time_factor_required,
     1 core_dt_suppression, core_dt_factor
      real (kind=R8), allocatable, save :: time_factor(:,:)
      character*256, save :: numerics_filename
      real (kind=R8), save :: numerics_time_mod
      real (kind=R8), save :: numerics_time_switch
      logical :: first_time_step
      logical, save :: b2mod_numerics_allocated = .false.
      integer, save :: iy_SOL_dtXX = 1

      private :: numerics
      logical, save :: write_nml_num = .true.
      namelist /numerics/ dtco, dtmo, dtee, dtei,
     & dtco_SOL, dtmo_SOL, dtee_SOL, dtei_SOL, iy_SOL_dtXX,
     & time_factor_required,
     & core_dt_suppression, core_dt_factor,
     & solveco, solvemo, solvemt, solvepo, solveee, solveei, solveet,       !srv 22.05.18
     & write_nml_num,
     & numerics_filename, numerics_time_mod, numerics_time_switch,
     & corr_core_dn, corr_core_dt, add_te_corr_to_po         ! IYS 19.12.2017

      contains

      subroutine alloc_b2mod_numerics(nxd,nyd,nsd)
      implicit none
      integer nxd,nyd,nsd

      if (b2mod_numerics_allocated) return
      write(*,*) 'Allocating space for b2mod_numerics'                !srv 25.02.06 11.06.08

      allocate(time_factor(-1:nxd,-1:nyd))
#ifdef ALLOCATE_AND_NAMELIST
      allocate(dtco(0:nsd-1,0:nnreg(0)))
      allocate(dtco_SOL(0:nsd-1))
      allocate(dtmo(0:nsd-1,0:nnreg(0)))
      allocate(dtmo_SOL(0:nsd-1))
      allocate(dtee(0:nnreg(0)))
      allocate(dtee_SOL)
      allocate(dtei_SOL)
      allocate(corr_core_dn(0:nsd-1))    ! IYS 19.12.2017
      allocate(solveco(0:nsd-1,0:nnreg(0)))
      allocate(solvemo(0:nsd-1,0:nnreg(0)))
      allocate(solvemt(0:nnreg(0)))
      allocate(solvepo(0:nnreg(0)))
      allocate(solveee(0:nnreg(0)))
      allocate(solveei(0:nnreg(0)))
      allocate(solveet(0:nnreg(0)))
      allocate(last_solve_5(0:nnreg(0)))
      allocate(last_solve_7(0:nnreg(0)))                           !srv 02.04.18
      allocate(last_solve_9(0:nnreg(0)))
      allocate(last_solve_11(0:nnreg(0)))                          !srv 02.04.18
      allocate(add_te_corr_to_po(0:nnreg(0)))                      ! IYS 21.09.2018, 21.09.2018
#else
      call xertst (nsd.le.nsdmax,
     1 'increase DEF_NSD in DIMENSIONS.F')
      call xertst (nnreg(0).le.nnregmax,
     1 'increase NNREGMAX in b2mod_indirect')
#endif
      dtco=1.0_R8
      dtmo=1.0_R8
      dtee=1.0_R8
      dtei=1.0_R8
      dtco_SOL=1.0_R8
      dtmo_SOL=1.0_R8
      dtee_SOL=1.0_R8
      dtei_SOL=1.0_R8
      corr_core_dn=1.0_R8     ! IYS 19.12.2017
      corr_core_dt=1.0_R8     ! IYS 19.12.2017
      time_factor_required=0.1_R8
      core_dt_suppression=1.0_R8
      core_dt_factor=1.0_R8
      solveco=.true.
      solvemo=.true.
      solvemt=.true.
      solvepo=.true.
      add_te_corr_to_po=.true.                                     ! IYS 21.09.2018
      solveee=.true.
      solveei=.true.
      solveet=.true.
      last_solve_5=.true.
      last_solve_7=.true.                                          !srv 02.04.18
      last_solve_9=.true.
      last_solve_11=.true.                                         !srv 02.04.18
      numerics_filename='b2.numerics.parameters'
      numerics_time_mod=0.0_R8
      numerics_time_switch=0.0_R8
      time_factor=1.0_R8
      write(*,*) 'Allocated space for b2mod_numerics'              !srv 11.11.09
      b2mod_numerics_allocated = .true.

      return
      end subroutine alloc_b2mod_numerics

      subroutine dealloc_b2mod_numerics
      implicit none

      if (.not.b2mod_numerics_allocated) return
      deallocate(time_factor)
#ifdef ALLOCATE_AND_NAMELIST
      deallocate(dtco, dtco_SOL)
      deallocate(dtmo, dtmo_SOL)
      deallocate(dtee, dtee_SOL)
      deallocate(dtei, dtei_SOL)
      deallocate(corr_core_dn)                                     ! IYS 19.12.2017
      deallocate(solveco)
      deallocate(solvemo)
      deallocate(solvemt)
      deallocate(solvepo)
      deallocate(add_te_corr_to_po)                                ! IYS 21.09.2018
      deallocate(solveee)
      deallocate(solveei)
      deallocate(solveet)
      deallocate(last_solve_5)
      deallocate(last_solve_7)                                     !srv 02.04.18
      deallocate(last_solve_9)
      deallocate(last_solve_11)                                    !srv 02.04.18
#endif
      b2mod_numerics_allocated = .false.
      return
      end subroutine dealloc_b2mod_numerics

      subroutine read_b2mod_numerics_namelist(nx,ny,ns,nsmin,nsmax)
      implicit none

      integer nx,ny,ns,nsmin,nsmax
      logical :: file_ok
      integer ix, iy, ireg, is
      integer, save :: style_b2mndt = 1
      integer, save :: ncall = 0
      real(kind=R8), save :: dtim = 1.0_R8
      real (kind=R8) :: ttf
      character*260 filename

      filename=numerics_filename
      numerics_time_mod = 0.0_R8
      numerics_time_switch = 0.0_R8
      if (ncall.eq.0) then
        call ipgeti('b2mndt_style', style_b2mndt)
        call ipgetr('b2mndr_dtim', dtim)
        dteff = dtim
      end if
      call find_file(filename,file_ok)
      if(file_ok) then
        open(99,file=filename)
        write(*,*) 'Reading namelist numerics from ',trim(filename)
#ifdef USE_MPI
        read(99,numerics,err=999)
#else
        read(99,numerics)
#endif
        close(99)
        write(*,*) 'Read namelist numerics from ',trim(filename)
        do ireg = 0, nnreg(0)
          solveet(ireg) = (solveee(ireg).and.solveei(ireg))
     &                                  .and.solveet(ireg)
          call xertst(style_b2mndt.le.1.or..not.solveet(ireg),
     &     'Time-dependent mode requires SOLVEET = .false.!')
          do is = nsmin, nsmax-1
            solvemt(ireg) = solvemt(ireg).and.solvemo(is,ireg)
          end do
          call xertst(style_b2mndt.le.1.or..not.solvemt(ireg),
     &     'Time-dependent mode requires SOLVEMT = .false.!')
        end do
      else
        call xerrab ('No '//trim(filename))
      endif

      call xertst (0.0_R8.le.numerics_time_mod,
     . 'faulty input boundary_time_mod')
      call xertst (0.0_R8.le.numerics_time_switch,
     . 'faulty input boundary_time_switch')
      call xertst (0.0_R8.lt.core_dt_suppression,
     & 'faulty parameter core_dt_suppression')
      call xertst (0.0_R8.lt.core_dt_factor,
     & 'faulty parameter core_dt_factor')
      call xertst (0.0_R8.lt.corr_core_dt,
     & 'faulty parameter corr_core_dt')
      do is = 0, ns-1
        call xertst (0.0_R8.lt.corr_core_dn(is),
     &   'faulty parameter corr_core_dn')
      end do
      call xertst (0.0_R8.le.time_factor_required,
     & 'faulty parameter time_factor_required')
      call xertst (0.0_R8.lt.min(minval(dtco), minval(dtco_SOL)),
     & 'faulty parameter dtco')
      call xertst (0.0_R8.lt.min(minval(dtmo), minval(dtmo_SOL)),
     & 'faulty parameter dtmo')
      call xertst (0.0_R8.lt.min(minval(dtee), dtee_SOL),
     & 'faulty parameter dtee')
      call xertst (0.0_R8.lt.min(minval(dtei), dtei_SOL),
     & 'faulty parameter dtei')

      ttf = min(minval(dtco(nsmin:nsmax-1,1:nnreg(0))),
     &                        minval(dtco_SOL(nsmin:nsmax-1)))
      if (style_b2mndt.ne.2.and.
     &    ttf.eq.max(maxval(dtco(nsmin:nsmax-1,1:nnreg(0))),
     &                        maxval(dtco_SOL(nsmin:nsmax-1))).and.
     &    ttf.eq.min(minval(dtmo(nsmin:nsmax-1,1:nnreg(0))),
     &                        minval(dtmo_SOL(nsmin:nsmax-1))).and.
     &    ttf.eq.max(maxval(dtmo(nsmin:nsmax-1,1:nnreg(0))),
     &                        maxval(dtmo_SOL(nsmin:nsmax-1))).and.
     &    ttf.eq.min(minval(dtee(1:nnreg(0))),dtee_SOL).and.
     &    ttf.eq.max(maxval(dtee(1:nnreg(0))),dtee_SOL).and.
     &    ttf.eq.min(minval(dtei(1:nnreg(0))),dtei_SOL).and.
     &    ttf.eq.max(maxval(dtei(1:nnreg(0))),dtei_SOL)) then
        dteff = dtim*ttf
        if (ttf.ne.1.0_R8) write(0,*)
     &   'Effective plasma timestep rescaled to ', dteff,' seconds'
      end if

c additional time modification factor
      ttf=core_dt_suppression
      do iy=-1,ny
        if(iy.gt.0) ttf=min(ttf*core_dt_factor,1.0_R8)
        do ix=-1,nx
          if(on_closed_surface(ix,iy)) then
            time_factor(ix,iy)=ttf
          else
            time_factor(ix,iy)=1.0_R8
          endif
        enddo
      enddo

      ncall = ncall + 1
      return
#ifdef USE_MPI
  999 call xerrab ('Error while reading NUMERICS namelist!')
#endif
      end subroutine read_b2mod_numerics_namelist

      subroutine write_b2mod_numerics_namelist
      implicit none

      if (write_nml_num) write(*,numerics)

      return
      end subroutine write_b2mod_numerics_namelist

      end module b2mod_numerics_namelist
