      module b2mod_numerics_namelist
      use b2mod_types
      use b2mod_ad
     & , only : nsdmax, cvregmax
      implicit none

#ifdef ALLOCATE_AND_NAMELIST
      real (kind=R8), allocatable, save ::
     *   dtco(:,:),
     *   dtmo(:,:),
     *   dtee(:),
     *   dtei(:),
     *   dten(:),
     *   corr_core_dn(:)                                   ! IYS 19.12.2017
      logical, allocatable, save ::
     *   solveco(:,:),
     *   solvemo(:,:),
     *   solvemt(:),
     *   solvepo(:),
     *   solveee(:),
     *   solveei(:),
     *   solveen(:),
     *   solveet(:),
     *   solvekt(:),
     *   solvezt(:),
     *   last_solve_5(:),
     *   last_solve_9(:),
     *   add_te_corr_to_po(:)                            ! IYS 21.09.2018
#else
#     include "DIMENSIONS.F"
      real (kind=R8), save :: dtco(0:nsdmax-1,0:CVREGMAX)
      real (kind=R8), save :: dtmo(0:nsdmax-1,0:CVREGMAX)
      real (kind=R8), save :: dtee(0:CVREGMAX)
      real (kind=R8), save :: dtei(0:CVREGMAX)
      real (kind=R8), save :: dten(0:CVREGMAX)
      real (kind=R8), save :: corr_core_dn(0:nsdmax-1)     ! IYS 19.12.2017
      logical, save :: solveco(0:nsdmax-1,0:CVREGMAX)
      logical, save :: solvemo(0:nsdmax-1,0:CVREGMAX)
      logical, save :: solvemt(0:CVREGMAX)
      logical, save :: solvepo(0:CVREGMAX)
      logical, save :: solveee(0:CVREGMAX)
      logical, save :: solveei(0:CVREGMAX)
      logical, save :: solveen(0:CVREGMAX)
      logical, save :: solveet(0:CVREGMAX)
      logical, save :: solvekt(0:CVREGMAX)
      logical, save :: solvezt(0:CVREGMAX)
      logical, save :: last_solve_5(0:CVREGMAX)
      logical, save :: last_solve_9(0:CVREGMAX)
      logical, save :: add_te_corr_to_po(0:CVREGMAX)       ! IYS 21.09.2018
#endif
      real (kind=R8), save :: corr_core_dt                 ! IYS 19.12.2017
      real (kind=R8), save :: time_factor_required,
     1 core_dt_suppression, core_dt_factor
      real (kind=R8), allocatable, save :: time_factor(:)
      character*256, save :: numerics_filename
      real (kind=R8), save :: numerics_time_mod
      real (kind=R8), save :: numerics_time_switch
      logical :: first_time_step
      logical, save :: b2mod_numerics_allocated = .false.

      private :: numerics
      logical, save :: write_nml_num = .true.
      namelist /numerics/ dtco, dtmo, dtee, dtei, dten,
     & time_factor_required,
     & core_dt_suppression, core_dt_factor,
     & solveco, solvemo, solvemt, solvepo,
     & solveee, solveei, solveen, solveet,                 !srv 22.05.18
     & write_nml_num,
     & numerics_filename, numerics_time_mod, numerics_time_switch,
     & solvekt, solvezt,
     & corr_core_dn, corr_core_dt, add_te_corr_to_po       ! IYS 19.12.2017

      contains

      subroutine alloc_b2mod_numerics(nCv,nsd,nnreg,switch)
      use b2mod_switches
      implicit none
      integer nCv,nsd,nnreg(0:1)
      type (switches), intent (in) :: switch

      if (allocated(time_factor)) return
      write(*,*) 'Allocating space for b2mod_numerics'                !srv 25.02.06 11.06.08

      allocate(time_factor(nCv))
#ifdef ALLOCATE_AND_NAMELIST
      allocate(dtco(0:nsd-1,0:nnreg(0)))
      allocate(dtmo(0:nsd-1,0:nnreg(0)))
      allocate(dtee(0:nnreg(0)))
      allocate(dtei(0:nnreg(0)))
      allocate(dten(0:nnreg(0)))
      allocate(corr_core_dn(0:nsd-1))    ! IYS 19.12.2017
      allocate(solveco(0:nsd-1,0:nnreg(0)))
      allocate(solvemo(0:nsd-1,0:nnreg(0)))
      allocate(solvemt(0:nnreg(0)))
      allocate(solvepo(0:nnreg(0)))
      allocate(solveee(0:nnreg(0)))
      allocate(solveei(0:nnreg(0)))
      allocate(solveen(0:nnreg(0)))
      allocate(solveet(0:nnreg(0)))
      allocate(solvekt(0:nnreg(0)))
      allocate(solvezt(0:nnreg(0)))
      allocate(last_solve_5(0:nnreg(0)))
      allocate(last_solve_9(0:nnreg(0)))
      allocate(add_te_corr_to_po(0:nnreg(0)))                      ! IYS 21.09.2018, 21.09.2018
#else
      call xertst (nsd.le.nsdmax,
     1 'increase DEF_NSD in DIMENSIONS.F')
      call xertst (nnreg(0).le.cvregmax,'increase CVREGMAX')
#endif
      dtco=1.0_R8
      dtmo=1.0_R8
      dtee=1.0_R8
      dtei=1.0_R8
      dten=1.0_R8
      corr_core_dn=1.0_R8     ! IYS 19.12.2017
      corr_core_dt=1.0_R8     ! IYS 19.12.2017
      time_factor_required=0.1_R8
      core_dt_suppression=1.0_R8
      core_dt_factor=1.0_R8
      solveco=.true.
      solvemo=.true.
      if (switch%recycle_afn.ne.1) then !! default, use total momentum equation
        solvemt=.true.
      else
        solvemt=.false.                 !! in case of KUL fluid neutrals, don't use total momentum equation
      endif
      solvepo=.true.
      add_te_corr_to_po=.true.                                     ! IYS 21.09.2018
      solveee=.true.
      solveei=.true.
      solveen=.true.
      solveet=.true.
      solvekt=.true.
      solvezt=.true.
      if (switch%b2mndt_style.eq.2) then !! don't solve total momentum and energy equation in time-dependent mode
        solveet=.false.
        solvemt=.false.
      endif
      last_solve_5=.true.
      last_solve_9=.true.
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
      deallocate(dtco)
      deallocate(dtmo)
      deallocate(dtee)
      deallocate(dtei)
      deallocate(dten)
      deallocate(corr_core_dn)                                     ! IYS 19.12.2017
      deallocate(solveco)
      deallocate(solvemo)
      deallocate(solvemt)
      deallocate(solvepo)
      deallocate(add_te_corr_to_po)                                ! IYS 21.09.2018
      deallocate(solveee)
      deallocate(solveei)
      deallocate(solveen)
      deallocate(solveet)
      deallocate(solvekt)
      deallocate(solvezt)
      deallocate(last_solve_5)
      deallocate(last_solve_9)
#endif
      b2mod_numerics_allocated = .false.
      return
      end subroutine dealloc_b2mod_numerics

      subroutine read_b2mod_numerics_namelist(nCv,ns,nsmin,nsmax,nnreg,
     .  cvOnClosedSurface,b2mndt_style)
      implicit none

      integer nCv,ns,nsmin,nsmax,nnreg(0:1)
      logical :: file_ok, cvOnClosedSurface(nCv)
      integer iCv, is, ireg, b2mndt_style
      real (kind=R8) :: ttf
      character*260 filename

      filename=numerics_filename
      numerics_time_mod = 0.0_R8
      numerics_time_switch = 0.0_R8
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
        do ireg = 0, nnreg(0)                                     !srv 22.05.18 {
          solveet(ireg) = (solveee(ireg).and.solveei(ireg))
     &                                  .and.solveet(ireg)
          if (solveet(ireg).eq..true.) then
            call xerrab('time-dependent mode requires solveet=false')
          endif
          do is = nsmin, nsmax-1
            solvemt(ireg) = solvemt(ireg).and.solvemo(is,ireg)
          end do
          if (solvemt(ireg).eq..true.) then
            call xerrab('time-dependent mode requires solvemt=false')
          endif
        end do                                                    !srv 22.05.18 }
      else
        call xerrab ('No '//trim(filename))
      endif

      call xertst (0.0_R8.le.numerics_time_mod,
     . 'faulty input numerics_time_mod')
      call xertst (0.0_R8.le.numerics_time_switch,
     . 'faulty input numerics_time_switch')
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
      call xertst (0.0_R8.lt.minval(dtco),
     & 'faulty parameter dtco')
      call xertst (0.0_R8.lt.minval(dtmo),
     & 'faulty parameter dtmo')
      call xertst (0.0_R8.lt.minval(dtee),
     & 'faulty parameter dtee')
      call xertst (0.0_R8.lt.minval(dtei),
     & 'faulty parameter dtei')
      call xertst (0.0_R8.lt.minval(dten),
     & 'faulty parameter dten')

c additional time modification factor
      ttf=core_dt_suppression
      do iCv = 1,nCv
        ttf=min(ttf*core_dt_factor,1.0_R8)
        if(cvOnClosedSurface(iCv)) then
          time_factor(iCv)=ttf
        else
          time_factor(iCv)=1.0_R8
        endif
      enddo

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
