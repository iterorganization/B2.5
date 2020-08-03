      module b2mod_numerics_namelist
      use b2mod_types
      implicit none

#ifdef ALLOCATE_AND_NAMELIST
      real (kind=R8), allocatable, save ::
     *   dtco(:,:),
     *   dtmo(:,:),
     *   dtee(:),
     *   dtei(:)
     *   dten(:)
      logical, allocatable, save ::
     *   solveco(:,:),
     *   solvemo(:,:),
     *   solvemt(:),
     *   solvepo(:),
     *   solveee(:),
     *   solveei(:),
     *   solveen(:),
     *   solveet(:),
     *   last_solve_5(:),
     *   last_solve_9(:)
#else
#     include "DIMENSIONS.F"
      integer NSDMAX, CVREGMAX
      parameter (NSDMAX=DEF_NSD)
      parameter (CVREGMAX=8)
      real (kind=R8), save :: dtco(0:nsdmax-1,0:CVREGMAX)
      real (kind=R8), save :: dtmo(0:nsdmax-1,0:CVREGMAX)
      real (kind=R8), save :: dtee(0:CVREGMAX)
      real (kind=R8), save :: dtei(0:CVREGMAX)
      real (kind=R8), save :: dten(0:CVREGMAX)
      logical, save :: solveco(0:nsdmax-1,0:CVREGMAX)
      logical, save :: solvemo(0:nsdmax-1,0:CVREGMAX)
      logical, save :: solvemt(0:CVREGMAX)
      logical, save :: solvepo(0:CVREGMAX)
      logical, save :: solveee(0:CVREGMAX)
      logical, save :: solveei(0:CVREGMAX)
      logical, save :: solveen(0:CVREGMAX)
      logical, save :: solveet(0:CVREGMAX)
      logical, save :: last_solve_5(0:CVREGMAX)
      logical, save :: last_solve_9(0:CVREGMAX)
#endif
      real (kind=R8), save :: time_factor_required,
     1 core_dt_suppression, core_dt_factor
      real (kind=R8), allocatable, save :: time_factor(:)
      character*256, save :: numerics_filename
      real (kind=R8), save :: numerics_time_mod
      real (kind=R8), save :: numerics_time_switch
      logical :: first_time_step

      private :: numerics
      logical, save :: write_nml_num = .true.
      namelist /numerics/ dtco,dtmo,dtee,dtei,dten,time_factor_required,
     1 core_dt_suppression,core_dt_factor,
     2 solveco,solvemo,solvemt,solvepo,solveee,solveei,solveen,solveet,        !srv 22.05.18
     3 write_nml_num,numerics_filename,numerics_time_mod,
     $ numerics_time_switch


      contains

      subroutine alloc_b2mod_numerics(nCv,nsd,nnreg)
      implicit none
      integer nCv,nsd,nnreg(0:1), kul_recycle

      if (allocated(time_factor)) return
      write(*,*) 'Allocating space for b2mod_numerics'                !srv 25.02.06 11.06.08

      kul_recycle = 0    
      call ipgeti ('b2stbr_kul_recycle', kul_recycle)

      allocate(time_factor(nCv))
#ifdef ALLOCATE_AND_NAMELIST
      allocate(dtco(0:nsd-1,0:nnreg(0)))
      allocate(dtmo(0:nsd-1,0:nnreg(0)))
      allocate(dtee(0:nnreg(0)))
      allocate(dtei(0:nnreg(0)))
      allocate(dten(0:nnreg(0)))
      allocate(solveco(0:nsd-1,0:nnreg(0)))
      allocate(solvemo(0:nsd-1,0:nnreg(0)))
      allocate(solvemt(0:nnreg(0)))
      allocate(solvepo(0:nnreg(0)))
      allocate(solveee(0:nnreg(0)))
      allocate(solveei(0:nnreg(0)))
      allocate(solveen(0:nnreg(0)))
      allocate(solveet(0:nnreg(0)))
      allocate(last_solve_5(0:nnreg(0)))
      allocate(last_solve_9(0:nnreg(0)))
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
      time_factor_required=0.1_R8
      core_dt_suppression=1.0_R8
      core_dt_factor=1.0_R8
      solveco=.true.
      solvemo=.true.
      if (kul_recycle.ne.1) then !! default, use total momentum equation
        solvemt=.true.
      else
        solvemt=.false.          !! in case of KUL fluid neutrals, don't use total momentum equation
      endif
      solvepo=.true.
      solveee=.true.
      solveei=.true.
      solveen=.true.
      solveet=.true.
      last_solve_5=.true.
      last_solve_9=.true.
      numerics_filename='b2.numerics.parameters'
      numerics_time_mod=0.0_R8
      numerics_time_switch=0.0_R8
      time_factor=1.0_R8
      write(*,*) 'Allocated space for b2mod_numerics'                !srv 11.11.09

      return
      end subroutine alloc_b2mod_numerics

      subroutine dealloc_b2mod_numerics
      implicit none

      if (.not.allocated(time_factor)) return
      deallocate(time_factor)
#ifdef ALLOCATE_AND_NAMELIST
      deallocate(dtco)
      deallocate(dtmo)
      deallocate(dtee)
      deallocate(dtei)
      deallocate(dten)
      deallocate(solveco)
      deallocate(solvemo)
      deallocate(solvemt)
      deallocate(solvepo)
      deallocate(solveee)
      deallocate(solveei)
      deallocate(solveen)
      deallocate(solveet)
      deallocate(last_solve_5)
      deallocate(last_solve_9)
#endif
      return
      end subroutine dealloc_b2mod_numerics

      subroutine read_b2mod_numerics_namelist(nCv,nsmin,nsmax,nnreg,
     .  cvOnClosedSurface)
      implicit none

      integer nCv,nsmin,nsmax,nnreg(0:1)
      logical :: file_ok, cvOnClosedSurface(nCv)
      integer iCv, is, ireg
      real (kind=R8) :: ttf
      character*260 filename


      filename=numerics_filename
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
!        do ireg = 0, nnreg(0)
!          solveet(ireg) = solveee(ireg).or.solveei(ireg)
!          solvemt(ireg) = .false.
!          do is = nsmin, nsmax-1
!            solvemt(ireg) = solvemt(ireg) .or. solvemo(is,ireg)
!          end do
!        end do
        do ireg = 0, nnreg(0)                                     !srv 22.05.18 {
          solveet(ireg) = (solveee(ireg).and.solveei(ireg))
     &                                  .and.solveet(ireg)
          do is = nsmin, nsmax-1
            solvemt(ireg) = solvemt(ireg).and.solvemo(is,ireg)
          end do
        end do                                                    !srv 22.05.18 }
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

c additional time-modification factor
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
 999  call xerrab ('Error while reading NUMERICS namelist!')
      end subroutine read_b2mod_numerics_namelist

      subroutine write_b2mod_numerics_namelist
      implicit none

      if (write_nml_num) write(*,numerics)

      return
      end subroutine write_b2mod_numerics_namelist

      end module b2mod_numerics_namelist
