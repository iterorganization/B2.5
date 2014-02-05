      module b2mod_numerics_namelist
      use b2mod_types
      use b2mod_indirect
      implicit none

#ifdef ALLOCATE_AND_NAMELIST
      real (kind=R8), allocatable, save ::
     *   dtco(:,:), 
     *   dtmo(:,:), 
     *   dtee(:), 
     *   dtei(:)
      logical, allocatable, save ::
     *   solveco(:,:), 
     *   solvemo(:,:),
     *   solvepo(:),
     *   solveee(:),
     *   solveei(:)
#else
#     include "DIMENSIONS.F"
      integer NSDMAX
      parameter (NSDMAX=DEF_NSD)
      real (kind=R8), save :: dtco(0:nsdmax-1,0:NNREGMAX)
      real (kind=R8), save :: dtmo(0:nsdmax-1,0:NNREGMAX)
      real (kind=R8), save :: dtee(0:NNREGMAX)
      real (kind=R8), save :: dtei(0:NNREGMAX)
      logical, save :: solveco(0:nsdmax-1,0:NNREGMAX)
      logical, save :: solvemo(0:nsdmax-1,0:NNREGMAX)
      logical, save :: solvepo(0:NNREGMAX)
      logical, save :: solveee(0:NNREGMAX)
      logical, save :: solveei(0:NNREGMAX)
#endif
      real (kind=R8), save :: time_factor_required,
     1 core_dt_suppression, core_dt_factor
      real (kind=R8), allocatable, save :: time_factor(:,:)
      character*256, save :: numerics_filename
      real (kind=R8), save :: numerics_time_mod
      real (kind=R8), save :: numerics_time_switch

      contains

      subroutine alloc_b2mod_numerics(nsd)
      implicit none
      integer nsd
      write(*,*) 'Allocating space for b2mod_numerics'                !srv 25.02.06 11.06.08

#ifdef ALLOCATE_AND_NAMELIST
      allocate(dtco(0:nsd-1,0:nnreg(0)))
      allocate(dtmo(0:nsd-1,0:nnreg(0)))
      allocate(dtee(0:nnreg(0)))
      allocate(dtei(0:nnreg(0)))
      allocate(solveco(0:nsd-1,0:nnreg(0)))
      allocate(solvemo(0:nsd-1,0:nnreg(0)))
      allocate(solvepo(0:nnreg(0)))
      allocate(solveee(0:nnreg(0)))
      allocate(solveei(0:nnreg(0)))
#else
      if(nsd.gt.nsdmax)
     1 stop 'increase NSDMAX in b2mod_numerics_namelist'
      if(nnreg(0).gt.nnregmax)
     1 stop 'increase NNREGMAX in b2mod_numerics_namelist'
#endif
      dtco=1.0_R8
      dtmo=1.0_R8
      dtee=1.0_R8
      dtei=1.0_R8
      time_factor_required=0.1_R8
      core_dt_suppression=1.0_R8
      core_dt_factor=1.0_R8
      solveco=.true.
      solvemo=.true.
      solvepo=.true.
      solveee=.true.
      solveei=.true.
      numerics_filename='b2.numerics.parameters'
      numerics_time_mod=0.0_R8
      numerics_time_switch=0.0_R8
      write(*,*) 'Allocated space for b2mod_numerics'                !srv 11.11.09

      return
      end subroutine alloc_b2mod_numerics

      subroutine read_b2mod_numerics_namelist
      implicit none

      namelist /numerics/ dtco,dtmo,dtee,dtei,time_factor_required,
     1 core_dt_suppression,core_dt_factor,
     2 solveco,solvemo,solvepo,solveee,solveei,
     3 numerics_filename,numerics_time_mod,numerics_time_switch
      logical :: file_ok
      character*260 filename

      filename='../' // numerics_filename
      inquire(file=filename,exist=file_ok)
      if(file_ok) then
        open(99,file=filename)
        write(*,*) 'Reading namelist numerics from ',trim(filename)
        read(99,numerics)
        write(*,*) 'Read namelist numerics from ',trim(filename)
      else
        write(*,*) 'No ',trim(filename)
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

      return
      end subroutine read_b2mod_numerics_namelist

      subroutine write_b2mod_numerics_namelist
      implicit none

      namelist /numerics/ dtco,dtmo,dtee,dtei,time_factor_required,
     1 core_dt_suppression,core_dt_factor,
     2 solveco,solvemo,solvepo,solveee,solveei,
     3 numerics_filename,numerics_time_mod,numerics_time_switch

      write(*,numerics)

      return
      end subroutine write_b2mod_numerics_namelist

      end module b2mod_numerics_namelist
