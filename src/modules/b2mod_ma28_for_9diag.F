      module b2mod_ma28_for_9diag                                !sv 13.10.06 {
      use b2mod_types
      implicit none

      integer, allocatable, save :: iwrk_9p(:)
      real (kind=R8), allocatable, save ::
     1 aaf_9p(:), wrk0_9p(:), wrk1_9p(:), wrk2_9p(:)
#ifdef _OPENMP
!$OMP THREADPRIVATE (iwrk_9p, aaf_9p, wrk0_9p, wrk1_9p, wrk2_9p)
#endif

      real (kind=R8), allocatable, save ::
     * r1_9p(:),
     * s1_9p(:),
     * a1h_9p(:),
     * w_9p(:),
     * res_9p(:),
     * a1h_save_9p(:)

#ifdef _OPENMP
!$OMP THREADPRIVATE(r1_9p, s1_9p, a1h_9p, w_9p, res_9p, a1h_save_9p)
#endif

      integer, allocatable, save ::
     * irn_9p(:),
     * icn_9p(:),
     * ish_9p(:),
     * jsh_9p(:),
     * ijsh_9p(:),
     * irnh_9p(:),
     * icnh_9p(:),
     * icnhkeep_9p(:),
     * ikeep_9p(:,:),
     * iw28_9p(:,:),
     * irnh_save_9p(:),
     * icnh_save_9p(:)
#ifdef _OPENMP
!$OMP THREADPRIVATE(irn_9p, icn_9p, ish_9p, jsh_9p, ijsh_9p)
!$OMP THREADPRIVATE(irnh_9p, icnh_9p, icnhkeep_9p)
!$OMP THREADPRIVATE(ikeep_9p, iw28_9p, irnh_save_9p, icnh_save_9p)
#endif

      integer, save :: nc0_9p,nnd_9p
      integer, save :: n_nonzeros_9p,solvdim_9p,solvdim1_9p
      integer, save :: nltrsol_9p
      logical, save :: lfirst_ilu9 = .true.
#ifdef _OPENMP
!$OMP THREADPRIVATE(nc0_9p, nltrsol_9p, nnd_9p, n_nonzeros_9p)
!$OMP THREADPRIVATE(solvdim_9p, solvdim1_9p)
!$OMP THREADPRIVATE(lfirst_ilu9)
#endif

      integer, allocatable, save ::
     * indexes_9p(:,:), iasave_9p(:), jasave_9p(:)
#ifdef _OPENMP
!$OMP THREADPRIVATE(indexes_9p, iasave_9p, jasave_9p)
#endif

      contains

      subroutine alloc_b2mod_ma28_for_9diag
      implicit none

      if (.not.allocated(r1_9p)) then

        write(*,*) 'Allocating space for b2mod_ma28_for_9diag'
        write(*,*) 'nnd        = ',nnd_9p
        write(*,*) 'n_nonzeros = ',n_nonzeros_9p
        write(*,*) 'solvdim    = ',solvdim_9p
        write(*,*) 'solvdim1   = ',solvdim1_9p

        allocate(r1_9p(1:nnd_9p))
        allocate(s1_9p(1:nnd_9p))
        allocate(a1h_9p(1:solvdim_9p))
        allocate(w_9p(1:nnd_9p))
        allocate(res_9p(1:nnd_9p))
        allocate(a1h_save_9p(1:n_nonzeros_9p))
        allocate(irn_9p(1:n_nonzeros_9p))
        allocate(icn_9p(1:n_nonzeros_9p))
        allocate(ish_9p(1:n_nonzeros_9p))
        allocate(jsh_9p(1:n_nonzeros_9p))
        allocate(ijsh_9p(1:n_nonzeros_9p))
        allocate(irnh_9p(1:solvdim1_9p))
        allocate(icnh_9p(1:solvdim_9p))
        allocate(icnhkeep_9p(1:solvdim_9p))
        allocate(ikeep_9p(1:nnd_9p,5))
        allocate(iw28_9p(1:nnd_9p,8))
        allocate(irnh_save_9p(1:n_nonzeros_9p))
        allocate(icnh_save_9p(1:n_nonzeros_9p))

        write(*,*) 'Allocated space for b2mod_ma28_for_9diag'
      end if
      return
      end subroutine alloc_b2mod_ma28_for_9diag

      subroutine dealloc_b2mod_ma28_for_9diag

#ifdef _OPENMP
!$OMP PARALLEL
#endif
      if (allocated(r1_9p)) then
        deallocate(r1_9p)
        deallocate(s1_9p)
        deallocate(a1h_9p)
        deallocate(w_9p)
        deallocate(res_9p)
        deallocate(a1h_save_9p)
        deallocate(irn_9p)
        deallocate(icn_9p)
        deallocate(ish_9p)
        deallocate(jsh_9p)
        deallocate(ijsh_9p)
        deallocate(irnh_9p)
        deallocate(icnh_9p)
        deallocate(icnhkeep_9p)
        deallocate(ikeep_9p)
        deallocate(iw28_9p)
        deallocate(irnh_save_9p)
        deallocate(icnh_save_9p)
      endif
#ifdef _OPENMP
!$OMP END PARALLEL
#endif
      return
      end subroutine dealloc_b2mod_ma28_for_9diag

      subroutine restart_ma28_for_9
      implicit none
      lfirst_ilu9=.true.
      write(*,*) 'Restarting ma28_for_9'
      return
      end subroutine restart_ma28_for_9

      subroutine solver_init_9p(nx, ny)
      use b2mod_indirect
      implicit none
      integer, intent(in) :: nx, ny
      integer, save :: style
      integer, save :: mult_nonzero, mult_solvdim, mult_solvdim1
      integer :: n2
      external ipgeti, xertst
      data style/2/

cxpb This routine to be called once to initialize the arrays used by
cxpb the matrix solver. Must be done once only and outside the iteration
cxpb loop to avoid errors when using OpenMP parallelization.

*   ..set n2
      n2 = (nx+2)*(ny+2)
*   ..test nx, ny, itcnt
      call xertst (0.le.nx.and.0.le.ny, 'faulty argument nx, ny')
*   ..alocate arrays
#ifdef _OPENMP
!$OMP PARALLEL
!$OMP& SHARED(nx, ny, n2)
#endif
      call ipgeti ('b2ux9p_style', style)
      call xertst (0.le.style.and.style.lt.3,
     &  'faulty internal parameter style')
      write(*,*)
     .  'Allocating space for solver for 9-point stencil matrix'  !srv 11.01.00
      if(style.eq.2) then
        mult_nonzero=10
        mult_solvdim=15
        mult_solvdim1=10
        call ipgeti ('b2ux9p_mult_nonzero', mult_nonzero)
        call ipgeti ('b2ux9p_mult_solvdim', mult_solvdim)
        call ipgeti ('b2ux9p_mult_solvdim1', mult_solvdim1)
        nnd_9p=(nx+2)*(ny+2)
        n_nonzeros_9p=mult_nonzero*nnd_9p
        solvdim_9p=mult_solvdim*n_nonzeros_9p
        solvdim1_9p=mult_solvdim1*n_nonzeros_9p
        call alloc_b2mod_ma28_for_9diag
        irnh_9p=0
        icnh_9p=0
        nltrsol_9p=2
        call ipgeti ('b2ux9p_nltrsol', nltrsol_9p)
        allocate(indexes_9p(1:3,1:solvdim_9p))
        allocate(iasave_9p(1:solvdim_9p),jasave_9p(1:solvdim_9p))
      else if (style.eq.0) then
        call xerrab ('b2ux9p: Style 0 is not implemented')
*    ..obtain workspace for iluter
        allocate(aaf_9p(1:n2*5))
        allocate(wrk0_9p(1:n2))
        allocate(wrk1_9p(1:n2))
        allocate(wrk2_9p(1:n2))
        allocate(iwrk_9p(1:1))
      else if (style.eq.1) then
        call xerrab ('b2ux9p: Style 1 is not implemented')
*    ..obtain workspace for slv5pt
        allocate(aaf_9p(1:n2*(3*min(nx+2,ny+2)+1)))
        allocate(wrk0_9p(1:n2))
        allocate(wrk1_9p(1:1))
        allocate(wrk2_9p(1:1))
        allocate(iwrk_9p(1:n2))
      endif
      write(*,*)
     .  'Allocated space for solver for 9-point stencil matrix'   !srv 11.01.00
#ifdef _OPENMP
!$OMP END PARALLEL
#endif
      return
      end subroutine solver_init_9p

      end module b2mod_ma28_for_9diag                             !sv 13.10.06 }
