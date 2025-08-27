      module b2mod_ma28_for_11diag                               !srv 02.04.18 {
      use b2mod_types
      implicit none

      integer, allocatable, save :: iwrk_11p(:)
      real (kind=R8), allocatable, save ::
     1 aaf_11p(:), wrk0_11p(:), wrk1_11p(:), wrk2_11p(:)
#ifdef _OPENMP
!$OMP THREADPRIVATE (iwrk_11p, aaf_11p, wrk0_11p, wrk1_11p, wrk2_11p)
#endif

      real (kind=R8), allocatable, save ::
     * r1_11p(:),
     * s1_11p(:),
     * a1h_11p(:),
     * w_11p(:),
     * res_11p(:),
     * a1h_save_11p(:)

#ifdef _OPENMP
!$OMP THREADPRIVATE(r1_11p, s1_11p)
!$OMP THREADPRIVATE(a1h_11p, w_11p, res_11p, a1h_save_11p)
#endif

      integer, allocatable, save ::
     * irn_11p(:),
     * icn_11p(:),
     * ish_11p(:),
     * jsh_11p(:),
     * ijsh_11p(:),
     * irnh_11p(:),
     * icnh_11p(:),
     * icnhkeep_11p(:),
     * ikeep_11p(:,:),
     * iw28_11p(:,:),
     * irnh_save_11p(:),
     * icnh_save_11p(:)
#ifdef _OPENMP
!$OMP THREADPRIVATE(irn_11p, icn_11p, ish_11p, jsh_11p, ijsh_11p)
!$OMP THREADPRIVATE(irnh_11p, icnh_11p, icnhkeep_11p)
!$OMP THREADPRIVATE(ikeep_11p, iw28_11p, irnh_save_11p, icnh_save_11p)
#endif

      integer, save :: nc0_11p,nnd_11p
      integer, save :: n_nonzeros_11p,solvdim_11p,solvdim1_11p
      integer, save :: nltrsol_11p
      logical, save :: lfirst_ilu11 = .true.
#ifdef _OPENMP
!$OMP THREADPRIVATE(nc0_11p, nltrsol_11p, nnd_11p, n_nonzeros_11p)
!$OMP THREADPRIVATE(solvdim_11p, solvdim1_11p)
!$OMP THREADPRIVATE(lfirst_ilu11)
#endif

      integer, allocatable, save ::
     * indexes_11p(:,:), iasave_11p(:), jasave_11p(:)
#ifdef _OPENMP
!$OMP THREADPRIVATE(indexes_11p, iasave_11p, jasave_11p)
#endif

      contains

      subroutine alloc_b2mod_ma28_for_11diag
      implicit none

      if (.not.allocated(r1_11p)) then

        write(*,*) 'Allocating space for b2mod_ma28_for_11diag'
        write(*,*) 'nnd        = ',nnd_11p
        write(*,*) 'n_nonzeros = ',n_nonzeros_11p
        write(*,*) 'solvdim    = ',solvdim_11p
        write(*,*) 'solvdim1   = ',solvdim1_11p

        allocate(r1_11p(1:nnd_11p))
        allocate(s1_11p(1:nnd_11p))
        allocate(a1h_11p(1:solvdim_11p))
        allocate(w_11p(1:nnd_11p))
        allocate(res_11p(1:nnd_11p))
        allocate(a1h_save_11p(1:n_nonzeros_11p))
        allocate(irn_11p(1:n_nonzeros_11p))
        allocate(icn_11p(1:n_nonzeros_11p))
        allocate(ish_11p(1:n_nonzeros_11p))
        allocate(jsh_11p(1:n_nonzeros_11p))
        allocate(ijsh_11p(1:n_nonzeros_11p))
        allocate(irnh_11p(1:solvdim1_11p))
        allocate(icnh_11p(1:solvdim_11p))
        allocate(icnhkeep_11p(1:solvdim_11p))
        allocate(ikeep_11p(1:nnd_11p,5))
        allocate(iw28_11p(1:nnd_11p,8))
        allocate(irnh_save_11p(1:n_nonzeros_11p))
        allocate(icnh_save_11p(1:n_nonzeros_11p))

        write(*,*) 'Allocated space for b2mod_ma28_for_11diag'
      endif
      return
      end subroutine alloc_b2mod_ma28_for_11diag

      subroutine dealloc_b2mod_ma28_for_11diag

#ifdef _OPENMP
!$OMP PARALLEL
#endif
      if (allocated(r1_11p)) then
        deallocate(r1_11p)
        deallocate(s1_11p)
        deallocate(a1h_11p)
        deallocate(w_11p)
        deallocate(res_11p)
        deallocate(a1h_save_11p)
        deallocate(irn_11p)
        deallocate(icn_11p)
        deallocate(ish_11p)
        deallocate(jsh_11p)
        deallocate(ijsh_11p)
        deallocate(irnh_11p)
        deallocate(icnh_11p)
        deallocate(icnhkeep_11p)
        deallocate(ikeep_11p)
        deallocate(iw28_11p)
        deallocate(irnh_save_11p)
        deallocate(icnh_save_11p)
      endif
#ifdef _OPENMP
!$OMP END PARALLEL
#endif
      return
      end subroutine dealloc_b2mod_ma28_for_11diag

      subroutine restart_ma28_for_11
      implicit none
      lfirst_ilu11=.true.
      write(*,*) 'Restarting ma28_for_11'
      return
      end subroutine restart_ma28_for_11

      subroutine solver_init_11p(nx, ny)
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
      call ipgeti ('b2ux11p_style', style)
      call xertst (0.le.style.and.style.lt.3,
     &  'faulty internal parameter style')
      write(*,*)
     .  'Allocating space for solver for 11-point stencil matrix'  !srv 11.01.00
      if(style.eq.2) then
        mult_nonzero=15
        mult_solvdim=15
        mult_solvdim1=10
        call ipgeti ('b2ux11p_mult_nonzero', mult_nonzero)
        call ipgeti ('b2ux11p_mult_solvdim', mult_solvdim)
        call ipgeti ('b2ux11p_mult_solvdim1', mult_solvdim1)
        nnd_11p=(nx+2)*(ny+2)
        n_nonzeros_11p=mult_nonzero*nnd_11p
        solvdim_11p=mult_solvdim*n_nonzeros_11p
        solvdim1_11p=mult_solvdim1*n_nonzeros_11p
        call alloc_b2mod_ma28_for_11diag
        irnh_11p=0
        icnh_11p=0
        nltrsol_11p=2
        call ipgeti ('b2ux11p_nltrsol', nltrsol_11p)
        allocate(indexes_11p(1:3,1:solvdim_11p))
        allocate(iasave_11p(1:solvdim_11p),jasave_11p(1:solvdim_11p))
      else if (style.eq.0) then
        call xerrab ('Style not applicable to 11-point-stencil!')
*    ..obtain workspace for iluter
        allocate(aaf_11p(1:n2*5))
        allocate(wrk0_11p(1:n2))
        allocate(wrk1_11p(1:n2))
        allocate(wrk2_11p(1:n2))
        allocate(iwrk_11p(1:1))
      else if (style.eq.1) then
        call xerrab ('Style not applicable to 11-point-stencil!')
*    ..obtain workspace for slv5pt
        allocate(aaf_11p(1:n2*(3*min(nx+2,ny+2)+1)))
        allocate(wrk0_11p(1:n2))
        allocate(wrk1_11p(1:1))
        allocate(wrk2_11p(1:1))
        allocate(iwrk_11p(1:n2))
      endif
      write(*,*)
     .  'Allocated space for solver for 11-point stencil matrix'   !srv 11.01.00
#ifdef _OPENMP
!$OMP END PARALLEL
#endif
      return
      end subroutine solver_init_11p

      end module b2mod_ma28_for_11diag                           !srv 02.04.18 }
