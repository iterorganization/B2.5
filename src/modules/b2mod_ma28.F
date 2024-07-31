      module b2mod_ma28
      use b2mod_types
      implicit none

      integer, allocatable, save :: iwrk_5p(:)
      real (kind=R8), allocatable, save ::
     1 aaf_5p(:), wrk0_5p(:), wrk1_5p(:), wrk2_5p(:)

      real (kind=R8), allocatable, save ::
     * r1_5p(:),
     * s1_5p(:),
     * a1h_5p(:),
     * w_5p(:),
     * res_5p(:),
     * a1h_save_5p(:)

      integer, allocatable, save ::
     * irn_5p(:),
     * icn_5p(:),
     * ish_5p(:),
     * jsh_5p(:),
     * ijsh_5p(:),
     * irnh_5p(:),
     * icnh_5p(:),
     * icnhkeep_5p(:),
     * ikeep_5p(:,:),
     * iw28_5p(:,:),
     * irnh_save_5p(:),
     * icnh_save_5p(:)

      integer, save :: nc0_5p,nnd_5p
      integer, save :: n_nonzeros_5p,solvdim_5p,solvdim1_5p
      integer, save :: nltrsol_5p
      logical, save :: lfirst_ilu = .true.

#ifdef _OPENMP
!!FR All the allocatable arrays need to be declared as THREADPRIVATE
!!FR as each thread needs its own copy of the arrays used in the solver.
!!FR All the integers also need to be threadprivate as they are used
!!FR in the ncall==0 block in b2ux5p that defines the array sizes.
!$OMP THREADPRIVATE (iwrk_5p, aaf_5p, wrk0_5p, wrk1_5p, wrk2_5p)
!$OMP THREADPRIVATE (r1_5p, s1_5p, a1h_5p, w_5p, res_5p, a1h_save_5p)
!$OMP THREADPRIVATE (irn_5p, icn_5p, irnh_5p, icnh_5p)
!$OMP THREADPRIVATE (ish_5p, jsh_5p, ijsh_5p, iw28_5p)
!$OMP THREADPRIVATE (icnhkeep_5p, ikeep_5p, irnh_save_5p, icnh_save_5p)
!$OMP THREADPRIVATE (nc0_5p, nltrsol_5p, nnd_5p, n_nonzeros_5p)
!$OMP THREADPRIVATE (solvdim_5p, solvdim1_5p)
!$OMP THREADPRIVATE (lfirst_ilu)
#endif

      contains

      subroutine alloc_b2mod_ma28
      implicit none

      if (.not.allocated(r1_5p)) then

        write(*,*) 'Allocating space for b2mod_ma28'
        write(*,*) 'nnd = ',nnd_5p
        write(*,*) 'n_nonzeros = ',n_nonzeros_5p
        write(*,*) 'solvdim = ',solvdim_5p
        write(*,*) 'solvdim1 = ',solvdim1_5p

        allocate(r1_5p(1:nnd_5p))
        allocate(s1_5p(1:nnd_5p))
        allocate(a1h_5p(1:solvdim_5p))
        allocate(w_5p(1:nnd_5p))
        allocate(res_5p(1:nnd_5p))
        allocate(a1h_save_5p(1:n_nonzeros_5p))
        allocate(irn_5p(1:n_nonzeros_5p))
        allocate(icn_5p(1:n_nonzeros_5p))
        allocate(ish_5p(1:n_nonzeros_5p))
        allocate(jsh_5p(1:n_nonzeros_5p))
        allocate(ijsh_5p(1:n_nonzeros_5p))
        allocate(irnh_5p(1:solvdim1_5p))
        allocate(icnh_5p(1:solvdim_5p))
        allocate(icnhkeep_5p(1:solvdim_5p))
        allocate(ikeep_5p(1:nnd_5p,5))
        allocate(iw28_5p(1:nnd_5p,8))
        allocate(irnh_save_5p(1:n_nonzeros_5p))
        allocate(icnh_save_5p(1:n_nonzeros_5p))

        write(*,*) 'Allocated space for b2mod_ma28'
      end if
      return
      end subroutine alloc_b2mod_ma28

      subroutine dealloc_b2mod_ma28
#ifdef _OPENMP
!$OMP PARALLEL
#endif
      if (allocated(r1_5p)) then
        deallocate(r1_5p)
        deallocate(s1_5p)
        deallocate(a1h_5p)
        deallocate(w_5p)
        deallocate(res_5p)
        deallocate(a1h_save_5p)
        deallocate(irn_5p)
        deallocate(icn_5p)
        deallocate(ish_5p)
        deallocate(jsh_5p)
        deallocate(ijsh_5p)
        deallocate(irnh_5p)
        deallocate(icnh_5p)
        deallocate(icnhkeep_5p)
        deallocate(ikeep_5p)
        deallocate(iw28_5p)
        deallocate(irnh_save_5p)
        deallocate(icnh_save_5p)
      endif
#ifdef _OPENMP
!$OMP END PARALLEL
#endif
      return
      end subroutine dealloc_b2mod_ma28

      subroutine restart_ma28
      implicit none
      lfirst_ilu=.true.
      write(*,*) 'Restarting ma28'
      return
      end subroutine restart_ma28

      subroutine solver_init_5p(nx, ny)
      use b2mod_indirect
      implicit none
      integer, intent(in) :: nx, ny
      integer, save :: style
      integer, save :: mult_nonzero, mult_solvdim, mult_solvdim1
      integer :: n2, ix, iy
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
!$OMP& PRIVATE(ix, iy)
#endif
      call ipgeti ('b2ux5p_style', style)
      call xertst (0.le.style.and.style.lt.3,
     &  'faulty internal parameter style')
      write(*,*) 'Allocating space for solver'
      if(style.eq.2) then
        mult_nonzero=10
        mult_solvdim=15
        mult_solvdim1=10
        call ipgeti ('b2ux5p_mult_nonzero', mult_nonzero)
        call ipgeti ('b2ux5p_mult_solvdim', mult_solvdim)
        call ipgeti ('b2ux5p_mult_solvdim1', mult_solvdim1)
        nnd_5p=(nx+2)*(ny+2)
        n_nonzeros_5p=mult_nonzero*nnd_5p
        solvdim_5p=mult_solvdim*n_nonzeros_5p
        solvdim1_5p=mult_solvdim1*n_nonzeros_5p
        call alloc_b2mod_ma28
        call iniilu(nx+2,ny+2)
        nltrsol_5p=2
        call ipgeti ('b2ux5p_nltrsol', nltrsol_5p)
      elseif (style.eq.0) then
*    ..obtain workspace for iluter
        allocate(aaf_5p(1:n2*5))
        allocate(wrk0_5p(1:n2))
        allocate(wrk1_5p(1:n2))
        allocate(wrk2_5p(1:n2))
        allocate(iwrk_5p(1:1))
      else if (style.eq.1) then
*    ..obtain workspace for slv5pt
        allocate(aaf_5p(1:n2*(3*min(nx+2,ny+2)+1)))
        allocate(wrk0_5p(1:n2))
        allocate(wrk1_5p(1:1))
        allocate(wrk2_5p(1:1))
        allocate(iwrk_5p(1:n2))
      endif
      write(*,*) 'Allocated space for solver'
#ifdef _OPENMP
!$OMP END PARALLEL
#endif
      if (style.lt.2) then
        do ix = -1, nx
          do iy = -1, ny
            call xertst (
     .    (leftix(ix,iy).eq.-2 .or. leftix(ix,iy).eq.ix-1) .and.
     .     leftiy(ix,iy).eq.iy .and.
     .    (rightix(ix,iy).eq.nx+1 .or. rightix(ix,iy).eq.ix+1) .and.
     .     rightiy(ix,iy).eq.iy .and.
     .     bottomix(ix,iy).eq.ix .and.
     .    (bottomiy(ix,iy).eq.-2 .or. bottomiy(ix,iy).eq.iy-1) .and.
     .     topix(ix,iy).eq.ix .and.
     .    (topiy(ix,iy).eq.ny+1 .or. topiy(ix,iy).eq.iy+1),
     .   'Solver type inappropriate for non-trivial neighbourhoods')
          enddo
        enddo
      endif

      return
      end subroutine solver_init_5p

      end module b2mod_ma28
