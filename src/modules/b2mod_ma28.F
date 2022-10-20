      module b2mod_ma28
      use b2mod_types

      real (kind=R8), allocatable, save ::
     * r1(:),
     * s1(:),
     * a1h(:),
     * w(:),
     * res(:),
     * a1h_save(:)

      integer, allocatable, save ::
     * irn(:),
     * icn(:),
     * ish(:),
     * jsh(:),
     * ijsh(:),
     * irnh(:),
     * icnh(:),
     * icnhkeep(:),
     * ikeep(:,:),
     * iw28(:,:),
     * irnh_save(:),
     * icnh_save(:)

      integer nc0,nltrsol,nnd,n_nonzeros,solvdim,solvdim1
      logical, save :: lfirst_ilu = .true.

!!FR All the allocatable arrays need to be declared as THREADPRIVATE
!!FR as each thread needs its own copy of the arrays used in the solver.
!!FR All the integers also need to be threadprivate as they are used
!!FR in the ncall==0 block in b2ux5p that defines the array sizes.
      !$OMP THREADPRIVATE (r1, s1, a1h, w, res, a1h_save)
      !$OMP THREADPRIVATE (irn, icn, ish, jsh, ijsh, irnh, icnh)
      !$OMP THREADPRIVATE (icnhkeep, ikeep, iw28, irnh_save, icnh_save)
      !$OMP THREADPRIVATE (nc0, nltrsol, nnd, n_nonzeros)
      !$OMP THREADPRIVATE (solvdim, solvdim1)

      contains

      subroutine alloc_b2mod_ma28

      write(*,*) 'Allocating space for b2mod_ma28'
      write(*,*) 'nnd = ',nnd
      write(*,*) 'n_nonzeros = ',n_nonzeros
      write(*,*) 'solvdim = ',solvdim
      write(*,*) 'solvdim1 = ',solvdim1

      allocate(r1(1:nnd))
      allocate(s1(1:nnd))
      allocate(a1h(1:solvdim))
      allocate(w(1:nnd))
      allocate(res(1:nnd))
      allocate(a1h_save(1:n_nonzeros))
      allocate(irn(1:n_nonzeros))
      allocate(icn(1:n_nonzeros))
      allocate(ish(1:n_nonzeros))
      allocate(jsh(1:n_nonzeros))
      allocate(ijsh(1:n_nonzeros))
      allocate(irnh(1:solvdim1))
      allocate(icnh(1:solvdim))
      allocate(icnhkeep(1:solvdim))
      allocate(ikeep(1:nnd,5))
      allocate(iw28(1:nnd,8))
      allocate(irnh_save(1:n_nonzeros))
      allocate(icnh_save(1:n_nonzeros))

      write(*,*) 'Allocated space for b2mod_ma28'

      end subroutine alloc_b2mod_ma28

      subroutine dealloc_b2mod_ma28

      if (.not.allocated(r1)) return
      deallocate(r1)
      deallocate(s1)
      deallocate(a1h)
      deallocate(w)
      deallocate(res)
      deallocate(a1h_save)
      deallocate(irn)
      deallocate(icn)
      deallocate(ish)
      deallocate(jsh)
      deallocate(ijsh)
      deallocate(irnh)
      deallocate(icnh)
      deallocate(icnhkeep)
      deallocate(ikeep)
      deallocate(iw28)
      deallocate(irnh_save)
      deallocate(icnh_save)

      end subroutine dealloc_b2mod_ma28

      subroutine restart_ma28
      implicit none
      lfirst_ilu=.true.
      write(*,*) 'Restarting ma28'
      return
      end subroutine restart_ma28

      end module b2mod_ma28
