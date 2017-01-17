      module b2mod_ma28_for_11diag                                     !sv 13.10.06 {

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

      contains

      subroutine alloc_b2mod_ma28_for_11diag

      write(*,*) 'Allocating space for b2mod_ma28_for_11diag'
      write(*,*) 'nnd        = ',nnd
      write(*,*) 'n_nonzeros = ',n_nonzeros
      write(*,*) 'solvdim    = ',solvdim
      write(*,*) 'solvdim1   = ',solvdim1

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
      allocate(irnh(1:solvdim))
      allocate(icnh(1:solvdim))
      allocate(icnhkeep(1:solvdim))
      allocate(ikeep(1:nnd,5))
      allocate(iw28(1:nnd,8))
      allocate(irnh_save(1:n_nonzeros))
      allocate(icnh_save(1:n_nonzeros))

      write(*,*) 'Allocated space for b2mod_ma28_for_11diag'

      end subroutine alloc_b2mod_ma28_for_11diag

      subroutine dealloc_b2mod_ma28_for_11diag

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

      end subroutine dealloc_b2mod_ma28_for_11diag

      end module b2mod_ma28_for_11diag                                 !sv 13.10.06 }
