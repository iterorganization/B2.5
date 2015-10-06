      module b2mod_geo_corner
      use b2mod_types
      integer, save :: nmdpl                                         !srv 22.05.09
      real (kind=R8), allocatable, save ::
     * OnedBsq(:,:),
     * OnedBsqc(:,:),
     * wbbl(:,:,:),
     * wbbr(:,:,:),
     * wbbv(:,:,:),
     * wbbc(:,:,:),
     * lengthy(:,:),                                                 !srv 24.04.02
     * lengthy1(:,:),                                                !srv 21.05.08
     * bzb(:,:),                                                     !srv 13.10.06
     * len_surf(:,:),                                                !srv 02.05.09
     * crdy_surf(:,:),
     * crdy_c_cell(:,:)                                              !srv 02.05.09

      contains

      subroutine alloc_b2mod_geo_corner(nxd,nyd)
      integer nxd, nyd                                               !srv 04.06.99 
      allocate (OnedBsq(-1:nxd,-1:nyd))
      allocate (OnedBsqc(-1:nxd+1,-1:nyd+1))
      allocate (wbbl(-1:nxd,-1:nyd,0:3))
      allocate (wbbr(-1:nxd,-1:nyd,0:3))
      allocate (wbbv(-1:nxd,-1:nyd,0:3))
      allocate (wbbc(-1:nxd+1,-1:nyd+1,0:3))
      allocate (lengthy(-1:nxd,-1:nyd))                              !srv 24.04.02
      allocate (lengthy1(-1:nxd,-1:nyd))                             !srv 21.05.08
      allocate (len_surf(-1:nxd,-1:nyd))                             !srv 02.05.09
      allocate (crdy_surf(-1:nxd,-1:nyd))
      allocate (crdy_c_cell(-1:nxd,-1:nyd))                          !srv 02.05.09
      allocate (bzb(-1:nxd,-1:nyd))                                  !srv 13.10.06

      write (*,*) 'Allocated space for b2mod_geo_corner' 

      end subroutine alloc_b2mod_geo_corner
	 
      subroutine dealloc_b2mod_geo_corner

      deallocate(OnedBsq)
      deallocate(OnedBsqc)
      deallocate(wbbl)
      deallocate(wbbr)
      deallocate(wbbv)
      deallocate(wbbc)
      deallocate(lengthy)                                            !srv 24.04.02
      deallocate(lengthy1)                                           !srv 21.05.08
      deallocate(len_surf)                                           !srv 02.05.09
      deallocate(crdy_surf)
      deallocate(crdy_c_cell)                                        !srv 02.05.09
      deallocate(bzb)                                                !srv 13.10.06

      end subroutine dealloc_b2mod_geo_corner

      end module b2mod_geo_corner
