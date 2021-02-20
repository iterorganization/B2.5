      module b2mod_b2zhco
      use b2mod_types
      real (kind=R8), allocatable ::
     & c_R_TA(:,:,:,:), c_R_wA(:,:,:,:)
      contains

      subroutine init_zhco(nnucl, nx, ny)
      integer :: nnucl, nx, ny
      
      allocate(c_R_TA(-1:nx,-1:ny,0:nnucl-1,0:nnucl-1))
      allocate(c_R_wA(-1:nx,-1:ny,0:nnucl-1,0:nnucl-1))

      end subroutine init_zhco

      end module b2mod_b2zhco
