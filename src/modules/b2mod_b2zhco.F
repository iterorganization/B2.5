      module b2mod_b2zhco
      use b2mod_types
      integer :: zhd_corr=1
      real (kind=R8) :: zhflcorr=1.0_R8
      real (kind=R8), allocatable, save ::
     & c_R_TA(:,:,:,:),c_R_TA_nofl(:,:,:,:),
     & c_R_wA(:,:,:,:),
     & c_R_w(:,:,:,:)

      contains

      subroutine init_zhco(nnucl, nx, ny)
      integer :: nnucl, nx, ny
      
      allocate(c_R_TA(-1:nx,-1:ny,0:nnucl-1,0:nnucl-1))
      allocate(c_R_TA_nofl(-1:nx,-1:ny,0:nnucl-1,0:nnucl-1))
      allocate(c_R_wA(-1:nx,-1:ny,0:nnucl-1,0:nnucl-1))
      allocate(c_R_w(-1:nx,-1:ny,0:nnucl-1,0:nnucl-1))

      end subroutine init_zhco

      end module b2mod_b2zhco
