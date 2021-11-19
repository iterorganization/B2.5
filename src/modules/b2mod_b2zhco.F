      module b2mod_b2zhco
      use b2mod_types
      implicit none
      integer :: zhd_corr=1, zhd_corrT=1, zhcscorr=1, zhcscorr_vq=1,
     &  zh_imp_analyt=0, zh_analyt_mdf=0, zh_dist_fun=0, zh_tf_toff=0
      real (kind=R8) :: zhflcorr=3.0_R8
      real (kind=R8) :: zhflcorrh=0.4_R8
      real (kind=R8), allocatable, save ::
     &  c_hw_save(:,:,:,:),
     &  c_R_TA(:,:,:,:), c_R_TA_nofl(:,:,:,:),
     &  c_R_wA(:,:,:,:),
     &  c_R_TB(:,:,:), c_R_TB_nofl(:,:,:),
     &  c_R_wB(:,:,:,:),
     &  c_R_w(:,:,:,:),
     &  tf_ton_mat(:,:)

*-----------------------------------------------------------------------
*.documentation
*
*  1. purpose
*
*     B2MOD_B2ZHCO stores the transport coefficients and the flags
*     used by the Zhdanov closure scheme.
*
*  3. description 
*
*     This module contains:
*     - c_R_TA: the thermal force transport coefficient
*               for an isonuclear sequence
*     - c_R_TA_nofl: the thermal force transport coefficient
*               for an isonuclear sequence without flux limiting
*     - c_R_TB: the thermal force transport coefficient correction
*               for each charge state
*     - c_R_TB_nofl: the thermal force transport coefficient correction
*               for each charge state without flux limiting
*     - c_R_wA: the friction force transport coefficient
*               for an isonuclear sequence
*     - c_R_wB: the friction force transport coefficient correction
*               for each charge state
*     - c_R_w: the total friction force transport coefficient 
*     - c_hwA: the velocity-dependent ion heat flux transport coefficient
*              for an isonuclear sequence
*
*-----------------------------------------------------------------------
      contains

      subroutine alloc_b2mod_zhco(nnucl, nx, ny)
      implicit none
      integer :: nnucl, nx, ny
      allocate(c_hw_save(-1:nx,-1:ny,0:nnucl-1,0:nnucl-1))
      allocate(c_R_TA(-1:nx,-1:ny,0:nnucl-1,0:nnucl-1))
      allocate(c_R_TA_nofl(-1:nx,-1:ny,0:nnucl-1,0:nnucl-1))
      allocate(c_R_wA(-1:nx,-1:ny,0:nnucl-1,0:nnucl-1))
      allocate(c_R_TB(-1:nx,-1:ny,0:nnucl-1))
      allocate(c_R_TB_nofl(-1:nx,-1:ny,0:nnucl-1))
      allocate(c_R_wB(-1:nx,-1:ny,0:nnucl-1,0:nnucl-1))
      allocate(c_R_w(-1:nx,-1:ny,0:nnucl-1,0:nnucl-1))
      allocate(tf_ton_mat(0:nnucl-1,0:nnucl-1))
      tf_ton_mat=1.0_R8

      end subroutine alloc_b2mod_zhco

      subroutine dealloc_b2mod_zhco
      implicit none
      if(.not.allocated(c_R_TA)) return
      deallocate(c_hw_save)
      deallocate(c_R_TA)
      deallocate(c_R_TA_nofl)
      deallocate(c_R_wA)
      deallocate(c_R_TB)
      deallocate(c_R_TB_nofl)
      deallocate(c_R_wB)
      deallocate(c_R_w)
      deallocate(tf_ton_mat)

      end subroutine dealloc_b2mod_zhco

      end module b2mod_b2zhco
