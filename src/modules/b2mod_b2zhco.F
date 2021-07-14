      module b2mod_b2zhco
      use b2mod_types
      implicit none
      integer :: zhd_corr=1,zhcscorr=1,zh_imp_analyt=0,
     & zh_analyt_mdf=0
      real (kind=R8) :: zhflcorr=1.0_R8
      real (kind=R8) :: zhflcorrh=1.0_R8
      real (kind=R8), allocatable, save ::
     & c_hwA(:,:,:,:),
     & c_R_TA(:,:,:,:),c_R_TA_nofl(:,:,:,:),
     & c_R_wA(:,:,:,:),
     & c_R_TB(:,:,:),c_R_TB_nofl(:,:,:),
     & c_R_wB(:,:,:,:),
     & c_R_w(:,:,:,:)

*-----------------------------------------------------------------------
*.documentation
*
*  1. purpose
*
*     b2mod_b2zhco store the transport coefficients and the flags 
*
*  3. description 
*
*     This routine computes:
*     - c_R_TA: the thermal force transport coefficient for isonuclear ions
*     - c_R_TA_nofl: the thermal force transport coefficient for isonuclear ions without flux limiting
*     - c_R_TB: the thermal force transport coefficient correction for each charge state
*     - c_R_TB_nofl: the thermal force transport coefficient correction for each charge state without flux limiting
*     - c_R_wA: the friction force transport coefficient for isonuclear ions
*     - c_R_wB: the friction force transport coefficient correction for each charge state
*     - c_R_w: the total friction force transport coefficient 
*     - c_hwA: the velocity dependent ion heat flux tansport coefficient for isonuclear ions
*
*-----------------------------------------------------------------------
      contains

      subroutine init_zhco(nnucl, nx, ny)
      implicit none
      integer :: nnucl, nx, ny
      allocate(c_hwA(-1:nx,-1:ny,0:nnucl-1,0:nnucl-1))
      allocate(c_R_TA(-1:nx,-1:ny,0:nnucl-1,0:nnucl-1))
      allocate(c_R_TA_nofl(-1:nx,-1:ny,0:nnucl-1,0:nnucl-1))
      allocate(c_R_wA(-1:nx,-1:ny,0:nnucl-1,0:nnucl-1))
      allocate(c_R_TB(-1:nx,-1:ny,0:nnucl-1))
      allocate(c_R_TB_nofl(-1:nx,-1:ny,0:nnucl-1))
      allocate(c_R_wB(-1:nx,-1:ny,0:nnucl-1,0:nnucl-1))
      allocate(c_R_w(-1:nx,-1:ny,0:nnucl-1,0:nnucl-1))

      end subroutine init_zhco

      end module b2mod_b2zhco
