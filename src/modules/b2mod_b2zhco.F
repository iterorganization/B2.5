      module b2mod_b2zhco
      use b2mod_types
      implicit none
      character*1, allocatable, save :: zh_tf_toff(:)
      real (kind=R8), allocatable, save ::
     &  c_hw_save(:,:,:),
     &  c_R_TA(:,:,:), c_R_TA_nofl(:,:,:),
     &  c_R_wA(:,:,:),
     &  c_R_TB(:,:), c_R_TB_nofl(:,:),
     &  c_R_wB(:,:,:),
     &  c_R_w(:,:,:),
     &  tf_ton_mat(:,:),
     &  c_hTA_an_fl_save(:,:),
     &  c_R_TA_an_save(:,:,:)
      real (kind=R8), allocatable, save :: amfact(:)

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

      subroutine alloc_b2mod_zhco(nnucl, nCv)
      implicit none
      integer :: nnucl, nCv
      allocate(c_hw_save(nCv,0:nnucl-1,0:nnucl-1))
      allocate(c_R_TA(nCv,0:nnucl-1,0:nnucl-1))
      allocate(c_R_TA_nofl(nCv,0:nnucl-1,0:nnucl-1))
      allocate(c_R_wA(nCv,0:nnucl-1,0:nnucl-1))
      allocate(c_R_TB(nCv,0:nnucl-1))
      allocate(c_R_TB_nofl(nCv,0:nnucl-1))
      allocate(c_R_wB(nCv,0:nnucl-1,0:nnucl-1))
      allocate(c_R_w(nCv,0:nnucl-1,0:nnucl-1))
      allocate(tf_ton_mat(0:nnucl-1,0:nnucl-1))
      allocate(c_hTA_an_fl_save(nCv,0:nnucl-1))
      allocate(c_R_TA_an_save(nCv,0:nnucl-1,0:nnucl-1))
      allocate(zh_tf_toff(max(1,(nnucl*(nnucl+1))/2-nnucl)))
      tf_ton_mat=1.0_R8

      return
      end subroutine alloc_b2mod_zhco

      subroutine alloc_b2mod_zhco_amfact(nnucl)
      implicit none
      integer :: nnucl
      allocate(amfact(0:nnucl-1))

      end subroutine alloc_b2mod_zhco_amfact

      subroutine dealloc_b2mod_zhco
      implicit none
      if (allocated(c_R_TA)) then
        deallocate(c_hw_save)
        deallocate(c_R_TA)
        deallocate(c_R_TA_nofl)
        deallocate(c_R_wA)
        deallocate(c_R_TB)
        deallocate(c_R_TB_nofl)
        deallocate(c_R_wB)
        deallocate(c_R_w)
        deallocate(tf_ton_mat)
        deallocate(zh_tf_toff)
        deallocate(c_hTA_an_fl_save)
        deallocate(c_R_TA_an_save)
      endif

      if(allocated(amfact)) then
        deallocate(amfact)
      end if

      return
      end subroutine dealloc_b2mod_zhco

      end module b2mod_b2zhco
