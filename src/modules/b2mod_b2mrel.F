      module b2mod_b2mrel
      use b2mod_types
      implicit none
      integer, allocatable, save ::
     &  nucl2s(:,:), znnucl_int(:), s2nucl(:)
      real (kind=R8), allocatable, save ::
     &  amnucl(:), znnucl(:),
     &  N2(:,:), N3(:,:),
     &  N4(:,:), N5(:,:),
     &  N6(:,:), N8(:,:),
     &  N9(:,:),
     &  N10(:,:), N11(:,:),
     &  N12(:,:),
     &  N13(:,:), N14(:,:),
     &  N15(:,:), N16(:,:),
     &  mu(:,:), mum(:,:),
     &  rmu(:,:), rmd2(:),
     &  mdm(:,:),
     &  amnucl_mdf(:),
     &  N2_mdf(:,:), N3_mdf(:,:),
     &  N4_mdf(:,:), N5_mdf(:,:),
     &  N6_mdf(:,:), N8_mdf(:,:),
     &  N9_mdf(:,:),
     &  N10_mdf(:,:), N11_mdf(:,:),
     &  N12_mdf(:,:),
     &  N13_mdf(:,:), N14_mdf(:,:),
     &  N15_mdf(:,:), N16_mdf(:,:),
     &  mu_mdf(:,:), mum_mdf(:,:),
     &  rmu_mdf(:,:), rmd2_mdf(:),
     &  mdm_mdf(:,:)
*-----------------------------------------------------------------------
*.documentation
*
*  1. purpose
*
*     B2MOD_B2MREL stores the N-coefficients and the mass-related variables.
*
*  3. description
*
*     This module contains:
*     - N*: the N-coefficients, which are used in the G-objects
*     - amnucl: mass of the nucleus of the isonuclear sequence
*     - znnucl: charge of the nucleus of the isonuclear sequence
*     - znnucl_int: the same as in integer array
*     - nucl2s: the matrix that connects the isonuclear sequence to the species sequence
*     - s2nucl: the array that connects the species sequence to the isonuclear sequence
*     - mu: mu_{ab} the reduced mass (in amu)
*     - rmu: sqrt{mu_{ab}} the square root of the reduced mass (in amu)
*     - mdm: m_a/m_b the mass ratio
*     - mum: mu_{ab}/m_a the reduced mass ratio
*
*-----------------------------------------------------------------------
      contains

      subroutine alloc_b2mod_mrel(nnucl,znmax,ns)
      implicit none
      integer :: nnucl, znmax, ns

      allocate(amnucl(0:nnucl-1))
      allocate(znnucl(0:nnucl-1))
      allocate(N2(0:nnucl-1,0:nnucl-1))
      allocate(N3(0:nnucl-1,0:nnucl-1))
      allocate(N4(0:nnucl-1,0:nnucl-1))
      allocate(N5(0:nnucl-1,0:nnucl-1))
      allocate(N6(0:nnucl-1,0:nnucl-1))
      allocate(N8(0:nnucl-1,0:nnucl-1))
      allocate(N9(0:nnucl-1,0:nnucl-1))
      allocate(N10(0:nnucl-1,0:nnucl-1))
      allocate(N11(0:nnucl-1,0:nnucl-1))
      allocate(N12(0:nnucl-1,0:nnucl-1))
      allocate(N13(0:nnucl-1,0:nnucl-1))
      allocate(N14(0:nnucl-1,0:nnucl-1))
      allocate(N15(0:nnucl-1,0:nnucl-1))
      allocate(N16(0:nnucl-1,0:nnucl-1))
      allocate(mu(0:nnucl-1,0:nnucl-1))
      allocate(mum(0:nnucl-1,0:nnucl-1))
      allocate(rmu(0:nnucl-1,0:nnucl-1))
      allocate(rmd2(0:nnucl-1))
      allocate(mdm(0:nnucl-1,0:nnucl-1))
      allocate(amnucl_mdf(0:nnucl-1))
      allocate(N2_mdf(0:nnucl-1,0:nnucl-1))
      allocate(N3_mdf(0:nnucl-1,0:nnucl-1))
      allocate(N4_mdf(0:nnucl-1,0:nnucl-1))
      allocate(N5_mdf(0:nnucl-1,0:nnucl-1))
      allocate(N6_mdf(0:nnucl-1,0:nnucl-1))
      allocate(N8_mdf(0:nnucl-1,0:nnucl-1))
      allocate(N9_mdf(0:nnucl-1,0:nnucl-1))
      allocate(N10_mdf(0:nnucl-1,0:nnucl-1))
      allocate(N11_mdf(0:nnucl-1,0:nnucl-1))
      allocate(N12_mdf(0:nnucl-1,0:nnucl-1))
      allocate(N13_mdf(0:nnucl-1,0:nnucl-1))
      allocate(N14_mdf(0:nnucl-1,0:nnucl-1))
      allocate(N15_mdf(0:nnucl-1,0:nnucl-1))
      allocate(N16_mdf(0:nnucl-1,0:nnucl-1))
      allocate(mu_mdf(0:nnucl-1,0:nnucl-1))
      allocate(mum_mdf(0:nnucl-1,0:nnucl-1))
      allocate(rmu_mdf(0:nnucl-1,0:nnucl-1))
      allocate(rmd2_mdf(0:nnucl-1))
      allocate(mdm_mdf(0:nnucl-1,0:nnucl-1))
      allocate(nucl2s(0:nnucl-1,0:znmax))
      allocate(znnucl_int(0:nnucl-1))
      allocate(s2nucl(0:ns-1))
      end subroutine alloc_b2mod_mrel

      subroutine dealloc_b2mod_mrel
      if (.not.allocated(amnucl)) return
      deallocate(amnucl)
      deallocate(znnucl)
      deallocate(N2)
      deallocate(N3)
      deallocate(N4)
      deallocate(N5)
      deallocate(N6)
      deallocate(N8)
      deallocate(N9)
      deallocate(N10)
      deallocate(N11)
      deallocate(N12)
      deallocate(N13)
      deallocate(N14)
      deallocate(N15)
      deallocate(N16)
      deallocate(mu)
      deallocate(mum)
      deallocate(rmu)
      deallocate(rmd2)
      deallocate(mdm)
      deallocate(amnucl_mdf)
      deallocate(N2_mdf)
      deallocate(N3_mdf)
      deallocate(N4_mdf)
      deallocate(N5_mdf)
      deallocate(N6_mdf)
      deallocate(N8_mdf)
      deallocate(N9_mdf)
      deallocate(N10_mdf)
      deallocate(N11_mdf)
      deallocate(N12_mdf)
      deallocate(N13_mdf)
      deallocate(N14_mdf)
      deallocate(N15_mdf)
      deallocate(N16_mdf)
      deallocate(mu_mdf)
      deallocate(mum_mdf)
      deallocate(rmu_mdf)
      deallocate(rmd2_mdf)
      deallocate(mdm_mdf)
      deallocate(nucl2s)
      deallocate(znnucl_int)
      deallocate(s2nucl)
      end subroutine dealloc_b2mod_mrel

      end module b2mod_b2mrel
