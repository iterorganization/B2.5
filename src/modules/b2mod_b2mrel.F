      module b2mod_b2mrel
      use b2mod_types
      implicit none
      integer, allocatable :: 
     & nucl2s(:,:), znnucl_int(:), s2nucl(:) 
      real (kind=R8), allocatable ::
     & amnucl(:), znnucl(:), 
     & N2(:,:), N3(:,:), 
     & N4(:,:), N5(:,:), 
     & N6(:,:), N8(:,:),
     & N9(:,:),
     & N10(:,:), N11(:,:),
     & N12(:,:),
     & N13(:,:), N14(:,:),
     & N15(:,:), N16(:,:),
     & mu(:,:), mum(:,:),
     & rmu(:,:),rmd2(:),
     & mdm(:,:)
*-----------------------------------------------------------------------
*.documentation
*
*  1. purpose
*
*     b2mod_b2mrel store the N-coefficients and the mass-related variables 
*
*  3. description 
*
*     This routine computes:
*     - N*: the N-coefficients, which are used in the G-objects
*     - amnucl: mass of the nucleus of the isonuclear species
*     - znnucl: charge of the nucleus of the isonuclear species
*     - znnucl_int: the same in integer values 
*     - nucl2s: the matrix that connect the isonuclear sequence to the species sequence 
*     - s2nucl: the array that connect the species sequence to the isonuclear sequence 
*     - mu: mu_{ab} the reduced mass (in atom mass numbers)
*     - rmu: sqrt{mu_{ab}} the square root of the reduced mass (in atom mass numbers)
*     - mdm: m_a/m_b the mass to mass ratio
*     - mum: mu_{ab}/m_a the reduced mass to mass ratio
*
*-----------------------------------------------------------------------
      contains

      subroutine init_mrel(nnucl,znmax,ns)
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
      allocate(nucl2s(0:nnucl-1,0:znmax))
      allocate(znnucl_int(0:nnucl-1))
      allocate(s2nucl(0:ns-1))
      end subroutine init_mrel

      end module b2mod_b2mrel
