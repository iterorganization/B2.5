      module b2mod_b2mrel
      use b2mod_types
      integer, allocatable :: 
     & nucl2s(:,:), znnucl_int(:), s2nucl(:) 
      real (kind=R8), allocatable ::
     & amnucl(:), znnucl(:), 
     & N2(:,:), N5(:,:), 
     & N6(:,:), N8(:,:),
     & N9(:,:),
     & N10(:,:), N11(:,:),
     & N12(:,:),
     & mu(:,:), mum(:,:),
     & rmu(:,:),
     & mdm(:,:)
      contains

      subroutine init_mrel(nnucl,znmax,ns)
      integer :: nnucl, znmax
      
      allocate(amnucl(0:nnucl-1))
      allocate(znnucl(0:nnucl-1))
      allocate(N2(0:nnucl-1,0:nnucl-1))
      allocate(N5(0:nnucl-1,0:nnucl-1)) 
      allocate(N6(0:nnucl-1,0:nnucl-1))
      allocate(N8(0:nnucl-1,0:nnucl-1))
      allocate(N9(0:nnucl-1,0:nnucl-1))
      allocate(N10(0:nnucl-1,0:nnucl-1))
      allocate(N11(0:nnucl-1,0:nnucl-1))
      allocate(N12(0:nnucl-1,0:nnucl-1))
      allocate(mu(0:nnucl-1,0:nnucl-1)) 
      allocate(mum(0:nnucl-1,0:nnucl-1))
      allocate(rmu(0:nnucl-1,0:nnucl-1))
      allocate(mdm(0:nnucl-1,0:nnucl-1))
      allocate(nucl2s(0:nnucl-1,0:znmax))
      allocate(znnucl_int(0:nnucl-1))
      allocate(s2nucl(0:ns-1))
      end subroutine init_mrel

      end module b2mod_b2mrel
