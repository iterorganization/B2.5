      module b2us_map

      use b2mod_types
      implicit none
      
      private
      public :: alloc_mapping, dealloc_mapping,
     .          alloc_mapping_bc, dealloc_mapping_bc,
     .          alloc_mapping_rc, dealloc_mapping_rc
      
      type, public :: mapping
        integer :: nCv, nFc, nVx, nCg, nCi, nCmxVx, nCmxFc, nVmxCv,
     &             nVmxFc, nFs, mxnFs, nBc, mxnBc, nRc, mxnRc, nCmxNv

        integer, allocatable :: cvFcP(:,:) ! (nCv,2) pointing for every cell to the first index number in the cvFc list and its number of faces
        integer, allocatable :: cvFc(:)    ! listing for every cell the corresponding face numbers 
        integer, allocatable :: fcCv(:,:)  ! (nFc,2) listing for every face the two corresponding cells
        integer, allocatable :: fcVx(:,:)  ! (nFc,2) Hajo used: xVx(ix,iy) and yVx(ix,iy) are used but for each cell
        integer, allocatable :: cvVxP(:,:) ! (nCv,2) pointing for every cell to the first index number in the cvVx list and its number of vertices 
        integer, allocatable :: cvVx(:)    ! listing for every cell the corresponding vertex numbers        
        integer, allocatable :: vxFcP(:,:) ! (nVx,2) pointing for every vertex to the first index number in the vxFc list and its number of faces  
        integer, allocatable :: vxFc(:)    ! listing for every vertex the corresponding face numbers
        integer, allocatable :: vxCvP(:,:) ! (nVx,2) pointing for every cell to the first index number in the vxCv list and its number of cells   
        integer, allocatable :: vxCv(:)    ! listing for every vertex the corresponding cell numbers
        integer, allocatable :: cvFsP(:,:) ! (nFs,2) pointing for every 
        integer, allocatable :: cvFs(:)
        integer, allocatable :: fcFsP(:,:) ! (nFs,2) pointing for every 
        integer, allocatable :: fcFs(:)
        integer, allocatable :: fcReg(:) 
        integer, allocatable :: cvReg(:)
        logical, allocatable :: cvOnClosedSurface(:)

        integer, allocatable :: bcCvP(:,:) ! (nFs,2) pointing for every 
        integer, allocatable :: bcCv(:,:)
        real (kind=R8), allocatable :: bcFcOr (:)  ! multiplier to indicate whether outward normal coincindes with the face normal (1.0_R8) or not (-1.0_R8)

        integer, allocatable :: rcCvP(:,:) ! (nFs,2) pointing for every 
        integer, allocatable :: rcCv(:,:)
        real (kind=R8), allocatable :: rcFcOr (:)  ! multiplier to indicate whether outward normal coincindes with the face normal (1.0_R8) or not (-1.0_R8)

!        real (kind=R8), allocatable :: intfaceV(:,:) ! volume factors for interpolation from centers to faces; currently: still geo%fcVol
!        real (kind=R8), allocatable :: intvertexV(:) ! volume factors for interpolation from centers to vertices; currently: still based on cvVol
        real (kind=R8), allocatable :: intcellP(:)   ! factors for interpolation from faces to centers in poloidal direction
        real (kind=R8), allocatable :: intcellR(:)   ! factors for interpolation from faces to centers in radial direction

        integer, allocatable :: fcNv(:,:)  !(nFc,2): indicating for each cell of a face the number of the neighboring cell in its stencil
        integer, allocatable :: aaP(:,:) !(nCv,2): pointing for each cell to the start in cvNv array, and the number of points in its stencil
        integer, allocatable :: cvNv(:)    !(nFc,2): indicating for each cell the numbers of the cells in its stencil

      end type mapping

      contains

**********************************************************************

      subroutine alloc_mapping (m)

      implicit none
      type(mapping), intent(inout) :: m

      allocate (m%cvFcP(m%nCv,2)) 
      allocate (m%cvFc(m%nCmxFc))
      allocate (m%fcCv(m%nFc,2))
      allocate (m%fcVx(m%nFc,2))
      allocate (m%cvVxP(m%nCv,2))
      allocate (m%cvVx(m%nCmxVx))
      allocate (m%vxFcP(m%nVx,2)) 
      allocate (m%vxFc(m%nVmxFc))
      allocate (m%vxCvP(m%nVx,2))  
      allocate (m%vxCv(m%nVmxCv)) 
      allocate (m%cvFs(m%nCv))
      allocate (m%cvFsP(m%mxnFs,2))
      allocate (m%fcFs(m%nCv))
      allocate (m%fcFsP(m%mxnFs,2))
      allocate (m%fcReg(m%nFC))
      allocate (m%cvReg(m%nCv))
      allocate (m%cvOnClosedSurface(m%nCv))
      
!      allocate (m%intfaceV(m%nFc,2))
!      allocate (m%intvertexV(m%nCmxVx * m%nCv))
      allocate (m%intcellP(m%nCmxFc))
      allocate (m%intcellR(m%nCmxFc))

      allocate (m%fcNv(m%nFc,2))
      allocate (m%aaP(m%nCv,2))
      allocate (m%cvNv(m%nCmxNv))

!      allocate (m%cvNv(m%nCmxVx*m%nCmxVx*m%nCv)) ! very conservative; reduce!

      write (*,*) 'Mapping allocated'

      return
      end subroutine alloc_mapping 

**********************************************************************

      subroutine dealloc_mapping (m)

      implicit none
      type(mapping), intent(in out) :: m

      deallocate (m%cvFcP) 
      deallocate (m%cvFc)
      deallocate (m%fcCv)
      deallocate (m%fcVx)
      deallocate (m%cvVxP)
      deallocate (m%cvVx)
      deallocate (m%vxFcP) 
      deallocate (m%vxFc) 
      deallocate (m%vxCvP)  
      deallocate (m%vxCv) 
      deallocate (m%cvFs)
      deallocate (m%cvFsP)
      deallocate (m%fcFs)
      deallocate (m%fcFsP)
      deallocate (m%fcReg)
      deallocate (m%cvReg)
      deallocate (m%cvOnClosedSurface)

!      deallocate (m%intfaceV) 
!      deallocate (m%intvertexV) 
      deallocate (m%intcellP) 
      deallocate (m%intcellR)   

      deallocate (m%fcNv)
      deallocate (m%aaP)
      deallocate (m%cvNv)    

      return
      end subroutine dealloc_mapping 

**********************************************************************

      subroutine alloc_mapping_bc (m)

      implicit none
      type(mapping), intent(in out) :: m

      allocate (m%bcCvP  ( m % nBc,   2 ) )
      allocate (m%bcCv   ( m % mxnBc, 2 ) )
      allocate (m%bcFcOr ( m % mxnBc    ) )

      return
      end subroutine alloc_mapping_bc 

**********************************************************************

      subroutine dealloc_mapping_bc (m)

      implicit none
      type(mapping), intent(in out) :: m

      deallocate ( m % bcCvP  )
      deallocate ( m % bcCv   )
      deallocate ( m % bcFcOr )
      
      return
      end subroutine dealloc_mapping_bc

**********************************************************************

      subroutine alloc_mapping_rc (m)

      implicit none
      type(mapping), intent(in out) :: m

      allocate (m%rcCvP  ( m % nRc,   2 ) )
      allocate (m%rcCv   ( m % mxnRc, 2 ) )
      allocate (m%rcFcOr ( m % mxnRc    ) )

      return      
      end subroutine alloc_mapping_rc 

**********************************************************************

      subroutine dealloc_mapping_rc (m)

      implicit none
      type(mapping), intent(in out) :: m

      deallocate ( m % rcCvP  )
      deallocate ( m % rcCv   )
      deallocate ( m % rcFcOr )

      return      
      end subroutine dealloc_mapping_rc 

**********************************************************************

      end module b2us_map
