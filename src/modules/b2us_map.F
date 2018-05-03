      module b2us_map

      use b2mod_types
      implicit none
      
      private
      public :: alloc_mapping, dealloc_mapping
      
      type, public :: mapping
        integer :: nCv, nFc, nVx, nCg, nCi, nCmxVx, nCmxFc, nFs, mxnFs
        integer, allocatable :: cvFcP(:,:) ! (2,nCv) pointing for every cell to the first index number in the cvFc list and its number of faces
        integer, allocatable :: cvFc(:)    ! listing for every cell the corresponding face numbers 
	integer, allocatable :: fcCv(:,:)  ! (2,nFc) listing for every face the two corresponding cells
	integer, allocatable :: fcVx(:,:)  ! (2, nFc) Hajo used: xVx(ix,iy) and yVx(ix,iy) are used but for each cell
	integer, allocatable :: cvVxP(:,:) ! (2,nCv) pointing for every cell to the first index number in the cvVx list and its number of vertices 
	integer, allocatable :: cvVx(:)    ! listing for every cell the corresponding vertex numbers
	integer, allocatable :: vxFcP(:,:) ! (2,nVx) pointing for every vertex to the first index number in the vxFc list and its number of faces  
	integer, allocatable :: vxFc(:)    ! listing for every vertex the corresponding face numbers
	integer, allocatable :: vxCvP(:,:) ! (2,nVx) pointing for every cell to the first index number in the vxCv list and its number of cells   
	integer, allocatable :: vxCv(:)    ! listing for every vertex the corresponding cell numbers
        integer, allocatable :: cvFsP(:,:) ! (2,nFs) pointing for every 
        integer, allocatable :: cvFs(:)
        integer, allocatable :: fcFsP(:,:) ! (2,nFs) pointing for every 
        integer, allocatable :: fcFs(:)
        integer, allocatable :: fcReg(:) 
        integer, allocatable :: cvReg(:)
      end type mapping

      contains

**********************************************************************

      subroutine alloc_mapping (m)

      implicit none
      type(mapping), intent(in out) :: m

      allocate (m%cvFcP(m%nCv,2)) 
      allocate (m%cvFc(m%nCmxFc * m%nCv))
      allocate (m%fcCv(m%nFc,2))
      allocate (m%fcVx(m%nFc,2))
      allocate (m%cvVxP(m%nCv,2))
      allocate (m%cvVx(m%nCmxVx * m%nCv))
      allocate (m%vxFcP(m%nVx,2)) 
      allocate (m%vxFc(m%nCmxFc * m%nVx))
      allocate (m%vxCvP(m%nVx,2))  
      allocate (m%vxCv(m%nCmxVx * m%nCv)) 
      allocate (m%cvFs(m%nCv))
      allocate (m%cvFsP(m%mxnFs,2))
      allocate (m%fcFs(m%nCv))
      allocate (m%fcFsP(m%mxnFs,2))
      allocate (m%fcReg(m%nFC))
      allocate (m%cvReg(m%nCv))

      
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

      
      end subroutine dealloc_mapping 

**********************************************************************

      end module b2us_map
