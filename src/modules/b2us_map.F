      module b2us_map

      use b2mod_types
      implicit none
      
      private
      public :: alloc_mapping, dealloc_mapping,
     .          alloc_mapping_bc, dealloc_mapping_bc,
     .          alloc_mapping_rc, dealloc_mapping_rc
      
      type, public :: mapping
        integer :: nCv, nFc, nVx, nCg, nCi, nCmxVx, nCmxFc, nVmxCv,
     &             nVmxFc, nFs, nBc, mxnBc, nRc, mxnRc, nCmxNv

        integer, allocatable :: cvFcP(:,:)             ! (nCv,2) pointing for every cell to the first index number in the cvFc list and its number of faces
        integer, allocatable :: cvFc(:)                !         listing for every cell the corresponding face numbers 
        integer, allocatable :: fcCv(:,:)              ! (nFc,2) listing for every face the two corresponding cells
        integer, allocatable :: fcVx(:,:)              ! (nFc,2) listing for every face the two corresponding vertices
        integer, allocatable :: cvVxP(:,:)             ! (nCv,2) pointing for every cell to the first index number in the cvVx list and its number of vertices 
        integer, allocatable :: cvVx(:)                !         listing for every cell the corresponding vertex numbers        
        integer, allocatable :: vxFcP(:,:)             ! (nVx,2) pointing for every vertex to the first index number in the vxFc list and its number of faces  
        integer, allocatable :: vxFc(:)                !         listing for every vertex the corresponding face numbers
        integer, allocatable :: vxCvP(:,:)             ! (nVx,2) pointing for every cell to the first index number in the vxCv list and its number of cells   
        integer, allocatable :: vxCv(:)                !         listing for every vertex the corresponding cell numbers
        integer, allocatable :: fsCvP(:,:)             ! (nFs,2) pointing for every flux surface to the first index number in the fsCv list and its number of cells
        integer, allocatable :: fsCv(:)                !         listing for every flux surface the corresponding cell numbers
        integer, allocatable :: fsFcP(:,:)             ! (nFs,2) pointing for every flux surface to the first index number in the fsFc list and its number of cell faces
        integer, allocatable :: fsFc(:)                !         listing for every flux surface the corresponding cell face numbers
        integer, allocatable :: fcReg(:)               ! Region numbers of the cell faces 
        integer, allocatable :: cvReg(:)               ! Region numbers of the cells
        logical, allocatable :: cvOnClosedSurface(:)   ! Logical to indicate whether a cell is on a closed flux surface

        integer, allocatable :: bcCvP(:,:)             ! (nBc,2) pointing for each boundary condition to the first index in the bcCv list, and its number of guard cells
        integer, allocatable :: bcCv(:,:)              !         listing for each boundary condition the guard cells, and the corresponding domain cells
        integer, allocatable :: bcFcP(:,:)             ! (nBc,2) pointing for each boundary condition to the first index in the bcFc list, and its number of boundary faces
        integer, allocatable :: bcFc(:)                !         listing for each boundary condition the boundary faces
        real (kind=R8), allocatable :: bcFcOr (:)      ! Multiplier to indicate whether outward normal coincindes with the face normal (1.0_R8) or not (-1.0_R8)

        integer, allocatable :: rcCvP(:,:)             ! (nRc,2) pointing for each stratum to the first index in the rcCv list, and its number of guard cells
        integer, allocatable :: rcCv(:,:)              !         listing for each stratum the guard cells, and the corresponding domain cells
        real (kind=R8), allocatable :: rcFcOr (:)      ! Multiplier to indicate whether outward normal coincindes with the face normal (1.0_R8) or not (-1.0_R8)

!        real (kind=R8), allocatable :: intfaceV(:,:)   ! volume factors for interpolation from centers to faces; currently: geo%fcVol
!        real (kind=R8), allocatable :: intvertexV(:)   ! volume factors for interpolation from centers to vertices; currently: geo%vxVol
        real (kind=R8), allocatable :: intcellP(:)     ! factors for interpolation from faces to centers in poloidal direction
        real (kind=R8), allocatable :: intcellR(:)     ! factors for interpolation from faces to centers in radial direction

        integer, allocatable :: cvNvP(:,:)             !(nCv,2): pointing for each cell to the start in cvNv array, and the number of points in its stencil
        integer, allocatable :: cvNv(:)                !         listing for each cell the numbers of the cells in its stencil

        integer :: isClassicalGrid                     ! indicates whether grid is derived from an "original", structured grid
        integer :: nx, ny, nncut                       ! some parameters related to the original grid

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
      allocate (m%fsCv(m%nCv))
      allocate (m%fsCvP(m%nFs,2))
      allocate (m%fsFc(m%nFc))
      allocate (m%fsFcP(m%nFs,2))
      allocate (m%fcReg(m%nFC))
      allocate (m%cvReg(m%nCv))
      allocate (m%cvOnClosedSurface(m%nCv))
      
!      allocate (m%intfaceV(m%nFc,2))
!      allocate (m%intvertexV(m%nCmxVx * m%nCv))
      allocate (m%intcellP(m%nCmxFc))
      allocate (m%intcellR(m%nCmxFc))

      allocate (m%cvNvP(m%nCv,2))
      allocate (m%cvNv(m%nCmxNv))

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
      deallocate (m%fsCv)
      deallocate (m%fsCvP)
      deallocate (m%fsFc)
      deallocate (m%fsFcP)
      deallocate (m%fcReg)
      deallocate (m%cvReg)
      deallocate (m%cvOnClosedSurface)

!      deallocate (m%intfaceV) 
!      deallocate (m%intvertexV) 
      deallocate (m%intcellP) 
      deallocate (m%intcellR)   

      deallocate (m%cvNvP)
      deallocate (m%cvNv)    

      return
      end subroutine dealloc_mapping 

**********************************************************************

      subroutine alloc_mapping_bc (m)

      implicit none
      type(mapping), intent(in out) :: m

      allocate (m%bcCvP  ( m % nBc,   2 ) )
      allocate (m%bcCv   ( m % mxnBc, 2 ) )
      allocate (m%bcFcP  ( m % mxnBc, 2 ) )
      allocate (m%bcFc   ( m % mxnBc    ) )
      allocate (m%bcFcOr ( m % mxnBc    ) )

      return
      end subroutine alloc_mapping_bc 

**********************************************************************

      subroutine dealloc_mapping_bc (m)

      implicit none
      type(mapping), intent(in out) :: m

      deallocate ( m % bcCvP  )
      deallocate ( m % bcCv   )
      deallocate ( m % bcFcP  )
      deallocate ( m % bcFc   )
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
