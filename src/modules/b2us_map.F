      module b2us_map

      use b2mod_types
      implicit none
      
      private
      public :: alloc_mapping, dealloc_mapping,
     .          alloc_mapping_bc, dealloc_mapping_bc,
     .          alloc_mapping_rc, dealloc_mapping_rc,
     .          read_mapping, write_mapping,
     .          init_mapping, get_version
      
      type, public :: mapping
        integer :: nCv, nFc, nVx, nCg, nCi, nCmxVx, nCmxFc, nVmxCv,
     &             nVmxFc, nFs, nBc, mxnBc, nRc, mxnRc, nCmxNv
        integer :: nnreg(0:1)

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
        integer, allocatable :: rcFcP(:,:)             ! (nRc,2) pointing for each stratum to the first index in the rcFc list, and its number of boundary faces
        integer, allocatable :: rcFc(:)                !         listing for each stratum the boundary faces
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

      if (allocated(m%cvFcP)) return

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

      write (*,*) '=================================='
      write (*,*) 'Mapping allocated'
      write (*,*) 'nCv, nCi, nCg', m%nCv, m%nCi, m%nCg
      write (*,*) 'nFc, nVx, nFs', m%nFc, m%nVx, m%nFs
      write (*,*) 'nCmxVx, nCmxFc', m%nCmxVx, m%nCmxFc
      write (*,*) 'nVmxCv, nVmxFc', m%nVmxCv, m%nVmxFc
      write (*,*) '=================================='

      return
      end subroutine alloc_mapping 

**********************************************************************

      subroutine dealloc_mapping (m)

      implicit none
      type(mapping), intent(in out) :: m

      call dealloc_mapping_bc(m)
      call dealloc_mapping_rc(m)

      if (.not.allocated(m%cvFcP)) return

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

      if (allocated(m%bcCvP)) return

      allocate (m%bcCvP  ( m % nBc,   2 ) )
      allocate (m%bcCv   ( m % mxnBc, 2 ) )
      allocate (m%bcFcP  ( m % nBc,   2 ) )
      allocate (m%bcFc   ( m % mxnBc    ) )
      allocate (m%bcFcOr ( m % mxnBc    ) )

      return
      end subroutine alloc_mapping_bc 

**********************************************************************

      subroutine dealloc_mapping_bc (m)

      implicit none
      type(mapping), intent(in out) :: m

      if (.not.allocated(m%bcCvP)) return

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

      if (allocated(m%rcCvP)) return

      allocate (m%rcCvP  ( m % nRc,   2 ) )
      allocate (m%rcCv   ( m % mxnRc, 2 ) )
      allocate (m%rcFcP  ( m % nRc,   2 ) )
      allocate (m%rcFc   ( m % mxnRc    ) )
      allocate (m%rcFcOr ( m % mxnRc    ) )

      return      
      end subroutine alloc_mapping_rc 

**********************************************************************

      subroutine dealloc_mapping_rc (m)

      implicit none
      type(mapping), intent(in out) :: m

      if (.not.allocated(m%rcCvP)) return

      deallocate ( m % rcCvP  )
      deallocate ( m % rcCv   )
      deallocate ( m % rcFcP  )
      deallocate ( m % rcFc   )
      deallocate ( m % rcFcOr )

      return      
      end subroutine dealloc_mapping_rc 

***********************************************************************


      subroutine read_mapping (iun, mpg)
#ifdef WG_DEBUG
      use b2us_debug
#endif
      implicit none
      integer, intent(in) :: iun
      integer :: idum(0:9)
      type (mapping), intent (inout) :: mpg

      call alloc_mapping (mpg)

      call cfruin (iun, 1, idum, 'isClassicalGrid')
      mpg%isClassicalGrid = idum(0)
      call cfruin (iun, 3, idum, 'nx,ny,nncut')
      mpg%nx    = idum(0)
      mpg%ny    = idum(1)
      mpg%nncut = idum(2)
#ifdef WG_DEBUG
      nx_st = mpg%nx
      ny_st = mpg%ny
      isClassicalGrid = mpg%isClassicalGrid
#endif

      call cfruin (iun, mpg%nCv*2,   mpg%cvFcP, 'cvFcP')
      call cfruin (iun, mpg%nCmxFc,  mpg%cvFc,  'cvFc')
      call cfruin (iun, mpg%nFc*2,   mpg%fcCv,  'fcCv')
      call cfruin (iun, mpg%nFc*2,   mpg%fcVx,  'fcVx')
      call cfruin (iun, mpg%nCv*2,   mpg%cvVxP, 'cvVxP')
      call cfruin (iun, mpg%nCmxVx,  mpg%cvVx,  'cvVx')
      call cfruin (iun, mpg%nVx*2,   mpg%vxFcP, 'vxFcP')
      call cfruin (iun, mpg%nVmxFc,  mpg%vxFc,  'vxFc')
      call cfruin (iun, mpg%nVx*2,   mpg%vxCvP, 'vxCvP')
      call cfruin (iun, mpg%nVmxCv,  mpg%vxCv,  'vxCv')
      call cfruin (iun, mpg%nFs*2,   mpg%fsCvP, 'fsCvP')
      call cfruin (iun, mpg%nCv,     mpg%fsCv,  'fsCv')
      call cfruin (iun, mpg%nFs*2,   mpg%fsFcP, 'fsFcP')
      call cfruin (iun, mpg%nFc,     mpg%fsFc,  'fsFc')
      call cfruin (iun, mpg%nFc,     mpg%fcReg, 'fcReg')
      call cfruin (iun, mpg%nCv,     mpg%cvReg, 'cvReg')
      call cfrure (iun, mpg%nCmxFc,  mpg%intcellP,'intcellP')
      call cfrure (iun, mpg%nCmxFc,  mpg%intcellR,'intcellR')

      call init_mapping (mpg)

      return
      end subroutine read_mapping

***********************************************************************

      subroutine write_mapping (nout, mpg)

      implicit none
      integer, intent(in) :: nout
      integer :: idum(0:9)
      type (mapping), intent (in) :: mpg

      idum(0) = mpg%isClassicalGrid
      call cfwuin (nout, 1, idum, 'isClassicalGrid')
      if (mpg%isClassicalGrid == 1) then
        idum(0) = mpg%nx
        idum(1) = mpg%ny
        idum(2) = mpg%nncut
      else
        ! dummy values for unstructured grid
        idum(0) = -1
        idum(1) = -1
        idum(2) = -1
      endif      
      call cfwuin (nout, 3, idum, 'nx,ny,nncut')

      call cfwuin (nout, mpg%nCv*2,   mpg%cvFcP, 'cvFcP')
      call cfwuin (nout, mpg%nCmxFc,  mpg%cvFc,  'cvFc')
      call cfwuin (nout, mpg%nFc*2,   mpg%fcCv,  'fcCv')
      call cfwuin (nout, mpg%nFc*2,   mpg%fcVx,  'fcVx')
      call cfwuin (nout, mpg%nCv*2,   mpg%cvVxP, 'cvVxP')
      call cfwuin (nout, mpg%nCmxVx,  mpg%cvVx,  'cvVx')
      call cfwuin (nout, mpg%nVx*2,   mpg%vxFcP, 'vxFcP')
      call cfwuin (nout, mpg%nVmxFc,  mpg%vxFc,  'vxFc')
      call cfwuin (nout, mpg%nVx*2,   mpg%vxCvP, 'vxCvP')
      call cfwuin (nout, mpg%nVmxCv,  mpg%vxCv,  'vxCv')
      call cfwuin (nout, mpg%nFs*2,   mpg%fsCvP, 'fsCvP')
      call cfwuin (nout, mpg%nCv,     mpg%fsCv,  'fsCv')
      call cfwuin (nout, mpg%nFs*2,   mpg%fsFcP, 'fsFcP')
      call cfwuin (nout, mpg%nFc,     mpg%fsFc,  'fsFc')
      call cfwuin (nout, mpg%nFc,     mpg%fcReg, 'fcReg')
      call cfwuin (nout, mpg%nCv,     mpg%cvReg, 'cvReg')
      call cfwure (nout, mpg%nCmxFc,  mpg%intcellP,'intcellP')
      call cfwure (nout, mpg%nCmxFc,  mpg%intcellR,'intcellR')

      return
      end subroutine write_mapping

**********************************************************************

      subroutine init_mapping(m)
      implicit none
      type (mapping), intent (inout) :: m
      integer :: i, iCv, iVx, iVx1, iNv, iNv1, iNv2, iNvmx
      integer, allocatable :: indCv(:)
cwdk  This initialization routine precomputes a number of arrays
cwdk  in the mapping that are not stored in the b2fgmtry file, but 
cwdk  can be derived directly from data in that file

*   ..regions
      m%nnreg(0) = maxval(m%cvReg)
      m%nnreg(1) = maxval(m%fcReg)

*   ..identify cells on closed surfaces
      if (m%nnreg(0).eq.1) then
        write (*,*) ('WG-warning: periodic_bc not yet supported.')
        write (*,*) ('WG-warning: assuming periodic_bc .eq. 0.')
        m%cvOnClosedSurface = .false.
!WG_TODO        m%cvOnClosedSurface = periodic_bc.eq.1
      else if (m%nnreg(0).eq.2 .or. m%nnreg(0).eq.4 .or.
     .         m%nnreg(0).eq.5 .or. m%nnreg(0).eq.8) then
        m%cvOnClosedSurface = mod(m%cvReg,4).eq.1
      else
        call xerrab ('Unknown number of regions!')
      end if

*   ..set up connectivity information for the matrix stencil

      allocate (indCv(4*4+4+1)) ! max 4 vertices per cell, max 8 cells for each vertex (very conservative)

      iNv = 1
      do iCv = 1, m%nCv
        m%cvNvP(iCv,1) = iNv
        ! init: each cell appears in its own stencil
        indCv = iCv
        iNvmx = 1
        ! list all cells surrounding cell iCv, by listing the
        ! cells connecting to all vertices of cell iCv
        do iVx = 1, m%cvVxP(iCv,2)
          iVx1 = m%cvVx(m%cvVxP(iCv,1)+iVx-1)     ! the number in the vertex list
          indCv(iNvmx+1:iNvmx+m%vxCvP(iVx1,2)) =  ! all cells connecting to this vertex
     .      m%vxCv(m%vxCvP(iVx1,1):m%vxCvP(iVx1,1)+m%vxCvP(iVx1,2)-1)
          iNvmx = iNvmx + m%vxCvP(iVx1,2)
        enddo
        ! remove duplicates
        m%cvNv(m%cvNvP(iCv,1)) = indCv(1)
        i = 1
        outer: do iNv1 = 2, iNvmx
          do iNv2 = 1, i 
            if (m%cvNv(m%cvNvP(iCv,1)+iNv2-1).eq.indCv(iNv1)) then
              ! duplicate element
              cycle outer
            endif
          enddo
          ! new element found
          i = i + 1
          m%cvNv(m%cvNvP(iCv,1)+i-1) = indCv(iNv1)
        end do outer
        ! TODO: sort by cell number, keeping first cell fixed (main diagonal)
        m%cvNvP(iCv,2) = i
        iNv = iNv + i
      enddo

      deallocate(indCv)



      return
      end subroutine init_mapping

**********************************************************************

      subroutine get_version(boundary_filename,version)
      implicit none
      logical :: file_ok
      character*260 boundary_filename
      character*10 version
      character*260 filename
      filename=trim(boundary_filename)
      call find_file(filename,file_ok)
      if(file_ok) then
        open(99,file=filename)
        call cfverr(99,version)
        close(99)
      else
        call xerrab ('No '//trim(filename))
      endif
      end subroutine get_version

      end module b2us_map
