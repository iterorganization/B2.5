      module b2mod_indirect
      use b2mod_types
      use carre_constants
      use tradui_constants
      use b2mod_cellhelper
      use b2mod_connectivity
      implicit none
#include "DIMENSIONS.F"

      integer, parameter :: NNREGMAX=8

      integer, allocatable, save ::
     *  leftix(:,:),
     *  leftiy(:,:),
     *  rightix(:,:),
     *  rightiy(:,:),
     *  topix(:,:),
     *  topiy(:,:),
     *  bottomix(:,:),
     *  bottomiy(:,:),
     *  region(:,:,:),
     *  resignore(:,:,:),
     *  leftcut(:),
     *  rightcut(:),
     *  topcut(:),
     *  bottomcut(:)
      integer, save ::
     *  nncut, nnreg(0:2), periodic_bc, nncutmax, redef_gmtry
      integer, allocatable, save ::
     *  cflags(:,:,:)
      logical, allocatable, save ::
     *  on_closed_surface(:,:)


      type, public :: mapping_ag

        integer :: nCv, nFc, nVx, nCg, nCi, nCmxVx, nCmxFc, nVmxCv,
     &             nVmxFc, nFs, nFt, nBc, mxnBc, nRc, mxnRc, nCmxNv, 
     &             nCf, mxnCf, nx, ny, nFmxCv
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
!       integer, allocatable :: ftCvP(:,:)             ! (nFt,2) pointing for every flux tube to the first index number in the ftCv list and its number of cells
!        integer, allocatable :: ftCv(:)                !         listing for every flux tube the corresponding cell numbers
!        integer, allocatable :: ftFcP(:,:)             ! (nFt,2) pointing for every flux tube to the first index number in the ftFc list and its number of cell faces
!        integer, allocatable :: ftFc(:)                !         listing for every flux tube the corresponding cell face numbers
!        integer, allocatable :: cvFt(:)                !         listing for every cell the corresponding flux tube
!        integer, allocatable :: fsFcP(:,:)             ! (nFs,2) pointing for every flux surface to the first index number in the fsFc list and its number of cell faces
!        integer, allocatable :: fsFc(:)                !         listing for every flux surface the corresponding cell face numbers
!        integer, allocatable :: fcReg(:)               ! Region numbers of the cell faces 
!        integer, allocatable :: cvReg(:)               ! Region numbers of the cells
!        integer, allocatable :: ftReg(:)               ! Region numbers of the flux tubes
!        logical, allocatable :: cvOnClosedSurface(:)   ! Logical to indicate whether a cell is on a closed flux surface

!        integer, allocatable :: bcCvP(:,:)             ! (nBc,2) pointing for each boundary condition to the first index in the bcCv list, and its number of guard cells
!        integer, allocatable :: bcCv(:,:)              !         listing for each boundary condition the guard cells, and the corresponding domain cells
!        integer, allocatable :: bcFcP(:,:)             ! (nBc,2) pointing for each boundary condition to the first index in the bcFc list, and its number of boundary faces
!        integer, allocatable :: bcFc(:)                !         listing for each boundary condition the boundary faces
!        real (kind=R8), allocatable :: bcFcOr (:)      ! Multiplier to indicate whether outward normal coincides with the face normal (1.0_R8) or not (-1.0_R8)

!        integer, allocatable :: rcCvP(:,:)             ! (nRc,2) pointing for each stratum to the first index in the rcCv list, and its number of guard cells
!        integer, allocatable :: rcCv(:,:)              !         listing for each stratum the guard cells, and the corresponding domain cells
!        integer, allocatable :: rcFcP(:,:)             ! (nRc,2) pointing for each stratum to the first index in the rcFc list, and its number of boundary faces
!        integer, allocatable :: rcFc(:)                !         listing for each stratum the boundary faces
!        real (kind=R8), allocatable :: rcFcOr (:)      ! Multiplier to indicate whether outward normal coincides with the face normal (1.0_R8) or not (-1.0_R8)

!        integer, allocatable :: cfCvP(:,:)             ! (nCf,2) pointing for each cost function to the first index in the cfCv list, and its number of cells
!        integer, allocatable :: cfCv(:)                !         listing for each cost function the corresponding domain cells
!        integer, allocatable :: cfFcP(:,:)             ! (nCf,2) pointing for each cost function to the first index in the cfFc list, and its number of faces
!        integer, allocatable :: cfFc(:)                !         listing for each cost function the corresponding faces
!        real (kind=R8), allocatable :: cfFcOr (:)      ! Multiplier to indicate whether outward normal coincides with the face normal (1.0_R8) or not (-1.0_R8)

!        real (kind=R8), allocatable :: intfaceV(:,:)   ! volume factors for interpolation from centers to faces; currently: geo%fcVol
!        real (kind=R8), allocatable :: intvertexV(:)   ! volume factors for interpolation from centers to vertices; currently: geo%vxVol
!        real (kind=R8), allocatable :: intcellP(:)     ! factors for interpolation from faces to centers in poloidal direction
!        real (kind=R8), allocatable :: intcellR(:)     ! factors for interpolation from faces to centers in radial direction

!        integer, allocatable :: cvNvP(:,:)             !(nCv,2): pointing for each cell to the start in cvNv array, and the number of points in its stencil
!        integer, allocatable :: cvNv(:)                !         listing for each cell the numbers of the cells in its stencil


        integer :: isClassicalGrid                     ! indicates whether grid is derived from an "original", structured grid
!        integer :: nncut, nctcnt                        ! some parameters related to the original grid
        integer, allocatable :: imapCv(:,:)            ! back projection of cv indices to structured grid
        integer, allocatable :: imapFcx(:,:)           ! back projection of fc indices for x-faces to structured grid
        integer, allocatable :: imapFcy(:,:)           ! back projection of fc indices for y-faces to structured grid
        integer, allocatable :: imapVx(:,:)            ! back projection of Vx indices to structured grid 

        integer, allocatable :: icornVx(:)             ! corner indices; to be used in converting target information
        integer, allocatable :: ictcntrl(:)            ! indices of central cuts; to be used in converting triangle information

!        integer, allocatable :: fcLbl(:)               ! (nFc) indicates for each face the corresponding face label
!        integer, allocatable :: ftLbl(:)               ! (nFt) indicates for each flux tube the corresponding label
!        integer, allocatable :: cvLbl(:)               ! (nCv) indicates for each control volume the corresponding label
        end type mapping_ag

      logical, save :: b2mod_indirect_allocated = .false.

      contains

*----------------------------------------------------------------------
      subroutine alloc_b2mod_indirect(m,nCv,nFc,nVx,nCmxVx,nCmxFc,
     &        nFmxCv, nVmxCv, nVmxFc, nx, ny)
      implicit none
      
      integer, intent(in) :: nCv, nFc, nVx, nCmxVx, nCmxFc, nFmxCv,
     &          nVmxCv, nVmxFc ,nx,ny
      type(mapping_ag), intent(inout) :: m

      if (b2mod_indirect_allocated) return
	   allocate (m%cvFcP(nCv,1:2)) 
      allocate (m%cvFc(nCmxFc))
      allocate (m%fcCv(nFc,1:2))
      allocate (m%fcVx(nFc,1:2))
      allocate (m%cvVxP(nCv,1:2))
      allocate (m%cvVx(nCmxVx))
      allocate (m%vxFcP(nVx,1:2)) 
      allocate (m%vxFc(nVmxFc))
      allocate (m%vxCvP(nVx,1:2))  
      allocate (m%vxCv(nVmxCv)) 
!      allocate (m%ftCv(nCv))
!      allocate (m%ftCvP(nFt,1:2))
!      allocate (m%ftFc(nFc))
!      allocate (m%ftFcP(nFt,1:2))
!      allocate (m%cvFt(nCv))
!      allocate (m%fsFc(nFc))
!      allocate (m%fsFcP(nFs,1:2))
!      allocate (m%fcReg(nFc))
!      allocate (m%cvReg(nCv))
!      allocate (m%ftReg(nFt))
!      allocate (m%cvOnClosedSurface(nCv))
      
!      allocate (m%intfaceV(nFc,1:2))
!      allocate (m%intvertexV(nCmxVx * nCv))
!      allocate (m%intcellP(nCmxFc))
!      allocate (m%intcellR(nCmxFc))

!      allocate (m%cvNvP(nCv,1:2))
!      allocate (m%cvNv(nCmxNv))

      if (m%isClassicalGrid.eq.1) then
        allocate (m%imapCv(-1:nx,-1:ny))
        allocate (m%imapFcx(-1:nx,-1:ny))
        allocate (m%imapFcy(-1:nx,-1:ny))
        allocate (m%imapVx(-1:nx,-1:ny))
        allocate (m%icornVx(1:nx))
#ifdef WG_DEBUG
!        allocate (imapCv(-1:nx,-1:ny))
!        allocate (imapFcx(-1:nx,-1:ny))
!        allocate (imapFcy(-1:nx,-1:ny))
!        allocate (imapVx(-1:nx,-1:ny))
#endif
      else
        allocate (m%imapCv(0,0))
        allocate (m%imapFcx(0,0))
        allocate (m%imapFcy(0,0))
        allocate (m%imapVx(0,0))
        allocate (m%icornVx(0))
#ifdef WG_DEBUG
!        allocate (imapCv(0,0))
!        allocate (imapFcx(0,0))
!        allocate (imapFcy(0,0))
!        allocate (imapVx(0,0))
#endif
      endif

!      allocate (m%fcLbl(nFc))
!      allocate (m%ftLbl(nFt))
!      allocate (m%cvLbl(nCv))


      call initialize_mapping_ag (m)

!      write (*,*) '=================================='
!      write (*,*) 'Mapping allocated'
!      write (*,*) 'nCv, nCi, nCg', m%nCv, m%nCi, m%nCg
!      write (*,*) 'nFc, nVx, nFt', m%nFc, m%nVx, m%nFt
!      write (*,*) 'nFs'          , m%nFs
!      write (*,*) 'nCmxVx, nCmxFc', m%nCmxVx, m%nCmxFc
!      write (*,*) 'nVmxCv, nVmxFc', m%nVmxCv, m%nVmxFc
!      write (*,*) '=================================='



      b2mod_indirect_allocated = .true.

      return
      end subroutine alloc_b2mod_indirect

*----------------------------------------------------------------------
      subroutine initialize_mapping_ag(m)

      implicit none
      type(mapping_ag), intent(inout) :: m

      m%cvFcP = 0 
      m%cvFc = 0
      m%fcCv = 0
      m%fcVx = 0
      m%cvVxP = 0
      m%cvVx = 0
      m%vxFcP = 0 
      m%vxFc = 0
      m%vxCvP = 0  
      m%vxCv = 0 
!      m%ftCv = 0
!      m%ftCvP = 0
!      m%ftFc = 0
!      m%ftFcP = 0
!      m%cvFt = 0
!      m%fcReg = 0
!      m%cvReg = 0
!      m%ftReg = 0
!      m%cvOnClosedSurface = .false.
      
!      m%intfaceV = 0._R8
!      m%intvertexV = 0._R8
!      m%intcellP = 0._R8
!      m%intcellR = 0_R8

!      m%cvNvP = 0
!      m%cvNv = 0

      if (m%isClassicalGrid.eq.1) then
        m%imapCv = 0
        m%imapFcx = 0
        m%imapFcy = 0
        m%imapVx = 0
        m%icornVx = 0
#ifdef WG_DEBUG
!        imapCv = 0
!        imapFcx = 0
!        imapFcy = 0
!        imapVx = 0
#endif
       endif

!      m%fcLbl = 0
!      m%ftLbl = 0
!      m%cvLbl = 0

!      write (*,*) '=================================='
!      write (*,*) 'Mapping initialized'
!      write (*,*) '=================================='

      return

      end subroutine initialize_mapping_ag

*----------------------------------------------------------------------    
      subroutine dealloc_b2mod_indirect(m)

      implicit none
      type(mapping_ag), intent(in out) :: m
      

      if (b2mod_indirect_allocated) then

      deallocate (m%cvFcP) 
      deallocate (m%cvFc)
      deallocate (m%fcCv)
      deallocate (m%fcVx)
      deallocate (m%cvVxP)
      deallocate (m%cvVx)
      deallocate (m%vxFcP) 
!      deallocate (m%vxFc) 
      deallocate (m%vxCvP)  
!      deallocate (m%vxCv) 
!      deallocate (m%ftCv)
!      deallocate (m%ftCvP)
!      deallocate (m%ftFc)
!      deallocate (m%ftFcP)
!      deallocate (m%cvFt)
!      deallocate (m%fcReg)
!      deallocate (m%cvReg)
!      deallocate (m%ftReg)
!      deallocate (m%cvOnClosedSurface)

!      deallocate (m%intfaceV) 
!      deallocate (m%intvertexV) 
!      deallocate (m%intcellP) 
!      deallocate (m%intcellR)   

!      deallocate (m%cvNvP)
!      deallocate (m%cvNv)  
    
      deallocate (m%imapCv)
      deallocate (m%imapFcx)  
      deallocate (m%imapFcy)  
      deallocate (m%imapVx)  
      deallocate (m%icornVx)
!      if (allocated(m%ictcntrl)) deallocate (m%ictcntrl)
#ifdef WG_DEBUG
!     deallocate (imapCv)
!     deallocate (imapFcx)
!      deallocate (imapFcy)
!     deallocate (imapVx)
#endif

!      deallocate (m%fcLbl)
!      deallocate (m%ftLbl)
!      deallocate (m%cvLbl)

	   b2mod_indirect_allocated = .false.

      endif

      return
      end subroutine dealloc_b2mod_indirect

*----------------------------------------------------------------------
      subroutine alloc_b2mod_indirect_st(nxd,nyd,nncutmax)
      implicit none
      integer nxd,nyd,nncutmax

      if (b2mod_indirect_allocated) return
      allocate(leftix(-1:nxd,-1:nyd))
      allocate(leftiy(-1:nxd,-1:nyd))
      allocate(rightix(-1:nxd,-1:nyd))
      allocate(rightiy(-1:nxd,-1:nyd))
      allocate(topix(-1:nxd,-1:nyd))
      allocate(topiy(-1:nxd,-1:nyd))
      allocate(bottomix(-1:nxd,-1:nyd))
      allocate(bottomiy(-1:nxd,-1:nyd))
      allocate(region(-1:nxd,-1:nyd,0:2))
      allocate(on_closed_surface(-1:nxd,-1:nyd))
      nnreg=0
      allocate(resignore(-1:nxd,-1:nyd,2))
      nncutmax=DEF_NCUT
      allocate(leftcut(1:nncutmax))
      allocate(rightcut(1:nncutmax))
      allocate(topcut(1:nncutmax))
      allocate(bottomcut(1:nncutmax))
      allocate(cflags(-1:nxd,-1:nyd,GRID_N_CELLFLAGS))
      bottomcut = nyd+1
      topcut    = -2
      rightcut  = -2
      leftcut   = nxd+1
      redef_gmtry = 0
      cflags    = 0
      b2mod_indirect_allocated = .true.

      return
      end subroutine alloc_b2mod_indirect_st


      subroutine dealloc_b2mod_indirect_st
      implicit none

      if (b2mod_indirect_allocated) then
        deallocate(leftix)
        deallocate(leftiy)
        deallocate(rightix)
        deallocate(rightiy)
        deallocate(topix)
        deallocate(topiy)
        deallocate(bottomix)
        deallocate(bottomiy)
        deallocate(region)
        deallocate(on_closed_surface)
        deallocate(resignore)
        deallocate(leftcut)
        deallocate(rightcut)
        deallocate(topcut)
        deallocate(bottomcut)
        deallocate(cflags)
        b2mod_indirect_allocated = .false.
      endif

      return
      end subroutine dealloc_b2mod_indirect_st

*----------------------------------------------------------------------
! For the ghost cell at position (ix,iy), get the face index
! on which the neighbour internal cell is connected
      integer function ghostGetBndFace(nx, ny, ix, iy)
      use b2mod_connectivity
      implicit none
      integer, intent(in) :: nx, ny, ix, iy
      
! internal
      integer :: nix, niy, iDir

      if (.not. isGhostCell(cflags(ix,iy,CELLFLAG_TYPE))) then
         ghostGetBndFace = NODIRECTION
         return
      end if

      do iDir = LEFT, TOP
         select case (iDir)
         case (LEFT)
            nix = leftix(ix,iy)
            niy = leftiy(ix,iy)
         case (BOTTOM)
            nix = bottomix(ix,iy)
            niy = bottomiy(ix,iy)
         case (RIGHT)
            nix = rightix(ix,iy)
            niy = rightiy(ix,iy)
         case (TOP)
            nix = topix(ix,iy)
            niy = topiy(ix,iy)
         end select

         if (isInDomain(nx, ny, nix, niy)) then
            if (cflags(nix, niy, CELLFLAG_TYPE) == GRID_BOUNDARY) then
               ghostGetBndFace = iDir
               return
            end if
         end if
      end do

! Ghost cell without internal neighbour (e.g. corner ghost cell)
      ghostGetBndFace = NODIRECTION
      return
      end function ghostGetBndFace

! Build a list of cell indices for all ghost cells located 
! at the given range of structures 
      subroutine build_boundary_list( structFrom, structTo,
     &   listLen, x, y, type )
      use b2mod_connectivity
      implicit none
      integer, intent(in) :: structFrom, structTo
      integer, intent(out) :: listLen
      integer, intent(out), dimension(:) :: x, y
      character, intent(out), dimension(:) :: type
      
! internal
      integer :: ix, iy, nx, ny, gix, giy, iFace, iStruct
      character :: ctype

      nx = ubound(cflags, 1)
      ny = ubound(cflags, 2)

      listLen = 0

      do ix = -1, nx
         do iy = -1, ny

            if (cflags(ix, iy, CELLFLAG_TYPE) /= GRID_BOUNDARY) cycle

            do iFace = CELLFLAG_LEFTFACE, CELLFLAG_TOPFACE

               iStruct = cflags(ix, iy, iFace)

               if ( iStruct == GRID_UNDEFINED ) cycle
               if ( iStruct < structFrom ) cycle
               if ( iStruct > structTo ) cycle

               select case (iFace)
               case(CELLFLAG_LEFTFACE)
                  gix = leftix(ix, iy)
                  giy = leftiy(ix, iy)
                  ctype = "W"
               case(CELLFLAG_BOTTOMFACE)
                  gix = bottomix(ix, iy)
                  giy = bottomiy(ix, iy)                
                  ctype = "S"
               case(CELLFLAG_RIGHTFACE)
                  gix = rightix(ix, iy)
                  giy = rightiy(ix, iy)                
                  ctype = "E"
               case(CELLFLAG_TOPFACE)
                  gix = topix(ix, iy)
                  giy = topiy(ix, iy)                
                  ctype = "N"
               end select

! A vanishing face has no ghost cell -> no neighbour on that side
               if (.not. isInDomain(nx, ny, gix, giy)) cycle

               call xertst(isGhostCell(cflags(gix, giy, CELLFLAG_TYPE)),
     &          "Broken cell type: expected ghost cell" )

               listLen = listLen + 1
               x(listLen) = gix
               y(listLen) = giy
               type(listLen) = ctype

            end do

         end do
      end do
      return
      end subroutine build_boundary_list

      integer function FindTopRealCell(ny,ix)
      implicit none
      integer ix, ny
      integer iy

      iy = ny
      do while (.not.isRealCell(cflags(ix, iy, CELLFLAG_TYPE)))
        iy = iy - 1
        if (iy.eq.-2) exit
      end do
      FindTopRealCell = iy
      return
      end function FindTopRealCell

      integer function FindTopCell(ny,ix)
      implicit none
      integer ix, ny
      integer iy
      logical found

      iy = ny
      found = .false.
      do while (.not.found)
        found = isRealCell(cflags(ix,iy,CELLFLAG_TYPE))
        found = found .or. 
     &   (isGhostCell(cflags(ix,iy,CELLFLAG_TYPE)) .and.
     &    cflags(ix,iy,CELLFLAG_BOTTOMFACE) /= GRID_UNDEFINED .and.
     &    bottomix(ix,iy).eq.ix .and. bottomiy(ix,iy).eq.iy-1)
        found = found .or. 
     &   (isGhostCell(cflags(ix,iy,CELLFLAG_TYPE)) .and.
     &    cflags(ix,iy,CELLFLAG_LEFTFACE) /= GRID_UNDEFINED .and.
     &    leftiy(ix,iy).eq.iy .and. leftix(ix,iy).eq.ix-1)
        found = found .or. 
     &   (isGhostCell(cflags(ix,iy,CELLFLAG_TYPE)) .and.
     &    cflags(ix,iy,CELLFLAG_RIGHTFACE) /= GRID_UNDEFINED .and.
     &    rightiy(ix,iy).eq.iy .and. rightix(ix,iy).eq.ix+1)
        if (.not.found) iy = iy - 1
        if (iy.eq.-2) found = .true.
      end do
      FindTopCell = iy
      return
      end function FindTopCell

      integer function FindBottomRealCell(ny,ix)
      implicit none
      integer ix, ny
      integer iy

      iy = -1
      do while (.not.isRealCell(cflags(ix, iy, CELLFLAG_TYPE)))
        iy = iy + 1
        if (iy.eq.ny+1) exit
      end do
      FindBottomRealCell = iy
      return
      end function FindBottomRealCell

      integer function FindBottomCell(ny,ix)
      implicit none
      integer ix, ny
      integer iy
      logical found

      iy = -1
      found = .false.
      do while (.not.found)
        found = isRealCell(cflags(ix,iy,CELLFLAG_TYPE))
        found = found .or. 
     &   (isGhostCell(cflags(ix,iy,CELLFLAG_TYPE)) .and.
     &    cflags(ix,iy,CELLFLAG_TOPFACE) /= GRID_UNDEFINED .and.
     &    topix(ix,iy).eq.ix .and. topiy(ix,iy).eq.iy+1)
        found = found .or. 
     &   (isGhostCell(cflags(ix,iy,CELLFLAG_TYPE)) .and.
     &    cflags(ix,iy,CELLFLAG_LEFTFACE) /= GRID_UNDEFINED .and.
     &    leftiy(ix,iy).eq.iy .and. leftix(ix,iy).eq.ix-1)
        found = found .or. 
     &   (isGhostCell(cflags(ix,iy,CELLFLAG_TYPE)) .and.
     &    cflags(ix,iy,CELLFLAG_RIGHTFACE) /= GRID_UNDEFINED .and.
     &    rightiy(ix,iy).eq.iy .and. rightix(ix,iy).eq.ix+1)
        if (.not.found) iy = iy + 1
        if (iy.eq.ny+1) found = .true.
      end do
      FindBottomCell = iy
      return
      end function FindBottomCell

      integer function FindRightRealCell(nx,iy)
      implicit none
      integer iy, nx
      integer ix

      ix = nx
      do while (.not.isRealCell(cflags(ix, iy, CELLFLAG_TYPE)))
        ix = ix - 1
        if (ix.eq.-2) exit
      end do
      FindRightRealCell = ix
      return
      end function FindRightRealCell

      integer function FindRightCell(nx,iy)
      implicit none
      integer iy, nx
      integer ix
      logical found

      ix = nx
      found = .false.
      do while (.not.found)
        found = isRealCell(cflags(ix,iy,CELLFLAG_TYPE))
        found = found .or. 
     &   (isGhostCell(cflags(ix,iy,CELLFLAG_TYPE)) .and.
     &    cflags(ix,iy,CELLFLAG_LEFTFACE) /= GRID_UNDEFINED .and.
     &    leftiy(ix,iy).eq.iy .and. leftix(ix,iy).eq.ix-1)
        found = found .or. 
     &   (isGhostCell(cflags(ix,iy,CELLFLAG_TYPE)) .and.
     &    cflags(ix,iy,CELLFLAG_BOTTOMFACE) /= GRID_UNDEFINED .and.
     &    bottomix(ix,iy).eq.ix .and. bottomiy(ix,iy).eq.iy-1)
        found = found .or. 
     &   (isGhostCell(cflags(ix,iy,CELLFLAG_TYPE)) .and.
     &    cflags(ix,iy,CELLFLAG_TOPFACE) /= GRID_UNDEFINED .and.
     &    topix(ix,iy).eq.ix .and. topiy(ix,iy).eq.iy+1)
        if (.not.found) ix = ix - 1
        if (ix.eq.-2) found = .true.
      end do
      FindRightCell = ix
      return
      end function FindRightCell

      integer function FindLeftRealCell(nx,iy)
      implicit none
      integer iy, nx
      integer ix

      ix = -1
      do while (.not.isRealCell(cflags(ix, iy, CELLFLAG_TYPE)))
        ix = ix + 1
        if (ix.eq.nx+1) exit 
      end do
      FindLeftRealCell = ix
      return
      end function FindLeftRealCell

      integer function FindLeftCell(nx,iy)
      implicit none
      integer iy, nx
      integer ix
      logical found

      ix = -1
      found = .false.
      do while (.not.found)
        found = isRealCell(cflags(ix,iy,CELLFLAG_TYPE))
        found = found .or. 
     &   (isGhostCell(cflags(ix,iy,CELLFLAG_TYPE)) .and.
     &    cflags(ix,iy,CELLFLAG_RIGHTFACE) /= GRID_UNDEFINED .and.
     &    rightiy(ix,iy).eq.iy .and. rightix(ix,iy).eq.ix+1)
        found = found .or. 
     &   (isGhostCell(cflags(ix,iy,CELLFLAG_TYPE)) .and.
     &    cflags(ix,iy,CELLFLAG_BOTTOMFACE) /= GRID_UNDEFINED .and.
     &    bottomix(ix,iy).eq.ix .and. bottomiy(ix,iy).eq.iy-1)
        found = found .or. 
     &   (isGhostCell(cflags(ix,iy,CELLFLAG_TYPE)) .and.
     &    cflags(ix,iy,CELLFLAG_TOPFACE) /= GRID_UNDEFINED .and.
     &    topix(ix,iy).eq.ix .and. topiy(ix,iy).eq.iy+1)
        if (.not.found) ix = ix + 1
        if (ix.eq.nx+1) found = .true.
      end do
      FindLeftCell = ix
      return
      end function FindLeftCell

  !> Area of triangle given by vertex coordinates
      function triangleArea( x0, y0, x1, y1, x2, y2 )
      implicit none
      real (kind=R8) :: triangleArea
      real (kind=R8) :: x0, y0, x1, y1, x2, y2
      intrinsic abs

      triangleArea = 0.5_R8 * abs( (x1-x0)*(y2-y0) - (y1-y0)*(x2-x0) )

      end function triangleArea

  !> Area of quadrangle given by vertex coordinates
  !> The vertices must be in sequence! (=> structured)
      function quadArea_st( x0, y0, x1, y1, x2, y2, x3, y3 )
      implicit none
      real (kind=R8) :: quadArea_st
      real (kind=R8) :: x0, y0, x1, y1, x2, y2, x3, y3
      intrinsic abs

      quadArea_st = TriangleArea(x0, y0, x1, y1, x2, y2) +
     .           TriangleArea(x2, y2, x3, y3, x0, y0)

      end function quadArea_st
      
  !> Area of quadrangle given by vertex coordinates
  !> The vertices can be in random sequence
  !> Method not 100% waterproof for very sqeeuzed and distorted cell
      function quadArea (x0, y0, x1, y1, x2, y2, x3, y3 )
      implicit none
      real (kind=R8) :: A1, A2, A3, A4, quadArea
      real (kind=R8) :: x0, y0, x1, y1, x2, y2, x3, y3
      intrinsic abs

      A1 = quadArea_sub(x0,y0,x1,y1,x2,y2,x3,y3)
      A2 = quadArea_sub(x1,y1,x2,y2,x3,y3,x0,y0)
      A3 = quadArea_sub(x2,y2,x3,y3,x0,y0,x1,y1)
      A4 = quadArea_sub(x3,y3,x0,y0,x1,y1,x2,y2)
         
      !Determine which two A's are equal, these are the correct ones
      if ((abs(A1 - A2)/A1 .lt. 1e-7) .or.
     &    (abs(A1 - A3)/A1 .lt. 1e-7) .or.
     &    (abs(A1 - A4)/A1 .lt. 1e-7)) then
               quadArea = A1;
      elseif ((abs(A2 - A3)/A2 .lt. 1e-7) .or.
     &    (abs(A2 - A4)/A2 .lt. 1e-7)) then
               quadArea = A2;
      elseif (abs(A3 - A4)/A3 .lt. 1e-7) then
               quadArea = A3;
      endif
        
      end function quadArea 


      function quadArea_sub (x0, y0, x1, y1, x2, y2, x3, y3 )
      implicit none
      real (kind=R8) :: quadArea_sub, Area1, Area2
      real (kind=R8) :: x0, y0, x1, y1, x2, y2, x3, y3
      real (kind=R8) :: fcx01, fcy01, fcx12, fcy12, fcx02, fcy02, 
     &    cvX, cvY
      real (kind=R8) :: dist01, dist12, dist02, min_dist
      intrinsic abs, min, sqrt

      Area1 = TriangleArea(x0, y0, x1, y1, x2, y2)
      
      cvX = 0.25_R8*(x0+x1+x2+x3)
      cvY = 0.25_R8*(y0+y1+y2+y3)
      ! Determine points of second triangle
      fcx01 = 0.5_R8 * (x0 + x1) 
      fcy01 = 0.5_R8 * (y0 + y1)

      fcx12 = 0.5_R8 * (x1 + x2)
      fcy12 = 0.5_R8 * (y1 + y2)

      fcx02 = 0.5_R8 * (x0 + x2)
      fcy02 = 0.5_R8 * (y0 + y2)

      dist01 = sqrt( (fcx01- x3)**2 + (fcy01- y3)**2 ) + 
     &      sqrt( (fcx01 - cvX)**2 + (fcy01 - cvY)**2 )
      dist12 = sqrt( (fcx12- x3)**2 + (fcy12- y3)**2 ) + 
     &      sqrt( (fcx12 - cvX)**2 + (fcy12 - cvY)**2 )     
      dist02 = sqrt( (fcx02- x3)**2 + (fcy02- y3)**2 ) + 
     &      sqrt( (fcx02 - cvX)**2 + (fcy02 - cvY)**2 )      

      min_dist = min(dist01,dist12,dist02)

      if (abs(min_dist - dist01)/min_dist .lt. 1e-7) then
         Area2 = TriangleArea(x0, y0, x1, y1, x3, y3)
         
      elseif (abs(min_dist - dist12)/min_dist .lt. 1e-7) then
         Area2 = TriangleArea(x1, y1, x2, y2, x3, y3)

      elseif (abs(min_dist - dist02)/min_dist .lt. 1e-7) then
         Area2 = TriangleArea(x0, y0, x2, y2, x3, y3)

      endif


      quadArea_sub = Area1 + Area2

      end function quadArea_sub



  !> Coordinates of centroid of a triangle given by vertex coordinates
      function triangleCentroid( x0, y0, x1, y1, x2, y2 )
      implicit none
      real (kind=R8) :: triangleCentroid(0:1)
      real (kind=R8) :: x0, y0, x1, y1, x2, y2

      triangleCentroid(0) = ( x0 + x1 + x2 ) / 3.0_R8
      triangleCentroid(1) = ( y0 + y1 + y2 ) / 3.0_R8
      return
      end function triangleCentroid

  !> Coordinates of centroid of a quadrangle given by vertex coordinates
  !> Corners are numbered in the B2 convention 2-3
  !>                                           0-1
      function quadCentroid( x0, y0, x1, y1, x2, y2, x3, y3 )
      implicit none
      real (kind=R8) :: quadCentroid(0:1)
      real (kind=R8) :: x0, y0, x1, y1, x2, y2, x3, y3

      quadCentroid = ( triangleCentroid(x0,y0,x1,y1,x3,y3) 
     &   * triangleArea(x0,y0,x1,y1,x3,y3) 
     &   + triangleCentroid(x0,y0,x2,y2,x3,y3) 
     &   * triangleArea(x0,y0,x2,y2,x3,y3) ) 
     &   / ( triangleArea(x0,y0,x1,y1,x3,y3) 
     &   + triangleArea(x0,y0,x2,y2,x3,y3) )
      return
      end function quadCentroid

  !> Average of a quantity given on the vertices of
  !> an arbitrary quadrangular cell
  !> Corners are numbered in the B2 convention 2-3
  !>                                           0-1
      function averageVertex( f0, x0, y0, f1, x1, y1, 
     &                        f2, x2, y2, f3, x3, y3 )
      implicit none
      real (kind=R8) :: averageVertex
      real (kind=R8) :: f0, f1, f2, f3
      real (kind=R8) :: x0, y0, x1, y1, x2, y2, x3, y3
      real (kind=R8) :: a0, a1, a2, a3
  
      a0 = TriangleArea( x2, y2, x0, y0, x1, y1 )
      a1 = TriangleArea( x0, y0, x1, y1, x3, y3 )
      a2 = TriangleArea( x0, y0, x2, y2, x3, y3 )
      a3 = TriangleArea( x2, y2, x3, y3, x1, y1 )
 
      averageVertex = (f0*a0 + f1*a1 + f2*a2 + f3*a3) /  
     &                (a0 + a1 + a2 + a3)
 
      return
      end function averageVertex


      end module b2mod_indirect
