      module b2mod_indirect
      use b2mod_types
      use carre_constants
      use tradui_constants
      use b2mod_cellhelper
      use b2mod_connectivity
      use b2mod_dimensions
      use b2mod_math
      implicit none

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
     *  nncut, ccut, nnreg(0:2), periodic_bc, nncutmax, redef_gmtry
      integer, allocatable, save ::
     *  cflags(:,:,:)
      logical, allocatable, save ::
     *  on_closed_surface(:,:)

      logical, save :: b2mod_indirect_allocated = .false.

      contains

*----------------------------------------------------------------------

      subroutine alloc_b2mod_indirect(nxd,nyd,nncutmax)
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
      nncutmax=max(1,DEF_NCUT)
      allocate(leftcut(1:nncutmax))
      allocate(rightcut(1:nncutmax))
      allocate(topcut(1:nncutmax))
      allocate(bottomcut(1:nncutmax))
      allocate(cflags(-1:nxd,-1:nyd,GRID_N_CELLFLAGS))
      bottomcut   = nyd+1
      topcut      = -2
      rightcut    = -2
      leftcut     = nxd+1
      ccut        = 0
      redef_gmtry = 0
      cflags      = 0
      b2mod_indirect_allocated = .true.

      return
      end subroutine alloc_b2mod_indirect

      subroutine dealloc_b2mod_indirect
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
      end subroutine dealloc_b2mod_indirect

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

! procedures
      external xertst

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

      end module b2mod_indirect
