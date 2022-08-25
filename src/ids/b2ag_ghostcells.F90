module b2ag_ghostcells

  use b2mod_types
  use carre_constants
  use tradui_constants
  use b2mod_connectivity
  use b2mod_cellhelper
  use Logging
  use Helper
  use b2mod_constants

  implicit none

  !integer :: nrid, nCi, nCg, nC

contains

  ! Given the grid dimensions without ghost cells and the number of cuts,
  ! figures out what the actual grid size nx, ny for the B2 data structures has to be
  ! accounting for ghost cells
  subroutine computeGridSizeWithGhostCells_st(nnx, nny, niso, nx, ny)
    integer, intent(in) :: nnx, nny, niso
    integer, intent(out) :: nx, ny

    nx = nnx + 2 * niso
    ny = nny
  end subroutine computeGridSizeWithGhostCells_st

  subroutine computeGridSizeWithGhostCells(nrid,nCi,nCg,nFc,nCmxVx0,nCmxFc0, &
              &  nFmxCv0, nVmxCv0, nVmxFc0, nCv, nCmxVx, nCmxFc, &
              &  nFmxCv, nVmxCv, nVmxFc, nCmxNv ) 
    integer, intent(in) :: nrid, nCi, nCg, nFc, nCmxVx0, nCmxFc0, nFmxCv0, &
              &  nVmxCv0, nVmxFc0
    integer, intent(out) :: nCv, nCmxVx, nCmxFc, nFmxCv, &
              &   nVmxCv, nVmxFc, nCmxNv
    
    
     
    !internal
    integer :: i, cflags(1:nCi), nBc, cvFc(1:nCmxFc0), iFc, nBF
    character id*32

    !With the assumption that each boundary cell has only one boundary face
    !read the cflag van cell information, count the number of 3's
    !read(nrid)
   
    call read_cvFc(nrid,nCv,cvFc,nCmxFc0)
    !count how many faces only appear once in cvFc, this are boundary faces (nBF)
    nBF = 0
    do iFc = 1,nFc
      if (count (cvFc == iFc) .eq. 1) then
         nBF = nBF + 1
      endif
    enddo
    
    if(nCg.ne.nBF) then        
       write(*,'(a,4i4)') 'b2ag_ghostcells: nCg and nBF disagree',nCg,nBF
       stop
    endif

    nCv = nCi + nCg

    !.. calculate nCmxVx, nCmxFc, nFmxCv, nVmxCv, nVmxFc
    nCmxVx = nCmxVx0 + 2*nCg  ! every ghost cell has two vertices
    nCmxFc = nCmxFc0 + nCg    ! every ghost cell has one face
    nFmxCv = nFmxCv0 + nCg    ! every boundary faces is now also used in a ghost cell
    nVmxCv = nVmxCv0 + 2*nCg  ! every vertex at the boundary get two more cells
    nVmxFc = nVmxFc0          ! no new faces created for ghost cells 
    nCmxNv = 9 * nCv          ! over-estimate for 9-point stencil


    return

96  call xerrab('error here')

  end subroutine computeGridSizeWithGhostCells


  !> Rearrange data to insert empty slots for ghost cells
  !> var is of size (-1:nx, -1:ny). 
  !> If withSourceOffset = .false., assumes var(-1:nnx-2, -1:nny-2) is filled with data.
  !> If withSourceOffset = .true.,  assumes var(0:nnx-1, 0:nny-1) is filled with data.
  subroutine rearrange_cells( nx, ny, nnx, nny, niso, nxiso, var, defaultValue, withSourceOffset )
    ! TODO: variable offset
    integer, intent(in) :: nx, ny, nnx, nny, niso, nxiso(:)
    real(R8), intent(inout) :: var(-1:nx,-1:ny)
    real(R8), intent(in) :: defaultValue
    logical, intent(in) :: withSourceOffset

    ! internal
    integer :: iBlock, ixl, iyl, ixu, iyu, xoffset, yoffset, sxoffset, syoffset
    !real(R8) :: tvar(-1:ubound(var,1),-1:ubound(var,2))
    real(R8) :: tvar(-1:nx,-1:ny)

    tvar = defaultValue
    
    ixl = -1
    iyl = -1
    iyu = nny-2

    sxoffset = 0
    syoffset = 0
    if ( withSourceOffset ) then
       sxoffset = 1
       syoffset = 1
!!$       ixl = ixl + 1
!!$       iyl = iyl + 1
!!$       iyu = iyu + 1 
    end if

    yoffset = 1
    xoffset = 1

    ! for every block (i.e. continous block of cells in computational space)
    ! (sn: 1 block, dn: 2 blocks)
    do iBlock = 1, niso+1

        ! figure out block range in x
        if (iBlock == niso + 1) then
           ixu = nnx-2
        else
            ixu = nxiso(iBlock) - 2
        end if
        
        !if ( withSourceOffset ) ixu = ixu + 1

        ! move block 
        tvar( ixl+xoffset:ixu+xoffset, iyl+yoffset:iyu+yoffset  ) = &
             & var( ixl+sxoffset:ixu+sxoffset, iyl+syoffset:iyu+syoffset )

        ! move to next block
        ixl = ixu + 1
        xoffset = xoffset + 2

    end do

    var = tvar
  end subroutine rearrange_cells


  !> Given a B2-style grid consisting out of multiple regions without ghost cells, 
  !> add the ghost cells
  !>
  !> Outline:
  !>  For every cell with a boundary face, reserve space for a ghost cell
  !>  -1st pass: try selecting the slot in the canonical position, i.e. adjacent to the cell
  !>  -2nd pass: for every corner (cell with two adjacent boundary faces), try selecting canonical slot
  !>  -3rd pass: for any faces not resolved up to now, look up alternative empty slot 
  !>                (should only be necessary for extended grid)
  !>  -4th pass: for any corners not resolved up to now, look up alternative empty slot
  !>
  !>  Corners: any cell with more than one boundary face needs "corner ghost cells" to connect
  !>           the ghost cells (which is required by code). 

  subroutine create_guard_cells_st(nnx, nny, nx, ny, niso, nxiso, &
       & crx, cry, psi, bp, bt, ffbz, cflags, withSourceOffset)
    use b2mod_cellhelper

    integer, intent(in) :: nnx, nny, nx, ny, niso, nxiso(:)
    real (kind=R8), intent(inout) :: crx(-1:nx,-1:ny,0:4), cry(-1:nx,-1:ny,0:4), &
         & psi(-1:nx,-1:ny,0:4,0:2), &
         & ffbz(-1:nx,-1:ny,0:3), &
         & bp(-1:nx,-1:ny), bt(-1:nx,-1:ny)
    integer, intent(inout) :: cflags(-1:nx,-1:ny,GRID_N_CELLFLAGS)
    logical, intent(in) :: withSourceOffset

    ! internal 
    ! neighbour connectivity arrays: neighbour on face iFace is at (nbix(ix, iy, iFace), nbiy(ix, iy, iFace))
    integer, dimension(-1:nx, -1:ny, 0:3) :: nbix, nbiy

    integer :: i, j, ix, iy, iPass, iFace, nextFace, gix, giy, nGhostCell, nNoCanonicalGhostCell
    real(R8) :: cflags_tmpreal(-1:nx,-1:ny)   
    logical :: classicalGrid
    
    call logmsg(LOGDEBUG, "create_guard_cells_st: entering")

    ! initialize connectivity
    nbix = NO_CONNECTIVITY
    nbiy = NO_CONNECTIVITY

    do i = 0, 3
        call rearrange_cells( nx, ny, nnx, nny, niso, nxiso, ffbz(:,:,i), &
             & INVALID_POSITION, withSourceOffset )
    end do

    ! rearrange cells to make space for ghost cells according to default B2 system
    do i = 0, 4
        call rearrange_cells( nx, ny, nnx, nny, niso, nxiso, crx(:,:,i), INVALID_POSITION, withSourceOffset )
        call rearrange_cells( nx, ny, nnx, nny, niso, nxiso, cry(:,:,i), INVALID_POSITION, withSourceOffset )
        do j = 0, 2
            call rearrange_cells( nx, ny, nnx, nny, niso, nxiso, psi(:,:,i, j), INVALID_DOUBLE, withSourceOffset )
        end do
    end do
    call rearrange_cells( nx, ny, nnx, nny, niso, nxiso, bp(:,:), INVALID_DOUBLE, withSourceOffset )
    call rearrange_cells( nx, ny, nnx, nny, niso, nxiso, bt(:,:), INVALID_DOUBLE, withSourceOffset )
    do i = 1,ubound(cflags,3)
        cflags_tmpreal(:,:) = real(cflags(:,:,i),R8)
        call rearrange_cells( nx, ny, nnx, nny, niso, nxiso, cflags_tmpreal(:,:), &
            & real(GRID_UNDEFINED, R8), withSourceOffset )
        cflags(:,:,i) = nint(cflags_tmpreal(:,:))
    end do

    ! Dirty hack: are we doing a classical grid, i.e. only covering a part of the vessel domain?
    ! -> Check whether there are any external cells.
    classicalGrid = isClassicalGrid(cflags)
    if (classicalGrid) &
        & call logmsg(LOGDEBUG, "create_guard_cells_st: assuming classical grid")

    ! find slots for ghost cells
    nGhostCell = 0
    nNoCanonicalGhostCell = 0
    do iPass = 1, 3
        do ix = -1, nx
            do iy = -1, ny

                ! only add guard cells for boundary cells
                if ( cflags(ix, iy, CELLFLAG_TYPE) /= GRID_BOUNDARY ) cycle

                do iFace = 0, 3 ! Order see above 
                    if (.not. isBoundary(ix, iy, iFace)) cycle
                    if (isZeroFace(ix, iy, iFace)) cycle                    

                    select case (iPass)
                    case (1)
                       ! ghost on that face: find slot and mark if available
                        nGhostCell = nGhostCell + 1
                        call populateCanonicalGhostSlot(ix, iy, iFace)                  
                    case (2)
                        ! We only do corner ghost cells for classical grids, where we aim to
                        ! exactly replicate the traditional grid setup.
                        if (.not. classicalGrid) cycle
                        ! try to find canonical position for corner ghost cell slot between
                        ! guard cells for faces iFace and iFace + 1
                        nextFace = mod(iFace + 1, 4)
                        if (isBoundary(ix, iy, nextFace)) then
                            nGhostCell = nGhostCell + 1
                            call populateCanonicalCornerGhostSlot(ix, iy, iFace, nextFace)
                        end if
                    case (3)
                        call getCanonicalGhostSlot(ix, iy, iFace, gix, giy)
                        if (.not. (gix == nbix(ix, iy, iFace) .and. giy == nbiy(ix, iy, iFace))) then
                        !if ( .not. isGhostCellLocal(nbix(gix, giy, iFace), nbiy(gix, giy, iFace)) ) then
                            write (*,*) "No canonical ghost cell for face ", ix, iy, iFace
                            nGhostCell = nGhostCell + 1
                            nNoCanonicalGhostCell = nNoCanonicalGhostCell + 1
                            call populateNonCanonicalGhostSlot(ix, iy, iFace)
                        end if
                    end select                    
                end do
                
            end do
        end do        
    end do

    ! psi(dx, dy), ffbz and geometry have been extrapolated to the ghost cells now
    ! However, we cannot compute the magnetic field for the ghost cells now, 
    ! for this we also require the metric information

    call logmsg(LOGDEBUG, "create_guard_cells_st: # of ghost cells "//int2str(nGhostCell) )    
    call logmsg(LOGDEBUG, "create_guard_cells_st: # of non-canonical ghost cells "//int2str(nNoCanonicalGhostCell) )
    call logmsg(LOGDEBUG, "create_guard_cells_st: leaving")

  contains

    ! Check whether a face is a boundary face
    logical function isBoundary(ix, iy, iFace)
      integer, intent(in) :: ix, iy, iFace

      isBoundary = (cflags(ix, iy, CELLFLAG_LEFTFACE + iFace) /= GRID_UNDEFINED)
    end function isBoundary   

    ! Check whether a face is a boundary face
    logical function isZeroFace(ix, iy, iFace)
      integer, intent(in) :: ix, iy, iFace

      integer :: vx1, vx2

      ! first entry of crx, cry is cell midpoint
      vx1 = FACE_VERTEX(VX_START, iFace) + 1
      vx2 = FACE_VERTEX(VX_END, iFace) + 1

      isZeroFace = points_match( &
           & crx(ix, iy, vx1), cry(ix, iy, vx1), &
           & crx(ix, iy, vx2), cry(ix, iy, vx2) )
    end function isZeroFace

    logical function slotAvailable(ix, iy)
      integer, intent(in) :: ix, iy

      slotAvailable = isUnusedCell(cflags(ix, iy, CELLFLAG_TYPE))  
    end function slotAvailable

    ! Set connectivity between two cells: gix, giy is on face iFace of ix, iy
    subroutine connect(ix, iy, iFace, gix, giy)
      integer, intent(in) :: ix, iy, iFace, gix, giy

      nbix(ix, iy, iFace) = gix
      nbiy(ix, iy, iFace) = giy
      nbix(gix, giy, mod(iFace + 2, 4)) = ix
      nbiy(gix, giy, mod(iFace + 2, 4)) = iy
      cflags(gix, giy, CELLFLAG_LEFTFACE + mod(iFace + 2, 4)) &
           & = cflags(ix, iy, CELLFLAG_LEFTFACE + iFace)
    end subroutine connect

    ! Check whether cell ix, iy has neighbour cell on face iFace
    logical function connected(ix, iy, iFace)
      integer, intent(in) :: ix, iy, iFace
      connected = (nbix(ix, iy, iFace) /= NO_CONNECTIVITY) &
             & .and. (nbiy(ix, iy, iFace) /= NO_CONNECTIVITY)
    end function connected
    
    ! Mark a cell slot to be a guard cell and establish connectivity between ghost and inner cell
    subroutine setGhostSlot(ix, iy, iFace, gix, giy, doGeo)
      integer, intent(in) :: ix, iy, iFace, gix, giy
      logical, intent(in) :: doGeo
      
      cflags(gix, giy, CELLFLAG_TYPE) = GRID_GUARD
      call connect(ix, iy, iFace, gix, giy)

      ! extrapolate geometry & other quantities for cell
      if (.not. doGeo) return

      select case (iFace)
      case (LEFT)
          call extendFace( ix, iy, 3, 2, 1, 0, gix, giy )
          call extendFace( ix, iy, 1, 0, 3, 2, gix, giy )
      case (RIGHT)
          call extendFace( ix, iy, 2, 3, 0, 1, gix, giy )
          call extendFace( ix, iy, 0, 1, 2, 3, gix, giy )
      case (BOTTOM)
          call extendFace( ix, iy, 2, 0, 3, 1, gix, giy )
          call extendFace( ix, iy, 3, 1, 2, 0, gix, giy )
      case (TOP)
          call extendFace( ix, iy, 0, 2, 1, 3, gix, giy )
          call extendFace( ix, iy, 1, 3, 0, 2, gix, giy )
      end select

    end subroutine setGhostSlot


    ! Extend the face nodeA1->nodeB1 from cell (oix, oiy)(o=original) to 
    ! to intialize the face nodeA1->nodeB1 in cell (tix, tiy)(t=target).
    ! 
    ! If the face nodeA1->nodeB1 has zero length, use data from face 
    ! nodeA2->nodeB2 instead.

    subroutine extendFace( oix, oiy, nodeA1, nodeB1, nodeA2, nodeB2, tix, tiy)
      integer, intent(in) :: oix, oiy, nodeA1, nodeB1, tix, tiy, nodeA2, nodeB2

      ! internal
      real(R8), parameter :: del = 1.0e-1_R8
      integer :: nodeAO, nodeBO

      ! Copy face endpoint to new face startpoint
      !crx(tix, tiy, nodeA1) = crx(oix, oiy, nodeB1)
      !cry(tix, tiy, nodeA1) = cry(oix, oiy, nodeB1)

      ! Extrapolate new face endpoint
      if (points_match(&
           & crx(ix, iy, nodeA1+1), cry(oix, oiy, nodeA1+1), &
           & crx(oix, oiy, nodeB1+1), cry(oix, oiy, nodeB1+1))) then
          ! Face connecting nodeA1 and nodeB1 in original cell has zero length. Use alternate face 
          ! connecting node nodeA2, nodeB2 for extrapolation.

          ! crx(tix, tiy, nodeB1) = crx(oix, oiy, nodeB1) + ( crx(oix, oiy, nodeB2) - crx(oix, oiy, nodeA2) ) * del
          ! cry(tix, tiy, nodeB1) = cry(oix, oiy, nodeB1) + ( cry(oix, oiy, nodeB2) - cry(oix, oiy, nodeA2) ) * del

          nodeAO = nodeA2
          nodeBO = nodeB2
      else
          ! crx(tix, tiy, nodeB1) = crx(oix, oiy, nodeB1) + ( crx(oix, oiy, nodeB1) - crx(oix, oiy, nodeA1) ) * del
          ! cry(tix, tiy, nodeB1) = cry(oix, oiy, nodeB1) + ( cry(oix, oiy, nodeB1) - cry(oix, oiy, nodeA1) ) * del

          nodeAO = nodeA1
          nodeBO = nodeB1
      end if

      call extendValue( oix, oiy, nodeA1, nodeB1, nodeAO, nodeBO, tix, tiy, crx(:,:,1:4) )
      call extendValue( oix, oiy, nodeA1, nodeB1, nodeAO, nodeBO, tix, tiy, cry(:,:,1:4) )
      call extendValue( oix, oiy, nodeA1, nodeB1, nodeAO, nodeBO, tix, tiy, psi(:,:,1:4,0) )
      call extendValue( oix, oiy, nodeA1, nodeB1, nodeAO, nodeBO, tix, tiy, psi(:,:,1:4,1) )
      call extendValue( oix, oiy, nodeA1, nodeB1, nodeAO, nodeBO, tix, tiy, psi(:,:,1:4,2) )
      call extendValue( oix, oiy, nodeA1, nodeB1, nodeAO, nodeBO, tix, tiy, ffbz)
    end subroutine extendFace        
    
    ! Fill data for face nodeAT->nodeBT in target cell (tix,tiy) by:
    ! -setting data at the start node nodeAT to value of end node nodeBT in original cell (oix,oiy)
    ! -extrapolating data to end node nodeBT in target cell (tix,tiy) by taking the slope from 
    !  the face nodeA0->nodeB0 in the original cell    
    subroutine extendValue( oix, oiy, nodeAT, nodeBT, nodeAO, nodeBO, tix, tiy, val)
      integer, intent(in) :: oix, oiy, nodeAT, nodeBT, nodeAO, nodeBO, tix, tiy
      double precision, intent(inout), dimension(-1:nx,-1:ny,0:3) :: val

      ! internal
      real(R8), parameter :: del = 1.0e-1_R8

      ! Copy value from original (o) face endpoint to target (t) face startpoint
      ! (using corner node indices AT, BT here because we want the actual value)
      val(tix, tiy, nodeAT) = val(oix, oiy, nodeBT)

      ! Extrapolate value to target face endpoint using suitable values from original cells for extrapolation
      ! (using corner node indices AO, BO here because we need a non-zero length face for extrapolation)
      val(tix, tiy, nodeBT) = val(oix, oiy, nodeBT) + ( val(oix, oiy, nodeBO) - val(oix, oiy, nodeAO) ) * del
    end subroutine extendValue


    ! Check whether cell has a ghost cell on the given face
    logical function hasGhostCell(ix, iy, iFace)
      integer, intent(in) :: ix, iy, iFace

      hasGhostCell = (nbix(ix, iy, iFace) /= NO_CONNECTIVITY) &
             & .and. (nbiy(ix, iy, iFace) /= NO_CONNECTIVITY)
      if (.not. hasGhostCell) return
      hasGhostCell = isGhostCellLocal( nbix(ix, iy, iFace), nbiy(ix, iy, iFace) )
    end function hasGhostCell


    ! Check whether cell is a ghost cell. Local wrapper deferring the decision to b2mod_connectivity.isGhostCell.
    logical function isGhostCellLocal(ix, iy)
      integer, intent(in) :: ix, iy
      
      isGhostCellLocal = isInGrid(ix, iy)
      if (.not. isGhostCellLocal) return
      isGhostCellLocal = isGhostCell(cflags(ix, iy, CELLFLAG_TYPE))
    end function isGhostCellLocal


    logical function isInGrid( ix, iy )
      integer, intent(in) :: ix, iy
      isInGrid = (ix >= -1) .and. (ix <= nx) .and. (iy >= -1) .and. (iy <= ny) 
    end function isInGrid

    ! Find the canonical ghost cell slot for a cell on the given face and mark if available
    subroutine populateCanonicalGhostSlot(ix, iy, iFace)
      integer, intent(in) :: ix, iy, iFace

      ! internal
      integer :: gix, giy
      logical :: avail

      call getCanonicalGhostSlot(ix, iy, iFace, gix, giy)
      avail = slotAvailable(gix, giy)
      if (.not. avail) then
          call logmsg(LOGDEBUG, "populateCanonicalGhostSlot: not available for cell "//int2str(ix)//", "//int2str(iy))
          return
      end if
      call setGhostSlot(ix, iy, iFace, gix, giy, .true.)
    end subroutine populateCanonicalGhostSlot
    
    ! Find the canonical ghost cell slot for a cell on the given face and mark if available
    subroutine getCanonicalGhostSlot(ix, iy, iFace, gix, giy)
      integer, intent(in) :: ix, iy, iFace
      integer, intent(out) :: gix, giy

      ! Canonical ghost cell slot: right next to the cell in the given direction
      gix = ix + dxDir(iFace)
      giy = iy + dyDir(iFace)
    end subroutine getCanonicalGhostSlot

    subroutine populateCanonicalCornerGhostSlot(ix, iy, faceA, faceB)
      integer, intent(in) :: ix, iy, faceA, faceB

      ! internal
      integer :: gix, giy

      call getCanonicalCornerGhostSlot(ix, iy, faceA, faceB, gix, giy)
      
      if (.not. slotAvailable(gix, giy)) then
          call logmsg(LOGDEBUG, "populateCanonicalCornerGhostSlot: not available for cell "//int2str(ix)//", "//int2str(iy))
          return
      end if
      
      if (nbix(ix, iy, faceA) /= NO_CONNECTIVITY) then 
          call setGhostSlot(nbix(ix, iy, faceA), nbiy(ix, iy, faceA), faceB, gix, giy, .true.)
      end if

      if (nbix(ix, iy, faceB) /= NO_CONNECTIVITY) then 
          call setGhostSlot(nbix(ix, iy, faceB), nbiy(ix, iy, faceB), faceA, gix, giy, .false. )
      end if

    end subroutine populateCanonicalCornerGhostSlot

    subroutine getCanonicalCornerGhostSlot(ix, iy, faceA, faceB, gix, giy)
      integer, intent(in) :: ix, iy, faceA, faceB
      integer, intent(out) :: gix, giy

      if (faceA == LEFT .and. faceB == BOTTOM) then
          gix = ix - 1
          giy = iy - 1
      else if (faceA == BOTTOM .and. faceB == RIGHT) then
          gix = ix + 1
          giy = iy - 1
      else if (faceA == RIGHT .and. faceB == TOP) then
          gix = ix + 1
          giy = iy + 1
      else if (faceA == TOP .and. faceB == LEFT) then
          gix = ix - 1
          giy = iy + 1
      else
          stop 'getCanonicalCornerGhostSlot: strange face combination'
      end if
    end subroutine getCanonicalCornerGhostSlot


    ! Find the canonical ghost cell slot for a cell on the given face and mark if available
    subroutine populateNoncanonicalGhostSlot(ix, iy, iFace)
      integer, intent(in) :: ix, iy, iFace

      ! internal
      integer :: gix, giy

      call getNoncanonicalGhostSlot(ix, iy, iFace, gix, giy)
      if (gix == NO_CONNECTIVITY) stop "populateNoncanonicalGhostSlot: no free slot left"
      call logmsg(LOGDEBUG, "populateNoncanonicalGhostSlot: face "&
          & //int2str(ix)//", "//int2str(iy)//", "//int2str(iFace)//" using slot "//int2str(gix)//", "//int2str(giy))
      call setGhostSlot(ix, iy, iFace, gix, giy, .true.)
    end subroutine populateNoncanonicalGhostSlot
    
    ! Find a free non-canonical ghost slot. If no slot can be found, 
    ! returns NO_CONNECTIVITY for both gix, giy
    subroutine getNoncanonicalGhostSlot(ix, iy, iFace, gix, giy)
      integer, intent(in) :: ix, iy, iFace
      integer, intent(out) :: gix, giy

      ! internal
      integer :: xstep, ystep, dx, dy
      logical :: cellFound

      select case (iFace)
      case(LEFT, RIGHT)

          if (iFace == LEFT) dx = -1
          if (iFace == RIGHT) dx = +1

          gix = ix
          giy = iy
          do ystep =  1, ny + 2
              do xstep = 1, nx + 2
                  gix = gix + dx
                  if (gix < -1) gix = nx
                  if (gix > nx) gix = -1
                  cellFound = slotAvailable(gix, giy)
                  if (cellFound) exit
              end do
              if (cellFound) exit
              giy = giy - 1
              if (giy < -1) giy = ny
          end do

      case(BOTTOM, TOP)

          if (iFace == BOTTOM) dy = -1
          if (iFace == TOP) dy = +1

          gix = ix
          giy = iy
          do xstep = 1, nx + 2
              do ystep =  1, ny + 2
                  giy = giy + dy
                  if (giy < -1) giy = ny
                  if (giy > ny) giy = -1
                  cellFound = slotAvailable(gix, giy)
                  if (cellFound) exit
              end do
              if (cellFound) exit
              gix = gix - 1
              if (gix < -1) gix = nx
          end do

      end select

      if (.not. cellFound) then
          gix = NO_CONNECTIVITY
          giy = NO_CONNECTIVITY
      end if

    end subroutine getNoncanonicalGhostSlot
    
  end subroutine create_guard_cells_st

  subroutine create_guard_cells(nnCv,nnFc,nnVx,g,m)
    use b2us_geo
    use b2us_map

    integer, intent(in) :: nnCv, nnFc, nnVx
    type(geometry), intent(inout) :: g
    type(mapping), intent(inout) :: m

    !internal

    integer :: i, iCv, cellnumber, cvFcpos, cvVxpos, nGhostCell
    !integer :: facenumbers(nGhostCell)
    integer :: facenumbers(m%nCg)  ! too big
    real(kind=R8) :: fcX,fcY 


    call logmsg(LOGDEBUG, "create_guard_cells:entering")

    ! initialize connectivity
    

    ! rearrange cells to make space for ghost cells according to default B2 system - not necessary
    ! cvFc and cvFcP changes here, as ghost cells are only one face cells
    ! cvVx and cvVxP changes here
    ! in cflags(:,1) guard_cell gets flag == GRID_GUARD
    nGhostCell = m%nCv - nnCv
    ! change cvFcP and cvVxP
    do iCv = nnCv+1, m%nCv
       m%cvFcP(iCv,2) = 1
       m%cvFcP(iCv,1) = m%cvFcP(iCv-1,1) + m%cvFcP(iCv-1,2)

       m%cvVxP(iCv,2) = 2
       m%cvVxP(iCv,1) = m%cvVxP(iCv-1,1) + m%cvVxP(iCv-1,2)
    enddo
    
    ! change cvFc
    ! create list of boundary faces
    call get_boundary_faces(facenumbers,g,m,nGhostCell,nnCv)

    do i = 1, nGhostCell
       cellnumber = nnCv + i 
       cvFcpos = m%cvFcP(cellnumber,1)
       m%cvFc(cvFcpos) = facenumbers(i) 

       cvVxpos = m%cvVxP(cellnumber,1)
       m%cvVx(cvVxpos) = m%fcVx(facenumbers(i),1)
       m%cvVx(cvVxpos+1) = m%fcVx(facenumbers(i),2)
    enddo

    ! change cvLbl(:,1), cvX, cvY
    do iCv = nnCv+1, m%nCv
       m%cvLbl(iCv) = GRID_GUARD
       !iface = m%cvFc(m%cvFc(i,1)) ! just one face
       ! temporary put the cell center of the guard cells of the boundary face
       ! in b2agmt the proper cell center is computed
       fcX = 0.5_R8*sum(g%vxX(m%cvVx(m%cvVxP(iCv,1): &
           &    m%cvVxP(iCv,1)+m%cvVxP(iCv,2)-1)))
       g%cvX(iCv) = fcX
       fcY = 0.5_R8*sum(g%vxY(m%cvVx(m%cvVxP(iCv,1): &
           &    m%cvVxP(iCv,1)+m%cvVxP(iCv,2)-1)))
       g%cvY(iCv) = fcY     
    enddo

    ! add flux tubes of guard cells
    do iCv = nnCv+1,m%nCv
       
    enddo
    

    ! calculate cvFpsi , cvBt and cvBp in routine b2agbb in b2agfs.f

    call logmsg(LOGDEBUG, "create_guard_cells: # of ghost cells "//int2str(nGhostCell))
    call logmsg(LOGDEBUG, "create_guard_cells: leaving")

  contains
    
    ! Give a list of the boundary faces
    subroutine get_boundary_faces(facenumbers,g,m,nGhostCell,nnCv)
      use b2us_geo
      use b2us_map

      integer, intent(in) :: nGhostCell, nnCv
      integer, intent(out) :: facenumbers(nGhostCell)
      type(mapping) :: m 
      type(geometry) :: g

      !internal
      integer :: iCv, j, counter, facesQuad(0:3), facesTria(0:2)
      logical :: Bface
      intrinsic count
      
      counter = 1
      do iCv = 1, nnCv
          if (m%cvLbl(iCv).eq.GRID_BOUNDARY) then
             if (m%cvFcP(iCv,2).eq.4) then  !Quad
               facesQuad = m%cvFc(m%cvFcP(iCv,1): &
                 &   m%cvFcP(iCv,1)+m%cvFcP(iCv,2)-1)
               do j = 0, 3
                  Bface = isBoundaryFace(facesQuad(j))
                  if (Bface) then
                     facenumbers(counter) = facesQuad(j)
                     counter = counter+1
                  endif
               enddo
             else !Triangle
               facesTria = m%cvFc(m%cvFcP(iCv,1): &
                 &   m%cvFcP(iCv,1)+m%cvFcP(iCv,2)-1)
               do j = 0,2
                 Bface = isBoundaryFace(facesTria(j))
                 if (Bface) then
                    facenumbers(counter) = facesTria(j)
                    counter = counter+1
                 endif
               enddo
             endif

          endif
      enddo

    end subroutine get_boundary_faces

    logical function isBoundaryFace(iface)
      integer, intent(in) :: iface

      integer :: counter

      counter = 0
      ! count how many times iface is present in cvFc
      !do i = 1, m%nCmxFc
      !   if ((m%cvFc(i)-1).eq.iface) then
      !       counter = counter +1
      !   endif 
      !enddo
     counter = count (m%cvFc == iface)


      ! if the face is present once, it is a boundary face
      if (counter.eq.1) then
         isBoundaryFace = 1
      else
         isBoundaryFace = 0
      endif 

    end function isBoundaryFace


  end subroutine create_guard_cells



  subroutine read_cvFc(nrid, nCv, cvFc,nCmxFc)
     
      use b2mod_types
      implicit none
      integer, intent(in) :: nrid, nCv, nCmxFc
      integer :: cflags1(nCv), &
       & cv(nCv), cvVxP1(nCv), & 
       & cvVxP2(nCv), cvreg(nCv), cvFc(nCmxFc)
      real (kind=R8) :: &
       & cvX(nCv), cvY(nCv), &
       & psi(nCv), bp(nCv), bt(nCv)
      integer, allocatable :: cvVx(:) 
      
      integer i, j
      character idcod*5, id0*32, id1*32
      logical streql
      external xertst, streql

      !..subprogram start-up calls
      call subini('read_cflags')

      j=0
      idcod = 'start'
      id1 = idcod

99    if  (j.eq.0) then
        if ((streql('*cf:',idcod)).and. &
       &   (streql('cv' ,id1))) then
           do i = 1, nCv
             read(nrid,*,err=97) cv(i), cvVxP1(i), cvVxP2(i), cvX(i), &
       &     cvY(i), psi(i), bp(i), bt(i), cflags1(i), cvreg(i)
           enddo
           j=1
           allocate(cvVx(cvVxP1(nCv)+cvVxP2(nCv)-1))
           !read cvVx
           call cfruin(nrid,cvVxP1(nCv)+cvVxP2(nCv)-1,cvVx,'cvVx')
           !read cvFc
           call cfruin(nrid,nCmxFc,cvFc,'cvFc')           
        else
         read(nrid,*,err=97) idcod, id1
!         read(nrid,'(2a8,a32)',err=96) idcod, id1 !read the next line
         goto 99
        end if
      end if 

      return

96 call xerrab('read_cflags - error at head')
97 call xerrab('read_cflags - error reading grid file')           


  end subroutine read_cvFc


end module b2ag_ghostcells
