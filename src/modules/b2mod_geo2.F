      module b2mod_geo2
      use b2mod_types
      use b2mod_indirect
     . , only: GRID_N_CELLFLAGS
      use b2mod_cellhelper
      use carre_constants
      use b2mod_grid_mapping

      real (kind=R8), allocatable, save ::
     *  bb1(:,:,:),
     *  crx1(:,:,:),
     *  cry1(:,:,:),
     *  ffbz1(:,:,:),
     *  fpsi1(:,:,:,:),
     *  gs1(:,:,:),
     *  hxx(:,:),
     *  hyy(:,:),
     *  hz1(:,:),
     *  hc1(:,:,:),
     *  ht1(:,:,:),
     *  qz1(:,:,:),
     *  qc1(:,:,:),
     *  qs1(:,:,:),
     *  qac1(:,:,:),
     *  qas1(:,:,:),
     *  pbs1(:,:,:),
     *  ebc1(:,:,:),
     *  wbbl1(:,:,:),
     *  wbbr1(:,:,:),
     *  wbbv1(:,:,:),
     *  wbbc1(:,:,:),
     *  vol1(:,:),
     *  cell_width1(:,:),
     *  cell_height1(:,:)

      integer, allocatable, save ::
     *  leftix1(:,:),
     *  leftiy1(:,:),
     *  rightix1(:,:),
     *  rightiy1(:,:),
     *  topix1(:,:),
     *  topiy1(:,:),
     *  bottomix1(:,:),
     *  bottomiy1(:,:),
     *  leftcut1(:),
     *  rightcut1(:),
     *  bottomcut1(:),
     *  topcut1(:),
     *  region1(:,:,:),
     *  resignore1(:,:,:),
     *  cflags1(:,:,:)

      type(B2GridMap), save :: gmap1

      logical, save :: b2mod_geo2_allocated = .false.

      contains

      subroutine alloc_b2mod_geo2(nxd,nyd,nncutmax)
      implicit none
      integer nxd, nyd, nncutmax
      allocate(bb1(-1:nxd,-1:nyd,0:3))
      allocate(crx1(-1:nxd,-1:nyd,0:3))
      allocate(cry1(-1:nxd,-1:nyd,0:3))
      allocate(ffbz1(-1:nxd,-1:nyd,0:3))
      allocate(fpsi1(-1:nxd,-1:nyd,0:3,0:2))
      allocate(gs1(-1:nxd,-1:nyd,0:2))
      allocate(hxx(-1:nxd,-1:nyd))
      allocate(hyy(-1:nxd,-1:nyd))
      allocate(hz1(-1:nxd,-1:nyd))
      allocate(hc1(-1:nxd,-1:nyd,0:3))
      allocate(ht1(-1:nxd,-1:nyd,0:1))
      allocate(qz1(-1:nxd,-1:nyd,0:1))
      allocate(qc1(-1:nxd,-1:nyd,0:1))
      allocate(qs1(-1:nxd,-1:nyd,0:1))
      allocate(qac1(-1:nxd,-1:nyd,0:1))
      allocate(qas1(-1:nxd,-1:nyd,0:1))
      allocate(pbs1(-1:nxd,-1:nyd,0:1))
      allocate(ebc1(-1:nxd,-1:nyd,0:2))
      allocate(wbbl1(-1:nxd,-1:nyd,0:3))
      allocate(wbbr1(-1:nxd,-1:nyd,0:3))
      allocate(wbbv1(-1:nxd,-1:nyd,0:3))
      allocate(wbbc1(-1:nxd+1,-1:nyd+1,0:3))
      allocate(vol1(-1:nxd,-1:nyd))
      allocate(cell_width1(-1:nxd,-1:nyd))
      allocate(cell_height1(-1:nxd,-1:nyd))
      allocate(leftix1(-1:nxd,-1:nyd))
      allocate(leftiy1(-1:nxd,-1:nyd))
      allocate(rightix1(-1:nxd,-1:nyd))
      allocate(rightiy1(-1:nxd,-1:nyd))
      allocate(topix1(-1:nxd,-1:nyd))
      allocate(topiy1(-1:nxd,-1:nyd))
      allocate(bottomix1(-1:nxd,-1:nyd))
      allocate(bottomiy1(-1:nxd,-1:nyd))
      allocate(leftcut1(1:nncutmax))
      allocate(rightcut1(1:nncutmax))
      allocate(topcut1(1:nncutmax))
      allocate(bottomcut1(1:nncutmax))
      allocate(region1(-1:nxd,-1:nyd,0:2))
      allocate(resignore1(-1:nxd,-1:nyd,1:2))
      allocate(cflags1(-1:nxd,-1:nyd,GRID_N_CELLFLAGS))
      b2mod_geo2_allocated = .true.
      end subroutine alloc_b2mod_geo2

      subroutine dealloc_b2mod_geo2
      implicit none

      if(b2mod_geo2_allocated) then
         deallocate(bb1)
         deallocate(crx1)
         deallocate(cry1)
         deallocate(ffbz1)
         deallocate(fpsi1)
         deallocate(wbbc1)
         deallocate(gs1)
         deallocate(hxx)
         deallocate(hyy)
         deallocate(hz1)
         deallocate(hc1)
         deallocate(ht1)
         deallocate(qz1)
         deallocate(qc1)
         deallocate(qs1)
         deallocate(qac1)
         deallocate(qas1)
         deallocate(pbs1)
         deallocate(ebc1)
         deallocate(wbbl1)
         deallocate(wbbr1)
         deallocate(wbbv1)
         deallocate(vol1)
         deallocate(cell_width1)
         deallocate(cell_height1)
         deallocate(leftix1)
         deallocate(leftiy1)
         deallocate(rightix1)
         deallocate(rightiy1)
         deallocate(topix1)
         deallocate(topiy1)
         deallocate(bottomix1)
         deallocate(bottomiy1)
         deallocate(leftcut1)
         deallocate(rightcut1)
         deallocate(topcut1)
         deallocate(bottomcut1)
         deallocate(region1)
         deallocate(resignore1)
         deallocate(cflags1)
         b2mod_geo2_allocated = .false.
      endif
      if (map1Initialized) call deallocateB2GridMap(gmap1)
      end subroutine dealloc_b2mod_geo2

  !> Computes the standard B2 connectivity information
  !! @note  A note on the cut arrays:
  !!        leftcut(i) holds the left boundary index of a region which is cut
  !!        rightcut(i) holds the right boundary index + 1 of a region which is cut
  !!        This means the region range in the x direction is
  !!        (leftcut(i):rightcut(i))
      subroutine init_connectivity( nx, ny, cflag,
     &     periodic_bc, nncut, inseltop, inselbot, istyle )

      use b2mod_types
      implicit none
#include <DIMENSIONS.F>

    !!   ..input arguments (unchanged on exit)
      integer, intent(in) :: nx, ny, istyle, periodic_bc
      integer cflag(-1:nx,-1:ny,CARREOUT_NCELLFLAGS)
    !!   .. output arguments
      integer, intent(out) ::
     & nncut, inseltop, inselbot

    !! internal

      integer :: ic, ix, iy, i, ixNb, iyNb, ixl, ixr
      integer :: xstep, ystep, nNb
      logical :: cutFound, cellFound, fullGrid

      nncut=0
      bottomcut1(:)=ny+1
      topcut1(:)=-2
      rightcut1(:)=-2
      leftcut1(:)=nx+1
      ic=0 !! cut counter

      rightix1 = NO_CONNECTIVITY
      rightiy1 = NO_CONNECTIVITY
      leftix1 = NO_CONNECTIVITY
      leftiy1 = NO_CONNECTIVITY
      topix1 = NO_CONNECTIVITY
      topiy1 = NO_CONNECTIVITY
      bottomix1 = NO_CONNECTIVITY
      bottomiy1 = NO_CONNECTIVITY
      inseltop = NO_CONNECTIVITY
      inselbot = NO_CONNECTIVITY

    !! First step: find cell connectivity

      do iy=-1,ny
        do ix=-1,nx

            !! unused cells have no connectivity, skip
            if (isUnusedCell(cflag(ix,iy,CELLFLAG_TYPE))) cycle

            !! Look for the left/right connectivity
            if (distc1(ix,iy,0,ix,iy,2) > 0.0_R8) then
              ixNb = ix
              iyNb = iy
              cellFound = .false.

              do ystep =  1, ny + 2
                do xstep = 1, nx + 2
                  ixNb = ixNb - 1
                  if (ixNb < -1) ixNb = nx

                  if (isUnusedCell(cflag(ixNb,iyNb,CELLFLAG_TYPE)))
     .             cycle
                  cellFound = matchLeft1(ixNb,iyNb,ix,iy)
                  if (cellFound) exit
                end do
                if (cellFound) exit
                iyNb = iyNb - 1
                if (iyNb < -1) iyNb = ny
              end do

              if (cellFound) then
                !! set connectivity in both directions
                leftix1(ix,iy)=ixNb
                leftiy1(ix,iy)=iyNb
                rightix1( ixNb, iyNb ) = ix
                rightiy1( ixNb, iyNb ) = iy
              end if
            end if

            !! Look for the top/bottom connectivity
            if (distc1(ix,iy,0,ix,iy,1) > 0.0_R8) then
              ixNb = ix
              iyNb = iy
              cellFound = .false.

              do xstep = 1, nx + 2
                do ystep =  1, ny + 2
                  iyNb = iyNb - 1
                  if (iyNb < -1) iyNb = ny

                  if (isUnusedCell(cflag(ixNb,iyNb,CELLFLAG_TYPE)))
     .             cycle
                  cellFound = matchBottom1(ixNb,iyNb,ix,iy)
                  if (cellFound) exit
                end do
                if (cellFound) exit
                ixNb = ixNb - 1
                if (ixNb < -1) ixNb = nx
              end do

              if (cellFound) then
                !! set connectivity in both directions
                bottomix1(ix,iy)=ixNb
                bottomiy1(ix,iy)=iyNb
                topix1(ixNb,iyNb)=ix
                topiy1(ixNb,iyNb)=iy
              else
                !write (*,*) "No bottom neighbour for ", ix, iy
              end if
            end if

        end do
      end do

    !! Fix connectivity to match B2 convention
      do iy=-1,ny
        do ix=-1,nx
          if (leftix1(ix, iy) == NO_CONNECTIVITY) then
            leftix1(ix, iy) = -2
            leftiy1(ix, iy) = iy
          end if
          if (rightix1(ix, iy) == NO_CONNECTIVITY) then
            rightix1(ix, iy) = nx + 1
            rightiy1(ix, iy) = iy
          end if
          if (bottomix1(ix, iy) == NO_CONNECTIVITY) then
            bottomix1(ix, iy) = ix
            bottomiy1(ix, iy) = -2
          end if
          if (topix1(ix, iy) == NO_CONNECTIVITY) then
            topix1(ix, iy) = ix
            topiy1(ix, iy) = ny + 1
          end if
        end do
      end do

    !! second step (only for "classical" grids with no cell type information):
    !! identify ghost cells

      if ( count(isGhostCell(cflag(:,:,CELLFLAG_TYPE))) == 0 ) then

       !! find guard cells
        do iy=-1,ny
          do ix=-1,nx

            nNb = 0
            if (isInDomain(nx, ny,
     .           leftix1(ix, iy), leftiy1(ix, iy))) nNb = nNb + 1
            if (isInDomain(nx, ny,
     .           rightix1(ix, iy), rightiy1(ix, iy))) nNb = nNb + 1
            if (isInDomain(nx, ny,
     .           topix1(ix, iy), topiy1(ix, iy))) nNb = nNb + 1
            if (isInDomain(nx, ny,
     .           bottomix1(ix, iy), bottomiy1(ix, iy))) nNb = nNb + 1

            if (nNb < 4) cflag(ix, iy, CELLFLAG_TYPE) = GRID_GUARD

          end do
        end do

       !! mark boundary cells
        do iy=-1,ny
          do ix=-1,nx

            if (.not. isGhostCell(cflag(ix, iy, CELLFLAG_TYPE))) cycle

            if (isInDomain(nx, ny,
     .           leftix1(ix, iy), leftiy1(ix, iy))) then
              if ( isRealCell(cflag(leftix1(ix, iy),
     .                        leftiy1(ix, iy), CELLFLAG_TYPE)) ) then
                cflag(leftix1(ix, iy),
     .                leftiy1(ix, iy), CELLFLAG_TYPE) = GRID_BOUNDARY
                cflag(ix, iy, CELLFLAG_LEFTFACE) = BOUNDARY_NOSTRUCTURE
                cflag(leftix1(ix, iy),
     .                leftiy1(ix, iy), CELLFLAG_RIGHTFACE) =
     .                                             BOUNDARY_NOSTRUCTURE
              end if
            end if

            if (isInDomain(nx, ny,
     .           rightix1(ix, iy), rightiy1(ix, iy))) then
              if ( isRealCell(cflag(rightix1(ix, iy),
     .                         rightiy1(ix, iy), CELLFLAG_TYPE)) ) then
                cflag(rightix1(ix, iy),
     .                rightiy1(ix, iy), CELLFLAG_TYPE) = GRID_BOUNDARY
                cflag(ix, iy, CELLFLAG_RIGHTFACE) = BOUNDARY_NOSTRUCTURE
                cflag(rightix1(ix, iy),
     .                rightiy1(ix, iy), CELLFLAG_LEFTFACE) =
     .                                              BOUNDARY_NOSTRUCTURE
              end if
            end if

            if (isInDomain(nx, ny,
     .           topix1(ix, iy), topiy1(ix, iy))) then
              if ( isRealCell(cflag(topix1(ix, iy),
     .                         topiy1(ix, iy), CELLFLAG_TYPE)) ) then
                cflag(topix1(ix, iy),
     .                topiy1(ix, iy), CELLFLAG_TYPE) = GRID_BOUNDARY
                cflag(ix, iy, CELLFLAG_TOPFACE) = BOUNDARY_NOSTRUCTURE
                cflag(topix1(ix, iy),
     .                topiy1(ix, iy), CELLFLAG_BOTTOMFACE) =
     .                                            BOUNDARY_NOSTRUCTURE
              end if
            end if

            if (isInDomain(nx, ny,
     .           bottomix1(ix, iy), bottomiy1(ix, iy))) then
              if ( isRealCell(cflag(bottomix1(ix, iy),
     .                         bottomiy1(ix, iy), CELLFLAG_TYPE)) ) then
                cflag(bottomix1(ix, iy),
     .                bottomiy1(ix, iy), CELLFLAG_TYPE) = GRID_BOUNDARY
                cflag(ix, iy, CELLFLAG_BOTTOMFACE) =
     .                                             BOUNDARY_NOSTRUCTURE
                cflag(bottomix1(ix, iy),
     .                bottomiy1(ix, iy), CELLFLAG_TOPFACE) =
     .                                             BOUNDARY_NOSTRUCTURE
              end if
            end if

          end do
        end do

      else

    !! second step (for "extended" grids with no cell face information):
    !! make sure ghost cells are not connected across different walls
    !! for full grids, we wish to allow corner cells to be connected

        fullGrid =
     .   (count(isUnusedCell(cflag(0:nx-1,0:ny-1,CELLFLAG_TYPE)))
     .     == 0)

        do ix = -1, nx
         do iy = -1, ny
          if(.not. isGhostCell(cflag(ix, iy, CELLFLAG_TYPE))) cycle
          if(cflag(ix,iy,CELLFLAG_LEFTFACE) /= GRID_UNDEFINED) then
            if(isInDomain(nx,ny,topix1(ix,iy),topiy1(ix,iy))) then
              if(cflag(ix,iy,CELLFLAG_LEFTFACE).ne.
     .           cflag(topix1(ix,iy),topiy1(ix,iy),CELLFLAG_LEFTFACE)
     .         .and. (((cflag(ix,iy,CELLFLAG_LEFTFACE).gt.0 .or.
     .                  cflag(topix1(ix,iy),topiy1(ix,iy),
     .                              CELLFLAG_LEFTFACE).gt.0) .and.
     .         .not.fullGrid) .or.
     .          (cflag(topix1(ix,iy),topiy1(ix,iy),CELLFLAG_LEFTFACE)
     .            /= GRID_UNDEFINED .and. fullGrid))) then
                bottomix1(topix1(ix,iy),topiy1(ix,iy)) = ix
                bottomiy1(topix1(ix,iy),topiy1(ix,iy)) = -2
                topix1(ix,iy) = ix
                topiy1(ix,iy) = ny+1
              end if
            end if
            if(isInDomain(nx,ny,bottomix1(ix,iy),
     .                          bottomiy1(ix,iy))) then
              if(cflag(ix,iy,CELLFLAG_LEFTFACE).ne.
     .           cflag(bottomix1(ix,iy),bottomiy1(ix,iy),
     .                       CELLFLAG_LEFTFACE) .and.
     .        (((cflag(ix,iy,CELLFLAG_LEFTFACE).gt.0 .or.
     .           cflag(bottomix1(ix,iy),bottomiy1(ix,iy),
     .                       CELLFLAG_LEFTFACE).gt.0) .and.
     .         .not.fullGrid) .or.
     .          (cflag(bottomix1(ix,iy),bottomiy1(ix,iy),
     .                       CELLFLAG_LEFTFACE) /= GRID_UNDEFINED
     .         .and. fullGrid))) then
                topix1(bottomix1(ix,iy),bottomiy1(ix,iy)) = ix
                topiy1(bottomix1(ix,iy),bottomiy1(ix,iy)) = ny+1
                bottomix1(ix,iy) = ix
                bottomiy1(ix,iy) = -2
              end if
            end if
          end if
          if(cflag(ix,iy,CELLFLAG_BOTTOMFACE) /= GRID_UNDEFINED) then
            if(isInDomain(nx,ny,leftix1(ix,iy),leftiy1(ix,iy))) then
              if(cflag(ix,iy,CELLFLAG_BOTTOMFACE).ne.
     .           cflag(leftix1(ix,iy),
     .                 leftiy1(ix,iy),CELLFLAG_BOTTOMFACE) .and.
     .        (((cflag(ix,iy,CELLFLAG_BOTTOMFACE).gt.0 .or.
     .           cflag(leftix1(ix,iy),leftiy1(ix,iy),
     .                       CELLFLAG_BOTTOMFACE).gt.0) .and.
     .          .not.fullGrid) .or.
     .          (cflag(leftix1(ix,iy),leftiy1(ix,iy),
     .                       CELLFLAG_BOTTOMFACE) /= GRID_UNDEFINED
     .          .and. fullGrid))) then
                rightix1(leftix1(ix,iy),leftiy1(ix,iy)) = nx+1
                rightiy1(leftix1(ix,iy),leftiy1(ix,iy)) = iy
                leftix1(ix,iy) = -2
                leftiy1(ix,iy) = iy
              end if
            end if
            if(isInDomain(nx,ny,rightix1(ix,iy),rightiy1(ix,iy))) then
              if(cflag(ix,iy,CELLFLAG_BOTTOMFACE).ne.
     .           cflag(rightix1(ix,iy),rightiy1(ix,iy),
     .                       CELLFLAG_BOTTOMFACE) .and.
     .        (((cflag(ix,iy,CELLFLAG_BOTTOMFACE).gt.0 .or.
     .           cflag(rightix1(ix,iy),rightiy1(ix,iy),
     .                       CELLFLAG_BOTTOMFACE).gt.0) .and.
     .          .not.fullGrid) .or.
     .          (cflag(rightix1(ix,iy),rightiy1(ix,iy),
     .                       CELLFLAG_BOTTOMFACE) /= GRID_UNDEFINED
     .          .and. fullGrid))) then
                leftix1(rightix1(ix,iy),rightiy1(ix,iy)) = -2
                leftiy1(rightix1(ix,iy),rightiy1(ix,iy)) = iy
                rightix1(ix,iy) = nx+1
                rightiy1(ix,iy) = -2
              end if
            end if
          end if
          if(cflag(ix,iy,CELLFLAG_RIGHTFACE) /= GRID_UNDEFINED) then
            if(isInDomain(nx,ny,topix1(ix,iy),topiy1(ix,iy))) then
              if(cflag(ix,iy,CELLFLAG_RIGHTFACE).ne.
     .           cflag(topix1(ix,iy),topiy1(ix,iy),
     .                       CELLFLAG_RIGHTFACE) .and.
     .        (((cflag(ix,iy,CELLFLAG_RIGHTFACE).gt.0 .or.
     .           cflag(topix1(ix,iy),topiy1(ix,iy),
     .                       CELLFLAG_RIGHTFACE).gt.0) .and.
     .           .not.fullGrid) .or.
     .          (cflag(topix1(ix,iy),topiy1(ix,iy),
     .                       CELLFLAG_RIGHTFACE) /= GRID_UNDEFINED
     .           .and. fullGrid))) then
                bottomix1(topix1(ix,iy),topiy1(ix,iy)) = ix
                bottomiy1(topix1(ix,iy),topiy1(ix,iy)) = -2
                topix1(ix,iy) = ix
                topiy1(ix,iy) = ny+1
              end if
            end if
            if(isInDomain(nx,ny,bottomix1(ix,iy),
     .                          bottomiy1(ix,iy))) then
              if(cflag(ix,iy,CELLFLAG_RIGHTFACE).ne.
     .           cflag(bottomix1(ix,iy),bottomiy1(ix,iy),
     .                       CELLFLAG_RIGHTFACE) .and.
     .        (((cflag(ix,iy,CELLFLAG_RIGHTFACE).gt.0 .or.
     .           cflag(bottomix1(ix,iy),bottomiy1(ix,iy),
     .                       CELLFLAG_RIGHTFACE).gt.0) .and.
     .           .not.fullGrid) .or.
     .          (cflag(bottomix1(ix,iy),bottomiy1(ix,iy),
     .                       CELLFLAG_RIGHTFACE) /= GRID_UNDEFINED
     .           .and. fullGrid))) then
                topix1(bottomix1(ix,iy),bottomiy1(ix,iy)) = ix
                topiy1(bottomix1(ix,iy),bottomiy1(ix,iy)) = ny+1
                bottomix1(ix,iy) = ix
                bottomiy1(ix,iy) = -2
              end if
            end if
          end if
          if(cflag(ix,iy,CELLFLAG_TOPFACE) /= GRID_UNDEFINED) then
            if(isInDomain(nx,ny,leftix1(ix,iy),leftiy1(ix,iy))) then
              if(cflag(ix,iy,CELLFLAG_TOPFACE).ne.
     .           cflag(leftix1(ix,iy),leftiy1(ix,iy),
     .                       CELLFLAG_TOPFACE) .and.
     .        (((cflag(ix,iy,CELLFLAG_TOPFACE).gt.0 .or.
     .           cflag(leftix1(ix,iy),leftiy1(ix,iy),
     .                       CELLFLAG_TOPFACE).gt.0) .and.
     .           .not.fullGrid) .or.
     .          (cflag(leftix1(ix,iy),leftiy1(ix,iy),
     .                       CELLFLAG_TOPFACE) /= GRID_UNDEFINED
     .           .and. fullGrid))) then
                rightix1(leftix1(ix,iy),leftiy1(ix,iy)) = nx+1
                rightiy1(leftix1(ix,iy),leftiy1(ix,iy)) = iy
                leftix1(ix,iy) = -2
                leftiy1(ix,iy) = iy
              end if
            end if
            if(isInDomain(nx,ny,rightix1(ix,iy),
     .                          rightiy1(ix,iy))) then
              if(cflag(ix,iy,CELLFLAG_TOPFACE).ne.
     .           cflag(rightix1(ix,iy),rightiy1(ix,iy),
     .                       CELLFLAG_TOPFACE) .and.
     .        (((cflag(ix,iy,CELLFLAG_TOPFACE).gt.0 .or.
     .           cflag(rightix1(ix,iy),rightiy1(ix,iy),
     .                       CELLFLAG_TOPFACE).gt.0) .and.
     .           .not.fullGrid) .or.
     .          (cflag(rightix1(ix,iy),rightiy1(ix,iy),
     .                       CELLFLAG_TOPFACE) /= GRID_UNDEFINED
     .           .and. fullGrid))) then
                leftix1(rightix1(ix,iy),rightiy1(ix,iy)) = -2
                leftiy1(rightix1(ix,iy),rightiy1(ix,iy)) = iy
                rightix1(ix,iy) = nx+1
                rightiy1(ix,iy) = -2
              end if
            end if
          end if
         end do
        end do

      end if


    !! third step: find cuts

      do iy=-1,ny
        do ix=-1,nx

          !! unused cells have no connectivity, skip
          if (isUnusedCell(cflag(ix,iy,CELLFLAG_TYPE))) cycle

          !! Get left neighbour
          ixNb = leftix1(ix,iy)
          iyNb = leftiy1(ix,iy)

          if (.not. isInDomain(nx, ny, ixNb, iyNb)) cycle

          if ((periodic_bc == 1) .and.
     &        (ixNb == nx-1) .and. (ix == 0)) cycle  !! Limiter connectivity case

          !! if neither this nor neighbour cell is a ghost cell
          !! and cells on same horizontal line but not next to each other
          !! do bookkeeping for cut
          if ( ( (iy /= iyNb) .or. (ix /= ixNb +1) ) .and.
     &           .not.isGhostCell(cflag(ix,iy,CELLFLAG_TYPE)) .and.
     &           .not.isGhostCell(cflag(ixNb,iyNb,CELLFLAG_TYPE)) ) then

             !! Found a left neighbour inside the domain on same line but not directly
             !! left of the cell -> there is a cut at the left face of the cell
             !! Figure out right and left end of cut region
             if (ixNb<ix) then
                !! jump to the left
                ixr = ix
                ixl = ixNb+1
             else
                !! jump to the right
                ixr = ixNb+1
                ixl = ix
             end if

             cutFound=.false.
             do i = 1, ic
                if(ixr.eq.rightcut1(i) .or.
     .             ixl.eq.leftcut1(i)) then
                   cutFound = .true.
                   exit
                endif
             end do
             if (.not.cutFound) then
                ic = ic + 1
                write(*,*) 'ic',ic
                i = ic
             endif
             bottomcut1(i)=min(iy,bottomcut1(i))
             topcut1(i)=max(iy+1,topcut1(i))
             leftcut1(i)=ixl
             rightcut1(i)=ixr
          endif

        end do

        nncut=max(nncut,ic)
        if(nncut.gt.DEF_NCUT) then
          write(*,*) ' Increase DEF_NCUT in DIMENSIONS.F!'
          write(*,*) ' nncut = ',nncut,' DEF_NCUT = ',DEF_NCUT
#ifdef BUILDING_CARRE
          stop 'faulty parameter DEF_NCUT'
#else
          call xertst (nncut.le.DEF_NCUT,'faulty parameter DEF_NCUT')
#endif

        end if
      end do

      write(*,*) 'istyle',istyle,'nncut',nncut
      if(nncut.eq.0) write(*,*) 'No cuts found'
      if(nncut.ge.1) write(*,'(1x,a,4i5)')
     & 'Calculated leftcut1, rightcut1, topcut1, bottomcut1 = ',
     & leftcut1(1), rightcut1(1), topcut1(1), bottomcut1(1)
      if(nncut.ge.2) write(*,'(1x,a,4i5)')
     & 'Calculated leftcut2, rightcut2, topcut2, bottomcut2 = ',
     & leftcut1(2), rightcut1(2), topcut1(2), bottomcut1(2)

      return
      end subroutine init_connectivity

      function b2gavg1 (nx1, ny1, fun)
cxpb  This performs an average over the whole domain of valid cells on the secondary grid
      use b2mod_types
      implicit none
      real (kind=R8) :: b2gavg1
      integer nx1, ny1
      real (kind=R8) :: fun(-1:nx1, -1:ny1), denominator
      integer ix, iy, celltype
      logical unused

      b2gavg1 = 0.0_R8
      denominator = 0.0_R8
      do ix = -1, nx1
        do iy = -1, ny1
          celltype = cflags1(ix,iy,CELLFLAG_TYPE)
          unused = (celltype == GRID_UNDEFINED)
     &        .or. (celltype == GRID_EXTERNAL)
     &        .or. (celltype == GRID_DEAD)
          if (.not.unused) then
            b2gavg1 = b2gavg1 + fun(ix,iy)
            denominator = denominator + 1.0_R8
          end if
        end do
      end do
      if (denominator.gt.0.0_R8) b2gavg1 = b2gavg1/denominator

      return
      end function b2gavg1

  ! Interpolate a cell quantity cv to a particular vertex
      function interpolateToVertices1(nx1,ny1,vx_index,cv) result(vx)
      use b2mod_b2cmfs
      use b2mod_indirect
      use b2mod_constants
      use b2mod_cellhelper
      implicit none
      integer, intent(in) :: nx1, ny1, vx_index
      real (kind=R8), intent(in), dimension(-1:nx1,-1:ny1) :: cv
      real (kind=R8), dimension(-1:nx1,-1:ny1) :: vx

      ! internal
      integer :: ix, iy, is, ixx, iyy
      real (kind=R8) :: av, wTot, minVal, maxVal, w, d, centroid(0:1)

      ! For every vertex, find connected cells and average values
      do ix = -1, nx1
        do iy = -1, ny1

          if(isUnusedCell(cflags1(ix,iy,CELLFLAG_TYPE))) then
             vx(ix,iy) = cv(ix,iy)
             cycle
          end if

          wTot = 0.0_R8
          av = 0.0_R8
          minVal = huge(1.0_R8)
          maxVal =-huge(1.0_R8)
          do is = 1, 8
            if ( gmap1%mapCvixVx( gmap1%mapVxI(ix,iy,vx_index), is ) ==
     .           B2_GRID_UNDEFINED .or.
     .           gmap1%mapCviyVx( gmap1%mapVxI(ix,iy,vx_index), is ) ==
     .           B2_GRID_UNDEFINED) cycle
            ixx = gmap1%mapCvixVx( gmap1%mapVxI(ix,iy,vx_index), is )
            iyy = gmap1%mapCviyVx( gmap1%mapVxI(ix,iy,vx_index), is )
            centroid = quadCentroid(crx1(ixx,iyy,0),cry1(ixx,iyy,0),
     .                              crx1(ixx,iyy,1),cry1(ixx,iyy,1),
     .                              crx1(ixx,iyy,2),cry1(ixx,iyy,2),
     .                              crx1(ixx,iyy,3),cry1(ixx,iyy,3))
            d = sqrt ( (centroid(0) - crx1(ix,iy,vx_index))**2 +
     .                 (centroid(1) - cry1(ix,iy,vx_index))**2 )
            w = 1.0_R8/d
            av = av + cv(ixx, iyy) * w
            wTot = wTot + w
            minVal = min(minVal, cv(ixx, iyy))
            maxVal = max(maxVal, cv(ixx, iyy))

          end do

          vx(ix, iy) = av / wTot
          d = vx(ix,iy) - minVal
          if (minVal.ne.0.0_R8) then
#ifdef BUILDING_CARRE
            if ( d/abs(minVal) < - 1.0e-15_R8 ) stop
     .       "Interpolated value smaller than minimum"
#else
            call xertst( d/abs(minVal) >= - 1.0e-15_R8,
     .       "Interpolated value smaller than minimum" )
#endif
          else
#ifdef BUILDING_CARRE
            if ( d.lt.minVal ) stop
     .       "Interpolated value smaller than minimum"
#else
            call xertst( d.ge.minVal,
     .       "Interpolated value smaller than minimum" )
#endif
          end if
          d = vx(ix,iy) - maxVal
          if (maxVal.ne.0.0_R8) then
#ifdef BUILDING_CARRE
            if ( d/abs(maxVal) > 1.0e-15_R8 ) stop
     .       "Interpolated value bigger than maximum"
#else
            call xertst( d/abs(maxVal) <=   1.0e-15_R8,
     .       "Interpolated value bigger than maximum" )
#endif
          else
#ifdef BUILDING_CARRE
            if ( d.gt.maxVal ) stop
     .       "Interpolated value bigger than maximum"
#else
            call xertst( d.le.maxVal,
     .       "Interpolated value bigger than maximum" )
#endif
          end if
        end do
      end do
      return

      end function interpolateToVertices1

      function cr1(ix,iy)
      implicit none
      integer, intent(in) :: ix, iy
      real(kind=R8) :: cr1
      cr1 = 0.25_R8*(crx1(ix,iy,0)+crx1(ix,iy,1)+
     .               crx1(ix,iy,2)+crx1(ix,iy,3))
      return
      end function cr1

      function cz1(ix,iy)
      implicit none
      integer, intent(in) :: ix, iy
      real(kind=R8) :: cz1
      cz1 = 0.25_R8*(cry1(ix,iy,0)+cry1(ix,iy,1)+
     .               cry1(ix,iy,2)+cry1(ix,iy,3))
      return
      end function cz1

      function dist1(ix,iy,i,j)
      use b2mod_math
      implicit none
      integer, intent(in) :: ix, iy, i, j
      real(kind=R8) :: dist1
      dist1 = norm (
     &  crx1(ix,iy,j)-crx1(ix,iy,i), cry1(ix,iy,j)-cry1(ix,iy,i))
      return
      end function dist1

      function distc1(ix1,iy1,ip1,ix2,iy2,ip2)
      use b2mod_math
      implicit none
      integer, intent(in) :: ix1, iy1, ip1, ix2, iy2, ip2
      real(kind=R8) :: distc1
      distc1 = norm(crx1(ix1,iy1,ip1)-crx1(ix2,iy2,ip2),
     .              cry1(ix1,iy1,ip1)-cry1(ix2,iy2,ip2))
      return
      end function distc1

    !! Matches right face of cell (ix1,iy1) to left face of cell (ix2,iy2)
      logical function matchLeft1(ix1,iy1,ix2,iy2)
      implicit none
      integer, intent(in) :: ix1, iy1, ix2, iy2
      matchLeft1 =
     & (distc1(ix1,iy1,1,ix2,iy2,0)+distc1(ix1,iy1,3,ix2,iy2,2)).lt.
     & geom_match_dist
      return
      end function matchLeft1

    !! Matches top face of cell (ix1,iy1) to bottom face of cell (ix2,iy2)
      logical function matchBottom1(ix1,iy1,ix2,iy2)
      implicit none
      integer, intent(in) :: ix1, iy1, ix2, iy2
      matchBottom1 =
     & (distc1(ix1,iy1,2,ix2,iy2,0)+distc1(ix1,iy1,3,ix2,iy2,1)).lt.
     & geom_match_dist
      return
      end function matchBottom1

      end module b2mod_geo2
