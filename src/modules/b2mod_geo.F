      module b2mod_geo
      use b2mod_types
      use b2mod_math
      use b2mod_version
      use b2mod_indirect
      use b2mod_geo_corner
      use b2mod_cellhelper
      use b2mod_grid_mapping
      implicit none

      real (kind=R8), allocatable, save ::
     *  bb(:,:,:),
     *  crx(:,:,:),
     *  cry(:,:,:),
     *  ffbz(:,:,:),
     *  fpsi(:,:,:),
     *  gs(:,:,:),
     *  hx(:,:),
     *  hy(:,:),
     *  hz(:,:),                                                     !srv 09.01.01
     *  hyl(:,:),                                                    ! IYS 04.03.2019
     *  hxb(:,:),                                                    ! IYS 04.03.2019
     *  hzl(:,:),                                                    ! IYS 04.04.2019
     *  hzb(:,:),                                                    ! IYS 04.04.2019
     *  qz(:,:,:),
     *  qc(:,:,:),
     *  qs(:,:,:),
     *  vol(:,:),
     *  pbs(:,:,:), pbshz(:,:,:),                                    !srv 09.01.01
     *  conn(:),
     *  cell_width(:,:),
     *  cell_height(:,:)

      type(B2GridMap), save :: gmap

      real (kind=R8) :: geom_match_dist
      data geom_match_dist/1.0e-6_R8/
      logical, save :: b2mod_geo_allocated = .false.

      public :: abs_diff, dist, disx, disy, match
      public :: divide_by_poloidal_areas, divide_by_contact_areas
      public :: hy1, pit, pitx, sy

      contains

      subroutine alloc_b2mod_geo(nxd,nyd)
      implicit none
      integer, intent(in) :: nxd, nyd

      allocate(bb(-1:nxd,-1:nyd,0:3))
      allocate(crx(-1:nxd,-1:nyd,0:3))
      allocate(cry(-1:nxd,-1:nyd,0:3))
      allocate(ffbz(-1:nxd,-1:nyd,0:3))
      allocate(fpsi(-1:nxd,-1:nyd,0:3))
      allocate(gs(-1:nxd,-1:nyd,0:2))
      allocate(hx(-1:nxd,-1:nyd))
      allocate(hyl(-1:nxd,-1:nyd))
      allocate(hxb(-1:nxd,-1:nyd))
      allocate(hzl(-1:nxd,-1:nyd))
      allocate(hzb(-1:nxd,-1:nyd))
      allocate(hy(-1:nxd,-1:nyd))
      allocate(hz(-1:nxd,-1:nyd))                                    !srv 09.01.01
      allocate(qz(-1:nxd,-1:nyd,0:1))
      allocate(qc(-1:nxd,-1:nyd,0:1))
      allocate(qs(-1:nxd,-1:nyd,0:1))
      allocate(vol(-1:nxd,-1:nyd))
      allocate(pbs(-1:nxd,-1:nyd,0:1))
      allocate(pbshz(-1:nxd,-1:nyd,0:1))                             !srv 09.01.01
      allocate(conn(-1:nyd))
      allocate(cell_width(-1:nxd,-1:nyd))
      allocate(cell_height(-1:nxd,-1:nyd))

      b2mod_geo_allocated = .true.

      call alloc_b2mod_geo_corner(nxd,nyd)

      return
      end subroutine alloc_b2mod_geo

      subroutine dealloc_b2mod_geo
      implicit none

      if(b2mod_geo_allocated) then

        deallocate(bb)
        deallocate(crx)
        deallocate(cry)
        deallocate(ffbz)
        deallocate(fpsi)
        deallocate(gs)
        deallocate(hx)
        deallocate(hy)
        deallocate(hyl)
        deallocate(hxb)
        deallocate(hzl)
        deallocate(hzb)
        deallocate(hz)
        deallocate(qz)
        deallocate(qc)
        deallocate(qs)
        deallocate(vol)
        deallocate(pbs)
        deallocate(pbshz)                                            !srv 09.01.01
        deallocate(conn)
        deallocate(cell_width)
        deallocate(cell_height)

        b2mod_geo_allocated = .false.

      endif

      if (mapInitialized) call deallocateB2GridMap(gmap)

      call dealloc_b2mod_geo_corner

      return
      end subroutine dealloc_b2mod_geo

      subroutine read_b2mod_geo(nx,ny,nread)
      use b2mod_math
      use b2mod_indirect
      implicit none
      integer nread, nx, ny, n2
      integer ix, iy
      external cfrure, cfrure_opt, xerrab

      n2 = (nx+1)*(ny+1)
      call cfrure (nread, size(bb), bb, 'bb')
      call cfrure (nread, size(crx), crx, 'crx')
      call cfrure (nread, size(cry), cry, 'cry')
      call cfrure (nread, size(ffbz), ffbz, 'ffbz')
      call cfrure (nread, size(fpsi), fpsi, 'fpsi')
      call cfrure (nread, size(gs), gs, 'gs')
      call cfrure (nread, size(hx), hx, 'hx')
      call cfrure (nread, size(hy), hy, 'hy')
      call cfrure (nread, size(qz), qz, 'qz')
      if (b2fplasma_version.ge.'03.001.000') then
        call cfrure (nread, size(qc), qc, 'qc')
        call cfrure (nread, size(qs), qs, 'qs')
      elseif (b2fplasma_version.gt.'01.001.021') then
       call cfrure (nread, n2, qc(-1,-1,0), 'qc')
      else
       do iy = -1, ny
        do ix = -1, nx
         qc(ix,iy,0) = cosang (
     &    (crx(ix,iy,1)+crx(ix,iy,3))-(crx(ix,iy,0)+crx(ix,iy,2)),
     &    (cry(ix,iy,1)+cry(ix,iy,3))-(cry(ix,iy,0)+cry(ix,iy,2)),
     &     crx(ix,iy,2)-crx(ix,iy,0),  cry(ix,iy,2)-cry(ix,iy,0))
         if(qc(ix,iy,0).le.0.0_R8) then
          write(*,*) 'Problem in b2mod_geo: qc(ix,iy) <= 0',
     1     qc(ix,iy,0),ix,iy
          write(*,'(a,4f12.6)') 'crx = ',crx(ix,iy,0:3)
          write(*,'(a,4f12.6)') 'cry = ',cry(ix,iy,0:3)
         endif
        enddo
       enddo
      endif
      if(b2fgmtry_version.lt.'03.001.000') then
       qc(:,:,1) = 1.0_R8
       qs(:,:,1) = 0.0_R8
       if(b2fgmtry_version.ge.'01.001.030') then
	 call cfrure_opt (nread, n2, qc(-1,-1,1), 'qcb')
       end if
       if (maxval(qc(:,:,1)).eq.1.0_R8 .and. 
     &     minval(qc(:,:,1)).eq.1.0_R8) then
        do iy = -1, ny
         do ix = -1, nx
          qc(ix,iy,1) = cosang (
     &      crx(ix,iy,1)-crx(ix,iy,0),  cry(ix,iy,1)-cry(ix,iy,0),
     &     (crx(ix,iy,2)+crx(ix,iy,3))-(crx(ix,iy,0)+crx(ix,iy,1)),
     &     (cry(ix,iy,2)+cry(ix,iy,3))-(cry(ix,iy,0)+cry(ix,iy,1)))
          if(qc(ix,iy,1).le.0.0_R8) then
           write(*,'(a,1f12.6,2i5)')
     1      'Problem in b2mod_geo: qc(ix,iy,1) <= 0',qc(ix,iy,1),ix,iy
           write(*,'(a,4f12.6)') 'crx = ', crx(ix,iy,0:3)
           write(*,'(a,4f12.6)') 'cry = ', cry(ix,iy,0:3)
          endif
         enddo
        enddo
       endif
       do iy = -1, ny
        do ix = -1, nx
         qs(ix,iy,0) = sinang (
     &    (crx(ix,iy,1)+crx(ix,iy,3))-(crx(ix,iy,0)+crx(ix,iy,2)),
     &    (cry(ix,iy,1)+cry(ix,iy,3))-(cry(ix,iy,0)+cry(ix,iy,2)),
     &    crx(ix,iy,2)-crx(ix,iy,0), cry(ix,iy,2)-cry(ix,iy,0))
        enddo
       enddo
cwdk   ensure qc(:,:,0)**2 + qs(:,:,0)**2 = 1, in possible event 
cwdk   of improved computation of qc
       qs(:,:,0) = sign(sqrt(1 - qc(:,:,0)**2),qs(:,:,0))
      endif
      call cfrure (nread, size(vol), vol, 'vol')
c conn
      if(b2fplasma_version.ge.'01.001.063') then        !srv 15.02.11
        call cfrure (nread, size(pbs), pbs, 'pbs')
      else if(b2fplasma_version.ge.'01.001.025') then   !xpb for 5.0 and 5.1 versions
        call cfrure_opt (nread, size(pbs), pbs, 'pbs')
      endif
      if(b2fplasma_version.ge.'01.001.033') then
        call cfrure_opt
     .   (nread, size(cell_width), cell_width, 'cell_width')
        call cfrure_opt
     .   (nread, size(cell_height), cell_height, 'cell_height')
      else
        if(isClassicalGrid(cflags)) then
          cell_width = hx
          cell_height(:,:) = hy(:,:)*qz(:,:,1)
        else
          call xerrab ('Cell widths and heights not set: re-run b2ag!')
        endif
      endif

      return
      end subroutine read_b2mod_geo

      subroutine write_b2mod_geo(nwrite)
      implicit none
      integer nwrite
      external cfwure

      call cfwure (nwrite, size(bb), bb, 'bb')
      call cfwure (nwrite, size(crx), crx, 'crx')
      call cfwure (nwrite, size(cry), cry, 'cry')
      call cfwure (nwrite, size(ffbz), ffbz, 'ffbz')
      call cfwure (nwrite, size(fpsi), fpsi, 'fpsi')
      call cfwure (nwrite, size(gs), gs, 'gs')
      call cfwure (nwrite, size(hx), hx, 'hx')
      call cfwure (nwrite, size(hy), hy, 'hy')
      call cfwure (nwrite, size(qz), qz, 'qz')
      if(newversion.gt.'01.001.021') then
        call cfwure (nwrite, size(qc), qc, 'qc')
      endif
      if(newversion.ge.'03.001.000') then
        call cfwure (nwrite, size(qs), qs, 'qs')
      endif

      call cfwure (nwrite, size(vol), vol, 'vol')
c conn
      if(newversion.ge.'01.001.063') then               !srv 15.02.11
        call cfwure (nwrite, size(pbs), pbs, 'pbs')
      endif
      return
      end subroutine write_b2mod_geo

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
      bottomcut(:)=ny+1
      topcut(:)=-2
      rightcut(:)=-2
      leftcut(:)=nx+1
      ic=0 !! cut counter

      rightix = NO_CONNECTIVITY
      rightiy = NO_CONNECTIVITY
      leftix = NO_CONNECTIVITY
      leftiy = NO_CONNECTIVITY
      topix = NO_CONNECTIVITY
      topiy = NO_CONNECTIVITY
      bottomix = NO_CONNECTIVITY
      bottomiy = NO_CONNECTIVITY
      inseltop = NO_CONNECTIVITY
      inselbot = NO_CONNECTIVITY

    !! First step: find cell connectivity

      do iy=-1,ny
        do ix=-1,nx

            !! unused cells have no connectivity, skip
            if (isUnusedCell(cflag(ix,iy,CELLFLAG_TYPE))) cycle

            !! Look for the left/right connectivity
            if (distc(ix,iy,0,ix,iy,2) > 0.0_R8) then
              ixNb = ix
              iyNb = iy
              cellFound = .false.

              do ystep =  1, ny + 2
                do xstep = 1, nx + 2
                  ixNb = ixNb - 1
                  if (ixNb < -1) ixNb = nx

                  if (isUnusedCell(cflag(ixNb,iyNb,CELLFLAG_TYPE)))
     .             cycle
                  cellFound = matchLeft(ixNb,iyNb,ix,iy)
                  if (cellFound) exit
                end do
                if (cellFound) exit
                iyNb = iyNb - 1
                if (iyNb < -1) iyNb = ny
              end do

              if (cellFound) then
                !! set connectivity in both directions
                leftix(ix,iy)=ixNb
                leftiy(ix,iy)=iyNb
                rightix( ixNb, iyNb ) = ix
                rightiy( ixNb, iyNb ) = iy
              end if
            end if

            !! Look for the top/bottom connectivity
            if (distc(ix,iy,0,ix,iy,1) > 0.0_R8) then
              ixNb = ix
              iyNb = iy
              cellFound = .false.

              do xstep = 1, nx + 2
                do ystep =  1, ny + 2
                  iyNb = iyNb - 1
                  if (iyNb < -1) iyNb = ny

                  if (isUnusedCell(cflag(ixNb,iyNb,CELLFLAG_TYPE)))
     .             cycle
                  cellFound = matchBottom(ixNb,iyNb,ix,iy)
                  if (cellFound) exit
                end do
                if (cellFound) exit
                ixNb = ixNb - 1
                if (ixNb < -1) ixNb = nx
              end do

              if (cellFound) then
                !! set connectivity in both directions
                bottomix(ix,iy)=ixNb
                bottomiy(ix,iy)=iyNb
                topix(ixNb,iyNb)=ix
                topiy(ixNb,iyNb)=iy
              else
                !write (*,*) "No bottom neighbour for ", ix, iy
              end if
            end if

        end do
      end do

    !! Fix connectivity to match B2 convention
      do iy=-1,ny
        do ix=-1,nx
          if (leftix(ix, iy) == NO_CONNECTIVITY) then
            leftix(ix, iy) = -2
            leftiy(ix, iy) = iy
          end if
          if (rightix(ix, iy) == NO_CONNECTIVITY) then
            rightix(ix, iy) = nx + 1
            rightiy(ix, iy) = iy
          end if
          if (bottomix(ix, iy) == NO_CONNECTIVITY) then
            bottomix(ix, iy) = ix
            bottomiy(ix, iy) = -2
          end if
          if (topix(ix, iy) == NO_CONNECTIVITY) then
            topix(ix, iy) = ix
            topiy(ix, iy) = ny + 1
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
     .           leftix(ix, iy), leftiy(ix, iy))) nNb = nNb + 1
            if (isInDomain(nx, ny,
     .           rightix(ix, iy), rightiy(ix, iy))) nNb = nNb + 1
            if (isInDomain(nx, ny,
     .           topix(ix, iy), topiy(ix, iy))) nNb = nNb + 1
            if (isInDomain(nx, ny,
     .           bottomix(ix, iy), bottomiy(ix, iy))) nNb = nNb + 1

            if (nNb < 4) cflag(ix, iy, CELLFLAG_TYPE) = GRID_GUARD

          end do
        end do

       !! mark boundary cells
        do iy=-1,ny
          do ix=-1,nx

            if (.not. isGhostCell(cflag(ix, iy, CELLFLAG_TYPE))) cycle

            if (isInDomain(nx, ny,
     .           leftix(ix, iy), leftiy(ix, iy))) then
              if ( isRealCell(cflag(leftix(ix, iy),
     .                        leftiy(ix, iy), CELLFLAG_TYPE)) ) then
                cflag(leftix(ix, iy),
     .                leftiy(ix, iy), CELLFLAG_TYPE) = GRID_BOUNDARY
                cflag(ix, iy, CELLFLAG_LEFTFACE) = BOUNDARY_NOSTRUCTURE
                cflag(leftix(ix, iy),
     .                leftiy(ix, iy), CELLFLAG_RIGHTFACE) =
     .                                             BOUNDARY_NOSTRUCTURE
              end if
            end if

            if (isInDomain(nx, ny,
     .           rightix(ix, iy), rightiy(ix, iy))) then
              if ( isRealCell(cflag(rightix(ix, iy),
     .                         rightiy(ix, iy), CELLFLAG_TYPE)) ) then
                cflag(rightix(ix, iy),
     .                rightiy(ix, iy), CELLFLAG_TYPE) = GRID_BOUNDARY
                cflag(ix, iy, CELLFLAG_RIGHTFACE) = BOUNDARY_NOSTRUCTURE
                cflag(rightix(ix, iy),
     .                rightiy(ix, iy), CELLFLAG_LEFTFACE) =
     .                                              BOUNDARY_NOSTRUCTURE
              end if
            end if

            if (isInDomain(nx, ny,
     .           topix(ix, iy), topiy(ix, iy))) then
              if ( isRealCell(cflag(topix(ix, iy),
     .                         topiy(ix, iy), CELLFLAG_TYPE)) ) then
                cflag(topix(ix, iy),
     .                topiy(ix, iy), CELLFLAG_TYPE) = GRID_BOUNDARY
                cflag(ix, iy, CELLFLAG_TOPFACE) = BOUNDARY_NOSTRUCTURE
                cflag(topix(ix, iy),
     .                topiy(ix, iy), CELLFLAG_BOTTOMFACE) =
     .                                            BOUNDARY_NOSTRUCTURE
              end if
            end if

            if (isInDomain(nx, ny,
     .           bottomix(ix, iy), bottomiy(ix, iy))) then
              if ( isRealCell(cflag(bottomix(ix, iy),
     .                         bottomiy(ix, iy), CELLFLAG_TYPE)) ) then
                cflag(bottomix(ix, iy),
     .                bottomiy(ix, iy), CELLFLAG_TYPE) = GRID_BOUNDARY
                cflag(ix, iy, CELLFLAG_BOTTOMFACE) =
     .                                             BOUNDARY_NOSTRUCTURE
                cflag(bottomix(ix, iy),
     .                bottomiy(ix, iy), CELLFLAG_TOPFACE) =
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
            if(isInDomain(nx,ny,topix(ix,iy),topiy(ix,iy))) then
              if(cflag(ix,iy,CELLFLAG_LEFTFACE).ne.
     .           cflag(topix(ix,iy),topiy(ix,iy),CELLFLAG_LEFTFACE)
     .         .and. (((cflag(ix,iy,CELLFLAG_LEFTFACE).gt.0 .or.
     .                  cflag(topix(ix,iy),topiy(ix,iy),
     .                              CELLFLAG_LEFTFACE).gt.0) .and.
     .         .not.fullGrid) .or.
     .          (cflag(topix(ix,iy),topiy(ix,iy),CELLFLAG_LEFTFACE)
     .            /= GRID_UNDEFINED .and. fullGrid))) then
                bottomix(topix(ix,iy),topiy(ix,iy)) = ix
                bottomiy(topix(ix,iy),topiy(ix,iy)) = -2
                topix(ix,iy) = ix
                topiy(ix,iy) = ny+1
              end if
            end if
            if(isInDomain(nx,ny,bottomix(ix,iy),
     .                          bottomiy(ix,iy))) then
              if(cflag(ix,iy,CELLFLAG_LEFTFACE).ne.
     .           cflag(bottomix(ix,iy),bottomiy(ix,iy),
     .                       CELLFLAG_LEFTFACE) .and.
     .        (((cflag(ix,iy,CELLFLAG_LEFTFACE).gt.0 .or.
     .           cflag(bottomix(ix,iy),bottomiy(ix,iy),
     .                       CELLFLAG_LEFTFACE).gt.0) .and.
     .         .not.fullGrid) .or.
     .          (cflag(bottomix(ix,iy),bottomiy(ix,iy),
     .                       CELLFLAG_LEFTFACE) /= GRID_UNDEFINED
     .         .and. fullGrid))) then
                topix(bottomix(ix,iy),bottomiy(ix,iy)) = ix
                topiy(bottomix(ix,iy),bottomiy(ix,iy)) = ny+1
                bottomix(ix,iy) = ix
                bottomiy(ix,iy) = -2
              end if
            end if
          end if
          if(cflag(ix,iy,CELLFLAG_BOTTOMFACE) /= GRID_UNDEFINED) then
            if(isInDomain(nx,ny,leftix(ix,iy),leftiy(ix,iy))) then
              if(cflag(ix,iy,CELLFLAG_BOTTOMFACE).ne.
     .           cflag(leftix(ix,iy),
     .                 leftiy(ix,iy),CELLFLAG_BOTTOMFACE) .and.
     .        (((cflag(ix,iy,CELLFLAG_BOTTOMFACE).gt.0 .or.
     .           cflag(leftix(ix,iy),leftiy(ix,iy),
     .                       CELLFLAG_BOTTOMFACE).gt.0) .and.
     .          .not.fullGrid) .or.
     .          (cflag(leftix(ix,iy),leftiy(ix,iy),
     .                       CELLFLAG_BOTTOMFACE) /= GRID_UNDEFINED
     .          .and. fullGrid))) then
                rightix(leftix(ix,iy),leftiy(ix,iy)) = nx+1
                rightiy(leftix(ix,iy),leftiy(ix,iy)) = iy
                leftix(ix,iy) = -2
                leftiy(ix,iy) = iy
              end if
            end if
            if(isInDomain(nx,ny,rightix(ix,iy),rightiy(ix,iy))) then
              if(cflag(ix,iy,CELLFLAG_BOTTOMFACE).ne.
     .           cflag(rightix(ix,iy),rightiy(ix,iy),
     .                       CELLFLAG_BOTTOMFACE) .and.
     .        (((cflag(ix,iy,CELLFLAG_BOTTOMFACE).gt.0 .or.
     .           cflag(rightix(ix,iy),rightiy(ix,iy),
     .                       CELLFLAG_BOTTOMFACE).gt.0) .and.
     .          .not.fullGrid) .or.
     .          (cflag(rightix(ix,iy),rightiy(ix,iy),
     .                       CELLFLAG_BOTTOMFACE) /= GRID_UNDEFINED
     .          .and. fullGrid))) then
                leftix(rightix(ix,iy),rightiy(ix,iy)) = -2
                leftiy(rightix(ix,iy),rightiy(ix,iy)) = iy
                rightix(ix,iy) = nx+1
                rightiy(ix,iy) = -2
              end if
            end if
          end if
          if(cflag(ix,iy,CELLFLAG_RIGHTFACE) /= GRID_UNDEFINED) then
            if(isInDomain(nx,ny,topix(ix,iy),topiy(ix,iy))) then
              if(cflag(ix,iy,CELLFLAG_RIGHTFACE).ne.
     .           cflag(topix(ix,iy),topiy(ix,iy),
     .                       CELLFLAG_RIGHTFACE) .and.
     .        (((cflag(ix,iy,CELLFLAG_RIGHTFACE).gt.0 .or.
     .           cflag(topix(ix,iy),topiy(ix,iy),
     .                       CELLFLAG_RIGHTFACE).gt.0) .and.
     .           .not.fullGrid) .or.
     .          (cflag(topix(ix,iy),topiy(ix,iy),
     .                       CELLFLAG_RIGHTFACE) /= GRID_UNDEFINED
     .           .and. fullGrid))) then
                bottomix(topix(ix,iy),topiy(ix,iy)) = ix
                bottomiy(topix(ix,iy),topiy(ix,iy)) = -2
                topix(ix,iy) = ix
                topiy(ix,iy) = ny+1
              end if
            end if
            if(isInDomain(nx,ny,bottomix(ix,iy),
     .                          bottomiy(ix,iy))) then
              if(cflag(ix,iy,CELLFLAG_RIGHTFACE).ne.
     .           cflag(bottomix(ix,iy),bottomiy(ix,iy),
     .                       CELLFLAG_RIGHTFACE) .and.
     .        (((cflag(ix,iy,CELLFLAG_RIGHTFACE).gt.0 .or.
     .           cflag(bottomix(ix,iy),bottomiy(ix,iy),
     .                       CELLFLAG_RIGHTFACE).gt.0) .and.
     .           .not.fullGrid) .or.
     .          (cflag(bottomix(ix,iy),bottomiy(ix,iy),
     .                       CELLFLAG_RIGHTFACE) /= GRID_UNDEFINED
     .           .and. fullGrid))) then
                topix(bottomix(ix,iy),bottomiy(ix,iy)) = ix
                topiy(bottomix(ix,iy),bottomiy(ix,iy)) = ny+1
                bottomix(ix,iy) = ix
                bottomiy(ix,iy) = -2
              end if
            end if
          end if
          if(cflag(ix,iy,CELLFLAG_TOPFACE) /= GRID_UNDEFINED) then
            if(isInDomain(nx,ny,leftix(ix,iy),leftiy(ix,iy))) then
              if(cflag(ix,iy,CELLFLAG_TOPFACE).ne.
     .           cflag(leftix(ix,iy),leftiy(ix,iy),
     .                       CELLFLAG_TOPFACE) .and.
     .        (((cflag(ix,iy,CELLFLAG_TOPFACE).gt.0 .or.
     .           cflag(leftix(ix,iy),leftiy(ix,iy),
     .                       CELLFLAG_TOPFACE).gt.0) .and.
     .           .not.fullGrid) .or.
     .          (cflag(leftix(ix,iy),leftiy(ix,iy),
     .                       CELLFLAG_TOPFACE) /= GRID_UNDEFINED
     .           .and. fullGrid))) then
                rightix(leftix(ix,iy),leftiy(ix,iy)) = nx+1
                rightiy(leftix(ix,iy),leftiy(ix,iy)) = iy
                leftix(ix,iy) = -2
                leftiy(ix,iy) = iy
              end if
            end if
            if(isInDomain(nx,ny,rightix(ix,iy),
     .                          rightiy(ix,iy))) then
              if(cflag(ix,iy,CELLFLAG_TOPFACE).ne.
     .           cflag(rightix(ix,iy),rightiy(ix,iy),
     .                       CELLFLAG_TOPFACE) .and.
     .        (((cflag(ix,iy,CELLFLAG_TOPFACE).gt.0 .or.
     .           cflag(rightix(ix,iy),rightiy(ix,iy),
     .                       CELLFLAG_TOPFACE).gt.0) .and.
     .           .not.fullGrid) .or.
     .          (cflag(rightix(ix,iy),rightiy(ix,iy),
     .                       CELLFLAG_TOPFACE) /= GRID_UNDEFINED
     .           .and. fullGrid))) then
                leftix(rightix(ix,iy),rightiy(ix,iy)) = -2
                leftiy(rightix(ix,iy),rightiy(ix,iy)) = iy
                rightix(ix,iy) = nx+1
                rightiy(ix,iy) = -2
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
          ixNb = leftix(ix,iy)
          iyNb = leftiy(ix,iy)

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
                if(ixr.eq.rightcut(i) .or.
     .             ixl.eq.leftcut(i)) then
                   cutFound = .true.
                   exit
                endif
             end do
             if (.not.cutFound) then
                ic = ic + 1
                write(*,*) 'ic',ic
                i = ic
             endif
             bottomcut(i)=min(iy,bottomcut(i))
             topcut(i)=max(iy+1,topcut(i))
             leftcut(i)=ixl
             rightcut(i)=ixr
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
     & leftcut(1), rightcut(1), topcut(1), bottomcut(1)
      if(nncut.ge.2) write(*,'(1x,a,4i5)')
     & 'Calculated leftcut2, rightcut2, topcut2, bottomcut2 = ',
     & leftcut(2), rightcut(2), topcut(2), bottomcut(2)

      end subroutine init_connectivity

      subroutine divide_by_contact_areas(nx,ny,field,result)
      implicit none
      integer, intent(in) :: nx, ny
      real(kind=R8), intent(in) :: field(-1:nx,-1:ny,0:1)
      real(kind=R8), intent(out) :: result(-1:nx,-1:ny,0:1)
      integer ix, iy
      do iy = -1, ny
        do ix = -1, nx
          if (gs(ix,iy,0).gt.0.0_R8) then
            result(ix,iy,0) = field(ix,iy,0)/gs(ix,iy,0)
          else
            result(ix,iy,0) = 0.0_R8
          end if
          if (gs(ix,iy,1).gt.0.0_R8) then
            result(ix,iy,1) = field(ix,iy,1)/gs(ix,iy,1)
          else
            result(ix,iy,1) = 0.0_R8
          end if
        end do
      end do
      return
      end subroutine divide_by_contact_areas

      subroutine divide_by_poloidal_areas(nx,ny,field,result)
      implicit none
      integer, intent(in) :: nx, ny
      real(kind=R8), intent(in) :: field(-1:nx,-1:ny,0:1)
      real(kind=R8), intent(out) :: result(-1:nx,-1:ny,0:1)
      integer ix, iy
      do iy = -1, ny
        do ix = -1, nx
          if (gs(ix,iy,0).gt.0.0_R8) then
            result(ix,iy,0) = field(ix,iy,0)/gs(ix,iy,0)/qc(ix,iy,0)
          else
            result(ix,iy,0) = 0.0_R8
          end if
          if (gs(ix,iy,1).gt.0.0_R8) then
            result(ix,iy,1) = field(ix,iy,1)/gs(ix,iy,1)
          else
            result(ix,iy,1) = 0.0_R8
          end if
        end do
      end do
      return
      end subroutine divide_by_poloidal_areas

      function abs_diff(p1,p2)
      implicit none
      real(kind=R8) :: abs_diff, p1, p2
      abs_diff = abs(p1-p2)
      return
      end function abs_diff

      function dist(ix,iy,i,j)
      implicit none
      integer, intent(in) :: ix, iy, i, j
      real(kind=R8) :: dist
      dist = norm(crx(ix,iy,j)-crx(ix,iy,i), cry(ix,iy,j)-cry(ix,iy,i))
      return
      end function dist

      function distc(ix1,iy1,ip1,ix2,iy2,ip2)
      implicit none
      integer, intent(in) :: ix1, iy1, ip1, ix2, iy2, ip2
      real(kind=R8) :: distc
      distc = norm(crx(ix1,iy1,ip1)-crx(ix2,iy2,ip2),
     .             cry(ix1,iy1,ip1)-cry(ix2,iy2,ip2))
      return
      end function distc

      function disx(ix,iy)
      implicit none
      integer, intent(in) :: ix, iy
      real(kind=R8) :: disx
      disx = (dist(bottomix(ix,iy),bottomiy(ix,iy),2,3)+
     1        dist(ix,iy,0,1))/2.0_R8
      return
      end function disx

      function disy(ix,iy)
      implicit none
      integer, intent(in) :: ix, iy
      real(kind=R8) :: disy
      disy = (dist(leftix(ix,iy),leftiy(ix,iy),1,3)+
     1        dist(ix,iy,0,2))/2.0_R8
      return
      end function disy

      function hy1(ix,iy)
      implicit none
      integer, intent(in) :: ix, iy
      real(kind=R8) :: hy1
      hy1 = hy(ix,iy)*qz(ix,iy,1)
      return
      end function hy1

      logical function match(ix1,iy1,ix2,iy2)
      implicit none
      integer, intent(in) :: ix1,iy1,ix2,iy2
      match =
     & (distc(ix1,iy1,1,ix2,iy2,0)+distc(ix1,iy1,3,ix2,iy2,2)).lt.
     & geom_match_dist
      return
      end function match

    !! Matches right face of cell (ix1,iy1) to left face of cell (ix2,iy2)
      logical function matchLeft(ix1,iy1,ix2,iy2)
      implicit none
      integer, intent(in) :: ix1, iy1, ix2, iy2
      matchLeft =
     & (distc(ix1,iy1,1,ix2,iy2,0)+distc(ix1,iy1,3,ix2,iy2,2)).lt.
     & geom_match_dist
      return
      end function matchLeft

    !! Matches top face of cell (ix1,iy1) to bottom face of cell (ix2,iy2)
      logical function matchBottom(ix1,iy1,ix2,iy2)
      implicit none
      integer, intent(in) :: ix1, iy1, ix2, iy2
      matchBottom =
     & (distc(ix1,iy1,2,ix2,iy2,0)+distc(ix1,iy1,3,ix2,iy2,1)).lt.
     & geom_match_dist
      return
      end function matchBottom

      function pit(ix,iy)
      implicit none
      integer, intent(in) :: ix, iy
      real(kind=R8) :: pit
      pit = sqrt(bb(ix,iy,0)**2+bb(ix,iy,1)**2)/bb(ix,iy,3)
      return
      end function pit

      function pitx(ix,iy)
      implicit none
      integer, intent(in) :: ix, iy
      real(kind=R8) :: pitx
      pitx = abs(bb(ix,iy,0))/bb(ix,iy,3)
      return
      end function pitx

      function sy(ix,iy)
      use b2mod_indirect
      implicit none
      integer, intent(in) :: ix, iy
      real(kind=R8) :: sy
      sy = gs(topix(ix,iy),topiy(ix,iy),1)
      return
      end function sy

      function syhz(ix,iy)
      use b2mod_indirect
      implicit none
      integer, intent(in) :: ix, iy
      real(kind=R8) :: syhz
      syhz = gs(topix(ix,iy),topiy(ix,iy),1)*                   !srv 24.08.08
     &       hz(ix,iy)
!     &       hz(topix(ix,iy),topiy(ix,iy))
      return
      end function syhz

      end module b2mod_geo
