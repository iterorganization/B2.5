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
     *  hx1(:,:),
     *  hy1(:,:),
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
      allocate(hx1(-1:nxd,-1:nyd))
      allocate(hy1(-1:nxd,-1:nyd))
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
         deallocate(hx1)
         deallocate(hy1)
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
!      intrinsic huge

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
            call xertst( d/abs(minVal) >= - 1.0e-15_R8,
     .       "Interpolated value smaller than minimum" )
          else
            call xertst( d.ge.minVal,
     .       "Interpolated value smaller than minimum" )
          end if
          d = vx(ix,iy) - maxVal
          if (maxVal.ne.0.0_R8) then
            call xertst( d/abs(maxVal) <=   1.0e-15_R8,
     .       "Interpolated value bigger than maximum" )
          else
            call xertst( d.le.minVal,
     .       "Interpolated value bigger than minimum" )
          end if
        end do
      end do
      return

      end function interpolateToVertices1

      end module b2mod_geo2
