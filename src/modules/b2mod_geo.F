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

      public :: abs_diff, dist, disx, disy, distc, match
      public :: divide_by_poloidal_areas, divide_by_contact_areas
      public :: cr, cz, hy1, pit, pitx, sy

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

      subroutine init_region_extended(nx,ny,nncut,nncutmax,
     & leftcut,rightcut,topcut,bottomcut,
     & leftix,leftiy,rightix,rightiy,topix,topiy,bottomix,bottomiy,
     & region,nnreg,resignore,
     & periodic_bc,cflag)
      use b2mod_types
      implicit none
      integer, intent(in) :: nx,ny,nncut,nncutmax
      integer, intent(in) :: leftcut(nncutmax),rightcut(nncutmax),
     & topcut(nncutmax),bottomcut(nncutmax),
     & leftix(-1:nx,-1:ny),leftiy(-1:nx,-1:ny),
     & rightix(-1:nx,-1:ny),rightiy(-1:nx,-1:ny),
     & topix(-1:nx,-1:ny),topiy(-1:nx,-1:ny),
     & bottomix(-1:nx,-1:ny),bottomiy(-1:nx,-1:ny), periodic_bc
      integer, intent(inout) :: cflag(-1:nx,-1:ny,CARREOUT_NCELLFLAGS)

      integer, intent(out) ::
     & region(-1:nx,-1:ny,0:2),nnreg(0:2),
     & resignore(-1:nx,-1:ny,1:2)
    !
    !      Cylindrical slab (core) and annular cases
    !    ++-----------------------------++
    !    ||                             ||
    !    ||                             ||
    !    ||              1              ||
    !    ||                             ||
    !    ||                             ||
    !    ++-----------------------------++
    !
    !    ++-----------------------------++
    !    ||                             ||
    !    ||                             ||
    !    ||1                            ||
    !    ||                             ||
    !    ||                             ||
    !    ++-----------------------------++
    !
    !    ++--------------2--------------++
    !    ||                             ||
    !    ||                             ||
    !    ||                             ||
    !    ||                             ||
    !    ||                             ||
    !    ++--------------1--------------++
    !
    !            Limiter case
    !    +-------------------------------+
    !    |                               |
    !    |               2               |
    !    |                               |
    !    +/-----------------------------\+
    !    ||                             ||
    !    ||              1              ||
    !    ||                             ||
    !    ++-----------------------------++
    !
    !    +-------------------------------+
    !    |                               |
    !    1                               2
    !    |                               |
    !    +/-----------------------------\+
    !    ||                             ||
    !    ||3                            ||
    !    ||                             ||
    !    ++-----------------------------++
    !
    !    +---------------3---------------+
    !    |                               |
    !    |                               |
    !    |                               |
    !    +/--------------2--------------\+
    !    ||                             ||
    !    ||                             ||
    !    ||                             ||
    !    ++--------------1--------------++
    !
    !
    ! Cases with no X-point:
    ! region 1 === core
    ! region 2 === SOL
    !
    !                           +++++++++++++++++++++++++
    !
    !                                       Periodic BC stellarator island case
    !    +-------+---------------+-------+  ++-----------------------------++
    !    |       :               :       |  ||              5              ||
    !    |       :       2       :       |  +\-----------------------------/+
    !    |       :               :       |  |               2               |
    !    |   3   +---------------+   4   |  +-------+---------------+-------+
    !    |       |               |       |  |       |               |       |
    !    |       |       1       |       |  |   3   |       1       |   4   |
    !    |       |               |       |  |       |               |       |
    !    +-------+---------------+-------+  +-------+---------------+-------+
    !
    !
    !    +-------+---------------+-------+  ++-----------------------------++
    !    |       :               :       |  ||7                            ||
    !    |       2               3       |  +\-----------------------------/+
    !    |       :               :       |  |       2               3       |
    !    1       +---------------+       4  1-------+---------------+-------4
    !    |       |               |       |  |       |               |       |
    !    |       |5              |6      |  |       |5              |6      |
    !    |       |               |       |  |       |               |       |
    !    +-------+---------------+-------+  +-------+---------------+-------+
    !
    !
    !    +---5---+-------6-------+---7---+  ++--------------6--------------++
    !    |       :               :       |  ||                             ||
    !    |       :               :       |  +\--------------8--------------/+
    !    |       :               :       |  |                               |
    !    |       +-------4-------+       |  +---5---+-------4-------+---7---+
    !    |       |               |       |  |       |               |       |
    !    |       |               |       |  |       |               |       |
    !    |       |               |       |  |       |               |       |
    !    +---1---+-------2-------+---3---+  +---1---+-------2-------+---3---+
    !
    ! Cases with one X-point:
    ! region 1 === core
    ! region 2 === SOL
    ! region 3 === Western divertor (inboard for LSN, outboard for USN)
    ! region 4 === Eastern divertor (outboard for LSN, inboard for USN)
    ! region 5 === island
    !
    !                           +++++++++++++++++++++++++
    !
    !                      Connected double-null geometry case
    !
    !    +-------+---------------+-------++-------+---------------+-------+
    !    |       :               :       ||       :               :       |
    !    |       :       2       :       ||       :       6       :       |
    !    |       :               :       ||       :               :       |
    !    |   3   +---------------+   4   ||   7   +---------------+   8   |
    !    |       |               |       ||       |               |       |
    !    |       |       1       |       ||       |       5       |       |
    !    |       |               |       ||       |               |       |
    !    +-------+---------------+-------++-------+---------------+-------+
    !
    !
    !    +-------+---------------+-------++-------+---------------+-------+
    !    |       :               :       ||       :               :       |
    !    |       2               3       ||       6               7       |
    !    |       :               :       ||       :               :       |
    !    1       +---------------+       45       +---------------+       8
    !    |       |               |       ||       |               |       |
    !    |       |9              |10     ||       |11             |12     |
    !    |       |               |       ||       |               |       |
    !    +-------+---------------+-------++-------+---------------+-------+
    !
    !
    !    +---5---+-------6-------+---7---++--12---+------13-------+--14---+
    !    |       :               :       ||       :               :       |
    !    |       :               :       ||       :               :       |
    !    |       :               :       ||       :               :       |
    !    |       +-------4-------+       ||       +------11-------+       |
    !    |       |               |       ||       |               |       |
    !    |       |               |       ||       |               |       |
    !    |       |               |       ||       |               |       |
    !    +---1---+-------2-------+---3---++---8---+-------9-------+--10---+
    !
    ! Cases with two X-points:
    ! region 1 === left core
    ! region 2 === left SOL
    ! region 3 === left inboard divertor
    ! region 4 === left outboard divertor
    ! region 5 === right core
    ! region 6 === right SOL
    ! region 7 === right inboard divertor
    ! region 8 === right outboard divertor
    !
    !                           +++++++++++++++++++++++++
    !
    !                      Disconnected double-null geometry case
    !                           (Bottom X-point is active)
    !
    !    +-------+---------------+-------++-------+---------------+-------+
    !    |       :               :       ||       :               :       |
    !    |       :       2       :       ||       :       6       :       |
    !    |       :               +...4...||...7...+               :       |
    !    |       :               |       ||       |               :       |
    !    |...3...+---------------+       ||       +---------------+...8...|
    !    |       |       1       |       ||       |       5       |       |
    !    |       |               |       ||       |               |       |
    !    +-------+---------------+-------++-------+---------------+-------+
    !
    !
    !    +-------+---------------+-------++-------+---------------+-------+
    !    |       :               :       ||       :               :       |
    !    |       2               3       ||       6               7       |
    !    |       :               +.......45.......+               :       |
    !    |       :               |       ||       |13             :       |
    !    1.......+---------------+10     ||       +---------------+.......8
    !    |       |9              |       ||       |11             |12     |
    !    |       |               |       ||       |               |       |
    !    +-------+---------------+-------++-------+---------------+-------+
    !
    !
    !    +---5---+-------6-------+---7---++--12---+------13-------+--14---+
    !    |       :               :       ||       :               :       |
    !    |       :               :       ||       :               :       |
    !    |       :               +.......||.......+               :       |
    !    |       :               |       ||       |               :       |
    !    |.......+-------4-------+       ||       +------11-------+.......|
    !    |       |               |       ||       |               |       |
    !    |       |               |       ||       |               |       |
    !    +---1---+-------2-------+---3---++---8---+-------9-------+--10---+
    !
    !
    ! For reference: CARRE region numbering referring to Fig. 3b in the CARRE96 paper:
    !
    !    +-------+---------------+-------++-------+---------------+-------+
    !    |4444444:444444444444444:4444444||2222222:222222222222222:2222222|
    !    |4444444:444444444444444:4444444||2222222:222222222222222:2222222|
    !    |.......:...............Y---b---||---a---Y...............:.......|
    !    |1111111:111111111111111|3333333||3333333|111111111111111:1111111|
    !    |---f---X-------d-------+3333333||3333333+-------c-------X---e---|
    !    |5555555|666666666666666|3333333||3333333|666666666666666|5555555|
    !    |5555555|666666666666666|3333333||3333333|666666666666666|5555555|
    !    +-------+---------------+-------++-------+---------------+-------+
    !
    !                           - - - - - - - - - - - - -
    !
    !                      Disconnected double-null geometry case
    !                            (Top X-point is active)
    !
    !    +-------+---------------+-------++-------+---------------+-------+
    !    |       :               :       ||       :               :       |
    !    |       :       2       :       ||       :       6       :       |
    !    |...3...+               :       ||       :               +...8...|
    !    |       |               :       ||       :               |       |
    !    |       +---------------+...4...||...7...+---------------+       |
    !    |       |       1       |       ||       |       5       |       |
    !    |       |               |       ||       |               |       |
    !    +-------+---------------+-------++-------+---------------+-------+
    !
    !
    !    +-------+---------------+-------++-------+---------------+-------+
    !    |       :               :       ||       :               :       |
    !    |       2               3       ||       6               7       |
    !    1.......+               :       ||       :               +.......8
    !    |       |13             :       ||       :               |       |
    !    |       +---------------+.......45.......+---------------+12     |
    !    |       |9              |10     ||       |11             |       |
    !    |       |               |       ||       |               |       |
    !    +-------+---------------+-------++-------+---------------+-------+
    !
    !
    !    +---5---+-------6-------+---7---++--12---+------13-------+--14---+
    !    |       :               :       ||       :               :       |
    !    |       :               :       ||       :               :       |
    !    |.......+               :       ||       :               +.......|
    !    |       |               :       ||       :               |       |
    !    |       +-------4-------+.......||.......+------11-------+       |
    !    |       |               |       ||       |               |       |
    !    |       |               |       ||       |               |       |
    !    +---1---+-------2-------+---3---++---8---+-------9-------+--10---+
    !
    ! Cases with two X-points:
    ! region 1 === left core
    ! region 2 === left SOL
    ! region 3 === left inboard divertor
    ! region 4 === left outboard divertor
    ! region 5 === right core
    ! region 6 === right SOL
    ! region 7 === right inboard divertor
    ! region 8 === right outboard divertor
    !
    ! For extended grid cases, the logic is to give the same number as in the
    ! rectangular grid case for the same topological boundary, although this
    ! means mixing x- and y- directed faces for the faces in contact with the
    ! targets and the vacuum vessel walls. To do so, when a y-directed face
    ! should be included in a traditionally x-directed surface (or vice versa),
    ! we assign it a negative number being minus the number of the x-directed
    ! region it should belong to.
    !

      integer ix,iy,ix1,ix2
      integer inseliy,inselix1,inselix2,iyt,geoType,iFace,offset
      integer ixbreak
      integer lefttargetindex(2), righttargetindex(2)
      logical CellToTest
    !
#ifndef BUILDING_CARRE
      call ipgetr ('b2agfs_geom_match_dist', geom_match_dist)
      call xertst (geom_match_dist.gt.0.0_R8,
     &  'faulty argument geom_match_dist')
#endif
    ! if stellarator island, limiter or cylindrical slab case, find periodicity domain coordinate
      if(periodic_bc.eq.1) then
        if (nncut.eq.1) then    ! stellarator island is on North side
            inseliy=ny+1
            inselix1=-2
            inselix2=-2
            do iy=ny,-1,-1
                do ix=2,0,-1
                    if(match(nx-ix,iy,ix-1,iy)) then
                        inselix1=ix-1
                        inselix2=nx-ix
                        inseliy=iy
                        exit
                    endif
                enddo
            enddo
        else if (nncut.eq.0) then   ! limited region is on South side
            inseliy=-2
            inselix1=-2
            inselix2=-2
            do iy = ny-1, 0, -1
              if (inseliy.gt.-2) cycle
              ix = -1
              ix1 = -2
              do while (ix1.eq.-2 .and. ix.lt.nx)
                geoType = cellGeoType(crx(ix,iy,:), cry(ix,iy,:))
                CellToTest =
     &            cflag(ix,iy,CELLFLAG_TYPE) == GRID_INTERNAL .or.
     &           (cflag(ix,iy,CELLFLAG_TYPE) == GRID_BOUNDARY .and.
     &            cflag(ix,iy,CELLFLAG_TOPFACE) /= GRID_UNDEFINED .and.
     &            geoType /= CGEO_TRIA_NOTOP)
                if (.not.CellToTest) then
                  ix = ix + 1
                else
                  ix1 = ix
                endif
              end do
              if(ix1.eq.-2) cycle
              ix = nx
              ix2 = -2
              do while (ix2.eq.-2 .and. ix.gt.-2)
                geoType = cellGeoType(crx(ix,iy,:), cry(ix,iy,:))
                CellToTest =
     &            cflag(ix,iy,CELLFLAG_TYPE) == GRID_INTERNAL .or.
     &           (cflag(ix,iy,CELLFLAG_TYPE) == GRID_BOUNDARY .and.
     &            cflag(ix,iy,CELLFLAG_TOPFACE) /= GRID_UNDEFINED .and.
     &            geoType /= CGEO_TRIA_NOTOP)
                if (.not.CellToTest) then
                  ix = ix - 1
                else
                  ix2 = ix
                endif
              end do
              if(ix2.eq.-2) cycle
              if(ix1.ne.ix2 .and. match(ix2,iy,ix1,iy)) then
                 inselix1=ix1
                 inselix2=ix2
                 if (distc(ix1,iy,0,ix1,iy,2).lt.geom_match_dist .and.
     &               distc(ix2,iy,1,ix2,iy,3).lt.geom_match_dist) then
                   inseliy=iy-1
                 else
                   inseliy=iy
                 endif
              endif
            enddo
            if (inseliy.eq.ny-1) then
            ! Test for cylindrical slab case : then the guard cells match as well
              if (match(inselix2,ny,inselix1,ny)) inseliy = ny
            endif
        else
#ifdef BUILDING_CARRE
          stop 'Unexpected geometry!'
#else
          call xerrab ('Unexpected geometry!')
#endif
        endif
      endif
    !
      region=0
    ! volume component
      if(nx.eq.1) then  ! 1d radial
        do iy=-1,ny
            do ix=-1,nx
              if (.not.isUnusedCell(cflag(ix,iy,CELLFLAG_TYPE)))
     >         region(ix,iy,0)=1
            enddo
        enddo
        nnreg(0)=1
      else if(ny.eq.1) then  ! 1d parallel
        do iy=-1,ny
            do ix=-1,nx
              if (.not.isUnusedCell(cflag(ix,iy,CELLFLAG_TYPE)))
     >         region(ix,iy,0)=1
            enddo
        enddo
        nnreg(0)=1
      else if (nncut.eq.0 .and. periodic_bc.eq.1) then   ! limiter and cylindrical slab cases
        do iy = -1, inseliy
            do ix = inselix1, inselix2
              if (.not.isUnusedCell(cflag(ix,iy,CELLFLAG_TYPE)))
     >         region(ix,iy,0) = 1
            enddo
        enddo
        do iy = inseliy+1, ny
            do ix = -1, nx
              if (.not.isUnusedCell(cflag(ix,iy,CELLFLAG_TYPE)))
     >         region(ix,iy,0) = 2
            enddo
        enddo
        if (inseliy.eq.ny) then
          nnreg(0) = 1   ! cylindrical slab
        else
          nnreg(0) = 2   ! limiter geometry
        end if
      else if (nncut.eq.0 .and. periodic_bc.eq.0) then   ! 2-D slab case
        do iy = -1, ny
            do ix = -1, nx
              if (.not.isUnusedCell(cflag(ix,iy,CELLFLAG_TYPE)))
     >         region(ix,iy,0) = 1
            end do
        end do
        nnreg(0) = 1
      else
        if(topcut(1).ge.-1) then
          if(nncut.eq.1) then
            do ix=-1,leftcut(1)-1
              do iy=-1,ny
                if (isUnusedCell(cflag(ix,iy,CELLFLAG_TYPE))) cycle
                if (periodic_bc.ne.1.or.iy.lt.topcut(1))
     >           region(ix,iy,0)=3
                if (periodic_bc.eq.1.and.
     >              iy.ge.topcut(1).and.iy.lt.inseliy) region(ix,iy,0)=2
                if (periodic_bc.eq.1.and.iy.ge.inseliy.and.
     >              ix.ge.inselix1.and.ix.le.inselix2) region(ix,iy,0)=5
              enddo
            enddo
            do ix=leftcut(1),rightcut(1)-1
              do iy=bottomcut(1),topcut(1)-1
                if (.not.isUnusedCell(cflag(ix,iy,CELLFLAG_TYPE)))
     >           region(ix,iy,0)=1
              enddo
              if(periodic_bc.eq.0) then
                do iy=topcut(1),ny
                  if (.not.isUnusedCell(cflag(ix,iy,CELLFLAG_TYPE)))
     >             region(ix,iy,0)=2
                enddo
              elseif(periodic_bc.eq.1) then
                do iy=topcut(1),inseliy-1
                  if (.not.isUnusedCell(cflag(ix,iy,CELLFLAG_TYPE)))
     >             region(ix,iy,0)=2
                enddo
                do iy=inseliy,ny
                  if (.not.isUnusedCell(cflag(ix,iy,CELLFLAG_TYPE)).and.
     &             ix.ge.inselix1.and.ix.le.inselix2) region(ix,iy,0)=5
                enddo
              endif
            enddo
            do ix=rightcut(1),nx
              do iy=-1,ny
                if (isUnusedCell(cflag(ix,iy,CELLFLAG_TYPE))) cycle
                if (periodic_bc.eq.0.or.iy.lt.topcut(1))
     >           region(ix,iy,0)=4
                if (periodic_bc.eq.1.and.
     &              iy.ge.topcut(1).and.iy.lt.inseliy) region(ix,iy,0)=2
                if (periodic_bc.eq.1.and.iy.ge.inseliy.and.
     &              ix.ge.inselix1.and.ix.le.inselix2) region(ix,iy,0)=5
              enddo
            enddo
            if(periodic_bc.eq.1) then
              nnreg(0)=5
            else
              nnreg(0)=4
            endif
          elseif(nncut.eq.2) then
            do ix = -1, leftcut(1)-1
              do iy = -1, ny
                if (.not.isUnusedCell(cflag(ix,iy,CELLFLAG_TYPE)))
     >           region(ix,iy,0)=3
              enddo
            enddo
            do ix = leftcut(1), leftcut(2)-1
              do iy = max(bottomcut(1),bottomcut(2)), ny
                if (isUnusedCell(cflag(ix,iy,CELLFLAG_TYPE))) cycle
                if (iy.lt.min(topcut(1),topcut(2))) then
                  region(ix,iy,0)=1
                else
                  region(ix,iy,0)=2
                endif
              enddo
            enddo
            do iy = -1, ny
              ix = leftcut(2)
              do while
     &         (isInDomain(nx,ny,leftix(ix,iy),leftiy(ix,iy)) .and.
     &          .not.isUnusedCell(cflag(ix,iy,CELLFLAG_TYPE)))
                region(ix,iy,0)=4
                ix=ix+1
              enddo
              ix = rightcut(2)-1
              do while
     &         (isInDomain(nx,ny,rightix(ix,iy),rightiy(ix,iy)) .and.
     &          .not.isUnusedCell(cflag(ix,iy,CELLFLAG_TYPE)))
                region(ix,iy,0)=7
                ix=ix-1
              enddo
            enddo
            do ix = rightcut(2), rightcut(1)-1
              do iy = max(bottomcut(1),bottomcut(2)), ny
                if (isUnusedCell(cflag(ix,iy,CELLFLAG_TYPE))) cycle
                if (iy.lt.min(topcut(1),topcut(2))) then
                  region(ix,iy,0)=5
                else
                  region(ix,iy,0)=6
                endif
              enddo
            enddo
            do ix = rightcut(1), nx
              do iy = -1, ny
                if (.not.isUnusedCell(cflag(ix,iy,CELLFLAG_TYPE)))
     >           region(ix,iy,0)=8
              enddo
            enddo
            nnreg(0)=8
          endif
        endif
      endif

    ! Ghost cells inherit region number from internal neighbour cell
      do ix=-1,nx
        do iy=-1,ny
          if ( cflag(ix, iy, CELLFLAG_TYPE) /= GRID_BOUNDARY ) cycle
          geoType = cellGeoType(crx(ix,iy,:), cry(ix,iy,:))

          if ( cflag(ix, iy, CELLFLAG_LEFTFACE) /= GRID_UNDEFINED
     &        .and. geoType /= CGEO_TRIA_NOLEFT )
     &     region(leftix(ix,iy), leftiy(ix,iy), 0) = region(ix,iy,0)
          if ( cflag(ix, iy, CELLFLAG_BOTTOMFACE) /= GRID_UNDEFINED
     &        .and. geoType /= CGEO_TRIA_NOBOT )
     &     region(bottomix(ix,iy), bottomiy(ix,iy), 0) = region(ix,iy,0)
          if ( cflag(ix, iy, CELLFLAG_RIGHTFACE) /= GRID_UNDEFINED
     &        .and. geoType /= CGEO_TRIA_NORIGHT )
     &     region(rightix(ix,iy), rightiy(ix,iy), 0) = region(ix,iy,0)
          if ( cflag(ix, iy, CELLFLAG_TOPFACE) /= GRID_UNDEFINED
     &        .and. geoType /= CGEO_TRIA_NOTOP )
     &     region(topix(ix,iy), topiy(ix,iy), 0) = region(ix,iy,0)
        end do
      end do

    ! Set unused cells to region 0
      do iy=-1,ny
        do ix=-1,nx
          if (isUnusedCell(cflag(ix,iy,CELLFLAG_TYPE)))
     >     region(ix,iy,0) = 0
        end do
      end do

    ! Check all valid cells have a nonzero region number
      do ix=-1,nx
        do iy=-1,ny
          if(.not.isUnusedCell(cflag(ix,iy,CELLFLAG_TYPE)) .and.
     &     region(ix,iy,0).eq.0)
     >     write(*,*) 'REGION not set for ',ix,iy,0
        enddo
      enddo

    ! After we know the volumetric region numbers, we fix the boundary indices for
    ! non-structure boundaries. CARRE2 sets them to -1 for all boundary faces not
    ! on a structure. Here they are set to -REGION.
    ! The corresponding face of the guard cell receives the same number.

      do iy=-1,ny
        do ix=-1,nx
          if (cflag(ix, iy, CELLFLAG_TYPE) /= GRID_BOUNDARY) cycle
          do iFace = CELLFLAG_LEFTFACE, CELLFLAG_TOPFACE
            if (cflag(ix, iy, iFace) == BOUNDARY_NOSTRUCTURE) then
              offset = (iFace - CELLFLAG_LEFTFACE + 1) * (-10)
              cflag(ix, iy, iFace) = offset - region(ix, iy, 0)
              geoType = cellGeoType(crx(ix,iy,:), cry(ix,iy,:))
              if (iFace == CELLFLAG_LEFTFACE .and.
     &            geoType /= CGEO_TRIA_NOLEFT) then
                cflag(leftix(ix,iy),leftiy(ix,iy),CELLFLAG_RIGHTFACE) =
     &           cflag(ix, iy, iFace)
              else if (iFace == CELLFLAG_RIGHTFACE .and.
     &                 geoType /= CGEO_TRIA_NORIGHT) then
                cflag(rightix(ix,iy),rightiy(ix,iy),CELLFLAG_LEFTFACE) =
     &           cflag(ix, iy, iFace)
              else if (iFace == CELLFLAG_TOPFACE .and.
     &                 geoType /= CGEO_TRIA_NOTOP) then
                cflag(topix(ix,iy),topiy(ix,iy),CELLFLAG_BOTTOMFACE) =
     &           cflag(ix, iy, iFace)
              else if (iFace == CELLFLAG_BOTTOMFACE .and.
     &                 geoType /= CGEO_TRIA_NOBOT) then
                cflag(bottomix(ix,iy),bottomiy(ix,iy),CELLFLAG_TOPFACE)=
     &           cflag(ix, iy, iFace)
              endif
            end if
          end do
        enddo
      enddo

      if(nncut.eq.0 .and. periodic_bc.eq.1 .and. nnreg(0).eq.2) then   ! limiter case
! We identify the structure index of the left target
        iy = inseliy+1
        ix = inselix1
        lefttargetindex(1)=GRID_UNDEFINED
        do while (ix.ge.-1 .and. lefttargetindex(1).eq.GRID_UNDEFINED)
          if (cflag(ix, iy, CELLFLAG_TYPE) == GRID_BOUNDARY) then
            geoType = cellGeoType(crx(ix,iy,:), cry(ix,iy,:))
            if (geoType /= CGEO_TRIA_NOLEFT .and.
     &          cflag(ix, iy, CELLFLAG_LEFTFACE) /= GRID_UNDEFINED) then
              lefttargetindex(1) = cflag(ix,iy,CELLFLAG_LEFTFACE)
            else if (geoType /= CGEO_TRIA_NOTOP .and.
     &          cflag(ix, iy, CELLFLAG_TOPFACE) /= GRID_UNDEFINED) then
              lefttargetindex(1) = cflag(ix,iy,CELLFLAG_TOPFACE)
            else if (geoType /= CGEO_TRIA_NOBOT .and.
     &          cflag(ix, iy, CELLFLAG_BOTTOMFACE) /= GRID_UNDEFINED)
     &       then
              lefttargetindex(1) = cflag(ix,iy,CELLFLAG_BOTTOMFACE)
            else
#ifdef BUILDING_CARRE
              stop
     &  'Invalid cflags settings (L), should not have gotten here!'
#else
              call xerrab
     & ('Invalid cflags settings (L), should not have gotten here!')
#endif
            endif
          else
            ix = ix-1
          endif
        end do
#ifdef BUILDING_CARRE
        if (.not. lefttargetindex(1).ne.GRID_UNDEFINED)
     >   stop 'Left target not found!'
#else
        call xertst(lefttargetindex(1).ne.GRID_UNDEFINED,
     .        'Left target not found!')
#endif
! We identify the structure index of the right target
        iy = inseliy+1
        ix = inselix2
        righttargetindex(1)=GRID_UNDEFINED
        do while (ix.le.nx .and. righttargetindex(1).eq.GRID_UNDEFINED)
          if (cflag(ix, iy, CELLFLAG_TYPE) == GRID_BOUNDARY) then
            geoType = cellGeoType(crx(ix,iy,:), cry(ix,iy,:))
            if (geoType /= CGEO_TRIA_NORIGHT .and.
     &       cflag(ix, iy, CELLFLAG_RIGHTFACE) /= GRID_UNDEFINED) then
              righttargetindex(1) = cflag(ix,iy,CELLFLAG_RIGHTFACE)
            else if (geoType /= CGEO_TRIA_NOTOP .and.
     &       cflag(ix, iy, CELLFLAG_TOPFACE) /= GRID_UNDEFINED) then
              righttargetindex(1) = cflag(ix,iy,CELLFLAG_TOPFACE)
            else if (geoType /= CGEO_TRIA_NOBOT .and.
     &       cflag(ix, iy, CELLFLAG_BOTTOMFACE) /= GRID_UNDEFINED) then
              righttargetindex(1) = cflag(ix,iy,CELLFLAG_BOTTOMFACE)
            else
#ifdef BUILDING_CARRE
              stop
     &  'Invalid cflags settings (R), should not have gotten here!'
#else
              call xerrab
     & ('Invalid cflags settings (R), should not have gotten here!')
#endif
            endif
          else
            ix = ix+1
          endif
        end do
#ifdef BUILDING_CARRE
        if (.not. righttargetindex(1).ne.GRID_UNDEFINED)
     >   stop 'Right target not found!'
#else
        call xertst(righttargetindex(1).ne.GRID_UNDEFINED,
     .        'Right target not found!')
#endif

! X-flux component
        do iy = inseliy+1, ny
! We find the cells contacting the left target: x-region 1
          ix = -1
          do while (ix.lt.nx
     &     .and. cflag(ix, iy, CELLFLAG_TYPE) /= GRID_INTERNAL)
            if (cflag(ix, iy, CELLFLAG_TYPE) == GRID_BOUNDARY) then
              geoType = cellGeoType(crx(ix,iy,:), cry(ix,iy,:))
              if (cflag(ix, iy, CELLFLAG_LEFTFACE) == lefttargetindex(1)
     &            .and. geoType /= CGEO_TRIA_NOLEFT) region(ix,iy,1)=1
              if (cflag(ix, iy, CELLFLAG_TOPFACE) == lefttargetindex(1)
     &            .and. geoType /= CGEO_TRIA_NOTOP)
     &         region(topix(ix,iy),topiy(ix,iy),2)=-1
              if (cflag(ix, iy, CELLFLAG_BOTTOMFACE) ==
     &                                               lefttargetindex(1)
     &            .and. geoType /= CGEO_TRIA_NOBOT) region(ix,iy,2)=-1
            endif
            ix = ix + 1
          enddo
! We find the cells contacting the right target: x-region 2
          ix = nx
          do while (ix.gt.-1
     &     .and. cflag(ix, iy, CELLFLAG_TYPE) /= GRID_INTERNAL)
            if (cflag(ix, iy, CELLFLAG_TYPE) == GRID_BOUNDARY) then
              geoType = cellGeoType(crx(ix,iy,:), cry(ix,iy,:))
              if (cflag(ix, iy, CELLFLAG_RIGHTFACE) ==
     &            righttargetindex(1) .and.
     &            geoType /= CGEO_TRIA_NORIGHT)
     &         region(rightix(ix,iy),rightiy(ix,iy),1)=2
              if (cflag(ix, iy, CELLFLAG_TOPFACE) == righttargetindex(1)
     &            .and. geoType /= CGEO_TRIA_NOTOP)
     &         region(topix(ix,iy),topiy(ix,iy),2)=-2
              if (cflag(ix, iy, CELLFLAG_BOTTOMFACE) ==
     &            righttargetindex(1) .and. geoType /= CGEO_TRIA_NOBOT)
     &         region(ix,iy,2)=-2
            endif
            ix = ix - 1
          enddo
        enddo
! We find the core periodicity line: x-region 3
        do iy = -1, inseliy
          if (.not.isUnusedCell(cflag(inselix1,iy,CELLFLAG_TYPE)))
     >     region(inselix1,iy,1)=3
        enddo
        nnreg(1)=3
! Y-flux component
        do ix = inselix1, inselix2
! We find the core boundary: y-region 1
          iy = -1
          do while (iy.le.inseliy
     &     .and. .not.isRealCell(cflag(ix, iy, CELLFLAG_TYPE)))
            iy = iy+1
          end do
          if (iy.le.inseliy) then
            region(ix,iy,2)=1
          endif
! We find the separatrix: y-region 2
          region(topix(ix,inseliy),topiy(ix,inseliy),2)=2
        enddo
! We find the contact to the main chamber wall: y-region 3
        do ix = -1, nx
          iy = ny
          do while (iy.gt.inseliy
     &     .and. cflag(ix,iy,CELLFLAG_TYPE) /= GRID_INTERNAL)
            if (cflag(ix, iy, CELLFLAG_TYPE) == GRID_BOUNDARY) then
              geoType = cellGeoType(crx(ix,iy,:), cry(ix,iy,:))
              if (cflag(ix, iy, CELLFLAG_TOPFACE) /= lefttargetindex(1)
     &         .and.
     &            cflag(ix, iy, CELLFLAG_TOPFACE) /= righttargetindex(1)
     &         .and.
     &            cflag(ix, iy, CELLFLAG_TOPFACE) /= GRID_UNDEFINED
     &         .and. geoType /= CGEO_TRIA_NOTOP)
     &         region(topix(ix,iy),topiy(ix,iy),2)=3
              if
     &         (cflag(ix, iy, CELLFLAG_RIGHTFACE) /= lefttargetindex(1)
     &          .and.
     &          cflag(ix, iy, CELLFLAG_RIGHTFACE) /= righttargetindex(1)
     &          .and.
     &          cflag(ix, iy, CELLFLAG_RIGHTFACE) /= GRID_UNDEFINED
     &          .and. geoType /= CGEO_TRIA_NORIGHT)
     &         region(rightix(ix,iy),rightiy(ix,iy),1)=-3
              if
     &         (cflag(ix, iy, CELLFLAG_LEFTFACE) /= lefttargetindex(1)
     &         .and.
     &          cflag(ix, iy, CELLFLAG_LEFTFACE) /= righttargetindex(1)
     &         .and.
     &          cflag(ix, iy, CELLFLAG_LEFTFACE) /= GRID_UNDEFINED
     &         .and. geoType /= CGEO_TRIA_NOLEFT) region(ix,iy,1)=-3
            end if
            iy = iy - 1
          end do
        end do
        nnreg(2)=3
      else if (nncut.eq.0 .and. periodic_bc.eq.1 .and. nnreg(0).eq.1)
     > then ! cylindrical slab case
        region(inselix1,:,1)=1
        nnreg(1)=1
        do ix = inselix1, inselix2
          region(topix(ix,-1),topiy(ix,-1),2)=1
          region(ix,ny,2)=1
        enddo
        nnreg(2)=2
      elseif(topcut(1).ge.-1) then   ! we have at least 1 X-point
        if(periodic_bc.eq.1) then
          iyt=inseliy-1
        else
          iyt=ny
        endif
! We identify the structure index of the left target
        iy = topcut(1)+1
        ix = leftcut(1)-1
        lefttargetindex(1)=GRID_UNDEFINED
        do while (ix.gt.-1 .and. lefttargetindex(1).eq.GRID_UNDEFINED)
          if (cflag(ix, iy, CELLFLAG_TYPE) == GRID_BOUNDARY) then
            geoType = cellGeoType(crx(ix,iy,:), cry(ix,iy,:))
            if (geoType /= CGEO_TRIA_NOLEFT .and.
     &       cflag(ix, iy, CELLFLAG_LEFTFACE) /= GRID_UNDEFINED) then
              lefttargetindex(1) = cflag(ix,iy,CELLFLAG_LEFTFACE)
            else if (geoType /= CGEO_TRIA_NOTOP .and.
     &       cflag(ix, iy, CELLFLAG_TOPFACE) /= GRID_UNDEFINED) then
              lefttargetindex(1) = cflag(ix,iy,CELLFLAG_TOPFACE)
            else if (geoType /= CGEO_TRIA_NOBOT .and.
     &       cflag(ix, iy, CELLFLAG_BOTTOMFACE) /= GRID_UNDEFINED) then
              lefttargetindex(1) = cflag(ix,iy,CELLFLAG_BOTTOMFACE)
            else
#ifdef BUILDING_CARRE
              stop
     &  'Invalid cflags settings (L), should not have gotten here!'
#else
              call xerrab
     & ('Invalid cflags settings (L), should not have gotten here!')
#endif
            endif
          else
            ix = ix-1
          endif
        end do
#ifdef BUILDING_CARRE
        if (.not. lefttargetindex(1).ne.GRID_UNDEFINED)
     &   stop 'Left target not found!'
#else
        call xertst(lefttargetindex(1).ne.GRID_UNDEFINED,
     &        'Left target not found!')
#endif
! We identify the structure index of the right target
        iy = topcut(1)+1
        ix = rightcut(1)
        righttargetindex(1)=GRID_UNDEFINED
        do while (ix.lt.nx .and. righttargetindex(1).eq.GRID_UNDEFINED)
          if (cflag(ix, iy, CELLFLAG_TYPE) == GRID_BOUNDARY) then
            geoType = cellGeoType(crx(ix,iy,:), cry(ix,iy,:))
            if (geoType /= CGEO_TRIA_NORIGHT .and.
     &       cflag(ix, iy, CELLFLAG_RIGHTFACE) /= GRID_UNDEFINED) then
              righttargetindex(1) = cflag(ix,iy,CELLFLAG_RIGHTFACE)
            else if (geoType /= CGEO_TRIA_NOTOP .and.
     &       cflag(ix, iy, CELLFLAG_TOPFACE) /= GRID_UNDEFINED) then
              righttargetindex(1) = cflag(ix,iy,CELLFLAG_TOPFACE)
            else if (geoType /= CGEO_TRIA_NOBOT .and.
     &       cflag(ix, iy, CELLFLAG_BOTTOMFACE) /= GRID_UNDEFINED) then
              righttargetindex(1) = cflag(ix,iy,CELLFLAG_BOTTOMFACE)
            else
#ifdef BUILDING_CARRE
              stop
     &  'Invalid cflags settings (R), should not have gotten here!'
#else
              call xerrab
     & ('Invalid cflags settings (R), should not have gotten here!')
#endif
            endif
          else
            ix = ix+1
          endif
        end do
#ifdef BUILDING_CARRE
        if (.not. righttargetindex(1).ne.GRID_UNDEFINED)
     >   stop 'Right target not found!'
#else
        call xertst(righttargetindex(1).ne.GRID_UNDEFINED,
     .        'Right target not found!')
#endif

! X-flux component
      do iy=-1,iyt
! We find the cells contacting the left target: x-region 1
        ix = -1
        do while (ix.lt.nx .and.
     &   cflag(ix, iy, CELLFLAG_TYPE) /= GRID_INTERNAL)
          if (cflag(ix, iy, CELLFLAG_TYPE) == GRID_BOUNDARY) then
            geoType = cellGeoType(crx(ix,iy,:), cry(ix,iy,:))
            if (cflag(ix, iy, CELLFLAG_LEFTFACE) == lefttargetindex(1)
     &       .and. geoType /= CGEO_TRIA_NOLEFT) region(ix,iy,1)=1
            if (cflag(ix, iy, CELLFLAG_TOPFACE) == lefttargetindex(1)
     &       .and. geoType /= CGEO_TRIA_NOTOP)
     &       region(topix(ix,iy),topiy(ix,iy),2)=-1
            if (cflag(ix, iy, CELLFLAG_BOTTOMFACE) == lefttargetindex(1)
     &       .and. geoType /= CGEO_TRIA_NOBOT) region(ix,iy,2)=-1
          endif
          ix = ix + 1
        enddo
! We find the cells contacting the right target: x-region 4 (1 X-point) or 8 (2 X-points)
        ix = nx
        do while (ix.gt.-1
     &   .and. cflag(ix, iy, CELLFLAG_TYPE) /= GRID_INTERNAL)
          if (cflag(ix, iy, CELLFLAG_TYPE) == GRID_BOUNDARY) then
            geoType = cellGeoType(crx(ix,iy,:), cry(ix,iy,:))
            if (cflag(ix, iy, CELLFLAG_RIGHTFACE) == righttargetindex(1)
     &       .and. geoType /= CGEO_TRIA_NORIGHT)
     &       region(rightix(ix,iy),rightiy(ix,iy),1)=4*nncut
            if (cflag(ix, iy, CELLFLAG_TOPFACE) == righttargetindex(1)
     &       .and. geoType /= CGEO_TRIA_NOTOP)
     &       region(topix(ix,iy),topiy(ix,iy),2)=-4*nncut
            if (cflag(ix, iy, CELLFLAG_BOTTOMFACE) ==
     &                                              righttargetindex(1)
     &       .and. geoType /= CGEO_TRIA_NOBOT) region(ix,iy,2)=-4*nncut
          endif
          ix = ix - 1
        enddo
      enddo
      if (nncut.eq.2) then  ! Additional data for DN cases
! We identify the location of the break between the two grid halves
        ixbreak = -2
        ix = leftcut(2)
        do while (ix.lt.rightcut(2) .and. ixbreak.eq.-2)
          iy = ny
          do while (.not.isRealCell(cflag(ix, iy, CELLFLAG_TYPE)))
            iy = iy - 1
            if (iy.eq.-2) exit
          end do
          if (iy.ne.-2) then
            ix = ix + 1
          else
            ixbreak = ix
          end if
        end do
#ifdef BUILDING_CARRE
        if (ixbreak.eq.-2) stop 'Did not find the break location!'
#else
        call xertst(ixbreak.ne.-2, 'Did not find the break location!')
#endif
! We identify the structure index of the second left target
        iy = topcut(2)+1
        ix = rightcut(2)-1
        do while (.not.isRealCell(cflag(ix, iy, CELLFLAG_TYPE)) .and.
     &             ix.lt.rightcut(1))
          ix = ix + 1
        end do
#ifdef BUILDING_CARRE
        if (ix.eq.rightcut(1))
     &   stop 'Second separatrix does not have a left strike point!'
#else
        call xertst(ix.lt.rightcut(1),
     &        'Second separatrix does not have a left strike point!')
#endif
        lefttargetindex(2)=GRID_UNDEFINED
        do while (ix.gt.ixbreak
     &   .and. lefttargetindex(2).eq.GRID_UNDEFINED)
          if (cflag(ix, iy, CELLFLAG_TYPE) == GRID_BOUNDARY) then
            geoType = cellGeoType(crx(ix,iy,:), cry(ix,iy,:))
            if (geoType /= CGEO_TRIA_NOLEFT .and.
     &       cflag(ix, iy, CELLFLAG_LEFTFACE) /= GRID_UNDEFINED) then
              lefttargetindex(2) = cflag(ix,iy,CELLFLAG_LEFTFACE)
            else if (geoType /= CGEO_TRIA_NOTOP .and.
     &       cflag(ix, iy, CELLFLAG_TOPFACE) /= GRID_UNDEFINED) then
              lefttargetindex(2) = cflag(ix,iy,CELLFLAG_TOPFACE)
            else if (geoType /= CGEO_TRIA_NOBOT .and.
     &       cflag(ix, iy, CELLFLAG_BOTTOMFACE) /= GRID_UNDEFINED) then
              lefttargetindex(2) = cflag(ix,iy,CELLFLAG_BOTTOMFACE)
            else
#ifdef BUILDING_CARRE
              stop
     >  'Invalid cflags settings (L2), should not have gotten here!'
#else
              call xerrab
     . ('Invalid cflags settings (L2), should not have gotten here!')
#endif
            endif
          else
            ix = ix-1
          endif
        end do
#ifdef BUILDING_CARRE
        if (.not. lefttargetindex(2).ne.GRID_UNDEFINED)
     &   stop 'Second left target not found!'
#else
        call xertst(lefttargetindex(2).ne.GRID_UNDEFINED,
     &        'Second left target not found!')
#endif
! We identify the structure index of the second right target
        iy = topcut(2)+1
        ix = leftcut(2)
        do while (.not.isRealCell(cflag(ix, iy, CELLFLAG_TYPE)) .and.
     &   ix.gt.leftcut(1) )
          ix = ix - 1
        end do
#ifdef BUILDING_CARRE
        if (ix.eq.leftcut(1))
     &   stop 'Second separatrix does not have a right strike point!'
#else
        call xertst(ix.gt.leftcut(1),
     &        'Second separatrix does not have a right strike point!')
#endif
        righttargetindex(2)=GRID_UNDEFINED
        do while (ix.lt.ixbreak
     &   .and. righttargetindex(2).eq.GRID_UNDEFINED)
          if (cflag(ix, iy, CELLFLAG_TYPE) == GRID_BOUNDARY) then
            geoType = cellGeoType(crx(ix,iy,:), cry(ix,iy,:))
            if (geoType /= CGEO_TRIA_NORIGHT .and.
     &       cflag(ix, iy, CELLFLAG_RIGHTFACE) /= GRID_UNDEFINED) then
              righttargetindex(2) = cflag(ix,iy,CELLFLAG_RIGHTFACE)
            else if (geoType /= CGEO_TRIA_NOTOP .and.
     &       cflag(ix, iy, CELLFLAG_TOPFACE) /= GRID_UNDEFINED) then
              righttargetindex(2) = cflag(ix,iy,CELLFLAG_TOPFACE)
            else if (geoType /= CGEO_TRIA_NOBOT .and.
     &       cflag(ix, iy, CELLFLAG_BOTTOMFACE) /= GRID_UNDEFINED) then
              righttargetindex(2) = cflag(ix,iy,CELLFLAG_BOTTOMFACE)
            else
#ifdef BUILDING_CARRE
              stop
     >  'Invalid cflags settings (R2), should not have gotten here!'
#else
              call xerrab
     . ('Invalid cflags settings (R2), should not have gotten here!')
#endif
            endif
          else
            ix = ix+1
          endif
        end do
#ifdef BUILDING_CARRE
        if (.not. righttargetindex(2).ne.GRID_UNDEFINED)
     &   stop 'Second right target not found!'
#else
        call xertst(righttargetindex(2).ne.GRID_UNDEFINED,
     &        'Second right target not found!')
#endif
        do iy=-1,iyt
! We find the cells contacting the second left target: x-region 5
          ix = nx
          do while (ix.gt.ixbreak)
            if (cflag(ix, iy, CELLFLAG_TYPE) == GRID_BOUNDARY) then
              geoType = cellGeoType(crx(ix,iy,:), cry(ix,iy,:))
              if (cflag(ix, iy, CELLFLAG_LEFTFACE) == lefttargetindex(2)
     &         .and. geoType /= CGEO_TRIA_NOLEFT) region(ix,iy,1)=5
              if (cflag(ix, iy, CELLFLAG_TOPFACE) == lefttargetindex(2)
     &         .and. geoType /= CGEO_TRIA_NOTOP)
     &         region(topix(ix,iy),topiy(ix,iy),2)=-5
              if (cflag(ix, iy, CELLFLAG_BOTTOMFACE) ==
     &                                               lefttargetindex(2)
     &         .and. geoType /= CGEO_TRIA_NOBOT) region(ix,iy,2)=-5
            endif
            ix = ix - 1
          enddo
! We find the cells contacting the second right target: x-region 4
          ix = -1
          do while (ix.lt.ixbreak)
            if (cflag(ix, iy, CELLFLAG_TYPE) == GRID_BOUNDARY) then
              geoType = cellGeoType(crx(ix,iy,:), cry(ix,iy,:))
              if (cflag(ix, iy, CELLFLAG_RIGHTFACE) ==
     &                                               righttargetindex(2)
     &         .and. geoType /= CGEO_TRIA_NORIGHT)
     &         region(rightix(ix,iy),rightiy(ix,iy),1)=4
              if (cflag(ix, iy, CELLFLAG_TOPFACE) == righttargetindex(2)
     &         .and. geoType /= CGEO_TRIA_NOTOP)
     &         region(topix(ix,iy),topiy(ix,iy),2)=-4
              if (cflag(ix, iy, CELLFLAG_BOTTOMFACE) ==
     &                                               righttargetindex(2)
     &         .and. geoType /= CGEO_TRIA_NOBOT) region(ix,iy,2)=-4
            endif
            ix = ix + 1
          enddo
        enddo
      end if  ! DN
      do iy=topcut(1),iyt
! We find the location of the Western divertor throat: x-region 2
        geoType = cellGeoType(crx(leftcut(1),iy,:),
     &                        cry(leftcut(1),iy,:))
        if (isRealCell(cflag(leftcut(1),iy,CELLFLAG_TYPE))
     &   .and. geoType /= CGEO_TRIA_NOLEFT) region(leftcut(1),iy,1) = 2
! We find the location of the Eastern divertor throat: x-region 3 (1 X-point) or 7 (2 X-points)
        geoType = cellGeoType(crx(rightcut(1),iy,:),
     &                        cry(rightcut(1),iy,:))
        if (isRealCell(cflag(rightcut(1),iy,CELLFLAG_TYPE))
     &   .and. geoType /= CGEO_TRIA_NOLEFT)
     &   region(rightcut(1),iy,1) = 3 + 4*(nncut-1)
      end do
      if (nncut.eq.2) then
        do iy=topcut(2),iyt
! We find the location of the inner upper divertor throat: x-region 3
          geoType = cellGeoType(crx(leftcut(2),iy,:),
     &                          cry(leftcut(2),iy,:))
          if (isRealCell(cflag(leftcut(2),iy,CELLFLAG_TYPE))
     &     .and. geoType /= CGEO_TRIA_NOLEFT)
     &     region(leftcut(2),iy,1) = 3
! We find the location of the outer upper divertor throat: x-region 6
          geoType = cellGeoType(crx(rightcut(2),iy,:),
     &                          cry(rightcut(2),iy,:))
          if (isRealCell(cflag(rightcut(2),iy,CELLFLAG_TYPE))
     &     .and. geoType /= CGEO_TRIA_NOLEFT)
     &     region(rightcut(2),iy,1) = 6
        end do
      end if
! We find the core and PFR bottom periodicity lines: x-regions 5 and 6 (1 X-point)
!                                                  : x-regions 9 and 12 (2 X-points)
! x-region 13 is the connection between the two halves of the SOL section between the two separatrices
      do iy=bottomcut(1),topcut(1)-1
        if(nncut.eq.1) then
          if (isRealCell(cflag(leftcut(1),iy,CELLFLAG_TYPE)))
     &     region(leftcut(1),iy,1)=5
          geoType = cellGeoType(crx(rightcut(1),iy,:),
     &                          cry(rightcut(1),iy,:))
          if (isRealCell(cflag(rightcut(1),iy,CELLFLAG_TYPE)) .and.
     &     cflag(rightcut(1), iy, CELLFLAG_LEFTFACE) == GRID_UNDEFINED
     &     .and. geoType /= CGEO_TRIA_NOLEFT)
     &     region(rightcut(1),iy,1) = 6
        elseif(nncut.eq.2) then
          if(iy.ge.bottomcut(2).and.iy.lt.topcut(2)) then
            if (isRealCell(cflag(leftcut(1),iy,CELLFLAG_TYPE)))
     >       region(leftcut(1),iy,1)=9
          else
            if (isRealCell(cflag(leftcut(1),iy,CELLFLAG_TYPE)))
     >       region(leftcut(1),iy,1)=13
          endif
          geoType = cellGeoType(crx(rightcut(1),iy,:),
     &                          cry(rightcut(1),iy,:))
          if (isRealCell(cflag(rightcut(1),iy,CELLFLAG_TYPE)) .and.
     &     cflag(rightcut(1), iy, CELLFLAG_LEFTFACE) == GRID_UNDEFINED)
     &     region(rightcut(1),iy,1)=12
        endif
      enddo
! We find the core and PFR top periodicity lines: x-regions 10 and 11 (2 X-points)
! x-region 13 is the connection between the two halves of the SOL section between the two separatrices
      if(nncut.eq.2) then
        do iy=bottomcut(2),topcut(2)-1
          geoType = cellGeoType(crx(leftcut(2),iy,:),
     &                          cry(leftcut(2),iy,:))
          if (isRealCell(cflag(leftcut(2),iy,CELLFLAG_TYPE)) .and.
     &     cflag(leftcut(2), iy, CELLFLAG_LEFTFACE) == GRID_UNDEFINED
     &     .and. geoType /= CGEO_TRIA_NOLEFT) region(leftcut(2),iy,1)=10
          if(iy.ge.bottomcut(1).and.iy.lt.topcut(1)) then
            if (isRealCell(cflag(rightcut(2),iy,CELLFLAG_TYPE)))
     >       region(rightcut(2),iy,1)=11
          else
            if (isRealCell(cflag(rightcut(2),iy,CELLFLAG_TYPE)))
     >       region(rightcut(2),iy,1)=13
          endif
        enddo
      endif
! We find the stellarator island periodicity line: x-region 7
      if(nncut.eq.1.and.periodic_bc.eq.1) then
        do iy=iyt+1,ny
          if (isRealCell(cflag(inselix1,iy,CELLFLAG_TYPE)))
     >     region(inselix1,iy,1)=7
        enddo
      endif
      if(nncut.eq.1.and.nnreg(0).eq.4) then    ! SN case
        nnreg(1)=6
      elseif(nncut.eq.1.and.nnreg(0).eq.5) then    ! Stellarator island case
        nnreg(1)=7
      elseif(nncut.eq.2.and.topcut(1).eq.topcut(2)) then   ! DNC case
        nnreg(1)=12
      elseif(nncut.eq.2.and.topcut(1).ne.topcut(2)) then   ! DND case
        nnreg(1)=13
      endif
! Y-flux component
      if(periodic_bc.eq.1) then
        iyt=topcut(1)
      else
        iyt=ny
      endif
! We find the contact to the left (bottom) PFR wall: y-region 1
      do ix = -1, leftcut(1)-1
        iy = -1
        do while (iy.lt.topcut(1)
     &   .and. cflag(ix,iy,CELLFLAG_TYPE) /= GRID_INTERNAL)
          if (cflag(ix, iy, CELLFLAG_TYPE) == GRID_BOUNDARY) then
            geoType = cellGeoType(crx(ix,iy,:), cry(ix,iy,:))
            if (cflag(ix, iy, CELLFLAG_BOTTOMFACE) /= lefttargetindex(1)
     &       .and.
     &       cflag(ix, iy, CELLFLAG_BOTTOMFACE) /= righttargetindex(1)
     &       .and.
     &       cflag(ix, iy, CELLFLAG_BOTTOMFACE) /= GRID_UNDEFINED
     &       .and. geoType /= CGEO_TRIA_NOBOT) region(ix,iy,2)=1
            if (cflag(ix, iy, CELLFLAG_RIGHTFACE) /= lefttargetindex(1)
     &       .and.
     &       cflag(ix, iy, CELLFLAG_RIGHTFACE) /= righttargetindex(1)
     &       .and.
     &       cflag(ix, iy, CELLFLAG_RIGHTFACE) /= GRID_UNDEFINED
     &       .and. geoType /= CGEO_TRIA_NORIGHT)
     &       region(rightix(ix,iy),rightiy(ix,iy),1)=-1
            if (cflag(ix, iy, CELLFLAG_LEFTFACE) /= lefttargetindex(1)
     &       .and.
     &       cflag(ix, iy, CELLFLAG_LEFTFACE) /= righttargetindex(1)
     &       .and.
     &       cflag(ix, iy, CELLFLAG_LEFTFACE) /= GRID_UNDEFINED
     &       .and. geoType /= CGEO_TRIA_NOLEFT) region(ix,iy,1)=-1
          end if
          iy = iy + 1
        end do
      end do
      if (periodic_bc.ne.1) then
! We find the contact to the left (bottom) divertor wall: y-region 5
        do ix = -1, leftcut(1)-1
          iy = ny
          do while (iy.ge.topcut(1)
     &     .and. cflag(ix,iy,CELLFLAG_TYPE) /= GRID_INTERNAL)
            if (cflag(ix, iy, CELLFLAG_TYPE) == GRID_BOUNDARY) then
              geoType = cellGeoType(crx(ix,iy,:), cry(ix,iy,:))
              if (cflag(ix, iy, CELLFLAG_TOPFACE) /= lefttargetindex(1)
     &         .and.
     &         cflag(ix, iy, CELLFLAG_TOPFACE) /= righttargetindex(1)
     &         .and.
     &         cflag(ix, iy, CELLFLAG_TOPFACE) /= GRID_UNDEFINED
     &         .and. geoType /= CGEO_TRIA_NOTOP)
     &         region(topix(ix,iy),topiy(ix,iy),2)=5
              if (cflag(ix, iy, CELLFLAG_RIGHTFACE) /=
     &                                              lefttargetindex(1)
     &         .and.
     &         cflag(ix, iy, CELLFLAG_RIGHTFACE) /= righttargetindex(1)
     &         .and.
     &         cflag(ix, iy, CELLFLAG_RIGHTFACE) /= GRID_UNDEFINED
     &         .and. geoType /= CGEO_TRIA_NORIGHT)
     &         region(rightix(ix,iy),rightiy(ix,iy),1)=-5
              if (cflag(ix, iy, CELLFLAG_LEFTFACE) /= lefttargetindex(1)
     &         .and.
     &         cflag(ix, iy, CELLFLAG_LEFTFACE) /= righttargetindex(1)
     &         .and.
     &         cflag(ix, iy, CELLFLAG_LEFTFACE) /= GRID_UNDEFINED
     &         .and. geoType /= CGEO_TRIA_NOLEFT) region(ix,iy,1)=-5
            end if
            iy = iy - 1
          end do
        end do
      elseif (periodic_bc.eq.1) then
! Stellarator island: We find the entrance to the left PFR: y-region 5
        do ix=-1,leftcut(1)-1
          if (isRealCell(cflag(ix,iyt,CELLFLAG_TYPE)))
     >     region(ix,iyt,2)=5
        enddo
      endif
      if(nncut.lt.2) then  ! SN case
        do ix=leftcut(1),rightcut(1)-1
! SN: We find the core boundary: y-region 2
          iy = -1
          do while (iy.le.topcut(1) .and.
     &     .not.isRealCell(cflag(ix, iy, CELLFLAG_TYPE)))
            iy = iy+1
          end do
          if (iy.le.topcut(1)) then
            region(ix,iy,2)=2
          endif
! SN: We find the separatrix: y-region 4
          region(ix,topcut(1),2)=4
        enddo
        if(periodic_bc.eq.1) then
          do ix=inselix1,inselix2
! We find the island inner boundary: y-region 6
            iy = ny
            do while (iy.ge.inseliy .and.
     &       .not.isRealCell(cflag(ix, iy, CELLFLAG_TYPE)))
              iy = iy-1
            end do
            if (iy.ge.inseliy) then
              region(topix(ix,iy),topiy(ix,iy),2)=6
            endif
! We find the island separatrix: y-region 8
            region(ix,inseliy,2)=8
          enddo
        else
! SN: We find the contact to the main chamber wall: y-region 6
          do ix = leftcut(1), rightcut(1)-1
            iy = ny
            do while (iy.ge.topcut(1)
     &       .and. cflag(ix,iy,CELLFLAG_TYPE) /= GRID_INTERNAL)
              if (cflag(ix, iy, CELLFLAG_TYPE) == GRID_BOUNDARY) then
                geoType = cellGeoType(crx(ix,iy,:), cry(ix,iy,:))
                if (cflag(ix, iy, CELLFLAG_TOPFACE) /=
     &                                              lefttargetindex(1)
     &           .and.
     &           cflag(ix, iy, CELLFLAG_TOPFACE) /= righttargetindex(1)
     &           .and.
     &           cflag(ix, iy, CELLFLAG_TOPFACE) /= GRID_UNDEFINED
     &           .and. geoType /= CGEO_TRIA_NOTOP)
     &           region(topix(ix,iy),topiy(ix,iy),2)=6
                if (cflag(ix, iy, CELLFLAG_RIGHTFACE) /=
     &                                                lefttargetindex(1)
     &           .and.
     &           cflag(ix, iy, CELLFLAG_RIGHTFACE) /=
     &                                               righttargetindex(1)
     &           .and.
     &           cflag(ix, iy, CELLFLAG_RIGHTFACE) /= GRID_UNDEFINED
     &           .and. geoType /= CGEO_TRIA_NORIGHT)
     &           region(rightix(ix,iy),rightiy(ix,iy),1)=-6
                if (cflag(ix, iy, CELLFLAG_LEFTFACE) /=
     &                                               lefttargetindex(1)
     &           .and.
     &           cflag(ix, iy, CELLFLAG_LEFTFACE) /= righttargetindex(1)
     &           .and.
     &           cflag(ix, iy, CELLFLAG_LEFTFACE) /= GRID_UNDEFINED
     &           .and. geoType /= CGEO_TRIA_NOLEFT) region(ix,iy,1)=-6
              end if
              iy = iy - 1
            end do
          end do
        endif
! SN: We find the contact to the right PFR wall: y-region 3
        do ix = rightcut(1), nx
          iy = -1
          do while (iy.le.topcut(1)
     &     .and. cflag(ix,iy,CELLFLAG_TYPE) /= GRID_INTERNAL)
            if (cflag(ix, iy, CELLFLAG_TYPE) == GRID_BOUNDARY) then
              geoType = cellGeoType(crx(ix,iy,:), cry(ix,iy,:))
              if (cflag(ix, iy, CELLFLAG_BOTTOMFACE) /=
     &                                               lefttargetindex(1)
     &         .and.
     &         cflag(ix, iy, CELLFLAG_BOTTOMFACE) /= righttargetindex(1)
     &         .and.
     &         cflag(ix, iy, CELLFLAG_BOTTOMFACE) /= GRID_UNDEFINED
     &         .and. geoType /= CGEO_TRIA_NOBOT) region(ix,iy,2)=3
              if (cflag(ix, iy, CELLFLAG_RIGHTFACE) /=
     &                                              lefttargetindex(1)
     &         .and.
     &         cflag(ix, iy, CELLFLAG_RIGHTFACE) /= righttargetindex(1)
     &         .and.
     &         cflag(ix, iy, CELLFLAG_RIGHTFACE) /= GRID_UNDEFINED
     &         .and. geoType /= CGEO_TRIA_NORIGHT)
     &         region(rightix(ix,iy),rightiy(ix,iy),1)=-3
              if (cflag(ix, iy, CELLFLAG_LEFTFACE) /=
     &                                             lefttargetindex(1)
     &         .and.
     &         cflag(ix, iy, CELLFLAG_LEFTFACE) /= righttargetindex(1)
     &         .and.
     &         cflag(ix, iy, CELLFLAG_LEFTFACE) /= GRID_UNDEFINED
     &         .and. geoType /= CGEO_TRIA_NOLEFT) region(ix,iy,1)=-3
            end if
            iy = iy + 1
          end do
! SN: We find the contact to the right divertor wall: y-region 7
          iy = ny
          do while (iy.gt.topcut(1)
     &     .and. cflag(ix,iy,CELLFLAG_TYPE) /= GRID_INTERNAL)
            if (cflag(ix, iy, CELLFLAG_TYPE) == GRID_BOUNDARY) then
              geoType = cellGeoType(crx(ix,iy,:), cry(ix,iy,:))
              if (cflag(ix, iy, CELLFLAG_TOPFACE) /= lefttargetindex(1)
     &         .and.
     &         cflag(ix, iy, CELLFLAG_TOPFACE) /= righttargetindex(1)
     &         .and.
     &         cflag(ix, iy, CELLFLAG_TOPFACE) /= GRID_UNDEFINED
     &         .and. geoType /= CGEO_TRIA_NOTOP)
     &         region(topix(ix,iy),topiy(ix,iy),2)=7
              if (cflag(ix, iy, CELLFLAG_RIGHTFACE) /=
     &                                              lefttargetindex(1)
     &         .and.
     &         cflag(ix, iy, CELLFLAG_RIGHTFACE) /= righttargetindex(1)
     &         .and.
     &         cflag(ix, iy, CELLFLAG_RIGHTFACE) /= GRID_UNDEFINED
     &         .and. geoType /= CGEO_TRIA_NORIGHT)
     &         region(rightix(ix,iy),rightiy(ix,iy),1)=-7
              if (cflag(ix, iy, CELLFLAG_LEFTFACE) /= lefttargetindex(1)
     &         .and.
     &         cflag(ix, iy, CELLFLAG_LEFTFACE) /= righttargetindex(1)
     &         .and.
     &         cflag(ix, iy, CELLFLAG_LEFTFACE) /= GRID_UNDEFINED
     &         .and. geoType /= CGEO_TRIA_NOLEFT) region(ix,iy,1)=-7
            end if
            iy = iy - 1
          end do
        end do
      elseif(nncut.eq.2) then  !DN case
        do ix=leftcut(1),leftcut(2)-1
! DN: We find the inner core boundary: y-region 2
          iy = -1
          do while (iy.le.topcut(1) .and. iy.le.topcut(2) .and.
     &     .not.isRealCell(cflag(ix, iy, CELLFLAG_TYPE)))
            iy = iy+1
          end do
          if (iy.le.topcut(1) .and. iy.le.topcut(2)) then
            region(ix,iy,2)=2
          endif
! DN: We find the inner separatrix: y-region 4
          region(ix,min(topcut(1),topcut(2)),2)=4
! DN: We find the contact to the inner main chamber wall: y-region 6
          iy = ny
          do while (iy.ge.topcut(1) .and. iy.ge.topcut(2) .and.
     &              cflag(ix,iy,CELLFLAG_TYPE) /= GRID_INTERNAL)
            if (cflag(ix, iy, CELLFLAG_TYPE) == GRID_BOUNDARY) then
              geoType = cellGeoType(crx(ix,iy,:), cry(ix,iy,:))
              if (cflag(ix, iy, CELLFLAG_TOPFACE) /= lefttargetindex(1)
     &         .and.
     &         cflag(ix, iy, CELLFLAG_TOPFACE) /= righttargetindex(1)
     &         .and.
     &         cflag(ix, iy, CELLFLAG_TOPFACE) /= lefttargetindex(2)
     &         .and.
     &         cflag(ix, iy, CELLFLAG_TOPFACE) /= righttargetindex(2)
     &         .and.
     &         cflag(ix, iy, CELLFLAG_TOPFACE) /= GRID_UNDEFINED
     &         .and. geoType /= CGEO_TRIA_NOTOP)
     &         region(topix(ix,iy),topiy(ix,iy),2)=6
              if (cflag(ix, iy, CELLFLAG_RIGHTFACE) /=
     &                                              lefttargetindex(1)
     &         .and.
     &         cflag(ix, iy, CELLFLAG_RIGHTFACE) /= righttargetindex(1)
     &         .and.
     &         cflag(ix, iy, CELLFLAG_RIGHTFACE) /= lefttargetindex(2)
     &         .and.
     &         cflag(ix, iy, CELLFLAG_RIGHTFACE) /= righttargetindex(2)
     &         .and.
     &         cflag(ix, iy, CELLFLAG_RIGHTFACE) /= GRID_UNDEFINED
     &         .and. geoType /= CGEO_TRIA_NORIGHT)
     &         region(rightix(ix,iy),rightiy(ix,iy),1)=-6
              if (cflag(ix, iy, CELLFLAG_LEFTFACE) /= lefttargetindex(1)
     &         .and.
     &         cflag(ix, iy, CELLFLAG_LEFTFACE) /= righttargetindex(1)
     &         .and.
     &         cflag(ix, iy, CELLFLAG_LEFTFACE) /= lefttargetindex(2)
     &         .and.
     &         cflag(ix, iy, CELLFLAG_LEFTFACE) /= righttargetindex(2)
     &         .and.
     &         cflag(ix, iy, CELLFLAG_LEFTFACE) /= GRID_UNDEFINED
     &         .and. geoType /= CGEO_TRIA_NOLEFT) region(ix,iy,1)=-6
            end if
            iy = iy - 1
          end do
        end do
! DN: We find the contact to the top left PFR wall: y-region 3
        do ix = leftcut(2), rightcut(2)-1
          iy = -1
          do while (iy.lt.topcut(2)
     &     .and. cflag(ix,iy,CELLFLAG_TYPE) /= GRID_INTERNAL)
            if (cflag(ix, iy, CELLFLAG_TYPE) == GRID_BOUNDARY) then
              geoType = cellGeoType(crx(ix,iy,:), cry(ix,iy,:))
              if (cflag(ix, iy, CELLFLAG_BOTTOMFACE) /=
     &                                               lefttargetindex(2)
     &         .and.
     &         cflag(ix, iy, CELLFLAG_BOTTOMFACE) /= righttargetindex(2)
     &         .and.
     &         cflag(ix, iy, CELLFLAG_BOTTOMFACE) /= GRID_UNDEFINED
     &         .and. geoType /= CGEO_TRIA_NOBOT) region(ix,iy,2)=3
              if (cflag(ix, iy, CELLFLAG_RIGHTFACE) /=
     &                                              lefttargetindex(2)
     &         .and.
     &         cflag(ix, iy, CELLFLAG_RIGHTFACE) /= righttargetindex(2)
     &         .and.
     &         cflag(ix, iy, CELLFLAG_RIGHTFACE) /= GRID_UNDEFINED
     &         .and. geoType /= CGEO_TRIA_NORIGHT)
     &         region(rightix(ix,iy),rightiy(ix,iy),1)=-3
              if (cflag(ix, iy, CELLFLAG_LEFTFACE) /= lefttargetindex(2)
     &         .and.
     &         cflag(ix, iy, CELLFLAG_LEFTFACE) /= righttargetindex(2)
     &         .and.
     &         cflag(ix, iy, CELLFLAG_LEFTFACE) /= GRID_UNDEFINED
     &         .and. geoType /= CGEO_TRIA_NOLEFT) region(ix,iy,1)=-3
            end if
            iy = iy + 1
          end do
          if (iy.eq.topcut(2)) exit
        end do
! DN: We find the contact to the top left divertor wall: y-region 7
        do ix = leftcut(2), rightcut(2)-1
          iy = ny
          do while (iy.ge.topcut(2)
     &     .and. cflag(ix,iy,CELLFLAG_TYPE) /= GRID_INTERNAL)
            if (cflag(ix, iy, CELLFLAG_TYPE) == GRID_BOUNDARY) then
              geoType = cellGeoType(crx(ix,iy,:), cry(ix,iy,:))
              if (cflag(ix, iy, CELLFLAG_TOPFACE) /= lefttargetindex(2)
     &         .and.
     &         cflag(ix, iy, CELLFLAG_TOPFACE) /= righttargetindex(2)
     &         .and.
     &         cflag(ix, iy, CELLFLAG_TOPFACE) /= GRID_UNDEFINED
     &         .and. geoType /= CGEO_TRIA_NOTOP)
     &         region(topix(ix,iy),topiy(ix,iy),2)=7
              if (cflag(ix, iy, CELLFLAG_RIGHTFACE) /=
     &                                              lefttargetindex(2)
     &         .and.
     &         cflag(ix, iy, CELLFLAG_RIGHTFACE) /= righttargetindex(2)
     &         .and.
     &         cflag(ix, iy, CELLFLAG_RIGHTFACE) /= GRID_UNDEFINED
     &         .and. geoType /= CGEO_TRIA_NORIGHT)
     &         region(rightix(ix,iy),rightiy(ix,iy),1)=-7
              if (cflag(ix, iy, CELLFLAG_LEFTFACE) /= lefttargetindex(2)
     &         .and.
     &         cflag(ix, iy, CELLFLAG_LEFTFACE) /= righttargetindex(2)
     &         .and.
     &         cflag(ix, iy, CELLFLAG_LEFTFACE) /= GRID_UNDEFINED
     &         .and. geoType /= CGEO_TRIA_NOLEFT) region(ix,iy,1)=-7
            end if
            iy = iy - 1
          end do
        end do
! DN: We find the contact to the right top PFR wall: y-region 8
        do ix = rightcut(2)-1, leftcut(2), -1
          iy = -1
          do while (iy.lt.topcut(2)
     &     .and. cflag(ix,iy,CELLFLAG_TYPE) /= GRID_INTERNAL)
            if (cflag(ix, iy, CELLFLAG_TYPE) == GRID_BOUNDARY) then
              geoType = cellGeoType(crx(ix,iy,:), cry(ix,iy,:))
              if (cflag(ix, iy, CELLFLAG_BOTTOMFACE) /=
     &                                               lefttargetindex(2)
     &         .and.
     &         cflag(ix, iy, CELLFLAG_BOTTOMFACE) /= righttargetindex(2)
     &         .and.
     &         cflag(ix, iy, CELLFLAG_BOTTOMFACE) /= GRID_UNDEFINED
     &         .and. geoType /= CGEO_TRIA_NOBOT) region(ix,iy,2)=8
              if (cflag(ix, iy, CELLFLAG_RIGHTFACE) /=
     &                                              lefttargetindex(2)
     &         .and.
     &         cflag(ix, iy, CELLFLAG_RIGHTFACE) /= righttargetindex(2)
     &         .and.
     &         cflag(ix, iy, CELLFLAG_RIGHTFACE) /= GRID_UNDEFINED
     &         .and. geoType /= CGEO_TRIA_NORIGHT)
     &         region(rightix(ix,iy),rightiy(ix,iy),1)=-8
              if (cflag(ix, iy, CELLFLAG_LEFTFACE) /= lefttargetindex(2)
     &         .and.
     &         cflag(ix, iy, CELLFLAG_LEFTFACE) /= righttargetindex(2)
     &         .and.
     &         cflag(ix, iy, CELLFLAG_LEFTFACE) /= GRID_UNDEFINED
     &         .and. geoType /= CGEO_TRIA_NOLEFT) region(ix,iy,1)=-8
            end if
            iy = iy + 1
          end do
! We find the contact to the right top divertor wall: y-region 12
          iy = ny
          do while (iy.ge.topcut(2)
     &     .and. cflag(ix,iy,CELLFLAG_TYPE) /= GRID_INTERNAL)
            if (cflag(ix, iy, CELLFLAG_TYPE) == GRID_BOUNDARY) then
              geoType = cellGeoType(crx(ix,iy,:), cry(ix,iy,:))
              if (cflag(ix, iy, CELLFLAG_TOPFACE) /= lefttargetindex(2)
     &         .and.
     &         cflag(ix, iy, CELLFLAG_TOPFACE) /= righttargetindex(2)
     &         .and.
     &         cflag(ix, iy, CELLFLAG_TOPFACE) /= GRID_UNDEFINED
     &         .and. geoType /= CGEO_TRIA_NOTOP)
     &         region(topix(ix,iy),topiy(ix,iy),2)=12
              if (cflag(ix, iy, CELLFLAG_RIGHTFACE) /=
     &                                              lefttargetindex(2)
     &         .and.
     &         cflag(ix, iy, CELLFLAG_RIGHTFACE) /= righttargetindex(2)
     &         .and.
     &         cflag(ix, iy, CELLFLAG_RIGHTFACE) /= GRID_UNDEFINED
     &         .and. geoType /= CGEO_TRIA_NORIGHT)
     &         region(rightix(ix,iy),rightiy(ix,iy),1)=-12
              if (cflag(ix, iy, CELLFLAG_LEFTFACE) /= lefttargetindex(2)
     &         .and.
     &         cflag(ix, iy, CELLFLAG_LEFTFACE) /= righttargetindex(2)
     &         .and.
     &         cflag(ix, iy, CELLFLAG_LEFTFACE) /= GRID_UNDEFINED
     &         .and. geoType /= CGEO_TRIA_NOLEFT) region(ix,iy,1)=-12
            end if
            iy = iy - 1
          end do
        end do
        do ix = rightcut(2),rightcut(1)-1
! DN: We find the outer core boundary: y-region 9
          iy = -1
          do while (iy.le.topcut(1) .and. iy.le.topcut(2) .and.
     &     .not.isRealCell(cflag(ix, iy, CELLFLAG_TYPE)))
            iy = iy+1
          end do
          if (iy.le.topcut(1) .and. iy.le.topcut(2)) then
            region(ix,iy,2)=9
          endif
! DN: We find the outer separatrix: y-region 11
          region(ix,min(topcut(1),topcut(2)),2)=11
! DN: We find the contact to the inner main chamber wall: y-region 13
          iy = ny
          do while (iy.ge.topcut(1) .and. iy.ge.topcut(2) .and.
     &              cflag(ix,iy,CELLFLAG_TYPE) /= GRID_INTERNAL)
            if (cflag(ix, iy, CELLFLAG_TYPE) == GRID_BOUNDARY) then
              geoType = cellGeoType(crx(ix,iy,:), cry(ix,iy,:))
              if (cflag(ix, iy, CELLFLAG_TOPFACE) /= lefttargetindex(1)
     &         .and.
     &         cflag(ix, iy, CELLFLAG_TOPFACE) /= righttargetindex(1)
     &         .and.
     &         cflag(ix, iy, CELLFLAG_TOPFACE) /= lefttargetindex(2)
     &         .and.
     &         cflag(ix, iy, CELLFLAG_TOPFACE) /= righttargetindex(2)
     &         .and.
     &         cflag(ix, iy, CELLFLAG_TOPFACE) /= GRID_UNDEFINED
     &         .and. geoType /= CGEO_TRIA_NOTOP)
     &         region(topix(ix,iy),topiy(ix,iy),2)=13
              if (cflag(ix, iy, CELLFLAG_RIGHTFACE) /=
     &                                              lefttargetindex(1)
     &         .and.
     &         cflag(ix, iy, CELLFLAG_RIGHTFACE) /= righttargetindex(1)
     &         .and.
     &         cflag(ix, iy, CELLFLAG_RIGHTFACE) /= lefttargetindex(2)
     &         .and.
     &         cflag(ix, iy, CELLFLAG_RIGHTFACE) /= righttargetindex(2)
     &         .and.
     &         cflag(ix, iy, CELLFLAG_RIGHTFACE) /= GRID_UNDEFINED
     &         .and. geoType /= CGEO_TRIA_NORIGHT)
     &         region(rightix(ix,iy),rightiy(ix,iy),1)=-13
              if (cflag(ix, iy, CELLFLAG_LEFTFACE) /= lefttargetindex(1)
     &         .and.
     &         cflag(ix, iy, CELLFLAG_LEFTFACE) /= righttargetindex(1)
     &         .and.
     &         cflag(ix, iy, CELLFLAG_LEFTFACE) /= lefttargetindex(2)
     &         .and.
     &         cflag(ix, iy, CELLFLAG_LEFTFACE) /= righttargetindex(2)
     &         .and.
     &         cflag(ix, iy, CELLFLAG_LEFTFACE) /= GRID_UNDEFINED
     &         .and. geoType /= CGEO_TRIA_NOLEFT) region(ix,iy,1)=-13
            end if
            iy = iy - 1
          end do
        end do
! DN: We find the contact to the bottom right PFR wall: y-region 10
        do ix = rightcut(1), nx
          iy = -1
          do while (iy.le.topcut(1)
     &     .and. cflag(ix,iy,CELLFLAG_TYPE) /= GRID_INTERNAL)
            if (cflag(ix, iy, CELLFLAG_TYPE) == GRID_BOUNDARY) then
              geoType = cellGeoType(crx(ix,iy,:), cry(ix,iy,:))
              if (cflag(ix, iy, CELLFLAG_BOTTOMFACE) /=
     &                                               lefttargetindex(1)
     &         .and.
     &         cflag(ix, iy, CELLFLAG_BOTTOMFACE) /= righttargetindex(1)
     &         .and.
     &         cflag(ix, iy, CELLFLAG_BOTTOMFACE) /= GRID_UNDEFINED
     &         .and. geoType /= CGEO_TRIA_NOBOT) region(ix,iy,2)=10
              if (cflag(ix, iy, CELLFLAG_RIGHTFACE) /=
     &                                              lefttargetindex(1)
     &         .and.
     &         cflag(ix, iy, CELLFLAG_RIGHTFACE) /= righttargetindex(1)
     &         .and.
     &         cflag(ix, iy, CELLFLAG_RIGHTFACE) /= GRID_UNDEFINED
     &         .and. geoType /= CGEO_TRIA_NORIGHT)
     &         region(rightix(ix,iy),rightiy(ix,iy),1)=-10
              if (cflag(ix, iy, CELLFLAG_LEFTFACE) /=
     &                                             lefttargetindex(1)
     &         .and.
     &         cflag(ix, iy, CELLFLAG_LEFTFACE) /= righttargetindex(1)
     &         .and.
     &         cflag(ix, iy, CELLFLAG_LEFTFACE) /= GRID_UNDEFINED
     &         .and. geoType /= CGEO_TRIA_NOLEFT) region(ix,iy,1)=-10
            end if
            iy = iy + 1
          end do
! DN: We find the contact to the right divertor wall: y-region 14
          iy = ny
          do while (iy.gt.topcut(1)
     &     .and. cflag(ix,iy,CELLFLAG_TYPE) /= GRID_INTERNAL)
            if (cflag(ix, iy, CELLFLAG_TYPE) == GRID_BOUNDARY) then
              geoType = cellGeoType(crx(ix,iy,:), cry(ix,iy,:))
              if (cflag(ix, iy, CELLFLAG_TOPFACE) /= lefttargetindex(1)
     &             .and.
     &            cflag(ix, iy, CELLFLAG_TOPFACE) /= righttargetindex(1)
     &             .and.
     &            cflag(ix, iy, CELLFLAG_TOPFACE) /= GRID_UNDEFINED
     &             .and. geoType /= CGEO_TRIA_NOTOP)
     &         region(topix(ix,iy),topiy(ix,iy),2)=14
              if (cflag(ix, iy, CELLFLAG_RIGHTFACE) /=
     &                                              lefttargetindex(1)
     &         .and.
     &         cflag(ix, iy, CELLFLAG_RIGHTFACE) /= righttargetindex(1)
     &         .and.
     &         cflag(ix, iy, CELLFLAG_RIGHTFACE) /= GRID_UNDEFINED
     &         .and. geoType /= CGEO_TRIA_NORIGHT)
     &         region(rightix(ix,iy),rightiy(ix,iy),1)=-14
              if (cflag(ix, iy, CELLFLAG_LEFTFACE) /= lefttargetindex(1)
     &         .and.
     &         cflag(ix, iy, CELLFLAG_LEFTFACE) /= righttargetindex(1)
     &         .and.
     &         cflag(ix, iy, CELLFLAG_LEFTFACE) /= GRID_UNDEFINED
     &         .and. geoType /= CGEO_TRIA_NOLEFT) region(ix,iy,1)=-14
              end if
              iy = iy - 1
            end do
          end do
        endif
        if(periodic_bc.eq.1) then
          nnreg(2)=8
        elseif(nncut.lt.2) then
          nnreg(2)=7
        elseif(nncut.eq.2) then
          nnreg(2)=14
        endif
      else ! 1-D and 2-D slab cases
! We identify the structure index of the left target
        iy = ny/2
        ix = -1
        lefttargetindex(1)=GRID_UNDEFINED
        do while (ix.le.nx .and. lefttargetindex(1).eq.GRID_UNDEFINED)
          if (cflag(ix, iy, CELLFLAG_TYPE) == GRID_BOUNDARY) then
            geoType = cellGeoType(crx(ix,iy,:), cry(ix,iy,:))
            if (geoType /= CGEO_TRIA_NOLEFT .and.
     &       cflag(ix, iy, CELLFLAG_LEFTFACE) /= GRID_UNDEFINED) then
              lefttargetindex(1) = cflag(ix,iy,CELLFLAG_LEFTFACE)
            else if (geoType /= CGEO_TRIA_NOTOP .and.
     &       cflag(ix, iy, CELLFLAG_TOPFACE) /= GRID_UNDEFINED) then
              lefttargetindex(1) = cflag(ix,iy,CELLFLAG_TOPFACE)
            else if (geoType /= CGEO_TRIA_NOBOT .and.
     &       cflag(ix, iy, CELLFLAG_BOTTOMFACE) /= GRID_UNDEFINED) then
              lefttargetindex(1) = cflag(ix,iy,CELLFLAG_BOTTOMFACE)
            else
#ifdef BUILDING_CARRE
              stop
     >  'Invalid cflags settings (L), should not have gotten here!'
#else
              call xerrab
     . ('Invalid cflags settings (L), should not have gotten here!')
#endif
            endif
          else
            ix = ix+1
          endif
        end do
#ifdef BUILDING_CARRE
        if (.not. lefttargetindex(1).ne.GRID_UNDEFINED)
     >   stop 'Left target not found!'
#else
        call xertst(lefttargetindex(1).ne.GRID_UNDEFINED,
     .        'Left target not found!')
#endif
! We identify the structure index of the right target
        ix = nx
        righttargetindex(1)=GRID_UNDEFINED
        do while (ix.ge.-1 .and. righttargetindex(1).eq.GRID_UNDEFINED)
          if (cflag(ix, iy, CELLFLAG_TYPE) == GRID_BOUNDARY) then
            geoType = cellGeoType(crx(ix,iy,:), cry(ix,iy,:))
            if (geoType /= CGEO_TRIA_NORIGHT .and.
     &       cflag(ix, iy, CELLFLAG_RIGHTFACE) /= GRID_UNDEFINED) then
              righttargetindex(1) = cflag(ix,iy,CELLFLAG_RIGHTFACE)
            else if (geoType /= CGEO_TRIA_NOTOP .and.
     &       cflag(ix, iy, CELLFLAG_TOPFACE) /= GRID_UNDEFINED) then
              righttargetindex(1) = cflag(ix,iy,CELLFLAG_TOPFACE)
            else if (geoType /= CGEO_TRIA_NOBOT .and.
     &       cflag(ix, iy, CELLFLAG_BOTTOMFACE) /= GRID_UNDEFINED) then
              righttargetindex(1) = cflag(ix,iy,CELLFLAG_BOTTOMFACE)
            else
#ifdef BUILDING_CARRE
              stop
     >   'Invalid cflags settings (R), should not have gotten here!'
#else
              call xerrab
     .  ('Invalid cflags settings (R), should not have gotten here!')
#endif
            endif
          else
            ix = ix-1
          endif
        end do
#ifdef BUILDING_CARRE
        if (.not. righttargetindex(1).ne.GRID_UNDEFINED)
     >   stop 'Right target not found!'
#else
        call xertst(righttargetindex(1).ne.GRID_UNDEFINED,
     .        'Right target not found!')
#endif
        do iy = -1, ny
! We identify the left edge
          ix = -1
          do while (ix.lt.nx
     &     .and. cflag(ix, iy, CELLFLAG_TYPE) /= GRID_INTERNAL)
            if (cflag(ix, iy, CELLFLAG_TYPE) == GRID_BOUNDARY) then
              geoType = cellGeoType(crx(ix,iy,:), cry(ix,iy,:))
              if (cflag(ix, iy, CELLFLAG_LEFTFACE) == lefttargetindex(1)
     &         .and. geoType /= CGEO_TRIA_NOLEFT) region(ix,iy,1)=1
              if (cflag(ix, iy, CELLFLAG_TOPFACE) == lefttargetindex(1)
     &         .and. geoType /= CGEO_TRIA_NOTOP)
     &         region(topix(ix,iy),topiy(ix,iy),2)=-1
              if (cflag(ix, iy, CELLFLAG_BOTTOMFACE) ==
     &                                               lefttargetindex(1)
     &         .and. geoType /= CGEO_TRIA_NOBOT) region(ix,iy,2)=-1
            endif
            ix = ix + 1
          enddo
! We identify the right edge
          ix = nx
          do while (ix.gt.-1
     &     .and. cflag(ix, iy, CELLFLAG_TYPE) /= GRID_INTERNAL)
            if (cflag(ix, iy, CELLFLAG_TYPE) == GRID_BOUNDARY) then
              geoType = cellGeoType(crx(ix,iy,:), cry(ix,iy,:))
              if (cflag(ix, iy, CELLFLAG_RIGHTFACE) ==
     &                                              righttargetindex(1)
     &         .and. geoType /= CGEO_TRIA_NORIGHT)
     &         region(rightix(ix,iy),rightiy(ix,iy),1)=2
              if (cflag(ix, iy, CELLFLAG_TOPFACE) == righttargetindex(1)
     &          .and. geoType /= CGEO_TRIA_NOTOP)
     &         region(topix(ix,iy),topiy(ix,iy),2)=-2
              if (cflag(ix, iy, CELLFLAG_BOTTOMFACE) ==
     &                                               righttargetindex(1)
     &         .and. geoType /= CGEO_TRIA_NOBOT) region(ix,iy,2)=-2
            endif
            ix = ix - 1
          enddo
        enddo
        nnreg(1)=2
! We identify the bottom edge
        do ix = -1, nx
          iy = -1
          do while (iy.lt.ny
     &     .and. cflag(ix,iy,CELLFLAG_TYPE) /= GRID_INTERNAL)
            if (cflag(ix, iy, CELLFLAG_TYPE) == GRID_BOUNDARY) then
              geoType = cellGeoType(crx(ix,iy,:), cry(ix,iy,:))
              if (cflag(ix, iy, CELLFLAG_BOTTOMFACE) /= GRID_UNDEFINED
     &         .and.
     &         cflag(ix, iy, CELLFLAG_BOTTOMFACE) /= lefttargetindex(1)
     &         .and.
     &         cflag(ix, iy, CELLFLAG_BOTTOMFACE) /= righttargetindex(1)
     &         .and. geoType /= CGEO_TRIA_NOBOT) region(ix,iy,2)=1
              if (cflag(ix, iy, CELLFLAG_RIGHTFACE) /= GRID_UNDEFINED
     &         .and.
     &         cflag(ix, iy, CELLFLAG_RIGHTFACE) /= lefttargetindex(1)
     &         .and.
     &         cflag(ix, iy, CELLFLAG_RIGHTFACE) /= righttargetindex(1)
     &         .and. geoType /= CGEO_TRIA_NORIGHT)
     &         region(rightix(ix,iy),rightiy(ix,iy),1)=-1
              if (cflag(ix, iy, CELLFLAG_LEFTFACE) /= GRID_UNDEFINED
     &         .and.
     &         cflag(ix, iy, CELLFLAG_LEFTFACE) /= lefttargetindex(1)
     &         .and.
     &         cflag(ix, iy, CELLFLAG_LEFTFACE) /= righttargetindex(1)
     &         .and. geoType /= CGEO_TRIA_NOLEFT) region(ix,iy,1)=-1
            end if
            iy = iy + 1
          end do
        end do
! We identify the top edge
        do ix = -1, nx
          iy = ny
          do while (iy.gt.-1
     &     .and. cflag(ix,iy,CELLFLAG_TYPE) /= GRID_INTERNAL)
            if (cflag(ix, iy, CELLFLAG_TYPE) == GRID_BOUNDARY) then
              geoType = cellGeoType(crx(ix,iy,:), cry(ix,iy,:))
              if (cflag(ix, iy, CELLFLAG_TOPFACE) /= GRID_UNDEFINED
     &         .and.
     &         cflag(ix, iy, CELLFLAG_TOPFACE) /= lefttargetindex(1)
     &         .and.
     &         cflag(ix, iy, CELLFLAG_TOPFACE) /= righttargetindex(1)
     &         .and. geoType /= CGEO_TRIA_NOTOP)
     &         region(topix(ix,iy),topiy(ix,iy),2)=2
              if (cflag(ix, iy, CELLFLAG_RIGHTFACE) /= GRID_UNDEFINED
     &         .and.
     &         cflag(ix, iy, CELLFLAG_RIGHTFACE) /= lefttargetindex(1)
     &         .and.
     &         cflag(ix, iy, CELLFLAG_RIGHTFACE) /= righttargetindex(1)
     &         .and. geoType /= CGEO_TRIA_NORIGHT)
     &         region(rightix(ix,iy),rightiy(ix,iy),1)=-2
              if (cflag(ix, iy, CELLFLAG_LEFTFACE) /= GRID_UNDEFINED
     &         .and.
     &         cflag(ix, iy, CELLFLAG_LEFTFACE) /= lefttargetindex(1)
     &         .and.
     &         cflag(ix, iy, CELLFLAG_LEFTFACE) /= righttargetindex(1)
     &         .and. geoType /= CGEO_TRIA_NOLEFT) region(ix,iy,1)=-2
            end if
            iy = iy - 1
          end do
        end do
        nnreg(2)=2
      endif

    ! set up region where the residuals should be ignored (none yet)
      do iy=-1,ny
        do ix=-1,nx
            if(region(ix,iy,0).eq.0 .or.
     &       isGhostCell(cflag(ix,iy,CELLFLAG_TYPE))) then
                resignore(ix,iy,1)=0
                resignore(ix,iy,2)=0
            else
                resignore(ix,iy,1)=1
                resignore(ix,iy,2)=1
            endif
        enddo
      enddo

      return
      end subroutine init_region_extended

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

      function cr(ix,iy)
      implicit none
      integer, intent(in) :: ix, iy
      real(kind=R8) :: cr
      cr = 0.25_R8*(crx(ix,iy,0)+crx(ix,iy,1)+crx(ix,iy,2)+crx(ix,iy,3))
      return
      end function cr

      function cz(ix,iy)
      implicit none
      integer, intent(in) :: ix, iy
      real(kind=R8) :: cz
      cz = 0.25_R8*(cry(ix,iy,0)+cry(ix,iy,1)+cry(ix,iy,2)+cry(ix,iy,3))
      return
      end function cz

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
