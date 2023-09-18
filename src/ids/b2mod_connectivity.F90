!!-----------------------------------------------------------------------------
!! DOCUMENTATION (doxygen 1.8.8)::
!>      @section b2mod_conn_desc Description
!!      Module providing Basic connectivity routines and routines for obtaining
!!      cell and grid characterization information.
!!
!!----------------------------------------------------------------------------
module b2mod_connectivity

    use b2mod_types
    use carre_constants
    use b2mod_cellhelper
    use logging
    use helper

    implicit none
#include "DIMENSIONS.F"

    integer, parameter :: NO_CONNECTIVITY = huge(0) !< Constant to mark in
        !< connectivity arrays that no connectivity available

contains

  !> Check the connectivity for errors.
  subroutine test_connectivity(nx,ny,crx,cry,cflag, &
      & leftix,leftiy,rightix,rightiy, &
      & topix,topiy,bottomix,bottomiy)

    use b2mod_types
    implicit none

    !!   ..input arguments (unchanged on exit)
    integer, intent(in) :: nx, ny
    integer cflag(-1:nx,-1:ny,CARREOUT_NCELLFLAGS)
    real (kind=R8), intent(in) :: &
        & crx(-1:nx,-1:ny,0:3), cry(-1:nx,-1:ny,0:3)
    integer, intent(in) :: &
        & leftix(-1:nx,-1:ny),leftiy(-1:nx,-1:ny), &
        & rightix(-1:nx,-1:ny),rightiy(-1:nx,-1:ny), &
        & topix(-1:nx,-1:ny),topiy(-1:nx,-1:ny), &
        & bottomix(-1:nx,-1:ny),bottomiy(-1:nx,-1:ny)

    !! internal
    integer :: i, ix, iy, counter
    logical :: rightFace, leftFace, topFace, botFace !! has ... face
    logical :: rightNb, leftNb, topNb, botNb !! has ... neighbour
    logical :: error, thisCellError !! error occurred
    character*80 :: message

    error = .false.
    counter = 0

    do ix = -1, nx
        do iy = -1, ny

            if ( isUnusedCell( cflag(ix, iy, CELLFLAG_TYPE) ) ) cycle

            leftFace = .not. points_match(crx(ix,iy,0), cry(ix,iy,0), &
                 & crx(ix,iy,2), cry(ix,iy,2))
            botFace = .not. points_match(crx(ix,iy,0), cry(ix,iy,0), &
                 & crx(ix,iy,1), cry(ix,iy,1))
            rightFace = .not. points_match(crx(ix,iy,1), cry(ix,iy,1), &
                 & crx(ix,iy,3), cry(ix,iy,3))
            topFace = .not. points_match(crx(ix,iy,2), cry(ix,iy,2), &
                 & crx(ix,iy,3), cry(ix,iy,3))

            leftNb = isInDomain(nx,ny,leftix(ix,iy),leftiy(ix,iy))
            botNb = isInDomain(nx,ny,bottomix(ix,iy),bottomiy(ix,iy))
            rightNb = isInDomain(nx,ny,rightix(ix,iy),rightiy(ix,iy))
            topNb = isInDomain(nx,ny,topix(ix,iy),topiy(ix,iy))

            thisCellError = .false.

            select case ( cflag(ix, iy, CELLFLAG_TYPE) )
            case ( GRID_INTERNAL, GRID_BOUNDARY )
                if (leftFace .and. .not. leftNb) then
                    write (*,'(a,2i4,a)') "test_connectivity: ix, iy = ", ix, iy," : inner cell missing left neighbour"
                    thisCellError = .true.
                end if
                if (botFace .and. .not. botNb) then
                    write (*,'(a,2i4,a)') "test_connectivity: ix, iy = ", ix, iy," : inner cell missing bottom neighbour"
                    thisCellError = .true.
                end if
                if (rightFace .and. .not. rightNb) then
                    write (*,'(a,2i4,a)') "test_connectivity: ix, iy = ", ix, iy," : inner cell missing right neighbour"
                    thisCellError = .true.
                end if
                if (topFace .and. .not. topNb) then
                    write (*,'(a,2i4,a)') "test_connectivity: ix, iy = ", ix, iy," : inner cell missing top neighbour"
                    thisCellError = .true.
                end if
            case ( GRID_GUARD )

                if (.not. (leftNb .or. botNb .or. rightNb .or. topNb)) then
                    write (*,'(a,2i4,a)') "test_connectivity: ix, iy = ", ix, iy," : guard cell without neighbour"
                    thisCellError = .true.
                end if
                if ( cflag(ix, iy, CELLFLAG_LEFTFACE) /= GRID_UNDEFINED .and. .not. leftNb ) then
                    write (*,'(a,2i4,a)') "test_connectivity: ix, iy = ", ix, iy," : guard cell, left nb broken"
                    thisCellError = .true.
                end if
                if ( cflag(ix, iy, CELLFLAG_BOTTOMFACE) /= GRID_UNDEFINED .and. .not. botNb ) then
                    write (*,'(a,2i4,a)') "test_connectivity: ix, iy = ", ix, iy," : guard cell, bottom nb broken"
                    thisCellError = .true.
                end if
                if ( cflag(ix, iy, CELLFLAG_RIGHTFACE) /= GRID_UNDEFINED .and. .not. rightNb ) then
                    write (*,'(a,2i4,a)') "test_connectivity: ix, iy = ", ix, iy," : guard cell, right nb broken"
                    thisCellError = .true.
                end if
                if ( cflag(ix, iy, CELLFLAG_TOPFACE) /= GRID_UNDEFINED .and. .not. topNb ) then
                    write (*,'(a,2i4,a)') "test_connectivity: ix, iy = ", ix, iy," : guard cell, top nb broken"
                    thisCellError = .true.
                end if
            end select

            if (thisCellError) then
                write(*,'(a,4f12.6)') 'crx = ',(crx(ix,iy,i),i=0,3)
                write(*,'(a,4f12.6)') 'cry = ',(cry(ix,iy,i),i=0,3)
                counter = counter + 1
                error = .true.
            end if

        end do
    end do

    if (error) then
      if (counter.eq.1) then
        message = "test_connectivity: 1 error found"
      else
        write(message,'(a,i4,a)') "test_connectivity: ",counter," errors found"
      end if
#ifdef BUILDING_CARRE
      stop "test_connectivity: error(s) found"
#else
      call xerrab ( message )
#endif
    end if

  return
  end subroutine test_connectivity



  subroutine init_region(nx,ny,nncut,nncutmax, &
      & leftcut,rightcut,topcut,bottomcut,ccut, &
      & leftix,leftiy,rightix,rightiy,topix,topiy,bottomix,bottomiy, &
      & region,nnreg,resignore, &
      & crx,cry,periodic_bc,cflag)
    use b2mod_types
    implicit none
    integer, intent(in) :: nx,ny,nncut,nncutmax,ccut
    integer, intent(in) :: leftcut(nncutmax),rightcut(nncutmax), &
        & topcut(nncutmax),bottomcut(nncutmax), &
        & leftix(-1:nx,-1:ny),leftiy(-1:nx,-1:ny),&
        & rightix(-1:nx,-1:ny),rightiy(-1:nx,-1:ny), &
        & topix(-1:nx,-1:ny),topiy(-1:nx,-1:ny),&
        & bottomix(-1:nx,-1:ny),bottomiy(-1:nx,-1:ny), periodic_bc
    integer, intent(inout) :: cflag(-1:nx,-1:ny,CARREOUT_NCELLFLAGS)
    real (kind=R8), intent(in) :: &
        & crx(-1:nx,-1:ny,0:3), cry(-1:nx,-1:ny,0:3)

    integer, intent(out) :: &
        & region(-1:nx,-1:ny,0:2),nnreg(0:2), &
        & resignore(-1:nx,-1:ny,1:2)
    !
    !      Cylindrical slab (core) case
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
    !            Limiter case (ID=1)
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
    ! region 3 === inboard divertor
    ! region 4 === outboard divertor
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

    integer ix,iy,inseliy,inselix1,inselix2,iyt,geoType,iFace,offset
    integer ix1,ix2,ixn,iyn,iregn,ncount
    integer ixbreak
    integer lefttargetindex(2), righttargetindex(2)
    logical CellToTest
    intrinsic min, max
    !
    real (kind=R8) :: Xvertices(0:3), Yvertices(0:3)
    !
#ifndef BUILDING_CARRE
    call ipgetr ('b2agfs_geom_match_dist', geom_match_dist)
    call xertst (geom_match_dist.gt.0.0_R8, &
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
                Xvertices(:) = crx(ix,iy,:)
                Yvertices(:) = cry(ix,iy,:)
                geoType = cellGeoType(Xvertices, Yvertices)
                CellToTest = cflag(ix,iy,CELLFLAG_TYPE) == GRID_INTERNAL .or.  &
            &               (cflag(ix,iy,CELLFLAG_TYPE) == GRID_BOUNDARY .and. &
            &                cflag(ix,iy,CELLFLAG_TOPFACE) /= GRID_UNDEFINED .and. &
            &                geoType /= CGEO_TRIA_NOTOP)
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
                Xvertices(:) = crx(ix,iy,:)
                Yvertices(:) = cry(ix,iy,:)
                geoType = cellGeoType(Xvertices, Yvertices)
                CellToTest = cflag(ix,iy,CELLFLAG_TYPE) == GRID_INTERNAL .or.  &
            &               (cflag(ix,iy,CELLFLAG_TYPE) == GRID_BOUNDARY .and. &
            &                cflag(ix,iy,CELLFLAG_TOPFACE) /= GRID_UNDEFINED .and. &
            &                geoType /= CGEO_TRIA_NOTOP)
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
                 if (dist(ix1,iy,0,ix1,iy,2).lt.geom_match_dist .and. &
                   & dist(ix2,iy,1,ix2,iy,3).lt.geom_match_dist) then
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
    ! first pass: set region for the real cells (internal + boundary cells)
    if(nx.eq.1) then  ! 1d-radial
        do iy=-1,ny
            do ix=-1,nx
              if (isRealCell(cflag(ix,iy,CELLFLAG_TYPE))) region(ix,iy,0)=1
            enddo
        enddo
        nnreg(0)=1
    else if(ny.eq.1) then  ! 1d-parallel
        do iy=-1,ny
            do ix=-1,nx
              if (isRealCell(cflag(ix,iy,CELLFLAG_TYPE))) region(ix,iy,0)=1
            enddo
        enddo
        nnreg(0)=1
    else if (nncut.eq.0 .and. periodic_bc.eq.1) then   ! limiter and cylindrical slab cases
        do iy = -1, inseliy
            do ix = inselix1, inselix2
              if (isRealCell(cflag(ix,iy,CELLFLAG_TYPE))) region(ix,iy,0) = 1
            enddo
        enddo
        do iy = inseliy+1, ny
            do ix = -1, nx
              if (isRealCell(cflag(ix,iy,CELLFLAG_TYPE))) region(ix,iy,0) = 2
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
              if (isRealCell(cflag(ix,iy,CELLFLAG_TYPE))) region(ix,iy,0) = 1
            end do
        end do
        nnreg(0) = 1
    else
        if(topcut(1).ge.-1) then
            if(nncut.eq.1) then
                do ix=-1,leftcut(1)-1
                    do iy=-1,ny
                        if (.not.isRealCell(cflag(ix,iy,CELLFLAG_TYPE))) cycle
                        if (periodic_bc.ne.1.or.iy.lt.topcut(1)) &
                            &           region(ix,iy,0)=3
                        if(periodic_bc.eq.1.and. &
                            &           iy.ge.topcut(1).and.iy.lt.inseliy) &
                            &           region(ix,iy,0)=2
                        if(periodic_bc.eq.1.and.iy.ge.inseliy.and. &
                            &           ix.ge.inselix1.and.ix.le.inselix2) &
                            &           region(ix,iy,0)=5
                    enddo
                enddo
                do ix=leftcut(1),rightcut(1)-1
                    do iy=bottomcut(1),topcut(1)-1
                      if (isRealCell(cflag(ix,iy,CELLFLAG_TYPE))) region(ix,iy,0)=1
                    enddo
                    if(periodic_bc.eq.0) then
                        do iy=topcut(1),ny
                          if (isRealCell(cflag(ix,iy,CELLFLAG_TYPE))) region(ix,iy,0)=2
                        enddo
                    elseif(periodic_bc.eq.1) then
                        do iy=topcut(1),inseliy-1
                          if (isRealCell(cflag(ix,iy,CELLFLAG_TYPE))) region(ix,iy,0)=2
                        enddo
                        do iy=inseliy,ny
                          if (isRealCell(cflag(ix,iy,CELLFLAG_TYPE)) .and. &
                          &   ix.ge.inselix1.and.ix.le.inselix2) region(ix,iy,0)=5
                        enddo
                    endif
                enddo
                do ix=rightcut(1),nx
                    do iy=-1,ny
                        if (.not.isRealCell(cflag(ix,iy,CELLFLAG_TYPE))) cycle
                        if (periodic_bc.eq.0.or.iy.lt.topcut(1)) &
                            &           region(ix,iy,0)=4
                        if(periodic_bc.eq.1.and. &
                            &           iy.ge.topcut(1).and.iy.lt.inseliy) &
                            &           region(ix,iy,0)=2
                        if(periodic_bc.eq.1.and.iy.ge.inseliy.and. &
                            &           ix.ge.inselix1.and.ix.le.inselix2) &
                            &           region(ix,iy,0)=5
                    enddo
                enddo
                if(periodic_bc.eq.1) then
                    nnreg(0)=5
                else
                    nnreg(0)=4
                endif
            elseif(nncut.eq.2) then
                if (ccut.eq.0) then
#ifdef BUILDING_CARRE
                   stop 'Found ccut = 0 : central cut for DN not specified.'
#else
                   call xerrab ('Found ccut = 0 : central cut for DN not specified.')
#endif
                endif
                do ix = -1, leftcut(1)-1
                    do iy = -1, ny
                      if (isRealCell(cflag(ix,iy,CELLFLAG_TYPE))) region(ix,iy,0)=3
                    enddo
                enddo
                do ix = leftcut(1), leftcut(2)-1
                    do iy = max(bottomcut(1),bottomcut(2)), ny
                        if (.not.isRealCell(cflag(ix,iy,CELLFLAG_TYPE))) cycle
                        if (iy.lt.min(topcut(1),topcut(2))) then
                            region(ix,iy,0)=1
                        else
                            region(ix,iy,0)=2
                        endif
                    enddo
                enddo
                do iy = -1, ny
                    do ix = leftcut(2), ccut-1
                        if (isRealCell(cflag(ix,iy,CELLFLAG_TYPE))) region(ix,iy,0)=4
                    enddo
                    do ix = ccut, rightcut(2)-1
                        if (isRealCell(cflag(ix,iy,CELLFLAG_TYPE))) region(ix,iy,0)=7
                    enddo
                enddo
                do ix = rightcut(2), rightcut(1)-1
                    do iy = max(bottomcut(1),bottomcut(2)), ny
                        if (.not.isRealCell(cflag(ix,iy,CELLFLAG_TYPE))) cycle
                        if (iy.lt.min(topcut(1),topcut(2))) then
                            region(ix,iy,0)=5
                        else
                            region(ix,iy,0)=6
                        endif
                    enddo
                enddo
                do ix = rightcut(1), nx
                    do iy = -1, ny
                      if (isRealCell(cflag(ix,iy,CELLFLAG_TYPE))) region(ix,iy,0)=8
                    enddo
                enddo
                nnreg(0)=8
            endif
        endif
    endif

    ! Ghost cells inherit region number from internal neighbour cell
    do ix=-1,nx
        do iy=-1,ny
            if (.not.isBoundaryCell(cflag(ix, iy, CELLFLAG_TYPE))) cycle
            Xvertices(:) = crx(ix,iy,:)
            Yvertices(:) = cry(ix,iy,:)
            geoType = cellGeoType(Xvertices, Yvertices)

            if ( cflag(ix, iy, CELLFLAG_LEFTFACE) /= GRID_UNDEFINED &
                 & .and. geoType /= CGEO_TRIA_NOLEFT ) &
                 & region(leftix(ix,iy), leftiy(ix,iy), 0) = region(ix,iy,0)
            if ( cflag(ix, iy, CELLFLAG_BOTTOMFACE) /= GRID_UNDEFINED &
                 & .and. geoType /= CGEO_TRIA_NOBOT ) &
                 & region(bottomix(ix,iy), bottomiy(ix,iy), 0) = region(ix,iy,0)
            if ( cflag(ix, iy, CELLFLAG_RIGHTFACE) /= GRID_UNDEFINED &
                 & .and. geoType /= CGEO_TRIA_NORIGHT ) &
                 & region(rightix(ix,iy), rightiy(ix,iy), 0) = region(ix,iy,0)
            if ( cflag(ix, iy, CELLFLAG_TOPFACE) /= GRID_UNDEFINED &
                 & .and. geoType /= CGEO_TRIA_NOTOP ) &
                 & region(topix(ix,iy), topiy(ix,iy), 0) = region(ix,iy,0)
        end do
    end do

    ! Final pass for corner cells surrounded by ghosts
    do iy=-1,ny
        do ix=-1,nx
            if (.not.isGhostCell(cflag(ix,iy,CELLFLAG_TYPE))) cycle
            if (region(ix,iy,0).ne.0) cycle
            ! check for presence of ghost neighbors
            ncount=0
            iregn=0
            ixn=leftix(ix,iy)
            iyn=leftiy(ix,iy)
            if (isInDomain(nx,ny,ixn,iyn)) then
                if (isGhostcell(cflag(ixn,iyn,CELLFLAG_TYPE)).and.region(ixn,iyn,0).ne.0) then
                    iregn=region(ixn,iyn,0)
                endif
            endif
            ixn=rightix(ix,iy)
            iyn=rightiy(ix,iy)
            if (isInDomain(nx,ny,ixn,iyn)) then
                if (isGhostcell(cflag(ixn,iyn,CELLFLAG_TYPE)).and.region(ixn,iyn,0).ne.0) then
                    if (iregn.ne.0.and.iregn.ne.region(ixn,iyn,0)) then
#ifdef BUILDING_CARRE
                      stop 'Problem with region definition'
#else
                      call xerrab('Problem with region definition')
#endif
                    end if
                    iregn=region(ixn,iyn,0)
                endif
            endif
            ixn=topix(ix,iy)
            iyn=topiy(ix,iy)
            if (isInDomain(nx,ny,ixn,iyn)) then
                if (isGhostcell(cflag(ixn,iyn,CELLFLAG_TYPE)).and.region(ixn,iyn,0).ne.0) then
                    if (iregn.ne.0.and.iregn.ne.region(ixn,iyn,0)) then
#ifdef BUILDING_CARRE
                      stop 'Problem with region definition'
#else
                      call xerrab('Problem with region definition')
#endif
                    end if
                    iregn=region(ixn,iyn,0)
                endif
            endif
            ixn=bottomix(ix,iy)
            iyn=bottomiy(ix,iy)
            if (isInDomain(nx,ny,ixn,iyn)) then
                if (isGhostcell(cflag(ixn,iyn,CELLFLAG_TYPE)).and.region(ixn,iyn,0).ne.0) then
                    if (iregn.ne.0.and.iregn.ne.region(ixn,iyn,0)) then
#ifdef BUILDING_CARRE
                      stop 'Problem with region definition'
#else
                      call xerrab('Problem with region definition')
#endif
                    end if
                    iregn=region(ixn,iyn,0)
                endif
            endif
            if (iregn.ne.0) region(ix,iy,0)=iregn
            if (iregn.eq.0) then
               write (*,*) 'Problem identifying region for corner cell ', ix,iy
            endif
        end do
    end do

    ! Set unused cells to region 0
    do iy=-1,ny
        do ix=-1,nx
            if (isUnusedCell(cflag(ix,iy,CELLFLAG_TYPE))) &
                 & region(ix,iy,0) = 0
        end do
    end do

    ! Check all valid cells have a nonzero region number
    do ix=-1,nx
        do iy=-1,ny
            if(.not.isUnusedCell(cflag(ix,iy,CELLFLAG_TYPE)) .and. &
             & region(ix,iy,0).eq.0) write(*,*) 'REGION not set for ',ix,iy,0
        enddo
    enddo

    write(*,*) 'region(:,:,0)'
    do iy=ny,-1,-1
        write(*,'(3276(10i1,1x))') (region(ix,iy,0),ix=-1,nx)
    enddo

    ! After we know the volumetric region numbers, we fix the boundary indices for
    ! non-structure boundaries. CARRE2 sets them to -1 for all boundary faces not
    ! on a structure. Here they are set to -REGION.
    ! The corresponding face of the guard cell receives the same number.

    ! We first identify the inner X-point
    ! For DN cases, only the inner X-point is guaranteed to lie within the
    ! valid domain!
    do iy=-1,ny
        do ix=-1,nx
            if (cflag(ix, iy, CELLFLAG_TYPE) /= GRID_BOUNDARY) cycle
            do iFace = CELLFLAG_LEFTFACE, CELLFLAG_TOPFACE
                if (cflag(ix, iy, iFace) == BOUNDARY_NOSTRUCTURE) then
                    offset = (iFace - CELLFLAG_LEFTFACE + 1) * (-10)
                    cflag(ix, iy, iFace) = offset - region(ix, iy, 0)
                    Xvertices(:) = crx(ix,iy,:)
                    Yvertices(:) = cry(ix,iy,:)
                    geoType = cellGeoType(Xvertices, Yvertices)
                    if (iFace == CELLFLAG_LEFTFACE .and. &
                      & geoType /= CGEO_TRIA_NOLEFT) then
                      cflag(leftix(ix,iy),leftiy(ix,iy),CELLFLAG_RIGHTFACE) =  &
                       & cflag(ix, iy, iFace)
                    else if (iFace == CELLFLAG_RIGHTFACE .and. &
                      & geoType /= CGEO_TRIA_NORIGHT) then
                      cflag(rightix(ix,iy),rightiy(ix,iy),CELLFLAG_LEFTFACE) =  &
                       & cflag(ix, iy, iFace)
                    else if (iFace == CELLFLAG_TOPFACE .and. &
                      & geoType /= CGEO_TRIA_NOTOP) then
                      cflag(topix(ix,iy),topiy(ix,iy),CELLFLAG_BOTTOMFACE) =  &
                       & cflag(ix, iy, iFace)
                    else if (iFace == CELLFLAG_BOTTOMFACE .and. &
                      & geoType /= CGEO_TRIA_NOBOT) then
                      cflag(bottomix(ix,iy),bottomiy(ix,iy),CELLFLAG_TOPFACE) =  &
                       & cflag(ix, iy, iFace)
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
          Xvertices(:) = crx(ix,iy,:)
          Yvertices(:) = cry(ix,iy,:)
          geoType = cellGeoType(Xvertices, Yvertices)
          if (geoType /= CGEO_TRIA_NOLEFT .and. &
            & cflag(ix, iy, CELLFLAG_LEFTFACE) /= GRID_UNDEFINED) then
            lefttargetindex(1) = cflag(ix,iy,CELLFLAG_LEFTFACE)
          else if (geoType /= CGEO_TRIA_NOTOP .and. &
            & cflag(ix, iy, CELLFLAG_TOPFACE) /= GRID_UNDEFINED) then
            lefttargetindex(1) = cflag(ix,iy,CELLFLAG_TOPFACE)
          else if (geoType /= CGEO_TRIA_NOBOT .and. &
            & cflag(ix, iy, CELLFLAG_BOTTOMFACE) /= GRID_UNDEFINED) then
            lefttargetindex(1) = cflag(ix,iy,CELLFLAG_BOTTOMFACE)
          endif
        else
          ix = ix-1
        endif
      end do
#ifdef BUILDING_CARRE
      if (.not. lefttargetindex(1).ne.GRID_UNDEFINED) stop 'Left target not found!'
#else
      call xertst(lefttargetindex(1).ne.GRID_UNDEFINED, 'Left target not found!')
#endif
! We identify the structure index of the right target
      iy = inseliy+1
      ix = inselix2
      righttargetindex(1)=GRID_UNDEFINED
      do while (ix.le.nx .and. righttargetindex(1).eq.GRID_UNDEFINED)
        if (cflag(ix, iy, CELLFLAG_TYPE) == GRID_BOUNDARY) then
          Xvertices(:) = crx(ix,iy,:)
          Yvertices(:) = cry(ix,iy,:)
          geoType = cellGeoType(Xvertices, Yvertices)
          if (geoType /= CGEO_TRIA_NORIGHT .and. &
            & cflag(ix, iy, CELLFLAG_RIGHTFACE) /= GRID_UNDEFINED) then
            righttargetindex(1) = cflag(ix,iy,CELLFLAG_RIGHTFACE)
          else if (geoType /= CGEO_TRIA_NOTOP .and. &
            & cflag(ix, iy, CELLFLAG_TOPFACE) /= GRID_UNDEFINED) then
            righttargetindex(1) = cflag(ix,iy,CELLFLAG_TOPFACE)
          else if (geoType /= CGEO_TRIA_NOBOT .and. &
            & cflag(ix, iy, CELLFLAG_BOTTOMFACE) /= GRID_UNDEFINED) then
            righttargetindex(1) = cflag(ix,iy,CELLFLAG_BOTTOMFACE)
          endif
        else
          ix = ix+1
        endif
      end do
#ifdef BUILDING_CARRE
      if (.not. righttargetindex(1).ne.GRID_UNDEFINED) stop 'Right target not found!'
#else
      call xertst(righttargetindex(1).ne.GRID_UNDEFINED, 'Right target not found!')
#endif

    ! X-flux component
      do iy = inseliy+1, ny
! We find the cells contacting the left target: x-region 1
        ix = -1
        do while (ix.lt.nx .and. cflag(ix, iy, CELLFLAG_TYPE) /= GRID_INTERNAL)
          if (cflag(ix, iy, CELLFLAG_TYPE) == GRID_BOUNDARY) then
            Xvertices(:) = crx(ix,iy,:)
            Yvertices(:) = cry(ix,iy,:)
            geoType = cellGeoType(Xvertices, Yvertices)
            if (cflag(ix, iy, CELLFLAG_LEFTFACE) == lefttargetindex(1) &
                 & .and. geoType /= CGEO_TRIA_NOLEFT) &
                 & region(ix,iy,1)=1
            if (cflag(ix, iy, CELLFLAG_TOPFACE) == lefttargetindex(1) &
                 & .and. geoType /= CGEO_TRIA_NOTOP) &
                 & region(topix(ix,iy),topiy(ix,iy),2)=-1
            if (cflag(ix, iy, CELLFLAG_BOTTOMFACE) == lefttargetindex(1) &
                 & .and. geoType /= CGEO_TRIA_NOBOT) &
                 & region(ix,iy,2)=-1
          endif
          ix = ix + 1
        enddo
! We find the cells contacting the right target: x-region 2
        ix = nx
        do while (ix.gt.-1 .and. cflag(ix, iy, CELLFLAG_TYPE) /= GRID_INTERNAL)
          if (cflag(ix, iy, CELLFLAG_TYPE) == GRID_BOUNDARY) then
            Xvertices(:) = crx(ix,iy,:)
            Yvertices(:) = cry(ix,iy,:)
            geoType = cellGeoType(Xvertices, Yvertices)
            if (cflag(ix, iy, CELLFLAG_RIGHTFACE) == righttargetindex(1) &
                 & .and. geoType /= CGEO_TRIA_NORIGHT) &
                 & region(rightix(ix,iy),rightiy(ix,iy),1)=2
            if (cflag(ix, iy, CELLFLAG_TOPFACE) == righttargetindex(1) &
                 & .and. geoType /= CGEO_TRIA_NOTOP) &
                 & region(topix(ix,iy),topiy(ix,iy),2)=-2
            if (cflag(ix, iy, CELLFLAG_BOTTOMFACE) == righttargetindex(1) &
                 & .and. geoType /= CGEO_TRIA_NOBOT) &
                 & region(ix,iy,2)=-2
          endif
          ix = ix - 1
        enddo
      enddo
! We find the core periodicity line: x-region 3
      do iy = -1, inseliy
        if (.not.isUnusedCell(cflag(inselix1,iy,CELLFLAG_TYPE))) region(inselix1,iy,1)=3
      enddo
      nnreg(1)=3
    ! Y-flux component
      do ix = inselix1, inselix2
! We find the core boundary: y-region 1
        iy = -1
        do while (iy.le.inseliy .and. .not.isRealCell(cflag(ix, iy, CELLFLAG_TYPE)))
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
        do while (iy.gt.inseliy .and. cflag(ix,iy,CELLFLAG_TYPE) /= GRID_INTERNAL)
          if (cflag(ix, iy, CELLFLAG_TYPE) == GRID_BOUNDARY) then
            Xvertices(:) = crx(ix,iy,:)
            Yvertices(:) = cry(ix,iy,:)
            geoType = cellGeoType(Xvertices, Yvertices)
            if (cflag(ix, iy, CELLFLAG_TOPFACE) /= lefttargetindex(1) .and. &
              & cflag(ix, iy, CELLFLAG_TOPFACE) /= righttargetindex(1) .and. &
              & cflag(ix, iy, CELLFLAG_TOPFACE) /= GRID_UNDEFINED .and. &
              & geoType /= CGEO_TRIA_NOTOP) &
                 & region(topix(ix,iy),topiy(ix,iy),2)=3
            if (cflag(ix, iy, CELLFLAG_RIGHTFACE) /= lefttargetindex(1) .and. &
              & cflag(ix, iy, CELLFLAG_RIGHTFACE) /= righttargetindex(1) .and. &
              & cflag(ix, iy, CELLFLAG_RIGHTFACE) /= GRID_UNDEFINED .and. &
              & geoType /= CGEO_TRIA_NORIGHT) &
                 & region(rightix(ix,iy),rightiy(ix,iy),1)=-3
            if (cflag(ix, iy, CELLFLAG_LEFTFACE) /= lefttargetindex(1) .and. &
              & cflag(ix, iy, CELLFLAG_LEFTFACE) /= righttargetindex(1) .and. &
              & cflag(ix, iy, CELLFLAG_LEFTFACE) /= GRID_UNDEFINED .and. &
              & geoType /= CGEO_TRIA_NOLEFT) &
                 & region(ix,iy,1)=-3
          end if
          iy = iy - 1
        end do
      end do
      nnreg(2)=3
    else if (nncut.eq.0 .and. periodic_bc.eq.1 .and. nnreg(0).eq.1) then ! cylindrical slab case
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
      if (.not.isUnusedCell(cflag(ix,iy,CELLFLAG_TYPE))) then
        do while (ix.gt.-1 .and. lefttargetindex(1).eq.GRID_UNDEFINED)
          if (cflag(ix, iy, CELLFLAG_TYPE) == GRID_BOUNDARY) then
            Xvertices(:) = crx(ix,iy,:)
            Yvertices(:) = cry(ix,iy,:)
            geoType = cellGeoType(Xvertices, Yvertices)
            if (geoType /= CGEO_TRIA_NOLEFT .and. &
              & cflag(ix, iy, CELLFLAG_LEFTFACE) /= GRID_UNDEFINED) then
              lefttargetindex(1) = cflag(ix,iy,CELLFLAG_LEFTFACE)
            else if (geoType /= CGEO_TRIA_NOTOP .and. &
              & cflag(ix, iy, CELLFLAG_TOPFACE) /= GRID_UNDEFINED) then
              lefttargetindex(1) = cflag(ix,iy,CELLFLAG_TOPFACE)
            else if (geoType /= CGEO_TRIA_NOBOT .and. &
              & cflag(ix, iy, CELLFLAG_BOTTOMFACE) /= GRID_UNDEFINED) then
              lefttargetindex(1) = cflag(ix,iy,CELLFLAG_BOTTOMFACE)
            endif
          else
            ix = ix-1
          endif
        end do
#ifdef BUILDING_CARRE
        if (.not. lefttargetindex(1).ne.GRID_UNDEFINED) stop 'Left target not found!'
#else
        call xertst(lefttargetindex(1).ne.GRID_UNDEFINED, 'Left target not found!')
#endif
      endif
! We identify the structure index of the right target
      iy = topcut(1)+1
      ix = rightcut(1)
      righttargetindex(1)=GRID_UNDEFINED
      if (.not.isUnusedCell(cflag(ix,iy,CELLFLAG_TYPE))) then
        do while (ix.lt.nx .and. righttargetindex(1).eq.GRID_UNDEFINED)
          if (cflag(ix, iy, CELLFLAG_TYPE) == GRID_BOUNDARY) then
            Xvertices(:) = crx(ix,iy,:)
            Yvertices(:) = cry(ix,iy,:)
            geoType = cellGeoType(Xvertices, Yvertices)
            if (geoType /= CGEO_TRIA_NORIGHT .and. &
              & cflag(ix, iy, CELLFLAG_RIGHTFACE) /= GRID_UNDEFINED) then
              righttargetindex(1) = cflag(ix,iy,CELLFLAG_RIGHTFACE)
            else if (geoType /= CGEO_TRIA_NOTOP .and. &
              & cflag(ix, iy, CELLFLAG_TOPFACE) /= GRID_UNDEFINED) then
              righttargetindex(1) = cflag(ix,iy,CELLFLAG_TOPFACE)
            else if (geoType /= CGEO_TRIA_NOBOT .and. &
              & cflag(ix, iy, CELLFLAG_BOTTOMFACE) /= GRID_UNDEFINED) then
              righttargetindex(1) = cflag(ix,iy,CELLFLAG_BOTTOMFACE)
            endif
          else
            ix = ix+1
          endif
        end do
#ifdef BUILDING_CARRE
        if (.not. righttargetindex(1).ne.GRID_UNDEFINED) stop 'Right target not found!'
#else
        call xertst(righttargetindex(1).ne.GRID_UNDEFINED, 'Right target not found!')
#endif
      endif

    ! X-flux component
      do iy=-1,iyt
! We find the cells contacting the left target: x-region 1
        ix = -1
        do while (ix.lt.nx .and. cflag(ix, iy, CELLFLAG_TYPE) /= GRID_INTERNAL)
          if (cflag(ix, iy, CELLFLAG_TYPE) == GRID_BOUNDARY) then
            Xvertices(:) = crx(ix,iy,:)
            Yvertices(:) = cry(ix,iy,:)
            geoType = cellGeoType(Xvertices, Yvertices)
            if (cflag(ix, iy, CELLFLAG_LEFTFACE) == lefttargetindex(1) &
                 & .and. lefttargetindex(1) /= GRID_UNDEFINED &
                 & .and. geoType /= CGEO_TRIA_NOLEFT) &
                 & region(ix,iy,1)=1
            if (cflag(ix, iy, CELLFLAG_TOPFACE) == lefttargetindex(1) &
                 & .and. lefttargetindex(1) /= GRID_UNDEFINED &
                 & .and. geoType /= CGEO_TRIA_NOTOP) &
                 & region(topix(ix,iy),topiy(ix,iy),2)=-1
            if (cflag(ix, iy, CELLFLAG_BOTTOMFACE) == lefttargetindex(1) &
                 & .and. lefttargetindex(1) /= GRID_UNDEFINED &
                 & .and. geoType /= CGEO_TRIA_NOBOT) &
                 & region(ix,iy,2)=-1
          endif
          ix = ix + 1
        enddo
! We find the cells contacting the right target: x-region 4 (1 X-point) or 8 (2 X-points)
        ix = nx
        do while (ix.gt.-1 .and. cflag(ix, iy, CELLFLAG_TYPE) /= GRID_INTERNAL)
          if (cflag(ix, iy, CELLFLAG_TYPE) == GRID_BOUNDARY) then
            Xvertices(:) = crx(ix,iy,:)
            Yvertices(:) = cry(ix,iy,:)
            geoType = cellGeoType(Xvertices, Yvertices)
            if (cflag(ix, iy, CELLFLAG_RIGHTFACE) == righttargetindex(1) &
                 & .and. righttargetindex(1) /= GRID_UNDEFINED &
                 & .and. geoType /= CGEO_TRIA_NORIGHT) &
                 & region(rightix(ix,iy),rightiy(ix,iy),1)=4*nncut
            if (cflag(ix, iy, CELLFLAG_TOPFACE) == righttargetindex(1) &
                 & .and. righttargetindex(1) /= GRID_UNDEFINED &
                 & .and. geoType /= CGEO_TRIA_NOTOP) &
                 & region(topix(ix,iy),topiy(ix,iy),2)=-4*nncut
            if (cflag(ix, iy, CELLFLAG_BOTTOMFACE) == righttargetindex(1) &
                 & .and. righttargetindex(1) /= GRID_UNDEFINED &
                 & .and. geoType /= CGEO_TRIA_NOBOT) &
                 & region(ix,iy,2)=-4*nncut
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
        if (ixbreak.eq.-2) &
         & stop 'Did not find the break location!'
#else
        call xertst(ixbreak.ne.-2, &
         & 'Did not find the break location!')
#endif
! We identify the structure index of the second left target
        iy = topcut(2)+1
        ix = rightcut(2)-1
        do while (.not.isRealCell(cflag(ix, iy, CELLFLAG_TYPE)) .and. &
         & ix.le.nx)
          ix = ix + 1
        end do
        lefttargetindex(2)=GRID_UNDEFINED
        if (ix.le.nx) then
          do while (ix.gt.ixbreak .and. lefttargetindex(2).eq.GRID_UNDEFINED)
            if (cflag(ix, iy, CELLFLAG_TYPE) == GRID_BOUNDARY) then
              Xvertices(:) = crx(ix,iy,:)
              Yvertices(:) = cry(ix,iy,:)
              geoType = cellGeoType(Xvertices, Yvertices)
              if (geoType /= CGEO_TRIA_NOLEFT .and. &
                & cflag(ix, iy, CELLFLAG_LEFTFACE) /= GRID_UNDEFINED) then
                lefttargetindex(2) = cflag(ix,iy,CELLFLAG_LEFTFACE)
              else if (geoType /= CGEO_TRIA_NOTOP .and. &
                & cflag(ix, iy, CELLFLAG_TOPFACE) /= GRID_UNDEFINED) then
                lefttargetindex(2) = cflag(ix,iy,CELLFLAG_TOPFACE)
              else if (geoType /= CGEO_TRIA_NOBOT .and. &
                & cflag(ix, iy, CELLFLAG_BOTTOMFACE) /= GRID_UNDEFINED) then
                lefttargetindex(2) = cflag(ix,iy,CELLFLAG_BOTTOMFACE)
              endif
            else
              ix = ix-1
            endif
          end do
#ifdef BUILDING_CARRE
          if (.not. lefttargetindex(2).ne.GRID_UNDEFINED) &
           & stop 'Second left target not found!'
#else
          call xertst(lefttargetindex(2).ne.GRID_UNDEFINED, &
           & 'Second left target not found!')
#endif
        endif
! We identify the structure index of the second right target
        iy = topcut(2)+1
        ix = leftcut(2)
        do while (.not.isRealCell(cflag(ix, iy, CELLFLAG_TYPE)) .and. &
         & ix.ge.0 )
          ix = ix - 1
        end do
        righttargetindex(2)=GRID_UNDEFINED
        if (ix.ge.0) then
          do while (ix.lt.ixbreak .and. righttargetindex(2).eq.GRID_UNDEFINED)
            if (cflag(ix, iy, CELLFLAG_TYPE) == GRID_BOUNDARY) then
              Xvertices(:) = crx(ix,iy,:)
              Yvertices(:) = cry(ix,iy,:)
              geoType = cellGeoType(Xvertices, Yvertices)
              if (geoType /= CGEO_TRIA_NORIGHT .and. &
                & cflag(ix, iy, CELLFLAG_RIGHTFACE) /= GRID_UNDEFINED) then
                righttargetindex(2) = cflag(ix,iy,CELLFLAG_RIGHTFACE)
              else if (geoType /= CGEO_TRIA_NOTOP .and. &
                & cflag(ix, iy, CELLFLAG_TOPFACE) /= GRID_UNDEFINED) then
                righttargetindex(2) = cflag(ix,iy,CELLFLAG_TOPFACE)
              else if (geoType /= CGEO_TRIA_NOBOT .and. &
                & cflag(ix, iy, CELLFLAG_BOTTOMFACE) /= GRID_UNDEFINED) then
                righttargetindex(2) = cflag(ix,iy,CELLFLAG_BOTTOMFACE)
              endif
            else
              ix = ix+1
            endif
          end do
#ifdef BUILDING_CARRE
          if (.not. righttargetindex(2).ne.GRID_UNDEFINED) &
      &    stop 'Second right target not found!'
#else
          call xertst(righttargetindex(2).ne.GRID_UNDEFINED, &
      &    'Second right target not found!')
#endif
        endif
        do iy=-1,iyt
! We find the cells contacting the second left target: x-region 5
          ix = nx
          do while (ix.gt.ixbreak)
            if (cflag(ix, iy, CELLFLAG_TYPE) == GRID_BOUNDARY) then
              Xvertices(:) = crx(ix,iy,:)
              Yvertices(:) = cry(ix,iy,:)
              geoType = cellGeoType(Xvertices, Yvertices)
              if (cflag(ix, iy, CELLFLAG_LEFTFACE) == lefttargetindex(2) &
                 & .and. lefttargetindex(2) /= GRID_UNDEFINED &
                 & .and. geoType /= CGEO_TRIA_NOLEFT) region(ix,iy,1)=5
              if (cflag(ix, iy, CELLFLAG_TOPFACE) == lefttargetindex(2) &
                 & .and. lefttargetindex(2) /= GRID_UNDEFINED &
                 & .and. geoType /= CGEO_TRIA_NOTOP) &
                 & region(topix(ix,iy),topiy(ix,iy),2)=-5
              if (cflag(ix, iy, CELLFLAG_BOTTOMFACE) == lefttargetindex(2) &
                 & .and. lefttargetindex(2) /= GRID_UNDEFINED &
                 & .and. geoType /= CGEO_TRIA_NOBOT) region(ix,iy,2)=-5
            endif
            ix = ix - 1
          enddo
! We find the cells contacting the second right target: x-region 4
          ix = -1
          do while (ix.lt.ixbreak)
            if (cflag(ix, iy, CELLFLAG_TYPE) == GRID_BOUNDARY) then
              Xvertices(:) = crx(ix,iy,:)
              Yvertices(:) = cry(ix,iy,:)
              geoType = cellGeoType(Xvertices, Yvertices)
              if (cflag(ix, iy, CELLFLAG_RIGHTFACE) == righttargetindex(2) &
                 & .and. righttargetindex(2) /= GRID_UNDEFINED &
                 & .and. geoType /= CGEO_TRIA_NORIGHT) &
                 & region(rightix(ix,iy),rightiy(ix,iy),1)=4
              if (cflag(ix, iy, CELLFLAG_TOPFACE) == righttargetindex(2) &
                 & .and. righttargetindex(2) /= GRID_UNDEFINED &
                 & .and. geoType /= CGEO_TRIA_NOTOP) &
                 & region(topix(ix,iy),topiy(ix,iy),2)=-4
              if (cflag(ix, iy, CELLFLAG_BOTTOMFACE) == righttargetindex(2) &
                 & .and. righttargetindex(2) /= GRID_UNDEFINED &
                 & .and. geoType /= CGEO_TRIA_NOBOT) &
                 & region(ix,iy,2)=-4
            endif
            ix = ix + 1
          enddo
        enddo
      end if  ! DN
      do iy=topcut(1),iyt
! We find the location of the inner divertor throat: x-region 2
        Xvertices(:) = crx(leftcut(1),iy,:)
        Yvertices(:) = cry(leftcut(1),iy,:)
        geoType = cellGeoType(Xvertices, Yvertices)
        if (isRealCell(cflag(leftcut(1),iy,CELLFLAG_TYPE)) &
                 & .and. geoType /= CGEO_TRIA_NOLEFT) &
     &   region(leftcut(1),iy,1) = 2
! We find the location of the outer divertor throat: x-region 3 (1 X-point) or 7 (2 X-points)
        Xvertices(:) = crx(rightcut(1),iy,:)
        Yvertices(:) = cry(rightcut(1),iy,:)
        geoType = cellGeoType(Xvertices, Yvertices)
        if (isRealCell(cflag(rightcut(1),iy,CELLFLAG_TYPE)) &
                 & .and. geoType /= CGEO_TRIA_NOLEFT) &
     &   region(rightcut(1),iy,1) = 3 + 4*(nncut-1)
      end do
      if (nncut.eq.2) then
        do iy=topcut(2),iyt
! We find the location of the inner upper divertor throat: x-region 3
          Xvertices(:) = crx(leftcut(2),iy,:)
          Yvertices(:) = cry(leftcut(2),iy,:)
          geoType = cellGeoType(Xvertices, Yvertices)
          if (isRealCell(cflag(leftcut(2),iy,CELLFLAG_TYPE)) &
                 & .and. geoType /= CGEO_TRIA_NOLEFT) &
     &     region(leftcut(2),iy,1) = 3
! We find the location of the outer upper divertor throat: x-region 6
          Xvertices(:) = crx(rightcut(2),iy,:)
          Yvertices(:) = cry(rightcut(2),iy,:)
          geoType = cellGeoType(Xvertices, Yvertices)
          if (isRealCell(cflag(rightcut(2),iy,CELLFLAG_TYPE)) &
                 & .and. geoType /= CGEO_TRIA_NOLEFT) &
     &     region(rightcut(2),iy,1) = 6
        end do
      end if
! We find the core and PFR bottom periodicity lines: x-regions 5 and 6 (1 X-point)
!                                                  : x-regions 9 and 12 (2 X-points)
! x-region 13 is the connection between the two halves of the SOL section between the two separatrices
      do iy=bottomcut(1),topcut(1)-1
        if(nncut.eq.1) then
          if (isRealCell(cflag(leftcut(1),iy,CELLFLAG_TYPE))) &
      &    region(leftcut(1),iy,1)=5
          Xvertices(:) = crx(rightcut(1),iy,:)
          Yvertices(:) = cry(rightcut(1),iy,:)
          geoType = cellGeoType(Xvertices, Yvertices)
          if (isRealCell(cflag(rightcut(1),iy,CELLFLAG_TYPE)) .and. &
            & cflag(rightcut(1), iy, CELLFLAG_LEFTFACE) == GRID_UNDEFINED .and. &
            & geoType /= CGEO_TRIA_NOLEFT) region(rightcut(1),iy,1)=6
        elseif(nncut.eq.2) then
          if(iy.ge.bottomcut(2).and.iy.lt.topcut(2)) then
            if (isRealCell(cflag(leftcut(1),iy,CELLFLAG_TYPE))) region(leftcut(1),iy,1)=9
          else
            if (isRealCell(cflag(leftcut(1),iy,CELLFLAG_TYPE))) region(leftcut(1),iy,1)=13
          endif
          Xvertices(:) = crx(rightcut(1),iy,:)
          Yvertices(:) = cry(rightcut(1),iy,:)
          geoType = cellGeoType(Xvertices, Yvertices)
          if (isRealCell(cflag(rightcut(1),iy,CELLFLAG_TYPE)) .and. &
            & cflag(rightcut(1), iy, CELLFLAG_LEFTFACE) == GRID_UNDEFINED) &
            & region(rightcut(1),iy,1)=12
        endif
      enddo
! We find the core and PFR top periodicity lines: x-regions 10 and 11 (2 X-points)
! x-region 13 is the connection between the two halves of the SOL section between the two separatrices
      if(nncut.eq.2) then
        do iy=bottomcut(2),topcut(2)-1
          Xvertices(:) = crx(leftcut(2),iy,:)
          Yvertices(:) = cry(leftcut(2),iy,:)
          geoType = cellGeoType(Xvertices, Yvertices)
          if (isRealCell(cflag(leftcut(2),iy,CELLFLAG_TYPE)) .and. &
            & cflag(leftcut(2), iy, CELLFLAG_LEFTFACE) == GRID_UNDEFINED .and. &
            & geoType /= CGEO_TRIA_NOLEFT) region(leftcut(2),iy,1)=10
          if(iy.ge.bottomcut(1).and.iy.lt.topcut(1)) then
            if (isRealCell(cflag(rightcut(2),iy,CELLFLAG_TYPE))) region(rightcut(2),iy,1)=11
          else
            if (isRealCell(cflag(rightcut(2),iy,CELLFLAG_TYPE))) region(rightcut(2),iy,1)=13
          endif
        enddo
      endif
! We find the stellarator island periodicity line: x-region 7
      if(nncut.eq.1.and.periodic_bc.eq.1) then
        do iy=iyt+1,ny
          if (isRealCell(cflag(inselix1,iy,CELLFLAG_TYPE))) region(inselix1,iy,1)=7
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
        do while (iy.lt.topcut(1) .and. cflag(ix,iy,CELLFLAG_TYPE) /= GRID_INTERNAL)
          if (cflag(ix, iy, CELLFLAG_TYPE) == GRID_BOUNDARY) then
            Xvertices(:) = crx(ix,iy,:)
            Yvertices(:) = cry(ix,iy,:)
            geoType = cellGeoType(Xvertices, Yvertices)
            if (cflag(ix, iy, CELLFLAG_BOTTOMFACE) /= lefttargetindex(1) .and. &
              & cflag(ix, iy, CELLFLAG_BOTTOMFACE) /= righttargetindex(1) .and. &
              & cflag(ix, iy, CELLFLAG_BOTTOMFACE) /= GRID_UNDEFINED .and. &
              & geoType /= CGEO_TRIA_NOBOT) &
                 & region(ix,iy,2)=1
            if (cflag(ix, iy, CELLFLAG_RIGHTFACE) /= lefttargetindex(1) .and. &
              & cflag(ix, iy, CELLFLAG_RIGHTFACE) /= righttargetindex(1) .and. &
              & cflag(ix, iy, CELLFLAG_RIGHTFACE) /= GRID_UNDEFINED .and. &
              & geoType /= CGEO_TRIA_NORIGHT) &
                 & region(rightix(ix,iy),rightiy(ix,iy),1)=-1
            if (cflag(ix, iy, CELLFLAG_LEFTFACE) /= lefttargetindex(1) .and. &
              & cflag(ix, iy, CELLFLAG_LEFTFACE) /= righttargetindex(1) .and. &
              & cflag(ix, iy, CELLFLAG_LEFTFACE) /= GRID_UNDEFINED .and. &
              & geoType /= CGEO_TRIA_NOLEFT) &
                 & region(ix,iy,1)=-1
          end if
          iy = iy + 1
        end do
      end do
      if (periodic_bc.ne.1) then
! We find the contact to the left (bottom) divertor wall: y-region 5
        do ix = -1, leftcut(1)-1
          iy = ny
          do while (iy.ge.topcut(1) .and. cflag(ix,iy,CELLFLAG_TYPE) /= GRID_INTERNAL)
            if (cflag(ix, iy, CELLFLAG_TYPE) == GRID_BOUNDARY) then
              Xvertices(:) = crx(ix,iy,:)
              Yvertices(:) = cry(ix,iy,:)
              geoType = cellGeoType(Xvertices, Yvertices)
              if (cflag(ix, iy, CELLFLAG_TOPFACE) /= lefttargetindex(1) .and. &
                & cflag(ix, iy, CELLFLAG_TOPFACE) /= righttargetindex(1) .and. &
                & cflag(ix, iy, CELLFLAG_TOPFACE) /= GRID_UNDEFINED .and. &
                & geoType /= CGEO_TRIA_NOTOP) &
                & region(topix(ix,iy),topiy(ix,iy),2)=5
              if (cflag(ix, iy, CELLFLAG_RIGHTFACE) /= lefttargetindex(1) .and. &
                & cflag(ix, iy, CELLFLAG_RIGHTFACE) /= righttargetindex(1) .and. &
                & cflag(ix, iy, CELLFLAG_RIGHTFACE) /= GRID_UNDEFINED .and. &
                & geoType /= CGEO_TRIA_NORIGHT) &
                & region(rightix(ix,iy),rightiy(ix,iy),1)=-5
              if (cflag(ix, iy, CELLFLAG_LEFTFACE) /= lefttargetindex(1) .and. &
                & cflag(ix, iy, CELLFLAG_LEFTFACE) /= righttargetindex(1) .and. &
                & cflag(ix, iy, CELLFLAG_LEFTFACE) /= GRID_UNDEFINED .and. &
                & geoType /= CGEO_TRIA_NOLEFT) &
                & region(ix,iy,1)=-5
            end if
            iy = iy - 1
          end do
        end do
      elseif (periodic_bc.eq.1) then
! Stellarator island: We find the entrance to the left PFR: y-region 5
        do ix=-1,leftcut(1)-1
          if (isRealCell(cflag(ix,iyt,CELLFLAG_TYPE))) region(ix,iyt,2)=5
        enddo
      endif
      if(nncut.lt.2) then  ! SN case
        do ix=leftcut(1),rightcut(1)-1
! SN: We find the core boundary: y-region 2
          iy = -1
          do while (iy.le.topcut(1) .and. .not.isRealCell(cflag(ix, iy, CELLFLAG_TYPE)))
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
            do while (iy.ge.inseliy .and. .not.isRealCell(cflag(ix, iy, CELLFLAG_TYPE)))
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
            do while (iy.ge.topcut(1) .and. cflag(ix,iy,CELLFLAG_TYPE) /= GRID_INTERNAL)
              if (cflag(ix, iy, CELLFLAG_TYPE) == GRID_BOUNDARY) then
                Xvertices(:) = crx(ix,iy,:)
                Yvertices(:) = cry(ix,iy,:)
                geoType = cellGeoType(Xvertices, Yvertices)
                if (cflag(ix, iy, CELLFLAG_TOPFACE) /= lefttargetindex(1) .and. &
                  & cflag(ix, iy, CELLFLAG_TOPFACE) /= righttargetindex(1) .and. &
                  & cflag(ix, iy, CELLFLAG_TOPFACE) /= GRID_UNDEFINED .and. &
                  & geoType /= CGEO_TRIA_NOTOP) &
                  & region(topix(ix,iy),topiy(ix,iy),2)=6
                if (cflag(ix, iy, CELLFLAG_RIGHTFACE) /= lefttargetindex(1) .and. &
                  & cflag(ix, iy, CELLFLAG_RIGHTFACE) /= righttargetindex(1) .and. &
                  & cflag(ix, iy, CELLFLAG_RIGHTFACE) /= GRID_UNDEFINED .and. &
                  & geoType /= CGEO_TRIA_NORIGHT) &
                  & region(rightix(ix,iy),rightiy(ix,iy),1)=-6
                if (cflag(ix, iy, CELLFLAG_LEFTFACE) /= lefttargetindex(1) .and. &
                  & cflag(ix, iy, CELLFLAG_LEFTFACE) /= righttargetindex(1) .and. &
                  & cflag(ix, iy, CELLFLAG_LEFTFACE) /= GRID_UNDEFINED .and. &
                  & geoType /= CGEO_TRIA_NOLEFT) &
                  & region(ix,iy,1)=-6
              end if
              iy = iy - 1
            end do
          end do
        endif
! SN: We find the contact to the right PFR wall: y-region 3
        do ix = rightcut(1), nx
          iy = -1
          do while (iy.le.topcut(1) .and. cflag(ix,iy,CELLFLAG_TYPE) /= GRID_INTERNAL)
            if (cflag(ix, iy, CELLFLAG_TYPE) == GRID_BOUNDARY) then
              Xvertices(:) = crx(ix,iy,:)
              Yvertices(:) = cry(ix,iy,:)
              geoType = cellGeoType(Xvertices, Yvertices)
              if (cflag(ix, iy, CELLFLAG_BOTTOMFACE) /= lefttargetindex(1) .and. &
                & cflag(ix, iy, CELLFLAG_BOTTOMFACE) /= righttargetindex(1) .and. &
                & cflag(ix, iy, CELLFLAG_BOTTOMFACE) /= GRID_UNDEFINED .and. &
                & geoType /= CGEO_TRIA_NOBOT) &
                & region(ix,iy,2)=3
              if (cflag(ix, iy, CELLFLAG_RIGHTFACE) /= lefttargetindex(1) .and. &
                & cflag(ix, iy, CELLFLAG_RIGHTFACE) /= righttargetindex(1) .and. &
                & cflag(ix, iy, CELLFLAG_RIGHTFACE) /= GRID_UNDEFINED .and. &
                & geoType /= CGEO_TRIA_NORIGHT) &
                & region(rightix(ix,iy),rightiy(ix,iy),1)=-3
              if (cflag(ix, iy, CELLFLAG_LEFTFACE) /= lefttargetindex(1) .and. &
                & cflag(ix, iy, CELLFLAG_LEFTFACE) /= righttargetindex(1) .and. &
                & cflag(ix, iy, CELLFLAG_LEFTFACE) /= GRID_UNDEFINED .and. &
                & geoType /= CGEO_TRIA_NOLEFT) &
                & region(ix,iy,1)=-3
            end if
            iy = iy + 1
          end do
! SN: We find the contact to the right divertor wall: y-region 7
          iy = ny
          do while (iy.gt.topcut(1) .and. cflag(ix,iy,CELLFLAG_TYPE) /= GRID_INTERNAL)
            if (cflag(ix, iy, CELLFLAG_TYPE) == GRID_BOUNDARY) then
              Xvertices(:) = crx(ix,iy,:)
              Yvertices(:) = cry(ix,iy,:)
              geoType = cellGeoType(Xvertices, Yvertices)
              if (cflag(ix, iy, CELLFLAG_TOPFACE) /= lefttargetindex(1) .and. &
                & cflag(ix, iy, CELLFLAG_TOPFACE) /= righttargetindex(1) .and. &
                & cflag(ix, iy, CELLFLAG_TOPFACE) /= GRID_UNDEFINED .and. &
                & geoType /= CGEO_TRIA_NOTOP) &
                & region(topix(ix,iy),topiy(ix,iy),2)=7
              if (cflag(ix, iy, CELLFLAG_RIGHTFACE) /= lefttargetindex(1) .and. &
                & cflag(ix, iy, CELLFLAG_RIGHTFACE) /= righttargetindex(1) .and. &
                & cflag(ix, iy, CELLFLAG_RIGHTFACE) /= GRID_UNDEFINED .and. &
                & geoType /= CGEO_TRIA_NORIGHT) &
                & region(rightix(ix,iy),rightiy(ix,iy),1)=-7
              if (cflag(ix, iy, CELLFLAG_LEFTFACE) /= lefttargetindex(1) .and. &
                & cflag(ix, iy, CELLFLAG_LEFTFACE) /= righttargetindex(1) .and. &
                & cflag(ix, iy, CELLFLAG_LEFTFACE) /= GRID_UNDEFINED .and. &
                & geoType /= CGEO_TRIA_NOLEFT) &
                & region(ix,iy,1)=-7
            end if
            iy = iy - 1
          end do
        end do
      elseif(nncut.eq.2) then  !DN case
        do ix=leftcut(1),leftcut(2)-1
! DN: We find the inner core boundary: y-region 2
          iy = -1
          do while (iy.le.topcut(1) .and. iy.le.topcut(2) .and. &
                  & .not.isRealCell(cflag(ix, iy, CELLFLAG_TYPE)))
            iy = iy+1
          end do
          if (iy.le.topcut(1) .and. iy.le.topcut(2)) then
            region(ix,iy,2)=2
          endif
! DN: We find the inner separatrix: y-region 4
          region(ix,min(topcut(1),topcut(2)),2)=4
! DN: We find the contact to the inner main chamber wall: y-region 6
          iy = ny
          do while (iy.ge.topcut(1) .and. iy.ge.topcut(2) .and. &
                  & cflag(ix,iy,CELLFLAG_TYPE) /= GRID_INTERNAL)
            if (cflag(ix, iy, CELLFLAG_TYPE) == GRID_BOUNDARY) then
              Xvertices(:) = crx(ix,iy,:)
              Yvertices(:) = cry(ix,iy,:)
              geoType = cellGeoType(Xvertices, Yvertices)
              if (cflag(ix, iy, CELLFLAG_TOPFACE) /= lefttargetindex(1) .and. &
                & cflag(ix, iy, CELLFLAG_TOPFACE) /= righttargetindex(1) .and. &
                & cflag(ix, iy, CELLFLAG_TOPFACE) /= lefttargetindex(2) .and. &
                & cflag(ix, iy, CELLFLAG_TOPFACE) /= righttargetindex(2) .and. &
                & cflag(ix, iy, CELLFLAG_TOPFACE) /= GRID_UNDEFINED .and. &
                & geoType /= CGEO_TRIA_NOTOP) &
                & region(topix(ix,iy),topiy(ix,iy),2)=6
              if (cflag(ix, iy, CELLFLAG_RIGHTFACE) /= lefttargetindex(1) .and. &
                & cflag(ix, iy, CELLFLAG_RIGHTFACE) /= righttargetindex(1) .and. &
                & cflag(ix, iy, CELLFLAG_RIGHTFACE) /= lefttargetindex(2) .and. &
                & cflag(ix, iy, CELLFLAG_RIGHTFACE) /= righttargetindex(2) .and. &
                & cflag(ix, iy, CELLFLAG_RIGHTFACE) /= GRID_UNDEFINED .and. &
                & geoType /= CGEO_TRIA_NORIGHT) &
                & region(rightix(ix,iy),rightiy(ix,iy),1)=-6
              if (cflag(ix, iy, CELLFLAG_LEFTFACE) /= lefttargetindex(1) .and. &
                & cflag(ix, iy, CELLFLAG_LEFTFACE) /= righttargetindex(1) .and. &
                & cflag(ix, iy, CELLFLAG_LEFTFACE) /= lefttargetindex(2) .and. &
                & cflag(ix, iy, CELLFLAG_LEFTFACE) /= righttargetindex(2) .and. &
                & cflag(ix, iy, CELLFLAG_LEFTFACE) /= GRID_UNDEFINED .and. &
                & geoType /= CGEO_TRIA_NOLEFT) &
                & region(ix,iy,1)=-6
            end if
            iy = iy - 1
          end do
        end do
! DN: We find the contact to the top left PFR wall: y-region 3
        do ix = leftcut(2), rightcut(2)-1
          iy = -1
          do while (iy.lt.topcut(2) .and. cflag(ix,iy,CELLFLAG_TYPE) /= GRID_INTERNAL)
            if (cflag(ix, iy, CELLFLAG_TYPE) == GRID_BOUNDARY) then
              Xvertices(:) = crx(ix,iy,:)
              Yvertices(:) = cry(ix,iy,:)
              geoType = cellGeoType(Xvertices, Yvertices)
              if (cflag(ix, iy, CELLFLAG_BOTTOMFACE) /= lefttargetindex(2) .and. &
                & cflag(ix, iy, CELLFLAG_BOTTOMFACE) /= righttargetindex(2) .and. &
                & cflag(ix, iy, CELLFLAG_BOTTOMFACE) /= GRID_UNDEFINED .and. &
                & geoType /= CGEO_TRIA_NOBOT) &
                & region(ix,iy,2)=3
              if (cflag(ix, iy, CELLFLAG_RIGHTFACE) /= lefttargetindex(2) .and. &
                & cflag(ix, iy, CELLFLAG_RIGHTFACE) /= righttargetindex(2) .and. &
                & cflag(ix, iy, CELLFLAG_RIGHTFACE) /= GRID_UNDEFINED .and. &
                & geoType /= CGEO_TRIA_NORIGHT) &
                & region(rightix(ix,iy),rightiy(ix,iy),1)=-3
              if (cflag(ix, iy, CELLFLAG_LEFTFACE) /= lefttargetindex(2) .and. &
                & cflag(ix, iy, CELLFLAG_LEFTFACE) /= righttargetindex(2) .and. &
                & cflag(ix, iy, CELLFLAG_LEFTFACE) /= GRID_UNDEFINED .and. &
                & geoType /= CGEO_TRIA_NOLEFT) &
                & region(ix,iy,1)=-3
            end if
            iy = iy + 1
          end do
          if (iy.eq.topcut(2)) exit
        end do
! DN: We find the contact to the top left divertor wall: y-region 7
        do ix = leftcut(2), rightcut(2)-1
          iy = ny
          do while (iy.ge.topcut(2) .and. cflag(ix,iy,CELLFLAG_TYPE) /= GRID_INTERNAL)
            if (cflag(ix, iy, CELLFLAG_TYPE) == GRID_BOUNDARY) then
              Xvertices(:) = crx(ix,iy,:)
              Yvertices(:) = cry(ix,iy,:)
              geoType = cellGeoType(Xvertices, Yvertices)
              if (cflag(ix, iy, CELLFLAG_TOPFACE) /= lefttargetindex(2) .and. &
                & cflag(ix, iy, CELLFLAG_TOPFACE) /= righttargetindex(2) .and. &
                & cflag(ix, iy, CELLFLAG_TOPFACE) /= GRID_UNDEFINED .and. &
                & geoType /= CGEO_TRIA_NOTOP) &
                & region(topix(ix,iy),topiy(ix,iy),2)=7
              if (cflag(ix, iy, CELLFLAG_RIGHTFACE) /= lefttargetindex(2) .and. &
                & cflag(ix, iy, CELLFLAG_RIGHTFACE) /= righttargetindex(2) .and. &
                & cflag(ix, iy, CELLFLAG_RIGHTFACE) /= GRID_UNDEFINED .and. &
                & geoType /= CGEO_TRIA_NORIGHT) &
                & region(rightix(ix,iy),rightiy(ix,iy),1)=-7
              if (cflag(ix, iy, CELLFLAG_LEFTFACE) /= lefttargetindex(2) .and. &
                & cflag(ix, iy, CELLFLAG_LEFTFACE) /= righttargetindex(2) .and. &
                & cflag(ix, iy, CELLFLAG_LEFTFACE) /= GRID_UNDEFINED .and. &
                & geoType /= CGEO_TRIA_NOLEFT) &
                & region(ix,iy,1)=-7
            end if
            iy = iy - 1
          end do
        end do
! DN:  We find the contact to the right top PFR wall: y-region 8
        do ix = rightcut(2)-1, leftcut(2), -1
          iy = -1
          do while (iy.lt.topcut(2) .and. cflag(ix,iy,CELLFLAG_TYPE) /= GRID_INTERNAL)
            if (cflag(ix, iy, CELLFLAG_TYPE) == GRID_BOUNDARY) then
              Xvertices(:) = crx(ix,iy,:)
              Yvertices(:) = cry(ix,iy,:)
              geoType = cellGeoType(Xvertices, Yvertices)
              if (cflag(ix, iy, CELLFLAG_BOTTOMFACE) /= lefttargetindex(2) .and. &
                & cflag(ix, iy, CELLFLAG_BOTTOMFACE) /= righttargetindex(2) .and. &
                & cflag(ix, iy, CELLFLAG_BOTTOMFACE) /= GRID_UNDEFINED .and. &
                & geoType /= CGEO_TRIA_NOBOT) &
                & region(ix,iy,2)=8
              if (cflag(ix, iy, CELLFLAG_RIGHTFACE) /= lefttargetindex(2) .and. &
                & cflag(ix, iy, CELLFLAG_RIGHTFACE) /= righttargetindex(2) .and. &
                & cflag(ix, iy, CELLFLAG_RIGHTFACE) /= GRID_UNDEFINED .and. &
                & geoType /= CGEO_TRIA_NORIGHT) &
                & region(rightix(ix,iy),rightiy(ix,iy),1)=-8
              if (cflag(ix, iy, CELLFLAG_LEFTFACE) /= lefttargetindex(2) .and. &
                & cflag(ix, iy, CELLFLAG_LEFTFACE) /= righttargetindex(2) .and. &
                & cflag(ix, iy, CELLFLAG_LEFTFACE) /= GRID_UNDEFINED .and. &
                & geoType /= CGEO_TRIA_NOLEFT) &
                & region(ix,iy,1)=-8
            end if
            iy = iy + 1
          end do
! We find the contact to the right top divertor wall: y-region 12
          iy = ny
          do while (iy.ge.topcut(2) .and. cflag(ix,iy,CELLFLAG_TYPE) /= GRID_INTERNAL)
            if (cflag(ix, iy, CELLFLAG_TYPE) == GRID_BOUNDARY) then
              Xvertices(:) = crx(ix,iy,:)
              Yvertices(:) = cry(ix,iy,:)
              geoType = cellGeoType(Xvertices, Yvertices)
              if (cflag(ix, iy, CELLFLAG_TOPFACE) /= lefttargetindex(2) .and. &
                & cflag(ix, iy, CELLFLAG_TOPFACE) /= righttargetindex(2) .and. &
                & cflag(ix, iy, CELLFLAG_TOPFACE) /= GRID_UNDEFINED .and. &
                & geoType /= CGEO_TRIA_NOTOP) &
                & region(topix(ix,iy),topiy(ix,iy),2)=12
              if (cflag(ix, iy, CELLFLAG_RIGHTFACE) /= lefttargetindex(2) .and. &
                & cflag(ix, iy, CELLFLAG_RIGHTFACE) /= righttargetindex(2) .and. &
                & cflag(ix, iy, CELLFLAG_RIGHTFACE) /= GRID_UNDEFINED .and. &
                & geoType /= CGEO_TRIA_NORIGHT) &
                & region(rightix(ix,iy),rightiy(ix,iy),1)=-12
              if (cflag(ix, iy, CELLFLAG_LEFTFACE) /= lefttargetindex(2) .and. &
                & cflag(ix, iy, CELLFLAG_LEFTFACE) /= righttargetindex(2) .and. &
                & cflag(ix, iy, CELLFLAG_LEFTFACE) /= GRID_UNDEFINED .and. &
                & geoType /= CGEO_TRIA_NOLEFT) &
                & region(ix,iy,1)=-12
            end if
            iy = iy - 1
          end do
        end do
        do ix = rightcut(2),rightcut(1)-1
! DN: We find the outer core boundary: y-region 9
          iy = -1
          do while (iy.le.topcut(1) .and. iy.le.topcut(2) .and. &
                  & .not.isRealCell(cflag(ix, iy, CELLFLAG_TYPE)))
            iy = iy+1
          end do
          if (iy.le.topcut(1) .and. iy.le.topcut(2)) then
            region(ix,iy,2)=9
          endif
! DN: We find the outer separatrix: y-region 11
          region(ix,min(topcut(1),topcut(2)),2)=11
! DN: We find the contact to the inner main chamber wall: y-region 13
          iy = ny
          do while (iy.ge.topcut(1) .and. iy.ge.topcut(2) .and. &
                  & cflag(ix,iy,CELLFLAG_TYPE) /= GRID_INTERNAL)
            if (cflag(ix, iy, CELLFLAG_TYPE) == GRID_BOUNDARY) then
              Xvertices(:) = crx(ix,iy,:)
              Yvertices(:) = cry(ix,iy,:)
              geoType = cellGeoType(Xvertices, Yvertices)
              if (cflag(ix, iy, CELLFLAG_TOPFACE) /= lefttargetindex(1) .and. &
                & cflag(ix, iy, CELLFLAG_TOPFACE) /= righttargetindex(1) .and. &
                & cflag(ix, iy, CELLFLAG_TOPFACE) /= lefttargetindex(2) .and. &
                & cflag(ix, iy, CELLFLAG_TOPFACE) /= righttargetindex(2) .and. &
                & cflag(ix, iy, CELLFLAG_TOPFACE) /= GRID_UNDEFINED .and. &
                & geoType /= CGEO_TRIA_NOTOP) &
                & region(topix(ix,iy),topiy(ix,iy),2)=13
              if (cflag(ix, iy, CELLFLAG_RIGHTFACE) /= lefttargetindex(1) .and. &
                & cflag(ix, iy, CELLFLAG_RIGHTFACE) /= righttargetindex(1) .and. &
                & cflag(ix, iy, CELLFLAG_RIGHTFACE) /= lefttargetindex(2) .and. &
                & cflag(ix, iy, CELLFLAG_RIGHTFACE) /= righttargetindex(2) .and. &
                & cflag(ix, iy, CELLFLAG_RIGHTFACE) /= GRID_UNDEFINED .and. &
                & geoType /= CGEO_TRIA_NORIGHT) &
                & region(rightix(ix,iy),rightiy(ix,iy),1)=-13
              if (cflag(ix, iy, CELLFLAG_LEFTFACE) /= lefttargetindex(1) .and. &
                & cflag(ix, iy, CELLFLAG_LEFTFACE) /= righttargetindex(1) .and. &
                & cflag(ix, iy, CELLFLAG_LEFTFACE) /= lefttargetindex(2) .and. &
                & cflag(ix, iy, CELLFLAG_LEFTFACE) /= righttargetindex(2) .and. &
                & cflag(ix, iy, CELLFLAG_LEFTFACE) /= GRID_UNDEFINED .and. &
                & geoType /= CGEO_TRIA_NOLEFT) &
                & region(ix,iy,1)=-13
            end if
            iy = iy - 1
          end do
        end do
! DN: We find the contact to the bottom right PFR wall: y-region 10
        do ix = rightcut(1), nx
          iy = -1
          do while (iy.le.topcut(1) .and. cflag(ix,iy,CELLFLAG_TYPE) /= GRID_INTERNAL)
            if (cflag(ix, iy, CELLFLAG_TYPE) == GRID_BOUNDARY) then
              Xvertices(:) = crx(ix,iy,:)
              Yvertices(:) = cry(ix,iy,:)
              geoType = cellGeoType(Xvertices, Yvertices)
              if (cflag(ix, iy, CELLFLAG_BOTTOMFACE) /= lefttargetindex(1) .and. &
                & cflag(ix, iy, CELLFLAG_BOTTOMFACE) /= righttargetindex(1) .and. &
                & cflag(ix, iy, CELLFLAG_BOTTOMFACE) /= GRID_UNDEFINED .and. &
                & geoType /= CGEO_TRIA_NOBOT) &
                & region(ix,iy,2)=10
              if (cflag(ix, iy, CELLFLAG_RIGHTFACE) /= lefttargetindex(1) .and. &
                & cflag(ix, iy, CELLFLAG_RIGHTFACE) /= righttargetindex(1) .and. &
                & cflag(ix, iy, CELLFLAG_RIGHTFACE) /= GRID_UNDEFINED .and. &
                & geoType /= CGEO_TRIA_NORIGHT) &
                & region(rightix(ix,iy),rightiy(ix,iy),1)=-10
              if (cflag(ix, iy, CELLFLAG_LEFTFACE) /= lefttargetindex(1) .and. &
                & cflag(ix, iy, CELLFLAG_LEFTFACE) /= righttargetindex(1) .and. &
                & cflag(ix, iy, CELLFLAG_LEFTFACE) /= GRID_UNDEFINED .and. &
                & geoType /= CGEO_TRIA_NOLEFT) &
                & region(ix,iy,1)=-10
            end if
            iy = iy + 1
          end do
! DN: We find the contact to the right divertor wall: y-region 14
          iy = ny
          do while (iy.gt.topcut(1) .and. cflag(ix,iy,CELLFLAG_TYPE) /= GRID_INTERNAL)
            if (cflag(ix, iy, CELLFLAG_TYPE) == GRID_BOUNDARY) then
              Xvertices(:) = crx(ix,iy,:)
              Yvertices(:) = cry(ix,iy,:)
              geoType = cellGeoType(Xvertices, Yvertices)
              if (cflag(ix, iy, CELLFLAG_TOPFACE) /= lefttargetindex(1) .and. &
                & cflag(ix, iy, CELLFLAG_TOPFACE) /= righttargetindex(1) .and. &
                & cflag(ix, iy, CELLFLAG_TOPFACE) /= GRID_UNDEFINED .and. &
                & geoType /= CGEO_TRIA_NOTOP) &
                & region(topix(ix,iy),topiy(ix,iy),2)=14
              if (cflag(ix, iy, CELLFLAG_RIGHTFACE) /= lefttargetindex(1) .and. &
                & cflag(ix, iy, CELLFLAG_RIGHTFACE) /= righttargetindex(1) .and. &
                & cflag(ix, iy, CELLFLAG_RIGHTFACE) /= GRID_UNDEFINED .and. &
                & geoType /= CGEO_TRIA_NORIGHT) &
                & region(rightix(ix,iy),rightiy(ix,iy),1)=-14
              if (cflag(ix, iy, CELLFLAG_LEFTFACE) /= lefttargetindex(1) .and. &
                & cflag(ix, iy, CELLFLAG_LEFTFACE) /= righttargetindex(1) .and. &
                & cflag(ix, iy, CELLFLAG_LEFTFACE) /= GRID_UNDEFINED .and. &
                & geoType /= CGEO_TRIA_NOLEFT) &
                & region(ix,iy,1)=-14
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
          Xvertices(:) = crx(ix,iy,:)
          Yvertices(:) = cry(ix,iy,:)
          geoType = cellGeoType(Xvertices, Yvertices)
          if (geoType /= CGEO_TRIA_NOLEFT .and. &
            & cflag(ix, iy, CELLFLAG_LEFTFACE) /= GRID_UNDEFINED) then
            lefttargetindex(1) = cflag(ix,iy,CELLFLAG_LEFTFACE)
          else if (geoType /= CGEO_TRIA_NOTOP .and. &
            & cflag(ix, iy, CELLFLAG_TOPFACE) /= GRID_UNDEFINED) then
            lefttargetindex(1) = cflag(ix,iy,CELLFLAG_TOPFACE)
          else if (geoType /= CGEO_TRIA_NOBOT .and. &
            & cflag(ix, iy, CELLFLAG_BOTTOMFACE) /= GRID_UNDEFINED) then
            lefttargetindex(1) = cflag(ix,iy,CELLFLAG_BOTTOMFACE)
          endif
        else
          ix = ix+1
        endif
      end do
#ifdef BUILDING_CARRE
      if (.not. lefttargetindex(1).ne.GRID_UNDEFINED) stop 'Left target not found!'
#else
      call xertst(lefttargetindex(1).ne.GRID_UNDEFINED, 'Left target not found!')
#endif
! We identify the structure index of the right target
      ix = nx
      righttargetindex(1)=GRID_UNDEFINED
      do while (ix.ge.-1 .and. righttargetindex(1).eq.GRID_UNDEFINED)
        if (cflag(ix, iy, CELLFLAG_TYPE) == GRID_BOUNDARY) then
          Xvertices(:) = crx(ix,iy,:)
          Yvertices(:) = cry(ix,iy,:)
          geoType = cellGeoType(Xvertices, Yvertices)
          if (geoType /= CGEO_TRIA_NORIGHT .and. &
            & cflag(ix, iy, CELLFLAG_RIGHTFACE) /= GRID_UNDEFINED) then
            righttargetindex(1) = cflag(ix,iy,CELLFLAG_RIGHTFACE)
          else if (geoType /= CGEO_TRIA_NOTOP .and. &
            & cflag(ix, iy, CELLFLAG_TOPFACE) /= GRID_UNDEFINED) then
            righttargetindex(1) = cflag(ix,iy,CELLFLAG_TOPFACE)
          else if (geoType /= CGEO_TRIA_NOBOT .and. &
            & cflag(ix, iy, CELLFLAG_BOTTOMFACE) /= GRID_UNDEFINED) then
            righttargetindex(1) = cflag(ix,iy,CELLFLAG_BOTTOMFACE)
          endif
        else
          ix = ix-1
        endif
      end do
#ifdef BUILDING_CARRE
      if (.not. righttargetindex(1).ne.GRID_UNDEFINED) stop 'Right target not found!'
#else
      call xertst(righttargetindex(1).ne.GRID_UNDEFINED, 'Right target not found!')
#endif
      do iy = -1, ny
! We identify the left edge
        ix = -1
        do while (ix.lt.nx .and. cflag(ix, iy, CELLFLAG_TYPE) /= GRID_INTERNAL)
          if (cflag(ix, iy, CELLFLAG_TYPE) == GRID_BOUNDARY) then
            Xvertices(:) = crx(ix,iy,:)
            Yvertices(:) = cry(ix,iy,:)
            geoType = cellGeoType(Xvertices, Yvertices)
            if (cflag(ix, iy, CELLFLAG_LEFTFACE) == lefttargetindex(1) .and. &
              & geoType /= CGEO_TRIA_NOLEFT) &
              & region(ix,iy,1)=1
            if (cflag(ix, iy, CELLFLAG_TOPFACE) == lefttargetindex(1) .and. &
              & geoType /= CGEO_TRIA_NOTOP) &
              & region(topix(ix,iy),topiy(ix,iy),2)=-1
            if (cflag(ix, iy, CELLFLAG_BOTTOMFACE) == lefttargetindex(1) .and. &
              & geoType /= CGEO_TRIA_NOBOT) &
              & region(ix,iy,2)=-1
          endif
          ix = ix + 1
        enddo
! We identify the right edge
        ix = nx
        do while (ix.gt.-1 .and. cflag(ix, iy, CELLFLAG_TYPE) /= GRID_INTERNAL)
          if (cflag(ix, iy, CELLFLAG_TYPE) == GRID_BOUNDARY) then
            Xvertices(:) = crx(ix,iy,:)
            Yvertices(:) = cry(ix,iy,:)
            geoType = cellGeoType(Xvertices, Yvertices)
            if (cflag(ix, iy, CELLFLAG_RIGHTFACE) == righttargetindex(1) .and. &
              & geoType /= CGEO_TRIA_NORIGHT) &
              & region(rightix(ix,iy),rightiy(ix,iy),1)=2
            if (cflag(ix, iy, CELLFLAG_TOPFACE) == righttargetindex(1) .and. &
              & geoType /= CGEO_TRIA_NOTOP) &
              & region(topix(ix,iy),topiy(ix,iy),2)=-2
            if (cflag(ix, iy, CELLFLAG_BOTTOMFACE) == righttargetindex(1) .and. &
              & geoType /= CGEO_TRIA_NOBOT) &
              & region(ix,iy,2)=-2
          endif
          ix = ix - 1
        enddo
      enddo
      nnreg(1)=2
! We identify the bottom edge
      do ix = -1, nx
        iy = -1
        do while (iy.lt.ny .and. cflag(ix,iy,CELLFLAG_TYPE) /= GRID_INTERNAL)
          if (cflag(ix, iy, CELLFLAG_TYPE) == GRID_BOUNDARY) then
            Xvertices(:) = crx(ix,iy,:)
            Yvertices(:) = cry(ix,iy,:)
            geoType = cellGeoType(Xvertices, Yvertices)
            if (cflag(ix, iy, CELLFLAG_BOTTOMFACE) /= GRID_UNDEFINED .and. &
              & cflag(ix, iy, CELLFLAG_BOTTOMFACE) /= lefttargetindex(1) .and. &
              & cflag(ix, iy, CELLFLAG_BOTTOMFACE) /= righttargetindex(1) .and. &
              & geoType /= CGEO_TRIA_NOBOT) &
              & region(ix,iy,2)=1
            if (cflag(ix, iy, CELLFLAG_RIGHTFACE) /= GRID_UNDEFINED .and. &
              & cflag(ix, iy, CELLFLAG_RIGHTFACE) /= lefttargetindex(1) .and. &
              & cflag(ix, iy, CELLFLAG_RIGHTFACE) /= righttargetindex(1) .and. &
              & geoType /= CGEO_TRIA_NORIGHT) &
              & region(rightix(ix,iy),rightiy(ix,iy),1)=-1
            if (cflag(ix, iy, CELLFLAG_LEFTFACE) /= GRID_UNDEFINED .and. &
              & cflag(ix, iy, CELLFLAG_LEFTFACE) /= lefttargetindex(1) .and. &
              & cflag(ix, iy, CELLFLAG_LEFTFACE) /= righttargetindex(1) .and. &
              & geoType /= CGEO_TRIA_NOLEFT) &
              & region(ix,iy,1)=-1
          end if
          iy = iy + 1
        end do
      end do
! We identify the top edge
      do ix = -1, nx
        iy = ny
        do while (iy.gt.-1 .and. cflag(ix,iy,CELLFLAG_TYPE) /= GRID_INTERNAL)
          if (cflag(ix, iy, CELLFLAG_TYPE) == GRID_BOUNDARY) then
            Xvertices(:) = crx(ix,iy,:)
            Yvertices(:) = cry(ix,iy,:)
            geoType = cellGeoType(Xvertices, Yvertices)
            if (cflag(ix, iy, CELLFLAG_TOPFACE) /= GRID_UNDEFINED .and. &
              & cflag(ix, iy, CELLFLAG_TOPFACE) /= lefttargetindex(1) .and. &
              & cflag(ix, iy, CELLFLAG_TOPFACE) /= righttargetindex(1) .and. &
              & geoType /= CGEO_TRIA_NOTOP) &
              & region(topix(ix,iy),topiy(ix,iy),2)=2
            if (cflag(ix, iy, CELLFLAG_RIGHTFACE) /= GRID_UNDEFINED .and. &
              & cflag(ix, iy, CELLFLAG_RIGHTFACE) /= lefttargetindex(1) .and. &
              & cflag(ix, iy, CELLFLAG_RIGHTFACE) /= righttargetindex(1) .and. &
              & geoType /= CGEO_TRIA_NORIGHT) &
              & region(rightix(ix,iy),rightiy(ix,iy),1)=-2
            if (cflag(ix, iy, CELLFLAG_LEFTFACE) /= GRID_UNDEFINED .and. &
              & cflag(ix, iy, CELLFLAG_LEFTFACE) /= lefttargetindex(1) .and. &
              & cflag(ix, iy, CELLFLAG_LEFTFACE) /= righttargetindex(1) .and. &
              & geoType /= CGEO_TRIA_NOLEFT) &
              & region(ix,iy,1)=-2
          end if
          iy = iy - 1
        end do
      end do
      nnreg(2)=2
    endif

    ! set up region where the residuals should be ignored (none yet)
    do iy=-1,ny
        do ix=-1,nx
            if(region(ix,iy,0).eq.0 .or. &
                &     isGhostCell(cflag(ix,iy,CELLFLAG_TYPE))) then
                resignore(ix,iy,1)=0
                resignore(ix,iy,2)=0
            else
                resignore(ix,iy,1)=1
                resignore(ix,iy,2)=1
            endif
        enddo
    enddo

    return

  contains

    function dist(ix1,iy1,ip1,ix2,iy2,ip2)
    implicit none
    integer, intent(in) :: ix1, iy1, ip1, ix2, iy2, ip2
    real (kind=R8) :: dist

    dist= sqrt((crx(ix1,iy1,ip1)-crx(ix2,iy2,ip2))**2+ &
        &      (cry(ix1,iy1,ip1)-cry(ix2,iy2,ip2))**2)
    return
    end function dist

    logical function match(ix1,iy1,ix2,iy2)
    implicit none
    integer, intent(in) :: ix1, iy1, ix2, iy2

    match = (dist(ix1,iy1,1,ix2,iy2,0)+dist(ix1,iy1,3,ix2,iy2,2)).lt. &
          &  geom_match_dist
    return
    end function match

  end subroutine init_region

  !> Cell categorizations

  !> Identify cells not used by the solver
  elemental logical function isUnusedCell(celltype)
    integer, intent(in) :: celltype

    isUnusedCell = (celltype == GRID_UNDEFINED) &
        & .or. (celltype == GRID_EXTERNAL) &
        & .or. (celltype == GRID_DEAD)
  end function isUnusedCell

  !> Identify cells not used by the solver
  elemental logical function isGhostCell(celltype)
    integer, intent(in) :: celltype

    isGhostCell = (celltype == GRID_GUARD)
  end function isGhostCell

  !> Identify cells not used by the solver
  elemental logical function isBoundaryCell(celltype)
    integer, intent(in) :: celltype

    isBoundaryCell = (celltype == GRID_BOUNDARY)
  end function isBoundaryCell

  !> Identify cells inside computational domain
  elemental logical function isRealCell(celltype)
    integer, intent(in) :: celltype

     isRealCell = (celltype == GRID_INTERNAL) .or. (celltype == GRID_BOUNDARY)
  end function isRealCell

  logical function isClassicalGrid(cflags)
    integer, dimension(:,:,:), intent(in) :: cflags

    isClassicalGrid = count( cflags(:, :, CELLFLAG_TYPE) == GRID_EXTERNAL ) == 0
  end function isClassicalGrid

  logical function isInDomain(nx, ny, ix, iy)
    integer, intent(in) :: nx, ny, ix, iy

    isInDomain = (ix >= -1 .and. ix <= nx &
          & .and. iy >= -1 .and. iy <= ny)

  end function isInDomain

!!$  !> Check if points (x1,y1) and (x2,y2) are identical
!!$  !> (i.e, very very close to each other)
!!$  logical function pointsIdentical( x1, y1, x2, y2, absTol )
!!$    real(R8), intent(in) :: x1, y1, x2, y2
!!$    real(R8), intent(in), optional :: absTol
!!$
!!$    ! internal
!!$    real(R8) :: lAbsTol
!!$
!!$    real(R8), parameter :: DEFAULTABSTOL = 1e-6
!!$
!!$    lAbsTol = DEFAULTABSTOL
!!$    if (present(absTol)) lAbsTol = absTol
!!$
!!$    pointsIdentical = ( points_dist(x1, y1, x2, y2) < lAbsTol )
!!$
!!$  end function pointsIdentical

end module b2mod_connectivity

!!!Local Variables:
!!! mode: f90
!!! End:
