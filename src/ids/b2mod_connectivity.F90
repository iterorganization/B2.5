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

    integer, parameter :: NO_CONNECTIVITY = huge(0) !< Constant to mark in
        !< connectivity arrays that no connectivity available

    !! Geometry/topology IDs (obtain using function geometryId(..:))

    integer, parameter :: GEOMETRY_COUNT = 13
        !< Number of different geometry/topology situations = max(GEOMETRY_*)

    !! The IDs, matching the IDS definitions of the GGD identifiers
    integer, parameter :: GEOMETRY_UNSPECIFIED = 0
    integer, parameter :: GEOMETRY_LINEAR = 1
    integer, parameter :: GEOMETRY_CYLINDER = 2
    integer, parameter :: GEOMETRY_LIMITER = 3
    integer, parameter :: GEOMETRY_SN = 4
    integer, parameter :: GEOMETRY_CDN = 5
    integer, parameter :: GEOMETRY_DDN_BOTTOM = 6
    integer, parameter :: GEOMETRY_DDN_TOP = 7
    integer, parameter :: GEOMETRY_ANNULUS = 8
    integer, parameter :: GEOMETRY_STELLARATORISLAND = 9
    integer, parameter :: GEOMETRY_STRUCTURED_SPACES = 10
    integer, parameter :: GEOMETRY_LFS_SNOWFLAKE_MINUS = 11
    integer, parameter :: GEOMETRY_LFS_SNOWFLAKE_PLUS = 12

    !! Region types
    !! Region type indices are the ones used in the B2 region array,
    !! i.e. zero-based.
    integer, parameter :: REGIONTYPE_COUNT = 3
        !< Number of different region types

    !! The types (indexing as in B2 region array, i.e. zero-based)
    integer, parameter :: REGIONTYPE_CELL = 0   !< First region type
    integer, parameter :: REGIONTYPE_XEDGE = 1  !< Second region type
    integer, parameter :: REGIONTYPE_YEDGE = 2  !< Third region type
    integer, parameter :: REGIONTYPE_EDGE = 3   !< Undifferentiated edge type

    !! Region counts and names

    !! Maximum number of regions of each type
    integer, parameter :: REGION_COUNT_MAX = 14

    !! Region counts
    !! First dimension: geometry type
    !! Second dimension: region type
    integer, dimension(0:REGIONTYPE_COUNT-1, 0:GEOMETRY_COUNT-1), parameter :: &
        &   regionCounts =  &
        &   reshape( (/     &
        &       1,  0,  0,  & !! GEOMETRY_UNSPECIFIED
        &       1,  2,  2,  & !! GEOMETRY_LINEAR
        &       1,  1,  2,  & !! GEOMETRY_CYLINDER
        &       2,  3,  3,  & !! GEOMETRY_LIMITER
        &       4,  6,  7,  & !! GEOMETRY_SN
        &       8, 12, 14,  & !! GEOMETRY_CDN
        &       8, 13, 14,  & !! GEOMETRY_DDN_BOTTOM
        &       8, 13, 14,  & !! GEOMETRY_DDN_TOP
        &       1,  2,  2,  & !! GEOMETRY_ANNULUS
        &       5,  7,  8,  & !! GEOMETRY_STELLARATORISLAND
        &       1,  0,  0,  & !! GEOMETRY_STRUCTURED_SPACES
        &       7, 13, 13,  & !! GEOMETRY_LFS_SNOWFLAKE_MINUS
        &       7, 13, 13   & !! GEOMETRY_LFS_SNOWFLAKE_PLUS
        &    /),            &
        &    (/ REGIONTYPE_COUNT, GEOMETRY_COUNT /) )   !< Region counts

    character(11), dimension(0:GEOMETRY_COUNT-1), parameter :: geometryName = &
        &   (/     &
        &    'UNSPECIFIED', &
        &    'LINEAR     ', &
        &    'CYLINDER   ', &
        &    'LIMITER    ', &
        &    'SN         ', &
        &    'CDN        ', &
        &    'DDN_BOTTOM ', &
        &    'DDN_TOP    ', &
        &    'ANNULUS    ', &
        &    'ISLAND     ', &
        &    'STRUCTURED ', &
        &    'LFS_SF-    ', &
        &    'LFS_SF+    '  &
        &   /)

    character(50), dimension(0:GEOMETRY_COUNT-1), parameter :: geometryDescription = &
        &   (/     &
        &    'Unspecified geometry                              ', &
        &    'Linear case                                       ', &
        &    'Cylinder geometry, straight in the third direction', &
        &    'Limiter geometry                                  ', &
        &    'Single null geometry                              ', &
        &    'Connected double null                             ', &
        &    'Disconnected double null, bottom X-point is active', &
        &    'Disconnected double null, top X-point is active   ', &
        &    'Annular geometry, curved in the third direction   ', &
        &    'Stellarator island geometry                       ', &
        &    'Structured grid built with multiple 1-D spaces    ', &
        &    'Low-field side Snowflake-minus geometry           ', &
        &    'Low-field side Snowflake-plus geometry            '  &
        &   /)

    !! Region names
    !! First dimension: geometry type (given in comments)
    !! Second dimension: region type
    !! Third dimension: region index

    character(32), parameter, private :: UU = repeat(' ', 32) !< UnUsed string

    character(32), dimension(REGION_COUNT_MAX, 0:REGIONTYPE_COUNT-1,                &
        &   0:GEOMETRY_COUNT-1) :: regionNames =                                    &
        &   reshape( (/                                                             &
        & & ! GEOMETRY_UNSPECIFIED
        &   'Plasma'//repeat(' ',26),                                               &
        &   UU, UU, UU, UU, UU, UU, UU, UU, UU, UU, UU, UU, UU,                     &
        & &
        &   'Left boundary                   ', 'Right boundary                  ', &
        &   UU, UU, UU, UU, UU, UU, UU, UU, UU, UU, UU, UU,                         &
        & &
        &   'Top boundary                    ', 'Bottom boundary                 ', &
        &   UU, UU, UU, UU, UU, UU, UU, UU, UU, UU, UU, UU,                         &
        & & ! GEOMETRY_LINEAR
        &   'Plasma'//repeat(' ',26),                                               &
        &   UU, UU, UU, UU, UU, UU, UU, UU, UU, UU, UU, UU, UU,                     &
        & &
        &   'Anti-clockwise boundary         ', 'Clockwise boundary              ', &
        &   UU, UU, UU, UU, UU, UU, UU, UU, UU, UU, UU, UU,                         &
        & &
        &   'Top boundary                    ', 'Bottom boundary                 ', &
        &   UU, UU, UU, UU, UU, UU, UU, UU, UU, UU, UU, UU,                         &
        & & ! GEOMETRY_CYLINDER
        &   'Plasma'//repeat(' ',26),                                               &
        &   UU, UU, UU, UU, UU, UU, UU, UU, UU, UU, UU, UU, UU,                     &
        & &
        &   'Periodicity boundary            ',                                     &
        &   UU, UU, UU, UU, UU, UU, UU, UU, UU, UU, UU, UU, UU,                     &
        & &
        &   'Top boundary                    ', 'Bottom boundary                 ', &
        &   UU, UU, UU, UU, UU, UU, UU, UU, UU, UU, UU, UU,                         &
        & & ! GEOMETRY_LIMITER
        &   'Core                            ', 'SOL                             ', &
        &   UU, UU, UU, UU, UU, UU, UU, UU, UU, UU, UU, UU,                         &
        & &
        &   'Anti-clockwise target           ', 'Clockwise target                ', &
        &   'Core cut                        ',                                     &
        & UU, UU, UU, UU, UU, UU, UU, UU, UU, UU, UU,                               &
        & &
        &   'Core boundary                   ', 'Separatrix                      ', &
        &   'Main chamber wall               ',                                     &
        &   UU, UU, UU, UU, UU, UU, UU, UU, UU, UU, UU,                             &
        & & ! GEOMETRY_SN Single null
        &   'Core                            ', 'SOL                             ', &
        &   'Western divertor                ', 'Eastern divertor                ', &
        &   UU, UU, UU, UU, UU, UU, UU, UU, UU, UU,                                 &
        & &
        &   'Western target                  ', 'Western throat                  ', &
        &   'Eastern throat                  ', 'Eastern target                  ', &
        &   'Core cut                        ', 'PFR cut                         ', &
        &   UU, UU, UU, UU, UU, UU, UU, UU,                                         &
        & &
        &   'Western PFR wall                ', 'Core boundary                   ', &
        &   'Eastern PFR wall                ', 'Separatrix                      ', &
        &   'Western baffle                  ', 'Main chamber wall               ', &
        &   'Eastern baffle                  ',                                     &
        &   UU, UU, UU, UU, UU, UU, UU,                                             &
        & & ! GEOMETRY_CDN Connected double null
        &   'Inner core                      ', 'Inner SOL                       ', &
        &   'Lower inner divertor            ', 'Upper inner divertor            ', &
        &   'Outer core                      ', 'Outer SOL                       ', &
        &   'Upper outer divertor            ', 'Lower outer divertor            ', &
        &   UU, UU, UU, UU, UU, UU,                                                 &
        & &
        &   'Lower inner target              ', 'Lower inner throat              ', &
        &   'Upper inner throat              ', 'Upper inner target              ', &
        &   'Upper outer target              ', 'Upper outer throat              ', &
        &   'Lower outer throat              ', 'Lower outer target              ', &
        &   'Lower core cut                  ', 'Upper PFR cut                   ', &
        &   'Upper core cut                  ', 'Lower PFR cut                   ', &
        &   UU, UU,                                                                 &
        & &
        &   'Lower inner PFR wall            ', 'Inner core boundary             ', &
        &   'Upper inner PFR wall            ', 'Inner separatrix                ', &
        &   'Lower inner baffle              ', 'Inner main chamber wall         ', &
        &   'Upper inner baffle              ', 'Upper outer PFR wall            ', &
        &   'Outer core boundary             ', 'Lower outer PFR wall            ', &
        &   'Outer separatrix                ', 'Upper outer baffle              ', &
        &   'Outer main chamber wall         ', 'Lower outer baffle              ', &
        & & ! GEOMETRY_DDN_BOTTOM
        &   'Inner core                      ', 'Inner SOL                       ', &
        &   'Lower inner divertor            ', 'Upper inner divertor            ', &
        &   'Outer core                      ', 'Outer SOL                       ', &
        &   'Upper outer divertor            ', 'Lower outer divertor            ', &
        &   UU, UU, UU, UU, UU, UU,                                                 &
        & &
        &   'Lower inner target              ', 'Lower inner throat              ', &
        &   'Upper inner throat              ', 'Upper inner target              ', &
        &   'Upper outer target              ', 'Upper outer throat              ', &
        &   'Lower outer throat              ', 'Lower outer target              ', &
        &   'Lower core cut                  ', 'Upper PFR cut                   ', &
        &   'Upper core cut                  ', 'Lower PFR cut                   ', &
        &   'Between separatrices core cut   ',                                     &
        &   UU,                                                                     &
        & &
        &   'Lower inner PFR wall            ', 'Inner core boundary             ', &
        &   'Upper inner PFR wall            ', 'Inner separatrix                ', &
        &   'Lower inner baffle              ', 'Inner main chamber wall         ', &
        &   'Upper inner baffle              ', 'Upper outer PFR wall            ', &
        &   'Outer core boundary             ', 'Lower outer PFR wall            ', &
        &   'Outer separatrix                ', 'Upper outer baffle              ', &
        &   'Outer main chamber wall         ', 'Lower outer baffle              ', &
        & & ! GEOMETRY_DDN_TOP
        &   'Inner core                      ', 'Inner SOL                       ', &
        &   'Lower inner divertor            ', 'Upper inner divertor            ', &
        &   'Outer core                      ', 'Outer SOL                       ', &
        &   'Upper outer divertor            ', 'Lower outer divertor            ', &
        &   UU, UU, UU, UU, UU, UU,                                                 &
        & &
        &   'Lower inner target              ', 'Lower inner throat              ', &
        &   'Upper inner throat              ', 'Upper inner target              ', &
        &   'Upper outer target              ', 'Upper outer throat              ', &
        &   'Lower outer throat              ', 'Lower outer target              ', &
        &   'Lower core cut                  ', 'Upper PFR cut                   ', &
        &   'Upper core cut                  ', 'Lower PFR cut                   ', &
        &   'Between separatrices core cut   ',                                     &
        &   UU,                                                                     &
        & &
        &   'Lower inner PFR wall            ', 'Inner core boundary             ', &
        &   'Upper inner PFR wall            ', 'Inner separatrix                ', &
        &   'Lower inner baffle              ', 'Inner main chamber wall         ', &
        &   'Upper inner baffle              ', 'Upper outer PFR wall            ', &
        &   'Outer core boundary             ', 'Lower outer PFR wall            ', &
        &   'Outer separatrix                ', 'Upper outer baffle              ', &
        &   'Outer main chamber wall         ', 'Lower outer baffle              ', &
        & & ! GEOMETRY_ANNULUS
        &   'Plasma'//repeat(' ',26),                                               &
        &   UU, UU, UU, UU, UU, UU, UU, UU, UU, UU, UU, UU, UU,                     &
        & &
        &   'Periodicity boundary            ',                                     &
        &   UU, UU, UU, UU, UU, UU, UU, UU, UU, UU, UU, UU, UU,                     &
        & &
        &   'Top boundary                    ', 'Bottom boundary                 ', &
        &   UU, UU, UU, UU, UU, UU, UU, UU, UU, UU, UU, UU,                         &
        & & ! GEOMETRY_STELLARATORISLAND
        &   'Core                            ', 'SOL                             ', &
        &   'Inner divertor                  ', 'Outer divertor                  ', &
        &   'Island                          ',                                     &
        &   UU, UU, UU, UU, UU, UU, UU, UU, UU,                                     &
        & &
        &   'Inner target                    ', 'Inner throat                    ', &
        &   'Outer throat                    ', 'Outer target                    ', &
        &   'Core cut                        ', 'PFR cut                         ', &
        &   'Island cut                      ',                                     &
        &   UU, UU, UU, UU, UU, UU, UU,                                             &
        & &
        &   'Inner PFR wall                  ', 'Core boundary                   ', &
        &   'Outer PFR wall                  ', 'Separatrix                      ', &
        &   'Entrance to inner PFR           ', 'Island center                   ', &
        &   'Entrance to outer PFR           ', 'Island boundary                 ', &
        &   UU, UU, UU, UU, UU, UU,                                                 &
        & & ! GEOMETRY_STRUCTURED_SPACES
        &   'Plasma'//repeat(' ',26),                                               &
        &   UU, UU, UU, UU, UU, UU, UU, UU, UU, UU, UU, UU, UU,                     &
        & &
        &   'Left boundary                   ', 'Right boundary                  ', &
        &   UU, UU, UU, UU, UU, UU, UU, UU, UU, UU, UU, UU,                         &
        & &
        &   'Top boundary                    ', 'Bottom boundary                 ', &
        &   UU, UU, UU, UU, UU, UU, UU, UU, UU, UU, UU, UU,                         &
        & & ! GEOMETRY_LFS_SNOWFLAKE_MINUS
        &   'Core                            ', 'SOL                             ', &
        &   'Inner divertor                  ', 'Outer divertor entrance         ', &
        &   'First outboard divertor leg     ', 'Second outboard divertor leg    ', &
        &   'Third outboard divertor leg     ',                                     &
        &   UU, UU, UU, UU, UU, UU, UU,                                             &
        & &
        &   'Inner target                    ', 'Inner divertor entrance         ', &
        &   'Outer divertor entrance         ', 'First outer divertor target     ', &
        &   'Second outer divertor target    ', 'SOL connection to outer leg 1   ', &
        &   'Connection of outer leg 2 and 3 ', 'Third outer divertor target     ', &
        &   'Core cut                        ', 'PFR cut                         ', &
        &   'Connection of outer leg 1 and 2 ', 'PFR connection to outer leg 3   ', &
        &   'CFR connection to outer leg 3   ',                                     &
        &   UU,                                                                     &
        & &
        &   'Inner PFR wall                  ', 'Core boundary                   ', &
        &   'Outer entrance PFR wall         ', 'Separatrix                      ', &
        &   'Inner baffle                    ', 'Main chamber wall               ', &
        &   'Outer entrance baffle           ', 'First outer leg PFR wall        ', &
        &   'Second outer leg PFR wall       ', 'Third outer leg PFR wall        ', &
        &   'First outer leg baffle          ', 'Second outer leg baffle         ', &
        &   'Third outer leg baffle          ',                                     &
        &   UU,                                                                     &
        & & ! GEOMETRY_LFS_SNOWFLAKE_PLUS
        &   'Core                            ', 'SOL                             ', &
        &   'Inner divertor                  ', 'Outer divertor entrance         ', &
        &   'First outboard divertor leg     ', 'Second outboard divertor leg    ', &
        &   'Third outboard divertor leg     ',                                     &
        &   UU, UU, UU, UU, UU, UU, UU,                                             &
        & &
        &   'Inner target                    ', 'Inner divertor entrance         ', &
        &   'Outer divertor entrance         ', 'First outer divertor target     ', &
        &   'Second outer divertor target    ', 'SOL connection to outer leg 1   ', &
        &   'Connection of outer leg 2 and 3 ', 'Third outer divertor target     ', &
        &   'Core cut                        ', 'PFR cut                         ', &
        &   'Connection of outer leg 1 and 2 ', 'PFR connection to outer leg 3   ', &
        &   'CFR connection to outer leg 3   ',                                     &
        &   UU,                                                                     &
        & &
        &   'Inner PFR wall                  ', 'Core boundary                   ', &
        &   'Outer entrance PFR wall         ', 'Separatrix                      ', &
        &   'Inner baffle                    ', 'Main chamber wall               ', &
        &   'Outer entrance baffle           ', 'First outer leg PFR wall        ', &
        &   'Second outer leg PFR wall       ', 'Third outer leg PFR wall        ', &
        &   'First outer leg baffle          ', 'Second outer leg baffle         ', &
        &   'Third outer leg baffle          ',                                     &
        &   UU                                                                      &
        & &
        &  /), &
        & (/REGION_COUNT_MAX, REGIONTYPE_COUNT, GEOMETRY_COUNT/) )

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
    integer :: ix, iy
    logical :: rightFace, leftFace, topFace, botFace !! has ... face
    logical :: rightNb, leftNb, topNb, botNb !! has ... neighbour
    logical :: error, thisCellError !! error occurred

    error = .false.

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
                write(*,'(a,4f12.6)') 'crx = ',crx(ix,iy,0:3)
                write(*,'(a,4f12.6)') 'cry = ',cry(ix,iy,0:3)
                error = .true.
            end if

        end do
    end do

    if (error) stop "test_connectivity: error(s) found"

  end subroutine test_connectivity


  !> Identify what geometry/topology is present from cut and periodicity data.
  !> Returns one of the GEOMETRY_* constants. Stops if unknown geometry.
  integer function geometryId( nnreg, isymm, periodic_bc, topcut )
    integer, intent(in) :: nnreg(0:2), isymm, periodic_bc, topcut(:)
    logical, save :: first
    data first/.true./

    if (nnreg(0) == 1 .and. periodic_bc.le.0) then
        geometryId = GEOMETRY_LINEAR
        if (first) then
            call logmsg( LOGDEBUG, "b2mod_connectivity.geometryId(): identified GEOMETRY_LINEAR")
            first = .false.
        end if
        return
    end if

    if (nnreg(0) == 1 .and. periodic_bc == 1) then
      if (isymm == 0) then
        geometryId = GEOMETRY_CYLINDER
        if (first) then
            call logmsg( LOGDEBUG, "b2mod_connectivity.geometryId(): identified GEOMETRY_CYLINDER")
            first = .false.
        end if
        return
      else
        geometryId = GEOMETRY_ANNULUS
        if (first) then
            call logmsg( LOGDEBUG, "b2mod_connectivity.geometryId(): identified GEOMETRY_ANNULUS")
            first = .false.
        end if
        return
      end if
    end if

    if (nnreg(0) == 2) then
        geometryId = GEOMETRY_LIMITER
        if (first) then
            call logmsg( LOGDEBUG, "b2mod_connectivity.geometryId(): identified GEOMETRY_LIMITER")
            first = .false.
        end if
        return
    end if

    if (nnreg(0) == 4) then
        geometryId = GEOMETRY_SN
        if (first) then
            call logmsg( LOGDEBUG, "b2mod_connectivity.geometryId(): identified GEOMETRY_SN")
            first = .false.
        end if
        return
    end if

    if (nnreg(0) == 5) then
        geometryId = GEOMETRY_STELLARATORISLAND
        if (first) then
            call logmsg( LOGDEBUG, "b2mod_connectivity.geometryId(): identified GEOMETRY_STELLARATORISLAND")
            first = .false.
        end if
        return
    end if

    if (nnreg(0) == 7) then
        if (topcut(1) < topcut(2)) then
            geometryId = GEOMETRY_LFS_SNOWFLAKE_MINUS
            if (first) then
                call logmsg( LOGDEBUG, "b2mod_connectivity.geometryId(): identified GEOMETRY_LFS_SNOWFLAKE_MINUS")
                first = .false.
            end if
            return
        end if
        if (topcut(1) > topcut(2)) then
            geometryId = GEOMETRY_LFS_SNOWFLAKE_PLUS
            if (first) then
                call logmsg( LOGDEBUG, "b2mod_connectivity.geometryId(): identified GEOMETRY_LFS_SNOWFLAKE_PLUS")
                first = .false.
            end if
            return
        end if
        if (topcut(1) == topcut(2)) then
            geometryId = GEOMETRY_UNSPECIFIED
            if (first) then
                call logmsg( LOGDEBUG, "b2mod_connectivity.geometryId(): unknown GEOMETRY_UNSPECIFIED")
                first = .false.
            end if
            return
        end if
    end if

    if (nnreg(0) == 8) then

        if (topcut(1) == topcut(2)) then
            geometryId = GEOMETRY_CDN
            if (first) then
                call logmsg( LOGDEBUG, "b2mod_connectivity.geometryId(): identified GEOMETRY_CDN")
                first = .false.
            end if
            return
        end if
        if (topcut(1) < topcut(2)) then
            geometryId = GEOMETRY_DDN_BOTTOM
            if (first) then
                call logmsg( LOGDEBUG, "b2mod_connectivity.geometryId(): identified GEOMETRY_DDN_BOTTOM")
                first = .false.
            end if
            return
        end if
        if (topcut(1) > topcut(2)) then
            geometryId = GEOMETRY_DDN_TOP
            if (first) then
                call logmsg( LOGDEBUG, "b2mod_connectivity.geometryId(): identified GEOMETRY_DDN_TOP")
                first = .false.
            end if
            return
        end if
    end if

    geometryId = GEOMETRY_UNSPECIFIED
    if (first) then
        call logmsg( LOGDEBUG, "b2mod_connectivity.geometryId(): unknown GEOMETRY_UNSPECIFIED")
        first = .false.
    end if

    stop 'b2mod_connectivity.geometryId: unknown geometry'

  end function geometryId


  !> Return number of regions of a given type for a given geometry.
  !> @param geometryId: see definition of GEOMETRY_* constants.
  !> @param regionType: see definition of REGIONTYPE_* constants.
  integer function regionCount( geometryId, regionType )
    integer, intent(in) :: geometryId, regionType

    !! TODO: maybe do some input parameter checking
    regionCount = regionCounts(regionType, geometryId)
  end function regionCount


  !> Return total number of regions (both cell and face regions) for
  !> the given geometry
  integer function regionCountTotal( geometryId )
    integer, intent(in) :: geometryId

    !! internal
    integer :: iType

    regionCountTotal = 0
    do iType = 0, REGIONTYPE_COUNT - 1
        regionCountTotal = regionCountTotal + regionCount( geometryId, iType )
    end do
  end function regionCountTotal


  !> Return name of a specific region.
  !> @param geometryId: see definition of GEOMETRY_* constants.
  !> @param regionType: see definition of REGIONTYPE_* constants.
  !> @param regionId: number of region. Should be between [1, regionCount(geometryId, regionType)].
  character(32) function regionName( geometryId, regionType, regionId )
    integer, intent(in) :: geometryId, regionType, regionId

    regionName = regionNames(regionId, regionType, geometryId)
  end function regionName


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
