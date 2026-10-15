module b2mod_geometry
    use logging &
     & , only: logmsg, LOGDEBUG, LOGWARNING
    use b2mod_cellhelper &
     & , only: LEFT, RIGHT, TOP, BOTTOM, NODIRECTION

    implicit none

    !! Geometry/topology IDs (obtain using function geometryId)

    integer, parameter :: GEOMETRY_COUNT = 14
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
    integer, parameter :: GEOMETRY_GENERAL = 13

    !! Region types
    !! Region type indices are the ones used in the B2 region array,
    !! i.e. zero-based.
    integer, parameter :: REGIONTYPE_COUNT = 4
        !< Number of different region types

    !! The types (indexing as in B2 region array, i.e. zero-based)
    integer, parameter :: REGIONTYPE_CELL = 0   !< First region type
    integer, parameter :: REGIONTYPE_XEDGE = 1  !< Second region type
    integer, parameter :: REGIONTYPE_YEDGE = 2  !< Third region type
    integer, parameter :: REGIONTYPE_EDGE = 3   !< Undifferentiated edge type

    !! Boundary direction indices
    integer, parameter :: WEST = LEFT, SOUTH = BOTTOM, EAST = RIGHT, NORTH = TOP

    !! Region counts and names

    !! Maximum number of regions of each type
    integer, parameter :: REGION_COUNT_MAX = 15

    !! Maximum number of regions over all geometry types
    integer, parameter :: REGION_NUMBER_MAX = 27

    !! Region counts
    !! First dimension: geometry type
    !! Second dimension: region type
    integer, dimension(0:REGIONTYPE_COUNT-1, 0:GEOMETRY_COUNT-1), parameter :: &
        &   regionCounts =     &
        &   reshape( (/        &
        &       1,  0,  0,  0, & !! GEOMETRY_UNSPECIFIED
        &       1,  0,  0,  4, & !! GEOMETRY_LINEAR
        &       1,  1,  2,  0, & !! GEOMETRY_CYLINDER
        &       2,  1,  2,  3, & !! GEOMETRY_LIMITER
        &       4,  4,  2,  7, & !! GEOMETRY_SN
        &       8,  8,  4, 14, & !! GEOMETRY_CDN
        &       8,  9,  4, 14, & !! GEOMETRY_DDN_BOTTOM
        &       8,  9,  4, 14, & !! GEOMETRY_DDN_TOP
        &       1,  2,  2,  0, & !! GEOMETRY_ANNULUS
        &       5,  5,  6,  4, & !! GEOMETRY_STELLARATORISLAND
        &       1,  0,  0,  0, & !! GEOMETRY_STRUCTURED_SPACES
        &       7,  9,  2, 15, & !! GEOMETRY_LFS_SNOWFLAKE_MINUS
        &       7,  9,  2, 15, & !! GEOMETRY_LFS_SNOWFLAKE_PLUS
        &       0,  0,  0,  0  & !! GEOMETRY_GENERAL
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
        &    'LFS_SF+    ', &
        &    'GENERAL    '  &
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
        &    'Low-field side Snowflake-plus geometry            ', &
        &    'General geometry and topology                     '  &
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
        &   UU, UU, UU, UU, UU, UU, UU, UU, UU, UU, UU, UU, UU, UU,                 &
        & &
        &   UU, UU, UU, UU, UU, UU, UU, UU, UU, UU, UU, UU, UU, UU, UU,             &
        & &
        &   UU, UU, UU, UU, UU, UU, UU, UU, UU, UU, UU, UU, UU, UU, UU,             &
        & &
        &   UU, UU, UU, UU, UU, UU, UU, UU, UU, UU, UU, UU, UU, UU, UU,             &
        & & ! GEOMETRY_LINEAR
        &   'Plasma'//repeat(' ',26),                                               &
        &   UU, UU, UU, UU, UU, UU, UU, UU, UU, UU, UU, UU, UU, UU,                 &
        & &
        &   UU, UU, UU, UU, UU, UU, UU, UU, UU, UU, UU, UU, UU, UU, UU,             &
        & &
        &   UU, UU, UU, UU, UU, UU, UU, UU, UU, UU, UU, UU, UU, UU, UU,             &
        & &
        &   'Anti-clockwise boundary         ', 'Clockwise boundary              ', &
        &   'Top boundary                    ', 'Bottom boundary                 ', &
        &   UU, UU, UU, UU, UU, UU, UU, UU, UU, UU, UU,                             &
        & & ! GEOMETRY_CYLINDER
        &   'Plasma'//repeat(' ',26),                                               &
        &   UU, UU, UU, UU, UU, UU, UU, UU, UU, UU, UU, UU, UU, UU,                 &
        & &
        &   'Periodicity boundary            ',                                     &
        &   UU, UU, UU, UU, UU, UU, UU, UU, UU, UU, UU, UU, UU, UU,                 &
        & &
        &   'Top boundary                    ', 'Bottom boundary                 ', &
        &   UU, UU, UU, UU, UU, UU, UU, UU, UU, UU, UU, UU, UU,                     &
        & &
        &   UU, UU, UU, UU, UU, UU, UU, UU, UU, UU, UU, UU, UU, UU, UU,             &
        & & ! GEOMETRY_LIMITER
        &   'Core                            ', 'Scrape-Off Layer                ', &
        &   UU, UU, UU, UU, UU, UU, UU, UU, UU, UU, UU, UU, UU,                     &
        & &
        &   'Core cut                        ',                                     &
        &   UU, UU, UU, UU, UU, UU, UU, UU, UU, UU, UU, UU, UU, UU,                 &
        & &
        &   'Core boundary                   ', 'Separatrix                      ', &
        &   UU, UU, UU, UU, UU, UU, UU, UU, UU, UU, UU, UU, UU,                     &
        & &
        &   'Anti-clockwise target           ', 'Clockwise target                ', &
        &   'Main chamber wall               ',                                     &
        &   UU, UU, UU, UU, UU, UU, UU, UU, UU, UU, UU, UU,                         &
        & & ! GEOMETRY_SN Single null
        &   'Core                            ', 'Scrape-Off Layer                ', &
        &   'Western divertor                ', 'Eastern divertor                ', &
        &   UU, UU, UU, UU, UU, UU, UU, UU, UU, UU, UU,                             &
        & &
        &   'Western throat                  ', 'Eastern throat                  ', &
        &   'Core cut                        ', 'PFR cut                         ', &
        &   UU, UU, UU, UU, UU, UU, UU, UU, UU, UU, UU,                             &
        & &
        &   'Core boundary                   ', 'Separatrix                      ', &
        &   UU, UU, UU, UU, UU, UU, UU, UU, UU, UU, UU, UU, UU,                     &
        & &
        &   'Western target                  ', 'Eastern target                  ', &
        &   'Western PFR wall                ', 'Eastern PFR wall                ', &
        &   'Western baffle                  ', 'Main chamber wall               ', &
        &   'Eastern baffle                  ',                                     &
        &   UU, UU, UU, UU, UU, UU, UU, UU,                                         &
        & & ! GEOMETRY_CDN Connected double null
        &   'Inner core                      ', 'Inner SOL                       ', &
        &   'Lower inner divertor            ', 'Upper inner divertor            ', &
        &   'Outer core                      ', 'Outer SOL                       ', &
        &   'Upper outer divertor            ', 'Lower outer divertor            ', &
        &   UU, UU, UU, UU, UU, UU, UU,                                             &
        & &
        &   'Lower inner throat              ', 'Upper inner throat              ', &
        &   'Upper outer throat              ', 'Lower outer throat              ', &
        &   'Lower core cut                  ', 'Upper PFR cut                   ', &
        &   'Upper core cut                  ', 'Lower PFR cut                   ', &
        &   UU, UU, UU, UU, UU, UU, UU,                                             &
        & &
        &   'Inner core boundary             ', 'Inner separatrix                ', &
        &   'Outer core boundary             ', 'Outer separatrix                ', &
        &   UU, UU, UU, UU, UU, UU, UU, UU, UU, UU, UU,                             &
        & &
        &   'Lower inner target              ', 'Upper inner target              ', &
        &   'Upper outer target              ', 'Lower outer target              ', &
        &   'Lower inner PFR wall            ', 'Upper inner PFR wall            ', &
        &   'Lower inner baffle              ', 'Inner main chamber wall         ', &
        &   'Upper inner baffle              ', 'Upper outer PFR wall            ', &
        &   'Lower outer PFR wall            ', 'Upper outer baffle              ', &
        &   'Outer main chamber wall         ', 'Lower outer baffle              ', &
        &   UU,                                                                     &
        & & ! GEOMETRY_DDN_BOTTOM
        &   'Inner core                      ', 'Inner SOL                       ', &
        &   'Lower inner divertor            ', 'Upper inner divertor            ', &
        &   'Outer core                      ', 'Outer SOL                       ', &
        &   'Upper outer divertor            ', 'Lower outer divertor            ', &
        &   UU, UU, UU, UU, UU, UU, UU,                                             &
        & &
        &   'Lower inner throat              ', 'Upper inner throat              ', &
        &   'Upper outer throat              ', 'Lower outer throat              ', &
        &   'Lower core cut                  ', 'Upper PFR cut                   ', &
        &   'Upper core cut                  ', 'Lower PFR cut                   ', &
        &   'Between separatrices SOL cut    ',                                     &
        &   UU, UU, UU, UU, UU, UU,                                                 &
        & &
        &   'Inner core boundary             ', 'Inner separatrix                ', &
        &   'Outer core boundary             ', 'Outer separatrix                ', &
        &   UU, UU, UU, UU, UU, UU, UU, UU, UU, UU, UU,                             &
        & &
        &   'Lower inner target              ', 'Upper inner target              ', &
        &   'Upper outer target              ', 'Lower outer target              ', &
        &   'Lower inner PFR wall            ', 'Upper inner PFR wall            ', &
        &   'Lower inner baffle              ', 'Inner main chamber wall         ', &
        &   'Upper inner baffle              ', 'Upper outer PFR wall            ', &
        &   'Lower outer PFR wall            ', 'Upper outer baffle              ', &
        &   'Outer main chamber wall         ', 'Lower outer baffle              ', &
        &   UU,                                                                     &
        & & ! GEOMETRY_DDN_TOP
        &   'Inner core                      ', 'Inner SOL                       ', &
        &   'Lower inner divertor            ', 'Upper inner divertor            ', &
        &   'Outer core                      ', 'Outer SOL                       ', &
        &   'Upper outer divertor            ', 'Lower outer divertor            ', &
        &   UU, UU, UU, UU, UU, UU, UU,                                             &
        & &
        &   'Lower inner throat              ', 'Upper inner throat              ', &
        &   'Upper outer throat              ', 'Lower outer throat              ', &
        &   'Lower core cut                  ', 'Upper PFR cut                   ', &
        &   'Upper core cut                  ', 'Lower PFR cut                   ', &
        &   'Between separatrices SOL cut    ',                                     &
        &   UU, UU, UU, UU, UU, UU,                                                 &
        & &
        &   'Inner core boundary             ', 'Inner separatrix                ', &
        &   'Outer core boundary             ', 'Outer separatrix                ', &
        &   UU, UU, UU, UU, UU, UU, UU, UU, UU, UU, UU,                             &
        & &
        &   'Lower inner target              ', 'Upper inner target              ', &
        &   'Upper outer target              ', 'Lower outer target              ', &
        &   'Lower inner PFR wall            ', 'Upper inner PFR wall            ', &
        &   'Lower inner baffle              ', 'Inner main chamber wall         ', &
        &   'Upper inner baffle              ', 'Upper outer PFR wall            ', &
        &   'Lower outer PFR wall            ', 'Upper outer baffle              ', &
        &   'Outer main chamber wall         ', 'Lower outer baffle              ', &
        &   UU,                                                                     &
        & & ! GEOMETRY_ANNULUS
        &   'Plasma'//repeat(' ',26),                                               &
        &   UU, UU, UU, UU, UU, UU, UU, UU, UU, UU, UU, UU, UU, UU,                 &
        & &
        &   'Periodicity boundary            ',                                     &
        &   UU, UU, UU, UU, UU, UU, UU, UU, UU, UU, UU, UU, UU, UU,                 &
        & &
        &   'Top boundary                    ', 'Bottom boundary                 ', &
        &   UU, UU, UU, UU, UU, UU, UU, UU, UU, UU, UU, UU, UU,                     &
        & &
        &   UU, UU, UU, UU, UU, UU, UU, UU, UU, UU, UU, UU, UU, UU, UU,             &
        & & ! GEOMETRY_STELLARATORISLAND
        &   'Core                            ', 'Scrape-Off Layer                ', &
        &   'Inner divertor                  ', 'Outer divertor                  ', &
        &   'Island                          ',                                     &
        &   UU, UU, UU, UU, UU, UU, UU, UU, UU, UU,                                 &
        & &
        &   'Inner throat                    ', 'Outer throat                    ', &
        &   'Core cut                        ', 'PFR cut                         ', &
        &   'Island cut                      ',                                     &
        &   UU, UU, UU, UU, UU, UU, UU, UU, UU, UU,                                 &
        & &
        &   'Core boundary                   ', 'Separatrix                      ', &
        &   'Entrance to inner PFR           ', 'Island center                   ', &
        &   'Entrance to outer PFR           ', 'Island boundary                 ', &
        &   UU, UU, UU, UU, UU, UU, UU, UU, UU,                                     &
        & &
        &   'Inner target                    ', 'Outer target                    ', &
        &   'Inner PFR wall                  ', 'Outer PFR wall                  ', &
        &   UU, UU, UU, UU, UU, UU, UU, UU, UU, UU, UU,                             &
        & & ! GEOMETRY_STRUCTURED_SPACES
        &   'Plasma'//repeat(' ',26),                                               &
        &   UU, UU, UU, UU, UU, UU, UU, UU, UU, UU, UU, UU, UU, UU,                 &
        & &
        &   'Left boundary                   ', 'Right boundary                  ', &
        &   UU, UU, UU, UU, UU, UU, UU, UU, UU, UU, UU, UU, UU,                     &
        & &
        &   'Top boundary                    ', 'Bottom boundary                 ', &
        &   UU, UU, UU, UU, UU, UU, UU, UU, UU, UU, UU, UU, UU,                     &
        & &
        &   UU, UU, UU, UU, UU, UU, UU, UU, UU, UU, UU, UU, UU, UU, UU,             &
        & & ! GEOMETRY_LFS_SNOWFLAKE_MINUS
        &   'Core                            ', 'Scrape-Off Layer                ', &
        &   'Inner divertor                  ', 'Outer divertor entrance         ', &
        &   'First outboard divertor leg     ', 'Second outboard divertor leg    ', &
        &   'Third outboard divertor leg     ',                                     &
        &   UU, UU, UU, UU, UU, UU, UU, UU,                                         &
        & &
        &   'Inner divertor entrance         ', 'Outer divertor entrance         ', &
        &   'SOL connection to outer leg 1   ', 'Connection of outer leg 2 and 3 ', &
        &   'Core cut                        ', 'PFR cut                         ', &
        &   'Connection of outer leg 1 and 2 ', 'PFR connection to outer leg 3   ', &
        &   'CFR connection to outer leg 3   ',                                     &
        &   UU, UU, UU, UU, UU, UU,                                                 &
        & &
        &   'Core boundary                   ', 'Separatrix                      ', &
        &   UU, UU, UU, UU, UU, UU, UU, UU, UU, UU, UU, UU, UU,                     &
        & &
        &   'Inner target                    ', 'First outer divertor target     ', &
        &   'Second outer divertor target    ', 'Third outer divertor target     ', &
        &   'Inner PFR wall                  ', 'Outer entrance PFR wall         ', &
        &   'Inner baffle                    ', 'Main chamber wall               ', &
        &   'Outer entrance baffle           ', 'First outer leg PFR wall        ', &
        &   'Second outer leg PFR wall       ', 'Third outer leg PFR wall        ', &
        &   'First outer leg baffle          ', 'Second outer leg baffle         ', &
        &   'Third outer leg baffle          ',                                     &
        & & ! GEOMETRY_LFS_SNOWFLAKE_PLUS
        &   'Core                            ', 'Scrape-Off Layer                ', &
        &   'Inner divertor                  ', 'Outer divertor entrance         ', &
        &   'First outboard divertor leg     ', 'Second outboard divertor leg    ', &
        &   'Third outboard divertor leg     ',                                     &
        &   UU, UU, UU, UU, UU, UU, UU, UU,                                         &
        & &
        &   'Inner divertor entrance         ', 'Outer divertor entrance         ', &
        &   'SOL connection to outer leg 1   ', 'Connection of outer leg 2 and 3 ', &
        &   'Core cut                        ', 'PFR cut                         ', &
        &   'Connection of outer leg 1 and 2 ', 'PFR connection to outer leg 3   ', &
        &   'CFR connection to outer leg 3   ',                                     &
        &   UU, UU, UU, UU, UU, UU,                                                 &
        & &
        &   'Core boundary                   ', 'Separatrix                      ', &
        &   UU, UU, UU, UU, UU, UU, UU, UU, UU, UU, UU, UU, UU,                     &
        & &
        &   'Inner target                    ', 'First outer divertor target     ', &
        &   'Second outer divertor target    ', 'Third outer divertor target     ', &
        &   'Inner PFR wall                  ', 'Outer entrance PFR wall         ', &
        &   'Inner baffle                    ', 'Main chamber wall               ', &
        &   'Outer entrance baffle           ', 'First outer leg PFR wall        ', &
        &   'Second outer leg PFR wall       ', 'Third outer leg PFR wall        ', &
        &   'First outer leg baffle          ', 'Second outer leg baffle         ', &
        &   'Third outer leg baffle          ',                                     &
        & & ! GEOMETRY_GENERAL
        &   'Plasma'//repeat(' ',26),                                               &
        &   UU, UU, UU, UU, UU, UU, UU, UU, UU, UU, UU, UU, UU, UU,                 &
        & &
        &   UU, UU, UU, UU, UU, UU, UU, UU, UU, UU, UU, UU, UU, UU, UU,             &
        & &
        &   UU, UU, UU, UU, UU, UU, UU, UU, UU, UU, UU, UU, UU, UU, UU,             &
        & &
        &   UU, UU, UU, UU, UU, UU, UU, UU, UU, UU, UU, UU, UU, UU, UU              &
        & &
        &  /), &
        & (/REGION_COUNT_MAX, REGIONTYPE_COUNT, GEOMETRY_COUNT/) )

   integer, dimension(REGION_COUNT_MAX, 0:REGIONTYPE_COUNT-1,           &
        &   0:GEOMETRY_COUNT-1) :: regionNumbers =                      &
        &   reshape( (/                                                 &
        & & ! GEOMETRY_UNSPECIFIED
        &    1,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0, &
        &    0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0, &
        &    0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0, &
        &    0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0, &
        & & ! GEOMETRY_LINEAR
        &    1,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0, &
        &    0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0, &
        &    0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0, &
        &    1,  2,  3,  4,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0, &
        & & ! GEOMETRY_CYLINDER
        &    1,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0, &
        &    1,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0, &
        &    2,  3,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0, &
        &    0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0, &
        & & ! GEOMETRY_LIMITER
        &    1,  2,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0, &
        &    3,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0, &
        &    4,  5,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0, &
        &    1,  2,  6,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0, &
        & & ! GEOMETRY_SN
        &    1,  2,  3,  4,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0, &
        &    2,  3,  5,  6,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0, &
        &    8, 10,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0, &
        &    1,  4,  7,  9, 11, 12, 13,  0,  0,  0,  0,  0,  0,  0,  0, &
        & & ! GEOMETRY_CDN
        &    1,  2,  3,  4,  5,  6,  7,  8,  0,  0,  0,  0,  0,  0,  0, &
        &    2,  3,  6,  7,  9, 10, 11, 12,  0,  0,  0,  0,  0,  0,  0, &
        &   14, 16, 21, 23,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0, &
        &    1,  4,  5,  8, 13, 15, 17, 18, 19, 20, 22, 24, 25, 26,  0, &
        & & ! GEOMETRY_DDN_BOTTOM
        &    1,  2,  3,  4,  5,  6,  7,  8,  0,  0,  0,  0,  0,  0,  0, &
        &    2,  3,  6,  7,  9, 10, 11, 12, 13,  0,  0,  0,  0,  0,  0, &
        &   15, 17, 22, 24,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0, &
        &    1,  4,  5,  8, 14, 16, 18, 19, 20, 21, 23, 25, 26, 27,  0, &
        & & ! GEOMETRY_DDN_TOP
        &    1,  2,  3,  4,  5,  6,  7,  8,  0,  0,  0,  0,  0,  0,  0, &
        &    2,  3,  6,  7,  9, 10, 11, 12, 13,  0,  0,  0,  0,  0,  0, &
        &   15, 17, 22, 24,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0, &
        &    1,  4,  5,  8, 14, 16, 18, 19, 20, 21, 23, 25, 26, 27,  0, &
        & & ! GEOMETRY_ANNULUS
        &    1,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0, &
        &    1,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0, &
        &    2,  3,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0, &
        &    0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0, &
        & & ! GEOMETRY_STELLARATORISLAND
        &    1,  2,  3,  4,  5,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0, &
        &    2,  3,  5,  6,  7,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0, &
        &    9, 11, 12, 13, 14, 15,  0,  0,  0,  0,  0,  0,  0,  0,  0, &
        &    1,  4,  8, 10,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0, &
        & & ! GEOMETRY_STRUCTURED_SPACES
        &    1,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0, &
        &    1,  2,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0, &
        &    3,  4,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0, &
        &    0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0, &
        & & ! GEOMETRY_LFS_SNOWFLAKE_MINUS
        &    1,  2,  3,  4,  5,  6,  7,  0,  0,  0,  0,  0,  0,  0,  0, &
        &    2,  3,  6,  7,  9, 10, 11, 12, 13,  0,  0,  0,  0,  0,  0, &
        &   15, 17,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0, &
        &    1,  4,  5,  8, 14, 16, 18, 19, 20, 21, 22, 23, 24, 25, 26, &
        & & ! GEOMETRY_LFS_SNOWFLAKE_PLUS
        &    1,  2,  3,  4,  5,  6,  7,  0,  0,  0,  0,  0,  0,  0,  0, &
        &    2,  3,  6,  7,  9, 10, 11, 12, 13,  0,  0,  0,  0,  0,  0, &
        &   15, 17,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0, &
        &    1,  4,  5,  8, 14, 16, 18, 19, 20, 21, 22, 23, 24, 25, 26, &
        & & ! GEOMETRY_GENERAL
        &    0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0, &
        &    0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0, &
        &    0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0, &
        &    0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0  &
        & &
        &  /), &
        & (/REGION_COUNT_MAX, REGIONTYPE_COUNT, GEOMETRY_COUNT/) )

   integer, dimension(REGION_NUMBER_MAX, 0:GEOMETRY_COUNT-1) :: boundaryAssignments = &
        &   reshape( (/                                                       &
        & & ! GEOMETRY_UNSPECIFIED
        &    NODIRECTION, NODIRECTION, NODIRECTION, NODIRECTION, NODIRECTION, &
        &    NODIRECTION, NODIRECTION, NODIRECTION, NODIRECTION, NODIRECTION, &
        &    NODIRECTION, NODIRECTION, NODIRECTION, NODIRECTION, NODIRECTION, &
        &    NODIRECTION, NODIRECTION, NODIRECTION, NODIRECTION, NODIRECTION, &
        &    NODIRECTION, NODIRECTION, NODIRECTION, NODIRECTION, NODIRECTION, &
        &    NODIRECTION, NODIRECTION,                                        &
        & & ! GEOMETRY_LINEAR
        &    WEST, EAST, SOUTH, NORTH,                           NODIRECTION, &
        &    NODIRECTION, NODIRECTION, NODIRECTION, NODIRECTION, NODIRECTION, &
        &    NODIRECTION, NODIRECTION, NODIRECTION, NODIRECTION, NODIRECTION, &
        &    NODIRECTION, NODIRECTION, NODIRECTION, NODIRECTION, NODIRECTION, &
        &    NODIRECTION, NODIRECTION, NODIRECTION, NODIRECTION, NODIRECTION, &
        &    NODIRECTION, NODIRECTION,                                        &
        & & ! GEOMETRY_CYLINDER
        &    NODIRECTION, SOUTH, NORTH,             NODIRECTION, NODIRECTION, &
        &    NODIRECTION, NODIRECTION, NODIRECTION, NODIRECTION, NODIRECTION, &
        &    NODIRECTION, NODIRECTION, NODIRECTION, NODIRECTION, NODIRECTION, &
        &    NODIRECTION, NODIRECTION, NODIRECTION, NODIRECTION, NODIRECTION, &
        &    NODIRECTION, NODIRECTION, NODIRECTION, NODIRECTION, NODIRECTION, &
        &    NODIRECTION, NODIRECTION,                                        &
        & & ! GEOMETRY_LIMITER
        &    WEST, EAST, NODIRECTION, SOUTH, NODIRECTION,                     &
        &    NORTH,       NODIRECTION, NODIRECTION, NODIRECTION, NODIRECTION, &
        &    NODIRECTION, NODIRECTION, NODIRECTION, NODIRECTION, NODIRECTION, &
        &    NODIRECTION, NODIRECTION, NODIRECTION, NODIRECTION, NODIRECTION, &
        &    NODIRECTION, NODIRECTION, NODIRECTION, NODIRECTION, NODIRECTION, &
        &    NODIRECTION, NODIRECTION,                                        &
        & & ! GEOMETRY_SN
        &    WEST, NODIRECTION, NODIRECTION, EAST, NODIRECTION,               &
        &    NODIRECTION, SOUTH, SOUTH, SOUTH, NODIRECTION,                   &
        &    NORTH, NORTH, NORTH,                   NODIRECTION, NODIRECTION, &
        &    NODIRECTION, NODIRECTION, NODIRECTION, NODIRECTION, NODIRECTION, &
        &    NODIRECTION, NODIRECTION, NODIRECTION, NODIRECTION, NODIRECTION, &
        &    NODIRECTION, NODIRECTION,                                        &
        & & ! GEOMETRY_CDN
        &    WEST, NODIRECTION, NODIRECTION, EAST, WEST,                      &
        &    NODIRECTION, NODIRECTION, EAST, NODIRECTION, NODIRECTION,        &
        &    NODIRECTION, NODIRECTION, SOUTH, SOUTH, SOUTH,                   &
        &    NODIRECTION, NORTH, NORTH, NORTH, SOUTH,                         &
        &    SOUTH, SOUTH, NODIRECTION, NORTH, NORTH,                         &
        &    NORTH,                                              NODIRECTION, &
        & & ! GEOMETRY_DDN_BOTTOM
        &    WEST, NODIRECTION, NODIRECTION, EAST, WEST,                      &
        &    NODIRECTION, NODIRECTION, EAST, NODIRECTION, NODIRECTION,        &
        &    NODIRECTION, NODIRECTION, NODIRECTION, SOUTH, SOUTH,             &
        &    SOUTH, NODIRECTION, NORTH, NORTH, NORTH,                         &
        &    SOUTH, SOUTH, SOUTH, NODIRECTION, NORTH,                         &
        &    NORTH, NORTH,                                                    &
        & & ! GEOMETRY_DDN_TOP
        &    WEST, NODIRECTION, NODIRECTION, EAST, WEST,                      &
        &    NODIRECTION, NODIRECTION, EAST, NODIRECTION, NODIRECTION,        &
        &    NODIRECTION, NODIRECTION, NODIRECTION, SOUTH, SOUTH,             &
        &    SOUTH, NODIRECTION, NORTH, NORTH, NORTH,                         &
        &    SOUTH, SOUTH, SOUTH, NODIRECTION, NORTH,                         &
        &    NORTH, NORTH,                                                    &
        & & ! GEOMETRY_ANNULUS
        &    NODIRECTION, SOUTH, NORTH,             NODIRECTION, NODIRECTION, &
        &    NODIRECTION, NODIRECTION, NODIRECTION, NODIRECTION, NODIRECTION, &
        &    NODIRECTION, NODIRECTION, NODIRECTION, NODIRECTION, NODIRECTION, &
        &    NODIRECTION, NODIRECTION, NODIRECTION, NODIRECTION, NODIRECTION, &
        &    NODIRECTION, NODIRECTION, NODIRECTION, NODIRECTION, NODIRECTION, &
        &    NODIRECTION, NODIRECTION,                                        &
        & & ! GEOMETRY_STELLARATORISLAND
        &    WEST, NODIRECTION, NODIRECTION, EAST, NODIRECTION,               &
        &    NODIRECTION, NODIRECTION, SOUTH, SOUTH, SOUTH,                   &
        &    NODIRECTION, NODIRECTION, NORTH,       NODIRECTION, NODIRECTION, &
        &    NODIRECTION, NODIRECTION, NODIRECTION, NODIRECTION, NODIRECTION, &
        &    NODIRECTION, NODIRECTION, NODIRECTION, NODIRECTION, NODIRECTION, &
        &    NODIRECTION, NODIRECTION,                                        &
        & & ! GEOMETRY_STRUCTURED_SPACES
        &    WEST, EAST, SOUTH, NORTH,                           NODIRECTION, &
        &    NODIRECTION, NODIRECTION, NODIRECTION, NODIRECTION, NODIRECTION, &
        &    NODIRECTION, NODIRECTION, NODIRECTION, NODIRECTION, NODIRECTION, &
        &    NODIRECTION, NODIRECTION, NODIRECTION, NODIRECTION, NODIRECTION, &
        &    NODIRECTION, NODIRECTION, NODIRECTION, NODIRECTION, NODIRECTION, &
        &    NODIRECTION, NODIRECTION,                                        &
        & & ! GEOMETRY_LFS_SNOWFLAKE_MINUS
        &    WEST, NODIRECTION, NODIRECTION, EAST, WEST,                      &
        &    NODIRECTION, NODIRECTION, WEST, NODIRECTION, NODIRECTION,        &
        &    NODIRECTION, NODIRECTION, NODIRECTION, SOUTH, SOUTH,             &
        &    SOUTH, NODIRECTION, NORTH, NORTH, NORTH,                         &
        &    SOUTH, SOUTH, SOUTH, NORTH, NORTH,                               &
        &    NORTH,                                              NODIRECTION, &
        & & ! GEOMETRY_LFS_SNOWFLAKE_PLUS
        &    WEST, NODIRECTION, NODIRECTION, EAST, WEST,                      &
        &    NODIRECTION, NODIRECTION, WEST, NODIRECTION, NODIRECTION,        &
        &    NODIRECTION, NODIRECTION, NODIRECTION, SOUTH, SOUTH,             &
        &    SOUTH, NODIRECTION, NORTH, NORTH, NORTH,                         &
        &    SOUTH, SOUTH, SOUTH, NORTH, NORTH,                               &
        &    NORTH,                                              NODIRECTION, &
        & & ! GEOMETRY_GENERAL
        &    NODIRECTION, NODIRECTION, NODIRECTION, NODIRECTION, NODIRECTION, &
        &    NODIRECTION, NODIRECTION, NODIRECTION, NODIRECTION, NODIRECTION, &
        &    NODIRECTION, NODIRECTION, NODIRECTION, NODIRECTION, NODIRECTION, &
        &    NODIRECTION, NODIRECTION, NODIRECTION, NODIRECTION, NODIRECTION, &
        &    NODIRECTION, NODIRECTION, NODIRECTION, NODIRECTION, NODIRECTION, &
        &    NODIRECTION, NODIRECTION                                         &
        & &
        &  /), &
        & (/REGION_NUMBER_MAX, GEOMETRY_COUNT/) )

contains

  !> Identify what geometry/topology is present from cut and periodicity data.
  !> Returns one of the GEOMETRY_* constants. Stops if unknown geometry.
  !> The object variable indicates for what object we are asking the topology of:
  !> 1: The GRID topology
  !> 2: The PLASMA topology
  !> The number of regions mpg%nnreg(0) and face indices mpg%fcReg refer to grid regions
  integer function geometryId( mpg, geo, object )
    use b2mod_types
    use b2us_map
    use b2us_geo
    use b2mod_ad &
     & , only: firstgmid
    implicit none
    type(mapping), intent(in) :: mpg
    type(geometry), intent(in) :: geo
    integer, intent(in) :: object
    integer :: i, j, iCv, nActive
    real(kind=R8) :: Xpsi_active, Xpsi_snowflake
    logical :: active, is_above_plasma, declaredOk
    logical, allocatable :: is_active(:)
    character(len=96) :: reason
    external xerrab, xertst

    call xertst ( object.eq.1.or.object.eq.2, 'incorrect object setting in geometryId')

    ! Explicit topology metadata is authoritative. Trust a declared
    ! standard family only after checking it against the serialized
    ! records. Reject malformed declarations rather than mixing the
    ! declared basic family with GEOMETRY_GENERAL behavior.
    if (mpg%hasTopologicalData) then
        geometryId = GEOMETRY_GENERAL
        if (mpg%geometryID.ge.GEOMETRY_LINEAR .and. &
          & mpg%geometryID.le.GEOMETRY_DDN_TOP) then
            call validateDeclaredGeometry(mpg, geo, declaredOk, reason)
            if (declaredOk) then
                geometryId = mpg%geometryID
                if (firstgmid(object)) call logmsg(LOGDEBUG, &
                  & "b2mod_geometry.geometryId(): using declared " // &
                  & "GEOMETRY_ID")
            else
                write(*,'(a,i0)') &
                  & "b2mod_geometry.geometryId(): invalid " // &
                  & "GEOMETRY_ID ", mpg%geometryID
                write(*,'(a,a)') &
                  & "  topology metadata check failed: ", trim(reason)
                call xerrab("Declared GEOMETRY_ID is inconsistent " // &
                  & "with explicit topology metadata")
            end if
        else if (mpg%geometryID.eq.GEOMETRY_GENERAL) then
            if (firstgmid(object)) call logmsg(LOGDEBUG, &
              & "b2mod_geometry.geometryId(): declared " // &
              & "GEOMETRY_GENERAL")
        else
            write(*,'(a,i0)') &
              & "b2mod_geometry.geometryId(): unsupported " // &
              & "GEOMETRY_ID ", mpg%geometryID
            call xerrab("Unsupported declared GEOMETRY_ID in " // &
              & "explicit topology metadata")
        end if
        firstgmid(object) = .false.
        return
    end if

    if (mpg%nnreg(0) == 1 .and. mpg%periodic_bc.le.0) then
        geometryId = GEOMETRY_LINEAR
        if (firstgmid(1)) then
            call logmsg( LOGDEBUG, "b2mod_connectivity.geometryId(): identified GEOMETRY_LINEAR")
            firstgmid(1) = .false.
        end if
        return
    end if

    if (mpg%nnreg(0) == 1 .and. mpg%periodic_bc == 1) then
        geometryId = GEOMETRY_CYLINDER
        if (firstgmid(1)) then
            call logmsg( LOGDEBUG, "b2mod_connectivity.geometryId(): identified GEOMETRY_CYLINDER")
            firstgmid(1) = .false.
        end if
        return
    end if

    if (mpg%nnreg(0) == 2) then
        geometryId = GEOMETRY_LIMITER
        if (firstgmid(1)) then
            call logmsg( LOGDEBUG, "b2mod_connectivity.geometryId(): identified GEOMETRY_LIMITER")
            firstgmid(1) = .false.
        end if
        return
    end if

    if (mpg%nnreg(0) == 4) then
        if (object.eq.1) then
          geometryId = GEOMETRY_SN
          if (firstgmid(1)) then
            call logmsg( LOGDEBUG, "b2mod_connectivity.geometryId(): identified grid GEOMETRY_SN")
            firstgmid(1) = .false.
          end if
        else if (object.eq.2) then
          if (maxval(mpg%fcReg(1:mpg%nFc)).le.6) then
            geometryId = GEOMETRY_LIMITER
            if (firstgmid(2)) then
              call logmsg( LOGDEBUG, "b2mod_connectivity.geometryId(): identified plasma GEOMETRY_LIMITER")
              firstgmid(2) = .false.
            end if
          else
            geometryId = GEOMETRY_SN
            if (firstgmid(2)) then
              call logmsg( LOGDEBUG, "b2mod_connectivity.geometryId(): identified plasma GEOMETRY_SN")
              firstgmid(2) = .false.
            end if
          end if
        end if
        return
    end if

    if (mpg%nnreg(0) == 5) then
        geometryId = GEOMETRY_STELLARATORISLAND
        if (firstgmid(1)) then
            call logmsg( LOGDEBUG, "b2mod_connectivity.geometryId(): identified GEOMETRY_STELLARATORISLAND")
            firstgmid(1) = .false.
        end if
        return
    end if

    if (mpg%nnreg(0) == 7) then
        Xpsi_active = geo%fsPsi(mpg%vxFs(mpg%Xpt(mpg%iactive)))
        Xpsi_snowflake = geo%fsPsi(mpg%vxFs(mpg%Xpt(mpg%iinactive)))
        if ((Xpsi_active.lt.Xpsi_snowflake.and.geo%psi_increasing).or. &
          & (Xpsi_active.gt.Xpsi_snowflake.and..not.geo%psi_increasing)) then
            geometryId = GEOMETRY_LFS_SNOWFLAKE_MINUS
            if (firstgmid(1)) then
                call logmsg( LOGDEBUG, "b2mod_connectivity.geometryId(): identified GEOMETRY_LFS_SNOWFLAKE_MINUS")
                firstgmid(1) = .false.
            end if
            return
        else if ((Xpsi_active.gt.Xpsi_snowflake.and.geo%psi_increasing).or. &
              &  (Xpsi_active.lt.Xpsi_snowflake.and..not.geo%psi_increasing)) then
            geometryId = GEOMETRY_LFS_SNOWFLAKE_PLUS
            if (firstgmid(1)) then
                call logmsg( LOGDEBUG, "b2mod_connectivity.geometryId(): identified GEOMETRY_LFS_SNOWFLAKE_PLUS")
                firstgmid(1) = .false.
            end if
            return
        end if
        if (mpg%vxFs(mpg%Xpt(mpg%iactive)) == mpg%vxFs(mpg%Xpt(mpg%iinactive))) then
            geometryId = GEOMETRY_UNSPECIFIED
            if (firstgmid(1)) then
                call logmsg( LOGDEBUG, "b2mod_connectivity.geometryId(): unknown GEOMETRY_UNSPECIFIED")
                firstgmid(1) = .false.
            end if
            return
        end if
    end if

    if (mpg%nnreg(0) == 8) then

      if (object.eq.1) then
        if (mpg%nXpt.le.1) then
          geometryId = GEOMETRY_DDN_BOTTOM ! Arbitrary ambiguous GEOMETRY_DDN assignment
          if (firstgmid(1)) then
            call logmsg( LOGDEBUG, "b2mod_connectivity.geometryId(): identified grid GEOMETRY_DDN(_BOTTOM?)")
            firstgmid(1) = .false.
          end if
          return
        else if (mpg%nXpt.eq.2) then ! We have a true DN grid
          if (mpg%vxFs(mpg%Xpt(1)) == mpg%vxFs(mpg%Xpt(2))) then
            geometryId = GEOMETRY_CDN
            if (firstgmid(1)) then
              call logmsg( LOGDEBUG, "b2mod_connectivity.geometryId(): identified grid GEOMETRY_CDN")
              firstgmid(1) = .false.
            end if
            return
          else
            active = .false.
            do i = mpg%vxCvP(mpg%Xpt(1),1), mpg%vxCvP(mpg%Xpt(1),1) + &
                                          & mpg%vxCvP(mpg%Xpt(1),2) - 1
              iCv = mpg%vxCv(i)
              if (mpg%cvReg(iCv).eq.1) active = .true.
            end do
            if ((geo%vxY(mpg%Xpt(1)) < geo%vxY(mpg%Xpt(2)).and.active).or. &
           &    (geo%vxY(mpg%Xpt(1)) > geo%vxY(mpg%Xpt(2)).and..not.active)) then
              geometryId = GEOMETRY_DDN_BOTTOM
              if (firstgmid(1)) then
                call logmsg( LOGDEBUG, "b2mod_connectivity.geometryId(): identified grid GEOMETRY_DDN_BOTTOM")
                firstgmid(1) = .false.
              end if
              return
            else
              geometryId = GEOMETRY_DDN_TOP
              if (firstgmid(1)) then
                call logmsg( LOGDEBUG, "b2mod_connectivity.geometryId(): identified grid GEOMETRY_DDN_TOP")
                firstgmid(1) = .false.
              end if
              return
            end if
          end if
        else ! We have more than 2 X-points, to avoid a crash,
             ! we base the plasma topology only on the active X-point
             ! and its location
          allocate(is_active(mpg%nXpt))
          is_active = .false.
          is_above_plasma = .false.
          do i = 1, mpg%nXpt
            do j = mpg%vxCvP(mpg%Xpt(i),1), mpg%vxCvP(mpg%Xpt(i),1) + &
                                          & mpg%vxCvP(mpg%Xpt(i),2) - 1
              iCv = mpg%vxCv(j)
              if (mpg%cvReg(iCv).eq.1) then
                is_active(i) = .true.
                is_above_plasma = geo%vxY(mpg%Xpt(i)).gt.geo%cvY(iCv)
              end if
            end do
          end do
          nActive = count( is_active(1:mpg%nXpt) .eqv. .true. )
          if (nActive.eq.2) then
            geometryId = GEOMETRY_CDN
            if (firstgmid(1)) then
              call logmsg( LOGDEBUG, "b2mod_connectivity.geometryId(): identified grid GEOMETRY_CDN")
              call logmsg( LOGWARNING, "This grid contains more than two X-points!")
              call logmsg( LOGWARNING, "Standard diagnostic output may fail.")
              firstgmid(1) = .false.
            end if
            deallocate( is_active )
            return
          else if (nActive.eq.1 .and. is_above_plasma) then
            geometryId = GEOMETRY_DDN_TOP
            if (firstgmid(1)) then
              call logmsg( LOGDEBUG, "b2mod_connectivity.geometryId(): identified grid GEOMETRY_DDN_TOP")
              call logmsg( LOGWARNING, "This grid contains more than two X-points!")
              call logmsg( LOGWARNING, "Standard diagnostic output may fail.")
              firstgmid(1) = .false.
            end if
            deallocate( is_active )
            return
          else if (nActive.eq.1 .and. .not.is_above_plasma) then
            geometryId = GEOMETRY_DDN_BOTTOM
            if (firstgmid(1)) then
              call logmsg( LOGDEBUG, "b2mod_connectivity.geometryId(): identified grid GEOMETRY_DDN_BOTTOM")
              call logmsg( LOGWARNING, "This grid contains more than two X-points!")
              call logmsg( LOGWARNING, "Standard diagnostic output may fail.")
              firstgmid(1) = .false.
            end if
            deallocate( is_active )
            return
          end if
        endif
      else if (object.eq.2) then
        if (maxval(mpg%fcReg(1:mpg%nFc)).le.6) then
          geometryId = GEOMETRY_LIMITER
          if (firstgmid(2)) then
            call logmsg( LOGDEBUG, "b2mod_connectivity.geometryId(): identified plasma GEOMETRY_LIMITER")
            firstgmid(2) = .false.
          end if
          return
        elseif (mpg%nXpt.eq.1) then !nh only 1 X-point for vessel mode grids
          geometryID = GEOMETRY_SN
          if (firstgmid(2)) then
            call logmsg( LOGDEBUG, "b2mod_connectivity.geometryId(): identified plasma GEOMETRY_SN")
            firstgmid(2) = .false.
          end if
          return
        else if (mpg%nXpt.eq.2) then ! The normal DN case
          if (mpg%vxFs(mpg%Xpt(1)) == mpg%vxFs(mpg%Xpt(2))) then
            geometryId = GEOMETRY_CDN
            if (firstgmid(2)) then
              call logmsg( LOGDEBUG, "b2mod_connectivity.geometryId(): identified plasma GEOMETRY_CDN")
              firstgmid(2) = .false.
            end if
            return
          else
            active = .false.
            do i = mpg%vxCvP(mpg%Xpt(1),1), mpg%vxCvP(mpg%Xpt(1),1) + &
                                          & mpg%vxCvP(mpg%Xpt(1),2) - 1
              iCv = mpg%vxCv(i)
              if (mpg%cvReg(iCv).eq.1) active = .true.
            end do
            if ((geo%vxY(mpg%Xpt(1)) < geo%vxY(mpg%Xpt(2)).and.active).or. &
             &  (geo%vxY(mpg%Xpt(1)) > geo%vxY(mpg%Xpt(2)).and..not.active)) then
              geometryId = GEOMETRY_DDN_BOTTOM
              if (firstgmid(2)) then
                call logmsg( LOGDEBUG, "b2mod_connectivity.geometryId(): identified plasma GEOMETRY_DDN_BOTTOM")
                firstgmid(2) = .false.
              end if
              return
            else
              geometryId = GEOMETRY_DDN_TOP
              if (firstgmid(2)) then
                call logmsg( LOGDEBUG, "b2mod_connectivity.geometryId(): identified plasma GEOMETRY_DDN_TOP")
                firstgmid(2) = .false.
              end if
              return
            end if
          end if
        else ! We have more than 2 X-points,
             ! we base the plasma topology only on the active X-point
             ! and its location
          allocate(is_active(mpg%nXpt))
          is_active = .false.
          is_above_plasma = .false.
          do i = 1, mpg%nXpt
            do j = mpg%vxCvP(mpg%Xpt(i),1), mpg%vxCvP(mpg%Xpt(i),1) + &
                                          & mpg%vxCvP(mpg%Xpt(i),2) - 1
              iCv = mpg%vxCv(j)
              if (mpg%cvReg(iCv).eq.1) then
                is_active(i) = .true.
                is_above_plasma = geo%vxY(mpg%Xpt(i)).gt.geo%cvY(iCv)
              end if
            end do
          end do
          nActive = count( is_active(1:mpg%nXpt) .eqv. .true. )
          if (nActive.eq.2) then
            geometryId = GEOMETRY_CDN
            if (firstgmid(2)) then
              call logmsg( LOGDEBUG, "b2mod_connectivity.geometryId(): identified plasma GEOMETRY_CDN")
              firstgmid(2) = .false.
            end if
            deallocate( is_active )
            return
          else if (nActive.eq.1 .and. is_above_plasma) then
            geometryId = GEOMETRY_DDN_TOP
            if (firstgmid(2)) then
              call logmsg( LOGDEBUG, "b2mod_connectivity.geometryId(): identified plasma GEOMETRY_DDN_TOP")
              firstgmid(2) = .false.
            end if
            deallocate( is_active )
            return
          else if (nActive.eq.1 .and. .not.is_above_plasma) then
            geometryId = GEOMETRY_DDN_BOTTOM
            if (firstgmid(2)) then
              call logmsg( LOGDEBUG, "b2mod_connectivity.geometryId(): identified plasma GEOMETRY_DDN_BOTTOM")
              firstgmid(2) = .false.
            end if
            deallocate( is_active )
            return
          end if
        end if
      end if
    end if

    geometryId = GEOMETRY_UNSPECIFIED
    if (firstgmid(1)) then
        call logmsg( LOGDEBUG, "b2mod_connectivity.geometryId(): unknown GEOMETRY_UNSPECIFIED")
        firstgmid(1) = .false.
    end if

    call xerrab ( 'b2mod_connectivity.geometryId: unknown geometry' )

  end function geometryId


  !> Validate a serialized basic-family GEOMETRY_ID against the
  !> topology records supplied by the grid producer.
  subroutine validateDeclaredGeometry( mpg, geo, ok, reason )
    use b2mod_types &
     & , only: R8
    use b2us_map &
     & , only: mapping
    use b2us_geo &
     & , only: geometry
    implicit none
    type(mapping), intent(in) :: mpg
    type(geometry), intent(in) :: geo
    logical, intent(out) :: ok
    character(len=*), intent(out) :: reason
    integer :: nPrimExp, nRegBase, iCv, iXpt, nCore
    integer :: nPrimBelow, nPrimAbove, nSecBelow, nSecAbove
    real(kind=R8) :: yRef

    ! GOAT primary means that the X-point is connected directly to the
    ! core: CDN has two primaries, while DDN has one primary and one
    ! secondary X-point.
    select case (mpg%geometryID)
    case (GEOMETRY_LINEAR, GEOMETRY_CYLINDER)
        nPrimExp = 0
        nRegBase = 1
    case (GEOMETRY_LIMITER)
        nPrimExp = 0
        nRegBase = 2
    case (GEOMETRY_SN)
        nPrimExp = 1
        nRegBase = 4
    case (GEOMETRY_CDN)
        nPrimExp = 2
        nRegBase = 8
    case default
        nPrimExp = 1
        nRegBase = 8
    end select

    ok = .false.
    reason = ' '
    if (mpg%nOpt.ge.2) then
        reason = 'more than one O-point'
        return
    end if
    if (maxval(mpg%cvReg).lt.nRegBase) then
        reason = 'too few volume regions'
        return
    end if
    if (mpg%geometryID.eq.GEOMETRY_LINEAR .and. &
      & mpg%periodic_bc.gt.0) then
        reason = 'linear grid is periodic'
        return
    end if
    if (mpg%geometryID.eq.GEOMETRY_CYLINDER .and. &
      & mpg%periodic_bc.ne.1) then
        reason = 'cylindrical grid is not periodic'
        return
    end if
    if (.not.allocated(mpg%isPrimaryXpt)) then
        ok = .true.
        return
    end if
    if (count(mpg%isPrimaryXpt.eq.1).ne.nPrimExp) then
        reason = 'primary X-point count'
        return
    end if
    if (nPrimExp.eq.0) then
        ok = .true.
        return
    end if

    ! Match the GOAT side test: use the O-point height when available,
    ! otherwise the centroid height of the core cells.
    if (mpg%nOpt.eq.1 .and. allocated(mpg%Opt)) then
        yRef = geo%vxY(mpg%Opt(1))
    else
        yRef = 0.0_R8
        nCore = 0
        do iCv = 1, mpg%nCi
            if (mpg%cvReg(iCv).ge.1 .and. &
              & mpg%cvReg(iCv).le.nRegBase .and. &
              & mod(mpg%cvReg(iCv)-1,4).eq.0) then
                yRef = yRef + geo%cvY(iCv)
                nCore = nCore + 1
            end if
        end do
        if (nCore.eq.0) then
            reason = 'no core cells'
            return
        end if
        yRef = yRef/real(nCore, R8)
    end if

    nPrimBelow = 0
    nPrimAbove = 0
    nSecBelow = 0
    nSecAbove = 0
    do iXpt = 1, mpg%nXpt
        if (geo%vxY(mpg%Xpt(iXpt)).ge.yRef) then
            if (mpg%isPrimaryXpt(iXpt).eq.1) then
                nPrimAbove = nPrimAbove + 1
            else
                nSecAbove = nSecAbove + 1
            end if
        else
            if (mpg%isPrimaryXpt(iXpt).eq.1) then
                nPrimBelow = nPrimBelow + 1
            else
                nSecBelow = nSecBelow + 1
            end if
        end if
    end do

    select case (mpg%geometryID)
    case (GEOMETRY_SN)
        ok = (nPrimBelow.eq.1 .and. &
           & nPrimAbove+nSecAbove.eq.0) .or. &
           & (nPrimAbove.eq.1 .and. nPrimBelow+nSecBelow.eq.0)
    case (GEOMETRY_CDN)
        ok = nPrimBelow.eq.1 .and. nPrimAbove.eq.1
    case (GEOMETRY_DDN_BOTTOM)
        if (mpg%isClassicalGrid.eq.1) then
            ok = nPrimBelow.eq.1
        else
            ok = nPrimBelow.eq.1 .and. nSecAbove.ge.1
        end if
    case (GEOMETRY_DDN_TOP)
        if (mpg%isClassicalGrid.eq.1) then
            ok = nPrimAbove.eq.1
        else
            ok = nPrimAbove.eq.1 .and. nSecBelow.ge.1
        end if
    end select
    if (.not.ok) reason = 'X-point positions relative to the core'

  end subroutine validateDeclaredGeometry


  !> Return number of regions of a given type for a given grid geometry.
  !> @param geometryId: see definition of GEOMETRY_* constants.
  !> @param regionType: see definition of REGIONTYPE_* constants.
  integer function regionCount( geometryId, regionType )
    integer, intent(in) :: geometryId, regionType

    !! TODO: maybe do some input parameter checking
    regionCount = regionCounts(regionType, geometryId)
  end function regionCount


  !> Return total number of regions (both cell and face regions) for
  !> the given grid geometry
  integer function regionCountTotal( geometryId )
    integer, intent(in) :: geometryId

    !! internal
    integer :: iType

    regionCountTotal = 0
    do iType = 0, REGIONTYPE_COUNT - 1
        regionCountTotal = regionCountTotal + regionCount( geometryId, iType )
    end do
  end function regionCountTotal


  !> Return name of a specific grid region.
  !> @param geometryId: see definition of GEOMETRY_* constants.
  !> @param regionType: see definition of REGIONTYPE_* constants.
  !> @param regionId: number of region. Should be between [1, regionCount(geometryId, regionType)].
  character(32) function regionName( geometryId, regionType, regionId )
    integer, intent(in) :: geometryId, regionType, regionId

    regionName = regionNames(regionId, regionType, geometryId)
  end function regionName

  !> Return number of a specific grid region as per the fcReg numbering
  !> @param geometryId: see definition of GEOMETRY_* constants.
  !> @param regionType: see definition of REGIONTYPE_* constants.
  !> @param regionId: number of region. Should be between [1, regionCount(geometryId, regionType)].
  integer function regionNumber( geometryId, regionType, regionId )
    integer, intent(in) :: geometryId, regionType, regionId

    regionNumber = regionNumbers(regionId, regionType, geometryId)
  end function regionNumber

end module b2mod_geometry

!!!Local Variables:
!!! mode: f90
!!! End:
