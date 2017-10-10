module b2mod_ual_io_grid

  !> This module provides:
  !>
  !> 2. A routine (b2ITMFillGridDescription/b2IMASFillGridDescription) to write 
  !> the B2 grid into an ITM/IMAS IDS grid description data structure (which 
  !> usually is part of a CPO/IDS). It also sets up the default subgrids for 
  !> the B2 grid.
  !> 
  !> 3. Routines to transform variables stored in the B2 data structure into
  !> the form expected CPO data structure.

    use b2mod_types , B2_R8 => R8, B2_R4 => R4
#ifdef IMAS
    use ids_schemas  ! IGNORE
    use ids_routines ! IGNORE
#else
# ifdef ITM
    use itm_types , ITM_R8 => R8, ITM_R4 => R4 ! IGNORE
    Use euITM_schemas ! IGNORE
    use itm_constants , pi => itm_pi ! IGNORE
# endif
#endif
    use ggd_common
    use helper
    use logging , only: logmsg, LOGDEBUG
    use ggd_assert
    use string
    use ggd , GGD_UNDEFINED => GRID_UNDEFINED
    use ggd_structured
    use b2mod_connectivity , REMOVED_B2_R8 => R8
    use carre_constants
    use b2mod_cellhelper

    use b2mod_grid_mapping
    use b2mod_indirect

implicit none

    !> Constants for use with the ITM grid description

    !> Space indices 
    integer, parameter :: SPACE_POLOIDALPLANE = 1
    integer, parameter :: SPACE_TOROIDALANGLE = 2

    !> Space setup: 
    !  SPACE_COUNT = SPACE_POLOIDALPLANE: do only do the poloidal plane space,
    !  SPACE_COUNT = SPACE_TOROIDALANGLE: will do the full 3d grid with two spaces.
    integer, parameter :: SPACE_COUNT = SPACE_POLOIDALPLANE

    !> Will have at maximum so many spaces
    integer, parameter :: SPACE_COUNT_MAX = 2

    !> Number of points in toroidal direction (only 1 makes sense here, this is for playing around)
    integer, parameter :: NNODES_TOROIDAL = 1

    !> Flag controlling whether the toroidal space is set up periodic or not.
    !> If periodic, the last node is connect to the first by an edge.
    !> If not periodic, an additional node at 2 pi is added
    logical, parameter :: TOROIDAL_PERIODIC = .false.


    !> Object class tuples
    integer, dimension(SPACE_COUNT_MAX), parameter :: CLASS_NODE = (/ 0, 0 /)

    integer, dimension(SPACE_COUNT_MAX), parameter :: CLASS_RZ_EDGE = (/ 1, 0 /)
    integer, dimension(SPACE_COUNT_MAX), parameter :: CLASS_PHI_EDGE = (/ 0, 1 /)

    integer, dimension(SPACE_COUNT_MAX), parameter :: CLASS_POLOIDALRADIAL_FACE = (/ 1, 1 /)
    integer, dimension(SPACE_COUNT_MAX), parameter :: CLASS_TOROIDAL_FACE = (/ 2, 0 /)

    integer, dimension(SPACE_COUNT_MAX), parameter :: CLASS_CELL = (/ 2, 1 /)


    !> Subgrid/Grid subset name constants

    !> Number of generic subgrids
    integer, parameter :: B2_GENERIC_SUBGRID_COUNT = 6

    ! Generic subgrids (all cells, all faces)
    ! Note: special subgrids (given by region ids) do not have specific constants (see also b2mod_connectivity.f90)

    !> Subgrid of all 2d cells
    integer, parameter :: B2_SUBGRID_CELLS = 1 
    !> Subgrid of all 0d nodes
    integer, parameter :: B2_SUBGRID_NODES = 2
    !> Subgrid of all faces in order given by grid map. First x-aligned faces, then y-aligned faces
    integer, parameter :: B2_SUBGRID_FACES = 3
    !> Subgrid of all x-aligned faces in order given by grid map. 
    integer, parameter :: B2_SUBGRID_FACES_X = 4
    !> Subgrid of all y-aligned faces in order given by grid map. 
    integer, parameter :: B2_SUBGRID_FACES_Y = 5
    !> Subgrid of all X-points
    integer, parameter :: B2_SUBGRID_XPOINTS = 6

!dpc  private :: R8

#ifdef IMAS

contains

    !> Routine that fills in a grid description which is part of a CPO
    !> using the given grid data and prepared mappings
    subroutine b2IMASFillGridDescription( gmap,ggd_grid, &
        & nx,ny,crx,cry, &
        & leftix,leftiy,rightix,rightiy, &
        & topix,topiy,bottomix,bottomiy,&
        & nnreg,topcut,region,cflag,includeGhostCells,vol,gs,qc )

        type(B2GridMap), intent(in) :: gmap
        type(ids_generic_grid_dynamic), intent(out) :: ggd_grid

        !> Size of grid arrays: (-1:nx, -1:ny) 
        integer, intent(in) :: nx, ny
        !>   .. output arguments
        !> vertex coordinates
        real (R8), intent(in) :: &
            & crx(-1:nx,-1:ny,0:3), cry(-1:nx,-1:ny,0:3)
        !> B2 connectivity array
        integer, intent(in) :: &
            & leftix(-1:nx,-1:ny),leftiy(-1:nx,-1:ny),&
            & rightix(-1:nx,-1:ny),rightiy(-1:nx,-1:ny),&
            & topix(-1:nx,-1:ny),topiy(-1:nx,-1:ny),&
            & bottomix(-1:nx,-1:ny),bottomiy(-1:nx,-1:ny)
        !> B2 region & cut information
        integer, intent(in) :: &
            & nnreg(0:2), topcut(:), &
            & region(-1:nx,-1:ny,0:2)
        !> Cell flags 
        integer cflag(-1:nx,-1:ny, CARREOUT_NCELLFLAGS)
        logical, intent(in) :: includeGhostCells
        !> Optional B2 measure information
        real(R8), intent(in), optional :: vol(-1:nx,-1:ny,0:4),     &
            &   gs(-1:nx,-1:ny,0:2), qc(-1:nx,-1:ny)

        !> internal
        integer, parameter :: NDIM = 2

        call assert( present(gs) .EQV. present(qc) )

        call fillInGridDescription()
#if 0
        call fillInSubGridDescription()
#endif

contains 

    !> Part 1: fill in grid description
    subroutine fillInGridDescription()

        !> internal
        integer :: ivx, ifc, icv, ix, iy, nix, niy, i, dir

        allocate( ggd_grid%space(SPACE_COUNT) )

        !> Coordinate types
        !> (dimension of space = NDIM = size( coordtype ) 
      
        ! ggd_grid%space(SPACE_POLOIDALPLANE)%geometry_type%name = 'Poloidal' 

        allocate( ggd_grid%space( SPACE_POLOIDALPLANE )%coordinates_type(NDIM) )
        ggd_grid%space( SPACE_POLOIDALPLANE )%coordinates_type(1) = COORDTYPE_R
        ggd_grid%space( SPACE_POLOIDALPLANE )%coordinates_type(2) = COORDTYPE_Z

        !> Have two types of objects: 0d nodes, 1d edges, 2d cells
        allocate( ggd_grid%space(SPACE_POLOIDALPLANE)%  &
            &   objects_per_dimension(NDIM + 1) )

        !> Allocate the number of objects of each type
        allocate( ggd_grid%space( SPACE_POLOIDALPLANE )%    &
            &   objects_per_dimension(1)%object( (nx+1)*(ny+1)) )       !> nodes
        allocate( ggd_grid%space( SPACE_POLOIDALPLANE )%    &
            &   objects_per_dimension(2)%object( nx*(ny+1)+(nx+1)*ny) ) !> edges
        allocate( ggd_grid%space( SPACE_POLOIDALPLANE )%    &
            &   objects_per_dimension(3)%object( nx*ny ) )              !> cells
      
        ! write(*,*) "* gmap%nvx: ", gmap%nvx
        !> Fill in node information
        allocate( ggd_grid%space( SPACE_POLOIDALPLANE )%    &
            &   objects_per_dimension(1)%object( gmap%nvx ) )
        do ivx = 1, gmap % nvx
            allocate( ggd_grid%space( SPACE_POLOIDALPLANE )%    &
                &   objects_per_dimension(1)%object( ivx )%geometry(2))
            allocate( ggd_grid%space( SPACE_POLOIDALPLANE )%    &
                &   objects_per_dimension(1)%object( ivx )%nodes(1))
            ! write(*,*) "** ivx: ", ivx
            ! write(*,*) "** gmap%mapVxix( ivx ): ", gmap%mapVxix( ivx )
            ! write(*,*) "** gmap%mapVxiy( ivx ): ", gmap%mapVxiy( ivx )
            ! write(*,*) "** gmap%mapVxIVx( ivx ): ", gmap%mapVxIVx( ivx )
            ! write(*,*) "** crx( ... ) : ",  &
                ! &   crx( gmap%mapVxix( ivx ), gmap%mapVxiy( ivx ), gmap%mapVxIVx( ivx ))
            ! write(*,*) "** cry( ... ) : ",  &
                ! &   cry( gmap%mapVxix( ivx ), gmap%mapVxiy( ivx ), gmap%mapVxIVx( ivx ))
            ggd_grid%space( SPACE_POLOIDALPLANE )%objects_per_dimension(1)%     &
                &   object( ivx )%nodes(1) = ivx             
            ggd_grid%space( SPACE_POLOIDALPLANE )%objects_per_dimension(1)%     &
                &   object( ivx )%geometry(1) = crx(    gmap%mapVxix( ivx ),    &
                                                    &   gmap%mapVxiy( ivx ),    &
                                                    &   gmap%mapVxIVx( ivx ))       
            ggd_grid%space( SPACE_POLOIDALPLANE )%objects_per_dimension(1)%     &
                &   object( ivx )%geometry(2) = cry(    gmap%mapVxix( ivx ),    &
                                                    &   gmap%mapVxiy( ivx ),    &
                                                    &   gmap%mapVxIVx( ivx ))       
        end do

# if 0

        ivx = 0
        do iy = 0, ny-1
            do ix = 0, nx-1
                write(*,*) "* ivx: ", ivx
                ivx = ivx + 1 !> Lower left corners
                write(*,*) "* crx(ix, iy, 0): ", crx(ix, iy, 0)
                write(*,*) "* cry(ix, iy, 0): ", cry(ix, iy, 0)
                ggd_grid%space(SPACE_POLOIDALPLANE)%    &
                    &   objects_per_dimension(1)%object(ivx)%geometry(1) =  &
                    &   crx(ix,iy,0)
                ggd_grid%space(SPACE_POLOIDALPLANE)%    &
                    &   objects_per_dimension(1)%object(ivx)%geometry(2) =  &
                    &   cry(ix,iy,0)
                if (ix.eq.nx-1) then
                    ivx = ivx + 1  !> Lower right corners
                    write(*,*) "* crx(ix, iy, 1): ", crx(ix, iy, 1)
                    write(*,*) "* cry(ix, iy, 1): ", cry(ix, iy, 1)
                    ggd_grid%space(SPACE_POLOIDALPLANE)%    &
                        &   objects_per_dimension(1)%object(ivx)%geometry(1) =  &
                        &   crx(ix,iy,1)
                    ggd_grid%space(SPACE_POLOIDALPLANE)%    &
                        &   objects_per_dimension(1)%object(ivx)%geometry(2) =  &
                        &   cry(ix,iy,1)
                end if
                if (iy.eq.ny-1) then
                    ivx = ivx + 1  !> Upper left corners
                    write(*,*) "* crx(ix, iy, 2): ", crx(ix, iy, 2)
                    write(*,*) "* cry(ix, iy, 2): ", cry(ix, iy, 2)
                    ggd_grid%space(SPACE_POLOIDALPLANE)%    &
                        &   objects_per_dimension(1)%object(ivx)%geometry(1) =  &
                        &   crx(ix,iy,2)
                    ggd_grid%space(SPACE_POLOIDALPLANE)%    &
                        &   objects_per_dimension(1)%object(ivx)%geometry(2) =  &
                        &   cry(ix,iy,2)
                    if (ix.eq.nx-1) then
                        ivx = ivx + 1  !> Upper right corners
                        write(*,*) "* crx(ix, iy, 3): ", crx(ix, iy, 3)
                        write(*,*) "* cry(ix, iy, 3): ", cry(ix, iy, 3)
                        ggd_grid%space(SPACE_POLOIDALPLANE)%    &
                            &   objects_per_dimension(1)%object(ivx)%   &
                            &   geometry(1) = crx(ix,iy,3)
                        ggd_grid%space(SPACE_POLOIDALPLANE)%    &
                            &   objects_per_dimension(1)%object(ivx)%   &
                            &   geometry(2) = cry(ix,iy,3)
                    end if
                end if
            end do
        end do



      ! Fill in object definitions (i.e. what objects compose an object)

      ! 1d objects: faces
      ! ...have two boundaries
      allocate( ggd_grid%space(SPACE_POLOIDALPLANE) % objects(2) % boundary( gmap%nfcx + gmap%nfcy, 2) )
      ! ...have two neighbours, in positive and negative coordinate direction, one on each side
      ! (for x-aligned faces: along flux surface, for y-aligned faces: orthogonal to flux surface)
      allocate( ggd_grid%space(SPACE_POLOIDALPLANE) % objects(2) % neighbour(gmap%nfcx + gmap%nfcy, 2, 1) )      
      ! 1d object measure: face area
      if (present(gs)) allocate( ggd_grid%space(SPACE_POLOIDALPLANE) % objects(2) % measure( gmap % nfcx + gmap % nfcy, 1 ) )
      ! first set all boundary & connectivity information to undefined
      ggd_grid%space(SPACE_POLOIDALPLANE) % objects(2) % boundary = GRID_UNDEFINED
      ggd_grid%space(SPACE_POLOIDALPLANE) % objects(2) % neighbour = GRID_UNDEFINED
      ! x-aligned faces    
      do ifc = 1, gmap % nfcx    
          ! get position of this face in the b2 grid
          ix = gmap % mapFcix( ifc )
          iy = gmap % mapFciy( ifc )
          ! get index of start vertex 
          ! objdef dims: index of face, 1=start node, 1=one-dimensional object
          select case ( gmap % mapFcIFace( ifc ) )
          case( BOTTOM )
             ! start index: 1=start node
             ggd_grid%space(SPACE_POLOIDALPLANE) % objects(2) % boundary( ifc, 1 ) = gmap % mapVxI( ix, iy, VX_LOWERLEFT )
             if (gmap % mapVxI( ix, iy, VX_LOWERLEFT ) == GRID_UNDEFINED) then
                call logmsg(LOGWARNING, "b2IMASFillGD: BOTTOM face at "//int2str(ix)//","//int2str(iy)//" has no start node")
             end if
             ! end vertex: 2=end node
             ggd_grid%space(SPACE_POLOIDALPLANE) % objects(2) % boundary( ifc, 2 ) = gmap % mapVxI( ix, iy, VX_LOWERRIGHT )
             if (gmap % mapVxI( ix, iy, VX_LOWERRIGHT ) == GRID_UNDEFINED) then
                call logmsg(LOGWARNING, "b2IMASFillGD: BOTTOM face at "//int2str(ix)//","//int2str(iy)//" has no end node")
             end if
          case( TOP )
             ! start index: 1=start node
             ggd_grid%space(SPACE_POLOIDALPLANE) % objects(2) % boundary( ifc, 1 ) = gmap % mapVxI( ix, iy, VX_UPPERLEFT )
             if (gmap % mapVxI( ix, iy, VX_UPPERLEFT ) == GRID_UNDEFINED) then
                call logmsg(LOGWARNING, "b2IMASFillGD: TOP face at "//int2str(ix)//","//int2str(iy)//" has no start node")
             end if
             ! end vertex: 2=end node
             ggd_grid%space(SPACE_POLOIDALPLANE) % objects(2) % boundary( ifc, 2 ) = gmap % mapVxI( ix, iy, VX_UPPERRIGHT )
             if (gmap % mapVxI( ix, iy, VX_UPPERRIGHT ) == GRID_UNDEFINED) then
                call logmsg(LOGWARNING, "b2IMASFillGD: TOP face at "//int2str(ix)//","//int2str(iy)//" has no end node")
             end if
          end select

          ! Neighbour faces of this face
          ! Left neighbour: face continuing to the left of this face
          nix = leftix( ix, iy )
          niy = leftiy( ix, iy )
          if ( .not. isUnneededCell( nx, ny, cflag, includeGhostCells, nix, niy ) ) then
             ggd_grid%space(SPACE_POLOIDALPLANE) % objects(2) % neighbour( ifc, 1, 1 ) = &
                  & gmap % mapFcI( nix, niy, gmap % mapFcIFace( ifc ) )
          end if
          ! Right neighbour: face continuing to the right of this face 
          nix = rightix( ix, iy )
          niy = rightiy( ix, iy )
          if ( .not. isUnneededCell( nx, ny, cflag, includeGhostCells, nix, niy ) ) then
             ggd_grid%space(SPACE_POLOIDALPLANE) % objects(2) % neighbour( ifc, 2, 1 ) = &
                  & gmap % mapFcI( nix, niy, gmap % mapFcIFace( ifc ) )
          end if

          ! measure: area
          if (present(gs)) ggd_grid%space(SPACE_POLOIDALPLANE) % objects(2) % measure( ifc, 1 ) = gs(ix, iy, ALIGNX)
      end do

      ! y-aligned faces    
      do ifc = gmap % nfcx + 1, gmap % nfcx + gmap % nfcy    
          ! get position of this face in the b2 grid
          ix = gmap % mapFcix( ifc )
          iy = gmap % mapFciy( ifc )
!!$          if (gmap%mapCvI(ix, iy) == GRID_UNDEFINED) then 
!!$                  call logmsg(LOGWARNING, "b2IMASFillGD: writing out faces for unused cell "//int2str(ix)//","//int2str(iy))
!!$          end if

          select case ( gmap % mapFcIFace( ifc ) )
          case( LEFT )
              ! start index: 1=start node
              ggd_grid%space(SPACE_POLOIDALPLANE) % objects(2) % boundary( ifc, 1 ) = gmap % mapVxI( ix, iy, VX_LOWERLEFT )
              if (gmap % mapVxI( ix, iy, VX_LOWERLEFT ) == GRID_UNDEFINED) then
                  call logmsg(LOGWARNING, "b2IMASFillGD: LEFT face at "//int2str(ix)//","//int2str(iy)//" has no start node")
              end if
          ! end vertex: 2=end node
              ggd_grid%space(SPACE_POLOIDALPLANE) % objects(2) % boundary( ifc, 2 ) = gmap % mapVxI( ix, iy, VX_UPPERLEFT )
              if (gmap % mapVxI( ix, iy, VX_UPPERLEFT ) == GRID_UNDEFINED) then
                  call logmsg(LOGWARNING, "b2IMASFillGD: LEFT face at "//int2str(ix)//","//int2str(iy)//" has no end node")
          end if
              !if (ggd_grid%space(SPACE_POLOIDALPLANE) % objects(2) % boundary( ifc, 1 )


          case( RIGHT )
              ! start index: 1=start node
              ggd_grid%space(SPACE_POLOIDALPLANE) % objects(2) % boundary( ifc, 1 ) = gmap % mapVxI( ix, iy, VX_LOWERRIGHT )
              if (gmap % mapVxI( ix, iy, VX_LOWERRIGHT ) == GRID_UNDEFINED) then
                  call logmsg(LOGWARNING, "b2IMASFillGD: RIGHT face at "//int2str(ix)//","//int2str(iy)//" has no start node")
          end if
              ! end vertex: 2=end node
              ggd_grid%space(SPACE_POLOIDALPLANE) % objects(2) % boundary( ifc, 2 ) = gmap % mapVxI( ix, iy, VX_UPPERRIGHT )
              if (gmap % mapVxI( ix, iy, VX_UPPERRIGHT ) == GRID_UNDEFINED) then
                  call logmsg(LOGWARNING, "b2IMASFillGD: RIGHT face at "//int2str(ix)//","//int2str(iy)//" has no end node")
              end if
          end select


          ! Neighbour faces of this face
          ! Bottom neighbour: face continuing to the bottom of this face
          nix = bottomix( ix, iy )
          niy = bottomiy( ix, iy )
          if ( .not. isUnneededCell( nx, ny, cflag, includeGhostCells, nix, niy ) ) then
             ggd_grid%space(SPACE_POLOIDALPLANE) % objects(2) % neighbour( ifc, 1, 1 ) = &
                  & gmap % mapFcI( nix, niy, gmap % mapFcIFace( ifc ) )
          end if
          ! Top neighbour: face continuing to the top of this face 
          nix = topix( ix, iy )
          niy = topiy( ix, iy )
          if ( .not. isUnneededCell( nx, ny, cflag, includeGhostCells, nix, niy ) ) then
             ggd_grid%space(SPACE_POLOIDALPLANE) % objects(2) % neighbour( ifc, 2, 1 ) = &
                  & gmap % mapFcI( nix, niy, gmap % mapFcIFace( ifc ) )
          end if
          ! measure: area
          if (present(gs)) ggd_grid%space(SPACE_POLOIDALPLANE) % objects(2) % measure( ifc, 1 ) = gs(ix, iy, ALIGNY)*qc(ix, iy)
      end do

      ! 2d objects: cells
      ! ...have four boundaries
      allocate( ggd_grid%space(SPACE_POLOIDALPLANE) % objects(3) % boundary( gmap%ncv, 4) )
      ! 2d object measure: cell volume
      if (present(vol)) allocate( ggd_grid%space(SPACE_POLOIDALPLANE) % objects(3) % measure( gmap % ncv, 1 ) )      
      ! Also store additional geometry information: position in computational space
      ! FIXME: this should go into alternate geometry, which is not available yet for grid objects
      allocate( ggd_grid%space(SPACE_POLOIDALPLANE) % objects(3) % geo(gmap%ncv, 2, 1, 1) )

      ! first set all boundary information to undefined
      ggd_grid%space(SPACE_POLOIDALPLANE) % objects(3) % boundary = GRID_UNDEFINED

      do icv = 1, gmap % ncv
          ix = gmap % mapCvix( icv )
          iy = gmap % mapCviy( icv )

          ! Set position in computational space
          ggd_grid%space(SPACE_POLOIDALPLANE) % objects(3) % geo(icv, 1, 1, 1) = ix
          ggd_grid%space(SPACE_POLOIDALPLANE) % objects(3) % geo(icv, 2, 1, 1) = iy
          ! put faces composing the quadliateral in the list: left face (y-aligned)
          ggd_grid%space(SPACE_POLOIDALPLANE) % objects(3) % boundary( icv, 1 ) = gmap % mapFcI( ix, iy, LEFT )
          ! bottom face (x-aligned
          ggd_grid%space(SPACE_POLOIDALPLANE) % objects(3) % boundary( icv, 2 ) = gmap % mapFcI( ix, iy, BOTTOM )
          ! right face (y-aligned)
          ggd_grid%space(SPACE_POLOIDALPLANE) % objects(3) % boundary( icv, 3 ) = gmap % mapFcI( ix, iy, RIGHT )
          ! top face (x-aligned)
          ggd_grid%space(SPACE_POLOIDALPLANE) % objects(3) % boundary( icv, 4 ) = gmap % mapFcI( ix, iy, TOP )
      end do

      ! Fill in connectivity information
      ! ...have one neighbour per boundary
      allocate( ggd_grid%space(SPACE_POLOIDALPLANE) % objects(3) % neighbour( gmap%ncv, 4, 1) )
      ! first set all to undefined
      ggd_grid%space(SPACE_POLOIDALPLANE) % objects(3) % neighbour = GRID_UNDEFINED

      do icv = 1, gmap % ncv
          ix = gmap % mapCvix( icv )
          iy = gmap % mapCviy( icv )

          do dir = LEFT, TOP
             call getNeighbour(nx, ny, leftix, leftiy, rightix, rightiy, topix, topiy, bottomix, bottomiy, &
                  & ix, iy, dir, nix, niy)             
             if ( .not. isUnneededCell( nx, ny, cflag, includeGhostCells, nix, niy ) ) then
                ggd_grid%space(SPACE_POLOIDALPLANE) % objects(3) % neighbour(icv, dir+1, 1) = gmap % mapCvI( nix, niy )
             end if
          end do

      end do

      ! Fill in x-point indices
      allocate( ggd_grid%space(SPACE_POLOIDALPLANE) % xpoints( gmap % nsv ) )
      ggd_grid%space(SPACE_POLOIDALPLANE) % xpoints = gmap % svi(1:gmap % nsv)

      ! If requested, add a second space for the toroidal angle
      if (SPACE_COUNT == SPACE_TOROIDALANGLE) then
          
          if ( TOROIDAL_PERIODIC ) then
          call gridSetupStruct1dSpace( ggd%spaces(SPACE_TOROIDALANGLE), &
              & COORDTYPE_PHI, &
                  & (/  ( ( 2*pi / NNODES_TOROIDAL ) * i, i = 0, NNODES_TOROIDAL - 1 ) /), &
                  & periodic = .true. )
          else
              call gridSetupStruct1dSpace( ggd%spaces(SPACE_TOROIDALANGLE), &
                  & COORDTYPE_PHI, &
                  & (/  ( ( 2*pi / NNODES_TOROIDAL ) * i, i = 0, NNODES_TOROIDAL ) /) )
          end if

      end if
#endif

    end subroutine fillInGridDescription

#if 0
    ! Part 2: define subgrids
    subroutine fillInSubGridDescription

      ! internal
      integer :: geoId, iRegion, subgridCount, iType, nSubgrid
      integer :: xIn, yIn, xOut, yOut, iCoreSg
      integer :: cls(SPACE_COUNT_MAX)
      integer, allocatable :: xpoints(:,:)

      geoId = geometryId(nnreg, periodic_bc, topcut)
    
      ! Figure out total number of subgrids
      ! Do generic subgrids + subgrids
      nSubgrid = B2_GENERIC_SUBGRID_COUNT + regionCountTotal(geoId)
      ! Inner/outer midplane subgrids
      nSubgrid = nSubgrid + 2

      call logmsg( LOGDEBUG, "b2IMASFillGridDescription: expecting total of "&
          &//Int2str(nSubgrid)//" subgrids" )
      allocate( ggd % subgrids( nSubgrid ) )

      ! Set up generic subgrids

      ! B2_SUBGRID_CELLS: all 2d cells, one implicit object list
      call createSubGridForClass( ggd, ggd % subgrids( B2_SUBGRID_CELLS ), &
          & CLASS_CELL(1:SPACE_COUNT), 'Cells' )

      ! B2_SUBGRID_NODES: all nodes, one implicit object list
      call createSubGridForClass( ggd, ggd % subgrids( B2_SUBGRID_NODES ), &
          & CLASS_NODE(1:SPACE_COUNT), 'Nodes' )

      ! B2_SUBGRID_FACES: all faces, one implicit object list
      call createSubGridForClass( ggd, ggd % subgrids( B2_SUBGRID_FACES ), &
          & CLASS_POLOIDALRADIAL_FACE(1:SPACE_COUNT), 'Faces' )

      ! B2_SUBGRID_FACES_X: x-aligned faces. One implicit object list, range over x faces
      ! Create subgrid with one object list
      call createSubGrid( ggd % subgrids( B2_SUBGRID_FACES_X ), 1, 'x-aligned faces' )
      ! Initialize implicit object list for faces (class (/1/) )
      call createImplicitObjectList( ggd, ggd % subgrids( B2_SUBGRID_FACES_X ) % list(1), &
          & CLASS_POLOIDALRADIAL_FACE(1:SPACE_COUNT) )
      ggd % subgrids( B2_SUBGRID_FACES_X ) % list(1) % indset(1) &
          & = createIndexListForRange( 1, gmap%nfcx )
      if ( SPACE_COUNT == SPACE_TOROIDALANGLE ) then
          ggd % subgrids( B2_SUBGRID_FACES_X ) % list(1) % indset(2) &
              & = createIndexListForRange( 1, 1 )
      end if

      ! B2_SUBGRID_FACES_Y: y-aligned faces. One implicit object list, range over y faces. Same procedure.
      call createSubGrid( ggd % subgrids( B2_SUBGRID_FACES_Y ), 1, 'y-aligned faces' )
      call createImplicitObjectList( ggd, ggd % subgrids( B2_SUBGRID_FACES_Y ) % list(1)&
          & , CLASS_POLOIDALRADIAL_FACE(1:SPACE_COUNT) )
      ggd % subgrids( B2_SUBGRID_FACES_Y ) % list(1) % indset(1) &
          & = createIndexListForRange( gmap%nfcx + 1, gmap%nfcx + gmap%nfcy )
      if ( SPACE_COUNT == SPACE_TOROIDALANGLE ) then
          ggd % subgrids( B2_SUBGRID_FACES_Y ) % list(1) % indset(2) &
              & = createIndexListForRange( 1, 1 )
      end if

      ! Subgrid of all x-points (in one poloidal plane at toroidal index 1)
      ! Assemble object descriptor for x-points
      allocate( xpoints(gmap%nsv, SPACE_COUNT) )
      xpoints = 1
      xpoints(:, SPACE_POLOIDALPLANE) = gmap%svi(1:gmap%nsv)
      call createSubGridForExplicitList( ggd, ggd % subgrids( B2_SUBGRID_XPOINTS ), &
          & CLASS_NODE(1:SPACE_COUNT), xpoints, 'x-points' )

      ! Set up specific subgrids by collection faces for regions

      ! Start counting from end of generic subgrids
      subgridCount = B2_GENERIC_SUBGRID_COUNT

      ! Cell + face subgrids
      do iType = REGIONTYPE_CELL, REGIONTYPE_YFACE
          
          select case(iType)
          case( REGIONTYPE_CELL )
              cls = CLASS_CELL
          case( REGIONTYPE_YFACE, REGIONTYPE_XFACE )
              cls = CLASS_POLOIDALRADIAL_FACE
          end select

          do iRegion = 1, regionCount(geoId, iType)              
              subgridCount = subgridCount + 1

              call logmsg( LOGDEBUG, "b2IMASFillGridDescription: add subgrid #"//int2str(subgridCount)//&
                  & " for iType "//int2str(iType)//&
                  &", iRegion "//int2str(iRegion)//": "//regionName(geoId, iType, iRegion) )

              call createSubGridForExplicitList( ggd, ggd % subgrids( subgridCount ), &
                  & cls(1:SPACE_COUNT), &
                  & collectIndexListForRegion(gmap, region, iType, iRegion), &
                  & regionName(geoId, iType, iRegion) )

          end do
      end do

      ! Add midplane node subgrids
      
      ! Find the core boundary subgrid by looking for its name as defined in b2mod_connectivity
      iCoreSg = gridFindSubGridByName(ggd, "Core boundary")
      ! For double null, we need the outer half of the core boundary
      if (iCoreSg == GRID_UNDEFINED) then 
          iCoreSg = gridFindSubGridByName(ggd, "Outer core boundary")          
      end if
      if (iCoreSg == GRID_UNDEFINED) stop "fillInSubGridDescription: &
          & did not find core boundary subgrid for assembling outer midplane subgrid"

      ! Figure out starting points for inner and outer midplane on core boundary
      call findMidplaneCells(ggd%subgrids(iCoreSg), gmap, crx, xIn, yIn, xOut, yOut)

      subgridCount = subgridCount + 1
      call createSubGridForExplicitList( ggd, ggd % subgrids( subgridCount ), &
          & CLASS_NODE(1:SPACE_COUNT), &
          & collectRadialVertexIndexList(gmap, cflag, xIn, yIn, topix, topiy), &
          & "Inner midplane" )
      
      subgridCount = subgridCount + 1
      call createSubGridForExplicitList( ggd, ggd % subgrids( subgridCount ), &
          & CLASS_NODE(1:SPACE_COUNT), &
          & collectRadialVertexIndexList(gmap, cflag, xOut, yOut, topix, topiy), &
          & "Outer midplane" )      

      call logmsg( LOGDEBUG, "b2IMASFillGridDescription: wrote total of "&
          &//int2str(subgridCount)//" subgrids (expected was "//int2str(size(ggd%subgrids))//')' )

      call assert( subgridCount == size(ggd%subgrids) )
    end subroutine fillInSubGridDescription

#endif

  end subroutine b2IMASFillGridDescription

#else
#ifdef ITM

contains

  !> Routine that fills in a grid description which is part of a CPO
  !> using the given grid data and prepared mappings
  subroutine b2ITMFillGridDescription( gmap,itmgrid, &
      & nx,ny,crx,cry, &
      & leftix,leftiy,rightix,rightiy, &
      & topix,topiy,bottomix,bottomiy,&
      & nnreg,topcut,region,cflag,includeGhostCells,vol,gs,qc )

    type(B2GridMap), intent(in) :: gmap
    type(type_complexgrid), intent(out) :: itmgrid

    ! Size of grid arrays: (-1:nx, -1:ny) 
    integer, intent(in) :: nx, ny
    !   .. output arguments
    ! vertex coordinates
    real (ITM_R8), intent(in) :: &
        & crx(-1:nx,-1:ny,0:3), cry(-1:nx,-1:ny,0:3)
    ! B2 connectivity array
    integer, intent(in) :: &
        & leftix(-1:nx,-1:ny),leftiy(-1:nx,-1:ny),&
        & rightix(-1:nx,-1:ny),rightiy(-1:nx,-1:ny),&
        & topix(-1:nx,-1:ny),topiy(-1:nx,-1:ny),&
        & bottomix(-1:nx,-1:ny),bottomiy(-1:nx,-1:ny)
    ! B2 region & cut information
    integer, intent(in) :: &
        & nnreg(0:2), topcut(:), &
        & region(-1:nx,-1:gmap%b2ny,0:2)
    ! Cell flags 
    integer cflag(-1:nx,-1:ny, CARREOUT_NCELLFLAGS)
    logical, intent(in) :: includeGhostCells
    ! Optional B2 measure information
    real(ITM_R8), intent(in), optional :: vol(-1:nx,-1:ny,0:4), gs(-1:nx,-1:ny,0:2), qc(-1:nx,-1:ny)

    ! internal
    integer, parameter :: NDIM = 2

    call assert( present(gs) .EQV. present(qc) )

    call fillInGridDescription()
    call fillInSubGridDescription()

  contains 

    ! Part 1: fill in grid description
    subroutine fillInGridDescription()

      ! internal
      integer :: ivx, ifc, icv, ix, iy, nix, niy, i, dir

      allocate( itmgrid % spaces(SPACE_COUNT) )

      ! Coordinate types
      ! (dimension of space = NDIM = size( coordtype )
      allocate( itmgrid % spaces(SPACE_POLOIDALPLANE) % coordtype(NDIM, 1) )    
      itmgrid % spaces(SPACE_POLOIDALPLANE) % coordtype(:, 1) &
           & = (/ COORDTYPE_R, COORDTYPE_Z /)

      ! Have two types of objects: 1d edges, 2d cells
      allocate( itmgrid % spaces(SPACE_POLOIDALPLANE) % objects(NDIM + 1) )

      ! Fill in node information
      allocate( itmgrid % spaces(SPACE_POLOIDALPLANE) % objects(1) % geo(gmap%nvx, NDIM, 1, 1) )
      do ivx = 1, gmap % nvx
          itmgrid % spaces(SPACE_POLOIDALPLANE) % objects(1) % geo(ivx, 1, 1, 1) = &
              & crx( gmap % mapVxix( ivx ), gmap % mapVxiy( ivx ), gmap % mapVxIVx( ivx ) )
          itmgrid % spaces(SPACE_POLOIDALPLANE) % objects(1) % geo(ivx, 2, 1, 1) = & 
              & cry( gmap % mapVxix( ivx ), gmap % mapVxiy( ivx ), gmap % mapVxIVx( ivx ) )
      end do

      ! Fill in object definitions (i.e. what objects compose an object)

      ! 1d objects: faces
      ! ...have two boundaries
      allocate( itmgrid % spaces(SPACE_POLOIDALPLANE) % objects(2) % boundary( gmap%nfcx + gmap%nfcy, 2) )
      ! ...have two neighbours, in positive and negative coordinate direction, one on each side
      ! (for x-aligned faces: along flux surface, for y-aligned faces: orthogonal to flux surface)
      allocate( itmgrid % spaces(SPACE_POLOIDALPLANE) % objects(2) % neighbour(gmap%nfcx + gmap%nfcy, 2, 1) )      
      ! 1d object measure: face area
      if (present(gs)) allocate( itmgrid % spaces(SPACE_POLOIDALPLANE) % objects(2) % measure( gmap % nfcx + gmap % nfcy, 1 ) )
      ! first set all boundary & connectivity information to undefined
      itmgrid % spaces(SPACE_POLOIDALPLANE) % objects(2) % boundary = GRID_UNDEFINED
      itmgrid % spaces(SPACE_POLOIDALPLANE) % objects(2) % neighbour = GRID_UNDEFINED
      ! x-aligned faces    
      do ifc = 1, gmap % nfcx    
          ! get position of this face in the b2 grid
          ix = gmap % mapFcix( ifc )
          iy = gmap % mapFciy( ifc )
          ! get index of start vertex 
          ! objdef dims: index of face, 1=start node, 1=one-dimensional object
          select case ( gmap % mapFcIFace( ifc ) )
          case( BOTTOM )
             ! start index: 1=start node
             itmgrid % spaces(SPACE_POLOIDALPLANE) % objects(2) % boundary( ifc, 1 ) = gmap % mapVxI( ix, iy, VX_LOWERLEFT )
             if (gmap % mapVxI( ix, iy, VX_LOWERLEFT ) == GRID_UNDEFINED) then
                call logmsg(LOGWARNING, "b2ITMFillGD: BOTTOM face at "//int2str(ix)//","//int2str(iy)//" has no start node")
             end if
             ! end vertex: 2=end node
             itmgrid % spaces(SPACE_POLOIDALPLANE) % objects(2) % boundary( ifc, 2 ) = gmap % mapVxI( ix, iy, VX_LOWERRIGHT )
             if (gmap % mapVxI( ix, iy, VX_LOWERRIGHT ) == GRID_UNDEFINED) then
                call logmsg(LOGWARNING, "b2ITMFillGD: BOTTOM face at "//int2str(ix)//","//int2str(iy)//" has no end node")
             end if
          case( TOP )
             ! start index: 1=start node
             itmgrid % spaces(SPACE_POLOIDALPLANE) % objects(2) % boundary( ifc, 1 ) = gmap % mapVxI( ix, iy, VX_UPPERLEFT )
             if (gmap % mapVxI( ix, iy, VX_UPPERLEFT ) == GRID_UNDEFINED) then
                call logmsg(LOGWARNING, "b2ITMFillGD: TOP face at "//int2str(ix)//","//int2str(iy)//" has no start node")
             end if
             ! end vertex: 2=end node
             itmgrid % spaces(SPACE_POLOIDALPLANE) % objects(2) % boundary( ifc, 2 ) = gmap % mapVxI( ix, iy, VX_UPPERRIGHT )
             if (gmap % mapVxI( ix, iy, VX_UPPERRIGHT ) == GRID_UNDEFINED) then
                call logmsg(LOGWARNING, "b2ITMFillGD: TOP face at "//int2str(ix)//","//int2str(iy)//" has no end node")
             end if
          end select

          ! Neighbour faces of this face
          ! Left neighbour: face continuing to the left of this face
          nix = leftix( ix, iy )
          niy = leftiy( ix, iy )
          if ( .not. isUnneededCell( nx, ny, cflag, includeGhostCells, nix, niy ) ) then
             itmgrid % spaces(SPACE_POLOIDALPLANE) % objects(2) % neighbour( ifc, 1, 1 ) = &
                  & gmap % mapFcI( nix, niy, gmap % mapFcIFace( ifc ) )
          end if
          ! Right neighbour: face continuing to the right of this face 
          nix = rightix( ix, iy )
          niy = rightiy( ix, iy )
          if ( .not. isUnneededCell( nx, ny, cflag, includeGhostCells, nix, niy ) ) then
             itmgrid % spaces(SPACE_POLOIDALPLANE) % objects(2) % neighbour( ifc, 2, 1 ) = &
                  & gmap % mapFcI( nix, niy, gmap % mapFcIFace( ifc ) )
          end if

          ! measure: area
          if (present(gs)) itmgrid % spaces(SPACE_POLOIDALPLANE) % objects(2) % measure( ifc, 1 ) = gs(ix, iy, ALIGNX)
      end do

      ! y-aligned faces    
      do ifc = gmap % nfcx + 1, gmap % nfcx + gmap % nfcy    
          ! get position of this face in the b2 grid
          ix = gmap % mapFcix( ifc )
          iy = gmap % mapFciy( ifc )
!!$          if (gmap%mapCvI(ix, iy) == GRID_UNDEFINED) then 
!!$                  call logmsg(LOGWARNING, "b2ITMFillGD: writing out faces for unused cell "//int2str(ix)//","//int2str(iy))
!!$          end if

          select case ( gmap % mapFcIFace( ifc ) )
          case( LEFT )
              ! start index: 1=start node
              itmgrid % spaces(SPACE_POLOIDALPLANE) % objects(2) % boundary( ifc, 1 ) = gmap % mapVxI( ix, iy, VX_LOWERLEFT )
              if (gmap % mapVxI( ix, iy, VX_LOWERLEFT ) == GRID_UNDEFINED) then
                  call logmsg(LOGWARNING, "b2ITMFillGD: LEFT face at "//int2str(ix)//","//int2str(iy)//" has no start node")
              end if
          ! end vertex: 2=end node
              itmgrid % spaces(SPACE_POLOIDALPLANE) % objects(2) % boundary( ifc, 2 ) = gmap % mapVxI( ix, iy, VX_UPPERLEFT )
              if (gmap % mapVxI( ix, iy, VX_UPPERLEFT ) == GRID_UNDEFINED) then
                  call logmsg(LOGWARNING, "b2ITMFillGD: LEFT face at "//int2str(ix)//","//int2str(iy)//" has no end node")
          end if
              !if (itmgrid % spaces(SPACE_POLOIDALPLANE) % objects(2) % boundary( ifc, 1 )


          case( RIGHT )
              ! start index: 1=start node
              itmgrid % spaces(SPACE_POLOIDALPLANE) % objects(2) % boundary( ifc, 1 ) = gmap % mapVxI( ix, iy, VX_LOWERRIGHT )
              if (gmap % mapVxI( ix, iy, VX_LOWERRIGHT ) == GRID_UNDEFINED) then
                  call logmsg(LOGWARNING, "b2ITMFillGD: RIGHT face at "//int2str(ix)//","//int2str(iy)//" has no start node")
          end if
              ! end vertex: 2=end node
              itmgrid % spaces(SPACE_POLOIDALPLANE) % objects(2) % boundary( ifc, 2 ) = gmap % mapVxI( ix, iy, VX_UPPERRIGHT )
              if (gmap % mapVxI( ix, iy, VX_UPPERRIGHT ) == GRID_UNDEFINED) then
                  call logmsg(LOGWARNING, "b2ITMFillGD: RIGHT face at "//int2str(ix)//","//int2str(iy)//" has no end node")
              end if
          end select


          ! Neighbour faces of this face
          ! Bottom neighbour: face continuing to the bottom of this face
          nix = bottomix( ix, iy )
          niy = bottomiy( ix, iy )
          if ( .not. isUnneededCell( nx, ny, cflag, includeGhostCells, nix, niy ) ) then
             itmgrid % spaces(SPACE_POLOIDALPLANE) % objects(2) % neighbour( ifc, 1, 1 ) = &
                  & gmap % mapFcI( nix, niy, gmap % mapFcIFace( ifc ) )
          end if
          ! Top neighbour: face continuing to the top of this face 
          nix = topix( ix, iy )
          niy = topiy( ix, iy )
          if ( .not. isUnneededCell( nx, ny, cflag, includeGhostCells, nix, niy ) ) then
             itmgrid % spaces(SPACE_POLOIDALPLANE) % objects(2) % neighbour( ifc, 2, 1 ) = &
                  & gmap % mapFcI( nix, niy, gmap % mapFcIFace( ifc ) )
          end if
          ! measure: area
          if (present(gs)) itmgrid % spaces(SPACE_POLOIDALPLANE) % objects(2) % measure( ifc, 1 ) = gs(ix, iy, ALIGNY)*qc(ix, iy)
      end do

      ! 2d objects: cells
      ! ...have four boundaries
      allocate( itmgrid % spaces(SPACE_POLOIDALPLANE) % objects(3) % boundary( gmap%ncv, 4) )
      ! 2d object measure: cell volume
      if (present(vol)) allocate( itmgrid % spaces(SPACE_POLOIDALPLANE) % objects(3) % measure( gmap % ncv, 1 ) )      
      ! Also store additional geometry information: position in computational space
      ! FIXME: this should go into alternate geometry, which is not available yet for grid objects
      allocate( itmgrid % spaces(SPACE_POLOIDALPLANE) % objects(3) % geo(gmap%ncv, 2, 1, 1) )

      ! first set all boundary information to undefined
      itmgrid % spaces(SPACE_POLOIDALPLANE) % objects(3) % boundary = GRID_UNDEFINED

      do icv = 1, gmap % ncv
          ix = gmap % mapCvix( icv )
          iy = gmap % mapCviy( icv )

          ! Set position in computational space
          itmgrid % spaces(SPACE_POLOIDALPLANE) % objects(3) % geo(icv, 1, 1, 1) = ix
          itmgrid % spaces(SPACE_POLOIDALPLANE) % objects(3) % geo(icv, 2, 1, 1) = iy
          ! put faces composing the quadliateral in the list: left face (y-aligned)
          itmgrid % spaces(SPACE_POLOIDALPLANE) % objects(3) % boundary( icv, 1 ) = gmap % mapFcI( ix, iy, LEFT )
          ! bottom face (x-aligned
          itmgrid % spaces(SPACE_POLOIDALPLANE) % objects(3) % boundary( icv, 2 ) = gmap % mapFcI( ix, iy, BOTTOM )
          ! right face (y-aligned)
          itmgrid % spaces(SPACE_POLOIDALPLANE) % objects(3) % boundary( icv, 3 ) = gmap % mapFcI( ix, iy, RIGHT )
          ! top face (x-aligned)
          itmgrid % spaces(SPACE_POLOIDALPLANE) % objects(3) % boundary( icv, 4 ) = gmap % mapFcI( ix, iy, TOP )
      end do

      ! Fill in connectivity information
      ! ...have one neighbour per boundary
      allocate( itmgrid % spaces(SPACE_POLOIDALPLANE) % objects(3) % neighbour( gmap%ncv, 4, 1) )
      ! first set all to undefined
      itmgrid % spaces(SPACE_POLOIDALPLANE) % objects(3) % neighbour = GRID_UNDEFINED

      do icv = 1, gmap % ncv
          ix = gmap % mapCvix( icv )
          iy = gmap % mapCviy( icv )

          do dir = LEFT, TOP
             call getNeighbour(nx, ny, leftix, leftiy, rightix, rightiy, topix, topiy, bottomix, bottomiy, &
                  & ix, iy, dir, nix, niy)             
             if ( .not. isUnneededCell( nx, ny, cflag, includeGhostCells, nix, niy ) ) then
                itmgrid % spaces(SPACE_POLOIDALPLANE) % objects(3) % neighbour(icv, dir+1, 1) = gmap % mapCvI( nix, niy )
             end if
          end do

      end do

      ! Fill in x-point indices
      allocate( itmgrid % spaces(SPACE_POLOIDALPLANE) % xpoints( gmap % nsv ) )
      itmgrid % spaces(SPACE_POLOIDALPLANE) % xpoints = gmap % svi(1:gmap % nsv)

      ! If requested, add a second space for the toroidal angle
      if (SPACE_COUNT == SPACE_TOROIDALANGLE) then
          
          if ( TOROIDAL_PERIODIC ) then
          call gridSetupStruct1dSpace( itmgrid%spaces(SPACE_TOROIDALANGLE), &
              & COORDTYPE_PHI, &
                  & (/  ( ( 2*pi / NNODES_TOROIDAL ) * i, i = 0, NNODES_TOROIDAL - 1 ) /), &
                  & periodic = .true. )
          else
              call gridSetupStruct1dSpace( itmgrid%spaces(SPACE_TOROIDALANGLE), &
                  & COORDTYPE_PHI, &
                  & (/  ( ( 2*pi / NNODES_TOROIDAL ) * i, i = 0, NNODES_TOROIDAL ) /) )
          end if

      end if

    end subroutine fillInGridDescription


    ! Part 2: define subgrids
    subroutine fillInSubGridDescription

      ! internal
      integer :: geoId, iRegion, subgridCount, iType, nSubgrid
      integer :: xIn, yIn, xOut, yOut, iCoreSg
      integer :: cls(SPACE_COUNT_MAX)
      integer, allocatable :: xpoints(:,:)

      geoId = geometryId(nnreg, periodic_bc, topcut)
    
      ! Figure out total number of subgrids
      ! Do generic subgrids + subgrids
      nSubgrid = B2_GENERIC_SUBGRID_COUNT + regionCountTotal(geoId)
      ! Inner/outer midplane subgrids
      nSubgrid = nSubgrid + 2

      call logmsg( LOGDEBUG, "b2ITMFillGridDescription: expecting total of "&
          &//Int2str(nSubgrid)//" subgrids" )
      allocate( itmgrid % subgrids( nSubgrid ) )

      ! Set up generic subgrids

      ! B2_SUBGRID_CELLS: all 2d cells, one implicit object list
      call createSubGridForClass( itmgrid, itmgrid % subgrids( B2_SUBGRID_CELLS ), &
          & CLASS_CELL(1:SPACE_COUNT), 'Cells' )

      ! B2_SUBGRID_NODES: all nodes, one implicit object list
      call createSubGridForClass( itmgrid, itmgrid % subgrids( B2_SUBGRID_NODES ), &
          & CLASS_NODE(1:SPACE_COUNT), 'Nodes' )

      ! B2_SUBGRID_FACES: all faces, one implicit object list
      call createSubGridForClass( itmgrid, itmgrid % subgrids( B2_SUBGRID_FACES ), &
          & CLASS_POLOIDALRADIAL_FACE(1:SPACE_COUNT), 'Faces' )

      ! B2_SUBGRID_FACES_X: x-aligned faces. One implicit object list, range over x faces
      ! Create subgrid with one object list
      call createSubGrid( itmgrid % subgrids( B2_SUBGRID_FACES_X ), 1, 'x-aligned faces' )
      ! Initialize implicit object list for faces (class (/1/) )
      call createImplicitObjectList( itmgrid, itmgrid % subgrids( B2_SUBGRID_FACES_X ) % list(1), &
          & CLASS_POLOIDALRADIAL_FACE(1:SPACE_COUNT) )
      itmgrid % subgrids( B2_SUBGRID_FACES_X ) % list(1) % indset(1) &
          & = createIndexListForRange( 1, gmap%nfcx )
      if ( SPACE_COUNT == SPACE_TOROIDALANGLE ) then
          itmgrid % subgrids( B2_SUBGRID_FACES_X ) % list(1) % indset(2) &
              & = createIndexListForRange( 1, 1 )
      end if

      ! B2_SUBGRID_FACES_Y: y-aligned faces. One implicit object list, range over y faces. Same procedure.
      call createSubGrid( itmgrid % subgrids( B2_SUBGRID_FACES_Y ), 1, 'y-aligned faces' )
      call createImplicitObjectList( itmgrid, itmgrid % subgrids( B2_SUBGRID_FACES_Y ) % list(1)&
          & , CLASS_POLOIDALRADIAL_FACE(1:SPACE_COUNT) )
      itmgrid % subgrids( B2_SUBGRID_FACES_Y ) % list(1) % indset(1) &
          & = createIndexListForRange( gmap%nfcx + 1, gmap%nfcx + gmap%nfcy )
      if ( SPACE_COUNT == SPACE_TOROIDALANGLE ) then
          itmgrid % subgrids( B2_SUBGRID_FACES_Y ) % list(1) % indset(2) &
              & = createIndexListForRange( 1, 1 )
      end if

      ! Subgrid of all x-points (in one poloidal plane at toroidal index 1)
      ! Assemble object descriptor for x-points
      allocate( xpoints(gmap%nsv, SPACE_COUNT) )
      xpoints = 1
      xpoints(:, SPACE_POLOIDALPLANE) = gmap%svi(1:gmap%nsv)
      call createSubGridForExplicitList( itmgrid, itmgrid % subgrids( B2_SUBGRID_XPOINTS ), &
          & CLASS_NODE(1:SPACE_COUNT), xpoints, 'x-points' )

      ! Set up specific subgrids by collection faces for regions

      ! Start counting from end of generic subgrids
      subgridCount = B2_GENERIC_SUBGRID_COUNT

      ! Cell + face subgrids
      do iType = REGIONTYPE_CELL, REGIONTYPE_YFACE
          
          select case(iType)
          case( REGIONTYPE_CELL )
              cls = CLASS_CELL
          case( REGIONTYPE_YFACE, REGIONTYPE_XFACE )
              cls = CLASS_POLOIDALRADIAL_FACE
          end select

          do iRegion = 1, regionCount(geoId, iType)              
              subgridCount = subgridCount + 1

              call logmsg( LOGDEBUG, "b2ITMFillGridDescription: add subgrid #"//int2str(subgridCount)//&
                  & " for iType "//int2str(iType)//&
                  &", iRegion "//int2str(iRegion)//": "//regionName(geoId, iType, iRegion) )

              call createSubGridForExplicitList( itmgrid, itmgrid % subgrids( subgridCount ), &
                  & cls(1:SPACE_COUNT), &
                  & collectIndexListForRegion(gmap, region, iType, iRegion), &
                  & regionName(geoId, iType, iRegion) )

          end do
      end do

      ! Add midplane node subgrids
      
      ! Find the core boundary subgrid by looking for its name as defined in b2mod_connectivity
      iCoreSg = gridFindSubGridByName(itmgrid, "Core boundary")
      ! For double null, we need the outer half of the core boundary
      if (iCoreSg == GRID_UNDEFINED) then 
          iCoreSg = gridFindSubGridByName(itmgrid, "Outer core boundary")          
      end if
      if (iCoreSg == GRID_UNDEFINED) stop "fillInSubGridDescription: &
          & did not find core boundary subgrid for assembling outer midplane subgrid"

      ! Figure out starting points for inner and outer midplane on core boundary
      call findMidplaneCells(itmgrid%subgrids(iCoreSg), gmap, crx, xIn, yIn, xOut, yOut)

      subgridCount = subgridCount + 1
      call createSubGridForExplicitList( itmgrid, itmgrid % subgrids( subgridCount ), &
          & CLASS_NODE(1:SPACE_COUNT), &
          & collectRadialVertexIndexList(gmap, cflag, xIn, yIn, topix, topiy), &
          & "Inner midplane" )
      
      subgridCount = subgridCount + 1
      call createSubGridForExplicitList( itmgrid, itmgrid % subgrids( subgridCount ), &
          & CLASS_NODE(1:SPACE_COUNT), &
          & collectRadialVertexIndexList(gmap, cflag, xOut, yOut, topix, topiy), &
          & "Outer midplane" )      

      call logmsg( LOGDEBUG, "b2ITMFillGridDescription: wrote total of "&
          &//int2str(subgridCount)//" subgrids (expected was "//int2str(size(itmgrid%subgrids))//')' )

      call assert( subgridCount == size(itmgrid%subgrids) )
    end subroutine fillInSubGridDescription

  end subroutine b2ITMFillGridDescription


  ! Figure out starting cells for inner and outer midplane on core boundary
  ! by finding the points on the core boundary with minimum and maximum r positions
  subroutine findMidplaneCells(coreBndSubgrid, gmap, crx, xIn, yIn, xOut, yOut)   
    type(type_complexgrid_subgrid), intent(in) :: coreBndSubgrid
    type(B2GridMap), intent(in) :: gmap
    ! x/radial vertex coordinates
    real (ITM_R8), intent(in) :: crx(-1:gmap%b2nx,-1:gmap%b2ny,0:3)
    integer, intent(out) :: xIn, yIn, xOut, yOut


    ! internal
    real(ITM_R8) :: rMin, rMax
    type(GridObject) :: obj
    integer :: ix, iy, iObj

    rMin = huge(rMin)
    rMax = -huge(rMax)

    xIn = huge(xIn)
    xOut = huge(xOut)

    ! Loop over all faces in core boundary subgrid
    do iObj = 1, gridSubGridSize(coreBndSubgrid)    
        obj = subGridGetObject(coreBndSubgrid, iObj)
        ! Expect a face
        call assert( all(obj%cls(1:SPACE_COUNT) == CLASS_POLOIDALRADIAL_FACE(1:SPACE_COUNT)) )        
        ! ...which is aligned along the x-direction
        call assert( gmap % mapFcIFace(obj%ind(SPACE_POLOIDALPLANE)) == BOTTOM )
        ix = gmap % mapFcix( obj%ind(SPACE_POLOIDALPLANE) )
        iy = gmap % mapFciy( obj%ind(SPACE_POLOIDALPLANE) )
        
        ! We want the vertex associated with the cell at ix, iy, which is number 0
        if ( crx(ix, iy, 0) < rMin ) then
            rMin = crx(ix, iy, 0)
            xIn = ix
            yIn = iy
        end if
        if ( crx(ix, iy, 0) > rMax ) then
            rMax = crx(ix, iy, 0)
            xOut = ix
            yOut = iy
        end if
    end do

    call assert(xIn /= huge(xIn), "findMidplaneCells: did not find inner midplane position")
    call assert(xOut /= huge(xOut), "findMidplaneCells: did not find outer midplane position")
  end subroutine


  !> Collect the grid indices of all vertices on the radial grid line outward 
  !> starting at the vertex at position six,siy in computational space.
  function collectRadialVertexIndexList(gmap, cflag, six, siy, topix, topiy) result( indexList )
    integer, allocatable, dimension(:,:) :: indexList

    type(B2GridMap), intent(in) :: gmap
    integer, intent(in) :: six, siy
    integer, intent(in) ::  cflag(-1:gmap%b2nx,-1:gmap%b2ny, CARREOUT_NCELLFLAGS)

    ! B2 connectivity array
    integer, intent(in) :: &
        & topix(-1:gmap%b2nx,-1:gmap%b2ny),topiy(-1:gmap%b2nx,-1:gmap%b2ny)

    ! internal
    integer :: ix, iy, nix, niy, nVx, iVx

    ! First figure out how many points we have: start at six, siy, 
    ! go towards top until running out of physical domain
    nVx = 1
    ix = six
    iy = siy
    do
        ! Take a step upwards
        nix = topix(ix, iy)
        niy = topiy(ix, iy)
        ! Stepped outside grid or into ghost cell?
        if (isUnneededCell( gmap%b2nx, gmap%b2ny, cflag, .true., nix, niy) ) then
            exit
        else
            nVx = nVx + 1
            ix = nix
            iy = niy
        end if
    end do
    
    allocate( indexList(nVx, SPACE_COUNT) )   
    if (SPACE_COUNT == SPACE_TOROIDALANGLE)&
        & indexList(:,SPACE_TOROIDALANGLE) = 1

    ! collect indices: repeat above loop with storing indices
    ix = six
    iy = siy
    ! Store starting point index
    iVx = 1    
    indexList(iVx, SPACE_POLOIDALPLANE) = gmap%mapVxI(ix, iy, VX_LOWERLEFT)
    do
        ! take a step
        nix = topix(ix, iy)
        niy = topiy(ix, iy)

        ! Stepped outside grid?
        if (isUnneededCell( gmap%b2nx, gmap%b2ny, cflag, .true., nix, niy) ) then          
            exit
        end if

        ! Store index for new point
        iVx = iVx + 1
        if (gmap%mapVxI(nix, niy, VX_LOWERLEFT) /= GRID_UNDEFINED) then
           indexList(iVx, SPACE_POLOIDALPLANE) = gmap%mapVxI(nix, niy, VX_LOWERLEFT)
        else if (gmap%mapVxI(ix, iy, VX_UPPERLEFT) /= GRID_UNDEFINED &
             & .and. iVx == nVx) then
           indexList(iVx, SPACE_POLOIDALPLANE) = gmap%mapVxI(ix, iy, VX_UPPERLEFT)
        else
           stop "collectRadialVertexIndexList: cannot find expected vertex index"
        end if

        ix = nix
        iy = niy
    end do

    call assert( iVx == nVx )

  end function collectRadialVertexIndexList


  !> Build an index list of all objects of a given region type (b2mod_connectivity.REGIONTYPE_*)
  !> for a given region id.
  !> @param gmap the B2<->CPO grid map, as built by b2CreateMap
  !> @param region the B2 region array
  !> @param iRegionType The region type
  !> @param iRegion The region 
  !> @result The list of indices for all objects that constitute this grid region. The array
  !>  has two dimensions because it is given as a list of object descriptors.
  function collectIndexListForRegion(gmap, region, iRegionType, iRegion) result( indexList )
    integer, allocatable, dimension(:,:) :: indexList

    type(B2GridMap), intent(in) :: gmap
    integer, intent(in) :: region(-1:gmap%b2nx,-1:gmap%b2ny,0:2)
    integer, intent(in) :: iRegionType, iRegion

    ! internal
    integer :: ix, iy, nInd, iInd, ind

    ! Figure out how many indices to expect. A simple count of the form
    ! nInd = count( region(:,:,iRegionType) == iRegion )
    ! will not do, because we have to account for removed objects (ghost cells/faces).

    ! search the relevant objects and count them
    nInd = 0
    do ix = -1, gmap%b2nx
        do iy = -1, gmap%b2ny

            if ( region(ix, iy, iRegionType) == iRegion ) then
                ! Get index depending on what object type we're looking at
                select case (iRegionType) 
                case (REGIONTYPE_CELL)
                    ind = gmap%mapCvI(ix, iy)
                case (REGIONTYPE_XFACE)
                    ind = gmap%mapFcI(ix, iy, LEFT)
                case (REGIONTYPE_YFACE)
                    ind = gmap%mapFcI(ix, iy, BOTTOM)
                end select

                ! Only count this index if not undefined
                if ( ind /= GRID_UNDEFINED ) nInd = nInd + 1
            end if

        end do
    end do

    allocate( indexList(nInd, SPACE_COUNT) )
    indexList = 1

    ! search the relevant objects and store their index consecutively
    iInd = 0
    do ix = -1, gmap%b2nx
        do iy = -1, gmap%b2ny

            if ( region(ix, iy, iRegionType) == iRegion ) then
                ! Get index depending on what object type we're looking at
                select case (iRegionType) 
                case (REGIONTYPE_CELL)
                    ind = gmap%mapCvI(ix, iy)
                case (REGIONTYPE_XFACE)
                    ind = gmap%mapFcI(ix, iy, LEFT)
                case (REGIONTYPE_YFACE)
                    ind = gmap%mapFcI(ix, iy, BOTTOM)
                end select

                if ( ind /= GRID_UNDEFINED ) then
                    iInd = iInd + 1                  
                    call assert(iInd <= nInd)
                    indexList( iInd, SPACE_POLOIDALPLANE ) = ind
                end if
            end if

        end do
    end do

    call assert( iInd == nInd )

  end function collectIndexListForRegion

#endif
#endif

end module b2mod_ual_io_grid

!!!Local Variables:
!!! mode: f90
!!! End:
