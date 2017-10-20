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
    use ids_schemas         ! IGNORE
    use ids_routines        ! IGNORE

    !> IMAS constants definitions (coordinate types identifiers, grid subset 
    !> identifiers ...)
    use ids_grid_common     ! IGNORE
    use ids_string          ! IGNORE
    use ids_grid_subgrid    ! IGNORE
#else
# ifdef ITM
    use itm_types , ITM_R8 => R8, ITM_R4 => R4 ! IGNORE
    Use euITM_schemas ! IGNORE
    use itm_constants , pi => itm_pi ! IGNORE
# endif
#endif
    ! use ggd_common
    use helper
    use logging , only: logmsg, LOGDEBUG
    ! use ggd_assert
    use string
    ! use ggd , GGD_UNDEFINED => GRID_UNDEFINED
    ! use ggd_structured
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

    !> Number of points in toroidal direction (only 1 makes sense here, this is 
    !> for playing around)
    integer, parameter :: NNODES_TOROIDAL = 1

    !> Flag controlling whether the toroidal space is set up periodic or not.
    !> If periodic, the last node is connect to the first by an edge.
    !> If not periodic, an additional node at 2 pi is added
    logical, parameter :: TOROIDAL_PERIODIC = .false.

    !> Object class tuples (ITM classes)
    integer, dimension(SPACE_COUNT_MAX), parameter :: CLASS_NODE = (/ 0, 0 /)
    integer, dimension(SPACE_COUNT_MAX), parameter :: CLASS_RZ_EDGE = (/ 1, 0 /)
    integer, dimension(SPACE_COUNT_MAX), parameter :: CLASS_PHI_EDGE = (/ 0, 1 /)
    integer, dimension(SPACE_COUNT_MAX), parameter ::   &
        &   CLASS_POLOIDALRADIAL_FACE = (/ 1, 1 /)
    integer, dimension(SPACE_COUNT_MAX), parameter ::   &
        &   CLASS_TOROIDAL_FACE = (/ 2, 0 /)
    integer, dimension(SPACE_COUNT_MAX), parameter :: CLASS_CELL = (/ 2, 1 /)  

    !> Object class tuples (one digit IDS classes, transformed from the above 
    !> ITM classes )
    !> Primary IDS classes are:
    !>      Class 1 - nodes/vertices (0D objects)
    !>      Class 2 - edges/faces (1D objects)
    !>      Class 2 - 2D cells (2D objects)
    !> Fortran90 does not allow initialization of constants using SUM. This is 
    !> permitted in newer Fortran 2003. Current workaraund is to directly 
    !> specify the primary IDS class constants
    
    ! integer, parameter :: IDS_CLASS_NODE = sum(CLASS_NODE) + 1 
    integer, parameter :: IDS_CLASS_NODE = 1
    ! integer, parameter :: IDS_CLASS_RZ_EDGE = sum(CLASS_RZ_EDGE) + 1
    integer, parameter :: IDS_CLASS_RZ_EDGE = 2
    ! integer, parameter :: IDS_CLASS_PHI_EDGE = sum(CLASS_PHI_EDGE) + 1
    integer, parameter :: IDS_CLASS_PHI_EDGE = 2
    ! integer, parameter :: IDS_CLASS_POLOIDALRADIAL_FACE =   &
        ! &   sum(CLASS_POLOIDALRADIAL_FACE) 
    integer, parameter :: IDS_CLASS_POLOIDALRADIAL_FACE = 2
    ! integer, parameter :: IDS_CLASS_TOROIDAL_FACE =         &
        ! &   sum(CLASS_TOROIDAL_FACE) 
    integer, parameter :: IDS_CLASS_TOROIDAL_FACE = 2
    ! integer, parameter :: IDS_CLASS_CELL = sum(IDS_CLASS_CELL) 
    integer, parameter :: IDS_CLASS_CELL = 3

    !> Subgrid/Grid subset name constants

    !> Number of generic subgrids
    integer, parameter :: B2_GENERIC_GSUBSET_COUNT = 6

    !> Generic grid subsets (all cells, all faces)
    !> Note: special grid subsets (given by region ids) do not have specific 
    !> constants (see also b2mod_connectivity.f90)

    !> IMAS uses GGD grid subset identifier definitions defined in GSL
    !> (in ids_grid_common)
#ifdef ITM
    !> Unspecified grid subset
    integer, parameter :: B2_GSUBSET_UNSPECIFIED = 0
    !> Grid subset containing all nodes (0D objects) belonging to associated 
    !> space
    integer, parameter :: B2_GSUBSET_NODES = 1
    !> Grid subset containing all faces (first x-aligned, then y-aligned) 
    !> belonging to the associated space (order given by grid map)
    integer, parameter :: B2_GSUBSET_FACES = 2            
    !> Grid subset containing all x-aligned (poloidally) aligned faces belonging
    !>  to the associated space (order given by grid map)
    integer, parameter :: B2_GSUBSET_X_ALIGNED_FACES = 3   
    !> Grid subset containing all y-aligned (radially) aligned faces belonging 
    !> to the associated space (order given by grid map)
    integer, parameter :: B2_GSUBSET_Y_ALIGNED_FACES = 4   
    !> Grid subset containing all cells belonging to the associated space 
    integer, parameter :: B2_GSUBSET_CELLS = 5              
    !> Grid subset containing nodes defining x-points      
    integer, parameter :: B2_GSUBSET_X_POINTS = 6 
#endif           

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
#if 1
        call fillInGridSubsetDescription()
#endif

contains 

    !> Part 1: fill in grid description
    subroutine fillInGridDescription()
        !> Description of some variables:
        !> ifc - Face/edge index 
        !> ivx - Vertex/node index
        !> nfc - Number of all faces/edges (x + y aligned)

        !> internal
        integer ::  ivx, ifc, icv, ix, iy, nix, niy, i, j, dir, nfc
        real(16), parameter :: PI_16 = 4 * atan (1.0_16)

        allocate( ggd_grid%space(SPACE_COUNT) )

        !> Coordinate types
        !> (dimension of space = NDIM = size( coordtype ) 
      
        ! ggd_grid%space(SPACE_POLOIDALPLANE)%geometry_type%name = 'Poloidal' 

        !> Set the space coordinates, also defining the dimension of the space
        allocate( ggd_grid%space( SPACE_POLOIDALPLANE )%coordinates_type(NDIM) )
        ggd_grid%space( SPACE_POLOIDALPLANE )%coordinates_type(1) = COORDTYPE_R
        ggd_grid%space( SPACE_POLOIDALPLANE )%coordinates_type(2) = COORDTYPE_Z

        !> Have two types of objects: 0d nodes, 1d edges, 2d cells
        allocate( ggd_grid%space(SPACE_POLOIDALPLANE)%  &
            &   objects_per_dimension(NDIM + 1) )

        !> Allocate the number of objects of each type
        allocate( ggd_grid%space( SPACE_POLOIDALPLANE )%    &
            &   objects_per_dimension(1)%object( (nx+1)*(ny+1) - 1) )   !> nodes
        allocate( ggd_grid%space( SPACE_POLOIDALPLANE )%    &
            &   objects_per_dimension(2)%object( nx*(ny+1)+(nx+1)*ny) ) !> edges
        allocate( ggd_grid%space( SPACE_POLOIDALPLANE )%    &
            &   objects_per_dimension(3)%object( nx*ny ) )              !> cells
      
        !> Fill in node information
        !> Set number of nodes (0D objects)
        ! allocate( ggd_grid%space( SPACE_POLOIDALPLANE )%    &
            ! &   objects_per_dimension(1)%object( gmap%nvx ) )
        do ivx = 1, gmap % nvx
            !> Allocate goometry leaf for each node
            allocate( ggd_grid%space( SPACE_POLOIDALPLANE )%    &
                &   objects_per_dimension(1)%object( ivx )%geometry(2)) 
#if 1       
            !> Set geometry (R and Z coordinates) of each node
            !> The way of writing nodes data identically to the to IDS converted 
            !> CPO 16151/1000 case, available on (written in 11th October 2017) 
            !> ITER HPC custer in directory
            !> /home/ITER/penkod/public/imasdb/solps-iter/3/0
            !> or 
            !> /home/ITER/kosl/public/imasdb/solps-iter/3/0
            !> while on GateWay Marconi
            !> /marconi_work/eufus_gw/work/g2penkod/imasdb/solps-iter/3/0
            !> or 
            !> /marconi_work/eufus_gw/work/g2kosl/imasdb/solps-iter/3/0
            ggd_grid%space( SPACE_POLOIDALPLANE )%objects_per_dimension(1)%     &
                &   object( ivx )%geometry(1) = crx(    gmap%mapVxix( ivx ),    &
                                                    &   gmap%mapVxiy( ivx ),    &
                                                    &   gmap%mapVxIVx( ivx ))       
            ggd_grid%space( SPACE_POLOIDALPLANE )%objects_per_dimension(1)%     &
                &   object( ivx )%geometry(2) = cry(    gmap%mapVxix( ivx ),    &
                                                    &   gmap%mapVxiy( ivx ),    &
                                                    &   gmap%mapVxIVx( ivx ))
#endif

            !> Set additional node index (REQUIRED!)
            allocate( ggd_grid%space( SPACE_POLOIDALPLANE )%    &
                &   objects_per_dimension(1)%object( ivx )%nodes(1))
            ggd_grid%space( SPACE_POLOIDALPLANE )%objects_per_dimension(1)%     &
                &   object( ivx )%nodes(1) = ivx  
        end do

#if 0
        !> Set geometry (R and Z coordinates) of each node
        !> (original b2mod. This node writing order does not work well with the
        !> code used to write 1D objects data (faces/edges). Thats why it is 
        !> currently not used (#if 0) and the the method above is used instead)
        ivx = 0
        do iy = 0, ny-1
            do ix = 0, nx-1
                ivx = ivx + 1 !> Lower left corners
                ggd_grid%space( SPACE_POLOIDALPLANE )%    &
                    &   objects_per_dimension(1)%object(ivx)%geometry(1) =  &
                    &   crx(ix,iy,0)
                ggd_grid%space( SPACE_POLOIDALPLANE )%    &
                    &   objects_per_dimension(1)%object(ivx)%geometry(2) =  &
                    &   cry(ix,iy,0)
                if (ivx .eq. ((nx+1)*(ny+1) - 1) ) exit
                if (ix.eq.nx-1) then
                    ivx = ivx + 1  !> Lower right corners
                    ggd_grid%space( SPACE_POLOIDALPLANE )%    &
                        &   objects_per_dimension(1)%object(ivx)%geometry(1) =  &
                        &   crx(ix,iy,1)
                    ggd_grid%space( SPACE_POLOIDALPLANE )%    &
                        &   objects_per_dimension(1)%object(ivx)%geometry(2) =  &
                        &   cry(ix,iy,1)
                    if (ivx .eq. ((nx+1)*(ny+1) - 1) ) exit
                end if
                if (iy.eq.ny-1) then
                    ivx = ivx + 1  !> Upper left corners
                    ggd_grid%space( SPACE_POLOIDALPLANE )%    &
                        &   objects_per_dimension(1)%object(ivx)%geometry(1) =  &
                        &   crx(ix,iy,2)
                    ggd_grid%space( SPACE_POLOIDALPLANE )%    &
                        &   objects_per_dimension(1)%object(ivx)%geometry(2) =  &
                        &   cry(ix,iy,2)
                    if (ivx .eq. ((nx+1)*(ny+1) - 1) ) exit
                    if (ix.eq.nx-1) then
                        ivx = ivx + 1  !> Upper right corners
                        ggd_grid%space( SPACE_POLOIDALPLANE )%    &
                            &   objects_per_dimension(1)%object(ivx)%   &
                            &   geometry(1) = crx(ix,iy,3)
                        ggd_grid%space( SPACE_POLOIDALPLANE )%    &
                            &   objects_per_dimension(1)%object(ivx)%   &
                            &   geometry(2) = cry(ix,iy,3)
                        if (ivx .eq. ((nx+1)*(ny+1) - 1) ) exit
                    end if
                end if
            end do
        end do
#endif 


        !> Fill in object definitions (i.e. what objects compose an object)
        !> 1D objects: faces/edges
        nfc = gmap%nfcx + gmap%nfcy !> Number of all faces/edges
        write(0,*) "num_obj_1D: ", nfc
        write(0,*) "num_obj_1D x-aligned: ", gmap%nfcx
        write(0,*) "num_obj_1D y-aligned: ", gmap%nfcy
        !> Allocate 1D objects 
        allocate( ggd_grid%space( SPACE_POLOIDALPLANE )%    &
            &   objects_per_dimension(2)%object( gmap%nfcx + gmap%nfcy) )
        !> Each 1D object has two boundaries and two neighbours, in positive 
        !> and negative coordinate direction, one on each side (for x-aligned 
        !> faces: along flux surface, for y-aligned faces: orthogonal to flux 
        !> surface)
        do ifc = 1, nfc
            !> Allocate and set all boundary & connectivity information to 
            !> undefined
            allocate( ggd_grid%space( SPACE_POLOIDALPLANE )%    &
                &   objects_per_dimension(2)%object( ifc )%boundary(2) )
            !> Allocate list of 0D objects forming the 1D object 
            !> Two 0D objects (vertices/nodes) form one 1D object (face/edge)
            allocate( ggd_grid%space( SPACE_POLOIDALPLANE )%    &
                &   objects_per_dimension(2)%object( ifc )%nodes(2) )
            do i = 1, 2
                !> Boundary to undefined
                ggd_grid%space( SPACE_POLOIDALPLANE )%objects_per_dimension(2)% &
                    &   object( ifc )%boundary(i)%index = B2_GRID_UNDEFINED
                !> Neighbours to undefined
                allocate(ggd_grid%space( SPACE_POLOIDALPLANE )% &
                    &   objects_per_dimension(2)%object( ifc )%boundary(i)% &
                    &   neighbours(1))
                ggd_grid%space( SPACE_POLOIDALPLANE )%objects_per_dimension(2)% &
                    &   object( ifc )%boundary(i)%neighbours(1) =   &
                    &   B2_GRID_UNDEFINED
                !> Nodes to undefined
                ggd_grid%space( SPACE_POLOIDALPLANE )%objects_per_dimension(2)% &
                    &   object( ifc )%nodes(i) = B2_GRID_UNDEFINED
            end do
        end do

        !> x-aligned faces    
        do ifc = 1, gmap%nfcx   
            !> get position of this face in the b2 grid
            ix = gmap%mapFcix( ifc )
            iy = gmap%mapFciy( ifc )
            !> get index of start vertex 
            !> objdef dims: index of face, 1=start node, 1=one-dimensional 
            !> object
            select case ( gmap%mapFcIFace( ifc ) )
            case( BOTTOM )
                !> start index: 1=start node
                ggd_grid%space( SPACE_POLOIDALPLANE )%objects_per_dimension(2)% &
                    &   object( ifc )%boundary(1)%index =     &
                    &   gmap%mapVxI( ix, iy, VX_LOWERLEFT )
                ggd_grid%space( SPACE_POLOIDALPLANE )%objects_per_dimension(2)% &
                    &   object( ifc )%nodes(1) = gmap%mapVxI( ix, iy, VX_LOWERLEFT )
                if (gmap%mapVxI( ix, iy, VX_LOWERLEFT ) ==  &
                    &   B2_GRID_UNDEFINED) then
                    call logmsg(LOGWARNING, "b2IMASFillGD: BOTTOM face at " &
                        &   //idsInt2str(ix)//","//idsInt2str(iy)//         &
                        &   " has no start node")
                end if
                !> end vertex: 2=end node
                ggd_grid%space( SPACE_POLOIDALPLANE )%objects_per_dimension(2)% &
                    &   object( ifc )%boundary(2)%index =     &
                    &   gmap%mapVxI( ix, iy, VX_LOWERRIGHT )
                ggd_grid%space( SPACE_POLOIDALPLANE )%objects_per_dimension(2)% &
                    &   object( ifc )%nodes(2) = gmap%mapVxI( ix, iy, VX_LOWERRIGHT )
                if (gmap%mapVxI( ix, iy, VX_LOWERRIGHT ) == &
                    &   B2_GRID_UNDEFINED) then
                    call logmsg(LOGWARNING, "b2IMASFillGD: BOTTOM face at " &
                        &   //idsInt2str(ix)//","//idsInt2str(iy)//         &
                        &   " has no end node")
                end if
            case( TOP )
                !> start index: 1=start node
                ggd_grid%space( SPACE_POLOIDALPLANE )%objects_per_dimension(2)% &
                    &   object( ifc )%boundary(1)%index =     &
                    &   gmap%mapVxI( ix, iy, VX_UPPERLEFT )
                ggd_grid%space( SPACE_POLOIDALPLANE )%objects_per_dimension(2)% &
                    &   object( ifc )%nodes(1) = gmap%mapVxI( ix, iy, VX_UPPERLEFT )
                if (gmap%mapVxI( ix, iy, VX_UPPERLEFT ) ==  &
                    &   B2_GRID_UNDEFINED) then
                    call logmsg(LOGWARNING, "b2IMASFillGD: TOP face at "    &
                        &   //idsInt2str(ix)//","//idsInt2str(iy)//         &
                        &   " has no start node")
                end if
                !> end vertex: 2=end node
                ggd_grid%space( SPACE_POLOIDALPLANE )%objects_per_dimension(2)% &
                    &   object( ifc )%boundary(2)%index =     &
                    &   gmap%mapVxI( ix, iy, VX_UPPERRIGHT )
                ggd_grid%space( SPACE_POLOIDALPLANE )%objects_per_dimension(2)% &
                    &   object( ifc )%nodes(2) = gmap%mapVxI( ix, iy, VX_UPPERRIGHT )
                if (gmap%mapVxI( ix, iy, VX_UPPERRIGHT ) == &
                    &   B2_GRID_UNDEFINED) then
                    call logmsg(LOGWARNING, "b2IMASFillGD: TOP face at "    &
                        &   //idsInt2str(ix)//","//idsInt2str(iy)//         &
                        &   " has no end node")
                end if
            end select

            !> Neighbour faces of this face
            !> Left neighbour: face continuing to the left of this face
            nix = leftix( ix, iy )
            niy = leftiy( ix, iy )
            if ( .not. isUnneededCell( nx, ny, cflag, includeGhostCells, nix, niy ) ) then
                ggd_grid%space( SPACE_POLOIDALPLANE )%objects_per_dimension(2)% &
                    &   object( ifc )%boundary(1)%neighbours(1) =    &
                    &   gmap%mapFcI( nix, niy, gmap%mapFcIFace( ifc ) )
            end if
            !> Right neighbour: face continuing to the right of this face 
            nix = rightix( ix, iy )
            niy = rightiy( ix, iy )
            if ( .not. isUnneededCell( nx, ny, cflag, includeGhostCells, nix, niy ) ) then
                ggd_grid%space( SPACE_POLOIDALPLANE )%objects_per_dimension(2)% &
                    &   object( ifc )%boundary(2)%neighbours(2) =     &
                    &   gmap%mapFcI( nix, niy, gmap%mapFcIFace( ifc ) )
            end if

            !> 1d object measure: face area
            if (present(gs)) then
                ggd_grid%space( SPACE_POLOIDALPLANE )%objects_per_dimension(2)% &
                    &   object( ifc )%measure = &
                    &   gs(ix, iy, ALIGNX)
            end if
        end do

        !> y-aligned faces    
        do ifc = gmap%nfcx + 1, gmap % nfcx + gmap % nfcy   
            !> get position of this face in the b2 grid
            ix = gmap % mapFcix( ifc )
            iy = gmap % mapFciy( ifc )
!!$            if (gmap%mapCvI(ix, iy) == B2_GRID_UNDEFINED) then 
!!$                call logmsg(LOGWARNING, "b2IMASFillGD: writing out faces 
!!$                for unused cell "//idsint2str(ix)//","//idsint2str(iy))
!!$            end if

            select case ( gmap % mapFcIFace( ifc ) )
            case( LEFT )
                !> start index: 1=start node
                ggd_grid%space( SPACE_POLOIDALPLANE )%objects_per_dimension(2)% &
                    &   object( ifc )%boundary(1)%index =                       &
                    &   gmap%mapVxI( ix, iy, VX_LOWERLEFT )
                ggd_grid%space( SPACE_POLOIDALPLANE )%objects_per_dimension(2)% &
                    &   object( ifc )%nodes(1) = gmap%mapVxI( ix, iy, VX_LOWERLEFT )
                if (gmap%mapVxI( ix, iy, VX_LOWERLEFT ) ==  &
                    &   B2_GRID_UNDEFINED) then
                    call logmsg(LOGWARNING, "b2IMASFillGD: LEFT face at "   &
                        &   //idsint2str(ix)//","//idsint2str(iy)//         &
                        &   " has no start node")
                end if
                !> end vertex: 2=end node
                ggd_grid%space( SPACE_POLOIDALPLANE )%objects_per_dimension(2)% &
                    &   object( ifc )%boundary(2)%index =                       &
                    &   gmap%mapVxI( ix, iy, VX_UPPERLEFT )
                ggd_grid%space( SPACE_POLOIDALPLANE )%objects_per_dimension(2)% &
                    &   object( ifc )%nodes(2) = gmap%mapVxI( ix, iy, VX_UPPERLEFT )
                if (gmap%mapVxI( ix, iy, VX_UPPERLEFT ) ==  &
                    &   B2_GRID_UNDEFINED) then
                    call logmsg(LOGWARNING, "b2IMASFillGD: LEFT face at "   &
                        &   //idsint2str(ix)//","//idsint2str(iy)//         &
                        &   " has no end node")
                end if

            case( RIGHT )
                !> start index: 1=start node
                ggd_grid%space( SPACE_POLOIDALPLANE )%objects_per_dimension(2)% &
                    &   object( ifc )%boundary(1)%index =                       &
                    &   gmap%mapVxI( ix, iy, VX_LOWERRIGHT )
                ggd_grid%space( SPACE_POLOIDALPLANE )%objects_per_dimension(2)% &
                    &   object( ifc )%nodes(1) = gmap%mapVxI( ix, iy, VX_LOWERRIGHT )
                if (gmap%mapVxI( ix, iy, VX_LOWERRIGHT ) == &
                    &   B2_GRID_UNDEFINED) then
                    call logmsg(LOGWARNING, "b2IMASFillGD: RIGHT face at "  &
                        &   //idsint2str(ix)//","//idsint2str(iy)//         &
                        &   " has no start node")
                end if
                !> end vertex: 2=end node
                ggd_grid%space( SPACE_POLOIDALPLANE )%objects_per_dimension(2)% &
                    &   object( ifc )%boundary(2)%index =                       &
                    &   gmap%mapVxI( ix, iy, VX_UPPERRIGHT )
                ggd_grid%space( SPACE_POLOIDALPLANE )%objects_per_dimension(2)% &
                    &   object( ifc )%nodes(2) = gmap%mapVxI( ix, iy, VX_UPPERRIGHT )
                if (gmap%mapVxI( ix, iy, VX_UPPERRIGHT ) == &
                    &    B2_GRID_UNDEFINED) then
                    call logmsg(LOGWARNING, "b2IMASFillGD: RIGHT face at "  &
                        &   //idsint2str(ix)//","//idsint2str(iy)//         &
                        &   " has no end node")
                end if
            end select

            !> Neighbour faces of this face
            !> Bottom neighbour: face continuing to the bottom of this face
            nix = bottomix( ix, iy )
            niy = bottomiy( ix, iy )
            if ( .not. isUnneededCell( nx, ny, cflag, includeGhostCells, nix, niy ) ) then
                ggd_grid%space( SPACE_POLOIDALPLANE )%objects_per_dimension(2)% &
                    &   object( ifc )%boundary(1)%neighbours(1) =    &
                    &   gmap%mapFcI( nix, niy, gmap%mapFcIFace( ifc ) )
            end if
            !> Top neighbour: face continuing to the top of this face 
            nix = topix( ix, iy )
            niy = topiy( ix, iy )
            if ( .not. isUnneededCell( nx, ny, cflag, includeGhostCells, nix, niy ) ) then
                ggd_grid%space( SPACE_POLOIDALPLANE )%objects_per_dimension(2)% &
                    &   object( ifc )%boundary(1)%neighbours(1) =    &
                    &   gmap%mapFcI( nix, niy, gmap%mapFcIFace( ifc ) )
            end if
            !> 1d object measure: edge length
            if (present(gs)) then
                ggd_grid%space( SPACE_POLOIDALPLANE )%objects_per_dimension(2)% &
                    &   object( ifc )%measure = &
                    &   gs(ix, iy, ALIGNY)*qc(ix, iy) 
            end if
        end do

        !> Fill in object definitions (i.e. what objects compose an object)
        !> 2D objects: Cells
        write(0,*) "num_obj_2D: ", gmap%ncv
        !> Allocate 2D objects
        allocate( ggd_grid%space(SPACE_POLOIDALPLANE)%objects_per_dimension(3)% &
            &   object( gmap%ncv ) )
        
        !> Each 2D object has four boundaries
        !> Each boundary has one neighbour
        do icv = 1, gmap%ncv
            !> Allocate and set all boundary & connectivity information to 
            !> undefined
            allocate( ggd_grid%space(SPACE_POLOIDALPLANE)%  &
                &   objects_per_dimension(3)%object( icv )%boundary(4) )
            !> Allocate list of 0D objects forming the 2D object 
            !> Four 0D objects (vertices/nodes) form one 2D object (2D cell)
            allocate( ggd_grid%space( SPACE_POLOIDALPLANE )%    &
                &   objects_per_dimension(3)%object( icv )%nodes(4) )
            !> Also store additional geometry information: position in 
            !> computational space
            !> FIXME:   this should go into alternate geometry, which is not 
            !>          available yet for grid objects
            allocate( ggd_grid%space(SPACE_POLOIDALPLANE)%  &
                &   objects_per_dimension(3)%object( icv )%geometry(2) )
            do i = 1, 4
                !> Boundary to undefined
                ggd_grid%space( SPACE_POLOIDALPLANE )%objects_per_dimension(3)% &
                    &   object( icv )%boundary(i)%index = B2_GRID_UNDEFINED
                !> Neighbours to undefined
                allocate(ggd_grid%space( SPACE_POLOIDALPLANE )% &
                    &   objects_per_dimension(3)%object( icv )%boundary(i)% &
                    &   neighbours(1))
                ggd_grid%space( SPACE_POLOIDALPLANE )%objects_per_dimension(3)% &
                    &   object( icv )%boundary(i)%neighbours(1) =   &
                    &   B2_GRID_UNDEFINED
                !> Nodes to undefined
                ggd_grid%space( SPACE_POLOIDALPLANE )%objects_per_dimension(3)% &
                    &   object( icv )%nodes(i) = B2_GRID_UNDEFINED
                !> Geometry to undefined
                if (i < 3) then
                    ggd_grid%space( SPACE_POLOIDALPLANE )%objects_per_dimension(3)% &
                        &   object( icv )%geometry(i) = B2_GRID_UNDEFINED
                end if
            end do
        end do

        do icv = 1, gmap%ncv
            !> Set position in computational space
            ix = gmap%mapCvix( icv )
            iy = gmap%mapCviy( icv )
            ggd_grid%space(SPACE_POLOIDALPLANE)%objects_per_dimension(3)%   &
                &   object( icv )%geometry(1) = ix
            ggd_grid%space(SPACE_POLOIDALPLANE)%objects_per_dimension(3)%   &
                &   object( icv )%geometry(2) = iy

            !> Set faces composing the quadliateral in the list: 
            !> left face (y-aligned)
            ggd_grid%space(SPACE_POLOIDALPLANE)%objects_per_dimension(3)%   &
                &   object( icv )%boundary(1)%index = gmap%mapFcI( ix, iy, LEFT )
            !> bottom face (x-aligned
            ggd_grid%space(SPACE_POLOIDALPLANE)%objects_per_dimension(3)%   &
                &   object( icv )%boundary(2)%index = gmap%mapFcI( ix, iy, BOTTOM )
            !> right face (y-aligned)
            ggd_grid%space(SPACE_POLOIDALPLANE)%objects_per_dimension(3)%   &
                &   object( icv )%boundary(3)%index = gmap%mapFcI( ix, iy, RIGHT )
            !> top face (x-aligned)
            ggd_grid%space(SPACE_POLOIDALPLANE)%objects_per_dimension(3)%   &
                &   object( icv )%boundary(4)%index = gmap%mapFcI( ix, iy, TOP )
            do dir = LEFT, TOP
                call getNeighbour(nx, ny, leftix, leftiy, rightix, rightiy,     &
                    &   topix, topiy, bottomix, bottomiy, ix, iy, dir, nix, niy)             
                if ( .not. isUnneededCell( nx, ny, cflag, includeGhostCells,    &
                    &   nix, niy ) ) then
                    ggd_grid%space(SPACE_POLOIDALPLANE)%            &
                        &   objects_per_dimension(3)%object (icv )% &
                        &   boundary( dir + 1 )%neighbours(1) =     &
                        &   gmap%mapCvI( nix, niy )
                end if
            end do
            !> 2d object measure: cell area
            if (present(vol)) then 
                ggd_grid%space(SPACE_POLOIDALPLANE)%objects_per_dimension(3)%   &
                    &   object( icv )%measure = vol(ix, iy, 1)
            end if 
        end do

        !> Set nodes list, composing the 2D objects - Cells, using a subroutine
        call setCellsConnectivityArrayNodes(ggd_grid)

#if 0
        !! TODO
        !> Fill in x-point indices
        !> In edge_profiles no node for data on x-points was found. There is 
        !> hovewer one in equilibrium%boundary%x_point
        allocate( ggd_grid%space(SPACE_POLOIDALPLANE)%xpoints( gmap%nsv ) )
        ggd_grid%space(SPACE_POLOIDALPLANE)%xpoints = gmap%svi(1:gmap%nsv)
#endif

        !> If requested, add a second space for the toroidal angle
        if (SPACE_COUNT == SPACE_TOROIDALANGLE) then
            if ( TOROIDAL_PERIODIC ) then
            call gridSetupStruct1dSpace( ggd_grid%space(SPACE_TOROIDALANGLE), &
                & COORDTYPE_PHI, &
                    & (/  ( ( 2*PI_16 / NNODES_TOROIDAL ) * i, i = 0, NNODES_TOROIDAL - 1 ) /), &
                    & .true. ) !> periodic = .true.
            else
                call gridSetupStruct1dSpace( ggd_grid%space(SPACE_TOROIDALANGLE), &
                    & COORDTYPE_PHI, &
                    & (/  ( ( 2*PI_16 / NNODES_TOROIDAL ) * i, i = 0, NNODES_TOROIDAL ) /) )
            end if
        end if

    end subroutine fillInGridDescription

    subroutine setCellsConnectivityArrayNodes(ggd_grid)
        type(ids_generic_grid_dynamic), intent(inout) :: ggd_grid
        !> internal
        integer, allocatable    ::  objects2Darray(:,:)
        integer ::  num_nodes_2D, num_boundary_2D, node1, node2
        integer, allocatable ::  node_idx(:)
        integer, allocatable    ::  free_edge(:)
        integer ::  edge_idx, last_idx
        integer ::  icv, ix, iy, m, loop_count


        !> Get the list of 0D objects forming the 2D objects and wirte it to IDS
        allocate( objects2Darray( gmap%ncv, 4 ) )
        ! objects2Darray = numpy.array([], dtype='int')

        ! ids_dim_2D.object.resize(ncv)

        ! Already done
        allocate( ggd_grid%space(SPACE_POLOIDALPLANE)%objects_per_dimension(3)% &
                &   object( gmap%ncv ) )

        !> Get 2D objects geometry (nodes) data from CPO, sort them into more 
        !> orderly form  and put them into the IDS
        num_nodes_2D    = 4
        num_boundary_2D = 4
        allocate( node_idx(4) )
        allocate( free_edge(3) )

        do icv = 1, gmap%ncv
            ix = gmap%mapCvix( icv )
            iy = gmap%mapCviy( icv )
            allocate( ggd_grid%space(SPACE_POLOIDALPLANE)%objects_per_dimension(3)% &
                &   object(icv)%nodes(num_nodes_2D) )
            
            edge_idx = gmap%mapFcI( ix, iy, LEFT )
            free_edge(1) = gmap%mapFcI( ix, iy, BOTTOM )
            free_edge(2) = gmap%mapFcI( ix, iy, RIGHT )
            free_edge(3) = gmap%mapFcI( ix, iy, TOP )

            node_idx(1) = ggd_grid%space( SPACE_POLOIDALPLANE )%   &
                &   objects_per_dimension(2)%object( edge_idx )%boundary(1)%index
            last_idx = 2
            node_idx(last_idx) = ggd_grid%space( SPACE_POLOIDALPLANE )%   &
                &   objects_per_dimension(2)%object( edge_idx )%boundary(2)%index

            do loop_count = 1, 4
                if (last_idx < 4) then
                    do m = 1, 3
                        edge_idx = free_edge(m)
                        if (edge_idx < 1) then
                            continue
                        end if
                        node1 = ggd_grid%space( SPACE_POLOIDALPLANE )%   &
                            &   objects_per_dimension(2)%object( edge_idx )% &
                            &   boundary(1)%index
                        node2 = ggd_grid%space( SPACE_POLOIDALPLANE )%   &
                            &   objects_per_dimension(2)%object( edge_idx )% &
                            &   boundary(2)%index

                        if (node_idx(last_idx) == node1) then
                            free_edge(m) = -1
                            last_idx = last_idx + 1
                            node_idx(last_idx)= node2
                            exit
                        end if
                        if (node_idx(last_idx) == node2) then
                            free_edge(m) = -1
                            last_idx = last_idx + 1
                            node_idx(last_idx) = node1
                            exit
                        end if
                    end do
                end if
            end do

            !> Set nodes list, composing the 2D objects - Cells
            ggd_grid%space( SPACE_POLOIDALPLANE )%objects_per_dimension(3)% &
                &   object( icv )%nodes(1) = node_idx(1)
            ggd_grid%space( SPACE_POLOIDALPLANE )%objects_per_dimension(3)% &
                &   object( icv )%nodes(2) = node_idx(4)
            ggd_grid%space( SPACE_POLOIDALPLANE )%objects_per_dimension(3)% &
                &   object( icv )%nodes(3) = node_idx(3)
            ggd_grid%space( SPACE_POLOIDALPLANE )%objects_per_dimension(3)% &
                &   object( icv )%nodes(4) = node_idx(2)
        end do

    end subroutine setCellsConnectivityArrayNodes


    ! Part 2: define grid subsets
    subroutine fillInGridSubsetDescription

        !> internal
        integer :: geoId, iRegion, GSubsetCount, iType, nGSubset
        integer :: xIn, yIn, xOut, yOut, iCoreGS
        integer :: cls(SPACE_COUNT_MAX)
        integer, allocatable :: xpoints(:,:)

        geoId = geometryId(nnreg, periodic_bc, topcut)
    
        !> Figure out total number of subgrids
        !> Do generic subgrids + subgrids
        nGSubset = B2_GENERIC_GSUBSET_COUNT + regionCountTotal(geoId)
        !> Inner/outer midplane subgrids
        ! nGSubset = nGSubset + 2

        call logmsg( LOGDEBUG, "b2IMASFillGridDescription: expecting total of " &
            &//idsInt2str(nGSubset)//" grid subsets" )
        allocate( ggd_grid%grid_subset( nGSubset ) )

        !> Set up generic grid subsets
        !> GRID_SUBSET_NODES: all nodes, one implicit object list
        !> The commented codes use CLASS_NODE, CLASS_POLOIDALRADIAL_FACE and 
        !> CLASS_CELL parameters, using CPO class definitions, which might not 
        !> work with the current IMAS GSL

        call createGridSubsetForClass(  ggd_grid,                   &
            &   ggd_grid%grid_subset( GRID_SUBSET_NODES ),          & 
            &   IDS_CLASS_NODE, 1, GRID_SUBSET_NODES, "Nodes",      &
            &   "All nodes (0D objects) in the domain."  )
        ! call createGridSubsetForClass(  ggd_grid,           &
        !     &   ggd_grid%grid_subset( GRID_SUBSET_NODES ),   &
        !     &   CLASS_NODE(1:SPACE_COUNT), 1, GRID_SUBSET_NODES, "Nodes B2.5",   &
        !     &   "All nodes (0D objects) in the domain."  )

        !> GRID_SUBSET_FACES: all faces, one implicit object list
        call createGridSubsetForClass(  ggd_grid,                       &
            &   ggd_grid%grid_subset( GRID_SUBSET_FACES ),              &
            &   IDS_CLASS_POLOIDALRADIAL_FACE, 1, GRID_SUBSET_FACES,    &
            &   "Faces", "All faces (1D objects) in the domain."  )
        ! call createGridSubsetForClass(  ggd_grid,         &
        !     &   ggd_grid%grid_subset( GRID_SUBSET_FACES ), &
        !     &   CLASS_POLOIDALRADIAL_FACE(1:SPACE_COUNT), 1, &
        !     &   GRID_SUBSET_FACES, "Faces",     &
        !     &   "All faces (1D objects) in the domain."  )

        !> GRID_SUBSET_X_ALIGNED_FACES: x-aligned faces. One implicit object 
        !> list, range over x faces
        !> Create grid subset with one object list
        call createEmptyGridSubset(                                     &
            &   ggd_grid%grid_subset( GRID_SUBSET_X_ALIGNED_FACES ),    &
            &   GRID_SUBSET_X_ALIGNED_FACES, 'x-aligned faces' )
        !> Initialize implicit object list for faces (class (/2/) )
        call createExplicitObjectListSingleSpace( ggd_grid,         &
            &   ggd_grid%grid_subset( GRID_SUBSET_X_ALIGNED_FACES), &
            &   IDS_CLASS_POLOIDALRADIAL_FACE, (/ 1:gmap%nfcx /),   &
            &   IDS_CLASS_POLOIDALRADIAL_FACE, 1)
        ! call createExplicitObjectListSingleSpace( ggd_grid,         &
        !     &   ggd_grid%grid_subset( GRID_SUBSET_X_ALIGNED_FACES),  &
        !     &   GRID_SUBSET_X_ALIGNED_FACES, (/ 1:gmap%nfcx /),      &
        !     &   CLASS_POLOIDALRADIAL_FACE(1:SPACE_COUNT), 1)

        if ( SPACE_COUNT == SPACE_TOROIDALANGLE ) then
            call createExplicitObjectListSingleSpace( ggd_grid,         &
                &   ggd_grid%grid_subset( GRID_SUBSET_X_ALIGNED_FACES), &
                &   IDS_CLASS_POLOIDALRADIAL_FACE, (/ 1:1 /),           &
                &   IDS_CLASS_POLOIDALRADIAL_FACE, 1)
        ! call createExplicitObjectListSingleSpace( ggd_grid,         &
        !     &   ggd_grid%grid_subset( GRID_SUBSET_X_ALIGNED_FACES),  &
        !     &   GRID_SUBSET_X_ALIGNED_FACES, (/ 1:1 /),      &
        !     &   CLASS_POLOIDALRADIAL_FACE(1:SPACE_COUNT), 1)
        end if

        !> GRID_SUBSET_Y_ALIGNED_FACES: y-aligned faces. One implicit object 
        !> list, range over y faces
        !> Create grid subset with one object list
        call createEmptyGridSubset(                                     &
            &   ggd_grid%grid_subset( GRID_SUBSET_Y_ALIGNED_FACES ),    & 
            &   GRID_SUBSET_Y_ALIGNED_FACES, 'y-aligned faces' )
        !> Initialize implicit object list for faces (class (/2/) )
        call createExplicitObjectListSingleSpace( ggd_grid,             &
            &   ggd_grid%grid_subset( GRID_SUBSET_Y_ALIGNED_FACES),     &
            &   IDS_CLASS_POLOIDALRADIAL_FACE,                          &
            &   (/ ( gmap%nfcx + 1 ) : ( gmap%nfcx + gmap%nfcy ) /),    &
            &   IDS_CLASS_POLOIDALRADIAL_FACE, 1)
        ! call createExplicitObjectListSingleSpace( ggd_grid,         &
        !     &   ggd_grid%grid_subset( GRID_SUBSET_Y_ALIGNED_FACES),  &
        !     &   GRID_SUBSET_Y_ALIGNED_FACES, (/ 1:gmap%nfcy /),      &
        !     &   CLASS_POLOIDALRADIAL_FACE(1:SPACE_COUNT), 1)

        if ( SPACE_COUNT == SPACE_TOROIDALANGLE ) then
        call createExplicitObjectListSingleSpace( ggd_grid,         &
            &   ggd_grid%grid_subset( GRID_SUBSET_Y_ALIGNED_FACES), &
            &   IDS_CLASS_POLOIDALRADIAL_FACE, (/ 1:1 /),           &
            &   IDS_CLASS_POLOIDALRADIAL_FACE, 1)
        ! call createExplicitObjectListSingleSpace( ggd_grid,         &
        !     &   ggd_grid%grid_subset( GRID_SUBSET_Y_ALIGNED_FACES),  &
        !     &   GRID_SUBSET_Y_ALIGNED_FACES, (/ 1:1 /),      &
        !     &   CLASS_POLOIDALRADIAL_FACE(1:SPACE_COUNT), 1)
        end if

        !> GRID_SUBSET_CELLS: all 2d cells, one implicit object list
        call createGridSubsetForClass(  ggd_grid,           &
            &   ggd_grid%grid_subset( GRID_SUBSET_CELLS ),  &
            &   IDS_CLASS_CELL, 1, GRID_SUBSET_CELLS, "Cells",  &
            &   "All faces (1D objects) in the domain."  )

        ! call createGridSubsetForClass(  ggd_grid,               &
        !     &   ggd_grid%grid_subset( GRID_SUBSET_CELLS ),       &
        !     &   CLASS_CELL(1:SPACE_COUNT), 1, GRID_SUBSET_CELLS, &
        !     &   "Cells", "All faces (1D objects) in the domain." )

        !> Grid subset of all x-points (in one poloidal plane at toroidal index 1)
        !> Assemble object descriptor for x-points
        allocate( xpoints(gmap%nsv, SPACE_COUNT) )
        xpoints = 1
        xpoints(:, SPACE_POLOIDALPLANE) = gmap%svi(1:gmap%nsv)
        !> Create grid subset with one object list
        call createEmptyGridSubset(                                     &
            &   ggd_grid%grid_subset( GRID_SUBSET_X_POINTS ),           &
            &   GRID_SUBSET_X_POINTS, 'x-points' )
        !> Initialize explicit object list for faces (class (/1/) )
        !> TODO: xpoints(:, 1 ) -> taking values for first space only. Set for 
        !> all spaces.
        call createExplicitObjectListSingleSpace( ggd_grid,         &
                &   ggd_grid%grid_subset( GRID_SUBSET_X_POINTS),    &
                &   IDS_CLASS_NODE , xpoints(:, 1), IDS_CLASS_NODE, 1)
        ! call createExplicitObjectListSingleSpace( ggd_grid,         &
        !       &   ggd_grid%grid_subset( GRID_SUBSET_X_POINTS),  &
        !       &   GRID_SUBSET_X_POINTS, xpoints, CLASS_NODE(1:SPACE_COUNT), 1)
#if 0
      !> Set up specific grid subset by collection faces for regions

      !> Start counting from end of generic grid subset
      GSubsetCount = B2_GENERIC_GSUBSET_COUNT

      !> Cell + face grid subset
      do iType = REGIONTYPE_CELL, REGIONTYPE_YFACE
          
          select case(iType)
          case( REGIONTYPE_CELL )
              cls = CLASS_CELL
          case( REGIONTYPE_YFACE, REGIONTYPE_XFACE )
              cls = CLASS_POLOIDALRADIAL_FACE
          end select

          do iRegion = 1, regionCount(geoId, iType)              
              GSubsetCount = GSubsetCount + 1

              call logmsg( LOGDEBUG, "b2IMASFillGridDescription: add subgrid #"//idsInt2str(GSubsetCount)//&
                  & " for iType "//idsInt2str(iType)//&
                  &", iRegion "//idsInt2str(iRegion)//": "//regionName(geoId, iType, iRegion) )

              call createSubGridForExplicitList( ggd_grid, ggd_grid%grid_subset( GSubsetCount ), &
                  & cls(1:SPACE_COUNT), &
                  & collectIndexListForRegion(gmap, region, iType, iRegion), &
                  & regionName(geoId, iType, iRegion) )

          end do
      end do

      !> Add midplane node subgrids
      
      !> Find the core boundary subgrid by looking for its name as defined in b2mod_connectivity
      iCoreGS = gridFindSubGridByName(ggd_grid, "Core boundary")
      !> For double null, we need the outer half of the core boundary
      if (iCoreGS == GRID_UNDEFINED) then 
          iCoreGS = gridFindSubGridByName(ggd_grid, "Outer core boundary")          
      end if
      if (iCoreGS == GRID_UNDEFINED) stop "fillInGridSubsetDescription: &
          & did not find core boundary subgrid for assembling outer midplane subgrid"

      !> Figure out starting points for inner and outer midplane on core boundary
      call findMidplaneCells(ggd_grid%grid_subset(iCoreGS), gmap, crx, xIn, yIn, xOut, yOut)

      GSubsetCount = GSubsetCount + 1
      call createSubGridForExplicitList( ggd_grid, ggd_grid%grid_subset( GSubsetCount ), &
          & CLASS_NODE(1:SPACE_COUNT), &
          & collectRadialVertexIndexList(gmap, cflag, xIn, yIn, topix, topiy), &
          & "Inner midplane" )
      
      GSubsetCount = GSubsetCount + 1
      call createSubGridForExplicitList( ggd_grid, ggd_grid%grid_subset( GSubsetCount ), &
          & CLASS_NODE(1:SPACE_COUNT), &
          & collectRadialVertexIndexList(gmap, cflag, xOut, yOut, topix, topiy), &
          & "Outer midplane" )      

      call logmsg( LOGDEBUG, "b2IMASFillGridDescription: wrote total of "&
          &//idsInt2str(GSubsetCount)//" subgrids (expected was "//idsInt2str(size(ggd_grid%grid_subset))//')' )

      call assert( GSubsetCount == size(ggd_grid%subgrids) )
#endif
    end subroutine fillInGridSubsetDescription



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
