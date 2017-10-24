module b2mod_ual_io

    !> B2 modules

    use b2mod_types
    !use b2mod_version
    use b2mod_b2cmpa
    use b2mod_geo
    use b2mod_plasma
    use b2mod_constants
    !use b2mod_rates
    !use b2mod_residuals
    use b2mod_sources
    use b2mod_transport
    use b2mod_anomalous_transport
    !use b2mod_work
    !use b2mod_ppout
    use b2mod_indirect
    !use b2mod_neutrals_namelist
    !use b2mod_neutr_src_scaling

    !> B2/CPO Mapping
    use b2mod_ual_io_data
    use b2mod_ual_io_grid

    use logging

    !> UAL Access
#ifdef IMAS
    use ids_schemas  ! IGNORE
    use ids_routines ! IGNORE
    use ids_assert   ! IGNORE
    use ids_grid_common, IDS_COORDTYPE_R => COORDTYPE_R,    &   ! IGNORE
        &   IDS_COORDTYPE_Z => COORDTYPE_Z ! IGNORE
    use ids_string              ! IGNORE
    use ids_grid_subgrid        ! IGNORE
    use ids_grid_objectlist     ! IGNORE
    use ids_grid_examples       ! IGNORE
    use ids_grid_unstructured   ! IGNORE
    use ids_grid_structured     ! IGNORE
    use ids_grid_data           ! IGNORE
#else
#ifdef ITM
    use euITM_schemas  ! IGNORE
    use euITM_routines ! IGNORE
#endif
#endif
  ! use ggd , GGD_UNDEFINED => GRID_UNDEFINED

  implicit none

  logical, parameter, private :: INCLUDE_GHOST_CELLS = .false.

#ifdef IMAS

contains

    subroutine write_ids(edge_profiles,edge_sources,edge_transport)
#       include <git_version_B25.h>
        type (ids_edge_profiles)    :: edge_profiles
        type (ids_edge_sources)     :: edge_sources
        type (ids_edge_transport)   :: edge_transport

        !! Internal
        type(B2GridMap) :: gmap  !! Where is this type defined?
        type(ids_generic_grid_dynamic_grid_subset) ::  gs_cell, gs_face,   &
            &   gs_bnd_core
        integer :: is, ns, nx, ny, i, j
        logical, parameter :: B2_WRITE_DATA = .true.
        real(IDS_real),   &
            &   dimension( -1:ubound( crx, 1 ), -1:ubound( crx, 2), 3, 3) :: e
        integer :: iGsCore, iGsInnerMidplane, iGsOuterMidplane

        real(IDS_real) :: tmpFace( -1:ubound( na, 1), -1:ubound( na, 2 ), 0:1)
        real(IDS_real) :: tmpVx( -1:ubound( na, 1), -1:ubound( na, 2 ))
        integer ::  homogeneous_time
        real(IDS_real)    ::  time
        integer :: ggd_slice, num_ggd_slice
        integer :: ne_slice, num_ne_slices
        integer :: eflux_slice, num_eflux_slices

        !> ===  SET UP IDS ===
        write(0,*) "Setting data for edge_profiles IDS"

        !> Preparing database for writing
        !> Through practice it was disclosed that there are some mandatory 
        !> steps to be done in order to assure for data to be successfully 
        !> written to IDS. Without going through those steps errors and failed 
        !> process of writing to IDS are to be expected.
        !> This can be done using exampleSetIDSFundamentals routine
        homogeneous_time = 1
        time = 0.0_IDS_real
        call exampleSetIDSFundamentals( edge_profiles, homogeneous_time, time) 

        !> Allocate ggd slice
        num_ggd_slice = 1
        ggd_slice = 1
        allocate( edge_profiles%ggd( num_ggd_slice ) )

        !> allocate and init the IDS
        allocate( edge_profiles%code%name(1) )
# ifdef B25_EIRENE
        edge_profiles%code%name = "SOLPS-ITER"
# else
        edge_profiles%code%name = "B2.5"
# endif
        allocate( edge_profiles%code%version(1) )
        edge_profiles%code%version = git_version_B25

        ns = size( na, 3 ) 
        nx = ubound( na, 1 )
        ny = ubound( na, 2 )

        !> List of species
        allocate( edge_profiles%ggd( ggd_slice )%ion( ns ) )
        do is = 0, ns-1
            allocate( edge_profiles%ggd( ggd_slice )%ion( is + 1 )%state(1) )
            allocate( edge_profiles%ggd( ggd_slice )%ion( is + 1 )%element(1) )

            call species( is, edge_profiles%ggd( ggd_slice )%ion( is + 1 )% &
                &   state(1)%label, .false.)
            edge_profiles%ggd( ggd_slice )%ion( is + 1 )%element(1)%a = am( is )
            edge_profiles%ggd( ggd_slice )%ion( is + 1 )%element(1)%z_n =   &
                &   zn( is )
            edge_profiles%ggd( ggd_slice )%ion( is + 1 )%state(1)%z_min =   &
                &   zamin( is )
            edge_profiles%ggd( ggd_slice )%ion( is + 1 )%state(1)%z_max =   &
                &   zamax( is )
 
        enddo

        write(*,*) "Running b2CreateMap subroutine"
        !> Set up the B2<->IDS mappings
        call b2CreateMap( nx, ny, crx( -1:nx, -1:ny, : ),             &
            &   cry( -1:nx, -1:ny, : ), cflags, leftix, leftiy,       &
            &   rightix, rightiy, topix, topiy, bottomix,bottomiy,    & 
            &   INCLUDE_GHOST_CELLS, gmap )
        mapInitialized = .true.

        !> Write grid & grid subsets/subgrids
        call b2IMASFillGridDescription( gmap, edge_profiles%ggd( ggd_slice )%grid,    &
            &   nx, ny, crx(-1:nx, -1:ny, :), cry(-1:nx, -1:ny, : ),        &
            &   leftix, leftiy, rightix, rightiy, topix, topiy, bottomix,   &
            &   bottomiy, nnreg, topcut, region, cflags,                    &
            &   INCLUDE_GHOST_CELLS, vol, gs, qc )

! #if 0

    call assert( geometryId( nnreg, periodic_bc, topcut ) ==    &
        &   GEOMETRY_SN, "write_ids: can only do single null" )

    !> Write plasma state 

    if ( B2_WRITE_DATA ) then

        write (*,*) "b2mod_ual_io.write_ids: writing plasma state"
        
        iGsCore = findGridSubsetByName( edge_profiles%ggd( ggd_slice )%grid,    &
            &   "Core boundary" )
        iGsInnerMidplane = findGridSubsetByName( edge_profiles% &
            &   ggd( ggd_slice )%grid, "Inner midplane" )
        iGsOuterMidplane = findGridSubsetByName( edge_profiles% &
            &   ggd( ggd_slice )%grid, "Outer midplane" )

        !> TODO: The fluxes are currently in the edge_transport. They are 
        !> supposed to be in the edge_profiles (24.10.2017)
        !> use electrons%particles or electrons%energy node? Currently using 
        !> electrons%energy node

        !> edge_transport fundamentals
        !> TODO: create it as a subroutine
        allocate( edge_transport%model(1) )
        allocate( edge_transport%model(1)%ggd( num_ggd_slice ) )

        !> ne
        num_ne_slices = 1
        ne_slice = 1
        num_eflux_slices = 1
        eflux_slice = 0
        allocate( edge_profiles%ggd( ggd_slice )%electrons% &
            &   density( num_ne_slices ) )
        call write_quantity( edge_profiles%ggd( ggd_slice )%electrons%  &
            &   density, edge_transport%model(1)%    &
            &   ggd( ggd_slice )%electrons%energy%flux, ne, fne, ggd_slice )
        ! call write_quantity( edge_profiles%ggd( ggd_slice )%electrons%  &
        !     &   density( ne_slice )%values, edge_transport%model(1)%    &
        !     &   ggd( ggd_slice )%electrons%energy%flux, ne, fne, ggd_slice )
#if 0
        call write_cell_scalar( edge_profiles%fluid%ne%source,  &
            &   sne(:,:,0) + sne(:,:,1) * ne )

        !> na
        allocate( edge_profiles%fluid%ni( ns ) )
        do is = 1, ns
            call write_quantity( edge_profiles%fluid%ni(is)%value,          &
                &   edge_profiles%fluid%ni( is )%flux, na(:,:, is - 1 ),    &
                &   fna(:,:,:, is - 1 ) )
            call write_cell_scalar( edge_profiles%fluid%ni( is )%source,    &
                &   sna(:,:,0, is - 1 ) + sna(:,:,1, is - 1 ) * na(:,:, is - 1 ) )
        end do

!!$    !> ue TODO: must be computed, refactor code from b2news into function
!!$    allocate(edge_profiles%fluid%ve%comps(1))
!!$    allocate(edge_profiles%fluid%ve%align(1))
!!$    allocate(edge_profiles%fluid%ve%alignid(1))
!!$    edge_profiles%fluid%ve%align(1) = VEC_ALIGN_PARALLEL
!!$    edge_profiles%fluid%ve%alignid(1) = VEC_ALIGN_PARALLEL_ID
!!$    call write_cell_scalar( edge_profiles%fluid%ve%comps(1)%value, ue(:,:) )

        !> ua
        allocate( edge_profiles%fluid%vi( ns ) )
        do is = 1, ns
            allocate( edge_profiles%fluid%vi( is )%comps(1) )
            allocate( edge_profiles%fluid%vi( is )%align(1) )
            allocate( edge_profiles%fluid%vi( is )%alignid(1) )
            edge_profiles%fluid%vi( is )%align(1) = VEC_ALIGN_PARALLEL
            edge_profiles%fluid%vi( is )%alignid(1) = VEC_ALIGN_PARALLEL_ID

            call write_cell_scalar( edge_profiles%fluid%vi( is )%comps(1)%value,&
                &   ua(:,:, is - 1 ) )
        end do

        !> te
        call write_quantity( edge_profiles%fluid%te%value,  &
            &   edge_profiles%fluid%te%flux, te/qe, fhe )

        !> ti
        allocate( edge_profiles%fluid%ti(1) )
        call write_quantity( edge_profiles%fluid%ti(1)%value,   &
            &   edge_profiles%fluid%ti(1)%flux, ti/qe, fhi )

        !> po
        call write_cell_scalar( edge_profiles%fluid%po%value, po )

        !> B (magnetic field vector)        
        allocate( edge_profiles%fluid%te_aniso%comps(4) )

        !> Compute unit basis vectors along the field directions
        call computeCoordinateUnitVectors(crx, cry, e(:,:,:,1), e(:,:,:,2),     &
            &   e(:,:,:,3))

        !> Write the three unit basis vectors
        do i = 1, 3
            allocate( edge_profiles%fluid%te_aniso%comps(i)%flux(1) )
            call write_cell_vector( &
                &   edge_profiles%fluid%te_aniso%comps(i)%flux(1), &
                & (/ VEC_ALIGN_DEFAULT, VEC_ALIGN_DEFAULT, VEC_ALIGN_DEFAULT /), &
                & (/ VEC_ALIGN_DEFAULT_ID, VEC_ALIGN_DEFAULT_ID, VEC_ALIGN_DEFAULT_ID /), &
                & e(:,:,:,i) )
        end do

        !> write the magnetic field vector in the b2 coordinate system
        allocate( edge_profiles%fluid%te_aniso%comps(4)%flux(1) )
        call write_cell_vector( edge_profiles%fluid%te_aniso%comps(4)%flux(1),  &
            & (/ VEC_ALIGN_POLOIDAL, VEC_ALIGN_RADIAL, VEC_ALIGN_TOROIDAL /),   &
            & (/ VEC_ALIGN_POLOIDAL_ID, VEC_ALIGN_RADIAL_ID,                    &
            &   VEC_ALIGN_TOROIDAL_ID /), bb(:,:,0:2) )
#endif
      
    end if

    call logmsg( LOGDEBUG, "b2mod_ual_io.write_cpo: done" )

    contains

        !> Write a scalar B2 cell quantity to a ids_generic_grid_scalar
        subroutine write_quantity( val, fluxes, value, flux, ggd_slice )
            use b2mod_interp    ! IGNORE
            type(ids_generic_grid_scalar), pointer, intent(inout) :: val(:)
            ! real(IDS_real), pointer, intent(inout)  :: val(:)
            type (ids_generic_grid_scalar), pointer, intent(inout) :: fluxes(:)
            real(IDS_real), intent(in) :: value( -1:gmap%b2nx, -1:gmap%b2ny )
            real(IDS_real), intent(in) :: flux( -1:gmap%b2nx, -1:gmap%b2ny, 0:1 )  
            integer, intent(in) ::  ggd_slice                   
            real(IDS_real), dimension(:), pointer :: idsdata

            allocate( val(5) )
#if 0
            idsdata => b2IMASTransformDataB2ToIDS( edge_profiles%ggd( ggd_slice )%grid,  &
                &   GRID_SUBSET_CELLS, gmap, value )
            call gridWriteDataScalar( val(1), GRID_SUBSET_CELLS, idsdata ) 
            deallocate( idsdata )
            tmpFace = 0.0_IDS_real
            call value_on_faces( nx, ny, vol, value, tmpFace)
            idsdata => b2IMASTransformDataB2ToIDS(  &
                &   edge_profiles%ggd( ggd_slice )%grid, iGsCore, gmap, tmpFace )
            call gridWriteData( val(2), iGsCore, idsdata )
            deallocate( idsdata )
            tmpVx = interpolateToVertices(  &
                &   gmap%b2nx, gmap%b2ny, VX_LOWERLEFT, value )
            idsdata => b2IMASTransformDataB2ToIDSVertex(    &
                &   edge_profiles%ggd( ggd_slice )%grid, iGsInnerMidplane, gmap, tmpVx  )
            call gridWriteData( val(3), iGsInnerMidplane, idsdata )
            deallocate( idsdata )
            idsdata => b2IMASTransformDataB2ToIDSVertex(    &
                &   edge_profiles%ggd( ggd_slice )%grid, iGsOuterMidplane, gmap, tmpVx )
            call gridWriteData( val(4), iGsOuterMidplane, idsdata )
            deallocate( idsdata )
            idsdata => b2IMASTransformDataB2ToIDSVertex(    &
                &   edge_profiles%ggd( ggd_slice )%grid, GRID_SUBSET_NODES, gmap, tmpVx )
            call gridWriteData( val(5), GRID_SUBSET_NODES, idsdata )
            deallocate( idsdata )
            allocate( fluxes(2) )
! #if 0
            call write_face_vector( fluxes(1), flux, ggd_slice)
            call write_face_vector( fluxes(2), flux, ggd_slice, gridSubsetId = iGsCore )
#endif 
        end subroutine write_quantity

#if 0
        ! Write a vector B2 face quantity to a complexgrid_vector
        subroutine write_face_vector(vector, b2FaceData, ggd_slice, gridSubsetId)
            type(ids_generic_grid_scalar), intent(inout) :: vector
            real(IDS_real), intent(in) :: b2FaceData(-1:gmap%b2nx, -1:gmap%b2ny, 0:1)
            integer, intent(in), optional :: gridSubsetId
            real(IDS_real), dimension(:), pointer :: idsdata
            integer, intent(in) :: ggd_slice

            ! if ( .not. present(gridSubsetId) ) then
            !     !> TODO: add checks whether already allocated
            !     allocate(vector%comp(2))
            !     allocate(vector%align(2))
            !     allocate(vector%alignid(2))
              
            !     !> Fill in alignment information for vector components
            !     vector%align(1) = VEC_ALIGN_POLOIDAL
            !     vector%alignid(1) = VEC_ALIGN_POLOIDAL_ID
              
            !     vector%align(2) = VEC_ALIGN_RADIAL
            !     vector%alignid(2) = VEC_ALIGN_RADIAL_ID          

            !     !> Fill in vector component data
            !     idsdata => b2IMASTransformDataB2ToIDS(edgecpo%ggd( ggd_slice )%grid, B2_SUBGRID_FACES_Y, gmap, b2FaceData)
            !     call gridWriteData( vector%comp(1), B2_SUBGRID_FACES_Y, idsdata ) 
            !     deallocate(idsdata)
            !     idsdata => b2IMASTransformDataB2ToIDS(edgecpo%ggd( ggd_slice )%grid, B2_SUBGRID_FACES_X, gmap, b2FaceData)
            !     call gridWriteData( vector%comp(2), B2_SUBGRID_FACES_X, idsdata ) 
            !     deallocate(idsdata)
            ! else
            !     allocate(vector%comp(1))
            !     allocate(vector%align(1))
            !     allocate(vector%alignid(1))

            !     vector%align(1) = GRID_UNDEFINED
            !     vector%alignid(1) = ""

            !     idsdata => b2IMASTransformDataB2ToIDS(edgecpo%ggd( ggd_slice )%grid, gridSubsetId, gmap, b2FaceData)
            !     call gridWriteData( vector%comp(1), gridSubsetId, idsdata ) 
            !     deallocate(idsdata)
            ! end if

        end subroutine write_face_vector
#endif
        ! return
    end subroutine write_ids

#else
# ifdef ITM

contains

  subroutine write_cpo(edgecpo)
    type (type_edge) :: edgecpo

    ! internal
    type(B2ITMGridMap) :: gmap
    type(type_complexgrid_subgrid) :: sg_cell, sg_face, sg_bnd_core
    integer :: is, ns, nx, ny, i
    logical, parameter :: B2_WRITE_DATA = .true.
    real(ITM_R8), dimension(-1:ubound(crx,1),-1:ubound(crx,2),3,3) :: e
    integer :: iSgCore, iSgInnerMidplane, iSgOuterMidplane

    real(ITM_R8) :: tmpFace(-1:ubound(na, 1), -1:ubound(na, 2), 0:1)
    real(ITM_R8) :: tmpVx(-1:ubound(na, 1), -1:ubound(na, 2))


    ! allocate and init the cpo
    allocate(edgecpo%datainfo%dataprovider(1))
    edgecpo%datainfo%dataprovider="IPP"
    allocate(edgecpo%codeparam%codename(1))
    edgecpo%codeparam%codename(1)="B2.5"
    edgecpo%time= 0.0D0


    ns = size(na, 3) 
    nx = ubound(na, 1)
    ny = ubound(na, 2)


! species block
    allocate(edgecpo%species(ns))
    do is = 0, ns-1
       allocate(edgecpo%species(is+1)%label(1))
       call species(is, edgecpo%species(is+1)%label, .false.)
       edgecpo%species(is+1)%amn = am(is)
       edgecpo%species(is+1)%zn = zn(is)
       edgecpo%species(is+1)%zmin = zamin(is)
       edgecpo%species(is+1)%zmax = zamax(is)
    enddo


    ! set up the B2<->CPO mappings
    call b2ITMCreateMap( nx,ny,crx(-1:nx,-1:ny,: ),cry(-1:nx,-1:ny,:),&
        & cflags,leftix,leftiy,rightix,rightiy, &
        & topix,topiy,bottomix,bottomiy, INCLUDE_GHOST_CELLS, gmap )
    mapInitialized = .true.

    ! write grid & subgrids
    call b2ITMFillGridDescription( gmap, edgecpo%grid, &
        & nx,ny,crx(-1:nx,-1:ny,:),cry(-1:nx,-1:ny,:), &
        & leftix,leftiy,rightix,rightiy, &
        & topix,topiy,bottomix,bottomiy, &
        & nnreg, topcut, region, cflags, INCLUDE_GHOST_CELLS, vol, gs, qc )


    call assert( geometryId( nnreg, periodic_bc, topcut ) == GEOMETRY_SN, "write_cpo: can only do single null" )

    ! Write plasma state 

    if ( B2_WRITE_DATA ) then

        write (*,*) "b2mod_ual_io.write_cpo: writing plasma state"
        
        iSgCore = gridFindSubGridByName( edgecpo%grid, "Core boundary" )
        iSgInnerMidplane = gridFindSubGridByName( edgecpo%grid, "Inner midplane" )
        iSgOuterMidplane = gridFindSubGridByName( edgecpo%grid, "Outer midplane" )

        ! ne
        call write_quantity( edgecpo%fluid%ne%value, edgecpo%fluid%ne%flux, ne, fne )
        call write_cell_scalar( edgecpo%fluid%ne%source, sne(:,:,0) + sne(:,:,1)*ne )

        ! na
        allocate(edgecpo%fluid%ni(ns))
        do is = 1, ns
            call write_quantity( edgecpo%fluid%ni(is)%value, edgecpo%fluid%ni(is)%flux, na(:,:,is-1), fna(:,:,:,is-1) )
            call write_cell_scalar( edgecpo%fluid%ni(is)%source, sna(:,:,0,is-1) + sna(:,:,1,is-1)*na(:,:,is-1) )
        end do

!!$    ! ue TODO: must be computed, refactor code from b2news into function
!!$    allocate(edgecpo%fluid%ve%comps(1))
!!$    allocate(edgecpo%fluid%ve%align(1))
!!$    allocate(edgecpo%fluid%ve%alignid(1))
!!$    edgecpo%fluid%ve%align(1) = VEC_ALIGN_PARALLEL
!!$    edgecpo%fluid%ve%alignid(1) = VEC_ALIGN_PARALLEL_ID
!!$    call write_cell_scalar( edgecpo%fluid%ve%comps(1)%value, ue(:,:) )

        ! ua
        allocate(edgecpo%fluid%vi(ns))
        do is = 1, ns
            allocate(edgecpo%fluid%vi(is)%comps(1))
            allocate(edgecpo%fluid%vi(is)%align(1))
            allocate(edgecpo%fluid%vi(is)%alignid(1))
            edgecpo%fluid%vi(is)%align(1) = VEC_ALIGN_PARALLEL
            edgecpo%fluid%vi(is)%alignid(1) = VEC_ALIGN_PARALLEL_ID

            call write_cell_scalar( edgecpo%fluid%vi(is)%comps(1)%value, ua(:,:,is-1) )
        end do

        ! te
        call write_quantity( edgecpo%fluid%te%value, edgecpo%fluid%te%flux, te/qe, fhe )

        ! ti
        allocate(edgecpo%fluid%ti(1))
        call write_quantity( edgecpo%fluid%ti(1)%value, edgecpo%fluid%ti(1)%flux, ti/qe, fhi )


        ! po
        call write_cell_scalar( edgecpo%fluid%po%value, po )

        ! B (magnetic field vector)        
        allocate(edgecpo%fluid%te_aniso%comps(4))

        ! Compute unit basis vectors along the field directions
        call computeCoordinateUnitVectors(crx, cry, e(:,:,:,1), e(:,:,:,2), e(:,:,:,3))

        ! Write the three unit basis vectors
        do i = 1, 3
            allocate(edgecpo%fluid%te_aniso%comps(i)%flux(1))
            call write_cell_vector( edgecpo%fluid%te_aniso%comps(i)%flux(1), &
                & (/ VEC_ALIGN_DEFAULT, VEC_ALIGN_DEFAULT, VEC_ALIGN_DEFAULT /), &
                & (/ VEC_ALIGN_DEFAULT_ID, VEC_ALIGN_DEFAULT_ID, VEC_ALIGN_DEFAULT_ID /), &
                & e(:,:,:,i) )
        end do

        ! write the magnetic field vector in the b2 coordinate system
        allocate(edgecpo%fluid%te_aniso%comps(4)%flux(1))
        call write_cell_vector( edgecpo%fluid%te_aniso%comps(4)%flux(1), &
            & (/ VEC_ALIGN_POLOIDAL, VEC_ALIGN_RADIAL, VEC_ALIGN_TOROIDAL /), &
            & (/ VEC_ALIGN_POLOIDAL_ID, VEC_ALIGN_RADIAL_ID, VEC_ALIGN_TOROIDAL_ID /), &
            &  bb(:,:,0:2) )
      
    end if

    call logmsg( LOGDEBUG, "b2mod_ual_io.write_cpo: done" )

  contains

    ! Write a scalar B2 cell quantity to a complexgrid_scalar
    subroutine write_quantity( values, fluxes, value, flux )
      use b2mod_interp
      type(type_complexgrid_scalar), pointer, intent(inout) :: values(:)
      type(type_complexgrid_vector), pointer, intent(inout) :: fluxes(:)
      real(ITM_R8), intent(in) :: value(-1:gmap%b2nx, -1:gmap%b2ny)
      real(ITM_R8), intent(in) :: flux(-1:gmap%b2nx, -1:gmap%b2ny, 0:1)                     
      real(ITM_R8), dimension(:), pointer :: cpodata

      allocate(values(5))
      cpodata => b2ITMTransformDataB2ToCpo( edgecpo%grid, B2_SUBGRID_CELLS, gmap, value )
      call gridWriteData( values(1), B2_SUBGRID_CELLS, cpodata ) 
      deallocate(cpodata)
      tmpFace = 0.0_ITM_R8
      call value_on_faces(nx,ny,vol,value,tmpFace)
      cpodata => b2ITMTransformDataB2ToCpo( edgecpo%grid, iSgCore, gmap, tmpFace )
      call gridWriteData( values(2), iSgCore, cpodata )
      deallocate(cpodata)
      tmpVx = interpolateToVertices( gmap%b2nx, gmap%b2ny, VX_LOWERLEFT, value )
      cpodata => b2ITMTransformDataB2ToCpoVertex( edgecpo%grid, iSgInnerMidplane, gmap, tmpVx  )
      call gridWriteData( values(3), iSgInnerMidplane, cpodata )
      deallocate(cpodata)
      cpodata => b2ITMTransformDataB2ToCpoVertex( edgecpo%grid, iSgOuterMidplane, gmap, tmpVx )
      call gridWriteData( values(4), iSgOuterMidplane, cpodata )
      deallocate(cpodata)
      cpodata => b2ITMTransformDataB2ToCpoVertex( edgecpo%grid, B2_SUBGRID_NODES, gmap, tmpVx )
      call gridWriteData( values(5), B2_SUBGRID_NODES, cpodata )
      deallocate(cpodata)
      allocate( fluxes(2) )
      call write_face_vector( fluxes(1), flux )
      call write_face_vector( fluxes(2), flux, subgridInd = iSgCore )
    end subroutine write_quantity


    ! Write a scalar B2 cell quantity to a complexgrid_scalar
    subroutine write_cell_scalar(scalar, b2CellData)
      type(type_complexgrid_scalar), intent(inout), pointer :: scalar(:)
      real(ITM_R8), intent(in) :: b2CellData(-1:gmap%b2nx, -1:gmap%b2ny)
      real(ITM_R8), dimension(:), pointer :: cpodata

      ! TODO: add checks whether already allocated
      allocate(scalar(1))
      cpodata => b2ITMTransformDataB2ToCpo( edgecpo%grid, B2_SUBGRID_CELLS, gmap, b2CellData )
      call gridWriteData( scalar(1), B2_SUBGRID_CELLS, cpodata ) 
      deallocate(cpodata)
    end subroutine write_cell_scalar


    ! Write a vector B2 cell quantity to a complexgrid_vector
    subroutine write_cell_vector(vector, align, alignid, vecdata)
      type(type_complexgrid_vector), intent(inout) :: vector
      real(ITM_R8), intent(in) :: vecdata(-1:gmap%b2nx, -1:gmap%b2ny, 0:2)
      integer, intent(in) :: align(3)
      character(LEN=132), intent(in) :: alignid(3)
      real(ITM_R8), dimension(:), pointer :: cpodata

      ! internal
      integer :: dim, i

      dim = size(vecdata, 3)

      ! TODO: add checks whether already allocated
      allocate(vector%comp(dim))
      allocate(vector%align(dim))
      allocate(vector%alignid(dim))

      ! Fill in alignment information for vector components
      vector%align = align
      vector%alignid = alignid

      ! Fill in vector component data
      do i = 1, dim
         cpodata => b2ITMTransformDataB2ToCpo(edgecpo%grid, B2_SUBGRID_CELLS, gmap, vecdata(:,:,i-1))
         call gridWriteData( vector%comp(i), B2_SUBGRID_CELLS, cpodata )
         deallocate(cpodata)
      end do

    end subroutine write_cell_vector


    ! Write a vector B2 face quantity to a complexgrid_vector
    subroutine write_face_vector(vector, b2FaceData, subgridInd)
      type(type_complexgrid_vector), intent(inout) :: vector
      real(ITM_R8), intent(in) :: b2FaceData(-1:gmap%b2nx, -1:gmap%b2ny, 0:1)
      integer, intent(in), optional :: subgridInd
      real(ITM_R8), dimension(:), pointer :: cpodata

!!$      if ( .not. present(subgridInd) ) then
!!$          ! TODO: add checks whether already allocated
!!$          allocate(vector%comp(2))
!!$          allocate(vector%align(2))
!!$          allocate(vector%alignid(2))
!!$          
!!$          ! Fill in alignment information for vector components
!!$          vector%align(1) = VEC_ALIGN_POLOIDAL
!!$          vector%alignid(1) = VEC_ALIGN_POLOIDAL_ID
!!$          
!!$          vector%align(2) = VEC_ALIGN_RADIAL
!!$          vector%alignid(2) = VEC_ALIGN_RADIAL_ID          
!!$
!!$          ! Fill in vector component data
!!$          cpodata => b2ITMTransformDataB2ToCpo(edgecpo%grid, B2_SUBGRID_FACES_Y, gmap, b2FaceData)
!!$          call gridWriteData( vector%comp(1), B2_SUBGRID_FACES_Y, cpodata ) 
!!$          deallocate(cpodata)
!!$          cpodata => b2ITMTransformDataB2ToCpo(edgecpo%grid, B2_SUBGRID_FACES_X, gmap, b2FaceData)
!!$          call gridWriteData( vector%comp(2), B2_SUBGRID_FACES_X, cpodata ) 
!!$          deallocate(cpodata)
!!$      else
!!$          allocate(vector%comp(1))
!!$          allocate(vector%align(1))
!!$          allocate(vector%alignid(1))
!!$
!!$          vector%align(1) = GRID_UNDEFINED
!!$          vector%alignid(1) = ""
!!$
!!$          cpodata => b2ITMTransformDataB2ToCpo(edgecpo%grid, subgridInd, gmap, b2FaceData)
!!$          call gridWriteData( vector%comp(1), subgridInd, cpodata ) 
!!$          deallocate(cpodata)
!!$      end if

    end subroutine write_face_vector

  end subroutine write_cpo

  !> From the B2 grid, compute the coordinate unit vectors (poloidal, radial. toroidal)
  subroutine computeCoordinateUnitVectors(crx, cry, e1, e2, e3)
    real(ITM_R8), intent(in), dimension(-1:,-1:,0:) :: crx, cry
    real(ITM_R8), intent(out), dimension(-1:ubound(crx,1),-1:ubound(crx,2),3) :: e1, e2, e3
    
    ! internal
    integer :: ix, iy, ixn, iyn, nx, ny
    real(ITM_R8), dimension(0:1) :: cC, cN
    real(ITM_R8) :: dir

    e1 = 0.0
    e2 = 0.0
    e3 = 0.0

    ! poloidal vectors

    nx = ubound(crx,1)
    ny = ubound(crx,2)
    do ix = -1, nx
        do iy = -1, ny

            cC = quadCentroid( &
                & crx(ix, iy, 0), cry(ix, iy, 0), &
                & crx(ix, iy, 1), cry(ix, iy, 1), &
                & crx(ix, iy, 2), cry(ix, iy, 2), &
                & crx(ix, iy, 3), cry(ix, iy, 3) )

            ! poloidal direction
            ! Try to find right neighbour
            dir = 1.0
            ixn = rightix( ix, iy )
            iyn = rightiy( ix, iy )

            if ( .not. isInDomain( nx, ny, ixn, iyn ) ) then
                ! If not found, try to find left neighbour
                ! ...and note to invert vector direction
                dir = -1.0
                ixn = leftix( ix, iy )
                iyn = leftiy( ix, iy )
            end if
            if ( .not. isInDomain( nx, ny, ixn, iyn ) ) then
                !stop "computeCoordinateUnitVectors: not able to find poloidal neighbour for cell"
                ! skip cell
                cycle
            end if
            
            cN = quadCentroid( &
                & crx(ixn, iyn, 0), cry(ixn, iyn, 0), &
                & crx(ixn, iyn, 1), cry(ixn, iyn, 1), &
                & crx(ixn, iyn, 2), cry(ixn, iyn, 2), &
                & crx(ixn, iyn, 3), cry(ixn, iyn, 3) )

            ! compute vector from one centroid to the other
            e1(ix,iy,1) = cN(0) - cC(0)   ! R
            e1(ix,iy,2) = 0.0             ! phi
            e1(ix,iy,3) = cN(1) - cC(1)   ! Z
            
            e1(ix,iy,:) = e1(ix,iy,:) * dir  ! fix direction


            ! radial direction
            ! Try to find top neighbour
            dir = 1.0
            ixn = topix( ix, iy )
            iyn = topiy( ix, iy )

            if ( .not. isInDomain( nx, ny, ixn, iyn ) ) then
                ! If not found, try to find bottom neighbour
                ! ...and note to invert vector direction
                dir = -1.0
                ixn = bottomix( ix, iy )
                iyn = bottomiy( ix, iy )
            end if
            if ( .not. isInDomain( nx, ny, ixn, iyn ) ) then
                !stop "computeCoordinateUnitVectors: not able to find toroidal neighbour for cell"
                ! skip cell
                cycle
            end if
            
            cN = quadCentroid( &
                & crx(ixn, iyn, 0), cry(ixn, iyn, 0), &
                & crx(ixn, iyn, 1), cry(ixn, iyn, 1), &
                & crx(ixn, iyn, 2), cry(ixn, iyn, 2), &
                & crx(ixn, iyn, 3), cry(ixn, iyn, 3) )

            ! compute vector from one centroid to the other
            e2(ix,iy,1) = cN(0) - cC(0)   ! R
            e2(ix,iy,2) = 0.0             ! phi
            e2(ix,iy,3) = cN(1) - cC(1)   ! Z
            
            e2(ix,iy,:) = e2(ix,iy,:) * dir  ! fix direction


            ! toroidal direction
            e3(ix,iy,1) = 0.0   ! R
            e3(ix,iy,2) = 1.0   ! phi
            e3(ix,iy,3) = 0.0   ! Z
            

            ! make unit vectors
            e1(ix,iy,:) = unitVector(e1(ix,iy,:))
            e2(ix,iy,:) = unitVector(e2(ix,iy,:))
            e3(ix,iy,:) = unitVector(e3(ix,iy,:))

        end do
    end do

  end subroutine computeCoordinateUnitVectors


  !> Return unit vector along direction of given vector
  function unitVector(v) result(unitV)
    real(ITM_R8), intent(in) :: v(:)
    real(ITM_R8) :: unitV(size(v))

    unitV = v / sqrt( sum( v**2 ) )
  end function unitVector

# endif
#endif
 
end module b2mod_ual_io

!!!Local Variables:
!!! mode: f90
!!! End:
