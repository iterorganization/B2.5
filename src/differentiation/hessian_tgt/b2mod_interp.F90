!!-----------------------------------------------------------------------------
!! DOCUMENTATION:
!>      @section b2mod_interp_desc Description
!!      Module providing routines for interpolation of a cell-centred quantity to
!!      cell faces, computing flow velocity quantities, providing values for
!!      cell-centered quantities etc.
!!
!!      @subsection b2mod_gmap_pv  Parameters/variables
!!      @param  BB_INVALID
!!
!!-----------------------------------------------------------------------------
module b2mod_interp

    use b2mod_types
    use b2mod_connectivity
    use carre_constants

    implicit none

    real(R8), parameter :: BB_INVALID = 1e5_R8

contains

    !> Routine for computing magnetic field at faces
    subroutine compute_b_at_faces( nx, ny, bb, wbbl, wbbr, wbbv,    &
        &   leftix, leftiy, rightix, rightiy, bottomix, bottomiy,   &
        &   topix, topiy, vol, gs, qc, cflags )

        implicit none
        integer, intent(in)  :: nx  !< Specifies the number of interior cells
                                    !< along the first coordinate
        integer, intent(in)  :: ny  !< Specifies the number of interior cells
                                    !< along the second coordinate
        real(R8), intent(in) :: bb( -1:nx, -1:ny, 0:3 ) !< Magnetic field
        real(R8), intent(in) :: vol( -1:nx ,-1:ny ,0:4 )    !< Cell volume
        real(R8), intent(in) :: gs( -1:nx, -1:ny, 0:2 )
        real(R8), intent(in) :: qc( -1:nx, -1:ny , 0:1) !< Cosine of the angle
            !< between flux line direction and left cell face
        real(R8), intent(out) :: wbbl( -1:nx, -1:ny, 0:3 )
        real(R8), intent(out) :: wbbr( -1:nx, -1:ny, 0:3 )
        real(R8), intent(out) :: wbbv( -1:nx, -1:ny, 0:3 )
        integer, intent(in) :: cflags( -1:nx, -1:ny, CARREOUT_NCELLFLAGS )
        integer, intent(in) :: leftix( -1:nx, -1:ny)    !< Left neighbour
            !< poloidal (first coordinate) index array
        integer, intent(in) :: leftiy( -1:nx, -1:ny )   !< Left neighbour radial
            !< (second coordinate) index
        integer, intent(in) :: rightix( -1:nx, -1:ny )  !< Right neighbour
            !< poloidal (first coordinate) index array
        integer, intent(in) :: rightiy( -1:nx, -1:ny )  !< Right neighbour
            !< radial (second coordinate) index
        integer, intent(in) :: topix( -1:nx, -1:ny )    !< Top neighbour
            !< poloidal (first coordinate) index array
        integer, intent(in) :: topiy( -1:nx, -1:ny )    !< Top neighbour radial
            !< (second coordinate) index
        integer, intent(in) :: bottomix( -1:nx, -1:ny ) !< Bottom neighbour
            !< poloidal (first coordinate) index array
        integer, intent(in) :: bottomiy( -1:nx, -1:ny ) !< Bottom neighbour
            !< radial (second coordinate) index

        !! Internal variables
        integer :: ix   !< x-aligned (poloidal) cell index
        integer :: iy   !< y-aligned (radial) cell index

        !! Interpolate magnetic field to cell faces
        call interp_magnetic_field( nx, ny, bb, wbbl, wbbr, wbbv,       &
            &   leftix, leftiy, rightix, rightiy, bottomix, bottomiy,   &
            &   topix, topiy, vol, gs, qc, cflags )

        !! Extrapolate magnetic field to edges
        do iy = -1, ny
            do ix = -1, nx
                if ( isUnusedCell( cflags( ix, iy, CELLFLAG_TYPE ))) cycle
                !! Left face
                if ( .not. isInDomain(  &
                    &   nx, ny, leftix( ix, iy ), leftiy( ix, iy ))) then
                    if ( .not. isInDomain(  &
                        &   nx, ny, rightix( ix, iy ), rightiy( ix, iy ))) then
                        !! Not really an extrapolation, have nothing to
                        !! work with.
                        wbbl( ix, iy, 0:2 ) = bb( ix, iy, 0:2 )
                    else
                        wbbl( ix, iy, 0:2 ) = bb( ix, iy, 0:2 ) -   &
                            &   ( wbbr( ix, iy, 0:2 )-bb( ix, iy, 0:2 ))
                    endif
                endif
                !! Right face
                if ( .not. isInDomain(  &
                    &   nx, ny, rightix( ix, iy ), rightiy( ix, iy ))) then
                    if ( .not. isInDomain(  &
                        &   nx, ny, leftix( ix, iy ), leftiy( ix, iy ))) then
                        wbbr( ix, iy, 0:2 ) = bb( ix, iy, 0:2 )
                    else
                        wbbr( ix, iy, 0:2 ) = bb( ix, iy, 0:2 ) -   &
                            &   ( wbbl( ix, iy, 0:2 ) - bb( ix, iy, 0:2 ))
                    endif
                endif
                !! Bottom face
                if ( .not. isInDomain(  &
                    &   nx, ny, bottomix( ix, iy), bottomiy( ix, iy ))) then
                    if ( .not. isInDomain(  &
                        &   nx, ny, topix( ix, iy), topiy( ix, iy ))) then
                        wbbv( ix, iy, 0:2 ) = bb( ix, iy, 0:2 )
                    else
                        wbbv( ix, iy, 0:2 ) = bb( ix, iy, 0:2 ) &
                            & -1- ( wbbv( topix( ix, iy ),      &
                            &   topiy( ix, iy ), 0:2 ) - bb( ix, iy, 0:2 ))
                    endif
                endif
                !! TODO: might have to introduce wbbu (upper face) here.
            enddo
        enddo

        !! Compute extrapolated wbb.(:,:,3) components
        wbbl(:,:,3) = sqrt( wbbl(:,:,0)**2 + wbbl(:,:,1)**2 + wbbl(:,:,2)**2)
        wbbr(:,:,3) = sqrt( wbbr(:,:,0)**2 + wbbr(:,:,1)**2 + wbbr(:,:,2)**2)
        wbbv(:,:,3) = sqrt( wbbv(:,:,0)**2 + wbbv(:,:,1)**2 + wbbv(:,:,2)**2)

        !! On all unused cells, set magnetic field to some big positive value
        where( isUnusedCell( cflags(:, :, CELLFLAG_TYPE )))
            wbbl(:,:,0) = BB_INVALID
            wbbl(:,:,1) = BB_INVALID
            wbbl(:,:,2) = BB_INVALID
            wbbl(:,:,3) = BB_INVALID
            wbbr(:,:,0) = BB_INVALID
            wbbr(:,:,1) = BB_INVALID
            wbbr(:,:,2) = BB_INVALID
            wbbr(:,:,3) = BB_INVALID
            wbbv(:,:,0) = BB_INVALID
            wbbv(:,:,1) = BB_INVALID
            wbbv(:,:,2) = BB_INVALID
            wbbv(:,:,3) = BB_INVALID
        end where

        return
    end subroutine compute_b_at_faces

    !> This routine performs interpolation of a cell-centered quantity to cell
    !! faces with special treatment for triangular cells at the edges.
    !! The interpolation is done according to the cell widths given by hx and hy.
    !! xory.eq.1: interpolation in the x-direction
    !! xory.eq.2: interpolation in the y-direction
    !! Use this instead of left-right and top-bottom interpolation
    !! Obsolete: did not take properly into account trapezoidal and triangle cells
    !! Use interp_volume below instead
    subroutine interp_width( xory, nx, ny, hx, hy, gs, qc, centre, face )
        use b2mod_indirect_diffv_diffv
        implicit none
        integer :: xory
        integer :: nx   !< Specifies the number of interior cells
                        !< along the first coordinate (poloidal)
        integer :: ny   !< Specifies the number of interior cells
                        !< along the second coordinate (radial)
        real (kind=R8) :: hx(-1:nx,-1:ny)
        real (kind=R8) :: hy(-1:nx,-1:ny)
        real (kind=R8) :: gs(-1:nx,-1:ny,0:2)
        real (kind=R8) :: qc(-1:nx,-1:ny,0:1) !< Cosine of the angle
            !< between flux line direction and left cell face
        real (kind=R8) :: centre(-1:nx,-1:ny)
        real (kind=R8) :: face(-1:nx,-1:ny)

        integer :: ix   !< x-aligned (poloidal) cell index
        integer :: iy   !< y-aligned (radial) cell index
        real (kind=R8) :: area_to_top
        real (kind=R8) :: area_to_bottom
        real (kind=R8) :: area_to_left
        real (kind=R8) :: area_to_right
        intrinsic sqrt

        face = INVALID_DOUBLE

        do ix = -1, nx
          do iy = -1, ny
            if(isUnusedCell(cflags(ix,iy,CELLFLAG_TYPE))) cycle
                if(xory.eq.1) then
                    if(isInDomain(nx,ny,leftix(ix,iy),leftiy(ix,iy))) then
                        face(ix,iy) = (centre(ix,iy)*hx(leftix(ix,iy),leftiy(ix,iy))+ &
                        &  centre(leftix(ix,iy),leftiy(ix,iy))*hx(ix,iy))/ &
                            & (hx(ix,iy)+hx(leftix(ix,iy),leftiy(ix,iy)))
                    elseif ((isBoundaryCell(cflags(ix,iy,CELLFLAG_TYPE)) .or. &
                        &   isGhostCell(cflags(ix,iy,CELLFLAG_TYPE))) .and. &
                        &  cflags(ix,iy,CELLFLAG_TOPFACE) /= GRID_UNDEFINED .and. &
                        &  cflags(ix,iy,CELLFLAG_BOTTOMFACE) == GRID_UNDEFINED) then
                        face(ix,iy) = (centre(ix,iy)*hy(topix(ix,iy),topiy(ix,iy))+ &
                            &  centre(topix(ix,iy),topiy(ix,iy))*hy(ix,iy))/ &
                            & (hy(ix,iy)+hy(topix(ix,iy),topiy(ix,iy)))
                    elseif ((isBoundaryCell(cflags(ix,iy,CELLFLAG_TYPE)) .or. &
                        &   isGhostCell(cflags(ix,iy,CELLFLAG_TYPE))) .and. &
                        &  cflags(ix,iy,CELLFLAG_TOPFACE) == GRID_UNDEFINED .and. &
                        &  cflags(ix,iy,CELLFLAG_BOTTOMFACE) /= GRID_UNDEFINED) then
                        face(ix,iy) = (centre(ix,iy)*hy(bottomix(ix,iy),bottomiy(ix,iy))+ &
                            &  centre(bottomix(ix,iy),bottomiy(ix,iy))*hy(ix,iy))/ &
                            & (hy(ix,iy)+hy(bottomix(ix,iy),bottomiy(ix,iy)))
                    elseif ((isBoundaryCell(cflags(ix,iy,CELLFLAG_TYPE)) .or. &
                        &   isGhostCell(cflags(ix,iy,CELLFLAG_TYPE))) .and. &
                        &  cflags(ix,iy,CELLFLAG_TOPFACE) /= GRID_UNDEFINED .and. &
                        &  cflags(ix,iy,CELLFLAG_BOTTOMFACE) /= GRID_UNDEFINED) then
                        area_to_top = gs(topix(ix,iy),topiy(ix,iy),1)* &
                            &  sqrt(1.0_R8 - qc(topix(ix,iy),topiy(ix,iy),1)**2)
                        area_to_bottom = gs(ix,iy,1)*sqrt(1.0_R8 - qc(ix,iy,1)**2)
                        face(ix,iy) = (centre(ix,iy)* &
                            &  (area_to_bottom * hy(bottomix(ix,iy),bottomiy(ix,iy)) + &
                            &   area_to_top * hy(topix(ix,iy),topiy(ix,iy))) + &
                            &  (centre(bottomix(ix,iy),bottomiy(ix,iy)) * area_to_bottom + &
                            &   centre(topix(ix,iy),topiy(ix,iy)) * area_to_top)*hy(ix,iy))/ &
                            &  (area_to_bottom * hy(bottomix(ix,iy),bottomiy(ix,iy)) + &
                            &   area_to_top * hy(topix(ix,iy),topiy(ix,iy)) + &
                            &  (area_to_bottom+area_to_top) * hy(ix,iy))
                    endif
                elseif (xory.eq.2) then
                    if(isInDomain(nx,ny,bottomix(ix,iy),bottomiy(ix,iy))) then
                        face(ix,iy) = (centre(ix,iy)*hy(bottomix(ix,iy),bottomiy(ix,iy))+ &
                            &  centre(bottomix(ix,iy),bottomiy(ix,iy))*hy(ix,iy))/ &
                            & (hy(ix,iy)+hy(bottomix(ix,iy),bottomiy(ix,iy)))
                    elseif ((isBoundaryCell(cflags(ix,iy,CELLFLAG_TYPE)) .or. &
                        &   isGhostCell(cflags(ix,iy,CELLFLAG_TYPE))) .and. &
                        &  cflags(ix,iy,CELLFLAG_LEFTFACE) /= GRID_UNDEFINED .and. &
                        &  cflags(ix,iy,CELLFLAG_RIGHTFACE) == GRID_UNDEFINED) then
                        face(ix,iy) = (centre(ix,iy)*hx(leftix(ix,iy),leftiy(ix,iy))+ &
                            &  centre(leftix(ix,iy),leftiy(ix,iy))*hx(ix,iy))/ &
                            & (hx(ix,iy)+hx(leftix(ix,iy),leftiy(ix,iy)))
                    elseif ((isBoundaryCell(cflags(ix,iy,CELLFLAG_TYPE)) .or. &
                        &   isGhostCell(cflags(ix,iy,CELLFLAG_TYPE))) .and. &
                        &  cflags(ix,iy,CELLFLAG_LEFTFACE) == GRID_UNDEFINED .and. &
                        &  cflags(ix,iy,CELLFLAG_RIGHTFACE) /= GRID_UNDEFINED) then
                        face(ix,iy) = (centre(ix,iy)*hx(rightix(ix,iy),rightiy(ix,iy))+ &
                            &  centre(rightix(ix,iy),rightiy(ix,iy))*hx(ix,iy))/ &
                            & (hx(ix,iy)+hx(rightix(ix,iy),rightiy(ix,iy)))
                    elseif ((isBoundaryCell(cflags(ix,iy,CELLFLAG_TYPE)) .or. &
                        &   isGhostCell(cflags(ix,iy,CELLFLAG_TYPE))) .and. &
                        &  cflags(ix,iy,CELLFLAG_LEFTFACE) /= GRID_UNDEFINED .and. &
                        &  cflags(ix,iy,CELLFLAG_RIGHTFACE) /= GRID_UNDEFINED) then
                        area_to_right = gs(rightix(ix,iy),rightiy(ix,iy),0)* &
                            &  sqrt(1.0_R8 - qc(rightix(ix,iy),rightiy(ix,iy),0)**2)
                        area_to_left = gs(ix,iy,0)*sqrt(1.0_R8 - qc(ix,iy,0)**2)
                        face(ix,iy) = (centre(ix,iy)* &
                            &  (area_to_left * hx(leftix(ix,iy),leftiy(ix,iy)) + &
                            &   area_to_right * hx(rightix(ix,iy),rightiy(ix,iy))) + &
                            &  (centre(leftix(ix,iy),leftiy(ix,iy)) * area_to_left + &
                            &   centre(rightix(ix,iy),rightiy(ix,iy)) * area_to_right ) * &
                            &   hx(ix,iy)) / &
                            &  (area_to_left * hx(leftix(ix,iy),leftiy(ix,iy)) + &
                            &   area_to_right * hx(rightix(ix,iy),rightiy(ix,iy)) + &
                            &  (area_to_left+area_to_right) * hx(ix,iy))
                endif
            endif
          end do
        end do

        return
    end subroutine interp_width

  subroutine interp_volume1( idir, nx, ny, ns, vol, gs, qc, centre, face )
  implicit none
  integer, intent(in) :: idir
  integer, intent(in) :: nx !< Specifies the number of interior cells
                            !< along the first coordinate (poloidal)
  integer, intent(in) :: ny !< Specifies the number of interior cells
                            !< along the second coordinate (radial)
  integer, intent(in) :: ns
  real (kind=R8), intent(in) :: vol(-1:nx,-1:ny,0:4)    !< Cell volume
  real (kind=R8), intent(in) :: gs(-1:nx,-1:ny,0:2)
  real (kind=R8), intent(in) :: qc(-1:nx,-1:ny,0:1) !< Cosine of the angle
    !< between flux line direction and left cell face
  real (kind=R8), intent(in) :: centre(-1:nx,-1:ny,0:ns-1)
  real (kind=R8), intent(out) :: face(-1:nx,-1:ny,0:ns-1)
  integer is

  do is = 0, ns-1
    call interp_volume( idir, nx, ny, vol, gs, qc, centre(-1,-1,is),   &
        &   face(-1,-1,is) )
  end do

  return
  end subroutine interp_volume1

  !> This routine performs interpolation of a cell-centred quantity to cell faces
  !! with special treatment for triangular cells at the edges.
  !! This routine is to be used when the result must also be present on
  !! non-existent 'virtual' faces that are missing in triangle cells.
  !! If one only needs the values on existing valid faces, then one should use
  !! values_on_faces instead.
  !! The interpolation is done according to the cell volumes given by vol.
  !! For trapezoidal cells, the volume used is that which corresponds to
  !! the projected area of contact between the cells being interpolated.
  !! idir.eq.TO_LEFT: interpolation to the LEFT
  !! idir.eq.TO_BOTTOM: interpolation in the BOTTOM
  !! idir.eq.TO_RIGHT: interpolation to the RIGHT
  !! idir.eq.TO_TOP: interpolation in the TOP
  subroutine interp_volume( idir, nx, ny, vol, gs, qc, centre, face )
  use b2mod_indirect_diffv_diffv
  use b2mod_cellhelper
  implicit none
  integer, intent(in) :: idir
  integer, intent(in) :: nx    !< Specifies the number of interior cells
                               !< along the first coordinate (poloidal)
  integer, intent(in) :: ny    !< Specifies the number of interior cells
                               !< along the second coordinate (radial)
  real (kind=R8), intent(in) :: vol(-1:nx,-1:ny,0:4)
  real (kind=R8), intent(in) :: gs(-1:nx,-1:ny,0:2)
  real (kind=R8), intent(in) :: qc(-1:nx,-1:ny,0:1) !< Cosine of the angle
    !< between flux line direction and left cell face
  real (kind=R8), intent(in) :: centre(-1:nx,-1:ny)
  real (kind=R8), intent(out) :: face(-1:nx,-1:ny)

  integer ix, iy, ixn, iyn
  real (kind=R8) :: area_to_top, area_to_bottom, area_to_left, area_to_right
  intrinsic sqrt

!! Do not forget identical code in interp_magnetic_field below !!

  face = INVALID_DOUBLE

  do ix = -1, nx
    do iy = -1, ny
      if(isUnusedCell(cflags(ix,iy,CELLFLAG_TYPE))) cycle
      if (idir.eq.TO_LEFT) then
        if(isInDomain(nx,ny,leftix(ix,iy),leftiy(ix,iy))) then
           ixn = leftix(ix,iy)
           iyn = leftiy(ix,iy)
           face(ix,iy) = (centre(ix,iy)*vol(ixn,iyn,TO_RIGHT) + &
                       &  centre(ixn,iyn)*vol(ix,iy,TO_LEFT))/ &
                       & (vol(ixn,iyn,TO_RIGHT)+vol(ix,iy,TO_LEFT))
        elseif ((isBoundaryCell(cflags(ix,iy,CELLFLAG_TYPE)) .or. &
             &   isGhostCell(cflags(ix,iy,CELLFLAG_TYPE))) .and. &
             &  cflags(ix,iy,CELLFLAG_BOTTOMFACE) == GRID_UNDEFINED .and. &
             &  cflags(ix,iy,CELLFLAG_TOPFACE) /= GRID_UNDEFINED) then
           ixn = topix(ix,iy)
           iyn = topiy(ix,iy)
           face(ix,iy) = (centre(ix,iy)*vol(ixn,iyn,TO_BOTTOM)+ &
                       &  centre(ixn,iyn)*vol(ix,iy,TO_TOP))/ &
                       & (vol(ix,iy,TO_TOP)+vol(ixn,iyn,TO_BOTTOM))
        elseif ((isBoundaryCell(cflags(ix,iy,CELLFLAG_TYPE)) .or. &
             &   isGhostCell(cflags(ix,iy,CELLFLAG_TYPE))) .and. &
             &  cflags(ix,iy,CELLFLAG_TOPFACE) == GRID_UNDEFINED .and. &
             &  cflags(ix,iy,CELLFLAG_BOTTOMFACE) /= GRID_UNDEFINED) then
           ixn = bottomix(ix,iy)
           iyn = bottomiy(ix,iy)
           face(ix,iy) = (centre(ix,iy)*vol(ixn,iyn,TO_TOP)+ &
                       &  centre(ixn,iyn)*vol(ix,iy,TO_BOTTOM))/ &
                       & (vol(ix,iy,TO_BOTTOM)+vol(ixn,iyn,TO_TOP))
        elseif ((isBoundaryCell(cflags(ix,iy,CELLFLAG_TYPE)) .or. &
             &   isGhostCell(cflags(ix,iy,CELLFLAG_TYPE))) .and. &
             &  cflags(ix,iy,CELLFLAG_BOTTOMFACE) /= GRID_UNDEFINED .and. &
             &  cflags(ix,iy,CELLFLAG_TOPFACE) /= GRID_UNDEFINED) then
           area_to_top = gs(topix(ix,iy),topiy(ix,iy),1)* &
                 &  sqrt(1.0_R8 - qc(topix(ix,iy),topiy(ix,iy),1)**2)
           area_to_bottom = gs(ix,iy,1)*sqrt(1.0_R8 - qc(ix,iy,1)**2)
           face(ix,iy) = (centre(ix,iy)* &
             &  (area_to_bottom * vol(bottomix(ix,iy),bottomiy(ix,iy),TO_TOP) + &
             &   area_to_top * vol(topix(ix,iy),topiy(ix,iy),TO_BOTTOM)) + &
             &  (centre(bottomix(ix,iy),bottomiy(ix,iy)) * area_to_bottom + &
             &   centre(topix(ix,iy),topiy(ix,iy)) * area_to_top)* vol(ix,iy,TO_SELF))/ &
             &  (area_to_bottom * vol(bottomix(ix,iy),bottomiy(ix,iy),TO_TOP) + &
             &   area_to_top * vol(topix(ix,iy),topiy(ix,iy),TO_BOTTOM) + &
             &  (area_to_bottom+area_to_top) * vol(ix,iy,TO_SELF))
        endif
      elseif (idir.eq.TO_BOTTOM) then
        if(isInDomain(nx,ny,bottomix(ix,iy),bottomiy(ix,iy))) then
           ixn = bottomix(ix,iy)
           iyn = bottomiy(ix,iy)
           face(ix,iy) = (centre(ix,iy)*vol(ixn,iyn,TO_TOP) + &
                       &  centre(ixn,iyn)*vol(ix,iy,TO_BOTTOM))/ &
                       & (vol(ixn,iyn,TO_TOP)+vol(ix,iy,TO_BOTTOM))
        elseif ((isBoundaryCell(cflags(ix,iy,CELLFLAG_TYPE)) .or. &
             &   isGhostCell(cflags(ix,iy,CELLFLAG_TYPE))) .and. &
             &  cflags(ix,iy,CELLFLAG_LEFTFACE) /= GRID_UNDEFINED .and. &
             &  cflags(ix,iy,CELLFLAG_RIGHTFACE) == GRID_UNDEFINED) then
           ixn = leftix(ix,iy)
           iyn = leftiy(ix,iy)
           face(ix,iy) = (centre(ix,iy)*vol(ixn,iyn,TO_RIGHT) + &
                       &  centre(ixn,iyn)*vol(ix,iy,TO_LEFT))/ &
                       & (vol(ixn,iyn,TO_RIGHT)+vol(ix,iy,TO_LEFT))
        elseif ((isBoundaryCell(cflags(ix,iy,CELLFLAG_TYPE)) .or. &
             &   isGhostCell(cflags(ix,iy,CELLFLAG_TYPE))) .and. &
             &  cflags(ix,iy,CELLFLAG_LEFTFACE) == GRID_UNDEFINED .and. &
             &  cflags(ix,iy,CELLFLAG_RIGHTFACE) /= GRID_UNDEFINED) then
           ixn = rightix(ix,iy)
           iyn = rightiy(ix,iy)
           face(ix,iy) = (centre(ix,iy)*vol(ixn,iyn,TO_LEFT) + &
                       &  centre(ixn,iyn)*vol(ix,iy,TO_RIGHT))/ &
                       & (vol(ixn,iyn,TO_LEFT)+vol(ix,iy,TO_RIGHT))
        elseif ((isBoundaryCell(cflags(ix,iy,CELLFLAG_TYPE)) .or. &
             &   isGhostCell(cflags(ix,iy,CELLFLAG_TYPE))) .and. &
             &  cflags(ix,iy,CELLFLAG_LEFTFACE) /= GRID_UNDEFINED .and. &
             &  cflags(ix,iy,CELLFLAG_RIGHTFACE) /= GRID_UNDEFINED) then
           area_to_right = gs(rightix(ix,iy),rightiy(ix,iy),0)* &
                 &  sqrt(1.0_R8 - qc(rightix(ix,iy),rightiy(ix,iy),0)**2)
           area_to_left = gs(ix,iy,0)*sqrt(1.0_R8 - qc(ix,iy,0)**2)
           face(ix,iy) = (centre(ix,iy)* &
             &  (area_to_left * vol(leftix(ix,iy),leftiy(ix,iy),TO_RIGHT) + &
             &   area_to_right * vol(rightix(ix,iy),rightiy(ix,iy),TO_LEFT)) + &
             &  (centre(leftix(ix,iy),leftiy(ix,iy)) * area_to_left + &
             &   centre(rightix(ix,iy),rightiy(ix,iy)) * area_to_right)* vol(ix,iy,TO_SELF))/ &
             &  (area_to_left * vol(leftix(ix,iy),leftiy(ix,iy),TO_RIGHT) + &
             &   area_to_right * vol(rightix(ix,iy),rightiy(ix,iy),TO_LEFT) + &
             &  (area_to_left+area_to_right) * vol(ix,iy,TO_SELF))
        endif
      elseif (idir.eq.TO_RIGHT) then
        if(isInDomain(nx,ny,rightix(ix,iy),rightiy(ix,iy))) then
           ixn = rightix(ix,iy)
           iyn = rightiy(ix,iy)
           face(ix,iy) = (centre(ix,iy)*vol(ixn,iyn,TO_LEFT) + &
                       &  centre(ixn,iyn)*vol(ix,iy,TO_RIGHT))/ &
                       & (vol(ixn,iyn,TO_LEFT)+vol(ix,iy,TO_RIGHT))
        elseif ((isBoundaryCell(cflags(ix,iy,CELLFLAG_TYPE)) .or. &
             &   isGhostCell(cflags(ix,iy,CELLFLAG_TYPE))) .and. &
             &  cflags(ix,iy,CELLFLAG_BOTTOMFACE) == GRID_UNDEFINED .and. &
             &  cflags(ix,iy,CELLFLAG_TOPFACE) /= GRID_UNDEFINED) then
           ixn = topix(ix,iy)
           iyn = topiy(ix,iy)
           face(ix,iy) = (centre(ix,iy)*vol(ixn,iyn,TO_BOTTOM)+ &
                       &  centre(ixn,iyn)*vol(ix,iy,TO_TOP))/ &
                       & (vol(ix,iy,TO_TOP)+vol(ixn,iyn,TO_BOTTOM))
        elseif ((isBoundaryCell(cflags(ix,iy,CELLFLAG_TYPE)) .or. &
             &   isGhostCell(cflags(ix,iy,CELLFLAG_TYPE))) .and. &
             &  cflags(ix,iy,CELLFLAG_TOPFACE) == GRID_UNDEFINED .and. &
             &  cflags(ix,iy,CELLFLAG_BOTTOMFACE) /= GRID_UNDEFINED) then
           ixn = bottomix(ix,iy)
           iyn = bottomiy(ix,iy)
           face(ix,iy) = (centre(ix,iy)*vol(ixn,iyn,TO_TOP)+ &
                       &  centre(ixn,iyn)*vol(ix,iy,TO_BOTTOM))/ &
                       & (vol(ix,iy,TO_BOTTOM)+vol(ixn,iyn,TO_TOP))
        elseif ((isBoundaryCell(cflags(ix,iy,CELLFLAG_TYPE)) .or. &
             &   isGhostCell(cflags(ix,iy,CELLFLAG_TYPE))) .and. &
             &  cflags(ix,iy,CELLFLAG_BOTTOMFACE) /= GRID_UNDEFINED .and. &
             &  cflags(ix,iy,CELLFLAG_TOPFACE) /= GRID_UNDEFINED) then
           area_to_top = gs(topix(ix,iy),topiy(ix,iy),1)* &
                 &  sqrt(1.0_R8 - qc(topix(ix,iy),topiy(ix,iy),1)**2)
           area_to_bottom = gs(ix,iy,1)*sqrt(1.0_R8 - qc(ix,iy,1)**2)
           face(ix,iy) = (centre(ix,iy)* &
             &  (area_to_bottom * vol(bottomix(ix,iy),bottomiy(ix,iy),TO_TOP) + &
             &   area_to_top * vol(topix(ix,iy),topiy(ix,iy),TO_BOTTOM)) + &
             &  (centre(bottomix(ix,iy),bottomiy(ix,iy)) * area_to_bottom + &
             &   centre(topix(ix,iy),topiy(ix,iy)) * area_to_top)* vol(ix,iy,TO_SELF))/ &
             &  (area_to_bottom * vol(bottomix(ix,iy),bottomiy(ix,iy),TO_TOP) + &
             &   area_to_top * vol(topix(ix,iy),topiy(ix,iy),TO_BOTTOM) + &
             &  (area_to_bottom+area_to_top) * vol(ix,iy,TO_SELF))
        endif
      elseif (idir.eq.TO_TOP) then
        if(isInDomain(nx,ny,topix(ix,iy),topiy(ix,iy))) then
           ixn = topix(ix,iy)
           iyn = topiy(ix,iy)
           face(ix,iy) = (centre(ix,iy)*vol(ixn,iyn,TO_BOTTOM) + &
                       &  centre(ixn,iyn)*vol(ix,iy,TO_TOP))/ &
                       & (vol(ixn,iyn,TO_BOTTOM)+vol(ix,iy,TO_TOP))
        elseif ((isBoundaryCell(cflags(ix,iy,CELLFLAG_TYPE)) .or. &
             &   isGhostCell(cflags(ix,iy,CELLFLAG_TYPE))) .and. &
             &  cflags(ix,iy,CELLFLAG_LEFTFACE) /= GRID_UNDEFINED .and. &
             &  cflags(ix,iy,CELLFLAG_RIGHTFACE) == GRID_UNDEFINED) then
           ixn = leftix(ix,iy)
           iyn = leftiy(ix,iy)
           face(ix,iy) = (centre(ix,iy)*vol(ixn,iyn,TO_RIGHT) + &
                       &  centre(ixn,iyn)*vol(ix,iy,TO_LEFT))/ &
                       & (vol(ixn,iyn,TO_RIGHT)+vol(ix,iy,TO_LEFT))
        elseif ((isBoundaryCell(cflags(ix,iy,CELLFLAG_TYPE)) .or. &
             &   isGhostCell(cflags(ix,iy,CELLFLAG_TYPE))) .and. &
             &  cflags(ix,iy,CELLFLAG_LEFTFACE) == GRID_UNDEFINED .and. &
             &  cflags(ix,iy,CELLFLAG_RIGHTFACE) /= GRID_UNDEFINED) then
           ixn = rightix(ix,iy)
           iyn = rightiy(ix,iy)
           face(ix,iy) = (centre(ix,iy)*vol(ixn,iyn,TO_RIGHT) + &
                       &  centre(ixn,iyn)*vol(ix,iy,TO_LEFT))/ &
                       & (vol(ixn,iyn,TO_LEFT)+vol(ix,iy,TO_RIGHT))
        elseif ((isBoundaryCell(cflags(ix,iy,CELLFLAG_TYPE)) .or. &
             &   isGhostCell(cflags(ix,iy,CELLFLAG_TYPE))) .and. &
             &  cflags(ix,iy,CELLFLAG_LEFTFACE) /= GRID_UNDEFINED .and. &
             &  cflags(ix,iy,CELLFLAG_RIGHTFACE) /= GRID_UNDEFINED) then
           area_to_right = gs(rightix(ix,iy),rightiy(ix,iy),0)* &
                 &  sqrt(1.0_R8 - qc(rightix(ix,iy),rightiy(ix,iy),0)**2)
           area_to_left = gs(ix,iy,0)*sqrt(1.0_R8 - qc(ix,iy,0)**2)
           face(ix,iy) = (centre(ix,iy)* &
             &  (area_to_left * vol(leftix(ix,iy),leftiy(ix,iy),TO_RIGHT) + &
             &   area_to_right * vol(rightix(ix,iy),rightiy(ix,iy),TO_LEFT)) + &
             &  (centre(leftix(ix,iy),leftiy(ix,iy)) * area_to_left + &
             &   centre(rightix(ix,iy),rightiy(ix,iy)) * area_to_right)* vol(ix,iy,TO_SELF))/ &
             &  (area_to_left * vol(leftix(ix,iy),leftiy(ix,iy),TO_RIGHT) + &
             &   area_to_right * vol(rightix(ix,iy),rightiy(ix,iy),TO_LEFT) + &
             &  (area_to_left+area_to_right) * vol(ix,iy,TO_SELF))
        endif
      endif
    end do
  end do

  return
  end subroutine interp_volume

  !> Identical code to interp_volume, but needs to refer to the
  !! geometric arrays explicitly
  subroutine interp_magnetic_field( nx, ny, bb, wbbl, wbbr, wbbv,   &
        &    leftix, leftiy, rightix, rightiy, bottomix, bottomiy,  &
        &    topix, topiy, vol, gs, qc, cflags)

    implicit none
    integer, intent(in) :: nx   !< Specifies the number of interior cells
                                !< along the first coordinate (poloidal)
    integer, intent(in) :: ny   !< Specifies the number of interior cells
                                !< along the second coordinate (radial)
    real(R8), intent(in) :: bb(-1:nx,-1:ny,0:3)
    real(R8), intent(in) :: vol(-1:nx,-1:ny,0:4)
    real(R8), intent(in) :: gs(-1:nx,-1:ny,0:2)
    real(R8), intent(in) :: qc(-1:nx,-1:ny,0:1)
    real(R8), intent(out) :: wbbl(-1:nx,-1:ny,0:3)
    real(R8), intent(out) :: wbbr(-1:nx,-1:ny,0:3)
    real(R8), intent(out) :: wbbv(-1:nx,-1:ny,0:3)
    integer, intent(in) :: cflags(-1:nx,-1:ny,CARREOUT_NCELLFLAGS)
    integer, intent(in) :: leftix( -1:nx, -1:ny)    !< Left neighbour
        !< poloidal (first coordinate) index array
    integer, intent(in) :: leftiy( -1:nx, -1:ny )   !< Left neighbour radial
        !< (second coordinate) index
    integer, intent(in) :: rightix( -1:nx, -1:ny )  !< Right neighbour
        !< poloidal (first coordinate) index array
    integer, intent(in) :: rightiy( -1:nx, -1:ny )  !< Right neighbour
        !< radial (second coordinate) index
    integer, intent(in) :: topix( -1:nx, -1:ny )    !< Top neighbour
        !< poloidal (first coordinate) index array
    integer, intent(in) :: topiy( -1:nx, -1:ny )    !< Top neighbour radial
        !< (second coordinate) index
    integer, intent(in) :: bottomix( -1:nx, -1:ny ) !< Bottom neighbour
        !< poloidal (first coordinate) index array
    integer, intent(in) :: bottomiy( -1:nx, -1:ny ) !< Bottom neighbour
        !< radial (second coordinate) index

    integer :: ix   !< x-aligned (poloidal) cell index
    integer :: iy   !< y-aligned (radial) cell index
    integer :: ixn
    integer :: iyn
    real (kind=R8) :: area_to_top, area_to_bottom, area_to_left, area_to_right
  intrinsic sqrt

  wbbl = BB_INVALID
  wbbr = BB_INVALID
  wbbv = BB_INVALID

  do ix = -1, nx
    do iy = -1, ny
      if(isUnusedCell(cflags(ix,iy,CELLFLAG_TYPE))) cycle
!! First compute wbbl, the interpolation to the left
      if(isInDomain(nx,ny,leftix(ix,iy),leftiy(ix,iy))) then
        ixn = leftix(ix,iy)
        iyn = leftiy(ix,iy)
        wbbl(ix,iy,:) = (bb(ix,iy,:)*vol(ixn,iyn,TO_RIGHT) + &
                      &  bb(ixn,iyn,:)*vol(ix,iy,TO_LEFT))/ &
                      & (vol(ixn,iyn,TO_RIGHT)+vol(ix,iy,TO_LEFT))
      elseif ((isBoundaryCell(cflags(ix,iy,CELLFLAG_TYPE)) .or. &
           &   isGhostCell(cflags(ix,iy,CELLFLAG_TYPE))) .and. &
           &  cflags(ix,iy,CELLFLAG_BOTTOMFACE) == GRID_UNDEFINED .and. &
           &  cflags(ix,iy,CELLFLAG_TOPFACE) /= GRID_UNDEFINED) then
        ixn = topix(ix,iy)
        iyn = topiy(ix,iy)
        wbbl(ix,iy,:) = (bb(ix,iy,:)*vol(ixn,iyn,TO_BOTTOM)+ &
                     &  bb(ixn,iyn,:)*vol(ix,iy,TO_TOP))/ &
                     & (vol(ix,iy,TO_TOP)+vol(ixn,iyn,TO_BOTTOM))
      elseif ((isBoundaryCell(cflags(ix,iy,CELLFLAG_TYPE)) .or. &
           &   isGhostCell(cflags(ix,iy,CELLFLAG_TYPE))) .and. &
           &  cflags(ix,iy,CELLFLAG_TOPFACE) == GRID_UNDEFINED .and. &
           &  cflags(ix,iy,CELLFLAG_BOTTOMFACE) /= GRID_UNDEFINED) then
        ixn = bottomix(ix,iy)
        iyn = bottomiy(ix,iy)
        wbbl(ix,iy,:) = (bb(ix,iy,:)*vol(ixn,iyn,TO_TOP)+ &
                      &  bb(ixn,iyn,:)*vol(ix,iy,TO_BOTTOM))/ &
                      & (vol(ix,iy,TO_BOTTOM)+vol(ixn,iyn,TO_TOP))
      elseif ((isBoundaryCell(cflags(ix,iy,CELLFLAG_TYPE)) .or. &
           &   isGhostCell(cflags(ix,iy,CELLFLAG_TYPE))) .and. &
           &  cflags(ix,iy,CELLFLAG_BOTTOMFACE) /= GRID_UNDEFINED .and. &
           &  cflags(ix,iy,CELLFLAG_TOPFACE) /= GRID_UNDEFINED) then
        area_to_top = gs(topix(ix,iy),topiy(ix,iy),1)* &
               &  sqrt(1.0_R8 - qc(topix(ix,iy),topiy(ix,iy),1)**2)
        area_to_bottom = gs(ix,iy,1)*sqrt(1.0_R8 - qc(ix,iy,1)**2)
        wbbl(ix,iy,:) = (bb(ix,iy,:)* &
           &  (area_to_bottom * vol(bottomix(ix,iy),bottomiy(ix,iy),TO_TOP) + &
           &   area_to_top * vol(topix(ix,iy),topiy(ix,iy),TO_BOTTOM)) + &
           &  (bb(bottomix(ix,iy),bottomiy(ix,iy),:) * area_to_bottom + &
           &   bb(topix(ix,iy),topiy(ix,iy),:) * area_to_top)* vol(ix,iy,TO_SELF))/ &
           &  (area_to_bottom * vol(bottomix(ix,iy),bottomiy(ix,iy),TO_TOP) + &
           &   area_to_top * vol(topix(ix,iy),topiy(ix,iy),TO_BOTTOM) + &
           &  (area_to_bottom+area_to_top) * vol(ix,iy,TO_SELF))
      endif
!! Second, compute wbbr, the interpolation to the right
      if(isInDomain(nx,ny,rightix(ix,iy),rightiy(ix,iy))) then
        ixn = rightix(ix,iy)
        iyn = rightiy(ix,iy)
        wbbr(ix,iy,:) = (bb(ix,iy,:)*vol(ixn,iyn,TO_LEFT) + &
                    &  bb(ixn,iyn,:)*vol(ix,iy,TO_RIGHT))/ &
                    & (vol(ixn,iyn,TO_LEFT)+vol(ix,iy,TO_RIGHT))
      elseif ((isBoundaryCell(cflags(ix,iy,CELLFLAG_TYPE)) .or. &
           &   isGhostCell(cflags(ix,iy,CELLFLAG_TYPE))) .and. &
           &  cflags(ix,iy,CELLFLAG_BOTTOMFACE) == GRID_UNDEFINED .and. &
           &  cflags(ix,iy,CELLFLAG_TOPFACE) /= GRID_UNDEFINED) then
        ixn = topix(ix,iy)
        iyn = topiy(ix,iy)
        wbbr(ix,iy,:) = (bb(ix,iy,:)*vol(ixn,iyn,TO_BOTTOM)+ &
                    &  bb(ixn,iyn,:)*vol(ix,iy,TO_TOP))/ &
                    & (vol(ix,iy,TO_TOP)+vol(ixn,iyn,TO_BOTTOM))
      elseif ((isBoundaryCell(cflags(ix,iy,CELLFLAG_TYPE)) .or. &
           &   isGhostCell(cflags(ix,iy,CELLFLAG_TYPE))) .and. &
           &  cflags(ix,iy,CELLFLAG_TOPFACE) == GRID_UNDEFINED .and. &
           &  cflags(ix,iy,CELLFLAG_BOTTOMFACE) /= GRID_UNDEFINED) then
        ixn = bottomix(ix,iy)
        iyn = bottomiy(ix,iy)
        wbbr(ix,iy,:) = (bb(ix,iy,:)*vol(ixn,iyn,TO_TOP)+ &
                    &  bb(ixn,iyn,:)*vol(ix,iy,TO_BOTTOM))/ &
                    & (vol(ix,iy,TO_BOTTOM)+vol(ixn,iyn,TO_TOP))
      elseif ((isBoundaryCell(cflags(ix,iy,CELLFLAG_TYPE)) .or. &
           &   isGhostCell(cflags(ix,iy,CELLFLAG_TYPE))) .and. &
           &  cflags(ix,iy,CELLFLAG_BOTTOMFACE) /= GRID_UNDEFINED .and. &
           &  cflags(ix,iy,CELLFLAG_TOPFACE) /= GRID_UNDEFINED) then
        area_to_top = gs(topix(ix,iy),topiy(ix,iy),1)* &
               &  sqrt(1.0_R8 - qc(topix(ix,iy),topiy(ix,iy),1)**2)
        area_to_bottom = gs(ix,iy,1)*sqrt(1.0_R8 - qc(ix,iy,1)**2)
        wbbr(ix,iy,:) = (bb(ix,iy,:)* &
           &  (area_to_bottom * vol(bottomix(ix,iy),bottomiy(ix,iy),TO_TOP) + &
           &   area_to_top * vol(topix(ix,iy),topiy(ix,iy),TO_BOTTOM)) + &
           &  (bb(bottomix(ix,iy),bottomiy(ix,iy),:) * area_to_bottom + &
           &   bb(topix(ix,iy),topiy(ix,iy),:) * area_to_top)* vol(ix,iy,TO_SELF))/ &
           &  (area_to_bottom * vol(bottomix(ix,iy),bottomiy(ix,iy),TO_TOP) + &
           &   area_to_top * vol(topix(ix,iy),topiy(ix,iy),TO_BOTTOM) + &
           &  (area_to_bottom+area_to_top) * vol(ix,iy,TO_SELF))
      endif
!! Third, compute wbbv, the interpolation to the bottom
      if(isInDomain(nx,ny,bottomix(ix,iy),bottomiy(ix,iy))) then
        ixn = bottomix(ix,iy)
        iyn = bottomiy(ix,iy)
        wbbv(ix,iy,:) = (bb(ix,iy,:)*vol(ixn,iyn,TO_TOP) + &
                    &  bb(ixn,iyn,:)*vol(ix,iy,TO_BOTTOM))/ &
                    & (vol(ixn,iyn,TO_TOP)+vol(ix,iy,TO_BOTTOM))
      elseif ((isBoundaryCell(cflags(ix,iy,CELLFLAG_TYPE)) .or. &
           &   isGhostCell(cflags(ix,iy,CELLFLAG_TYPE))) .and. &
           &  cflags(ix,iy,CELLFLAG_LEFTFACE) /= GRID_UNDEFINED .and. &
           &  cflags(ix,iy,CELLFLAG_RIGHTFACE) == GRID_UNDEFINED) then
        ixn = leftix(ix,iy)
        iyn = leftiy(ix,iy)
        wbbv(ix,iy,:) = (bb(ix,iy,:)*vol(ixn,iyn,TO_RIGHT) + &
                    &  bb(ixn,iyn,:)*vol(ix,iy,TO_LEFT))/ &
                    & (vol(ixn,iyn,TO_RIGHT)+vol(ix,iy,TO_LEFT))
      elseif ((isBoundaryCell(cflags(ix,iy,CELLFLAG_TYPE)) .or. &
           &   isGhostCell(cflags(ix,iy,CELLFLAG_TYPE))) .and. &
           &  cflags(ix,iy,CELLFLAG_LEFTFACE) == GRID_UNDEFINED .and. &
           &  cflags(ix,iy,CELLFLAG_RIGHTFACE) /= GRID_UNDEFINED) then
        ixn = rightix(ix,iy)
        iyn = rightiy(ix,iy)
        wbbv(ix,iy,:) = (bb(ix,iy,:)*vol(ixn,iyn,TO_LEFT) + &
                    &  bb(ixn,iyn,:)*vol(ix,iy,TO_RIGHT))/ &
                    & (vol(ixn,iyn,TO_LEFT)+vol(ix,iy,TO_RIGHT))
      elseif ((isBoundaryCell(cflags(ix,iy,CELLFLAG_TYPE)) .or. &
           &   isGhostCell(cflags(ix,iy,CELLFLAG_TYPE))) .and. &
           &  cflags(ix,iy,CELLFLAG_LEFTFACE) /= GRID_UNDEFINED .and. &
           &  cflags(ix,iy,CELLFLAG_RIGHTFACE) /= GRID_UNDEFINED) then
        area_to_right = gs(rightix(ix,iy),rightiy(ix,iy),0)* &
           &  sqrt(1.0_R8 - qc(rightix(ix,iy),rightiy(ix,iy),0)**2)
        area_to_left = gs(ix,iy,0)*sqrt(1.0_R8 - qc(ix,iy,0)**2)
        wbbv(ix,iy,:) = (bb(ix,iy,:)* &
           &  (area_to_left * vol(leftix(ix,iy),leftiy(ix,iy),TO_RIGHT) + &
           &   area_to_right * vol(rightix(ix,iy),rightiy(ix,iy),TO_LEFT)) + &
           &  (bb(leftix(ix,iy),leftiy(ix,iy),:) * area_to_left + &
           &   bb(rightix(ix,iy),rightiy(ix,iy),:) * area_to_right)* vol(ix,iy,TO_SELF))/ &
           &  (area_to_left * vol(leftix(ix,iy),leftiy(ix,iy),TO_RIGHT) + &
           &   area_to_right * vol(rightix(ix,iy),rightiy(ix,iy),TO_LEFT) + &
           &  (area_to_left+area_to_right) * vol(ix,iy,TO_SELF))
      endif
    end do
  end do

  return
  end subroutine interp_magnetic_field

#ifdef ITM_ENVIRONMENT_LOADED
  subroutine value_on_faces1(nx,ny,ns,weight,centre,face)
  implicit none
  integer, intent(in) :: nx, ny, ns
  real (kind=R8), intent(in) :: centre(-1:nx,-1:ny,0:ns-1), &
                              & weight(-1:nx,-1:ny,0:4)
  real (kind=R8), intent(out) :: face(-1:nx,-1:ny,0:1,0:ns-1)
  integer is

  do is = 0, ns-1
    call value_on_faces(nx,ny,weight,centre(-1,-1,is),face(-1,-1,0,is))
  end do

  return
  end subroutine value_on_faces1

  !> This subroutine computes an interpolated value on the existing faces only
  subroutine value_on_faces(nx,ny,weight,centre,face)
  use b2mod_geo_diffv_diffv
  use b2mod_indirect_diffv_diffv
  implicit none
!! input arguments
  integer, intent(in) :: nx, ny
  real (kind=R8), intent(in) :: centre(-1:nx,-1:ny), &
                              & weight(-1:nx,-1:ny,0:4)
!! output arguments
  real (kind=R8), intent(out) :: face(-1:nx,-1:ny,0:1)
!! local variables
  integer ix, iy, cgeo
  real (kind=R8) :: cxx(0:3), cyy(0:3)

  face = 0.0_R8

  do iy = -1, ny
    do ix = -1, nx
      if (isUnusedCell(cflags(ix,iy,CELLFLAG_TYPE))) cycle
      cxx(0:3) = crx(ix,iy,0:3)
      cyy(0:3) = cry(ix,iy,0:3)
      cgeo = cellGeoType( cxx, cyy )
      if (cgeo == CGEO_TRIA_NOLEFT) then
        face(ix,iy,0) = 0.0_R8
      else if (.not.isInDomain(nx,ny,leftix(ix,iy),leftiy(ix,iy))) then
        face(ix,iy,0) = centre(ix,iy)
      else
        face(ix,iy,0) = &
         &  (centre(ix,iy)*weight(leftix(ix,iy),leftiy(ix,iy),TO_RIGHT) + &
         &   centre(leftix(ix,iy),leftiy(ix,iy))*weight(ix,iy,TO_LEFT)) / &
         &  (weight(leftix(ix,iy),leftiy(ix,iy),TO_RIGHT)+weight(ix,iy,TO_LEFT))
      end if
      if (cgeo == CGEO_TRIA_NOBOT) then
        face(ix,iy,1) = 0.0_R8
      else if (.not.isInDomain(nx,ny,bottomix(ix,iy),bottomiy(ix,iy))) then
        face(ix,iy,1) = centre(ix,iy)
      else
        face(ix,iy,1) = &
         &  (centre(ix,iy)*weight(bottomix(ix,iy),bottomiy(ix,iy),TO_TOP) + &
         &   centre(bottomix(ix,iy),bottomiy(ix,iy))*weight(ix,iy,TO_BOTTOM)) / &
         &  (weight(bottomix(ix,iy),bottomiy(ix,iy),TO_TOP)+ &
         &   weight(ix,iy,TO_BOTTOM))
      end if
    end do
  end do

  return
  end subroutine value_on_faces

  !> Interpolate a cell quantity cv to a particular vertex
  function interpolateToVertices( nx, ny, vx_index, cv ) result( vx )
    use b2mod_geo_diffv_diffv
    use b2mod_b2cmfs
    use b2mod_indirect_diffv_diffv
    use b2mod_constants
    use b2mod_cellhelper
    implicit none
    integer, intent(in) :: nx, ny, vx_index
    real (kind=R8), intent(in), dimension(-1:nx,-1:ny) :: cv
    real (kind=R8), dimension(-1:nx,-1:ny) :: vx

    ! internal
    integer :: ix, iy, is, ixx, iyy
    real (kind=R8) :: av, wTot, minVal, maxVal, w, d, centroid(0:1)

    ! procedures
    external xertst

    ! For every vertex, find connected cells and average values
    do ix = -1, nx
       do iy = -1, ny

          if(isUnusedCell(cflags(ix,iy,CELLFLAG_TYPE))) then
             vx(ix,iy) = cv(ix,iy)
             cycle
          end if

          wTot = 0.0_R8
          av = 0.0_R8
          minVal = huge(1.0_R8)
          maxVal =-huge(1.0_R8)
          do is = 1, 8
            if ( gmap%mapCvixVx( gmap%mapVxI(ix,iy,vx_index), is ) == B2_GRID_UNDEFINED .or. &
               & gmap%mapCviyVx( gmap%mapVxI(ix,iy,vx_index), is ) == B2_GRID_UNDEFINED) cycle
            ixx = gmap%mapCvixVx( gmap%mapVxI(ix,iy,vx_index), is )
            iyy = gmap%mapCviyVx( gmap%mapVxI(ix,iy,vx_index), is )
            centroid = quadCentroid(crx(ixx,iyy,0),cry(ixx,iyy,0), &
                                  & crx(ixx,iyy,1),cry(ixx,iyy,1), &
                                  & crx(ixx,iyy,2),cry(ixx,iyy,2), &
                                  & crx(ixx,iyy,3),cry(ixx,iyy,3))
            d = sqrt ( (centroid(0) - crx(ix,iy,vx_index))**2 + &
                     & (centroid(1) - cry(ix,iy,vx_index))**2 )
            w = 1.0_R8/d
            av = av + cv(ixx, iyy) * w
            wTot = wTot + w
            minVal = min(minVal, cv(ixx, iyy))
            maxVal = max(maxVal, cv(ixx, iyy))

          end do

          vx(ix, iy) = av / wTot
          d = vx(ix,iy) - minVal
          if (minVal.ne.0.0_R8) then
            call xertst ( d/abs(minVal) >= - 1.0e-15_R8, "Interpolated value smaller than minimum" )
          else
            call xertst ( d.ge.minVal, "Interpolated value smaller than minimum" )
          end if
          d = vx(ix,iy) - maxVal
          if (maxVal.ne.0.0_R8) then
            call xertst ( d/abs(maxVal) <=   1.0e-15_R8, "Interpolated value bigger than maximum" )
          else
            call xertst ( d.le.maxVal, "Interpolated value bigger than maximum" )
          end if
       end do
    end do
    return

  end function interpolateToVertices
#endif

end module b2mod_interp

!!!Local Variables:
!!! mode: f90
!!! End:
