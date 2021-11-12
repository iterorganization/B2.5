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
        &   topix, topiy, vol, gs, qc, qcb, cflags )

        implicit none
        integer, intent(in)  :: nx  !< Specifies the number of interior cells
                                    !< along the first coordinate
        integer, intent(in)  :: ny  !< Specifies the number of interior cells
                                    !< along the second coordinate
        real(R8), intent(in) :: bb( -1:nx, -1:ny, 0:3 ) !< Magnetic field
        real(R8), intent(in) :: vol( -1:nx ,-1:ny ,0:4 )    !< Cell volume
        real(R8), intent(in) :: gs( -1:nx, -1:ny, 0:2 )
        real(R8), intent(in) :: qc( -1:nx, -1:ny ) !< Cosine of the angle
            !< between flux line direction and left cell face
        real(R8), intent(in) :: qcb( -1:nx, -1:ny )
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
            &   topix, topiy, vol, gs, qc, qcb, cflags )

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
    subroutine interp_width( xory, nx, ny, hx, hy, gs, qc, qcb, centre, face )
        use b2mod_indirect
        implicit none
        integer :: xory
        integer :: nx   !< Specifies the number of interior cells
                        !< along the first coordinate (poloidal)
        integer :: ny   !< Specifies the number of interior cells
                        !< along the second coordinate (radial)
        real (kind=R8) :: hx(-1:nx,-1:ny)
        real (kind=R8) :: hy(-1:nx,-1:ny)
        real (kind=R8) :: gs(-1:nx,-1:ny,0:2)
        real (kind=R8) :: qc(-1:nx,-1:ny) !< Cosine of the angle
            !< between flux line direction and left cell face
        real (kind=R8) :: qcb(-1:nx,-1:ny)
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
                            &  sqrt(1.0_R8 - qcb(topix(ix,iy),topiy(ix,iy))**2)
                        area_to_bottom = gs(ix,iy,1)*sqrt(1.0_R8 - qcb(ix,iy)**2)
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
                            &  sqrt(1.0_R8 - qc(rightix(ix,iy),rightiy(ix,iy))**2)
                        area_to_left = gs(ix,iy,0)*sqrt(1.0_R8 - qc(ix,iy)**2)
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

  subroutine interp_volume1( idir, nx, ny, ns, vol, gs, qc, qcb, centre, face )
  implicit none
  integer, intent(in) :: idir
  integer, intent(in) :: nx !< Specifies the number of interior cells
                            !< along the first coordinate (poloidal)
  integer, intent(in) :: ny !< Specifies the number of interior cells
                            !< along the second coordinate (radial)
  integer, intent(in) :: ns
  real (kind=R8), intent(in) :: vol(-1:nx,-1:ny,0:4)    !< Cell volume
  real (kind=R8), intent(in) :: gs(-1:nx,-1:ny,0:2)
  real (kind=R8), intent(in) :: qc(-1:nx,-1:ny) !< Cosine of the angle
    !< between flux line direction and left cell face
  real (kind=R8), intent(in) :: qcb(-1:nx,-1:ny)
  real (kind=R8), intent(in) :: centre(-1:nx,-1:ny,0:ns-1)
  real (kind=R8), intent(out) :: face(-1:nx,-1:ny,0:ns-1)
  integer is

  do is = 0, ns-1
    call interp_volume( idir, nx, ny, vol, gs, qc, qcb, centre(-1,-1,is),   &
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
  subroutine interp_volume( idir, nx, ny, vol, gs, qc, qcb, centre, face )
  use b2mod_indirect
  use b2mod_cellhelper
  implicit none
  integer, intent(in) :: idir
  integer, intent(in) :: nx    !< Specifies the number of interior cells
                               !< along the first coordinate (poloidal)
  integer, intent(in) :: ny    !< Specifies the number of interior cells
                               !< along the second coordinate (radial)
  real (kind=R8), intent(in) :: vol(-1:nx,-1:ny,0:4)
  real (kind=R8), intent(in) :: gs(-1:nx,-1:ny,0:2)
  real (kind=R8), intent(in) :: qc(-1:nx,-1:ny) !< Cosine of the angle
    !< between flux line direction and left cell face
  real (kind=R8), intent(in) :: qcb(-1:nx,-1:ny)
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
                 &  sqrt(1.0_R8 - qcb(topix(ix,iy),topiy(ix,iy))**2)
           area_to_bottom = gs(ix,iy,1)*sqrt(1.0_R8 - qcb(ix,iy)**2)
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
                 &  sqrt(1.0_R8 - qc(rightix(ix,iy),rightiy(ix,iy))**2)
           area_to_left = gs(ix,iy,0)*sqrt(1.0_R8 - qc(ix,iy)**2)
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
                 &  sqrt(1.0_R8 - qcb(topix(ix,iy),topiy(ix,iy))**2)
           area_to_bottom = gs(ix,iy,1)*sqrt(1.0_R8 - qcb(ix,iy)**2)
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
                 &  sqrt(1.0_R8 - qc(rightix(ix,iy),rightiy(ix,iy))**2)
           area_to_left = gs(ix,iy,0)*sqrt(1.0_R8 - qc(ix,iy)**2)
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
        &    topix, topiy, vol, gs, qc, qcb, cflags)

    implicit none
    integer, intent(in) :: nx   !< Specifies the number of interior cells
                                !< along the first coordinate (poloidal)
    integer, intent(in) :: ny   !< Specifies the number of interior cells
                                !< along the second coordinate (radial)
    real(R8), intent(in) :: bb(-1:nx,-1:ny,0:3)
    real(R8), intent(in) :: vol(-1:nx,-1:ny,0:4)
    real(R8), intent(in) :: gs(-1:nx,-1:ny,0:2)
    real(R8), intent(in) :: qc(-1:nx,-1:ny)
    real(R8), intent(in) :: qcb(-1:nx,-1:ny)
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
               &  sqrt(1.0_R8 - qcb(topix(ix,iy),topiy(ix,iy))**2)
        area_to_bottom = gs(ix,iy,1)*sqrt(1.0_R8 - qcb(ix,iy)**2)
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
               &  sqrt(1.0_R8 - qcb(topix(ix,iy),topiy(ix,iy))**2)
        area_to_bottom = gs(ix,iy,1)*sqrt(1.0_R8 - qcb(ix,iy)**2)
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
           &  sqrt(1.0_R8 - qc(rightix(ix,iy),rightiy(ix,iy))**2)
        area_to_left = gs(ix,iy,0)*sqrt(1.0_R8 - qc(ix,iy)**2)
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

  !> This routine performs an interpolation of cell-faced quantities to cell centres
  !! with special treatment for non-rectangular cells and attention to the direction
  !! of slanted faces.
  !! If the cell-faced quantity is a flux where the directional sign matters, we
  !! take it into account
  !! The cell-centered is bi-directional (0: poloidal, 1: radial)
  !! For a positive slant value, a positive poloidal flux in a slanted face is
  !! directed towards the top
  subroutine interp_from_face(isflux,isparallel,nx,ny,flux,centre)
  use b2mod_geo , only: crx, cry, gs, qz, qc, qcb, pbs, vol
  use b2mod_math
  use b2mod_indirect

  implicit none
  logical, intent(in) :: isflux, isparallel
  integer, intent(in) :: nx, ny
  real (kind=R8), intent(in) :: flux(-1:nx,-1:ny,0:1)
  real (kind=R8), intent(out) :: centre(-1:nx,-1:ny,0:1)

  integer i, ix, iy
  integer cgeo
  real (kind=R8) :: area_to_top, area_to_bottom, area_to_left, area_to_right
  real (kind=R8) :: face(-1:nx,-1:ny,0:1,0:1), slant
  logical classical
  real (kind=R8) :: weight(-1:nx,-1:ny,TO_SELF:TO_TOP)

  classical = isClassicalGrid(cflags)
!! If isparallel is .true., then the background flow is in the parallel direction
!! face has dimensions (ix,iy,iFace,iDir)
  do i = TO_SELF, TO_TOP
    weight(:,:,i) = vol(:,:)
  end do
  face = 0.0_R8
  do ix = -1, nx
    do iy = -1, ny
      if (isUnusedCell(cflags(ix,iy,CELLFLAG_TYPE))) cycle
      cgeo = cellGeoType( crx(ix,iy,:), cry(ix,iy,:) )
      if (classical .or. cflags(ix,iy,CELLFLAG_LEFTFACE) == GRID_UNDEFINED) then
        slant = 0.0_R8
      elseif (cflags(ix,iy,CELLFLAG_BOTTOMFACE) == GRID_UNDEFINED .and. &
            &  cgeo /= CGEO_TRIA_NOBOT .and. cgeo /= CGEO_TRIA_NOTOP .and. &
            &  cflags(ix,iy,CELLFLAG_TOPFACE) == GRID_UNDEFINED ) then
        slant = sinang ( crx(ix,iy,2) - crx(ix,iy,0), cry(ix,iy,2) - cry(ix,iy,0), &
                      & (crx(ix,iy,1)+crx(ix,iy,3)) - (crx(ix,iy,0)+crx(ix,iy,2)), &
                      & (cry(ix,iy,1)+cry(ix,iy,3)) - (cry(ix,iy,0)+cry(ix,iy,2)))
      elseif ((cgeo == CGEO_TRIA_NOBOT .or. &
            &  cflags(ix,iy,CELLFLAG_BOTTOMFACE) /= GRID_UNDEFINED) .and. &
            &  cgeo /= CGEO_TRIA_NOTOP .and. &
            &  cflags(ix,iy,CELLFLAG_TOPFACE) == GRID_UNDEFINED ) then
        slant = sinang ( crx(ix,iy,2) - crx(ix,iy,0), cry(ix,iy,2) - cry(ix,iy,0), &
                       & crx(ix,iy,3) - crx(ix,iy,2), cry(ix,iy,3) - cry(ix,iy,2))
      elseif (cflags(ix,iy,CELLFLAG_BOTTOMFACE) == GRID_UNDEFINED .and. &
           &  cgeo /= CGEO_TRIA_NOBOT .and. (cgeo == CGEO_TRIA_NOTOP .or. &
           &  cflags(ix,iy,CELLFLAG_TOPFACE) /= GRID_UNDEFINED) ) then
        slant = sinang ( crx(ix,iy,2) - crx(ix,iy,0), cry(ix,iy,2) - cry(ix,iy,0), &
                       & crx(ix,iy,1) - crx(ix,iy,0), cry(ix,iy,1) - cry(ix,iy,0))
      elseif ((cgeo == CGEO_TRIA_NOBOT .or. &
            &  cflags(ix,iy,CELLFLAG_BOTTOMFACE) /= GRID_UNDEFINED) .and. &
            &  cgeo /= CGEO_TRIA_NORIGHT .and. &
            & (cgeo == CGEO_TRIA_NOTOP .or. &
            &  cflags(ix,iy,CELLFLAG_TOPFACE) /= GRID_UNDEFINED) ) then
        slant = sinang ( crx(ix,iy,2) - crx(ix,iy,0), cry(ix,iy,2) - cry(ix,iy,0), &
                       & cry(ix,iy,3) - cry(ix,iy,1), crx(ix,iy,1) - crx(ix,iy,3))
      endif
      if (classical .or. &
       & (cflags(ix,iy,CELLFLAG_LEFTFACE) == GRID_UNDEFINED .and. cgeo /= CGEO_TRIA_NOLEFT)) then
        face(ix,iy,0,0) = flux(ix,iy,0)
      else if (cflags(ix,iy,CELLFLAG_LEFTFACE) /= GRID_UNDEFINED .and. cgeo /= CGEO_TRIA_NOLEFT) then
        if (isparallel) then
          face(ix,iy,0,0) = flux(ix,iy,0)
        else
          face(ix,iy,0,0) = flux(ix,iy,0)*sqrt(1.0_R8-slant**2)
          if (slant.ge.0.0_R8 .or. .not.isflux) then
            face(ix,iy,0,1) = flux(ix,iy,0)*(1.0_R8-sqrt(1.0_R8-slant**2))
          else
            face(ix,iy,0,1) =-flux(ix,iy,0)*(1.0_R8-sqrt(1.0_R8-slant**2))
          endif
        endif
      endif
      if (classical .or. &
       & (cflags(ix,iy,CELLFLAG_BOTTOMFACE) == GRID_UNDEFINED .and. cgeo /= CGEO_TRIA_NOBOT)) then
        face(ix,iy,1,1) = flux(ix,iy,1)
      else if (cflags(ix,iy,CELLFLAG_BOTTOMFACE) /= GRID_UNDEFINED .and. cgeo /= CGEO_TRIA_NOBOT) then
        if (isparallel) then
          if (pbs(ix,iy,1).ge.0.0_R8 .or. .not.isflux) then
            face(ix,iy,1,0) = flux(ix,iy,1)
          else
            face(ix,iy,1,0) =-flux(ix,iy,1)
          endif
        else
          if (pbs(ix,iy,1).ge.0.0_R8 .or. .not.isflux) then
            face(ix,iy,1,0) = flux(ix,iy,1)*(1.0_R8-sqrt(1.0_R8-qcb(ix,iy)**2))
          else
            face(ix,iy,1,0) =-flux(ix,iy,1)*(1.0_R8-sqrt(1.0_R8-qcb(ix,iy)**2))
          endif
          face(ix,iy,1,1) = flux(ix,iy,1)*sqrt(1.0_R8-qcb(ix,iy)**2)
        endif
      endif
    enddo
  enddo

  centre = 0.0_R8

  do ix = -1, nx
    do iy = -1, ny
      if(isUnusedCell(cflags(ix,iy,CELLFLAG_TYPE))) cycle
!! We treat the guard cells at the end so they can inherit real cell value when
!! necessary
      if(isGhostCell(cflags(ix,iy,CELLFLAG_TYPE))) cycle
      cgeo = cellGeoType( crx(ix,iy,:), cry(ix,iy,:) )
      if (isInDomain(nx,ny,topix(ix,iy),topiy(ix,iy))) then
        area_to_top = gs(topix(ix,iy),topiy(ix,iy),1) * &
                    & sqrt(1.0_R8 - qcb(topix(ix,iy),topiy(ix,iy))**2)
      else
        area_to_top = 0.0_R8
      endif
      if (isInDomain(nx,ny,bottomix(ix,iy),bottomiy(ix,iy))) then
        area_to_bottom = gs(ix,iy,1) * sqrt(1.0_R8 - qcb(ix,iy)**2)
      else
        area_to_bottom = 0.0_R8
      endif
      if (isInDomain(nx,ny,leftix(ix,iy),leftiy(ix,iy))) then
        area_to_left = gs(ix,iy,0) * qc(ix,iy)
      else
        area_to_left = 0.0_R8
      endif
      if (isInDomain(nx,ny,rightix(ix,iy),rightiy(ix,iy))) then
        area_to_right = gs(rightix(ix,iy),rightiy(ix,iy),0) * &
                      & qc(rightix(ix,iy),rightiy(ix,iy))
      else
        area_to_right = 0.0_R8
      endif
      if (cgeo == CGEO_TRIA_NOLEFT) then
        if (cflags(ix,iy,CELLFLAG_TOPFACE) /= GRID_UNDEFINED .and. &
         &  cflags(ix,iy,CELLFLAG_BOTTOMFACE) /= GRID_UNDEFINED) then
          centre(ix,iy,0) = 0.5_R8*( face(rightix(ix,iy),rightiy(ix,iy),0,0) + &
               & ( face(ix,iy,1,0)*area_to_bottom + &
               &   face(topix(ix,iy),topiy(ix,iy),1,0)*area_to_top ) / &
               & ( area_to_bottom + area_to_top) )
        elseif (cflags(ix,iy,CELLFLAG_TOPFACE) == GRID_UNDEFINED .and. &
               &  cflags(ix,iy,CELLFLAG_BOTTOMFACE) /= GRID_UNDEFINED) then
          centre(ix,iy,0) = 0.5_R8*( face(rightix(ix,iy),rightiy(ix,iy),0,0) + &
               &   face(ix,iy,1,0) )
        elseif (cflags(ix,iy,CELLFLAG_TOPFACE) /= GRID_UNDEFINED .and. &
               &  cflags(ix,iy,CELLFLAG_BOTTOMFACE) == GRID_UNDEFINED) then
          centre(ix,iy,0) = 0.5_R8*( face(rightix(ix,iy),rightiy(ix,iy),0,0) + &
               &                     face(topix(ix,iy),topiy(ix,iy),1,0) )
        endif
      else if (cgeo == CGEO_TRIA_NORIGHT) then
        if (cflags(ix,iy,CELLFLAG_TOPFACE) /= GRID_UNDEFINED .and. &
         &  cflags(ix,iy,CELLFLAG_BOTTOMFACE) /= GRID_UNDEFINED) then
          centre(ix,iy,0) = 0.5_R8*( face(ix,iy,0,0) + &
               & ( face(ix,iy,1,0)*area_to_bottom + &
               &   face(topix(ix,iy),topiy(ix,iy),1,0)*area_to_top ) / &
               & ( area_to_bottom + area_to_top) )
        elseif (cflags(ix,iy,CELLFLAG_TOPFACE) == GRID_UNDEFINED .and. &
             &  cflags(ix,iy,CELLFLAG_BOTTOMFACE) /= GRID_UNDEFINED) then
          centre(ix,iy,0) = 0.5_R8*( face(ix,iy,0,0) + face(ix,iy,1,0) )
        elseif (cflags(ix,iy,CELLFLAG_TOPFACE) /= GRID_UNDEFINED .and. &
             &  cflags(ix,iy,CELLFLAG_BOTTOMFACE) == GRID_UNDEFINED) then
          centre(ix,iy,0) = 0.5_R8*( face(ix,iy,0,0) + &
               &   face(topix(ix,iy),topiy(ix,iy),1,0) )
        endif
      else if (cgeo == CGEO_TRIA_NOTOP) then
        if (cflags(ix,iy,CELLFLAG_BOTTOMFACE) /= GRID_UNDEFINED) then
          if (cflags(ix,iy,CELLFLAG_LEFTFACE) /= GRID_UNDEFINED) then
            centre(ix,iy,0) = 0.5_R8*( face(rightix(ix,iy),rightiy(ix,iy),0,0) + &
               & ( face(ix,iy,1,0)*area_to_bottom + face(ix,iy,0,0)*area_to_left ) / &
               & ( area_to_bottom + area_to_left) )
          else if (cflags(ix,iy,CELLFLAG_RIGHTFACE) /= GRID_UNDEFINED) then
            centre(ix,iy,0) = 0.5_R8*( face(ix,iy,0,0) + &
               & ( face(ix,iy,1,0)*area_to_bottom + &
               &   face(rightix(ix,iy),rightiy(ix,iy),0,0)*area_to_right ) / &
               & ( area_to_bottom + area_to_right) )
          else
            centre(ix,iy,0) = 0.5_R8*( face(ix,iy,0,0) + &
               &   face(rightix(ix,iy),rightiy(ix,iy),0,0) )
          endif
        else
          centre(ix,iy,0) = 0.5_R8*( face(ix,iy,0,0) + face(rightix(ix,iy),rightiy(ix,iy),0,0) )
        endif
      else if (cgeo == CGEO_TRIA_NOBOT) then
        if (cflags(ix,iy,CELLFLAG_TOPFACE) /= GRID_UNDEFINED) then
          if (cflags(ix,iy,CELLFLAG_LEFTFACE) /= GRID_UNDEFINED) then
            centre(ix,iy,0) = 0.5_R8*( face(rightix(ix,iy),rightiy(ix,iy),0,0) + &
               & ( face(topix(ix,iy),topiy(ix,iy),1,0)*area_to_top + face(ix,iy,0,0)*area_to_left ) / &
               & ( area_to_top + area_to_left) )
          else if (cflags(ix,iy,CELLFLAG_RIGHTFACE) /= GRID_UNDEFINED) then
            centre(ix,iy,0) = 0.5_R8*( face(ix,iy,0,0) + &
               & ( face(topix(ix,iy),topiy(ix,iy),1,0)*area_to_top + &
               &   face(rightix(ix,iy),rightiy(ix,iy),0,0)*area_to_right ) / &
               & ( area_to_top + area_to_right) )
          else
            centre(ix,iy,0) = 0.5_R8*( face(ix,iy,0,0) + face(rightix(ix,iy),rightiy(ix,iy),0,0) )
          endif
        else
          centre(ix,iy,0) = 0.5_R8*( face(ix,iy,0,0) + face(rightix(ix,iy),rightiy(ix,iy),0,0) )
        endif
      else if (cgeo == CGEO_QUAD) then
        if (classical .or. (cflags(ix,iy,CELLFLAG_TOPFACE) == GRID_UNDEFINED .and. &
                          & cflags(ix,iy,CELLFLAG_BOTTOMFACE) == GRID_UNDEFINED) ) then
          centre(ix,iy,0) = 0.5_R8*( face(ix,iy,0,0) + face(rightix(ix,iy),rightiy(ix,iy),0,0) )
        elseif (cflags(ix,iy,CELLFLAG_TOPFACE) /= GRID_UNDEFINED .and. &
              & cflags(ix,iy,CELLFLAG_BOTTOMFACE) /= GRID_UNDEFINED) then
 !! both top and bottom sides are misaligned
          if (pbs(ix,iy,1).gt.0.0_R8 .and. & !! the bottom edge is being shaved going to the left
            & pbs(topix(ix,iy),topiy(ix,iy),1).gt.0.0_R8) then  !! the top edge is being shaved going to the right
            centre(ix,iy,0) = 0.5_R8* ( &
                & (face(ix,iy,0,0)*area_to_left + face(ix,iy,1,0)*area_to_bottom) / &
                & (area_to_left + area_to_bottom) + &
                & (face(topix(ix,iy),topiy(ix,iy),1,0)*area_to_top + &
                &  face(rightix(ix,iy),rightiy(ix,iy),0,0)*area_to_right) / &
                & (area_to_top + area_to_right) )
          elseif (pbs(ix,iy,1).lt.0.0_R8 .and. & !! the bottom edge is being shaved going to the right
                & pbs(topix(ix,iy),topiy(ix,iy),1).gt.0.0_R8) then  !! the top edge is being shaved going to the right
            centre(ix,iy,0) = 0.5_R8* ( face(ix,iy,0,0) + &
                & (face(topix(ix,iy),topiy(ix,iy),1,0)*area_to_top + &
                &  face(rightix(ix,iy),rightiy(ix,iy),0,0)*area_to_right + &
                &  face(ix,iy,1,0)*area_to_bottom) / &
                & (area_to_top + area_to_right + area_to_bottom) )
          elseif (pbs(ix,iy,1).lt.0.0_R8 .and. & !! the bottom edge is being shaved going to the right
                & pbs(topix(ix,iy),topiy(ix,iy),1).lt.0.0_R8) then  !! the top edge is being shaved going to the left
            centre(ix,iy,0) = 0.5_R8* ( &
                & (face(ix,iy,0,0)*area_to_left + &
                &  face(topix(ix,iy),topiy(ix,iy),1,0)*area_to_top) / &
                & (area_to_left + area_to_top) + &
                & (face(ix,iy,1,0)*area_to_bottom + &
                &  face(rightix(ix,iy),rightiy(ix,iy),0,0)*area_to_right) / &
                & (area_to_bottom + area_to_right) )
          elseif (pbs(ix,iy,1).gt.0.0_R8 .and. & !! the bottom edge is being shaved going to the left
                & pbs(topix(ix,iy),topiy(ix,iy),1).lt.0.0_R8) then  !! the top edge is being shaved going to the left
            centre(ix,iy,0) = 0.5_R8* ( face(rightix(ix,iy),rightiy(ix,iy),0,0) + &
                & (face(topix(ix,iy),topiy(ix,iy),1,0)*area_to_top + &
                &  face(ix,iy,0,0)*area_to_left + &
                &  face(ix,iy,1,0)*area_to_bottom) / &
                & (area_to_top + area_to_left + area_to_bottom) )
          elseif (pbs(ix,iy,1).eq.0.0_R8 .and. & !! the bottom edge is approximately aligned
                & pbs(topix(ix,iy),topiy(ix,iy),1).lt.0.0_R8) then  !! the top edge is being shaved going to the left
            centre(ix,iy,0) = 0.5_R8* ( face(rightix(ix,iy),rightiy(ix,iy),0,0) + &
                & (face(topix(ix,iy),topiy(ix,iy),1,0)*area_to_top + &
                &  face(ix,iy,0,0)*area_to_left) / &
                & (area_to_top + area_to_left) )
          elseif (pbs(ix,iy,1).eq.0.0_R8 .and. & !! the bottom edge is approximately aligned
                & pbs(topix(ix,iy),topiy(ix,iy),1).gt.0.0_R8) then  !! the top edge is being shaved going to the right
            centre(ix,iy,0) = 0.5_R8* ( face(ix,iy,0,0) + &
                & (face(topix(ix,iy),topiy(ix,iy),1,0)*area_to_top + &
                &  face(rightix(ix,iy),rightiy(ix,iy),0,0)*area_to_right) / &
                & (area_to_top + area_to_right) )
          elseif (pbs(ix,iy,1).lt.0.0_R8 .and. & !! the bottom edge is being shaved going to the right
                & pbs(topix(ix,iy),topiy(ix,iy),1).eq.0.0_R8) then  !! the top edge is approximately aligned
            centre(ix,iy,0) = 0.5_R8* ( face(ix,iy,0,0) + &
                & (face(rightix(ix,iy),rightiy(ix,iy),0,0)*area_to_right + &
                &  face(ix,iy,1,0)*area_to_bottom) / &
                & (area_to_right + area_to_bottom) )
          elseif (pbs(ix,iy,1).gt.0.0_R8 .and. & !! the bottom edge is being shaved going to the left
                & pbs(topix(ix,iy),topiy(ix,iy),1).eq.0.0_R8) then  !! the top edge is approximately aligned
            centre(ix,iy,0) = 0.5_R8* ( &
                & (face(ix,iy,0,0)*area_to_left + face(ix,iy,1,0)*area_to_bottom) / &
                & (area_to_left + area_to_bottom) + &
                &  face(rightix(ix,iy),rightiy(ix,iy),0,0) )
          elseif (pbs(ix,iy,1).eq.0.0_R8 .and. &
                & pbs(topix(ix,iy),topiy(ix,iy),1).eq.0.0_R8) then  !! both edges are approximately aligned
            centre(ix,iy,0) = 0.5_R8*( face(ix,iy,0,0) + face(rightix(ix,iy),rightiy(ix,iy),0,0) )
          endif
        elseif (cflags(ix,iy,CELLFLAG_TOPFACE) == GRID_UNDEFINED .and. &
              & cflags(ix,iy,CELLFLAG_BOTTOMFACE) /= GRID_UNDEFINED) then
          if (pbs(ix,iy,1).lt.0.0_R8) then  !! the bottom edge is being shaved going to the right
            centre(ix,iy,0) = 0.5_R8* ( face(ix,iy,0,0) + &
                & (face(rightix(ix,iy),rightiy(ix,iy),0,0)*area_to_right + &
                &  face(ix,iy,1,0)*area_to_bottom) / &
                & (area_to_right + area_to_bottom) )
          elseif (pbs(ix,iy,1).gt.0.0_R8) then  !! the bottom edge is being shaved going to the left
            centre(ix,iy,0) = 0.5_R8* ( &
                & (face(ix,iy,0,0)*area_to_left + face(ix,iy,1,0)*area_to_bottom) / &
                & (area_to_left + area_to_bottom) + &
                &  face(rightix(ix,iy),rightiy(ix,iy),0,0) )
          elseif (pbs(ix,iy,1).eq.0.0_R8) then  !! the bottom edge is approximately aligned
            centre(ix,iy,0) = 0.5_R8*( face(ix,iy,0,0) + face(rightix(ix,iy),rightiy(ix,iy),0,0) )
          endif
        elseif (cflags(ix,iy,CELLFLAG_TOPFACE) /= GRID_UNDEFINED .and. &
              & cflags(ix,iy,CELLFLAG_BOTTOMFACE) == GRID_UNDEFINED) then
          if (pbs(topix(ix,iy),topiy(ix,iy),1).lt.0.0_R8) then  !! the top edge is being shaved going to the left
            centre(ix,iy,0) = 0.5_R8* ( face(rightix(ix,iy),rightiy(ix,iy),0,0) + &
                & (face(topix(ix,iy),topiy(ix,iy),1,0)*area_to_top + &
                &  face(ix,iy,0,0)*area_to_left) / &
                & (area_to_top + area_to_left) )
          elseif (pbs(topix(ix,iy),topiy(ix,iy),1).gt.0.0_R8) then  !! the top edge is being shaved going to the right
            centre(ix,iy,0) = 0.5_R8* ( face(ix,iy,0,0) + &
                & (face(topix(ix,iy),topiy(ix,iy),1,0)*area_to_top + &
                &  face(rightix(ix,iy),rightiy(ix,iy),0,0)*area_to_right) / &
                & (area_to_top + area_to_right) )
          elseif (pbs(topix(ix,iy),topiy(ix,iy),1).eq.0.0_R8) then  !! the top edge is approximately aligned
            centre(ix,iy,0) = 0.5_R8*( face(ix,iy,0,0) + face(rightix(ix,iy),rightiy(ix,iy),0,0) )
          endif
        endif
      endif
      if (isInDomain(nx,ny,topix(ix,iy),topiy(ix,iy))) then
        area_to_top = gs(topix(ix,iy),topiy(ix,iy),1) * qcb(topix(ix,iy),topiy(ix,iy))
      else
        area_to_top = 0.0_R8
      endif
      if (isInDomain(nx,ny,bottomix(ix,iy),bottomiy(ix,iy))) then
        area_to_bottom = gs(ix,iy,1) * qcb(ix,iy)
      else
        area_to_bottom = 0.0_R8
      endif
      if (isInDomain(nx,ny,leftix(ix,iy),leftiy(ix,iy))) then
        area_to_left = gs(ix,iy,0) * sqrt(1.0_R8 - qc(ix,iy)**2 )
      else
        area_to_left = 0.0_R8
      endif
      if (isInDomain(nx,ny,rightix(ix,iy),rightiy(ix,iy))) then
        area_to_right = gs(rightix(ix,iy),rightiy(ix,iy),0) * &
                      & sqrt(1.0_R8 - qc(rightix(ix,iy),rightiy(ix,iy))**2 )
      else
        area_to_right = 0.0_R8
      endif
      if (cgeo == CGEO_TRIA_NOBOT) then
        if (cflags(ix,iy,CELLFLAG_LEFTFACE) /= GRID_UNDEFINED .and. &
         &  cflags(ix,iy,CELLFLAG_RIGHTFACE) /= GRID_UNDEFINED) then
          centre(ix,iy,1) = 0.5_R8 * ( face(topix(ix,iy),topiy(ix,iy),1,1) + &
                 & (face(ix,iy,0,1) * area_to_left + &
                 &  face(rightix(ix,iy),rightiy(ix,iy),0,1) * area_to_right) / &
                 & (area_to_left + area_to_right) )
        elseif (cflags(ix,iy,CELLFLAG_LEFTFACE) == GRID_UNDEFINED .and. &
             &  cflags(ix,iy,CELLFLAG_RIGHTFACE) /= GRID_UNDEFINED) then
          centre(ix,iy,1) = 0.5_R8 * ( face(topix(ix,iy),topiy(ix,iy),1,1) + &
                 &  face(rightix(ix,iy),rightiy(ix,iy),0,1) )
        elseif (cflags(ix,iy,CELLFLAG_LEFTFACE) /= GRID_UNDEFINED .and. &
         &      cflags(ix,iy,CELLFLAG_RIGHTFACE) == GRID_UNDEFINED) then
          centre(ix,iy,1) = 0.5_R8 * ( face(topix(ix,iy),topiy(ix,iy),1,1) + &
                 &  face(ix,iy,0,1) )
        endif
      else if (cgeo == CGEO_TRIA_NOTOP) then
        if (cflags(ix,iy,CELLFLAG_LEFTFACE) /= GRID_UNDEFINED .and. &
         &  cflags(ix,iy,CELLFLAG_RIGHTFACE) /= GRID_UNDEFINED) then
          centre(ix,iy,1) = 0.5_R8 * ( face(ix,iy,1,1) + &
                 & (face(ix,iy,0,1) * area_to_left + &
                 &  face(rightix(ix,iy),rightiy(ix,iy),0,1) * area_to_right) / &
                 & (area_to_left + area_to_right) )
        elseif (cflags(ix,iy,CELLFLAG_LEFTFACE) == GRID_UNDEFINED .and. &
             &  cflags(ix,iy,CELLFLAG_RIGHTFACE) /= GRID_UNDEFINED) then
          centre(ix,iy,1) = 0.5_R8 * ( face(ix,iy,1,1) + &
                 &  face(rightix(ix,iy),rightiy(ix,iy),0,1) )
        elseif (cflags(ix,iy,CELLFLAG_LEFTFACE) /= GRID_UNDEFINED .and. &
             &  cflags(ix,iy,CELLFLAG_RIGHTFACE) == GRID_UNDEFINED) then
          centre(ix,iy,1) = 0.5_R8 * ( face(ix,iy,1,1) + face(ix,iy,0,1) )
        endif
      else if (cgeo == CGEO_TRIA_NOLEFT) then
        if (cflags(ix,iy,CELLFLAG_RIGHTFACE) /= GRID_UNDEFINED) then
          if (cflags(ix,iy,CELLFLAG_TOPFACE) /= GRID_UNDEFINED) then
            centre(ix,iy,1) = 0.5_R8 * ( face(ix,iy,1,1) + &
                  & ( face(topix(ix,iy),topiy(ix,iy),1,1)*area_to_top + &
                  &   face(rightix(ix,iy),rightiy(ix,iy),0,1)*area_to_right ) / &
                  & ( area_to_top + area_to_right ) )
          elseif (cflags(ix,iy,CELLFLAG_BOTTOMFACE) /= GRID_UNDEFINED) then
            centre(ix,iy,1) = 0.5_R8 * ( face(topix(ix,iy),topiy(ix,iy),1,1) + &
                  & ( face(ix,iy,1,1)*area_to_bottom + &
                  &   face(rightix(ix,iy),rightiy(ix,iy),0,1)*area_to_right ) / &
                  & ( area_to_bottom + area_to_right ) )
          else
            centre(ix,iy,1) = 0.5_R8 * ( face(ix,iy,1,1) + &
                  &   face(topix(ix,iy),topiy(ix,iy),1,1) )
          endif
        else
          centre(ix,iy,1) = 0.5_R8 * ( face(ix,iy,1,1) + &
                &   face(topix(ix,iy),topiy(ix,iy),1,1) )
        endif
      else if (cgeo == CGEO_TRIA_NORIGHT) then
        if (cflags(ix,iy,CELLFLAG_LEFTFACE) /= GRID_UNDEFINED) then
          if (cflags(ix,iy,CELLFLAG_TOPFACE) /= GRID_UNDEFINED) then
            centre(ix,iy,1) = 0.5_R8 * ( face(ix,iy,1,1) + &
                  & ( face(topix(ix,iy),topiy(ix,iy),1,1)*area_to_top + &
                  &   face(ix,iy,0,1)*area_to_left ) / &
                  & ( area_to_top + area_to_left ) )
          elseif (cflags(ix,iy,CELLFLAG_BOTTOMFACE) /= GRID_UNDEFINED) then
            centre(ix,iy,1) = 0.5_R8 * ( face(topix(ix,iy),topiy(ix,iy),1,1) + &
                  & ( face(ix,iy,1,1)*area_to_bottom + &
                  &   face(ix,iy,0,1)*area_to_left ) / &
                  & ( area_to_bottom + area_to_left ) )
          else
            centre(ix,iy,1) = 0.5_R8 * ( face(ix,iy,1,1) + &
                  &   face(topix(ix,iy),topiy(ix,iy),1,1) )
          endif
        else
          centre(ix,iy,1) = 0.5_R8 * ( face(ix,iy,1,1) + &
                  &   face(topix(ix,iy),topiy(ix,iy),1,1) )
        endif
      else if (cgeo == CGEO_QUAD) then
        if (classical .or. (cflags(ix,iy,CELLFLAG_LEFTFACE) == GRID_UNDEFINED .and. &
                          & cflags(ix,iy,CELLFLAG_RIGHTFACE) == GRID_UNDEFINED) ) then
          centre(ix,iy,1) = 0.5_R8*( face(ix,iy,1,1) + face(topix(ix,iy),topiy(ix,iy),1,1) )
        elseif (cflags(ix,iy,CELLFLAG_LEFTFACE) /= GRID_UNDEFINED .and. &
              & cflags(ix,iy,CELLFLAG_RIGHTFACE) /= GRID_UNDEFINED) then
 !! both left and right sides are misaligned
          if (weight(ix,iy,TO_TOP).gt.weight(ix,iy,TO_BOTTOM) .and. &
           &  weight(ix,iy,TO_TOP).eq.weight(ix,iy,TO_SELF)) then  !! the cell is keystone-shaped
            centre(ix,iy,1) = 0.5_R8 * ( face(topix(ix,iy),topiy(ix,iy),1,1) + &
                   & ( face(ix,iy,1,1)*area_to_bottom + &
                   &   face(ix,iy,0,1)*area_to_left + &
                   &   face(rightix(ix,iy),rightiy(ix,iy),0,1)*area_to_right ) / &
                   & ( area_to_bottom + area_to_left + area_to_right ) )
          else if (weight(ix,iy,TO_TOP).lt.weight(ix,iy,TO_BOTTOM) .and. &
              &    weight(ix,iy,TO_BOTTOM).eq.weight(ix,iy,TO_SELF)) then  !! the cell is an inverted keystone
            centre(ix,iy,1) = 0.5_R8 * ( face(ix,iy,1,1) + &
                   & ( face(topix(ix,iy),topiy(ix,iy),1,1)*area_to_top + &
                   &   face(ix,iy,0,1)*area_to_left + &
                   &   face(rightix(ix,iy),rightiy(ix,iy),0,1)*area_to_right ) / &
                   & ( area_to_top + area_to_left + area_to_right ) )
          else if (weight(ix,iy,TO_TOP).lt.weight(ix,iy,TO_SELF).and. &
                &  weight(ix,iy,TO_BOTTOM).lt.weight(ix,iy,TO_SELF)) then
 !! left and right edges are being shaved in the same direction
            if (qz(ix,iy,0).lt.0.0_R8) then  !! the cell is slanted to the right
              centre(ix,iy,1) = 0.5_R8 * ( &
                     & ( face(ix,iy,1,1)*area_to_bottom + &
                     &   face(rightix(ix,iy),rightiy(ix,iy),0,1)*area_to_right ) / &
                     & ( area_to_bottom + area_to_right ) + &
                     & ( face(ix,iy,0,1)*area_to_left + &
                     &   face(topix(ix,iy),topiy(ix,iy),1,1)*area_to_top ) / &
                     & ( area_to_top + area_to_left ) )
            elseif (qz(ix,iy,0).gt.0.0_R8) then  !! the cell is slanted to the left
              centre(ix,iy,1) = 0.5_R8 * ( &
                     & ( face(ix,iy,1,1)*area_to_bottom + &
                     &   face(ix,iy,0,1)*area_to_left ) / &
                     & ( area_to_bottom + area_to_left ) + &
                     & ( face(rightix(ix,iy),rightiy(ix,iy),0,1)*area_to_right + &
                     &   face(topix(ix,iy),topiy(ix,iy),1,1)*area_to_top ) / &
                     & ( area_to_top + area_to_right ) )
            elseif (qz(ix,iy,0).eq.0.0_R8) then  !! the cell is approximately aligned
              centre(ix,iy,1) = 0.5_R8 * ( face(ix,iy,1,1) + &
                     &   face(topix(ix,iy),topiy(ix,iy),1,1) )
            endif
          else !! both edges are approximately aligned
            centre(ix,iy,1) = 0.5_R8 * ( face(ix,iy,1,1) + &
                     &   face(topix(ix,iy),topiy(ix,iy),1,1) )
          endif
        elseif (cflags(ix,iy,CELLFLAG_LEFTFACE) /= GRID_UNDEFINED .and. &
              & cflags(ix,iy,CELLFLAG_RIGHTFACE) == GRID_UNDEFINED) then
!! the left side is misaligned
          if (weight(ix,iy,TO_TOP).lt.weight(ix,iy,TO_BOTTOM)) then  !! the cell is being shaved going up
            centre(ix,iy,1) = 0.5_R8 * ( face(ix,iy,1,1) + &
                     & ( face(ix,iy,0,1)*area_to_left + &
                     &   face(topix(ix,iy),topiy(ix,iy),1,1)*area_to_top ) / &
                     & ( area_to_top + area_to_left ) )
          else if (weight(ix,iy,TO_TOP).gt.weight(ix,iy,TO_BOTTOM)) then  !! the cell is being shaved going down
            centre(ix,iy,1) = 0.5_R8 * ( &
                     & ( face(ix,iy,1,1)*area_to_bottom + &
                     &   face(ix,iy,0,1)*area_to_left ) / &
                     & ( area_to_bottom + area_to_left ) + &
                     &   face(topix(ix,iy),topiy(ix,iy),1,1) )
          else !! case of approximately rectangular boundary cell
            centre(ix,iy,1) = 0.5_R8 * ( face(ix,iy,1,1) + &
                     &   face(topix(ix,iy),topiy(ix,iy),1,1) )
          endif
        elseif (cflags(ix,iy,CELLFLAG_LEFTFACE) == GRID_UNDEFINED .and. &
              & cflags(ix,iy,CELLFLAG_RIGHTFACE) /= GRID_UNDEFINED) then
 !! the right side is misaligned
          if (weight(ix,iy,TO_TOP).lt.weight(ix,iy,TO_BOTTOM)) then  !! the cell is being shaved going up
            centre(ix,iy,1) = 0.5_R8 * ( face(ix,iy,1,1) + &
                     & ( face(rightix(ix,iy),rightiy(ix,iy),0,1)*area_to_right + &
                     &   face(topix(ix,iy),topiy(ix,iy),1,1)*area_to_top ) / &
                     & ( area_to_top + area_to_right ) )
          else if (weight(ix,iy,TO_TOP).gt.weight(ix,iy,TO_BOTTOM)) then  !! the cell is being shaved going down
            centre(ix,iy,1) = 0.5_R8 * ( &
                     & ( face(ix,iy,1,1)*area_to_bottom + &
                     &   face(rightix(ix,iy),rightiy(ix,iy),0,1)*area_to_right ) / &
                     & ( area_to_bottom + area_to_right ) + &
                     & face(topix(ix,iy),topiy(ix,iy),1,1)*area_to_top )
          else !! case of approximately rectangular boundary cell
            centre(ix,iy,1) = 0.5_R8 * ( face(ix,iy,1,1) + &
                     &   face(topix(ix,iy),topiy(ix,iy),1,1) )
          endif
        endif
      endif
    end do
  end do

!! We now treat the special case of guard cells
!! that inherit directly the value from the valid face when available
!! If not, and neighbours exist, we do a simple interpolation
!! If not again, we recopy the value from the neighbouring real cell
  do ix = -1, nx
    do iy = -1, ny
      if(.not.isGhostCell(cflags(ix,iy,CELLFLAG_TYPE))) cycle
      if(ghostGetBndFace(nx,ny,ix,iy) == RIGHT) then
        centre(ix,iy,0) = face(rightix(ix,iy),rightiy(ix,iy),0,0)
        if (isInDomain(nx,ny,bottomix(ix,iy),bottomiy(ix,iy)) .and. &
          & isInDomain(nx,ny,topix(ix,iy),topiy(ix,iy))) then
          centre(ix,iy,1) = 0.5_R8* ( face(ix,iy,1,1) + &
          & face(topix(ix,iy),topiy(ix,iy),1,1) )
        elseif (isInDomain(nx,ny,bottomix(ix,iy),bottomiy(ix,iy)) .and. &
         & .not.isInDomain(nx,ny,topix(ix,iy),topiy(ix,iy))) then
          centre(ix,iy,1) = face(ix,iy,1,1)
        elseif (.not.isInDomain(nx,ny,bottomix(ix,iy),bottomiy(ix,iy)) .and. &
          &          isInDomain(nx,ny,topix(ix,iy),topiy(ix,iy))) then
          centre(ix,iy,1) = face(topix(ix,iy),topiy(ix,iy),1,1)
        endif
      else if(ghostGetBndFace(nx,ny,ix,iy) == LEFT) then
        centre(ix,iy,0) = face(ix,iy,0,0)
        if (isInDomain(nx,ny,bottomix(ix,iy),bottomiy(ix,iy)) .and. &
          & isInDomain(nx,ny,topix(ix,iy),topiy(ix,iy))) then
          centre(ix,iy,1) = 0.5_R8* ( face(ix,iy,0,1) + &
          & face(topix(ix,iy),topiy(ix,iy),0,1) )
        elseif (isInDomain(nx,ny,bottomix(ix,iy),bottomiy(ix,iy)) .and. &
         & .not.isInDomain(nx,ny,topix(ix,iy),topiy(ix,iy))) then
          centre(ix,iy,1) = face(ix,iy,0,1)
        elseif (.not.isInDomain(nx,ny,bottomix(ix,iy),bottomiy(ix,iy)) .and. &
          &          isInDomain(nx,ny,topix(ix,iy),topiy(ix,iy))) then
          centre(ix,iy,1) = face(topix(ix,iy),topiy(ix,iy),0,1)
        endif
      else if(ghostGetBndFace(nx,ny,ix,iy) == BOTTOM) then
        if (isparallel) then
          centre(ix,iy,0) = face(ix,iy,1,0)
        else
          if (isInDomain(nx,ny,leftix(ix,iy),leftiy(ix,iy)) .and. &
            & isInDomain(nx,ny,rightix(ix,iy),rightiy(ix,iy))) then
            centre(ix,iy,0) = 0.5_R8* ( face(ix,iy,0,0) + &
            &   face(rightix(ix,iy),rightiy(ix,iy),0,0) )
          else if (.not.isInDomain(nx,ny,leftix(ix,iy),leftiy(ix,iy)) .and. &
            &           isInDomain(nx,ny,rightix(ix,iy),rightiy(ix,iy))) then
            centre(ix,iy,0) = face(rightix(ix,iy),rightiy(ix,iy),0,0)
          else if (isInDomain(nx,ny,leftix(ix,iy),leftiy(ix,iy)) .and. &
            & .not.isInDomain(nx,ny,rightix(ix,iy),rightiy(ix,iy))) then
            centre(ix,iy,0) = face(ix,iy,0,0)
          end if
          centre(ix,iy,1) = face(ix,iy,1,1)
        endif
      else if(ghostGetBndFace(nx,ny,ix,iy) == TOP) then
        if (isparallel) then
          centre(ix,iy,0) = face(topix(ix,iy),topiy(ix,iy),1,0)
        else
          if (isInDomain(nx,ny,leftix(ix,iy),leftiy(ix,iy)) .and. &
            & isInDomain(nx,ny,rightix(ix,iy),rightiy(ix,iy))) then
            centre(ix,iy,0) = 0.5_R8* ( face(ix,iy,0,0) + &
            &   face(rightix(ix,iy),rightiy(ix,iy),0,0) )
          else if (.not.isInDomain(nx,ny,leftix(ix,iy),leftiy(ix,iy)) .and. &
            &           isInDomain(nx,ny,rightix(ix,iy),rightiy(ix,iy))) then
            centre(ix,iy,0) = face(rightix(ix,iy),rightiy(ix,iy),0,0)
          else if (isInDomain(nx,ny,leftix(ix,iy),leftiy(ix,iy)) .and. &
            & .not.isInDomain(nx,ny,rightix(ix,iy),rightiy(ix,iy))) then
            centre(ix,iy,0) = face(ix,iy,0,0)
          end if
          centre(ix,iy,1) = face(ix,iy,1,1)
        end if
      end if
    end do
  end do

  return
  end subroutine interp_from_face

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
  use b2mod_geo
  use b2mod_indirect
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

  !> This subroutine computes a FACE-CENTERED velocity by dividing a flow
  !> by a projected area and a FACE-CENTERED density
  subroutine face_velocity_from_flow(nx,ny,ns,flow,density,velocity)
  use b2mod_geo
  use b2mod_indirect
  implicit none
!! input arguments
  integer, intent(in) :: nx, ny, ns
  real (kind=R8), intent(in) :: flow(-1:nx,-1:ny,0:1,0:ns-1), &
                           & density(-1:nx,-1:ny,0:1,0:ns-1)
!! output arguments
  real (kind=R8), intent(out) :: velocity(-1:nx,-1:ny,0:1,0:ns-1)
!! local variables
  integer ix, iy, is
  real (kind=R8) :: area_to_top, area_to_bottom

  velocity = 0.0_R8

  do is = 0, ns-1
    do iy = -1, ny
      do ix = -1, nx
        if (isUnusedCell(cflags(ix,iy,CELLFLAG_TYPE))) cycle
        if (density(ix,iy,0,is).eq.0.0_R8) cycle
        if (density(ix,iy,1,is).eq.0.0_R8) cycle
        if (isInDomain(nx,ny,leftix(ix,iy),leftiy(ix,iy))) then
          velocity(ix,iy,0,is) = flow(ix,iy,0,is)/ &
           & (gs(ix,iy,0)*qc(ix,iy))/density(ix,iy,0,is)
        elseif ((isBoundaryCell(cflags(ix,iy,CELLFLAG_TYPE)) .or. &
              &  isGhostCell(cflags(ix,iy,CELLFLAG_TYPE))) .and. &
              &  cflags(ix,iy,CELLFLAG_TOPFACE) /= GRID_UNDEFINED .and. &
              &  cflags(ix,iy,CELLFLAG_BOTTOMFACE) == GRID_UNDEFINED) then
          if (pbs(topix(ix,iy),topiy(ix,iy),1).ne.0.0_R8) &
          & velocity(ix,iy,0,is) = flow(topix(ix,iy),topiy(ix,iy),1,is)/ &
              &  (gs(topix(ix,iy),topiy(ix,iy),1)* &
              &   sqrt(1.0_R8-qcb(topix(ix,iy),topiy(ix,iy))**2))/ &
              &  density(topix(ix,iy),topiy(ix,iy),1,is)
        elseif ((isBoundaryCell(cflags(ix,iy,CELLFLAG_TYPE)) .or. &
              &  isGhostCell(cflags(ix,iy,CELLFLAG_TYPE))) .and. &
              &  cflags(ix,iy,CELLFLAG_BOTTOMFACE) /= GRID_UNDEFINED .and. &
              &  cflags(ix,iy,CELLFLAG_TOPFACE) == GRID_UNDEFINED) then
          if (pbs(ix,iy,1).ne.0.0_R8) &
          & velocity(ix,iy,0,is) = flow(ix,iy,1,is)/density(ix,iy,1,is)/ &
              &  (gs(ix,iy,1)*sqrt(1.0_R8-qcb(ix,iy)**2))
        elseif ((isBoundaryCell(cflags(ix,iy,CELLFLAG_TYPE)) .or. &
              &  isGhostCell(cflags(ix,iy,CELLFLAG_TYPE))) .and. &
              &  cflags(ix,iy,CELLFLAG_BOTTOMFACE) /= GRID_UNDEFINED .and. &
              &  cflags(ix,iy,CELLFLAG_TOPFACE) /= GRID_UNDEFINED) then
          area_to_top = gs(topix(ix,iy),topiy(ix,iy),1)* &
              &   sqrt(1.0_R8-qcb(topix(ix,iy),topiy(ix,iy))**2)
          area_to_bottom = gs(ix,iy,1)*sqrt(1.0_R8-qcb(ix,iy)**2)
          if (area_to_top.ne.0.0_R8) velocity(ix,iy,0,is) = &
              &   flow(topix(ix,iy),topiy(ix,iy),1,is)/area_to_top/ &
              &   density(topix(ix,iy),topiy(ix,iy),1,is)
          if (area_to_bottom.ne.0.0_R8) velocity(ix,iy,0,is) = &
              &   velocity(ix,iy,0,is) + &
              &   flow(ix,iy,1,is)/area_to_bottom/density(ix,iy,1,is)
        end if
        if (isInDomain(nx,ny,bottomix(ix,iy),bottomiy(ix,iy))) then
          if (pbs(ix,iy,1).eq.0.0_R8) then
            velocity(ix,iy,1,is) = flow(ix,iy,1,is)/gs(ix,iy,1)/ &
                              & density(ix,iy,1,is)
          else
            velocity(ix,iy,1,is) = flow(ix,iy,1,is)/(gs(ix,iy,1)*qcb(ix,iy))/ &
                              & density(ix,iy,1,is)
          endif
        elseif ((isBoundaryCell(cflags(ix,iy,CELLFLAG_TYPE)) .or. &
               & isGhostCell(cflags(ix,iy,CELLFLAG_TYPE))) .and. &
               & cflags(ix,iy,CELLFLAG_LEFTFACE) /= GRID_UNDEFINED .and. &
               & cflags(ix,iy,CELLFLAG_RIGHTFACE) == GRID_UNDEFINED) then
          if (qc(ix,iy).ne.1.0_R8) &
           & velocity(ix,iy,1,is) = flow(ix,iy,0,is)/density(ix,iy,0,is) &
               & /(gs(ix,iy,0)*sqrt(1.0_R8-qc(ix,iy)**2))
        elseif ((isBoundaryCell(cflags(ix,iy,CELLFLAG_TYPE)) .or. &
               & isGhostCell(cflags(ix,iy,CELLFLAG_TYPE))) .and. &
               & cflags(ix,iy,CELLFLAG_RIGHTFACE) /= GRID_UNDEFINED .and. &
               & cflags(ix,iy,CELLFLAG_LEFTFACE) == GRID_UNDEFINED) then
          if (qc(rightix(ix,iy),rightiy(ix,iy)).ne.1.0_R8) &
           & velocity(ix,iy,1,is) = &
               & flow(rightix(ix,iy),rightiy(ix,iy),0,is)/ &
               & density(rightix(ix,iy),rightiy(ix,iy),0,is)/ &
               & (gs(rightix(ix,iy),rightiy(ix,iy),0)* &
               & sqrt(1.0_R8-qc(rightix(ix,iy),rightiy(ix,iy))**2))
        elseif ((isBoundaryCell(cflags(ix,iy,CELLFLAG_TYPE)) .or. &
               & isGhostCell(cflags(ix,iy,CELLFLAG_TYPE))) .and. &
               & cflags(ix,iy,CELLFLAG_RIGHTFACE) /= GRID_UNDEFINED .and. &
               & cflags(ix,iy,CELLFLAG_LEFTFACE) /= GRID_UNDEFINED) then
          if (qc(ix,iy).ne.1.0_R8) &
           & velocity(ix,iy,1,is) = flow(ix,iy,0,is)/density(ix,iy,0,is) &
               & /(gs(ix,iy,0)*sqrt(1.0_R8-qc(ix,iy)**2))
          if (qc(rightix(ix,iy),rightiy(ix,iy)).ne.1.0_R8) &
           & velocity(ix,iy,1,is) = velocity(ix,iy,1,is) + &
               & flow(rightix(ix,iy),rightiy(ix,iy),0,is)/ &
               & density(rightix(ix,iy),rightiy(ix,iy),0,is)/ &
               & (gs(rightix(ix,iy),rightiy(ix,iy),0)* &
               & sqrt(1.0_R8-qc(rightix(ix,iy),rightiy(ix,iy))**2))
        endif
      end do
    end do
  end do

  return
  end subroutine face_velocity_from_flow

  !> This subroutine computes a CELL-CENTERED 2-component velocity by dividing
  !> FACE-CENTERED flows by a cross-sectional area and a CELL-CENTERED density
  !> If isparallel is .true., then the flow is purely poloidal
  subroutine cell_velocity_from_flow(nx,ny,ns,isparallel,flow,density,velocity)
  use b2mod_geo
  use b2mod_indirect
  implicit none
!! input arguments
  integer, intent(in) :: nx, ny, ns
  logical, intent(in) :: isparallel
  real (kind=R8), intent(in) :: flow(-1:nx,-1:ny,0:1,0:ns-1), &
                           & density(-1:nx,-1:ny,0:ns-1)
!! output arguments
  real (kind=R8), intent(out) :: velocity(-1:nx,-1:ny,0:1,0:ns-1)
!! local variables
  integer ix, iy, is
  real (kind=R8) :: flux(-1:nx,-1:ny,0:1,0:ns-1)

  do is = 0, ns-1
    do iy = -1, ny
      do ix = -1, nx
        if (isInDomain(nx,ny,leftix(ix,iy),leftiy(ix,iy))) then
          flux(ix,iy,0,is) = flow(ix,iy,0,is)/gs(ix,iy,0)/qc(ix,iy)
        else
          flux(ix,iy,0,is) = 0.0_R8
        end if
        if (isInDomain(nx,ny,bottomix(ix,iy),bottomiy(ix,iy))) then
          if (pbs(ix,iy,1).ne.0.0_R8) then
            flux(ix,iy,1,is) = flow(ix,iy,1,is)/gs(ix,iy,1)/qcb(ix,iy)
          else
            flux(ix,iy,1,is) = flow(ix,iy,1,is)/gs(ix,iy,1)
          end if
        else
          flux(ix,iy,1,is) = 0.0_R8
        end if
      end do
    end do
  end do

  velocity = 0.0_R8

  do is = 0, ns-1
    call interp_from_face(.true.,isparallel,nx,ny,flux(-1,-1,0,is),velocity(-1,-1,0,is))
    do iy = -1, ny
      do ix = -1, nx
        if (isUnusedCell(cflags(ix,iy,CELLFLAG_TYPE))) then
          velocity(ix,iy,0,is) = 0.0_R8
          velocity(ix,iy,1,is) = 0.0_R8
          cycle
        end if
        if (density(ix,iy,is).eq.0.0_R8) cycle
        velocity(ix,iy,0,is) = velocity(ix,iy,0,is)/density(ix,iy,is)
        velocity(ix,iy,1,is) = velocity(ix,iy,1,is)/density(ix,iy,is)
      end do
    end do
  end do

  return
  end subroutine cell_velocity_from_flow

  !> This subroutine computes a flow by multiplying CELL-CENTERED velocities and densities
  !> through a projected contact area
  !> The flow produced is a cell-faced quantity
  !> The flow is positive in the usual directions, thus the sign of pbs appears
  !> We interpolate internally to the cell faces
  subroutine flow_from_velocity(nx,ny,ns,velocity,density,flow)
  use b2mod_geo
  use b2mod_indirect
  implicit none
!! input arguments
  integer, intent(in) :: nx, ny, ns
  real (kind=R8), intent(in) :: velocity(-1:nx,-1:ny,0:1,0:ns-1), &
                              & density(-1:nx,-1:ny,0:ns-1)
!! output arguments
  real (kind=R8), intent(out) :: flow(-1:nx,-1:ny,0:1,0:ns-1)
!! local variables
  integer ix, iy, is, cgeo
  real (kind=R8) :: den(-1:nx,-1:ny,0:1,0:ns-1), vv(-1:nx,-1:ny,0:1,0:ns-1)
  real (kind=R8) :: weight(-1:nx,-1:ny,TO_SELF:TO_TOP)
!! procedures
  intrinsic sqrt
  real (kind=R8) :: b2sign
  external b2sign

  flow = 0.0_R8
  do is = TO_SELF, TO_TOP
    weight(:,:,is) = vol(:,:)
  end do

  call value_on_faces1(nx,ny,ns,weight,density,den)
  do is = 0, ns-1
    call interp_volume(TO_LEFT,nx,ny,vol,gs,qc,qcb,velocity(-1,-1,0,is),vv(-1,-1,0,is))
    call interp_volume(TO_BOTTOM,nx,ny,vol,gs,qc,qcb,velocity(-1,-1,1,is),vv(-1,-1,1,is))
    do iy = -1, ny
      do ix = -1, nx
        if (isUnusedCell(cflags(ix,iy,CELLFLAG_TYPE))) cycle
!! poloidal contributions
        if (isInDomain(nx,ny,leftix(ix,iy),leftiy(ix,iy))) &
             & flow(ix,iy,0,is) = flow(ix,iy,0,is) + &
             & vv(ix,iy,0,is)*den(ix,iy,0,is)*gs(ix,iy,0)*qc(ix,iy)
        if ((isBoundaryCell(cflags(ix,iy,CELLFLAG_TYPE)) .or. &
             &  isGhostCell(cflags(ix,iy,CELLFLAG_TYPE))) .and. &
             &  cflags(ix,iy,CELLFLAG_TOPFACE) /= GRID_UNDEFINED .and. &
             &  isInDomain(nx,ny,topix(ix,iy),topiy(ix,iy))) &
             & flow(topix(ix,iy),topiy(ix,iy),1,is) = &
             & flow(topix(ix,iy),topiy(ix,iy),1,is) + &
             & b2sign(1.0_R8,pbs(topix(ix,iy),topiy(ix,iy),1))* &
             & vv(ix,iy,0,is)*den(ix,iy,0,is)* &
             & gs(topix(ix,iy),topiy(ix,iy),1)* &
             & sqrt(1.0_R8 - qcb(topix(ix,iy),topiy(ix,iy))**2)
        if ((isBoundaryCell(cflags(ix,iy,CELLFLAG_TYPE)) .or. &
             &  isGhostCell(cflags(ix,iy,CELLFLAG_TYPE))) .and. &
             &  cflags(ix,iy,CELLFLAG_BOTTOMFACE) /= GRID_UNDEFINED .and. &
             &  isInDomain(nx,ny,bottomix(ix,iy),bottomiy(ix,iy)) ) &
             & flow(ix,iy,1,is) = flow(ix,iy,1,is) + &
             & b2sign(1.0_R8,pbs(ix,iy,1))* &
             & vv(ix,iy,0,is)*den(ix,iy,0,is)* &
             & gs(ix,iy,1)*sqrt(1.0_R8 - qcb(ix,iy)**2)

!! radial contributions
        if (isInDomain(nx,ny,bottomix(ix,iy),bottomiy(ix,iy))) then
          if (pbs(ix,iy,1).eq.0.0_R8) then
            flow(ix,iy,1,is) = flow(ix,iy,1,is) + &
             & vv(ix,iy,1,is)*den(ix,iy,1,is)*gs(ix,iy,1)
          else
            flow(ix,iy,1,is) = flow(ix,iy,1,is) + b2sign(1.0_R8,pbs(ix,iy,1))* &
             & vv(ix,iy,1,is)*den(ix,iy,1,is)*gs(ix,iy,1)*qcb(ix,iy)
          endif
        endif
        cgeo = cellGeoType(crx(ix,iy,:),cry(ix,iy,:))
        if (cgeo == CGEO_TRIA_NOLEFT) then
          if (cflags(ix,iy,CELLFLAG_RIGHTFACE) /= GRID_UNDEFINED .and. &
            & cflags(ix,iy,CELLFLAG_TOPFACE) /= GRID_UNDEFINED ) &
             & flow(rightix(ix,iy),rightiy(ix,iy),0,is) = &
             & flow(rightix(ix,iy),rightiy(ix,iy),0,is) + &
             & vv(ix,iy,1,is)*den(ix,iy,1,is)* &
             & gs(rightix(ix,iy),rightiy(ix,iy),0)* &
             & sqrt(1.0_R8 - qc(rightix(ix,iy),rightiy(ix,iy))**2)
          if (cflags(ix,iy,CELLFLAG_RIGHTFACE) /= GRID_UNDEFINED .and. &
            & cflags(ix,iy,CELLFLAG_BOTTOMFACE) /= GRID_UNDEFINED ) &
             & flow(rightix(ix,iy),rightiy(ix,iy),0,is) = &
             & flow(rightix(ix,iy),rightiy(ix,iy),0,is) - &
             & vv(ix,iy,1,is)*den(ix,iy,1,is)* &
             & gs(rightix(ix,iy),rightiy(ix,iy),0)* &
             & sqrt(1.0_R8 - qc(rightix(ix,iy),rightiy(ix,iy))**2)
        else if (cgeo == CGEO_TRIA_NORIGHT) then
          if (cflags(ix,iy,CELLFLAG_LEFTFACE) /= GRID_UNDEFINED .and. &
            & cflags(ix,iy,CELLFLAG_TOPFACE) /= GRID_UNDEFINED) &
             & flow(ix,iy,0,is) = flow(ix,iy,0,is) - &
             & vv(ix,iy,1,is)*den(ix,iy,1,is)* &
             & gs(ix,iy,0)* sqrt(1.0_R8 - qc(ix,iy)**2)
          if (cflags(ix,iy,CELLFLAG_LEFTFACE) /= GRID_UNDEFINED .and. &
            & cflags(ix,iy,CELLFLAG_BOTTOMFACE) /= GRID_UNDEFINED ) &
             & flow(ix,iy,0,is) = flow(ix,iy,0,is) + &
             & vv(ix,iy,1,is)*den(ix,iy,1,is)* &
             & gs(ix,iy,0)* sqrt(1.0_R8 - qc(ix,iy)**2)
        else if (cgeo == CGEO_TRIA_NOBOT) then
          if (cflags(ix,iy,CELLFLAG_LEFTFACE) /= GRID_UNDEFINED) &
             & flow(ix,iy,0,is) = flow(ix,iy,0,is) + &
             & vv(ix,iy,1,is)*den(ix,iy,1,is)* &
             & gs(ix,iy,0)* sqrt(1.0_R8 - qc(ix,iy)**2)
          if (cflags(ix,iy,CELLFLAG_RIGHTFACE) /= GRID_UNDEFINED) &
             & flow(rightix(ix,iy),rightiy(ix,iy),0,is) = &
             & flow(rightix(ix,iy),rightiy(ix,iy),0,is) - &
             & vv(ix,iy,1,is)*den(ix,iy,1,is)* &
             & gs(rightix(ix,iy),rightiy(ix,iy),0)* &
             & sqrt(1.0_R8 - qc(rightix(ix,iy),rightiy(ix,iy))**2)
        else if (cgeo == CGEO_TRIA_NOTOP) then
          if (cflags(ix,iy,CELLFLAG_LEFTFACE) /= GRID_UNDEFINED) &
             & flow(ix,iy,0,is) = flow(ix,iy,0,is) - &
             & vv(ix,iy,1,is)*den(ix,iy,1,is)* &
             & gs(ix,iy,0)* sqrt(1.0_R8 - qc(ix,iy)**2)
          if (cflags(ix,iy,CELLFLAG_RIGHTFACE) /= GRID_UNDEFINED) &
             & flow(rightix(ix,iy),rightiy(ix,iy),0,is) = &
             & flow(rightix(ix,iy),rightiy(ix,iy),0,is) + &
             & vv(ix,iy,1,is)*den(ix,iy,1,is)* &
             & gs(rightix(ix,iy),rightiy(ix,iy),0)* &
             & sqrt(1.0_R8 - qc(rightix(ix,iy),rightiy(ix,iy))**2)
        end if
      end do
    end do
  end do

  return
  end subroutine flow_from_velocity

  !> Interpolate a cell quantity cv to all cell vertices
  function interpolateToAllVertices( nx, ny, cv ) result ( vx )
    use b2mod_cellhelper
    implicit none
    integer, intent(in) :: nx, ny
    real (kind=R8), intent(in), dimension(-1:nx,-1:ny) :: cv
    real (kind=R8), dimension(-1:nx,-1:ny,0:3) :: vx

    !! internal
    integer ivertex

    do ivertex = VX_LOWERLEFT, VX_UPPERRIGHT
      vx(:,:,ivertex) = interpolateToVertices( nx, ny, ivertex, cv)
    end do
    return
  end function interpolateToAllVertices

  !> Interpolate a cell quantity cv to a particular vertex
  function interpolateToVertices( nx, ny, vx_index, cv ) result( vx )
    use b2mod_geo
    use b2mod_b2cmfs
    use b2mod_indirect
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

  !> The purpose of this routine is to provide values of cell-centered quantities
  !> projected to the 'sides' of the cells, i.e. taking into account the fact that,
  !> on each side, more that one neighbor may contribute, or, in the case of
  !> triangular cells, the contribution does necessarily come from the obvious
  !> place.
  !> The projected values will be stored in side(,,TO_LEFT:TO_TOP).
  !> The interpolation to the sides is weighted by the 'weight' function.
  subroutine value_to_side(nx, ny, weight, centre, side)
  use b2mod_geo , only: crx, cry, gs, qc, qcb, pbs, vol
  use b2mod_indirect
  use b2mod_cellhelper
  implicit none
  integer, intent(in) :: nx, ny
  real(R8), intent(in) :: centre(-1:nx,-1:ny)
  real(R8), intent(in) :: weight(-1:nx,-1:ny,0:4)
  real(R8), intent(out) :: side(-1:nx,-1:ny,1:4)
  real (kind=R8) :: vol1(-1:nx,-1:ny,0:4)

  integer i, ix, iy, cgeo
  real (kind=R8) :: area_to_top, area_to_bottom, area_to_left, area_to_right
  logical classical, rectangular
  intrinsic sqrt

  side = INVALID_DOUBLE

  classical = isClassicalGrid(cflags)

  do i = TO_SELF, TO_TOP
    vol1(:,:,i) = vol(:,:)
  end do
  do ix = -1, nx
    do iy = -1, ny
      if (isUnusedCell(cflags(ix,iy,CELLFLAG_TYPE))) cycle
      rectangular = classical .or. cflags(ix,iy,CELLFLAG_TYPE) == GRID_INTERNAL
      if (rectangular .or. cflags(ix,iy,CELLFLAG_TYPE) == GRID_GUARD) then
        if (isInDomain(nx,ny,leftix(ix,iy),leftiy(ix,iy))) then
          side(ix,iy,TO_LEFT) = &
           & (centre(ix,iy)*weight(leftix(ix,iy),leftiy(ix,iy),TO_RIGHT)+ &
           &  centre(leftix(ix,iy),leftiy(ix,iy))*weight(ix,iy,TO_LEFT))/ &
           & (weight(leftix(ix,iy),leftiy(ix,iy),TO_RIGHT)+ &
           &  weight(ix,iy,TO_LEFT))
        else
          side(ix,iy,TO_LEFT) = centre(ix,iy)
        end if
        if (isInDomain(nx,ny,bottomix(ix,iy),bottomiy(ix,iy))) then
          side(ix,iy,TO_BOTTOM) = &
           & (centre(ix,iy)*weight(bottomix(ix,iy),bottomiy(ix,iy),TO_TOP)+ &
           &  centre(bottomix(ix,iy),bottomiy(ix,iy))*weight(ix,iy,TO_BOTTOM))/ &
           & (weight(bottomix(ix,iy),bottomiy(ix,iy),TO_TOP)+ &
           &  weight(ix,iy,TO_BOTTOM))
        else
          side(ix,iy,TO_BOTTOM) = centre(ix,iy)
        end if
        if (isInDomain(nx,ny,rightix(ix,iy),rightiy(ix,iy))) then
          side(ix,iy,TO_RIGHT) = &
           & (centre(ix,iy)*weight(rightix(ix,iy),rightiy(ix,iy),TO_LEFT)+ &
           &  centre(rightix(ix,iy),rightiy(ix,iy))*weight(ix,iy,TO_RIGHT))/ &
           & (weight(rightix(ix,iy),rightiy(ix,iy),TO_LEFT)+ &
           &  weight(ix,iy,TO_RIGHT))
        else
          side(ix,iy,TO_RIGHT) = centre(ix,iy)
        end if
        if (isInDomain(nx,ny,topix(ix,iy),topiy(ix,iy))) then
          side(ix,iy,TO_TOP) = &
           & (centre(ix,iy)*weight(topix(ix,iy),topiy(ix,iy),TO_BOTTOM)+ &
           &  centre(topix(ix,iy),topiy(ix,iy))*weight(ix,iy,TO_TOP))/ &
           & (weight(topix(ix,iy),topiy(ix,iy),TO_BOTTOM)+ &
           &  weight(ix,iy,TO_TOP))
        else
          side(ix,iy,TO_TOP) = centre(ix,iy)
        end if
      else
        cgeo = cellGeoType( crx(ix,iy,:), cry(ix,iy,:) )
!! Do left side first
        if (isInDomain(nx,ny,leftix(ix,iy),leftiy(ix,iy))) then
          area_to_left = gs(ix,iy,0) * qc(ix,iy)
        else
          area_to_left = 0.0_R8
        endif
        if (isInDomain(nx,ny,bottomix(ix,iy),bottomiy(ix,iy))) then
          area_to_bottom = gs(ix,iy,1)*sqrt(1.0_R8 - qcb(ix,iy)**2)
        else
          area_to_bottom = 0.0_R8
        end if
        if (cgeo == CGEO_TRIA_NOLEFT) then
!! Triangles with no left face
          area_to_top = gs(topix(ix,iy),topiy(ix,iy),1)* &
                 &  sqrt(1.0_R8 - qcb(topix(ix,iy),topiy(ix,iy))**2)
          if (pbs(topix(ix,iy),topiy(ix,iy),1).lt.0.0_R8 .and. &
            & pbs(ix,iy,1).gt.0.0_R8) then  !! cell is a left-pointing arrow
            side(ix,iy,TO_LEFT) = &
             & (area_to_top * &
             &  (centre(ix,iy)*weight(topix(ix,iy),topiy(ix,iy),TO_BOTTOM) + &
             &   centre(topix(ix,iy),topiy(ix,iy))*weight(ix,iy,TO_TOP))/ &
             &  (weight(topix(ix,iy),topiy(ix,iy),TO_BOTTOM) + &
             &   weight(ix,iy,TO_TOP)) + &
             &  area_to_bottom * &
             &  (centre(ix,iy)*weight(bottomix(ix,iy),bottomiy(ix,iy),TO_TOP) + &
             &   centre(bottomix(ix,iy),bottomiy(ix,iy))*weight(ix,iy,TO_BOTTOM))/ &
             &  (weight(bottomix(ix,iy),bottomiy(ix,iy),TO_TOP) + &
             &   weight(ix,iy,TO_BOTTOM))) / &
             & (area_to_top + area_to_bottom)
          else if (pbs(topix(ix,iy),topiy(ix,iy),1).lt.0.0_R8 .and. &
                &  pbs(ix,iy,1).le.0.0_R8) then
            side(ix,iy,TO_LEFT) = &
             &  (centre(ix,iy)*weight(topix(ix,iy),topiy(ix,iy),TO_BOTTOM) + &
             &   centre(topix(ix,iy),topiy(ix,iy))*weight(ix,iy,TO_TOP))/ &
             &  (weight(topix(ix,iy),topiy(ix,iy),TO_BOTTOM) + &
             &   weight(ix,iy,TO_TOP))
          else if (pbs(topix(ix,iy),topiy(ix,iy),1).ge.0.0_R8 .and. &
                &  pbs(ix,iy,1).gt.0.0_R8) then
            side(ix,iy,TO_LEFT) = &
             &  (centre(ix,iy)*weight(bottomix(ix,iy),bottomiy(ix,iy),TO_TOP) + &
             &   centre(bottomix(ix,iy),bottomiy(ix,iy))*weight(ix,iy,TO_BOTTOM))/ &
             &  (weight(bottomix(ix,iy),bottomiy(ix,iy),TO_TOP) + &
             &   weight(ix,iy,TO_BOTTOM))
          end if
        else if (vol1(ix,iy,TO_LEFT).eq.vol1(ix,iy,TO_SELF)) then
!! Other triangles and trapezoids with a full left side
          side(ix,iy,TO_LEFT) = &
           & (centre(ix,iy)*weight(leftix(ix,iy),leftiy(ix,iy),TO_RIGHT)+ &
           &  centre(leftix(ix,iy),leftiy(ix,iy))*weight(ix,iy,TO_LEFT))/ &
           & (weight(leftix(ix,iy),leftiy(ix,iy),TO_RIGHT)+ &
           &  weight(ix,iy,TO_LEFT))
        else
!! Other trapezoids
          area_to_top = gs(topix(ix,iy),topiy(ix,iy),1)* &
                 &  sqrt(1.0_R8 - qcb(topix(ix,iy),topiy(ix,iy))**2)
          if (pbs(topix(ix,iy),topiy(ix,iy),1).lt.0.0_R8 .and. &
            & pbs(ix,iy,1).gt.0.0_R8) then  !! contribution from all three neighbors
            side(ix,iy,TO_LEFT) = &
             & (area_to_top * &
             &  (centre(ix,iy)*weight(topix(ix,iy),topiy(ix,iy),TO_BOTTOM) + &
             &   centre(topix(ix,iy),topiy(ix,iy))*weight(ix,iy,TO_TOP))/ &
             &  (weight(topix(ix,iy),topiy(ix,iy),TO_BOTTOM) + &
             &   weight(ix,iy,TO_TOP)) + &
             &  area_to_left * &
             &  (centre(ix,iy)*weight(leftix(ix,iy),leftiy(ix,iy),TO_RIGHT) + &
             &   centre(leftix(ix,iy),leftiy(ix,iy))*weight(ix,iy,TO_LEFT))/ &
             &  (weight(leftix(ix,iy),leftiy(ix,iy),TO_RIGHT) + &
             &   weight(ix,iy,TO_LEFT)) + &
             &  area_to_bottom * &
             &  (centre(ix,iy)*weight(bottomix(ix,iy),bottomiy(ix,iy),TO_TOP) + &
             &   centre(bottomix(ix,iy),bottomiy(ix,iy))*weight(ix,iy,TO_BOTTOM))/ &
             &  (weight(bottomix(ix,iy),bottomiy(ix,iy),TO_TOP) + &
             &   weight(ix,iy,TO_BOTTOM))) / &
             & (area_to_top + area_to_left + area_to_bottom)
          else if (pbs(topix(ix,iy),topiy(ix,iy),1).lt.0.0_R8 .and. &
                &  pbs(ix,iy,1).le.0.0_R8) then  !! contribution from top and left neighbors
            side(ix,iy,TO_LEFT) = &
             & (area_to_top * &
             &  (centre(ix,iy)*weight(topix(ix,iy),topiy(ix,iy),TO_BOTTOM) + &
             &   centre(topix(ix,iy),topiy(ix,iy))*weight(ix,iy,TO_TOP))/ &
             &  (weight(topix(ix,iy),topiy(ix,iy),TO_BOTTOM) + &
             &   weight(ix,iy,TO_TOP)) + &
             &  area_to_left * &
             &  (centre(ix,iy)*weight(leftix(ix,iy),leftiy(ix,iy),TO_RIGHT) + &
             &   centre(leftix(ix,iy),leftiy(ix,iy))*weight(ix,iy,TO_LEFT))/ &
             &  (weight(leftix(ix,iy),leftiy(ix,iy),TO_RIGHT) + &
             &   weight(ix,iy,TO_LEFT))) / &
             & (area_to_top + area_to_left)
          else if (pbs(topix(ix,iy),topiy(ix,iy),1).ge.0.0_R8 .and. &
                &  pbs(ix,iy,1).gt.0.0_R8) then  !! contribution from bottom and left neighbors
            side(ix,iy,TO_LEFT) = &
             & (area_to_left * &
             &  (centre(ix,iy)*weight(leftix(ix,iy),leftiy(ix,iy),TO_RIGHT) + &
             &   centre(leftix(ix,iy),leftiy(ix,iy))*weight(ix,iy,TO_LEFT))/ &
             &  (weight(leftix(ix,iy),leftiy(ix,iy),TO_RIGHT) + &
             &   weight(ix,iy,TO_LEFT)) + &
             &  area_to_bottom * &
             &  (centre(ix,iy)*weight(bottomix(ix,iy),bottomiy(ix,iy),TO_TOP) + &
             &   centre(bottomix(ix,iy),bottomiy(ix,iy))*weight(ix,iy,TO_BOTTOM))/ &
             &  (weight(bottomix(ix,iy),bottomiy(ix,iy),TO_TOP) + &
             &   weight(ix,iy,TO_BOTTOM))) / &
             & (area_to_left + area_to_bottom)
          end if
        end if
!! Do bottom side second
        area_to_bottom = gs(ix,iy,1) * qcb(ix,iy)
        area_to_left = gs(ix,iy,0)*sqrt(1.0_R8 - qc(ix,iy)**2)
        if (cgeo == CGEO_TRIA_NOBOT) then
!! Triangles with no bottom face
          area_to_right = gs(rightix(ix,iy),rightiy(ix,iy),0)* &
                 &  sqrt(1.0_R8 - qc(rightix(ix,iy),rightiy(ix,iy))**2)
          if (cflags(ix,iy,CELLFLAG_LEFTFACE) /= GRID_UNDEFINED .and. &
           &  cflags(ix,iy,CELLFLAG_RIGHTFACE) /= GRID_UNDEFINED) then
            side(ix,iy,TO_BOTTOM) = &
             & (area_to_left * &
             &  (centre(ix,iy)*weight(leftix(ix,iy),leftiy(ix,iy),TO_RIGHT) + &
             &   centre(leftix(ix,iy),leftiy(ix,iy))*weight(ix,iy,TO_LEFT))/ &
             &  (weight(leftix(ix,iy),leftiy(ix,iy),TO_RIGHT) + &
             &   weight(ix,iy,TO_LEFT)) + &
             &  area_to_right * &
             &  (centre(ix,iy)*weight(rightix(ix,iy),rightiy(ix,iy),TO_LEFT) + &
             &   centre(rightix(ix,iy),rightiy(ix,iy))*weight(ix,iy,TO_RIGHT))/ &
             &  (weight(rightix(ix,iy),rightiy(ix,iy),TO_LEFT) + &
             &   weight(ix,iy,TO_RIGHT))) / &
             & (area_to_left + area_to_right)
          else if (cflags(ix,iy,CELLFLAG_LEFTFACE) == GRID_UNDEFINED .and. &
                &  cflags(ix,iy,CELLFLAG_RIGHTFACE) /= GRID_UNDEFINED) then
            side(ix,iy,TO_BOTTOM) = &
             &  (centre(ix,iy)*weight(rightix(ix,iy),rightiy(ix,iy),TO_LEFT) + &
             &   centre(rightix(ix,iy),rightiy(ix,iy))*weight(ix,iy,TO_RIGHT))/ &
             &  (weight(rightix(ix,iy),rightiy(ix,iy),TO_LEFT) + &
             &   weight(ix,iy,TO_RIGHT))
          else if (cflags(ix,iy,CELLFLAG_LEFTFACE) /= GRID_UNDEFINED .and. &
                &  cflags(ix,iy,CELLFLAG_RIGHTFACE) == GRID_UNDEFINED) then
            side(ix,iy,TO_BOTTOM) = &
             &  (centre(ix,iy)*weight(leftix(ix,iy),leftiy(ix,iy),TO_RIGHT) + &
             &   centre(leftix(ix,iy),leftiy(ix,iy))*weight(ix,iy,TO_LEFT))/ &
             &  (weight(leftix(ix,iy),leftiy(ix,iy),TO_RIGHT) + &
             &   weight(ix,iy,TO_LEFT))
          end if
        else if (vol1(ix,iy,TO_BOTTOM).eq.vol1(ix,iy,TO_SELF)) then
!! Other triangles and trapezoids with a full bottom side
          side(ix,iy,TO_BOTTOM) = &
           & (centre(ix,iy)*weight(bottomix(ix,iy),bottomiy(ix,iy),TO_TOP)+ &
           &  centre(bottomix(ix,iy),bottomiy(ix,iy))*weight(ix,iy,TO_BOTTOM))/ &
           & (weight(bottomix(ix,iy),bottomiy(ix,iy),TO_TOP)+ &
           &  weight(ix,iy,TO_BOTTOM))
        else
!! Other trapezoids
          area_to_right = gs(rightix(ix,iy),rightiy(ix,iy),0)* &
                 &  sqrt(1.0_R8 - qc(rightix(ix,iy),rightiy(ix,iy))**2)
          if (cflags(ix,iy,CELLFLAG_LEFTFACE) /= GRID_UNDEFINED .and. &
           &  cflags(ix,iy,CELLFLAG_RIGHTFACE) /= GRID_UNDEFINED) then
            side(ix,iy,TO_BOTTOM) = &
             & (area_to_left * &
             &  (centre(ix,iy)*weight(leftix(ix,iy),leftiy(ix,iy),TO_RIGHT) + &
             &   centre(leftix(ix,iy),leftiy(ix,iy))*weight(ix,iy,TO_LEFT))/ &
             &  (weight(leftix(ix,iy),leftiy(ix,iy),TO_RIGHT) + &
             &   weight(ix,iy,TO_LEFT)) + &
             &  area_to_bottom * &
             &  (centre(ix,iy)*weight(bottomix(ix,iy),bottomiy(ix,iy),TO_TOP)+ &
             &   centre(bottomix(ix,iy),bottomiy(ix,iy))*weight(ix,iy,TO_BOTTOM))/ &
             &  (weight(bottomix(ix,iy),bottomiy(ix,iy),TO_TOP)+ &
             &   weight(ix,iy,TO_BOTTOM)) + &
             &  area_to_right * &
             &  (centre(ix,iy)*weight(rightix(ix,iy),rightiy(ix,iy),TO_LEFT) + &
             &   centre(rightix(ix,iy),rightiy(ix,iy))*weight(ix,iy,TO_RIGHT))/ &
             &  (weight(rightix(ix,iy),rightiy(ix,iy),TO_LEFT) + &
             &   weight(ix,iy,TO_RIGHT))) / &
             & (area_to_left + area_to_bottom + area_to_right)
          elseif (cflags(ix,iy,CELLFLAG_LEFTFACE) == GRID_UNDEFINED .and. &
           &  cflags(ix,iy,CELLFLAG_RIGHTFACE) /= GRID_UNDEFINED) then
            side(ix,iy,TO_BOTTOM) = &
             & (area_to_bottom * &
             &  (centre(ix,iy)*weight(bottomix(ix,iy),bottomiy(ix,iy),TO_TOP)+ &
             &   centre(bottomix(ix,iy),bottomiy(ix,iy))*weight(ix,iy,TO_BOTTOM))/ &
             &  (weight(bottomix(ix,iy),bottomiy(ix,iy),TO_TOP)+ &
             &   weight(ix,iy,TO_BOTTOM)) + &
             &  area_to_right * &
             &  (centre(ix,iy)*weight(rightix(ix,iy),rightiy(ix,iy),TO_LEFT) + &
             &   centre(rightix(ix,iy),rightiy(ix,iy))*weight(ix,iy,TO_RIGHT))/ &
             &  (weight(rightix(ix,iy),rightiy(ix,iy),TO_LEFT) + &
             &   weight(ix,iy,TO_RIGHT))) / &
             & (area_to_bottom + area_to_right)
          elseif (cflags(ix,iy,CELLFLAG_LEFTFACE) /= GRID_UNDEFINED .and. &
           &  cflags(ix,iy,CELLFLAG_RIGHTFACE) == GRID_UNDEFINED) then
            side(ix,iy,TO_BOTTOM) = &
             & (area_to_left * &
             &  (centre(ix,iy)*weight(leftix(ix,iy),leftiy(ix,iy),TO_RIGHT) + &
             &   centre(leftix(ix,iy),leftiy(ix,iy))*weight(ix,iy,TO_LEFT))/ &
             &  (weight(leftix(ix,iy),leftiy(ix,iy),TO_RIGHT) + &
             &   weight(ix,iy,TO_LEFT)) + &
             &  area_to_bottom * &
             &  (centre(ix,iy)*weight(bottomix(ix,iy),bottomiy(ix,iy),TO_TOP)+ &
             &   centre(bottomix(ix,iy),bottomiy(ix,iy))*weight(ix,iy,TO_BOTTOM))/ &
             &  (weight(bottomix(ix,iy),bottomiy(ix,iy),TO_TOP)+ &
             &   weight(ix,iy,TO_BOTTOM))) / &
             & (area_to_bottom + area_to_left)
          end if
        end if
!! Do right side third
        if (cgeo == CGEO_TRIA_NORIGHT) then
!! Triangles with no right face
          area_to_bottom = gs(ix,iy,1)*sqrt(1.0_R8 - qcb(ix,iy)**2)
          area_to_top = gs(topix(ix,iy),topiy(ix,iy),1)* &
                 &  sqrt(1.0_R8 - qcb(topix(ix,iy),topiy(ix,iy))**2)
          if (pbs(topix(ix,iy),topiy(ix,iy),1).gt.0.0_R8 .and. &
            & pbs(ix,iy,1).lt.0.0_R8) then  !! cell is a right-pointing arrow
            side(ix,iy,TO_RIGHT) = &
             & (area_to_top * &
             &  (centre(ix,iy)*weight(topix(ix,iy),topiy(ix,iy),TO_BOTTOM) + &
             &   centre(topix(ix,iy),topiy(ix,iy))*weight(ix,iy,TO_TOP))/ &
             &  (weight(topix(ix,iy),topiy(ix,iy),TO_BOTTOM) + &
             &   weight(ix,iy,TO_TOP)) + &
             &  area_to_bottom * &
             &  (centre(ix,iy)*weight(bottomix(ix,iy),bottomiy(ix,iy),TO_TOP) + &
             &   centre(bottomix(ix,iy),bottomiy(ix,iy))*weight(ix,iy,TO_BOTTOM))/ &
             &  (weight(bottomix(ix,iy),bottomiy(ix,iy),TO_TOP) + &
             &   weight(ix,iy,TO_BOTTOM))) / &
             & (area_to_top + area_to_bottom)
          else if (pbs(topix(ix,iy),topiy(ix,iy),1).gt.0.0_R8 .and. &
                &  pbs(ix,iy,1).ge.0.0_R8) then
            side(ix,iy,TO_RIGHT) = &
             &  (centre(ix,iy)*weight(topix(ix,iy),topiy(ix,iy),TO_BOTTOM) + &
             &   centre(topix(ix,iy),topiy(ix,iy))*weight(ix,iy,TO_TOP))/ &
             &  (weight(topix(ix,iy),topiy(ix,iy),TO_BOTTOM) + &
             &   weight(ix,iy,TO_TOP))
          else if (pbs(topix(ix,iy),topiy(ix,iy),1).le.0.0_R8 .and. &
                &  pbs(ix,iy,1).lt.0.0_R8) then
            side(ix,iy,TO_RIGHT) = &
             &  (centre(ix,iy)*weight(bottomix(ix,iy),bottomiy(ix,iy),TO_TOP) + &
             &   centre(bottomix(ix,iy),bottomiy(ix,iy))*weight(ix,iy,TO_BOTTOM))/ &
             &  (weight(bottomix(ix,iy),bottomiy(ix,iy),TO_TOP) + &
             &   weight(ix,iy,TO_BOTTOM))
          end if
        else if (vol1(ix,iy,TO_RIGHT).eq.vol1(ix,iy,TO_SELF)) then
!! Other triangles and trapezoids with a full right side
          side(ix,iy,TO_RIGHT) = &
           & (centre(ix,iy)*weight(rightix(ix,iy),rightiy(ix,iy),TO_LEFT)+ &
           &  centre(rightix(ix,iy),rightiy(ix,iy))*weight(ix,iy,TO_RIGHT))/ &
           & (weight(rightix(ix,iy),rightiy(ix,iy),TO_LEFT)+ &
           &  weight(ix,iy,TO_RIGHT))
        else
!! Other trapezoids
          area_to_right = gs(rightix(ix,iy),rightiy(ix,iy),0) * &
                        & qc(rightix(ix,iy),rightiy(ix,iy))
          area_to_bottom = gs(ix,iy,1)*sqrt(1.0_R8 - qcb(ix,iy)**2)
          area_to_top = gs(topix(ix,iy),topiy(ix,iy),1)* &
                 &  sqrt(1.0_R8 - qcb(topix(ix,iy),topiy(ix,iy))**2)
          if (pbs(topix(ix,iy),topiy(ix,iy),1).gt.0.0_R8 .and. &
            & pbs(ix,iy,1).lt.0.0_R8) then  !! contribution from all three neighbors
            side(ix,iy,TO_RIGHT) = &
             & (area_to_top * &
             &  (centre(ix,iy)*weight(topix(ix,iy),topiy(ix,iy),TO_BOTTOM) + &
             &   centre(topix(ix,iy),topiy(ix,iy))*weight(ix,iy,TO_TOP))/ &
             &  (weight(topix(ix,iy),topiy(ix,iy),TO_BOTTOM) + &
             &   weight(ix,iy,TO_TOP)) + &
             &  area_to_right * &
             &  (centre(ix,iy)*weight(rightix(ix,iy),rightiy(ix,iy),TO_LEFT) + &
             &   centre(rightix(ix,iy),rightiy(ix,iy))*weight(ix,iy,TO_RIGHT))/ &
             &  (weight(rightix(ix,iy),rightiy(ix,iy),TO_LEFT) + &
             &   weight(ix,iy,TO_RIGHT)) + &
             &  area_to_bottom * &
             &  (centre(ix,iy)*weight(bottomix(ix,iy),bottomiy(ix,iy),TO_TOP) + &
             &   centre(bottomix(ix,iy),bottomiy(ix,iy))*weight(ix,iy,TO_BOTTOM))/ &
             &  (weight(bottomix(ix,iy),bottomiy(ix,iy),TO_TOP) + &
             &   weight(ix,iy,TO_BOTTOM))) / &
             & (area_to_top + area_to_right + area_to_bottom)
          else if (pbs(topix(ix,iy),topiy(ix,iy),1).gt.0.0_R8 .and. &
                &  pbs(ix,iy,1).ge.0.0_R8) then  !! contribution from top and right neighbors
            side(ix,iy,TO_RIGHT) = &
             & (area_to_top * &
             &  (centre(ix,iy)*weight(topix(ix,iy),topiy(ix,iy),TO_BOTTOM) + &
             &   centre(topix(ix,iy),topiy(ix,iy))*weight(ix,iy,TO_TOP))/ &
             &  (weight(topix(ix,iy),topiy(ix,iy),TO_BOTTOM) + &
             &   weight(ix,iy,TO_TOP)) + &
             &  area_to_right * &
             &  (centre(ix,iy)*weight(rightix(ix,iy),rightiy(ix,iy),TO_LEFT) + &
             &   centre(rightix(ix,iy),rightiy(ix,iy))*weight(ix,iy,TO_RIGHT))/ &
             &  (weight(rightix(ix,iy),rightiy(ix,iy),TO_LEFT) + &
             &   weight(ix,iy,TO_RIGHT))) / &
             & (area_to_top + area_to_right)
          else if (pbs(topix(ix,iy),topiy(ix,iy),1).le.0.0_R8 .and. &
                &  pbs(ix,iy,1).lt.0.0_R8) then  !! contribution from bottom and right neighbors
            side(ix,iy,TO_RIGHT) = &
             & (area_to_right * &
             &  (centre(ix,iy)*weight(rightix(ix,iy),rightiy(ix,iy),TO_LEFT) + &
             &   centre(rightix(ix,iy),rightiy(ix,iy))*weight(ix,iy,TO_RIGHT))/ &
             &  (weight(rightix(ix,iy),rightiy(ix,iy),TO_LEFT) + &
             &   weight(ix,iy,TO_RIGHT)) + &
             &  area_to_bottom * &
             &  (centre(ix,iy)*weight(bottomix(ix,iy),bottomiy(ix,iy),TO_TOP) + &
             &   centre(bottomix(ix,iy),bottomiy(ix,iy))*weight(ix,iy,TO_BOTTOM))/ &
             &  (weight(bottomix(ix,iy),bottomiy(ix,iy),TO_TOP) + &
             &   weight(ix,iy,TO_BOTTOM))) / &
             & (area_to_right + area_to_bottom)
          end if
        end if
!! Do top side fourth
        if (cgeo == CGEO_TRIA_NOTOP) then
!! Triangles with no top face
          area_to_left = gs(ix,iy,0)*sqrt(1.0_R8 - qc(ix,iy)**2)
          area_to_right = gs(rightix(ix,iy),rightiy(ix,iy),0)* &
                 &  sqrt(1.0_R8 - qc(rightix(ix,iy),rightiy(ix,iy))**2)
          if (cflags(ix,iy,CELLFLAG_LEFTFACE) /= GRID_UNDEFINED .and. &
           &  cflags(ix,iy,CELLFLAG_RIGHTFACE) /= GRID_UNDEFINED) then
            side(ix,iy,TO_TOP) = &
             & (area_to_left * &
             &  (centre(ix,iy)*weight(leftix(ix,iy),leftiy(ix,iy),TO_RIGHT) + &
             &   centre(leftix(ix,iy),leftiy(ix,iy))*weight(ix,iy,TO_LEFT))/ &
             &  (weight(leftix(ix,iy),leftiy(ix,iy),TO_RIGHT) + &
             &   weight(ix,iy,TO_LEFT)) + &
             &  area_to_right * &
             &  (centre(ix,iy)*weight(rightix(ix,iy),rightiy(ix,iy),TO_LEFT) + &
             &   centre(rightix(ix,iy),rightiy(ix,iy))*weight(ix,iy,TO_RIGHT))/ &
             &  (weight(rightix(ix,iy),rightiy(ix,iy),TO_LEFT) + &
             &   weight(ix,iy,TO_RIGHT))) / &
             & (area_to_left + area_to_right)
          else if (cflags(ix,iy,CELLFLAG_LEFTFACE) == GRID_UNDEFINED .and. &
                &  cflags(ix,iy,CELLFLAG_RIGHTFACE) /= GRID_UNDEFINED) then
            side(ix,iy,TO_TOP) = &
             &  (centre(ix,iy)*weight(rightix(ix,iy),rightiy(ix,iy),TO_LEFT) + &
             &   centre(rightix(ix,iy),rightiy(ix,iy))*weight(ix,iy,TO_RIGHT))/ &
             &  (weight(rightix(ix,iy),rightiy(ix,iy),TO_LEFT) + &
             &   weight(ix,iy,TO_RIGHT))
          else if (cflags(ix,iy,CELLFLAG_LEFTFACE) /= GRID_UNDEFINED .and. &
                &  cflags(ix,iy,CELLFLAG_RIGHTFACE) == GRID_UNDEFINED) then
            side(ix,iy,TO_TOP) = &
             &  (centre(ix,iy)*weight(leftix(ix,iy),leftiy(ix,iy),TO_RIGHT) + &
             &   centre(leftix(ix,iy),leftiy(ix,iy))*weight(ix,iy,TO_LEFT))/ &
             &  (weight(leftix(ix,iy),leftiy(ix,iy),TO_RIGHT) + &
             &   weight(ix,iy,TO_LEFT))
          end if
        else if (vol1(ix,iy,TO_TOP).eq.vol1(ix,iy,TO_SELF)) then
!! Other triangles and trapezoids with a full top side
          side(ix,iy,TO_TOP) = &
           & (centre(ix,iy)*weight(topix(ix,iy),topiy(ix,iy),TO_BOTTOM)+ &
           &  centre(topix(ix,iy),topiy(ix,iy))*weight(ix,iy,TO_TOP))/ &
           & (weight(topix(ix,iy),topiy(ix,iy),TO_BOTTOM)+ &
           &  weight(ix,iy,TO_TOP))
        else
!! Other trapezoids
          area_to_top = gs(topix(ix,iy),topiy(ix,iy),1) * &
                     & qcb(topix(ix,iy),topiy(ix,iy))
          area_to_left = gs(ix,iy,0)*sqrt(1.0_R8 - qc(ix,iy)**2)
          area_to_right = gs(rightix(ix,iy),rightiy(ix,iy),0)* &
                 &  sqrt(1.0_R8 - qc(rightix(ix,iy),rightiy(ix,iy))**2)
          if (cflags(ix,iy,CELLFLAG_LEFTFACE) /= GRID_UNDEFINED .and. &
           &  cflags(ix,iy,CELLFLAG_RIGHTFACE) /= GRID_UNDEFINED) then
            side(ix,iy,TO_TOP) = &
             & (area_to_left * &
             &  (centre(ix,iy)*weight(leftix(ix,iy),leftiy(ix,iy),TO_RIGHT) + &
             &   centre(leftix(ix,iy),leftiy(ix,iy))*weight(ix,iy,TO_LEFT))/ &
             &  (weight(leftix(ix,iy),leftiy(ix,iy),TO_RIGHT) + &
             &   weight(ix,iy,TO_LEFT)) + &
             &  area_to_top * &
             &  (centre(ix,iy)*weight(topix(ix,iy),topiy(ix,iy),TO_BOTTOM)+ &
             &   centre(topix(ix,iy),topiy(ix,iy))*weight(ix,iy,TO_TOP))/ &
             &  (weight(topix(ix,iy),topiy(ix,iy),TO_BOTTOM)+ &
             &   weight(ix,iy,TO_TOP)) + &
             &  area_to_right * &
             &  (centre(ix,iy)*weight(rightix(ix,iy),rightiy(ix,iy),TO_LEFT) + &
             &   centre(rightix(ix,iy),rightiy(ix,iy))*weight(ix,iy,TO_RIGHT))/ &
             &  (weight(rightix(ix,iy),rightiy(ix,iy),TO_LEFT) + &
             &   weight(ix,iy,TO_RIGHT))) / &
             & (area_to_left + area_to_top + area_to_right)
          else if (cflags(ix,iy,CELLFLAG_LEFTFACE) == GRID_UNDEFINED .and. &
                &  cflags(ix,iy,CELLFLAG_RIGHTFACE) /= GRID_UNDEFINED) then
            side(ix,iy,TO_TOP) = &
             & (area_to_top * &
             &  (centre(ix,iy)*weight(topix(ix,iy),topiy(ix,iy),TO_BOTTOM)+ &
             &   centre(topix(ix,iy),topiy(ix,iy))*weight(ix,iy,TO_TOP))/ &
             &  (weight(topix(ix,iy),topiy(ix,iy),TO_BOTTOM)+ &
             &   weight(ix,iy,TO_TOP)) + &
             &  area_to_right * &
             &  (centre(ix,iy)*weight(rightix(ix,iy),rightiy(ix,iy),TO_LEFT) + &
             &   centre(rightix(ix,iy),rightiy(ix,iy))*weight(ix,iy,TO_RIGHT))/ &
             &  (weight(rightix(ix,iy),rightiy(ix,iy),TO_LEFT) + &
             &   weight(ix,iy,TO_RIGHT))) / &
             & (area_to_top + area_to_right)
          else if (cflags(ix,iy,CELLFLAG_LEFTFACE) /= GRID_UNDEFINED .and. &
                &  cflags(ix,iy,CELLFLAG_RIGHTFACE) == GRID_UNDEFINED) then
            side(ix,iy,TO_TOP) = &
             & (area_to_left * &
             &  (centre(ix,iy)*weight(leftix(ix,iy),leftiy(ix,iy),TO_RIGHT) + &
             &   centre(leftix(ix,iy),leftiy(ix,iy))*weight(ix,iy,TO_LEFT))/ &
             &  (weight(leftix(ix,iy),leftiy(ix,iy),TO_RIGHT) + &
             &   weight(ix,iy,TO_LEFT)) + &
             &  area_to_top * &
             &  (centre(ix,iy)*weight(topix(ix,iy),topiy(ix,iy),TO_BOTTOM)+ &
             &   centre(topix(ix,iy),topiy(ix,iy))*weight(ix,iy,TO_TOP))/ &
             &  (weight(topix(ix,iy),topiy(ix,iy),TO_BOTTOM)+ &
             &   weight(ix,iy,TO_TOP))) / &
             & (area_to_left + area_to_top)
          end if
        end if
      end if
    end do
  end do

  return
  end subroutine value_to_side

  function get_connection_length( iclass, imode, nx, ny, ix, iy, ivx )
  use b2mod_types
  use b2mod_geo
  use b2mod_indirect
#ifdef IMAS
#if IMAS_MINOR_VERSION > 8
  use ids_schemas  & ! IGNORE
     & , only : IDS_REAL_INVALID
#endif
#endif
  implicit none
  integer, intent(in) :: iclass, imode
  integer, intent(in) :: nx, ny, ix, iy
  integer, intent(in), optional :: ivx
  real(kind=R8) :: get_connection_length
  external xertst, xerrab

!  This function computes the connection length
!  For nodes (iclass.eq.1), we compute the length along a field line
!  For cells (iclass.eq.3), we compute the length of the flux tube
!  If imode.eq.1, the full length
!  If imode.eq.2, the length to closest solid surface

  integer ix1, ixl, ixr, iyl, iop
  real(kind=R8) :: t1, t2

#ifdef IMAS
#if IMAS_MINOR_VERSION > 8
  get_connection_length = IDS_REAL_INVALID
#else
  get_connection_length = 0.0_R8
#endif
#else
  get_connection_length = 0.0_R8
#endif

  call xertst(0.lt.nx, 'faulty input nx')
  call xertst(0.lt.ny, 'faulty input ny')
  call xertst(-1.le.ix .and. ix.le.nx, 'faulty input ix')
  call xertst(-1.le.iy .and. iy.le.ny, 'faulty input iy')
  call xertst(imode.eq.1 .or. imode.eq.2, 'faulty input imode')
  call xertst(iclass.eq.3 .or. present(ivx), &
     & 'iclass.eq.1 requires a vertex index !')
  if (region(ix,iy,0).eq.0) return

  select case ( iclass )
     case ( 1 )  ! IDS_CLASS_NODE
        select case ( ivx )
        case( 0 )
          iop = 1
          iyl = iy
        case( 1 )
          iop = 0
          iyl = iy
        case( 2 )
          iop = 3
          iyl = min(ny,topiy(ix,iy))
        case( 3 )
          iop = 2
          iyl = min(ny,topiy(ix,iy))
        case default
          call xerrab ('Unrecognized ivx value !')
        end select
        ixl=leftix(ix,iy)
        t1=0.0_R8
        if (ixl.ne.-2.and.(ivx.eq.1.or.ivx.eq.3)) then
          if (region(ix,iy,0).ne.0) t1=dist(ix,iy,ivx,iop)/ &
     &       (sqrt(wbbv(ix,iyl,0)**2+wbbv(ix,iyl,1)**2)/wbbv(ix,iyl,3))
        endif
        ix1=ixl
        if (ixl.ne.-2) ixl=leftix(ix1,iy)
        do while(ix1.ne.ix.and.ixl.ne.-2)
          if (region(ix1,iy,0).ne.0) t1=t1+dist(ix1,iy,ivx,iop)/ &
     &       (sqrt(wbbv(ix1,iyl,0)**2+wbbv(ix1,iyl,1)**2)/wbbv(ix1,iyl,3))
          ix1=ixl
          ixl=leftix(ix1,iy)
        end do
        if(ix1.eq.ix.and.(ivx.eq.0.or.ivx.eq.2)) t1=t1+dist(ix,iy,ivx,iop)/ &
     &    (sqrt(wbbv(ix1,iyl,0)**2+wbbv(ix1,iyl,1)**2)/wbbv(ix1,iyl,3))
        ixr=rightix(ix,iy)
        t2=0.0_R8
        if (ixr.ne.nx+1.and.(ivx.eq.0.or.ivx.eq.2)) then
          if (region(ix,iy,0).ne.0) t2=dist(ix,iy,ivx,iop)/ &
     &       (sqrt(wbbv(ix,iyl,0)**2+wbbv(ix,iyl,1)**2)/wbbv(ix,iyl,3))
        endif
        ix1=ixr
        if (ixr.ne.nx+1) ixr=rightix(ix1,iy)
        do while(ix1.ne.ix.and.ixr.ne.nx+1)
          if (region(ix1,iy,0).ne.0) t2=t2+dist(ix1,iy,ivx,iop)/ &
     &     (sqrt(wbbv(ix1,iyl,0)**2+wbbv(ix1,iyl,1)**2)/wbbv(ix1,iyl,3))
          ix1=ixr
          ixr=rightix(ix1,iy)
        end do
        if(ix1.eq.ix.and.(ivx.eq.1.or.ivx.eq.3)) t2=t2+dist(ix,iy,ivx,iop)/ &
     &    (sqrt(wbbv(ix1,iyl,0)**2+wbbv(ix1,iyl,1)**2)/wbbv(ix1,iyl,3))
        if(imode.eq.1.and.ix1.ne.ix) then
          get_connection_length = t1+t2
        else
          get_connection_length = min(t1,t2)
        end if
     case ( 3 )  ! IDS_CLASS_CELL
        ixl=leftix(ix,iy)
        t1=0.0_R8
        if (ixl.ne.-2) then
          if (region(ix,iy,0).ne.0) t1=hx(ix,iy)/2.0_R8/ &
     &       (sqrt(bb(ix,iy,0)**2+bb(ix,iy,1)**2)/bb(ix,iy,3))
        endif
        ix1=ixl
        if (ixl.ne.-2) ixl=leftix(ix1,iy)
        do while(ix1.ne.ix.and.ixl.ne.-2)
          if (region(ix1,iy,0).ne.0) t1=t1+hx(ix1,iy)/ &
     &       (sqrt(bb(ix1,iy,0)**2+bb(ix1,iy,1)**2)/bb(ix1,iy,3))
          ix1=ixl
          ixl=leftix(ix1,iy)
        end do
        if(ix1.eq.ix) t1=t1+hx(ix1,iy)/2.0_R8/ &
     &    (sqrt(bb(ix1,iy,0)**2+bb(ix1,iy,1)**2)/bb(ix1,iy,3))
        ixr=rightix(ix,iy)
        t2=0.0_R8
        if (ixr.ne.nx+1) then
          if (region(ix,iy,0).ne.0) t2=hx(ix,iy)/2.0_R8/ &
     &       (sqrt(bb(ix,iy,0)**2+bb(ix,iy,1)**2)/bb(ix,iy,3))
        endif
        ix1=ixr
        if (ixr.ne.nx+1) ixr=rightix(ix1,iy)
        do while(ix1.ne.ix.and.ixr.ne.nx+1)
          if (region(ix1,iy,0).ne.0) t2=t2+hx(ix1,iy)/ &
     &       (sqrt(bb(ix1,iy,0)**2+bb(ix1,iy,1)**2)/bb(ix1,iy,3))
          ix1=ixr
          ixr=rightix(ix1,iy)
        end do
        if(ix1.eq.ix) t2=t2+hx(ix1,iy)/2.0_R8/ &
     &    (sqrt(bb(ix1,iy,0)**2+bb(ix1,iy,1)**2)/bb(ix1,iy,3))
        if(imode.eq.1.and.ix1.ne.ix) then
          get_connection_length = t1+t2
        else
          get_connection_length = min(t1,t2)
        end if
     case default
        call xerrab ( 'Unsupported class in get_connection_length !')
  end select

  return
  end function get_connection_length

end module b2mod_interp

!!!Local Variables:
!!! mode: f90
!!! End:
