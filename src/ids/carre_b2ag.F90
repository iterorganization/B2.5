module carre_b2ag

  use b2mod_indirect
  use b2mod_interp
  implicit none

contains

  subroutine carre_b2agbb (nx,ny,cflag,fpsi,ffbz,wbbc,bb, & 
       &  crx,cry,psidx,psidy,isymm)

    !======================================================================
    !   ..input arguments (unchanged on exit)
    integer nx, ny, isymm
    integer cflag(-1:nx,-1:ny,CARREOUT_NCELLFLAGS)
    real*8 fpsi(-1:nx,-1:ny,0:3), ffbz(-1:nx,-1:ny,0:3)
    real*8 crx(-1:nx,-1:ny,0:3),cry(-1:nx,-1:ny,0:3), & 
         &  psidx(-1:nx,-1:ny,0:3),psidy(-1:nx,-1:ny,0:3)
    !   ..output arguments (unspecified on entry)
    real*8, intent(out) :: bb(-1:nx,-1:ny,0:3), wbbc(-1:nx+1,-1:ny+1,0:3)

    !-----------------------------------------------------------------------
    !.documentation
    !
    !  1. purpose
    !
    !     B2AGBB computes the magnetic field from the flux functions.
    !
    !
    !  3. description (see also routine b2cdca)
    !
    !     This routine computes the three components of the magnetic field
    !     and the total magnetic field strength, all at cell centers. The
    !     components (bx,by) are obtained by differencing the potential fpsi
    !     from its values at cell corners. The component bz is obtained
    !     directly from the field function ffbz.
    !
    !
    !  5. parameters (see also routine b2cdcv)
    !
    !     nx, ny - integer, input.
    !     See the calling routine for description.
    !
    !     vol, hx, hy, qz - real array, input.
    !     See the calling routine for description.
    !
    !     fpsi - (-1:nx,-1:ny,0:3) real array, input.
    !     For (ix,iy,i) in (-1:nx,-1:ny,0:3), fpsi(ix,iy,i) specifies the
    !     potential for the magnetic field components (bx,by) at the
    !     (ix,iy,i) position.
    !     (In the case of straight symmetry, bx=dpsi/dy and by=-dpsi/dx.
    !     In the case of toroidal symmetry, bx=(dpsi/dy)/(2*pi*R) and
    !     by=-(dpsi/dx)/(2*pi*R).)
    !
    !     ffbz - (-1:nx,-1:ny,0:3) real array, input.
    !     For (ix,iy,i) in (-1:nx,-1:ny,0:3), ffbz(ix,iy,i) specifies the
    !     field function for the magnetic field component bz at the
    !     (ix,iy,i) position.
    !     (In the case of straight symmetry, bz=ffbz. In the case of
    !     toroidal symmetry, bz=ffbz/(2*pi*R).)
    !
    !     bb - (-1:nx,-1:ny,0:3) real array, output.
    !     For (ix,iy) in (-1:nx,-1:ny), bb(ix,iy,0:3) specifies the
    !     magnetic field at the center of the (ix,iy) cell.
    !     bb(ix,iy,0:2) are the (x,y,z)-components and bb(ix,iy,3) is
    !     the absolute magnetic field strength.
    !     (bb(,,3)=sqrt(bb(,,0)**2+bb(,,1)**2+bb(,,2)**2.)
    !     It will hold that 0.lt.bb(,,3).
    !
    !-----------------------------------------------------------------------
    !.declarations

    !   ..local variables
    integer ix, iy
    real*8 t0,psdx,psdy,babs,pi
    real*8 centroid(0:1)
    !   ..procedures
    intrinsic sqrt,sign
    !-----------------------------------------------------------------------
    !.computation

    ! ..compute magnetic field
    pi=4.*atan(1.)
    !   ..loop over cells

    bb = BB_INVALID
    wbbc = BB_INVALID

    do iy = -1, ny
        do ix = -1, nx
            if (cflag(ix,iy,CELLFLAG_TYPE) == GRID_UNDEFINED) cycle
            !     ..compute magnetic field at cell center
            if (isymm.ne.0) then
!WG         WG-TMP back to simple averaging for consistency later on (b2agmt)  
              centroid(0) = 0.25_R8*(crx(ix,iy,0) + crx(ix,iy,1) + &
     &                               crx(ix,iy,2) + crx(ix,iy,3))
              centroid(1) = 0.25_R8*(cry(ix,iy,0) + cry(ix,iy,1) + &
     &                               cry(ix,iy,2) + cry(ix,iy,3))
!WG              centroid=quadCentroid(crx(ix,iy,0),cry(ix,iy,0),  &
!WG     &                              crx(ix,iy,1),cry(ix,iy,1),  &
!WG     &                              crx(ix,iy,2),cry(ix,iy,2),  &
!WG     &                              crx(ix,iy,3),cry(ix,iy,3))
              if (isymm.eq.1 .or. isymm.eq.2) then
                t0=2.0*pi*centroid(0)
              elseif (isymm.eq.3 .or. isymm.eq.4) then
                t0=2.0*pi*centroid(1)
              endif
            else
              t0=1.0
            endif

            if (t0 == 0.0) cycle

            bb(ix,iy,0) = & 
                 &     (fpsi(ix,iy,2)-fpsi(ix,iy,0)+fpsi(ix,iy,3)-fpsi(ix,iy,1))

            !       psdx=0.25*(psidx(ix,iy,0)+psidx(ix,iy,1)+psidx(ix,iy,2)
            !    .    +psidx(ix,iy,3))
            !       psdy=0.25*(psidy(ix,iy,0)+psidy(ix,iy,1)+psidy(ix,iy,2)
            !    .    +psidy(ix,iy,3))
            !       babs=sqrt(psdx*psdx+psdy*psdy)/t0
!WG            babs= averageVertex( &
!WG             & sqrt(psidx(ix,iy,0)**2+psidy(ix,iy,0)**2),crx(ix,iy,0),cry(ix,iy,0), &
!WG             & sqrt(psidx(ix,iy,1)**2+psidy(ix,iy,1)**2),crx(ix,iy,1),cry(ix,iy,1), &
!WG             & sqrt(psidx(ix,iy,2)**2+psidy(ix,iy,2)**2),crx(ix,iy,2),cry(ix,iy,2), &
!WG             & sqrt(psidx(ix,iy,3)**2+psidy(ix,iy,3)**2),crx(ix,iy,2),cry(ix,iy,3) )/t0
!WG         WG-TMP back to simple averaging for consistency later on (b2agmt)  
            babs = 0.25_R8*( &
     &        sqrt(psidx(ix,iy,0)**2+psidy(ix,iy,0)**2) + &
     &        sqrt(psidx(ix,iy,1)**2+psidy(ix,iy,1)**2) + &
     &        sqrt(psidx(ix,iy,2)**2+psidy(ix,iy,2)**2) + &
     &        sqrt(psidx(ix,iy,3)**2+psidy(ix,iy,3)**2) )/t0
            bb(ix,iy,0)=sign(babs,bb(ix,iy,0))

            if (isymm.eq.1 .or. isymm.eq.2) then
              babs=sqrt(psidx(ix,iy,0)**2+psidy(ix,iy,0)**2)/ &
                 &    (2.0*pi*crx(ix,iy,0))
            elseif (isymm.eq.3 .or. isymm.eq.4) then
              babs=sqrt(psidx(ix,iy,0)**2+psidy(ix,iy,0)**2)/ &
                 &    (2.0*pi*cry(ix,iy,0))
            else
              babs=sqrt(psidx(ix,iy,0)**2+psidy(ix,iy,0)**2)
            endif
            wbbc(ix,iy,0)=sign(babs,bb(ix,iy,0))

            bb(ix,iy,1) = 0.0e0
            wbbc(ix,iy,1) = 0.0e0

!WG         WG-TMP back to simple averaging for consistency later on (b2agmt)  
            bb(ix,iy,2) = & 
     &            0.25_R8*(ffbz(ix,iy,0) + ffbz(ix,iy,1) + &
     &                     ffbz(ix,iy,2) + ffbz(ix,iy,3) ) / t0
!WG            bb(ix,iy,2) = & 
!WG     &       averageVertex(ffbz(ix,iy,0),crx(ix,iy,0),cry(ix,iy,0),  &
!WG     &                     ffbz(ix,iy,1),crx(ix,iy,1),cry(ix,iy,1),  &
!WG     &                     ffbz(ix,iy,2),crx(ix,iy,2),cry(ix,iy,2),  &
!WG     &                     ffbz(ix,iy,3),crx(ix,iy,3),cry(ix,iy,3)) / t0
            if (isymm.eq.1 .or. isymm.eq.2) then
              wbbc(ix,iy,2) = ffbz(ix,iy,0)/(2.0*pi*crx(ix,iy,0))
            elseif (isymm.eq.3 .or. isymm.eq.4) then
              wbbc(ix,iy,2) = ffbz(ix,iy,0)/(2.0*pi*cry(ix,iy,0))
            else
              wbbc(ix,iy,2) = ffbz(ix,iy,0)
            endif
            bb(ix,iy,3) = & 
                 &     sqrt(bb(ix,iy,0)**2+bb(ix,iy,1)**2+bb(ix,iy,2)**2)
            wbbc(ix,iy,3) = & 
                 &     sqrt(wbbc(ix,iy,0)**2+wbbc(ix,iy,1)**2+wbbc(ix,iy,2)**2)

            if (ix.eq.nx .and. iy.lt.ny) then
              if (isymm.eq.1 .or. isymm.eq.2) then
                babs=sqrt(psidx(nx,iy,1)**2+psidy(nx,iy,1)**2)/ &
                 &    (2.0*pi*crx(nx,iy,1))
              elseif (isymm.eq.3 .or. isymm.eq.4) then
                babs=sqrt(psidx(nx,iy,1)**2+psidy(nx,iy,1)**2)/ &
                 &    (2.0*pi*cry(nx,iy,1))
              else
                babs=sqrt(psidx(nx,iy,1)**2+psidy(nx,iy,1)**2)
              endif
              wbbc(nx+1,iy,0) = sign(babs,wbbc(nx,iy,0))
              wbbc(nx+1,iy,1) = 0.0e0
              if (isymm.eq.1 .or. isymm.eq.2) then
                wbbc(nx+1,iy,2) = ffbz(nx,iy,1)/(2.0*pi*crx(nx,iy,1))
              elseif (isymm.eq.3 .or. isymm.eq.4) then
                wbbc(nx+1,iy,2) = ffbz(nx,iy,1)/(2.0*pi*cry(nx,iy,1))
              else
                wbbc(nx+1,iy,2) = ffbz(nx,iy,1)
              endif
              wbbc(nx+1,iy,3) = & 
                 &     sqrt(wbbc(nx+1,iy,0)**2+wbbc(nx+1,iy,1)**2+ &
                 &          wbbc(nx+1,iy,2)**2)
            else if (ix.lt.nx .and. iy.eq.ny) then
              if (isymm.eq.1 .or. isymm.eq.2) then
                babs=sqrt(psidx(ix,ny,2)**2+psidy(ix,ny,2)**2)/ &
                 &    (2.0*pi*crx(ix,ny,2))
              elseif (isymm.eq.3 .or. isymm.eq.4) then
                babs=sqrt(psidx(ix,ny,2)**2+psidy(ix,ny,2)**2)/ &
                 &    (2.0*pi*cry(ix,ny,2))
              else
                babs=sqrt(psidx(ix,ny,2)**2+psidy(ix,ny,2)**2)
              endif
              wbbc(ix,ny+1,0) = sign(babs,wbbc(ix,ny,0))
              wbbc(ix,ny+1,1) = 0.0e0
              if (isymm.eq.1 .or. isymm.eq.2) then
                wbbc(ix,ny+1,2) = ffbz(ix,ny,2)/(2.0*pi*crx(ix,ny,2))
              elseif (isymm.eq.3 .or. isymm.eq.4) then
                wbbc(ix,ny+1,2) = ffbz(ix,ny,2)/(2.0*pi*cry(ix,ny,2))
              else
                wbbc(ix,ny+1,2) = ffbz(ix,ny,2)
              endif
              wbbc(ix,ny+1,3) = & 
                 &     sqrt(wbbc(ix,ny+1,0)**2+wbbc(ix,ny+1,1)**2+ &
                 &          wbbc(ix,ny+1,2)**2)
            else if (ix.eq.nx .and. iy.eq.ny) then
              if (isymm.eq.1 .or. isymm.eq.2) then
                babs=sqrt(psidx(nx,ny,3)**2+psidy(nx,ny,3)**2)/ &
                 &    (2.0*pi*crx(nx,ny,3))
              elseif (isymm.eq.3 .or. isymm.eq.4) then
                babs=sqrt(psidx(nx,ny,3)**2+psidy(nx,ny,3)**2)/ &
                 &    (2.0*pi*cry(nx,ny,3))
              else
                babs=sqrt(psidx(nx,ny,3)**2+psidy(nx,ny,3)**2)
              endif
              wbbc(nx+1,ny+1,0) = sign(babs,wbbc(nx,ny,0))
              wbbc(nx+1,ny+1,1) = 0.0e0
              if (isymm.eq.1 .or. isymm.eq.2) then
                wbbc(nx+1,ny+1,2) = ffbz(nx,ny,3)/(2.0*pi*crx(nx,ny,3))
              elseif (isymm.eq.3 .or. isymm.eq.4) then
                wbbc(nx+1,ny+1,2) = ffbz(nx,ny,3)/(2.0*pi*cry(nx,ny,3))
              else
                wbbc(nx+1,ny+1,2) = ffbz(nx,ny,3)
              endif
              wbbc(nx+1,ny+1,3) = & 
                 &     sqrt(wbbc(nx+1,ny+1,0)**2+wbbc(nx+1,ny+1,1)**2+ &
                 &          wbbc(nx+1,ny+1,2)**2)
            end if

      enddo
    enddo

  end subroutine carre_b2agbb

end module carre_b2ag
