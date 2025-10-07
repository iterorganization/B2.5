      subroutine intflux1(nFc,nCv,nz,mpg,geo,faceflux,density,cvflux)
      use b2mod_types
      use b2us_map
      use b2us_geo
      implicit none
      integer, intent (in) :: nCv,nFc,nz
      type (mapping), intent (in) :: mpg
      type (geometry), intent (in) :: geo
      real (kind=R8), intent (in) :: faceflux(nFc,0:1,0:nz-1)
      real (kind=R8), intent (in) :: density(nCv,0:nz-1)
      real (kind=R8), intent (out) :: cvflux(nCv,0:1,0:nz-1)
      integer is
      do is=0,nz-1
        call intflux(nFc,nCv,mpg,geo,faceflux(1,0,is),density(1,is),
     .   cvflux(1,0,is))
      enddo
      return
      end subroutine intflux1
c
! $AD NOCHECKPOINT
      subroutine intflux(nFc,nCv,mpg,geo,faceflux,density,cvflux)
c Provides an interpolation of face flux quantities to cell centers
c The flux carries the cell-centered density of the quantity being transported
c We first interpolate the density to cell faces
c We then compute an effective velocity at cell faces
C We interpolate this velocity to cell centers
c We multiply the interpolated velocity by the density and the projected area
      use b2mod_types
      use b2us_map
      use b2us_geo
      implicit none
      integer, intent (in) :: nCv,nFc
      type (mapping), intent (in) :: mpg
      type (geometry), intent (in) :: geo
      real (kind=R8), intent (in) :: faceflux(nFc,0:1)
      real (kind=R8), intent (in) :: density(nCv)
      real (kind=R8), intent (out) :: cvflux(nCv,0:1)

      integer idir,iFc
      real (kind=R8) :: qf(nFc), wrkf(nFc,0:1), wrkc(nCv,0:1)

c ..interpolate the quantity to cell faces
      call intface(nCv,nFc,mpg%fcCv,geo%fcVol,density,qf)

c ..compute the effective velocity
      do idir = 0, 1
        do iFc = 1, nFc
          if (geo%fcQalf(iFc,idir).ne.0.0_R8.and.qf(iFc).ne.0.0_R8) then
            wrkf(iFc,idir) = faceflux(iFc,idir)/qf(iFc)/
     .                      (geo%fcS(iFc)*abs(geo%fcQalf(iFc,idir)))
          else
            wrkf(iFc,idir) = 0.0_R8
          end if
        end do
      end do

c ..interpolate the effective velocity to cell centers
      call intcell(nFc,nCv,mpg,mpg%intcellP,wrkf(:,0),wrkc(:,0))
      call intcell(nFc,nCv,mpg,mpg%intcellR,wrkf(:,1),wrkc(:,1))

c ..multiply the interpolated velocity
c ..by the projected cell areas and the cell-centered quantity density
      cvflux(:,0) = wrkc(:,0)*density(:)*geo%cvHy(:)*geo%cvHz(:)
      cvflux(:,1) = wrkc(:,1)*density(:)*geo%cvHx(:)*geo%cvHz(:)

      return
      end subroutine intflux
