      subroutine intface_us1(nCv,nFc,nz,fcCv,vol,centre,face)
      use b2mod_types
      implicit none
      integer nCv,nFc,nz
      integer fcCv(nFc,2)
      real (kind=R8) ::
     *   vol(nCv),centre(nCv,0:nz),face(nFc,0:nz)
      integer is
      do is=0,nz
        call intface_us(nCv,nFc,fcCv,vol,centre(1,is),face(1,is))
      enddo
      return
      end subroutine intface_us1
c
      subroutine intface_us(nCv,nFc,fcCv,vol,centre,face)
      use b2mod_types
      implicit none

      integer nCv,nFc
      integer fcCv(nFc,2)
      real (kind=R8) ::
     *   vol(nCv),centre(nCv),face(nFc)
      
      integer iFc
      do iFc = 1, nFc
        face(iFc) =
     &    (vol(fcCv(iFc,1))*centre(fcCv(iFc,2))+
     &     vol(fcCv(iFc,2))*centre(fcCv(iFc,1)))/
     &    (vol(fcCv(iFc,1))+vol(fcCv(iFc,2)))
      enddo
!wdk   to be reactivated after adaptation of b2xxgs...
!      call b2xxgs (nx, ny, 0.0_R8, face, 1)
     
      return
      end subroutine intface_us
