      subroutine intface1(nCv,nFc,nz,fcCv,fcVol,centre,face)
      use b2mod_types
      implicit none
      integer nCv,nFc,nz
      integer fcCv(nFc,2)
      real (kind=R8) ::
     *   fcVol(nFc,2),centre(nCv,0:nz),face(nFc,0:nz)
      integer is
      external intface
      do is=0,nz
        call intface(nCv,nFc,fcCv,fcVol,centre(1,is),face(1,is))
      enddo
      return
      end subroutine intface1
c
! $AD NOCHECKPOINT
      subroutine intface(nCv,nFc,fcCv,fcVol,centre,face)
      use b2mod_types
      implicit none

      integer, intent(in) :: nCv,nFc
      integer, intent(in) :: fcCv(nFc,2)
      real (kind=R8), intent(in) ::
     *   fcVol(nFc,2),centre(nCv)
      real(kind=R8), intent(out) :: face(nFc)
      integer iFc
      do iFc = 1, nFc
        face(iFc) =
     &    (fcVol(iFc,1)*centre(fcCv(iFc,2))+
     &     fcVol(iFc,2)*centre(fcCv(iFc,1)))/
     &    (fcVol(iFc,1)+fcVol(iFc,2))
      enddo
     
      return
      end subroutine intface
c
      function intface_s(iFc,nCv,nFc,fcCv,fcVol,centre) result(face)
      !same as intface, but interpolation to a single face iFc
      use b2mod_types
      implicit none

      integer, intent(in):: nCv,nFc, iFc
      integer, intent(in) :: fcCv(nFc,2)
      real (kind=R8), intent(in) ::
     *   fcVol(nFc,2),centre(nCv)
      real (kind=R8) :: face

        face =
     &    (fcVol(iFc,1)*centre(fcCv(iFc,2))+
     &     fcVol(iFc,2)*centre(fcCv(iFc,1)))/
     &    (fcVol(iFc,1)+fcVol(iFc,2))
     
      return
      end function intface_s
