      subroutine intfaceu1(nCv,nFc,nz,fcCv,flo,centre,face)
      use b2mod_types
      implicit none
      integer nCv,nFc,nz
      integer fcCv(nFc,2)
      real (kind=R8) ::
     *   flo(nFc),centre(nCv,0:nz),face(nFc,0:nz)
      integer is
      do is=0,nz
        call intfaceu(nCv,nFc,fcCv,flo,centre(1,is),face(1,is))
      enddo
      return
      end subroutine intfaceu1
c
      subroutine intfaceu(nCv,nFc,fcCv,flo,centre,face)
      use b2mod_types
      implicit none

      integer nCv,nFc
      integer fcCv(nFc,2)
      real (kind=R8) ::
     *   flo(nFc),centre(nCv),face(nFc)
      integer iFc

      do iFc = 1, nFc
        if (flo(iFc).gt.0.0_R8) face(iFc) = centre(fcCv(iFc,1))
        if (flo(iFc).lt.0.0_R8) face(iFc) = centre(fcCv(iFc,2))
        if (flo(iFc).eq.0.0_R8) 
     &    face(iFc) = 0.5_R8*(centre(fcCv(iFc,1)) + centre(fcCv(iFc,2)))
      enddo
     
      return
      end subroutine intfaceu
c
      function intfaceu_s(iFc,nCv,nFc,fcCv,flo,centre) result(face)
      use b2mod_types
      implicit none

      integer, intent(in):: nCv,nFc, iFc
      integer, intent(in) :: fcCv(nFc,2)
      real (kind=R8), intent(in) ::
     *   flo(nFc),centre(nCv)
      real (kind=R8) :: face

      if (flo(iFc).gt.0.0_R8) face = centre(fcCv(iFc,1))
      if (flo(iFc).lt.0.0_R8) face = centre(fcCv(iFc,2))
      if (flo(iFc).eq.0.0_R8) 
     &  face = 0.5_R8*(centre(fcCv(iFc,1)) + centre(fcCv(iFc,2)))

      return
      end function intfaceu_s
