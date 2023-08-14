! $AD NOCHECKPOINT
      subroutine intvertex1(nCv,nVx,nz,mpg,vxvol,centre,vertex)
      use b2mod_types
      use b2us_map
      implicit none
      type (mapping), intent (in) :: mpg
      integer nCv,nVx,nz
      real (kind=R8) ::
     *   vxvol(mpg%nVmxCv),centre(nCv,0:nz),vertex(nVx,0:nz)
      integer is
      external intvertex
      do is=0,nz
        call intvertex(nCv,nVx,mpg,vxvol,centre(1,is),vertex(1,is))
      enddo
      return
      end subroutine intvertex1
c
! $AD NOCHECKPOINT
      subroutine intvertex(nCv,nVx,mpg,vxvol,centre,vertex)
      use b2mod_types
      use b2us_map
      implicit none
      integer nCv,nVx
      type (mapping), intent (in) :: mpg
      real (kind=R8) ::
     *   vxvol(mpg%nVmxCv),centre(nCv),vertex(nVx)
  
      integer iVx,iCv
      real (kind=R8) :: volsum

      vertex = 0.0_R8
      do iVx = 1, nVx
        volsum = 0.0_R8
        do iCv = 1,mpg%vxCvP(iVx,2)
          volsum = volsum + 1._R8/vxvol(mpg%vxCvP(iVx,1) + iCv - 1)
          vertex(iVx) = vertex(iVx) + 
     &      centre(mpg%vxCv(mpg%vxCvP(iVx,1) + iCv - 1))/
     &      vxvol(mpg%vxCvP(iVx,1) + iCv - 1)
        enddo
        vertex(iVx) = vertex(iVx)/(volsum + 1.0e-30_R8)
      enddo

      return
      end

      function intvertex_s(iVx,nCv,nVx,mpg,vxvol,centre) result(vertex)
      !same as intevertex, but interpolation to a single vertex iVx
      use b2mod_types
      use b2us_map
      implicit none
      integer nCv,nVx,iVx
      type (mapping), intent (in) :: mpg
      real (kind=R8), intent(in) ::
     *   vxvol(mpg%nVmxCv),centre(nCv)
      real (kind=R8) :: vertex
  
      integer iCv
      real (kind=R8) :: volsum

      vertex = 0.0_R8
      volsum = 0.0_R8
      do iCv = 1,mpg%vxCvP(iVx,2)
        volsum = volsum + 1._R8/vxvol(mpg%vxCvP(iVx,1) + iCv - 1)
        vertex = vertex + 
     &    centre(mpg%vxCv(mpg%vxCvP(iVx,1) + iCv - 1))/
     &    vxvol(mpg%vxCvP(iVx,1) + iCv - 1)
      enddo
      vertex = vertex/(volsum + 1.0e-30_R8)

      return
      end function intvertex_s

