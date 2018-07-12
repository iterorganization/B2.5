      subroutine intvertex1(nCv,nVx,nz,mpg,vol,centre,vertex)
      use b2mod_types
      use b2us_map
      implicit none
      type (mapping), intent (in) :: mpg
      integer nCv,nVx,nz
      real (kind=R8) ::
     *   vol(nCv),centre(nCv,0:nz),vertex(nVx,0:nz)
      integer is
      do is=0,nz
        call intvertex(nCv,nVx,mpg,vol,centre(1,is),vertex(1,is))
      enddo
      return
      end subroutine intvertex1
c
      subroutine intvertex(nCv,nVx,mpg,vol,centre,vertex)
      use b2mod_types
      use b2us_map
      implicit none
      integer nCv,nVx
      type (mapping), intent (in) :: mpg
      real (kind=R8) ::
     *   vol(nCv),centre(nCv),vertex(nVx)
  
      integer iVx,iCv
      real (kind=R8) :: volsum

      vertex=0.0_R8
      do iVx = 1, nVx
        volsum = 0.0_R8
        do iCv = 1,mpg%vxCvP(iVx,2)
          volsum = volsum + 1/vol(mpg%vxCv(mpg%vxCvP(iVx,1) + iCv - 1))
          vertex(iVx) = vertex(iVx) + 
     &      centre(mpg%vxCv(mpg%vxCvP(iVx,1) + iCv - 1))/
     &      vol(mpg%vxCv(mpg%vxCvP(iVx,1) + iCv - 1))
        enddo
        vertex(iVx) = vertex(iVx)/volsum
      enddo

      return
      end
