      subroutine intfacev1(nx,ny,nz,bottomix,bottomiy,vol,centre,face)
      use b2mod_types
      implicit none
      integer nx,ny,nz
      integer bottomix(-1:nx,-1:ny), bottomiy(-1:nx,-1:ny)
      real (kind=R8) ::
     *   vol(-1:nx,-1:ny),centre(-1:nx,-1:ny,0:nz),
     1 face(-1:nx,-1:ny,0:nz)
      integer is
      do is=0,nz
        call intfacev(nx,ny,bottomix,bottomiy,vol,centre(-1,-1,is),
     .   face(-1,-1,is))
      enddo
      return
      end subroutine intfacev1
c
      subroutine intfacev(nx,ny,bottomix,bottomiy,vol,centre,face)
      use b2mod_types
      implicit none

      integer nx,ny
      integer bottomix(-1:nx,-1:ny), bottomiy(-1:nx,-1:ny)
      real (kind=R8) ::
     *   vol(-1:nx,-1:ny),centre(-1:nx,-1:ny),face(-1:nx,-1:ny)
      
      integer ix,iy

      do ix = -1, nx
        do iy = -1, ny
          if (bottomiy(ix,iy).eq.-2) then
            face(ix,iy)=0.
          else
            face(ix,iy) =
     &       (vol(ix,iy)*centre(bottomix(ix,iy),bottomiy(ix,iy))+
     &        vol(bottomix(ix,iy),bottomiy(ix,iy))*centre(ix,iy))/
     &       (vol(bottomix(ix,iy),bottomiy(ix,iy))+vol(ix,iy))
          endif
        enddo
      enddo
     
      return
      end subroutine intfacev
