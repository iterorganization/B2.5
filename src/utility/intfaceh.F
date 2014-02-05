      subroutine intfaceh1(nx,ny,nz,leftix,leftiy,vol,centre,leftface)
      use b2mod_types
      implicit none
      integer nx,ny,nz
      integer leftix(-1:nx,-1:ny), leftiy(-1:nx,-1:ny)
      real (kind=R8) ::
     *   vol(-1:nx,-1:ny),centre(-1:nx,-1:ny,0:nz),
     1 leftface(-1:nx,-1:ny,0:nz)
      integer is
      do is=0,nz
        call intfaceh(nx,ny,leftix,leftiy,vol,centre(-1,-1,is),
     .   leftface(-1,-1,is))
      enddo
      return
      end subroutine intfaceh1
c
      subroutine intfaceh(nx,ny,leftix,leftiy,vol,centre,leftface)
      use b2mod_types
      implicit none

      integer nx,ny
      integer leftix(-1:nx,-1:ny), leftiy(-1:nx,-1:ny)
      real (kind=R8) ::
     *   vol(-1:nx,-1:ny),centre(-1:nx,-1:ny),leftface(-1:nx,-1:ny)
      
      integer ix,iy
      do iy = -1, ny
       do ix = -1, nx
        if(leftix(ix,iy).eq.-2) then
          leftface(ix,iy) = 0.0_R8
        else
          leftface(ix,iy) =
     &     (vol(ix,iy)*centre(leftix(ix,iy),leftiy(ix,iy))+
     &      vol(leftix(ix,iy),leftiy(ix,iy))*centre(ix,iy))/
     &     (vol(leftix(ix,iy),leftiy(ix,iy))+vol(ix,iy))
        endif
       enddo
      enddo
      call b2xxgs (nx, ny, 0.0_R8, leftface, 1)
     
      return
      end subroutine intfaceh
c
      subroutine intfacer1(nx,ny,nz,
     . rightix,rightiy,vol,centre,rightface)
      use b2mod_types
      implicit none
      integer nx,ny,nz
      integer rightix(-1:nx,-1:ny), rightiy(-1:nx,-1:ny)
      real (kind=R8) ::
     *   vol(-1:nx,-1:ny),centre(-1:nx,-1:ny,0:nz),
     1 rightface(-1:nx,-1:ny,0:nz)
      integer is
      do is=0,nz
        call intfacer(nx,ny,rightix,rightiy,vol,centre(-1,-1,is),
     .   rightface(-1,-1,is))
      enddo
      return
      end subroutine intfacer1
c
      subroutine intfacer(nx,ny,rightix,rightiy,vol,centre,rightface)
      use b2mod_types
      implicit none

      integer nx,ny
      integer rightix(-1:nx,-1:ny), rightiy(-1:nx,-1:ny)
      real (kind=R8) ::
     *   vol(-1:nx,-1:ny),centre(-1:nx,-1:ny),rightface(-1:nx,-1:ny)
      
      integer ix,iy
      do iy = -1, ny
       do ix = -1, nx
        if(rightix(ix,iy).eq.nx+1) then
          rightface(ix,iy) = 0.0_R8
        else
          rightface(ix,iy) =
     &     (vol(ix,iy)*centre(rightix(ix,iy),rightiy(ix,iy))+
     &      vol(rightix(ix,iy),rightiy(ix,iy))*centre(ix,iy))/
     &     (vol(rightix(ix,iy),rightiy(ix,iy))+vol(ix,iy))
        endif
       enddo
      enddo
      call b2xxgs (nx, ny, 0.0_R8, rightface, 1)
     
      return
      end subroutine intfacer
