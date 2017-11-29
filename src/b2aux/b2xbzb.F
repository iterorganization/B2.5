      subroutine b2xbzb (nx, ny, bb, vol, bzb)              !srv 13.10.06 23.07.10 {
      use b2mod_types
      use b2mod_indirect
      implicit none
      integer nx, ny
      real (kind=R8) ::
     *   bb(-1:nx,-1:ny,0:3), vol(-1:nx,-1:ny), bzb(-1:nx,-1:ny)
*     ------------------------------------------------------------------
*     B2XBZB calculates Bz/B**2 - Bz/<B**2>.
*     ------------------------------------------------------------------
      integer ix, iy
      real (kind=R8) ::
     *   bbsum, volsum
*     ------------------------------------------------------------------
*$$$  call subini ('b2xbzb')

      bbsum  = 0.0e0_R8
      volsum = 0.0e0_R8
      do iy = -1,ny
       do ix = -1,nx
        if(on_closed_surface(ix,iy)) then      ! for CORE
         bbsum = bbsum + bb(ix,iy,3)**2*vol(ix,iy)
         volsum = volsum + vol(ix,iy)
        endif
       enddo
      enddo
      if(volsum.gt.0.0_R8) then
        bzb = bbsum/volsum
        bzb = bb(:,:,2)/bb(:,:,3)**2 - bb(:,:,2)/bzb
      else
        bzb = 0.0_R8
      endif

*$$$  call subend ()
      return
*     ------------------------------------------------------------------
      end                                                   !srv 13.10.06 23.07.10 }
