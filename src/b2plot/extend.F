
      subroutine extend(field,nxd,nyd,nsd,nx,ny,ns)
      use b2mod_types
      use b2mod_indirect
      integer nxd,nyd,nsd,nx,ny,ns
      real (kind=R8) ::
     *   field(-1:nxd,-1:nyd,0:nsd-1)
      integer ix,iy,is
      integer ixn,iyn
      do is=0,ns-1
        do iy=-1,ny
          do ix=-1,nx
            ixn=ix
            iyn=iy
            if (bottomiy(ix,iy).eq.-2 .and. leftix(ix,iy).eq.-2) then
              ixn=rightix(topix(ix,iy),topiy(ix,iy))
              iyn=rightiy(topix(ix,iy),topiy(ix,iy))
            elseif (bottomiy(ix,iy).eq.-2 .and. rightix(ix,iy).eq.nx+1) 
     >       then
              ixn=leftix(topix(ix,iy),topiy(ix,iy))
              iyn=leftiy(topix(ix,iy),topiy(ix,iy))
            elseif (bottomiy(ix,iy).eq.-2) then
              ixn=topix(ix,iy)
              iyn=topiy(ix,iy)
            elseif (topiy(ix,iy).eq.ny+1 .and. leftix(ix,iy).eq.-2) then
              ixn=rightix(bottomix(ix,iy),bottomiy(ix,iy))
              iyn=rightiy(bottomix(ix,iy),bottomiy(ix,iy))
            elseif (topiy(ix,iy).eq.ny+1 .and. rightix(ix,iy).eq.nx+1) 
     >       then
              ixn=leftix(bottomix(ix,iy),bottomiy(ix,iy))
              iyn=leftiy(bottomix(ix,iy),bottomiy(ix,iy))
            elseif (topiy(ix,iy).eq.ny+1) then
              ixn=bottomix(ix,iy)
              iyn=bottomiy(ix,iy)
            elseif (leftix(ix,iy).eq.-2) then
              ixn=rightix(ix,iy)
              iyn=rightiy(ix,iy)
            elseif (rightix(ix,iy).eq.nx+1) then
              ixn=leftix(ix,iy)
              iyn=leftiy(ix,iy)
            endif
            if (ixn.ne.ix .or. iyn.ne.iy) 
     >        field(ix,iy,is)=field(ixn,iyn,is)
          enddo
        enddo
      enddo
      return
      end
