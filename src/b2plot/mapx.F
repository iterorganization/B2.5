c
      subroutine mapx(nx,ny,ns,field,work)
      use b2mod_types
      use b2mod_indirect
      implicit none
      integer nx,ny,ns,ix,iy,is
      real (kind=R8) ::
     1 field(-1:nx,-1:ny,0:ns-1),work(-1:nx,-1:ny,0:ns-1)

      do is=0,ns-1
        do iy=-1,ny
          do ix=-1,nx
            if(leftix(ix,iy).eq.-2) then
              work(ix,iy,is)=field(rightix(ix,iy),rightiy(ix,iy),is)
            elseif(rightix(ix,iy).eq.nx+1) then
              work(ix,iy,is)=field(ix,iy,is)
            else
              work(ix,iy,is)=0.5_R8*
     1         (field(rightix(ix,iy),rightiy(ix,iy),is)+field(ix,iy,is))
            endif
          enddo
        enddo
      enddo
      do is=0,ns-1
        do iy=-1,ny
          do ix=-1,nx
            field(ix,iy,is)=work(ix,iy,is)
          enddo
        enddo
      enddo
      return
      end
