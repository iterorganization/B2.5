c
      subroutine mapy(nx, ny, ns, copy, field, work)
      use b2mod_types
      use b2mod_indirect
	logical copy
      integer nx, ny, ns
      real (kind=R8) ::
     1 field(-1:nx,-1:ny,0:ns-1), work(-1:nx,-1:ny,0:ns-1)
      integer ix, iy, is

      do is=0,ns-1
         do iy=-1,ny
           do ix=-1,nx
             if(iy.eq.-1) then
               work(ix,iy,is)=field(ix,iy+1,is)
             elseif(iy.eq.ny) then
               work(ix,iy,is)=field(ix,iy,is)
             else
               work(ix,iy,is)=0.5_R8*
     1	       (field(topix(ix,iy),topiy(ix,iy),is)+field(ix,iy,is))
             endif
           enddo
         enddo
      enddo
	if(copy)then                                                   !srv 07.07.09
       do is=0,ns-1
         do iy=-1,ny
           do ix=-1,nx
             field(ix,iy,is)=work(ix,iy,is)
           enddo
         enddo
       enddo
	endif
      return
      end
