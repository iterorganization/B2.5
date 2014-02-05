
      subroutine extend(field,nxd,nyd,nsd,nx,ny,ns)
      use b2mod_types
      integer nxd,nyd,nsd,nx,ny,ns
      real (kind=R8) ::
     *   field(-1:nxd,-1:nyd,0:nsd-1)
      integer ix,iy,is
      do is=0,ns-1
        do iy=0,ny-1
          field(-1,iy,is)=field(0,iy,is)
          field(nx,iy,is)=field(nx-1,iy,is)
        enddo
        do ix=-1,nx
          field(ix,-1,is)=field(ix,0,is)
          field(ix,ny,is)=field(ix,ny-1,is)
        enddo
      enddo
      return
      end
