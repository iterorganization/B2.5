      module b2mod_geo
      use b2mod_types
      use b2mod_version
      use b2mod_geo_corner
      implicit none

      real (kind=R8), allocatable, save ::
     *  bb(:,:,:),
     *  crx(:,:,:),
     *  cry(:,:,:),
     *  ffbz(:,:,:),
     *  fpsi(:,:,:),
     *  gs(:,:,:),
     *  hx(:,:),
     *  hy(:,:),
     *  hz(:,:),                                                     !srv 09.01.01
     *  qz(:,:,:),
     *  qc(:,:),
     *  vol(:,:),
     *  pbs(:,:,:), pbshz(:,:,:),                                    !srv 09.01.01
     *  conn(:)

      contains

      subroutine alloc_b2mod_geo(nxd,nyd)
      implicit none
      integer nxd,nyd

      allocate(bb(-1:nxd,-1:nyd,0:3))
      allocate(crx(-1:nxd,-1:nyd,0:3))
      allocate(cry(-1:nxd,-1:nyd,0:3))
      allocate(ffbz(-1:nxd,-1:nyd,0:3))
      allocate(fpsi(-1:nxd,-1:nyd,0:3))
      allocate(gs(-1:nxd,-1:nyd,0:2))
      allocate(hx(-1:nxd,-1:nyd))
      allocate(hy(-1:nxd,-1:nyd))
      allocate(hz(-1:nxd,-1:nyd))                                    !srv 09.01.01
      allocate(qz(-1:nxd,-1:nyd,0:1))
      allocate(qc(-1:nxd,-1:nyd))
      allocate(vol(-1:nxd,-1:nyd))
      allocate(pbs(-1:nxd,-1:nyd,0:1))
      allocate(pbshz(-1:nxd,-1:nyd,0:1))                             !srv 09.01.01
      allocate(conn(-1:nyd))

      call alloc_b2mod_geo_corner(nxd,nyd)

      return
      end subroutine alloc_b2mod_geo

      subroutine dealloc_b2mod_geo
      implicit none

      deallocate(bb)
      deallocate(crx)
      deallocate(cry)
      deallocate(ffbz)
      deallocate(fpsi)
      deallocate(gs)
      deallocate(hx)
      deallocate(hy)
      deallocate(hz)
      deallocate(qz)
      deallocate(qc)
      deallocate(vol)
      deallocate(pbs)
      deallocate(pbshz)                                              !srv 09.01.01
      deallocate(conn)

      call dealloc_b2mod_geo_corner

      return
      end subroutine dealloc_b2mod_geo

      subroutine read_b2mod_geo(nx,ny,nread)
      implicit none
      integer nread, nx, ny
      integer ix, iy
      real (kind=R8) ::
     *   p0x, p0y, p1x, p1y
*   ..procedures
      intrinsic min, max, sqrt
      real (kind=R8) ::
     *   trim1, norm, cosang
      norm(p0x,p0y) = sqrt(p0x**2+p0y**2)
      trim1(p0x) = max(-1.0_R8,min(1.0_R8,p0x))
      cosang(p0x,p0y,p1x,p1y) = trim1(
     &  (p0x*p1y-p0y*p1x)/(norm(p0x,p0y)*norm(p1x,p1y)))

      call cfrure (nread, size(bb), bb, 'bb')
      call cfrure (nread, size(crx), crx, 'crx')
      call cfrure (nread, size(cry), cry, 'cry')
      call cfrure (nread, size(ffbz), ffbz, 'ffbz')
      call cfrure (nread, size(fpsi), fpsi, 'fpsi')
      call cfrure (nread, size(gs), gs, 'gs')
      call cfrure (nread, size(hx), hx, 'hx')
      call cfrure (nread, size(hy), hy, 'hy')
      call cfrure (nread, size(qz), qz, 'qz')
      if (b2fplasma_version.gt.'01.001.021') then
       call cfrure (nread, size(qc), qc, 'qc')
      else
       do iy = -1, ny
        do ix = -1, nx
         qc(ix,iy) = cosang (
     &    (crx(ix,iy,1)+crx(ix,iy,3))-(crx(ix,iy,0)+crx(ix,iy,2)),
     &    (cry(ix,iy,1)+cry(ix,iy,3))-(cry(ix,iy,0)+cry(ix,iy,2)),
     &    crx(ix,iy,2)-crx(ix,iy,0), cry(ix,iy,2)-cry(ix,iy,0))
         if(qc(ix,iy).le.0.0_R8) then
          write(*,*) 'Problem in b2mod_geo: qc(ix,iy) <= 0',
     1     qc(ix,iy),ix,iy
          write(*,'(a,4f12.6)') 'crx = ',crx(ix,iy,0:3)
          write(*,'(a,4f12.6)') 'cry = ',cry(ix,iy,0:3)
         endif
        enddo
       enddo
      endif
      call cfrure (nread, size(vol), vol, 'vol')
c conn
      if(b2fplasma_version.ge.'01.001.063') then        !srv 15.02.11
        call cfrure (nread, size(pbs), pbs, 'pbs')
      endif

      return
      end subroutine read_b2mod_geo

      subroutine write_b2mod_geo(nwrite)
      implicit none
      integer nwrite

      call cfwure (nwrite, size(bb), bb, 'bb')
      call cfwure (nwrite, size(crx), crx, 'crx')
      call cfwure (nwrite, size(cry), cry, 'cry')
      call cfwure (nwrite, size(ffbz), ffbz, 'ffbz')
      call cfwure (nwrite, size(fpsi), fpsi, 'fpsi')
      call cfwure (nwrite, size(gs), gs, 'gs')
      call cfwure (nwrite, size(hx), hx, 'hx')
      call cfwure (nwrite, size(hy), hy, 'hy')
      call cfwure (nwrite, size(qz), qz, 'qz')
      if(newversion.gt.'01.001.021') then
        call cfwure (nwrite, size(qc), qc, 'qc')
      endif
      call cfwure (nwrite, size(vol), vol, 'vol')
c conn
      if(newversion.ge.'01.001.063') then               !srv 15.02.11
        call cfwure (nwrite, size(pbs), pbs, 'pbs')
      endif
      return
      end subroutine write_b2mod_geo

      end module b2mod_geo
