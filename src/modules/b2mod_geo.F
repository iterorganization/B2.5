      module b2mod_geo
      use b2mod_types
      use b2mod_version
      use b2mod_geo_corner
      use b2mod_cellhelper
      use b2mod_grid_mapping
      implicit none

      real (kind=R8), allocatable, save ::
     *  bb(:,:,:),
     *  crx(:,:,:),
     *  cry(:,:,:),
     *  ffbz(:,:,:),
     *  fpsi(:,:,:,:),
     *  gs(:,:,:),
     *  hx(:,:),
     *  hy(:,:),
     *  hz(:,:),                                                     !srv 09.01.01
     *  hc(:,:,:),
     *  ht(:,:,:),
     *  qz(:,:,:),
     *  qc(:,:,:),
     *  qs(:,:,:),
     *  qac(:,:,:),
     *  qas(:,:,:),
     *  vol(:,:),
     *  pbs(:,:,:), pbshz(:,:,:),                                    !srv 09.01.01
     *  ebc(:,:,:),
     *  conn(:),
     *  cell_width(:,:),
     *  cell_height(:,:)

      type(B2GridMap), save :: gmap

      logical, save :: b2mod_geo_allocated = .false.

            type, public :: geometry_ag
!  cell quantites
        real (kind=R8), allocatable  :: cvBb(:,:)      ! Magnetic field: (nCv,0:3) Bx, By, Bz, |B|
        real (kind=R8), allocatable  :: cvBp(:)        ! Magnetic field in poloidal direction: (nCv)
        real (kind=R8), allocatable  :: cvBt(:)        ! Magnetic field in toroidal direction: (nCv)
        real (kind=R8), allocatable  :: cvX(:)         ! X-co of cell center
        real (kind=R8), allocatable  :: cvY(:)         ! Y-co of cell center
        real (kind=R8), allocatable  :: cvSz(:)        ! Area perp. to third (z) direction
        real (kind=R8), allocatable  :: cvHz(:)        ! Length in third (z) direction
        real (kind=R8), allocatable  :: cvHx(:)        ! Temp: characteristic poloidal length; eventually remove, because meaning unclear in WG context
        real (kind=R8), allocatable  :: cvQgam(:,:)    ! Temp: complementary angle between poloidal and radial directions; remove?
        real (kind=R8), allocatable  :: cvVol(:)       ! Cell volume
        real (kind=R8), allocatable  :: cvOnedBsq(:)   ! 1/B**2 in cell centers
        real (kind=R8), allocatable  :: cvBzb(:)       ! Factor for P-S/mdf flows
        real (kind=R8), allocatable  :: cvEb(:,:)      ! Unit vector pointing along magnetic field (Cartesian X,Y,Z coordinates)
        integer (kind=R8), allocatable  :: cvCflags(:,:)  ! Cellflags 
        integer (kind=R8), allocatable  :: cvCtype(:)     ! Celltype
        real (kind=R8), allocatable  :: cvFpsi(:)      ! (nCv) psi at cell center

!  face based quantities
        real (kind=R8), allocatable  :: fcBb(:,:)      ! Magnetic field: (nFc,0:3) to be recomputed at cell faces (Bx,By,Bz,|B|)
        real (kind=R8), allocatable  :: fcS(:)         ! Total cell face area
        real (kind=R8), allocatable  :: fcHc(:,:)      ! Lengths between cell face center and neighboring cell centers ("connector" direction)
        real (kind=R8), allocatable  :: fcHt(:)        ! Lengths along cell face ("tangential" direction)
        real (kind=R8), allocatable  :: fcHz(:)        ! Length in third (z) direction
        real (kind=R8), allocatable  :: fcVol(:,:)     ! Part of volume of neighboring cell centers associated with the face; for interpolation
        real (kind=R8), allocatable  :: fcQgam(:,:)    ! Cosine (:,0) and sine (:,1) of gamma
        real (kind=R8), allocatable  :: fcQalf(:,:)    ! Cosine (:,0) and sine (:,1) of alpha
        real (kind=R8), allocatable  :: fcQbet(:,:)    ! Cosine (:,0) and sine (:,1) of beta
        real (kind=R8), allocatable  :: fcPbs(:)       ! Magnetic contact area
        real (kind=R8), allocatable  :: fcPbshz(:)     ! Magnetic contact area * hz (at faces)
        real (kind=R8), allocatable  :: fcBzb(:)       ! Factor for P-S/mdf flows
        real (kind=R8), allocatable  :: fcEb(:,:)      ! Unit vector pointing along magnetic field (Cartesian X,Y,Z coordinates)

!  vertex quantities
        real (kind=R8), allocatable  :: vxBb(:,:)      ! Magnetic field: (nVx,0:3) to be recomputed at vertices (Bx,By,Bz,|B|)
        real (kind=R8), allocatable  :: vxBx(:)        ! Magnetic field Bx = dpsi/dx
        real (kind=R8), allocatable  :: vxBy(:)        ! Magnetic field By = dpsi/dy
        real (kind=R8), allocatable  :: vxX(:)         ! X-co of vertex
        real (kind=R8), allocatable  :: vxY(:)         ! Y-co of vertex
        real (kind=R8), allocatable  :: vxHz(:)        ! Length in third (z) direction
        !real (kind=R8), allocatable  :: vxVol(:)       ! Part of volume of surrounding cell centers associated with the vertex; for interpolation
        real (kind=R8), allocatable  :: vxFfbz(:)      !
        real (kind=R8), allocatable  :: vxFpsi(:)      ! 
        real (kind=R8), allocatable  :: vxOnedBsq(:)   ! 1/B**2 at vertices
        real (kind=R8), allocatable  :: vxBzb(:)       ! Factor for P-S/mdf flows

!  flux tube/surface quantities
        real (kind=R8), allocatable  :: cvConn(:)      ! (nCv) Connection length for each CV
        !real (kind=R8), allocatable  :: ftConn(:)      ! (nFt) Connection length for each flux tube
        !real (kind=R8), allocatable  :: ftEps(:)       ! (nFt) Inverse aspect ratio
        !real (kind=R8), allocatable  :: ftBbav2(:)     ! (nFt) Average B**2 on the flux tube


      end type geometry_ag

      contains

*********************************************************************** 

      subroutine alloc_b2mod_geo(g,nCv,nFc,nVx)

      implicit none
      integer, intent(in) :: nCv, nFc, nVx
      type(geometry_ag), intent(out) :: g

      !cell quantities
      allocate (g%cvBb(0:nCv-1,0:3))
      allocate (g%cvBp(0:nCv-1))
      allocate (g%cvBt(0:nCv-1))
      allocate (g%cvX(0:nCv-1))
      allocate (g%cvY(0:nCv-1))
      allocate (g%cvSz(0:nCv-1))
      allocate (g%cvHz(0:nCv-1))
      allocate (g%cvHx(0:nCv-1))
      allocate (g%cvQgam(0:nCv-1,0:1))
      allocate (g%cvVol(0:nCv-1))
      allocate (g%cvOnedBsq(0:nCv-1))
      allocate (g%cvBzb(0:nCv-1))
      allocate (g%cvEb(0:nCv-1,0:2))
      allocate (g%cvCflags(0:nCv-1,0:4))
      allocate (g%cvCtype(0:nCv-1))
      allocate (g%cvFpsi(0:nCv-1))

      !face quantites
      allocate (g%fcBb(0:nFc-1,0:3)) 
      allocate (g%fcS(0:nFc-1))
      allocate (g%fcHc(0:nFc-1,1:2))
      allocate (g%fcHt(0:nFc-1))
      allocate (g%fcHz(0:nFc-1))
      allocate (g%fcVol(0:nFc-1,1:2))      
      allocate (g%fcQgam(0:nFc-1,0:1))
      allocate (g%fcQalf(0:nFc-1,0:1))
      allocate (g%fcQbet(0:nFc-1,0:1))
      allocate (g%fcPbs(0:nFc-1))
      allocate (g%fcPbshz(0:nFc-1))
      allocate (g%fcBzb(0:nFc-1))
      allocate (g%fcEb(0:nFc-1,0:2))

      !vertex quantites
      allocate (g%vxBb(0:nVx-1,0:3))
      allocate (g%vxBx(0:nVx-1))
      allocate (g%vxBy(0:nVx-1))
      allocate (g%vxX(0:nVx-1))
      allocate (g%vxY(0:nVx-1))
      allocate (g%vxHz(0:nVx-1))
      !allocate (g%vxVol(nVmxCv)) nVmxCv not known here
      allocate (g%vxFfbz(0:nVx-1))
      allocate (g%vxFpsi(0:nVx-1))
      allocate (g%vxOnedBsq(0:nVx-1))
      allocate (g%vxBzb(0:nVx-1))
      
      !flux tube quantites
      allocate (g%cvConn(0:nCv-1))
      !allocate (g%ftConn(0:nFt-1))
      !allocate (g%ftEps(0:nFt-1))
      !allocate (g%ftBbav2(0:nFt-1))

      !  flux surface quantities
      !allocate (g%fsPsi(0:nFs-1))

      call initialize_b2mod_geo (g)


      b2mod_geo_allocated = .true.

      return
      end subroutine alloc_b2mod_geo

***********************************************************************

      subroutine initialize_b2mod_geo (g)

      implicit none
      type(geometry_ag), intent(in out) :: g

      !  cell quantites
      g%cvBb = 0._R8
      g%cvBp = 0._R8
      g%cvBt = 0._R8
      g%cvX = 0._R8
      g%cvY = 0._R8
      g%cvSz = 0._R8
      g%cvHz = 0._R8
      g%cvHx = 0._R8
      g%cvQgam = 0._R8
      g%cvVol = 0._R8
      g%cvOnedBsq = 0._R8
      g%cvBzb = 0._R8
      g%cvEb = 0._R8
      g%cvCflags = 0
      g%cvCtype = 0
      g%cvFpsi = 0._R8

      !  face quantites
      g%fcBb = 0._R8
      g%fcS = 0._R8
      g%fcHc = 0._R8
      g%fcHt = 0._R8
      g%fcHz = 0._R8
      g%fcVol = 0._R8    
      g%fcQgam = 0._R8
      g%fcQalf = 0._R8
      g%fcQbet = 0._R8
      g%fcPbs = 0._R8
      g%fcPbshz = 0._R8
      g%fcBzb = 0._R8
      g%fcEb = 0._R8

      !  vertex quantites
      g%vxBb = 0._R8
      g%vxBx = 0._R8
      g%vxBy = 0._R8
      g%vxX = 0._R8
      g%vxY = 0._R8
      g%vxHz = 0._R8
      !g%vxVol = 0._R8
      g%vxFfbz = 0._R8
      g%vxFpsi = 0._R8
      g%vxOnedBsq = 0._R8
      g%vxBzb = 0._R8

      !  flux surface quantites
      g%cvConn = 0._R8
      !g%ftConn = 0._R8
      !g%ftEps = 0._R8
      !g%ftBbav2 = 0._R8

      !g%fsPsi = 0._R8

      return
      end subroutine initialize_b2mod_geo

***********************************************************************

      subroutine dealloc_b2mod_geo (g)

      implicit none
      type(geometry_ag), intent(in out)  :: g


      if(b2mod_geo_allocated) then

      !  cell quantites
        deallocate (g%cvBb) 
        deallocate (g%cvX)
        deallocate (g%cvY)
        deallocate (g%cvSz)
        deallocate (g%cvHz)
        deallocate (g%cvHx)
        deallocate (g%cvQgam)
        deallocate (g%cvVol)
        deallocate (g%cvOnedBsq)
        deallocate (g%cvBzb)
        deallocate (g%cvEb) 
      
      !  face quantites
        deallocate (g%fcBb) 
        deallocate (g%fcS)
        deallocate (g%fcHc)
        deallocate (g%fcHt)
        deallocate (g%fcHz)
        deallocate (g%fcVol)      
        deallocate (g%fcQgam)
        deallocate (g%fcQalf)
        deallocate (g%fcQbet)
        deallocate (g%fcPbs)
        deallocate (g%fcPbshz)
        deallocate (g%fcBzb)
        deallocate (g%fcEb) 
      
      !  vertex quantites
        deallocate (g%vxBb)
        deallocate (g%vxX)
        deallocate (g%vxY)
        deallocate (g%vxHz)
        !deallocate (g%vxVol)
        deallocate (g%vxFfbz)
        deallocate (g%vxFpsi)
        deallocate (g%vxOnedBsq)
        deallocate (g%vxBzb)

      !  flux surface quantites
        deallocate (g%cvConn)
        !deallocate (g%ftConn)
        !deallocate (g%ftEps)
        !deallocate (g%ftBbav2)
        !deallocate (g%fsPsi)


        b2mod_geo_allocated = .false.

      endif

      return
      end subroutine dealloc_b2mod_geo

***********************************************************************

      subroutine alloc_b2mod_geo_st(nxd,nyd)
      implicit none
      integer nxd,nyd

      allocate(bb(-1:nxd,-1:nyd,0:3))
      allocate(crx(-1:nxd,-1:nyd,0:3))
      allocate(cry(-1:nxd,-1:nyd,0:3))
      allocate(ffbz(-1:nxd,-1:nyd,0:3))
      allocate(fpsi(-1:nxd,-1:nyd,0:3,0:2))
      allocate(gs(-1:nxd,-1:nyd,0:2))
      allocate(hx(-1:nxd,-1:nyd))
      allocate(hy(-1:nxd,-1:nyd))
      allocate(hz(-1:nxd,-1:nyd))                                    !srv 09.01.01
      allocate(hc(-1:nxd,-1:nyd,0:3))
      allocate(ht(-1:nxd,-1:nyd,0:1))
      allocate(qz(-1:nxd,-1:nyd,0:1))
      allocate(qc(-1:nxd,-1:nyd,0:1))
      allocate(qs(-1:nxd,-1:nyd,0:1))
      allocate(qac(-1:nxd,-1:nyd,0:1))
      allocate(qas(-1:nxd,-1:nyd,0:1))
      allocate(vol(-1:nxd,-1:nyd))
      allocate(pbs(-1:nxd,-1:nyd,0:1))
      allocate(pbshz(-1:nxd,-1:nyd,0:1))                             !srv 09.01.01
      allocate(ebc(-1:nxd,-1:nyd,0:2))
      allocate(conn(-1:nyd))
      allocate(cell_width(-1:nxd,-1:nyd))
      allocate(cell_height(-1:nxd,-1:nyd))

      b2mod_geo_allocated = .true.

      call alloc_b2mod_geo_corner(nxd,nyd)

      return
      end subroutine alloc_b2mod_geo_st

      subroutine dealloc_b2mod_geo_st
      implicit none

      if(b2mod_geo_allocated) then

        deallocate(bb)
        deallocate(crx)
        deallocate(cry)
        deallocate(ffbz)
        deallocate(fpsi)
        deallocate(gs)
        deallocate(hx)
        deallocate(hy)
        deallocate(hz)
        deallocate(hc)
        deallocate(ht)
        deallocate(qz)
        deallocate(qc)
        deallocate(qs)
        deallocate(qac)
        deallocate(qas)
        deallocate(vol)
        deallocate(pbs)
        deallocate(pbshz)                                            !srv 09.01.01
        deallocate(ebc)
        deallocate(conn)
        deallocate(cell_width)
        deallocate(cell_height)

        b2mod_geo_allocated = .false.

      endif

      if (mapInitialized) call deallocateB2GridMap(gmap)

      call dealloc_b2mod_geo_corner

      return
      end subroutine dealloc_b2mod_geo_st

***********************************************************************

      subroutine read_b2mod_geo(nx,ny,nread)
      use b2mod_indirect
      implicit none
      integer nread, nx, ny, n2
      integer ix, iy
      real (kind=R8) ::
     *   p0x, p0y, p1x, p1y
*   ..procedures
      intrinsic min, max, sqrt
      real (kind=R8) ::
     *   trim1, norm, cosang, sinang
      norm(p0x,p0y) = sqrt(p0x**2+p0y**2)
      trim1(p0x) = max(-1.0_R8,min(1.0_R8,p0x))
      sinang(p0x,p0y,p1x,p1y) = trim1(
     &  (p0x*p1x+p0y*p1y)/(norm(p0x,p0y)*norm(p1x,p1y)))
      cosang(p0x,p0y,p1x,p1y) = trim1(
     &  (p0x*p1y-p0y*p1x)/(norm(p0x,p0y)*norm(p1x,p1y)))

      n2 = (nx+1)*(ny+1)
      call cfrure (nread, size(bb), bb, 'bb')
      call cfrure (nread, size(crx), crx, 'crx')
      call cfrure (nread, size(cry), cry, 'cry')
      call cfrure (nread, size(ffbz), ffbz, 'ffbz')
      if (b2fplasma_version.ge.'03.002.000') then
       call cfrure (nread, size(fpsi), fpsi, 'fpsi')
      else
       call cfrure (nread, n2*4, fpsi(-1,-1,0,0), 'fpsi')
       fpsi(:,:,:,1:2) = 0.0_R8
      endif
      call cfrure (nread, size(gs), gs, 'gs')
      call cfrure (nread, size(hx), hx, 'hx')
      call cfrure (nread, size(hy), hy, 'hy')
      call cfrure (nread, size(qz), qz, 'qz')
      if (b2fplasma_version.ge.'03.001.000') then
        call cfrure (nread, size(qc), qc, 'qc')
        call cfrure (nread, size(qs), qs, 'qs')
      elseif (b2fplasma_version.gt.'01.001.021') then
       call cfrure (nread, n2, qc(-1,-1,0), 'qc')
      else
       do iy = -1, ny
        do ix = -1, nx
         qc(ix,iy,0) = cosang (
     &    (crx(ix,iy,1)+crx(ix,iy,3))-(crx(ix,iy,0)+crx(ix,iy,2)),
     &    (cry(ix,iy,1)+cry(ix,iy,3))-(cry(ix,iy,0)+cry(ix,iy,2)),
     &     crx(ix,iy,2)-crx(ix,iy,0),  cry(ix,iy,2)-cry(ix,iy,0))
         if(qc(ix,iy,0).le.0.0_R8) then
          write(*,*) 'Problem in b2mod_geo: qc(ix,iy) <= 0',
     1     qc(ix,iy,0),ix,iy
          write(*,'(a,4f12.6)') 'crx = ',crx(ix,iy,0:3)
          write(*,'(a,4f12.6)') 'cry = ',cry(ix,iy,0:3)
         endif
        enddo
       enddo
      endif
      if(b2fgmtry_version.lt.'03.001.000') then
       qc(:,:,1) = 1.0_R8
       qs(:,:,1) = 0.0_R8
       if(b2fgmtry_version.ge.'01.001.030') then
	 call cfrure_opt (nread, n2, qc(-1,-1,1), 'qcb')
       end if
       if (maxval(qc(:,:,1)).eq.1.0_R8 .and. 
     &     minval(qc(:,:,1)).eq.1.0_R8) then
        do iy = -1, ny
         do ix = -1, nx
          qc(ix,iy,1) = cosang (
     &      crx(ix,iy,1)-crx(ix,iy,0),  cry(ix,iy,1)-cry(ix,iy,0),
     &     (crx(ix,iy,2)+crx(ix,iy,3))-(crx(ix,iy,0)+crx(ix,iy,1)),
     &     (cry(ix,iy,2)+cry(ix,iy,3))-(cry(ix,iy,0)+cry(ix,iy,1)))
          if(qc(ix,iy,1).le.0.0_R8) then
           write(*,'(a,1f12.6,2i5)')
     1      'Problem in b2mod_geo: qc(ix,iy,1) <= 0',qc(ix,iy,1),ix,iy
           write(*,'(a,4f12.6)') 'crx = ', crx(ix,iy,0:3)
           write(*,'(a,4f12.6)') 'cry = ', cry(ix,iy,0:3)
          endif
         enddo
        enddo
       endif
       do iy = -1, ny
        do ix = -1, nx
         qs(ix,iy,0) = sinang (
     &    (crx(ix,iy,1)+crx(ix,iy,3))-(crx(ix,iy,0)+crx(ix,iy,2)),
     &    (cry(ix,iy,1)+cry(ix,iy,3))-(cry(ix,iy,0)+cry(ix,iy,2)),
     &    crx(ix,iy,2)-crx(ix,iy,0), cry(ix,iy,2)-cry(ix,iy,0))
        enddo
       enddo
cwdk   ensure qc(:,:,0)**2 + qs(:,:,0)**2 = 1, in possible event 
cwdk   of improved computation of qc
       qs(:,:,0) = sign(sqrt(1 - qc(:,:,0)**2),qs(:,:,0))
      endif
      call cfrure (nread, size(vol), vol, 'vol')
c conn
      if(b2fplasma_version.ge.'01.001.063') then        !srv 15.02.11
        call cfrure (nread, size(pbs), pbs, 'pbs')
      endif
      if(b2fplasma_version.ge.'01.001.033') then
        call cfrure_opt
     .   (nread, size(cell_width), cell_width, 'cell_width')
        call cfrure_opt
     .   (nread, size(cell_height), cell_height, 'cell_height')
      else
        if(isClassicalGrid(cflags)) then
          cell_width = hx
          cell_height(:,:) = hy(:,:)*qz(:,:,1)
        else
          call xerrab ('Cell widths and heights not set: re-run b2ag!')
        endif
      endif

      return
      end subroutine read_b2mod_geo

***********************************************************************

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
      if(newversion.ge.'03.001.000') then
        call cfwure (nwrite, size(qs), qs, 'qs')
      endif

      call cfwure (nwrite, size(vol), vol, 'vol')
c conn
      if(newversion.ge.'01.001.063') then               !srv 15.02.11
        call cfwure (nwrite, size(pbs), pbs, 'pbs')
      endif
      return
      end subroutine write_b2mod_geo

      end module b2mod_geo
