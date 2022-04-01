      module b2mod_boundary_sources                          !srv 26.02.18 {
      use b2mod_types
      implicit none
      real (kind=R8), allocatable, save ::
     *  bv_po(:,:),
     *  bv_te(:,:),
     *  bv_ti(:,:),
     *  bv_ua(:,:,:),
     *  bv_na(:,:,:)

      contains

      subroutine alloc_b2mod_boundary_sources(nxd,nyd,nsd)
      implicit none
      integer nxd, nyd, nsd

      allocate(bv_po(-1:nxd,-1:nyd))
      allocate(bv_te(-1:nxd,-1:nyd))
      allocate(bv_ti(-1:nxd,-1:nyd))
      allocate(bv_ua(-1:nxd,-1:nyd,0:nsd-1))
      allocate(bv_na(-1:nxd,-1:nyd,0:nsd-1))

      bv_po = 0.0_R8
      bv_te = 0.0_R8
      bv_ti = 0.0_R8
      bv_ua = 0.0_R8
      bv_na = 0.0_R8

      write(*,*) 'Allocated space for b2mod_boundary_sources'
      return

      end subroutine alloc_b2mod_boundary_sources

      subroutine dealloc_b2mod_boundary_sources
      implicit none

      deallocate(bv_po)
      deallocate(bv_te)
      deallocate(bv_ti)
      deallocate(bv_ua)
      deallocate(bv_na)

      return
      end subroutine dealloc_b2mod_boundary_sources

      subroutine read_boundary_sources(nx, ny, ns)
      use b2mod_types
      use b2mod_boundary_namelist
      implicit none
      integer nx, ny, ns

      integer is, ix, iy
      character*36 filename
      logical lexist

      filename='bv_na.dat'
      call find_file(filename,lexist)
      call xertst(lexist,'File '//trim(filename)//' is not found')
      open(70,file=trim(filename))
      read (70,*) (((bv_na(ix,iy,is),ix=-1,nx),iy=-1,ny),is=0,ns-1)
      close(70)

      filename='bv_ua.dat'
      call find_file(filename,lexist)
      call xertst(lexist,'File '//trim(filename)//' is not found')
      open(70,file=trim(filename))
      read (70,*) (((bv_ua(ix,iy,is),ix=-1,nx),iy=-1,ny),is=0,ns-1)
      close(70)

      filename='bv_te.dat'
      call find_file(filename,lexist)
      call xertst(lexist,'File '//trim(filename)//' is not found')
      open(70,file=trim(filename))
      read (70,*) ((bv_te(ix,iy),ix=-1,nx),iy=-1,ny)
      close(70)

      filename='bv_ti.dat'
      call find_file(filename,lexist)
      call xertst(lexist,'File '//trim(filename)//' is not found')
      open(70,file=trim(filename))
      read (70,*) ((bv_ti(ix,iy),ix=-1,nx),iy=-1,ny)
      close(70)

      filename='bv_po.dat'
      call find_file(filename,lexist)
      call xertst(lexist,'File '//trim(filename)//' is not found')
      open(70,file=trim(filename))
      read (70,*) ((bv_po(ix,iy),ix=-1,nx),iy=-1,ny)
      close(70)

      return
      end subroutine read_boundary_sources

      end module b2mod_boundary_sources                      !srv 26.02.18 }
