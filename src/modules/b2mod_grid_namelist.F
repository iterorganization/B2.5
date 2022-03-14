      module b2mod_grid_namelist
      use b2mod_types
      use b2mod_indirect
      implicit none

#ifdef ALLOCATE_AND_NAMELIST
      real (kind=R8), allocatable, save ::
     *     R1(;), 
     *     Z1(;), 
     *     R2(;), 
     *     Z2(:),      
     *     CURRENT(:)
      integer, allocatable, save ::
     *     WINDINGS(:),
     *     startCellR(:),
     *     startCellZ(:),
     *     endCellR(:),
     *     endCellZ(:)
      real (kind=R8), allocatable, save :: 
     *     regionR(:, :),
     *     regionZ(:, :)
      integer, allocatable, save :: nNodesInExternalRegion(:)
      character*260, allocatable, save :: area(:), angle(:)
#else
#     include "DIMENSIONS.F"
      integer NCOIL_MAX
      parameter (NCOIL_MAX=20)
      real (kind=R8), save :: R1(0:NCOIL_MAX - 1)
      real (kind=R8), save :: R2(0:NCOIL_MAX - 1)
      real (kind=R8), save :: Z1(0:NCOIL_MAX - 1)
      real (kind=R8), save :: Z2(0:NCOIL_MAX - 1)
      real (kind=R8), save :: CURRENT(0:NCOIL_MAX - 1)
      integer, save :: WINDINGS(0:NCOIL_MAX - 1)
	  integer, save :: nlayers(0:NCOIL_MAX - 1)										

      integer NEXTERNAL_REGIONS_MAX
      parameter (NEXTERNAL_REGIONS_MAX = 10)
      integer NNODES_IN_EXTERNAL_REGION_MAX
      parameter (NNODES_IN_EXTERNAL_REGION_MAX = 50)
      integer, save :: startCellR(0:NEXTERNAL_REGIONS_MAX - 1)
      integer, save :: endCellR(0:NEXTERNAL_REGIONS_MAX - 1)
      integer, save :: startCellZ(0:NEXTERNAL_REGIONS_MAX - 1)
      integer, save :: endCellZ(0:NEXTERNAL_REGIONS_MAX - 1)

      real (kind=R8), save :: 
     *     regionR(0:NNODES_IN_EXTERNAL_REGION_MAX - 1, 
     *     0:NEXTERNAL_REGIONS_MAX - 1)
      real (kind=R8), save :: 
     *     regionZ(0:NNODES_IN_EXTERNAL_REGION_MAX - 1,
     *     0:NEXTERNAL_REGIONS_MAX - 1)
      integer, save :: nNodesInRegion(0:NEXTERNAL_REGIONS_MAX - 1)
      character*260, save :: 
     *     area(0:NEXTERNAL_REGIONS_MAX - 1),
     *     angle(0:NEXTERNAL_REGIONS_MAX - 1)
#endif

      integer, save :: NX, NY, NCOIL, nExternalRegions,
     *     triangulationOption
      real (kind=R8), save :: r_ne, r_nw, r_se, r_sw, 
     *     z_ne, z_nw, z_se, z_sw,
     *     r_start, r_end, z_start, z_end,
     .     bTor, stepsize, gc_n, gc_s, gc_e, gc_w
      character*260, save :: fn_z, fn_r
      logical, save :: triangularGrid, b25Grid

      contains


      subroutine alloc_b2mod_grid_namelist(nx_in,ny_in,ncoil_in,
     *     nExternalRegions_in, nNodesInExternalRegion_in)
      implicit none

      integer nx_in, ny_in, ncoil_in, nExternalRegions_in, 
     *     nNodesInExternalRegion_in

#ifndef ALLOCATE_AND_NAMELIST
      external xertst
#endif

#ifdef ALLOCATE_AND_NAMELIST
      nx=nx_in
      ny=ny_in
      ncoil=NCOIL_IN
      allocate(R1(0:ncoil - 1))
      allocate(R2(0:ncoil - 1))
      allocate(Z1(0:ncoil - 1))
      allocate(Z2(0:ncoil - 1))
      allocate(CURRENT(0:ncoil - 1))
      allocate(WINDINGS(0:ncoil - 1))
      allocate(regionR(0:nNodesInExternalRegion_in - 1,
     *     0:nExternalRegions_in - 1))
      allocate(regionZ(0:nNodesInExternalRegion_in - 1,
     *     0:nExternalRegions_in - 1))
      allocate(nNodesInExternalRegion(0:nExternalRegions_in - 1))

      allocate(startCellR(0:ExternalRegions_in - 1))
      allocate(endCellR(0:ExternalRegions_in - 1))
      allocate(startCellZ(0:ExternalRegions_in - 1))
      allocate(endCellZ(0:ExternalRegions_in - 1))

#endif
      bTor = 0.0_R8
      nlayers(:) = 5 ! default value is 5
      R1(:) = 0.0_R8 
      Z1(:) = 0.0_R8 
      R2(:) = 0.0_R8 
      Z2(:) = 0.0_R8      
      CURRENT(:) = 0.0_R8 
      r_ne = 0.0_R8 
      r_nw = 0.0_R8 
      r_se = 0.0_R8 
      r_sw = 0.0_R8 
      z_ne = 0.0_R8 
      z_nw = 0.0_R8 
      z_se = 0.0_R8 
      z_sw = 0.0_R8 
      r_start = 0.0_R8 
      r_end = 0.0_R8 
      z_start = 0.0_R8 
      z_end = 0.0_R8 
      stepsize = 1.0e-4_R8
      windings(:) = 0
      fn_z = 'z.dat'
      fn_r = 'r.dat'
      gc_n = 1.0e-4_R8
      gc_s = 1.0e-4_R8
      gc_w = 1.0e-4_R8
      gc_e = 1.0e-4_R8

      angle(:) = ''
      area(:) = ''
      triangularGrid = .FALSE.
      b25Grid = .FALSE.
      triangulationOption = 0
      regionR(:, :) = 0.0_R8
      regionZ(:, :) = 0.0_R8
      nNodesInRegion(:) = 0
      startCellR(:) = -1
      endCellR(:) = -1
      startCellZ(:) = -1
      endCellZ(:) = -1

      end subroutine alloc_b2mod_grid_namelist


      subroutine read_b2mod_grid_namelist(grid_filename, dotest)
      implicit none

      integer is
      real (kind=R8) :: eps, tmp
      external xertst, xerrab

      namelist /grid/ R1, R2, Z1, Z2, current, windings,
     1     ncoil, nx, ny, nlayers, fn_r, fn_z,
     2     r_ne, r_nw, r_se, r_sw, z_ne, z_nw, z_se, z_sw,
     3     r_start, r_end, z_start, z_end,      
     4     gc_n, gc_s, gc_w, gc_e, stepsize, bTor,
     5     b25Grid, triangularGrid, triangulationOption,
     6     nExternalRegions, regionR, regionZ, nNodesInRegion,
     7     startCellR, endCellR, startCellZ, endCellZ,
     8     angle, area

      logical :: file_ok, dotest
      character*260 filename, grid_filename

      eps=1.0e-7_R8
      filename=grid_filename
      inquire(file=filename,exist=file_ok)
      if(file_ok) then
        open(99,file=filename)
        write(*,*) 'Reading namelist grid from ',trim(filename)
        read(99,grid)
        write(*,*) 'Read namelist grid ',trim(filename)
      else 
        write (*,*) 'read_b2mod_grid_namelist: file =', filename
        call xerrab('grid input file did not exist')
      end if

      if ((.not.file_ok).or.(.not.dotest)) then
        if (.not.file_ok) write(*,*) 'No ',trim(filename)
        write(*,*) 'Testing skipped!'
        return
      end if

*   ..test input
      if (b25Grid) then
        ! create a B2.5 grid  

        call xertst (0.lt.nx, 'faulty input nx')
        call xertst (0.lt.ny, 'faulty input ny')
        call xertst (0.lt.ncoil, 'faulty input ncoil')
        do is = 0, ncoil-1
          call xertst (0.lt.nlayers(is), 'faulty input nlayers')														  
          call xertst (eps.le.current(is), 'faulty input current')
          call xertst (0.le.windings(is), 'faulty input windings')
          call xertst (0.0_R8.le.R1(is), 'faulty input R_in')
          call xertst (0.0_R8.le.R2(is), 'faulty input R_out')
          if(R1(is).gt.R2(is)) then
            tmp=R1(is)
            R1(is) = R2(is)
            R2(is) = R1(is)
          end if
          if(Z1(is).gt.Z2(is)) then
            tmp=Z1(is)
            Z1(is) = Z2(is)
            Z2(is) = Z1(is)
          end if
          call xertst ((R2(is)-R1(is)).gt.eps, 
     .         'difference R1, R2 too small')
          call xertst ((Z2(is)-Z1(is)).gt.eps, 
     ,         'difference Z1, Z2 too small')
          call xertst (r_ne.ge.0.0_R8, 'r_ne should be positive')
          call xertst (r_nw.ge.0.0_R8, 'r_nw should be positive')
          call xertst (r_se.ge.0.0_R8, 'r_se should be positive')
          call xertst (r_sw.ge.0.0_R8, 'r_sw should be positive')
          call xertst (stepsize.gt.0.0_R8, 
     ,         'stepsize negative, should be a small positive number')
          call xertst (stepsize.lt.0.1_R8, 
     ,         'stepsize too large, should be a small positive number')
        end do

      else
        ! read in existing grid
        ! check for OUT.ASCII to be existing
        inquire(file='OUT.ASCII',exist=file_ok)
        call xertst (file_ok, 
     ,         'OUT.ASCII file required, but not available')
      end if

      if (triangularGrid) then
        select case (triangulationOption)
        case (0, 1)
          ! 0: Split B2.5 cell in 2 triangles
          ! 1: split B2.5 cell into 4 triangles (centroid as additional node)
          ! nothing to be done

        case default
          call xerrab('Undefined triangulation option')
        end select
        call xertst(nExternalRegions.ge.0,
     *       'nExternalRegions should be 0 or positive.')
        do is = 0, nExternalRegions - 1
          call xertst(nNodesInRegion(is).ge.1, ! changed by JAS
     *         'nNodesInRegion should be minimally 1')
          call xertst(nNodesInRegion(is).le.
     *         NNODES_IN_EXTERNAL_REGION_MAX,
     *         'nNodesInRegion exceeds the maximum value')
          call xertst(startCellR(is) == endCellR(is) .OR.
     *         startCellZ(is) == endCellZ(is),
     *         'External region must connect to single B2.5 boundary')
          call xertst(startCellR(is) /= endCellR(is) .OR.
     *         startCellZ(is) /= endCellZ(is),
     *         'External region must connect to B2.5 mesh by >1 node')
        end do
      end if

      return

      end subroutine read_b2mod_grid_namelist


      subroutine write_b2mod_grid_namelist
      implicit none

      namelist /grid/ R1, R2, Z1, Z2, current, windings,
     1     ncoil, nx, ny, nlayers, fn_r, fn_z,
     2     r_ne, r_nw, r_se, r_sw, z_ne, z_nw, z_se, z_sw,  
     3     r_start, r_end, z_start, z_end,      
     4     gc_n, gc_s, gc_w, gc_e, stepsize, bTor,
     5     b25Grid, triangularGrid, triangulationOption,
     6     nExternalRegions, regionR, regionZ, nNodesInRegion,
     7     startCellR, endCellR, startCellZ, endCellZ,
     8     angle, area

      write(*,grid)

      end subroutine write_b2mod_grid_namelist

      end module b2mod_grid_namelist
