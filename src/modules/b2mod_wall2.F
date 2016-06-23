      module b2mod_wall2
      use b2mod_types
      use b2mod_wall
     . , only : ndepth, ntrack, nwall_segments
      use b2mod_layer
     . , only : nwall
#     include "DIMENSIONS.F"
      
      integer, save :: nwall1, nbnd1
      real (kind=R8), save, allocatable :: contr_area(:,:)
      real (kind=R8), save, allocatable :: acc_area(:), acc_area1(:)
cxpb  Namelist data
      integer, allocatable, save :: imapx1(:), imapy1(:)
      logical, allocatable :: inertial_cooling1(:)
      real (kind=R8), allocatable, save ::
     *  target_temp1(:,:),
     *  backplate_temp1(:),
     *  plate_thickness1(:),
     *  plate_time_factor1(:),
     *  coating_thickness1(:),
     *  deposition1(:,:),
     *  monolayer_deposition1(:,:),
     *  erosion1(:,:),
     *  monolayer_erosion1(:,:),
     *  plate_time1(:),
     *  plate_area1(:),
     *  chemical_sputtering1(:,:),
     *  physical_sputtering1(:,:),
     *  physical_sputtering_energy1(:,:),
     *  RES_sputtering1(:,:),
     *  thermal_evaporation1(:,:),
     *  backscattering1(:,:),
     *  backscattering_energy1(:,:)

cxpb Data file material properties
      character*6, allocatable, save ::
     *  surface_material_name1(:),
     *  coating_material_name1(:),
     *  bulk_material_name1(:)
      integer, allocatable, save ::
     *  layer_nconstituents1(:),
     *  layer_nzconstituents1(:,:)
      real (kind=R8), allocatable, save ::
     *  layer_nrelconstituents1(:,:)
      
      integer, allocatable, save ::
     *  xymap1(:,:), wall_begin1(:), wall_end1(:)

cxpb Mapping arrays:
cxpb A wall element numbered iwall can be found in the new B2 grid at the position
cxpb   ( imapx1(iwall) , imapy1(iwall) )
cxpb Reversely, the B2 boundary cell ( ix, iy ) corresponds to the boundary
cxpb element: xymap1(ix,iy)

      contains

      subroutine alloc_b2mod_wall2_init(nncut1,nx1,ny1,nlim,nsts)
      implicit none
      integer nncut1, nx1, ny1, nlim, nsts

      allocate(xymap1(-1:nx1,-1:ny1))
      nwall1 = 2*(nx1+2)+max(nncut1,2)*(ny1+2)+nlim+nsts

      allocate(imapx1(nwall1), imapy1(nwall1))
      allocate(acc_area(nwall),acc_area1(nwall1))
      allocate(contr_area(nwall,nwall1))
      allocate(plate_area1(nwall1))
      allocate(wall_begin1(nwall_segments))
      allocate(wall_end1(nwall_segments))

      xymap1 = 0
      imapx1 = 0
      imapy1 = 0
      wall_begin1 = 0
      wall_end1 = 0
      contr_area = 0.0_R8
      plate_area1 = 0.0_R8
      acc_area = 0.0_R8
      acc_area1 = 0.0_R8

      return
      end subroutine alloc_b2mod_wall2_init

      subroutine alloc_b2mod_wall2(ns1,natm1,nmol,nion,ndepth1)
      implicit none
      integer ns1, natm1, nmol, nion, ndepth1
      external subini, subend, xertst

      call subini('b2mod_wall2')
      call xertst(1.le.ndepth1 .and. ndepth1.le.ndepth,
     . 'invalid dimension ndepth1')

cxpb  Namelist data
      allocate(inertial_cooling1(nwall1))
      allocate(
     *  target_temp1(nwall1,ndepth),
     *  backplate_temp1(nwall1),
     *  plate_thickness1(nwall1),
     *  plate_time_factor1(nwall1),
     *  coating_thickness1(nwall1),
     *  deposition1(nwall1,ntrack),
     *  monolayer_deposition1(nwall1,ntrack),
     *  erosion1(nwall1,ntrack),
     *  monolayer_erosion1(nwall1,ntrack),
     *  plate_time1(nwall1),
     *  chemical_sputtering1(nwall1,0:ns1-1+natm1+nmol+nion),
     *  physical_sputtering1(nwall1,0:ns1-1+natm1+nmol+nion),
     *  physical_sputtering_energy1(nwall1,0:ns1-1+natm1+nmol+nion),
     *  RES_sputtering1(nwall1,0:ns1-1+natm1+nmol+nion),
     *  thermal_evaporation1(nwall1,0:ns1-1+natm1+nmol+nion),
     *  backscattering1(nwall1,0:ns1-1+natm1+nmol+nion),
     *  backscattering_energy1(nwall1,0:ns1-1+natm1+nmol+nion))

cxpb Data file material properties
      allocate(
     *  surface_material_name1(nwall1),
     *  coating_material_name1(nwall1),
     *  bulk_material_name1(nwall1))
      allocate(
     *  layer_nconstituents1(nwall1),
     *  layer_nzconstituents1(nwall1,6+ntrack))
      allocate(
     *  layer_nrelconstituents1(nwall1,6+ntrack))
      
      inertial_cooling1 = .false.
      target_temp1 = 0.0_R8
      backplate_temp1 = 0.0_R8
      plate_thickness1 = 0.0_R8
      coating_thickness1 = 0.0_R8
      plate_time_factor1 = 0.0_R8
      deposition1 = 0.0_R8
      monolayer_deposition1 = 0.0_R8
      erosion1 = 0.0_R8
      monolayer_erosion1 = 0.0_R8
      plate_time1 = 0.0_R8
      chemical_sputtering1 = 0.0_R8
      physical_sputtering1 = 0.0_R8
      physical_sputtering_energy1 = 0.0_R8
      RES_sputtering1 = 0.0_R8
      thermal_evaporation1 = 0.0_R8
      backscattering1 = 0.0_R8
      backscattering_energy1 = 0.0_R8

      call subend()
      return

      end subroutine alloc_b2mod_wall2

      subroutine dealloc_b2mod_wall2_init
      implicit none

      if (.not.allocated(xymap1)) return
      deallocate(xymap1)
      deallocate(imapx1, imapy1)
      deallocate(acc_area, acc_area1)
      deallocate(contr_area)
      deallocate(plate_area1)
      deallocate(wall_begin1)
      deallocate(wall_end1)
      return
      end subroutine dealloc_b2mod_wall2_init

      subroutine dealloc_b2mod_wall2
      implicit none

      if (.not.allocated(inertial_cooling1)) return
      deallocate(inertial_cooling1)
      deallocate(
     *  target_temp1,
     *  backplate_temp1,
     *  plate_thickness1,
     *  plate_time_factor1,
     *  coating_thickness1,
     *  deposition1,
     *  monolayer_deposition1,
     *  erosion1,
     *  monolayer_erosion1,
     *  plate_time1,
     *  chemical_sputtering1,
     *  physical_sputtering1,
     *  physical_sputtering_energy1,
     *  RES_sputtering1,
     *  thermal_evaporation1,
     *  backscattering1,
     *  backscattering_energy1)
      deallocate(
     *  surface_material_name1,
     *  coating_material_name1,
     *  bulk_material_name1)
      deallocate(
     *  layer_nconstituents1,
     *  layer_nzconstituents1)
      deallocate(
     *  layer_nrelconstituents1)
      
      return
      end subroutine dealloc_b2mod_wall2

      end module b2mod_wall2
