      module b2mod_wall2
      use b2mod_types
      use b2mod_wall
     . , only : nwall, ndepth, ntrack
#     include "DIMENSIONS.F"
      
cxpb  Namelist data
      integer, save :: imapx1(nwall), imapy1(nwall)
      logical inertial_cooling1(nwall)
      real (kind=R8), save ::
     *  target_temp1(nwall,ndepth),
     *  backplate_temp1(nwall),
     *  plate_thickness1(nwall),
     *  plate_time_factor1(nwall),
     *  coating_thickness1(nwall),
     *  deposition1(nwall,ntrack),
     *  monolayer_deposition1(nwall,ntrack),
     *  erosion1(nwall,ntrack),
     *  monolayer_erosion1(nwall,ntrack),
     *  plate_time1(nwall),
     *  plate_area1(nwall),
     *  chemical_sputtering1
     .   (nwall,0:DEF_NSD-1+DEF_NATM+DEF_NMOL+DEF_NION),
     *  physical_sputtering1
     .   (nwall,0:DEF_NSD-1+DEF_NATM+DEF_NMOL+DEF_NION),
     *  physical_sputtering_energy1
     .   (nwall,0:DEF_NSD-1+DEF_NATM+DEF_NMOL+DEF_NION),
     *  RES_sputtering1(nwall,0:DEF_NSD-1+DEF_NATM+DEF_NMOL+DEF_NION),
     *  thermal_evaporation1
     .   (nwall,0:DEF_NSD-1+DEF_NATM+DEF_NMOL+DEF_NION),
     *  backscattering1(nwall,0:DEF_NSD-1+DEF_NATM+DEF_NMOL+DEF_NION),
     *  backscattering_energy1
     .   (nwall,0:DEF_NSD-1+DEF_NATM+DEF_NMOL+DEF_NION)

cxpb Data file material properties
      character*6, save ::
     *  surface_material_name1(nwall),
     *  coating_material_name1(nwall),
     *  bulk_material_name1(nwall)
      integer, save ::
     *  layer_nconstituents1(nwall),
     *  layer_nzconstituents1(nwall,6+ntrack)
      real (kind=R8), save ::
     *  layer_nrelconstituents1(nwall,6+ntrack)
      
#ifdef ALLOCATE_AND_NAMELIST
      integer, allocatable, save ::
     *  xymap1(:,:)
#else
      integer xymap1(-1:DEF_NXD,-1:DEF_NYD)
#endif

cxpb Mapping arrays:
cxpb A wall element numbered iwall can be found in the new B2 grid at the position
cxpb   ( imapx1(iwall) , imapy1(iwall) )
cxpb Reversely, the B2 boundary cell ( ix, iy ) corresponds to the boundary
cxpb element: xymap1(ix,iy)

      contains

      subroutine alloc_b2mod_wall2(nx1,ny1,ndepth1)
      implicit none
      integer nx1, ny1, ndepth1
      logical errdim
      external subini, subend, xertst

      call subini('b2mod_wall2')
#ifdef ALLOCATE_AND_NAMELIST
      allocate(xymap1(-1:nx1,-1:ny1))
#else
      errdim=.false.
      if(nx1.gt.DEF_NXD) then
        write(*,*) 'alloc_b2mod_wall2: nxd > def_nxd',nx1,DEF_NXD
        errdim=.true.
      endif
      if(ny1.gt.DEF_NYD) then
        write(*,*) 'alloc_b2mod_wall2: nyd > def_nyd',ny1,DEF_NYD
        errdim=.true.
      endif
      call xertst (.not.errdim,
     1 'faulty def_nxd, def_nyd')
#endif
      call xertst(1.le.ndepth1 .and. ndepth1.le.ndepth,
     . 'invalid dimension ndepth1')

      imapx1 = 0
      imapy1 = 0
      xymap1 = 0
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
      plate_area1 = 0.0_R8
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

      end module b2mod_wall2
