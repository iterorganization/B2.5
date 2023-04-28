      module b2mod_layer
      use b2mod_types
      implicit none

      integer nwall, ntrack, nalloys
      parameter(ntrack=4,nalloys=6)

#ifdef ALLOCATE_AND_NAMELIST
      integer, allocatable :: layer_nconstituents(:)
      integer, allocatable :: layer_nzconstituents(:,:)
      real (kind=R8), allocatable :: layer_nrelconstituents(:,:)
#else
#     include "DIMENSIONS.F"
      parameter(
     . nwall = 2*(DEF_NXD+2)+(DEF_NCUT+2)*(DEF_NYD+2)+DEF_NLIM+DEF_NSTS)
      integer layer_nconstituents(nwall)
      integer layer_nzconstituents(nwall,6+ntrack)
      real (kind=R8) :: layer_nrelconstituents(nwall,6+ntrack)
#endif
      real (kind=R8), allocatable :: layer_thickness(:)
      real (kind=R8), allocatable :: layer_mass_density(:)
      real (kind=R8), allocatable :: layer_depth(:)
      real (kind=R8), allocatable :: layer_composition(:,:)
      real (kind=R8), allocatable :: particle_inventory(:,:)
      real (kind=R8), save :: layer_thermal_contact

      character*6 layer_alloys(nalloys)
      integer alloy_composition(nalloys,ntrack)
      character*32, save :: alloy_long_name(nalloys)
      character*1, save :: alloy_temp_unit(nalloys)
      integer, save :: alloy_nconstituents(nalloys)
      integer, save ::
     *  alloy_format(nalloys),
     *  alloy_nzconstituents(nalloys,ntrack),
     *  alloy_nrelconstituents(nalloys,ntrack)
      real (kind=R8), save ::
     *  alloy_chemsput_erel(nalloys),
     *  alloy_blackbody_emiss(nalloys),
     *  alloy_vapprs_A(nalloys),
     *  alloy_vapprs_B(nalloys),
     *  alloy_stick_coef(nalloys),
     *  alloy_vapor_temp(nalloys),
     *  alloy_heat_sublim(nalloys)
      real (kind=R8), save ::
     *  alloy_target_density(nalloys),
     *  alloy_Kappa_KA(nalloys),
     *  alloy_Kappa_KB(nalloys),
     *  alloy_Kappa_A1(nalloys),
     *  alloy_Kappa_A2(nalloys),
     *  alloy_Kappa_A3(nalloys),
     *  alloy_Kappa_A4(nalloys),
     *  alloy_Kappa_A5(nalloys),
     *  alloy_Cp_B1(nalloys),
     *  alloy_Cp_B2(nalloys),
     *  alloy_Cp_B3(nalloys),
     *  alloy_Cp_B4(nalloys)

      character*32, save :: species_long_name(ntrack)
      character*1, save :: species_temp_unit(ntrack)
      integer, save :: species_cp_format(ntrack)
      real (kind=R8), save ::
     *  species_target_density(ntrack),
     *  species_Kappa_KA(ntrack),
     *  species_Kappa_KB(ntrack),
     *  species_Kappa_A1(ntrack),
     *  species_Kappa_A2(ntrack),
     *  species_Kappa_A3(ntrack),
     *  species_Kappa_A4(ntrack),
     *  species_Kappa_A5(ntrack),
     *  species_Cp_B1(ntrack),
     *  species_Cp_B2(ntrack),
     *  species_Cp_B3(ntrack),
     *  species_Cp_B4(ntrack),
     *  species_blackbody_emiss(ntrack),
     *  species_vapprs_A(ntrack),
     *  species_vapprs_B(ntrack),
     *  species_stick_coef(ntrack),
     *  species_vapor_temp(ntrack),
     *  species_heat_sublim(ntrack)

      logical, save :: layer_model
      data layer_model /.false./

      contains

      subroutine alloc_b2mod_layer(nx,ny,nncut,nlim,nsts)
      implicit none
      integer nx,ny,nncut,nlim,nsts
      external xertst

*   ..test arguments
      call xertst (0.le.nx.and.0.le.ny, 'faulty argument nx, ny')
      call xertst (0.le.nncut, 'faulty argument nncut')
      call xertst (0.le.nlim, 'faulty argument nlim')
      call xertst (0.le.nsts, 'faulty argument nsts')

#ifdef ALLOCATE_AND_NAMELIST
      nwall = 2*(nx+2)+max(nncut,2)*(ny+2)+nlim+nsts
      allocate(layer_nconstituents(nwall))
      allocate(layer_nzconstituents(nwall,6+ntrack))
      allocate(layer_nrelconstituents(nwall,6+ntrack))
#endif
      allocate(layer_thickness(nwall))
      allocate(layer_mass_density(nwall))
      allocate(layer_depth(nwall))
      allocate(layer_composition(nwall,6+ntrack+nalloys))
      allocate(particle_inventory(nwall,ntrack))

      layer_alloys = ' '
      layer_nconstituents = 0
      layer_nzconstituents = 0
      layer_nrelconstituents = 0.0_R8

      return

      end subroutine alloc_b2mod_layer

      subroutine dealloc_b2mod_layer
      if (.not.allocated(layer_thickness)) return
#ifdef ALLOCATE_AND_NAMELIST
      deallocate(layer_nconstituents)
      deallocate(layer_nzconstituents)
      deallocate(layer_nrelconstituents)
#endif
      deallocate(layer_thickness)
      deallocate(layer_mass_density)
      deallocate(layer_depth)
      deallocate(layer_composition)
      deallocate(particle_inventory)
      end subroutine dealloc_b2mod_layer

      end module b2mod_layer
