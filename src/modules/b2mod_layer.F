      module b2mod_layer
      use b2mod_types
#     include "DIMENSIONS.F"
      implicit none

      integer nwall, ntrack, nalloys
      parameter(
     *  nwall=2*DEF_NCUT*(DEF_NXD+2)+2*(DEF_NYD+2)+DEF_NLIM+DEF_NSTS,
     *  ntrack=4,
     *  nalloys=6)

      integer layer_nconstituents(nwall)
      integer layer_nzconstituents(nwall,6+ntrack)
      real (kind=R8) :: layer_nrelconstituents(nwall,6+ntrack)
      real (kind=R8) :: layer_thickness(nwall)
      real (kind=R8) :: layer_mass_density(nwall)
      real (kind=R8) :: layer_depth(nwall)
      real (kind=R8) :: layer_composition(nwall,6+ntrack+nalloys)
      real (kind=R8) :: particle_inventory(nwall,ntrack)
      real (kind=R8) :: layer_thermal_contact

      character*6 layer_alloys(nalloys)
      integer alloy_composition(nalloys,ntrack)
      character*32, save :: alloy_long_name(nalloys)
      integer, save :: alloy_nconstituents(nalloys)
      integer, save ::
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
     *  alloy_Cp_B1(nalloys),
     *  alloy_Cp_B2(nalloys),
     *  alloy_Cp_B3(nalloys)

      character*32, save :: species_long_name(ntrack)
      real (kind=R8), save ::
     *  species_target_density(ntrack),
     *  species_Kappa_KA(ntrack),
     *  species_Kappa_KB(ntrack),
     *  species_Kappa_A1(ntrack),
     *  species_Kappa_A2(ntrack),
     *  species_Kappa_A3(ntrack),
     *  species_Kappa_A4(ntrack),
     *  species_Cp_B1(ntrack),
     *  species_Cp_B2(ntrack),
     *  species_Cp_B3(ntrack),
     *  species_blackbody_emiss(ntrack),
     *  species_vapprs_A(ntrack),
     *  species_vapprs_B(ntrack),
     *  species_stick_coef(ntrack),
     *  species_vapor_temp(ntrack),
     *  species_heat_sublim(ntrack)

      end module b2mod_layer
