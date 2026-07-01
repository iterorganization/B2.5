#ifdef CONSTANTS_PROVIDED
include 'mathematical_constants.f90'   ! IGNORE
#if CONSTANTS_MINOR_VERSION > 1
include 'codata_2022.f90'              ! IGNORE
#else
include 'codata_2018.f90'              ! IGNORE
#endif
#endif

  module b2mod_constants
  use b2mod_types
!     -------------------------------------------------------------------------
!     09/12/2007 xpb : source CODATA 2006 (http://www.nist.gov/)
!     08/19/2011 xpb : source CODATA 2010 (http://www.nist.gov/)
!     03/19/2018 xpb : source CODATA 2014 (http://www.nist.gov/)
!     03/06/2019 dpc : source CODATA 2018 (http://www.nist.gov/)
!     26/03/2026 xpb : source CODATA 2022 (http://www.nist.gov/)
!     -------------------------------------------------------------------------
#ifdef CONSTANTS_PROVIDED
  use mathematical_constants  ! IGNORE
  use codata                  ! IGNORE
#endif
  real (kind=R8) :: pi
  real (kind=R8) :: c, me, mp, ev, qe, eps0, mu0, avogr, KBolt
#ifdef CONSTANTS_PROVIDED
  parameter (pi=m_pi)
  parameter (c=speed_of_light_in_vacuum)
  parameter (me=electron_mass)
  parameter (mp=proton_mass)
  parameter (ev=electron_volt, qe=ev)
  parameter (mu0=vacuum_mag_permeability)
  parameter (eps0=vacuum_electric_permittivity)
  parameter (avogr=Avogadro_constant)
  parameter (KBolt=Boltzmann_constant)
  real (kind=R8) :: const_h
  parameter (const_h=Planck_constant)
#else
  parameter (pi=3.141592653589793238462643383280_R8)
  parameter (c=2.99792458e8_R8, me=9.1093837139e-31_R8,    &
 &  mp=1.67262192595e-27_R8, ev=1.602176634e-19_R8, qe=ev, &
 &  mu0=1.25663706127e-6_R8, eps0=8.8541878188e-12_R8,     &
 &  avogr=6.02214076e23_R8, KBolt=1.380649e-23_R8)
#endif
  end module b2mod_constants

!!!Local Variables:
!!! mode: f90
!!! End:
