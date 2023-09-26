c**************************************************************************************
c********************  NEW KU Leuven module for ECH/ECWC plasmas  *********************
c**************************************************************************************
! mb  01/08/23
! main author: Anthony Piras
c   This module is intended to simulate power depositions in Electron Cyclotron Heated
c     plasmas. To turn on the model, the flag b2stbr_ech has to be set to '1'. This
c     module computes an additional heating source due to EC heating injected either
c     at the cyclotron frequency or its harmonics
      module b2mod_ech
*   ..other modules called 
      use b2mod_types
      use b2mod_constants
      use b2mod_wall
      use b2mod_math
      use b2mod_types
      use b2mod_math
      use b2mod_transport_namelist
      use b2mod_input_profile
      use b2mod_constants
      use b2mod_time
      use b2mod_diag
      use b2mod_b2cmpa
      use b2mod_b2cmpt
      use b2mod_b2cmpb
      use b2mod_ranges
      use b2mod_neoclassical
      use b2mod_geo_corner
      use b2mod_user_namelist
      use b2mod_switches
      use b2us_geo
      use b2us_map
      use b2us_plasma

      implicit none
      save

      real (kind=R8), allocatable, save ::
     &  sECH(:) 

      contains

c Alloc
      subroutine alloc_b2mod_ech(nFcd)
      implicit none
      integer,intent(in) :: nFcd

      allocate(sECH(nFcd))

      sECH = 0.0_R8

      return
      end subroutine alloc_b2mod_ech
c End alloc

c ECH power computation
      subroutine CalcECH(iCv,nCi,iFc,
     &           isn,istra,iwall,pl,geo,mpg,sECH)
c     This routine calculates the ECH power
      use b2mod_b2cmpa
      use b2mod_neutrals_namelist , only:E_fc
      use b2us_plasma
      use b2us_geo
      use b2us_map
      use b2mod_subsys
      implicit none
      type (B2Plasma), intent (in) :: pl
      type (geometry), intent (in) :: geo
      type (mapping), intent (in) :: mpg
      integer,intent(in):: iCv, nCi, iFc, isn, iwall, istra
c     Output variables
      real (kind=R8),intent(out):: sECH

c     Internal variables
      real (kind=R8):: power
      integer :: is

      call subini ('CalcECH')
      power = 1.0E+06 ! MW
      call sfill(sECH, power)

      end subroutine CalcECH

c Dealloc
      subroutine dealloc_b2mod_ech
      implicit none

      deallocate(sECH)

      return
      end subroutine dealloc_b2mod_ech
c End dealloc

      end module b2mod_ech
