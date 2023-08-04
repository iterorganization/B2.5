c**************************************************************************************
c********************  NEW KU Leuven module for ECH/ECWC plasmas  *********************
c**************************************************************************************
!mb   06/07/2017
c   This module is intended to simulate power depositions in Electron Cyclotron Heated
c     plasmas. To turn on the model, the flag b2stbr_ech has to be set to '1'. This
c     module computes an additional heating source due to EC heating injected at the
c     cyclotron frequency
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

      logical, save :: recycle_initialised
      data recycle_initialised /.false./

      real (kind=R8),save::uold=0.0_R8,told=0.0_R8,mold=0.0_R8
      integer, save :: nm, nt, nshpottau,nx1,ny1,ns1

c     Recycle_index and target_index refer for each species and boundary in which integer
c     the corresponding recycling coefficients are stored.
      integer, allocatable, save ::
     *  recycle_index(:),target_index(:)
c     Boolean to indicate whether trim data is present for species/material combination
      logical, allocatable, save ::
     *  trimrecyc_present(:,:)
c      real (kind=R8), allocatable, save ::

      real (kind=R8),save::int0l,int0r,
     &    int1l,int1r,int2l,int2r,
     &    int3l,int3r,int4l,int4r,
     &    int5l,int5r,int6l,int6r

      contains

      subroutine alloc_b2mod_recycle(nFcd,nsd)
      implicit none
      integer,intent(in) :: nFcd,nsd

!      allocate()

      return
      end subroutine alloc_b2mod_recycle

      end module b2mod_ech
