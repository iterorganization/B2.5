      module b2mod_ad
      use b2mod_types
      use b2mod_numerics_namelist
     & , only : cvregmax
      use b2mod_user_namelist
     & , only : nncf
      implicit none
!     This module is intended to collect variables within B2.5 which need to be 'visible'
!     to calling routines as required by adjoint AD with Tapenade

      character(7), save :: my_out_folder
#ifdef WG_DEBUG
      character(8), save :: my_out_sf(0:1)
#endif

      integer, save :: ntstep_b2wall !CHECK!
      real (kind=R8), allocatable, save ::
     *  old_erosion(:,:), old_deposition(:,:)

      integer, save :: b2sral_elm_count = 0

      real (kind=R8), save :: b2sian_eps=1.0e-60_R8

      real (kind=R8), save :: b2tqna_keps_eps = 1.0e-10

      logical, save :: b2npmo_solvedum(0:CVREGMAX) = .true.

      logical, save :: b2news_solving(4) = .true.    !srv 06.07.00

      real (kind=R8), save :: cfnorm(nncf), vold(nncf)
      real (kind=R8), allocatable, save :: ffdata(:,:)

      logical, save :: test_residual = .false. 
      integer, save :: b2mndt_itcnt = 0
      end module b2mod_ad
