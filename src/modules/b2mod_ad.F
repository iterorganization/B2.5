      module b2mod_ad
      use b2mod_types
#     include "DIMENSIONS.F"
      implicit none
!     This module is intended to collect variables within B2.5 which need to be 'visible'
!     to calling routines as required by adjoint AD with Tapenade

      integer NSDMAX
      parameter (NSDMAX=DEF_NSD)
      integer CVREGMAX
      parameter (CVREGMAX=8)
      integer nncf
      parameter (nncf=8)

      character(7), save :: my_out_folder
#ifdef WG_DEBUG
      character(8), save :: my_out_sf(0:1)
#endif

      integer, save :: ntstep_b2wall = 0 
      real (kind=R8), allocatable, save ::
     *  old_erosion(:,:), old_deposition(:,:)

      integer, save :: b2sral_elm_count = 0

      real (kind=R8), save :: b2sian_eps=1.0e-60_R8

      real (kind=R8), save :: b2tqna_keps_eps = 1.0e-10

      logical, save :: b2npmo_solvedum(0:CVREGMAX) = .true.

      logical, save :: b2news_solving(4) = .true.    !srv 06.07.00

      real (kind=R8), save :: cfnorm(nncf), vold(nncf)
      real (kind=R8), allocatable, save :: TMP1(:), TMP2(:)
      real (kind=R8), allocatable, save :: rrdata(:,:), vvdata(:,:)

      logical, save :: test_residual = .false. 
      integer, save :: b2mndt_itcnt = 0

      character, save :: fmt(0:5)*32
      data fmt /'(1p,5(1x,e25.17e3))', '(1p,6(1x,e21.13e3))',
     &  '(1p,7(1x,e17.9e3))', '(1p,8(1x,e15.7e3))',
     &  '(1p,9(1x,e13.5e3))', '(1p,10(1x,e12.4e3))'/

! ncall variables from transport routines
      integer, save :: ncall_b2tcpa = 0, ncall_b2tfcc = 0,
     & ncall_b2tfch = 0, ncall_b2tfed = 0, ncall_b2tfhe = 0,
     & ncall_b2tfhi = 0, ncall_b2tfnb = 0, ncall_b2tfrn = 0,
     & ncall_b2tiner = 0, ncall_b2tlh0 = 0, ncall_b2tlhe = 0,
     & ncall_b2tlmv = 0, ncall_b2tlnl = 0, ncall_b2tqca = 0,
     & ncall_b2tqce = 0, ncall_b2tqin = 0, ncall_b2tqna = 0,
     & ncall_b2tral = 0, ncall_transp_keps = 0, ncall_b2trcl = 0,
     & ncall_b2treq = 0, ncall_b2trno = 0, ncall_b2trql = 0,
! ncall variables from sources routines
     & ncall_b2siav = 0, ncall_b2sigp = 0, ncall_b2sihs_ = 0,
     & ncall_b2sikt = 0, ncall_b2spcx = 0, ncall_b2spel = 0,
     & ncall_b2sqcx = 0, ncall_b2sqel = 0, ncall_b2sral = 0,
     & ncall_b2srdt = 0, ncall_b2srsm = 0, ncall_b2srst = 0,
     & ncall_b2stbc = 0, ncall_b2stbc_phys = 0, ncall_b2stbr = 0,
     & ncall_b2stbr_phys = 0, ncall_b2stcx = 0, ncall_b2stel = 0,
     & ncall_b2sicf,
! ncall variables from equation routines
     & ncall_b2npco = 0, ncall_b2news_ = 0, ncall_b2npht = 0,
     & ncall_b2npmo = 0, ncall_b2nppo = 0, ncall_b2nxdp = 0,
     & ncall_b2nxdu = 0, ncall_b2nxfc = 0,
! ncall variables from solvers routines
     & ncall_b2upco = 0, ncall_b2upht = 0, ncall_b2uppo = 0, 
     & ncall_b2urmo = 0,ncall_b2ursc = 0, ncall_b2ursd = 0, 
     & ncall_b2usco = 0, ncall_b2usht = 0, ncall_b2usmo = 0,
     & ncall_b2uspo = 0,
! ncall variables from other routines
     & ncall_b2xehx = 0, ncall_b2xehy = 0, ncall_b2mndt = 0

      end module b2mod_ad
