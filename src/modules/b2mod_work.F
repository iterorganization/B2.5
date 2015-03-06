      module b2mod_work
      use b2mod_types
      real (kind=R8), allocatable, save ::
     *  conb(:,:,:,:),  conc(:,:,:),                                  !srv 13.10.06
     *  cone(:,:,:,:),  coni(:,:,:,:),                                !srv 02.01.07
     *  corpa(:,:,:),  corpo(:,:),
     *  corte(:,:),  corti(:,:),
     *  corua(:,:,:),  corut(:,:),
     *  delne(:,:), delni(:,:,:), delna(:,:,:),  delpo(:,:),  
     *  delte(:,:),  delti(:,:),  delua(:,:,:),
     *  flob(:,:,:),  floe(:,:,:),  floi(:,:,:),
     *  floe_noc(:,:,:),                                              !srv 28.10.02
     *  floi_noc(:,:,:),                                              !srv 28.10.02
     *  pcca(:,:,:,:),  pccm(:,:,:),
     *  avisper(:,:,:),                                               !srv 05.07.00
     *  bvisper(:,:,:),                                               !srv 05.07.00
     *  cvisper(:,:,:),                                               !srv 05.07.00
     *  wrk0(:,:),  wrk1(:,:)
      real (kind=R8), allocatable, save :: 
     &  facdrift_profile(:,:), fac_ExB_profile(:,:), ds(:),
     1	fac_vis_profile(:,:)
#ifdef USE_EIRENE
      real (kind=R8), allocatable, save ::                            !IYS 28.03.12 {
     &  sna0_eir_h(:,:,:,:,:),
     &  smo0_eir_h(:,:,:,:,:),
     &  she0_eir_h(:,:,:,:), shi0_eir_h(:,:,:,:),
     &  rfluxa_h(:,:,:,:,:), rfluxm_h(:,:,:,:,:)
! _h is an abbreviation from 'history'
      integer, save :: max_eir_hist_slc = 100                          !IYS 28.03.12 }
#endif
      integer, allocatable, save ::
     .  iffboundx(:,:), iffboundy(:,:),
     .  ifenebound(:,:), iftimbound(:,:)

      contains

      subroutine alloc_b2mod_work(nxd,nyd,nsd)
      integer nxd,nyd,nsd

      allocate(conb(-1:nxd,-1:nyd,0:1,0:4))                           !srv 13.10.06
      allocate(conc(-1:nxd,-1:nyd,0:1))
      allocate(cone(-1:nxd,-1:nyd,0:1,0:2))                           !srv 02.01.07
      allocate(coni(-1:nxd,-1:nyd,0:1,0:2))                           !srv 02.01.07
      allocate(corpa(-1:nxd,-1:nyd,0:nsd-1))
      allocate(corpo(-1:nxd,-1:nyd))
      allocate(corte(-1:nxd,-1:nyd))
      allocate(corti(-1:nxd,-1:nyd))
      allocate(corua(-1:nxd,-1:nyd,0:nsd-1))
      allocate(corut(-1:nxd,-1:nyd))
      allocate(delne(-1:nxd,-1:nyd))
      allocate(delni(-1:nxd,-1:nyd,0:1))
      allocate(delna(-1:nxd,-1:nyd,0:nsd-1))
      allocate(delpo(-1:nxd,-1:nyd))
      allocate(delte(-1:nxd,-1:nyd))
      allocate(delti(-1:nxd,-1:nyd))
      allocate(delua(-1:nxd,-1:nyd,0:nsd-1))
      allocate(flob(-1:nxd,-1:nyd,0:1))
      allocate(floe(-1:nxd,-1:nyd,0:1))
      allocate(floi(-1:nxd,-1:nyd,0:1))
      allocate(floe_noc(-1:nxd,-1:nyd,0:1))                           !srv 28.10.02
      allocate(floi_noc(-1:nxd,-1:nyd,0:1))                           !srv 28.10.02
      allocate(pcca(-1:nxd,-1:nyd,0:1,0:nsd-1))
      allocate(pccm(-1:nxd,-1:nyd,0:1))
      allocate(avisper(-1:nxd,-1:nyd,0:1))                            !srv 05.07.00
      allocate(bvisper(-1:nxd,-1:nyd,0:1))                            !srv 05.07.00
      allocate(cvisper(-1:nxd,-1:nyd,0:1))                            !srv 05.07.00
      allocate(wrk0(-1:nxd,-1:nyd))
      allocate(wrk1(-1:nxd,-1:nyd))
      allocate(iffboundx(-1:nxd,-1:nyd))
      allocate(iffboundy(-1:nxd,-1:nyd))
      allocate(ifenebound(-1:nxd,-1:nyd))
      allocate(iftimbound(-1:nxd,-1:nyd))

      avisper=0.0_R8                                                       !srv 05.07.00
      bvisper=0.0_R8                                                       !srv 05.07.00
      cvisper=0.0_R8                                                       !srv 05.07.00
      conb=0.0_R8
      conc=0.0_R8
      cone=0.0_R8
      coni=0.0_R8
      corpa=0.0_R8
      corpo=0.0_R8
      corte=0.0_R8
      corti=0.0_R8
      corua=0.0_R8
      corut=0.0_R8
      delne=0.0_R8
      delni=0.0_R8
      delna=0.0_R8
      delpo=0.0_R8
      delte=0.0_R8
      delti=0.0_R8
      delua=0.0_R8
      flob=0.0_R8
      floe=0.0_R8
      floi=0.0_R8
      floe_noc=0.0_R8                                                     !srv 28.10.02
      floi_noc=0.0_R8                                                     !srv 28.10.02
      pcca=0.0_R8
      pccm=0.0_R8
      wrk0=0.0_R8
      wrk1=0.0_R8
      iffboundx=0
      iffboundy=0
      ifenebound=0
      iftimbound=0

      end subroutine alloc_b2mod_work

#ifdef USE_EIRENE  
      subroutine alloc_hist_arrays(nxd,nyd,nsd_local,nhistd)      !IYS 28.03.12 {

      use eirmod_parmmod
      use eirmod_braeir
      implicit none
      integer nxd,nyd,nsd_local,nhistd
#     include "DIMENSIONS.F"

! to store the time history of sources from EIRENE
      allocate(sna0_eir_h(-1:nxd,-1:nyd,0:1,0:nsd_local-1,nhistd))
      allocate(smo0_eir_h(-1:nxd,-1:nyd,0:3,0:nsd_local-1,nhistd))
      allocate(she0_eir_h(-1:nxd,-1:nyd,0:3,nhistd))
      allocate(shi0_eir_h(-1:nxd,-1:nyd,0:3,nhistd))
      allocate(rfluxa_h(0:ndxp,0:ndyp,natm,1,nhistd))
      allocate(rfluxm_h(0:ndxp,0:ndyp,natm,1,nhistd))
      call clear_hist_arrays
      end subroutine alloc_hist_arrays

      subroutine clear_hist_arrays                                 !csw 24feb2014
      implicit none
      sna0_eir_h = 0.0_R8
      smo0_eir_h = 0.0_R8
      she0_eir_h = 0.0_R8
      shi0_eir_h = 0.0_R8 
      rfluxa_h = 0.0_R8  
      rfluxm_h = 0.0_R8    
      end subroutine clear_hist_arrays                             !IYS 28.03.12 }      
#endif            

      end module b2mod_work
