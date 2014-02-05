      module b2mod_work
      use b2mod_types
      real (kind=R8), allocatable, save ::
     *  conb(:,:,:,:),  conc(:,:,:),                                  !srv 13.10.06
     *  cone(:,:,:,:),  coni(:,:,:,:),                                !srv 02.01.07
     *  corpa(:,:,:),  corpo(:,:),
     *  corte(:,:),  corti(:,:),
     *  corua(:,:,:),  corut(:,:),
     *  delna(:,:,:),  delpo(:,:),  delte(:,:),  delti(:,:),
     *  delua(:,:,:),
     *  flob(:,:,:),  floe(:,:,:),  floi(:,:,:),
     *  floe_noc(:,:,:),                                              !srv 28.10.02
     *  floi_noc(:,:,:),                                              !srv 28.10.02
     *  pcca(:,:,:,:),  pccm(:,:,:),
     *  avisper(:,:,:),                                               !srv 05.07.00
     *  bvisper(:,:,:),                                               !srv 05.07.00
     *  cvisper(:,:,:),                                               !srv 05.07.00
     *  wrk0(:,:),  wrk1(:,:)
#ifdef USE_EIRENE
      real (kind=R8), allocatable, save ::                            !IYS 28.03.12 {
     &  sna0_eir_h(:,:,:,:,:),
     &  smo0_eir_h(:,:,:,:,:),
     &  she0_eir_h(:,:,:,:), shi0_eir_h(:,:,:,:),
     &  rfluxa_h(:,:,:,:,:), rfluxm_h(:,:,:,:,:)
! _h is an abbreviation from 'history'
      integer, save :: max_eir_hist_slc = 100                          !IYS 28.03.12 }
#endif

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

      avisper=0                                                       !srv 05.07.00
      bvisper=0                                                       !srv 05.07.00
      cvisper=0                                                       !srv 05.07.00
      conb=0
      conc=0
      cone=0
      coni=0
      corpa=0
      corpo=0
      corte=0
      corti=0
      corua=0
      corut=0
      delna=0
      delpo=0
      delte=0
      delti=0
      delua=0
      flob=0
      floe=0
      floi=0
      floe_noc=0                                                     !srv 28.10.02
      floi_noc=0                                                     !srv 28.10.02
      pcca=0
      pccm=0
      wrk0=0
      wrk1=0

      end subroutine alloc_b2mod_work

#ifdef USE_EIRENE  
      subroutine alloc_hist_arrays(nxd,nyd,nsd)                     !IYS 28.03.12 {

      implicit none
      integer nxd,nyd,nsd
#     include "DIMENSIONS.F"
#     include <COUPLE/KOPPLDIM.F>

! to store the time history of sources from EIRENE
      allocate(sna0_eir_h(-1:nxd,-1:nyd,0:1,0:nsd-1,max_eir_hist_slc))
      allocate(smo0_eir_h(-1:nxd,-1:nyd,0:3,0:nsd-1,max_eir_hist_slc))
      allocate(she0_eir_h(-1:nxd,-1:nyd,0:3,max_eir_hist_slc))
      allocate(shi0_eir_h(-1:nxd,-1:nyd,0:3,max_eir_hist_slc))
      allocate(rfluxa_h(0:ndxp,0:ndyp,natm,1,max_eir_hist_slc))
      allocate(rfluxm_h(0:ndxp,0:ndyp,natm,1,max_eir_hist_slc))
      sna0_eir_h = 0.0_R8
      smo0_eir_h = 0.0_R8
      she0_eir_h = 0.0_R8
      shi0_eir_h = 0.0_R8 
      rfluxa_h = 0.0_R8  
      rfluxm_h = 0.0_R8    

      end subroutine alloc_hist_arrays                             !IYS 28.03.12 }      
#endif            

      end module b2mod_work
