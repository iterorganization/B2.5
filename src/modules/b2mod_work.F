      module b2mod_work
      use b2mod_types
      real (kind=R8), allocatable, save ::
     *  conb(:,:,:,:),  conc(:,:,:),                              !srv 13.10.06
     *  cone(:,:,:,:),  coni(:,:,:,:),                            !srv 02.01.07
     *  corpa(:,:,:),  corpo(:,:),
     *  corte(:,:),  corti(:,:),
     *  corua(:,:,:),  corut(:,:),
     *  delne(:,:), delni(:,:,:), delna(:,:,:),  delpo(:,:),  
     *  delte(:,:),  delti(:,:),  delua(:,:,:),
     *  flob(:,:,:),  floe(:,:,:),  floi(:,:,:),
     *  floe_noc(:,:,:),                                          !srv 28.10.02
     *  floi_noc(:,:,:),                                          !srv 28.10.02
     *  pcca(:,:,:,:),  pccm(:,:,:),
     *  avisper(:,:,:),                                           !srv 05.07.00
     *  bvisper(:,:,:),                                           !srv 05.07.00
     *  cvisper(:,:,:),                                           !srv 05.07.00
     *  wrk0(:,:),  wrk1(:,:)
      real (kind=R8), allocatable, save :: 
     &  facdrift_profile(:,:), fac_ExB_profile(:,:), ds(:),
     1	fac_vis_profile(:,:)
      integer, allocatable, save ::
     .  iffboundx(:,:), iffboundy(:,:),
     .  ifenebound(:,:), iftimbound(:,:)

      contains

      subroutine alloc_b2mod_work(nxd,nyd,nsd)
      integer nxd,nyd,nsd

      allocate(conb(-1:nxd,-1:nyd,0:1,0:4))                       !srv 13.10.06
      allocate(conc(-1:nxd,-1:nyd,0:1))
      allocate(cone(-1:nxd,-1:nyd,0:1,0:2))                       !srv 02.01.07
      allocate(coni(-1:nxd,-1:nyd,0:1,0:2))                       !srv 02.01.07
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
      allocate(floe_noc(-1:nxd,-1:nyd,0:1))                       !srv 28.10.02
      allocate(floi_noc(-1:nxd,-1:nyd,0:1))                       !srv 28.10.02
      allocate(pcca(-1:nxd,-1:nyd,0:1,0:nsd-1))
      allocate(pccm(-1:nxd,-1:nyd,0:1))
      allocate(avisper(-1:nxd,-1:nyd,0:1))                        !srv 05.07.00
      allocate(bvisper(-1:nxd,-1:nyd,0:1))                        !srv 05.07.00
      allocate(cvisper(-1:nxd,-1:nyd,0:1))                        !srv 05.07.00
      allocate(wrk0(-1:nxd,-1:nyd))
      allocate(wrk1(-1:nxd,-1:nyd))
      allocate(iffboundx(-1:nxd,-1:nyd))
      allocate(iffboundy(-1:nxd,-1:nyd))
      allocate(ifenebound(-1:nxd,-1:nyd))
      allocate(iftimbound(-1:nxd,-1:nyd))
      allocate(facdrift_profile(-1:nxd,-1:nyd))
      allocate(fac_ExB_profile(-1:nxd,-1:nyd))
      allocate(fac_vis_profile(-1:nxd,-1:nyd))
      allocate(ds(-1:nyd))

      avisper=0.0_R8                                              !srv 05.07.00
      bvisper=0.0_R8                                              !srv 05.07.00
      cvisper=0.0_R8                                              !srv 05.07.00
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
      floe_noc=0.0_R8                                             !srv 28.10.02
      floi_noc=0.0_R8                                             !srv 28.10.02
      pcca=0.0_R8
      pccm=0.0_R8
      wrk0=0.0_R8
      wrk1=0.0_R8
      facdrift_profile=1.0_R8
      fac_ExB_profile=1.0_R8
      fac_vis_profile=1.0_R8
      ds=0.0_R8
      iffboundx=0
      iffboundy=0
      ifenebound=0
      iftimbound=0

      end subroutine alloc_b2mod_work

      end module b2mod_work
