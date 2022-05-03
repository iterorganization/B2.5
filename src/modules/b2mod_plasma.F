      module b2mod_plasma
      use b2mod_types
      use b2mod_version
      implicit none

      real (kind=R8), allocatable, save ::
     *  fch(:,:,:,:),  fch0(:,:,:,:),  fchp(:,:,:,:),  fch_32(:,:,:,:),
     *  fch_52(:,:,:,:),
     *  fna_32(:,:,:,:,:),  fna_52(:,:,:,:,:),
     *  fni_32(:,:,:,:),    fni_52(:,:,:,:),
     *  fne_32(:,:,:,:),    fne_52(:,:,:,:),
     *  fhe(:,:,:,:),  fhe0(:,:,:,:),  fhep(:,:,:,:),  fhet(:,:,:,:),
     *  fhi(:,:,:,:),  fhi0(:,:,:,:),  fhip(:,:,:,:),  fhit(:,:,:,:),
     *  fna(:,:,:,:,:),  fna0(:,:,:,:,:),  fnap(:,:,:,:,:),
     *  fne(:,:,:,:),    fni(:,:,:,:),
     *  na(:,:,:),  na0(:,:,:),  nap(:,:,:),
     *  ne(:,:),  ne0(:,:),  ne2(:,:),  nep(:,:),
     *  ni(:,:,:),  ni0(:,:,:), nip(:,:,:),                              !srv 07.07.21
     *  pb(:,:),
     *  po(:,:),  po0(:,:),  pop(:,:),
     *  te(:,:),  te0(:,:),  tep(:,:),
     *  ti(:,:),  ti0(:,:),  tip(:,:),
     *  ua(:,:,:),  ua0(:,:,:),  uap(:,:,:), ue(:,:),
     *  uadia(:,:,:,:,:),  fchdia(:,:,:,:),
     *  fch_p(:,:,:,:), fch_pi_c(:,:,:,:), fch_pi_f(:,:,:,:),            !srv 10.06.08 { 13.01.17
     *  fhm(:,:,:,:,:), fhp(:,:,:,:,:), fhj(:,:,:,:), fht(:,:,:,:),
     *  fhe_mdf(:,:,:,:), fhePSch(:,:,:,:),
     *  fhi_mdf(:,:,:,:), fhiPSch(:,:,:,:),                              !srv 13.10.06
     *  fna_mdf(:,:,:,:,:), fnaPSch(:,:,:,:,:),
     *  fna_fcor(:,:,:,:,:), fna_he(:,:,:,:,:),                          !srv 01.10.99 21.01.00
     *  fna_eir(:,:,:,:,:), fne_eir(:,:,:,:), fhe_eir(:,:,:,:),          !srv 14.07.10
     *  fhi_eir(:,:,:,:),
     *  fna_32nd(:,:,:,:,:), fna_52nd(:,:,:,:,:)                         !srv 18.06.08

      real (kind=R8), allocatable, save ::
     *  fchvispar(:,:,:,:), fchvisper(:,:,:,:), fchvisq(:,:,:,:),        !srv 11.07.99
     *  fchinert(:,:,:,:), fchin(:,:,:,:), fchanml(:,:,:,:),             !srv 07.07.99 24.03.07 10.06.08 }
     *  fchstoch(:,:,:),                                                 !srv 24.08.21
     *  vadia(:,:,:,:,:), vedia(:,:,:,:), wadia(:,:,:,:,:),              !xpb
     *  vaecrb(:,:,:,:,:), veecrb(:,:,:,:),
     *  wedia(:,:,:,:), ve(:,:),                                         !srv 14.07.10
     *  fna_nodrift(:,:,:,:,:),
     *  kinrgy(:,:,:), kinrgy0(:,:,:), kinrgyp(:,:,:),                   !xpb
     *  fmo(:,:,:,:,:)
      real (kind=R8), allocatable, save ::                               !srv 08.09.21 {
     &  fchanml_a(:,:,:,:,:),
     &  fchinert_a(:,:,:,:,:),
     &  fchvispar_a(:,:,:,:,:),
     &  fchvisper_a(:,:,:,:,:),
     &  fchvisq_a(:,:,:,:,:)                                             !srv 08.09.21 }
      real (kind=R8), allocatable, save ::
     *  facdrift(:,:), fac_ExB(:,:), fac_vis(:,:)

      integer, save :: facdrift_ncall
      real (kind=R8), allocatable, save ::
     &  facdrift_profile(:,:),
     &  fac_ExB_profile(:,:), fac_vis_profile(:,:)
      real (kind=R8), save ::
     *  facdrift_scalar, fac_ExB_scalar, fac_vis_scalar,
     *  facdrift_start, facExB_start, facvis_start,
     *  facdrift_target, facExB_target, facvis_target,
     *  facdrift_inc, facExB_inc, facvis_inc,
     *  facdrift_dec, facExB_dec, facvis_dec
      data facdrift_ncall/0/, facdrift_scalar/0.0_R8/,
     *  fac_ExB_scalar/0.0_R8/, fac_vis_scalar/0.0_R8/,
     *  facdrift_start/0.0_R8/, facdrift_target/0.0_R8/,
     *  facExB_start/0.0_R8/, facExB_target/0.0_R8/,
     *  facvis_start/0.0_R8/, facvis_target/0.0_R8/,
     *  facdrift_inc/1.0_R8/, facdrift_dec/0.0_R8/,
     *  facExB_inc/1.0_R8/, facExB_dec/0.0_R8/,
     *  facvis_inc/1.0_R8/, facvis_dec/0.0_R8/

      contains

      subroutine alloc_b2mod_plasma(nxd,nyd,nsd)
      implicit none
      integer nxd,nyd,nsd

      allocate(fch(-1:nxd,-1:nyd,0:1,0:1))
      allocate(fch0(-1:nxd,-1:nyd,0:1,0:1))
      allocate(fchp(-1:nxd,-1:nyd,0:1,0:1))
      allocate(fch_32(-1:nxd,-1:nyd,0:1,0:1))
      allocate(fch_52(-1:nxd,-1:nyd,0:1,0:1))
      allocate(fna_32(-1:nxd,-1:nyd,0:1,0:1,0:nsd-1))
      allocate(fna_52(-1:nxd,-1:nyd,0:1,0:1,0:nsd-1))
      allocate(fni_32(-1:nxd,-1:nyd,0:1,0:1))
      allocate(fni_52(-1:nxd,-1:nyd,0:1,0:1))
      allocate(fne_32(-1:nxd,-1:nyd,0:1,0:1))
      allocate(fne_52(-1:nxd,-1:nyd,0:1,0:1))
      allocate(fna_32nd(-1:nxd,-1:nyd,0:1,0:1,0:nsd-1))            !srv 18.06.08
      allocate(fna_52nd(-1:nxd,-1:nyd,0:1,0:1,0:nsd-1))            !srv 18.06.08
      allocate(fhe(-1:nxd,-1:nyd,0:1,0:1))
      allocate(fhe0(-1:nxd,-1:nyd,0:1,0:1))
      allocate(fhep(-1:nxd,-1:nyd,0:1,0:1))
      allocate(fhet(-1:nxd,-1:nyd,0:1,0:1))
      allocate(fhe_mdf(-1:nxd,-1:nyd,0:1,0:1))                     !srv 13.10.06
      allocate(fhePSch(-1:nxd,-1:nyd,0:1,0:1))                     !srv 13.10.06
      allocate(fhi(-1:nxd,-1:nyd,0:1,0:1))
      allocate(fhi0(-1:nxd,-1:nyd,0:1,0:1))
      allocate(fhip(-1:nxd,-1:nyd,0:1,0:1))
      allocate(fhit(-1:nxd,-1:nyd,0:1,0:1))
      allocate(fhi_mdf(-1:nxd,-1:nyd,0:1,0:1))                     !srv 13.10.06
      allocate(fhiPSch(-1:nxd,-1:nyd,0:1,0:1))                     !srv 13.10.06
      allocate(fna_he(-1:nxd,-1:nyd,0:1,0:1,0:nsd-1))              !srv 21.01.00
      allocate(fna_fcor(-1:nxd,-1:nyd,0:1,0:1,0:nsd-1))            !srv 01.10.99
      allocate(fnaPSch(-1:nxd,-1:nyd,0:1,0:1,0:nsd-1))             !srv 13.10.06
      allocate(fna_mdf(-1:nxd,-1:nyd,0:1,0:1,0:nsd-1))             !srv 13.10.06
      allocate(fna_eir(-1:nxd,-1:nyd,0:1,0:1,0:nsd-1))             !srv 14.07.10 {
      allocate(fne_eir(-1:nxd,-1:nyd,0:1,0:1))
      allocate(fhe_eir(-1:nxd,-1:nyd,0:1,0:1))
      allocate(fhi_eir(-1:nxd,-1:nyd,0:1,0:1))
      fna_he = 0.0_R8
      fnaPSch = 0.0_R8
      fhePSch = 0.0_R8
      fhiPSch = 0.0_R8
      fna_mdf = 0.0_R8
      fhe_mdf = 0.0_R8
      fhi_mdf = 0.0_R8
      fna_eir = 0.0_R8
      fne_eir = 0.0_R8
      fhe_eir = 0.0_R8
      fhi_eir = 0.0_R8
      fna_fcor = 0.0_R8                                            !srv 14.07.10 }
      allocate(fna(-1:nxd,-1:nyd,0:1,0:1,0:nsd-1))
      allocate(fna0(-1:nxd,-1:nyd,0:1,0:1,0:nsd-1))
      allocate(fnap(-1:nxd,-1:nyd,0:1,0:1,0:nsd-1))
      allocate(fne(-1:nxd,-1:nyd,0:1,0:1))
      allocate(fni(-1:nxd,-1:nyd,0:1,0:1))
      allocate(fhm(-1:nxd,-1:nyd,0:1,0:1,0:nsd-1))
      allocate(fhp(-1:nxd,-1:nyd,0:1,0:1,0:nsd-1))
      allocate(fhj(-1:nxd,-1:nyd,0:1,0:1))
      allocate(fht(-1:nxd,-1:nyd,0:1,0:1))
      allocate(na(-1:nxd,-1:nyd,0:nsd-1))
      allocate(na0(-1:nxd,-1:nyd,0:nsd-1))
      allocate(nap(-1:nxd,-1:nyd,0:nsd-1))
      allocate(ne(-1:nxd,-1:nyd))
      allocate(ne0(-1:nxd,-1:nyd))
      allocate(ne2(-1:nxd,-1:nyd))
      allocate(nep(-1:nxd,-1:nyd))
      allocate(ni(-1:nxd,-1:nyd,0:1))
      allocate(ni0(-1:nxd,-1:nyd,0:1))
      allocate(nip(-1:nxd,-1:nyd,0:1))                               !srv 07.07.21
      allocate(pb(-1:nxd,-1:nyd))
      allocate(po(-1:nxd,-1:nyd))
      allocate(po0(-1:nxd,-1:nyd))
      allocate(pop(-1:nxd,-1:nyd))
      allocate(te(-1:nxd,-1:nyd))
      allocate(te0(-1:nxd,-1:nyd))
      allocate(tep(-1:nxd,-1:nyd))
      allocate(ti(-1:nxd,-1:nyd))
      allocate(ti0(-1:nxd,-1:nyd))
      allocate(tip(-1:nxd,-1:nyd))
      allocate(ue(-1:nxd,-1:nyd))
      allocate(ua(-1:nxd,-1:nyd,0:nsd-1))
      allocate(ua0(-1:nxd,-1:nyd,0:nsd-1))
      allocate(uap(-1:nxd,-1:nyd,0:nsd-1))
      allocate(uadia(-1:nxd,-1:nyd,0:1,0:1,0:nsd-1))
      allocate(fch_pi_c(-1:nxd,-1:nyd,0:1,0:1))                      !srv 13.01.17
      allocate(fch_pi_f(-1:nxd,-1:nyd,0:1,0:1))                      !srv 13.01.17
      allocate(fchdia(-1:nxd,-1:nyd,0:1,0:1))
      fchdia=0.0_R8
      allocate(fch_p(-1:nxd,-1:nyd,0:1,0:1))                         !xpb srv 10.06.08
      allocate(fchvispar(-1:nxd,-1:nyd,0:1,0:1))                     !xpb
      fchvispar=0.0_R8
      allocate(fchvisper(-1:nxd,-1:nyd,0:1,0:1))                     !xpb
      fchvisper=0.0_R8
      allocate(fchvisq(-1:nxd,-1:nyd,0:1,0:1))                       !srv 11.07.99
      fchvisq=0.0_R8
      allocate(fchinert(-1:nxd,-1:nyd,0:1,0:1))                      !srv 07.07.99
      fchinert=0.0_R8
      allocate(fchin(-1:nxd,-1:nyd,0:1,0:1))                         !srv 30.09.02
      fchin=0.0_R8
      allocate(fchanml(-1:nxd,-1:nyd,0:1,0:1))                       !srv 24.03.07
      allocate(fchstoch(-1:nxd,-1:nyd,0:1))                          !srv 24.08.21
      fchstoch = 0.0_R8                                              !srv 24.08.21
      allocate(fchanml_a(-1:nxd,-1:nyd,0:1,0:1,0:nsd-1))             !srv 08.09.21 {
      allocate(fchinert_a(-1:nxd,-1:nyd,0:1,0:1,0:nsd-1))
      allocate(fchvispar_a(-1:nxd,-1:nyd,0:1,0:1,0:nsd-1))
      allocate(fchvisper_a(-1:nxd,-1:nyd,0:1,0:1,0:nsd-1))
      allocate(fchvisq_a(-1:nxd,-1:nyd,0:1,0:1,0:nsd-1))             !srv 08.09.21 }
      allocate(vadia(-1:nxd,-1:nyd,0:1,0:1,0:nsd-1))                 !xpb
      vadia=0.0_R8
      allocate(vedia(-1:nxd,-1:nyd,0:1,0:1))                         !xpb
      allocate(wadia(-1:nxd,-1:nyd,0:1,0:1,0:nsd-1))                 !xpb
      allocate(wedia(-1:nxd,-1:nyd,0:1,0:1))                         !srv 14.07.10
      allocate(ve(-1:nxd,-1:nyd))                                    !srv 14.07.10
      allocate(vaecrb(-1:nxd,-1:nyd,0:1,0:1,0:nsd-1))                !xpb
      vaecrb=0.0_R8
      allocate(veecrb(-1:nxd,-1:nyd,0:1,0:1))                        !xpb
      allocate(fna_nodrift(-1:nxd,-1:nyd,0:1,0:1,0:nsd-1))           !xpb
      fna_nodrift=0.0_R8
      allocate(kinrgy(-1:nxd,-1:nyd,0:nsd-1))                        !xpb
      allocate(kinrgy0(-1:nxd,-1:nyd,0:nsd-1))                       !xpb
      allocate(kinrgyp(-1:nxd,-1:nyd,0:nsd-1))                       !xpb
      allocate(fmo(-1:nxd,-1:nyd,0:1,0:1,0:nsd-1))
      allocate(facdrift(-1:nxd,-1:nyd))                              !dpc
      facdrift=0.0_R8
      allocate(fac_ExB(-1:nxd,-1:nyd))                               !dpc
      fac_ExB=0.0_R8
      allocate(fac_vis(-1:nxd,-1:nyd))                               !dpc
      fac_vis=0.0_R8
      allocate(facdrift_profile(-1:nxd,-1:nyd))
      allocate(fac_ExB_profile(-1:nxd,-1:nyd))
      allocate(fac_vis_profile(-1:nxd,-1:nyd))
      facdrift_profile=1.0_R8
      fac_ExB_profile=1.0_R8
      fac_vis_profile=1.0_R8

      return
      end subroutine alloc_b2mod_plasma

      subroutine dealloc_b2mod_plasma
      implicit none

      deallocate(fch)
      deallocate(fch0)
      deallocate(fchp)
      deallocate(fch_32)
      deallocate(fch_52)
      deallocate(fna_32)
      deallocate(fna_52)
      deallocate(fna_32nd)                                         !srv 18.06.08
      deallocate(fna_52nd)                                         !srv 18.06.08
      deallocate(fni_32)
      deallocate(fni_52)
      deallocate(fne_32)
      deallocate(fne_52)
      deallocate(fhe)
      deallocate(fhe0)
      deallocate(fhep)
      deallocate(fhet)
      deallocate(fhi)
      deallocate(fhi0)
      deallocate(fhip)
      deallocate(fhit)
      deallocate(fna)
      deallocate(fna0)
      deallocate(fnap)
      deallocate(fne)
      deallocate(fni)
      deallocate(fhm)
      deallocate(fhp)
      deallocate(fhj)
      deallocate(fht)
      deallocate(na)
      deallocate(na0)
      deallocate(nap)
      deallocate(ne)
      deallocate(ne0)
      deallocate(ne2)
      deallocate(nep)
      deallocate(ni)
      deallocate(ni0)
      deallocate(nip)                                              !srv 07.07.21
      deallocate(pb)
      deallocate(po)
      deallocate(po0)
      deallocate(pop)
      deallocate(te)
      deallocate(te0)
      deallocate(tep)
      deallocate(ti)
      deallocate(ti0)
      deallocate(tip)
      deallocate(ue)
      deallocate(ua)
      deallocate(ua0)
      deallocate(uap)
      deallocate(uadia)
      deallocate(fchdia)
      deallocate(fch_p)
      deallocate(fchvispar)
      deallocate(fchvisper)
      deallocate(fchvisq)                           !srv 11.07.99
      deallocate(fchinert)                          !srv 07.07.99
      deallocate(fchin)                             !srv 30.09.02
      deallocate(fchanml)                           !srv 24.03.07
      deallocate(fchstoch)                          !srv 24.08.21
      deallocate(fchanml_a)                         !srv 08.09.21 {
      deallocate(fchinert_a)
      deallocate(fchvispar_a)
      deallocate(fchvisper_a)
      deallocate(fchvisq_a)                         !srv 08.09.21 }
      deallocate(vadia)
      deallocate(vedia)
      deallocate(wadia)
      deallocate(vaecrb)
      deallocate(veecrb)
      deallocate(fna_nodrift)
      deallocate(kinrgy)
      deallocate(kinrgy0)
      deallocate(kinrgyp)
      deallocate(fmo)
      deallocate(facdrift)
      deallocate(fac_ExB)
      deallocate(fac_vis)
      deallocate(fhe_mdf)                           !srv 13.10.06
      deallocate(fhePSch)                           !srv 13.10.06
      deallocate(fhi_mdf)                           !srv 13.10.06
      deallocate(fhiPSch)                           !srv 13.10.06
      deallocate(fna_he)                            !srv 21.01.00
      deallocate(fna_fcor)                          !srv 01.10.99
      deallocate(fnaPSch)                           !srv 13.10.06
      deallocate(fna_mdf)                           !srv 13.10.06
      deallocate(fch_pi_c)                          !srv 13.01.17
      deallocate(fch_pi_f)                          !srv 13.01.17
      deallocate(fna_eir)                           !srv 14.07.10 {
      deallocate(fne_eir)
      deallocate(fhe_eir)
      deallocate(fhi_eir)
      deallocate(wedia)
      deallocate(ve)                                !srv 14.07.10 }
      deallocate(facdrift_profile)
      deallocate(fac_ExB_profile)
      deallocate(fac_vis_profile)

      return
      end subroutine dealloc_b2mod_plasma

      subroutine read_b2mod_plasma(nxd,nyd,nsd,nread)
      use b2mod_constants
      use b2mod_b2cmpa
      implicit none
      integer nxd,nyd,nsd,nread
      integer is
*   ..local arrays for reading fluxes of old dimensions
      real (kind=R8) :: floc(-1:nxd,-1:nyd,0:1),
     &  flocs(-1:nxd,-1:nyd,0:1,0:nsd-1)
      external b2xpni, cfrure, cfrure_opt

      if(b2fplasma_version.ge.'03.001.000') then
        call cfrure (nread, size(fch),  fch,  'fch')
        call cfrure (nread, size(fch0), fch0, 'fch0')
        call cfrure (nread, size(fchp), fchp, 'fchp')
        call cfrure (nread, size(fhe),  fhe,  'fhe')
        call cfrure (nread, size(fhe0), fhe0, 'fhe0')
        call cfrure (nread, size(fhep), fhep, 'fhep')
        call cfrure (nread, size(fhet), fhet, 'fhet')
        call cfrure (nread, size(fhi),  fhi,  'fhi')
        call cfrure (nread, size(fhi0), fhi0, 'fhi0')
        call cfrure (nread, size(fhip), fhip, 'fhip')
        call cfrure (nread, size(fhit), fhit, 'fhit')
        call cfrure (nread, size(fna),  fna,  'fna' )
        call cfrure (nread, size(fna0), fna0, 'fna0')
        call cfrure (nread, size(fnap), fnap, 'fnap')
        call cfrure (nread, size(fne),  fne,  'fne' )
        call cfrure (nread, size(fni),  fni,  'fni' )
      else
        call cfrure (nread, size(floc), floc, 'fch')
        fch = 0.0_R8
        fch(:,:,0,0) = floc(:,:,0)
        fch(:,:,1,1) = floc(:,:,1)
        call cfrure (nread, size(floc), floc, 'fch0')
        fch0 = 0.0_R8
        fch0(:,:,0,0) = floc(:,:,0)
        fch0(:,:,1,1) = floc(:,:,1)
        call cfrure (nread, size(floc), floc, 'fchp')
        fchp = 0.0_R8
        fchp(:,:,0,0) = floc(:,:,0)
        fchp(:,:,1,1) = floc(:,:,1)
        call cfrure (nread, size(floc), floc, 'fhe')
        fhe = 0.0_R8
        fhe(:,:,0,0) = floc(:,:,0)
        fhe(:,:,1,1) = floc(:,:,1)
        call cfrure (nread, size(floc), floc, 'fhe0')
        fhe0 = 0.0_R8
        fhe0(:,:,0,0) = floc(:,:,0)
        fhe0(:,:,1,1) = floc(:,:,1)
        call cfrure (nread, size(floc), floc, 'fhep')
        fhep = 0.0_R8
        fhep(:,:,0,0) = floc(:,:,0)
        fhep(:,:,1,1) = floc(:,:,1)
        if(b2fplasma_version.gt.'01.001.003') then
          call cfrure (nread, size(fhet), fhet, 'fhet')
        else
          fhet = fhe
        endif
        call cfrure (nread, size(floc), floc, 'fhi')
        fhi = 0.0_R8
        fhi(:,:,0,0) = floc(:,:,0)
        fhi(:,:,1,1) = floc(:,:,1)
        call cfrure (nread, size(floc), floc, 'fhi0')
        fhi0 = 0.0_R8
        fhi0(:,:,0,0) = floc(:,:,0)
        fhi0(:,:,1,1) = floc(:,:,1)
        call cfrure (nread, size(floc), floc, 'fhi')
        fhi = 0.0_R8
        fhi(:,:,0,0) = floc(:,:,0)
        fhi(:,:,1,1) = floc(:,:,1)
        if(b2fplasma_version.gt.'01.001.003') then
          call cfrure (nread, size(fhit), fhit, 'fhit')
        else
          fhit = fhi
        endif
        call cfrure (nread, size(flocs), flocs, 'fna')
        fna = 0.0_R8
        fna(:,:,0,0,:) = flocs(:,:,0,:)
        fna(:,:,1,1,:) = flocs(:,:,1,:)
        call cfrure (nread, size(flocs), flocs, 'fna0')
        fna0 = 0.0_R8
        fna0(:,:,0,0,:) = flocs(:,:,0,:)
        fna0(:,:,1,1,:) = flocs(:,:,1,:)
        call cfrure (nread, size(flocs), flocs, 'fnap')
        fnap = 0.0_R8
        fnap(:,:,0,0,:) = flocs(:,:,0,:)
        fnap(:,:,1,1,:) = flocs(:,:,1,:)
        call cfrure (nread, size(floc), floc, 'fne')
        fne = 0.0_R8
        fne(:,:,0,0) = floc(:,:,0)
        fne(:,:,1,1) = floc(:,:,1)
        call cfrure (nread, size(floc), floc, 'fni')
        fni = 0.0_R8
        fni(:,:,0,0) = floc(:,:,0)
        fni(:,:,1,1) = floc(:,:,1)
      endif
      call cfrure (nread, size(na), na, 'na')
      call cfrure (nread, size(na0), na0, 'na0')
      call cfrure (nread, size(nap), nap, 'nap')
      call cfrure (nread, size(ne), ne, 'ne')
      call cfrure (nread, size(ne0), ne0, 'ne0')
      call cfrure (nread, size(ne2), ne2, 'ne2')
      call cfrure (nread, size(nep), nep, 'nep')
      call cfrure (nread, size(ni), ni, 'ni')
      call cfrure (nread, size(ni0), ni0, 'ni0')
      call b2xpni (nxd, nyd, nsd, nap, nip)
      call cfrure (nread, size(pb), pb, 'pb')
      call cfrure (nread, size(po), po, 'po')
      call cfrure (nread, size(po0), po0, 'po0')
      call cfrure (nread, size(pop), pop, 'pop')
      call cfrure (nread, size(te), te, 'te')
      call cfrure (nread, size(te0), te0, 'te0')
      call cfrure (nread, size(tep), tep, 'tep')
      call cfrure (nread, size(ti), ti, 'ti')
      call cfrure (nread, size(ti0), ti0, 'ti0')
      call cfrure (nread, size(tip), tip, 'tip')
      call cfrure (nread, size(ua), ua, 'ua')
      call cfrure (nread, size(ua0), ua0, 'ua0')
      call cfrure (nread, size(uap), uap, 'uap')
      if (b2fplasma_version.ge.'03.001.000') then
        call cfrure (nread, size(uadia), uadia, 'uadia')
        call cfrure (nread, size(fchdia), fchdia, 'fchdia')
        call cfrure (nread, size(fmo),    fmo,    'fmo')
        call cfrure (nread, size(fna_32), fna_32, 'fna_32')
        call cfrure (nread, size(fna_52), fna_52, 'fna_52')
        call cfrure (nread, size(fni_32), fni_32, 'fni_32')
        call cfrure (nread, size(fni_52), fni_52, 'fni_52')
        call cfrure (nread, size(fne_32), fne_32, 'fne_32')
        call cfrure (nread, size(fne_52), fne_52, 'fne_52')
      else
        call cfrure (nread, size(flocs), uadia, 'uadia')
        uadia = 0.0_R8
        uadia(:,:,0,0,:) = flocs(:,:,0,:)
        uadia(:,:,1,1,:) = flocs(:,:,1,:)
        call cfrure (nread, size(floc), floc, 'fchdia')
        fchdia = 0.0_R8
        fchdia(:,:,0,0) = floc(:,:,0)
        fchdia(:,:,1,1) = floc(:,:,1)
        if(b2fplasma_version.ge.'01.001.008') then
          call cfrure (nread, size(flocs), flocs, 'fmo')
          fmo = 0.0_R8
          fmo(:,:,0,0,:) = flocs(:,:,0,:)
          fmo(:,:,1,1,:) = flocs(:,:,1,:)
          call cfrure (nread, size(flocs), flocs, 'fna_32')
          fna_32 = 0.0_R8
          fna_32(:,:,0,0,:) = flocs(:,:,0,:)
          fna_32(:,:,1,1,:) = flocs(:,:,1,:)
          call cfrure (nread, size(flocs), flocs, 'fna_52')
          fna_52 = 0.0_R8
          fna_52(:,:,0,0,:) = flocs(:,:,0,:)
          fna_52(:,:,1,1,:) = flocs(:,:,1,:)
          call cfrure (nread, size(floc), floc, 'fni_32')
          fni_32 = 0.0_R8
          fni_32(:,:,0,0) = floc(:,:,0)
          fni_32(:,:,1,1) = floc(:,:,1)
          call cfrure (nread, size(floc), floc, 'fni_52')
          fni_52 = 0.0_R8
          fni_52(:,:,0,0) = floc(:,:,0)
          fni_52(:,:,1,1) = floc(:,:,1)
          call cfrure (nread, size(floc), floc, 'fne_32')
          fne_32 = 0.0_R8
          fne_32(:,:,0,0) = floc(:,:,0)
          fne_32(:,:,1,1) = floc(:,:,1)
          call cfrure (nread, size(floc), floc, 'fne_52')
          fne_52 = 0.0_R8
          fne_52(:,:,0,0) = floc(:,:,0)
          fne_52(:,:,1,1) = floc(:,:,1)
        else
          do is = 0, nsd-1
            fmo(-1:nxd,-1:nyd,0,0,is) = am(is)*mp*ua(-1:nxd,-1:nyd,is)*
     .       fna(-1:nxd,-1:nyd,0,0,is)
            fmo(-1:nxd,-1:nyd,1,1,is) = am(is)*mp*ua(-1:nxd,-1:nyd,is)*
     .       fna(-1:nxd,-1:nyd,1,1,is)
          enddo
          fna_32(:,:,0,0,:) = fna(:,:,0,0,:)
          fna_32(:,:,1,1,:) = fna(:,:,1,1,:)
          fna_52 = 0.0_R8
          fni_32 = fni
          fni_52 = 0.0_R8
          fne_32 = fne
          fne_52 = 0.0_R8
        endif
      endif
      if(b2fplasma_version.ge.'03.001.000') then
        call cfrure (nread, size(wadia), wadia, 'wadia')
        call cfrure (nread, size(vaecrb), vaecrb, 'vaecrb')
      elseif(b2fplasma_version.ge.'01.001.010') then
        call cfrure (nread, size(flocs), wadia, 'wadia')
        wadia = 0.0_R8
        wadia(:,:,0,0,:) = flocs(:,:,0,:)
        wadia(:,:,1,1,:) = flocs(:,:,1,:)
        call cfrure (nread, size(flocs), vaecrb, 'vaecrb')
        vaecrb = 0.0_R8
        vaecrb(:,:,0,0,:) = flocs(:,:,0,:)
        vaecrb(:,:,1,1,:) = flocs(:,:,1,:)
      else
        wadia = uadia
        vaecrb = 0.0_R8
      endif
      if(b2fplasma_version.ge.'01.001.012') then
        call cfrure (nread, size(facdrift), facdrift, 'facdrift')
        call cfrure (nread, size(fac_ExB), fac_ExB, 'fac_ExB')
      else
        facdrift = 0.0_R8
        fac_ExB = 0.0_R8
      endif
      if(b2fplasma_version.ge.'03.001.000') then
        call cfrure (nread, size(fchvispar), fchvispar, 'fchvispar')
        call cfrure (nread, size(fchvisper), fchvisper, 'fchvisper')
        call cfrure (nread, size(fchin), fchin, 'fchin')
        call cfrure (nread, size(fna_nodrift), fna_nodrift,
     1   'fna_nodrift')
      elseif(b2fplasma_version.ge.'01.001.014') then
        call cfrure (nread, size(floc), floc, 'fchvispar')
        fchvispar = 0.0_R8
        fchvispar(:,:,0,0) = floc(:,:,0)
        fchvispar(:,:,1,1) = floc(:,:,1)
        call cfrure (nread, size(floc), floc, 'fchvisper')
        fchvisper = 0.0_R8
        fchvisper(:,:,0,0) = floc(:,:,0)
        fchvisper(:,:,1,1) = floc(:,:,1)
        call cfrure (nread, size(floc), floc, 'fchin')
        fchin = 0.0_R8
        fchin(:,:,0,0) = floc(:,:,0)
        fchin(:,:,1,1) = floc(:,:,1)
        call cfrure (nread, size(flocs), flocs,'fna_nodrift')
        fna_nodrift = 0.0_R8
        fna_nodrift(:,:,0,0,:) = flocs(:,:,0,:)
        fna_nodrift(:,:,1,1,:) = flocs(:,:,1,:)
      else
        fchvispar = 0.0_R8
        fchvisper = 0.0_R8
        fchin = 0.0_R8
        fna_nodrift(:,:,0,0,:) = fna(:,:,0,0,:)
        fna_nodrift(:,:,1,1,:) = fna(:,:,1,1,:)
      endif
      if(b2fplasma_version.ge.'01.001.015') then
        call cfrure (nread, size(fac_vis), fac_vis, 'fac_vis')
      else
        fac_vis = 0.0_R8
      endif
      if(b2fplasma_version.ge.'03.001.000') then
        call cfrure (nread, size(fna_mdf), fna_mdf, 'fna_mdf')
        call cfrure (nread, size(fhe_mdf), fhe_mdf, 'fhe_mdf')
        call cfrure (nread, size(fhi_mdf), fhi_mdf, 'fhi_mdf')
        call cfrure (nread, size(fnaPSch), fnaPSch, 'fnaPSch')
        call cfrure (nread, size(fhePSch), fhePSch, 'fhePSch')
        call cfrure (nread, size(fhiPSch), fhiPSch, 'fhiPSch')
        call cfrure (nread, size(fna_fcor), fna_fcor, 'fna_fcor')
        call cfrure (nread, size(fna_he), fna_he, 'fna_he')
        call cfrure (nread, size(fchvisq), fchvisq, 'fchvisq')
        call cfrure (nread, size(fchinert), fchinert, 'fchinert')
      elseif(b2fplasma_version.ge.'03.000.003') then
        call cfrure (nread, size(flocs), flocs, 'fna_mdf')
        fna_mdf = 0.0_R8
        fna_mdf(:,:,0,0,:) = flocs(:,:,0,:)
        fna_mdf(:,:,1,1,:) = flocs(:,:,1,:)
        call cfrure (nread, size(floc), floc, 'fhe_mdf')
        fhe_mdf = 0.0_R8
        fhe_mdf(:,:,0,0) = floc(:,:,0)
        fhe_mdf(:,:,1,1) = floc(:,:,1)
        call cfrure (nread, size(floc), floc, 'fhi_mdf')
        fhi_mdf = 0.0_R8
        fhi_mdf(:,:,0,0) = floc(:,:,0)
        fhi_mdf(:,:,1,1) = floc(:,:,1)
        call cfrure (nread, size(flocs), flocs, 'fnaPSch')
        fnaPSch = 0.0_R8
        fnaPSch(:,:,0,0,:) = flocs(:,:,0,:)
        fnaPSch(:,:,1,1,:) = flocs(:,:,1,:)
        call cfrure (nread, size(floc), floc, 'fhePSch')
        fhePSch = 0.0_R8
        fhePSch(:,:,0,0) = floc(:,:,0)
        fhePSch(:,:,1,1) = floc(:,:,1)
        call cfrure (nread, size(floc), floc, 'fhiPSch')
        fhiPSch = 0.0_R8
        fhiPSch(:,:,0,0) = floc(:,:,0)
        fhiPSch(:,:,1,1) = floc(:,:,1)
        call cfrure (nread, size(flocs), flocs, 'fna_fcor')
        fna_fcor = 0.0_R8
        fna_fcor(:,:,0,0,:) = flocs(:,:,0,:)
        fna_fcor(:,:,1,1,:) = flocs(:,:,1,:)
        call cfrure (nread, size(flocs), flocs, 'fna_he')
        fna_he = 0.0_R8
        fna_he(:,:,0,0,:) = flocs(:,:,0,:)
        fna_he(:,:,1,1,:) = flocs(:,:,1,:)
        call cfrure (nread, size(floc), floc, 'fchvisq')
        fchvisq = 0.0_R8
        fchvisq(:,:,0,0) = floc(:,:,0)
        fchvisq(:,:,1,1) = floc(:,:,1)
        call cfrure (nread, size(floc), floc, 'fchinert')
        fchinert = 0.0_R8
        fchinert(:,:,0,0) = floc(:,:,0)
        fchinert(:,:,1,1) = floc(:,:,1)
      endif
      if(b2fplasma_version.ge.'03.001.000') then
        call cfrure_opt (nread, size(floc), floc, 'fht')
        fht = 0.0_R8
        fht(:,:,0,0) = floc(:,:,0)
        fht(:,:,1,1) = floc(:,:,1)
        call cfrure_opt (nread, size(floc), floc, 'fhj')
        fhj = 0.0_R8
        fhj(:,:,0,0) = floc(:,:,0)
        fhj(:,:,1,1) = floc(:,:,1)
        call cfrure_opt (nread, size(flocs), flocs, 'fhm')
        fhm = 0.0_R8
        fhm(:,:,0,0,:) = flocs(:,:,0,:)
        fhm(:,:,1,1,:) = flocs(:,:,1,:)
        call cfrure_opt (nread, size(flocs), flocs, 'fhp')
        fhp = 0.0_R8
        fhp(:,:,0,0,:) = flocs(:,:,0,:)
        fhp(:,:,1,1,:) = flocs(:,:,1,:)
      elseif(b2fplasma_version.ge.'03.000.007') then
        call cfrure_opt (nread, size(fht), fht, 'fht')
        call cfrure_opt (nread, size(fhj), fhj, 'fhj')
        call cfrure_opt (nread, size(fhm), fhm, 'fhm')
        call cfrure_opt (nread, size(fhp), fhp, 'fhp')
      endif

      return
      end subroutine read_b2mod_plasma

      subroutine write_b2mod_plasma(nwrite)
      implicit none
      integer nwrite
      external cfwure

      call cfwure (nwrite, size(fch),  fch,  'fch')
      call cfwure (nwrite, size(fch0), fch0, 'fch0')
      call cfwure (nwrite, size(fchp), fchp, 'fchp')
      call cfwure (nwrite, size(fhe),  fhe,  'fhe')
      call cfwure (nwrite, size(fhe0), fhe0, 'fhe0')
      call cfwure (nwrite, size(fhep), fhep, 'fhep')
      call cfwure (nwrite, size(fhet), fhet, 'fhet')
      call cfwure (nwrite, size(fhi),  fhi,  'fhi')
      call cfwure (nwrite, size(fhi0), fhi0, 'fhi0')
      call cfwure (nwrite, size(fhip), fhip, 'fhip')
      call cfwure (nwrite, size(fhit), fhit, 'fhit')
      call cfwure (nwrite, size(fna),  fna,  'fna')
      call cfwure (nwrite, size(fna0), fna0, 'fna0')
      call cfwure (nwrite, size(fnap), fnap, 'fnap')
      call cfwure (nwrite, size(fne),  fne,  'fne')
      call cfwure (nwrite, size(fni),  fni,  'fni')
      call cfwure (nwrite, size(na),   na,   'na')
      call cfwure (nwrite, size(na0),  na0,  'na0')
      call cfwure (nwrite, size(nap),  nap,  'nap')
      call cfwure (nwrite, size(ne),   ne,   'ne')
      call cfwure (nwrite, size(ne0),  ne0,  'ne0')
      call cfwure (nwrite, size(ne2),  ne2,  'ne2')
      call cfwure (nwrite, size(nep),  nep,  'nep')
      call cfwure (nwrite, size(ni),   ni,   'ni')
      call cfwure (nwrite, size(ni0),  ni0,  'ni0')
      call cfwure (nwrite, size(pb),   pb,   'pb')
      call cfwure (nwrite, size(po),   po,   'po')
      call cfwure (nwrite, size(po0),  po0,  'po0')
      call cfwure (nwrite, size(pop),  pop,  'pop')
      call cfwure (nwrite, size(te),   te,   'te')
      call cfwure (nwrite, size(te0),  te0,  'te0')
      call cfwure (nwrite, size(tep),  tep,  'tep')
      call cfwure (nwrite, size(ti),   ti,   'ti')
      call cfwure (nwrite, size(ti0),  ti0,  'ti0')
      call cfwure (nwrite, size(tip),  tip,  'tip')
      call cfwure (nwrite, size(ua),   ua,   'ua')
      call cfwure (nwrite, size(ua0),  ua0,  'ua0')
      call cfwure (nwrite, size(uap),  uap,  'uap')
      call cfwure (nwrite, size(uadia), uadia, 'uadia')
      call cfwure (nwrite, size(fchdia), fchdia, 'fchdia')
      call cfwure (nwrite, size(fmo), fmo, 'fmo')
      call cfwure (nwrite, size(fna_32), fna_32, 'fna_32')
      call cfwure (nwrite, size(fna_52), fna_52, 'fna_52')
      call cfwure (nwrite, size(fni_32), fni_32, 'fni_32')
      call cfwure (nwrite, size(fni_52), fni_52, 'fni_52')
      call cfwure (nwrite, size(fne_32), fne_32, 'fne_32')
      call cfwure (nwrite, size(fne_52), fne_52, 'fne_52')
      call cfwure (nwrite, size(wadia),  wadia,  'wadia')
      call cfwure (nwrite, size(vaecrb), vaecrb, 'vaecrb')
      call cfwure (nwrite, size(facdrift), facdrift, 'facdrift')
      call cfwure (nwrite, size(fac_ExB), fac_ExB, 'fac_ExB')
      call cfwure (nwrite, size(fchvispar), fchvispar, 'fchvispar')
      call cfwure (nwrite, size(fchvisper), fchvisper, 'fchvisper')
      call cfwure (nwrite, size(fchin), fchin, 'fchin')
      call cfwure (nwrite, size(fna_nodrift), fna_nodrift,'fna_nodrift')
      call cfwure (nwrite, size(fac_vis), fac_vis, 'fac_vis')
      call cfwure (nwrite, size(fna_mdf), fna_mdf, 'fna_mdf')
      call cfwure (nwrite, size(fhe_mdf), fhe_mdf, 'fhe_mdf')
      call cfwure (nwrite, size(fhi_mdf), fhi_mdf, 'fhi_mdf')
      call cfwure (nwrite, size(fnaPSch), fnaPSch, 'fnaPSch')
      call cfwure (nwrite, size(fhePSch), fhePSch, 'fhePSch')
      call cfwure (nwrite, size(fhiPSch), fhiPSch, 'fhiPSch')
      call cfwure (nwrite, size(fna_fcor), fna_fcor, 'fna_fcor')
      call cfwure (nwrite, size(fna_he), fna_he, 'fna_he')
      call cfwure (nwrite, size(fchvisq), fchvisq, 'fchvisq')
      call cfwure (nwrite, size(fchinert), fchinert, 'fchinert')
      call cfwure (nwrite, size(fht), fht, 'fht')
      call cfwure (nwrite, size(fhj), fhj, 'fhj')
      call cfwure (nwrite, size(fhm), fhm, 'fhm')
      call cfwure (nwrite, size(fhp), fhp, 'fhp')

      return
      end subroutine write_b2mod_plasma
c
      subroutine set_facdrift_fac_ExB_init(nx, ny, ramp_slow)
      use b2mod_types
      use b2mod_indirect
      use b2mod_geo
      use b2mod_work
      implicit none
      integer nx, ny, ramp_slow
      real (kind=R8) ::
     & facdrift_tanh_a, facExB_tanh_a, facvis_tanh_a,
     & facdrift_tanh_b, facExB_tanh_b, facvis_tanh_b
      integer iy_nocoreExB, fac_ref
      integer ix, iy, jxi, jxa, jsep
      real (kind=R8) dssep
      external get_jsep, ipgeti, ipgetr, xertst
      save
     & iy_nocoreExB, fac_ref,
     & facdrift_tanh_a, facdrift_tanh_b,
     & facExB_tanh_a, facExB_tanh_b,
     & facvis_tanh_a, facvis_tanh_b
      data
     & iy_nocoreExB/-2/,
     & facdrift_tanh_a/0.0_R8/, facdrift_tanh_b/0.0_R8/,
     & facExB_tanh_a/0.0_R8/, facExB_tanh_b/0.0_R8/,
     & facvis_tanh_a/0.0_R8/, facvis_tanh_b/0.0_R8/

      if(facdrift_ncall.eq.0) then
       call get_jsep(nx,ny,jxi,jxa,jsep)
       fac_ref=jxa
       call ipgeti ('b2news_fac_ref',fac_ref)
       call xertst (-1.le.fac_ref.and.fac_ref.le.nx,
     2      'faulty internal parameter fac_ref')
       call ipgetr ('b2news_facdrift_start', facdrift_start)
       call ipgetr ('b2news_facdrift_target', facdrift_target)
       call ipgetr ('b2news_facdrift_inc', facdrift_inc)
       call ipgetr ('b2news_facdrift_dec', facdrift_dec)
       call xertst (0.0_R8.le.facdrift_start .and.                       !xpb
     &              1.0_R8.ge.facdrift_start,                            !xpb
     &   'faulty argument facdrift_start')                               !xpb
       call xertst (0.0_R8.le.facdrift_target .and.                      !xpb
     &              1.0_R8.ge.facdrift_target,                           !xpb
     &   'faulty argument facdrift_target')                              !xpb
       call xertst (1.0_R8.le.facdrift_inc,                              !xpb
     &   'faulty argument facdrift_inc')                                 !xpb
       call xertst (0.0_R8.le.facdrift_dec .and. 1.0_R8.ge.facdrift_dec,
     &   'faulty argument facdrift_dec')                                 !xpb
       facdrift_scalar=facdrift_start
       if(facdrift_start.eq.facdrift_target .and.
     &  facdrift_start.ne.0.0_R8) then
         write(*,*) 'facdrift set to a fixed value :', facdrift_start
       elseif(facdrift_start.ne.facdrift_target.and.ramp_slow.eq.0) then
         write(*,'(a,1p,1g9.2,a,1g9.2)')
     &    'facdrift ramped normally from ',facdrift_start,
     &    ' towards ',facdrift_target
       elseif(facdrift_start.ne.facdrift_target.and.ramp_slow.ne.0) then
         write(*,'(a,1p,1g9.2,a,1g9.2)')
     &    'facdrift ramped slowly from ',facdrift_start,
     &    ' towards ',facdrift_target
       else
         write(*,*) 'No diamagnetic drifts'
       endif
       call ipgetr ('b2news_facExB_start', facExB_start)
       call ipgetr ('b2news_facExB_target', facExB_target)
       call ipgetr ('b2news_facExB_inc', facExB_inc)
       call ipgetr ('b2news_facExB_dec', facExB_dec)
       call xertst (0.0_R8.le.facExB_start .and.
     &              1.0_R8.ge.facExB_start,
     &   'faulty argument facExB_start')
       call xertst (0.0_R8.le.facExB_target .and.
     &              1.0_R8.ge.facExB_target,
     &   'faulty argument facExB_target')
       call xertst (1.0_R8.le.facExB_inc,
     &   'faulty argument facExB_inc')
       call xertst (0.0_R8.le.facExB_dec .and. 1.0_R8.ge.facExB_dec,
     &   'faulty argument facExB_dec')
       call ipgetr ('b2news_ExB', fac_ExB_scalar)                        !xpb
       call xertst (0.0_R8.le.fac_ExB_scalar .and.
     &              1.0_R8.ge.fac_ExB_scalar,
     &  'faulty parameter b2news_ExB')                                   !xpb
       if(fac_ExB_scalar.ne.0.0_R8) then
c old style
         facExB_start=fac_ExB_scalar
         facExB_target=fac_ExB_scalar
         facExB_inc=1.0_R8
         facExB_dec=1.0_R8
         write(*,*) 'fac_ExB set to a fixed value: ', fac_ExB_scalar
       elseif(facExB_start.ne.0.0_R8.and.facExB_start.eq.facExB_target)
     >  then
         fac_ExB_scalar=facExB_start
         facExB_inc=1.0_R8
         facExB_dec=1.0_R8
         write(*,*) 'fac_ExB set to a fixed value: ', fac_ExB_scalar
       elseif(facExB_start.ne.0.0_R8.and.facExB_start.ne.facExB_target)
     >  then
         fac_ExB_scalar=facExB_start
         if (ramp_slow.eq.0) then
           write(*,'(a,1p,1g9.2,a,1g9.2)')
     &      'fac_ExB ramped normally from ',facExB_start,
     &      ' towards ',facExB_target
         else
           write(*,'(a,1p,1g9.2,a,1g9.2)')
     &      'fac_ExB ramped slowly from ',facExB_start,
     &      ' towards ',facExB_target
         endif
       else
         write(*,*) 'No ExB'
       endif
       call ipgeti ('b2news_iy_nocoreExB', iy_nocoreExB)                 !xpb
       call xertst (-2.le.iy_nocoreExB .and. ny.ge.iy_nocoreExB,         !xpb
     &   'faulty parameter iy_nocoreExB')                                !xpb
       do iy=-1,iy_nocoreExB
         do ix=-1,nx
           if(on_closed_surface(ix,iy)) then
             fac_ExB_profile(ix,iy)=0.0_R8
           endif
         enddo
       enddo
       call ipgetr ('b2news_facvis_start', facvis_start)
       call ipgetr ('b2news_facvis_target', facvis_target)
       call ipgetr ('b2news_facvis_inc', facvis_inc)
       call ipgetr ('b2news_facvis_dec', facvis_dec)
       call xertst (0.0_R8.le.facvis_start .and.
     &              1.0_R8.ge.facvis_start,
     &   'faulty argument facvis_start')
       call xertst (0.0_R8.le.facvis_target .and.
     &              1.0_R8.ge.facvis_target,
     &   'faulty argument facvis_target')
       call xertst (1.0_R8.le.facvis_inc,
     &   'faulty argument facvis_inc')
       call xertst (0.0_R8.le.facvis_dec .and. 1.0_R8.ge.facvis_dec,
     &   'faulty argument facvis_dec')
       call ipgetr ('b2news_vis', fac_vis_scalar)                        !xpb
       call xertst (0.0_R8.le.fac_vis_scalar .and.
     &              1.0_R8.ge.fac_vis_scalar,
     &  'faulty parameter b2news_vis')                                   !xpb
       if(fac_vis_scalar.ne.0.0_R8) then
c old style
         facvis_start=fac_vis_scalar
         facvis_target=fac_vis_scalar
         facvis_inc=1.0_R8
         facvis_dec=1.0_R8
         write(*,*) 'fac_vis set to a fixed value: ', fac_vis_scalar
       elseif(facvis_start.ne.0.0_R8.and.facvis_start.eq.facvis_target)
     >  then
         fac_vis_scalar=facvis_start
         facvis_inc=1.0_R8
         facvis_dec=1.0_R8
         write(*,*) 'fac_vis set to a fixed value: ', fac_vis_scalar
       elseif(facvis_start.ne.0.0_R8.and.facvis_start.ne.facvis_target)
     >  then
         fac_vis_scalar=facvis_start
         if (ramp_slow.eq.0) then
           write(*,'(a,1p,1g9.2,a,1g9.2)')
     &      'fac_vis ramped normally from ',facvis_start,
     &      ' towards ',facvis_target
         else
           write(*,'(a,1p,1g9.2,a,1g9.2)')
     &      'fac_vis ramped slowly from ',facvis_start,
     &      ' towards ',facvis_target
         endif
       else
         write(*,*) 'No vis'
       endif
       ds(-1)=
     1   sqrt((cr(fac_ref,-1)-
     1      0.5_R8*(crx(fac_ref,-1,0)+crx(fac_ref,-1,1)))**2+
     2      (cz(fac_ref,-1)-
     3      0.5_R8*(cry(fac_ref,-1,0)+cry(fac_ref,-1,1)))**2)
       do iy=0,ny
          ds(iy)=ds(iy-1)+
     1        sqrt((cr(fac_ref,iy)-cr(fac_ref,iy-1))**2+
     2        (cz(fac_ref,iy)-cz(fac_ref,iy-1))**2)
       enddo
       dssep=(ds(jsep)+ds(jsep+1))/2.0_R8
       do iy=-1,ny
         ds(iy)=ds(iy)-dssep
       enddo

       call ipgetr ('b2news_facdrift_tanh_a', facdrift_tanh_a)
       call ipgetr ('b2news_facdrift_tanh_b', facdrift_tanh_b)
       if(facdrift_tanh_b.ne.0.0_R8) then
          do iy=-1,ny
             do ix=-1,nx
                if(on_closed_surface(ix,iy)) then
                   facdrift_profile(ix,iy)=facdrift_profile(ix,iy)*
     &                  0.5_R8*(1.0_R8+
     &                    tanh((ds(iy)-facdrift_tanh_a)/
     &                   facdrift_tanh_b))
                endif
             enddo
          enddo
       endif
       call ipgetr ('b2news_fac_ExB_tanh_a', facExB_tanh_a)
       if (facExB_tanh_a.eq.0.0_R8)
     >  call ipgetr ('b2news_facExB_tanh_a', facExB_tanh_a)
       call ipgetr ('b2news_fac_ExB_tanh_b', facExB_tanh_b)
       if (facExB_tanh_b.eq.0.0_R8)
     >  call ipgetr ('b2news_facExB_tanh_b', facExB_tanh_b)
       if(facExB_tanh_b.ne.0.0_R8) then
          do iy=-1,ny
             do ix=-1,nx
                if(on_closed_surface(ix,iy)) then
                   fac_ExB_profile(ix,iy)=fac_ExB_profile(ix,iy)*
     &                  0.5_R8*(1.0_R8+
     &                    tanh((ds(iy)-facExB_tanh_a)/
     &                   facExB_tanh_b))
                endif
             enddo
          enddo
       endif
       call ipgetr ('b2news_fac_vis_tanh_a', facvis_tanh_a)
       if (facvis_tanh_a.eq.0.0_R8)
     >  call ipgetr ('b2news_facvis_tanh_a', facvis_tanh_a)
       call ipgetr ('b2news_fac_vis_tanh_b', facvis_tanh_b)
       if (facvis_tanh_b.eq.0.0_R8)
     >  call ipgetr ('b2news_facvis_tanh_b', facvis_tanh_b)
       if(facvis_tanh_b.ne.0.0_R8) then
          do iy=-1,ny
             do ix=-1,nx
                if(on_closed_surface(ix,iy)) then
                   fac_vis_profile(ix,iy)=fac_vis_profile(ix,iy)*
     &                  0.5_R8*(1.0_R8+
     &                    tanh((ds(iy)-facvis_tanh_a)/
     &                   facvis_tanh_b))
                endif
             enddo
          enddo
       endif
      endif
      call set_facdrift_fac_ExB_scale
      return
      end subroutine set_facdrift_fac_ExB_init

      subroutine set_facdrift_fac_ExB_inc
      implicit none
      external xertst
      call xertst (0.lt.facdrift_ncall,
     &     'set_facdrift_fac_ExB_inc called before initialisation')
      facdrift_scalar=min(facdrift_scalar*facdrift_inc,facdrift_target)
      fac_ExB_scalar=min(fac_ExB_scalar*facExB_inc,facExB_target)
      fac_vis_scalar=min(fac_vis_scalar*facvis_inc,facvis_target)
      call set_facdrift_fac_ExB_scale
      return
      end subroutine set_facdrift_fac_ExB_inc

      subroutine set_facdrift_fac_ExB_dec
      implicit none
      external xertst
      call xertst (0.lt.facdrift_ncall,
     &     'set_facdrift_fac_ExB_dec called before initialisation')
      facdrift_scalar=max(facdrift_scalar*facdrift_dec,facdrift_start)
      fac_ExB_scalar=max(fac_ExB_scalar*facExB_dec,facExB_start)
      fac_vis_scalar=max(fac_vis_scalar*facvis_dec,facvis_start)
      call set_facdrift_fac_ExB_scale
      return
      end subroutine set_facdrift_fac_ExB_dec

      subroutine set_facdrift_fac_ExB_scale
      implicit none
      facdrift=facdrift_scalar*facdrift_profile
      fac_ExB=fac_ExB_scalar*fac_ExB_profile
      fac_vis=fac_vis_scalar*fac_vis_profile
      write(*,'(1x,a,1p,6g12.5)') 'facdrift, fac_ExB, fac_vis',
     &   minval(facdrift),maxval(facdrift),
     &   minval(fac_ExB),maxval(fac_ExB),
     &   minval(fac_vis),maxval(fac_vis)
      facdrift_ncall=facdrift_ncall+1
      return
      end subroutine set_facdrift_fac_ExB_scale

      function rfchx(nx,ny,ix)
      implicit none
      integer, intent(in) :: nx, ny
      integer, intent(in) :: ix
      real(kind=R8) :: rfchx
      real(kind=R8) :: b2ssum
      external b2ssum
      rfchx = b2ssum(ny+2,fch(ix,-1,0,0),nx+2)
      return
      end function rfchx

      function rfchy(nx,iy)
      implicit none
      integer, intent(in) :: nx
      integer, intent(in) :: iy
      real(kind=R8) :: rfchy
      real(kind=R8) :: b2ssum
      external b2ssum
      rfchy = b2ssum(nx+2,fch(-1,iy,1,1),1)
      return
      end function rfchy

      function rfhex(nx,ny,ix)
      implicit none
      integer, intent(in) :: nx, ny
      integer, intent(in) :: ix
      real(kind=R8) :: rfhex
      real(kind=R8) :: b2ssum
      external b2ssum
      rfhex = b2ssum(ny+2,fhe(ix,-1,0,0),nx+2)
      return
      end function rfhex

      function rfhey(nx,iy)
      implicit none
      integer, intent(in) :: nx
      integer, intent(in) :: iy
      real(kind=R8) :: rfhey
      real(kind=R8) :: b2ssum
      external b2ssum
      rfhey = b2ssum(nx+2,fhe(-1,iy,1,1),1)
      return
      end function rfhey

      function rfhix(nx,ny,ix)
      implicit none
      integer, intent(in) :: nx, ny
      integer, intent(in) :: ix
      real(kind=R8) :: rfhix
      real(kind=R8) :: b2ssum
      external b2ssum
      rfhix = b2ssum(ny+2,fhi(ix,-1,0,0),nx+2)
      return
      end function rfhix

      function rfhiy(nx,iy)
      implicit none
      integer, intent(in) :: nx
      integer, intent(in) :: iy
      real(kind=R8) :: rfhiy
      real(kind=R8) :: b2ssum
      external b2ssum
      rfhiy = b2ssum(nx+2,fhi(-1,iy,1,1),1)
      return
      end function rfhiy

      function rfnax(nx,ny,ix,is)
      implicit none
      integer, intent(in) :: nx, ny
      integer, intent(in) :: ix, is
      real(kind=R8) :: rfnax
      real(kind=R8) :: b2ssum
      external b2ssum
      rfnax = b2ssum(ny+2,fna(ix,-1,0,0,is),nx+2)
      return
      end function rfnax

      function rfnay(nx,iy,is)
      implicit none
      integer, intent(in) :: nx
      integer, intent(in) :: iy, is
      real(kind=R8) :: rfnay
      real(kind=R8) :: b2ssum
      external b2ssum
      rfnay = b2ssum(nx+2,fna(-1,iy,1,1,is),1)
      return
      end function rfnay

      function rfnex(nx,ny,ix)
      implicit none
      integer, intent(in) :: nx, ny
      integer, intent(in) :: ix
      real(kind=R8) :: rfnex
      real(kind=R8) :: b2ssum
      external b2ssum
      rfnex = b2ssum(ny+2,fne(ix,-1,0,0),nx+2)
      return
      end function rfnex

      function rfney(nx,iy)
      implicit none
      integer, intent(in) :: nx
      integer, intent(in) :: iy
      real(kind=R8) :: rfney
      real(kind=R8) :: b2ssum
      external b2ssum
      rfney = b2ssum(nx+2,fne(-1,iy,1,1),1)
      return
      end function rfney

      function rfnix(nx,ny,ix)
      implicit none
      integer, intent(in) :: nx, ny
      integer, intent(in) :: ix
      real(kind=R8) :: rfnix
      real(kind=R8) :: b2ssum
      external b2ssum
      rfnix = b2ssum(ny+2,fni(ix,-1,0,0),nx+2)
      return
      end function rfnix

      function rfniy(nx,iy)
      implicit none
      integer, intent(in) :: nx
      integer, intent(in) :: iy
      real(kind=R8) :: rfniy
      real(kind=R8) :: b2ssum
      external b2ssum
      rfniy = b2ssum(nx+2,fni(-1,iy,1,1),1)
      return
      end function rfniy

      function roxa(ix,iy,is)
      use b2mod_constants
      use b2mod_b2cmpa
      implicit none
      integer, intent(in) :: ix, iy, is
      real(kind=R8) :: roxa
      roxa = am(is)*mp*na(ix,iy,is)
      return
      end function roxa

      function uasx(ix,iy,is)
      use b2mod_geo
      use b2mod_indirect
      implicit none
      integer, intent(in) :: ix, iy, is
      real(kind=R8) :: uasx
      uasx =
     . (hx(ix,iy)*ua(leftix(ix,iy),leftiy(ix,iy),is)+
     1  hx(leftix(ix,iy),leftiy(ix,iy))*ua(ix,iy,is))/
     2 (hx(leftix(ix,iy),leftiy(ix,iy))+hx(ix,iy))*pbs(ix,iy,0)
      return
      end function uasx

      function uasy(ix,iy,is)
      use b2mod_geo
      use b2mod_indirect
      implicit none
      integer, intent(in) :: ix, iy, is
      real(kind=R8) :: uasy
      uasy =
     . (hy1(ix,iy)*ua(bottomix(ix,iy),bottomiy(ix,iy),is)+
     1  hy1(bottomix(ix,iy),bottomiy(ix,iy))*ua(ix,iy,is))/
     2 (hy1(bottomix(ix,iy),bottomiy(ix,iy))+hy1(ix,iy))*pbs(ix,iy,1)
      return
      end function uasy

      function uax(ix,iy,is)
      use b2mod_geo
      use b2mod_indirect
      implicit none
      integer, intent(in) :: ix, iy, is
      real(kind=R8) :: uax
      uax = fna(ix,iy,0,0,is)*
     &  (hx(ix,iy)/na(leftix(ix,iy),leftiy(ix,iy),is)+
     1   hx(leftix(ix,iy),leftiy(ix,iy))/na(ix,iy,is))/
     &  (vol(leftix(ix,iy),leftiy(ix,iy))+vol(ix,iy))
      return
      end function uax

      function uay(ix,iy,is)
      use b2mod_geo
      use b2mod_indirect
      implicit none
      integer, intent(in) :: ix, iy, is
      real(kind=R8) :: uay
      uay = fna(ix,iy,1,1,is)*
     &  (hy1(ix,iy)/na(bottomix(ix,iy),bottomiy(ix,iy),is)+
     1   hy1(bottomix(ix,iy),bottomiy(ix,iy))/na(ix,iy,is))/
     &  (vol(bottomix(ix,iy),bottomiy(ix,iy))+vol(ix,iy))
      return
      end function uay

      function uesx(ix,iy)
      use b2mod_geo
      use b2mod_indirect
      implicit none
      integer, intent(in) :: ix, iy
      real(kind=R8) :: uesx
      uesx =
     . (hx(ix,iy)*ue(leftix(ix,iy),leftiy(ix,iy))+
     1  hx(leftix(ix,iy),leftiy(ix,iy))*ue(ix,iy))/
     2 (hx(leftix(ix,iy),leftiy(ix,iy))+hx(ix,iy))*pbs(ix,iy,0)
      return
      end function uesx

      end module b2mod_plasma
