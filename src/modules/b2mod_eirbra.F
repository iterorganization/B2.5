      module b2mod_eirbra
      use b2mod_types
      implicit none
c
C     EIRENE TO BRAAMS
c
      real (kind=R8), allocatable ::
     ,  sni(:,:,:), smo(:,:,:), see(:,:), sei(:,:,:),
     ,  sna_paat(:,:,:), sna_pmat(:,:,:),
     ,  sna_piat(:,:,:), sna_paml(:,:,:),
     ,  sna_pmml(:,:,:), sna_piml(:,:,:),
     ,  sna_paio(:,:,:), sna_pmio(:,:,:),
     ,  sna_piio(:,:,:),
     ,  sni_papl(:,:,:), sni_pmpl(:,:,:),
     ,  sni_pipl(:,:,:), sni_pppl(:,:,:),
     ,  smo_mapl(:,:,:), smo_mmpl(:,:,:),
     ,  smo_mipl(:,:,:), smo_cppv(:,:,:),
     ,  sne_pael(:,:), sne_pmel(:,:), sne_piel(:,:),
     ,  see_eael(:,:), see_emel(:,:),
     ,  see_eiel(:,:), see_epel(:,:),
     ,  sei_eapl(:,:,:), sei_empl(:,:,:),
     ,  sei_eipl(:,:,:), sei_eppl(:,:,:),
     ,  fluxh(:,:),
     ,  srcstrn(:),
     ,  volsumn(:), volsumm(:), volsumei(:), volsumee(:),
     ,  flxspci(:,:), srccrfc(:,:)
c
c*** fluxh  : not used
c*** volsumn: plasma particle source due to volumetric recombination
c*** volsum[m,ee,ei]: the same for momentum and energy
c*** srcstrn: intensities of different neutral sources (fluxt)
c*** flxspci: fluxes of plasma ions to the recycling surfaces
c*** srccrfc: source correction factors (from infcop)
C

C> Default output file prefix
#ifndef HPUX
      CHARACTER(5), PUBLIC, PARAMETER :: FORT    = 'FORT.'
      CHARACTER(5), PUBLIC, PARAMETER :: FORT_LC = 'fort.'
#else
      CHARACTER(3), PUBLIC, PARAMETER :: FORT    = 'FTN'
      CHARACTER(3), PUBLIC, PARAMETER :: FORT_LC = 'ftn'
#endif

      contains

      subroutine alloc_b2mod_eirbra(ncv,nsd,natmd,nmold,niond,nplsd,
     .  nstrad)
      implicit none
      integer ncv, nsd, natmd, nmold, niond, nplsd, nstrad

      allocate(sni(ncv,nsd,nstrad))
      allocate(smo(ncv,nsd,nstrad))
      allocate(see(ncv,nstrad))
      allocate(sei(ncv,nsd,nstrad))
      allocate(sna_paat(ncv,max(1,natmd),nstrad))
      allocate(sna_pmat(ncv,max(1,natmd),nstrad))
      allocate(sna_piat(ncv,max(1,natmd),nstrad))
      allocate(sna_paml(ncv,max(1,nmold),nstrad))
      allocate(sna_pmml(ncv,max(1,nmold),nstrad))
      allocate(sna_piml(ncv,max(1,nmold),nstrad))
      allocate(sna_paio(ncv,max(1,niond),nstrad))
      allocate(sna_pmio(ncv,max(1,niond),nstrad))
      allocate(sna_piio(ncv,max(1,niond),nstrad))
      allocate(sni_papl(ncv,nsd,nstrad))
      allocate(sni_pmpl(ncv,nsd,nstrad))
      allocate(sni_pipl(ncv,nsd,nstrad))
      allocate(sni_pppl(ncv,nsd,nstrad))
      allocate(smo_mapl(ncv,nsd,nstrad))
      allocate(smo_mmpl(ncv,nsd,nstrad))
      allocate(smo_mipl(ncv,nsd,nstrad))
      allocate(smo_cppv(ncv,nsd,nstrad))
      allocate(sne_pael(ncv,nstrad))
      allocate(sne_pmel(ncv,nstrad))
      allocate(sne_piel(ncv,nstrad))
      allocate(see_eael(ncv,nstrad))
      allocate(see_emel(ncv,nstrad))
      allocate(see_eiel(ncv,nstrad))
      allocate(see_epel(ncv,nstrad))
      allocate(sei_eapl(ncv,nsd,nstrad))
      allocate(sei_empl(ncv,nsd,nstrad))
      allocate(sei_eipl(ncv,nsd,nstrad))
      allocate(sei_eppl(ncv,nsd,nstrad))
      allocate(fluxh(10,nstrad))
      allocate(srcstrn(nstrad))
      allocate(volsumn(nstrad))
      allocate(volsumm(nstrad))
      allocate(volsumei(nstrad))
      allocate(volsumee(nstrad))
      allocate(flxspci(nplsd,nstrad))
      allocate(srccrfc(max(1,natmd),nstrad))

      call init_b2mod_eirbra

      return
      end subroutine alloc_b2mod_eirbra

      subroutine init_b2mod_eirbra
      implicit none

      volsumn  = 0.0_R8
      volsumm  = 0.0_R8
      volsumei = 0.0_R8
      volsumee = 0.0_R8

      return
      end subroutine init_b2mod_eirbra

      subroutine dealloc_b2mod_eirbra
      implicit none

      deallocate(sni)
      deallocate(smo)
      deallocate(see)
      deallocate(sei)
      deallocate(sna_paat)
      deallocate(sna_pmat)
      deallocate(sna_piat)
      deallocate(sna_paml)
      deallocate(sna_pmml)
      deallocate(sna_piml)
      deallocate(sna_paio)
      deallocate(sna_pmio)
      deallocate(sna_piio)
      deallocate(sni_papl)
      deallocate(sni_pmpl)
      deallocate(sni_pipl)
      deallocate(sni_pppl)
      deallocate(smo_mapl)
      deallocate(smo_mmpl)
      deallocate(smo_mipl)
      deallocate(smo_cppv)
      deallocate(sne_pael)
      deallocate(sne_pmel)
      deallocate(sne_piel)
      deallocate(see_eael)
      deallocate(see_emel)
      deallocate(see_eiel)
      deallocate(see_epel)
      deallocate(sei_eapl)
      deallocate(sei_empl)
      deallocate(sei_eipl)
      deallocate(sei_eppl)
      deallocate(fluxh)
      deallocate(srcstrn)
      deallocate(volsumn)
      deallocate(volsumm)
      deallocate(volsumei)
      deallocate(volsumee)
      deallocate(flxspci)
      deallocate(srccrfc)

      return
      end subroutine dealloc_b2mod_eirbra

      end module b2mod_eirbra
