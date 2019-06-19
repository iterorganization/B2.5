      module b2mod_eirbra
      use b2mod_types
      implicit none
c
C     EIRENE TO BRAAMS
c
      real (kind=R8), allocatable ::
     ,  sni(:,:,:,:), smo(:,:,:,:), see(:,:,:), sei(:,:,:),
     ,  sni_papl(:,:,:,:), sni_pmpl(:,:,:,:),
     ,  sni_pipl(:,:,:,:), sni_pppl(:,:,:,:),
     ,  smo_mapl(:,:,:,:), smo_mmpl(:,:,:,:),
     ,  smo_mipl(:,:,:,:), smo_cppv(:,:,:,:),
     ,  see_eael(:,:,:), see_emel(:,:,:),
     ,  see_eiel(:,:,:), see_epel(:,:,:),
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

      contains

      subroutine alloc_b2mod_eirbra(nxd,nyd,nsd,natmd,nplsd,nstrad)
      implicit none
      integer nxd, nyd, nsd, natmd, nplsd, nstrad

      allocate(sni(0:nxd+1,0:nyd+1,nsd,nstrad))
      allocate(smo(0:nxd+1,0:nyd+1,nsd,nstrad))
      allocate(see(0:nxd+1,0:nyd+1,nstrad))
      allocate(sei(0:nxd+1,0:nyd+1,nstrad))
      allocate(sni_papl(0:nxd+1,0:nyd+1,nsd,nstrad))
      allocate(sni_pmpl(0:nxd+1,0:nyd+1,nsd,nstrad))
      allocate(sni_pipl(0:nxd+1,0:nyd+1,nsd,nstrad))
      allocate(sni_pppl(0:nxd+1,0:nyd+1,nsd,nstrad))
      allocate(smo_mapl(0:nxd+1,0:nyd+1,nsd,nstrad))
      allocate(smo_mmpl(0:nxd+1,0:nyd+1,nsd,nstrad))
      allocate(smo_mipl(0:nxd+1,0:nyd+1,nsd,nstrad))
      allocate(smo_cppv(0:nxd+1,0:nyd+1,nsd,nstrad))
      allocate(see_eael(0:nxd+1,0:nyd+1,nstrad))
      allocate(see_emel(0:nxd+1,0:nyd+1,nstrad))
      allocate(see_eiel(0:nxd+1,0:nyd+1,nstrad))
      allocate(see_epel(0:nxd+1,0:nyd+1,nstrad))
      allocate(sei_eapl(0:nxd+1,0:nyd+1,nstrad))
      allocate(sei_empl(0:nxd+1,0:nyd+1,nstrad))
      allocate(sei_eipl(0:nxd+1,0:nyd+1,nstrad))
      allocate(sei_eppl(0:nxd+1,0:nyd+1,nstrad))
      allocate(fluxh(10,nstrad))
      allocate(srcstrn(nstrad))
      allocate(volsumn(nstrad))
      allocate(volsumm(nstrad))
      allocate(volsumei(nstrad))
      allocate(volsumee(nstrad))
      allocate(flxspci(nplsd,nstrad))
      allocate(srccrfc(natmd,nstrad))

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
      deallocate(sni_papl)
      deallocate(sni_pmpl)
      deallocate(sni_pipl)
      deallocate(sni_pppl)
      deallocate(smo_mapl)
      deallocate(smo_mmpl)
      deallocate(smo_mipl)
      deallocate(smo_cppv)
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
