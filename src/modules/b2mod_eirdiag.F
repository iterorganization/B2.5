      module b2mod_eirdiag
      use b2mod_types
      implicit none

      real (kind=R8), allocatable ::
     &  dab2(:,:,:), dmb2(:,:,:), dib2(:,:,:),
     &  tab2(:,:,:), tmb2(:,:,:), tib2(:,:,:),
     &  rfluxa(:,:,:), rfluxm(:,:,:),
     &  refluxa(:,:,:), refluxm(:,:,:),
     &  pfluxa(:,:,:), pfluxm(:,:,:),
     &  pefluxa(:,:,:), pefluxm(:,:,:),
     &  tfluxa(:,:,:), tfluxm(:,:,:),
     &  tefluxa(:,:,:), tefluxm(:,:,:),
     &  emiss(:,:,:), emissmol(:,:,:),
     &  srcml(:,:,:), edissml(:,:,:),
     &  wldnek(:,:), wldnep(:,:),
     &  wldna(:,:,:), ewlda(:,:,:),
     &  wldnm(:,:,:), ewldm(:,:,:),
     &  wldra(:,:,:), wldrm(:,:,:),
     &  wldpp(:,:,:), wldpa(:,:,:),
     &  wldpm(:,:,:), wldpeb(:,:),
     &  wldspt(:,:), wldspta(:,:,:), wldsptm(:,:,:),
     &  eneutrad(:,:,:),
     &  wlarea(:), wlabsrp(:,:),
     &  wldna_res(:,:),     !incident atom flux
     &  wldnm_res(:,:),     !incident molecule flux
     &  ewlda_res(:,:),     !incident energy flux of atoms
     &  ewldm_res(:,:),     !incident energy flux of molecules
     &  ewldea_res(:,:),    !emitted energy flux of atoms
     &  ewldem_res(:,:),    !net incident energy flux due of molecules
     &  ewldrp_res(:),      !kinetic energy of reflected neutrals originated from ions
     &  ewldmr_res(:,:),    !energy due to recombination of atoms into molecules
     &  wldspt_res(:),      !flux of sputtered wall material
     &  wldspta_res(:,:),   !sputtered flux for each type of emitted atom
     &  wldsptm_res(:,:),   !sputtered flux for each type of emitted molecule
     &  wlpump(:,:),        !pumped flux for each type of pumped particle
     &  wlpump_res(:,:),    !pumped flux for each type of pumped particle
     &  sarea_res(:)        !area of the surface elements

CVK NEUTRAL FLUXES, SPATIALLY RESOLVED ON NON-DEFAULT STANDARD SURFACES (NDS)
C First dimension: index of atom or molecule. Second dimension is the index of
C surface elements, controlled by array 'eirdiag_nds_ind'
C

! See WNEUTRAL_FLUXES for the description
      INTEGER, allocatable :: eirdiag_nds_ind(:), fnmti_b2(:)

C  eirdiag_nds_ind(IS) : is the index AFTER WHICH the data for surface IS is started in arrays '*_res'.
C                        IF _ind(IS)<0 then this surface is skipped
C  fnmti_b2            : the list of the triangle faces making up the non-standard surfaces
C
C
      CHARACTER(8), allocatable :: eirtxt(:)

      integer nnlimi,nnstsi,nnatmi,nnmoli,nnioni,nnplsi
      integer, allocatable :: isrftype(:)
      integer, allocatable :: lkindp(:), lkindm(:), lkindi(:)
      logical lhalpha,lvib
c
c*** Volume data:
c***    srcml   :   power loss due to molecules, including
c***    edissml :   power loss due to molecule dissociation
c*** Surface data:
c***    wldnek  :   net kinetic energy deposited by neutrals
c***    wldnep  :   energy released by recombination into molecules
c***    wldna   :   flux of atoms impinging onto the surface
c***    ewlda   :   their average energy
c***    wldnm   :   flux of molecules impinging onto the surface
c***    ewldm   :   their average energy
c***    wldra   :   flux of reflected atoms
c***    wldrm   :   flux of reflected molecules
c***    wldpp   :   flux of plasma ions impinging onto the surface
c***    wldpa   :   flux of resulting atoms
c***    wldpm   :   flux of resulting molecules
c***    wldpeb  :   power carried away by these atoms and molecules
c***    wldspt  :   flux of sputtered wall material
c***    isrftype:   surface type (iliin in Eirene)
c***    eneutrad:   power radiated due to neutrals (only atoms at the moment)
c

      contains

      subroutine alloc_b2mod_eirdiag(nCv,
     ,  natmd,nmold,niond,nplsd,nlimps,nstrad,nstsd,nsgmx)
      implicit none
      integer nCv, natmd, nmold, niond, nplsd,
     &        nlimps, nstrad, nstsd, nsgmx

      allocate(dab2(nCv,max(1,natmd),1))
      allocate(dmb2(nCv,max(1,nmold),1))
      allocate(dib2(nCv,max(1,niond),1))
      allocate(tab2(nCv,max(1,natmd),1))
      allocate(tmb2(nCv,max(1,nmold),1))
      allocate(tib2(nCv,max(1,niond),1))
      allocate(rfluxa(nCv,max(1,natmd),1))
      allocate(rfluxm(nCv,max(1,nmold),1))
      allocate(refluxa(nCv,max(1,natmd),1))
      allocate(refluxm(nCv,max(1,nmold),1))
      allocate(pfluxa(nCv,max(1,natmd),1))
      allocate(pfluxm(nCv,max(1,nmold),1))
      allocate(pefluxa(nCv,max(1,natmd),1))
      allocate(pefluxm(nCv,max(1,nmold),1))
      allocate(tfluxa(nCv,max(1,natmd),1))
      allocate(tfluxm(nCv,max(1,nmold),1))
      allocate(tefluxa(nCv,max(1,natmd),1))
      allocate(tefluxm(nCv,max(1,nmold),1))
      allocate(emiss(nCv,1,1))
      allocate(emissmol(nCv,1,1))
      allocate(srcml(nCv,max(1,nmold),1))
      allocate(edissml(nCv,max(1,nmold),0:nstrad))
      allocate(wldnek(nlimps,0:nstrad))
      allocate(wldnep(nlimps,0:nstrad))
      allocate(wldna(nlimps,max(1,natmd),0:nstrad))
      allocate(ewlda(nlimps,max(1,natmd),0:nstrad))
      allocate(wldnm(nlimps,max(1,nmold),0:nstrad))
      allocate(ewldm(nlimps,max(1,nmold),0:nstrad))
      allocate(wldra(nlimps,max(1,natmd),0:nstrad))
      allocate(wldrm(nlimps,max(1,nmold),0:nstrad))
      allocate(wldpp(nlimps,nplsd,0:nstrad))
      allocate(wldpa(nlimps,max(1,natmd),0:nstrad))
      allocate(wldpm(nlimps,max(1,nmold),0:nstrad))
      allocate(wldpeb(nlimps,0:nstrad))
      allocate(wldspt(nlimps,0:nstrad))
      allocate(wldspta(nlimps,max(1,natmd),0:nstrad))
      allocate(wldsptm(nlimps,max(1,nmold),0:nstrad))
      allocate(eneutrad(nCv,max(1,natmd),0:nstrad))
      allocate(wlarea(nlimps))
      allocate(wlabsrp(natmd+nmold+niond+nplsd,nlimps))
      allocate(wldna_res(max(1,natmd),nsgmx))
      allocate(wldnm_res(max(1,nmold),nsgmx))
      allocate(ewlda_res(max(1,natmd),nsgmx))
      allocate(ewldm_res(max(1,nmold),nsgmx))
      allocate(ewldea_res(max(1,natmd),nsgmx))
      allocate(ewldem_res(max(1,nmold),nsgmx))
      allocate(ewldrp_res(nsgmx))
      allocate(ewldmr_res(max(1,nmold),nsgmx))
      allocate(wldspt_res(nsgmx))
      allocate(wldspta_res(max(1,natmd),nsgmx))
      allocate(wldsptm_res(max(1,nmold),nsgmx))
      allocate(wlpump(natmd+nmold+niond+nplsd,nsgmx))
      allocate(wlpump_res(natmd+nmold+niond+nplsd,nsgmx))
      allocate(sarea_res(nsgmx))
      allocate(eirdiag_nds_ind(nstsd+1))
      allocate(fnmti_b2(nstsd+1))
      allocate(isrftype(nlimps))
      allocate(lkindp(nplsd))
      allocate(lkindm(max(1,nmold)))
      allocate(lkindi(max(1,niond)))

      return
      end subroutine alloc_b2mod_eirdiag

      subroutine dealloc_b2mod_eirdiag
      implicit none

      deallocate(dab2)
      deallocate(dmb2)
      deallocate(dib2)
      deallocate(tab2)
      deallocate(tmb2)
      deallocate(tib2)
      deallocate(rfluxa)
      deallocate(rfluxm)
      deallocate(refluxa)
      deallocate(refluxm)
      deallocate(pfluxa)
      deallocate(pfluxm)
      deallocate(pefluxa)
      deallocate(pefluxm)
      deallocate(tfluxa)
      deallocate(tfluxm)
      deallocate(tefluxa)
      deallocate(tefluxm)
      deallocate(emiss)
      deallocate(emissmol)
      deallocate(srcml)
      deallocate(edissml)
      deallocate(wldnek)
      deallocate(wldnep)
      deallocate(wldna)
      deallocate(ewlda)
      deallocate(wldnm)
      deallocate(ewldm)
      deallocate(wldra)
      deallocate(wldrm)
      deallocate(wldpp)
      deallocate(wldpa)
      deallocate(wldpm)
      deallocate(wldpeb)
      deallocate(wldspt)
      deallocate(wldspta)
      deallocate(wldsptm)
      deallocate(eneutrad)
      deallocate(wlarea)
      deallocate(wlabsrp)
      deallocate(wldna_res)
      deallocate(wldnm_res)
      deallocate(ewlda_res)
      deallocate(ewldm_res)
      deallocate(ewldea_res)
      deallocate(ewldem_res)
      deallocate(ewldrp_res)
      deallocate(ewldmr_res)
      deallocate(wldspt_res)
      deallocate(wldspta_res)
      deallocate(wldsptm_res)
      deallocate(wlpump)
      deallocate(wlpump_res)
      deallocate(sarea_res)
      deallocate(eirdiag_nds_ind)
      deallocate(fnmti_b2)
      deallocate(isrftype)
      deallocate(lkindp)
      deallocate(lkindm)
      deallocate(lkindi)

      return
      end subroutine dealloc_b2mod_eirdiag

      end module b2mod_eirdiag
