      subroutine eirene_mc_init(nx,ny)
      use b2mod_types
#ifdef B25_EIRENE
      use eirmod_parmmod
      use eirmod_eirbra
      use eirmod_comusr
     , , only : lkindm, lkindp, lkindi
      use eirmod_wneutrals
      use eirmod_extrab25
#else
      use b2mod_eirbra
      use b2mod_eirdiag
#endif
      implicit none
      integer nx, ny
#ifndef B25_EIRENE
#     include <COUPLE/KOPPLDIM.F>
#endif

*   ..test nx, ny
      call xertst (0.le.nx.and.0.le.ny, 'faulty argument nx, ny')
cwdk  copy-paste from eirene_extrab25 allocations (Eirene_dummy)
c     avoids double module definition in case code is compiled with Eirene 
      if(.not.allocated(dab2)) then
#ifdef B25_EIRENE
        allocate(dab2(0:nxdd,0:nydd,natm,1))
        allocate(dmb2(0:nxdd,0:nydd,nmol,1))
        allocate(dib2(0:nxdd,0:nydd,nion,1))
        allocate(tab2(0:nxdd,0:nydd,natm,1))
        allocate(tmb2(0:nxdd,0:nydd,nmol,1))
        allocate(tib2(0:nxdd,0:nydd,nion,1))
        allocate(rfluxa(0:nxdd,0:nydd,natm,1))
        allocate(rfluxm(0:nxdd,0:nydd,nmol,1))
        allocate(refluxa(0:nxdd,0:nydd,natm,1))
        allocate(refluxm(0:nxdd,0:nydd,nmol,1))
        allocate(pfluxa(0:nxdd,0:nydd,natm,1))
        allocate(pfluxm(0:nxdd,0:nydd,nmol,1))
        allocate(pefluxa(0:nxdd,0:nydd,natm,1))
        allocate(pefluxm(0:nxdd,0:nydd,nmol,1))
        allocate(emiss(0:nxdd,0:nydd,1,1))
        allocate(emissmol(0:nxdd,0:nydd,1,1))
        allocate(srcml(0:nxdd,0:nydd,nmol,1))
        allocate(edissml(0:nxdd,0:nydd,nmol,0:nstra+1))
        allocate(wldnek(nlimps,0:nstra+1)) 
        allocate(wldnep(nlimps,0:nstra+1))
        allocate(wldna(nlimps,natm,0:nstra+1))
        allocate(ewlda(nlimps,natm,0:nstra+1))
        allocate(wldnm(nlimps,nmol,0:nstra+1))
        allocate(ewldm(nlimps,nmol,0:nstra+1))
        allocate(wldra(nlimps,natm,0:nstra+1))
        allocate(wldrm(nlimps,nmol,0:nstra+1))
        allocate(wldpp(nlimps,npls,0:nstra+1))
        allocate(wldpa(nlimps,natm,0:nstra+1))
        allocate(wldpm(nlimps,nmol,0:nstra+1))
        allocate(wldpeb(nlimps,0:nstra+1))
        allocate(wldspt(nlimps,0:nstra+1))
        allocate(eneutrad(0:nxdd,0:nydd,natm,0:nstra+1))
        allocate(emolrad(0:nxdd,0:nydd,nmol,0:nstra+1))
        allocate(eionrad(0:nxdd,0:nydd,nion,0:nstra+1))
        allocate(isrftype(nlimps))
#else
        call alloc_b2mod_eirdiag(nx,ny,
     .   natm,nmol,nion,npls,nlimps,nstra,nsts,nsgmx)
#endif
!       write(*,*) ' Allocated arrays in eirene_extrab25_eirdiag_alloc ',  &
!        ' with nxdd,nydd,natm,nmol,nion,nlimps,nstra,npls',             &
!               nxdd,nydd,natm,nmol,nion,nlimps,nstra,npls
        dab2=0.0_R8
        dmb2=0.0_R8
        dib2=0.0_R8
        tab2=0.0_R8
        tmb2=0.0_R8
        tib2=0.0_R8
        rfluxa=0.0_R8
        rfluxm=0.0_R8
        refluxa=0.0_R8
        refluxm=0.0_R8
        pfluxa=0.0_R8
        pfluxm=0.0_R8
        pefluxa=0.0_R8
        pefluxm=0.0_R8
        emiss=0.0_R8
        emissmol=0.0_R8
        srcml=0.0_R8
        edissml=0.0_R8
        eneutrad=0.0_R8
        emolrad=0.0_R8
        eionrad=0.0_R8
        wldnek=0.0_R8
        wldnep=0.0_R8
        wldna=0.0_R8
        ewlda=0.0_R8
        wldnm=0.0_R8
        ewldm=0.0_R8
        wldra=0.0_R8
        wldrm=0.0_R8
        wldpp=0.0_R8
        wldpa=0.0_R8
        wldpm=0.0_R8
        wldpeb=0.0_R8
        wldspt=0.0_R8
        isrftype=0
      endif

cwdk  copy-paste from eirene_extrab25_eirpbls allocations (Eirene_dummy, only case im=1)
c     avoids double module definition in case code is compiled with Eirene  
c     Note: all arrays allocated in eirene_extrab25_eirpbls_init (Eirene) are allocated here,
c     so the standard cleanup routine can still be used (checks on allocation of lkindp)
      if(.not.allocated(lkindm)) then 
#ifdef B25_EIRENE
        if (.not.allocated(lkindm)) allocate(lkindm(nmol))
        if (.not.allocated(lkindp)) allocate(lkindp(npls))
        if (.not.allocated(lkindi)) allocate(lkindi(nion))
        if (.not.allocated(srccrfc)) allocate(srccrfc(natm,nstra))
        if (.not.allocated(srcstrn)) allocate(srcstrn(nstra))
#endif
        lkindm=0
        lkindp=0
        lkindi=0
        srccrfc=0.0_R8
        srcstrn=0.0_R8
        flxspci=0.0_R8
      endif 

      return
      end subroutine eirene_mc_init
