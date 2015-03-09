      subroutine eirene_mc_init
      use b2mod_types
#ifdef B25_EIRENE
      use eirmod_parmmod
      use eirmod_extrab25
#endif
      implicit none
#ifndef B25_EIRENE
#     include <COUPLE/KOPPLDIM.F>
#     include <COUPLE/EIRDIAG.F>
#     include <COUPLE/EIRBRA.F>
#endif
cwdk
c#ifndef B25_EIRENE                                                   !srv (mail 15.01.08)
c      call eirene_extrab25_eirdiag_alloc(nxdd,nydd,
c     .                     natm,nmol,nion,nlmpgs,nstra,npls)          !srv 14.10.13
c      call eirene_extrab25_eirpbls_alloc(1,natm,nmol,nion,npls,       !srv 14.10.13
c     &                                       nstra,nfl)               !srv 03.04.13
c#endif

cwdk  copy-paste from eirene_extrab25_eirdiag_alloc (Eirene_dummy)
c     avoids double module definition in case code is compiled with Eirene 
#ifdef B25_EIRENE  
      if(.not.allocated(dab2)) then
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
        allocate(edissml(0:nxdd,0:nydd,nmol,1))
        allocate(wldnek(nlmpgs,0:nstra+1)) 
        allocate(wldnep(nlmpgs,0:nstra+1))
        allocate(wldna(nlmpgs,natm,0:nstra+1))
        allocate(ewlda(nlmpgs,natm,0:nstra+1))
        allocate(wldnm(nlmpgs,nmol,0:nstra+1))
        allocate(ewldm(nlmpgs,nmol,0:nstra+1))
        allocate(wldra(nlmpgs,natm,0:nstra+1))
        allocate(wldrm(nlmpgs,nmol,0:nstra+1))
        allocate(wldpp(nlmpgs,npls,0:nstra+1))
        allocate(wldpa(nlmpgs,natm,0:nstra+1))
        allocate(wldpm(nlmpgs,nmol,0:nstra+1))
        allocate(wldpeb(nlmpgs,0:nstra+1))
        allocate(wldspt(nlmpgs,0:nstra+1))
        allocate(eneutrad(0:nxdd,0:nydd,natm,1))
        allocate(isrftype(nlmpgs))
#endif
!       write(*,*) ' Allocated arrays in eirene_extrab25_eirdiag_alloc ',  &
!        ' with nxdd,nydd,natm,nmol,nion,nlmpgs,nstra,npls',             &
!               nxdd,nydd,natm,nmol,nion,nlmpgs,nstra,npls
        dab2=0
        dmb2=0
        dib2=0
        tab2=0
        tmb2=0
        tib2=0
        rfluxa=0
        rfluxm=0
        refluxa=0
        refluxm=0
        pfluxa=0
        pfluxm=0
        pefluxa=0
        pefluxm=0
        emiss=0
        emissmol=0
        srcml=0
        edissml=0
        wldnek=0
        wldnep=0
        wldna=0
        ewlda=0
        wldnm=0
        ewldm=0
        wldra=0
        wldrm=0
        wldpp=0
        wldpa=0
        wldpm=0
        wldpeb=0
        wldspt=0
        eneutrad=0
        isrftype=0
#ifdef B25_EIRENE
      endif
#endif

cwdk  copy-paste from eirene_extrab25_eirpbls_alloc (Eirene_dummy, only case im=1)
c     avoids double module definition in case code is compiled with Eirene  
c     Note: all arrays allocated in eirene_extrab25_eirpbls_init (Eirene) are allocated here,
c     so the standard cleanup routine can still be used (checks on allocation of lkindp)
#ifdef B25_EIRENE
      if(.not.allocated(lkindm)) then 
        allocate(lkindm(nmol))
        allocate(lkindp(npls))
        allocate(lkindi(nion))
        allocate(srccrfc(natm,nstra))
        allocate(srcstrn(nstra))
#endif
        lkindm=0
        lkindp=0
        lkindi=0
        srccrfc=0.
        srcstrn=0.
#ifdef B25_EIRENE
      endif 
#endif

      return
      end
