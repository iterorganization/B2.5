      module b2mod_sources
      use b2mod_types
      use b2mod_version
      implicit none
      real (kind=R8), allocatable, save ::
     *  sch(:,:,:),
     *  she(:,:,:),
     *  shi(:,:,:),
     *  sne(:,:,:),
     *  smo(:,:,:,:),
     *  shef(:,:,:),                                           !srv 10.11.02                
     *  shif(:,:,:),                                           !srv 10.11.02 
     *  smof(:,:,:,:),                                         !srv 10.11.02 
     *  smq(:,:,:,:),
     *  sna(:,:,:,:),
     *  rsana(:,:,:),
     *  rsahi(:,:,:),
     *  rsamo(:,:,:),
     *  rrana(:,:,:),
     *  rrahi(:,:,:),
     *  rramo(:,:,:),
     *  rqahe(:,:,:),
     *  rqrad(:,:,:),
     *  rqbrm(:,:,:),
     *  rcxna(:,:,:),
     *  rcxhi(:,:,:),
     *  rcxmo(:,:,:)
      real (kind=R8), allocatable, save ::
     *  b2stbr_sna(:,:,:),
     *  b2stbr_smo(:,:,:),
     *  b2stbr_sne(:,:),
     *  b2stbr_she(:,:),
     *  b2stbr_shi(:,:),
     *  b2stbr_sch(:,:),
     *  b2stbc_sna(:,:,:),
     *  b2stbc_smo(:,:,:),
     *  b2stbc_sne(:,:),
     *  b2stbc_she(:,:),
     *  b2stbc_shi(:,:),
     *  b2stbc_sch(:,:),
     *  b2stbm_sna(:,:,:),
     *  b2stbm_smo(:,:,:),
     *  b2stbm_sne(:,:),
     *  b2stbm_she(:,:),
     *  b2stbm_shi(:,:),
     *  b2stbm_sch(:,:),
     *  b2sihs_divue(:,:),
     *  b2sihs_divua(:,:),
     *  b2sihs_exbe(:,:),
     *  b2sihs_exba(:,:),
     *  b2sihs_visa(:,:),
     *  b2sihs_joule(:,:),
     *  b2sihs_fraa(:,:),
     *  b2sihs_str(:,:),
     1  b2npmo_smaf(:,:,:,:),
     1  b2npmo_smag(:,:,:,:),
     1  b2npmo_smav(:,:,:,:),
     1  smpr(:,:,:),
     1  smpt(:,:,:),
     1  smfr(:,:,:),
     1  smcf(:,:,:),
     *  ext_sna(:,:,:),
     *  ext_smo(:,:,:),
     *  ext_sne(:,:),
     *  ext_she(:,:),
     *  ext_shi(:,:),
     *  ext_sch(:,:)
cdjm Jan2017
      real (kind=R8), allocatable, save ::
     *  b2stbr_she0to3(:,:,:),
     *  b2stbr_shi0to3(:,:,:),
     *  b2stbc_she0to3(:,:,:),
     *  b2stbc_shi0to3(:,:,:),
     *  b2stel_sna_ion(:,:,:),
     *  b2stel_sna_rec(:,:,:),
     *  b2stel_smq_ion(:,:,:),
     *  b2stel_smq_rec(:,:,:),
     *  b2stel_she(:,:),
     *  b2stel_shi_ion(:,:),
     *  b2stel_shi_rec(:,:),
     *  b2stel_she0to3(:,:,:),
     *  b2stel_shi_ion0to3(:,:,:),
     *  b2stel_shi_rec0to3(:,:,:),
     *  b2stcx_sna(:,:,:),
     *  b2stcx_smq(:,:,:),
     *  b2stcx_shi(:,:),
     *  b2stcx_shi0to3(:,:,:),
     *  b2srsm_sna(:,:,:),
     *  b2srsm_smo(:,:,:),
     *  b2srsm_she(:,:),
     *  b2srsm_shi(:,:),
     *  b2srsm_she0to3(:,:,:),
     *  b2srsm_shi0to3(:,:,:),
     *  b2srsm_sch(:,:),
     *  b2srdt_sna(:,:,:),
     *  b2srdt_smo(:,:,:),
     *  b2srdt_she(:,:),
     *  b2srdt_shi(:,:),
     *  b2srdt_she0to3(:,:,:),
     *  b2srdt_shi0to3(:,:,:),
     *  b2srdt_sch(:,:),
     *  b2srst_sna(:,:,:),
     *  b2srst_smo(:,:,:),
     *  b2srst_she(:,:),
     *  b2srst_shi(:,:),
     *  b2srst_she0to3(:,:,:),
     *  b2srst_shi0to3(:,:,:),
     *  b2srst_sch(:,:),
     *  eirene_mc_t0_papl(:,:,:,:),
     *  eirene_mc_t0_pmpl(:,:,:,:),
     *  eirene_mc_t0_pipl(:,:,:,:),
     *  eirene_mc_t0_pppl(:,:,:,:),
     *  eirene_mc_t0_core(:,:,:),
     *  eirene_mc_t0_mapl(:,:,:,:),
     *  eirene_mc_t0_mmpl(:,:,:,:),
     *  eirene_mc_t0_mipl(:,:,:,:),
     *  eirene_mc_t0_cppv(:,:,:,:),
     *  eirene_mc_t0_eael(:,:,:),
     *  eirene_mc_t0_emel(:,:,:),
     *  eirene_mc_t0_eiel(:,:,:),
     *  eirene_mc_t0_epel(:,:,:),
     *  eirene_mc_t0_eapl(:,:,:),
     *  eirene_mc_t0_empl(:,:,:),
     *  eirene_mc_t0_eipl(:,:,:),
     *  eirene_mc_t0_eppl(:,:,:),
     *  b2sihs_diae(:,:),
     *  b2sihs_diaa(:,:),
     *  b2sifr_smoch(:,:,:),
     *  b2sifr_smotf(:,:,:),
     *  b2sigp_smogp(:,:,:),
     *  b2sigp_pstat(:,:,:),
     *  b2siav_smovh(:,:,:),
     *  b2siav_smovv(:,:,:),
     *  b2nxdv_smo(:,:,:),
     *  b2npht_shei(:,:)
      real (kind=R8), allocatable, save ::                             !srv 11.09.09 16.07.10
     *  snasm(:,:,:,:), smosm(:,:,:,:), 
     &  shesm(:,:,:), shism(:,:,:),
     *  snadt(:,:,:,:), smodt(:,:,:,:), 
     &  shedt(:,:,:), shidt(:,:,:),
     &  schdt(:,:,:), snedt(:,:,:)

      contains

      subroutine alloc_b2mod_sources(nxd,nyd,nsd)
      implicit none
      integer nxd,nyd,nsd

      allocate(sch(-1:nxd,-1:nyd,0:3))
      allocate(she(-1:nxd,-1:nyd,0:3))
      allocate(sne(-1:nxd,-1:nyd,0:1))
      allocate(shi(-1:nxd,-1:nyd,0:3))
      allocate(smo(-1:nxd,-1:nyd,0:3,0:nsd-1))
      allocate(sna(-1:nxd,-1:nyd,0:1,0:nsd-1))
      allocate(shef(-1:nxd,-1:nyd,0:3))                        !sv 10.11.02
      allocate(shif(-1:nxd,-1:nyd,0:3))                        !sv 10.11.02
      allocate(smof(-1:nxd,-1:nyd,0:3,0:nsd-1))                !sv 10.11.02
      allocate(smq(-1:nxd,-1:nyd,0:3,0:nsd-1))
      allocate(b2npmo_smag(-1:nxd,-1:nyd,0:3,0:nsd-1))
      allocate(b2npmo_smaf(-1:nxd,-1:nyd,0:3,0:nsd-1))
      allocate(b2npmo_smav(-1:nxd,-1:nyd,0:3,0:nsd-1))
      allocate(smpr(-1:nxd,-1:nyd,0:nsd-1))
      allocate(smpt(-1:nxd,-1:nyd,0:nsd-1))
      allocate(smfr(-1:nxd,-1:nyd,0:nsd-1))
      allocate(smcf(-1:nxd,-1:nyd,0:nsd-1))
      allocate(rsana(-1:nxd,-1:nyd,0:nsd-1))
      allocate(rsahi(-1:nxd,-1:nyd,0:nsd-1))
      allocate(rsamo(-1:nxd,-1:nyd,0:nsd-1))
      allocate(rrana(-1:nxd,-1:nyd,0:nsd-1))
      allocate(rrahi(-1:nxd,-1:nyd,0:nsd-1))
      allocate(rramo(-1:nxd,-1:nyd,0:nsd-1))
      allocate(rqahe(-1:nxd,-1:nyd,0:nsd-1))
      allocate(rqrad(-1:nxd,-1:nyd,0:nsd-1))
      allocate(rqbrm(-1:nxd,-1:nyd,0:nsd-1))
      allocate(rcxna(-1:nxd,-1:nyd,0:nsd-1))
      allocate(rcxhi(-1:nxd,-1:nyd,0:nsd-1))
      allocate(rcxmo(-1:nxd,-1:nyd,0:nsd-1))
      allocate(b2stbr_sna(-1:nxd,-1:nyd,0:nsd-1))
      allocate(b2stbr_smo(-1:nxd,-1:nyd,0:nsd-1))
      allocate(b2stbr_sne(-1:nxd,-1:nyd))
      allocate(b2stbr_she(-1:nxd,-1:nyd))
      allocate(b2stbr_shi(-1:nxd,-1:nyd))
      allocate(b2stbr_sch(-1:nxd,-1:nyd))
      allocate(b2stbc_sna(-1:nxd,-1:nyd,0:nsd-1))
      allocate(b2stbc_smo(-1:nxd,-1:nyd,0:nsd-1))
      allocate(b2stbc_sne(-1:nxd,-1:nyd))
      allocate(b2stbc_she(-1:nxd,-1:nyd))
      allocate(b2stbc_shi(-1:nxd,-1:nyd))
      allocate(b2stbc_sch(-1:nxd,-1:nyd))
      allocate(b2stbm_sna(-1:nxd,-1:nyd,0:nsd-1))
      allocate(b2stbm_smo(-1:nxd,-1:nyd,0:nsd-1))
      allocate(b2stbm_sne(-1:nxd,-1:nyd))
      allocate(b2stbm_she(-1:nxd,-1:nyd))
      allocate(b2stbm_shi(-1:nxd,-1:nyd))
      allocate(b2stbm_sch(-1:nxd,-1:nyd))
      allocate(b2sihs_divue(-1:nxd,-1:nyd))
      allocate(b2sihs_divua(-1:nxd,-1:nyd))
      allocate(b2sihs_exbe(-1:nxd,-1:nyd))
      allocate(b2sihs_exba(-1:nxd,-1:nyd))
      allocate(b2sihs_visa(-1:nxd,-1:nyd))
      allocate(b2sihs_joule(-1:nxd,-1:nyd))
      allocate(b2sihs_fraa(-1:nxd,-1:nyd))
      allocate(b2sihs_str(-1:nxd,-1:nyd))
      allocate(ext_sna(-1:nxd,-1:nyd,0:nsd-1))
      allocate(ext_smo(-1:nxd,-1:nyd,0:nsd-1))
      allocate(ext_sne(-1:nxd,-1:nyd))
      allocate(ext_she(-1:nxd,-1:nyd))
      allocate(ext_shi(-1:nxd,-1:nyd))
      allocate(ext_sch(-1:nxd,-1:nyd))
      allocate(snedt(-1:nxd,-1:nyd,0:1))
      allocate(schdt(-1:nxd,-1:nyd,0:3))                        !srv 16.07.10 {
      allocate(shedt(-1:nxd,-1:nyd,0:3))
      allocate(shidt(-1:nxd,-1:nyd,0:3))
      allocate(smodt(-1:nxd,-1:nyd,0:3,0:nsd-1))
      allocate(snadt(-1:nxd,-1:nyd,0:1,0:nsd-1))
      allocate(shesm(-1:nxd,-1:nyd,0:3))
      allocate(shism(-1:nxd,-1:nyd,0:3))
      allocate(smosm(-1:nxd,-1:nyd,0:3,0:nsd-1))
      allocate(snasm(-1:nxd,-1:nyd,0:1,0:nsd-1))                 !srv 16.07.10 }
cdjm Jan2017
      allocate(b2sifr_smoch(-1:nxd,-1:nyd,0:nsd-1))
      allocate(b2sifr_smotf(-1:nxd,-1:nyd,0:nsd-1))
      allocate(b2sigp_smogp(-1:nxd,-1:nyd,0:nsd-1))
      allocate(b2sigp_pstat(-1:nxd,-1:nyd,0:nsd-1))
      allocate(b2siav_smovh(-1:nxd,-1:nyd,0:nsd-1))
      allocate(b2siav_smovv(-1:nxd,-1:nyd,0:nsd-1))
      allocate(b2nxdv_smo(-1:nxd,-1:nyd,0:nsd-1))
      allocate(b2stbr_she0to3(-1:nxd,-1:nyd,0:3))
      allocate(b2stbr_shi0to3(-1:nxd,-1:nyd,0:3))
      allocate(b2stbc_she0to3(-1:nxd,-1:nyd,0:3))
      allocate(b2stbc_shi0to3(-1:nxd,-1:nyd,0:3))
      allocate(b2stel_sna_ion(-1:nxd,-1:nyd,0:nsd-1))
      allocate(b2stel_sna_rec(-1:nxd,-1:nyd,0:nsd-1))
      allocate(b2stel_smq_ion(-1:nxd,-1:nyd,0:nsd-1))
      allocate(b2stel_smq_rec(-1:nxd,-1:nyd,0:nsd-1))
      allocate(b2stel_she(-1:nxd,-1:nyd))
      allocate(b2stel_shi_ion(-1:nxd,-1:nyd))
      allocate(b2stel_shi_rec(-1:nxd,-1:nyd))
      allocate(b2stel_she0to3(-1:nxd,-1:nyd,0:3))
      allocate(b2stel_shi_ion0to3(-1:nxd,-1:nyd,0:3))
      allocate(b2stel_shi_rec0to3(-1:nxd,-1:nyd,0:3))
      allocate(b2stcx_sna(-1:nxd,-1:nyd,0:nsd-1))
      allocate(b2stcx_smq(-1:nxd,-1:nyd,0:nsd-1))
      allocate(b2stcx_shi(-1:nxd,-1:nyd))
      allocate(b2stcx_shi0to3(-1:nxd,-1:nyd,0:3))
      allocate(b2srsm_sna(-1:nxd,-1:nyd,0:nsd-1))
      allocate(b2srsm_smo(-1:nxd,-1:nyd,0:nsd-1))
      allocate(b2srsm_she(-1:nxd,-1:nyd))
      allocate(b2srsm_shi(-1:nxd,-1:nyd))
      allocate(b2srsm_she0to3(-1:nxd,-1:nyd,0:3))
      allocate(b2srsm_shi0to3(-1:nxd,-1:nyd,0:3))
      allocate(b2srsm_sch(-1:nxd,-1:nyd))
      allocate(b2srdt_sna(-1:nxd,-1:nyd,0:nsd-1))
      allocate(b2srdt_smo(-1:nxd,-1:nyd,0:nsd-1))
      allocate(b2srdt_she(-1:nxd,-1:nyd))
      allocate(b2srdt_shi(-1:nxd,-1:nyd))
      allocate(b2srdt_she0to3(-1:nxd,-1:nyd,0:3))
      allocate(b2srdt_shi0to3(-1:nxd,-1:nyd,0:3))
      allocate(b2srdt_sch(-1:nxd,-1:nyd))     
      allocate(b2srst_sna(-1:nxd,-1:nyd,0:nsd-1))
      allocate(b2srst_smo(-1:nxd,-1:nyd,0:nsd-1))
      allocate(b2srst_she(-1:nxd,-1:nyd))
      allocate(b2srst_shi(-1:nxd,-1:nyd))
      allocate(b2srst_she0to3(-1:nxd,-1:nyd,0:3))
      allocate(b2srst_shi0to3(-1:nxd,-1:nyd,0:3))
      allocate(b2srst_sch(-1:nxd,-1:nyd))     
      allocate(b2sihs_diae(-1:nxd,-1:nyd))
      allocate(b2sihs_diaa(-1:nxd,-1:nyd))
      allocate(b2npht_shei(-1:nxd,-1:nyd))

      rsana=0.0_R8
      rsahi=0.0_R8
      rsamo=0.0_R8
      rrana=0.0_R8
      rrahi=0.0_R8
      rramo=0.0_R8
      rqahe=0.0_R8
      rqrad=0.0_R8
      rqbrm=0.0_R8
      rcxna=0.0_R8
      rcxhi=0.0_R8
      rcxmo=0.0_R8
      b2stbr_sna=0.0_R8
      b2stbr_smo=0.0_R8
      b2stbr_sne=0.0_R8
      b2stbr_she=0.0_R8
      b2stbr_shi=0.0_R8
      b2stbr_sch=0.0_R8
      b2stbc_sna=0.0_R8
      b2stbc_smo=0.0_R8
      b2stbc_sne=0.0_R8
      b2stbc_she=0.0_R8
      b2stbc_shi=0.0_R8
      b2stbc_sch=0.0_R8
      b2stbm_sna=0.0_R8
      b2stbm_smo=0.0_R8
      b2stbm_sne=0.0_R8
      b2stbm_she=0.0_R8
      b2stbm_shi=0.0_R8
      b2stbm_sch=0.0_R8
      b2sihs_divue=0.0_R8
      b2sihs_divua=0.0_R8
      b2sihs_exbe=0.0_R8
      b2sihs_exba=0.0_R8
      b2sihs_visa=0.0_R8
      b2sihs_joule=0.0_R8
      b2sihs_fraa=0.0_R8
      b2sihs_str=0.0_R8
      b2npmo_smaf=0.0_R8
      b2npmo_smag=0.0_R8
      b2npmo_smav=0.0_R8
      smpr=0.0_R8
      smpt=0.0_R8
      smfr=0.0_R8
      smcf=0.0_R8
      ext_sna=0.0_R8
      ext_smo=0.0_R8
      ext_sne=0.0_R8
      ext_she=0.0_R8
      ext_shi=0.0_R8
      ext_sch=0.0_R8
cdjm Jan2017
      b2stbr_she0to3=0.0_R8
      b2stbr_shi0to3=0.0_R8
      b2stbc_she0to3=0.0_R8
      b2stbc_shi0to3=0.0_R8
      b2stel_sna_ion=0.0_R8
      b2stel_sna_rec=0.0_R8
      b2stel_smq_ion=0.0_R8
      b2stel_smq_rec=0.0_R8
      b2stel_she=0.0_R8
      b2stel_shi_ion=0.0_R8
      b2stel_shi_rec=0.0_R8
      b2stel_she0to3=0.0_R8
      b2stel_shi_ion0to3=0.0_R8
      b2stel_shi_rec0to3=0.0_R8
      b2stcx_sna=0.0_R8
      b2stcx_smq=0.0_R8
      b2stcx_shi=0.0_R8
      b2stcx_shi0to3=0.0_R8
      b2srsm_sna=0.0_R8
      b2srsm_smo=0.0_R8
      b2srsm_she=0.0_R8
      b2srsm_shi=0.0_R8
      b2srsm_she0to3=0.0_R8
      b2srsm_shi0to3=0.0_R8
      b2srsm_sch=0.0_R8
      b2srdt_sna=0.0_R8
      b2srdt_smo=0.0_R8
      b2srdt_she=0.0_R8
      b2srdt_shi=0.0_R8
      b2srdt_she0to3=0.0_R8
      b2srdt_shi0to3=0.0_R8
      b2srdt_sch=0.0_R8
      b2srst_sna=0.0_R8
      b2srst_smo=0.0_R8
      b2srst_she=0.0_R8
      b2srst_shi=0.0_R8
      b2srst_she0to3=0.0_R8
      b2srst_shi0to3=0.0_R8
      b2srst_sch=0.0_R8
      b2sihs_diae=0.0_R8
      b2sihs_diaa=0.0_R8
      b2sifr_smoch=0.0_R8
      b2sifr_smotf=0.0_R8
      b2sigp_smogp=0.0_R8
      b2sigp_pstat=0.0_R8
      b2siav_smovh=0.0_R8
      b2siav_smovv=0.0_R8
      b2nxdv_smo=0.0_R8
      b2npht_shei=0.0_R8

      return
      end subroutine alloc_b2mod_sources

      subroutine dealloc_b2mod_sources
      implicit none

      deallocate(sch)
      deallocate(she)
      deallocate(sne)
      deallocate(shi)
      deallocate(smo)
      deallocate(shef)                        !sv 10.11.02
      deallocate(shif)                        !sv 10.11.02
      deallocate(smof)                        !sv 10.11.02
      deallocate(smq)
      deallocate(b2npmo_smag)
      deallocate(b2npmo_smaf)
      deallocate(b2npmo_smav)
      deallocate(smpr)
      deallocate(smpt)
      deallocate(smfr)
      deallocate(smcf)
      deallocate(sna)
      deallocate(rsana)
      deallocate(rsahi)
      deallocate(rsamo)
      deallocate(rrana)
      deallocate(rrahi)
      deallocate(rramo)
      deallocate(rqahe)
      deallocate(rqrad)
      deallocate(rqbrm)
      deallocate(rcxna)
      deallocate(rcxhi)
      deallocate(rcxmo)
      deallocate(b2stbr_sna)
      deallocate(b2stbr_smo)
      deallocate(b2stbr_sne)
      deallocate(b2stbr_she)
      deallocate(b2stbr_shi)
      deallocate(b2stbr_sch)
      deallocate(b2stbc_sna)
      deallocate(b2stbc_smo)
      deallocate(b2stbc_sne)
      deallocate(b2stbc_she)
      deallocate(b2stbc_shi)
      deallocate(b2stbc_sch)
      deallocate(b2stbm_sna)
      deallocate(b2stbm_smo)
      deallocate(b2stbm_sne)
      deallocate(b2stbm_she)
      deallocate(b2stbm_shi)
      deallocate(b2stbm_sch)
      deallocate(b2sihs_divue)
      deallocate(b2sihs_divua)
      deallocate(b2sihs_exbe)
      deallocate(b2sihs_exba)
      deallocate(b2sihs_visa)
      deallocate(b2sihs_joule)
      deallocate(b2sihs_fraa)
      deallocate(b2sihs_str)
      deallocate(ext_sna)
      deallocate(ext_smo)
      deallocate(ext_sne)
      deallocate(ext_she)
      deallocate(ext_shi)
      deallocate(ext_sch)
      deallocate(schdt)                        !srv 16.07.10 {
      deallocate(snedt)
      deallocate(shedt)
      deallocate(shidt)
      deallocate(smodt)
      deallocate(snadt)
      deallocate(shesm)
      deallocate(shism)
      deallocate(smosm)
      deallocate(snasm)                         !srv 16.07.10 }
cdjm Jan2017
      deallocate(b2sifr_smoch)
      deallocate(b2sifr_smotf)
      deallocate(b2sigp_smogp)
      deallocate(b2sigp_pstat)
      deallocate(b2siav_smovh)
      deallocate(b2siav_smovv)
      deallocate(b2nxdv_smo)
      deallocate(b2stbr_she0to3)
      deallocate(b2stbr_shi0to3)
      deallocate(b2stbc_she0to3)
      deallocate(b2stbc_shi0to3)
      deallocate(b2stel_sna_ion)
      deallocate(b2stel_sna_rec)
      deallocate(b2stel_smq_ion)
      deallocate(b2stel_smq_rec)
      deallocate(b2stel_she)
      deallocate(b2stel_shi_ion)
      deallocate(b2stel_shi_rec)
      deallocate(b2stel_she0to3)
      deallocate(b2stel_shi_ion0to3)
      deallocate(b2stel_shi_rec0to3)
      deallocate(b2stcx_sna)
      deallocate(b2stcx_smq)
      deallocate(b2stcx_shi)
      deallocate(b2stcx_shi0to3)
      deallocate(b2srsm_sna)
      deallocate(b2srsm_smo)
      deallocate(b2srsm_she)
      deallocate(b2srsm_shi)
      deallocate(b2srsm_she0to3)
      deallocate(b2srsm_shi0to3)
      deallocate(b2srsm_sch)
      deallocate(b2srdt_sna)
      deallocate(b2srdt_smo)
      deallocate(b2srdt_she)
      deallocate(b2srdt_shi)
      deallocate(b2srdt_she0to3)
      deallocate(b2srdt_shi0to3)
      deallocate(b2srdt_sch)
      deallocate(b2srst_sna)
      deallocate(b2srst_smo)
      deallocate(b2srst_she)
      deallocate(b2srst_shi)
      deallocate(b2srst_she0to3)
      deallocate(b2srst_shi0to3)
      deallocate(b2srst_sch)
      deallocate(b2sihs_diae)
      deallocate(b2sihs_diaa)
      deallocate(b2npht_shei)

      return
      end subroutine dealloc_b2mod_sources

      subroutine alloc_b2mod_eirene_sources(nxd,nyd,nsd,nstrad)
      implicit none
      integer nxd,nyd,nsd,nstrad

      allocate(eirene_mc_t0_papl(-1:nxd,-1:nyd,0:nsd-1,1:nstrad)) ! djm Jan2017
      allocate(eirene_mc_t0_pmpl(-1:nxd,-1:nyd,0:nsd-1,1:nstrad))
      allocate(eirene_mc_t0_pipl(-1:nxd,-1:nyd,0:nsd-1,1:nstrad))
      allocate(eirene_mc_t0_pppl(-1:nxd,-1:nyd,0:nsd-1,1:nstrad))
      allocate(eirene_mc_t0_core(-1:nxd,-1:nyd,0:nsd-1))
      allocate(eirene_mc_t0_mapl(-1:nxd,-1:nyd,0:nsd-1,1:nstrad))
      allocate(eirene_mc_t0_mmpl(-1:nxd,-1:nyd,0:nsd-1,1:nstrad))
      allocate(eirene_mc_t0_mipl(-1:nxd,-1:nyd,0:nsd-1,1:nstrad))
      allocate(eirene_mc_t0_cppv(-1:nxd,-1:nyd,0:nsd-1,1:nstrad))
      allocate(eirene_mc_t0_eael(-1:nxd,-1:nyd,1:nstrad))
      allocate(eirene_mc_t0_emel(-1:nxd,-1:nyd,1:nstrad))
      allocate(eirene_mc_t0_eiel(-1:nxd,-1:nyd,1:nstrad))
      allocate(eirene_mc_t0_epel(-1:nxd,-1:nyd,1:nstrad))
      allocate(eirene_mc_t0_eapl(-1:nxd,-1:nyd,1:nstrad))
      allocate(eirene_mc_t0_empl(-1:nxd,-1:nyd,1:nstrad))
      allocate(eirene_mc_t0_eipl(-1:nxd,-1:nyd,1:nstrad))
      allocate(eirene_mc_t0_eppl(-1:nxd,-1:nyd,1:nstrad))

      eirene_mc_t0_papl=0.0_R8 !djm Jan2017
      eirene_mc_t0_pmpl=0.0_R8
      eirene_mc_t0_pipl=0.0_R8
      eirene_mc_t0_pppl=0.0_R8
      eirene_mc_t0_core=0.0_R8
      eirene_mc_t0_mapl=0.0_R8
      eirene_mc_t0_mmpl=0.0_R8
      eirene_mc_t0_mipl=0.0_R8
      eirene_mc_t0_cppv=0.0_R8
      eirene_mc_t0_eael=0.0_R8
      eirene_mc_t0_emel=0.0_R8
      eirene_mc_t0_eiel=0.0_R8
      eirene_mc_t0_epel=0.0_R8
      eirene_mc_t0_eapl=0.0_R8
      eirene_mc_t0_empl=0.0_R8
      eirene_mc_t0_eipl=0.0_R8
      eirene_mc_t0_eppl=0.0_R8

      return
      end subroutine alloc_b2mod_eirene_sources

      subroutine dealloc_b2mod_eirene_sources
      implicit none

      if(.not.allocated(eirene_mc_t0_papl)) return
      deallocate(eirene_mc_t0_papl)
      deallocate(eirene_mc_t0_pmpl)
      deallocate(eirene_mc_t0_pipl)
      deallocate(eirene_mc_t0_pppl)
      deallocate(eirene_mc_t0_core)
      deallocate(eirene_mc_t0_mapl)
      deallocate(eirene_mc_t0_mmpl)
      deallocate(eirene_mc_t0_mipl)
      deallocate(eirene_mc_t0_cppv)
      deallocate(eirene_mc_t0_eael)
      deallocate(eirene_mc_t0_emel)
      deallocate(eirene_mc_t0_eiel)
      deallocate(eirene_mc_t0_epel)
      deallocate(eirene_mc_t0_eapl)
      deallocate(eirene_mc_t0_empl)
      deallocate(eirene_mc_t0_eipl)
      deallocate(eirene_mc_t0_eppl)

      return
      end subroutine dealloc_b2mod_eirene_sources

      subroutine alloc_b2mod_eirene_sources(nxd,nyd,nsd,nstrad)
      implicit none
      integer nxd,nyd,nsd,nstrad

      allocate(eirene_mc_t0_papl(-1:nxd,-1:nyd,0:nsd-1,1:nstrad)) ! djm Jan2017
      allocate(eirene_mc_t0_pmpl(-1:nxd,-1:nyd,0:nsd-1,1:nstrad))
      allocate(eirene_mc_t0_pipl(-1:nxd,-1:nyd,0:nsd-1,1:nstrad))
      allocate(eirene_mc_t0_pppl(-1:nxd,-1:nyd,0:nsd-1,1:nstrad))
      allocate(eirene_mc_t0_core(-1:nxd,-1:nyd,0:nsd-1))
      allocate(eirene_mc_t0_mapl(-1:nxd,-1:nyd,0:nsd-1,1:nstrad))
      allocate(eirene_mc_t0_mmpl(-1:nxd,-1:nyd,0:nsd-1,1:nstrad))
      allocate(eirene_mc_t0_mipl(-1:nxd,-1:nyd,0:nsd-1,1:nstrad))
      allocate(eirene_mc_t0_cppv(-1:nxd,-1:nyd,0:nsd-1,1:nstrad))
      allocate(eirene_mc_t0_eael(-1:nxd,-1:nyd,1:nstrad))
      allocate(eirene_mc_t0_emel(-1:nxd,-1:nyd,1:nstrad))
      allocate(eirene_mc_t0_eiel(-1:nxd,-1:nyd,1:nstrad))
      allocate(eirene_mc_t0_epel(-1:nxd,-1:nyd,1:nstrad))
      allocate(eirene_mc_t0_eapl(-1:nxd,-1:nyd,1:nstrad))
      allocate(eirene_mc_t0_empl(-1:nxd,-1:nyd,1:nstrad))
      allocate(eirene_mc_t0_eipl(-1:nxd,-1:nyd,1:nstrad))
      allocate(eirene_mc_t0_eppl(-1:nxd,-1:nyd,1:nstrad))

      eirene_mc_t0_papl=0.0_R8 !djm Jan2017
      eirene_mc_t0_pmpl=0.0_R8
      eirene_mc_t0_pipl=0.0_R8
      eirene_mc_t0_pppl=0.0_R8
      eirene_mc_t0_core=0.0_R8
      eirene_mc_t0_mapl=0.0_R8
      eirene_mc_t0_mmpl=0.0_R8
      eirene_mc_t0_mipl=0.0_R8
      eirene_mc_t0_cppv=0.0_R8
      eirene_mc_t0_eael=0.0_R8
      eirene_mc_t0_emel=0.0_R8
      eirene_mc_t0_eiel=0.0_R8
      eirene_mc_t0_epel=0.0_R8
      eirene_mc_t0_eapl=0.0_R8
      eirene_mc_t0_empl=0.0_R8
      eirene_mc_t0_eipl=0.0_R8
      eirene_mc_t0_eppl=0.0_R8

      return
      end subroutine alloc_b2mod_eirene_sources

      subroutine dealloc_b2mod_eirene_sources
      implicit none

      if(.not.allocated(eirene_mc_t0_papl)) return
      deallocate(eirene_mc_t0_papl)
      deallocate(eirene_mc_t0_pmpl)
      deallocate(eirene_mc_t0_pipl)
      deallocate(eirene_mc_t0_pppl)
      deallocate(eirene_mc_t0_core)
      deallocate(eirene_mc_t0_mapl)
      deallocate(eirene_mc_t0_mmpl)
      deallocate(eirene_mc_t0_mipl)
      deallocate(eirene_mc_t0_cppv)
      deallocate(eirene_mc_t0_eael)
      deallocate(eirene_mc_t0_emel)
      deallocate(eirene_mc_t0_eiel)
      deallocate(eirene_mc_t0_epel)
      deallocate(eirene_mc_t0_eapl)
      deallocate(eirene_mc_t0_empl)
      deallocate(eirene_mc_t0_eipl)
      deallocate(eirene_mc_t0_eppl)

      return
      end subroutine dealloc_b2mod_eirene_sources

      subroutine read_b2mod_sources(nread)
      implicit none
      integer nread

      call cfrure (nread, size(sch), sch, 'sch')
      call cfrure (nread, size(she), she, 'she')
      call cfrure (nread, size(shi), shi, 'shi')
      call cfrure (nread, size(smo), smo, 'smo')
      call cfrure (nread, size(smq), smq, 'smq')
      call cfrure (nread, size(sna), sna, 'sna')
      if(b2fplasma_version.ge.'01.001.026')
     >  call cfrure_opt (nread, size(sne), sne, 'sne')
      if(b2fplasma_version.gt.'00.000.000') then
        call cfrure (nread, size(rsana), rsana, 'rsana')
        call cfrure (nread, size(rsahi), rsahi, 'rsahi')
        if(b2fplasma_version.gt.'01.001.003') then
          call cfrure (nread, size(rsamo), rsamo, 'rsamo')
        endif
        call cfrure (nread, size(rrana), rrana, 'rrana')
        call cfrure (nread, size(rrahi), rrahi, 'rrahi')
        if(b2fplasma_version.gt.'01.001.003') then
          call cfrure (nread, size(rramo), rramo, 'rramo')
        endif
        call cfrure (nread, size(rqahe), rqahe, 'rqahe')
        if(b2fplasma_version.ge.'01.001.016') then
          call cfrure (nread, size(rqrad), rqrad, 'rqrad')
        endif
        if(b2fplasma_version.ge.'01.001.019') then
          call cfrure (nread, size(rqbrm), rqbrm, 'rqbrm')
        endif
        call cfrure (nread, size(rcxna), rcxna, 'rcxna')
        call cfrure (nread, size(rcxhi), rcxhi, 'rcxhi')
        if(b2fplasma_version.gt.'01.001.003') then
          call cfrure (nread, size(rcxmo), rcxmo, 'rcxmo')
        endif
      endif
      if(b2fplasma_version.gt.'01.001.000') then
        call cfrure (nread, size(b2stbr_sna), b2stbr_sna, 'b2stbr_sna')
        call cfrure (nread, size(b2stbr_smo), b2stbr_smo, 'b2stbr_smo')
        call cfrure (nread, size(b2stbr_she), b2stbr_she, 'b2stbr_she')
        call cfrure (nread, size(b2stbr_shi), b2stbr_shi, 'b2stbr_shi')
        call cfrure (nread, size(b2stbr_sch), b2stbr_sch, 'b2stbr_sch')
      endif
      if(b2fplasma_version.ge.'01.001.026') call cfrure_opt 
     >  (nread, size(b2stbr_sne), b2stbr_sne, 'b2stbr_sne')
      if(b2fplasma_version.gt.'01.001.010') then
        call cfrure (nread, size(b2stbc_sna), b2stbc_sna, 'b2stbc_sna')
        call cfrure (nread, size(b2stbc_smo), b2stbc_smo, 'b2stbc_smo')
        call cfrure (nread, size(b2stbc_she), b2stbc_she, 'b2stbc_she')
        call cfrure (nread, size(b2stbc_shi), b2stbc_shi, 'b2stbc_shi')
        call cfrure (nread, size(b2stbc_sch), b2stbc_sch, 'b2stbc_sch')
      endif
      if(b2fplasma_version.ge.'03.000.006') then !djm Jan2017
        call cfrure (nread, size(b2stel_sna_ion),
     .                b2stel_sna_ion, 'b2stel_sna_ion')
        call cfrure (nread, size(b2stel_sna_rec),
     .                b2stel_sna_rec, 'b2stel_sna_rec')
        call cfrure (nread, size(b2stel_smq_ion),
     .                b2stel_smq_ion, 'b2stel_smq_ion')
        call cfrure (nread, size(b2stel_smq_rec),
     .                b2stel_smq_rec, 'b2stel_smq_rec')
        call cfrure (nread, size(b2stel_she),
     .                b2stel_she, 'b2stel_she')
        call cfrure (nread, size(b2stel_shi_ion),
     .                b2stel_shi_ion, 'b2stel_shi_ion')
        call cfrure (nread, size(b2stel_shi_rec),
     .                b2stel_shi_rec, 'b2stel_shi_rec')
        call cfrure (nread, size(b2stcx_sna), b2stcx_sna, 'b2stcx_sna')
        call cfrure (nread, size(b2stcx_smq), b2stcx_smq, 'b2stcx_smq')
        call cfrure (nread, size(b2stcx_shi), b2stcx_shi, 'b2stcx_shi')
        call cfrure (nread, size(b2srsm_sna), b2srsm_sna, 'b2srsm_sna')
        call cfrure (nread, size(b2srsm_smo), b2srsm_smo, 'b2srsm_smo')
        call cfrure (nread, size(b2srsm_she), b2srsm_she, 'b2srsm_she')
        call cfrure (nread, size(b2srsm_shi), b2srsm_shi, 'b2srsm_shi')
        call cfrure (nread, size(b2srsm_sch), b2srsm_sch, 'b2srsm_sch')
        call cfrure (nread, size(b2srdt_sna), b2srdt_sna, 'b2srdt_sna')
        call cfrure (nread, size(b2srdt_smo), b2srdt_smo, 'b2srdt_smo')
        call cfrure (nread, size(b2srdt_she), b2srdt_she, 'b2srdt_she')
        call cfrure (nread, size(b2srdt_shi), b2srdt_shi, 'b2srdt_shi')
        call cfrure (nread, size(b2srdt_sch), b2srdt_sch, 'b2srdt_sch')
        call cfrure (nread, size(b2srst_sna), b2srst_sna, 'b2srst_sna')
        call cfrure (nread, size(b2srst_smo), b2srst_smo, 'b2srst_smo')
        call cfrure (nread, size(b2srst_she), b2srst_she, 'b2srst_she')
        call cfrure (nread, size(b2srst_shi), b2srst_shi, 'b2srst_shi')
        call cfrure (nread, size(b2srst_sch), b2srst_sch, 'b2srst_sch')
        call cfrure (nread, size(b2sifr_smoch),
     &               b2sifr_smoch, 'b2sifr_smoch')
        call cfrure (nread, size(b2sifr_smotf),
     &               b2sifr_smotf, 'b2sifr_smotf')
        call cfrure (nread, size(b2sigp_smogp),
     &               b2sigp_smogp, 'b2sigp_smogp')
        call cfrure (nread, size(b2sigp_pstat),
     &               b2sigp_pstat, 'b2sigp_pstat')
        call cfrure (nread, size(b2siav_smovh),
     &               b2siav_smovh, 'b2siav_smovh')
        call cfrure (nread, size(b2siav_smovv),
     &               b2siav_smovv, 'b2siav_smovv')
        call cfrure (nread, size(b2nxdv_smo),
     &               b2nxdv_smo, 'b2nxdv_smo')
        call cfrure (nread, size(b2sihs_diae),
     &               b2sihs_diae, 'b2sihs_diae')
        call cfrure (nread, size(b2sihs_diaa),
     &               b2sihs_diaa, 'b2sihs_diaa')
        call cfrure (nread, size(b2npht_shei),
     &               b2npht_shei, 'b2npht_shei')
      endif
      if(b2fplasma_version.ge.'01.001.026') call cfrure_opt 
     >  (nread, size(b2stbc_sne), b2stbc_sne, 'b2stbc_sne')
      if(b2fplasma_version.ge.'01.001.027') then
        call cfrure_opt 
     >   (nread, size(b2stbm_sna), b2stbm_sna, 'b2stbm_sna')
        call cfrure_opt 
     >   (nread, size(b2stbm_smo), b2stbm_smo, 'b2stbm_smo')
        call cfrure_opt 
     >   (nread, size(b2stbm_she), b2stbm_she, 'b2stbm_she')
        call cfrure_opt 
     >   (nread, size(b2stbm_shi), b2stbm_shi, 'b2stbm_shi')
        call cfrure_opt 
     >   (nread, size(b2stbm_sch), b2stbm_sch, 'b2stbm_sch')
        call cfrure_opt 
     >   (nread, size(b2stbm_sne), b2stbm_sne, 'b2stbm_sne')
      endif
      if(b2fplasma_version.gt.'01.001.002') then
        call cfrure (nread, size(b2sihs_divue), b2sihs_divue,
     1   'b2sihs_divue') 
        call cfrure (nread, size(b2sihs_divua), b2sihs_divua,
     1   'b2sihs_divua') 
        call cfrure (nread, size(b2sihs_exbe),  b2sihs_exbe,
     1   'b2sihs_exbe')  
        call cfrure (nread, size(b2sihs_exba),  b2sihs_exba,
     1   'b2sihs_exba')  
        call cfrure (nread, size(b2sihs_visa),  b2sihs_visa,
     1   'b2sihs_visa')  
        call cfrure (nread, size(b2sihs_joule), b2sihs_joule,
     1   'b2sihs_joule') 
        call cfrure (nread, size(b2sihs_fraa),  b2sihs_fraa,
     1   'b2sihs_fraa')  
        call cfrure (nread, size(b2sihs_str),   b2sihs_str,
     1   'b2sihs_str')   
      endif
      if(b2fplasma_version.gt.'01.001.005') then
        call cfrure (nread, size(b2npmo_smaf),  b2npmo_smaf,
     1   'b2npmo_smaf')
        call cfrure (nread, size(b2npmo_smag),  b2npmo_smag,
     1   'b2npmo_smag')
        call cfrure (nread, size(b2npmo_smav),  b2npmo_smav,
     1   'b2npmo_smav')
      endif
      if(b2fplasma_version.gt.'01.001.030') then
         call cfrure_opt (nread, size(smpr),  smpr, 'smpr')
         call cfrure_opt (nread, size(smpt),  smpt, 'smpt')
         call cfrure_opt (nread, size(smfr),  smfr, 'smfr')
      endif
      if(b2fplasma_version.gt.'03.000.000') then
         call cfrure (nread, size(smcf),  smcf, 'smcf')
      endif
      if(b2fplasma_version.gt.'01.001.012') then
        call cfrure (nread, size(ext_sna), ext_sna, 'ext_sna')
        call cfrure (nread, size(ext_smo), ext_smo, 'ext_smo')
        call cfrure (nread, size(ext_she), ext_she, 'ext_she')
        call cfrure (nread, size(ext_shi), ext_shi, 'ext_shi')
        call cfrure (nread, size(ext_sch), ext_sch, 'ext_sch')
      endif
      if(b2fplasma_version.ge.'01.001.026')
     >  call cfrure_opt (nread, size(ext_sne), ext_sne, 'ext_sne')

      return
      end subroutine read_b2mod_sources

      subroutine read_b2mod_eirene_sources(nxd,nyd,nsd,nstrad,nread)
      implicit none
      integer nxd,nyd,nsd,nstrad,nread

      if(b2fplasma_version.ge.'03.000.006') then
        if (.not.allocated(eirene_mc_t0_papl))
     >   call alloc_b2mod_eirene_sources(nxd,nyd,nsd,nstrad)
        call cfrure_opt (nread, size(eirene_mc_t0_papl),
     &                   eirene_mc_t0_papl, 'eirene_mc_t0_papl')
        call cfrure_opt (nread, size(eirene_mc_t0_pmpl),
     &                   eirene_mc_t0_pmpl, 'eirene_mc_t0_pmpl')
        call cfrure_opt (nread, size(eirene_mc_t0_pipl),
     &                   eirene_mc_t0_pipl, 'eirene_mc_t0_pipl')
        call cfrure_opt (nread, size(eirene_mc_t0_pppl),
     &                   eirene_mc_t0_pppl, 'eirene_mc_t0_pppl')
        call cfrure_opt (nread, size(eirene_mc_t0_core),
     &                   eirene_mc_t0_core, 'eirene_mc_t0_core')
        call cfrure_opt (nread, size(eirene_mc_t0_mapl),
     &                   eirene_mc_t0_mapl, 'eirene_mc_t0_mapl')
        call cfrure_opt (nread, size(eirene_mc_t0_mmpl),
     &                   eirene_mc_t0_mmpl, 'eirene_mc_t0_mmpl')
        call cfrure_opt (nread, size(eirene_mc_t0_mipl),
     &                   eirene_mc_t0_mipl, 'eirene_mc_t0_mipl')
        call cfrure_opt (nread, size(eirene_mc_t0_cppv),
     &                   eirene_mc_t0_cppv, 'eirene_mc_t0_cppv')
        call cfrure_opt (nread, size(eirene_mc_t0_eael),
     &                   eirene_mc_t0_eael, 'eirene_mc_t0_eael')
        call cfrure_opt (nread, size(eirene_mc_t0_emel),
     &                   eirene_mc_t0_emel, 'eirene_mc_t0_emel')
        call cfrure_opt (nread, size(eirene_mc_t0_eiel),
     &                   eirene_mc_t0_eiel, 'eirene_mc_t0_eiel')
        call cfrure_opt (nread, size(eirene_mc_t0_epel),
     &                   eirene_mc_t0_epel, 'eirene_mc_t0_epel')
        call cfrure_opt (nread, size(eirene_mc_t0_eapl),
     &                   eirene_mc_t0_eapl, 'eirene_mc_t0_eapl')
        call cfrure_opt (nread, size(eirene_mc_t0_empl),
     &                   eirene_mc_t0_empl, 'eirene_mc_t0_empl')
        call cfrure_opt (nread, size(eirene_mc_t0_eipl),
     &                   eirene_mc_t0_eipl, 'eirene_mc_t0_eipl')
        call cfrure_opt (nread, size(eirene_mc_t0_eppl),
     &                   eirene_mc_t0_eppl, 'eirene_mc_t0_eppl')
      end if

      return
      end subroutine read_b2mod_eirene_sources

      subroutine write_b2mod_sources(nwrite)
      implicit none
      integer nwrite
      character*3 source_name

      if (newversion.ge.'01.001.026') then
        source_name = 'sne'
      else if (newversion.gt.'01.001.010') then
        source_name = 'sch'
      end if

      call cfwure (nwrite, size(sch), sch, 'sch')
      call cfwure (nwrite, size(she), she, 'she')
      call cfwure (nwrite, size(shi), shi, 'shi')
      call cfwure (nwrite, size(smo), smo, 'smo')
      call cfwure (nwrite, size(smq), smq, 'smq')
      call cfwure (nwrite, size(sna), sna, 'sna')
      if(newversion.ge.'01.001.026')
     >  call cfwure (nwrite, size(sne), sne, 'sne')
      if(newversion.gt.'00.000.000') then
        write(*,*) newversion,'rsana --- rcxhi'
        call cfwure (nwrite, size(rsana), rsana, 'rsana')
        call cfwure (nwrite, size(rsahi), rsahi, 'rsahi')
        if(newversion.gt.'01.001.003') then
          write(*,*) newversion,'rsamo --- rcxmo'
          call cfwure (nwrite, size(rsamo), rsamo, 'rsamo')
        endif
        call cfwure (nwrite, size(rrana), rrana, 'rrana')
        call cfwure (nwrite, size(rrahi), rrahi, 'rrahi')
        if(newversion.gt.'01.001.003') then
          call cfwure (nwrite, size(rramo), rramo, 'rramo')
        endif
        call cfwure (nwrite, size(rqahe), rqahe, 'rqahe')
        if(newversion.ge.'01.001.016') then
          call cfwure (nwrite, size(rqrad), rqrad, 'rqrad')
        endif
        if(newversion.ge.'01.001.019') then
          call cfwure (nwrite, size(rqbrm), rqbrm, 'rqbrm')
        endif
        call cfwure (nwrite, size(rcxna), rcxna, 'rcxna')
        call cfwure (nwrite, size(rcxhi), rcxhi, 'rcxhi')
        if(newversion.gt.'01.001.003') then
          call cfwure (nwrite, size(rcxmo), rcxmo, 'rcxmo')
        endif
      endif
      if(newversion.gt.'01.001.000') then
        write(*,*) newversion,'b2stbr_sna --- _'//source_name
        call cfwure (nwrite, size(b2stbr_sna), b2stbr_sna, 'b2stbr_sna')
        call cfwure (nwrite, size(b2stbr_smo), b2stbr_smo, 'b2stbr_smo')
        call cfwure (nwrite, size(b2stbr_she), b2stbr_she, 'b2stbr_she')
        call cfwure (nwrite, size(b2stbr_shi), b2stbr_shi, 'b2stbr_shi')
        call cfwure (nwrite, size(b2stbr_sch), b2stbr_sch, 'b2stbr_sch')
      endif
      if(newversion.ge.'01.001.026')
     >  call cfwure (nwrite, size(b2stbr_sne), b2stbr_sne, 'b2stbr_sne')
      if(newversion.gt.'01.001.010') then
        write(*,*) newversion,'b2stbc_sna --- _'//source_name
        call cfwure (nwrite, size(b2stbc_sna), b2stbc_sna, 'b2stbc_sna')
        call cfwure (nwrite, size(b2stbc_smo), b2stbc_smo, 'b2stbc_smo')
        call cfwure (nwrite, size(b2stbc_she), b2stbc_she, 'b2stbc_she')
        call cfwure (nwrite, size(b2stbc_shi), b2stbc_shi, 'b2stbc_shi')
        call cfwure (nwrite, size(b2stbc_sch), b2stbc_sch, 'b2stbc_sch')
      endif
      if(newversion.ge.'03.000.006') then !djm Jan2017
        call cfwure (nwrite, size(b2stel_sna_ion),
     .                b2stel_sna_ion, 'b2stel_sna_ion')
        call cfwure (nwrite, size(b2stel_sna_rec),
     .                b2stel_sna_rec, 'b2stel_sna_rec')
        call cfwure (nwrite, size(b2stel_smq_ion),
     .                b2stel_smq_ion, 'b2stel_smq_ion')
        call cfwure (nwrite, size(b2stel_smq_rec),
     .                b2stel_smq_rec, 'b2stel_smq_rec')
        call cfwure (nwrite, size(b2stel_she),
     .                b2stel_she, 'b2stel_she')
        call cfwure (nwrite, size(b2stel_shi_ion),
     .                b2stel_shi_ion, 'b2stel_shi_ion')
        call cfwure (nwrite, size(b2stel_shi_rec),
     .                b2stel_shi_rec, 'b2stel_shi_rec')
        call cfwure (nwrite, size(b2stcx_sna), b2stcx_sna, 'b2stcx_sna')
        call cfwure (nwrite, size(b2stcx_smq), b2stcx_smq, 'b2stcx_smq')
        call cfwure (nwrite, size(b2stcx_shi), b2stcx_shi, 'b2stcx_shi')
        call cfwure (nwrite, size(b2srsm_sna), b2srsm_sna, 'b2srsm_sna')
        call cfwure (nwrite, size(b2srsm_smo), b2srsm_smo, 'b2srsm_smo')
        call cfwure (nwrite, size(b2srsm_she), b2srsm_she, 'b2srsm_she')
        call cfwure (nwrite, size(b2srsm_shi), b2srsm_shi, 'b2srsm_shi')
        call cfwure (nwrite, size(b2srsm_sch), b2srsm_sch, 'b2srsm_sch')
        call cfwure (nwrite, size(b2srdt_sna), b2srdt_sna, 'b2srdt_sna')
        call cfwure (nwrite, size(b2srdt_smo), b2srdt_smo, 'b2srdt_smo')
        call cfwure (nwrite, size(b2srdt_she), b2srdt_she, 'b2srdt_she')
        call cfwure (nwrite, size(b2srdt_shi), b2srdt_shi, 'b2srdt_shi')
        call cfwure (nwrite, size(b2srdt_sch), b2srdt_sch, 'b2srdt_sch')
        call cfwure (nwrite, size(b2srst_sna), b2srst_sna, 'b2srst_sna')
        call cfwure (nwrite, size(b2srst_smo), b2srst_smo, 'b2srst_smo')
        call cfwure (nwrite, size(b2srst_she), b2srst_she, 'b2srst_she')
        call cfwure (nwrite, size(b2srst_shi), b2srst_shi, 'b2srst_shi')
        call cfwure (nwrite, size(b2srst_sch), b2srst_sch, 'b2srst_sch')
        call cfwure (nwrite, size(b2sifr_smoch),
     &               b2sifr_smoch, 'b2sifr_smoch')
        call cfwure (nwrite, size(b2sifr_smotf),
     &               b2sifr_smotf, 'b2sifr_smotf')
        call cfwure (nwrite, size(b2sigp_smogp),
     &               b2sigp_smogp, 'b2sigp_smogp')
        call cfwure (nwrite, size(b2sigp_pstat),
     &               b2sigp_pstat, 'b2sigp_pstat')
        call cfwure (nwrite, size(b2siav_smovh),
     &               b2siav_smovh, 'b2siav_smovh')
        call cfwure (nwrite, size(b2siav_smovv),
     &               b2siav_smovv, 'b2siav_smovv')
        call cfwure (nwrite, size(b2nxdv_smo),
     &               b2nxdv_smo, 'b2nxdv_smo')
        call cfwure (nwrite, size(b2sihs_diae),
     &               b2sihs_diae, 'b2sihs_diae')
        call cfwure (nwrite, size(b2sihs_diaa),
     &               b2sihs_diaa, 'b2sihs_diaa')
        call cfwure (nwrite, size(b2npht_shei),
     &               b2npht_shei, 'b2npht_shei')
      endif
      if(newversion.ge.'01.001.026')
     >  call cfwure (nwrite, size(b2stbc_sne), b2stbc_sne, 'b2stbc_sne')
      if(newversion.ge.'01.001.027') then
        write(*,*) newversion,'b2stbm_sna --- _'//source_name
        call cfwure (nwrite, size(b2stbm_sna), b2stbm_sna, 'b2stbm_sna')
        call cfwure (nwrite, size(b2stbm_smo), b2stbm_smo, 'b2stbm_smo')
        call cfwure (nwrite, size(b2stbm_she), b2stbm_she, 'b2stbm_she')
        call cfwure (nwrite, size(b2stbm_shi), b2stbm_shi, 'b2stbm_shi')
        call cfwure (nwrite, size(b2stbm_sch), b2stbm_sch, 'b2stbm_sch')
        call cfwure (nwrite, size(b2stbm_sne), b2stbm_sne, 'b2stbm_sne')
      endif
       if(newversion.gt.'01.001.002') then
        write(*,*) newversion,'b2sihs_divue --- _str'
        call cfwure (nwrite, size(b2sihs_divue),b2sihs_divue,
     1   'b2sihs_divue') 
        call cfwure (nwrite, size(b2sihs_divua),b2sihs_divua,
     1   'b2sihs_divua') 
        call cfwure (nwrite, size(b2sihs_exbe), b2sihs_exbe,
     1   'b2sihs_exbe')  
        call cfwure (nwrite, size(b2sihs_exba), b2sihs_exba,
     1   'b2sihs_exba')  
        call cfwure (nwrite, size(b2sihs_visa), b2sihs_visa,
     1   'b2sihs_visa')  
        call cfwure (nwrite, size(b2sihs_joule),b2sihs_joule,
     1   'b2sihs_joule') 
        call cfwure (nwrite, size(b2sihs_fraa), b2sihs_fraa,
     1   'b2sihs_fraa')  
        call cfwure (nwrite, size(b2sihs_str),  b2sihs_str,
     1   'b2sihs_str')   
      endif
      if(newversion.gt.'01.001.005') then
        call cfwure (nwrite, size(b2npmo_smaf), b2npmo_smaf,
     1   'b2npmo_smaf')
        call cfwure (nwrite, size(b2npmo_smag), b2npmo_smag,
     1   'b2npmo_smag')
        call cfwure (nwrite, size(b2npmo_smav), b2npmo_smav,
     1   'b2npmo_smav')
      endif
      if(newversion.gt.'01.001.030') then
         call cfwure (nwrite, size(smpr),  smpr, 'smpr')
         call cfwure (nwrite, size(smpt),  smpt, 'smpt')
         call cfwure (nwrite, size(smfr),  smfr, 'smfr')
      endif
      if(newversion.gt.'03.000.000') then
         call cfwure (nwrite, size(smcf),  smcf, 'smcf')
      endif
      if(newversion.gt.'01.001.012') then
        call cfwure (nwrite, size(ext_sna), ext_sna, 'ext_sna')
        call cfwure (nwrite, size(ext_smo), ext_smo, 'ext_smo')
        call cfwure (nwrite, size(ext_she), ext_she, 'ext_she')
        call cfwure (nwrite, size(ext_shi), ext_shi, 'ext_shi')
        call cfwure (nwrite, size(ext_sch), ext_sch, 'ext_sch')
      endif
      if(newversion.ge.'01.001.026')
     >  call cfwure (nwrite, size(ext_sne), ext_sne, 'ext_sne')

      return
      end subroutine write_b2mod_sources

      subroutine write_b2mod_eirene_sources(nwrite)
      implicit none
      integer nwrite

      if (.not.allocated(eirene_mc_t0_papl)) return
      call cfwure (nwrite, size(eirene_mc_t0_papl),
     &             eirene_mc_t0_papl, 'eirene_mc_t0_papl')
      call cfwure (nwrite, size(eirene_mc_t0_pmpl),
     &             eirene_mc_t0_pmpl, 'eirene_mc_t0_pmpl')
      call cfwure (nwrite, size(eirene_mc_t0_pipl),
     &             eirene_mc_t0_pipl, 'eirene_mc_t0_pipl')
      call cfwure (nwrite, size(eirene_mc_t0_pppl),
     &             eirene_mc_t0_pppl, 'eirene_mc_t0_pppl')
      call cfwure (nwrite, size(eirene_mc_t0_core),
     &             eirene_mc_t0_core, 'eirene_mc_t0_core')
      call cfwure (nwrite, size(eirene_mc_t0_mapl),
     &             eirene_mc_t0_mapl, 'eirene_mc_t0_mapl')
      call cfwure (nwrite, size(eirene_mc_t0_mmpl),
     &             eirene_mc_t0_mmpl, 'eirene_mc_t0_mmpl')
      call cfwure (nwrite, size(eirene_mc_t0_mipl),
     &             eirene_mc_t0_mipl, 'eirene_mc_t0_mipl')
      call cfwure (nwrite, size(eirene_mc_t0_cppv),
     &             eirene_mc_t0_cppv, 'eirene_mc_t0_cppv')
      call cfwure (nwrite, size(eirene_mc_t0_eael),
     &             eirene_mc_t0_eael, 'eirene_mc_t0_eael')
      call cfwure (nwrite, size(eirene_mc_t0_emel),
     &             eirene_mc_t0_emel, 'eirene_mc_t0_emel')
      call cfwure (nwrite, size(eirene_mc_t0_eiel),
     &             eirene_mc_t0_eiel, 'eirene_mc_t0_eiel')
      call cfwure (nwrite, size(eirene_mc_t0_epel),
     &             eirene_mc_t0_epel, 'eirene_mc_t0_epel')
      call cfwure (nwrite, size(eirene_mc_t0_eapl),
     &             eirene_mc_t0_eapl, 'eirene_mc_t0_eapl')
      call cfwure (nwrite, size(eirene_mc_t0_empl),
     &             eirene_mc_t0_empl, 'eirene_mc_t0_empl')
      call cfwure (nwrite, size(eirene_mc_t0_eipl),
     &             eirene_mc_t0_eipl, 'eirene_mc_t0_eipl')
      call cfwure (nwrite, size(eirene_mc_t0_eppl),
     &             eirene_mc_t0_eppl, 'eirene_mc_t0_eppl')

      return
      end subroutine write_b2mod_eirene_sources

      end module b2mod_sources
