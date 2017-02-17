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
      real (kind=R8), allocatable, save ::                             !srv 11.09.09 16.07.10
     *  snasm(:,:,:,:), smosm(:,:,:,:), 
     &  shesm(:,:,:), shism(:,:,:),
     *  snadt(:,:,:,:), smodt(:,:,:,:), 
     &  shedt(:,:,:), shidt(:,:,:),
     &  schdt(:,:,:), snedt(:,:,:)
cdjm Jan2017 Source arrays required for balance. The linearised components are stored
c            and then the total source is updated just before each residual calculation.
c            Only this total source is output to b2fplasma.
cdjm Jan2017 Arrays for particle balance:
      real (kind=R8), allocatable, save ::
     *  b2stbr_sna0to1(:,:,:,:),
     *  b2stbr_sna_bal(:,:,:),
     *  b2stbc_sna0to1(:,:,:,:),
     *  b2stbc_sna_bal(:,:,:),
     *  eirene_mc_papl_sna0to1(:,:,:,:,:),
     *  eirene_mc_pmpl_sna0to1(:,:,:,:,:),
     *  eirene_mc_pipl_sna0to1(:,:,:,:,:),
     *  eirene_mc_pppl_sna0to1(:,:,:,:,:),
     *  eirene_mc_core_sna0to1(:,:,:,:),
     *  eirene_mc_papl_sna(:,:,:,:),
     *  eirene_mc_pmpl_sna(:,:,:,:),
     *  eirene_mc_pipl_sna(:,:,:,:),
     *  eirene_mc_pppl_sna(:,:,:,:),
     *  eirene_mc_core_sna(:,:,:),
     *  b2stbm_sna0to1(:,:,:,:),
     *  b2stbm_sna_bal(:,:,:),
     *  ext_sna0to1(:,:,:,:),
     *  ext_sna_bal(:,:,:),
     *  b2stel_sna_ion0to1(:,:,:,:),
     *  b2stel_sna_ion_bal(:,:,:),
     *  b2stel_sna_rec0to1(:,:,:,:),
     *  b2stel_sna_rec_bal(:,:,:),
     *  b2stcx_sna0to1(:,:,:,:),
     *  b2stcx_sna_bal(:,:,:),
     *  b2srsm_sna0to1(:,:,:,:),
     *  b2srsm_sna_bal(:,:,:),
     *  b2srdt_sna0to1(:,:,:,:),
     *  b2srdt_sna_bal(:,:,:),
     *  b2srst_sna0to1(:,:,:,:),
     *  b2srst_sna_bal(:,:,:)
cdjm Jan2017 Arrays for momentum balance:
      real (kind=R8), allocatable, save ::
     *  b2stbr_smo0to3(:,:,:,:),
     *  b2stbr_smo_bal(:,:,:),
     *  b2stbc_smo0to3(:,:,:,:),
     *  b2stbc_smo_bal(:,:,:),
     *  eirene_mc_mapl_smo0to3(:,:,:,:,:),
     *  eirene_mc_mmpl_smo0to3(:,:,:,:,:),
     *  eirene_mc_mipl_smo0to3(:,:,:,:,:),
     *  eirene_mc_cppv_smo0to3(:,:,:,:,:),
     *  eirene_mc_mapl_smo(:,:,:,:),
     *  eirene_mc_mmpl_smo(:,:,:,:),
     *  eirene_mc_mipl_smo(:,:,:,:),
     *  eirene_mc_cppv_smo(:,:,:,:),
     *  b2stbm_smo0to3(:,:,:,:),
     *  b2stbm_smo_bal(:,:,:),
     *  ext_smo0to3(:,:,:,:),
     *  ext_smo_bal(:,:,:),
     *  b2stel_smq_ion0to3(:,:,:,:),
     *  b2stel_smq_ion_bal(:,:,:),
     *  b2stel_smq_rec0to3(:,:,:,:),
     *  b2stel_smq_rec_bal(:,:,:),
     *  b2stcx_smq0to3(:,:,:,:),
     *  b2stcx_smq_bal(:,:,:),
     *  b2srsm_smo0to3(:,:,:,:),
     *  b2srsm_smo_bal(:,:,:),
     *  b2srdt_smo0to3(:,:,:,:),
     *  b2srdt_smo_bal(:,:,:),
     *  b2srst_smo0to3(:,:,:,:),
     *  b2srst_smo_bal(:,:,:),
     *  b2sifr_smoch0to3(:,:,:,:),
     *  b2sifr_smoch_bal(:,:,:),
     *  b2sifr_smotf0to3(:,:,:,:),
     *  b2sifr_smotf_bal(:,:,:),
     *  b2siav_smovh0to3(:,:,:,:),
     *  b2siav_smovh_bal(:,:,:),
     *  b2siav_smovv0to3(:,:,:,:),
     *  b2siav_smovv_bal(:,:,:),
     *  b2sicf_smo0to3(:,:,:,:),
     *  b2sicf_smo_bal(:,:,:),
     *  b2nxdv_smo0to3(:,:,:,:),
     *  b2nxdv_smo_bal(:,:,:),
     *  b2sigp_smogp0to3(:,:,:,:),
     *  b2sigp_smogp_bal(:,:,:),
     *  b2sigp_pstat0to3(:,:,:,:),
     *  b2sigp_pstat_bal(:,:,:)
cdjm Jan2017 Arrays for electron heat balance:
      real (kind=R8), allocatable, save ::
     *  b2stbr_she0to3(:,:,:),
     *  b2stbr_she_bal(:,:),
     *  b2stbc_she0to3(:,:,:),
     *  b2stbc_she_bal(:,:),
     *  eirene_mc_eael_she0to3(:,:,:,:),
     *  eirene_mc_emel_she0to3(:,:,:,:),
     *  eirene_mc_eiel_she0to3(:,:,:,:),
     *  eirene_mc_epel_she0to3(:,:,:,:),
     *  eirene_mc_eael_she(:,:,:),
     *  eirene_mc_emel_she(:,:,:),
     *  eirene_mc_eiel_she(:,:,:),
     *  eirene_mc_epel_she(:,:,:),
     *  b2stbm_she0to3(:,:,:),
     *  b2stbm_she_bal(:,:),
     *  ext_she0to3(:,:,:),
     *  ext_she_bal(:,:),
     *  b2stel_she0to3(:,:,:),
     *  b2stel_she_bal(:,:),
     *  b2srsm_she0to3(:,:,:),
     *  b2srsm_she_bal(:,:),
     *  b2srdt_she0to3(:,:,:),
     *  b2srdt_she_bal(:,:),
     *  b2srst_she0to3(:,:,:),
     *  b2srst_she_bal(:,:),
     *  b2sihs_diae0to3(:,:,:),
     *  b2sihs_diae_bal(:,:),
     *  b2sihs_divue0to3(:,:,:),
     *  b2sihs_divue_bal(:,:),
     *  b2sihs_exbe0to3(:,:,:),
     *  b2sihs_exbe_bal(:,:),
     *  b2sihs_joule0to3(:,:,:),
     *  b2sihs_joule_bal(:,:),
     *  b2npht_shei_bal(:,:) !djm equipartition is not linearised
cdjm Jan2017 Arrays for ion heat balance:
      real (kind=R8), allocatable, save ::
     *  b2stbr_shi0to3(:,:,:),
     *  b2stbr_shi_bal(:,:),
     *  b2stbc_shi0to3(:,:,:),
     *  b2stbc_shi_bal(:,:),
     *  eirene_mc_eapl_shi0to3(:,:,:,:),
     *  eirene_mc_empl_shi0to3(:,:,:,:),
     *  eirene_mc_eipl_shi0to3(:,:,:,:),
     *  eirene_mc_eppl_shi0to3(:,:,:,:),
     *  eirene_mc_eapl_shi(:,:,:),
     *  eirene_mc_empl_shi(:,:,:),
     *  eirene_mc_eipl_shi(:,:,:),
     *  eirene_mc_eppl_shi(:,:,:),
     *  b2stbm_shi0to3(:,:,:),
     *  b2stbm_shi_bal(:,:),
     *  ext_shi0to3(:,:,:),
     *  ext_shi_bal(:,:),
     *  b2stel_shi_ion0to3(:,:,:),
     *  b2stel_shi_ion_bal(:,:),
     *  b2stel_shi_rec0to3(:,:,:),
     *  b2stel_shi_rec_bal(:,:),
     *  b2stcx_shi0to3(:,:,:),
     *  b2stcx_shi_bal(:,:),
     *  b2srsm_shi0to3(:,:,:),
     *  b2srsm_shi_bal(:,:),
     *  b2srdt_shi0to3(:,:,:),
     *  b2srdt_shi_bal(:,:),
     *  b2srst_shi0to3(:,:,:),
     *  b2srst_shi_bal(:,:),
     *  b2sihs_diaa0to3(:,:,:),
     *  b2sihs_diaa_bal(:,:),
     *  b2sihs_divua0to3(:,:,:),
     *  b2sihs_divua_bal(:,:),
     *  b2sihs_exba0to3(:,:,:),
     *  b2sihs_exba_bal(:,:),
     *  b2sihs_visa0to3(:,:,:),
     *  b2sihs_visa_bal(:,:),
     *  b2sihs_fraa0to3(:,:,:),
     *  b2sihs_fraa_bal(:,:)

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

cdjm Jan2017 Arrays for particle balance:
      allocate(b2stbr_sna0to1(-1:nxd,-1:nyd,0:1,0:nsd-1))
      allocate(b2stbr_sna_bal(-1:nxd,-1:nyd,0:nsd-1))
      allocate(b2stbc_sna0to1(-1:nxd,-1:nyd,0:1,0:nsd-1))
      allocate(b2stbc_sna_bal(-1:nxd,-1:nyd,0:nsd-1))
      allocate(b2stbm_sna0to1(-1:nxd,-1:nyd,0:1,0:nsd-1))
      allocate(b2stbm_sna_bal(-1:nxd,-1:nyd,0:nsd-1))
      allocate(ext_sna0to1(-1:nxd,-1:nyd,0:1,0:nsd-1))
      allocate(ext_sna_bal(-1:nxd,-1:nyd,0:nsd-1))
      allocate(b2stel_sna_ion0to1(-1:nxd,-1:nyd,0:1,0:nsd-1))
      allocate(b2stel_sna_ion_bal(-1:nxd,-1:nyd,0:nsd-1))
      allocate(b2stel_sna_rec0to1(-1:nxd,-1:nyd,0:1,0:nsd-1))
      allocate(b2stel_sna_rec_bal(-1:nxd,-1:nyd,0:nsd-1))
      allocate(b2stcx_sna0to1(-1:nxd,-1:nyd,0:1,0:nsd-1))
      allocate(b2stcx_sna_bal(-1:nxd,-1:nyd,0:nsd-1))
      allocate(b2srsm_sna0to1(-1:nxd,-1:nyd,0:1,0:nsd-1))
      allocate(b2srsm_sna_bal(-1:nxd,-1:nyd,0:nsd-1))
      allocate(b2srdt_sna0to1(-1:nxd,-1:nyd,0:1,0:nsd-1))
      allocate(b2srdt_sna_bal(-1:nxd,-1:nyd,0:nsd-1))
      allocate(b2srst_sna0to1(-1:nxd,-1:nyd,0:1,0:nsd-1))
      allocate(b2srst_sna_bal(-1:nxd,-1:nyd,0:nsd-1))
      b2stbr_sna0to1=0.0_R8
      b2stbr_sna_bal=0.0_R8
      b2stbc_sna0to1=0.0_R8
      b2stbc_sna_bal=0.0_R8
      b2stbm_sna0to1=0.0_R8
      b2stbm_sna_bal=0.0_R8
      ext_sna0to1=0.0_R8
      ext_sna_bal=0.0_R8
      b2stel_sna_ion0to1=0.0_R8
      b2stel_sna_ion_bal=0.0_R8
      b2stel_sna_rec0to1=0.0_R8
      b2stel_sna_rec_bal=0.0_R8
      b2stcx_sna0to1=0.0_R8
      b2stcx_sna_bal=0.0_R8
      b2srsm_sna0to1=0.0_R8
      b2srsm_sna_bal=0.0_R8
      b2srdt_sna0to1=0.0_R8
      b2srdt_sna_bal=0.0_R8
      b2srst_sna0to1=0.0_R8
      b2srst_sna_bal=0.0_R8
cdjm Arrays for momentum balance:
      allocate(b2stbr_smo0to3(-1:nxd,-1:nyd,0:3,0:nsd-1))
      allocate(b2stbr_smo_bal(-1:nxd,-1:nyd,0:nsd-1))
      allocate(b2stbc_smo0to3(-1:nxd,-1:nyd,0:3,0:nsd-1))
      allocate(b2stbc_smo_bal(-1:nxd,-1:nyd,0:nsd-1))      
      allocate(b2stbm_smo0to3(-1:nxd,-1:nyd,0:3,0:nsd-1))
      allocate(b2stbm_smo_bal(-1:nxd,-1:nyd,0:nsd-1))
      allocate(ext_smo0to3(-1:nxd,-1:nyd,0:3,0:nsd-1))
      allocate(ext_smo_bal(-1:nxd,-1:nyd,0:nsd-1))
      allocate(b2stel_smq_ion0to3(-1:nxd,-1:nyd,0:3,0:nsd-1))
      allocate(b2stel_smq_ion_bal(-1:nxd,-1:nyd,0:nsd-1))
      allocate(b2stel_smq_rec0to3(-1:nxd,-1:nyd,0:3,0:nsd-1))
      allocate(b2stel_smq_rec_bal(-1:nxd,-1:nyd,0:nsd-1))
      allocate(b2stcx_smq0to3(-1:nxd,-1:nyd,0:3,0:nsd-1))
      allocate(b2stcx_smq_bal(-1:nxd,-1:nyd,0:nsd-1))
      allocate(b2srsm_smo0to3(-1:nxd,-1:nyd,0:3,0:nsd-1))
      allocate(b2srsm_smo_bal(-1:nxd,-1:nyd,0:nsd-1))
      allocate(b2srdt_smo0to3(-1:nxd,-1:nyd,0:3,0:nsd-1))
      allocate(b2srdt_smo_bal(-1:nxd,-1:nyd,0:nsd-1))
      allocate(b2srst_smo0to3(-1:nxd,-1:nyd,0:3,0:nsd-1))
      allocate(b2srst_smo_bal(-1:nxd,-1:nyd,0:nsd-1))
      allocate(b2sifr_smoch0to3(-1:nxd,-1:nyd,0:3,0:nsd-1))
      allocate(b2sifr_smoch_bal(-1:nxd,-1:nyd,0:nsd-1))
      allocate(b2sifr_smotf0to3(-1:nxd,-1:nyd,0:3,0:nsd-1))
      allocate(b2sifr_smotf_bal(-1:nxd,-1:nyd,0:nsd-1))
      allocate(b2siav_smovh0to3(-1:nxd,-1:nyd,0:3,0:nsd-1))
      allocate(b2siav_smovh_bal(-1:nxd,-1:nyd,0:nsd-1))
      allocate(b2siav_smovv0to3(-1:nxd,-1:nyd,0:3,0:nsd-1))
      allocate(b2siav_smovv_bal(-1:nxd,-1:nyd,0:nsd-1))
      allocate(b2sicf_smo0to3(-1:nxd,-1:nyd,0:3,0:nsd-1))
      allocate(b2sicf_smo_bal(-1:nxd,-1:nyd,0:nsd-1))
      allocate(b2nxdv_smo0to3(-1:nxd,-1:nyd,0:3,0:nsd-1))
      allocate(b2nxdv_smo_bal(-1:nxd,-1:nyd,0:nsd-1))
      allocate(b2sigp_smogp0to3(-1:nxd,-1:nyd,0:3,0:nsd-1))
      allocate(b2sigp_smogp_bal(-1:nxd,-1:nyd,0:nsd-1))
      allocate(b2sigp_pstat0to3(-1:nxd,-1:nyd,0:3,0:nsd-1))
      allocate(b2sigp_pstat_bal(-1:nxd,-1:nyd,0:nsd-1))
      b2stbr_smo0to3=0.0_R8
      b2stbr_smo_bal=0.0_R8
      b2stbc_smo0to3=0.0_R8
      b2stbc_smo_bal=0.0_R8
      b2stbm_smo0to3=0.0_R8
      b2stbm_smo_bal=0.0_R8
      ext_smo0to3=0.0_R8
      ext_smo_bal=0.0_R8
      b2stel_smq_ion0to3=0.0_R8
      b2stel_smq_ion_bal=0.0_R8
      b2stel_smq_rec0to3=0.0_R8
      b2stel_smq_rec_bal=0.0_R8
      b2stcx_smq0to3=0.0_R8
      b2stcx_smq_bal=0.0_R8
      b2srsm_smo0to3=0.0_R8
      b2srsm_smo_bal=0.0_R8
      b2srdt_smo0to3=0.0_R8
      b2srdt_smo_bal=0.0_R8
      b2srst_smo0to3=0.0_R8
      b2srst_smo_bal=0.0_R8
      b2sifr_smoch0to3=0.0_R8
      b2sifr_smoch_bal=0.0_R8
      b2sifr_smotf0to3=0.0_R8
      b2sifr_smotf_bal=0.0_R8
      b2siav_smovh0to3=0.0_R8
      b2siav_smovh_bal=0.0_R8
      b2siav_smovv0to3=0.0_R8
      b2siav_smovv_bal=0.0_R8
      b2sicf_smo0to3=0.0_R8
      b2sicf_smo_bal=0.0_R8
      b2nxdv_smo0to3=0.0_R8
      b2nxdv_smo_bal=0.0_R8
      b2sigp_smogp0to3=0.0_R8
      b2sigp_smogp_bal=0.0_R8
      b2sigp_pstat0to3=0.0_R8
      b2sigp_pstat_bal=0.0_R8
cdjm Arrays for electron heat balance:
      allocate(b2stbr_she0to3(-1:nxd,-1:nyd,0:3))
      allocate(b2stbr_she_bal(-1:nxd,-1:nyd))
      allocate(b2stbc_she0to3(-1:nxd,-1:nyd,0:3))
      allocate(b2stbc_she_bal(-1:nxd,-1:nyd))
      allocate(b2stbm_she0to3(-1:nxd,-1:nyd,0:3))
      allocate(b2stbm_she_bal(-1:nxd,-1:nyd))
      allocate(ext_she0to3(-1:nxd,-1:nyd,0:3))
      allocate(ext_she_bal(-1:nxd,-1:nyd))
      allocate(b2stel_she0to3(-1:nxd,-1:nyd,0:3))
      allocate(b2stel_she_bal(-1:nxd,-1:nyd))
      allocate(b2srsm_she0to3(-1:nxd,-1:nyd,0:3))
      allocate(b2srsm_she_bal(-1:nxd,-1:nyd))
      allocate(b2srdt_she0to3(-1:nxd,-1:nyd,0:3))
      allocate(b2srdt_she_bal(-1:nxd,-1:nyd))
      allocate(b2srst_she0to3(-1:nxd,-1:nyd,0:3))
      allocate(b2srst_she_bal(-1:nxd,-1:nyd))
      allocate(b2sihs_diae0to3(-1:nxd,-1:nyd,0:3))
      allocate(b2sihs_diae_bal(-1:nxd,-1:nyd))
      allocate(b2sihs_divue0to3(-1:nxd,-1:nyd,0:3))
      allocate(b2sihs_divue_bal(-1:nxd,-1:nyd))
      allocate(b2sihs_exbe0to3(-1:nxd,-1:nyd,0:3))
      allocate(b2sihs_exbe_bal(-1:nxd,-1:nyd))
      allocate(b2sihs_joule0to3(-1:nxd,-1:nyd,0:3))
      allocate(b2sihs_joule_bal(-1:nxd,-1:nyd))
      allocate(b2npht_shei_bal(-1:nxd,-1:nyd))
      b2stbr_she0to3=0.0_R8
      b2stbr_she_bal=0.0_R8
      b2stbc_she0to3=0.0_R8
      b2stbc_she_bal=0.0_R8
      b2stbm_she0to3=0.0_R8
      b2stbm_she_bal=0.0_R8
      ext_she0to3=0.0_R8
      ext_she_bal=0.0_R8
      b2stel_she0to3=0.0_R8
      b2stel_she_bal=0.0_R8
      b2srsm_she0to3=0.0_R8
      b2srsm_she_bal=0.0_R8
      b2srdt_she0to3=0.0_R8
      b2srdt_she_bal=0.0_R8
      b2srst_she0to3=0.0_R8
      b2srst_she_bal=0.0_R8
      b2sihs_diae0to3=0.0_R8
      b2sihs_diae_bal=0.0_R8
      b2sihs_divue0to3=0.0_R8
      b2sihs_divue_bal=0.0_R8
      b2sihs_exbe0to3=0.0_R8
      b2sihs_exbe_bal=0.0_R8
      b2sihs_joule0to3=0.0_R8
      b2sihs_joule_bal=0.0_R8
      b2npht_shei_bal=0.0_R8
cdjm Arrays for ion heat balance:
      allocate(b2stbr_shi0to3(-1:nxd,-1:nyd,0:3))
      allocate(b2stbr_shi_bal(-1:nxd,-1:nyd))
      allocate(b2stbc_shi0to3(-1:nxd,-1:nyd,0:3))
      allocate(b2stbc_shi_bal(-1:nxd,-1:nyd))
      allocate(b2stbm_shi0to3(-1:nxd,-1:nyd,0:3))
      allocate(b2stbm_shi_bal(-1:nxd,-1:nyd))
      allocate(ext_shi0to3(-1:nxd,-1:nyd,0:3))
      allocate(ext_shi_bal(-1:nxd,-1:nyd))
      allocate(b2stel_shi_ion0to3(-1:nxd,-1:nyd,0:3))
      allocate(b2stel_shi_ion_bal(-1:nxd,-1:nyd))
      allocate(b2stel_shi_rec0to3(-1:nxd,-1:nyd,0:3))
      allocate(b2stel_shi_rec_bal(-1:nxd,-1:nyd))
      allocate(b2stcx_shi0to3(-1:nxd,-1:nyd,0:3))
      allocate(b2stcx_shi_bal(-1:nxd,-1:nyd))
      allocate(b2srsm_shi0to3(-1:nxd,-1:nyd,0:3))
      allocate(b2srsm_shi_bal(-1:nxd,-1:nyd))
      allocate(b2srdt_shi0to3(-1:nxd,-1:nyd,0:3))
      allocate(b2srdt_shi_bal(-1:nxd,-1:nyd))
      allocate(b2srst_shi0to3(-1:nxd,-1:nyd,0:3))
      allocate(b2srst_shi_bal(-1:nxd,-1:nyd))
      allocate(b2sihs_diaa0to3(-1:nxd,-1:nyd,0:3))
      allocate(b2sihs_diaa_bal(-1:nxd,-1:nyd))
      allocate(b2sihs_divua0to3(-1:nxd,-1:nyd,0:3))
      allocate(b2sihs_divua_bal(-1:nxd,-1:nyd))
      allocate(b2sihs_exba0to3(-1:nxd,-1:nyd,0:3))
      allocate(b2sihs_exba_bal(-1:nxd,-1:nyd))
      allocate(b2sihs_visa0to3(-1:nxd,-1:nyd,0:3))
      allocate(b2sihs_visa_bal(-1:nxd,-1:nyd))
      allocate(b2sihs_fraa0to3(-1:nxd,-1:nyd,0:3))
      allocate(b2sihs_fraa_bal(-1:nxd,-1:nyd))
      b2stbr_shi0to3=0.0_R8
      b2stbr_shi_bal=0.0_R8
      b2stbc_shi0to3=0.0_R8
      b2stbc_shi_bal=0.0_R8
      b2stbm_shi0to3=0.0_R8
      b2stbm_shi_bal=0.0_R8
      ext_shi0to3=0.0_R8
      ext_shi_bal=0.0_R8
      b2stel_shi_ion0to3=0.0_R8
      b2stel_shi_ion_bal=0.0_R8
      b2stel_shi_rec0to3=0.0_R8
      b2stel_shi_rec_bal=0.0_R8
      b2stcx_shi0to3=0.0_R8
      b2stcx_shi_bal=0.0_R8
      b2srsm_shi0to3=0.0_R8
      b2srsm_shi_bal=0.0_R8
      b2srdt_shi0to3=0.0_R8
      b2srdt_shi_bal=0.0_R8
      b2srst_shi0to3=0.0_R8
      b2srst_shi_bal=0.0_R8
      b2sihs_diaa0to3=0.0_R8
      b2sihs_diaa_bal=0.0_R8
      b2sihs_divua0to3=0.0_R8
      b2sihs_divua_bal=0.0_R8
      b2sihs_exba0to3=0.0_R8
      b2sihs_exba_bal=0.0_R8
      b2sihs_visa0to3=0.0_R8
      b2sihs_visa_bal=0.0_R8
      b2sihs_fraa0to3=0.0_R8
      b2sihs_fraa_bal=0.0_R8
      
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

cdjm Jan2017 Arrays for particle balance:
      deallocate(b2stbr_sna0to1)
      deallocate(b2stbr_sna_bal)
      deallocate(b2stbc_sna0to1)
      deallocate(b2stbc_sna_bal)
      deallocate(b2stbm_sna0to1)
      deallocate(b2stbm_sna_bal)
      deallocate(ext_sna0to1)
      deallocate(ext_sna_bal)
      deallocate(b2stel_sna_ion0to1)
      deallocate(b2stel_sna_ion_bal)
      deallocate(b2stel_sna_rec0to1)
      deallocate(b2stel_sna_rec_bal)
      deallocate(b2stcx_sna0to1)
      deallocate(b2stcx_sna_bal)
      deallocate(b2srsm_sna0to1)
      deallocate(b2srsm_sna_bal)
      deallocate(b2srdt_sna0to1)
      deallocate(b2srdt_sna_bal)
      deallocate(b2srst_sna0to1)
      deallocate(b2srst_sna_bal)
cdjm Arrays for momentum balance:
      deallocate(b2stbr_smo0to3)
      deallocate(b2stbr_smo_bal)
      deallocate(b2stbc_smo0to3)
      deallocate(b2stbc_smo_bal)      
      deallocate(b2stbm_smo0to3)
      deallocate(b2stbm_smo_bal)
      deallocate(ext_smo0to3)
      deallocate(ext_smo_bal)
      deallocate(b2stel_smq_ion0to3)
      deallocate(b2stel_smq_ion_bal)
      deallocate(b2stel_smq_rec0to3)
      deallocate(b2stel_smq_rec_bal)
      deallocate(b2stcx_smq0to3)
      deallocate(b2stcx_smq_bal)
      deallocate(b2srsm_smo0to3)
      deallocate(b2srsm_smo_bal)
      deallocate(b2srdt_smo0to3)
      deallocate(b2srdt_smo_bal)
      deallocate(b2srst_smo0to3)
      deallocate(b2srst_smo_bal)
      deallocate(b2sifr_smoch0to3)
      deallocate(b2sifr_smoch_bal)
      deallocate(b2sifr_smotf0to3)
      deallocate(b2sifr_smotf_bal)
      deallocate(b2siav_smovh0to3)
      deallocate(b2siav_smovh_bal)
      deallocate(b2siav_smovv0to3)
      deallocate(b2siav_smovv_bal)
      deallocate(b2sicf_smo0to3)
      deallocate(b2sicf_smo_bal)
      deallocate(b2nxdv_smo0to3)
      deallocate(b2nxdv_smo_bal)
      deallocate(b2sigp_smogp0to3)
      deallocate(b2sigp_smogp_bal)
      deallocate(b2sigp_pstat0to3)
      deallocate(b2sigp_pstat_bal)
cdjm Arrays for electron heat balance:
      deallocate(b2stbr_she0to3)
      deallocate(b2stbr_she_bal)
      deallocate(b2stbc_she0to3)
      deallocate(b2stbc_she_bal)
      deallocate(b2stbm_she0to3)
      deallocate(b2stbm_she_bal)
      deallocate(ext_she0to3)
      deallocate(ext_she_bal)
      deallocate(b2stel_she0to3)
      deallocate(b2stel_she_bal)
      deallocate(b2srsm_she0to3)
      deallocate(b2srsm_she_bal)
      deallocate(b2srdt_she0to3)
      deallocate(b2srdt_she_bal)
      deallocate(b2srst_she0to3)
      deallocate(b2srst_she_bal)
      deallocate(b2sihs_diae0to3)
      deallocate(b2sihs_diae_bal)
      deallocate(b2sihs_divue0to3)
      deallocate(b2sihs_divue_bal)
      deallocate(b2sihs_exbe0to3)
      deallocate(b2sihs_exbe_bal)
      deallocate(b2sihs_joule0to3)
      deallocate(b2sihs_joule_bal)
      deallocate(b2npht_shei_bal)
cdjm Arrays for ion heat balance:
      deallocate(b2stbr_shi0to3)
      deallocate(b2stbr_shi_bal)
      deallocate(b2stbc_shi0to3)
      deallocate(b2stbc_shi_bal)
      deallocate(b2stbm_shi0to3)
      deallocate(b2stbm_shi_bal)
      deallocate(ext_shi0to3)
      deallocate(ext_shi_bal)
      deallocate(b2stel_shi_ion0to3)
      deallocate(b2stel_shi_ion_bal)
      deallocate(b2stel_shi_rec0to3)
      deallocate(b2stel_shi_rec_bal)
      deallocate(b2stcx_shi0to3)
      deallocate(b2stcx_shi_bal)
      deallocate(b2srsm_shi0to3)
      deallocate(b2srsm_shi_bal)
      deallocate(b2srdt_shi0to3)
      deallocate(b2srdt_shi_bal)
      deallocate(b2srst_shi0to3)
      deallocate(b2srst_shi_bal)
      deallocate(b2sihs_diaa0to3)
      deallocate(b2sihs_diaa_bal)
      deallocate(b2sihs_divua0to3)
      deallocate(b2sihs_divua_bal)
      deallocate(b2sihs_exba0to3)
      deallocate(b2sihs_exba_bal)
      deallocate(b2sihs_visa0to3)
      deallocate(b2sihs_visa_bal)
      deallocate(b2sihs_fraa0to3)
      deallocate(b2sihs_fraa_bal)

      return
      end subroutine dealloc_b2mod_sources

      subroutine alloc_b2mod_eirene_sources(nxd,nyd,nsd,nstrad)
      implicit none
      integer nxd,nyd,nsd,nstrad

cdjm Jan2017 Eirene arrays for particle balance:
      allocate(
     & eirene_mc_papl_sna0to1(-1:nxd,-1:nyd,0:1,0:nsd-1,1:nstrad))
      allocate(
     & eirene_mc_pmpl_sna0to1(-1:nxd,-1:nyd,0:1,0:nsd-1,1:nstrad))
      allocate(
     & eirene_mc_pipl_sna0to1(-1:nxd,-1:nyd,0:1,0:nsd-1,1:nstrad))
      allocate(
     & eirene_mc_pppl_sna0to1(-1:nxd,-1:nyd,0:1,0:nsd-1,1:nstrad))
      allocate(
     & eirene_mc_core_sna0to1(-1:nxd,-1:nyd,0:1,0:nsd-1))
      allocate(eirene_mc_papl_sna(-1:nxd,-1:nyd,0:nsd-1,1:nstrad))
      allocate(eirene_mc_pmpl_sna(-1:nxd,-1:nyd,0:nsd-1,1:nstrad))
      allocate(eirene_mc_pipl_sna(-1:nxd,-1:nyd,0:nsd-1,1:nstrad))
      allocate(eirene_mc_pppl_sna(-1:nxd,-1:nyd,0:nsd-1,1:nstrad))
      allocate(eirene_mc_core_sna(-1:nxd,-1:nyd,0:nsd-1))
      eirene_mc_papl_sna0to1=0.0_R8
      eirene_mc_pmpl_sna0to1=0.0_R8
      eirene_mc_pipl_sna0to1=0.0_R8
      eirene_mc_pppl_sna0to1=0.0_R8
      eirene_mc_core_sna0to1=0.0_R8
      eirene_mc_papl_sna=0.0_R8
      eirene_mc_pmpl_sna=0.0_R8
      eirene_mc_pipl_sna=0.0_R8
      eirene_mc_pppl_sna=0.0_R8
      eirene_mc_core_sna=0.0_R8
cdjm Jan2017 Eirene arrays for momentum balance:
      allocate(
     & eirene_mc_mapl_smo0to3(-1:nxd,-1:nyd,0:3,0:nsd-1,1:nstrad))
      allocate(
     & eirene_mc_mmpl_smo0to3(-1:nxd,-1:nyd,0:3,0:nsd-1,1:nstrad))
      allocate(
     & eirene_mc_mipl_smo0to3(-1:nxd,-1:nyd,0:3,0:nsd-1,1:nstrad))
      allocate(
     & eirene_mc_cppv_smo0to3(-1:nxd,-1:nyd,0:3,0:nsd-1,1:nstrad))
      allocate(eirene_mc_mapl_smo(-1:nxd,-1:nyd,0:nsd-1,1:nstrad))
      allocate(eirene_mc_mmpl_smo(-1:nxd,-1:nyd,0:nsd-1,1:nstrad))
      allocate(eirene_mc_mipl_smo(-1:nxd,-1:nyd,0:nsd-1,1:nstrad))
      allocate(eirene_mc_cppv_smo(-1:nxd,-1:nyd,0:nsd-1,1:nstrad))
      eirene_mc_mapl_smo0to3=0.0_R8
      eirene_mc_mmpl_smo0to3=0.0_R8
      eirene_mc_mipl_smo0to3=0.0_R8
      eirene_mc_cppv_smo0to3=0.0_R8
      eirene_mc_mapl_smo=0.0_R8
      eirene_mc_mmpl_smo=0.0_R8
      eirene_mc_mipl_smo=0.0_R8
      eirene_mc_cppv_smo=0.0_R8
cdjm Jan2017 Eirene arrays for electron heat balance:
      allocate(eirene_mc_eael_she0to3(-1:nxd,-1:nyd,0:3,1:nstrad))
      allocate(eirene_mc_emel_she0to3(-1:nxd,-1:nyd,0:3,1:nstrad))
      allocate(eirene_mc_eiel_she0to3(-1:nxd,-1:nyd,0:3,1:nstrad))
      allocate(eirene_mc_epel_she0to3(-1:nxd,-1:nyd,0:3,1:nstrad))
      allocate(eirene_mc_eael_she(-1:nxd,-1:nyd,1:nstrad))
      allocate(eirene_mc_emel_she(-1:nxd,-1:nyd,1:nstrad))
      allocate(eirene_mc_eiel_she(-1:nxd,-1:nyd,1:nstrad))
      allocate(eirene_mc_epel_she(-1:nxd,-1:nyd,1:nstrad))
      eirene_mc_eael_she0to3=0.0_R8
      eirene_mc_emel_she0to3=0.0_R8
      eirene_mc_eiel_she0to3=0.0_R8
      eirene_mc_epel_she0to3=0.0_R8
      eirene_mc_eael_she=0.0_R8
      eirene_mc_emel_she=0.0_R8
      eirene_mc_eiel_she=0.0_R8
      eirene_mc_epel_she=0.0_R8      
cdjm Jan2017 Eirene arrays for ion heat balance:
      allocate(eirene_mc_eapl_shi0to3(-1:nxd,-1:nyd,0:3,1:nstrad))
      allocate(eirene_mc_empl_shi0to3(-1:nxd,-1:nyd,0:3,1:nstrad))
      allocate(eirene_mc_eipl_shi0to3(-1:nxd,-1:nyd,0:3,1:nstrad))
      allocate(eirene_mc_eppl_shi0to3(-1:nxd,-1:nyd,0:3,1:nstrad))
      allocate(eirene_mc_eapl_shi(-1:nxd,-1:nyd,1:nstrad))
      allocate(eirene_mc_empl_shi(-1:nxd,-1:nyd,1:nstrad))
      allocate(eirene_mc_eipl_shi(-1:nxd,-1:nyd,1:nstrad))
      allocate(eirene_mc_eppl_shi(-1:nxd,-1:nyd,1:nstrad))
      eirene_mc_eapl_shi0to3=0.0_R8
      eirene_mc_empl_shi0to3=0.0_R8
      eirene_mc_eipl_shi0to3=0.0_R8
      eirene_mc_eppl_shi0to3=0.0_R8
      eirene_mc_eapl_shi=0.0_R8
      eirene_mc_empl_shi=0.0_R8
      eirene_mc_eipl_shi=0.0_R8
      eirene_mc_eppl_shi=0.0_R8

      return
      end subroutine alloc_b2mod_eirene_sources

      subroutine dealloc_b2mod_eirene_sources
      implicit none

      if(.not.allocated(eirene_mc_papl_sna)) return
cdjm Jan2017 Eirene arrays for particle balance:
      deallocate(eirene_mc_papl_sna0to1)
      deallocate(eirene_mc_pmpl_sna0to1)
      deallocate(eirene_mc_pipl_sna0to1)
      deallocate(eirene_mc_pppl_sna0to1)
      deallocate(eirene_mc_core_sna0to1)
      deallocate(eirene_mc_papl_sna)
      deallocate(eirene_mc_pmpl_sna)
      deallocate(eirene_mc_pipl_sna)
      deallocate(eirene_mc_pppl_sna)
      deallocate(eirene_mc_core_sna)
cdjm Jan2017 Eirene arrays for momentum balance:
      deallocate(eirene_mc_mapl_smo0to3)
      deallocate(eirene_mc_mmpl_smo0to3)
      deallocate(eirene_mc_mipl_smo0to3)
      deallocate(eirene_mc_cppv_smo0to3)
      deallocate(eirene_mc_mapl_smo)
      deallocate(eirene_mc_mmpl_smo)
      deallocate(eirene_mc_mipl_smo)
      deallocate(eirene_mc_cppv_smo)
cdjm Jan2017 Eirene arrays for electron heat balance:
      deallocate(eirene_mc_eael_she0to3)
      deallocate(eirene_mc_emel_she0to3)
      deallocate(eirene_mc_eiel_she0to3)
      deallocate(eirene_mc_epel_she0to3)
      deallocate(eirene_mc_eael_she)
      deallocate(eirene_mc_emel_she)
      deallocate(eirene_mc_eiel_she)
      deallocate(eirene_mc_epel_she)
cdjm Jan2017 Eirene arrays for ion heat balance:      
      deallocate(eirene_mc_eapl_shi0to3)
      deallocate(eirene_mc_empl_shi0to3)
      deallocate(eirene_mc_eipl_shi0to3)
      deallocate(eirene_mc_eppl_shi0to3)
      deallocate(eirene_mc_eapl_shi)
      deallocate(eirene_mc_empl_shi)
      deallocate(eirene_mc_eipl_shi)
      deallocate(eirene_mc_eppl_shi)

      return
      end subroutine dealloc_b2mod_eirene_sources

      subroutine update_balance_sources_particles(nx,ny,isb,nb) !djm Jan2017
      implicit none
*   ..input arguments (unchanged on exit)
      integer nx,ny,isb
      real (kind=R8):: nb(-1:nx,-1:ny)
*   ..local variables
      integer ix,iy,istrai,nstrad

      nstrad = size(eirene_mc_papl_sna,4)
      
      do ix=-1,nx
        do iy=-1,ny
          b2stbr_sna_bal(ix,iy,isb)=b2stbr_sna0to1(ix,iy,0,isb)+
     &                     b2stbr_sna0to1(ix,iy,1,isb)*nb(ix,iy)
          b2stbc_sna_bal(ix,iy,isb)=b2stbc_sna0to1(ix,iy,0,isb)+
     &                     b2stbc_sna0to1(ix,iy,1,isb)*nb(ix,iy)
          if (allocated(eirene_mc_papl_sna)) then
            do istrai=1,nstrad
              eirene_mc_papl_sna(ix,iy,isb,istrai)=
     &         eirene_mc_papl_sna0to1(ix,iy,0,isb,istrai)+
     &         eirene_mc_papl_sna0to1(ix,iy,1,isb,istrai)*nb(ix,iy)
              eirene_mc_pmpl_sna(ix,iy,isb,istrai)=
     &         eirene_mc_pmpl_sna0to1(ix,iy,0,isb,istrai)+
     &         eirene_mc_pmpl_sna0to1(ix,iy,1,isb,istrai)*nb(ix,iy)
              eirene_mc_pipl_sna(ix,iy,isb,istrai)=
     &         eirene_mc_pipl_sna0to1(ix,iy,0,isb,istrai)+
     &         eirene_mc_pipl_sna0to1(ix,iy,1,isb,istrai)*nb(ix,iy)
              eirene_mc_pppl_sna(ix,iy,isb,istrai)=
     &         eirene_mc_pppl_sna0to1(ix,iy,0,isb,istrai)+
     &         eirene_mc_pppl_sna0to1(ix,iy,1,isb,istrai)*nb(ix,iy)
            enddo
            eirene_mc_core_sna(ix,iy,isb)=
     &       eirene_mc_core_sna0to1(ix,iy,0,isb)+
     &       eirene_mc_core_sna0to1(ix,iy,1,isb)*nb(ix,iy)
          endif
          b2stbm_sna_bal(ix,iy,isb)=b2stbm_sna0to1(ix,iy,0,isb)+
     &                     b2stbm_sna0to1(ix,iy,1,isb)*nb(ix,iy)
          ext_sna_bal(ix,iy,isb)=ext_sna0to1(ix,iy,0,isb)+
     &                     ext_sna0to1(ix,iy,1,isb)*nb(ix,iy)
          b2stel_sna_ion_bal(ix,iy,isb)=b2stel_sna_ion0to1(ix,iy,0,isb)+
     &                     b2stel_sna_ion0to1(ix,iy,1,isb)*nb(ix,iy)
          b2stel_sna_rec_bal(ix,iy,isb)=b2stel_sna_rec0to1(ix,iy,0,isb)+
     &                     b2stel_sna_rec0to1(ix,iy,1,isb)*nb(ix,iy)
          b2stcx_sna_bal(ix,iy,isb)=b2stcx_sna0to1(ix,iy,0,isb)+
     &                     b2stcx_sna0to1(ix,iy,1,isb)*nb(ix,iy)
          b2srsm_sna_bal(ix,iy,isb)=b2srsm_sna0to1(ix,iy,0,isb)+
     &                     b2srsm_sna0to1(ix,iy,1,isb)*nb(ix,iy)
          b2srdt_sna_bal(ix,iy,isb)=b2srdt_sna0to1(ix,iy,0,isb)+
     &                     b2srdt_sna0to1(ix,iy,1,isb)*nb(ix,iy)
          b2srst_sna_bal(ix,iy,isb)=b2srst_sna0to1(ix,iy,0,isb)+
     &                     b2srst_sna0to1(ix,iy,1,isb)*nb(ix,iy)
        enddo
      enddo
      return
      end subroutine update_balance_sources_particles

      subroutine update_balance_sources_momentum(nx,ny,isb,ub,rob) !djm Jan2017
      use b2mod_constants
     & , only : mp
      use b2mod_b2cmpa
     & , only : am
      implicit none
*   ..input arguments (unchanged on exit)
      integer nx,ny,isb
      real (kind=R8):: ub(-1:nx,-1:ny), rob(-1:nx,-1:ny)
*   ..local variables
      integer ix,iy,istrai,nstrad

      nstrad = size(eirene_mc_papl_sna,4)

      do ix=-1,nx
        do iy=-1,ny
          b2stbr_smo_bal(ix,iy,isb)=b2stbr_smo0to3(ix,iy,0,isb)+
     &            b2stbr_smo0to3(ix,iy,1,isb)*ub(ix,iy)+
     &            b2stbr_smo0to3(ix,iy,2,isb)*rob(ix,iy)+
     &            b2stbr_smo0to3(ix,iy,3,isb)*rob(ix,iy)*ub(ix,iy)
          b2stbc_smo_bal(ix,iy,isb)=b2stbc_smo0to3(ix,iy,0,isb)+
     &            b2stbc_smo0to3(ix,iy,1,isb)*ub(ix,iy)+
     &            b2stbc_smo0to3(ix,iy,2,isb)*rob(ix,iy)+
     &            b2stbc_smo0to3(ix,iy,3,isb)*rob(ix,iy)*ub(ix,iy)
          if (allocated(eirene_mc_mapl_smo)) then
            do istrai=1,nstrad
              eirene_mc_mapl_smo(ix,iy,isb,istrai)=
     &         eirene_mc_mapl_smo0to3(ix,iy,0,isb,istrai)+
     &         eirene_mc_mapl_smo0to3(ix,iy,1,isb,istrai)*ub(ix,iy)+
     &         eirene_mc_mapl_smo0to3(ix,iy,2,isb,istrai)*rob(ix,iy)+
     &         eirene_mc_mapl_smo0to3(ix,iy,3,isb,istrai)*rob(ix,iy)*
     &                                                    ub(ix,iy)
              eirene_mc_mmpl_smo(ix,iy,isb,istrai)=
     &         eirene_mc_mmpl_smo0to3(ix,iy,0,isb,istrai)+
     &         eirene_mc_mmpl_smo0to3(ix,iy,1,isb,istrai)*ub(ix,iy)+
     &         eirene_mc_mmpl_smo0to3(ix,iy,2,isb,istrai)*rob(ix,iy)+
     &         eirene_mc_mmpl_smo0to3(ix,iy,3,isb,istrai)*rob(ix,iy)*
     &                                                    ub(ix,iy)
              eirene_mc_mipl_smo(ix,iy,isb,istrai)=
     &         eirene_mc_mipl_smo0to3(ix,iy,0,isb,istrai)+
     &         eirene_mc_mipl_smo0to3(ix,iy,1,isb,istrai)*ub(ix,iy)+
     &         eirene_mc_mipl_smo0to3(ix,iy,2,isb,istrai)*rob(ix,iy)+
     &         eirene_mc_mipl_smo0to3(ix,iy,3,isb,istrai)*rob(ix,iy)*
     &                                                    ub(ix,iy)
              eirene_mc_cppv_smo(ix,iy,isb,istrai)=
     &         eirene_mc_cppv_smo0to3(ix,iy,0,isb,istrai)+
     &         eirene_mc_cppv_smo0to3(ix,iy,1,isb,istrai)*ub(ix,iy)+
     &         eirene_mc_cppv_smo0to3(ix,iy,2,isb,istrai)*rob(ix,iy)+
     &         eirene_mc_cppv_smo0to3(ix,iy,3,isb,istrai)*rob(ix,iy)*
     &                                                    ub(ix,iy)
            enddo
          endif
          b2stbm_smo_bal(ix,iy,isb)=b2stbm_smo0to3(ix,iy,0,isb)+
     &            b2stbm_smo0to3(ix,iy,1,isb)*ub(ix,iy)+
     &            b2stbm_smo0to3(ix,iy,2,isb)*rob(ix,iy)+
     &            b2stbm_smo0to3(ix,iy,3,isb)*rob(ix,iy)*ub(ix,iy)
          ext_smo_bal(ix,iy,isb)=ext_smo0to3(ix,iy,0,isb)+
     &            ext_smo0to3(ix,iy,1,isb)*ub(ix,iy)+
     &            ext_smo0to3(ix,iy,2,isb)*rob(ix,iy)+
     &            ext_smo0to3(ix,iy,3,isb)*rob(ix,iy)*ub(ix,iy)
          b2stel_smq_ion_bal(ix,iy,isb)=b2stel_smq_ion0to3(ix,iy,0,isb)+
     &            b2stel_smq_ion0to3(ix,iy,1,isb)*ub(ix,iy)+
     &            b2stel_smq_ion0to3(ix,iy,2,isb)*rob(ix,iy)+
     &            b2stel_smq_ion0to3(ix,iy,3,isb)*rob(ix,iy)*ub(ix,iy)
          b2stel_smq_rec_bal(ix,iy,isb)=b2stel_smq_rec0to3(ix,iy,0,isb)+
     &            b2stel_smq_rec0to3(ix,iy,1,isb)*ub(ix,iy)+
     &            b2stel_smq_rec0to3(ix,iy,2,isb)*rob(ix,iy)+
     &            b2stel_smq_rec0to3(ix,iy,3,isb)*rob(ix,iy)*ub(ix,iy)
          b2stcx_smq_bal(ix,iy,isb)=b2stcx_smq0to3(ix,iy,0,isb)+
     &            b2stcx_smq0to3(ix,iy,1,isb)*ub(ix,iy)+
     &            b2stcx_smq0to3(ix,iy,2,isb)*rob(ix,iy)+
     &            b2stcx_smq0to3(ix,iy,3,isb)*rob(ix,iy)*ub(ix,iy)
          b2srsm_smo_bal(ix,iy,isb)=b2srsm_smo0to3(ix,iy,0,isb)+
     &            b2srsm_smo0to3(ix,iy,1,isb)*ub(ix,iy)+
     &            b2srsm_smo0to3(ix,iy,2,isb)*rob(ix,iy)+
     &            b2srsm_smo0to3(ix,iy,3,isb)*rob(ix,iy)*ub(ix,iy)
          b2srdt_smo_bal(ix,iy,isb)=b2srdt_smo0to3(ix,iy,0,isb)+
     &            b2srdt_smo0to3(ix,iy,1,isb)*ub(ix,iy)+
     &            b2srdt_smo0to3(ix,iy,2,isb)*rob(ix,iy)+
     &            b2srdt_smo0to3(ix,iy,3,isb)*rob(ix,iy)*ub(ix,iy)
          b2srst_smo_bal(ix,iy,isb)=b2srst_smo0to3(ix,iy,0,isb)+
     &            b2srst_smo0to3(ix,iy,1,isb)*ub(ix,iy)+
     &            b2srst_smo0to3(ix,iy,2,isb)*rob(ix,iy)+
     &            b2srst_smo0to3(ix,iy,3,isb)*rob(ix,iy)*ub(ix,iy)
          b2sifr_smoch_bal(ix,iy,isb)=b2sifr_smoch0to3(ix,iy,0,isb)+
     &            b2sifr_smoch0to3(ix,iy,1,isb)*ub(ix,iy)+
     &            b2sifr_smoch0to3(ix,iy,2,isb)*rob(ix,iy)+
     &            b2sifr_smoch0to3(ix,iy,3,isb)*rob(ix,iy)*ub(ix,iy)
          b2sifr_smotf_bal(ix,iy,isb)=b2sifr_smotf0to3(ix,iy,0,isb)+
     &            b2sifr_smotf0to3(ix,iy,1,isb)*ub(ix,iy)+
     &            b2sifr_smotf0to3(ix,iy,2,isb)*rob(ix,iy)+
     &            b2sifr_smotf0to3(ix,iy,3,isb)*rob(ix,iy)*ub(ix,iy)
          b2siav_smovh_bal(ix,iy,isb)=b2siav_smovh0to3(ix,iy,0,isb)+
     &            b2siav_smovh0to3(ix,iy,1,isb)*ub(ix,iy)+
     &            b2siav_smovh0to3(ix,iy,2,isb)*rob(ix,iy)+
     &            b2siav_smovh0to3(ix,iy,3,isb)*rob(ix,iy)*ub(ix,iy)
          b2siav_smovv_bal(ix,iy,isb)=b2siav_smovv0to3(ix,iy,0,isb)+
     &            b2siav_smovv0to3(ix,iy,1,isb)*ub(ix,iy)+
     &            b2siav_smovv0to3(ix,iy,2,isb)*rob(ix,iy)+
     &            b2siav_smovv0to3(ix,iy,3,isb)*rob(ix,iy)*ub(ix,iy)
          b2sicf_smo_bal(ix,iy,isb)=b2sicf_smo0to3(ix,iy,0,isb)+
     &            b2sicf_smo0to3(ix,iy,1,isb)*ub(ix,iy)+
     &            b2sicf_smo0to3(ix,iy,2,isb)*rob(ix,iy)+
     &            b2sicf_smo0to3(ix,iy,3,isb)*rob(ix,iy)*ub(ix,iy)
          b2nxdv_smo_bal(ix,iy,isb)=b2nxdv_smo0to3(ix,iy,0,isb)+
     &            b2nxdv_smo0to3(ix,iy,1,isb)*ub(ix,iy)+
     &            b2nxdv_smo0to3(ix,iy,2,isb)*rob(ix,iy)+
     &            b2nxdv_smo0to3(ix,iy,3,isb)*rob(ix,iy)*ub(ix,iy)
          b2sigp_smogp_bal(ix,iy,isb)=b2sigp_smogp0to3(ix,iy,0,isb)+
     &            b2sigp_smogp0to3(ix,iy,1,isb)*ub(ix,iy)+
     &            b2sigp_smogp0to3(ix,iy,2,isb)*rob(ix,iy)+
     &            b2sigp_smogp0to3(ix,iy,3,isb)*rob(ix,iy)*ub(ix,iy)
          b2sigp_pstat_bal(ix,iy,isb)=b2sigp_pstat0to3(ix,iy,0,isb)+
     &            b2sigp_pstat0to3(ix,iy,1,isb)*ub(ix,iy)+
     &            b2sigp_pstat0to3(ix,iy,2,isb)*rob(ix,iy)+
     &            b2sigp_pstat0to3(ix,iy,3,isb)*rob(ix,iy)*ub(ix,iy)
        enddo
      enddo
      return
      end subroutine update_balance_sources_momentum

      subroutine update_balance_sources_heat(nx,ny,ne,ni,te,ti)
      implicit none
*   ..input arguments (unchanged on exit)
      integer nx,ny,isb
      real (kind=R8):: ne(-1:nx,-1:ny), ni(-1:nx,-1:ny),
     &                 te(-1:nx,-1:ny), ti(-1:nx,-1:ny)
*   ..local variables
      integer ix,iy,istrai,nstrad

      nstrad = size(eirene_mc_papl_sna,4)

cdjm Update electron heat sources
      do ix=-1,nx
        do iy=-1,ny
          b2stbr_she_bal(ix,iy)=
     &     b2stbr_she0to3(ix,iy,0)+
     &     b2stbr_she0to3(ix,iy,1)*te(ix,iy)+
     &     b2stbr_she0to3(ix,iy,2)*ne(ix,iy)+
     &     b2stbr_she0to3(ix,iy,3)*ne(ix,iy)*te(ix,iy)
          b2stbc_she_bal(ix,iy)=
     &     b2stbc_she0to3(ix,iy,0)+
     &     b2stbc_she0to3(ix,iy,1)*te(ix,iy)+
     &     b2stbc_she0to3(ix,iy,2)*ne(ix,iy)+
     &     b2stbc_she0to3(ix,iy,3)*ne(ix,iy)*te(ix,iy)
          if (allocated(eirene_mc_eael_she)) then
            do istrai=1,nstrad
              eirene_mc_eael_she(ix,iy,istrai)=
     &         eirene_mc_eael_she0to3(ix,iy,0,istrai)+
     &         eirene_mc_eael_she0to3(ix,iy,1,istrai)*te(ix,iy)+
     &         eirene_mc_eael_she0to3(ix,iy,2,istrai)*ne(ix,iy)+
     &         eirene_mc_eael_she0to3(ix,iy,3,istrai)*ne(ix,iy)*
     &                                                te(ix,iy)
              eirene_mc_emel_she(ix,iy,istrai)=
     &         eirene_mc_emel_she0to3(ix,iy,0,istrai)+
     &         eirene_mc_emel_she0to3(ix,iy,1,istrai)*te(ix,iy)+
     &         eirene_mc_emel_she0to3(ix,iy,2,istrai)*ne(ix,iy)+
     &         eirene_mc_emel_she0to3(ix,iy,3,istrai)*ne(ix,iy)*
     &                                                te(ix,iy)
              eirene_mc_eiel_she(ix,iy,istrai)=
     &         eirene_mc_eiel_she0to3(ix,iy,0,istrai)+
     &         eirene_mc_eiel_she0to3(ix,iy,1,istrai)*te(ix,iy)+
     &         eirene_mc_eiel_she0to3(ix,iy,2,istrai)*ne(ix,iy)+
     &         eirene_mc_eiel_she0to3(ix,iy,3,istrai)*ne(ix,iy)*
     &                                                te(ix,iy)
              eirene_mc_epel_she(ix,iy,istrai)=
     &         eirene_mc_epel_she0to3(ix,iy,0,istrai)+
     &         eirene_mc_epel_she0to3(ix,iy,1,istrai)*te(ix,iy)+
     &         eirene_mc_epel_she0to3(ix,iy,2,istrai)*ne(ix,iy)+
     &         eirene_mc_epel_she0to3(ix,iy,3,istrai)*ne(ix,iy)*
     &                                                te(ix,iy)
            enddo
          endif
          b2stbm_she_bal(ix,iy)=
     &     b2stbm_she0to3(ix,iy,0)+
     &     b2stbm_she0to3(ix,iy,1)*te(ix,iy)+
     &     b2stbm_she0to3(ix,iy,2)*ne(ix,iy)+
     &     b2stbm_she0to3(ix,iy,3)*ne(ix,iy)*te(ix,iy)
          ext_she_bal(ix,iy)=
     &     ext_she0to3(ix,iy,0)+
     &     ext_she0to3(ix,iy,1)*te(ix,iy)+
     &     ext_she0to3(ix,iy,2)*ne(ix,iy)+
     &     ext_she0to3(ix,iy,3)*ne(ix,iy)*te(ix,iy)
          b2stel_she_bal(ix,iy)=
     &     b2stel_she0to3(ix,iy,0)+
     &     b2stel_she0to3(ix,iy,1)*te(ix,iy)+
     &     b2stel_she0to3(ix,iy,2)*ne(ix,iy)+
     &     b2stel_she0to3(ix,iy,3)*ne(ix,iy)*te(ix,iy)
          b2srsm_she_bal(ix,iy)=
     &     b2srsm_she0to3(ix,iy,0)+
     &     b2srsm_she0to3(ix,iy,1)*te(ix,iy)+
     &     b2srsm_she0to3(ix,iy,2)*ne(ix,iy)+
     &     b2srsm_she0to3(ix,iy,3)*ne(ix,iy)*te(ix,iy)
          b2srdt_she_bal(ix,iy)=
     &     b2srdt_she0to3(ix,iy,0)+
     &     b2srdt_she0to3(ix,iy,1)*te(ix,iy)+
     &     b2srdt_she0to3(ix,iy,2)*ne(ix,iy)+
     &     b2srdt_she0to3(ix,iy,3)*ne(ix,iy)*te(ix,iy)
          b2srst_she_bal(ix,iy)=
     &     b2srst_she0to3(ix,iy,0)+
     &     b2srst_she0to3(ix,iy,1)*te(ix,iy)+
     &     b2srst_she0to3(ix,iy,2)*ne(ix,iy)+
     &     b2srst_she0to3(ix,iy,3)*ne(ix,iy)*te(ix,iy)
          b2sihs_diae_bal(ix,iy)=
     &     b2sihs_diae0to3(ix,iy,0)+
     &     b2sihs_diae0to3(ix,iy,1)*te(ix,iy)+
     &     b2sihs_diae0to3(ix,iy,2)*ne(ix,iy)+
     &     b2sihs_diae0to3(ix,iy,3)*ne(ix,iy)*te(ix,iy)
          b2sihs_divue_bal(ix,iy)=
     &     b2sihs_divue0to3(ix,iy,0)+
     &     b2sihs_divue0to3(ix,iy,1)*te(ix,iy)+
     &     b2sihs_divue0to3(ix,iy,2)*ne(ix,iy)+
     &     b2sihs_divue0to3(ix,iy,3)*ne(ix,iy)*te(ix,iy)
          b2sihs_exbe_bal(ix,iy)=
     &     b2sihs_exbe0to3(ix,iy,0)+
     &     b2sihs_exbe0to3(ix,iy,1)*te(ix,iy)+
     &     b2sihs_exbe0to3(ix,iy,2)*ne(ix,iy)+
     &     b2sihs_exbe0to3(ix,iy,3)*ne(ix,iy)*te(ix,iy)
          b2sihs_joule_bal(ix,iy)=
     &     b2sihs_joule0to3(ix,iy,0)+
     &     b2sihs_joule0to3(ix,iy,1)*te(ix,iy)+
     &     b2sihs_joule0to3(ix,iy,2)*ne(ix,iy)+
     &     b2sihs_joule0to3(ix,iy,3)*ne(ix,iy)*te(ix,iy)
        enddo
      enddo
cdjm Update ion heat sources
      do ix=-1,nx
        do iy=-1,ny
          b2stbr_shi_bal(ix,iy)=
     &     b2stbr_shi0to3(ix,iy,0)+
     &     b2stbr_shi0to3(ix,iy,1)*ti(ix,iy)+
     &     b2stbr_shi0to3(ix,iy,2)*ni(ix,iy)+
     &     b2stbr_shi0to3(ix,iy,3)*ni(ix,iy)*ti(ix,iy)
          b2stbc_shi_bal(ix,iy)=
     &     b2stbc_shi0to3(ix,iy,0)+
     &     b2stbc_shi0to3(ix,iy,1)*ti(ix,iy)+
     &     b2stbc_shi0to3(ix,iy,2)*ni(ix,iy)+
     &     b2stbc_shi0to3(ix,iy,3)*ni(ix,iy)*ti(ix,iy)
          if (allocated(eirene_mc_eapl_shi)) then
            do istrai=1,nstrad
              eirene_mc_eapl_shi(ix,iy,istrai)=
     &         eirene_mc_eapl_shi0to3(ix,iy,0,istrai)+
     &         eirene_mc_eapl_shi0to3(ix,iy,1,istrai)*ti(ix,iy)+
     &         eirene_mc_eapl_shi0to3(ix,iy,2,istrai)*ni(ix,iy)+
     &         eirene_mc_eapl_shi0to3(ix,iy,3,istrai)*ni(ix,iy)*
     &                                                ti(ix,iy)
              eirene_mc_empl_shi(ix,iy,istrai)=
     &         eirene_mc_empl_shi0to3(ix,iy,0,istrai)+
     &         eirene_mc_empl_shi0to3(ix,iy,1,istrai)*ti(ix,iy)+
     &         eirene_mc_empl_shi0to3(ix,iy,2,istrai)*ni(ix,iy)+
     &         eirene_mc_empl_shi0to3(ix,iy,3,istrai)*ni(ix,iy)*
     &                                                ti(ix,iy)
              eirene_mc_eipl_shi(ix,iy,istrai)=
     &         eirene_mc_eipl_shi0to3(ix,iy,0,istrai)+
     &         eirene_mc_eipl_shi0to3(ix,iy,1,istrai)*ti(ix,iy)+
     &         eirene_mc_eipl_shi0to3(ix,iy,2,istrai)*ni(ix,iy)+
     &         eirene_mc_eipl_shi0to3(ix,iy,3,istrai)*ni(ix,iy)*
     &                                                ti(ix,iy)
              eirene_mc_eppl_shi(ix,iy,istrai)=
     &         eirene_mc_eppl_shi0to3(ix,iy,0,istrai)+
     &         eirene_mc_eppl_shi0to3(ix,iy,1,istrai)*ti(ix,iy)+
     &         eirene_mc_eppl_shi0to3(ix,iy,2,istrai)*ni(ix,iy)+
     &         eirene_mc_eppl_shi0to3(ix,iy,3,istrai)*ni(ix,iy)*
     &                                                ti(ix,iy)
            enddo
          endif
          b2stbm_shi_bal(ix,iy)=
     &     b2stbm_shi0to3(ix,iy,0)+
     &     b2stbm_shi0to3(ix,iy,1)*ti(ix,iy)+
     &     b2stbm_shi0to3(ix,iy,2)*ni(ix,iy)+
     &     b2stbm_shi0to3(ix,iy,3)*ni(ix,iy)*ti(ix,iy)
          ext_shi_bal(ix,iy)=
     &     ext_shi0to3(ix,iy,0)+
     &     ext_shi0to3(ix,iy,1)*ti(ix,iy)+
     &     ext_shi0to3(ix,iy,2)*ni(ix,iy)+
     &     ext_shi0to3(ix,iy,3)*ni(ix,iy)*ti(ix,iy)
          b2stel_shi_ion_bal(ix,iy)=
     &     b2stel_shi_ion0to3(ix,iy,0)+
     &     b2stel_shi_ion0to3(ix,iy,1)*ti(ix,iy)+
     &     b2stel_shi_ion0to3(ix,iy,2)*ni(ix,iy)+
     &     b2stel_shi_ion0to3(ix,iy,3)*ni(ix,iy)*ti(ix,iy)
          b2stel_shi_rec_bal(ix,iy)=
     &     b2stel_shi_rec0to3(ix,iy,0)+
     &     b2stel_shi_rec0to3(ix,iy,1)*ti(ix,iy)+
     &     b2stel_shi_rec0to3(ix,iy,2)*ni(ix,iy)+
     &     b2stel_shi_rec0to3(ix,iy,3)*ni(ix,iy)*ti(ix,iy)
          b2stcx_shi_bal(ix,iy)=
     &     b2stcx_shi0to3(ix,iy,0)+
     &     b2stcx_shi0to3(ix,iy,1)*ti(ix,iy)+
     &     b2stcx_shi0to3(ix,iy,2)*ni(ix,iy)+
     &     b2stcx_shi0to3(ix,iy,3)*ni(ix,iy)*ti(ix,iy)
          b2srsm_shi_bal(ix,iy)=
     &     b2srsm_shi0to3(ix,iy,0)+
     &     b2srsm_shi0to3(ix,iy,1)*ti(ix,iy)+
     &     b2srsm_shi0to3(ix,iy,2)*ni(ix,iy)+
     &     b2srsm_shi0to3(ix,iy,3)*ni(ix,iy)*ti(ix,iy)
          b2srdt_shi_bal(ix,iy)=
     &     b2srdt_shi0to3(ix,iy,0)+
     &     b2srdt_shi0to3(ix,iy,1)*ti(ix,iy)+
     &     b2srdt_shi0to3(ix,iy,2)*ni(ix,iy)+
     &     b2srdt_shi0to3(ix,iy,3)*ni(ix,iy)*ti(ix,iy)
          b2srst_shi_bal(ix,iy)=
     &     b2srst_shi0to3(ix,iy,0)+
     &     b2srst_shi0to3(ix,iy,1)*ti(ix,iy)+
     &     b2srst_shi0to3(ix,iy,2)*ni(ix,iy)+
     &     b2srst_shi0to3(ix,iy,3)*ni(ix,iy)*ti(ix,iy)
          b2sihs_diaa_bal(ix,iy)=
     &     b2sihs_diaa0to3(ix,iy,0)+
     &     b2sihs_diaa0to3(ix,iy,1)*ti(ix,iy)+
     &     b2sihs_diaa0to3(ix,iy,2)*ni(ix,iy)+
     &     b2sihs_diaa0to3(ix,iy,3)*ni(ix,iy)*ti(ix,iy)
          b2sihs_divua_bal(ix,iy)=
     &     b2sihs_divua0to3(ix,iy,0)+
     &     b2sihs_divua0to3(ix,iy,1)*ti(ix,iy)+
     &     b2sihs_divua0to3(ix,iy,2)*ni(ix,iy)+
     &     b2sihs_divua0to3(ix,iy,3)*ni(ix,iy)*ti(ix,iy)
          b2sihs_exba_bal(ix,iy)=
     &     b2sihs_exba0to3(ix,iy,0)+
     &     b2sihs_exba0to3(ix,iy,1)*ti(ix,iy)+
     &     b2sihs_exba0to3(ix,iy,2)*ni(ix,iy)+
     &     b2sihs_exba0to3(ix,iy,3)*ni(ix,iy)*ti(ix,iy)
          b2sihs_visa_bal(ix,iy)=
     &     b2sihs_visa0to3(ix,iy,0)+
     &     b2sihs_visa0to3(ix,iy,1)*ti(ix,iy)+
     &     b2sihs_visa0to3(ix,iy,2)*ni(ix,iy)+
     &     b2sihs_visa0to3(ix,iy,3)*ni(ix,iy)*ti(ix,iy)
          b2sihs_fraa_bal(ix,iy)=
     &     b2sihs_fraa0to3(ix,iy,0)+
     &     b2sihs_fraa0to3(ix,iy,1)*ti(ix,iy)+
     &     b2sihs_fraa0to3(ix,iy,2)*ni(ix,iy)+
     &     b2sihs_fraa0to3(ix,iy,3)*ni(ix,iy)*ti(ix,iy)
        enddo
      enddo
      return
      end subroutine update_balance_sources_heat

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
cdjm Arrays for particle balance
        call cfrure (nread, size(b2stbr_sna_bal),
     .               b2stbr_sna_bal, 'b2stbr_sna_bal')
        call cfrure (nread, size(b2stbc_sna_bal),
     .               b2stbc_sna_bal, 'b2stbc_sna_bal')
        call cfrure (nread, size(b2stbm_sna_bal),
     .               b2stbm_sna_bal, 'b2stbm_sna_bal')
        call cfrure (nread, size(ext_sna_bal),
     .               ext_sna_bal, 'ext_sna_bal')
        call cfrure (nread, size(b2stel_sna_ion_bal),
     .                b2stel_sna_ion_bal, 'b2stel_sna_ion_bal')
        call cfrure (nread, size(b2stel_sna_rec_bal),
     .                b2stel_sna_rec_bal, 'b2stel_sna_rec_bal')
        call cfrure (nread, size(b2stcx_sna_bal),
     .                b2stcx_sna_bal, 'b2stcx_sna_bal')
        call cfrure (nread, size(b2srsm_sna_bal),
     .                b2srsm_sna_bal, 'b2srsm_sna_bal')
        call cfrure (nread, size(b2srdt_sna_bal),
     .                b2srdt_sna_bal, 'b2srdt_sna_bal')
        call cfrure (nread, size(b2srst_sna_bal),
     .                b2srst_sna_bal, 'b2srst_sna_bal')
cdjm Arrays for momentum balance
        call cfrure (nread, size(b2stbr_smo_bal),
     .                b2stbr_smo_bal, 'b2stbr_smo_bal')
        call cfrure (nread, size(b2stbc_smo_bal),
     .                b2stbc_smo_bal, 'b2stbc_smo_bal')
        call cfrure (nread, size(b2stbm_smo_bal),
     .                b2stbm_smo_bal, 'b2stbm_smo_bal')
        call cfrure (nread, size(ext_smo_bal),
     .                ext_smo_bal, 'ext_smo_bal')
        call cfrure (nread, size(b2stel_smq_ion_bal),
     .                b2stel_smq_ion_bal, 'b2stel_smq_ion_bal')
        call cfrure (nread, size(b2stel_smq_rec_bal),
     .                b2stel_smq_rec_bal, 'b2stel_smq_rec_bal')
        call cfrure (nread, size(b2stcx_smq_bal),
     .                b2stcx_smq_bal, 'b2stcx_smq_bal')
        call cfrure (nread, size(b2srsm_smo_bal),
     .                b2srsm_smo_bal, 'b2srsm_smo_bal')
        call cfrure (nread, size(b2srdt_smo_bal),
     .                b2srdt_smo_bal, 'b2srdt_smo_bal')
        call cfrure (nread, size(b2srst_smo_bal),
     .                b2srst_smo_bal, 'b2srst_smo_bal')
        call cfrure (nread, size(b2sifr_smoch_bal),
     .               b2sifr_smoch_bal, 'b2sifr_smoch_bal')
        call cfrure (nread, size(b2sifr_smotf_bal),
     .               b2sifr_smotf_bal, 'b2sifr_smotf_bal')
        call cfrure (nread, size(b2siav_smovh_bal),
     .               b2siav_smovh_bal, 'b2siav_smovh_bal')
        call cfrure (nread, size(b2siav_smovv_bal),
     .               b2siav_smovv_bal, 'b2siav_smovv_bal')
        call cfrure (nread, size(b2sicf_smo_bal),
     .               b2sicf_smo_bal, 'b2sicf_smo_bal')
        call cfrure (nread, size(b2nxdv_smo_bal),
     .               b2nxdv_smo_bal, 'b2nxdv_smo_bal')
        call cfrure (nread, size(b2sigp_smogp_bal),
     .               b2sigp_smogp_bal, 'b2sigp_smogp_bal')
        call cfrure (nread, size(b2sigp_pstat_bal),
     .               b2sigp_pstat_bal, 'b2sigp_pstat_bal')
cdjm Arrays for electron heat balance
        call cfrure (nread, size(b2stbr_she_bal),
     .                b2stbr_she_bal, 'b2stbr_she_bal')
        call cfrure (nread, size(b2stbc_she_bal),
     .                b2stbc_she_bal, 'b2stbc_she_bal')
        call cfrure (nread, size(b2stbm_she_bal),
     .                b2stbm_she_bal, 'b2stbm_she_bal')
        call cfrure (nread, size(ext_she_bal),
     .                ext_she_bal, 'ext_she_bal')
        call cfrure (nread, size(b2stel_she_bal),
     .                b2stel_she_bal, 'b2stel_she_bal')
        call cfrure (nread, size(b2srsm_she_bal),
     .                b2srsm_she_bal, 'b2srsm_she_bal')
        call cfrure (nread, size(b2srdt_she_bal),
     .                b2srdt_she_bal, 'b2srdt_she_bal')
        call cfrure (nread, size(b2srst_she_bal),
     .                b2srst_she_bal, 'b2srst_she_bal')
        call cfrure (nread, size(b2sihs_diae_bal),
     .                b2sihs_diae_bal, 'b2sihs_diae_bal')
        call cfrure (nread, size(b2sihs_divue_bal),
     .                b2sihs_divue_bal, 'b2sihs_divue_bal')
        call cfrure (nread, size(b2sihs_exbe_bal),
     .                b2sihs_exbe_bal, 'b2sihs_exbe_bal')
        call cfrure (nread, size(b2sihs_joule_bal),
     .                b2sihs_joule_bal, 'b2sihs_joule_bal')
        call cfrure (nread, size(b2npht_shei_bal),
     .               b2npht_shei_bal, 'b2npht_shei_bal')
cdjm Arrays for ion heat balance
        call cfrure (nread, size(b2stbr_shi_bal),
     .                b2stbr_shi_bal, 'b2stbr_shi_bal')
        call cfrure (nread, size(b2stbc_shi_bal),
     .                b2stbc_shi_bal, 'b2stbc_shi_bal')
        call cfrure (nread, size(b2stbm_shi_bal),
     .                b2stbm_shi_bal, 'b2stbm_shi_bal')
        call cfrure (nread, size(ext_shi_bal),
     .                ext_shi_bal, 'ext_shi_bal')
        call cfrure (nread, size(b2stel_shi_ion_bal),
     .                b2stel_shi_ion_bal, 'b2stel_shi_ion_bal')
        call cfrure (nread, size(b2stel_shi_rec_bal),
     .                b2stel_shi_rec_bal, 'b2stel_shi_rec_bal')
        call cfrure (nread, size(b2stcx_shi_bal),
     .                b2stcx_shi_bal, 'b2stcx_shi_bal')
        call cfrure (nread, size(b2srsm_shi_bal),
     .                b2srsm_shi_bal, 'b2srsm_shi_bal')
        call cfrure (nread, size(b2srdt_shi_bal),
     .                b2srdt_shi_bal, 'b2srdt_shi_bal')
        call cfrure (nread, size(b2srst_shi_bal),
     .                b2srst_shi_bal, 'b2srst_shi_bal')
        call cfrure (nread, size(b2sihs_diaa_bal),
     .               b2sihs_diaa_bal, 'b2sihs_diaa_bal')
        call cfrure (nread, size(b2sihs_divua_bal),
     .               b2sihs_divua_bal, 'b2sihs_divua_bal')
        call cfrure (nread, size(b2sihs_exba_bal),
     .               b2sihs_exba_bal, 'b2sihs_exba_bal')
        call cfrure (nread, size(b2sihs_visa_bal),
     .               b2sihs_visa_bal, 'b2sihs_visa_bal')
        call cfrure (nread, size(b2sihs_fraa_bal),
     .               b2sihs_fraa_bal, 'b2sihs_fraa_bal')
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
        if (.not.allocated(eirene_mc_papl_sna))
     >   call alloc_b2mod_eirene_sources(nxd,nyd,nsd,nstrad)
        call cfrure_opt (nread, size(eirene_mc_papl_sna),
     &                   eirene_mc_papl_sna, 'eirene_mc_papl_sna')
        call cfrure_opt (nread, size(eirene_mc_pmpl_sna),
     &                   eirene_mc_pmpl_sna, 'eirene_mc_pmpl_sna')
        call cfrure_opt (nread, size(eirene_mc_pipl_sna),
     &                   eirene_mc_pipl_sna, 'eirene_mc_pipl_sna')
        call cfrure_opt (nread, size(eirene_mc_pppl_sna),
     &                   eirene_mc_pppl_sna, 'eirene_mc_pppl_sna')
        call cfrure_opt (nread, size(eirene_mc_core_sna),
     &                   eirene_mc_core_sna, 'eirene_mc_core_sna')
        call cfrure_opt (nread, size(eirene_mc_mapl_smo),
     &                   eirene_mc_mapl_smo, 'eirene_mc_mapl_smo')
        call cfrure_opt (nread, size(eirene_mc_mmpl_smo),
     &                   eirene_mc_mmpl_smo, 'eirene_mc_mmpl_smo')
        call cfrure_opt (nread, size(eirene_mc_mipl_smo),
     &                   eirene_mc_mipl_smo, 'eirene_mc_mipl_smo')
        call cfrure_opt (nread, size(eirene_mc_cppv_smo),
     &                   eirene_mc_cppv_smo, 'eirene_mc_cppv_smo')
        call cfrure_opt (nread, size(eirene_mc_eael_she),
     &                   eirene_mc_eael_she, 'eirene_mc_eael_she')
        call cfrure_opt (nread, size(eirene_mc_emel_she),
     &                   eirene_mc_emel_she, 'eirene_mc_emel_she')
        call cfrure_opt (nread, size(eirene_mc_eiel_she),
     &                   eirene_mc_eiel_she, 'eirene_mc_eiel_she')
        call cfrure_opt (nread, size(eirene_mc_epel_she),
     &                   eirene_mc_epel_she, 'eirene_mc_epel_she')
        call cfrure_opt (nread, size(eirene_mc_eapl_shi),
     &                   eirene_mc_eapl_shi, 'eirene_mc_eapl_shi')
        call cfrure_opt (nread, size(eirene_mc_empl_shi),
     &                   eirene_mc_empl_shi, 'eirene_mc_empl_shi')
        call cfrure_opt (nread, size(eirene_mc_eipl_shi),
     &                   eirene_mc_eipl_shi, 'eirene_mc_eipl_shi')
        call cfrure_opt (nread, size(eirene_mc_eppl_shi),
     &                   eirene_mc_eppl_shi, 'eirene_mc_eppl_shi')
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
cdjm Arrays for particle balance
        call cfwure (nwrite, size(b2stbr_sna_bal),
     .               b2stbr_sna_bal, 'b2stbr_sna_bal')
        call cfwure (nwrite, size(b2stbc_sna_bal),
     .               b2stbc_sna_bal, 'b2stbc_sna_bal')
        call cfwure (nwrite, size(b2stbm_sna_bal),
     .               b2stbm_sna_bal, 'b2stbm_sna_bal')
        call cfwure (nwrite, size(ext_sna_bal),
     .               ext_sna_bal, 'ext_sna_bal')
        call cfwure (nwrite, size(b2stel_sna_ion_bal),
     .                b2stel_sna_ion_bal, 'b2stel_sna_ion_bal')
        call cfwure (nwrite, size(b2stel_sna_rec_bal),
     .                b2stel_sna_rec_bal, 'b2stel_sna_rec_bal')
        call cfwure (nwrite, size(b2stcx_sna_bal),
     .                b2stcx_sna_bal, 'b2stcx_sna_bal')
        call cfwure (nwrite, size(b2srsm_sna_bal),
     .                b2srsm_sna_bal, 'b2srsm_sna_bal')
        call cfwure (nwrite, size(b2srdt_sna_bal),
     .                b2srdt_sna_bal, 'b2srdt_sna_bal')
        call cfwure (nwrite, size(b2srst_sna_bal),
     .                b2srst_sna_bal, 'b2srst_sna_bal')
cdjm Arrays for momentum balance
        call cfwure (nwrite, size(b2stbr_smo_bal),
     .                b2stbr_smo_bal, 'b2stbr_smo_bal')
        call cfwure (nwrite, size(b2stbc_smo_bal),
     .                b2stbc_smo_bal, 'b2stbc_smo_bal')
        call cfwure (nwrite, size(b2stbm_smo_bal),
     .                b2stbm_smo_bal, 'b2stbm_smo_bal')
        call cfwure (nwrite, size(ext_smo_bal),
     .                ext_smo_bal, 'ext_smo_bal')
        call cfwure (nwrite, size(b2stel_smq_ion_bal),
     .                b2stel_smq_ion_bal, 'b2stel_smq_ion_bal')
        call cfwure (nwrite, size(b2stel_smq_rec_bal),
     .                b2stel_smq_rec_bal, 'b2stel_smq_rec_bal')
        call cfwure (nwrite, size(b2stcx_smq_bal),
     .                b2stcx_smq_bal, 'b2stcx_smq_bal')
        call cfwure (nwrite, size(b2srsm_smo_bal),
     .                b2srsm_smo_bal, 'b2srsm_smo_bal')
        call cfwure (nwrite, size(b2srdt_smo_bal),
     .                b2srdt_smo_bal, 'b2srdt_smo_bal')
        call cfwure (nwrite, size(b2srst_smo_bal),
     .                b2srst_smo_bal, 'b2srst_smo_bal')
        call cfwure (nwrite, size(b2sifr_smoch_bal),
     .               b2sifr_smoch_bal, 'b2sifr_smoch_bal')
        call cfwure (nwrite, size(b2sifr_smotf_bal),
     .               b2sifr_smotf_bal, 'b2sifr_smotf_bal')
        call cfwure (nwrite, size(b2siav_smovh_bal),
     .               b2siav_smovh_bal, 'b2siav_smovh_bal')
        call cfwure (nwrite, size(b2siav_smovv_bal),
     .               b2siav_smovv_bal, 'b2siav_smovv_bal')
        call cfwure (nwrite, size(b2sicf_smo_bal),
     .               b2sicf_smo_bal, 'b2sicf_smo_bal')
        call cfwure (nwrite, size(b2nxdv_smo_bal),
     .               b2nxdv_smo_bal, 'b2nxdv_smo_bal')
        call cfwure (nwrite, size(b2sigp_smogp_bal),
     .               b2sigp_smogp_bal, 'b2sigp_smogp_bal')
        call cfwure (nwrite, size(b2sigp_pstat_bal),
     .               b2sigp_pstat_bal, 'b2sigp_pstat_bal')
cdjm Arrays for electron heat balance
        call cfwure (nwrite, size(b2stbr_she_bal),
     .                b2stbr_she_bal, 'b2stbr_she_bal')
        call cfwure (nwrite, size(b2stbc_she_bal),
     .                b2stbc_she_bal, 'b2stbc_she_bal')
        call cfwure (nwrite, size(b2stbm_she_bal),
     .                b2stbm_she_bal, 'b2stbm_she_bal')
        call cfwure (nwrite, size(ext_she_bal),
     .                ext_she_bal, 'ext_she_bal')
        call cfwure (nwrite, size(b2stel_she_bal),
     .                b2stel_she_bal, 'b2stel_she_bal')
        call cfwure (nwrite, size(b2srsm_she_bal),
     .                b2srsm_she_bal, 'b2srsm_she_bal')
        call cfwure (nwrite, size(b2srdt_she_bal),
     .                b2srdt_she_bal, 'b2srdt_she_bal')
        call cfwure (nwrite, size(b2srst_she_bal),
     .                b2srst_she_bal, 'b2srst_she_bal')
        call cfwure (nwrite, size(b2sihs_diae_bal),
     .                b2sihs_diae_bal, 'b2sihs_diae_bal')
        call cfwure (nwrite, size(b2sihs_divue_bal),
     .                b2sihs_divue_bal, 'b2sihs_divue_bal')
        call cfwure (nwrite, size(b2sihs_exbe_bal),
     .                b2sihs_exbe_bal, 'b2sihs_exbe_bal')
        call cfwure (nwrite, size(b2sihs_joule_bal),
     .                b2sihs_joule_bal, 'b2sihs_joule_bal')
        call cfwure (nwrite, size(b2npht_shei_bal),
     .               b2npht_shei_bal, 'b2npht_shei_bal')
cdjm Arrays for ion heat balance
        call cfwure (nwrite, size(b2stbr_shi_bal),
     .                b2stbr_shi_bal, 'b2stbr_shi_bal')
        call cfwure (nwrite, size(b2stbc_shi_bal),
     .                b2stbc_shi_bal, 'b2stbc_shi_bal')
        call cfwure (nwrite, size(b2stbm_shi_bal),
     .                b2stbm_shi_bal, 'b2stbm_shi_bal')
        call cfwure (nwrite, size(ext_shi_bal),
     .                ext_shi_bal, 'ext_shi_bal')
        call cfwure (nwrite, size(b2stel_shi_ion_bal),
     .                b2stel_shi_ion_bal, 'b2stel_shi_ion_bal')
        call cfwure (nwrite, size(b2stel_shi_rec_bal),
     .                b2stel_shi_rec_bal, 'b2stel_shi_rec_bal')
        call cfwure (nwrite, size(b2stcx_shi_bal),
     .                b2stcx_shi_bal, 'b2stcx_shi_bal')
        call cfwure (nwrite, size(b2srsm_shi_bal),
     .                b2srsm_shi_bal, 'b2srsm_shi_bal')
        call cfwure (nwrite, size(b2srdt_shi_bal),
     .                b2srdt_shi_bal, 'b2srdt_shi_bal')
        call cfwure (nwrite, size(b2srst_shi_bal),
     .                b2srst_shi_bal, 'b2srst_shi_bal')
        call cfwure (nwrite, size(b2sihs_diaa_bal),
     .               b2sihs_diaa_bal, 'b2sihs_diaa_bal')
        call cfwure (nwrite, size(b2sihs_divua_bal),
     .               b2sihs_divua_bal, 'b2sihs_divua_bal')
        call cfwure (nwrite, size(b2sihs_exba_bal),
     .               b2sihs_exba_bal, 'b2sihs_exba_bal')
        call cfwure (nwrite, size(b2sihs_visa_bal),
     .               b2sihs_visa_bal, 'b2sihs_visa_bal')
        call cfwure (nwrite, size(b2sihs_fraa_bal),
     .               b2sihs_fraa_bal, 'b2sihs_fraa_bal')
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

      if (.not.allocated(eirene_mc_papl_sna)) return
      call cfwure (nwrite, size(eirene_mc_papl_sna),
     &             eirene_mc_papl_sna, 'eirene_mc_papl_sna')
      call cfwure (nwrite, size(eirene_mc_pmpl_sna),
     &             eirene_mc_pmpl_sna, 'eirene_mc_pmpl_sna')
      call cfwure (nwrite, size(eirene_mc_pipl_sna),
     &             eirene_mc_pipl_sna, 'eirene_mc_pipl_sna')
      call cfwure (nwrite, size(eirene_mc_pppl_sna),
     &             eirene_mc_pppl_sna, 'eirene_mc_pppl_sna')
      call cfwure (nwrite, size(eirene_mc_core_sna),
     &             eirene_mc_core_sna, 'eirene_mc_core_sna')
      call cfwure (nwrite, size(eirene_mc_mapl_smo),
     &             eirene_mc_mapl_smo, 'eirene_mc_mapl_smo')
      call cfwure (nwrite, size(eirene_mc_mmpl_smo),
     &             eirene_mc_mmpl_smo, 'eirene_mc_mmpl_smo')
      call cfwure (nwrite, size(eirene_mc_mipl_smo),
     &             eirene_mc_mipl_smo, 'eirene_mc_mipl_smo')
      call cfwure (nwrite, size(eirene_mc_cppv_smo),
     &             eirene_mc_cppv_smo, 'eirene_mc_cppv_smo')
      call cfwure (nwrite, size(eirene_mc_eael_she),
     &             eirene_mc_eael_she, 'eirene_mc_eael_she')
      call cfwure (nwrite, size(eirene_mc_emel_she),
     &             eirene_mc_emel_she, 'eirene_mc_emel_she')
      call cfwure (nwrite, size(eirene_mc_eiel_she),
     &             eirene_mc_eiel_she, 'eirene_mc_eiel_she')
      call cfwure (nwrite, size(eirene_mc_epel_she),
     &             eirene_mc_epel_she, 'eirene_mc_epel_she')
      call cfwure (nwrite, size(eirene_mc_eapl_shi),
     &             eirene_mc_eapl_shi, 'eirene_mc_eapl_shi')
      call cfwure (nwrite, size(eirene_mc_empl_shi),
     &             eirene_mc_empl_shi, 'eirene_mc_empl_shi')
      call cfwure (nwrite, size(eirene_mc_eipl_shi),
     &             eirene_mc_eipl_shi, 'eirene_mc_eipl_shi')
      call cfwure (nwrite, size(eirene_mc_eppl_shi),
     &             eirene_mc_eppl_shi, 'eirene_mc_eppl_shi')

      return
      end subroutine write_b2mod_eirene_sources

      end module b2mod_sources
