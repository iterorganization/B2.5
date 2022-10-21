      module b2mod_balance
      use b2mod_types
      implicit none

*   ..Flux and source arrays required for balance postprocessing routines. For the
*     sources, the linearised components are stored and then the total source is updated
*     just before each residual calculation.
*     Only this total source is output to balance.nc
*   ..Arrays for particle balance:
*     Fluxes
      real (kind=R8), allocatable, save ::
     *  fna_pinch(:,:,:),
     *  fna_pll(:,:,:),
     *  fna_drift(:,:,:),
     *  fna_ch(:,:,:),
     *  fna_nanom(:,:,:),
     *  fna_panom(:,:,:),
     *  fna_pschused(:,:,:)
*     Sources
      real (kind=R8), allocatable, save ::
     *  b2stbr_phys_sna0to1(:,:,:),
     *  b2stbr_phys_sna_bal(:,:),
     *  b2stbr_bas_sna0to1(:,:,:),
     *  b2stbr_bas_sna_bal(:,:),
     *  b2stbr_first_flight_sna0to1(:,:,:),
     *  b2stbr_first_flight_sna_bal(:,:),
     *  b2stbc_sna0to1(:,:,:),
     *  b2stbc_sna_bal(:,:),
     *  eirene_mc_papl_sna0to1(:,:,:,:),
     *  eirene_mc_pmpl_sna0to1(:,:,:,:),
     *  eirene_mc_pipl_sna0to1(:,:,:,:),
     *  eirene_mc_pppl_sna0to1(:,:,:,:),
     *  eirene_mc_core_sna0to1(:,:,:),
     *  eirene_mc_papl_sna_bal(:,:,:),
     *  eirene_mc_pmpl_sna_bal(:,:,:),
     *  eirene_mc_pipl_sna_bal(:,:,:),
     *  eirene_mc_pppl_sna_bal(:,:,:),
     *  eirene_mc_core_sna_bal(:,:),
     *  b2stbm_sna0to1(:,:,:),
     *  b2stbm_sna_bal(:,:),
     *  ext_sna0to1(:,:,:),
     *  ext_sna_bal(:,:),
     *  b2stel_sna_ion0to1(:,:,:),
     *  b2stel_sna_ion_bal(:,:),
     *  b2stel_sna_rec0to1(:,:,:),
     *  b2stel_sna_rec_bal(:,:),
     *  b2stcx_sna0to1(:,:,:),
     *  b2stcx_sna_bal(:,:),
     *  b2srsm_sna0to1(:,:,:),
     *  b2srsm_sna_bal(:,:),
     *  b2srdt_sna0to1(:,:,:),
     *  b2srdt_sna_bal(:,:),
     *  b2srst_sna0to1(:,:,:),
     *  b2srst_sna_bal(:,:)
*     Time-averaged quantities
      real (kind=R8), allocatable, save ::
     *  fna_pinch_av(:,:,:),
     *  fna_pll_av(:,:,:),
     *  fna_drift_av(:,:,:),
     *  fna_ch_av(:,:,:),
     *  fna_nanom_av(:,:,:),
     *  fna_panom_av(:,:,:),
     *  fna_pschused_av(:,:,:),
     *  b2stbr_phys_sna_bal_av(:,:),
     *  b2stbr_bas_sna_bal_av(:,:),
     *  b2stbr_first_flight_sna_bal_av(:,:),
     *  b2stbc_sna_bal_av(:,:),
     *  eirene_mc_papl_sna_bal_av(:,:,:),
     *  eirene_mc_pmpl_sna_bal_av(:,:,:),
     *  eirene_mc_pipl_sna_bal_av(:,:,:),
     *  eirene_mc_pppl_sna_bal_av(:,:,:),
     *  eirene_mc_core_sna_bal_av(:,:),
     *  b2stbm_sna_bal_av(:,:),
     *  ext_sna_bal_av(:,:),
     *  b2stel_sna_ion_bal_av(:,:),
     *  b2stel_sna_rec_bal_av(:,:),
     *  b2stcx_sna_bal_av(:,:),
     *  b2srsm_sna_bal_av(:,:),
     *  b2srdt_sna_bal_av(:,:),
     *  b2srst_sna_bal_av(:,:),
     *  resco_av(:,:)
*   ..Arrays for momentum balance:
*     Fluxes
      real (kind=R8), allocatable, save ::
     *  fmo_flua(:,:,:),
     *  fmo_cvsa(:,:,:),
     *  fmo_b2nxfv(:,:,:)
*     Sources
      real (kind=R8), allocatable, save ::
     *  b2stbr_phys_smo0to3(:,:,:),
     *  b2stbr_phys_smo_bal(:,:),
     *  b2stbr_bas_smo0to3(:,:,:),
     *  b2stbr_bas_smo_bal(:,:),
     *  b2stbc_smo0to3(:,:,:),
     *  b2stbc_smo_bal(:,:),
     *  eirene_mc_mapl_smo0to3(:,:,:,:),
     *  eirene_mc_mmpl_smo0to3(:,:,:,:),
     *  eirene_mc_mipl_smo0to3(:,:,:,:),
     *  eirene_mc_cppv_smo0to3(:,:,:,:),
     *  eirene_mc_mapl_smo_bal(:,:,:),
     *  eirene_mc_mmpl_smo_bal(:,:,:),
     *  eirene_mc_mipl_smo_bal(:,:,:),
     *  eirene_mc_cppv_smo_bal(:,:,:),
     *  b2stbm_smo0to3(:,:,:),
     *  b2stbm_smo_bal(:,:),
     *  ext_smo0to3(:,:,:),
     *  ext_smo_bal(:,:),
     *  b2stel_smq_ion0to3(:,:,:),
     *  b2stel_smq_ion_bal(:,:),
     *  b2stel_smq_rec0to3(:,:,:),
     *  b2stel_smq_rec_bal(:,:),
     *  b2stcx_smq0to3(:,:,:),
     *  b2stcx_smq_bal(:,:),
     *  b2srsm_smo0to3(:,:,:),
     *  b2srsm_smo_bal(:,:),
     *  b2srdt_smo0to3(:,:,:),
     *  b2srdt_smo_bal(:,:),
     *  b2srst_smo0to3(:,:,:),
     *  b2srst_smo_bal(:,:),
     *  b2sifr_smoch0to3(:,:,:),
     *  b2sifr_smoch_bal(:,:),
     *  b2sifr_smotf0to3(:,:,:),
     *  b2sifr_smotf_bal(:,:),
     *  b2siav_smovh0to3(:,:,:),
     *  b2siav_smovh_bal(:,:),
     *  b2siav_smovv0to3(:,:,:),
     *  b2siav_smovv_bal(:,:),
     *  b2sicf_smo0to3(:,:,:),
     *  b2sicf_smo_bal(:,:),
     *  b2sian_smo0to3(:,:,:),               !srv 25.09.17
     *  b2sian_smo_bal(:,:),                 !srv 25.09.17
     *  b2nxdv_smo0to3(:,:,:),
     *  b2nxdv_smo_bal(:,:),
     *  b2sigp_smogp0to3(:,:,:),
     *  b2sigp_smogp_bal(:,:),
     *  b2sigp_pstat0to3(:,:,:),
     *  b2sigp_pstat_bal(:,:)
*     Time-averaged quantities
      real (kind=R8), allocatable, save ::
     *  fmo_flua_av(:,:,:),
     *  fmo_cvsa_av(:,:,:),
     *  fmo_b2nxfv_av(:,:,:),
     *  b2stbr_phys_smo_bal_av(:,:),
     *  b2stbr_bas_smo_bal_av(:,:),
     *  b2stbc_smo_bal_av(:,:),
     *  eirene_mc_mapl_smo_bal_av(:,:,:),
     *  eirene_mc_mmpl_smo_bal_av(:,:,:),
     *  eirene_mc_mipl_smo_bal_av(:,:,:),
     *  eirene_mc_cppv_smo_bal_av(:,:,:),
     *  b2stbm_smo_bal_av(:,:),
     *  ext_smo_bal_av(:,:),
     *  b2stel_smq_ion_bal_av(:,:),
     *  b2stel_smq_rec_bal_av(:,:),
     *  b2stcx_smq_bal_av(:,:),
     *  b2srsm_smo_bal_av(:,:),
     *  b2srdt_smo_bal_av(:,:),
     *  b2srst_smo_bal_av(:,:),
     *  b2sifr_smoch_bal_av(:,:),
     *  b2sifr_smotf_bal_av(:,:),
     *  b2siav_smovh_bal_av(:,:),
     *  b2siav_smovv_bal_av(:,:),
     *  b2sicf_smo_bal_av(:,:),
     *  b2sian_smo_bal_av(:,:),              !srv 25.07.17
     *  b2nxdv_smo_bal_av(:,:),
     *  b2sigp_smogp_bal_av(:,:),
     *  b2sigp_pstat_bal_av(:,:),
     *  resmo_av(:,:)
*   ..Arrays for electron heat balance:
*     Fluxes
      real (kind=R8), allocatable, save ::
     *  fhe_32(:,:),
     *  fhe_52(:,:),
     *  fhe_thermj(:,:),
     *  fhe_cond(:,:),
     *  fhe_dia(:,:),
     *  fhe_ecrb(:,:),
     *  fhe_strange(:,:),
     *  fhe_pschused(:,:)
*     Sources
      real (kind=R8), allocatable, save ::
     *  b2stbr_phys_she0to3(:,:),
     *  b2stbr_phys_she_bal(:),
     *  b2stbr_bas_she0to3(:,:),
     *  b2stbr_bas_she_bal(:),
     *  b2stbr_first_flight_she0to3(:,:),
     *  b2stbr_first_flight_she_bal(:),
     *  b2stbc_she0to3(:,:),
     *  b2stbc_she_bal(:),
     *  eirene_mc_eael_she0to3(:,:,:),
     *  eirene_mc_emel_she0to3(:,:,:),
     *  eirene_mc_eiel_she0to3(:,:,:),
     *  eirene_mc_epel_she0to3(:,:,:),
     *  eirene_mc_eael_she_bal(:,:),
     *  eirene_mc_emel_she_bal(:,:),
     *  eirene_mc_eiel_she_bal(:,:),
     *  eirene_mc_epel_she_bal(:,:),
     *  b2stbm_she0to3(:,:),
     *  b2stbm_she_bal(:),
     *  ext_she0to3(:,:),
     *  ext_she_bal(:),
     *  b2stel_she0to3(:,:),
     *  b2stel_she_bal(:),
     *  b2srsm_she0to3(:,:),
     *  b2srsm_she_bal(:),
     *  b2srdt_she0to3(:,:),
     *  b2srdt_she_bal(:),
     *  b2srst_she0to3(:,:),
     *  b2srst_she_bal(:),
     *  b2sihs_diae0to3(:,:),
     *  b2sihs_diae_bal(:),
     *  b2sihs_divue0to3(:,:),
     *  b2sihs_divue_bal(:),
     *  b2sihs_exbe0to3(:,:),
     *  b2sihs_exbe_bal(:),
     *  b2sihs_joule0to3(:,:),
     *  b2sihs_joule_bal(:),
     *  b2npht_shei_bal(:) !djm equipartition is not linearised
*     Time-averaged quantities
      real (kind=R8), allocatable, save ::
     *  fhe_32_av(:,:),
     *  fhe_52_av(:,:),
     *  fhe_thermj_av(:,:),
     *  fhe_cond_av(:,:),
     *  fhe_dia_av(:,:),
     *  fhe_ecrb_av(:,:),
     *  fhe_strange_av(:,:),
     *  fhe_pschused_av(:,:),
     *  b2stbr_phys_she_bal_av(:),
     *  b2stbr_bas_she_bal_av(:),
     *  b2stbr_first_flight_she_bal_av(:),
     *  b2stbc_she_bal_av(:),
     *  eirene_mc_eael_she_bal_av(:,:),
     *  eirene_mc_emel_she_bal_av(:,:),
     *  eirene_mc_eiel_she_bal_av(:,:),
     *  eirene_mc_epel_she_bal_av(:,:),
     *  b2stbm_she_bal_av(:),
     *  ext_she_bal_av(:),
     *  b2stel_she_bal_av(:),
     *  b2srsm_she_bal_av(:),
     *  b2srdt_she_bal_av(:),
     *  b2srst_she_bal_av(:),
     *  b2sihs_diae_bal_av(:),
     *  b2sihs_divue_bal_av(:),
     *  b2sihs_exbe_bal_av(:),
     *  b2sihs_joule_bal_av(:),
     *  b2npht_shei_bal_av(:),
     *  reshe_av(:)
*   ..Arrays for ion heat balance:
*     Fluxes
      real (kind=R8), allocatable, save ::
     *  fhi_32(:,:),
     *  fhi_52(:,:),
     *  fhi_cond(:,:),
     *  fhi_dia(:,:),
     *  fhi_ecrb(:,:),
     *  fhi_strange(:,:),
     *  fhi_pschused(:,:),
     *  fhi_inert(:,:),
     *  fhi_vispar(:,:),
     *  fhi_visper(:,:),
     *  fhi_visq(:,:),
     *  fhi_anml(:,:),
     *  fhi_kevis(:,:)
*     Source
      real (kind=R8), allocatable, save ::
     *  b2stbr_phys_shi0to3(:,:),
     *  b2stbr_phys_shi_bal(:),
     *  b2stbr_bas_shi0to3(:,:),
     *  b2stbr_bas_shi_bal(:),
     *  b2stbr_first_flight_shi0to3(:,:),
     *  b2stbr_first_flight_shi_bal(:),
     *  b2stbc_shi0to3(:,:),
     *  b2stbc_shi_bal(:),
     *  eirene_mc_eapl_shi0to3(:,:,:),
     *  eirene_mc_empl_shi0to3(:,:,:),
     *  eirene_mc_eipl_shi0to3(:,:,:),
     *  eirene_mc_eppl_shi0to3(:,:,:),
     *  eirene_mc_eapl_shi_bal(:,:),
     *  eirene_mc_empl_shi_bal(:,:),
     *  eirene_mc_eipl_shi_bal(:,:),
     *  eirene_mc_eppl_shi_bal(:,:),
     *  b2stbm_shi0to3(:,:),
     *  b2stbm_shi_bal(:),
     *  ext_shi0to3(:,:),
     *  ext_shi_bal(:),
     *  b2stel_shi_ion0to3(:,:),
     *  b2stel_shi_ion_bal(:),
     *  b2stel_shi_rec0to3(:,:),
     *  b2stel_shi_rec_bal(:),
     *  b2stcx_shi0to3(:,:),
     *  b2stcx_shi_bal(:),
     *  b2srsm_shi0to3(:,:),
     *  b2srsm_shi_bal(:),
     *  b2srdt_shi0to3(:,:),
     *  b2srdt_shi_bal(:),
     *  b2srst_shi0to3(:,:),
     *  b2srst_shi_bal(:),
     *  b2sihs_diaa0to3(:,:),
     *  b2sihs_diaa_bal(:),
     *  b2sihs_divua0to3(:,:),
     *  b2sihs_divua_bal(:),
     *  b2sihs_exba0to3(:,:),
     *  b2sihs_exba_bal(:),
     *  b2sihs_visa0to3(:,:),
     *  b2sihs_visa_bal(:),
     *  b2sihs_fraa0to3(:,:),
     *  b2sihs_fraa_bal(:)
*     Time-averaged quantities
      real (kind=R8), allocatable, save ::
     *  fhi_32_av(:,:),
     *  fhi_52_av(:,:),
     *  fhi_cond_av(:,:),
     *  fhi_dia_av(:,:),
     *  fhi_ecrb_av(:,:),
     *  fhi_strange_av(:,:),
     *  fhi_pschused_av(:,:),
     *  fhi_inert_av(:,:),
     *  fhi_vispar_av(:,:),
     *  fhi_visper_av(:,:),
     *  fhi_visq_av(:,:),
     *  fhi_anml_av(:,:),
     *  fhi_kevis_av(:,:),
     *  b2stbr_phys_shi_bal_av(:),
     *  b2stbr_bas_shi_bal_av(:),
     *  b2stbr_first_flight_shi_bal_av(:),
     *  b2stbc_shi_bal_av(:),
     *  eirene_mc_eapl_shi_bal_av(:,:),
     *  eirene_mc_empl_shi_bal_av(:,:),
     *  eirene_mc_eipl_shi_bal_av(:,:),
     *  eirene_mc_eppl_shi_bal_av(:,:),
     *  b2stbm_shi_bal_av(:),
     *  ext_shi_bal_av(:),
     *  b2stel_shi_ion_bal_av(:),
     *  b2stel_shi_rec_bal_av(:),
     *  b2stcx_shi_bal_av(:),
     *  b2srsm_shi_bal_av(:),
     *  b2srdt_shi_bal_av(:),
     *  b2srst_shi_bal_av(:),
     *  b2sihs_diaa_bal_av(:),
     *  b2sihs_divua_bal_av(:),
     *  b2sihs_exba_bal_av(:),
     *  b2sihs_visa_bal_av(:),
     *  b2sihs_fraa_bal_av(:),
     *  reshi_av(:)
*     Additional averaged quantities
      real (kind=R8), allocatable, save ::
     *  ne_av(:),
     *  na_av(:,:),
     *  ua_av(:,:),
     *  po_av(:),
     *  te_av(:),
     *  ti_av(:),
     *  rpt_av(:,:),
     *  kinrgy_av(:,:),
     *  fne_av(:,:),
     *  dab2_av(:,:,:),
     *  dmb2_av(:,:,:),
     *  tab2_av(:,:,:),
     *  tmb2_av(:,:,:),
     *  rfluxa_av(:,:,:),
     *  rfluxm_av(:,:,:),
     *  refluxa_av(:,:,:),
     *  refluxm_av(:,:,:),
     *  pfluxa_av(:,:,:),
     *  pfluxm_av(:,:,:),
     *  pefluxa_av(:,:,:),
     *  pefluxm_av(:,:,:)

*   ..balance_netcdf decides whether to initialise, track, and write balance arrays.
*     Default: do not create balance arrays (balance_netcdf.eq.0)
      integer, save :: balance_netcdf

*   ..balance_average decides whether to average balance arrays over time.
*     Default: do not average (balance_average.eq.0)
      integer, save :: balance_average

      contains

      subroutine alloc_b2mod_balance(nCv,nFc,nsd)
      implicit none
*   ..Input arguments (unchanged on exit)
      integer nCv,nFc,nsd

*   ..Arrays for particle balance:
*     Fluxes:
      allocate(fna_pinch(nFc,0:1,0:nsd-1))
      allocate(fna_pll(nFc,0:1,0:nsd-1))
      allocate(fna_drift(nFc,0:1,0:nsd-1))
      allocate(fna_ch(nFc,0:1,0:nsd-1))
      allocate(fna_nanom(nFc,0:1,0:nsd-1))
      allocate(fna_panom(nFc,0:1,0:nsd-1))
      allocate(fna_pschused(nFc,0:1,0:nsd-1))
*     Sources
      allocate(b2stbr_phys_sna0to1(nCv,0:1,0:nsd-1))
      allocate(b2stbr_phys_sna_bal(nCv,0:nsd-1))
      allocate(b2stbr_bas_sna0to1(nCv,0:1,0:nsd-1))
      allocate(b2stbr_bas_sna_bal(nCv,0:nsd-1))
      allocate(b2stbr_first_flight_sna0to1(nCv,0:1,0:nsd-1))
      allocate(b2stbr_first_flight_sna_bal(nCv,0:nsd-1))
      allocate(b2stbc_sna0to1(nCv,0:1,0:nsd-1))
      allocate(b2stbc_sna_bal(nCv,0:nsd-1))
      allocate(b2stbm_sna0to1(nCv,0:1,0:nsd-1))
      allocate(b2stbm_sna_bal(nCv,0:nsd-1))
      allocate(ext_sna0to1(nCv,0:1,0:nsd-1))
      allocate(ext_sna_bal(nCv,0:nsd-1))
      allocate(b2stel_sna_ion0to1(nCv,0:1,0:nsd-1))
      allocate(b2stel_sna_ion_bal(nCv,0:nsd-1))
      allocate(b2stel_sna_rec0to1(nCv,0:1,0:nsd-1))
      allocate(b2stel_sna_rec_bal(nCv,0:nsd-1))
      allocate(b2stcx_sna0to1(nCv,0:1,0:nsd-1))
      allocate(b2stcx_sna_bal(nCv,0:nsd-1))
      allocate(b2srsm_sna0to1(nCv,0:1,0:nsd-1))
      allocate(b2srsm_sna_bal(nCv,0:nsd-1))
      allocate(b2srdt_sna0to1(nCv,0:1,0:nsd-1))
      allocate(b2srdt_sna_bal(nCv,0:nsd-1))
      allocate(b2srst_sna0to1(nCv,0:1,0:nsd-1))
      allocate(b2srst_sna_bal(nCv,0:nsd-1))
      b2stbr_phys_sna0to1=0.0_R8
      b2stbr_phys_sna_bal=0.0_R8
      b2stbr_bas_sna0to1=0.0_R8
      b2stbr_bas_sna_bal=0.0_R8
      b2stbr_first_flight_sna0to1=0.0_R8
      b2stbr_first_flight_sna_bal=0.0_R8
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
*     Time-averaged quantities
      allocate(fna_pinch_av(nFc,0:1,0:nsd-1))
      allocate(fna_pll_av(nFc,0:1,0:nsd-1))
      allocate(fna_drift_av(nFc,0:1,0:nsd-1))
      allocate(fna_ch_av(nFc,0:1,0:nsd-1))
      allocate(fna_nanom_av(nFc,0:1,0:nsd-1))
      allocate(fna_panom_av(nFc,0:1,0:nsd-1))
      allocate(fna_pschused_av(nFc,0:1,0:nsd-1))
      allocate(b2stbr_phys_sna_bal_av(nCv,0:nsd-1))
      allocate(b2stbr_bas_sna_bal_av(nCv,0:nsd-1))
      allocate(b2stbr_first_flight_sna_bal_av(nCv,0:nsd-1))
      allocate(b2stbc_sna_bal_av(nCv,0:nsd-1))
      allocate(b2stbm_sna_bal_av(nCv,0:nsd-1))
      allocate(ext_sna_bal_av(nCv,0:nsd-1))
      allocate(b2stel_sna_ion_bal_av(nCv,0:nsd-1))
      allocate(b2stel_sna_rec_bal_av(nCv,0:nsd-1))
      allocate(b2stcx_sna_bal_av(nCv,0:nsd-1))
      allocate(b2srsm_sna_bal_av(nCv,0:nsd-1))
      allocate(b2srdt_sna_bal_av(nCv,0:nsd-1))
      allocate(b2srst_sna_bal_av(nCv,0:nsd-1))
      allocate(resco_av(nCv,0:nsd-1))
*   ..Arrays for momentum balance:
*     Fluxes
      allocate(fmo_flua(nFc,0:1,0:nsd-1))
      allocate(fmo_cvsa(nFc,0:1,0:nsd-1))
      allocate(fmo_b2nxfv(nFc,0:1,0:nsd-1))
*     Sources
      allocate(b2stbr_phys_smo0to3(nCv,0:3,0:nsd-1))
      allocate(b2stbr_phys_smo_bal(nCv,0:nsd-1))
      allocate(b2stbr_bas_smo0to3(nCv,0:3,0:nsd-1))
      allocate(b2stbr_bas_smo_bal(nCv,0:nsd-1))
      allocate(b2stbc_smo0to3(nCv,0:3,0:nsd-1))
      allocate(b2stbc_smo_bal(nCv,0:nsd-1))      
      allocate(b2stbm_smo0to3(nCv,0:3,0:nsd-1))
      allocate(b2stbm_smo_bal(nCv,0:nsd-1))
      allocate(ext_smo0to3(nCv,0:3,0:nsd-1))
      allocate(ext_smo_bal(nCv,0:nsd-1))
      allocate(b2stel_smq_ion0to3(nCv,0:3,0:nsd-1))
      allocate(b2stel_smq_ion_bal(nCv,0:nsd-1))
      allocate(b2stel_smq_rec0to3(nCv,0:3,0:nsd-1))
      allocate(b2stel_smq_rec_bal(nCv,0:nsd-1))
      allocate(b2stcx_smq0to3(nCv,0:3,0:nsd-1))
      allocate(b2stcx_smq_bal(nCv,0:nsd-1))
      allocate(b2srsm_smo0to3(nCv,0:3,0:nsd-1))
      allocate(b2srsm_smo_bal(nCv,0:nsd-1))
      allocate(b2srdt_smo0to3(nCv,0:3,0:nsd-1))
      allocate(b2srdt_smo_bal(nCv,0:nsd-1))
      allocate(b2srst_smo0to3(nCv,0:3,0:nsd-1))
      allocate(b2srst_smo_bal(nCv,0:nsd-1))
      allocate(b2sifr_smoch0to3(nCv,0:3,0:nsd-1))
      allocate(b2sifr_smoch_bal(nCv,0:nsd-1))
      allocate(b2sifr_smotf0to3(nCv,0:3,0:nsd-1))
      allocate(b2sifr_smotf_bal(nCv,0:nsd-1))
      allocate(b2siav_smovh0to3(nCv,0:3,0:nsd-1))
      allocate(b2siav_smovh_bal(nCv,0:nsd-1))
      allocate(b2siav_smovv0to3(nCv,0:3,0:nsd-1))
      allocate(b2siav_smovv_bal(nCv,0:nsd-1))
      allocate(b2sicf_smo0to3(nCv,0:3,0:nsd-1))
      allocate(b2sicf_smo_bal(nCv,0:nsd-1))
      allocate(b2sian_smo0to3(nCv,0:3,0:nsd-1))              !srv 25.09.17
      allocate(b2sian_smo_bal(nCv,0:nsd-1))                  !srv 25.09.17
      allocate(b2nxdv_smo0to3(nCv,0:3,0:nsd-1))
      allocate(b2nxdv_smo_bal(nCv,0:nsd-1))
      allocate(b2sigp_smogp0to3(nCv,0:3,0:nsd-1))
      allocate(b2sigp_smogp_bal(nCv,0:nsd-1))
      allocate(b2sigp_pstat0to3(nCv,0:3,0:nsd-1))
      allocate(b2sigp_pstat_bal(nCv,0:nsd-1))
      b2stbr_phys_smo0to3=0.0_R8
      b2stbr_phys_smo_bal=0.0_R8
      b2stbr_bas_smo0to3=0.0_R8
      b2stbr_bas_smo_bal=0.0_R8
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
      b2sian_smo0to3=0.0_R8                    !srv 25.09.17
      b2sian_smo_bal=0.0_R8                    !srv 25.09.17
      b2nxdv_smo0to3=0.0_R8
      b2nxdv_smo_bal=0.0_R8
      b2sigp_smogp0to3=0.0_R8
      b2sigp_smogp_bal=0.0_R8
      b2sigp_pstat0to3=0.0_R8
      b2sigp_pstat_bal=0.0_R8
*     Time-averaged quantities
      allocate(fmo_flua_av(nFc,0:1,0:nsd-1))
      allocate(fmo_cvsa_av(nFc,0:1,0:nsd-1))
      allocate(fmo_b2nxfv_av(nFc,0:1,0:nsd-1))
      allocate(b2stbr_phys_smo_bal_av(nCv,0:nsd-1))
      allocate(b2stbr_bas_smo_bal_av(nCv,0:nsd-1))
      allocate(b2stbc_smo_bal_av(nCv,0:nsd-1))
      allocate(b2stbm_smo_bal_av(nCv,0:nsd-1))
      allocate(ext_smo_bal_av(nCv,0:nsd-1))
      allocate(b2stel_smq_ion_bal_av(nCv,0:nsd-1))
      allocate(b2stel_smq_rec_bal_av(nCv,0:nsd-1))
      allocate(b2stcx_smq_bal_av(nCv,0:nsd-1))
      allocate(b2srsm_smo_bal_av(nCv,0:nsd-1))
      allocate(b2srdt_smo_bal_av(nCv,0:nsd-1))
      allocate(b2srst_smo_bal_av(nCv,0:nsd-1))
      allocate(b2sifr_smoch_bal_av(nCv,0:nsd-1))
      allocate(b2sifr_smotf_bal_av(nCv,0:nsd-1))
      allocate(b2siav_smovh_bal_av(nCv,0:nsd-1))
      allocate(b2siav_smovv_bal_av(nCv,0:nsd-1))
      allocate(b2sicf_smo_bal_av(nCv,0:nsd-1))
      allocate(b2sian_smo_bal_av(nCv,0:nsd-1))          !srv 25.09.17
      allocate(b2nxdv_smo_bal_av(nCv,0:nsd-1))
      allocate(b2sigp_smogp_bal_av(nCv,0:nsd-1))
      allocate(b2sigp_pstat_bal_av(nCv,0:nsd-1))
      allocate(resmo_av(nCv,0:nsd-1))
*   ..Arrays for electron heat balance:
*     Fluxes
      allocate(fhe_32(nFc,0:1))
      allocate(fhe_52(nFc,0:1))
      allocate(fhe_thermj(nFc,0:1))
      allocate(fhe_cond(nFc,0:1))
      allocate(fhe_dia(nFc,0:1))
      allocate(fhe_ecrb(nFc,0:1))
      allocate(fhe_strange(nFc,0:1))
      allocate(fhe_pschused(nFc,0:1))
*     Sources
      allocate(b2stbr_phys_she0to3(nCv,0:3))
      allocate(b2stbr_phys_she_bal(nCv))
      allocate(b2stbr_bas_she0to3(nCv,0:3))
      allocate(b2stbr_bas_she_bal(nCv))
      allocate(b2stbr_first_flight_she0to3(nCv,0:3))
      allocate(b2stbr_first_flight_she_bal(nCv))
      allocate(b2stbc_she0to3(nCv,0:3))
      allocate(b2stbc_she_bal(nCv))
      allocate(b2stbm_she0to3(nCv,0:3))
      allocate(b2stbm_she_bal(nCv))
      allocate(ext_she0to3(nCv,0:3))
      allocate(ext_she_bal(nCv))
      allocate(b2stel_she0to3(nCv,0:3))
      allocate(b2stel_she_bal(nCv))
      allocate(b2srsm_she0to3(nCv,0:3))
      allocate(b2srsm_she_bal(nCv))
      allocate(b2srdt_she0to3(nCv,0:3))
      allocate(b2srdt_she_bal(nCv))
      allocate(b2srst_she0to3(nCv,0:3))
      allocate(b2srst_she_bal(nCv))
      allocate(b2sihs_diae0to3(nCv,0:3))
      allocate(b2sihs_diae_bal(nCv))
      allocate(b2sihs_divue0to3(nCv,0:3))
      allocate(b2sihs_divue_bal(nCv))
      allocate(b2sihs_exbe0to3(nCv,0:3))
      allocate(b2sihs_exbe_bal(nCv))
      allocate(b2sihs_joule0to3(nCv,0:3))
      allocate(b2sihs_joule_bal(nCv))
      allocate(b2npht_shei_bal(nCv))
      b2stbr_phys_she0to3=0.0_R8
      b2stbr_phys_she_bal=0.0_R8
      b2stbr_bas_she0to3=0.0_R8
      b2stbr_bas_she_bal=0.0_R8
      b2stbr_first_flight_she0to3=0.0_R8
      b2stbr_first_flight_she_bal=0.0_R8
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
*     Time-averaged quantities
      allocate(fhe_32_av(nFc,0:1))
      allocate(fhe_52_av(nFc,0:1))
      allocate(fhe_thermj_av(nFc,0:1))
      allocate(fhe_cond_av(nFc,0:1))
      allocate(fhe_dia_av(nFc,0:1))
      allocate(fhe_ecrb_av(nFc,0:1))
      allocate(fhe_strange_av(nFc,0:1))
      allocate(fhe_pschused_av(nFc,0:1))
      allocate(b2stbr_phys_she_bal_av(nCv))
      allocate(b2stbr_bas_she_bal_av(nCv))
      allocate(b2stbr_first_flight_she_bal_av(nCv))
      allocate(b2stbc_she_bal_av(nCv))
      allocate(b2stbm_she_bal_av(nCv))
      allocate(ext_she_bal_av(nCv))
      allocate(b2stel_she_bal_av(nCv))
      allocate(b2srsm_she_bal_av(nCv))
      allocate(b2srdt_she_bal_av(nCv))
      allocate(b2srst_she_bal_av(nCv))
      allocate(b2sihs_diae_bal_av(nCv))
      allocate(b2sihs_divue_bal_av(nCv))
      allocate(b2sihs_exbe_bal_av(nCv))
      allocate(b2sihs_joule_bal_av(nCv))
      allocate(b2npht_shei_bal_av(nCv))
      allocate(reshe_av(nCv))
*   ..Arrays for ion heat balance:
*     Fluxes
      allocate(fhi_32(nFc,0:1))
      allocate(fhi_52(nFc,0:1))
      allocate(fhi_cond(nFc,0:1))
      allocate(fhi_dia(nFc,0:1))
      allocate(fhi_ecrb(nFc,0:1))
      allocate(fhi_strange(nFc,0:1))
      allocate(fhi_pschused(nFc,0:1))
      allocate(fhi_inert(nFc,0:1))
      allocate(fhi_vispar(nFc,0:1))
      allocate(fhi_visper(nFc,0:1))
      allocate(fhi_visq(nFc,0:1))
      allocate(fhi_anml(nFc,0:1))
      allocate(fhi_kevis(nFc,0:1))
*     Sources
      allocate(b2stbr_phys_shi0to3(nCv,0:3))
      allocate(b2stbr_phys_shi_bal(nCv))
      allocate(b2stbr_bas_shi0to3(nCv,0:3))
      allocate(b2stbr_bas_shi_bal(nCv))
      allocate(b2stbr_first_flight_shi0to3(nCv,0:3))
      allocate(b2stbr_first_flight_shi_bal(nCv))
      allocate(b2stbc_shi0to3(nCv,0:3))
      allocate(b2stbc_shi_bal(nCv))
      allocate(b2stbm_shi0to3(nCv,0:3))
      allocate(b2stbm_shi_bal(nCv))
      allocate(ext_shi0to3(nCv,0:3))
      allocate(ext_shi_bal(nCv))
      allocate(b2stel_shi_ion0to3(nCv,0:3))
      allocate(b2stel_shi_ion_bal(nCv))
      allocate(b2stel_shi_rec0to3(nCv,0:3))
      allocate(b2stel_shi_rec_bal(nCv))
      allocate(b2stcx_shi0to3(nCv,0:3))
      allocate(b2stcx_shi_bal(nCv))
      allocate(b2srsm_shi0to3(nCv,0:3))
      allocate(b2srsm_shi_bal(nCv))
      allocate(b2srdt_shi0to3(nCv,0:3))
      allocate(b2srdt_shi_bal(nCv))
      allocate(b2srst_shi0to3(nCv,0:3))
      allocate(b2srst_shi_bal(nCv))
      allocate(b2sihs_diaa0to3(nCv,0:3))
      allocate(b2sihs_diaa_bal(nCv))
      allocate(b2sihs_divua0to3(nCv,0:3))
      allocate(b2sihs_divua_bal(nCv))
      allocate(b2sihs_exba0to3(nCv,0:3))
      allocate(b2sihs_exba_bal(nCv))
      allocate(b2sihs_visa0to3(nCv,0:3))
      allocate(b2sihs_visa_bal(nCv))
      allocate(b2sihs_fraa0to3(nCv,0:3))
      allocate(b2sihs_fraa_bal(nCv))
      b2stbr_phys_shi0to3=0.0_R8
      b2stbr_phys_shi_bal=0.0_R8
      b2stbr_bas_shi0to3=0.0_R8
      b2stbr_bas_shi_bal=0.0_R8
      b2stbr_first_flight_shi0to3=0.0_R8
      b2stbr_first_flight_shi_bal=0.0_R8
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
*     Time-averaged quantities
      allocate(fhi_32_av(nFc,0:1))
      allocate(fhi_52_av(nFc,0:1))
      allocate(fhi_cond_av(nFc,0:1))
      allocate(fhi_dia_av(nFc,0:1))
      allocate(fhi_ecrb_av(nFc,0:1))
      allocate(fhi_strange_av(nFc,0:1))
      allocate(fhi_pschused_av(nFc,0:1))
      allocate(fhi_inert_av(nFc,0:1))
      allocate(fhi_vispar_av(nFc,0:1))
      allocate(fhi_visper_av(nFc,0:1))
      allocate(fhi_visq_av(nFc,0:1))
      allocate(fhi_anml_av(nFc,0:1))
      allocate(fhi_kevis_av(nFc,0:1))
      allocate(b2stbr_phys_shi_bal_av(nCv))
      allocate(b2stbr_bas_shi_bal_av(nCv))
      allocate(b2stbr_first_flight_shi_bal_av(nCv))
      allocate(b2stbc_shi_bal_av(nCv))
      allocate(b2stbm_shi_bal_av(nCv))
      allocate(ext_shi_bal_av(nCv))
      allocate(b2stel_shi_ion_bal_av(nCv))
      allocate(b2stel_shi_rec_bal_av(nCv))
      allocate(b2stcx_shi_bal_av(nCv))
      allocate(b2srsm_shi_bal_av(nCv))
      allocate(b2srdt_shi_bal_av(nCv))
      allocate(b2srst_shi_bal_av(nCv))
      allocate(b2sihs_diaa_bal_av(nCv))
      allocate(b2sihs_divua_bal_av(nCv))
      allocate(b2sihs_exba_bal_av(nCv))
      allocate(b2sihs_visa_bal_av(nCv))
      allocate(b2sihs_fraa_bal_av(nCv))
      allocate(reshi_av(nCv))
*     Additional averaged quantities
      allocate(ne_av(nCv))
      allocate(na_av(nCv,0:nsd-1))
      allocate(ua_av(nCv,0:nsd-1))
      allocate(po_av(nCv))
      allocate(te_av(nCv))
      allocate(ti_av(nCv))
      allocate(rpt_av(nCv,0:nsd-1))
      allocate(kinrgy_av(nCv,0:nsd-1))
      allocate(fne_av(nFc,0:1))

      return
      end subroutine alloc_b2mod_balance

      subroutine dealloc_b2mod_balance
      implicit none

      if(.not.allocated(fna_pinch)) return
cdjm Jan2017 Flux arrays:
*   ..Arrays for particle balance:
*     Fluxes
      deallocate(fna_pinch)
      deallocate(fna_pll)
      deallocate(fna_drift)
      deallocate(fna_ch)
      deallocate(fna_nanom)
      deallocate(fna_panom)
      deallocate(fna_pschused)
*     Sources
      deallocate(b2stbr_phys_sna0to1)
      deallocate(b2stbr_phys_sna_bal)
      deallocate(b2stbr_bas_sna0to1)
      deallocate(b2stbr_bas_sna_bal)
      deallocate(b2stbr_first_flight_sna0to1)
      deallocate(b2stbr_first_flight_sna_bal)
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
*     Time-averaged quantities
      deallocate(fna_pinch_av)
      deallocate(fna_pll_av)
      deallocate(fna_drift_av)
      deallocate(fna_ch_av)
      deallocate(fna_nanom_av)
      deallocate(fna_panom_av)
      deallocate(fna_pschused_av)
      deallocate(b2stbr_phys_sna_bal_av)
      deallocate(b2stbr_bas_sna_bal_av)
      deallocate(b2stbr_first_flight_sna_bal_av)
      deallocate(b2stbc_sna_bal_av)
      deallocate(b2stbm_sna_bal_av)
      deallocate(ext_sna_bal_av)
      deallocate(b2stel_sna_ion_bal_av)
      deallocate(b2stel_sna_rec_bal_av)
      deallocate(b2stcx_sna_bal_av)
      deallocate(b2srsm_sna_bal_av)
      deallocate(b2srdt_sna_bal_av)
      deallocate(b2srst_sna_bal_av)
      deallocate(resco_av)
*   ..Arrays for momentum balance:
*     Fluxes
      deallocate(fmo_flua)
      deallocate(fmo_cvsa)
      deallocate(fmo_b2nxfv)
*     Sources
      deallocate(b2stbr_phys_smo0to3)
      deallocate(b2stbr_phys_smo_bal)
      deallocate(b2stbr_bas_smo0to3)
      deallocate(b2stbr_bas_smo_bal)
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
      deallocate(b2sian_smo0to3)              !srv 25.09.17
      deallocate(b2sian_smo_bal)              !srv 25.09.17
      deallocate(b2nxdv_smo0to3)
      deallocate(b2nxdv_smo_bal)
      deallocate(b2sigp_smogp0to3)
      deallocate(b2sigp_smogp_bal)
      deallocate(b2sigp_pstat0to3)
      deallocate(b2sigp_pstat_bal)
*     Time-averaged quantities
      deallocate(fmo_flua_av)
      deallocate(fmo_cvsa_av)
      deallocate(fmo_b2nxfv_av)
      deallocate(b2stbr_phys_smo_bal_av)
      deallocate(b2stbr_bas_smo_bal_av)
      deallocate(b2stbc_smo_bal_av)
      deallocate(b2stbm_smo_bal_av)
      deallocate(ext_smo_bal_av)
      deallocate(b2stel_smq_ion_bal_av)
      deallocate(b2stel_smq_rec_bal_av)
      deallocate(b2stcx_smq_bal_av)
      deallocate(b2srsm_smo_bal_av)
      deallocate(b2srdt_smo_bal_av)
      deallocate(b2srst_smo_bal_av)
      deallocate(b2sifr_smoch_bal_av)
      deallocate(b2sifr_smotf_bal_av)
      deallocate(b2siav_smovh_bal_av)
      deallocate(b2siav_smovv_bal_av)
      deallocate(b2sicf_smo_bal_av)
      deallocate(b2sian_smo_bal_av)              !srv 25.09.17
      deallocate(b2nxdv_smo_bal_av)
      deallocate(b2sigp_smogp_bal_av)
      deallocate(b2sigp_pstat_bal_av)
      deallocate(resmo_av)
*   ..Arrays for electron heat balance:
*     Fluxes
      deallocate(fhe_32)
      deallocate(fhe_52)
      deallocate(fhe_thermj)
      deallocate(fhe_cond)
      deallocate(fhe_dia)
      deallocate(fhe_ecrb)
      deallocate(fhe_strange)
      deallocate(fhe_pschused)
*     Sources
      deallocate(b2stbr_phys_she0to3)
      deallocate(b2stbr_phys_she_bal)
      deallocate(b2stbr_bas_she0to3)
      deallocate(b2stbr_bas_she_bal)
      deallocate(b2stbr_first_flight_she0to3)
      deallocate(b2stbr_first_flight_she_bal)
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
*     Time-averaged quantities
      deallocate(fhe_32_av)
      deallocate(fhe_52_av)
      deallocate(fhe_thermj_av)
      deallocate(fhe_cond_av)
      deallocate(fhe_dia_av)
      deallocate(fhe_ecrb_av)
      deallocate(fhe_strange_av)
      deallocate(fhe_pschused_av)
      deallocate(b2stbr_phys_she_bal_av)
      deallocate(b2stbr_bas_she_bal_av)
      deallocate(b2stbr_first_flight_she_bal_av)
      deallocate(b2stbc_she_bal_av)
      deallocate(b2stbm_she_bal_av)
      deallocate(ext_she_bal_av)
      deallocate(b2stel_she_bal_av)
      deallocate(b2srsm_she_bal_av)
      deallocate(b2srdt_she_bal_av)
      deallocate(b2srst_she_bal_av)
      deallocate(b2sihs_diae_bal_av)
      deallocate(b2sihs_divue_bal_av)
      deallocate(b2sihs_exbe_bal_av)
      deallocate(b2sihs_joule_bal_av)
      deallocate(b2npht_shei_bal_av)
      deallocate(reshe_av)
*   ..Arrays for ion heat balance:
*     Fluxes
      deallocate(fhi_32)
      deallocate(fhi_52)
      deallocate(fhi_cond)
      deallocate(fhi_dia)
      deallocate(fhi_ecrb)
      deallocate(fhi_strange)
      deallocate(fhi_pschused)
      deallocate(fhi_inert)
      deallocate(fhi_vispar)
      deallocate(fhi_visper)
      deallocate(fhi_visq)
      deallocate(fhi_anml)
      deallocate(fhi_kevis)
*     Sources
      deallocate(b2stbr_phys_shi0to3)
      deallocate(b2stbr_phys_shi_bal)
      deallocate(b2stbr_bas_shi0to3)
      deallocate(b2stbr_bas_shi_bal)
      deallocate(b2stbr_first_flight_shi0to3)
      deallocate(b2stbr_first_flight_shi_bal)
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
*     Time-averaged quantities
      deallocate(fhi_32_av)
      deallocate(fhi_52_av)
      deallocate(fhi_cond_av)
      deallocate(fhi_dia_av)
      deallocate(fhi_ecrb_av)
      deallocate(fhi_strange_av)
      deallocate(fhi_pschused_av)
      deallocate(fhi_inert_av)
      deallocate(fhi_vispar_av)
      deallocate(fhi_visper_av)
      deallocate(fhi_visq_av)
      deallocate(fhi_anml_av)
      deallocate(fhi_kevis_av)
      deallocate(b2stbr_phys_shi_bal_av)
      deallocate(b2stbr_bas_shi_bal_av)
      deallocate(b2stbr_first_flight_shi_bal_av)
      deallocate(b2stbc_shi_bal_av)
      deallocate(b2stbm_shi_bal_av)
      deallocate(ext_shi_bal_av)
      deallocate(b2stel_shi_ion_bal_av)
      deallocate(b2stel_shi_rec_bal_av)
      deallocate(b2stcx_shi_bal_av)
      deallocate(b2srsm_shi_bal_av)
      deallocate(b2srdt_shi_bal_av)
      deallocate(b2srst_shi_bal_av)
      deallocate(b2sihs_diaa_bal_av)
      deallocate(b2sihs_divua_bal_av)
      deallocate(b2sihs_exba_bal_av)
      deallocate(b2sihs_visa_bal_av)
      deallocate(b2sihs_fraa_bal_av)
      deallocate(reshi_av)
*     Additional averaged quantities
      deallocate(ne_av)
      deallocate(na_av)
      deallocate(ua_av)
      deallocate(po_av)
      deallocate(te_av)
      deallocate(ti_av)

      return
      end subroutine dealloc_b2mod_balance

      subroutine alloc_b2mod_eirene_sources(nCv,nFc,nsd,nstrad)
#ifdef B25_EIRENE
#ifndef SOLPS4_3
      use eirmod_extrab25
     & , only : dab2,dmb2,tab2,tmb2,rfluxa,rfluxm,refluxa,refluxm,
     &          pfluxa,pfluxm,pefluxa,pefluxm
#else
      use eirdiag ! IGNORE
     & , only : dab2,dmb2,tab2,tmb2,rfluxa,rfluxm,refluxa,refluxm,
     &          pfluxa,pfluxm,pefluxa,pefluxm
#endif
#else
      use b2mod_eirdiag
#endif
      implicit none
      integer nCv,nFc,nsd,nstrad
*   ..procedures
      intrinsic size

cdjm Jan2017 Eirene arrays for particle balance:
      allocate(
     & eirene_mc_papl_sna0to1(nCv,0:1,0:nsd-1,1:nstrad))
      allocate(
     & eirene_mc_pmpl_sna0to1(nCv,0:1,0:nsd-1,1:nstrad))
      allocate(
     & eirene_mc_pipl_sna0to1(nCv,0:1,0:nsd-1,1:nstrad))
      allocate(
     & eirene_mc_pppl_sna0to1(nCv,0:1,0:nsd-1,1:nstrad))
      allocate(
     & eirene_mc_core_sna0to1(nCv,0:1,0:nsd-1))
      allocate(eirene_mc_papl_sna_bal(nCv,0:nsd-1,1:nstrad))
      allocate(eirene_mc_pmpl_sna_bal(nCv,0:nsd-1,1:nstrad))
      allocate(eirene_mc_pipl_sna_bal(nCv,0:nsd-1,1:nstrad))
      allocate(eirene_mc_pppl_sna_bal(nCv,0:nsd-1,1:nstrad))
      allocate(eirene_mc_core_sna_bal(nCv,0:nsd-1))
      eirene_mc_papl_sna0to1=0.0_R8
      eirene_mc_pmpl_sna0to1=0.0_R8
      eirene_mc_pipl_sna0to1=0.0_R8
      eirene_mc_pppl_sna0to1=0.0_R8
      eirene_mc_core_sna0to1=0.0_R8
      eirene_mc_papl_sna_bal=0.0_R8
      eirene_mc_pmpl_sna_bal=0.0_R8
      eirene_mc_pipl_sna_bal=0.0_R8
      eirene_mc_pppl_sna_bal=0.0_R8
      eirene_mc_core_sna_bal=0.0_R8
cdjm Jun2017 Time-averaged arrays for particle balance:
      allocate(
     & eirene_mc_papl_sna_bal_av(nCv,0:nsd-1,1:nstrad))
      allocate(
     & eirene_mc_pmpl_sna_bal_av(nCv,0:nsd-1,1:nstrad))
      allocate(
     & eirene_mc_pipl_sna_bal_av(nCv,0:nsd-1,1:nstrad))
      allocate(
     & eirene_mc_pppl_sna_bal_av(nCv,0:nsd-1,1:nstrad))
      allocate(
     & eirene_mc_core_sna_bal_av(nCv,0:nsd-1))
cdjm Jan2017 Eirene arrays for momentum balance:
      allocate(
     & eirene_mc_mapl_smo0to3(nCv,0:3,0:nsd-1,1:nstrad))
      allocate(
     & eirene_mc_mmpl_smo0to3(nCv,0:3,0:nsd-1,1:nstrad))
      allocate(
     & eirene_mc_mipl_smo0to3(nCv,0:3,0:nsd-1,1:nstrad))
      allocate(
     & eirene_mc_cppv_smo0to3(nCv,0:3,0:nsd-1,1:nstrad))
      allocate(eirene_mc_mapl_smo_bal(nCv,0:nsd-1,1:nstrad))
      allocate(eirene_mc_mmpl_smo_bal(nCv,0:nsd-1,1:nstrad))
      allocate(eirene_mc_mipl_smo_bal(nCv,0:nsd-1,1:nstrad))
      allocate(eirene_mc_cppv_smo_bal(nCv,0:nsd-1,1:nstrad))
      eirene_mc_mapl_smo0to3=0.0_R8
      eirene_mc_mmpl_smo0to3=0.0_R8
      eirene_mc_mipl_smo0to3=0.0_R8
      eirene_mc_cppv_smo0to3=0.0_R8
      eirene_mc_mapl_smo_bal=0.0_R8
      eirene_mc_mmpl_smo_bal=0.0_R8
      eirene_mc_mipl_smo_bal=0.0_R8
      eirene_mc_cppv_smo_bal=0.0_R8
cdjm Jun2017 Time-averaged arrays for momentum balance:
      allocate(
     & eirene_mc_mapl_smo_bal_av(nCv,0:nsd-1,1:nstrad))
      allocate(
     & eirene_mc_mmpl_smo_bal_av(nCv,0:nsd-1,1:nstrad))
      allocate(
     & eirene_mc_mipl_smo_bal_av(nCv,0:nsd-1,1:nstrad))
      allocate(
     & eirene_mc_cppv_smo_bal_av(nCv,0:nsd-1,1:nstrad))
cdjm Jan2017 Eirene arrays for electron heat balance:
      allocate(eirene_mc_eael_she0to3(nCv,0:3,1:nstrad))
      allocate(eirene_mc_emel_she0to3(nCv,0:3,1:nstrad))
      allocate(eirene_mc_eiel_she0to3(nCv,0:3,1:nstrad))
      allocate(eirene_mc_epel_she0to3(nCv,0:3,1:nstrad))
      allocate(eirene_mc_eael_she_bal(nCv,1:nstrad))
      allocate(eirene_mc_emel_she_bal(nCv,1:nstrad))
      allocate(eirene_mc_eiel_she_bal(nCv,1:nstrad))
      allocate(eirene_mc_epel_she_bal(nCv,1:nstrad))
      eirene_mc_eael_she0to3=0.0_R8
      eirene_mc_emel_she0to3=0.0_R8
      eirene_mc_eiel_she0to3=0.0_R8
      eirene_mc_epel_she0to3=0.0_R8
      eirene_mc_eael_she_bal=0.0_R8
      eirene_mc_emel_she_bal=0.0_R8
      eirene_mc_eiel_she_bal=0.0_R8
      eirene_mc_epel_she_bal=0.0_R8      
cdjm Jun2017 Time-averaged arrays for electron heat balance:
      allocate(eirene_mc_eael_she_bal_av(nCv,1:nstrad))
      allocate(eirene_mc_emel_she_bal_av(nCv,1:nstrad))
      allocate(eirene_mc_eiel_she_bal_av(nCv,1:nstrad))
      allocate(eirene_mc_epel_she_bal_av(nCv,1:nstrad))
cdjm Jan2017 Eirene arrays for ion heat balance:
      allocate(eirene_mc_eapl_shi0to3(nCv,0:3,1:nstrad))
      allocate(eirene_mc_empl_shi0to3(nCv,0:3,1:nstrad))
      allocate(eirene_mc_eipl_shi0to3(nCv,0:3,1:nstrad))
      allocate(eirene_mc_eppl_shi0to3(nCv,0:3,1:nstrad))
      allocate(eirene_mc_eapl_shi_bal(nCv,1:nstrad))
      allocate(eirene_mc_empl_shi_bal(nCv,1:nstrad))
      allocate(eirene_mc_eipl_shi_bal(nCv,1:nstrad))
      allocate(eirene_mc_eppl_shi_bal(nCv,1:nstrad))
      eirene_mc_eapl_shi0to3=0.0_R8
      eirene_mc_empl_shi0to3=0.0_R8
      eirene_mc_eipl_shi0to3=0.0_R8
      eirene_mc_eppl_shi0to3=0.0_R8
      eirene_mc_eapl_shi_bal=0.0_R8
      eirene_mc_empl_shi_bal=0.0_R8
      eirene_mc_eipl_shi_bal=0.0_R8
      eirene_mc_eppl_shi_bal=0.0_R8
cdjm Jun2017 Time-averaged arrays for ion heat balance:
      allocate(eirene_mc_eapl_shi_bal_av(nCv,1:nstrad))
      allocate(eirene_mc_empl_shi_bal_av(nCv,1:nstrad))
      allocate(eirene_mc_eipl_shi_bal_av(nCv,1:nstrad))
      allocate(eirene_mc_eppl_shi_bal_av(nCv,1:nstrad))
cdjm Jun2017 Time-averaged arrays for fluxes and densities:
      allocate(dab2_av(size(dab2,1),size(dab2,2),size(dab2,3)))
      allocate(dmb2_av(size(dmb2,1),size(dmb2,2),size(dmb2,3)))
      allocate(tab2_av(size(tab2,1),size(tab2,2),size(tab2,3)))
      allocate(tmb2_av(size(tmb2,1),size(tmb2,2),size(tmb2,3)))
      allocate(rfluxa_av(size(rfluxa,1),size(rfluxa,2),size(rfluxa,3)))
      allocate(rfluxm_av(size(rfluxm,1),size(rfluxm,2),size(rfluxm,3)))
      allocate(refluxa_av(size(refluxa,1),size(refluxa,2),
     &                    size(refluxa,3)))
      allocate(refluxm_av(size(refluxm,1),size(refluxm,2),
     &                    size(refluxm,3)))
      allocate(pfluxa_av(size(pfluxa,1),size(pfluxa,2),size(pfluxa,3)))
      allocate(pfluxm_av(size(pfluxm,1),size(pfluxm,2),size(pfluxm,3)))
      allocate(pefluxa_av(size(pefluxa,1),size(pefluxa,2),
     &                    size(pefluxa,3)))
      allocate(pefluxm_av(size(pefluxm,1),size(pefluxm,2),
     &                    size(pefluxm,3)))

      dab2_av=0.0_R8
      dmb2_av=0.0_R8
      tab2_av=0.0_R8
      tmb2_av=0.0_R8
      rfluxa_av=0.0_R8
      rfluxm_av=0.0_R8
      refluxa_av=0.0_R8
      refluxm_av=0.0_R8
      pfluxa_av=0.0_R8
      pfluxm_av=0.0_R8
      pefluxa_av=0.0_R8
      pefluxm_av=0.0_R8

      return
      end subroutine alloc_b2mod_eirene_sources

      subroutine dealloc_b2mod_eirene_sources
      implicit none

      if(.not.allocated(eirene_mc_papl_sna_bal)) return
cdjm Jan2017 Eirene arrays for particle balance:
      deallocate(eirene_mc_papl_sna0to1)
      deallocate(eirene_mc_pmpl_sna0to1)
      deallocate(eirene_mc_pipl_sna0to1)
      deallocate(eirene_mc_pppl_sna0to1)
      deallocate(eirene_mc_core_sna0to1)
      deallocate(eirene_mc_papl_sna_bal)
      deallocate(eirene_mc_pmpl_sna_bal)
      deallocate(eirene_mc_pipl_sna_bal)
      deallocate(eirene_mc_pppl_sna_bal)
      deallocate(eirene_mc_core_sna_bal)
      deallocate(eirene_mc_papl_sna_bal_av)
      deallocate(eirene_mc_pmpl_sna_bal_av)
      deallocate(eirene_mc_pipl_sna_bal_av)
      deallocate(eirene_mc_pppl_sna_bal_av)
      deallocate(eirene_mc_core_sna_bal_av)
cdjm Jan2017 Eirene arrays for momentum balance:
      deallocate(eirene_mc_mapl_smo0to3)
      deallocate(eirene_mc_mmpl_smo0to3)
      deallocate(eirene_mc_mipl_smo0to3)
      deallocate(eirene_mc_cppv_smo0to3)
      deallocate(eirene_mc_mapl_smo_bal)
      deallocate(eirene_mc_mmpl_smo_bal)
      deallocate(eirene_mc_mipl_smo_bal)
      deallocate(eirene_mc_cppv_smo_bal)
cdjm Jan2017 Eirene arrays for electron heat balance:
      deallocate(eirene_mc_eael_she0to3)
      deallocate(eirene_mc_emel_she0to3)
      deallocate(eirene_mc_eiel_she0to3)
      deallocate(eirene_mc_epel_she0to3)
      deallocate(eirene_mc_eael_she_bal)
      deallocate(eirene_mc_emel_she_bal)
      deallocate(eirene_mc_eiel_she_bal)
      deallocate(eirene_mc_epel_she_bal)
cdjm Jan2017 Eirene arrays for ion heat balance:
      deallocate(eirene_mc_eapl_shi0to3)
      deallocate(eirene_mc_empl_shi0to3)
      deallocate(eirene_mc_eipl_shi0to3)
      deallocate(eirene_mc_eppl_shi0to3)
      deallocate(eirene_mc_eapl_shi_bal)
      deallocate(eirene_mc_empl_shi_bal)
      deallocate(eirene_mc_eipl_shi_bal)
      deallocate(eirene_mc_eppl_shi_bal)
cdjm Jen2017 Time-averaged arrays for fluxes and densities:
      deallocate(dab2_av)
      deallocate(dmb2_av)
      deallocate(tab2_av)
      deallocate(tmb2_av)
      deallocate(rfluxa_av)
      deallocate(rfluxm_av)
      deallocate(refluxa_av)
      deallocate(refluxm_av)
      deallocate(pfluxa_av)
      deallocate(pfluxm_av)
      deallocate(pefluxa_av)
      deallocate(pefluxm_av)

      return
      end subroutine dealloc_b2mod_eirene_sources

      subroutine update_balance_sources_particles(nCv,isb,nb) !djm Jan2017
      implicit none
*   ..input arguments (unchanged on exit)
      integer nCv,isb
      real (kind=R8):: nb(nCv)
*   ..local variables
      integer iCv,istrai,nstrad
*   ..procedures
      intrinsic size

      nstrad = size(eirene_mc_papl_sna_bal,3)
      
      do iCv = 1, nCv
          b2stbr_phys_sna_bal(iCv,isb)=
     &                b2stbr_phys_sna0to1(iCv,0,isb)+
     &                b2stbr_phys_sna0to1(iCv,1,isb)*nb(iCv)
          b2stbr_bas_sna_bal(iCv,isb)=
     &                b2stbr_bas_sna0to1(iCv,0,isb)+
     &                b2stbr_bas_sna0to1(iCv,1,isb)*nb(iCv)
          b2stbr_first_flight_sna_bal(iCv,isb)=
     &                b2stbr_first_flight_sna0to1(iCv,0,isb)+
     &                b2stbr_first_flight_sna0to1(iCv,1,isb)*nb(iCv)
          b2stbc_sna_bal(iCv,isb)=b2stbc_sna0to1(iCv,0,isb)+
     &                     b2stbc_sna0to1(iCv,1,isb)*nb(iCv)
          if (allocated(eirene_mc_papl_sna_bal)) then
            do istrai=1,nstrad
              eirene_mc_papl_sna_bal(iCv,isb,istrai)=
     &         eirene_mc_papl_sna0to1(iCv,0,isb,istrai)+
     &         eirene_mc_papl_sna0to1(iCv,1,isb,istrai)*nb(iCv)
              eirene_mc_pmpl_sna_bal(iCv,isb,istrai)=
     &         eirene_mc_pmpl_sna0to1(iCv,0,isb,istrai)+
     &         eirene_mc_pmpl_sna0to1(iCv,1,isb,istrai)*nb(iCv)
              eirene_mc_pipl_sna_bal(iCv,isb,istrai)=
     &         eirene_mc_pipl_sna0to1(iCv,0,isb,istrai)+
     &         eirene_mc_pipl_sna0to1(iCv,1,isb,istrai)*nb(iCv)
              eirene_mc_pppl_sna_bal(iCv,isb,istrai)=
     &         eirene_mc_pppl_sna0to1(iCv,0,isb,istrai)+
     &         eirene_mc_pppl_sna0to1(iCv,1,isb,istrai)*nb(iCv)
            enddo
            eirene_mc_core_sna_bal(iCv,isb)=
     &       eirene_mc_core_sna0to1(iCv,0,isb)+
     &       eirene_mc_core_sna0to1(iCv,1,isb)*nb(iCv)
          endif
          b2stbm_sna_bal(iCv,isb)=b2stbm_sna0to1(iCv,0,isb)+
     &                     b2stbm_sna0to1(iCv,1,isb)*nb(iCv)
          ext_sna_bal(iCv,isb)=ext_sna0to1(iCv,0,isb)+
     &                     ext_sna0to1(iCv,1,isb)*nb(iCv)
          b2stel_sna_ion_bal(iCv,isb)=b2stel_sna_ion0to1(iCv,0,isb)+
     &                     b2stel_sna_ion0to1(iCv,1,isb)*nb(iCv)
          b2stel_sna_rec_bal(iCv,isb)=b2stel_sna_rec0to1(iCv,0,isb)+
     &                     b2stel_sna_rec0to1(iCv,1,isb)*nb(iCv)
          b2stcx_sna_bal(iCv,isb)=b2stcx_sna0to1(iCv,0,isb)+
     &                     b2stcx_sna0to1(iCv,1,isb)*nb(iCv)
          b2srsm_sna_bal(iCv,isb)=b2srsm_sna0to1(iCv,0,isb)+
     &                     b2srsm_sna0to1(iCv,1,isb)*nb(iCv)
          b2srdt_sna_bal(iCv,isb)=b2srdt_sna0to1(iCv,0,isb)+
     &                     b2srdt_sna0to1(iCv,1,isb)*nb(iCv)
          b2srst_sna_bal(iCv,isb)=b2srst_sna0to1(iCv,0,isb)+
     &                     b2srst_sna0to1(iCv,1,isb)*nb(iCv)
      enddo
      return
      end subroutine update_balance_sources_particles

      subroutine update_balance_sources_momentum(nCv,isb,ub,rob) !djm Jan2017
      use b2mod_constants
     & , only : mp
      use b2mod_b2cmpa
     & , only : am
      implicit none
*   ..input arguments (unchanged on exit)
      integer nCv,isb
      real (kind=R8):: ub(nCv), rob(nCv)
*   ..local variables
      integer iCv,istrai,nstrad
*   ..procedures
      intrinsic size

      nstrad = size(eirene_mc_mapl_smo_bal,3)

      do iCv = 1, nCv
          b2stbr_phys_smo_bal(iCv,isb)=
     &            b2stbr_phys_smo0to3(iCv,0,isb)+
     &            b2stbr_phys_smo0to3(iCv,1,isb)*ub(iCv)+
     &            b2stbr_phys_smo0to3(iCv,2,isb)*rob(iCv)+
     &            b2stbr_phys_smo0to3(iCv,3,isb)*rob(iCv)*ub(iCv)
          b2stbr_bas_smo_bal(iCv,isb)=
     &            b2stbr_bas_smo0to3(iCv,0,isb)+
     &            b2stbr_bas_smo0to3(iCv,1,isb)*ub(iCv)+
     &            b2stbr_bas_smo0to3(iCv,2,isb)*rob(iCv)+
     &            b2stbr_bas_smo0to3(iCv,3,isb)*rob(iCv)*ub(iCv)
          b2stbc_smo_bal(iCv,isb)=b2stbc_smo0to3(iCv,0,isb)+
     &            b2stbc_smo0to3(iCv,1,isb)*ub(iCv)+
     &            b2stbc_smo0to3(iCv,2,isb)*rob(iCv)+
     &            b2stbc_smo0to3(iCv,3,isb)*rob(iCv)*ub(iCv)
          if (allocated(eirene_mc_mapl_smo_bal)) then
            do istrai=1,nstrad
              eirene_mc_mapl_smo_bal(iCv,isb,istrai)=
     &         eirene_mc_mapl_smo0to3(iCv,0,isb,istrai)+
     &         eirene_mc_mapl_smo0to3(iCv,1,isb,istrai)*ub(iCv)+
     &         eirene_mc_mapl_smo0to3(iCv,2,isb,istrai)*rob(iCv)+
     &         eirene_mc_mapl_smo0to3(iCv,3,isb,istrai)*rob(iCv)*
     &                                                    ub(iCv)
              eirene_mc_mmpl_smo_bal(iCv,isb,istrai)=
     &         eirene_mc_mmpl_smo0to3(iCv,0,isb,istrai)+
     &         eirene_mc_mmpl_smo0to3(iCv,1,isb,istrai)*ub(iCv)+
     &         eirene_mc_mmpl_smo0to3(iCv,2,isb,istrai)*rob(iCv)+
     &         eirene_mc_mmpl_smo0to3(iCv,3,isb,istrai)*rob(iCv)*
     &                                                    ub(iCv)
              eirene_mc_mipl_smo_bal(iCv,isb,istrai)=
     &         eirene_mc_mipl_smo0to3(iCv,0,isb,istrai)+
     &         eirene_mc_mipl_smo0to3(iCv,1,isb,istrai)*ub(iCv)+
     &         eirene_mc_mipl_smo0to3(iCv,2,isb,istrai)*rob(iCv)+
     &         eirene_mc_mipl_smo0to3(iCv,3,isb,istrai)*rob(iCv)*
     &                                                    ub(iCv)
              eirene_mc_cppv_smo_bal(iCv,isb,istrai)=
     &         eirene_mc_cppv_smo0to3(iCv,0,isb,istrai)+
     &         eirene_mc_cppv_smo0to3(iCv,1,isb,istrai)*ub(iCv)+
     &         eirene_mc_cppv_smo0to3(iCv,2,isb,istrai)*rob(iCv)+
     &         eirene_mc_cppv_smo0to3(iCv,3,isb,istrai)*rob(iCv)*
     &                                                    ub(iCv)
            enddo
          endif
          b2stbm_smo_bal(iCv,isb)=b2stbm_smo0to3(iCv,0,isb)+
     &            b2stbm_smo0to3(iCv,1,isb)*ub(iCv)+
     &            b2stbm_smo0to3(iCv,2,isb)*rob(iCv)+
     &            b2stbm_smo0to3(iCv,3,isb)*rob(iCv)*ub(iCv)
          ext_smo_bal(iCv,isb)=ext_smo0to3(iCv,0,isb)+
     &            ext_smo0to3(iCv,1,isb)*ub(iCv)+
     &            ext_smo0to3(iCv,2,isb)*rob(iCv)+
     &            ext_smo0to3(iCv,3,isb)*rob(iCv)*ub(iCv)
          b2stel_smq_ion_bal(iCv,isb)=b2stel_smq_ion0to3(iCv,0,isb)+
     &            b2stel_smq_ion0to3(iCv,1,isb)*ub(iCv)+
     &            b2stel_smq_ion0to3(iCv,2,isb)*rob(iCv)+
     &            b2stel_smq_ion0to3(iCv,3,isb)*rob(iCv)*ub(iCv)
          b2stel_smq_rec_bal(iCv,isb)=b2stel_smq_rec0to3(iCv,0,isb)+
     &            b2stel_smq_rec0to3(iCv,1,isb)*ub(iCv)+
     &            b2stel_smq_rec0to3(iCv,2,isb)*rob(iCv)+
     &            b2stel_smq_rec0to3(iCv,3,isb)*rob(iCv)*ub(iCv)
          b2stcx_smq_bal(iCv,isb)=b2stcx_smq0to3(iCv,0,isb)+
     &            b2stcx_smq0to3(iCv,1,isb)*ub(iCv)+
     &            b2stcx_smq0to3(iCv,2,isb)*rob(iCv)+
     &            b2stcx_smq0to3(iCv,3,isb)*rob(iCv)*ub(iCv)
          b2srsm_smo_bal(iCv,isb)=b2srsm_smo0to3(iCv,0,isb)+
     &            b2srsm_smo0to3(iCv,1,isb)*ub(iCv)+
     &            b2srsm_smo0to3(iCv,2,isb)*rob(iCv)+
     &            b2srsm_smo0to3(iCv,3,isb)*rob(iCv)*ub(iCv)
          b2srdt_smo_bal(iCv,isb)=b2srdt_smo0to3(iCv,0,isb)+
     &            b2srdt_smo0to3(iCv,1,isb)*ub(iCv)+
     &            b2srdt_smo0to3(iCv,2,isb)*rob(iCv)+
     &            b2srdt_smo0to3(iCv,3,isb)*rob(iCv)*ub(iCv)
          b2srst_smo_bal(iCv,isb)=b2srst_smo0to3(iCv,0,isb)+
     &            b2srst_smo0to3(iCv,1,isb)*ub(iCv)+
     &            b2srst_smo0to3(iCv,2,isb)*rob(iCv)+
     &            b2srst_smo0to3(iCv,3,isb)*rob(iCv)*ub(iCv)
          b2sifr_smoch_bal(iCv,isb)=b2sifr_smoch0to3(iCv,0,isb)+
     &            b2sifr_smoch0to3(iCv,1,isb)*ub(iCv)+
     &            b2sifr_smoch0to3(iCv,2,isb)*rob(iCv)+
     &            b2sifr_smoch0to3(iCv,3,isb)*rob(iCv)*ub(iCv)
          b2sifr_smotf_bal(iCv,isb)=b2sifr_smotf0to3(iCv,0,isb)+
     &            b2sifr_smotf0to3(iCv,1,isb)*ub(iCv)+
     &            b2sifr_smotf0to3(iCv,2,isb)*rob(iCv)+
     &            b2sifr_smotf0to3(iCv,3,isb)*rob(iCv)*ub(iCv)
          b2siav_smovh_bal(iCv,isb)=b2siav_smovh0to3(iCv,0,isb)+
     &            b2siav_smovh0to3(iCv,1,isb)*ub(iCv)+
     &            b2siav_smovh0to3(iCv,2,isb)*rob(iCv)+
     &            b2siav_smovh0to3(iCv,3,isb)*rob(iCv)*ub(iCv)
          b2siav_smovv_bal(iCv,isb)=b2siav_smovv0to3(iCv,0,isb)+
     &            b2siav_smovv0to3(iCv,1,isb)*ub(iCv)+
     &            b2siav_smovv0to3(iCv,2,isb)*rob(iCv)+
     &            b2siav_smovv0to3(iCv,3,isb)*rob(iCv)*ub(iCv)
          b2sicf_smo_bal(iCv,isb)=b2sicf_smo0to3(iCv,0,isb)+
     &            b2sicf_smo0to3(iCv,1,isb)*ub(iCv)+
     &            b2sicf_smo0to3(iCv,2,isb)*rob(iCv)+
     &            b2sicf_smo0to3(iCv,3,isb)*rob(iCv)*ub(iCv)
          b2sian_smo_bal(iCv,isb)=b2sian_smo0to3(iCv,0,isb)+
     &            b2sian_smo0to3(iCv,1,isb)*ub(iCv)+
     &            b2sian_smo0to3(iCv,2,isb)*rob(iCv)+
     &            b2sian_smo0to3(iCv,3,isb)*rob(iCv)*ub(iCv)      !srv 25.09.17
          b2nxdv_smo_bal(iCv,isb)=b2nxdv_smo0to3(iCv,0,isb)+
     &            b2nxdv_smo0to3(iCv,1,isb)*ub(iCv)+
     &            b2nxdv_smo0to3(iCv,2,isb)*rob(iCv)+
     &            b2nxdv_smo0to3(iCv,3,isb)*rob(iCv)*ub(iCv)
          b2sigp_smogp_bal(iCv,isb)=b2sigp_smogp0to3(iCv,0,isb)+
     &            b2sigp_smogp0to3(iCv,1,isb)*ub(iCv)+
     &            b2sigp_smogp0to3(iCv,2,isb)*rob(iCv)+
     &            b2sigp_smogp0to3(iCv,3,isb)*rob(iCv)*ub(iCv)
          b2sigp_pstat_bal(iCv,isb)=b2sigp_pstat0to3(iCv,0,isb)+
     &            b2sigp_pstat0to3(iCv,1,isb)*ub(iCv)+
     &            b2sigp_pstat0to3(iCv,2,isb)*rob(iCv)+
     &            b2sigp_pstat0to3(iCv,3,isb)*rob(iCv)*ub(iCv)
      enddo
      return
      end subroutine update_balance_sources_momentum

      subroutine update_balance_sources_heat(nCv,ne,ni,te,ti)
      implicit none
*   ..input arguments (unchanged on exit)
      integer nCv
      real (kind=R8):: ne(nCv), ni(nCv),
     &                 te(nCv), ti(nCv)
*   ..local variables
      integer iCv,istrai,nstrad
*   ..procedures
      intrinsic size

      nstrad = size(eirene_mc_eael_she_bal,2)

cdjm Update electron heat sources
      do iCv = 1, nCv
          b2stbr_phys_she_bal(iCv)=
     &     b2stbr_phys_she0to3(iCv,0)+
     &     b2stbr_phys_she0to3(iCv,1)*te(iCv)+
     &     b2stbr_phys_she0to3(iCv,2)*ne(iCv)+
     &     b2stbr_phys_she0to3(iCv,3)*ne(iCv)*te(iCv)
          b2stbr_bas_she_bal(iCv)=
     &     b2stbr_bas_she0to3(iCv,0)+
     &     b2stbr_bas_she0to3(iCv,1)*te(iCv)+
     &     b2stbr_bas_she0to3(iCv,2)*ne(iCv)+
     &     b2stbr_bas_she0to3(iCv,3)*ne(iCv)*te(iCv)
          b2stbr_first_flight_she_bal(iCv)=
     &     b2stbr_first_flight_she0to3(iCv,0)+
     &     b2stbr_first_flight_she0to3(iCv,1)*te(iCv)+
     &     b2stbr_first_flight_she0to3(iCv,2)*ne(iCv)+
     &     b2stbr_first_flight_she0to3(iCv,3)*ne(iCv)*te(iCv)
          b2stbc_she_bal(iCv)=
     &     b2stbc_she0to3(iCv,0)+
     &     b2stbc_she0to3(iCv,1)*te(iCv)+
     &     b2stbc_she0to3(iCv,2)*ne(iCv)+
     &     b2stbc_she0to3(iCv,3)*ne(iCv)*te(iCv)
          if (allocated(eirene_mc_eael_she_bal)) then
            do istrai=1,nstrad
              eirene_mc_eael_she_bal(iCv,istrai)=
     &         eirene_mc_eael_she0to3(iCv,0,istrai)+
     &         eirene_mc_eael_she0to3(iCv,1,istrai)*te(iCv)+
     &         eirene_mc_eael_she0to3(iCv,2,istrai)*ne(iCv)+
     &         eirene_mc_eael_she0to3(iCv,3,istrai)*ne(iCv)*
     &                                                te(iCv)
              eirene_mc_emel_she_bal(iCv,istrai)=
     &         eirene_mc_emel_she0to3(iCv,0,istrai)+
     &         eirene_mc_emel_she0to3(iCv,1,istrai)*te(iCv)+
     &         eirene_mc_emel_she0to3(iCv,2,istrai)*ne(iCv)+
     &         eirene_mc_emel_she0to3(iCv,3,istrai)*ne(iCv)*
     &                                                te(iCv)
              eirene_mc_eiel_she_bal(iCv,istrai)=
     &         eirene_mc_eiel_she0to3(iCv,0,istrai)+
     &         eirene_mc_eiel_she0to3(iCv,1,istrai)*te(iCv)+
     &         eirene_mc_eiel_she0to3(iCv,2,istrai)*ne(iCv)+
     &         eirene_mc_eiel_she0to3(iCv,3,istrai)*ne(iCv)*
     &                                                te(iCv)
              eirene_mc_epel_she_bal(iCv,istrai)=
     &         eirene_mc_epel_she0to3(iCv,0,istrai)+
     &         eirene_mc_epel_she0to3(iCv,1,istrai)*te(iCv)+
     &         eirene_mc_epel_she0to3(iCv,2,istrai)*ne(iCv)+
     &         eirene_mc_epel_she0to3(iCv,3,istrai)*ne(iCv)*
     &                                                te(iCv)
            enddo
          endif
          b2stbm_she_bal(iCv)=
     &     b2stbm_she0to3(iCv,0)+
     &     b2stbm_she0to3(iCv,1)*te(iCv)+
     &     b2stbm_she0to3(iCv,2)*ne(iCv)+
     &     b2stbm_she0to3(iCv,3)*ne(iCv)*te(iCv)
          ext_she_bal(iCv)=
     &     ext_she0to3(iCv,0)+
     &     ext_she0to3(iCv,1)*te(iCv)+
     &     ext_she0to3(iCv,2)*ne(iCv)+
     &     ext_she0to3(iCv,3)*ne(iCv)*te(iCv)
          b2stel_she_bal(iCv)=
     &     b2stel_she0to3(iCv,0)+
     &     b2stel_she0to3(iCv,1)*te(iCv)+
     &     b2stel_she0to3(iCv,2)*ne(iCv)+
     &     b2stel_she0to3(iCv,3)*ne(iCv)*te(iCv)
          b2srsm_she_bal(iCv)=
     &     b2srsm_she0to3(iCv,0)+
     &     b2srsm_she0to3(iCv,1)*te(iCv)+
     &     b2srsm_she0to3(iCv,2)*ne(iCv)+
     &     b2srsm_she0to3(iCv,3)*ne(iCv)*te(iCv)
          b2srdt_she_bal(iCv)=
     &     b2srdt_she0to3(iCv,0)+
     &     b2srdt_she0to3(iCv,1)*te(iCv)+
     &     b2srdt_she0to3(iCv,2)*ne(iCv)+
     &     b2srdt_she0to3(iCv,3)*ne(iCv)*te(iCv)
          b2srst_she_bal(iCv)=
     &     b2srst_she0to3(iCv,0)+
     &     b2srst_she0to3(iCv,1)*te(iCv)+
     &     b2srst_she0to3(iCv,2)*ne(iCv)+
     &     b2srst_she0to3(iCv,3)*ne(iCv)*te(iCv)
          b2sihs_diae_bal(iCv)=
     &     b2sihs_diae0to3(iCv,0)+
     &     b2sihs_diae0to3(iCv,1)*te(iCv)+
     &     b2sihs_diae0to3(iCv,2)*ne(iCv)+
     &     b2sihs_diae0to3(iCv,3)*ne(iCv)*te(iCv)
          b2sihs_divue_bal(iCv)=
     &     b2sihs_divue0to3(iCv,0)+
     &     b2sihs_divue0to3(iCv,1)*te(iCv)+
     &     b2sihs_divue0to3(iCv,2)*ne(iCv)+
     &     b2sihs_divue0to3(iCv,3)*ne(iCv)*te(iCv)
          b2sihs_exbe_bal(iCv)=
     &     b2sihs_exbe0to3(iCv,0)+
     &     b2sihs_exbe0to3(iCv,1)*te(iCv)+
     &     b2sihs_exbe0to3(iCv,2)*ne(iCv)+
     &     b2sihs_exbe0to3(iCv,3)*ne(iCv)*te(iCv)
          b2sihs_joule_bal(iCv)=
     &     b2sihs_joule0to3(iCv,0)+
     &     b2sihs_joule0to3(iCv,1)*te(iCv)+
     &     b2sihs_joule0to3(iCv,2)*ne(iCv)+
     &     b2sihs_joule0to3(iCv,3)*ne(iCv)*te(iCv)
      enddo
cdjm Update ion heat sources
      do iCv = 1, nCv
          b2stbr_phys_shi_bal(iCv)=
     &     b2stbr_phys_shi0to3(iCv,0)+
     &     b2stbr_phys_shi0to3(iCv,1)*ti(iCv)+
     &     b2stbr_phys_shi0to3(iCv,2)*ni(iCv)+
     &     b2stbr_phys_shi0to3(iCv,3)*ni(iCv)*ti(iCv)
          b2stbr_bas_shi_bal(iCv)=
     &     b2stbr_bas_shi0to3(iCv,0)+
     &     b2stbr_bas_shi0to3(iCv,1)*ti(iCv)+
     &     b2stbr_bas_shi0to3(iCv,2)*ni(iCv)+
     &     b2stbr_bas_shi0to3(iCv,3)*ni(iCv)*ti(iCv)
          b2stbr_first_flight_shi_bal(iCv)=
     &     b2stbr_first_flight_shi0to3(iCv,0)+
     &     b2stbr_first_flight_shi0to3(iCv,1)*ti(iCv)+
     &     b2stbr_first_flight_shi0to3(iCv,2)*ni(iCv)+
     &     b2stbr_first_flight_shi0to3(iCv,3)*ni(iCv)*ti(iCv)
          b2stbc_shi_bal(iCv)=
     &     b2stbc_shi0to3(iCv,0)+
     &     b2stbc_shi0to3(iCv,1)*ti(iCv)+
     &     b2stbc_shi0to3(iCv,2)*ni(iCv)+
     &     b2stbc_shi0to3(iCv,3)*ni(iCv)*ti(iCv)
          if (allocated(eirene_mc_eapl_shi_bal)) then
            do istrai=1,nstrad
              eirene_mc_eapl_shi_bal(iCv,istrai)=
     &         eirene_mc_eapl_shi0to3(iCv,0,istrai)+
     &         eirene_mc_eapl_shi0to3(iCv,1,istrai)*ti(iCv)+
     &         eirene_mc_eapl_shi0to3(iCv,2,istrai)*ni(iCv)+
     &         eirene_mc_eapl_shi0to3(iCv,3,istrai)*ni(iCv)*
     &                                                ti(iCv)
              eirene_mc_empl_shi_bal(iCv,istrai)=
     &         eirene_mc_empl_shi0to3(iCv,0,istrai)+
     &         eirene_mc_empl_shi0to3(iCv,1,istrai)*ti(iCv)+
     &         eirene_mc_empl_shi0to3(iCv,2,istrai)*ni(iCv)+
     &         eirene_mc_empl_shi0to3(iCv,3,istrai)*ni(iCv)*
     &                                                ti(iCv)
              eirene_mc_eipl_shi_bal(iCv,istrai)=
     &         eirene_mc_eipl_shi0to3(iCv,0,istrai)+
     &         eirene_mc_eipl_shi0to3(iCv,1,istrai)*ti(iCv)+
     &         eirene_mc_eipl_shi0to3(iCv,2,istrai)*ni(iCv)+
     &         eirene_mc_eipl_shi0to3(iCv,3,istrai)*ni(iCv)*
     &                                                ti(iCv)
              eirene_mc_eppl_shi_bal(iCv,istrai)=
     &         eirene_mc_eppl_shi0to3(iCv,0,istrai)+
     &         eirene_mc_eppl_shi0to3(iCv,1,istrai)*ti(iCv)+
     &         eirene_mc_eppl_shi0to3(iCv,2,istrai)*ni(iCv)+
     &         eirene_mc_eppl_shi0to3(iCv,3,istrai)*ni(iCv)*
     &                                                ti(iCv)
            enddo
          endif
          b2stbm_shi_bal(iCv)=
     &     b2stbm_shi0to3(iCv,0)+
     &     b2stbm_shi0to3(iCv,1)*ti(iCv)+
     &     b2stbm_shi0to3(iCv,2)*ni(iCv)+
     &     b2stbm_shi0to3(iCv,3)*ni(iCv)*ti(iCv)
          ext_shi_bal(iCv)=
     &     ext_shi0to3(iCv,0)+
     &     ext_shi0to3(iCv,1)*ti(iCv)+
     &     ext_shi0to3(iCv,2)*ni(iCv)+
     &     ext_shi0to3(iCv,3)*ni(iCv)*ti(iCv)
          b2stel_shi_ion_bal(iCv)=
     &     b2stel_shi_ion0to3(iCv,0)+
     &     b2stel_shi_ion0to3(iCv,1)*ti(iCv)+
     &     b2stel_shi_ion0to3(iCv,2)*ni(iCv)+
     &     b2stel_shi_ion0to3(iCv,3)*ni(iCv)*ti(iCv)
          b2stel_shi_rec_bal(iCv)=
     &     b2stel_shi_rec0to3(iCv,0)+
     &     b2stel_shi_rec0to3(iCv,1)*ti(iCv)+
     &     b2stel_shi_rec0to3(iCv,2)*ni(iCv)+
     &     b2stel_shi_rec0to3(iCv,3)*ni(iCv)*ti(iCv)
          b2stcx_shi_bal(iCv)=
     &     b2stcx_shi0to3(iCv,0)+
     &     b2stcx_shi0to3(iCv,1)*ti(iCv)+
     &     b2stcx_shi0to3(iCv,2)*ni(iCv)+
     &     b2stcx_shi0to3(iCv,3)*ni(iCv)*ti(iCv)
          b2srsm_shi_bal(iCv)=
     &     b2srsm_shi0to3(iCv,0)+
     &     b2srsm_shi0to3(iCv,1)*ti(iCv)+
     &     b2srsm_shi0to3(iCv,2)*ni(iCv)+
     &     b2srsm_shi0to3(iCv,3)*ni(iCv)*ti(iCv)
          b2srdt_shi_bal(iCv)=
     &     b2srdt_shi0to3(iCv,0)+
     &     b2srdt_shi0to3(iCv,1)*ti(iCv)+
     &     b2srdt_shi0to3(iCv,2)*ni(iCv)+
     &     b2srdt_shi0to3(iCv,3)*ni(iCv)*ti(iCv)
          b2srst_shi_bal(iCv)=
     &     b2srst_shi0to3(iCv,0)+
     &     b2srst_shi0to3(iCv,1)*ti(iCv)+
     &     b2srst_shi0to3(iCv,2)*ni(iCv)+
     &     b2srst_shi0to3(iCv,3)*ni(iCv)*ti(iCv)
          b2sihs_diaa_bal(iCv)=
     &     b2sihs_diaa0to3(iCv,0)+
     &     b2sihs_diaa0to3(iCv,1)*ti(iCv)+
     &     b2sihs_diaa0to3(iCv,2)*ni(iCv)+
     &     b2sihs_diaa0to3(iCv,3)*ni(iCv)*ti(iCv)
          b2sihs_divua_bal(iCv)=
     &     b2sihs_divua0to3(iCv,0)+
     &     b2sihs_divua0to3(iCv,1)*ti(iCv)+
     &     b2sihs_divua0to3(iCv,2)*ni(iCv)+
     &     b2sihs_divua0to3(iCv,3)*ni(iCv)*ti(iCv)
          b2sihs_exba_bal(iCv)=
     &     b2sihs_exba0to3(iCv,0)+
     &     b2sihs_exba0to3(iCv,1)*ti(iCv)+
     &     b2sihs_exba0to3(iCv,2)*ni(iCv)+
     &     b2sihs_exba0to3(iCv,3)*ni(iCv)*ti(iCv)
          b2sihs_visa_bal(iCv)=
     &     b2sihs_visa0to3(iCv,0)+
     &     b2sihs_visa0to3(iCv,1)*ti(iCv)+
     &     b2sihs_visa0to3(iCv,2)*ni(iCv)+
     &     b2sihs_visa0to3(iCv,3)*ni(iCv)*ti(iCv)
          b2sihs_fraa_bal(iCv)=
     &     b2sihs_fraa0to3(iCv,0)+
     &     b2sihs_fraa0to3(iCv,1)*ti(iCv)+
     &     b2sihs_fraa0to3(iCv,2)*ni(iCv)+
     &     b2sihs_fraa0to3(iCv,3)*ni(iCv)*ti(iCv)
      enddo
      return
      end subroutine update_balance_sources_heat

      subroutine update_average_balance(state,itim,ntim)
      use b2us_plasma
#ifdef B25_EIRENE
#ifndef SOLPS4_3
      use eirmod_extrab25
     & , only : dab2,dmb2,tab2,tmb2,rfluxa,rfluxm,refluxa,refluxm,
     &          pfluxa,pfluxm,pefluxa,pefluxm
#else
      use eirdiag ! IGNORE
     & , only : dab2,dmb2,tab2,tmb2,rfluxa,rfluxm,refluxa,refluxm,
     &          pfluxa,pfluxm,pefluxa,pefluxm
#endif
#else
      use b2mod_eirdiag
#endif
      implicit none
*   ..input arguments (unchanged on exit)
      integer, intent(in) :: itim,ntim
      type(B2State), intent(in) :: state
*   ..local variables
      real (kind=R8) :: avmlt, curmlt
*   ..procedures
      real (kind=R8) :: rratio
      external rratio

      if (itim.eq.0) then
*     ..Arrays for particle balance:
        fna_pinch_av=0.0_R8
        fna_pll_av=0.0_R8
        fna_drift_av=0.0_R8
        fna_ch_av=0.0_R8
        fna_nanom_av=0.0_R8
        fna_panom_av=0.0_R8
        fna_pschused_av=0.0_R8
        b2stbr_phys_sna_bal_av=0.0_R8
        b2stbr_bas_sna_bal_av=0.0_R8
        b2stbr_first_flight_sna_bal_av=0.0_R8
        b2stbc_sna_bal_av=0.0_R8
        b2stbm_sna_bal_av=0.0_R8
        ext_sna_bal_av=0.0_R8
        b2stel_sna_ion_bal_av=0.0_R8
        b2stel_sna_rec_bal_av=0.0_R8
        b2stcx_sna_bal_av=0.0_R8
        b2srsm_sna_bal_av=0.0_R8
        b2srdt_sna_bal_av=0.0_R8
        b2srst_sna_bal_av=0.0_R8
        resco_av=0.0_R8
        if (allocated(eirene_mc_papl_sna_bal)) then
          eirene_mc_papl_sna_bal_av=0.0_R8
          eirene_mc_pmpl_sna_bal_av=0.0_R8
          eirene_mc_pipl_sna_bal_av=0.0_R8
          eirene_mc_pppl_sna_bal_av=0.0_R8
          eirene_mc_core_sna_bal_av=0.0_R8
        endif
*     ..Arrays for momentum balance:
        fmo_flua_av=0.0_R8
        fmo_cvsa_av=0.0_R8
        fmo_b2nxfv_av=0.0_R8
        b2stbr_phys_smo_bal_av=0.0_R8
        b2stbr_bas_smo_bal_av=0.0_R8
        b2stbc_smo_bal_av=0.0_R8
        if (allocated(eirene_mc_mapl_smo_bal)) then
          eirene_mc_mapl_smo_bal_av=0.0_R8
          eirene_mc_mmpl_smo_bal_av=0.0_R8
          eirene_mc_mipl_smo_bal_av=0.0_R8
          eirene_mc_cppv_smo_bal_av=0.0_R8
        endif
        b2stbm_smo_bal_av=0.0_R8
        ext_smo_bal_av=0.0_R8
        b2stel_smq_ion_bal_av=0.0_R8
        b2stel_smq_rec_bal_av=0.0_R8
        b2stcx_smq_bal_av=0.0_R8
        b2srsm_smo_bal_av=0.0_R8
        b2srdt_smo_bal_av=0.0_R8
        b2srst_smo_bal_av=0.0_R8
        b2sifr_smoch_bal_av=0.0_R8
        b2sifr_smotf_bal_av=0.0_R8
        b2siav_smovh_bal_av=0.0_R8
        b2siav_smovv_bal_av=0.0_R8
        b2sicf_smo_bal_av=0.0_R8
        b2sian_smo_bal_av=0.0_R8           !srv 25.09.17
        b2nxdv_smo_bal_av=0.0_R8
        b2sigp_smogp_bal_av=0.0_R8
        b2sigp_pstat_bal_av=0.0_R8
        resmo_av=0.0_R8
*     ..Arrays for electron heat balance:
        fhe_32_av=0.0_R8
        fhe_52_av=0.0_R8
        fhe_thermj_av=0.0_R8
        fhe_cond_av=0.0_R8
        fhe_dia_av=0.0_R8
        fhe_ecrb_av=0.0_R8
        fhe_strange_av=0.0_R8
        fhe_pschused_av=0.0_R8
        b2stbr_phys_she_bal_av=0.0_R8
        b2stbr_bas_she_bal_av=0.0_R8
        b2stbr_first_flight_she_bal_av=0.0_R8
        b2stbc_she_bal_av=0.0_R8
        if (allocated(eirene_mc_eael_she_bal)) then
          eirene_mc_eael_she_bal_av=0.0_R8
          eirene_mc_emel_she_bal_av=0.0_R8
          eirene_mc_eiel_she_bal_av=0.0_R8
          eirene_mc_epel_she_bal_av=0.0_R8
        endif
        b2stbm_she_bal_av=0.0_R8
        ext_she_bal_av=0.0_R8
        b2stel_she_bal_av=0.0_R8
        b2srsm_she_bal_av=0.0_R8
        b2srdt_she_bal_av=0.0_R8
        b2srst_she_bal_av=0.0_R8
        b2sihs_diae_bal_av=0.0_R8
        b2sihs_divue_bal_av=0.0_R8
        b2sihs_exbe_bal_av=0.0_R8
        b2sihs_joule_bal_av=0.0_R8
        b2npht_shei_bal_av=0.0_R8
        reshe_av=0.0_R8
*     ..Arrays for ion heat balance:
        fhi_32_av=0.0_R8
        fhi_52_av=0.0_R8
        fhi_cond_av=0.0_R8
        fhi_dia_av=0.0_R8
        fhi_ecrb_av=0.0_R8
        fhi_strange_av=0.0_R8
        fhi_pschused_av=0.0_R8
        fhi_inert_av=0.0_R8
        fhi_vispar_av=0.0_R8
        fhi_visper_av=0.0_R8
        fhi_visq_av=0.0_R8
        fhi_anml_av=0.0_R8
        fhi_kevis_av=0.0_R8
        b2stbr_phys_shi_bal_av=0.0_R8
        b2stbr_bas_shi_bal_av=0.0_R8
        b2stbr_first_flight_shi_bal_av=0.0_R8
        b2stbc_shi_bal_av=0.0_R8
        if (allocated(eirene_mc_eapl_shi_bal)) then
          eirene_mc_eapl_shi_bal_av=0.0_R8
          eirene_mc_empl_shi_bal_av=0.0_R8
          eirene_mc_eipl_shi_bal_av=0.0_R8
          eirene_mc_eppl_shi_bal_av=0.0_R8
        endif
        b2stbm_shi_bal_av=0.0_R8
        ext_shi_bal_av=0.0_R8
        b2stel_shi_ion_bal_av=0.0_R8
        b2stel_shi_rec_bal_av=0.0_R8
        b2stcx_shi_bal_av=0.0_R8
        b2srsm_shi_bal_av=0.0_R8
        b2srdt_shi_bal_av=0.0_R8
        b2srst_shi_bal_av=0.0_R8
        b2sihs_diaa_bal_av=0.0_R8
        b2sihs_divua_bal_av=0.0_R8
        b2sihs_exba_bal_av=0.0_R8
        b2sihs_visa_bal_av=0.0_R8
        b2sihs_fraa_bal_av=0.0_R8
        reshi_av=0.0_R8
*     ..Additional averaged quantities
        ne_av=0.0_R8
        na_av=0.0_R8
        ua_av=0.0_R8
        po_av=0.0_R8
        te_av=0.0_R8
        ti_av=0.0_R8
        rpt_av=0.0_R8
        kinrgy_av=0.0_R8
        fne_av=0.0_R8
        if (allocated(dab2_av)) then
          dab2_av=0.0_R8
          dmb2_av=0.0_R8
          tab2_av=0.0_R8
          tmb2_av=0.0_R8
          rfluxa_av=0.0_R8
          rfluxm_av=0.0_R8
          refluxa_av=0.0_R8
          refluxm_av=0.0_R8
          pfluxa_av=0.0_R8
          pfluxm_av=0.0_R8
          pefluxa_av=0.0_R8
          pefluxm_av=0.0_R8
        endif
      endif

      if (balance_average.ne.0) then
        avmlt = 1.0_R8
        curmlt = rratio(1,ntim)
      else
        avmlt = 0.0_R8
        curmlt = 1.0_R8
      endif

*   ..Arrays for particle balance:
      fna_pinch_av = avmlt*fna_pinch_av +
     &               curmlt*fna_pinch
      fna_pll_av = avmlt*fna_pll_av +
     &             curmlt*fna_pll
      fna_drift_av = avmlt*fna_drift_av +
     &               curmlt*fna_drift
      fna_ch_av = avmlt*fna_ch_av +
     &            curmlt*fna_ch
      fna_nanom_av = avmlt*fna_nanom_av +
     &               curmlt*fna_nanom
      fna_panom_av = avmlt*fna_panom_av +
     &               curmlt*fna_panom
      fna_pschused_av = avmlt*fna_pschused_av +
     &                  curmlt*fna_pschused
      b2stbr_phys_sna_bal_av = avmlt*b2stbr_phys_sna_bal_av +
     &                         curmlt*b2stbr_phys_sna_bal
      b2stbr_bas_sna_bal_av = avmlt*b2stbr_bas_sna_bal_av +
     &                        curmlt*b2stbr_bas_sna_bal
      b2stbr_first_flight_sna_bal_av =
     &                  avmlt*b2stbr_first_flight_sna_bal_av +
     &                  curmlt*b2stbr_first_flight_sna_bal
      b2stbc_sna_bal_av = avmlt*b2stbc_sna_bal_av +
     &                    curmlt*b2stbc_sna_bal
      if (allocated(eirene_mc_papl_sna_bal)) then
        eirene_mc_papl_sna_bal_av = avmlt*eirene_mc_papl_sna_bal_av +
     &                              curmlt*eirene_mc_papl_sna_bal
        eirene_mc_pmpl_sna_bal_av = avmlt*eirene_mc_pmpl_sna_bal_av +
     &                              curmlt*eirene_mc_pmpl_sna_bal
        eirene_mc_pipl_sna_bal_av = avmlt*eirene_mc_pipl_sna_bal_av +
     &                              curmlt*eirene_mc_pipl_sna_bal
        eirene_mc_pppl_sna_bal_av = avmlt*eirene_mc_pppl_sna_bal_av +
     &                              curmlt*eirene_mc_pppl_sna_bal
        eirene_mc_core_sna_bal_av = avmlt*eirene_mc_core_sna_bal_av +
     &                              curmlt*eirene_mc_core_sna_bal
      endif
      b2stbm_sna_bal_av = avmlt*b2stbm_sna_bal_av +
     &                    curmlt*b2stbm_sna_bal
      ext_sna_bal_av = avmlt*ext_sna_bal_av +
     &                 curmlt*ext_sna_bal
      b2stel_sna_ion_bal_av = avmlt*b2stel_sna_ion_bal_av +
     &                        curmlt*b2stel_sna_ion_bal
      b2stel_sna_rec_bal_av = avmlt*b2stel_sna_rec_bal_av +
     &                        curmlt*b2stel_sna_rec_bal
      b2stcx_sna_bal_av = avmlt*b2stcx_sna_bal_av +
     &                    curmlt*b2stcx_sna_bal
      b2srsm_sna_bal_av = avmlt*b2srsm_sna_bal_av +
     &                    curmlt*b2srsm_sna_bal
      b2srdt_sna_bal_av = avmlt*b2srdt_sna_bal_av +
     &                    curmlt*b2srdt_sna_bal
      b2srst_sna_bal_av = avmlt*b2srst_sna_bal_av +
     &                    curmlt*b2srst_sna_bal
      resco_av = avmlt*resco_av +
     &           curmlt*state%dv%resco
*   ..Arrays for momentum balance:
      fmo_flua_av = avmlt*fmo_flua_av +
     &              curmlt*fmo_flua
      fmo_cvsa_av = avmlt*fmo_cvsa_av +
     &              curmlt*fmo_cvsa
      fmo_b2nxfv_av = avmlt*fmo_b2nxfv_av +
     &                curmlt*fmo_b2nxfv
      b2stbr_phys_smo_bal_av = avmlt*b2stbr_phys_smo_bal_av +
     &                         curmlt*b2stbr_phys_smo_bal
      b2stbr_bas_smo_bal_av = avmlt*b2stbr_bas_smo_bal_av +
     &                        curmlt*b2stbr_bas_smo_bal
      b2stbc_smo_bal_av = avmlt*b2stbc_smo_bal_av +
     &                    curmlt*b2stbc_smo_bal
      if (allocated(eirene_mc_mapl_smo_bal)) then
        eirene_mc_mapl_smo_bal_av = avmlt*eirene_mc_mapl_smo_bal_av +
     &                              curmlt*eirene_mc_mapl_smo_bal
        eirene_mc_mmpl_smo_bal_av = avmlt*eirene_mc_mmpl_smo_bal_av +
     &                              curmlt*eirene_mc_mmpl_smo_bal
        eirene_mc_mipl_smo_bal_av = avmlt*eirene_mc_mipl_smo_bal_av +
     &                              curmlt*eirene_mc_mipl_smo_bal
        eirene_mc_cppv_smo_bal_av = avmlt*eirene_mc_cppv_smo_bal_av +
     &                              curmlt*eirene_mc_cppv_smo_bal
      endif
      b2stbm_smo_bal_av = avmlt*b2stbm_smo_bal_av +
     &                    curmlt*b2stbm_smo_bal
      ext_smo_bal_av = avmlt*ext_smo_bal_av +
     &                 curmlt*ext_smo_bal
      b2stel_smq_ion_bal_av = avmlt*b2stel_smq_ion_bal_av +
     &                        curmlt*b2stel_smq_ion_bal
      b2stel_smq_rec_bal_av = avmlt*b2stel_smq_rec_bal_av +
     &                        curmlt*b2stel_smq_rec_bal
      b2stcx_smq_bal_av = avmlt*b2stcx_smq_bal_av +
     &                    curmlt*b2stcx_smq_bal
      b2srsm_smo_bal_av = avmlt*b2srsm_smo_bal_av +
     &                    curmlt*b2srsm_smo_bal
      b2srdt_smo_bal_av = avmlt*b2srdt_smo_bal_av +
     &                    curmlt*b2srdt_smo_bal
      b2srst_smo_bal_av = avmlt*b2srst_smo_bal_av +
     &                    curmlt*b2srst_smo_bal
      b2sifr_smoch_bal_av = avmlt*b2sifr_smoch_bal_av +
     &                      curmlt*b2sifr_smoch_bal
      b2sifr_smotf_bal_av = avmlt*b2sifr_smotf_bal_av +
     &                      curmlt*b2sifr_smotf_bal
      b2siav_smovh_bal_av = avmlt*b2siav_smovh_bal_av +
     &                      curmlt*b2siav_smovh_bal
      b2siav_smovv_bal_av = avmlt*b2siav_smovv_bal_av +
     &                      curmlt*b2siav_smovv_bal
      b2sicf_smo_bal_av = avmlt*b2sicf_smo_bal_av +
     &                    curmlt*b2sicf_smo_bal
      b2sian_smo_bal_av = avmlt*b2sian_smo_bal_av +
     &                    curmlt*b2sian_smo_bal           !srv 25.09.17
      b2nxdv_smo_bal_av = avmlt*b2nxdv_smo_bal_av +
     &                    curmlt*b2nxdv_smo_bal
      b2sigp_smogp_bal_av = avmlt*b2sigp_smogp_bal_av +
     &                      curmlt*b2sigp_smogp_bal
      b2sigp_pstat_bal_av = avmlt*b2sigp_pstat_bal_av +
     &                      curmlt*b2sigp_pstat_bal
      resmo_av = avmlt*resmo_av +
     &           curmlt*state%dv%resmo
*   ..Arrays for electron heat balance:
      fhe_32_av = avmlt*fhe_32_av +
     &            curmlt*fhe_32
      fhe_52_av = avmlt*fhe_52_av +
     &            curmlt*fhe_52
      fhe_thermj_av = avmlt*fhe_thermj_av +
     &                curmlt*fhe_thermj
      fhe_cond_av = avmlt*fhe_cond_av +
     &              curmlt*fhe_cond
      fhe_dia_av = avmlt*fhe_dia_av +
     &             curmlt*fhe_dia
      fhe_ecrb_av = avmlt*fhe_ecrb_av +
     &              curmlt*fhe_ecrb
      fhe_strange_av = avmlt*fhe_strange_av +
     &                 curmlt*fhe_strange
      fhe_pschused_av = avmlt*fhe_pschused_av +
     &                  curmlt*fhe_pschused
      b2stbr_phys_she_bal_av = avmlt*b2stbr_phys_she_bal_av +
     &                         curmlt*b2stbr_phys_she_bal
      b2stbr_bas_she_bal_av = avmlt*b2stbr_bas_she_bal_av +
     &                        curmlt*b2stbr_bas_she_bal
      b2stbr_first_flight_she_bal_av = 
     &            avmlt*b2stbr_first_flight_she_bal_av +
     &            curmlt*b2stbr_first_flight_she_bal
      b2stbc_she_bal_av = avmlt*b2stbc_she_bal_av +
     &                    curmlt*b2stbc_she_bal
      if (allocated(eirene_mc_eael_she_bal)) then
        eirene_mc_eael_she_bal_av = avmlt*eirene_mc_eael_she_bal_av +
     &                              curmlt*eirene_mc_eael_she_bal
        eirene_mc_emel_she_bal_av = avmlt*eirene_mc_emel_she_bal_av +
     &                              curmlt*eirene_mc_emel_she_bal
        eirene_mc_eiel_she_bal_av = avmlt*eirene_mc_eiel_she_bal_av +
     &                              curmlt*eirene_mc_eiel_she_bal
        eirene_mc_epel_she_bal_av = avmlt*eirene_mc_epel_she_bal_av +
     &                              curmlt*eirene_mc_epel_she_bal
      endif
      b2stbm_she_bal_av = avmlt*b2stbm_she_bal_av +
     &                    curmlt*b2stbm_she_bal
      ext_she_bal_av = avmlt*ext_she_bal_av +
     &                 curmlt*ext_she_bal
      b2stel_she_bal_av = avmlt*b2stel_she_bal_av +
     &                    curmlt*b2stel_she_bal
      b2srsm_she_bal_av = avmlt*b2srsm_she_bal_av +
     &                    curmlt*b2srsm_she_bal
      b2srdt_she_bal_av = avmlt*b2srdt_she_bal_av +
     &                    curmlt*b2srdt_she_bal
      b2srst_she_bal_av = avmlt*b2srst_she_bal_av +
     &                    curmlt*b2srst_she_bal
      b2sihs_diae_bal_av = avmlt*b2sihs_diae_bal_av +
     &                     curmlt*b2sihs_diae_bal
      b2sihs_divue_bal_av = avmlt*b2sihs_divue_bal_av +
     &                      curmlt*b2sihs_divue_bal
      b2sihs_exbe_bal_av = avmlt*b2sihs_exbe_bal_av +
     &                     curmlt*b2sihs_exbe_bal
      b2sihs_joule_bal_av = avmlt*b2sihs_joule_bal_av +
     &                      curmlt*b2sihs_joule_bal
      b2npht_shei_bal_av = avmlt*b2npht_shei_bal_av +
     &                     curmlt*b2npht_shei_bal
      reshe_av = avmlt*reshe_av +
     &           curmlt*state%dv%reshe
*   ..Arrays for ion heat balance:
      fhi_32_av = avmlt*fhi_32_av +
     &            curmlt*fhi_32
      fhi_52_av = avmlt*fhi_52_av +
     &            curmlt*fhi_52
      fhi_cond_av = avmlt*fhi_cond_av +
     &              curmlt*fhi_cond
      fhi_dia_av = avmlt*fhi_dia_av +
     &             curmlt*fhi_dia
      fhi_ecrb_av = avmlt*fhi_ecrb_av +
     &              curmlt*fhi_ecrb
      fhi_strange_av = avmlt*fhi_strange_av +
     &                 curmlt*fhi_strange
      fhi_pschused_av = avmlt*fhi_pschused_av +
     &                  curmlt*fhi_pschused
      fhi_inert_av = avmlt*fhi_inert_av +
     &               curmlt*fhi_inert
      fhi_vispar_av = avmlt*fhi_vispar_av +
     &                curmlt*fhi_vispar
      fhi_visper_av = avmlt*fhi_visper_av +
     &                curmlt*fhi_visper
      fhi_visq_av = avmlt*fhi_visq_av +
     &                curmlt*fhi_visq
      fhi_anml_av = avmlt*fhi_anml_av +
     &              curmlt*fhi_anml
      fhi_kevis_av = avmlt*fhi_kevis_av +
     &               curmlt*fhi_kevis
      b2stbr_phys_shi_bal_av = avmlt*b2stbr_phys_shi_bal_av +
     &                         curmlt*b2stbr_phys_shi_bal
      b2stbr_bas_shi_bal_av = avmlt*b2stbr_bas_shi_bal_av +
     &                        curmlt*b2stbr_bas_shi_bal
      b2stbr_first_flight_shi_bal_av = 
     &            avmlt*b2stbr_first_flight_shi_bal_av +
     &            curmlt*b2stbr_first_flight_shi_bal
      b2stbc_shi_bal_av = avmlt*b2stbc_shi_bal_av +
     &                    curmlt*b2stbc_shi_bal
      if (allocated(eirene_mc_eapl_shi_bal)) then
        eirene_mc_eapl_shi_bal_av = avmlt*eirene_mc_eapl_shi_bal_av +
     &                              curmlt*eirene_mc_eapl_shi_bal
        eirene_mc_empl_shi_bal_av = avmlt*eirene_mc_empl_shi_bal_av +
     &                              curmlt*eirene_mc_empl_shi_bal
        eirene_mc_eipl_shi_bal_av = avmlt*eirene_mc_eipl_shi_bal_av +
     &                              curmlt*eirene_mc_eipl_shi_bal
        eirene_mc_eppl_shi_bal_av = avmlt*eirene_mc_eppl_shi_bal_av +
     &                              curmlt*eirene_mc_eppl_shi_bal
      endif
      b2stbm_shi_bal_av = avmlt*b2stbm_shi_bal_av +
     &                    curmlt*b2stbm_shi_bal
      ext_shi_bal_av = avmlt*ext_shi_bal_av +
     &                 curmlt*ext_shi_bal
      b2stel_shi_ion_bal_av = avmlt*b2stel_shi_ion_bal_av +
     &                        curmlt*b2stel_shi_ion_bal
      b2stel_shi_rec_bal_av = avmlt*b2stel_shi_rec_bal_av +
     &                        curmlt*b2stel_shi_rec_bal
      b2stcx_shi_bal_av = avmlt*b2stcx_shi_bal_av +
     &                    curmlt*b2stcx_shi_bal
      b2srsm_shi_bal_av = avmlt*b2srsm_shi_bal_av +
     &                    curmlt*b2srsm_shi_bal
      b2srdt_shi_bal_av = avmlt*b2srdt_shi_bal_av +
     &                    curmlt*b2srdt_shi_bal
      b2srst_shi_bal_av = avmlt*b2srst_shi_bal_av +
     &                    curmlt*b2srst_shi_bal
      b2sihs_diaa_bal_av = avmlt*b2sihs_diaa_bal_av +
     &                     curmlt*b2sihs_diaa_bal
      b2sihs_divua_bal_av = avmlt*b2sihs_divua_bal_av +
     &                      curmlt*b2sihs_divua_bal
      b2sihs_exba_bal_av = avmlt*b2sihs_exba_bal_av +
     &                     curmlt*b2sihs_exba_bal
      b2sihs_visa_bal_av = avmlt*b2sihs_visa_bal_av +
     &                     curmlt*b2sihs_visa_bal
      b2sihs_fraa_bal_av = avmlt*b2sihs_fraa_bal_av +
     &                     curmlt*b2sihs_fraa_bal
      reshi_av = avmlt*reshi_av +
     &           curmlt*state%dv%reshi
*   ..Additional averaged quantities
      ne_av=avmlt*ne_av+curmlt*state%dv%ne
      na_av=avmlt*na_av+curmlt*state%pl%na
      ua_av=avmlt*ua_av+curmlt*state%pl%ua
      po_av=avmlt*po_av+curmlt*state%pl%po
      te_av=avmlt*te_av+curmlt*state%pl%te
      ti_av=avmlt*ti_av+curmlt*state%pl%ti
      rpt_av=avmlt*rpt_av+curmlt*state%rt%rpt
      kinrgy_av=avmlt*kinrgy_av+curmlt*state%dv%kinrgy
      fne_av=avmlt*fne_av+curmlt*state%dv%fne
      if (allocated(dab2_av)) then
#ifndef SOLPS4_3
        dab2_av=avmlt*dab2_av+curmlt*dab2
        dmb2_av=avmlt*dmb2_av+curmlt*dmb2
        tab2_av=avmlt*tab2_av+curmlt*tab2
        tmb2_av=avmlt*tmb2_av+curmlt*tmb2
        rfluxa_av=avmlt*rfluxa_av+curmlt*rfluxa
        rfluxm_av=avmlt*rfluxm_av+curmlt*rfluxm
        refluxa_av=avmlt*refluxa_av+curmlt*refluxa
        refluxm_av=avmlt*refluxm_av+curmlt*refluxm
        pfluxa_av=avmlt*pfluxa_av+curmlt*pfluxa
        pfluxm_av=avmlt*pfluxm_av+curmlt*pfluxm
        pefluxa_av=avmlt*pefluxa_av+curmlt*pefluxa
        pefluxm_av=avmlt*pefluxm_av+curmlt*pefluxm
#else
        dab2_av(:,:,:)=avmlt*dab2_av(:,:,:)+curmlt*dab2(:,:,:,1)
        dmb2_av(:,:,:)=avmlt*dmb2_av(:,:,:)+curmlt*dmb2(:,:,:,1)
        tab2_av(:,:,:)=avmlt*tab2_av(:,:,:)+curmlt*tab2(:,:,:,1)
        tmb2_av(:,:,:)=avmlt*tmb2_av(:,:,:)+curmlt*tmb2(:,:,:,1)
        rfluxa_av(:,:,:)=avmlt*rfluxa_av(:,:,:)+curmlt*rfluxa(:,:,:,1)
        rfluxm_av(:,:,:)=avmlt*rfluxm_av(:,:,:)+curmlt*rfluxm(:,:,:,1)
        refluxa_av(:,:,:)=avmlt*refluxa_av(:,:,:)+
     .                                         curmlt*refluxa(:,:,:,1)
        refluxm_av(:,:,:)=avmlt*refluxm_av(:,:,:)+
     .                                         curmlt*refluxm(:,:,:,1)
        pfluxa_av(:,:,:)=avmlt*pfluxa_av(:,:,:)+curmlt*pfluxa(:,:,:,1)
        pfluxm_av(:,:,:)=avmlt*pfluxm_av(:,:,:)+curmlt*pfluxm(:,:,:,1)
        pefluxa_av(:,:,:)=avmlt*pefluxa_av(:,:,:)+
     .                                         curmlt*pefluxa(:,:,:,1)
        pefluxm_av(:,:,:)=avmlt*pefluxm_av(:,:,:)+
     .                                         curmlt*pefluxm(:,:,:,1)
#endif
      endif

      end subroutine update_average_balance

#ifdef WG_TODO
      subroutine write_balance(nxd,nyd,nsd)
#ifndef NO_CDF
      use b2mod_geo
     & , only : crx,cry,hx,hz,bb,vol,gs
      use b2mod_indirect
     & , only : leftix,leftiy,rightix,rightiy,
     &          topix,topiy,bottomix,bottomiy
      use b2mod_b2cmpa
     & , only : zamin,zamax
      use b2mod_residuals
     & , only : resco, resmo, reshe, reshi
      use b2mod_b2cmpa
     & , only : am
      use b2mod_constants
     & , only : mp,ev
#else
#     include "DIMENSIONS.F"
#endif
      implicit none
*   ..input arguments (unchanged on exit)
      integer nxd,nyd,nsd
#ifndef NO_CDF
#include <netcdf.inc>
*   ..local variables
      integer nstrad
      integer status
      integer ncid,nxid,nyid,n1id,n2id,n3id,n4id,nsid,nstraid,specid,
     .        eirnxid,eirnyid,eirnatmid,eirnmolid
      integer vdims(5),idum(1)
      integer jxi,jxa,jsep
      integer is
      character*13, allocatable :: species_name(:)
      integer use_eirene
      integer, allocatable :: za(:)
      real (kind=R8), save :: fac_hz = 1.0e0_R8
      real (kind=R8) :: rdum(1)
*   ..procedures
      intrinsic size
      external ipgeti, get_jsep, ipgetr, species, check_cdf_status

      write(*,*) 'Begin writing balance.nc file'

*   ..create the balance.nc file
      status=nf_create('balance.nc',NF_CLOBBER,ncid)
      call check_cdf_status(status)

*   ..define all of the dimensions that will be used
      status=nf_def_dim(ncid,'nx_plus2',nxd+2,nxid)
      call check_cdf_status(status)
      status=nf_def_dim(ncid,'ny_plus2',nyd+2,nyid)
      call check_cdf_status(status)
      status=nf_def_dim(ncid,'1',1,n1id)
      call check_cdf_status(status)
      status=nf_def_dim(ncid,'2',2,n2id)
      call check_cdf_status(status)
      status=nf_def_dim(ncid,'3',3,n3id)
      call check_cdf_status(status)
      status=nf_def_dim(ncid,'4',4,n4id)
      call check_cdf_status(status)
      status=nf_def_dim(ncid,'ns',nsd,nsid)
      call check_cdf_status(status)
      if (allocated(eirene_mc_papl_sna_bal)) then
        nstrad=size(eirene_mc_papl_sna_bal,3)
      else
        nstrad=1
      endif
      status=nf_def_dim(ncid,'nstra',nstrad,nstraid)
      call check_cdf_status(status)
      status=nf_def_dim(ncid,'speclng',13,specid)
      call check_cdf_status(status)
      if (allocated(dab2_av)) then
        status=nf_def_dim(ncid,'eirnx',size(pfluxa_av,1),eirnxid)
        call check_cdf_status(status)
        status=nf_def_dim(ncid,'eirny',size(pfluxa_av,2),eirnyid)
        call check_cdf_status(status)
        status=nf_def_dim(ncid,'natm',size(pfluxa_av,3),eirnatmid)
        call check_cdf_status(status)
        status=nf_def_dim(ncid,'nmol',size(pfluxm_av,3),eirnmolid)
        call check_cdf_status(status)
      endif

*   ..take dataset out of define mode
      status=nf_enddef(ncid)
      call check_cdf_status(status)

*   ..crx, cry:
      vdims(1)=nxid
      vdims(2)=nyid
      vdims(3)=n4id
      call write_cdf_double(ncid,'crx',crx,vdims,3)
      call write_cdf_double(ncid,'cry',cry,vdims,3)
*   ..bb,hx,hz,vol,gs,am:
      call write_cdf_double(ncid,'bb',bb,vdims,3)
      call write_cdf_double(ncid,'hx',hx,vdims,2)
      call write_cdf_double(ncid,'hz',hz,vdims,2)
      call write_cdf_double(ncid,'vol',vol,vdims,2)
      vdims(3)=n3id
      call write_cdf_double(ncid,'gs',gs,vdims,3)
      vdims(1)=nsid
      call write_cdf_double(ncid,'am',am,vdims,1)
*   ..mp, ev:
      vdims(1)=n1id
      rdum(1)=mp
      call write_cdf_double_single(ncid,'mp',rdum,vdims(1))
      rdum(1)=ev
      call write_cdf_double_single(ncid,'ev',rdum,vdims(1))
*   ..rightix, rightiy etc:
      vdims(1)=nxid
      call write_cdf_int(ncid,'leftix',leftix,vdims,2)
      call write_cdf_int(ncid,'leftiy',leftiy,vdims,2)
      call write_cdf_int(ncid,'rightix',rightix,vdims,2)
      call write_cdf_int(ncid,'rightiy',rightiy,vdims,2)
      call write_cdf_int(ncid,'topix',topix,vdims,2)
      call write_cdf_int(ncid,'topiy',topiy,vdims,2)
      call write_cdf_int(ncid,'bottomix',bottomix,vdims,2)
      call write_cdf_int(ncid,'bottomiy',bottomiy,vdims,2)
*   ..jxi, jxa, jsep:
      call get_jsep(nxd,nyd,jxi,jxa,jsep)
      vdims(1)=n1id
      idum(1)=jxi
      call write_cdf_int_single(ncid,'jxi',idum,vdims(1))
      idum(1)=jxa
      call write_cdf_int_single(ncid,'jxa',idum,vdims(1))
      idum(1)=jsep
      call write_cdf_int_single(ncid,'jsep',idum,vdims(1))
*   ..species string:
      vdims(1)=specid
      vdims(2)=nsid
      allocate(species_name(0:nsd-1))
      do is=0,nsd-1
        call species(is,species_name(is),.false.)
      enddo
      call write_cdf_text(ncid,'species',species_name(0:nsd-1),vdims,2)
*   ..b2mndr_eirene:
      use_eirene=0
      call ipgeti ('b2mndr_eirene', use_eirene)
      vdims(1)=n1id
      idum(1)=use_eirene
      call write_cdf_int_single(ncid,'b2mndr_eirene',idum,vdims(1))
*   ..fac_hz:
      call ipgetr ('b2mndr_hz', fac_hz)
      vdims(1)=n1id
      rdum(1)=fac_hz
      call write_cdf_double_single(ncid,'b2mndr_hz',rdum,vdims(1))
*   ..za:
      vdims(1)=nsid
      allocate(za(0:nsd-1))
      za(0:nsd-1)=nint((zamin(0:nsd-1)+zamax(0:nsd-1))/2.0_R8)
      call write_cdf_int(ncid,'za',za(0:nsd-1),vdims,1)
*   ..Arrays for particle balance:
*     Fluxes
      vdims(1)=nxid
      vdims(2)=nyid
      vdims(3)=n2id
      vdims(4)=n2id
      vdims(5)=nsid
      call write_cdf_double(ncid,'fna_pinch',
     &                   fna_pinch_av,vdims,5)
      call write_cdf_double(ncid,'fna_pll',
     &                   fna_pll_av,vdims,5)
      call write_cdf_double(ncid,'fna_drift',
     &                   fna_drift_av,vdims,5)
      call write_cdf_double(ncid,'fna_ch',
     &                   fna_ch_av,vdims,5)
      call write_cdf_double(ncid,'fna_nanom',
     &                   fna_nanom_av,vdims,5)
      call write_cdf_double(ncid,'fna_panom',
     &                   fna_panom_av,vdims,5)
      call write_cdf_double(ncid,'fna_pschused',
     &                   fna_pschused_av,vdims,5)
*     Sources
      vdims(3)=nsid
      call write_cdf_double(ncid,'b2stbr_phys_sna_bal',
     &                   b2stbr_phys_sna_bal_av,vdims,3)
      call write_cdf_double(ncid,'b2stbr_bas_sna_bal',
     &                   b2stbr_bas_sna_bal_av,vdims,3)
      call write_cdf_double(ncid,'b2stbr_first_flight_sna_bal',
     &                   b2stbr_first_flight_sna_bal_av,vdims,3)
      call write_cdf_double(ncid,'b2stbc_sna_bal',
     &                   b2stbc_sna_bal_av,vdims,3)
      call write_cdf_double(ncid,'b2stbm_sna_bal',
     &                   b2stbm_sna_bal_av,vdims,3)
      call write_cdf_double(ncid,'ext_sna_bal',
     &                   ext_sna_bal_av,vdims,3)
      call write_cdf_double(ncid,'b2stel_sna_ion_bal',
     &                   b2stel_sna_ion_bal_av,vdims,3)
      call write_cdf_double(ncid,'b2stel_sna_rec_bal',
     &                   b2stel_sna_rec_bal_av,vdims,3)
      call write_cdf_double(ncid,'b2stcx_sna_bal',
     &                   b2stcx_sna_bal_av,vdims,3)
      call write_cdf_double(ncid,'b2srsm_sna_bal',
     &                   b2srsm_sna_bal_av,vdims,3)
      call write_cdf_double(ncid,'b2srdt_sna_bal',
     &                   b2srdt_sna_bal_av,vdims,3)
      call write_cdf_double(ncid,'b2srst_sna_bal',
     &                   b2srst_sna_bal_av,vdims,3)
*     Residual
      call write_cdf_double(ncid,'resco',resco_av,vdims,3)
*   ..Arrays for momentum balance:
*     Fluxes
      vdims(3)=n2id
      vdims(4)=n2id
      vdims(5)=nsid
      call write_cdf_double(ncid,'fmo_flua',fmo_flua_av,vdims,5)
      call write_cdf_double(ncid,'fmo_cvsa',fmo_cvsa_av,vdims,5)
      call write_cdf_double(ncid,'fmo_b2nxfv',fmo_b2nxfv_av,vdims,5)
*     Sources
      vdims(3)=nsid
      call write_cdf_double(ncid,'b2stbr_phys_smo_bal',
     &                   b2stbr_phys_smo_bal_av,vdims,3)
      call write_cdf_double(ncid,'b2stbr_bas_smo_bal',
     &                   b2stbr_bas_smo_bal_av,vdims,3)
      call write_cdf_double(ncid,'b2stbc_smo_bal',
     &                   b2stbc_smo_bal_av,vdims,3)
      call write_cdf_double(ncid,'b2stbm_smo_bal',
     &                   b2stbm_smo_bal_av,vdims,3)
      call write_cdf_double(ncid,'ext_smo_bal',
     &                   ext_smo_bal_av,vdims,3)
      call write_cdf_double(ncid,'b2stel_smq_ion_bal',
     &                   b2stel_smq_ion_bal_av,vdims,3)
      call write_cdf_double(ncid,'b2stel_smq_rec_bal',
     &                   b2stel_smq_rec_bal_av,vdims,3)
      call write_cdf_double(ncid,'b2stcx_smq_bal',
     &                   b2stcx_smq_bal_av,vdims,3)
      call write_cdf_double(ncid,'b2srsm_smo_bal',
     &                   b2srsm_smo_bal_av,vdims,3)
      call write_cdf_double(ncid,'b2srdt_smo_bal',
     &                   b2srdt_smo_bal_av,vdims,3)
      call write_cdf_double(ncid,'b2srst_smo_bal',
     &                   b2srst_smo_bal_av,vdims,3)
      call write_cdf_double(ncid,'b2sifr_smoch_bal',
     &                   b2sifr_smoch_bal_av,vdims,3)
      call write_cdf_double(ncid,'b2sifr_smotf_bal',
     &                   b2sifr_smotf_bal_av,vdims,3)
      call write_cdf_double(ncid,'b2siav_smovh_bal',
     &                   b2siav_smovh_bal_av,vdims,3)
      call write_cdf_double(ncid,'b2siav_smovv_bal',
     &                   b2siav_smovv_bal_av,vdims,3)
      call write_cdf_double(ncid,'b2sicf_smo_bal',
     &                   b2sicf_smo_bal_av,vdims,3)
      call write_cdf_double(ncid,'b2sian_smo_bal',
     &                   b2sian_smo_bal_av,vdims,3)            !srv 25.09.17
      call write_cdf_double(ncid,'b2nxdv_smo_bal',
     &                   b2nxdv_smo_bal_av,vdims,3)
      call write_cdf_double(ncid,'b2sigp_smogp_bal',
     &                   b2sigp_smogp_bal_av,vdims,3)
      call write_cdf_double(ncid,'b2sigp_pstat_bal',
     &                   b2sigp_pstat_bal_av,vdims,3)
*     Residual
      call write_cdf_double(ncid,'resmo',resmo_av,vdims,3)
*   ..Arrays for electron heat balance:
*     Fluxes
      vdims(3)=n2id
      vdims(4)=n2id
      call write_cdf_double(ncid,'fhe_32',
     &                   fhe_32_av,vdims,4)
      call write_cdf_double(ncid,'fhe_52',
     &                   fhe_52_av,vdims,4)
      call write_cdf_double(ncid,'fhe_thermj',
     &                   fhe_thermj_av,vdims,4)
      call write_cdf_double(ncid,'fhe_cond',
     &                   fhe_cond_av,vdims,4)
      call write_cdf_double(ncid,'fhe_dia',
     &                   fhe_dia_av,vdims,4)
      call write_cdf_double(ncid,'fhe_ecrb',
     &                   fhe_ecrb_av,vdims,4)
      call write_cdf_double(ncid,'fhe_strange',
     &                   fhe_strange_av,vdims,4)
      call write_cdf_double(ncid,'fhe_pschused',
     &                   fhe_pschused_av,vdims,4)
*     Sources
      call write_cdf_double(ncid,'b2stbr_phys_she_bal',
     &                   b2stbr_phys_she_bal_av,vdims,2)
      call write_cdf_double(ncid,'b2stbr_bas_she_bal',
     &                   b2stbr_bas_she_bal_av,vdims,2)
      call write_cdf_double(ncid,'b2stbr_first_flight_she_bal',
     &                   b2stbr_first_flight_she_bal_av,vdims,2)
      call write_cdf_double(ncid,'b2stbc_she_bal',
     &                   b2stbc_she_bal_av,vdims,2)
      call write_cdf_double(ncid,'b2stbm_she_bal',
     &                   b2stbm_she_bal_av,vdims,2)
      call write_cdf_double(ncid,'ext_she_bal',
     &                   ext_she_bal_av,vdims,2)
      call write_cdf_double(ncid,'b2stel_she_bal',
     &                   b2stel_she_bal_av,vdims,2)
      call write_cdf_double(ncid,'b2srsm_she_bal',
     &                   b2srsm_she_bal_av,vdims,2)
      call write_cdf_double(ncid,'b2srdt_she_bal',
     &                   b2srdt_she_bal_av,vdims,2)
      call write_cdf_double(ncid,'b2srst_she_bal',
     &                   b2srst_she_bal_av,vdims,2)
      call write_cdf_double(ncid,'b2sihs_diae_bal',
     &                   b2sihs_diae_bal_av,vdims,2)
      call write_cdf_double(ncid,'b2sihs_divue_bal',
     &                   b2sihs_divue_bal_av,vdims,2)
      call write_cdf_double(ncid,'b2sihs_exbe_bal',
     &                   b2sihs_exbe_bal_av,vdims,2)
      call write_cdf_double(ncid,'b2sihs_joule_bal',
     &                   b2sihs_joule_bal_av,vdims,2)
      call write_cdf_double(ncid,'b2npht_shei_bal',
     &                   b2npht_shei_bal_av,vdims,2)
*     Residual
      call write_cdf_double(ncid,'reshe',reshe_av,vdims,2)
*   ..Arrays for ion heat balance:
*     Fluxes
      vdims(3)=n2id
      vdims(4)=n2id
      call write_cdf_double(ncid,'fhi_32',
     &                      fhi_32_av,vdims,4)
      call write_cdf_double(ncid,'fhi_52',
     &                      fhi_52_av,vdims,4)
      call write_cdf_double(ncid,'fhi_cond',
     &                      fhi_cond_av,vdims,4)
      call write_cdf_double(ncid,'fhi_dia',
     &                   fhi_dia_av,vdims,4)
      call write_cdf_double(ncid,'fhi_ecrb',
     &                   fhi_ecrb_av,vdims,4)
      call write_cdf_double(ncid,'fhi_strange',
     &                   fhi_strange_av,vdims,4)
      call write_cdf_double(ncid,'fhi_pschused',
     &                   fhi_pschused_av,vdims,4)
      call write_cdf_double(ncid,'fhi_inert',
     &                   fhi_inert_av,vdims,4)
      call write_cdf_double(ncid,'fhi_vispar',
     &                   fhi_vispar_av,vdims,4)
      call write_cdf_double(ncid,'fhi_visper',
     &                   fhi_visper_av,vdims,4)
      call write_cdf_double(ncid,'fhi_visq',
     &                   fhi_visq_av,vdims,4)
      call write_cdf_double(ncid,'fhi_anml',
     &                   fhi_anml_av,vdims,4)
      call write_cdf_double(ncid,'fhi_kevis',
     &                   fhi_kevis_av,vdims,4)
*     Sources
      call write_cdf_double(ncid,'b2stbr_phys_shi_bal',
     &                      b2stbr_phys_shi_bal_av,vdims,2)
      call write_cdf_double(ncid,'b2stbr_bas_shi_bal',
     &                      b2stbr_bas_shi_bal_av,vdims,2)
      call write_cdf_double(ncid,'b2stbr_first_flight_shi_bal',
     &                      b2stbr_first_flight_shi_bal_av,vdims,2)
      call write_cdf_double(ncid,'b2stbc_shi_bal',
     &                      b2stbc_shi_bal_av,vdims,2)
      call write_cdf_double(ncid,'b2stbm_shi_bal',
     &                      b2stbm_shi_bal_av,vdims,2)
      call write_cdf_double(ncid,'ext_shi_bal',
     &                      ext_shi_bal_av,vdims,2)
      call write_cdf_double(ncid,'b2stel_shi_ion_bal',
     &                      b2stel_shi_ion_bal_av,vdims,2)
      call write_cdf_double(ncid,'b2stel_shi_rec_bal',
     &                      b2stel_shi_rec_bal_av,vdims,2)
      call write_cdf_double(ncid,'b2stcx_shi_bal',
     &                      b2stcx_shi_bal_av,vdims,2)
      call write_cdf_double(ncid,'b2srsm_shi_bal',
     &                      b2srsm_shi_bal_av,vdims,2)
      call write_cdf_double(ncid,'b2srdt_shi_bal',
     &                      b2srdt_shi_bal_av,vdims,2)
      call write_cdf_double(ncid,'b2srst_shi_bal',
     &                      b2srst_shi_bal_av,vdims,2)
      call write_cdf_double(ncid,'b2sihs_diaa_bal',
     &                      b2sihs_diaa_bal_av,vdims,2)
      call write_cdf_double(ncid,'b2sihs_divua_bal',
     &                      b2sihs_divua_bal_av,vdims,2)
      call write_cdf_double(ncid,'b2sihs_exba_bal',
     &                      b2sihs_exba_bal_av,vdims,2)
      call write_cdf_double(ncid,'b2sihs_visa_bal',
     &                      b2sihs_visa_bal_av,vdims,2)
      call write_cdf_double(ncid,'b2sihs_fraa_bal',
     &                      b2sihs_fraa_bal_av,vdims,2)
*     Residual
      call write_cdf_double(ncid,'reshi',reshi_av,vdims,2)
*   ..Eirene sources
      if (allocated(eirene_mc_papl_sna_bal)) then
        vdims(3)=nsid
        vdims(4)=nstraid
        call write_cdf_double(ncid,'eirene_mc_papl_sna_bal',
     &                     eirene_mc_papl_sna_bal_av,vdims,4)
        call write_cdf_double(ncid,'eirene_mc_pmpl_sna_bal',
     &                     eirene_mc_pmpl_sna_bal_av,vdims,4)
        call write_cdf_double(ncid,'eirene_mc_pipl_sna_bal',
     &                     eirene_mc_pipl_sna_bal_av,vdims,4)
        call write_cdf_double(ncid,'eirene_mc_pppl_sna_bal',
     &                     eirene_mc_pppl_sna_bal_av,vdims,4)
        call write_cdf_double(ncid,'eirene_mc_core_sna_bal',
     &                     eirene_mc_core_sna_bal_av,vdims,4)
        call write_cdf_double(ncid,'eirene_mc_mapl_smo_bal',
     &                     eirene_mc_mapl_smo_bal_av,vdims,4)
        call write_cdf_double(ncid,'eirene_mc_mmpl_smo_bal',
     &                     eirene_mc_mmpl_smo_bal_av,vdims,4)
        call write_cdf_double(ncid,'eirene_mc_mipl_smo_bal',
     &                     eirene_mc_mipl_smo_bal_av,vdims,4)
        call write_cdf_double(ncid,'eirene_mc_cppv_smo_bal',
     &                     eirene_mc_cppv_smo_bal_av,vdims,4)
        vdims(3)=nstraid
        call write_cdf_double(ncid,'eirene_mc_eael_she_bal',
     &                     eirene_mc_eael_she_bal_av,vdims,3)
        call write_cdf_double(ncid,'eirene_mc_emel_she_bal',
     &                     eirene_mc_emel_she_bal_av,vdims,3)
        call write_cdf_double(ncid,'eirene_mc_eiel_she_bal',
     &                     eirene_mc_eiel_she_bal_av,vdims,3)
        call write_cdf_double(ncid,'eirene_mc_epel_she_bal',
     &                     eirene_mc_epel_she_bal_av,vdims,3)
        call write_cdf_double(ncid,'eirene_mc_eapl_shi_bal',
     &                     eirene_mc_eapl_shi_bal_av,vdims,3)
        call write_cdf_double(ncid,'eirene_mc_empl_shi_bal',
     &                     eirene_mc_empl_shi_bal_av,vdims,3)
        call write_cdf_double(ncid,'eirene_mc_eipl_shi_bal',
     &                     eirene_mc_eipl_shi_bal_av,vdims,3)
        call write_cdf_double(ncid,'eirene_mc_eppl_shi_bal',
     &                     eirene_mc_eppl_shi_bal_av,vdims,3)
      else
        write(*,*) 'skipping Eirene sources'
      endif
*   ..Additional averaged quantities
      vdims(3)=nsid
      call write_cdf_double(ncid,'ne',ne_av,vdims,2)
      call write_cdf_double(ncid,'na',na_av,vdims,3)
      call write_cdf_double(ncid,'ua',ua_av,vdims,3)
      call write_cdf_double(ncid,'po',po_av,vdims,2)
      call write_cdf_double(ncid,'te',te_av,vdims,2)
      call write_cdf_double(ncid,'ti',ti_av,vdims,2)
      call write_cdf_double(ncid,'rpt',rpt_av,vdims,3)
      call write_cdf_double(ncid,'kinrgy',kinrgy_av,vdims,3)
      vdims(3)=n2id
      call write_cdf_double(ncid,'fne',fne_av,vdims,3)
      if (allocated(dab2_av)) then
        vdims(1)=eirnxid
        vdims(2)=eirnyid
        vdims(3)=eirnatmid
        call write_cdf_double(ncid,'dab2',dab2_av,vdims,3)
        call write_cdf_double(ncid,'tab2',tab2_av,vdims,3)
        call write_cdf_double(ncid,'rfluxa',rfluxa_av,vdims,3)
        call write_cdf_double(ncid,'refluxa',refluxa_av,vdims,3)
        call write_cdf_double(ncid,'pfluxa',pfluxa_av,vdims,3)
        call write_cdf_double(ncid,'pefluxa',pefluxa_av,vdims,3)
        vdims(3)=eirnmolid
        call write_cdf_double(ncid,'dmb2',dmb2_av,vdims,3)
        call write_cdf_double(ncid,'tmb2',tmb2_av,vdims,3)
        call write_cdf_double(ncid,'rfluxm',rfluxm_av,vdims,3)
        call write_cdf_double(ncid,'refluxm',refluxm_av,vdims,3)
        call write_cdf_double(ncid,'pfluxm',pfluxm_av,vdims,3)
        call write_cdf_double(ncid,'pefluxm',pefluxm_av,vdims,3)
      endif

      status=nf_close(ncid)
      call check_cdf_status(status)

      write(*,*) 'balance.nc written'
#else
*   ..procedures
      external ipgeti, xertst

      call xertst (nxd.le.DEF_NXD, 'faulty parameter nxd')
      call xertst (nyd.le.DEF_NYD, 'faulty parameter nyd')
      call xertst (nsd.le.DEF_NSD, 'faulty parameter nsd')
      write(*,*) 'Compilation with NO_CDF: balance.nc file not written'
#endif

      return
      end subroutine write_balance

      subroutine read_balance
      implicit none
#ifndef NO_CDF
#     include <netcdf.inc>
      integer :: ncid, ndims, nvars, natts, unlimid, status, xtype
      integer :: ivar, idim
      integer :: len, vdims(NF_MAX_VAR_DIMS), start(NF_MAX_VAR_DIMS)
      real (kind=R8), allocatable :: data2(:,:), data3(:,:,:),
     ,                               data4(:,:,:,:)
      character*(NF_MAX_NAME) :: varname, dimname
      character*256 filename
      logical :: exists
      external find_file, check_cdf_status

      write(*,*) 'Begin reading balance.nc file'

c access the balance.nc file
      filename='balance.nc'
      call find_file(filename,exists)
      status=nf_open(filename,NF_NOWRITE,ncid)
      call check_cdf_status(status)

      status=nf_inq(ncid,ndims,nvars,natts,unlimid)
      call check_cdf_status(status)

c loop over all variables in the file
c only read a subset of the data as needed by the calling routine(s)
      do ivar=1,nvars
        status=nf_inq_var(ncid,ivar,varname,xtype,ndims,vdims,natts)
        call check_cdf_status(status)
        do idim=1,ndims
          status=nf_inq_dim(ncid,vdims(idim),dimname,len)
          call check_cdf_status(status)
          vdims(idim)=len
          start(idim)=1
        enddo
        select case (varname)
        case ('b2stel_she_ion_bal')
          allocate(data2(vdims(1),vdims(2)))
          status=nf_get_vara_double(ncid,ivar,start,vdims,data2)
          call check_cdf_status(status)
          b2stel_she_ion_bal = data2
          deallocate(data2)
        case ('b2stel_she_rec_bal')
          allocate(data2(vdims(1),vdims(2)))
          status=nf_get_vara_double(ncid,ivar,start,vdims,data2)
          call check_cdf_status(status)
          b2stel_she_rec_bal = data2
          deallocate(data2)
        case ('b2stel_shi_ion_bal')
          allocate(data2(vdims(1),vdims(2)))
          status=nf_get_vara_double(ncid,ivar,start,vdims,data2)
          call check_cdf_status(status)
          b2stel_shi_ion_bal = data2
          deallocate(data2)
        case ('b2stel_shi_rec_bal')
          allocate(data2(vdims(1),vdims(2)))
          status=nf_get_vara_double(ncid,ivar,start,vdims,data2)
          call check_cdf_status(status)
          b2stel_shi_rec_bal = data2
          deallocate(data2)
        case ('b2npht_shei_bal')
          allocate(data2(vdims(1),vdims(2)))
          status=nf_get_vara_double(ncid,ivar,start,vdims,data2)
          call check_cdf_status(status)
          b2npht_shei_bal = data2
          deallocate(data2)
        case ('eirene_mc_pael_sne_bal')
          allocate(data3(vdims(1),vdims(2),vdims(3)))
          status=nf_get_vara_double(ncid,ivar,start,vdims,data3)
          call check_cdf_status(status)
          eirene_mc_pael_sne_bal = data3
          deallocate(data3)
        case ('eirene_mc_pmel_sne_bal')
          allocate(data3(vdims(1),vdims(2),vdims(3)))
          status=nf_get_vara_double(ncid,ivar,start,vdims,data3)
          call check_cdf_status(status)
          eirene_mc_pmel_sne_bal = data3
          deallocate(data3)
        case ('eirene_mc_eael_she_bal')
          allocate(data3(vdims(1),vdims(2),vdims(3)))
          status=nf_get_vara_double(ncid,ivar,start,vdims,data3)
          call check_cdf_status(status)
          eirene_mc_eael_she_bal = data3
          deallocate(data3)
        case ('eirene_mc_emel_she_bal')
          allocate(data3(vdims(1),vdims(2),vdims(3)))
          status=nf_get_vara_double(ncid,ivar,start,vdims,data3)
          call check_cdf_status(status)
          eirene_mc_emel_she_bal = data3
          deallocate(data3)
        case ('eirene_mc_eapl_shi_bal')
          allocate(data3(vdims(1),vdims(2),vdims(3)))
          status=nf_get_vara_double(ncid,ivar,start,vdims,data3)
          call check_cdf_status(status)
          eirene_mc_eapl_shi_bal = data3
          deallocate(data3)
        case ('eirene_mc_empl_shi_bal')
          allocate(data3(vdims(1),vdims(2),vdims(3)))
          status=nf_get_vara_double(ncid,ivar,start,vdims,data3)
          call check_cdf_status(status)
          eirene_mc_empl_shi_bal = data3
          deallocate(data3)
        case ('eirene_mc_papl_sna_bal')
          allocate(data4(vdims(1),vdims(2),vdims(3),vdims(4)))
          status=nf_get_vara_double(ncid,ivar,start,vdims,data4)
          call check_cdf_status(status)
          eirene_mc_papl_sna_bal = data4
          deallocate(data4)
        case ('eirene_mc_pmpl_sna_bal')
          allocate(data4(vdims(1),vdims(2),vdims(3),vdims(4)))
          status=nf_get_vara_double(ncid,ivar,start,vdims,data4)
          call check_cdf_status(status)
          eirene_mc_pmpl_sna_bal = data4
          deallocate(data4)
        case ('eirene_mc_paat_sna_bal')
          allocate(data4(vdims(1),vdims(2),vdims(3),vdims(4)))
          status=nf_get_vara_double(ncid,ivar,start,vdims,data4)
          call check_cdf_status(status)
          eirene_mc_paat_sna_bal = data4
          deallocate(data4)
        case ('eirene_mc_pmat_sna_bal')
          allocate(data4(vdims(1),vdims(2),vdims(3),vdims(4)))
          status=nf_get_vara_double(ncid,ivar,start,vdims,data4)
          call check_cdf_status(status)
          eirene_mc_pmat_sna_bal = data4
          deallocate(data4)
        case ('eirene_mc_piat_sna_bal')
          allocate(data4(vdims(1),vdims(2),vdims(3),vdims(4)))
          status=nf_get_vara_double(ncid,ivar,start,vdims,data4)
          call check_cdf_status(status)
          eirene_mc_piat_sna_bal = data4
          deallocate(data4)
        case ('eirene_mc_paml_sna_bal')
          allocate(data4(vdims(1),vdims(2),vdims(3),vdims(4)))
          status=nf_get_vara_double(ncid,ivar,start,vdims,data4)
          call check_cdf_status(status)
          eirene_mc_paml_sna_bal = data4
          deallocate(data4)
        case ('eirene_mc_pmml_sna_bal')
          allocate(data4(vdims(1),vdims(2),vdims(3),vdims(4)))
          status=nf_get_vara_double(ncid,ivar,start,vdims,data4)
          call check_cdf_status(status)
          eirene_mc_pmml_sna_bal = data4
          deallocate(data4)
        case ('eirene_mc_piml_sna_bal')
          allocate(data4(vdims(1),vdims(2),vdims(3),vdims(4)))
          status=nf_get_vara_double(ncid,ivar,start,vdims,data4)
          call check_cdf_status(status)
          eirene_mc_piml_sna_bal = data4
          deallocate(data4)
        case ('eirene_mc_paio_sna_bal')
          allocate(data4(vdims(1),vdims(2),vdims(3),vdims(4)))
          status=nf_get_vara_double(ncid,ivar,start,vdims,data4)
          call check_cdf_status(status)
          eirene_mc_paio_sna_bal = data4
          deallocate(data4)
        case ('eirene_mc_pmio_sna_bal')
          allocate(data4(vdims(1),vdims(2),vdims(3),vdims(4)))
          status=nf_get_vara_double(ncid,ivar,start,vdims,data4)
          call check_cdf_status(status)
          eirene_mc_pmio_sna_bal = data4
          deallocate(data4)
        case ('eirene_mc_piio_sna_bal')
          allocate(data4(vdims(1),vdims(2),vdims(3),vdims(4)))
          status=nf_get_vara_double(ncid,ivar,start,vdims,data4)
          call check_cdf_status(status)
          eirene_mc_piio_sna_bal = data4
          deallocate(data4)
        case ('eirene_mc_mapl_smo_bal')
          allocate(data4(vdims(1),vdims(2),vdims(3),vdims(4)))
          status=nf_get_vara_double(ncid,ivar,start,vdims,data4)
          call check_cdf_status(status)
          eirene_mc_mapl_smo_bal = data4
          deallocate(data4)
        case ('eirene_mc_mmpl_smo_bal')
          allocate(data4(vdims(1),vdims(2),vdims(3),vdims(4)))
          status=nf_get_vara_double(ncid,ivar,start,vdims,data4)
          call check_cdf_status(status)
          eirene_mc_mmpl_smo_bal = data4
          deallocate(data4)
        case default
          continue
        end select
      enddo

      status=nf_close(ncid)
      call check_cdf_status(status)
#endif

      return
      end subroutine read_balance
#endif

#ifndef NO_CDF
      subroutine write_cdf_double(ncid,varname,var,vdims,ndim)
      implicit none
#include <netcdf.inc>
*   ..input arguments
      integer :: ncid,vdims(*),ndim
      character*(*) :: varname
      real (kind=R8) :: var(*)
*   ..local variables
      integer status,varid
*   ..procedures
      external check_cdf_status

*   ..put NetCDF dataset into define mode
      status=nf_redef(ncid)
      call check_cdf_status(status)
*   ..define this variable
      status=nf_def_var(ncid,varname,NF_DOUBLE,ndim,vdims,varid)
      call check_cdf_status(status)
*   ..take dataset out of define mode
      status=nf_enddef(ncid)
      call check_cdf_status(status)
*   ..before writing the data
      status=nf_put_var_double(ncid,varid,var)
      call check_cdf_status(status)

      return
      end subroutine write_cdf_double

      subroutine write_cdf_double_single(ncid,varname,var,vdim)
      implicit none
#include <netcdf.inc>
*   ..input arguments
      integer :: ncid
      integer :: vdim(1)
      character*(*) :: varname
      real (kind=R8) :: var(1)
*   ..local variables
      integer status,varid
*   ..procedures
      external check_cdf_status

*   ..put NetCDF dataset into define mode
      status=nf_redef(ncid)
      call check_cdf_status(status)
*   ..define this variable
      status=nf_def_var(ncid,varname,NF_DOUBLE,1,vdim,varid)
      call check_cdf_status(status)
*   ..take dataset out of define mode
      status=nf_enddef(ncid)
      call check_cdf_status(status)
*   ..before writing the data
      status=nf_put_var_double(ncid,varid,var)
      call check_cdf_status(status)

      return
      end subroutine write_cdf_double_single

      subroutine write_cdf_int(ncid,varname,var,vdims,ndim)
      implicit none
#include <netcdf.inc>
*   ..input arguments
      integer :: ncid,vdims(*),ndim
      character*(*) :: varname
      integer :: var(*)
*   ..local variables
      integer status,varid
*   ..procedures
      external check_cdf_status

*   ..put NetCDF dataset into define mode
      status=nf_redef(ncid)
      call check_cdf_status(status)
*   ..define this variable
      status=nf_def_var(ncid,varname,NF_INT,ndim,vdims,varid)
      call check_cdf_status(status)
*   ..take dataset out of define mode
      status=nf_enddef(ncid)
      call check_cdf_status(status)
*   ..before writing the data
      status=nf_put_var_int(ncid,varid,var)
      call check_cdf_status(status)

      return
      end subroutine write_cdf_int

      subroutine write_cdf_int_single(ncid,varname,var,vdim)
      implicit none
#include <netcdf.inc>
*   ..input arguments
      integer :: ncid
      integer :: vdim(1),var(1)
      character*(*) :: varname
*   ..local variables
      integer status,varid
*   ..procedures
      external check_cdf_status

*   ..put NetCDF dataset into define mode
      status=nf_redef(ncid)
      call check_cdf_status(status)
*   ..define this variable
      status=nf_def_var(ncid,varname,NF_INT,1,vdim,varid)
      call check_cdf_status(status)
*   ..take dataset out of define mode
      status=nf_enddef(ncid)
      call check_cdf_status(status)
*   ..before writing the data
      status=nf_put_var_int(ncid,varid,var)
      call check_cdf_status(status)

      return
      end subroutine write_cdf_int_single

      subroutine write_cdf_text(ncid,varname,var,vdims,ndim)
      implicit none
#include <netcdf.inc>
*   ..input arguments
      integer :: ncid,vdims(*),ndim
      character*(*) :: varname
      character*(*) :: var(*)
*   ..local variables
      integer status,varid
*   ..procedures
      external check_cdf_status

*   ..put NetCDF dataset into define mode
      status=nf_redef(ncid)
      call check_cdf_status(status)
*   ..define this variable
      status=nf_def_var(ncid,varname,NF_CHAR,ndim,vdims,varid)
      call check_cdf_status(status)
*   ..take dataset out of define mode
      status=nf_enddef(ncid)
      call check_cdf_status(status)
*   ..before writing the data
      status=nf_put_var_text(ncid,varid,var)
      call check_cdf_status(status)

      return
      end subroutine write_cdf_text
#endif

      end module b2mod_balance
