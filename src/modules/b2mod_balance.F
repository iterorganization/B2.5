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
     *  fna_pinch(:,:,:,:),
     *  fna_pll(:,:,:,:),
     *  fna_drift(:,:,:,:),
     *  fna_ch(:,:,:,:),
     *  fna_nanom(:,:,:,:),
     *  fna_panom(:,:,:,:),
     *  fna_pschused(:,:,:,:),
     *  fna_tot(:,:,:,:)
*     Sources
      real (kind=R8), allocatable, save ::
     *  b2stbr_phys_sna0to1(:,:,:,:),
     *  b2stbr_phys_sna_bal(:,:,:),
     *  b2stbr_bas_sna0to1(:,:,:,:),
     *  b2stbr_bas_sna_bal(:,:,:),
     *  b2stbr_first_flight_sna0to1(:,:,:,:),
     *  b2stbr_first_flight_sna_bal(:,:,:),
     *  b2stbc_sna0to1(:,:,:,:),
     *  b2stbc_sna_bal(:,:,:),
     *  eirene_mc_pael_sne0to1(:,:,:,:),
     *  eirene_mc_papl_sna0to1(:,:,:,:,:),
     *  eirene_mc_pmel_sne0to1(:,:,:,:),
     *  eirene_mc_pmpl_sna0to1(:,:,:,:,:),
     *  eirene_mc_pipl_sna0to1(:,:,:,:,:),
     *  eirene_mc_pppl_sna0to1(:,:,:,:,:),
     *  eirene_mc_core_sna0to1(:,:,:,:),
     *  eirene_mc_pael_sne_bal(:,:,:),
     *  eirene_mc_papl_sna_bal(:,:,:,:),
     *  eirene_mc_pmel_sne_bal(:,:,:),
     *  eirene_mc_pmpl_sna_bal(:,:,:,:),
     *  eirene_mc_pipl_sna_bal(:,:,:,:),
     *  eirene_mc_pppl_sna_bal(:,:,:,:),
     *  eirene_mc_paat_sna_bal(:,:,:,:),
     *  eirene_mc_pmat_sna_bal(:,:,:,:),
     *  eirene_mc_piat_sna_bal(:,:,:,:),
     *  eirene_mc_paml_sna_bal(:,:,:,:),
     *  eirene_mc_pmml_sna_bal(:,:,:,:),
     *  eirene_mc_piml_sna_bal(:,:,:,:),
     *  eirene_mc_paio_sna_bal(:,:,:,:),
     *  eirene_mc_pmio_sna_bal(:,:,:,:),
     *  eirene_mc_piio_sna_bal(:,:,:,:),
     *  eirene_mc_core_sna_bal(:,:,:),
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
     *  b2srst_sna_bal(:,:,:),
     *  tot_sna0to1(:,:,:,:),
     *  tot_sna_bal(:,:,:)
*     Time-averaged quantities
      real (kind=R8), allocatable, save ::
     *  fna_pinch_av(:,:,:,:),
     *  fna_pll_av(:,:,:,:),
     *  fna_drift_av(:,:,:,:),
     *  fna_ch_av(:,:,:,:),
     *  fna_nanom_av(:,:,:,:),
     *  fna_panom_av(:,:,:,:),
     *  fna_pschused_av(:,:,:,:),
     *  fna_tot_av(:,:,:,:),
     *  b2stbr_phys_sna_bal_av(:,:,:),
     *  b2stbr_bas_sna_bal_av(:,:,:),
     *  b2stbr_first_flight_sna_bal_av(:,:,:),
     *  b2stbc_sna_bal_av(:,:,:),
     *  eirene_mc_pael_sne_bal_av(:,:,:),
     *  eirene_mc_papl_sna_bal_av(:,:,:,:),
     *  eirene_mc_pmel_sne_bal_av(:,:,:),
     *  eirene_mc_pmpl_sna_bal_av(:,:,:,:),
     *  eirene_mc_pipl_sna_bal_av(:,:,:,:),
     *  eirene_mc_pppl_sna_bal_av(:,:,:,:),
     *  eirene_mc_paat_sna_bal_av(:,:,:,:),
     *  eirene_mc_pmat_sna_bal_av(:,:,:,:),
     *  eirene_mc_piat_sna_bal_av(:,:,:,:),
     *  eirene_mc_paml_sna_bal_av(:,:,:,:),
     *  eirene_mc_pmml_sna_bal_av(:,:,:,:),
     *  eirene_mc_piml_sna_bal_av(:,:,:,:),
     *  eirene_mc_paio_sna_bal_av(:,:,:,:),
     *  eirene_mc_pmio_sna_bal_av(:,:,:,:),
     *  eirene_mc_piio_sna_bal_av(:,:,:,:),
     *  eirene_mc_core_sna_bal_av(:,:,:),
     *  b2stbm_sna_bal_av(:,:,:),
     *  ext_sna_bal_av(:,:,:),
     *  b2stel_sna_ion_bal_av(:,:,:),
     *  b2stel_sna_rec_bal_av(:,:,:),
     *  b2stcx_sna_bal_av(:,:,:),
     *  b2srsm_sna_bal_av(:,:,:),
     *  b2srdt_sna_bal_av(:,:,:),
     *  b2srst_sna_bal_av(:,:,:),
     *  tot_sna_bal_av(:,:,:),
     *  resco_av(:,:,:)
*   ..Arrays for momentum balance:
*     Fluxes
      real (kind=R8), allocatable, save ::
     *  fmo_flua(:,:,:,:),
     *  fmo_cvsa(:,:,:,:),
     *  fmo_hybr(:,:,:,:),
     *  fmo_b2nxfv(:,:,:,:),
     *  fmo_tot(:,:,:,:)
*     Sources
      real (kind=R8), allocatable, save ::
     *  b2stbr_phys_smo0to3(:,:,:,:),
     *  b2stbr_phys_smo_bal(:,:,:),
     *  b2stbr_bas_smo0to3(:,:,:,:),
     *  b2stbr_bas_smo_bal(:,:,:),
     *  b2stbc_smo0to3(:,:,:,:),
     *  b2stbc_smo_bal(:,:,:),
     *  eirene_mc_mapl_smo0to3(:,:,:,:,:),
     *  eirene_mc_mmpl_smo0to3(:,:,:,:,:),
     *  eirene_mc_mipl_smo0to3(:,:,:,:,:),
     *  eirene_mc_mppl_smo0to3(:,:,:,:,:),
     *  eirene_mc_mapl_smo_bal(:,:,:,:),
     *  eirene_mc_mmpl_smo_bal(:,:,:,:),
     *  eirene_mc_mipl_smo_bal(:,:,:,:),
     *  eirene_mc_mppl_smo_bal(:,:,:,:),
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
     *  b2sifr_smotf_ehxp0to3(:,:,:,:),
     *  b2sifr_smotf_cthe0to3(:,:,:,:),
     *  b2sifr_smotf_cthi0to3(:,:,:,:),
     *  b2sifr_smotf_ehxp_bal(:,:,:),
     *  b2sifr_smotf_cthe_bal(:,:,:),
     *  b2sifr_smotf_cthi_bal(:,:,:),
     *  b2sifr_smofrea0to3(:,:,:,:),
     *  b2sifr_smofrea_bal(:,:,:),
     *  b2sifr_smofria0to3(:,:,:,:),
     *  b2sifr_smofria_bal(:,:,:),
     *  b2sifr_smotfea0to3(:,:,:,:),
     *  b2sifr_smotfea_bal(:,:,:),
     *  b2sifr_smotfia0to3(:,:,:,:),
     *  b2sifr_smotfia_bal(:,:,:),
     *  b2siav_smovh0to3(:,:,:,:),
     *  b2siav_smovh_bal(:,:,:),
     *  b2siav_smovv0to3(:,:,:,:),
     *  b2siav_smovv_bal(:,:,:),
     *  b2siav_smovi0to3(:,:,:,:),             !som 13.08.21
     *  b2siav_smovi_bal(:,:,:),               !som 13.08.21
     *  b2sicf_smo0to3(:,:,:,:),
     *  b2sicf_smo_bal(:,:,:),
     *  b2sian_smo0to3(:,:,:,:),               !srv 25.09.17
     *  b2sian_smo_bal(:,:,:),                 !srv 25.09.17
     *  b2nxdv_smo0to3(:,:,:,:),
     *  b2nxdv_smo_bal(:,:,:),
     *  b2sigp_smogp0to3(:,:,:,:),
     *  b2sigp_smogp_bal(:,:,:),
     *  b2sigp_smogpi0to3(:,:,:,:),
     *  b2sigp_smogpi_bal(:,:,:),
     *  b2sigp_smogpe0to3(:,:,:,:),
     *  b2sigp_smogpe_bal(:,:,:),
     *  b2sigp_smogpgr0to3(:,:,:,:),
     *  b2sigp_smogpgr_bal(:,:,:),
     *  b2sigp_pstat0to3(:,:,:,:),
     *  b2sigp_pstat_bal(:,:,:),
     *  b2sigp_pstati0to3(:,:,:,:),
     *  b2sigp_pstati_bal(:,:,:),
     *  b2sigp_pstate0to3(:,:,:),
     *  b2sigp_pstate_bal(:,:),
     *  tot_smo0to3(:,:,:,:),
     *  tot_smo_bal(:,:,:)
*     Time-averaged quantities
      real (kind=R8), allocatable, save ::
     *  fmo_flua_av(:,:,:,:),
     *  fmo_cvsa_av(:,:,:,:),
     *  fmo_hybr_av(:,:,:,:),
     *  fmo_b2nxfv_av(:,:,:,:),
     *  fmo_tot_av(:,:,:,:),
     *  b2stbr_phys_smo_bal_av(:,:,:),
     *  b2stbr_bas_smo_bal_av(:,:,:),
     *  b2stbc_smo_bal_av(:,:,:),
     *  eirene_mc_mapl_smo_bal_av(:,:,:,:),
     *  eirene_mc_mmpl_smo_bal_av(:,:,:,:),
     *  eirene_mc_mipl_smo_bal_av(:,:,:,:),
     *  eirene_mc_mppl_smo_bal_av(:,:,:,:),
     *  b2stbm_smo_bal_av(:,:,:),
     *  ext_smo_bal_av(:,:,:),
     *  b2stel_smq_ion_bal_av(:,:,:),
     *  b2stel_smq_rec_bal_av(:,:,:),
     *  b2stcx_smq_bal_av(:,:,:),
     *  b2srsm_smo_bal_av(:,:,:),
     *  b2srdt_smo_bal_av(:,:,:),
     *  b2srst_smo_bal_av(:,:,:),
     *  b2sifr_smoch_bal_av(:,:,:),
     *  b2sifr_smotf_ehxp_bal_av(:,:,:),
     *  b2sifr_smotf_cthe_bal_av(:,:,:),
     *  b2sifr_smotf_cthi_bal_av(:,:,:),
     *  b2sifr_smofrea_bal_av(:,:,:),
     *  b2sifr_smofria_bal_av(:,:,:),
     *  b2sifr_smotfea_bal_av(:,:,:),
     *  b2sifr_smotfia_bal_av(:,:,:),
     *  b2siav_smovh_bal_av(:,:,:),
     *  b2siav_smovv_bal_av(:,:,:),
     *  b2siav_smovi_bal_av(:,:,:),             !som 13.08.21
     *  b2sicf_smo_bal_av(:,:,:),
     *  b2sian_smo_bal_av(:,:,:),              !srv 25.07.17
     *  b2nxdv_smo_bal_av(:,:,:),
     *  b2sigp_smogp_bal_av(:,:,:),
     *  b2sigp_smogpi_bal_av(:,:,:),
     *  b2sigp_smogpe_bal_av(:,:,:),
     *  b2sigp_smogpgr_bal_av(:,:,:),
     *  b2sigp_pstat_bal_av(:,:,:),
     *  b2sigp_pstati_bal_av(:,:,:),
     *  b2sigp_pstate_bal_av(:,:),
     *  tot_smo_bal_av(:,:,:),
     *  resmo_av(:,:,:)
*   ..Arrays for electron heat balance:
*     Fluxes
      real (kind=R8), allocatable, save ::
     *  fhe_32(:,:,:),
     *  fhe_52(:,:,:),
     *  fhe_thermj(:,:,:),
     *  fhe_cond(:,:,:),
     *  fhe_dia(:,:,:),
     *  fhe_ecrb(:,:,:),
     *  fhe_strange(:,:,:),
     *  fhe_pschused(:,:,:)
*     Sources
      real (kind=R8), allocatable, save ::
     *  b2stbr_phys_she0to3(:,:,:),
     *  b2stbr_phys_she_bal(:,:),
     *  b2stbr_bas_she0to3(:,:,:),
     *  b2stbr_bas_she_bal(:,:),
     *  b2stbr_first_flight_she0to3(:,:,:),
     *  b2stbr_first_flight_she_bal(:,:),
     *  b2stbc_she0to3(:,:,:),
     *  b2stbc_she_bal(:,:),
     *  eirene_mc_eael_she0to3(:,:,:,:),
     *  eirene_mc_emel_she0to3(:,:,:,:),
     *  eirene_mc_eiel_she0to3(:,:,:,:),
     *  eirene_mc_epel_she0to3(:,:,:,:),
     *  eirene_mc_eael_she_bal(:,:,:),
     *  eirene_mc_emel_she_bal(:,:,:),
     *  eirene_mc_eiel_she_bal(:,:,:),
     *  eirene_mc_epel_she_bal(:,:,:),
     *  b2stbm_she0to3(:,:,:),
     *  b2stbm_she_bal(:,:),
     *  ext_she0to3(:,:,:),
     *  ext_she_bal(:,:),
     *  b2stel_she0to3(:,:,:,:),
     *  b2stel_she_bal(:,:,:),
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
*     Time-averaged quantities
      real (kind=R8), allocatable, save ::
     *  fhe_32_av(:,:,:),
     *  fhe_52_av(:,:,:),
     *  fhe_thermj_av(:,:,:),
     *  fhe_cond_av(:,:,:),
     *  fhe_dia_av(:,:,:),
     *  fhe_ecrb_av(:,:,:),
     *  fhe_strange_av(:,:,:),
     *  fhe_pschused_av(:,:,:),
     *  b2stbr_phys_she_bal_av(:,:),
     *  b2stbr_bas_she_bal_av(:,:),
     *  b2stbr_first_flight_she_bal_av(:,:),
     *  b2stbc_she_bal_av(:,:),
     *  eirene_mc_eael_she_bal_av(:,:,:),
     *  eirene_mc_emel_she_bal_av(:,:,:),
     *  eirene_mc_eiel_she_bal_av(:,:,:),
     *  eirene_mc_epel_she_bal_av(:,:,:),
     *  b2stbm_she_bal_av(:,:),
     *  ext_she_bal_av(:,:),
     *  b2stel_she_bal_av(:,:,:),
     *  b2srsm_she_bal_av(:,:),
     *  b2srdt_she_bal_av(:,:),
     *  b2srst_she_bal_av(:,:),
     *  b2sihs_diae_bal_av(:,:),
     *  b2sihs_divue_bal_av(:,:),
     *  b2sihs_exbe_bal_av(:,:),
     *  b2sihs_joule_bal_av(:,:),
     *  b2npht_shei_bal_av(:,:),
     *  reshe_av(:,:)
*   ..Arrays for ion heat balance:
*     Fluxes
      real (kind=R8), allocatable, save ::
     *  fhi_32(:,:,:),
     *  fhi_52(:,:,:),
     *  fhi_cond(:,:,:),
     *  fhi_dia(:,:,:),
     *  fhi_ecrb(:,:,:),
     *  fhi_strange(:,:,:),
     *  fhi_pschused(:,:,:),
     *  fhi_inert(:,:,:),
     *  fhi_vispar(:,:,:),
     *  fhi_visper(:,:,:),
     *  fhi_visq(:,:,:),
     *  fhi_anml(:,:,:),
     *  fhi_kevis(:,:,:),
     *  fhi_vzh(:,:,:)
*     Source
      real (kind=R8), allocatable, save ::
     *  b2stbr_phys_shi0to3(:,:,:),
     *  b2stbr_phys_shi_bal(:,:),
     *  b2stbr_bas_shi0to3(:,:,:),
     *  b2stbr_bas_shi_bal(:,:),
     *  b2stbr_first_flight_shi0to3(:,:,:),
     *  b2stbr_first_flight_shi_bal(:,:),
     *  b2stbc_shi0to3(:,:,:),
     *  b2stbc_shi_bal(:,:),
     *  eirene_mc_eapl_shi0to3(:,:,:,:),
     *  eirene_mc_empl_shi0to3(:,:,:,:),
     *  eirene_mc_eipl_shi0to3(:,:,:,:),
     *  eirene_mc_eppl_shi0to3(:,:,:,:),
     *  eirene_mc_eapl_shi_bal(:,:,:),
     *  eirene_mc_empl_shi_bal(:,:,:),
     *  eirene_mc_eipl_shi_bal(:,:,:),
     *  eirene_mc_eppl_shi_bal(:,:,:),
     *  b2stbm_shi0to3(:,:,:),
     *  b2stbm_shi_bal(:,:),
     *  ext_shi0to3(:,:,:),
     *  ext_shi_bal(:,:),
     *  b2stel_she_ion0to3(:,:,:),
     *  b2stel_she_ion_bal(:,:),
     *  b2stel_she_rec0to3(:,:,:),
     *  b2stel_she_rec_bal(:,:),
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
*     Time-averaged quantities
      real (kind=R8), allocatable, save ::
     *  fhi_32_av(:,:,:),
     *  fhi_52_av(:,:,:),
     *  fhi_cond_av(:,:,:),
     *  fhi_dia_av(:,:,:),
     *  fhi_ecrb_av(:,:,:),
     *  fhi_strange_av(:,:,:),
     *  fhi_pschused_av(:,:,:),
     *  fhi_inert_av(:,:,:),
     *  fhi_vispar_av(:,:,:),
     *  fhi_visper_av(:,:,:),
     *  fhi_visq_av(:,:,:),
     *  fhi_anml_av(:,:,:),
     *  fhi_kevis_av(:,:,:),
     *  fhi_vzh_av(:,:,:),
     *  b2stbr_phys_shi_bal_av(:,:),
     *  b2stbr_bas_shi_bal_av(:,:),
     *  b2stbr_first_flight_shi_bal_av(:,:),
     *  b2stbc_shi_bal_av(:,:),
     *  eirene_mc_eapl_shi_bal_av(:,:,:),
     *  eirene_mc_empl_shi_bal_av(:,:,:),
     *  eirene_mc_eipl_shi_bal_av(:,:,:),
     *  eirene_mc_eppl_shi_bal_av(:,:,:),
     *  b2stbm_shi_bal_av(:,:),
     *  ext_shi_bal_av(:,:),
     *  b2stel_she_ion_bal_av(:,:),
     *  b2stel_she_rec_bal_av(:,:),
     *  b2stel_shi_ion_bal_av(:,:),
     *  b2stel_shi_rec_bal_av(:,:),
     *  b2stcx_shi_bal_av(:,:),
     *  b2srsm_shi_bal_av(:,:),
     *  b2srdt_shi_bal_av(:,:),
     *  b2srst_shi_bal_av(:,:),
     *  b2sihs_diaa_bal_av(:,:),
     *  b2sihs_divua_bal_av(:,:),
     *  b2sihs_exba_bal_av(:,:),
     *  b2sihs_visa_bal_av(:,:),
     *  b2sihs_fraa_bal_av(:,:),
     *  reshi_av(:,:)
*     Additional averaged quantities
      real (kind=R8), allocatable, save ::
     *  ne_av(:,:),
     *  na_av(:,:,:),
     *  ua_av(:,:,:),
     *  po_av(:,:),
     *  te_av(:,:),
     *  ti_av(:,:),
     *  rpt_av(:,:,:),
     *  kinrgy_av(:,:,:),
     *  fne_av(:,:,:),
     *  dab2_av(:,:,:,:),
     *  dmb2_av(:,:,:,:),
     *  tab2_av(:,:,:,:),
     *  tmb2_av(:,:,:,:),
     *  rfluxa_av(:,:,:,:),
     *  rfluxm_av(:,:,:,:),
     *  refluxa_av(:,:,:,:),
     *  refluxm_av(:,:,:,:),
     *  pfluxa_av(:,:,:,:),
     *  pfluxm_av(:,:,:,:),
     *  pefluxa_av(:,:,:,:),
     *  pefluxm_av(:,:,:,:)

*   ..balance_netcdf decides whether to initialise, track, and write balance arrays.
*     Default: do not create balance arrays (balance_netcdf.eq.0)
      integer, save :: balance_netcdf

*   ..balance_average decides whether to average balance arrays over time.
*     Default: do not average (balance_average.eq.0)
      integer, save :: balance_average

      contains

      subroutine alloc_b2mod_balance(nxd,nyd,nsd)
      implicit none
*   ..Input arguments (unchanged on exit)
      integer nxd,nyd,nsd

*   ..Arrays for particle balance:
*     Fluxes:
      allocate(fna_pinch(-1:nxd,-1:nyd,0:1,0:nsd-1))
      allocate(fna_pll(-1:nxd,-1:nyd,0:1,0:nsd-1))
      allocate(fna_drift(-1:nxd,-1:nyd,0:1,0:nsd-1))
      allocate(fna_ch(-1:nxd,-1:nyd,0:1,0:nsd-1))
      allocate(fna_nanom(-1:nxd,-1:nyd,0:1,0:nsd-1))
      allocate(fna_panom(-1:nxd,-1:nyd,0:1,0:nsd-1))
      allocate(fna_pschused(-1:nxd,-1:nyd,0:1,0:nsd-1))
      allocate(fna_tot(-1:nxd,-1:nyd,0:1,0:nsd-1))
*     Sources
      allocate(b2stbr_phys_sna0to1(-1:nxd,-1:nyd,0:1,0:nsd-1))
      allocate(b2stbr_phys_sna_bal(-1:nxd,-1:nyd,0:nsd-1))
      allocate(b2stbr_bas_sna0to1(-1:nxd,-1:nyd,0:1,0:nsd-1))
      allocate(b2stbr_bas_sna_bal(-1:nxd,-1:nyd,0:nsd-1))
      allocate(b2stbr_first_flight_sna0to1(-1:nxd,-1:nyd,0:1,0:nsd-1))
      allocate(b2stbr_first_flight_sna_bal(-1:nxd,-1:nyd,0:nsd-1))
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
      allocate(tot_sna0to1(-1:nxd,-1:nyd,0:1,0:nsd-1))
      allocate(tot_sna_bal(-1:nxd,-1:nyd,0:nsd-1))
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
      tot_sna0to1=0.0_R8
      tot_sna_bal=0.0_R8
*     Time-averaged quantities
      allocate(fna_pinch_av(-1:nxd,-1:nyd,0:1,0:nsd-1))
      allocate(fna_pll_av(-1:nxd,-1:nyd,0:1,0:nsd-1))
      allocate(fna_drift_av(-1:nxd,-1:nyd,0:1,0:nsd-1))
      allocate(fna_ch_av(-1:nxd,-1:nyd,0:1,0:nsd-1))
      allocate(fna_nanom_av(-1:nxd,-1:nyd,0:1,0:nsd-1))
      allocate(fna_panom_av(-1:nxd,-1:nyd,0:1,0:nsd-1))
      allocate(fna_pschused_av(-1:nxd,-1:nyd,0:1,0:nsd-1))
      allocate(fna_tot_av(-1:nxd,-1:nyd,0:1,0:nsd-1))
      allocate(b2stbr_phys_sna_bal_av(-1:nxd,-1:nyd,0:nsd-1))
      allocate(b2stbr_bas_sna_bal_av(-1:nxd,-1:nyd,0:nsd-1))
      allocate(b2stbr_first_flight_sna_bal_av(-1:nxd,-1:nyd,0:nsd-1))
      allocate(b2stbc_sna_bal_av(-1:nxd,-1:nyd,0:nsd-1))
      allocate(b2stbm_sna_bal_av(-1:nxd,-1:nyd,0:nsd-1))
      allocate(ext_sna_bal_av(-1:nxd,-1:nyd,0:nsd-1))
      allocate(b2stel_sna_ion_bal_av(-1:nxd,-1:nyd,0:nsd-1))
      allocate(b2stel_sna_rec_bal_av(-1:nxd,-1:nyd,0:nsd-1))
      allocate(b2stcx_sna_bal_av(-1:nxd,-1:nyd,0:nsd-1))
      allocate(b2srsm_sna_bal_av(-1:nxd,-1:nyd,0:nsd-1))
      allocate(b2srdt_sna_bal_av(-1:nxd,-1:nyd,0:nsd-1))
      allocate(b2srst_sna_bal_av(-1:nxd,-1:nyd,0:nsd-1))
      allocate(tot_sna_bal_av(-1:nxd,-1:nyd,0:nsd-1))
      allocate(resco_av(-1:nxd,-1:nyd,0:nsd-1))
*   ..Arrays for momentum balance:
*     Fluxes
      allocate(fmo_flua(-1:nxd,-1:nyd,0:1,0:nsd-1))
      allocate(fmo_cvsa(-1:nxd,-1:nyd,0:1,0:nsd-1))
      allocate(fmo_hybr(-1:nxd,-1:nyd,0:1,0:nsd-1))
      allocate(fmo_b2nxfv(-1:nxd,-1:nyd,0:1,0:nsd-1))
      allocate(fmo_tot(-1:nxd,-1:nyd,0:1,0:nsd-1))
*     Sources
      allocate(b2stbr_phys_smo0to3(-1:nxd,-1:nyd,0:3,0:nsd-1))
      allocate(b2stbr_phys_smo_bal(-1:nxd,-1:nyd,0:nsd-1))
      allocate(b2stbr_bas_smo0to3(-1:nxd,-1:nyd,0:3,0:nsd-1))
      allocate(b2stbr_bas_smo_bal(-1:nxd,-1:nyd,0:nsd-1))
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
      allocate(b2sifr_smotf_ehxp0to3(-1:nxd,-1:nyd,0:3,0:nsd-1))
      allocate(b2sifr_smotf_ehxp_bal(-1:nxd,-1:nyd,0:nsd-1))
      allocate(b2sifr_smotf_cthe0to3(-1:nxd,-1:nyd,0:3,0:nsd-1))
      allocate(b2sifr_smotf_cthe_bal(-1:nxd,-1:nyd,0:nsd-1))
      allocate(b2sifr_smotf_cthi0to3(-1:nxd,-1:nyd,0:3,0:nsd-1))
      allocate(b2sifr_smotf_cthi_bal(-1:nxd,-1:nyd,0:nsd-1))
      allocate(b2sifr_smofrea0to3(-1:nxd,-1:nyd,0:3,0:nsd-1))
      allocate(b2sifr_smofrea_bal(-1:nxd,-1:nyd,0:nsd-1))
      allocate(b2sifr_smofria0to3(-1:nxd,-1:nyd,0:3,0:nsd-1))
      allocate(b2sifr_smofria_bal(-1:nxd,-1:nyd,0:nsd-1))
      allocate(b2sifr_smotfea0to3(-1:nxd,-1:nyd,0:3,0:nsd-1))
      allocate(b2sifr_smotfea_bal(-1:nxd,-1:nyd,0:nsd-1))
      allocate(b2sifr_smotfia0to3(-1:nxd,-1:nyd,0:3,0:nsd-1))
      allocate(b2sifr_smotfia_bal(-1:nxd,-1:nyd,0:nsd-1))
      allocate(b2siav_smovh0to3(-1:nxd,-1:nyd,0:3,0:nsd-1))
      allocate(b2siav_smovh_bal(-1:nxd,-1:nyd,0:nsd-1))
      allocate(b2siav_smovv0to3(-1:nxd,-1:nyd,0:3,0:nsd-1))
      allocate(b2siav_smovv_bal(-1:nxd,-1:nyd,0:nsd-1))
      allocate(b2siav_smovi0to3(-1:nxd,-1:nyd,0:3,0:nsd-1))            !som 13.08.21
      allocate(b2siav_smovi_bal(-1:nxd,-1:nyd,0:nsd-1))                !som 13.08.21
      allocate(b2sicf_smo0to3(-1:nxd,-1:nyd,0:3,0:nsd-1))
      allocate(b2sicf_smo_bal(-1:nxd,-1:nyd,0:nsd-1))
      allocate(b2sian_smo0to3(-1:nxd,-1:nyd,0:3,0:nsd-1))              !srv 25.09.17
      allocate(b2sian_smo_bal(-1:nxd,-1:nyd,0:nsd-1))                  !srv 25.09.17
      allocate(b2nxdv_smo0to3(-1:nxd,-1:nyd,0:3,0:nsd-1))
      allocate(b2nxdv_smo_bal(-1:nxd,-1:nyd,0:nsd-1))
      allocate(b2sigp_smogp0to3(-1:nxd,-1:nyd,0:3,0:nsd-1))
      allocate(b2sigp_smogp_bal(-1:nxd,-1:nyd,0:nsd-1))
      allocate(b2sigp_smogpi0to3(-1:nxd,-1:nyd,0:3,0:nsd-1))
      allocate(b2sigp_smogpi_bal(-1:nxd,-1:nyd,0:nsd-1))
      allocate(b2sigp_smogpe0to3(-1:nxd,-1:nyd,0:3,0:nsd-1))
      allocate(b2sigp_smogpe_bal(-1:nxd,-1:nyd,0:nsd-1))
      allocate(b2sigp_smogpgr0to3(-1:nxd,-1:nyd,0:3,0:nsd-1))
      allocate(b2sigp_smogpgr_bal(-1:nxd,-1:nyd,0:nsd-1))
      allocate(b2sigp_pstat0to3(-1:nxd,-1:nyd,0:3,0:nsd-1))
      allocate(b2sigp_pstat_bal(-1:nxd,-1:nyd,0:nsd-1))
      allocate(b2sigp_pstati0to3(-1:nxd,-1:nyd,0:3,0:nsd-1))
      allocate(b2sigp_pstati_bal(-1:nxd,-1:nyd,0:nsd-1))
      allocate(b2sigp_pstate0to3(-1:nxd,-1:nyd,0:3))
      allocate(b2sigp_pstate_bal(-1:nxd,-1:nyd))
      allocate(tot_smo0to3(-1:nxd,-1:nyd,0:3,0:nsd-1))
      allocate(tot_smo_bal(-1:nxd,-1:nyd,0:nsd-1))
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
      b2sifr_smotf_ehxp0to3=0.0_R8
      b2sifr_smotf_ehxp_bal=0.0_R8
      b2sifr_smotf_cthe0to3=0.0_R8
      b2sifr_smotf_cthe_bal=0.0_R8
      b2sifr_smotf_cthi0to3=0.0_R8
      b2sifr_smotf_cthi_bal=0.0_R8
      b2sifr_smofrea0to3=0.0_R8
      b2sifr_smofrea_bal=0.0_R8
      b2sifr_smofria0to3=0.0_R8
      b2sifr_smofria_bal=0.0_R8
      b2sifr_smotfea0to3=0.0_R8
      b2sifr_smotfea_bal=0.0_R8
      b2sifr_smotfia0to3=0.0_R8
      b2sifr_smotfia_bal=0.0_R8
      b2siav_smovh0to3=0.0_R8
      b2siav_smovh_bal=0.0_R8
      b2siav_smovv0to3=0.0_R8
      b2siav_smovv_bal=0.0_R8
      b2siav_smovi0to3=0.0_R8                  !som 13.08.21
      b2siav_smovi_bal=0.0_R8                  !som 13.08.21
      b2sicf_smo0to3=0.0_R8
      b2sicf_smo_bal=0.0_R8
      b2sian_smo0to3=0.0_R8                    !srv 25.09.17
      b2sian_smo_bal=0.0_R8                    !srv 25.09.17
      b2nxdv_smo0to3=0.0_R8
      b2nxdv_smo_bal=0.0_R8
      b2sigp_smogp0to3=0.0_R8
      b2sigp_smogp_bal=0.0_R8
      b2sigp_smogpi0to3=0.0_R8
      b2sigp_smogpi_bal=0.0_R8
      b2sigp_smogpe0to3=0.0_R8
      b2sigp_smogpe_bal=0.0_R8
      b2sigp_smogpgr0to3=0.0_R8
      b2sigp_smogpgr_bal=0.0_R8
      b2sigp_pstat0to3=0.0_R8
      b2sigp_pstat_bal=0.0_R8
      b2sigp_pstati0to3=0.0_R8
      b2sigp_pstati_bal=0.0_R8
      b2sigp_pstate0to3=0.0_R8
      b2sigp_pstate_bal=0.0_R8
      tot_smo0to3=0.0_R8
      tot_smo_bal=0.0_R8
*     Time-averaged quantities
      allocate(fmo_flua_av(-1:nxd,-1:nyd,0:1,0:nsd-1))
      allocate(fmo_cvsa_av(-1:nxd,-1:nyd,0:1,0:nsd-1))
      allocate(fmo_hybr_av(-1:nxd,-1:nyd,0:1,0:nsd-1))
      allocate(fmo_b2nxfv_av(-1:nxd,-1:nyd,0:1,0:nsd-1))
      allocate(fmo_tot_av(-1:nxd,-1:nyd,0:1,0:nsd-1))
      allocate(b2stbr_phys_smo_bal_av(-1:nxd,-1:nyd,0:nsd-1))
      allocate(b2stbr_bas_smo_bal_av(-1:nxd,-1:nyd,0:nsd-1))
      allocate(b2stbc_smo_bal_av(-1:nxd,-1:nyd,0:nsd-1))
      allocate(b2stbm_smo_bal_av(-1:nxd,-1:nyd,0:nsd-1))
      allocate(ext_smo_bal_av(-1:nxd,-1:nyd,0:nsd-1))
      allocate(b2stel_smq_ion_bal_av(-1:nxd,-1:nyd,0:nsd-1))
      allocate(b2stel_smq_rec_bal_av(-1:nxd,-1:nyd,0:nsd-1))
      allocate(b2stcx_smq_bal_av(-1:nxd,-1:nyd,0:nsd-1))
      allocate(b2srsm_smo_bal_av(-1:nxd,-1:nyd,0:nsd-1))
      allocate(b2srdt_smo_bal_av(-1:nxd,-1:nyd,0:nsd-1))
      allocate(b2srst_smo_bal_av(-1:nxd,-1:nyd,0:nsd-1))
      allocate(b2sifr_smoch_bal_av(-1:nxd,-1:nyd,0:nsd-1))
      allocate(b2sifr_smotf_ehxp_bal_av(-1:nxd,-1:nyd,0:nsd-1))
      allocate(b2sifr_smotf_cthe_bal_av(-1:nxd,-1:nyd,0:nsd-1))
      allocate(b2sifr_smotf_cthi_bal_av(-1:nxd,-1:nyd,0:nsd-1))
      allocate(b2sifr_smofrea_bal_av(-1:nxd,-1:nyd,0:nsd-1))
      allocate(b2sifr_smofria_bal_av(-1:nxd,-1:nyd,0:nsd-1))
      allocate(b2sifr_smotfea_bal_av(-1:nxd,-1:nyd,0:nsd-1))
      allocate(b2sifr_smotfia_bal_av(-1:nxd,-1:nyd,0:nsd-1))
      allocate(b2siav_smovh_bal_av(-1:nxd,-1:nyd,0:nsd-1))
      allocate(b2siav_smovv_bal_av(-1:nxd,-1:nyd,0:nsd-1))
      allocate(b2siav_smovi_bal_av(-1:nxd,-1:nyd,0:nsd-1))        !som 13.08.21
      allocate(b2sicf_smo_bal_av(-1:nxd,-1:nyd,0:nsd-1))
      allocate(b2sian_smo_bal_av(-1:nxd,-1:nyd,0:nsd-1))          !srv 25.09.17
      allocate(b2nxdv_smo_bal_av(-1:nxd,-1:nyd,0:nsd-1))
      allocate(b2sigp_smogp_bal_av(-1:nxd,-1:nyd,0:nsd-1))
      allocate(b2sigp_smogpi_bal_av(-1:nxd,-1:nyd,0:nsd-1))
      allocate(b2sigp_smogpe_bal_av(-1:nxd,-1:nyd,0:nsd-1))
      allocate(b2sigp_smogpgr_bal_av(-1:nxd,-1:nyd,0:nsd-1))
      allocate(b2sigp_pstat_bal_av(-1:nxd,-1:nyd,0:nsd-1))
      allocate(b2sigp_pstati_bal_av(-1:nxd,-1:nyd,0:nsd-1))
      allocate(b2sigp_pstate_bal_av(-1:nxd,-1:nyd))
      allocate(tot_smo_bal_av(-1:nxd,-1:nyd,0:nsd-1))
      allocate(resmo_av(-1:nxd,-1:nyd,0:nsd-1))
*   ..Arrays for electron heat balance:
*     Fluxes
      allocate(fhe_32(-1:nxd,-1:nyd,0:1))
      allocate(fhe_52(-1:nxd,-1:nyd,0:1))
      allocate(fhe_thermj(-1:nxd,-1:nyd,0:1))
      allocate(fhe_cond(-1:nxd,-1:nyd,0:1))
      allocate(fhe_dia(-1:nxd,-1:nyd,0:1))
      allocate(fhe_ecrb(-1:nxd,-1:nyd,0:1))
      allocate(fhe_strange(-1:nxd,-1:nyd,0:1))
      allocate(fhe_pschused(-1:nxd,-1:nyd,0:1))
*     Sources
      allocate(b2stbr_phys_she0to3(-1:nxd,-1:nyd,0:3))
      allocate(b2stbr_phys_she_bal(-1:nxd,-1:nyd))
      allocate(b2stbr_bas_she0to3(-1:nxd,-1:nyd,0:3))
      allocate(b2stbr_bas_she_bal(-1:nxd,-1:nyd))
      allocate(b2stbr_first_flight_she0to3(-1:nxd,-1:nyd,0:3))
      allocate(b2stbr_first_flight_she_bal(-1:nxd,-1:nyd))
      allocate(b2stbc_she0to3(-1:nxd,-1:nyd,0:3))
      allocate(b2stbc_she_bal(-1:nxd,-1:nyd))
      allocate(b2stbm_she0to3(-1:nxd,-1:nyd,0:3))
      allocate(b2stbm_she_bal(-1:nxd,-1:nyd))
      allocate(ext_she0to3(-1:nxd,-1:nyd,0:3))
      allocate(ext_she_bal(-1:nxd,-1:nyd))
      allocate(b2stel_she0to3(-1:nxd,-1:nyd,0:3,0:nsd-1))
      allocate(b2stel_she_bal(-1:nxd,-1:nyd,0:nsd-1))
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
      allocate(fhe_32_av(-1:nxd,-1:nyd,0:1))
      allocate(fhe_52_av(-1:nxd,-1:nyd,0:1))
      allocate(fhe_thermj_av(-1:nxd,-1:nyd,0:1))
      allocate(fhe_cond_av(-1:nxd,-1:nyd,0:1))
      allocate(fhe_dia_av(-1:nxd,-1:nyd,0:1))
      allocate(fhe_ecrb_av(-1:nxd,-1:nyd,0:1))
      allocate(fhe_strange_av(-1:nxd,-1:nyd,0:1))
      allocate(fhe_pschused_av(-1:nxd,-1:nyd,0:1))
      allocate(b2stbr_phys_she_bal_av(-1:nxd,-1:nyd))
      allocate(b2stbr_bas_she_bal_av(-1:nxd,-1:nyd))
      allocate(b2stbr_first_flight_she_bal_av(-1:nxd,-1:nyd))
      allocate(b2stbc_she_bal_av(-1:nxd,-1:nyd))
      allocate(b2stbm_she_bal_av(-1:nxd,-1:nyd))
      allocate(ext_she_bal_av(-1:nxd,-1:nyd))
      allocate(b2stel_she_bal_av(-1:nxd,-1:nyd,0:nsd-1))
      allocate(b2srsm_she_bal_av(-1:nxd,-1:nyd))
      allocate(b2srdt_she_bal_av(-1:nxd,-1:nyd))
      allocate(b2srst_she_bal_av(-1:nxd,-1:nyd))
      allocate(b2sihs_diae_bal_av(-1:nxd,-1:nyd))
      allocate(b2sihs_divue_bal_av(-1:nxd,-1:nyd))
      allocate(b2sihs_exbe_bal_av(-1:nxd,-1:nyd))
      allocate(b2sihs_joule_bal_av(-1:nxd,-1:nyd))
      allocate(b2npht_shei_bal_av(-1:nxd,-1:nyd))
      allocate(reshe_av(-1:nxd,-1:nyd))
*   ..Arrays for ion heat balance:
*     Fluxes
      allocate(fhi_32(-1:nxd,-1:nyd,0:1))
      allocate(fhi_52(-1:nxd,-1:nyd,0:1))
      allocate(fhi_cond(-1:nxd,-1:nyd,0:1))
      allocate(fhi_dia(-1:nxd,-1:nyd,0:1))
      allocate(fhi_ecrb(-1:nxd,-1:nyd,0:1))
      allocate(fhi_strange(-1:nxd,-1:nyd,0:1))
      allocate(fhi_pschused(-1:nxd,-1:nyd,0:1))
      allocate(fhi_inert(-1:nxd,-1:nyd,0:1))
      allocate(fhi_vispar(-1:nxd,-1:nyd,0:1))
      allocate(fhi_visper(-1:nxd,-1:nyd,0:1))
      allocate(fhi_visq(-1:nxd,-1:nyd,0:1))
      allocate(fhi_anml(-1:nxd,-1:nyd,0:1))
      allocate(fhi_kevis(-1:nxd,-1:nyd,0:1))
      allocate(fhi_vzh(-1:nxd,-1:nyd,0:1))
*     Sources
      allocate(b2stbr_phys_shi0to3(-1:nxd,-1:nyd,0:3))
      allocate(b2stbr_phys_shi_bal(-1:nxd,-1:nyd))
      allocate(b2stbr_bas_shi0to3(-1:nxd,-1:nyd,0:3))
      allocate(b2stbr_bas_shi_bal(-1:nxd,-1:nyd))
      allocate(b2stbr_first_flight_shi0to3(-1:nxd,-1:nyd,0:3))
      allocate(b2stbr_first_flight_shi_bal(-1:nxd,-1:nyd))
      allocate(b2stbc_shi0to3(-1:nxd,-1:nyd,0:3))
      allocate(b2stbc_shi_bal(-1:nxd,-1:nyd))
      allocate(b2stbm_shi0to3(-1:nxd,-1:nyd,0:3))
      allocate(b2stbm_shi_bal(-1:nxd,-1:nyd))
      allocate(ext_shi0to3(-1:nxd,-1:nyd,0:3))
      allocate(ext_shi_bal(-1:nxd,-1:nyd))
      allocate(b2stel_she_ion0to3(-1:nxd,-1:nyd,0:3))
      allocate(b2stel_she_ion_bal(-1:nxd,-1:nyd))
      allocate(b2stel_she_rec0to3(-1:nxd,-1:nyd,0:3))
      allocate(b2stel_she_rec_bal(-1:nxd,-1:nyd))
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
      b2stel_she_ion0to3=0.0_R8
      b2stel_she_ion_bal=0.0_R8
      b2stel_she_rec0to3=0.0_R8
      b2stel_she_rec_bal=0.0_R8
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
      allocate(fhi_32_av(-1:nxd,-1:nyd,0:1))
      allocate(fhi_52_av(-1:nxd,-1:nyd,0:1))
      allocate(fhi_cond_av(-1:nxd,-1:nyd,0:1))
      allocate(fhi_dia_av(-1:nxd,-1:nyd,0:1))
      allocate(fhi_ecrb_av(-1:nxd,-1:nyd,0:1))
      allocate(fhi_strange_av(-1:nxd,-1:nyd,0:1))
      allocate(fhi_pschused_av(-1:nxd,-1:nyd,0:1))
      allocate(fhi_inert_av(-1:nxd,-1:nyd,0:1))
      allocate(fhi_vispar_av(-1:nxd,-1:nyd,0:1))
      allocate(fhi_visper_av(-1:nxd,-1:nyd,0:1))
      allocate(fhi_visq_av(-1:nxd,-1:nyd,0:1))
      allocate(fhi_anml_av(-1:nxd,-1:nyd,0:1))
      allocate(fhi_kevis_av(-1:nxd,-1:nyd,0:1))
      allocate(fhi_vzh_av(-1:nxd,-1:nyd,0:1))
      allocate(b2stbr_phys_shi_bal_av(-1:nxd,-1:nyd))
      allocate(b2stbr_bas_shi_bal_av(-1:nxd,-1:nyd))
      allocate(b2stbr_first_flight_shi_bal_av(-1:nxd,-1:nyd))
      allocate(b2stbc_shi_bal_av(-1:nxd,-1:nyd))
      allocate(b2stbm_shi_bal_av(-1:nxd,-1:nyd))
      allocate(ext_shi_bal_av(-1:nxd,-1:nyd))
      allocate(b2stel_she_ion_bal_av(-1:nxd,-1:nyd))
      allocate(b2stel_she_rec_bal_av(-1:nxd,-1:nyd))
      allocate(b2stel_shi_ion_bal_av(-1:nxd,-1:nyd))
      allocate(b2stel_shi_rec_bal_av(-1:nxd,-1:nyd))
      allocate(b2stcx_shi_bal_av(-1:nxd,-1:nyd))
      allocate(b2srsm_shi_bal_av(-1:nxd,-1:nyd))
      allocate(b2srdt_shi_bal_av(-1:nxd,-1:nyd))
      allocate(b2srst_shi_bal_av(-1:nxd,-1:nyd))
      allocate(b2sihs_diaa_bal_av(-1:nxd,-1:nyd))
      allocate(b2sihs_divua_bal_av(-1:nxd,-1:nyd))
      allocate(b2sihs_exba_bal_av(-1:nxd,-1:nyd))
      allocate(b2sihs_visa_bal_av(-1:nxd,-1:nyd))
      allocate(b2sihs_fraa_bal_av(-1:nxd,-1:nyd))
      allocate(reshi_av(-1:nxd,-1:nyd))
*     Additional averaged quantities
      allocate(ne_av(-1:nxd,-1:nyd))
      allocate(na_av(-1:nxd,-1:nyd,0:nsd-1))
      allocate(ua_av(-1:nxd,-1:nyd,0:nsd-1))
      allocate(po_av(-1:nxd,-1:nyd))
      allocate(te_av(-1:nxd,-1:nyd))
      allocate(ti_av(-1:nxd,-1:nyd))
      allocate(rpt_av(-1:nxd,-1:nyd,0:nsd-1))
      allocate(kinrgy_av(-1:nxd,-1:nyd,0:nsd-1))
      allocate(fne_av(-1:nxd,-1:nyd,0:1))

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
      deallocate(fna_tot)
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
      deallocate(tot_sna0to1)
      deallocate(tot_sna_bal)
*     Time-averaged quantities
      deallocate(fna_pinch_av)
      deallocate(fna_pll_av)
      deallocate(fna_drift_av)
      deallocate(fna_ch_av)
      deallocate(fna_nanom_av)
      deallocate(fna_panom_av)
      deallocate(fna_pschused_av)
      deallocate(fna_tot_av)
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
      deallocate(tot_sna_bal_av)
      deallocate(resco_av)
*   ..Arrays for momentum balance:
*     Fluxes
      deallocate(fmo_flua)
      deallocate(fmo_cvsa)
      deallocate(fmo_hybr)
      deallocate(fmo_b2nxfv)
      deallocate(fmo_tot)
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
      deallocate(b2sifr_smotf_ehxp0to3)
      deallocate(b2sifr_smotf_ehxp_bal)
      deallocate(b2sifr_smotf_cthe0to3)
      deallocate(b2sifr_smotf_cthe_bal)
      deallocate(b2sifr_smotf_cthi0to3)
      deallocate(b2sifr_smotf_cthi_bal)
      deallocate(b2sifr_smofrea0to3)
      deallocate(b2sifr_smofrea_bal)
      deallocate(b2sifr_smofria0to3)
      deallocate(b2sifr_smofria_bal)
      deallocate(b2sifr_smotfea0to3)
      deallocate(b2sifr_smotfea_bal)
      deallocate(b2sifr_smotfia0to3)
      deallocate(b2sifr_smotfia_bal)
      deallocate(b2siav_smovh0to3)
      deallocate(b2siav_smovh_bal)
      deallocate(b2siav_smovv0to3)
      deallocate(b2siav_smovv_bal)
      deallocate(b2siav_smovi0to3)            !som 13.08.21
      deallocate(b2siav_smovi_bal)            !som 13.08.21
      deallocate(b2sicf_smo0to3)
      deallocate(b2sicf_smo_bal)
      deallocate(b2sian_smo0to3)              !srv 25.09.17
      deallocate(b2sian_smo_bal)              !srv 25.09.17
      deallocate(b2nxdv_smo0to3)
      deallocate(b2nxdv_smo_bal)
      deallocate(b2sigp_smogp0to3)
      deallocate(b2sigp_smogp_bal)
      deallocate(b2sigp_smogpi0to3)
      deallocate(b2sigp_smogpi_bal)
      deallocate(b2sigp_smogpe0to3)
      deallocate(b2sigp_smogpe_bal)
      deallocate(b2sigp_smogpgr0to3)
      deallocate(b2sigp_smogpgr_bal)
      deallocate(b2sigp_pstat0to3)
      deallocate(b2sigp_pstat_bal)
      deallocate(b2sigp_pstati0to3)
      deallocate(b2sigp_pstati_bal)
      deallocate(b2sigp_pstate0to3)
      deallocate(b2sigp_pstate_bal)
      deallocate(tot_smo0to3)
      deallocate(tot_smo_bal)
*     Time-averaged quantities
      deallocate(fmo_flua_av)
      deallocate(fmo_cvsa_av)
      deallocate(fmo_hybr_av)
      deallocate(fmo_b2nxfv_av)
      deallocate(fmo_tot_av)
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
      deallocate(b2sifr_smotf_ehxp_bal_av)
      deallocate(b2sifr_smotf_cthe_bal_av)
      deallocate(b2sifr_smotf_cthi_bal_av)
      deallocate(b2sifr_smofrea_bal_av)
      deallocate(b2sifr_smofria_bal_av)
      deallocate(b2sifr_smotfea_bal_av)
      deallocate(b2sifr_smotfia_bal_av)
      deallocate(b2siav_smovh_bal_av)
      deallocate(b2siav_smovv_bal_av)
      deallocate(b2siav_smovi_bal_av)            !som 13.08.21
      deallocate(b2sicf_smo_bal_av)
      deallocate(b2sian_smo_bal_av)              !srv 25.09.17
      deallocate(b2nxdv_smo_bal_av)
      deallocate(b2sigp_smogp_bal_av)
      deallocate(b2sigp_smogpi_bal_av)
      deallocate(b2sigp_smogpe_bal_av)
      deallocate(b2sigp_smogpgr_bal_av)
      deallocate(b2sigp_pstat_bal_av)
      deallocate(b2sigp_pstati_bal_av)
      deallocate(b2sigp_pstate_bal_av)
      deallocate(tot_smo_bal_av)
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
      deallocate(fhi_vzh)
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
      deallocate(b2stel_she_ion0to3)
      deallocate(b2stel_she_ion_bal)
      deallocate(b2stel_she_rec0to3)
      deallocate(b2stel_she_rec_bal)
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
      deallocate(fhi_vzh_av)
      deallocate(b2stbr_phys_shi_bal_av)
      deallocate(b2stbr_bas_shi_bal_av)
      deallocate(b2stbr_first_flight_shi_bal_av)
      deallocate(b2stbc_shi_bal_av)
      deallocate(b2stbm_shi_bal_av)
      deallocate(ext_shi_bal_av)
      deallocate(b2stel_she_ion_bal_av)
      deallocate(b2stel_she_rec_bal_av)
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

      subroutine alloc_b2mod_eirene_sources
     &          (nxd,nyd,nsd,natmd,nmold,niond,nstrad)
#ifdef B25_EIRENE
#ifndef SOLPS4_3
      use eirmod_wneutrals ! IGNORE
     & , only : dab2,dmb2,tab2,tmb2,rfluxa,rfluxm,refluxa,refluxm,
     &          pfluxa,pfluxm,pefluxa,pefluxm
#else
      use eirdiag          ! IGNORE
     & , only : dab2,dmb2,tab2,tmb2,rfluxa,rfluxm,refluxa,refluxm,
     &          pfluxa,pfluxm,pefluxa,pefluxm
#endif
#else
      use b2mod_eirdiag
#endif
      implicit none
      integer nxd,nyd,nsd,natmd,nmold,niond,nstrad
*   ..procedures
      intrinsic size

cdjm Jan2017 Eirene arrays for particle balance:
      allocate(eirene_mc_pael_sne0to1(-1:nxd,-1:nyd,0:1,1:nstrad))
      allocate(
     & eirene_mc_papl_sna0to1(-1:nxd,-1:nyd,0:1,0:nsd-1,1:nstrad))
      allocate(eirene_mc_pmel_sne0to1(-1:nxd,-1:nyd,0:1,1:nstrad))
      allocate(
     & eirene_mc_pmpl_sna0to1(-1:nxd,-1:nyd,0:1,0:nsd-1,1:nstrad))
      allocate(
     & eirene_mc_pipl_sna0to1(-1:nxd,-1:nyd,0:1,0:nsd-1,1:nstrad))
      allocate(
     & eirene_mc_pppl_sna0to1(-1:nxd,-1:nyd,0:1,0:nsd-1,1:nstrad))
      allocate(
     & eirene_mc_core_sna0to1(-1:nxd,-1:nyd,0:1,0:nsd-1))
      allocate(eirene_mc_pael_sne_bal(-1:nxd,-1:nyd,1:nstrad))
      allocate(eirene_mc_papl_sna_bal(-1:nxd,-1:nyd,0:nsd-1,1:nstrad))
      allocate(eirene_mc_pmel_sne_bal(-1:nxd,-1:nyd,1:nstrad))
      allocate(eirene_mc_pmpl_sna_bal(-1:nxd,-1:nyd,0:nsd-1,1:nstrad))
      allocate(eirene_mc_pipl_sna_bal(-1:nxd,-1:nyd,0:nsd-1,1:nstrad))
      allocate(eirene_mc_pppl_sna_bal(-1:nxd,-1:nyd,0:nsd-1,1:nstrad))
      if (natmd.gt.0) then
        allocate(eirene_mc_paat_sna_bal(-1:nxd,-1:nyd,natmd,1:nstrad))
        allocate(eirene_mc_pmat_sna_bal(-1:nxd,-1:nyd,natmd,1:nstrad))
        allocate(eirene_mc_piat_sna_bal(-1:nxd,-1:nyd,natmd,1:nstrad))
      else
        allocate(eirene_mc_paat_sna_bal(-1:nxd,-1:nyd,1,1:nstrad))
        allocate(eirene_mc_pmat_sna_bal(-1:nxd,-1:nyd,1,1:nstrad))
        allocate(eirene_mc_piat_sna_bal(-1:nxd,-1:nyd,1,1:nstrad))
      end if
      if (nmold.gt.0) then
        allocate(eirene_mc_paml_sna_bal(-1:nxd,-1:nyd,nmold,1:nstrad))
        allocate(eirene_mc_pmml_sna_bal(-1:nxd,-1:nyd,nmold,1:nstrad))
        allocate(eirene_mc_piml_sna_bal(-1:nxd,-1:nyd,nmold,1:nstrad))
      else
        allocate(eirene_mc_paml_sna_bal(-1:nxd,-1:nyd,1,1:nstrad))
        allocate(eirene_mc_pmml_sna_bal(-1:nxd,-1:nyd,1,1:nstrad))
        allocate(eirene_mc_piml_sna_bal(-1:nxd,-1:nyd,1,1:nstrad))
      end if
      if (niond.gt.0) then
        allocate(eirene_mc_paio_sna_bal(-1:nxd,-1:nyd,niond,1:nstrad))
        allocate(eirene_mc_pmio_sna_bal(-1:nxd,-1:nyd,niond,1:nstrad))
        allocate(eirene_mc_piio_sna_bal(-1:nxd,-1:nyd,niond,1:nstrad))
      else
        allocate(eirene_mc_paio_sna_bal(-1:nxd,-1:nyd,1,1:nstrad))
        allocate(eirene_mc_pmio_sna_bal(-1:nxd,-1:nyd,1,1:nstrad))
        allocate(eirene_mc_piio_sna_bal(-1:nxd,-1:nyd,1,1:nstrad))
      end if
      allocate(eirene_mc_core_sna_bal(-1:nxd,-1:nyd,0:nsd-1))
      eirene_mc_pael_sne0to1=0.0_R8
      eirene_mc_papl_sna0to1=0.0_R8
      eirene_mc_pmel_sne0to1=0.0_R8
      eirene_mc_pmpl_sna0to1=0.0_R8
      eirene_mc_pipl_sna0to1=0.0_R8
      eirene_mc_pppl_sna0to1=0.0_R8
      eirene_mc_core_sna0to1=0.0_R8
      eirene_mc_pael_sne_bal=0.0_R8
      eirene_mc_papl_sna_bal=0.0_R8
      eirene_mc_pmel_sne_bal=0.0_R8
      eirene_mc_pmpl_sna_bal=0.0_R8
      eirene_mc_pipl_sna_bal=0.0_R8
      eirene_mc_pppl_sna_bal=0.0_R8
      eirene_mc_paat_sna_bal=0.0_R8
      eirene_mc_pmat_sna_bal=0.0_R8
      eirene_mc_piat_sna_bal=0.0_R8
      eirene_mc_paml_sna_bal=0.0_R8
      eirene_mc_pmml_sna_bal=0.0_R8
      eirene_mc_piml_sna_bal=0.0_R8
      eirene_mc_paio_sna_bal=0.0_R8
      eirene_mc_pmio_sna_bal=0.0_R8
      eirene_mc_piio_sna_bal=0.0_R8
      eirene_mc_core_sna_bal=0.0_R8
cdjm Jun2017 Time-averaged arrays for particle balance:
      allocate(eirene_mc_pael_sne_bal_av(-1:nxd,-1:nyd,1:nstrad))
      allocate(
     & eirene_mc_papl_sna_bal_av(-1:nxd,-1:nyd,0:nsd-1,1:nstrad))
      allocate(eirene_mc_pmel_sne_bal_av(-1:nxd,-1:nyd,1:nstrad))
      allocate(
     & eirene_mc_pmpl_sna_bal_av(-1:nxd,-1:nyd,0:nsd-1,1:nstrad))
      allocate(
     & eirene_mc_pipl_sna_bal_av(-1:nxd,-1:nyd,0:nsd-1,1:nstrad))
      allocate(
     & eirene_mc_pppl_sna_bal_av(-1:nxd,-1:nyd,0:nsd-1,1:nstrad))
      if (natmd.gt.0) then
        allocate(
     &   eirene_mc_paat_sna_bal_av(-1:nxd,-1:nyd,natmd,1:nstrad))
        allocate(
     &   eirene_mc_pmat_sna_bal_av(-1:nxd,-1:nyd,natmd,1:nstrad))
        allocate(
     &   eirene_mc_piat_sna_bal_av(-1:nxd,-1:nyd,natmd,1:nstrad))
      else
        allocate(eirene_mc_paat_sna_bal_av(-1:nxd,-1:nyd,1,1:nstrad))
        allocate(eirene_mc_pmat_sna_bal_av(-1:nxd,-1:nyd,1,1:nstrad))
        allocate(eirene_mc_piat_sna_bal_av(-1:nxd,-1:nyd,1,1:nstrad))
      end if
      if (nmold.gt.0) then
        allocate(
     &   eirene_mc_paml_sna_bal_av(-1:nxd,-1:nyd,nmold,1:nstrad))
        allocate(
     &   eirene_mc_pmml_sna_bal_av(-1:nxd,-1:nyd,nmold,1:nstrad))
        allocate(
     &   eirene_mc_piml_sna_bal_av(-1:nxd,-1:nyd,nmold,1:nstrad))
      else
        allocate(eirene_mc_paml_sna_bal_av(-1:nxd,-1:nyd,1,1:nstrad))
        allocate(eirene_mc_pmml_sna_bal_av(-1:nxd,-1:nyd,1,1:nstrad))
        allocate(eirene_mc_piml_sna_bal_av(-1:nxd,-1:nyd,1,1:nstrad))
      end if
      if (niond.gt.0) then
        allocate(
     &   eirene_mc_paio_sna_bal_av(-1:nxd,-1:nyd,niond,1:nstrad))
        allocate(
     &   eirene_mc_pmio_sna_bal_av(-1:nxd,-1:nyd,niond,1:nstrad))
        allocate(
     &   eirene_mc_piio_sna_bal_av(-1:nxd,-1:nyd,niond,1:nstrad))
      else
        allocate(eirene_mc_paio_sna_bal_av(-1:nxd,-1:nyd,1,1:nstrad))
        allocate(eirene_mc_pmio_sna_bal_av(-1:nxd,-1:nyd,1,1:nstrad))
        allocate(eirene_mc_piio_sna_bal_av(-1:nxd,-1:nyd,1,1:nstrad))
      end if
      allocate(
     & eirene_mc_core_sna_bal_av(-1:nxd,-1:nyd,0:nsd-1))
cdjm Jan2017 Eirene arrays for momentum balance:
      allocate(
     & eirene_mc_mapl_smo0to3(-1:nxd,-1:nyd,0:3,0:nsd-1,1:nstrad))
      allocate(
     & eirene_mc_mmpl_smo0to3(-1:nxd,-1:nyd,0:3,0:nsd-1,1:nstrad))
      allocate(
     & eirene_mc_mipl_smo0to3(-1:nxd,-1:nyd,0:3,0:nsd-1,1:nstrad))
      allocate(
     & eirene_mc_mppl_smo0to3(-1:nxd,-1:nyd,0:3,0:nsd-1,1:nstrad))
      allocate(eirene_mc_mapl_smo_bal(-1:nxd,-1:nyd,0:nsd-1,1:nstrad))
      allocate(eirene_mc_mmpl_smo_bal(-1:nxd,-1:nyd,0:nsd-1,1:nstrad))
      allocate(eirene_mc_mipl_smo_bal(-1:nxd,-1:nyd,0:nsd-1,1:nstrad))
      allocate(eirene_mc_mppl_smo_bal(-1:nxd,-1:nyd,0:nsd-1,1:nstrad))
      eirene_mc_mapl_smo0to3=0.0_R8
      eirene_mc_mmpl_smo0to3=0.0_R8
      eirene_mc_mipl_smo0to3=0.0_R8
      eirene_mc_mppl_smo0to3=0.0_R8
      eirene_mc_mapl_smo_bal=0.0_R8
      eirene_mc_mmpl_smo_bal=0.0_R8
      eirene_mc_mipl_smo_bal=0.0_R8
      eirene_mc_mppl_smo_bal=0.0_R8
cdjm Jun2017 Time-averaged arrays for momentum balance:
      allocate(
     & eirene_mc_mapl_smo_bal_av(-1:nxd,-1:nyd,0:nsd-1,1:nstrad))
      allocate(
     & eirene_mc_mmpl_smo_bal_av(-1:nxd,-1:nyd,0:nsd-1,1:nstrad))
      allocate(
     & eirene_mc_mipl_smo_bal_av(-1:nxd,-1:nyd,0:nsd-1,1:nstrad))
      allocate(
     & eirene_mc_mppl_smo_bal_av(-1:nxd,-1:nyd,0:nsd-1,1:nstrad))
cdjm Jan2017 Eirene arrays for electron heat balance:
      allocate(eirene_mc_eael_she0to3(-1:nxd,-1:nyd,0:3,1:nstrad))
      allocate(eirene_mc_emel_she0to3(-1:nxd,-1:nyd,0:3,1:nstrad))
      allocate(eirene_mc_eiel_she0to3(-1:nxd,-1:nyd,0:3,1:nstrad))
      allocate(eirene_mc_epel_she0to3(-1:nxd,-1:nyd,0:3,1:nstrad))
      allocate(eirene_mc_eael_she_bal(-1:nxd,-1:nyd,1:nstrad))
      allocate(eirene_mc_emel_she_bal(-1:nxd,-1:nyd,1:nstrad))
      allocate(eirene_mc_eiel_she_bal(-1:nxd,-1:nyd,1:nstrad))
      allocate(eirene_mc_epel_she_bal(-1:nxd,-1:nyd,1:nstrad))
      eirene_mc_eael_she0to3=0.0_R8
      eirene_mc_emel_she0to3=0.0_R8
      eirene_mc_eiel_she0to3=0.0_R8
      eirene_mc_epel_she0to3=0.0_R8
      eirene_mc_eael_she_bal=0.0_R8
      eirene_mc_emel_she_bal=0.0_R8
      eirene_mc_eiel_she_bal=0.0_R8
      eirene_mc_epel_she_bal=0.0_R8      
cdjm Jun2017 Time-averaged arrays for electron heat balance:
      allocate(eirene_mc_eael_she_bal_av(-1:nxd,-1:nyd,1:nstrad))
      allocate(eirene_mc_emel_she_bal_av(-1:nxd,-1:nyd,1:nstrad))
      allocate(eirene_mc_eiel_she_bal_av(-1:nxd,-1:nyd,1:nstrad))
      allocate(eirene_mc_epel_she_bal_av(-1:nxd,-1:nyd,1:nstrad))
cdjm Jan2017 Eirene arrays for ion heat balance:
      allocate(eirene_mc_eapl_shi0to3(-1:nxd,-1:nyd,0:3,1:nstrad))
      allocate(eirene_mc_empl_shi0to3(-1:nxd,-1:nyd,0:3,1:nstrad))
      allocate(eirene_mc_eipl_shi0to3(-1:nxd,-1:nyd,0:3,1:nstrad))
      allocate(eirene_mc_eppl_shi0to3(-1:nxd,-1:nyd,0:3,1:nstrad))
      allocate(eirene_mc_eapl_shi_bal(-1:nxd,-1:nyd,1:nstrad))
      allocate(eirene_mc_empl_shi_bal(-1:nxd,-1:nyd,1:nstrad))
      allocate(eirene_mc_eipl_shi_bal(-1:nxd,-1:nyd,1:nstrad))
      allocate(eirene_mc_eppl_shi_bal(-1:nxd,-1:nyd,1:nstrad))
      eirene_mc_eapl_shi0to3=0.0_R8
      eirene_mc_empl_shi0to3=0.0_R8
      eirene_mc_eipl_shi0to3=0.0_R8
      eirene_mc_eppl_shi0to3=0.0_R8
      eirene_mc_eapl_shi_bal=0.0_R8
      eirene_mc_empl_shi_bal=0.0_R8
      eirene_mc_eipl_shi_bal=0.0_R8
      eirene_mc_eppl_shi_bal=0.0_R8
cdjm Jun2017 Time-averaged arrays for ion heat balance:
      allocate(eirene_mc_eapl_shi_bal_av(-1:nxd,-1:nyd,1:nstrad))
      allocate(eirene_mc_empl_shi_bal_av(-1:nxd,-1:nyd,1:nstrad))
      allocate(eirene_mc_eipl_shi_bal_av(-1:nxd,-1:nyd,1:nstrad))
      allocate(eirene_mc_eppl_shi_bal_av(-1:nxd,-1:nyd,1:nstrad))
cdjm Jun2017 Time-averaged arrays for fluxes and densities:
      allocate(dab2_av(0:(size(dab2,1)-1),0:(size(dab2,2)-1),
     &                   size(dab2,3),1))
      allocate(dmb2_av(0:(size(dmb2,1)-1),0:(size(dmb2,2)-1),
     &                   size(dmb2,3),1))
      allocate(tab2_av(0:(size(tab2,1)-1),0:(size(tab2,2)-1),
     &                   size(tab2,3),1))
      allocate(tmb2_av(0:(size(tmb2,1)-1),0:(size(tmb2,2)-1),
     &                   size(tmb2,3),1))
      allocate(rfluxa_av(0:(size(rfluxa,1)-1),0:(size(rfluxa,2)-1),
     &                   size(rfluxa,3),1))
      allocate(rfluxm_av(0:(size(rfluxm,1)-1),0:(size(rfluxm,2)-1),
     &                   size(rfluxm,3),1))
      allocate(refluxa_av(0:(size(refluxa,1)-1),0:(size(refluxa,2)-1),
     &                   size(refluxa,3),1))
      allocate(refluxm_av(0:(size(refluxm,1)-1),0:(size(refluxm,2)-1),
     &                   size(refluxm,3),1))
      allocate(pfluxa_av(0:(size(pfluxa,1)-1),0:(size(pfluxa,2)-1),
     &                   size(pfluxa,3),1))
      allocate(pfluxm_av(0:(size(pfluxm,1)-1),0:(size(pfluxm,2)-1),
     &                   size(pfluxm,3),1))
      allocate(pefluxa_av(0:(size(pefluxa,1)-1),0:(size(pefluxa,2)-1),
     &                   size(pefluxa,3),1))
      allocate(pefluxm_av(0:(size(pefluxm,1)-1),0:(size(pefluxm,2)-1),
     &                   size(pefluxm,3),1))

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
      deallocate(eirene_mc_pael_sne0to1)
      deallocate(eirene_mc_papl_sna0to1)
      deallocate(eirene_mc_pmel_sne0to1)
      deallocate(eirene_mc_pmpl_sna0to1)
      deallocate(eirene_mc_pipl_sna0to1)
      deallocate(eirene_mc_pppl_sna0to1)
      deallocate(eirene_mc_core_sna0to1)
      deallocate(eirene_mc_pael_sne_bal)
      deallocate(eirene_mc_papl_sna_bal)
      deallocate(eirene_mc_pmel_sne_bal)
      deallocate(eirene_mc_pmpl_sna_bal)
      deallocate(eirene_mc_pipl_sna_bal)
      deallocate(eirene_mc_pppl_sna_bal)
      deallocate(eirene_mc_paat_sna_bal)
      deallocate(eirene_mc_pmat_sna_bal)
      deallocate(eirene_mc_piat_sna_bal)
      deallocate(eirene_mc_paml_sna_bal)
      deallocate(eirene_mc_pmml_sna_bal)
      deallocate(eirene_mc_piml_sna_bal)
      deallocate(eirene_mc_paio_sna_bal)
      deallocate(eirene_mc_pmio_sna_bal)
      deallocate(eirene_mc_piio_sna_bal)
      deallocate(eirene_mc_core_sna_bal)
      deallocate(eirene_mc_pael_sne_bal_av)
      deallocate(eirene_mc_papl_sna_bal_av)
      deallocate(eirene_mc_pmel_sne_bal_av)
      deallocate(eirene_mc_pmpl_sna_bal_av)
      deallocate(eirene_mc_pipl_sna_bal_av)
      deallocate(eirene_mc_pppl_sna_bal_av)
      deallocate(eirene_mc_paat_sna_bal_av)
      deallocate(eirene_mc_pmat_sna_bal_av)
      deallocate(eirene_mc_piat_sna_bal_av)
      deallocate(eirene_mc_paml_sna_bal_av)
      deallocate(eirene_mc_pmml_sna_bal_av)
      deallocate(eirene_mc_piml_sna_bal_av)
      deallocate(eirene_mc_paio_sna_bal_av)
      deallocate(eirene_mc_pmio_sna_bal_av)
      deallocate(eirene_mc_piio_sna_bal_av)
      deallocate(eirene_mc_core_sna_bal_av)
cdjm Jan2017 Eirene arrays for momentum balance:
      deallocate(eirene_mc_mapl_smo0to3)
      deallocate(eirene_mc_mmpl_smo0to3)
      deallocate(eirene_mc_mipl_smo0to3)
      deallocate(eirene_mc_mppl_smo0to3)
      deallocate(eirene_mc_mapl_smo_bal)
      deallocate(eirene_mc_mmpl_smo_bal)
      deallocate(eirene_mc_mipl_smo_bal)
      deallocate(eirene_mc_mppl_smo_bal)
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

      subroutine update_balance_sources_particles(nx,ny,isb,nb) !djm Jan2017
      implicit none
*   ..input arguments (unchanged on exit)
      integer nx,ny,isb
      real (kind=R8):: nb(-1:nx,-1:ny)
*   ..local variables
      integer ix,iy,istrai,nstrad
*   ..procedures
      intrinsic size

      nstrad = size(eirene_mc_papl_sna_bal,4)
      
      do ix=-1,nx
        do iy=-1,ny
          b2stbr_phys_sna_bal(ix,iy,isb)=
     &                b2stbr_phys_sna0to1(ix,iy,0,isb)+
     &                b2stbr_phys_sna0to1(ix,iy,1,isb)*nb(ix,iy)
          b2stbr_bas_sna_bal(ix,iy,isb)=
     &                b2stbr_bas_sna0to1(ix,iy,0,isb)+
     &                b2stbr_bas_sna0to1(ix,iy,1,isb)*nb(ix,iy)
          b2stbr_first_flight_sna_bal(ix,iy,isb)=
     &                b2stbr_first_flight_sna0to1(ix,iy,0,isb)+
     &                b2stbr_first_flight_sna0to1(ix,iy,1,isb)*nb(ix,iy)
          b2stbc_sna_bal(ix,iy,isb)=b2stbc_sna0to1(ix,iy,0,isb)+
     &                     b2stbc_sna0to1(ix,iy,1,isb)*nb(ix,iy)
          if (allocated(eirene_mc_papl_sna_bal)) then
            do istrai=1,nstrad
              eirene_mc_papl_sna_bal(ix,iy,isb,istrai)=
     &         eirene_mc_papl_sna0to1(ix,iy,0,isb,istrai)+
     &         eirene_mc_papl_sna0to1(ix,iy,1,isb,istrai)*nb(ix,iy)
              eirene_mc_pmpl_sna_bal(ix,iy,isb,istrai)=
     &         eirene_mc_pmpl_sna0to1(ix,iy,0,isb,istrai)+
     &         eirene_mc_pmpl_sna0to1(ix,iy,1,isb,istrai)*nb(ix,iy)
              eirene_mc_pipl_sna_bal(ix,iy,isb,istrai)=
     &         eirene_mc_pipl_sna0to1(ix,iy,0,isb,istrai)+
     &         eirene_mc_pipl_sna0to1(ix,iy,1,isb,istrai)*nb(ix,iy)
              eirene_mc_pppl_sna_bal(ix,iy,isb,istrai)=
     &         eirene_mc_pppl_sna0to1(ix,iy,0,isb,istrai)+
     &         eirene_mc_pppl_sna0to1(ix,iy,1,isb,istrai)*nb(ix,iy)
            enddo
            eirene_mc_core_sna_bal(ix,iy,isb)=
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
          tot_sna_bal(ix,iy,isb)=tot_sna0to1(ix,iy,0,isb)+
     &                     tot_sna0to1(ix,iy,1,isb)*nb(ix,iy)
        enddo
      enddo
      return
      end subroutine update_balance_sources_particles

      subroutine update_balance_sources_momentum(nx,ny,isb,ub,rob) !djm Jan2017
      implicit none
*   ..input arguments (unchanged on exit)
      integer nx,ny,isb
      real (kind=R8):: ub(-1:nx,-1:ny), rob(-1:nx,-1:ny)
*   ..local variables
      integer ix,iy,istrai,nstrad
*   ..procedures
      intrinsic size

      nstrad = size(eirene_mc_mapl_smo_bal,4)

      do ix=-1,nx
        do iy=-1,ny
          b2stbr_phys_smo_bal(ix,iy,isb)=
     &            b2stbr_phys_smo0to3(ix,iy,0,isb)+
     &            b2stbr_phys_smo0to3(ix,iy,1,isb)*ub(ix,iy)+
     &            b2stbr_phys_smo0to3(ix,iy,2,isb)*rob(ix,iy)+
     &            b2stbr_phys_smo0to3(ix,iy,3,isb)*rob(ix,iy)*ub(ix,iy)
          b2stbr_bas_smo_bal(ix,iy,isb)=
     &            b2stbr_bas_smo0to3(ix,iy,0,isb)+
     &            b2stbr_bas_smo0to3(ix,iy,1,isb)*ub(ix,iy)+
     &            b2stbr_bas_smo0to3(ix,iy,2,isb)*rob(ix,iy)+
     &            b2stbr_bas_smo0to3(ix,iy,3,isb)*rob(ix,iy)*ub(ix,iy)
          b2stbc_smo_bal(ix,iy,isb)=b2stbc_smo0to3(ix,iy,0,isb)+
     &            b2stbc_smo0to3(ix,iy,1,isb)*ub(ix,iy)+
     &            b2stbc_smo0to3(ix,iy,2,isb)*rob(ix,iy)+
     &            b2stbc_smo0to3(ix,iy,3,isb)*rob(ix,iy)*ub(ix,iy)
          if (allocated(eirene_mc_mapl_smo_bal)) then
            do istrai=1,nstrad
              eirene_mc_mapl_smo_bal(ix,iy,isb,istrai)=
     &         eirene_mc_mapl_smo0to3(ix,iy,0,isb,istrai)+
     &         eirene_mc_mapl_smo0to3(ix,iy,1,isb,istrai)*ub(ix,iy)+
     &         eirene_mc_mapl_smo0to3(ix,iy,2,isb,istrai)*rob(ix,iy)+
     &         eirene_mc_mapl_smo0to3(ix,iy,3,isb,istrai)*rob(ix,iy)*
     &                                                    ub(ix,iy)
              eirene_mc_mmpl_smo_bal(ix,iy,isb,istrai)=
     &         eirene_mc_mmpl_smo0to3(ix,iy,0,isb,istrai)+
     &         eirene_mc_mmpl_smo0to3(ix,iy,1,isb,istrai)*ub(ix,iy)+
     &         eirene_mc_mmpl_smo0to3(ix,iy,2,isb,istrai)*rob(ix,iy)+
     &         eirene_mc_mmpl_smo0to3(ix,iy,3,isb,istrai)*rob(ix,iy)*
     &                                                    ub(ix,iy)
              eirene_mc_mipl_smo_bal(ix,iy,isb,istrai)=
     &         eirene_mc_mipl_smo0to3(ix,iy,0,isb,istrai)+
     &         eirene_mc_mipl_smo0to3(ix,iy,1,isb,istrai)*ub(ix,iy)+
     &         eirene_mc_mipl_smo0to3(ix,iy,2,isb,istrai)*rob(ix,iy)+
     &         eirene_mc_mipl_smo0to3(ix,iy,3,isb,istrai)*rob(ix,iy)*
     &                                                    ub(ix,iy)
              eirene_mc_mppl_smo_bal(ix,iy,isb,istrai)=
     &         eirene_mc_mppl_smo0to3(ix,iy,0,isb,istrai)+
     &         eirene_mc_mppl_smo0to3(ix,iy,1,isb,istrai)*ub(ix,iy)+
     &         eirene_mc_mppl_smo0to3(ix,iy,2,isb,istrai)*rob(ix,iy)+
     &         eirene_mc_mppl_smo0to3(ix,iy,3,isb,istrai)*rob(ix,iy)*
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
          b2sifr_smotf_ehxp_bal(ix,iy,isb)=
     &           b2sifr_smotf_ehxp0to3(ix,iy,0,isb)+
     &           b2sifr_smotf_ehxp0to3(ix,iy,1,isb)*ub(ix,iy)+
     &           b2sifr_smotf_ehxp0to3(ix,iy,2,isb)*rob(ix,iy)+
     &           b2sifr_smotf_ehxp0to3(ix,iy,3,isb)*rob(ix,iy)*ub(ix,iy)
          b2sifr_smotf_cthe_bal(ix,iy,isb)=
     &           b2sifr_smotf_cthe0to3(ix,iy,0,isb)+
     &           b2sifr_smotf_cthe0to3(ix,iy,1,isb)*ub(ix,iy)+
     &           b2sifr_smotf_cthe0to3(ix,iy,2,isb)*rob(ix,iy)+
     &           b2sifr_smotf_cthe0to3(ix,iy,3,isb)*rob(ix,iy)*ub(ix,iy)
          b2sifr_smotf_cthi_bal(ix,iy,isb)=
     &           b2sifr_smotf_cthi0to3(ix,iy,0,isb)+
     &           b2sifr_smotf_cthi0to3(ix,iy,1,isb)*ub(ix,iy)+
     &           b2sifr_smotf_cthi0to3(ix,iy,2,isb)*rob(ix,iy)+
     &           b2sifr_smotf_cthi0to3(ix,iy,3,isb)*rob(ix,iy)*ub(ix,iy)
          b2sifr_smofrea_bal(ix,iy,isb)=
     &           b2sifr_smofrea0to3(ix,iy,0,isb)+
     &           b2sifr_smofrea0to3(ix,iy,1,isb)*ub(ix,iy)+
     &           b2sifr_smofrea0to3(ix,iy,2,isb)*rob(ix,iy)+
     &           b2sifr_smofrea0to3(ix,iy,3,isb)*rob(ix,iy)*ub(ix,iy)
          b2sifr_smofria_bal(ix,iy,isb)=
     &           b2sifr_smofria0to3(ix,iy,0,isb)+
     &           b2sifr_smofria0to3(ix,iy,1,isb)*ub(ix,iy)+
     &           b2sifr_smofria0to3(ix,iy,2,isb)*rob(ix,iy)+
     &           b2sifr_smofria0to3(ix,iy,3,isb)*rob(ix,iy)*ub(ix,iy)
          b2sifr_smotfea_bal(ix,iy,isb)=
     &           b2sifr_smotfea0to3(ix,iy,0,isb)+
     &           b2sifr_smotfea0to3(ix,iy,1,isb)*ub(ix,iy)+
     &           b2sifr_smotfea0to3(ix,iy,2,isb)*rob(ix,iy)+
     &           b2sifr_smotfea0to3(ix,iy,3,isb)*rob(ix,iy)*ub(ix,iy)
          b2sifr_smotfia_bal(ix,iy,isb)=
     &           b2sifr_smotfia0to3(ix,iy,0,isb)+
     &           b2sifr_smotfia0to3(ix,iy,1,isb)*ub(ix,iy)+
     &           b2sifr_smotfia0to3(ix,iy,2,isb)*rob(ix,iy)+
     &           b2sifr_smotfia0to3(ix,iy,3,isb)*rob(ix,iy)*ub(ix,iy)
          b2siav_smovh_bal(ix,iy,isb)=b2siav_smovh0to3(ix,iy,0,isb)+
     &            b2siav_smovh0to3(ix,iy,1,isb)*ub(ix,iy)+
     &            b2siav_smovh0to3(ix,iy,2,isb)*rob(ix,iy)+
     &            b2siav_smovh0to3(ix,iy,3,isb)*rob(ix,iy)*ub(ix,iy)
          b2siav_smovv_bal(ix,iy,isb)=b2siav_smovv0to3(ix,iy,0,isb)+
     &            b2siav_smovv0to3(ix,iy,1,isb)*ub(ix,iy)+
     &            b2siav_smovv0to3(ix,iy,2,isb)*rob(ix,iy)+
     &            b2siav_smovv0to3(ix,iy,3,isb)*rob(ix,iy)*ub(ix,iy)
          b2siav_smovi_bal(ix,iy,isb)=b2siav_smovi0to3(ix,iy,0,isb)+
     &            b2siav_smovi0to3(ix,iy,1,isb)*ub(ix,iy)+
     &            b2siav_smovi0to3(ix,iy,2,isb)*rob(ix,iy)+
     &            b2siav_smovi0to3(ix,iy,3,isb)*rob(ix,iy)*ub(ix,iy)    !som 13.08.21
          b2sicf_smo_bal(ix,iy,isb)=b2sicf_smo0to3(ix,iy,0,isb)+
     &            b2sicf_smo0to3(ix,iy,1,isb)*ub(ix,iy)+
     &            b2sicf_smo0to3(ix,iy,2,isb)*rob(ix,iy)+
     &            b2sicf_smo0to3(ix,iy,3,isb)*rob(ix,iy)*ub(ix,iy)
          b2sian_smo_bal(ix,iy,isb)=b2sian_smo0to3(ix,iy,0,isb)+
     &            b2sian_smo0to3(ix,iy,1,isb)*ub(ix,iy)+
     &            b2sian_smo0to3(ix,iy,2,isb)*rob(ix,iy)+
     &            b2sian_smo0to3(ix,iy,3,isb)*rob(ix,iy)*ub(ix,iy)      !srv 25.09.17
          b2nxdv_smo_bal(ix,iy,isb)=b2nxdv_smo0to3(ix,iy,0,isb)+
     &            b2nxdv_smo0to3(ix,iy,1,isb)*ub(ix,iy)+
     &            b2nxdv_smo0to3(ix,iy,2,isb)*rob(ix,iy)+
     &            b2nxdv_smo0to3(ix,iy,3,isb)*rob(ix,iy)*ub(ix,iy)
          b2sigp_smogp_bal(ix,iy,isb)=b2sigp_smogp0to3(ix,iy,0,isb)+
     &            b2sigp_smogp0to3(ix,iy,1,isb)*ub(ix,iy)+
     &            b2sigp_smogp0to3(ix,iy,2,isb)*rob(ix,iy)+
     &            b2sigp_smogp0to3(ix,iy,3,isb)*rob(ix,iy)*ub(ix,iy)
          b2sigp_smogpi_bal(ix,iy,isb)=b2sigp_smogpi0to3(ix,iy,0,isb)+
     &            b2sigp_smogpi0to3(ix,iy,1,isb)*ub(ix,iy)+
     &            b2sigp_smogpi0to3(ix,iy,2,isb)*rob(ix,iy)+
     &            b2sigp_smogpi0to3(ix,iy,3,isb)*rob(ix,iy)*ub(ix,iy)
          b2sigp_smogpe_bal(ix,iy,isb)=b2sigp_smogpe0to3(ix,iy,0,isb)+
     &            b2sigp_smogpe0to3(ix,iy,1,isb)*ub(ix,iy)+
     &            b2sigp_smogpe0to3(ix,iy,2,isb)*rob(ix,iy)+
     &            b2sigp_smogpe0to3(ix,iy,3,isb)*rob(ix,iy)*ub(ix,iy)
          b2sigp_smogpgr_bal(ix,iy,isb)=b2sigp_smogpgr0to3(ix,iy,0,isb)+
     &            b2sigp_smogpgr0to3(ix,iy,1,isb)*ub(ix,iy)+
     &            b2sigp_smogpgr0to3(ix,iy,2,isb)*rob(ix,iy)+
     &            b2sigp_smogpgr0to3(ix,iy,3,isb)*rob(ix,iy)*ub(ix,iy)
          b2sigp_pstat_bal(ix,iy,isb)=b2sigp_pstat0to3(ix,iy,0,isb)+
     &            b2sigp_pstat0to3(ix,iy,1,isb)*ub(ix,iy)+
     &            b2sigp_pstat0to3(ix,iy,2,isb)*rob(ix,iy)+
     &            b2sigp_pstat0to3(ix,iy,3,isb)*rob(ix,iy)*ub(ix,iy)
          b2sigp_pstati_bal(ix,iy,isb)=b2sigp_pstati0to3(ix,iy,0,isb)+
     &            b2sigp_pstati0to3(ix,iy,1,isb)*ub(ix,iy)+
     &            b2sigp_pstati0to3(ix,iy,2,isb)*rob(ix,iy)+
     &            b2sigp_pstati0to3(ix,iy,3,isb)*rob(ix,iy)*ub(ix,iy)
          b2sigp_pstate_bal(ix,iy)=b2sigp_pstate0to3(ix,iy,0)+
     &            b2sigp_pstate0to3(ix,iy,1)*ub(ix,iy)+
     &            b2sigp_pstate0to3(ix,iy,2)*rob(ix,iy)+
     &            b2sigp_pstate0to3(ix,iy,3)*rob(ix,iy)*ub(ix,iy)
          tot_smo_bal(ix,iy,isb)=tot_smo0to3(ix,iy,0,isb)+
     &            tot_smo0to3(ix,iy,1,isb)*ub(ix,iy)+
     &            tot_smo0to3(ix,iy,2,isb)*rob(ix,iy)+
     &            tot_smo0to3(ix,iy,3,isb)*rob(ix,iy)*ub(ix,iy)
        enddo
      enddo
      return
      end subroutine update_balance_sources_momentum

      subroutine update_balance_sources_heat(nx,ny,ne,ni,te,ti,ns)
      implicit none
*   ..input arguments (unchanged on exit)
      integer nx,ny,ns
      real (kind=R8):: ne(-1:nx,-1:ny), ni(-1:nx,-1:ny),
     &                 te(-1:nx,-1:ny), ti(-1:nx,-1:ny)
*   ..local variables
      integer ix,iy,istrai,nstrad,is
*   ..procedures
      intrinsic size

      nstrad = size(eirene_mc_eael_she_bal,3)

cdjm Update electron heat sources
cxpb We also treat the electron particle sources
      do ix=-1,nx
        do iy=-1,ny
          b2stbr_phys_she_bal(ix,iy)=
     &     b2stbr_phys_she0to3(ix,iy,0)+
     &     b2stbr_phys_she0to3(ix,iy,1)*te(ix,iy)+
     &     b2stbr_phys_she0to3(ix,iy,2)*ne(ix,iy)+
     &     b2stbr_phys_she0to3(ix,iy,3)*ne(ix,iy)*te(ix,iy)
          b2stbr_bas_she_bal(ix,iy)=
     &     b2stbr_bas_she0to3(ix,iy,0)+
     &     b2stbr_bas_she0to3(ix,iy,1)*te(ix,iy)+
     &     b2stbr_bas_she0to3(ix,iy,2)*ne(ix,iy)+
     &     b2stbr_bas_she0to3(ix,iy,3)*ne(ix,iy)*te(ix,iy)
          b2stbr_first_flight_she_bal(ix,iy)=
     &     b2stbr_first_flight_she0to3(ix,iy,0)+
     &     b2stbr_first_flight_she0to3(ix,iy,1)*te(ix,iy)+
     &     b2stbr_first_flight_she0to3(ix,iy,2)*ne(ix,iy)+
     &     b2stbr_first_flight_she0to3(ix,iy,3)*ne(ix,iy)*te(ix,iy)
          b2stbc_she_bal(ix,iy)=
     &     b2stbc_she0to3(ix,iy,0)+
     &     b2stbc_she0to3(ix,iy,1)*te(ix,iy)+
     &     b2stbc_she0to3(ix,iy,2)*ne(ix,iy)+
     &     b2stbc_she0to3(ix,iy,3)*ne(ix,iy)*te(ix,iy)
          if (allocated(eirene_mc_eael_she_bal)) then
            do istrai=1,nstrad
              eirene_mc_pael_sne_bal(ix,iy,istrai)=
     &         eirene_mc_pael_sne0to1(ix,iy,0,istrai)+
     &         eirene_mc_pael_sne0to1(ix,iy,1,istrai)*ne(ix,iy)
              eirene_mc_pmel_sne_bal(ix,iy,istrai)=
     &         eirene_mc_pmel_sne0to1(ix,iy,0,istrai)+
     &         eirene_mc_pmel_sne0to1(ix,iy,1,istrai)*ne(ix,iy)
              eirene_mc_eael_she_bal(ix,iy,istrai)=
     &         eirene_mc_eael_she0to3(ix,iy,0,istrai)+
     &         eirene_mc_eael_she0to3(ix,iy,1,istrai)*te(ix,iy)+
     &         eirene_mc_eael_she0to3(ix,iy,2,istrai)*ne(ix,iy)+
     &         eirene_mc_eael_she0to3(ix,iy,3,istrai)*ne(ix,iy)*
     &                                                te(ix,iy)
              eirene_mc_emel_she_bal(ix,iy,istrai)=
     &         eirene_mc_emel_she0to3(ix,iy,0,istrai)+
     &         eirene_mc_emel_she0to3(ix,iy,1,istrai)*te(ix,iy)+
     &         eirene_mc_emel_she0to3(ix,iy,2,istrai)*ne(ix,iy)+
     &         eirene_mc_emel_she0to3(ix,iy,3,istrai)*ne(ix,iy)*
     &                                                te(ix,iy)
              eirene_mc_eiel_she_bal(ix,iy,istrai)=
     &         eirene_mc_eiel_she0to3(ix,iy,0,istrai)+
     &         eirene_mc_eiel_she0to3(ix,iy,1,istrai)*te(ix,iy)+
     &         eirene_mc_eiel_she0to3(ix,iy,2,istrai)*ne(ix,iy)+
     &         eirene_mc_eiel_she0to3(ix,iy,3,istrai)*ne(ix,iy)*
     &                                                te(ix,iy)
              eirene_mc_epel_she_bal(ix,iy,istrai)=
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
          do is = 0, ns-1
            b2stel_she_bal(ix,iy,is)=
     &       b2stel_she0to3(ix,iy,0,is)+
     &       b2stel_she0to3(ix,iy,1,is)*te(ix,iy)+
     &       b2stel_she0to3(ix,iy,2,is)*ne(ix,iy)+
     &       b2stel_she0to3(ix,iy,3,is)*ne(ix,iy)*te(ix,iy)
          enddo
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
          b2stbr_phys_shi_bal(ix,iy)=
     &     b2stbr_phys_shi0to3(ix,iy,0)+
     &     b2stbr_phys_shi0to3(ix,iy,1)*ti(ix,iy)+
     &     b2stbr_phys_shi0to3(ix,iy,2)*ni(ix,iy)+
     &     b2stbr_phys_shi0to3(ix,iy,3)*ni(ix,iy)*ti(ix,iy)
          b2stbr_bas_shi_bal(ix,iy)=
     &     b2stbr_bas_shi0to3(ix,iy,0)+
     &     b2stbr_bas_shi0to3(ix,iy,1)*ti(ix,iy)+
     &     b2stbr_bas_shi0to3(ix,iy,2)*ni(ix,iy)+
     &     b2stbr_bas_shi0to3(ix,iy,3)*ni(ix,iy)*ti(ix,iy)
          b2stbr_first_flight_shi_bal(ix,iy)=
     &     b2stbr_first_flight_shi0to3(ix,iy,0)+
     &     b2stbr_first_flight_shi0to3(ix,iy,1)*ti(ix,iy)+
     &     b2stbr_first_flight_shi0to3(ix,iy,2)*ni(ix,iy)+
     &     b2stbr_first_flight_shi0to3(ix,iy,3)*ni(ix,iy)*ti(ix,iy)
          b2stbc_shi_bal(ix,iy)=
     &     b2stbc_shi0to3(ix,iy,0)+
     &     b2stbc_shi0to3(ix,iy,1)*ti(ix,iy)+
     &     b2stbc_shi0to3(ix,iy,2)*ni(ix,iy)+
     &     b2stbc_shi0to3(ix,iy,3)*ni(ix,iy)*ti(ix,iy)
          if (allocated(eirene_mc_eapl_shi_bal)) then
            do istrai=1,nstrad
              eirene_mc_eapl_shi_bal(ix,iy,istrai)=
     &         eirene_mc_eapl_shi0to3(ix,iy,0,istrai)+
     &         eirene_mc_eapl_shi0to3(ix,iy,1,istrai)*ti(ix,iy)+
     &         eirene_mc_eapl_shi0to3(ix,iy,2,istrai)*ni(ix,iy)+
     &         eirene_mc_eapl_shi0to3(ix,iy,3,istrai)*ni(ix,iy)*
     &                                                ti(ix,iy)
              eirene_mc_empl_shi_bal(ix,iy,istrai)=
     &         eirene_mc_empl_shi0to3(ix,iy,0,istrai)+
     &         eirene_mc_empl_shi0to3(ix,iy,1,istrai)*ti(ix,iy)+
     &         eirene_mc_empl_shi0to3(ix,iy,2,istrai)*ni(ix,iy)+
     &         eirene_mc_empl_shi0to3(ix,iy,3,istrai)*ni(ix,iy)*
     &                                                ti(ix,iy)
              eirene_mc_eipl_shi_bal(ix,iy,istrai)=
     &         eirene_mc_eipl_shi0to3(ix,iy,0,istrai)+
     &         eirene_mc_eipl_shi0to3(ix,iy,1,istrai)*ti(ix,iy)+
     &         eirene_mc_eipl_shi0to3(ix,iy,2,istrai)*ni(ix,iy)+
     &         eirene_mc_eipl_shi0to3(ix,iy,3,istrai)*ni(ix,iy)*
     &                                                ti(ix,iy)
              eirene_mc_eppl_shi_bal(ix,iy,istrai)=
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
          b2stel_she_ion_bal(ix,iy)=
     &     b2stel_she_ion0to3(ix,iy,0)+
     &     b2stel_she_ion0to3(ix,iy,1)*te(ix,iy)+
     &     b2stel_she_ion0to3(ix,iy,2)*ne(ix,iy)+
     &     b2stel_she_ion0to3(ix,iy,3)*ne(ix,iy)*te(ix,iy)
          b2stel_she_rec_bal(ix,iy)=
     &     b2stel_she_rec0to3(ix,iy,0)+
     &     b2stel_she_rec0to3(ix,iy,1)*te(ix,iy)+
     &     b2stel_she_rec0to3(ix,iy,2)*ne(ix,iy)+
     &     b2stel_she_rec0to3(ix,iy,3)*ne(ix,iy)*te(ix,iy)
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

      subroutine update_average_balance(itim,ntim)
      use b2mod_residuals
     & , only : resco,resmo,reshe,reshi
      use b2mod_plasma
     & , only : ne,na,ua,po,te,ti,kinrgy,fne
      use b2mod_rates
     & , only : rpt
#ifdef B25_EIRENE
#ifndef SOLPS4_3
      use eirmod_wneutrals ! IGNORE
     & , only : dab2,dmb2,tab2,tmb2,rfluxa,rfluxm,refluxa,refluxm,
     &          pfluxa,pfluxm,pefluxa,pefluxm
#else
      use eirdiag          ! IGNORE
     & , only : dab2,dmb2,tab2,tmb2,rfluxa,rfluxm,refluxa,refluxm,
     &          pfluxa,pfluxm,pefluxa,pefluxm
#endif
#else
      use b2mod_eirdiag
#endif
      implicit none
*   ..input arguments (unchanged on exit)
      integer, intent(in) :: itim,ntim
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
        fna_tot_av=0.0_R8
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
        tot_sna_bal_av=0.0_R8
        resco_av=0.0_R8
        if (allocated(eirene_mc_papl_sna_bal)) then
          eirene_mc_pael_sne_bal_av=0.0_R8
          eirene_mc_papl_sna_bal_av=0.0_R8
          eirene_mc_pmel_sne_bal_av=0.0_R8
          eirene_mc_pmpl_sna_bal_av=0.0_R8
          eirene_mc_pipl_sna_bal_av=0.0_R8
          eirene_mc_pppl_sna_bal_av=0.0_R8
          eirene_mc_paat_sna_bal_av=0.0_R8
          eirene_mc_pmat_sna_bal_av=0.0_R8
          eirene_mc_piat_sna_bal_av=0.0_R8
          eirene_mc_paml_sna_bal_av=0.0_R8
          eirene_mc_pmml_sna_bal_av=0.0_R8
          eirene_mc_piml_sna_bal_av=0.0_R8
          eirene_mc_paio_sna_bal_av=0.0_R8
          eirene_mc_pmio_sna_bal_av=0.0_R8
          eirene_mc_piio_sna_bal_av=0.0_R8
          eirene_mc_core_sna_bal_av=0.0_R8
        endif
*     ..Arrays for momentum balance:
        fmo_flua_av=0.0_R8
        fmo_cvsa_av=0.0_R8
        fmo_hybr_av=0.0_R8
        fmo_b2nxfv_av=0.0_R8
        fmo_tot_av=0.0_R8
        b2stbr_phys_smo_bal_av=0.0_R8
        b2stbr_bas_smo_bal_av=0.0_R8
        b2stbc_smo_bal_av=0.0_R8
        if (allocated(eirene_mc_mapl_smo_bal)) then
          eirene_mc_mapl_smo_bal_av=0.0_R8
          eirene_mc_mmpl_smo_bal_av=0.0_R8
          eirene_mc_mipl_smo_bal_av=0.0_R8
          eirene_mc_mppl_smo_bal_av=0.0_R8
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
        b2sifr_smotf_ehxp_bal_av=0.0_R8
        b2sifr_smotf_cthe_bal_av=0.0_R8
        b2sifr_smotf_cthi_bal_av=0.0_R8
        b2sifr_smofrea_bal_av=0.0_R8
        b2sifr_smofria_bal_av=0.0_R8
        b2sifr_smotfea_bal_av=0.0_R8
        b2sifr_smotfia_bal_av=0.0_R8
        b2siav_smovh_bal_av=0.0_R8
        b2siav_smovv_bal_av=0.0_R8
        b2siav_smovi_bal_av=0.0_R8         !som 13.08.21
        b2sicf_smo_bal_av=0.0_R8
        b2sian_smo_bal_av=0.0_R8           !srv 25.09.17
        b2nxdv_smo_bal_av=0.0_R8
        b2sigp_smogp_bal_av=0.0_R8
        b2sigp_smogpi_bal_av=0.0_R8
        b2sigp_smogpe_bal_av=0.0_R8
        b2sigp_smogpgr_bal_av=0.0_R8
        b2sigp_pstati_bal_av=0.0_R8
        b2sigp_pstate_bal_av=0.0_R8
        b2sigp_pstat_bal_av=0.0_R8
        tot_smo_bal_av=0.0_R8
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
        fhi_vzh_av=0.0_R8
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
        b2stel_she_ion_bal_av=0.0_R8
        b2stel_she_rec_bal_av=0.0_R8
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

*     Average now accumulates total to be rescaled before write
*     using ntim < 0
      if (balance_average.ne.0) then
        if (ntim < 0) then
           avmlt = rratio(1,itim)
           curmlt = 0._R8
        else
           avmlt = 1.0_R8
           curmlt = 1.0_R8
        endif
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
      fna_tot_av = avmlt*fna_tot_av +
     &                  curmlt*fna_tot
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
        eirene_mc_pael_sne_bal_av = avmlt*eirene_mc_pael_sne_bal_av +
     &                              curmlt*eirene_mc_pael_sne_bal
        eirene_mc_papl_sna_bal_av = avmlt*eirene_mc_papl_sna_bal_av +
     &                              curmlt*eirene_mc_papl_sna_bal
        eirene_mc_pmel_sne_bal_av = avmlt*eirene_mc_pmel_sne_bal_av +
     &                              curmlt*eirene_mc_pmel_sne_bal
        eirene_mc_pmpl_sna_bal_av = avmlt*eirene_mc_pmpl_sna_bal_av +
     &                              curmlt*eirene_mc_pmpl_sna_bal
        eirene_mc_pipl_sna_bal_av = avmlt*eirene_mc_pipl_sna_bal_av +
     &                              curmlt*eirene_mc_pipl_sna_bal
        eirene_mc_pppl_sna_bal_av = avmlt*eirene_mc_pppl_sna_bal_av +
     &                              curmlt*eirene_mc_pppl_sna_bal
        eirene_mc_core_sna_bal_av = avmlt*eirene_mc_core_sna_bal_av +
     &                              curmlt*eirene_mc_core_sna_bal
        eirene_mc_paat_sna_bal_av = avmlt*eirene_mc_paat_sna_bal_av +
     &                              curmlt*eirene_mc_paat_sna_bal
        eirene_mc_pmat_sna_bal_av = avmlt*eirene_mc_pmat_sna_bal_av +
     &                              curmlt*eirene_mc_pmat_sna_bal
        eirene_mc_piat_sna_bal_av = avmlt*eirene_mc_piat_sna_bal_av +
     &                              curmlt*eirene_mc_piat_sna_bal
        eirene_mc_paml_sna_bal_av = avmlt*eirene_mc_paml_sna_bal_av +
     &                              curmlt*eirene_mc_paml_sna_bal
        eirene_mc_pmml_sna_bal_av = avmlt*eirene_mc_pmml_sna_bal_av +
     &                              curmlt*eirene_mc_pmml_sna_bal
        eirene_mc_piml_sna_bal_av = avmlt*eirene_mc_piml_sna_bal_av +
     &                              curmlt*eirene_mc_piml_sna_bal
        eirene_mc_paio_sna_bal_av = avmlt*eirene_mc_paio_sna_bal_av +
     &                              curmlt*eirene_mc_paio_sna_bal
        eirene_mc_pmio_sna_bal_av = avmlt*eirene_mc_pmio_sna_bal_av +
     &                              curmlt*eirene_mc_pmio_sna_bal
        eirene_mc_piio_sna_bal_av = avmlt*eirene_mc_piio_sna_bal_av +
     &                              curmlt*eirene_mc_piio_sna_bal
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
      tot_sna_bal_av = avmlt*tot_sna_bal_av +
     &                    curmlt*tot_sna_bal
      resco_av = avmlt*resco_av +
     &           curmlt*resco
*   ..Arrays for momentum balance:
      fmo_flua_av = avmlt*fmo_flua_av +
     &              curmlt*fmo_flua
      fmo_cvsa_av = avmlt*fmo_cvsa_av +
     &              curmlt*fmo_cvsa
      fmo_hybr_av = avmlt*fmo_hybr_av +
     &              curmlt*fmo_hybr
      fmo_b2nxfv_av = avmlt*fmo_b2nxfv_av +
     &                curmlt*fmo_b2nxfv
      fmo_tot_av = avmlt*fmo_tot_av +
     &                curmlt*fmo_tot
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
        eirene_mc_mppl_smo_bal_av = avmlt*eirene_mc_mppl_smo_bal_av +
     &                              curmlt*eirene_mc_mppl_smo_bal
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
      b2sifr_smotf_ehxp_bal_av = avmlt*b2sifr_smotf_ehxp_bal_av +
     &                      curmlt*b2sifr_smotf_ehxp_bal
      b2sifr_smotf_cthe_bal_av = avmlt*b2sifr_smotf_cthe_bal_av +
     &                      curmlt*b2sifr_smotf_cthe_bal
      b2sifr_smotf_cthi_bal_av = avmlt*b2sifr_smotf_cthi_bal_av +
     &                      curmlt*b2sifr_smotf_cthi_bal
      b2sifr_smofrea_bal_av = avmlt*b2sifr_smofrea_bal_av +
     &                      curmlt*b2sifr_smofrea_bal
      b2sifr_smofria_bal_av = avmlt*b2sifr_smofria_bal_av +
     &                      curmlt*b2sifr_smofria_bal
      b2sifr_smotfea_bal_av = avmlt*b2sifr_smotfea_bal_av +
     &                      curmlt*b2sifr_smotfea_bal
      b2sifr_smotfia_bal_av = avmlt*b2sifr_smotfia_bal_av +
     &                      curmlt*b2sifr_smotfia_bal
      b2siav_smovh_bal_av = avmlt*b2siav_smovh_bal_av +
     &                      curmlt*b2siav_smovh_bal
      b2siav_smovv_bal_av = avmlt*b2siav_smovv_bal_av +
     &                      curmlt*b2siav_smovv_bal
      b2siav_smovi_bal_av = avmlt*b2siav_smovi_bal_av +
     &                      curmlt*b2siav_smovi_bal       !som 13.08.21
      b2sicf_smo_bal_av = avmlt*b2sicf_smo_bal_av +
     &                    curmlt*b2sicf_smo_bal
      b2sian_smo_bal_av = avmlt*b2sian_smo_bal_av +
     &                    curmlt*b2sian_smo_bal           !srv 25.09.17
      b2nxdv_smo_bal_av = avmlt*b2nxdv_smo_bal_av +
     &                    curmlt*b2nxdv_smo_bal
      b2sigp_smogp_bal_av = avmlt*b2sigp_smogp_bal_av +
     &                      curmlt*b2sigp_smogp_bal
      b2sigp_smogpi_bal_av = avmlt*b2sigp_smogpi_bal_av +
     &                       curmlt*b2sigp_smogpi_bal
      b2sigp_smogpe_bal_av = avmlt*b2sigp_smogpe_bal_av +
     &                       curmlt*b2sigp_smogpe_bal
      b2sigp_smogpgr_bal_av = avmlt*b2sigp_smogpgr_bal_av +
     &                       curmlt*b2sigp_smogpgr_bal
      b2sigp_pstat_bal_av = avmlt*b2sigp_pstat_bal_av +
     &                      curmlt*b2sigp_pstat_bal
      b2sigp_pstati_bal_av = avmlt*b2sigp_pstati_bal_av +
     &                       curmlt*b2sigp_pstati_bal
      b2sigp_pstate_bal_av = avmlt*b2sigp_pstate_bal_av +
     &                       curmlt*b2sigp_pstate_bal
      tot_smo_bal_av = avmlt*tot_smo_bal_av +
     &                 curmlt*tot_smo_bal
      resmo_av = avmlt*resmo_av +
     &           curmlt*resmo
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
     &           curmlt*reshe
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
     &              curmlt*fhi_visq
      fhi_anml_av = avmlt*fhi_anml_av +
     &              curmlt*fhi_anml
      fhi_kevis_av = avmlt*fhi_kevis_av +
     &               curmlt*fhi_kevis
      fhi_vzh_av = avmlt*fhi_vzh_av +
     &             curmlt*fhi_vzh
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
      b2stel_she_ion_bal_av = avmlt*b2stel_she_ion_bal_av +
     &                        curmlt*b2stel_she_ion_bal
      b2stel_she_rec_bal_av = avmlt*b2stel_she_rec_bal_av +
     &                        curmlt*b2stel_she_rec_bal
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
     &           curmlt*reshi
*   ..Additional averaged quantities
      ne_av=avmlt*ne_av+curmlt*ne
      na_av=avmlt*na_av+curmlt*na
      ua_av=avmlt*ua_av+curmlt*ua
      po_av=avmlt*po_av+curmlt*po
      te_av=avmlt*te_av+curmlt*te
      ti_av=avmlt*ti_av+curmlt*ti
      rpt_av=avmlt*rpt_av+curmlt*rpt
      kinrgy_av=avmlt*kinrgy_av+curmlt*kinrgy
      fne_av=avmlt*fne_av+curmlt*fne
      if (allocated(dab2_av)) then
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
      endif

      end subroutine update_average_balance

      subroutine write_balance(nxd,nyd,nsd)
#ifndef NO_CDF
      use b2mod_geo
     & , only : crx,cry,hx,hy,hz,bb,vol,gs
      use b2mod_indirect
     & , only : leftix,leftiy,rightix,rightiy,
     &          topix,topiy,bottomix,bottomiy
      use b2mod_b2cmpa
     & , only : zamin,zamax
      use b2mod_b2cmpa
     & , only : am
      use b2mod_constants
     & , only : mp,ev
#else
      use b2mod_dimensions
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
     .        eirnxid,eirnyid,eirnatmid,eirnmolid,eirnionid
      integer vdims(4),idum(1)
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
        nstrad=size(eirene_mc_papl_sna_bal,4)
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
        status=nf_def_dim(ncid,'nion',size(eirene_mc_paio_sna_bal,3),
     &                                                  eirnionid)
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
*   ..bb,hx,hy,hz,vol,gs,am:
      call write_cdf_double(ncid,'bb',bb,vdims,3)
      call write_cdf_double(ncid,'hx',hx,vdims,2)
      call write_cdf_double(ncid,'hy',hy,vdims,2)
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
      vdims(4)=nsid
      call write_cdf_double(ncid,'fna_pinch',
     &                   fna_pinch_av,vdims,4)
      call write_cdf_double(ncid,'fna_pll',
     &                   fna_pll_av,vdims,4)
      call write_cdf_double(ncid,'fna_drift',
     &                   fna_drift_av,vdims,4)
      call write_cdf_double(ncid,'fna_ch',
     &                   fna_ch_av,vdims,4)
      call write_cdf_double(ncid,'fna_nanom',
     &                   fna_nanom_av,vdims,4)
      call write_cdf_double(ncid,'fna_panom',
     &                   fna_panom_av,vdims,4)
      call write_cdf_double(ncid,'fna_pschused',
     &                   fna_pschused_av,vdims,4)
      call write_cdf_double(ncid,'fna_tot',
     &                   fna_tot_av,vdims,4)
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
      call write_cdf_double(ncid,'tot_sna_bal',
     &                   tot_sna_bal_av,vdims,3)
*     Residual
      call write_cdf_double(ncid,'resco',resco_av,vdims,3)
*   ..Arrays for momentum balance:
*     Fluxes
      vdims(3)=n2id
      vdims(4)=nsid
      call write_cdf_double(ncid,'fmo_flua',fmo_flua_av,vdims,4)
      call write_cdf_double(ncid,'fmo_cvsa',fmo_cvsa_av,vdims,4)
      call write_cdf_double(ncid,'fmo_hybr',fmo_hybr_av,vdims,4)
      call write_cdf_double(ncid,'fmo_b2nxfv',fmo_b2nxfv_av,vdims,4)
      call write_cdf_double(ncid,'fmo_tot',fmo_tot_av,vdims,4)
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
      call write_cdf_double(ncid,'b2sifr_smotf_ehxp_bal',
     &                   b2sifr_smotf_ehxp_bal_av,vdims,3)
      call write_cdf_double(ncid,'b2sifr_smotf_cthe_bal',
     &                   b2sifr_smotf_cthe_bal_av,vdims,3)
      call write_cdf_double(ncid,'b2sifr_smotf_cthi_bal',
     &                   b2sifr_smotf_cthi_bal_av,vdims,3)
      call write_cdf_double(ncid,'b2sifr_smofrea_bal',
     &                      b2sifr_smofrea_bal_av,vdims,3)
      call write_cdf_double(ncid,'b2sifr_smofria_bal',
     &                      b2sifr_smofria_bal_av,vdims,3)
      call write_cdf_double(ncid,'b2sifr_smotfea_bal',
     &                      b2sifr_smotfea_bal_av,vdims,3)
      call write_cdf_double(ncid,'b2sifr_smotfia_bal',
     &                      b2sifr_smotfia_bal_av,vdims,3)
      call write_cdf_double(ncid,'b2siav_smovh_bal',
     &                   b2siav_smovh_bal_av,vdims,3)
      call write_cdf_double(ncid,'b2siav_smovv_bal',
     &                   b2siav_smovv_bal_av,vdims,3)
      call write_cdf_double(ncid,'b2siav_smovi_bal',
     &                   b2siav_smovi_bal_av,vdims,3)          !som 13.08.21
      call write_cdf_double(ncid,'b2sicf_smo_bal',
     &                   b2sicf_smo_bal_av,vdims,3)
      call write_cdf_double(ncid,'b2sian_smo_bal',
     &                   b2sian_smo_bal_av,vdims,3)            !srv 25.09.17
      call write_cdf_double(ncid,'b2nxdv_smo_bal',
     &                   b2nxdv_smo_bal_av,vdims,3)
      call write_cdf_double(ncid,'b2sigp_smogp_bal',
     &                   b2sigp_smogp_bal_av,vdims,3)
      call write_cdf_double(ncid,'b2sigp_smogpi_bal',
     &                   b2sigp_smogpi_bal_av,vdims,3)
      call write_cdf_double(ncid,'b2sigp_smogpe_bal',
     &                   b2sigp_smogpe_bal_av,vdims,3)
      call write_cdf_double(ncid,'b2sigp_smogpgr_bal',
     &                   b2sigp_smogpgr_bal_av,vdims,3)
      call write_cdf_double(ncid,'b2sigp_pstat_bal',
     &                   b2sigp_pstat_bal_av,vdims,3)
      call write_cdf_double(ncid,'b2sigp_pstati_bal',
     &                   b2sigp_pstati_bal_av,vdims,3)
      call write_cdf_double(ncid,'b2sigp_pstate_bal',
     &                   b2sigp_pstate_bal_av,vdims,2)
      call write_cdf_double(ncid,'tot_smo_bal',
     &                   tot_smo_bal_av,vdims,3)
*     Residual
      call write_cdf_double(ncid,'resmo',resmo_av,vdims,3)
*   ..Arrays for electron heat balance:
*     Fluxes
      vdims(3)=n2id
      call write_cdf_double(ncid,'fhe_32',
     &                   fhe_32_av,vdims,3)
      call write_cdf_double(ncid,'fhe_52',
     &                   fhe_52_av,vdims,3)
      call write_cdf_double(ncid,'fhe_thermj',
     &                   fhe_thermj_av,vdims,3)
      call write_cdf_double(ncid,'fhe_cond',
     &                   fhe_cond_av,vdims,3)
      call write_cdf_double(ncid,'fhe_dia',
     &                   fhe_dia_av,vdims,3)
      call write_cdf_double(ncid,'fhe_ecrb',
     &                   fhe_ecrb_av,vdims,3)
      call write_cdf_double(ncid,'fhe_strange',
     &                   fhe_strange_av,vdims,3)
      call write_cdf_double(ncid,'fhe_pschused',
     &                   fhe_pschused_av,vdims,3)
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
      vdims(3)=nsid
      call write_cdf_double(ncid,'b2stel_she_bal',
     &                   b2stel_she_bal_av,vdims,3)
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
      call write_cdf_double(ncid,'fhi_32',
     &                      fhi_32_av,vdims,3)
      call write_cdf_double(ncid,'fhi_52',
     &                      fhi_52_av,vdims,3)
      call write_cdf_double(ncid,'fhi_cond',
     &                      fhi_cond_av,vdims,3)
      call write_cdf_double(ncid,'fhi_dia',
     &                   fhi_dia_av,vdims,3)
      call write_cdf_double(ncid,'fhi_ecrb',
     &                   fhi_ecrb_av,vdims,3)
      call write_cdf_double(ncid,'fhi_strange',
     &                   fhi_strange_av,vdims,3)
      call write_cdf_double(ncid,'fhi_pschused',
     &                   fhi_pschused_av,vdims,3)
      call write_cdf_double(ncid,'fhi_inert',
     &                   fhi_inert_av,vdims,3)
      call write_cdf_double(ncid,'fhi_vispar',
     &                   fhi_vispar_av,vdims,3)
      call write_cdf_double(ncid,'fhi_visper',
     &                   fhi_visper_av,vdims,3)
      call write_cdf_double(ncid,'fhi_visq',
     &                   fhi_visq_av,vdims,3)
      call write_cdf_double(ncid,'fhi_anml',
     &                   fhi_anml_av,vdims,3)
      call write_cdf_double(ncid,'fhi_kevis',
     &                   fhi_kevis_av,vdims,3)
      call write_cdf_double(ncid,'fhi_vzh',
     &                   fhi_vzh_av,vdims,3)
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
      call write_cdf_double(ncid,'b2stel_she_ion_bal',
     &                      b2stel_she_ion_bal_av,vdims,2)
      call write_cdf_double(ncid,'b2stel_she_rec_bal',
     &                      b2stel_she_rec_bal_av,vdims,2)
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
        call write_cdf_double(ncid,'eirene_mc_mppl_smo_bal',
     &                     eirene_mc_mppl_smo_bal_av,vdims,4)
        vdims(3)=eirnatmid
        call write_cdf_double(ncid,'eirene_mc_paat_sna_bal',
     &                     eirene_mc_paat_sna_bal_av,vdims,4)
        call write_cdf_double(ncid,'eirene_mc_pmat_sna_bal',
     &                     eirene_mc_pmat_sna_bal_av,vdims,4)
        call write_cdf_double(ncid,'eirene_mc_piat_sna_bal',
     &                     eirene_mc_piat_sna_bal_av,vdims,4)
        vdims(3)=eirnmolid
        call write_cdf_double(ncid,'eirene_mc_paml_sna_bal',
     &                     eirene_mc_paml_sna_bal_av,vdims,4)
        call write_cdf_double(ncid,'eirene_mc_pmml_sna_bal',
     &                     eirene_mc_pmml_sna_bal_av,vdims,4)
        call write_cdf_double(ncid,'eirene_mc_piml_sna_bal',
     &                     eirene_mc_piml_sna_bal_av,vdims,4)
        vdims(3)=eirnionid
        call write_cdf_double(ncid,'eirene_mc_paio_sna_bal',
     &                     eirene_mc_paio_sna_bal_av,vdims,4)
        call write_cdf_double(ncid,'eirene_mc_pmio_sna_bal',
     &                     eirene_mc_pmio_sna_bal_av,vdims,4)
        call write_cdf_double(ncid,'eirene_mc_piio_sna_bal',
     &                     eirene_mc_piio_sna_bal_av,vdims,4)
        vdims(3)=nstraid
        call write_cdf_double(ncid,'eirene_mc_pael_sne_bal',
     &                     eirene_mc_pael_sne_bal_av,vdims,3)
        call write_cdf_double(ncid,'eirene_mc_pmel_sne_bal',
     &                     eirene_mc_pmel_sne_bal_av,vdims,3)
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
