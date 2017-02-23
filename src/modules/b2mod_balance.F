      module b2mod_balance
      use b2mod_types
      implicit none

*   ..Flux and source arrays required for balance post-processing routines. For the
*     sources, the linearised components are stored and then the total source is updated
*     just before each residual calculation. Only this total source is output to balance.nc
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
*     Sources
     *  b2stbr_sna0to1(:,:,:,:),
     *  b2stbr_sna_bal(:,:,:),
     *  b2stbc_sna0to1(:,:,:,:),
     *  b2stbc_sna_bal(:,:,:),
     *  eirene_mc_papl_sna0to1(:,:,:,:,:),
     *  eirene_mc_pmpl_sna0to1(:,:,:,:,:),
     *  eirene_mc_pipl_sna0to1(:,:,:,:,:),
     *  eirene_mc_pppl_sna0to1(:,:,:,:,:),
     *  eirene_mc_core_sna0to1(:,:,:,:),
     *  eirene_mc_papl_sna_bal(:,:,:,:),
     *  eirene_mc_pmpl_sna_bal(:,:,:,:),
     *  eirene_mc_pipl_sna_bal(:,:,:,:),
     *  eirene_mc_pppl_sna_bal(:,:,:,:),
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
*   ..Arrays for momentum balance:
*     Fluxes
     *  fmo_flua(:,:,:,:),
     *  fmo_cvsa(:,:,:,:),
     *  fmo_b2nxfv(:,:,:,:),
*     Sources
     *  b2stbr_smo0to3(:,:,:,:),
     *  b2stbr_smo_bal(:,:,:),
     *  b2stbc_smo0to3(:,:,:,:),
     *  b2stbc_smo_bal(:,:,:),
     *  eirene_mc_mapl_smo0to3(:,:,:,:,:),
     *  eirene_mc_mmpl_smo0to3(:,:,:,:,:),
     *  eirene_mc_mipl_smo0to3(:,:,:,:,:),
     *  eirene_mc_cppv_smo0to3(:,:,:,:,:),
     *  eirene_mc_mapl_smo_bal(:,:,:,:),
     *  eirene_mc_mmpl_smo_bal(:,:,:,:),
     *  eirene_mc_mipl_smo_bal(:,:,:,:),
     *  eirene_mc_cppv_smo_bal(:,:,:,:),
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
     *  b2sigp_pstat_bal(:,:,:),
*   ..Arrays for electron heat balance:
*     Fluxes
     *  fhe_32(:,:,:),
     *  fhe_52(:,:,:),
     *  fhe_thermj(:,:,:),
     *  fhe_cond(:,:,:),
     *  fhe_dia(:,:,:),
     *  fhe_ecrb(:,:,:),
     *  fhe_strange(:,:,:),
     *  fhe_pschused(:,:,:),
*     Sources
     *  b2stbr_she0to3(:,:,:),
     *  b2stbr_she_bal(:,:),
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
     *  b2npht_shei_bal(:,:), !djm equipartition is not linearised
*   ..Arrays for ion heat balance:
*     Fluxes
     *  fhi_32(:,:,:),
     *  fhi_52(:,:,:),
     *  fhi_cond(:,:,:),
     *  fhi_dia(:,:,:),
     *  fhi_ecrb(:,:,:),
     *  fhi_strange(:,:,:),
     *  fhi_pschused(:,:,:),
     *  fhi_inert(:,:,:),
     *  fhi_vispar(:,:,:),
     *  fhi_anml(:,:,:),
     *  fhi_kevis(:,:,:),
*     Sources
     *  b2stbr_shi0to3(:,:,:),
     *  b2stbr_shi_bal(:,:),
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

*   ..Balance_netcdf decides whether to initialise, track, and write balance arrays.
*     Default: do not create balance arrays (balance_netcdf.eq.0)
      integer, save :: balance_netcdf
      data balance_netcdf/0/

      contains

      subroutine alloc_b2mod_balance(nxd,nyd,nsd)
      implicit none
*   ..Input arguments (unchanged on exit)
      integer nxd,nyd,nsd
*   ..Procedures
      external ipgeti

*   ..Read input variable balance_netcdf
      call ipgeti ('balance_netcdf', balance_netcdf)

      if (balance_netcdf.ne.0) then
        write(*,*) 'balance_netcdf.ne.0. Allocating balance arrays'
*     ..Arrays for particle balance:
*       Fluxes:
        allocate(fna_pinch(-1:nxd,-1:nyd,0:1,0:nsd-1))
        allocate(fna_pll(-1:nxd,-1:nyd,0:1,0:nsd-1))
        allocate(fna_drift(-1:nxd,-1:nyd,0:1,0:nsd-1))
        allocate(fna_ch(-1:nxd,-1:nyd,0:1,0:nsd-1))
        allocate(fna_nanom(-1:nxd,-1:nyd,0:1,0:nsd-1))
        allocate(fna_panom(-1:nxd,-1:nyd,0:1,0:nsd-1))
        allocate(fna_pschused(-1:nxd,-1:nyd,0:1,0:nsd-1))
*       Sources
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
*     ..Arrays for momentum balance:
*       Fluxes
        allocate(fmo_flua(-1:nxd,-1:nyd,0:1,0:nsd-1))
        allocate(fmo_cvsa(-1:nxd,-1:nyd,0:1,0:nsd-1))
        allocate(fmo_b2nxfv(-1:nxd,-1:nyd,0:1,0:nsd-1))
*       Sources
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
*     ..Arrays for electron heat balance:
*       Fluxes
        allocate(fhe_32(-1:nxd,-1:nyd,0:1))
        allocate(fhe_52(-1:nxd,-1:nyd,0:1))
        allocate(fhe_thermj(-1:nxd,-1:nyd,0:1))
        allocate(fhe_cond(-1:nxd,-1:nyd,0:1))
        allocate(fhe_dia(-1:nxd,-1:nyd,0:1))
        allocate(fhe_ecrb(-1:nxd,-1:nyd,0:1))
        allocate(fhe_strange(-1:nxd,-1:nyd,0:1))
        allocate(fhe_pschused(-1:nxd,-1:nyd,0:1))
*       Sources
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
*     ..Arrays for ion heat balance:
*       Fluxes
        allocate(fhi_32(-1:nxd,-1:nyd,0:1))
        allocate(fhi_52(-1:nxd,-1:nyd,0:1))
        allocate(fhi_cond(-1:nxd,-1:nyd,0:1))
        allocate(fhi_dia(-1:nxd,-1:nyd,0:1))
        allocate(fhi_ecrb(-1:nxd,-1:nyd,0:1))
        allocate(fhi_strange(-1:nxd,-1:nyd,0:1))
        allocate(fhi_pschused(-1:nxd,-1:nyd,0:1))
        allocate(fhi_inert(-1:nxd,-1:nyd,0:1))
        allocate(fhi_vispar(-1:nxd,-1:nyd,0:1))
        allocate(fhi_anml(-1:nxd,-1:nyd,0:1))
        allocate(fhi_kevis(-1:nxd,-1:nyd,0:1))
*       Sources
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
      endif

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
*   ..Arrays for momentum balance:
*     Fluxes
      deallocate(fmo_flua)
      deallocate(fmo_cvsa)
      deallocate(fmo_b2nxfv)
*     Sources
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
      deallocate(fhi_anml)
      deallocate(fhi_kevis)
*     Sources
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
      end subroutine dealloc_b2mod_balance

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
      allocate(eirene_mc_papl_sna_bal(-1:nxd,-1:nyd,0:nsd-1,1:nstrad))
      allocate(eirene_mc_pmpl_sna_bal(-1:nxd,-1:nyd,0:nsd-1,1:nstrad))
      allocate(eirene_mc_pipl_sna_bal(-1:nxd,-1:nyd,0:nsd-1,1:nstrad))
      allocate(eirene_mc_pppl_sna_bal(-1:nxd,-1:nyd,0:nsd-1,1:nstrad))
      allocate(eirene_mc_core_sna_bal(-1:nxd,-1:nyd,0:nsd-1))
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
cdjm Jan2017 Eirene arrays for momentum balance:
      allocate(
     & eirene_mc_mapl_smo0to3(-1:nxd,-1:nyd,0:3,0:nsd-1,1:nstrad))
      allocate(
     & eirene_mc_mmpl_smo0to3(-1:nxd,-1:nyd,0:3,0:nsd-1,1:nstrad))
      allocate(
     & eirene_mc_mipl_smo0to3(-1:nxd,-1:nyd,0:3,0:nsd-1,1:nstrad))
      allocate(
     & eirene_mc_cppv_smo0to3(-1:nxd,-1:nyd,0:3,0:nsd-1,1:nstrad))
      allocate(eirene_mc_mapl_smo_bal(-1:nxd,-1:nyd,0:nsd-1,1:nstrad))
      allocate(eirene_mc_mmpl_smo_bal(-1:nxd,-1:nyd,0:nsd-1,1:nstrad))
      allocate(eirene_mc_mipl_smo_bal(-1:nxd,-1:nyd,0:nsd-1,1:nstrad))
      allocate(eirene_mc_cppv_smo_bal(-1:nxd,-1:nyd,0:nsd-1,1:nstrad))
      eirene_mc_mapl_smo0to3=0.0_R8
      eirene_mc_mmpl_smo0to3=0.0_R8
      eirene_mc_mipl_smo0to3=0.0_R8
      eirene_mc_cppv_smo0to3=0.0_R8
      eirene_mc_mapl_smo_bal=0.0_R8
      eirene_mc_mmpl_smo_bal=0.0_R8
      eirene_mc_mipl_smo_bal=0.0_R8
      eirene_mc_cppv_smo_bal=0.0_R8
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
          b2stbr_sna_bal(ix,iy,isb)=b2stbr_sna0to1(ix,iy,0,isb)+
     &                     b2stbr_sna0to1(ix,iy,1,isb)*nb(ix,iy)
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
*   ..procedures
      intrinsic size

      nstrad = size(eirene_mc_mapl_smo_bal,4)

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
              eirene_mc_cppv_smo_bal(ix,iy,isb,istrai)=
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
      integer nx,ny
      real (kind=R8):: ne(-1:nx,-1:ny), ni(-1:nx,-1:ny),
     &                 te(-1:nx,-1:ny), ti(-1:nx,-1:ny)
*   ..local variables
      integer ix,iy,istrai,nstrad
*   ..procedures
      intrinsic size

      nstrad = size(eirene_mc_eael_she_bal,3)

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
          if (allocated(eirene_mc_eael_she_bal)) then
            do istrai=1,nstrad
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

#ifndef NO_CDF
      subroutine write_balance(nxd,nyd,nsd)
      use b2mod_geo
     & , only : crx,cry,hx,bb,vol
      use b2mod_indirect
     & , only : leftix,leftiy,rightix,rightiy,
     &          topix,topiy,bottomix,bottomiy
      use b2mod_b2cmpa
     & , only : zamin,zamax
      use b2mod_residuals
     & , only : resco, resmo, reshe, reshi
      implicit none
#include <netcdf.inc>
*   ..input arguments (unchanged on exit)
      integer nxd,nyd,nsd
*   ..local variables
      integer nstrad
      integer status
      integer ncid,nxid,nyid,n1id,n2id,n4id,nsid,nstraid,specid,varid
      integer vdims(4)
      integer jxi,jxa,jsep
      integer is
      character*13, allocatable :: species_name(:)
      integer use_eirene
      integer, allocatable :: za(:)
*   ..procedures
      intrinsic size

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

*   ..take dataset out of define mode
      status=nf_enddef(ncid)
      call check_cdf_status(status)

*   ..crx, cry:
      vdims(1)=nxid
      vdims(2)=nyid
      vdims(3)=n4id
      call write_cdf_double(ncid,'crx',crx,vdims,3)
      call write_cdf_double(ncid,'cry',cry,vdims,3)
*   ..bb,hx,vol:
      call write_cdf_double(ncid,'bb',bb,vdims,3)
      call write_cdf_double(ncid,'hx',hx,vdims,2)
      call write_cdf_double(ncid,'vol',vol,vdims,2)
*   ..rightix, rightiy etc:
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
      call write_cdf_int_single(ncid,'jxi',jxi,n1id)
      call write_cdf_int_single(ncid,'jxa',jxa,n1id)
      call write_cdf_int_single(ncid,'jsep',jsep,n1id)
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
      call ipgeti ('b2mndr_eirene',use_eirene)
      call write_cdf_int_single(ncid,'b2mndr_eirene',use_eirene,n1id)
*   ..za:
      vdims(1)=nsid
      allocate(za(0:nsd-1))
      za(0:nsd-1)=(zamin(0:nsd-1)+zamax(0:nsd-1))/2.0_R8
      call write_cdf_int(ncid,'za',za(0:nsd-1),vdims,1)
*   ..Arrays for particle balance:
*     Fluxes
      vdims(1)=nxid
      vdims(2)=nyid
      vdims(3)=n2id
      vdims(4)=nsid
      call write_cdf_double(ncid,'fna_pinch',
     &                   fna_pll,vdims,4)
      call write_cdf_double(ncid,'fna_pll',
     &                   fna_drift,vdims,4)
      call write_cdf_double(ncid,'fna_drift',
     &                   fna_drift,vdims,4)
      call write_cdf_double(ncid,'fna_ch',
     &                   fna_ch,vdims,4)
      call write_cdf_double(ncid,'fna_nanom',
     &                   fna_nanom,vdims,4)
      call write_cdf_double(ncid,'fna_panom',
     &                   fna_panom,vdims,4)
      call write_cdf_double(ncid,'fna_pschused',
     &                   fna_pschused,vdims,4)
*     Sources
      vdims(3)=nsid
      vdims(4)=nstraid
      call write_cdf_double(ncid,'b2stbr_sna_bal',
     &                   b2stbr_sna_bal,vdims,3)
      call write_cdf_double(ncid,'b2stbc_sna_bal',
     &                   b2stbc_sna_bal,vdims,3)
      call write_cdf_double(ncid,'b2stbm_sna_bal',
     &                   b2stbm_sna_bal,vdims,3)
      call write_cdf_double(ncid,'ext_sna_bal',
     &                   ext_sna_bal,vdims,3)
      call write_cdf_double(ncid,'b2stel_sna_ion_bal',
     &                   b2stel_sna_ion_bal,vdims,3)
      call write_cdf_double(ncid,'b2stel_sna_rec_bal',
     &                   b2stel_sna_rec_bal,vdims,3)
      call write_cdf_double(ncid,'b2stcx_sna_bal',
     &                   b2stcx_sna_bal,vdims,3)
      call write_cdf_double(ncid,'b2srsm_sna_bal',
     &                   b2srsm_sna_bal,vdims,3)
      call write_cdf_double(ncid,'b2srdt_sna_bal',
     &                   b2srdt_sna_bal,vdims,3)
      call write_cdf_double(ncid,'b2srst_sna_bal',
     &                   b2srst_sna_bal,vdims,3)
*     Residual
      call write_cdf_double(ncid,'resco',resco,vdims,3)
*   ..Arrays for momentum balance:
*     Fluxes
      vdims(3)=n2id
      vdims(4)=nsid
      call write_cdf_double(ncid,'fmo_flua',fmo_flua,vdims,4)
      call write_cdf_double(ncid,'fmo_cvsa',fmo_cvsa,vdims,4)
      call write_cdf_double(ncid,'fmo_b2nxfv',fmo_b2nxfv,vdims,4)
*     Sources
      vdims(3)=nsid
      vdims(4)=nstraid
      call write_cdf_double(ncid,'b2stbr_smo_bal',
     &                   b2stbr_smo_bal,vdims,3)
      call write_cdf_double(ncid,'b2stbc_smo_bal',
     &                   b2stbc_smo_bal,vdims,3)
      call write_cdf_double(ncid,'b2stbm_smo_bal',
     &                   b2stbm_smo_bal,vdims,3)
      call write_cdf_double(ncid,'ext_smo_bal',
     &                   ext_smo_bal,vdims,3)
      call write_cdf_double(ncid,'b2stel_smq_ion_bal',
     &                   b2stel_smq_ion_bal,vdims,3)
      call write_cdf_double(ncid,'b2stel_smq_rec_bal',
     &                   b2stel_smq_rec_bal,vdims,3)
      call write_cdf_double(ncid,'b2stcx_smq_bal',
     &                   b2stcx_smq_bal,vdims,3)
      call write_cdf_double(ncid,'b2srsm_smo_bal',
     &                   b2srsm_smo_bal,vdims,3)
      call write_cdf_double(ncid,'b2srdt_smo_bal',
     &                   b2srdt_smo_bal,vdims,3)
      call write_cdf_double(ncid,'b2srst_smo_bal',
     &                   b2srst_smo_bal,vdims,3)
      call write_cdf_double(ncid,'b2sifr_smoch_bal',
     &                   b2sifr_smoch_bal,vdims,3)
      call write_cdf_double(ncid,'b2sifr_smotf_bal',
     &                   b2sifr_smotf_bal,vdims,3)
      call write_cdf_double(ncid,'b2siav_smovh_bal',
     &                   b2siav_smovh_bal,vdims,3)
      call write_cdf_double(ncid,'b2siav_smovv_bal',
     &                   b2siav_smovv_bal,vdims,3)
      call write_cdf_double(ncid,'b2sicf_smo_bal',
     &                   b2sicf_smo_bal,vdims,3)
      call write_cdf_double(ncid,'b2nxdv_smo_bal',
     &                   b2nxdv_smo_bal,vdims,3)
      call write_cdf_double(ncid,'b2sigp_smogp_bal',
     &                   b2sigp_smogp_bal,vdims,3)
      call write_cdf_double(ncid,'b2sigp_pstat_bal',
     &                   b2sigp_pstat_bal,vdims,3)
*     Residual
      call write_cdf_double(ncid,'resmo',resmo,vdims,3)
*   ..Arrays for electron heat balance:
*     Fluxes
      vdims(3)=n2id
      call write_cdf_double(ncid,'fhe_32',
     &                   fhe_32,vdims,3)
      call write_cdf_double(ncid,'fhe_52',
     &                   fhe_52,vdims,3)
      call write_cdf_double(ncid,'fhe_thermj',
     &                   fhe_thermj,vdims,3)
      call write_cdf_double(ncid,'fhe_cond',
     &                   fhe_cond,vdims,3)
      call write_cdf_double(ncid,'fhe_dia',
     &                   fhe_dia,vdims,3)
      call write_cdf_double(ncid,'fhe_ecrb',
     &                   fhe_ecrb,vdims,3)
      call write_cdf_double(ncid,'fhe_strange',
     &                   fhe_strange,vdims,3)
      call write_cdf_double(ncid,'fhe_pschused',
     &                   fhe_pschused,vdims,3)
*     Sources
      vdims(3)=nstraid
      call write_cdf_double(ncid,'b2stbr_she_bal',
     &                   b2stbr_she_bal,vdims,2)
      call write_cdf_double(ncid,'b2stbc_she_bal',
     &                   b2stbc_she_bal,vdims,2)
      call write_cdf_double(ncid,'b2stbm_she_bal',
     &                   b2stbm_she_bal,vdims,2)
      call write_cdf_double(ncid,'ext_she_bal',
     &                   ext_she_bal,vdims,2)
      call write_cdf_double(ncid,'b2stel_she_bal',
     &                   b2stel_she_bal,vdims,2)
      call write_cdf_double(ncid,'b2srsm_she_bal',
     &                   b2srsm_she_bal,vdims,2)
      call write_cdf_double(ncid,'b2srdt_she_bal',
     &                   b2srdt_she_bal,vdims,2)
      call write_cdf_double(ncid,'b2srst_she_bal',
     &                   b2srst_she_bal,vdims,2)
      call write_cdf_double(ncid,'b2sihs_diae_bal',
     &                   b2sihs_diae_bal,vdims,2)
      call write_cdf_double(ncid,'b2sihs_divue_bal',
     &                   b2sihs_divue_bal,vdims,2)
      call write_cdf_double(ncid,'b2sihs_exbe_bal',
     &                   b2sihs_exbe_bal,vdims,2)
      call write_cdf_double(ncid,'b2sihs_joule_bal',
     &                   b2sihs_joule_bal,vdims,2)
      call write_cdf_double(ncid,'b2npht_shei_bal',
     &                   b2npht_shei_bal,vdims,2)
*     Residual
      call write_cdf_double(ncid,'reshe',reshe,vdims,2)
*   ..Arrays for ion heat balance:
*     Fluxes
      vdims(3)=n2id
      call write_cdf_double(ncid,'fhi_32',
     &                      fhi_32,vdims,3)
      call write_cdf_double(ncid,'fhi_52',
     &                      fhi_52,vdims,3)
      call write_cdf_double(ncid,'fhi_cond',
     &                      fhi_cond,vdims,3)
      call write_cdf_double(ncid,'fhi_dia',
     &                   fhi_dia,vdims,3)
      call write_cdf_double(ncid,'fhi_ecrb',
     &                   fhi_ecrb,vdims,3)
      call write_cdf_double(ncid,'fhi_strange',
     &                   fhi_strange,vdims,3)
      call write_cdf_double(ncid,'fhi_pschused',
     &                   fhi_pschused,vdims,3)
      call write_cdf_double(ncid,'fhi_inert',
     &                   fhi_inert,vdims,3)
      call write_cdf_double(ncid,'fhi_vispar',
     &                   fhi_vispar,vdims,3)
      call write_cdf_double(ncid,'fhi_anml',
     &                   fhi_anml,vdims,3)
      call write_cdf_double(ncid,'fhi_kevis',
     &                   fhi_kevis,vdims,3)
*     Sources
      vdims(3)=nstraid
      call write_cdf_double(ncid,'b2stbr_shi_bal',
     &                      b2stbr_shi_bal,vdims,2)
      call write_cdf_double(ncid,'b2stbc_shi_bal',
     &                      b2stbc_shi_bal,vdims,2)
      call write_cdf_double(ncid,'b2stbm_shi_bal',
     &                      b2stbm_shi_bal,vdims,2)
      call write_cdf_double(ncid,'ext_shi_bal',
     &                      ext_shi_bal,vdims,2)
      call write_cdf_double(ncid,'b2stel_shi_ion_bal',
     &                      b2stel_shi_ion_bal,vdims,2)
      call write_cdf_double(ncid,'b2stel_shi_rec_bal',
     &                      b2stel_shi_rec_bal,vdims,2)
      call write_cdf_double(ncid,'b2stcx_shi_bal',
     &                      b2stcx_shi_bal,vdims,2)
      call write_cdf_double(ncid,'b2srsm_shi_bal',
     &                      b2srsm_shi_bal,vdims,2)
      call write_cdf_double(ncid,'b2srdt_shi_bal',
     &                      b2srdt_shi_bal,vdims,2)
      call write_cdf_double(ncid,'b2srst_shi_bal',
     &                      b2srst_shi_bal,vdims,2)
      call write_cdf_double(ncid,'b2sihs_diaa_bal',
     &                      b2sihs_diaa_bal,vdims,2)
      call write_cdf_double(ncid,'b2sihs_divua_bal',
     &                      b2sihs_divua_bal,vdims,2)
      call write_cdf_double(ncid,'b2sihs_exba_bal',
     &                      b2sihs_exba_bal,vdims,2)
      call write_cdf_double(ncid,'b2sihs_visa_bal',
     &                      b2sihs_visa_bal,vdims,2)
      call write_cdf_double(ncid,'b2sihs_fraa_bal',
     &                      b2sihs_fraa_bal,vdims,2)
*     Residual
      call write_cdf_double(ncid,'reshi',reshi,vdims,2)
*   ..Eirene sources
      if (allocated(eirene_mc_papl_sna_bal)) then
        call write_cdf_double(ncid,'eirene_mc_papl_sna_bal',
     &                     eirene_mc_papl_sna_bal,vdims,4)
        call write_cdf_double(ncid,'eirene_mc_pmpl_sna_bal',
     &                     eirene_mc_pmpl_sna_bal,vdims,4)
        call write_cdf_double(ncid,'eirene_mc_pipl_sna_bal',
     &                     eirene_mc_pipl_sna_bal,vdims,4)
        call write_cdf_double(ncid,'eirene_mc_pppl_sna_bal',
     &                     eirene_mc_pppl_sna_bal,vdims,4)
        call write_cdf_double(ncid,'eirene_mc_core_sna_bal',
     &                     eirene_mc_core_sna_bal,vdims,4)
        call write_cdf_double(ncid,'eirene_mc_mapl_smo_bal',
     &                     eirene_mc_mapl_smo_bal,vdims,4)
        call write_cdf_double(ncid,'eirene_mc_mmpl_smo_bal',
     &                     eirene_mc_mmpl_smo_bal,vdims,4)
        call write_cdf_double(ncid,'eirene_mc_mipl_smo_bal',
     &                     eirene_mc_mipl_smo_bal,vdims,4)
        call write_cdf_double(ncid,'eirene_mc_cppv_smo_bal',
     &                     eirene_mc_cppv_smo_bal,vdims,4)
        call write_cdf_double(ncid,'eirene_mc_eael_she_bal',
     &                     eirene_mc_eael_she_bal,vdims,3)
        call write_cdf_double(ncid,'eirene_mc_emel_she_bal',
     &                     eirene_mc_emel_she_bal,vdims,3)
        call write_cdf_double(ncid,'eirene_mc_eiel_she_bal',
     &                     eirene_mc_eiel_she_bal,vdims,3)
        call write_cdf_double(ncid,'eirene_mc_epel_she_bal',
     &                     eirene_mc_epel_she_bal,vdims,3)
        call write_cdf_double(ncid,'eirene_mc_eapl_shi_bal',
     &                     eirene_mc_eapl_shi_bal,vdims,3)
        call write_cdf_double(ncid,'eirene_mc_empl_shi_bal',
     &                     eirene_mc_empl_shi_bal,vdims,3)
        call write_cdf_double(ncid,'eirene_mc_eipl_shi_bal',
     &                     eirene_mc_eipl_shi_bal,vdims,3)
        call write_cdf_double(ncid,'eirene_mc_eppl_shi_bal',
     &                     eirene_mc_eppl_shi_bal,vdims,3)
      else
        write(*,*) 'skipping Eirene sources'
      endif

      status=nf_close(ncid)
      call check_cdf_status(status)

      write(*,*) 'balance.nc written'

      return
      end subroutine write_balance

      subroutine write_cdf_double(ncid,varname,var,vdims,ndim)
      implicit none
#include <netcdf.inc>
*   ..input arguments
      integer :: ncid,vdims(*),ndim
      character*(*) :: varname
      real (kind=R8) :: var(*)
*   ..local variables
      integer status,varid

*   ..put netcdf dataset into define mode
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

      subroutine write_cdf_int(ncid,varname,var,vdims,ndim)
      implicit none
#include <netcdf.inc>
*   ..input arguments
      integer :: ncid,vdims(*),ndim
      character*(*) :: varname
      integer :: var(*)
*   ..local variables
      integer status,varid

*   ..put netcdf dataset into define mode
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
      integer :: ncid,var,vdim
      character*(*) :: varname
*   ..local variables
      integer status,varid

*   ..put netcdf dataset into define mode
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

*   ..put netcdf dataset into define mode
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
