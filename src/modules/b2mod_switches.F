      module b2mod_switches

c This module centralizes all the switches available in the B2.5 side of SOLPS-ITER provided by b2mn.dat
c First, the defaults are loaded (call set_defaults_switches in b2mod_driver).
c Next, the defaults are overwritten by the values provided in b2mn.dat (call read_switches in b2mod_driver).
c Finally, consistency checks on the provided switches are performed. Both the values themselves (check_values_switches) as well as
c consistency checks for combinations of different switches (check_consistency_switches) are performed here.

      use b2mod_types
      use b2mod_constants
     & , only : ev, c
      implicit none

      private
      public :: alloc_switches, dealloc_switches,
     .          read_switches, set_defaults_switches,
     .          check_switches, check_values_switches,
     .          check_consistency_switches

      type, public :: switches

      integer :: get_residuals !srv 22.05.18
      integer :: mdf_fhe       !srv 13.10.06
      integer :: mdf_fhi
      integer :: mdf_fnb
      integer :: pflux_cor
      integer :: pot_eq          !srv 13.07.05
      integer :: no_current      !srv 29.01.20
      integer :: sheath_drift_fix
      integer :: secmodel
      integer :: b2stbc_fchy_dia_coreonly
      integer :: nsmin
      integer :: nsmax
      integer :: icase_ii
      integer :: icase_ee
      integer :: icase_ei
      integer :: stylec
      integer :: zero_dead_and_core
      integer :: fke_Zhdanov
      integer :: conductive_limit
      integer :: core_cond_limit
      integer :: lluciani
      integer :: no_solve
      integer :: no_b2sral_call
      integer :: do_2nd_b2npco_call
      integer :: re_eval_prtls_fluxes
      integer :: istyle_cur_contr_on_S_and_N
      integer :: use_astra
      integer :: lin_shi0
      integer :: flux_limit_to_dpa
      integer :: flux_limit_to_vsa
      integer :: csig_an_style
      integer :: user_transport

c Advanced fluid neutral models
      integer :: tn_style
      integer :: recycle_afn
      integer :: use_uy_uz_0
      integer :: use_uy_uz_1
      integer :: remove_FC_el
      integer :: spatial_hybrid
      integer :: use_auto_spatial_hyb
      integer :: diffusion_bc_safeguard
      integer :: do_volrec
      integer :: transport_afn !! use advanced fluid neutral (afn) transport coefficients

      real (kind=R8) :: fhe_vis_per  !srv 06.07.00
      real (kind=R8) :: fhe_vis_par  !srv 22.03.02
      real (kind=R8) :: fnb_vis_per  !srv 02.10.17
      real (kind=R8) :: fnb_vis_q    !srv 02.10.17
      real (kind=R8) :: fhe_vis_q
      real (kind=R8) :: PSch         !srv 17.01.07
      real (kind=R8) :: fhePSch      !srv 17.01.07
      real (kind=R8) :: fhiPSch
      real (kind=R8) :: b2tfnb_xcur         !xpb 27.09.17
      real (kind=R8) :: b2tfnb_ycur         !srv 07.07.08
      real (kind=R8) :: fch_ion_neutral
      real (kind=R8) :: art_rad
      real (kind=R8) :: sna0ep
      real (kind=R8) :: she0ep
      real (kind=R8) :: shi0ep
      real (kind=R8) :: b2stbc_neoclassical
      real (kind=R8) :: integral_current
      real (kind=R8) :: b2stbc_fchy_dia
      real (kind=R8) :: b2mndr_na_min
      real (kind=R8) :: b2mndr_ua_max !! in absolute sense: abs(ua) < ua_max
      real (kind=R8) :: b2mndr_na_new
      real (kind=R8) :: cfc0
      real (kind=R8) :: prl_cur
      real (kind=R8) :: dia_cur
      real (kind=R8) :: fch_pTe
      real (kind=R8) :: alfTeEh
      real (kind=R8) :: fhe_vdia_par
      real (kind=R8) :: fhi_vdia_par
      real (kind=R8) :: xvecrb
      real (kind=R8) :: xwdia
      real (kind=R8) :: fnb_drift_hyb
      real (kind=R8) :: xfac
      real (kind=R8) :: fch_inert
      real (kind=R8) :: cthiv
      real (kind=R8) :: cthev
      real (kind=R8) :: b2trcl_lambda
      real (kind=R8) :: poleldr
      real (kind=R8) :: radeldr
      real (kind=R8) :: poldidr
      real (kind=R8) :: raddidr
      real (kind=R8) :: b2tlv0_alpha
      real (kind=R8) :: b2tlv0_gamma
      real (kind=R8) :: b2tlh0_alpha
      real (kind=R8) :: b2tlh0_gamma
      real (kind=R8) :: b2tlc0_alpha
      real (kind=R8) :: b2tlc0_gamma
      real (kind=R8) :: b2tfnb_alpha
      real (kind=R8) :: b2tfnb_gamma
      real (kind=R8) :: fch_stochastic
      real (kind=R8) :: Qalfmin
      real (kind=R8) :: delpo
      real (kind=R8) :: Kn_b1
      real (kind=R8) :: Kn_b2
      real (kind=R8) :: pol_anom_scale
      real (kind=R8) :: auto_spatial_hyb_Kn_1
      real (kind=R8) :: auto_spatial_hyb_Kn_2
      real (kind=R8) :: L_macro_AFN
      real (kind=R8) :: neutral_sources_rescale

c iout switches
      integer :: iout_b2wdat
      integer :: b2sian_iout
      integer :: b2sicf_iout
      integer :: b2sigp_iout
      integer :: b2sihs_iout
      integer :: b2srdt_iout
      integer :: b2srsm_iout
      integer :: b2srst_iout
      integer :: b2stcx_iout
      integer :: b2stel_iout
      integer :: b2stbr_iout
      integer :: b2sral_iout
      integer :: b2stbr_output
      integer :: b2npht_iout
      integer :: b2npmo_iout
      integer :: b2npco_iout
      integer :: b2npp7_iout
      integer :: b2nppo_iout
!      integer :: b2nph9_iout       ! No longer needed for unstructured grids
      integer :: b2news_iout
      integer :: b2nxfx_iout
      integer :: b2nxfc_iout        !srv 17.06.02
      integer :: b2usco_iout        !srv 13.10.06
      integer :: b2uspo_iout
      integer :: b2usht_iout
      integer :: b2usmo_iout         !srv 17.06.02
      integer :: b2tfhi_iout
      integer :: b2tfhe_iout
      integer :: b2tfrn_iout
      integer :: b2tfch_iout
      integer :: b2tlmv_iout
      integer :: b2tqca_iout
      integer :: b2tqce_iout
      integer :: b2trcl_iout
      integer :: b2tfnb_iout
      integer :: b2tral_iout
      integer :: b2tlv0_iout
      integer :: b2tlh0_iout
      integer :: b2tlc0_iout
      integer :: b2tqna_iout
      integer :: b2tvspr_iout
      integer :: b2tvskt_iout
      integer :: b2stbc_iout
      integer :: b2trno_iout
      integer :: b2stbr_phys_iout
      integer :: b2wdat_append
      integer :: AFN_out
      integer :: copy_background_iout
      integer :: eirene_mc_iout
      integer :: b2siav_iout
      integer :: b2mndr_iout

c style  switches
      integer :: b2sicf_style
      integer :: b2sigp_style
      integer :: b2sigp_pressure_restriction
      integer :: b2sihs_style
      integer :: b2sral_style
      integer :: style_visc
      integer :: b2tfnb_style_int_vel
      integer :: b2sihs_style_int_vel
      integer :: istyle_Joule_heating
      integer :: b2sqcx_styl0
      integer :: b2stcx_styl0
      integer :: b2stel_styl0
      integer :: b2nxdv_style         !srv 26.10.08
      integer :: b2nxfx_style
      integer :: b2nxfc_style
      integer :: mode_ehy
      integer :: mode_ehx
      integer :: b2tlmv_style
      integer :: b2tqca_model
      integer :: b2tqce_model
      integer :: b2tqce_style_guard_cells
      integer :: b2tral_mode
      integer :: b2tlv0_style
      integer :: b2tlh0_style
      integer :: b2tlc0_style
      integer :: b2tinnt_style
      integer :: model_sig
      integer :: b2mwqt_style
      integer :: b2npht_style
      integer :: fnb_drift_style
      integer :: b2sifr_styl0

c phm switches
      real (kind=R8) :: b2sian_phm0
      real (kind=R8) :: b2sicf_phm0
      real (kind=R8) :: b2sicf_phm1
      real (kind=R8) :: b2sigp_phm0
      real (kind=R8) :: b2sihs_phm0  !! srv 27.12.01 28.08.02
      real (kind=R8) :: b2sihs_phm1
      real (kind=R8) :: b2sihs_phm2
      real (kind=R8) :: b2sihs_phm3
      real (kind=R8) :: b2sihs_phm4
      real (kind=R8) :: b2sihs_phm5
      real (kind=R8) :: b2sihs_phm6
      real (kind=R8) :: b2sihs_phm7
      real (kind=R8) :: b2sihs_phm8
      real (kind=R8) :: b2sqcx_phm0
      real (kind=R8) :: b2sqel_phm0
      real (kind=R8) :: b2sqel_phm1
      real (kind=R8) :: b2sqel_phm2
      real (kind=R8) :: b2stel_phm0
      real (kind=R8) :: b2srdt_phm0
      real (kind=R8) :: b2srdt_phm1
      real (kind=R8) :: b2srdt_phm3
      real (kind=R8) :: b2srdt_phm4
      real (kind=R8) :: b2srdt_phm5
      real (kind=R8) :: b2nxfv_phm0
      real (kind=R8) :: b2nxfv_phm1
      real (kind=R8) :: b2sifr_limthee
      real (kind=R8) :: b2sifr_limthii
      real (kind=R8) :: b2sifr_phm0
      real (kind=R8) :: b2sifr_phm1
      real (kind=R8) :: b2sifr_phm2
      real (kind=R8) :: b2sifr_phm3
      real (kind=R8) :: b2tqca_phm0
      real (kind=R8) :: b2treq_phm0
      real (kind=R8) :: b2stbc_phm0
      real (kind=R8) :: b2stbc_phm1
      real (kind=R8) :: b2stbc_phm2

c pcm switches
      real (kind=R8) :: b2npco_pcm0
      real (kind=R8) :: b2npco_pcm1
      real (kind=R8) :: b2npht_pcm0
      real (kind=R8) :: b2npht_pcm1

c relaxation factors
      real (kind=R8) :: b2stcx_rg0
      real (kind=R8) :: b2npco_rxg
      real (kind=R8) :: b2npht_rxg
      real (kind=R8) :: b2npmo_rxg
      real (kind=R8) :: xfm0
      real (kind=R8) :: xfm1
      real (kind=R8) :: xfm2
      real (kind=R8) :: xfm3
      real (kind=R8) :: b2stel_rxm0
      real (kind=R8) :: b2stel_rxm1
      real (kind=R8) :: b2stel_rxm2
      real (kind=R8) :: b2stel_rxm3
      real (kind=R8) :: b2stel_rxm4
      real (kind=R8) :: b2stel_rg0
      real (kind=R8) :: b2stel_rg1

c discretization methods
      integer :: b2npmo_discr_meth
      integer :: b2tfhe_discr_meth
      integer :: b2tfhi_discr_meth
      integer :: b2tfnb_discr_meth

c namelist switches
      integer :: b2stbc_boundary_namelist
      integer :: b2stbr_neutrals_namelist
      integer :: b2tqna_transport_namelist
      integer :: b2srdt_numerics_namelist
      integer :: b2optim_namelist

c diagnostics switches
      integer :: b2stbc_diagno
      integer :: b2npmo_diagno
      integer :: ion_vlct_restrict
      integer :: vlct_diagno
      integer :: b2tqna_diagno

      real(kind=R8) :: ion_vlct_restrict_M

c MMS
      integer:: set_na_numerical
      integer:: set_ua_numerical
      integer:: set_ti_numerical
      integer:: set_te_numerical
      integer:: set_po_numerical

c geometry
      integer :: fcHz_style
      integer :: vxHz_style
      integer :: use_9pt_stencil
      integer :: fcVol_style
      integer :: vxVol_style

c Eirene coupling
      integer :: use_eirene
      integer :: eirene_sheath_pot
      integer :: ird_aver
      integer :: iav_run
      integer :: eirene_extrap
      integer :: eirene_uub_style
      integer :: eirene_background
      integer :: eirene_print_minmax
      integer :: save_f30
      integer :: save_f31
      integer :: ionising_core
      integer :: ank_mods

      real (kind=R8) :: eir_te_min
      real (kind=R8) :: eir_te_max
      real (kind=R8) :: eir_ti_min
      real (kind=R8) :: eir_ti_max
      real (kind=R8) :: eir_na_min
      real (kind=R8) :: eir_na_max
      real (kind=R8) :: eir_ne_min
      real (kind=R8) :: eir_ne_max
      real (kind=R8) :: eir_ua_min
      real (kind=R8) :: eir_ua_max

c k-enstrophy model
      integer :: transport_keps
      integer :: solve_keps
      integer :: keps_local
      integer :: keps_iout
      real(kind=R8) :: keps_cd
      real(kind=R8) :: keps_heat
      real(kind=R8) :: dna_min, vsa_min, hce_min, hci_min
      real(kind=R8) :: keps_heat_i
      real(kind=R8) :: keps_sig
      real(kind=R8) :: keps_alf
      real(kind=R8) :: keps_visc
      real(kind=R8) :: keps_dkt
      real(kind=R8) :: keps_dzt
      real(kind=R8) :: keps_inc
      real(kind=R8) :: keps_fac
      real(kind=R8) :: keps_shear
      integer :: b2sikt_iout
      integer :: b2sikt_style
      integer :: b2sikt_model
      integer :: b2sikt_local
      integer :: b2sikt_sheath_local
      integer :: b2sikt_fac_diss_core_mode
      integer :: b2sikt_min_source
      integer :: b2sikt_kt_source_stab
      real(kind=R8) :: b2sikt_fac_aniso
      real(kind=R8) :: b2sikt_fac_sheath
      real(kind=R8) :: b2sikt_fac_sheath_core
      real(kind=R8) :: b2sikt_fac_diss
      real(kind=R8) :: b2sikt_fac_diss_core
      real(kind=R8) :: b2sikt_fac_vis_RS
      real (kind=R8) :: b2tfhi_fflokt
      real (kind=R8) :: b2tfhi_fconkt
      real (kind=R8) :: b2tfhi_fflozt
      real (kind=R8) :: b2tfhi_fconzt
      real (kind=R8) :: b2tfhi_fsigkt
      real (kind=R8) :: b2tfhi_fkt_hie
      real (kind=R8) :: b2tfhe_vis_kt

c b2tqna
      real (kind=R8) :: b2tqna_ballooning
      real (kind=R8) :: b2tqna_bb_ref
      real (kind=R8) :: b2tqna_ballooning_rescale
      real (kind=R8) :: b2tqna_pfr_rescale
      real (kind=R8) :: b2tqna_divsol_rescale
      real (kind=R8) :: b2tqna_min_df0
      real (kind=R8) :: b2tqna_max_df0
      real (kind=R8) :: b2tqna_cfvma
      integer :: b2tqna_new_df0
      integer :: b2tqna_ixref
      integer :: b2tqna_transport_inputfile
      integer :: b2tqna_ballooning_sig
      integer :: b2tqna_limit_coeff

c b2tanml
      real (kind=R8) :: fch_anomalous

c b2stbm
      integer :: b2stbm_impgyro_mod
      integer :: b2stbm_internal_energy_sources
      real (kind=R8) :: b2stbm_linearisation

c b2stbr
      integer :: b2stbr_potential_at_guard_cell
      integer :: b2stbr_temperature_at_guard_cell
      integer :: b2stbr_sput_frac_flag
      integer :: b2stbr_b2wall_netcdf
      integer :: b2stbr_coreregno
      integer :: b2stbr_first_flight
      real (kind=R8) :: b2stbr_core_sources_rescale
      real (kind=R8) :: b2stbr_sput_frc
      real (kind=R8) :: b2stbr_alpha
      real (kind=R8) :: b2stbr_redep_alpha
      real (kind=R8) :: b2stbr_sput_mixed_beta
      real (kind=R8) :: b2stbr_sput_mixed_alpha
      real (kind=R8) :: b2stbr_sput_chem_cutoff_alpha
      real (kind=R8) :: b2stbr_sput_chem_cutoff_beta
      real (kind=R8) :: b2stbr_therm_evap
      real (kind=R8) :: b2stbr_sput_res
      real (kind=R8) :: b2stbr_sput_phys
      real (kind=R8) :: b2stbr_sput_chem_alpha
      real (kind=R8) :: recycled_neutrals_contr       !srv 12.04.08
      real (kind=R8) :: sput_phys_alpha
      real (kind=R8) :: b2stbr_in_dl

c b2sihs_
      real (kind=R8) :: b2sihs__rf0
      real (kind=R8) :: b2sihs__rf1
      real (kind=R8) :: b2sihs__rf2
      real (kind=R8) :: b2sihs__rf3
      real (kind=R8) :: b2sihs__rf4

c b2news_
      integer :: b2news_ncallout

c b2srst
      real (kind=R8) :: b2srst_rf0
      real (kind=R8) :: b2srst_rf1
      real (kind=R8) :: b2srst_rf2
      real (kind=R8) :: b2srst_rf3

c b2tfnb
      integer :: b2tfnb_no_hybr
      real (kind=R8) :: b2tfnb_flux_limit_min_ti

c b2srsm
      integer :: b2srsm_enable

c b2mndt
      integer :: density_control
      integer :: nstg(0:2)
      integer :: b2mndt_moqtlv
      integer :: b2mndt_moitlv
      real (kind=R8) :: nstg_aresco
      real (kind=R8) :: nstg_areshe
      real (kind=R8) :: nstg_areshi
      real (kind=R8) :: BoRiS
      real (kind=R8) :: b2mndt_rxf
      real (kind=R8) :: facvis_start
      real (kind=R8) :: facdrift_start
      real (kind=R8) :: facExB_start

c b2mndr
      real (kind=R8) :: b2mndr_hz
      real (kind=R8) :: b2mndr_stim
      real (kind=R8) :: neutral_rescale

c b2ag
      real (kind=R8) :: b2agmt_width
      real (kind=R8) :: geom_match_dist

c b2siav
      integer :: b2siav_qstyle
      real (kind=R8) :: addvis
      real (kind=R8) :: addvis1
      real (kind=R8) :: cqip1

c b2trcl
      real (kind=R8) :: cvsa_mltpl

c b2stbc
      integer :: bcene_16_style
      integer :: bceni_16_style

      real (kind=R8) :: b2stbc_cbc
      real (kind=R8) :: b2stbc_bc_ref
      real (kind=R8) :: b2stbc_bc_ref_te
      real (kind=R8) :: b2stbc_bc_ref_ti
      real (kind=R8) :: b2stbc_corr_flux
      real (kind=R8) :: bc_type13_norm
      real (kind=R8) :: bc_type13_fac
      real (kind=R8) :: nesepm
      real (kind=R8) :: ndes_sol
      real (kind=R8) :: nesepm_overshoot

c b2trno
      real (kind=R8) :: b2trno_alpha_stoch
c Feedback
      integer :: b2stbc_feedback

c Med
      integer :: med_style

      end type switches

      contains

***********************************************************************

      subroutine alloc_switches (s)

      implicit none
      type(switches), intent(in out) :: s

      !! TODO

      return
      end subroutine alloc_switches

***********************************************************************

      subroutine set_defaults_switches (s,ns)
      use b2mod_types
      implicit none
      integer, intent(in) :: ns
      type(switches), intent(in out) :: s

      s % get_residuals = 0
      s % mdf_fhe = 0
      s % mdf_fhi = 0
      s % mdf_fnb = 0
      s % pflux_cor = 1
      s % pot_eq = 1
      s % no_current = 0
      s % sheath_drift_fix = 1
      s % secmodel = 0
      s % b2stbc_fchy_dia_coreonly = 1
      s % nsmin = 0
      s % nsmax = ns
      s % icase_ii = 0
      s % icase_ee = 0
      s % icase_ei = 0
      s % stylec = 0
      s % zero_dead_and_core = 1
      s % fke_Zhdanov = 1
      s % conductive_limit = 1
      s % core_cond_limit = 0
      s % lluciani = 3
      s % no_solve = 0
      s % no_b2sral_call = 1
      s % do_2nd_b2npco_call = 0
      s % re_eval_prtls_fluxes = 0
      s % istyle_cur_contr_on_S_and_N = 2
      s % use_astra = 0
      s % lin_shi0 = 0
      s % flux_limit_to_dpa = 1
      s % flux_limit_to_vsa = 1
      s % csig_an_style = 1
      s % user_transport = 0

c Advanced fluid neutral models
      s % tn_style = 0
      s % recycle_afn = 0
      s % use_uy_uz_0 = 0
      s % use_uy_uz_1 = 0
      s % remove_FC_el = 0
      s % spatial_hybrid = 0
      s % use_auto_spatial_hyb = 0
      s % diffusion_bc_safeguard = 1
      s % do_volrec = 1
      s % transport_afn = 0

      s % fhe_vis_per = 0.0e0_R8
      s % fhe_vis_par = 0.0e0_R8
      s % fnb_vis_per = 0.0e0_R8
      s % fnb_vis_q   = 1.0e0_R8
      s % fhe_vis_q   = 1.0e0_R8
      s % PSch        = 1.0e0_R8
      s % fhePSch        = 1.0e0_R8
      s % fhiPSch        = 1.0e0_R8
      s % b2tfnb_xcur        = 0.0e0_R8
      s % b2tfnb_ycur        = 1.0e0_R8
      s % fch_ion_neutral = 0.0e0_R8
      s % art_rad = 0.0e0_R8
      s % sna0ep = 1.0e-36_R8
      s % she0ep = 1.0e-36_R8
      s % shi0ep = 1.0e-36_R8
      s % b2stbc_neoclassical = 0.0e0_R8
      s % integral_current = 0.0e0_R8
      s % b2stbc_fchy_dia = 0.0e0_R8
      s % b2mndr_na_min = 1.0e4_R8
      s % b2mndr_ua_max = 1.0e6_R8
      s % b2mndr_na_new = 1.0e14_R8
      s % cfc0 = 1.0e0_R8
      s % prl_cur = 1.0e0_R8
      s % dia_cur = 1.0e0_R8
      s % fch_pTe = 1.0e0_R8
      s % alfTeEh = 0.0_R8
      s % fhe_vdia_par = 1.0e0_R8
      s % fhi_vdia_par = 1.0e0_R8
      s % xvecrb = 1.0e0_R8
      s % xwdia = 1.0e0_R8
      s % fnb_drift_hyb = 1.0e0_R8
      s % xfac = 1.0e0_R8
      s % fch_inert = 1.0e0_R8
      s % cthiv = 2.65_R8
      s % cthev = 0.0e0_R8
      s % b2trcl_lambda = -5.0e0_R8
      s % poleldr = 1.0_R8
      s % radeldr = 1.0_R8
      s % poldidr = 1.0_R8
      s % raddidr = 1.0_R8
      s % b2tlv0_alpha = 0.0_R8
      s % b2tlv0_gamma = 2.0_R8
      s % b2tlh0_alpha = 0.0_R8
      s % b2tlh0_gamma = 2.0_R8
      s % b2tlc0_alpha = 0.0_R8
      s % b2tlc0_gamma = 2.0_R8
      s % b2tfnb_alpha = 0.0_R8
      s % b2tfnb_gamma = 2.0_R8
      s % fch_stochastic = 0.0_R8
      s % Qalfmin        = 1.0e-3_R8
      s % delpo          = 3.1e0_R8
      s % Kn_b1         = 0.01e0_R8
      s % Kn_b2         = 0.10e0_R8
      s % pol_anom_scale = 1.0_R8
      s % auto_spatial_hyb_Kn_1 = 0.0_R8
      s % auto_spatial_hyb_Kn_2 = 100.0_R8
      s % L_macro_AFN = 0.10_R8

      s % iout_b2wdat = 0
      s % b2sian_iout = 0
      s % b2sicf_iout = 0
      s % b2sigp_iout = 0
      s % b2sihs_iout = 0
      s % b2srdt_iout = 0
      s % b2srsm_iout = 0
      s % b2srst_iout = 0
      s % b2stcx_iout = 0
      s % b2stel_iout = 0
      s % b2stbr_iout = 0
      s % b2sral_iout = 0
      s % b2stbr_output = 0
      s % b2npht_iout = 0
      s % b2npmo_iout = 0
      s % b2npco_iout = 0
      s % b2npp7_iout = 0
      s % b2nppo_iout = 0
!      s % b2nph9_iout = 0
      s % b2news_iout = 0
      s % b2nxfx_iout = 0
      s % b2nxfc_iout = 0
      s % b2usco_iout = 0
      s % b2uspo_iout = 0
      s % b2usht_iout = 0
      s % b2usmo_iout = 0
      s % b2tfhi_iout = 0
      s % b2tfhe_iout = 0
      s % b2tfrn_iout = 0
      s % b2tfch_iout = 0
      s % b2tlmv_iout = 0
      s % b2tqca_iout = 0
      s % b2tqce_iout = 0
      s % b2trcl_iout = 0
      s % b2tfnb_iout = 0
      s % b2tral_iout = 0
      s % b2tlv0_iout = 0
      s % b2tlh0_iout = 0
      s % b2tlc0_iout = 0
      s % b2tqna_iout = 0
      s % b2tvspr_iout = 0
      s % b2tvskt_iout = 0
      s % b2stbc_iout = 0
      s % b2trno_iout = 0
      s % b2stbr_phys_iout = 0
      s % b2wdat_append = 0
      s % AFN_out = 0
      s % copy_background_iout =  0 !srv 14.07.10
      s % eirene_mc_iout = 0
      s % b2siav_iout = 0
      s % b2mndr_iout = 0

      s % b2sicf_style = 0
      s % b2sigp_style = 2
      s % b2sigp_pressure_restriction = 1
      s % b2sihs_style = 0
      s % b2sral_style = 2
      s % style_visc = 0
      s % b2tfnb_style_int_vel = 0
      s % b2sihs_style_int_vel = 0
      s % istyle_Joule_heating = 1
      s % b2sqcx_styl0 = 0
      s % b2stcx_styl0 = 0
      s % b2stel_styl0 = 0
      s % b2nxdv_style = 1
      s % b2nxfx_style = 1
      s % b2nxfc_style = 1
      s % mode_ehy = 0
      s % mode_ehx = 0
      s % b2tlmv_style = 1
      s % b2tqca_model = 1
      s % b2tqce_model = 1
      s % b2tqce_style_guard_cells = 0
      s % b2tral_mode = 1
      s % b2tlv0_style = 1
      s % b2tlh0_style = 1
      s % b2tlc0_style = 1
      s % b2tinnt_style = 1
      s % model_sig = 0
      s % b2mwqt_style = 1
      s % b2npht_style = 1
      s % fnb_drift_style = 1
      s % b2sifr_styl0 = 0

      s % b2sian_phm0 = 1.0e0_R8
      s % b2sicf_phm0 = 1.0e0_R8
      s % b2sicf_phm1 = 1.0e0_R8
      s % b2sigp_phm0 = 1.0e0_R8
      s % b2sihs_phm0 = 1.0e0_R8
      s % b2sihs_phm1 = 1.0e0_R8
      s % b2sihs_phm2 = 1.0e0_R8
      s % b2sihs_phm3 = 1.0e0_R8
      s % b2sihs_phm4 = 1.0e0_R8
      s % b2sihs_phm5 = 1.0e0_R8
      s % b2sihs_phm6 = 1.0e0_R8
      s % b2sihs_phm7 = 0.0e0_R8
      s % b2sihs_phm8 = 1.0e0_R8
      s % b2sqcx_phm0 = 1.0e0_R8
      s % b2sqel_phm0 = 1.0e0_R8
      s % b2sqel_phm1 = 1.0e0_R8
      s % b2sqel_phm2 = 1.0e0_R8
      s % b2stel_phm0 = 1.0e0_R8
      s % b2srdt_phm0 = 1.0e0_R8
      s % b2srdt_phm1 = 1.0e0_R8
      s % b2srdt_phm3 = 1.0e0_R8
      s % b2srdt_phm4 = 0.0e0_R8
      s % b2srdt_phm5 = 1.0e0_R8
      s % b2nxfv_phm0 = 1.0e0_R8
      s % b2nxfv_phm1 = 1.0e0_R8
      s % b2sifr_limthii = 0.3e0_R8
      s % b2sifr_limthee = 0.3e0_R8
      s % b2sifr_phm0 = 1.0e0_R8
      s % b2sifr_phm1 = 1.0e0_R8
      s % b2sifr_phm2 = 1.0e0_R8
      s % b2sifr_phm3 = 1.0e0_R8
      s % b2tqca_phm0 = 1.0e0_R8
      s % b2treq_phm0 = 1.0e0_R8
      s % b2stbc_phm0 = 1.0e0_R8
      s % b2stbc_phm1 = 1.0e0_R8
      s % b2stbc_phm2 = 1.0e0_R8

      s% b2npco_pcm0 = 1.0e0_R8
      s% b2npco_pcm1 = 1.0e0_R8
      s% b2npht_pcm0 = 1.0e0_R8
      s% b2npht_pcm1 = 1.0e0_R8

      s % b2stcx_rg0  = 1.0e0_R8
      s % b2npco_rxg =  1.0e0_R8
      s % b2npht_rxg =  1.0e0_R8
      s % b2npmo_rxg =  1.0e6_R8 !! 10^6 !!
      s % xfm0       =  1.0e0_R8
      s % xfm1       =  1.0e0_R8
      s % xfm2       =  1.0e0_R8
      s % xfm3       =  1.0e0_R8
      s % b2stel_rxm0=  0.0e0_R8
      s % b2stel_rxm1=  0.0e0_R8
      s % b2stel_rxm2=  0.0e0_R8
      s % b2stel_rxm3=  0.0e0_R8
      s % b2stel_rxm4=  0.0e0_R8
      s % b2stel_rg0 =  1.0e0_R8
      s % b2stel_rg1 =  1.0e0_R8

      s % b2npmo_discr_meth = 2
      s % b2tfhe_discr_meth = 2
      s % b2tfhi_discr_meth = 2
      s % b2tfnb_discr_meth = 2

      s % b2stbc_boundary_namelist = 0
      s % b2stbr_neutrals_namelist = 0
      s % b2tqna_transport_namelist = 0
      s % b2srdt_numerics_namelist = 0
      s % b2optim_namelist = 0

      s % b2stbc_diagno = 0
      s % b2npmo_diagno = 0
      s % ion_vlct_restrict = 0
      s % ion_vlct_restrict_M = 3.0_R8
      s % vlct_diagno = 1
      s % b2tqna_diagno = 0

      s % set_na_numerical = 0
      s % set_ua_numerical = 0
      s % set_te_numerical = 0
      s % set_ti_numerical = 0
      s % set_po_numerical = 0

c geometry
      s % fcHz_style = 1
      s % vxHz_style = 1
      s % use_9pt_stencil = 1
      s % fcVol_style = 0
      s % vxVol_style = 2

c eirene coupling
      s % use_eirene = 0
      s % eirene_sheath_pot = 1
      s % ird_aver = 0
      s % iav_run = 0
      s % eirene_extrap = 0
      s % eirene_uub_style = 0
      s % eirene_background = 1
      s % eirene_print_minmax = 0
      s % save_f30 = 0
      s % save_f31 = 0
      s % ionising_core = 0
      s % ank_mods = 0

      s % eir_te_min = 0.0e0_R8
      s % eir_te_max = 1.0e30_R8
      s % eir_ti_min = 0.0e0_R8
      s % eir_ti_max = 1.0e30_R8
      s % eir_na_min = 0.0e0_R8
      s % eir_na_max = 1.0e30_R8
      s % eir_ne_min = 0.0e0_R8
      s % eir_ne_max = 1.0e30_R8
      s % eir_ua_min = -c
      s % eir_ua_max = c
      s % neutral_sources_rescale = 1.0_R8

c k-enstrophy model
      s % transport_keps = 0
      s % solve_keps = 0
      s % keps_local = 1
      s % keps_iout = 0
      s % keps_cd = 1.0_R8
      s % keps_heat = 2.0_R8
      s % keps_heat_i = 2.0_R8
      s % keps_sig = 1.0e-4_R8
      s % keps_alf = 1.0e-4_R8
      s % dna_min = 1.0e-2_R8
      s % vsa_min = 1.0e-2_R8
      s % hce_min = 1.0e-2_R8
      s % hci_min = 1.0e-2_R8
      s % keps_visc = 0.2_R8
      s % keps_dkt = 0.01_R8
      s % keps_dzt = 0.01_R8
      s % keps_inc = 1.0_R8
      s % keps_fac = 1.0_R8
      s % keps_shear = 1.0e1_R8
      s % b2sikt_iout = 0
      s % b2sikt_style = 0
      s % b2sikt_model = 1
      s % b2sikt_local = 1
      s % b2sikt_fac_diss_core_mode = 0
      s % b2sikt_sheath_local = 1
      s % b2sikt_min_source = 0
      s % b2sikt_kt_source_stab = 1
      s % b2sikt_fac_aniso = 1.0_R8
      s % b2sikt_fac_sheath = 1.0_R8
      s % b2sikt_fac_sheath_core = 0.0_R8
      s % b2sikt_fac_diss = 0.0_R8
      s % b2sikt_fac_diss_core = 0.0_R8
      s % b2sikt_fac_vis_RS = 0.0_R8
      s % b2tfhi_fflokt = 0.0_R8
      s % b2tfhi_fconkt = 0.0_R8
      s % b2tfhi_fflozt = 0.0_R8
      s % b2tfhi_fconzt = 0.0_R8
      s % b2tfhi_fsigkt = 0.0_R8
      s % b2tfhi_fkt_hie = 0.0_R8
      s % b2tfhe_vis_kt = 0.0e0_R8

c b2tqna
      s % b2tqna_ballooning = 0.0_R8
      s % b2tqna_bb_ref = 0.0_R8
      s % b2tqna_ballooning_rescale = 1.0_R8
      s % b2tqna_pfr_rescale = 1.0_R8
      s % b2tqna_divsol_rescale = 1.0_R8
      s % b2tqna_min_df0 = 0.0_R8
      s % b2tqna_max_df0 = 1.0e30_R8
      s % b2tqna_cfvma = 0.0e0_R8                   !srv 15.12.05
      s % b2tqna_new_df0 = 0
      s % b2tqna_ixref = 1
      s % b2tqna_transport_inputfile = 0
      s % b2tqna_ballooning_sig = 1
      s % b2tqna_limit_coeff = 0

c fch_anomalous
      s % fch_anomalous = 1.0e0_R8

c b2stbm
      s % b2stbm_impgyro_mod = 0
      s % b2stbm_internal_energy_sources = 0
      s % b2stbm_linearisation = 1.0_R8

c b2stbr
      s % b2stbr_potential_at_guard_cell = 1
      s % b2stbr_temperature_at_guard_cell = 1
      s % b2stbr_sput_frac_flag = 0
      s % b2stbr_b2wall_netcdf = 0
      s % b2stbr_coreregno = 1
      s % b2stbr_first_flight = 0
      s % b2stbr_core_sources_rescale = 1.0_R8
      s % b2stbr_sput_frc = 0.0_R8
      s % b2stbr_alpha = 0.25_R8
      s % b2stbr_redep_alpha = 0.0_R8
      s % b2stbr_sput_mixed_beta = 1.0_R8
      s % b2stbr_sput_mixed_alpha = 1.0_R8
      s % b2stbr_sput_chem_cutoff_alpha = 1.0_R8
      s % b2stbr_sput_chem_cutoff_beta = 3.0_R8
      s % b2stbr_therm_evap = 0.0_R8
      s % b2stbr_sput_res = 0.0_R8
      s % b2stbr_sput_phys = 0.0_R8
      s % b2stbr_sput_chem_alpha = 1.0e-2_R8
      s % recycled_neutrals_contr = 1.0_R8
      s % sput_phys_alpha = 1.0e-2_R8
      s % b2stbr_in_dl = 0.001_R8

c b2sihs_
      s % b2sihs__rf0 = 1.0_R8
      s % b2sihs__rf1 = 1.0_R8
      s % b2sihs__rf2 = 1.0_R8
      s % b2sihs__rf3 = 1.0_R8
      s % b2sihs__rf4 = 1.0_R8

c b2news_
      s % b2news_ncallout = -1

c b2sihs_
      s % b2srst_rf0 = 1.0_R8
      s % b2srst_rf1 = 1.0_R8
      s % b2srst_rf2 = 1.0_R8
      s % b2srst_rf3 = 1.0_R8

c b2tfnb
      s % b2tfnb_no_hybr = 0
      s % b2tfnb_flux_limit_min_ti = 0.0_R8

c b2srsm
      s % b2srsm_enable = 0

c b2mndt
      s % density_control = 0
      s % nstg(0) = 1
      s % nstg(1) = 1
      s % nstg(2) = 1
      s % b2mndt_moqtlv = 3
      s % b2mndt_moitlv = -1
      s % nstg_aresco = 0.0_R8
      s % nstg_areshe = 0.0_R8
      s % nstg_areshi = 0.0_R8
      s % BoRiS = 0.0_R8
      s % b2mndt_rxf = 0.5_R8
      s % facvis_start = 0.0_R8
      s % facdrift_start = 0.0_R8
      s % facExB_start = 0.0_R8

c Feedback
      s % b2stbc_feedback = 0

c b2mndr
      s % b2mndr_hz = 1.0_R8
      s % b2mndr_stim = 0.0_R8
      s % neutral_rescale = 1.0_R8

c b2ag
      s % b2agmt_width = 1.0_R8
      s % geom_match_dist = 1.0e-6_R8

c b2siav
      s % b2siav_qstyle = 0
      s % addvis = 0.0_R8
      s % addvis1 = 1.0_R8
      s % cqip1 = 0.5_R8
      s % cvsa_mltpl = 1.0_R8

c b2stbc_phys
      s % bcene_16_style = 0
      s % bceni_16_style = 0

      s % b2stbc_cbc = 1.0e0_R8              !srv 14.10.99
      s % b2stbc_bc_ref = 0.01_R8            !srv 14.10.99
      s % b2stbc_bc_ref_te = 0.01         !IYS 07.03.2013 11.04.13
      s % b2stbc_bc_ref_ti=0.01         !IYS 07.03.2013 11.04.13
      s % b2stbc_corr_flux   = 0.0e0_R8                    !srv 25.07.05
      s % bc_type13_norm = 1.0e15_R8
      s % bc_type13_fac = 1.0_R8
      s % nesepm = 0.0_R8
      s % nesepm_overshoot = 0.0_R8
      s % ndes_sol = 0.0_R8

c b2trno
      s % b2trno_alpha_stoch =  1.0e0_R8                  !srv 17.12.13

c med
      s % med_style = 0

      return
      end subroutine set_defaults_switches

***********************************************************************

      subroutine dealloc_switches (s)

      implicit none
      type(switches), intent(in out) :: s

      !! TODO

      return
      end subroutine dealloc_switches

***********************************************************************


      subroutine read_switches (s)

      implicit none
      type (switches), intent (inout) :: s

      call alloc_switches (s)

      call ipgeti ('get_residuals',s % get_residuals)
      call ipgeti ('b2tfhe_mdf', s % mdf_fhe)
      call ipgeti ('b2tfhi_mdf', s % mdf_fhi)
      call ipgeti ('b2tfnb_mdf', s % mdf_fnb)
      call ipgeti ('b2tfnb_pflux_cor', s % pflux_cor)
      call ipgeti ('b2news_poteq', s % pot_eq)
      call ipgeti ('b2tfhe_no_current',s % no_current)
      call ipgeti ('b2stbc_sheath_drift_fix', s % sheath_drift_fix)
      call ipgeti ('b2stbc_secmodel', s % secmodel)
      call ipgeti ('b2stbc_fchy_dia_coreonly',
     .               s % b2stbc_fchy_dia_coreonly)
      call ipgeti ('b2news_nsmin', s % nsmin)
      call ipgeti ('b2news_nsmax', s % nsmax)
      call ipgeti ('b2tlnl_ii', s %  icase_ii)
      call ipgeti ('b2tlnl_ee', s %  icase_ee)
      call ipgeti ('b2tlnl_ei', s %  icase_ei)
      call ipgeti ('b2upht_stylec', s % stylec)
      call ipgeti('b2xzdd_zero_dead_and_core', s % zero_dead_and_core)
      call ipgeti ('b2tqce_fke_Zhdanov', s % fke_Zhdanov)                  !srv 02.12.10
      call ipgeti ('b2trcl_conductive_limit', s % conductive_limit)        !srv 16.06.08
      call ipgeti ('b2trcl_core_cond_limit', s % core_cond_limit)          !csw 13mar2012
      call ipgeti ('b2trcl_lluciani', s % lluciani)
      call ipgeti ('b2trcl_core_cond_limit', s % core_cond_limit)          !csw 13mar2012
      call ipgeti ('b2news_no_solve', s % no_solve)                        !srv 02.07.08
      call ipgeti ('b2news_no_b2sral_call', s%no_b2sral_call)
      call ipgeti ('b2news_do_2nd_b2npco_call', s%do_2nd_b2npco_call)
      call ipgeti ('b2news_re_eval_prtls_fluxes',s%re_eval_prtls_fluxes)   !srv 28.06.13
      call ipgeti ('b2stbc_istyle_cur_contr_on_S_and_N',
     &                       s%istyle_cur_contr_on_S_and_N)                !srv 01.02.09
      call ipgeti ('b2mndr_astra', s%use_astra)
      call ipgeti ('b2stbr_phys_lin_shi0', s % lin_shi0)
      call ipgeti ('b2trno_flux_limit_to_dpa', s % flux_limit_to_dpa)
      call ipgeti ('b2trno_flux_limit_to_vsa', s % flux_limit_to_vsa)
      call ipgeti ('b2trno_csig_an_style', s % csig_an_style)               !srv 15.06.08
      call ipgeti ('b2tqna_user_transport', s % user_transport)             !mb  11.10.17

c Advanced fluid neutral models
      call ipgeti ('b2mn_tn_style',s % tn_style)
      call ipgeti ('b2stbr_recycle_afn', s % recycle_afn)
      call ipgeti ('b2stbr_use_uy_uz_0', s % use_uy_uz_0)
      call ipgeti ('b2stbr_use_uy_uz_1', s % use_uy_uz_1)
      call ipgeti ('b2stbr_remove_FC_el', s % remove_FC_el)
      call ipgeti ('b2mn_spatial_hybrid', s % spatial_hybrid)
      call ipgeti ('use_auto_spatial_hyb', s % use_auto_spatial_hyb)
      call ipgeti ('b2mn_do_volrec', s % do_volrec)
      call ipgeti ('b2tqna_transport_afn', s % transport_afn)
      call ipgeti ('diffusion_bc_safeguard', s % diffusion_bc_safeguard)

      call ipgetr ('b2tfhe_vis_per', s % fhe_vis_per)
      call ipgetr ('b2tfhe_vis_par', s % fhe_vis_par)
      call ipgetr ('b2tfnb_vis_per', s % fnb_vis_per)
      call ipgetr ('b2tfnb_vis_q', s % fnb_vis_q)
      call ipgetr ('b2tfhe_vis_q', s % fhe_vis_q)
      call ipgetr ('b2tfnb_PSch', s % PSch)
      call ipgetr ('b2tfhe_PSch', s % fhePSch)
      call ipgetr ('b2tfhi_PSch', s % fhiPSch)
      call ipgetr ('b2tfnb_xcur', s % b2tfnb_xcur)
      call ipgetr ('b2tfnb_ycur', s % b2tfnb_ycur)
      call ipgetr ('b2tfhe_neutral', s % fch_ion_neutral)
      call ipgetr ('b2sqel_artificial_radiation', s % art_rad)
      call ipgetr ('b2stbc_sna0ep', s % sna0ep)
      call ipgetr ('b2stbc_shi0ep', s % shi0ep)
      call ipgetr ('b2stbc_she0ep', s % she0ep)
      call ipgetr ('b2stbc_neoclassical', s % b2stbc_neoclassical)          !xpb
      call ipgetr ('b2stbc_integral_current', s % integral_current)         !dpc
      call ipgetr ('b2stbc_fchy_dia', s % b2stbc_fchy_dia)                  !xpb
      call ipgetr ('b2mndr_na_min', s%b2mndr_na_min)
      call ipgetr ('b2mndr_ua_max', s%b2mndr_ua_max)
      call ipgetr ('b2mndr_na_new', s%b2mndr_na_new)
      call ipgetr ('b2usmo_cfc0', s % cfc0)
      call ipgetr ('b2tfhe_prl_cur', s % prl_cur)
      call ipgetr ('b2tfhe_dia_cur', s % dia_cur)                           !srv 17.05.99
      call ipgetr ('b2tfhe_fch_pTe', s % fch_pTe)                           !srv 30.12.09
      call ipgetr ('b2tfhe_alfTeEh', s % alfTeEh)                           !srv 28.04.08
      call ipgetr ('b2tfhe_vdia_par', s % fhe_vdia_par)                     !srv 17.05.99
      call ipgetr ('b2tfhi_vdia_par', s % fhi_vdia_par)                     !srv 17.05.99
      call ipgetr ('b2tfrn_xvecrb', s % xvecrb)                             !srv 27.01.12
      call ipgetr ('b2tfrn_xwdia', s % xwdia)                               !srv 27.01.12
      call ipgetr ('b2tfnb_drift_hyb',s % fnb_drift_hyb)
      call ipgetr ('b2tfcc_xfac', s % xfac)
      call ipgetr ('b2tiner_inert', s % fch_inert)                          !srv 07.07.99
      call ipgetr ('b2trcl_cthi', s % cthiv)
      call ipgetr ('b2trcl_cthe', s % cthev)
      call ipgetr ('b2trcl_lambda', s % b2trcl_lambda)
      call ipgetr ('b2tfnb_poleldr', s % poleldr)                         !srv 25.04.00 {
      call ipgetr ('b2tfnb_radeldr', s % radeldr)
      call ipgetr ('b2tfnb_poldidr', s % poldidr)
      call ipgetr ('b2tfnb_raddidr', s % raddidr)
      call ipgetr ('b2tlv0_alpha', s % b2tlv0_alpha)
      call ipgetr ('b2tlv0_gamma', s % b2tlv0_gamma)
      call ipgetr ('b2tlh0_alpha', s % b2tlh0_alpha)
      call ipgetr ('b2tlh0_gamma', s % b2tlh0_gamma)
      call ipgetr ('b2tlc0_alpha', s % b2tlc0_alpha)
      call ipgetr ('b2tlc0_gamma', s % b2tlc0_gamma)
      call ipgetr ('b2tfnb_alpha', s % b2tfnb_alpha)
      call ipgetr ('b2tfnb_gamma', s % b2tfnb_gamma)
      call ipgetr ('b2tfhe_stochastic', s%fch_stochastic)           !lk 31.01.08
      call ipgetr ('b2stbc_Qalfmin', s % Qalfmin)
      call ipgetr ('b2stbc_delpo', s % delpo)
      call ipgetr ('b2stbr_Kn_b1',s % Kn_b1)
      call ipgetr ('b2stbr_Kn_b2',s % Kn_b2)
      call ipgetr ('b2trno_pol_anom_scale', s % pol_anom_scale)
      call ipgetr ('auto_spatial_hyb_Kn_1', s % auto_spatial_hyb_Kn_1)
      call ipgetr ('auto_spatial_hyb_Kn_2', s % auto_spatial_hyb_Kn_2)
      call ipgetr ('L_macro_AFN', s % L_macro_AFN)

c iout switches
      call ipgeti ('b2wdat_iout', s % iout_b2wdat)
      call ipgeti ('b2sian_iout', s % b2sian_iout)
      call ipgeti ('b2sicf_iout', s % b2sicf_iout)
      call ipgeti ('b2sigp_iout', s % b2sigp_iout)
      call ipgeti ('b2sihs_iout', s % b2sihs_iout)
      call ipgeti ('b2srdt_iout', s % b2srdt_iout)
      call ipgeti ('b2srsm_iout', s % b2srsm_iout)
      call ipgeti ('b2srst_iout', s % b2srst_iout)
      call ipgeti ('b2stcx_iout', s % b2stcx_iout)
      call ipgeti ('b2stel_iout', s % b2stel_iout)
      call ipgeti ('b2stbr_iout', s % b2stbr_iout)
      call ipgeti ('b2sral_iout', s % b2sral_iout)                 !srv 30.06.08
      call ipgeti ('b2stbr_output', s % b2stbr_output)
      call ipgeti ('b2npht_iout', s % b2npht_iout)
      call ipgeti ('b2npmo_iout', s % b2npmo_iout)
      call ipgeti ('b2npco_iout', s % b2npco_iout)
      call ipgeti ('b2npp7_iout', s % b2npp7_iout)
      call ipgeti ('b2nppo_iout', s % b2nppo_iout)
!      call ipgeti ('b2nph9_iout', s % b2nph9_iout)
      call ipgeti ('b2news_iout', s % b2news_iout)
      call ipgeti ('b2nxfx_iout', s % b2nxfx_iout)
      call ipgeti ('b2nxfc_iout', s % b2nxfc_iout)
      call ipgeti ('b2usco_iout', s % b2usco_iout)
      call ipgeti ('b2uspo_iout', s % b2uspo_iout)
      call ipgeti ('b2usht_iout', s % b2usht_iout)
      call ipgeti ('b2usmo_iout', s % b2usmo_iout)
      call ipgeti ('b2tfhi_iout', s % b2tfhi_iout)
      call ipgeti ('b2tfhe_iout', s % b2tfhe_iout)
      call ipgeti ('b2tfrn_iout', s % b2tfrn_iout)
      call ipgeti ('b2tfch_iout', s % b2tfch_iout)
      call ipgeti ('b2tlmv_iout', s % b2tlmv_iout)
      call ipgeti ('b2tqca_iout', s % b2tqca_iout)                            !srv 16.06.08
      call ipgeti ('b2tqce_iout', s % b2tqce_iout)                            !srv 16.06.08
      call ipgeti ('b2trcl_iout', s % b2trcl_iout)                            !srv 17.06.08
      call ipgeti ('b2tfnb_iout', s % b2tfnb_iout)                            !srv 16.09.02
      call ipgeti ('b2tral_iout', s % b2tral_iout)
      call ipgeti ('b2tlv0_iout', s % b2tlv0_iout)
      call ipgeti ('b2tlh0_iout', s % b2tlh0_iout)
      call ipgeti ('b2tlc0_iout', s % b2tlc0_iout)
      call ipgeti ('b2tqna_iout', s % b2tqna_iout)                            !srv 10.11.09
      call ipgeti ('b2tvspr_iout', s % b2tvspr_iout)
      call ipgeti ('b2tvskt_iout', s % b2tvskt_iout)
      call ipgeti ('b2stbc_iout', s % b2stbc_iout)                            !srv 11.07.05
      call ipgeti ('b2trno_iout', s % b2trno_iout)                            !srv 17.06.08
      call ipgeti ('b2stbr_phys_iout', s % b2stbr_phys_iout)                  !srv 11.07.05
      call ipgeti('b2wdat_append', s % b2wdat_append)
      call ipgeti ('b2stbr_AFN_out',s % AFN_out)
      call ipgeti('copy_background_iout', s%copy_background_iout)   !srv 14.07.10
      call ipgeti ('eirene_mc_iout', s % eirene_mc_iout)
      call ipgeti ('b2mndr_iout', s % b2mndr_iout)

c style/mode switches
      call ipgeti ('b2sicf_style', s % b2sicf_style)
      call ipgeti ('b2sigp_style', s % b2sigp_style)
      call ipgeti ('b2sigp_pressure_restriction',
     &  s % b2sigp_pressure_restriction)
      call ipgeti ('b2sihs_style', s % b2sihs_style)
      call ipgeti ('b2sral_style', s % b2sral_style)                       !srv 28.07.08
      call ipgeti ('b2sihs_style_visc', s % style_visc)
      call ipgeti ('b2tfnb_style_int_vel', s % b2tfnb_style_int_vel)
      call ipgeti ('b2sihs_style_int_vel', s % b2sihs_style_int_vel)
      call ipgeti ('b2sihs_istyle_Joule_heating',
     &                      s % istyle_Joule_heating)
      call ipgeti ('b2sqcx_styl0', s % b2sqcx_styl0)
      call ipgeti ('b2stcx_styl0', s % b2stcx_styl0)
      call ipgeti ('b2stel_styl0', s % b2stel_styl0)
      call ipgeti ('b2nxdv_style', s % b2nxdv_style)
      call ipgeti ('b2nxfx_style', s % b2nxfx_style)
      call ipgeti ('b2nxfc_style', s % b2nxfc_style)
      call ipgeti ('b2tfhe_mode_ehy',s% mode_ehy)                          !srv 06.04.15
      call ipgeti ('b2tfhe_mode_ehx',s% mode_ehx)                          !srv 06.04.15
      call ipgeti ('b2tlmv_style', s% b2tlmv_style)                        !srv 30.05.08
      call ipgeti ('b2tqca_model', s% b2tqca_model)
      call ipgeti ('b2tqce_model', s% b2tqce_model)
      call ipgeti ('b2tqce_style_guard_cells',
     &                                    s%b2tqce_style_guard_cells)      !srv 19.04.13
      call ipgeti ('b2tral_mode',s% b2tral_mode)
      call ipgeti ('b2tlv0_style', s % b2tlv0_style)
      call ipgeti ('b2tlh0_style', s % b2tlh0_style)
      call ipgeti ('b2tlc0_style', s % b2tlc0_style)
      call ipgeti ('b2tinnt_style', s % b2tinnt_style)
      call ipgeti ('b2tqna_model_sig', s % model_sig)               !srv 02.12.04
      call ipgeti ('b2mwqt_style', s % b2mwqt_style)
      call ipgeti ('b2npht_style', s % b2npht_style)                !srv 09.07.08
      call ipgeti ('b2tfnb_drift_style', s % fnb_drift_style)       !srv 21.06.08f
      call ipgeti ('b2sifr_styl0', s % b2sifr_styl0)

c phm switches
      call ipgetr ('b2sian_phm0', s % b2sian_phm0)
      call ipgetr ('b2sicf_phm0', s % b2sicf_phm0)
      call ipgetr ('b2sicf_phm1', s % b2sicf_phm1)
      call ipgetr ('b2sigp_phm0', s % b2sigp_phm0)
      call ipgetr ('b2sihs_phm0', s % b2sihs_phm0)
      call ipgetr ('b2sihs_phm1', s % b2sihs_phm1)
      call ipgetr ('b2sihs_phm2', s % b2sihs_phm2)
      call ipgetr ('b2sihs_phm3', s % b2sihs_phm3)
      call ipgetr ('b2sihs_phm4', s % b2sihs_phm4)
      call ipgetr ('b2sihs_phm5', s % b2sihs_phm5)
      call ipgetr ('b2sihs_phm6', s % b2sihs_phm6)
      call ipgetr ('b2sihs_phm7', s % b2sihs_phm7)
      call ipgetr ('b2sihs_phm8', s % b2sihs_phm8)
      call ipgetr ('b2sqcx_phm0', s % b2sqcx_phm0)
      call ipgetr ('b2sqel_phm0', s % b2sqel_phm0)
      call ipgetr ('b2sqel_phm1', s % b2sqel_phm1)
      call ipgetr ('b2sqel_phm2', s % b2sqel_phm2)
      call ipgetr ('b2stel_phm0', s % b2stel_phm0)
      call ipgetr ('b2srdt_phm0', s % b2srdt_phm0)
      call ipgetr ('b2srdt_phm1', s % b2srdt_phm1)
      call ipgetr ('b2srdt_phm3', s % b2srdt_phm3)
      call ipgetr ('b2srdt_phm4', s % b2srdt_phm4)
      call ipgetr ('b2srdt_phm5', s % b2srdt_phm5)
      call ipgetr ('b2nxfv_phm0', s % b2nxfv_phm0)
      call ipgetr ('b2nxfv_phm1', s % b2nxfv_phm1)
      call ipgetr ('b2sifr_limthii', s % b2sifr_limthii)
      call ipgetr ('b2sifr_limthee', s % b2sifr_limthee)
      call ipgetr ('b2sifr_phm0', s % b2sifr_phm0)
      call ipgetr ('b2sifr_phm1', s % b2sifr_phm1)
      call ipgetr ('b2sifr_phm2', s % b2sifr_phm2)
      call ipgetr ('b2sifr_phm3', s % b2sifr_phm3)
      call ipgetr ('b2tqca_phm0', s % b2tqca_phm0)
      call ipgetr ('b2treq_phm0', s % b2treq_phm0)
      call ipgetr ('b2stbc_phm0', s % b2stbc_phm0)  !srv 17.01.07
      call ipgetr ('b2stbc_phm1', s % b2stbc_phm1)
      call ipgetr ('b2stbc_phm2', s % b2stbc_phm2)


c pcm switches
      call ipgetr ('b2npco_pcm0', s% b2npco_pcm0)
      call ipgetr ('b2npco_pcm1', s% b2npco_pcm1)
      call ipgetr ('b2npht_pcm0', s% b2npht_pcm0)
      call ipgetr ('b2npht_pcm1', s% b2npht_pcm1)

c addition numerical relaxation
      call ipgetr ('b2stcx_rg0', s % b2stcx_rg0)
      call ipgetr ('b2npco_rxg', s % b2npco_rxg)
      call ipgetr ('b2npht_rxg', s % b2npht_rxg)
      call ipgetr ('b2npmo_rxg', s % b2npmo_rxg)
      call ipgetr ('b2news_xfm0', s % xfm0)
      call ipgetr ('b2news_xfm1', s % xfm1)
      call ipgetr ('b2news_xfm2', s % xfm2)
      call ipgetr ('b2news_xfm3', s % xfm3)
      call ipgetr ('b2stel_rxm0', s % b2stel_rxm0)
      call ipgetr ('b2stel_rxm1', s % b2stel_rxm1)
      call ipgetr ('b2stel_rxm2', s % b2stel_rxm2)
      call ipgetr ('b2stel_rxm3', s % b2stel_rxm3)
      call ipgetr ('b2stel_rxm4', s % b2stel_rxm4)
      call ipgetr ('b2stel_rg0', s % b2stel_rg0)
      call ipgetr ('b2stel_rg1', s % b2stel_rg1)

c discretization methods
      call ipgeti ('b2npmo_discr_meth', s % b2npmo_discr_meth)
      call ipgeti ('b2tfhe_discr_meth', s % b2tfhe_discr_meth)
      call ipgeti ('b2tfhi_discr_meth', s % b2tfhi_discr_meth)
      call ipgeti ('b2tfnb_discr_meth', s % b2tfnb_discr_meth)

c namelist switches
      call ipgeti('b2stbc_boundary_namelist',s%b2stbc_boundary_namelist)
      call ipgeti('b2stbr_neutrals_namelist',s%b2stbr_neutrals_namelist)
      call ipgeti('b2tqna_transport_namelist',
     &                                      s%b2tqna_transport_namelist)
      call ipgeti ('b2srdt_numerics_namelist',
     &  s%b2srdt_numerics_namelist)
      call ipgeti ('b2optim_namelist', s%b2optim_namelist)

c diagnostic switches
      call ipgeti ('b2stbc_diagno', s % b2stbc_diagno)                     !xpb
      call ipgeti ('b2npmo_diagno', s % b2npmo_diagno)
      call ipgeti ('b2npmo_ion_vlct_restrict', s % ion_vlct_restrict)
      call ipgetr ('b2npmo_ion_vlct_restrict_M', s%ion_vlct_restrict_M)
      call ipgeti ('b2npmo_vlct_diagno', s % vlct_diagno)
      call ipgeti ('b2tqna_diagno', s % b2tqna_diagno)

c MMS
      call ipgeti ('b2mndr_set_na_numerical', s % set_na_numerical)
      call ipgeti ('b2mndr_set_ua_numerical', s % set_ua_numerical)
      call ipgeti ('b2mndr_set_te_numerical', s % set_te_numerical)
      call ipgeti ('b2mndr_set_ti_numerical', s % set_ti_numerical)
      call ipgeti ('b2mndr_set_po_numerical', s % set_po_numerical)

c geometry
      call ipgeti ('b2mndr_fcHz_style',s % fcHz_style)
      call ipgeti ('b2mndr_vxHz_style',s % vxHz_style)
      call ipgeti ('b2mndr_use_9pt_stencil',s % use_9pt_stencil)
      call ipgeti ('b2mndr_fcVol_style',s % fcVol_style)
      call ipgeti ('b2mndr_vxVol_style', s % vxVol_style)

c Eirene coupling
      call ipgeti ('b2mndr_eirene', s % use_eirene)
      call ipgeti ('eirene_sheath_pot',s % eirene_sheath_pot)
      call ipgeti ('b2mndr_av_read', s % ird_aver)                !pb 29.09.16
      call ipgeti ('b2mndt_av', s % iav_run)
      call ipgeti('eirene_extrap', s % eirene_extrap)
      call ipgeti('eirene_uub_style', s % eirene_uub_style)
      call ipgeti('eirene_background', s % eirene_background)
      call ipgeti('eirene_print_minmax', s % eirene_print_minmax)
      call ipgeti ('eirene_savef30', s % save_f30)
      call ipgeti ('eirene_savef31', s % save_f31)                !srv 14.07.10
      call ipgeti ('eirene_ionising_core', s % ionising_core)
      if (s%ionising_core.eq.0)     ! for backward compatibility
     > call ipgeti ('eirene_ionizing_core', s % ionising_core)
      call ipgeti('eirene_ank_mods', s % ank_mods)

      call ipgetr('eirene_te_min',s % eir_te_min)
      call ipgetr('eirene_te_max',s % eir_te_max)
      call ipgetr('eirene_ti_min',s % eir_ti_min)
      call ipgetr('eirene_ti_max',s % eir_ti_max)
      call ipgetr('eirene_na_min',s % eir_na_min)
      call ipgetr('eirene_na_max',s % eir_na_max)
      call ipgetr('eirene_ne_min',s % eir_ne_min)
      call ipgetr('eirene_ne_max',s % eir_ne_max)
      call ipgetr('eirene_ua_min',s % eir_ua_min)
      call ipgetr('eirene_ua_max',s % eir_ua_max)
      call ipgetr ('b2mndr_rescale_neutrals_sources',
     &                s % neutral_sources_rescale)


c k-enstrophy model
      call ipgeti ('b2mndr_solve_keps',s % solve_keps)
      call ipgeti ('b2tqna_transport_keps', s % transport_keps)
      call ipgetr ('b2tqna_keps_cd', s % keps_cd)
      call ipgetr ('b2tqna_keps_visc', s % keps_visc)
      call ipgetr ('b2tqna_keps_dkt', s % keps_dkt)
      call ipgetr ('b2tqna_keps_dzt', s % keps_dzt)
      call ipgetr ('b2tqna_keps_heat', s % keps_heat)
      call ipgetr ('b2tqna_keps_heat_i', s % keps_heat_i)
      call ipgetr ('b2tqna_keps_sig', s % keps_sig)
      call ipgetr ('b2tqna_keps_alf', s % keps_alf)
      call ipgeti ('b2tqna_keps_local', s % keps_local)
      call ipgeti ('b2tqna_keps_iout', s % keps_iout)
      call ipgetr ('b2tqna_keps_dna_min', s % dna_min)
      call ipgetr ('b2tqna_keps_vsa_min', s % vsa_min)
      call ipgetr ('b2tqna_keps_hce_min', s % hce_min)
      call ipgetr ('b2tqna_keps_hci_min', s % hci_min)
      call ipgetr ('b2tqna_keps_inc', s % keps_inc)
      call ipgetr ('b2tqna_keps_init', s % keps_fac)
      call ipgetr ('b2tqna_keps_shear', s % keps_shear)
      call ipgeti ('b2sikt_iout', s % b2sikt_iout)
      call ipgeti ('b2sikt_style', s % b2sikt_style)
      call ipgeti ('b2sikt_model', s % b2sikt_model)
      call ipgeti ('b2sikt_keps_local', s % b2sikt_local)
      call ipgeti ('b2sikt_fac_sheath_local', s % b2sikt_sheath_local)
      call ipgeti ('b2sikt_min_source', s % b2sikt_min_source)
      call ipgeti ('b2sikt_kt_source_stab', s % b2sikt_kt_source_stab)
      call ipgetr ('b2sikt_aniso', s % b2sikt_fac_aniso)
      call ipgetr ('b2sikt_fac_sheath', s % b2sikt_fac_sheath)
      call ipgetr ('b2sikt_fac_sheath_core', s % b2sikt_fac_sheath_core)
      call ipgetr ('b2sikt_fac_diss', s % b2sikt_fac_diss)
      call ipgetr ('b2sikt_fac_diss_core', s % b2sikt_fac_diss_core)
      call ipgeti ('b2sikt_fac_diss_core_mode',
     .          s % b2sikt_fac_diss_core_mode)
      call ipgetr ('b2sikt_fac_vis_RS', s % b2sikt_fac_vis_RS)
      call ipgetr ('b2tfhi_fflokt', s % b2tfhi_fflokt)
      call ipgetr ('b2tfhi_fconkt', s % b2tfhi_fconkt)
      call ipgetr ('b2tfhi_fflozt', s % b2tfhi_fflozt)
      call ipgetr ('b2tfhi_fconzt', s % b2tfhi_fconzt)
      call ipgetr ('b2tfhi_fsigkt', s % b2tfhi_fsigkt)
      call ipgetr ('b2tfhi_fkt_hie', s % b2tfhi_fkt_hie)
      call ipgetr ('b2tfhe_vis_kt', s % b2tfhe_vis_kt)

c b2tqna
      call ipgeti ('b2tqna_new_df0', s%b2tqna_new_df0)
      call ipgetr ('b2tqna_min_df0', s%b2tqna_min_df0)
      call ipgetr ('b2tqna_max_df0', s%b2tqna_max_df0)
      call ipgeti ('b2tqna_ixref', s%b2tqna_ixref)
      call ipgeti ('b2tqna_inputfile', s%b2tqna_transport_inputfile)
      call ipgetr ('b2tqna_ballooning', s%b2tqna_ballooning)
      call ipgetr ('b2tqna_bb_ref', s%b2tqna_bb_ref)
      call ipgetr ('b2tqna_ballooning_rescale',
     & s%b2tqna_ballooning_rescale)
      call ipgeti ('b2tqna_ballooning_sig', s%b2tqna_ballooning_sig)
      call ipgeti ('b2tqna_limit_coeff', s%b2tqna_limit_coeff)
      call ipgetr ('b2tqna_pfr_rescale', s%b2tqna_pfr_rescale)
      call ipgetr ('b2tqna_divsol_rescale', s%b2tqna_divsol_rescale)
      call ipgetr ('b2tqna_cfvma', s%b2tqna_cfvma)                !srv 15.12.05

c fch anomalous
      call ipgetr ('b2tfhe_anomalous', s%fch_anomalous)
       if (s%fch_anomalous.eq.1.0_R8)
     &  call ipgetr ('b2tanml_anomalous', s%fch_anomalous) ! for backward compatibility

c b2stbm
      call ipgetr ('b2stbm_linearisation', s % b2stbm_linearisation)
      if (s%b2stbm_linearisation.eq.1.0_R8)   ! for backward compatibility
     &  call ipgetr ('b2stbm_linearization', s % b2stbm_linearisation)
      call ipgeti ('b2stbm_impgyro_mod', s % b2stbm_impgyro_mod)
      call ipgeti ('b2stbm_internal_energy_sources',
     &  s % b2stbm_internal_energy_sources)

c b2stbr
      call ipgeti ('b2stbr_potential_at_guard_cell',
     &  s % b2stbr_potential_at_guard_cell)
      call ipgeti ('b2stbr_temperature_at_guard_cell',
     &  s % b2stbr_temperature_at_guard_cell)
      call ipgeti ('b2stbc_coreregno', s % b2stbr_coreregno)
      call ipgeti ('b2stbr_first_flight', s % b2stbr_first_flight)
      call ipgeti ('b2stbr_sput_frac_flag', s % b2stbr_sput_frac_flag)
      call ipgeti ('b2stbr_b2wall_netcdf', s % b2stbr_b2wall_netcdf)
      call ipgetr ('b2stbr_core_sources_rescale',
     &  s % b2stbr_core_sources_rescale)
      call ipgetr ('b2stbr_sput_frc', s % b2stbr_sput_frc)
      call ipgetr ('b2stbr_alpha', s % b2stbr_alpha)
      if(s%b2stbr_sput_frac_flag.ne.0) then
        call ipgetr ('b2stbr_sput_mixed_alpha',
     &    s % b2stbr_sput_mixed_alpha)
        call ipgetr ('b2stbr_sput_mixed_beta',
     &    s % b2stbr_sput_mixed_beta)
      endif
      if(s%b2stbr_sput_frc.ne.0.0_R8) then
        call ipgetr ('b2stbr_sput_chem_alpha', s%b2stbr_sput_chem_alpha)
        call ipgetr ('b2stbr_sput_chem_cutoff_alpha',
     &   s%b2stbr_sput_chem_cutoff_alpha)
        call ipgetr ('b2stbr_sput_chem_cutoff_beta',
     &   s%b2stbr_sput_chem_cutoff_beta)
      endif
      call ipgetr ('b2stbr_therm_evap', s%b2stbr_therm_evap)
      call ipgetr ('b2stbr_sput_res', s%b2stbr_sput_res)
      call ipgetr ('b2stbr_sput_phys', s%b2stbr_sput_phys)
      call ipgetr ('b2stbr_redep_alpha', s%b2stbr_redep_alpha)
      call ipgetr ('b2stbr_recycled_neutrals_contr',
     & s%recycled_neutrals_contr)
      call ipgetr ('b2stbr_sput_phys_alpha', s%sput_phys_alpha)
      call ipgetr ('b2stbr_first_flight_dl', s%b2stbr_in_dl)

c b2sihs_
      call ipgetr ('b2sihs_rf0', s%b2sihs__rf0)
      call ipgetr ('b2sihs_rf1', s%b2sihs__rf1)
      call ipgetr ('b2sihs_rf2', s%b2sihs__rf2)
      call ipgetr ('b2sihs_rf3', s%b2sihs__rf3)
      call ipgetr ('b2sihs_rf4', s%b2sihs__rf4)

c b2news_
      call ipgeti ('b2news_ncallout', s%b2news_ncallout)

c b2srst
      call ipgetr ('b2srst_rf0', s%b2srst_rf0)
      call ipgetr ('b2srst_rf1', s%b2srst_rf1)
      call ipgetr ('b2srst_rf2', s%b2srst_rf2)
      call ipgetr ('b2srst_rf3', s%b2srst_rf3)

c b2tfnb
      call ipgetr ('b2tfnb_flux_limit_min_ti',
     &  s%b2tfnb_flux_limit_min_ti)
      s%b2tfnb_flux_limit_min_ti=s%b2tfnb_flux_limit_min_ti*ev
      call ipgeti ('b2tfnb_no_hybr', s%b2tfnb_no_hybr)

c b2srsm
      call ipgeti ('b2srsm_enable', s%b2srsm_enable)

c b2mndt
      call ipgeti ('b2mndt_density_control',s%density_control)
      call ipgeti ('b2mndt_nstg0', s%nstg(0))
      call ipgeti ('b2mndt_nstg1', s%nstg(1))
      call ipgeti ('b2mndt_nstg2', s%nstg(2))
      call ipgeti ('b2mndt_moqtlv', s%b2mndt_moqtlv)
      call ipgeti ('b2mndt_moitlv', s%b2mndt_moitlv)
      call ipgetr ('b2mndt_nstg_areshe', s%nstg_areshe)
      call ipgetr ('b2mndt_nstg_areshi', s%nstg_areshi)
      call ipgetr ('b2mndt_nstg_aresco', s%nstg_aresco)
      call ipgetr ('b2news_BoRiS', s%BoRiS)
      call ipgetr ('b2mndt_rxf', s%b2mndt_rxf)
      call ipgetr ('b2news_facExB_start',   s%facExB_start)
      call ipgetr ('b2news_facdrift_start', s%facdrift_start)
      call ipgetr ('b2news_facvis_start',   s%facvis_start)

c b2mndr
      call ipgetr ('b2mndr_hz',   s%b2mndr_hz)                       !srv 09.01.01 09.06.08
      call ipgetr ('b2mndr_stim', s%b2mndr_stim)
      call ipgetr ('b2mndr_rescale_neutrals', s%neutral_rescale)


c b2ag
      call ipgetr ('b2agmt_1d_width', s%b2agmt_width)
      call ipgetr ('b2agfs_geom_match_dist', s%geom_match_dist)

c b2siav
      call ipgetr ('b2siav_addvis', s%addvis)                        !srv 23.07.99
      call ipgetr ('b2siav_addvis1', s%addvis1)                      !srv 04.10.99
      call ipgetr ('b2siav_cqip1', s%cqip1)
      call ipgeti ('b2siav_iout', s%b2siav_iout)
      call ipgeti ('b2siav_style_qip', s%b2siav_qstyle)                          !srv 28.03.06

c b2trcl
      call ipgetr ('b2trcl_cvsa_mltpl', s%cvsa_mltpl)                            !iyv 06.12.13

c b2stbc
        call ipgeti ('b2stbc_bcene_16_style', s%bcene_16_style)
        call ipgeti ('b2stbc_bceni_16_style', s%bceni_16_style)

        call ipgetr ('b2stbc_cbc', s%b2stbc_cbc)                              !srv 14.10.99
        call ipgetr ('b2stbc_bc_ref_ti', s%b2stbc_bc_ref_ti)                  !srv 11.04.13
        call ipgetr ('b2stbc_bc_ref', s%b2stbc_bc_ref)
        call ipgetr ('b2stbc_bc_ref_te', s%b2stbc_bc_ref_te)                  !srv 11.04.13
        call ipgetr ('b2stbc_corr_flux_in_core', s%b2stbc_corr_flux)          !srv 25.07.05
        call ipgetr ('b2stbc_type13_norm', s%bc_type13_norm)
        call ipgetr ('b2stbc_type13_fac', s%bc_type13_fac)
        call ipgetr ('b2stbc_nesepm', s%nesepm)
        call ipgetr ('b2stbc_ndes_sol', s%ndes_sol)
        call ipgetr ('b2stbc_nesepm_overshoot', s%nesepm_overshoot)

c b2trno
        call ipgetr ('b2trno_con_e_stochastic', s%b2trno_alpha_stoch)         !lk 31.01.08

c Feedback
      call ipgeti ('b2stbc_feedback', s%b2stbc_feedback)

c Med
      call ipgeti ('b2mndr_med_style', s%med_style)

      return
      end subroutine read_switches

***********************************************************************

      subroutine check_switches (s,ns)

      implicit none
      integer, intent(in) :: ns
      type(switches), intent(in out) :: s

c check the bounds of the switches themselves
      call check_values_switches(s,ns)
c check the inter-consistency of different switches (e.g. certain combinations of swithces may not (yet) be allowed or developed)
      call check_consistency_switches(s)

      return
      end subroutine check_switches

***********************************************************************

      subroutine check_values_switches (s,ns)

c check the bounds of the switches themselves

      implicit none
      integer, intent(in) :: ns
      type(switches), intent(in out) :: s

c tn_style
c tn_style = 0: default combined ion-neutral energy equation
c tn_style = 1: pure-ion energy equation;
c               after solving ti-equation, tn is set equal to ti
c tn_style = 2: separate neutral energy equation for hydrogenic atoms (H/D/T) solved for tn
c               all ions and non-hydrogenic atoms combined into ti-equation
      if (s%tn_style.lt.0.or.s%tn_style.gt.2) then
        call xerrab('tn_tyle should be 0, 1 or 2')
      endif

      if (s%tn_style.ne.0.and.s%BoRiS.ne.0.0_R8) then
        call xerrab('tn_tyle.ne.0 not compatible with BoRiS.ne.0.0_R8')
      endif

c get_residuals
      if (.not.(s%get_residuals.eq.0.or.s%get_residuals.eq.1)) then
        call xerrab('get_residuals is expected to be 0 or 1')
      endif

c mdf_fhe
      if (.not.(s%mdf_fhe.eq.0.or.s%mdf_fhe.eq.1)) then
        call xerrab('mdf_fhe is expected to be 0 or 1')
      endif

c mdf_fhi
      if (.not.(s%mdf_fhi.eq.0.or.s%mdf_fhi.eq.1)) then
        call xerrab('mdf_fhi is expected to be 0 or 1')
      endif

c mdf_fnb
      if (.not.(s%mdf_fnb.eq.0.or.s%mdf_fnb.eq.1)) then
        call xerrab('mdf_fnb is expected to be 0 or 1')
      endif

c pot_eq
c pot_eq = 0: don't solve pot. eq.;
c pot_eq = 1: solve (default)
c pot_eq = 2: pot. = 3.1*Te/qe
      if (s%pot_eq.lt.0.or.s%pot_eq.gt.2) then  !srv 16.07.12 27.02.13
        call xerrab('faulty argument pot_eq')
      endif

c b2tfhe_no_current
      if (.not.(s%no_current.eq.0.or.s%no_current.eq.1)) then
        call xerrab('b2tfhe_no_current is expected to be 0 or 1')
      endif

c fch_anomalous
      call xertst(s%fch_anomalous.ge.0.0_R8 .and.
     &  s%fch_anomalous.le.1.0_R8,' Invalid fch_anomalous')

c use_eirene
      if (.not.(s%use_eirene.eq.0.or.s%use_eirene.eq.1)) then
        call xerrab('use_eirene is expected to be 0 or 1')
      endif

c b2sicf_style
       call xertst (0.le.s%b2sicf_style.and.s%b2sicf_style.lt.2,
     &   'faulty internal parameter b2sicf_style')

c b2sigp_style
       call xertst (0.le.s%b2sigp_style.and.s%b2sigp_style.lt.3,
     &   'b2sigp_style should be 0, 1 or 2')

c phm factors in b2sihs                       !srv 02.06.99
       call xertst (0.0_R8.le.min(s%b2sihs_phm0,s%b2sihs_phm1,
     &  s%b2sihs_phm2,s%b2sihs_phm3,s%b2sihs_phm4,
     &  s%b2sihs_phm5,s%b2sihs_phm6),
     &   'faulty parameter b2sihs phm0/1/2/3/4/5/6')

c b2sqcx
       call xertst (0.le.s%b2sqcx_styl0,
     &          'faulty parameter b2sqcx_styl0')
       call xertst (0.0_R8.le.s%b2sqcx_phm0,
     &           'faulty parameter b2sqcx_phm0')

c b2sqel
       call xertst (0.0_R8.le.min(s%b2sqel_phm0,
     &               s%b2sqel_phm1,s%b2sqel_phm2),
     &  'faulty parameter b2sqel phm0, phm1, phm2')

c b2stel
       call xertst (0.le.s%b2stel_styl0.and.s%b2stel_styl0.lt.4,
     &   'faulty parameter b2stel styl0')

       call xertst (0.0_R8.le.s%b2stel_rxm0.and.
     &   0.0_R8.le.s%b2stel_rxm1.and.0.0_R8.le.s%b2stel_rxm2.and.
     &   0.0_R8.le.s%b2stel_rxm3.and.0.0_R8.le.s%b2stel_rxm4,
     &   'faulty parameter b2stel rxm0/1/2/3/4')

       call xertst (0.0_R8.le.s%b2stel_phm0,
     &                          'faulty parameter b2stel phm0')

       call xertst (0.0_R8.le.s%b2stel_rg0.and.0.0_R8.le.s%b2stel_rg1,
     &   'faulty parameter b2stel rg0/1')

c b2srdt
        call xertst (0.0_R8.le.s%b2srdt_phm0.and.
     &        0.0_R8.le.s%b2srdt_phm1.and.0.0_R8.le.s%b2srdt_phm3.and.
     &        0.0_R8.le.s%b2srdt_phm4.and.0.0_R8.le.s%b2srdt_phm5,    !xpb
     &    'faulty internal parameter phm0, phm1, phm3, phm4, phm5')   !xpb

c b2stcx
       call xertst (0.0_R8.le.s%b2stcx_rg0,
     &            'faulty parameter b2stcx rg0')

       call xertst (0.le.s%b2stcx_styl0.and.s%b2stcx_styl0.lt.4,
     &   'faulty parameter b2stcx styl0')

c b2stbc
       call xertst (0.0_R8.le.s%sna0ep,
     .  'faulty parameter b2stbc: sna0ep')
       call xertst (0.0_R8.le.s%she0ep,
     .  'faulty parameter b2stbc: she0ep')
       call xertst (0.0_R8.le.s%shi0ep,
     .  'faulty parameter b2stbc: shi0ep')
       call xertst (0.0_R8.le.s%b2stbc_neoclassical .and.                  !xpb
     &              s%b2stbc_neoclassical.le.1.0_R8,                       !xpb
     1              'faulty iparameter b2stbc_neoclassical')               !xpb

       call xertst (0.0_R8.le.s%b2stbc_fchy_dia .and.                      !xpb
     &              s%b2stbc_fchy_dia.le.1.0_R8,                           !xpb
     1              'faulty internal parameter b2stbc_fchy_dia')           !xpb

c b2sihs
       call xertst (0.le.s%b2sihs_style,
     &                        'faulty parameter b2sihs style')

c b2sral
       call xertst (1.ne.s%b2sral_style,
     .  'b2sral_style = 1 not yet developed for WG')

       call xertst (0.le.s%b2sral_style.and.s%b2sral_style.le.2,           !srv 03.11.17
     .  'b2sral_style should be 0, 1 or 2')


c b2npco
       call xertst (0.lt.s%b2npco_rxg, 'faulty parameter b2npco rxg')
       call xertst (0.le.s%b2npco_pcm0, 'faulty parameter b2npco pcm0')
       call xertst (0.le.s%b2npco_pcm1, 'faulty parameter b2npco pcm1')

c b2npht
       call xertst (0.le.s%b2npht_pcm0, 'faulty parameter b2npht pcm0')
       call xertst (0.le.s%b2npht_pcm1, 'faulty parameter b2npht pcm1')
       call xertst (0.lt.s%b2npht_rxg, 'faulty parameter b2npht rxg')

c b2nxfc
       call xertst (0.le.s%b2nxfc_style.and.s%b2nxfc_style.lt.2,
     &   'b2nxfc style should be 0 or 1')

c b2npmo
        call xertst(0.0_R8.lt.s%b2npmo_rxg,
     &                                'faulty parameter b2npmo rxg')

c nsmin and nxmax
       if (s%nsmin.gt.ns) then
         s%nsmin=ns
       elseif (s%nsmin.lt.0) then
         s%nsmin=0
       endif

       if (s%nsmax.gt.ns) then
         s%nsmax=ns
       elseif( s%nsmax.lt.0) then
         s%nsmax=0
       endif

c b2news
        call xertst (0.eq.s%icase_ii.or.1.eq.s%icase_ii,
     &                      'b2tlnl_ii should be 0 or 1')

        call xertst (0.eq.s%icase_ee .or. 2.eq.s%icase_ee,
     &               'faulty input b2tlnl_ee')

        call xertst (0.eq.s%icase_ei .or. 3.eq.s%icase_ei,
     &              'faulty input b2tlnl_ei')

       call xertst (s%b2mndr_na_min.gt.0.0_R8,
     &  'faulty parameter b2mndr na_min')
       call xertst (s%b2mndr_na_new.gt.0.0_R8,
     &  'faulty parameter b2mndr na_new')
       call xertst (s%b2mndr_na_new.gt.s%b2mndr_na_min,
     &  'na_new must be larger than na_min!')

       call xertst (0.le.s%stylec, 'faulty parameter stylec')

       call xertst (0.0_R8.le.s%xfm0.and.s%xfm0.le.1.0_R8,
     &   'xfm0 should be between 0 and 1')
       call xertst (0.0_R8.le.s%xfm1.and.s%xfm1.le.1.0_R8,
     &   'xfm1 should be between 0 and 1')
       call xertst (0.0_R8.le.s%xfm2.and.s%xfm2.le.1.0_R8,
     &   'xfm2 should be between 0 and 1')
       call xertst (0.0_R8.le.s%xfm3.and.s%xfm3.le.1.0_R8,
     &   'xfm3 should be between 0 and 1')

       call xertst (-15.le.s%no_solve.and.s%no_solve.le.1,
     &   'faulty input no_solve')

c b2usmo
       call xertst (0.0_R8.le.s%cfc0, 'faulty parameter cfc0')

c b2tfcc
       call xertst (0.0_R8.le.s%xfac, 'faulty parameter xfac')

c b2tqca
       call xertst (0.0_R8.le.s%b2tqca_phm0,'faulty input b2tqca phm0')

c b2treq
       call xertst (0.0_R8.le.s%b2treq_phm0,'faulty input b2treq phm0')

c b2trcl
        call xertst (0.0_R8.le.s%cthiv, 'faulty parameter cthiv')
        call xertst (0.0_R8.le.s%cthev, 'faulty parameter cthev')

c b2lnl
        call xertst(0.0_R8.ne.s%b2trcl_lambda,
     &                           'faulty internal parameter lambda') !srv 02.07.13

c neutral flux limiters
        call xertst (0.0_R8.le.s%b2tlv0_alpha,
     &                             'invalid parameter b2tlv0 alpha')
        call xertst (0.0_R8.le.s%b2tlv0_gamma,
     &                             'invalid parameter b2tlv0 gamma')
        call xertst (0.0_R8.le.s%b2tlh0_alpha,
     &                             'invalid parameter b2tlh0 alpha')
        call xertst (0.0_R8.le.s%b2tlh0_gamma,
     &                             'invalid parameter b2tlh0 gamma')
        call xertst (0.0_R8.le.s%b2tlc0_alpha,
     &                             'invalid parameter b2tlc0 alpha')
        call xertst (0.0_R8.le.s%b2tlc0_gamma,
     &                             'invalid parameter b2tlc0 gamma')

c b2mwqt
       call xertst (0.le.s%b2mwqt_style,
     &                              'faulty parameter b2mwqt style')

c k-enstrophy model
      call xertst (0.le.s%solve_keps .and. s%solve_keps.le.2,
     &  'b2mndr_solve_keps value not coded')
      call xertst (0.le.s%keps_local .and. s%keps_local.le.1,
     &  'b2tqna_keps_local value not coded')
      call xertst (0.le.s%transport_keps .and. s%transport_keps.le.1,
     &  'b2tqna_transport_keps value not coded')
      call xertst (0.0_R8.lt.s%keps_dkt,
     &  'keps_dkt must be .gt. 0.0')
      call xertst (0.0_R8.lt.s%keps_dzt,
     &  'keps_dzt must be .gt. 0.0')
      call xertst (0.0_R8.lt.s%keps_cd,
     &  'keps_cd must be .gt. 0.0')
      call xertst (0.0_R8.lt.s%keps_heat,
     &  'keps_heat must be .gt. 0.0')
      call xertst (0.0_R8.lt.s%keps_heat_i,
     &  'keps_heat_i must be .gt. 0.0')
      call xertst (0.0_R8.le.s%keps_sig,
     &  'keps_sig must be .ge. 0.0')
      call xertst (0.0_R8.le.s%keps_alf,
     &  'keps_alf must be .ge. 0.0')
      call xertst (0.0_R8.lt.s%dna_min,
     &  'dna_min must be .gt. 0.0')
      call xertst (0.0_R8.lt.s%vsa_min,
     &  'vsa_min must be .gt. 0.0')
      call xertst (0.0_R8.lt.s%hce_min,
     &  'hce_min must be .gt. 0.0')
      call xertst (0.0_R8.lt.s%hci_min,
     &  'hci_min must be .gt. 0.0')
      call xertst (0.0_R8.lt.s%keps_visc,
     &  'keps_visc must be .gt. 0.0')
      call xertst (0.0_R8.le.s%keps_shear,
     &  'keps_shear must be .ge. 0.0')
      call xertst (0.0_R8.le.s%b2sikt_fac_aniso,
     &  'b2sikt_fac_aniso must be .ge. 0.0')
      call xertst (0.0_R8.le.s%b2sikt_fac_sheath,
     &  'b2sikt_fac_sheath must be .ge. 0.0')
      call xertst (0.0_R8.le.s%b2sikt_fac_sheath_core,
     &  'b2sikt_fac_sheath_core must be .ge. 0.0')
      call xertst (0.0_R8.le.s%b2sikt_fac_diss,
     &  'b2sikt_fac_diss must be .ge. 0.0')
      call xertst (0.0_R8.le.s%b2sikt_fac_diss_core,
     &  'b2sikt_fac_diss_core must be .ge. 0.0')
      call xertst (0.0_R8.le.s%b2tfhi_fflokt,
     &  'b2tfhi_fflokt must be .ge. 0.0')
      call xertst (0.0_R8.le.s%b2tfhi_fconkt,
     &  'b2tfhi_fconkt must be .ge. 0.0')
      call xertst (0.0_R8.le.s%b2tfhi_fflozt,
     &  'b2tfhi_fflozt must be .ge. 0.0')
      call xertst (0.0_R8.le.s%b2tfhi_fconzt,
     &  'b2tfhi_fconzt must be .ge. 0.0')
      call xertst (0.0_R8.le.s%b2tfhi_fsigkt,
     &  'b2tfhi_fsigkt must be .ge. 0.0')
      call xertst (0.0_R8.le.s%b2tfhi_fkt_hie,
     &  'b2tfhi_fkt_hie must be .ge. 0.0')
      call xertst (0.0_R8.le.s%keps_inc,
     &  'keps_inc must be .ge. 0.0')
      call xertst (0.0_R8.le.s%keps_fac,
     &  'keps_fac must be .ge. 0.0')
      call xertst (0.le.s%b2sikt_local .and. s%b2sikt_local.le.1,
     &  'b2sikt_local must be either 0 or 1')
      call xertst (0.le.s%b2sikt_sheath_local .and.
     &  s%b2sikt_sheath_local.le.1,
     &  'b2sikt_sheath_local must be either 0 or 1')
      call xertst (0.le.s%b2sikt_min_source .and.
     &  s%b2sikt_min_source.le.2,
     &  'b2sikt_min_source must be either 0, 1 or 2')

c eirene_mc
      call xertst (0.0_R8.le.s%neutral_sources_rescale,
     1   'faulty input neutral_sources_rescale')

c b2tqna
      call xertst (0.0_R8.le.s%b2tqna_min_df0.and.s%b2tqna_min_df0.le.
     &   s%b2tqna_max_df0.and.s%b2tqna_max_df0.gt.0.0_R8,
     &   'Invalid range for b2tqna_min_df0 and/or b2tqna_max_df0')
      call xertst (s%b2tqna_pfr_rescale.gt.0.0_R8,
     &  'faulty input pfr_rescale')
      call xertst (s%b2tqna_divsol_rescale.gt.0.0_R8,
     &   'faulty input divsol_rescale')
      if (s%user_transport.eq.10) then
        write(*,*) 'user_transport = 10 is now obsolete, because it',
     &      'did not allow to combine advanced neutral transport',
     &      'with non-default plasma transport.'
        call xerrab('use "b2tqna_transport_afn = 1" instead')
      endif

c b2stbm
      call xertst (s%b2stbm_linearisation.le.1.0_R8 .and.
     &  s%b2stbm_linearisation.ge.0.0_R8,
     &  'faulty input parameter b2stbm_linearisation')
      call xertst (s%b2stbm_internal_energy_sources.eq.0,
     &  'b2stbm internal_energy_sources =/ 0 not coded yet')

c b2stbr
      call xertst (0.0_R8.le.s%b2stbr_core_sources_rescale,
     &  'faulty internal parameter core_sources_rescale')
      call xertst (0.0_R8.le.s%b2stbr_alpha,'faulty b2stbr_alpha')
      call xertst (0.0_R8.le.s%b2stbr_sput_frc,'faulty b2stbr_sput_frc')
      call xertst (s%b2stbr_sput_mixed_alpha.ge.0.0_R8,
     &  'faulty parameter sput_mixed_alpha')
      call xertst (s%b2stbr_sput_mixed_beta.gt.0.0_R8,
     &  'faulty parameter sput_mixed_beta')
      call xertst (0.0_R8.le.s%b2stbr_sput_chem_alpha,
     &  'faulty b2stbr_sput_chem_alpha')
      call xertst (0.0_R8.le.s%b2stbr_sput_chem_cutoff_alpha,
     &  'faulty b2stbr_sput_chem_cutoff_alpha')
      call xertst (0.0_R8.le.s%b2stbr_sput_chem_cutoff_beta,
     &  'faulty b2stbr_sput_chem_cutoff_beta')
      call xertst (0.0_R8.le.s%b2stbr_therm_evap,
     &  'faulty b2stbr_therm_evap')
      call xertst (0.0_R8.le.s%b2stbr_sput_res,
     &  'faulty b2stbr_sput_res')
      call xertst (0.0_R8.le.s%b2stbr_sput_phys,
     &  'faulty input b2stbr_sput_phys')
      call xertst (0.0_R8.le.s%b2stbr_redep_alpha,
     &  'faulty internal parameter redep_alpha')
      call xertst (0.0_R8.le.s%recycled_neutrals_contr,
     &  'faulty internal parameter recycled_neutrals_contr')
      call xertst (0.0_R8.le.s%sput_phys_alpha,
     &  'faulty b2stbr_sput_phys_alpha')

c b2sihs_
      call xertst (0.0_R8.le.min(s%b2sihs__rf0,s%b2sihs__rf1,
     &  s%b2sihs__rf2,s%b2sihs__rf3,s%b2sihs__rf4),
     &  'faulty internal parameter rf0/1/2/3/4')

c b2srst
      call xertst (0.0_R8.le.s%b2srst_rf0.and.0.0_R8.le.s%b2srst_rf1
     & .and.0.0_R8.le.s%b2srst_rf2.and.0.0_R8.le.s%b2srst_rf3,
     & 'faulty internal parameter rf0/1/2/3')

c b2tfnb
      call xertst (0.0_R8.le.s%b2tfnb_flux_limit_min_ti,
     &  'faulty internal parameter flux_limit_min_ti')

c b2mndt
      call xertst (1.le.s%nstg(0).and.1.le.s%nstg(1).and.1.le.s%nstg(2),
     &  'faulty internal parameter nstg(0:2)')
      call xertst (-1.le.s%b2mndt_moqtlv.and.s%b2mndt_moqtlv.le.3,
     &  'faulty internal parameter moqtlv')
      call xertst (-1.le.s%b2mndt_moitlv.and.s%b2mndt_moitlv.le.3,
     &   'faulty internal parameter moitlv')
      call xertst (s%nstg_areshe.ge.0.0_R8,
     &  'faulty internal parameter nstg_areshe')
      call xertst (s%nstg_areshi.ge.0.0_R8,
     &  'faulty internal parameter nstg_areshi')
      call xertst (s%nstg_aresco.ge.0.0_R8,
     &  'faulty internal parameter nstg_aresco')
      call xertst (0.0_R8.le.s%b2mndt_rxf.and.s%b2mndt_rxf.le.1.0_R8,
     &  'faulty internal parameter rxf')

c b2mndr
       call xertst (0.0.le.s%b2mndr_hz .and. 1.0.ge.s%b2mndr_hz,
     &  'faulty parameter b2mndr_hz')                          !srv 09.01.01 09.06.08 }

c b2ag
      call xertst (s%geom_match_dist.gt.0.0_R8,
     &  'faulty argument geom_match_dist')

c med
      call xertst (s%med_style.ge.0 .and. s%med_style.le.1,
     &  'faulty parameter b2mndr_med_style')

      return
      end subroutine check_values_switches

***********************************************************************

      subroutine check_consistency_switches (s)
c      use b2mod_boundary_namelist
c     & , only : nbc
c check the inter-consistency of different switches (e.g. certain combinations of swithces may not (yet) be allowed or developed)

      implicit none
      type(switches), intent(in out) :: s

      if (s%pot_eq.ne.1.and.s%no_current.ne.1) then
        call xerrab ('If b2news_poteq.ne.1,'//
     &                  'then b2tfhe_no_current must be 1!')
      endif

      if (s%tn_style.ne.0.and.s%stylec.ne.0) then
        call xerrab('b2upht stylec =/= 0 not yet available for'//
     &                  'non-default energy equation')
      endif

      if(s%fch_pTe .eq. 1.0_R8) s%alfTeEh = 0.0_R8                          !srv 16.04.10
      call xertst (.not.(0.0_R8.ne.s%alfTeEh.and.0.0_R8.ne.s%fch_pTe),
     &        ' Do not set 0.0_R8.ne.alfTeEh.and.0.0_R8.ne.fch_pTe')

      call xertst (.not.(0.ne.s%b2tfnb_alpha.and.0.ne.s%b2tlc0_alpha),      !srv 14.10.09
     &   'b2tfnb_alpha and b2tlc0_alpha cannot be both different '//
     &   'from zero ')

      if (s%pot_eq .eq. 1) call xertst(s%prl_cur .gt. 0.0e0_R8,
     &  'No parallel current!')

      if (s%pot_eq .eq. 1) then
       call xertst(s%fch_ion_neutral .gt. 0.0e0_R8 .or.
     *             s%fch_anomalous .gt. 0.0e0_R8 .or.
     *             s%fhe_vis_q .gt. 0.0_R8 .or.  !fnb/fhe?
     *             s%fch_inert .gt. 0.0_R8 .or.
     *             s%fhe_vis_per .gt. 0.0_R8 .or. !fnb/fhe?
     *             s%fch_stochastic .gt. 0.0_R8 .or.
     *             s%dia_cur .gt. 0.0e0_R8 .or.
     *             s%fhe_vis_par .gt. 0.0_R8,
     *             'No radial current terms!')
      endif

c AFN neutral model
      if(s%AFN_out.ne.0.and.s%recycle_afn.eq.0) then
        call xerrab('b2stbr_recycle_afn should be 1 when writing '//
     &   'AFN output')
      endif

c spatially hybrid neutrals
      if(s%spatial_hybrid.ne.0.and.s%use_eirene.eq.0) then
        call xerrab('spatial hybrid neutrals must be with eirene')
      endif

c      if(s%spatial_hybrid.ne.0.and.s%pot_eq.eq.1) then
c        call xerrab('spatial hybrid neutrals not yet developed'//
c     &               'with drifts and currents')
c      endif

      if(s%use_auto_spatial_hyb.ne.0.and.s%spatial_hybrid.eq.0) then
        call xerrab('automatic spatial hybrid needs spatial hybrid')
      endif

      if(s%use_auto_spatial_hyb.ne.0.and.s%recycle_afn.eq.0) then
        call xerrab('automatic spatial hybrid needs KUL neutral BCs')
      endif

c k-enstrophy model
      if (s%solve_keps.eq.0 .and. s%transport_keps.gt.0) then
        write(*,*) 'WARNING! b2mndr_solve_keps = ',
     &    s%solve_keps,' and b2tqna_transport_keps =',
     &    s%transport_keps,' should go together'
      endif

c b2tfnb b2tfrn, this needs to be here otherwise clashes with Tapenade adjoint
      if (s%b2tfnb_no_hybr .eq. 1) s%b2tfnb_discr_meth = 0

      return
      end subroutine check_consistency_switches

***********************************************************************



      end module b2mod_switches
