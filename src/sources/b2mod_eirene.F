      module b2mod_eirene
      use b2mod_types
      use b2mod_plasma
      use b2mod_indirect
      use b2mod_geo
      use b2mod_user_namelist
      use b2mod_neutrals_namelist
      use b2mod_numerics_namelist
      use b2mod_boundary_namelist
      use b2mod_constants
      use b2mod_eirene_globals
      use b2mod_tallies
      use b2mod_indirect
      use b2mod_neutr_src_scaling
      use b2mod_feedback
      use b2mod_diag
      use b2mod_time
      use b2mod_rates
      use b2mod_trace
      use b2mod_b2cmpa
      use b2mod_b2cmpb
      use b2mod_balance !djm Jan2017
     & ,only : eirene_mc_papl_sna0to1, eirene_mc_pmpl_sna0to1,
     &         eirene_mc_pael_sne0to1, eirene_mc_pmel_sne0to1,
     &         eirene_mc_pipl_sna0to1, eirene_mc_pppl_sna0to1,
     &         eirene_mc_core_sna0to1, eirene_mc_mapl_smo0to3,
     &         eirene_mc_mmpl_smo0to3, eirene_mc_mipl_smo0to3,
     &         eirene_mc_cppv_smo0to3, eirene_mc_eael_she0to3,
     &         eirene_mc_emel_she0to3, eirene_mc_eiel_she0to3,
     &         eirene_mc_epel_she0to3, eirene_mc_eapl_shi0to3,
     &         eirene_mc_empl_shi0to3, eirene_mc_eipl_shi0to3,
     &         eirene_mc_eppl_shi0to3, eirene_mc_paat_sna_bal,
     &         eirene_mc_pmat_sna_bal, eirene_mc_piat_sna_bal,
     &         eirene_mc_paml_sna_bal, eirene_mc_pmml_sna_bal,
     &         eirene_mc_piml_sna_bal, eirene_mc_paio_sna_bal,
     &         eirene_mc_pmio_sna_bal, eirene_mc_piio_sna_bal,
     &         alloc_b2mod_eirene_sources,
     &         balance_netcdf
#ifdef B25_EIRENE
      use b2mod_work
      use eirmod_braeir
      use eirmod_eirbra
      use eirmod_ctrcei
      use eirmod_comusr
     , , only: texts, nsph, nspa, nspam
      use eirmod_cinit
     , , only: fort_lc
      use eirmod_wneutrals
      use eirmod_clgin
     ,   , only : iliin
#else
      use b2mod_braeir
      use b2mod_eirbra
#endif
      implicit none

      integer, save :: ank_mods = 0, use_recyceir = 1
      integer, save :: nhist = 1, nrelax = 0
      real (kind=R8), save, allocatable ::
     .  sni_av(:,:,:,:), smo_av(:,:,:,:), sei_av(:,:,:), see_av(:,:,:)

      real (kind=R8), save, allocatable ::                        !djm Jan2017
     .  sne_av_pael(:,:,:), sne_av_pmel(:,:,:),
     .  sna_av_paat(:,:,:,:), sna_av_pmat(:,:,:,:),
     .  sna_av_piat(:,:,:,:), sna_av_paml(:,:,:,:),
     .  sna_av_pmml(:,:,:,:), sna_av_piml(:,:,:,:),
     .  sna_av_paio(:,:,:,:), sna_av_pmio(:,:,:,:),
     .  sna_av_piio(:,:,:,:),
     .  sni_av_papl(:,:,:,:), sni_av_pmpl(:,:,:,:),
     .  sni_av_pipl(:,:,:,:), sni_av_pppl(:,:,:,:),
     .  smo_av_mapl(:,:,:,:), smo_av_mmpl(:,:,:,:),
     .  smo_av_mipl(:,:,:,:), smo_av_cppv(:,:,:,:),
     .  see_av_eael(:,:,:), see_av_emel(:,:,:),
     .  see_av_eiel(:,:,:), see_av_epel(:,:,:),
     .  sei_av_eapl(:,:,:), sei_av_empl(:,:,:),
     .  sei_av_eipl(:,:,:), sei_av_eppl(:,:,:)
      real (kind=R8), allocatable, save ::                        !IYS 28.03.12 {
     &  sni_h(:,:,:,:,:),
     &  smo_h(:,:,:,:,:),
     &  see_h(:,:,:,:), sei_h(:,:,:,:),
     &  pfluxa_h(:,:,:,:,:), rfluxa_h(:,:,:,:,:),
     &  pfluxm_h(:,:,:,:,:), rfluxm_h(:,:,:,:,:),
! _h is an abbreviation from 'history'
     &  sni_s(:,:,:,:),
     &  smo_s(:,:,:,:),
     &  see_s(:,:,:), sei_s(:,:,:),
     &  pfluxa_s(:,:,:,:), rfluxa_s(:,:,:,:),
     &  pfluxm_s(:,:,:,:), rfluxm_s(:,:,:,:),
! _s is an abbreviation from 'sum'
     &  dab2_h(:,:,:,:,:), dab2_s(:,:,:,:)
      real (kind=R8), save, allocatable ::                        !djm Jan2017
     &  sna_h_paat(:,:,:,:,:), sna_h_pmat(:,:,:,:,:),
     &  sna_h_piat(:,:,:,:,:), sna_h_paml(:,:,:,:,:),
     &  sna_h_pmml(:,:,:,:,:), sna_h_piml(:,:,:,:,:),
     &  sna_h_paio(:,:,:,:,:), sna_h_pmio(:,:,:,:,:),
     &  sna_h_piio(:,:,:,:,:),
     &  sni_h_papl(:,:,:,:,:), sni_h_pmpl(:,:,:,:,:),
     &  sni_h_pipl(:,:,:,:,:), sni_h_pppl(:,:,:,:,:),
     &  smo_h_mapl(:,:,:,:,:), smo_h_mmpl(:,:,:,:,:),
     &  smo_h_mipl(:,:,:,:,:), smo_h_cppv(:,:,:,:,:),
     &  sne_h_pael(:,:,:,:), sne_h_pmel(:,:,:,:),
     &  see_h_eael(:,:,:,:), see_h_emel(:,:,:,:),
     &  see_h_eiel(:,:,:,:), see_h_epel(:,:,:,:),
     &  sei_h_eapl(:,:,:,:), sei_h_empl(:,:,:,:),
     &  sei_h_eipl(:,:,:,:), sei_h_eppl(:,:,:,:),
     &  sna_s_paat(:,:,:,:), sna_s_pmat(:,:,:,:),
     &  sna_s_piat(:,:,:,:), sna_s_paml(:,:,:,:),
     &  sna_s_pmml(:,:,:,:), sna_s_piml(:,:,:,:),
     &  sna_s_paio(:,:,:,:), sna_s_pmio(:,:,:,:),
     &  sna_s_piio(:,:,:,:),
     &  sni_s_papl(:,:,:,:), sni_s_pmpl(:,:,:,:),
     &  sni_s_pipl(:,:,:,:), sni_s_pppl(:,:,:,:),
     &  smo_s_mapl(:,:,:,:), smo_s_mmpl(:,:,:,:),
     &  smo_s_mipl(:,:,:,:), smo_s_cppv(:,:,:,:),
     &  sne_s_pael(:,:,:), sne_s_pmel(:,:,:),
     &  see_s_eael(:,:,:), see_s_emel(:,:,:),
     &  see_s_eiel(:,:,:), see_s_epel(:,:,:),
     &  sei_s_eapl(:,:,:), sei_s_empl(:,:,:),
     &  sei_s_eipl(:,:,:), sei_s_eppl(:,:,:)

      real(kind=R8), allocatable, save :: hlp_stra(:)
      logical, allocatable, save :: hlp_pff(:,:)
      logical, save :: saved_volrec_used
      data saved_volrec_used/.false./

      contains

      subroutine eirene_mc(nx, ny, ns, dtim, BoRiS,
     1 sna0, smo0, she0, shi0, main_call)
      use b2mod_subsys
      implicit none
      integer nx,ny,ns
      real (kind=R8) ::
     * dtim, BoRiS,
     * sna0(-1:nx,-1:ny,0:1,0:ns-1), smo0(-1:nx,-1:ny,0:3,0:ns-1),
     & she0(-1:nx,-1:ny,0:3), shi0(-1:nx,-1:ny,0:3)
      logical main_call                                           !dpc
csw 19oct2012, added save attribute!!!
      logical, save :: my_ltime

      integer ix,iy,is,ib,ih,it,ifl,istrai,ncall,iatm,imol,iion
      integer pfrregno1, pfrregno2, solregno, isfeedback, ncallfeedback,
     1 nesepm_istra, b2stbc_feedback
      logical feedback_stratum
      real (kind=R8) ::
     * recyc_hlp,vthe,t0,t1,rf0,dua,t2,t3,t4,f,
     * neutral_sources_rescale, ndes, private_flux_puff, nesepm_sol,
     * recyctmp, nesepm_pfr, volrec_sol, ndes_sol, nepedm_sol
      real (kind=R8), save, allocatable ::
     1 SSNI(:,:), SSMO(:,:), SSEE(:), SSEI(:), SSEICORR(:),
     1 SSEISTAB(:),
     1 tfl(:,:), volsumnold(:)
      real (kind=R8), save, allocatable ::
     * SMO_TOT(:,:,:), SNI_TOT(:,:,:)
      real (kind=R8) ::                                           !djm Jan2017
     * t0_a, t0_m, t0_i, t0_p,
     * t1_apl, t1_mpl, t1_ipl, t1_ppl,
     * t2_apl, t2_mpl, t2_ipl, t2_ppl,
     * t3_apl, t3_mpl, t3_ipl, t3_ppl,
     * t4_apl, t4_mpl, t4_ipl, t4_ppl

cank-20020724 {
c*** The type of source scaling is determined by the lstrascl parameter
c*** for each stratum, each atom species:
c***  < 0   - prehistoric treatment (left for compatibility)
c***  = 0   - no source scaling at all
c***  =< natmi - full treatment, on each sub-iteration
c***  > natmi  - partial treatment, after Eirene call only
c***
c*** In the Eirene particle balance, we distinguish the following
c*** components for each stratum, each atomic species:
c***
c***  hlp_pi  corresponding plasma species onto the surfaces
c***  hlp_po  resulting neutrals from the surfaces (recycling)
c***  hlp_ni  neutrals onto the surfaces
c***  hlp_no  reflected neutrals
c***  hlp_g   neutrals from the gas puffing sources
c***  hlp_s   neutrals ionising in the plasma (<0 for recombination)
c***  hlp_t   neutrals absorbed in the "time layer"
c***  tfl     plasma onto the surfaces in B2
c***
c*** For a recycling stratum, first scale all the sources and fluxes by
c*** the ratio of tfl/hlp_pi to have exact match between B2 & Eirene.
c*** Then scale hlp_s, hlp_t, hlp_ni, hlp_no by x ~ 1:
c***
c***    (hlp_s + hlp_t + hlp_ni - hlp_no) * x = hlp_po
c***
c*** (this should only be done for the species scaled with themselves)
c***
c*** For a non-recycling stratum, scale hlp_s & hlp_t by x, and
c*** hlp_ni & hlp_no by 1/x:
c***
c***    (hlp_s + hlp_t) * x + (hlp_ni - hlp_no) / x = hlp_g
c***
c*** selecting the x value positive, closest to 1.
c***
c*** No scaling for the "time layer"
c***
c*** 10.02.2013: If l_neutrad<0, then old (incorrect) treatment of the neutral radiation.
c*** Otherwise: take the radiation explicitly calculated in EIRENE.

      real (kind=R8) :: s,x1,x2,t00,rec_tol
      real (kind=R8), save, allocatable ::
     1 hlp_pi(:),hlp_po(:),hlp_ni(:),
     2 hlp_no(:),hlp_g(:),hlp_s(:),hlp_t(:),
     1 hlp_scl(:),hlp_scp(:,:),
     2 core_flux(:)
      real (kind=R8), save :: core_surface
      logical hlp_bad_isp, hlp_log, hlp_pr, hlp_tt, hlp_cor
      integer i,j,k,l,js,ls,kks,n,icount
      integer, save :: l_ddt
      integer, save, allocatable :: ispc_ref(:), ishigh(:)
      character*13, save, allocatable :: ion_name(:)
      real (kind=R8), allocatable :: sum_wld(:,:)
cank }
      character*24 hlp_frm, hlp_frm1, hlp_frm2, hlp_frm3, hlp_frm4
      character*256 filename
c      integer lns_hlp
c      parameter (lns_hlp=4)
c      integer lis_hlp(lns_hlp)
c      data lis_hlp / 110, 111, 112, 49 /
      logical ex, last, found, na_feedback_on
      integer :: igpuff(NSTRAT)

      integer output_style, linearisation, save_f31,
     1 repeat_first_call, ionising_core, dpc_fix, eircount,
     2 no_sseicorr, save_f30, stab_corr
      save ncall, output_style,
     * pfrregno1, pfrregno2, solregno, isfeedback, b2stbc_feedback,
     * nesepm_sol, nesepm_pfr, ndes, private_flux_puff, nesepm_istra,
     * neutral_sources_rescale, ncallfeedback, save_f31,
     1 repeat_first_call, ionising_core, dpc_fix,
     2 volrec_sol, ndes_sol, no_sseicorr, nepedm_sol, save_f30,
     3 stab_corr
csw 26mar2013
      integer, save :: iscltst=0
csw 26feb2014
      integer, save :: ihist=0, neutr_avg=0
      integer, save :: ih_end
      logical, save :: full_history = .false.
      real (kind=R8) :: relax
      integer, save :: jxi, jxa, jsep
      integer, save :: nxtl, nxtr   ! IYS 29.03.2022

*   ..procedures
      real (kind=R8) :: rratio
      intrinsic min, mod, abs, sqrt, log
      external xertst, ipgeti, ipgetr, rratio
#ifdef B25_EIRENE
      logical streql
      external streql
#endif
*   ..initialisation
      data ncall/0/, output_style/1/,
     .  linearisation/1/, solregno/3/, b2stbc_feedback/0/,
     1  pfrregno1/0/, pfrregno2/2/, isfeedback/0/, ncallfeedback/0/,
     2  neutral_sources_rescale/1.0_R8/, rec_tol/1.e-4_R8/,
     3  nesepm_sol/0.0_R8/, nesepm_pfr/0.0_R8/, nesepm_istra/-1/,
     4  private_flux_puff/0.0_R8/, ionising_core/0/, dpc_fix/1/,
     5  ndes/0.0_R8/, save_f31/0/, repeat_first_call/1/,
     6  volrec_sol/0.0_R8/, ndes_sol/0.0_R8/,
     7  no_sseicorr/0/, nepedm_sol/0.0_R8/, save_f30/0/, stab_corr/1/

*-----------------------------------------------------------------------
*.computation

* ..preliminaries
*   ..subprogram start-up calls
      call subini('eirene_mc')
*   ..set internal parameters on first call
      if(ncall.eq.0) then !{
c not sure this is the right place
        allocate (hlp_pi(NATM))
        allocate (hlp_po(NATM))
        allocate (hlp_ni(NATM))
        allocate (hlp_no(NATM))
        allocate (hlp_g(NATM))
        allocate (hlp_s(NATM))
        allocate (hlp_t(NATM))
        allocate (hlp_scl(NATM))
        allocate (hlp_scp(NATM,NSTRAT))
        allocate (hlp_stra(NATM))
        allocate (core_flux(NATM))
        allocate (hlp_pff(NATM,NSTRAT))
        allocate (ion_name(NATM))
        allocate (tfl(NSTRAT,NATM))
        allocate (SSNI(NSTRAT,0:NS-1))
        allocate (SSMO(NSTRAT,0:NS-1))
        allocate (SSEE(NSTRAT))
        allocate (SSEI(NSTRAT))
        allocate (SSEICORR(NSTRAT))
        allocate (SSEISTAB(NSTRAT))
        allocate (SMO_TOT(0:NDXP,0:NDYP,NFL))
        allocate (SNI_TOT(0:NDXP,0:NDYP,NFL))
        allocate (volsumnold(NSTRAT))
        volsumnold=-1.0_R8
        allocate(ispc_ref(NATM))
        allocate(ishigh(NATM))
        allocate(sum_wld(4,2+NATM))
        allocate(feedback_strata(0:nspecies-1))
        if (balance_netcdf.ne.0) then !djm Jan2017
         call alloc_b2mod_eirene_sources(nx,ny,ns,natm,nmol,nion,nstrat)
        endif

        write(*,*) 'Calling test_couple'
        call test_couple(ndx,ndy,nstra,nfl,natm,nmol,nion,nlimps)
        write(*,*) 'Called test_couple'
        call ipgeti ('eirene_mc_output_style', output_style)
        call ipgeti ('eirene_repeat_first_call', repeat_first_call)
        call ipgeti ('eirene_mc_linearisation', linearisation)
        if (linearisation.eq.1)   ! for backward compatibility
     >   call ipgeti ('eirene_mc_linearization', linearisation)
        call xertst (linearisation.eq.1 .or. linearisation.eq.0,
     .   'faulty input parameter linearisation')
        call ipgeti ('b2stbc_feedback', b2stbc_feedback)
        lfeedback = b2stbc_feedback.ne.0
        call ipgetr ('b2stbc_ndes', ndes)
        call xertst (0.0_R8.le.ndes, 'faulty internal parameter ndes')
        call ipgetr ('b2stbc_nesepm_sol', nesepm_sol)
        call xertst (0.0_R8.le.nesepm_sol,
     1   'faulty internal parameter nesepm_sol')
        call ipgetr ('b2stbc_nepedm_sol', nepedm_sol)
        call xertst (0.0_R8.le.nepedm_sol,
     1   'faulty internal parameter nepedm_sol')
        call ipgetr ('b2stbc_nesepm_pfr', nesepm_pfr)
        call xertst (0.0_R8.le.nesepm_pfr,
     1   'faulty internal parameter nesepm_pfr')
        call ipgetr ('b2stbc_volrec_sol', volrec_sol)
        call xertst (0.0_R8.le.volrec_sol,
     1   'faulty internal parameter volrec_sol')
        call ipgetr ('b2stbc_ndes_sol', ndes_sol)
        call xertst (0.0_R8.le.ndes_sol,
     1   'faulty internal parameter ndes_sol')
        call ipgeti ('b2stbc_isfeedback', isfeedback)
        call xertst (0.le.isfeedback.and.isfeedback.lt.ns,
     1   'faulty internal parameter isfeedback')
        feedback_strata=0
        call ipgeti ('eirene_nesepm_istra', nesepm_istra)
        if(nesepm_istra.gt.0.and.lfeedback) then
          call xertst (nesepm_istra.le.nstrat,
     1     'faulty internal parameter nesepm_istra')
          call xertst (crcstra(nesepm_istra).eq.'C',
     1     'Stratum nesepm_istra is not declared as a puff stratum!')
        else
          do istrai = 1, nstrat
            if(nesepm_istra.le.0.and.CRCSTRA(istrai).EQ.'C') then
              nesepm_istra=istrai
            endif
          enddo
        endif
        if((ndes.gt.0.0_R8.or.nesepm_pfr.gt.0.0_R8.or.
     &    private_flux_puff.gt.0.0_R8.or.nesepm_sol.gt.0.0_R8.or.
     *    volrec_sol.gt.0.0_R8.or.ndes_sol.gt.0.0_R8.or.
     &    nepedm_sol.gt.0.0_R8).and.lfeedback) then
          call xertst (nesepm_istra.gt.0,
     &     'No "C" stratum for feedback scheme specified !')
          icount = 0
          do istrai = 1, nstrat
            if (crcstra(istrai).eq.'C' .and. lfeedback .and.
     .          species_start(istrai).le.isfeedback .and.
     .          species_end(istrai).ge.isfeedback) icount = icount+1
          end do
          call xertst (icount.eq.1 .or. nesepm_istra.gt.0,
     .     'Multiple definitions of puff strata, use nesepm_istra')
          feedback_strata(b2espcr(isfeedback)-1)=nesepm_istra
        endif
        do is = 0, nspecies-1
          if (na_feedback_actuator(is).eq.1) then
            do istrai = 1, nstrat
              if (crcstra(istrai).eq.'C' .and.
     .          species_start(istrai).le.eb2spcr(is+1) .and.
     .          species_end(istrai).ge.eb2spcr(is+1)) then
                  feedback_strata(is) = istrai
              end if
            end do
          end if
        end do
        call ipgetr ('b2stbc_private_flux_puff', private_flux_puff)
        call xertst (0.0_R8.le.private_flux_puff,
     1   'faulty internal parameter private_flux_puff')
        call ipgeti ('b2stbc_pfrregno1', pfrregno1)
        call ipgeti ('b2stbc_pfrregno2', pfrregno2)
        call ipgeti ('b2stbc_solregno', solregno)
        if (lfeedback) then
          if (ndes.gt.0.0_R8.or.private_flux_puff.gt.0.0_R8
     1       .or.nesepm_pfr.gt.0.0_R8) then !{
            call xertst (cbirso.le.pfrregno1.and.
     1       pfrregno1.le.cbirso+cbnrso-1,
     2       'faulty internal parameter pfrregno1')
            call xertst (cbirso.le.pfrregno2.and.
     1       pfrregno2.le.cbirso+cbnrso-1,
     2       'faulty internal parameter pfrregno2')
          endif !}
          if (nesepm_sol.gt.0.0_R8.or.volrec_sol.gt.0.0_R8.or.
     1        nepedm_sol.gt.0.0_R8.or.ndes_sol.gt.0.0_R8) then
            call xertst (cbirno.le.solregno.and.
     1       solregno.le.cbirno+cbnrno-1,
     2       'faulty internal parameter solregno')
          endif
        endif
        call ipgetr ('b2mndr_rescale_neutrals_sources',
     *       neutral_sources_rescale)
        call xertst (0.0_R8.le.neutral_sources_rescale,
     1   'faulty input neutral_sources_rescale')
        call ipgeti ('b2stbc_ncallfeedback', ncallfeedback)
        call ipgeti ('eirene_ionising_core', ionising_core)
        if (ionising_core.eq.0)     ! for backward compatibility
     >   call ipgeti ('eirene_ionizing_core', ionising_core)
        call xertst (ionising_core.ge.-nbc.and.ionising_core.le.1,
     1   'faulty input ionising_core, abs(ionising_core) > NBC')
        call ipgeti ('eirene_ank_mods', ank_mods)
        call ipgeti ('eirene_dpc_fix', dpc_fix)
        call ipgeti ('eirene_savef30', save_f30)
        call ipgeti ('eirene_savef31', save_f31)
        call ipgeti ('eirene_use_recyceir', use_recyceir)
csw 28jan2011, call to eirene_mc_init not used anymore, use new modules
c        call eirene_mc_init(nx,ny)
        call eirene_extrab25_alloc_mods
!pb 13dec2013 eirpbls removed - not needed
!pb        call eirene_extrab25_eirpbls_init(nmol,nion,npls)
csw 22jul2011
        call eirene_extrab25_iniusr_init(n_spcsrf,l_spcsrf,i_spcsrf,
     .           j_spcsrf,sps_sgrp,sps_absr,sps_trno,sps_trni,
     .           sps_mtri,sps_tmpr,sps_spph,sps_spch,
     .           sps_mtrl,sps_id)
#ifndef NO_B2_CHEM_SPUT
        call eirene_extrab25_braeir_init(chemical_sputter_yield)
#endif
csw 26mar2013
        call ipgeti ('eirene_iscltst', iscltst)

csw 26feb2013
        call ipgeti ('eirene_neutr_avg', neutr_avg)
        call ipgeti ('b2stbr_eir_src_nhist', nhist)
        call ipgeti ('eirene_underrelax', nrelax)
        call xertst (.not.(
     .    (neutr_avg.gt.1 .and. (nhist.gt.1 .or. nrelax.gt.1)) .or.
     .    (nhist.gt.1 .and. (neutr_avg.gt.1 .or. nrelax.gt.1)) .or.
     .    (nrelax.gt.1 .and. (neutr_avg.gt.1 .or. nhist.gt.1))),
     .   'Cannot have simultaneous multiple Eirene averaging schemes!')
        if(nhist.gt.1) call alloc_hist_arrays(nhist)
        if(nrelax.gt.1 .or. nhist.gt.1) call alloc_sum_arrays
        call xertst(ank_mods.eq.0.or.ionising_core.le.0,
     1   'contradictory inputs ank_mods and ionising_core!')
        if(eirene_mod.le.0) eirene_mod=1
        call ipgeti ('eirene_no_sseicorr', no_sseicorr)
        call ipgeti ('eirene_stab_corr', stab_corr)
      end if !}

      if(ncall.eq.0) then
        call write_f30(nx, ny)
        call get_nxt(nx, nxtl, nxtr)
        call get_jsep(nx,ny,jxi,jxa,jsep)
      endif
      if(ncall.eq.0.or.(main_call.and.mod(ncall,eirene_mod).eq.0)) then
        call copy_background(nx,ny,ns,BoRiS)
        if(ncall.eq.0.or.save_f31.ne.0) then
          call write_f31(nx, ny)
        endif
      endif

      if(ncall.eq.0) then !{
        call alloc_b2mod_neutrals_save(nstraid)
        call read_b2mod_neutrals_save(99,'    ',saved_volrec_used)
        call xertst(nstrat.le.nstraid,
     &   'Increase DEF_NSTRA in DIMENSIONS.F')
        if(maxval(saved_volrec).gt.0.0_R8) then !{
          volsumnold(1:nstrat) = saved_volrec(1:nstrat)
        else !}{
          do istrai = 1, nstrat !{
            if (CRCSTRA(istrai).EQ.'V') then !{
              if(volsumnold(istrai).lt.0.0_R8) then
                volsumnold(istrai)=volrecstart(istrai)
                saved_volrec_used =.true.
              endif
            else !}{
              volsumnold(istrai)=0.0_R8
            end if !}
          end do !}
        end if !}

c*** Check the b2espcr array
        hlp_bad_isp=.false.
        ispc_ref=-1
        do i=0,ns-1 !{
cwd exclude fluid neutrals from this check: their recycling coefficient does not
c   necessarily have to match that of the ions, since they are not treated in the
c   same Eirene stratum
          if (is_neutral(i)) cycle
          j=b2espcr(i)
          if(ispc_ref(j).lt.0) then !{
            ispc_ref(j)=i
            do k=1,nstrat !{
              if(use_recyceir.gt.0) then
                recyctmp=recyceir(k)
              else
                recyctmp=recyc(i,k)
              endif
              hlp_scp(j,k)=recyctmp
            end do !}
          else !}{
            hlp_tt=.false.
            do k=1,nstrat !{
              if(use_recyceir.gt.0) then
                recyctmp=recyceir(k)
              else
                recyctmp=recyc(i,k)
              endif
              hlp_tt=hlp_tt.or.
     1         abs(hlp_scp(j,k)-recyctmp).gt.rec_tol
            end do !}
            if(hlp_tt) then !{
              write (0,*)
     1           'eirene_mc: different values of recyceir ',
     ,                                  'species,is= ',j,i
              hlp_bad_isp=.true.
            end if !}
          end if !}
        end do !}
        call xertst(.not.hlp_bad_isp,
     1     'eirene_mc: error in recyceir specification')

csw 08feb2011 l_ddt definition moved down after Eirene call when it is
c             defined
c c*** "time layer"
c       l_ddt=0
c       if(crcstra(nstrat).eq.'T') l_ddt=nnlimi+nnstsi
        hlp_tt=.false.
      end if !}

      if (ncall.eq.0 .and. ionising_core.ne.0) then
cxpb compute core boundary surface area
        core_surface = 0.0_R8
        if (boundary_namelist.eq.1) then
          do ib = 1, nbc
            if (core_boundary(ib)) then
              do ix = 1, bc_list_size(ib)
                core_surface = core_surface +
     .            gs(topix(bc_list_x(ix,ib),bc_list_y(ix,ib)),
     .               topiy(bc_list_x(ix,ib),bc_list_y(ix,ib)),1)
              enddo
            endif
          enddo
        else
          do ix = -1, nx
            if (on_closed_surface(ix,0)) then
              core_surface = core_surface + gs(ix,0,1)
            endif
          enddo
        endif
cxpb compute ishigh and ion_name arrays
        ishigh = -1
        ion_name = ' '
        do j = 1, nnatmi
          is = eb2atcr(j)
          do js = is+1, ns-1
            if (lnext(is,js)) ishigh(j) = js
          enddo
          if (ishigh(j).eq.-1) ishigh(j) = is
          call species (ishigh(j), ion_name(j), .false.)
        enddo
      endif

*  -- COMPUTE TOTAL FLUXES TO RECYCLING BOUNDARIES --

c*** tflux(i)  : source strength for stratum i as given to Eirene
c*** tfl(i,j)  : the same for species j only
c*** flux_scale: conversion factor Eirene-B2
c*** recyc     : fudge factors for the sources in B2
cwdk only take outgoing fluxes into account in tflux and tfl, since
cwdk ingoing fluxes will be cut off on the EIRENE side
      na_feedback_on = .false.
      do is = 0, nspecies-1
        if (na_feedback_on) cycle
        if (na_feedback_actuator(is).eq.1)
     >   call xertst(nesepm_istra.gt.0.or.maxval(feedback_strata).gt.0,
     &    'No gas puff feedback strata set!')
        na_feedback_on = na_feedback_target(is).gt.0.0_R8.and.
     &                   na_feedback_actuator(is).eq.1.and.
     &                   feedback_strata(is).gt.0
        if (na_feedback_on) then
          na_feedback_on = na_feedback_on .and.
     &      species_start(feedback_strata(is)).le.eb2spcr(is+1).and.
     &      species_end(feedback_strata(is)).ge.eb2spcr(is+1)
        end if
      end do

      tfl=0.0_R8
      tflux=0.0_R8            !! the actual flux passed to EIRENE
      flux_scale=0.0_R8
      do istrai = 1, nstrat !{
        feedback_stratum=.false.
        do is=0,nspecies-1
          feedback_stratum=feedback_stratum.or.
     .                     istrai.eq.feedback_strata(is)
        end do
        flux_scale(istrai)=1.0_R8
        IF(crcstra(istrai).eq.'W') then !{
          do is = species_start(istrai),species_end(istrai) !{
            if(.not.is_neutral(is)) then !{
              do iy=1,rc_list_size(istrai) !{
                tflux(istrai)=tflux(istrai)-min(0.0_R8,
     1     fna(rightix(rc_list_x(iy,istrai),rc_list_y(iy,istrai)),
     2         rightiy(rc_list_x(iy,istrai),rc_list_y(iy,istrai)),0,is))
                tfl(istrai,b2eatcr(is))=tfl(istrai,b2eatcr(is))-
     .      min(0.0_R8,
     1     fna(rightix(rc_list_x(iy,istrai),rc_list_y(iy,istrai)),
     2         rightiy(rc_list_x(iy,istrai),rc_list_y(iy,istrai)),0,is))
              end do !}
            end if !}
          end do !}
        else if(crcstra(istrai).eq.'E') then !}{
          do is = species_start(istrai),species_end(istrai) !{
            if(.not.is_neutral(is)) then !{
              do iy=1,rc_list_size(istrai) !{
                tflux(istrai)=tflux(istrai)+max(0.0_R8,
     .           fna(rc_list_x(iy,istrai),rc_list_y(iy,istrai),0,is))
                tfl(istrai,b2eatcr(is))=tfl(istrai,b2eatcr(is))+
     .           max(0.0_R8,
     1           fna(rc_list_x(iy,istrai),rc_list_y(iy,istrai),0,is))
              end do !}
            end if !}
          end do !}
        else if(crcstra(istrai).eq.'S') then !}{
          do is = species_start(istrai),species_end(istrai) !{
            if(.not.is_neutral(is)) then !{
              do ix=1,rc_list_size(istrai) !{
                tflux(istrai)=tflux(istrai)-min(0.0_R8,
     1      fna(topix(rc_list_x(ix,istrai),rc_list_y(ix,istrai)),
     2          topiy(rc_list_x(ix,istrai),rc_list_y(ix,istrai)),1,is))
                tfl(istrai,b2eatcr(is))=tfl(istrai,b2eatcr(is))-
     .      min(0.0_R8,
     1      fna(topix(rc_list_x(ix,istrai),rc_list_y(ix,istrai)),
     2          topiy(rc_list_x(ix,istrai),rc_list_y(ix,istrai)),1,is))
              end do !}
            end if !}
          end do !}
        else if(crcstra(istrai).eq.'N') then !}{
          do is = species_start(istrai),species_end(istrai) !{
            if(.not.is_neutral(is)) then !{
              do ix=1,rc_list_size(istrai) !{
                tflux(istrai)=tflux(istrai)+max(0.0_R8,
     .           fna(rc_list_x(ix,istrai),rc_list_y(ix,istrai),1,is))
                tfl(istrai,b2eatcr(is))=tfl(istrai,b2eatcr(is))+
     .           max(0.0_R8,
     1           fna(rc_list_x(ix,istrai),rc_list_y(ix,istrai),1,is))
              end do !}
            end if !}
          end do !}
        else if(crcstra(istrai).eq.'C') then !}{
          if (ncall.ge.ncallfeedback.and.istrai.eq.nesepm_istra.and.
     &     (ndes.gt.0.0_R8.or.nesepm_pfr.gt.0.0_R8.or.
     &      private_flux_puff.gt.0.0_R8).and.lfeedback) then !{
            if (pfrregno1.ne.pfrregno2) then !{
              tflux(istrai)=(cbsna(0,isfeedback,pfrregno1)+
     1         cbsna(0,isfeedback,pfrregno2))/neutral_sources_rescale
            else !}{
              tflux(istrai)=
     1         cbsna(0,isfeedback,pfrregno1)/neutral_sources_rescale
            end if !}
            if(ncall.eq.ncallfeedback) then
              if(nesepm_pfr.gt.0.0_R8) then
                write(*,*) '"nesepm_pfr" feedback on stratum ',istrai
              else if(ndes.gt.0.0_R8) then
                write(*,*) '"ndes" feedback on stratum ',istrai
              else if(private_flux_puff.gt.0.0_R8) then
                write(*,*) '"private_flux_puff" feedback on stratum ',
     1           istrai
              endif
            endif
          else if (ncall.ge.ncallfeedback.and.nesepm_sol.gt.0.0_R8.and.
     1     istrai.eq.nesepm_istra.and.lfeedback) then !}{
            TFLUX(istrai)=
     1       (cbsna(0,isfeedback,solregno))/neutral_sources_rescale
            if(ncall.eq.ncallfeedback)
     1       write(*,*) '"nesepm_sol" feedback on stratum ',istrai
          else if (ncall.ge.ncallfeedback.and.volrec_sol.gt.0.0_R8.and.
     1     istrai.eq.nesepm_istra.and.lfeedback) then !}{
            TFLUX(istrai)=
     1       (cbsna(0,isfeedback,solregno))/neutral_sources_rescale
            if(ncall.eq.ncallfeedback)
     1       write(*,*) '"volrec_sol" feedback on stratum ',istrai
          else if (ncall.ge.ncallfeedback.and.ndes_sol.gt.0.0_R8.and.
     1     istrai.eq.nesepm_istra.and.lfeedback) then !}{
            TFLUX(istrai)=
     1       (cbsna(0,isfeedback,solregno))/neutral_sources_rescale
            if(ncall.eq.ncallfeedback)
     1       write(*,*) '"ndes_sol" feedback on stratum ',istrai
          else if (ncall.ge.ncallfeedback.and.nepedm_sol.gt.0.0_R8.and.
     1     istrai.eq.nesepm_istra.and.lfeedback) then !}{
            TFLUX(istrai)=
     1       (cbsna(0,isfeedback,solregno))/neutral_sources_rescale
            if(ncall.eq.ncallfeedback)
     1       write(*,*) '"nepedm_sol" feedback on stratum ',istrai
          else if (time_dep_puff(istrai)) then !}{
            if (time_dep_puff_func(istrai)) then
             call xertst(time_dep_puff_case(istrai).ge.0
     1               .and.time_dep_puff_case(istrai).le.5,
     1         'bad value for time_dep_puff_case')
             if (time_dep_puff_case(istrai) .eq. 0) then ! sinusoidal
               tflux(istrai) = time_dep_puff_param(1,istrai)
     &               +time_dep_puff_param(2,istrai)*(sin(2.d0*pi
     &               *time_dep_puff_param(3,istrai)*tim
     &               +time_dep_puff_param(4,istrai)))
             elseif (time_dep_puff_case(istrai) .eq. 1) then  ! stair+ramp
               if ((tim-INT(tim/2.d0/time_dep_puff_param(3,istrai))
     &            *2.*time_dep_puff_param(3,istrai)).ge.0. .and.
     &             (tim-INT(tim/2.d0/time_dep_puff_param(3,istrai))
     &            *2.*time_dep_puff_param(3,istrai)).lt.
     &             time_dep_puff_param(3,istrai)) then ! mod of tim/step size
                 tflux(istrai) = time_dep_puff_param(1,istrai)
     &                 +INT(tim/2.d0/time_dep_puff_param(3,istrai))
     &                 *time_dep_puff_param(2,istrai)
               else
                 tflux(istrai) = time_dep_puff_param(1,istrai)
     &                 +INT(tim/2.d0/time_dep_puff_param(3,istrai))
     &                 *time_dep_pufF_param(2,istrai)
     &                 +(tim-INT(tim/2.d0/time_dep_puff_param(3,istrai))
     &                 *2.*time_dep_puff_param(3,istrai)
     &                 -time_dep_puff_param(3,istrai))
     &                 *time_dep_puff_param(2,istrai)
     &                 /time_dep_puff_param(3,istrai)
               endif
             elseif (time_dep_puff_case(istrai) .eq. 2) then  ! stair+ramp
               tflux(istrai) = time_dep_puff_param(1,istrai)
     &               +tim*time_dep_puff_param(2,istrai)
             elseif (time_dep_puff_case(istrai) .eq. 3) then  ! linear + sinusoidal
               tflux(istrai) = time_dep_puff_param(1,istrai)
     &               + time_dep_puff_param(2,istrai)*(sin(2.d0*3.14159
     &               *time_dep_puff_param(3,istrai)*tim
     &               +time_dep_puff_param(4,istrai)))
     &               +tim*time_dep_puff_param(5,istrai)
             elseif (time_dep_puff_case(istrai) .eq. 4) then  ! step function
                 tflux(istrai) = time_dep_puff_param(1,istrai)
     &                 +INT(tim/time_dep_puff_param(3,istrai))
     &                 *time_dep_puff_param(2,istrai)
             elseif (time_dep_puff_case(istrai) .eq. 5) then  ! triangle wave
                 tflux(istrai) = (time_dep_puff_param(1,istrai)
     &                 -time_dep_puff_param(2,istrai))
     &                 +(4.d0*time_dep_puff_param(2,istrai)
     &                 /time_dep_puff_param(3,istrai))
     &                 * ((time_dep_puff_param(3,istrai)/2.d0)
     &                 - abs(mod((tim+time_dep_puff_param(3,istrai)/4.d0
     &                 +time_dep_puff_param(4,istrai))
     &                  ,(2.d0*(time_dep_puff_param(3,istrai)/2.d0)))
     &                 - time_dep_puff_param(3,istrai)/2.d0));
             endif
            else !if time_dep_puff_func
             if (tim.lt.gpdata(1,1,istrai)) then
              tflux(istrai)=userfluxparm(istrai,1)
             else if (tim.gt.gpdata(ngpdata(istrai),1,istrai)) then
              tflux(istrai)=gpdata(ngpdata(istrai),2,istrai)
             else
              ib=1
              it=ngpdata(istrai)
              do while (it-ib.gt.1)
                ih=(ib+it)/2
                if(tim.lt.gpdata(ih,1,istrai)) then
                  it=ih
                elseif(tim.gt.gpdata(ih,1,istrai)) then
                  ib=ih
                else
                  it=ih
                  ib=ih
                endif
              enddo
              if (ib.lt.it) then
                f=(tim-gpdata(ib,1,istrai))/
     1                (gpdata(it,1,istrai)-gpdata(ib,1,istrai))
                tflux(istrai)=f*
     1                (gpdata(it,2,istrai)-gpdata(ib,2,istrai))+
     2                 gpdata(ib,2,istrai)
              else
                tflux(istrai)=gpdata(ih,2,istrai)
              endif
             endif
           endif
          else if (ncall.ge.ncallfeedback.and.
     &             feedback_stratum.and.na_feedback_on) then !}{
            do is = 0, nspecies-1 !{
              if (feedback_strata(is).gt.0) then !{
                if (na_feedback_target(is).gt.0.0_R8.and.
     &              na_feedback_actuator(is).eq.1.and.
     &              species_start(feedback_strata(is)).le.
     &              eb2spcr(is+1)+nfluids(is+1).and.
     &            species_end(feedback_strata(is)).ge.eb2spcr(is+1))
     >           then !{
                  tflux(istrai)=tflux(istrai)+
     .             conpar(eb2spcr(is+1),na_feedback_ib(is),1)
                end if !}
              end if !}
            end do !}
          else !}{
            tflux(istrai)=userfluxparm(istrai,1)
          end if !}
          do is = species_start(istrai),species_end(istrai) !{
            tfl(istrai,b2eatcr(is))=1.0_R8
          end do !}
          if(tflux(istrai).gt.0.0_R8) then !{
            flux_scale(istrai)=1.0_R8/(tflux(istrai)*qe)
          else !}{
            flux_scale(istrai)=0.0_R8
          end if !}
        else if(crcstra(istrai).eq.'V') then !}{
          tflux(istrai)=1.0_R8/qe
          do is = species_start(istrai),species_end(istrai) !{
            tfl(istrai,b2eatcr(is))=1.0_R8
          end do !}
        else if(crcstra(istrai).eq.'T') then !}{
          tflux(istrai)=1.0_R8/qe
          do is = species_start(istrai),species_end(istrai) !{
            tfl(istrai,b2eatcr(is))=1.0_R8
          end do !}
c        else if(crcstra(istrai).eq.'U') then !}{
cc         tflux(istrai)=userfluxparm(istrai,1)
c          write(*,*) 'Option "U" not yet implemented'
c          tflux(istrai)=0.0_R8
c          do is = species_start(istrai),species_end(istrai) !{
c            tfl(istrai,b2eatcr(is))=1.0_R8
c          end do !}
c          if(tflux(istrai).gt.0.0_R8) then !{
cc warning: the following looks strange but is apparently correct
c             flux_scale(istrai)=1.0_R8/(tflux(istrai)*qe)
c          else !}{
c             flux_scale(istrai)=0.0_R8
c          end if !}
        else !}{
          WRITE(*,*) 'OPTION ',CRCSTRA(istrai),' NOT AVAILABLE'
        end if !}
      end do !}

      write(hlp_frm,'(a,i2,a)') '(a,1p,',nstrat,'(1x,1e14.7))'
      if(main_call.or.ncall.eq.0.or.output_style.gt.1) then !{
        write(*,hlp_frm) 'tflux: ',tflux
        write(*,hlp_frm) 'flux_scale: ',flux_scale
      end if !}

      last=.false.

cxpb !!!!!!!!!!!!!!!!!!!!!!!!!!! EIRENE CALL HERE !!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      eircount=0
      if(main_call .and. mod(ncall,eirene_mod).eq.0) then !{

csw 25jan2011 LTIME=.false. in eirsrt call
        my_ltime=.false.

        filename=fort_lc//'44'
        call find_file(filename,ex)
        call eirene_eirsrt(last,my_ltime,eirene_step_dt,tflux,
     *   b2brem,b2rad,b2qie(0),b2vdp,
     1   eirene_step_cpu)
#ifdef B25_EIRENE
        if(.not.ex) then
          do is=1,nnatmi
            if(streql(textan(is-1),' ')) textan(is-1) = texts(is+nsph)
          end do
          do is=1,nnmoli
            if(streql(textmn(is-1),' ')) textmn(is-1) = texts(is+nspa)
          end do
          do is=1,nnioni
            if(streql(textin(is-1),' ')) textin(is-1) = texts(is+nspam)
          end do
          call map_ion_components
        end if
#endif
        eircount=1
        if(ncall.eq.0.and.repeat_first_call.gt.0) then
          do while (eircount.le.repeat_first_call)
           call copy_background(nx,ny,ns,BoRiS)
csw 25jan2011 LTIME=.false. in eirsrt call
           my_ltime=.false.
           call eirene_eirsrt(last,my_ltime,eirene_step_dt,tflux,
     *       b2brem,b2rad,b2qie(0),b2vdp,
     1       eirene_step_cpu)
           eircount=eircount+1
          enddo
        endif
        if(ncall.eq.0) then
          write(*,*) 'Calling post_couple'
          call post_couple(nx,ny,ns,nstrat)
          write(*,*) 'Called post_couple'
          write(*,*) 'Checking strata correspondence'
          call strata_check
          write(*,*) 'Strata correspondence checked!'
        endif
      elseif(.not.main_call) then !}{
        if (ncall.lt.3)
     1   write(*,*) 'Skipping call to EIRENE because not main_call'
      elseif(mod(ncall,eirene_mod).ne.0) then !}{
        write(*,*) 'Skipping call to EIRENE because of mod condition'
      end if !}
cxpb !!!!!!!!!!!!!!!!!!!!!!!!!!! EIRENE CALL ENDED !!!!!!!!!!!!!!!!!!!!!!!!!!!!

! IYS 20.01.2017 to extend Eirene quantities to guard cells
#ifdef B25_EIRENE
      call extend(dab2(0:ndxp,0:ndyp,1:nnatmi,1),
     &                   ndxp-1,ndyp-1,natm,nx,ny,nnatmi)
      call extend(tab2(0:ndxp,0:ndyp,1:nnatmi,1),
     &                   ndxp-1,ndyp-1,natm,nx,ny,nnatmi)
      call extend(dmb2(0:ndxp,0:ndyp,1:nnmoli,1),
     &                   ndxp-1,ndyp-1,nmol,nx,ny,nnmoli)
      call extend(tmb2(0:ndxp,0:ndyp,1:nnmoli,1),
     &                   ndxp-1,ndyp-1,nmol,nx,ny,nnmoli)
      call extend(dib2(0:ndxp,0:ndyp,1:nnioni,1),
     &                   ndxp-1,ndyp-1,nion,nx,ny,nnioni)
      call extend(tib2(0:ndxp,0:ndyp,1:nnioni,1),
     &                   ndxp-1,ndyp-1,nion,nx,ny,nnioni)
      call extend(rfluxa(0:ndxp,0:ndyp,1:nnatmi,1),
     &                   ndxp-1,ndyp-1,natm,nx,ny,nnatmi)
      call extend(rfluxm(0:ndxp,0:ndyp,1:nnmoli,1),
     &                   ndxp-1,ndyp-1,nmol,nx,ny,nnmoli)
      call extend(pfluxa(0:ndxp,0:ndyp,1:nnatmi,1),
     &                   ndxp-1,ndyp-1,natm,nx,ny,nnatmi)
      call extend(pfluxm(0:ndxp,0:ndyp,1:nnmoli,1),
     &                   ndxp-1,ndyp-1,nmol,nx,ny,nnmoli)
      call extend(tfluxa(0:ndxp,0:ndyp,1:nnatmi,1),
     &                   ndxp-1,ndyp-1,natm,nx,ny,nnatmi)
      call extend(tfluxm(0:ndxp,0:ndyp,1:nnmoli,1),
     &                   ndxp-1,ndyp-1,nmol,nx,ny,nnmoli)
#endif
! IYS 20.01.2017 to extend Eirene quantities to guard cells


csw 08feb2011 put l_ddt definition here since nnlimi(nlim) is known only after Eirene call
      if(ank_mods.ne.0 .or. ionising_core.eq.1) then
c*** "time layer"
        l_ddt=0
        if(crcstra(nstrai).eq.'T') then
          l_ddt=nnlimi+nnstsi
csw 2011 in case of no time layer in Eirene
          if(.not.my_ltime) then
csw 2011
c           l_ddt=l_ddt+1
            l_ddt=0
          endif
        endif
      end if !}

      if(main_call.or.ncall.eq.0.or.output_style.gt.1) then !{
        write(*,hlp_frm) 'tflux: ',tflux
        write(*,hlp_frm) 'flux_scale: ',flux_scale
      end if !}

cxpb  Additional diagnostic for comparison with SOLPS4.0
      if (output_style.gt.1) then
        write(*,*) 'in b2mod_eirene: ncall = ', ncall
        write(*,*) 'in b2mod_eirene: main_call = ', main_call
      end if
      if (ncall.eq.0) then
        call alloc_b2mod_user(nx,ns,nlim,nmol,ntns)
        if (main_call) call b2trci(nx,ny,ns,nxtl,nxtr,jsep,nnreg(0))
        if (write_nml_user) call write_b2mod_user_namelist
      endif
      if (main_call.and.mod(ncall,eirene_mod).eq.0
     & .and.ank_tracing.ne.0.and.output_style.ge.2) then
        write(*,*) 'Eirene output in eirene_mc :'
        sum_wld(1:4,1:2+nnatmi)=0.0_R8
        do istrai = 1, nstrat
          do i=1,ntns
            sum_wld(ktns(i),1)=sum_wld(ktns(i),1)+
     .                         wldnek(nnlimi+ltns(i),istrai)
            sum_wld(ktns(i),2)=sum_wld(ktns(i),2)+
     .                         wldnep(nnlimi+ltns(i),istrai)
            do is = 1, nnatmi
              sum_wld(ktns(i),2+is)=sum_wld(ktns(i),2+is)+
     .                                  wldna(nnlimi+ltns(i),is,istrai)
            end do
          end do
        end do
        write(*,*) 'Summed target data : '
        write(hlp_frm,'(a,i2,a)') '(a,i1,a,a,',2+nnatmi,'e14.6)'
        do i = 1, 2*max(1,nncut)
          write(*,hlp_frm) 'Target ',i,' : ',
     .     'wldnek, wldnep, wldna = ', sum_wld(i,1:2+nnatmi)
        end do
      endif
cxpb

cxpb  Compute recombination sources rescaling factor if necessary
      do istrai = 1, nstrat !{
        if(volsumnold(istrai).gt.0.0_R8) then !{
          if(output_style.gt.1.or.main_call)
     .     write(*,'(a,i2,1p,2(1x,1e14.7))')
     .      'stratum, volsumnold, abs(volsumn) ',
     .      istrai, volsumnold(istrai), abs(volsumn(istrai))/qe
        end if !}
        if(volsumnold(istrai).gt.0.0_R8.and.volrecinc.gt.0.0_R8) then !{
          if(abs(volsumn(istrai)/qe).gt.
     &           volsumnold(istrai)*volrecinc) then !{
            if(main_call) then !{
              flux_scale(istrai)=volsumnold(istrai)*volrecinc/
     1         abs(volsumn(istrai)/qe)
              if(output_style.ge.1) write(*,'(a,i2,a,1p,g12.5)')
     .         'Rescaling volume recombination stratum ',
     1         istrai,' by ',flux_scale(istrai)
            else !}{
              flux_scale(istrai)=volsumnold(istrai)/
     .         abs(volsumn(istrai)/qe)
            end if !}
          end if !}
          volsumnold(istrai)=(1.0_R8-volrecwt)*volsumnold(istrai)+
     1     volrecwt*abs(volsumn(istrai)/qe)*flux_scale(istrai)
          saved_volrec(istrai) = volsumnold(istrai)
        end if !}
      end do !}
cxpb

      if(main_call.and.mod(ncall,eirene_mod).eq.0.and.ank_mods.ne.0)
     > then
        if(ncall.eq.0) then !{
          hlp_tt=.false.
          hlp_bad_isp=.false.
          do i=0,ns-1 !{
            hlp_tt=b2eatcr(i).lt.0 .or. b2eatcr(i).gt.nnatmi
            hlp_bad_isp= hlp_bad_isp .or. hlp_tt
            if(hlp_tt) write (0,*) 'is,b2eatcr,nnatmi: ',
     ,         i,b2eatcr(i),nnatmi
          end do !}
          call xertst (.not.hlp_bad_isp,
     1       'eirene_mc: wrong value found in b2eatcr')
        end if !}

        if(dbg_eir_mc.gt.9 .and.
     &       main_call .and. mod(ncall,eirene_mod).eq.0) then !{
          write(*,'()')
          print *,'%%% Eirene data in B2 %%%',nlimps
          print '(a,5(i3,1x))','nnatmi,nnmoli,nnplsi,nnstsi,nnlimi = ',
     ,       nnatmi,nnmoli,nnplsi,nnstsi,nnlimi
          write(hlp_frm1,'(a,i2,a)')
     ,     '(1x,',2*nnatmi+2*nnmoli,'(6x,a2,i3.2))'
          write(hlp_frm2,'(a,i2,a)')
     ,     '(1p,i4,',2*nnatmi+2*nnmoli,'e11.4)'
          write(hlp_frm3,'(a,i2,a)')
     ,     '(1p,i4,',nnplsi,'e11.4)'
          write(hlp_frm4,'(a,i2,a)')
     ,     '(1p,i4,',1+nnatmi+nnmoli,'e11.4)'
          do k=0,nstrat !{
            write(*,'()')
            if (nnmoli.gt.0) then
              print *,'wldna, wldnm, wldra, wldrm:',k
            else
              print *,'wldna, wldra:',k
            end if
            write (*,hlp_frm1)
     ,         ('na',i,i=1,nnatmi),('nm',i,i=1,nnmoli),
     ,         ('ra',i,i=1,nnatmi),('rm',i,i=1,nnmoli)
            do j=1,nnlimi+nnstsi !{
              if(j.le.nnlimi .or. j.gt.nnlimi) then !{
                hlp_pr=.false.
                do i=1,nnatmi !{
                  hlp_pr= hlp_pr .or. wldna(j,i,k).ne.0.0_R8
     .               .or. wldra(j,i,k).ne.0.0_R8
                end do !}
                do i=1,nnmoli !{
                  hlp_pr= hlp_pr .or. wldnm(j,i,k).ne.0.0_R8
     .               .or. wldrm(j,i,k).ne.0.0_R8
                end do !}
                if(hlp_pr) then !{
                  write(*,hlp_frm2) j,
     ,               (wldna(j,i,k),i=1,nnatmi),
     ,               (wldnm(j,i,k),i=1,nnmoli),
     ,               (wldra(j,i,k),i=1,nnatmi),
     ,               (wldrm(j,i,k),i=1,nnmoli)
                end if !}
              end if !}
            end do !}

            write(*,'()')
            print *,'wldpp:',k
            do j=nnlimi+1,nnlimi+nnstsi !{
              hlp_pr=.false.
              do i=1,nnplsi !{
                hlp_pr = hlp_pr .or. wldpp(j,i,k).ne.0.0_R8
              end do !}
              if (hlp_pr)
     >         write(*,hlp_frm3) j,(wldpp(j,i,k),i=1,nnplsi)
            end do !}
            write(*,'()')
            if (nnmoli.gt.0) then
              print *,'wldpeb, wldpa, wldpm:',k
            else
              print *,'wldpeb, wldpa:',k
            end if
            do j=nnlimi+1,nnlimi+nnstsi !{
              hlp_pr=wldpeb(j,k).ne.0.0_R8
              do i=1,nnatmi !{
                hlp_pr= hlp_pr .or. wldpa(j,i,k).ne.0.0_R8
              end do !}
              do i=1,nnmoli !{
                hlp_pr= hlp_pr .or. wldpm(j,i,k).ne.0.0_R8
              end do !}
              if (hlp_pr) write(*,hlp_frm4) j,wldpeb(j,k),
     ,           (wldpa(j,i,k),i=1,nnatmi),
     ,           (wldpm(j,i,k),i=1,nnmoli)
            end do !}
          end do !}
          print '(a)','%%%'
        end if !}
        if(dbg_eir_mc.gt.0 .and.
     &       main_call .and. mod(ncall,eirene_mod).eq.0) then !{
c*** Print the particle fluxes for each stratum and species
          print '(a)',
     ,     '%%%=============== total fluxes by stratum and '//
     ,     'species before correction ===============%%%'
          print '(a25,20(4x,a5,4x))',' Stratum Species lstrascl',
     ,       ' tfl ','wldpi','wldpo','wldni','wldno',
     ,       'Srcpl','puff ','d/dt ','crsrc','error'
          do istrai=1,nstrat !{
            do ls=1,nnatmi !{
              hlp_pi(ls)=0.0_R8
              hlp_po(ls)=0.0_R8
              hlp_ni(ls)=0.0_R8
              hlp_no(ls)=0.0_R8
              hlp_g(ls) =0.0_R8
              hlp_s(ls) =0.0_R8
              hlp_t(ls) =0.0_R8
            end do !}
            do j=1,nnlimi+nnstsi !{
              if(j.gt.nnlimi .and. j.ne.l_ddt) then !{
                do is=0,ns-1 !{
                  if(.not.is_neutral(is)) hlp_pi(b2eatcr(is))=
     +             hlp_pi(b2eatcr(is))+wldpp(j,fluids_list(is),istrai)
                end do !}
                do i=1,nnatmi !{
                  hlp_po(i)=hlp_po(i)+wldpa(j,i,istrai)
                  do is=1,nnmoli !{
                    s=mlcmp(i,is)
                    if(s.gt.0.0_R8)
     -                 hlp_po(i)=hlp_po(i)+wldpm(j,is,istrai)*s
                  end do !}
                end do !}
              end if !}
              hlp_log=j.le.nnlimi
#ifdef B25_EIRENE
     &                .and. iliin(j) .gt. 0                ! IYS 03.12.2019 Consider non-transparent surfaces only,
                                                           ! transparent ones should not contribute into balance
#endif
              if(hlp_log .or. (j.gt.nnlimi .and. j.ne.l_ddt)) then !{
                if(.not.hlp_log) then !{
                  do i=1,ncns !{
                    hlp_log=hlp_log .or. j.eq.lcns(i)+nnlimi
                  end do !}
                  do i=1,ntns !{
                    hlp_log=hlp_log .or. j.eq.ltns(i)+nnlimi
                  end do !}
                end if !}

                if(hlp_log) then !{
                  do i=1,nnatmi !{
                    hlp_ni(i)=hlp_ni(i)+wldna(j,i,istrai)
                    hlp_no(i)=hlp_no(i)+wldra(j,i,istrai)
                    do is=1,nnmoli !{
                      s=mlcmp(i,is)
                      hlp_ni(i)=hlp_ni(i)+wldnm(j,is,istrai)*s
                      hlp_no(i)=hlp_no(i)+wldrm(j,is,istrai)*s
                    end do !}
                  end do !}
                end if !}
              else if(j.eq.l_ddt) then !}{
                do i=1,nnatmi !{
                  hlp_t(i)=hlp_t(i)+wldna(j,i,istrai)
                  do is=1,nnmoli !{
                    s=mlcmp(i,is)
                    hlp_t(i)=hlp_t(i)+wldnm(j,is,istrai)*s
                  end do !}
                end do !}
              end if !}
            end do !}

            do ix=1,nx !{
              do iy=1,ny !{
                if(resignore(ix-1,iy-1,1).ne.0) then !{
                  do is=0,ns-1 !{
                    if(.not.is_neutral(is)) then !{
                      hlp_s(b2eatcr(is))=hlp_s(b2eatcr(is))+
     +                   sni(ix,iy,fluids_list(is),istrai)
                    end if !}
                  end do !}
                end if !}
              end do !}
            end do !}

            do ls=1,nnatmi !{
              kks=ispc_ref(ls)
              hlp_s(ls)=hlp_s(ls)*tflux(istrai)*flux_scale(istrai)
              if(hlp_pi(ls).gt.0.0_R8) then !{
                t0=tfl(istrai,ls)/hlp_pi(ls)
              else !}{
                t0=1.0_R8
              end if !}
              if(use_recyceir.gt.0) then
                recyctmp=recyceir(istrai)
              else
                recyctmp=recyc(kks,istrai)
              endif
              if(crcstra(istrai).eq.'C')
     -           hlp_g(ls)=tflux(istrai)*
     *           gpfc(ls,istrai)*recyctmp
              s=hlp_po(ls)+hlp_g(ls)+hlp_no(ls)-
     -           (hlp_ni(ls)+hlp_s(ls)+hlp_t(ls))
              print '(i6,2i8,3x,1p,10e13.5)',istrai,ls,
     ,           lstrascl(istrai,ls),tfl(istrai,ls),
     ,           hlp_pi(ls),hlp_po(ls),hlp_ni(ls),hlp_no(ls),
     ,           hlp_s(ls),hlp_g(ls),hlp_t(ls),t0,s
            end do !}
          end do !}
          print '(a)',
     ,     '%%%============================================'//
     ,     '=========================================%%%'
        end if !}
      endif

c carlis adding gas puff output for particle balance routine
      do istrai=1,nstrat
        igpuff(istrai)=0
        if (crcstra(istrai).eq.'C') then
          igpuff(istrai)=1
        end if
      end do
      print *,'Number of gas puff sources: ',sum(igpuff)
      if (sum(igpuff).gt.0) then
        do istrai=1,nstrat
          do ls=1,nnatmi
            if (igpuff(istrai).gt.0.0_R8 .and.
     .          gpfc(ls,istrai).gt.0.0_R8) then
              if(use_recyceir.gt.0) then
                recyctmp=recyceir(istrai)
              else
                kks=ispc_ref(ls)
                recyctmp=recyc(kks,istrai)
              endif
              hlp_g(ls)=tflux(istrai)*
     *           gpfc(ls,istrai)*recyctmp
              print *,'Gas puff strength, atom: ',hlp_g(ls),ls
            end if
          end do
        end do
      end if

c*** Re-normalisation of the sources and fluxes in order to enforce
c*** the particle balance for recycling and recombining species

csw 23nov2012
      hlp_tt=.false.

csw 2012 fixed rescaling procedure
c        (otherwise rescaling procedure does not correctly work with internal iterations nstg2>1)
c        this is phase 1 of rescaling process
c
      if(main_call.and.mod(ncall,eirene_mod).eq.0
     .            .and.ank_mods.ne.0) then
        write(hlp_frm1,'(a,i2,a)') '(a,1p,i3,',nnatmi,'e13.5)'
        do istrai=1,nstrat !{
          if((crcstra(istrai).eq.'W' .or.
     .        crcstra(istrai).eq.'E' .or.
     .        crcstra(istrai).eq.'N' .or.
     .        crcstra(istrai).eq.'S')) then

c*** Recycling stratum
c*** First, correct the balance between B2 and Eirene sources for the
c*** species to be scaled with themselves (recycling species)

            do ls=1,nnatmi !{
              hlp_pi(ls)=0.0_R8
              hlp_po(ls)=0.0_R8
              hlp_ni(ls)=0.0_R8
              hlp_no(ls)=0.0_R8
              hlp_g(ls) =0.0_R8
              hlp_s(ls) =0.0_R8
              hlp_t(ls) =0.0_R8
            end do !}
            n=nnstsi
            if(l_ddt.gt.0) n=n-1
            do i=1,n !{
              j=nnlimi+i
              do is=0,ns-1 !{
                if(.not.is_neutral(is)) then !{
                  ls=b2eatcr(is)
                  if(lstrascl(istrai,ls).eq.ls) then !{
                    hlp_pi(ls)=hlp_pi(ls)+
     .               wldpp(j,fluids_list(is),istrai)
                  end if !}
                end if !}
              end do !}
            end do !}
            do ls=1,nnatmi !{
              if(hlp_pi(ls).gt.0.0_R8) then !{
                hlp_scl(ls)=tfl(istrai,ls)/hlp_pi(ls)
              else !}{
                hlp_scl(ls)=1.0_R8
              end if !}
            end do !}
c!###{
            if(dbg_eir_mc.gt.0) then
              write(*,hlp_frm1) ' @@@ istrai,pi: ',
     ,                            istrai,(hlp_pi(ls),ls=1,nnatmi)
              write(*,hlp_frm1) ' @@@       tfl: ',
     ,                        istrai,(tfl(istrai,ls),ls=1,nnatmi)
              write(*,hlp_frm1) ' @@@       scl: ',
     ,                           istrai,(hlp_scl(ls),ls=1,nnatmi)
            endif
c!###}

            do ix=1,nx !{
              do iy=1,ny !{
                do is=0,ns-1 !{
                  ls=b2eatcr(is)
                  if(.not.is_neutral(is)) then !{
                    if(lstrascl(istrai,ls).gt.0) then !{
                      sni(ix,iy,fluids_list(is),istrai)=
     *                   sni(ix,iy,fluids_list(is),istrai)*hlp_scl(ls)
                      smo(ix,iy,fluids_list(is),istrai)=
     *                   smo(ix,iy,fluids_list(is),istrai)*hlp_scl(ls)
                      if (balance_netcdf.ne.0) then !djm Jan2017
                        sni_papl(ix,iy,fluids_list(is),istrai)=
     *                     sni_papl(ix,iy,fluids_list(is),istrai)*
     *                     hlp_scl(ls)
                        sni_pmpl(ix,iy,fluids_list(is),istrai)=
     *                     sni_pmpl(ix,iy,fluids_list(is),istrai)*
     *                     hlp_scl(ls)
                        sni_pipl(ix,iy,fluids_list(is),istrai)=
     *                     sni_pipl(ix,iy,fluids_list(is),istrai)*
     *                     hlp_scl(ls)
                        sni_pppl(ix,iy,fluids_list(is),istrai)=
     *                     sni_pppl(ix,iy,fluids_list(is),istrai)*
     *                     hlp_scl(ls)
                        smo_mapl(ix,iy,fluids_list(is),istrai)=
     *                     smo_mapl(ix,iy,fluids_list(is),istrai)*
     *                     hlp_scl(ls)
                        smo_mmpl(ix,iy,fluids_list(is),istrai)=
     *                     smo_mmpl(ix,iy,fluids_list(is),istrai)*
     *                     hlp_scl(ls)
                        smo_mipl(ix,iy,fluids_list(is),istrai)=
     *                     smo_mipl(ix,iy,fluids_list(is),istrai)*
     *                     hlp_scl(ls)
                        smo_cppv(ix,iy,fluids_list(is),istrai)=
     *                     smo_cppv(ix,iy,fluids_list(is),istrai)*
     *                     hlp_scl(ls)
                      endif
                    end if !}
                    if(resignore(ix-1,iy-1,1).ne.0) then !{
                      hlp_s(ls)=hlp_s(ls)+
     .                   sni(ix,iy,fluids_list(is),istrai)
                    end if !}
                  else
                    sna_paat(ix,iy,ls,istrai)=
     *                 sna_paat(ix,iy,ls,istrai)*hlp_scl(ls)
                    sna_pmat(ix,iy,ls,istrai)=
     *                 sna_pmat(ix,iy,ls,istrai)*hlp_scl(ls)
                    sna_piat(ix,iy,ls,istrai)=
     *                 sna_piat(ix,iy,ls,istrai)*hlp_scl(ls)
                  end if !}
                end do !}
                do is = 1, nnmoli
                  ls=lmolscl(is)
                  sna_paml(ix,iy,is,istrai)=
     *               sna_paml(ix,iy,is,istrai)*hlp_scl(ls)
                  sna_pmml(ix,iy,is,istrai)=
     *               sna_pmml(ix,iy,is,istrai)*hlp_scl(ls)
                  sna_piml(ix,iy,is,istrai)=
     *               sna_piml(ix,iy,is,istrai)*hlp_scl(ls)
                end do
                do is = 1, nnioni
                  ls=lionscl(is)
                  sna_paio(ix,iy,is,istrai)=
     *               sna_paio(ix,iy,is,istrai)*hlp_scl(ls)
                  sna_pmio(ix,iy,is,istrai)=
     *               sna_pmio(ix,iy,is,istrai)*hlp_scl(ls)
                  sna_piio(ix,iy,is,istrai)=
     *               sna_piio(ix,iy,is,istrai)*hlp_scl(ls)
                end do
              end do !}
            end do !}
            do ls=1,nnatmi !{
              hlp_s(ls)=hlp_s(ls)*tflux(istrai)*flux_scale(istrai)
            end do !}

            do j=1,nnlimi+nnstsi !{
              if(j.gt.nnlimi .and. j.ne.l_ddt) then !{
                do is=0,ns-1 !{
                  if(.not.is_neutral(is)) then !{
                    ls=b2eatcr(is)
                    if(lstrascl(istrai,ls).gt.0) then !{
                      wldpp(j,fluids_list(is),istrai)=
     =                   wldpp(j,fluids_list(is),istrai)*hlp_scl(ls)
                    end if !}
                  end if !}
                end do !}
                do ls=1,nnatmi !{
                  wldpa(j,ls,istrai)=wldpa(j,ls,istrai)*hlp_scl(ls)
                  hlp_po(ls)=hlp_po(ls)+wldpa(j,ls,istrai)
                end do !}

                do i=1,nnmoli !{
                  ls=lmolscl(i)
                  if(ls.gt.0 .and. ls.le.nnatmi) then !{
                    wldpm(j,i,istrai)=wldpm(j,i,istrai)*hlp_scl(ls)
                  end if !}
                  do ls=1,nnatmi !{
                    hlp_po(ls)=hlp_po(ls)+
     +               wldpm(j,i,istrai)*mlcmp(ls,i)
                  end do !}
                end do !}
              end if !}

              hlp_log=j.le.nnlimi
#ifdef B25_EIRENE
     &                .and. iliin(j) .gt. 0                ! IYS 03.12.2019 Consider non-transparent surfaces only,
                                                           ! transparent ones should not contribute into balance
#endif
              if(.not.hlp_log) then !{
                do i=1,ncns !{
                  hlp_log=hlp_log .or. j.eq.lcns(i)+nnlimi
                end do !}
                if(.not.hlp_log) then !{
                  do i=1,ntns !{
                    hlp_log=hlp_log .or. j.eq.ltns(i)+nnlimi
                  end do !}
                end if !}
              end if !}

              do ls=1,nnatmi !{
                if(lstrascl(istrai,ls).gt.0) then !{
                  wldna(j,ls,istrai)=wldna(j,ls,istrai)*hlp_scl(ls)
                  wldra(j,ls,istrai)=wldra(j,ls,istrai)*hlp_scl(ls)
                end if !}
c!###{
c               if(hlp_tt) print '(a,1p,3i4,l3,4e13.5)',
c     ,          ' @@@ istrai,j,ls,log,scl,wldpa,wldna,wldra',
c     ,          istrai,j,ls,hlp_log,hlp_scl(ls),wldpa(j,ls,istrai),
c     ,                       wldna(j,ls,istrai),wldra(j,ls,istrai)
c!###}
                if(hlp_log) then !{
                  hlp_ni(ls)=hlp_ni(ls)+wldna(j,ls,istrai)
                  hlp_no(ls)=hlp_no(ls)+wldra(j,ls,istrai)
                else if(j.eq.l_ddt) then !}{
                  hlp_t(ls)=hlp_t(ls)+wldna(j,ls,istrai)
                end if !}
              end do !}

              do i=1,nnmoli !{
                ls=lmolscl(i)
                if(ls.gt.0 .and. ls.le.nnatmi) then !{
                  wldnm(j,i,istrai)=wldnm(j,i,istrai)*hlp_scl(ls)
                  wldrm(j,i,istrai)=wldrm(j,i,istrai)*hlp_scl(ls)
                end if !}
                do ls=1,nnatmi !{
                  if(hlp_log) then !{
                    hlp_ni(ls)=hlp_ni(ls)+
     +               wldnm(j,i,istrai)*mlcmp(ls,i)
                    hlp_no(ls)=hlp_no(ls)+
     +               wldrm(j,i,istrai)*mlcmp(ls,i)
                  else if(j.eq.l_ddt) then !}{
                    hlp_t(ls)=hlp_t(ls)+
     +               wldnm(j,i,istrai)*mlcmp(ls,i)
                  end if !}
                end do !}
              end do !}
            end do !}
c!###{
            if(dbg_eir_mc.gt.0) then
              do ls=1,nnatmi !{
                s=hlp_po(ls)+hlp_no(ls)-hlp_s(ls)-hlp_ni(ls)
                print '(i6,2i8,3x,1p,8e13.5)',istrai,ls,
     ,                    lstrascl(istrai,ls),tfl(istrai,ls),
     ,                    hlp_pi(ls),hlp_po(ls),hlp_ni(ls),
     ,                    hlp_no(ls),hlp_s(ls),hlp_scl(ls),s
              end do !}
            endif
c!###}

c*** Now rectify the Eirene particle balance for each species -
c*** simply correcting the neutral sources and sinks

            do ls=1,nnatmi !{
c<<<
              if(dbg_eir_mc.gt.0) then
                print *,'@@@ istrai,ls,lstrascl: ',
     ,                       istrai,ls,lstrascl(istrai,ls)
              endif
c>>>
              hlp_scl(ls)=1.0_R8
              if(lstrascl(istrai,ls).eq.ls) then !{

c*** Consider flux of neutrals from plasma recombination on the target
c*** fixed - and adjust other terms to keep the balance

                s=hlp_s(ls)+hlp_t(ls)+hlp_ni(ls)-hlp_no(ls)
                if(s.gt.0.0_R8) hlp_scl(ls)=hlp_po(ls)/s

              else if(lstrascl(istrai,ls).gt.0) then

*** This species in this stratum is born from some interactions of
c*** other species. Try to scale the fluxes to keep the same source
c*** terms in balance.

                s=hlp_po(ls)+hlp_t(ls)+hlp_no(ls)-hlp_ni(ls)
               if(abs(s).gt.0.0_R8 .and.
     .            sign(1.0_R8,s).eq.sign(1.0_R8,hlp_s(ls)))
     >               hlp_scl(ls)=hlp_s(ls)/s
c<<<
               if(dbg_eir_mc.gt.0) then
                 print '(1x,a,1p,5e13.5)','@@@ scl,s,po,no,ni,t: ',
     ,          hlp_scl(ls),s,hlp_po(ls),hlp_no(ls),hlp_ni(ls),hlp_t(ls)
               endif
c>>>

              end if !}
              hlp_scp(ls,istrai)=hlp_scl(ls)
            end do !}

            do ix=1,nx !{
              do iy=1,ny !{
                do is=0,ns-1 !{
                  ls=b2eatcr(is)
                  if(.not.is_neutral(is)) then !{
                    if(lstrascl(istrai,ls).gt.0) then !{
                      sni(ix,iy,fluids_list(is),istrai)=
     *                 sni(ix,iy,fluids_list(is),istrai)*hlp_scl(ls)
                      smo(ix,iy,fluids_list(is),istrai)=
     *                 smo(ix,iy,fluids_list(is),istrai)*hlp_scl(ls)
                      if (balance_netcdf.ne.0) then !djm Jan2017
                        sni_papl(ix,iy,fluids_list(is),istrai)=
     *                   sni_papl(ix,iy,fluids_list(is),istrai)*
     *                   hlp_scl(ls)
                        sni_pmpl(ix,iy,fluids_list(is),istrai)=
     *                   sni_pmpl(ix,iy,fluids_list(is),istrai)*
     *                   hlp_scl(ls)
                        sni_pipl(ix,iy,fluids_list(is),istrai)=
     *                   sni_pipl(ix,iy,fluids_list(is),istrai)*
     *                   hlp_scl(ls)
                        sni_pppl(ix,iy,fluids_list(is),istrai)=
     *                   sni_pppl(ix,iy,fluids_list(is),istrai)*
     *                   hlp_scl(ls)
                        smo_mapl(ix,iy,fluids_list(is),istrai)=
     *                   smo_mapl(ix,iy,fluids_list(is),istrai)*
     *                   hlp_scl(ls)
                        smo_mmpl(ix,iy,fluids_list(is),istrai)=
     *                   smo_mmpl(ix,iy,fluids_list(is),istrai)*
     *                   hlp_scl(ls)
                        smo_mipl(ix,iy,fluids_list(is),istrai)=
     *                   smo_mipl(ix,iy,fluids_list(is),istrai)*
     *                   hlp_scl(ls)
                        smo_cppv(ix,iy,fluids_list(is),istrai)=
     *                   smo_cppv(ix,iy,fluids_list(is),istrai)*
     *                   hlp_scl(ls)
                      endif
                    end if !}
                  else !}{
                    sna_paat(ix,iy,ls,istrai)=
     *                 sna_paat(ix,iy,ls,istrai)*hlp_scl(ls)
                    sna_pmat(ix,iy,ls,istrai)=
     *                 sna_pmat(ix,iy,ls,istrai)*hlp_scl(ls)
                    sna_piat(ix,iy,ls,istrai)=
     *                 sna_piat(ix,iy,ls,istrai)*hlp_scl(ls)
                  end if !}
                end do !}
                do is = 1, nnmoli
                  ls = lmolscl(is)
                  sna_paml(ix,iy,is,istrai)=
     *               sna_paml(ix,iy,is,istrai)*hlp_scl(ls)
                  sna_pmml(ix,iy,is,istrai)=
     *               sna_pmml(ix,iy,is,istrai)*hlp_scl(ls)
                  sna_piml(ix,iy,is,istrai)=
     *               sna_piml(ix,iy,is,istrai)*hlp_scl(ls)
                end do
                do is = 1, nnioni
                  ls = lionscl(is)
                  sna_paio(ix,iy,is,istrai)=
     *               sna_paio(ix,iy,is,istrai)*hlp_scl(ls)
                  sna_pmio(ix,iy,is,istrai)=
     *               sna_pmio(ix,iy,is,istrai)*hlp_scl(ls)
                  sna_piio(ix,iy,is,istrai)=
     *               sna_piio(ix,iy,is,istrai)*hlp_scl(ls)
                end do
              end do !}
            end do !}

            do j=1,nnlimi+nnstsi !{
              do ls=1,nnatmi !{
                wldna(j,ls,istrai)=wldna(j,ls,istrai)*hlp_scl(ls)
                wldra(j,ls,istrai)=wldra(j,ls,istrai)*hlp_scl(ls)
                if(lstrascl(istrai,ls).ne.ls) then !{
                  wldpa(j,ls,istrai)=wldpa(j,ls,istrai)*hlp_scl(ls)
                end if !}
              end do !}
              do i=1,nnmoli !{
                ls=lmolscl(i)
                if(ls.gt.0 .and. ls.le.nnatmi) then !{
                  wldnm(j,i,istrai)=wldnm(j,i,istrai)*hlp_scl(ls)
                  wldrm(j,i,istrai)=wldrm(j,i,istrai)*hlp_scl(ls)
                  if(lstrascl(istrai,ls).ne.ls) then !{
                    wldpm(j,i,istrai)=wldpm(j,i,istrai)*hlp_scl(ls)
                  end if !}
                end if !}
              end do !}
            end do !}
          else if(crcstra(istrai).eq.'V' .or.
     .            crcstra(istrai).eq.'C') then !}{

c*** Recombination or gas puffing stratum.
c*** Here only the second pass of correction

            do ls=1,nnatmi !{
              hlp_pi(ls)=0.0_R8
              hlp_po(ls)=0.0_R8
              hlp_ni(ls)=0.0_R8
              hlp_no(ls)=0.0_R8
              hlp_g(ls) =0.0_R8
              hlp_s(ls) =0.0_R8
              hlp_t(ls) =0.0_R8
            end do !}

            do ix=1,nx !{
              do iy=1,ny !{
                if(resignore(ix-1,iy-1,1).ne.0) then !{
                  do is=0,ns-1 !{
                    if(.not.is_neutral(is)) then !{
                      ls=b2eatcr(is)
                      hlp_s(ls)=hlp_s(ls)+
     .                 sni(ix,iy,fluids_list(is),istrai)
                    end if !}
                  end do !}
                end if !}
              end do !}
            end do !}
            do ls=1,nnatmi !{
              hlp_s(ls)=hlp_s(ls)*tflux(istrai)*flux_scale(istrai)
            end do !}

            do j=1,nnlimi+nnstsi !{

              hlp_log=j.le.nnlimi
#ifdef B25_EIRENE
     &                .and. iliin(j) .gt. 0                ! IYS 03.12.2019 Consider non-transparent surfaces only,
                                                           ! transparent ones should not contribute into balance
#endif
              if(.not.hlp_log) then !{
                do i=1,ncns !{
                  hlp_log=hlp_log .or. j.eq.lcns(i)+nnlimi
                end do !}
                do i=1,ntns !{
                  hlp_log=hlp_log .or. j.eq.ltns(i)+nnlimi
                end do !}
              end if !}

              do ls=1,nnatmi !{
                if(hlp_log) then !{
                  hlp_ni(ls)=hlp_ni(ls)+wldna(j,ls,istrai)
                  hlp_no(ls)=hlp_no(ls)+wldra(j,ls,istrai)
                else if(j.eq.l_ddt) then !}{
                  hlp_t(ls)=hlp_t(ls)+wldna(j,ls,istrai)
                end if !}
              end do !}
              do i=1,nnmoli !{
                ls=lmolscl(i)
                if(ls.gt.0 .and. ls.le.nnatmi) then !{
                  if(hlp_log) then !{
                    hlp_ni(ls)=hlp_ni(ls)+
     +               wldnm(j,i,istrai)*mlcmp(ls,i)
                    hlp_no(ls)=hlp_no(ls)+
     +               wldrm(j,i,istrai)*mlcmp(ls,i)
                  else if(j.eq.l_ddt) then !}{
                    hlp_t(ls)=hlp_t(ls)+
     +               wldnm(j,i,istrai)*mlcmp(ls,i)
                  end if !}
                end if !}
              end do !}
            end do !}

            if(crcstra(istrai).eq.'C') then !{
              do ls=1,nnatmi !{
                kks=ispc_ref(ls)
                if(use_recyceir.gt.0) then
                  recyctmp=recyceir(istrai)
                else
                  recyctmp=recyc(kks,istrai)
                endif
                hlp_g(ls)=tflux(istrai)*gpfc(ls,istrai)*
     *             recyctmp
              end do !}
            end if !}

            do ls=1,nnatmi !{
              hlp_scl(ls)=1.0_R8
              if(lstrascl(istrai,ls).eq.ls) then !{
                s=hlp_s(ls)+hlp_t(ls)
                if(s.eq.0.0_R8) then !{
                  if(hlp_g(ls).gt.0.0_R8) then !{
                    hlp_scl(ls)=(hlp_ni(ls)-hlp_no(ls))/hlp_g(ls)
                  else !}{
                    if(hlp_ni(ls)-hlp_no(ls).ne.0.0_R8)
     -               write(*,'(a,a,a,i3,1x,i3)')
     ,                'eirene_mc: error in scaling - no ',
     ,                'sources and nonzero wall flux. ',
     ,                'stratum, species = ',istrai,ls
                  end if !}
                else !}{
                  s=hlp_g(ls)**2-4.0_R8*(hlp_ni(ls)-hlp_no(ls))*s
                  if(s.gt.0.0_R8) then !{
                    s=sqrt(s)
                    x1=0.5_R8*(hlp_g(ls)+s)/(hlp_s(ls)+hlp_t(ls))
                    x2=0.5_R8*(hlp_g(ls)-s)/(hlp_s(ls)+hlp_t(ls))
                    if(x2.le.0.0_R8) then !{
                      hlp_scl(ls)=x1
                    else if(x1.le.0.0_R8) then !}{
                      hlp_scl(ls)=x2
                    else !}{
                      if(abs(log(x1)).le.abs(log(x2))) then !{
                        hlp_scl(ls)=x1
                      else !}{
                        hlp_scl(ls)=x2
                      end if !}
                    end if !}
                    print '(a,i3,1p,8(1x,e12.5))',
     .                'ank+xpb: ls,s,g,t,ni,no,x1,x2,scl',ls,
     1                 hlp_s(ls),hlp_g(ls),hlp_t(ls),
     .                 hlp_ni(ls),hlp_no(ls),
     2                 x1,x2,hlp_scl(ls)
                  else !}{
                    write(*,'(a,a,i3,1x,i3)')
     ,               'eirene_mc: error in scaling - zero ',
     ,               'discriminant.  stratum, species = ',istrai,ls
                  end if !}
                end if !}
              end if !}
              if(hlp_scl(ls).le.0.0_R8) then !{
                write(*,'(a,a,a,i3,1x,i3,1x,1p,e11.4)')
     ,           'eirene_mc: error in scaling - non-',
     ,           'positive scaling factor. stratum, species, ',
     ,           'factor = ',istrai,ls,hlp_scl(ls)
                hlp_scl(ls)=1.0_R8
              else if(abs(log(hlp_scl(ls))).gt.
     1         abs(log(neut_scl_lim))) then !}{
                write(*,'(a,a,a,i3,1x,i3,1x,1p,e11.4)')
     ,           'eirene_mc: error in scaling - too ',
     ,           'strong variation. stratum, species, ',
     ,           'factor = ',istrai,ls,hlp_scl(ls)
                hlp_scl(ls)=1.0_R8
              end if !}
              hlp_scp(ls,istrai)=hlp_scl(ls)
            end do !}

            do j=1,nnlimi+nnstsi !{
              do ls=1,nnatmi !{
                if(j.ne.l_ddt) then !{
                  wldna(j,ls,istrai)=wldna(j,ls,istrai)/hlp_scl(ls)
                  wldra(j,ls,istrai)=wldra(j,ls,istrai)/hlp_scl(ls)
                else !}{
                  wldna(j,ls,istrai)=wldna(j,ls,istrai)*hlp_scl(ls)
                  wldra(j,ls,istrai)=wldra(j,ls,istrai)*hlp_scl(ls)
                end if !}
              end do !}
              do i=1,nnmoli !{
                ls=lmolscl(i)
                if(ls.gt.0 .and. ls.le.nnatmi) then !{
                  if(j.ne.l_ddt) then !{
                    wldnm(j,i,istrai)=wldnm(j,i,istrai)/hlp_scl(ls)
                    wldrm(j,i,istrai)=wldrm(j,i,istrai)/hlp_scl(ls)
                  else !}{
                    wldnm(j,i,istrai)=wldnm(j,i,istrai)*hlp_scl(ls)
                    wldrm(j,i,istrai)=wldrm(j,i,istrai)*hlp_scl(ls)
                  end if !}
                end if !}
              end do !}
            end do !}

            do is=0,ns-1 !{
              ls=b2eatcr(is)
              if(.not.is_neutral(is)) then !{
                if(lstrascl(istrai,ls).gt.0) then !{
                  do ix=1,nx !{
                    do iy=1,ny !{
                      sni(ix,iy,fluids_list(is),istrai)=
     *                 sni(ix,iy,fluids_list(is),istrai)*hlp_scl(ls)
                      smo(ix,iy,fluids_list(is),istrai)=
     *                 smo(ix,iy,fluids_list(is),istrai)*hlp_scl(ls)
                      if (balance_netcdf.ne.0) then !djm Jan2017
                        sni_papl(ix,iy,fluids_list(is),istrai)=
     *                   sni_papl(ix,iy,fluids_list(is),istrai)*
     *                   hlp_scl(ls)
                        sni_pmpl(ix,iy,fluids_list(is),istrai)=
     *                   sni_pmpl(ix,iy,fluids_list(is),istrai)*
     *                   hlp_scl(ls)
                        sni_pipl(ix,iy,fluids_list(is),istrai)=
     *                   sni_pipl(ix,iy,fluids_list(is),istrai)*
     *                   hlp_scl(ls)
                        sni_pppl(ix,iy,fluids_list(is),istrai)=
     *                   sni_pppl(ix,iy,fluids_list(is),istrai)*
     *                   hlp_scl(ls)
                        smo_mapl(ix,iy,fluids_list(is),istrai)=
     *                   smo_mapl(ix,iy,fluids_list(is),istrai)*
     *                   hlp_scl(ls)
                        smo_mmpl(ix,iy,fluids_list(is),istrai)=
     *                   smo_mmpl(ix,iy,fluids_list(is),istrai)*
     *                   hlp_scl(ls)
                        smo_mipl(ix,iy,fluids_list(is),istrai)=
     *                   smo_mipl(ix,iy,fluids_list(is),istrai)*
     *                   hlp_scl(ls)
                        smo_cppv(ix,iy,fluids_list(is),istrai)=
     *                   smo_cppv(ix,iy,fluids_list(is),istrai)*
     *                   hlp_scl(ls)
                      endif
                    end do !}
                  end do !}
                end if !}
              else !}{
                do ix=1,nx !{
                  do iy=1,ny !{
                    sna_paat(ix,iy,ls,istrai)=
     *                 sna_paat(ix,iy,ls,istrai)*hlp_scl(ls)
                    sna_pmat(ix,iy,ls,istrai)=
     *                 sna_pmat(ix,iy,ls,istrai)*hlp_scl(ls)
                    sna_piat(ix,iy,ls,istrai)=
     *                 sna_piat(ix,iy,ls,istrai)*hlp_scl(ls)
                  end do !}
                end do !}
              end if !}
            end do !}
            do is = 1, nnmoli
              ls = lmolscl(is)
              do ix=1,nx
                do iy=1,ny
                  sna_paml(ix,iy,is,istrai)=
     *               sna_paml(ix,iy,is,istrai)*hlp_scl(ls)
                  sna_pmml(ix,iy,is,istrai)=
     *               sna_pmml(ix,iy,is,istrai)*hlp_scl(ls)
                  sna_piml(ix,iy,is,istrai)=
     *               sna_piml(ix,iy,is,istrai)*hlp_scl(ls)
                end do
              end do
            end do
            do is = 1, nnioni
              ls = lionscl(is)
              do ix=1,nx
                do iy=1,ny
                  sna_paio(ix,iy,is,istrai)=
     *               sna_paio(ix,iy,is,istrai)*hlp_scl(ls)
                  sna_pmio(ix,iy,is,istrai)=
     *               sna_pmio(ix,iy,is,istrai)*hlp_scl(ls)
                  sna_piio(ix,iy,is,istrai)=
     *               sna_piio(ix,iy,is,istrai)*hlp_scl(ls)
                end do
              end do
            end do
          else if(crcstra(istrai).eq.'T') then !}{
            do ls = 1,nnatmi !{
              hlp_scp(ls,istrai)=1.0_R8
            end do !}
          end if !}
c!###{
          if(dbg_eir_mc.gt.0) then !{
            do ls=1,nnatmi !{
              print '(i6,2i8,3x,1p,7e13.5)',istrai,ls,
     ,         lstrascl(istrai,ls),tfl(istrai,ls),
     ,         hlp_pi(ls),hlp_po(ls),
     ,         hlp_ni(ls),hlp_no(ls),hlp_s(ls),hlp_scl(ls)
            end do !}
          end if !}
c!###}
        end do !}

        if(dbg_eir_mc.gt.0 .and.
     &       main_call .and. mod(ncall,eirene_mod).eq.0) then !{
c*** Print the particle fluxes for each stratum and species
          print '(a)',
     ,     '%%%=============== total fluxes by stratum and '//
     ,     'species after correction ================%%%'
          print '(a25,10(4x,a5,4x))',' Stratum Species lstrascl',
     ,       ' tfl ','wldpi','wldpo','wldni','wldno',
     ,       'Srcpl','puff ',' d/dt','crsrc','error'
          do istrai=1,nstrat !{
            do ls=1,nnatmi !{
              hlp_pi(ls)=0.0_R8
              hlp_po(ls)=0.0_R8
              hlp_ni(ls)=0.0_R8
              hlp_no(ls)=0.0_R8
              hlp_g(ls) =0.0_R8
              hlp_s(ls) =0.0_R8
              hlp_t(ls) =0.0_R8
            end do !}
            do j=1,nnlimi+nnstsi !{
              if(j.gt.nnlimi .and. j.ne.l_ddt) then !{
                do is=0,ns-1 !{
                  if(.not.is_neutral(is)) then !{
                    hlp_pi(b2eatcr(is))=hlp_pi(b2eatcr(is))+
     +               wldpp(j,fluids_list(is),istrai)
                  end if !}
                end do !}
                do i=1,nnatmi !{
                  hlp_po(i)=hlp_po(i)+wldpa(j,i,istrai)
                  do is=1,nnmoli !{
                    s=mlcmp(i,is)
                    if(s.gt.0.0_R8)
     -               hlp_po(i)=hlp_po(i)+wldpm(j,is,istrai)*s
                  end do !}
                end do !}
              end if !}

              hlp_log=j.le.nnlimi
#ifdef B25_EIRENE
     &                .and. iliin(j) .gt. 0                ! IYS 03.12.2019 Consider non-transparent surfaces only,
                                                           ! transparent ones should not contribute into balance
#endif
              if(hlp_log .or. (j.gt.nnlimi .and. j.ne.l_ddt)) then !{
                if(.not.hlp_log) then !{
                  do i=1,ncns !{
                    hlp_log=hlp_log .or. j.eq.lcns(i)+nnlimi
                  end do !}
                  do i=1,ntns !{
                    hlp_log=hlp_log .or. j.eq.ltns(i)+nnlimi
                  end do !}
                end if !}

                if(hlp_log) then !{
                  do i=1,nnatmi !{
                    hlp_ni(i)=hlp_ni(i)+wldna(j,i,istrai)
                    hlp_no(i)=hlp_no(i)+wldra(j,i,istrai)
                    do is=1,nnmoli !{
                      s=mlcmp(i,is)
                      hlp_ni(i)=hlp_ni(i)+wldnm(j,is,istrai)*s
                      hlp_no(i)=hlp_no(i)+wldrm(j,is,istrai)*s
                    end do !}
                  end do !}
                end if !}
              else if(j.eq.l_ddt) then !}{
                do i=1,nnatmi !{
                  hlp_t(i)=hlp_t(i)+wldna(j,i,istrai)
                  do is=1,nnmoli !{
                    s=mlcmp(i,is)
                    hlp_t(i)=hlp_t(i)+wldnm(j,is,istrai)*s
                  end do !}
                end do !}
              end if !}
            end do !}

            do ix=1,nx !{
              do iy=1,ny !{
                if(resignore(ix-1,iy-1,1).ne.0) then !{
                  do is=0,ns-1 !{
                    if(.not.is_neutral(is)) then !{
                      hlp_s(b2eatcr(is))=hlp_s(b2eatcr(is))+
     +                 sni(ix,iy,fluids_list(is),istrai)
                    end if !}
                  end do !}
                end if !}
              end do !}
            end do !}

            do ls=1,nnatmi !{
              kks=ispc_ref(ls)
              hlp_s(ls)=hlp_s(ls)*tflux(istrai)*flux_scale(istrai)
              t0=hlp_scp(ls,istrai)
              if(use_recyceir.gt.0) then
                recyctmp=recyceir(istrai)
              else
                recyctmp=recyc(kks,istrai)
              endif
              if(crcstra(istrai).eq.'C')
     -           hlp_g(ls)=tflux(istrai)*
     *           gpfc(ls,istrai)*recyctmp
              s=hlp_po(ls)+hlp_g(ls)+hlp_no(ls)-
     -         (hlp_ni(ls)+hlp_s(ls)+hlp_t(ls))
              print '(i6,2i8,3x,1p,10e13.5)',istrai,ls,
     ,         lstrascl(istrai,ls),tfl(istrai,ls),
     ,         hlp_pi(ls),hlp_po(ls),hlp_ni(ls),hlp_no(ls),
     ,         hlp_s(ls),hlp_g(ls),hlp_t(ls),t0,s
              end do !}
          end do !}
          print '(a)',
     ,     '%%%============================================'//
     ,     '=========================================%%%'
        end if !}

c!###{
c        write(hlp_frm1,'(a,i2,a)') '(3x,',nstrat,'i12)'
c        write(hlp_frm2,'(a,i2,a)') '(1x,a5,2x,1p,',nstrat,'e12.4)'
c        do j=1,lns_hlp !{
c          l=lis_hlp(j)
c          print *,'%%%=============== surface ',l
c          write(*,hlp_frm1) (i,i=1,nstrat)
c          write(*,hlp_frm2) 'wldnm',
c     ,                                    (wldnm(l,1,i),i=1,nstrat)
c          write(*,hlp_frm2) 'wldrm',
c     ,                                    (wldrm(l,1,i),i=1,nstrat)
c          write(*,hlp_frm2) 'wldna',
c     ,                                    (wldna(l,1,i),i=1,nstrat)
c          write(*,hlp_frm2) 'wldra',
c     ,                                    (wldra(l,1,i),i=1,nstrat)
c          print *
c        end do !}
c!###}
      endif !ank_mods


csw 26feb2014 accumulate source terms for SOLPS4.3-like averaging procedure
c             (this is new and similar to SOLPS4 style from Kotov,
c              it is NOT IYSs style)
      if(ncall .eq. 0 .and. .not.allocated(sni_av)) then
        allocate(sni_av(0:(size(sni,1)-1),0:(size(sni,2)-1),
     .                     size(sni,3),size(sni,4)))
        allocate(smo_av(0:(size(smo,1)-1),0:(size(smo,2)-1),
     .                     size(smo,3),size(smo,4)))
        allocate(sei_av(0:(size(sei,1)-1),0:(size(sei,2)-1),
     .                     size(sei,3)))
        allocate(see_av(0:(size(see,1)-1),0:(size(see,2)-1),
     .                     size(see,3)))
        sni_av=0.0_R8
        smo_av=0.0_R8
        see_av=0.0_R8
        sei_av=0.0_R8
        if (balance_netcdf.ne.0) then !djm Jan2017
          allocate(sne_av_pael(0:(size(see,1)-1),0:(size(see,2)-1),
     .                    size(see,3)))
          allocate(sne_av_pmel(0:(size(see,1)-1),0:(size(see,2)-1),
     .                    size(see,3)))
          allocate(sna_av_paat(0:(size(sni,1)-1),0:(size(sni,2)-1),
     .                    nnatmi, size(sni,4)))
          allocate(sna_av_pmat(0:(size(sni,1)-1),0:(size(sni,2)-1),
     .                    nnatmi, size(sni,4)))
          allocate(sna_av_piat(0:(size(sni,1)-1),0:(size(sni,2)-1),
     .                    nnatmi, size(sni,4)))
          allocate(sna_av_paml(0:(size(sni,1)-1),0:(size(sni,2)-1),
     .                    nnmoli, size(sni,4)))
          allocate(sna_av_pmml(0:(size(sni,1)-1),0:(size(sni,2)-1),
     .                    nnmoli, size(sni,4)))
          allocate(sna_av_piml(0:(size(sni,1)-1),0:(size(sni,2)-1),
     .                    nnmoli, size(sni,4)))
          allocate(sna_av_paio(0:(size(sni,1)-1),0:(size(sni,2)-1),
     .                    nnioni, size(sni,4)))
          allocate(sna_av_pmio(0:(size(sni,1)-1),0:(size(sni,2)-1),
     .                    nnioni, size(sni,4)))
          allocate(sna_av_piio(0:(size(sni,1)-1),0:(size(sni,2)-1),
     .                    nnioni, size(sni,4)))
          allocate(sni_av_papl(0:(size(sni,1)-1),0:(size(sni,2)-1),
     .                    size(sni,3), size(sni,4)))
          allocate(sni_av_pmpl(0:(size(sni,1)-1),0:(size(sni,2)-1),
     .                    size(sni,3), size(sni,4)))
          allocate(sni_av_pipl(0:(size(sni,1)-1),0:(size(sni,2)-1),
     .                    size(sni,3), size(sni,4)))
          allocate(sni_av_pppl(0:(size(sni,1)-1),0:(size(sni,2)-1),
     .                    size(sni,3), size(sni,4)))
          allocate(smo_av_mapl(0:(size(smo,1)-1),0:(size(smo,2)-1),
     .                    size(smo,3), size(smo,4)))
          allocate(smo_av_mmpl(0:(size(smo,1)-1),0:(size(smo,2)-1),
     .                    size(smo,3), size(smo,4)))
          allocate(smo_av_mipl(0:(size(smo,1)-1),0:(size(smo,2)-1),
     .                    size(smo,3), size(smo,4)))
          allocate(smo_av_cppv(0:(size(smo,1)-1),0:(size(smo,2)-1),
     .                    size(smo,3), size(smo,4)))
          allocate(sei_av_eapl(0:(size(sei,1)-1),0:(size(sei,2)-1),
     .                    size(sei,3)))
          allocate(sei_av_empl(0:(size(sei,1)-1),0:(size(sei,2)-1),
     .                    size(sei,3)))
          allocate(sei_av_eipl(0:(size(sei,1)-1),0:(size(sei,2)-1),
     .                    size(sei,3)))
          allocate(sei_av_eppl(0:(size(sei,1)-1),0:(size(sei,2)-1),
     .                    size(sei,3)))
          allocate(see_av_eael(0:(size(see,1)-1),0:(size(see,2)-1),
     .                    size(see,3)))
          allocate(see_av_emel(0:(size(see,1)-1),0:(size(see,2)-1),
     .                    size(see,3)))
          allocate(see_av_eiel(0:(size(see,1)-1),0:(size(see,2)-1),
     .                    size(see,3)))
          allocate(see_av_epel(0:(size(see,1)-1),0:(size(see,2)-1),
     .                    size(see,3)))
          sne_av_pael=0.0_R8
          sne_av_pmel=0.0_R8
          sna_av_paat=0.0_R8
          sna_av_pmat=0.0_R8
          sna_av_piat=0.0_R8
          sna_av_paml=0.0_R8
          sna_av_pmml=0.0_R8
          sna_av_piml=0.0_R8
          sni_av_papl=0.0_R8
          sni_av_pmpl=0.0_R8
          sni_av_pipl=0.0_R8
          sni_av_pppl=0.0_R8
          smo_av_mapl=0.0_R8
          smo_av_mmpl=0.0_R8
          smo_av_mipl=0.0_R8
          smo_av_cppv=0.0_R8
          see_av_eael=0.0_R8
          see_av_emel=0.0_R8
          see_av_eiel=0.0_R8
          see_av_epel=0.0_R8
          sei_av_eapl=0.0_R8
          sei_av_empl=0.0_R8
          sei_av_eipl=0.0_R8
          sei_av_eppl=0.0_R8
        endif
      endif

      if(neutr_avg .gt. 0) then

       if(main_call.and.mod(ncall,eirene_mod).eq.0) then
        ihist=ihist+1
        if(ihist.eq.neutr_avg+1) then
         sni_av=0.0_R8
         smo_av=0.0_R8
         see_av=0.0_R8
         sei_av=0.0_R8
         if (balance_netcdf.ne.0) then !djm Jan2017
           sne_av_pael=0.0_R8
           sne_av_pmel=0.0_R8
           sna_av_paat=0.0_R8
           sna_av_pmat=0.0_R8
           sna_av_piat=0.0_R8
           sna_av_paml=0.0_R8
           sna_av_pmml=0.0_R8
           sna_av_piml=0.0_R8
           sni_av_papl=0.0_R8
           sni_av_pmpl=0.0_R8
           sni_av_pipl=0.0_R8
           sni_av_pppl=0.0_R8
           smo_av_mapl=0.0_R8
           smo_av_mmpl=0.0_R8
           smo_av_mipl=0.0_R8
           smo_av_cppv=0.0_R8
           see_av_eael=0.0_R8
           see_av_emel=0.0_R8
           see_av_eiel=0.0_R8
           see_av_epel=0.0_R8
           sei_av_eapl=0.0_R8
           sei_av_empl=0.0_R8
           sei_av_eipl=0.0_R8
           sei_av_eppl=0.0_R8
         endif
         dab2_s=0.0_R8
         pfluxa_s=0.0_R8
         pfluxm_s=0.0_R8
         rfluxa_s=0.0_R8
         rfluxm_s=0.0_R8
         ihist = 1
        endif
        sni_av = sni_av + sni
        smo_av = smo_av + smo
        sei_av = sei_av + sei
        see_av = see_av + see
        if (balance_netcdf.ne.0) then !djm Jan2017
          sne_av_pael = sne_av_pael + sne_pael
          sne_av_pmel = sne_av_pmel + sne_pmel
          sna_av_paat = sna_av_paat + sna_paat
          sna_av_pmat = sna_av_pmat + sna_pmat
          sna_av_piat = sna_av_piat + sna_piat
          sna_av_paml = sna_av_paml + sna_paml
          sna_av_pmml = sna_av_pmml + sna_pmml
          sna_av_piml = sna_av_piml + sna_piml
          sna_av_paio = sna_av_paio + sna_paio
          sna_av_pmio = sna_av_pmio + sna_pmio
          sna_av_piio = sna_av_piio + sna_piio
          sni_av_papl = sni_av_papl + sni_papl
          sni_av_pmpl = sni_av_pmpl + sni_pmpl
          sni_av_pipl = sni_av_pipl + sni_pipl
          sni_av_pppl = sni_av_pppl + sni_pppl
          smo_av_mapl = smo_av_mapl + smo_mapl
          smo_av_mmpl = smo_av_mmpl + smo_mmpl
          smo_av_mipl = smo_av_mipl + smo_mipl
          smo_av_cppv = smo_av_cppv + smo_cppv
          sei_av_eapl = sei_av_eapl + sei_eapl
          sei_av_empl = sei_av_empl + sei_empl
          sei_av_eipl = sei_av_eipl + sei_eipl
          sei_av_eppl = sei_av_eppl + sei_eppl
          see_av_eael = see_av_eael + see_eael
          see_av_emel = see_av_emel + see_emel
          see_av_eiel = see_av_eiel + see_eiel
          see_av_epel = see_av_epel + see_epel
        endif
        dab2_s = dab2_s + dab2
        pfluxa_s = pfluxa_s + pfluxa
        pfluxm_s = pfluxm_s + pfluxm
        rfluxa_s = rfluxa_s + rfluxa
        rfluxm_s = rfluxm_s + rfluxm
        dab2 = dab2_s/dble(ihist)
        pfluxa = pfluxa_s/dble(ihist)
        pfluxm = pfluxm_s/dble(ihist)
        rfluxa = rfluxa_s/dble(ihist)
        rfluxm = rfluxm_s/dble(ihist)
       endif
       ih_end=ihist
      else if (nhist.gt.1) then
csw 25feb2014  fix IYS version of source term averaging (running average)
c fixed SMOOTH SOURCES FROM EIRENE, averaging procedure
       if(main_call.and.mod(ncall,eirene_mod).eq.0) then
        ihist=ihist+1
        if (ihist.eq.nhist) full_history = .true.
        if (ihist.eq.nhist+1) ihist = 1
        write(*,*) 'eirene_mc: EIRENE source averaging ', nhist, ihist
     &      ,full_history
        sni_h(:,:,:,:,ihist) = sni
        smo_h(:,:,:,:,ihist) = smo
        see_h(:,:,:,ihist) = see
        sei_h(:,:,:,ihist) = sei
        if (balance_netcdf.ne.0) then !djm Jan2017
          sna_h_paat(:,:,:,:,ihist) = sna_paat
          sna_h_pmat(:,:,:,:,ihist) = sna_pmat
          sna_h_piat(:,:,:,:,ihist) = sna_piat
          sna_h_paml(:,:,:,:,ihist) = sna_paml
          sna_h_pmml(:,:,:,:,ihist) = sna_pmml
          sna_h_piml(:,:,:,:,ihist) = sna_piml
          sna_h_paio(:,:,:,:,ihist) = sna_paio
          sna_h_pmio(:,:,:,:,ihist) = sna_pmio
          sna_h_piio(:,:,:,:,ihist) = sna_piio
          sni_h_papl(:,:,:,:,ihist) = sni_papl
          sni_h_pmpl(:,:,:,:,ihist) = sni_pmpl
          sni_h_pipl(:,:,:,:,ihist) = sni_pipl
          sni_h_pppl(:,:,:,:,ihist) = sni_pppl
          smo_h_mapl(:,:,:,:,ihist) = smo_mapl
          smo_h_mmpl(:,:,:,:,ihist) = smo_mmpl
          smo_h_mipl(:,:,:,:,ihist) = smo_mipl
          smo_h_cppv(:,:,:,:,ihist) = smo_cppv
          sne_h_pael(:,:,:,ihist) = sne_pael
          sne_h_pmel(:,:,:,ihist) = sne_pmel
          see_h_eael(:,:,:,ihist) = see_eael
          see_h_emel(:,:,:,ihist) = see_emel
          see_h_eiel(:,:,:,ihist) = see_eiel
          see_h_epel(:,:,:,ihist) = see_epel
          sei_h_eapl(:,:,:,ihist) = sei_eapl
          sei_h_empl(:,:,:,ihist) = sei_empl
          sei_h_eipl(:,:,:,ihist) = sei_eipl
          sei_h_eppl(:,:,:,ihist) = sei_eppl
        endif
        do iatm = 1,nnatmi
         pfluxa_h(:,:,iatm,1,ihist) = pfluxa(:,:,iatm,1)
         rfluxa_h(:,:,iatm,1,ihist) = rfluxa(:,:,iatm,1)
         dab2_h(:,:,iatm,1,ihist) = dab2(:,:,iatm,1)
        enddo
        do imol = 1,nnmoli
         pfluxm_h(:,:,imol,1,ihist) = pfluxm(:,:,imol,1)
         rfluxm_h(:,:,imol,1,ihist) = rfluxm(:,:,imol,1)
        enddo
        if (full_history) then
         ih_end = nhist
        else
         ih_end = ihist
        end if
        sni_av = 0.0_R8
        smo_av = 0.0_R8
        see_av = 0.0_R8
        sei_av = 0.0_R8
        if (balance_netcdf.ne.0) then !djm Jan2017
          sne_av_pael = 0.0_R8
          sne_av_pmel = 0.0_R8
          sna_av_paat = 0.0_R8
          sna_av_pmat = 0.0_R8
          sna_av_piat = 0.0_R8
          sna_av_paml = 0.0_R8
          sna_av_pmml = 0.0_R8
          sna_av_piml = 0.0_R8
          sni_av_papl = 0.0_R8
          sni_av_pmpl = 0.0_R8
          sni_av_pipl = 0.0_R8
          sni_av_pppl = 0.0_R8
          smo_av_mapl = 0.0_R8
          smo_av_mmpl = 0.0_R8
          smo_av_mipl = 0.0_R8
          smo_av_cppv = 0.0_R8
          see_av_eael = 0.0_R8
          see_av_emel = 0.0_R8
          see_av_eiel = 0.0_R8
          see_av_epel = 0.0_R8
          sei_av_eapl = 0.0_R8
          sei_av_empl = 0.0_R8
          sei_av_eipl = 0.0_R8
          sei_av_eppl = 0.0_R8
        endif
        dab2_s = 0.0_R8
        pfluxa_s = 0.0_R8
        pfluxm_s = 0.0_R8
        rfluxa_s = 0.0_R8
        rfluxm_s = 0.0_R8
        do ih = 1, ih_end
         sni_av = sni_av + sni_h(:,:,:,:,ih)
         smo_av = smo_av + smo_h(:,:,:,:,ih)
         see_av = see_av + see_h(:,:,:,ih)
         sei_av = sei_av + sei_h(:,:,:,ih)
         if (balance_netcdf.ne.0) then !djm Jan2017
           sne_av_pael = sne_av_pael + sne_h_pael(:,:,:,ih)
           sne_av_pmel = sne_av_pmel + sne_h_pmel(:,:,:,ih)
           sna_av_paat = sna_av_paat + sna_h_paat(:,:,:,:,ih)
           sna_av_pmat = sna_av_pmat + sna_h_pmat(:,:,:,:,ih)
           sna_av_piat = sna_av_piat + sna_h_piat(:,:,:,:,ih)
           sna_av_paml = sna_av_paml + sna_h_paml(:,:,:,:,ih)
           sna_av_pmml = sna_av_pmml + sna_h_pmml(:,:,:,:,ih)
           sna_av_piml = sna_av_piml + sna_h_piml(:,:,:,:,ih)
           sni_av_papl = sni_av_papl + sni_h_papl(:,:,:,:,ih)
           sni_av_pmpl = sni_av_pmpl + sni_h_pmpl(:,:,:,:,ih)
           sni_av_pipl = sni_av_pipl + sni_h_pipl(:,:,:,:,ih)
           sni_av_pppl = sni_av_pppl + sni_h_pppl(:,:,:,:,ih)
           smo_av_mapl = smo_av_mapl + smo_h_mapl(:,:,:,:,ih)
           smo_av_mmpl = smo_av_mmpl + smo_h_mmpl(:,:,:,:,ih)
           smo_av_mipl = smo_av_mipl + smo_h_mipl(:,:,:,:,ih)
           smo_av_cppv = smo_av_cppv + smo_h_cppv(:,:,:,:,ih)
           see_av_eael = see_av_eael + see_h_eael(:,:,:,ih)
           see_av_emel = see_av_emel + see_h_emel(:,:,:,ih)
           see_av_eiel = see_av_eiel + see_h_eiel(:,:,:,ih)
           see_av_epel = see_av_epel + see_h_epel(:,:,:,ih)
           sei_av_eapl = sei_av_eapl + sei_h_eapl(:,:,:,ih)
           sei_av_empl = sei_av_empl + sei_h_empl(:,:,:,ih)
           sei_av_eipl = sei_av_eipl + sei_h_eipl(:,:,:,ih)
           sei_av_eppl = sei_av_eppl + sei_h_eppl(:,:,:,ih)
         endif
         do iatm = 1,nnatmi
          dab2_s(:,:,iatm,1) = dab2_s(:,:,iatm,1)
     .                       + dab2_h(:,:,iatm,1,ih)
          pfluxa_s(:,:,iatm,1) = pfluxa_s(:,:,iatm,1)
     .                         + pfluxa_h(:,:,iatm,1,ih)
          rfluxa_s(:,:,iatm,1) = rfluxa_s(:,:,iatm,1)
     .                         + rfluxa_h(:,:,iatm,1,ih)
         enddo
         do imol = 1,nnmoli
          pfluxm_s(:,:,imol,1) = pfluxm_s(:,:,imol,1)
     .                         + pfluxm_h(:,:,imol,1,ih)
          rfluxm_s(:,:,imol,1) = rfluxm_s(:,:,imol,1)
     .                         + rfluxm_h(:,:,imol,1,ih)
         enddo
        enddo
        dab2 = dab2_s / dble(ih_end)
        pfluxa = pfluxa_s / dble(ih_end)
        pfluxm = pfluxm_s / dble(ih_end)
        rfluxa = rfluxa_s / dble(ih_end)
        rfluxm = rfluxm_s / dble(ih_end)
       end if
      else if (nrelax.gt.1) then
cxpb Under-relaxation scheme for smoothing of Eirene sources
       if(main_call.and.mod(ncall,eirene_mod).eq.0) then
        ihist=ihist+1
        if (ihist.eq.nrelax) full_history = .true.
        if (ihist.eq.nrelax+1) ihist = 1
        write(*,*)
     .   'eirene_mc: EIRENE source under-relaxation ', nrelax, ihist
        if (full_history) then
         ih_end = nrelax
        else
         ih_end = ihist
        end if
        relax = rratio(ih_end-1,ih_end)
        sni_s = sni_s*relax + sni/dble(ih_end)
        smo_s = smo_s*relax + smo/dble(ih_end)
        see_s = see_s*relax + see/dble(ih_end)
        sei_s = sei_s*relax + sei/dble(ih_end)
        if (balance_netcdf.ne.0) then !djm Jan2017
          sne_s_pael = sne_s_pael*relax + sne_pael/dble(ih_end)
          sne_s_pmel = sne_s_pmel*relax + sne_pmel/dble(ih_end)
          sna_s_paat = sna_s_paat*relax + sna_paat/dble(ih_end)
          sna_s_pmat = sna_s_pmat*relax + sna_pmat/dble(ih_end)
          sna_s_piat = sna_s_piat*relax + sna_piat/dble(ih_end)
          sna_s_paml = sna_s_paml*relax + sna_paml/dble(ih_end)
          sna_s_pmml = sna_s_pmml*relax + sna_pmml/dble(ih_end)
          sna_s_piml = sna_s_piml*relax + sna_piml/dble(ih_end)
          sna_s_paio = sna_s_paio*relax + sna_paio/dble(ih_end)
          sna_s_pmio = sna_s_pmio*relax + sna_pmio/dble(ih_end)
          sna_s_piio = sna_s_piio*relax + sna_piio/dble(ih_end)
          sni_s_papl = sni_s_papl*relax + sni_papl/dble(ih_end)
          sni_s_pmpl = sni_s_pmpl*relax + sni_pmpl/dble(ih_end)
          sni_s_pipl = sni_s_pipl*relax + sni_pipl/dble(ih_end)
          sni_s_pppl = sni_s_pppl*relax + sni_pppl/dble(ih_end)
          smo_s_mapl = smo_s_mapl*relax + smo_mapl/dble(ih_end)
          smo_s_mmpl = smo_s_mmpl*relax + smo_mmpl/dble(ih_end)
          smo_s_mipl = smo_s_mipl*relax + smo_mipl/dble(ih_end)
          smo_s_cppv = smo_s_cppv*relax + smo_cppv/dble(ih_end)
          see_s_eael = see_s_eael*relax + see_eael/dble(ih_end)
          see_s_emel = see_s_emel*relax + see_emel/dble(ih_end)
          see_s_eiel = see_s_eiel*relax + see_eiel/dble(ih_end)
          see_s_epel = see_s_epel*relax + see_epel/dble(ih_end)
          sei_s_eapl = sei_s_eapl*relax + sei_eapl/dble(ih_end)
          sei_s_empl = sei_s_empl*relax + sei_empl/dble(ih_end)
          sei_s_eipl = sei_s_eipl*relax + sei_eipl/dble(ih_end)
          sei_s_eppl = sei_s_eppl*relax + sei_eppl/dble(ih_end)
        endif
        do iatm = 1,nnatmi
         dab2_s(:,:,iatm,1) = dab2_s(:,:,iatm,1)*relax
     .                      + dab2(:,:,iatm,1)/dble(ih_end)
         pfluxa_s(:,:,iatm,1) = pfluxa_s(:,:,iatm,1)*relax
     .                        + pfluxa(:,:,iatm,1)/dble(ih_end)
         rfluxa_s(:,:,iatm,1) = rfluxa_s(:,:,iatm,1)*relax
     .                        + rfluxa(:,:,iatm,1)/dble(ih_end)
        enddo
        do imol = 1,nnmoli
         pfluxm_s(:,:,imol,1) = pfluxm_s(:,:,imol,1)*relax
     .                        + pfluxm(:,:,imol,1)/dble(ih_end)
         rfluxm_s(:,:,imol,1) = rfluxm_s(:,:,imol,1)*relax
     .                        + rfluxm(:,:,imol,1)/dble(ih_end)
        enddo
        sni_av = sni_s*dble(ih_end)
        smo_av = smo_s*dble(ih_end)
        sei_av = sei_s*dble(ih_end)
        see_av = see_s*dble(ih_end)
        if (balance_netcdf.ne.0) then !djm Jan2017
          sne_av_pael = sne_s_pael*dble(ih_end)
          sne_av_pmel = sne_s_pmel*dble(ih_end)
          sna_av_paat = sna_s_paat*dble(ih_end)
          sna_av_pmat = sna_s_pmat*dble(ih_end)
          sna_av_piat = sna_s_piat*dble(ih_end)
          sna_av_paml = sna_s_paml*dble(ih_end)
          sna_av_pmml = sna_s_pmml*dble(ih_end)
          sna_av_piml = sna_s_piml*dble(ih_end)
          sni_av_papl = sni_s_papl*dble(ih_end)
          sni_av_pmpl = sni_s_pmpl*dble(ih_end)
          sni_av_pipl = sni_s_pipl*dble(ih_end)
          sni_av_pppl = sni_s_pppl*dble(ih_end)
          smo_av_mapl = smo_s_mapl*dble(ih_end)
          smo_av_mmpl = smo_s_mmpl*dble(ih_end)
          smo_av_mipl = smo_s_mipl*dble(ih_end)
          smo_av_cppv = smo_s_cppv*dble(ih_end)
          sei_av_eapl = sei_s_eapl*dble(ih_end)
          sei_av_empl = sei_s_empl*dble(ih_end)
          sei_av_eipl = sei_s_eipl*dble(ih_end)
          sei_av_eppl = sei_s_eppl*dble(ih_end)
          see_av_eael = see_s_eael*dble(ih_end)
          see_av_emel = see_s_emel*dble(ih_end)
          see_av_eiel = see_s_eiel*dble(ih_end)
          see_av_epel = see_s_epel*dble(ih_end)
        endif
        dab2 = dab2_s
        pfluxa = pfluxa_s
        pfluxm = pfluxm_s
        rfluxa = rfluxa_s
        rfluxm = rfluxm_s
       end if
      else
       ih_end=1
       sni_av = sni
       smo_av = smo
       sei_av = sei
       see_av = see
       if (balance_netcdf.ne.0) then !djm Jan2017
         sne_av_pael = sne_pael
         sne_av_pmel = sne_pmel
         sna_av_paat = sna_paat
         sna_av_pmat = sna_pmat
         sna_av_piat = sna_piat
         sna_av_paml = sna_paml
         sna_av_pmml = sna_pmml
         sna_av_piml = sna_piml
         sna_av_paio = sna_paio
         sna_av_pmio = sna_pmio
         sna_av_piio = sna_piio
         sni_av_papl = sni_papl
         sni_av_pmpl = sni_pmpl
         sni_av_pipl = sni_pipl
         sni_av_pppl = sni_pppl
         smo_av_mapl = smo_mapl
         smo_av_mmpl = smo_mmpl
         smo_av_mipl = smo_mipl
         smo_av_cppv = smo_cppv
         sei_av_eapl = sei_eapl
         sei_av_empl = sei_empl
         sei_av_eipl = sei_eipl
         sei_av_eppl = sei_eppl
         see_av_eael = see_eael
         see_av_emel = see_emel
         see_av_eiel = see_eiel
         see_av_epel = see_epel
       endif
      endif

csw 2012 fixed rescaling procedure
c  otherwise rescaling procedure does not work with internal iterations nstg2>1
c  (sclrtio/sclalph will be reinitialised which should not happen at istg2 > 1)
c  this is part of phase 2 of rescaling process
c
      if(main_call.and.mod(ncall,eirene_mod).eq.0.and.ank_mods.ne.0)
     .   then

c*** Calculate the basic quantities - i.e. the neutral source strengths
c*** and total plasma sources for each stratum for subsequent rescaling

        do istrai=1, nstrat !{

c*** Accumulators:
c***  sclsrcn : the neutral sources (source strengths)
c***  sclbase : the plasma sources

          do j=0,nnatmi !{
            sclsrcn(j,istrai)=0.0_R8
            sclbase(j,istrai)=0.0_R8
          end do !}
          do is=0,ns-1 !{
            if(.not.is_neutral(is)) then !{
              j=b2eatcr(is)
              if(use_recyceir.gt.0) then
                recyctmp=recyceir(istrai)
              else
                recyctmp=recyc(is,istrai)
              endif
              do iy=1,ny !{
                do ix=1,nx !{
                  if(resignore(ix-1,iy-1,1).ne.0) then !{
                    s=sni(ix,iy,fluids_list(is),istrai)*tflux(istrai)*
     *               recyctmp*flux_scale(istrai)
                    if(crcstra(istrai).eq.'V' .and. s.lt.0.0_R8) then !{
                      sclsrcn(j,istrai)=sclsrcn(j,istrai)-s
                      sclsrcn(0,istrai)=sclsrcn(0,istrai)-s
                    else !}{
                      sclbase(j,istrai)=sclbase(j,istrai)+s
                      sclbase(0,istrai)=sclbase(0,istrai)+s
                    end if !}
                  end if !}
                end do !}
              end do !}
            end if !}
          end do !}

          do j=1,nnatmi !{
            is=ispc_ref(j)
            if(use_recyceir.gt.0) then
              recyctmp=recyceir(istrai)
            else
              recyctmp=recyc(is,istrai)
            endif
            if(crcstra(istrai).eq.'E' .or.
     .         crcstra(istrai).eq.'W' .or.
     .         crcstra(istrai).eq.'N' .or.
     .         crcstra(istrai).eq.'S' ) then !{
              sclsrcn(j,istrai)=sclsrcn(j,istrai)+
     +         tfl(istrai,j)*recyctmp
              sclsrcn(0,istrai)=sclsrcn(0,istrai)+
     +         tfl(istrai,j)*recyctmp
            else if(crcstra(istrai).eq.'C'
c     .                          .or. crcstra(istrai).eq.'U'
     .       ) then !}{
              sclsrcn(j,istrai)=sclsrcn(j,istrai)+
     +         tflux(istrai)*gpfc(j,istrai)*recyctmp
            end if !}

            l=lstrascl(istrai,j)
            if(l.gt.nnatmi .or. l.le.0) l=j
            if(sclsrcn(l,istrai).gt.0.0_R8) then !{
              sclbase(j,istrai)=sclbase(j,istrai)/sclsrcn(l,istrai)
            else !}{
              if(sclsrcn(l,istrai).lt.0.0_R8) write(*,*)
     -         '@ neutral source < 0. species,stratum =',l,istrai
              sclbase(j,istrai)=0.0_R8
            end if !}

            sclalph(j,istrai)=1.0_R8
            sclstra(j,istrai)=1.0_R8
          end do !}

          sclsrcn(0,istrai)=0.0_R8
          do j=1,nnatmi !{
            hlp_log=.false.
            do k=1,nnatmi !{
              hlp_log = hlp_log .or. lstrascl(istrai,k).eq.j
            end do !}
            if(hlp_log) sclsrcn(0,istrai)=sclsrcn(0,istrai)+
     +       sclsrcn(j,istrai)
          end do !}
          if(sclsrcn(0,istrai).gt.0.0_R8) then !{
            sclbase(0,istrai)=sclbase(0,istrai)/sclsrcn(0,istrai)
          else !}{
            if(sclsrcn(0,istrai).lt.0.0_R8) write(*,*)
     -       '@ neutral source < 0. species,stratum =',0,istrai
            sclbase(0,istrai)=0.0_R8
          end if !}

          sclalph(0,istrai)=1.0_R8
          sclstra(0,istrai)=1.0_R8
        end do !}
      endif ! ank_mods

      core_flux = 0.0_R8
      if (ionising_core.ne.0) then
        do j = 1, nnlimi+nnstsi
          hlp_cor = .false.
          do i=1,ncns !{
            hlp_cor = hlp_cor .or. j.eq.lcns(i)+nnlimi
          end do !}
          if(hlp_cor) then !{
            do ls=1,nnatmi !{
              core_flux(ls)=core_flux(ls)+
     .         (wldna(j,ls,0)-wldra(j,ls,0)-wldpa(j,ls,0))
              do is=1,nnmoli !{
                s=mlcmp(ls,is)
                core_flux(ls)=core_flux(ls)+
     .           (wldnm(j,is,0)-wldrm(j,is,0)-wldpm(j,is,0))*s
                if (dbg_eir_mc.ge.12 .and. main_call) then
                  write(*,'(a,5i4)') 'J, ISTRAI, IATM, IMOL, S',
     .             j, istrai, ls, is, mlcmp(ls,is)
                  write(*,'(a,1p,5(1x,1e14.7))')
     .             'wldna, wldra, wldnm, wldrm, core_flux',
     .              wldna(j,ls,0), wldra(j,ls,0)+wldpa(j,ls,0),
     .              wldnm(j,is,0), wldrm(j,is,0)+wldpm(j,is,0),
     .              core_flux(ls)
                endif
              end do !}
            end do !}
          end if !}
        end do !}
      endif

csw 2012 now transfer (rescaled) sources into sna0, smo0, she0, shi0,
c   and apply linearisation

      ssni=0.0_R8
      ssmo=0.0_R8
      ssee=0.0_R8
      ssei=0.0_R8
      radarr=0.0_R8
      rad_neutr=0.0_R8
      rdneureg=0.0_R8
      recycireg=0.0_R8
      recycereg=0.0_R8
      sseereg=0.0_R8
      sseireg=0.0_R8
      ssnireg=0.0_R8
      ssnastra=0.0_R8
      sseicorr=0.0_R8
      sseistab=0.0_R8
      src_ni(:,3,:)=0.0_R8
      src_mo(:,3,:)=0.0_R8
      src_ne(:,3)=0.0_R8
      src_ms(:,3)=0.0_R8
      src_ee(:,3)=0.0_R8
      src_ei(:,3)=0.0_R8
      src_es(:,3)=0.0_R8
      src_ne(:,3)=0.0_R8
      smo_tot=0.0_R8
      sni_tot=0.0_R8
      do is = 0, ns-1 !{
        if(.not.is_neutral(is)) then !{
          j=b2eatcr(is)
          ifl=fluids_list(is)
          do istrai=1,nstrat !{
            if(use_recyceir.gt.0) then !{
              recyctmp=recyceir(istrai)
            else !}{
              recyctmp=recyc(is,istrai)
            endif !}
            do iy = -1, ny !{
              do ix = -1, nx !{
                smo_tot(ix+1,iy+1,ifl)=smo_tot(ix+1,iy+1,ifl)+
     .           smo_av(ix+1,iy+1,ifl,istrai)*tflux(istrai)*
     *           recyctmp*flux_scale(istrai)*sclstra(j,istrai)/
     *           dble(ih_end)
                t00 =
     .           sni_av(ix+1,iy+1,ifl,istrai)*tflux(istrai)*
     *           recyctmp*flux_scale(istrai)/dble(ih_end)
                if(crcstra(istrai).ne.'V' .or. t00.gt.0.0_R8) then
                  if (iscltst.eq.0) then
                    t00=t00*sclalph(j,istrai)
                  else
                    t00=t00*sclstra(j,istrai)
                  endif
                endif
                sni_tot(ix+1,iy+1,ifl)=sni_tot(ix+1,iy+1,ifl)+t00
              enddo !}
            enddo !}
          enddo !}
        endif !}
      enddo !}
      if (balance_netcdf.ne.0) then
        eirene_mc_pael_sne0to1 = 0.0_R8 !djm Jan2017
        eirene_mc_papl_sna0to1 = 0.0_R8
        eirene_mc_pmel_sne0to1 = 0.0_R8
        eirene_mc_pmpl_sna0to1 = 0.0_R8
        eirene_mc_pipl_sna0to1 = 0.0_R8
        eirene_mc_pppl_sna0to1 = 0.0_R8
        eirene_mc_mapl_smo0to3 = 0.0_R8
        eirene_mc_mmpl_smo0to3 = 0.0_R8
        eirene_mc_mipl_smo0to3 = 0.0_R8
        eirene_mc_cppv_smo0to3 = 0.0_R8
        eirene_mc_eael_she0to3 = 0.0_R8
        eirene_mc_emel_she0to3 = 0.0_R8
        eirene_mc_eiel_she0to3 = 0.0_R8
        eirene_mc_epel_she0to3 = 0.0_R8
        eirene_mc_eapl_shi0to3 = 0.0_R8
        eirene_mc_empl_shi0to3 = 0.0_R8
        eirene_mc_eipl_shi0to3 = 0.0_R8
        eirene_mc_eppl_shi0to3 = 0.0_R8
      endif
      do istrai=1,nstrat !{
* by ns
        recyc_hlp=0.0_R8
        do is=0,ns-1 !{
          j=b2eatcr(is)
          if(use_recyceir.gt.0) then
            recyctmp=recyceir(istrai)
          else
            recyctmp=recyc(is,istrai)
          endif
          recyc_hlp=max(recyc_hlp,recyctmp)
          if(.not.is_neutral(is)) then !{
            ifl=fluids_list(is)
            do iy=-1,ny !{
              do ix=-1,nx !{
                if(leftix(ix,iy).ne.-2 .and.
     &             rightix(ix,iy).ne.nx+1 .and.
     &             bottomiy(ix,iy).ne.-2 .and.
     &             topiy(ix,iy).ne.ny+1 .and.
     &             region(ix,iy,0).ne.0) then !{
* continuity
                  rf0=RF_NEUT(1)

csw 26feb2014 averaging
c                  t00=sni(ix+1,iy+1,ifl,istrai)*tflux(istrai)*
c     *             recyctmp*flux_scale(istrai)
                  t00=sni_av(ix+1,iy+1,ifl,istrai)*tflux(istrai)*
     *             recyctmp*flux_scale(istrai)/
     *             dble(ih_end)

                  if(dpc_fix.lt.0) then
                    write(*,*) 'You are running the code with ' //
     1               'a default option that has been determined ' //
     2               'to be wrong!'
                    write(*,*) 'If you really want to do this, ' //
     1               'you will need to set "eirene_dpc_fix" ' //
     2               'in "b2mn.dat" to "0"'
                    call xertst(dpc_fix.ge.0,
     1               'Invalid setting of "eirene_dpc_fix"')
                    call xerrab ('Error in value of eirene_dpc_fix!')
                  else if(dpc_fix.eq.0) then
                                                  ! original wrong code
                  else if(dpc_fix.eq.1) then
                    t0=t00                        ! dpc fix!
                    if (balance_netcdf.ne.0) then
                      t0_a = sni_av_papl(ix+1,iy+1,ifl,istrai)*  !djm Jan2017
     *                       tflux(istrai)*
     *                       recyctmp*flux_scale(istrai)/
     *                       dble(ih_end)
                      t0_m = sni_av_pmpl(ix+1,iy+1,ifl,istrai)*
     *                       tflux(istrai)*
     *                       recyctmp*flux_scale(istrai)/
     *                       dble(ih_end)
                      t0_i = sni_av_pipl(ix+1,iy+1,ifl,istrai)*
     *                       tflux(istrai)*
     *                       recyctmp*flux_scale(istrai)/
     *                       dble(ih_end)
                      t0_p = sni_av_pppl(ix+1,iy+1,ifl,istrai)*
     *                       tflux(istrai)*
     *                       recyctmp*flux_scale(istrai)/
     *                       dble(ih_end)
                    endif
                  else if(dpc_fix.eq.2) then
                    t0=0.0_R8                     ! at least get the wrong code to go wrong safely!
                    if (balance_netcdf.ne.0) then
                      t0_a=0.0_R8 !djm Jan2017
                      t0_m=0.0_R8
                      t0_i=0.0_R8
                      t0_p=0.0_R8
                    endif
                  endif
                  if(crcstra(istrai).ne.'V' .or. t00.gt.0.0_R8) then
                    t0=t00*sclalph(j,istrai)
                    if (balance_netcdf.ne.0) then !djm Jan2017
                      t0_a=t0_a*sclalph(j,istrai)
                      t0_m=t0_m*sclalph(j,istrai)
                      t0_i=t0_i*sclalph(j,istrai)
                      t0_p=t0_p*sclalph(j,istrai)
                    endif
                  endif

c                  if (t0.lt.0.0_R8.and.main_call.and.
c     &                crcstra(istrai).ne.'V'.and.output_style.ge.1) then !{
c                    write(6,*) 'neutral particle source < 0! ',t0
c                    write(6,*) 'istrai,is,ix,iy ',istrai,is,ix,iy
c                  else if (t0.gt.0.0_R8.and.main_call.and.
c     &                crcstra(istrai).eq.'V'.and.output_style.ge.1) then !{
c                    write(6,*) 'recombination particle source > 0! ',t0
c                    write(6,*) 'istrai,is,ix,iy ',istrai,is,ix,iy
c                  end if !}

                  if(linearisation.eq.1) then !{
                    if(t0.ge.0.0_R8) then !{
                      sna0(ix,iy,0,is)=sna0(ix,iy,0,is)+t0
                    else !}{
                      sna0(ix,iy,1,is)=sna0(ix,iy,1,is)+t0/na(ix,iy,is)
                    end if !}
                    sna0(ix,iy,0,is)=sna0(ix,iy,0,is)+rf0*abs(t0)
                    sna0(ix,iy,1,is)=sna0(ix,iy,1,is)-
     1               rf0*abs(t0)/na(ix,iy,is)
                    if (balance_netcdf.ne.0) then
cdjm Jan2017 Calculate linearised components of source
                      if (t0.ge.0.0_R8) then
                        eirene_mc_papl_sna0to1(ix,iy,0,is,istrai)=
     &                    eirene_mc_papl_sna0to1(ix,iy,0,is,istrai)+
     &                    t0_a*(1.0_R8+rf0)
                        eirene_mc_pmpl_sna0to1(ix,iy,0,is,istrai)=
     &                    eirene_mc_pmpl_sna0to1(ix,iy,0,is,istrai)+
     &                    t0_m*(1.0_R8+rf0)
                        eirene_mc_pipl_sna0to1(ix,iy,0,is,istrai)=
     &                    eirene_mc_pipl_sna0to1(ix,iy,0,is,istrai)+
     &                    t0_i*(1.0_R8+rf0)
                        eirene_mc_pppl_sna0to1(ix,iy,0,is,istrai)=
     &                    eirene_mc_pppl_sna0to1(ix,iy,0,is,istrai)+
     &                    t0_p*(1.0_R8+rf0)
                        eirene_mc_papl_sna0to1(ix,iy,1,is,istrai)=
     &                    eirene_mc_papl_sna0to1(ix,iy,1,is,istrai)-
     &                    rf0*t0_a/na(ix,iy,is)
                        eirene_mc_pmpl_sna0to1(ix,iy,1,is,istrai)=
     &                    eirene_mc_pmpl_sna0to1(ix,iy,1,is,istrai)-
     &                    rf0*t0_m/na(ix,iy,is)
                        eirene_mc_pipl_sna0to1(ix,iy,1,is,istrai)=
     &                    eirene_mc_pipl_sna0to1(ix,iy,1,is,istrai)-
     &                    rf0*t0_i/na(ix,iy,is)
                        eirene_mc_pppl_sna0to1(ix,iy,1,is,istrai)=
     &                    eirene_mc_pppl_sna0to1(ix,iy,1,is,istrai)-
     &                    rf0*t0_p/na(ix,iy,is)
c$$$                    if ((ix.eq.90).and.(iy.eq.20)) then !DJM DEBUG
c$$$                      write(*,*) 'sna0',t0+rf0*abs(t0),
c$$$     &                           (t0_a+t0_m+t0_i+t0_p)*(1.0_R8+rf0)
c$$$                      write(*,*) 'sna1',-rf0*abs(t0)/na(ix,iy,is),
c$$$     &                          -(t0_a+t0_m+t0_i+t0_p)*rf0/na(ix,iy,is)
c$$$                    endif
                      else
                        eirene_mc_papl_sna0to1(ix,iy,0,is,istrai)=
     &                    eirene_mc_papl_sna0to1(ix,iy,0,is,istrai)-
     &                    rf0*t0_a
                        eirene_mc_pmpl_sna0to1(ix,iy,0,is,istrai)=
     &                    eirene_mc_pmpl_sna0to1(ix,iy,0,is,istrai)-
     &                    rf0*t0_m
                        eirene_mc_pipl_sna0to1(ix,iy,0,is,istrai)=
     &                    eirene_mc_pipl_sna0to1(ix,iy,0,is,istrai)-
     &                    rf0*t0_i
                        eirene_mc_pppl_sna0to1(ix,iy,0,is,istrai)=
     &                    eirene_mc_pppl_sna0to1(ix,iy,0,is,istrai)-
     &                    rf0*t0_p
                        eirene_mc_papl_sna0to1(ix,iy,1,is,istrai)=
     &                    eirene_mc_papl_sna0to1(ix,iy,1,is,istrai)+
     &                    (1.0_R8+rf0)*t0_a/na(ix,iy,is)
                        eirene_mc_pmpl_sna0to1(ix,iy,1,is,istrai)=
     &                    eirene_mc_pmpl_sna0to1(ix,iy,1,is,istrai)+
     &                    (1.0_R8+rf0)*t0_m/na(ix,iy,is)
                        eirene_mc_pipl_sna0to1(ix,iy,1,is,istrai)=
     &                    eirene_mc_pipl_sna0to1(ix,iy,1,is,istrai)+
     &                    (1.0_R8+rf0)*t0_i/na(ix,iy,is)
                        eirene_mc_pppl_sna0to1(ix,iy,1,is,istrai)=
     &                    eirene_mc_pppl_sna0to1(ix,iy,1,is,istrai)+
     &                    (1.0_R8+rf0)*t0_p/na(ix,iy,is)
c$$$                    if ((ix.eq.90).and.(iy.eq.20)) then !DJM DEBUG
c$$$                      write(*,*) 'sna0',rf0*abs(t0),
c$$$     &                           -(t0_a+t0_m+t0_i+t0_p)*rf0
c$$$                      write(*,*) 'sna1',t0/na(ix,iy,is)-
c$$$     &                                  rf0*abs(t0)/na(ix,iy,is),
c$$$     &                   (t0_a+t0_m+t0_i+t0_p)*(1.0_R8+rf0)/na(ix,iy,is)
c$$$                    endif
                      endif
                    endif
                  elseif(linearisation.eq.0) then !}{
                    sna0(ix,iy,0,is)=sna0(ix,iy,0,is)+t0
                    if (balance_netcdf.ne.0) then
                      eirene_mc_papl_sna0to1(ix,iy,0,is,istrai)= !djm Jan2017
     &                 eirene_mc_papl_sna0to1(ix,iy,0,is,istrai)+t0_a
                      eirene_mc_pmpl_sna0to1(ix,iy,0,is,istrai)=
     &                 eirene_mc_pmpl_sna0to1(ix,iy,0,is,istrai)+t0_m
                      eirene_mc_pipl_sna0to1(ix,iy,0,is,istrai)=
     &                 eirene_mc_pipl_sna0to1(ix,iy,0,is,istrai)+t0_i
                      eirene_mc_pppl_sna0to1(ix,iy,0,is,istrai)=
     &                 eirene_mc_pppl_sna0to1(ix,iy,0,is,istrai)+t0_p
                    endif
                  end if !}
                  if(crcstra(istrai).ne.'V' .or. t00.gt.0.0_R8)
     >             t0=t00*sclstra(j,istrai)
                  ssni(istrai,is)=ssni(istrai,is)+t0
c*** calculate the integral sources for nregion regions, 3 instants:
c*** before correction (r), after correction (c), after iteration (i)
                  if (main_call.and.mod(ncall,eirene_mod).eq.0) then
                    src_ni(region(ix,iy,0),1,b2espcr(is))=
     .               src_ni(region(ix,iy,0),1,b2espcr(is))+t0
                    src_ni(region(ix,iy,0),2,b2espcr(is))=
     .               src_ni(region(ix,iy,0),2,b2espcr(is))+t0
                    cor_ni(ix,iy,b2espcr(is))=
     .               cor_ni(ix,iy,b2espcr(is))+t0
                    src_ni(0,1,b2espcr(is))=src_ni(0,1,b2espcr(is))+t0
                    src_ni(0,2,b2espcr(is))=src_ni(0,2,b2espcr(is))+t0
                  else
                    src_ni(region(ix,iy,0),3,b2espcr(is))=
     .               src_ni(region(ix,iy,0),3,b2espcr(is))+t0
                    src_ni(0,3,b2espcr(is))=src_ni(0,3,b2espcr(is))+t0
                  endif

csw 21feb2012, add l_neutrad (SOLPS4.3)
                  if(l_neutrad.lt.0) then
                    radarr(ix,iy)=radarr(ix,iy)+t0*rpt(ix,iy,is)*ev
                    rad_neutr=rad_neutr+t0*rpt(ix,iy,is)*ev
                    rdneureg(region(ix,iy,0))=
     .               rdneureg(region(ix,iy,0))+t0*rpt(ix,iy,is)*ev
                    recycereg(region(ix,iy,0))=
     .               recycereg(region(ix,iy,0))+t0*rpt(ix,iy,is)*ev
                  endif
                  ssnireg(b2espcr(is),region(ix,iy,0))=
     .             ssnireg(b2espcr(is),region(ix,iy,0))+t0
                  ssnastra(j,istrai)=
     .             ssnastra(j,istrai)+t0
* parallel momentum
                  rf0=RF_NEUT(2)
csw 26feb2014 averaging
c                  t0=smo(ix+1,iy+1,ifl,istrai)*tflux(istrai)*
c     *             recyctmp*flux_scale(istrai)*sclstra(j,istrai)
                  t0=smo_av(ix+1,iy+1,ifl,istrai)*tflux(istrai)*
     *             recyctmp*flux_scale(istrai)*sclstra(j,istrai)/
     *             dble(ih_end)*hz(ix,iy)
                  if (balance_netcdf.ne.0) then !djm Jan2017
                    t0_a=smo_av_mapl(ix+1,iy+1,ifl,istrai)*
     *                   tflux(istrai)*
     *                   recyctmp*flux_scale(istrai)*sclstra(j,istrai)/
     *                   dble(ih_end)*hz(ix,iy)
                    t0_m=smo_av_mmpl(ix+1,iy+1,ifl,istrai)*
     *                   tflux(istrai)*
     *                   recyctmp*flux_scale(istrai)*sclstra(j,istrai)/
     *                   dble(ih_end)*hz(ix,iy)
                    t0_i=smo_av_mipl(ix+1,iy+1,ifl,istrai)*
     *                   tflux(istrai)*
     *                   recyctmp*flux_scale(istrai)*sclstra(j,istrai)/
     *                   dble(ih_end)*hz(ix,iy)
                    t0_p=smo_av_cppv(ix+1,iy+1,ifl,istrai)*
     *                   tflux(istrai)*
     *                   recyctmp*flux_scale(istrai)*sclstra(j,istrai)/
     *                   dble(ih_end)*hz(ix,iy)
                  endif

                  VTHE=SQRT(TE(IX,IY)/me)
c SCALE VTHE TO GET A VELOCITY CLOSER TO THE ION VELOCITIES
                  VTHE=VTHE/sqrt(am(is)*mp/me)
c
                  if (linearisation.eq.1) then !{
cxpb 2002.01.28 moved to cell-centered velocities
                    smo0(IX,IY,0,IS)=smo0(IX,IY,0,IS)+
     .               t0+RF0*ABS(t0)*UA(IX,IY,IS)/VTHE
                    smo0(IX,IY,1,IS)=smo0(IX,IY,1,IS)-
     .               RF0*ABS(t0)/VTHE
cxpb was                    smo0(IX,IY,0,IS)=smo0(IX,IY,0,IS)+
cxpb was     .               (t0/2.0_R8)+RF0*ABS(t0/2.0_R8)*UA(IX,IY,IS)/VTHE
cxpb was                    smo0(IX,IY,1,IS)=smo0(IX,IY,1,IS)-
cxpb was     .               RF0*ABS((t0/2.0_R8))/VTHE

cxpb was                    smo0(RIGHTIX(IX,IY),RIGHTIY(IX,IY),0,IS)=
cxpb was     1               smo0(RIGHTIX(IX,IY),RIGHTIY(IX,IY),0,IS)+
cxpb was     .               (t0/2.0_R8)+RF0*ABS(t0/2.0_R8)*
cxpb was     1               UA(RIGHTIX(IX,IY),RIGHTIY(IX,IY),IS)/VTHE
cxpb was                    smo0(RIGHTIX(IX,IY),RIGHTIY(IX,IY),1,IS)=
cxpb was     1               smo0(RIGHTIX(IX,IY),RIGHTIY(IX,IY),1,IS)-
cxpb was     .               RF0*ABS((t0/2.0_R8))/VTHE
                    if (balance_netcdf.ne.0) then
cdjm Jan2017 Calculate linearised components of source
                      if (t0.ge.0.0_R8) then
                        eirene_mc_mapl_smo0to3(ix,iy,0,is,istrai)=
     &                   eirene_mc_mapl_smo0to3(ix,iy,0,is,istrai)+
     &                   (1.0_R8+rf0*ua(ix,iy,is)/vthe)*t0_a
                        eirene_mc_mmpl_smo0to3(ix,iy,0,is,istrai)=
     &                   eirene_mc_mmpl_smo0to3(ix,iy,0,is,istrai)+
     &                   (1.0_R8+rf0*ua(ix,iy,is)/vthe)*t0_m
                        eirene_mc_mipl_smo0to3(ix,iy,0,is,istrai)=
     &                   eirene_mc_mipl_smo0to3(ix,iy,0,is,istrai)+
     &                   (1.0_R8+rf0*ua(ix,iy,is)/vthe)*t0_i
                        eirene_mc_cppv_smo0to3(ix,iy,0,is,istrai)=
     &                   eirene_mc_cppv_smo0to3(ix,iy,0,is,istrai)+
     &                   (1.0_R8+rf0*ua(ix,iy,is)/vthe)*t0_p
                        eirene_mc_mapl_smo0to3(ix,iy,1,is,istrai)=
     &                   eirene_mc_mapl_smo0to3(ix,iy,1,is,istrai)-
     &                   rf0*t0_a/vthe
                        eirene_mc_mmpl_smo0to3(ix,iy,1,is,istrai)=
     &                   eirene_mc_mmpl_smo0to3(ix,iy,1,is,istrai)-
     &                   rf0*t0_m/vthe
                        eirene_mc_mipl_smo0to3(ix,iy,1,is,istrai)=
     &                   eirene_mc_mipl_smo0to3(ix,iy,1,is,istrai)-
     &                   rf0*t0_i/vthe
                        eirene_mc_cppv_smo0to3(ix,iy,1,is,istrai)=
     &                   eirene_mc_cppv_smo0to3(ix,iy,1,is,istrai)-
     &                   rf0*t0_p/vthe
c$$$                    if ((ix.eq.90).and.(iy.eq.20)) then !DJM DEBUG
c$$$                      write(*,*) 'smo0',
c$$$     &                           t0+RF0*ABS(t0)*UA(IX,IY,IS)/VTHE,
c$$$     &              (t0_a+t0_m+t0_i+t0_p)*(1.0_R8+rf0*ua(ix,iy,is)/vthe)
c$$$                      write(*,*) 'smo1',-RF0*ABS(t0)/VTHE,
c$$$     &                          -(t0_a+t0_m+t0_i+t0_p)*rf0/vthe
c$$$                    endif
                      else
                        eirene_mc_mapl_smo0to3(ix,iy,0,is,istrai)=
     &                   eirene_mc_mapl_smo0to3(ix,iy,0,is,istrai)+
     &                   (1.0_R8-rf0*ua(ix,iy,is)/vthe)*t0_a
                        eirene_mc_mmpl_smo0to3(ix,iy,0,is,istrai)=
     &                   eirene_mc_mmpl_smo0to3(ix,iy,0,is,istrai)+
     &                   (1.0_R8-rf0*ua(ix,iy,is)/vthe)*t0_m
                        eirene_mc_mipl_smo0to3(ix,iy,0,is,istrai)=
     &                   eirene_mc_mipl_smo0to3(ix,iy,0,is,istrai)+
     &                   (1.0_R8-rf0*ua(ix,iy,is)/vthe)*t0_i
                        eirene_mc_cppv_smo0to3(ix,iy,0,is,istrai)=
     &                   eirene_mc_cppv_smo0to3(ix,iy,0,is,istrai)+
     &                   (1.0_R8-rf0*ua(ix,iy,is)/vthe)*t0_p
                        eirene_mc_mapl_smo0to3(ix,iy,1,is,istrai)=
     &                   eirene_mc_mapl_smo0to3(ix,iy,1,is,istrai)+
     &                   rf0*t0_a/vthe
                        eirene_mc_mmpl_smo0to3(ix,iy,1,is,istrai)=
     &                   eirene_mc_mmpl_smo0to3(ix,iy,1,is,istrai)+
     &                   rf0*t0_m/vthe
                        eirene_mc_mipl_smo0to3(ix,iy,1,is,istrai)=
     &                   eirene_mc_mipl_smo0to3(ix,iy,1,is,istrai)+
     &                   rf0*t0_i/vthe
                        eirene_mc_cppv_smo0to3(ix,iy,1,is,istrai)=
     &                   eirene_mc_cppv_smo0to3(ix,iy,1,is,istrai)+
     &                   rf0*t0_p/vthe
c$$$                    if ((ix.eq.90).and.(iy.eq.20)) then !DJM DEBUG
c$$$                      write(*,*) 'smo0',
c$$$     &                           t0+RF0*ABS(t0)*UA(IX,IY,IS)/VTHE,
c$$$     &              (t0_a+t0_m+t0_i+t0_p)*(1.0_R8-rf0*ua(ix,iy,is)/vthe)
c$$$                      write(*,*) 'smo1',-RF0*ABS(t0)/VTHE,
c$$$     &                          (t0_a+t0_m+t0_i+t0_p)*rf0/vthe
c$$$                    endif
                      endif
                    endif
                  elseif (linearisation.eq.0) then !}{
                    smo0(IX,IY,0,IS)=smo0(IX,IY,0,IS)+t0
                    if (balance_netcdf.ne.0) then
                      eirene_mc_mapl_smo0to3(ix,iy,0,is,istrai)= !djm Jan2017
     &                 eirene_mc_mapl_smo0to3(ix,iy,0,is,istrai)+t0_a
                      eirene_mc_mmpl_smo0to3(ix,iy,0,is,istrai)=
     &                 eirene_mc_mmpl_smo0to3(ix,iy,0,is,istrai)+t0_m
                      eirene_mc_mipl_smo0to3(ix,iy,0,is,istrai)=
     &                 eirene_mc_mipl_smo0to3(ix,iy,0,is,istrai)+t0_i
                      eirene_mc_cppv_smo0to3(ix,iy,0,is,istrai)=
     &                 eirene_mc_cppv_smo0to3(ix,iy,0,is,istrai)+t0_p
                    endif
                  end if !}
                  SSMO(istrai,IS)=SSMO(istrai,IS)+t0
c*** calculate the integral sources for nregion regions, 3 instants:
c*** before correction (r), after correction (c), after iteration (i)
                  if (main_call.and.mod(ncall,eirene_mod).eq.0) then
                    src_mo(region(ix,iy,0),1,b2espcr(is))=
     .               src_mo(region(ix,iy,0),1,b2espcr(is))+t0
                    src_mo(region(ix,iy,0),2,b2espcr(is))=
     .               src_mo(region(ix,iy,0),2,b2espcr(is))+t0
                    cor_mo(ix,iy,b2espcr(is))=
     .               cor_mo(ix,iy,b2espcr(is))+t0
                    src_mo(0,1,b2espcr(is))=src_mo(0,1,b2espcr(is))+t0
                    src_mo(0,2,b2espcr(is))=src_mo(0,2,b2espcr(is))+t0
                    src_ms(region(ix,iy,0),1)=
     .               src_ms(region(ix,iy,0),1)+t0
                    src_ms(region(ix,iy,0),2)=
     .               src_ms(region(ix,iy,0),2)+t0
                    src_ms(0,1)=src_ms(0,1)+t0
                    src_ms(0,2)=src_ms(0,2)+t0
                  else
                    src_mo(region(ix,iy,0),3,b2espcr(is))=
     .               src_mo(region(ix,iy,0),3,b2espcr(is))+t0
                    src_ms(region(ix,iy,0),3)=
     .               src_ms(region(ix,iy,0),3)+t0
                    src_mo(0,3,b2espcr(is))=src_mo(0,3,b2espcr(is))+t0
                    src_ms(0,3)=src_ms(0,3)+t0
                  endif
                end if !}
              end do !}  ! ix
            end do !}  ! iy
          end if !}
        end do !}  ! is
        if (balance_netcdf.ne.0) then !{
          rf0=RF_NEUT(1)
          do is=1,nnatmi !{
            do iy=-1,ny !{
              do ix=-1,nx !{
                if(leftix(ix,iy).ne.-2 .and. rightix(ix,iy).ne.nx+1.and.
     &             bottomiy(ix,iy).ne.-2 .and. topiy(ix,iy).ne.ny+1.and.
     &             region(ix,iy,0).ne.0) then !{
                  t0_a = sna_av_paat(ix+1,iy+1,is,istrai)*
     *                   tflux(istrai)*
     *                   recyc_hlp*flux_scale(istrai)/
     *                   dble(ih_end)
                  t0_m = sna_av_pmat(ix+1,iy+1,is,istrai)*
     *                   tflux(istrai)*
     *                   recyc_hlp*flux_scale(istrai)/
     *                   dble(ih_end)
                  t0_i = sna_av_piat(ix+1,iy+1,is,istrai)*
     *                   tflux(istrai)*
     *                   recyc_hlp*flux_scale(istrai)/
     *                   dble(ih_end)
                  eirene_mc_paat_sna_bal(ix,iy,is,istrai)=
     *                   eirene_mc_paat_sna_bal(ix,iy,is,istrai)+t0_a
                  eirene_mc_pmat_sna_bal(ix,iy,is,istrai)=
     *                   eirene_mc_pmat_sna_bal(ix,iy,is,istrai)+t0_m
                  eirene_mc_piat_sna_bal(ix,iy,is,istrai)=
     *                   eirene_mc_piat_sna_bal(ix,iy,is,istrai)+t0_i
                end if !}
              end do !}  ! ix
            end do !}  ! iy
          end do !}  ! natmi
          do is=1,nnmoli !{
            do iy=-1,ny !{
              do ix=-1,nx !{
                if(leftix(ix,iy).ne.-2 .and. rightix(ix,iy).ne.nx+1.and.
     &             bottomiy(ix,iy).ne.-2 .and. topiy(ix,iy).ne.ny+1.and.
     &             region(ix,iy,0).ne.0) then !{
                  t0_a = sna_av_paml(ix+1,iy+1,is,istrai)*
     *                   tflux(istrai)*
     *                   recyc_hlp*flux_scale(istrai)/
     *                   dble(ih_end)
                  t0_m = sna_av_pmml(ix+1,iy+1,is,istrai)*
     *                   tflux(istrai)*
     *                   recyc_hlp*flux_scale(istrai)/
     *                   dble(ih_end)
                  t0_i = sna_av_piml(ix+1,iy+1,is,istrai)*
     *                   tflux(istrai)*
     *                   recyc_hlp*flux_scale(istrai)/
     *                   dble(ih_end)
                  eirene_mc_paml_sna_bal(ix,iy,is,istrai)=
     *                   eirene_mc_paml_sna_bal(ix,iy,is,istrai)+t0_a
                  eirene_mc_pmml_sna_bal(ix,iy,is,istrai)=
     *                   eirene_mc_pmml_sna_bal(ix,iy,is,istrai)+t0_m
                  eirene_mc_piml_sna_bal(ix,iy,is,istrai)=
     *                   eirene_mc_piml_sna_bal(ix,iy,is,istrai)+t0_i
                end if !}
              end do !}  ! ix
            end do !}  ! iy
          end do !}  ! nmoli
          do is=1,nnioni !{
            do iy=-1,ny !{
              do ix=-1,nx !{
                if(leftix(ix,iy).ne.-2 .and. rightix(ix,iy).ne.nx+1.and.
     &             bottomiy(ix,iy).ne.-2 .and. topiy(ix,iy).ne.ny+1.and.
     &             region(ix,iy,0).ne.0) then !{
                  t0_a = sna_av_paio(ix+1,iy+1,is,istrai)*
     *                   tflux(istrai)*
     *                   recyc_hlp*flux_scale(istrai)/
     *                   dble(ih_end)
                  t0_m = sna_av_pmio(ix+1,iy+1,is,istrai)*
     *                   tflux(istrai)*
     *                   recyc_hlp*flux_scale(istrai)/
     *                   dble(ih_end)
                  t0_i = sna_av_piio(ix+1,iy+1,is,istrai)*
     *                   tflux(istrai)*
     *                   recyc_hlp*flux_scale(istrai)/
     *                   dble(ih_end)
                  eirene_mc_paio_sna_bal(ix,iy,is,istrai)=
     *                   eirene_mc_paio_sna_bal(ix,iy,is,istrai)+t0_a
                  eirene_mc_pmio_sna_bal(ix,iy,is,istrai)=
     *                   eirene_mc_pmio_sna_bal(ix,iy,is,istrai)+t0_m
                  eirene_mc_piio_sna_bal(ix,iy,is,istrai)=
     *                   eirene_mc_piio_sna_bal(ix,iy,is,istrai)+t0_i
                end if !}
              end do !}  ! ix
            end do !}  ! iy
          end do !}  ! nioni
        end if !}
        do iy=-1,ny !{
          do ix=-1,nx !{
            if(leftix(ix,iy).ne.-2 .and. rightix(ix,iy).ne.nx+1 .and.
     &         bottomiy(ix,iy).ne.-2 .and. topiy(ix,iy).ne.ny+1 .and.
     &         region(ix,iy,0).ne.0) then !{
* electron particle source
              rf0=RF_NEUT(1)
              if (balance_netcdf.ne.0) then
                t0_a = sne_av_pael(ix+1,iy+1,istrai)*tflux(istrai)* !djm Jan2017
     *           recyc_hlp*flux_scale(istrai)*sclstra(0,istrai)/
     *           dble(ih_end)
                t0_m = sne_av_pmel(ix+1,iy+1,istrai)*tflux(istrai)*
     *           recyc_hlp*flux_scale(istrai)*sclstra(0,istrai)/
     *           dble(ih_end)
              endif
              if (linearisation.eq.1) then !{
                if (balance_netcdf.ne.0) then
cdjm Jan2017 Calculate linearised components of source
                  if (t0_a.ge.0.0_R8) then
                    eirene_mc_pael_sne0to1(ix,iy,0,istrai)=
     &               eirene_mc_pael_sne0to1(ix,iy,0,istrai)+
     &               (1.0_R8+rf0)*t0_a
                    eirene_mc_pael_sne0to1(ix,iy,1,istrai)=
     &               eirene_mc_pael_sne0to1(ix,iy,1,istrai)-
     &               rf0*t0_a/ne(ix,iy)
                  else
                    eirene_mc_pael_sne0to1(ix,iy,0,istrai)=
     &               eirene_mc_pael_sne0to1(ix,iy,0,istrai)-
     &               rf0*t0_a
                    eirene_mc_pael_sne0to1(ix,iy,1,istrai)=
     &               eirene_mc_pael_sne0to1(ix,iy,1,istrai)+
     &               (1.0_R8+rf0)*t0_a/ne(ix,iy)
                  endif
                  if (t0_m.ge.0.0_R8) then
                    eirene_mc_pmel_sne0to1(ix,iy,0,istrai)=
     &               eirene_mc_pmel_sne0to1(ix,iy,0,istrai)+
     &               (1.0_R8+rf0)*t0_m
                    eirene_mc_pmel_sne0to1(ix,iy,1,istrai)=
     &               eirene_mc_pmel_sne0to1(ix,iy,1,istrai)-
     &               rf0*t0_m/ne(ix,iy)
                  else
                    eirene_mc_pmel_sne0to1(ix,iy,0,istrai)=
     &               eirene_mc_pmel_sne0to1(ix,iy,0,istrai)-
     &               rf0*t0_m
                    eirene_mc_pmel_sne0to1(ix,iy,1,istrai)=
     &               eirene_mc_pmel_sne0to1(ix,iy,1,istrai)+
     &               (1.0_R8+rf0)*t0_m/ne(ix,iy)
                  endif
                endif
              elseif (linearisation.eq.0) then !}{
                if (balance_netcdf.ne.0) then
                  eirene_mc_pael_sne0to1(ix,iy,0,istrai)= !djm Jan2017
     &             eirene_mc_pael_sne0to1(ix,iy,0,istrai)+t0_a
                  eirene_mc_pmel_sne0to1(ix,iy,0,istrai)=
     &             eirene_mc_pmel_sne0to1(ix,iy,0,istrai)+t0_m
                endif
              end if !}
* electron energy
              rf0=RF_NEUT(3)
csw 26feb2014 averaging
c              t0=see(ix+1,iy+1,istrai)*tflux(istrai)*
c     *         recyc_hlp*flux_scale(istrai)*sclstra(0,istrai)
              t0=see_av(ix+1,iy+1,istrai)*tflux(istrai)*
     *         recyc_hlp*flux_scale(istrai)*sclstra(0,istrai)/
     *         dble(ih_end)
              if (balance_netcdf.ne.0) then
                t0_a = see_av_eael(ix+1,iy+1,istrai)*tflux(istrai)* !djm Jan2017
     *           recyc_hlp*flux_scale(istrai)*sclstra(0,istrai)/
     *           dble(ih_end)
                t0_m = see_av_emel(ix+1,iy+1,istrai)*tflux(istrai)*
     *           recyc_hlp*flux_scale(istrai)*sclstra(0,istrai)/
     *           dble(ih_end)
                t0_i = see_av_eiel(ix+1,iy+1,istrai)*tflux(istrai)*
     *           recyc_hlp*flux_scale(istrai)*sclstra(0,istrai)/
     *           dble(ih_end)
                t0_p = see_av_epel(ix+1,iy+1,istrai)*tflux(istrai)*
     *           recyc_hlp*flux_scale(istrai)*sclstra(0,istrai)/
     *           dble(ih_end)
              endif

              if (linearisation.eq.1) then !{
                she0(ix,iy,0)=she0(ix,iy,0)+max(0.0_R8,t0)+abs(t0)*rf0
                she0(ix,iy,1)=she0(ix,iy,1)+
     1           (min(0.0_R8,t0)-abs(t0)*rf0)/te(ix,iy)
                if (balance_netcdf.ne.0) then
cdjm Jan2017 Calculate linearised components of source
                  if (t0.ge.0.0_R8) then
                    eirene_mc_eael_she0to3(ix,iy,0,istrai)=
     &               eirene_mc_eael_she0to3(ix,iy,0,istrai)+
     &               (1.0_R8+rf0)*t0_a
                    eirene_mc_emel_she0to3(ix,iy,0,istrai)=
     &               eirene_mc_emel_she0to3(ix,iy,0,istrai)+
     &               (1.0_R8+rf0)*t0_m
                    eirene_mc_eiel_she0to3(ix,iy,0,istrai)=
     &               eirene_mc_eiel_she0to3(ix,iy,0,istrai)+
     &               (1.0_R8+rf0)*t0_i
                    eirene_mc_epel_she0to3(ix,iy,0,istrai)=
     &               eirene_mc_epel_she0to3(ix,iy,0,istrai)+
     &               (1.0_R8+rf0)*t0_p
                    eirene_mc_eael_she0to3(ix,iy,1,istrai)=
     &               eirene_mc_eael_she0to3(ix,iy,1,istrai)-
     &               t0_a*rf0/te(ix,iy)
                    eirene_mc_emel_she0to3(ix,iy,1,istrai)=
     &               eirene_mc_emel_she0to3(ix,iy,1,istrai)-
     &               t0_m*rf0/te(ix,iy)
                    eirene_mc_eiel_she0to3(ix,iy,1,istrai)=
     &               eirene_mc_eiel_she0to3(ix,iy,1,istrai)-
     &               t0_i*rf0/te(ix,iy)
                    eirene_mc_epel_she0to3(ix,iy,1,istrai)=
     &               eirene_mc_epel_she0to3(ix,iy,1,istrai)-
     &               t0_p*rf0/te(ix,iy)
c$$$                    if ((ix.eq.90).and.(iy.eq.20)) then !DJM DEBUG
c$$$                      write(*,*) 'she0',max(0.0_R8,t0)+abs(t0)*rf0,
c$$$     &              (t0_a+t0_m+t0_i+t0_p)*(1.0_R8+rf0)
c$$$                      write(*,*) 'she1',
c$$$     &                   (min(0.0_R8,t0)-abs(t0)*rf0)/te(ix,iy),
c$$$     &                         -(t0_a+t0_m+t0_i+t0_p)*rf0/te(ix,iy)
c$$$                    endif
                  else
                    eirene_mc_eael_she0to3(ix,iy,0,istrai)=
     &               eirene_mc_eael_she0to3(ix,iy,0,istrai)-
     &               rf0*t0_a
                    eirene_mc_emel_she0to3(ix,iy,0,istrai)=
     &               eirene_mc_emel_she0to3(ix,iy,0,istrai)-
     &               rf0*t0_m
                    eirene_mc_eiel_she0to3(ix,iy,0,istrai)=
     &               eirene_mc_eiel_she0to3(ix,iy,0,istrai)-
     &               rf0*t0_i
                    eirene_mc_epel_she0to3(ix,iy,0,istrai)=
     &               eirene_mc_epel_she0to3(ix,iy,0,istrai)-
     &               rf0*t0_p
                    eirene_mc_eael_she0to3(ix,iy,1,istrai)=
     &               eirene_mc_eael_she0to3(ix,iy,1,istrai)+
     &               (1.0_R8+rf0)*t0_a/te(ix,iy)
                    eirene_mc_emel_she0to3(ix,iy,1,istrai)=
     &               eirene_mc_emel_she0to3(ix,iy,1,istrai)+
     &               (1.0_R8+rf0)*t0_m/te(ix,iy)
                    eirene_mc_eiel_she0to3(ix,iy,1,istrai)=
     &               eirene_mc_eiel_she0to3(ix,iy,1,istrai)+
     &               (1.0_R8+rf0)*t0_i/te(ix,iy)
                    eirene_mc_epel_she0to3(ix,iy,1,istrai)=
     &               eirene_mc_epel_she0to3(ix,iy,1,istrai)+
     &               (1.0_R8+rf0)*t0_p/te(ix,iy)
c$$$                    if ((ix.eq.90).and.(iy.eq.20)) then !DJM DEBUG
c$$$                      write(*,*) 'she0',max(0.0_R8,t0)+abs(t0)*rf0,
c$$$     &              -(t0_a+t0_m+t0_i+t0_p)*rf0
c$$$                      write(*,*) 'she1',
c$$$     &                   (min(0.0_R8,t0)-abs(t0)*rf0)/te(ix,iy),
c$$$     &                    (t0_a+t0_m+t0_i+t0_p)*(1.0_R8+rf0)/te(ix,iy)
c$$$                    endif
                  endif
                endif
              elseif (linearisation.eq.0) then !}{
                she0(ix,iy,0)=she0(ix,iy,0)+t0
                if (balance_netcdf.ne.0) then
                  eirene_mc_eael_she0to3(ix,iy,0,istrai)= !djm Jan2017
     &             eirene_mc_eael_she0to3(ix,iy,0,istrai)+t0_a
                  eirene_mc_emel_she0to3(ix,iy,0,istrai)=
     &             eirene_mc_emel_she0to3(ix,iy,0,istrai)+t0_m
                  eirene_mc_eiel_she0to3(ix,iy,0,istrai)=
     &             eirene_mc_eiel_she0to3(ix,iy,0,istrai)+t0_i
                  eirene_mc_epel_she0to3(ix,iy,0,istrai)=
     &             eirene_mc_epel_she0to3(ix,iy,0,istrai)+t0_p
                endif
              end if !}
              ssee(istrai)=ssee(istrai)+t0
c*** calculate the integral sources for nregion regions, 3 instants:
c*** before correction (r), after correction (c), after iteration (i)
              if (main_call.and.mod(ncall,eirene_mod).eq.0) then
                src_ee(region(ix,iy,0),1)=
     .           src_ee(region(ix,iy,0),1)+t0
                src_ee(region(ix,iy,0),2)=
     .           src_ee(region(ix,iy,0),2)+t0
                cor_ee(ix,iy)=cor_ee(ix,iy)+t0
                src_es(region(ix,iy,0),1)=
     .           src_es(region(ix,iy,0),1)+t0
                src_es(region(ix,iy,0),2)=
     .           src_es(region(ix,iy,0),2)+t0
                src_ee(0,1)=src_ee(0,1)+t0
                src_ee(0,2)=src_ee(0,2)+t0
                src_es(0,1)=src_es(0,1)+t0
                src_es(0,2)=src_es(0,2)+t0
              else
                src_ee(region(ix,iy,0),3)=
     .           src_ee(region(ix,iy,0),3)+t0
                src_es(region(ix,iy,0),3)=
     .           src_es(region(ix,iy,0),3)+t0
                src_ee(0,3)=src_ee(0,3)+t0
                src_es(0,3)=src_es(0,3)+t0
              endif
csw 21feb2012, add l_neutrad (SOLPS4.3)
              if(l_neutrad.lt.0) then
                radarr(ix,iy)=radarr(ix,iy)+t0
                rad_neutr=rad_neutr+t0
                rdneureg(region(ix,iy,0))=rdneureg(region(ix,iy,0))+t0
                recycereg(region(ix,iy,0))=recycereg(region(ix,iy,0))
     .                                     +t0
              endif
              sseereg(region(ix,iy,0))=sseereg(region(ix,iy,0))+t0
* ion energy
              rf0=RF_NEUT(4)
csw 26feb2014 averaging
c              t0=sei(ix+1,iy+1,istrai)*tflux(istrai)*
c     1         recyc_hlp*flux_scale(istrai)*sclstra(0,istrai)
              t0=sei_av(ix+1,iy+1,istrai)*tflux(istrai)*
     1         recyc_hlp*flux_scale(istrai)*sclstra(0,istrai)/
     1         dble(ih_end)
              if (balance_netcdf.ne.0) then
                t0_a=sei_av_eapl(ix+1,iy+1,istrai)*tflux(istrai)* !djm Jan2017
     1           recyc_hlp*flux_scale(istrai)*sclstra(0,istrai)/
     1           dble(ih_end)
                t0_m=sei_av_empl(ix+1,iy+1,istrai)*tflux(istrai)*
     1           recyc_hlp*flux_scale(istrai)*sclstra(0,istrai)/
     1           dble(ih_end)
                t0_i=sei_av_eipl(ix+1,iy+1,istrai)*tflux(istrai)*
     1           recyc_hlp*flux_scale(istrai)*sclstra(0,istrai)/
     1           dble(ih_end)
                t0_p=sei_av_eppl(ix+1,iy+1,istrai)*tflux(istrai)*
     1           recyc_hlp*flux_scale(istrai)*sclstra(0,istrai)/
     1           dble(ih_end)
              endif

              ssei(istrai)=ssei(istrai)+t0
c ? RADARR ?
csw 06mar2012 add l_neutrad
              if(l_neutrad.lt.0) then
                recycireg(region(ix,iy,0))=recycireg(region(ix,iy,0))
     .                                     +t0
              endif
              sseireg(region(ix,iy,0))=sseireg(region(ix,iy,0))+t0
c*** calculate the integral sources for nregion regions, 3 instants:
c*** before correction (r), after correction (c), after iteration (i)
              if (main_call.and.mod(ncall,eirene_mod).eq.0) then
                src_ei(region(ix,iy,0),1)=
     .           src_ei(region(ix,iy,0),1)+t0
                src_ei(region(ix,iy,0),2)=
     .           src_ei(region(ix,iy,0),2)+t0
                cor_ei(ix,iy)=cor_ei(ix,iy)+t0
                src_es(region(ix,iy,0),1)=
     .           src_es(region(ix,iy,0),1)+t0
                src_es(region(ix,iy,0),2)=
     .           src_es(region(ix,iy,0),2)+t0
                src_ei(0,1)=src_ei(0,1)+t0
                src_es(0,1)=src_es(0,1)+t0
                src_ei(0,2)=src_ei(0,2)+t0
                src_es(0,2)=src_es(0,2)+t0
              else
                src_ei(region(ix,iy,0),3)=
     .           src_ei(region(ix,iy,0),3)+t0
                src_es(region(ix,iy,0),3)=
     .           src_es(region(ix,iy,0),3)+t0
                src_ei(0,3)=src_ei(0,3)+t0
                src_es(0,3)=src_es(0,3)+t0
              endif
              t1=0.0_R8
              if (balance_netcdf.ne.0) then
                t1_apl=0.0_R8 !djm Jan2017
                t1_mpl=0.0_R8
                t1_ipl=0.0_R8
                t1_ppl=0.0_R8
              endif
cxpb Correction of Eirene energy sources to B2 heat sources
cxpb Only necessary for SOLPS5.0 equation formulation
              if (BoRiS.eq.0.0_R8.and.no_sseicorr.eq.0) then !{
                do is=0,ns-1 !{
                  if(use_recyceir.gt.0) then
                    recyctmp=recyceir(istrai)
                  else
                    recyctmp=recyc(is,istrai)
                  endif
                  if(.not.is_neutral(is)) then !{
                    ifl=fluids_list(is)
                    j=b2eatcr(is)
                    if (smo_tot(ix+1,iy+1,ifl).ne.0.0_R8 .and.
     .                  sni_tot(ix+1,iy+1,ifl).ne.0.0_R8) then
                      dua = (smo_tot(ix+1,iy+1,ifl)*
     .                       dtim*dtmo(is,region(ix,iy,0))/(am(is)*mp) -
     .                       ua(ix,iy,is)*sni_tot(ix+1,iy+1,ifl)*
     .                       dtim*dtco(is,region(ix,iy,0)) )/
     .                      (na(ix,iy,is)*vol(ix,iy) +
     .                       sni_tot(ix+1,iy+1,ifl)*
     .                       dtim*dtco(is,region(ix,iy,0)))
                    else
                      dua = 0.0_R8
                    end if

csw 06mar2012 momentum source
csw 26feb2014 averaging
c                    t2=smo(ix+1,iy+1,ifl,istrai)*tflux(istrai)*
c     *               recyctmp*flux_scale(istrai)
                    t2=smo_av(ix+1,iy+1,ifl,istrai)*tflux(istrai)*
     *               recyctmp*flux_scale(istrai)/
     *               dble(ih_end)
                    if (balance_netcdf.ne.0) then
                      t2_apl= !djm Jan2017
     *                 smo_av_mapl(ix+1,iy+1,ifl,istrai)*tflux(istrai)*
     *                 recyctmp*flux_scale(istrai)/
     *                 dble(ih_end)
                      t2_mpl=
     *                 smo_av_mmpl(ix+1,iy+1,ifl,istrai)*tflux(istrai)*
     *                 recyctmp*flux_scale(istrai)/
     *                 dble(ih_end)
                      t2_ipl=
     *                 smo_av_mipl(ix+1,iy+1,ifl,istrai)*tflux(istrai)*
     *                 recyctmp*flux_scale(istrai)/
     *                 dble(ih_end)
                      t2_ppl=
     *                 smo_av_cppv(ix+1,iy+1,ifl,istrai)*tflux(istrai)*
     *                 recyctmp*flux_scale(istrai)/
     *                 dble(ih_end)
                    endif

                    t2=t2*sclstra(j,istrai)
                    if (balance_netcdf.ne.0) then
                      t2_apl=t2_apl*sclstra(j,istrai) !djm Jan2017
                      t2_mpl=t2_mpl*sclstra(j,istrai)
                      t2_ipl=t2_ipl*sclstra(j,istrai)
                      t2_ppl=t2_ppl*sclstra(j,istrai)
                    endif

csw 06mar2012 particle source
csw 26feb2014 averaging
c                    t3=sni(ix+1,iy+1,ifl,istrai)*tflux(istrai)*
c     *               recyctmp*flux_scale(istrai)
                    t3=sni_av(ix+1,iy+1,ifl,istrai)*tflux(istrai)*
     *               recyctmp*flux_scale(istrai)/
     *               dble(ih_end)
                    if (balance_netcdf.ne.0) then
                      t3_apl= !djm Jan2017
     *                 sni_av_papl(ix+1,iy+1,ifl,istrai)*tflux(istrai)*
     *                 recyctmp*flux_scale(istrai)/
     *                 dble(ih_end)
                      t3_mpl=
     *                 sni_av_pmpl(ix+1,iy+1,ifl,istrai)*tflux(istrai)*
     *                 recyctmp*flux_scale(istrai)/
     *                 dble(ih_end)
                      t3_ipl=
     *                 sni_av_pipl(ix+1,iy+1,ifl,istrai)*tflux(istrai)*
     *                 recyctmp*flux_scale(istrai)/
     *                 dble(ih_end)
                      t3_ppl=
     *                 sni_av_pppl(ix+1,iy+1,ifl,istrai)*tflux(istrai)*
     *                 recyctmp*flux_scale(istrai)/
     *                 dble(ih_end)
                    endif

csw 26mar2013 (for TESTING)
                    if(iscltst .eq. 0) then
                      if(crcstra(istrai).ne.'V' .or. t3.gt.0.0_R8)
     -                 t3=t3*sclalph(j,istrai)
                       if (balance_netcdf.ne.0) then
                         t3_apl=t3_apl*sclalph(j,istrai) !djm Jan2017
                         t3_mpl=t3_mpl*sclalph(j,istrai)
                         t3_ipl=t3_ipl*sclalph(j,istrai)
                         t3_ppl=t3_ppl*sclalph(j,istrai)
                       endif
                    else
                      if(crcstra(istrai).ne.'V' .or. t3.gt.0.0_R8)
     -                 t3=t3*sclstra(j,istrai)
                       if (balance_netcdf.ne.0) then
                         t3_apl=t3_apl*sclstra(j,istrai) !djm Jan2017
                         t3_mpl=t3_mpl*sclstra(j,istrai)
                         t3_ipl=t3_ipl*sclstra(j,istrai)
                         t3_ppl=t3_ppl*sclstra(j,istrai)
                       endif
                    endif

csw 06mar2012 total correction
                    t1=t1-ua(ix,iy,is)*t2+kinrgy(ix,iy,is)*t3
                    if (balance_netcdf.ne.0) then
                      t1_apl=t1_apl-ua(ix,iy,is)*t2_apl+ !djm Jan2017
     &                       kinrgy(ix,iy,is)*t3_apl
                      t1_mpl=t1_mpl-ua(ix,iy,is)*t2_mpl+
     &                       kinrgy(ix,iy,is)*t3_mpl
                      t1_ipl=t1_ipl-ua(ix,iy,is)*t2_ipl+
     &                       kinrgy(ix,iy,is)*t3_ipl
                      t1_ppl=t1_ppl-ua(ix,iy,is)*t2_ppl+
     &                       kinrgy(ix,iy,is)*t3_ppl
                    endif
cxpb 13feb2015 additional correction with second-order effects for large sources
cxpb 12oct2016 multiply na term by fractional stratum weight to avoid double-counting
                    if (smo_tot(ix+1,iy+1,ifl).ne.0.0_R8.and.
     .                  t3.ne.0.0_R8 .and.
     .                  stab_corr.eq.1) then
                      t4=(na(ix,iy,is)*vol(ix,iy)*
     .                    (t2/smo_tot(ix+1,iy+1,ifl))+
     .                    t3*dtim*dtco(is,region(ix,iy,0)))/
     .                      (dtim*dtei(region(ix,iy,0)))*
     .                    0.5_R8*mp*am(is)*dua**2
                      if (balance_netcdf.ne.0) then
                        t4_apl=(na(ix,iy,is)*vol(ix,iy)* !djm Jan2017
     .                      (t2_apl/smo_tot(ix+1,iy+1,ifl))+
     .                      t3_apl*dtim*dtco(is,region(ix,iy,0)))/
     .                        (dtim*dtei(region(ix,iy,0)))*
     .                      0.5_R8*mp*am(is)*dua**2
                        t4_mpl=(na(ix,iy,is)*vol(ix,iy)*
     .                      (t2_mpl/smo_tot(ix+1,iy+1,ifl))+
     .                      t3_mpl*dtim*dtco(is,region(ix,iy,0)))/
     .                        (dtim*dtei(region(ix,iy,0)))*
     .                      0.5_R8*mp*am(is)*dua**2
                        t4_ipl=(na(ix,iy,is)*vol(ix,iy)*
     .                      (t2_ipl/smo_tot(ix+1,iy+1,ifl))+
     .                      t3_ipl*dtim*dtco(is,region(ix,iy,0)))/
     .                        (dtim*dtei(region(ix,iy,0)))*
     .                      0.5_R8*mp*am(is)*dua**2
                        t4_ppl=(na(ix,iy,is)*vol(ix,iy)*
     .                      (t2_ppl/smo_tot(ix+1,iy+1,ifl))+
     .                      t3_ppl*dtim*dtco(is,region(ix,iy,0)))/
     .                        (dtim*dtei(region(ix,iy,0)))*
     .                      0.5_R8*mp*am(is)*dua**2
                      endif
                    else
                      t4=0.0_R8
                      if (balance_netcdf.ne.0) then
                        t4_apl=0.0_R8 !djm Jan2017
                        t4_mpl=0.0_R8
                        t4_ipl=0.0_R8
                        t4_ppl=0.0_R8
                      endif
                    endif
                    t1=t1-t4
                    sseistab(istrai)=sseistab(istrai)-t4
                    if (balance_netcdf.ne.0) then
                      t1_apl=t1_apl-t4_apl !djm Jan2017
                      t1_mpl=t1_mpl-t4_mpl
                      t1_ipl=t1_ipl-t4_ipl
                      t1_ppl=t1_ppl-t4_ppl
                    endif
                  end if !}
                end do !} !is
              end if !}
              sseicorr(istrai)=sseicorr(istrai)+t1
              t0=t0+t1
              if (balance_netcdf.ne.0) then
                t0_a=t0_a+t1_apl !djm Jan2017
                t0_m=t0_m+t1_mpl
                t0_i=t0_i+t1_ipl
                t0_p=t0_p+t1_ppl
              endif
              if (linearisation.eq.1) then !{
                shi0(ix,iy,0)=shi0(ix,iy,0)+max(0.0_R8,t0)+abs(t0)*rf0
                shi0(ix,iy,1)=shi0(ix,iy,1)+
     1           (min(0.0_R8,t0)-abs(t0)*rf0)/ti(ix,iy)
                if (balance_netcdf.ne.0) then
cdjm Jan2017 Calculate linearised components of source
                  if (t0.ge.0.0_R8) then
                    eirene_mc_eapl_shi0to3(ix,iy,0,istrai)=
     &               eirene_mc_eapl_shi0to3(ix,iy,0,istrai)+
     &               (1.0_R8+rf0)*t0_a
                    eirene_mc_empl_shi0to3(ix,iy,0,istrai)=
     &               eirene_mc_empl_shi0to3(ix,iy,0,istrai)+
     &               (1.0_R8+rf0)*t0_m
                    eirene_mc_eipl_shi0to3(ix,iy,0,istrai)=
     &               eirene_mc_eipl_shi0to3(ix,iy,0,istrai)+
     &               (1.0_R8+rf0)*t0_i
                    eirene_mc_eppl_shi0to3(ix,iy,0,istrai)=
     &               eirene_mc_eppl_shi0to3(ix,iy,0,istrai)+
     &               (1.0_R8+rf0)*t0_p
                    eirene_mc_eapl_shi0to3(ix,iy,1,istrai)=
     &               eirene_mc_eapl_shi0to3(ix,iy,1,istrai)-
     &               t0_a*rf0/ti(ix,iy)
                    eirene_mc_empl_shi0to3(ix,iy,1,istrai)=
     &               eirene_mc_empl_shi0to3(ix,iy,1,istrai)-
     &               t0_m*rf0/ti(ix,iy)
                    eirene_mc_eipl_shi0to3(ix,iy,1,istrai)=
     &               eirene_mc_eipl_shi0to3(ix,iy,1,istrai)-
     &               t0_i*rf0/ti(ix,iy)
                    eirene_mc_eppl_shi0to3(ix,iy,1,istrai)=
     &               eirene_mc_eppl_shi0to3(ix,iy,1,istrai)-
     &               t0_p*rf0/ti(ix,iy)
c$$$                    if ((ix.eq.90).and.(iy.eq.20)) then !DJM DEBUG
c$$$                      write(*,*) 'shi0',max(0.0_R8,t0)+abs(t0)*rf0,
c$$$     &              (t0_a+t0_m+t0_i+t0_p)*(1.0_R8+rf0)
c$$$                      write(*,*) 'shi1',
c$$$     &                   (min(0.0_R8,t0)-abs(t0)*rf0)/ti(ix,iy),
c$$$     &                         -(t0_a+t0_m+t0_i+t0_p)*rf0/ti(ix,iy)
c$$$                    endif
                  else
                    eirene_mc_eapl_shi0to3(ix,iy,0,istrai)=
     &               eirene_mc_eapl_shi0to3(ix,iy,0,istrai)-
     &               rf0*t0_a
                    eirene_mc_empl_shi0to3(ix,iy,0,istrai)=
     &               eirene_mc_empl_shi0to3(ix,iy,0,istrai)-
     &               rf0*t0_m
                    eirene_mc_eipl_shi0to3(ix,iy,0,istrai)=
     &               eirene_mc_eipl_shi0to3(ix,iy,0,istrai)-
     &               rf0*t0_i
                    eirene_mc_eppl_shi0to3(ix,iy,0,istrai)=
     &               eirene_mc_eppl_shi0to3(ix,iy,0,istrai)-
     &               rf0*t0_p
                    eirene_mc_eapl_shi0to3(ix,iy,1,istrai)=
     &               eirene_mc_eapl_shi0to3(ix,iy,1,istrai)+
     &               (1.0_R8+rf0)*t0_a/ti(ix,iy)
                    eirene_mc_empl_shi0to3(ix,iy,1,istrai)=
     &               eirene_mc_empl_shi0to3(ix,iy,1,istrai)+
     &               (1.0_R8+rf0)*t0_m/ti(ix,iy)
                    eirene_mc_eipl_shi0to3(ix,iy,1,istrai)=
     &               eirene_mc_eipl_shi0to3(ix,iy,1,istrai)+
     &               (1.0_R8+rf0)*t0_i/ti(ix,iy)
                    eirene_mc_eppl_shi0to3(ix,iy,1,istrai)=
     &               eirene_mc_eppl_shi0to3(ix,iy,1,istrai)+
     &               (1.0_R8+rf0)*t0_p/ti(ix,iy)
c$$$                    if ((ix.eq.90).and.(iy.eq.20)) then !DJM DEBUG
c$$$                      write(*,*) 'shi0',max(0.0_R8,t0)+abs(t0)*rf0,
c$$$     &              -(t0_a+t0_m+t0_i+t0_p)*rf0
c$$$                      write(*,*) 'shi1',
c$$$     &                   (min(0.0_R8,t0)-abs(t0)*rf0)/ti(ix,iy),
c$$$     &                    (t0_a+t0_m+t0_i+t0_p)*(1.0_R8+rf0)/ti(ix,iy)
c$$$                    endif
                  endif
                endif
              elseif (linearisation.eq.0) then !}{
                shi0(ix,iy,0)=shi0(ix,iy,0)+t0
                if (balance_netcdf.ne.0) then !djm Jan2017
                  eirene_mc_eapl_shi0to3(ix,iy,0,istrai)=
     &             eirene_mc_eapl_shi0to3(ix,iy,0,istrai)+t0_a
                  eirene_mc_empl_shi0to3(ix,iy,0,istrai)=
     &             eirene_mc_empl_shi0to3(ix,iy,0,istrai)+t0_m
                  eirene_mc_eipl_shi0to3(ix,iy,0,istrai)=
     &             eirene_mc_eipl_shi0to3(ix,iy,0,istrai)+t0_i
                  eirene_mc_eppl_shi0to3(ix,iy,0,istrai)=
     &             eirene_mc_eppl_shi0to3(ix,iy,0,istrai)+t0_p
                endif
              end if !}
            end if !}
          end do !} ix
        end do !} iy
      end do !} istrai
c
c*** correct radiation losses for the molecule dissociation
c
csw 21feb2012 add test l_neutrad (SOLPS4.3)
      if(l_neutrad.lt.0) then
        do iy=-1,ny
          do ix=-1,nx
            if (leftix(ix,iy).ne.-2 .and. rightix(ix,iy).ne.nx+1 .and.
     &          topiy(ix,iy).ne.ny+1 .and. bottomiy(ix,iy).ne.-2 .and.
     &          region(ix,iy,0).ne.0) then
              k=region(ix,iy,0)
              do imol = 1, nnmoli
!pb  eneutrad dimensioned per stratum to allow for short cycling
!pb  only sum over strata used here
!                radarr(ix,iy)=radarr(ix,iy)-edissml(ix+1,iy+1,imol,1)
!                rad_neutr=rad_neutr-edissml(ix+1,iy+1,imol,1)
!                rdneureg(k)=rdneureg(k)-edissml(ix+1,iy+1,imol,1)
!                recycereg(k)=recycereg(k)-edissml(ix+1,iy+1,imol,1)
                radarr(ix,iy)=radarr(ix,iy)-edissml(ix+1,iy+1,imol,0)
                rad_neutr=rad_neutr-edissml(ix+1,iy+1,imol,0)
                rdneureg(k)=rdneureg(k)-edissml(ix+1,iy+1,imol,0)
                recycereg(k)=recycereg(k)-edissml(ix+1,iy+1,imol,0)
              end do
            end if
          end do
        end do
        write (hlp_frm,'(a,i2,a)') '(a,1p,',nnmoli,'(1x,e14.7))'
        if (dbg_eir_mc.ge.2 .and. main_call) write(*,hlp_frm)
     .   'EIRENE_MC. EDISSML : ',
     .   (SUM(edissml(1:nx,1:ny,imol,0)),imol=1,nnmoli)
      else

csw 21feb2012, calculate neutral radiation directly from EIRENE (SOLPS4.3)
c*** calculate radiation losses directly (at the moment only atoms + edissml)
c
        rad_neutr=0.0_R8
        rdneureg=0.0_R8
        recycereg=0.0_R8
        do iy=-1,ny
          do ix=-1,nx
            if (leftix(ix,iy).ne.-2 .and. rightix(ix,iy).ne.nx+1 .and.
     &          topiy(ix,iy).ne.ny+1 .and. bottomiy(ix,iy).ne.-2 .and.
     &          region(ix,iy,0).ne.0) then
              k=region(ix,iy,0)
              do iatm = 1, nnatmi
!pb  18.12.2013
!pb  eneutrad dimensioned per stratum to allow for short cycling
!pb  only sum over strata used here
!                radarr(ix,iy)=radarr(ix,iy)+eneutrad(ix+1,iy+1,iatm,1)
!                rad_neutr=rad_neutr+eneutrad(ix+1,iy+1,iatm,1)
!                rdneureg(k)=rdneureg(k)+eneutrad(ix+1,iy+1,iatm,1)
!                recycereg(k)=recycereg(k)+eneutrad(ix+1,iy+1,iatm,1)
!pb
                 radarr(ix,iy)=radarr(ix,iy)+eneutrad(ix+1,iy+1,iatm,0)
                 rad_neutr=rad_neutr+eneutrad(ix+1,iy+1,iatm,0)
                 rdneureg(k)=rdneureg(k)+eneutrad(ix+1,iy+1,iatm,0)
                 recycereg(k)=recycereg(k)+eneutrad(ix+1,iy+1,iatm,0)
              enddo
              do imol = 1, nnmoli
!pb  edissml dimensioned per stratum to allow for short cycling
!pb  only sum over strata used here
!                radarr(ix,iy)=radarr(ix,iy)-edissml(ix+1,iy+1,imol,1)
!                rad_neutr=rad_neutr-edissml(ix+1,iy+1,imol,1)
!                rdneureg(k)=rdneureg(k)-edissml(ix+1,iy+1,imol,1)
!                recycereg(k)=recycereg(k)-edissml(ix+1,iy+1,imol,1)
cxpb Replacing edissml contribution with emolrad
c$$$                radarr(ix,iy)=radarr(ix,iy)-edissml(ix+1,iy+1,imol,0)
c$$$                rad_neutr=rad_neutr-edissml(ix+1,iy+1,imol,0)
c$$$                rdneureg(k)=rdneureg(k)-edissml(ix+1,iy+1,imol,0)
c$$$                recycereg(k)=recycereg(k)-edissml(ix+1,iy+1,imol,0)
                radarr(ix,iy)=radarr(ix,iy)+emolrad(ix+1,iy+1,imol,0)
                rad_neutr=rad_neutr+emolrad(ix+1,iy+1,imol,0)
                rdneureg(k)=rdneureg(k)+emolrad(ix+1,iy+1,imol,0)
                recycereg(k)=recycereg(k)+emolrad(ix+1,iy+1,imol,0)
              end do
              do iion = 1, nnioni
                radarr(ix,iy)=radarr(ix,iy)+eionrad(ix+1,iy+1,iion,0)
                rad_neutr=rad_neutr+eionrad(ix+1,iy+1,iion,0)
                rdneureg(k)=rdneureg(k)+eionrad(ix+1,iy+1,iion,0)
                recycereg(k)=recycereg(k)+eionrad(ix+1,iy+1,iion,0)
              end do
            endif
          enddo
        enddo
        write (hlp_frm,'(a,i2,a)') '(a,1p,',nnatmi,'(1x,e14.7))'
        if (dbg_eir_mc.ge.2 .and. main_call) write(*,hlp_frm)
     .   'EIRENE_MC. ENEUTRAD : ',
     .   (SUM(eneutrad(1:nx,1:ny,iatm,0)),iatm=1,nnatmi)
        if (dbg_eir_mc.ge.2 .and. main_call) write(*,hlp_frm)
     .   'EIRENE_MC. EMOLRAD : ',
     .   (SUM(emolrad(1:nx,1:ny,imol,0)),imol=1,nnmoli)
        if (dbg_eir_mc.ge.2 .and. main_call) write(*,hlp_frm)
     .   'EIRENE_MC. EIONRAD : ',
     .   (SUM(eionrad(1:nx,1:ny,iion,0)),iion=1,nnioni)
      endif
      if (dbg_eir_mc.ge.2 .and. main_call) write(*,'(a,1p,2(1x,e14.7))')
     . 'EIRENE_MC. RADARR, RAD_NEUTR : ',
     .  SUM(radarr(0:nx-1,0:ny-1)), rad_neutr

      if(.not.filedata(ntrc_int)) then
        write(*,'(a,1p,1000g12.4)')
     -      'recycecore,recycesol,recycedivin,recycedivout= ',recycereg
        write(*,'(a,1p,1000g12.4)')
     -      'recycicore,recycisol,recycidivin,recycidivout= ',recycireg
      end if

C CORRECTING THE POTENTIAL ENERGY SPENT FOR RECOMBINATION OF ATOMS AND ATOMIC IONS INTO MOLECULES
      IF(L_NEUTFLUX.GE.0) THEN
        DO I=1,NNSTSI
          k=eirdiag_nds_ind(I)
          if(k.le.0) cycle
          IS=I+nlim
          wldnep(IS,0)=0.0_R8
          DO J=1,NNMOLI
            wldnep(IS,0)=wldnep(IS,0)+ewldmr_res(j,k)
          END DO
        END DO
      END IF

cxpb  core recycling into ions of the highest ionisation stage
      if (ionising_core.gt.0) then
        rf0=RF_NEUT(1)
        if (balance_netcdf.ne.0) then
          eirene_mc_core_sna0to1=0.0_R8 !djm Jan2017
        endif
        do j = 1, nnatmi
          ! compute the average density on the core boundary (t1)
          t1 = 0.0_R8
          if (boundary_namelist.eq.1) then
            do ib = 1, nbc
              if (core_boundary(ib)) then
                do ix = 1, bc_list_size(ib)
                  t1 = t1 +
     &             na(topix(bc_list_x(ix,ib),bc_list_y(ix,ib)),
     ,                topiy(bc_list_x(ix,ib),bc_list_y(ix,ib)),
     ,                ishigh(j))*
     &             gs(topix(bc_list_x(ix,ib),bc_list_y(ix,ib)),
     ,                topiy(bc_list_x(ix,ib),bc_list_y(ix,ib)),1)
                end do
              end if
            end do
          else
            do ix = -1, nx
              if (on_closed_surface(ix,0)) then
                t1 = t1 + na(ix,0,ishigh(j))*gs(ix,0,1)
              end if
            end do
          end if
          t1 = t1/core_surface
          if (boundary_namelist.eq.1) then
            do ib = 1, nbc
              if (core_boundary(ib)) then
                do ix = 1, bc_list_size(ib)
! the recycling flux is modulated as per the poloidal density distribution
                  t0=core_flux(j)*
     *             (na(topix(bc_list_x(ix,ib),bc_list_y(ix,ib)),
     ,                 topiy(bc_list_x(ix,ib),bc_list_y(ix,ib)),
     ,                 ishigh(j))
     /              /t1)*
     *             (gs(topix(bc_list_x(ix,ib),bc_list_y(ix,ib)),
     ,                 topiy(bc_list_x(ix,ib),bc_list_y(ix,ib)),1)/
     /              core_surface)
                  if(linearisation.eq.1) then !{
                    if(t0.gt.0.0_R8) then !{
                      sna0(bc_list_x(ix,ib),
     ,                     bc_list_y(ix,ib),0,ishigh(j))=
     &                 sna0(bc_list_x(ix,ib),
     ,                      bc_list_y(ix,ib),0,ishigh(j))+t0
                      if (balance_netcdf.ne.0) then
                        eirene_mc_core_sna0to1(bc_list_x(ix,ib),
     ,                                         bc_list_y(ix,ib),0,
     ,                                         ishigh(j))
     .                  = eirene_mc_core_sna0to1(bc_list_x(ix,ib),
     ,                                           bc_list_y(ix,ib),0,
     ,                                           ishigh(j))+t0
                      endif
                    else !}{
                      sna0(bc_list_x(ix,ib),
     ,                     bc_list_y(ix,ib),1,ishigh(j))=
     .                 sna0(bc_list_x(ix,ib),
     ,                      bc_list_y(ix,ib),1,ishigh(j))
     1                +t0/na(bc_list_x(ix,ib),
     ,                       bc_list_y(ix,ib),ishigh(j))
                      if (balance_netcdf.ne.0) then
                        eirene_mc_core_sna0to1(bc_list_x(ix,ib),
     ,                                         bc_list_y(ix,ib),1,
     ,                                         ishigh(j))
     .                   = eirene_mc_core_sna0to1(bc_list_x(ix,ib),
     ,                                            bc_list_y(ix,ib),1,
     ,                                            ishigh(j))
     &                   + t0
     .                   / na(bc_list_x(ix,ib),
     ,                        bc_list_y(ix,ib),ishigh(j))
                      endif
                    end if !}
                    sna0(bc_list_x(ix,ib),bc_list_y(ix,ib),0,ishigh(j))=
     .               sna0(bc_list_x(ix,ib),bc_list_y(ix,ib),0,ishigh(j))
     1               +rf0*abs(t0)
                    if (balance_netcdf.ne.0) then
                      eirene_mc_core_sna0to1(bc_list_x(ix,ib),
     .                                       bc_list_y(ix,ib),0,
     .                                       ishigh(j))= !djm Jan2017
     &                eirene_mc_core_sna0to1(bc_list_x(ix,ib),
     ,                                       bc_list_y(ix,ib),0,
     ,                                       ishigh(j))+rf0*abs(t0)
                    endif
                    sna0(bc_list_x(ix,ib),bc_list_y(ix,ib),1,ishigh(j))=
     .               sna0(bc_list_x(ix,ib),bc_list_y(ix,ib),1,ishigh(j))
     1               -rf0*abs(t0)/
     .                na(bc_list_x(ix,ib),bc_list_y(ix,ib),ishigh(j))
                    if (balance_netcdf.ne.0) then
                      eirene_mc_core_sna0to1(bc_list_x(ix,ib),
     ,                                       bc_list_y(ix,ib),1,
     ,                                       ishigh(j))= !djm Jan2017
     &                eirene_mc_core_sna0to1(bc_list_x(ix,ib),
     ,                                       bc_list_y(ix,ib),0,
     ,                                       ishigh(j))-rf0*abs(t0)/
     .                  na(bc_list_x(ix,ib),bc_list_y(ix,ib),ishigh(j))
                    endif
                  elseif(linearisation.eq.0) then !}{
                    sna0(bc_list_x(ix,ib),bc_list_y(ix,ib),0,ishigh(j))=
     .               sna0(bc_list_x(ix,ib),bc_list_y(ix,ib),0,ishigh(j))
     .               +t0
                    if (balance_netcdf.ne.0) then
                      eirene_mc_core_sna0to1(bc_list_x(ix,ib),
     .                                       bc_list_y(ix,ib),0,
     .                                       ishigh(j))= !djm Jan2017
     &                eirene_mc_core_sna0to1(bc_list_x(ix,ib),
     .                                       bc_list_y(ix,ib),0,
     .                                       ishigh(j))+t0
                    endif
                  end if !}
                enddo  !ix}
              endif  !core}
            enddo  !ic}
          else
            do ix = -1, nx
              if (on_closed_surface(ix,0)) then
                t0=core_flux(j)*(na(ix,0,ishigh(j))/t1)
     .                         *(gs(ix,0,1)/core_surface)
                if(linearisation.eq.1) then !{
                  if(t0.gt.0.0_R8) then !{
                    sna0(ix,-1,0,ishigh(j))=sna0(ix,-1,0,ishigh(j))+t0
                    if (balance_netcdf.ne.0) then
                      eirene_mc_core_sna0to1(ix,-1,0,ishigh(j))=
     .                 eirene_mc_core_sna0to1(ix,-1,0,ishigh(j))+t0
                    endif
                  else !}{
                    sna0(ix,-1,1,ishigh(j))=sna0(ix,-1,1,ishigh(j))+
     .               t0/na(ix,-1,ishigh(j))
                    if (balance_netcdf.ne.0) then
                      eirene_mc_core_sna0to1(ix,-1,1,ishigh(j))=
     .                 eirene_mc_core_sna0to1(ix,-1,1,ishigh(j))+
     .                 t0/na(ix,-1,ishigh(j))
                    endif
                  end if !}
                  sna0(ix,-1,0,ishigh(j))=sna0(ix,-1,0,ishigh(j))+
     .             rf0*abs(t0)
                  if (balance_netcdf.ne.0) then
                    eirene_mc_core_sna0to1(ix,-1,0,ishigh(j))=
     .               eirene_mc_core_sna0to1(ix,-1,0,ishigh(j))+
     .               rf0*abs(t0)
                  endif
                  sna0(ix,-1,1,ishigh(j))=sna0(ix,-1,1,ishigh(j))-
     1             rf0*abs(t0)/na(ix,-1,ishigh(j))
                  if (balance_netcdf.ne.0) then
                    eirene_mc_core_sna0to1(ix,-1,1,ishigh(j))=
     .               eirene_mc_core_sna0to1(ix,-1,0,ishigh(j))-
     .               rf0*abs(t0)/na(ix,-1,ishigh(j))
                  endif
                elseif(linearisation.eq.0) then !}{
                  sna0(ix,-1,0,ishigh(j))=sna0(ix,-1,0,ishigh(j))+t0
                  if (balance_netcdf.ne.0) then
                    eirene_mc_core_sna0to1(ix,-1,0,ishigh(j))=
     .               eirene_mc_core_sna0to1(ix,-1,0,ishigh(j))+t0
                  endif
                end if !}
              endif  !core}
            enddo  !ix}
          endif  !boundary_namelist}
        enddo  !j}
      else if (ionising_core.lt.0) then
        ib=abs(ionising_core)
        call xertst (ib.le.nbc,
     1   'faulty input ionising_core, abs(ionising_core) > NBC')
        found = .false.
        conpar(:,ib,3) = 0.0_R8
        do j = 1, nnatmi
          is = ishigh(j)
          if(bccon(is,ib).eq.13) then
            found = .true.
            conpar(is,ib,3) = conpar(is,ib,3) + core_flux(j)
            if(output_style.ge.2 .or. ncall.lt.3)
     1       write(*,'(a,1p,1e12.5,a,i3,1x,i3,1x,i3)')
     .        'eirene_mc: conpar(is,ib,3) = ',
     2         conpar(is,ib,3),' for j,is,ib = ',j,is,ib
          end if
        end do
        if(.not.found .and. (output_style.ge.2 .or. ncall.lt.3)) then
          write(hlp_frm,'(a,i3,a)') '(a,i3,a,',ns,'(i3,1x))'
          write(*,hlp_frm)
     1     'eirene_mc: bccon for ib = ',ib,
     2     ' <> 13, current value = ',(bccon(is,ib),is=0,ns-1)
c          call xerrab
c     .     ('Negative eirene_ionising_core must be used with BCCON=13')
        endif
      endif

      cum_volrec=0.0_R8
      do istrai=1, nstrat
        cum_volrec=cum_volrec+VOLSUMN(istrai)
      enddo
      cum_volrec=-cum_volrec/qe

! { IYS 12.01.2018
      sna_eir = 0.0_R8
      if (ank_tracing .gt. 0) then
!        do is = 1, nspecies
!          do iy = -1, ny
!            do ix = -1, nx
!              sna_eir(is) = sna_eir(is) +
!     &           sna0(ix,iy,0,eb2spcr(is)+1) +
!     &           sna0(ix,iy,1,eb2spcr(is)+1) * na(ix,iy,eb2spcr(is)+1)
!            enddo
!          enddo
!        enddo
        do ix=-1, nx
          do iy=-1,ny
            do is = 0, ns-1
              if( is_neutral(is) ) cycle
              k = b2espcr(is)
              sna_eir(k) = sna_eir(k) +
     &           sna0(ix,iy,0,is) +
     &           sna0(ix,iy,1,is)*na(ix,iy,is)
            enddo
          enddo
        enddo
      endif     ! IYS 12.01.2018 }

      if((main_call.and.output_style.eq.1).or.output_style.gt.1) then !{
        do istrai=1, nstrat !{
          write (*,*) 'ISTRATUM = ',istrai
          do is=0,ns-1 !{
            write (*,'(a,i3,1p,2(1x,1e14.7))') 'IS,SSNI,SSMO = ',
     1       is,ssni(istrai,is)*qe,ssmo(istrai,is)
          end do !}
          write (*,'(a,1p,5(1x,1e14.7))') 'SSEE,SSEI = ',
     1     ssee(istrai),ssei(istrai)+sseicorr(istrai),
     2     ssei(istrai),sseicorr(istrai),sseistab(istrai)
        end do !}
        write(*,'(a,1p,7(1x,1e14.7))')
     .  'CUM:SSEE,SSEI[1234],SSNI,SSMO = ',
     1   sum(ssee(1:nstrat)),
     2   sum(ssei(1:nstrat))+sum(sseicorr(1:nstrat)),
     3   sum(ssei(1:nstrat)),sum(sseicorr(1:nstrat)),
     4   sum(sseistab(1:nstrat)),
     2   sum(ssni(1:nstrat,0:ns-1)),sum(ssmo(1:nstrat,0:ns-1))
        if (ionising_core.ne.0) then
          do j = 1, nnatmi
            write(*,*) 'Recycled',trim(ion_name(j)),
     1       ' ion flux through core boundary', core_flux(j)
          enddo
        endif
        write(*,*) 'Cumulative_volume_recombination ',cum_volrec
      end if !}

      if (main_call.or.ncall.eq.0) ncall=ncall+1

* ..return
      call subend ()
      return
      end subroutine eirene_mc

* ..finalise save
      subroutine eirene_mc_final(nx,ny,ns,BoRiS)
      implicit none
      integer, intent(in) :: nx, ny, ns
      real(kind=R8), intent(in) :: BoRiS
      if(saved_volrec_used) then
        call write_b2mod_neutrals_save(99,'    ')
        call dealloc_b2mod_neutrals_save
      end if
      call copy_background(nx, ny, ns, BoRiS)
      call write_f30(nx, ny)
      call write_f31(nx, ny)
csw 26jan2011
      call eirene_extrab25_cleanup
      call eirene_wneutrals_dealloc
#ifndef B25_EIRENE
      call dealloc_b2mod_eirdiag
#endif
csw 26feb2014
      if(allocated(sni_av)) then
       deallocate(sni_av)
       deallocate(smo_av)
       deallocate(see_av)
       deallocate(sei_av)
       if (allocated(sni_av_papl)) then !djm Jan2017
         deallocate(sne_av_pael)
         deallocate(sne_av_pmel)
         deallocate(sna_av_paat)
         deallocate(sna_av_pmat)
         deallocate(sna_av_piat)
         deallocate(sna_av_paml)
         deallocate(sna_av_pmml)
         deallocate(sna_av_piml)
         deallocate(sna_av_paio)
         deallocate(sna_av_pmio)
         deallocate(sna_av_piio)
         deallocate(sni_av_papl)
         deallocate(sni_av_pmpl)
         deallocate(sni_av_pipl)
         deallocate(sni_av_pppl)
         deallocate(smo_av_mapl)
         deallocate(smo_av_mmpl)
         deallocate(smo_av_mipl)
         deallocate(smo_av_cppv)
         deallocate(see_av_eael)
         deallocate(see_av_emel)
         deallocate(see_av_eiel)
         deallocate(see_av_epel)
         deallocate(sei_av_eapl)
         deallocate(sei_av_empl)
         deallocate(sei_av_eipl)
         deallocate(sei_av_eppl)
       endif
      endif
      if(nhist.gt.1) call dealloc_sum_arrays
      if(nrelax.gt.1 .or. nhist.gt.1) call dealloc_hist_arrays
      deallocate(feedback_strata)
      deallocate(hlp_pff,hlp_stra)

      return
      end subroutine eirene_mc_final

c======================================================================
      subroutine neutr_scl(nx, ny, ns)
      implicit none
      integer, intent(in) :: nx, ny, ns
      integer :: ix, iy, is, istrai, j, k, l
      real(kind=R8) :: s, t0, recyctmp, hlp_st
      logical hlp_log
c======================================================================
      if(dbg_eir_mc.gt.0) write(*,*) 'entry neutr_scl'
      if(ank_mods.ne.0) then
c*** Calculate the basic quantities - i.e. the neutral source strengths
c*** and total plasma sources for each stratum for subsequent rescaling

c*** First, initialise the accumulators.
c***  sclstra: neutral sources
c***  sclrtio: plasma sources

        do istrai=1, nstrat !{
          do j=0,nnatmi !{
            sclstra(j,istrai)=0.0_R8
            sclrtio(j,istrai)=0.0_R8
          end do !}
        end do !}

c*** Evaluate the new volumetric sources cell by cell
c*** Here, the calculation of t0 and s must be consistent with
c*** SNIC and SNIV evaluation above

        do iy=0,ny-1 !{
          do ix=0,nx-1 !{
            if(resignore(ix,iy,1).ne.0) then !{
              do istrai=1, nstrat !{
                do is=0,ns-1 !{
                  if(use_recyceir.gt.0) then
                    recyctmp=recyceir(istrai)
                  else
                    recyctmp=recyc(is,istrai)
                  endif
                  if(.not.is_neutral(is)) then !{
                    j=b2eatcr(is)
                    t0=sni(ix+1,iy+1,fluids_list(is),istrai)*recyctmp*
     *               tflux(istrai)*flux_scale(istrai)
c                s=t0+rf0*abs(t0)*(1.0_R8-na(ix,iy,is)/na0(ix,iy,is))
                    if(t0.ge.0.0_R8) then !{
                      s=t0
                    else !}{
                      s=t0*na(ix,iy,is)/na0(ix,iy,is)
                    end if !}
                    if(crcstra(istrai).eq.'V' .and. t0.lt.0.0_R8) then !{
                      sclstra(j,istrai)=sclstra(j,istrai)-s
c                  sclstra(0,istrai)=sclstra(0,istrai)-s
                    else !}{
                      sclrtio(j,istrai)=sclrtio(j,istrai)+
     +                 s*sclalph(j,istrai)
                      sclrtio(0,istrai)=sclrtio(0,istrai)+
     +                 s*sclalph(j,istrai)
                    end if !}
                  end if !}
                end do !}
              end do !}
            end if !}
          end do !}
        end do !}

c*** Calculate the new neutral sources for each stratum. Only the
c*** surface sources should be counted here - the volumetric
c*** recombination is already evaluated above.

        hlp_pff=.true.
        do istrai=1, nstrat !{
          do is=0,ns-1 !{
            if(.not.is_neutral(is)) then !{
              j=b2eatcr(is)
              if(use_recyceir.gt.0) then
                recyctmp=recyceir(istrai)
              else
                recyctmp=recyc(is,istrai)
              endif
              hlp_st=0.0_R8
              if(crcstra(istrai).eq.'W') then !{
                do iy=1,rc_list_size(istrai) !{
                  hlp_st=hlp_st-
     .   fna(rightix(rc_list_x(iy,istrai),rc_list_y(iy,istrai)),
     .       rightiy(rc_list_x(iy,istrai),rc_list_y(iy,istrai)),0,is)
                end do !}
              else if(crcstra(istrai).eq.'E') then !}{
                do iy=1,rc_list_size(istrai) !{
                  hlp_st=hlp_st+
     .             fna(rc_list_x(iy,istrai),rc_list_y(iy,istrai),0,is)
                end do !}
              else if(crcstra(istrai).eq.'N') then !}{
                do ix=1,rc_list_size(istrai) !{
                  hlp_st=hlp_st+
     .                 fna(rc_list_x(ix,istrai),
     .                     rc_list_y(ix,istrai),1,is)
                end do !}
              else if(crcstra(istrai).eq.'S') then !}{
                do ix=1,rc_list_size(istrai) !{
                  hlp_st=hlp_st-
     .                 fna(topix(rc_list_x(ix,istrai),
     .                           rc_list_y(ix,istrai)),
     .                     topiy(rc_list_x(ix,istrai),
     .                           rc_list_y(ix,istrai)),1,is)
                end do !}
              else if(crcstra(istrai).eq.'C') then !}{
                if(hlp_pff(j,istrai)) then !{
c     .           if(gpfc(j,istrai).gt.0.0_R8)                          !###
c     .            write(*,*) 'istrai,j,is,tflux,gpfc= ',               !###
c     .                        istrai,j,is,tflux(istrai),gpfc(j,istrai) !###
                  hlp_st=tflux(istrai)*gpfc(j,istrai)
                  hlp_pff(j,istrai)=.false.
                end if !}
c             else if(crcstra(istrai).eq.'U') then !}{
c               hlp_st=userfluxparm(istrai,1)   ! this is possibly wrong!
              end if !}
              sclstra(j,istrai)=sclstra(j,istrai)+hlp_st*recyctmp
c             sclstra(0,istrai)=sclstra(0,istrai)+hlp_st*recyctmp
            end if !}
          end do !}
        end do !}
c!###{
      if (dbg_eir_mc.gt.0) then
        print *,'@@=============== before correction factors =========='
        print '(a16,20(4x,a7,2x))',
     ,   ' Stratum Species','sclbase','sclstra',
     ,                      'sclrtio','sclalph','sclsrcn'
        do istrai=1, nstrat !{
          do j=0,nnatmi !{
            print '(i6,i8,2x,1p,20e13.5)',istrai,j,sclbase(j,istrai),
     ,         sclstra(j,istrai),sclrtio(j,istrai),sclalph(j,istrai),
     ,                                             sclsrcn(j,istrai)
          end do !}
        end do !}
        print *,'@@===================================================='
      endif
c!###}

c*** Finally, calculate the scaling factors for the plasma sources and
c*** neutral fluxes.
c*** Indeed, and equation F(sclalph)=0 should be solved for each
c*** species and stratum, where F = sclrtio - 1.

        do istrai=1, nstrat !{
          sclstra(0,istrai)=0.0_R8
          do j=1,nnatmi !{
            hlp_stra(j)=sclstra(j,istrai)
c!###{
c          hlp_pr= istrai.eq.7 .and. j.eq.3
c          if(hlp_pr) print *,'%# species,stratum,lstrascl= ',
c     ,                                      j,istrai,lstrascl(istrai,j)
c!###}
            hlp_log=.false.
            do k=1,nnatmi !{
              hlp_log = hlp_log .or. lstrascl(istrai,k).eq.j
            end do !}
            if(hlp_log) sclstra(0,istrai)=sclstra(0,istrai)+
     +       sclstra(j,istrai)
          end do !}
          do j=0,nnatmi !{
c          hlp_pr= istrai.eq.7 .and. j.eq.3 !###
            if(j.gt.0) then !{
              l=lstrascl(istrai,j)
              if(l.le.0 .or. l.gt.nnatmi) l=j
            else !}{
              l=0
            end if !}
            if(sclrtio(j,istrai).gt.0.0_R8 .and.
     &         sclbase(j,istrai).gt.0.0_R8) then !{
              sclrtio(j,istrai)=sclbase(j,istrai)*
     *         sclstra(l,istrai)/sclrtio(j,istrai)
            else !}{
              if(sclrtio(j,istrai).lt.0.0_R8)
     -         write(0,'(a,a,i3,1x,i3,1x,1p,1e14.7)')
     ,          '@eirene_mc: negative plasma source. ',
     ,          'istrai,ispecies,source = ',istrai,j,sclrtio(j,istrai)
              sclrtio(j,istrai)=1.0_R8
            end if !}
          end do !}

          do j=0,nnatmi !{
            if(j.gt.0 .and. lstrascl(istrai,j).gt.0 .and.
     .       lstrascl(istrai,j).le.nnatmi) then !{

c*** A primitive iteration scheme - to be modified later if necessary

              l=lstrascl(istrai,j)
              if(l.le.0 .or. l.gt.nnatmi) l=j
              if(sclsrcn(l,istrai).gt.0.0_R8) then !{
                sclstra(j,istrai)=hlp_stra(l)/sclsrcn(l,istrai)
              else !}{
                sclstra(j,istrai)=1.0_R8
              end if !}
              sclstra(j,istrai)=max(1.0_R8/neut_scl_lim,
     1         min(neut_scl_lim,sclstra(j,istrai)))
              if(sclbase(j,istrai).gt.0.0_R8) then !{
                sclalph(j,istrai)=sclalph(j,istrai)*sclrtio(j,istrai)
                sclalph(j,istrai)=max(1.0_R8/neut_scl_lim,
     1           min(neut_scl_lim,sclalph(j,istrai)))
              else !}{
                sclalph(j,istrai)=sclstra(j,istrai)
              end if !}
            else if(j.eq.0 .or. lstrascl(istrai,j).lt.0) then !}{
              if(sclsrcn(j,istrai).gt.0.0_R8) then !{
                sclstra(j,istrai)=sclstra(j,istrai)/sclsrcn(j,istrai)
              else !}{
                sclstra(j,istrai)=1.0_R8
              end if !}
              sclalph(j,istrai)=sclstra(j,istrai)
            else !}{
              sclstra(j,istrai)=1.0_R8
              sclalph(j,istrai)=1.0_R8
            end if !}
          end do !}
        end do !}

c!###{
      if (dbg_eir_mc.gt.0) then
        print *,'@@=============== after correction factors ==========='
        print '(a16,20(4x,a7,2x))',
     ,   ' Stratum Species','sclbase','sclstra',
     ,                      'sclrtio','sclalph','sclsrcn'
        do istrai=1, nstrat !{
          do j=0,nnatmi !{
            print '(i6,i8,2x,1p,20e13.5)',istrai,j,sclbase(j,istrai),
     ,         sclstra(j,istrai),sclrtio(j,istrai),sclalph(j,istrai),
     ,                                             sclsrcn(j,istrai)
          end do !}
        end do !}
        print *,'@@===================================================='
      endif
c!###}

      endif ! ank_mods
      return
      end subroutine neutr_scl

      subroutine alloc_hist_arrays(nhistd)      !IYS 28.03.12 {
#ifdef B25_EIRENE
      use eirmod_parmmod
      use eirmod_braeir
#else
      use b2mod_braeir
#endif
      implicit none
      integer nhistd
#ifndef B25_EIRENE
#include "COUPLE/KOPPLDIM.F"
#endif

! to store the time history of sources from EIRENE
      allocate(sni_h(0:(size(sni,1)-1), 0:(size(sni,2)-1),
     .                    size(sni,3), size(sni,4), nhistd))
      allocate(smo_h(0:(size(smo,1)-1), 0:(size(smo,2)-1),
     .                    size(smo,3), size(smo,4), nhistd))
      allocate(see_h(0:(size(see,1)-1), 0:(size(see,2)-1),
     .                    size(see,3), nhistd))
      allocate(sei_h(0:(size(sei,1)-1), 0:(size(sei,2)-1),
     .                    size(sei,3), nhistd))
      if (balance_netcdf.ne.0) then !djm Jan2017
        allocate(sna_h_paat(0:(size(sni,1)-1), 0:(size(sni,2)-1),
     .                    nnatmi, size(sni,4), nhistd))
        allocate(sna_h_pmat(0:(size(sni,1)-1), 0:(size(sni,2)-1),
     .                    nnatmi, size(sni,4), nhistd))
        allocate(sna_h_piat(0:(size(sni,1)-1), 0:(size(sni,2)-1),
     .                    nnatmi, size(sni,4), nhistd))
        allocate(sna_h_paml(0:(size(sni,1)-1), 0:(size(sni,2)-1),
     .                    nnmoli, size(sni,4), nhistd))
        allocate(sna_h_pmml(0:(size(sni,1)-1), 0:(size(sni,2)-1),
     .                    nnmoli, size(sni,4), nhistd))
        allocate(sna_h_piml(0:(size(sni,1)-1), 0:(size(sni,2)-1),
     .                    nnmoli, size(sni,4), nhistd))
        allocate(sna_h_paio(0:(size(sni,1)-1), 0:(size(sni,2)-1),
     .                    nnioni, size(sni,4), nhistd))
        allocate(sna_h_pmio(0:(size(sni,1)-1), 0:(size(sni,2)-1),
     .                    nnioni, size(sni,4), nhistd))
        allocate(sna_h_piio(0:(size(sni,1)-1), 0:(size(sni,2)-1),
     .                    nnioni, size(sni,4), nhistd))
        allocate(sni_h_papl(0:(size(sni,1)-1), 0:(size(sni,2)-1),
     .                    size(sni,3), size(sni,4), nhistd))
        allocate(sni_h_pmpl(0:(size(sni,1)-1), 0:(size(sni,2)-1),
     .                    size(sni,3), size(sni,4), nhistd))
        allocate(sni_h_pipl(0:(size(sni,1)-1), 0:(size(sni,2)-1),
     .                    size(sni,3), size(sni,4), nhistd))
        allocate(sni_h_pppl(0:(size(sni,1)-1), 0:(size(sni,2)-1),
     .                    size(sni,3), size(sni,4), nhistd))
        allocate(smo_h_mapl(0:(size(smo,1)-1), 0:(size(smo,2)-1),
     .                    size(smo,3), size(smo,4), nhistd))
        allocate(smo_h_mmpl(0:(size(smo,1)-1), 0:(size(smo,2)-1),
     .                    size(smo,3), size(smo,4), nhistd))
        allocate(smo_h_mipl(0:(size(smo,1)-1), 0:(size(smo,2)-1),
     .                    size(smo,3), size(smo,4), nhistd))
        allocate(smo_h_cppv(0:(size(smo,1)-1), 0:(size(smo,2)-1),
     .                    size(smo,3), size(smo,4), nhistd))
        allocate(sne_h_pael(0:(size(see,1)-1), 0:(size(see,2)-1),
     .                    size(see,3), nhistd))
        allocate(sne_h_pmel(0:(size(see,1)-1), 0:(size(see,2)-1),
     .                    size(see,3), nhistd))
        allocate(see_h_eael(0:(size(see,1)-1), 0:(size(see,2)-1),
     .                    size(see,3), nhistd))
        allocate(see_h_emel(0:(size(see,1)-1), 0:(size(see,2)-1),
     .                    size(see,3), nhistd))
        allocate(see_h_eiel(0:(size(see,1)-1), 0:(size(see,2)-1),
     .                    size(see,3), nhistd))
        allocate(see_h_epel(0:(size(see,1)-1), 0:(size(see,2)-1),
     .                    size(see,3), nhistd))
        allocate(sei_h_eapl(0:(size(sei,1)-1), 0:(size(sei,2)-1),
     .                    size(sei,3), nhistd))
        allocate(sei_h_empl(0:(size(sei,1)-1), 0:(size(sei,2)-1),
     .                    size(sei,3), nhistd))
        allocate(sei_h_eipl(0:(size(sei,1)-1), 0:(size(sei,2)-1),
     .                    size(sei,3), nhistd))
        allocate(sei_h_eppl(0:(size(sei,1)-1), 0:(size(sei,2)-1),
     .                    size(sei,3), nhistd))
      endif
      allocate(pfluxa_h(0:ndxp,0:ndyp,natm,1,nhistd))
      allocate(pfluxm_h(0:ndxp,0:ndyp,nmol,1,nhistd))
      allocate(rfluxa_h(0:ndxp,0:ndyp,natm,1,nhistd))
      allocate(rfluxm_h(0:ndxp,0:ndyp,nmol,1,nhistd))
      allocate(dab2_h(0:ndxp,0:ndyp,natm,1,nhistd))
      call clear_hist_arrays
      return
      end subroutine alloc_hist_arrays

      subroutine dealloc_hist_arrays
      implicit none

      if (.not.allocated(sni_h)) return
      deallocate(sni_h)
      deallocate(smo_h)
      deallocate(see_h)
      deallocate(sei_h)
      deallocate(pfluxa_h)
      deallocate(pfluxm_h)
      deallocate(rfluxa_h)
      deallocate(rfluxm_h)
      deallocate(dab2_h)
      if (.not.allocated(sni_h_papl)) return !djm Jan2017
      deallocate(sne_h_pael)
      deallocate(sne_h_pmel)
      deallocate(sna_h_paat)
      deallocate(sna_h_pmat)
      deallocate(sna_h_piat)
      deallocate(sna_h_paml)
      deallocate(sna_h_pmml)
      deallocate(sna_h_piml)
      deallocate(sna_h_paio)
      deallocate(sna_h_pmio)
      deallocate(sna_h_piio)
      deallocate(sni_h_papl)
      deallocate(sni_h_pmpl)
      deallocate(sni_h_pipl)
      deallocate(sni_h_pppl)
      deallocate(smo_h_mapl)
      deallocate(smo_h_mmpl)
      deallocate(smo_h_mipl)
      deallocate(smo_h_cppv)
      deallocate(see_h_eael)
      deallocate(see_h_emel)
      deallocate(see_h_eiel)
      deallocate(see_h_epel)
      deallocate(sei_h_eapl)
      deallocate(sei_h_empl)
      deallocate(sei_h_eipl)
      deallocate(sei_h_eppl)
      return
      end subroutine dealloc_hist_arrays

      subroutine clear_hist_arrays                                !csw 24feb2014
      implicit none
      sni_h = 0.0_R8
      smo_h = 0.0_R8
      see_h = 0.0_R8
      sei_h = 0.0_R8
      pfluxa_h = 0.0_R8
      pfluxm_h = 0.0_R8
      rfluxa_h = 0.0_R8
      rfluxm_h = 0.0_R8
      dab2_h = 0.0_R8
      if (balance_netcdf.ne.0) then !djm Jan2017
        sne_h_pael = 0.0_R8
        sne_h_pmel = 0.0_R8
        sna_h_paat = 0.0_R8
        sna_h_pmat = 0.0_R8
        sna_h_piat = 0.0_R8
        sna_h_paml = 0.0_R8
        sna_h_pmml = 0.0_R8
        sna_h_piml = 0.0_R8
        sni_h_papl = 0.0_R8
        sni_h_pmpl = 0.0_R8
        sni_h_pipl = 0.0_R8
        sni_h_pppl = 0.0_R8
        smo_h_mapl = 0.0_R8
        smo_h_mmpl = 0.0_R8
        smo_h_mipl = 0.0_R8
        smo_h_cppv = 0.0_R8
        see_h_eael = 0.0_R8
        see_h_emel = 0.0_R8
        see_h_eiel = 0.0_R8
        see_h_epel = 0.0_R8
        sei_h_eapl = 0.0_R8
        sei_h_empl = 0.0_R8
        sei_h_eipl = 0.0_R8
        sei_h_eppl = 0.0_R8
      endif
      return
      end subroutine clear_hist_arrays                            !IYS 28.03.12 }

      subroutine alloc_sum_arrays
#ifdef B25_EIRENE
      use eirmod_parmmod
      use eirmod_braeir
#else
      use b2mod_braeir
#endif
      implicit none
#ifndef B25_EIRENE
#include "COUPLE/KOPPLDIM.F"
#endif
      allocate(sni_s(0:(size(sni,1)-1), 0:(size(sni,2)-1),
     .                    size(sni,3), size(sni,4)))
      allocate(smo_s(0:(size(smo,1)-1), 0:(size(smo,2)-1),
     .                    size(smo,3), size(smo,4)))
      allocate(see_s(0:(size(see,1)-1), 0:(size(see,2)-1),
     .                    size(see,3)))
      allocate(sei_s(0:(size(sei,1)-1), 0:(size(sei,2)-1),
     .                    size(sei,3)))
      if (balance_netcdf.ne.0) then !djm Jan2017
        allocate(sna_s_paat(0:(size(sni,1)-1), 0:(size(sni,2)-1),
     .                      nnatmi, size(sni,4)))
        allocate(sna_s_pmat(0:(size(sni,1)-1), 0:(size(sni,2)-1),
     .                      nnatmi, size(sni,4)))
        allocate(sna_s_piat(0:(size(sni,1)-1), 0:(size(sni,2)-1),
     .                      nnatmi, size(sni,4)))
        allocate(sna_s_paml(0:(size(sni,1)-1), 0:(size(sni,2)-1),
     .                      nnmoli, size(sni,4)))
        allocate(sna_s_pmml(0:(size(sni,1)-1), 0:(size(sni,2)-1),
     .                      nnmoli, size(sni,4)))
        allocate(sna_s_piml(0:(size(sni,1)-1), 0:(size(sni,2)-1),
     .                      nnmoli, size(sni,4)))
        allocate(sna_s_paio(0:(size(sni,1)-1), 0:(size(sni,2)-1),
     .                      nnioni, size(sni,4)))
        allocate(sna_s_pmio(0:(size(sni,1)-1), 0:(size(sni,2)-1),
     .                      nnioni, size(sni,4)))
        allocate(sna_s_piio(0:(size(sni,1)-1), 0:(size(sni,2)-1),
     .                      nnioni, size(sni,4)))
        allocate(sni_s_papl(0:(size(sni,1)-1), 0:(size(sni,2)-1),
     .                      size(sni,3), size(sni,4)))
        allocate(sni_s_pmpl(0:(size(sni,1)-1), 0:(size(sni,2)-1),
     .                      size(sni,3), size(sni,4)))
        allocate(sni_s_pipl(0:(size(sni,1)-1), 0:(size(sni,2)-1),
     .                      size(sni,3), size(sni,4)))
        allocate(sni_s_pppl(0:(size(sni,1)-1), 0:(size(sni,2)-1),
     .                      size(sni,3), size(sni,4)))
        allocate(smo_s_mapl(0:(size(smo,1)-1), 0:(size(smo,2)-1),
     .                      size(smo,3), size(smo,4)))
        allocate(smo_s_mmpl(0:(size(smo,1)-1), 0:(size(smo,2)-1),
     .                      size(smo,3), size(smo,4)))
        allocate(smo_s_mipl(0:(size(smo,1)-1), 0:(size(smo,2)-1),
     .                      size(smo,3), size(smo,4)))
        allocate(smo_s_cppv(0:(size(smo,1)-1), 0:(size(smo,2)-1),
     .                      size(smo,3), size(smo,4)))
        allocate(sne_s_pael(0:(size(see,1)-1), 0:(size(see,2)-1),
     .                      size(see,3)))
        allocate(sne_s_pmel(0:(size(see,1)-1), 0:(size(see,2)-1),
     .                      size(see,3)))
        allocate(see_s_eael(0:(size(see,1)-1), 0:(size(see,2)-1),
     .                      size(see,3)))
        allocate(see_s_emel(0:(size(see,1)-1), 0:(size(see,2)-1),
     .                      size(see,3)))
        allocate(see_s_eiel(0:(size(see,1)-1), 0:(size(see,2)-1),
     .                      size(see,3)))
        allocate(see_s_epel(0:(size(see,1)-1), 0:(size(see,2)-1),
     .                      size(see,3)))
        allocate(sei_s_eapl(0:(size(sei,1)-1), 0:(size(sei,2)-1),
     .                      size(sei,3)))
        allocate(sei_s_empl(0:(size(sei,1)-1), 0:(size(sei,2)-1),
     .                      size(sei,3)))
        allocate(sei_s_eipl(0:(size(sei,1)-1), 0:(size(sei,2)-1),
     .                      size(sei,3)))
        allocate(sei_s_eppl(0:(size(sei,1)-1), 0:(size(sei,2)-1),
     .                      size(sei,3)))
      endif
      allocate(pfluxa_s(0:ndxp,0:ndyp,natm,1))
      allocate(pfluxm_s(0:ndxp,0:ndyp,nmol,1))
      allocate(rfluxa_s(0:ndxp,0:ndyp,natm,1))
      allocate(rfluxm_s(0:ndxp,0:ndyp,nmol,1))
      allocate(dab2_s(0:ndxp,0:ndyp,natm,1))
      call clear_sum_arrays
      return
      end subroutine alloc_sum_arrays

      subroutine dealloc_sum_arrays
      implicit none
      deallocate(sni_s)
      deallocate(smo_s)
      deallocate(see_s)
      deallocate(sei_s)
      if (allocated(sni_s_papl)) then !djm Jan2017
        deallocate(sne_s_pael)
        deallocate(sne_s_pmel)
        deallocate(sna_s_paat)
        deallocate(sna_s_pmat)
        deallocate(sna_s_piat)
        deallocate(sna_s_paml)
        deallocate(sna_s_pmml)
        deallocate(sna_s_piml)
        deallocate(sna_s_paio)
        deallocate(sna_s_pmio)
        deallocate(sna_s_piio)
        deallocate(sni_s_papl)
        deallocate(sni_s_pmpl)
        deallocate(sni_s_pipl)
        deallocate(sni_s_pppl)
        deallocate(smo_s_mapl)
        deallocate(smo_s_mmpl)
        deallocate(smo_s_mipl)
        deallocate(smo_s_cppv)
        deallocate(see_s_eael)
        deallocate(see_s_emel)
        deallocate(see_s_eiel)
        deallocate(see_s_epel)
        deallocate(sei_s_eapl)
        deallocate(sei_s_empl)
        deallocate(sei_s_eipl)
        deallocate(sei_s_eppl)
      endif
      deallocate(pfluxa_s)
      deallocate(pfluxm_s)
      deallocate(rfluxa_s)
      deallocate(rfluxm_s)
      deallocate(dab2_s)
      return
      end subroutine dealloc_sum_arrays

      subroutine clear_sum_arrays                                !csw 24feb2014
      implicit none
      sni_s = 0.0_R8
      smo_s = 0.0_R8
      see_s = 0.0_R8
      sei_s = 0.0_R8
      if (balance_netcdf.ne.0) then !djm Jan2017
        sne_s_pael = 0.0_R8
        sne_s_pmel = 0.0_R8
        sni_s_papl = 0.0_R8
        sni_s_pmpl = 0.0_R8
        sni_s_pipl = 0.0_R8
        sni_s_pppl = 0.0_R8
        smo_s_mapl = 0.0_R8
        smo_s_mmpl = 0.0_R8
        smo_s_mipl = 0.0_R8
        smo_s_cppv = 0.0_R8
        see_s_eael = 0.0_R8
        see_s_emel = 0.0_R8
        see_s_eiel = 0.0_R8
        see_s_epel = 0.0_R8
        sei_s_eapl = 0.0_R8
        sei_s_empl = 0.0_R8
        sei_s_eipl = 0.0_R8
        sei_s_eppl = 0.0_R8
      endif
      pfluxa_s = 0.0_R8
      pfluxm_s = 0.0_R8
      rfluxa_s = 0.0_R8
      rfluxm_s = 0.0_R8
      dab2_s = 0.0_R8
      return
      end subroutine clear_sum_arrays                            !IYS 28.03.12 }

#ifdef B25_EIRENE
      subroutine strata_check
      implicit none
      integer itrg, jtrg, ntrg, ierror, icnt, is1, is2
      character*10 hlp_frm
      logical first

      ntrg = size(rcpos_eir)
      write(hlp_frm,'(a,i3,a)') '(a,',ntrg,'i5)'
      ierror = 0
      jtrg = 0
      icnt = 0
      first = .true.
      do itrg = 1, ntrg
        if (rcprt_eir(itrg).eq.1) then
          icnt = 1
          first = .true.
          jtrg = jtrg + 1
          if (is_neutral(species_start(jtrg))) then
            is1 = fluids_list(species_start(jtrg)+1)
          else
            is1 = fluids_list(species_start(jtrg))
          end if
          is2 = fluids_list(species_end(jtrg))
          if (crcstra(jtrg).ne.rcchr_eir(itrg) .or.
     &        rcpos(jtrg)  .ne.rcpos_eir(itrg) .or.
     &        rcstart(jtrg).ne.rcbeg_eir(itrg) .or.
     &        rcend(jtrg)  .ne.rcend_eir(itrg) .or.
     &        is1          .ne.rcspi_eir(itrg) .or.
     &        is2          .ne.rcspe_eir(itrg)) ierror = ierror + 1
        else if (rcprt_eir(itrg).gt.1) then
          if (first) icnt = 1
          if (first) jtrg = jtrg + 1
          if (is_neutral(species_start(jtrg))) then
            is1 = fluids_list(species_start(jtrg)+1)
          else
            is1 = fluids_list(species_start(jtrg))
          end if
          is2 = fluids_list(species_end(jtrg))
          if (crcstra(jtrg).ne.rcchr_eir(itrg) .or.
     &        rcpos(jtrg)  .ne.rcpos_eir(itrg) .or.
     &        rcstart(jtrg).gt.rcbeg_eir(itrg) .or.
     &        rcend(jtrg)  .lt.rcend_eir(itrg) .or.
     &        is1          .ne.rcspi_eir(itrg) .or.
     &        is2          .ne.rcspe_eir(itrg)) ierror = ierror + 1
          if (.not.first) icnt = icnt + 1
          first = icnt.eq.rcprt_eir(itrg)
        end if
      end do
      if (ierror.gt.0) write(*,*)
     & 'Mismatch between Eirene and B2 recycling strata !'
      if (trcint.or.ierror.gt.0) then
        write(*,*)
        write(*,*) 'From B2:'
        write(*,*) 'RCCHR = ',crcstra(1:jtrg)
        write(*,hlp_frm) 'RCPOS = ',rcpos(1:jtrg)
        write(*,hlp_frm) 'RCSTRT= ',rcstart(1:jtrg)
        write(*,hlp_frm) 'RCEND = ',rcend(1:jtrg)
        write(*,hlp_frm) 'SPEC. = ',species_start(1:jtrg)
        write(*,hlp_frm) ' TO   = ',species_end(1:jtrg)
        write(*,*)
        write(*,*) 'From EIRENE:'
        write(*,*) 'RCCHR = ',rcchr_eir(1:ntrg)
        write(*,hlp_frm) 'RCPOS = ',rcpos_eir(1:ntrg)
        write(*,hlp_frm) 'RCSTRT= ',rcbeg_eir(1:ntrg)
        write(*,hlp_frm) 'RCEND = ',rcend_eir(1:ntrg)
        write(*,hlp_frm) 'NSPZI = ',rcspi_eir(1:ntrg)
        write(*,hlp_frm) 'NSPZE = ',rcspe_eir(1:ntrg)
      end if
      call xertst(ierror.eq.0,
     & 'Mismatch between Eirene and B2 recycling strata !')

      return
      end subroutine strata_check
#endif

      end module b2mod_eirene
