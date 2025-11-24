      module b2mod_usrtrc
      use b2mod_types
      use b2mod_diag_diffv
      use b2mod_time
      use b2mod_b2cmpa_diffv
      use b2mod_user_namelist_diffv
      use b2mod_b2plot_wall_loading
#ifdef B25_EIRENE
      use b2mod_b2plot
     .   , only : b2_fnmti
      use eirmod_parmmod
     .   , only : nstra, natm, nmol
      use eirmod_wneutrals
      use eirmod_eirbra
      use eirmod_clgin
     ,   , only : iliin
#else
      use b2mod_eirbra
#endif
      use b2mod_dimensions
      use b2mod_geometry_diffv
      implicit none

c=======================================================================
c*** Puts user-specific quantities on the output file
c*** "tracing/user.trc"
c-----------------------------------------------------------------------
      integer lun, ltrc
      parameter (ltrc=ntrc_usr)

c-----------------------------------------------------------------------
c*** lv    :  version number to appear in the file
c*** lvmj  :  version number referring to the last major change

      integer lv,lvmj
      parameter (lv=20000727, lvmj=970212)

c----------------------------------------------------------------------
c*** debug_flags usage:
c***  xx  >0 : debugging printout on
c***  81     : target fluxes from b2usr_loads printed to targb2*.dat
c***  99     : enables debug printouts in b2mod_usrtrc.F
c-----------------------------------------------------------------------
c*** User-specific data:
c***
c***     header     var       Legend
c***                                                          (global)
c***    Pneu_PFR  hlp_ppfr   neutral pressure in PFR - mean
c***    Pneu_PFi  hlp_ppfi   neutral pressure in PFR - inner
c***    Pneu_PFo  hlp_ppfo   neutral pressure in PFR - outer
c***    He_fract  hlp_hefr   helium fraction in the pump
c***    He_enrch  hlp_henr   helium "enrichment"
c***    He_p_prs  hlp_hepr   helium partial pressure in PFR
c***    pk_pwr_i  hlp_ppwi   peak power load onto the inner target
c***    pk_pwr_o  hlp_ppwo   peak power load onto the outer target
c***    pk_pw_iu  hlp_pwiu   peak power load on the inner upper target
c***    pk_pw_ou  hlp_pwou   peak power load on the outer upper target
c***                         (the peaks are calculated in a more simple
c***                         way compared to user_spb however they allow
c**                          for an additional multiplicator to the plasma loads)
c***    Sp_DT     hlp_spdt   effective pumping speed for DT
c***    Sp_He     hlp_sphe   effective pumping speed for He
c***    Sp_He_sc  hlp_sshe   re-scaled pumping speed for He
c***    He_src    hlp_hesr   helium source
c***    He_sc_cr  hlp_hesc   rescaled average nHe/ne at the CEI
c***    He_sc_sp  hlp_hess   rescaled average nHe/ne at the separatrix
c***    He_og_cr  hlp_heoc   original average nHe/ne at the CEI
c***    He_og_sp  hlp_heos   original average nHe/ne at the separatrix
c***    He_enrcs  hlp_hens   helium enrichment to ns from throughputs
c***    He_enrcc  hlp_henc   helium enrichment to nc from throughputs
c***    DT_thrpt  hlp_dttp   DT throughput
c***    He_thrpt  hlp_hetp   He throughput
c***    pk_pwp_w  hlp_ppfw   peak parallel heat flux to main chamber wall
c***                                       (PFR - by group and species)
c***    PFRflx    hlp_pfrf   net nuclei flux through each surface group
c***                              (neutral data - by group and species)
c***    n_        hlp_dndts  neutral density                      [1/m**3]
c***    P_        hlp_prdts  neutral pressure                     [Pa]
c***    T_        hlp_tmdts  neutral temperature                  [eV]
c***    j_        hlp_fxdts  impinging neutral flux               [1/m**2/s]
c***    G_        hlp_fldts  neutral flux pumped out              [1/s]
c***    S_        hlp_stdts  technical pumping speed (albedo)     [m**3/s]
c***    s_        hlp_sedts  effective pumping speed (G/n)        [m**3/s]
c***                              (PFR interface - by side and species)
c***    PFRi_area hlp_w      PFR interface area, inner            [m**2]
c***    PFRi_n_   hlp_dnpfs  average neutral density, inner       [1/m**3]
c***    PFRi_T_   hlp_tmpfs  average neutral temperature, inner   [eV]
c***    PFRi_v_   hlp_vlpfs  average neutral flow velocity, inner [m/s]
c***    PFRo_area hlp_w      PFR interface area, outer            [m**2]
c***    PFRo_n_   hlp_dnpfs  average neutral density, outer       [1/m**3]
c***    PFRo_T_   hlp_tmpfs  average neutral temperature, outer   [eV]
c***    PFRo_v_   hlp_vlpfs  average neutral flow velocity, outer [m/s]
c***    DoD_il    hlp_dodil  min. DoD, inner lower divertor
c***    DoD_ol    hlp_dodol  min. DoD, outer lower divertor
c***    DoD_iu    hlp_dodiu  min. DoD, inner upper divertor
c***    DoD_ou    hlp_dodou  min. DoD, outer upper divertor
c***    isatt_il  hlp_istil  total Isat to inner lower divertor
c***    isatt_ol  hlp_istol  total Isat to outer lower divertor
c***    isatt_iu  hlp_istiu  total Isat to inner upper divertor
c***    isatt_ou  hlp_istou  total Isat to inner upper divertor
c***    albll_mx  al_bll_mx  ballooning alpha parameter maximum value
c***    albll_ic  al_bll_ic  ballooning alpha parameter maximum iCv
c-----------------------------------------------------------------------

c*** Global values
      integer nglb_tmp, ng_tmp
      parameter (nglb_tmp=23)

c*** File markers: routine version, number of global values, number of
c***               PFR values, number of species, number of PFR surface
c***               groups, number of diagnostic surface groups
      integer nrkn_hlp
      parameter(nrkn_hlp=7)
      integer mrkn_hlp(nrkn_hlp), mrks_hlp(nrkn_hlp)
c***    number of tallies, particle index (<0 - mol.), species index
      integer ndta_tmp(nntrgrpd,DEF_NATM+DEF_NMOL),
     *        ndtm_tmp(DEF_NATM+DEF_NMOL), ndtx_tmp(DEF_NATM+DEF_NMOL)
c-----------------------------------------------------------------------
c*** Number of positions in the output line and temporary arrays for
c*** the data and header lines
      integer ntrc_tmp
      real (kind=R8), allocatable :: hlp_out(:)
      character*(mhdr), allocatable :: hlp_hdr(:)
c-----------------------------------------------------------------------

c*** Neutral-related data on surfaces
c***                          - by diagnostic surface groups and species
c***    max. and min. number of tallies, number of species
      integer ndtk_tmp, nmnt_tmp, ndts_tmp
      parameter (ndtk_tmp=7, nmnt_tmp=4)
c*** Neutral-related data on PFR interfaces
c***                                       - by side and neutral species
c***    number of PFR interfaces and of tallies
      integer ngl2_tmp
      parameter (ngl2_tmp=4)
      integer ngl3_tmp
      parameter (ngl3_tmp=4)
      integer ngl4_tmp
      parameter (ngl4_tmp=2)
      integer npff_tmp, npfi_tmp, npft_tmp
      parameter (npff_tmp=1, npfi_tmp=2, npft_tmp=3)
      character, save :: hlp_hg(nglb_tmp)*8
      character, save :: hlp_hg2(ngl2_tmp)*8
      character, save :: hlp_hg3(ngl3_tmp)*10
      character, save :: hlp_hg4(ngl4_tmp)*10
      character, save :: hlp_hn(ndtk_tmp)
      character, save :: hlp_hp(npff_tmp)*8
      character, save :: hlp_hr(npft_tmp), hlp_hs(npfi_tmp)
c-----------------------------------------------------------------------
      integer, save :: itg_tmp(4)
      real (kind=R8), save :: hlp_ksdts(nntrgrpd,DEF_NATM+DEF_NMOL)
c*** Column headers:
      data hlp_hg /
     1      'Pneu_PFR', 'Pneu_PFi', 'Pneu_PFo', 'He_fract', 'He_enrch',
     2      'He_p_prs', 'pk_pwr_i', 'pk_pwr_o', 'pk_pw_iu', 'pk_pw_ou',
     3      'Sp_DT   ', 'Sp_He   ', 'Sp_He_sc', 'He_src  ', 'He_sc_cr',
     4      'He_sc_sp', 'He_og_cr', 'He_og_sp', 'He_enrcs', 'He_enrcc',
     5      'DT_thrpt', 'He_thrpt', 'pk_ppf_w' /
      data hlp_hg2 /
     1      'DoD_il  ', 'DoD_ol  ', 'DoD_iu  ', 'DoD_ou  ' /
      data hlp_hg3 /
     1      'isatt_il', 'isatt_ol', 'isatt_iu', 'isatt_ou' /
      data hlp_hg4 /
     1      'albll_mx', 'albll_ic' /
      data hlp_hp / 'PFR_flow' /
      data hlp_hn / 'n','p','T','j','G','S','s' /
      data hlp_hr / 'n', 'T', 'v' /
      data hlp_hs / 'i', 'o' /

      integer, private :: gridGeometry, plasmaGeometry

      integer, save :: n_pad_str
      real (kind=R8), save :: xdum(1)
      character, save :: tab, frm*16, hlp_frm*120

      external fix_user

      contains

      subroutine b2usrtrc(mpg,geo,ns,st,st_ext,switch)
      use b2us_geo_diffv
      use b2us_plasma_diffv
      use b2mod_constants
      use b2mod_boundary_namelist_diffv
      use b2mod_numerics_namelist_diffv
#ifdef B25_EIRENE
      use eirmod_extrab25
#endif

      implicit none
      type(mapping), intent(in) :: mpg
      type(geometry), intent(in) :: geo
      type(B2state), intent(in) :: st
      type(B2StateExt), intent(in) :: st_ext
      type(switches), intent(in) :: switch
      integer, intent(in) :: ns

      real (kind=R8) ::
     *     equg_tmp(nglb_tmp)
      real (kind=R8) ::
     *     hlp_ppfr, hlp_ppfi, hlp_ppfo, hlp_hefr, hlp_henr, hlp_hepr,
     ,     hlp_ppwi, hlp_ppwo, hlp_pwiu, hlp_pwou, hlp_spdt, hlp_sphe,
     ,     hlp_sshe, hlp_hesr, hlp_hesc, hlp_hess, hlp_heoc, hlp_heos,
     ,     hlp_hens, hlp_henc, hlp_dttp, hlp_hetp, hlp_ppfw
      equivalence
     ,     (equg_tmp(1),  hlp_ppfr), (equg_tmp(2),  hlp_ppfi),
     ,     (equg_tmp(3),  hlp_ppfo), (equg_tmp(4),  hlp_hefr),
     ,     (equg_tmp(5),  hlp_henr), (equg_tmp(6),  hlp_hepr),
     ,     (equg_tmp(7),  hlp_ppwi), (equg_tmp(8),  hlp_ppwo),
     ,     (equg_tmp(9),  hlp_pwiu), (equg_tmp(10), hlp_pwou),
     ,     (equg_tmp(11), hlp_spdt), (equg_tmp(12), hlp_sphe),
     ,     (equg_tmp(13), hlp_sshe), (equg_tmp(14), hlp_hesr),
     ,     (equg_tmp(15), hlp_hesc), (equg_tmp(16), hlp_hess),
     ,     (equg_tmp(17), hlp_heoc), (equg_tmp(18), hlp_heos),
     ,     (equg_tmp(19), hlp_hens), (equg_tmp(20), hlp_henc),
     ,     (equg_tmp(21), hlp_dttp), (equg_tmp(22), hlp_hetp),
     ,     (equg_tmp(23), hlp_ppfw)
      real (kind=R8) ::
     *     equg2_tmp(ngl2_tmp)
      real (kind=R8) ::
     *     hlp_dod(4), hlp_dodil, hlp_dodol, hlp_dodiu, hlp_dodou
      equivalence (equg2_tmp(1), hlp_dod),
     ,     (equg2_tmp(1), hlp_dodil), (equg2_tmp(2), hlp_dodol),
     ,     (equg2_tmp(3), hlp_dodiu), (equg2_tmp(4), hlp_dodou)
      integer mgl3_tmp
      real (kind=R8) ::
     *     equg3_tmp(ngl3_tmp)
      real (kind=R8) ::
     *     hlp_isat(4), hlp_istil, hlp_istol, hlp_istiu, hlp_istou
      equivalence (equg3_tmp(1), hlp_isat),
     ,     (equg3_tmp(1), hlp_istil), (equg3_tmp(2), hlp_istol),
     ,     (equg3_tmp(3), hlp_istiu), (equg3_tmp(4), hlp_istou)
c     integer kgl4_tmp(2)
      real (kind=R8) ::
     *     equg4_tmp(ngl4_tmp)
      real (kind=R8) ::
     *     hlp_albll(ngl4_tmp), al_bll_mx, al_bll_ic
      equivalence (equg4_tmp(1), hlp_albll),
     ,     (equg4_tmp(1), al_bll_mx), (equg4_tmp(2), al_bll_ic)
      real (kind=R8) :: wtarg(DEF_NYD,2*DEF_NCUT),
     ,     wnutrg(DEF_NYD,2*DEF_NCUT),wnutpr(DEF_NYD,2*DEF_NCUT)
      real(kind=R8) :: wrad(DEF_NYD,2*DEF_NCUT,0:DEF_NATM)

c*** PFR flow values
      real (kind=R8) ::
     *     equp_tmp(0:npfrgrpd,DEF_NATM,npff_tmp)
      real (kind=R8) ::
     *     hlp_pfrf(0:npfrgrpd,DEF_NATM,npff_tmp)
      equivalence (equp_tmp(0,1,1), hlp_pfrf)
      character hlp_hb*8

c***  atomic/molecular weight, group area, absorption
      real (kind=R8) ::
     *  hlp_awdts(DEF_NATM+DEF_NMOL),hlp_ardts(nntrgrpd)
      real (kind=R8) ::
     *   equn_tmp(nntrgrpd,DEF_NATM+DEF_NMOL,ndtk_tmp),
     ,  hlp_dndts(nntrgrpd,DEF_NATM+DEF_NMOL),
     ,  hlp_prdts(nntrgrpd,DEF_NATM+DEF_NMOL),
     ,  hlp_tmdts(nntrgrpd,DEF_NATM+DEF_NMOL),
     ,  hlp_fxdts(nntrgrpd,DEF_NATM+DEF_NMOL),
     ,  hlp_fldts(nntrgrpd,DEF_NATM+DEF_NMOL),
     ,  hlp_stdts(nntrgrpd,DEF_NATM+DEF_NMOL),
     ,  hlp_sedts(nntrgrpd,DEF_NATM+DEF_NMOL)
      equivalence
     ,  (equn_tmp(1,1,1),hlp_dndts), (equn_tmp(1,1,2),hlp_prdts),
     ,  (equn_tmp(1,1,3),hlp_tmdts), (equn_tmp(1,1,4),hlp_fxdts),
     ,  (equn_tmp(1,1,5),hlp_fldts), (equn_tmp(1,1,6),hlp_stdts),
     ,  (equn_tmp(1,1,7),hlp_sedts)
      character hlp_m*16

c***  surface area
      real (kind=R8) ::
     *  hlp_w(npfi_tmp)
      real (kind=R8) ::
     *  equr_tmp(npfi_tmp,DEF_NATM+DEF_NMOL,npft_tmp),
     ,  hlp_dnpfs(npfi_tmp,DEF_NATM+DEF_NMOL),
     ,  hlp_tmpfs(npfi_tmp,DEF_NATM+DEF_NMOL),
     ,  hlp_vlpfs(npfi_tmp,DEF_NATM+DEF_NMOL)
      equivalence
     ,  (equr_tmp(1,1,1),hlp_dnpfs), (equr_tmp(1,1,2),hlp_tmpfs),
     ,  (equr_tmp(1,1,3),hlp_vlpfs)

      real (kind=R8) ::
     *  hlp, hlp_pi, hlp_si, hlp_po, hlp_so,
     ,  hlp_he, hlp_hep,
     ,  hlp_ne_he_c, hlp_ne_he_s,
     ,  hlp_alcff
      real (kind=R8) ::
     *  d, e, s, x, kef, taf, uaf
      real (kind=R8) :: pz(1:mpg%nCv), pv(1:mpg%nVx), grpz(1:mpg%nCv),
     *  albll(1:mpg%nCi)
      integer i, j, k, l, m, n, i0, is, iCv, iFc, iFt, irt,
     ,  j_he_i, l_pad_str
#ifdef B25_EIRENE
      integer ia, ias, ic, im, iyp, js, ks, ls
      real (kind=R8) :: a, g, t, v, hlp_h, hlp_p, hlp_s, hlp_tot, hlpnw
c     real (kind=R8) :: e0
      logical is_target
#endif
c*** radial position of the separatrix and poloidal positions of
c*** targets and midplane for DoD calculation for each target
      integer use_eirene
      logical ex, first, first_out, hlp_dbg
      character f*(lenfnm)
      data use_eirene / 0 /
      data first,first_out / 2*.true. /
      real (kind=R8), external :: intface_s
c-----------------------------------------------------------------------

      save
      if(first_time_step) return
      if(.not.filedata(ltrc)) return
c=======================================================================

      hlp_dbg=debug_flags(99).gt.0

      if (first) then !{

        if(hlp_dbg) then !{
          print *
          print *,'b2usrtrc: before initialisation. ',
     ,     'natm,nmol,nnatmi,nnmoli,nntrgrp=',
     ,      natm,nmol,nnatmi,nnmoli,nntrgrp
          if(nntrgrp.gt.0) then !{
            print '(1x,a9,2x,20i4/(12x,20i4))','lntrgrp:',
     ,                                (lntrgrp(i),i=1,jntrgrp(nntrgrp))
            print '(1x,a9,2x,16a9)','gntrgrp:',
     ,                                (gntrgrp(i),i=1,nntrgrp)
            print '(1x,a9,1x,i4,15i9)','intrgrp:',
     ,                                (intrgrp(i),i=1,nntrgrp)
            print '(1x,a9,1x,i4,15i9)','jntrgrp:',
     ,                                (jntrgrp(i),i=1,nntrgrp)
            print *
          end if !}
        endif !}

        tab=char(9)
        call ipgeti ('b2mndr_eirene', use_eirene)
        call alloc_b2plot_wall_loading(nlim,nsgmx)
c-----------------------------------------------------------------------
c*** Check the topology

c*** Check the input data for consistency and prepare indices

        j_he_i=0
        do j=1,nspecies !{
          if(nfluids(j).eq.1) then !{
            do i=1,3 !{
              if(nint(am(eb2spcr(j))).eq.i) then !{
                if(j_h_at(i).le.0 .or. j_h_at(i).gt.nnatmi) then !{
                  j_h_at(i)=b2eatcr(eb2spcr(j))
                end if !}
              end if !}
            end do !}
          end if !}
          if(nfluids(j).eq.2) then !{
            j_he_i=eb2spcr(j)
            if(j_he_at.le.0 .or. j_he_at.gt.nnatmi) then !{
              j_he_at=b2eatcr(eb2spcr(j))
            end if !}
          end if !}
        end do !}

        if(hlp_dbg) then !{
          print '(a,5(i3,1x),2a)',
     ,     '%% b2usrtrc: j_he_at,j_h_at(1:3),j_he_i,'//
     ,                  'gridGeometry,plasmaGeometry = ', !###
     ,                   j_he_at,j_h_at     ,j_he_i,
     ,                   geometryName(gridGeometry),
     ,                   geometryName(plasmaGeometry)     !###
        endif !}

c*** Prepare for neutral data evaluation at diagnostic surfaces

c*** Check whether the specified surfaces are real

c*** Calculate the basic quantities for each surface/species

        ndts_tmp=nnatmi+nnmoli
        do i=1,ndts_tmp !{
          if(i.le.nnatmi) then !{
            ndtm_tmp(i)=i
            hlp_awdts(i)=am(eb2atcr(i))
          else !}{
            ndtm_tmp(i)=nnatmi-i
            hlp_awdts(i)=0.0_R8
            do j=1,nnatmi !{
              hlp_awdts(i)=hlp_awdts(i)+
     +                                mlcmp(j,-ndtm_tmp(i))*hlp_awdts(j)
            end do !}
          end if !}
        end do !}

        do j=1,ndts_tmp !{
          if(ndtm_tmp(j).gt.0) then !{
            m=ndtm_tmp(j)
          else !}{
            m=natm-ndtm_tmp(j)
          end if !}
          ndtx_tmp(j)=m
caap tbd add nphoti shift
        end do !}

        hlp_ardts=0.0_R8
        ndta_tmp=nmnt_tmp
        if (use_eirene.ne.0) then
          do k=1,nntrgrp !{
            do i=intrgrp(k),jntrgrp(k) !{
              l=lntrgrp(i)
              hlp_ardts(k)=hlp_ardts(k)+wlarea(l)
              do j=1,ndts_tmp !{
                m=ndtx_tmp(j)
                x=wlabsrp(m,lntrgrp(intrgrp(k)))
                if(i.eq.intrgrp(k)) then !{
                  if(x.gt.0.0_R8) ndta_tmp(k,j)=ndtk_tmp
                  hlp_ksdts(k,j)=x
                else !}{
                  if(x.ge.0.0_R8 .and. wlabsrp(m,l).ne.x) then !{
                    ndta_tmp(k,j)=nmnt_tmp
                    print '(4a,i3)',
     ,                 'b2usrtrc: different absorption coefficients',
     ,                 ' for diagnostic surface group ',gntrgrp(k),
     ,                 ', species ',j
                    hlp_ksdts(k,j)=0.0_R8
                    exit
                  end if !}
                end if !}
              end do !}
            end do !}
          end do !}

          if(hlp_dbg) then !{
           print *
           print *,'b2usrtrc: after initialisation. ',
     ,             'nntrgrp,ndtk_tmp,nmnt_tmp,ndts_tmp=',
     ,              nntrgrp,ndtk_tmp,nmnt_tmp,ndts_tmp
           write(f,'(a,i2,a)') '(1x,a9,2x,',ndts_tmp,'i3)'
           print *,'tallies (ndta_tmp):'
           print f,(gntrgrp(i),(ndta_tmp(i,j),j=1,ndts_tmp),i=1,nntrgrp)
           print *
           print f,'ndtm_tmp',(ndtm_tmp(j),j=1,ndts_tmp)
           print f,'ndtx_tmp',(ndtx_tmp(j),j=1,ndts_tmp)
           write(f,'(a,i2,a)') '(1x,a9,2x,1p,',ndts_tmp,'e10.2)'
           print f,'hlp_awdts',(hlp_awdts(j),j=1,ndts_tmp)
           print *,'area/absorption:'
           write(f,'(a,i2,a)') '(1x,a9,2x,1p,',ndts_tmp+1,'e10.2)'
           print f,(gntrgrp(i),hlp_ardts(i),
     ,                        (hlp_ksdts(i,j),j=1,ndts_tmp),i=1,nntrgrp)
           print *
          endif !}

        end if

c*** Find the B2 cell ranges corresponding to PFR interfaces
c***  (ipfc_tmp, jpfc_tmp)

        hlp_w=0
        if(lpfrs_pmp.ne.0.and.npfr_cvs.gt.0) then !{
          do i=1,npfr_cvs !{
            k=mod(mpg%cvReg(pfr_cvs(i)),4)
            if (k.eq.3) hlp_w(1)=hlp_w(1)+geo%cvSz(pfr_cvs(i))
            if (k.eq.0) hlp_w(2)=hlp_w(2)+geo%cvSz(pfr_cvs(i))
          end do !}
        end if !}

        if(hlp_dbg) then !{
          print *
          print *,'b2usrtrc: PFR initialisation. npfi=',npfi_tmp
          do k=1,npfi_tmp !{
            print *,'PFR areaB2 ',k,hlp_w(k)
          end do !}
          print *
        endif !}

c-----------------------------------------------------------------------
c*** Prepare constants for ballooning alpha calculation

        if(R0_albll.le.0.0_R8.or.q95_albll.le.0.0_R8) then
          e=1.0e20_R8
          d=0.0_R8
          s=0.0_R8
          do iFc = 1, mpg%nFc
            if (mpg%fcReg(iFc).eq.0) cycle
            select case ( gridGeometry )
            case ( GEOMETRY_LINEAR, GEOMETRY_STRUCTURED_SPACES )
              if (mpg%fcReg(iFc).eq.3) then
                e=min(e,geo%vxX(mpg%fcVx(iFc,1)),
     .                  geo%vxX(mpg%fcVx(iFc,2)))
                d=max(d,geo%vxX(mpg%fcVx(iFc,1)),
     .                  geo%vxX(mpg%fcVx(iFc,2)))
                s=s+geo%fcHt(iFc)*geo%fcBb(iFc,0)/geo%fcBb(iFc,3)
              end if
            case ( GEOMETRY_CYLINDER, GEOMETRY_ANNULUS )
              if (mpg%fcReg(iFc).eq.2) then
                e=min(e,geo%vxX(mpg%fcVx(iFc,1)),
     .                  geo%vxX(mpg%fcVx(iFc,2)))
                d=max(d,geo%vxX(mpg%fcVx(iFc,1)),
     .                  geo%vxX(mpg%fcVx(iFc,2)))
                s=s+geo%fcHt(iFc)*geo%fcBb(iFc,0)/geo%fcBb(iFc,3)
              end if
            case ( GEOMETRY_LIMITER )
              if (mpg%fcReg(iFc).eq.4) then
                e=min(e,geo%vxX(mpg%fcVx(iFc,1)),
     .                  geo%vxX(mpg%fcVx(iFc,2)))
                d=max(d,geo%vxX(mpg%fcVx(iFc,1)),
     .                  geo%vxX(mpg%fcVx(iFc,2)))
                s=s+geo%fcHt(iFc)*geo%fcBb(iFc,0)/geo%fcBb(iFc,3)
              end if
            case ( GEOMETRY_SN )
              if (mpg%fcReg(iFc).eq.8) then
                e=min(e,geo%vxX(mpg%fcVx(iFc,1)),
     .                  geo%vxX(mpg%fcVx(iFc,2)))
                d=max(d,geo%vxX(mpg%fcVx(iFc,1)),
     .                  geo%vxX(mpg%fcVx(iFc,2)))
                s=s+geo%fcHt(iFc)*geo%fcBb(iFc,0)/geo%fcBb(iFc,3)
              end if
            case ( GEOMETRY_STELLARATORISLAND )
              if (mpg%fcReg(iFc).eq.9) then
                e=min(e,geo%vxX(mpg%fcVx(iFc,1)),
     .                  geo%vxX(mpg%fcVx(iFc,2)))
                d=max(d,geo%vxX(mpg%fcVx(iFc,1)),
     .                  geo%vxX(mpg%fcVx(iFc,2)))
                s=s+geo%fcHt(iFc)*geo%fcBb(iFc,0)/geo%fcBb(iFc,3)
              end if
            case ( GEOMETRY_CDN )
              if (mpg%fcReg(iFc).eq.14 .or. mpg%fcReg(iFc).eq.21) then
                e=min(e,geo%vxX(mpg%fcVx(iFc,1)),
     .                  geo%vxX(mpg%fcVx(iFc,2)))
                d=max(d,geo%vxX(mpg%fcVx(iFc,1)),
     .                  geo%vxX(mpg%fcVx(iFc,2)))
                s=s+geo%fcHt(iFc)*geo%fcBb(iFc,0)/geo%fcBb(iFc,3)
              end if
            case ( GEOMETRY_DDN_BOTTOM, GEOMETRY_DDN_TOP )
              if (mpg%fcReg(iFc).eq.15 .or. mpg%fcReg(iFc).eq.22) then
                e=min(e,geo%vxX(mpg%fcVx(iFc,1)),
     .                  geo%vxX(mpg%fcVx(iFc,2)))
                d=max(d,geo%vxX(mpg%fcVx(iFc,1)),
     .                  geo%vxX(mpg%fcVx(iFc,2)))
                s=s+geo%fcHt(iFc)*geo%fcBb(iFc,0)/geo%fcBb(iFc,3)
              end if
            case ( GEOMETRY_LFS_SNOWFLAKE_MINUS,
     &             GEOMETRY_LFS_SNOWFLAKE_PLUS )
              if (mpg%fcReg(iFc).eq.15) then
                e=min(e,geo%vxX(mpg%fcVx(iFc,1)),
     .                  geo%vxX(mpg%fcVx(iFc,2)))
                d=max(d,geo%vxX(mpg%fcVx(iFc,1)),
     .                  geo%vxX(mpg%fcVx(iFc,2)))
                s=s+geo%fcHt(iFc)*geo%fcBb(iFc,0)/geo%fcBb(iFc,3)
              end if
            end select
          end do
          if(R0_albll.le.0.0_R8) R0_albll=0.5_R8*(e+d)
          if(q95_albll.le.0.0_R8) q95_albll=s/(2.0_R8*pi*R0_albll)
        endif
        if(B0_albll.le.0.0_R8) B0_albll=
     &   geo%cvBb(omp(icsepomp),3)*geo%cvX(omp(icsepomp))/R0_albll
        if(min(q95_albll,B0_albll,R0_albll).gt.0.0_R8) then
          hlp_alcff=(q95_albll/B0_albll)**2*R0_albll*2.0_R8*mu0
        else
          hlp_alcff=0.0_R8
        endif
        if(hlp_dbg) print '(1x,a,4e12.4)', "### albll: B,R,q,hlp =>",
     .                            B0_albll,R0_albll,q95_albll,hlp_alcff

c-----------------------------------------------------------------------
c*** Allocate the buffer arrays for header and table lines

        ntrc_tmp=1+nglb_tmp+nntrgrp+
     +    npfi_tmp*(1+npft_tmp*ndts_tmp)+ngl2_tmp+ngl4_tmp
        k=1+npfrgrp
        if(npfrgrp.le.1) k=max(0,npfrgrp)
        ntrc_tmp=ntrc_tmp+nspecies*npff_tmp*k
        do k=1,nntrgrp !{
          do j=1,ndts_tmp !{
            ntrc_tmp=ntrc_tmp+ndta_tmp(k,j)
            if(hlp_dbg) then !{
              print '(a,4i6)','ntrc_tmp,k,j,ndta_tmp(k,j):',
     ,          ntrc_tmp,k,j,ndta_tmp(k,j)
            endif !}
          end do !}
        end do !}
        if(print_isat) then !{
          if(plasmaGeometry.ne.GEOMETRY_ANNULUS) ntrc_tmp=ntrc_tmp+2
          if(plasmaGeometry.eq.GEOMETRY_CDN .or.
     &       plasmaGeometry.eq.GEOMETRY_DDN_BOTTOM .or.
     &       plasmaGeometry.eq.GEOMETRY_DDN_TOP .or.
     &       plasmaGeometry.eq.GEOMETRY_LFS_SNOWFLAKE_MINUS .or.
     &       plasmaGeometry.eq.GEOMETRY_LFS_SNOWFLAKE_PLUS)
     >      ntrc_tmp=ntrc_tmp+2
        end if !}

        if(hlp_dbg) then !{
          write(6,'(3x,11a10)') 'ntrc_tmp','nglb_tmp','nspecies',
     ,     'npff_tmp','npfrgrp','nntrgrp','npfi_tmp','npft_tmp',
     ,     'ndts_tmp','ngl2_tmp','ngl4_tmp'
          write(6,'(11i10)') ntrc_tmp,nglb_tmp,nspecies,
     ,      npff_tmp,npfrgrp,nntrgrp,npfi_tmp,npft_tmp,
     ,      ndts_tmp,ngl2_tmp,ngl4_tmp
          endif !}

        if(.not.allocated(hlp_out)) then !{
          allocate(hlp_out(ntrc_tmp), hlp_hdr(ntrc_tmp))
        end if !}

        first=.false.
      end if !}
c=======================================================================
c*** Initialise the arrays

      equg_tmp=0.0_R8
      equp_tmp=0.0_R8
      equn_tmp=0.0_R8
      equr_tmp=0.0_R8

#ifdef B25_EIRENE
c*** Re-scale the neutral fluxes

      do i=1,nnlimi+nnstsi !{
        do j=1,nnatmi !{
          wldna(i,j,0)=0.0_R8
          wldra(i,j,0)=0.0_R8
          wldpa(i,j,0)=0.0_R8
          do k=1,nstrai !{
            wldna(i,j,0)=wldna(i,j,0)+wldna(i,j,k)*sclstra(j,k)
            wldra(i,j,0)=wldra(i,j,0)+wldra(i,j,k)*sclstra(j,k)
            wldpa(i,j,0)=wldpa(i,j,0)+wldpa(i,j,k)*sclstra(j,k)
          end do !}
        end do !}
        do j=1,nnmoli !{
          l=lmolscl(j)
          wldnm(i,j,0)=0.0_R8
          wldrm(i,j,0)=0.0_R8
          wldpm(i,j,0)=0.0_R8
          do k=1,nstrai !{
            if(l.gt.0 .and. l.le.nnatmi) then !{
              wldnm(i,j,0)=wldnm(i,j,0)+wldnm(i,j,k)*sclstra(l,k)
              wldrm(i,j,0)=wldrm(i,j,0)+wldrm(i,j,k)*sclstra(l,k)
              wldpm(i,j,0)=wldpm(i,j,0)+wldpm(i,j,k)*sclstra(l,k)
            else !}{
              wldnm(i,j,0)=wldnm(i,j,0)+wldnm(i,j,k)
              wldrm(i,j,0)=wldrm(i,j,0)+wldrm(i,j,k)
              wldpm(i,j,0)=wldpm(i,j,0)+wldpm(i,j,k)
            end if !}
          end do !}
        end do !}
      end do !}
#endif
c=======================================================================
      if(j_he_at.gt.0) then !{

c*** Calculate the helium-related data. First, average densities
c*** for the CEI...

        hlp_pi=0.0_R8
        hlp_he=0.0_R8
        hlp_si=0.0_R8
        hlp_ne_he_c=0.0_R8

c*** ... and separatrix
        hlp_po=0.0_R8
        hlp_hep=0.0_R8
        hlp_so=0.0_R8
        hlp_ne_he_s=0.0_R8

        if(hlp_dbg) print *,'%% b2usrtrc: j_he_i = ',j_he_i !###

        do iFc=1,mpg%nFc !{
          if (mpg%fcReg(iFc).eq.0) cycle
          select case ( gridGeometry )
          case ( GEOMETRY_LINEAR, GEOMETRY_STRUCTURED_SPACES )
            if (mpg%fcReg(iFc).eq.3) then
              iCv=min(mpg%fcCv(iFc,1),mpg%fcCv(iFc,2))
              hlp_si=hlp_si+geo%fcS(iFc)
              hlp_pi=hlp_pi+geo%fcS(iFc)*st%dv%ne(iCv)
              hlp=0.0_R8
              do j=0,2 !{
                hlp=hlp+st%pl%na(iCv,j_he_i+j)
                hlp_ne_he_c=
     .           hlp_ne_he_c+st%pl%na(iCv,j_he_i+j)*j*geo%fcS(iFc)
              end do !}
              hlp_he=hlp_he+geo%fcS(iFc)*hlp
            end if
            if (mpg%fcFs(iFc).eq.mpg%ifssep.and.mpg%ifssep.ne.0) then
              iCv=min(mpg%fcCv(iFc,1),mpg%fcCv(iFc,2))
              hlp_so=hlp_so+geo%fcS(iFc)
              hlp_po=hlp_po+geo%fcS(iFc)*st%dv%ne(iCv)
              do j=0,2 !{
                hlp=hlp+st%pl%na(iCv,j_he_i+j)
                hlp_ne_he_s=hlp_ne_he_s+
     .           st%pl%na(iCv,j_he_i+j)*j*geo%fcS(iFc)
              end do !}
              hlp_hep=hlp_hep+geo%fcS(iFc)*hlp
            end if
          case ( GEOMETRY_CYLINDER, GEOMETRY_ANNULUS )
            if (mpg%fcReg(iFc).eq.2) then
              iCv=min(mpg%fcCv(iFc,1),mpg%fcCv(iFc,2))
              hlp_si=hlp_si+geo%fcS(iFc)
              hlp_pi=hlp_pi+geo%fcS(iFc)*st%dv%ne(iCv)
              hlp=0.0_R8
              do j=0,2 !{
                hlp=hlp+st%pl%na(iCv,j_he_i+j)
                hlp_ne_he_c=
     .           hlp_ne_he_c+st%pl%na(iCv,j_he_i+j)*j*geo%fcS(iFc)
              end do !}
              hlp_he=hlp_he+geo%fcS(iFc)*hlp
            end if
            if (mpg%fcFs(iFc).eq.mpg%ifssep.and.mpg%ifssep.ne.0) then
              iCv=min(mpg%fcCv(iFc,1),mpg%fcCv(iFc,2))
              hlp_so=hlp_so+geo%fcS(iFc)
              hlp_po=hlp_po+geo%fcS(iFc)*st%dv%ne(iCv)
              do j=0,2 !{
                hlp=hlp+st%pl%na(iCv,j_he_i+j)
                hlp_ne_he_s=hlp_ne_he_s+
     .           st%pl%na(iCv,j_he_i+j)*j*geo%fcS(iFc)
              end do !}
              hlp_hep=hlp_hep+geo%fcS(iFc)*hlp
            end if
          case ( GEOMETRY_LIMITER )
            if (mpg%fcReg(iFc).eq.4) then
              iCv=min(mpg%fcCv(iFc,1),mpg%fcCv(iFc,2))
              hlp_si=hlp_si+geo%fcS(iFc)
              hlp_pi=hlp_pi+geo%fcS(iFc)*st%dv%ne(iCv)
              hlp=0.0_R8
              do j=0,2 !{
                hlp=hlp+st%pl%na(iCv,j_he_i+j)
                hlp_ne_he_c=
     .           hlp_ne_he_c+st%pl%na(iCv,j_he_i+j)*j*geo%fcS(iFc)
              end do !}
              hlp_he=hlp_he+geo%fcS(iFc)*hlp
            end if
            if (mpg%fcReg(iFc).eq.5) then
              if (mpg%cvOnClosedSurface(mpg%fcCv(iFc,1))) then
                iCv=mpg%fcCv(iFc,1)
              else
                iCv=mpg%fcCv(iFc,2)
              end if
              hlp_so=hlp_so+geo%fcS(iFc)
              hlp_po=hlp_po+geo%fcS(iFc)*st%dv%ne(iCv)
              do j=0,2 !{
                hlp=hlp+st%pl%na(iCv,j_he_i+j)
                hlp_ne_he_s=hlp_ne_he_s+
     .           st%pl%na(iCv,j_he_i+j)*j*geo%fcS(iFc)
              end do !}
              hlp_hep=hlp_hep+geo%fcS(iFc)*hlp
            end if
          case ( GEOMETRY_SN )
            if (mpg%fcReg(iFc).eq.8) then
              iCv=min(mpg%fcCv(iFc,1),mpg%fcCv(iFc,2))
              hlp_si=hlp_si+geo%fcS(iFc)
              hlp_pi=hlp_pi+geo%fcS(iFc)*st%dv%ne(iCv)
              hlp=0.0_R8
              do j=0,2 !{
                hlp=hlp+st%pl%na(iCv,j_he_i+j)
                hlp_ne_he_c=
     .           hlp_ne_he_c+st%pl%na(iCv,j_he_i+j)*j*geo%fcS(iFc)
              end do !}
              hlp_he=hlp_he+geo%fcS(iFc)*hlp
            end if
            if (mpg%fcReg(iFc).eq.10) then
              if (mpg%cvOnClosedSurface(mpg%fcCv(iFc,1))) then
                iCv=mpg%fcCv(iFc,1)
              else
                iCv=mpg%fcCv(iFc,2)
              end if
              hlp_so=hlp_so+geo%fcS(iFc)
              hlp_po=hlp_po+geo%fcS(iFc)*st%dv%ne(iCv)
              do j=0,2 !{
                hlp=hlp+st%pl%na(iCv,j_he_i+j)
                hlp_ne_he_s=hlp_ne_he_s+
     .           st%pl%na(iCv,j_he_i+j)*j*geo%fcS(iFc)
              end do !}
              hlp_hep=hlp_hep+geo%fcS(iFc)*hlp
            end if
          case ( GEOMETRY_CDN )
            if (mpg%fcReg(iFc).eq.14 .or. mpg%fcReg(iFc).eq.21) then
              iCv=min(mpg%fcCv(iFc,1),mpg%fcCv(iFc,2))
              hlp_si=hlp_si+geo%fcS(iFc)
              hlp_pi=hlp_pi+geo%fcS(iFc)*st%dv%ne(iCv)
              hlp=0.0_R8
              do j=0,2 !{
                hlp=hlp+st%pl%na(iCv,j_he_i+j)
                hlp_ne_he_c=
     .           hlp_ne_he_c+st%pl%na(iCv,j_he_i+j)*j*geo%fcS(iFc)
              end do !}
              hlp_he=hlp_he+geo%fcS(iFc)*hlp
            end if
            if (mpg%fcReg(iFc).eq.16 .or. mpg%fcReg(iFc).eq.23) then
              if (mpg%cvOnClosedSurface(mpg%fcCv(iFc,1))) then
                iCv=mpg%fcCv(iFc,1)
              else
                iCv=mpg%fcCv(iFc,2)
              end if
              hlp_so=hlp_so+geo%fcS(iFc)
              hlp_po=hlp_po+geo%fcS(iFc)*st%dv%ne(iCv)
              do j=0,2 !{
                hlp=hlp+st%pl%na(iCv,j_he_i+j)
                hlp_ne_he_s=hlp_ne_he_s+
     .           st%pl%na(iCv,j_he_i+j)*j*geo%fcS(iFc)
              end do !}
              hlp_hep=hlp_hep+geo%fcS(iFc)*hlp
            end if
          case ( GEOMETRY_DDN_TOP, GEOMETRY_DDN_BOTTOM )
            if (mpg%fcReg(iFc).eq.14 .or. mpg%fcReg(iFc).eq.21) then
              iCv=min(mpg%fcCv(iFc,1),mpg%fcCv(iFc,2))
              hlp_si=hlp_si+geo%fcS(iFc)
              hlp_pi=hlp_pi+geo%fcS(iFc)*st%dv%ne(iCv)
              hlp=0.0_R8
              do j=0,2 !{
                hlp=hlp+st%pl%na(iCv,j_he_i+j)
                hlp_ne_he_c=
     .           hlp_ne_he_c+st%pl%na(iCv,j_he_i+j)*j*geo%fcS(iFc)
              end do !}
              hlp_he=hlp_he+geo%fcS(iFc)*hlp
            end if
            if (mpg%fcReg(iFc).eq.16 .or. mpg%fcReg(iFc).eq.23) then
              if (mpg%cvOnClosedSurface(mpg%fcCv(iFc,1))) then
                iCv=mpg%fcCv(iFc,1)
              else
                iCv=mpg%fcCv(iFc,2)
              end if
              hlp_so=hlp_so+geo%fcS(iFc)
              hlp_po=hlp_po+geo%fcS(iFc)*st%dv%ne(iCv)
              do j=0,2 !{
                hlp=hlp+st%pl%na(iCv,j_he_i+j)
                hlp_ne_he_s=hlp_ne_he_s+
     .           st%pl%na(iCv,j_he_i+j)*j*geo%fcS(iFc)
              end do !}
              hlp_hep=hlp_hep+geo%fcS(iFc)*hlp
            end if
          case ( GEOMETRY_STELLARATORISLAND )
            if (mpg%fcReg(iFc).eq.9) then
              iCv=min(mpg%fcCv(iFc,1),mpg%fcCv(iFc,2))
              hlp_si=hlp_si+geo%fcS(iFc)
              hlp_pi=hlp_pi+geo%fcS(iFc)*st%dv%ne(iCv)
              hlp=0.0_R8
              do j=0,2 !{
                hlp=hlp+st%pl%na(iCv,j_he_i+j)
                hlp_ne_he_c=
     .           hlp_ne_he_c+st%pl%na(iCv,j_he_i+j)*j*geo%fcS(iFc)
              end do !}
              hlp_he=hlp_he+geo%fcS(iFc)*hlp
            end if
            if (mpg%fcReg(iFc).eq.11) then
              if (mpg%cvOnClosedSurface(mpg%fcCv(iFc,1))) then
                iCv=mpg%fcCv(iFc,1)
              else
                iCv=mpg%fcCv(iFc,2)
              end if
              hlp_so=hlp_so+geo%fcS(iFc)
              hlp_po=hlp_po+geo%fcS(iFc)*st%dv%ne(iCv)
              do j=0,2 !{
                hlp=hlp+st%pl%na(iCv,j_he_i+j)
                hlp_ne_he_s=hlp_ne_he_s+
     .           st%pl%na(iCv,j_he_i+j)*j*geo%fcS(iFc)
              end do !}
              hlp_hep=hlp_hep+geo%fcS(iFc)*hlp
            end if
          case ( GEOMETRY_LFS_SNOWFLAKE_MINUS,
     &           GEOMETRY_LFS_SNOWFLAKE_PLUS )
            if (mpg%fcReg(iFc).eq.15) then
              iCv=min(mpg%fcCv(iFc,1),mpg%fcCv(iFc,2))
              hlp_si=hlp_si+geo%fcS(iFc)
              hlp_pi=hlp_pi+geo%fcS(iFc)*st%dv%ne(iCv)
              hlp=0.0_R8
              do j=0,2 !{
                hlp=hlp+st%pl%na(iCv,j_he_i+j)
                hlp_ne_he_c=
     .           hlp_ne_he_c+st%pl%na(iCv,j_he_i+j)*j*geo%fcS(iFc)
              end do !}
              hlp_he=hlp_he+geo%fcS(iFc)*hlp
            end if
            if (mpg%fcReg(iFc).eq.17) then
              if (mpg%cvOnClosedSurface(mpg%fcCv(iFc,1))) then
                iCv=mpg%fcCv(iFc,1)
              else
                iCv=mpg%fcCv(iFc,2)
              end if
              hlp_so=hlp_so+geo%fcS(iFc)
              hlp_po=hlp_po+geo%fcS(iFc)*st%dv%ne(iCv)
              do j=0,2 !{
                hlp=hlp+st%pl%na(iCv,j_he_i+j)
                hlp_ne_he_s=hlp_ne_he_s+
     .           st%pl%na(iCv,j_he_i+j)*j*geo%fcS(iFc)
              end do !}
              hlp_hep=hlp_hep+geo%fcS(iFc)*hlp
            end if
          end select
        end do !}
        if(hlp_pi.gt.0.0_R8) then !{
          hlp_heoc=hlp_he/hlp_pi
          hlp_ne_he_c=hlp_ne_he_c/hlp_pi
        else !}{
          hlp_heoc=0.0_R8
          hlp_ne_he_c=0.0_R8
        end if !}
        if(hlp_po.gt.0.0_R8) then !{
          hlp_heos=hlp_hep/hlp_po
          hlp_ne_he_s=hlp_ne_he_s/hlp_po
        else !}{
          hlp_heos=0.0_R8
          hlp_ne_he_s=0.0_R8
        end if !}
      end if !}
c=======================================================================
#ifdef B25_EIRENE
      if(lpfrs_pmp.gt.0 .and. lpfrs_pmp.le.4) then !{

c-----------------------------------------------------------------------
c*** Calculate the pumping-related data.
c***
c*** First, estimate the neutral pressure in the duct - as the average
c*** over the plasma cells adjacent to the duct, weighted with the
c*** surface area. The pump location is specified by parameters
c***
c***  lpfrs_pmp:  1 - in PFR, 2 - on the outer wall,
c***              3 - in upper PFR, 4 - on the inner wall
c***  pfr_cvs  :  list of control volumes used for the neutral averaging
c***
c***  j_he_at  is the helium index in the Eirene atoms list. If 0, then
c***           the B2 species is taken assuming the same order.
c***  j_h_at   is the same for the hydrogen isotopes, and
c***  l_h_mol  gives the number of hydrogen isotope nuclei in each molecule
c***  lhetrgts is the list of Eirene surfaces used to estimate the
c***           He enrichment from the single-sided fluxes in Eirene.
c-----------------------------------------------------------------------

        hlp_p=0.0_R8
        hlp_s=0.0_R8
        hlp_hep=0.0_R8

        if(hlp_dbg) then !{
          print '(1x,a,6i4)',
     ,          '%%## j_he_at,nnatmi',!###
     ,                j_he_at,nnatmi  !###
        endif !}

        hlp_pi=0.0_R8
        hlp_si=0.0_R8
        hlp_po=0.0_R8
        hlp_so=0.0_R8
        if(lpfrs_pmp.ne.0) then !{

c*** Western divertor

          do j=1,npfr_cvs !{
            if(mod(mpg%cvReg(pfr_cvs(j)),4).ne.0) cycle
            hlp_si=hlp_si+geo%cvSz(pfr_cvs(j))
            hlp=0.0_R8
            do k=1,nnatmi !{
              hlp=hlp+dab2(pfr_cvs(j),k,1)*tab2(pfr_cvs(j),k,1)
              if(k.eq.j_he_at) hlp_hep=hlp_hep+
     +          dab2(pfr_cvs(j),k,1)*tab2(pfr_cvs(j),k,1)*
     .          geo%cvSz(pfr_cvs(j))
            end do !}
            do k=1,nnmoli !{
              hlp=hlp+dmb2(pfr_cvs(j),k,1)*tmb2(pfr_cvs(j),k,1)
            end do !}
            hlp_pi=hlp_pi+hlp*geo%cvSz(pfr_cvs(j))
          end do !}

          if(hlp_dbg) then !{
            print '(a,1p,2i4,2e12.3)',' %%## m,n,hlp_si,hlp_pi: ',  !###
     ,                                       m,n,hlp_si,hlp_pi      !###
          endif !}

c*** Eastern divertor

          do j=1,npfr_cvs !{
            if(mod(mpg%cvReg(pfr_cvs(j)),4).ne.3) cycle
            hlp_so=hlp_so+geo%cvSz(pfr_cvs(j))
            hlp=0.0_R8
            do k=1,nnatmi !{
              hlp=hlp+dab2(pfr_cvs(j),k,1)*tab2(pfr_cvs(j),k,1)
              if(k.eq.j_he_at) hlp_hep=hlp_hep+
     +          dab2(pfr_cvs(j),k,1)*tab2(pfr_cvs(j),k,1)*
     .          geo%cvSz(pfr_cvs(j))
            end do !}
            do k=1,nnmoli !{
              hlp=hlp+dmb2(pfr_cvs(j),k,1)*tmb2(pfr_cvs(j),k,1)
            end do !}
            hlp_po=hlp_po+hlp*geo%cvSz(pfr_cvs(j))
          end do !}

          if(hlp_dbg) then !{
            print '(a,1p,2i4,2e12.3)',' %%## m,n,hlp_so,hlp_po: ',  !###
     ,                                       m,n,hlp_so,hlp_po      !###
          endif !}

          hlp_s=hlp_si+hlp_so
          hlp_p=hlp_pi+hlp_po
        end if !}

        if(hlp_dbg)
     .    print '(a,1p,2i4,3e12.3)',' %%## m,n,hlp_s,hlp_hep,hlp_p: ',
     ,                                     m,n,hlp_s,hlp_hep,hlp_p

        if(hlp_s.gt.0.0_R8) then !{
          hlp_ppfr=hlp_p/hlp_s
          hlp_hepr=hlp_hep/hlp_s
        end if !}
        if(hlp_si.gt.0.0_R8) hlp_ppfi=hlp_pi/hlp_si
        if(hlp_so.gt.0.0_R8) hlp_ppfo=hlp_po/hlp_so

c*** Calculate throughputs, pumping speed, and re-scale He if necessary

        if(hlp_dbg) print *,' %%## l,js,hlp_wa,hlp_wm,hlpnw: '
        do l=1,4 !{
          if (l.eq.4) then !{
            ls=j_he_at
          else !}{
            ls=j_h_at(l)
          end if !}
          hlpnw=0.0_R8
          if(ls.gt.0 .and. ls.le.nnatmi) then !{

c*** Neutral flux across the side walls

            do js=1,nnlimi !{
              if(iliin(js).le.0) cycle
              hlpnw=hlpnw+wldna(js,ls,0)-wldra(js,ls,0)
              if(hlp_dbg) hlp_tot=0.0_R8
              do j=1,nnmoli !{
                if(mlcmp(ls,j).gt.0) hlpnw=hlpnw+
     +           mlcmp(ls,j)*(wldnm(js,j,0)-wldrm(js,j,0))
                if(mlcmp(ls,j).gt.0.and.hlp_dbg)
     .           hlp_tot = hlp_tot +
     +                     mlcmp(ls,j)*(wldnm(js,j,0)-wldrm(js,j,0))
              end do !}
              if(hlp_dbg) print '(1x,2i4,3e12.3)',
     .                       l,js,wldna(js,ls,0)-wldra(js,ls,0),
     ,                       hlp_tot,hlpnw
            end do !}

c*** flux across other recycling surfaces for vessel mode grids
c*** we must remove targets and transparent surfaces from that list
            do js=1,nsts !{
              is_target = .false.
              do ic=1,ntns !{
                is_target = is_target.or.js.eq.ltns(ic)
              end do !}
              if (.not.is_target) then !{
                ks = js+nnlimi
                if (iliin(js).gt.0) then !{
                  hlpnw=hlpnw+wldna(ks,ls,0)-wldra(ks,ls,0)
                  do j=1,nnmoli !{
                    if(mlcmp(ls,j).gt.0) hlpnw=hlpnw+
     .                mlcmp(ls,j)*(wldnm(ks,j,0)-wldrm(ks,j,0))
                    if(mlcmp(ls,j).gt.0.and.hlp_dbg)
     .               hlp_tot = hlp_tot +
     +                mlcmp(ls,j)*(wldnm(ks,j,0)-wldrm(ks,j,0))
                  end do !}
                  if(hlp_dbg) print '(1x,2i4,3e12.3)',
     .                       l,js,wldna(ks,ls,0)-wldra(ks,ls,0),
     ,                       hlp_tot,hlpnw
                end if !}
              end if !}
            end do !}
          end if !}
          if (l.le.3) then !{
            hlp_dttp=hlp_dttp+hlpnw/2.7e20_R8
          else if (l.eq.4) then !}{
            hlp_hetp=hlpnw/2.7e20_R8
          end if !}
        end do !}
        hlp_dttp=0.5_R8*hlp_dttp     ! assume molecular DT pumping
        if(hlp_ppfr.gt.0.0_R8) hlp_spdt=hlp_dttp/hlp_ppfr
        if(hlp_hepr.gt.0.0_R8) hlp_sphe=hlp_hetp/hlp_hepr
        if(scale_he) then !{
          hlp_sshe=spmp_he_to_d*hlp_spdt
        else !}{
          hlp_sshe=hlp_spdt
        end if !}

c*** Calculate the helium-related data

        if(hlp_dbg) then !{
          print '(a,1p,2i4,1e12.3)',
     ,           '%%## j_he_at,nhetrgtx,fusion_power =',
     ,                 j_he_at,nhetrgtx,fusion_power      !###
        endif !}

c*** Re-scale the helium density upstream to meet both throughput and
c*** pumping speed

        if(fusion_power.gt.0.0_R8) then !{
          hlp_hesr=3.55e17_R8/2.7e20_R8*fusion_power
          hlp=hlp_sshe*hlp_hepr
          if(hlp.gt.0.0_R8) then !{
            hlp=hlp_hesr/hlp
            if(hlp_ne_he_c.lt.1.0_R8) then !{
              hlp_hesc=
     .         1.0_R8-1.0_R8/(1.0_R8+hlp_heoc*hlp/(1.0_R8-hlp_ne_he_c))
            else !}{
              hlp_hesc=1
            end if !}
            if(hlp_ne_he_s.lt.1.0_R8) then !{
              hlp_hess=
     .         1.0_R8-1.0_R8/(1.0_R8+hlp_heos*hlp/(1.0_R8-hlp_ne_he_s))
            else !}{
              hlp_hess=1.0_R8
            end if !}
          end if !}
          hlp=2.0_R8*hlp_dttp*hlp_hesc
          if(hlp.gt.0.0_R8) hlp_henc=hlp_hesr/hlp
          hlp=2.0_R8*hlp_dttp*hlp_hess
          if(hlp.gt.0.0_R8) hlp_hens=hlp_hesr/hlp
        end if !}
      end if !}

      if(j_he_at.gt.0) then !{

c*** Helium fraction and related data - flux estimation

        hlp_tot=0.0_R8
        hlp_he=0.0_R8
        hlp_h=0.0_R8
        do i=1,nhetrgtx !{
          k=lhetrgts(i)
          if(k.eq.0) go to 10
          if(k.lt.0) k=nnlimi-k
          if(j_he_at.le.nnatmi) hlp_he=hlp_he+wldna(k,j_he_at,0)
          do j=1,3 !{
            if(j_h_at(j).le.nnatmi.and.j_h_at(j).gt.0)
     .       hlp_h=hlp_h+wldna(k,j_h_at(j),0)
          end do !}
          do j=1,nnatmi !{
            hlp_tot=hlp_tot+wldna(k,j,0)
          end do !}
          do j=1,nnmoli !{
            hlp_tot=hlp_tot+wldnm(k,j,0)
            hlp_h=hlp_h+wldnm(k,j,0)*
     .       (l_h_mol(j,1)+l_h_mol(j,2)+l_h_mol(j,3))
          end do !}
        end do !}
   10   if(hlp_h.gt.0.0_R8) hlp_hefr=hlp_he/(hlp_h+hlp_he)
        if(hlp_heoc.gt.0.0_R8) hlp_henr=hlp_hefr/hlp_heoc
      end if !}
#endif

c*** Calculate the peak powers to the targets

      wtarg=0.0_R8
      wnutpr=0.0_R8
      wnutrg=0.0_R8
      call b2ptrdl(DEF_NYD,2*DEF_NCUT,wrad,.true.)
c West target
      hlp_ppwi=0.0_R8
      do i=1,mpg%divFcP(1,2) !{
        iFc=mpg%divFc(mpg%divFcP(1,1)+i-1)
        iCv=min(mpg%fcCv(iFc,1),mpg%fcCv(iFc,2))
        s=geo%fcS(iFc)
        hlp=abs(st%dv%fht(iFc,0)+st%dv%fht(iFc,1)+
     .          st_ext%fhi(iFc,0)+st_ext%fhi(iFc,1))
        do is=0,st_ext%ns-1 !{
          taf=intface_s(iFc,mpg%nCv,mpg%nFc,mpg%fcCv,geo%fcVol,
     .     st_ext%ta(1,is))
          uaf=intface_s(iFc,mpg%nCv,mpg%nFc,mpg%fcCv,geo%fcVol,
     .     st_ext%ua(1,is))
          hlp=hlp+abs((st_ext%fa(iFc,0,is)+st_ext%fa(iFc,1,is))*(
     .     st_ext%pt(iCv,is)*ev+(1.0_R8-switch%BoRiS)*
     .    (0.5_R8*st_ext%am(is)*mp*uaf**2+taf)))
        end do !}
        hlp=hlp/s
#ifdef B25_EIRENE
        ias=b2_fnmti(iFc)
        hlp=hlp+ewldt_res(ias)
        do is=1,nnatmi
          wnutrg(i,1)=wnutrg(i,1)+ewlda_res(is,ias)-ewldea_res(is,ias)
        end do
        do is=1,nnmoli
          wnutrg(i,1)=wnutrg(i,1)+ewldm_res(is,ias)-ewldem_res(is,ias)+
     .                            ewldmr_res(is,ias)
        end do
        wnutpr(i,1)=ewldrp_res(ias)
#endif
        wtarg(i,1)=(hlp-wnutpr(i,1))*trgshp(1)
     .                 +wnutrg(i,1)+wrad(i,1,0)
        hlp_ppwi=max(hlp_ppwi,wtarg(i,1))
      end do !}
c East target
      hlp_ppwo=0.0_R8
      do i=1,mpg%divFcP(maxval(mpg%strDiv),2) !{
        iFc=mpg%divFc(mpg%divFcP(maxval(mpg%strDiv),1)+i-1)
        iCv=min(mpg%fcCv(iFc,1),mpg%fcCv(iFc,2))
        s=geo%fcS(iFc)
        hlp=abs(st%dv%fht(iFc,0)+st%dv%fht(iFc,1)+
     .          st_ext%fhi(iFc,0)+st_ext%fhi(iFc,1))
        do is=0,st_ext%ns-1 !{
          taf=intface_s(iFc,mpg%nCv,mpg%nFc,mpg%fcCv,geo%fcVol,
     .     st_ext%ta(1,is))
          uaf=intface_s(iFc,mpg%nCv,mpg%nFc,mpg%fcCv,geo%fcVol,
     .     st_ext%ua(1,is))
          hlp=hlp+abs((st_ext%fa(iFc,0,is)+st_ext%fa(iFc,1,is))*(
     .     st_ext%pt(iCv,is)*ev+(1.0_R8-switch%BoRiS)*
     .    (0.5_R8*st_ext%am(is)*mp*uaf**2+taf)))
        end do !}
        hlp=hlp/s
#ifdef B25_EIRENE
        ias=b2_fnmti(iFc)
        hlp=hlp+ewldt_res(ias)
        do is=1,nnatmi
          wnutrg(i,maxval(mpg%strDiv))=wnutrg(i,maxval(mpg%strDiv))+
     .     ewlda_res(is,ias)-ewldea_res(is,ias)
        end do
        do is=1,nnmoli
          wnutrg(i,maxval(mpg%strDiv))=wnutrg(i,maxval(mpg%strDiv))+
     .     ewldm_res(is,ias)-ewldem_res(is,ias)+ewldmr_res(is,ias)
        end do
        wnutpr(i,maxval(mpg%strDiv))=ewldrp_res(ias)
#endif
        wtarg(i,maxval(mpg%strDiv))=
     .   (hlp-wnutpr(i,maxval(mpg%strDiv)))*trgshp(maxval(mpg%strDiv))
     .       +wnutrg(i,maxval(mpg%strDiv))+wrad(i,maxval(mpg%strDiv),0)
        hlp_ppwo=max(hlp_ppwo,wtarg(i,maxval(mpg%strDiv)))
      end do !}

      if(maxval(mpg%strDiv).gt.2) then !{
        hlp_pwou=0.0_R8
        do i=1,mpg%divFcP(2,2) !{
          iFc=mpg%divFc(mpg%divFcP(2,1)+i-1)
          iCv=min(mpg%fcCv(iFc,1),mpg%fcCv(iFc,2))
          s=geo%fcS(iFc)
          hlp=abs(st%dv%fht(iFc,0)+st%dv%fht(iFc,1)+
     .            st_ext%fhi(iFc,0)+st_ext%fhi(iFc,1))
          do is=0,st_ext%ns-1 !{
            taf=intface_s(iFc,mpg%nCv,mpg%nFc,mpg%fcCv,geo%fcVol,
     .       st_ext%ta(1,is))
            uaf=intface_s(iFc,mpg%nCv,mpg%nFc,mpg%fcCv,geo%fcVol,
     .       st_ext%ua(1,is))
            hlp=hlp+abs((st_ext%fa(iFc,0,is)+st_ext%fa(iFc,1,is))*(
     .       st_ext%pt(iCv,is)*ev+(1.0_R8-switch%BoRiS)*
     .      (0.5_R8*st_ext%am(is)*mp*uaf**2+taf)))
          end do !}
          hlp=hlp/s
#ifdef B25_EIRENE
          ias=b2_fnmti(iFc)
          hlp=hlp+ewldt_res(ias)
          do is=1,nnatmi
            wnutrg(i,2)=wnutrg(i,2)+
     .      ewlda_res(is,ias)-ewldea_res(is,ias)
          end do
          do is=1,nnmoli
            wnutrg(i,2)=wnutrg(i,2)+
     .       ewldm_res(is,ias)-ewldem_res(is,ias)+ewldmr_res(is,ias)
          end do
          wnutpr(i,2)=ewldrp_res(ias)
#endif
          wtarg(i,2)=(hlp-wnutpr(i,2))*trgshp(2)
     .                   +wnutrg(i,2)+wrad(i,2,0)
          hlp_pwou=max(hlp_pwou,wtarg(i,2))
        end do !}

        hlp_pwiu=0.0_R8
        do i=1,mpg%divFcP(3,2) !{
          iFc=mpg%divFc(mpg%divFcP(2,1)+i-1)
          iCv=min(mpg%fcCv(iFc,1),mpg%fcCv(iFc,2))
          s=geo%fcS(iFc)
          hlp=abs(st%dv%fht(iFc,0)+st%dv%fht(iFc,1)+
     .            st_ext%fhi(iFc,0)+st_ext%fhi(iFc,1))
          do is=0,st_ext%ns-1 !{
            taf=intface_s(iFc,mpg%nCv,mpg%nFc,mpg%fcCv,geo%fcVol,
     .       st_ext%ta(1,is))
            uaf=intface_s(iFc,mpg%nCv,mpg%nFc,mpg%fcCv,geo%fcVol,
     .       st_ext%ua(1,is))
            hlp=hlp+abs((st_ext%fa(iFc,0,is)+st_ext%fa(iFc,1,is))*(
     .       st_ext%pt(iCv,is)*ev+(1.0_R8-switch%BoRiS)*
     .      (0.5_R8*st_ext%am(is)*mp*uaf**2+taf)))
          end do !}
          hlp=hlp/s
#ifdef B25_EIRENE
          ias=b2_fnmti(iFc)
          hlp=hlp+ewldt_res(ias)
          do is=1,nnatmi
            wnutrg(i,3)=wnutrg(i,3)+
     .      ewlda_res(is,ias)-ewldea_res(is,ias)
          end do
          do is=1,nnmoli
            wnutrg(i,3)=wnutrg(i,3)+
     .       ewldm_res(is,ias)-ewldem_res(is,ias)+ewldmr_res(is,ias)
          end do
          wnutpr(i,3)=ewldrp_res(ias)
#endif
          wtarg(i,3)=(hlp-wnutpr(i,3))*trgshp(3)
     .                   +wnutrg(i,3)+wrad(i,3,0)
          hlp_pwiu=max(hlp_pwiu,wtarg(i,3))
        end do !}
      end if !}

      if(debug_flags(81).gt.0) then
        write(frm,700) maxval(mpg%strDiv)
  700   format('(',i3,'(2x,e15.7))')
C total power
        open(3,file='targb2.datv')
        do i=1,maxval(mpg%divFcP(1:maxval(mpg%strDiv),2))
          write(3,frm) (wtarg(i,itg_tmp(j)),j=1,maxval(mpg%strDiv))
        end do
        close(3)
C neutrals
        open(3,file='targb2n.datv')
        do i=1,maxval(mpg%divFcP(1:maxval(mpg%strDiv),2))
          write(3,frm) (wnutrg(i,itg_tmp(j)),j=1,maxval(mpg%strDiv))
        end do
        close(3)
C plasma
        do j=1,maxval(mpg%strDiv)
          do i=1,mpg%divFcP(j,2)
            wtarg(i,itg_tmp(j))=
     .       (wtarg(i,itg_tmp(j))-wnutrg(i,itg_tmp(j))
     .       -wrad(i,itg_tmp(j),0))/trgshp(itg_tmp(j))
          end do
        end do
        open(3,file='targb2pl.datv')
        do i=1,maxval(mpg%divFcP(1:maxval(mpg%strDiv),2))
          write(3,frm) (wtarg(i,itg_tmp(j)),j=1,maxval(mpg%strDiv))
        end do
        close(3)
C radiation
        open(3,file='targb2r.datv')
        do i=1,maxval(mpg%divFcP(1:maxval(mpg%strDiv),2))
          write(3,frm) (wrad(i,itg_tmp(j),0),j=1,maxval(mpg%strDiv))
        end do
        close(3)
      end if

c*** Calculate the peak parallel power flux to the main chamber wall

      hlp_ppfw=0.0_R8
      if(hlp_dbg) print *,'### b2usrtrc: calculating max. parallel '//
     .                    'heat flux to the main chamber wall'
      do i=1,mpg%nCg !{
        iFc=mpg%fcs_wall(i)
        iCv=min(mpg%fcCv(iFc,1),mpg%fcCv(iFc,2))
        if (mpg%fcReg(iFc).eq.0) cycle
        select case ( gridGeometry )
        case ( GEOMETRY_LINEAR )
          if (mpg%fcReg(iFc).ne.4) cycle
        case ( GEOMETRY_LIMITER )
          if (mpg%fcReg(iFc).ne.6) cycle
        case ( GEOMETRY_SN )
          if (mpg%fcReg(iFc).ne.12) cycle
        case ( GEOMETRY_CDN )
          if (mpg%fcReg(iFc).ne.18 .and. mpg%fcReg(iFc).ne.25) cycle
        case ( GEOMETRY_DDN_BOTTOM, GEOMETRY_DDN_TOP )
          if (mpg%fcReg(iFc).ne.19 .and. mpg%fcReg(iFc).ne.26) cycle
        case ( GEOMETRY_LFS_SNOWFLAKE_MINUS,
     &         GEOMETRY_LFS_SNOWFLAKE_PLUS )
          if (mpg%fcReg(iFc).ne.19) cycle
        case default
          cycle
        end select
        s=abs(geo%fcPbs(iFc))
        hlp=abs(st%dv%fhe(iFc,0)+st%dv%fhi(iFc,0)+st_ext%fhi(iFc,0)+
     .          0.5_R8*(1.0_R8-switch%BoRiS)*st%dv%fne_53(iFc,0)*
     .          (st%pl%te(mpg%fcCv(iFc,1))+st%pl%te(mpg%fcCv(iFc,2))))
        do is=0,ns-1 !{
          kef=intface_s(iFc,mpg%nCv,mpg%nFc,mpg%fcCv,geo%fcVol,
     .          st%dv%kinrgy(1,is))
          hlp=hlp+abs((st%dv%fna(iFc,0,is))*
     .     (st%rt%rpt(iCv,is)*ev+(1.0_R8-switch%BoRiS)*kef)+
     .      0.5_R8*(1.0_R8-switch%BoRiS)*st%dv%fna_53(iFc,0,is)*
     .      (st%pl%ti(mpg%fcCv(iFc,1))+st%pl%ti(mpg%fcCv(iFc,2))))
        end do !}
        hlp = abs(hlp)
        do is=0,st_ext%ns-1 !{
          taf=intface_s(iFc,mpg%nCv,mpg%nFc,mpg%fcCv,geo%fcVol,
     .     st_ext%ta(1,is))
          uaf=intface_s(iFc,mpg%nCv,mpg%nFc,mpg%fcCv,geo%fcVol,
     .     st_ext%ua(1,is))
          hlp=hlp+abs((st_ext%fa(iFc,0,is))*
     .    (st_ext%pt(iCv,is)*ev+(1.0_R8-switch%BoRiS)*
     .    (0.5_R8*st_ext%am(is)*mp*uaf**2+taf)))
        end do !}
        hlp=hlp/s
        hlp_ppfw=max(hlp_ppfw,hlp)
          if(hlp_dbg) print '(1x,a,i6,3e12.4)',
     .     "iCv,pbs,hlp_ppfw,hlp_ppfw_max =>",iCv,s,hlp,hlp_ppfw
      end do !}

#ifdef B25_EIRENE
c*** Calculate the neutral flows in PFR
      if (use_eirene.eq.0) goto 30
      do k=1,npfrgrp !{
        do j=1,nspecies !{
caap          hlp=1.0_R8/2.7e20_R8
caap          if(species_list(j).eq.'H' .or. species_list(j).eq.'D' .or.
caap     ,       species_list(j).eq.'T' .or. species_list(j).eq.'DT')
caap     ,      hlp=0.5_R8*hlp
caap no need to convert to molecular throughput here
          hlp=1.0_R8
          do i=ipfrgrp(k),jpfrgrp(k) !{
            m=lpfrgrp(i)
            if(m.lt.0) m=nnlimi-m
            if(isrftype(m).ne.-3) then !{
              if(isrftype(m).eq.0) then !}
                write (*,'(1x,a,i4,a,i4)')
     .           '#b2usrtrc: missing type of diagnostic surface. '//
     .           'Surface:',m,',  type: ',isrftype(m)
                write (*,*) 'Check for missing Eirene data!!!'
              else ! }{
                write (*,'(1x,a,i4,a,i4)')
     .           '#b2usrtrc: wrong type of diagnostic surface. '//
     .           'Surface:',m,',  type: ',isrftype(m)
                write (*,*) 'Check the DG model!!!'
              endif !}
              hlp_pfrf(k,j,1)=0.0_R8
              go to 30
            end if !}
            do ia=1,nnatmi !{
              if (latmscl(ia).eq.j) then !{
                hlp_pfrf(k,j,1)=hlp_pfrf(k,j,1)+wldna(m,ia,0)
                do im=1,nnmoli !{
                  hlp_pfrf(k,j,1)=hlp_pfrf(k,j,1)+
     +                                      wldnm(m,im,0)*mlcmp(ia,im)
                end do !}
              end if !}
            end do !}
          end do !}
          hlp_pfrf(k,j,1)=hlp_pfrf(k,j,1)*hlp
        end do !}
      end do !}
   30 continue
#endif
      do k=1,npfrgrp !{
        do j=1,nspecies !{
          hlp_pfrf(0,j,1)=hlp_pfrf(0,j,1)+hlp_pfrf(k,j,1)
        end do !}
      end do !}

#ifdef B25_EIRENE
c-----------------------------------------------------------------------
c*** Evaluate the neutral data at diagnostic surface groups
c*** Here s is the total area of a diagnostic surface group,
c*** g total flux of particles of a certain sort hitting these surfaces,
c*** e average energy of these particles,
c*** a their atomic mass in AMU, x the absorption probability.

      n=nnatmi+nnmoli
      do k=1,nntrgrp !{
        s=hlp_ardts(k)
        do j=1,n !{
          e=0.0_R8
          g=0.0_R8
          do i=intrgrp(k),jntrgrp(k) !{
            l=lntrgrp(i)
            if(j.le.nnatmi) then !{
              g=g+wldna(l,j,0)
              e=e+wldna(l,j,0)*ewlda(l,j,0)
            else !}{
              m=j-nnatmi
              g=g+wldnm(l,m,0)
              e=e+wldnm(l,m,0)*ewldm(l,m,0)
            end if !}
          end do !}
          if(g.gt.0.0_R8) e=e/g
          if(e.le.0.0_R8) e=1.0_R8
          a=hlp_awdts(j)
          x=hlp_ksdts(k,j)
          hlp_tmdts(k,j)=0.5_R8*e
          hlp_dndts(k,j)=2.6e-4_R8*g*sqrt(a/(e*s*s))
          hlp_prdts(k,j)=ev*hlp_tmdts(k,j)*hlp_dndts(k,j)
          hlp_fxdts(k,j)=g/s
          if(ndta_tmp(k,j).gt.nmnt_tmp) then !{
c            e0=0.03_R8 ! room temperature in eV ==> 0.025
            x=hlp_ksdts(k,j)
            hlp_fldts(k,j)=x*g
            hlp_stdts(k,j)=622._R8*x*s/sqrt(a)
            hlp_sedts(k,j)=hlp_fldts(k,j)/max(hlp_dndts(k,j),1.0_R8)
          end if !}
        end do !}
      end do !}
c-----------------------------------------------------------------------
c*** evaluate neutral data on PFR interface surfaces

      n=nnatmi+nnmoli
      do k=1,npfi_tmp !{
        do j=1,n !{
c*** here s accumulates the area, d the density*area,
c*** e the temperature*density*area, v the flux*area
          s=0.0_R8
          d=0.0_R8
          e=0.0_R8
          v=0.0_R8
          do i=1,npfr_cvs !{
            if (k.eq.1.and.mod(mpg%cvReg(pfr_cvs(i)),4).ne.3) cycle
            if (k.eq.2.and.mod(mpg%cvReg(pfr_cvs(i)),4).ne.0) cycle
            s=s+geo%cvSz(pfr_cvs(i))
            if(ndtm_tmp(j).gt.0) then !{
              g=dab2(pfr_cvs(i),ndtm_tmp(j),1)
              t=tab2(pfr_cvs(i),ndtm_tmp(j),1)/ev
            else !}{
              g=dmb2(pfr_cvs(i),-ndtm_tmp(j),1)
              t=tmb2(pfr_cvs(i),-ndtm_tmp(j),1)/ev
            end if !}
            d=d+geo%cvSz(pfr_cvs(i))*g
            e=e+geo%cvSz(pfr_cvs(i))*g*t
            v=v+geo%cvSz(pfr_cvs(i))*g*sqrt(max(t,0.0_R8))
          end do !}
          if(s.gt.0.0_R8 .and. g.gt.0.0_R8) then !{
            e=e/d
            v=v/d
            d=d/s
          else !}{
            e=0.0_R8
            v=0.0_R8
            d=0.0_R8
          end if !}
          hlp_dnpfs(k,j)=d
          hlp_tmpfs(k,j)=e
          hlp_vlpfs(k,j)=v*0.39e4_R8/sqrt(hlp_awdts(j))
        end do !}
      end do !}
      if(hlp_dbg) then !{
        d=0.0_R8
        e=0.0_R8
        do j=1,n !{
          d=d+hlp_dnpfs(1,j)*hlp_tmpfs(1,j)
          e=e+hlp_dnpfs(2,j)*hlp_tmpfs(2,j)
        end do !}
        print *,'prfi,prfo=',ev*d,ev*e
      endif !}
#endif
c-----------------------------------------------------------------------
c*** evaluate "Degree of detachment" and Isat for each divertor

      hlp_dod=0.0_R8
      hlp_isat=0.0_R8
      mgl3_tmp=0
      call b2xppz_nodiff(mpg%nCv, ns, st%dv%ne, st%pl%na,
     .  st%pl%te, st%pl%ti, pz, st_ext)

      if(hlp_dbg) print *,'### print_isat=',print_isat    !###

      do k=1,maxval(mpg%strDiv) !{
        d=1.0e18_R8
        do i=mpg%ifdiv(k),mpg%divFcP(k,2) !{
          iFc=mpg%divFc(mpg%divFcP(k,1)+i-1)
          iCv=min(mpg%fcCv(iFc,1),mpg%fcCv(iFc,2))
          iFt=mpg%cvFt(iCv)
          e=2.0_R8*pz(iCv)
          if(e.gt.0.0_R8.and.iFt.gt.0) then
            if (nimp.gt.0.and.k.le.maxval(mpg%strDiv)/2) then
              d=min(d,pz(imp(ft_imp(iFt)))/e)
            else
              d=min(d,pz(omp(ft_omp(iFt)))/e)
            end if
          end if
          if(hlp_dbg) print *,'e,d=',e,d  !###
        end do !}
        hlp_dod(itg_tmp(k))=d
        if(hlp_dbg) print *, 'hpl_dod(k)',hlp_dod(itg_tmp(k)) !###
        if(print_isat) then !{
          mgl3_tmp=mgl3_tmp+1
csw 03jan2017 corrected do loop to exclude neutrals from isat calculation
          do is=0,ns-1 !{
            if(.not.is_neutral(is)) then !{
              do i=1,mpg%divFcP(k,2) !{
                iFc=mpg%divFc(mpg%divFcP(k,1)+i-1)
                iCv=min(mpg%fcCv(iFc,1),mpg%fcCv(iFc,2))
                hlp_isat(itg_tmp(k))=hlp_isat(itg_tmp(k))+
     +           (st%dv%fna(iFc,0,is)+st%dv%fna(iFc,1,is))*
     *            mpg%divFcOr(mpg%divFcP(k,1)+i-1)*
     *            st%rt%rza(iCv,is)
              end do !}
            end if !}
          end do !}
          hlp_isat(itg_tmp(k))=hlp_isat(itg_tmp(k))*qe
          if(hlp_dbg) print *,'### mgl3_tmp=',mgl3_tmp,itg_tmp(k)  !###
        end if !}
      end do !}
c-----------------------------------------------------------------------
c*** find the maximum ballooning factor and its location

      al_bll_mx=0.0_R8
      al_bll_ic=0.0_R8
      call gradc_r_nodiff(mpg%nCv, mpg%nFc, mpg%nVx, 0, geo,
     . mpg, pz, pv, grpz)
      do iCv=1,mpg%nCi
        if(mpg%cvOnClosedSurface(iCv) .and. mpg%cvReg(iCv).ne.0) then
          albll(iCv)=abs(grpz(iCv))*hlp_alcff
          if(hlp_dbg) print '(1x,a,i6,2e12.4)',
     .     "iCv,al_bll_mx,albll =>",iCv,al_bll_mx,albll(iCv)
          if(al_bll_mx.lt.albll(iCv)) then
            al_bll_mx=albll(iCv)
            al_bll_ic=real(iCv,R8)
          endif
        endif
      enddo
      if(hlp_dbg) print '(1x,a,1i6,1e12.4)',
     . "all_bll_ic,albll_mx, =>",nint(al_bll_ic),al_bll_mx

c=======================================================================
      if(filedata(ltrc)) then !{

c*** Write the data for the table line on the output array

        ex=nglb_tmp.lt.ntrc_tmp
        if(ex) hlp_out(1)=tim
        if(ex) hlp_out(2:nglb_tmp+1)=equg_tmp(1:nglb_tmp)
        n=nglb_tmp+1
        if(npfrgrp.gt.0) then !{
          i0=0
          if(npfrgrp.le.1) i0=1
          do i=i0,npfrgrp !{
            do j=1,nspecies !{
              do k=1,npff_tmp !{
                ex = ex .and. n.lt.ntrc_tmp
                n=n+1
                if(ex) hlp_out(n)=equp_tmp(i,j,k)
              end do !}
            end do !}
          end do !}
        end if !}

        if(hlp_dbg) print *,'eqp_tmp:',n  !###

        do k=1,nntrgrp !{
          n=n+1
          ex=ex .and. n.le.ntrc_tmp
          if(ex) hlp_out(n)=hlp_ardts(k)
          do j=1,ndts_tmp !{
            ex = ex .and. n+ndta_tmp(k,j).le.ntrc_tmp
            do i=1,ndta_tmp(k,j) !{
              if(ex) hlp_out(n+i)=equn_tmp(k,j,i)
            end do !}
            n=n+ndta_tmp(k,j)
          end do !}
        end do !}

        if(hlp_dbg) print *,'ndta_tmp:',n  !###

        do k=1,npfi_tmp !{
          n=n+1
          ex=ex .and. n.le.ntrc_tmp
          if(ex) hlp_out(n)=hlp_w(k)
          do j=1,ndts_tmp !{
            ex = ex .and. n+npft_tmp.le.ntrc_tmp
            do i=1,npft_tmp !{
              if(ex) hlp_out(n+i)=equr_tmp(k,j,i)
            end do !}
            n=n+npft_tmp
          end do !}
        end do !}

        if(hlp_dbg) print *,'npft_tmp:',n  !###

        ex=ex .and. n+ngl2_tmp.le.ntrc_tmp
        if(ex) then !{
          do i=1,ngl2_tmp !{
            hlp_out(n+i)=equg2_tmp(i)
          end do !}
        end if !}
        n=n+ngl2_tmp

        if(hlp_dbg) print *,'ngl2_tmp:',n  !###

        ex=ex .and. n+mgl3_tmp.le.ntrc_tmp
        if(hlp_dbg) print *,'###a ex,mgl3_tmp,n,ntrc_tmp=',
     .                            ex,mgl3_tmp,n,ntrc_tmp !###
        if(ex) then !{
          do i=1,mgl3_tmp !{
            hlp_out(n+i)=equg3_tmp(i)
          end do !}
        end if !}
        n=n+mgl3_tmp
        if(hlp_dbg) print '(a,4i8)','n,mgl3_tmp,ngl4_tmp,ntrc_tmp:',
     ,                               n,mgl3_tmp,ngl4_tmp,ntrc_tmp
        ex=ex .and. n+ngl4_tmp.le.ntrc_tmp
        if(ex) then !{
          do i=1,ngl4_tmp !{
            hlp_out(n+i)=equg4_tmp(i)
          end do !}
        end if !}
        n=n+ngl4_tmp

        if(hlp_dbg) print *,'ngl4_tmp:',n  !###

        if(.not.ex .or. n.ne.ntrc_tmp) go to 96
      end if !}

      if(first_out.and.hlp_dbg) then !{
        print *
        print *,'b2usrtc: Eirene data. nnatmi,nnmoli=',nnatmi,nnmoli
        print *,'eirtxt:'
        print '(10a10)',(eirtxt(i),i=1,nnatmi+nnmoli)
caap tbd add nphoti
        print *,'Species: ',ndts_tmp,(ndtm_tmp(j),j=1,ndts_tmp)
        do j=1,ndts_tmp !{
          do i=1,ndta_tmp(k,j) !{
            if(ndtm_tmp(j).gt.0) then !{
              hlp_hdr(j)=adjustl(eirtxt(ndtm_tmp(j)))
caap tbd add nphoti
            else !}{
              hlp_hdr(j)=adjustl(eirtxt(nnatmi-ndtm_tmp(j)))
caap tbd add nphoti
            end if !}
          end do !}
        end do !}
        print *,'Diagnostic surface groups:'
        do k=1,nntrgrp !{
          print *
          print '(1x,a8,i9,9i10/(8x,10i10))',gntrgrp(k),
     ,                            (lntrgrp(i),i=intrgrp(k),jntrgrp(k))

          print '(3x,a5,4x,1p,10e10.2)','area',
     ,                    (wlarea(lntrgrp(i)),i=intrgrp(k),jntrgrp(k))
          do j=1,ndts_tmp !{
            if(j.le.nnatmi) then !{
              print '(3x,a5,1x,a2,1x,1p,10e10.2)','wldna',hlp_hdr(j),
     ,                  (wldna(lntrgrp(i),j,0),i=intrgrp(k),jntrgrp(k))
              print '(3x,a5,1x,a2,1x,1p,10e10.2)','ewlda',hlp_hdr(j),
     ,                  (ewlda(lntrgrp(i),j,0),i=intrgrp(k),jntrgrp(k))
            else !}{
              print '(3x,a5,1x,a2,1x,1p,10e10.2)','wldnm',hlp_hdr(j),
     ,           (wldnm(lntrgrp(i),j-nnatmi,0),i=intrgrp(k),jntrgrp(k))
              print '(3x,a5,1x,a2,1x,1p,10e10.2)','ewldm',hlp_hdr(j),
     ,           (ewldm(lntrgrp(i),j-nnatmi,0),i=intrgrp(k),jntrgrp(k))
            end if !}
          end do !}
        end do !}
      end if !}
c=======================================================================
      if(first_out .and. filedata(ltrc)) then !{

c*** Prepare the headers, markers, and padding data for the output file

        ex=nglb_tmp.lt.ntrc_tmp
        if(ex) hlp_hdr(1)='time'
        if(ex) hlp_hdr(2:nglb_tmp+1)=hlp_hg
        n=nglb_tmp+1
        i0=0
        if(npfrgrp.le.1) i0=1
        do i=i0,npfrgrp !{
          hlp_hb='total'
          if(i.gt.0) hlp_hb=gpfrgrp(i)
          do j=1,nspecies !{
            do k=1,npff_tmp !{
              ex = ex .and. n.lt.ntrc_tmp
              n=n+1
              if(ex) hlp_hdr(n)=adjustl(trim(hlp_hp(k))//'_'//
     /                             trim(hlp_hb)//'_'//species_list(j))
            end do !}
          end do !}
        end do !}

        n_pad_str=0
        l_pad_str=n+1

        do k=1,nntrgrp !{
          n=n+1
          hlp_hdr(n)='Area_'//gntrgrp(k)
          do j=1,ndts_tmp !{
            ex = ex .and. n+ndta_tmp(k,j).le.ntrc_tmp
            do i=1,ndta_tmp(k,j) !{
              if(ndtm_tmp(j).gt.0) then !{
                hlp_m=adjustl(eirtxt(ndtm_tmp(j)))
              else !}{
                hlp_m=adjustl(eirtxt(nnatmi-ndtm_tmp(j)))
              end if !}
              if(ex) hlp_hdr(n+i)=adjustl(trim(hlp_hn(i))//'_'//
     /                                    trim(hlp_m)//'_'//gntrgrp(k))
            end do !}
            n=n+ndta_tmp(k,j)
          end do !}
        end do !}
        do k=1,npfi_tmp !{
          n=n+1
          hlp_hdr(n)='PFR'//trim(hlp_hs(k))//'_area'
          do j=1,ndts_tmp !{
            ex = ex .and. n+npft_tmp.le.ntrc_tmp
            do i=1,npft_tmp !{
              if(ndtm_tmp(j).gt.0) then !{
                hlp_m=adjustl(eirtxt(ndtm_tmp(j)))
caap tbd add nphoti
              else !}{
                hlp_m=adjustl(eirtxt(nnatmi-ndtm_tmp(j)))
caap tbd add nphoti
              end if !}
              if(ex) hlp_hdr(n+i)='PFR'//trim(hlp_hs(k))//'_'//
     /             trim(hlp_hr(i))//'_'//trim(hlp_m)
            end do !}
            n=n+npft_tmp
          end do !}
        end do !}
        ex=ex .and. n+ngl2_tmp.le.ntrc_tmp
        if(ex) hlp_hdr(n+1:n+ngl2_tmp)=hlp_hg2
        n=n+ngl2_tmp
        ex=ex .and. n+mgl3_tmp.le.ntrc_tmp
        if(hlp_dbg) print *,'###b ex,mgl3_tmp,n,ntrc_tmp=',
     .                            ex,mgl3_tmp,n,ntrc_tmp !###
        if(ex .and. mgl3_tmp.gt.0) then !{
          i=0
          do k=1,maxval(mpg%strDiv) !{
            n=n+1
            i=i+1
            hlp_hdr(n)=hlp_hg3(itg_tmp(k))
          end do !}
          if(i.ne.mgl3_tmp) then !{
            print *,'### b2usrtrc: i.ne.mgl3_tmp: ',i,mgl3_tmp,
     ,         '  ==> no Isat data in the file'
            print_isat=.false.
            mgl3_tmp=0
            n=n-i
            ntrc_tmp=n
          end if !}
        end if !}
        if(ex .and. ngl4_tmp.gt.0) then !{
          ex=ex .and. n+ngl4_tmp.le.ntrc_tmp
          hlp_hdr(n+1:n+ngl4_tmp)=hlp_hg4
          n=n+ngl4_tmp
        end if !}

        n_pad_str=n-l_pad_str+1
        if(n_pad_str.eq.0) l_pad_str=1

        if(hlp_dbg) print *,' n,hlp_hdr= ',
     .                        n,('  ',trim(hlp_hdr(i)),i=1,n) !###

        if(.not.ex .or. n.ne.ntrc_tmp) go to 98

        mrks_hlp(1)=lv
        mrks_hlp(2)=nglb_tmp
        mrks_hlp(3)=npff_tmp
        mrks_hlp(4)=nspecies
        mrks_hlp(5)=npfrgrp
        mrks_hlp(6)=nntrgrp
        mrks_hlp(7)=ngl2_tmp
        mrkn_hlp=mrks_hlp

c*** Initialize the file for tracing data.

c*** Moved here from b2usrtri to ensure being called after the first
c*** call to Eirene.

c*** lun   :  logical unit number for the output file (=ntrc(ltrc))
c*** f     :  the associated file name

        lun=ntrc(ltrc)
        f=fn_ntrc(ltrc)

        n=nrkn_hlp
        m=0
        xdum(1)=0.0_R8
        call combfile(lun,f,lenfnm,hlp_hdr,ntrc_tmp,label,.false.,
     ,                fix_user,hlp_hdr(l_pad_str),hlp_out(l_pad_str),
     ,                n_pad_str,7,4,tab,mrks_hlp,n,xdum,m,hlp_frm,irt)

        if(hlp_dbg) print *,'b2usrtrc: after combfile. irt,m,f=',  !###
     ,                                         irt,m,'  ',trim(f)  !###

        if(irt.lt.0) then !{
          go to 92
        else if(irt.gt.0) then !}{
          if(irt.ge.8) then !{
            if(mrks_hlp(4).ne.mrkn_hlp(4)) go to 93
          end if !}
          if(irt.ge.16) then !{
            print *,'b2usrtrc: new output file opened - ',trim(f)
          end if !}
          fn_ntrc(ltrc)=f
        end if !}
        if(hlp_dbg) print *,'hlp_frm=',trim(hlp_frm)  !###
        first_out=.false.
      end if !}
c=======================================================================

      if(ank_tracing.eq.0) return
      if(filedata(ltrc)) then !{

c*** Write the data in a separate file for plotting...

        if(hlp_dbg) then !{
          inquire(lun,opened=ex,name=f)
          print *,'b2usrtrc: unit ',lun,', opened=',ex,
     ,                                           ', file ',trim(f)
        endif !}

        write(lun,hlp_frm,err=94) (hlp_out(j),j=1,ntrc_tmp)
caap        if(npfrgrp.gt.0) then !{
caap          i0=0
caap          if(npfrgrp.le.1) i0=1
caap          write(lun,hlp_frm,err=94) tim,(equg_tmp(j),j=1,ng_tmp),
caap     ,      (((equp_tmp(i,j,k),i=i0,npfrgrp),j=1,nspecies),k=1,npff_tmp)
caap        else !}{
caap          write(lun,hlp_frm,err=94) tim,(equg_tmp(j),j=1,ng_tmp)
caap        end if !}
      end if !}
      return
c=======================================================================
   92 write(*,*) 'b2usrtrc: wrong format of the file ',trim(f)
      filedata(ltrc)=.false.
      return

   93 write(*,'(1x,4a)') 'b2usrtrc: number of species changed - ',
     ,                'remove the ',trim(f),' tracing file first!'
      filedata(ltrc)=.false.
      return

   94 write(*,*) 'b2usrtrc: cannot write in the file ',trim(f)
      filedata(ltrc)=.false.
      return

 96   write(*,'(2a,2i8)') 'b2usrtrc: number of columns ',
     ,  'in the header and data differ:',ntrc_tmp,n
      filedata(ltrc)=.false.
      return

   98 write(*,'(1x,2a,2i8)') 'b2usrtrc: internal error. ',
     ,  'Number of columns in the header - ',ntrc_tmp,n
      filedata(ltrc)=.false.
      return
      end subroutine b2usrtrc
c=======================================================================

      subroutine b2usrtri(mpg,geo,ns,switch)
      use b2mod_switches_diffv
      use b2us_map_diffv
      use b2us_geo_diffv
      implicit none
      type (mapping), intent (in) :: mpg
      type (geometry), intent (in) :: geo
      type (switches), intent (in) :: switch
      integer, intent(in) :: ns
      integer i, j, k, l, i0, nv, na_tmp, np_tmp, ns_tmp
      character hlp_c, f*(lenfnm), hlp_hb*8
      integer lnh_tmp
      parameter (lnh_tmp=20)
      character hlp_ha(0:npfrgrpd,ns,npff_tmp)*(lnh_tmp)
      logical ex, warn_tmp
c-----------------------------------------------------------------------
c*** Initialise the file for tracing data
c-----------------------------------------------------------------------

      tab=char(9)
      gridGeometry = geometryId( mpg, geo, 1 )
      plasmaGeometry = geometryId( mpg, geo, 2 )

c*** lun   :  logical unit number for the output file (=ntrc(ltrc))
c*** f     :  the associated file name

      lun=ntrc(ltrc)
      f=fn_ntrc(ltrc)
      i0=0                                                           !srv 05.11.10
      if(npfrgrp.le.1) i0=1                                          !srv 05.11.10

c*** assign target indices
      itg_tmp = 0
      select case ( plasmaGeometry )
      case ( GEOMETRY_LINEAR, GEOMETRY_LIMITER, GEOMETRY_SN,
     &       GEOMETRY_STELLARATORISLAND )
        itg_tmp(1) = 1
        itg_tmp(2) = 2
      case ( GEOMETRY_CDN, GEOMETRY_DDN_BOTTOM, GEOMETRY_DDN_TOP )
        itg_tmp(1) = 1
        itg_tmp(2) = 3
        itg_tmp(3) = 4
        itg_tmp(4) = 2
      case ( GEOMETRY_LFS_SNOWFLAKE_MINUS, GEOMETRY_LFS_SNOWFLAKE_PLUS )
        itg_tmp(1) = 1
        itg_tmp(2) = 4
        itg_tmp(3) = 3
        itg_tmp(4) = 2
      end select

      call find_file(f,ex)
      open(lun,file=f,err=90)
      rewind(lun)
      if(ex) read(lun,*,end=91) hlp_c
      goto 95
   91 ex = .false.
   95 continue
      rewind(lun)
      if(ex.and.switch%b2mndr_stim.lt.0.0_R8) then !{
        na_tmp=0
        np_tmp=0
        ns_tmp=0
        read(lun,*,err=92) hlp_c,nv,ng_tmp
        lv_trc(ltrc)=nv
        if(nv.ge.20000727) then !{
          rewind(lun)
          read(lun,*,err=92) hlp_c,nv,ng_tmp,np_tmp,ns_tmp,na_tmp
        end if !}
        warn_tmp=.false.
        if(nv.ne.lv) then !{
          warn_tmp=.true.
          write(*,'(1x,a,i8,3a,i8)')
     ,     '!!! b2usrtrc: current version ',lv,
     ,     ' differs from that in the file ',trim(f),' : ',nv
          if(nv .lt. lvmj) then !{
c*** It was a major change of the file format - impossible to continue
            write(*,*) 'No tracing file - ',
     ,                 'remove the existing file first: ',trim(f)
            filedata(ltrc)=.false.
            return
          end if !}
          if(nv .lt. 20000727 .and. npfrgrp.gt.0) then !{
            write(*,*) 'No tracing file - ',
     ,                 'remove the existing file first: ',trim(f)
            filedata(ltrc)=.false.
            return
          end if !}
        end if !}
        if(na_tmp.ne.npfrgrp) then !{
          write(*,'(1x,a,i2,3a,i2)')
     ,     '!!! b2usrtrc: current number of PFR surface groups ',
     ,     npfrgrp,
     ,     ' differs from that in the file ',trim(f),' : ',na_tmp
          filedata(ltrc)=.false.
        end if !}
        if(npfrgrp.gt.0) then !{
          if(ng_tmp.ne.nglb_tmp) then !{
            write(*,'(1x,a,i4,3a,i4)')
     ,       '!!! b2usrtrc: current number of global values ',nglb_tmp,
     ,       ' differs from that in the file ',trim(f),' : ',ng_tmp
            filedata(ltrc)=.false.
          end if !}
          if(np_tmp.ne.npff_tmp) then !{
            write(*,'(1x,a,i2,3a,i2)')
     ,       '!!! b2usrtrc: current number of PFR values ',npff_tmp,
     ,       ' differs from that in the file ',trim(f),' : ',np_tmp
            filedata(ltrc)=.false.
          end if !}
          if(ns_tmp.ne.nspecies) then !{
            write(*,'(1x,a,i2,3a,i2)')
     ,       '!!! b2usrtrc: current number of species ',nspecies,
     ,       ' differs from that in the file ',trim(f),' : ',ns_tmp
            filedata(ltrc)=.false.
          end if !}
        end if !}

        if(.not.filedata(ltrc)) then !{
          write(*,*) 'No ',trim(f),' tracing file - ',
     ,                                'remove the existing file first!'
          return
        end if !}
        if(warn_tmp) write(*,*) ' You should remove the file!'
c*** Read to the end of the file
   50   read(lun,'(a1)',end=55) hlp_c
        go to 50
   55   backspace lun
      else !}{
c*** New file
        ng_tmp=nglb_tmp
        np_tmp=npff_tmp
        ns_tmp=nspecies
        na_tmp=npfrgrp
        do i=0,npfrgrp !{
          hlp_hb='total'
          if(i.gt.0) hlp_hb=gpfrgrp(i)
          do j=1,nspecies !{
            do k=1,npff_tmp !{
              hlp_ha(i,j,k)=hlp_hp(k)//'_'//hlp_hb//'_'//species_list(j)
              call lefta(hlp_ha(i,j,k),lnh_tmp,l)
            end do !}
          end do !}
        end do !}
        write(lun,700) lv,ng_tmp,np_tmp,ns_tmp,na_tmp
        write(lun,'(2a)') '#  ',label
        write(hlp_frm,'(a,i3,a)')
     ,   '(''# time'',8x,',ng_tmp+np_tmp*(na_tmp+1-i0)*ns_tmp,'(3x,a))'
#ifdef DBG
        print *,'b2usrtrc frm:',trim(hlp_frm)  !###
#endif
        if(npfrgrp.gt.0) then !{
          write(lun,hlp_frm) (hlp_hg(i),i=1,ng_tmp),
     ,        (((hlp_ha(i,j,k),i=i0,npfrgrp),j=1,nspecies),k=1,npff_tmp)
        else !}{
          write(lun,hlp_frm) (hlp_hg(i),i=1,ng_tmp)
        end if !}
      end if !}
      write(hlp_frm,'(a,i3,a)')
     - '(1p,e16.7,',ng_tmp+np_tmp*(na_tmp+1-i0)*ns_tmp,'(2x,e12.4e3))'
      return
c=======================================================================
   90 write(*,*) 'b2usrtrc: cannot open file ',trim(f)
      filedata(ltrc)=.false.
      return

   92 write(*,*) 'b2usrtrc: wrong format of the file ',trim(f)
      filedata(ltrc)=.false.
      return

c=======================================================================
  700 format('#   ',5i8)
      end subroutine b2usrtri
c======================================================================
c
      subroutine b2usrtrcc(mpg,iline)
      use b2us_map_diffv
      implicit none
      type(mapping), intent(in) :: mpg
      integer, intent(out) :: iline
      integer :: i, j, k, l, m, n, i0, irt, l_pad_str, mgl3_tmp
      real (kind=R8) :: tdummy, x
      character hlp_c, hlp_hb*8, f*(lenfnm), hlp_m*16
      logical ex
c----------------------------------------------------------------------
c*** - to be called by b2co for locating truncation for savefile
c----------------------------------------------------------------------
c
      if(.not.filedata(ntrc_usr)) return
      ndta_tmp=nmnt_tmp
      ndts_tmp=nnatmi+nnmoli
      do i=1,ndts_tmp !{
        if(i.le.nnatmi) then !{
          ndtm_tmp(i)=i
        else !}{
          ndtm_tmp(i)=nnatmi-i
        end if !}
      end do !}
      do j=1,ndts_tmp !{
        if(ndtm_tmp(j).gt.0) then !{
          m=ndtm_tmp(j)
        else !}{
          m=natm-ndtm_tmp(j)
        end if !}
        ndtx_tmp(j)=m
      end do !}
      do k=1,nntrgrp !{
        do i=intrgrp(k),jntrgrp(k) !{
          l=lntrgrp(i)
          do j=1,ndts_tmp !{
            m=ndtx_tmp(j)
            x=wlabsrp(m,lntrgrp(intrgrp(k)))
            if(i.eq.intrgrp(k)) then !{
              if(x.gt.0.0_R8) ndta_tmp(k,j)=ndtk_tmp
              hlp_ksdts(k,j)=x
            else !}{
              if(x.ge.0.0_R8 .and. wlabsrp(m,l).ne.x) then !{
                ndta_tmp(k,j)=nmnt_tmp
              end if !}
            end if !}
          end do !}
        end do !}
      end do !}
      lun=ntrc(ltrc)
      f=fn_ntrc(ltrc)
      ntrc_tmp=1+nglb_tmp+nntrgrp+
     + npfi_tmp*(1+npft_tmp*ndts_tmp)+ngl2_tmp+ngl4_tmp
      k=1+npfrgrp
      if(npfrgrp.le.1) k=max(0,npfrgrp)
      ntrc_tmp=ntrc_tmp+nspecies*npff_tmp*k
      do k=1,nntrgrp !{
        do j=1,ndts_tmp !{
          ntrc_tmp=ntrc_tmp+ndta_tmp(k,j)
        end do !}
      end do !}
      if(print_isat) then !{
        if(plasmaGeometry.ne.GEOMETRY_ANNULUS) ntrc_tmp=ntrc_tmp+2
        if(plasmaGeometry.eq.GEOMETRY_CDN .or.
     &       plasmaGeometry.eq.GEOMETRY_DDN_BOTTOM .or.
     &       plasmaGeometry.eq.GEOMETRY_DDN_TOP .or.
     &       plasmaGeometry.eq.GEOMETRY_LFS_SNOWFLAKE_MINUS .or.
     &       plasmaGeometry.eq.GEOMETRY_LFS_SNOWFLAKE_PLUS)
     &   ntrc_tmp=ntrc_tmp+2
      end if !}
      allocate(hlp_hdr(ntrc_tmp))
      allocate(hlp_out(ntrc_tmp))
      ex=nglb_tmp.lt.ntrc_tmp
      if(ex) hlp_hdr(1)='time'
      if(ex) hlp_hdr(2:nglb_tmp+1)=hlp_hg
      n=nglb_tmp+1
      i0=0
      if(npfrgrp.le.1) i0=1
      do i=i0,npfrgrp !{
        hlp_hb='total'
        if(i.gt.0) hlp_hb=gpfrgrp(i)
        do j=1,nspecies !{
          do k=1,npff_tmp !{
            ex = ex .and. n.lt.ntrc_tmp
            n=n+1
            if(ex) hlp_hdr(n)=adjustl(trim(hlp_hp(k))//'_'//
     /                             trim(hlp_hb)//'_'//species_list(j))
          end do !}
        end do !}
      end do !}
      do i=1,ndts_tmp !{
        if(i.le.nnatmi) then !{
          ndtm_tmp(i)=i
        else !}{
          ndtm_tmp(i)=nnatmi-i
        end if !}
      end do !}
      n_pad_str=0
      l_pad_str=n+1
      do k=1,nntrgrp !{
        n=n+1
        hlp_hdr(n)='Area_'//gntrgrp(k)
        do j=1,ndts_tmp !{
          ex = ex .and. n+nmnt_tmp.le.ntrc_tmp
          do i=1,nmnt_tmp !{
            if(ndtm_tmp(j).gt.0) then !{
              hlp_m=adjustl(textan(ndtm_tmp(j)-1))
            else !}{
              hlp_m=adjustl(textmn(-ndtm_tmp(j)-1))
            end if !}
            if(ex) hlp_hdr(n+i)=adjustl(trim(hlp_hn(i))//'_'//
     /                                  trim(hlp_m)//'_'//gntrgrp(k))
          end do !}
          n=n+nmnt_tmp
        end do !}
      end do !}
      do k=1,npfi_tmp !{
        n=n+1
        hlp_hdr(n)='PFR'//trim(hlp_hs(k))//'_area'
        do j=1,ndts_tmp !{
          ex = ex .and. n+npft_tmp.le.ntrc_tmp
          do i=1,npft_tmp !{
            if(ndtm_tmp(j).gt.0) then !{
              hlp_m=adjustl(textan(ndtm_tmp(j)-1))
            else !}{
              hlp_m=adjustl(textmn(-ndtm_tmp(j)-1))
            end if !}
            if(ex) hlp_hdr(n+i)='PFR'//trim(hlp_hs(k))//'_'//
     /                                 trim(hlp_hr(i))//'_'//trim(hlp_m)
          end do !}
          n=n+npft_tmp
        end do !}
      end do !}
      ex=ex .and. n+ngl2_tmp.le.ntrc_tmp
      if(ex) hlp_hdr(n+1:n+ngl2_tmp)=hlp_hg2
      n=n+ngl2_tmp
      mgl3_tmp=0
      do k=1,maxval(mpg%strDiv) !{
        if(print_isat) mgl3_tmp=mgl3_tmp+1
      end do !}
      ex=ex .and. n+mgl3_tmp.le.ntrc_tmp
      if(ex .and. mgl3_tmp.gt.0) then !{
        i=0
        do k=1,maxval(mpg%strDiv) !{
          n=n+1
          i=i+1
          hlp_hdr(n)=hlp_hg3(itg_tmp(k))
        end do !}
        if(i.ne.mgl3_tmp) then !{
          print *,'### b2usrtrc: i.ne.mgl3_tmp: ',i,mgl3_tmp,
     ,            '  ==> no Isat data in the file'
          print_isat=.false.
          mgl3_tmp=0
          n=n-i
          ntrc_tmp=n
        end if !}
      end if !}
      if(ex .and. ngl4_tmp.gt.0) then !{
        ex=ex .and. n+ngl4_tmp.le.ntrc_tmp
        hlp_hdr(n+1:n+ngl4_tmp)=hlp_hg4
        n=n+ngl4_tmp
      end if !}
      n_pad_str=n-l_pad_str+1
      if(n_pad_str.eq.0) l_pad_str=1
      n=nrkn_hlp
      m=0
      xdum(1)=0.0_R8
      mrks_hlp(1)=lv
      mrks_hlp(2)=nglb_tmp
      mrks_hlp(3)=npff_tmp
      mrks_hlp(4)=nspecies
      mrks_hlp(5)=npfrgrp
      mrks_hlp(6)=nntrgrp
      mrks_hlp(7)=ngl2_tmp
      call combfile(lun,f,lenfnm,hlp_hdr,ntrc_tmp,label,.false.,
     ,              fix_user,hlp_hdr(l_pad_str),hlp_out(l_pad_str),
     ,              n_pad_str,7,4,tab,mrks_hlp,n,xdum,m,hlp_frm,irt)
      print *,'b2usrtrcc: after combfile. irt,m,f=',          !###
     ,                                    irt,m,'  ',trim(f)  !###
      if(irt.lt.0) then !{
        go to 92
      else if(irt.gt.0) then !}{
        if(irt.ge.8) then !{
          if(mrks_hlp(4).ne.mrkn_hlp(4)) go to 93
        end if !}
        if(irt.ge.16) then !{
          print *,'b2usrtrc: new output file opened - ',trim(f)
        end if !}
        fn_ntrc(ltrc)=f
      end if !}
      call find_file(f,ex)
      open(lun,file=f,err=190)
      rewind(lun)
      iline = 0
      read(lun,'(a1)',err=155,end=155) hlp_c
      do while (hlp_c.eq.'#')
        iline = iline+1
        read(lun,'(a1)',err=155,end=155) hlp_c
      end do
      backspace lun
      read(lun,*,err=155,end=155) tdummy
      do while (tdummy.le.tim)
        iline = iline+1
        read(lun,*,err=155,end=155) tdummy
      end do
  155 continue
!      close(lun)
      return
c======================================================================
   92 write(*,*) 'b2usrtrc: wrong format of the file ',trim(f)
      filedata(ltrc)=.false.
      return

   93 write(*,'(1x,4a)') 'b2usrtrc: number of species changed - ',
     ,                'remove the ',trim(f),' tracing file first!'
      filedata(ltrc)=.false.
      return

  190 write(*,*) 'b2usrtrc: cannot open file ',trim(f)
      filedata(ltrc)=.false.
      iline=-1
      return
      end subroutine b2usrtrcc
c=======================================================================
      subroutine b2usrtrs

      return
      end subroutine b2usrtrs
c=======================================================================
      subroutine b2usrtrf
c     Finalization step
      logical file_open

      lun=ntrc(ltrc)
      inquire(unit=lun,opened=file_open)
      if (file_open) close(lun)

      if (.not.allocated(hlp_hdr)) return
      deallocate(hlp_hdr)
      deallocate(hlp_out)

      return
      end subroutine b2usrtrf
c=======================================================================
      end module b2mod_usrtrc
