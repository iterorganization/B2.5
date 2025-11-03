      module b2mod_usrtrc
      use b2mod_types
      use b2mod_diag_diffv_diffv
      use b2mod_time
      use b2mod_b2cmpa_diffv_diffv
      use b2mod_user_namelist_diffv_diffv
#ifdef B25_EIRENE
      use eirmod_parmmod
     .   , only : nstra, natm
      use eirmod_wneutrals
      use eirmod_eirbra
#else
      use b2mod_eirbra
#endif
      use b2mod_dimensions
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
c***                         (the peaks do not include radiation load)
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
c-----------------------------------------------------------------------

c*** Global values
      integer nglb_tmp, ng_tmp
      parameter (nglb_tmp=22)

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
      parameter (ngl4_tmp=3)
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
      data itg_tmp / 4*-2 /
      real (kind=R8), save :: hlp_ksdts(nntrgrpd,DEF_NATM+DEF_NMOL)
c*** Column headers:
      data hlp_hg /
     1      'Pneu_PFR', 'Pneu_PFi', 'Pneu_PFo', 'He_fract', 'He_enrch',
     2      'He_p_prs', 'pk_pwr_i', 'pk_pwr_o', 'pk_pw_iu', 'pk_pw_ou',
     3      'Sp_DT   ', 'Sp_He   ', 'Sp_He_sc', 'He_src  ', 'He_sc_cr',
     4      'He_sc_sp', 'He_og_cr', 'He_og_sp', 'He_enrcs', 'He_enrcc',
     4      'DT_thrpt', 'He_thrpt'/
      data hlp_hg2 /
     1      'DoD_il', 'DoD_ol', 'DoD_iu', 'DoD_ou' /
      data hlp_hg3 /
     1      'isat_tt_il', 'isat_tt_ol', 'isat_tt_iu', 'isat_tt_ou' /
      data hlp_hg4 /
     1      'al_bll_max', 'al_bll_ix', 'al_bll_iy' /
      data hlp_hp / 'PFR_flow' /
      data hlp_hn / 'n','p','T','j','G','S','s' /
      data hlp_hr / 'n', 'T', 'v' /
      data hlp_hs / 'i', 'o' /

      logical double_null
      data double_null / .false. /

      integer, save :: jxi, jxa, jsep, n_pad_str
      real (kind=R8), save :: xdum(1)
      character, save :: tab, hlp_frm*120

      external fix_user

      contains

      subroutine b2usrtrc(nx,ny,ns,nxtl,nxtr,BoRiS)
      use b2mod_geo_diffv_diffv
      use b2mod_rates
      use b2mod_plasma_diffv_diffv
      use b2mod_external_diffv_diffv
      use b2mod_indirect_diffv_diffv
      use b2mod_constants
      use b2mod_boundary_namelist_diffv_diffv
      use b2mod_numerics_namelist_diffv_diffv
#ifdef B25_EIRENE
      use eirmod_extrab25
#endif

      implicit none

      integer nx, ny, ns, nxtl, nxtr
      real (kind=R8) :: BoRiS
      real (kind=R8) ::
     *     equg_tmp(nglb_tmp)
      real (kind=R8) ::
     *     hlp_ppfr, hlp_ppfi, hlp_ppfo, hlp_hefr, hlp_henr, hlp_hepr,
     ,     hlp_ppwi, hlp_ppwo, hlp_pwiu, hlp_pwou, hlp_spdt, hlp_sphe,
     ,     hlp_sshe, hlp_hesr, hlp_hesc, hlp_hess, hlp_heoc, hlp_heos,
     ,     hlp_hens, hlp_henc, hlp_dttp, hlp_hetp
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
     ,     (equg_tmp(21), hlp_dttp), (equg_tmp(22), hlp_hetp)
      real (kind=R8) ::
     *     equg2_tmp(ngl2_tmp)
      real (kind=R8) ::
     *     hlp_dod(4), hlp_dodil, hlp_dodol,hlp_dodiu, hlp_dodou
      equivalence (equg2_tmp(1), hlp_dod),
     ,     (equg2_tmp(1),  hlp_dodil), (equg2_tmp(2),  hlp_dodol),
     ,     (equg2_tmp(3),  hlp_dodiu), (equg2_tmp(4),  hlp_dodou)
      integer mgl3_tmp
      integer kgl3_tmp(ngl3_tmp)
      data kgl3_tmp / -1, 0, 0, -1 /
      real (kind=R8) ::
     *     equg3_tmp(ngl3_tmp)
      real (kind=R8) ::
     *     hlp_isat(4), hlp_istil, hlp_istol, hlp_istiu, hlp_istou
      equivalence (equg3_tmp(1), hlp_isat),
     ,     (equg3_tmp(1),  hlp_istil), (equg3_tmp(2),  hlp_istol),
     ,     (equg3_tmp(3),  hlp_istiu), (equg3_tmp(4),  hlp_istou)
c     integer kgl4_tmp(2)
      real (kind=R8) ::
     *     equg4_tmp(ngl4_tmp)
      real (kind=R8) ::
     *     hlp_albll(ngl4_tmp),al_bll_mx,al_bll_ix,al_bll_iy
      equivalence (equg4_tmp(1), hlp_albll),
     ,     (equg4_tmp(1),  al_bll_mx), (equg4_tmp(2),  al_bll_ix),
     ,     (equg4_tmp(3),  al_bll_iy)

c*** PFR flow values
      real (kind=R8) ::
     *     equp_tmp(0:npfrgrpd,DEF_NATM,npff_tmp)
      real (kind=R8) ::
     *     hlp_pfrf(0:npfrgrpd,DEF_NATM,npff_tmp)
      real (kind=R8) ::
     *     wtarg_max(max(4,DEF_NCUT)) !vk, max target heat flux density
      equivalence (equp_tmp(0,1,1), hlp_pfrf)
      character hlp_hb*8

c***    atomic/molecular weight, group area, absorption
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

c***    number of Eirene surfaces per interface, list of these surfaces,
c***    starting and ending B2 cell for interface
      integer npfs_tmp(npfi_tmp),lpfs_tmp(nstra,npfi_tmp),
     ,  ipfc_tmp(npfi_tmp),jpfc_tmp(npfi_tmp)
c***    surface area,cell interface areas
      real (kind=R8) ::
     *  hlp_arpfs(npfi_tmp)
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
     *  hlp, hlp_s, hlp_pi, hlp_si, hlp_po, hlp_so,
     ,  hlp_he, hlp_hep, hlp_neu_peak,
     ,  hlp_ne_he_c, hlp_ne_he_s
      real (kind=R8) ::
     *  d, e, s, x
      real (kind=R8) :: pz(-1:nx,-1:ny)
      integer i, j, k, l, m, n, i0, is, ix, iy, js, irt,
     ,  j_he_i, l_pad_str
#ifdef B25_EIRENE
      integer ia, im, iyp, ls
      real (kind=R8) :: a, g, v, hlp_h, hlp_p, hlp_tot, hlpnw
c     real (kind=R8) :: e0
#endif
c*** radial position of the separatrix and poloidal positions of
c*** targets and midplane for DoD calculation for each target
      integer isp_tmp(4),imp_tmp(4)
      integer use_eirene
      logical ex, first, first_out
      character f*(lenfnm)
      character hlp_tpl*2
      data hlp_neu_peak / 5.0_R8 /
      data use_eirene / 0 /
      data first,first_out / 2*.true. /
c-----------------------------------------------------------------------

      save
c=======================================================================
      if(.not.filedata(ltrc)) return
c=======================================================================

      if (first) then !{

c!###{
c        print *
c        print *,'b2usrtrc: before initialisation. ',
c     ,    'natm,nmol,nnatmi,nnmoli,nntrgrp=',
c     ,     natm,nmol,nnatmi,nnmoli,nntrgrp
c        if(nntrgrp.gt.0) then !{
c          print '(1x,a9,2x,20i4/(12x,20i4))','lntrgrp:',
c     ,                                (lntrgrp(i),i=1,jntrgrp(nntrgrp))
c          print '(1x,a9,2x,16a9)','gntrgrp:',
c     ,                                (gntrgrp(i),i=1,nntrgrp)
c          print '(1x,a9,1x,i4,15i9)','intrgrp:',
c     ,                                (intrgrp(i),i=1,nntrgrp)
c          print '(1x,a9,1x,i4,15i9)','jntrgrp:',
c     ,                                (jntrgrp(i),i=1,nntrgrp)
c          print *
c        end if !}
c!###}
        tab=char(9)
        call ipgeti ('b2mndr_eirene', use_eirene)
c-----------------------------------------------------------------------
c*** Check the topology

        hlp_tpl='UD'
        if(nncut.eq.1) then !{
c*** Single null
          hlp_tpl='SN'
        else if(nncut.eq.2) then !}{
c*** Double null
          hlp_tpl='DN'
        else if(nncut.eq.0) then !}{
c*** Limiter case
          hlp_tpl='LM'
        end if !}
        if(hlp_tpl.eq.'SN') then !{
          double_null=.false.
        else if(hlp_tpl.eq.'DN') then !}{
          double_null=.true.
        else !}{
          write(0,*) 'b2usrtrc: only works for SN or DN topology.'
          filedata(ltrc)=.false.
          return
        end if !}

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

        print '(a,5(i3,1x),l1)',
     ,   '%% b2usrtrc: j_he_at,j_h_at(1:3),j_he_i,double_null = ', !###
     ,                 j_he_at,j_h_at     ,j_he_i,double_null      !###

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
c!###{
c          print *
c          print *,'b2usrtrc: after initialisation. ',
c     ,            'nntrgrp,ndtk_tmp,nmnt_tmp,ndts_tmp=',
c     ,             nntrgrp,ndtk_tmp,nmnt_tmp,ndts_tmp
c          write(f,'(10h(1x,a9,2x,,i2,3hi3))') ndts_tmp
c          print *,'tallies (ndta_tmp):'
c          print f,(gntrgrp(i),(ndta_tmp(i,j),j=1,ndts_tmp),i=1,nntrgrp)
c          print *
c          print f,'ndtm_tmp',(ndtm_tmp(j),j=1,ndts_tmp)
c          print f,'ndtx_tmp',(ndtx_tmp(j),j=1,ndts_tmp)
c          write(f,'(13h(1x,a9,2x,1p,,i2,6he10.2))') ndts_tmp
c          print f,'hlp_awdts',(hlp_awdts(j),j=1,ndts_tmp)
c          print *,'area/absorption:'
c          write(f,'(13h(1x,a9,2x,1p,,i2,6he10.2))') ndts_tmp+1
c          print f,(gntrgrp(i),hlp_ardts(i),
c     ,                       (hlp_ksdts(i,j),j=1,ndts_tmp),i=1,nntrgrp)
c          print *
c!###}
        end if

c*** Find the Eirene surfaces corresponding to PFR interfaces
c***  (npfs_tmp, lpfs_tmp, hlp_arpfs)

        npfs_tmp=0
        lpfs_tmp=0
        hlp_arpfs=0
        if (use_eirene.ne.0) then
          do i=1,nstrai !{
            k=0
            if(crcstra(i).eq.'S') then !{
              if(rcend(i).le.leftcut(1)-1) then !{
                k=1
              else if(rcstart(i).gt.leftcut(nncut)-1) then !}{
                k=2
              end if !}
              if(k.gt.0) then !{
c                print *,'i,k,msns=',i,k,msns(1,i),msns(2,i) !###
                do j=msns(1,i),msns(1,i)+msns(2,i)-1 !{
                  npfs_tmp(k)=npfs_tmp(k)+1
                  lpfs_tmp(npfs_tmp(k),k)=lsns(j)
                  hlp_arpfs(k)=hlp_arpfs(k)+wlarea(lsns(j))
                end do !}
              end if !}
            end if !}
          end do !}
        end if

c*** Find the B2 cell ranges corresponding to PFR interfaces
c***  (ipfc_tmp, jpfc_tmp)

        hlp_w=0
        if(lpfrs_pmp.ne.0) then !{

c*** Western divertor

          m=lpfrb_i
          if(m.le.0) then !{
            if(lpfrs_pmp.eq.1) then !{
              m=leftcut(1)-1+m
            elseif(lpfrs_pmp.eq.2) then !}{
              m=nx-1+m
            elseif(lpfrs_pmp.eq.3 .or. lpfrs_pmp.eq.4) then !}{
              m=nxtl-1+m
            end if !}
          end if !}
          n=lpfrt_i
          if(n.le.0) then !{
            if(lpfrs_pmp.eq.1) then !{
              n=leftcut(1)-1+n
            elseif(lpfrs_pmp.eq.2) then !}{
              n=nx-1+n
            elseif(lpfrs_pmp.eq.3) then !}{
              n=leftcut(2)-n
            elseif(lpfrs_pmp.eq.4) then !}{
              n=nxtl-1+n
            end if !}
          end if !}
          ipfc_tmp(1)=m
          jpfc_tmp(1)=n

c*** Eastern divertor

          m=lpfrt_o
          if(m.le.0) then !{
            if(lpfrs_pmp.eq.1) then !{
              m=rightcut(1)-m
            elseif(lpfrs_pmp.eq.2) then !}{
              m=nx-1+m
            elseif(lpfrs_pmp.eq.3) then !}{
              m=rightcut(2)-1+m
            elseif(lpfrs_pmp.eq.4) then !}{
              m=nxtr-1-m
            end if !}
          end if !}
          n=lpfrb_o
          if(n.le.0) then !{
            if(lpfrs_pmp.eq.1 .or. lpfrs_pmp.eq.2) then !{
              n=nx-1+n
            else if(lpfrs_pmp.eq.3) then !}{
              n=nxtr+1-n
            else if(lpfrs_pmp.eq.4) then !}{
              n=nxtl-1+n
            end if !}
          end if !}
          ipfc_tmp(2)=m
          jpfc_tmp(2)=n

        else !}{
          ipfc_tmp(1)=1
          jpfc_tmp(1)=0
          ipfc_tmp(2)=1
          jpfc_tmp(2)=0
        end if !}

        iy=0
        if(lpfrs_pmp.eq.2 .or. lpfrs_pmp.eq.4) iy=ny
        do k=1,npfi_tmp !{
          do i=ipfc_tmp(k),jpfc_tmp(k) !{
            hlp_w(k)=hlp_w(k)+gs(i,iy,1)
          end do !}
        end do !}
c!###{
c        print *
c        print *,'b2usrtrc: PFR initialisation. npfi=',npfi_tmp
c        do k=1,npfi_tmp !{
c          print *,'PFR areaEir, areaB2 ',k,hlp_arpfs(k),hlp_w(k)
c          print *,'npfs,lpfs=',npfs_tmp(k),
c     ,                                  (lpfs_tmp(i,k),i=1,npfs_tmp(k))
c          print *,'ipfc,jpfc=',ipfc_tmp(k),jpfc_tmp(k)
c        end do !}
c        print *,'nlim,lsns: ',nlim,lsns(1:sum(ksns(1:nstrai)))
c        print *
c!###}
c-----------------------------------------------------------------------
c*** Prepare constants for DoD calculation

        isp_tmp(1)=topcut(1)
        isp_tmp(2)=isp_tmp(1)
        itg_tmp(1)=0
        itg_tmp(2)=nx-1
        imp_tmp(1)=jxa
        imp_tmp(2)=jxa
        if(double_null) then !{
          isp_tmp(3)=topcut(2)
          isp_tmp(4)=isp_tmp(3)
          itg_tmp(3)=nxtl-1
          itg_tmp(4)=nxtr+1
          imp_tmp(1)=jxi
          imp_tmp(3)=jxi
          imp_tmp(4)=jxa
        end if !}
c!###{
c        print '(10x,4i6)',1,2,3,4
c        print '(2x,a8,4i6)','isp_tmp',(isp_tmp(i),i=1,4)
c        print '(2x,a8,4i6)','itg_tmp',(itg_tmp(i),i=1,4)
c        print '(2x,a8,4i6)','imp_tmp',(imp_tmp(i),i=1,4)
c!###}
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
c!###{
c            print '(a,4i6)','ntrc_tmp,k,j,ndta_tmp(k,j):',
c     ,        ntrc_tmp,k,j,ndta_tmp(k,j)
c!###}
          end do !}
        end do !}
        if(print_isat) then !{
          ntrc_tmp=ntrc_tmp+2
          if(double_null) ntrc_tmp=ntrc_tmp+2
        end if !}

!###{
        write(6,'(3x,11a10)') 'ntrc_tmp','nglb_tmp','nspecies',
     ,    'npff_tmp','npfrgrp','nntrgrp','npfi_tmp','npft_tmp',
     ,    'ndts_tmp','ngl2_tmp','ngl4_tmp'
        write(6,'(11i10)') ntrc_tmp,nglb_tmp,nspecies,
     ,    npff_tmp,npfrgrp,nntrgrp,npfi_tmp,npft_tmp,
     ,    ndts_tmp,ngl2_tmp,ngl4_tmp
!###}
        if(.not.allocated(hlp_out)) then !{
          allocate(hlp_out(ntrc_tmp), hlp_hdr(ntrc_tmp))
c          print *,'b2usrtrc: allocated hlp_out,hlp_hdr:',ntrc_tmp !###
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
        iy=jsep+1

c        print *,'%% b2usrtrc: j_he_i,iy = ',j_he_i,iy !###

        do ix=0,nx-1 !{
          if (on_closed_surface(ix,iy-1)) then  !{
            hlp_si=hlp_si+gs(ix,0,1)
            hlp_pi=hlp_pi+gs(ix,0,1)*ne(ix,0)
            hlp=0.0_R8
            do j=0,2 !{
              hlp=hlp+na(ix,0,j_he_i+j)
              hlp_ne_he_c=hlp_ne_he_c+na(ix,0,j_he_i+j)*j*gs(ix,0,1)
            end do !}
            hlp_he=hlp_he+gs(ix,0,1)*hlp

            hlp_so=hlp_so+gs(ix,iy,1)
            hlp_po=hlp_po+gs(ix,iy,1)*ne(ix,iy)
            hlp=0.0_R8
            do j=0,2 !{
              hlp=hlp+na(ix,iy,j_he_i+j)
              hlp_ne_he_s=hlp_ne_he_s+na(ix,iy,j_he_i+j)*j*gs(ix,iy,1)
            end do !}
            hlp_hep=hlp_hep+gs(ix,iy,1)*hlp
          end if !}
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
c***  lpfrb_[io]: cell indices corresponding to the bottom (PFR) or
c***              just beginning of the region open to the duct
c***  lpfrt_[io]: cell indices corresponding to the top (PFR) or just
c***              end of the region open to the duct
c***
c*** In case of PFR, the bottom cells are counted back from the targets
c*** and the top cells down from the x-point if the corresponding value
c*** is negative.
c*** In case of SOL, negative values of lpfrb_[io] mean counting
c*** back from the next target.
c*** In case of DN configuration, negative values of lpfrt_[io] mean
c*** counting from the top targets.
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
#ifdef DBG
        print '(1x,a,6i4)',
     ,          '%%## j_he_at,nnatmi,lpfrb_i,lpfrt_i,lpfrb_o,lpfrt_o',!###
     ,                j_he_at,nnatmi,lpfrb_i,lpfrt_i,lpfrb_o,lpfrt_o  !###
#endif
        hlp_pi=0.0_R8
        hlp_si=0.0_R8
        hlp_po=0.0_R8
        hlp_so=0.0_R8
        if(lpfrs_pmp.ne.0) then !{

c*** Western divertor

          m=lpfrb_i
          if(m.le.0) then !{
            if(lpfrs_pmp.eq.1) then !{
              m=leftcut(1)-1+m
            elseif(lpfrs_pmp.eq.2) then !}{
              m=nx-1+m
            elseif(lpfrs_pmp.eq.3 .or. lpfrs_pmp.eq.4) then !}{
              m=nxtl-1+m
            end if !}
          end if !}
          n=lpfrt_i
          if(n.le.0) then !{
            if(lpfrs_pmp.eq.1) then !{
              n=leftcut(1)-1+n
            elseif(lpfrs_pmp.eq.2) then !}{
              n=nx-1+n
            elseif(lpfrs_pmp.eq.3) then !}{
              n=leftcut(2)-n
            elseif(lpfrs_pmp.eq.4) then !}{
              n=nxtl-1+n
            end if !}
          end if !}
          iy=0
          iyp=1
          if(lpfrs_pmp.eq.2 .or. lpfrs_pmp.eq.2) then
            iy=ny
            iyp=ny
          end if
          do j=m,n !{
            hlp_si=hlp_si+gs(j,iy,1)
            hlp=0.0_R8
#ifdef WG_TODO
            do k=1,nnatmi !{
              hlp=hlp+dab2(j+1,iyp,k,1)*tab2(j+1,iyp,k,1)
              if(k.eq.j_he_at) hlp_hep=hlp_hep+
     +          dab2(j+1,iyp,k,1)*tab2(j+1,iyp,k,1)*gs(j,iy,1)
            end do !}
            do k=1,nnmoli !{
              hlp=hlp+dmb2(j+1,iyp,k,1)*tmb2(j+1,iyp,k,1)
            end do !}
#endif
            hlp_pi=hlp_pi+hlp*gs(j,iy,1)
          end do !}

c          print '(a,1p,2i4,2e12.3)',' %%## m,n,hlp_si,hlp_hep: ',  !###
c     ,                                     m,n,hlp_si,hlp_hep      !###

c*** Eastern divertor

          m=lpfrt_o
          if(m.le.0) then !{
            if(lpfrs_pmp.eq.1) then !{
              m=rightcut(1)-m
            elseif(lpfrs_pmp.eq.2) then !}{
              m=nx-1+m
            elseif(lpfrs_pmp.eq.3) then !}{
              m=rightcut(2)-1+m
            elseif(lpfrs_pmp.eq.4) then !}{
              m=nxtr-1-m
            end if !}
          end if !}
          n=lpfrb_o
          if(n.le.0) then !{
            if(lpfrs_pmp.eq.1 .or. lpfrs_pmp.eq.2) then !{
              n=nx-1+n
            else if(lpfrs_pmp.eq.3) then !}{
              n=nxtr+1-n
            else if(lpfrs_pmp.eq.4) then !}{
              n=nxtl-1+n
            end if !}
          end if !}
          iy=0
          iyp=1
          if(lpfrs_pmp.eq.2 .or. lpfrs_pmp.eq.4) then !{
            iy=ny
            iyp=ny
          endif !}
          do j=m,n !{
            hlp_so=hlp_so+gs(j,iy,1)
            hlp=0.0_R8
#ifdef WG_TODO
            do k=1,nnatmi !{
              hlp=hlp+dab2(j+1,iyp,k,1)*tab2(j+1,iyp,k,1)
              if(k.eq.j_he_at) hlp_hep=hlp_hep+
     +            dab2(j+1,iyp,k,1)*tab2(j+1,iyp,k,1)*gs(j,iy,1)
            end do !}
            do k=1,nnmoli !{
              hlp=hlp+dmb2(j+1,iyp,k,1)*tmb2(j+1,iyp,k,1)
            end do !}
#endif
            hlp_po=hlp_po+hlp*gs(j,iy,1)
          end do !}

          hlp_s=hlp_si+hlp_so
          hlp_p=hlp_pi+hlp_po
        end if !}
c<<<
c        print '(a,1p,2i4,3e12.3)',' %%## m,n,hlp_so,hlp_hep,hlp_p: ',
c     ,                                   m,n,hlp_so,hlp_hep,hlp_p
c>>>
        if(hlp_s.gt.0.0_R8) then !{
          hlp_ppfr=hlp_p/hlp_s
          hlp_hepr=hlp_hep/hlp_s
        end if !}
        if(hlp_si.gt.0.0_R8) hlp_ppfi=hlp_pi/hlp_si
        if(hlp_so.gt.0.0_R8) hlp_ppfo=hlp_po/hlp_so

c*** Calculate throughputs, pumping speed, and re-scale He if necessary

        do l=1,4 !{
          if (l.eq.4) then !{
            ls=j_he_at
          else !}{
            ls=j_h_at(l)
          end if !}
          if(ls.gt.0 .and. ls.le.nnatmi) then !{

c*** Neutral flux across the side walls

            hlpnw=0.0_R8
            do js=1,nnlimi !{
              hlpnw=hlpnw+wldna(js,ls,0)-wldra(js,ls,0)-wldpa(js,ls,0)
              do j=1,nnmoli !{
                if(mlcmp(ls,j).gt.0) hlpnw=hlpnw+
     +           mlcmp(ls,j)*(wldnm(js,j,0)-wldrm(js,j,0)-wldpm(js,j,0))
              end do !}
            end do !}
          end if !}
          if (l.eq.2 .or. l.eq.3) then !{
            equg_tmp(21)=equg_tmp(21)+hlpnw/2.7e20_R8
          else if (l.eq.4) then !}{
            equg_tmp(22)=hlpnw/2.7e20_R8
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

c       print *,'%%## j_he_at,nhetrgtx,fusion_power =',  !###
c     ,               j_he_at,nhetrgtx,fusion_power      !###
        if(fusion_power.gt.0.0_R8) then !{

c*** Re-scale the helium density upstream to meet both throughput and
c*** pumping speed

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
c***  (assuming a peaking factor of hlp_neu_peak for neutrals)

      if(l_neutflux.lt.0 .or. use_eirene.eq.0) then !{ VK
c West target
        hlp_ppwi=0.0_R8
        hlp_s=0.0_R8
        do iy=0,ny-1 !{
          s=gs(0,iy,0)
          hlp=abs(fhe(0,iy,0,0)+fhi(0,iy,0,0))
          do is=0,ns-1 !{
            hlp=hlp+abs(fna(0,iy,0,0,is)*(
     .       (rpt(0,iy,is)*hx(-1,iy)+rpt(-1,iy,is)*hx(0,iy))/
     .       (hx(-1,iy)+hx(0,iy))*ev+(1.0_R8-BoRiS)*
     .       (kinrgy(0,iy,is)*hx(-1,iy)+kinrgy(-1,iy,is)*hx(0,iy))/
     .       (hx(-1,iy)+hx(0,iy))))
          end do !}
          do is=0,ns_ext-1 !{
            hlp=hlp+abs(fa_ext(0,iy,0,0,is)*(
     .       (pt_ext(0,iy,is)*hx(-1,iy)+pt_ext(-1,iy,is)*hx(0,iy))/
     .       (hx(-1,iy)+hx(0,iy))*ev+
     .       (1.0_R8-BoRiS)*0.5_R8*am_ext(is)*mp*
     .       (ua_ext(0,iy,is)**2*hx(-1,iy)+
     .        ua_ext(-1,iy,is)**2*hx(0,iy))/
     .       (hx(-1,iy)+hx(0,iy))))
          end do !}
          hlp=hlp/s
          hlp_ppwi=max(hlp_ppwi,hlp)
          hlp_s=hlp_s+s
        end do !}
        js=ltns(1)+nnlimi
        hlp_ppwi=hlp_ppwi+(wldnek(js,0)+wldnep(js,0))*hlp_neu_peak/hlp_s
c East target
        hlp_ppwo=0.0_R8
        hlp_s=0.0_R8
        do iy=0,ny-1 !{
          hlp=abs(fhe(nx,iy,0,0)+fhi(nx,iy,0,0))
          s=gs(nx,iy,0)
          do is=0,ns-1 !{
            hlp=hlp+abs(fna(nx,iy,0,0,is)*(
     .       (rpt(nx,iy,is)*hx(nx-1,iy)+rpt(nx-1,iy,is)*hx(nx,iy))/
     .       (hx(nx-1,iy)+hx(nx,iy))*ev+(1.0_R8-BoRiS)*
     .       (kinrgy(nx,iy,is)*hx(nx-1,iy)+
     .        kinrgy(nx-1,iy,is)*hx(nx,iy))/
     .       (hx(nx-1,iy)+hx(nx,iy))))
          end do !}
          do is=0,ns_ext-1 !{
            hlp=hlp+abs(fa_ext(nx,iy,0,0,is)*(
     .       (pt_ext(nx,iy,is)*hx(nx-1,iy)+
     .        pt_ext(nx-1,iy,is)*hx(nx,iy))/
     .       (hx(nx-1,iy)+hx(nx,iy))*ev+
     .       (1.0_R8-BoRiS)*0.5_R8*am_ext(is)*mp*
     .       (ua(nx,iy,is)**2*hx(nx-1,iy)+
     .        ua(nx-1,iy,is)**2*hx(nx,iy))/
     .       (hx(nx-1,iy)+hx(nx,iy))))
          end do !}
          hlp=hlp/s
          hlp_ppwo=max(hlp_ppwo,hlp)
          hlp_s=hlp_s+s
        end do !}
        js=ltns(nncut)+nnlimi
        hlp=(wldnek(js,0)+wldnep(js,0))*hlp_neu_peak/hlp_s
        hlp_ppwo=hlp_ppwo+(wldnek(js,0)+wldnep(js,0))*hlp_neu_peak/hlp_s

        if(double_null) then !{
          hlp_pwou=0.0_R8
          ix=nxtr+1
          hlp_s=0.0_R8
          do iy=0,ny-1 !{
            s=gs(ix,iy,0)
            hlp=abs(fhe(ix,iy,0,0)+fhi(ix,iy,0,0))
            do is=0,ns-1 !{
              hlp=hlp+abs(fna(ix,iy,0,0,is)*(
     .         (rpt(ix,iy,is)*hx(ix-1,iy)+
     .          rpt(ix-1,iy,is)*hx(ix,iy))/
     .         (hx(ix-1,iy)+hx(ix,iy))*ev+(1.0_R8-BoRiS)*
     .         (kinrgy(ix,iy,is)*hx(ix-1,iy)+
     .          kinrgy(ix-1,iy,is)*hx(ix,iy))/
     .         (hx(ix-1,iy)+hx(ix,iy))))
            end do !}
            do is=0,ns_ext-1 !{
              hlp=hlp+abs(fa_ext(ix,iy,0,0,is)*(
     .         (pt_ext(ix,iy,is)*hx(ix-1,iy)+
     .          pt_ext(ix-1,iy,is)*hx(ix,iy))/
     .         (hx(ix-1,iy)+hx(ix,iy))*ev+(1.0_R8-BoRiS)*
     .         (0.5_R8*am_ext(is)*mp)*
     .         (ua(ix,iy,is)**2*hx(ix-1,iy)+
     .          ua(ix-1,iy,is)**2*hx(ix,iy))/
     .         (hx(ix-1,iy)+hx(ix,iy))))
            end do !}
            hlp=hlp/s
            hlp_pwou=max(hlp_pwou,hlp)
            hlp_s=hlp_s+s
          end do !}
          js=ltns(3)+nnlimi
          hlp_pwou=hlp_pwou+(wldnek(js,0)+
     +                                 wldnep(js,0))*hlp_neu_peak/hlp_s

          hlp_pwiu=0.0_R8
          ix=nxtl
          hlp_s=0.0_R8
          do iy=1,ny !{
            hlp=abs(fhe(ix,iy,0,0)+fhi(ix,iy,0,0))
            s=gs(ix,iy,0)
            do is=0,ns-1 !{
              hlp=hlp+abs(fna(ix,iy,0,0,is)*(
     .         (rpt(ix,iy,is)*hx(ix-1,iy)+
     .          rpt(ix-1,iy,is)*hx(ix,iy))/
     .         (hx(ix-1,iy)+hx(ix,iy))*ev+(1.0_R8-BoRiS)*
     .         (kinrgy(ix,iy,is)*hx(ix-1,iy)+
     .          kinrgy(ix-1,iy,is)*hx(ix,iy))/
     .         (hx(ix-1,iy)+hx(ix,iy))))
            end do !}
            do is=0,ns_ext-1 !{
              hlp=hlp+abs(fa_ext(ix,iy,0,0,is)*(
     .         (pt_ext(ix,iy,is)*hx(ix-1,iy)+
     .          pt_ext(ix-1,iy,is)*hx(ix,iy))/
     .         (hx(ix-1,iy)+hx(ix,iy))*ev+(1.0_R8-BoRiS)*
     .         (0.5_R8*am_ext(is)*mp)*
     .         (ua(ix,iy,is)**2*hx(ix-1,iy)+
     .          ua(ix-1,iy,is)**2*hx(ix,iy))/
     .         (hx(ix-1,iy)+hx(ix,iy))))
            end do !}
            hlp=hlp/s
            hlp_pwiu=max(hlp_pwiu,hlp)
            hlp_s=hlp_s+s
          end do !}
          js=ltns(2)+nnlimi
          hlp=(wldnek(js,0)+wldnep(js,0))*hlp_neu_peak/hlp_s
          hlp_pwiu=hlp_pwiu+(wldnek(js,0)+
     +                                 wldnep(js,0))*hlp_neu_peak/hlp_s
        end if !}
      else !}{VK
#ifdef DBG
        write(*,*) 'In b2mod_usrtrc: nxtl, nxtr = ', nxtl, nxtr
#endif
        call b2usr_loads_nodiff_nodiff(nx,ny,ns,nxtl,nxtr,BoRiS,
     &    wtarg_max)
#ifdef DBG
        write(*,*) 'b2usr_loads finished'
#endif
        hlp_ppwi=wtarg_max(1)
csw 02jan2017 corrected hlp_ppwo
c       hlp_ppwo=wtarg_max(nncut)
        hlp_ppwo=wtarg_max(nncut*2)
        if(double_null) then
          hlp_pwou=wtarg_max(2)
          hlp_pwiu=wtarg_max(3)
        end if
      end if !}

#ifdef B25_EIRENE
c*** Calculate the neutral flows in PFR
      if (use_eirene.eq.0) goto 30
      do k=1,npfrgrp !{
        do j=1,nspecies !{
          hlp=1.0_R8/2.7e20_R8
          if(species_list(j).eq.'H' .or. species_list(j).eq.'D' .or.
     ,       species_list(j).eq.'T') hlp=0.5_R8*hlp
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
          hlp_dndts(k,j)=3.6e-4_R8*g*sqrt(a/(e*s*s))
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

      iy=1
      is=0
      if(lpfrs_pmp.eq.2 .or. lpfrs_pmp.eq.4) then !{
        iy=ny
        is=ny
      end if !}
      n=nnatmi+nnmoli
      do k=1,npfi_tmp !{
        do j=1,n !{
c*** here s accumulates the area, d the density*area,
c*** e the temperature*density*area, v the flux*area
          s=0.0_R8
          d=0.0_R8
          e=0.0_R8
          v=0.0_R8
          do i=ipfc_tmp(k),jpfc_tmp(k) !{
            s=s+gs(i,is,1)
#ifdef WG_TODO
            if(ndtm_tmp(j).gt.0) then !{
              g=dab2(i,iy,ndtm_tmp(j),1)
              t=tab2(i,iy,ndtm_tmp(j),1)/ev
            else !}{
              g=dmb2(i,iy,-ndtm_tmp(j),1)
              t=tmb2(i,iy,-ndtm_tmp(j),1)/ev
            end if !}
            d=d+gs(i,is,1)*g
            e=e+gs(i,is,1)*g*t
            v=v+gs(i,is,1)*g*sqrt(max(t,0.0_R8))
#endif
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
c!###{
c      d=0.0_R8
c      e=0.0_R8
c      do j=1,n !{
c        d=d+hlp_dnpfs(1,j)*hlp_tmpfs(1,j)
c        e=e+hlp_dnpfs(2,j)*hlp_tmpfs(2,j)
c      end do !}
c      print *,'prfi,prfo=',ev*d,ev*e
c!###}
#endif
c-----------------------------------------------------------------------
c*** evaluate "Degree of detachment" and Isat for each divertor

      hlp_dod=0.0_R8
      hlp_isat=0.0_R8
      mgl3_tmp=0
      call b2xppz_st_nodiff_nodiff(nx,ny,ns,ne,na,te,ti,pz)
c      print *,'### print_isat=',print_isat    !###
      do k=1,4 !{
        if(itg_tmp(k).ne.-2) then !{
c          print * !###
          d=1.e18_R8
          do iy=isp_tmp(k),ny-1 !{
            e=2.0_R8*pz(itg_tmp(k),iy)
            if(e.gt.0.0_R8) d=min(d,pz(imp_tmp(k),iy)/e)
c            print *,'e,d=',e,d  !###
          end do !}
          hlp_dod(k)=d
c          print *, k,hlp_dod(k) !###
c          print * !###
          if(print_isat) then !{
            mgl3_tmp=mgl3_tmp+1
csw 03jan2017 corrected do loop to avoid exceed of array limits and to
csw           exclude neutrals from isat calculation
            do is=0,ns-1 !{
              if(.not.is_neutral(is)) then !{
                do iy=0,ny-1 !{
                  hlp_isat(mgl3_tmp)=hlp_isat(mgl3_tmp)+
     +                      fna(itg_tmp(k)+kgl3_tmp(k)+1,iy,0,0,is)*
     *                      rza(itg_tmp(k),iy,is)
                end do !}
              end if !}
            end do !}
            hlp_isat(mgl3_tmp)=hlp_isat(mgl3_tmp)*
     *                      qe*(2.0_R8*kgl3_tmp(k)+1)
c            print *,'### mgl3_tmp=',mgl3_tmp  !###
          end if !}
        end if !}
      end do !}
c-----------------------------------------------------------------------
c*  Find the maximum ballooning factor and its location

c$$$      m=topcut(1)
c$$$      n=ny
c$$$      k=0
c$$$      if(nncut.eq.2) then !{
c$$$        if(topcut(2).gt.topcut(1)) then !{
c$$$          n=topcut(2)-1
c$$$        else if(topcut(2).eq.topcut(1)) then !}{
c$$$          k=leftcut(2)-1
c$$$        else !}{
c$$$          n=m-1
c$$$          m=topcut(2)
c$$$        end if !}
c$$$      end if !}
c$$$c      print *,'b2usrtrc: m,n,k=',m,n,k  !###
c$$$      if(k.eq.0) then !{
c$$$        al_bll_mx=maxval(aalbll(leftcut(1):leftcut(nncut),m:n))
c$$$        kgl4_tmp=maxloc(aalbll(leftcut(1):leftcut(nncut),m:n))
c$$$        al_bll_ix=kgl4_tmp(1)+leftcut(1)-1
c$$$        al_bll_iy=kgl4_tmp(2)+m-1
c$$$      else !}{
c$$$        al_bll_mx=maxval(aalbll(rightcut(2):leftcut(nncut),m:n))
c$$$        kgl4_tmp=maxloc(aalbll(rightcut(2):leftcut(nncut),m:n))
c$$$        al_bll_ix=kgl4_tmp(1)+rightcut(1)-1
c$$$        al_bll_iy=kgl4_tmp(2)+m-1
c$$$        t=maxval(aalbll(leftcut(1):rightcut(1),m:n))
c$$$        if(t.gt.al_bll_mx) then !{
c$$$          al_bll_mx=t
c$$$          kgl4_tmp=maxloc(aalbll(leftcut(1):rightcut(1),m:n))
c$$$          al_bll_ix=kgl4_tmp(1)+leftcut(1)-1
c$$$          al_bll_iy=kgl4_tmp(2)+m-1
c$$$        end if !}
c$$$      end if !}
c!###{
c      i=al_bll_ix
c      j=al_bll_iy
c      print '(2a,1p,4i5,2e12.3)','al_bll_ix,al_bll_iy,',
c     ,  'kgl4_tmp,al_bll_mx,aalbll',i,j,
c     ,  kgl4_tmp,al_bll_mx,aalbll(i,j)
c!###}
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
              else !}{
                hlp_m=adjustl(eirtxt(nnatmi-ndtm_tmp(j)))
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
c        print *,'###b ex,mgl3_tmp,n,ntrc_tmp=',ex,mgl3_tmp,n,ntrc_tmp !###
        if(ex .and. mgl3_tmp.gt.0) then !{
          i=0
          do k=1,4 !{
            if(itg_tmp(k).ne.-2) then !{
              n=n+1
              i=i+1
              hlp_hdr(n)=hlp_hg3(k)
            end if !}
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

c        print *,'n,hlp_hdr= ',n,('  ',trim(hlp_hdr(i)),i=1,n) !###

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

c        print *,'b2usrtrc: after combfile. irt,m,f=',          !###
c     ,                                     irt,m,'  ',trim(f)  !###
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
#ifdef DBG
        print *,'hlp_frm=',trim(hlp_frm)  !###
#endif
        first_out=.false.
      end if !}
c=======================================================================

      if(ank_tracing.eq.0) return
      if(filedata(ltrc)) then !{

c*** Write the data in a separate file for plotting...

        if(npfrgrp.gt.0) then !{
          i0=0
          if(npfrgrp.le.1) i0=1
          write(lun,hlp_frm,err=94) tim,(equg_tmp(j),j=1,ng_tmp),
     ,      (((equp_tmp(i,j,k),i=i0,npfrgrp),j=1,nspecies),k=1,npff_tmp)
        else !}{
          write(lun,hlp_frm,err=94) tim,(equg_tmp(j),j=1,ng_tmp)
        end if !}
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

   98 write(*,'(1x,2a,2i8)') 'b2usrtrc: internal error. ',
     ,  'Number of columns in the header - ',ntrc_tmp,n
      filedata(ltrc)=.false.
      return
      end subroutine b2usrtrc
c=======================================================================

      subroutine b2usrtri(nx,ny,ns,switch)
      use b2mod_switches_diffv_diffv
      implicit none
      integer, intent(in) :: nx, ny, ns
      integer i, j, k, l, i0, nv, na_tmp, np_tmp, ns_tmp
      type (switches), intent (in) :: switch
      character hlp_c, f*(lenfnm), hlp_hb*8
      integer lnh_tmp
      parameter (lnh_tmp=20)
      character hlp_ha(0:npfrgrpd,ns,npff_tmp)*(lnh_tmp)
      logical ex, warn_tmp
c-----------------------------------------------------------------------
c*** Initialise the file for tracing data
c-----------------------------------------------------------------------

      tab=char(9)
      call get_jsep(nx,ny,jxi,jxa,jsep)

c*** lun   :  logical unit number for the output file (=ntrc(ltrc))
c*** f     :  the associated file name

      lun=ntrc(ltrc)
      f=fn_ntrc(ltrc)
      i0=0                                                           !srv 05.11.10
      if(npfrgrp.le.1) i0=1                                          !srv 05.11.10

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
     ,   '(''# time'',',ng_tmp+np_tmp*(na_tmp+1-i0)*ns_tmp,'(2x,a))'
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
     - '(1p,e16.7,',ng_tmp+np_tmp*(na_tmp+1-i0)*ns_tmp,'(2x,e12.4))'
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
      subroutine b2usrtrcc(nx,nxtl,nxtr,iline)
      use b2mod_indirect_diffv_diffv
      implicit none
      integer, intent(in) :: nx, nxtl, nxtr
      integer, intent(out) :: iline
      integer :: i, j, k, l, m, n, i0, irt, l_pad_str, mgl3_tmp
      real (kind=R8) :: tdummy, x
      character hlp_c, hlp_tpl*2, hlp_hb*8, f*(lenfnm), hlp_m*16
      logical ex
c----------------------------------------------------------------------
c*** - to be called by b2co for locating truncation for savefile
c----------------------------------------------------------------------
c
      if(.not.filedata(ntrc_usr)) return
      hlp_tpl='UD'
      if(nncut.eq.1) then !{
c*** Single null
        hlp_tpl='SN'
      else if(nncut.eq.2) then !}{
c*** Double null
        hlp_tpl='DN'
      else if(nncut.eq.0) then !}{
c*** Limiter case
        hlp_tpl='LM'
      end if !}
      if(hlp_tpl.eq.'SN') then !{
         double_null=.false.
      else if(hlp_tpl.eq.'DN') then !}{
        double_null=.true.
      else !}{
        write(0,*) 'b2usrtrc: only works for SN or DN topology.'
        filedata(ltrc)=.false.
        return
      end if !}
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
        ntrc_tmp=ntrc_tmp+2
        if(double_null) ntrc_tmp=ntrc_tmp+2
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
      itg_tmp(1)=0
      itg_tmp(2)=nx-1
      if(double_null) then !{
        itg_tmp(3)=nxtl-1
        itg_tmp(4)=nxtr+1
      end if !}
      do k=1,4 !{
        if(itg_tmp(k).ne.-2) then !{
          if(print_isat) mgl3_tmp=mgl3_tmp+1
        end if !}
      end do !}
      ex=ex .and. n+mgl3_tmp.le.ntrc_tmp
      if(ex .and. mgl3_tmp.gt.0) then !{
        i=0
        do k=1,4 !{
          if(itg_tmp(k).ne.-2) then !{
            n=n+1
            i=i+1
            hlp_hdr(n)=hlp_hg3(k)
          end if !}
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
