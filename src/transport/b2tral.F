*-----------------------------------------------------------------------
*.specification

      subroutine b2tral (nCv, nFc, nVx, ns, nscx, nscxmax, iscx, ismain,
     1 switch, geo, mpg, pl, dv, rt, st_ext, co, co_ns)
      use b2mod_types
      use b2mod_constants
      use b2mod_b2cmpa
      use b2mod_b2zhco
     & , only: zhflcorr, zhcscorr, zhcscorr_vq,                         !som 21.03.21
     &         zh_imp_analyt, zh_analyt_mdf,
     &         zhflcorrh, zh_dist_fun, zh_tf_toff
      use b2mod_zhfrtf
     & , only: avm_u, w_out, hTdp_out, hwdp_out,
     &         rTdpgamma_out, rwdpgamma_out
!WG_TBD      use b2mod_neoclassical
      use b2mod_switches
      use b2us_geo
      use b2us_map
      use b2us_plasma
      use b2mod_input_profile
     & , only : alloc_input_profile
      use b2mod_ad
     & , only : b2news_solving, ncall_b2tral
#ifndef NAGFOR
! csc The following are not necessary for computation but are needed
!     for adjoint AD to avoid side-effect variables
      use b2mod_b2cmpt
     &  , only : cfvla, cflim, cfvsa, cfalf, cfdpa, cfsig, cfdna, 
     &           cfhce, cfhci
      use b2mod_transport_namelist
     & , only : parm_hce, parm_hci, parm_dna, parm_dpa, parm_vla,
     &          parm_vsa, parm_alf, parm_sig
      use b2mod_input_profile
     & , only : ndata, no_pflux, tr_ip_new_files, region_flags, tdata,
     &          tr_ip_elm_count, no_div, addspec,
     &          nss, nrr, nkind_coeff, nkind_data, poloidal_scaling,
     &          transport_ip_filename, nscale, transport_ip_time_switch,
     &          transport_ip_time_mod
      use b2mod_ad
     & , only : my_out_folder,
#ifdef WG_DEBUG
     &          my_out_sf,                                              !srv 27.06.21
#endif
     &          ncall_b2tlmv, ncall_b2tqin, ncall_b2tqce, ncall_b2tlnl,
     &          ncall_b2tqca, ncall_b2tlh0, ncall_b2tqna, ncall_b2trno, 
     &          ncall_transp_keps, nsdmax, cvregmax, ncall_b2trcl,
     &          ncall_b2ttia
      use b2mod_math
     & , only : cutlo, cutll, b2mod_math_initialised, small_r4_constant
#endif
      use b2mod_subsys
      use b2mod_ipmain
      implicit none
*   ..input arguments (unchanged on exit)
      integer, intent (in) :: nCv, nFc, nVx, ns, nscx, nscxmax,
     &                        iscx(0:nscxmax-1), ismain
      type (switches), intent (inout) :: switch
      type (geometry), intent (in) :: geo
      type (mapping), intent (in) :: mpg
      type (B2Plasma), intent (in) :: pl
      type (B2Derivatives), intent (inout) :: dv
      type (B2Rates), intent (in) :: rt
      type (B2StateExt), intent (in) :: st_ext
*   ..input/output arguments
      type (B2Coeff), intent (inout) :: co
      type (B2Coeff_nspecies), intent (inout) :: co_ns
*-----------------------------------------------------------------------
*.documentation
*
*  1. purpose
*
*     B2TRAL computes the transport coefficients, including geometric
*     factors such as vol/hx**2 and vol/hy**2.
*
*
*  3. description (see also routine b2cdca)
*
*     This routine computes the following discrete transport
*     coefficients.
*     - cdna: particle transport proportional to grad na.
*     - cdpa: particle transport proportional to grad pa.
*     - cvla: anomalous velocity.
*     - cvsa: viscosity.
*     - chce: electron thermal conductivity.
*     - chve: electron thermal strange velocity.
*     - chci: atom thermal conductivity.
*     - chvi: atom thermal strange velocity.
*     - csig: electrical conductivity.
*     - csigin: ion-neutral conductivity.
*     - calf: thermo-electric coefficient.
*
*     The routine operates as a simple interface to routines b2trcl
*     and b2trno, for the classical and the anomalous contributions,
*     respectively.
*
*
*-----------------------------------------------------------------------
*.declarations

*   ..local variables
      integer iFc, is, is0
      character(3) chns, chns0                                          !srv 17.06.08
      integer, save :: nnucl, zmaxmax                                   !som 09.02.21
      real (kind=R8), save :: gamma = 1.0_R8                            !som 20.03.21
      real (kind=R8) ::
     1  chvi0(nFc,0:1), csig0(nFc,0:1), calf0(nFc,0:1), chci0(nFc,0:1),
     3  cvsa0(nFc,0:1,0:ns-1), cvsahz0(nFc,0:1,0:ns-1),                 !srv 09.01.01
     4  chce0(nFc,0:1), chve0(nFc,0:1), cdkt0(nFc,0:1),
     &  chci0_fl(nFc,0:1), chci0_nofl(nFc,0:1), chci0_luc(nFc,0:1),     !som 16.06.21
     &  cvsahz0_fl(nFc,0:1,0:ns-1), cvsahz0_nofl(nFc,0:1,0:ns-1),       !som 16.06.21
     &  cvsahz0_luc(nFc,0:1,0:ns-1),
     &  cvsa0_luc(nFc,0:1,0:ns-1),
     &  chci0_al(nFc,0:1,0:nspecies-1),                                 !som 08.10.21
     &  chci0_a(nFc,0:1,0:ns-1),                                        !som 08.10.21
     &  cvsa0_hAdp(nFc,0:1,0:nspecies-1,0:nspecies-1),                  !som 20.10.21
     &  cvsa0_hBdp(nFc,0:1,0:nspecies-1),                               !som 20.10.21
     &  cvsahz0_hAdp(nFc,0:1,0:nspecies-1,0:nspecies-1),                !som 08.10.21
     &  cvsahz0_hBdp(nFc,0:1,0:nspecies-1),                             !som 08.10.21
     &  vsaf0(nFc,0:1,0:ns-1),
     &  vsaf0_luc(nFc,0:1,0:ns-1),
     1  ft,                                                             !srv 16.10.17
     *  wrkc(nCv), wrkf(nFc), wrkmin, wght(nFc,2)
        
      character*72 message
*   ..procedures
      external xertst, ipgetr, b2saxpy
      external b2trcl, b2trno, b2trql
*   ..initialisation

*-----------------------------------------------------------------------
*.computation

* ..preliminaries
*   ..subprogram start-up calls
      call subini ('b2tral')
*   ..test nCv, nFc, ns
      call xertst (0.lt.nCv.and.0.lt.nFc, 'faulty argument nCv, nFc')
      call xertst (1.le.ns, 'faulty argument ns')
*   ..read data on first call
      if (ncall_b2tral.eq.0) then
        b2news_solving(1:4)=switch%no_solve.eq.0
        do is=1,4
          b2news_solving(is)=b2news_solving(is).or.
     &    (switch%no_solve.lt.0.and.
     &     mod(abs(switch%no_solve),2**is).ge.2**(is-1))
        enddo
        call alloc_input_profile(nCv,ns)
      end if

* ..compute transport coefficients
*   ..compute flux limit coefficients
      call b2trql (nCv, nFc, ns, switch, geo, mpg, pl, dv,
     .  st_ext, co%chvemx, co%chvimx)
*   ..compute anomalous terms
      call b2trno (nCv, nFc, nVx, ns, nscx, nscxmax, iscx, ismain,
     &  switch, geo, mpg, pl, dv, rt, st_ext, co)
*   ..compute classical terms
      call b2trcl (nCv, nFc, nVx, ns, nscx, iscx, ismain,               !srv 23.11.10
     1  switch, geo, mpg, pl, dv, rt, co, st_ext, cvsahz0,              !srv 09.01.01 03.06.03
     3  cvsa0, chce0, chve0, chci0, chvi0, csig0, calf0, cdkt0,         !xpb
     4  co%csigin, co%cthe, co%cthi, co%vsaf_cl, co%cvsa_cl, 
     5  co%cvsahz_cl, co%fllimvisc, co%csig_cl, co%calf_cl)
*   ..compute Zhdanov closure coefficients
      if (switch%Zhdanov_closure.eq.1) then                             !som 09.02.21
        if (ncall_b2tral.eq.0) then
          call b2trnu (ns, nnucl, zmaxmax)
          call xertst (nnucl.eq.nspecies,
     &     'Error while building isonuclear sequences!')
        endif
        call b2tlnl (nCv, switch, switch%icase_ii, pl%te, pl%ti, dv%ne,
     &    dv%lnlam)                                                     !srv 01.07.09  20.09.11
        call b2trzh (nCv, nFc, nVx, ns, zmaxmax, switch,
     &    geo, mpg, pl, dv, rt, co, co_ns,
     &    co_ns%hci_al, co%hci_a,
     &    chci0_nofl, chci0_fl, chci0_luc,
     &    chci0_al, co_ns%chci_al, chci0_a, co%chci_a,
     &    vsaf0_luc, co%vsaf_drho,
     &    cvsahz0_nofl, cvsahz0_fl,cvsahz0_luc, co%cvsahz_drho,
     &    cvsa0_luc, co%cvsa_drho,
     &    cvsa0_hAdp, cvsa0_hBdp,
     &    co_ns%cvsa_hAdp_albe, co_ns%cvsa_hBdp_al,
     &    cvsahz0_hAdp, cvsahz0_hBdp,
     &    co_ns%cvsahz_hAdp_albe, co_ns%cvsahz_hBdp_al)                 !som 08.10.21

        if (switch%Zhdanov_test.eq.0 .and.
     &      switch%Zhdanov_cond.eq.1) then                              !som 10.08.21
          chci0(:,0:1)=chci0_luc(:,0:1)
        endif

        if (switch%Zhdanov_test.eq.0 .and.
     &      switch%Zhdanov_visu.eq.1) then                              !som 10.08.21
          if (switch%zh_vis_type.eq.0) then                             !som 13.01.22
            vsaf0(:,0:1,0:ns-1) = co%vsaf_drho(:,0:1,0:ns-1)
            cvsahz0(:,0:1,0:ns-1) = co%cvsahz_drho(:,0:1,0:ns-1)
            cvsa0(:,0:1,0:ns-1) = co%cvsa_drho(:,0:1,0:ns-1)
            co%vsaf_cl(:,0:1,0:ns-1) = co%vsaf_drho(:,0:1,0:ns-1)
            co%cvsahz_cl(:,0:1,0:ns-1) = co%cvsahz_drho(:,0:1,0:ns-1)
            co%cvsa_cl(:,0:1,0:ns-1) = co%cvsa_drho(:,0:1,0:ns-1)
          else if (switch%zh_vis_type.eq.1) then
            vsaf0(:,0:1,0:ns-1) = vsaf0_luc(:,0:1,0:ns-1)
            cvsahz0(:,0:1,0:ns-1) = cvsahz0_luc(:,0:1,0:ns-1)
            cvsa0(:,0:1,0:ns-1)= cvsa0_luc(:,0:1,0:ns-1)
            co%vsaf_cl(:,0:1,0:ns-1) = vsaf0_luc(:,0:1,0:ns-1)
            co%cvsahz_cl(:,0:1,0:ns-1) = cvsahz0_luc(:,0:1,0:ns-1)
            co%cvsa_cl(:,0:1,0:ns-1) = cvsa0_luc(:,0:1,0:ns-1)
          endif
        endif

        if(switch%b2npht_iout.ne.0 .or. switch%iout_b2wdat.eq.4) then
          call my_out_us(70,nFc,1,chci0_nofl(1,0),'b2trzh_chci_th')
          call my_out_us(70,nFc,1,chci0_fl(1,0),'b2trzh_fllim_chci_th')
          call my_out_us(70,nFc,1,chci0_luc(1,0),
     &     'b2trzh_luciani_fllim_chci_th')
          call my_out_us(70,nFc,1,chci0(1,0),
     &     'b2trzh_applied_chci_th')
          do is = 0, ns-1
            write (chns,'(i3.3)') is
            call my_out_us(70,nFc,1,chci0_a(1,0,is),
     &       'b2trzh_fllim_chci_a_th'//chns)
            call my_out_us(70,nFc,1,co_ns%chci_al(1,0,is),
     &       'b2trzh_luciani_fllim_chci_a_th'//chns)
            call my_out_us(70,nFc,1,cvsahz0_nofl(1,0,is),
     &       'b2trzh_cvsahz_th'//chns)
            call my_out_us(70,nFc,1,cvsahz0_fl(1,0,is),
     &       'b2trzh_fllim_cvsahz_th'//chns)
            call my_out_us(70,nFc,1,cvsahz0_luc(1,0,is),
     &       'b2trzh_luciani_fllim_cvsahz_th'//chns)
            call my_out_us(70,nFc,1,co%cvsahz_drho(1,0,is),
     &       'b2trzh_luciani_fllim_cvsahz_drho_th'//chns)
            call my_out_us(70,nFc,1,co%cvsahz_cl(1,0,is),
     &       'b2trzh_applied_cvsahz_th'//chns)
          enddo
          do is = 0, nspecies-1
            write (chns,'(i3.3)') is
            call my_out_us(70,nFc,1,chci0_al(1,0,is),
     &       'b2trzh_fllim_chci_al_th'//chns)
            call my_out_us(70,nFc,1,co_ns%chci_al(1,0,is),
     &       'b2trzh_luciani_fllim_chci_al_th'//chns)
            do is0 = 0, nspecies-1
             write (chns0,'(i3.3)') is0
              call my_out_us(70,nFc,1,cvsa0_hAdp(1,0,is,is0),
     &         'b2trzh_cvsa_hAdp_th'//chns//'_'//chns0)
              call my_out_us(70,nFc,1,co_ns%cvsa_hAdp_albe(1,0,is,is0),
     &         'b2trzh_luciani_cvsa_hAdp_th'//chns//'_'//chns0)
              call my_out_us(70,nFc,1,cvsahz0_hAdp(1,0,is,is0),
     &         'b2trzh_cvsahz_hAdp_th'//chns//'_'//chns0)
              call 
     &          my_out_us(70,nFc,1,co_ns%cvsahz_hAdp_albe(1,0,is,is0),
     &         'b2trzh_luciani_cvsahz_hAdp_th'//chns//'_'//chns0)
            enddo
            call my_out_us(70,nFc,1,cvsa0_hBdp(1,0,is),
     &       'b2trzh_cvsa_hBdp_th'//chns)
            call my_out_us(70,nFc,1,co_ns%cvsa_hBdp_al(1,0,is),
     &       'b2trzh_luciani_cvsa_hBdp_th'//chns)
            call my_out_us(70,nFc,1,cvsahz0_hBdp(1,0,is),
     &       'b2trzh_cvsahz_hBdp_th'//chns)
            call my_out_us(70,nFc,1,co_ns%cvsahz_hBdp_al(1,0,is),
     &       'b2trzh_luciani_cvsahz_hBdp_th'//chns)
          enddo
        endif

        if((switch%b2npht_iout.ne.0 .or.
     &      switch%iout_b2wdat.eq.4) .and. zh_dist_fun.eq.1) then
          call my_out_us(70,nCv,0,avm_u(1), 'b2trzh_avm_u')
          do is = 0, ns-1
            write (chns,'(i3.3)') is
            call my_out_us(70,nCv,0,w_out(1,is),
     &       'b2trzh_w'//chns)
            call my_out_us(70,nCv,0,hTdp_out(1,is),
     &       'b2trzh_hTdp'//chns)
            call my_out_us(70,nCv,0,hwdp_out(1,is),
     &       'b2trzh_hwdp'//chns)
            call my_out_us(70,nCv,0,rTdpgamma_out(1,is),
     &       'b2trzh_rTdpgamma'//chns)
            call my_out_us(70,nCv,0,rwdpgamma_out(1,is),
     &       'b2trzh_rwdpgamma'//chns)
          enddo
        endif
      else
        co%cvsahz_drho=0.0e0_R8
        co%cvsa_drho=0.0e0_R8
      endif

*   ..add to previous contributions=
      ft = 4.0_R8/3.0_R8                                          !srv 16.10.17 {
      do is = 0, ns-1
       if (.not.is_neutral(is)) then
        co%cvsa(:,0,is)   = ft*co%cvsa_cl(:,0,is)   + co%cvsa(:,0,is)
        co%cvsa(:,1,is)   =    co%cvsa_cl(:,1,is)   + co%cvsa(:,1,is)
        co%cvsahz(:,0,is) = ft*co%cvsahz_cl(:,0,is) + co%cvsahz(:,0,is)
        co%cvsahz(:,1,is) =    co%cvsahz_cl(:,1,is) + co%cvsahz(:,1,is)
       else
        co%cvsa(:,:,is)   =    co%cvsa_cl(:,:,is) + co%cvsa(:,:,is)
        co%cvsahz(:,:,is) =  co%cvsahz_cl(:,:,is) + co%cvsahz(:,:,is)
       end if
      end do                                                      !srv 16.10.17 }
      if (switch%fch_anomalous.ne.1.0_R8) then
        co%csig = co%csig * switch%fch_anomalous
        co%calf = co%calf * switch%fch_anomalous
        co%csig_an = co%csig_an * switch%fch_anomalous
        co%calf_an = co%calf_an * switch%fch_anomalous
      end if
      call b2saxpy (nFc*2, 1.0_R8, chce0, 1, co%chce, 1)
      call b2saxpy (nFc*2, 1.0_R8, chve0, 1, co%chve, 1)
      call b2saxpy (nFc*2, 1.0_R8, chci0, 1, co%chci, 1)
      call b2saxpy (nFc*2, 1.0_R8, chvi0, 1, co%chvi, 1)
      call b2saxpy (nFc*2, 1.0_R8, csig0, 1, co%csig, 1)
      call b2saxpy (nFc*2, 1.0_R8, calf0, 1, co%calf, 1)
      call b2saxpy (nFc*2, 1.0_R8, cdkt0, 1, co%cdkt, 1)
*   ..compute terms from diamagnetic parts                        !srv 13.10.06 {
      !wdk note: treatment made consistent with the computation of
      !          PSch-flows in b2tfnb and b2tfhe/i
      !wdk note: for correct treatment of flux-limiting of the PSch flows
      !          in b2tfnb, we take out the fcQalf factor here (and reapply it later).
      !          This effectively makes cdde/i normal conductances (not poloidal/radial).
      !          This renders 0 and 1 components identical at the moment.
      co%cdde(:,0) = geo%fcBzb*geo%fcS/
     &  (geo%fcQgam(:,0)*(geo%fcHc(:,1) + geo%fcHc(:,2)))/qe
      co%cdde(:,1) = geo%fcBzb*geo%fcS/
     &  (geo%fcQgam(:,0)*(geo%fcHc(:,1) + geo%fcHc(:,2)))/qe
      wght = 1.0_R8
      do is = 0, ns-1
       if(.not.is_neutral(is)) then
        wrkc = 1.0_R8/rt%rza(:,is)
        call intface(nCv,nFc,mpg%fcCv,wght,wrkc,wrkf)
        co%cddi(:,0,is) = wrkf*co%cdde(:,0)
        co%cddi(:,1,is) = wrkf*co%cdde(:,1)
       else
        co%cddi(:,:,is) = 0.0_R8
       endif 
      enddo                                                           !srv 13.10.06 }

*  ..check that we have some particle diffusion
      if (ncall_b2tral.le.3) then                                            !srv 12.03.12 {
        do is = 0, ns-1
cwdk combined check for poloidal or radial diffusion
          wrkf(:) = co%cdna(:,0,is)*geo%fcQalf(:,0) + 
     .              co%cdpa(:,0,is)*geo%fcQalf(:,0) + 
     .              co%cdna(:,1,is)*geo%fcQalf(:,1) + 
     .              co%cdpa(:,1,is)*geo%fcQalf(:,1)
          wrkmin = huge(1.0_R8)
          do iFc = 1, nFc
            wrkmin = min(wrkmin, wrkf(iFc))
          enddo
          write(message,'(a, i3)')
     .     'No particle diffusion for species ',is
          call xertst (wrkmin.gt.0.0_R8 .or. 
     .     switch%pol_anom_scale.eq.0.0_R8 .or. .not.b2news_solving(2),
     .     message)
*  ..check that we have some anomalous ion heat transport
          wrkmin = minval(co%hcib(:,is))
          write(message,'(a, i3)')
     .     'No ion heat diffusivity for species ',is
          call xertst (wrkmin.gt.0.0_R8 .or. .not.b2news_solving(4),
     .      message)
        enddo
*  ..check that we have some anomalous perpendicular currents
        wrkf(:) = co%calf_an(:,1)*geo%fcQalf(:,1) +
     .            co%csig_an(:,1)*geo%fcQalf(:,1) +
     .            co%csig(:,0)*geo%fcQalf(:,0)      ! on 'poloidal' faces, classical conduction is sufficient
        wrkmin = huge(1.0_R8)
        do iFc = 1, nFc
              wrkmin = min(wrkmin, wrkf(iFc))
        enddo
        call xertst (wrkmin.gt.0.0_R8 .or. switch%pot_eq.ne.1 .or.
     .   switch%fch_anomalous.eq.0.0_R8 .or. .not.b2news_solving(3) .or.
     .   switch%no_current.eq.1, 'No anomalous perpendicular currents')
      endif                                                           !srv 12.03.12 }

      if(switch%b2tral_iout.ne.0) then                                              !srv 17.06.08 {
       do is = 0, ns-1
        write (chns,'(i3.3)') is
        call my_out_us(70,nFc,1,co%cdna(1,0,is),'cdna_th'//chns)
        call my_out_us(70,nFc,1,co%cdna(1,1,is),'cdna_r'//chns)
        call my_out_us(70,nFc,1,co%cdpa(1,0,is),'cdpa_th'//chns)
        call my_out_us(70,nFc,1,co%cdpa(1,1,is),'cdpa_r'//chns)
        call my_out_us(70,nFc,1,co%cddi(1,0,is),'cddi_th'//chns)
        call my_out_us(70,nFc,1,co%cddi(1,1,is),'cddi_r'//chns)
        call my_out_us(70,nFc,1,co%cvla(1,0,is),'cvla_th'//chns)
        call my_out_us(70,nFc,1,co%cvla(1,1,is),'cvla_r'//chns)
        call my_out_us(70,nFc,1,co%cvsa(1,0,is),'cvsa_th'//chns)
        call my_out_us(70,nFc,1,co%cvsa(1,1,is),'cvsa_r'//chns)
        call my_out_us(70,nFc,1,co%cdnahz(1,0,is),'cdnahz_th'//chns)
        call my_out_us(70,nFc,1,co%cdnahz(1,1,is),'cdnahz_r'//chns)
        call my_out_us(70,nFc,1,co%cdpahz(1,0,is),'cdpahz_th'//chns)
        call my_out_us(70,nFc,1,co%cdpahz(1,1,is),'cdpahz_r'//chns)
        call my_out_us(70,nFc,1,co%cvlahz(1,0,is),'cvlahz_th'//chns)
        call my_out_us(70,nFc,1,co%cvlahz(1,1,is),'cvlahz_r'//chns)
        call my_out_us(70,nFc,1,co%cvsahz(1,0,is),'cvsahz_th'//chns)
        call my_out_us(70,nFc,1,co%cvsahz(1,1,is),'cvsahz_r'//chns)
        call my_out_us(70,nFc,1,co%cvmahz(1,0,is),'cvmahz_th'//chns)
        call my_out_us(70,nFc,1,co%cvmahz(1,1,is),'cvmahz_r'//chns)
        call my_out_us(70,nFc,1,co%cvsa_cl(1,0,is),'cvsa_cl_th'//chns)
        call my_out_us(70,nFc,1,co%cvsa_cl(1,1,is),'cvsa_cl_r'//chns)
        call my_out_us(70,nFc,1,co%cvsahz_cl(1,0,is),
     &    'cvsahz_cl_th'//chns)
        call my_out_us(70,nFc,1,co%cvsahz_cl(1,1,is),
     &    'cvsahz_cl_r'//chns)
        call my_out_us(70,nCv,0,co%cthe(1,is),'cthe'//chns)
        call my_out_us(70,nCv,0,co%cthi(1,is),'cthi'//chns)
        do is0=0,ns-1
         write (chns0,'(i3.3)') is0
         if(is0.ne.is .and. .not.is_neutral(is)) then 
          call my_out_us(70,nFc,1,co%csigin(1,0,is,is0),'csigin_th'//
     &        chns//'_'//chns0)
          call my_out_us(70,nFc,1,co%csigin(1,1,is,is0),'csigin_r'//
     &        chns//'_'//chns0)
         endif
        enddo
       enddo
       call my_out_us(70,nFc,1,co%chce(1,0),'chce_th')
       call my_out_us(70,nFc,1,co%chce(1,1),'chce_r')
       call my_out_us(70,nFc,1,co%chci(1,0),'chci_th')
       call my_out_us(70,nFc,1,co%chci(1,1),'chci_r')
       call my_out_us(70,nFc,1,co%chve(1,0),'chve_th')
       call my_out_us(70,nFc,1,co%chve(1,1),'chve_r')
       call my_out_us(70,nFc,1,co%chvi(1,0),'chvi_th')
       call my_out_us(70,nFc,1,co%chvi(1,1),'chvi_r')
       call my_out_us(70,nFc,1,co%cdde(1,0),'cdde_th')
       call my_out_us(70,nFc,1,co%cdde(1,1),'cdde_r')
       call my_out_us(70,nFc,1,co%csig(1,0),'csig_th')
       call my_out_us(70,nFc,1,co%csig(1,1),'csig_r')
       call my_out_us(70,nFc,1,co%calf(1,0),'calf_th')
       call my_out_us(70,nFc,1,co%calf(1,1),'calf_r')
       call my_out_us(70,nFc,1,co%cdkt(1,0),'cdkt_th')
       call my_out_us(70,nFc,1,co%cdkt(1,1),'cdkt_r')
      endif                                                          !srv 17.06.08 }

      if(switch%b2npco_iout.ne.0) then                               !srv 11.09.09 {
       do is = 0, ns-1
        write (chns,'(i3.3)') is 
        call my_out_us(70,nFc,1,co%cdna(1,0,is),'b2tral_cdna_th'//chns)
        call my_out_us(70,nFc,1,co%cdna(1,1,is),'b2tral_cdna_r'//chns)
        call my_out_us(70,nFc,1,co%cdpa(1,0,is),'b2tral_cdpa_th'//chns)
        call my_out_us(70,nFc,1,co%cdpa(1,1,is),'b2tral_cdpa_r'//chns)
        if(.not.is_neutral(is)) then
         call my_out_us(70,nFc,1,co%cddi(1,0,is),'b2tral_cddi_th'//chns)
         call my_out_us(70,nFc,1,co%cddi(1,1,is),'b2tral_cddi_r'//chns)
        endif
!        call my_out_us(70,nFc,1,co%cvla(1,0,is),'b2tral_cvlax'//chns)
!        call my_out_us(70,nFc,1,co%cvla(1,1,is),'b2tral_cvlay'//chns)
       enddo
      endif

      if(switch%b2npmo_iout.ne.0) then
       do is = 0, ns-1
        write (chns,'(i3.3)') is
        call my_out_us(70,nFc,1,co%cvsahz(1,0,is),'b2tral_cvsahz_th'//
     &              chns)
        call my_out_us(70,nFc,1,co%cvsahz(1,1,is),'b2tral_cvsahz_r'//
     &              chns)
!        call my_out_us(70,nFc,1,co%cvsa_cl(1,0,is),'b2tral_cvsa_clx'//
!     &              chns)
!        call my_out_us(70,nFc,1,co%cvsa_cl(1,1,is),'b2tral_cvsa_cly'//
!     &              chns)
!        call my_out_us(70,nFc,1,co%cvsahz_cl(1,0,is),'b2tral_cvsahz_clx'//
!     &             chns)
!        call my_out_us(70,nFc,1,co%cvsahz_cl(1,1,is),'b2tral_cvsahz_cly'//
!     &             chns)
        call my_out_us(70,nCv,0,co%cthe(1,is),'b2tral_cthe'//chns)
        call my_out_us(70,nCv,0,co%cthi(1,is),'b2tral_cthi'//chns)
       enddo
      endif
      if(switch%b2npht_iout.ne.0) then
       call my_out_us(70,nFc,1,co%chce(1,0),'b2tral_chce_th')
       call my_out_us(70,nFc,1,co%chce(1,1),'b2tral_chce_r')
       call my_out_us(70,nFc,1,co%chci(1,0),'b2tral_chci_th')
       call my_out_us(70,nFc,1,co%chci(1,1),'b2tral_chci_r')
!       call my_out_us(70,nFc,1,co%chve(1,0),'b2tral_chvex')
!       call my_out_us(70,nFc,1,co%chve(1,1),'b2tral_chvey')
!       call my_out_us(70,nFc,1,co%chvi(1,0),'b2tral_chvix')
!       call my_out_us(70,nFc,1,co%chvi(1,1),'b2tral_chviy')
      endif
      if(switch%b2npp7_iout.ne.0 .or. switch%b2nppo_iout.ne.0) then
       do is = 0, ns-1
        write (chns,'(i3.3)') is
        do is0 = 0, ns-1
         write (chns0,'(i3.3)') is0
         if(is0.ne.is .and. .not.is_neutral(is)) then
          call my_out_us(70,nFc,1,co%csigin(1,0,is,is0),
     &        'b2tral_csigin_th'//chns//'_'//chns0)
          call my_out_us(70,nFc,1,co%csigin(1,1,is,is0),
     &        'b2tral_csigin_r'//chns//'_'//chns0)
         endif
        enddo
       enddo
       call my_out_us(70,nFc,1,co%csig(1,0),'b2tral_csig_th')
       call my_out_us(70,nFc,1,co%csig(1,1),'b2tral_csig_r')
      endif                                                          !srv 11.09.09 }

* ..return
      ncall_b2tral = ncall_b2tral+1
      call subend ()
      return

*-----------------------------------------------------------------------
*.end b2tral

      end subroutine b2tral
