*-----------------------------------------------------------------------
*.specification

      subroutine b2trzh (nCv, nFc, nVx, ns, znmax, switch,
     &  geo, mpg, pl, dv, rt, co, co_ns,
     &  hci0_al_luc, hci0_a_luc,
     &  chci0_nofl, chci0, chci0_luc,
     &  chci0_al, chci0_al_luc, chci0_a, chci0_a_luc,
     &  vsaf0_luc, vsaf_drho_luc,
     &  cvsahz0_nofl, cvsahz0, cvsahz0_luc, cvsahz_drho_luc,
     &  cvsa0_luc, cvsa_drho_luc,
     &  cvsa_hAdp, cvsa_hBdp, cvsa_hAdp_luc, cvsa_hBdp_luc,
     &  cvsahz_hAdp, cvsahz_hBdp, cvsahz_hAdp_luc, cvsahz_hBdp_luc)     !som 08.10.2021
      use b2mod_types
      use b2mod_b2cmpa
      use b2mod_b2cmpt
      use b2mod_constants
      use b2mod_subsys
      use b2mod_b2mrel
      use b2mod_b2zhco
      use b2mod_zhfrtf
      use b2mod_transport_namelist
      use b2mod_switches
      use b2us_geo
      use b2us_map
      use b2us_plasma
      implicit none
*   ..input arguments (unchanged on exit)
      integer nCv, nFc, nVx, ns, znmax
      type (switches), intent (inout) :: switch
      type (geometry), intent (in) :: geo
      type (mapping), intent (in) :: mpg
      type (B2Plasma), intent (in) :: pl
      type (B2Derivatives), intent (inout) :: dv
      type (B2Rates), intent (in) :: rt
*   ..input/output arguments
      type (B2Coeff), intent (inout) :: co
      type (B2Coeff_nspecies), intent (inout) :: co_ns
*   ..output arguments (unspecified on entry)
      real (kind=R8) ::
     &  chci0(nFc,0:1), chci0_nofl(nFc,0:1),
     &  chci0_luc(nFc,0:1),
     &  chci0_al(nFc,0:1,0:nspecies-1),
     &  chci0_a(nFc,0:1,0:ns-1),
     &  hci0_al_luc(nCv,0:nspecies-1),
     &  hci0_a_luc(nCv,0:ns-1),
     &  chci0_al_luc(nFc,0:1,0:nspecies-1),
     &  chci0_a_luc(nFc,0:1,0:ns-1),
     &  cvsahz0(nFc,0:1,0:ns-1),
     &  cvsahz0_nofl(nFc,0:1,0:ns-1),
     &  vsaf0_luc(nFc,0:1,0:ns-1),
     &  vsaf_drho_luc(nFc,0:1,0:ns-1),
     &  cvsahz0_luc(nFc,0:1,0:ns-1),
     &  cvsahz_drho_luc(nFc,0:1,0:ns-1),
     &  cvsa0_luc(nFc,0:1,0:ns-1),
     &  cvsa_drho_luc(nFc,0:1,0:ns-1),
     &  cvsa_hAdp(nFc,0:1,0:nspecies-1,0:nspecies-1),
     &  cvsa_hBdp(nFc,0:1,0:nspecies-1),
     &  cvsa_hAdp_luc(nFc,0:1,0:nspecies-1,0:nspecies-1),
     &  cvsa_hBdp_luc(nFc,0:1,0:nspecies-1),
     &  cvsahz_hAdp(nFc,0:1,0:nspecies-1,0:nspecies-1),
     &  cvsahz_hBdp(nFc,0:1,0:nspecies-1),
     &  cvsahz_hAdp_luc(nFc,0:1,0:nspecies-1,0:nspecies-1),
     &  cvsahz_hBdp_luc(nFc,0:1,0:nspecies-1)

*-----------------------------------------------------------------------
*.documentation
*
*  1. purpose
*
*     B2TRZH computes the transport coefficients using the 21N-moment
*     approach, as detailed in [V. M. Zhdanov. Transport Processes in
*     Multicomponent Plasma, vol. 44. 10 2002].
*
*  3. description
*
*     This routine computes:
*     - c_R_TA: the thermal force transport coefficient
*               for an isonuclear sequence
*     - c_R_TA_nofl: the thermal force transport coefficient
*               for an isonuclear sequence without flux limiting
*     - c_R_TB: the thermal force transport coefficient correction
*               for each charge state
*     - c_R_TB_nofl: the thermal force transport coefficient correction
*               for each charge state without flux limiting
*     - c_R_wA: the friction force transport coefficient
*               for an isonuclear sequence
*     - c_R_wB: the friction force transport coefficient correction
*               for each charge state
*     - c_R_w: the total friction force transport coefficient
*     - c_hw_save: the velocity-dependent ion heat flux transport coefficient
*     - temperature gradients
*     - mass-related variables
*
*-----------------------------------------------------------------------
*.declarations

*   ..local variables
      integer ncall, mode
      integer iCv, iFc, is, inucl, iz, i, j
      real (kind=R8) :: temp, gamma
      real (kind=R8) :: wrk(nCv), fcSHz(nFc), fcVolHz(nFc,2)
      real (kind=R8) :: cfl(nCv)
      real (kind=R8) :: zn_temp, am_temp, amti(0:nspecies-1),
     &  zeffne, nz2n(0:nspecies-1), nati(0:ns-1), eta_cl_xy(0:ns-1),
     &  G2(0:nspecies-1,0:nspecies-1), G8(0:nspecies-1,0:nspecies-1),
     &  Ainv(0:2*nspecies-1,0:2*nspecies-1),
     &  Ainv_mdf(0:nspecies-1,0:2*nspecies-1),
     &  Ainv_fl(0:2*nspecies-1,0:2*nspecies-1),
     &  Ainv_fl_h(0:2*nspecies-1,0:2*nspecies-1),
     &  Ainv_mdf_fl(0:nspecies-1,0:2*nspecies-1),
     &  Cinv(0:2*nspecies-1,0:2*nspecies-1),
     &  c_rwA(0:nspecies-1,0:nspecies-1),
     &  c_R_TA_an_nofl(0:nspecies-1,0:nspecies-1),
     &  c_R_TA_an(0:nspecies-1,0:nspecies-1), c_rTA_an(0:nspecies-1),
     &  c_hTA_an(0:nspecies-1), c_hTA_an_fl(0:nspecies-1),
     &  c_pi_uA(0:nspecies-1,0:nspecies-1),
     &  c_pi_uB(0:nspecies-1),
     &  c_pi_uA_an(0:nspecies-1),
     &  c_pi_hA(0:nspecies-1,0:nspecies-1),
     &  c_pi_hB(0:nspecies-1),
     &  c_pi_hA_an(0:nspecies-1),
     &  tau_Zh(0:nspecies-1,0:nspecies-1),
     &  eta_cl(nCv,0:ns-1),
     &  eta_cl_fl(nCv,0:ns-1),
     &  eta_hAdp(nCv,0:nspecies-1,0:nspecies-1),
     &  eta_hBdp(nCv,0:nspecies-1),
     &  t0, t1, czetap, zetap(nCv), I_al(0:nspecies-1),
     &  cF_limh(0:nspecies-1), cF_limBh(0:nspecies-1),
     &  cF_lim(0:nspecies-1), cF_limB(0:nspecies-1),
     &  cF_limv(0:ns-1),
     &  wrkc(0:nspecies-1), wrkcc(0:nspecies-1,0:nspecies-1),
     &  wrkcA(0:nspecies-1,0:nspecies-1),
     &  wrkcB(0:nspecies-1,0:nspecies-1),
     &  wrkG2(0:nspecies-1,0:nspecies-1),
     &  wrkG8(0:nspecies-1,0:nspecies-1),
     &  hci0_nofl(nCv,0:1), hci0(nCv,0:1),
     &  hci0_al_nofl(nCv,0:nspecies-1),
     &  hci0_al(nCv,0:nspecies-1),
     &  hciB_al_nofl(0:nspecies-1),
     &  hciB_al(0:nspecies-1),
     &  hci0_a(nCv,0:ns-1),
     &  c_hw(0:nspecies-1,0:nspecies-1),c_rw(0:nspecies-1,0:nspecies-1),
     &  c_hwA(0:nspecies-1,0:nspecies-1),
     &  c_hwB(0:nspecies-1),c_rwB(0:nspecies-1),
     &  c_hTB(0:nspecies-1),c_rTB(0:nspecies-1),
     &  kappa_hTA(0:nspecies-1,0:nspecies-1),kappa_hTB(0:nspecies-1),
     &  kappa_rTA(0:nspecies-1,0:nspecies-1),kappa_rTB(0:nspecies-1),
     &  cvsahz_drho(nFc,0:1,0:ns-1),
     &  cvsa0(nFc,0:1,0:ns-1),
     &  cvsa_drho(nFc,0:1,0:ns-1)
*   ..initialisation
      save ncall
      data ncall/0/

*-----------------------------------------------------------------------
*.computation

*   ..subprogram start-up calls
      call subini ('b2trzh')
*   ..set internal parameters on first call
      if (ncall.eq.0) then
*   ..allocate
        call alloc_b2mod_mrel(nspecies,znmax,ns)
        call alloc_b2mod_zhco(nspecies, nCv)
        call alloc_b2mod_zhfrtf(nCv, ns)
        if (zh_dist_fun.eq.1) then
          call alloc_b2mod_zhfrtf_df(nCv, ns)
        endif

*   ..calculate the atom number and the atomic mass of each isonuclear sequence
        inucl=1
        znnucl(0)=zn(0)
        amnucl(0)=am(0)
        zn_temp=zn(0)
        am_temp=am(0)
        do is = 0, ns-1
          if (nint(zn(is)).ne.nint(zn_temp) .or.
     &        nint(am(is)).ne.nint(am_temp)) then
            inucl=inucl+1
            znnucl(inucl-1)=zn(is)
            amnucl(inucl-1)=am(is)
            zn_temp=zn(is)
            am_temp=am(is)
          endif
*   ..correspondence matrix
          nucl2s(inucl-1,nint(zamin(is)):nint(zamax(is)))=is
          s2nucl(is)=inucl
          znnucl_int(inucl-1)=nint(znnucl(inucl-1))
        enddo
*   ..calculate the N-coefficients and all mass-related variables
        call b2trcn (nspecies, amnucl, N2, N3, N4, N5, N6, N8, N9,
     &    N10, N11, N12, N13, N14, N15, N16, mu, rmu, rmd2, mum, mdm)   !som 02.08.2021

*   ..calculate the matrix of thermal force switches
        if (nspecies.gt.1) then
          i=nspecies-1
          j=i-1
          do is = (nspecies*(nspecies+1))/2-nspecies, 1, -1
            if (zh_tf_toff(is:is).eq.'1') then
              tf_ton_mat(i,j)=0.0_R8
              tf_ton_mat(j,i)=0.0_R8
            else
              tf_ton_mat(i,j)=1.0_R8
              tf_ton_mat(j,i)=1.0_R8
            endif
            j=j-1
            if (j.eq.-1) then
              i=i-1
              j=i-1
            endif
          enddo
        endif
      endif

      gamma = switch%b2tlc0_gamma

*   ..main part
*   ..calculate the average Z**2*n and density of each isonuclear sequence
      call b2mod_zhfrtf_prep(nCv, ns, pl%na, pl%ua)

*   ..calculate temperature gradients
      call b2mod_zhfrtf_tgta(nCv, nFc, nVx, ns, geo, mpg, pl%na, pl%ua,
     &  pl%ti)

*   ..compute zetap
      t1 = (4.0_R8*pi*eps0/qe)/qe
      t0 = (3.0_R8/(4.0_R8*sqrt(2.0_R8*pi)))*t1*(t1*ev*sqrt(ev))
      czetap = t0*sqrt(mp)
      zetap = czetap*(pl%ti/ev)*sqrt(pl%ti/ev)/dv%lnlam

*   ..compute viscous flux limit profile
      call cfl_tanh_profile(nCv, geo, mpg, cflim(2),
     .  cflmv_core, cflmv_sol, cflmv_tanh_a, cflmv_tanh_b, cfl)

*   ..calculate transport coefficietns
      do iCv = 1, nCv
        amti=amnucl*mp/pl%ti(iCv)
        zeffne=sum(z2n_xy(iCv,0:nspecies-1))
*    ..normalized average Z**2*n
        nz2n=z2n_xy(iCv,0:nspecies-1)/zeffne
*    ..calculate the correction for each charge state (heat flux/R-moment/friction)
        call b2trcz(nspecies, mum, mdm, nz2n,
     &    wrkc, wrkcB, c_hwB, c_rwB, c_hTB, c_rTB)
        c_R_TB_nofl(iCv,0:nspecies-1) = wrkc
        c_R_wB(iCv,0:nspecies-1,0:nspecies-1) = wrkcB
*    ..calculate the correction for each charge state (viscous stress tensor)
        call b2trcw(nspecies, mum, mdm, nz2n,
     &    c_pi_uB(0:nspecies-1), c_pi_hB(0:nspecies-1))

        if (zh_imp_analyt.eq.0) then
*    ..calculate the inverted matrix for the heat flux/R-moment
          call b2trca(nspecies, amti, N5, N6, N9,
     &      N10, N11, N12, nz2n, Ainv)
*    ..calculate the inverted matrix for the viscous stress tensor
          call b2trcv(nspecies, amti, mdm, N3, N4, N13,
     &      N14, N15, N16, nz2n, Cinv)

*    ..re-normalize part of the inverted matrix
          do j = 0, 2*nspecies-1
            Ainv_mdf(0:nspecies-1,j) =
     &        amti(0:nspecies-1)*Ainv(nspecies:2*nspecies-1,j)
          enddo

*    ..flux-limiting coefficient
          call b2trfl(iCv, Ainv, zeffne, amti, zetap(iCv), gamma,
     &      cF_lim, cF_limh, cF_limB, cF_limBh)
          do i = 0, nspecies-1
            Ainv_fl(i,0:nspecies-1) =
     &        cF_lim(i)*Ainv(i,0:nspecies-1)
            Ainv_mdf_fl(i,0:nspecies-1) =
     &        cF_lim(i)*Ainv_mdf(i,0:nspecies-1)
            Ainv_fl_h(i,0:nspecies-1) =
     &        cF_limh(i)*Ainv(i,0:nspecies-1)
          enddo
*    ..G-objects for the thermal and friction forces
          do j = 0, nspecies-1
            G2(:,j) = nz2n(:)*nz2n(j)*N2(:,j)
            G8(:,j) = nz2n(:)*nz2n(j)*N8(:,j)
          enddo

*    ..calculate the classical ion conductivity without flux limiting
          call b2trht(nspecies, iCv, geo%cvBb(iCv,0)/geo%cvBb(iCv,3),
     &      zetap(iCv), zeffne, pl%ti(iCv), Ainv, wrkc)
          hci0_al_nofl(iCv,0:nspecies-1) = wrkc
*    ..calculate the classical ion conductivity with flux limiting
          call b2trht(nspecies, iCv, geo%cvBb(iCv,0)/geo%cvBb(iCv,3),
     &      zetap(iCv), zeffne, pl%ti(iCv), Ainv_fl_h, wrkc)
          hci0_al(iCv,0:nspecies-1) = wrkc
*    ..calculate the velocity parts of heat flux and R-transport coefficient
          call b2trhw(nspecies, G2, G8, mum,
     &      Ainv, Ainv_mdf,
     &      c_hwA(0:nspecies-1,0:nspecies-1),
     &      c_rwA(0:nspecies-1,0:nspecies-1))
*    ..calculate the thermal force transport coefficient without flux limiting
          call b2trtf(nspecies, G2, G8, mum, mdm, Ainv, Ainv_mdf,
     &       wrkcc)
          c_R_TA_nofl(iCv,0:nspecies-1,0:nspecies-1) = wrkcc
*    ..calculate the thermal force transport coefficient with flux limiting
          wrkG2 = G2*tf_ton_mat
          wrkG8 = G8*tf_ton_mat
          call b2trtf(nspecies, wrkG2, wrkG8,
     &       mum, mdm, Ainv_fl, Ainv_mdf_fl,
     &       wrkcc)
          c_R_TA(iCv,0:nspecies-1,0:nspecies-1) = wrkcc

          do j = 0, nspecies-1
*    ..calculate the transport coefficient of the velocity-dependent part
*      of the viscosity stress
            c_pi_uA(0:nspecies-1,j) = -nz2n(0:nspecies-1)*nz2n(j)*
     &        Cinv(0:nspecies-1,j)/rmu(0:nspecies-1,j)
*    ..calculate the transport coefficient of the heat flux dependent part
*      of the viscosity stress
            c_pi_hA(0:nspecies-1,j) = -nz2n(0:nspecies-1)*nz2n(j)*
     &        (Cinv(0:nspecies-1,j)+
     &         3.5_R8*Cinv(0:nspecies-1,j+nspecies))
     &        /rmu(0:nspecies-1,j)
          enddo

*    ..for the distribution function
          if (zh_dist_fun.eq.1) then
            do i = 0, nspecies-1
              kappa_hTA(i,0:nspecies-1) =
     &          2.5_R8*(zetap(iCv)/(zeffne**2*mp))*
     &          Ainv(i,0:nspecies-1)*nal(iCv,0:nspecies-1)
              kappa_hTB(i) =
     &         -(zetap(iCv)/(zeffne**2*mp))*
     &          (c_hTB(i)/rmd2(i)*nal(iCv,i)/nz2n(i)**2)

              kappa_rTA(i,0:nspecies-1) =
     &          2.5_R8*(zetap(iCv)/(zeffne**2*mp))*
     &          Ainv_mdf(i,0:nspecies-1)*nal(iCv,0:nspecies-1)
              kappa_rTB(i) =
     &         -(zetap(iCv)/(zeffne**2*mp))*
     &          (c_rTB(i)/rmd2(i)*nal(iCv,i)/nz2n(i)**2)
            enddo
          endif

        else
*    ..the analytical thermal force transport coefficient
          call b2trhrtf_an(nspecies, mum, mdm,nz2n,
     &      c_hwA(0:nspecies-1,0:nspecies-1),
     &      c_rwA(0:nspecies-1,0:nspecies-1),
     &      c_hTA_an, c_rTA_an, c_R_TA_an_nofl)
*    ..calculate the classical ion conductivity without flux limiting
          call b2trht_an(nspecies, iCv, geo%cvBb(iCv,0)/geo%cvBb(iCv,3),
     &      zetap(iCv), zeffne, pl%ti(iCv), c_hTA_an, nz2n, rmd2,
     &      wrkc)
          hci0_al_nofl(iCv,0:nspecies-1) = wrkc
*    ..the thermal force transport coefficient without flux limiting
          do i = 0, nspecies-1
            temp = sum(c_R_TA_an_nofl(i,0:nspecies-1)*
     &        (nz2n(0:nspecies-1)/nz2n(i)))
            c_R_TA_nofl(iCv,i,0:nspecies-1) =
     &        -c_R_TA_an_nofl(0:nspecies-1,i)*
     &        (nz2n(i)/nz2n(0:nspecies-1))
            c_R_TA_nofl(iCv,i,i) = c_R_TA_nofl(iCv,i,i)+temp
          enddo

*    ..flux-limiting coefficient (analytical)
          call b2trfl_an(iCv, c_hTA_an,
     &      amnucl, amti, zetap(iCv), gamma,
     &      cF_lim, cF_limh, cF_limB, cF_limBh)
          do i = 0, nspecies-1
            c_R_TA_an(i,0:nspecies-1) =
     &        cF_lim(i)*c_R_TA_an_nofl(i,0:nspecies-1)
     &       *tf_ton_mat(i,0:nspecies-1)
             c_hTA_an_fl(i) = cF_limh(i)*c_hTA_an(i)
          enddo
*    ..calculate the classical ion conductivity with flux limiting
          call b2trht_an(nspecies, iCv, geo%cvBb(iCv,0)/geo%cvBb(iCv,3),
     &      zetap(iCv), zeffne, pl%ti(iCv), c_hTA_an_fl, nz2n, rmd2,
     &      wrkc)
          hci0_al(iCv,0:nspecies-1) = wrkc
*    ..the thermal force transport coefficient with flux limiting
          do i = 0, nspecies-1
            temp = sum(c_R_TA_an(i,0:nspecies-1)*
     &        (nz2n(0:nspecies-1)/nz2n(i)))
            c_R_TA(iCv,i,0:nspecies-1) =
     &        -c_R_TA_an(0:nspecies-1,i)*(nz2n(i)/nz2n(0:nspecies-1))
            c_R_TA(iCv,i,i) = c_R_TA(iCv,i,i)+temp
          enddo

*    ..the analytical viscous stress transport coefficient
          call b2trvis_an(nspecies, mum, mdm,nz2n,
     &      c_pi_uA_an, c_pi_hA_an)

          c_pi_uA = 0.0_R8
          c_pi_hA = 0.0_R8
          do j = 0, nspecies-1
*    ..calculate the transport coefficient of the velocity-dependent part
*      of the viscosity stress
            c_pi_uA(j,j) = 2.0_R8*c_pi_uA_an(j)
*    ..calculate the transport coefficient of the heat flux dependent part
*      of the viscosity stress
            c_pi_hA(j,j) = 2.0_R8*c_pi_hA_an(j)
          enddo

*    ..for the distribution function
          if (zh_dist_fun.eq.1) then
            kappa_hTA = 0.0_R8
            kappa_rTA = 0.0_R8
            do i = 0, nspecies-1
              kappa_hTA(i,i) =
     &         -(zetap(iCv)/(zeffne**2*mp))*
     &          (c_hTA_an(i)/rmd2(i)*nal(iCv,i)/nz2n(i)**2)
              kappa_hTB(i) =
     &         -(zetap(iCv)/(zeffne**2*mp))*
     &          (c_hTB(i)/rmd2(i)*nal(iCv,i)/nz2n(i)**2)

              kappa_rTA(i,i) =
     &         -(zetap(iCv)/(zeffne**2*mp))*
     &          (c_rTA_an(i)/rmd2(i)*nal(iCv,i)/nz2n(i)**2)
              kappa_rTB(i) =
     &         -(zetap(iCv)/(zeffne**2*mp))*
     &          (c_rTB(i)/rmd2(i)*nal(iCv,i)/nz2n(i)**2)
            enddo
          endif

        endif

*    ..calculate the thermal force correction for each charge state
        c_R_TB(iCv,0:nspecies-1) = cF_limB(0:nspecies-1)*
     &    c_R_TB_nofl(iCv,0:nspecies-1)

*    ..calculate the friction force transport coefficient
        wrkc = z2n_xy(iCv,0:nspecies-1)
        call b2trfr(nspecies, mum, rmu, wrkc,
     &    c_hwA(0:nspecies-1,0:nspecies-1),
     &    c_rwA(0:nspecies-1,0:nspecies-1),
     &    wrkcB, wrkcA, wrkcc)
        c_R_wA(iCv,0:nspecies-1,0:nspecies-1) = wrkcA
        c_R_wB(iCv,0:nspecies-1,0:nspecies-1) = wrkcB
        c_R_w (iCv,0:nspecies-1,0:nspecies-1) = wrkcc

*    ..calculate the classical eta coefficient including charge state correction
        eta_cl(iCv,0:ns-1) = 0.0_R8
        do i = 0, nspecies-1
          tau_Zh(i,0:nspecies-1) = rmu(i,0:nspecies-1)*((zetap(iCv)*
     &      nal(iCv,i))/(z2n_xy(iCv,i)*z2n_xy(iCv,0:nspecies-1)))

*    ..the velocity-dependent part
          temp = sum(c_pi_uA(i,0:nspecies-1)*tau_Zh(i,0:nspecies-1)*
     &      (nal(iCv,0:nspecies-1)/nal(iCv,i)))

          do iz= 1, znnucl_int(i)
            eta_cl(iCv,nucl2s(i,iz)) =
     &        pl%na(iCv,nucl2s(i,iz))*pl%ti(iCv)*
     &       (temp+zhcscorr*2.0_R8*c_pi_uB(i)*tau_Zh(i,i)*
     &       (z2n_xy(iCv,i)/(nal(iCv,i)*iz**2)-1.0_R8))
          enddo

*    ..the heat flux dependent part
          do j = 0, nspecies-1
            eta_hAdp(iCv,i,j) = c_pi_hA(i,j)*tau_Zh(i,j)
          enddo

          eta_hBdp(iCv,i) = 2.0_R8*c_pi_hB(i)*tau_Zh(i,i)
        enddo

*    ..calculate flux limiting for the viscous stress
        nati(0:ns-1) = pl%na(iCv,0:ns-1)*pl%ti(iCv)
        eta_cl_xy(0:ns-1) = eta_cl(iCv,0:ns-1)
        call b2trfl_vis (iCv, ns,
     &    eta_cl_xy, nati, gamma, cfl(iCv), cF_limv)

*    ..the velocity-dependent part
        do is = 0, ns-1
*    ..apply bx^2
          eta_cl(iCv,is) = eta_cl(iCv,is)*
     &      (geo%cvBb(iCv,0)/geo%cvBb(iCv,3))**2
*    ..apply flux limiting
          eta_cl_fl(iCv,is) = eta_cl(iCv,is)*
     &      cF_limv(is)
        enddo

*    ..the heat flux dependent part
        do i = 0, nspecies-1
*    ..apply bx^2
          do j = 0, nspecies-1
            eta_hAdp(iCv,i,j) = eta_hAdp(iCv,i,j)*
     &        (geo%cvBb(iCv,0)/geo%cvBb(iCv,3))**2
          enddo

          eta_hBdp(iCv,i) = eta_hBdp(iCv,i)*
     &      (geo%cvBb(iCv,0)/geo%cvBb(iCv,3))**2
        enddo

*   ..tilde I for each isonuclear sequence
        do inucl = 0, nspecies-1
          I_al(inucl)=0.0_R8
          do iz = 1, znnucl_int(inucl)
            I_al(inucl) = I_al(inucl)+
     &        (pl%na(iCv,nucl2s(inucl,iz))/nal(iCv,inucl))**2/
     &         Ia(iCv,nucl2s(inucl,iz))
          enddo
        enddo

*    ..calculate the classical ion conductivity correction
        call b2trht_an(nspecies, iCv, geo%cvBb(iCv,0)/geo%cvBb(iCv,3),
     &    zetap(iCv), zeffne, pl%ti(iCv), c_hTB, nz2n, rmd2,
     &    hciB_al_nofl(0:nspecies-1))
*    ..calculate the classical ion conductivity correction with flux limiting
        wrkc = cF_limBh*c_hTB
        call b2trht_an(nspecies, iCv, geo%cvBb(iCv,0)/geo%cvBb(iCv,3),
     &    zetap(iCv), zeffne, pl%ti(iCv), wrkc, nz2n, rmd2,
     &    hciB_al(0:nspecies-1))

*   ..apply charge state correction for the heat flux
        hci0_nofl(iCv,0) =
     &    sum(hci0_al_nofl(iCv,0:nspecies-1)+zhcscorr*
     &    hciB_al_nofl(0:nspecies-1)*(I_al(0:nspecies-1)-1.0_R8))
        hci0(iCv,0) =
     &    sum(hci0_al(iCv,0:nspecies-1)+zhcscorr*
     &    hciB_al(0:nspecies-1)*(I_al(0:nspecies-1)-1.0_R8))

*   ..ion heat conductivity for each charge state
        hci0_a(iCv,0:ns-1) = 0.0_R8
        do i = 0, nspecies-1
          do iz= 1, znnucl_int(i)
            hci0_a(iCv,nucl2s(i,iz)) =
     &        pl%na(iCv,nucl2s(i,iz))/nal(iCv,i)*
     &       (hci0_al(iCv,i)+zhcscorr*hciB_al(i)*
     &        (z2n_xy(iCv,i)/(nal(iCv,i)*iz**2)-1.0_R8))
          enddo
        enddo

        c_hw(0:nspecies-1,0:nspecies-1) =
     &    c_hwA(0:nspecies-1,0:nspecies-1)
        do inucl = 0, nspecies-1
          c_hw(inucl,inucl) = c_hw(inucl,inucl)+zhcscorr*
     &      c_hwB(inucl)-sum(c_hwA(0:nspecies-1,inucl))
        enddo
*   ..save heat flux velocity-dependent part transport coefficient
        c_hw_save(iCv,0:nspecies-1,0:nspecies-1) =
     &    c_hw(0:nspecies-1,0:nspecies-1)

        if (zh_dist_fun.eq.1) then
*   ..apply charge state correction for the R-moment
          c_rw(0:nspecies-1,0:nspecies-1) =
     &      c_rwA(0:nspecies-1,0:nspecies-1)
          do inucl = 0, nspecies-1
            c_rw(inucl,inucl) = c_rw(inucl,inucl)+zhcscorr*
     &        c_rwB(inucl)-sum(c_rwA(0:nspecies-1,inucl))
          enddo

*    ..moments of the distribution function
          do inucl = 0, nspecies-1
            do iz= 1, znnucl_int(inucl)
              w_out(iCv,nucl2s(inucl,iz)) =
     &          pl%ua(iCv,nucl2s(inucl,iz))-avm_u(iCv)

              hTdp_out(iCv,nucl2s(inucl,iz)) =
     &          sum(kappa_hTA(inucl,0:nspecies-1)*
     &              gtalc(iCv,0:nspecies-1))+
     &          zhcscorr*kappa_hTB(inucl)*
     &         (pl%na(iCv,nucl2s(inucl,iz))/nal(iCv,inucl)/
     &          Ia(iCv,nucl2s(inucl,iz))*
     &          gt_Ac(iCv,nucl2s(inucl,iz))-
     &          gtalc(iCv,inucl))

              hwdp_out(iCv,nucl2s(inucl,iz)) =
     &          sum(c_hw(0:nspecies-1,inucl)*
     &              (pl%ua(iCv,nucl2s(inucl,iz))-
     &               av_ualpha(iCv,0:nspecies-1)))

              rTdpgamma_out(iCv,nucl2s(inucl,iz)) =
     &          sum(kappa_rTA(inucl,0:nspecies-1)*
     &              gtalc(iCv,0:nspecies-1))+
     &          zhcscorr*kappa_rTB(inucl)*
     &         (pl%na(iCv,nucl2s(inucl,iz))/nal(iCv,inucl)/
     &          Ia(iCv,nucl2s(inucl,iz))*
     &          gt_Ac(iCv,nucl2s(inucl,iz))-
     &          gtalc(iCv,inucl))

              rwdpgamma_out(iCv,nucl2s(inucl,iz)) =
     &          sum(c_rw(0:nspecies-1,inucl)*
     &             (pl%ua(iCv,nucl2s(inucl,iz))-
     &              av_ualpha(iCv,0:nspecies-1)))
            enddo
          enddo
        endif

      enddo

      mode = switch%b2tral_mode

*    ..ion conductivity on the cell boundary
      call b2txcx (nCv, nFc, mode, geo, mpg, geo%fcVol, geo%fcS,
     &  hci0(:,0), chci0(:,0))
      call b2txcx (nCv, nFc, mode, geo, mpg, geo%fcVol, geo%fcS,
     & hci0_nofl(:,0), chci0_nofl(:,0))
      do i = 0, nspecies-1
        call b2txcx (nCv, nFc, mode, geo, mpg, geo%fcVol, geo%fcS,
     &   hci0_al(:,i), chci0_al(:,0,i))
      enddo
      do is = 0, ns-1
        call b2txcx (nCv, nFc, mode, geo, mpg, geo%fcVol, geo%fcS,
     &   hci0_a(:,is), chci0_a(:,0,is))
      enddo

*    ..ion velocity-dependent viscosity on the cell boundary
      fcSHz = geo%fcS*geo%fcHz
      do iFc = 1, nFc
        fcVolHz(iFc,1) = geo%fcVol(iFc,1)*geo%cvHz(mpg%fcCv(iFc,1))
        fcVolHz(iFc,2) = geo%fcVol(iFc,2)*geo%cvHz(mpg%fcCv(iFc,2))
      enddo
      cvsa_drho   = 0.0_R8                                              !srv 30.01.23
      cvsahz_drho = 0.0_R8                                              !srv 30.01.23
      do is = 0, ns-1
        wrk = eta_cl_fl(:,is)*rho_A_rel(:,is)
        do iFc = 1, nFc
          vsaf0_luc(iFc,0,is) =
     &      (geo%fcVol(iFc,1)*wrk(mpg%fcCv(iFc,1))+
     &       geo%fcVol(iFc,2)*wrk(mpg%fcCv(iFc,2)))/
     &      (geo%fcVol(iFc,1)+geo%fcVol(iFc,2))
        enddo
        vsaf0_luc(:,1,is) = 0.0_R8
        call b2txcx (nCv, nFc, mode, geo, mpg, fcVolHz, fcSHz,
     &   wrk, cvsahz0(:,0,is))
        call b2txcx (nCv, nFc, mode, geo, mpg, geo%fcVol, geo%fcS,
     &   wrk, cvsa0(:,0,is))
        call b2txcx (nCv, nFc, mode, geo, mpg, fcVolHz, fcSHz,
     &   eta_cl_fl(:,is), cvsahz_drho(:,0,is))
        wrk = eta_cl(:,is)*rho_A_rel(:,is)
        do iFc = 1, nFc
          vsaf_drho_luc(iFc,0,is) =
     &      (geo%fcVol(iFc,1)*wrk(mpg%fcCv(iFc,1))+
     &       geo%fcVol(iFc,2)*wrk(mpg%fcCv(iFc,2)))/
     &      (geo%fcVol(iFc,1)+geo%fcVol(iFc,2))
        enddo
        vsaf_drho_luc(:,1,is) = 0.0_R8
        call b2txcx (nCv, nFc, mode, geo, mpg, fcVolHz, fcSHz,
     &   wrk, cvsahz0_nofl(:,0,is))
        call b2txcx (nCv, nFc, mode, geo, mpg, geo%fcVol, geo%fcS,
     &   eta_cl_fl(:,is), cvsa_drho(:,0,is))
      enddo

*    ..ion heat flux dependent viscosity on the cell boundary
      do i = 0, nspecies-1
        do j = 0, nspecies-1
          call b2txcx (nCv, nFc, mode, geo, mpg, fcVolHz, fcSHz,
     &     eta_hAdp(:,i,j), cvsahz_hAdp(:,0,i,j))
          call b2txcx (nCv, nFc, mode, geo, mpg, geo%fcVol, geo%fcS,
     &     eta_hAdp(:,i,j), cvsa_hAdp(:,0,i,j))
        enddo

        call b2txcx (nCv, nFc, mode, geo, mpg, fcVolHz, fcSHz,
     &   eta_hBdp(:,i), cvsahz_hBdp(:,0,i))
        call b2txcx (nCv, nFc, mode, geo, mpg, geo%fcVol, geo%fcS,
     &   eta_hBdp(:,i), cvsa_hBdp(:,0,i))
      enddo

*    ..apply Luciani correction for the ion conductivity
      chci0_luc=0.0_R8
      chci0_luc(:,0) = chci0_luc(:,0)*co%f_luc_ki
      hci0(:,0) = hci0(:,0)*co%f_luc_ki_c
      do i = 0, nspecies-1
        chci0_al_luc(:,0,i) = chci0_al(:,0,i)*co%f_luc_ki
        hci0_al_luc(:,i) = hci0_al(:,i)*co%f_luc_ki_c
      enddo
      do is = 0, ns-1
        chci0_a_luc(:,0,is) = chci0_a(:,0,is)*co%f_luc_ki
        hci0_a_luc(:,is) = hci0_a(:,is)*co%f_luc_ki_c
      enddo

*    ..apply Luciani correction for the ion velocity-dependent viscosity
      cvsahz0_luc     = cvsahz0                                         !srv 30.01.23 {
      cvsahz_drho_luc = cvsahz_drho
      cvsa0_luc       = cvsa0
      cvsa_drho_luc   = cvsa_drho                                       !srv 30.01.23 }
      do is = 0, ns-1
        if (.not.is_neutral(is)) then
          vsaf0_luc(:,0,is) = vsaf0_luc(:,0,is)*co%f_luc_et
          vsaf_drho_luc(:,0,is) = vsaf_drho_luc(:,0,is)*co%f_luc_et
          cvsahz0_luc(:,0,is) = cvsahz0(:,0,is)*co%f_luc_et
          cvsahz_drho_luc(:,0,is) = cvsahz_drho(:,0,is)*co%f_luc_et
          cvsa0_luc(:,0,is) = cvsa0(:,0,is)*co%f_luc_et
          cvsa_drho_luc(:,0,is) = cvsa_drho(:,0,is)*co%f_luc_et
        endif
      enddo

*    ..apply Luciani and neoclassical corrections for the ion heat flux
*      dependent viscosity
      do i = 0, nspecies-1
        do j = 0, nspecies-1
          cvsahz_hAdp_luc(:,0,i,j) =
     &      cvsahz_hAdp(:,0,i,j)*co%f_luc_et*co%alpha1f
          cvsa_hAdp_luc(:,0,i,j) =
     &      cvsa_hAdp(:,0,i,j)*co%f_luc_et*co%alpha1f
        enddo
        cvsahz_hBdp_luc(:,0,i) = cvsahz_hBdp(:,0,i)*
     &    co%f_luc_et*co%alpha1f
        cvsa_hBdp_luc(:,0,i) = cvsa_hBdp(:,0,i)*
     &    co%f_luc_et*co%alpha1f
      enddo

      ncall = ncall+1
      call subend ()
      return

      end subroutine b2trzh
