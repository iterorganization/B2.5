*-----------------------------------------------------------------------
*.specification

      subroutine b2trzh (nx, ny, ns, znmax,
     &  na, ua, ti, lnlam, hx, hz, bb, vol, gs, qz,
     &  mode, gamma, chci0_nofl, chci0, chci0_luc,
     &  chci0_al, chci0_al_luc, chci0_a, chci0_a_luc,
     &  cvsahz0_nofl, cvsahz0, cvsahz0_luc, cvsahz_drho_luc,
     &  cvsa0_luc, cvsa_drho_luc,
     &  cvsa_hAdp, cvsa_hBdp, cvsa_hAdp_luc, cvsa_hBdp_luc,
     &  cvsahz_hAdp, cvsahz_hBdp, cvsahz_hAdp_luc, cvsahz_hBdp_luc)     !som 08.10.2021
      use b2mod_types
      use b2mod_b2cmpa
      use b2mod_constants
      use b2mod_subsys
      use b2mod_b2mrel
      use b2mod_b2zhco
      use b2mod_zhfrtf
      use b2mod_transport
     , , only: f_luc_ki, f_luc_et, alpha1
      use b2mod_indirect
      implicit none
*   ..input arguments (unchanged on exit)
      integer nx, ny, ns, znmax, mode
      real (kind=R8) :: gamma                             !som 20.03.21
      real (kind=R8) ::
     &  na(-1:nx,-1:ny,0:ns-1), ti(-1:nx,-1:ny), ua(-1:nx,-1:ny,0:ns-1),
     &  lnlam(-1:nx,-1:ny), hx(-1:nx,-1:ny), hz(-1:nx,-1:ny),
     &  bb(-1:nx,-1:ny,0:3), vol(-1:nx,-1:ny), gs(-1:nx,-1:ny,0:2),
     &  qz(-1:nx,-1:ny,0:1)
*   ..output arguments (unspecified on entry)
      real (kind=R8) ::
     &  chci0(-1:nx,-1:ny,0:1), chci0_nofl(-1:nx,-1:ny,0:1),
     &  chci0_luc(-1:nx,-1:ny,0:1),
     &  chci0_al(-1:nx,-1:ny,0:1,0:nspecies-1),
     &  chci0_a(-1:nx,-1:ny,0:1,0:ns-1),
     &  chci0_al_luc(-1:nx,-1:ny,0:1,0:nspecies-1),
     &  chci0_a_luc(-1:nx,-1:ny,0:1,0:ns-1),
     &  cvsahz0(-1:nx,-1:ny,0:1,0:ns-1),
     &  cvsahz_drho(-1:nx,-1:ny,0:1,0:ns-1),
     &  cvsahz0_nofl(-1:nx,-1:ny,0:1,0:ns-1),
     &  cvsahz0_luc(-1:nx,-1:ny,0:1,0:ns-1),
     &  cvsahz_drho_luc(-1:nx,-1:ny,0:1,0:ns-1),
     &  cvsa0(-1:nx,-1:ny,0:1,0:ns-1),
     &  cvsa_drho(-1:nx,-1:ny,0:1,0:ns-1),
     &  cvsa0_luc(-1:nx,-1:ny,0:1,0:ns-1),
     &  cvsa_drho_luc(-1:nx,-1:ny,0:1,0:ns-1),
     &  cvsa_hAdp(-1:nx,-1:ny,0:1,0:nspecies-1,0:nspecies-1),
     &  cvsa_hBdp(-1:nx,-1:ny,0:1,0:nspecies-1),
     &  cvsa_hAdp_luc(-1:nx,-1:ny,0:1,0:nspecies-1,0:nspecies-1),
     &  cvsa_hBdp_luc(-1:nx,-1:ny,0:1,0:nspecies-1),
     &  cvsahz_hAdp(-1:nx,-1:ny,0:1,0:nspecies-1,0:nspecies-1),
     &  cvsahz_hBdp(-1:nx,-1:ny,0:1,0:nspecies-1),
     &  cvsahz_hAdp_luc(-1:nx,-1:ny,0:1,0:nspecies-1,0:nspecies-1),
     &  cvsahz_hBdp_luc(-1:nx,-1:ny,0:1,0:nspecies-1)

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
      integer ncall
      integer ix, iy, is, inucl, iz,i, j
      real (kind=R8) :: temp
      real (kind=R8) :: zn_temp, am_temp, amti(0:nspecies-1),
     &  zeffne, nz2n(0:nspecies-1),
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
     &  eta_cl(-1:nx,-1:ny,0:ns-1),
     &  eta_cl_fl(-1:nx,-1:ny,0:ns-1),
     &  eta_hAdp(-1:nx,-1:ny,0:nspecies-1,0:nspecies-1),
     &  eta_hBdp(-1:nx,-1:ny,0:nspecies-1),
     &  t0, t1, czetap,zetap(-1:nx,-1:ny),I_al(0:nspecies-1),
     &  cF_limh(0:nspecies-1),cF_limBh(0:nspecies-1),
     &  cF_lim(0:nspecies-1),cF_limB(0:nspecies-1),
     &  cF_limv(0:ns-1),
     &  hci0_nofl(-1:nx,-1:ny,0:1),hci0(-1:nx,-1:ny,0:1),
     &  hci0_al_nofl(-1:nx,-1:ny,0:nspecies-1),
     &  hci0_al(-1:nx,-1:ny,0:nspecies-1),
     &  hciB_al_nofl(0:nspecies-1),
     &  hciB_al(0:nspecies-1),
     &  hci0_a_nofl(-1:nx,-1:ny,0:ns-1),
     &  hci0_a(-1:nx,-1:ny,0:ns-1),
     &  c_hw(0:nspecies-1,0:nspecies-1),c_rw(0:nspecies-1,0:nspecies-1),
     &  c_hwA(0:nspecies-1,0:nspecies-1),
     &  c_hwB(0:nspecies-1),c_rwB(0:nspecies-1),
     &  c_hTB(0:nspecies-1),c_rTB(0:nspecies-1),
     &  kappa_hTA(0:nspecies-1,0:nspecies-1),kappa_hTB(0:nspecies-1),
     &  kappa_rTA(0:nspecies-1,0:nspecies-1),kappa_rTB(0:nspecies-1)
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
        call alloc_b2mod_zhco(nspecies, nx, ny)
        call alloc_b2mod_zhfrtf(nx, ny, ns)
        if (zh_dist_fun.eq.1) then
          call alloc_b2mod_zhfrtf_df(nx, ny, ns)
        endif

*   ..calculate the atom number and the atomic mass of each isonuclear sequence
        inucl=1
        znnucl(0)=zn(0)
        amnucl(0)=am(0)
        zn_temp=zn(0)
        am_temp=am(0)
        do is = 0, ns-1
          if (zn(is).ne.zn_temp .and. am(is).ne.am_temp) then
            inucl=inucl+1
            znnucl(inucl-1)=zn(is)
            amnucl(inucl-1)=am(is)
            zn_temp=zn(is)
            am_temp=am(is)
          endif
*   ..correspondence matrix
          nucl2s(inucl-1,nint(zamax(is)))=is
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

*   ..main part
*   ..calculate the average Z**2*n and density of each isonuclear sequence
      call b2mod_zhfrtf_prep(nx,ny,ns,na,ua)

*   ..calculate temperature gradients
      call b2mod_zhfrtf_tgta(nx, ny, ns, na, ua, ti, hx, bb)

*   ..compute zetap
      t1 = (4.0_R8*pi*eps0/qe)/qe
      t0 = (3.0_R8/(4.0_R8*sqrt(2.0_R8*pi)))*t1*(t1*ev*sqrt(ev))
      czetap = t0*sqrt(mp)
      zetap = czetap*(ti/ev)*sqrt(ti/ev)/lnlam

*   ..clean output coefficients
      chci0_luc=0.0_R8
      cvsahz0_luc=0.0_R8
      cvsa0_luc=0.0_R8

*   ..calculate transport coefficietns
      do iy = -1, ny
        do ix = -1, nx
          amti=amnucl*mp/ti(ix,iy)
          zeffne=sum(z2n_xy(ix,iy,0:nspecies-1))
*    ..normalized average Z**2*n
          nz2n=z2n_xy(ix,iy,0:nspecies-1)/zeffne
*    ..calculate the correction for each charge state (heat flux/R-moment/friction)
          call b2trcc(nspecies, mum, mdm, nz2n,
     &     c_R_TB_nofl(ix,iy,0:nspecies-1),
     &     c_R_wB(ix,iy,0:nspecies-1,0:nspecies-1),
     &     c_hwB, c_rwB, c_hTB, c_rTB)
*    ..calculate the correction for each charge state (viscous stress tensor)
          call b2trcs(nspecies, mum, mdm, nz2n,
     &     c_pi_uB(0:nspecies-1), c_pi_hB(0:nspecies-1))

          if (zh_imp_analyt.eq.0) then
*    ..calculate the inverted matrix for the heat flux/R-moment
            call b2trca(nspecies, amti, N5, N6, N9,
     &       N10, N11, N12, nz2n, Ainv)
*    ..calculate the inverted matrix for the viscous stress tensor
            call b2trcv(nspecies, amti, mdm, N3, N4, N13,
     &       N14, N15, N16, nz2n, Cinv)

*    ..re-normalize part of the inverted matrix
            do j = 0, 2*nspecies-1
              Ainv_mdf(0:nspecies-1,j) =
     &         amti(0:nspecies-1)*Ainv(nspecies:2*nspecies-1,j)
            enddo

*    ..flux-limiting coefficient
            call b2trfl(ix, iy, Ainv, zeffne, amti, zetap(ix,iy), gamma,
     &       cF_lim, cF_limh, cF_limB, cF_limBh)
            do i = 0, nspecies-1
              Ainv_fl(i,0:nspecies-1) =
     &         cF_lim(i)*Ainv(i,0:nspecies-1)
              Ainv_mdf_fl(i,0:nspecies-1) =
     &         cF_lim(i)*Ainv_mdf(i,0:nspecies-1)
              Ainv_fl_h(i,0:nspecies-1) =
     &         cF_limh(i)*Ainv(i,0:nspecies-1)
            enddo
*    ..G-objects for the thermal and friction forces
            do j = 0, nspecies-1
              G2(:,j) = nz2n(:)*nz2n(j)*N2(:,j)
              G8(:,j) = nz2n(:)*nz2n(j)*N8(:,j)
            enddo

*    ..calculate the classical ion conductivity without flux limiting
            call b2trht(nspecies, ix, iy, bb(ix,iy,0)/bb(ix,iy,3),
     &       zetap(ix,iy), zeffne, ti(ix,iy), Ainv,
     &       hci0_al_nofl(ix,iy,0:nspecies-1))
*    ..calculate the classical ion conductivity with flux limiting
            call b2trht(nspecies, ix, iy, bb(ix,iy,0)/bb(ix,iy,3),
     &       zetap(ix,iy), zeffne, ti(ix,iy), Ainv_fl_h,
     &       hci0_al(ix,iy,0:nspecies-1))
*    ..calculate the velocity parts of heat flux and R-transport coefficient
            call b2trhw(nspecies, G2, G8, mum,
     &       Ainv, Ainv_mdf,
     &       c_hwA(0:nspecies-1,0:nspecies-1),
     &       c_rwA(0:nspecies-1,0:nspecies-1))
*    ..calculate the thermal force transport coefficient without flux limiting
            call b2trtf(nspecies, G2, G8, mum, mdm, Ainv, Ainv_mdf,
     &       c_R_TA_nofl(ix,iy,0:nspecies-1,0:nspecies-1))
*    ..calculate the thermal force transport coefficient with flux limiting
            call b2trtf(nspecies, G2*tf_ton_mat, G8*tf_ton_mat,
     &        mum, mdm, Ainv_fl, Ainv_mdf_fl,
     &        c_R_TA(ix,iy,0:nspecies-1,0:nspecies-1))

            do j = 0, nspecies-1
*    ..calculate the transport coefficient of the velocity-dependent part
*      of the viscosity stress
              c_pi_uA(0:nspecies-1,j) = -nz2n(0:nspecies-1)*nz2n(j)*
     &         Cinv(0:nspecies-1,j)/rmu(0:nspecies-1,j)
*    ..calculate the transport coefficient of the heat flux dependent part
*      of the viscosity stress
              c_pi_hA(0:nspecies-1,j) = -nz2n(0:nspecies-1)*nz2n(j)*
     &         (Cinv(0:nspecies-1,j)+
     &          3.5_R8*Cinv(0:nspecies-1,j+nspecies))
     &         /rmu(0:nspecies-1,j)
            enddo

*    ..for the distribution function
            if (zh_dist_fun.eq.1) then
              do i = 0, nspecies-1
                kappa_hTA(i,0:nspecies-1) =
     &           2.5_R8*(zetap(ix,iy)/(zeffne**2*mp))*
     &           Ainv(i,0:nspecies-1)*nal(ix,iy,0:nspecies-1)
                kappa_hTB(i) =
     &          -(zetap(ix,iy)/(zeffne**2*mp))*
     &           (c_hTB(i)/rmd2(i)*nal(ix,iy,i)/nz2n(i)**2)

                kappa_rTA(i,0:nspecies-1) =
     &           2.5_R8*(zetap(ix,iy)/(zeffne**2*mp))*
     &           Ainv_mdf(i,0:nspecies-1)*nal(ix,iy,0:nspecies-1)
                kappa_rTB(i) =
     &          -(zetap(ix,iy)/(zeffne**2*mp))*
     &           (c_rTB(i)/rmd2(i)*nal(ix,iy,i)/nz2n(i)**2)
              enddo
            endif

          else
*    ..the analytical thermal force transport coefficient
            call b2trhrtf_an(nspecies, mum, mdm,nz2n,
     &       c_hwA(0:nspecies-1,0:nspecies-1),
     &       c_rwA(0:nspecies-1,0:nspecies-1),
     &       c_hTA_an, c_rTA_an, c_R_TA_an_nofl)
*    ..calculate the classical ion conductivity without flux limiting
            call b2trht_an(nspecies, ix, iy, bb(ix,iy,0)/bb(ix,iy,3),
     &       zetap(ix,iy),zeffne,ti(ix,iy),c_hTA_an,nz2n,rmd2,
     &       hci0_al_nofl(ix,iy,0:nspecies-1))
*    ..the thermal force transport coefficient without flux limiting
            do i = 0, nspecies-1
              temp = sum(c_R_TA_an_nofl(i,0:nspecies-1)*
     &         (nz2n(0:nspecies-1)/nz2n(i)))
              c_R_TA_nofl(ix,iy,i,0:nspecies-1) =
     &         -c_R_TA_an_nofl(0:nspecies-1,i)*
     &         (nz2n(i)/nz2n(0:nspecies-1))
              c_R_TA_nofl(ix,iy,i,i) = c_R_TA_nofl(ix,iy,i,i)+temp
            enddo

*    ..flux-limiting coefficient (analytical)
            call b2trfl_an(ix, iy, c_hTA_an,
     &       amnucl, amti, zetap(ix,iy), gamma,
     &       cF_lim, cF_limh, cF_limB, cF_limBh)
            do i = 0, nspecies-1
              c_R_TA_an(i,0:nspecies-1) =
     &         cF_lim(i)*c_R_TA_an_nofl(i,0:nspecies-1)
     &        *tf_ton_mat(i,0:nspecies-1)
              c_hTA_an_fl(i) = cF_limh(i)*c_hTA_an(i)
            enddo
*    ..calculate the classical ion conductivity with flux limiting
            call b2trht_an(nspecies, ix, iy, bb(ix,iy,0)/bb(ix,iy,3),
     &       zetap(ix,iy), zeffne, ti(ix,iy), c_hTA_an_fl, nz2n, rmd2,
     &       hci0_al(ix,iy,0:nspecies-1))
*    ..the thermal force transport coefficient with flux limiting
            do i = 0, nspecies-1
              temp = sum(c_R_TA_an(i,0:nspecies-1)*
     &         (nz2n(0:nspecies-1)/nz2n(i)))
              c_R_TA(ix,iy,i,0:nspecies-1) =
     &         -c_R_TA_an(0:nspecies-1,i)*(nz2n(i)/nz2n(0:nspecies-1))
              c_R_TA(ix,iy,i,i) = c_R_TA(ix,iy,i,i)+temp
            enddo

*    ..the analytical viscous stress transport coefficient
            call b2trvis_an(nspecies, mum, mdm,nz2n,
     &       c_pi_uA_an, c_pi_hA_an)

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
     &          -(zetap(ix,iy)/(zeffne**2*mp))*
     &           (c_hTA_an(i)/rmd2(i)*nal(ix,iy,i)/nz2n(i)**2)
                kappa_hTB(i) =
     &          -(zetap(ix,iy)/(zeffne**2*mp))*
     &           (c_hTB(i)/rmd2(i)*nal(ix,iy,i)/nz2n(i)**2)

                kappa_rTA(i,i) =
     &          -(zetap(ix,iy)/(zeffne**2*mp))*
     &           (c_rTA_an(i)/rmd2(i)*nal(ix,iy,i)/nz2n(i)**2)
                kappa_rTB(i) =
     &          -(zetap(ix,iy)/(zeffne**2*mp))*
     &           (c_rTB(i)/rmd2(i)*nal(ix,iy,i)/nz2n(i)**2)
              enddo
            endif

          endif

*    ..calculate the thermal force correction for each charge state
          c_R_TB(ix,iy,0:nspecies-1) = cF_limB(0:nspecies-1)*
     &     c_R_TB_nofl(ix,iy,0:nspecies-1)

*    ..calculate the friction force transport coefficient
          call b2trfr(nspecies, mum, rmu,
     &     z2n_xy(ix,iy,0:nspecies-1),
     &     c_hwA(0:nspecies-1,0:nspecies-1),
     &     c_rwA(0:nspecies-1,0:nspecies-1),
     &     c_R_wB(ix,iy,0:nspecies-1,0:nspecies-1),
     &     c_R_wA(ix,iy,0:nspecies-1,0:nspecies-1),
     &     c_R_w(ix,iy,0:nspecies-1,0:nspecies-1))

*    ..calculate the classical eta coefficient including charge state correction
          eta_cl(ix,iy,0:ns-1) = 0.0_R8
          do i = 0, nspecies-1
            tau_Zh(i,0:nspecies-1) = rmu(i,0:nspecies-1)*((zetap(ix,iy)*
     &       nal(ix,iy,i))/(z2n_xy(ix,iy,i)*z2n_xy(ix,iy,0:nspecies-1)))

*    ..the velocity-dependent part
            temp = sum(c_pi_uA(i,0:nspecies-1)*tau_Zh(i,0:nspecies-1)*
     &       (nal(ix,iy,0:nspecies-1)/nal(ix,iy,i)))

            do iz= 1, znnucl_int(i)
              eta_cl(ix,iy,nucl2s(i,iz)) =
     &         na(ix,iy,nucl2s(i,iz))*ti(ix,iy)*
     &        (temp+zhcscorr*2.0_R8*c_pi_uB(i)*tau_Zh(i,i)*
     &        (z2n_xy(ix,iy,i)/(nal(ix,iy,i)*iz**2)-1.0_R8))
            enddo

*    ..the heat flux dependent part
            do j = 0, nspecies-1
              eta_hAdp(ix,iy,i,j) = c_pi_hA(i,j)*tau_Zh(i,j)
            enddo

            eta_hBdp(ix,iy,i) = 2.0_R8*c_pi_hB(i)*tau_Zh(i,i)
          enddo

*    ..calculate flux limiting for the viscous stress
          call b2trfl_vis (ix, iy,ns, eta_cl(ix,iy,0:ns-1),
     &     na(ix,iy,0:ns-1)*ti(ix,iy), gamma, cF_limv)

*    ..the velocity-dependent part
          do is = 0, ns-1
*    ..apply bx^2
            eta_cl(ix,iy,is) = eta_cl(ix,iy,is)*
     &       (bb(ix,iy,0)/bb(ix,iy,3))**2
*    ..apply flux limiting
            eta_cl_fl(ix,iy,is) = eta_cl(ix,iy,is)*
     &       cF_limv(is)
          enddo

*    ..the heat flux dependent part
          do i = 0, nspecies-1
*    ..apply bx^2
            do j = 0, nspecies-1
              eta_hAdp(ix,iy,i,j) = eta_hAdp(ix,iy,i,j)*
     &         (bb(ix,iy,0)/bb(ix,iy,3))**2
            enddo

            eta_hBdp(ix,iy,i) = eta_hBdp(ix,iy,i)*
     &       (bb(ix,iy,0)/bb(ix,iy,3))**2
          enddo

*   ..tilde I for each isonuclear sequence
          do inucl = 0, nspecies-1
            I_al(inucl)=0.0_R8
            do iz = 1, znnucl_int(inucl)
              I_al(inucl) = I_al(inucl)+
     &         (na(ix,iy,nucl2s(inucl,iz))/nal(ix,iy,inucl))**2/
     &          Ia(ix,iy,nucl2s(inucl,iz))
            enddo
          enddo

*    ..calculate the classical ion conductivity correction
          call b2trht_an(nspecies, ix, iy, bb(ix,iy,0)/bb(ix,iy,3),
     &     zetap(ix,iy), zeffne, ti(ix,iy), c_hTB, nz2n, rmd2,
     &     hciB_al_nofl(0:nspecies-1))
*    ..calculate the classical ion conductivity correction with flux limiting
          call b2trht_an(nspecies, ix, iy, bb(ix,iy,0)/bb(ix,iy,3),
     &     zetap(ix,iy), zeffne, ti(ix,iy), cF_limBh*c_hTB, nz2n, rmd2,
     &     hciB_al(0:nspecies-1))

*   ..apply charge state correction for the heat flux
          hci0_nofl(ix,iy,0) =
     &     sum(hci0_al_nofl(ix,iy,0:nspecies-1)+zhcscorr*
     &     hciB_al_nofl(0:nspecies-1)*(I_al(0:nspecies-1)-1.0_R8))
          hci0(ix,iy,0) =
     &     sum(hci0_al(ix,iy,0:nspecies-1)+zhcscorr*
     &     hciB_al(0:nspecies-1)*(I_al(0:nspecies-1)-1.0_R8))

*   ..ion heat conductivity for each charge state
          hci0_a_nofl(ix,iy,0:ns-1) = 0.0_R8
          hci0_a(ix,iy,0:ns-1) = 0.0_R8
          do i = 0, nspecies-1
            do iz= 1, znnucl_int(i)
              hci0_a_nofl(ix,iy,nucl2s(i,iz)) =
     &         na(ix,iy,nucl2s(i,iz))/nal(ix,iy,i)*
     &        (hci0_al_nofl(ix,iy,i)+zhcscorr*hciB_al_nofl(i)*
     &        (z2n_xy(ix,iy,i)/(nal(ix,iy,i)*iz**2)-1.0_R8))
              hci0_a(ix,iy,nucl2s(i,iz)) =
     &         na(ix,iy,nucl2s(i,iz))/nal(ix,iy,i)*
     &        (hci0_al(ix,iy,i)+zhcscorr*hciB_al(i)*
     &         (z2n_xy(ix,iy,i)/(nal(ix,iy,i)*iz**2)-1.0_R8))
            enddo
          enddo

          c_hw(0:nspecies-1,0:nspecies-1) =
     &     c_hwA(0:nspecies-1,0:nspecies-1)
          do inucl = 0, nspecies-1
            c_hw(inucl,inucl) = c_hw(inucl,inucl)+zhcscorr*
     &       c_hwB(inucl)-sum(c_hwA(0:nspecies-1,inucl))
          enddo
*   ..save heat flux velocity-dependent part transport coefficient
          c_hw_save(ix,iy,0:nspecies-1,0:nspecies-1) =
     &     c_hw(0:nspecies-1,0:nspecies-1)

          if (zh_dist_fun.eq.1) then
*   ..apply charge state correction for the R-moment
            c_rw(0:nspecies-1,0:nspecies-1) =
     &       c_rwA(0:nspecies-1,0:nspecies-1)
            do inucl = 0, nspecies-1
              c_rw(inucl,inucl) = c_rw(inucl,inucl)+zhcscorr*
     &         c_rwB(inucl)-sum(c_rwA(0:nspecies-1,inucl))
            enddo

*    ..moments of the distribution function
            do inucl = 0, nspecies-1
              do iz= 1, znnucl_int(inucl)
                w_out(ix,iy,nucl2s(inucl,iz)) =
     &           ua(ix,iy,nucl2s(inucl,iz))-avm_u(ix,iy)

                hTdp_out(ix,iy,nucl2s(inucl,iz)) =
     &           sum(kappa_hTA(inucl,0:nspecies-1)*
     &               gtalc(ix,iy,0:nspecies-1))+
     &           zhcscorr*kappa_hTB(inucl)*
     &          (na(ix,iy,nucl2s(inucl,iz))/nal(ix,iy,inucl)/
     &           Ia(ix,iy,nucl2s(inucl,iz))*
     &           gt_Ac(ix,iy,nucl2s(inucl,iz))-
     &           gtalc(ix,iy,inucl))

                hwdp_out(ix,iy,nucl2s(inucl,iz)) =
     &           sum(c_hw(0:nspecies-1,inucl)*
     &               (ua(ix,iy,nucl2s(inucl,iz))-
     &                av_ualpha(ix,iy,0:nspecies-1)))

                rTdpgamma_out(ix,iy,nucl2s(inucl,iz)) =
     &           sum(kappa_rTA(inucl,0:nspecies-1)*
     &               gtalc(ix,iy,0:nspecies-1))+
     &           zhcscorr*kappa_rTB(inucl)*
     &          (na(ix,iy,nucl2s(inucl,iz))/nal(ix,iy,inucl)/
     &           Ia(ix,iy,nucl2s(inucl,iz))*
     &           gt_Ac(ix,iy,nucl2s(inucl,iz))-
     &           gtalc(ix,iy,inucl))

                rwdpgamma_out(ix,iy,nucl2s(inucl,iz)) =
     &           sum(c_rw(0:nspecies-1,inucl)*
     &              (ua(ix,iy,nucl2s(inucl,iz))-
     &               av_ualpha(ix,iy,0:nspecies-1)))
              enddo
            enddo
          endif

        enddo
      enddo

*    ..ion conductivity on the cell boundary
      call b2txcx (nx, ny, mode, vol, gs, hx, qz,
     & hci0(-1,-1,0), chci0(-1,-1,0))
      call b2txcx (nx, ny, mode, vol, gs, hx, qz,
     & hci0_nofl(-1,-1,0), chci0_nofl(-1,-1,0))
      do i = 0, nspecies-1
        call b2txcx (nx, ny, mode, vol, gs, hx, qz,
     &   hci0_al(-1,-1,i), chci0_al(-1,-1,0,i))
      enddo
      do is = 0, ns-1
        call b2txcx (nx, ny, mode, vol, gs, hx, qz,
     &   hci0_a(-1,-1,is), chci0_a(-1,-1,0,is))
      enddo

*    ..ion velocity-dependent viscosity on the cell boundary
      do is = 0, ns-1
       call b2txcx (nx, ny, mode, hz*vol, gs, hx, qz,
     &   eta_cl_fl(-1:nx,-1:ny,is)*rho_A_rel(-1:nx,-1:ny,is),
     &   cvsahz0(-1,-1,0,is))
       call b2txcx (nx, ny, mode, hz*vol, gs, hx, qz,
     &   eta_cl_fl(-1,-1,is), cvsahz_drho(-1,-1,0,is))
       call b2txcx (nx, ny, mode, hz*vol, gs, hx, qz,
     &   eta_cl(-1:nx,-1:ny,is)*rho_A_rel(-1:nx,-1:ny,is),
     &   cvsahz0_nofl(-1,-1,0,is))
       call b2txcx (nx, ny, mode, vol, gs, hx, qz,
     &   eta_cl_fl(-1:nx,-1:ny,is)*rho_A_rel(-1:nx,-1:ny,is),
     &   cvsa0(-1,-1,0,is))
       call b2txcx (nx, ny, mode, vol, gs, hx, qz,
     &   eta_cl_fl(-1,-1,is), cvsa_drho(-1,-1,0,is))
      enddo

*    ..ion heat flux dependent viscosity on the cell boundary
      do i = 0, nspecies-1
        do j = 0, nspecies-1
          call b2txcx (nx, ny, mode, hz*vol, gs, hx, qz,
     &     eta_hAdp(-1,-1,i,j), cvsahz_hAdp(-1,-1,0,i,j))
          call b2txcx (nx, ny, mode, vol, gs, hx, qz,
     &     eta_hAdp(-1,-1,i,j), cvsa_hAdp(-1,-1,0,i,j))
        enddo

        call b2txcx (nx, ny, mode, hz*vol, gs, hx, qz,
     &   eta_hBdp(-1,-1,i), cvsahz_hBdp(-1,-1,0,i))
        call b2txcx (nx, ny, mode, vol, gs, hx, qz,
     &   eta_hBdp(-1,-1,i), cvsa_hBdp(-1,-1,0,i))
      enddo

*    ..apply Luciani correction for the ion conductivity
      do iy = -1, ny
        do ix = 0, nx
          if (on_closed_surface(ix,iy)) then
            chci0_luc(ix,iy,0) = chci0(ix,iy,0)*f_luc_ki(ix,iy)
          else
            chci0_luc(ix,iy,0) = chci0(ix,iy,0)
          endif
        enddo
      enddo
      do i = 0, nspecies-1
        do iy = -1, ny
          do ix = 0, nx
            if (on_closed_surface(ix,iy)) then
              chci0_al_luc(ix,iy,0,i) = chci0_al(ix,iy,0,i)*
     &         f_luc_ki(ix,iy)
            else
              chci0_al_luc(ix,iy,0,i) = chci0_al(ix,iy,0,i)
            endif
          enddo
        enddo
      enddo
      do is = 0, ns-1
        do iy = -1, ny
          do ix = 0, nx
            if (on_closed_surface(ix,iy)) then
              chci0_a_luc(ix,iy,0,is) = chci0_a(ix,iy,0,is)*
     &         f_luc_ki(ix,iy)
            else
              chci0_a_luc(ix,iy,0,is) = chci0_a(ix,iy,0,is)
            endif
          enddo
        enddo
      enddo

*    ..apply Luciani correction for the ion velocity-dependent viscosity
      do is = 0, ns-1
        if(.not.is_neutral(is)) then
          do iy = -1, ny
            do ix = 0, nx
              if (on_closed_surface(ix,iy)) then
                cvsahz0_luc(ix,iy,0,is) = cvsahz0(ix,iy,0,is)*
     &           f_luc_et(ix,iy)
                cvsahz_drho_luc(ix,iy,0,is) = cvsahz_drho(ix,iy,0,is)*
     &           f_luc_et(ix,iy)
                cvsa0_luc(ix,iy,0,is) = cvsa0(ix,iy,0,is)*
     &           f_luc_et(ix,iy)
                cvsa_drho_luc(ix,iy,0,is) = cvsa_drho(ix,iy,0,is)*
     &           f_luc_et(ix,iy)
              else
                cvsahz0_luc(ix,iy,0,is) = cvsahz0(ix,iy,0,is)
                cvsahz_drho_luc(ix,iy,0,is) = cvsahz_drho(ix,iy,0,is)
                cvsa0_luc(ix,iy,0,is) = cvsa0(ix,iy,0,is)
                cvsa_drho_luc(ix,iy,0,is) = cvsa_drho(ix,iy,0,is)
              endif
            enddo
          enddo
        endif
      enddo

*    ..apply Luciani and neoclassical corrections for the ion heat flux
*      dependent viscosity
      do i = 0, nspecies-1
        do j = 0, nspecies-1
          do iy = -1, ny
            do ix = 0, nx
              if (on_closed_surface(ix,iy)) then
                cvsahz_hAdp_luc(ix,iy,0,i,j) = cvsahz_hAdp(ix,iy,0,i,j)*
     &           f_luc_et(ix,iy)*alpha1(ix,iy)
                cvsa_hAdp_luc(ix,iy,0,i,j) = cvsa_hAdp(ix,iy,0,i,j)*
     &           f_luc_et(ix,iy)*alpha1(ix,iy)
              else
                cvsahz_hAdp_luc(ix,iy,0,i,j) = cvsahz_hAdp(ix,iy,0,i,j)
                cvsa_hAdp_luc(ix,iy,0,i,j) = cvsa_hAdp(ix,iy,0,i,j)
              endif
            enddo
          enddo
        enddo

        do iy = -1, ny
          do ix = 0, nx
            if (on_closed_surface(ix,iy)) then
              cvsahz_hBdp_luc(ix,iy,0,i) = cvsahz_hBdp(ix,iy,0,i)*
     &         f_luc_et(ix,iy)*alpha1(ix,iy)
              cvsa_hBdp_luc(ix,iy,0,i) = cvsa_hBdp(ix,iy,0,i)*
     &         f_luc_et(ix,iy)*alpha1(ix,iy)
            else
              cvsahz_hBdp_luc(ix,iy,0,i) = cvsahz_hBdp(ix,iy,0,i)
              cvsa_hBdp_luc(ix,iy,0,i) = cvsa_hBdp(ix,iy,0,i)
            endif
          enddo
        enddo
      enddo

      ncall = ncall+1
      call subend ()
      return

      end subroutine b2trzh
