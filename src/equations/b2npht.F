*-----------------------------------------------------------------------
*.specification

      subroutine b2npht (nCv, nFc, nVx, ns, switch, geo, mpg,           !srv 05.07.17
     &  itcnt, ismain, solving, solvereg, solvireg, solvnreg, solvtreg, !srv 22.05.18
     &  solvpreg, solvmreg, solvkreg, solvzreg, rxf, dtim,
     &  pl, dv, co, rt, sr, srw, psnc, psnl, st_ext, ierr)
      use b2mod_types
      use b2mod_constants
      use b2mod_b2cmpa
      use b2mod_switches
      use b2us_geo
      use b2us_map
      use b2us_plasma
      use b2mod_numerics_namelist                                     !srv 07.07.21 {
     & , only : time_factor, dtei, dtee, dten                         !srv 07.07.21 }
!WG_TODO      use b2mod_balance !djm Jan2017
!WG_TODO     & , only : b2npht_shei_bal, update_balance_sources_heat,
!WG_TODO     &          balance_netcdf
      use b2mod_ad
     & , only : ncall_b2npht
#ifndef NAGFOR
! csc The following are not necessary for computation but are needed
!     for adjoint AD to avoid side-effect variables
      use b2mod_numerics_namelist
     & , only : last_solve_9
      use b2mod_ad
     & , only : my_out_folder,
#ifdef WG_DEBUG
     &          my_out_sf,                                              !srv 27.06.21
#endif
     &          ncall_b2sikt, ncall_b2sifrtf,
     &          ncall_b2xehx, ncall_b2xehy, ncall_b2upht, ncall_b2usht,
     &          ncall_b2ursd, ncall_b2sihs_, ncall_b2tlnl, cvregmax
#endif
      use b2mod_subsys
      implicit none
*   ..input arguments (unchanged on exit)
      integer nCv, nFc, nVx, ns, itcnt, ismain                          !srv 05.07.17
      type (switches), intent (in) :: switch
      type (geometry), intent (in) :: geo
      type (mapping), intent (in) :: mpg
      type (B2PlasmaSnapshot), intent (inout) :: psnc, psnl
      real (kind=R8) :: rxf, dtim
      logical solving(4),
     *  solvereg(0:mpg%nnreg(0)), solvireg(0:mpg%nnreg(0)),
     *  solvnreg(0:mpg%nnreg(0)), solvtreg(0:mpg%nnreg(0)),             !srv 22.05.18
     *  solvpreg(0:mpg%nnreg(0)), solvmreg(0:ns-1,0:mpg%nnreg(0)),
     *  solvkreg(0:mpg%nnreg(0)), solvzreg(0:mpg%nnreg(0))
*   ..input/output arguments
      type (B2Plasma), intent (inout) :: pl
      type (B2Derivatives), intent (inout) :: dv
      type (B2Coeff), intent (inout) :: co
      type (B2Source), intent (inout) :: sr
      type (B2SourceWork), intent (inout) :: srw
      type (B2Rates), intent (in) :: rt
      type (B2StateExt), intent (in) :: st_ext
*   ..output arguments (unspecified on entry)
      integer ierr
*   ..common blocks
*-----------------------------------------------------------------------
*.documentation
*
*  1. purpose
*
*     B2NPHT performs one iteration of a relaxation process for the
*     solution of the equations for internal electron and atom energy.
*
*
*  3. description (see also routine b2cdca)
*
*     For a detailed description of the equations and the method of
*     solution, see the references cited in the main program unit.
*
*     This routine performs one inner iteration for the equations that
*     govern the temperatures te and ti.
*     An outline of the procedure follows.
*     ..compute auxiliaries
*     ..compute residuals
*     ..compute corrections
*     ..apply corrections
*
*-----------------------------------------------------------------------
*.declarations

*   ..local variables
      integer iCv, is
      real (kind=R8) ::
     * aa(mpg%nCmxNv), wrkf(nFc,0:1),                                   !srv 02.01.07
     1 wrkc(nFc,0:1,0:2),                                               !srv 02.01.07
     2 she0(nCv,0:3), shi0(nCv,0:3), shn0(nCv,0:3),
     3 wrk0(nCv), shei(nCv),                                            !srv 18.09.09
     4 shef(nCv,0:3), shif(nCv,0:3),                                    !srv 10.11.02
     5 shekt(nCv,0:3), shikt(nCv,0:3), skt0(nCv,0:3), szt0(nCv,0:3),
     6 skt_prod(nCv), skt_diss(nCv)
      character chns*3                                                  !srv 17.07.05
      real (kind=R8), save :: rxg
      logical :: warning_i, warning_e, warning_n, warning_kt, warning_zt
*   ..procedures
      intrinsic sqrt
      external xertst, ipgetr, b2saxpy
      external b2xvsg, b2xvff, b2xvfx, b2ursd, b2usht, b2upht           !srv 02.01.07
*   ..initialisation

*-----------------------------------------------------------------------
*.documentation-internal
*
*     Description of some local variables:
*
*     rxg - real.
*     rxg specifies a special under-relaxation parameter. A term of the
*     form "abs(residual)/rxg" is added to the diagonal of the matrix
*     of the correction equation, with the effect of limiting the
*     computed correction. rxg is dimensionless and of order unity;
*     smaller values of rxg imply stronger damping.
*
*     pcm0 - real.
*     ...
*
*     pcm1 - real.
*     ...
*
*-----------------------------------------------------------------------
*.computation

* ..preliminaries
*   ..subprogram start-up calls
      call subini ('b2npht')

*   ..set internal parameters on first call
      if (ncall_b2npht.eq.0) then
        if(switch%b2mndt_style.eq.2) then
          rxg = huge(1.0_R8)
          if (switch%b2npht_pcm0.gt.0.0_R8) then
            call xerrab('b2npht_pcm0 must be 0 for time-dependent mode')
          endif
          if (switch%b2npht_pcm1.gt.0.0_R8) then
            call xerrab('b2npht_pcm1 must be 0 for time-dependent mode')
          endif
          warning_e = .false.
          warning_i = .false.
          warning_n = .false.
          warning_kt = .false.
          warning_zt = .false.
          do iCv = 1, nCv
            warning_e = warning_e .or.
     .        rxf*dtee(mpg%cvReg(iCv))*time_factor(iCv).ne.1.0_R8
            warning_i = warning_i .or.
     .        rxf*dtei(mpg%cvReg(iCv))*time_factor(iCv).ne.1.0_R8
            warning_n = warning_n .or.
     .        rxf*dten(mpg%cvReg(iCv))*time_factor(iCv).ne.1.0_R8
            if (switch%solve_keps.gt.0) then
              warning_kt = warning_kt.or.
     .        rxf*time_factor(iCv).ne.1.0_R8
            endif
            if (switch%solve_keps.gt.1) then
              warning_zt = warning_zt.or.
     .        rxf*time_factor(iCv).ne.1.0_R8
            endif
          enddo
          if (warning_e) write(6,*) 'rxf*dtee*time_factor should be 1 !'
          if (warning_i) write(6,*) 'rxf*dtei*time_factor should be 1 !'
          if (warning_n) write(6,*) 'rxf*dten*time_factor should be 1 !'
          if (warning_kt.or.warning_zt) then
            write(6,*) 'rxf*time_factor should be 1 for k-eps !'
          endif


          if (warning_e .or. warning_i .or. warning_n)
     .   write(6,*) 'Time derivative terms may be inacurrate !'
          if (warning_kt .or. warning_zt)
     .   write(6,*) 'Time derivative terms may be inacurrate for k-eps!'
        else
          rxg = 1.0_R8
        endif                                                           !srv 07.07.21 }
      endif
*   ..test nCv, nFc
      call xertst (0.le.nCv.and.0.le.nFc, 'faulty argument nCv, nFc')
*   ..test rxf
      call xertst (0.le.rxf.and.rxf.le.1, 'faulty argument rxf')
*   ..extensive tests on first few calls
      if (ncall_b2npht.lt.3) then
*    ..test sign of transport coefficients
       call b2xvsg (nCv, co%ceqp, 1, 'ceqp', '.ge.')
*    ..test sign of she, shi
       call b2xvsg (nCv, sr%she(1,0), 1, 'she0', '.ge.')
       call b2xvsg (nCv, sr%she(1,1), 1, 'she1', '.le.')
       call b2xvsg (nCv, sr%she(1,2), 1, 'she2', '.ge.')
       call b2xvsg (nCv, sr%she(1,3), 1, 'she3', '.le.')
       call b2xvsg (nCv, sr%shi(1,0), 1, 'shi0', '.ge.')
       call b2xvsg (nCv, sr%shi(1,1), 1, 'shi1', '.le.')
       call b2xvsg (nCv, sr%shi(1,2), 1, 'shi2', '.ge.')
       call b2xvsg (nCv, sr%shi(1,3), 1, 'shi3', '.le.')
       call b2xvsg (nCv, sr%shn(1,0), 1, 'shn0', '.ge.')
       call b2xvsg (nCv, sr%shn(1,1), 1, 'shn1', '.le.')
       call b2xvsg (nCv, sr%shn(1,2), 1, 'shn2', '.ge.')
       call b2xvsg (nCv, sr%shn(1,3), 1, 'shn3', '.le.')
*    ..test sign of te, ti
       call b2xvsg (nCv, pl%te, 1, 'te', '.gt.')
       call b2xvsg (nCv, pl%ti, 1, 'ti', '.gt.')
       call b2xvsg (nCv, pl%tn, 1, 'tn', '.gt.')
      endif
*    ..throw warning for use of time-dependent with new k-eps equations
      if (switch%solve_keps.gt.0.and.switch%b2mndt_style.eq.2) then
        write(6,*) 'WARNING: the time-dependent treatment of k-eps',
     &              'equations was implemented similarly as the energy',
     &              'equations but lacks theoretical development'
      endif


* ..main computation

c make sure fhn is zero when not solving a separate neutral energy equation
      if (switch % tn_style.ne.2) then
        dv%fhn = 0.0_R8
      endif

c time-dependent treatment                                              !srv 07.07.21 {
      if(switch%b2mndt_style.eq.2) then
*  .. calculate time-dependent source
        call b2shdt (nCv, mpg%nCi, ns, dtim, switch%BoRiS,
     &    geo%cvVol, psnl%na, psnl%te, psnl%ti, psnl%tn,
     &    psnl%kt, psnl%zt, dv%ne, dv%ni, dv%nn,
     &    psnl%kinrgy, switch, mpg, sr%shedt, sr%shidt, sr%shndt,
     &    sr%sktdt, sr%sztdt)
*  .. add contribution from time-dependent source
        sr%she = sr%she+sr%shedt
        sr%she(:,0)=sr%she(:,0)-1.5_R8*(dv%ne-psnl%ne)*
     &                                  psnl%te/dtim*geo%cvVol
        sr%shi = sr%shi+sr%shidt
        sr%shi(:,0)=sr%shi(:,0)-1.5_R8*(dv%ni(:,0)-psnl%ni(:,0))* !! TODO B2MNDT depends on tn_style ?
     &                                  psnl%ti/dtim*geo%cvVol
        if (switch%tn_style.eq.2) then
          sr%shn = sr%shn+sr%shndt
          sr%shn(:,0)=sr%shn(:,0)-1.5_R8*(dv%nn(:)-psnl%nn(:))*
     &                                    psnl%tn/dtim*geo%cvVol
        endif
        if (switch%solve_keps.gt.0) then
          sr%skt = sr%skt + sr%sktdt
          sr%skt(:,0)=sr%skt(:,0)-1.5_R8*(dv%ni(:,1)-psnl%ni(:,1))*
     &                                    psnl%kt/dtim*geo%cvVol
        endif
        if (switch%solve_keps.gt.1) then
          sr%szt = sr%szt + sr%sztdt
          sr%szt(:,0)=sr%szt(:,0)-1.5_R8*(dv%ni(:,1)-psnl%ni(:,1))*
     &                                    psnl%zt/dtim*geo%cvVol
        endif
      endif                                                              !srv 207.07.21 }

*   ..compute final heat sources
      if(switch%b2npht_style.eq.0) then
       call xerrab('b2npht -- style 0 not available for WG')
#ifdef WG_TODO
       call b2sihs (nCv, nFc, nVx, ns, geo, mpg, switch%BoRiS, ismain,
     &  dv%fac_ExB, dv%lnlam, pl%na, pl%ua, dv%ue, pl%te, pl%ti, pl%po,
     &  dv%ne, dv%ni, dv%fch, co%cvsa, dv%fchdia, dv%uadia, dv%wadia,
     &  dv%vaecrb, she0, shi0)
#endif
      else
       call b2sihs_ (nCv, nFc, nVx, ns, switch, geo, mpg, ismain,
     &  dv%fac_ExB, dv%lnlam, pl%na, pl%ua, dv%ue, pl%te, pl%ti, pl%tn,
     &  pl%po, dv%ne, dv%ni, dv%nn, dv%fna_fcor, dv%fne, dv%fch,
     &  co%cvsa, co%cvsahz_eff, co%f_luc_sg, co%alfx_c, co%sigx_c,
     &  rt%rza, rt%rz2, st_ext, she0, shi0, shn0, srw)
      endif
      call b2saxpy (nCv*4, 1.0e0_R8, sr%she, 1, she0, 1)
      call b2saxpy (nCv*4, 1.0e0_R8, sr%shi, 1, shi0, 1)
      call b2saxpy (nCv*4, 1.0e0_R8, sr%shn, 1, shn0, 1)

*   ..sources due to k-eps model
      sr%skt_prod = 0.0_R8
      sr%skt_diss = 0.0_R8
      if (switch%solve_keps.gt.0) then
        call b2sikt (nCv, nFc, nVx, ns, ismain, switch, geo, mpg,
     &    pl, dv, rt%rza, co, st_ext, shekt, shikt, skt0, skt_prod,
     &    skt_diss)
        call b2saxpy (nCv*4, 1.0e0_R8, shekt, 1, she0, 1)
        call b2saxpy (nCv*4, 1.0e0_R8, shikt, 1, shi0, 1)
        call b2saxpy (nCv*4, 1.0e0_R8, sr%skt, 1, skt0, 1)
        call b2saxpy (nCv, 1.0e0_R8, skt_prod, 1, sr%skt_prod, 1)
        call b2saxpy (nCv, 1.0e0_R8, skt_diss, 1, sr%skt_diss, 1)
      endif

#ifdef WG_TODO
cdjm Jan2017 update particle source components for balance
      if (balance_netcdf.ne.0) then
        call update_balance_sources_heat(nx,ny,ne,ni,te,ti)
      endif
#endif
*   ..compute provisional reshe, reshi
*     (first ignore equipartition)
*   ..save final electron heat source ignoring equipartition            !srv 10.11.02 {
      shef = she0
*   ..save final ion heat source ignoring equipartition
      shif = shi0                                                       !srv 10.11.02 }
      if ((solving(4).and.any(solvereg(0:mpg%nnreg(0)))).or.
     &     switch%get_residuals.eq.1) then                              !srv 22.05.18
        if (switch%mdf_fhe.eq.0) then
          call b2ursd (nCv, nFc, mpg, dv%ne, pl%te, she0, dv%fhe,
     &     dv%reshe)
        else
          call b2ursd (nCv, nFc, mpg, dv%ne, pl%te, she0,
     &     dv%fhe_mdf, dv%reshe)
        endif
        if(switch%b2mndt_style.eq.2) then                               !srv 10.10.22 {
          dv%reshe=dv%reshe+
     &    (pl%te-psnl%te)*dv%ne*sr%shedt(:,3)
        endif                                                           !srv 10.10.22 }
      endif                                                             !srv 22.05.18
      if ((solving(4).and.any(solvireg(0:mpg%nnreg(0)))).or.
     &     switch%get_residuals.eq.1) then    !srv 22.05.18
        if (switch%mdf_fhi.eq.0) then
          call b2ursd (nCv, nFc, mpg, dv%ni(:,0), pl%ti, shi0,
     &     dv%fhi, dv%reshi)
        else
          call b2ursd (nCv, nFc, mpg, dv%ni(:,0), pl%ti, shi0,
     &     dv%fhi_mdf, dv%reshi)
        endif
        if(switch%b2mndt_style.eq.2) then                               !srv 10.10.22 {
          dv%reshi=dv%reshi+
     &    (pl%ti-psnl%ti)*dv%ni(:,0)*sr%shidt(:,3)
        endif                                                           !srv 10.10.22 }
      endif                                                             !srv 22.05.18
      if ((switch%solve_keps.gt.0.and.any(solvkreg(0:mpg%nnreg(0)))).or.
     &     switch%get_residuals.eq.1) then
        call b2ursd (nCv, nFc, mpg, dv%ni(:,1), pl%kt, skt0,
     &   dv%fkt, dv%reskt)
      endif
      if ((switch%solve_keps.gt.1.and.any(solvzreg(0:mpg%nnreg(0)))).or.
     &     switch%get_residuals.eq.1) then
        call b2ursd (nCv, nFc, mpg, dv%ni(:,1), pl%zt, szt0,
     &   dv%fzt, dv%reszt)
      endif

c residuals for separate neutral energy equation
      if (switch % tn_style.eq.2) then   !! separate neutral energy equation
        if (any(solvnreg(0:mpg%nnreg(0))).or.
     &      switch%get_residuals.eq.1) then
                                !! AFN_TODO add neutral energy equation in solving switch?
          call b2ursd (nCv, nFc, mpg, dv%nn, pl%tn, shn0, dv%fhn,
     &     dv%reshn)
        endif
        if(switch%b2mndt_style.eq.2) then                               !srv 10.10.22 {
          dv%reshn=dv%reshn+
     &    (pl%tn-psnl%tn)*dv%ni(:,0)*sr%shndt(:,3) !! TODO_B2MNDT double check source scaling?
        endif                                                           !srv 10.10.22 }
      else
        dv%reshn = 0.0_R8         !! this is needed because reshn is added to the total energy equation
      endif

      do iCv = 1, mpg%nCi
        shei(iCv) = co%ceqp(iCv)*dv%ne(iCv)*dv%ni(iCv,0)/
     &     (pl%te(iCv)*sqrt(pl%te(iCv)))*(pl%te(iCv)-pl%ti(iCv))
      enddo
      do iCv = mpg%nCi+1, nCv
        shei(iCv) = 0.0_R8
      enddo
#ifdef WG_TODO
      if (balance_netcdf.ne.0) then
        b2npht_shei_bal(-1:nx,-1:ny) = shei(-1:nx,-1:ny)                !djm Jan2017
      endif
#endif
*   ..compute final resht, reshe, reshi
!      if (solving(4).or.get_residuals.eq.1) then                        !srv 22.05.18
!        dv%resht = dv%reshe + dv%reshi
!        dv%reshe = dv%reshe - shei
!        dv%reshi = dv%reshi + shei
!      endif                                                             !srv 22.05.18
      if ((solving(4).and.any(solvereg(0:mpg%nnreg(0))))                     !srv 22.05.18 {
     &    .or.switch%get_residuals.eq.1) dv%reshe = dv%reshe - shei
      if ((solving(4).and.any(solvireg(0:mpg%nnreg(0))))
     &    .or.switch%get_residuals.eq.1) dv%reshi = dv%reshi + shei
      if (solving(4).and.any(solvereg(0:mpg%nnreg(0))).and.            ! if total equation is not solved then resmt = 0
     &   any(solvireg(0:mpg%nnreg(0))).and.
     &   any(solvtreg(0:mpg%nnreg(0))).or.
     &   switch%get_residuals.eq.1)
     &   dv%resht = dv%reshe + dv%reshi + dv%reshn                       !srv 22.05.18 }

      if(switch%b2npht_iout .eq. 1 .or. switch%iout_b2wdat .eq. 4) then  !srv 17.07.05 01.11.10 {
        call my_out_us(70,nCv,0,shei,'b2npht_shei')
        call my_out_us(70,nCv,0,co%ceqp,'b2npht_ceqp')
!        do is = 0,ns-1
!         write (chns,'(i3.3)') is
!         call my_out_us(70,nCv,0,na(1,is),'b2npht_na'//chns//'_input')
!         call my_out_us(70,nCv,0,ua(1,is),'b2npht_ua'//chns//'_input')
!        enddo
!        call my_out_us(70,nCv,0,te,'b2npht_te_input')
!        call my_out_us(70,nCv,0,ti,'b2npht_ti_input')
!        call my_out_us(70,nCv,0,po,'b2npht_po_input')
        call my_out_us(70,nCv,0,dv%reshe,'b2npht_reshe')
        call my_out_us(70,nFc,1,dv%fhe(1,0),'b2npht_fhe_th')
        call my_out_us(70,nFc,1,dv%fhe(1,1),'b2npht_fhe_r')
        call my_out_us(70,nCv,0,she0(1,0),'b2npht_she0')
        call my_out_us(70,nCv,0,she0(1,1),'b2npht_she1')
        call my_out_us(70,nCv,0,she0(1,2),'b2npht_she2')
        call my_out_us(70,nCv,0,she0(1,3),'b2npht_she3')
        wrk0=she0(:,0)+she0(:,1)*pl%te+she0(:,2)*dv%ne+
     &       she0(:,3)*pl%te*dv%ne
        call my_out_us(70,nCv,0,wrk0,'b2npht_she')
        call my_out_us(70,nCv,0,dv%reshi,'b2npht_reshi')
        call my_out_us(70,nCv,0,dv%reshn,'b2npht_reshn')
        call my_out_us(70,nCv,0,dv%resht,'b2npht_resht')                  !srv 22.05.18
        call my_out_us(70,nFc,1,dv%fhi(1,0),'b2npht_fhi_th')
        call my_out_us(70,nFc,1,dv%fhi(1,1),'b2npht_fhi_r')
        call my_out_us(70,nFc,1,dv%fhn(1,0),'b2npht_fhn_th')
        call my_out_us(70,nFc,1,dv%fhn(1,1),'b2npht_fhn_r')
        call my_out_us(70,nFc,1,dv%fkt(1,0),'b2npht_fkt_th')
        call my_out_us(70,nFc,1,dv%fkt(1,1),'b2npht_fkt_r')
        call my_out_us(70,nFc,1,dv%fzt(1,0),'b2npht_fzt_th')
        call my_out_us(70,nFc,1,dv%fzt(1,1),'b2npht_fzt_r')
        call my_out_us(70,nCv,0,shi0(1,0),'b2npht_shi0')
        call my_out_us(70,nCv,0,shi0(1,1),'b2npht_shi1')
        call my_out_us(70,nCv,0,shi0(1,2),'b2npht_shi2')
        call my_out_us(70,nCv,0,shi0(1,3),'b2npht_shi3')
        wrk0=shi0(:,0)+shi0(:,1)*pl%ti+
     &       shi0(:,2)*dv%ni(:,0)+shi0(:,3)*pl%ti*dv%ni(:,0)
        call my_out_us(70,nCv,0,wrk0,'b2npht_shi')
        wrk0=shn0(:,0)+shn0(:,1)*pl%tn+
     &       shn0(:,2)*dv%nn+
     &       shn0(:,3)*pl%tn*dv%nn
        call my_out_us(70,nCv,0,wrk0,'b2npht_shn')
        wrk0=skt0(:,0)+skt0(:,1)*pl%kt+
     &       skt0(:,2)*dv%ni(:,1)+skt0(:,3)*pl%kt*dv%ni(:,1)
        call my_out_us(70,nCv,0,wrk0,'b2npht_skt')
        wrk0=szt0(:,0)+szt0(:,1)*pl%zt+
     &       szt0(:,2)*dv%ni(:,1)+szt0(:,3)*pl%zt*dv%ni(:,1)
        call my_out_us(70,nCv,0,wrk0,'b2npht_szt')
      endif                                                             !srv 17.07.05 01.11.10 }
*   ..compute correction
      if (solving(4)) then                                              !srv 11.01.13
        call b2usht (nCv, nFc, nVx, ns, switch, geo, mpg,
     &  itcnt, dv%ne, dv%ni, dv%nn, pl%te, pl%ti,
     &  pl%tn, pl%kt, pl%zt, dv%floe, dv%cone, dv%floi, dv%coni,
     &  dv%flon, dv%conn, dv%flokt, dv%conkt, dv%flozt, dv%conzt,
     &  she0, shi0, shn0, skt0, szt0,
     &  co%ceqp, dv%resht, dv%reshe, dv%reshi, dv%reshn, dv%reskt,
     &  dv%reszt, dv%cortt, dv%corte,
     &  dv%corti, dv%cortn, dv%corkt, dv%corzt, wrkf, wrkc, aa)
*   ..apply correction
        call b2upht (nCv, nFc, nVx, ns, switch, geo, mpg,
     &  solving(3), solving(1), solvereg, solvireg, solvnreg, solvtreg,  !srv 22.05.18
     &  solvpreg, solvmreg, solvkreg, solvzreg, rxf,
     &  dv%cortt, dv%corte, dv%corti, dv%cortn, dv%corkt, dv%corzt,
     &  dv%pccm(1:nCv,0), pl%na, pl%ua, pl%te, pl%ti, pl%tn, pl%kt,
     &  pl%zt, pl%po, dv%ne, dv%ni)
      endif
*   ..set error parameter and return
*     (at present all errors cause an abort through xerrab)
      ierr = 0

      if (switch % tn_style.le.1) then  !! if we solve a combined ion-neutral (tn_style.eq.0) or pure ion (tn_style.eq.1) energy equation,
        pl%tn = pl%ti                   !! tn must be exactly equal to ti.
      endif

      if (switch%b2npht_iout.eq.1) then
        do is = 0, ns-1
          write (chns,'(i3.3)') is
          call my_out_us(70,nCv,0,pl%na(1,is),'b2npht_na'//chns)
          call my_out_us(70,nCv,0,pl%ua(1,is),'b2npht_ua'//chns)
        end do
        call my_out_us(70,nCv,0,pl%te,'b2npht_te')
        call my_out_us(70,nCv,0,pl%ti,'b2npht_ti')
        call my_out_us(70,nCv,0,pl%tn,'b2npht_tn')
        call my_out_us(70,nCv,0,pl%po,'b2npht_po')
        call my_out_us(70,nCv,0,dv%corte,'b2npht_corte')
        call my_out_us(70,nCv,0,dv%corti,'b2npht_corti')
        call my_out_us(70,nCv,0,dv%cortt,'b2npht_cortt')
        call my_out_us(70,nCv,0,dv%cortn,'b2npht_cortn')
      elseif (switch%iout_b2wdat.eq.4) then
        call my_out_us(70,nCv,0,pl%te,'b2npht_te')
        call my_out_us(70,nCv,0,pl%ti,'b2npht_ti')
      end if


* ..return
      ncall_b2npht = ncall_b2npht+1
      call subend ()
      return

*-----------------------------------------------------------------------
*.end b2npht

      end subroutine b2npht
