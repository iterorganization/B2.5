*-----------------------------------------------------------------------
*.specification

      subroutine b2npmo (nCv, nFc, nVx, ns, ismain,
     &  switch, geo, mpg,
     &  itcnt, solving, rxf, dtim, pl, dv, rt, co, co_ns, sr, srw,
     &  psnl, st_ext)
      use b2mod_types
      use b2mod_math
      use b2mod_constants
#ifndef FIXED_POINT
      use b2mod_balance !djm Jan2017
     1 , only : update_balance_sources_momentum,
     2          update_balance_eirene_sources_momentum,
     3          balance_netcdf, fmo_b2nxfv, tot_smo0to3, fmo_tot
#endif
      use b2mod_numerics_namelist
      use b2mod_transport_fun                                          !srv 13.01.17
      use b2mod_b2cmpa
      use b2mod_switches
      use b2us_geo
      use b2us_map
      use b2us_plasma
      use b2mod_openmp
      use b2mod_ma28_for_us
     , , only : restart_ma28_for_us
      use b2mod_ad
     & , only : b2npmo_solvedum, ncall_b2npmo, rxg_npmo
#ifndef NAGFOR
! csc The following are not necessary for computation but are needed
!     for adjoint AD to avoid side-effect variables
      use b2mod_ad
     & , only : ncall_b2usmo, nsmooth_siav_zh,
     &          ncall_b2urmo, ncall_b2sicf, ncall_b2nxfc, ncall_b2sigp,
     &          ncall_b2nxdu, ncall_b2ttia, ncall_b2sifr, ncall_b2siav,
#ifdef WG_DEBUG
     &          my_out_sf,                                             !srv 27.06.21
#endif
     &          my_out_folder, ncall_b2tlnl, ncall_b2sifrtf, icase_sifr
#endif
      use b2mod_zhfrtf                                                 !som 13.08.21
     & , only : rho_A_rel, avm_u, z2n_cv
      use b2mod_b2zhco
     & , only : c_hTA_an_fl_save, c_R_TA_an_save
      use b2mod_frtf_nccorr
     & , only : b2mod_frtf_nccorr_coefs
      use b2mod_subsys
      implicit none
*   ..input arguments (unchanged on exit)
      integer nCv, nFc, nVx, ns, ismain, itcnt                         !srv 13.01.17
      type (switches), intent (in) :: switch
      type (geometry), intent (in) :: geo
      type (mapping), intent (in) :: mpg
      type (B2Coeff_nspecies), intent (in) :: co_ns
      type (B2Rates), intent (in) :: rt
      type (B2Source), intent (inout) :: sr
      type (B2SourceWork), intent (inout) :: srw
      type (B2PlasmaSnapshot), intent (inout) :: psnl
      type (B2StateExt), intent (in) :: st_ext
      ! type(B2State), intent(inout)   :: st
      real (kind=R8) ::
     *  rxf, dtim
      logical solving
*   ..input/output arguments
      type (B2Plasma), intent (inout) :: pl
      type (B2Derivatives), intent (inout) :: dv
      type (B2Coeff), intent (inout) :: co ! due to kt_neo

*-----------------------------------------------------------------------
*.documentation
*
*  1. purpose
*
*     B2NPMO performs one iteration of a relaxation process for the
*     solution of the B2 equations for parallel velocity.
*
*  3. description (see also routine b2cdca)
*
*     For a detailed description of the equations and the method of
*     solution, see the references cited in the main program unit.
*
*     This routine performs one inner iteration for the equations that
*     govern the parallel velocity ua.
*     An outline of the procedure follows.
*     ..compute auxiliaries
*     ..do is=0,ns-1
*      ..compute residual for species (is)
*      ..compute correction for species (is)
*      ..compute contribution to total residual
*     ..enddo
*     ..compute simultaneous correction for all species
*     ..apply correction to all species
*
*-----------------------------------------------------------------------
*.declarations

*   ..local variables
      integer iCv, is, ireg, k
      integer :: vcount(0:ns-1)
      character chns*3                                                 !srv 04.07.08
      real (kind=R8) :: mm(0:ns-1), am1(0:ns-1)                        !som 13.11.20
      real (kind=R8) ::
     *  t0, t1, czetap, czetae                                         !srv 13.01.17
      real (kind=R8) ::
     *  smbch(nCv,0:3), smbtf(nCv,0:3),                                !srv 11.09.09
     *  smbgp(nCv,0:3),
     *  smbvh(nCv,0:3), smbvv(nCv,0:3), smbvi(nCv,0:3),                !srv 23.07.99 04.10.99
     *  smbvtmp(nCv,0:3),
     *  smbcf(nCv,0:3), smban(nCv,0:3),                                !srv 09.07.01 25.09.17
     1  cvct(nFc,0:1), pb(nCv),
     2  flub(nFc,0:1), lnlam(nCv), kt_neo(nCv), rob(nCv),
     *  flubv(nFc,0:1),                                                !srv 11.09.09
     3  wrkf(nFc,0:1), wrkff(nFc,0:1), wrkfff(nFc,0:1),
     *  wrkv(nVx), aa(mpg%nCmxNv),
     4  smb(nCv,0:3), flcb(nFc,0:1),
     5  cvcb(nFc,0:1), ro(nCv), rz(nCv),
     6  pr(nCv), pz(nCv), rob1(nCv), pb1(nCv),
     8  ehxp(nCv), smt(nCv,0:3), cst,                                  !srv 17.01.14
     9  flct(nFc,0:1), wrk1(nCv), wrks(nCv),
     1  smqdu(nCv), wrk0(nCv), wrk2(nCv),
     2  ctcfb(nCv), ctcft(nCv),                                        !srv 06.04.07
     3  zeff(nCv), zmff(nCv),                                          !srv 13.01.17 {
     4  ce1(nCv), ce2(nCv),
     5  cimp1(nCv,0:ns-1), cimp2(nCv,0:ns-1),
     6  zetap(nCv), zetae(nCv),
     7  gti(nCv), gte(nCv), gtn(nCv),                                  !srv 13.01.17 }
     8  wrkf0(nFc), wrkf1(nFc),
     9  zmc(nCv), z11(nCv),                                            !som 13.11.20
     1  z12(nCv), z22(nCv),
     2  zf1(nCv), zf2(nCv),
     3  rho_b_rel(nCv), resmt(nCv)

      real (kind=R8) ::
     *  smbfria(nCv,0:3), smbfrea(nCv,0:3),
     *  smbtfia(nCv,0:3), smbtfea(nCv,0:3)                             !srv 13.01.17 05.07.17
#ifndef FIXED_POINT
      real (kind=R8) ::
     *  ubf(nFc), robf(nFc)
#endif
      logical new_matrix, main_lightest, warning
*   ..procedures
      external xertst, ipgetr, sfill, b2scopy, b2sscal, b2saxpy
      external
     &  b2xvsg, b2xpne, b2xprz, b2xpro, b2xppr, b2xppz,
     &  b2xppb, b2nxdu, b2nxdv, b2nxfx, b2nxfc, b2nxst,
     &  b2sigp, b2siav, b2sicf, b2nxcm, b2urmo, b2usmo

*-----------------------------------------------------------------------
*.computation

* ..preliminaries
*   ..subprogram start-up calls
      call subini ('b2npmo')
*   ..test nCv, nFc, ns
      call xertst (0.lt.nCv.and.0.lt.nFc, 'faulty argument nCv, nFc')
      call xertst (1.le.ns, 'faulty argument ns')
*   ..test itcnt
      call xertst (0.le.itcnt, 'faulty argument itcnt')
*   ..test rxf
      call xertst (0.0_R8.le.rxf.and.rxf.le.1.0_R8,
     . 'faulty argument rxf')
*   ..set internal parameters on first call
      if (ncall_b2npmo.eq.0) then
        write(*,*) 'nsmin, nsmax = ', switch%nsmin, switch%nsmax
        main_lightest = .true.
        do is = 0, ns-1
          main_lightest = main_lightest .and.
     .     nint(am(is)).ge.nint(am(ismain))
        end do
        if (.not.main_lightest .and. switch%b2sigp_style.eq.2 .and.
     .      switch%Zhdanov_closure.eq.0) then
         write(0,'(2a)') 'Warning: b2sigp_style.eq.2 not appropriate ',
     .    'for plasmas where the main species is not the lightest one!'
         write(6,'(2a)') 'Warning: b2sigp_style.eq.2 not appropriate ',
     .    'for plasmas where the main species is not the lightest one!'
        endif
c time-dependent treatment
        if(switch%b2mndt_style.eq.2) then
         rxg_npmo = huge(1.0_R8)
         do is = 0,ns-1
          warning = .false.
          do iCv = 1, nCv
            warning = warning .or.
     .       rxf*dtmo(is,mpg%cvReg(iCv))*time_factor(iCv).ne.1.0_R8
          enddo
          if (warning) then
           write(6,'(a,i3,a)')
     .      'rxf*dtmo*time_factor for species ',is,' should be 1 !'
           write(6,*) 'Time derivative terms may be inacurrate !'
          end if
         enddo
        else
         rxg_npmo = switch%b2npmo_rxg
        endif
        call xertst (0.0_R8.lt.rxg_npmo,
     &    'faulty internal parameter rxg')
      endif
*   ..extensive tests on first few calls
      if (ncall_b2npmo.lt.3) then
*    ..test sign of na, te, ti
       do is = 0, ns-1
        call b2xvsg (nCv, pl%na(:,is), 1, 'na', '.gt.')
       enddo
       call b2xvsg (nCv, pl%te, 1, 'te', '.gt.')
       call b2xvsg (nCv, pl%ti, 1, 'ti', '.gt.')
       call b2xvsg (nCv, pl%tn, 1, 'tn', '.gt.')
*    ..test sign of cvsa
       do is = 0, ns-1
        wrkf0(:) = co%cvsa(:,0,is)*geo%fcQalf(:,0)
        wrkf1(:) = co%cvsa(:,1,is)*geo%fcQalf(:,1)
        call b2xvsg (nFc, wrkf0, 1, 'cvsa0', '.ge.')
        call b2xvsg (nFc, wrkf1, 1, 'cvsa1', '.ge.')
       enddo
*    ..test sign of smo, smq
       do is = 0, ns-1
        call b2xvsg (nCv, sr%smo(:,1,is), 1, 'smo 1', '.le.')
        call b2xvsg (nCv, sr%smo(:,3,is), 1, 'smo 3', '.le.')
        call b2xvsg (nCv, sr%smq(:,1,is), 1, 'smq 1', '.le.')
        call b2xvsg (nCv, sr%smq(:,3,is), 1, 'smq 3', '.le.')
       enddo
      endif

*   ..compute ne, ne2, ro, rz, pr, pz
      call b2xpni (nCv, ns, pl%na, dv%ni)
      call b2xpnn (nCv, ns, pl%na, dv%nn)
      call b2xpne (nCv, ns, rt%rza, pl%na, st_ext%ne, dv%ne)
      call b2xpne (nCv, ns, rt%rz2, pl%na, st_ext%ne2, dv%ne2)
      call b2xpro (nCv, ns, mp, am, pl%na, ro, st_ext)
      call b2xprz (nCv, ns, mp, am, pl%na, rz, st_ext)
      call b2xppr (nCv, ns, dv%ne, pl%na, pl%te, pl%ti, pl%tn, pr,
     &             st_ext)
      call b2xppz (nCv, ns, dv%ne, pl%na, pl%te, pl%ti, pz, st_ext)
      call b2xzef (nCv, ns, rt%rz2, pl%na, dv%ne, zeff, st_ext)        !srv 13.01.17 {
      do iCv = 1, nCv
        ce1(iCv)=fce1(zeff(iCv))
        ce2(iCv)=fce2(zeff(iCv))
      enddo                                                            !srv 13.01.17 }
      if (switch%Zhdanov_closure.eq.0.or.switch%Zhdanov_test.eq.1 .or.
     &    switch%Zhdanov_tf_fr.eq.0) then
        mm(0:ns-1) = switch%mass_multiplicator
        mm(ismain) = 1.0_R8
        am1(0:ns-1)= mm(0:ns-1)*am(0:ns-1)
        call b2xzmf (nCv, ns, ismain, rt%rz2, pl%na, zmff)
        call b2xzmc (nCv, ns, ismain, rt%rz2, pl%na, am1, zmc)
        call b2xzf1 (nCv, ns, ismain, rt%rz2, pl%na, am1, zf1)
        call b2xzf2 (nCv, ns, ismain, rt%rz2, pl%na, am1, zf2)
        call b2xz11 (nCv, ns, ismain, rt%rz2, pl%na, am1, z11)
        call b2xz12 (nCv, ns, ismain, rt%rz2, pl%na, am1, z12)
        call b2xz22 (nCv, ns, ismain, rt%rz2, pl%na, am1, z22)
        if (switch%impr_form_fr.eq.0) then
          do iCv = 1, nCv
            do is = 0, ns-1
              cimp1(iCv,is)=fce1(zmff(iCv))
            end do
          enddo
        else
          do is = 0, ns-1
            do iCv = 1, nCv
              cimp1(iCv,is)=fcim1(ismain,is,ns,am1,z11(iCv),
     &          z12(iCv),z22(iCv),zf1(iCv),zf2(iCv))
            enddo
          enddo
        endif
        if (switch%impr_form_tf.eq.0) then
          do iCv = 1, nCv
            do is = 0, ns-1
              cimp2(iCv,is)=fce2(zmff(iCv))
            enddo
          enddo
        else
          do is = 0, ns-1
            do iCv = 1, nCv
              cimp2(iCv,is)=fcim2(ismain,is,ns,am1,zmc(iCv),
     &          z11(iCv),z12(iCv),z22(iCv))
            enddo
          enddo
        endif
      endif
*   ..compute the Coulomb logarithm
      call b2tlnl (nCv, switch, switch%icase_ii, pl%te, pl%ti,
     &             dv%ne, dv%lnlam)
*   ..compute zetap and zetae
      t1 = (4.0_R8*pi*eps0/qe)/qe
      t0 = (3.0_R8/(4.0_R8*sqrt(2.0_R8*pi)))*t1*(t1*ev*sqrt(ev))
      czetap = t0*sqrt(mp)
      czetae = t0*sqrt(me)
      zetap = czetap*(pl%ti/ev)*sqrt(pl%ti/ev)/dv%lnlam
      zetae = czetae*(pl%te/ev)*sqrt(pl%te/ev)/dv%lnlam
*   ..compute gti, gtn, gte
      !wdk gti, gtn and gte are now computed directly in cell centers
      !wdk (instead of faces as in structured code)
      call gradc_p(nCv, nFc, nVx, 0, geo, mpg, pl%ti, wrkv, gti)
      call gradc_p(nCv, nFc, nVx, 0, geo, mpg, pl%tn, wrkv, gtn)
      call gradc_p(nCv, nFc, nVx, 0, geo, mpg, pl%te, wrkv, gte)
*   ..compute ehxp
      if(switch%b2sigp_style.ne.2) then
        call b2nxfx (nCv, nFc, nVx, switch, geo, mpg, qe, dv%ne,
     &     pl%te, pl%po, ehxp)
      endif                                                            !srv 13.01.17 }
*   ..compute smqdu
      call b2nxdu (nCv, ns, pl%na, pl%ua, sr%smq, smqdu, wrk0, wrk1)
      if (switch%b2npmo_iout.ne.0) then
       call my_out_us(70,nCv,0,smqdu,'b2npmo_smqdu')
      endif
*   ..initialise
      call sfill (nCv*4, 0.0_R8, smt, 1)
      call sfill (nFc*2, 0.0_R8, flct, 1)
      call sfill (nFc*2, 0.0_R8, cvct, 1)
      call sfill (nCv, 0.0_R8, dv%resmt, 1)
      call sfill (nCv, 0.0_R8, resmt, 1)
      call sfill (nCv, 0.0_R8, ctcft, 1)                               !srv 06.04.07

* ..main computation
      lnlam = dv%lnlam ! copy
*   ..compute individual species corrections
#if ( defined(_OPENMP) && !defined(NO_OPENMP_B2NPMO) )
#ifdef WITH_FTIMINGS
      call timer%start("parallel")
#endif
!!FR smt, resmt, cvct & flct accumulate by species and thus need
!!FR to be reduction variables.
!!FR rob, pb, flub, flcb, cvcb, wrkf, smb are all
!!FR local to the parallel region and thus are private
!!FR wrk0 is used outside parallel region but is zeroed before use
!!FR aa5, wrk1 are used outside parallel region but are output arguments
!!FR and thus their contents should not matter
!!FR all other variables are shared
!$omp parallel do default(none)
!$omp& shared(geo, nCv, nFc, nVx, switch, mpg, psnl)
!$omp& shared(co, co_ns, dv, pl, rt, srw, sr)
!$omp& shared(st_ext, ehxp)
!$omp& reduction(+:resmt, flct, cvct, smt)
!$omp& reduction(+:ctcft)
!$omp& private(rob, pb, wrk0, flub, flcb, b2npmo_solvedum)
!$omp& private(pb1, rob1, wrkf0)
!$omp& firstprivate(lnlam)
!$omp& private(is, k)
!$omp& lastprivate(kt_neo)
!$omp& private(wrkf, wrk1, smb, ro)
!$omp& private(smbfrea, smbfria, smbtfea, smbtfia)
!$omp& private(smbch, smbtf, smbvtmp)
!$omp& private(smbgp, flubv, wrk2, smbvh,smbvv)
!$omp& private(smbcf, ctcfb, smban, cvcb)
!$omp& private(smbvi, rho_b_rel)
!$omp& private(chns, wrks, aa, wrkff, wrkfff)
!$omp& firstprivate(last_solve_9) lastprivate(last_solve_9)
!$omp& shared(solvemt, solvemo, ns, am, solving, ismain, smqdu, itcnt)
#ifndef FIXED_POINT
!$omp& shared(balance_netcdf, fmo_b2nxfv)
!$omp& shared(fmo_tot, tot_smo0to3)
!$omp& private(robf, ubf)
#endif
!$omp& shared(is_neutral)
!$omp& shared(zetap, zetae, gti, gte, ce1, ce2, cimp1)
!$omp& shared(cimp2, zeff, z2n_cv, rxg_npmo, dtim)
!$omp& shared(c_hTA_an_fl_save, c_R_TA_an_save)
!$omp& shared(rho_A_rel, avm_u, ts_factor)
#endif
!      do is = 0, ns-1
      do is = switch%nsmin, switch%nsmax-1
*    ..compute rob, pb
       call b2scopy (nCv, pl%na(:,is), 1, rob, 1)
       call b2sscal (nCv, am(is)*mp, rob, 1)
       call b2xppb (nCv, rt%rza(:,is), pl%na(:,is), pl%te, pl%ti,
     &               pl%tn, is, pb)
*    ..compute extended friction term
       if(switch%b2sigp_style.eq.2) then                               !srv 13.01.17
         if(switch%Zhdanov_nc.eq.1 .and. switch%Zhdanov_closure.eq.1)
     >    call b2mod_frtf_nccorr_coefs (nCv, mpg%nFt,
     &     mpg, geo, pl, co, rt, switch,
     &     c_hTA_an_fl_save, c_R_TA_an_save, z2n_cv, zeff)
          call b2sifrtf (nCv, nFc, nVx, ns, is, ismain, switch, geo,
     &   mpg, pl%na, pl%ua, dv%ue, pl%te, pl%ti, dv%ni, dv%ne, dv%ne2,
     &   rt%rza, rt%rz2, pl%po,zetap, zetae, gti, gte, ce1, ce2,
     &   cimp1, cimp2, zeff, co%f_luc_sg, co%alfx_c, co%sigx_c, st_ext,
     &   smbfrea, smbfria, smbtfea, smbtfia, smbch,
     &   smbtf)             ! smbch and smbtf will be multiplied by vol*hz below !srv 05.07.17
        do k = 0, 3
         smbch(:,k) = geo%cvVol*geo%cvHz*smbch(:,k)
         smbtf(:,k) = geo%cvVol*geo%cvHz*smbtf(:,k)
        enddo
        srw%smfr(:,is) = geo%cvVol*geo%cvHz*(
     &   (smbfrea(:,0) + smbfria(:,0)) +
     &   (smbfrea(:,1) + smbfria(:,1))*pl%ua(:,is) +
     &   (smbfrea(:,2) + smbfria(:,2))*am(is)*mp*pl%na(:,is) +
     &   (smbfrea(:,3) + smbfria(:,3))*am(is)*mp*pl%na(:,is)*pl%ua(:,is)
     &   )
       else
        smbfria = 0.0_R8
        smbfrea = 0.0_R8
        smbtfia = 0.0_R8
        smbtfea = 0.0_R8
        call b2sifr_ (nCv, nFc, nVx, ns, is, switch, geo, mpg,
     &   pl%na, pl%ua, pl%te, pl%ti, dv%ni, dv%ne, dv%ne2, rt%rza,
     &   rt%rz2, ehxp, lnlam, co%cthe, co%cthi, st_ext,
     &   smbch, smbtf, srw%smpr(:,is), srw%smpt(:,is), srw%smfr(:,is),
     &   wrk0, wrk1, wrk2)
       endif
*    ..compute pressure gradient term
       call b2sigp (nCv, nFc, nVx, is, switch, geo, mpg, rt%rza(:,is),
     &   pl%na(:,is), pl%ti, pl%tn, dv%ne, pl%te, pl%po, smbgp, wrk0,
     &   wrk1)
*    ..compute flub
       call b2scopy (nFc*2, dv%fna_fcor(:,:,is), 1, flub, 1)
       call b2sscal (nFc*2, am(is)*mp, flub, 1)
*    ..add convective term due to new form of ion viscosity term only
       if (switch%Zhdanov_closure.ne.0 .and. ! can be put outside of loop?
     &      switch%Zhdanov_test.ne.1 .and.
     &      switch%Zhdanov_visu.ne.0 .and.
     &      switch%zh_vis_type.eq.1) then
         rho_b_rel = rho_A_rel(:,is)
       else
         rho_b_rel = 0.0_R8
       endif
       call b2nxfv (nCv, nFc, nVx, is, ismain, switch, geo, mpg,
     &    rho_b_rel, co%cvsahz_cl(:,:,is), co%cvsahz(:,:,is), flubv)
#ifndef FIXED_POINT
       if (balance_netcdf.ne.0) then
         fmo_b2nxfv(:,:,is)=flubv                                       !djm Jan2017
       endif
#endif
       call b2saxpy (nFc*2, 1.0_R8, flubv, 1, flub, 1)                  !srv 11.09.09
*    ..compute flcb, cvcb
       call b2nxfc (nCv, nFc, nVx, is, switch, geo, mpg, rob,
     &  pl%ua(:,is), flub, co%cvsahz(:,:,is), flcb, cvcb,
     &  wrkf0)

*    ..compute additional viscosity term
       if (switch%Zhdanov_closure.eq.0.or.switch%Zhdanov_test.eq.1 .or.
     &    (switch%Zhdanov_visu.eq.0.and.switch%Zhdanov_vish.eq.0)) then !som 13.08.21
         call b2siav (nCv, nFc, nVx, ns, is, ismain, switch, geo, mpg,
     &    am(is)*mp, pl%ua(:,is), pl%te, pl%ti, pl%na(:,is),
     &    dv%ni(:,1),
     &    dv%ne2, rt%rz2, lnlam, co%hcix_c, co%cvsahz_cl(:,:,is),
     &    smbvh, smbvv, kt_neo, wrk0, wrk1,
     &    wrk2)
          smbvi = 0.0_R8
       elseif (switch%Zhdanov_visu.eq.1 .and.
     &         switch%Zhdanov_vish.eq.1) then
         call b2siav_zh (nCv, nFc, nVx, ns, is, switch, geo, mpg,
     &    am(is)*mp, rt%rza(:,is), pl%te, pl%ti, pl%na, pl%ua,
     &    dv%ne2, rt%rz2, dv%lnlam, co%hci_a,
     &    co_ns%hci_al_ast, co%cvsahz_cl(:,:,is),
     &    co%cvsahz_drho(:,:,is), co%vsaf_drho(:,:,is),
     &    co_ns%vsaf_hAdp_albe, co_ns%vsaf_hBdp_al,
     &    co_ns%vsaf_uAdp_albe, co_ns%vsaf_uBdp_al,
     &    smbvh, smbvv, smbvi, rho_A_rel(:,is), avm_u, wrk0, wrk1, wrk2)
       elseif (switch%Zhdanov_visu.eq.1 .and.
     &         switch%Zhdanov_vish.eq.0) then
         call b2siav (nCv, nFc, nVx, ns, is, ismain, switch, geo, mpg,
     &    am(is)*mp, pl%ua(:,is), pl%te, pl%ti, pl%na(:,is),
     &    dv%ni(:,1),
     &    dv%ne2, rt%rz2, lnlam, co%hcix_c, co%cvsahz_cl(:,:,is),
     &    smbvh, smbvtmp, kt_neo, wrk0, wrk1,
     &    wrk2)
         call b2siav_zh (nCv, nFc, nVx, ns, is, switch, geo, mpg,
     &    am(is)*mp, rt%rza(:,is), pl%te, pl%ti, pl%na, pl%ua,
     &    dv%ne2, rt%rz2, dv%lnlam, co%hci_a,
     &    co_ns%hci_al_ast, co%cvsahz_cl(:,:,is),
     &    co%cvsahz_drho(:,:,is), co%vsaf_drho(:,:,is),
     &    co_ns%vsaf_hAdp_albe, co_ns%vsaf_hBdp_al,
     &    co_ns%vsaf_uAdp_albe, co_ns%vsaf_uBdp_al,
     &    smbvtmp, smbvv, smbvi, rho_A_rel(:,is), avm_u,
     &    wrk0, wrk1, wrk2)
       elseif (switch%Zhdanov_visu.eq.0 .and.
     &         switch%Zhdanov_vish.eq.1) then
         call b2siav (nCv, nFc, nVx, ns, is, ismain, switch, geo, mpg,
     &    am(is)*mp, pl%ua(:,is), pl%te, pl%ti, pl%na(:,is),
     &    dv%ni(:,1),
     &    dv%ne2, rt%rz2, lnlam, co%hcix_c, co%cvsahz_cl(:,:,is),
     &    smbvtmp, smbvv, kt_neo, wrk0, wrk1,
     &    wrk2)
         call b2siav_zh (nCv, nFc, nVx, ns, is, switch, geo, mpg,
     &    am(is)*mp, rt%rza(:,is), pl%te, pl%ti, pl%na, pl%ua,
     &    dv%ne2, rt%rz2, dv%lnlam, co%hci_a,
     &    co_ns%hci_al_ast, co%cvsahz_cl(:,:,is),
     &    co%cvsahz_drho(:,:,is), co%vsaf_drho(:,:,is),
     &    co_ns%vsaf_hAdp_albe, co_ns%vsaf_hBdp_al,
     &    co_ns%vsaf_uAdp_albe, co_ns%vsaf_uBdp_al,
     &    smbvh, smbvtmp, smbvi, rho_A_rel(:,is), avm_u,
     &    wrk0, wrk1, wrk2)
        smbvi = 0.0_R8
       else
        kt_neo = 0.0_R8
       endif
       smbvi = smbvi*switch%smbvi_factor
*    ..compute centrifugal force term                                   !srv 09.07.01
       call b2sicf (nCv, nFc, nVx, is, am(is)*mp, switch,  geo, mpg,    !srv 09.07.01 28.01.02 22.07.05
     *   pl%na(:,is), pl%ua(:,is), smbcf, ctcfb, wrk0)
       srw%smcf(:,is)=smbcf(:,0)+smbcf(:,1)*pl%ua(:,is)+
     &                smbcf(:,2)*pl%na(:,is)*mp*am(is)+
     &                smbcf(:,3)*pl%ua(:,is)*pl%na(:,is)*mp*am(is)
*    ..compute correcting source due to artificial anomalous electric conductivity
       call b2sian (nCv, nFc, is, switch, geo, mpg,
     &   dv%fchanml_a(:,:,is), smban)
*    ..compute total parallel momentum source
       call b2scopy (nCv*4, sr%smo(:,:,is), 1, smb, 1)
       call b2saxpy (nCv*4, 1.0_R8, sr%smq(:,:,is), 1, smb, 1)
       srw%b2npmo_smaf(:,:,is)=smbch(:,:)+smbtf(:,:)                    !srv 11.09.09 {
       call b2saxpy (nCv*4, 1.0_R8, smbch, 1, smb, 1)
       call b2saxpy (nCv*4, 1.0_R8, smbtf, 1, smb, 1)                   !srv 11.09.09 }
       srw%b2npmo_smag(:,:,is)=smbgp(:,:)
       call b2saxpy (nCv*4, 1.0_R8, smbgp, 1, smb, 1)
       srw%b2npmo_smav(:,:,is)=smbvh(:,:)+smbvv(:,:)+smbvi(:,:)        !srv 03.07.08 {
       call b2saxpy (nCv*4, 1.0_R8, smbvh,  1, smb, 1)
       call b2saxpy (nCv*4, 1.0_R8, smbvv,  1, smb, 1)
       call b2saxpy (nCv*4, 1.0_R8, smbvi,  1, smb, 1)
       call b2saxpy (nCv*4, 1.0_R8, smbcf,  1, smb, 1)                  !srv 03.07.08 }
       call b2saxpy (nCv*4, 1.0_R8, smban,  1, smb, 1)                  !srv 25.09.17
       call b2nxdv (nCv,is,switch,mpg,rob,sr%smq(:,:,is),smqdu,smb)     !srv
#ifndef FIXED_POINT
cdjm Jan2017 update momentum source components for balance
       if (balance_netcdf.ne.0) then
         tot_smo0to3(:,:,is)=smb !djm Mar2019
        call intface(nCv,nFc,mpg%fcCv,geo%fcVol,pl%ua(:,is),ubf)
        call intface(nCv,nFc,mpg%fcCv,geo%fcVol,rob,robf)
        call update_balance_sources_momentum(nCv,nFc,is,pl%ua(:,is),rob,
     &    ubf,robf)
       endif
       if (balance_netcdf.ne.0 .or.
     &    switch%write_eirene_sources.ne.0) then
        call update_balance_eirene_sources_momentum(
     &   nCv,is,pl%ua(1,is),rob)
       endif
#endif

*    ..time-dependent treatment
*    ..compute residual without time-stepping sources (needed at next time-step for ckn scheme)
       if (switch%b2mndt_ckn.eq.2) then
         if((solving.and.any(solvemo(is,0:mpg%nnreg(0))))
     &            .or.switch%get_residuals.eq.1) then
           call b2urmo (nCv, nFc, nVx, switch, geo, mpg, is,
     &       pl%ua(:,is), rob, smb, flub, co%cvsahz(:,:,is),
     &       dv%resmo0(:,is), dv%fmo(:,:,is), dv%fhm(:,:,is),
     &       co%cvsahz_eff(:,:,is))
         endif
       endif
*    .. calculate time-dependent source
*    .. ts_factor = time-stepping multiplication factor
       if (switch%b2mndt_style.eq.2) then
         call b2smdt (nCv, mpg%nCi, ns, dtim, geo%cvHz, geo%cvVol,
     &                pl%na, psnl%ua, switch, mpg, sr%smodt, is)
*    .. add contribution from time-dependent source
         sr%smodt(:,0,is) = sr%smodt(:,0,is)-ts_factor*
     &    am(is)*mp*geo%cvVol*geo%cvHz*(pl%na(:,is)-psnl%na(:,is))*
     &    psnl%ua(:,is)/dtim
         smb = smb + sr%smodt(:,:,is)
       endif

*    ..compute residual
      if((solving.and.any(solvemo(is,0:mpg%nnreg(0))))
     &            .or.switch%get_residuals.eq.1) then                   !srv 22.05.18
        call b2urmo (nCv, nFc, nVx, switch, geo, mpg, is, pl%ua(:,is),
     &    rob, smb, flub, co%cvsahz(:,:,is), dv%resmo(:,is),
     &    dv%fmo(:,:,is), dv%fhm(:,:,is), co%cvsahz_eff(:,:,is))        !srv 04.07.08
      endif                                                             !srv 22.05.18
*    ..modify residual for ckn scheme
      if (switch%b2mndt_ckn.eq.2) then
        if((solving.and.any(solvemo(is,0:mpg%nnreg(0))))
     &      .or.switch%get_residuals.eq.1) then
          dv%resmo(1:mpg%nCi,is) = dv%resmo(1:mpg%nCi,is) +
     &     psnl%resmo0(1:mpg%nCi,is)
        endif
      elseif (switch%b2mndt_bdf.eq.2) then
        if((solving.and.any(solvemo(is,0:mpg%nnreg(0))))
     &     .or.switch%get_residuals.eq.1) then
          dv%resmo(1:mpg%nCi,is) = dv%resmo(1:mpg%nCi,is) +
     &     0.5_R8*psnl%dmodt(1:mpg%nCi,is) ! +1/2*m*V*[d(nu)/dt]_i (where i denotes the previous time-step)
        endif
      endif

#ifndef FIXED_POINT
cdjm Mar2019 output the total momentum flux and source for balance
      if (balance_netcdf.ne.0) then
        fmo_tot(:,:,is)=dv%fmo(:,:,is)
      endif
#endif

      write (chns,'(i3.3)') is                                          !srv 18.03.10

      if (switch%b2npmo_iout.ne.0 .or. switch%iout_b2wdat .eq. 4) then  !srv 12.07.02 {
        call my_out_us(70,nCv,0,dv%resmo(1,is),'b2npmo_resmo'//chns)
        call my_out_us(70,nFc,1,dv%fmo(1,0,is),'b2npmo_fmo_th'//chns)
        call my_out_us(70,nFc,1,dv%fmo(1,1,is),'b2npmo_fmo_r'//chns)
        if(.not. is_neutral(is)) then
          call my_out_us(70,nFc,1,flubv(1,0),'b2npmo_flubv_th'//chns)
          wrks=smbvv(:,0)+smbvv(:,1)*pl%ua(:,is)+
     &         smbvv(:,2)*pl%na(:,is)*mp*am(is)+
     &         smbvv(:,3)*pl%ua(:,is)*pl%na(:,is)*mp*am(is)
          call my_out_us(70,nCv,0,wrks,'b2npmo_smovv'//chns)
          wrks=smbvi(:,0)+smbvi(:,1)*pl%ua(:,is)+
     &         smbvi(:,2)*pl%na(:,is)*mp*am(is)+
     &         smbvi(:,3)*pl%ua(:,is)*pl%na(:,is)*mp*am(is)
          call my_out_us(70,nCv,0,wrks,'b2npmo_smovi'//chns)
          wrks=smbvh(:,0)+smbvh(:,1)*pl%ua(:,is)+
     &         smbvh(:,2)*pl%na(:,is)*mp*am(is)+
     &         smbvh(:,3)*pl%ua(:,is)*pl%na(:,is)*mp*am(is)
          call my_out_us(70,nCv,0,wrks,'b2npmo_smovh'//chns)
        endif
        wrks=smb(:,0)+smb(:,1)*pl%ua(:,is)+
     &       smb(:,2)*pl%na(:,is)*mp*am(is)+
     &       smb(:,3)*pl%ua(:,is)*pl%na(:,is)*mp*am(is)
        call my_out_us(70,nCv,0,wrks,'b2npmo_smb'//chns)                !srv 28.02.17
        call my_out_us(70,nCv,0,srw%smcf(1,is),'b2npmo_smocf'//chns)
        wrks = geo%cvVol*geo%cvHz*(smbfria(:,0)+
     &        smbfria(:,1)*pl%ua(:,is)+      !srv 05.07.17 {
     &        smbfria(:,2)*pl%na(:,is)*mp*am(is)+
     &        smbfria(:,3)*pl%ua(:,is)*pl%na(:,is)*mp*am(is))
        call my_out_us(70,nCv,0,wrks,'b2npmo_smofria'//chns)
        wrks = geo%cvVol*geo%cvHz*(smbfrea(:,0)+
     &                             smbfrea(:,1)*pl%ua(:,is))
        call my_out_us(70,nCv,0,wrks,'b2npmo_smofrea'//chns)
        wrks = geo%cvVol*geo%cvHz*smbtfia(:,0)
        call my_out_us(70,nCv,0,wrks,'b2npmo_smotfia'//chns)
        wrks = geo%cvVol*geo%cvHz*smbtfea(:,0)
        call my_out_us(70,nCv,0,wrks,'b2npmo_smotfea'//chns)
        call my_out_us(70,nCv,0,dv%ue,'b2npmo_ue')                      !srv 05.07.17 }
        wrks=smban(:,0)+smban(:,1)*pl%ua(:,is)+
     &       smban(:,2)*pl%na(:,is)*mp*am(is)+
     &       smban(:,3)*pl%ua(:,is)*pl%na(:,is)*mp*am(is)
        call my_out_us(70,nCv,0,wrks,'b2npmo_smoan'//chns)              !srv 25.09.17
        wrks=sr%smodt(:,0,is)+
     &       sr%smodt(:,1,is)*pl%ua(:,is)+
     &       sr%smodt(:,2,is)*pl%na(:,is)*mp*am(is)+
     &       sr%smodt(:,3,is)*pl%ua(:,is)*pl%na(:,is)*mp*am(is)
         call my_out_us(70,nCv,0,wrks,'b2npmo_smodt'//chns)
      endif                                                             !srv 05.07.17

      if (switch%b2npmo_iout .ne. 0 .or. switch%iout_b2wdat .eq. 4) then
        call my_out_us(70,nCv,0,pl%ua(1,is),'b2npmo_ua_input'//chns)
        wrks=sr%smo(:,0,is)+sr%smo(:,1,is)*pl%ua(:,is)+
     &       sr%smo(:,2,is)*pl%na(:,is)*mp*am(is)+
     &       sr%smo(:,3,is)*pl%ua(:,is)*pl%na(:,is)*mp*am(is)
        call my_out_us(70,nCv,0,wrks,'b2npmo_smo'//chns)                !srv 28.02.17
        wrks=smbch(:,0)+smbch(:,1)*pl%ua(:,is)+
     &       smbch(:,2)*pl%na(:,is)*mp*am(is)+
     &       smbch(:,3)*pl%ua(:,is)*pl%na(:,is)*mp*am(is)
        call my_out_us(70,nCv,0,wrks,'b2npmo_smoch'//chns)
        wrks=smbtf(:,0)+smbtf(:,1)*pl%ua(:,is)+
     &       smbtf(:,2)*pl%na(:,is)*mp*am(is)+
     &       smbtf(:,3)*pl%ua(:,is)*pl%na(:,is)*mp*am(is)
        call my_out_us(70,nCv,0,wrks,'b2npmo_smotf'//chns)
        wrks=smbgp(:,0)+smbgp(:,1)*pl%ua(:,is)+
     &       smbgp(:,2)*pl%na(:,is)*mp*am(is)+
     &       smbgp(:,3)*pl%ua(:,is)*pl%na(:,is)*mp*am(is)
        call my_out_us(70,nCv,0,wrks,'b2npmo_smogp'//chns)
       endif                                                            !srv 12.07.02 }

       if (solving.and.any(solvemo(is,0:mpg%nnreg(0)))) then            !srv 06.07.10 22.05.18
         b2npmo_solvedum(0:mpg%nnreg(0))=solvemo(is,0:mpg%nnreg(0))
c    ..compare matrices
#if ( defined(_OPENMP) && !defined(NO_OPENMP_B2NPMO) )
         call restart_ma28_for_us
#else
         new_matrix = .false.
         do ireg = 0, mpg%nnreg(0)
          new_matrix = new_matrix.or.
     &     solvemo(is,ireg).neqv.last_solve_9(ireg)
         end do
         if (new_matrix.or..true.) call restart_ma28_for_us
#endif
*    ..compute correction
         pb1 = pb !csc to avoid aliasing in Tapenade
         rob1 = rob !csc to avoid aliasing in Tapenade
         call b2usmo (nCv, nFc, nVx, switch, geo, mpg, mpg%nnreg(0),
     &    b2npmo_solvedum, itcnt, rxg_npmo, rob, rob1, pb,
     &    pb1, pl%ua(:,is), smb,
     &    flcb, cvcb, dv%resmo(:,is), ctcfb, dv%corua(:,is),
     &    dv%pcca(:,:,is), aa, wrk1, 'b2usmo_'//chns)
         last_solve_9(0:mpg%nnreg(0)) = solvemo(is,0:mpg%nnreg(0))
csw 21feb2014 added this else branch
       else
         dv%pcca(:,:,is)  = 0.0_R8
!         dv%resmo(:,is) = 0.0_R8       ! it is useful to know resmo    !srv 22.05.18
         dv%corua(:,is) = 0.0_R8
       endif

*    ..compute contribution to total residual
csw 21feb2014 and do this only for the solved species, i.e. added next branch
      if((solving.and.any(solvemo(0:ns-1,0:mpg%nnreg(0)))
     &   .and.any(solvemt(0:mpg%nnreg(0))))) then !srv 22.05.18 resmt is non-zero when all equations are solved
!     &   .and.all(solvemt(0:mpg%nnreg(0)))).or.get_residuals.eq.1) then !srv 22.05.18
        call b2scopy (nCv*4, sr%smo(:,:,is), 1, smb, 1)
        if(switch%b2sigp_style.eq.2) then                               !srv 13.01.17 {
*       (in this case thermal force and the friction are taken into account)
         call b2saxpy (nCv*4, 1.0_R8, smbtf, 1, smb, 1)                 !srv 08.10.14
         call b2saxpy (nCv*4, 1.0_R8, smbch, 1, smb, 1)
        else
*       (now the friction cancels out)
        endif                                                           !srv 13.01.17 }
        call b2saxpy (nCv*4, 1.0_R8, smbgp, 1, smb, 1)
        call b2saxpy (nCv*4, 1.0_R8, smbvh, 1, smb, 1)
        call b2saxpy (nCv*4, 1.0_R8, smbvv, 1, smb, 1)                  !srv 03.07.08
        call b2saxpy (nCv*4, 1.0_R8, smbvi, 1, smb, 1)
        call b2saxpy (nCv*4, 1.0_R8, smbcf, 1, smb, 1)                  !srv 09.07.01
        call b2saxpy (nCv*4, 1.0_R8, smban, 1, smb, 1)                  !srv 25.09.17
        call b2urmo (nCv, nFc, nVx, switch, geo, mpg,  is, pl%ua(:,is),
     .     rob, smb, flub, co%cvsahz(:,:,is), wrk0, wrkf, wrkff, wrkfff)
*    ..accumulate residual
        call b2saxpy (nCv, 1.0_R8, wrk0, 1, resmt, 1)
*    ..accumulate smb
        call b2nxst (nCv, rob, smb, ro, smt)
*    ..accumulate transport coefficients
        call b2saxpy (nFc*2, 1.0_R8, flcb, 1, flct, 1)
        call b2saxpy (nFc*2, 1.0_R8, cvcb, 1, cvct, 1)
        call b2saxpy (nCv, 1.0_R8, ctcfb, 1, ctcft, 1)                  !srv 06.04.07
       endif                                                            !srv 22.05.18
      enddo
#if ( defined(_OPENMP) && !defined(NO_OPENMP_B2NPMO) )
!$omp end parallel do
#endif

*     .. attribute kt_neo, resmt
      dv%resmt  = resmt
      co%kt_neo = kt_neo

csw 21feb2014, execute simultaneous corrections always
      if (solving.and.any(solvemo(0:ns-1,0:mpg%nnreg(0)))
     &                  .and.any(solvemt(0:mpg%nnreg(0)))) then         !srv 06.07.10 22.05.18
*     ..compute simultaneous correction
        call sfill (nCv, 0.0_R8, wrk0, 1)
c     ..compare 9-pt stencil matrix
#if ( defined(_OPENMP) && !defined(NO_OPENMP_B2NPMO) )
        call restart_ma28_for_us
#else
        new_matrix = .false.
        do ireg = 0, mpg%nnreg(0)
         new_matrix = new_matrix.or.
     &    solvemt(ireg).neqv.last_solve_9(ireg)
        end do
        if (new_matrix.or..true.) call restart_ma28_for_us
#endif
*     ..compute correction
        call b2usmo (nCv, nFc, nVx, switch, geo, mpg, mpg%nnreg(0),
     &    solvemt, itcnt, rxg_npmo, ro, rz, pr, pz, wrk0,
     &    smt, flct, cvct,
     &    dv%resmt, ctcft, dv%corut, dv%pccm, aa, wrk1, 'b2usmo_total')
        last_solve_9(0:mpg%nnreg(0)) = solvemt(0:mpg%nnreg(0))
      else
        dv%pccm  = 0.0_R8
!        dv%resmt = 0.0_R8            ! it is useful to know resmt      !srv 22.05.18
        dv%corut = 0.0_R8
      endif                                                             !srv 06.07.10

*   ..apply correction
      if (rxf.ne.0.0_R8) then
       do is = 0, ns-1
        do iCv = 1, mpg%nCi ! internal cells
         ireg = mpg%cvReg(iCv)
         if (solvemo(is,ireg)) then
          if(solvemt(ireg)) then                                        !srv 22.05.18
            pl%ua(iCv,is) = pl%ua(iCv,is)+
     &        rxf*med(dv%corua(iCv,is),dv%corut(iCv),
     &        dv%corua(iCv,is)+dv%corut(iCv))
          else                                                          !srv 22.05.18 {
            pl%ua(iCv,is) = pl%ua(iCv,is) + rxf*dv%corua(iCv,is)
          endif                                                         !srv 22.05.18 }
         endif
        enddo
        do iCv = mpg%nCi+1, nCv ! guard cells
         ireg = mpg%cvReg(iCv)
         if (solvemo(is,ireg)) then
           pl%ua(iCv,is) = pl%ua(iCv,is) + rxf*dv%corua(iCv,is)
         endif
        enddo
       enddo
      endif

      vcount = 0
      if( switch%ion_vlct_restrict .ne. 0 ) then                        !iyv 19.07.12 { !srv 17.01.14
        do is = 0, ns-1
          do iCv = 1, nCv
            if (.not.is_neutral(is).or.switch%tn_style.le.1.or.
     .         (switch%tn_style.eq.2.and.nint(zn(is)).ne.1)) then
              cst=sqrt(pz(iCv)/rz(iCv))                                 !iyv 03.02.14
            else
              cst=sqrt(pl%tn(iCv)/(mp*am(is)))
            endif
            if( abs(pl%ua(iCv,is)).gt.
     .         switch%ion_vlct_restrict_M*cst ) then
              vcount(is) = vcount(is) + 1
              if (switch%vlct_diagno.ge.2) then
                write(6,'(a,a,i3,a,i4,a,i4,a)')
     .             'Parallel velocity restriction applied ',
     .             'to species ',is,' in cell (',iCv,'):'
                write(6,'(a,1p,1e12.5,a,1e12.5,a,1e12.5)')
     .             'Mach ',switch%ion_vlct_restrict_M,': ',
     .                     switch%ion_vlct_restrict_M*cst,
     .              ' Ua:',pl%ua(iCv,is)
              end if
              pl%ua(iCv,is)= sign(switch%ion_vlct_restrict_M*cst,
     .                              pl%ua(iCv,is))                      !iyv 03.02.14
            endif
          enddo
        enddo
        if (switch%vlct_diagno.ge.1) then
          do is = 0, ns-1
            if (vcount(is).gt.1) then
              write(6,'(a,a,1e12.5,a,i5,a,i3)')
     .         'In b2npmo, parallel velocity restriction ',
     .         'to Mach ', switch%ion_vlct_restrict_M,
     .         ' was applied ',vcount(is),
     .         ' times for species ',is
            else if (vcount(is).eq.1) then
              write(6,'(a,a,1e12.5,a,i5,a,i3)')
     .         'In b2npmo, parallel velocity restriction ',
     .         'to Mach ', switch%ion_vlct_restrict_M,
     .         ' was applied ',vcount(is),
     .         ' time for species ',is
            end if
          end do
        endif
      endif

!! do the same but with possibly a lower max M for the core region.
      vcount = 0
      if( switch%ion_vlct_restrict_core .ne. 0 ) then
        do is = 0, ns-1
          do iCv = 1, nCv
            if (mpg%cvOnClosedSurface(iCv)) then
              if (.not.is_neutral(is).or.switch%tn_style.le.1.or.
     .           (switch%tn_style.eq.2.and.nint(zn(is)).ne.1)) then
                cst=sqrt(pz(iCv)/rz(iCv))
              else
                cst=sqrt(pl%tn(iCv)/(mp*am(is)))
              endif
              if( abs(pl%ua(iCv,is)).gt.
     .            switch%ion_vlct_restrict_core_M*cst ) then
                vcount(is) = vcount(is) + 1
                if (switch%vlct_diagno.ge.2) then
                  write(6,'(a,a,i3,a,i4,a,i4,a)')
     .             'Core parallel velocity restriction applied ',
     .             'to species ',is,' in cell (',iCv,'):'
                  write(6,'(a,1p,1e12.5,a,1e12.5,a,1e12.5)')
     .             'Mach ',switch%ion_vlct_restrict_core_M,': ',
     .                     switch%ion_vlct_restrict_core_M*cst,
     .              ' Ua:',pl%ua(iCv,is)
                end if
                pl%ua(iCv,is)= sign(switch%ion_vlct_restrict_core_M*cst,
     .                              pl%ua(iCv,is))
              endif
            endif
          enddo
        enddo
        if (switch%vlct_diagno.ge.1) then
          do is = 0, ns-1
            if (vcount(is).gt.1) then
              write(6,'(a,a,1e12.5,a,i5,a,i3)')
     .         'In b2npmo, core parallel velocity restriction ',
     .         'to Mach ', switch%ion_vlct_restrict_core_M,
     .         ' was applied ',vcount(is),
     .         ' times for species ',is
            else if (vcount(is).eq.1) then
              write(6,'(a,a,1e12.5,a,i5,a,i3)')
     .         'In b2npmo, parallel velocity restriction ',
     .         'to Mach ', switch%ion_vlct_restrict_core_M,
     .         ' was applied ',vcount(is),
     .         ' time for species ',is
            end if
          end do
        endif
      endif

!! do the same but with a value for the parallel velocity rather than the Mach number.
!! Limiting only Mach number may not be enough if there is a feedback loop on viscous heating -> increased T -> increased allowable ua and sound speed -> increased ua -> increased viscous heating -> ...
      vcount = 0
      if( switch%ion_vlct_restrict_abs .ne. 0 ) then
        do is = 0, ns-1
          do iCv = 1, nCv
              if( abs(pl%ua(iCv,is)).gt.
     .         switch%ion_vlct_restrict_ua ) then
                vcount(is) = vcount(is) + 1
                if (switch%vlct_diagno.ge.2) then
                  write(6,'(a,a,i3,a,i4,a,i4,a)')
     .             'Parallel velocity restriction applied on abs(ua)',
     .             'to species ',is,' in cell (',iCv,'):'
                  write(6,'(a,1p,1e12.5,a,1e12.5,a,1e12.5)')
     .             'Max ua ',switch%ion_vlct_restrict_ua,
     .              ' Ua:',pl%ua(iCv,is)
                end if
                pl%ua(iCv,is)= sign(switch%ion_vlct_restrict_ua,
     .                              pl%ua(iCv,is))
              endif
          enddo
        enddo
        if (switch%vlct_diagno.ge.1) then
          do is = 0, ns-1
            if (vcount(is).gt.1) then
              write(6,'(a,a,1e12.5,a,i5,a,i3)')
     .         'In b2npmo, parallel velocity restriction ',
     .         'to ua ', switch%ion_vlct_restrict_ua,
     .         ' was applied ',vcount(is),
     .         ' times for species ',is
            else if (vcount(is).eq.1) then
              write(6,'(a,a,1e12.5,a,i5,a,i3)')
     .         'In b2npmo, parallel velocity restriction ',
     .         'to ua ', switch%ion_vlct_restrict_ua,
     .         ' was applied ',vcount(is),
     .         ' time for species ',is
            end if
          end do
        endif
      endif

*   ..compute physical time-derivative after correction has been applied
*   .. m*d(nu)/dt = m*[u*(dn/dt) + n*(du/dt)]
! to obtain the true physical derivative, we need to divide by ts_factor, which was previously multiplied to smodt
      do is = 0, ns-1
        dv%dmodt(:,is) = ( -sr%smodt(:,0,is)                          ! = -ma*mp*V/dt*h_z*n_(i+1)*u_i +ma*mp*V/dt*h_z*u_i*(n_(i+1)-n_i)
     &    -sr%smodt(:,1,is)*pl%ua(:,is)                             ! = 0
     &    -sr%smodt(:,2,is)*pl%na(:,is)*mp*am(is)                   ! = 0
     &    -sr%smodt(:,3,is)*pl%ua(:,is)*pl%na(:,is)*mp*am(is) )      ! = +ma*mp*V/dt*h_z*n_(i+1)*u_(i+1)
     &    /ts_factor
      enddo

*   ..correct pcca
      call b2nxcm (nCv, ns, dv%pccm, dv%pcca)
      if (switch%b2npmo_iout.ne.0 .or. switch%iout_b2wdat .eq. 4) then  !srv 12.07.02 {
       do is = 0, ns-1
        write (chns,'(i3.3)') is                                        !srv 18.05.02
        call my_out_us(70,nCv,0,pl%ua(1,is),'b2npmo_ua'//chns)
        call my_out_us(70,nCv,0,dv%dmodt(1,is),'b2npmo_dmodt'//chns)
       enddo
      endif                                                             !srv 12.07.02 }
      if (switch%b2npmo_iout.ne.0) then
       call my_out_us(70,nCv,0,dv%pccm(1,0),'b2npmo_pccm_th')           !srv 11.01.13
       call my_out_us(70,nCv,0,dv%pccm(1,1),'b2npmo_pccm_r')            !srv 11.01.13
      endif                                                             !srv 12.07.02 }

* ..return
      ncall_b2npmo = ncall_b2npmo+1
      call subend ()
      return
*-----------------------------------------------------------------------
*.end b2npmo

      end subroutine b2npmo
