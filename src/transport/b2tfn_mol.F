*-----------------------------------------------------------------------
*.specification

      subroutine b2tfn_mol (nCv, nFc, nVx,
     &  switch, geo, mpg, pl, dv, co)
      use b2mod_types
      use b2mod_constants
      use b2mod_b2cmpa
      use b2mod_b2cmfs
      use b2mod_switches
      use b2us_geo
      use b2us_map
      use b2us_plasma
      use b2mod_subsys
      use b2mod_fluid_molecules
      implicit none
*   ..input arguments (unchanged on exit)
      integer nCv, nFc, nVx, isb, ismain
      type (switches), intent (in) :: switch
      type (geometry), intent (in) :: geo
      type (mapping), intent (in) :: mpg
      type(B2Plasma), intent (in) :: pl
*   ..input/output arguments
      type(B2Coeff), intent (inout) :: co
      type(B2Derivatives), intent (inout) :: dv
*-----------------------------------------------------------------------
*.documentation
*
*  1. purpose
*
*     B2TFMOL computes the molecular fluxes and the coefficients that 
*     specify its linearisation.
*
*-----------------------------------------------------------------------
*.declarations

*   ..local variables
      integer :: iFc, meth, is_main_ion                                
      character chns*3                                                 
      real (kind=R8) :: nmf(nFc), pm(nCv), fn_mol0(nFc,0:1),
     *   flo(nFc,0:1), term(nFc), flob0(nFc,0:1), conb0(nFc,0:1),
     *   cdpm(nFc,0:1), cdpm_RC(nFc,0:1),  dpbm(nFc), dpm(nFc,0:1),
     *   flo0(nFc,0:1), con0(nFc,0:1), flo1(nFc,0:1), con1(nFc,0:1),
     *   flo_tmp(nFc, 0:1), con_tmp(nFc, 0:1, 0:4),
     *   nif(nFc), fnan(nFc,0:1), vif(nFc,0:1), nm_eq(nCv), nm_eqf(nFc),
     *   uif(nFc)
      real (kind=R8) ::
     *   umf(nFc), wrk(nCv), wrkf(nFc,0:1), wrkVx(nVx), conv(nFc),
     *   dummy_fc0(nFc,0:1), dummy_fc1(nFc,0:1), zeroes_fc0(nFc,0:1),
     *   weight(nFc,2), kmpr(nCv), kmprv(nVx), pit(nFc), pit2(nFc),
     *   coeff_vif(nFc), fn_dp(nFc,0:1)
      real (kind=R8) ::                                     
     *   dpmc(nCv)           !! poloidal pressure gradient at cell center
      real (kind=R8) ::      !! only for postprocessing
     *  ww(nFc), vv(nFc), uu(nFc), 
     *  wwc(nCv), vvc(nCv), uuc(nCv), pitc(nCv), pit2c(nCv)  
      real (kind=R8) ::      !! test
     *  tmf(nFc), dtm(nFc,0:1), dp_mol_f(nFc)              
      real (kind=R8), parameter :: eps = 1.0e-10_R8
*   ..procedures
      external xertst
      external b2xvsg, intface, diff

*-----------------------------------------------------------------------
*.computation

* ..preliminaries
*   ..subprogram start-up calls
      call subini ('b2tfn_mol')
*   ..test nCv, nFc
      call xertst (0.lt.nCv.and.0.lt.nFc, 'faulty argument nCv, nFc')
*   ..set internal parameters on first call

*   ..initializations
      zeroes_fc0 = 0.0_R8
      flo_tmp = 0.0_R8
      con_tmp = 0.0_R8
      vif = 0.0_R8
      fnan = 0.0_R8
      meth = switch%b2tfnb_discr_meth

*   ..interpolate nm and pm to cell faces
      weight = 1.0_R8
      call intface(nCv,nFc,mpg%fcCv,weight,n_mol(:,0),nmf)
      call intface(nCv,nFc,mpg%fcCv,weight,t_mol(:,0),tmf)

*   ..compute differences of pm at cell faces and poloidal gradients
*     at cell centers
      pm = n_mol(:,0)*t_mol(:,0)
      call diff(nCv, nFc, nVx, 0, geo, mpg, pm, wrkVx, dpm)
      call gradc_p(nCv, nFc, nVx, 1, geo, mpg, pm, wrkVx,
     &  dpmc)

*   ..precompute kmpr and its differences
      kmpr = n_mol(:,0) / pm
      call intvertex(nCv,nVx,mpg,geo%vxVol,kmpr,kmprv)

*   ..interpolate um to cell faces
      if(switch%b2tfnb_style_int_vel.eq.0) then
        call intface(nCv,nFc,mpg%fcCv,geo%fcVol,u_mol(:,0),umf)
      else
        call intface(nCv,nFc,mpg%fcCv,geo%fcHc,u_mol(:,0),umf)
      endif

*   ..quantities needed for the ion-convective part in the pressure-diffusion equation
      !! equilibrium density
      nm_eq = n_mol(:,0)*beta_eq
      call intface(nCv,nFc,mpg%fcCv,weight,nm_eq,nm_eqf)
      call intface(nCv,nFc,mpg%fcCv,weight,coeff_vi,coeff_vif)

      !! flux per square meter
      is_main_ion = 1 !! TODO generalize
      fnan(:,0) = dv%fna_eir(:,0,is_main_ion)/                !! poloidal
     .                    (geo%fcS*geo%fcQalf(:,0)+eps)
      fnan(:,1) = dv%fna_eir(:,1,is_main_ion)/                !! radial
     .                    (geo%fcS*geo%fcQalf(:,1)+eps)
      !! interpolate main ion density to faces
      call intface(nCv,nFc,mpg%fcCv,weight,pl%na(:,is_main_ion),nif)
      call intface(nCv,nFc,mpg%fcCv,weight,pl%ua(:,is_main_ion),uif)

      !! poloidal and radial velocities of main ion on faces
      vif(:,0) = fnan(:,0)/nif
      vif(:,1) = fnan(:,1)/nif
      !! correct from poloidal to diamagnetic
      pit = dabs(geo%fcBb(:,0))/geo%fcBb(:,3)  !! pitch, b_theta / B
      pit2 = dabs(geo%fcBb(:,2))/geo%fcBb(:,3) !! b_phi / B
      vif(:,0) = (vif(:,0)-uif*pit)/pit2

      if (switch%b2tfn_mol_conv.eq.1) then
        conv = geo%fcPbs*umf
      else
        conv = 0.0_R8
      endif


*   ..compute convective coefficients (poloidal and radial)
      if (switch%b2tfn_mol_ionv.ne.0) then
        flo(:,0) = conv + 
     .             vif(:,0)*coeff_vif*geo%fcS*geo%fcQalf(:,0) 
        flo(:,1) = vif(:,1)*coeff_vif*geo%fcS*geo%fcQalf(:,1) 
      else
        flo(:,0) = conv
        flo(:,1) = 0.0_R8
      endif

*   ..apply discretization scheme for density-part of fna_mdf
      call calcflow(nCv, nFc, nVx, meth, geo, mpg, n_mol(:,0),
     &     flo, zeroes_fc0, fn_mol0,
     &     dummy_fc0, dummy_fc1)

*   ..compute Rhie-Chow term
      do iFc = 1, nFc
        cdpm_RC(iFc,0) = switch%b2npco_pcm1*0.5_R8*geo%fcPbs(iFc)*
     &   (dv%pcc_mol(mpg%fcCv(iFc,1),1)*n_mol(mpg%fcCv(iFc,1),0) +
     &    dv%pcc_mol(mpg%fcCv(iFc,2),1)*n_mol(mpg%fcCv(iFc,2),0))/
     &   (geo%fcQgam(iFc,0)*(geo%fcHc(iFc,1) + geo%fcHc(iFc,2)))
        dpbm(iFc) = 0.5_R8*geo%fcQgam(iFc,0)*
     &   (dpmc(mpg%fcCv(iFc,1)) + dpmc(mpg%fcCv(iFc,2)))*
     &   (geo%fcHc(iFc,1) + geo%fcHc(iFc,2))
      enddo
      cdpm_RC(:,1) = 0.0_R8
      term   = -cdpm_RC(:,0)*(dpm(:,0) - dpbm)

*   ..pressure-diffusion contribution
      fn_dp(:,0) = - co%cdp_mol(:,0)*dpm(:,0)
      fn_dp(:,1) = - co%cdp_mol(:,1)*dpm(:,1)


*   ..compute other fluxes
      dv%fn_mol(:,0) = fn_mol0(:,0) + fn_dp(:,0) + term 
      dv%fn_mol(:,1) = fn_mol0(:,1) + fn_dp(:,1)

      dv%fn_mol_fcor(:,0) = fn_mol0(:,0)*geo%fcHz - 
     .                      co%cdphz_mol(:,0)*dpm(:,0)
      dv%fn_mol_fcor(:,1) = fn_mol0(:,1)*geo%fcHz - 
     .                      co%cdphz_mol(:,1)*dpm(:,1)

*   ..compute flo0, con0
      call calccoef(nCv,nFc,nVx,meth,geo,
     &    flo,zeroes_fc0,flo0,con0)                  !! no density-driven part for fluid molecules -> zeroes_fc0
      do iFc = 1, nFc
        flo_tmp(iFc,0)   = flo_tmp(iFc,0) +
     &     0.5_R8*flo0(iFc,0)*
     &     (kmpr(mpg%fcCv(iFc,1))+kmpr(mpg%fcCv(iFc,2))) -
     &     con0(iFc,0)*(kmpr(mpg%fcCv(iFc,2))-kmpr(mpg%fcCv(iFc,1)))
        flo_tmp(iFc,1)   = flo_tmp(iFc,1) +
     &     0.5_R8*flo0(iFc,1)*
     &     (kmprv(mpg%fcVx(iFc,1))+kmprv(mpg%fcVx(iFc,2))) -
     &     con0(iFc,1)*(kmprv(mpg%fcVx(iFc,2))-kmprv(mpg%fcVx(iFc,1)))
        con_tmp(iFc,0,0) = con_tmp(iFc,0,0) +
     &     0.5_R8*con0(iFc,0)*
     &     (kmpr(mpg%fcCv(iFc,1))+kmpr(mpg%fcCv(iFc,2))) - 0.25_R8*
     &     flo0(iFc,0)*(kmpr(mpg%fcCv(iFc,2))-kmpr(mpg%fcCv(iFc,1)))
        con_tmp(iFc,1,0) = con_tmp(iFc,1,0) +
     &     0.5_R8*con0(iFc,1)*
     &     (kmprv(mpg%fcVx(iFc,1))+kmprv(mpg%fcVx(iFc,2))) - 0.25_R8*
     &     flo0(iFc,1)*(kmprv(mpg%fcVx(iFc,2))-kmprv(mpg%fcVx(iFc,1)))
      enddo


*   ..contribution from pressure diffusion and Rhie-Chow term
      cdpm = cdpm_RC + co%cdp_mol(:,:) 
      call calccoef(nCv,nFc,nVx,0,geo,zeroes_fc0,cdpm,flo1,con1)

      con_tmp(:,:,0) = con_tmp(:,:,0) + con1
      con_tmp(:,0,1) = 0.0_R8                       
      con_tmp(:,1,1) = 0.0_R8                       
      con_tmp(:,:,2) = con_tmp(:,:,1)

      flo_tmp = flo_tmp + flo1

      dv%flo_mol(:,:)    = flo_tmp
      dv%con_mol(:,:,:)  = con_tmp

c compute vv, ww and uu
       vv = dv%fn_mol(:,1)/(geo%fcS*geo%fcQalf(:,1)+eps)/nmf
       uu = dv%fn_mol(:,0)/(geo%fcS*geo%fcQalf(:,0)+eps)/nmf
       call intcell(nFc, nCv, mpg, mpg%intcellP,
     &               uu,uuc)
       call intcell(nFc, nCv, mpg, mpg%intcellR,
     &               vv,vvc)
       pitc = abs(geo%cvBb(:,0))/geo%cvBb(:,3)
       pit2c = abs(geo%cvBb(:,2))/geo%cvBb(:,3)
       wwc = n_moli(:,0)/n_mol(:,0)*pit2c*pl%ua(:,is_main_ion)

c limit also these velocity components to avoid numerical issues
       uuc =sign(min(abs(uuc),um_max),uuc)
       vvc =sign(min(abs(vvc),um_max),vvc)
       wwc =sign(min(abs(wwc),um_max),wwc)

       if (switch%b2stfm_recompute_perp_vels.ne.0) then
         u_pol_mol(:,0) = uuc
         u_rad_mol(:,0) = vvc
         u_tor_mol(:,0) = wwc !! TODO?
         u_dia_mol(:,0) = (u_pol_mol(:,0)-u_mol(:,0)*pitc)/pit2c
       else 
cc molecule velocities are kept at the value read from the input mol_bg
       endif

      if(switch%b2tfmol_iout .eq. 1 ) then      
       call my_out_us(70,nFc,1,dv%fn_mol(:,0),'b2tfn_mol_th')
       call my_out_us(70,nFc,1,dv%fn_mol(:,1),'b2tfn_mol_r')
       call my_out_us(70,nFc,1,co%cdp_mol(:,0),'b2tfn_cdp_mol_th')
       call my_out_us(70,nFc,1,co%cdp_mol(:,1),'b2tfn_cdp_mol_r')
       call my_out_us(70,nFc,1,dpm(:,0),'b2tfn_dpm_mol_th')
       call my_out_us(70,nFc,1,dpm(:,1),'b2tfn_dpm_mol_r')
       call my_out_us(70,nCv,0,uuc,'b2tfn_mol_uuc') ! poloidal
       call my_out_us(70,nCv,0,vvc,'b2tfn_mol_vvc') ! radial
       call my_out_us(70,nCv,0,wwc,'b2tfn_mol_wwc') ! toroidal
       call my_out_us(70,nFc,1,fn_dp(:,0),'b2tfn_dp_mol_th')
       call my_out_us(70,nFc,1,fn_dp(:,1),'b2tfn_dp_mol_r')
       call my_out_us(70,nFc,1,term,'b2tfn_term_mol_th')
       call my_out_us(70,nFc,1,fn_mol0(:,0),'b2tfn_mol_conv_th')
       call my_out_us(70,nFc,1,fn_mol0(:,1),'b2tfn_mol_conv_r')
       call my_out_us(70,nCv,0,u_pol_mol(:,0),'b2tfn_mol_u_pol') ! poloidal
       call my_out_us(70,nCv,0,u_rad_mol(:,0),'b2tfn_mol_u_rad') ! radial
       call my_out_us(70,nCv,0,u_tor_mol(:,0),'b2tfn_mol_u_tor') ! toroidal
      endif

* ..return
      call subend ()
      return

*-----------------------------------------------------------------------
*.end b2tfn_mol

      end subroutine b2tfn_mol
