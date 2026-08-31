*-----------------------------------------------------------------------
*.specification

      subroutine b2npmo_mol_3D (nCv, nFc, nVx,
     &  switch, geo, mpg,
     &  itcnt, rxf, dtim, pl, dv, co, rt, sr_mol,
     &  psnl)
      use b2mod_types
      use b2mod_math
      use b2mod_constants
      use b2mod_numerics_namelist
      use b2mod_transport_fun                                         
      use b2mod_b2cmpa
      use b2mod_switches
      use b2us_geo
      use b2us_map
      use b2us_plasma
      use b2mod_openmp
      use b2mod_ma28_for_us
     , , only : restart_ma28_for_us
      use b2mod_fluid_molecules
      use b2mod_subsys
      use b2mod_ad
     & , only : b2npmo_solvedum
      implicit none
*   ..input arguments (unchanged on exit)
      integer nCv, nFc, nVx, itcnt                        
      type (switches), intent (in) :: switch
      type (geometry), intent (in) :: geo
      type (mapping), intent (in) :: mpg
      type (B2Rates), intent (in) :: rt
      type (B2Coeff), intent (in) :: co
      real (kind=R8) ::
     *  rxf, dtim
*   ..input/output arguments
      type (B2Plasma), intent (inout) :: pl
      type (B2Derivatives), intent (inout) :: dv
      type (B2SourceMol), intent (inout) :: sr_mol
      type (B2PlasmaSnapshot), intent (inout) :: psnl

*-----------------------------------------------------------------------
*.documentation
*
*  1. purpose
*
*     B2NPMO_MOL_3D performs one iteration for the solution of
*     the 3D velocities of fluid molecules.
*
*-----------------------------------------------------------------------
*.declarations

*   ..local variables
      integer :: iCv, ireg
      real (kind=R8) ::
     .  smgp_pol(nCv), smgp_rad(nCv),
     .  pm(nCv),
     .  flub(nFc,0:1), rob(nCv), 
     .  flubv(nFc,0:1), flub_temp(nFc,0:1),
     .  flcb_rad(nFc,0:1), cvcb_rad(nFc,0:1),      
     .  flcb_pol(nFc,0:1), cvcb_pol(nFc,0:1),  
     .  flcb_tor(nFc,0:1), cvcb_tor(nFc,0:1),                                                                                                                                                                          
     .  wrkf0(nFc,0:1), 
     .  wrkv(nVx), aa(mpg%nCmxNv),
     .  sm_mol_rad(nCv,0:3), sm_mol_pol(nCv,0:3), sm_mol_tor(nCv,0:3), 
     .  ro(nCv), cst, cvs_temp(nFc,0:1),                                
     .  wrk0(nCv),wrk1(nCv), mass, dummy_out(nCv,0:1),
     .  zeroes_f(nFc,0:1), zeroes_c(nFc), zeroes_nCv(nCv),
     .  pit(nCv), pit2(nCv), smcu_pol(nCv), smcu_rad(nCv), smcu_tor(nCv)
      real (kind=R8), save :: rxg_npmo_mol_perp = 0.0_R8
      logical new_matrix, warning
*   ..procedures
      external xertst, ipgetr, sfill, b2scopy, b2sscal, b2saxpy
      external
     &  b2xvsg, b2nxfc,
     &  b2sicf, b2nxcm, b2urmo, b2usmo_perp2

*-----------------------------------------------------------------------
*.computation

* ..preliminaries
*   ..subprogram start-up calls
      call subini ('b2npmo_mol_perp2')
*   ..test nCv, nFc
      call xertst (0.lt.nCv.and.0.lt.nFc, 'faulty argument nCv, nFc')
*   ..test itcnt
      call xertst (0.le.itcnt, 'faulty argument itcnt')
*   ..test rxf
      call xertst (0.0_R8.le.rxf.and.rxf.le.1.0_R8,
     . 'faulty argument rxf')

*   ..set internal parameters on first call
      rxg_npmo_mol_perp = switch%b2npmo_rxg

      pit = abs(geo%cvBb(:,0))/geo%cvBb(:,3)
      pit2 = abs(geo%cvBb(:,2))/geo%cvBb(:,3)

*   ..compute ro
      mass = mp*am_d2
      ro = n_mol(:,0)*mass
 
*   ..initialise
      call sfill (nCv, 0.0_R8, dv%resmo_mol_pol, 1)
      call sfill (nCv, 0.0_R8, dv%resmo_mol_rad, 1)
      call sfill (nCv, 0.0_R8, dv%resmo_mol_tor, 1)
      zeroes_f = 0.0_R8
      zeroes_c = 0.0_R8
      zeroes_nCv = 0.0_R8
* ..main computation

*    ..compute pm
      pm = t_mol(:,0)*n_mol(:,0)

*    ..compute pressure gradient term
      call b2smgp_3D (nCv, nFc, nVx, switch, geo, mpg,
     &   n_mol(:,0), t_mol(:,0), smgp_rad, smgp_pol)
     !! toroidal pressure gradient is zero

*    ..compute flub
       flub = am_d2*mp*dv%fn_mol

      !! FM_TODO ??
!      call b2nxfv (nCv, nFc, nVx, 0, 0, switch, geo, mpg,
!     &    zeroes_c, zeroes_f, co%cvshz_mol(:,:), flubv)
!      call b2saxpy (nFc*2, 1.0_R8, flubv, 1, flub, 1)  

      if (switch%b2tfmol_3D_style1.eq.1) then
        flub_temp = flub
      else
        flub_temp = 0.0_R8 
      endif

      if (switch%b2tfmol_3D_style2.eq.1) then
        cvs_temp = co%cvs_mol 
      else
        cvs_temp = 0.0_R8 
      endif
               


     !! FM_TODO, neglected for now
*    ..compute centrifugal force term                                  
!      call b2sicf (nCv, nFc, nVx, 0, am_d2*mp, switch,  geo, mpg,   
!     *   n_mol(:,0), u_mol(:,0), sm_cf, ctcfb, wrk0)

      call b2smcu_3D (nCv, nFc, nVx, am_d2*mp, switch,  geo, mpg,   
     *   n_mol(:,0), u_pol_mol(:,0), u_rad_mol(:,0), u_tor_mol(:,0),
     *   smcu_pol, smcu_rad, smcu_tor)

*    ..compute total perpendicular momentum sources
      sm_mol_rad(:,0) = sr_mol%sm_mol_rad(:,0) + smgp_rad + smcu_rad
      sm_mol_rad(:,1) = sr_mol%sm_mol_rad(:,1)
      sm_mol_rad(:,2) = 0.0_R8
      sm_mol_rad(:,3) = 0.0_R8 

      sm_mol_pol(:,0) = sr_mol%sm_mol_pol(:,0) + smgp_pol + smcu_pol
      sm_mol_pol(:,1) = sr_mol%sm_mol_pol(:,1)
      sm_mol_pol(:,2) = 0.0_R8
      sm_mol_pol(:,3) = 0.0_R8 

      sm_mol_tor(:,0) = sr_mol%sm_mol_tor(:,0) + smcu_tor
      sm_mol_tor(:,1) = sr_mol%sm_mol_tor(:,1)
      sm_mol_tor(:,2) = 0.0_R8
      sm_mol_tor(:,3) = 0.0_R8 

*    ..compute residual
      call b2urmo (nCv, nFc, nVx, switch, geo, mpg, 0, u_rad_mol(:,0),
     &    ro, sm_mol_rad, flub_temp, cvs_temp, dv%resmo_mol_rad(:),
     &    dv%fm_mol_rad(:,:), zeroes_f, wrkf0)  

      call b2urmo (nCv, nFc, nVx, switch, geo, mpg, 0, u_pol_mol(:,0),
     &    ro, sm_mol_pol, flub_temp, cvs_temp, dv%resmo_mol_pol(:),
     &    dv%fm_mol_pol(:,:), zeroes_f, wrkf0)  

      call b2urmo (nCv, nFc, nVx, switch, geo, mpg, 0, u_tor_mol(:,0),
     &    ro, sm_mol_tor, flub_temp, cvs_temp, dv%resmo_mol_tor(:),
     &    dv%fm_mol_tor(:,:), zeroes_f, wrkf0)  

      if (switch%b2npmo_iout.ne.0) then 
        call my_out_us(70,nCv,0,dv%resmo_mol_rad,'b2npmo_resmo_mol_rad')
        call my_out_us(70,nCv,0,dv%resmo_mol_dia,'b2npmo_resmo_mol_pol')
        call my_out_us(70,nCv,0,dv%resmo_mol_dia,'b2npmo_resmo_mol_tor')
      endif  
      
      
*    ..compute flcb, cvcb
      call b2nxfc_3D (nCv, nFc, nVx, 0, switch, geo, mpg, ro,
     &  u_rad_mol(:,0), flub_temp, cvs_temp, flcb_rad, cvcb_rad,
     &  wrkf0)
      call b2nxfc_3D (nCv, nFc, nVx, 0, switch, geo, mpg, ro,
     &  u_pol_mol(:,0), flub_temp, cvs_temp, flcb_pol, cvcb_pol,
     &  wrkf0)
      call b2nxfc_3D (nCv, nFc, nVx, 0, switch, geo, mpg, ro,
     &  u_tor_mol(:,0), flub_temp, cvs_temp, flcb_tor, cvcb_tor,
     &  wrkf0)

*    ..compute correction
      call b2usmo_3D (nCv, nFc, nVx, switch, geo, mpg, mpg%nnreg(0),
     &    b2npmo_solvedum, itcnt, rxg_npmo_mol_perp, ro, pm,
     &    u_rad_mol(:,0), sm_mol_rad,
     &    flcb_rad, cvcb_rad, dv%resmo_mol_rad, zeroes_nCv, 
     &    dv%coru_mol_rad, dv%pcc_mol_rad, aa, wrk1, 'b2usmo_mol_rad') 

      call b2usmo_3D (nCv, nFc, nVx, switch, geo, mpg, mpg%nnreg(0),
     &    b2npmo_solvedum, itcnt, rxg_npmo_mol_perp, ro, pm,
     &    u_pol_mol(:,0), sm_mol_pol,
     &    flcb_pol, cvcb_pol, dv%resmo_mol_pol, zeroes_nCv, 
     &    dv%coru_mol_pol, dv%pcc_mol_pol, aa, wrk1, 'b2usmo_mol_pol') 

      call b2usmo_3D (nCv, nFc, nVx, switch, geo, mpg, mpg%nnreg(0),
     &    b2npmo_solvedum, itcnt, rxg_npmo_mol_perp, ro, pm,
     &    u_tor_mol(:,0), sm_mol_tor,
     &    flcb_tor, cvcb_tor, dv%resmo_mol_tor, zeroes_nCv, 
     &    dv%coru_mol_tor, zeroes_f, aa, wrk1, 'b2usmo_mol_tor') 


cc print out sources before update
      if(switch%b2tfmol_iout.ne.0) then     
        call my_out_us(70,nCv,0,smgp_rad,'b2npmo_mol_smgp_rad')
        call my_out_us(70,nCv,0,smgp_pol,'b2npmo_mol_smgp_pol')

        wrk0 = sm_mol_rad(:,0) + sm_mol_rad(:,1)*u_rad_mol(:,0) + 
     &         sm_mol_rad(:,2)*rob + sm_mol_rad(:,3)*rob*u_rad_mol(:,0)  
        call my_out_us(70,nCv,0,wrk0,'b2npmo_mol_sm_rad_tot')
        wrk0 = sm_mol_pol(:,0) + sm_mol_pol(:,1)*u_pol_mol(:,0) + 
     &         sm_mol_pol(:,2)*rob + sm_mol_pol(:,3)*rob*u_pol_mol(:,0)  
        call my_out_us(70,nCv,0,wrk0,'b2npmo_mol_sm_pol_tot')
        wrk0 = sm_mol_tor(:,0) + sm_mol_tor(:,1)*u_tor_mol(:,0) + 
     &         sm_mol_tor(:,2)*rob + sm_mol_tor(:,3)*rob*u_tor_mol(:,0)  
        call my_out_us(70,nCv,0,wrk0,'b2npmo_mol_sm_tor_tot')

        call my_out_us(70,nFc,1,dv%fm_mol_rad(:,0),
     .                   'b2npmo_mol_fm_rad_th')
        call my_out_us(70,nFc,1,dv%fm_mol_rad(:,1),
     .                   'b2npmo_mol_fm_rad_r')
        call my_out_us(70,nFc,1,dv%fm_mol_pol(:,0),
     .                   'b2npmo_mol_fm_pol_th')
        call my_out_us(70,nFc,1,dv%fm_mol_pol(:,1),
     .                   'b2npmo_mol_fm_pol_r')
        call my_out_us(70,nFc,1,dv%fm_mol_tor(:,0),
     .                   'b2npmo_mol_fm_tor_th')
        call my_out_us(70,nFc,1,dv%fm_mol_tor(:,1),
     .                   'b2npmo_mol_fm_tor_r') 
        call my_out_us(70,nCv,0,dv%coru_mol_rad, 'b2npmo_coru_mol_rad')
        call my_out_us(70,nCv,0,dv%coru_mol_pol, 'b2npmo_coru_mol_pol')
        call my_out_us(70,nCv,0,dv%coru_mol_tor, 'b2npmo_coru_mol_tor')
        call my_out_us(70,nCv,0,dv%resmo_mol_rad,'b2npmo_resmo_mol_rad')
        call my_out_us(70,nCv,0,dv%resmo_mol_pol,'b2npmo_resmo_mol_pol')
        call my_out_us(70,nCv,0,dv%resmo_mol_tor,'b2npmo_resmo_mol_tor')
        !! sources in parts
       call my_out_us(70,nCv,0,sm_mol_rad(:,0),'b2npmo_mol_sm_rad0_tot')
       call my_out_us(70,nCv,0,sm_mol_rad(:,1),'b2npmo_mol_sm_rad1_tot')
       call my_out_us(70,nCv,0,sm_mol_rad(:,2),'b2npmo_mol_sm_rad2_tot')
       call my_out_us(70,nCv,0,sm_mol_rad(:,3),'b2npmo_mol_sm_rad3_tot')
       call my_out_us(70,nCv,0,sm_mol_pol(:,0),'b2npmo_mol_sm_pol0_tot')
       call my_out_us(70,nCv,0,sm_mol_pol(:,1),'b2npmo_mol_sm_pol1_tot')
       call my_out_us(70,nCv,0,sm_mol_pol(:,2),'b2npmo_mol_sm_pol2_tot')
       call my_out_us(70,nCv,0,sm_mol_pol(:,3),'b2npmo_mol_sm_pol3_tot')
       call my_out_us(70,nCv,0,sm_mol_tor(:,0),'b2npmo_mol_sm_tor0_tot')
       call my_out_us(70,nCv,0,sm_mol_tor(:,1),'b2npmo_mol_sm_tor1_tot')
       call my_out_us(70,nCv,0,sm_mol_tor(:,2),'b2npmo_mol_sm_tor2_tot')
       call my_out_us(70,nCv,0,sm_mol_tor(:,3),'b2npmo_mol_sm_tor3_tot')
      end if


*   ..apply correction
      if (rxf.ne.0.0_R8) then
       do iCv = 1, nCv
        if (switch%b2tfmol_3D_rad.eq.1) then
          u_rad_mol(iCv,0) = u_rad_mol(iCv,0) + rxf*dv%coru_mol_rad(iCv)
        else
          dv%resmo_mol_rad = 0.0_R8
        endif
        if (switch%b2tfmol_3D_pol.eq.1) then
          u_pol_mol(iCv,0) = u_pol_mol(iCv,0) + rxf*dv%coru_mol_pol(iCv)
        else
          dv%resmo_mol_pol = 0.0_R8
        endif
        if (switch%b2tfmol_3D_tor.eq.1) then
          u_tor_mol(iCv,0) = u_tor_mol(iCv,0) + rxf*dv%coru_mol_tor(iCv)
        else 
          dv%resmo_mol_tor = 0.0_R8
        endif
       enddo
      endif

*   ..recompute the other velocity components
      u_mol(:,0) = pit*u_pol_mol(:,0) + pit2*u_tor_mol(:,0)
      u_dia_mol(:,0) = pit2*u_pol_mol(:,0) - pit*u_tor_mol(:,0)

c      if( switch%ion_vlct_restrict .ne. 0 ) then                       
c        do iCv = 1, nCv
c          cst=sqrt(t_mol(iCv,0)/(mp*am_d2))
c          if( abs(u_rad_mol(iCv,0)).gt.
c     .         switch%ion_vlct_restrict_M*cst ) then
c              u_rad_mol(iCv,0)= sign(switch%ion_vlct_restrict_M*cst,
c     .                              u_rad_mol(iCv,0))                     
c          endif
c          if( abs(u_dia_mol(iCv,0)).gt.
c     .         switch%ion_vlct_restrict_M*cst ) then
c              u_dia_mol(iCv,0)= sign(switch%ion_vlct_restrict_M*cst,
c     .                              u_dia_mol(iCv,0))                     
c          endif
c        enddo
c      endif



*   ..compute physical time-derivative after correction has been applied
*   .. m*d(nu)/dt = m*[u*(dn/dt) + n*(du/dt)]
! to obtain the true physical derivative, we need to divide by ts_factor, which was previously multiplied to smodt
!        dv%dmdt_mol = ( -sr_mol%sm_dt_mol(:,0)                        
!     &    - sr_mol%sm_dt_mol(:,1)*u_mol(:,0))/ts_factor
! FM_TODO add later if needed

*   ..correct pcca
cc      call b2nxcm (nCv, ns, dv%pccm, dv%pcca) TODO ?? Not needed for a single momentum equation (only needed for total momentum equation for all species)
cc FM_TODO to be included later ...

* ..return
      call subend ()
      return
*-----------------------------------------------------------------------
*.end b2npmo_mol

      end subroutine b2npmo_mol_3D
