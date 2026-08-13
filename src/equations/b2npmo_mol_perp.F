*-----------------------------------------------------------------------
*.specification

      subroutine b2npmo_mol_perp (nCv, nFc, nVx,
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
*     B2NPMO_MOL_PERP performs one iteration for the solution of
*     the perpendicular velocities of fluid molecules.
*
*-----------------------------------------------------------------------
*.declarations

*   ..local variables
      integer :: iCv, ireg
      real (kind=R8) ::
     .  smgp_rad(nCv), smgp_dia(nCv),
     .  pm(nCv),
     .  flub(nFc,0:1), rob(nCv), 
     .  flubv(nFc,0:1), 
     .  flcb_rad(nFc,0:1), cvcb_rad(nFc,0:1),      
     .  flcb_dia(nFc,0:1), cvcb_dia(nFc,0:1),                                                                                      
     .  wrkf0(nFc,0:1), 
     .  wrkv(nVx), aa(mpg%nCmxNv),
     .  sm_mol_rad(nCv,0:3), sm_mol_dia(nCv,0:3), ro(nCv),
     .  cst, cvs_temp(nFc,0:1),                                
     .  wrk0(nCv),wrk1(nCv), mass, dummy_out(nCv,0:1),
     .  zeroes_f(nFc,0:1), zeroes_c(nFc), zeroes_nCv(nCv)
      real (kind=R8), save :: rxg_npmo_mol_perp = 0.0_R8
      logical new_matrix, warning
*   ..procedures
      external xertst, ipgetr, sfill, b2scopy, b2sscal, b2saxpy
      external
     &  b2xvsg, b2nxfc,
     &  b2sigp, b2sicf, b2nxcm, b2urmo, b2usmo_perp

*-----------------------------------------------------------------------
*.computation

* ..preliminaries
*   ..subprogram start-up calls
      call subini ('b2npmo_mol_perp')
*   ..test nCv, nFc
      call xertst (0.lt.nCv.and.0.lt.nFc, 'faulty argument nCv, nFc')
*   ..test itcnt
      call xertst (0.le.itcnt, 'faulty argument itcnt')
*   ..test rxf
      call xertst (0.0_R8.le.rxf.and.rxf.le.1.0_R8,
     . 'faulty argument rxf')

*   ..set internal parameters on first call
      rxg_npmo_mol_perp = switch%b2npmo_rxg

*   ..compute ro
      mass = mp*am_d2
      ro = n_mol(:,0)*mass
 
*   ..initialise
      call sfill (nCv, 0.0_R8, dv%resmo_mol_dia, 1)
      call sfill (nCv, 0.0_R8, dv%resmo_mol_rad, 1)
      zeroes_f = 0.0_R8
      zeroes_c = 0.0_R8
      zeroes_nCv = 0.0_R8
* ..main computation

*    ..compute pm
      pm = t_mol(:,0)*n_mol(:,0)

*    ..compute pressure gradient term
      call b2smgp_perp (nCv, nFc, nVx, switch, geo, mpg,
     &   n_mol(:,0), t_mol(:,0), smgp_rad, smgp_dia)

*    ..compute flub
       flub = am_d2*mp*dv%fn_mol

      !! FM_TODO ??
!      call b2nxfv (nCv, nFc, nVx, 0, 0, switch, geo, mpg,
!     &    zeroes_c, zeroes_f, co%cvshz_mol(:,:), flubv)
!      call b2saxpy (nFc*2, 1.0_R8, flubv, 1, flub, 1)  

      if (switch%b2tfmol_3D_style2.eq.1) then
        cvs_temp = co%cvs_mol 
      else
        cvs_temp = 0.0_R8 
      endif
               
*    ..compute flcb, cvcb
      call b2nxfc (nCv, nFc, nVx, 0, switch, geo, mpg, ro,
     &  u_rad_mol(:,0), flub, cvs_temp, flcb_rad, cvcb_rad,
     &  wrkf0)
      call b2nxfc (nCv, nFc, nVx, 0, switch, geo, mpg, ro,
     &  u_dia_mol(:,0), flub, cvs_temp, flcb_dia, cvcb_dia,
     &  wrkf0)

      if (switch%b2tfmol_3D_style1.eq.0) then
        flcb_dia = 0.0_R8
        flcb_rad = 0.0_R8
        cvcb_dia = 0.0_R8
        cvcb_rad = 0.0_R8
      endif


     !! FM_TODO, neglected for now
*    ..compute centrifugal force term                                  
!      call b2sicf (nCv, nFc, nVx, 0, am_d2*mp, switch,  geo, mpg,   
!     *   n_mol(:,0), u_mol(:,0), sm_cf, ctcfb, wrk0)

*    ..compute total perpendicular momentum sources
      sm_mol_rad(:,0) = sr_mol%sm_mol_rad(:,0) + smgp_rad
      sm_mol_rad(:,1) = sr_mol%sm_mol_rad(:,1)
      sm_mol_rad(:,2) = 0.0_R8
      sm_mol_rad(:,3) = 0.0_R8 

      sm_mol_dia(:,0) = sr_mol%sm_mol_dia(:,0) + smgp_dia
      sm_mol_dia(:,1) = sr_mol%sm_mol_dia(:,1)
      sm_mol_dia(:,2) = 0.0_R8
      sm_mol_dia(:,3) = 0.0_R8 

*    ..compute residual
      call b2urmo (nCv, nFc, nVx, switch, geo, mpg, 0, u_rad_mol(:,0),
     &    ro, sm_mol_rad, flub, cvs_temp, dv%resmo_mol_rad(:),
     &    dv%fm_mol_rad(:,:), zeroes_f, wrkf0)  

      call b2urmo (nCv, nFc, nVx, switch, geo, mpg, 0, u_dia_mol(:,0),
     &    ro, sm_mol_dia, flub, cvs_temp, dv%resmo_mol_dia(:),
     &    dv%fm_mol_dia(:,:), zeroes_f, wrkf0)  

      if (switch%b2npmo_iout.ne.0) then 
        call my_out_us(70,nCv,0,dv%resmo_mol_rad,'b2npmo_resmo_mol_rad')
        call my_out_us(70,nCv,0,dv%resmo_mol_dia,'b2npmo_resmo_mol_dia')
      endif                                                          

*    ..compute correction
      call b2usmo_perp (nCv, nFc, nVx, switch, geo, mpg, mpg%nnreg(0),
     &    b2npmo_solvedum, itcnt, rxg_npmo_mol_perp, ro, pm,
     &    u_rad_mol(:,0), sm_mol_rad,
     &    flcb_rad, cvcb_rad, dv%resmo_mol_rad, zeroes_nCv, 
     &    dv%coru_mol_rad, aa, wrk1, 'b2usmo_mol_rad') 

      call b2usmo_perp (nCv, nFc, nVx, switch, geo, mpg, mpg%nnreg(0),
     &    b2npmo_solvedum, itcnt, rxg_npmo_mol_perp, ro, pm,
     &    u_dia_mol(:,0), sm_mol_dia,
     &    flcb_dia, cvcb_dia, dv%resmo_mol_dia, zeroes_nCv, 
     &    dv%coru_mol_dia, aa, wrk1, 'b2usmo_mol_dia') 


cc print out sources before update
      if(switch%b2tfmol_iout.ne.0) then     
        call my_out_us(70,nCv,0,smgp_rad,'b2npmo_mol_smgp_rad')
        call my_out_us(70,nCv,0,smgp_dia,'b2npmo_mol_smgp_dia')
        wrk0 = sm_mol_rad(:,0) + sm_mol_rad(:,1)*u_rad_mol(:,0) + 
     &         sm_mol_rad(:,2)*rob + sm_mol_rad(:,3)*rob*u_rad_mol(:,0)  
        call my_out_us(70,nCv,0,wrk0,'b2npmo_mol_sm_rad_tot')
        wrk0 = sm_mol_dia(:,0) + sm_mol_dia(:,1)*u_dia_mol(:,0) + 
     &         sm_mol_dia(:,2)*rob + sm_mol_dia(:,3)*rob*u_dia_mol(:,0)  
        call my_out_us(70,nCv,0,wrk0,'b2npmo_mol_sm_dia_tot')
        call my_out_us(70,nFc,1,dv%fm_mol_rad(:,0),
     .                   'b2npmo_mol_fm_rad_th')
        call my_out_us(70,nFc,1,dv%fm_mol_rad(:,1),
     .                   'b2npmo_mol_fm_rad_r')
        call my_out_us(70,nFc,1,dv%fm_mol_dia(:,0),
     .                   'b2npmo_mol_fm_dia_th')
        call my_out_us(70,nFc,1,dv%fm_mol_dia(:,1),
     .                   'b2npmo_mol_fm_dia_r')
        call my_out_us(70,nCv,0,dv%coru_mol_rad, 'b2npmo_coru_mol_rad')
        call my_out_us(70,nCv,0,dv%coru_mol_dia, 'b2npmo_coru_mol_dia')
        call my_out_us(70,nCv,0,dv%resmo_mol_rad,'b2npmo_resmo_mol_rad')
        call my_out_us(70,nCv,0,dv%resmo_mol_dia,'b2npmo_resmo_mol_dia')
        !! sources in parts
       call my_out_us(70,nCv,0,sm_mol_dia(:,0),'b2npmo_mol_sm_dia0_tot')
       call my_out_us(70,nCv,0,sm_mol_dia(:,1),'b2npmo_mol_sm_dia1_tot')
       call my_out_us(70,nCv,0,sm_mol_dia(:,2),'b2npmo_mol_sm_dia2_tot')
       call my_out_us(70,nCv,0,sm_mol_dia(:,3),'b2npmo_mol_sm_dia3_tot')
       call my_out_us(70,nCv,0,sm_mol_rad(:,0),'b2npmo_mol_sm_rad0_tot')
       call my_out_us(70,nCv,0,sm_mol_rad(:,1),'b2npmo_mol_sm_rad1_tot')
       call my_out_us(70,nCv,0,sm_mol_rad(:,2),'b2npmo_mol_sm_rad2_tot')
       call my_out_us(70,nCv,0,sm_mol_rad(:,3),'b2npmo_mol_sm_rad3_tot')

      end if


*   ..apply correction
      if (rxf.ne.0.0_R8) then
       do iCv = 1, nCv
        u_rad_mol(iCv,0) = u_rad_mol(iCv,0) + rxf*dv%coru_mol_rad(iCv)
        u_dia_mol(iCv,0) = u_dia_mol(iCv,0) + rxf*dv%coru_mol_dia(iCv)
       enddo
      endif

      if( switch%ion_vlct_restrict .ne. 0 ) then                       
        do iCv = 1, nCv
          cst=sqrt(t_mol(iCv,0)/(mp*am_d2))
          if( abs(u_rad_mol(iCv,0)).gt.
     .         switch%ion_vlct_restrict_M*cst ) then
              u_rad_mol(iCv,0)= sign(switch%ion_vlct_restrict_M*cst,
     .                              u_rad_mol(iCv,0))                     
          endif
          if( abs(u_dia_mol(iCv,0)).gt.
     .         switch%ion_vlct_restrict_M*cst ) then
              u_dia_mol(iCv,0)= sign(switch%ion_vlct_restrict_M*cst,
     .                              u_dia_mol(iCv,0))                     
          endif
        enddo
      endif



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

      end subroutine b2npmo_mol_perp
