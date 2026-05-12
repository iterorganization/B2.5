*-----------------------------------------------------------------------
*.specification

      subroutine b2npmo_mol (nCv, nFc, nVx,
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
*     B2NPMO_MOL performs one iteration for the solution of
*     the B2.5 equations for parallel velocity for fluid molecules.
*
*-----------------------------------------------------------------------
*.declarations

*   ..local variables
      integer :: iCv, ireg
      integer, save :: ncall_b2npmo_mol = 0
      real (kind=R8) ::
     .  sm_gp(nCv,0:3), sm_cf(nCv,0:3), 
     .  pm(nCv),
     .  flub(nFc,0:1), rob(nCv), ctcfb(nCv),
     .  flubv(nFc,0:1), flcb(nFc,0:1), cvcb(nFc,0:1),                                              
     .  wrkf0(nFc,0:1), 
     .  wrkv(nVx), aa(mpg%nCmxNv),
     .  sm_mol(nCv,0:3), ro(nCv),
     .  cst, ro1(nCv), pm1(nCv),                                
     .  wrk0(nCv),wrk1(nCv), mass,
     .  gtm(nCv), zeroes_f(nFc,0:1), zeroes_c(nFc)
      real (kind=R8), save :: rxg_npmo_mol = 0.0_R8
      logical new_matrix, warning
*   ..procedures
      external xertst, ipgetr, sfill, b2scopy, b2sscal, b2saxpy
      external
     &  b2xvsg, b2nxfc,
     &  b2sigp, b2sicf, b2nxcm, b2urmo, b2usmo

*-----------------------------------------------------------------------
*.computation

* ..preliminaries
*   ..subprogram start-up calls
      call subini ('b2npmo_mol')
*   ..test nCv, nFc
      call xertst (0.lt.nCv.and.0.lt.nFc, 'faulty argument nCv, nFc')
*   ..test itcnt
      call xertst (0.le.itcnt, 'faulty argument itcnt')
*   ..test rxf
      call xertst (0.0_R8.le.rxf.and.rxf.le.1.0_R8,
     . 'faulty argument rxf')

*   ..set internal parameters on first call
      rxg_npmo_mol = switch%b2npmo_rxg

*   ..compute ro
      mass = mp*am_d2
      ro = nm_ext(:,0)*mass
 
*   ..compute gtm
      call gradc_p(nCv, nFc, nVx, 0, geo, mpg, tm_ext(:,0), wrkv, gtm)


*   ..initialise
      call sfill (nCv, 0.0_R8, dv%resmo_mol, 1)
      zeroes_f = 0.0_R8
      zeroes_c = 0.0_R8
      sm_gp = 0.0_R8
      sm_cf = 0.0_R8

* ..main computation
*   ..compute individual species corrections

*    ..compute pm
      pm = tm_ext(:,0)*nm_ext(:,0)

*    ..compute pressure gradient term
      call b2sigp (nCv, nFc, nVx, 0, switch, geo, mpg, zeroes_c,
     &   nm_ext(:,0), tm_ext(:,0), tm_ext(:,0), dv%ne, pl%te, pl%po, 
     &   sm_gp, wrk0, wrk1)
*    ..compute flub
       flub = am_d2*mp*dv%fn_mol_fcor

      call b2nxfv (nCv, nFc, nVx, 0, 0, switch, geo, mpg,
     &    zeroes_c, zeroes_f, co%cvshz_mol(:,:), flubv)

      call b2saxpy (nFc*2, 1.0_R8, flubv, 1, flub, 1)  
               
*    ..compute flcb, cvcb
      call b2nxfc (nCv, nFc, nVx, 0, switch, geo, mpg, ro,
     &  um_ext(:,0), flub, co%cvshz_mol(:,:), flcb, cvcb,
     &  wrkf0)

*    ..compute centrifugal force term                                  
      call b2sicf (nCv, nFc, nVx, 0, am_d2*mp, switch,  geo, mpg,   
     *   nm_ext(:,0), um_ext(:,0), sm_cf, ctcfb, wrk0)

*    ..compute total parallel momentum source
      sm_mol(:,0) = sr_mol%sm_mol(:,0)
      sm_mol(:,1) = sr_mol%sm_mol(:,1)
      sm_mol(:,2) = 0.0_R8
      sm_mol(:,3) = 0.0_R8
      call b2saxpy (nCv*4, 1.0_R8, sm_gp, 1, sm_mol, 1)
      call b2saxpy (nCv*4, 1.0_R8, sm_cf, 1, sm_mol, 1)   

*    ..compute residual
      call b2urmo (nCv, nFc, nVx, switch, geo, mpg, 0, um_ext(:,0),
     &    ro, sm_mol, flub, co%cvshz_mol(:,:), dv%resmo_mol(:),
     &    dv%fm_mol(:,:), zeroes_f, co%cvshz_mol_eff(:,:))      


      if (switch%b2npmo_iout.ne.0) then 
        call my_out_us(70,nCv,0,dv%resmo_mol,'b2npmo_resmo_mol')
      endif                                                          

*    ..compute correction
      pm1 = pm !csc to avoid aliasing in Tapenade
      ro1 = ro !csc to avoid aliasing in Tapenade
      call b2usmo (nCv, nFc, nVx, switch, geo, mpg, mpg%nnreg(0),
     &    b2npmo_solvedum, itcnt, rxg_npmo_mol, ro, ro1, pm,
     &    pm1, um_ext(:,0), sm_mol,
     &    flcb, cvcb, dv%resmo_mol, ctcfb, dv%coru_mol,
     &    dv%pcc_mol(:,:), aa, wrk1, 'b2usmo_mol')


cc print out sources before update
      if(switch%b2tfmol_iout.ne.0) then     
        wrk0 = sm_gp(:,0) + sm_gp(:,1)*um_ext(:,0) + 
     &         sm_gp(:,2)*rob + sm_gp(:,3)*rob*um_ext(:,0)  
        call my_out_us(70,nCv,0,wrk0,'b2npmo_mol_smgp')
        wrk0 = sm_cf(:,0) + sm_cf(:,1)*um_ext(:,0) + 
     &         sm_cf(:,2)*rob + sm_cf(:,3)*rob*um_ext(:,0)  
        call my_out_us(70,nCv,0,wrk0,'b2npmo_mol_smcf')
        wrk0 = sm_mol(:,0) + sm_mol(:,1)*um_ext(:,0) + 
     &         sm_mol(:,2)*rob + sm_mol(:,3)*rob*um_ext(:,0)  
        call my_out_us(70,nCv,0,wrk0,'b2npmo_mol_tot')
        call my_out_us(70,nFc,1,dv%fm_mol(:,0),
     .                   'b2npmo_mol_fm_th')
        call my_out_us(70,nFc,1,dv%fm_mol(:,1),
     .                   'b2npmo_mol_fm_r')
        call my_out_us(70,nCv,0,dv%coru_mol, 'b2npmo_coru_mol')
        call my_out_us(70,nCv,0,dv%resmo_mol,'b2npmo_resmo_mol')
      end if


*   ..apply correction
      if (rxf.ne.0.0_R8) then
       do iCv = 1, nCv
        um_ext(iCv,0) = um_ext(iCv,0) + rxf*dv%coru_mol(iCv)
       enddo
      endif

      if( switch%ion_vlct_restrict .ne. 0 ) then                       
        do iCv = 1, nCv
          cst=sqrt(tm_ext(iCv,0)/(mp*am_d2))
          if( abs(um_ext(iCv,0)).gt.
     .         switch%ion_vlct_restrict_M*cst ) then
              um_ext(iCv,0)= sign(switch%ion_vlct_restrict_M*cst,
     .                              um_ext(iCv,0))                     
          endif
        enddo
      endif

      if( switch%ion_vlct_restrict_core .ne. 0 ) then                       
        do iCv = 1, nCv
          if (mpg%cvOnClosedSurface(iCv)) then
            cst=sqrt(tm_ext(iCv,0)/(mp*am_d2))
            if( abs(um_ext(iCv,0)).gt.
     .         switch%ion_vlct_restrict_core_M*cst ) then
              um_ext(iCv,0)= sign(switch%ion_vlct_restrict_core_M*cst,
     .                              um_ext(iCv,0))                     
            endif
          endif
        enddo
      endif


*   ..compute physical time-derivative after correction has been applied
*   .. m*d(nu)/dt = m*[u*(dn/dt) + n*(du/dt)]
! to obtain the true physical derivative, we need to divide by ts_factor, which was previously multiplied to smodt
        dv%dmdt_mol = ( -sr_mol%sm_dt_mol(:,0)                        
     &    - sr_mol%sm_dt_mol(:,1)*um_ext(:,0))/ts_factor

*   ..correct pcca
cc      call b2nxcm (nCv, ns, dv%pccm, dv%pcca) TODO ?? Not needed for a single momentum equation (only needed for total momentum equation for all species)
cc TODO to be included later ...




* ..return
      ncall_b2npmo_mol = ncall_b2npmo_mol+1
      call subend ()
      return
*-----------------------------------------------------------------------
*.end b2npmo_mol

      end subroutine b2npmo_mol
