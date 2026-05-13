*-----------------------------------------------------------------------
*.specification

      subroutine b2npht_mol (nCv, nFc, nVx, switch, geo, mpg,           
     &  itcnt,  ua_solve, rxf, dtim,
     &  pl, dv, co, rt, sr_mol, psnl)
      use b2mod_types
      use b2mod_constants
      use b2mod_b2cmpa
      use b2mod_switches
      use b2us_geo
      use b2us_map
      use b2us_plasma
      use b2mod_ma28_for_us
     , , only : restart_ma28_for_us
      use b2mod_numerics_namelist                                       
     & , only : time_factor, dtei, dtee, dten, ts_factor              
      use b2mod_fluid_molecules
      use b2mod_subsys
      implicit none
*   ..input arguments (unchanged on exit)
      integer nCv, nFc, nVx, itcnt 
      integer, save :: ncall_b2npht_mol = 0                       
      type (switches) :: switch
      type (geometry), intent (in) :: geo
      type (mapping), intent (in) :: mpg
      type (B2PlasmaSnapshot), intent (inout) :: psnl
      real (kind=R8) :: rxf, dtim
      logical :: ua_solve
*   ..input/output arguments
      type (B2Plasma), intent (inout) :: pl
      type (B2Derivatives), intent (inout) :: dv
      type (B2Coeff), intent (inout) :: co
      type (B2SourceMol), intent (inout) :: sr_mol
      type (B2Rates), intent (in) :: rt
*-----------------------------------------------------------------------
*.documentation
*
*  1. purpose
*
*     B2NPHT_MOL performs one iteration of a relaxation process for the
*     solution of the equations for internal energy of fluid molecules.
*
*-----------------------------------------------------------------------
*.declarations

*   ..local variables
      integer iCv
      real (kind=R8) ::
     * aa(mpg%nCmxNv), wrkf(nFc,0:1), wrkc(nFc,0:1,0:2), tm0(nCv),
     * nm0(nCv), um0(nCv), sh_final(nCv,0:3), sh0_mol(nCv,0:3), wrk(nCv)                                         
*   ..procedures
      intrinsic sqrt
      external xertst, ipgetr, b2saxpy
      external b2xvsg, b2ursd, b2usht, b2upht                           
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
*-----------------------------------------------------------------------
*.computation

* ..preliminaries
*   ..subprogram start-up calls
      call subini ('b2npht_mol')
*   ..test nCv, nFc
      call xertst (0.lt.nCv.and.0.lt.nFc, 'faulty argument nCv, nFc')
*   ..test rxf
      call xertst (0.0_R8.le.rxf.and.rxf.le.1.0_R8,
     . 'faulty argument rxf')

*   ..extensive tests on first few calls
      if (ncall_b2npht_mol.lt.3) then

*    ..test sign of sh_mol
       call b2xvsg (nCv, sr_mol%sh_mol(1,0), 1, 'sh_mol0', '.ge.')
       call b2xvsg (nCv, sr_mol%sh_mol(1,1), 1, 'sh_mol1', '.le.')

*    ..test sign of tm
       call b2xvsg (nCv, t_mol(:,0), 1, 'tm', '.gt.')

      endif

* ..initizalization
      tm0 = t_mol(:,0)
      nm0 = n_mol(:,0)
      um0 = u_mol(:,0)
      sh_final(:,0) = sr_mol%sh_mol(:,0)
      sh_final(:,1) = sr_mol%sh_mol(:,1)
      sh_final(:,2) = 0.0_R8
      sh_final(:,3) = 0.0_R8

* ..main computation

*   ..compute final heat sources
      call b2sihs_mol (nCv, nFc, nVx, switch, geo, mpg,
     &  co%cvs_mol,  co%cvshz_mol_eff, sh0_mol)

      call b2saxpy (nCv*4, 1.0_R8, sh0_mol, 1, sh_final, 1)

*   ..compute resh_mol
      call b2ursd (nCv, nFc, mpg, nm0, tm0, sh_final, dv%fh_mol,
     &   dv%resh_mol)

cc        switch%b2uxus_style = 2        

*   ..compute correction
      call restart_ma28_for_us
      call b2usht_mol (nCv, nFc, nVx, geo, mpg,
     &  itcnt, nm0, tm0, dv%floh_mol, dv%conh_mol,
     &  sh_final, dv%resh_mol,
     &  dv%cort_mol, wrkf, wrkc, aa, switch)

*   ..apply correction
      call b2upht_mol (nCv, nFc, nVx, switch, geo, mpg,
     &  ua_solve, rxf, dv%cort_mol,
     &  dv%pcc_mol(1:nCv,0),nm0, um0, tm0)

*   ..compute physical time-derivative after correction has been applied
*   .. m*d(nT)/dt = m*[n*(dT/dt) + T*(dn/dt)]
      dv%dhdt_mol = - ( sr_mol%sh_dt_mol(:,0)
     &       +sr_mol%sh_dt_mol(:,1)*tm0)/ts_factor

      t_mol(:,0) = tm0

      if(switch%b2tfmol_iout.ne.0) then       
        call my_out_us(70,nCv,0,tm0,'b2npht_mol_tm')
        call my_out_us(70,nCv,0,dv%resh_mol,'b2npht_resh_mol')
        wrk = sh0_mol(:,0) + sh0_mol(:,1)*tm0 + sh0_mol(:,2)*nm0 +
     &        sh0_mol(:,3)*nm0*tm0
        call my_out_us(70,nCv,0,wrk,'b2sihs_sh0_mol')
      end if

* ..return
      ncall_b2npht_mol = ncall_b2npht_mol+1
      call subend ()
      return

*-----------------------------------------------------------------------
*.end b2npht_mol

      end subroutine b2npht_mol
