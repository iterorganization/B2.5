*-----------------------------------------------------------------------
*.specification

      subroutine b2npco_mol (nCv, nFc, nVx, nregionv, solve_ua, 
     &  solvcreg, rxf, dtim, switch, geo, mpg, pl, dv, sr_mol, psnl)
      use b2mod_types
      use b2mod_b2cmpa
      use b2mod_switches
      use b2us_geo
      use b2us_map
      use b2mod_openmp
      use b2us_plasma
      use b2mod_fluid_molecules
      use b2mod_numerics_namelist                            
     & , only : dtco, time_factor, ts_factor                             
      use b2mod_subsys
      implicit none
*   ..input arguments (unchanged on exit)
      integer, intent(in) :: nCv, nFc, nVx, nregionv        
      type (switches), intent (in) :: switch
      type (geometry), intent (in) :: geo
      type (mapping), intent (in) :: mpg
      type (B2SourceMol), intent (inout) :: sr_mol
      logical solvcreg(0:nregionv), solve_ua
      real (kind=R8) ::
     *  rxf, dtim
*   ..input/output arguments
      type (B2Plasma), intent (inout) :: pl
      type (B2Derivatives), intent (inout) :: dv
      type (B2PlasmaSnapshot), intent (inout) :: psnl
*   ..common blocks

*-----------------------------------------------------------------------
*.documentation
*
*  1. purpose
*
*     B2NPCO_MOL performs one iteration of a relaxation process for the
*     solution of the continuity equation for one molecular species.
*
*-----------------------------------------------------------------------
*.declarations

*   ..local variables
      real (kind=R8) ::
     *   aa(mpg%nCmxNv), wrk0(nCv), pm(nCv), um(nCv), nm(nCv),
     *   rxg_npco_mol
      logical warning
      integer :: iCv
      integer, save :: ncall_b2npco_mol = 0
*   ..procedures
      external xerrab, ipgetr
      external b2xvsg, b2ursc, b2usco, b2upco
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
      call subini ('b2npco_mol')
*   ..test nCv, nFc
      call xertst (0.lt.nCv.and.0.lt.nFc, 'faulty argument nCv, nFc')
*   ..test rxf
      call xertst (0.0_R8.le.rxf.and.rxf.le.1.0_R8,
     . 'faulty argument rxf')

*   ..initizalization
      rxg_npco_mol = switch%b2npco_rxg
      nm = nm_ext(:,0)
      pm = nm*tm_ext(:,0)
      um = um_ext(:,0)

*   ..extensive tests on first few calls
      if (ncall_b2npco_mol.lt.3) then
*    ..test sign of con_mol
       call b2xvsg (nFc, dv%con_mol(1,0,0), 1, 'con00', '.ge.')
       call b2xvsg (nFc, dv%con_mol(1,0,1), 1, 'con01', '.ge.')
       call b2xvsg (nFc, dv%con_mol(1,0,2), 1, 'con02', '.ge.')
*    ..test sign of snc, snv
       call b2xvsg (nCv, sr_mol%sn_mol(1,0), 1, 'sn_c', '.ge.')
       call b2xvsg (nCv, sr_mol%sn_mol(1,1), 1, 'sn_v', '.le.')
*    ..test sign of pm, nm
       call b2xvsg (nCv, pm, 1, 'pm', '.gt.')
       call b2xvsg (nCv, nm, 1, 'nm', '.gt.')
      endif

      if(switch%b2mndt_style.ne.1) then                            
        call xerrab('b2mndt_style.ne.1 not yet for molecules')
      endif                                                          

* ..particle balance
*   ..compute residual
      call b2ursc (nCv, nFc, mpg, nm, sr_mol%sn_mol(:,:),
     &    dv%fn_mol, dv%resco_mol)
                                
      solvcreg = .true. !! TODO introduce new molecule-specific variables                   
*   ..compute correction
      call b2usco (nCv, nFc, nVx, 1, switch, geo, mpg, nregionv,    
     &   solvcreg, 1, rxg_npco_mol, nm, pm,
     &   sr_mol%sn_mol(:,1), dv%flo_mol(:,:), dv%con_mol(:,:,:),
     &   dv%resco_mol(:), dv%corp_mol, aa, 'b2usco__mol')
*   ..apply correction
      call b2upco (nCv, nFc, nVx, geo, mpg, switch, nregionv,
     &   solve_ua, solvcreg, rxf, switch%b2npco_pcm0, dv%corp_mol,
     &   dv%pcc_mol(:,0),
     &   pm, nm, um, 1, psnl)

*    ..compute time derivative
      ! to obtain the true physical derivative V(dn/dt), we need to divide by ts_factor, which was previously multiplied to snadt inside b2scdt
      dv%dndt_mol = (-sr_mol%sn_dt_mol(:,1)*nm
     &                   -sr_mol%sn_dt_mol(:,0))

*   ..update states
      nm_ext(:,0) = nm
      um_ext(:,0) = um 

*    ..output quantities for continuity equation
      if(switch%b2tfmol_iout.ne.0) then       
        call my_out_us(70,nCv,0,pm,'b2npco__mol_pm')
      endif

* ..return
      ncall_b2npco_mol = ncall_b2npco_mol+1
      call subend ()
      return

*-----------------------------------------------------------------------
*.end b2npco_mol

      end subroutine b2npco_mol
