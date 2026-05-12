*-----------------------------------------------------------------------
*.specification
      subroutine b2tfh_mol (nCv, nFc, nVx, switch, geo, mpg,  
     &   pl, dv, co)                                           
      use b2mod_types
      use b2mod_math
      use b2mod_constants
      use b2mod_b2cmfs
      use b2mod_b2cmpa
      use b2mod_switches
      use b2us_geo
      use b2us_map
      use b2us_plasma
      use b2mod_subsys
      use b2mod_fluid_molecules
      implicit none
*   ..input arguments (unchanged on exit)
      integer nCv, nFc, nVx, ns, ismain
      type (switches), intent (in) :: switch
      type (geometry), intent (in) :: geo
      type (mapping), intent (in) :: mpg
      type(B2Plasma), intent (in) :: pl
      logical balance                                                   !djm Jan2017
*   ..input/output arguments (unspecified on entry)
      type(B2Coeff), intent (inout) :: co
      type(B2Derivatives), intent (inout) :: dv
*   ..common blocks
*-----------------------------------------------------------------------
*.documentation
*
*  1. purpose
*
*     B2TFH_MOL computes the molecular heat flux between neighboring cells,
*     and the coefficients that specify its linearization.
*
*-----------------------------------------------------------------------
*.declarations

*   ..local variables
      integer  meth
      real (kind=R8) ::
     *   flo0(nFc,0:1), con0(nFc,0:1), 
     *   wght(nFc,2), dumm0(nFc,0:1),
     *   dumm1(nFc,0:1), dumm2(nFc,0:1), dumm3(nFc,0:1)
      real (kind=R8) ::
     *   nmf(nFc), tmf(nFc), wrkf(nFc,0:1)
*   ..procedures
      external xertst
      external b2xvsg
*   ..initialization

*-----------------------------------------------------------------------
*.computation

* ..preliminaries
*   ..subprogram start-up calls
      call subini ('b2tfh_mol')
*   ..set internal parameters on first call                                                                
*   ..test nCv, nFc
      call xertst (0.lt.nCv.and.0.lt.nFc, 'faulty argument nCv, nFc')

      meth = switch%b2tfhi_discr_meth

*   ..set dummy arrays
      dumm0 = 0.0_R8
      dumm1 = 0.0_R8

*   ..interpolate nm and tm to cell faces

      call intface(nCv,nFc,mpg%fcCv,wght,tm_ext(:,0),tmf)
      call intface(nCv,nFc,mpg%fcCv,wght,nm_ext(:,0),nmf)


*   ..convective coefficients
      !flo0
      flo0(:,0) = 1.5_R8*dv%fn_mol(:,0)  
      flo0(:,1) = 1.5_R8*dv%fn_mol(:,1)

*   ..apply discretization scheme
      call calcflow(nCv,nFc,nVx,meth,geo,mpg,tm_ext(:,0),flo0,
     &  co%chc_mol,dv%fh_mol,dumm1,dumm2)

*   ..compute flo, con
      call calccoef(nCv,nFc,nVx,meth,geo,flo0,co%chc_mol,
     &  flo0,con0)

!! TODO add linearization coefficients for molecules?
      dv%floh_mol = dv%floh_mol + flo0
      dv%conh_mol(:,:,0) = dv%conh_mol(:,:,0) + con0


      if(switch%b2tfmol_iout.ne.0) then                                  
       call my_out_us(70,nFc,1,dv%fh_mol(:,0),'b2tfh_mol_th')
       call my_out_us(70,nFc,1,dv%fh_mol(:,1),'b2tfh_mol_r')
      endif                                                    

* ..return
      call subend ()
      return


*.end b2tfh_mol

      end subroutine b2tfh_mol
