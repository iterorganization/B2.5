*-----------------------------------------------------------------------
*.specification

      subroutine b2srdt_mol (nCv, ns, dtim, switch, geo, mpg, 
     .                       sr_mol)             
      use b2mod_types
      use b2mod_numerics_namelist
      use b2mod_constants
      use b2mod_b2cmpa
      use b2mod_switches
      use b2us_geo
      use b2us_map
      use b2us_plasma
      use b2mod_fluid_molecules
      use b2mod_subsys
      implicit none
*   ..input arguments (unchanged on exit)
      logical :: lout                                                   !srv 16.07.10
      integer nCv, ns
      type (switches), intent (in) :: switch
      type (geometry), intent (in) :: geo
      type (mapping), intent (in) :: mpg
      real (kind=R8) ::
     *  dtim
*   ..input/output arguments
      type (B2SourceMol), intent (inout) :: sr_mol
*   ..common blocks
*-----------------------------------------------------------------------
*.documentation
*
*  1. purpose
*
*     B2SRDT_MOL modifies the source coefficients to account for the
*     implicit time discretisation for fluid molecules.
*
*
*-----------------------------------------------------------------------
*.declarations

*   ..local variables
      integer iCv
      character chk*1                                       
      real (kind=R8) :: t0, ttim
*   ..procedures
      external xertst
      external b2xvsg
*   ..initialisation

*-----------------------------------------------------------------------
*.computation

* ..preliminaries
*   ..subprogram start-up calls
      call subini ('b2srdt_mol')
*   ..test nCv
      call xertst (0.lt.nCv, 'faulty argument nCv')
*   ..test sign of dtim
      call xertst (0.0_R8.lt.dtim, 'faulty argument dtim')


* ..include the contributions from dtim
*   ..modify sn_dt_mol

      do iCv = 1, mpg%nCi
        ttim=dtim*time_factor(iCv)
        t0=(switch%b2srdt_phm0/ttim)*geo%cvVol(iCv)*ts_factor    
        sr_mol%sn_dt_mol(iCv,0) = t0*nm_ext(iCv,0)
        sr_mol%sn_dt_mol(iCv,1) = -t0
      enddo

*   ..modify sm_dt_mol

      do iCv = 1, mpg%nCi
        ttim=dtim*time_factor(iCv)
        t0 = (switch%b2srdt_phm1/ttim)*geo%cvVol(iCv)*geo%cvHz(iCv) 
     &         *ts_factor
        sr_mol%sm_dt_mol(iCv,0) = t0*am_d2*mp*nm_ext(iCv,0)*
     &                            um_ext(iCv,0)
        sr_mol%sm_dt_mol(iCv,1) =-t0*am_d2*mp*nm_ext(iCv,0)
      enddo

*   ..modify sh_dt_mol

      do iCv = 1, mpg%nCi
        ttim=dtim*time_factor(iCv)
        t0 = (switch%b2srdt_phm3/ttim)*geo%cvVol(iCv)*ts_factor
        sr_mol%sh_dt_mol(iCv,0) = 1.5_R8*t0*nm_ext(iCv,0)*tm_ext(iCv,0)
        sr_mol%sh_dt_mol(iCv,1) =-1.5_R8*t0*nm_ext(iCv,0)
      enddo

*   ..add to total sources
       sr_mol%sn_mol = sr_mol%sn_mol + sr_mol%sn_dt_mol
       sr_mol%sm_mol = sr_mol%sm_mol + sr_mol%sm_dt_mol
       sr_mol%sh_mol = sr_mol%sh_mol + sr_mol%sh_dt_mol

* ..return
      call subend ()
      return

*-----------------------------------------------------------------------
*.end b2srdt_mol

      end subroutine b2srdt_mol
