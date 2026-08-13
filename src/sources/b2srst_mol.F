*-----------------------------------------------------------------------
*.specification

      subroutine b2srst_mol (nCv, switch, sr_mol)
      use b2mod_types
      use b2mod_constants
      use b2mod_b2cmpa
      use b2mod_switches
      use b2us_plasma
      use b2mod_subsys
      use b2mod_fluid_molecules
      implicit none
*   ..input arguments (unchanged on exit)
      integer nCv                               
      type (switches), intent (in) :: switch
*   ..input/output arguments
      type (B2SourceMol), intent (inout) :: sr_mol
*   ..common blocks
*-----------------------------------------------------------------------
*.documentation
*
*  1. purpose
*
*     B2SRST_MOL modifies the source coefficients to provide numerical
*     stabilisation.
*
*
*-----------------------------------------------------------------------
*.declarations

*   ..local variables
      integer iCv
      real (kind=R8) :: t0
*   ..procedures
      intrinsic min, max
      external xertst, b2xvsg
*   ..initialisation
*-----------------------------------------------------------------------
*.computation

* ..preliminaries
*   ..subprogram start-up calls
      call subini ('b2srst_mol')
*   ..set internal parameters on first call
*   ..test nCv
      call xertst (0.lt.nCv, 'faulty argument nCv')

*  ..stabilise
*   ..modify sn

      do iCv = 1, nCv
        t0 = sr_mol%sn_mol(iCv,0)+sr_mol%sn_mol(iCv,1)*n_mol(iCv,0)
        sr_mol%sn_mol(iCv,0) =
     &     max(sr_mol%sn_mol(iCv,0),(1.0_R8+switch%b2srst_rf0)*t0,
     &     -switch%b2srst_rf0*t0)
        sr_mol%sn_mol(iCv,1) = (1.0_R8/n_mol(iCv,0))*
     &     min(sr_mol%sn_mol(iCv,1)*n_mol(iCv,0),-switch%b2srst_rf0*t0,
     &     (1.0_R8+switch%b2srst_rf0)*t0)
      enddo

*   ..modify sm

      do iCv = 1, nCv
*      ..make sm(,,1) negative
        sr_mol%sm_mol(iCv,0) = sr_mol%sm_mol(iCv,0)+
     &     max(0.0_R8,sr_mol%sm_mol(iCv,1))*u_mol(iCv,0)
        sr_mol%sm_mol(iCv,1) = min(sr_mol%sm_mol(iCv,1),0.0_R8)

        sr_mol%sm_mol_rad(iCv,0) = sr_mol%sm_mol_rad(iCv,0)+
     &     max(0.0_R8,sr_mol%sm_mol_rad(iCv,1))*u_rad_mol(iCv,0)
        sr_mol%sm_mol_rad(iCv,1) = min(sr_mol%sm_mol_rad(iCv,1),0.0_R8)

        sr_mol%sm_mol_dia(iCv,0) = sr_mol%sm_mol_dia(iCv,0)+
     &     max(0.0_R8,sr_mol%sm_mol_dia(iCv,1))*u_dia_mol(iCv,0)
        sr_mol%sm_mol_dia(iCv,1) = min(sr_mol%sm_mol_dia(iCv,1),0.0_R8)
      enddo


*   ..modify sh

      do iCv = 1, nCv
*     ..modify sh(,0:1)
        t0 = sr_mol%sh_mol(iCv,0)+sr_mol%sh_mol(iCv,1)*t_mol(iCv,0)
        sr_mol%sh_mol(iCv,0) = 
     &        max(sr_mol%sh_mol(iCv,0),(1.0_R8+switch%b2srst_rf1)*t0,
     &                  -switch%b2srst_rf1*t0)
        sr_mol%sh_mol(iCv,1) = (1.0_R8/t_mol(iCv,0))*
     &    min(sr_mol%sh_mol(iCv,1)*t_mol(iCv,0),-switch%b2srst_rf1*t0,
     &    (1.0_R8+switch%b2srst_rf1)*t0)
      enddo


* ..return
      call subend ()
      return

*-----------------------------------------------------------------------
*.end b2srst_mol

      end subroutine b2srst_mol
