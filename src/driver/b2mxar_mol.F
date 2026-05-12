*-----------------------------------------------------------------------
*.specification

      subroutine b2mxar_mol (nCv, switch, geo, mpg, pl, dv, diag)
      use b2mod_types
      use b2mod_constants
      use b2mod_b2cmpa
      use b2mod_switches
      use b2us_geo
      use b2us_map
      use b2us_plasma
      use b2mod_subsys
      use b2mod_fluid_molecules
      implicit none
*   ..input arguments (unchanged on exit)
      integer nCv, ns
      type (switches), intent (in) :: switch
      type (geometry), intent (in) :: geo
      type (mapping), intent (in) :: mpg
      type (B2Plasma), intent (in) :: pl
      type (B2Derivatives), intent (in) :: dv
*   ..output arguments (unspecified on entry)
      type (B2Diagnostic), intent (inout) :: diag
*   ..common blocks
*-----------------------------------------------------------------------
*.documentation
*
*  1. purpose
*
*     B2MXAR_MOL computes scaled (L1-)norms of the residuals for 
*     fluid molecules.
*
*-----------------------------------------------------------------------
*.declarations

*   ..local variables
      integer iCv
      real (kind=R8) ::
     *  navi, movi, prvm, t0, t1
*   ..procedures
      intrinsic sqrt, abs
      external xertst
*   ..initialisation

*-----------------------------------------------------------------------
*.computation

* ..preliminaries
*   ..subprogram start-up calls
      call subini ('b2mxar_mol')
*   ..set internal parameters on first call
*   ..test nCv
      call xertst (0.lt.nCv, 'faulty argument nCv')

* ..compute scaled norms
*   ..compute aresco, aresmo
      navi = 0.0_R8
      movi = 0.0_R8
      do iCv = 1, mpg%nCi
          navi = navi+geo%cvVol(iCv)*nm_ext(iCv,0)
          movi = movi+geo%cvVol(iCv)*
     &      nm_ext(iCv,0)*sqrt(tm_ext(iCv,0)*am_d2*mp)

      enddo
      t0 = 0.0_R8
      t1 = 0.0_R8
      do iCv = 1, mpg%nCi
       t0 = t0+abs(dv%resco_mol(iCv))
       t1 = t1+abs(dv%resmo_mol(iCv))
      enddo
      diag%aresco_mol = t0/navi
      diag%aresmo_mol = t1/movi


*   ..compute prvm
      prvm = 0.0_R8
      do iCv = 1, mpg%nCi
        prvm = prvm+geo%cvVol(iCv)*nm_ext(iCv,0)*tm_ext(iCv,0)
      enddo
*   ..compute aresh_mol
      t0 = 0.0_R8

      do iCv = 1, mpg%nCi
       t0 = t0+abs(dv%resh_mol(iCv))
      enddo
      diag%aresh_mol = t0/(1.5_R8*prvm)

* ..return
      call subend ()
      return

*-----------------------------------------------------------------------
*.end b2mxar_mol

      end subroutine b2mxar_mol
