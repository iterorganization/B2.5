*-----------------------------------------------------------------------
*.specification

      subroutine b2mwqt_mol (nput, nCv, itim, itcnt, switch,
     &  geo, pl, dv, diag)
      use b2mod_types
      use b2mod_constants
      use b2mod_b2cmpa
      use b2mod_switches
      use b2us_geo
      use b2us_plasma
      use b2mod_subsys
      use b2mod_fluid_molecules
      implicit none
*   ..input arguments (unchanged on exit)
      integer nput, nCv, ns, itim, itcnt
      type (switches), intent (in) :: switch
      type (geometry), intent (in) :: geo
      type (B2Plasma), intent (in) :: pl
      type (B2Derivatives), intent (in) :: dv
      type (B2Diagnostic), intent (in) :: diag
*   ..common blocks
*-----------------------------------------------------------------------
*.documentation
*
*  1. purpose
*
*     B2MWQT_MOL writes brief data concerning the numerical evolution.
*
*
*  3. description (see also routine b2cdca)
*
*     This routine writes to unit (nput) the 'quick trace' evolution
*     data. Unit (nput) must be connected to an un*formatted file.
*
*-----------------------------------------------------------------------
*.declarations

*   ..local variables
      integer nfldmx, nfld
      parameter (nfldmx=1000)
      integer k, is
      real (kind=R8) ::
     *   z(0:nfldmx-1)
*   ..procedures
      real (kind=R8) ::
     *   smin
      external xertst, ipgeti, smin, cfwure
*-----------------------------------------------------------------------
*.documentation-internal
*
*     Description of some local variables:
*
*     style - integer.
*     style specifies the amount of output produced on each call.
*     Legitimate values are in the range 0..4. style.eq.0 implies
*     no output, and the amount of output is an ascending function
*     of the value of style. See the code for precise information.

*-----------------------------------------------------------------------
*.computation

* ..preliminaries
*   ..subprogram start-up calls
      call subini ('b2mwqt_mol')
*   ..test nput
      call xertst (1.le.nput, 'faulty argument nput')
*   ..test itim, itcnt
      call xertst (0.le.itim.and.0.le.itcnt,
     &  'faulty argument itim, itcnt')
*   ..test residuals
      call xertst (0.0_R8.le.diag%aresco_mol.and.
     &  0.0_R8.le.diag%aresh_mol,
     & 'faulty argument aresco_mol / aresh_mol')
*   ..set and test nfld (number of fields)
      nfld = 7
      call xertst (nfld.le.nfldmx, 'faulty parameter nfldmx')

* ..get data for output
      k = 0
*    ..iteration counters
      z(k) = itim
      k = k+1
      z(k) = itcnt
      k = k+1
*    ..norms of residuals
      z(k) = diag%aresco_mol
      k = k+1
      z(k) = diag%aresmo_mol
      k = k+1
      z(k) = diag%aresmo_mol_rad
      k = k+1
      z(k) = diag%aresmo_mol_dia
      k = k+1
      z(k) = diag%aresh_mol
      k = k+1

*   ..test k
      call xertst (k.eq.nfld, 'faulty parameter nfld')

* ..produce output
      call cfwure (nput, nfld, z, 'data')

* ..return
      call subend ()
      return

*-----------------------------------------------------------------------
*.end b2mwqt_mol

      end subroutine b2mwqt_mol
