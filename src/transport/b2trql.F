*-----------------------------------------------------------------------
*.specification

      subroutine b2trql (nCv, nFc, ns, switch, geo, mpg, pl, dv, st_ext, 
     .  chvemx, chvimx)
      use b2mod_types
      use b2mod_constants
      use b2mod_b2cmpa
      use b2mod_b2cmpt
      use b2mod_switches
      use b2us_geo
      use b2us_map
      use b2us_plasma
! csc The following are not necessary for computation but are needed
!     for adjoint AD to avoid side-effect variables
      use b2mod_ad
     & , only : ncall_b2trql, ncall_b2tlhe
      use b2mod_subsys
      implicit none
*   ..input arguments (unchanged on exit)
      integer nCv, nFc, ns
      type (switches), intent (in) :: switch
      type (geometry), intent (in) :: geo
      type (mapping), intent (in) :: mpg
      type (B2Plasma), intent (in) :: pl
      type (B2Derivatives), intent (in) :: dv
      type (B2StateExt), intent (in) :: st_ext
*   ..output arguments (unspecified on entry)
      real (kind=R8) ::
     *  chvemx(nFc), chvimx(nFc)
*   ..common blocks
*-----------------------------------------------------------------------
*.documentation
*
*  1. purpose
*
*     B2TRQL computes the flux limit velocities.
*
*
*  3. description (see also routine b2cdca)
*
*     This routine computes the following flux limit velocities:
*     - chvemx: electron heat transport.
*     - chvimx: ion heat transport.
*
*-----------------------------------------------------------------------
*.declarations

*   ..local variables
      real (kind=R8) ::
     *   cflme, cflmi
      real (kind=R8) ::
     *   nirm(nCv)
*   ..procedures
      real (kind=R8) ::
     *   smin, smax
      external xertst, smin, smax
      external b2tlhe, b2tlhi, b2tlmv, b2xpnr, b2xvsg
*   ..initialisation

*-----------------------------------------------------------------------
*.documentation-internal
*
*     The workspace array w is used as a stack for storage of dynamic
*     real*8 arrays. An array 'name' in w(0:*) will be stored starting at
*     location w(l'name').
*
*     Description of some local variables (including arrays in w(0:*):
*
*     nirm - (nCv) real*8 array.
*     For (iCv) in (1:nCv), nirm(iCv) specifies the atom
*     density weighted by 1/sqrt(am) on the (iCv) cell. It will hold
*     that 0.lt.nirm().
*     The defining expression for nirm is
*       nirm() = (sum is : 0.le.is.lt.ns : na(,is)/sqrt(am(is))).
*
*-----------------------------------------------------------------------
*.computation

* ..preliminaries
*   ..subprogram start-up calls
      call subini ('b2trql')
*   ..test nCv, nFc, ns
      call xertst (0.le.nCv.and.0.le.nFc, 'faulty argument nCv, nFc')
      call xertst (1.le.ns, 'faulty argument ns')
*   ..extensive tests on first few calls
      if (ncall_b2trql.lt.3) then
*    ..test bb
       call b2xvsg (nCv, geo%cvBb(1,3), 1, 'bb3', '.gt.')
*    ..test sign of cvVol
       call b2xvsg (nCv, geo%cvVol, 1, 'vol', '.gt.')
*    ..test range of cvQgam, fcQgam
       call xertst (-1.0_R8.lt.smin(nCv,geo%cvQgam(1,1),1).and.
     &   smax(nCv,geo%cvQgam(1,1),1).lt.1.0_R8, 
     &   'faulty argument range cvQgam(,1)')
       call xertst (0.0_R8.lt.smin(nCv,geo%cvQgam(1,0),1).and.
     &   smax(nCv,geo%cvQgam(1,0),1).le.1.0_R8, 
     &   'faulty argument range cvQgam(,0)')
       call xertst (-1.0_R8.lt.smin(nFc,geo%fcQgam(1,1),1).and.
     &   smax(nFc,geo%fcQgam(1,1),1).lt.1.0_R8, 
     &   'faulty argument range fcQgam(,1)')
       call xertst (0.0_R8.lt.smin(nFc,geo%fcQgam(1,0),1).and.
     &   smax(nFc,geo%fcQgam(1,0),1).le.1.0_R8, 
     &   'faulty argument range fcQgam(,0)')
*    ..test sign of na, te, ti, ne2
       call b2xvsg (nCv*ns, pl%na, 1, 'na', '.gt.')
       call b2xvsg (nCv, pl%te, 1, 'te', '.gt.')
       call b2xvsg (nCv, pl%ti, 1, 'ti', '.gt.')
       call b2xvsg (nCv, dv%ne2, 1, 'ne2', '.gt.')
      endif
*   ..compute nirm (atom density weighted by 1/sqrt(am)
      call b2xpnr (nCv, ns, pl%na, nirm)

* ..compute flux limit velocities
*   ..compute electron heat flux limit velocity
      cflme = cflim(0)
      call b2tlhe (nCv, nFc, me, cflme, switch, geo, mpg, dv%ne, 
     &  pl%te, chvemx)
*   ..compute all atom heat flux limit velocity
      cflmi = cflim(1)
      call b2tlhi (nCv, nFc, mp, cflmi, switch, geo, mpg, dv%ni,
     &  nirm, st_ext%ni, pl%ti, pl%tn, chvimx)

* ..return
      ncall_b2trql = ncall_b2trql+1
      call subend ()
      return

*-----------------------------------------------------------------------
*.end b2trql

      end subroutine b2trql
