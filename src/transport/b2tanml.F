*-----------------------------------------------------------------------
*.specification

      subroutine b2tanml (nCv, nFc, ns, switch, geo, mpg,
     & csig_an, po, fchanml)                                            !srv 20.09.06 {
      use b2mod_types
      use b2mod_constants
      use b2mod_b2cmpa
      use b2mod_switches
      use b2us_geo
      use b2us_map
      use b2mod_subsys
      implicit none
*   ..input arguments (unchanged on exit)
      integer, intent (in) :: nCv, nFc, ns
      type (switches), intent (in) :: switch
      type (geometry), intent (in) :: geo
      type (mapping), intent (in) :: mpg
      real (kind=R8) ::
     & po(nCv), csig_an(nFc,0:1)
*   ..output arguments (unspecified on entry)
      real (kind=R8) ::
     &    fchanml(nFc,0:1)
*   ..common blocks
*-----------------------------------------------------------------------
*.documentation
*
*  1. purpose
*
*     B2TANML computes the poloidal and radial components of the anomalous current.
*
*
*  8. further comments
*
*-----------------------------------------------------------------------
*.declarations

*   ..global named constants

*   ..local variables
      real (kind=R8) :: dpo(nFc,0:1), pov(mpg%nVx)
*   ..procedures
      external xertst
      external b2xvsg, diff_r
*-----------------------------------------------------------------------
*.computation

* ..preliminaries
*   ..subprogram start-up calls
      call subini ('b2tanml')
*   ..test nCv, nFc, ns
      call xertst (0.lt.nCv.and.0.lt.nFc, 'faulty argument nCv, nFc')
      call xertst (1.le.ns, 'faulty argument ns')

      if(switch%no_current.ne.0.or.
     &   switch%fch_anomalous.eq.0.0_R8) then                           !srv 27.02.13
        fchanml = 0.0_R8
      else                                                              !srv 27.02.13

*   ..compute differences of po
        call diff(nCv, nFc, mpg%nVx, 0, geo, mpg, po, pov, dpo)

*   ..compute the current
        fchanml = -csig_an*dpo

      endif                                                             !srv 27.02.13

* ..return
      call subend ()
      return

*-----------------------------------------------------------------------
*.end b2tanml

      end subroutine b2tanml                                             !srv 20.09.06 }
