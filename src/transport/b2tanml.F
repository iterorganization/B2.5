*-----------------------------------------------------------------------
*.specification

      subroutine b2tanml (nCv, nFc, ns, geo, mpg, csig_an, po,       !srv 20.09.06 {
     &                    fchanml)
      use b2mod_types
      use b2mod_constants
      use b2mod_b2cmpa
      use b2us_geo
      use b2us_map
      implicit none
*   ..input arguments (unchanged on exit)
      integer nCv, nFc, ns
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
*     B2TANML computes the x,y-components of anomalous current
*
*
*  8. further comments
*
*-----------------------------------------------------------------------
*.declarations

*   ..global named constants
      
*   ..local variables
      integer ncall, iCv, iFc
      integer, save :: no_current = 0                          !srv 27.02.13
      integer, save :: pot_eq
      integer, save :: guard_flows
      real (kind=R8) :: fch_anomalous, dpo(nFc,0:1),
     *  wght(nCv), pov(mpg%nVx)
      save ncall, fch_anomalous
*   ..procedures
      intrinsic max, abs
      external subini, subend, xertst
      external b2xvsg, b2xvff, b2xvfx, diff_r
*   ..initialization
      data ncall/0/, fch_anomalous/1.0e0_R8/, pot_eq/1/, guard_flows/2/
*-----------------------------------------------------------------------
*.computation

* ..preliminaries
*   ..subprogram start-up calls
      call subini ('b2tanml')
*   ..set internal parameters on first call
      if (ncall.eq.0) then
       call ipgetr ('b2tfhe_anomalous', fch_anomalous)
       if (fch_anomalous.eq.1.0_R8) 
     >  call ipgetr ('b2tanml_anomalous', fch_anomalous) ! for backward compatibility
       call ipgeti ('b2tfhe_no_current', no_current)                !srv 27.02.13
       call ipgeti ('b2news_poteq', pot_eq)                         !srv 27.02.13
       call xertst (pot_eq.eq.1.or.no_current.eq.1,
     &  'If b2news_poteq.ne.1, then b2tfhe_no_current must be 1!')
      endif
*   ..test nCv, nFc
      call xertst (0.le.nCv.and.0.le.nFc, 'faulty argument nCv, nFc')
*   ..extensive tests on first few calls
      if (ncall.lt.3) then
      endif

      if(no_current.ne.0.or.fch_anomalous.eq.0.0e0_R8) then                                !srv 27.02.13
        fchanml = 0.0e0_R8
      else                                                    !srv 27.02.13

*   ..compute differences of po
        wght = 1.0_R8
        call diff(nCv, nFc, mpg%nVx, 0, geo, mpg, po, pov, dpo)

*   ..compute the current
        fchanml = -csig_an*dpo

      endif                                                   !srv 27.02.13

* ..return
      ncall = ncall+1
      call subend ()
      return

*-----------------------------------------------------------------------
*.end b2tanml

      end                                                             !srv 20.09.06 }
