*-----------------------------------------------------------------------
*.specification

      subroutine b2tcpa (nCv, nFc, nVx, geo, mpg, te, po, ne, csig,
     *    calf, ehx, fch_p)  !srv 13.01.17
      use b2mod_types
      use b2mod_constants
      use b2mod_b2cmpa
      use b2us_geo
      use b2us_map
      implicit none
*   ..input arguments (unchanged on exit)
      integer nCv, nFc, nVx
      type (geometry), intent (in) :: geo
      type (mapping), intent (in) :: mpg
      real (kind=R8) ::
     *   te(nCv), po(nCv), ne(nCv),
     *   csig(nFc,0:1), calf(nFc,0:1), ehx(nFc)
*   ..output arguments (unspecified on entry)
      real (kind=R8) ::
     &   fch_p(nFc,0:1)
*   ..common blocks
*-----------------------------------------------------------------------
*.documentation
*
*  1. purpose
*
*     B2TCPA computes the parallel current between neighboring cells.
*
*  8. further comments
*
*-----------------------------------------------------------------------
*.declarations

*   ..global named constants
      
*   ..local variables
      integer ncall
      integer, save :: pot_eq = 1                                     !srv 26.07.06
      integer, save :: no_current=0                                   !srv 27.02.13
      real (kind=R8), save :: prl_cur = 1.0_R8, phm2 = 1.0_R8         !srv 02.04.15 08.10.14
      real (kind=R8) ::
     *   tev(nVx), dtep(nFc)
      save ncall
*   ..procedures
      external subini, subend, xertst
      external b2xvsg, b2xvff
*   ..initialization
      data ncall/0/
*-----------------------------------------------------------------------
*.computation

* ..preliminaries
*   ..subprogram start-up calls
      call subini ('b2tcpa')
*   ..set internal parameters on first call                              !srv
      if (ncall.eq.0) then                                               !srv
       call ipgeti ('b2news_poteq', pot_eq)                              !srv 26.07.06
       call ipgeti ('b2tfhe_no_current',no_current)                      !srv 27.02.13
       call ipgetr ('b2sifr_phm2', phm2)                                 !srv 08.10.14
       call ipgetr ('b2tfhe_prl_cur', prl_cur)                           !srv 02.04.15
       call xertst (pot_eq.eq.1.or.no_current.eq.1,
     &  'If b2news_poteq.ne.1, then b2tfhe_no_current must be 1!')
      endif                                                              !srv
*   ..test nCv, nFc
      call xertst (0.le.nCv.and.0.le.nFc, 'faulty argument nCv, nFc')
*   ..extensive tests on first few calls
      if (ncall.lt.3) then
*    ..test sign of ne, te, csig, calf
       call b2xvsg (nCv, ne, 1, 'ne', '.gt.')
       call b2xvsg (nCv, te, 1, 'te', '.gt.')
       call b2xvsg (nFc*2, csig*geo%fcQalf, 1, 'csig', '.ge.')
      endif

* ..compute parallel current
      if(pot_eq.eq.1 .and. no_current.eq.0) then                       !srv 26.07.06 27.02.13
*   ..compute poloidal differences of te
        call diff_p(nCv, nFc, nVx, 0, geo, mpg, te, tev, dtep)
*   ..compute fch_p
        fch_p(:,0) = (csig(:,0)*ehx - phm2*calf(:,0)*dtep)*prl_cur
        fch_p(:,1) = 0.0e0_R8
      else                                                             !srv 26.07.06 {
        fch_p = 0.0e0_R8
      endif                                                            !srv 26.07.06 }

* ..return
      ncall = ncall+1
      call subend ()
      return

*-----------------------------------------------------------------------
*.end b2tcpa

      end
