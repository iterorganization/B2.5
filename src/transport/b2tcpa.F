*-----------------------------------------------------------------------
*.specification

      subroutine b2tcpa (nCv, nFc, nVx, ns, switch, geo, mpg, te, po,
     *    ne, na, ua, rza, rz2, zeff, csig, calf, ehx, st_ext,
     &    fch_p, fch_pi_c, fch_pi_f)                                    !srv 13.01.17
      use b2mod_types
      use b2mod_constants
      use b2mod_b2cmpa
      use b2mod_switches
      use b2us_geo
      use b2us_map
      use b2us_plasma
! csc The following are not necessary for computation but are needed
!     for adjoint AD to avoid side-effect variables
      use b2mod_ad
     & , only : ncall_b2tcpa, ncall_b2xehx
      use b2mod_subsys
      implicit none
*   ..input arguments (unchanged on exit)
      integer nCv, nFc, nVx, ns
      type (switches), intent (in) :: switch
      type (geometry), intent (in) :: geo
      type (mapping), intent (in) :: mpg
      type (B2StateExt), intent (in) :: st_ext
      real (kind=R8) ::
     *   te(nCv), po(nCv), ne(nCv), na(nCv,0:ns-1), ua(nCv,0:ns-1),
     *   rza(nCv,0:ns-1), rz2(nCv,0:ns-1), zeff(nCv), csig(nFc,0:1),
     *   calf(nFc,0:1)
*   ..output arguments (unspecified on entry)
      real (kind=R8) ::
     &   fch_p(nFc,0:1), ehx(nFc), fch_pi_c(nCv,0:1), fch_pi_f(nFc,0:1) !srv 29.01.20
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
      real (kind=R8) ::
     *   tev(nVx), dtep(nFc), wrk(nFc,0:1)
*   ..procedures
      external xertst
      external b2xvsg, b2xvff
*   ..initialization
*-----------------------------------------------------------------------
*.computation

* ..preliminaries
*   ..subprogram start-up calls
      call subini ('b2tcpa')
*   ..test nCv, nFc
      call xertst (0.le.nCv.and.0.le.nFc, 'faulty argument nCv, nFc')
*   ..extensive tests on first few calls
      if (ncall_b2tcpa.lt.3) then
*    ..test sign of ne, te, csig, calf
       call b2xvsg (nCv, ne, 1, 'ne', '.gt.')
       call b2xvsg (nCv, te, 1, 'te', '.gt.')
       wrk = csig*geo%fcQalf
       call b2xvsg (nFc*2, wrk, 1, 'csig', '.ge.')
      endif

* ..compute parallel current
      if(switch%pot_eq.eq.1 .and. switch%no_current.eq.0) then         !srv 26.07.06 27.02.13
        call b2xehx (nCv, nFc, nVx, switch, geo, mpg, te, po, ne, ehx) !srv 13.01.17 29.01.20
*   ..compute poloidal differences of te
        call diff_p(nCv, nFc, nVx, 0, geo, mpg, te, tev, dtep)
*   ..compute fch_p
        fch_p(:,0) = (csig(:,0)*ehx - switch%b2sifr_phm2*
     &                               calf(:,0)*dtep)*switch%prl_cur
        if(switch%b2sigp_style.eq.2) then
          call b2xzef (nCv, ns, rz2, na, ne, zeff, st_ext)
          call b2xpic (nCv, ns, mpg, qe, rza, geo%cvBb, na, ua, zeff, 
     &      st_ext, fch_pi_c)
          call b2xpif (nCv, nFc, switch, geo, mpg, fch_pi_c, fch_pi_f)
        else
          fch_pi_f = 0.0_R8                                            !srv 04.04.17
        endif
        fch_p = fch_p - fch_pi_f                                       !srv 13.01.17 }
        fch_p(:,1) = 0.0e0_R8
      else                                                             !srv 26.07.06 {
        fch_p = 0.0e0_R8
        ehx   = 0.0_R8
      endif                                                            !srv 26.07.06 }

* ..return
      ncall_b2tcpa = ncall_b2tcpa+1
      call subend ()
      return

*-----------------------------------------------------------------------
*.end b2tcpa

      end subroutine b2tcpa
