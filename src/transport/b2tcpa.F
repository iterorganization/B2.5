*-----------------------------------------------------------------------
*.specification

      subroutine b2tcpa (nx, ny, ns, te, po, ne, na, ua, csig, calf,
     &          ehx, zeff, fch_p)                                        !srv 13.01.17 27.12.19
      use b2mod_types
      use b2mod_indirect
      use b2mod_constants
      use b2mod_b2cmpa
      use b2mod_geo
     & , only : hx, hy, vol, bb
      use b2mod_plasma
     & , only : fch_pi_f, fch_pi_c                                       !srv 13.01.17
      use b2mod_rates
      implicit none
*   ..input arguments (unchanged on exit)
      integer nx, ny, ns
      real (kind=R8) ::
     *   te(-1:nx,-1:ny), po(-1:nx,-1:ny), ne(-1:nx,-1:ny),
     &   na(-1:nx,-1:ny,0:ns-1), ua(-1:nx,-1:ny,0:ns-1),
     &   csig(-1:nx,-1:ny,0:1), calf(-1:nx,-1:ny,0:1)
*   ..output arguments (unspecified on entry)
      real (kind=R8) ::
     &   ehx(-1:nx,-1:ny), zeff(-1:nx,-1:ny),
     &   fch_p(-1:nx,-1:ny,0:1)
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
      integer ncall, n2, ix, iy
      integer, save :: pot_eq = 1                                        !srv 26.07.06
      integer, save :: no_current=0                                      !srv 27.02.13
      integer, save :: style = 2                                         !srv 13.01.17
      real (kind=R8), save :: prl_cur = 1.0_R8, phm2 = 1.0_R8            !srv 02.04.15 08.10.14
      real (kind=R8) ::
     *   tex
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
       call ipgeti ('b2sigp_style', style)                               !srv 13.01.17
      endif                                                              !srv
*   ..set n2
      n2 = (nx+2)*(ny+2)
*   ..test nx, ny
      call xertst (0.le.nx.and.0.le.ny, 'faulty argument nx, ny')
*   ..extensive tests on first few calls
      if (ncall.lt.3) then
*    ..test sign of ne, te, csig
       call b2xvsg (n2, ne, 1, 'ne', '.gt.')
       call b2xvsg (n2, te, 1, 'te', '.gt.')
       call b2xvsg (n2*2, csig, 1, 'csig', '.ge.')
*    ..test edge values of csig, calf
       call b2xvff (nx, ny, csig, 'csig')
       call b2xvff (nx, ny, calf, 'calf')
      endif

* ..compute parallel current
      if(pot_eq.eq.1 .and. no_current.eq.0) then                         !srv 26.07.06 27.02.13
        call b2xehx (nx, ny, hx, te, po, ne, ehx)                        !srv 13.01.17 {
*   ..compute fch_p
        do iy = -1, ny
          do ix = -1, nx                                                 !srv 26.11.02 {
            if(leftix(ix,iy).eq.-2)then                                  !srv 26.11.02
              fch_p(ix,iy,0)  = 0.0e0_R8
            else                                                         !srv 26.11.02 }
              tex = te(ix,iy)-te(leftix(ix,iy),leftiy(ix,iy))
              fch_p(ix,iy,0) = (csig(ix,iy,0)*ehx(ix,iy)-
     &                     phm2*calf(ix,iy,0)*tex)*prl_cur               !srv 03.02.03 08.10.14
            endif                                                        !srv 26.11.02
          enddo
        enddo
        if(style.eq.2) then
          call b2xzef (nx, ny, ns, rz2, na, ne, zeff)
          call b2xpic (nx, ny, ns, qe, rza, bb, na, ua, zeff, fch_pi_c)
          call b2xpif (nx, ny, hx, hy, vol, fch_pi_c, fch_pi_f)
        else
          fch_pi_f = 0.0_R8                                              !srv 04.04.17
        endif
        fch_p = fch_p - fch_pi_f                                         !srv 13.01.17 }
        fch_p(:,:,1) = 0.0e0_R8
      else                                                               !srv 26.07.06 {
        fch_p = 0.0e0_R8
      endif                                                              !srv 26.07.06 }

* ..return
      ncall = ncall+1
      call subend ()
      return

*-----------------------------------------------------------------------
*.end b2tcpa

      end subroutine b2tcpa
