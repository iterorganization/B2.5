*-----------------------------------------------------------------------
*.specification

      subroutine b2xehx (nCv, nFc, nVx, geo, mpg, te, po, ne, ehx)       !srv 13.01.17
      use b2mod_types
      use b2mod_constants
      use b2us_geo
      use b2us_map
      implicit none
*   ..input arguments (unchanged on exit)
      integer nCv, nFc, nVx
      type (geometry), intent (in) :: geo
      type (mapping), intent (in) :: mpg
      real (kind=R8) ::
     &   te(nCv), po(nCv), ne(nCv)
*   ..output arguments (unspecified on entry)
      real (kind=R8) ::
     &   ehx(nFc)
*   ..common blocks
*-----------------------------------------------------------------------
*.documentation
*
*  1. purpose
*
*     B2XEHX computes the modified poloidal electrical field on all cell
*     faces.
*
*  8. further comments
*
*-----------------------------------------------------------------------
*.declarations

*   ..global named constants
      
*   ..local variables
      integer ncall
      integer, save :: mode_ehx = 0                                     !srv 06.04.15
      integer, save :: pot_eq = 1, no_current= 0                        !srv 29.01.20
      real (kind=R8) ::
     *   dpop(nFc), wrk0(nFc), 
     *   wrk1(nFc), wrk2(nFc), wrk(nCv),
     *   wrkv(nVx), ones(nFc,2)
      save ncall
*   ..procedures
      external subini, subend, xertst
      external b2xvsg
*   ..initialization
      data ncall/0/
*-----------------------------------------------------------------------
*.computation

* ..preliminaries
*   ..subprogram start-up calls
      call subini ('b2xehx')
*   ..set internal parameters on first call
      if (ncall.eq.0) then
       call ipgeti ('b2tfhe_mode_ehx',mode_ehx)                          !srv 06.04.15
       call ipgeti ('b2news_poteq', pot_eq)                              !srv 29.01.20 {
       call ipgeti ('b2tfhe_no_current',no_current)               
       call xertst (pot_eq.eq.1.or.no_current.eq.1,
     &  'If b2news_poteq.ne.1, then b2tfhe_no_current must be 1!')       !srv 29.01.20 }
      endif
*   ..test nCv, nFc
      call xertst (0.le.nCv.and.0.le.nFc, 'faulty argument nCv, nFc')
*   ..extensive tests on first few calls
      if (ncall.lt.3) then
*    ..test sign of ne, te, csig
       call b2xvsg (nCv, ne, 1, 'ne', '.gt.')
       call b2xvsg (nCv, te, 1, 'te', '.gt.')
      endif

*   ..compute ehx
      if(pot_eq.ne.1 .and. no_current.ne.0) then                         !srv 29.01.20 {
        ehx = 0.0_R8
      else
        if(mode_ehx.eq.0) then
          ones = 1.0_R8
*       ..compute poloidal differences of po and pe
          call diff_p(nCv, nFc, nVx, 0, geo, mpg, po, wrkv, dpop)
          call diff_p(nCv, nFc, nVx, 0, geo, mpg, ne*te, wrkv, wrk0)
*       ..interpolate ne to faces
          call intface(nCv,nFc,mpg%fcCv,ones,ne,wrk1)
*       ..compute ehx      
          ehx = -dpop + wrk0/(qe*wrk1)
        else if (mode_ehx.eq.1) then
*       ..compute poloidal differences of po and pe
          call diff_p(nCv, nFc, nVx, 0, geo, mpg, po, wrkv, dpop)
          call diff_p(nCv, nFc, nVx, 0, geo, mpg, ne*te, wrkv, wrk0)
*       ..interpolate ne to faces
          call intface(nCv,nFc,mpg%fcCv,geo%fcHc,ne,wrk1)
*       ..compute ehx      
          ehx = -dpop + wrk0/(qe*wrk1)
        else
          ones = 1.0_R8
*       ..compute poloidal differences of po, te and log(ne) on faces
          call diff_p(nCv, nFc, nVx, 0, geo, mpg, po, wrkv, dpop)
          call diff_p(nCv, nFc, nVx, 0, geo, mpg, te, wrkv, wrk0)
          call diff_p(nCv, nFc, nVx, 0, geo, mpg, log(ne), wrkv, wrk1)
*       ..interpolate te to faces
          call intface(nCv,nFc,mpg%fcCv,ones,te,wrk2)
*       ..compute ehx      
          ehx = -dpop + (wrk0 + wrk1*wrk2)/qe
        endif
      endif                                                              !srv 29.01.20 }

* ..return
      ncall = ncall+1
      call subend ()
      return

*-----------------------------------------------------------------------
*.end b2xehx

      end
