*-----------------------------------------------------------------------
*.specification

      subroutine b2xehx (nCv, nFc, nVx, switch, geo, mpg, te, po, 
     &                   ne, ehx)       !srv 13.01.17
      use b2mod_types
      use b2mod_constants
      use b2us_geo
      use b2us_map
      use b2mod_switches
! csc The following are not necessary for computation but are needed
!     for adjoint AD to avoid side-effect variables
      use b2mod_ad
     & , only : ncall_b2xehx
      implicit none
*   ..input arguments (unchanged on exit)
      integer nCv, nFc, nVx
      type (switches), intent (in) :: switch
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
      real (kind=R8) ::
     *   dpop(nFc), wrk0(nFc), 
     *   wrk1(nFc), wrk2(nFc), wrk(nCv),
     *   wrkv(nVx), ones(nFc,2)
*   ..procedures
      external subini, subend, xertst
      external b2xvsg
*   ..initialization
*-----------------------------------------------------------------------
*.computation

* ..preliminaries
*   ..subprogram start-up calls
      call subini ('b2xehx')
*   ..test nCv, nFc
      call xertst (0.le.nCv.and.0.le.nFc, 'faulty argument nCv, nFc')
*   ..extensive tests on first few calls
      if (ncall_b2xehx.lt.3) then
*    ..test sign of ne, te, csig
       call b2xvsg (nCv, ne, 1, 'ne', '.gt.')
       call b2xvsg (nCv, te, 1, 'te', '.gt.')
      endif

*   ..compute ehx
      if(switch%pot_eq.ne.1 .and. switch%no_current.ne.0) then    !srv 29.01.20 {
        ehx = 0.0_R8
      else
        if(switch%mode_ehx.eq.0) then
          ones = 1.0_R8
*       ..compute poloidal differences of po and pe
          call diff_p(nCv, nFc, nVx, 0, geo, mpg, po, wrkv, dpop)
          call diff_p(nCv, nFc, nVx, 0, geo, mpg, ne*te, wrkv, wrk0)
*       ..interpolate ne to faces
          call intface(nCv,nFc,mpg%fcCv,ones,ne,wrk1)
*       ..compute ehx      
          ehx = -dpop + wrk0/(qe*wrk1)
        else if (switch%mode_ehx.eq.1) then
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
      ncall_b2xehx = ncall_b2xehx+1
      call subend ()
      return

*-----------------------------------------------------------------------
*.end b2xehx

      end
