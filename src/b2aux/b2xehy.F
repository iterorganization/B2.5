*-----------------------------------------------------------------------
*.specification

      subroutine b2xehy (nCv, nFc, nVx, switch, geo, mpg, te, po,
     &                   ne, ehy)       !srv 13.01.17
      use b2mod_types
      use b2mod_constants
      use b2mod_switches
      use b2us_geo
      use b2us_map
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
     &   ehy(nFc)
*   ..common blocks
*-----------------------------------------------------------------------
*.documentation
*
*  1. purpose
*
*     B2XEHY computes the modified radial electrical field on all cell
*     faces.
*
*  8. further comments
*
*-----------------------------------------------------------------------
*.declarations

*   ..global named constants
      
*   ..local variables
      integer ncall
      real (kind=R8) ::
     *   dpor(nFc), wrk0(nFc), 
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
      call subini ('b2xehy')
*   ..test nCv, nFc
      call xertst (0.le.nCv.and.0.le.nFc, 'faulty argument nCv, nFc')
*   ..extensive tests on first few calls
      if (ncall.lt.3) then
*    ..test sign of ne, te, csig
       call b2xvsg (nCv, ne, 1, 'ne', '.gt.')
       call b2xvsg (nCv, te, 1, 'te', '.gt.')
      endif

*   ..compute ehy
      if(switch%pot_eq.ne.1 .and. switch%no_current.ne.0) then           !srv 29.01.20 {
        ehy = 0.0_R8
      else
        if(switch%mode_ehy.eq.0) then
          ones = 1.0_R8
*       ..compute radial differences of po and pe
          call diff_r(nCv, nFc, nVx, 0, geo, mpg, po, wrkv, dpor)
          call diff_r(nCv, nFc, nVx, 0, geo, mpg, ne*te, wrkv, wrk0)
*       ..interpolate ne to faces
          call intface(nCv,nFc,mpg%fcCv,ones,ne,wrk1)
*       ..compute ehy      
          ehy = -dpor + wrk0/(qe*wrk1)
        else if (switch%mode_ehy.eq.1) then
          call diff_r(nCv, nFc, nVx, 0, geo, mpg, po, wrkv, dpor)
          call diff_r(nCv, nFc, nVx, 0, geo, mpg, ne*te, wrkv, wrk0)
*       ..interpolate ne to faces
cwdk small loss of backwards compatibility:
cwdk now inverse distance weighting based on full distance cell center to face center ("hy")
cwdk instead of radial distance ("hy1_")
          call intface(nCv,nFc,mpg%fcCv,geo%fcHc,ne,wrk1)
*       ..compute ehy      
          ehy = -dpor + wrk0/(qe*wrk1)
        else
          ones = 1.0_R8
*       ..compute radial differences of po, te and log(ne) on faces
          call diff_r(nCv, nFc, nVx, 0, geo, mpg, po, wrkv, dpor)
          call diff_r(nCv, nFc, nVx, 0, geo, mpg, te, wrkv, wrk0)
          call diff_r(nCv, nFc, nVx, 0, geo, mpg, log(ne), wrkv, wrk1)
*       ..interpolate te to faces
          call intface(nCv,nFc,mpg%fcCv,ones,ne,wrk1)
*       ..compute ehy      
          ehy = -dpor + (wrk0 + wrk1*wrk2)/qe
        endif
      endif                                                              !srv 29.01.20 }

* ..return
      ncall = ncall+1
      call subend ()
      return

*-----------------------------------------------------------------------
*.end b2xehy

      end
