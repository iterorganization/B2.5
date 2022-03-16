*-----------------------------------------------------------------------
*.specification

      subroutine b2xehx (nx, ny, hx, hy, qz, qs, te, po, ne, ehx)        !srv 13.01.17
      use b2mod_types
      use b2mod_indirect
      use b2mod_constants
      use b2mod_subsys
      implicit none
*   ..input arguments (unchanged on exit)
      integer nx, ny
      real (kind=R8) ::
     &   hx(-1:nx,-1:ny), hy(-1:nx,-1:ny), qz(-1:nx,-1:ny,0:1),
     &   te(-1:nx,-1:ny), po(-1:nx,-1:ny), ne(-1:nx,-1:ny),
     &   qs(-1:nx,-1:ny,0:1)
*   ..output arguments (unspecified on entry)
      real (kind=R8) ::
     &   ehx(-1:nx,-1:ny,0:1)
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
      integer ncall, n2, ix, iy
      integer, save :: mode_ehx = 0                                      !srv 06.04.15
      integer, save :: pot_eq = 1, no_current= 0                         !srv 29.01.20
      real (kind=R8) ::
     *   t0, hy1_(-1:nx,-1:ny),
     *   dpop(-1:nx,-1:ny,0:1), wrk0(-1:nx,-1:ny,0:1), 
     *   wrk1(-1:nx,-1:ny,0:1), wrk2(-1:nx,-1:ny,0:1), wrk(-1:nx,-1:ny),
     *   dum(-1:nx,-1:ny)
      save ncall
*   ..procedures
      external xertst
      external b2xvsg
      real (kind=R8) ::
     *   hy1
      hy1(ix,iy) = hy(ix,iy)*qz(ix,iy,1)
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
*   ..set n2
      n2 = (nx+2)*(ny+2)
*   ..test nx, ny
      call xertst (0.le.nx.and.0.le.ny, 'faulty argument nx, ny')
*   ..extensive tests on first few calls
      if (ncall.lt.3) then
*    ..test sign of ne, te, csig
       call b2xvsg (n2, ne, 1, 'ne', '.gt.')
       call b2xvsg (n2, te, 1, 'te', '.gt.')
      endif

      if(pot_eq.ne.1 .and. no_current.ne.0) then                         !srv 29.01.20 {
        ehx = 0.0_R8
      else
        do iy=-1,ny
          do ix=-1,nx
            hy1_(ix,iy)=hy1(ix,iy)
          enddo
        enddo

cwdk  wrk will be used for interpolation to vertices; now: just 1.0_R8
c     to be considered: volume-weighted interpolation
        wrk = 1.0_R8

*   ..compute ehx
        if(mode_ehx.eq.0) then
*       ..compute poloidal differences of po and pe
          call diff_p(nx, ny, 0, hx, hy, wrk, qs, po, dum, dpop)
          call diff_p(nx, ny, 0, hx, hy, wrk, qs, ne*te, dum, wrk0)
*       ..interpolate ne to faces
          call intfaceh(nx,ny,leftix,leftiy,wrk,ne,wrk1(-1,-1,0))
          call intfacev(nx,ny,bottomix,bottomiy,wrk,ne,wrk1(-1,-1,1))
*       ..compute ehx      
          ehx = -dpop + wrk0/(qe*wrk1)
        else if (mode_ehx.eq.1) then
*       ..compute poloidal differences of po and pe
          call diff_p(nx, ny, 0, hx, hy, wrk, qs, po, dum, dpop)
          call diff_p(nx, ny, 0, hx, hy, wrk, qs, ne*te, dum, wrk0)
*       ..interpolate ne to faces
          call intfaceh(nx,ny,leftix,leftiy,hx,ne,wrk1(-1,-1,0))
          call intfacev(nx,ny,bottomix,bottomiy,hy1_,ne,wrk1(-1,-1,1))
*       ..compute ehx      
          ehx = -dpop + wrk0/(qe*wrk1)
        else
*       ..compute poloidal differences of po, te and log(ne) on faces
          call diff_p(nx, ny, 0, hx, hy, wrk, qs, po, dum, dpop)
          call diff_p(nx, ny, 0, hx, hy, wrk, qs, te, dum, wrk0)
          call diff_p(nx, ny, 0, hx, hy, wrk, qs, log(ne), dum, wrk1)
*       ..interpolate te to faces
          wrk = 1.0_R8
          call intfaceh(nx,ny,leftix,leftiy,wrk,te,wrk2(-1,-1,0))
          call intfacev(nx,ny,bottomix,bottomiy,wrk,te,wrk2(-1,-1,1))
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
