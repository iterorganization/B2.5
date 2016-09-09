*-----------------------------------------------------------------------
*.specification

      subroutine b2xehy (nx, ny, hx, hy, qz, te, po, ne, ehy)           !srv 13.01.17
      use b2mod_types
      use b2mod_indirect
      use b2mod_constants
      implicit none
*   ..input arguments (unchanged on exit)
      integer nx, ny
      real (kind=R8) ::
     &   hx(-1:nx,-1:ny), hy(-1:nx,-1:ny), qz(-1:nx,-1:ny,0:1),
     &   te(-1:nx,-1:ny), po(-1:nx,-1:ny), ne(-1:nx,-1:ny)
*   ..output arguments (unspecified on entry)
      real (kind=R8) ::
     &   ehy(-1:nx,-1:ny,0:1)
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
      integer ncall, n2, ix, iy
      integer, save :: mode_ehy = 0                                     !srv 06.04.15
      real (kind=R8) ::
     *   t0, hy1_(-1:nx,-1:ny),
     *   dpor(-1:nx,-1:ny,0:1), wrk0(-1:nx,-1:nx,0:1), 
     *   wrk1(-1:nx,-1:ny,0:1), wrk2(-1:nx,-1:ny,0:1), wrk(-1:nx,-1:ny)
      save ncall
*   ..procedures
      external subini, subend, xertst
      external b2xvsg
      real (kind=R8) ::
     *   hy1
      hy1(ix,iy) = hy(ix,iy)*qz(ix,iy,1)                                !srv 07.05.08 04.07.14
*   ..initialization
      data ncall/0/
*-----------------------------------------------------------------------
*.computation

* ..preliminaries
*   ..subprogram start-up calls
      call subini ('b2xehy')
*   ..set internal parameters on first call                             !srv
      if (ncall.eq.0) then                                              !srv
       call ipgeti ('b2tfhe_mode_ehy',mode_ehy)                         !srv 06.04.15
      endif                                                             !srv
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

      do iy=-1,ny
        do ix=-1,nx
          hy1_(ix,iy)=hy1(ix,iy)
        enddo
      enddo

cwdk  wrk will be used for interpolation to vertices; now: just 1.0_R8
c     to be considered: volume-weighted interpolation
      wrk = 1.0_R8

*   ..compute ehy
      if(mode_ehy.eq.0) then
*     ..compute radial differences of po and pe
        call diff_r(nx, ny, hx, hy, wrk, qz, po, dpor)
        call diff_r(nx, ny, hx, hy, wrk, qz, ne*te, wrk0)
*     ..interpolate ne to faces
        call intfaceh(nx,ny,leftix,leftiy,wrk,ne,wrk1(-1,-1,0))
        call intfacev(nx,ny,bottomix,bottomiy,wrk,ne,wrk1(-1,-1,1))
*     ..compute ehy      
        ehy = -dpor + wrk0/(qe*wrk1)
      else if (mode_ehy.eq.1) then
*     ..compute radial differences of po and pe
        call diff_r(nx, ny, hx, hy, wrk, qz, po, dpor)
        call diff_r(nx, ny, hx, hy, wrk, qz, ne*te, wrk0)
*     ..interpolate ne to faces
cwdk TODO: implement this as a standard interpolation option...
        do iy = -1, ny
          do ix = -1, nx
            if(leftix(ix,iy).eq.-2) then
              wrk1(ix,iy,0) = 0.0_R8
            else
              wrk1(ix,iy,0) = 
     &            (ne(leftix(ix,iy),leftiy(ix,iy))*hx(ix,iy)+
     &             ne(ix,iy)*hx(leftix(ix,iy),leftiy(ix,iy)))/
     &            (hx(leftix(ix,iy),leftiy(ix,iy))+hx(ix,iy))
            endif
            if(bottomiy(ix,iy).eq.-2) then
              wrk1(ix,iy,1) = 0.0_R8
	    else
              wrk1(ix,iy,1) = 
     &            (ne(bottomix(ix,iy),bottomiy(ix,iy))*hy1_(ix,iy)+
     &             ne(ix,iy)*hy1_(bottomix(ix,iy),bottomiy(ix,iy)))/
     &            (hy1_(bottomix(ix,iy),bottomiy(ix,iy))+hy1_(ix,iy))
            endif
          enddo
        enddo
*     ..compute ehy      
        ehy = -dpor + wrk0/(qe*wrk1)
      else
*     ..compute radial differences of po, te and log(ne) on faces
        call diff_r(nx, ny, hx, hy, wrk, qz, po, dpor)
        call diff_r(nx, ny, hx, hy, wrk, qz, te, wrk0)
        call diff_r(nx, ny, hx, hy, wrk, qz, log(ne), wrk1)
*     ..interpolate te to faces
        wrk = 1.0_R8
        call intfaceh(nx,ny,leftix,leftiy,wrk,ne,wrk2(-1,-1,0))
        call intfacev(nx,ny,bottomix,bottomiy,wrk,ne,wrk2(-1,-1,1))
*     ..compute ehy      
        ehy = -dpor + (wrk0 + wrk1*wrk2)/qe
      endif

* ..return
      ncall = ncall+1
      call subend ()
      return

*-----------------------------------------------------------------------
*.end b2xehy

      end
