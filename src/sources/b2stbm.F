*-----------------------------------------------------------------------
*.specification

      subroutine b2stbm (nx, ny, ns, vol, hx, hy, qz, pbs, na, ua,
     &  te, ti, po, ne, ni, fna, fhe, fhi, fch, fne, fni,
     &  sna0, smo0, she0, shi0, sch0, wrk0, wrk1)
      use b2mod_types
      use b2mod_constants
      use b2mod_b2cmpa
      use b2mod_b2cmpb
      implicit none
*   ..input arguments (unchanged on exit)
      integer nx, ny, ns
      real (kind=R8) ::
     *   vol(-1:nx,-1:ny), hx(-1:nx,-1:ny), hy(-1:nx,-1:ny),
     &  qz(-1:nx,-1:ny,0:1), pbs(-1:nx,-1:ny,0:1)
      real (kind=R8) ::
     *   na(-1:nx,-1:ny,0:ns-1), ua(-1:nx,-1:ny,0:ns-1),
     &  te(-1:nx,-1:ny), ti(-1:nx,-1:ny), po(-1:nx,-1:ny),
     &  ne(-1:nx,-1:ny), ni(-1:nx,-1:ny,0:1)
      real (kind=R8) ::
     *   fna(-1:nx,-1:ny,0:1,0:ns-1), fhe(-1:nx,-1:ny,0:1),
     &  fhi(-1:nx,-1:ny,0:1), fch(-1:nx,-1:ny,0:1),
     &  fne(-1:nx,-1:ny,0:1), fni(-1:nx,-1:ny,0:1)
*   ..output arguments (unspecified on entry)
      real (kind=R8) ::
     *   sna0(-1:nx,-1:ny,0:1,0:ns-1), smo0(-1:nx,-1:ny,0:3,0:ns-1),
     &  she0(-1:nx,-1:ny,0:3), shi0(-1:nx,-1:ny,0:3),
     &  sch0(-1:nx,-1:ny,0:3)
*   ..workspace arguments (unspecified on entry and on exit)
      real (kind=R8) ::
     *   wrk0(-1:nx,-1:ny), wrk1(-1:nx,-1:ny)
*   ..common blocks
*-----------------------------------------------------------------------
*.documentation
*
*  1. purpose
*
*     B2STBM provides linearised expressions for the sources of
*     density, momentum, and heat, due to miscellaneous application
*     specific processes.
*
*
*-----------------------------------------------------------------------
*.declarations

*   ..local variables
      integer n2, ncall, is
*   ..procedures
      real (kind=R8) ::
     * smin, smax
      external subini, subend, xertst, sfill,
     * smin, smax, b2xvsg, b2xvff
*   ..initialisation
      data ncall/0/

*-----------------------------------------------------------------------
*.documentation-internal
*
*     This version of b2stbm just computes zero sources.
*
*-----------------------------------------------------------------------
*.computation

* ..preliminaries
*   ..subprogram start-up calls
      call subini ('b2stbm')
*   ..set n2
      n2 = (nx+2)*(ny+2)
*   ..test nx, ny, ns
      call xertst (0.le.nx.and.0.le.ny, 'faulty argument nx, ny')
      call xertst (1.le.ns, 'faulty argument ns')
*   ..extensive tests on first few calls
      if(ncall.lt.3) then
*    ..test sign of vol, hx, hy
       call b2xvsg (n2, vol, 1, 'vol', '.gt.')
       call b2xvsg (n2, hx, 1, 'hx', '.gt.')
       call b2xvsg (n2, hy, 1, 'hy', '.gt.')
*    ..test range of qz
       call xertst (-1.lt.smin(n2,qz(-1,-1,0),1).and.
     &   smax(n2,qz(-1,-1,0),1).lt.1, 'faulty argument range qz(,,0)')
       call xertst (0.lt.smin(n2,qz(-1,-1,1),1).and.
     &   smax(n2,qz(-1,-1,1),1).le.1, 'faulty argument range qz(,,1)')
*    ..test sign of na, ni, ne, te, ti
       call b2xvsg (n2*ns, na, 1, 'na', '.gt.')
       call b2xvsg (n2*2, ni, 1, 'ni', '.gt.')
       call b2xvsg (n2, ne, 1, 'ne', '.gt.')
       call b2xvsg (n2, te, 1, 'te', '.gt.')
       call b2xvsg (n2, ti, 1, 'ti', '.gt.')
*   ..test edge values of fluxes
       do is = 0, ns-1
        call b2xvff (nx, ny, fna(-1,-1,0,is), 'fna')
       enddo
       call b2xvff (nx, ny, fhe, 'fhe')
       call b2xvff (nx, ny, fhi, 'fhi')
       call b2xvff (nx, ny, fch, 'fch')
       call b2xvff (nx, ny, fne, 'fne')
       call b2xvff (nx, ny, fni, 'fni')
      endif

*   ..initialise sources to 0
      call sfill (n2*2*ns, 0.0_R8, sna0, 1)
      call sfill (n2*4*ns, 0.0_R8, smo0, 1)
      call sfill (n2*4, 0.0_R8, she0, 1)
      call sfill (n2*4, 0.0_R8, shi0, 1)
      call sfill (n2*4, 0.0_R8, sch0, 1)

* ..return
      ncall = ncall+1
      call subend ()
      return

*-----------------------------------------------------------------------
*.end b2stbm

      end subroutine b2stbm
