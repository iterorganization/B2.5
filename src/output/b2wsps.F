*-----------------------------------------------------------------------
*.specification

      subroutine b2wsps (nspc, nx, ny, ns,
     &  crx, cry, bb, vol,
     &  hx, hy, qz, qc, qs, gs, pbs, cell_width, cell_height,
     &  ne, na, ua, te, ti, po, fna, fhe, fhi, fch)
      use b2mod_types
      use b2mod_rates
      use b2mod_indirect
      use b2mod_constants
      use b2mod_external
      use b2mod_b2cmpa
      implicit none
*   ..input arguments (unchanged on exit)
      integer nspc, nx, ny, ns
      real (kind=R8) ::
     *  crx(-1:nx,-1:ny,0:3), cry(-1:nx,-1:ny,0:3),
     &  bb(-1:nx,-1:ny,0:3), vol(-1:nx,-1:ny), hx(-1:nx,-1:ny),
     &  hy(-1:nx,-1:ny), qz(-1:nx,-1:ny,0:1), gs(-1:nx,-1:ny,0:2),
     &  qc(-1:nx,-1:ny,0:1), qs(-1:nx,-1:ny,0:1), pbs(-1:nx,-1:ny,0:1),
     &  cell_width(-1:nx,-1:ny), cell_height(-1:nx,-1:ny)
      real (kind=R8) ::
     *  na(-1:nx,-1:ny,0:ns-1), ua(-1:nx,-1:ny,0:ns-1),
     &  te(-1:nx,-1:ny), ti(-1:nx,-1:ny), po(-1:nx,-1:ny),
     &  ne(-1:nx,-1:ny)
      real (kind=R8) ::
     *  fna(-1:nx,-1:ny,0:1,0:1,0:ns-1), fhe(-1:nx,-1:ny,0:1,0:1),
     &  fhi(-1:nx,-1:ny,0:1,0:1), fch(-1:nx,-1:ny,0:1,0:1)
*   ..common blocks
*-----------------------------------------------------------------------
*.documentation
*
*  1. purpose
*
*     B2WSPS provides special purpose printed output concerning the
*     plasma state. This routine is an interface to b2wsst.
*
*
*  3. description (see also routine b2cdca)
*
*     This routine is available for use as the occasion may require;
*     it is a suitable location for code to produce output that is of
*     interest only for a specific application of b2.
*     The routine is set up as an interface to b2wsst.
*
*
*  5. parameters (see also routine b2cdcv)
*
*     Refer to the calling routine for a description of all parameters.
*
*-----------------------------------------------------------------------
*.declarations

*   ..local variables
      real (kind=R8) ::
     * fni(-1:nx,-1:ny,0:1,0:1), fne(-1:nx,-1:ny,0:1,0:1),
     1 wrk1(-1:nx,-1:ny), wrk0(-1:nx,-1:ny),
     2 ni(-1:nx,-1:ny,0:1), 
     3 ro(-1:nx,-1:ny), pr(-1:nx,-1:ny)
*   ..procedures
      external subini, subend, xertst
      external b2xpni, b2xppr, b2xpro, b2xpfe, b2xpfi, b2wsst

*-----------------------------------------------------------------------
*.computation

* ..preliminaries
*   ..subprogram start-up calls
      call subini ('b2wsps')
*   ..test nspc, nx, ny, ns
      call xertst (1.le.nspc, 'faulty argument nspc')
      call xertst (0.le.nx.and.0.le.ny, 'faulty argument nx, ny')
      call xertst (1.le.ns, 'faulty argument ns')
*     (further tests are performed inside b2wsst)

* ..special output
*   ..compute auxiliaries
      call b2xpni (nx, ny, ns, na, ni)
      call b2xppr (nx, ny, ns, ne, na, te, ti, pr)
      call b2xpro (nx, ny, ns, mp, am, na, ro)
      call b2xpfe (nx, ny, ns, ns_ext, qe, 
     &  rza, za_ext, fna, fa_ext, fch, fne)
      call b2xpfi (nx, ny, ns, fna, fni)
*   ..call output routine
      call b2wsst (nspc, nx, ny, ns,
     &  crx, cry, bb, vol,
     &  hx, hy, qz, qc, qs, gs, pbs, cell_width, cell_height,
     &  na, ua, te, ti, po, fna, fhe, fhi, fch, ne, ni,
     &  pr, ro, fne, fni, wrk0, wrk1)

* ..return
      call subend ()
      return

*-----------------------------------------------------------------------
*.end b2wsps

      end subroutine b2wsps
