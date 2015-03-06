*-----------------------------------------------------------------------
*.specification

      subroutine b2wfpi (nprt, idout, nx, ny, ns, 
     &  ne, na, ua, te, ti, po, fna, fhe, fhi, fch)
      use b2mod_types
      use b2mod_rates
      use b2mod_indirect
      use b2mod_constants
      use b2mod_external
      use b2mod_b2cmpa
      implicit none
*   ..input arguments (unchanged on exit)
      integer nprt, nx, ny, ns
      character idout*(*)
      real (kind=R8) ::
     *  na(-1:nx,-1:ny,0:ns-1), ua(-1:nx,-1:ny,0:ns-1),
     &  te(-1:nx,-1:ny), ti(-1:nx,-1:ny), po(-1:nx,-1:ny),
     &  fna(-1:nx,-1:ny,0:1,0:ns-1), fhe(-1:nx,-1:ny,0:1),
     &  fhi(-1:nx,-1:ny,0:1), fch(-1:nx,-1:ny,0:1),
     &  ne(-1:nx,-1:ny)
*   ..common blocks
*-----------------------------------------------------------------------
*.documentation
*
*  1. purpose
*
*     B2WFPI writes a summary of the plasma state to a formatted file.
*     This routine is an interface to b2wfpj.
*
*
*  3. description (see also routine b2cdca)
*
*     This routine writes to a formatted file summary information
*     describing the plasma state. The routine operates as an interface
*     to b2wfpj.
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
     * fni(-1:nx,-1:ny,0:1), fne(-1:nx,-1:ny,0:1),
     1 wrk0(-1:nx,-1:ny),
     2 ni(-1:nx,-1:ny,0:1), 
     3 ro(-1:nx,-1:ny), pr(-1:nx,-1:ny)
*   ..procedures
      external subini, subend, xertst
      external b2xpni, b2xppr, b2xpro, b2xpfe, b2xpfi, b2wfpj

*-----------------------------------------------------------------------
*.computation

* ..preliminaries
*   ..subprogram start-up calls
      call subini ('b2wfpi')
*   ..test nprt, nx, ny
      call xertst (1.le.nprt, 'faulty argument nprt')
      call xertst (0.le.nx.and.0.le.ny, 'faulty argument nx, ny')
*     (further tests are performed inside b2wfpj)

* ..write state
*   ..compute auxiliaries
      call b2xpni (nx, ny, ns, na, ni)
      call b2xppr (nx, ny, ns, ne, na, te, ti, pr)
      call b2xpro (nx, ny, ns, mp, am, na, ro)
      call b2xpfe (nx, ny, ns, ns_ext, qe, 
     &  rza, za_ext, fna, fa_ext, fch, fne)
      call b2xpfi (nx, ny, ns, fna, fni)
*   ..write state
      call b2wfpj (nprt, idout, nx, ny, ns, na, ua, te, ti, po,
     &  fna, fhe, fhi, fch, ne, ni, pr, ro,
     &  fne, fni, wrk0)

* ..return
      call subend ()
      return

*-----------------------------------------------------------------------
*.end b2wfpi

      end subroutine b2wfpi
