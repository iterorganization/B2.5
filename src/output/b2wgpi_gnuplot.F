*-----------------------------------------------------------------------
*.specification

      subroutine b2wgpi_gnuplot (nprt, nx, ny, ns, 
     &  ne, na, ua, te, ti, po, fna, fhe, fhi, fch,
     &  kinrgy, gs, qc, ni)
      use b2mod_types
      use b2mod_rates
      use b2mod_indirect
      use b2mod_constants
      use b2mod_external
      use b2mod_subsys
      implicit none
*   ..input arguments (unchanged on exit)
      integer nprt, nx, ny, ns
      real (kind=R8) ::
     *  na(-1:nx,-1:ny,0:ns-1), ua(-1:nx,-1:ny,0:ns-1),
     &  te(-1:nx,-1:ny), ti(-1:nx,-1:ny), po(-1:nx,-1:ny),
     &  fna(-1:nx,-1:ny,0:1,0:1,0:ns-1), fhe(-1:nx,-1:ny,0:1,0:1),
     &  fhi(-1:nx,-1:ny,0:1,0:1), fch(-1:nx,-1:ny,0:1,0:1),
     &  ne(-1:nx,-1:ny),ni(-1:nx,-1:ny,0:1),
     &  kinrgy(-1:nx,-1:ny,0:ns-1),
     &  gs(-1:nx,-1:ny, 0:2), qc(-1:nx,-1:ny,0:1)
*   ..common blocks
*-----------------------------------------------------------------------
*.documentation
*
*  1. purpose
*
*     B2WGPI_GNUPLOT provides graphical output of the plasma state on a 
*     uniform rectangular mesh. This routine is an interface to 
*     b2wgpj_gnuplot.
*
*
*  3. description (see also routine b2cdca)
*
*     This routine produces plots of the plasma state that results
*     from a B2 calculation. The routine operates as an interface
*     to b2wgpj_gnuplot.
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
     *  fne(-1:nx,-1:ny,0:1,0:1)
*   ..procedures
      external xertst
      external b2xpfe, b2wgpj_gnuplot 

*-----------------------------------------------------------------------
*.computation

* ..preliminaries
*   ..subprogram start-up calls
      call subini ('b2wgpi_gnuplot')
*   ..test nx, ny
      call xertst (0.le.nx.and.0.le.ny, 'faulty argument nx, ny')

* ..make plots
*   ..compute auxiliaries
#ifdef WG_TODO
      call b2xpfe (nx, ny, ns, ns_ext, qe, 
     &  rza, za_ext, fna, fa_ext, fch, fne)
#endif
*   ..plot state
      call b2wgpj_gnuplot (nprt, nx, ny, ns, na, ua, te, ti, po,
     &  fna, fhe, fhi, fch, ne, fne, kinrgy, gs, qc, ni)

* ..return
      call subend ()
      return

*-----------------------------------------------------------------------
*.end b2wgpi_gnuplot

      end subroutine b2wgpi_gnuplot
