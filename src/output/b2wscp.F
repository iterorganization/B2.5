*-----------------------------------------------------------------------
*.specification

      subroutine b2wscp (nspc, nx, ny, ns)
      use b2mod_types
      use b2mod_constants
      use b2mod_b2cmgs
      use b2mod_b2cmpa
      use b2mod_b2cmpb
      use b2mod_b2cmpt
      use b2mod_subsys
      implicit none
*   ..input arguments (unchanged on exit)
      integer nspc, nx, ny, ns
*   ..output arguments (unspecified on entry)
*     (none)
*   ..common blocks
*-----------------------------------------------------------------------
*.documentation
*
*  1. purpose
*
*     B2WSCP provides special purpose printed output of data in the
*     physics common blocks.
*
*
*  3. description (see also routine b2cdca)
*
*     This routine is available for use as the occasion may require;
*     it is a suitable location for code to produce output that is of
*     interest only for a specific application of b2.
*     The basic version of b2wscp produces no output.
*
*
*  5. parameters (see also routine b2cdcv)
*
*     nspc - integer, input.
*     nspc specifies the unit number for special purpose
*     output.
*
*     nx, ny, ns - integer, input.
*     See the calling routine for description.
*
*     (parameters in common:)
*
*     /b2cmgs/ - input.
*     /b2cmgs/ specifies the singular regions in the geometry.
*     See the calling routine for description.
*
*     /b2cmpa/ - input.
*     /b2cmpa/ specifies basic physics parameters.
*     See the calling routine for description.
*
*     /b2cmpb/ - input.
*     /b2cmpb/ specifies physics parameters for boundary conditions.
*     See the calling routine for description.
*
*     /b2cmpt/ - input.
*     /b2cmpt/ specifies physics parameters for transport.
*     See the calling routine for description.
*
*-----------------------------------------------------------------------
*.declarations

*   ..local variables
*     (none)
*   ..procedures
      external xertst

*-----------------------------------------------------------------------
*.computation

* ..preliminaries
*   ..subprogram start-up calls
      call subini ('b2wscp')
*   ..test nspc
      call xertst (1.le.nspc, 'faulty argument nspc')
*   ..test nx, ny, ns
      call xertst (0.le.nx.and.0.le.ny, 'faulty argument nx, ny')
      call xertst (1.le.ns, 'faulty argument ns')

* ..produce output
*     (nothing)

* ..return
      call subend ()
      return

*-----------------------------------------------------------------------
*.formats

*-----------------------------------------------------------------------
*.end b2wscp

      end
