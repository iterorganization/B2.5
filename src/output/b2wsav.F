*-----------------------------------------------------------------------
*.specification

      subroutine b2wsav (nsav, nx, ny, ns, lblgm, crx, cry, bb, vol,
     &  hx, hy, qz, qc, qcb, gs, pbs, cell_width, cell_height,
     &  na, ua, te, ti, po, fna, fhe, fhi, fch)
      use b2mod_types
      use b2mod_rates
      use b2mod_indirect
      use b2mod_constants
      use b2mod_external
      implicit none
*   ..input arguments (unchanged on exit)
      integer nsav, nx, ny, ns
      character lblgm*120
      real (kind=R8) ::
     *  crx(-1:nx,-1:ny,0:3), cry(-1:nx,-1:ny,0:3), qc(-1:nx,-1:ny),
     &  bb(-1:nx,-1:ny,0:3), vol(-1:nx,-1:ny), hx(-1:nx,-1:ny),
     &  hy(-1:nx,-1:ny), qz(-1:nx,-1:ny,0:1), gs(-1:nx,-1:ny,0:2),
     &  pbs(-1:nx,-1:ny,0:1), qcb(-1:nx,-1:ny),
     &  cell_width(-1:nx,-1:ny), cell_height(-1:nx,-1:ny)
      real (kind=R8) ::
     *  na(-1:nx,-1:ny,0:ns-1), ua(-1:nx,-1:ny,0:ns-1),
     &  te(-1:nx,-1:ny), ti(-1:nx,-1:ny), po(-1:nx,-1:ny),
     &  fna(-1:nx,-1:ny,0:1,0:ns-1), fhe(-1:nx,-1:ny,0:1),
     &  fhi(-1:nx,-1:ny,0:1), fch(-1:nx,-1:ny,0:1)
*   ..output arguments (unspecified on entry)
*     (none)
*   ..common blocks
*-----------------------------------------------------------------------
*.documentation
*
*  1. purpose
*
*     B2WSAV writes the state of the plasma to an unformatted file,
*     in the style required for input into the 'old b2' code.
*     The routine operates as an interface to b2wsax.
*
*
*  3. description (see also routine b2cdca)
*
*     See the routine b2wsax.
*
*
*  5. parameters (see also routine b2cdcv)
*
*     Refer to the calling routine for a description of all parameters.
*
*-----------------------------------------------------------------------
*.declarations

*   ..local variables
      integer is
      real (kind=R8) ::
     *  pr(-1:nx,-1:ny), uac(-1:nx,-1:ny,0:ns-1),
     1  ne(-1:nx,-1:ny), ne2(-1:nx,-1:ny)
*   ..procedures
      external subini, subend
      external b2xpne, b2xppr, b2xpxc, b2wsax

*-----------------------------------------------------------------------
*.computation

* ..preliminaries
*   ..subprogram start-up calls
      call subini ('b2wsav')

* ..produce output
*   ..compute auxiliaries
      call b2xpne (nx, ny, ns, rza, na, ne_ext, ne)
      call b2xpne (nx, ny, ns, rz2, na, ne2_ext, ne2)
      call b2xppr (nx, ny, ns, ne, na, te, ti, pr)
      do is = 0, ns-1
       call b2xpxc (nx, ny, ua(-1,-1,is), uac(-1,-1,is))
      enddo
*   ..produce output
      call b2wsax (nsav, nx, ny, ns, lblgm, crx, cry, bb, vol,
     &  hx, hy, qz, qc, qcb, gs, pbs, cell_width, cell_height,
     &  na, ua, te, ti, po, fna, fhe, fhi, fch, ne, ne2, rz2, pr, uac)

* ..return
      call subend ()
      return

*-----------------------------------------------------------------------
*.end b2wsav

      end subroutine b2wsav
