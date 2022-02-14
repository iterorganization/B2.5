*-----------------------------------------------------------------------
*.specification

      subroutine b2wsav (nsav, nx, ny, ns, lblgm)
      use b2mod_types
      use b2mod_geo
      use b2mod_rates
      use b2mod_plasma
      use b2mod_indirect
      use b2mod_constants
      use b2mod_external
      implicit none
*   ..input arguments (unchanged on exit)
      integer nsav, nx, ny, ns
      character lblgm*120
*   ..output arguments (unspecified on entry)
*     (none)
*   ..common blocks
*-----------------------------------------------------------------------
*.documentation
*
*  1. purpose
*
*     B2WSAV writes the state of the plasma to an unformatted file,
*     in the style required for input into the 'old B2' code.
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
     *  pr(-1:nx,-1:ny), uac(-1:nx,-1:ny,0:ns-1)
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
      call b2xppr (nx, ny, ns, ne, na, te, ti, pr)
      do is = 0, ns-1
       call b2xpxc (nx, ny, ua(-1,-1,is), uac(-1,-1,is))
      enddo
*   ..produce output
      call b2wsax (nsav, nx, ny, ns, lblgm, rz2, pr, uac)

* ..return
      call subend ()
      return

*-----------------------------------------------------------------------
*.end b2wsav

      end subroutine b2wsav
