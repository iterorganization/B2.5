*-----------------------------------------------------------------------
*.specification

      subroutine b2mwit (nput, nx, ny, ns, itim, itcnt, resco, resmo,
     &  resmt, reshe, reshi, respo, corpa, corua, corut, corte, corti,
     &  corpo)
      use b2mod_types
      implicit none
*   ..input arguments (unchanged on exit)
      integer nput, nx, ny, ns, itim, itcnt
      real (kind=R8) ::
     *   resco(-1:nx,-1:ny,0:ns-1), resmo(-1:nx,-1:ny,0:ns-1),
     &  resmt(-1:nx,-1:ny), reshe(-1:nx,-1:ny), reshi(-1:nx,-1:ny),
     &  respo(-1:nx,-1:ny)
      real (kind=R8) ::
     *   corpa(-1:nx,-1:ny,0:ns-1), corua(-1:nx,-1:ny,0:ns-1),
     &  corut(-1:nx,-1:ny), corte(-1:nx,-1:ny), corti(-1:nx,-1:ny),
     &  corpo(-1:nx,-1:ny)
*   ..output arguments (unspecified on entry)
*     (none)

*-----------------------------------------------------------------------
*.documentation
*
*  1. purpose
*
*     B2MWIT writes to an un*formatted file detailed monitoring output
*     for the inner iterations of the b2 solution process.
*
*
*  3. description (see also routine b2cdca)
*
*     This routine writes to unit (nput) data describing the numerical
*     convergence. Unit nput must be connected to an un*formatted file.
*
*
*  5. parameters (see also routine b2cdcv)
*
*     nput - integer, input.
*     nput specifies a unit number for un*formatted output.
*
*     nx, ny, ns - integer, input.
*     See the calling routine for description.
*
*     itim - integer, input.
*     itim specifies the timestep counter. It is initialised to 0 and
*     incremented by 1 following each successful call to b2mndt.
*
*     itcnt - integer, input.
*     itcnt counts the number of inner iterations. It is initialised to
*     0 and incremented by 1 following each call to b2news.
*
*     resco - (-1:nx,-1:ny,0:ns-1) real*8 array, input.
*     For (ix,iy,is) in (-1:nx,-1:ny,0:ns-1), resco(ix,iy,is) specifies
*     the residual of the continuity equation on the (ix,iy) cell for
*     species (is).
*
*     For (ix,iy,is) in (-1:nx,-1:ny,0:ns-1), resmo(ix,iy,is) specifies
*     the residual of the momentum equation on the (ix,iy) cell for
*     species (is).
*
*
*     resmt - (-1:nx,-1:ny) real*8 array.
*     For (ix,iy) in (-1:nx,-1:ny), resmt(ix,iy) specifies the residual
*     of the total momentum equation on the (ix,iy) cell.
*
*
*     reshe - (-1:nx,-1:ny) real*8 array, input.
*     For (ix,iy) in (-1:nx,-1:ny), reshe(ix,iy) specifies the residual
*     of the electron heat balance equation on the (ix,iy) cell.
*
*     reshi - (-1:nx,-1:ny) real*8 array, input.
*     For (ix,iy) in (-1:nx,-1:ny), reshi(ix,iy) specifies the residual
*     of the all atom heat balance equation on the (ix,iy) cell.
*
*     respo - (-1:nx,-1:ny) real*8 array, input.
*     For (ix,iy) in (-1:nx,-1:ny), respo(ix,iy) specifies the
*     residual of the electric potential equation on the (ix,iy) cell.
*
*     corpa - (-1:nx,-1:ny,0:ns-1) real*8 array, input.
*     For (ix,iy,is) in (-1:nx,-1:ny,0:ns-1), corpa(ix,iy,is) specifies
*     the last computed correction to the pressure of species (is) on
*     the (ix,iy) cell. corpa is computed by solving an equation that
*     has resco as its right hand side. It is not necessarily the same
*     as the pressure correction over one complete iteration.
*
*
*     corua - (-1:nx,-1:ny,0:ns-1) real*8 array.
*     For (ix,iy,is) in (0:nx,-1:ny,0:ns-1), corua(ix,iy,is) specifies
*     the last computed correction to the parallel velocity ua of
*     species (is) on the (ix,iy) cell.
*     corua is computed by solving an equation that has resmo as its
*     right hand side. It is not necessarily the same as the change in
*     the parallel velocity over one complete iteration.
*
*
*
*     corut - (-1:nx,-1:ny) real*8 array.
*     For (ix,iy) in (0:nx,-1:ny), corut(ix,iy) specifies the last
*     computed correction to the parallel velocity ua of all species
*     simultaneously on the (ix,iy) cell.
*     corut is computed by solving an equation that has resmt as its
*     right hand side. It is not necessarily the same as the change in
*     the parallel velocity over one complete iteration.
*
*
*     corte - (-1:nx,-1:ny) real*8 array, input.
*     For (ix,iy) in (-1:nx,-1:ny), corte(ix,iy) specifies the last
*     computed correction to the electron temperature on the (ix,iy)
*     cell. corte is computed by solving an equation that has reshe as
*     its right hand side. It is not necessarily the same as the change
*     in the electron temperature over one complete iteration.
*
*     corti - (-1:nx,-1:ny) real*8 array, input.
*     For (ix,iy) in (-1:nx,-1:ny), corti(ix,iy) specifies the last
*     computed correction to the atom temperature on the (ix,iy) cell.
*     corti is computed by solving an equation that has reshi as its
*     right hand side. It is not necessarily the same as the change in
*     the atom temperature over one complete iteration.
*
*     corpo - (-1:nx,-1:ny) real*8 array, input.
*     For (ix,iy) in (-1:nx,-1:ny), corpo(ix,iy) specifies the last
*     computed correction to the electric potential on the (ix,iy)
*     cell. corpo is computed by solving an equation that has respo as
*     its right hand side. It is not necessarily the same as the change
*     in the electric potential over one complete iteration.
*
*
*  8. further comments
*
*     See routine b2cdca for explanation of 'un*formatted file'.
*
*-----------------------------------------------------------------------
*.declarations

*   ..local variables
      integer ncall, n2, idum(0:9)
      save ncall
*   ..procedures
      external subini, subend, xertst, cfwuin, cfwure
      external b2xvfx
*   ..initialisation
      data ncall/0/

*-----------------------------------------------------------------------
*.computation

* ..preliminaries
*   ..subprogram start-up calls
      call subini ('b2mwit')
*   ..set n2
      n2 = (nx+2)*(ny+2)
*   ..test nput
      call xertst (1.le.nput, 'faulty argument nput')
*   ..test nx, ny, ns
      call xertst (0.le.nx.and.0.le.ny, 'faulty argument nx, ny')
      call xertst (1.le.ns, 'faulty argument ns')
*   ..test itim, itcnt
      call xertst (0.le.itim.and.0.le.itcnt,
     &  'faulty argument itim, itcnt')

* ..save detailed convergence data
*   ..save itim, itcnt
      idum(0) = itim
      idum(1) = itcnt
      call cfwuin (nput, 2, idum, 'itim,itcnt')
*   ..save residuals and corrections
      call cfwure (nput, n2*ns, resco, 'resco')
      call cfwure (nput, n2*ns, corpa, 'corpa')
      call cfwure (nput, n2*ns, resmo, 'resmo')
      call cfwure (nput, n2*ns, corua, 'corua')
      call cfwure (nput, n2, resmt, 'resmt')
      call cfwure (nput, n2, corut, 'corut')
      call cfwure (nput, n2, reshe, 'reshe')
      call cfwure (nput, n2, corte, 'corte')
      call cfwure (nput, n2, reshi, 'reshi')
      call cfwure (nput, n2, corti, 'corti')
      call cfwure (nput, n2, respo, 'respo')
      call cfwure (nput, n2, corpo, 'corpo')

* ..return
      ncall = ncall+1
      call subend ()
      return

*-----------------------------------------------------------------------
*.end b2mwit

      end subroutine b2mwit
