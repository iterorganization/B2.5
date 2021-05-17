*-----------------------------------------------------------------------
*.specification

      subroutine b2xvps (nx, ny, ns, ne, na, ua, te, ti, po,
     &  fna, fhe, fhi, fch)
      use b2mod_types
      implicit none
*   ..input arguments (unchanged on exit)
      integer nx, ny, ns
      real (kind=R8) ::
     *  na(-1:nx,-1:ny,0:ns-1), ua(-1:nx,-1:ny,0:ns-1),
     &  te(-1:nx,-1:ny), ti(-1:nx,-1:ny), po(-1:nx,-1:ny),
     &  ne(-1:nx,-1:ny)
      real (kind=R8) ::
     *  fna(-1:nx,-1:ny,0:1,0:ns-1), fhe(-1:nx,-1:ny,0:1),
     &  fhi(-1:nx,-1:ny,0:1), fch(-1:nx,-1:ny,0:1)
*   ..output arguments (unspecified on entry)
*     (none)

*-----------------------------------------------------------------------
*.documentation
*
*  1. purpose
*
*     B2XVPS tests the validity of the state of the plasma.
*
*
*  3. description (see also routine b2cdca)
*
*     This routine tests the validity of the state of the plasma:
*     densities and temperatures must be positive, and appropriate
*     edge values of the fluxes must vanish.
*
*
*  5. parameters (see also routine b2cdcv)
*
*     nx, ny, ns - integer, input.
*     ne, na, ua, te, ti, po - real array, input.
*     fna, fhe, fhi, fch - real array, input.
*     See the calling routine for description.
*
*-----------------------------------------------------------------------
*.declarations

*   ..local variables
      integer n2, is
*   ..procedures
      external subini, subend, xertst
      external b2xvsg, b2xvff

*-----------------------------------------------------------------------
*.computation

* ..preliminaries
*   ..subprogram start-up calls
      call subini ('b2xvps')
*   ..set n2
      n2 = (nx+2)*(ny+2)
*   ..test nx, ny, ns
      call xertst (0.le.nx.and.0.le.ny, 'faulty argument nx, ny')
      call xertst (1.le.ns, 'faulty argument ns')

* ..test state
*   ..test sign of ne, na, te, ti
      call b2xvsg (n2, ne, 1, 'ne', '.gt.')
      call b2xvsg (n2*ns, na, 1, 'na', '.gt.')
      call b2xvsg (n2, te, 1, 'te', '.gt.')
      call b2xvsg (n2, ti, 1, 'ti', '.gt.')

*   ..test edge values of fna, fhe, fhi, fch
      do is = 0, ns-1
       call b2xvff (nx, ny, fna(-1,-1,0,is), 'fna')
      enddo
      call b2xvff (nx, ny, fhe, 'fhe')
      call b2xvff (nx, ny, fhi, 'fhi')
      call b2xvff (nx, ny, fch, 'fch')

* ..return
      call subend ()
      return

*-----------------------------------------------------------------------
*.end b2xvps

      end subroutine b2xvps
