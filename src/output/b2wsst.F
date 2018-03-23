*-----------------------------------------------------------------------
*.specification

      subroutine b2wsst (nspc, nx, ny, ns, crx, cry, bb, vol, hx, hy,
     &  qz, qc, qcb, gs, pbs, cell_width, cell_height,
     &  na, ua, te, ti, po, fna, fhe, fhi, fch,
     &  ne, ni, pr, ro, fne, fni, wrk0, wrk1)
      use b2mod_types
      use b2mod_constants
      use b2mod_b2cmpa
      implicit none
*   ..input arguments (unchanged on exit)
      integer nspc, nx, ny, ns
      real (kind=R8) ::
     *  crx(-1:nx,-1:ny,0:3), cry(-1:nx,-1:ny,0:3),
     &  bb(-1:nx,-1:ny,0:3), vol(-1:nx,-1:ny), hx(-1:nx,-1:ny),
     &  hy(-1:nx,-1:ny), qz(-1:nx,-1:ny,0:1), qc(-1:nx,-1:ny),
     &  qcb(-1:nx,-1:ny), gs(-1:nx,-1:ny,0:2), pbs(-1:nx,-1:ny,0:1),
     &  cell_width(-1:nx,-1:ny), cell_height(-1:nx,-1:ny)
      real (kind=R8) ::
     *  na(-1:nx,-1:ny,0:ns-1), ua(-1:nx,-1:ny,0:ns-1),
     &  te(-1:nx,-1:ny), ti(-1:nx,-1:ny), po(-1:nx,-1:ny),
     &  fna(-1:nx,-1:ny,0:1,0:1,0:ns-1), fhe(-1:nx,-1:ny,0:1,0:1),
     &  fhi(-1:nx,-1:ny,0:1,0:1), fch(-1:nx,-1:ny,0:1,0:1)
      real (kind=R8) ::
     *  ne(-1:nx,-1:ny), ni(-1:nx,-1:ny,0:1), pr(-1:nx,-1:ny),
     &  ro(-1:nx,-1:ny), fne(-1:nx,-1:ny,0:1,0:1), 
     &  fni(-1:nx,-1:ny,0:1,0:1)
*   ..workspace arguments (unspecified on entry and on exit)
      real (kind=R8) ::
     *  wrk0(-1:nx,-1:ny), wrk1(-1:nx,-1:ny)
*   ..common blocks
*-----------------------------------------------------------------------
*.documentation
*
*  1. purpose
*
*     B2WSST provides special purpose printed output concerning the
*     plasma state.
*
*
*  3. description (see also routine b2cdca)
*
*     This routine is available for use as the occasion may require;
*     it is a suitable location for code to produce output that is of
*     interest only for a specific application of b2.
*     The basic version of b2wsst produces no output.
*
*
*  5. parameters (see also routine b2cdcv)
*
*     nspc - integer, input.
*     nspc specifies the unit number for special purpose output.
*
*     nx, ny, ns - integer, input.
*     See the calling routine for description.
*
*     crx, cry, bb, vol, hx, hy, qz, qc, qcb, gs, pbs - real*8 array, input.
*     See the calling routine for description.
*
*     na, ua, te, ti, po, fna, fhe, fhi, fch - real*8 array, input.
*     See the calling routine for description.
*
*     ne, ni, pr, ro, fne, fni - (-1:nx,-1:ny) real*8 array, input.
*     See the calling routine for description.
*
*     wrk0, wrk1 - (-1:nx,-1:ny) real*8 array, workspace.
*
*     (parameters in common:)
*
*     /b2cmpa/ - input.
*     /b2cmpa/ specifies basic physics parameters.
*     See the calling routine for description.
*
*-----------------------------------------------------------------------
*.declarations

*   ..local variables
      integer n2, ix, iy, is
*   ..procedures
      real (kind=R8) ::
     *  b2ssum, rfnax, rfnay, rfhex, rfhey,
     &  rfhix, rfhiy, rfchx, rfchy, rfnex, rfney, rfnix, rfniy
      external subini, subend, xertst, b2ssum
      external b2xvsg, b2xvff, b2xvgm, b2xvps
      rfnax(ix,is) = b2ssum(ny+2,fna(ix,-1,0,0,is),nx+2)
      rfnay(iy,is) = b2ssum(nx+2,fna(-1,iy,1,1,is),1)
      rfhex(ix) = b2ssum(ny+2,fhe(ix,-1,0,0),nx+2)
      rfhey(iy) = b2ssum(nx+2,fhe(-1,iy,1,1),1)
      rfhix(ix) = b2ssum(ny+2,fhi(ix,-1,0,0),nx+2)
      rfhiy(iy) = b2ssum(nx+2,fhi(-1,iy,1,1),1)
      rfchx(ix) = b2ssum(ny+2,fch(ix,-1,0,0),nx+2)
      rfchy(iy) = b2ssum(nx+2,fch(-1,iy,1,1),1)
      rfnex(ix) = b2ssum(ny+2,fne(ix,-1,0,0),nx+2)
      rfney(iy) = b2ssum(nx+2,fne(-1,iy,1,1),1)
      rfnix(ix) = b2ssum(ny+2,fni(ix,-1,0,0),nx+2)
      rfniy(iy) = b2ssum(nx+2,fni(-1,iy,1,1),1)

*-----------------------------------------------------------------------
*.computation

* ..preliminaries
*   ..subprogram start-up calls
      call subini ('b2wsst')
*   ..set n2
      n2 = (nx+2)*(ny+2)
*   ..test nspc, nx, ny, ns
      call xertst (1.le.nspc, 'faulty argument nspc')
      call xertst (0.le.nx.and.0.le.ny, 'faulty argument nx, ny')
      call xertst (1.le.ns, 'faulty argument ns')
*   ..test sign of bb
      call b2xvsg (n2, bb(-1,-1,3), 1, 'bb3', '.gt.')
*   ..test metric coefficients
      call b2xvgm (nx, ny, vol, hx, hy, qz, qc, qcb, gs, pbs,
     &  cell_width, cell_height)
*   ..test plasma state
      call b2xvps (nx, ny, ns, 
     &  ne, na, ua, te, ti, po, fna, fhe, fhi, fch)
*   ..test sign of ni, pr, ro
      call b2xvsg (n2*2, ni, 1, 'ni', '.gt.')
      call b2xvsg (n2, pr, 1, 'pr', '.gt.')
      call b2xvsg (n2, ro, 1, 'ro', '.gt.')
*   ..test edge values of fne, fni
      call b2xvff (nx, ny, fne, 'fne')
      call b2xvff (nx, ny, fni, 'fni')

* ..produce output
*     (nothing)

* ..return
      call subend ()
      return

*-----------------------------------------------------------------------
*.formats

*-----------------------------------------------------------------------
*.end b2wsst

      end subroutine b2wsst
