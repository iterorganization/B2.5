*-----------------------------------------------------------------------
*.specification

      subroutine b2wsax (nsav, nx, ny, ns, lblgm, 
     &  crx, cry, bb, vol,
     &  hx, hy, qz, qc, qs, gs, pbs, cell_width, cell_height,
     &  na, ua, te, ti, po, fna, fhe, fhi, fch, ne, ne2, rz2, pr, uac)
      use b2mod_types
      use b2mod_indirect
      use b2mod_constants
      use b2mod_cellhelper
      use b2mod_b2cmpa
      implicit none
*   ..input arguments (unchanged on exit)
      integer nsav, nx, ny, ns
      character lblgm*120
      real (kind=R8) ::
     *  crx(-1:nx,-1:ny,0:3), cry(-1:nx,-1:ny,0:3),  
     &  bb(-1:nx,-1:ny,0:3), vol(-1:nx,-1:ny), hx(-1:nx,-1:ny),
     &  hy(-1:nx,-1:ny), qz(-1:nx,-1:ny,0:1), qc(-1:nx,-1:ny,0:1), 
     &  qs(-1:nx,-1:ny,0:1), gs(-1:nx,-1:ny,0:2), pbs(-1:nx,-1:ny,0:1),
     &  cell_width(-1:nx,-1:ny), cell_height(-1:nx,-1:ny)
      real (kind=R8) ::
     *  na(-1:nx,-1:ny,0:ns-1), ua(-1:nx,-1:ny,0:ns-1),
     &  te(-1:nx,-1:ny), ti(-1:nx,-1:ny), po(-1:nx,-1:ny),
     &  fna(-1:nx,-1:ny,0:1,0:1,0:ns-1), fhe(-1:nx,-1:ny,0:1,0:1),
     &  fhi(-1:nx,-1:ny,0:1,0:1), fch(-1:nx,-1:ny,0:1,0:1)
      real (kind=R8) ::
     *  ne(-1:nx,-1:ny), ne2(-1:nx,-1:ny), pr(-1:nx,-1:ny),
     &  rz2(-1:nx,-1:ny,0:ns-1), uac(-1:nx,-1:ny,0:ns-1)
*   ..output arguments (unspecified on entry)
*     (none)
*   ..common blocks
*-----------------------------------------------------------------------
*.documentation
*
*  1. purpose
*
*     B2WSAX writes the state of the plasma to an unformatted file,
*     in the style required for input into the 'old B2' code.
*
*
*  3. description (see also routine b2cdca)
*
*     This routine writes the state of the plasma to an unformatted
*     file, in the style required for input into the 'old B2' code.
*     One difference between the two codes is in the dimensioning
*     of the arrays: (-1:nx,-1:ny) now versus (0:nx+1,0:ny+1) in the
*     old code; and for quantities defined on cell faces, in the
*     new code it is the first column or the first row of the array
*     that must contain zeroes whereas in the old code it is the last
*     column or row.
*
*
*  5. parameters (see also routine b2cdcv)
*
*     nsav - integer, input.
*     nsav specifies a unit number for unformatted output.
*
*     nx, ny, ns - integer, input.
*     See the calling routine for description.
*
*     lblgm - character*120, input.
*     lblgm identifies the instance of the program employed to prepare
*     the geometry.
*
*     crx, cry, bb, vol, hx, hy, qz, qc, qs, gs, pbs - real*8 array, input.
*     na, ua, te, ti, po, fna, fhe, fhi, fch - real*8 array, input.
*     ne, ne2, pr, uac - real*8 array, input.
*     See the calling routine for description.
*
*
*  8. further comments
*
*!!   Beware: The new code works with heat fluxes whereas the old code
*     used energy fluxes. That conversion is not made in this output
*     routine.
*
*-----------------------------------------------------------------------
*.declarations

*   ..local variables
      integer nmax, n2, ix, iy, i, is
      parameter (nmax=300)
      real (kind=R8) ::
     *   t0, t1, zff(1:6), zxx(-1:nmax,1:6), zyy(-1:nmax,1:6),
     *   zeff(-1:nx,-1:ny)
*   ..procedures
      real (kind=R8) ::
     *   b2ssum, pit, hy1, uax, uay, uaz, crxc, cryc
      external subini, subend, xertst, b2ssum
      external b2xvsg, b2xvgm, b2xvps, b2xgxy
      pit(ix,iy) = bb(ix,iy,0)/bb(ix,iy,3)
      hy1(ix,iy) = hy(ix,iy)*qz(ix,iy,1)
      uax(ix,iy,is) = fna(ix,iy,0,0,is)*
     &  (hx(ix,iy)/na(leftix(ix,iy),leftiy(ix,iy),is)+
     1   hx(leftix(ix,iy),leftiy(ix,iy))/na(ix,iy,is))/
     &  (vol(leftix(ix,iy),leftiy(ix,iy))+vol(ix,iy))
      uay(ix,iy,is) = fna(ix,iy,1,1,is)*
     &  (hy1(ix,iy)/na(bottomix(ix,iy),bottomiy(ix,iy),is)+
     1   hy1(bottomix(ix,iy),bottomiy(ix,iy))/na(ix,iy,is))/
     &  (vol(bottomix(ix,iy),bottomiy(ix,iy))+vol(ix,iy))
      uaz(ix,iy,is) = uac(ix,iy,is)*abs(bb(ix,iy,2))/bb(ix,iy,3)
      crxc(ix,iy) = (crx(ix,iy,0)+crx(ix,iy,1)+crx(ix,iy,2)+
     &  crx(ix,iy,3))/4
      cryc(ix,iy) = (cry(ix,iy,0)+cry(ix,iy,1)+cry(ix,iy,2)+
     &  cry(ix,iy,3))/4

*-----------------------------------------------------------------------
*.computation

* ..preliminaries
*   ..subprogram start-up calls
      call subini ('b2wsax')
*   ..set n2
      n2 = (nx+2)*(ny+2)
*   ..test nsav
      call xertst (1.le.nsav, 'faulty argument nsav')
*   ..test nx, ny, ns
      call xertst (0.le.nx.and.0.le.ny, 'faulty argument nx, ny')
      call xertst (nx.le.nmax.and.ny.le.nmax, 'faulty parameter nmax')
      call xertst (1.le.ns, 'faulty argument ns')
*   ..test sign of bb
      call b2xvsg (n2, bb(-1,-1,3), 1, 'bb3', '.gt.')
*   ..test metric coefficients
      call b2xvgm (nx, ny, vol, hx, hy, qz, qc, qs, gs, pbs,
     &  cell_width, cell_height)
*   ..test plasma state
#ifdef WG_TODO
      call b2xvps (nx, ny, ns, 
     &  ne, na, ua, te, ti, po, fna, fhe, fhi, fch)
      call b2xzef (nx, ny, ns, rz2, na, ne, zeff)
#endif

* ..compute zxx and zyy
*   ..compute zxx(1:2) and zyy(1:2)
      call b2xgxy (nx, ny, crx, cry, 0, 0, zxx(-1,1), zxx(-1,2),
     &  zyy(-1,1), zyy(-1,2))
*   ..compute zxx(3:6)
      do ix = -1, nx
*    ..'vol(ix)'
       zxx(ix,3) = b2ssum(ny+2,vol(ix,-1),nx+2)
*    ..'hx(ix)'
       t0 = 0.0_R8
       do iy = -1, ny
        t0 = t0+vol(ix,iy)/hx(ix,iy)
       enddo
       zxx(ix,4) = zxx(ix,3)/t0
*    ..'sx(ix)'
       zxx(ix,5) = b2ssum(ny+2,gs(ix,-1,0),nx+2)
*    ..'pit(ix)'
       t0 = 0.0_R8
       do iy = -1, ny
        t0 = t0+vol(ix,iy)*pit(ix,iy)
       enddo
       zxx(ix,6) = t0/zxx(ix,3)
      enddo
*   ..compute zyy(3:6)
      do iy = -1, ny
*    ..'vol(iy)'
       zyy(iy,3) = b2ssum(nx+2,vol(-1,iy),1)
*    ..'hy(iy)'
       t0 = 0.0_R8
       do ix = -1, nx
        t0 = t0+vol(ix,iy)/hy(ix,iy)
       enddo
       zyy(iy,4) = zyy(iy,3)/t0
*    ..'sy(iy)'
       zyy(iy,5) = b2ssum(nx+2,gs(-1,iy,1),1)
*    ..'pit(iy)'
       t0 = 0.0_R8
       t1 = 0.0_R8
       do ix = -1, nx
        t0 = t0+hx(ix,iy)
        t1 = t1+hx(ix,iy)/pit(ix,iy)
       enddo
       zyy(iy,6) = t0/t1
      enddo
*   ..compute zff(1:6)
*   ..total volume
      zff(1) = b2ssum(n2,vol,1)
*   ..global hx
      zff(2) = b2ssum(nx+2,zxx(-1,4),1)
*   ..global hy
      zff(3) = b2ssum(ny+2,zyy(-1,4),1)
*   ..global sx
      zff(4) = zff(1)/zff(2)
*   ..global sy
      zff(5) = zff(1)/zff(3)
*   ..global pit
      t0 = 0.0_R8
      t1 = 0.0_R8
      do iy = -1, ny
       t0 = t0+zyy(iy,6)*zyy(iy,3)
       t1 = t1+zyy(iy,3)
      enddo
      zff(6) = t0/t1

* ..produce unformatted output
*   ..save label, nx, ny, ns
      write (nsav) lblgm(61:120), nx, ny, ns
*   ..save reduced geometry
      write (nsav) (zff(i),i=1,6),
     &  (zxx(ix,1),ix=-1,nx),
     &  (zxx(ix,2),ix=0,nx),0.0_R8,
     &  (zxx(ix,3),ix=-1,nx),
     &  (zxx(ix,4),ix=-1,nx),
     &  (zxx(ix,5),ix=0,nx),0.0_R8,
     &  (zxx(ix,6),ix=-1,nx),
     &  (zyy(iy,1),iy=-1,ny),
     &  (zyy(iy,2),iy=0,ny),0.0_R8,
     &  (zyy(iy,3),iy=-1,ny),
     &  (zyy(iy,4),iy=-1,ny),
     &  (zyy(iy,5),iy=0,ny),0.0_R8,
     &  (zyy(iy,6),iy=-1,ny)
*   ..save geometry
      write (nsav)
     &  ((vol(ix,iy),ix=-1,nx),iy=-1,ny),
     &  ((hx(ix,iy),ix=-1,nx),iy=-1,ny),
     &  ((hy(ix,iy),ix=-1,nx),iy=-1,ny),
     &  ((gs(ix,iy,0),ix=0,nx),0.0_R8,iy=-1,ny),
     &  ((gs(ix,iy,1),ix=-1,nx),iy=0,ny),(0.0_R8,ix=-1,nx),
     &  ((pit(ix,iy),ix=-1,nx),iy=-1,ny),
     &  ((crxc(ix,iy),ix=-1,nx),iy=-1,ny),
     & (((crx(ix,iy,i),ix=-1,nx),iy=-1,ny),i=0,3),
     &  ((cryc(ix,iy),ix=-1,nx),iy=-1,ny),
     & (((cry(ix,iy,i),ix=-1,nx),iy=-1,ny),i=0,3)
*   ..save mi, zamin, zamax
      write (nsav) (mp*am(is),zamin(is),zamax(is),is=0,ns-1)
*   ..save plasma state
      write (nsav)
     &  (((na(ix,iy,is),ix=-1,nx),iy=-1,ny),is=0,ns-1),
     &  (((uax(ix,iy,is),ix=0,nx),0.0_R8,iy=-1,ny),is=0,ns-1),
     &  (((uay(ix,iy,is),ix=-1,nx),iy=0,ny),(0.0_R8,ix=-1,nx),
     &    is=0,ns-1),
     &  (((uaz(ix,iy,is),ix=-1,nx),iy=-1,ny),is=0,ns-1),
     &  ((te(ix,iy),ix=-1,nx),iy=-1,ny),
     &  ((ti(ix,iy),ix=-1,nx),iy=-1,ny)
      write (nsav)
     &  ((pr(ix,iy),ix=-1,nx),iy=-1,ny),
     &  ((ne(ix,iy),ix=-1,nx),iy=-1,ny),
     &  ((po(ix,iy),ix=-1,nx),iy=-1,ny),
     &  ((zeff(ix,iy),ix=-1,nx),iy=-1,ny),
     &  (((ua(ix,iy,is),ix=-1,nx),iy=-1,ny),is=0,ns-1)
*   ..save fluxes
      write (nsav)
     &  (((fna(ix,iy,0,0,is),ix=0,nx),0.0_R8,iy=-1,ny),is=0,ns-1),
     &  (((fna(ix,iy,1,1,is),ix=-1,nx),iy=0,ny),(0.0_R8,ix=-1,nx),
     &    is=0,ns-1),
     &  ((fhe(ix,iy,0,0),ix=0,nx),0.0_R8,iy=-1,ny),
     &  ((fhe(ix,iy,1,1),ix=-1,nx),iy=0,ny),(0.0_R8,ix=-1,nx),
     &  ((fhi(ix,iy,0,0),ix=0,nx),0.0_R8,iy=-1,ny),
     &  ((fhi(ix,iy,1,1),ix=-1,nx),iy=0,ny),(0.0_R8,ix=-1,nx)
*!!   (heat fluxes should be energy fluxes.)
*   ..save zeroes for residuals
      write (nsav)
     &  (0.0_R8,i=0,ns-1),0.0_R8,0.0_R8
*   ..save zeroes for common blocks
      write (nsav)
     &  (0.0_R8,i=0,8*ns-1),
     &  (0.0_R8,i=0,10),
     &  (0.0_R8,i=0,2*ns-1)
      write (nsav)
     &  (0.0_R8,i=0,4*ns-1),
     &  0.0_R8,0.0_R8

* ..return
      call subend ()
      return

*-----------------------------------------------------------------------
*.end b2wsax

      end subroutine b2wsax
