
*-----------------------------------------------------------------------
*.specification

      subroutine eseec0 (tef,tif,nef,pof,bpf,bbf,material,ix,iy,nx,ny,
     &  ns,ismain,sec,diagno)
      use b2mod_types
      use b2mod_constants
      use b2mod_indirect
      use b2mod_geo
      use b2mod_b2cmpa
      implicit none
*   ..input arguments (unchanged on exit)
      integer ix, iy, nx, ny, ns, ismain, material, diagno
      real (kind=R8) ::
     *  tef, tif, nef, pof, bpf, bbf
*   ..output argument (unspecified on entry)
      real (kind=R8) ::
     *  sec
*   ..common blocks
*-----------------------------------------------------------------------
*.documentation
*
*  1. purpose
*     
*     This subroutine is an interface between the b2 code and 
*     the -average- subroutine.  It receives data from b2 code and 
*     returns the average value over the emission angle
*     of the ratio mu = g/g0, where g0 is the secondary electron emission
*     coefficient for a given surface.
*     
*     Input parameters:
*     tef      - Electron temperature.
*     tif      - Ion temperature.
*     nef      - Electron density
*     pof      - Plasma-magnetic presheath potential drop
*     bpf      - Poloidal magnetic field
*     bbf      - Total magnetic field
*     material - Index of material for the first wall (see below)
*     ix       - poloidal mesh index.
*     iy       - radial mesh index.
*     nx       - poloidal mesh dimension - 2
*     ny       - radial mesh dimension - 2
*     ns       - number of species
*     ismain   - index of main plasma ion species
*     diagno   - output flag
*     
*     Output parameters:
*     sec - effective secondary electron emission coefficient.
*     
*     Internal variables
*     nmat - maximum number of considered materials - 1
*     imat - material index
*     wrkf - work function for different materials
*     eref - reference energy for different materials
*     ggmm - reference ESEEC for different materials
*     w0   - work function for selected material
*     u0   - reference energy for selected material
*     s0   - reference ESEEC for selected material
*     
*     Follows a table of the possible materials now considered.
*     
*     Number   |  Name
*     -----------------
*        0     |  Graphite
*     
*-----------------------------------------------------------------------
*.declarations

*   ..local variables
      integer nmat,imat,ncall
      parameter (nmat=9)
      real (kind=R8) :: wrkf(0:nmat), eref(0:nmat), ggmm(0:nmat)
      real (kind=R8) :: ef,w0,u0,s0,phi,g0,mu,debye,sphic
      real (kind=R8) ::
     *  r0,r1,r2,r3,z0,z1,z2,z3,nhatr,nhatz,bphatr,bphatz,bdotn
      save ncall
*   ..procedures
      intrinsic exp, log, sqrt, abs, asin, sin
      external subini, subend, xertst, average
*   ..initialisation
      data ncall/0/
*     
*     Assign table of work functions. Negative values for empty entries.
*     
      data     wrkf(0),      wrkf(1),      wrkf(2),      wrkf(3),  
     &         wrkf(4),      wrkf(5),      wrkf(6),      wrkf(7),  
     &         wrkf(8),      wrkf(9)       
     &
     &    /8.01e-19_R8,      -1.0_R8,      -1.0_R8,      -1.0_R8,
     &         -1.0_R8,      -1.0_R8,      -1.0_R8,      -1.0_R8,
     &         -1.0_R8,      -1.0_R8/  
*     
*     Assign table of reference energies. Negative values for empty entries.
*     
      data     eref(0),      eref(1),      eref(2),      eref(3),  
     &         eref(4),      eref(5),      eref(6),      eref(7),  
     &         eref(8),      eref(9)       
     &
     &    /4.81e-17_R8,      -1.0_R8,      -1.0_R8,      -1.0_R8,
     &         -1.0_R8,      -1.0_R8,      -1.0_R8,      -1.0_R8,
     &         -1.0_R8,      -1.0_R8/
*     
*     Assign table of reference ESEEC. Negative values for empty entries.
*     
      data     ggmm(0),      ggmm(1),      ggmm(2),      ggmm(3),  
     &         ggmm(4),      ggmm(5),      ggmm(6),      ggmm(7),  
     &         ggmm(8),      ggmm(9)       
     &
     &         /1.0_R8,      -1.0_R8,      -1.0_R8,      -1.0_R8,
     &         -1.0_R8,      -1.0_R8,      -1.0_R8,      -1.0_R8,
     &         -1.0_R8,      -1.0_R8/
*-----------------------------------------------------------------------
*.computation

* ..preliminaries
*   ..subprogram start-up calls
      call subini ('eseec0')
*   ..test input parameters
      call xertst (0.le.nx.and.0.le.ny, 'faulty argument nx, ny')
      call xertst (-1.le.ix.and.ix.le.nx, 'faulty argument ix')
      call xertst (-1.le.iy.and.iy.le.ny, 'faulty argument iy')
      call xertst (1.le.ns, 'faulty argument ns')
      call xertst (0.le.ismain.and.ismain.lt.ns,
     & 'faulty argument ismain')
      call xertst (.not.is_neutral(ismain),
     & 'invalid main plasma species ismain; must not be neutral')
      call xertst (0.le.material.and.material.le.nmat, 
     & 'faulty argument material')
      call xertst (0.lt.tef, 'faulty argument tef')
      call xertst (0.lt.tif, 'faulty argument tif')
      call xertst (0.lt.nef, 'faulty argument nef')
      call xertst (0.lt.bbf, 'faulty argument bbf')
*     
*     Test spatial position
*     
      call xertst (leftix(ix,iy).eq.-2.or.rightix(ix,iy).eq.nx+1,
     & 'eseec0 called for non-target plate point')
*     
*     Check material. Exit if not found.
*     
      imat = material
      call xertst (0.le.imat.and.imat.le.0,
     & 'Wall material index not supported')

      w0 = wrkf(imat)
      u0 = eref(imat)
      s0 = ggmm(imat)
*     
*     Compute angle between the magnetic field and the surface 
*     (i.e. the complementary to the incidence angle)
*      and the magnetic field
*     
      r0 = crx(ix,iy,0)
      z0 = cry(ix,iy,0)
      r1 = crx(ix,iy,1)
      z1 = cry(ix,iy,1)
      r2 = crx(ix,iy,2)
      z2 = cry(ix,iy,2)
      r3 = crx(ix,iy,3)
      z3 = cry(ix,iy,3)
      if (leftix(ix,iy).eq.-2) then
cxpb Western edge
        nhatr = (z3-z1)/sqrt((z3-z1)**2+(r1-r3)**2)
        nhatz = (r1-r3)/sqrt((z3-z1)**2+(r1-r3)**2)
        bphatr = (r1-r0+r3-r2)/sqrt((r1-r0+r3-r2)**2+(z1-z0+z3-z2)**2)
        bphatz = (z1-z0+z3-z2)/sqrt((r1-r0+r3-r2)**2+(z1-z0+z3-z2)**2)
      elseif (rightix(ix,iy).eq.nx+1) then
cxpb Eastern edge
        nhatr = (z0-z2)/sqrt((z0-z2)**2+(r2-r0)**2)
        nhatz = (r2-r0)/sqrt((z0-z2)**2+(r2-r0)**2)
        bphatr = (r0-r1+r2-r3)/sqrt((r0-r1+r2-r3)**2+(z0-z1+z2-z3)**2)
        bphatz = (z0-z1+z2-z3)/sqrt((r0-r1+r2-r3)**2+(z0-z1+z2-z3)**2)
      endif
      bdotn = (nhatr*bphatr+nhatz*bphatz)*bpf
      phi = (180.0_R8/pi)*abs(asin(bdotn/bbf))
*
*     Compute the electric field at the wall
*     The Debye length is computed according to NRL formulary (in MKS)
*
      debye = sqrt(eps0*tef/(nef*qe**2))
      sphic = sqrt(me*tef/mp/tif)
!     ef = (pof+tef/ev*log(sqrt(2*pi)*sin(phi*pi/180.0_R8)))/
!     & (3.0_R8*debye)
      ef = pof/(3.0_R8*debye)*
     &   (1.0_R8-(log(sin(phi*pi/180.0_R8)))/
     &    log(sphic))
*     
*     Compute ESEEC not averaged over distribution
*     
      g0 = s0*(2.0_R8*tef/u0)*
     *     exp(2.0_R8*(1.0_R8-sqrt(2.0_R8*tef/u0)))       !Tentative formula.
      g0 = 2.0_R8 * g0
*     
*     Compute electric field and averaged ratio gg/gg0
*     
      call average (ef,bbf,w0,phi,mu)
*     
*     Assign ESEEC
*      
      sec = mu * g0
*
*     Output
*
      if (diagno.ge.2) write (*,*) ix,iy,ef,phi,g0,mu
*
*   ..return
      ncall = ncall + 1
      call subend ()
      return

*-----------------------------------------------------------------------
*.end eseec0

      end subroutine eseec0
