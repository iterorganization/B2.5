*-----------------------------------------------------------------------
*.specification

      subroutine b2xtps (nx, ny, ns, nx1, ny1, ns1, nncutmax, 
     &  non_commensurate, nnreg, nncut, periodic_bc, jsep, jsep1, 
     &  region, region1,
     &  leftcut, leftcut1, rightcut, rightcut1,
     &  bottomcut, bottomcut1, topcut, topcut1,
     &  hx, hy, qc, hx1, hy1, qc1, 
     &  bundle_frac, bottomiy1, ne, na, ua, uadia, te, ti,	
     &  po, fna, fhe, fhi, fch, fch_32, fch_52, fch_p, kinrgy, 
     &  ne1, na1, ua1, uadia1, te1, ti1, po1, fna1, fhe1, fhi1, fch1, 
     &  fch1_32, fch1_52, fch1_p, kinrgy1)
      use b2mod_types
      implicit none
*   ..input arguments (unchanged on exit)
      integer nx, ny, ns, nx1, ny1, ns1, nncutmax, non_commensurate
      integer nnreg(0:2), nncut, periodic_bc, jsep, jsep1
      integer region(-1:nx,-1:ny,0:2), region1(-1:nx1,-1:ny1,0:2)
      integer leftcut(nncutmax), leftcut1(nncutmax)
      integer rightcut(nncutmax), rightcut1(nncutmax)
      integer bottomcut(nncutmax), bottomcut1(nncutmax)
      integer topcut(nncutmax), topcut1(nncutmax)
      real (kind=R8) ::
     *  hx(-1:nx,-1:ny), hy(-1:nx,-1:ny), qc(-1:nx,-1:ny),
     *  hx1(-1:nx1,-1:ny1), hy1(-1:nx1,-1:ny1), qc1(-1:nx1,-1:ny1)
      integer bottomiy1(-1:nx1,-1:ny1)
      real (kind=R8) ::
     *  bundle_frac(0:ns1-1,0:ns-1)
      real (kind=R8) ::
     *  na(-1:nx,-1:ny,0:ns-1), ua(-1:nx,-1:ny,0:ns-1),
     &  uadia(-1:nx,-1:ny,0:1,0:ns-1), ne(-1:nx,-1:ny),
     &  te(-1:nx,-1:ny), ti(-1:nx,-1:ny), po(-1:nx,-1:ny)
      real (kind=R8) ::
     *  fna(-1:nx,-1:ny,0:1,0:ns-1), fhe(-1:nx,-1:ny,0:1),
     &  fhi(-1:nx,-1:ny,0:1), fch(-1:nx,-1:ny,0:1),
     6  fch_32(-1:nx,-1:ny,0:1), fch_52(-1:nx,-1:ny,0:1), 
     &  fch_p(-1:nx,-1:ny,0:1), kinrgy(-1:nx,-1:ny,0:ns-1)
*   ..output arguments (unspecified on entry)
      real (kind=R8) ::
     *  na1(-1:nx1,-1:ny1,0:ns1-1), ua1(-1:nx1,-1:ny1,0:ns1-1),
     *  uadia1(-1:nx1,-1:ny1,0:1,0:ns1-1), ne1(-1:nx1,-1:ny1),
     &  te1(-1:nx1,-1:ny1), ti1(-1:nx1,-1:ny1), po1(-1:nx1,-1:ny1)
      real (kind=R8) ::
     *  fna1(-1:nx1,-1:ny1,0:1,0:ns1-1), fhe1(-1:nx1,-1:ny1,0:1),
     &  fhi1(-1:nx1,-1:ny1,0:1), fch1(-1:nx1,-1:ny1,0:1),
     1	fch1_32(-1:nx1,-1:ny1,0:1), fch1_52(-1:nx1,-1:ny1,0:1),
     2  fch1_p(-1:nx1,-1:ny1,0:1), kinrgy1(-1:nx1,-1:ny1,0:ns1-1)

*-----------------------------------------------------------------------
*.documentation
*
*  1. purpose
*
*     B2XTPS transfers the plasma state to a mesh of different size.
*
*
*  3. description (see also routine b2cdca)
*
*     This routine transfers the plasma state to a mesh of different
*     size.  Because of possible cuts in the geometry and associated
*     discontinuities, the lowest order (nearest-neighbour) form of
*     interpolation is employed.
*     In the case of a non-conformal mapping, then we do a 2-D
*     interpolation of values in the computational domain, per region.
*
*
*  5. parameters (see also routine b2cdcv)
*
*     nx, ny, ns, nx1, ny1 - integers, input.
*     ne, na, ua, uadia, te, ti, po, fna, fhe, fhi, fch - real*8 arrays, input.
*     See the calling routine for description.
*
*     ne1, na1, ua1, uadia1, te1, ti1, po1 - real*8 arrays, output.
*     fna1, fhe1, fhi1, fch1 - real*8 arrays, output.
*     Transformed plasma state.
*
*     (parameters in common:)
*
*-----------------------------------------------------------------------
*.declarations

*   ..local variables
      integer is, is1, n2
      real (kind=R8) ::
     *   wrk0(-1:nx,-1:ny1), wrk1(-1:nx1,-1:ny1), na0eps, sum
      save na0eps
      data na0eps/1.0e4_R8/
*   ..procedures
      external subini, subend, xertst, ipgetr, sfill, b2saxpy
      external b2xtvc, b2xtfx, b2xtfy

*-----------------------------------------------------------------------
*.computation

* ..preliminaries
*   ..subprogram start-up calls
      call subini ('b2xtps')
*   ..test nx, ny, ns, nx1, ny1, ns1
      call xertst (0.le.nx.and.0.le.ny, 'faulty argument nx, ny')
      call xertst (1.le.ns, 'faulty argument ns')
      call xertst (0.le.nx1.and.0.le.ny1, 'faulty argument nx1, ny1')
      call xertst (1.le.ns1, 'faulty argument ns1')
*   ..read input parameters
      call ipgetr ('b2mndr_na0eps', na0eps)
      call xertst (na0eps.gt.0.0_R8,
     & 'faulty internal parameter na0eps')
* ..transform state
      n2 = (nx1+2)*(ny1+2)
      call sfill (n2*ns1, 0.0_R8, na1, 1) 
      call sfill (n2*ns1, 0.0_R8, ua1, 1) 
      call sfill (n2*ns1*2, 0.0_R8, uadia1, 1) 
      call sfill (n2*ns1, 0.0_R8, kinrgy1, 1) 
*   ..transform field variables
      do is1 = 0, ns1-1
        sum = 0.0_R8
        do is = 0, ns-1
          sum = sum + bundle_frac(is1,is)
          if (bundle_frac(is1,is).ge.1.0_R8) then  ! simple case
            call b2xtvc (nx, ny, nx1, ny1, na(-1,-1,is), na1(-1,-1,is1),
     &        nncutmax, non_commensurate, 
     &        nnreg, nncut, periodic_bc, jsep, jsep1, region, region1,
     &        leftcut, leftcut1, rightcut, rightcut1,
     &        bottomcut, bottomcut1, topcut, topcut1,
     &        wrk0)
            call b2xtvc (nx, ny, nx1, ny1, ua(-1,-1,is), ua1(-1,-1,is1),
     &        nncutmax, non_commensurate, 
     &        nnreg, nncut, periodic_bc, jsep, jsep1, region, region1,
     &        leftcut, leftcut1, rightcut, rightcut1,
     &        bottomcut, bottomcut1, topcut, topcut1,
     &        wrk0)
            call b2xtvc (nx, ny, nx1, ny1, uadia(-1,-1,0,is), 
     &        uadia1(-1,-1,0,is1),
     &        nncutmax, non_commensurate, 
     &        nnreg, nncut, periodic_bc, jsep, jsep1, region, region1,
     &        leftcut, leftcut1, rightcut, rightcut1,
     &        bottomcut, bottomcut1, topcut, topcut1, wrk0)
            call b2xtvc (nx, ny, nx1, ny1, uadia(-1,-1,1,is), 
     &        uadia1(-1,-1,1,is1),
     &        nncutmax, non_commensurate, 
     &        nnreg, nncut, periodic_bc, jsep, jsep1, region, region1,
     &        leftcut, leftcut1, rightcut, rightcut1,
     &        bottomcut, bottomcut1, topcut, topcut1, wrk0)
            call b2xtvc (nx, ny, nx1, ny1, kinrgy(-1,-1,is), 
     &        kinrgy1(-1,-1,is1),
     &        nncutmax, non_commensurate, 
     &        nnreg, nncut, periodic_bc, jsep, jsep1, region, region1,
     &        leftcut, leftcut1, rightcut, rightcut1,
     &        bottomcut, bottomcut1, topcut, topcut1, wrk0)
          else if (bundle_frac(is1,is).gt.0.0_R8) then
            call b2xtvc (nx, ny, nx1, ny1, na(-1,-1,is), wrk1,
     &        nncutmax, non_commensurate, 
     &        nnreg, nncut, periodic_bc, jsep, jsep1, region, region1,
     &        leftcut, leftcut1, rightcut, rightcut1,
     &        bottomcut, bottomcut1, topcut, topcut1, wrk0)
            call b2saxpy (n2, bundle_frac(is1,is), wrk1, 1, 
     &        na1(-1,-1,is1), 1)
            call b2xtvc (nx, ny, nx1, ny1, ua(-1,-1,is), wrk1,
     &        nncutmax, non_commensurate, 
     &        nnreg, nncut, periodic_bc, jsep, jsep1, region, region1,
     &        leftcut, leftcut1, rightcut, rightcut1,
     &        bottomcut, bottomcut1, topcut, topcut1, wrk0)
            call b2saxpy (n2, bundle_frac(is1,is), wrk1, 1, 
     &        ua1(-1,-1,is1), 1)
            call b2xtvc (nx, ny, nx1, ny1, uadia(-1,-1,0,is), wrk1,
     &        nncutmax, non_commensurate, 
     &        nnreg, nncut, periodic_bc, jsep, jsep1, region, region1,
     &        leftcut, leftcut1, rightcut, rightcut1,
     &        bottomcut, bottomcut1, topcut, topcut1, wrk0)
            call b2saxpy (n2, bundle_frac(is1,is), wrk1, 1, 
     &        uadia1(-1,-1,0,is1), 1)
            call b2xtvc (nx, ny, nx1, ny1, uadia(-1,-1,1,is), wrk1,
     &        nncutmax, non_commensurate, 
     &        nnreg, nncut, periodic_bc, jsep, jsep1, region, region1,
     &        leftcut, leftcut1, rightcut, rightcut1,
     &        bottomcut, bottomcut1, topcut, topcut1, wrk0)
            call b2saxpy (n2, bundle_frac(is1,is), wrk1, 1, 
     &        uadia1(-1,-1,1,is1), 1)
            call b2xtvc (nx, ny, nx1, ny1, kinrgy(-1,-1,is), wrk1,
     &        nncutmax, non_commensurate, 
     &        nnreg, nncut, periodic_bc, jsep, jsep1, region, region1,
     &        leftcut, leftcut1, rightcut, rightcut1,
     &        bottomcut, bottomcut1, topcut, topcut1, wrk0)
            call b2saxpy (n2, bundle_frac(is1,is), wrk1, 1, 
     &        kinrgy1(-1,-1,is1), 1)
          endif
        enddo
        if (sum.eq.0.0_R8) then  ! this is a new species 
          call sfill (n2, na0eps, na1(-1,-1,is1), 1) 
        endif
      enddo
      call b2xtvc (nx, ny, nx1, ny1, ne, ne1,
     &        nncutmax, non_commensurate, 
     &        nnreg, nncut, periodic_bc, jsep, jsep1, region, region1,
     &        leftcut, leftcut1, rightcut, rightcut1,
     &        bottomcut, bottomcut1, topcut, topcut1, wrk0)
      call b2xtvc (nx, ny, nx1, ny1, te, te1,
     &        nncutmax, non_commensurate, 
     &        nnreg, nncut, periodic_bc, jsep, jsep1, region, region1,
     &        leftcut, leftcut1, rightcut, rightcut1,
     &        bottomcut, bottomcut1, topcut, topcut1, wrk0)
      call b2xtvc (nx, ny, nx1, ny1, ti, ti1,
     &        nncutmax, non_commensurate, 
     &        nnreg, nncut, periodic_bc, jsep, jsep1, region, region1,
     &        leftcut, leftcut1, rightcut, rightcut1,
     &        bottomcut, bottomcut1, topcut, topcut1, wrk0)
      call b2xtvc (nx, ny, nx1, ny1, po, po1,
     &        nncutmax, non_commensurate, 
     &        nnreg, nncut, periodic_bc, jsep, jsep1, region, region1,
     &        leftcut, leftcut1, rightcut, rightcut1,
     &        bottomcut, bottomcut1, topcut, topcut1, wrk0)
*   ..transform fluxes
      call sfill (n2*ns1*2, 0.0_R8, fna1, 1) 
      do is = 0, ns-1
        do is1 = 0, ns1-1
          if (bundle_frac(is1,is).ge.1.0_R8) then  ! simple case
            call b2xtfx (nx, ny, nx1, ny1, fna(-1,-1,0,is), 
     &        fna1(-1,-1,0,is1), nncutmax, non_commensurate,
     &        nnreg, nncut, periodic_bc, jsep, jsep1, region, region1,
     &        leftcut, leftcut1, rightcut, rightcut1,
     &        bottomcut, bottomcut1, topcut, topcut1,
     &        hy, qc, hy1, qc1)
            call b2xtfy (nx, ny, nx1, ny1, fna(-1,-1,1,is), 
     &        fna1(-1,-1,1,is1), 
     &        nncutmax, non_commensurate,
     &        nnreg, nncut, periodic_bc, jsep, jsep1, region, region1,
     &        leftcut, leftcut1, rightcut, rightcut1,
     &        bottomcut, bottomcut1, topcut, topcut1, 
     &        hx, hx1, bottomiy1)
          else if (bundle_frac(is1,is).gt.0.0_R8) then
            call b2xtfx (nx, ny, nx1, ny1, fna(-1,-1,0,is), wrk1, 
     &        nncutmax, non_commensurate,
     &        nnreg, nncut, periodic_bc, jsep, jsep1, region, region1,
     &        leftcut, leftcut1, rightcut, rightcut1,
     &        bottomcut, bottomcut1, topcut, topcut1,
     &        hy, qc, hy1, qc1)
            call b2saxpy (n2, bundle_frac(is1,is), wrk1, 1, 
     &        fna1(-1,-1,0,is1), 1)
            call b2xtfy (nx, ny, nx1, ny1, fna(-1,-1,1,is), wrk1, 
     &        nncutmax, non_commensurate,
     &        nnreg, nncut, periodic_bc, jsep, jsep1, region, region1,
     &        leftcut, leftcut1, rightcut, rightcut1,
     &        bottomcut, bottomcut1, topcut, topcut1, 
     .        hx, hx1, bottomiy1)
            call b2saxpy (n2, bundle_frac(is1,is), wrk1, 1, 
     &        fna1(-1,-1,1,is1), 1)
          endif
        enddo
      enddo
      call b2xtfx (nx, ny, nx1, ny1, fhe(-1,-1,0), fhe1(-1,-1,0), 
     &        nncutmax, non_commensurate,
     &        nnreg, nncut, periodic_bc, jsep, jsep1, region, region1,
     &        leftcut, leftcut1, rightcut, rightcut1,
     &        bottomcut, bottomcut1, topcut, topcut1,
     &        hy, qc, hy1, qc1)
      call b2xtfy (nx, ny, nx1, ny1, fhe(-1,-1,1), fhe1(-1,-1,1), 
     &        nncutmax, non_commensurate,
     &        nnreg, nncut, periodic_bc, jsep, jsep1, region, region1,
     &        leftcut, leftcut1, rightcut, rightcut1,
     &        bottomcut, bottomcut1, topcut, topcut1, 
     .        hx, hx1, bottomiy1)
      call b2xtfx (nx, ny, nx1, ny1, fhi(-1,-1,0), fhi1(-1,-1,0), 
     &        nncutmax, non_commensurate,
     &        nnreg, nncut, periodic_bc, jsep, jsep1, region, region1,
     &        leftcut, leftcut1, rightcut, rightcut1,
     &        bottomcut, bottomcut1, topcut, topcut1,
     &        hy, qc, hy1, qc1)
      call b2xtfy (nx, ny, nx1, ny1, fhi(-1,-1,1), fhi1(-1,-1,1), 
     &        nncutmax, non_commensurate,
     &        nnreg, nncut, periodic_bc, jsep, jsep1, region, region1,
     &        leftcut, leftcut1, rightcut, rightcut1,
     &        bottomcut, bottomcut1, topcut, topcut1, 
     .        hx, hx1, bottomiy1)
      call b2xtfx (nx, ny, nx1, ny1, fch(-1,-1,0), fch1(-1,-1,0), 
     &        nncutmax, non_commensurate,
     &        nnreg, nncut, periodic_bc, jsep, jsep1, region, region1,
     &        leftcut, leftcut1, rightcut, rightcut1,
     &        bottomcut, bottomcut1, topcut, topcut1,
     &        hy, qc, hy1, qc1)
      call b2xtfy (nx, ny, nx1, ny1, fch(-1,-1,1), fch1(-1,-1,1), 
     &        nncutmax, non_commensurate,
     &        nnreg, nncut, periodic_bc, jsep, jsep1, region, region1,
     &        leftcut, leftcut1, rightcut, rightcut1,
     &        bottomcut, bottomcut1, topcut, topcut1, 
     .        hx, hx1, bottomiy1)
      call b2xtfx (nx, ny, nx1, ny1, fch_32(-1,-1,0), fch1_32(-1,-1,0), 
     &        nncutmax, non_commensurate,
     &        nnreg, nncut, periodic_bc, jsep, jsep1, region, region1,
     &        leftcut, leftcut1, rightcut, rightcut1,
     &        bottomcut, bottomcut1, topcut, topcut1,
     &        hy, qc, hy1, qc1)
      call b2xtfy (nx, ny, nx1, ny1, fch_32(-1,-1,1), fch1_32(-1,-1,1), 
     &        nncutmax, non_commensurate,
     &        nnreg, nncut, periodic_bc, jsep, jsep1, region, region1,
     &        leftcut, leftcut1, rightcut, rightcut1,
     &        bottomcut, bottomcut1, topcut, topcut1, 
     .        hx, hx1, bottomiy1)
      call b2xtfx (nx, ny, nx1, ny1, fch_52(-1,-1,0), fch1_52(-1,-1,0), 
     &        nncutmax, non_commensurate,
     &        nnreg, nncut, periodic_bc, jsep, jsep1, region, region1,
     &        leftcut, leftcut1, rightcut, rightcut1,
     &        bottomcut, bottomcut1, topcut, topcut1,
     &        hy, qc, hy1, qc1)
      call b2xtfy (nx, ny, nx1, ny1, fch_52(-1,-1,1), fch1_52(-1,-1,1), 
     &        nncutmax, non_commensurate,
     &        nnreg, nncut, periodic_bc, jsep, jsep1, region, region1,
     &        leftcut, leftcut1, rightcut, rightcut1,
     &        bottomcut, bottomcut1, topcut, topcut1, 
     .        hx, hx1, bottomiy1)
      call b2xtfx (nx, ny, nx1, ny1, fch_p(-1,-1,0), fch1_p(-1,-1,0), 
     &        nncutmax, non_commensurate, 
     &        nnreg, nncut, periodic_bc, jsep, jsep1, region, region1,
     &        leftcut, leftcut1, rightcut, rightcut1,
     &        bottomcut, bottomcut1, topcut, topcut1,
     &        hy, qc, hy1, qc1)
      call b2xtfy (nx, ny, nx1, ny1, fch_p(-1,-1,1), fch1_p(-1,-1,1), 
     &        nncutmax, non_commensurate,
     &        nnreg, nncut, periodic_bc, jsep, jsep1, region, region1,
     &        leftcut, leftcut1, rightcut, rightcut1,
     &        bottomcut, bottomcut1, topcut, topcut1, 
     .        hx, hx1, bottomiy1)

* ..return
      call subend ()
      return

*-----------------------------------------------------------------------
*.end b2xtps

      end subroutine b2xtps
