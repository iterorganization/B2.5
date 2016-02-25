*-----------------------------------------------------------------------
*.specification

      subroutine b2xvgm (nx, ny, vol, hx, hy, qz, qc, qcb, gs, pbs,
     &  cell_width, cell_height)
      use b2mod_types
      use b2mod_b2cmgs
      use b2mod_cellhelper
      implicit none
*   ..input arguments (unchanged on exit)
      integer nx, ny
      real (kind=R8) ::
     *  vol(-1:nx,-1:ny), hx(-1:nx,-1:ny), hy(-1:nx,-1:ny),
     &  qz(-1:nx,-1:ny,0:1), qc(-1:nx,-1:ny,0:1), qcb(-1:nx,-1:ny),
     &  gs(-1:nx,-1:ny,0:2), pbs(-1:nx,-1:ny,0:1),
     &  cell_width(-1:nx,-1:ny), cell_height(-1:nx,-1:ny)
*   ..output arguments (unspecified on entry)
*     (none)
*   ..variables in common

*-----------------------------------------------------------------------
*.documentation
*
*  1. purpose
*
*     B2XVGM tests the validity of the geometric parameters.
*
*
*  3. description (see also routine b2cdca)
*
*     This routine tests the validity of the geometric parameters:
*     vol, hx and hy must be positive, qz and qc must be in range, 
*     gs(-1,,0), gs(,-1,1), pbs(-1,,0) and pbs(,-1,1) must vanish, 
*     gs and pbs must otherwise be positive, 
*     and the singular regions in /b2cmgs/ must be valid.
*
*
*  5. parameters (see also routine b2cdcv)
*
*     nx, ny - integer, input.
*     vol, hx, hy, qz, qc, qcb, gs, pbs - real array, input.
*     See the calling routine for description.
*
*     (parameters in common:)
*
*     /b2cmgs/ - input.
*     /b2cmgs/ specifies the singular regions in the geometry.
*     See the calling routine for description.
*
*-----------------------------------------------------------------------
*.declarations

*   ..local variables
      integer n2, i
*   ..procedures
      real (kind=R8) ::
     *   smin, smax
      external subini, subend, xertst, smin, smax
      external b2xvsg, b2xvfx, b2xvfy

*-----------------------------------------------------------------------
*.computation

* ..preliminaries
*   ..subprogram start-up calls
      call subini ('b2xvgm')
*   ..set n2
      n2 = (nx+2)*(ny+2)
*   ..test nx, ny
      call xertst (0.le.nx.and.0.le.ny, 'faulty argument nx, ny')

* ..test the geometry
*   ..test sign of vol, hx, hy
      call b2xvsg (n2, vol, 1, 'vol', '.gt.')
      call b2xvsg (n2, hx, 1, 'hx', '.gt.')
      call b2xvsg (n2, hy, 1, 'hy', '.gt.')
      call b2xvsg (n2, cell_width, 1, 'cell_width', '.gt.')
      call b2xvsg (n2, cell_height, 1, 'cell_height', '.gt.')
*   ..test range of qz, qc
      call xertst (smin(n2,qz(-1,-1,0),1).gt.-1.0_R8.and.
     &  smax(n2,qz(-1,-1,0),1).lt.1.0_R8,
     &  'faulty argument range qz(,,0)')
      call xertst (smin(n2,qz(-1,-1,1),1).gt.0.0_R8.and.
     &  smax(n2,qz(-1,-1,1),1).le.1.0_R8,
     &  'faulty argument range qz(,,1)')
      call xertst (smin(n2,qc,1).gt.0.0_R8.and.
     &  smax(n2,qc,1).le.1.0_R8,
     &  'faulty argument range qc')
      call xertst (smin(n2,qcb,1).gt.0.0_R8.and.
     &  smax(n2,qcb,1).le.1.0_R8,
     &  'faulty argument range qcb')
*   ..test edge values of gs
      call b2xvfx (nx, ny, gs(-1,-1,0), 'gsx')
      call b2xvfy (nx, ny, gs(-1,-1,1), 'gsy')
*   ..test sign of gs
      call b2xvsg (n2, gs(-1,-1,0), 1, 'gsx', '.ge.')
      call b2xvsg (n2, gs(-1,-1,1), 1, 'gsy', '.ge.')
      call b2xvsg (n2, gs(-1,-1,2), 1, 'gsz', '.gt.')
*   ..test edge values of pbs
      call b2xvfx (nx, ny, pbs(-1,-1,0), 'pbsx')
      call b2xvfy (nx, ny, pbs(-1,-1,1), 'pbsy')
*   ..test sign of pbs
!      call b2xvsg (n2, pbs(-1,-1,0), 1, 'pbsx', '.ge.')       !srv 11.01.13
!      call b2xvsg (n2, pbs(-1,-1,1), 1, 'pbsy', '.ge.')       !srv 11.01.13
*   ..test /b2cmgs/
      call xertst (0.le.nlreg.and.nlreg.le.nlrmx,
     &  'faulty argument /b2cmgs/, nlreg')
      do i = 0, nlreg-1
       call xertst (-1.le.nlxlo(i).and.nlxlo(i).le.nlxhi(i).and.
     &   nlxhi(i).le.nx,
     &   'faulty argument /b2cmgs/, range nlxlo, nlxhi')
       call xertst (-1.le.nlylo(i).and.nlylo(i).le.nlyhi(i).and.
     &   nlyhi(i).le.ny,
     &   'faulty argument /b2cmgs/, range nlylo, nlyhi')
       call xertst (0.le.nlloc(i).and.nlloc(i).lt.4,
     &   'faulty argument /b2cmgs/, nlloc')
      enddo

* ..return
      call subend ()
      return

*-----------------------------------------------------------------------
*.end b2xvgm

      end subroutine b2xvgm
*-----------------------------------------------------------------------
*.specification

      subroutine b2xvgm1(nx, ny, vol, hx, hy, qz, qc, qcb, gs, pbs,
     . cell_width, cell_height, bottomiy1)
      use b2mod_types
      use b2mod_b2cmgs
      use b2mod_cellhelper
      implicit none
*   ..input arguments (unchanged on exit)
      integer nx, ny, bottomiy1(-1:nx,-1:ny)
      real (kind=R8) ::
     *  vol(-1:nx,-1:ny), hx(-1:nx,-1:ny), hy(-1:nx,-1:ny),
     &  qz(-1:nx,-1:ny,0:1), qc(-1:nx,-1:ny,0:1), qcb(-1:nx,-1:ny),
     &  gs(-1:nx,-1:ny,0:2), pbs(-1:nx,-1:ny,0:1),
     &  cell_width(-1:nx,-1:ny), cell_height(-1:nx,-1:ny)
*   ..output arguments (unspecified on entry)
*     (none)
*   ..variables in common

*-----------------------------------------------------------------------
*.documentation
*
*  1. purpose
*
*     B2XVGM1 tests the validity of the geometric parameters on the 
*     coarse mesh.
*
*
*  3. description (see also routine b2cdca)
*
*     This routine tests the validity of the geometric parameters:
*     vol, hx and hy must be positive, qz and qc must be in range, 
*     gs(-1,,0), gs(,-1,1), pbs(-1,,0) and pbs(,-1,1) must vanish, 
*     gs and pbs must otherwise be positive, 
*     and the singular regions in /b2cmgs/ must be valid.
*
*
*  5. parameters (see also routine b2cdcv)
*
*     nx, ny - integer, input.
*     vol, hx, hy, qz, qc, qcb, gs, pbs - real array, input.
*     See the calling routine for description.
*
*     (parameters in common:)
*
*     /b2cmgs/ - input.
*     /b2cmgs/ specifies the singular regions in the geometry.
*     See the calling routine for description.
*
*-----------------------------------------------------------------------
*.declarations

*   ..local variables
      integer n2, i
*   ..procedures
      real (kind=R8) ::
     *   smin, smax
      external subini, subend, xertst, smin, smax
      external b2xvsg, b2xvfx, b2xvfy1

*-----------------------------------------------------------------------
*.computation

* ..preliminaries
*   ..subprogram start-up calls
      call subini ('b2xvgm1')
*   ..set n2
      n2 = (nx+2)*(ny+2)
*   ..test nx, ny
      call xertst (0.le.nx.and.0.le.ny, 'faulty argument nx, ny')

* ..test the geometry
*   ..test sign of vol, hx, hy
      call b2xvsg (n2, vol, 1, 'vol', '.gt.')
      call b2xvsg (n2, hx, 1, 'hx', '.gt.')
      call b2xvsg (n2, hy, 1, 'hy', '.gt.')
      call b2xvsg (n2, cell_width, 1, 'cell_width', '.gt.')
      call b2xvsg (n2, cell_height, 1, 'cell_height', '.gt.')
*   ..test range of qz, qc
      call xertst (smin(n2,qz(-1,-1,0),1).gt.-1.0_R8.and.
     &  smax(n2,qz(-1,-1,0),1).lt.1.0_R8,
     &  'faulty argument range qz(,,0)')
      call xertst (smin(n2,qz(-1,-1,1),1).gt.0.0_R8.and.
     &  smax(n2,qz(-1,-1,1),1).le.1.0_R8,
     &  'faulty argument range qz(,,1)')
      call xertst (smin(n2,qc,1).gt.0.0_R8.and.
     &  smax(n2,qc,1).le.1.0_R8,
     &  'faulty argument range qc')
      call xertst (smin(n2,qcb,1).gt.0.0_R8.and.
     &  smax(n2,qcb,1).le.1.0_R8,
     &  'faulty argument range qcb')
*   ..test edge values of gs
      call b2xvfx (nx, ny, gs(-1,-1,0), 'gsx')
      call b2xvfy1(nx, ny, gs(-1,-1,1), 'gsy', bottomiy1)
*   ..test sign of gs
      call b2xvsg (n2, gs(-1,-1,0), 1, 'gsx', '.ge.')
      call b2xvsg (n2, gs(-1,-1,1), 1, 'gsy', '.ge.')
      call b2xvsg (n2, gs(-1,-1,2), 1, 'gsz', '.gt.')
*   ..test edge values of pbs
      call b2xvfx (nx, ny, pbs(-1,-1,0), 'pbsx')
      call b2xvfy1(nx, ny, pbs(-1,-1,1), 'pbsy', bottomiy1)
*   ..test sign of pbs
!      call b2xvsg (n2, pbs(-1,-1,0), 1, 'pbsx', '.ge.')         !srv 11.01.13
!      call b2xvsg (n2, pbs(-1,-1,1), 1, 'pbsy', '.ge.')         !srv 11.01.13
*   ..test /b2cmgs/
      call xertst (0.le.nlreg.and.nlreg.le.nlrmx,
     &  'faulty argument /b2cmgs/, nlreg')
      do i = 0, nlreg-1
       call xertst (-1.le.nlxlo(i).and.nlxlo(i).le.nlxhi(i).and.
     &   nlxhi(i).le.nx,
     &   'faulty argument /b2cmgs/, range nlxlo, nlxhi')
       call xertst (-1.le.nlylo(i).and.nlylo(i).le.nlyhi(i).and.
     &   nlyhi(i).le.ny,
     &   'faulty argument /b2cmgs/, range nlylo, nlyhi')
       call xertst (0.le.nlloc(i).and.nlloc(i).lt.4,
     &   'faulty argument /b2cmgs/, nlloc')
      enddo

* ..return
      call subend ()
      return

*-----------------------------------------------------------------------
*.end b2xvgm1

      end subroutine b2xvgm1
