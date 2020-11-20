*-----------------------------------------------------------------------
*.specification

      subroutine b2agmx (nx1, ny1, nx, ny, 
     & vol1, hxx, hyy, qz1, qc1, qcb1,
     & gs1, crx1, cry1, bb1, pbs1, fpsi1, ffbz1,
     & wbbl1, wbbr1, wbbv1, wbbc1,
     & cell_width1, cell_height1,
     & vol, hx, hy, qz, qc, qcb,
     & gs, crx, cry, bb, pbs, fpsi, ffbz,
     & wbbl, wbbr, wbbv, wbbc,
     & cell_width, cell_height,
     1 leftix1, leftiy1, rightix1, rightiy1,
     2 topix1, topiy1, bottomix1, bottomiy1,
     3 leftcut1, rightcut1, bottomcut1, topcut1,
     $ cflags1)
      use b2mod_types
      use b2mod_math
      use b2mod_b2cmfs
      use b2mod_indirect
      use b2mod_constants
      use b2mod_cellhelper
      use b2mod_connectivity
      use b2mod_geo2
     , , only : dist1
      implicit none
*   ..input arguments (unchanged on exit)
      integer nx1, ny1, nx, ny
      real (kind=R8) ::
     *  vol1(-1:nx1,-1:ny1), hxx(-1:nx1,-1:ny1), hyy(-1:nx1,-1:ny1),
     &  qz1(-1:nx1,-1:ny1,0:1), qc1(-1:nx1,-1:ny1), qcb1(-1:nx1,-1:ny1),
     &  gs1(-1:nx1,-1:ny1,0:2),
     1	crx1(-1:nx1,-1:ny1,0:3), cry1(-1:nx1,-1:ny1,0:3),
     2	bb1(-1:nx1,-1:ny1,0:3), pbs1(-1:nx1,-1:ny1,0:1),
     3  fpsi1(-1:nx1,-1:ny1,0:3), ffbz1(-1:nx1,-1:ny1,0:3),
     4  wbbl1(-1:nx1,-1:ny1,0:3), wbbr1(-1:nx1,-1:ny1,0:3),
     &  wbbv1(-1:nx1,-1:ny1,0:3), wbbc1(-1:nx1+1,-1:ny1+1,0:3),
     5  cell_width1(-1:nx1,-1:ny1), cell_height1(-1:nx1,-1:ny1)
      integer
     *  leftix1(-1:nx1,-1:ny1), leftiy1(-1:nx1,-1:ny1),
     1	rightix1(-1:nx1,-1:ny1), rightiy1(-1:nx1,-1:ny1),
     2	topix1(-1:nx1,-1:ny1), topiy1(-1:nx1,-1:ny1),
     3	bottomix1(-1:nx1,-1:ny1), bottomiy1(-1:nx1,-1:ny1),
     4	leftcut1(1:nncutmax), rightcut1(1:nncutmax),
     5  bottomcut1(1:nncutmax), topcut1(1:nncutmax)
      integer, intent(in) :: cflags1(-1:nx1,-1:ny1, CARREOUT_NCELLFLAGS)
*   ..output arguments (unspecified on entry)
      real (kind=R8) ::
     *  vol(-1:nx,-1:ny), hx(-1:nx,-1:ny), hy(-1:nx,-1:ny),
     &  qz(-1:nx,-1:ny,0:1), qc(-1:nx,-1:ny), qcb(-1:nx,-1:ny),
     1	gs(-1:nx,-1:ny,0:2), crx(-1:nx,-1:ny,0:3), cry(-1:nx,-1:ny,0:3),
     2	bb(-1:nx,-1:ny,0:3), pbs(-1:nx,-1:ny,0:1),
     3  fpsi(-1:nx,-1:ny,0:3), ffbz(-1:nx,-1:ny,0:3),
     4  wbbl(-1:nx,-1:ny,0:3), wbbr(-1:nx,-1:ny,0:3),
     &  wbbv(-1:nx,-1:ny,0:3), wbbc(-1:nx+1,-1:ny+1,0:3),
     5  cell_width(-1:nx,-1:ny), cell_height(-1:nx,-1:ny)

*-----------------------------------------------------------------------
*.documentation
*
*  1. purpose
*
*     B2AGMX projects the metric coefficients from a fine mesh to a
*     coarsening of that mesh.
*
*
*  3. description (see also routine b2cdca)
*
*     This routine takes as input the metric coefficients on a fine
*     mesh: vol1, hxx, hyy, qz1, qc1 and gs1.  It computes the metric
*     coefficients on a corresponding coarsened mesh.
*
*
*  5. parameters (see also routine b2cdcv)
*
*     nx1, ny1, nx, ny - integer, input.
*     (nx1,ny1) are the dimensions of the fine mesh; (nx,ny) are the
*     dimensions of the coarse mesh.
*     See the calling routine for description.
*
*     vol1, hxx, hyy, qz1, qc1, gs1, pbs1 - real*8 array, input.
*     These arrays hold the metric coefficients on the fine mesh.
*     See below for detailed description.
*
*     vol - (-1:nx,-1:ny) real*8 array, output.
*     For (ix,iy) in (-1:nx,-1:ny), vol(ix,iy) specifies the volume of
*     the (ix,iy) cell. (In the case of toroidal symmetry this volume
*     includes a factor 2*pi*crx, while in the case of cylindrical
*     symmetry it includes a factor 2*pi*cry.)
*     It will hold that 0.lt.vol(,).
*
*     hx, hy - (-1:nx,-1:ny) real*8 array, output.
*     For (ix,iy) in (-1:nx,-1:ny), hx(ix,iy) and hy(ix,iy) specify the
*     x-diameter and the y-diameter, respectively, of the (ix,iy) cell.
*     It will hold that 0.lt.hx(,) and 0.lt.hy(,).
*
*     cell_width, cell_height - (-1:nx,-1:ny) real*8 array, output.
*     For (ix,iy) in (-1:nx,-1:ny), cell_width(ix,iy) and cell_height(ix,iy)
*     specify the average width and height, respectively, of the (ix,iy) cell.
*     It will hold that 0.lt.cell_width(,) and 0.lt.cell_height(,).
*
*     qz - (-1:nx,-1:ny,0:1) real*8 array, output.
*     For (ix,iy) in (-1:nx,-1:ny), qz(ix,iy,0) specifies sin(t) and
*     qz(ix,iy,1) specifies cos(t), where t is the complementary angle
*     between the x-direction and the y-direction on the (ix,iy) cell.
*     t is in the range (-pi/2 to pi/2) and vanishes for an orthogonal
*     geometry. (qz(ix,iy,0) is the inner product of the x and y unit
*     vectors.)
*     It will hold that -1.lt.qz(,,0).lt.1 and 0.lt.qz(,,1).le.1.
*
*     qc - (-1:nx,-1:ny) real*8 array, output.
*     For (ix,iy) in (-1:nx,-1:ny), qc(ix,iy) specifies cos(t), where 
*     t is the complementary angle between the x-direction of the 
*     (ix,iy) cell and its left face.
*     t is in the range (-pi/2 to pi/2) and vanishes for an orthogonal
*     geometry. 
*     It will hold that 0.lt.qc(,).le.1.
*
*     qcb - (-1:nx,-1:ny) real*8 array, output.
*     For (ix,iy) in (-1:nx,-1:ny), qcb(ix,iy) specifies cos(t), where
*     t is the complementary angle between the y-direction of the
*     (ix,iy) cell and its bottom face.
*     t is in the range (-pi/2 to pi/2) and vanishes for an orthogonal
*     geometry.
*     It will hold that 0.lt.qcb(,).le.1.
*
*     gs - (-1:nx,-1:ny,0:2) real*8 array, output.
*     (Let gsx(,)=gs(,,0), gsy(,)=gs(,,1) and gsz(,)=gs(,,2).)
*     For (ix,iy) in (0:nx,-1:ny), gsx(ix,iy) specifies the area of the
*     face between the (ix,iy) cell and its left neighbor. For ix.eq.-1,
*     gsx(ix,-1:ny) will hold 0.
*     For (ix,iy) in (-1:nx,0:ny), gsy(ix,iy) specifies the area of the
*     face between the (ix,iy) cell and its bottom neighbor. For iy.eq.-1,
*     gsy(-1:nx,iy) will hold 0.
*     For (ix,iy) in (-1:nx,-1:ny), gsz(ix,iy) specifies the area of the
*     (ix,iy) cell normal to the third (ignorable) coordinate.
*
*     pbs - (-1:nx,-1:ny,0:1) real*8 array.
*     (Let pbsx(,)=pbs(,,0) and pbsy(,)=pbs(,,1).)
*     For (ix,iy) in (-1:nx,-1:ny), pbsx(ix,iy) specifies the product
*     (bx/bb)*sx and pbsy(ix,iy) specifies (by/bb)*sy, or (magnetic 
*     field pitch)*(surface area), namely the parallel contact area, 
*     on the (ix,iy) cell left and bottom faces respectively.  
*     See routine b2xgbs for the form of interpolation employed.
*
*-----------------------------------------------------------------------
*.declarations

*   ..local variables
      integer ic, isox, isoy, ix, iy, ix1, iy1, jx1, jy1, kx1, ky1, 
     *   nxeff, nx1eff, ixiso(0:nncutmax), ix1iso(0:nncutmax),
     *   nyeff, ny1eff, iyiso(0:nncutmax), iy1iso(0:nncutmax)
      real (kind=R8) ::
     *   t0, t1, t2, t3, t4, t5
      logical iniso, error
      integer pbs_from_basis_mesh
*   ..procedures
      intrinsic sqrt
      external subini, subend, xertst, ipgeti
      external b2xvgm1
*   ..initialisation
      save pbs_from_basis_mesh
      data pbs_from_basis_mesh /1/

*-----------------------------------------------------------------------
*.computation

* ..preliminaries
*   ..subprogram start-up calls
      call subini ('b2agmx')
*   ..test nx1, ny1, nx, ny
      call xertst (1.le.nx.and.1.le.ny, 'faulty argument nx, ny')
      call xertst (1.le.nx1.and.1.le.ny1, 'faulty input nx1, ny1')
*   ..test input metric
      call b2xvgm1(nx1, ny1, vol1, hxx, hyy, qz1, qc1, qcb1, gs1, pbs1,
     . cell_width1, cell_height1, bottomiy1)
*   ..read input switch
      call ipgeti ('b2agmx_pbs_from_basis_mesh', pbs_from_basis_mesh)

* ..compute metric coefficients on coarse mesh
*   ..for multiple-null case, identify inner isolating cuts
      ix1iso(:) = nx1+1
      ix1iso(0) = -1
      isox = 0
      do ix1 = 0, nx1-1
        iniso = .false.
        do iy1 = 0, ny1-1
          if (.not.isInDomain(nx1,ny1,leftix1(ix1,iy1),leftiy1(ix1,iy1))
     &        .and. .not. iniso) then
            iniso = .true.
            isox = isox+1
            ix1iso(isox) = ix1
          endif
        enddo
      enddo
*   ..obtain nxeff, nx1eff
      nxeff = nx - 2*isox
      nx1eff= nx1- 2*isox
*   ..test nxeff, nx1eff
      call xertst (1.le.nxeff.and.1.le.nx1eff, 
     &  'faulty values nxeff, nx1eff')
      call xertst (nxeff*(nx1eff/nxeff).eq.nx1eff, 
     & 'nx1eff/nxeff not integer')
*   ..for multiple-null case, inner isolating cuts, mapped to new mesh
      ixiso(:) = nx+1
      ixiso(0) = -1
      do ic = 1, isox
        ixiso(ic) = floor(1.0_R8*(ix1iso(ic)-(2*ic-1))*nxeff/nx1eff) 
     .   + 2*ic-1
      enddo
*   ..for limiter case, identify horizontal inner boundaries
      iy1iso(:) = ny1+1
      iy1iso(0) = -1
      isoy = 0
      do iy1 = 0, ny1-1
        iniso = .false.
        do ix1 = 0, nx1-1
          if ((.not.isInDomain(nx1,ny1,topix1(ix1,iy1),topiy1(ix1,iy1))
     &     .or. .not.isInDomain(nx1,ny1,bottomix1(ix1,iy1),
     &                                  bottomiy1(ix1,iy1)))
     &     .and. (.not. iniso)) then  
            iniso = .true.
            isoy = isoy+1
            iy1iso(isoy) = iy1
          endif
        enddo
      enddo
*   ..obtain nyeff, ny1eff
      nyeff = ny - isoy
      ny1eff= ny1- isoy
*   ..test nyeff, ny1eff
      call xertst (1.le.nyeff.and.1.le.ny1eff, 
     &  'faulty values nyeff, ny1eff')
      call xertst (nyeff*(ny1eff/nyeff).eq.ny1eff, 
     & 'ny1eff/nyeff not integer')
*   ..for limiter case, horizontal boundaries, mapped to new mesh
      iyiso(:) = ny+1
      iyiso(0) = -1
      do ic = 1, isoy
        iyiso(ic) = floor(1.0_R8*(iy1iso(ic)-(ic-1))*nyeff/ny1eff) 
     .   + ic
      enddo
*   ..loop over cells
      isoy = 0
      do iy = -1, ny
        if (iy.gt.iyiso(isoy).and.iy.lt.iyiso(isoy+1)-1) then   ! inner cells
          ky1 = ny1eff/nyeff
          iy1 = (iy-isoy)*ny1eff/nyeff + isoy
        else if (iy.eq.iyiso(isoy)) then          ! lower edges
          ky1 = 1
          iy1 = (iy-(isoy-1))*ny1eff/nyeff + (isoy-1)
        else if (iy.eq.iyiso(isoy+1)-1) then      ! top edges
          ky1 = 1
          iy1 = (iy-isoy)*ny1eff/nyeff + isoy
        endif
        isox = 0
        do ix = -1, nx
          if (ix.gt.ixiso(isox).and.ix.lt.ixiso(isox+1)-1) then   ! inner cells
            kx1 = nx1eff/nxeff
            ix1 = (ix-2*isox)*nx1eff/nxeff + 2*isox 
          elseif (ix.eq.ixiso(isox)) then           ! left edges
            kx1 = 1
            ix1 = (ix-(2*isox-1))*nx1eff/nxeff + (2*isox-1)
          else if (ix.eq.ixiso(isox+1)-1) then      ! right edges
            kx1 = 1
            ix1 = (ix-2*isox)*nx1eff/nxeff + 2*isox
          endif
*     ..project gsx, pbsx and qc to cell (ix,iy)
          if (.not.isInDomain(nx1,ny1,leftix1(ix1,iy1),
     .                                leftiy1(ix1,iy1))) then
            gs(ix,iy,0) = 0.0_R8
            pbs(ix,iy,0) = 0.0_R8
            t0 = 0.0_R8
            t1 = 0.0_R8
            do jy1 = 0, ky1-1
              t0 = t0+pi*dist1(ix1,iy1+jy1,0,2)*
     .         (crx1(ix1,iy1+jy1,0)+crx1(ix1,iy1+jy1,2))
              t1 = t1+qc1(ix1,iy1+jy1)*pi*dist1(ix1,iy1+jy1,0,2)*
     .         (crx1(ix1,iy1+jy1,0)+crx1(ix1,iy1+jy1,2))
            enddo
            qc(ix,iy) = t1/t0
          else
            t0 = 0.0_R8
            t1 = 0.0_R8
            t2 = 0.0_R8
            do jy1 = 0, ky1-1
              t0 = t0+gs1(ix1,iy1+jy1,0)
              t1 = t1+qc1(ix1,iy1+jy1)*gs1(ix1,iy1+jy1,0)
              t2 = t2+pbs1(ix1,iy1+jy1,0)
            enddo
            gs(ix,iy,0) = t0
            qc(ix,iy) = t1/t0
            pbs(ix,iy,0) = t2
          endif
*     ..project gsy, pbsy and qcb to cell (ix,iy)
          if (.not.isInDomain(nx1,ny1,bottomix1(ix1,iy1),
     .                                bottomiy1(ix1,iy1))) then
            gs(ix,iy,1) = 0.0_R8
            pbs(ix,iy,1) = 0.0_R8
            t0 = 0.0_R8
            t1 = 0.0_R8
            do jx1 = 0, kx1-1
              t0 = t0+pi*dist1(ix1+jx1,iy1,0,1)*
     .         (crx1(ix1+jx1,iy1,0)+crx1(ix1+jx1,iy1,1))
              t1 = t1+qcb1(ix1+jx1,iy1)*pi*dist1(ix1+jx1,iy1,0,1)*
     .         (crx1(ix1+jx1,iy1,0)+crx1(ix1+jx1,iy1,1))
            enddo
            qcb(ix,iy) = t1/t0
          else
            t0 = 0.0_R8
            t1 = 0.0_R8
            t2 = 0.0_R8
            do jx1 = 0, kx1-1
              t0 = t0+gs1(ix1+jx1,iy1,1)
              t1 = t1+qcb1(ix1+jx1,iy1)*gs1(ix1+jx1,iy1,1)
              t2 = t2+pbs1(ix1+jx1,iy1,1)
            enddo
            gs(ix,iy,1) = t0
            qcb(ix,iy) = t1/t0
            pbs(ix,iy,1) = t2
          endif
*     ..project vol, gsz, and qz to cell (ix,iy)
          t0 = 0.0_R8
          t1 = 0.0_R8
          t2 = 0.0_R8
cdpc
          t3 = 0.0_R8
          t4 = 0.0_R8
          t5 = 0.0_R8
cdpc-end
          do jy1 = 0, ky1-1
            do jx1 = 0, kx1-1
              t0 = t0+vol1(ix1+jx1,iy1+jy1)
              t1 = t1+gs1(ix1+jx1,iy1+jy1,2)
              t2 = t2+gs1(ix1+jx1,iy1+jy1,2)*qz1(ix1+jx1,iy1+jy1,0)
cdpc
              t3 = t3+bb1(ix1+jx1,iy1+jy1,0)*vol1(ix1+jx1,iy1+jy1)
              t4 = t4+bb1(ix1+jx1,iy1+jy1,1)*vol1(ix1+jx1,iy1+jy1)
              t5 = t5+bb1(ix1+jx1,iy1+jy1,2)*vol1(ix1+jx1,iy1+jy1)
cdpc-end
            enddo
          enddo
          vol(ix,iy) = t0
          gs(ix,iy,2) = t1
          qz(ix,iy,0) = t2/t1
          qz(ix,iy,1) = sqrt(1.0_R8-qz(ix,iy,0)**2)
cdpc
          bb(ix,iy,0)=t3/t0
          bb(ix,iy,1)=t4/t0
          bb(ix,iy,2)=t5/t0
          bb(ix,iy,3)=sqrt(bb(ix,iy,0)**2+bb(ix,iy,1)**2+bb(ix,iy,2)**2)
          crx(ix,iy,0)=crx1(ix1,iy1,0)
          crx(ix,iy,1)=crx1(ix1+kx1-1,iy1,1)
          crx(ix,iy,2)=crx1(ix1,iy1+ky1-1,2)
          crx(ix,iy,3)=crx1(ix1+kx1-1,iy1+ky1-1,3)
          cry(ix,iy,0)=cry1(ix1,iy1,0)
          cry(ix,iy,1)=cry1(ix1+kx1-1,iy1,1)
          cry(ix,iy,2)=cry1(ix1,iy1+ky1-1,2)
          cry(ix,iy,3)=cry1(ix1+kx1-1,iy1+ky1-1,3)
          fpsi(ix,iy,0)=fpsi1(ix1,iy1,0)
          fpsi(ix,iy,1)=fpsi1(ix1+kx1-1,iy1,1)
          fpsi(ix,iy,2)=fpsi1(ix1,iy1+ky1-1,2)
          fpsi(ix,iy,3)=fpsi1(ix1+kx1-1,iy1+ky1-1,3)
          ffbz(ix,iy,0)=ffbz1(ix1,iy1,0)
          ffbz(ix,iy,1)=ffbz1(ix1+kx1-1,iy1,1)
          ffbz(ix,iy,2)=ffbz1(ix1,iy1+ky1-1,2)
          ffbz(ix,iy,3)=ffbz1(ix1+kx1-1,iy1+ky1-1,3)
          wbbc(ix,iy,0:3)=wbbc1(ix1,iy1,0:3)
          if (ix.eq.nx .and. iy.lt.ny) then
            wbbc(nx+1,iy,0:3)=wbbc1(nx1+1,iy1,0:3)
          else if (ix.lt.nx .and. iy.eq.ny) then
            wbbc(ix,ny+1,0:3)=wbbc1(ix1,ny1+1,0:3)
          else if (ix.eq.nx .and. iy.eq.ny) then
            wbbc(nx+1,ny+1,0:3)=wbbc1(nx1+1,ny1+1,0:3)
          endif
cxpb Initial import of cflags
          cflags(ix,iy,CELLFLAG_TYPE)=cflags1(ix1,iy1,CELLFLAG_TYPE)
cxpb Indicate we have done some coarsening
          if (cflags(ix,iy,CELLFLAG_TYPE) == GRID_INTERNAL_COARSEN)
     >     cflags(ix,iy,CELLFLAG_TYPE) = GRID_INTERNAL
cxpb Now check boundary cells, if they are in the right direction,
cxpb and if they needed coarsening
          if (cflags1(ix1,iy1,CELLFLAG_TYPE) == GRID_BOUNDARY_COARSEN
     &       .and. cflags1(ix1,iy1,CELLFLAG_BOTTOMFACE).ne.0 .and.
     &       kx1.gt.1)
     >     cflags(ix,iy,CELLFLAG_TYPE) = GRID_BOUNDARY
          if (cflags1(ix1,iy1+ky1-1,CELLFLAG_TYPE) ==
     &                                     GRID_BOUNDARY_COARSEN
     &       .and. cflags1(ix1,iy1,CELLFLAG_TOPFACE).ne.0 .and.
     &       kx1.gt.1)
     >     cflags(ix,iy,CELLFLAG_TYPE) = GRID_BOUNDARY
          if (cflags1(ix1,iy1,CELLFLAG_TYPE) == GRID_BOUNDARY_COARSEN
     &       .and. cflags1(ix1,iy1,CELLFLAG_LEFTFACE).ne.0 .and.
     &       ky1.gt.1)
     >     cflags(ix,iy,CELLFLAG_TYPE) = GRID_BOUNDARY
          if (cflags1(ix1+kx1-1,iy1,CELLFLAG_TYPE) ==
     &                                     GRID_BOUNDARY_COARSEN
     &       .and. cflags1(ix1,iy1,CELLFLAG_RIGHTFACE).ne.0 .and. 
     &       ky1.gt.1)
     >     cflags(ix,iy,CELLFLAG_TYPE) = GRID_BOUNDARY    
cxpb We now check type compatibilities
          error = .false.
          do jx1 = 0, kx1-1
            if (error) exit  ! no need to look further
            do jy1 = 0, ky1-1
              if (error) exit  ! no need to look further
              if (jx1.eq.0 .and. jy1.eq.0) cycle  ! initial data
              if (cflags(ix,iy,CELLFLAG_TYPE) == 
     .            cflags1(ix1+jx1,iy1+jy1,CELLFLAG_TYPE)) cycle  ! OK case
              TYPE: select case (cflags1(ix1+jx1,iy1+jy1,CELLFLAG_TYPE))
                case(GRID_UNDEFINED)
                  error = .not. (
     &                    cflags(ix,iy,CELLFLAG_TYPE) == GRID_EXTERNAL
     &               .or. cflags(ix,iy,CELLFLAG_TYPE) == GRID_DEAD)
                case(GRID_INTERNAL)
                  error = cflags(ix,iy,CELLFLAG_TYPE) == GRID_UNDEFINED
     &               .or. cflags(ix,iy,CELLFLAG_TYPE) == GRID_EXTERNAL
     &               .or. cflags(ix,iy,CELLFLAG_TYPE) == GRID_GUARD
     &               .or. cflags(ix,iy,CELLFLAG_TYPE) == GRID_DEAD
                case(GRID_EXTERNAL)
                  error = .not.(
     &                    cflags(ix,iy,CELLFLAG_TYPE) == GRID_UNDEFINED
     &               .or. cflags(ix,iy,CELLFLAG_TYPE) == GRID_GUARD
     &               .or. cflags(ix,iy,CELLFLAG_TYPE) == GRID_DEAD)
                  if (cflags(ix,iy,CELLFLAG_TYPE) == GRID_UNDEFINED)
     >                cflags(ix,iy,CELLFLAG_TYPE) = GRID_EXTERNAL
                case(GRID_BOUNDARY)
                  error = cflags(ix,iy,CELLFLAG_TYPE) == GRID_UNDEFINED
     &               .or. cflags(ix,iy,CELLFLAG_TYPE) == GRID_EXTERNAL
     &               .or. cflags(ix,iy,CELLFLAG_TYPE) == GRID_GUARD
     &               .or. cflags(ix,iy,CELLFLAG_TYPE) == GRID_DEAD
                  if (cflags(ix,iy,CELLFLAG_TYPE) == GRID_INTERNAL)
     >                cflags(ix,iy,CELLFLAG_TYPE) = GRID_BOUNDARY
                case(GRID_BOUNDARY_REFINE)
                  error = cflags(ix,iy,CELLFLAG_TYPE) == GRID_UNDEFINED
     &               .or. cflags(ix,iy,CELLFLAG_TYPE) == GRID_EXTERNAL
     &               .or. cflags(ix,iy,CELLFLAG_TYPE) == GRID_GUARD
     &               .or. cflags(ix,iy,CELLFLAG_TYPE) == GRID_DEAD
                  if (cflags(ix,iy,CELLFLAG_TYPE) == GRID_INTERNAL .or.
     &        cflags(ix,iy,CELLFLAG_TYPE) == GRID_INTERNAL_COARSEN .or.
     &                cflags(ix,iy,CELLFLAG_TYPE) == GRID_BOUNDARY .or.
     &        cflags(ix,iy,CELLFLAG_TYPE) == GRID_BOUNDARY_COARSEN .or.
     &                cflags(ix,iy,CELLFLAG_TYPE) == GRID_REFINE)
     >                cflags(ix,iy,CELLFLAG_TYPE) = GRID_BOUNDARY_REFINE
                case(GRID_BOUNDARY_REFINE_FIX)
                  error = cflags(ix,iy,CELLFLAG_TYPE) == GRID_UNDEFINED
     &               .or. cflags(ix,iy,CELLFLAG_TYPE) == GRID_EXTERNAL
     &               .or. cflags(ix,iy,CELLFLAG_TYPE) == GRID_GUARD
     &               .or. cflags(ix,iy,CELLFLAG_TYPE) == GRID_DEAD
                  if (cflags(ix,iy,CELLFLAG_TYPE) == GRID_INTERNAL .or.
     &        cflags(ix,iy,CELLFLAG_TYPE) == GRID_INTERNAL_COARSEN .or.
     &                cflags(ix,iy,CELLFLAG_TYPE) == GRID_BOUNDARY .or.
     &        cflags(ix,iy,CELLFLAG_TYPE) == GRID_BOUNDARY_COARSEN .or.
     &         cflags(ix,iy,CELLFLAG_TYPE) == GRID_BOUNDARY_REFINE .or.
     &                cflags(ix,iy,CELLFLAG_TYPE) == GRID_REFINE)
     >            cflags(ix,iy,CELLFLAG_TYPE) = GRID_BOUNDARY_REFINE_FIX
                case(GRID_INTERNAL_COARSEN)
                  error = cflags(ix,iy,CELLFLAG_TYPE) == GRID_UNDEFINED
     &               .or. cflags(ix,iy,CELLFLAG_TYPE) == GRID_EXTERNAL
     &               .or. cflags(ix,iy,CELLFLAG_TYPE) == GRID_GUARD
     &               .or. cflags(ix,iy,CELLFLAG_TYPE) == GRID_DEAD
                case(GRID_BOUNDARY_COARSEN)
                  error = cflags(ix,iy,CELLFLAG_TYPE) == GRID_UNDEFINED
     &               .or. cflags(ix,iy,CELLFLAG_TYPE) == GRID_EXTERNAL
     &               .or. cflags(ix,iy,CELLFLAG_TYPE) == GRID_GUARD
     &               .or. cflags(ix,iy,CELLFLAG_TYPE) == GRID_DEAD
                  if (cflags(ix,iy,CELLFLAG_TYPE) == GRID_INTERNAL .or.
     &             cflags(ix,iy,CELLFLAG_TYPE) == GRID_INTERNAL_COARSEN)
     >             cflags(ix,iy,CELLFLAG_TYPE) = GRID_BOUNDARY
                case(GRID_REFINE)
                  error = cflags(ix,iy,CELLFLAG_TYPE) == GRID_UNDEFINED
                  if (cflags(ix,iy,CELLFLAG_TYPE) == GRID_INTERNAL .or.
     &                cflags(ix,iy,CELLFLAG_TYPE) == GRID_EXTERNAL .or.
     &        cflags(ix,iy,CELLFLAG_TYPE) == GRID_INTERNAL_COARSEN .or.
     &                cflags(ix,iy,CELLFLAG_TYPE) == GRID_GUARD .or.
     &                cflags(ix,iy,CELLFLAG_TYPE) == GRID_DEAD)
     >             cflags(ix,iy,CELLFLAG_TYPE) = GRID_REFINE
                  if (cflags(ix,iy,CELLFLAG_TYPE) == GRID_BOUNDARY .or.
     &             cflags(ix,iy,CELLFLAG_TYPE) == GRID_BOUNDARY_COARSEN)
     >             cflags(ix,iy,CELLFLAG_TYPE) = GRID_BOUNDARY_REFINE
                case(GRID_GUARD)
                  error = cflags(ix,iy,CELLFLAG_TYPE) == GRID_UNDEFINED
     &               .or. cflags(ix,iy,CELLFLAG_TYPE) == GRID_INTERNAL
     &               .or. cflags(ix,iy,CELLFLAG_TYPE) == GRID_BOUNDARY
     &        .or. cflags(ix,iy,CELLFLAG_TYPE) == GRID_BOUNDARY_REFINE
     &    .or. cflags(ix,iy,CELLFLAG_TYPE) == GRID_BOUNDARY_REFINE_FIX
     &       .or. cflags(ix,iy,CELLFLAG_TYPE) == GRID_INTERNAL_COARSEN
     &       .or. cflags(ix,iy,CELLFLAG_TYPE) == GRID_BOUNDARY_COARSEN
                  if (cflags(ix,iy,CELLFLAG_TYPE) == GRID_EXTERNAL .or.
     &                cflags(ix,iy,CELLFLAG_TYPE) == GRID_DEAD)
     >             cflags(ix,iy,CELLFLAG_TYPE) = GRID_GUARD
                case(GRID_DEAD)
                  error = cflags(ix,iy,CELLFLAG_TYPE) == GRID_INTERNAL
     &               .or. cflags(ix,iy,CELLFLAG_TYPE) == GRID_BOUNDARY
     &        .or. cflags(ix,iy,CELLFLAG_TYPE) == GRID_BOUNDARY_REFINE
     &    .or. cflags(ix,iy,CELLFLAG_TYPE) == GRID_BOUNDARY_REFINE_FIX
     &       .or. cflags(ix,iy,CELLFLAG_TYPE) == GRID_INTERNAL_COARSEN
     &       .or. cflags(ix,iy,CELLFLAG_TYPE) == GRID_BOUNDARY_COARSEN
                  if (cflags(ix,iy,CELLFLAG_TYPE) == GRID_UNDEFINED)
     >             cflags(ix,iy,CELLFLAG_TYPE) = GRID_DEAD
                case default TYPE
                  error = .true.  ! we should not get here
              end select TYPE
            end do
          end do
          call xertst(.not.error,
     .     'Cell merging failed because of mismatched types')
cxpb  We import face data, which must match
          cflags(ix,iy,CELLFLAG_LEFTFACE) =
     .     cflags1(ix1,iy1,CELLFLAG_LEFTFACE)
          do jy1 = 1, ky1-1
            call xertst(
     .       cflags1(ix1,iy1+jy1,CELLFLAG_LEFTFACE).eq.
     .       cflags(ix,iy,CELLFLAG_LEFTFACE),
     .       'Boundary number mismatch on left face')
          end do
          cflags(ix,iy,CELLFLAG_BOTTOMFACE) =
     .     cflags1(ix1,iy1,CELLFLAG_BOTTOMFACE)
          do jx1 = 1, kx1-1
            call xertst(
     .       cflags1(ix1+jx1,iy1,CELLFLAG_BOTTOMFACE).eq.
     .       cflags(ix,iy,CELLFLAG_BOTTOMFACE),
     .       'Boundary number mismatch on bottom face')
          end do
          cflags(ix,iy,CELLFLAG_RIGHTFACE) =
     .     cflags1(ix1+kx1-1,iy1,CELLFLAG_RIGHTFACE)
          do jy1 = 1, ky1-1
            call xertst(
     .       cflags1(ix1+kx1-1,iy1+jy1,CELLFLAG_RIGHTFACE).eq. 
     .       cflags(ix,iy,CELLFLAG_RIGHTFACE),
     .       'Boundary number mismatch on right face')
          end do
          cflags(ix,iy,CELLFLAG_TOPFACE) =
     .     cflags1(ix1,iy1+ky1-1,CELLFLAG_TOPFACE)
          do jx1 = 1, kx1-1
            call xertst(
     .       cflags1(ix1+jx1,iy1+ky1-1,CELLFLAG_TOPFACE).eq. 
     .       cflags(ix,iy,CELLFLAG_TOPFACE),
     .       'Boundary number mismatch on top face')
          end do
cxpb Moving on to connectivity arrays
          if(isInDomain(nx1,ny1,leftix1(ix1,iy1),leftiy1(ix1,iy1))) then
            ic = 0
            do while (leftix1(ix1,iy1).ge.ix1iso(ic+1))
              ic = ic + 1
            enddo
            if (leftix1(ix1,iy1).eq.ix1iso(ic)) then
              leftix(ix,iy)=
     .         floor(1.0_R8*(leftix1(ix1,iy1)-(2*ic-1))*nxeff/nx1eff)+
     .        (2*ic-1)
            else
              leftix(ix,iy)=
     .         floor(1.0_R8*(leftix1(ix1,iy1)-(2*ic))*nxeff/nx1eff)+
     .        (2*ic)
            endif
          else
            leftix(ix,iy)=-2
          endif
          if(rightix1(ix1+kx1-1,iy1).ne.nx1+1) then
            ic = 0
            do while (rightix1(ix1+kx1-1,iy1).ge.ix1iso(ic+1))
              ic = ic + 1
            enddo
            rightix(ix,iy)=
     .         floor(1.0_R8*(rightix1(ix1+kx1-1,iy1)-2*ic)*nxeff/nx1eff)
     .         +2*ic
          else
            rightix(ix,iy)=nx+1
          endif
          ic = 0
          do while (topix1(ix1+kx1-1,iy1+ky1-1).ge.ix1iso(ic+1))
            ic = ic + 1
          enddo
          if(topix1(ix1+kx1-1,iy1+ky1-1).eq.ix1iso(ic)) then
            topix(ix,iy) = floor(1.0_R8*
     .       (topix1(ix1+kx1-1,iy1+ky1-1) - (2*ic-1)) *nxeff/nx1eff) + 
     .       (2*ic-1)
          else
            topix(ix,iy) = floor(1.0_R8*
     .       (topix1(ix1+kx1-1,iy1+ky1-1) - (2*ic)) *nxeff/nx1eff) + 
     .       (2*ic)
          endif
          ic = 0
          do while (bottomix1(ix1,iy1).ge.ix1iso(ic+1))
            ic = ic + 1
          enddo
          if (bottomix1(ix1,iy1).eq.ix1iso(ic)) then
            bottomix(ix,iy) = floor(1.0_R8*
     .       (bottomix1(ix1,iy1) - (2*ic-1)) *nxeff/nx1eff) + (2*ic-1)
          else
            bottomix(ix,iy) = floor(1.0_R8*
     .       (bottomix1(ix1,iy1) - (2*ic)) *nxeff/nx1eff) + (2*ic)
          endif
          ic = 0
          do while (leftiy1(ix1,iy1).ge.iy1iso(ic+1))
            ic = ic+1
          enddo
          if (leftiy1(ix1,iy1).eq.iy1iso(ic+1)) then
            leftiy(ix,iy) =
     .       floor(1.0_R8*(leftiy1(ix1,iy1)-(ic-1))*nyeff/ny1eff)+(ic-1)
          else
            leftiy(ix,iy) = 
     .       floor(1.0_R8*(leftiy1(ix1,iy1)-ic)*nyeff/ny1eff)+ic
          endif
          ic = 0
          do while (rightiy1(ix1+kx1-1,iy1).ge.iy1iso(ic+1))
            ic = ic+1
          enddo
          if (rightiy1(ix1+kx1-1,iy1).eq.iy1iso(ic+1)) then
            rightiy(ix,iy)=
     .       floor(1.0_R8*(rightiy1(ix1+kx1-1,iy1)-(ic-1))*nyeff/ny1eff)
     .      +(ic-1)
          else
            rightiy(ix,iy)=
     .       floor(1.0_R8*(rightiy1(ix1+kx1-1,iy1)-ic)*nyeff/ny1eff)+ic
          endif
          if(topiy1(ix1+kx1-1,iy1+ky1-1).ne.ny1+1) then
            ic = 0
            do while (topiy1(ix1+kx1-1,iy1+ky1-1).ge.iy1iso(ic+1))
              ic = ic + 1
            enddo
            if (topiy1(ix1+kx1-1,iy1+ky1-1).eq.iy1iso(ic+1)) then
              topiy(ix,iy) = floor(1.0_R8*
     .         (topiy1(ix1+kx1-1,iy1+ky1-1)-(ic-1))*nyeff/ny1eff)+(ic-1)
            else
              topiy(ix,iy) = floor(1.0_R8*
     .         (topiy1(ix1+kx1-1,iy1+ky1-1)-ic)*nyeff/ny1eff) + ic
            endif
          else
            topiy(ix,iy)=ny+1
          endif
          if(bottomiy1(ix1,iy1).ne.-2) then
            ic = 0
            do while (bottomiy1(ix1,iy1).ge.iy1iso(ic+1)) 
              ic = ic + 1
            enddo
            if (bottomiy1(ix1,iy1).eq.iy1iso(ic+1)) then
              bottomiy(ix,iy)=
     .         floor(1.0_R8*(bottomiy1(ix1,iy1)-(ic-1))*nyeff/ny1eff)
     .        +(ic-1)
            else
              bottomiy(ix,iy)=
     .         floor(1.0_R8*(bottomiy1(ix1,iy1)-ic)*nyeff/ny1eff)+ic
            endif
          else
            bottomiy(ix,iy)=-2
          endif
cdpc-end
*     ..project hx to cell (ix,iy)
          t0 = 0.0_R8
          do jy1 = 0, ky1-1
            t1 = 0.0_R8
            do jx1 = 0, kx1-1
              t1 = t1+hxx(ix1+jx1,iy1+jy1)**2/vol1(ix1+jx1,iy1+jy1)
            enddo
            t0 = t0+1/t1
          enddo
          hx(ix,iy) = sqrt(vol(ix,iy)/t0)
          cell_width(ix,iy) = hx(ix,iy)
*     ..project hy to cell (ix,iy)
          t0 = 0.0_R8
          do jx1 = 0, kx1-1
            t1 = 0.0_R8
            do jy1 = 0, ky1-1
              t1 = t1+(hyy(ix1+jx1,iy1+jy1)*qz1(ix1+jx1,iy1+jy1,1))**2/
     &	       vol1(ix1+jx1,iy1+jy1)
            enddo
            t0 = t0+1.0_R8/t1
          enddo
          hy(ix,iy) = sqrt(vol(ix,iy)/t0)/qz(ix,iy,1)
          cell_height(ix,iy) = hy(ix,iy)*qz(ix,iy,1)
*     ..update cut counters
          if (ix.eq.ixiso(isox+1)-1) then
            isox = isox+1
          endif
        enddo
        if (iy.eq.iyiso(isoy+1)-1) then
          isoy = isoy+1
        endif
      enddo

c dpc 2000.06.02 begin   xpb 11.08.00
      do iy = -1, ny
        if(leftix(-1,iy).eq.-1) leftix(-1,iy)=-2
        if(rightix(nx,iy).eq.nx) rightix(nx,iy)=nx+1
      enddo
      do ix = -1, nx
        if(bottomiy(ix,-1).eq.-1) bottomiy(ix,-1)=-2
        if(topiy(ix,ny).eq.ny) topiy(ix,ny)=ny+1
      enddo
c dpc 2000.06.02 end    xpb 11.08.00

      do ic=1,nncutmax
        bottomcut(ic)=ny+1
        topcut(ic)=-2
        rightcut(ic)=-2
        leftcut(ic)=nx+1
      enddo
      do ic = 1, nncut
        isox = 0
        do while (leftcut1(ic).ge.ix1iso(isox+1))
          isox = isox + 1
        enddo
        if(leftcut1(ic).ne.nx1+1) then
          if (leftcut1(ic).eq.ix1iso(isox+1)) then
            leftcut(ic) = 
     .       floor(1.0_R8*(leftcut1(ic)-(2*isox-1))*nxeff/nx1eff) 
     .      + (2*isox-1)
          else
            leftcut(ic) = 
     .       floor(1.0_R8*(leftcut1(ic)-2*isox)*nxeff/nx1eff) + (2*isox)
          endif
        endif
        isox = 0
        do while (rightcut1(ic).ge.ix1iso(isox+1))
          isox = isox + 1
        enddo
        if(rightcut1(ic).ne.-2) then
          if (rightcut1(ic).eq.ix1iso(isox+1)) then
            rightcut(ic) = 
     .       floor(1.0_R8*(rightcut1(ic)-(2*isox-1))*nxeff/nx1eff)+
     .       (2*isox-1)
          else
            rightcut(ic) = 
     .       floor(1.0_R8*(rightcut1(ic)-2*isox)*nxeff/nx1eff)+(2*isox)
          endif
        endif
        isoy = 0
        do while (topcut1(ic).ge.iy1iso(isoy+1))
          isoy = isoy + 1
        enddo
        if(topcut1(ic).ne.-2) then
          if (topcut1(ic).eq.iy1iso(isoy+1)) then
            topcut(ic)=floor(1.0_R8*(topcut1(ic)-(isoy-1))*nyeff/ny1eff)
     .       +(isoy-1)
          else
            topcut(ic)=floor(1.0_R8*(topcut1(ic)-isoy)*nyeff/ny1eff)+
     .       isoy
          endif
        endif
        isoy = 0
        do while (bottomcut1(ic).ge.iy1iso(isoy+1))
          isoy = isoy + 1
        enddo
        if(bottomcut1(ic).ne.ny1+1) then
          if (bottomcut1(ic).eq.iy1iso(isoy+1)) then
            bottomcut(ic)=floor(1.0_R8*
     .       (bottomcut1(ic)-(isoy-1))*nyeff/ny1eff)+(isoy-1)
          else
            bottomcut(ic)=floor(1.0_R8*
     .       (bottomcut1(ic)-isoy)*nyeff/ny1eff)+isoy
          endif
        endif
      enddo

      call ipgeti ('b2agfs_nncut', nncut)
      call xertst (0.le.nncut.and.nncut.le.nncutmax,
     . 'faulty argument nncut')
      if(nncut.ge.1) then
        call ipgeti ('b2agfs_leftcut', leftcut(1))
        call xertst (-1.le.leftcut(1).and.leftcut(1).lt.nx,
     .     'faulty argument leftcut')
        call ipgeti ('b2agfs_rightcut', rightcut(1))
        call xertst (leftcut(1).lt.rightcut(1).and.rightcut(1).le.nx,
     .   'faulty argument rightcut')
        call ipgeti ('b2agfs_topcut', topcut(1))
        call xertst (-2.le.topcut(1).and.topcut(1).le.ny,
     .   'faulty argument topcut')
        call ipgeti ('b2agfs_bottomcut', bottomcut(1))
        call xertst (bottomcut(1).le.topcut(1) .or. bottomcut(1).eq.-1,
     .   'faulty argument bottomcut')
      endif
      if(nncut.ge.2) then
        call ipgeti ('b2agfs_leftcut2', leftcut(2))
        call xertst (leftcut(1).lt.leftcut(2).and.
     .               leftcut(2).lt.rightcut(1),
     .              'faulty argument leftcut2')
        call ipgeti ('b2agfs_rightcut2', rightcut(2))
        call xertst (leftcut(2).lt.rightcut(2).and.
     .               rightcut(2).lt.rightcut(1),
     .              'faulty argument rightcut2')
        call ipgeti ('b2agfs_topcut2', topcut(2))
        call xertst (-2.le.topcut(2).and.topcut(2).le.ny,
     .              'faulty argument topcut2')
        call ipgeti ('b2agfs_bottomcut2', bottomcut(2))
        call xertst (bottomcut(2).le.topcut(2),
     .              'faulty argument bottomcut2')
      endif
cxpb  extra test since cases beyond double-null are not yet coded
      call xertst (nncut.le.2, 'faulty argument nncut')

      if(nncut.eq.0) write(*,'(a)') 'No cuts found'
      if(nncut.ge.1) write(*,'(a,4i4)') 
     > 'Using leftcut , rightcut , topcut , bottomcut  = ',
     1 leftcut(1), rightcut(1), topcut(1), bottomcut(1)
      if(nncut.ge.2) write(*,'(a,4i4)') 
     > 'Using leftcut2, rightcut2, topcut2, bottomcut2 = ',
     1 leftcut(2), rightcut(2), topcut(2), bottomcut(2)

      call init_region(nx,ny,nncut,nncutmax,
     1 leftcut,rightcut,topcut,bottomcut,
     1 leftix,rightix,rightiy,topix,topiy,bottomiy,
     2 region,nnreg,resignore,crx,cry,periodic_bc)

cxpb  extra code to allow for overwriting the pbs array with the 
cxpb  working mesh values
      if (pbs_from_basis_mesh.eq.0 .or. redef_gmtry.eq.1) then
*  ..interpolate magnetic field to cell faces
c!!!This interpolation scheme does not guarantee that bb(3)=sqrt(bb(0:2)**2)!!!
        call intfaceh1(nx,ny,2,leftix,leftiy,vol,bb,wbbl)
        call intfacer1(nx,ny,2,rightix,rightiy,vol,bb,wbbr)
        call intfacev1(nx,ny,2,bottomix,bottomiy,vol,bb,wbbv)
*   ..extrapolate magnetic field to edges
        do iy = -1, ny
          do ix = -1, nx
            if (.not.isInDomain(nx,ny,leftix(ix,iy),leftiy(ix,iy))) 
     >       wbbl(ix,iy,0:2) =
     &        bb(ix,iy,0:2) - (wbbr(ix,iy,0:2)-bb(ix,iy,0:2))
            if (.not.isInDomain(nx,ny,rightix(ix,iy),rightiy(ix,iy))) 
     >        wbbr(ix,iy,0:2) =
     &         bb(ix,iy,0:2) - (wbbl(ix,iy,0:2)-bb(ix,iy,0:2))
            if (.not.isInDomain(nx,ny,bottomix(ix,iy),bottomiy(ix,iy)))
     >         wbbv(ix,iy,0:2) = bb(ix,iy,0:2) -
     &       (wbbv(topix(ix,iy),topiy(ix,iy),0:2)-bb(ix,iy,0:2))
          enddo
        enddo
*   .. correct interpolation of magnetic field for diverted island case
*      as B is not a periodic quantity
        if(periodic_bc.eq.1 .and. nncut.ge.1) then
          do iy = -1, ny
            do ix = -1, nx
              if(leftix(ix,iy).ne.ix-1) then
                if(region(ix,iy,0).eq.1 .or. region(ix,iy,0).eq.4) then
                  wbbl(ix,iy,0)=bb(ix,iy,0)-(wbbr(ix,iy,0)-bb(ix,iy,0))
                  wbbl(ix,iy,1)=bb(ix,iy,1)-(wbbr(ix,iy,1)-bb(ix,iy,1))
                  wbbl(ix,iy,2)=bb(ix,iy,2)-(wbbr(ix,iy,2)-bb(ix,iy,2))
                endif
              else if(rightix(ix,iy).ne.ix+1) then
                if(region(ix,iy,0).eq.1 .or. region(ix,iy,0).eq.3) then
                  wbbr(ix,iy,0)=bb(ix,iy,0)-(wbbl(ix,iy,0)-bb(ix,iy,0))
                  wbbr(ix,iy,1)=bb(ix,iy,1)-(wbbl(ix,iy,1)-bb(ix,iy,1))
                  wbbr(ix,iy,2)=bb(ix,iy,2)-(wbbl(ix,iy,2)-bb(ix,iy,2))
                endif
              endif
            enddo
          enddo
        endif
*   ..compute extrapolated wbb.(:,:,3) components
        wbbl(:,:,3)=sqrt(wbbl(:,:,0)**2+wbbl(:,:,1)**2+wbbl(:,:,2)**2)
        wbbr(:,:,3)=sqrt(wbbr(:,:,0)**2+wbbr(:,:,1)**2+wbbr(:,:,2)**2)
        wbbv(:,:,3)=sqrt(wbbv(:,:,0)**2+wbbv(:,:,1)**2+wbbv(:,:,2)**2)
*   ..compute pbs (magnetic field pitch * surface area) at cell faces
        call b2xgbs (nx, ny, wbbl, wbbv, gs, qc, crx, cry, pbs)      
      endif

* ..return
      call subend ()
      return

*-----------------------------------------------------------------------
*.end b2agmx

      end subroutine b2agmx
