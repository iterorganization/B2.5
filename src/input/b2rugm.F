*-----------------------------------------------------------------------
*.specification

      subroutine b2rugm (nget, nx, ny, isymm)
      use b2mod_types
      use b2mod_geo
      use b2mod_math
      use b2mod_indirect
      use b2mod_version
      use b2mod_b2cmgs
      use b2mod_grid_mapping
      use b2mod_subsys
      implicit none
*   ..input arguments (unchanged on exit)
      integer nget, nx, ny, isymm
*-----------------------------------------------------------------------
*.documentation
*
*  1. purpose
*
*     B2RUGM reads from an un*formatted file the geometry and magnetic
*     field.
*
*
*  3. description (see also routine b2cdca)
*
*     This routine reads from unit (nget) the specification of a B2
*     geometry and magnetic field. Unit nget must be connected to an
*     un*formatted file.
*     The data consists of the coordinates of the four corners of all
*     the cells in the mesh, the magnetic flux and field functions at
*     those positions, the components of the magnetic field at cell
*     centers, a collection of metric coefficients, and a specification
*     of singular regions in the geometry. See the code to discern the
*     required format of the input file.
*
*
*  5. parameters (see also routine b2cdcv)
*
*     nget - integer, input.
*     nget specifies a unit number for un*formatted input.
*
*     nx, ny - integer, input.
*     See the calling routine for description.
*
*     crx, cry - (-1:nx,-1:ny,0:3) real*8 array, output.
*     For (ix,iy) in (-1:nx,-1:ny), crx(ix,iy,0:3) specify the
*     horizontal coordinates and cry(ix,iy,0:3) specify the vertical
*     coordinates (in the physical plane) of the four corners of the
*     (ix,iy) cell, mapped as follows:
*     (,,0): lower left corner,  (,,1): lower right corner,
*     (,,2): upper left corner,  (,,3): upper right corner.
*     One expects the following identities to hold away from the
*     boundaries:
*     crx(ix,iy,0)=crx(ix-1,iy,1)=crx(ix,iy-1,2)=crx(ix-1,iy-1,3),
*     cry(ix,iy,0)=cry(ix-1,iy,1)=cry(ix,iy-1,2)=cry(ix-1,iy-1,3),
*     However, there may be cuts in the geometry, in which case these
*     identities need not hold.
*
*     fpsi - (-1:nx,-1:ny,0:3,0:2) real*8 array, output.
*     For (ix,iy,i) in (-1:nx,-1:ny,0:3), fpsi(ix,iy,i,0) specifies the
*     potential for the magnetic field components (bx,by) at the
*     (ix,iy,i) position.
*     fpsi(ix,iy,i,1) specifies the (dpsi/dx) derivative.
*     fpsi(ix,iy,i,2) specifies the (dpsi/dy) derivative.
*     (In the case of straight symmetry, bx=dpsi/dy and by=-dpsi/dx.
*     In the case of toroidal or cylindrical symmetry,
*       bx=(dpsi/dy)/(2*pi*R) and by=-(dpsi/dx)/(2*pi*R),
*     where R is the distance to the axis of symmetry.)
*
*     ffbz - (-1:nx,-1:ny,0:3) real*8 array, output.
*     For (ix,iy,i) in (-1:nx,-1:ny,0:3), ffbz(ix,iy,i) specifies the
*     field function for the magnetic field component bz at the
*     (ix,iy,i) position.
*     (In the case of straight symmetry, bz=ffbz. In the case of
*     toroidal or cylindrical symmetry, bz=ffbz/(2*pi*R), where R is
*     the distance to the axis of symmetry.)
*
*     bb - (-1:nx,-1:ny,0:3) real*8 array, output.
*     For (ix,iy) in (-1:nx,-1:ny), bb(ix,iy,0:3) specifies the
*     magnetic field at the center of the (ix,iy) cell.
*     bb(ix,iy,0:2) are the (x,y,z)-components and bb(ix,iy,3) is
*     the absolute magnetic field strength.
*     (bb(,,3)=sqrt(bb(,,0)**2+bb(,,1)**2+bb(,,2)**2.)
*     It will hold that 0.lt.bb(,,3).
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
*     (parameters in common:)
*
*     /b2cmgs/ - output.
*     /b2cmgs/ specifies the singular regions in the geometry.
*      nlrmx - integer constant.
*     nlrmx is an upper bound for nlreg, and is used in dimension
*     specifications.
*      nlreg - integer.
*     nlreg specifies the number of singular regions. It will hold that
*     0.le.nlreg.le.nlrmx.
*      nlxlo, nlxhi, nlylo, nlyhi - (0:nlrmx-1) integer array.
*     For 0.le.i.lt.nlreg, nlxlo(i), nlxhi(i), nlylo(i) and nlyhi(i) are
*     the bounds of the i-th singular region. Mesh locations (ix,iy)
*     in (nlxlo(i):nlxhi(i),nlylo(i):nlyhi(i)) make up the region, which
*     will be included in (-1:nx,-1:ny).
*     For nlreg.le.i.lt.nlrmx, nlxlo(i), nlxhi(i), nlylo(i) and nlyhi(i)
*     are unspecified.
*      nlloc - (0:nlrmx-1) integer array.
*     For 0.le.i.lt.4, nlloc(i) specifies the class of quantities for
*     which the i-th region is special, as follows:
*       nlloc(i) = 0: cell centers, indices (-1:nx,-1:ny).
*       nlloc(i) = 1: cell x-faces, indices (0:nx,-1:ny).
*       nlloc(i) = 2: cell y-faces, indices (-1:nx,0:ny).
*       nlloc(i) = 3: cell corners, indices (0:nx,0:ny).
*     For nlreg.le.i.lt.nlrmx, nlloc(i) is unspecified.
*
*
*  8. further comments
*
*     See routine b2cdca for explanation of 'un*formatted file'.
*
*-----------------------------------------------------------------------
*.declarations

*   ..local variables
      integer n2, ix, iy, idum(0:9)
      integer iFace, offset
      real (kind=R8) ::
     *   t0, t1, t2
      integer, save :: ht_style = 1
*   ..procedures
      external xertst, cfrure, cfruin
      external intfaceh1, intfacer1, intfacev1
      external b2xgbs, b2xvsg, b2xvgm

*-----------------------------------------------------------------------
*.computation

* ..preliminaries
*   ..subprogram start-up calls
      call subini ('b2rugm')
*   ..set n2
      n2 = (nx+2)*(ny+2)
*   ..test nget, nx, ny
      call xertst (1.le.nget, 'faulty argument nget')
      call xertst (0.le.nx.and.0.le.ny, 'faulty argument nx, ny')

* ..read geometry and field
*   ..read coordinates
      call cfrure (nget, n2*4, crx, 'crx')
      call cfrure (nget, n2*4, cry, 'cry')
*   ..read flux functions
      if(b2fgmtry_version.ge.'03.002.000') then
       call cfrure (nget, n2*12, fpsi, 'fpsi')
      else
       call cfrure (nget, n2*4, fpsi(-1,-1,0,0), 'fpsi')
       fpsi(:,:,:,1:2) = 0.0_R8
      endif
      call cfrure (nget, n2*4, ffbz, 'ffbz')
*   ..read magnetic field
      call cfrure (nget, n2*4, bb, 'bb')
*   ..read vol, hx, hy
      call cfrure (nget, n2, vol, 'vol')
      call cfrure (nget, n2, hx, 'hx')
      call cfrure (nget, n2, hy, 'hy')
*   ..read qz, qc, gs
      call cfrure (nget, n2*2, qz, 'qz')
cxpb
      if(b2fgmtry_version.ge.'03.001.000') then
       call cfrure (nget, n2*2, qc, 'qc')
       call cfrure (nget, n2*2, qs, 'qs')
      elseif(b2fgmtry_version.ge.'01.001.022') then
       call cfrure (nget, n2, qc(-1,-1,0), 'qc')
      else
       do iy = -1, ny
        do ix = -1, nx
         qc(ix,iy,0) = cosang (
     &    (crx(ix,iy,1)+crx(ix,iy,3))-(crx(ix,iy,0)+crx(ix,iy,2)),
     &    (cry(ix,iy,1)+cry(ix,iy,3))-(cry(ix,iy,0)+cry(ix,iy,2)),
     &     crx(ix,iy,2)-crx(ix,iy,0),  cry(ix,iy,2)-cry(ix,iy,0))
         if(qc(ix,iy,0).le.0.0_R8) then
          write(*,'(a,1f12.6,2i5)')
     1     'Problem in b2rugm: qc(ix,iy) <= 0',qc(ix,iy,0),ix,iy
          write(*,'(a,4f12.6)') 'crx = ', crx(ix,iy,0:3)
          write(*,'(a,4f12.6)') 'cry = ', cry(ix,iy,0:3)
         endif
        enddo
       enddo
      endif
cxpb-end
      call cfrure (nget, n2*3, gs, 'gs')
*   ..read /b2cmgs/
      call cfruin (nget, 1, idum, 'nlreg')
      nlreg = idum(0)
      call xertst (0.le.nlreg.and.nlreg.le.nlrmx, 'faulty input nlreg')
      call cfruin (nget, nlreg, nlxlo, 'nlxlo')
      call cfruin (nget, nlreg, nlxhi, 'nlxhi')
      call cfruin (nget, nlreg, nlylo, 'nlylo')
      call cfruin (nget, nlreg, nlyhi, 'nlyhi')
      call cfruin (nget, nlreg, nlloc, 'nlloc')
cdpc
      call cfruin (nget, 1, idum, 'nncut')
      nncut=idum(0)
      call xertst (0.le.nncut.and.nncut.le.nncutmax,
     . 'faulty input nncut')
      call cfruin (nget, nncut, leftcut, 'leftcut')
      call cfruin (nget, nncut, rightcut, 'rightcut')
      call cfruin (nget, nncut, topcut, 'topcut')
      call cfruin (nget, nncut, bottomcut, 'bottomcut')
      call cfruin (nget, n2, leftix, 'leftix')
      call cfruin (nget, n2, rightix, 'rightix')
      call cfruin (nget, n2, topix, 'topix')
      call cfruin (nget, n2, bottomix, 'bottomix')
      call cfruin (nget, n2, leftiy, 'leftiy')
      call cfruin (nget, n2, rightiy, 'rightiy')
      call cfruin (nget, n2, topiy, 'topiy')
      call cfruin (nget, n2, bottomiy, 'bottomiy')
      call cfruin (nget, n2*3, region, 'region')
      call cfruin (nget, 3, nnreg, 'nnreg')
      call cfruin (nget, n2*2, resignore, 'resignore')
cdpc-end
      if(b2fgmtry_version.lt.'03.001.000') then
*    ..compute qc and qs at south faces
cwdk   For 9-point stencil, compute south-face values of qc (cos) and qs (sin)
cwdk   temporary computation; should go into geometry routines!
cwk:   this interpolation for backwards compatibility in computation of radial conductances etc.
       do iy = -1, ny
        do ix = -1, nx
         if (bottomiy(ix,iy).eq.-2) then
          qc(ix,iy,1)=0.0_R8
         else
          qc(ix,iy,1) =
     .    (hy(ix,iy)*qz(ix,iy,1) + hy(bottomix(ix,iy),bottomiy(ix,iy))*
     .      qz(bottomix(ix,iy),bottomiy(ix,iy),1))/ 
     .     (hy(ix,iy)  + hy(bottomix(ix,iy),bottomiy(ix,iy)))
          qs(ix,iy,1) =
     .     (hy(ix,iy)*qz(ix,iy,0) + hy(bottomix(ix,iy),bottomiy(ix,iy))*
     .      qz(bottomix(ix,iy),bottomiy(ix,iy),0))/ 
     .     (hy(ix,iy)  + hy(bottomix(ix,iy),bottomiy(ix,iy)))
         endif
        enddo
       enddo
cwdk   now ensure sin^2 + cos^2 = 1, cos > 0, and keep sign of sin
       qs(:,:,1) = sign(sqrt(1.0_R8 - qc(:,:,1)**2),qs(:,:,1))
       ! Init qs(:,:,0)
       do iy = -1, ny
        do ix = -1, nx
         qs(ix,iy,0) = sinang (
     &    (crx(ix,iy,1)+crx(ix,iy,3))-(crx(ix,iy,0)+crx(ix,iy,2)),
     &    (cry(ix,iy,1)+cry(ix,iy,3))-(cry(ix,iy,0)+cry(ix,iy,2)),
     &    crx(ix,iy,2)-crx(ix,iy,0), cry(ix,iy,2)-cry(ix,iy,0))
        enddo
       enddo
cwdk   ensure qc(:,:,0)**2 + qs(:,:,0)**2 = 1, in possible event 
cwdk   of improved computation of qc
       qs(:,:,0) = sign(sqrt(1 - qc(:,:,0)**2),qs(:,:,0))
      endif
cxpb
      if(b2fgmtry_version.ge.'01.001.000') then
        call cfruin (nget, 1, idum, 'periodic_bc')
        periodic_bc=idum(0)
      else
        periodic_bc=0
      endif
cxpb-end
      do iy = -1, ny
        do ix = -1, nx
          if (nnreg(0).eq.1) then
            on_closed_surface(ix,iy) = periodic_bc.eq.1
          else if (nnreg(0).eq.2 .or. nnreg(0).eq.4) then
            on_closed_surface(ix,iy) = region(ix,iy,0).eq.1
          else if (nnreg(0).eq.5 .or. nnreg(0).eq.8) then
            on_closed_surface(ix,iy) = mod(region(ix,iy,0),4).eq.1
          else
            call xerrab ('Unknown number of regions!')
          end if
        end do
      end do
*   ..compute or read pbs (magnetic field pitch * surface area) at cell faces
      if(b2fgmtry_version.ge.'01.001.063') then                      !srv 15.02.11
        call cfrure (nget, n2*2, pbs, 'pbs')
      endif
      if(b2fgmtry_version.eq.'03.000.006') then
        call cfruin (nget, 1, idum, 'redef_pbs')
        redef_gmtry=idum(0)
      elseif(b2fgmtry_version.ge.'03.000.007') then
        call cfruin (nget, 1, idum, 'redef_gmtry')
        redef_gmtry=idum(0)
        call cfrure (nget, n2*4, wbbl, 'wbbl')
        if (redef_gmtry.eq.0) write(*,'(a)')
     .   'Warning: wbbl value read from b2fgmtry will be overwritten '//
     .   'because redef_gmtry.eq.0 !'
      else
        redef_gmtry=0
      endif
      if(b2fgmtry_version.ge.'03.002.000') then
        call cfrure (nget, (nx+3)*(ny+3)*4, wbbc, 'wbbc')
      else
*    ..interpolate magnetic field to cell corners
        call my_intcor1(nx,ny,2,vol,periodic_bc,nncut,
     1     leftix,leftiy,rightix,rightiy,
     1     bottomix,bottomiy,topix,topiy,region,bb,wbbc)
        wbbc(:,:,3)=sqrt(wbbc(:,:,0)**2+wbbc(:,:,1)**2+wbbc(:,:,2)**2)
      endif

*   ..read extended grid cell flags
      cflags = 0
      if(b2fgmtry_version.ge.'01.001.028') then
        call cfruin_opt (nget, n2*GRID_N_CELLFLAGS, cflags, 'cflags')
      endif
      if(maxval(cflags).eq.0) then
        cflags(1:nx-2,1:ny-2,CELLFLAG_TYPE) = GRID_INTERNAL
        cflags(0,0:ny-1,CELLFLAG_TYPE) = GRID_BOUNDARY
        cflags(nx-1,0:ny-1,CELLFLAG_TYPE) = GRID_BOUNDARY
        cflags(0:nx-1,0,CELLFLAG_TYPE) = GRID_BOUNDARY
        cflags(0:nx-1,ny-1,CELLFLAG_TYPE) = GRID_BOUNDARY
        cflags(-1,-1:ny,CELLFLAG_TYPE) = GRID_GUARD
        cflags(nx,-1:ny,CELLFLAG_TYPE) = GRID_GUARD
        cflags(-1:nx,-1,CELLFLAG_TYPE) = GRID_GUARD
        cflags(-1:nx,ny,CELLFLAG_TYPE) = GRID_GUARD
        if (nncut.ge.2) then
          iy = ny/2
          do ix = 1, nx-2
            if (rightix(ix,iy).eq.nx+1) then
              cflags(ix-1,0:ny-1,CELLFLAG_TYPE) = GRID_BOUNDARY
              cflags(ix,-1:ny,CELLFLAG_TYPE) = GRID_GUARD
            end if
            if (leftix(ix,iy).eq.-2) then
              cflags(ix,-1:ny,CELLFLAG_TYPE) = GRID_GUARD
              cflags(ix+1,0:ny-1,CELLFLAG_TYPE) = GRID_BOUNDARY
            end if
          end do
        end if
        do iy=-1,ny
          do ix=-1,nx
            if (cflags(ix, iy, CELLFLAG_TYPE) /= GRID_BOUNDARY) cycle
            if (cflags(leftix(ix,iy),leftiy(ix,iy),CELLFLAG_TYPE) ==
     .          GRID_GUARD) then
              iFace = CELLFLAG_LEFTFACE
              offset = (iFace - CELLFLAG_LEFTFACE + 1) * (-10)
              cflags(ix, iy, iFace) = offset - region(ix, iy, 0)
              cflags(leftix(ix,iy),leftiy(ix,iy),CELLFLAG_RIGHTFACE) =
     .         cflags(ix, iy, iFace)
            end if
            if (cflags(rightix(ix,iy),rightiy(ix,iy),CELLFLAG_TYPE) ==
     .          GRID_GUARD) then
              iFace = CELLFLAG_RIGHTFACE
              offset = (iFace - CELLFLAG_LEFTFACE + 1) * (-10)
              cflags(ix, iy, iFace) = offset - region(ix, iy, 0)
              cflags(rightix(ix,iy),rightiy(ix,iy),CELLFLAG_LEFTFACE) =
     .         cflags(ix, iy, iFace)
            end if
            if (cflags(topix(ix,iy),topiy(ix,iy),CELLFLAG_TYPE) ==
     .          GRID_GUARD) then
              iFace = CELLFLAG_TOPFACE
              offset = (iFace - CELLFLAG_LEFTFACE + 1) * (-10)
              cflags(ix, iy, iFace) = offset - region(ix, iy, 0)
              cflags(topix(ix,iy),topiy(ix,iy),CELLFLAG_BOTTOMFACE) =
     .         cflags(ix, iy, iFace)
            end if
            if (cflags(bottomix(ix,iy),bottomiy(ix,iy),CELLFLAG_TYPE) ==
     .          GRID_GUARD) then
              iFace = CELLFLAG_BOTTOMFACE
              offset = (iFace - CELLFLAG_LEFTFACE + 1) * (-10)
              cflags(ix, iy, iFace) = offset - region(ix, iy, 0)
              cflags(bottomix(ix,iy),bottomiy(ix,iy),CELLFLAG_TOPFACE) =
     .         cflags(ix, iy, iFace)
            end if
          end do
        end do
      end if
      cell_width = hx
      cell_height(:,:) = hy(:,:)*qz(:,:,1)
      if(b2fgmtry_version.ge.'01.001.033') then
        call cfrure_opt (nget, n2, cell_width, 'cell_width')
        call cfrure_opt (nget, n2, cell_height, 'cell_height')
      endif
*  ..interpolate magnetic field to cell faces
      call b2agbf (nx, ny, isymm, periodic_bc, nncut,
     & crx, cry, fpsi, ffbz, vol, gs, qz, qc, qs, bb, wbbl, wbbr, wbbv,
     & region, leftix, leftiy, bottomix, bottomiy, rightix, rightiy,
     & topix, topiy, cflags)
*  ..init some arrays for classical grids
*   ..read  additional geometric parameters for extended grids
      if(b2fgmtry_version.ge.'03.002.000') then
        call cfrure(nget, n2*4, hc, 'hc')
        call cfrure(nget, n2*2, ht, 'ht')
        call cfrure(nget, n2*2, qac, 'qac')
        call cfrure(nget, n2*2, qas, 'qas')
        call cfrure(nget, n2*3, ebc, 'ebc')
      else
        hc(:,:,0) = 0.5*hx(:,:)
        hc(:,:,1) = 0.5*hy(:,:)
        hc(:,:,2) = 0.5*hx(:,:)
        hc(:,:,3) = 0.5*hy(:,:)
        qac(:,:,0) =  qc(:,:,0)*sign(1._R8,wbbl(:,:,0))
        qas(:,:,0) = -qs(:,:,0)*sign(1._R8,wbbl(:,:,0))
        qac(:,:,1) = 0.0_R8
        qas(:,:,1) = 1.0_R8*sign(1._R8,wbbv(:,:,0))
        ! Initialize ebc and ht for 'classical' grids
        do iy=-1,ny
          do ix=-1,nx

            ! Unit vectors of line between W/E face centers
            t0 = (crx(ix,iy,1)+crx(ix,iy,3))-(crx(ix,iy,0)+crx(ix,iy,2)) 
            t1 = (cry(ix,iy,1)+cry(ix,iy,3))-(cry(ix,iy,0)+cry(ix,iy,2))
            t2 = norm(t0,t1)

            ! Unit vector of magnetic field in Cartesian X,Y,Z coordinates
            ebc(ix,iy,0) = bb(ix,iy,0)/bb(ix,iy,3)*t0/t2
            ebc(ix,iy,1) = bb(ix,iy,0)/bb(ix,iy,3)*t1/t2
            ebc(ix,iy,2) = bb(ix,iy,2)/bb(ix,iy,3)

            ! Length along faces, ht
            if (leftix(ix,iy).ne.-2) then
              if (ht_style.eq.0) then
                ! Implementation intended for backwards compatibility
                ! testing
                ht(ix,iy,0) =
     &            0.5_R8*(hy(leftix(ix,iy),leftiy(ix,iy))+hy(ix,iy))
              else
                ! Precise computation
                ht(ix,iy,0) = dist(ix,iy,0,2)
              endif
            else
              ht(ix,iy,0) = 0.0_R8
            endif
            if (bottomiy(ix,iy).ne.-2) then
              if (ht_style.eq.0) then
                ! Implementation intended for backwards compatibility
                ! testing
                ht(ix,iy,1) =
     &            0.5_R8*(hx(bottomix(ix,iy),bottomiy(ix,iy))+hx(ix,iy))
              else
                ! Precise computation
                ht(ix,iy,1) = dist(ix,iy,0,1)
              endif
            else
              ht(ix,iy,1) = 0.0_R8
            endif
          enddo
        enddo
      endif
      if(b2fgmtry_version.lt.'01.001.063') then                      !srv 15.02.11
        ! Moved here; b2xgbs now requires cflags and qac
        call b2xgbs (nx, ny, wbbl, wbbv, gs, qac, pbs,
     &  leftix, leftiy, bottomix, bottomiy, cflags)
      endif
*   ..test sign of bb
      call b2xvsg (n2, bb(-1,-1,3), 1, 'bb3', '.ge.')
      call b2xvsg (n2, wbbl(-1,-1,3), 1, 'wbbl3', '.ge.')
      call b2xvsg (n2, wbbr(-1,-1,3), 1, 'wbbr3', '.ge.')
      call b2xvsg (n2, wbbv(-1,-1,3), 1, 'wbbv3', '.ge.')
      call b2xvsg ((nx+3)*(ny+3), wbbc(-1,-1,3), 1, 'wbbc3', '.ge.')
*   ..test geometry
      call xertst(maxval(region(-1:nx,-1:ny,1)).gt.0 .or.
     &            minval(region(-1:nx,-1:ny,2)).lt.0,
     $            'X-regions not set, re-run b2ag!')
      call xertst(maxval(region(-1:nx,-1:ny,2)).gt.0 .or.
     &            minval(region(-1:nx,-1:ny,1)).lt.0,
     $            'Y-regions not set, re-run b2ag!')
      call b2xvgm (nx, ny, vol, hx, hy, qz, qc, qs, gs, pbs,
     &     cell_width, cell_height)
*   ..create B2GridMap data structure
      if (nget.ne.70) then ! skip case when reading b2fgmtry.old from b2yt
        call b2CreateMap(nx, ny, crx, cry, cflags,
     &     leftix, leftiy, rightix, rightiy,
     &     topix, topiy, bottomix, bottomiy, .true., gmap, .false.)
        mapInitialized = .true.
      endif

* ..return
      call subend ()
      return

*-----------------------------------------------------------------------
*.end b2rugm

      end subroutine b2rugm
