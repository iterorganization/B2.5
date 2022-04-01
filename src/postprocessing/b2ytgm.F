*-----------------------------------------------------------------------
*.specification

      subroutine b2ytgm (nget, nx1, ny1, nncut, nnreg)
      use b2mod_types
      use b2mod_geo2
      use b2mod_math
      use b2mod_version
      use b2mod_b2cmgs
      use b2mod_subsys
      implicit none
#include "DIMENSIONS.F"
*   ..input arguments (unchanged on exit)
      integer nget, nx1, ny1
*   ..output arguments (unspecified on entry)
      integer nncut, nnreg(0:2)
*   ..common blocks
*-----------------------------------------------------------------------
*.documentation
*
*  1. purpose
*
*     B2YTGM reads from an un*formatted file the new geometry and magnetic
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
*     nx1, ny1 - integer, input.
*     See the calling routine for description.
*
*     crx1, cry1 - (-1:nx1,-1:ny1,0:3) real*8 array, output.
*     For (ix,iy) in (-1:nx1,-1:ny1), crx1(ix,iy,0:3) specify the
*     horizontal coordinates and cry1(ix,iy,0:3) specify the vertical
*     coordinates (in the physical plane) of the four corners of the
*     (ix,iy) cell, mapped as follows:
*     (,,0): lower left corner,  (,,1): lower right corner,
*     (,,2): upper left corner,  (,,3): upper right corner.
*     One expects the following identities to hold away from the
*     boundaries:
*     crx1(ix,iy,0)=crx1(ix-1,iy,1)=crx1(ix,iy-1,2)=crx1(ix-1,iy-1,3),
*     cry1(ix,iy,0)=cry1(ix-1,iy,1)=cry1(ix,iy-1,2)=cry1(ix-1,iy-1,3),
*     However, there may be cuts in the geometry, in which case these
*     identities need not hold.
*
*     fpsi1 - (-1:nx1,-1:ny1,0:3) real*8 array, output.
*     For (ix,iy,i) in (-1:nx1,-1:ny1,0:3), fpsi1(ix,iy,i) specifies the
*     potential for the magnetic field components (bx,by) at the
*     (ix,iy,i) position.
*     (In the case of straight symmetry, bx=dpsi/dy and by=-dpsi/dx.
*     In the case of toroidal or cylindrical symmetry,
*       bx=(dpsi/dy)/(2*pi*R) and by=-(dpsi/dx)/(2*pi*R),
*     where R is the distance to the axis of symmetry.)
*
*     ffbz1 - (-1:nx1,-1:ny1,0:3) real*8 array, output.
*     For (ix,iy,i) in (-1:nx1,-1:ny1,0:3), ffbz1(ix,iy,i) specifies the
*     field function for the magnetic field component bz at the
*     (ix,iy,i) position.
*     (In the case of straight symmetry, bz=ffbz. In the case of
*     toroidal or cylindrical symmetry, bz=ffbz/(2*pi*R), where R is
*     the distance to the axis of symmetry.)
*
*     bb1 - (-1:nx1,-1:ny1,0:3) real*8 array, output.
*     For (ix,iy) in (-1:nx1,-1:ny1), bb1(ix,iy,0:3) specifies the
*     magnetic field at the center of the (ix,iy) cell.
*     bb1(ix,iy,0:2) are the (x,y,z)-components and bb1(ix,iy,3) is
*     the absolute magnetic field strength.
*     (bb1(,,3)=sqrt(bb1(,,0)**2+bb1(,,1)**2+bb1(,,2)**2.)
*     It will hold that 0.lt.bb1(,,3).
*
*     vol1 - (-1:nx1,-1:ny1) real*8 array, output.
*     For (ix,iy) in (-1:nx1,-1:ny1), vol1(ix,iy) specifies the volume of
*     the (ix,iy) cell. (In the case of toroidal symmetry this volume
*     includes a factor 2*pi*crx1, while in the case of cylindrical
*     symmetry it includes a factor 2*pi*cry1.)
*     It will hold that 0.lt.vol1(,).
*
*     hxx, hyy - (-1:nx1,-1:ny1) real*8 array, output.
*     For (ix,iy) in (-1:nx1,-1:ny1), hxx(ix,iy) and hyy(ix,iy) specify the
*     x-diameter and the y-diameter, respectively, of the (ix,iy) cell.
*     It will hold that 0.lt.hxx(,) and 0.lt.hyy(,).
*
*     hz1 - (-1:nx1,-1:ny1) real*8 array, output.
*     For (ix,iy) in (-1:nx1,-1:ny1), hz1(ix,iy) specifies the
*     length along the toroidal direction of the (ix,iy) cell.
*     It will hold that 0.lt.hz1(,).
*
*     qz1 - (-1:nx1,-1:ny1,0:1) real*8 array, output.
*     For (ix,iy) in (-1:nx1,-1:ny1), qz1(ix,iy,0) specifies sin(t) and
*     qz1(ix,iy,1) specifies cos(t), where t is the complementary angle
*     between the x-direction and the y-direction on the (ix,iy) cell.
*     t is in the range (-pi/2 to pi/2) and vanishes for an orthogonal
*     geometry. (qz1(ix,iy,0) is the inner product of the x and y unit
*     vectors.)
*     It will hold that -1.lt.qz1(,,0).lt.1 and 0.lt.qz1(,,1).le.1.
*
*     qc1 - (-1:nx,-1:ny) real*8 array, output.
*     For (ix,iy) in (-1:nx,-1:ny), qc1(ix,iy) specifies cos(t), where 
*     t is the complementary angle between the x-direction of the 
*     (ix,iy) cell and its left face.
*     t is in the range (-pi/2 to pi/2) and vanishes for an orthogonal
*     geometry. 
*     It will hold that 0.lt.qc1(,).le.1.
*
*     gs1 - (-1:nx1,-1:ny1,0:2) real*8 array, output.
*     (Let gsx(,)=gs1(,,0), gsy(,)=gs1(,,1) and gsz(,)=gs1(,,2).)
*     For (ix,iy) in (0:nx1,-1:ny1), gsx(ix,iy) specifies the area of the
*     face between the (ix,iy) cell and its left neighbor. For ix.eq.-1,
*     gsx(ix,-1:ny1) will hold 0.
*     For (ix,iy) in (-1:nx1,0:ny1), gsy(ix,iy) specifies the area of the
*     face between the (ix,iy-1) and the (ix,iy) cell. For iy.eq.-1,
*     gsy(-1:nx1,iy) will hold 0.
*     For (ix,iy) in (-1:nx1,-1:ny1), gsz(ix,iy) specifies the area of the
*     (ix,iy) cell normal to the third (ignorable) coordinate.
*
*     pbs1 - (-1:nx1,-1:ny1,0:1) real*8 array.
*     (Let pbsx(,)=pbs1(,,0) and pbsy(,)=pbs1(,,1).)
*     For (ix,iy) in (-1:nx1,-1:ny1), pbsx(ix,iy) specifies the product
*     (bx1/bb1)*sx1 and pbsy(ix,iy) specifies (by1/bb1)*sy1, or (magnetic 
*     field pitch)*(surface area), namely the parallel contact area, 
*     on the (ix,iy) cell left and bottom faces respectively.  
*     See routine b2xgbs1 for the form of interpolation employed.
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
      integer n2, ix, iy, idum(0:9), periodic_bc
      integer iFace, offset
*   ..procedures
      external xertst, cfrure, cfruin
      external b2xvsg, b2xvgm1

*-----------------------------------------------------------------------
*.computation

* ..preliminaries
*   ..subprogram start-up calls
      call subini ('b2ytgm')
*   ..set n2
      n2 = (nx1+2)*(ny1+2)
*   ..test nget, nx, ny
      call xertst (1.le.nget, 'faulty argument nget')
      call xertst (0.le.nx1.and.0.le.ny1, 'faulty argument nx, ny')

* ..read geometry and field
*   ..read coordinates
      call cfrure (nget, n2*4, crx1, 'crx')
      call cfrure (nget, n2*4, cry1, 'cry')
*   ..read flux functions
      call cfrure (nget, n2*4, fpsi1, 'fpsi')
      call cfrure (nget, n2*4, ffbz1, 'ffbz')
*   ..read magnetic field
      call cfrure (nget, n2*4, bb1, 'bb')
*   ..read vol, hx, hy
      call cfrure (nget, n2, vol1, 'vol')
      call cfrure (nget, n2, hxx, 'hx')
      call cfrure (nget, n2, hyy, 'hy')
*   ..read qz, qc, gs
      call cfrure (nget, n2*2, qz1, 'qz')
cxpb
      if(b2fgmtry_version.ge.'03.001.000') then
        call cfrure (nget, n2*2, qc1, 'qc')
        call cfrure (nget, n2*2, qs1, 'qs')
      elseif(b2fgmtry_version.ge.'01.001.022') then
        call cfrure (nget, n2, qc1(-1,-1,0), 'qc')
      else
       do iy = -1, ny1
        do ix = -1, nx1
         qc1(ix,iy,0) = cosang (
     &    (crx1(ix,iy,1)+crx1(ix,iy,3))-(crx1(ix,iy,0)+crx1(ix,iy,2)),
     &    (cry1(ix,iy,1)+cry1(ix,iy,3))-(cry1(ix,iy,0)+cry1(ix,iy,2)),
     &     crx1(ix,iy,2)-crx1(ix,iy,0),  cry1(ix,iy,2)-cry1(ix,iy,0))
         if(qc1(ix,iy,0).le.0.0_R8) then
          write(*,*) 'Problem in b2ytgm: qc1(ix,iy) <= 0',
     1     qc1(ix,iy,0),ix,iy
          write(*,'(a,4f12.6)') 'crx1 = ', crx1(ix,iy,0:3)
          write(*,'(a,4f12.6)') 'cry1 = ', cry1(ix,iy,0:3)
         endif
        enddo
       enddo
      endif
      if(b2fgmtry_version.lt.'03.001.000') then
       qc1(:,:,1) = 1.0_R8
       qs1(:,:,1) = 0.0_R8
       if(b2fgmtry_version.ge.'01.001.030') then
         call cfrure_opt (nget, n2, qc1(-1,-1,1), 'qcb')
       end if
       if (maxval(qc1(:,:,1)).eq.1.0_R8 .and. 
     &     minval(qc1(:,:,1)).eq.1.0_R8) then
        do iy = -1, ny1
         do ix = -1, nx1
          qc1(ix,iy,1) = cosang (
     &      crx1(ix,iy,1)-crx1(ix,iy,0),  cry1(ix,iy,1)-cry1(ix,iy,0),
     &     (crx1(ix,iy,2)+crx1(ix,iy,3))-(crx1(ix,iy,0)+crx1(ix,iy,1)),
     &     (cry1(ix,iy,2)+cry1(ix,iy,3))-(cry1(ix,iy,0)+cry1(ix,iy,1)))
          if(qc1(ix,iy,1).le.0.0_R8) then
           write(*,'(a,1f12.6,2i5)')
     1      'Problem in b2rugm: qc1(ix,iy,1) <= 0',qc1(ix,iy,1),ix,iy
           write(*,'(a,4f12.6)') 'crx1 = ', crx1(ix,iy,0:3)
           write(*,'(a,4f12.6)') 'cry1 = ', cry1(ix,iy,0:3)
          endif
         enddo
        enddo
       endif
       do iy = -1, ny1
        do ix = -1, nx1
         qs1(ix,iy,0) = sinang (
     &    (crx1(ix,iy,1)+crx1(ix,iy,3))-(crx1(ix,iy,0)+crx1(ix,iy,2)),
     &    (cry1(ix,iy,1)+cry1(ix,iy,3))-(cry1(ix,iy,0)+cry1(ix,iy,2)),
     &     crx1(ix,iy,2)-crx1(ix,iy,0),  cry1(ix,iy,2)-cry1(ix,iy,0))
        enddo
       enddo
      endif
cxpb-end
      call cfrure (nget, n2*3, gs1, 'gs')
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
      nncut = idum(0)
      call xertst (0.le.nncut,'faulty input nncut')
      call xertst (nncut.le.DEF_NCUT,
     . 'increase DEF_NCUT in DIMENSIONS.F!')
      call cfruin (nget, nncut, leftcut1, 'leftcut')
      call cfruin (nget, nncut, rightcut1, 'rightcut')
      call cfruin (nget, nncut, topcut1, 'topcut')
      call cfruin (nget, nncut, bottomcut1, 'bottomcut')
      call cfruin (nget, n2, leftix1, 'leftix')
      call cfruin (nget, n2, rightix1, 'rightix')
      call cfruin (nget, n2, topix1, 'topix')
      call cfruin (nget, n2, bottomix1, 'bottomix')
      call cfruin (nget, n2, leftiy1, 'leftiy')
      call cfruin (nget, n2, rightiy1, 'rightiy')
      call cfruin (nget, n2, topiy1, 'topiy')
      call cfruin (nget, n2, bottomiy1, 'bottomiy')
      call cfruin (nget, n2*3, region1, 'region')
      call cfruin (nget, 3, nnreg, 'nnreg')
      call cfruin (nget, n2*2, resignore1, 'resignore')
cdpc-end
cxpb
      if(b2fgmtry_version.ge.'01.001.000') then
        call cfruin (nget, 1, idum, 'periodic_bc')
        periodic_bc=idum(0)
      else
        periodic_bc=0
      endif
cxpb-end
*  ..interpolate magnetic field to cell faces
c!!!This interpolation scheme does not guarantee that bb1(3)=sqrt(bb1(0:2)**2)!!!
      call intfaceh1(nx1,ny1,2,leftix1,leftiy1,vol1,bb1,wbbl1)
      call intfacer1(nx1,ny1,2,rightix1,rightiy1,vol1,bb1,wbbr1)
      call intfacev1(nx1,ny1,2,bottomix1,bottomiy1,vol1,bb1,wbbv1)
      call my_intcor1(nx1,ny1,2,vol1,periodic_bc,nncut,
     1     leftix1,leftiy1,rightix1,rightiy1,
     1     bottomix1,bottomiy1,topix1,topiy1,region1,bb1,wbbc1)
*   ..extrapolate magnetic field to edges
      do iy = -1, ny1
        do ix = -1, nx1
          if (leftix1(ix,iy).eq.-2) wbbl1(ix,iy,0:2) =
     &     bb1(ix,iy,0:2) - (wbbr1(ix,iy,0:2)-bb1(ix,iy,0:2))
          if (rightix1(ix,iy).eq.nx1+1) wbbr1(ix,iy,0:2) =
     &     bb1(ix,iy,0:2) - (wbbl1(ix,iy,0:2)-bb1(ix,iy,0:2))
          if (bottomiy1(ix,iy).eq.-2)
     &     wbbv1(ix,iy,0:2) = bb1(ix,iy,0:2) -
     &     (wbbv1(topix1(ix,iy),topiy1(ix,iy),0:2)-bb1(ix,iy,0:2))
        enddo
      enddo
*   .. correct interpolation of magnetic field for diverted island case
*      as B is not a periodic quantity
      if(periodic_bc.eq.1 .and. nncut.ge.1) then
        do iy = -1, ny1
          do ix = -1, nx1
            if(leftix1(ix,iy).ne.ix-1) then
              if(region1(ix,iy,0).eq.1 .or. region1(ix,iy,0).eq.4) then
                wbbl1(ix,iy,0)=
     .           bb1(ix,iy,0)-(wbbr1(ix,iy,0)-bb1(ix,iy,0))
                wbbl1(ix,iy,1)=
     .           bb1(ix,iy,1)-(wbbr1(ix,iy,1)-bb1(ix,iy,1))
                wbbl1(ix,iy,2)=
     .           bb1(ix,iy,2)-(wbbr1(ix,iy,2)-bb1(ix,iy,2))
              endif
            else if(rightix1(ix,iy).ne.ix+1) then
              if(region1(ix,iy,0).eq.1 .or. region1(ix,iy,0).eq.3) then
                wbbr1(ix,iy,0)=
     .           bb1(ix,iy,0)-(wbbl1(ix,iy,0)-bb1(ix,iy,0))
                wbbr1(ix,iy,1)=
     .           bb1(ix,iy,1)-(wbbl1(ix,iy,1)-bb1(ix,iy,1))
                wbbr1(ix,iy,2)=
     .           bb1(ix,iy,2)-(wbbl1(ix,iy,2)-bb1(ix,iy,2))
              endif
            endif
          enddo
        enddo
      endif
*   ..compute extrapolated wbb.(:,:,3) components
      wbbl1(:,:,3)=
     . sqrt(wbbl1(:,:,0)**2+wbbl1(:,:,1)**2+wbbl1(:,:,2)**2)
      wbbr1(:,:,3)=
     . sqrt(wbbr1(:,:,0)**2+wbbr1(:,:,1)**2+wbbr1(:,:,2)**2)
      wbbv1(:,:,3)=
     . sqrt(wbbv1(:,:,0)**2+wbbv1(:,:,1)**2+wbbv1(:,:,2)**2)
      wbbc1(:,:,3)=
     . sqrt(wbbc1(:,:,0)**2+wbbc1(:,:,1)**2+wbbc1(:,:,2)**2)
*   ..compute or read pbs (magnetic field pitch * surface area) at cell faces
      if(b2fgmtry_version.ge.'01.001.063') then                !srv 15.02.11
        call cfrure (nget, n2*2, pbs1, 'pbs')
      else
        call b2xgbs1 (nx1, ny1,
     .   bottomiy1, wbbl1, wbbv1, gs1, qc1, crx1, cry1, pbs1)
      endif
*   ..read extended grid cell flags
      cflags1 = 0
      if(b2fgmtry_version.ge.'01.001.028') then
        call cfruin_opt (nget, n2*GRID_N_CELLFLAGS, cflags1, 'cflags')
      endif
      if(maxval(cflags1).eq.0) then
        cflags1(1:nx1-2,1:ny1-2,CELLFLAG_TYPE) = GRID_INTERNAL
        cflags1(0,0:ny1-1,CELLFLAG_TYPE) = GRID_BOUNDARY
        cflags1(nx1-1,0:ny1-1,CELLFLAG_TYPE) = GRID_BOUNDARY
        cflags1(0:nx1-1,0,CELLFLAG_TYPE) = GRID_BOUNDARY
        cflags1(0:nx1-1,ny1-1,CELLFLAG_TYPE) = GRID_BOUNDARY
        cflags1(-1,-1:ny1,CELLFLAG_TYPE) = GRID_GUARD
        cflags1(nx1,-1:ny1,CELLFLAG_TYPE) = GRID_GUARD
        cflags1(-1:nx1,-1,CELLFLAG_TYPE) = GRID_GUARD
        cflags1(-1:nx1,ny1,CELLFLAG_TYPE) = GRID_GUARD
        if (nncut.ge.2) then
          iy = ny1/2
          do ix = 1, nx1-2
            if (rightix1(ix,iy).eq.nx1+1) then
              cflags1(ix-1,1:ny1-2,CELLFLAG_TYPE) = GRID_BOUNDARY
              cflags1(ix,0:ny1-1,CELLFLAG_TYPE) = GRID_GUARD
            end if
            if (leftix1(ix,iy).eq.-2) then
              cflags1(ix,0:ny1-1,CELLFLAG_TYPE) = GRID_GUARD
              cflags1(ix+1,1:ny1-2,CELLFLAG_TYPE) = GRID_BOUNDARY
            end if
          end do
        end if
        do iy=-1,ny1
          do ix=-1,nx1
            if (cflags1(ix, iy, CELLFLAG_TYPE) /= GRID_BOUNDARY) cycle
            if (cflags1(leftix1(ix,iy),leftiy1(ix,iy),CELLFLAG_TYPE) ==
     .          GRID_GUARD) then
              iFace = CELLFLAG_LEFTFACE
              offset = (iFace - CELLFLAG_LEFTFACE + 1) * (-10)
              cflags1(ix, iy, iFace) = offset - region1(ix, iy, 0)
              cflags1(leftix1(ix,iy),
     .                leftiy1(ix,iy),CELLFLAG_RIGHTFACE) =
     .         cflags1(ix, iy, iFace)
            end if
            if (cflags1(rightix1(ix,iy),
     .                  rightiy1(ix,iy),CELLFLAG_TYPE) ==
     .          GRID_GUARD) then
              iFace = CELLFLAG_RIGHTFACE
              offset = (iFace - CELLFLAG_LEFTFACE + 1) * (-10)
              cflags1(ix, iy, iFace) = offset - region1(ix, iy, 0)
              cflags1(rightix1(ix,iy),
     .                rightiy1(ix,iy),CELLFLAG_LEFTFACE) =
     .         cflags1(ix, iy, iFace)
            end if
            if (cflags1(topix1(ix,iy),topiy1(ix,iy),CELLFLAG_TYPE) ==
     .          GRID_GUARD) then
              iFace = CELLFLAG_TOPFACE
              offset = (iFace - CELLFLAG_LEFTFACE + 1) * (-10)
              cflags1(ix, iy, iFace) = offset - region1(ix, iy, 0)
              cflags1(topix1(ix,iy),topiy1(ix,iy),CELLFLAG_BOTTOMFACE) =
     .         cflags1(ix, iy, iFace)
            end if
            if (cflags1(bottomix1(ix,iy),
     .                  bottomiy1(ix,iy),CELLFLAG_TYPE) ==
     .          GRID_GUARD) then
              iFace = CELLFLAG_BOTTOMFACE
              offset = (iFace - CELLFLAG_LEFTFACE + 1) * (-10)
              cflags1(ix, iy, iFace) = offset - region1(ix, iy, 0)
              cflags1(bottomix1(ix,iy),
     .                bottomiy1(ix,iy),CELLFLAG_TOPFACE) =
     .         cflags1(ix, iy, iFace)
            end if
          end do
        end do
      end if
      cell_width1 = hxx
      cell_height1(:,:) = hyy(:,:)*qz1(:,:,1)
      if(b2fgmtry_version.ge.'01.001.033') then
        call cfrure_opt (nget, n2, cell_width1, 'cell_width')
        call cfrure_opt (nget, n2, cell_height1, 'cell_height')
      endif
*   ..test sign of bb
      call b2xvsg (n2, bb1(-1,-1,3), 1, 'bb3', '.gt.')
      call b2xvsg (n2, wbbl1(-1,-1,3), 1, 'wbbl3', '.gt.')
      call b2xvsg (n2, wbbr1(-1,-1,3), 1, 'wbbr3', '.gt.')
      call b2xvsg (n2, wbbv1(-1,-1,3), 1, 'wbbv3', '.gt.')
      call b2xvsg ((nx1+3)*(ny1+3), wbbc1(-1,-1,3), 1, 'wbbc3', '.gt.')
*   ..test geometry
      call xertst(maxval(region1(-1:nx1,-1:ny1,1)).gt.0 .or.
     &            minval(region1(-1:nx1,-1:ny1,2)).lt.0,
     $            'X-regions not set, re-run b2ag!')
      call xertst(maxval(region1(-1:nx1,-1:ny1,2)).gt.0 .or.
     &            minval(region1(-1:nx1,-1:ny1,1)).lt.0,
     $            'Y-regions not set, re-run b2ag!')
      call b2xvgm1(nx1, ny1, vol1, hxx, hyy, qz1, qc1, qs1, gs1, pbs1,
     . cell_width1, cell_height1, bottomiy1)
*   ..create B2GridMap data structure
      call b2CreateMap(nx1,ny1,crx1,cry1,cflags1,
     &     leftix1,leftiy1,rightix1,rightiy1,
     &     topix1,topiy1,bottomix1,bottomiy1,.true.,gmap1,.false.)
      map1Initialized = .true.

* ..return
      call subend ()
      return

*-----------------------------------------------------------------------
*.end b2ytgm

      end subroutine b2ytgm
