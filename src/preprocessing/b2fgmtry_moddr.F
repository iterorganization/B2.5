*-----------------------------------------------------------------------
*.specification

      subroutine b2fgmtry_moddr (ninp, nout)
      use b2mod_types
      use b2mod_geo
      use b2mod_geo2
      use b2mod_version
      use b2mod_indirect
      use b2mod_b2cmfs
      use b2mod_b2cmgs
      use b2mod_subsys
      implicit none
*   ..input arguments (unchanged on exit)
      integer ninp(0:1), nout(0:1), nx, ny
*   ..output arguments (unspecified on entry)
*     (none)

*-----------------------------------------------------------------------
*.documentation
*
*  1. purpose
*
*
*
*  3. description (see also routine b2cdca)
*
*
*
*  5. parameters (see also routine b2cdcv)
*
*-----------------------------------------------------------------------
*.declarations

*   ..common blocks
*   ..local variables
      integer idum(0:9), ix, iy, ix1
      character lblgm*120, idout*80
      real (kind=R8) ::
     *   parg(0:99)
*   ..procedures
      external xertst, xerrab, cfwuin, cfwuch, cfwure
      external ipgeti, intfaceh1, intfacer1, intfacev1
      external b2xxid, b2agsi, b2agfs, b2agcf, b2agmt, b2agmx, b2agbb,
     &  b2xgbs1, b2wfgm, b2wugm

      real (kind=R8), allocatable :: br(:,:,:), bz(:,:,:), qcdiff(:,:),
     .           vldiff(:,:), sxdiff(:,:), hxdiff(:,:), hydiff(:,:)
      real (kind=R8) :: pi, s213, r213, s143, r143, t0

*-----------------------------------------------------------------------
*.documentation-internal
*
*     Description of some local variables (including arrays in w(0:*):
*
*     crx, cry - (-1:nx,-1:ny,0:3) real*8 array.
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
*     fpsi - (-1:nx,-1:ny,0:3) real*8 array.
*     For (ix,iy,i) in (-1:nx,-1:ny,0:3), fpsi(ix,iy,i) specifies the
*     potential for the magnetic field components (bx,by) at the
*     (ix,iy,i) position.
*     (In the case of straight symmetry, bx=dpsi/dy and by=-dpsi/dx.
*     In the case of toroidal or cylindrical symmetry,
*       bx=(dpsi/dy)/(2*pi*R) and by=-(dpsi/dx)/(2*pi*R),
*     where R is the distance to the axis of symmetry.)
*
*     ffbz - (-1:nx,-1:ny,0:3) real*8 array.
*     For (ix,iy,i) in (-1:nx,-1:ny,0:3), ffbz(ix,iy,i) specifies the
*     field function for the magnetic field component bz at the
*     (ix,iy,i) position.
*     (In the case of straight symmetry, bz=ffbz. In the case of
*     toroidal or cylindrical symmetry, bz=ffbz/(2*pi*R), where R is
*     the distance to the axis of symmetry.)
*
*     bb - (-1:nx,-1:ny,0:3) real*8 array.
*     For (ix,iy) in (-1:nx,-1:ny), bb(ix,iy,0:3) specifies the
*     magnetic field at the center of the (ix,iy) cell.
*     bb(ix,iy,0:2) are the (x,y,z)-components and bb(ix,iy,3) is
*     the absolute magnetic field strength.
*     (bb(,,3)=sqrt(bb(,,0)**2+bb(,,1)**2+bb(,,2)**2.)
*     It will hold that 0.lt.bb(,,3).
*
*     vol - (-1:nx,-1:ny) real*8 array.
*     For (ix,iy) in (-1:nx,-1:ny), vol(ix,iy) specifies the volume of
*     the (ix,iy) cell. (In the case of toroidal symmetry this volume
*     includes a factor 2*pi*crx, while in the case of cylindrical
*     symmetry it includes a factor 2*pi*cry.)
*     It will hold that 0.lt.vol(,).
*
*     hx, hy - (-1:nx,-1:ny) real*8 array.
*     For (ix,iy) in (-1:nx,-1:ny), hx(ix,iy) and hy(ix,iy) specify the
*     x-diameter and the y-diameter, respectively, of the (ix,iy) cell.
*     It will hold that 0.lt.hx(,) and 0.lt.hy(,).
*
*     qz - (-1:nx,-1:ny,0:1) real*8 array.
*     For (ix,iy) in (-1:nx,-1:ny), qz(ix,iy,0) specifies sin(t) and
*     qz(ix,iy,1) specifies cos(t), where t is the complementary angle
*     between the x-direction and the y-direction on the (ix,iy) cell.
*     t is in the range (-pi/2 to pi/2) and vanishes for an orthogonal
*     geometry. (qz(ix,iy,0) is the inner product of the x and y unit
*     vectors.)
*     It will hold that -1.lt.qz(,,0).lt.1 and 0.lt.qz(,,1).le.1.
*
*     qc - (-1:nx,-1:ny,0:1) real*8 array, output.
*     For (ix,iy) in (-1:nx,-1:ny), qc(ix,iy,0) specifies cos(t), where
*     t is the complementary angle between the x-direction of the
*     (ix,iy) cell and its left face.
*     For (ix,iy) in (-1:nx,-1:ny), qc(ix,iy,1) specifies cos(u), where
*     u is the complementary angle between the y-direction of the
*     (ix,iy) cell and its bottom face.
*     t and u are in the range (-pi/2 to pi/2) and vanish for an orthogonal
*     geometry.
*     It will hold that 0.lt.qc(,,).le.1.
*
*     qs - (-1:nx,-1:ny,0:1) real*8 array, output.
*     For (ix,iy) in (-1:nx,-1:ny), qs(ix,iy,0) specifies sin(t), where
*     t is the complementary angle between the x-direction of the
*     (ix,iy) cell and its left face.
*     For (ix,iy) in (-1:nx,-1:ny), qs(ix,iy,1) specifies sin(u), where
*     u is the complementary angle between the y-direction of the
*     (ix,iy) cell and its bottom face.
*     t and u are in the range (-pi/2 to pi/2) and vanish for an orthogonal
*     geometry.
*     It will hold that -1.lt.qs(,,).lt.1.
*
*     gs - (-1:nx,-1:ny,0:2) real*8 array.
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
*     lblgm - character*120.
*     lblgm identifies the instance of the program employed to prepare
*     the geometry.
*
*     idout - character*80.
*     idout specifies the desired selection of output segments for
*     printed and graphical output. Associated with each segment is a
*     four-character word, which is identified on the printed output.
*     (The same word may be associated with several segments.)
*     A particular segment will be produced only if its associated word
*     has a match, as defined by the routine strmas, in idout.
*
*     isymm - integer.
*     isymm specifies the type of symmetry of the geometry: isymm.eq.0
*     implies a slab geometry, isymm.eq.1.or.isymm.eq.2 imply rotational 
*     symmetry about the crx=0 axis, and isymm.eq.3.or.isymm.eq.4 
*     indicate rotational symmetry about the cry=0 axis. 
*     Other values are not allowed. 
*     isymm.eq.1.or.isymm.eq.3 mean a toroidal geometry case, i.e. 
*     with a toroidal (out-of-plane) magnetic field component.
*     isymm.eq.2.or.isymm.eq.4 mean a cylindrical geometry case, i.e. 
*     no toroidal magnetic field component.
*
*     parg - (0:99) real*8 array.
*     parg(0:*) contains a set of numerical parameters that specify
*     details of the geometry. The interpretation of these parameters
*     is problem-dependent; see the code.
*
*     periodic_bc - integer
*     periodic_bc tells whether the cuts correspond to periodic boundary
*     conditions (for example in the case of limited plasmas or
*     a diverted stellarator island in a toroidally averaged geometry). 
*     periodic_bc.eq.0 is the default.
*     periodic_bc.eq.1 indicates that periodic boundary conditions are
*     being used.
*
*-----------------------------------------------------------------------
*.computation

* ..preliminaries
*   ..subprogram start-up calls
      call subini ('b2fgmtry_moddr')
*   ..test ninp, nout
      call xertst (1.le.ninp(0).and.1.le.ninp(1).and.
     &             1.le.nout(0).and.1.le.nout(1),
     &  'faulty argument ninp, nout')
*   ..set output control
      idout = '....'

* ..read the geometry
* ..obtain version number
      call cfverr (ninp(1), b2fgmtry_version)
*   ..obtain nx, ny from geometry file
      call cfruin (ninp(1), 2, idum, 'nx,ny')
      nx = idum(0)
      ny = idum(1)
      call xertst (0.le.nx.and.0.le.ny,
     &  'faulty input nx, ny from geometry file')
*   ..read geometry and metric
      call alloc_b2mod_geo(nx,ny)
      call alloc_b2mod_indirect(nx,ny,nncutmax)
      call cfruch (ninp(1), 120, lblgm, 'label')
      call cfruin (ninp(1), 1, idum, 'isymm')
      isymm=idum(0)
      call b2rugm (ninp(1), nx, ny, crx, cry, fpsi, ffbz,
     &  bb, vol, hx, hy, hz, qz, qc, qs, gs, pbs,
     &  hyl, hxb, hzl, hzb, wbbl, wbbr, wbbv, wbbc, gmap)
      call cfrure (ninp(1), 100, parg, 'parg')

* ..set constants
      pi = 4.d0*atan(1.d0)

* ..calculate temporary br/bz from crx/cry
      allocate (br(-1:nx, -1:ny, 4))
      allocate (bz(-1:nx, -1:ny, 4))

      do ix=-1,nx
        do iy=-1,ny
          br(ix,iy,1) = crx(ix,iy,1)
          br(ix,iy,2) = crx(ix,iy,3)
          br(ix,iy,3) = crx(ix,iy,2)
          br(ix,iy,4) = crx(ix,iy,0)

          bz(ix,iy,1) = cry(ix,iy,1)
          bz(ix,iy,2) = cry(ix,iy,3)
          bz(ix,iy,3) = cry(ix,iy,2)
          bz(ix,iy,4) = cry(ix,iy,0)
        enddo
      enddo
       

* .. allocate some hel arrays
      allocate(qcdiff(-1:nx, -1:ny))
      allocate(sxdiff(-1:nx, -1:ny))
      allocate(vldiff(-1:nx, -1:ny))
      allocate(hxdiff(-1:nx, -1:ny))
      allocate(hydiff(-1:nx, -1:ny))
      do ix=-1,nx
       do iy=-1,ny
         qcdiff(ix,iy) = qc(ix,iy,0)
         sxdiff(ix,iy) = qc(ix,iy,0)*gs(ix,iy,0)
         vldiff(ix,iy) = vol(ix,iy)
         hxdiff(ix,iy) = hx(ix,iy)
         hydiff(ix,iy) = hy(ix,iy)
       enddo
      enddo

* ..modify hx,hy
      do ix = -1, nx
        do iy = -1, ny
         hx(ix,iy) = sqrt( (bz(ix,iy,1) + bz(ix,iy,2) - 
     -                      bz(ix,iy,4) - bz(ix,iy,3))**2 
     +                   + (br(ix,iy,1) + br(ix,iy,2) -
     -                      br(ix,iy,4) - br(ix,iy,3))**2 )
         hx(ix,iy) = 0.5*hx(ix,iy)
         hy(ix,iy) = sqrt( (bz(ix,iy,1) + bz(ix,iy,4) - 
     -                      bz(ix,iy,2) - bz(ix,iy,3))**2 
     +                   + (br(ix,iy,1) + br(ix,iy,4) - 
     -                      br(ix,iy,2) - br(ix,iy,3))**2 )
         hy(ix,iy) = 0.5*hy(ix,iy)
        end do
      end do

* ..modify vol
      do ix = -1, nx
        do iy = -1, ny
         S213=br(ix,iy,2)*(bz(ix,iy,1)-bz(ix,iy,3))+
     +       br(ix,iy,1)*(bz(ix,iy,3)-bz(ix,iy,2))+
     +       br(ix,iy,3)*(bz(ix,iy,2)-bz(ix,iy,1))
         S143=br(ix,iy,1)*(bz(ix,iy,4)-bz(ix,iy,3))+
     +       br(ix,iy,4)*(bz(ix,iy,3)-bz(ix,iy,1))+
     +       br(ix,iy,3)*(bz(ix,iy,1)-bz(ix,iy,4))
         r213=br(ix,iy,2)+br(ix,iy,1)+br(ix,iy,3)
         r143=br(ix,iy,1)+br(ix,iy,4)+br(ix,iy,3)
         vol(ix,iy) = pi*abs((r213*S213+r143*S143))/3. 
        enddo
      enddo

* ..modify gsx (and pbsx) by modifying qc in case of area_fix=3
      do iy=-1,ny
        qc(nx,iy,0) = qz(nx,iy,1)
        do ix=-1,nx-1
          ix1=ix+1
          t0 = (vol(ix1,iy)*qz(ix1,iy,1)+vol(ix,iy)*qz(ix,iy,1)) /
     .         (vol(ix,iy)+vol(ix1,iy))
          qc(ix,iy,0) = t0
        enddo
      enddo


* ..calc diffs and output
      open(unit=4711, file='diff.prt', status='replace',
     .     form='formatted')
      do ix=-1,nx
        do iy=-1,ny
          qcdiff(ix,iy) = qcdiff(ix,iy) / qc(ix,iy,0) - 1.0_R8
          sxdiff(ix,iy) = sxdiff(ix,iy) /
     .                       (qc(ix,iy,0)*gs(ix,iy,0))- 1.0_R8
          vldiff(ix,iy) = vldiff(ix,iy) / vol(ix,iy)- 1.0_R8
          hxdiff(ix,iy) = hxdiff(ix,iy) / hx(ix,iy) - 1.0_R8
          hydiff(ix,iy) = hydiff(ix,iy) / hy(ix,iy) - 1.0_R8
          write(4711,'(2i6,5(1x,e13.6))') ix,iy,vldiff(ix,iy), 
     .           qcdiff(ix,iy), sxdiff(ix,iy), hxdiff(ix,iy),
     .           hydiff(ix,iy)
        enddo
      enddo
      close(4711)

* ..output
*   ..write summary of geometry
      call b2wfgm (nout(0), idout, nx, ny)
*   ..write version number
      call cfverw (nout(1), newversion)
*   ..save dimensions
      idum(0) = nx
      idum(1) = ny
      call cfwuin (nout(1), 2, idum, 'nx,ny')
*   ..save label
      call cfwuch (nout(1), 120, lblgm, 'label')
*   ..save isymm
      idum(0) = isymm
      call cfwuin (nout(1), 1, idum, 'isymm')
*   ..save geometry
      call b2wugm (nout(1), nx, ny)
*   ..save parameters
      call cfwure (nout(1), 100, parg, 'parg')

* ..clean up
      deallocate(br)
      deallocate(bz)
      deallocate(vldiff)
      deallocate(qcdiff)
      deallocate(sxdiff)
      call dealloc_b2mod_geo
      call dealloc_b2mod_indirect
* ..return
      call subend ()
      return

*-----------------------------------------------------------------------
*.end b2fgmtry_moddr

      end subroutine b2fgmtry_moddr
