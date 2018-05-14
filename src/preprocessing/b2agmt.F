*-----------------------------------------------------------------------
*.specification

      subroutine b2agmt (nx, ny, isymm, 
     & crx, cry, vol, hx, hy, qz, qc, qs, gs,
     & cell_width, cell_height, leftix, bottomiy)
      use b2mod_types
      use b2mod_constants
      implicit none
*   ..input arguments (unchanged on exit)
      integer nx, ny, isymm
      integer leftix(-1:nx,-1:ny), bottomiy(-1:nx,-1:ny) 
      real (kind=R8) ::
     *  crx(-1:nx,-1:ny,0:3), cry(-1:nx,-1:ny,0:3)
*   ..output arguments (unspecified on entry)
      real (kind=R8) ::
     *  vol(-1:nx,-1:ny), hx(-1:nx,-1:ny), hy(-1:nx,-1:ny),
     &  qz(-1:nx,-1:ny,0:1), qc(-1:nx,-1:ny,0:1),
     &  gs(-1:nx,-1:ny,0:2), qs(-1:nx,-1:ny,0:1), 
     &  cell_width(-1:nx,-1:ny), cell_height(-1:nx,-1:ny)

*-----------------------------------------------------------------------
*.documentation
*
*  1. purpose
*
*     B2AGMT computes the metric coefficients given the locations of
*     cell vertices.
*
*
*  3. description (see also routine b2cdca)
*
*     This routine takes as input the coordinates of the four corners
*     of all the cells in the mesh. It computes metric coefficients of
*     the mesh: vol, hx, hy, qz, and gs.
*
*
*  5. parameters (see also routine b2cdcv)
*
*     nx, ny - integer, input.
*     See the calling routine for description.
*
*     isymm - integer, input.
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
*     crx, cry - (-1:nx,-1:ny,0:3) real*8 array, input.
*     See the calling routine for description.
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
*     qc - (-1:nx,-1:ny,0:1) real*8 array, output.
*     For (ix,iy) in (-1:nx,-1:ny), qc(ix,iy) specifies cos(t), where 
*     t is the complementary angle between the x-direction of the 
*     (ix,iy) cell and its left face.
*     t is in the range (-pi/2 to pi/2) and vanishes for an orthogonal
*     geometry. 
*     It will hold that 0.lt.qc(,,).le.1.
*
*     qs - (-1:nx,-1:ny,0:1) real*8 array, output.
*     For (ix,iy) in (-1:nx,-1:ny), qs(ix,iy) specifies sin(t), where 
*     t is the complementary angle between the x-direction of the 
*     (ix,iy) cell and its left face.
*     t is in the range (-pi/2 to pi/2) and vanishes for an orthogonal
*     geometry. 
*     It will hold that -1.lt.qs(,,).le.1.
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
*-----------------------------------------------------------------------
*.declarations

*   ..local variables
      integer n2, ix, iy, i, j
      integer volhmk
      real (kind=R8) ::
     *   p0x, p0y, p1x, p1y, width, centroid(0:1)
*   ..procedures
      intrinsic min, max, sqrt
      real (kind=R8) ::
     *   trim1, norm, dist, sinang, cosang
      external subini, subend, xertst, ipgetr, ipgeti
      external b2xvsg
      norm(p0x,p0y) = sqrt(p0x**2+p0y**2)
      trim1(p0x) = max(-1.0_R8,min(1.0_R8,p0x))
      dist(ix,iy,i,j) = norm (
     &  crx(ix,iy,j)-crx(ix,iy,i), cry(ix,iy,j)-cry(ix,iy,i))
      sinang(p0x,p0y,p1x,p1y) = trim1(
     &  (p0x*p1x+p0y*p1y)/(norm(p0x,p0y)*norm(p1x,p1y)))
      cosang(p0x,p0y,p1x,p1y) = trim1(
     &  (p0x*p1y-p0y*p1x)/(norm(p0x,p0y)*norm(p1x,p1y)))
*   ..initialisation
      data volhmk/1/, width /1.0_R8/
      save volhmk, width

*-----------------------------------------------------------------------
*.documentation-internal
*
*     Let (p0x,p0y) specify a vector p0, and (p1x,p1y) specify
*     a vector p1. Then norm(p0x,p0y) is the standard norm of p0,
*     sinang(p0x,p0y,p1x,p1y) is sin(t), and cosang(p0x,p0y,p1x,p1y)
*     is cos(t), where t is the complementary angle between p0 and p1.
*     t=pi/2 when p0 and p1 are parallel, t=0 or t=pi when p0 and p1
*     are orthogonal, and t=-pi/2 when p0 and p1 are anti-parallel.
*
*-----------------------------------------------------------------------
*.computation

* ..preliminaries
*   ..subprogram start-up calls
      call subini ('b2agmt')
*   ..read input data
      call ipgeti ('b2agmt_volhmk', volhmk)
      call ipgetr ('b2agmt_1d_width', width)
      call xertst (0.0_R8.lt.width, 'faulty input width')
*   ..set n2
      n2 = (nx+2)*(ny+2)
*   ..test nx, ny, isymm
      call xertst (0.le.nx.and.0.le.ny, 'faulty argument nx, ny')
      call xertst (0.le.isymm.and.isymm.le.4,
     &  'faulty argument isymm')

* ..compute metric coefficients
*   ..compute hx, hy, qz, qc
      do iy = -1, ny
       do ix = -1, nx
        hx(ix,iy) = norm((crx(ix,iy,1)+crx(ix,iy,3))/2.0_R8-
     &                   (crx(ix,iy,0)+crx(ix,iy,2))/2.0_R8,
     &                   (cry(ix,iy,1)+cry(ix,iy,3))/2.0_R8-
     &                   (cry(ix,iy,0)+cry(ix,iy,2))/2.0_R8)
        hy(ix,iy) = norm((crx(ix,iy,2)+crx(ix,iy,3))/2.0_R8-
     &                   (crx(ix,iy,0)+crx(ix,iy,1))/2.0_R8,
     &                   (cry(ix,iy,2)+cry(ix,iy,3))/2.0_R8-
     &                   (cry(ix,iy,0)+cry(ix,iy,1))/2.0_R8) 
        cell_width(ix,iy) = hx(ix,iy)
        cell_height(ix,iy) = hy(ix,iy)
        qz(ix,iy,0) = sinang (
     &    (crx(ix,iy,1)+crx(ix,iy,3))-(crx(ix,iy,0)+crx(ix,iy,2)),
     &    (cry(ix,iy,1)+cry(ix,iy,3))-(cry(ix,iy,0)+cry(ix,iy,2)),
     &    (crx(ix,iy,2)+crx(ix,iy,3))-(crx(ix,iy,0)+crx(ix,iy,1)),
     &    (cry(ix,iy,2)+cry(ix,iy,3))-(cry(ix,iy,0)+cry(ix,iy,1)))
        qz(ix,iy,1) = cosang (
     &    (crx(ix,iy,1)+crx(ix,iy,3))-(crx(ix,iy,0)+crx(ix,iy,2)),
     &    (cry(ix,iy,1)+cry(ix,iy,3))-(cry(ix,iy,0)+cry(ix,iy,2)),
     &    (crx(ix,iy,2)+crx(ix,iy,3))-(crx(ix,iy,0)+crx(ix,iy,1)),
     &    (cry(ix,iy,2)+cry(ix,iy,3))-(cry(ix,iy,0)+cry(ix,iy,1)))
        qc(ix,iy,0) = cosang (
     &    (crx(ix,iy,1)+crx(ix,iy,3))-(crx(ix,iy,0)+crx(ix,iy,2)),
     &    (cry(ix,iy,1)+cry(ix,iy,3))-(cry(ix,iy,0)+cry(ix,iy,2)),
     &    crx(ix,iy,2)-crx(ix,iy,0), cry(ix,iy,2)-cry(ix,iy,0))
        qc(ix,iy,1) = 1.0_R8
        qs(ix,iy,0) = sinang (
     &    (crx(ix,iy,1)+crx(ix,iy,3))-(crx(ix,iy,0)+crx(ix,iy,2)),
     &    (cry(ix,iy,1)+cry(ix,iy,3))-(cry(ix,iy,0)+cry(ix,iy,2)),
     &    crx(ix,iy,2)-crx(ix,iy,0), cry(ix,iy,2)-cry(ix,iy,0))
        qs(ix,iy,1) = 0.0_R8
        if(qz(ix,iy,1).le.0.0_R8) then
          write(*,'(a,1x,1f9.6,2i4)')
     1     'Problem in b2agmt: qz(ix,iy,1) <= 0',qz(ix,iy,1),ix,iy
          write(*,'(a,1p,4g13.6)') 'crx = ',(crx(ix,iy,i),i=0,3)
          write(*,'(a,1p,4g13.6)') 'cry = ',(cry(ix,iy,i),i=0,3)
        endif
        if(qc(ix,iy,0).le.0.0_R8) then
          write(*,'(a,1x,1f9.6,2i4)')
     1     'Problem in b2agmt: qc(ix,iy) <= 0',qc(ix,iy,0),ix,iy
          write(*,'(a,1p,4g13.6)') 'crx = ',(crx(ix,iy,i),i=0,3)
          write(*,'(a,1p,4g13.6)') 'cry = ',(cry(ix,iy,i),i=0,3)
        endif
        if(qc(ix,iy,1).le.0.0_R8) then
          write(*,'(a,1f12.6,2i4)')
     1     'Problem in b2agmt: qc(ix,iy,1) <= 0',qc(ix,iy,1),ix,iy
          write(*,'(a,4f13.6)') 'crx = ',(crx(ix,iy,i),i=0,3)
          write(*,'(a,4f13.6)') 'cry = ',(cry(ix,iy,i),i=0,3)
        endif
       enddo
      enddo
*   ..test sign of qz(,,1) and qc(,)
*     (A negative entry in qz(,,1) corresponds to a left-handed
*     orientation of the coordinate system, which is not allowed in
*     this code.)
      call b2xvsg (n2, qz(-1,-1,1), 1, 'qz1', '.gt.')
      call b2xvsg (n2*2, qc, 1, 'qc', '.gt.')
*   ..compute vol, gs
      if (isymm.eq.0) then
*    ..slab geometry
       do iy = -1, ny
        do ix = -1, nx
         gs(ix,iy,2) = 0.5_R8*
     &     ((crx(ix,iy,3)-crx(ix,iy,0))*(cry(ix,iy,2)-cry(ix,iy,1))-
     &      (crx(ix,iy,2)-crx(ix,iy,1))*(cry(ix,iy,3)-cry(ix,iy,0)))
         vol(ix,iy) = gs(ix,iy,2)*width
        enddo
       enddo
       do iy = -1, ny
        do ix = -1, nx
         if(leftix(ix,iy).eq.-2) then
          gs(ix,iy,0) = 0.0_R8
         else
          gs(ix,iy,0) = dist(ix,iy,0,2)*width
         endif
         if(bottomiy(ix,iy).eq.-2) then
          gs(ix,iy,1) = 0.0_R8
         else
          gs(ix,iy,1) = dist(ix,iy,0,1)*width
         endif
        enddo
       enddo
      else if (isymm.eq.1 .or. isymm.eq.2) then
*    ..case of rotational symmetry about crx=0 axis
       do iy = -1, ny
        do ix = -1, nx
         gs(ix,iy,2) = 0.5_R8*
     &     ((crx(ix,iy,3)-crx(ix,iy,0))*(cry(ix,iy,2)-cry(ix,iy,1))-
     &      (crx(ix,iy,2)-crx(ix,iy,1))*(cry(ix,iy,3)-cry(ix,iy,0)))
         if (volhmk.eq.0) then
          centroid(0) = 0.25_R8*
     &     (crx(ix,iy,0)+crx(ix,iy,1)+crx(ix,iy,2)+crx(ix,iy,3))
          centroid(1) = 0.25_R8*
     &     (cry(ix,iy,0)+cry(ix,iy,1)+cry(ix,iy,2)+cry(ix,iy,3))
         else
          centroid(0) = XquadCentroid(crx(ix,iy,0),cry(ix,iy,0),
     &                                crx(ix,iy,1),cry(ix,iy,1),
     &                                crx(ix,iy,2),cry(ix,iy,2),
     &                                crx(ix,iy,3),cry(ix,iy,3))
          centroid(1) = YquadCentroid(crx(ix,iy,0),cry(ix,iy,0),
     &                                crx(ix,iy,1),cry(ix,iy,1),
     &                                crx(ix,iy,2),cry(ix,iy,2),
     &                                crx(ix,iy,3),cry(ix,iy,3))
         endif
         vol(ix,iy) = gs(ix,iy,2)*2.0_R8*pi*centroid(0)
        enddo
       enddo
       do iy = -1, ny
        do ix = -1, nx
         if(leftix(ix,iy).eq.-2) then
          gs(ix,iy,0) = 0.0_R8
         else
          gs(ix,iy,0) = pi*dist(ix,iy,0,2)*(crx(ix,iy,0)+crx(ix,iy,2))
         endif
         if(bottomiy(ix,iy).eq.-2) then
          gs(ix,iy,1) = 0.0_R8
         else
          gs(ix,iy,1) = pi*dist(ix,iy,0,1)*(crx(ix,iy,0)+crx(ix,iy,1))
         endif
        enddo
       enddo
      else if (isymm.eq.3 .or. isymm.eq.4) then
*    ..case of rotational symmetry about cry=0 axis
       do iy = -1, ny
        do ix = -1, nx
         gs(ix,iy,2) = 0.5_R8*
     &     ((crx(ix,iy,3)-crx(ix,iy,0))*(cry(ix,iy,2)-cry(ix,iy,1))-
     &      (crx(ix,iy,2)-crx(ix,iy,1))*(cry(ix,iy,3)-cry(ix,iy,0)))
         if (volhmk.eq.0) then
          centroid(0) = 0.25_R8*
     &     (crx(ix,iy,0)+crx(ix,iy,1)+crx(ix,iy,2)+crx(ix,iy,3))
          centroid(1) = 0.25_R8*
     &     (cry(ix,iy,0)+cry(ix,iy,1)+cry(ix,iy,2)+cry(ix,iy,3))
         else
          centroid(0) = XquadCentroid(crx(ix,iy,0),cry(ix,iy,0),
     &                                crx(ix,iy,1),cry(ix,iy,1),
     &                                crx(ix,iy,2),cry(ix,iy,2),
     &                                crx(ix,iy,3),cry(ix,iy,3))
          centroid(1) = YquadCentroid(crx(ix,iy,0),cry(ix,iy,0),
     &                                crx(ix,iy,1),cry(ix,iy,1),
     &                                crx(ix,iy,2),cry(ix,iy,2),
     &                                crx(ix,iy,3),cry(ix,iy,3))
         endif
         vol(ix,iy) = gs(ix,iy,2)*2.0_R8*pi*centroid(1)
        enddo
       enddo
       do iy = -1, ny
        do ix = -1, nx
         if(leftix(ix,iy).eq.-2) then
          gs(ix,iy,0) = 0.0_R8
         else
          gs(ix,iy,0) = pi*dist(ix,iy,0,2)*(cry(ix,iy,0)+cry(ix,iy,2))
         endif
         if(bottomiy(ix,iy).eq.-2) then
          gs(ix,iy,1) = 0.0_R8
         else
          gs(ix,iy,1) = pi*dist(ix,iy,0,1)*(cry(ix,iy,0)+cry(ix,iy,1))
         endif
        enddo
       enddo
      endif
*   ..test sign of vol
      call b2xvsg (n2, vol, 1, 'vol', '.gt.')

* ..return
      call subend ()
      return

*-----------------------------------------------------------------------
*.end b2agmt

      contains

C Coordinates of centroid of a quadrangle given by vertex coordinates
C Corners are numbered in the B2 convention 2-3
C                                           0-1
      function XquadCentroid( x0, y0, x1, y1, x2, y2, x3, y3 )
      use b2mod_types
      implicit none
      real (kind=R8) :: XquadCentroid
      real (kind=R8) :: triangleArea
      external triangleArea
      real (kind=R8) :: x0, y0, x1, y1, x2, y2, x3, y3
  
      XquadCentroid =
     & ( XtriangleCentroid(x0,y0,x1,y1,x3,y3)
     & * triangleArea(x0,y0,x1,y1,x3,y3)
     & + XtriangleCentroid(x0,y0,x2,y2,x3,y3)
     & * triangleArea(x0,y0,x2,y2,x3,y3) )
     & / ( triangleArea(x0,y0,x1,y1,x3,y3)
     &   + triangleArea(x0,y0,x2,y2,x3,y3) )
      return
      end function XquadCentroid

      function YquadCentroid( x0, y0, x1, y1, x2, y2, x3, y3 )
      use b2mod_types
      implicit none
      real (kind=R8) :: YquadCentroid
      real (kind=R8) :: triangleArea
      external triangleArea
      real (kind=R8) :: x0, y0, x1, y1, x2, y2, x3, y3
  
      YquadCentroid =
     & ( YtriangleCentroid(x0,y0,x1,y1,x3,y3)
     & * triangleArea(x0,y0,x1,y1,x3,y3)
     & + YtriangleCentroid(x0,y0,x2,y2,x3,y3)
     & * triangleArea(x0,y0,x2,y2,x3,y3) )
     & / ( triangleArea(x0,y0,x1,y1,x3,y3)
     &   + triangleArea(x0,y0,x2,y2,x3,y3) )
      return
      end function YquadCentroid

C Coordinates of centroid of a triangle given by vertex coordinates
      function XtriangleCentroid( x0, y0, x1, y1, x2, y2 )
      use b2mod_types
      implicit none
      real (kind=R8) :: XtriangleCentroid
      real (kind=R8) :: x0, y0, x1, y1, x2, y2

      XtriangleCentroid = ( x0 + x1 + x2 ) / 3.0_R8
      return
      end function XtriangleCentroid

      function YtriangleCentroid( x0, y0, x1, y1, x2, y2 )
      use b2mod_types
      implicit none
      real (kind=R8) :: YtriangleCentroid
      real (kind=R8) :: x0, y0, x1, y1, x2, y2

      YtriangleCentroid = ( y0 + y1 + y2 ) / 3.0_R8
      return
      end function YtriangleCentroid

      end subroutine b2agmt

C Area of triangle given by vertex coordinates
      function triangleArea( x0, y0, x1, y1, x2, y2 )
      use b2mod_types
      implicit none
      real (kind=R8) :: triangleArea
      real (kind=R8) :: x0, y0, x1, y1, x2, y2
      intrinsic abs
  
      triangleArea = 0.5_R8 * abs( (x1-x0)*(y2-y0)
     & - (y1-y0)*(x2-x0) )
      return
      end function triangleArea
