*-----------------------------------------------------------------------
*.specification

*      subroutine b2agmt (nx, ny, isymm, crx, cry,
*     & fpsi, bb, vol, hx, hy, hc, ht, qz, qc, qs, qac, qas, ebc,
*     & gs, cell_width, cell_height, leftix, leftiy, bottomix, bottomiy, 
*     & rightix, rightiy, topix, topiy, cflag)
      subroutine b2agmt (g, m)
      use b2mod_types
      use b2mod_constants
      use tradui_constants
      use b2mod_connectivity
      use b2mod_cellhelper
!      use b2mod_indirect
!     , ,only : TriangleArea, TriangleCentroid, quadCentroid, quadArea
!     use b2mod_geo
      use b2us_geo
      use b2us_map
 
      implicit none
*   ..input arguments (unchanged on exit)
      integer isymm
!      integer nx, ny, isymm
!      integer leftix(-1:nx,-1:ny), leftiy(-1:nx,-1:ny), 
!     &  bottomix(-1:nx,-1:ny), bottomiy(-1:nx,-1:ny),
!     &  rightix(-1:nx,-1:ny), rightiy(-1:nx,-1:ny),
!     &  topix(-1:nx,-1:ny), topiy(-1:nx,-1:ny),
!     $  cflag(-1:nx,-1:ny,GRID_N_CELLFLAGS) 
!      real (kind=R8) ::
!     *  crx(-1:nx,-1:ny,0:3), cry(-1:nx,-1:ny,0:3),
!     &  fpsi(-1:nx,-1:ny,0:3,0:2), bb(-1:nx,-1:ny,0:3)
*   ..output arguments (unspecified on entry)
!      real (kind=R8) ::
!     *  vol(-1:nx,-1:ny), hx(-1:nx,-1:ny), hy(-1:nx,-1:ny),
!     &  qz(-1:nx,-1:ny,0:1), qc(-1:nx,-1:ny,0:1), qs(-1:nx,-1:ny,0:1),
!     &  gs(-1:nx,-1:ny,0:2), qac(-1:nx,-1:ny,0:1), qas(-1:nx,-1:ny,0:1),
!     &  ebc(-1:nx,-1:ny,0:2), 
!     &  cell_width(-1:nx,-1:ny), cell_height(-1:nx,-1:ny),
!     &  hc(-1:nx,-1:ny,0:3), ht(-1:nx,-1:ny,0:1)

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
      integer ix, iy, i, j, iCv, iFc
      integer verts4(0:3), verts3(0:2), faces3(0:2)
      type(geometry) :: g
      type(mapping) :: m
      integer i1, i2, ii
      integer volhmk, cgeo, cgeon
      real (kind=R8) ::
     *   p0x, p0y, p1x, p1y, width, centroid(0:1), t0, t1, psic
      logical classical
!      real (kind=R8) ::
!     *   ecc(-1:nx,-1:ny,0:3,0:1), ecf(0:1), ebf(0:1),
!     *   isx(0:1), isy(0:1)
      real (kind=R8) :: 
     *   fcX(1:m%nFc), fcY(1:m%nFc) 
      integer, save :: style = 0
*   ..procedures
      intrinsic min, max, sqrt
      real (kind=R8) :: norm
!     *   trim1, norm
!     * sinang, cosang
!     *   dist, 
      real (kind=R8) ::
     &  sin_t, cos_t, sin_tl, cos_tl, tt
      external subini, subend, xertst, ipgetr, ipgeti
      external b2xvsg
      norm(p0x,p0y) = sqrt(p0x**2+p0y**2)
!      trim1(p0x) = max(-1.0_R8,min(1.0_R8,p0x))
!      dist(ix,iy,i,j) = norm (
!     &  crx(ix,iy,j)-crx(ix,iy,i), cry(ix,iy,j)-cry(ix,iy,i))
!      sinang(p0x,p0y,p1x,p1y) = trim1(
!     &  (p0x*p1x+p0y*p1y)/(norm(p0x,p0y)*norm(p1x,p1y)))
!      cosang(p0x,p0y,p1x,p1y) = trim1(
!     &  (p0x*p1y-p0y*p1x)/(norm(p0x,p0y)*norm(p1x,p1y)))
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
      call ipgeti ('b2agmt_style', style)
      call xertst (0.0_R8.lt.width, 'faulty input width')
*   ..test nx, ny, isymm
      call xertst (0.le.m%nCv.and.0.le.m%nFc.and.0.le.m%nVx,
     &  'faulty argument nCv, nFc, nVx')
      call xertst (0.le.isymm.and.isymm.le.4,
     &  'faulty argument isymm')

* ..compute metric coefficients
!      classical = isClassicalGrid(cflag)

cnv  Exploit the structuredness of a classical grid in unstructured format
cnv  similar methods as in b2agmt_st, but first the ording of the vertices 
cnv  need to be checked and adapted if necessary. For non-classical grid,
cnv  new methods are employed
!      if (m%isClassicalGrid.eq.1) then
!         call order_vertices(g,m)
!      endif 


*   ..compute hx, hy, hc, ht, ecc, ebc
      !hx = INVALID_DOUBLE
      !hy = INVALID_DOUBLE
      !cell_width = INVALID_DOUBLE
      !cell_height = INVALID_DOUBLE
      !qz(:,:,0) = 0.0_R8
      !qz(:,:,1) = 1.0_R8



! loop over cells
      do iCv = 1,m%nCv

       
       select case (m%cvFcP(iCv,2))

       case(4)  !Quad

       !ebc unit vector of poloidal field in cell center
       verts4 = m%cvVx(m%cvVxP(iCv,1):m%cvVxP(iCv,1) + 
     &          m%cvVxP(iCv,2)-1 )
       ! b_y
       t0 = 0.25_R8 * ( g%vxBb(verts4(0),2) + g%vxBb(verts4(1),2) +   
     &        g%vxBb(verts4(2),2) + g%vxBb(verts4(3),2) )
       ! b_x
       t1 = - 0.25_R8 * ( g%vxBb(verts4(0),1) + g%vxBb(verts4(1),1) +
     &        g%vxBb(verts4(2),1) + g%vxBb(verts4(3),1) )

       g%cvEb(iCv,0) = t0 / norm(t0,t1)
       g%cvEb(iCv,1) = t1 / norm(t0,t1)

*    ..unit vector of total magnetic field, defined in cell center
      ! until this point, cvEb was unit vector of poloidal field
      ! now, correct for toroidal field
       if (g%cvBb(iCv,3) == 0.0_R8) Then
         g%cvEb(iCv,:) = 0.0_R8
       else
         g%cvEb(iCv,0) = abs(g%cvBb(iCv,0))/g%cvBb(iCv,3)*g%cvEb(iCv,0)
         g%cvEb(iCv,1) = abs(g%cvBb(iCv,1))/g%cvBb(iCv,3)*g%cvEb(iCv,1)
         g%cvEb(iCv,2) = g%cvBb(iCv,2)/g%cvBb(iCv,3)
       endif 

*    ..compute volume, gs
      ! g%cvSz = gs(:,:,2)
      ! g%fcS = gs(:,:,0:1)

      ! quadArea compute area of quads with random order of vertices
      !g%cvSz(iCv) = quadArea(g%vxX(verts4(0)), g%vxY(verts4(0)),  
      !&     g%vxX(verts4(1)), g%vxY(verts4(1)), g%vxX(verts4(2)),
      !&     g%vxY(verts4(2)), g%vxX(verts4(3)), g%vxY(verts4(3)) )

      ! we assume the vertices are ordered (done in Make_traduit_b2us)
      ! This should be introduced in b2ag!!
      ! Because of this we can employ the quadArea_st which assumes the
      ! vertices in the correct sequence
       g%cvSz(iCv) = quadArea_st(g%vxX(verts4(0)), g%vxY(verts4(0)),  
     &     g%vxX(verts4(1)), g%vxY(verts4(1)), g%vxX(verts4(2)),
     &     g%vxY(verts4(2)), g%vxX(verts4(3)), g%vxY(verts4(3)) )


      if  (isymm.eq.0) then  !slab geometry
      ! volume is area times width
        g%cvVol = g%cvSz *width     
      elseif (isymm.eq.1 .or. isymm.eq.2) then ! rotational symmertry about crx=0 axis
        if (volhmk.eq.0) then 
          centroid(0) = 0.25*R8*
     &      (g%vxX(verts4(0)) + g%vxX(verts4(1)) + g%vxX(verts4(2)) 
     &       + g%vxX(verts4(3)) )
          centroid(1) = 0.25*R8*
     &      (g%vxY(verts4(0)) + g%vxY(verts4(1)) + g%vxY(verts4(2)) 
     &       + g%vxY(verts4(3)) )
        else
          ! QUADCENTROID NOT YET ADAPTED for unstructured - WIP
          centroid = quadCentroid_us(g%vxX(verts4(0)),g%vxY(verts4(0)),
     &                            g%vxX(verts4(1)),g%vxY(verts4(1)),
     &                            g%vxX(verts4(2)),g%vxY(verts4(2)),
     &                            g%vxX(verts4(3)),g%vxY(verts4(3)))
        endif
        ! volume is area time circumference of circle with radius equal to the x-coordinate of centroid 
        g%cvVol(iCv) = g%cvSz(iCv)*2.0_R8*pi*centroid(0)
      elseif  (isymm.eq.3 .or. isymm.eq.4) then ! rotational symmetry about cry=0
        if (volhmk.eq.0) then
          centroid(0) = 0.25*R8*
     &      (g%vxX(verts4(0)) + g%vxX(verts4(1)) + g%vxX(verts4(2)) 
     &       + g%vxX(verts4(3)) )
          centroid(1) = 0.25*R8*
     &      (g%vxY(verts4(0)) + g%vxY(verts4(1)) + g%vxY(verts4(2)) 
     &       + g%vxY(verts4(3)) )
        else 
          ! old function quadCentroid can be used if vertices are ordered correctly
          centroid = quadCentroid_us(g%vxX(verts4(0)),g%vxY(verts4(0)),
     &                            g%vxX(verts4(1)),g%vxY(verts4(1)),
     &                            g%vxX(verts4(2)),g%vxY(verts4(2)),
     &                            g%vxX(verts4(3)),g%vxY(verts4(3)))
        endif
        ! volume is area time circumference of circle with radius equal to the y-coordinate of centroid 
        g%cvVol(iCv) = g%cvSz(iCv)*2.0_R8*pi*centroid(1)       

      endif
      ! compute ecc : Unit vectors of connector lines cell center to face centers
      if (m%isClassicalGrid) then
         ! vertices need to be ordered as follows
         ! verts4(0) = left under
         ! verts4(1) = right under
         ! verts4(2) = right top
         ! verts4(3) = left top

      else


      endif
        



 

       case(3)  !Triangle
       ! Method Wouter to determine cell center for triangles!!! line 290 in structured version
cwdk      General approach for arbitrary triangle cases:
cwdk      1) identify vertices with max and min psi values psi_min, psi_max
cwdk      2) assume the triangle lies on the flux surface with psic = 0.5*(psi_min + psi_max)
cwdk         => through middle of the corresponding face; defines 1st intersection
cwdk      3) identify 2nd intersection of this surface with one of the other two faces
cwdk      4) centroid: average between the intersections with faces       

       ! Identify vertices with min and max psi values
       verts3 = m%cvVx(m%cvVxP(iCv,1):m%cvVxP(iCv,1) + 
     &          m%cvVxP(iCv,2)-1 )
       i1 = minloc(g%vxFpsi(verts3(0:2)),1)-1
       i2 = maxloc(g%vxFpsi(verts3(0:2)),1)-1
       psic = 0.5_R8*(g%vxFpsi(verts3(i1)) + g%vxFpsi(verts3(i2)))

      ! Loop over faces to find which ones are intersected by the psic surface
      ! (Loop includes degenerate face; should not cause problems)
      faces3 = m%cvFc(m%cvFcP(iCv,1):m%cvFcP(iCv,1) + 
     &          m%cvFcP(iCv,2)-1 )
      ii = 0 
      


       end select
       
         

       

       


      enddo

* ..compute qc, qs, qac, qas at cell faces
      ! Requirements for US solver:
      ! qc:  cosine of angle between connector line and face normal
      ! qs:    sine of angle between connector line and face normal
      !      => qc and qs will translate into fcQgam
      ! qac: cosine of angle between poloidal field and face normal
      ! qas:   sine of angle between poloidal field and face normal
      !     => qac and qas will translate into fcQalf
      !qc  = 1.0_R8
      !qs  = 0.0_R8
      !qac = 0.0_R8
      !qas = 1.0_R8

! loop over faces (in st first poloidal faces, then radial faces)
      do iFc = 1, m%nFc
        !compute face centers
        fcX(iFc) = (g%vxX(m%fcVx(iFc,1)) + g%vxX(m%fcVx(iFc,2)) )/2.0_R8
        fcY(iFc) = (g%vxY(m%fcVx(iFc,1)) + g%vxY(m%fcVx(iFc,2)) )/2.0_R8

        !fcHc(m%nFc,2) : Lengths between cell face center and neighboring cell centers ("connector" direction)
        g%fcHc(iFc,1) = sqrt(  ( g%cvX(m%fcCv(iFc,1)) - fcX(iFc) )**2
     &      + ( g%cvY(m%fcCV(iFc,1)) - fcY(iFc) )**2 ) 
        g%fcHc(iFc,2) = sqrt(  ( g%cvX(m%fcCv(iFc,2)) - fcX(iFc) )**2
     &      + ( g%cvY(m%fcCV(iFc,2)) - fcY(iFc) )**2 ) 

        !fcHt(m%nFc) Lengths along cell face ("tangential" direction)
        g%fcHt(iFc) = sqrt( ( g%vxX(m%fcVx(iFc,1)) - 
     &     g%vxX(m%fcVx(iFc,2)) )**2 + ( g%vxY(m%fcVx(iFc,1)) - 
     &     g%vxY(m%fcVx(iFc,2)) )**2  )


       ! g%fcS = gs(:,:,0:1)
       if (isymm.eq.0) then
          g%fcS(iFc) = g%fcHt(iFc) * width
       elseif (isymm.eq.1 .or. isymm.eq.2) then ! rotational symmetry about crx=0 axis
          g%fcS(iFc) = pi*g%fcHt(iFc) * (g%vxX(m%fcVx(iFc,1)) +
     &                                   g%vxX(m%fcVx(iFc,2)) )
       elseif (isymm.eq.3 .or. isymm.eq.4) then
          g%fcS(iFc) = pi*g%fcHt(iFc) * (g%vxY(m%fcVx(iFc,1)) +
     &                                   g%vxY(m%fcVx(iFc,2)) )

       endif 
      enddo
    


*    ..test sign of g%cvVol, g%cvS, g%fcS
       ! based on b2xvs1





* ..return
      call subend ()
      return

*-----------------------------------------------------------------------
*.end b2agmt

      end subroutine b2agmt


