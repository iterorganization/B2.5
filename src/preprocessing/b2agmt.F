*-----------------------------------------------------------------------
*.specification

*      subroutine b2agmt (nx, ny, isymm, crx, cry,
*     & fpsi, bb, vol, hx, hy, hc, ht, qz, qc, qs, qac, qas, ebc,
*     & gs, cell_width, cell_height, leftix, leftiy, bottomix, bottomiy, 
*     & rightix, rightiy, topix, topiy, cflag)
      subroutine b2agmt (g, m, isymm)
      use b2mod_types
      use b2mod_constants
      use tradui_constants
      use carre_constants
      !use b2mod_connectivity
      !use b2mod_cellhelper
      use b2us_geo
      use b2us_map
 
      implicit none
*   ..input arguments (unchanged on exit)
      integer isymm
      type(geometry) :: g
      type(mapping) :: m
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
      integer ix, iy, i, j, iCv, iFc, iF, face, gc, bc, v1, v2
      integer verts4(0:3), verts3(0:2), faces3(0:2), cells(1:2), 
     * fcType(m%nFc)
      integer i1, i2, ii
      integer volhmk, cgeo, cgeon
      real (kind=R8) ::
     *   p0x, p0y, p1x, p1y, width, centroid(m%nCv,0:1), bcFpsi,
     *   t0, t1, psic, dx, dy, cf(0:1), isx(0:1), isy(0:1),
     *   guard_rescale, fvec(0:1), dg, fcPsi, ec1fc(0:1), efcc2(0:1),
     *   econf(0:1), fn(0:1), tol_angles, tol_aligned
      logical classical
!      real (kind=R8) ::
!     *   ecc(-1:nx,-1:ny,0:3,0:1), ecf(0:1), ebf(0:1),
      real (kind=R8) :: 
     *   fcX(1:m%nFc), fcY(1:m%nFc) 
      integer, save :: style = 0
*   ..procedures
      intrinsic min, max, sqrt
      real (kind=R8) :: norm, trim1,sinang,cosang
      real (kind=R8) ::
     &  sin_t, cos_t, sin_tl, cos_tl, tt
      external subini, subend, xertst, ipgetr, ipgeti
      external b2xvsg
      norm(p0x,p0y) = sqrt(p0x**2+p0y**2)
      trim1(p0x) = max(-1.0_R8,min(1.0_R8,p0x))
!      dist(ix,iy,i,j) = norm (
!     &  crx(ix,iy,j)-crx(ix,iy,i), cry(ix,iy,j)-cry(ix,iy,i))
      sinang(p0x,p0y,p1x,p1y) = trim1(
     &  (p0x*p1x+p0y*p1y)/(norm(p0x,p0y)*norm(p1x,p1y)))
      cosang(p0x,p0y,p1x,p1y) = trim1(
     &  (p0x*p1y-p0y*p1x)/(norm(p0x,p0y)*norm(p1x,p1y)))
*   ..initialisation
      data volhmk/1/, width /1.0_R8/, guard_rescale /1.0e-3_R8/,
     &  tol_aligned /1.0e-4_R8/
      save volhmk, width, guard_rescale, tol_aligned

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
      call ipgetr ('b2agmt_guard_rescale',guard_rescale)
      call xertst (0.0_R8.lt.width, 'faulty input width')
      call xertst (0.0_R8.lt.guard_rescale, 
     &   'faulty input guard rescale')
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
       t0 = 0.25_R8 * ( g%vxBb(verts4(0),1) + g%vxBb(verts4(1),1) +   
     &        g%vxBb(verts4(2),1) + g%vxBb(verts4(3),1) )
       ! b_x
       t1 = - 0.25_R8 * ( g%vxBb(verts4(0),0) + g%vxBb(verts4(1),0) +
     &        g%vxBb(verts4(2),0) + g%vxBb(verts4(3),0) )

       g%cvEb(iCv,0) = t0 / norm(t0,t1)
       g%cvEb(iCv,1) = t1 / norm(t0,t1)

*    ..unit vector of total magnetic field, defined in cell center
      ! until this point, cvEb was unit vector of poloidal field
      ! now, correct for toroidal field
       if (g%cvBb(iCv,3) == 0.0_R8) then
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

      ! we assume the vertices are ordered 
      ! Because of this we can employ the quadArea_st which assumes the
      ! vertices in the correct sequence
       g%cvSz(iCv) = quadArea_st(g%vxX(verts4(0)), g%vxY(verts4(0)),  
     &     g%vxX(verts4(1)), g%vxY(verts4(1)), g%vxX(verts4(3)),
     &     g%vxY(verts4(3)), g%vxX(verts4(2)), g%vxY(verts4(2)) )


      if  (isymm.eq.0) then  !slab geometry
      ! volume is area times width
        g%cvVol(iCv) = g%cvSz(iCv) *width     
      elseif (isymm.eq.1 .or. isymm.eq.2) then ! rotational symmertry about crx=0 axis
        if (volhmk.eq.0) then 
          centroid(iCv,0) = 0.25_R8*
     &      (g%vxX(verts4(0)) + g%vxX(verts4(1)) + g%vxX(verts4(2)) 
     &       + g%vxX(verts4(3)) )
          centroid(iCv,1) = 0.25_R8*
     &      (g%vxY(verts4(0)) + g%vxY(verts4(1)) + g%vxY(verts4(2)) 
     &       + g%vxY(verts4(3)) )
        else
          centroid(iCv,:) = quadCentroid_us(g%vxX(verts4(0)),
     &          g%vxY(verts4(0)),g%vxX(verts4(1)),g%vxY(verts4(1)),           
     &                           g%vxX(verts4(2)),g%vxY(verts4(2)),
     &                           g%vxX(verts4(3)),g%vxY(verts4(3)))

          ! Calculate the centroid similar to the method for the triangle centroid
          call calc_psi_centroid(g,m,centroid,iCv)   

        endif
        ! volume is area time circumference of circle with radius equal to the x-coordinate of centroid 
        g%cvVol(iCv) = g%cvSz(iCv)*2.0_R8*pi*centroid(iCv,0)
      elseif  (isymm.eq.3 .or. isymm.eq.4) then ! rotational symmetry about cry=0
        if (volhmk.eq.0) then
          centroid(iCv,0) = 0.25*R8*
     &      (g%vxX(verts4(0)) + g%vxX(verts4(1)) + g%vxX(verts4(2)) 
     &       + g%vxX(verts4(3)) )
          centroid(iCv,1) = 0.25*R8*
     &      (g%vxY(verts4(0)) + g%vxY(verts4(1)) + g%vxY(verts4(2)) 
     &       + g%vxY(verts4(3)) )
        else ! default
          ! old function quadCentroid can be used if vertices are ordered correctly
          centroid(iCv,0:1) = quadCentroid_us(g%vxX(verts4(0)),
     &           g%vxY(verts4(0)),g%vxX(verts4(1)),g%vxY(verts4(1)),       
     &                            g%vxX(verts4(3)),g%vxY(verts4(3)),
     &                            g%vxX(verts4(2)),g%vxY(verts4(2)))


          ! Calculate the centroid similar to the method for the triangle centroid
          call calc_psi_centroid(g,m,centroid,iCv)
          
        endif
        ! volume is area time circumference of circle with radius equal to the y-coordinate of centroid 
        g%cvVol(iCv) = g%cvSz(iCv)*2.0_R8*pi*centroid(iCv,1)       

      endif

      !Compute qz which is cvQgam
      ! g%cvQgam(iCv,0) = qz(iCv,1) = cosang()
      ! g%cvQgam(iCv,1) = -qz(iCv,0) = - sinang()



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
      faces3 = m%cvFc(m%cvFcP(iCv,1):m%cvFcP(iCv,1) + 
     &          m%cvFcP(iCv,2)-1 )
      ii = 0 
      do iF = 0,2 
      if (psic.ge.min(g%vxFpsi(m%fcVx(faces3(iF),1)) , 
     &             g%vxFpsi(m%fcVx(faces3(iF),2))) .and.
     &   (psic.lt.max(g%vxFpsi(m%fcVx(faces3(iF),1)),  
     &             g%vxFpsi(m%fcVx(faces3(iF),2)))) ) then !intersected  
        t0 = (psic - g%vxFpsi(m%fcVx(faces3(iF),1))) / 
     &        (g%vxFpsi(m%fcVx(faces3(iF),2)) - 
     &         g%vxFpsi(m%fcVx(faces3(iF),1)))
        isx(ii) = g%vxX(m%fcVx(faces3(iF),1)) + t0 * 
     &  ( g%vxX(m%fcVx(faces3(iF),2)) - g%vxX(m%fcVx(faces3(iF),1)))
        isy(ii) = g%vxY(m%fcVx(faces3(iF),1)) + t0 * 
     &  ( g%vxY(m%fcVx(faces3(iF),2)) - g%vxY(m%fcVx(faces3(iF),1)))     
        ii = ii+1
      endif
      enddo  
      if (ii.ne.2) then
         write(*,*) 'Problem with intersections in triangle ', iCv
         write(*,*) 'Number of intersections found: ' , ii
         call xerrab('Intersections in triangle not found!')
      endif         
      
      !Compute centroid
      centroid(iCv,0) = 0.5_R8*sum(isx)
      centroid(iCv,1) = 0.5_R8*sum(isy)

      !Compute cvSz and cvVol
      g%cvSz = TriangleArea_us(g%vxX(verts3(0)), g%vxY(verts3(0)),
     &                        g%vxX(verts3(1)), g%vxY(verts3(1)),
     &                        g%vxX(verts3(2)), g%vxY(verts3(2)))

      if  (isymm.eq.0) then  !slab geometry
      ! volume is area times width
        g%cvVol(iCv) = g%cvSz(iCv) *width      
      elseif (isymm.eq.1 .or. isymm.eq.2) then ! rotational symmertry about crx=0 axis
        ! volume is area time circumference of circle with radius equal to the x-coordinate of centroid 
        g%cvVol(iCv) = g%cvSz(iCv)*2.0_R8*pi*centroid(iCv,0)
      elseif  (isymm.eq.3 .or. isymm.eq.4) then ! rotational symmetry about cry=0      
        ! volume is area time circumference of circle with radius equal to the y-coordinate of centroid 
        g%cvVol(iCv) = g%cvSz(iCv)*2.0_R8*pi*centroid(iCv,1) 
      endif


      ! then in b2agmt_st => hc, ht, ecc = > now in face-loop

      !Compute ebc



       case(1)  !Guard Cells (are at the end of the cells array so can us information 
                ! on triangles and quad calculated above)
         face = m%cvFc(m%cvFcP(iCv,1)) ! face of the guard cell
         v1 = m%fcVx(face,1)
         v2 = m%fcVx(face,2)
         cells = m%fcCv(face,:)

         fcPsi = 0.5_R8 * (g%vxFpsi(v1) + g%vxFpsi(v2))

         if ((m%cvLbl(cells(1)).eq.GRID_GUARD) .and. 
     &      (m%cvLbl(cells(2)).eq.GRID_BOUNDARY)) then
            gc = cells(1)   !guard cell 
            bc = cells(2)   ! boundary cell
         elseif ((m%cvLbl(cells(1)).eq.GRID_BOUNDARY) .and. 
     &      (m%cvLbl(cells(2)).eq.GRID_GUARD)) then
            gc = cells(2)
            bc = cells(1)
         else
            call xerrab
     & ('b2agmt cv-loop, guard: no guard and boundary cell assigned')
         endif 

         if (m%cvFcP(bc,2).eq.4) then !quad
         !psi value at cell center - not correct
         bcFpsi = 0.25_R8*sum(g%vxFpsi(m%cvVx(m%cvVxP(bc,1):
     &        m%cvVxP(bc,1)+m%cvVxP(bc,2)-1)))

         elseif (m%cvFcP(bc,2).eq.3) then ! triangle
         !psi value at cell center according to b2agmt
         ! Identify vertices with min and max psi values
         verts3 = m%cvVx(m%cvVxP(bc,1):m%cvVxP(bc,1) + 
     &            m%cvVxP(bc,2)-1 )
         i1 = minloc(g%vxFpsi(verts3(0:2)),1)-1
         i2 = maxloc(g%vxFpsi(verts3(0:2)),1)-1
         bcFpsi = 0.5_R8*(g%vxFpsi(verts3(i1)) + g%vxFpsi(verts3(i2)))
         
         endif

         !Determine intersection of bcFpsi with face
         if (bcFpsi.gt.min(g%vxFpsi(v1),g%vxFpsi(v2) ) .and.
     &       bcFpsi.lt.max(g%vxFpsi(v1),g%vxFpsi(v2)) )  then
            !Calculate intersection point
            t0 = (bcFpsi - g%vxFpsi(v1))/(g%vxFpsi(v2)-g%vxFpsi(v1))      
            isx(1) =  g%vxX(v1) + t0 *(g%vxX(v2)- g%vxX(v1))
            isy(1) =  g%vxY(v1) + t0 *(g%vxY(v2)- g%vxY(v1))

            !Introduce some cut-off!!
            !Check angle between vector tangential to face and vector between intersection and centroid
            ! Calc vector from cell center to intersection
            cf(0) = isx(1) - centroid(bc,0)
            cf(1) = isy(1) - centroid(bc,1)
            ! Calc vector tangential to boundary face
            fvec(0) = g%vxX(v2) - g%vxX(v1)
            fvec(1) = g%vxY(v2) - g%vxY(v1)
            
            !Calculate the cosine of the inclosed angle
            ! Cut-off is 0.8
            if (abs(sinang(fvec(0),fvec(1),cf(0),cf(1)) ) 
     &               .gt. 0.8_R8 ) then
               ! just use the face center
               isx(1) = 0.5_R8*(g%vxX(v1) + g%vxX(v2))
               isy(1) = 0.5_R8*(g%vxY(v1) + g%vxY(v2))
            endif


         else 
            !no intersection point, so boundary face is close to aligned with 
            ! magnetic field line
            ! just use face center
            isx(1) = 0.5_R8*(g%vxX(v1) + g%vxX(v2))
            isy(1) = 0.5_R8*(g%vxY(v1) + g%vxY(v2))
         endif


          !Determine cell center of guard cell
          !Vector from cell center to face center
          cf(0) = isx(1) - centroid(bc,0)
          cf(1) = isy(1) - centroid(bc,1)

          !Determine dx, dy of guard cell where sqrt(dx**2 + dy**2) is the
          ! the cell width of the guard cell 
          ! assumption: guard_cell is 1e-3 times smaller than boundary cells in outward direction
          dx = 2* cf(0) / guard_rescale
          dy = 2* cf(1) / guard_rescale
            
          dg = sqrt(dx**2+dy**2)
            
          !cell area is dg time face Area
          g%cvSz = dg * sqrt( (g%vxX(v1) - g%vxX(v2))**2 + 
     &                        (g%vxY(v1) - g%vxY(v2))**2)         

          !Cell center of guard is face center or intersection point plus cf/guard_rescale
          ! Vector from face center is direction of (dx/2, dy/2)
          ! (Because cell center of guard cell should be in the middle of the cell)
          centroid(gc,0) = isx(1) + cf(0) / guard_rescale
          centroid(gc,1) = isy(1) + cf(1) / guard_rescale 

          g%cvX(iCv) = centroid(iCv,0)  ! gc is equal to icv
          g%cvY(iCv) = centroid(iCv,1)


         !Determine cvVol for the guard cell
         if  (isymm.eq.0) then  !slab geometry
         ! volume is area times width
         g%cvVol(iCv) = g%cvSz(iCv) *width      
         elseif (isymm.eq.1 .or. isymm.eq.2) then ! rotational symmertry about crx=0 axis
         ! volume is area time circumference of circle with radius equal to the x-coordinate of centroid 
         g%cvVol(iCv) = g%cvSz(iCv)*2.0_R8*pi*centroid(iCv,0)
         elseif  (isymm.eq.3 .or. isymm.eq.4) then ! rotational symmetry about cry=0      
         ! volume is area time circumference of circle with radius equal to the y-coordinate of centroid 
         g%cvVol(iCv) = g%cvSz(iCv)*2.0_R8*pi*centroid(iCv,1) 
         endif



       end select
       
    


      enddo  !cell loop

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


* ..For the angle calculation, the connectivity data fcVx and fcCv need to comply
*  to a right turning scheme. For cell fcCv(iFc,1), the vertices in fcVx need to
*  ordered right turning. Only here (in the code) this is possible to perform 
*  as the guard cells, only now have a proper cell center. When the cell center 
*  coincides with the face center, this would complicate the ordering significantly.

      call order_face_inf(g,m)

      call build_fcType(g,m,fcType,tol_aligned)

! loop over faces (in st first poloidal faces, then radial faces)
      do iFc = 1, m%nFc
        !compute face centers
        fcX(iFc) = (g%vxX(m%fcVx(iFc,1)) + g%vxX(m%fcVx(iFc,2)) )/2.0_R8
        fcY(iFc) = (g%vxY(m%fcVx(iFc,1)) + g%vxY(m%fcVx(iFc,2)) )/2.0_R8

        !fcHc(m%nFc,2) : Lengths between cell face center and neighboring cell centers ("connector" direction)
        g%fcHc(iFc,1) = sqrt(  ( centroid(m%fcCv(iFc,1),0) - fcX(iFc) )**2
     &      + ( centroid(m%fcCv(iFc,1),1) - fcY(iFc) )**2 ) 
        g%fcHc(iFc,2) = sqrt(  ( centroid(m%fcCv(iFc,2),1) - fcX(iFc) )**2
     &      + ( centroid(m%fcCv(iFc,2),1) - fcY(iFc) )**2 ) 

        !fcHt(m%nFc) Lengths along cell face ("tangential" direction)
        g%fcHt(iFc) = sqrt( ( g%vxX(m%fcVx(iFc,1)) - 
     &     g%vxX(m%fcVx(iFc,2)) )**2 + ( g%vxY(m%fcVx(iFc,1)) - 
     &     g%vxY(m%fcVx(iFc,2)) )**2  )


       ! g%fcS = gs(:,:,0:1)
       if (isymm.eq.0) then
          g%fcS(iFc) = g%fcHt(iFc) * width
       elseif (isymm.eq.1 .or. isymm.eq.2) then ! rotational symmetry about crx=0 axis
          ! face surface =  distance between two vertices * 2 pi * x-coordinate of face center
          g%fcS(iFc) = pi*g%fcHt(iFc) * (g%vxX(m%fcVx(iFc,1)) +
     &                                   g%vxX(m%fcVx(iFc,2)) )
       elseif (isymm.eq.3 .or. isymm.eq.4) then ! rotational symmetry about cry=0 axis
          ! face surface =  distance between two vertices * 2 pi * y-coordinate of face center
          g%fcS(iFc) = pi*g%fcHt(iFc) * (g%vxY(m%fcVx(iFc,1)) +
     &                                   g%vxY(m%fcVx(iFc,2)) )

       endif


       !Calculate fcEb, unit vector of magnetic field at face center
       ! calculated as the weighted average vector of cvEb of the neighboring cells
       ! using the distance between cell center and face center fcHc as weighting factors
       ! This is called inverse distance interpolation
       g%fcEb(iFc,0) = (g%fcHc(iFc,2) * g%cvEb(m%fcCv(iFc,1),0) + 
     &           g%fcHc(iFc,1) * g%cvEb(m%fcCv(iFc,2),0)) /
     &           (g%fcHc(iFc,1) + g%fcHc(iFc,2))
       g%fcEb(iFc,1) = (g%fcHc(iFc,2) * g%cvEb(m%fcCv(iFc,1),1) + 
     &          g%fcHc(iFc,1) * g%cvEb(m%fcCv(iFc,2),1))/ 
     &        (g%fcHc(iFc,1) + g%fcHc(iFc,2))



       !Calculate unit connector econf on face center
       ! is inverse distance average of ecc (unit vector 
       ! between cell center and face center)
       t0 = fcX(iFc) - centroid(m%fcCv(iFc,1),0)  
       t1 = fcY(iFc) - centroid(m%fcCv(iFc,1),1)
       ec1fc(0) = t0 / norm(t0,t1)
       ec1fc(1) = t1 / norm(t0,t1)

       t0 = centroid(m%fcCv(iFc,2),0) - fcX(iFc)
       t1 = centroid(m%fcCv(iFc,2),1) - fcY(iFc)
       efcc2(0) = t0 / norm(t0,t1)
       efcc2(1) = t1 / norm(t0,t1)
       
       econf(0) = (g%fcHc(iFc,2)*ec1fc(0) + g%fcHc(iFc,1)*efcc2(0)) 
     &             /  (g%fcHc(iFc,1) + g%fcHc(iFc,2))
       econf(1) = (g%fcHc(iFc,2)*ec1fc(1) + g%fcHc(iFc,1)*efcc2(1)) 
     &            /  (g%fcHc(iFc,1) + g%fcHc(iFc,2))

       !Calculate face normal ! NOG NIET CORRECT
       ! algemeen : dx_normaal = dy_verts  (= y2 - y1)
       !            dy_normaal = -dx_verts (= x1 - x2)
       t0 = g%vxY(m%fcVx(iFc,2)) - g%vxY(m%fcVx(iFc,1))
       t1 = - ( g%vxX(m%fcVx(iFc,2)) - g%vxX(m%fcVx(iFc,1)) )
       fn(0) = t0 /norm(t0,t1)
       fn(1) = t1 /norm(t0,t1)     

       
*    .. compute angles
       !gamma - angle between cell connector (on face) and the face normal
       g%fcQgam(iFc,0) = cosang(econf(0), econf(1),fn(0), fn(1) )
       g%fcQgam(iFc,1) = sinang(econf(0), econf(1),fn(0), fn(1) )
       
       !alpha - angle between poloidal field vector on face and 
       g%fcQalf(iFc,0) = 
     &      cosang(g%fcEb(iFc,0),g%fcEb(iFc,1),fn(0),fn(1))
       g%fcQalf(iFc,0) = 
     &      sinang(g%fcEb(iFc,0),g%fcEb(iFc,0),fn(0),fn(1))

       ! beta - angle between connector and poloidal magnetic field vector
       ! or is equal to gamma minus alpha
       ! use sin(x+-y) = sin(x)*cos(y) +- cos(x)*sin(y)
       !     cos(x+-y) = cos(x)*cos(y) -+ sin(x)*sin(y)

!        g%fcQbet(iFc,0) = cosang(econf(0),econf(1),
!     &                  g%fcEb(iFc,0) , g%fcEb(iFc,1))
!       g%fcQbet(iFc,1) = sinang(econf(0),econf(1),
!     &                  g%fcEb(iFc,0) , g%fcEb(iFc,1))

       g%fcQbet(iFc,0) = g%fcQgam(iFc,0)*g%fcQalf(iFc,0) + 
     &                   g%fcQgam(iFc,1)*g%fcQalf(iFc,1)  
       g%fcQbet(iFc,1) = g%fcQgam(ifc,1)*g%fcQalf(ifc,0) - 
     &                   g%fcQgam(ifc,0)*g%fcQalf(ifc,1)     
            
       
*    .. correction of angles for negative toroidal field  
       g%fcQalf(iFc,1) = g%fcQalf(iFc,1)*sign(1.0_R8,g%fcBb(iFc,2))
       g%fcQbet(iFc,0) = g%fcQbet(iFc,0)*sign(1.0_R8,g%fcBb(iFc,2))       
        

*    .. correct angles in special cases (aligned, poloidal,...) 
       if (fcType(iFc).eq.1) then ! radial face
           g%fcQalf(iFc,0) = 0.0_R8
           !g%fcQalf(iFc,1) = +- 1.0_R8 (iets met sign van oorspronkelijke hoek ofzo?)
           g%fcQalf(iFc,1) = sign(1.0_R8,g%fcQalf(iFc,1))

       elseif (fcType(iFc).eq.2) then ! poloidal face
           !beta = 0
           !g%fcQbet(iFc,0) = +-1.0_R8 !cos is +-1
           g%fcQbet(iFc,0) = sign(1.0_R8,g%fcQbet(iFc,0))
           g%fcQbet(iFc,1) = 0.0_R8 ! sin is 0

       endif




*    .. angles checks: cos2 + sin2 = 1
        tol_angles = 1e-5_R8
        if (( g%fcQalf(iFc,0)**2 + g%fcQalf(iFc,1)**2 - 1) 
     &         .gt. tol_angles) then
              write(*,*) "b2agmt: fcQalf not conform at", iFc
        endif
        if (( g%fcQgam(iFc,0)**2 + g%fcQgam(iFc,1)**2 - 1) 
     &         .gt. tol_angles) then
              write(*,*) "b2agmt: fcQgam not conform at", iFc
        endif
        if (( g%fcQbet(iFc,0)**2 + g%fcQbet(iFc,1)**2 - 1) 
     &         .gt. tol_angles) then
              write(*,*) "b2agmt: fcQbet not conform at", iFc
        endif       
        

       ! all angles calculated 

       !Calculate fcPbs and fcPbshz - pitch * magnetic surface
       !structured code
       !left faces


        !unstructured code
        g%fcPbs(iFc) = g%fcS(iFc) * g%fcQalf(iFc,0) * 
     &   abs(g%fcBb(iFc,0))/g%fcBb(iFc,3)

      



      enddo
    


*    ..test sign of g%cvVol, g%cvS, g%fcS
       ! based on b2xvs1





* ..return
      call subend ()
      return

*-----------------------------------------------------------------------
*.end b2agmt

      end subroutine b2agmt


*-----------------------------------------------------------------------
      subroutine calc_psi_centroid(g,m,centroid,iCv)
      use b2mod_types
      use b2us_map
      use b2us_geo
      implicit none

      type(mapping) :: m
      type(geometry) :: g
      real(kind=R8) ::  centroid(m%nCv,0:1), psic, isx(0:1), isy(0:1),
     * t0
      integer iCv, verts(0:3), faces(0:3), i1, i2, ii, iF

cnv   General approach similar to calculation of cell center of triangles
cnv   but now done for quads
      
      !Identify vertices with min and max psi values
      verts = m%cvVx(m%cvVxP(iCv,1):m%cvVxP(iCv,1)+m%cvVxP(iCv,2)-1)
      faces = m%cvFc(m%cvFcP(iCv,1):m%cvFcP(iCv,1)+m%cvFcP(iCv,2)-1)
      
      i1 = minloc(g%vxFpsi(verts),1) -1
      i2 = maxloc(g%vxFpsi(verts),1) -1
      psic = 0.5_R8*(g%vxFpsi(verts(i1)) + g%vxFpsi(verts(i2)))
      
      ! Loop over faces to find which ones are intersected by the psic surface
      ii = 0      
      do iF = 0,3 
      if (psic.ge.min(g%vxFpsi(m%fcVx(faces(iF),1)) , 
     &             g%vxFpsi(m%fcVx(faces(iF),2))) .and.
     &   (psic.lt.max(g%vxFpsi(m%fcVx(faces(iF),1)),  
     &             g%vxFpsi(m%fcVx(faces(iF),2)))) ) then !intersected  
        t0 = (psic - g%vxFpsi(m%fcVx(faces(iF),1))) / 
     &        (g%vxFpsi(m%fcVx(faces(iF),2)) - 
     &         g%vxFpsi(m%fcVx(faces(iF),1)))
        isx(ii) = g%vxX(m%fcVx(faces(iF),1)) + t0 * 
     &  ( g%vxX(m%fcVx(faces(iF),2)) - g%vxX(m%fcVx(faces(iF),1)))
        isy(ii) = g%vxY(m%fcVx(faces(iF),1)) + t0 * 
     &  ( g%vxY(m%fcVx(faces(iF),2)) - g%vxY(m%fcVx(faces(iF),1)))     
        ii = ii+1
      endif
      enddo  
      if (ii.ne.2) then
         write(*,*) 'Problem with intersections in quad', iCv
         write(*,*) 'Number of intersections found: ' , ii
         call xerrab('Intersections in quad not found!')
      endif         
      
      !Compute centroid
      centroid(iCv,0) = 0.5_R8*sum(isx)
      centroid(iCv,1) = 0.5_R8*sum(isy)



      end subroutine calc_psi_centroid

*-----------------------------------------------------------------------
      subroutine order_face_inf(g,m)
      use b2mod_types
      use b2us_map
      use b2us_geo
      implicit none

      type(mapping) :: m
      type(geometry) :: g

      integer verts(1:2), cells (1:2), iFc
      real(kind=R8) ::  fcX, fcY, vec_verts(0:1), vec_cell(0:1), sin

      do iFc = 1,m%nFc
        verts = m%fcVx(iFc,:)
        cells = m%fcCv(iFc,:)

        !vector from V1 to V2
        vec_verts(0) = g%vxX(verts(2)) - g%vxX(verts(1)) 
        vec_verts(1) = g%vxY(verts(2)) - g%vxY(verts(1))

        !vector from geometric face center to cell center of c1
        fcX = 0.5_R8*(g%vxX(verts(2)) + g%vxX(verts(1)))
        fcY = 0.5_R8*(g%vxY(verts(2)) + g%vxY(verts(1)))

        vec_cell(0) = g%cvX(cells(1)) - fcX
        vec_cell(1) = g%cvY(cells(1)) - fcY
        
        !Calculate sinus of vec_verts and vec_cell
        ! sin = | a x b | / norm(a) * norm(b)
        ! with | a x b | = ax*by - bx*ay
        !only the sign of sin matters to no norms are used in the denominator
        sin = vec_verts(0) * vec_cell(1) - vec_cell(0) * vec_verts(1)

        if (sin.lt.0) then  !reordering needed: change V1 and V2
           m%fcVx(iFc,1) = verts(2)
           m%fcVx(iFc,2) = verts(1)
        endif 


      enddo


      end subroutine order_face_inf


*-----------------------------------------------------------------------

      subroutine build_fcType(g,m,fcType,tol)
      use b2mod_types
      use b2us_map
      use b2us_geo
      use b2us_cellhelper
      implicit none

      type(mapping) :: m
      type(geometry) :: g  
      integer :: fcType(m%nFc), iFc
      real(kind=R8) :: tol

      !tol = 1.0e-4_R8

      do iFc = 1,m%nFc
          if (isAligned(iFc,m,g,tol)) then
             fcType(iFc) = 1
          elseif (isPoloidal(iFc,m,g,tol)) then
             fcType(iFc) = 2
          else 
             fcType(iFc) = 0
          endif
      enddo


      end subroutine build_fcType     


