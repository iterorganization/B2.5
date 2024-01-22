*-----------------------------------------------------------------------
*.specification

      subroutine b2agmt (g, m, isymm,artificial_slab)
      use b2mod_types
      use b2mod_constants
      use tradui_constants
      use carre_constants
      use b2us_geo
      use b2us_map
      use b2mod_math
      use b2mod_subsys
      implicit none
*   ..input arguments (unchanged on exit)
      integer isymm
      type(geometry) :: g
      type(mapping) :: m

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
*     This routine takes as input the coordinates and potential 
*     of the all vertices in the mesh, the potential field and all
*     connectivity data of the mesh. It computes metric coefficients of
*     the mesh: cvVol, cvSz, cvBb (of guard cell only, for other cell it
*     was already determined in b2agbb), fcQalf, fcQbet, fcQgam at cell
*     faces.
*
*
*  5. parameters (see also routine b2cdcv)
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
*     cvX, cvY - (1:nCv) real*8 array, input.
*     See the calling routine for description.
*
*     cvVol - (1:nCv) real*8 array.
*     cvVol(iCv) specifies the volume of the cell iCv. (In the 
*     case of toroidal symmetry this volume
*     includes a factor 2*pi*cvX, while in the case of cylindrical
*     symmetry it includes a factor 2*pi*cvY.)
*     It will hold that 0.lt.cvVol(,).
*
*     fcQgam - (1:nFc,0:1) real*8 array, output.
*     For iFc in (1:nFc), fcQgam(iFc,0) specifies cos(t) and 
*     fcQgam(iFc,1) specifies sin(t), where t is the included angle
*     between the the cell-connector vector and the normal vector
*     of face iFc. It will hold that 0.lt.fcQgam(,).le.1.
*
*     fcQalf - (1:nFc,0:1) real*8 array.
*     For iFc in (1:nFc), fcQalf(iFc,0) specifies cos(t) and
*     fcQalf(iFc,1) specifies sin(t), where t is the include angle
*     between the magnetic field vector and the face normal.
*
*     fcQbet - (1:nFc,0:1) real*8 array.
*     For iFc in (1:nFc), fcQbet(iFc,0) specifies cos(t) and
*     fcQbet(iFc,1) specifies sin(t), where t is the include angle
*     between the cell-connector and the magnetic field vector.
*   
*     fcHc - (1:nFc,1:2) real*8 array.
*     For iFc in (1:nFc), fcHc(iFc,:) are the distances from
*     the cell centers fcCv(iFc,:) to the face center of face iFc.
*
*     fcHt - (1:nFc) real*8 array.
*     For iFc in (1:nFc), fcHt is the face length of face iFc.
*     
*   
*-----------------------------------------------------------------------
*.declarations

*   ..local variables
      integer i, j, iCv, iFc, iF, face, gc, bc, v1, v2, iVx,
     * artificial_slab 
      integer verts3(0:2), faces3(0:2), cells(1:2),
     * fcType(m%nFc), c1, c2,  verts2(0:1), intersectF(0:1), 
     * faces4_1(0:3), faces4_2(0:3)
      integer, allocatable :: verts(:)
      integer i1, i2, ii, i1_1, i1_2, i2_1, i2_2, nv
      integer volhmk, psi_flag
      real (kind=R8) ::
     *   p0x, p0y, p1x, p1y, width, centroid(m%nCv,0:1), bcFpsi,
     *   t0, t1, psic, dx, dy, cf(0:1), isx(0:1), isy(0:1),
     *   guard_rescale, fvec(0:1), dg, fcPsi, ec1fc(0:1), efcc2(0:1),
     *   econf(0:1), fn(0:1), tol_angles, cvFpsi(1:m%nCv),
     *   fcH, fcW, fcX(1:m%nFc), fcY(1:m%nFc), bcBb(0:1), 
     *   fcQgam1(1:m%nFc,0:1), fcQgam2(1:m%nFc,0:1)
      integer, save :: style = 0
*   ..procedures
      intrinsic min, max, sqrt
      real (kind=R8) ::
     &  sin_t, cos_t, sin_tl, cos_tl, tt
      external xertst, ipgetr, ipgeti
      external b2xvsg
*   ..initialisation
      data volhmk/1/, width /1.0_R8/, guard_rescale /1.0e-1_R8/
      save volhmk, width, guard_rescale


*-----------------------------------------------------------------------
*.documentation-internal
*
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

! loop over cells
      do iCv = 1,m%nCv

      !Default values
      !Compute qz which is cvQgam
      g%cvQgam(iCv,0) = 1.0_R8 ! default value to get through tests of solver
      g%cvQgam(iCv,1) = 0.0_R8 ! default value to get through tests of solver

       
       select case (m%cvFcP(iCv,2))

       case(3,4,5)  !Triangles, quadrilateral and pentagonal

       nv = m%cvVxP(iCv,2)
       allocate(verts(0:nv-1))

       !ebc unit vector of poloidal field in cell center
       verts = m%cvVx(m%cvVxP(iCv,1):m%cvVxP(iCv,1)+nv-1 )
       ! b_y
       t0 = sum(g%vxBb(verts,1)) / nv
       ! b_x
       t1 = - sum(g%vxBb(verts,0)) / nv

       g%cvEb(iCv,0) = t0 / norm(t0,t1)
       g%cvEb(iCv,1) = t1 / norm(t0,t1)

*    ..unit vector of total magnetic field, defined in cell center
      ! until this point, cvEb was unit vector of poloidal field
      ! now, correct for toroidal field
       if (g%cvBb(iCv,3) == 0.0_R8) then
         g%cvEb(iCv,:) = 0.0_R8
       else
         g%cvEb(iCv,0) = 
     &       abs(g%cvBb(iCv,0))/g%cvBb(iCv,3)*g%cvEb(iCv,0)
         g%cvEb(iCv,1) = 
     &       abs(g%cvBb(iCv,0))/g%cvBb(iCv,3)*g%cvEb(iCv,1)
         g%cvEb(iCv,2) = g%cvBb(iCv,2)/g%cvBb(iCv,3)
       endif 

*    ..compute volume

      ! we assume the vertices are ordered 
      ! Because of this we can employ the quadArea for quads
      ! which assumes the vertices in the correct sequence - right turning
      if (nv.eq.4) then
       g%cvSz(iCv) = quadArea(g%vxX(verts(0)), g%vxY(verts(0)),  
     &     g%vxX(verts(1)), g%vxY(verts(1)), g%vxX(verts(2)),
     &     g%vxY(verts(2)), g%vxX(verts(3)), g%vxY(verts(3)) )

       centroid(iCv,:) = quadCentroid(g%vxX(verts(0)), 
     &      g%vxY(verts(0)),g%vxX(verts(1)), g%vxY(verts(1)), 
     &      g%vxX(verts(3)),g%vxY(verts(3)), g%vxX(verts(2)),
     &      g%vxY(verts(2)))

      elseif (nv.eq.5) then
       g%cvSz(iCv) = pentArea(g%vxX(verts(0)), g%vxY(verts(0)),  
     &     g%vxX(verts(1)), g%vxY(verts(1)), g%vxX(verts(2)),
     &     g%vxY(verts(2)), g%vxX(verts(3)), g%vxY(verts(3)), 
     &     g%vxX(verts(4)), g%vxY(verts(4)) )

       centroid(iCv,:) = pentCentroid(g%vxX(verts(0)), g%vxY(verts(0)),  
     &     g%vxX(verts(1)), g%vxY(verts(1)), g%vxX(verts(2)),
     &     g%vxY(verts(2)), g%vxX(verts(3)), g%vxY(verts(3)), 
     &     g%vxX(verts(4)), g%vxY(verts(4)) )

      elseif (nv.eq.3) then
       g%cvSz(iCv) = TriangleArea(g%vxX(verts(0)), g%vxY(verts(0)),
     &                        g%vxX(verts(1)), g%vxY(verts(1)),
     &                        g%vxX(verts(2)), g%vxY(verts(2)))
       centroid(iCv,:) = triangleCentroid(g%vxX(verts(0)), 
     &                        g%vxY(verts(0)),
     &                        g%vxX(verts(1)), g%vxY(verts(1)),
     &                        g%vxX(verts(2)), g%vxY(verts(2)))
      else 
         call xerrab('Area computation not implemented for this type')
      endif

      if  (isymm.eq.0) then  !slab geometry
      ! volume is area times width
        g%cvVol(iCv) = g%cvSz(iCv) *width 
        call calc_psi_cellcenter(g,m,iCv,psic,isymm)
        cvFpsi(iCv) = psic          
      elseif (isymm.eq.1 .or. isymm.eq.2) then ! rotational symmertry about x=0 axis
        if (volhmk.eq.0) then 
          g%cvX(iCv) = sum(g%vxX(verts)) / nv
          g%cvY(iCv) = sum(g%vxY(verts)) / nv
          cvFpsi(iCv) = sum(g%vxFpsi(verts)) / nv
        else
          ! Calculate the cellcenter similar to the method for the triangle cellcenter
          call calc_psi_cellcenter(g,m,iCv,psic,isymm)
          cvFpsi(iCv) = psic   

        endif
        ! volume is area time circumference of circle with radius equal to the x-coordinate of centroid 
        g%cvVol(iCv) = g%cvSz(iCv)*2.0_R8*pi*centroid(iCv,0)
      elseif  (isymm.eq.3 .or. isymm.eq.4) then ! rotational symmetry about y=0
        if (volhmk.eq.0) then
          g%cvX(iCv) = sum(g%vxX(verts)) / nv
          g%cvY(iCv) = sum(g%vxY(verts)) / nv
          cvFpsi(iCv) = sum(g%vxFpsi(verts)) / nv
        else 
          ! Calculate the cellcenter similar to the method for the triangle cellcenter
          call calc_psi_cellcenter(g,m,iCv,psic,isymm)
          cvFpsi(iCv) = psic

        endif
        ! volume is area time circumference of circle with radius equal to the y-coordinate of centroid 
        g%cvVol(iCv) = g%cvSz(iCv)*2.0_R8*pi*centroid(iCv,1)       

      endif
       
       deallocate(verts)


       case(1)  ! Guard Cells (are at the end of the cells array so can use information 
                ! on other cells calculated above)
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

          ! set cvReg of guard cell equal to cvReg of boundary cell
          m%cvReg(gc) = m%cvReg(bc)

          ! set bcFpsi equal to cvFpsi value 
          bcFpsi = cvFpsi(bc)

          ! Determine cell center of guard cell
          ! If it is possible to put the cell center of the guard cell on the same psi
          ! as the boundary cell, this is done. Otherwise the cell center is placed
          ! in the extension of the vector from the boundary cell center and the face center.
          if (abs(fcPsi - bcFpsi)/abs((g%vxFpsi(v2)-g%vxFpsi(v1)))
     &       .lt. 1e-1_R8) then !get the guard cell on the same fpsi
            !Calculate intersection point - find point on face where psi = fcPsi
            t0 = (bcFpsi - g%vxFpsi(v1))/(g%vxFpsi(v2)-g%vxFpsi(v1))      
            isx(1) =  g%vxX(v1) + t0 *(g%vxX(v2)- g%vxX(v1))
            isy(1) =  g%vxY(v1) + t0 *(g%vxY(v2)- g%vxY(v1))
            
            !Set psi value
            psi_flag = 1       
            cvFpsi(iCv) = bcFpsi
                   
            
          else ! use geometric face center
            isx(1) = 0.5_R8 * (g%vxX(v1) + g%vxX(v2))
            isy(1) = 0.5_R8 * (g%vxY(v1) + g%vxY(v2))

            psi_flag = 0
          
          endif

          
          ! Vector from cell center to face center
          cf(0) = isx(1) - g%cvX(bc)
          cf(1) = isy(1) - g%cvY(bc) 

          ! Determine length of face (height of cell)
          fcH = sqrt( (g%vxX(v2)- g%vxX(v1))**2 + 
     &                (g%vxY(v2)- g%vxY(v1))**2 )

          ! Determine width of cell along fcPsi (width of cell)
          fcW = 2.0_R8 * sqrt(cf(0)**2 + cf(1)**2)
          
          ! Determine width (dx) and height (dy) of the guard cell depending on the
          ! width and height of the boundary cell
          if (fcH .lt. fcW) then
             dx = 2.0_R8*cf(0)* guard_rescale * fcH/fcW 
             dy = 2.0_R8*cf(1)* guard_rescale * fcH/fcW 
          else
             dx = 2.0_R8*cf(0) * guard_rescale
             dy = 2.0_R8*cf(1) * guard_rescale
          endif          
        
          dg = sqrt(dx**2+dy**2)
          g%cvHx(iCv) = dg * 2
            
          ! cell area is dg time face Area (times cosine of dg and face normal, but neglected here)
          g%cvSz(iCv) = dg * fcH        

          ! Cell center of guard is face center or intersection point plus (dx,dy)
          ! Vector from face center is direction of (dx/2, dy/2)
          ! (Because cell center of guard cell should be in the middle of the cell)
          ! For a guard cell, the cellcenter is put equal to the centroid
          centroid(gc,0) = isx(1) + 0.5_R8 * dx
          centroid(gc,1) = isy(1) + 0.5_R8 * dy 

          g%cvX(iCv) = centroid(iCv,0)  ! gc is equal to icv
          g%cvY(iCv) = centroid(iCv,1)  

         !Correct the magnetic field of guard cell 
         if (artificial_slab.eq.0) then
          call recalc_guard_bb(m,g,isymm,bc,iCv,v1,v2)
         else
          g%cvBb(iCv,0) = g%cvBb(1,0)
          g%cvBb(iCv,1) = g%cvBb(1,1)
          g%cvBb(iCv,2) = g%cvBb(1,2)
          g%cvBb(iCv,3) = g%cvBb(1,3)
         endif


          !Compute ebc - same method as quad - NEEDS REVIEW
          !ebc unit vector of poloidal field in cell center
          verts2 = m%cvVx(m%cvVxP(iCv,1):m%cvVxP(iCv,1) + 
     &          m%cvVxP(iCv,2)-1 )
          ! b_y
          t0 = ( g%vxBb(verts2(0),1) + g%vxBb(verts2(1),1)    
     &           ) / 2.0_R8
          ! b_x
          t1 = - ( g%vxBb(verts2(0),0) + g%vxBb(verts2(1),0) 
     &      ) / 2.0_R8

          g%cvEb(iCv,0) = t0 / norm(t0,t1)
          g%cvEb(iCv,1) = t1 / norm(t0,t1)

*    ..unit vector of total magnetic field, defined in cell center
          ! until this point, cvEb was unit vector of poloidal field
          ! now, correct for toroidal field
          if (g%cvBb(iCv,3) == 0.0_R8) then
            g%cvEb(iCv,:) = 0.0_R8
          else
            g%cvEb(iCv,0) = 
     &          abs(g%cvBb(iCv,0))/g%cvBb(iCv,3)*g%cvEb(iCv,0)
            g%cvEb(iCv,1) = 
     &          abs(g%cvBb(iCv,0))/g%cvBb(iCv,3)*g%cvEb(iCv,1) 
            g%cvEb(iCv,2) = g%cvBb(iCv,2)/g%cvBb(iCv,3)
          endif

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

       case default
          call xerrab('b2agmt: Type of cell not implemented')

       end select
       
 
      enddo  !cell loop

* ..compute fcQalf, fcQbet, fcQgam at cell faces

* ..For the angle calculation, the connectivity data fcVx and fcCv need to comply
*  to a right turning scheme. For cell fcCv(iFc,1), the vertices in fcVx need to
*  ordered right turning. Only here (in the code) this is possible to perform 
*  as the guard cells, only now have a proper cell center. When the cell center 
*  coincides with the face center (as set as dummy cell-center for the guard cell
*  determined in b2ag_ghostcells), this would complicate the ordering significantly.

      call order_face_inf(g,m)

* ..Also for angle calculation, it is necessary to distinguish between radial, poloidal
*  and other faces, indicated by the array fcType. fcType.EQ.1 indicates a radial face,
*  fcType.EQ.2 indicates a poloidal face. A face is a poloidal face when it is not aligned
*  with the magnetic field and when the cell centers of its cells lay on the same psi.      

      call build_fcType(g,m,fcType,cvFpsi)

! loop over faces 
      !pre compute face centers
      do iFc = 1,m%nFc
         v1 = m%fcVx(iFc,1)
         v2 = m%fcVx(iFc,2)

        ! compute face centers
        fcX(iFc) = (g%vxX(v1) + g%vxX(v2) )/2.0_R8
        fcY(iFc) = (g%vxY(v1) + g%vxY(v2) )/2.0_R8   

      enddo
      do iFc = 1, m%nFc
        v1 = m%fcVx(iFc,1)
        v2 = m%fcVx(iFc,2)
        c1 = m%fcCv(iFc,1)
        c2 = m%fcCv(iFc,2)

        ! fcHc(m%nFc,2) : Lengths between cell face center and neighboring cell centers ("connector" direction)
        g%fcHc(iFc,1) = sqrt(  ( g%cvX(c1) - fcX(iFc) )**2
     &      + ( g%cvY(c1) - fcY(iFc) )**2 ) 
        g%fcHc(iFc,2) = sqrt(  ( g%cvX(c2) - fcX(iFc) )**2
     &      + ( g%cvY(c2) - fcY(iFc) )**2 ) 

        !fcHt(m%nFc) Lengths along cell face ("tangential" direction)
        g%fcHt(iFc) = sqrt( ( g%vxX(v1) - g%vxX(v2) )**2 + 
     &                    ( g%vxY(v1) -  g%vxY(v2) )**2  )

       ! g%fcS : face surface
       if (isymm.eq.0) then
          g%fcS(iFc) = g%fcHt(iFc) * width
       elseif (isymm.eq.1 .or. isymm.eq.2) then ! rotational symmetry about crx=0 axis
          ! face surface =  face length * 2 pi * mean x-coordinate of face center
          g%fcS(iFc) = pi*g%fcHt(iFc) * (g%vxX(v1) + g%vxX(v2) )
       elseif (isymm.eq.3 .or. isymm.eq.4) then ! rotational symmetry about cry=0 axis
          ! face surface =  face length * 2 pi * mean y-coordinate of face center
          g%fcS(iFc) = pi*g%fcHt(iFc) * (g%vxY(v1) + g%vxY(v2) )
       endif


       ! Calculate fcEb, unit vector of magnetic field at face center
       ! calculated as the weighted average vector of cvEb of the neighboring cells
       ! using the distance between cell center and face center fcHc as weighting factors
       ! This is called inverse distance interpolation
       !OLD METHOD
       !g%fcEb(iFc,0) = (g%fcHc(iFc,2) * g%cvEb(c1,0) + 
       !&           g%fcHc(iFc,1) * g%cvEb(c2,0)) /
       !&           (g%fcHc(iFc,1) + g%fcHc(iFc,2))
       !g%fcEb(iFc,1) = (g%fcHc(iFc,2) * g%cvEb(c1,1) + 
       !&          g%fcHc(iFc,1) * g%cvEb(c2,1))/ 
       !&        (g%fcHc(iFc,1) + g%fcHc(iFc,2))
       
       ! similar to linear interpolating Bb to the face centers
       ! b_y
       t0 = 0.5_R8 * ( g%vxBb(v1,1) + g%vxBb(v2,1) )
       ! b_x
       t1 = - 0.5_R8 * ( g%vxBb(v1,0) + g%vxBb(v2,0))

       g%fcEb(iFc,0) = t0 / norm(t0,t1)
       g%fcEb(iFc,1) = t1 / norm(t0,t1)



       ! Calculate unit connector econf on face center
       ! is inverse distance average of ecc (unit vector 
       ! between cell center and face center)
       t0 = fcX(iFc) - g%cvX(m%fcCv(iFc,1))  
       t1 = fcY(iFc) - g%cvY(m%fcCv(iFc,1))
       ec1fc(0) = t0 / norm(t0,t1)
       ec1fc(1) = t1 / norm(t0,t1)

       t0 = g%cvX(m%fcCv(iFc,2)) - fcX(iFc)
       t1 = g%cvY(m%fcCv(iFc,2)) - fcY(iFc)
       efcc2(0) = t0 / norm(t0,t1)
       efcc2(1) = t1 / norm(t0,t1)
       
       econf(0) = (g%fcHc(iFc,2)*ec1fc(0) + g%fcHc(iFc,1)*efcc2(0)) 
     &             /  (g%fcHc(iFc,1) + g%fcHc(iFc,2))
       econf(1) = (g%fcHc(iFc,2)*ec1fc(1) + g%fcHc(iFc,1)*efcc2(1)) 
     &            /  (g%fcHc(iFc,1) + g%fcHc(iFc,2))

       ! Corrections for specific geometries
       ! poloidal faces
       ! =  is a face across the poloidal direction in a cell with aligned radial faces
       if (fcType(iFc).eq.2) then 
          
          t0 = econf(0)*g%fcEb(iFc,0) + econf(1)*g%fcEb(iFc,1)
          if (t0.ge.0.0_R8) then ! econf in same direction as fcEb
              econf(0) = g%fcEb(iFc,0)
              econf(1) = g%fcEb(iFc,1)
              ec1fc(0) = g%fcEb(iFc,0)
              ec1fc(1) = g%fcEb(iFc,1)
              efcc2(0) = g%fcEb(iFc,0)
              efcc2(1) = g%fcEb(iFc,1)              
          else 
              econf(0) = -g%fcEb(iFc,0)
              econf(1) = -g%fcEb(iFc,1)
              ec1fc(0) = -g%fcEb(iFc,0)
              ec1fc(1) = -g%fcEb(iFc,1)
              efcc2(0) = -g%fcEb(iFc,0)
              efcc2(1) = -g%fcEb(iFc,1)              
          endif
    
       endif

       ! Calculate face normal (not normalized!)
       ! general: dx_normal = dy_verts  (= y2 - y1)
       !          dy_normal = -dx_verts (= x1 - x2)
       t0 = g%vxY(m%fcVx(iFc,2)) - g%vxY(m%fcVx(iFc,1))
       t1 = - ( g%vxX(m%fcVx(iFc,2)) - g%vxX(m%fcVx(iFc,1)) )
       fn(0) = t0 
       fn(1) = t1    

       
*    .. compute angles
       ! gamma - angle between cell connector (on face) and the face normal
       ! NEW IMPLEMENTATION 28/10: to mitigate numerical error on fluxes due to 
       ! cell connector approximation (inverse distance averaging)
       ! gamma is also computed with distance averaging between 
       ! gamma1 and gamma2 which are gamma's calculated between on the 
       ! unapproximated connector between cell center and face center.
       ! With this method, numerical error on the 5-point contribution on 
       ! gradient calculation is significantly reduced for general grids

       fcQgam1(iFc,0) = sinang(ec1fc(0),ec1fc(1),fn(0),fn(1))
       fcQgam1(iFc,1) = cosang(ec1fc(0),ec1fc(1),fn(0),fn(1))

       fcQgam2(iFc,0) = sinang(efcc2(0),efcc2(1),fn(0),fn(1))
       fcQgam2(iFc,1) = cosang(efcc2(0),efcc2(1),fn(0),fn(1))

       g%fcQgam(iFc,0) = ( fcQgam1(iFc,0)*g%fcHc(iFc,1) + 
     &    fcQgam2(iFc,0)*g%fcHc(iFc,2)) /
     &    (g%fcHc(iFc,1) + g%fcHc(iFc,2))
       g%fcQgam(iFc,1) = ( fcQgam1(iFc,1)*g%fcHc(iFc,1) + 
     &    fcQgam2(iFc,1)*g%fcHc(iFc,2)) /
     &    (g%fcHc(iFc,1) + g%fcHc(iFc,2))

       ! Correction to comply with sin**2 + cos**2 = 1       
       g%fcQgam(iFc,1) = sign(1.0,g%fcQgam(iFc,1))* 
     &    sqrt(1 - g%fcQgam(iFc,0)**2) 

       !OLD METHOD to compute gamma
       !g%fcQgam(iFc,0) = sinang(econf(0), econf(1),fn(0), fn(1) )
       !g%fcQgam(iFc,1) = cosang(econf(0), econf(1),fn(0), fn(1) )
       
       !alpha - angle between magnetic field vector and normal vector
       g%fcQalf(iFc,0) =
     &         sinang(g%fcEb(iFc,0),g%fcEb(iFc,1),fn(0),fn(1))  
       g%fcQalf(iFc,1) = 
     &         cosang(g%fcEb(iFc,0),g%fcEb(iFc,1),fn(0),fn(1))         


       ! beta - angle between connector and poloidal magnetic field vector
       ! or is equal to gamma minus alpha
       ! use sin(x-y) = sin(x)*cos(y) - cos(x)*sin(y)
       !     cos(x-y) = cos(x)*cos(y) + sin(x)*sin(y)

       g%fcQbet(iFc,0) = g%fcQgam(iFc,0)*g%fcQalf(iFc,0) + 
     &                   g%fcQgam(iFc,1)*g%fcQalf(iFc,1)  
       g%fcQbet(iFc,1) = g%fcQgam(iFc,1)*g%fcQalf(iFc,0) - 
     &                   g%fcQgam(iFc,0)*g%fcQalf(iFc,1)     
            
            

*    .. correct angles in special cases (aligned, poloidal,...) 
       if (fcType(iFc).eq.1) then ! radial face
           g%fcQalf(iFc,0) = 0.0_R8
           g%fcQalf(iFc,1) = sign(1.0_R8,g%fcQalf(iFc,1)) ! rounds to -1 or 1

           !recalculation of fcQbet
           g%fcQbet(iFc,0) = g%fcQgam(iFc,0)*g%fcQalf(iFc,0) + 
     &                   g%fcQgam(iFc,1)*g%fcQalf(iFc,1)  
           g%fcQbet(iFc,1) = g%fcQgam(iFc,1)*g%fcQalf(iFc,0) - 
     &                   g%fcQgam(iFc,0)*g%fcQalf(iFc,1)           

       elseif (fcType(iFc).eq.2) then ! poloidal face
           
           !beta = 0
           !g%fcQbet(iFc,0) = +-1.0_R8 !cos is +-1

           g%fcQbet(iFc,0) = sign(1.0_R8,g%fcQbet(iFc,0))
           g%fcQbet(iFc,1) = 0.0_R8 ! sin is 0

           ! gamma = +- alfa
           g%fcQgam(iFc,0) = abs(g%fcQalf(iFc,0))
           if (g%fcQbet(iFc,0).eq.-1.0_R8) then
           g%fcQgam(iFc,1) = sign(g%fcQalf(iFc,1),g%fcQgam(iFc,1))
           !for sin gamma, it has the value of fcQalf(:,1) and the sign of fcQgam(:,1)
           else
              g%fcQgam(iFc,1) = g%fcQalf(iFc,1)
           endif
       endif

*    .. correction of angles for negative toroidal field  
       g%fcQalf(iFc,1) = g%fcQalf(iFc,1)*sign(1.0_R8,g%fcBb(iFc,2))
       g%fcQbet(iFc,0) = g%fcQbet(iFc,0)*sign(1.0_R8,g%fcBb(iFc,2))


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
        !unstructured code
        g%fcPbs(iFc) = g%fcS(iFc) * g%fcQalf(iFc,0) * 
     &   abs(g%fcBb(iFc,0))/g%fcBb(iFc,3)

      
      enddo


      ! Give vxBb the correct magnetic field value, till now it was
      ! g%vxBb(:,0) = dpsi/dy
      ! g%vxBb(:,1) = -dpsi/dx, for directional purposes

      do iVx = 1,m%nVx
         if (isymm.eq.0) then
            t0 = 1.0_R8
         elseif (isymm.eq.1 .or. isymm.eq.2) then
            t0 = 2.0_R8 * pi * g%vxX(iVx)
         elseif (isymm.eq.3 .or. isymm.eq.4) then
            t0 = 2.0_R8 * pi * g%vxY(iVx)
         endif

         g%vxBb(iVx,0) = sqrt(g%vxBb(iVx,0)**2 + g%vxBb(iVx,1)**2)/t0
         g%vxBb(iVx,1) = 0.0_R8
         g%vxBb(iVx,2) = abs(g%vxFfbz(iVx))/t0
         g%vxBb(iVx,3) = sqrt(g%vxBb(iVx,0)**2 + g%vxBb(iVx,1)**2 + 
     &                      g%vxBb(iVx,2)**2)
         

      enddo
       
*    ..test sign of g%cvVol, g%cvS, g%fcS     
       call b2xvs1(m%nCv, g%cvVol, 1,'cvVol' ,'.gt.')
       call b2xvs1(m%nCv, g%cvSz,   1,'cvS'   ,'.gt.')
       call b2xvs1(m%nFc, g%fcS,   1,'fcS'   ,'.gt.')

* ..return
      call subend ()
      return

*-----------------------------------------------------------------------
*.end b2agmt

      end subroutine b2agmt


*-----------------------------------------------------------------------
      subroutine calc_psi_cellcenter(g,m,iCv,psic,isymm)
      use b2mod_types
      use b2us_map
      use b2us_geo
      implicit none

      type(mapping) :: m
      type(geometry) :: g
      real(kind=R8) ::  psic, isx(0:1), isy(0:1),
     * t0
      integer iCv, i1, i2, ii, iF,
     *  intersectF(0:1), isymm, nv
      integer, allocatable :: verts(:), faces(:)

cnv   General approach similar to calculation of cell center of triangles
cnv   but now for quads

      nv = m%cvVxP(iCv,2)
      allocate(verts(0:nv-1))
      allocate(faces(0:nv-1))
      
      !Identify vertices with min and max psi values
      verts = m%cvVx(m%cvVxP(iCv,1):m%cvVxP(iCv,1)+nv-1)
      faces = m%cvFc(m%cvFcP(iCv,1):m%cvFcP(iCv,1)+nv-1)
      
      i1 = minloc(g%vxFpsi(verts),1) -1
      i2 = maxloc(g%vxFpsi(verts),1) -1
      psic = 0.5_R8*(g%vxFpsi(verts(i1)) + g%vxFpsi(verts(i2)))
      
      ! Loop over faces to find which ones are intersected by the psic surface
      ii = 0      
      do iF = 0,nv-1 
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
        intersectF(ii) = faces(iF)     
        ii = ii+1
      endif
      enddo  
      if (ii.ne.2) then
         write(*,*) 'Problem with intersections in cell', iCv
         write(*,*) 'Number of intersections found: ' , ii
         call xerrab('Intersections in cell not found!')
      endif         
      
      ! Compute cellcenter
      g%cvX(iCv) = 0.5_R8*sum(isx)
      g%cvY(iCv) = 0.5_R8*sum(isy)
      
      ! Compute poloidal length
      g%cvHx(iCv) = sqrt( (isx(1) - isx(0))**2 +  (isy(1) - isy(0))**2 ) 

      ! Correct the magnetic field in the new cell center
      if (isymm.ne.0) then
        call recalc_cvBb(iCv,g,m,isx,isy,intersectF,isymm)
      endif

      ! 

      deallocate(verts)
      deallocate(faces)


      end subroutine calc_psi_cellcenter

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

      subroutine build_fcType(g,m,fcType,cvFpsi)
      use b2mod_types
      use b2us_map
      use b2us_geo
      use b2us_cellhelper
      implicit none

      type(mapping) :: m
      type(geometry) :: g  
      integer :: fcType(m%nFc), iFc
      real(kind=R8) :: cvFpsi(m%nCv)

      do iFc = 1,m%nFc
          if (isAligned(iFc,m)) then
             fcType(iFc) = 1
          elseif (isPoloidal(iFc,m,g,cvFpsi)) then
             fcType(iFc) = 2
          else 
             fcType(iFc) = 0
          endif
      enddo


      end subroutine build_fcType     

*-----------------------------------------------------------------------

      subroutine recalc_guard_bb(m,g,isymm,bc,iCv, v1, v2)
      use b2mod_types
      use b2us_map
      use b2us_geo
      use b2mod_constants
      implicit none

      type(mapping) :: m
      type(geometry) :: g
      integer isymm, bc, iCv, v1, v2, nv
      real(kind=R8) :: vec1(0:2), vec2(0:2), p(0:2), d, bcBb(0:2), 
     &   gcBb(0:2), t0
      integer, allocatable :: verts(:)

      !1. interpolate vxBb(0:1) to the boundary cell center
      nv = m%cvVxP(bc,2)
      allocate(verts(nv))

      verts = m%cvVx(m%cvVxP(bc,1):m%cvVxP(bc,1)+nv-1)
      bcBb(0) = 0.25_R8 * sum( g%vxBb(verts,0) )
      bcBb(1) = 0.25_R8 * sum( g%vxBb(verts,1) )
      bcBb(2) = 0.25_R8 * sum( g%vxFfbz(verts))


      !2. construct 3d-plane with two guard cell vertices and cell center
      ! p(0)*x+ p(1)*y + p(2)*z = d met x = x-coordinate, y = y-coordinate, z = bb-value
      !first the dpsi/dy  = bx value (vxBb(:,0))
      vec1(0) = g%vxX(v1) - g%cvX(bc)
      vec1(1) = g%vxY(v1) - g%cvY(bc)
      vec1(2) = g%vxBb(v1,0) - bcBb(0)
      vec2(0) = g%vxX(v2) - g%cvX(bc)
      vec2(1) = g%vxY(v2) - g%cvY(bc)
      vec2(2) = g%vxBb(v2,0) - bcBb(0)

      !calculate cross to get the normal of the 3d-plane vector p(0:2)
      p(0) = vec1(1)*vec2(2) - vec1(2)*vec2(1) 
      p(1) = vec1(0)*vec2(2) - vec1(2)*vec1(0)
      p(2) = vec1(0)*vec2(1) - vec1(1)*vec2(0)

      !determine 3d-plane with cell center point
      d = p(0)*g%cvX(bc) + p(1)*g%cvY(bc) + p(2)*bcBb(0)

      ! fill guard cell center coordinate in into the plane
      gcBb(0) = (d - p(0)*g%cvX(iCv) - p(1)*g%cvY(iCv)) / p(2)

      !same for vxBb(:,1)
      vec1(2) = g%vxBb(v1,1) - bcBb(1)
      vec2(2) = g%vxBb(v2,1) - bcBb(1)

      p(0) = vec1(1)*vec2(2) - vec1(2)*vec2(1) 
      p(1) = vec1(0)*vec2(2) - vec1(2)*vec1(0)
      p(2) = vec1(0)*vec2(1) - vec1(1)*vec2(0)           
      d = p(0)*g%cvX(bc) + p(1)*g%cvY(bc) + p(2)*bcBb(1)

      gcBb(1) = (d - p(0)*g%cvX(iCv) - p(1)*g%cvY(iCv)) / p(2)

      !same for vxFfbz
      vec1(2) = g%vxFfbz(v1) - bcBb(2)
      vec2(2) = g%vxFfbz(v2) - bcBb(2)

      p(0) = vec1(1)*vec2(2) - vec1(2)*vec2(1) 
      p(1) = vec1(0)*vec2(2) - vec1(2)*vec1(0)
      p(2) = vec1(0)*vec2(1) - vec1(1)*vec2(0)           
      d = p(0)*g%cvX(bc) + p(1)*g%cvY(bc) + p(2)*bcBb(2)

      gcBb(2) = (d - p(0)*g%cvX(iCv) - p(1)*g%cvY(iCv)) / p(2)

      !3. compute magnetic field in guard cell center - same method as in b2agbb but now with only one point
      if (isymm.ne.0) then
        if (isymm.eq.1 .or. isymm.eq.2) then
           t0 = 2.0*pi*g%cvX(bc)
        elseif (isymm.eq.3 .or. isymm.eq.4) then
           t0 = 2.0*pi*g%cvY(bc)
        endif
      else 
        t0 = 1.0
      endif

      g%cvBb(iCv,0) = sqrt(gcBb(0)**2 + gcBb(1)**2)/t0
      g%cvBb(iCv,1) = 0.0_R8
      g%cvBb(iCv,2) = gcBb(2)/t0
      g%cvBb(iCv,3) = sqrt( g%cvBb(iCv,0)**2 + 
     &        g%cvBb(iCv,1)**2 + g%cvBb(iCv,2)**2)

      deallocate(verts)


      end subroutine recalc_guard_bb


*-----------------------------------------------------------------------

      subroutine recalc_cvBb(iCv,g,m,isx,isy,intersectF,isymm)
      use b2mod_types
      use b2us_map
      use b2us_geo
      use b2mod_constants
      implicit none

      type(mapping) :: m
      type(geometry) :: g
      integer :: iCv, isymm, v1, v2, v3, v4, intersectF(0:1)
      real(kind=R8) :: isx(0:1), isy(0:1), dist_v1v2, dist_is0v1,
     & dist_is0v2 , dist_is1v3, dist_is1v4, dist_v3v4,
     & isbb(0:1,0:2), t0

            
      !Recompute magnetic field in the new cell-center
      !magnetic field on intersection points
      v1 = m%fcVx(intersectF(0),1)
      v2 = m%fcVx(intersectF(0),2)
      v3 = m%fcVx(intersectF(1),1)
      v4 = m%fcVx(intersectF(1),2)
     
      !magnetic field on intersection point 0
      ! via inverse distance weighting 
      dist_is0v1 = sqrt((isx(0) - g%vxX(v1) )**2 + 
     &                  (isy(0) - g%vxY(v1)   )**2)
      dist_is0v2 = sqrt((isx(0) - g%vxX(v2)  )**2 + 
     &                  (isy(0) - g%vxY(v2)   )**2)
      dist_v1v2 = sqrt((g%vxX(v2) - g%vxX(v1)  )**2 + 
     &                 (g%vxY(v2) - g%vxY(v1)   )**2)
      

      isbb(0,0) = (dist_is0v2 * g%vxBb(v1,0) + 
     &              dist_is0v1* g%vxBb(v2,0))/dist_v1v2
      isbb(0,1) = (dist_is0v2 * g%vxBb(v1,1) + 
     &              dist_is0v1* g%vxBb(v2,1))/dist_v1v2
      isbb(0,2) = (dist_is0v2 * g%vxFfBz(v1) + 
     &              dist_is0v1* g%vxFfBz(v2))/dist_v1v2

      !magnetic field at intersection point 1
      dist_is1v3 = sqrt( (isx(1) - g%vxX(v3) )**2 +
     &                   (isy(1) - g%vxY(v3) )**2  )
      dist_is1v4 = sqrt( (isx(1) - g%vxX(v4) )**2 + 
     &                   (isy(1) - g%vxY(v4) )**2  )
      dist_v3v4 = sqrt(   (g%vxX(v3) - g%vxX(v4) )**2 +
     &                    (g%vxY(v3) - g%vxY(v4) )**2  )

      isbb(1,0) = (dist_is1v4 * g%vxBb(v3,0) + 
     &            dist_is1v3 * g%vxBb(v4,0) ) / dist_v3v4
      isbb(1,1) = (dist_is1v4 * g%vxBb(v3,1) + 
     &            dist_is1v3 * g%vxBb(v4,1)) /dist_v3v4
      isbb(1,2) = (dist_is1v4 * g%vxFfBz(v3) + 
     &            dist_is1v3 * g%vxFfBz(v4)) /dist_v3v4  

      if (isymm.eq.0) then
          t0 = 1.0_R8
      elseif (isymm.eq.1 .or. isymm.eq.2) then
          t0 = 2.0_R8 * pi * g%cvX(iCv)
      elseif (isymm.eq.3 .or. isymm.eq.4) then
          t0 = 2.0_R8 * pi * g%cvY(iCv)
      endif    

      !same method to calculate magnetic field as for the faces
      ! here cell center is the mean of the intersection points

      g%cvBb(iCv,0) = 0.5_R8 * 
     &   ( sqrt( isbb(0,0)**2 + isbb(0,1)**2 ) +
     &     sqrt( isbb(1,0)**2 + isbb(1,1)**2) ) / t0 
      g%cvBb(iCv,1) = 0.0_R8
      g%cvBb(iCv,2) = 0.5_R8 * ( isbb(0,2)+isbb(1,2)) /t0
      g%cvBb(iCv,3) = sqrt(g%cvBb(iCv,0)**2 + g%cvBb(iCv,1)**2 +  
     &                  g%cvBb(iCv,2)**2)


      end subroutine recalc_cvBb
