*-----------------------------------------------------------------------
*.specification

      subroutine b2agbb (g,m,isymm)
      !use b2mod_geo
      !use b2mod_indirect
      use b2mod_types
      use b2us_geo
      use b2us_map
      use carre_constants
      implicit none
*   ..input arguments
      integer isymm
      type(mapping) :: m 

*   ..input and output arguments
      type(geometry) :: g

*-----------------------------------------------------------------------
*.documentation
*
*  1. purpose
*     
*     B2AGBB computes the magnetic field from the flux functions.
*
*  3. description
*  
*     This routine computes the three components of the magnetic field
*     and the total magnetic field strength, at the cell centers and 
*     the face centers. The componets (bx,by) are obtained by 
*     differencing the potential fpsi from its values at cell corners.
*     The component bz is obtained directly from the field function ffbz.
*
*
*  5. parameters
*  WIP-NV - see b2agbb_st.F as an example
*-----------------------------------------------------------------------
*.declarations

*   ..local variables
      integer iCv, iFc, vx4(0:3), vx3(0:2), vx2(0:1), vx0, iCvBnd,  
     &  ghostGetBndCell
      real (kind=R8) ::
     *  t0, pi, babs, psdx, psdy, width, hzf,
     *  centroid(0:1) ! 0 = x-coord, 1 = y-coord
*   ..procedures
      intrinsic sqrt, sign

      !data width /1.0_R8/
      !save width

*-----------------------------------------------------------------------
*.computation
 

* ..preliminaries
*   ..subprogram start-up calls
      call subini ('b2agbb')
*   ..test nCv, nFc, nVx
      call xertst (0.le.m%nCv.and.0.le.m%nFc.and.0.le.m%nVx, 
     &   'faulty argument nx, ny')

      width = 1.0_R8
      call ipgetr('b2agmt_1d_width',width)

      pi = 4.*atan(1.)

*   ..compute magnetic field

      do iCv = 1, m%nCv
        if (m%cvLbl(iCv) .ne. GRID_GUARD)  then !guard cells are skipped
        !   ..compute magnetic field at cell center
        if (isymm.ne.0) then
          if (m%cvFcP(iCv,2).eq.4) then !Quad 
           vx4 = m%cvVx(m%cvVxP(iCv,1):m%cvVxP(iCv,1)+m%cvVxP(iCv,2)-1)
           centroid(0) = 0.25_R8*(sum(g%vxX(vx4)))
           centroid(1) = 0.25_R8*(sum(g%vxY(vx4)))
          elseif (m%cvFcP(iCv,2).eq.1) then !Guard cell
           vx2 = m%cvVx(m%cvVxP(iCv,1):m%cvVxP(iCv,1)+m%cvVxP(iCv,2)-1)
           centroid(0) = 0.5_R8*(sum(g%vxX(vx2)))
           centroid(1) = 0.5_R8*(sum(g%vxY(vx2)))
          else ! Triangles
           vx3 = m%cvVx(m%cvVxP(iCv,1):m%cvVxP(iCv,1)+m%cvVxP(iCv,2)-1)
           centroid(0) = 0.3333333_R8*(sum(g%vxX(vx3)))
           centroid(1) = 0.3333333_R8*(sum(g%vxY(vx3)))
          endif

          if (isymm.eq.1 .or. isymm.eq.2) then
           t0 = 2.0*pi*centroid(0)
          elseif (isymm.eq.3 .or. isymm.eq.4) then
           t0 = 2.0*pi*centroid(1)
          endif

        else
          t0 = 1.0
        endif


        if (t0 == 0.0) cycle
          if (m%cvFcP(iCv,2).eq.4) then
           ! structured code
           !bb(ix,iy,0) = & 
           !       &     (fpsi(ix,iy,2)-fpsi(ix,iy,0)+fpsi(ix,iy,3)-fpsi(ix,iy,1))

            vx4 = m%cvVx(m%cvVxP(iCv,1):m%cvVxP(iCv,1)+m%cvVxP(iCv,2)-1)
           ! not sure always correct but only the sign of cvBb 
           ! is important here -ISSUE
           !g%cvBb(iCv,0) = ( g%vxFpsi(vx4(2)) - g%vxFpsi(vx4(0)) 
           !&      +  g%vxFpsi(vx4(3)) -  g%vxFpsi(vx4(1)) )
          
           babs = 0.25_R8*(  
     &      sqrt( g%vxBb(vx4(0),0)**2 + g%vxBb(vx4(0),1)**2 ) + 
     &      sqrt( g%vxBb(vx4(1),0)**2 + g%vxBb(vx4(1),1)**2 ) +      
     &      sqrt( g%vxBb(vx4(2),0)**2 + g%vxBb(vx4(2),1)**2 ) + 
     &      sqrt( g%vxBb(vx4(3),0)**2 + g%vxBb(vx4(3),1)**2 ) )/t0
        
           ! sign takes value of first argument with sign of second argument
           !g%cvBb(iCv,0) = sign(babs,g%cvBb(iCv,0))
           g%cvBb(iCv,0) = babs  !just absolute value
        
          elseif (m%cvFcP(iCv,2).eq.1) then !at the end of array, exploit this
           vx2 = m%cvVx(m%cvVxP(iCv,1):m%cvVxP(iCv,1)+m%cvVxP(iCv,2)-1)
           ! assum sign of cvBb(1) is everywhere the same

           babs = 0.5_R8*(  
     &      sqrt( g%vxBb(vx2(0),0)**2 + g%vxBb(vx2(0),1)**2 ) + 
     &      sqrt( g%vxBb(vx2(1),0)**2 + g%vxBb(vx2(1),1)**2 ) )/t0
           ! cell 3 is high probably a quad (PF region)
           !g%cvBb(iCv,0) = sign(babs,g%cvBb(3,0)) !MAKE THIS GENERAL
           g%cvBb(iCv,0) = babs !just absolute value

          else !Triangle
           vx3 = m%cvVx(m%cvVxP(iCv,1):m%cvVxP(iCv,1)+m%cvVxP(iCv,2)-1)
           ! assum sign of cvBb(1) is everywhere the same

           babs = 0.3333333_R8*( 
     &      sqrt( g%vxBb(vx3(0),0)**2 + g%vxBb(vx3(0),1)**2 ) + 
     &      sqrt( g%vxBb(vx3(1),0)**2 + g%vxBb(vx3(1),1)**2 ) + 
     &      sqrt( g%vxBb(vx3(2),0)**2 + g%vxBb(vx3(2),1)**2 ) )/t0

           !g%cvBb(iCv,0) = sign(babs,g%cvBb(3,0)) !MAKE THIS GENERAL
           g%cvBb(iCv,0) = babs

          endif  
        

          ! operations with single vertex
          vx0 = m%cvVx(m%cvVxP(iCv,1)) !first vertex of the cell
          ! choose here right vertex (same as corrections for
          ! right and upper boundary in structured version)


          !if (isymm.eq.1.or.isymm.eq.2) then
          !  babs = sqrt( g%vxBb(vx0,0)**2 + g%vxBb(vx0,1)**2 )/  
            !&         (2.0*pi*g%vxX(vx0))
          ! elseif (isymm.eq.3.or.isymm.eq.4) then
          !   babs = sqrt( g%vxBb(vx0,0)**2 + g%vxBb(vx0,1)**2 )/ 
          ! &         (2.0*pi*g%vxY(vx0))
          !else
          !  babs=sqrt( g%vxBb(vx0,0)**2 + g%vxBb(vx0,1)**2 )
          !endif


          g%cvBb(iCv,1) = 0.0e0
          if (m%cvFcP(iCv,2).eq.4) then !Quad
            g%cvBb(iCv,2) = 0.25_R8*(sum(g%vxFfbz(vx4)))/t0
          elseif (m%cvFcP(iCv,2).eq.1) then
            g%cvBb(iCv,2) = 0.5_R8*(sum(g%vxFfbz(vx2)))/t0
          else 
            g%cvBb(iCv,2) = 0.3333333_R8 * sum(g%vxFfbz(vx3))/t0
          endif 
        
          g%cvBb(iCv,3) = sqrt( g%cvBb(iCv,0)**2 + 
     &     g%cvBb(iCv,1)**2 + g%cvBb(iCv,2)**2)


      
        !endif ! from if

        !else   ! Guard cells
          ! calculation of cvBb at the guard cells

          ! 28/4/22 : take the cvBb of a guard cell equal the Bb in the face center of the guard cell
          ! Magnetic field on face center is the average of the magnetic field on the two vertices
          ! -Later improve this 
        
         ! iCvBnd = ghostGetBndCell(iCv,g,m)  ! get the boundary cell neighboring the ghost cell
        
          !Temporary equal to zero
          !g%cvBb(iCv,0) = 0.0e0
          !g%cvBb(iCv,1) = 0.0e0
          !g%cvBb(iCv,2) = 0.0e0
          !g%cvBb(iCv,3) = 0.0e0


        endif  !from iCv ne GRID_GUARD


      enddo

      !  ..loop over faces to immediatelly compute magnetic field on the face centers
      !  .. trying analog computations as in b2agbf

      do iFc = 1,m%nFc
        
        vx2 = m%fcVx(iFc,:)
        if (isymm.eq.0) then
           hzf = width
        else if (isymm.eq.1 .or. isymm.eq.2) then
           hzf = pi * ( g%vxX(vx2(0)) + g%vxX(vx2(1)) )
        else if (isymm.eq.3 .or. isymm.eq.4) then
           hzf = pi * ( g%vxY(vx2(0)) + g%vxY(vx2(1)) )
        endif
        
        g%fcBb(iFc,0) = 0.5_R8 *   
     &   ( sqrt( g%vxBb(vx2(0),0)**2  + g%vxBb(vx2(0),1)**2) + 
     &     sqrt( g%vxBb(vx2(1),0)**2  + g%vxBb(vx2(1),1)**2) ) / hzf

        g%fcBb(iFc,1) = 0.0_R8
        g%fcBb(iFc,2) = 0.5_R8 *(g%vxFfbz(vx2(0))  + 
     &                    g%vxFfbz(vx2(1)) ) /hzf

      enddo

      g%fcBb(:,3) = sqrt( g%fcBb(:,0)**2 + 
     &   g%fcBb(:,1)**2 + g%fcBb(:,2)**2)

      end subroutine b2agbb



