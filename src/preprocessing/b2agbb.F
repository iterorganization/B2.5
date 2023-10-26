*-----------------------------------------------------------------------
*.specification

      subroutine b2agbb (g,m,isymm)
      use b2mod_types
      use b2us_geo
      use b2us_map
      use carre_constants
      use b2mod_subsys
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
      integer iCv, iFc, vx2(0:1), vx0, iCvBnd, nv
      integer, allocatable :: vx(:)
      real (kind=R8) ::
     *  t0, pi, babs, psdx, psdy, width, hzf, babs_test,
     *  centroid(0:1) ! 0 = x-coord, 1 = y-coord
*   ..procedures
      intrinsic sqrt, sign, sum

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
        !if (m%cvLbl(iCv) .ne. GRID_GUARD)  then !guard cells are skipped
        !   ..compute magnetic field at cell center - is approximated here
        nv = m%cvVxP(iCv,2)
        allocate(vx(0:nv-1))
        if (isymm.ne.0) then
           vx = m%cvVx(m%cvVxP(iCv,1):m%cvVxP(iCv,1)+nv-1)
           centroid(0) = (sum(g%vxX(vx))) / nv
           centroid(1) = (sum(g%vxY(vx))) / nv

          if (isymm.eq.1 .or. isymm.eq.2) then
           t0 = 2.0*pi*centroid(0)
          elseif (isymm.eq.3 .or. isymm.eq.4) then
           t0 = 2.0*pi*centroid(1)
          endif

        else
          t0 = 1.0
        endif


        if (t0 == 0.0) cycle
           ! structured code
           !bb(ix,iy,0) = & 
           !       &     (fpsi(ix,iy,2)-fpsi(ix,iy,0)+fpsi(ix,iy,3)-fpsi(ix,iy,1))
           ! not sure always correct but only the sign of cvBb 
           ! is important here -ISSUE
           !g%cvBb(iCv,0) = ( g%vxFpsi(vx4(2)) - g%vxFpsi(vx4(0)) 
           !&      +  g%vxFpsi(vx4(3)) -  g%vxFpsi(vx4(1)) )
          
           !babs = 0.25_R8*(  
           !&      sqrt( g%vxBb(vx4(0),0)**2 + g%vxBb(vx4(0),1)**2 ) + 
           !&      sqrt( g%vxBb(vx4(1),0)**2 + g%vxBb(vx4(1),1)**2 ) +      
           !&      sqrt( g%vxBb(vx4(2),0)**2 + g%vxBb(vx4(2),1)**2 ) + 
           !&      sqrt( g%vxBb(vx4(3),0)**2 + g%vxBb(vx4(3),1)**2 ) )/t0

           vx = m%cvVx(m%cvVxP(iCv,1):m%cvVxP(iCv,1)+nv-1)
           babs = 
     &      sum( sqrt(g%vxBb(vx,0)**2 + g%vxBb(vx,1)**2) ) /
     &       (t0 * nv)
        
           g%cvBb(iCv,0) = babs  !just absolute value
        
          g%cvBb(iCv,1) = 0.0e0
          g%cvBb(iCv,2) = (sum(g%vxFfbz(vx)))/(t0*nv)        
          g%cvBb(iCv,3) = sqrt( g%cvBb(iCv,0)**2 + 
     &     g%cvBb(iCv,1)**2 + g%cvBb(iCv,2)**2)

          deallocate(vx)


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
       ! here direction is lost

       ! other way
       ! fcBx = (g%vxBb(vx2(0),0) + g%vxBb(vx2(1),0)) / 2.0_R8
       ! fcBy = (g%vxBb(vx2(0),1) + g%vxBb(vx2(1),1)) / 2.0_R8
       ! fcBb(iFc,0) = sqrt(fcBx**2 + fcBy**2) / hzf
       ! tan = fcBy / fcBx 
       ! cos = fcBx / sqrt(fcBx**2 + fcBy**2)
       ! sin = fcBy / sqrt(fcBx**2 + fcBy**2)


        g%fcBb(iFc,1) = 0.0_R8
        g%fcBb(iFc,2) = 0.5_R8 *(g%vxFfbz(vx2(0))  + 
     &                    g%vxFfbz(vx2(1)) ) /hzf

      enddo

      g%fcBb(:,3) = sqrt( g%fcBb(:,0)**2 + 
     &   g%fcBb(:,1)**2 + g%fcBb(:,2)**2)

      end subroutine b2agbb



