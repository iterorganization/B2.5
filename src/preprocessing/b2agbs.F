      subroutine b2agbs (g,m)
      use b2mod_types
      use b2mod_b2cmfs
      use b2mod_constants
      use b2us_map
      use b2us_geo
      implicit none
*   ..input arguments (unchanged on exit)
      type(mapping) :: m
*   ..input/ouput arguments (defined on entry, changed upon exit)
      type(geometry) :: g
*   ..local variables
      integer :: iCv, sv, nv, verts(0:100), ref_ind_cv(m%nFt), iFt,
     .  ref_ind_fc(m%nFt), sc, nc, cells(0:m%nCv), faces(0;m%nFc), iFc,
     .  i, omp_cells(0:m%nCv), imp_cells(0:m%nCv)
      data width /1.0_R8/
      save width
*     ------------------------------------------------------------------
*     B2AGBS computes a correction to the pbsx array to ensure the
*     constancy of magnetic flux along flux tubes for unstructured 
*     grids.
*     ------------------------------------------------------------------
*.computation

* ..preliminaries
*   ..subprogram start-up calls
      call subini ('b2agbs')
      call ipgetr ('b2agmt_1d_width', width)

*   ..preparatory quantities
      g%cvHz = g%cvVol/g%cvSz
      if (isymm.eq.0) then
        g%fcHz(:) = width
      else if (isymm.eq.1 .or. isymm.eq.2) then
        g%fcHz(:) = pi * ( g%vxX(m%fcVx(:,1)) + g%vxX(m%fcVx(:,2)) )
      else if (isymm.eq.3 .or. isymm.eq.4) then
        g%fcHz(:) = pi * ( g%vxY(m%fcVx(:,1)) + g%vxY(m%fcVx(:,2)) )
      endif 

*   ..toroidal magnetic field
      if (isymm.eq.1 .or. isymm.eq.3) then
        if (maxval(abs(g%vxFfbz(:))).gt.0.0_R8) then
          ! Toroidal field in cell center
          do iCv = 1,m%nCv
            sv = m%cvVxP(iCv,1)
            nv = m%cvVxP(iCv,2)
            verts(0:nv-1) = m%cvVx(sv:sv+nv-1)
            g%cvBb(iCv,2) = sum(g%vxFfbz(verts(0:nv-1)))/ 
     &         ( nv * g%cvHz(iCv))
          enddo
          ! Toroidal field in face center
          g%fcBb(:,2) = (g%vxFfbz(m%fcVx(:,1))+g%vxFfbz(m%fcVx(:,2))/
     &        (2.0_R8 * g%fcHz(:))
        endif
      endif
      g%fcBb(:,1) = 0.0_R8

*   ..identify flux tubes - already done in GOAT

*   ..select reference location for each flux tube
*   ..we use the outer midplane whenever possible, then the inner midplane,
*   ..otherwise the cell to the right of the left end of the flux tube.
      ! Get omp cells
      omp_cells
      ! Get imp cells
      imp_cells
 
      do iFt = 1,m%nFt
        ! Get cells
        sc = m%ftCvP(iFt,1)
        nc = m%ftCvP(iFt,2)
        cells(0:nc-1) = m%ftCv(sc:sc+nc-1)
   
        ! If there is a common cell between omp cells and tube cells
           ! Take that cell and one of its poloidal faces
        ! Elseif there is common cell between imp cells and tube cells
           ! Take that cell and one of its poloidal faces
        ! Else  do something else
        
         ! dummy for the moment
         ref_ind_fc(iFt) = 1
         ref_ind_cv(iFt) = 1
      enddo

*   ..correct fcPbs (pbsx), cvHx (hx), and cvHy (hy) to enforce:
*   ..B_{pol}*hy1*hz constant (on western cell faces and at cell centers)
*   ..vol = hy1*hx*hz - (already determined by cvHz = cvVol/cvSz?)
      do iFt = 1,m%nFt
        ! Get cells
        sc = m%ftCvP(iFt,1)
        nc = m%ftCvP(iFt,2)
        cells(0:nc-1) = m%ftCv(sc:sc+nc-1)

        ! Get faces
        sf = m%ftFcP(iFt,1)
        nf = m%ftFcP(iFt,2)
        faces(0:nf-1) = m%ftFc(sf:sf+nf-1)

        ! Get reference
        pol_flux = g%fcBb(ref_ind_fc(iFt),0) 
     &           * g%fcS(ref_ind_fc(iFt)) 
     &           * g%fcQalf(ref_ind_fc(iFt),0)

        ! Apply on the tube - faces
        do i = 0,nf-1
           iFc = faces(i)

           ! Adjust fcBb(:,0)
           g%fcBb(iFc,0) = pol_flux / ( g%fcS(iFc) * g%fcQalf(iFc,0))

           ! Recompute fcBb(:,3)
           g%fcBb(iFc,3) = sqrt(g%fcBb(iFc,0)**2 + g%fcBb(iFc,1)**2 
     &                     + g%fcBb(iFc,2)**2)
            
	   ! Recompute fcPbs(:)
           g%fcPbs(iFc) = pol_flux / g%fcBb(iFc,3)
        enddo

       ! HOW TO APPLY FOR CELLS BECAUSE vol = hx*hy*hz already enforced
      enddo

*   ..copying over modified data
       
      g%cvBb(:,3) = sqrt(g%cvBb(:,0)**2 + g%cvBb(:,1)**2 
     &        + g%cvBb(:,2)**2)

*   ..return
      call subend ()
      return
      end subroutine b2agbs
