      subroutine b2agbs (g,m)
      use b2mod_types
      use b2mod_b2cmfs
      use b2mod_constants
      use b2us_map
      use b2us_geo
      use b2us_cellhelper
      implicit none
*   ..input arguments (unchanged on exit)
      type(mapping) :: m
*   ..input/ouput arguments (defined on entry, changed upon exit)
      type(geometry) :: g
*   ..local variables
      integer :: iCv, sv, nv, sf, nf, verts(0:100), ref_ind_cv(m%nFt), 
     .  iFt, ref_ind_fc(m%nFt), sc, nc, cells(0:m%nCv), 
     .  faces(0:m%nFc), iFc, i, cv
      real (kind=R8) :: pol_flux, width
      intrinsic int
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
          g%fcBb(:,2) = (g%vxFfbz(m%fcVx(:,1))+g%vxFfbz(m%fcVx(:,2)))/
     &        (2.0_R8 * g%fcHz(:))
        endif
      endif
      g%fcBb(:,1) = 0.0_R8

*   ..identify flux tubes - already done in GOAT

*   ..select reference location for each flux tube
*   ..we use the outer midplane whenever possible, then the inner midplane,
*   ..otherwise the cell to the right of the left end of the flux tube.
    
      ref_ind_cv = 0
      do iFt = 1,m%nFt
        ! Get cells
        sc = m%ftCvP(iFt,1)
        nc = m%ftCvP(iFt,2)
        cells(0:nc-1) = m%ftCv(sc:sc+nc-1)
        
        ! Get faces
        sf = m%ftFcP(iFt,1)
        nf = m%ftFcP(iFt,2)
        faces(0:nf-1) = m%ftFc(sf:sf+nf-1)
   
        ! If there is a common cell between omp cells and tube cells
        do i = 0,nc-1
           if (isMember(cells(i),m%cvOMP,m%ncvOMP)) then
              ref_ind_cv(iFt) = cells(i)
              exit
           endif
        enddo
        ! Elseif there is common cell between imp cells and tube cells
        if (ref_ind_cv(iFt).eq.0) then
          do i = 0,nc-1
           if (isMember(cells(i),m%cvIMP,m%ncvIMP)) then
              ref_ind_cv(iFt) = cells(i)
              exit
           endif   
          enddo 
        endif
        ! Else  take the middle of the tube
        if (ref_ind_cv(iFt).eq.0) then
             ref_ind_cv(iFt) = cells(int(0.5_R8*nc-1))
        endif

        ! Get the poloidal face of that cell
        cv = ref_ind_cv(iFt);
        
        ! Get a poloidal faces
        do i = 1,m%cvFcP(cv,2)
          iFc = m%cvFc(m%cvFcP(cv,1)+i-1)
          if (m%fcAligned(iFc).eq.0) then
             ref_ind_fc(iFt) = iFc
             exit
          endif
        enddo
      
      enddo

*   ..correct fcPbs (pbsx), cvHx (hx), and cvHy (hy) to enforce:
*   ..B_{pol}*hy1*hz constant (on western cell faces and at cell centers)
*   ..vol = hy1*hx*hz - (already determined by cvHz = cvVol/cvSz?)
      do iFt = 1,m%nFt
        ! Get cells
        sc = m%ftCvP(iFt,1)
        nc = m%ftCvP(iFt,2)
        cells(0:nc-1) = m%ftCv(sc:sc+nc-1)

        ! Get poloidal faces in the tube
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
           g%fcBb(iFc,0) = sign(1.0,g%fcBb(iFc,0))
     &        * abs(pol_flux / ( g%fcS(iFc) * g%fcQalf(iFc,0)))

           ! Recompute fcBb(:,3)
           g%fcBb(iFc,3) = sqrt(g%fcBb(iFc,0)**2 + g%fcBb(iFc,1)**2 
     &                     + g%fcBb(iFc,2)**2)
            
	   ! Recompute fcPbs(:)
           g%fcPbs(iFc) = sign(1.0,g%fcPbs(iFc))
     &       * abs(pol_flux / g%fcBb(iFc,3))
        enddo

      enddo

*   ..copying over modified data
       
      g%cvBb(:,3) = sqrt(g%cvBb(:,0)**2 + g%cvBb(:,1)**2 
     &        + g%cvBb(:,2)**2)

*   ..return
      call subend ()
      return
      end subroutine b2agbs

*----------------------------------------------------------------------
      
      subroutine AdjustPitch(g,iFc,pol_flux_ref)
      use b2mod_types
      use b2us_geo

*   ..input and output argument      
      type(geometry) :: g
*   ..input
      integer :: iFc
      real(kind=R8) :: pol_flux_ref
      intrinsic sign, abs, sqrt

      ! Adjust fcBb
      g%fcBb(iFc,0) = sign(1.0_R8,g%fcBb(iFc,0)) 
     &  * abs(pol_flux_ref/(g%fcS(iFc)*g%fcQalf(iFc,0)))

      ! Recompute fcBb(:,3)
      g%fcBb(iFc,3) = sqrt( g%fcBb(iFc,0)**2 + g%fcBb(iFc,1)**2
     &  + g%fcBb(iFc,2)**2 ) 

      ! Recompute pitch
      g%fcPbs(iFc) = sign(1.0_R8,g%fcPbs(iFc)) 
     &  * abs(pol_flux_ref/ g%fcBb(iFc,3))

      end subroutine AdjustPitch

*----------------------------------------------------------------------
    
      subroutine CompPolflux(g,iFc,pol_flux)
      use b2mod_types
      
      end subroutine 

*----------------------------------------------------------------------


      subroutine propagate
      &   (g,cells,indc,faces,indf,cv_done,fc_done,psif_max,nct,nft)
      use b2mod_types
      use b2us_geo
      
*   ..input and output arguments
      type(geometry) :: g
      integer :: ncpt, nfpt, indc, indf
      logical :: cv_done(0:ncpt), fc_done(0:nfpt)
      real(kind=R8) :: psif_max(0:nfpt)


      ! Local variables
      integer :: k, l, icv, nf, s, fcs_s(0:3), fcs_incl, i, fc
      real(kind=R8) :: pol_flux_ref, pol_flux_ref_k2, pol_flux_ref_k1,
     &                 pol_flux_ref_incl


      ! Get right variable name
      k = indc
      l = indf
      icv = cells(k)
      nf = g%cvFcP(icv,2)

      if ( nf .eq. 3 ) then

        ! Just propagate forward or backward
        if (k .eq. l ) then

          ! Reference is faces(k+1)
          pol_flux_ref = g%fcBb(faces(k+1),0)*g%fcS(faces(k+1))
     &           g%fcQalf(faces(k+1),0)
          call AdjustPitch(g,faces(k),pol_flux_ref) 
    
        elseif (k .lt. l ) then

          ! Reference is faces(k)
          pol_flux_ref = g%fcBb(faces(k),0)*g%fcS(faces(k))
     &           g%fcQalf(faces(k),0)          
          call AdjustPitch(g,faces(k),pol_flux_ref) 
         
        else

          call xerrab('b2agbs: something wrong')

        endif
       
        ! Indicate done faces and cells  
        fc_done(l) = .true.
        cv_done(k) = .true.

      elseif (nf .eq. 4 ) then

        ! Get inclined face
        s = g%cvFcP(icv,1)
        fcs_s = g%cvFc(s:s+nf-1)
        fcs_incl = 0
        do i = 0, nf-1
          fc = fcs_s(i)
          if ((fc .ne. faces(k)) .and. (fc .ne. faces(k+1))
     &        .and. (g%fcQalf(fc,0) .ne. 0.0_R8)) then
            fcs_incl = fc
            if (fcs_incl .ne. 0 ) exit
          endif 
        enddo

        if fcs_incl .eq. 0  then
          ! Just propagate
          if (k .eq. l ) then
            pol_flux_ref_k2 = abs(g%fcBb(faces(k+1),0)
     &        *g%fcS(faces(k+1))*g%fcQalf(faces(k+1),0))
            call AdjustPitch(g,faces(k),pol_flux_ref_k2)
          elseif (l .eq. k + 1) then
            pol_flux_ref_k1 = abs(g%fcBb(faces(k),0)
     &        *g%fcS(faces(k))*g%fcQalf(faces(k),0))
            call AdjustPitch(g,faces(k+1),pol_flux_ref_k1)   
          endif

	  ! Indicate done faces and cells 
          fc_done(l) = .true.
          cv_done(k) = .true.          
        
        else 
          ! Complex propagation
          if (k .eq. l) then
            ! Reference is faces(k+1) 
            if ( fc_done(l) .eq. .true. ) then
              ! Adjust the incline face
              ! Compression or expansion
                pol_flux_ref_k1 = abs(g%fcBb(faces(k),0)
     &           * g%fcS(faces(k)) * g%fcQalf(faces(k),0))
                pol_flux_ref_k2 = abs(g%fcBb(faces(k+1),0) 
     &           *g%fcS(faces(k+1)) * g%fcQalf(faces(k+1),0))

              if ( psif_max(k+1) .gt. psif_max(k)) then

                ! Compression
                pol_flux_ref = abs(pol_flux_ref_k2 - 
     &             pol_flux_ref_k1)
                
              
              else
                ! Expansion
                 pol_flux_ref = abs(pol_flux_ref_k1 - pol_flux_ref_k2)
                 
              endif

              ! Adjust the pitch
              call AdjustPitch(g,fcs_incl,pol_flux_ref)

            else

	      ! Adjust fc_done(l)
              ! compression or expansion
              pol_flux_ref_incl = abs(g%fcBb(fcs_incl,0)
     &          *g%fcS(fcs_incl)*g%fcQalf(fcs_incl,0))
              pol_flux_ref_k2 = abs(g%fcBb(faces(k+1),0)
     &          *g%fcS(faces(k+1))*g%fcQalf(faces(k+1),0))

              if (psif_max(k+1) .gt. psif_max(k)) then
                ! Compression
                pol_flux_ref_k1 = abs(g%fcBb(faces(k),0)
     &           * g%fcS(faces(k))*g%fcQalf(faces(k),0))

                pol_flux_ref_incl = pol_flux_ref_k2 
     &           *pol_flux_ref_incl
     &             /(pol_flux_ref_k1 + pol_flux_ref_incl)
                pol_flux_ref_k1 = pol_flux_ref_k2
     &           *pol_flux_ref_k1
     &            /(pol_flux_ref_k1 + pol_flux_ref_incl)

                call AdjustPitch(g,faces(k),pol_flux_ref_k1)

              else
                ! Expansion
                 pol_flux_ref_incl = abs(pol_flux_ref_k2 
     &             pol_flux_ref_incl)

              endif 
              
              call AdjustPitch(g,fcs_incl,pol_flux_ref_incl)
              
            endif

          elseif ( l .eq. k + 1) then ! forward
            ! Reference is faces(k)
            if (fc_done(k+1) .eq. .true. ) then
              ! Adjust inclined face
              pol_flux_ref_k1 = abs(g%fcBb(faces(k),0)
     &         *g%fcS(faces(k))*g%fcQalf(faces(k),0))
              pol_flux_ref_k2 = abs(g%fcBb(faces(k+1),0)
     &         *g%fcS(faces(k+1))*g%fcQalf(faces(k+1),0))

              if ( psif_max(k) .gt. psif_max(k+1) ) then
                ! compression - situation 3
                pol_flux_ref_incl = abs(pol_flux_k1 - pol_flux_k2)
              else 
                ! expansion - situation 2
                pol_flux_ref_incl = abs(pol_flux_ref_k2 
     &            - pol_flux_ref_k1)
              endif
              
              ! Adjust pitch
              call AdjustPitch(g,fcs_incl,pol_flux_ref_incl)
              
            else
              !  Adjust fc_done(k+1) and inclined face
              pol_flux_ref_incl = = abs(g%fcBb(fcs_incl,0)
     &          *g%fcS(fcs_incl)*g%fcQalf(fcs_incl,0))
              pol_flux_ref_k1 = abs(g%fcBb(faces(k),0)
     &         *g%fcS(faces(k))*g%fcQalf(faces(k),0))

              if (psif_max(k) .lt. psif_max(k+1)) then
                ! expansion - situation 4
                pol_flux_ref_k2 = abs(pol_flux_rel_inc 
     &              + pol_flux_ref_k1)

              else
                ! compression - situation 1
                pol_flux_ref_k2 = abs(g%fcBb(faces(k+1),0)
     &          *g%fcS(faces(k+1))*g%fcQalf(faces(k+1),0))   

                pol_flux_ref_incl = pol_flux_ref_k1
     &           *pol_flux_ref_incl
     &           /(pol_flux_ref_k2 + pol_flux_ref_incl)

                call AdjustPitch(g,fcs_incl,pol_flux_ref_incl)             
 	      endif
             
              call AdjustPitch(g,faces(k+1),pol_flux_ref_incl)

            endif

          else
            call xerrab('b2agbs: something wrong in propagation')
          endif

          ! Update fc_done and cv_done
          fc_done(l) = .true.
          cv_done(k) = .true.

	endif 

      else 
        call xerrab('b2agbs: pentagons or higher not supported')
      endif 

      
      end subroutine propagate





