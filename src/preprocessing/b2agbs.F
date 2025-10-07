      subroutine b2agbs (g,m)
      use b2mod_types
      use b2mod_b2cmfs
      use b2mod_constants
      use b2us_map
      use b2us_geo
      use b2us_cellhelper
      use b2mod_subsys
      implicit none
*   ..input arguments (unchanged on exit)
      type(mapping) :: m
*   ..input/output arguments (defined on entry, changed upon exit)
      type(geometry) :: g
*   ..local variables
      integer :: iCv, sv, nv, sf, nf, verts(0:100), ref_ind_cv(m%nFt),
     .  iFt, ref_ind_fc(m%nFt), sc, nc, cells(0:m%nCv), typef,
     .  faces(0:m%nFc), iFc, i, cv, v1, v2,
     .  scf, ncf, fcs(0:100), c_al, j, ref_fc(m%nFt), ref_cv(m%nFt),
     .  counter, pfaces(0:100), ref_face, ind,
     .  face1, face2, nf2, k, l, ncvOMP, ncvIMP, cvOMP(1:m%nCv),
     .  cvIMP(1:m%nCv)
!      real (kind=R8) :: pol_flux
      real (kind=R8) :: width, vpsi(0:100), psi_max(0:m%nCv),
     .  psif_max(0:m%nFc), psi_ref, pol_flux_ref, cutoff_width ,
     .  pol_flux_ref1, pol_flux_ref2, pol_flux1, pol_flux2, psif(m%nFc),
     .  rzomp(2,2), rzimp(2,2)
      logical :: fc_done(0:m%nFc), cv_done(0:m%nCv),
     .  cv_aligned(0:m%nCv), ifail
      intrinsic int, maxval, maxloc
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
      end if

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

*   ..First determine cells from OMP and IMP
      rzomp(1,:) = m%OMPr
      rzomp(2,:) = m%OMPz
      rzimp(1,:) = m%IMPr
      rzimp(2,:) = m%IMPz
      ncvOMP = 0
      cvOMP = 0
      call cv_intersections (g, m, rzomp,
     &  cvOMP, ncvOMP, m%nCv, 1, ifail)
      call xertst(.not.ifail,'Getting OMP cells failed in b2agbs')
      ncvIMP = 0
      cvIMP = 0
      call cv_intersections (g, m, rzimp,
     &  cvIMP, ncvIMP, m%nCv, 1, ifail)
      call xertst(.not.ifail,'Getting OMP cells failed in b2agbs')
      ref_ind_cv = 0
      do iFt = 1,m%nFt
        ! Get cells
        sc = m%ftCvP(iFt,1)
        nc = m%ftCvP(iFt,2)
        cells = 0
        cells(0:nc-1) = m%ftCv(sc:sc+nc-1)

        ! Get faces
        sf = m%ftFcP(iFt,1)
        nf = m%ftFcP(iFt,2)
        faces = 0
        faces(0:nf-1) = m%ftFc(sf:sf+nf-1)

        ! If there is a common cell between omp cells and tube cells
        do i = 0,nc-1
           if (isMember(cells(i),cvOMP(1:ncvOMP),ncvOMP)) then
              ref_ind_cv(iFt) = i
              exit
           endif
        enddo
        ! Elseif there is common cell between imp cells and tube cells
        if (ref_ind_cv(iFt).eq.0) then
          do i = 0,nc-1
           if (isMember(cells(i),cvIMP(1:ncvIMP),ncvIMP)) then
              ref_ind_cv(iFt) = i
              exit
           endif
          enddo
        endif
        ! Else  take the middle of the tube
        if (ref_ind_cv(iFt).eq.0) then
             ref_ind_cv(iFt) = int(0.5_R8*nc-1)
        endif

        ref_cv(iFt) = cells(ref_ind_cv(iFt))

        ! Get the poloidal face of that cell
        cv = ref_cv(iFt);

        ! Get a poloidal face
        ind = -1
        do i = 1,m%cvFcP(cv,2)
          iFc = m%cvFc(m%cvFcP(cv,1)+i-1)
          if (m%fcAligned(iFc).eq.0) then
             ! Get index of reference
             !ind = pack( wrkf , faces == iFc )
             do j = 0,nf-1
                 if ( faces(j) .eq. iFc ) then
                    ind = j
                    ref_ind_fc(iFt) = ind
                    ref_fc(iFt) = iFc
                    exit
                 endif
             enddo
             if (ind .ne. -1 ) exit
          endif
        enddo

        ! Check if face is assigned
        if ( ind .eq. -1 ) then
          call xerrab('b2agbs: no reference face assigned')
        endif

      enddo

*   ..correct fcPbs (pbsx), cvHx (hx), and cvHy (hy) to enforce:
*   ..B_{pol}*hy1*hz constant (on western cell faces and at cell centers)
*   ..vol = hy1*hx*hz - (already determined by cvHz = cvVol/cvSz?)

      ! Compute psi width of faces
      do iFc = 1,m%nFc
        v1 = m%fcVx(iFc,1)
        v2 = m%fcVx(iFc,2)

        psif(iFc) = abs(g%vxFpsi(v1) - g%vxFpsi(v2))
      enddo

      do iFt = 1,m%nFt
        ! Get cells
        cells = 0
        sc = m%ftCvP(iFt,1)
        nc = m%ftCvP(iFt,2)
        cells = 0
        cells(0:nc-1) = m%ftCv(sc:sc+nc-1)

        ! Get poloidal faces in the tube
        sf = m%ftFcP(iFt,1)
        nf = m%ftFcP(iFt,2)
        faces  = 0
        faces(0:nf-1) = m%ftFc(sf:sf+nf-1)

        ! Determine type of tube
        if ( nc .eq. nf ) then
          typef = 1 ! closed
        else
          typef = 2 ! open
        endif

        psi_max = 0.0_R8
        cv_aligned = .false.
        cv_done = .false.
        fc_done = .false.
        do i = 0,nc-1
          ! Get cell
          iCv = cells(i)

          ! Psi max
          verts = 0
          vpsi = 0
          sv = m%cvVxP(iCv,1)
          nv = m%cvVxP(iCv,2)
          verts(0:nv-1) = m%cvVx(sv:sv+nv-1)
          vpsi(0:nv-1) = g%vxFpsi(verts(0:nv-1))
          psi_max(i) = abs(maxval(vpsi(0:nv-1))
     &       - minval(vpsi(0:nv-1)))

          ! Determine aligned
          fcs = 0
          scf = m%cvFcP(iCv,1)
          ncf = m%cvFcP(iCv,2)
          fcs(0:ncf-1) = m%cvFc(scf:scf+ncf-1)

          c_al = 0
          do j = 0,ncf-1
            if (g%fcQalf(fcs(j),0) .eq. 0.0_R8) then
              c_al = c_al + 1
            endif
          enddo

          if  ((( ncf .eq. 4 ) .and. ( c_al .eq. 2 )) .or.
     &     ( ( ncf .eq. 3 ) .and. ( c_al .eq. 1)  )) then
            cv_aligned(i) = .true.
          elseif ( ncf .gt. 4 ) then
            call xerrab('b2agbs: Pentagons or higher not supported')
          endif

        enddo

        ! Get psi width of face in the tube
        psif_max = 0
        do j = 0,nf-1
          psif_max(j) = psif(faces(j))
        enddo

        ! Do the adjustement for the cells with proper width and alignment
        ! Loop over cells with aligned properties
        ! Get refernce of psi width and pol_flux
        psi_ref = psi_max(ref_ind_cv(iFt))
        pol_flux_ref = g%fcBb(ref_fc(iFt),0)*g%fcS(ref_fc(iFt))
     &    *g%fcQalf(ref_fc(iFt),0)
        cutoff_width = 0.05_R8

        ! Check if reference cell is aligned
        if ( cv_aligned(ref_ind_cv(iFt)) ) then

          ! Loop over cells
          do k = 0,nc-1
            if ( cv_aligned(k) ) then
               if ((( m%cvLbl(cells(k)) .eq. 1)
     &          .and. (m%cvFcP(cells(k),2) .eq. 4))
     &          .or. ( abs( (psi_ref - psi_max(k))/psi_ref)
     &          .lt. cutoff_width )) then
                cv_done(k) = .true.

                 ! Adjust pitch of faces(k) and next face
                 call AdjustPitch(g,faces(k),pol_flux_ref)
                 fc_done(k) = .true.
                 if ( typef .eq. 1 ) then
                   ! Closed tube
                   if ( k .ne. nc-1 ) then

                     call AdjustPitch(g,faces(k+1),pol_flux_ref)
                     fc_done(k+1) = .true.

                   else

                     call AdjustPitch(g,faces(0),pol_flux_ref)
                     fc_done(0) = .true.

                   endif
                 elseif ( typef .eq. 2 ) then
                   ! Open tube
                   call AdjustPitch(g,faces(k+1),pol_flux_ref)
                   fc_done(k+1) = .true.

                 else

                   call xerrab('b2agbs: type not defined')

                 endif
               endif
            endif

          enddo

        else

          ! Get reference cell
          cv = ref_cv(iFt)
          sf = m%cvFcP(cv,1)
          nf2 = m%cvFcP(cv,2)
          fcs = 0
          fcs(0:nf2-1) = m%cvFc(sf:sf+nf2-1)

          ! Not aligned to a quad
          if ( nf2 .ne. 4 ) then
            call xerrab('b2agbs: Something wrong')
          endif

          ! Excluded the aligned face
          counter = 0
          pfaces = 0
          do i = 0,nf2-1
            if (g%fcQalf(fcs(i),0) .ne. 0.0_R8) then
              pfaces(counter) = fcs(i)
              counter = counter + 1
            endif
          enddo

          if ( counter .ne. 3 )  then
            call xerrab('b2agbs: not supported')
          else
            ! Get face with max psi
            ind = maxloc(psif(pfaces(0:counter-1)),1) - 1
            ref_face = pfaces(ind)

            pol_flux_ref = abs(g%fcBb(ref_face,0)*g%fcS(ref_face)
     &       *g%fcQalf(ref_face,0))

            ! Two other faces need to the sum of the reference
            do i = 0,2
              if ( pfaces(i) .ne. ref_face ) then
                face1 = pfaces(i)
              endif
            enddo

            do i = 0,2
              if (( pfaces(i).ne. ref_face)
     &           .and. (pfaces(i) .ne. face1 )) then
                face2 = pfaces(i)
              endif
            enddo
          endif

          pol_flux1 = abs(g%fcBb(face1,0)*g%fcS(face1)
     &       *g%fcQalf(face1,0))
          pol_flux2 = abs(g%fcBb(face2,0)*g%fcS(face2)
     &       *g%fcQalf(face2,0))

          pol_flux_ref1 = pol_flux_ref * pol_flux1
     &      /(pol_flux1 + pol_flux2)
          pol_flux_ref2 = pol_flux_ref * pol_flux2
     &      /(pol_flux1 + pol_flux2)

          call AdjustPitch(g,face1,pol_flux_ref1)
          call AdjustPitch(g,face2,pol_flux_ref2)

          ! Update fc_done and cv_done
          cv_done(ref_ind_cv(iFt)) = .true.
          fc_done(ref_ind_cv(iFt)) = .true.
          fc_done(ref_ind_cv(iFt)+1) = .true.

        endif

        ! Complete cells which are not of full width
        ! Forward propagation
        do j = 1,nc-1 ! Starting from second cell
          if (( .not. cv_done(j) )
     &     .and. ( cv_done(j-1) )) then

             if ( fc_done(j) ) then

                call propagate(g,m,cells(0:nc-1),j,faces(0:nf-1),
     &       j+1,cv_done(0:nc-1),fc_done(0:nf-1),psif_max(0:nf-1),
     &        nc,nf)
             else

               call xerrab
     &    ('b2agbs: forward prop: fc_done of cv_done cells not done')

             endif

          endif

        enddo

        ! Backward propagation
        do  j = 1,nc-1
           l = nc - 1 - j
           if  (( .not. cv_done(l) )
     &        .and. ( cv_done(l+1) ) ) then

             if ( fc_done(l+1) ) then

                call propagate(g,m,cells(0:nc-1),l,faces(0:nf-1),
     &     l,cv_done(0:nc-1),fc_done(0:nf-1),psif_max(0:nf-1),nc,nf)

             else

               call xerrab
     &    ('b2agbs: backward prop: fc_done of cv_done cells not done')

             endif
           endif
        enddo

!        ! Get reference
!        pol_flux = g%fcBb(ref_ind_fc(iFt),0)
!     &           * g%fcS(ref_ind_fc(iFt))
!     &           * g%fcQalf(ref_ind_fc(iFt),0)
!
!        ! Apply on the tube - faces
!        do i = 0,nf-1
!           iFc = faces(i)
!
!          ! Adjust fcBb(:,0)
!          g%fcBb(iFc,0) = sign(1.0,g%fcBb(iFc,0))
!    &        * abs(pol_flux / ( g%fcS(iFc) * g%fcQalf(iFc,0)))

!           ! Recompute fcBb(:,3)
!           g%fcBb(iFc,3) = sqrt(g%fcBb(iFc,0)**2 + g%fcBb(iFc,1)**2
!     &                     + g%fcBb(iFc,2)**2)
!
            ! Recompute fcPbs(:)
!           g%fcPbs(iFc) = sign(1.0,g%fcPbs(iFc))
!     &       * abs(pol_flux / g%fcBb(iFc,3))
!        enddo

      enddo ! Loop over tubes

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
      implicit none

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


      subroutine propagate(g,m,cells,indc,faces,
     &   indf,cv_done,fc_done,psif_max,ncpt,nfpt)
      use b2mod_types
      use b2us_geo
      use b2us_map
      implicit none

*   ..input and output arguments
      type(geometry) :: g
      type(mapping) :: m
      integer :: ncpt, nfpt, indc, indf, cells(0:ncpt-1),
     &     faces(0:nfpt-1)
      logical :: cv_done(0:ncpt-1), fc_done(0:nfpt-1)
      real(kind=R8) :: psif_max(0:nfpt-1)

      ! Local variables
      integer :: k, l, icv, nf, s, fcs_s(0:3), fcs_incl, i, fc,
     &         fcs_incl_a(0:100), n_incl
      real(kind=R8) :: pol_flux_ref, pol_flux_ref_k2, pol_flux_ref_k1,
     &                 pol_flux_ref_incl, pol_flux_k1, pol_flux_k2,
     &                 pol_flux_incl
      intrinsic abs

      ! Get right variable name
      k = indc
      l = indf
      icv = cells(k)
      nf = m%cvFcP(icv,2)

      if ( nf .eq. 3 ) then

        ! Just propagate forward or backward
        if (k .eq. l ) then

          ! Reference is faces(k+1)
          pol_flux_ref = abs(g%fcBb(faces(k+1),0)*g%fcS(faces(k+1))
     &           *g%fcQalf(faces(k+1),0))
          call AdjustPitch(g,faces(k),pol_flux_ref)

        elseif (k .lt. l ) then

          ! Reference is faces(k)
          pol_flux_ref = abs(g%fcBb(faces(k),0)*g%fcS(faces(k))
     &           *g%fcQalf(faces(k),0))
          call AdjustPitch(g,faces(k+1),pol_flux_ref)

        else

          call xerrab('b2agbs: k and l indices not conform')

        endif

        ! Indicate done faces and cells
        fc_done(l) = .true.
        cv_done(k) = .true.

      elseif (nf .eq. 4 ) then

        ! Get inclined face
        s = m%cvFcP(icv,1)
        fcs_s = m%cvFc(s:s+nf-1)
        fcs_incl = 0
        fcs_incl_a = 0
        n_incl = 0
        do i = 0, nf-1
          fc = fcs_s(i)
          if ((fc .ne. faces(k)) .and. (fc .ne. faces(k+1))
     &        .and. (g%fcQalf(fc,0) .ne. 0.0_R8)) then
            fcs_incl_a(n_incl) = fc
            n_incl = n_incl + 1
            !if (fcs_incl .ne. 0 ) exit
          endif
        enddo

        if ( n_incl .eq. 0 )  then
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

        elseif ( n_incl .eq. 1 ) then
          ! Get inclined face
          fcs_incl = fcs_incl_a(0)

          ! Complex propagation
          if (k .eq. l) then
            ! Reference is faces(k+1)
            if ( fc_done(l) ) then
              ! Adjust the incline face
              ! Compression or expansion
                pol_flux_k1 = abs(g%fcBb(faces(k),0)
     &           * g%fcS(faces(k)) * g%fcQalf(faces(k),0))
                pol_flux_k2 = abs(g%fcBb(faces(k+1),0)
     &           *g%fcS(faces(k+1)) * g%fcQalf(faces(k+1),0))

              if ( psif_max(k+1) .gt. psif_max(k)) then

                ! Compression
                pol_flux_ref = abs(pol_flux_k2 -
     &             pol_flux_k1)

              else
                ! Expansion
                 pol_flux_ref = abs(pol_flux_k1 - pol_flux_k2)

              endif

              ! Adjust the pitch
              call AdjustPitch(g,fcs_incl,pol_flux_ref)

            else

              ! Adjust fc_done(l)
              ! compression or expansion
              pol_flux_incl = abs(g%fcBb(fcs_incl,0)
     &          *g%fcS(fcs_incl)*g%fcQalf(fcs_incl,0))
              pol_flux_k2 = abs(g%fcBb(faces(k+1),0)
     &          *g%fcS(faces(k+1))*g%fcQalf(faces(k+1),0))

              if (psif_max(k+1) .gt. psif_max(k)) then
                ! Compression
                pol_flux_k1 = abs(g%fcBb(faces(k),0)
     &           * g%fcS(faces(k))*g%fcQalf(faces(k),0))

                pol_flux_ref_incl = pol_flux_k2
     &           *pol_flux_incl
     &             /(pol_flux_k1 + pol_flux_incl)
                pol_flux_ref_k1 = pol_flux_k2
     &           *pol_flux_k1
     &            /(pol_flux_k1 + pol_flux_incl)

                call AdjustPitch(g,faces(k),pol_flux_ref_k1)
                call AdjustPitch(g,fcs_incl,pol_flux_ref_incl)

              else
                ! Expansion
                 pol_flux_ref_k1 = abs(pol_flux_k2
     &             + pol_flux_incl)
                call AdjustPitch(g,faces(k),pol_flux_ref_k1)

              endif

            endif

          elseif ( l .eq. k + 1) then ! forward
            ! Reference is faces(k)
            if (fc_done(k+1)) then
              ! Adjust inclined face
              pol_flux_k1 = abs(g%fcBb(faces(k),0)
     &         *g%fcS(faces(k))*g%fcQalf(faces(k),0))
              pol_flux_k2 = abs(g%fcBb(faces(k+1),0)
     &         *g%fcS(faces(k+1))*g%fcQalf(faces(k+1),0))

              if ( psif_max(k) .gt. psif_max(k+1) ) then
                ! compression - situation 3
                pol_flux_ref_incl = abs(pol_flux_k1
     &            - pol_flux_k2)
              else
                ! expansion - situation 2
                pol_flux_ref_incl = abs(pol_flux_k2
     &            - pol_flux_k1)
              endif

              ! Adjust pitch
              call AdjustPitch(g,fcs_incl,pol_flux_ref_incl)

            else
              !  Adjust fc_done(k+1) and inclined face
              pol_flux_incl = abs(g%fcBb(fcs_incl,0)
     &          *g%fcS(fcs_incl)*g%fcQalf(fcs_incl,0))
              pol_flux_k1 = abs(g%fcBb(faces(k),0)
     &         *g%fcS(faces(k))*g%fcQalf(faces(k),0))

              if (psif_max(k) .lt. psif_max(k+1)) then
                ! expansion - situation 4
                pol_flux_ref_k2 = abs(pol_flux_incl
     &              + pol_flux_k1)

              else
                ! compression - situation 1
                pol_flux_k2 = abs(g%fcBb(faces(k+1),0)
     &          *g%fcS(faces(k+1))*g%fcQalf(faces(k+1),0))

                pol_flux_ref_incl = pol_flux_k1
     &           *pol_flux_incl
     &           /(pol_flux_k2 + pol_flux_incl)
                pol_flux_ref_k2 = pol_flux_k1 * pol_flux_k2
     &           /(pol_flux_k2 + pol_flux_incl)

                call AdjustPitch(g,fcs_incl,pol_flux_ref_incl)
              endif

              call AdjustPitch(g,faces(k+1),pol_flux_ref_k2)

            endif

          else
            call xerrab
     &   ('b2agbs: propagate: something wrong in forward propagation')
          endif

          ! Update fc_done and cv_done
          fc_done(l) = .true.
          cv_done(k) = .true.
        elseif ( n_incl .gt. 1 ) then
           write(*,*)'Warning: b2agbs: cells without aligned faces '//
     &     'are not yet supported'
        endif

      else
        call xerrab
     &   ('b2agbs: propagate: pentagons or higher not supported')
      endif

      return
      end subroutine propagate
