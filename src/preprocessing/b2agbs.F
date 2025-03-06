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
      integer i, ix, iy, ixl, ixr, iyl, iyr, jxi, jxa, jsep
      integer ntubes(-1:ny), ixref(nncutmax,-1:ny), iseg(-1:nx,-1:ny),
     *  ixleft(nncutmax,-1:ny),
     *  isnum(nncutmax,-1:ny)
      real (kind=R8) :: hz(-1:nx,-1:ny),
     *                  hzl(-1:nx,-1:ny),    ! hz at the left cell face
     *                  pol_flux, width, wrk
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
      ixref = -2
      call get_jsep(nx,ny,jxi,jxa,jsep)
      do iy = -1, ny
        do i = 1, ntubes(iy)
          if (i.eq.iseg(jxa,iy)) then
            ixref(i,iy) = jxa
          else if (i.eq.iseg(jxi,iy)) then
            ixref(i,iy) = jxi
          else
            do ix=1,nx
              if (isRealCell(cflags(ix,iy,CELLFLAG_TYPE)) .and.
     &          (.not.isBoundaryCell(cflags(ix,iy,CELLFLAG_TYPE)))) exit
            enddo
            if (ix.ge.nx) then
              ixref(i,iy) = rightix(ixleft(i,iy),iy)
            else
              ixref(i,iy) = ix
            endif
          end if
        end do
      end do

*   ..First determine cells from OMP and IMP
      rzomp(1,:) = m%OMPr
      rzomp(2,:) = m%OMPz
      rzimp(1,:) = m%IMPr
      rzimp(2,:) = m%IMPz
      ncvOMP = 0
      cvOMP = 0
      call cv_intersections (g, m, rzimp,
     &  cvOMP, ciOMP, fcOMP, ncvOMP, nciOMP,
     &  nfcOMP, m%nCv, xiOMP, yiOMP, 1, ifail)
      call xertst(.not.ifail,'Getting OMP cells failed in b2agbs')
      ncvIMP = 0
      cvIMP = 0
      call cv_intersections (g, m, rzimp,
     &  cvIMP, ciIMP, fcIMP, ncvIMP, nciIMP,
     &  nfcIMP, m%nCv, xiIMP, yiIMP, 1, ifail)

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
*   ..vol = hy1*hx*hz
      do iy = -1, ny
        do ix = -1, nx
          if (region(ix,iy,0).eq.0) cycle
          if (iseg(ix,iy).eq.0) cycle
          if (ixref(iseg(ix,iy),iy).eq.-2) cycle

          if (leftix(ixref(iseg(ix,iy),iy),iy).ne.-2) then
            pol_flux = wbbl(ixref(iseg(ix,iy),iy),iy,0) *
     &                   gs(ixref(iseg(ix,iy),iy),iy,0) *
     &                   qc(ixref(iseg(ix,iy),iy),iy,0)
          else
            pol_flux = bb(ixref(iseg(ix,iy),iy),iy,0) *
     &                 hy(ixref(iseg(ix,iy),iy),iy)
          end if
          wrk = wbbl(ix,iy,0)*gs(ix,iy,0)*qc(ix,iy,0)
          if ((abs(wrk).lt.abs(pol_flux)) .and. 
     &        isBoundaryCell(cflags(ix,iy,CELLFLAG_TYPE)) ) cycle
          if (leftix(ix,iy).ne.-2) then ! do not adjust pbs if leftix=-2
                                        ! since there is no western face
            wbbl(ix,iy,0) = pol_flux / (gs(ix,iy,0) * qc(ix,iy,0))
            wbbl(ix,iy,3) =
     &       sqrt(wbbl(ix,iy,0)**2+wbbl(ix,iy,1)**2+wbbl(ix,iy,2)**2)
            pbs(ix,iy,0) = pol_flux / wbbl(ix,iy,3)
          else
            wbbl(ix,iy,0) = pol_flux /
     &       (gs(rightix(ix,iy),rightiy(ix,iy),0) *
     &        qc(rightix(ix,iy),rightiy(ix,iy),0))
            wbbl(ix,iy,3) =
     &       sqrt(wbbl(ix,iy,0)**2+wbbl(ix,iy,1)**2+wbbl(ix,iy,2)**2)
          end if
          wrk = pol_flux /
     .               (hz(ix,iy)*qz(ix,iy,1)*bb(ix,iy,0))
          if (wrk .lt. 0.0_R8) then
            write(*,'(a,i2,a,i2,a)') "For cell (ix,iy) (", ix,",", iy, 
     &        ") redefined hy < 0!"
            cycle
          endif
          hy(ix,iy) = wrk
          hx(ix,iy) = vol(ix,iy)/
     .               (hz(ix,iy)*hy(ix,iy)*qz(ix,iy,1))
        end do
      end do

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
