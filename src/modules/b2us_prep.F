      module b2us_prep

      use b2mod_types
      use b2mod_constants
      use b2mod_geo
      use b2mod_geo_corner
      use b2mod_indirect
      use b2mod_plasma
      use b2mod_transport
      use b2mod_anomalous_transport
      use b2mod_sources
      use b2mod_work
      use b2mod_residuals
      use b2mod_input_profile
      use b2mod_rates
      use b2mod_diag
      use b2mod_time
      use b2mod_version
      use b2mod_b2cmpa
      use b2us_map
      use b2us_geo
      use b2us_plasma
      use b2mod_b2cmfs
      use b2mod_connectivity
      use carre_constants

      implicit none

      private

      public :: calc_mapping,
     .          prep_geometry,
     .          prep_plasma, prep_plasma_ext, prep_state,
     .          to_struct_plasma,
     .          conv_regions, conv_midplane                             !srv 16.11.22

      contains

***********************************************************************
*     Module containing basic routines to convert geomtry and state
*     information from structured to unstructured format, and vice
*     versa.
***********************************************************************

      subroutine calc_mapping(nx, ny, m)
#ifdef WG_DEBUG
      use b2us_debug
     .  , isStructured => isClassicalGrid
#endif
      implicit none
      type (mapping), intent (inout) :: m
      integer, intent(in) :: nx, ny
      integer :: ix, iy, i, if1, if2, if3, if4, ic, ifc,
     .           vxFc_end, vxCv_end, icut, ix1, ix2, iy1, iy2,
     .           ivx, ivx1, ivx2, mxreg, ift, icv, ifs,
     .           ifcc, icount, n2, ix2o, iy2o,
     .           ncorneri, nghostfcxi, nghostfcyi, ifcy,
     .           nghostvxi, ici, icg, ict, j, itmpx, itmpy, lix, liy,
     .           ixb, iyb, bix, biy, tix, tiy, rix, riy, ifccy,
     .           nyreg, ixreg, iyreg, is, nxtl, nxtr,
     .           iunused, nunused, ivxe, nvxe, lbl, ier, nfcxReg,
     .           istrtc, istrtfc, istrtfy, ilbl
      integer, allocatable :: indFc(:), indCv(:), fh(:), ch(:),
     .           icorner(:,:), ighostcv(:,:),
     .           ixghostcv(:), iyghostcv(:), ixcorner(:), iycorner(:),
     .           ighostfcx(:,:), ixghostfcx(:), iyghostfcx(:),
     .           ighostfcy(:,:), ixghostfcy(:), iyghostfcy(:),
     .           ighostvx(:,:), ixghostvx(:), iyghostvx(:),
!     .           locleftix(:,:), locleftiy(:,:),
     .           locrightix(:,:), locrightiy(:,:),
!     .           locbottomix(:,:), locbottomiy(:,:),
     .           loctopix(:,:), loctopiy(:,:),
     .           ixmapghostvx(:,:), iymapghostvx(:,:),
     .           ixvxe(:), iyvxe(:), ftFcP(:,:), ftCvP(:,:), fsFcP(:,:),
     .           ftReg(:), nxreg(:), reg(:,:), ixstart(:,:,:), ftLbl(:),
     .           ixend(:,:,:), iystart(:), iyend(:), nseg(:,:)
      real(R8) :: xco(4), yco(4), babs
      logical lfound, swapped, classical, bound, cellError
      intrinsic pack

!  assume a structured quadrilateral grid

!  store information for easy back-transformation to structured format
!WG-TODO: rename m%isClassicalGrid to "m%isStructuredGrid"?
      m%isClassicalGrid = 1
      m%nx    = nx
      m%ny    = ny
      m%nncut = nncut

      classical = isClassicalGrid(cflags)

      write (*,*) 'Creating mapping ...'

#ifdef WG_DEBUG
      nx_st = nx
      ny_st = ny
#endif

!  no. of cells
! Count internal ghost cells and "corner" points
! Remark in this mapping routine the old B2-numbering is used; do loop over internal cells run from 1->nx, etc.

! Initialize the mapping matrices
      n2 = (nx+2)*(ny+2)
      allocate (icorner(-1:nx,-1:ny))
      allocate (ighostcv(-1:nx,-1:ny))
      allocate (ixghostcv(n2))
      allocate (iyghostcv(n2))
      allocate (ixcorner(n2))
      allocate (iycorner(n2))
      allocate (ighostfcx(-1:nx,-1:ny))
      allocate (ixghostfcx(n2))
      allocate (iyghostfcx(n2))
      allocate (ighostfcy(-1:nx,-1:ny))
      allocate (ixghostfcy(n2))
      allocate (iyghostfcy(n2))
      allocate (ighostvx(-1:nx,-1:ny))
      allocate (ixghostvx(n2))
      allocate (iyghostvx(n2))
!      allocate (locleftix(-1:nx,-1:ny))
!      allocate (locleftiy(-1:nx,-1:ny))
      allocate (locrightix(-1:nx,-1:ny))
      allocate (locrightiy(-1:nx,-1:ny))
!      allocate (locbottomix(-1:nx,-1:ny))
!      allocate (locbottomiy(-1:nx,-1:ny))
      allocate (loctopix(-1:nx,-1:ny))
      allocate (loctopiy(-1:nx,-1:ny))
      allocate (ixmapghostvx(-1:nx,-1:ny))
      allocate (iymapghostvx(-1:nx,-1:ny))
      allocate (ixvxe(n2))
      allocate (iyvxe(n2))
      ighostcv = 0
      icorner = 0
      ixghostcv = 0
      iyghostcv = 0
      ixcorner = 0
      iycorner = 0
      ighostfcx = 0
      ixghostfcx = 0
      iyghostfcx = 0
      ighostfcy = 0
      ixghostfcy = 0
      iyghostfcy = 0
      ighostvx = 0
      ixghostvx = 0
      iyghostvx = 0
!      locleftix = leftix
!      locleftiy = leftiy
      locrightix = rightix
      locrightiy = rightiy
!      locbottomix = bottomix
!      locbottomiy = bottomiy
      loctopix = topix
      loctopiy = topiy
      ixmapghostvx = -2
      iymapghostvx = -2
      ixvxe = -2
      iyvxe = -2

! no. of cells
! Count internal ghost cells
      icount = 0
! Count "corner" cells, by checking whether it is surrounded by ghosts and/or
! non-existing cells
      do iy = -1, ny
       do ix = -1, nx
        if (.not.isGhostCell(cflags(ix, iy, CELLFLAG_TYPE))) cycle
        ict = 0
! Count inexisting neighbors
        if (rightix(ix,iy)==nx+1) ict = ict+1
        if (leftix(ix,iy)==-2  )  ict = ict+1
        if (topiy(ix,iy)==ny+1)   ict = ict+1
        if (bottomiy(ix,iy)==-2)  ict = ict+1
        ! Count guard neighbors
        if (rightix(ix,iy)/=nx+1) then
         if (cflags(rightix(ix,iy),rightiy(ix,iy),CELLFLAG_TYPE)==
     .       GRID_GUARD) ict = ict+1
        endif
        if (leftix(ix,iy)/=-2) then
         if (cflags(leftix(ix,iy),leftiy(ix,iy),CELLFLAG_TYPE)==
     .       GRID_GUARD) ict = ict+1
        endif
        if (topiy(ix,iy)/=ny+1) then
         if (cflags(topix(ix,iy),topiy(ix,iy),CELLFLAG_TYPE)==
     .       GRID_GUARD) ict = ict+1
        endif
        if (bottomiy(ix,iy)/=-2) then
         if (cflags(bottomix(ix,iy),bottomiy(ix,iy),CELLFLAG_TYPE)==
     .       GRID_GUARD) ict = ict+1
        endif
        ! Corner if all neighbors are either guard or non-existing
        if (ict == 4) then
         icorner(ix,iy)=1
         icount=icount+1
         ixcorner(icount) = ix
         iycorner(icount) = iy
        endif
       enddo
      enddo
      ncorneri = icount
      m%nCi = count( cflags(:, :, CELLFLAG_TYPE) == GRID_INTERNAL)
     .       + count( cflags(:, :, CELLFLAG_TYPE) == GRID_BOUNDARY)
      m%nCg = count( cflags(:, :, CELLFLAG_TYPE) == GRID_GUARD )
     .       - ncorneri
      m%nCv = m%nCi + m%nCg             ! total number of cells

      write (*,*) 'nCv, nCi, nCg', m%nCv, m%nCi, m%nCg
      write (*,*) 'ncorneri     ', ncorneri

! find position of central cut
      if (classical) then
       iy = 1
       m%nctcnt = count(leftix(0:nx-1,iy-1) == -2)
       allocate (m%ictcntrl(m%nctcnt))
       m%ictcntrl = 0
       ict = 0
       do ix = 1, nx
         if (leftix(ix-1,iy-1) == -2) then
           ict = ict+1
           m%ictcntrl(ict) = ix-1
         end if
       end do
      endif

! no. of faces
! Count internal faces in between ghost cells (first x-faces then y-faces) and degenerate faces from triangles
      icount = 0
      do iy = 0, ny
        do ix = 0, nx
          if (leftix(ix,iy).eq.-2) then
            ! W guard cells, unused cells, degenerate triangle faces
            ighostfcx(ix,iy)=1
            icount=icount+1
            ixghostfcx(icount) = ix
            iyghostfcx(icount) = iy
          elseif(.not.isRealCell(cflags(ix,iy,CELLFLAG_TYPE)).and.
     .           .not.isRealCell(cflags(leftix(ix,iy),leftiy(ix,iy),
     .                                 CELLFLAG_TYPE)))then
            ! Faces between guard cells
            ighostfcx(ix,iy)=1
            icount=icount+1
            ixghostfcx(icount) = ix
            iyghostfcx(icount) = iy
          endif
        enddo
      enddo
      nghostfcxi = icount

      icount = 0
      do ix = 0, nx
        do iy = 0, ny
          if (bottomiy(ix,iy).eq.-2) then
            ! S guard cells, unused cells, degenerate triangle faces
            ighostfcy(ix,iy)=1
            icount=icount+1
            ixghostfcy(icount) = ix
            iyghostfcy(icount) = iy
          elseif(.not.isRealCell(cflags(ix,iy,CELLFLAG_TYPE)).and.
     .           .not.isRealCell(cflags(bottomix(ix,iy),bottomiy(ix,iy),
     .                                 CELLFLAG_TYPE)))then
            ! Faces between guard cells
            ighostfcy(ix,iy)=1
            icount=icount+1
            ixghostfcy(icount) = ix
            iyghostfcy(icount) = iy
          endif
        enddo
      enddo
      nghostfcyi = icount
      m%nFc = 2*(ny+1)*(nx+1) - nghostfcxi - nghostfcyi

! no. of vertices
! Count vertices to be removed
      iunused = 0
      icount  = 0
      nvxe    = 0

! First pass: remove vertices of all unused cells
      do iy = -1, ny
        do ix = -1, nx
          if (isUnusedCell(cflags(ix,iy,CELLFLAG_TYPE))) then
            iunused = iunused + 1
            ighostvx(ix,iy) = 1
            ixmapghostvx(ix,iy) = -2
            iymapghostvx(ix,iy) = -2
          endif
        enddo
      enddo
      nunused = iunused

! Second pass: add external vertices if needed
      do iy = -1, ny
        do ix = -1, nx
          if (cflags(ix,iy,CELLFLAG_TYPE) == GRID_INTERNAL) then
            ! Regular internal cell
            ! Need to add a new, external vertex if both
            ! E a boundary cell with degenerate N face and
            ! N a boundary cell with degenerate E face
            tix = topix(ix,iy)
            tiy = topiy(ix,iy)
            rix = rightix(ix,iy)
            riy = rightiy(ix,iy)
            if (isBoundaryCell(cflags(tix,tiy,CELLFLAG_TYPE)) .and.
     .           points_match(crx(tix,tiy,1), cry(tix,tiy,1),
     .                        crx(tix,tiy,3), cry(tix,tiy,3)) .and.
     .          isBoundaryCell(cflags(rix,riy,CELLFLAG_TYPE)) .and.
     .           points_match(crx(rix,riy,2), cry(rix,riy,2),
     .                        crx(rix,riy,3), cry(rix,riy,3))) then
              ! Check whether an external vertex has already been created here
              lfound = .false.
              ivxe   = 0
              do while (.not.lfound .and. ivxe.lt.nvxe)
                ivxe = ivxe + 1
                lfound = lfound.or.
     .            points_match(crx(ix,iy,3), cry(ix,iy,3),
     .                         crx(ixvxe(ivxe),iyvxe(ivxe),0),
     .                         cry(ixvxe(ivxe),iyvxe(ivxe),0))
              enddo
              if (.not.lfound) then
                ! Store in an available external cell
                nvxe = nvxe + 1
                ix2  = ix + 1
                iy2  = iy
                do while ((ighostvx(ix2,iy2) == 0) .or.
     .           .not.isUnusedCell(cflags(ix2,iy2,CELLFLAG_TYPE)))
                  ix2 = ix2 + 1
                  if (ix2 == nx+1) then
                    ix2 = ix2-1
                    do while ((ighostvx(ix2,iy2) == 0) .or.
     .                .not.isUnusedCell
     .                   (cflags(ix2,iy2,CELLFLAG_TYPE)))
                      ix2 = ix2 - 1
                      if (ix2 == 0) then
                        write (*,*) 'ix,  iy  = ', ix, iy
                        write (*,*) 'ix2, iy2 = ', ix2, iy2
                        call xerrab
     .                    ('No suitable unused cell found (0).')
                      endif
                    enddo
                  endif
                enddo
                ! Create vertex, assign relevant data
                ixvxe(nvxe) = ix2
                iyvxe(nvxe) = iy2
                nunused = nunused - 1
                ighostvx(ix2,iy2) = 0
                crx(ix2,iy2,0)    = crx(ix,iy,3)
                cry(ix2,iy2,0)    = cry(ix,iy,3)
                ffbz(ix2,iy2,0)   = ffbz(ix,iy,3)
                fpsi(ix2,iy2,0,:) = fpsi(ix,iy,3,:)
!WG_TODO need to recompute field at corner here;
!WG_TODO tmp: assume isymm = 1 or 2 treatment... (see carre_b2ag.F90)
!WG_TODO eventually: do this in separate routine later on?
                if (isymm.ne.1.and.isymm.ne.2)
     .            call xerrab('isymm value not yet available')
                babs=sqrt(fpsi(ix,iy,3,1)**2+fpsi(ix,iy,3,2)**2)/
     .            (2.0*pi*crx(ix,iy,3))
                wbbc(ix2,iy2,0)   = sign(babs,bb(ix,iy,0))
                wbbc(ix2,iy2,1)   = 0.0_R8
                wbbc(ix2,iy2,2)   = ffbz(ix,iy,3)/
     .            (2.0*pi*crx(ix,iy,3))
                wbbc(ix2,iy2,3)   = sqrt(wbbc(ix2,iy2,0)**2 +
     .            wbbc(ix2,iy2,1)**2 + wbbc(ix2,iy2,2)**2)
                ! adapt local connectivity arrays
                locrightix(tix,tiy) = ixvxe(nvxe)
                locrightiy(tix,tiy) = iyvxe(nvxe)
                loctopix(rix,riy)   = ixvxe(nvxe)
                loctopiy(rix,riy)   = iyvxe(nvxe)
                write (*,*)
     .           '(1) Created vertex at ix, iy, crx, cry = ',
     .            nvxe, ix, iy, crx(ix2,iy2,0), cry(ix2,iy2,0)
              else
                ! adapt local connectivity arrays
                locrightix(tix,tiy) = ixvxe(ivxe)
                locrightiy(tix,tiy) = iyvxe(ivxe)
                loctopix(rix,riy)   = ixvxe(ivxe)
                loctopiy(rix,riy)   = iyvxe(ivxe)
              endif
            endif
          elseif (isBoundaryCell(cflags(ix,iy,CELLFLAG_TYPE))) then
            ! For boundary cells: keep vertex, except if points coincide
            if (points_match(crx(ix,iy,0), cry(ix,iy,0),
     .                       crx(ix,iy,2), cry(ix,iy,2))) then
              ! Boundary cell with degenerate W face
              ! Merge vertex with top neighbor
              icount = icount + 1
              ighostvx(ix,iy)   = 1
              ixmapghostvx(ix,iy) = topix(ix,iy)
              iymapghostvx(ix,iy) = topiy(ix,iy)
            elseif (points_match(crx(ix,iy,0), cry(ix,iy,0),
     .                           crx(ix,iy,1), cry(ix,iy,1))) then
              ! Boundary cell with degenerate S face
              ! Merge vertex with right neighbor
              icount = icount + 1
              ighostvx(ix,iy)   = 1
              ixmapghostvx(ix,iy) = rightix(ix,iy)
              iymapghostvx(ix,iy) = rightiy(ix,iy)
            elseif (points_match(crx(ix,iy,1), cry(ix,iy,1),
     .                           crx(ix,iy,3), cry(ix,iy,3))) then
              ! Boundary cell with degenerate E face
              ! In 3 cases, need to add a new, external vertex:
              ! 1) if both N+S neighbors are guard
              ! 2) if S is guard, and N a boundary cell with degenerate E face
              ! 3) if N is guard, and S a boundary cell with degenerate E face
              tix = topix(ix,iy)
              tiy = topiy(ix,iy)
              bix = bottomix(ix,iy)
              biy = bottomiy(ix,iy)
              if ((isGhostCell(cflags(tix,tiy,CELLFLAG_TYPE)) .and.
     .             isGhostCell(cflags(bix,biy,CELLFLAG_TYPE))) .or.
     .            (isGhostCell(cflags(tix,tiy,CELLFLAG_TYPE)) .and.
     .             isBoundaryCell(cflags(bix,biy,CELLFLAG_TYPE)) .and.
     .             points_match(crx(bix,biy,1), cry(bix,biy,1),
     .                          crx(bix,biy,3), cry(bix,biy,3))) .or.
     .            (isGhostCell(cflags(bix,biy,CELLFLAG_TYPE)) .and.
     .             isBoundaryCell(cflags(tix,tiy,CELLFLAG_TYPE)) .and.
     .             points_match(crx(tix,tiy,1), cry(tix,tiy,1),
     .                          crx(tix,tiy,3), cry(tix,tiy,3)))) then
                ! Check whether an external vertex has already been created here
                lfound = .false.
                ivxe   = 0
                do while (.not.lfound .and. ivxe.lt.nvxe)
                  ivxe = ivxe + 1
                  lfound = lfound.or.
     .              points_match(crx(ix,iy,3), cry(ix,iy,3),
     .                           crx(ixvxe(ivxe),iyvxe(ivxe),0),
     .                           cry(ixvxe(ivxe),iyvxe(ivxe),0))
                enddo
                if (.not.lfound) then
                  ! Store in an available external cell
                  nvxe = nvxe + 1
                  ix2  = ix + 1
                  iy2  = iy
                  do while ((ighostvx(ix2,iy2) == 0) .or.
     .             .not.isUnusedCell(cflags(ix2,iy2,CELLFLAG_TYPE)))
                    ix2 = ix2 + 1
                    if (ix2 == nx+1) then
                      ix2 = ix2-1
                      do while ((ighostvx(ix2,iy2) == 0) .or.
     .                  .not.isUnusedCell
     .                     (cflags(ix2,iy2,CELLFLAG_TYPE)))
                        ix2 = ix2 - 1
                        if (ix2 == 0) then
                          write (*,*) 'ix,  iy  = ', ix, iy
                          write (*,*) 'ix2, iy2 = ', ix2, iy2
                          call xerrab
     .                      ('No suitable unused cell found (1).')
                        endif
                      enddo
                    endif
                  enddo
                  ! Create vertex, assign relevant data
                  ixvxe(nvxe) = ix2
                  iyvxe(nvxe) = iy2
                  nunused = nunused - 1
                  ighostvx(ix2,iy2) = 0
                  crx(ix2,iy2,0)    = crx(ix,iy,3)
                  cry(ix2,iy2,0)    = cry(ix,iy,3)
                  ffbz(ix2,iy2,0)   = ffbz(ix,iy,3)
                  fpsi(ix2,iy2,0,:) = fpsi(ix,iy,3,:)
!WG_TODO need to recompute field at corner here;
!WG_TODO tmp: assume isymm = 1 or 2 treatment... (see carre_b2ag.F90)
!WG_TODO eventually: do this in separate routine later on?
                  if (isymm.ne.1.and.isymm.ne.2)
     .              call xerrab('isymm value not yet available')
                  babs=sqrt(fpsi(ix,iy,3,1)**2+fpsi(ix,iy,3,2)**2)/
     .              (2.0*pi*crx(ix,iy,3))
                  wbbc(ix2,iy2,0)   = sign(babs,bb(ix,iy,0))
                  wbbc(ix2,iy2,1)   = 0.0_R8
                  wbbc(ix2,iy2,2)   = ffbz(ix,iy,3)/
     .              (2.0*pi*crx(ix,iy,3))
                  wbbc(ix2,iy2,3)   = sqrt(wbbc(ix2,iy2,0)**2 +
     .              wbbc(ix2,iy2,1)**2 + wbbc(ix2,iy2,2)**2)
                  ! adapt local connectivity arrays
                  locrightix(ix,iy) = ixvxe(nvxe)
                  locrightiy(ix,iy) = iyvxe(nvxe)
                  if (isGhostCell(cflags(tix,tiy,CELLFLAG_TYPE))) then
                    locrightix(tix,tiy) = ixvxe(nvxe)
                    locrightiy(tix,tiy) = iyvxe(nvxe)
                  endif
                  write (*,*)
     .             '(1) Created vertex at ix, iy, crx, cry = ',
     .              nvxe, ix, iy, crx(ix2,iy2,0), cry(ix2,iy2,0)
                else
                  ! adapt local connectivity arrays
                  locrightix(ix,iy) = ixvxe(ivxe)
                  locrightiy(ix,iy) = iyvxe(ivxe)
                  if (isGhostCell(cflags(tix,tiy,CELLFLAG_TYPE))) then
                    locrightix(tix,tiy) = ixvxe(ivxe)
                    locrightiy(tix,tiy) = iyvxe(ivxe)
                  endif
                endif
              elseif (isGhostCell(cflags(tix,tiy,CELLFLAG_TYPE))) then
                ! adapt local connectivity arrays
                locrightix(ix,iy) = loctopix(locrightix(bix,biy),
     .                                        locrightiy(bix,biy))
                locrightiy(ix,iy) = loctopiy(locrightix(bix,biy),
     .                                       locrightiy(bix,biy))
                locrightix(tix,tiy) = loctopix(locrightix(bix,biy),
     .                                         locrightiy(bix,biy))
                locrightiy(tix,tiy) = loctopiy(locrightix(bix,biy),
     .                                         locrightiy(bix,biy))
              elseif (isGhostCell(cflags(bix,biy,CELLFLAG_TYPE))) then
                ! adapt local connectivity arrays
                locrightix(ix,iy) = locrightix(tix,tiy)
                locrightiy(ix,iy) = locrightiy(tix,tiy)
              endif
            elseif (points_match(crx(ix,iy,2), cry(ix,iy,2),
     .                           crx(ix,iy,3), cry(ix,iy,3))) then
              ! Boundary cell with degenerate N face
              ! In 3 cases, need to add a new, external vertex:
              ! 1) if both E+W neighbors are guard
              ! 2) if E is guard, and W a boundary cell with degenerate N face
              ! 3) if W is guard, and E a boundary cell with degenerate N face
              lix = leftix(ix,iy)
              liy = leftiy(ix,iy)
              rix = rightix(ix,iy)
              riy = rightiy(ix,iy)
              if ((isGhostCell(cflags(lix,liy,CELLFLAG_TYPE)) .and.
     .             isGhostCell(cflags(rix,riy,CELLFLAG_TYPE))) .or.
     .            (isGhostCell(cflags(lix,liy,CELLFLAG_TYPE)) .and.
     .             isBoundaryCell(cflags(rix,riy,CELLFLAG_TYPE)) .and.
     .             points_match(crx(rix,riy,2), cry(rix,riy,2),
     .                          crx(rix,riy,3), cry(rix,riy,3))) .or.
     .            (isGhostCell(cflags(rix,riy,CELLFLAG_TYPE)) .and.
     .             isBoundaryCell(cflags(lix,liy,CELLFLAG_TYPE)) .and.
     .             points_match(crx(lix,liy,2), cry(lix,liy,2),
     .                          crx(lix,liy,3), cry(lix,liy,3)))) then
                ! Check whether an external vertex has already been created here
                lfound = .false.
                ivxe   = 0
                do while (.not.lfound .and. ivxe.lt.nvxe)
                  ivxe = ivxe + 1
                  lfound = lfound.or.
     .              points_match(crx(ix,iy,3), cry(ix,iy,3),
     .                           crx(ixvxe(ivxe),iyvxe(ivxe),0),
     .                           cry(ixvxe(ivxe),iyvxe(ivxe),0))
                enddo
                if (.not.lfound) then
                  ! Store in an available external cell
                  nvxe = nvxe + 1
                  ix2  = ix
                  iy2  = iy + 1
                  do while ((ighostvx(ix2,iy2) == 0) .or.
     .             .not.isUnusedCell(cflags(ix2,iy2,CELLFLAG_TYPE)))
                    iy2 = iy2 + 1
                    if (iy2 == ny+1) then
                      iy2 = iy2 - 1
                      do while ((ighostvx(ix2,iy2) == 0) .or.
     .                  .not.isUnusedCell
     .                      (cflags(ix2,iy2,CELLFLAG_TYPE)))
                        iy2 = iy2 - 1
                        if (iy2 == 0) then
                          write (*,*) 'ix,  iy  = ', ix, iy
                          write (*,*) 'ix2, iy2 = ', ix2, iy2
                          call xerrab
     .                      ('No suitable unused cell found (2).')
                        endif
                      enddo
                    endif
                  enddo
                  ! Create vertex, assign relevant data
                  ixvxe(nvxe) = ix2
                  iyvxe(nvxe) = iy2
                  nunused = nunused - 1
                  ighostvx(ix2,iy2) = 0
                  crx(ix2,iy2,0)    = crx(ix,iy,3)
                  cry(ix2,iy2,0)    = cry(ix,iy,3)
                  ffbz(ix2,iy2,0)   = ffbz(ix,iy,3)
                  fpsi(ix2,iy2,0,:) = fpsi(ix,iy,3,:)
!WG_TODO need to recompute field at corner here;
!WG_TODO tmp: assume isymm = 1 or 2 treatment... (see carre_b2ag.F90)
!WG_TODO eventually: do this in separate routine later on?
                  if (isymm.ne.1.and.isymm.ne.2)
     .              call xerrab('isymm value not yet available')
                  babs=sqrt(fpsi(ix,iy,3,1)**2+fpsi(ix,iy,3,2)**2)/
     .              (2.0*pi*crx(ix,iy,3))
                  wbbc(ix2,iy2,0)   = sign(babs,bb(ix,iy,0))
                  wbbc(ix2,iy2,1)   = 0.0_R8
                  wbbc(ix2,iy2,2)   = ffbz(ix,iy,3)/
     .              (2.0*pi*crx(ix,iy,3))
                  wbbc(ix2,iy2,3)   = sqrt(wbbc(ix2,iy2,0)**2 +
     .              wbbc(ix2,iy2,1)**2 + wbbc(ix2,iy2,2)**2)
                  ! adapt local connectivity arrays
                  loctopix(ix,iy) = ixvxe(nvxe)
                  loctopiy(ix,iy) = iyvxe(nvxe)
                  if (isGhostCell(cflags(rix,riy,CELLFLAG_TYPE))) then
                    loctopix(rix,riy) = ixvxe(nvxe)
                    loctopiy(rix,riy) = iyvxe(nvxe)
                  endif
                  write (*,*)
     .             '(2) Created vertex at ix, iy, crx, cry = ',
     .              nvxe, ix, iy, crx(ix2,iy2,0), cry(ix2,iy2,0)
                else
                  ! adapt local connectivity arrays
                  loctopix(ix,iy) = ixvxe(ivxe)
                  loctopiy(ix,iy) = iyvxe(ivxe)
                  if (isGhostCell(cflags(rix,riy,CELLFLAG_TYPE))) then
                    loctopix(rix,riy) = ixvxe(ivxe)
                    loctopiy(rix,riy) = iyvxe(ivxe)
                  endif
                endif
              elseif (isGhostCell(cflags(rix,riy,CELLFLAG_TYPE))) then
                ! adapt local connectivity arrays
                loctopix(ix,iy) = locrightix(loctopix(lix,liy),
     .                                       loctopiy(lix,liy))
                loctopiy(ix,iy) = locrightiy(loctopix(lix,liy),
     .                                       loctopiy(lix,liy))
                loctopix(rix,riy) = locrightix(loctopix(lix,liy),
     .                                         loctopiy(lix,liy))
                loctopiy(rix,riy) = locrightiy(loctopix(lix,liy),
     .                                         loctopiy(lix,liy))
              elseif (isGhostCell(cflags(lix,liy,CELLFLAG_TYPE))) then
                ! adapt local connectivity arrays
                loctopix(ix,iy) = loctopix(rix,riy)
                loctopiy(ix,iy) = loctopiy(rix,riy)
              endif
            else
              tix = topix(ix,iy)
              tiy = topiy(ix,iy)
              rix = rightix(ix,iy)
              riy = rightiy(ix,iy)
              if (isGhostCell(cflags(tix,tiy,CELLFLAG_TYPE)) .and.
     .            isGhostCell(cflags(rix,riy,CELLFLAG_TYPE))) then
                ! Corner; check if top-right vertex exists
                ! If not, add external vertex
                if (rightix(tix,tiy)==nx+1) then
                    ! Check whether an external vertex has already been created here
                    lfound = .false.
                    ivxe   = 0
                    do while (.not.lfound .and. ixvxe(ivxe+1) /= -2)
                      ivxe = ivxe + 1
                      lfound = lfound.or.
     .                  points_match(crx(ix,iy,3), cry(ix,iy,3),
     .                               crx(ixvxe(ivxe),iyvxe(ivxe),0),
     .                               cry(ixvxe(ivxe),iyvxe(ivxe),0))
                    enddo
                    if (.not.lfound) then
                      ! Store in an available external cell
                      nvxe = nvxe + 1
                      ix2  = tix + 1
                      iy2  = tiy
                      do while ((ighostvx(ix2,iy2) == 0) .or.
     .                 .not.isUnusedCell(cflags(ix2,iy2,CELLFLAG_TYPE)))
                        iy2 = iy2 + 1
                        if (iy2 == ny+1) then
                          iy2 = iy2 - 1
                          do while ((ighostvx(ix2,iy2) == 0) .or.
     .                      .not.isUnusedCell
     .                          (cflags(ix2,iy2,CELLFLAG_TYPE)))
                            iy2 = iy2 - 1
                            if (iy2 == 0) then
                              ix2 = tix + 1
                              iy2 = tiy
                              do while ((ighostvx(ix2,iy2) == 0) .or.
     .                          .not.isUnusedCell
     .                              (cflags(ix2,iy2,CELLFLAG_TYPE)))
                                ix2 = ix2 + 1
                                if (ix2 == nx+1) then
                                  write (*,*) 'ix,  iy  = ', ix, iy
                                  write (*,*) 'ix2, iy2 = ', ix2, iy2
                                  call xerrab
     .                          ('No suitable unused cell found (3).')
                                endif
                              enddo
                            endif
                          enddo
                        endif
                      enddo
                      ! Create vertex, assign relevant data
                      ixvxe(nvxe) = ix2
                      iyvxe(nvxe) = iy2
                      nunused = nunused - 1
                      ighostvx(ix2,iy2) = 0
                      crx(ix2,iy2,0)    = crx(ix,iy,2)
                      cry(ix2,iy2,0)    = cry(ix,iy,2)
                      ffbz(ix2,iy2,0)   = ffbz(ix,iy,2)
                      fpsi(ix2,iy2,0,:) = fpsi(ix,iy,2,:)
!WG_TODO need to recompute field at corner here;
!WG_TODO tmp: assume isymm = 1 or 2 treatment... (see carre_b2ag.F90)
!WG_TODO eventually: do this in separate routine later on?
                      if (isymm.ne.1.and.isymm.ne.2)
     .                  call xerrab('isymm value not yet available')
                      babs=sqrt(fpsi(ix,iy,2,1)**2+fpsi(ix,iy,2,2)**2)/
     .                  (2.0*pi*crx(ix,iy,2))
                      wbbc(ix2,iy2,0)   = sign(babs,bb(ix,iy,0))
                      wbbc(ix2,iy2,1)   = 0.0_R8
                      wbbc(ix2,iy2,2)   = ffbz(ix,iy,2)/
     .                  (2.0*pi*crx(ix,iy,2))
                      wbbc(ix2,iy2,3)   = sqrt(wbbc(ix2,iy2,0)**2 +
     .                  wbbc(ix2,iy2,1)**2 + wbbc(ix2,iy2,2)**2)
                      ! Adapt local connectivity arrays
                      locrightix(tix,tiy) = ixvxe(nvxe)
                      locrightiy(tix,tiy) = iyvxe(nvxe)
                      loctopix(rix,riy)   = ixvxe(nvxe)
                      loctopiy(rix,riy)   = iyvxe(nvxe)
                      write (*,*)
     .                 '(3) Created vertex at ix, iy, crx, cry = ',
     .                  nvxe, ix, iy, crx(ix2,iy2,0), cry(ix2,iy2,0)
                    else
                      ! Adapt local connectivity arrays
                      locrightix(tix,tiy) = ixvxe(ivxe)
                      locrightiy(tix,tiy) = iyvxe(ivxe)
                      loctopix(rix,riy)   = ixvxe(ivxe)
                      loctopiy(rix,riy)   = iyvxe(ivxe)
                    endif
                else
                  ! Adapt local connectivity arrays
                  loctopix(rix,riy)   = rightix(tix,tiy)
                  loctopiy(rix,riy)   = rightiy(tix,tiy)
                endif
              endif
            endif
          elseif (isGhostCell(cflags(ix,iy,CELLFLAG_TYPE))) then
            if(icorner(ix,iy) == 0) then
              ! Remove vertex of S/W cells
              tix = topix(ix,iy)
              tiy = topiy(ix,iy)
              if (tiy/=ny+1) then
                if (isBoundaryCell(cflags(tix,tiy,CELLFLAG_TYPE))) then
                  ! (ix,iy) is S guard cell; remove
                  icount = icount + 1
                  ighostvx(ix,iy)   = 1
                  ixghostvx(icount) = ix
                  iyghostvx(icount) = iy
                  ixmapghostvx(ix,iy) = -2
                  iymapghostvx(ix,iy) = -2
                endif
              endif
              rix = rightix(ix,iy)
              riy = rightiy(ix,iy)
              if (rix/=nx+1) then
                if (isBoundaryCell(cflags(rix,riy,CELLFLAG_TYPE))) then
                  ! (ix,iy) is W guard cell; remove
                  icount = icount + 1
                  ighostvx(ix,iy)   = 1
                  ixghostvx(icount) = ix
                  iyghostvx(icount) = iy
                  ixmapghostvx(ix,iy) = -2
                  iymapghostvx(ix,iy) = -2
                endif
              endif
              bix = bottomix(ix,iy)
              biy = bottomiy(ix,iy)
              if (biy/=-2) then
                if (isBoundaryCell(cflags(bix,biy,CELLFLAG_TYPE))) then
                  ! For N guard cells, check if vertex not already added as an E
                  ! (special case of corners protruding into the plasma)
                  ! If so, link to the other vertex
                  lix = leftix(bix,biy)
                  liy = leftiy(bix,biy)
                  if (lix/=-2) then
                    tix = topix(lix,liy)
                    tiy = topiy(lix,liy)
                    if (tiy/=ny+1) then
                      if (.not.isGhostCell(cflags(tix,tiy,
     .                         CELLFLAG_TYPE))) then
                      rix = rightix(tix,tiy)
                      riy = rightiy(tix,tiy)
                      if (rix/=nx+1) then
                        icount = icount + 1
                        ighostvx(ix,iy)   = 1
                        ixghostvx(icount) = ix
                        iyghostvx(icount) = iy
                        ixmapghostvx(ix,iy) = rix
                        iymapghostvx(ix,iy) = riy
                      endif
                      endif
                    endif
                  endif
                endif
              endif
              lix = leftix(ix,iy)
              liy = leftiy(ix,iy)
              if (lix/=-2) then
                if (isBoundaryCell(cflags(lix,liy,CELLFLAG_TYPE))) then
                  ! For E guard cells, check if top vertex exists (corners in the geometry)
                  ! If not, add external vertex
                  lfound=.false.
                  tix = topix(lix,liy)
                  tiy = topiy(lix,liy)
                  if (tiy/=ny+1) then
                    rix = rightix(tix,tiy)
                    riy = rightiy(tix,tiy)
                    if (rix/=nx+1) lfound=.true.
                  elseif (topiy(leftix(lix,liy),
     .                          leftiy(lix,liy))/=ny+1) then
                    tix = topix(leftix(lix,liy),leftiy(lix,liy))
                    tiy = topiy(leftix(lix,liy),leftiy(lix,liy))
                    rix = rightix(tix,tiy)
                    riy = rightiy(tix,tiy)
                    if (rix/=nx+1) lfound=.true.
                  endif
                  if (.not.lfound) then
                    ! Check whether an external vertex has already been created here
                    ivxe   = 0
                    do while (.not.lfound .and. ivxe.lt.nvxe)
                      ivxe = ivxe + 1
                      lfound = lfound.or.
     .                  points_match(crx(ix,iy,3), cry(ix,iy,3),
     .                               crx(ixvxe(ivxe),iyvxe(ivxe),0),
     .                               cry(ixvxe(ivxe),iyvxe(ivxe),0))
                    enddo
                    if (.not.lfound) then
                      ! Store in an available external cell
                      nvxe = nvxe + 1
                      ix2  = ix
                      iy2  = iy + 1
                      do while ((ighostvx(ix2,iy2) == 0) .or.
     .                 .not.isUnusedCell(cflags(ix2,iy2,CELLFLAG_TYPE)))
                        iy2 = iy2 + 1
                        if (iy2 == ny+1) then
                          iy2 = iy2 - 1
                          do while ((ighostvx(ix2,iy2) == 0) .or.
     .                      .not.isUnusedCell
     .                          (cflags(ix2,iy2,CELLFLAG_TYPE)))
                            iy2 = iy2 - 1
                            if (iy2 == 0) then
                              ix2 = ix
                              iy2 = iy + 1
                              do while ((ighostvx(ix2,iy2) == 0) .or.
     .                         .not.isUnusedCell
     .                             (cflags(ix2,iy2,CELLFLAG_TYPE)))
                                ix2 = ix2 + 1
                                if (ix2 == nx+1) then
                                  write (*,*) 'ix,  iy  = ', ix, iy
                                  write (*,*) 'ix2, iy2 = ', ix2, iy2
                                  call xerrab
     .                          ('No suitable unused cell found (4).')
                                endif
                              enddo
                            endif
                          enddo
                        endif
                      enddo
                      ! Create vertex, assign relevant data
                      ixvxe(nvxe) = ix2
                      iyvxe(nvxe) = iy2
                      nunused = nunused - 1
                      ighostvx(ix2,iy2) = 0
                      crx(ix2,iy2,0)    = crx(ix,iy,2)
                      cry(ix2,iy2,0)    = cry(ix,iy,2)
                      ffbz(ix2,iy2,0)   = ffbz(ix,iy,2)
                      fpsi(ix2,iy2,0,:) = fpsi(ix,iy,2,:)
!WG_TODO need to recompute field at corner here;
!WG_TODO tmp: assume isymm = 1 or 2 treatment... (see carre_b2ag.F90)
!WG_TODO eventually: do this in separate routine later on?
                      if (isymm.ne.1.and.isymm.ne.2)
     .                  call xerrab('isymm value not yet available')
                      babs=sqrt(fpsi(ix,iy,2,1)**2+fpsi(ix,iy,2,2)**2)/
     .                  (2.0*pi*crx(ix,iy,2))
                      wbbc(ix2,iy2,0)   = sign(babs,bb(ix,iy,0))
                      wbbc(ix2,iy2,1)   = 0.0_R8
                      wbbc(ix2,iy2,2)   = ffbz(ix,iy,2)/
     .                  (2.0*pi*crx(ix,iy,2))
                      wbbc(ix2,iy2,3)   = sqrt(wbbc(ix2,iy2,0)**2 +
     .                  wbbc(ix2,iy2,1)**2 + wbbc(ix2,iy2,2)**2)
                      ! Adapt local connectivity arrays
                      loctopix(ix,iy) = ixvxe(nvxe)
                      loctopiy(ix,iy) = iyvxe(nvxe)
                      write (*,*)
     .                 '(4) Created vertex at ix, iy, crx, cry = ',
     .                  nvxe, ix, iy, crx(ix2,iy2,0), cry(ix2,iy2,0)
                    endif
                  endif
                endif
              endif
            else
              ! For corners, only keep vertex of NE corners
              ! Current assumption: CARRE2/b2ag will not create corners for extended cases
              if ((leftix(ix,iy)   == -2)  .or.
     .            (bottomiy(ix,iy) == -2)) then
                ! (ix,iy) is not a NE corner cell; remove
                icount = icount + 1
                ighostvx(ix,iy)   = 1
                ixghostvx(icount) = ix
                iyghostvx(icount) = iy
                ixmapghostvx(ix,iy) = -2
                iymapghostvx(ix,iy) = -2
              endif
            endif
          endif
        enddo
      enddo
      nghostvxi = nunused + icount

! Third pass: adapt local connectivity arrays of boundary cells
      do iy = -1, ny
        do ix = -1, nx
          if (isBoundaryCell(cflags(ix,iy,CELLFLAG_TYPE))) then
            if (points_match(crx(ix,iy,1), cry(ix,iy,1),
     .                       crx(ix,iy,3), cry(ix,iy,3))) then
              ! Boundary cell with degenerate E face
              tix = topix(ix,iy)
              tiy = topiy(ix,iy)
              bix = bottomix(ix,iy)
              biy = bottomiy(ix,iy)
              if ((isGhostCell(cflags(tix,tiy,CELLFLAG_TYPE)) .and.
     .             isGhostCell(cflags(bix,biy,CELLFLAG_TYPE))) .or.
     .            (isGhostCell(cflags(tix,tiy,CELLFLAG_TYPE)) .and.
     .             isBoundaryCell(cflags(bix,biy,CELLFLAG_TYPE)) .and.
     .             points_match(crx(bix,biy,1), cry(bix,biy,1),
     .                          crx(bix,biy,3), cry(bix,biy,3))) .or.
     .            (isGhostCell(cflags(bix,biy,CELLFLAG_TYPE)) .and.
     .             isBoundaryCell(cflags(tix,tiy,CELLFLAG_TYPE)) .and.
     .             points_match(crx(tix,tiy,1), cry(tix,tiy,1),
     .                          crx(tix,tiy,3), cry(tix,tiy,3)))) then
              elseif (isGhostCell(cflags(tix,tiy,CELLFLAG_TYPE))) then
                ! adapt local connectivity arrays
                locrightix(ix,iy) = loctopix(rightix(bix,biy),
     .                                       rightiy(bix,biy))
                locrightiy(ix,iy) = loctopiy(rightix(bix,biy),
     .                                       rightiy(bix,biy))
                locrightix(tix,tiy) = loctopix(rightix(bix,biy),
     .                                         rightiy(bix,biy))
                locrightiy(tix,tiy) = loctopiy(rightix(bix,biy),
     .                                         rightiy(bix,biy))
              elseif (isGhostCell(cflags(bix,biy,CELLFLAG_TYPE))) then
                ! adapt local connectivity arrays
                locrightix(ix,iy) = rightix(tix,tiy)
                locrightiy(ix,iy) = rightiy(tix,tiy)
              endif
            elseif (points_match(crx(ix,iy,2), cry(ix,iy,2),
     .                           crx(ix,iy,3), cry(ix,iy,3))) then
              ! Boundary cell with degenerate N face
              lix = leftix(ix,iy)
              liy = leftiy(ix,iy)
              rix = rightix(ix,iy)
              riy = rightiy(ix,iy)
              if ((isGhostCell(cflags(lix,liy,CELLFLAG_TYPE)) .and.
     .             isGhostCell(cflags(rix,riy,CELLFLAG_TYPE))) .or.
     .            (isGhostCell(cflags(lix,liy,CELLFLAG_TYPE)) .and.
     .             isBoundaryCell(cflags(rix,riy,CELLFLAG_TYPE)) .and.
     .             points_match(crx(rix,riy,2), cry(rix,riy,2),
     .                          crx(rix,riy,3), cry(rix,riy,3))) .or.
     .            (isGhostCell(cflags(rix,riy,CELLFLAG_TYPE)) .and.
     .             isBoundaryCell(cflags(lix,liy,CELLFLAG_TYPE)) .and.
     .             points_match(crx(lix,liy,2), cry(lix,liy,2),
     .                          crx(lix,liy,3), cry(lix,liy,3)))) then
              elseif (isGhostCell(cflags(rix,riy,CELLFLAG_TYPE))) then
                ! adapt local connectivity arrays
                loctopix(ix,iy) = locrightix(topix(lix,liy),
     .                                       topiy(lix,liy))
                loctopiy(ix,iy) = locrightiy(topix(lix,liy),
     .                                       topiy(lix,liy))
                loctopix(rix,riy) = locrightix(topix(lix,liy),
     .                                         topiy(lix,liy))
                loctopiy(rix,riy) = locrightiy(topix(lix,liy),
     .                                         topiy(lix,liy))
              elseif (isGhostCell(cflags(lix,liy,CELLFLAG_TYPE))) then
                ! adapt local connectivity arrays
                loctopix(ix,iy) = topix(rix,riy)
                loctopiy(ix,iy) = topiy(rix,riy)
              endif
            endif
          endif
        enddo
      enddo

! Final pass: adapt local connectivity arrays of guard cells
      do iy = -1, ny
        do ix = -1, nx
          if (isGhostCell(cflags(ix,iy,CELLFLAG_TYPE))) then
            if(icorner(ix,iy) == 0) then
              bix = bottomix(ix,iy)
              biy = bottomiy(ix,iy)
              if (biy/=-2) then
                if (isBoundaryCell(cflags(bix,biy,CELLFLAG_TYPE))) then
                  ! N guard cell
                  if (locrightix(ix,iy)==nx+1) then
                    locrightix(ix,iy) = loctopix(rightix(bix,biy),
     .                                           rightiy(bix,biy))
                    locrightiy(ix,iy) = loctopiy(rightix(bix,biy),
     .                                           rightiy(bix,biy))
                  endif
                endif
              endif
              lix = leftix(ix,iy)
              liy = leftiy(ix,iy)
              if (lix/=-2) then
                if (isBoundaryCell(cflags(lix,liy,CELLFLAG_TYPE))) then
                  ! For E guard cells
                  if (loctopiy(ix,iy)==ny+1) then
                    loctopix(ix,iy) = locrightix(topix(lix,liy),
     .                                           topiy(lix,liy))
                    loctopiy(ix,iy) = locrightiy(topix(lix,liy),
     .                                           topiy(lix,liy))
                  endif
                endif
              endif
            else
            endif
          endif
        enddo
      enddo


! Sanity check
      cellError=.false.
      do iy = -1, ny
        do ix = -1, nx
          if (isRealCell(cflags(ix,iy,CELLFLAG_TYPE))) then
            if ((locrightix(ix,iy)==nx+1) .or.
     .          (locrightiy(ix,iy)==ny+1)) then
              write (*,*) 'Real cell ix, iy = ', ix, iy
              write (*,*) 'locrightix, locrightiy = ',
     .          locrightix(ix,iy),
     .          locrightiy(ix,iy)
              cellError=.true.
!              call xerrab ('Problem in definition locrightix.')
            endif
            if ((loctopix(ix,iy)==nx+1) .or.
     .          (loctopiy(ix,iy)==ny+1)) then
              write (*,*) 'Real cell ix, iy = ', ix, iy
              write (*,*) 'loctopix, loctopiy = ',
     .          loctopix(ix,iy),
     .          loctopiy(ix,iy)
              cellError=.true.
!              call xerrab ('Problem in definition loctopiy.')
            endif
          elseif (isGhostCell(cflags(ix,iy,CELLFLAG_TYPE))) then
            if(icorner(ix,iy) == 0) then
              bix = bottomix(ix,iy)
              biy = bottomiy(ix,iy)
              if (biy/=-2) then
                if (isBoundaryCell(cflags(bix,biy,CELLFLAG_TYPE))) then
                  ! N guard cell
                  if ((locrightix(ix,iy)==nx+1) .or.
     .                (locrightiy(ix,iy)==ny+1)) then
                   write (*,*) 'Ghost cell ix, iy = ', ix, iy
                   write (*,*) 'locrightix, locrightiy = ',
     .               locrightix(ix,iy),
     .               locrightiy(ix,iy)
                   cellError=.true.
!                    call xerrab ('Problem in definition locrightix.')
                  endif
                endif
              endif
              lix = leftix(ix,iy)
              liy = leftiy(ix,iy)
              if (lix/=-2) then
                if (isBoundaryCell(cflags(lix,liy,CELLFLAG_TYPE))) then
                  ! For E guard cells
                  if ((loctopix(ix,iy)==nx+1) .or.
     .                (loctopiy(ix,iy)==ny+1)) then
                   write (*,*) 'Ghost cell ix, iy = ', ix, iy
                   write (*,*) 'loctopix, loctopiy = ',
     .               loctopix(ix,iy),
     .               loctopiy(ix,iy)
                   cellError=.true.
!                    call xerrab ('Problem in definition loctopiy.')
                  endif
                endif
              endif
            endif
          endif
        enddo
      enddo

      if (cellError)
     .  call xerrab ('Problem in definition local connectivity.')

      m%nVx = (nx+2)*(ny+2) - nghostvxi

! no. of flux surfaces (conservative estimate)
      m%nFs = max(1,(1 + 1*nncut)*(ny + 2))
      m%nFt = m%nFs

! max. no. of entries in cvVx (vertex list of each cell)
      m%nCmxVx = 4 * m%nCi + 2*m%nCg

! max. no. of entries in cvFc (face list of each cell)
      m%nCmxFc = 4 * m%nCi + m%nCg

! max. no. of entries in vxCv (cell list of each vertex)
! corner vertices have only 3;
! each X-point vertex has 8; duplicate X-point vertices for slabs
!      m%nVmxCv = 4 * m%nVx - 4*nncut + 2*nncut*4
      m%nVmxCv = 4 * (nx+1) * (ny+1) + 2*nncut*4

! max. no. of entries is vxFc (face list of each vertex)
! corner vertices have only 2;
! other vertices on boundaries have only 3;
! each X-point vertex has 8; duplicate X-point vertices for slabs
!      m%nVmxFc = 4 * m%nVx - 4*2*nncut - 2*((nx-1) + (ny-1)) + 2*nncut*4
      m%nVmxFc = 4 * (nx+1) * (ny+1) - 2*((nx-1) + (ny-1))

! max no. of entries in matrix
!      m%nCmxNv = 9*m%nCi - 4*nncut + ! internal cells (corners only 8)
!     .           6*m%nCg - 8*nncut + ! guard cells (at corners only 5)
!     .           2*8*nncut*4         ! correction for cells around X-points
!      m%nCmxNv = 9 * (nx+1) * (ny+1)
      m%nCmxNv = 9 * (nx+2) * (ny+2)  ! there are 2 guard cells
      write (*,*) 'Mapping defined, with'
      write (*,*) 'nCv, nCi, nCg  = ', m%nCv, m%nCi, m%nCg
      write (*,*) 'nFc, nVx       = ', m%nFc, m%nVx
      write (*,*) 'nFs, nFt       = ', m%nFs, m%nFt
      write (*,*) 'nCmxFc, nCmxVx = ', m%nCmxFc, m%nCmxVx
      write (*,*) 'nVmxFc, nVmxCv = ', m%nVmxFc, m%nVmxCv
      write (*,*) 'nCmxNv         = ', m%nCmxNv

      call alloc_mapping (m)

! Compose the mapping matrix from structured to unstructured index
! for Cv's
      m%imapCv= 0 ! 0 means not in use in unstructured
      ici =0
      icg = m%nCi
      do iy = 0, ny+1
        ixloop: do ix = 0, nx+1
          if(isRealCell(cflags(ix-1,iy-1,CELLFLAG_TYPE))) then
            ici = ici + 1
            m%imapCv(ix-1,iy-1) = ici
          elseif(isGhostCell(cflags(ix-1,iy-1,CELLFLAG_TYPE))) then
            if (icorner(ix-1,iy-1)/=0) cycle ixloop
            icg = icg+1
            m%imapCv(ix-1,iy-1) = icg
          endif
        enddo ixloop
      enddo

      if (ici /= m%nCi) then
        call xerrab
     .    ('Inconsistency between counting and mapping in nCi!')
      endif
      if (icg /= m%nCv) then
        call xerrab
     .    ('Inconsistency between counting and mapping in nCv!')
      endif

! for Fcx's and Fcy's
      m%imapFcx= 0 ! 0 means not in use in unstructured
      ici = 0
      do iy = 0, ny
        do ix = 0, nx
          if(ighostfcx(ix,iy)==0) then
            ici = ici + 1
            m%imapFcx(ix,iy) = ici
           endif
        enddo
      enddo

      m%imapFcy= 0 ! 0 means not in use in unstructured
      do ix = 0, nx
        do iy = 0, ny
          if(ighostfcy(ix,iy)==0) then
            ici = ici + 1
            m%imapFcy(ix,iy) = ici
           endif
        enddo
      enddo

      if (ici /= m%nFc) then
        call xerrab
     .  ('Inconsistency between counting and mapping in nFc!')
      endif

      ier = 0
      do iy = -1, ny
        do ix = -1, nx
          if (any(cflags(ix,iy,CELLFLAG_LEFTFACE:CELLFLAG_TOPFACE)/=0))
     .      then
! leftface
            if (cflags(ix,iy,CELLFLAG_LEFTFACE) /= 0) then
              lbl = cflags(ix,iy,CELLFLAG_LEFTFACE)
              iFc = m%imapFcx(ix,iy)
              if (iFc > 0) then
                if (m%fcLbl(iFc) == 0) then
                  m%fcLbl(iFc) = lbl
                else if (m%fcLbl(iFc) /= lbl) then
                  write (*,'(a,5i6)')
     .                 'LEFTFACE: ix,iy,iFc,fcLbl,cflags ',
     .                 ix,iy,iFc,m%fcLbl(iFc),lbl
                  ier = ier + 1
                end if
              end if
            end if
! rightface
            if (cflags(ix,iy,CELLFLAG_RIGHTFACE) /= 0) then
              lbl = cflags(ix,iy,CELLFLAG_RIGHTFACE)
              if (rightix(ix,iy) < nx+1) then
                iFc = m%imapFcx(rightix(ix,iy),rightiy(ix,iy))
              else
                iFc = 0
              end if
              if (iFc > 0) then
                if (m%fcLbl(iFc) == 0) then
                  m%fcLbl(iFc) = lbl
                else if (m%fcLbl(iFc) /= lbl) then
                  write (*,'(a,5i6)')
     .                 'RIGHTFACE: ix,iy,iFc,fcLbl,cflags ',
     .                 ix,iy,iFc,m%fcLbl(iFc),lbl
                  ier = ier + 1
                end if
              end if
            end if
! bottomface
            if (cflags(ix,iy,CELLFLAG_BOTTOMFACE) /= 0) then
              lbl = cflags(ix,iy,CELLFLAG_BOTTOMFACE)
              iFc = m%imapFcy(ix,iy)
              if (iFc > 0) then
                if (m%fcLbl(iFc) == 0) then
                  m%fcLbl(iFc) = lbl
                else if (m%fcLbl(iFc) /= lbl) then
                  write (*,'(a,5i6)')
     .                 'BOTTOMFACE: ix,iy,iFc,fcLbl,cflags ',
     .                 ix,iy,iFc,m%fcLbl(iFc),lbl
                  ier = ier + 1
                end if
              end if
            end if
! topface
            if (cflags(ix,iy,CELLFLAG_TOPFACE) /= 0) then
              lbl = cflags(ix,iy,CELLFLAG_TOPFACE)
              if (topiy(ix,iy) < ny+1) then
                iFc = m%imapFcy(topix(ix,iy),topiy(ix,iy))
              else
                iFc = 0
              end if
              if (iFc > 0) then
                if (m%fcLbl(iFc) == 0) then
                  m%fcLbl(iFc) = lbl
                else if (m%fcLbl(iFc) /= lbl) then
                  write (*,'(a,5i6)')
     .                 'TOPFACE: ix,iy,iFc,fcLbl,cflags ',
     .                 ix,iy,iFc,m%fcLbl(iFc),lbl
                  ier = ier + 1
                end if
              end if
            end if
          end if
        end do
      end do

      write (*,*) 'face labels set ',count(m%fcLbl /= 0)

      write (*,*) 'face labels '
      do iFc = 1, m%nFc
        if (m%fcLbl(iFc) /= 0) write (*,*) iFc, m%fcLbl(iFc)
      end do

      if (ier > 0)
     .  call xerrab ('calc_mapping -- problem with fcLbl')

! for Vx
      m%imapVx = 0 ! 0 means not in use in unstructured
      ici = 0
      do iy = -1, ny
        do ix = -1, nx
          if(ighostvx(ix,iy)==0) then
            ici = ici + 1
            m%imapVx(ix,iy) = ici
           endif
        enddo
      enddo

      if (ici /= m%nVx) then
        call xerrab
     .    ('Inconsistency between counting and mapping in nVx!')
      endif

!     assign same vertex number to coinciding vertices
      write (*,*) 'Merging vertices'
      do iy = -1, ny
        do ix = -1, nx
          ix2 = ix
          iy2 = iy
          do while (ighostvx(ix2,iy2)==1)
            ix2o = ix2
            iy2o = iy2
            ix2= ixmapghostvx(ix2o,iy2)
            iy2= iymapghostvx(ix2o,iy2)
            if ((ix2 == -2) .or. (ix2 == nx+1) .or.
     .          (iy2 == -2) .or. (iy2 == ny+1)) exit
          enddo
          m%imapVx(ix,iy)= m%imapVx(min(max(ix2,-1),nx),
     .                              min(max(iy2,-1),ny))
        enddo
      enddo

! consistency check
      do iy = -1, ny
        do ix = -1, nx
          if(ighostvx(ix,iy) == 0 .and. m%imapVx(ix,iy) == 0) then
            write (*,*) 'ix, iy = ', ix, iy
            call xerrab ('Problem with vertex mapping')
           endif
        enddo
      enddo

      write (*,*) 'Merged vertices'

! store indices of corner coordinates
      m%icornVx = 0
      if (ncorneri > 0) then
!  sort corner points per ix position first
        DO j = ncorneri-1, 1, -1
          swapped = .FALSE.
          DO i = 1, j
            IF (ixcorner(i) > ixcorner(i+1)) THEN
              itmpx = ixcorner(i)
              itmpy = iycorner(i)
              ixcorner(i) = ixcorner(i+1)
              iycorner(i) = iycorner(i+1)
              ixcorner(i+1) = itmpx
              iycorner(i+1) = itmpy
              swapped = .TRUE.
            END IF
          END DO
          IF (.NOT. swapped) EXIT
        END DO

!  now sort corner points in iy position
        do j=ncorneri-1, 1, -1
          swapped = .FALSE.
          do i = 1, j
            if (ixcorner(i) /= ixcorner(i+1)) cycle
            if (iycorner(i) > iycorner(i+1)) then
              itmpx = ixcorner(i)
              itmpy = iycorner(i)
              ixcorner(i) = ixcorner(i+1)
              iycorner(i) = iycorner(i+1)
              ixcorner(i+1) = itmpx
              iycorner(i+1) = itmpy
              swapped = .TRUE.
            end if
          end do
          IF (.NOT. swapped) EXIT
        end do

        icount = 0
        do ic = 1, ncorneri
          icount = icount + 1
          ix = ixcorner(ic)
          if (count(m%ictcntrl(1:m%nctcnt) ==ix) > 0) ix = ix + 1
          ix = min(max(ix,0),nx)
          iy = iycorner(ic)
          iy = min(max(iy,0),ny)
          m%icornVx(icount) = m%imapVx(ix,iy)
        end do
      end if

#ifdef WG_DEBUG
      imapCv  = m%imapCv
      imapFcx = m%imapFcx
      imapFcy = m%imapFcy
      imapVx  = m%imapVx
#endif

! no. of faces and vertices to be considered for this cell; only consider all internal and boundary cells are either quadrilateral or triangle
      m%cvFcP(m%nCi+1:m%nCv,2) = 1
      m%cvVxP(m%nCi+1:m%nCv,2) = 2

      do iy = -1, ny
        do ix = -1, nx
          if ((m%imapCv(ix,iy) /= 0).and.
     .        (m%imapCv(ix,iy) .le. m%nCi)) then
            xco = crx(ix,iy,0:3)
            yco = cry(ix,iy,0:3)
            if(cellGeoType(xco,yco) == CGEO_QUAD) then
              m%cvFcP(m%imapCv(ix,iy),2) = 4
              m%cvVxP(m%imapCv(ix,iy),2) = 4
            else
              m%cvFcP(m%imapCv(ix,iy),2) = 3
              m%cvVxP(m%imapCv(ix,iy),2) = 3
            endif
          endif
        enddo
      enddo

! starting index in the list of cell faces
! starting index for every cell in the cvVx list and its number of vertices
      m%cvFcP(1,1) = 1
      m%cvVxP(1,1) = 1
      do i=2, m%nCv
        m%cvFcP(i,1) = m%cvFcP(i-1,1) + m%cvFcP(i-1,2)
        m%cvVxP(i,1) = m%cvVxP(i-1,1) + m%cvVxP(i-1,2)
      end do

! list of cell faces and vertices for cells
! interpolation coefficients for interpolation to cell centers
      do iy = 0, ny+1
        do ix = 0, nx+1
          if ((m%imapCv(ix-1,iy-1) /= 0).and.
     .        (m%imapCv(ix-1,iy-1) .le. m%nCi)) then
            ic = m%imapCv(ix-1,iy-1)
            xco = crx(ix-1,iy-1,0:3)
            yco = cry(ix-1,iy-1,0:3)
            select case
     .        ( cellGeoType(xco,yco) )
            case ( CGEO_QUAD )
              if1 = m%cvFcP(ic,1)
              if2 = if1 + 1
              if3 = if2 + 1
              if4 = if3 + 1
              m%cvFc(if1) = m%imapFcx(ix-1,iy-1)
              m%cvFc(if2) = m%imapFcx(rightix(ix-1,iy-1),
     .                                rightiy(ix-1,iy-1))
              m%cvFc(if3) = m%imapFcy(ix-1,iy-1)
              m%cvFc(if4) = m%imapFcy(topix(ix-1,iy-1),
     .                                topiy(ix-1,iy-1))
              m%intcellP(if1) = 1.0_R8
              m%intcellP(if2) = 1.0_R8
              m%intcellP(if3) = 0.0_R8
              m%intcellP(if4) = 0.0_R8
              m%intcellR(if1) = 0.0_R8
              m%intcellR(if2) = 0.0_R8
              m%intcellR(if3) = 1.0_R8
              m%intcellR(if4) = 1.0_R8
              m%cvVx(if1) = m%imapVx(ix-1,iy-1)
              m%cvVx(if2) = m%imapVx(rightix(ix-1,iy-1),
     .                               rightiy(ix-1,iy-1))
              m%cvVx(if3) = m%imapVx(topix(ix-1,iy-1),
     .                               topiy(ix-1,iy-1))
              m%cvVx(if4) = m%imapVx(locrightix(loctopix(ix-1,iy-1),
     .                                          loctopiy(ix-1,iy-1)),
     .                               locrightiy(loctopix(ix-1,iy-1),
     .                                          loctopiy(ix-1,iy-1)))
           case ( CGEO_TRIA_NOLEFT )
              if1 = m%cvFcP(ic,1)
              if2 = if1 + 1
              if3 = if2 + 1
              m%cvFc(if1) = m%imapFcx(rightix(ix-1,iy-1),
     .                                rightiy(ix-1,iy-1))
              m%cvFc(if2) = m%imapFcy(ix-1,iy-1)
              m%cvFc(if3) = m%imapFcy(topix(ix-1,iy-1),
     .                                topiy(ix-1,iy-1))
              m%intcellP(if1) = 1.0_R8
              m%intcellP(if2) =
     .          abs(fpsi(ix-1,iy-1,2,0)-fpsi(ix-1,iy-1,1,0))/
     .          abs(fpsi(ix-1,iy-1,3,0)-fpsi(ix-1,iy-1,1,0))
              m%intcellP(if3) =
     .          abs(fpsi(ix-1,iy-1,2,0)-fpsi(ix-1,iy-1,3,0))/
     .          abs(fpsi(ix-1,iy-1,3,0)-fpsi(ix-1,iy-1,1,0))
              m%intcellR(if1) = 1.0_R8
              m%intcellR(if2) =
     .          abs(fpsi(ix-1,iy-1,2,0)-fpsi(ix-1,iy-1,1,0))/
     .          abs(fpsi(ix-1,iy-1,3,0)-fpsi(ix-1,iy-1,1,0))
              m%intcellR(if3) =
     .          abs(fpsi(ix-1,iy-1,2,0)-fpsi(ix-1,iy-1,3,0))/
     .          abs(fpsi(ix-1,iy-1,3,0)-fpsi(ix-1,iy-1,1,0))
              m%cvVx(if1) = m%imapVx(ix-1,iy-1)
              m%cvVx(if2) = m%imapVx(rightix(ix-1,iy-1),
     .                               rightiy(ix-1,iy-1))
              m%cvVx(if3) = m%imapVx(locrightix(loctopix(ix-1,iy-1),
     .                                          loctopiy(ix-1,iy-1)),
     .                               locrightiy(loctopix(ix-1,iy-1),
     .                                          loctopiy(ix-1,iy-1)))
            case ( CGEO_TRIA_NOBOT )
              if1 = m%cvFcP(ic,1)
              if2 = if1 + 1
              if3 = if2 + 1
              m%cvFc(if1) = m%imapFcx(ix-1,iy-1)
              m%cvFc(if2) = m%imapFcx(rightix(ix-1,iy-1),
     .                                rightiy(ix-1,iy-1))
              m%cvFc(if3) = m%imapFcy(topix(ix-1,iy-1),
     .                                topiy(ix-1,iy-1))
              m%intcellP(if1) = 1.0_R8
              m%intcellP(if2) = 1.0_R8
              m%intcellP(if3) = 0.0_R8
              m%intcellR(if1) = 1.0_R8
              m%intcellR(if2) = 1.0_R8
              m%intcellR(if3) = 0.0_R8
              m%cvVx(if1) = m%imapVx(ix-1,iy-1)
              m%cvVx(if2) = m%imapVx(topix(ix-1,iy-1),
     .                               topiy(ix-1,iy-1))
              m%cvVx(if3) = m%imapVx(locrightix(loctopix(ix-1,iy-1),
     .                                          loctopiy(ix-1,iy-1)),
     .                               locrightiy(loctopix(ix-1,iy-1),
     .                                          loctopiy(ix-1,iy-1)))
            case ( CGEO_TRIA_NORIGHT )
              if1 = m%cvFcP(ic,1)
              if2 = if1 + 1
              if3 = if2 + 1
              m%cvFc(if1) = m%imapFcx(ix-1,iy-1)
              m%cvFc(if2) = m%imapFcy(ix-1,iy-1)
              m%cvFc(if3) = m%imapFcy(topix(ix-1,iy-1),
     .                                topiy(ix-1,iy-1))
              m%intcellP(if1) = 1.0_R8
              m%intcellP(if2) =
     .          abs(fpsi(ix-1,iy-1,0,0)-fpsi(ix-1,iy-1,3,0))/
     .          abs(fpsi(ix-1,iy-1,2,0)-fpsi(ix-1,iy-1,0,0))
              m%intcellP(if3) =
     .          abs(fpsi(ix-1,iy-1,2,0)-fpsi(ix-1,iy-1,3,0))/
     .          abs(fpsi(ix-1,iy-1,2,0)-fpsi(ix-1,iy-1,0,0))
              m%intcellR(if1) = 1.0_R8
              m%intcellR(if2) =
     .          abs(fpsi(ix-1,iy-1,0,0)-fpsi(ix-1,iy-1,3,0))/
     .          abs(fpsi(ix-1,iy-1,2,0)-fpsi(ix-1,iy-1,0,0))
              m%intcellR(if3) =
     .          abs(fpsi(ix-1,iy-1,2,0)-fpsi(ix-1,iy-1,3,0))/
     .          abs(fpsi(ix-1,iy-1,2,0)-fpsi(ix-1,iy-1,0,0))
              m%cvVx(if1) = m%imapVx(ix-1,iy-1)
              m%cvVx(if2) = m%imapVx(locrightix(ix-1,iy-1),!cwdk is this not the vertex that was thrown out??ctb solved?
     .                               locrightiy(ix-1,iy-1))
              m%cvVx(if3) = m%imapVx(loctopix(ix-1,iy-1),
     .                               loctopiy(ix-1,iy-1))
            case ( CGEO_TRIA_NOTOP )
              if1 = m%cvFcP(ic,1)
              if2 = if1 + 1
              if3 = if2 + 1
              m%cvFc(if1) = m%imapFcx(ix-1,iy-1)
              m%cvFc(if2) = m%imapFcx(rightix(ix-1,iy-1),
     .                                rightiy(ix-1,iy-1))
              m%cvFc(if3) = m%imapFcy(ix-1,iy-1)
              m%intcellP(if1) = 1.0_R8
              m%intcellP(if2) = 1.0_R8
              m%intcellP(if3) = 0.0_R8
              m%intcellR(if1) = 1.0_R8
              m%intcellR(if2) = 1.0_R8
              m%intcellR(if3) = 0.0_R8
              m%cvVx(if1) = m%imapVx(ix-1,iy-1)
              m%cvVx(if2) = m%imapVx(locrightix(ix-1,iy-1),
     .                               locrightiy(ix-1,iy-1))
              m%cvVx(if3) = m%imapVx(loctopix(ix-1,iy-1),
     .                               loctopiy(ix-1,iy-1))!cwdk is this not the vertex that was thrown out?? ctb solved?
            end select
          elseif ((m%imapCv(ix-1,iy-1) /= 0).and.
     .            (m%imapCv(ix-1,iy-1) .gt. m%nCi)) then
            bound = .false.

!tb  the vertex numbering might not be consistent for fcVx, to be checked
            ! N guard cell
            ixb = bottomix(ix-1,iy-1)
            iyb = bottomiy(ix-1,iy-1)
            if (iyb /= -2) then
              if (isBoundaryCell(cflags(ixb,iyb,CELLFLAG_TYPE))) then
                ic = m%imapCv(ix-1,iy-1)
                if1 = m%cvFcP(ic,1)
                m%cvFc(if1) = m%imapFcy(ix-1,iy-1)
                m%intcellP(if1) = 1.0_R8
                m%intcellR(if1) = 1.0_R8
                if1 = m%cvVxP(ic,1)
                m%cvVx(if1)   = m%imapVx(ix-1,iy-1)
                m%cvVx(if1+1) = m%imapVx(locrightix(ix-1,iy-1),
     .                                   locrightiy(ix-1,iy-1))
                bound = .true.
              end if
            end if

            ! S guard cell
            ixb = topix(ix-1,iy-1)
            iyb = topiy(ix-1,iy-1)
            if (iyb /= ny+1) then
              if (isBoundaryCell(cflags(ixb,iyb,CELLFLAG_TYPE))) then
                if (bound) call xerrab('calc_mapping - double boundary')
                ic = m%imapCv(ix-1,iy-1)
                if1 = m%cvFcP(ic,1)
                m%cvFc(if1) = m%imapFcy(ixb,iyb)
                m%intcellP(if1) = 1.0_R8
                m%intcellR(if1) = 1.0_R8
                if1 = m%cvVxP(ic,1)
                m%cvVx(if1)   = m%imapVx(ixb,iyb)
                m%cvVx(if1+1) = m%imapVx(locrightix(ixb,iyb),
     .                                   locrightiy(ixb,iyb))
                bound = .true.
              end if
            end if

            ! E guard cell
            ixb = leftix(ix-1,iy-1)
            iyb = leftiy(ix-1,iy-1)
            if (ixb /= -2) then
              if (isBoundaryCell(cflags(ixb,iyb,CELLFLAG_TYPE))) then
                if (bound) call xerrab('calc_mapping - double boundary')
                ic = m%imapCv(ix-1,iy-1)
                if1 = m%cvFcP(ic,1)
                m%cvFc(if1) = m%imapFcx(ix-1,iy-1)
                m%intcellP(if1) = 1.0_R8
                m%intcellR(if1) = 1.0_R8
                if1 = m%cvVxP(ic,1)
                m%cvVx(if1)   = m%imapVx(ix-1,iy-1)
                m%cvVx(if1+1) = m%imapVx(loctopix(ix-1,iy-1),
     .                                   loctopiy(ix-1,iy-1))
                bound = .true.
              end if
            end if

            ! W guard cell
            ixb = rightix(ix-1,iy-1)
            iyb = rightiy(ix-1,iy-1)
            if (ixb /= nx+1) then
              if (isBoundaryCell(cflags(ixb,iyb,CELLFLAG_TYPE))) then
                if (bound) call xerrab('calc_mapping - double boundary')
                ic = m%imapCv(ix-1,iy-1)
                if1 = m%cvFcP(ic,1)
                m%cvFc(if1) = m%imapFcx(ixb,iyb)
                m%intcellP(if1) = 1.0_R8
                m%intcellR(if1) = 1.0_R8
                if1 = m%cvVxP(ic,1)
                m%cvVx(if1)   = m%imapVx(ixb,iyb)
                m%cvVx(if1+1) = m%imapVx(loctopix(ixb,iyb),
     .                                   loctopiy(ixb,iyb))
                bound = .true.
              end if
            end if

            if (.not.bound) then
              call xerrab ('Problem with guard cell definition.')
            endif
          endif
        end do
      end do


! listing for every face the two corresponding cells
! and the corresponding vertices
! x faces
       do iy = 1, ny
         do ix = 1, nx+1
           if(m%imapFcx(ix-1,iy-1)/=0) then
             ifc =m%imapFcx(ix-1,iy-1)
             m%fcCv(ifc,1) = m%imapCv(leftix(ix-1,iy-1),
     .                                leftiy(ix-1,iy-1))
             m%fcCv(ifc,2) = m%imapCv(ix-1,iy-1)
             m%fcVx(ifc,1) = m%imapVx(ix-1,iy-1)
             m%fcVx(ifc,2) = m%imapVx(loctopix(ix-1,iy-1),
     .                                loctopiy(ix-1,iy-1))
           endif
         end do
       end do


! y faces

       do ix = 1, nx
         do iy = 1, ny+1
           if(m%imapFcy(ix-1,iy-1)/=0) then
             ifc =m%imapFcy(ix-1,iy-1)
             m%fcCv(ifc,2) = m%imapCv(ix-1,iy-1)
             m%fcCv(ifc,1) = m%imapCv(bottomix(ix-1,iy-1),
     .                                bottomiy(ix-1,iy-1))
             m%fcVx(ifc,2) = m%imapVx(ix-1,iy-1)
             m%fcVx(ifc,1) = m%imapVx(locrightix(ix-1,iy-1),
     .                                locrightiy(ix-1,iy-1))
           endif
         end do
       end do

! merging x-point vertex !!! TB not changed... seems not to be necessary
       do icut = 1, nncut
         ix1 = leftcut(icut)
         iy1 = topcut(icut)
         ix2 = rightcut(icut)
         iy2 = topcut(icut)
         ivx1 = m%imapVx(ix1,iy1)
         ivx2 = m%imapVx(ix2,iy2)

         !wdk Check whether X-point vertices actually match
         !    (One or both of the vertices might fall outside of domain)
         if (points_match(crx(ix1,iy1,0), cry(ix1,iy1,0),
     .                    crx(ix2,iy2,0), cry(ix2,iy2,0))) then
           where (m%imapVx == ivx2)
             m%imapVx = ivx1
           end where

           where (m%cvVx == ivx2)
             m%cvVx = ivx1
           end where

           where (m%fcVx == ivx2)
             m%fcVx = ivx1
           end where

           write (*,*) 'X-point ', icut , ' merged.'
         else
           write (*,*) 'X-point ', icut , ' not merged!'
         endif
       end do

      ! sanity check
      do ic = 1,m%nCv
        do ifc = m%cvFcP(ic,1),m%cvFcP(ic,1)+m%cvFcP(ic,2)-1
          if(m%cvFc(ifc) == 0) then
            write (*,*) 'ic, ifc = ', ic, ifc
            call xerrab ('Problem with mapping: undefined face')
           endif
        enddo
        do ivx = m%cvVxP(ic,1),m%cvVxP(ic,1)+m%cvVxP(ic,2)-1
          if(m%cvVx(ivx) == 0) then
            write (*,*) 'ic, ivx = ', ic, ivx
            call xerrab ('Problem with mapping: undefined vertex')
           endif
        enddo
      enddo

      do ifc = 1,m%nFc
        if(m%fcCv(ifc,1) == 0) then
          write (*,*) 'ifc = ', ifc
          call xerrab ('Problem with mapping: undefined cv1')
        endif
        if(m%fcCv(ifc,2) == 0) then
          write (*,*) 'ifc = ', ifc
          call xerrab ('Problem with mapping: undefined cv2')
        endif
        if(m%fcVx(ifc,1) == 0) then
          write (*,*) 'ifc = ', ifc
          call xerrab ('Problem with mapping: undefined vx1')
        endif
        if(m%fcVx(ifc,2) == 0) then
          write (*,*) 'ifc = ', ifc
          call xerrab ('Problem with mapping: undefined vx2')
        endif
      enddo


! listing for every vertex the cell and face numbers
! and the starting points in the list
       allocate (indFc(2*m%nFc))
       indFc(1:m%nFc)         = (/ (i, i=1,m%nFc) /)
       indFc(m%nFc+1:2*m%nFc) = (/ (i, i=1,m%nFc) /)

       allocate (fh(2*m%nFc))
       fh(1:m%nFc)         = m%fcVx(1:m%nFc,1)
       fh(m%nFc+1:2*m%nFc) = m%fcVx(1:m%nFc,2)


       allocate (indCv(m%nCmxVx))
       do ic = 1, m%nCv
          indCv(m%cvVxP(ic,1):m%cvVxP(ic,1)+m%cvVxP(ic,2)-1) = ic
       enddo

       allocate (ch(m%nCmxVx))
       ch(1:m%nCmxVx) = m%cvVx(1:m%nCmxVx)

       m%vxFcP(1,1) = 1
       m%vxCvP(1,1) = 1
       do ivx = 1, m%nVx-1
         m%vxFcP(ivx,2) = count(fh == ivx)
         m%vxFcP(ivx+1,1) = m%vxFcP(ivx,1) + m%vxFcP(ivx,2)
         if (m%vxFcP(ivx,2) > 0) then
           m%vxFc(m%vxFcP(ivx,1):m%vxFcP(ivx+1,1)-1) =
     .          pack( indFc, fh == ivx )
         end if
         m%vxCvP(ivx,2) = count(ch == ivx)
         m%vxCvP(ivx+1,1) = m%vxCvP(ivx,1) + m%vxCvP(ivx,2)
         if (m%vxCvP(ivx,2) > 0) then
           m%vxCv(m%vxCvP(ivx,1):m%vxCvP(ivx+1,1)-1) =
     .          pack( indCv, ch == ivx )
         end if
       end do

       ivx = m%nVx
       m%vxFcP(ivx,2) = count(fh == ivx)
       vxFc_end = m%vxFcP(ivx,1) + m%vxFcP(ivx,2) -1
       if (m%vxFcP(ivx,2) > 0) then
         m%vxFc(m%vxFcP(ivx,1):vxFc_end) =
     .          pack( indFc, fh == ivx )
       end if
       m%vxCvP(ivx,2) = count(ch == ivx)
       vxCv_end = m%vxCvP(ivx,1) + m%vxCvP(ivx,2) - 1
       if (m%vxCvP(ivx,2) > 0) then
         m%vxCv(m%vxCvP(ivx,1):vxCv_end) =
     .          pack( indCv, ch == ivx )
       end if

       deallocate (indFc)
       deallocate (indCv)
       deallocate (fh)
       deallocate (ch)

! regions (cell volumes)
      do iy = -1, ny
        do ix = -1, nx
          if(m%imapCv(ix,iy)/=0) then
            ic = m%imapCv(ix,iy)
            m%cvReg(ic) = region(ix,iy,0)
            m%cvOnClosedSurface(ic) = on_closed_surface(ix,iy)
          endif
        end do
      end do

! regions (faces)
! x faces
       nfcxReg = maxval(region(-1:nx,-1:ny,1))
       do iy = 0, ny-1
         do ix = 0, nx
           if(m%imapFcx(ix,iy)/=0) then
             ifc = m%imapFcx(ix,iy)
             if (region(ix,iy,1) .gt. 0) then
               m%fcReg(ifc) = region(ix,iy,1)
             elseif (region(ix,iy,1) .lt. 0) then
               ! actually a y-face region!
               m%fcReg(ifc) = abs(region(ix,iy,1)) + nfcxReg
             endif
           endif
         end do
       end do

! y faces
       do ix = 0, nx-1
         do iy = 0, ny
           if(m%imapFcy(ix,iy)/=0) then
             ifc = m%imapFcy(ix,iy)
             if (region(ix,iy,2) .gt. 0) then
               m%fcReg(ifc) = region(ix,iy,2) + nfcxReg
             elseif (region(ix,iy,2) .lt. 0) then
               ! actually a x-face region!
               m%fcReg(ifc) = abs(region(ix,iy,2))
             endif
           end if
         end do
      end do

      mxReg = maxval(m%cvReg(1:m%nCv))

! identify the various flux tubes and surfaces
      m%ftCv = 0 ! init to 0
      m%ftFc = 0 ! init to 0
      m%fsFc = 0 ! init to 0
      m%cvFt = 0 ! init to 0
      m%ftCvP = 0
      m%ftFcP = 0
      m%fsFcP = 0
      m%ftLbl = 0
      if (m%isClassicalGrid.eq.1) then
        call conv_regions(nnreg)
! csc these arrays are needed because otherwise the nFt with which the mapping is written
!     is different wrt the nFt with which they have been allocated
        allocate (ftCvP(m%nFt,2))
        allocate (ftFcP(m%nFt,2))
        allocate (fsFcP(m%nFs,2))
        allocate (ftReg(m%nFt))
        allocate (ftLbl(m%nFt))
        ftCvP = 0
        ftFcP = 0
        fsFcP = 0
        ftReg = 0
        ftLbl = 0
        ift = 0
        ifs = 0
        icv = 0
        ifcc = 0
        ifccy = 0
        ilbl = 0
        if (nncut.eq.0) then
          !Linear case
          nyreg = 1
          allocate(nxreg(nyreg),reg(1,nyreg),nseg(1,nyreg))
          allocate(ixstart(1,nyreg,1),ixend(1,nyreg,1))
          allocate(iystart(nyreg),iyend(nyreg))
          iystart(1) = -1
          iyend(1) = ny
          nxreg(1) = 1
          nseg(1,1) = 1 !
          ixstart(1,1,1) = -1
          ixend  (1,1,1) = nx
          reg(1,1) = 1
        else if (nncut.eq.1) then
          nyreg = 2
          allocate(nxreg(nyreg),reg(2,nyreg),nseg(2,nyreg))
          allocate(ixstart(2,nyreg,2),ixend(2,nyreg,2))
          allocate(iystart(nyreg),iyend(nyreg))
          iystart(1) = -1
          iyend(1) = topcut(1)-1
          iystart(2) = topcut(1)
          iyend(2) = ny
          nxreg(1) = 2
          nxreg(2) = 1
          nseg(1,1) = 2 !segments for pfr
          nseg(2,1) = 1 !segments for core
          nseg(1,2) = 1 !segments for SOL
          ixstart(1,1,1) = -1
          ixend  (1,1,1) = leftcut(1)-1
          ixstart(1,1,2) = rightcut(1)
          ixend  (1,1,2) = nx
          ixstart(2,1,1) = leftcut(1)
          ixend  (2,1,1) = rightcut(1)-1
          ixstart(1,2,1) = -1
          ixend  (1,2,1) = nx
          reg(1,1) = 3
          reg(2,1) = 1
          reg(1,2) = 2
        else
          nyreg = 3
          call get_nxt(nx, nxtl, nxtr)
          allocate(nxreg(nyreg),reg(3,nyreg),nseg(3,nyreg))
          allocate(ixstart(3,nyreg,2),ixend(3,nyreg,2))
          allocate(iystart(nyreg),iyend(nyreg))
          nxreg(1) = 3
          nxreg(2) = 2 !both DDN & DDNU
          nxreg(3) = 2
          iystart(1) = -1
          iyend(1) = min(topcut(1),topcut(2))-1
          iystart(2) = min(topcut(1),topcut(2))
          iyend(2) = max(topcut(1),topcut(2))-1
          iystart(3) = max(topcut(1),topcut(2))
          iyend(3) = ny
          nseg(:,1) = 2 !segments for pfr1/core/pfr2
          nseg(:,2) = 2 !segments for sol1/pfr2
          nseg(:,3) = 1 !segments for sol2/sol3
          ! pfr1
          ixstart(1,1,1) = -1
          ixend  (1,1,1) = leftcut(1)-1
          ixstart(1,1,2) = rightcut(1)
          ixend  (1,1,2) = nx
          reg(1,1) = 3
          ! core
          ixstart(2,1,1) = leftcut(1)
          ixend  (2,1,1) = leftcut(2)-1
          ixstart(2,1,2) = rightcut(2)
          ixend  (2,1,2) = rightcut(1)-1
          reg(2,1) = 1
          ! pfr2
          ixstart(3,1,1) = nxtr -1
          ixend  (3,1,1) = rightcut(2)-1
          ixstart(3,1,2) = leftcut(2)
          ixend  (3,1,2) = nxtl -1
          reg(3,1) = 3
          if (topcut(1).gt.topcut(2)) then !DDNU
            ix = 1
            reg(1,2) = 2
            reg(2,2) = 3
          elseif (topcut(1).lt.topcut(2)) then !DDND
            ix = 2
            reg(1,2) = 3
            reg(2,2) = 2
          else
            ix = 1
            nxreg(2) = 0 !CDN
          endif
          ! sol1(ddnu)/pfr2(ddnd)
          ixstart(1,2,1) = nxtr-1
          ixend  (1,2,1) = rightcut(ix)-1
          ixstart(1,2,2) = leftcut(ix)
          ixend  (1,2,2) = nxtl-1
          ! pfr2(ddnu)/sol1(ddnd)
          ixstart(2,2,1) = -1
          ixend  (2,2,1) = leftcut(ix)-1
          ixstart(2,2,2) = rightcut(ix)
          ixend  (2,2,2) = nx
          ! sol2
          ixstart(1,3,1) = -1
          ixend  (1,3,1) = nxtl-1
          reg(1,3) = 2
          ! sol3
          ixstart(2,3,1) = nxtr-1
          ixend  (2,3,1) = nx
          reg(2,3) = 2
        endif
        do iyreg = 1, nyreg
          do iy = iystart(iyreg), iyend(iyreg)
            ilbl = ilbl + 1
            do ixreg = 1, nxreg(iyreg)
              istrtc = -1
              istrtfc = -1
              istrtfy = -1
              do is = 1, nseg(ixreg,iyreg)
                do ix = ixstart(ixreg,iyreg,is), ixend(ixreg,iyreg,is)
                  if((m%imapCv(ix,iy)/=0)) then
                    ic = m%imapCv(ix,iy)
                    ifc = m%imapFcx(ix,iy)
                    ifcy = m%imapFcy(ix,iy)
                    if (istrtc.eq.-1) then
                      ift = ift + 1
                      icv = icv + 1
                      m%ftCv(icv) = ic
                      ftCvP(ift,1) = icv
                      m%cvFt(ic) = ift
                      ftReg(ift) = reg(ixreg,iyreg)
                      ftLbl(ift) = ilbl
                      istrtc = 1
                    else
                      icv = icv + 1
                      m%ftCv(icv) = ic
                      m%cvFt(ic) = ift
                      ftCvP(ift,2) = icv - ftCvP(ift,1) + 1
                    end if
                    if (ifc.ne.0) then
                      if (istrtfc.eq.-1) then
                        ifcc = ifcc + 1
                        m%ftFc(ifcc) = ifc
                        ftFcP(ift,1) = ifcc
                        istrtfc = 1
                      else
                        ifcc = ifcc + 1
                        m%ftFc(ifcc) = ifc
                        ftFcP(ift,2) = ifcc - ftFcP(ift,1) + 1
                      end if
                    endif
                    if (ifcy.ne.0) then
                      if (istrtfy.eq.-1) then
                        ifs = ifs + 1
                        ifccy = ifccy + 1
                        m%fsFc(ifccy) = ifcy
                        fsFcP(ifs,1) = ifccy
                        istrtfy = 1
                      else
                        ifccy = ifccy + 1
                        m%fsFc(ifccy) = ifcy
                        fsFcP(ifs,2) = ifccy - fsFcP(ifs,1) + 1
                      end if
                    endif
                  end if
                end do !ix
              end do !is
            end do !ixreg
          end do !iy
        end do !iyreg
! csc last flux surface needs ny + 1
        do ixreg = 1, nxreg(nyreg)
          istrtfy = -1
          do is = 1, nseg(ixreg,nyreg)
            do ix = ixstart(ixreg,nyreg,is), ixend(ixreg,nyreg,is)
              if((m%imapCv(ix,ny)/=0)) then
                ifcy = m%imapFcy(ix,ny)
                if (ifcy.ne.0) then
                  if (istrtfy == -1) then
                    ifs = ifs + 1
                    ifccy = ifccy + 1
                    m%fsFc(ifccy) = ifcy
                    fsFcP(ifs,1) = ifccy
                    istrtfy = 1
                  else
                    ifccy = ifccy + 1
                    m%fsFc(ifccy) = ifcy
                    fsFcP(ifs,2) = ifccy - fsFcP(ifs,1) + 1
                  endif
                endif
              endif
            enddo !ix
          enddo !is
        enddo !ixreg
        m%nFt = ift
        m%nFs = ifs
! csc Now that real number of flux surfaces is known re-allocate and assign ftFcP and ftCvP
        deallocate(m%ftFcP,m%ftCvP,m%fsFcP,m%ftReg,m%ftLbl)
        allocate (m%ftFcP(m%nFt,2))
        allocate (m%ftCvP(m%nFt,2))
        allocate (m%fsFcP(m%nFs,2))
        allocate (m%ftReg(m%nFt))
        allocate (m%ftLbl(m%nFt))
        m%ftFcP = ftFcP(1:m%nFt,:)
        m%ftCvP = ftCvP(1:m%nFt,:)
        m%fsFcP = fsFcP(1:m%nFs,:)
        m%ftReg = ftReg(1:m%nFt)
        m%ftLbl = ftLbl(1:m%nFt)
        deallocate(ftFcP,ftCvP,fsFcP,ftReg,ftLbl)
        deallocate(ixstart, ixend, iystart, iyend, reg, nxreg, nseg)
      else
        write (*,*) 'WARNING flux surfaces not yet available'//
     .              ' for non-structured grids'
      endif

!     Create b2.user.parameters file with definition of outer/inner midplanes
      call conv_midplane(nx,ny)

      deallocate (icorner)
      deallocate (ighostcv)
      deallocate (ixghostcv)
      deallocate (iyghostcv)
      deallocate (ixcorner)
      deallocate (iycorner)
      deallocate (ighostfcx)
      deallocate (ixghostfcx)
      deallocate (iyghostfcx)
      deallocate (ighostfcy)
      deallocate (ixghostfcy)
      deallocate (iyghostfcy)
      deallocate (ighostvx)
      deallocate (ixghostvx)
      deallocate (iyghostvx)
!      deallocate (locleftix)
!      deallocate (locleftiy)
      deallocate (locrightix)
      deallocate (locrightiy)
!      deallocate (locbottomix)
!      deallocate (locbottomiy)
      deallocate (loctopix)
      deallocate (loctopiy)

      return
      end subroutine calc_mapping

***********************************************************************

      subroutine prep_geometry(nx, ny, nCv, nCi, nFc, nVx, nFs, nFt,
     .  nVmxCv, gm, mpg)

      implicit none
      integer, intent(in) :: nx, ny, nCv, nCi, nFc, nVx, nFs, nFt,
     .                       nVmxCv
      integer :: ix, iy, ic, ifc, ivx, ift
#ifdef WG_TODO
      integer :: i, j, idb, ivx1, ivx2
#endif
      type (mapping) , intent (in) :: mpg
      type (geometry), intent (out) :: gm
#ifndef F2003
      integer newunit
      external newunit
#endif

      call alloc_geometry (gm, nCv, nFc, nVx, nFs, nFt, nVmxCv)

! cell volumes
      do iy = -1, ny
        do ix = -1, nx
           if(mpg%imapCv(ix,iy)/=0) then
            ic = mpg%imapCv(ix,iy)
! b field
            ! wdk sign/orientation of field now included in angles/unit vector cvEb
            gm%cvBb(ic,0:3) = abs(bb(ix,iy,0:3))
            gm%cvEb(ic,0:2) = ebc(ix,iy,0:2)
! cell center
            gm%cvX(ic) = 0.25_R8*sum(crx(ix,iy,0:3))
            gm%cvY(ic) = 0.25_R8*sum(cry(ix,iy,0:3))
            gm%cvSz(ic) = gs(ix,iy,2)
            gm%cvHz(ic) = hz(ix,iy)
            gm%cvHx(ic) = hx(ix,iy)
            gm%cvQgam(ic,0) = qz(ix,iy,1)    ! cos
            gm%cvQgam(ic,1) = -qz(ix,iy,0)   ! sin
            gm%cvVol(ic) = vol(ix,iy)
            gm%cvOnedBsq(ic) = OnedBsq(ix,iy)
            gm%cvBzb(ic) = bzb(ix,iy) !wdk to check: need abs here?
          end if
        end do
      end do

! x faces
      do iy = 0, ny
        do ix = 0, nx
          if(mpg%imapFcx(ix,iy)/=0) then
            ifc = mpg%imapFcx(ix,iy)
            ! wdk sign of field now included in angles
            gm%fcBb(ifc,0:3) = abs(wbbl(ix,iy,0:3))
            gm%fcS(ifc) = gs(ix,iy,0)
            gm%fcHc(ifc,1) = hc(leftix(ix,iy),leftiy(ix,iy),2)
            gm%fcHc(ifc,2) = hc(ix,iy,0)
            gm%fcHt(ifc)   = ht(ix,iy,0)
            gm%fcVol(ifc,1) = 0.5_R8*vol(leftix(ix,iy),leftiy(ix,iy))
            gm%fcVol(ifc,2) = 0.5_R8*vol(ix,iy)
            gm%fcQgam(ifc,0) = qc(ix,iy,0)
            gm%fcQgam(ifc,1) = -qs(ix,iy,0)
            gm%fcQalf(ifc,0) = qac(ix,iy,0)
            gm%fcQalf(ifc,1) = qas(ix,iy,0)
            gm%fcQbet(ifc,0) = gm%fcQgam(ifc,0)*gm%fcQalf(ifc,0) +
     .                         gm%fcQgam(ifc,1)*gm%fcQalf(ifc,1)
            gm%fcQbet(ifc,1) = gm%fcQgam(ifc,1)*gm%fcQalf(ifc,0) -
     .                         gm%fcQgam(ifc,0)*gm%fcQalf(ifc,1)
            ! wdk account for negative toroidal field through sign of angles
            gm%fcQalf(ifc,1) = gm%fcQalf(ifc,1)*
     .        sign(1.0_R8,wbbl(ix,iy,2))
            gm%fcQbet(ifc,0) = gm%fcQbet(ifc,0)*
     .        sign(1.0_R8,wbbl(ix,iy,2))
            gm%fcPbs(ifc) = pbs(ix,iy,0)
            gm%fcPbshz(ifc) = pbshz(ix,iy,0)
          end if
        end do
      end do

! y faces
      do ix = 0, nx
        do iy = 0, ny
          if(mpg%imapFcy(ix,iy)/=0) then
            ifc = mpg%imapFcy(ix,iy)
            ! wdk sign of field now included in angles
            gm%fcBb(ifc,0:3) = abs(wbbv(ix,iy,0:3))
            gm%fcS(ifc) = gs(ix,iy,1)
            gm%fcHc(ifc,1) = hc(bottomix(ix,iy),bottomiy(ix,iy),3)
            gm%fcHc(ifc,2) = hc(ix,iy,1)
            gm%fcHt(ifc)   = ht(ix,iy,1)
            gm%fcVol(ifc,1) = 0.5_R8*vol(bottomix(ix,iy),
     .                                    bottomiy(ix,iy))
            gm%fcVol(ifc,2) = 0.5_R8*vol(ix,iy)
            gm%fcQgam(ifc,0) = qc(ix,iy,1)
            gm%fcQgam(ifc,1) = qs(ix,iy,1)
            gm%fcQalf(ifc,0) = qac(ix,iy,1)
            gm%fcQalf(ifc,1) = qas(ix,iy,1)
            gm%fcQbet(ifc,0) = gm%fcQgam(ifc,0)*gm%fcQalf(ifc,0) +
     .                         gm%fcQgam(ifc,1)*gm%fcQalf(ifc,1)
            gm%fcQbet(ifc,1) = gm%fcQgam(ifc,1)*gm%fcQalf(ifc,0) -
     .                         gm%fcQgam(ifc,0)*gm%fcQalf(ifc,1)
            ! wdk account for negative toroidal field through sign of angles
            gm%fcQalf(ifc,1) = gm%fcQalf(ifc,1)*
     .        sign(1.0_R8,wbbv(ix,iy,2))
            gm%fcQbet(ifc,0) = gm%fcQbet(ifc,0)*
     .        sign(1.0_R8,wbbv(ix,iy,2))
            gm%fcPbs(ifc) = pbs(ix,iy,1)
            gm%fcPbshz(ifc) = pbshz(ix,iy,1)
          end if
        end do
      end do

! vertices
      gm%vxX = 0._R8
      gm%vxY = 0._R8
      do iy = 0, ny
        do ix = 0, nx
          if(mpg%imapVx(ix,iy)/=0) then
            ivx = mpg%imapVx(ix,iy)
            gm%vxBb(ivx,0:3) = abs(wbbc(ix,iy,0:3))
            gm%vxX(ivx) = crx(ix,iy,0)
            gm%vxY(ivx) = cry(ix,iy,0)
            gm%vxFfbz(ivx) = ffbz(ix,iy,0)
            gm%vxFpsi(ivx) = fpsi(ix,iy,0,0)
          end if
        end do
      end do

#ifdef WG_TODO
#ifdef F2003
      open (newunit=idb,file='debug.out')
#else
      idb = newunit()
      open (unit=idb,file='debug.out')
#endif
      do ix = 236, 238
        do iy = 8, 12
          write (idb,'(2i6)') ix, iy
          write (idb,'(2es17.6)') crx(ix,iy,0),cry(ix,iy,0)
          write (idb,'(2es17.6)') crx(ix,iy,1),cry(ix,iy,1)
          write (idb,'(2es17.6)') crx(ix,iy,3),cry(ix,iy,3)
          write (idb,'(2es17.6)') crx(ix,iy,2),cry(ix,iy,2)
          write (idb,'(2es17.6)') crx(ix,iy,0),cry(ix,iy,0)
       end do
      end do

      write (idb,*)
      do ix = 235, 250
        do iy = 5, 12
          ic = mpg%imapcv(ix,iy)
          write (idb,'(a,3i6)') 'ix,iy,ic ',ix,iy,ic
          write (idb,'(a,5i6)') 'cflags ',cflags(ix,iy,:)
          if (ic == 0) cycle
          write (idb,*) 'faces '
          do i = 1,mpg%cvFcP(ic,2)
            ifc=mpg%cvFc(mpg%cvFcP(ic,1)+i-1)
            ivx1 = mpg%fcVx(ifc,1)
            ivx2 = mpg%fcVx(ifc,2)
            write (idb,'(3i6,4es16.7)') ifc, ivx1, ivx2,
     .        gm%vxX(ivx1), gm%vxY(ivx1),
     .        gm%vxX(ivx2), gm%vxY(ivx2)
          end do
          do i = 1, 7
            write (idb,*)
          end do
        end do
      end do

      write (idb,*)
      write (idb,*)
      write (idb,*) ' cell coordinates '
      do ic = 1, nCv
        write (idb,'(a,3i6)') 'ic ',ic
        do i = 1, mpg%cvVxP(ic,2)
          ivx = mpg%cvVx(mpg%cvVxP(ic,1)+i-1)
          if (ivx == 0) then
            write (idb,'(i6)') ivx
            cycle
          end if
          write (idb,'(i6,2es16.7)') ivx, gm%vxX(ivx),gm%vxY(ivx)
        end do
        write (idb,*)
      end do

      write (idb,*)
      write (idb,*)
      write (idb,*) ' indices of unset vertices '
      do i = 1, nVx
        if (gm%vxX(i)**2+gm%vxY(i)**2 == 2._R8*666._R8*666._R8)
     .    write (idb,'(i6,2es16.7)') i, gm%vxX(i),gm%vxY(i)
      end do
      close (idb)
#endif
! Psi for each flux surface
      deallocate(gm%fsPsi)
      allocate (gm%fsPsi(mpg%nFs))
      gm%fsPsi = 0.0_R8
      do ift = 1, nFs
        ifc = mpg%fsFc(mpg%fsFcP(ift,1))
        gm%fsPsi(ift) = gm%vxFpsi(mpg%fcVx(ifc,1))
      end do
! calculate connection length
      call calc_conn(mpg, gm)

      end subroutine prep_geometry

**********************************************************************

      subroutine prep_plasma (nx, ny, ns, nCv, nCi, nFc, mpg, st)

      implicit none
      integer, intent(in) :: nx, ny, ns, nCi, nCv, nFc
      type(mapping), intent(in) :: mpg
      type(B2State), intent(inout) :: st
      integer :: ix, iy, ic, ifc!, nfx

      call createB2State (nCv, nFc, ns, st)
      call createB2Diagnostic( nCv, ns, nnreg, st % diag )

!      nfx = (nx+1)*ny

      st % ns = ns
      st % ismain = 1

! cell volumes
      do iy = 0, ny-1
        do ix = 0, nx-1
          if(mpg%imapCv(ix,iy)/=0) then
            ic =  mpg%imapCv(ix,iy)
!          ic = (ix+1) + iy*nx
! B2plasma
            st % pl % na(ic,0:ns-1) = na(ix,iy,0:ns-1)
            st % pl % ua(ic,0:ns-1) = ua(ix,iy,0:ns-1)
            st % pl % po(ic) = po(ix,iy)
            st % pl % te(ic) = te(ix,iy)
            st % pl % ti(ic) = ti(ix,iy)
            st % pl % tn(ic) = tn(ix,iy)
            st % pl % kt(ic) = kt(ix,iy)
            st % pl % zt(ic) = zt(ix,iy)

! B2coeff
            st % co % ceqp(ic) = ceqp(ix,iy)

! from b2coeffwork
            st % co % dna0(ic,0:ns-1) = dna0(ix,iy,0:ns-1)
            st % co % dpa0(ic,0:ns-1) = dpa0(ix,iy,0:ns-1)
            st % co % vla0(ic,0:1,0:ns-1) = vla0(ix,iy,0:1,0:ns-1)
            st % co % vsa0(ic,0:ns-1) = vsa0(ix,iy,0:ns-1)
            st % co % hci0(ic) = hci0(ix,iy)
            st % co % hcib(ic,0:ns-1) = hcib(ix,iy,0:ns-1)
            st % co % hce0(ic) = hce0(ix,iy)
            st % co % sig0(ic) = sig0(ix,iy)
            st % co % alf0(ic) = alf0(ix,iy)
            st % co % dkt0(ic) = dkt0(ix,iy)
            st % co % dzt0(ic) = dzt0(ix,iy)
            st % co % dna_ExB(ic) = dna_ExB(ix,iy)
            st % co % hce_ExB(ic) = hce_ExB(ix,iy)
            st % co % hci_ExB(ic) = hci_ExB(ix,iy)
            st % co % cthe(ic,0:ns-1) = cthe(ix,iy,0:ns-1)
            st % co % cthi(ic,0:ns-1) = cthi(ix,iy,0:ns-1)

! kinetic energy
            st % dv % kinrgy(ic,0:ns-1) = kinrgy(ix,iy,0:ns-1)

! derived state quantities
            st % dv % ne(ic) = ne(ix,iy)
            st % dv % ne2(ic) = ne2(ix,iy)
            st % dv % ni(ic,0:1) = ni(ix,iy,0:1)
!          st % dv % ue(ic) = ue(ix,iy)

! residuals
            st % dv % resmo(ic,0:ns-1) = resmo(ix,iy,0:ns-1)
            st % dv % resco(ic,0:ns-1) = resco(ix,iy,0:ns-1)
            st % dv % respo(ic) = respo(ix,iy)
            st % dv % reshe(ic) = reshe(ix,iy)
            st % dv % reshi(ic) = reshi(ix,iy)
!          st % dv % resht(ic) = resht(ix,iy)
            st % dv % resmt(ic) = resmt(ix,iy)
            st % dv % reskt(ic) = reskt(ix,iy)
            st % dv % reszt(ic) = reszt(ix,iy)

! corrections
            st % dv % corua(ic,0:ns-1) = corua(ix,iy,0:ns-1)
            st % dv % corpa(ic,0:ns-1) = corpa(ix,iy,0:ns-1)
            st % dv % corut(ic) = corut(ix,iy)
            st % dv % corpo(ic) = corpo(ix,iy)
!          st % dv % cortt(ic) = cortt(ix,iy)
            st % dv % corte(ic) = corte(ix,iy)
            st % dv % corti(ic) = corti(ix,iy)

! coupling coefficients
            st % dv % pcca(ic,0,0:ns-1) = pcca(ix,iy,0,0,0:ns-1)
            st % dv % pccm(ic,0) = pccm(ix,iy,0,0)
            st % dv % pcca(ic,1,0:ns-1) = pcca(ix,iy,1,1,0:ns-1)
            st % dv % pccm(ic,1) = pccm(ix,iy,1,1)

! B2SourceWork
!          st % srw % sch0(ic,0:3) = sch0(ix,iy,0:3)
!          st % srw % she0(ic,0:3) = she0(ix,iy,0:3)
!          st % srw % shi0(ic,0:3) = shi0(ix,iy,0:3)
!          st % srw % smo0(ic,0:3,0:ns-1) = smo0(ix,iy,0:3,0:ns-1)
!          st % srw % smq0(ic,0:3,0:ns-1) = smq0(ix,iy,0:3,0:ns-1)
!          st % srw % sna0(ic,0:3,0:ns-1) = sna0(ix,iy,0:3,0:ns-1)

! B2rates
            st % rt % rlqa(ic,0:1,0:ns-1) = rlqa(ix,iy,0:1,0:ns-1)
            st % rt % rlrd(ic,0:1,0:ns-1) = rlrd(ix,iy,0:1,0:ns-1)
            st % rt % rlqa(ic,0:1,0:ns-1) = rlqa(ix,iy,0:1,0:ns-1)
            st % rt % rlra(ic,0:1,0:ns-1) = rlra(ix,iy,0:1,0:ns-1)
            st % rt % rlsa(ic,0:1,0:ns-1) = rlsa(ix,iy,0:1,0:ns-1)
            st % rt % rlza(ic,0:1,0:ns-1) = rlza(ix,iy,0:1,0:ns-1)
            st % rt % rlz2(ic,0:1,0:ns-1) = rlz2(ix,iy,0:1,0:ns-1)
            st % rt % rlpt(ic,0:1,0:ns-1) = rlpt(ix,iy,0:1,0:ns-1)
            st % rt % rlpi(ic,0:1,0:ns-1) = rlpi(ix,iy,0:1,0:ns-1)

            st % rt % rza(ic,0:ns-1) = rza(ix,iy,0:ns-1)
            st % rt % rz2(ic,0:ns-1) = rz2(ix,iy,0:ns-1)
            st % rt % rpt(ic,0:ns-1) = rpt(ix,iy,0:ns-1)
            st % rt % rpi(ic,0:ns-1) = rpi(ix,iy,0:ns-1)
          end if
        end do
      end do

! guard cells iy = -1 and iy = ny
      do ix = 0, nx-1
!        ic = nCi + ix+1
        iy = -1
        if(mpg%imapCv(ix,iy)/=0) then
          ic=mpg%imapCv(ix,iy)
! B2plasma
          st % pl % na(ic,0:ns-1) = na(ix,iy,0:ns-1)
          st % pl % ua(ic,0:ns-1) = ua(ix,iy,0:ns-1)
          st % pl % po(ic) = po(ix,iy)
          st % pl % te(ic) = te(ix,iy)
          st % pl % ti(ic) = ti(ix,iy)
          st % pl % tn(ic) = tn(ix,iy)
          st % pl % kt(ic) = kt(ix,iy)
          st % pl % zt(ic) = zt(ix,iy)

! B2coeff
          st % co % ceqp(ic) = ceqp(ix,iy)

! from b2coeffwork
          st % co % dna0(ic,0:ns-1) = dna0(ix,iy,0:ns-1)
          st % co % dpa0(ic,0:ns-1) = dpa0(ix,iy,0:ns-1)
          st % co % vla0(ic,0:1,0:ns-1) = vla0(ix,iy,0:1,0:ns-1)
          st % co % vsa0(ic,0:ns-1) = vsa0(ix,iy,0:ns-1)
          st % co % hci0(ic) = hci0(ix,iy)
          st % co % hcib(ic,0:ns-1) = hcib(ix,iy,0:ns-1)
          st % co % hce0(ic) = hce0(ix,iy)
          st % co % sig0(ic) = sig0(ix,iy)
          st % co % alf0(ic) = alf0(ix,iy)
          st % co % dkt0(ic) = dkt0(ix,iy)
          st % co % dzt0(ic) = dzt0(ix,iy)
          st % co % dna_ExB(ic) = dna_ExB(ix,iy)
          st % co % hce_ExB(ic) = hce_ExB(ix,iy)
          st % co % hci_ExB(ic) = hci_ExB(ix,iy)
          st % co % cthe(ic,0:ns-1) = cthe(ix,iy,0:ns-1)
          st % co % cthi(ic,0:ns-1) = cthi(ix,iy,0:ns-1)

! kinetic energy
          st % dv % kinrgy(ic,0:ns-1) = kinrgy(ix,iy,0:ns-1)

! derived state quantities
          st % dv % ne(ic) = ne(ix,iy)
          st % dv % ne2(ic) = ne2(ix,iy)
          st % dv % ni(ic,0:1) = ni(ix,iy,0:1)
!        st % dv % ue(ic) = ue(ix,iy)

! residuals
          st % dv % resmo(ic,0:ns-1) = resmo(ix,iy,0:ns-1)
          st % dv % resco(ic,0:ns-1) = resco(ix,iy,0:ns-1)
          st % dv % respo(ic) = respo(ix,iy)
          st % dv % reshe(ic) = reshe(ix,iy)
          st % dv % reshi(ic) = reshi(ix,iy)
!        st % dv % resht(ic) = resht(ix,iy)
          st % dv % resmt(ic) = resmt(ix,iy)
          st % dv % reskt(ic) = reskt(ix,iy)
          st % dv % reszt(ic) = reszt(ix,iy)

! corrections
          st % dv % corua(ic,0:ns-1) = corua(ix,iy,0:ns-1)
          st % dv % corpa(ic,0:ns-1) = corpa(ix,iy,0:ns-1)
          st % dv % corut(ic) = corut(ix,iy)
          st % dv % corpo(ic) = corpo(ix,iy)
!        st % dv % cortt(ic) = cortt(ix,iy)
          st % dv % corte(ic) = corte(ix,iy)
          st % dv % corti(ic) = corti(ix,iy)

! coupling coefficients
          st % dv % pcca(ic,0,0:ns-1) = pcca(ix,iy,0,0,0:ns-1)
          st % dv % pccm(ic,0) = pccm(ix,iy,0,0)
          st % dv % pcca(ic,1,0:ns-1) = pcca(ix,iy,1,1,0:ns-1)
          st % dv % pccm(ic,1) = pccm(ix,iy,1,1)

! B2SourceWork
!        st % srw % sch0(ic,0:3) = sch0(ix,iy,0:3)
!        st % srw % she0(ic,0:3) = she0(ix,iy,0:3)
!        st % srw % shi0(ic,0:3) = shi0(ix,iy,0:3)
!        st % srw % smo0(ic,0:3,0:ns-1) = smo0(ix,iy,0:3,0:ns-1)
!        st % srw % smq0(ic,0:3,0:ns-1) = smq0(ix,iy,0:3,0:ns-1)
!        st % srw % sna0(ic,0:3,0:ns-1) = sna0(ix,iy,0:3,0:ns-1)

! B2rates
          st % rt % rlqa(ic,0:1,0:ns-1) = rlqa(ix,iy,0:1,0:ns-1)
          st % rt % rlrd(ic,0:1,0:ns-1) = rlrd(ix,iy,0:1,0:ns-1)
          st % rt % rlqa(ic,0:1,0:ns-1) = rlqa(ix,iy,0:1,0:ns-1)
          st % rt % rlqr(ic,0:1,0:ns-1) = rlqr(ix,iy,0:1,0:ns-1)
          st % rt % rlra(ic,0:1,0:ns-1) = rlra(ix,iy,0:1,0:ns-1)
          st % rt % rlsa(ic,0:1,0:ns-1) = rlsa(ix,iy,0:1,0:ns-1)
          st % rt % rlza(ic,0:1,0:ns-1) = rlza(ix,iy,0:1,0:ns-1)
          st % rt % rlz2(ic,0:1,0:ns-1) = rlz2(ix,iy,0:1,0:ns-1)
          st % rt % rlpt(ic,0:1,0:ns-1) = rlpt(ix,iy,0:1,0:ns-1)
          st % rt % rlpi(ic,0:1,0:ns-1) = rlpi(ix,iy,0:1,0:ns-1)

          st % rt % rza(ic,0:ns-1) = rza(ix,iy,0:ns-1)
          st % rt % rz2(ic,0:ns-1) = rz2(ix,iy,0:ns-1)
          st % rt % rpt(ic,0:ns-1) = rpt(ix,iy,0:ns-1)
          st % rt % rpi(ic,0:ns-1) = rpi(ix,iy,0:ns-1)
        endif
 !       ic = nCi + nx + ix+1
        iy = ny
        if(mpg%imapCv(ix,iy)/=0) then
          ic=mpg%imapCv(ix,iy)
! B2plasma
          st % pl % na(ic,0:ns-1) = na(ix,iy,0:ns-1)
          st % pl % ua(ic,0:ns-1) = ua(ix,iy,0:ns-1)
          st % pl % po(ic) = po(ix,iy)
          st % pl % te(ic) = te(ix,iy)
          st % pl % ti(ic) = ti(ix,iy)
          st % pl % tn(ic) = tn(ix,iy)
          st % pl % kt(ic) = kt(ix,iy)
          st % pl % zt(ic) = zt(ix,iy)

! B2coeff
          st % co % ceqp(ic) = ceqp(ix,iy)

! from b2coeffwork
          st % co % dna0(ic,0:ns-1) = dna0(ix,iy,0:ns-1)
          st % co % dpa0(ic,0:ns-1) = dpa0(ix,iy,0:ns-1)
          st % co % vla0(ic,0:1,0:ns-1) = vla0(ix,iy,0:1,0:ns-1)
          st % co % vsa0(ic,0:ns-1) = vsa0(ix,iy,0:ns-1)
          st % co % hci0(ic) = hci0(ix,iy)
          st % co % hcib(ic,0:ns-1) = hcib(ix,iy,0:ns-1)
          st % co % hce0(ic) = hce0(ix,iy)
          st % co % sig0(ic) = sig0(ix,iy)
          st % co % alf0(ic) = alf0(ix,iy)
          st % co % dkt0(ic) = dkt0(ix,iy)
          st % co % dzt0(ic) = dzt0(ix,iy)
          st % co % dna_ExB(ic) = dna_ExB(ix,iy)
          st % co % hce_ExB(ic) = hce_ExB(ix,iy)
          st % co % hci_ExB(ic) = hci_ExB(ix,iy)
          st % co % cthe(ic,0:ns-1) = cthe(ix,iy,0:ns-1)
          st % co % cthi(ic,0:ns-1) = cthi(ix,iy,0:ns-1)

! kinetic energy
          st % dv % kinrgy(ic,0:ns-1) = kinrgy(ix,iy,0:ns-1)

! derived state quantities
          st % dv % ne(ic) = ne(ix,iy)
          st % dv % ne2(ic) = ne2(ix,iy)
          st % dv % ni(ic,0:1) = ni(ix,iy,0:1)
!        st % dv % ue(ic) = ue(ix,iy)

! residuals
          st % dv % resmo(ic,0:ns-1) = resmo(ix,iy,0:ns-1)
          st % dv % resco(ic,0:ns-1) = resco(ix,iy,0:ns-1)
          st % dv % respo(ic) = respo(ix,iy)
          st % dv % reshe(ic) = reshe(ix,iy)
          st % dv % reshi(ic) = reshi(ix,iy)
!        st % dv % resht(ic) = resht(ix,iy)
          st % dv % resmt(ic) = resmt(ix,iy)
          st % dv % reskt(ic) = reskt(ix,iy)
          st % dv % reszt(ic) = reszt(ix,iy)

! corrections
          st % dv % corua(ic,0:ns-1) = corua(ix,iy,0:ns-1)
          st % dv % corpa(ic,0:ns-1) = corpa(ix,iy,0:ns-1)
          st % dv % corut(ic) = corut(ix,iy)
          st % dv % corpo(ic) = corpo(ix,iy)
!        st % dv % cortt(ic) = cortt(ix,iy)
          st % dv % corte(ic) = corte(ix,iy)
          st % dv % corti(ic) = corti(ix,iy)

! coupling coefficients
          st % dv % pcca(ic,0,0:ns-1) = pcca(ix,iy,0,0,0:ns-1)
          st % dv % pccm(ic,0) = pccm(ix,iy,0,0)
          st % dv % pcca(ic,1,0:ns-1) = pcca(ix,iy,1,1,0:ns-1)
          st % dv % pccm(ic,1) = pccm(ix,iy,1,1)

! B2SourceWork
!        st % srw % sch0(ic,0:3) = sch0(ix,iy,0:3)
!        st % srw % she0(ic,0:3) = she0(ix,iy,0:3)
!        st % srw % shi0(ic,0:3) = shi0(ix,iy,0:3)
!        st % srw % smo0(ic,0:3,0:ns-1) = smo0(ix,iy,0:3,0:ns-1)
!        st % srw % smq0(ic,0:3,0:ns-1) = smq0(ix,iy,0:3,0:ns-1)
!        st % srw % sna0(ic,0:3,0:ns-1) = sna0(ix,iy,0:3,0:ns-1)

! B2rates
          st % rt % rlqa(ic,0:1,0:ns-1) = rlqa(ix,iy,0:1,0:ns-1)
          st % rt % rlrd(ic,0:1,0:ns-1) = rlrd(ix,iy,0:1,0:ns-1)
          st % rt % rlqa(ic,0:1,0:ns-1) = rlqa(ix,iy,0:1,0:ns-1)
          st % rt % rlqr(ic,0:1,0:ns-1) = rlqr(ix,iy,0:1,0:ns-1)
          st % rt % rlra(ic,0:1,0:ns-1) = rlra(ix,iy,0:1,0:ns-1)
          st % rt % rlsa(ic,0:1,0:ns-1) = rlsa(ix,iy,0:1,0:ns-1)
          st % rt % rlza(ic,0:1,0:ns-1) = rlza(ix,iy,0:1,0:ns-1)
          st % rt % rlz2(ic,0:1,0:ns-1) = rlz2(ix,iy,0:1,0:ns-1)
          st % rt % rlpt(ic,0:1,0:ns-1) = rlpt(ix,iy,0:1,0:ns-1)
          st % rt % rlpi(ic,0:1,0:ns-1) = rlpi(ix,iy,0:1,0:ns-1)

          st % rt % rza(ic,0:ns-1) = rza(ix,iy,0:ns-1)
          st % rt % rz2(ic,0:ns-1) = rz2(ix,iy,0:ns-1)
          st % rt % rpt(ic,0:ns-1) = rpt(ix,iy,0:ns-1)
          st % rt % rpi(ic,0:ns-1) = rpi(ix,iy,0:ns-1)
        endif
      end do

! guard cells ix = -1 and ix = nx
      do iy = 0, ny-1
!        ic = nCi + 2*nx + iy+1
        ix = -1
        if(mpg%imapCv(ix,iy)/=0) then
          ic=mpg%imapCv(ix,iy)
! B2plasma
          st % pl % na(ic,0:ns-1) = na(ix,iy,0:ns-1)
          st % pl % ua(ic,0:ns-1) = ua(ix,iy,0:ns-1)
          st % pl % po(ic) = po(ix,iy)
          st % pl % te(ic) = te(ix,iy)
          st % pl % ti(ic) = ti(ix,iy)
          st % pl % tn(ic) = tn(ix,iy)
          st % pl % kt(ic) = kt(ix,iy)
          st % pl % zt(ic) = zt(ix,iy)

! B2coeff
          st % co % ceqp(ic) = ceqp(ix,iy)

! from b2coeffwork
          st % co % dna0(ic,0:ns-1) = dna0(ix,iy,0:ns-1)
          st % co % dpa0(ic,0:ns-1) = dpa0(ix,iy,0:ns-1)
          st % co % vla0(ic,0:1,0:ns-1) = vla0(ix,iy,0:1,0:ns-1)
          st % co % vsa0(ic,0:ns-1) = vsa0(ix,iy,0:ns-1)
          st % co % hci0(ic) = hci0(ix,iy)
          st % co % hcib(ic,0:ns-1) = hcib(ix,iy,0:ns-1)
          st % co % hce0(ic) = hce0(ix,iy)
          st % co % sig0(ic) = sig0(ix,iy)
          st % co % alf0(ic) = alf0(ix,iy)
          st % co % dkt0(ic) = dkt0(ix,iy)
          st % co % dzt0(ic) = dzt0(ix,iy)
          st % co % dna_ExB(ic) = dna_ExB(ix,iy)
          st % co % hce_ExB(ic) = hce_ExB(ix,iy)
          st % co % hci_ExB(ic) = hci_ExB(ix,iy)
          st % co % cthe(ic,0:ns-1) = cthe(ix,iy,0:ns-1)
          st % co % cthi(ic,0:ns-1) = cthi(ix,iy,0:ns-1)

! kinetic energy
          st % dv % kinrgy(ic,0:ns-1) = kinrgy(ix,iy,0:ns-1)

! derived state quantities
          st % dv % ne(ic) = ne(ix,iy)
          st % dv % ne2(ic) = ne2(ix,iy)
          st % dv % ni(ic,0:1) = ni(ix,iy,0:1)
!        st % dv % ue(ic) = ue(ix,iy)

! residuals
          st % dv % resmo(ic,0:ns-1) = resmo(ix,iy,0:ns-1)
          st % dv % resco(ic,0:ns-1) = resco(ix,iy,0:ns-1)
          st % dv % respo(ic) = respo(ix,iy)
          st % dv % reshe(ic) = reshe(ix,iy)
          st % dv % reshi(ic) = reshi(ix,iy)
!        st % dv % resht(ic) = resht(ix,iy)
          st % dv % resmt(ic) = resmt(ix,iy)
          st % dv % reskt(ic) = reskt(ix,iy)
          st % dv % reszt(ic) = reszt(ix,iy)

! corrections
          st % dv % corua(ic,0:ns-1) = corua(ix,iy,0:ns-1)
          st % dv % corpa(ic,0:ns-1) = corpa(ix,iy,0:ns-1)
          st % dv % corut(ic) = corut(ix,iy)
          st % dv % corpo(ic) = corpo(ix,iy)
!        st % dv % cortt(ic) = cortt(ix,iy)
          st % dv % corte(ic) = corte(ix,iy)
          st % dv % corti(ic) = corti(ix,iy)

! coupling coefficients
          st % dv % pcca(ic,0,0:ns-1) = pcca(ix,iy,0,0,0:ns-1)
          st % dv % pccm(ic,0) = pccm(ix,iy,0,0)
          st % dv % pcca(ic,1,0:ns-1) = pcca(ix,iy,1,1,0:ns-1)
          st % dv % pccm(ic,1) = pccm(ix,iy,1,1)

! B2SourceWork
!        st % srw % sch0(ic,0:3) = sch0(ix,iy,0:3)
!        st % srw % she0(ic,0:3) = she0(ix,iy,0:3)
!        st % srw % shi0(ic,0:3) = shi0(ix,iy,0:3)
!        st % srw % smo0(ic,0:3,0:ns-1) = smo0(ix,iy,0:3,0:ns-1)
!        st % srw % smq0(ic,0:3,0:ns-1) = smq0(ix,iy,0:3,0:ns-1)
!        st % srw % sna0(ic,0:3,0:ns-1) = sna0(ix,iy,0:3,0:ns-1)

! B2rates
          st % rt % rlqa(ic,0:1,0:ns-1) = rlqa(ix,iy,0:1,0:ns-1)
          st % rt % rlrd(ic,0:1,0:ns-1) = rlrd(ix,iy,0:1,0:ns-1)
          st % rt % rlqa(ic,0:1,0:ns-1) = rlqa(ix,iy,0:1,0:ns-1)
          st % rt % rlqr(ic,0:1,0:ns-1) = rlqr(ix,iy,0:1,0:ns-1)
          st % rt % rlra(ic,0:1,0:ns-1) = rlra(ix,iy,0:1,0:ns-1)
          st % rt % rlsa(ic,0:1,0:ns-1) = rlsa(ix,iy,0:1,0:ns-1)
          st % rt % rlza(ic,0:1,0:ns-1) = rlza(ix,iy,0:1,0:ns-1)
          st % rt % rlz2(ic,0:1,0:ns-1) = rlz2(ix,iy,0:1,0:ns-1)
          st % rt % rlpt(ic,0:1,0:ns-1) = rlpt(ix,iy,0:1,0:ns-1)
          st % rt % rlpi(ic,0:1,0:ns-1) = rlpi(ix,iy,0:1,0:ns-1)

          st % rt % rza(ic,0:ns-1) = rza(ix,iy,0:ns-1)
          st % rt % rz2(ic,0:ns-1) = rz2(ix,iy,0:ns-1)
          st % rt % rpt(ic,0:ns-1) = rpt(ix,iy,0:ns-1)
          st % rt % rpi(ic,0:ns-1) = rpi(ix,iy,0:ns-1)
        end if

!        ic = nCi + 2*nx + ny + iy+1
        ix = nx
        if(mpg%imapCv(ix,iy)/=0) then
          ic = mpg%imapCv(ix,iy)
! B2plasma
          st % pl % na(ic,0:ns-1) = na(ix,iy,0:ns-1)
          st % pl % ua(ic,0:ns-1) = ua(ix,iy,0:ns-1)
          st % pl % po(ic) = po(ix,iy)
          st % pl % te(ic) = te(ix,iy)
          st % pl % ti(ic) = ti(ix,iy)
          st % pl % tn(ic) = tn(ix,iy)
          st % pl % kt(ic) = kt(ix,iy)
          st % pl % zt(ic) = zt(ix,iy)

! B2coeff
          st % co % ceqp(ic) = ceqp(ix,iy)

! from b2coeffwork
          st % co % dna0(ic,0:ns-1) = dna0(ix,iy,0:ns-1)
          st % co % dpa0(ic,0:ns-1) = dpa0(ix,iy,0:ns-1)
          st % co % vla0(ic,0:1,0:ns-1) = vla0(ix,iy,0:1,0:ns-1)
          st % co % vsa0(ic,0:ns-1) = vsa0(ix,iy,0:ns-1)
          st % co % hci0(ic) = hci0(ix,iy)
          st % co % hcib(ic,0:ns-1) = hcib(ix,iy,0:ns-1)
          st % co % hce0(ic) = hce0(ix,iy)
          st % co % sig0(ic) = sig0(ix,iy)
          st % co % alf0(ic) = alf0(ix,iy)
          st % co % dkt0(ic) = dkt0(ix,iy)
          st % co % dzt0(ic) = dzt0(ix,iy)
          st % co % dna_ExB(ic) = dna_ExB(ix,iy)
          st % co % hce_ExB(ic) = hce_ExB(ix,iy)
          st % co % hci_ExB(ic) = hci_ExB(ix,iy)
          st % co % cthe(ic,0:ns-1) = cthe(ix,iy,0:ns-1)
          st % co % cthi(ic,0:ns-1) = cthi(ix,iy,0:ns-1)

! kinetic energy
          st % dv % kinrgy(ic,0:ns-1) = kinrgy(ix,iy,0:ns-1)

! derived state quantities
          st % dv % ne(ic) = ne(ix,iy)
          st % dv % ne2(ic) = ne2(ix,iy)
          st % dv % ni(ic,0:1) = ni(ix,iy,0:1)
!        st % dv % ue(ic) = ue(ix,iy)

! residuals
          st % dv % resmo(ic,0:ns-1) = resmo(ix,iy,0:ns-1)
          st % dv % resco(ic,0:ns-1) = resco(ix,iy,0:ns-1)
          st % dv % respo(ic) = respo(ix,iy)
          st % dv % reshe(ic) = reshe(ix,iy)
          st % dv % reshi(ic) = reshi(ix,iy)
!        st % dv % resht(ic) = resht(ix,iy)
          st % dv % resmt(ic) = resmt(ix,iy)
          st % dv % reskt(ic) = reskt(ix,iy)
          st % dv % reszt(ic) = reszt(ix,iy)

! corrections
          st % dv % corua(ic,0:ns-1) = corua(ix,iy,0:ns-1)
          st % dv % corpa(ic,0:ns-1) = corpa(ix,iy,0:ns-1)
          st % dv % corut(ic) = corut(ix,iy)
          st % dv % corpo(ic) = corpo(ix,iy)
!        st % dv % cortt(ic) = cortt(ix,iy)
          st % dv % corte(ic) = corte(ix,iy)
          st % dv % corti(ic) = corti(ix,iy)

! coupling coefficients
          st % dv % pcca(ic,0,0:ns-1) = pcca(ix,iy,0,0,0:ns-1)
          st % dv % pccm(ic,0) = pccm(ix,iy,0,0)
          st % dv % pcca(ic,1,0:ns-1) = pcca(ix,iy,1,1,0:ns-1)
          st % dv % pccm(ic,1) = pccm(ix,iy,1,1)

! B2SourceWork
!        st % srw % sch0(ic,0:3) = sch0(ix,iy,0:3)
!        st % srw % she0(ic,0:3) = she0(ix,iy,0:3)
!        st % srw % shi0(ic,0:3) = shi0(ix,iy,0:3)
!        st % srw % smo0(ic,0:3,0:ns-1) = smo0(ix,iy,0:3,0:ns-1)
!        st % srw % smq0(ic,0:3,0:ns-1) = smq0(ix,iy,0:3,0:ns-1)
!        st % srw % sna0(ic,0:3,0:ns-1) = sna0(ix,iy,0:3,0:ns-1)

! B2rates
          st % rt % rlqa(ic,0:1,0:ns-1) = rlqa(ix,iy,0:1,0:ns-1)
          st % rt % rlrd(ic,0:1,0:ns-1) = rlrd(ix,iy,0:1,0:ns-1)
          st % rt % rlqa(ic,0:1,0:ns-1) = rlqa(ix,iy,0:1,0:ns-1)
          st % rt % rlqr(ic,0:1,0:ns-1) = rlqr(ix,iy,0:1,0:ns-1)
          st % rt % rlra(ic,0:1,0:ns-1) = rlra(ix,iy,0:1,0:ns-1)
          st % rt % rlsa(ic,0:1,0:ns-1) = rlsa(ix,iy,0:1,0:ns-1)
          st % rt % rlza(ic,0:1,0:ns-1) = rlza(ix,iy,0:1,0:ns-1)
          st % rt % rlz2(ic,0:1,0:ns-1) = rlz2(ix,iy,0:1,0:ns-1)
          st % rt % rlpt(ic,0:1,0:ns-1) = rlpt(ix,iy,0:1,0:ns-1)
          st % rt % rlpi(ic,0:1,0:ns-1) = rlpi(ix,iy,0:1,0:ns-1)

          st % rt % rza(ic,0:ns-1) = rza(ix,iy,0:ns-1)
          st % rt % rz2(ic,0:ns-1) = rz2(ix,iy,0:ns-1)
          st % rt % rpt(ic,0:ns-1) = rpt(ix,iy,0:ns-1)
          st % rt % rpi(ic,0:ns-1) = rpi(ix,iy,0:ns-1)
        end if
      end do

! x faces
      do iy = 0, ny-1
        do ix = 0, nx
!          ifc = ix+1 + iy*(nx+1)
          if(mpg%imapFcx(ix,iy)/=0) then
            ifc = mpg%imapFcx(ix,iy)
! B2coeff
            st % co % csig_an(ifc,0:1) = csig_an(ix,iy,0,0:1)
            st % co % calf(ifc,0:1) = calf(ix,iy,0,0:1)
            st % co % cdkt(ifc,0:1) = cdkt(ix,iy,0,0:1)
            st % co % cdzt(ifc,0:1) = cdzt(ix,iy,0,0:1)

            st % co % csig(ifc,0:1) = csig(ix,iy,0,0:1)
            st % co % csigin(ifc,0:1,0:ns-1,0:ns-1) =
     .              csigin(ix,iy,0,0:1,0:ns-1,0:ns-1)

            st % co % chve(ifc,0:1) = chve(ix,iy,0,0:1)
            st % co % chce(ifc,0:1) = chce(ix,iy,0,0:1)
            st % co % chvemx(ifc) = chvemx(ix,iy)

            st % co % chvi(ifc,0:1) = chvi(ix,iy,0,0:1)
            st % co % chci(ifc,0:1) = chci(ix,iy,0,0:1)
            st % co % chvimx(ifc) = chvimx(ix,iy)

            st % co % cvla(ifc,0:1,0:ns-1) = cvla(ix,iy,0,0:1,0:ns-1)
            st % co % cdna(ifc,0:1,0:ns-1) = cdna(ix,iy,0,0:1,0:ns-1)
            st % co % cdpa(ifc,0:1,0:ns-1) = cdpa(ix,iy,0,0:1,0:ns-1)

            st % co % cvsa(ifc,0:1,0:ns-1) = cvsa(ix,iy,0,0:1,0:ns-1)
            st % co % cvsa_cl(ifc,0:1,0:ns-1) =
     .                                 cvsa_cl(ix,iy,0,0:1,0:ns-1)

            st % co % cddi(ifc,0:1,0:ns-1) = cddi(ix,iy,0,0:1,0:ns-1)

            st % co % fllim0fna(ifc,0:1,0:ns-1) =
     &              fllim0fna(ix,iy,0,0:1,0:ns-1)
            st % co % fllim0fhi(ifc,0:1,0:ns-1) =
     &              fllim0fhi(ix,iy,0,0:1,0:ns-1)
            st % co % fllimvisc(ifc,0:ns-1) = fllimvisc(ix,iy,0:ns-1)

! fluxes
            st % dv % fch(ifc,0:1) = fch(ix,iy,0,0:1)
            st % dv % fch_32(ifc,0:1) = fch_32(ix,iy,0,0:1)
            st % dv % fch_52(ifc,0:1) = fch_52(ix,iy,0,0:1)
            st % dv % fch_p(ifc,0:1) = fch_p(ix,iy,0,0:1)
            st % dv % fchdia(ifc,0:1) = fchdia(ix,iy,0,0:1)
            st % dv % fchin(ifc,0:1) = fchin(ix,iy,0,0:1)
            st % dv % fchvispar(ifc,0:1) = fchvispar(ix,iy,0,0:1)
            st % dv % fchvisper(ifc,0:1) = fchvisper(ix,iy,0,0:1)
            st % dv % fchvisq(ifc,0:1) = fchvisq(ix,iy,0,0:1)
            st % dv % fchviskt(ifc,0:1) = fchviskt(ix,iy,0,0:1)
            st % dv % fchinert(ifc,0:1) = fchinert(ix,iy,0,0:1)
            st % dv % fchanml(ifc,0:1) = fchanml(ix,iy,0,0:1)

            st % dv % fne(ifc,0:1) = fne(ix,iy,0,0:1)
            st % dv % fne_32(ifc,0:1) = fne_32(ix,iy,0,0:1)
            st % dv % fne_52(ifc,0:1) = fne_52(ix,iy,0,0:1)

            st % dv % fni(ifc,0:1) = fni(ix,iy,0,0:1)
            st % dv % fni_32(ifc,0:1) = fni_32(ix,iy,0,0:1)
            st % dv % fni_52(ifc,0:1) = fni_52(ix,iy,0,0:1)

            st % dv % fna(ifc,0:1,0:ns-1) = fna(ix,iy,0,0:1,0:ns-1)
            st % dv % fna_mdf(ifc,0:1,0:ns-1) =
     .              fna_mdf(ix,iy,0,0:1,0:ns-1)
            st % dv % fna_he(ifc,0:1,0:ns-1)=fna_he(ix,iy,0,0:1,0:ns-1)
            st % dv % fna_fcor(ifc,0:1,0:ns-1) =
     .              fna_fcor(ix,iy,0,0:1,0:ns-1)
            st % dv % fnaPSch(ifc,0:1,0:ns-1) =
     .              fnaPSch(ix,iy,0,0:1,0:ns-1)
            st % dv % fna_32(ifc,0:1,0:ns-1)=fna_32(ix,iy,0,0:1,0:ns-1)
            st % dv % fna_52(ifc,0:1,0:ns-1)=fna_52(ix,iy,0,0:1,0:ns-1)
            st % dv % fna_32nd(ifc,0:1,0:ns-1) =
     .              fna_32nd(ix,iy,0,0:1,0:ns-1)
            st % dv % fna_52nd(ifc,0:1,0:ns-1) =
     .              fna_52nd(ix,iy,0,0:1,0:ns-1)
            st % dv % fna_eir(ifc,0:1,0:ns-1) =
     .              fna_eir(ix,iy,0,0:1,0:ns-1)
            st % dv % fne_eir(ifc,0:1) =
     .              fne_eir(ix,iy,0,0:1)
           st % dv % fna_nodrift(ifc,0:1,0:ns-1) =
     .              fna_nodrift(ix,iy,0,0:1,0:ns-1)

            st % dv % fhe(ifc,0:1) = fhe(ix,iy,0,0:1)
            st % dv % fkt(ifc,0:1) = fkt(ix,iy,0,0:1)
            st % dv % fzt(ifc,0:1) = fzt(ix,iy,0,0:1)
            st % dv % fhe_mdf(ifc,0:1) = fhe_mdf(ix,iy,0,0:1)
            st % dv % fhePSch(ifc,0:1) = fhePSch(ix,iy,0,0:1)
            st % dv % fhet(ifc,0:1) = fhet(ix,iy,0,0:1)
            st % dv % fhe_eir(ifc,0:1) = fhe_eir(ix,iy,0,0:1)
            st % dv % fhi(ifc,0:1) = fhi(ix,iy,0,0:1)
            st % dv % fhi_mdf(ifc,0:1) = fhi_mdf(ix,iy,0,0:1)
            st % dv % fhiPSch(ifc,0:1) = fhiPSch(ix,iy,0,0:1)
            st % dv % fhit(ifc,0:1) = fhit(ix,iy,0,0:1)
            st % dv % fhi_eir(ifc,0:1) = fhi_eir(ix,iy,0,0:1)

! derived quantities
            st % dv % uadia(ifc,0:1,0:ns-1) =
     .              uadia(ix,iy,0,0:1,0:ns-1)
            st % dv % vadia(ifc,0:1,0:ns-1) =
     .               vadia(ix,iy,0,0:1,0:ns-1)
            st % dv % wadia(ifc,0:1,0:ns-1) =
     .               wadia(ix,iy,0,0:1,0:ns-1)
            st % dv % vaecrb(ifc,0:1,0:ns-1) =
     .              vaecrb(ix,iy,0,0:1,0:ns-1)
            st % dv % wedia(ifc,0:1) =
     .              wedia(ix,iy,0,0:1)
            st % dv % vedia(ifc,0:1) =
     .              vedia(ix,iy,0,0:1)
            st % dv % veecrb(ifc,0:1) =
     .              veecrb(ix,iy,0,0:1)
! note that facdrift, fac_ExB, fac_vis used to be defined as cell center quantities,
! but actually used as face quantities in case of drift_flows = 1
            st % dv % facdrift(ifc) = facdrift(ix,iy)
            st % dv % fac_ExB(ifc) = fac_ExB(ix,iy)
            st % dv % fac_vis(ifc) = fac_vis(ix,iy)

! numerical coefficients
            st % dv % conc(ifc,0:1) = conc(ix,iy,0,0:1)
            st % dv % fllime(ifc) = fllime(ix,iy)
            st % dv % floe(ifc,0:1) = floe(ix,iy,0,0:1)
            st % dv % cone(ifc,0:1,:) = cone(ix,iy,0,0:1,:)
            st % dv % fllimi(ifc) = fllimi(ix,iy)
            st % dv % floi(ifc,0:1) = floi(ix,iy,0,0:1)
            st % dv % coni(ifc,0:1,:) = coni(ix,iy,0,0:1,:)
            st % dv % floe_noc(ifc,0:1) = floe_noc(ix,iy,0,0:1)
            st % dv % floi_noc(ifc,0:1) = floi_noc(ix,iy,0,0:1)

!B2PlasmaSnapshot
!          st % psnl % fch(ifc,0:1) = fch(ix,iy,0,0:1)
!          st % psnl % fna(ifc,0:1,0:ns-1) = fna(ix,iy,0,0:1,0:ns-1)
!          st % psnl % fhe(ifc,0:1) = fhe(ix,iy,0,0:1)
!          st % psnl % fhi(ifc,0:1) = fhi(ix,iy,0,0:1)
!          st % psnc % fch(ifc,0:1) = fch(ix,iy,0,0:1)
!          st % psnc % fna(ifc,0:1,0:ns-1) = fna(ix,iy,0,0:1,0:ns-1)
!          st % psnc % fhe(ifc,0:1) = fhe(ix,iy,0,0:1)
!          st % psnc % fhi(ifc,0:1) = fhi(ix,iy,0,0:1)
          end if
        end do
      end do

! y faces
      do ix = 0, nx-1
        do iy = 0, ny
!          ifc = nfx + iy+1 + ix*(ny+1)
          if(mpg%imapFcy(ix,iy)/=0) then
            ifc = mpg%imapFcy(ix,iy)

! B2coeff
          st % co % csig_an(ifc,0:1) = csig_an(ix,iy,1,0:1)
          st % co % calf(ifc,0:1) = calf(ix,iy,1,0:1)
          st % co % cdkt(ifc,0:1) = cdkt(ix,iy,1,0:1)
          st % co % cdzt(ifc,0:1) = cdzt(ix,iy,1,0:1)

          st % co % csig(ifc,0:1) = csig(ix,iy,1,0:1)
          st % co % csigin(ifc,0:1,0:ns-1,0:ns-1) =
     .              csigin(ix,iy,1,0:1,0:ns-1,0:ns-1)

          st % co % chve(ifc,0:1) = chve(ix,iy,1,0:1)
          st % co % chce(ifc,0:1) = chce(ix,iy,1,0:1)
          st % co % chvemx(ifc) = chvemx(ix,iy)

          st % co % chvi(ifc,0:1) = chvi(ix,iy,1,0:1)
          st % co % chci(ifc,0:1) = chci(ix,iy,1,0:1)
          st % co % chvimx(ifc) = chvimx(ix,iy)

          st % co % cvla(ifc,0:1,0:ns-1) = cvla(ix,iy,1,0:1,0:ns-1)
          st % co % cdna(ifc,0:1,0:ns-1) = cdna(ix,iy,1,0:1,0:ns-1)
          st % co % cdpa(ifc,0:1,0:ns-1) = cdpa(ix,iy,1,0:1,0:ns-1)

          st % co % cvsa(ifc,0:1,0:ns-1) = cvsa(ix,iy,1,0:1,0:ns-1)
          st % co % cvsa_cl(ifc,0:1,0:ns-1) =
     .              cvsa_cl(ix,iy,1,0:1,0:ns-1)

          st % co % cddi(ifc,0:1,0:ns-1) = cddi(ix,iy,1,0:1,0:ns-1)

          st % co % fllim0fhi(ifc,0:1,0:ns-1) =
     &              fllim0fhi(ix,iy,1,0:1,0:ns-1)
          st % co % fllimvisc(ifc,0:ns-1) = fllimvisc(ix,iy,0:ns-1)

! fluxes
          st % dv % fch(ifc,0:1) = fch(ix,iy,1,0:1)
          st % dv % fch_32(ifc,0:1) = fch_32(ix,iy,1,0:1)
          st % dv % fch_52(ifc,0:1) = fch_52(ix,iy,1,0:1)
          st % dv % fch_p(ifc,0:1) = fch_p(ix,iy,1,0:1)
          st % dv % fchdia(ifc,0:1) = fchdia(ix,iy,1,0:1)
          st % dv % fchin(ifc,0:1) = fchin(ix,iy,1,0:1)
          st % dv % fchvispar(ifc,0:1) = fchvispar(ix,iy,1,0:1)
          st % dv % fchvisper(ifc,0:1) = fchvisper(ix,iy,1,0:1)
          st % dv % fchvisq(ifc,0:1) = fchvisq(ix,iy,1,0:1)
          st % dv % fchviskt(ifc,0:1) = fchviskt(ix,iy,1,0:1)
          st % dv % fchinert(ifc,0:1) = fchinert(ix,iy,1,0:1)
          st % dv % fchanml(ifc,0:1) = fchanml(ix,iy,1,0:1)

          st % dv % fne(ifc,0:1) = fne(ix,iy,1,0:1)
          st % dv % fne_32(ifc,0:1) = fne_32(ix,iy,1,0:1)
          st % dv % fne_52(ifc,0:1) = fne_52(ix,iy,1,0:1)

          st % dv % fni(ifc,0:1) = fni(ix,iy,1,0:1)
          st % dv % fni_32(ifc,0:1) = fni_32(ix,iy,1,0:1)
          st % dv % fni_52(ifc,0:1) = fni_52(ix,iy,1,0:1)

          st % dv % fna(ifc,0:1,0:ns-1) = fna(ix,iy,1,0:1,0:ns-1)
          st % dv % fna_mdf(ifc,0:1,0:ns-1) =
     .              fna_mdf(ix,iy,1,0:1,0:ns-1)
          st % dv % fna_he(ifc,0:1,0:ns-1) = fna_he(ix,iy,1,0:1,0:ns-1)
          st % dv % fna_fcor(ifc,0:1,0:ns-1) =
     .              fna_fcor(ix,iy,1,0:1,0:ns-1)
          st % dv % fnaPSch(ifc,0:1,0:ns-1) =
     .              fnaPSch(ix,iy,1,0:1,0:ns-1)
          st % dv % fna_32(ifc,0:1,0:ns-1) = fna_32(ix,iy,1,0:1,0:ns-1)
          st % dv % fna_52(ifc,0:1,0:ns-1) = fna_52(ix,iy,1,0:1,0:ns-1)
          st % dv % fna_32nd(ifc,0:1,0:ns-1) =
     .              fna_32nd(ix,iy,1,0:1,0:ns-1)
          st % dv % fna_52nd(ifc,0:1,0:ns-1) =
     .              fna_52nd(ix,iy,1,0:1,0:ns-1)
          st % dv % fna_eir(ifc,0:1,0:ns-1) =
     .              fna_eir(ix,iy,1,0:1,0:ns-1)
          st % dv % fne_eir(ifc,0:1) =
     .              fne_eir(ix,iy,1,0:1)
          st % dv % fna_nodrift(ifc,0:1,0:ns-1) =
     .              fna_nodrift(ix,iy,1,0:1,0:ns-1)

          st % dv % fhe(ifc,0:1) = fhe(ix,iy,1,0:1)
          st % dv % fkt(ifc,0:1) = fkt(ix,iy,1,0:1)
          st % dv % fzt(ifc,0:1) = fzt(ix,iy,1,0:1)
          st % dv % fhe_mdf(ifc,0:1) = fhe_mdf(ix,iy,1,0:1)
          st % dv % fhePSch(ifc,0:1) = fhePSch(ix,iy,1,0:1)
          st % dv % fhet(ifc,0:1) = fhet(ix,iy,1,0:1)
          st % dv % fhe_eir(ifc,0:1) = fhe_eir(ix,iy,1,0:1)
          st % dv % fhi(ifc,0:1) = fhi(ix,iy,1,0:1)
          st % dv % fhi_mdf(ifc,0:1) = fhi_mdf(ix,iy,1,0:1)
          st % dv % fhiPSch(ifc,0:1) = fhiPSch(ix,iy,1,0:1)
          st % dv % fhit(ifc,0:1) = fhit(ix,iy,1,0:1)
          st % dv % fhi_eir(ifc,0:1) = fhi_eir(ix,iy,1,0:1)

! derived quantities
          st % dv % uadia(ifc,0:1,0:ns-1) =
     .              uadia(ix,iy,1,0:1,0:ns-1)
          st % dv % vadia(ifc,0:1,0:ns-1) =
     .              vadia(ix,iy,1,0:1,0:ns-1)
          st % dv % wadia(ifc,0:1,0:ns-1) =
     .              wadia(ix,iy,1,0:1,0:ns-1)
          st % dv % vaecrb(ifc,0:1,0:ns-1) =
     .              vaecrb(ix,iy,1,0:1,0:ns-1)
          st % dv % wedia(ifc,0:1) =
     .              wedia(ix,iy,1,0:1)
          st % dv % vedia(ifc,0:1) =
     .              vedia(ix,iy,1,0:1)
          st % dv % veecrb(ifc,0:1) =
     .              veecrb(ix,iy,1,0:1)
! note that facdrift, fac_ExB, fac_vis used to be defined as cell center quantities,
! but actually used as face quantities in case of drift_flows = 1
          st % dv % facdrift(ifc) = facdrift(ix,iy)
          st % dv % fac_ExB(ifc) = fac_ExB(ix,iy)
          st % dv % fac_vis(ifc) = fac_vis(ix,iy)

! numerical coefficients
          st % dv % conc(ifc,0:1) = conc(ix,iy,1,0:1)
          st % dv % fllime(ifc) = fllime(ix,iy)
          st % dv % floe(ifc,0:1) = floe(ix,iy,1,0:1)
          st % dv % cone(ifc,0:1,:) = cone(ix,iy,1,0:1,:)
          st % dv % fllimi(ifc) = fllimi(ix,iy)
          st % dv % floi(ifc,0:1) = floi(ix,iy,1,0:1)
          st % dv % coni(ifc,0:1,:) = coni(ix,iy,1,0:1,:)
          st % dv % floe_noc(ifc,0:1) = floe_noc(ix,iy,1,0:1)
          st % dv % floi_noc(ifc,0:1) = floi_noc(ix,iy,1,0:1)

! B2PlasmaSnapshot
!          st % psnl % fch(ifc,0:1) = fch(ix,iy,1,0:1)
!          st % psnl % fna(ifc,0:1,0:ns-1) = fna(ix,iy,1,0:1,0:ns-1)
!          st % psnl % fhe(ifc,0:1) = fhe(ix,iy,1,0:1)
!          st % psnl % fhi(ifc,0:1) = fhi(ix,iy,1,0:1)
!          st % psnc % fch(ifc,0:1) = fch(ix,iy,1,0:1)
!         st % psnc % fna(ifc,0:1,0:ns-1) = fna(ix,iy,1,0:1,0:ns-1)
!          st % psnc % fhe(ifc,0:1) = fhe(ix,iy,1,0:1)
!          st % psnc % fhi(ifc,0:1) = fhi(ix,iy,1,0:1)
          end if
        end do
      end do

      st % diag % aresco(0:ns-1) = aresco(0:ns-1)
      st % diag % aresmo(0:ns-1) = aresmo(0:ns-1)

!      st % diag % acorpa(0:ns-1) = acorpa(0:ns-1)
!      st % diag % acorua(0:ns-1) = acorua(0:ns-1)

!      st % diag % rescoreg(0:nnreg(0),0:ns-1) =
!     .            rescoreg(0:nnreg(0),0:ns-1)
!      st % diag % resmoreg(0:nnreg(0),0:ns-1) =
!     .            resmoreg(0:nnreg(0),0:ns-1)
!      st % diag % reshereg(0:nnreg(0)) = reshereg(0:nnreg(0))
!      st % diag % reshireg(0:nnreg(0)) = reshireg(0:nnreg(0))

      end subroutine prep_plasma

**********************************************************************

      subroutine prep_state (nx, ny, ns, nCv, nCi, nFc, mpg, st)
      !! convert only a limited set of quantities as needed for the
      !! creation of the b2fstate file
      implicit none
      integer, intent(in) :: nx, ny, ns, nCi, nCv, nFc
      type(mapping), intent(in) :: mpg
      type(B2State), intent(inout) :: st
      integer :: ix, iy, ic, ifc!, nfx

      call createB2Plasma( nCv, ns, st % pl )
      call createB2Derivatives( nCv, nFc, ns, st % dv )

!      nfx = (nx+1)*ny

      st % ns = ns
      st % ismain = 1

! cell volumes
      do iy = -1, ny
        do ix = -1, nx
          if(mpg%imapCv(ix,iy)/=0) then
            ic =  mpg%imapCv(ix,iy)
! B2plasma
            st % pl % na(ic,0:ns-1) = na(ix,iy,0:ns-1)
            st % pl % ua(ic,0:ns-1) = ua(ix,iy,0:ns-1)
            st % pl % po(ic) = po(ix,iy)
            st % pl % te(ic) = te(ix,iy)
            st % pl % ti(ic) = ti(ix,iy)
            st % pl % tn(ic) = tn(ix,iy)
            st % pl % kt(ic) = kt(ix,iy)
            st % pl % zt(ic) = zt(ix,iy)

! derived state quantities
            st % dv % ne(ic) = ne(ix,iy)

          end if
        end do
      end do

! x faces
      do iy = 0, ny-1
        do ix = 0, nx
          if(mpg%imapFcx(ix,iy)/=0) then
            ifc = mpg%imapFcx(ix,iy)
! fluxes
            st % dv % fch(ifc,0:1) = fch(ix,iy,0,0:1)
            st % dv % fch_32(ifc,0:1) = fch_32(ix,iy,0,0:1)
            st % dv % fch_52(ifc,0:1) = fch_52(ix,iy,0,0:1)
            st % dv % fch_p(ifc,0:1) = fch_p(ix,iy,0,0:1)
            st % dv % fchdia(ifc,0:1) = fchdia(ix,iy,0,0:1)
            st % dv % fchin(ifc,0:1) = fchin(ix,iy,0,0:1)
            st % dv % fchvispar(ifc,0:1) = fchvispar(ix,iy,0,0:1)
            st % dv % fchvisper(ifc,0:1) = fchvisper(ix,iy,0,0:1)
            st % dv % fchvisq(ifc,0:1) = fchvisq(ix,iy,0,0:1)
            st % dv % fchviskt(ifc,0:1) = fchviskt(ix,iy,0,0:1)
            st % dv % fchinert(ifc,0:1) = fchinert(ix,iy,0,0:1)
            st % dv % fchanml(ifc,0:1) = fchanml(ix,iy,0,0:1)

            st % dv % fne(ifc,0:1) = fne(ix,iy,0,0:1)
            st % dv % fne_32(ifc,0:1) = fne_32(ix,iy,0,0:1)
            st % dv % fne_52(ifc,0:1) = fne_52(ix,iy,0,0:1)

            st % dv % fni(ifc,0:1) = fni(ix,iy,0,0:1)
            st % dv % fni_32(ifc,0:1) = fni_32(ix,iy,0,0:1)
            st % dv % fni_52(ifc,0:1) = fni_52(ix,iy,0,0:1)

            st % dv % fna(ifc,0:1,0:ns-1) = fna(ix,iy,0,0:1,0:ns-1)
            st % dv % fna_mdf(ifc,0:1,0:ns-1) =
     .              fna_mdf(ix,iy,0,0:1,0:ns-1)
            st % dv % fna_he(ifc,0:1,0:ns-1)=fna_he(ix,iy,0,0:1,0:ns-1)
            st % dv % fna_fcor(ifc,0:1,0:ns-1) =
     .              fna_fcor(ix,iy,0,0:1,0:ns-1)
            st % dv % fnaPSch(ifc,0:1,0:ns-1) =
     .              fnaPSch(ix,iy,0,0:1,0:ns-1)
            st % dv % fna_32(ifc,0:1,0:ns-1)=fna_32(ix,iy,0,0:1,0:ns-1)
            st % dv % fna_52(ifc,0:1,0:ns-1)=fna_52(ix,iy,0,0:1,0:ns-1)
            st % dv % fna_32nd(ifc,0:1,0:ns-1) =
     .              fna_32nd(ix,iy,0,0:1,0:ns-1)
            st % dv % fna_52nd(ifc,0:1,0:ns-1) =
     .              fna_52nd(ix,iy,0,0:1,0:ns-1)
            st % dv % fna_eir(ifc,0:1,0:ns-1) =
     .              fna_eir(ix,iy,0,0:1,0:ns-1)
            st % dv % fne_eir(ifc,0:1) =
     .              fne_eir(ix,iy,0,0:1)
            st % dv % fna_nodrift(ifc,0:1,0:ns-1) =
     .              fna_nodrift(ix,iy,0,0:1,0:ns-1)

            st % dv % fhe(ifc,0:1) = fhe(ix,iy,0,0:1)
            st % dv % fkt(ifc,0:1) = fkt(ix,iy,0,0:1)
            st % dv % fzt(ifc,0:1) = fzt(ix,iy,0,0:1)
            st % dv % fhe_mdf(ifc,0:1) = fhe_mdf(ix,iy,0,0:1)
            st % dv % fhePSch(ifc,0:1) = fhePSch(ix,iy,0,0:1)
            st % dv % fhet(ifc,0:1) = fhet(ix,iy,0,0:1)
            st % dv % fhe_eir(ifc,0:1) = fhe_eir(ix,iy,0,0:1)
            st % dv % fhi(ifc,0:1) = fhi(ix,iy,0,0:1)
            st % dv % fhi_mdf(ifc,0:1) = fhi_mdf(ix,iy,0,0:1)
            st % dv % fhiPSch(ifc,0:1) = fhiPSch(ix,iy,0,0:1)
            st % dv % fhit(ifc,0:1) = fhit(ix,iy,0,0:1)
            st % dv % fhi_eir(ifc,0:1) = fhi_eir(ix,iy,0,0:1)

! derived quantities
            st % dv % uadia(ifc,0:1,0:ns-1) =
     .              uadia(ix,iy,0,0:1,0:ns-1)
            st % dv % vadia(ifc,0:1,0:ns-1) =
     .               vadia(ix,iy,0,0:1,0:ns-1)
            st % dv % wadia(ifc,0:1,0:ns-1) =
     .               wadia(ix,iy,0,0:1,0:ns-1)
            st % dv % vaecrb(ifc,0:1,0:ns-1) =
     .              vaecrb(ix,iy,0,0:1,0:ns-1)
            st % dv % wedia(ifc,0:1) =
     .              wedia(ix,iy,0,0:1)
            st % dv % vedia(ifc,0:1) =
     .              vedia(ix,iy,0,0:1)
            st % dv % veecrb(ifc,0:1) =
     .              veecrb(ix,iy,0,0:1)
! note that facdrift, fac_ExB, fac_vis used to be defined as cell center quantities,
! but actually used as face quantities in case of drift_flows = 1
            st % dv % facdrift(ifc) = facdrift(ix,iy)
            st % dv % fac_ExB(ifc) = fac_ExB(ix,iy)
            st % dv % fac_vis(ifc) = fac_vis(ix,iy)

! numerical coefficients
            st % dv % floe_noc(ifc,0:1) = floe_noc(ix,iy,0,0:1)
            st % dv % floi_noc(ifc,0:1) = floi_noc(ix,iy,0,0:1)

          end if
        end do
      end do

! y faces
      do ix = 0, nx-1
        do iy = 0, ny
!          ifc = nfx + iy+1 + ix*(ny+1)
          if(mpg%imapFcy(ix,iy)/=0) then
            ifc = mpg%imapFcy(ix,iy)


! fluxes
          st % dv % fch(ifc,0:1) = fch(ix,iy,1,0:1)
          st % dv % fch_32(ifc,0:1) = fch_32(ix,iy,1,0:1)
          st % dv % fch_52(ifc,0:1) = fch_52(ix,iy,1,0:1)
          st % dv % fch_p(ifc,0:1) = fch_p(ix,iy,1,0:1)
          st % dv % fchdia(ifc,0:1) = fchdia(ix,iy,1,0:1)
          st % dv % fchin(ifc,0:1) = fchin(ix,iy,1,0:1)
          st % dv % fchvispar(ifc,0:1) = fchvispar(ix,iy,1,0:1)
          st % dv % fchvisper(ifc,0:1) = fchvisper(ix,iy,1,0:1)
          st % dv % fchvisq(ifc,0:1) = fchvisq(ix,iy,1,0:1)
          st % dv % fchviskt(ifc,0:1) = fchviskt(ix,iy,1,0:1)
          st % dv % fchinert(ifc,0:1) = fchinert(ix,iy,1,0:1)
          st % dv % fchanml(ifc,0:1) = fchanml(ix,iy,1,0:1)

          st % dv % fne(ifc,0:1) = fne(ix,iy,1,0:1)
          st % dv % fne_32(ifc,0:1) = fne_32(ix,iy,1,0:1)
          st % dv % fne_52(ifc,0:1) = fne_52(ix,iy,1,0:1)

          st % dv % fni(ifc,0:1) = fni(ix,iy,1,0:1)
          st % dv % fni_32(ifc,0:1) = fni_32(ix,iy,1,0:1)
          st % dv % fni_52(ifc,0:1) = fni_52(ix,iy,1,0:1)

          st % dv % fna(ifc,0:1,0:ns-1) = fna(ix,iy,1,0:1,0:ns-1)
          st % dv % fna_mdf(ifc,0:1,0:ns-1) =
     .              fna_mdf(ix,iy,1,0:1,0:ns-1)
          st % dv % fna_he(ifc,0:1,0:ns-1) = fna_he(ix,iy,1,0:1,0:ns-1)
          st % dv % fna_fcor(ifc,0:1,0:ns-1) =
     .              fna_fcor(ix,iy,1,0:1,0:ns-1)
          st % dv % fnaPSch(ifc,0:1,0:ns-1) =
     .              fnaPSch(ix,iy,1,0:1,0:ns-1)
          st % dv % fna_32(ifc,0:1,0:ns-1) = fna_32(ix,iy,1,0:1,0:ns-1)
          st % dv % fna_52(ifc,0:1,0:ns-1) = fna_52(ix,iy,1,0:1,0:ns-1)
          st % dv % fna_32nd(ifc,0:1,0:ns-1) =
     .              fna_32nd(ix,iy,1,0:1,0:ns-1)
          st % dv % fna_52nd(ifc,0:1,0:ns-1) =
     .              fna_52nd(ix,iy,1,0:1,0:ns-1)
          st % dv % fna_eir(ifc,0:1,0:ns-1) =
     .              fna_eir(ix,iy,1,0:1,0:ns-1)
          st % dv % fne_eir(ifc,0:1) =
     .              fne_eir(ix,iy,1,0:1)
          st % dv % fna_nodrift(ifc,0:1,0:ns-1) =
     .              fna_nodrift(ix,iy,1,0:1,0:ns-1)

          st % dv % fhe(ifc,0:1) = fhe(ix,iy,1,0:1)
          st % dv % fkt(ifc,0:1) = fkt(ix,iy,1,0:1)
          st % dv % fzt(ifc,0:1) = fzt(ix,iy,1,0:1)
          st % dv % fhe_mdf(ifc,0:1) = fhe_mdf(ix,iy,1,0:1)
          st % dv % fhePSch(ifc,0:1) = fhePSch(ix,iy,1,0:1)
          st % dv % fhet(ifc,0:1) = fhet(ix,iy,1,0:1)
          st % dv % fhe_eir(ifc,0:1) = fhe_eir(ix,iy,1,0:1)
          st % dv % fhi(ifc,0:1) = fhi(ix,iy,1,0:1)
          st % dv % fhi_mdf(ifc,0:1) = fhi_mdf(ix,iy,1,0:1)
          st % dv % fhiPSch(ifc,0:1) = fhiPSch(ix,iy,1,0:1)
          st % dv % fhit(ifc,0:1) = fhit(ix,iy,1,0:1)
          st % dv % fhi_eir(ifc,0:1) = fhi_eir(ix,iy,1,0:1)

! derived quantities
          st % dv % uadia(ifc,0:1,0:ns-1) =
     .              uadia(ix,iy,1,0:1,0:ns-1)
          st % dv % vadia(ifc,0:1,0:ns-1) =
     .              vadia(ix,iy,1,0:1,0:ns-1)
          st % dv % wadia(ifc,0:1,0:ns-1) =
     .              wadia(ix,iy,1,0:1,0:ns-1)
          st % dv % vaecrb(ifc,0:1,0:ns-1) =
     .              vaecrb(ix,iy,1,0:1,0:ns-1)
          st % dv % wedia(ifc,0:1) =
     .              wedia(ix,iy,1,0:1)
          st % dv % vedia(ifc,0:1) =
     .              vedia(ix,iy,1,0:1)
          st % dv % veecrb(ifc,0:1) =
     .              veecrb(ix,iy,1,0:1)
! note that facdrift, fac_ExB, fac_vis used to be defined as cell center quantities,
! but actually used as face quantities in case of drift_flows = 1
          st % dv % facdrift(ifc) = facdrift(ix,iy)
          st % dv % fac_ExB(ifc) = fac_ExB(ix,iy)
          st % dv % fac_vis(ifc) = fac_vis(ix,iy)

! numerical coefficients
          st % dv % floe_noc(ifc,0:1) = floe_noc(ix,iy,1,0:1)
          st % dv % floi_noc(ifc,0:1) = floi_noc(ix,iy,1,0:1)

          end if
        end do
      end do

      end subroutine prep_state

**********************************************************************

      subroutine to_struct_plasma (ns, mpg, st)

      implicit none
      integer, intent(in) :: ns
      type(mapping), intent(in) :: mpg
      type(B2State), intent(in) :: st
      integer :: nx, ny, nCi, nCv, nFc
      integer :: ix, iy, ic, ifc
      integer :: nout_st, io, idum(0:9)
      character lblst*120, lblps*120
#ifndef F2003
      integer newunit
      external newunit
#endif

      nCv = mpg%nCv
      nCi = mpg%nCi
      nFc = mpg%nFc

      if (mpg%isClassicalGrid == 1) then
        nx = mpg%nx
        ny = mpg%ny
      else
        call xerrab ('Not a classical grid!')
      endif

*   ..allocate modules for structured grid
      if (.not.allocated(fch))  call alloc_b2mod_plasma(nx,ny,ns)
      if (.not.allocated(ceqp)) call alloc_b2mod_transport(nx,ny,ns)
      if (.not.allocated(dna0))
     .    call alloc_b2mod_anomalous_transport(nx,ny,ns)
      if (.not.allocated(conb))  call alloc_b2mod_work(nx,ny,ns)
      if (.not.allocated(resco)) call alloc_b2mod_residuals(nx,ny,ns)
      if (.not.allocated(rlqa))  call alloc_b2mod_rates(nx,ny,ns)
      if (.not.allocated(sch))   call alloc_b2mod_sources(nx,ny,ns)

! cell volumes
      do iy = 0, ny-1
        do ix = 0, nx-1
          if(mpg%imapCv(ix,iy)/=0) then
            ic =  mpg%imapCv(ix,iy)
! B2plasma
          na(ix,iy,0:ns-1) = st % pl % na(ic,0:ns-1)
          ua(ix,iy,0:ns-1) = st % pl % ua(ic,0:ns-1)
          po(ix,iy)        = st % pl % po(ic)
          te(ix,iy)        = st % pl % te(ic)
          ti(ix,iy)        = st % pl % ti(ic)
          tn(ix,iy)        = st % pl % tn(ic)
          kt(ix,iy)        = st % pl % kt(ic)
          zt(ix,iy)        = st % pl % zt(ic)

! B2coeff
          ceqp(ix,iy) = st % co % ceqp(ic)

! from b2coeffwork
          dna0(ix,iy,0:ns-1)     = st % co % dna0(ic,0:ns-1)
          dpa0(ix,iy,0:ns-1)     = st % co % dpa0(ic,0:ns-1)
          vla0(ix,iy,0:1,0:ns-1) = st % co % vla0(ic,0:1,0:ns-1)
          vsa0(ix,iy,0:ns-1)     = st % co % vsa0(ic,0:ns-1)
          hci0(ix,iy)            = st % co % hci0(ic)
          hcn0(ix,iy)            = st % co % hcn0(ic)
          hcib(ix,iy,0:ns-1)     = st % co % hcib(ic,0:ns-1)
          hce0(ix,iy)            = st % co % hce0(ic)
          sig0(ix,iy)            = st % co % sig0(ic)
          alf0(ix,iy)            = st % co % alf0(ic)
          dkt0(ix,iy)            = st % co % dkt0(ic)
          dzt0(ix,iy)            = st % co % dzt0(ic)
          dna_ExB(ix,iy)         = st % co % dna_ExB(ic)
          hce_ExB(ix,iy)         = st % co % hce_ExB(ic)
          hci_ExB(ix,iy)         = st % co % hci_ExB(ic)
          cthe(ix,iy,0:ns-1)     = st % co % cthe(ic,0:ns-1)
          cthi(ix,iy,0:ns-1)     = st % co % cthi(ic,0:ns-1)


! kinetic energy
          kinrgy(ix,iy,0:ns-1) = st % dv % kinrgy(ic,0:ns-1)

! derived state quantities
          ne(ix,iy)     = st % dv % ne(ic)
          ne2(ix,iy)    = st % dv % ne2(ic)
          ni(ix,iy,0:1) = st % dv % ni(ic,0:1)
!         ue(ix,iy)     = st % dv % ue(ic)

! residuals
          resmo(ix,iy,0:ns-1) = st % dv % resmo(ic,0:ns-1)
          resco(ix,iy,0:ns-1) = st % dv % resco(ic,0:ns-1)
          respo(ix,iy)        = st % dv % respo(ic)
          reshe(ix,iy)        = st % dv % reshe(ic)
          reshi(ix,iy)        = st % dv % reshi(ic)
          reshn(ix,iy)        = st % dv % reshn(ic)
!         resht(ix,iy)        = st % dv % resht(ic)
          resmt(ix,iy)        = st % dv % resmt(ic)
          reskt(ix,iy)        = st % dv % reskt(ic)
          reszt(ix,iy)        = st % dv % reszt(ic)

! corrections
          corua(ix,iy,0:ns-1) = st % dv % corua(ic,0:ns-1)
          corpa(ix,iy,0:ns-1) = st % dv % corpa(ic,0:ns-1)
          corut(ix,iy)        = st % dv % corut(ic)
          corpo(ix,iy)        = st % dv % corpo(ic)
!         cortt(ix,iy)        = st % dv % cortt(ic)
          corte(ix,iy)        = st % dv % corte(ic)
          corti(ix,iy)        = st % dv % corti(ic)
          cortn(ix,iy)        = st % dv % cortn(ic)

! coupling coefficients
          pcca(ix,iy,0,0,0:ns-1) = st % dv % pcca(ic,0,0:ns-1)
          pccm(ix,iy,0,0)        = st % dv % pccm(ic,0)
          pcca(ix,iy,1,1,0:ns-1) = st % dv % pcca(ic,1,0:ns-1)
          pccm(ix,iy,1,1)        = st % dv % pccm(ic,1)

! B2Source
          sch(ix,iy,0:3)        = st % sr % sch(ic,0:3)
          she(ix,iy,0:3)        = st % sr % she(ic,0:3)
          shi(ix,iy,0:3)        = st % sr % shi(ic,0:3)
          shn(ix,iy,0:3)        = st % sr % shn(ic,0:3)
          skt(ix,iy,0:3)        = st % sr % skt(ic,0:3)
          skt_prod(ix,iy)       = st % sr % skt_prod(ic)
          skt_diss(ix,iy)       = st % sr % skt_diss(ic)
          szt(ix,iy,0:3)        = st % sr % szt(ic,0:3)
          smo(ix,iy,0:3,0:ns-1) = st % sr % smo(ic,0:3,0:ns-1)
          smq(ix,iy,0:3,0:ns-1) = st % sr % smq(ic,0:3,0:ns-1)
          sna(ix,iy,0:1,0:ns-1) = st % sr % sna(ic,0:1,0:ns-1)

! B2rates
          rlqa(ix,iy,0:1,0:ns-1) = st % rt % rlqa(ic,0:1,0:ns-1)
          rlrd(ix,iy,0:1,0:ns-1) = st % rt % rlrd(ic,0:1,0:ns-1)
          rlqa(ix,iy,0:1,0:ns-1) = st % rt % rlqa(ic,0:1,0:ns-1)
          rlqr(ix,iy,0:1,0:ns-1) = st % rt % rlqr(ic,0:1,0:ns-1)
          rlra(ix,iy,0:1,0:ns-1) = st % rt % rlra(ic,0:1,0:ns-1)
          rlsa(ix,iy,0:1,0:ns-1) = st % rt % rlsa(ic,0:1,0:ns-1)
          rlza(ix,iy,0:1,0:ns-1) = st % rt % rlza(ic,0:1,0:ns-1)
          rlz2(ix,iy,0:1,0:ns-1) = st % rt % rlz2(ic,0:1,0:ns-1)
          rlpt(ix,iy,0:1,0:ns-1) = st % rt % rlpt(ic,0:1,0:ns-1)
          rlpi(ix,iy,0:1,0:ns-1) = st % rt % rlpi(ic,0:1,0:ns-1)

          rza(ix,iy,0:ns-1) = st % rt % rza(ic,0:ns-1)
          rz2(ix,iy,0:ns-1) = st % rt % rz2(ic,0:ns-1)
          rpt(ix,iy,0:ns-1) = st % rt % rpt(ic,0:ns-1)
          rpi(ix,iy,0:ns-1) = st % rt % rpi(ic,0:ns-1)
          endif
        end do
      end do

! gueard cells iy = -1 and iy = ny
      do ix = 0, nx-1
        iy = -1
        if(mpg%imapCv(ix,iy)/=0) then
          ic =  mpg%imapCv(ix,iy)

! B2plasma
        na(ix,iy,0:ns-1) = st % pl % na(ic,0:ns-1)
        ua(ix,iy,0:ns-1) = st % pl % ua(ic,0:ns-1)
        po(ix,iy)        = st % pl % po(ic)
        te(ix,iy)        = st % pl % te(ic)
        ti(ix,iy)        = st % pl % ti(ic)
        tn(ix,iy)        = st % pl % tn(ic)
        kt(ix,iy)        = st % pl % kt(ic)
        zt(ix,iy)        = st % pl % zt(ic)

! B2coeff
        ceqp(ix,iy) = st % co % ceqp(ic)

! from b2coeffwork
        dna0(ix,iy,0:ns-1)     = st % co % dna0(ic,0:ns-1)
        dpa0(ix,iy,0:ns-1)     = st % co % dpa0(ic,0:ns-1)
        vla0(ix,iy,0:1,0:ns-1) = st % co % vla0(ic,0:1,0:ns-1)
        vsa0(ix,iy,0:ns-1)     = st % co % vsa0(ic,0:ns-1)
        hci0(ix,iy)            = st % co % hci0(ic)
        hcn0(ix,iy)            = st % co % hcn0(ic)
        hcib(ix,iy,0:ns-1)     = st % co % hcib(ic,0:ns-1)
        hce0(ix,iy)            = st % co % hce0(ic)
        sig0(ix,iy)            = st % co % sig0(ic)
        alf0(ix,iy)            = st % co % alf0(ic)
        dkt0(ix,iy)            = st % co % dkt0(ic)
        dzt0(ix,iy)            = st % co % dzt0(ic)
        dna_ExB(ix,iy)         = st % co % dna_ExB(ic)
        hce_ExB(ix,iy)         = st % co % hce_ExB(ic)
        hci_ExB(ix,iy)         = st % co % hci_ExB(ic)
        cthe(ix,iy,0:ns-1)     = st % co % cthe(ic,0:ns-1)
        cthi(ix,iy,0:ns-1)     = st % co % cthi(ic,0:ns-1)

! kinetic energy
        kinrgy(ix,iy,0:ns-1) = st % dv % kinrgy(ic,0:ns-1)

! derived state quantities
        ne(ix,iy)     = st % dv % ne(ic)
        ne2(ix,iy)    = st % dv % ne2(ic)
        ni(ix,iy,0:1) = st % dv % ni(ic,0:1)
!       ue(ix,iy)     = st % dv % ue(ic)

! residuals
        resmo(ix,iy,0:ns-1) = st % dv % resmo(ic,0:ns-1)
        resco(ix,iy,0:ns-1) = st % dv % resco(ic,0:ns-1)
        respo(ix,iy)        = st % dv % respo(ic)
        reshe(ix,iy)        = st % dv % reshe(ic)
        reshi(ix,iy)        = st % dv % reshi(ic)
        reshn(ix,iy)        = st % dv % reshn(ic)
!       resht(ix,iy)        = st % dv % resht(ic)
        resmt(ix,iy)        = st % dv % resmt(ic)
        reskt(ix,iy)        = st % dv % reskt(ic)
        reszt(ix,iy)        = st % dv % reszt(ic)

! corrections
        corua(ix,iy,0:ns-1) = st % dv % corua(ic,0:ns-1)
        corpa(ix,iy,0:ns-1) = st % dv % corpa(ic,0:ns-1)
        corut(ix,iy)        = st % dv % corut(ic)
        corpo(ix,iy)        = st % dv % corpo(ic)
!       cortt(ix,iy)        = st % dv % cortt(ic)
        corte(ix,iy)        = st % dv % corte(ic)
        corti(ix,iy)        = st % dv % corti(ic)
        cortn(ix,iy)        = st % dv % cortn(ic)

! coupling coefficients
        pcca(ix,iy,0,0,0:ns-1) = st % dv % pcca(ic,0,0:ns-1)
        pccm(ix,iy,0,0)        = st % dv % pccm(ic,0)
        pcca(ix,iy,1,1,0:ns-1) = st % dv % pcca(ic,1,0:ns-1)
        pccm(ix,iy,1,1)        = st % dv % pccm(ic,1)

! B2SourceWork
        sch(ix,iy,0:3)        = st % sr % sch(ic,0:3)
        she(ix,iy,0:3)        = st % sr % she(ic,0:3)
        shi(ix,iy,0:3)        = st % sr % shi(ic,0:3)
        shn(ix,iy,0:3)        = st % sr % shn(ic,0:3)
        skt(ix,iy,0:3)        = st % sr % skt(ic,0:3)
        skt_prod(ix,iy)       = st % sr % skt_prod(ic)
        skt_diss(ix,iy)       = st % sr % skt_diss(ic)
        szt(ix,iy,0:3)        = st % sr % szt(ic,0:3)
        smo(ix,iy,0:3,0:ns-1) = st % sr % smo(ic,0:3,0:ns-1)
        smq(ix,iy,0:3,0:ns-1) = st % sr % smq(ic,0:3,0:ns-1)
        sna(ix,iy,0:1,0:ns-1) = st % sr % sna(ic,0:1,0:ns-1)

! B2rates
        rlqa(ix,iy,0:1,0:ns-1) = st % rt % rlqa(ic,0:1,0:ns-1)
        rlrd(ix,iy,0:1,0:ns-1) = st % rt % rlrd(ic,0:1,0:ns-1)
        rlqa(ix,iy,0:1,0:ns-1) = st % rt % rlqa(ic,0:1,0:ns-1)
        rlqr(ix,iy,0:1,0:ns-1) = st % rt % rlqr(ic,0:1,0:ns-1)
        rlra(ix,iy,0:1,0:ns-1) = st % rt % rlra(ic,0:1,0:ns-1)
        rlsa(ix,iy,0:1,0:ns-1) = st % rt % rlsa(ic,0:1,0:ns-1)
        rlza(ix,iy,0:1,0:ns-1) = st % rt % rlza(ic,0:1,0:ns-1)
        rlz2(ix,iy,0:1,0:ns-1) = st % rt % rlz2(ic,0:1,0:ns-1)
        rlpt(ix,iy,0:1,0:ns-1) = st % rt % rlpt(ic,0:1,0:ns-1)
        rlpi(ix,iy,0:1,0:ns-1) = st % rt % rlpi(ic,0:1,0:ns-1)

        rza(ix,iy,0:ns-1) = st % rt % rza(ic,0:ns-1)
        rz2(ix,iy,0:ns-1) = st % rt % rz2(ic,0:ns-1)
        rpt(ix,iy,0:ns-1) = st % rt % rpt(ic,0:ns-1)
        rpi(ix,iy,0:ns-1) = st % rt % rpi(ic,0:ns-1)
        endif

        iy = ny
        if(mpg%imapCv(ix,iy)/=0) then
          ic =  mpg%imapCv(ix,iy)

! B2plasma
        na(ix,iy,0:ns-1) = st % pl % na(ic,0:ns-1)
        ua(ix,iy,0:ns-1) = st % pl % ua(ic,0:ns-1)
        po(ix,iy)        = st % pl % po(ic)
        te(ix,iy)        = st % pl % te(ic)
        ti(ix,iy)        = st % pl % ti(ic)
        tn(ix,iy)        = st % pl % tn(ic)
        kt(ix,iy)        = st % pl % kt(ic)
        zt(ix,iy)        = st % pl % zt(ic)

! B2coeff
        ceqp(ix,iy) = st % co % ceqp(ic)

! from b2coeffwork
        dna0(ix,iy,0:ns-1)     = st % co % dna0(ic,0:ns-1)
        dpa0(ix,iy,0:ns-1)     = st % co % dpa0(ic,0:ns-1)
        vla0(ix,iy,0:1,0:ns-1) = st % co % vla0(ic,0:1,0:ns-1)
        vsa0(ix,iy,0:ns-1)     = st % co % vsa0(ic,0:ns-1)
        hci0(ix,iy)            = st % co % hci0(ic)
        hcn0(ix,iy)            = st % co % hcn0(ic)
        hcib(ix,iy,0:ns-1)     = st % co % hcib(ic,0:ns-1)
        hce0(ix,iy)            = st % co % hce0(ic)
        sig0(ix,iy)            = st % co % sig0(ic)
        alf0(ix,iy)            = st % co % alf0(ic)
        dkt0(ix,iy)            = st % co % dkt0(ic)
        dzt0(ix,iy)            = st % co % dzt0(ic)
        dna_ExB(ix,iy)         = st % co % dna_ExB(ic)
        hce_ExB(ix,iy)         = st % co % hce_ExB(ic)
        hci_ExB(ix,iy)         = st % co % hci_ExB(ic)
        cthe(ix,iy,0:ns-1)     = st % co % cthe(ic,0:ns-1)
        cthi(ix,iy,0:ns-1)     = st % co % cthi(ic,0:ns-1)

! kinetic energy
        kinrgy(ix,iy,0:ns-1) = st % dv % kinrgy(ic,0:ns-1)

! derived state quantities
        ne(ix,iy)     = st % dv % ne(ic)
        ne2(ix,iy)    = st % dv % ne2(ic)
        ni(ix,iy,0:1) = st % dv % ni(ic,0:1)
!       ue(ix,iy)     = st % dv % ue(ic)

! residuals
        resmo(ix,iy,0:ns-1) = st % dv % resmo(ic,0:ns-1)
        resco(ix,iy,0:ns-1) = st % dv % resco(ic,0:ns-1)
        respo(ix,iy)        = st % dv % respo(ic)
        reshe(ix,iy)        = st % dv % reshe(ic)
        reshi(ix,iy)        = st % dv % reshi(ic)
        reshn(ix,iy)        = st % dv % reshn(ic)
!       resht(ix,iy)        = st % dv % resht(ic)
        resmt(ix,iy)        = st % dv % resmt(ic)
        reskt(ix,iy)        = st % dv % reskt(ic)
        reszt(ix,iy)        = st % dv % reszt(ic)

! corrections
        corua(ix,iy,0:ns-1) = st % dv % corua(ic,0:ns-1)
        corpa(ix,iy,0:ns-1) = st % dv % corpa(ic,0:ns-1)
        corut(ix,iy)        = st % dv % corut(ic)
        corpo(ix,iy)        = st % dv % corpo(ic)
!       cortt(ix,iy)        = st % dv % cortt(ic)
        corte(ix,iy)        = st % dv % corte(ic)
        corti(ix,iy)        = st % dv % corti(ic)
        cortn(ix,iy)        = st % dv % cortn(ic)

! coupling coefficients
        pcca(ix,iy,0,0,0:ns-1) = st % dv % pcca(ic,0,0:ns-1)
        pccm(ix,iy,0,0)        = st % dv % pccm(ic,0)
        pcca(ix,iy,1,1,0:ns-1) = st % dv % pcca(ic,1,0:ns-1)
        pccm(ix,iy,1,1)        = st % dv % pccm(ic,1)

! B2SourceWork
        sch(ix,iy,0:3)        = st % sr % sch(ic,0:3)
        she(ix,iy,0:3)        = st % sr % she(ic,0:3)
        shi(ix,iy,0:3)        = st % sr % shi(ic,0:3)
        shn(ix,iy,0:3)        = st % sr % shn(ic,0:3)
        skt(ix,iy,0:3)        = st % sr % skt(ic,0:3)
        skt_prod(ix,iy)       = st % sr % skt_prod(ic)
        skt_diss(ix,iy)       = st % sr % skt_diss(ic)
        szt(ix,iy,0:3)        = st % sr % szt(ic,0:3)
        smo(ix,iy,0:3,0:ns-1) = st % sr % smo(ic,0:3,0:ns-1)
        smq(ix,iy,0:3,0:ns-1) = st % sr % smq(ic,0:3,0:ns-1)
        sna(ix,iy,0:1,0:ns-1) = st % sr % sna(ic,0:1,0:ns-1)

! B2rates
        rlqa(ix,iy,0:1,0:ns-1) = st % rt % rlqa(ic,0:1,0:ns-1)
        rlrd(ix,iy,0:1,0:ns-1) = st % rt % rlrd(ic,0:1,0:ns-1)
        rlqa(ix,iy,0:1,0:ns-1) = st % rt % rlqa(ic,0:1,0:ns-1)
        rlqr(ix,iy,0:1,0:ns-1) = st % rt % rlqr(ic,0:1,0:ns-1)
        rlra(ix,iy,0:1,0:ns-1) = st % rt % rlra(ic,0:1,0:ns-1)
        rlsa(ix,iy,0:1,0:ns-1) = st % rt % rlsa(ic,0:1,0:ns-1)
        rlza(ix,iy,0:1,0:ns-1) = st % rt % rlza(ic,0:1,0:ns-1)
        rlz2(ix,iy,0:1,0:ns-1) = st % rt % rlz2(ic,0:1,0:ns-1)
        rlpt(ix,iy,0:1,0:ns-1) = st % rt % rlpt(ic,0:1,0:ns-1)
        rlpi(ix,iy,0:1,0:ns-1) = st % rt % rlpi(ic,0:1,0:ns-1)

        rza(ix,iy,0:ns-1) = st % rt % rza(ic,0:ns-1)
        rz2(ix,iy,0:ns-1) = st % rt % rz2(ic,0:ns-1)
        rpt(ix,iy,0:ns-1) = st % rt % rpt(ic,0:ns-1)
        rpi(ix,iy,0:ns-1) = st % rt % rpi(ic,0:ns-1)
        endif
      end do

! guard cells ix = -1 and ix = nx
      do iy = 0, ny-1
        ix = -1
        if(mpg%imapCv(ix,iy)/=0) then
          ic =  mpg%imapCv(ix,iy)

! B2plasma
        na(ix,iy,0:ns-1) = st % pl % na(ic,0:ns-1)
        ua(ix,iy,0:ns-1) = st % pl % ua(ic,0:ns-1)
        po(ix,iy)        = st % pl % po(ic)
        te(ix,iy)        = st % pl % te(ic)
        ti(ix,iy)        = st % pl % ti(ic)
        tn(ix,iy)        = st % pl % tn(ic)
        kt(ix,iy)        = st % pl % kt(ic)
        zt(ix,iy)        = st % pl % zt(ic)

! B2coeff
        ceqp(ix,iy) = st % co % ceqp(ic)

! from b2coeffwork
        dna0(ix,iy,0:ns-1)     = st % co % dna0(ic,0:ns-1)
        dpa0(ix,iy,0:ns-1)     = st % co % dpa0(ic,0:ns-1)
        vla0(ix,iy,0:1,0:ns-1) = st % co % vla0(ic,0:1,0:ns-1)
        vsa0(ix,iy,0:ns-1)     = st % co % vsa0(ic,0:ns-1)
        hci0(ix,iy)            = st % co % hci0(ic)
        hcn0(ix,iy)            = st % co % hcn0(ic)
        hcib(ix,iy,0:ns-1)     = st % co % hcib(ic,0:ns-1)
        hce0(ix,iy)            = st % co % hce0(ic)
        sig0(ix,iy)            = st % co % sig0(ic)
        alf0(ix,iy)            = st % co % alf0(ic)
        dkt0(ix,iy)            = st % co % dkt0(ic)
        dzt0(ix,iy)            = st % co % dzt0(ic)
        dna_ExB(ix,iy)         = st % co % dna_ExB(ic)
        hce_ExB(ix,iy)         = st % co % hce_ExB(ic)
        hci_ExB(ix,iy)         = st % co % hci_ExB(ic)
        cthe(ix,iy,0:ns-1)     = st % co % cthe(ic,0:ns-1)
        cthi(ix,iy,0:ns-1)     = st % co % cthi(ic,0:ns-1)

! kinetic energy
        kinrgy(ix,iy,0:ns-1) = st % dv % kinrgy(ic,0:ns-1)

! derived state quantities
        ne(ix,iy)     = st % dv % ne(ic)
        ne2(ix,iy)    = st % dv % ne2(ic)
        ni(ix,iy,0:1) = st % dv % ni(ic,0:1)
!       ue(ix,iy)     = st % dv % ue(ic)

! residuals
        resmo(ix,iy,0:ns-1) = st % dv % resmo(ic,0:ns-1)
        resco(ix,iy,0:ns-1) = st % dv % resco(ic,0:ns-1)
        respo(ix,iy)        = st % dv % respo(ic)
        reshe(ix,iy)        = st % dv % reshe(ic)
        reshi(ix,iy)        = st % dv % reshi(ic)
        reshn(ix,iy)        = st % dv % reshn(ic)
!       resht(ix,iy)        = st % dv % resht(ic)
        resmt(ix,iy)        = st % dv % resmt(ic)
        reskt(ix,iy)        = st % dv % reskt(ic)
        reszt(ix,iy)        = st % dv % reszt(ic)

! corrections
        corua(ix,iy,0:ns-1) = st % dv % corua(ic,0:ns-1)
        corpa(ix,iy,0:ns-1) = st % dv % corpa(ic,0:ns-1)
        corut(ix,iy)        = st % dv % corut(ic)
        corpo(ix,iy)        = st % dv % corpo(ic)
!       cortt(ix,iy)        = st % dv % cortt(ic)
        corte(ix,iy)        = st % dv % corte(ic)
        corti(ix,iy)        = st % dv % corti(ic)
        cortn(ix,iy)        = st % dv % cortn(ic)

! coupling coefficients
        pcca(ix,iy,0,0,0:ns-1) = st % dv % pcca(ic,0,0:ns-1)
        pccm(ix,iy,0,0)        = st % dv % pccm(ic,0)
        pcca(ix,iy,1,1,0:ns-1) = st % dv % pcca(ic,1,0:ns-1)
        pccm(ix,iy,1,1)        = st % dv % pccm(ic,1)

! B2SourceWork
        sch(ix,iy,0:3)        = st % sr % sch(ic,0:3)
        she(ix,iy,0:3)        = st % sr % she(ic,0:3)
        shi(ix,iy,0:3)        = st % sr % shi(ic,0:3)
        shn(ix,iy,0:3)        = st % sr % shn(ic,0:3)
        skt(ix,iy,0:3)        = st % sr % skt(ic,0:3)
        skt_prod(ix,iy)       = st % sr % skt_prod(ic)
        skt_diss(ix,iy)       = st % sr % skt_diss(ic)
        szt(ix,iy,0:3)        = st % sr % szt(ic,0:3)
        smo(ix,iy,0:3,0:ns-1) = st % sr % smo(ic,0:3,0:ns-1)
        smq(ix,iy,0:3,0:ns-1) = st % sr % smq(ic,0:3,0:ns-1)
        sna(ix,iy,0:1,0:ns-1) = st % sr % sna(ic,0:1,0:ns-1)

! B2rates
        rlqa(ix,iy,0:1,0:ns-1) = st % rt % rlqa(ic,0:1,0:ns-1)
        rlrd(ix,iy,0:1,0:ns-1) = st % rt % rlrd(ic,0:1,0:ns-1)
        rlqa(ix,iy,0:1,0:ns-1) = st % rt % rlqa(ic,0:1,0:ns-1)
        rlqr(ix,iy,0:1,0:ns-1) = st % rt % rlqr(ic,0:1,0:ns-1)
        rlra(ix,iy,0:1,0:ns-1) = st % rt % rlra(ic,0:1,0:ns-1)
        rlsa(ix,iy,0:1,0:ns-1) = st % rt % rlsa(ic,0:1,0:ns-1)
        rlza(ix,iy,0:1,0:ns-1) = st % rt % rlza(ic,0:1,0:ns-1)
        rlz2(ix,iy,0:1,0:ns-1) = st % rt % rlz2(ic,0:1,0:ns-1)
        rlpt(ix,iy,0:1,0:ns-1) = st % rt % rlpt(ic,0:1,0:ns-1)
        rlpi(ix,iy,0:1,0:ns-1) = st % rt % rlpi(ic,0:1,0:ns-1)

        rza(ix,iy,0:ns-1) = st % rt % rza(ic,0:ns-1)
        rz2(ix,iy,0:ns-1) = st % rt % rz2(ic,0:ns-1)
        rpt(ix,iy,0:ns-1) = st % rt % rpt(ic,0:ns-1)
        rpi(ix,iy,0:ns-1) = st % rt % rpi(ic,0:ns-1)
        endif

        ix = nx
        if(mpg%imapCv(ix,iy)/=0) then
          ic =  mpg%imapCv(ix,iy)

! B2plasma
        na(ix,iy,0:ns-1) = st % pl % na(ic,0:ns-1)
        ua(ix,iy,0:ns-1) = st % pl % ua(ic,0:ns-1)
        po(ix,iy)        = st % pl % po(ic)
        te(ix,iy)        = st % pl % te(ic)
        ti(ix,iy)        = st % pl % ti(ic)
        tn(ix,iy)        = st % pl % tn(ic)
        kt(ix,iy)        = st % pl % kt(ic)
        zt(ix,iy)        = st % pl % zt(ic)

! B2coeff
        ceqp(ix,iy) = st % co % ceqp(ic)

! from b2coeffwork
        dna0(ix,iy,0:ns-1)     = st % co % dna0(ic,0:ns-1)
        dpa0(ix,iy,0:ns-1)     = st % co % dpa0(ic,0:ns-1)
        vla0(ix,iy,0:1,0:ns-1) = st % co % vla0(ic,0:1,0:ns-1)
        vsa0(ix,iy,0:ns-1)     = st % co % vsa0(ic,0:ns-1)
        hci0(ix,iy)            = st % co % hci0(ic)
        hcn0(ix,iy)            = st % co % hcn0(ic)
        hcib(ix,iy,0:ns-1)     = st % co % hcib(ic,0:ns-1)
        hce0(ix,iy)            = st % co % hce0(ic)
        sig0(ix,iy)            = st % co % sig0(ic)
        alf0(ix,iy)            = st % co % alf0(ic)
        dkt0(ix,iy)            = st % co % dkt0(ic)
        dzt0(ix,iy)            = st % co % dzt0(ic)
        dna_ExB(ix,iy)         = st % co % dna_ExB(ic)
        hce_ExB(ix,iy)         = st % co % hce_ExB(ic)
        hci_ExB(ix,iy)         = st % co % hci_ExB(ic)
        cthe(ix,iy,0:ns-1)     = st % co % cthe(ic,0:ns-1)
        cthi(ix,iy,0:ns-1)     = st % co % cthi(ic,0:ns-1)

! kinetic energy
        kinrgy(ix,iy,0:ns-1) = st % dv % kinrgy(ic,0:ns-1)

! derived state quantities
        ne(ix,iy)     = st % dv % ne(ic)
        ne2(ix,iy)    = st % dv % ne2(ic)
        ni(ix,iy,0:1) = st % dv % ni(ic,0:1)
!       ue(ix,iy)     = st % dv % ue(ic)

! residuals
        resmo(ix,iy,0:ns-1) = st % dv % resmo(ic,0:ns-1)
        resco(ix,iy,0:ns-1) = st % dv % resco(ic,0:ns-1)
        respo(ix,iy)        = st % dv % respo(ic)
        reshe(ix,iy)        = st % dv % reshe(ic)
        reshi(ix,iy)        = st % dv % reshi(ic)
        reshn(ix,iy)        = st % dv % reshn(ic)
!       resht(ix,iy)        = st % dv % resht(ic)
        resmt(ix,iy)        = st % dv % resmt(ic)
        reskt(ix,iy)        = st % dv % reskt(ic)
        reszt(ix,iy)        = st % dv % reszt(ic)

! corrections
        corua(ix,iy,0:ns-1) = st % dv % corua(ic,0:ns-1)
        corpa(ix,iy,0:ns-1) = st % dv % corpa(ic,0:ns-1)
        corut(ix,iy)        = st % dv % corut(ic)
        corpo(ix,iy)        = st % dv % corpo(ic)
!       cortt(ix,iy)        = st % dv % cortt(ic)
        corte(ix,iy)        = st % dv % corte(ic)
        corti(ix,iy)        = st % dv % corti(ic)
        cortn(ix,iy)        = st % dv % cortn(ic)

! coupling coefficients
        pcca(ix,iy,0,0,0:ns-1) = st % dv % pcca(ic,0,0:ns-1)
        pccm(ix,iy,0,0)        = st % dv % pccm(ic,0)
        pcca(ix,iy,1,1,0:ns-1) = st % dv % pcca(ic,1,0:ns-1)
        pccm(ix,iy,1,1)        = st % dv % pccm(ic,1)

! B2SourceWork
        sch(ix,iy,0:3)        = st % sr % sch(ic,0:3)
        she(ix,iy,0:3)        = st % sr % she(ic,0:3)
        shi(ix,iy,0:3)        = st % sr % shi(ic,0:3)
        shn(ix,iy,0:3)        = st % sr % shn(ic,0:3)
        skt(ix,iy,0:3)        = st % sr % skt(ic,0:3)
        skt_prod(ix,iy)       = st % sr % skt_prod(ic)
        skt_diss(ix,iy)       = st % sr % skt_diss(ic)
        szt(ix,iy,0:3)        = st % sr % szt(ic,0:3)
        smo(ix,iy,0:3,0:ns-1) = st % sr % smo(ic,0:3,0:ns-1)
        smq(ix,iy,0:3,0:ns-1) = st % sr % smq(ic,0:3,0:ns-1)
        sna(ix,iy,0:1,0:ns-1) = st % sr % sna(ic,0:1,0:ns-1)

! B2rates
        rlqa(ix,iy,0:1,0:ns-1) = st % rt % rlqa(ic,0:1,0:ns-1)
        rlrd(ix,iy,0:1,0:ns-1) = st % rt % rlrd(ic,0:1,0:ns-1)
        rlqa(ix,iy,0:1,0:ns-1) = st % rt % rlqa(ic,0:1,0:ns-1)
        rlqr(ix,iy,0:1,0:ns-1) = st % rt % rlqr(ic,0:1,0:ns-1)
        rlra(ix,iy,0:1,0:ns-1) = st % rt % rlra(ic,0:1,0:ns-1)
        rlsa(ix,iy,0:1,0:ns-1) = st % rt % rlsa(ic,0:1,0:ns-1)
        rlza(ix,iy,0:1,0:ns-1) = st % rt % rlza(ic,0:1,0:ns-1)
        rlz2(ix,iy,0:1,0:ns-1) = st % rt % rlz2(ic,0:1,0:ns-1)
        rlpt(ix,iy,0:1,0:ns-1) = st % rt % rlpt(ic,0:1,0:ns-1)
        rlpi(ix,iy,0:1,0:ns-1) = st % rt % rlpi(ic,0:1,0:ns-1)

        rza(ix,iy,0:ns-1) = st % rt % rza(ic,0:ns-1)
        rz2(ix,iy,0:ns-1) = st % rt % rz2(ic,0:ns-1)
        rpt(ix,iy,0:ns-1) = st % rt % rpt(ic,0:ns-1)
        rpi(ix,iy,0:ns-1) = st % rt % rpi(ic,0:ns-1)
        endif
      end do

#ifdef WG_TODO
! Cannot do this anymore => b2mod_indirect not available in general
! Need a back-transformation to_struct_geometry first !?!
! Corner cells: add-hoc interpolation

!!!! CTB: I think this needs still to be adapted for Double null?

        iy=-1                                                         !srv 27.06.08 {
        do ix = -1, nx
          if(leftix(ix,iy).eq.-2 ) then
           na(ix,iy,:) = exp(
     .      log(na(rightix(ix,iy),rightiy(ix,iy),:)) +
     .      log(na(topix(ix,iy),topiy(ix,iy),:)) -
     &      log(na(topix(rightix(ix,iy),rightiy(ix,iy)),
     &             topiy(rightix(ix,iy),rightiy(ix,iy)),:)))
           ua(ix,iy,:) =
     .      ua(rightix(ix,iy),rightiy(ix,iy),:) +
     .      ua(topix(ix,iy),topiy(ix,iy),:) -
     &      ua(topix(rightix(ix,iy),rightiy(ix,iy)),
     &             topiy(rightix(ix,iy),rightiy(ix,iy)),:)
           te(ix,iy) = exp(
     .      log(te(rightix(ix,iy),rightiy(ix,iy))) +
     .      log(te(topix(ix,iy),topiy(ix,iy))) -
     &      log(te(topix(rightix(ix,iy),rightiy(ix,iy)),
     &             topiy(rightix(ix,iy),rightiy(ix,iy)))))
           ti(ix,iy) = exp(
     .      log(ti(rightix(ix,iy),rightiy(ix,iy))) +
     .      log(ti(topix(ix,iy),topiy(ix,iy))) -
     &      log(ti(topix(rightix(ix,iy),rightiy(ix,iy)),
     &             topiy(rightix(ix,iy),rightiy(ix,iy)))))
           po(ix,iy) =
     .      po(rightix(ix,iy),rightiy(ix,iy)) +
     .      po(topix(ix,iy),topiy(ix,iy)) -
     &      po(topix(rightix(ix,iy),rightiy(ix,iy)),
     &             topiy(rightix(ix,iy),rightiy(ix,iy)))
          endif
          if(rightix(ix,iy).eq.nx+1 ) then
           na(ix,iy,:) = exp(
     .      log(na(leftix(ix,iy),leftiy(ix,iy),:)) +
     .      log(na(topix(ix,iy),topiy(ix,iy),:)) -
     &      log(na(topix(leftix(ix,iy),leftiy(ix,iy)),
     &             topiy(leftix(ix,iy),leftiy(ix,iy)),:)))
           ua(ix,iy,:) =
     .      ua(leftix(ix,iy),leftiy(ix,iy),:) +
     .      ua(topix(ix,iy),topiy(ix,iy),:) -
     &      ua(topix(leftix(ix,iy),leftiy(ix,iy)),
     &             topiy(leftix(ix,iy),leftiy(ix,iy)),:)
           te(ix,iy) = exp(
     .      log(te(leftix(ix,iy),leftiy(ix,iy))) +
     .      log(te(topix(ix,iy),topiy(ix,iy))) -
     &      log(te(topix(leftix(ix,iy),leftiy(ix,iy)),
     &             topiy(leftix(ix,iy),leftiy(ix,iy)))))
           ti(ix,iy) = exp(
     .      log(ti(leftix(ix,iy),leftiy(ix,iy))) +
     .      log(ti(topix(ix,iy),topiy(ix,iy))) -
     &      log(ti(topix(leftix(ix,iy),leftiy(ix,iy)),
     &             topiy(leftix(ix,iy),leftiy(ix,iy)))))
           po(ix,iy) =
     .      po(leftix(ix,iy),leftiy(ix,iy)) +
     .      po(topix(ix,iy),topiy(ix,iy)) -
     &      po(topix(leftix(ix,iy),leftiy(ix,iy)),
     &             topiy(leftix(ix,iy),leftiy(ix,iy)))
          endif
        enddo
        iy=ny
        do ix = -1, nx
          if(leftix(ix,iy).eq.-2 ) then
           na(ix,iy,:) = exp(
     .      log(na(rightix(ix,iy),rightiy(ix,iy),:)) +
     .      log(na(bottomix(ix,iy),bottomiy(ix,iy),:)) -
     &      log(na(bottomix(rightix(ix,iy),rightiy(ix,iy)),
     &             bottomiy(rightix(ix,iy),rightiy(ix,iy)),:)))
           ua(ix,iy,:) =
     .      ua(rightix(ix,iy),rightiy(ix,iy),:) +
     .      ua(bottomix(ix,iy),bottomiy(ix,iy),:) -
     &      ua(bottomix(rightix(ix,iy),rightiy(ix,iy)),
     &             bottomiy(rightix(ix,iy),rightiy(ix,iy)),:)
           te(ix,iy) = exp(
     .      log(te(rightix(ix,iy),rightiy(ix,iy))) +
     .      log(te(bottomix(ix,iy),bottomiy(ix,iy))) -
     &      log(te(bottomix(rightix(ix,iy),rightiy(ix,iy)),
     &             bottomiy(rightix(ix,iy),rightiy(ix,iy)))))
           ti(ix,iy) = exp(
     .      log(ti(rightix(ix,iy),rightiy(ix,iy))) +
     .      log(ti(bottomix(ix,iy),bottomiy(ix,iy))) -
     &      log(ti(bottomix(rightix(ix,iy),rightiy(ix,iy)),
     &             bottomiy(rightix(ix,iy),rightiy(ix,iy)))))
           po(ix,iy) =
     .      po(rightix(ix,iy),rightiy(ix,iy)) +
     .      po(bottomix(ix,iy),bottomiy(ix,iy)) -
     &      po(bottomix(rightix(ix,iy),rightiy(ix,iy)),
     &             bottomiy(rightix(ix,iy),rightiy(ix,iy)))
          endif
          if(rightix(ix,iy).eq.nx+1 ) then
           na(ix,iy,:) = exp(
     .      log(na(leftix(ix,iy),leftiy(ix,iy),:)) +
     .      log(na(bottomix(ix,iy),bottomiy(ix,iy),:)) -
     &      log(na(bottomix(leftix(ix,iy),leftiy(ix,iy)),
     &             bottomiy(leftix(ix,iy),leftiy(ix,iy)),:)))
           ua(ix,iy,:) =
     .      ua(leftix(ix,iy),leftiy(ix,iy),:) +
     .      ua(bottomix(ix,iy),bottomiy(ix,iy),:) -
     &      ua(bottomix(leftix(ix,iy),leftiy(ix,iy)),
     &             bottomiy(leftix(ix,iy),leftiy(ix,iy)),:)
           te(ix,iy) = exp(
     .      log(te(leftix(ix,iy),leftiy(ix,iy))) +
     .      log(te(bottomix(ix,iy),bottomiy(ix,iy))) -
     &      log(te(bottomix(leftix(ix,iy),leftiy(ix,iy)),
     &             bottomiy(leftix(ix,iy),leftiy(ix,iy)))))
           ti(ix,iy) = exp(
     .      log(ti(leftix(ix,iy),leftiy(ix,iy))) +
     .      log(ti(bottomix(ix,iy),bottomiy(ix,iy))) -
     &      log(ti(bottomix(leftix(ix,iy),leftiy(ix,iy)),
     &             bottomiy(leftix(ix,iy),leftiy(ix,iy)))))
           po(ix,iy) =
     .      po(leftix(ix,iy),leftiy(ix,iy)) +
     .      po(bottomix(ix,iy),bottomiy(ix,iy)) -
     &      po(bottomix(leftix(ix,iy),leftiy(ix,iy)),
     &             bottomiy(leftix(ix,iy),leftiy(ix,iy)))
          endif
        enddo                                                         !srv 27.06.08 }
#endif

! x faces
      do iy = 0, ny-1
        do ix = 0, nx
          if(mpg%imapFcx(ix,iy)/=0) then
            ifc = mpg%imapFcx(ix,iy)

! B2coeff
          csig_an(ix,iy,0,0:1) = st % co % csig_an(ifc,0:1)
          calf(ix,iy,0,0:1)    = st % co % calf(ifc,0:1)
          cdkt(ix,iy,0,0:1)    = st % co % cdkt(ifc,0:1)
          cdzt(ix,iy,0,0:1)    = st % co % cdzt(ifc,0:1)

          csig(ix,iy,0,0:1) = st % co % csig(ifc,0:1)
          csigin(ix,iy,0,0:1,0:ns-1,0:ns-1) =
     .           st % co % csigin(ifc,0:1,0:ns-1,0:ns-1)

          chve(ix,iy,0,0:1) = st % co % chve(ifc,0:1)
          chce(ix,iy,0,0:1) = st % co % chce(ifc,0:1)
          chvemx(ix,iy)     = st % co % chvemx(ifc)

          chvi(ix,iy,0,0:1) = st % co % chvi(ifc,0:1)
          chci(ix,iy,0,0:1) = st % co % chci(ifc,0:1)
          chvimx(ix,iy)     = st % co % chvimx(ifc)

          cvla(ix,iy,0,0:1,0:ns-1) = st % co % cvla(ifc,0:1,0:ns-1)
          cdna(ix,iy,0,0:1,0:ns-1) = st % co % cdna(ifc,0:1,0:ns-1)
          cdpa(ix,iy,0,0:1,0:ns-1) = st % co % cdpa(ifc,0:1,0:ns-1)

          cvsa(ix,iy,0,0:1,0:ns-1) = st % co % cvsa(ifc,0:1,0:ns-1)
          cvsa_cl(ix,iy,0,0:1,0:ns-1) =
     .            st % co % cvsa_cl(ifc,0:1,0:ns-1)

          cddi(ix,iy,0,0:1,0:ns-1) = st % co % cddi(ifc,0:1,0:ns-1)

          fllim0fna(ix,iy,0,0:1,0:ns-1) =
     .              st % co % fllim0fna(ifc,0:1,0:ns-1)
          fllim0fhi(ix,iy,0,0:1,0:ns-1) =
     .              st % co % fllim0fhi(ifc,0:1,0:ns-1)
          fllimvisc(ix,iy,0:ns-1) = st % co % fllimvisc(ifc,0:ns-1)

! fluxes
          fch(ix,iy,0,0:1)       = st % dv % fch(ifc,0:1)
          fch_32(ix,iy,0,0:1)    = st % dv % fch_32(ifc,0:1)
          fch_52(ix,iy,0,0:1)    = st % dv % fch_52(ifc,0:1)
          fch_p(ix,iy,0,0:1)     = st % dv % fch_p(ifc,0:1)
          fchdia(ix,iy,0,0:1)    = st % dv % fchdia(ifc,0:1)
          fchin(ix,iy,0,0:1)     = st % dv % fchin(ifc,0:1)
          fchvispar(ix,iy,0,0:1) = st % dv % fchvispar(ifc,0:1)
          fchvisper(ix,iy,0,0:1) = st % dv % fchvisper(ifc,0:1)
          fchvisq(ix,iy,0,0:1)   = st % dv % fchvisq(ifc,0:1)
          fchviskt(ix,iy,0,0:1)  = st % dv % fchviskt(ifc,0:1)
          fchinert(ix,iy,0,0:1)  = st % dv % fchinert(ifc,0:1)
          fchanml(ix,iy,0,0:1)   = st % dv % fchanml(ifc,0:1)

          fne(ix,iy,0,0:1)    = st % dv % fne(ifc,0:1)
          fne_32(ix,iy,0,0:1) = st % dv % fne_32(ifc,0:1)
          fne_52(ix,iy,0,0:1) = st % dv % fne_52(ifc,0:1)

          fni(ix,iy,0,0:1)    = st % dv % fni(ifc,0:1)
          fni_32(ix,iy,0,0:1) = st % dv % fni_32(ifc,0:1)
          fni_52(ix,iy,0,0:1) = st % dv % fni_52(ifc,0:1)

          fna(ix,iy,0,0:1,0:ns-1)     = st % dv % fna(ifc,0:1,0:ns-1)
          fna_mdf(ix,iy,0,0:1,0:ns-1) =
     .            st % dv % fna_mdf(ifc,0:1,0:ns-1)
          fna_he(ix,iy,0,0:1,0:ns-1)  = st % dv % fna_he(ifc,0:1,0:ns-1)
          fna_fcor(ix,iy,0,0:1,0:ns-1) =
     .            st % dv % fna_fcor(ifc,0:1,0:ns-1)
          fnaPSch(ix,iy,0,0:1,0:ns-1) =
     .            st % dv % fnaPSch(ifc,0:1,0:ns-1)
          fna_32(ix,iy,0,0:1,0:ns-1)  = st % dv % fna_32(ifc,0:1,0:ns-1)
          fna_52(ix,iy,0,0:1,0:ns-1)  = st % dv % fna_52(ifc,0:1,0:ns-1)
          fna_32nd(ix,iy,0,0:1,0:ns-1) =
     .           st % dv % fna_32nd(ifc,0:1,0:ns-1)
          fna_52nd(ix,iy,0,0:1,0:ns-1) =
     .           st % dv % fna_52nd(ifc,0:1,0:ns-1)
          fna_eir(ix,iy,0,0:1,0:ns-1) =
     .            st % dv % fna_eir(ifc,0:1,0:ns-1)
          fne_eir(ix,iy,0,0:1) =
     .            st % dv % fne_eir(ifc,0:1)
          fna_nodrift(ix,iy,0,0:1,0:ns-1) =
     .           st % dv % fna_nodrift(ifc,0:1,0:ns-1)

          fmo(ix,iy,0,0:1,0:ns-1)     = st % dv % fmo(ifc,0:1,0:ns-1)

          fhe(ix,iy,0,0:1)     = st % dv % fhe(ifc,0:1)
          fkt(ix,iy,0,0:1)     = st % dv % fkt(ifc,0:1)
          fzt(ix,iy,0,0:1)     = st % dv % fzt(ifc,0:1)
          fhe_mdf(ix,iy,0,0:1) = st % dv % fhe_mdf(ifc,0:1)
          fhePSch(ix,iy,0,0:1) = st % dv % fhePSch(ifc,0:1)
          fhet(ix,iy,0,0:1)    = st % dv % fhet(ifc,0:1)
          fhe_eir(ix,iy,0,0:1) = st % dv % fhe_eir(ifc,0:1)
          fhi(ix,iy,0,0:1)     = st % dv % fhi(ifc,0:1)
          fhn(ix,iy,0,0:1)     = st % dv % fhn(ifc,0:1)
          fhi_mdf(ix,iy,0,0:1) = st % dv % fhi_mdf(ifc,0:1)
          fhiPSch(ix,iy,0,0:1) = st % dv % fhiPSch(ifc,0:1)
          fhit(ix,iy,0,0:1)    = st % dv % fhit(ifc,0:1)
          fhi_eir(ix,iy,0,0:1) = st % dv % fhi_eir(ifc,0:1)

! derived quantities
          uadia(ix,iy,0,0:1,0:ns-1) =
     .          st % dv % uadia(ifc,0:1,0:ns-1)
          vadia(ix,iy,0,0:1,0:ns-1) =
     .          st % dv % vadia(ifc,0:1,0:ns-1)
          wadia(ix,iy,0,0:1,0:ns-1) =
     .          st % dv % wadia(ifc,0:1,0:ns-1)
          vaecrb(ix,iy,0,0:1,0:ns-1) =
     .          st % dv % vaecrb(ifc,0:1,0:ns-1)
          wedia(ix,iy,0,0:1) =
     .          st % dv % wedia(ifc,0:1)
          vedia(ix,iy,0,0:1) =
     .          st % dv % vedia(ifc,0:1)
          veecrb(ix,iy,0,0:1) =
     .          st % dv % veecrb(ifc,0:1)
! note that facdrift, fac_ExB, fac_vis used to be defined as cell center quantities,
! but actually used as face quantities in case of drift_flows = 1
          facdrift(ix,iy) = st % dv % facdrift(ifc)
          fac_ExB(ix,iy)  = st % dv % fac_ExB(ifc)
          fac_vis(ix,iy)  = st % dv % fac_vis(ifc)

! numerical coefficients
          conc(ix,iy,0,0:1)   = st % dv % conc(ifc,0:1)
          fllime(ix,iy)       = st % dv % fllime(ifc)
          floe(ix,iy,0,0:1)   = st % dv % floe(ifc,0:1)
          cone(ix,iy,0,0:1,:) = st % dv % cone(ifc,0:1,:)
          fllimi(ix,iy)       = st % dv % fllimi(ifc)
          floi(ix,iy,0,0:1)   = st % dv % floi(ifc,0:1)
          coni(ix,iy,0,0:1,:) = st % dv % coni(ifc,0:1,:)
          floe_noc(ix,iy,0,0:1) = st % dv % floe_noc(ifc,0:1)
          floi_noc(ix,iy,0,0:1) = st % dv % floi_noc(ifc,0:1)
          endif
        end do
      end do

! y faces
      do ix = 0, nx-1
        do iy = 0, ny
          if(mpg%imapFcy(ix,iy)/=0) then
            ifc = mpg%imapFcy(ix,iy)

! B2coeff
          csig_an(ix,iy,1,0:1) = st % co % csig_an(ifc,0:1)
          calf(ix,iy,1,0:1)    = st % co % calf(ifc,0:1)
          cdkt(ix,iy,1,0:1)    = st % co % cdkt(ifc,0:1)
          cdzt(ix,iy,1,0:1)    = st % co % cdzt(ifc,0:1)

          csig(ix,iy,1,0:1) = st % co % csig(ifc,0:1)
          csigin(ix,iy,1,0:1,0:ns-1,0:ns-1) =
     .           st % co % csigin(ifc,0:1,0:ns-1,0:ns-1)

          chve(ix,iy,1,0:1) = st % co % chve(ifc,0:1)
          chce(ix,iy,1,0:1) = st % co % chce(ifc,0:1)
          chvemx(ix,iy)     = st % co % chvemx(ifc)

          chvi(ix,iy,1,0:1) = st % co % chvi(ifc,0:1)
          chci(ix,iy,1,0:1) = st % co % chci(ifc,0:1)
          chvimx(ix,iy)     = st % co % chvimx(ifc)

          cvla(ix,iy,1,0:1,0:ns-1) = st % co % cvla(ifc,0:1,0:ns-1)
          cdna(ix,iy,1,0:1,0:ns-1) = st % co % cdna(ifc,0:1,0:ns-1)
          cdpa(ix,iy,1,0:1,0:ns-1) = st % co % cdpa(ifc,0:1,0:ns-1)

          cvsa(ix,iy,1,0:1,0:ns-1)    = st % co % cvsa(ifc,0:1,0:ns-1)
          cvsa_cl(ix,iy,1,0:1,0:ns-1) =
     .            st % co % cvsa_cl(ifc,0:1,0:ns-1)

          cddi(ix,iy,1,0:1,0:ns-1) = st % co % cddi(ifc,0:1,0:ns-1)

          fllim0fhi(ix,iy,1,0:1,0:ns-1) =
     .              st % co % fllim0fhi(ifc,0:1,0:ns-1)
          fllimvisc(ix,iy,0:ns-1) = st % co % fllimvisc(ifc,0:ns-1)

! fluxes
          fch(ix,iy,1,0:1)       = st % dv % fch(ifc,0:1)
          fch_32(ix,iy,1,0:1)    = st % dv % fch_32(ifc,0:1)
          fch_52(ix,iy,1,0:1)    = st % dv % fch_52(ifc,0:1)
          fch_p(ix,iy,1,0:1)     = st % dv % fch_p(ifc,0:1)
          fchdia(ix,iy,1,0:1)    = st % dv % fchdia(ifc,0:1)
          fchin(ix,iy,1,0:1)     = st % dv % fchin(ifc,0:1)
          fchvispar(ix,iy,1,0:1) = st % dv % fchvispar(ifc,0:1)
          fchvisper(ix,iy,1,0:1) = st % dv % fchvisper(ifc,0:1)
          fchvisq(ix,iy,1,0:1)   = st % dv % fchvisq(ifc,0:1)
          fchviskt(ix,iy,1,0:1)  = st % dv % fchviskt(ifc,0:1)
          fchinert(ix,iy,1,0:1)  = st % dv % fchinert(ifc,0:1)
          fchanml(ix,iy,1,0:1)   = st % dv % fchanml(ifc,0:1)

          fne(ix,iy,1,0:1)    = st % dv % fne(ifc,0:1)
          fne_32(ix,iy,1,0:1) = st % dv % fne_32(ifc,0:1)
          fne_52(ix,iy,1,0:1) = st % dv % fne_52(ifc,0:1)

          fni(ix,iy,1,0:1)    = st % dv % fni(ifc,0:1)
          fni_32(ix,iy,1,0:1) = st % dv % fni_32(ifc,0:1)
          fni_52(ix,iy,1,0:1) = st % dv % fni_52(ifc,0:1)

          fna(ix,iy,1,0:1,0:ns-1) = st % dv % fna(ifc,0:1,0:ns-1)
          fna_mdf(ix,iy,1,0:1,0:ns-1) =
     .            st % dv % fna_mdf(ifc,0:1,0:ns-1)
          fna_he(ix,iy,1,0:1,0:ns-1) = st % dv % fna_he(ifc,0:1,0:ns-1)
          fna_fcor(ix,iy,1,0:1,0:ns-1) =
     .            st % dv % fna_fcor(ifc,0:1,0:ns-1)
          fnaPSch(ix,iy,1,0:1,0:ns-1) =
     .            st % dv % fnaPSch(ifc,0:1,0:ns-1)
          fna_32(ix,iy,1,0:1,0:ns-1) = st % dv % fna_32(ifc,0:1,0:ns-1)
          fna_52(ix,iy,1,0:1,0:ns-1) = st % dv % fna_52(ifc,0:1,0:ns-1)
          fna_32nd(ix,iy,1,0:1,0:ns-1) =
     .            st % dv % fna_32nd(ifc,0:1,0:ns-1)
          fna_52nd(ix,iy,1,0:1,0:ns-1) =
     .            st % dv % fna_52nd(ifc,0:1,0:ns-1)
          fna_eir(ix,iy,1,0:1,0:ns-1) =
     .            st % dv % fna_eir(ifc,0:1,0:ns-1)
          fne_eir(ix,iy,1,0:1) =
     .            st % dv % fne_eir(ifc,0:1)
          fna_nodrift(ix,iy,1,0:1,0:ns-1) =
     .            st % dv % fna_nodrift(ifc,0:1,0:ns-1)

          fmo(ix,iy,1,0:1,0:ns-1) = st % dv % fmo(ifc,0:1,0:ns-1)

          fhe(ix,iy,1,0:1)     = st % dv % fhe(ifc,0:1)
          fkt(ix,iy,1,0:1)     = st % dv % fkt(ifc,0:1)
          fzt(ix,iy,1,0:1)     = st % dv % fzt(ifc,0:1)
          fhe_mdf(ix,iy,1,0:1) = st % dv % fhe_mdf(ifc,0:1)
          fhePSch(ix,iy,1,0:1) = st % dv % fhePSch(ifc,0:1)
          fhet(ix,iy,1,0:1)    = st % dv % fhet(ifc,0:1)
          fhe_eir(ix,iy,1,0:1) = st % dv % fhe_eir(ifc,0:1)
          fhi(ix,iy,1,0:1)     = st % dv % fhi(ifc,0:1)
          fhn(ix,iy,1,0:1)     = st % dv % fhn(ifc,0:1)
          fhi_mdf(ix,iy,1,0:1) = st % dv % fhi_mdf(ifc,0:1)
          fhiPSch(ix,iy,1,0:1) = st % dv % fhiPSch(ifc,0:1)
          fhit(ix,iy,1,0:1)    = st % dv % fhit(ifc,0:1)
          fhi_eir(ix,iy,1,0:1) = st % dv % fhi_eir(ifc,0:1)

! derived quantities
          uadia(ix,iy,1,0:1,0:ns-1) =
     .          st % dv % uadia(ifc,0:1,0:ns-1)
          vadia(ix,iy,1,0:1,0:ns-1) =
     .          st % dv % vadia(ifc,0:1,0:ns-1)
          wadia(ix,iy,1,0:1,0:ns-1) =
     .          st % dv % wadia(ifc,0:1,0:ns-1)
          vaecrb(ix,iy,1,0:1,0:ns-1) =
     .          st % dv % vaecrb(ifc,0:1,0:ns-1)
          wedia(ix,iy,1,0:1) =
     .          st % dv % wedia(ifc,0:1)
          vedia(ix,iy,1,0:1) =
     .          st % dv % vedia(ifc,0:1)
          veecrb(ix,iy,1,0:1) =
     .          st % dv % veecrb(ifc,0:1)
! note that facdrift, fac_ExB, fac_vis used to be defined as cell center quantities,
! but actually used as face quantities in case of drift_flows = 1
          facdrift(ix,iy) = st % dv % facdrift(ifc)
          fac_ExB(ix,iy)  = st % dv % fac_ExB(ifc)
          fac_vis(ix,iy)  = st % dv % fac_vis(ifc)

! numerical coefficients
          conc(ix,iy,1,0:1)   = st % dv % conc(ifc,0:1)
          fllime(ix,iy)       = st % dv % fllime(ifc)
          floe(ix,iy,1,0:1)   = st % dv % floe(ifc,0:1)
          cone(ix,iy,1,0:1,:) = st % dv % cone(ifc,0:1,:)
          fllimi(ix,iy)       = st % dv % fllimi(ifc)
          floi(ix,iy,1,0:1)   = st % dv % floi(ifc,0:1)
          coni(ix,iy,1,0:1,:) = st % dv % coni(ifc,0:1,:)
          floe_noc(ix,iy,1,0:1) = st % dv % floe_noc(ifc,0:1)
          floi_noc(ix,iy,1,0:1) = st % dv % floi_noc(ifc,0:1)
          endif
        end do
      end do

!      aresco(0:ns-1) = st % diag % aresco(0:ns-1)
!      aresmo(0:ns-1) = st % diag % aresmo(0:ns-1)

!      rescoreg(0:nnreg(0),0:ns-1) =
!     .         st % diag % rescoreg(0:nnreg(0),0:ns-1)
!      resmoreg(0:nnreg(0),0:ns-1)
!     .         st % diag % resmoreg(0:nnreg(0),0:ns-1)
!      reshereg(0:nnreg(0)) = st % diag % reshereg(0:nnreg(0))
!      reshireg(0:nnreg(0)) = st % diag % reshireg(0:nnreg(0))

*   ..write b2fstate_st in structured format
#ifdef F2003
      open (newunit=nout_st,file='b2fstate_st',iostat=io)
#else
      nout_st = newunit()
      open (unit=nout_st,file='b2fstate_st',iostat=io)
#endif
      if (io /= 0) then
        write (*,*) ' File b2fstate_st not writable'
        call xerrab ('to_struct_plasma -- no state file ')
      end if
      lblst = 'Plasma state in structured format'

*   ..write file header
      call cfverw (nout_st, newversion)
      idum(0) = nx
      idum(1) = ny
      idum(2) = ns
      call cfwuin (nout_st, 3, idum, 'nx,ny,ns')
      call cfwuch (nout_st, 120, lblst, 'label')
      call b2wuzd (nout_st, newversion, ns, zamin, zamax, zn, am)

*   ..write plasma state
      call b2wups (nout_st, nx, ny, ns, ne, na, ua, uadia, te,
     &  ti, tn, po, kt, zt, fna, fhe, fhi, fhn, fch, fch_32, fch_52,
     &  fch_p, kinrgy,
     &  fna_mdf, fhe_mdf, fhi_mdf, fna_fcor, fna_nodrift, fna_he,
     &  fnaPSch, fhePSch, fhiPSch, fna_eir, fne_eir, fhe_eir, fhi_eir,
     &  fna_32, fna_52, fni_32, fni_52, fne_32, fne_52,
     &  fchdia, fchin, fchvispar, fchvisper, fchvisq, fchinert,
     &  vaecrb, vadia, wadia, veecrb, vedia, floe_noc, floi_noc)

      close(nout_st)

*   ..write b2fplasmf_st in structured format
#ifdef F2003
      open (newunit=nout_st,file='b2fplasmf_st',iostat=io)
#else
      nout_st = newunit()
      open (unit=nout_st,file='b2fplasmf_st',iostat=io)
#endif
      if (io /= 0) then
        write (*,*) ' File b2fplasmf_st not writable'
        call xerrab ('to_struct_plasma -- no plasma file ')
      end if
      lblps = 'Plasma state in structured format'

      call cfverw(nout_st, newversion)
!WG_TODO      call write_b2mod_geo(nout_st)
      call write_b2mod_plasma(nout_st)
      call write_b2mod_residuals(nout_st)
      call write_b2mod_sources(nout_st)
      call write_b2mod_transport(nout_st)
      call write_b2mod_anomalous_transport(nout_st)
      call write_b2mod_neutr_src_scaling(nout_st,nstrat)

      close(nout_st)

*   ..deallocate modules for structured grid
      call dealloc_b2mod_plasma
      call dealloc_b2mod_transport
      call dealloc_anomalous_transport
      call dealloc_b2mod_work
      call dealloc_b2mod_residuals
      call dealloc_b2mod_rates
      call dealloc_b2mod_sources

      return
      end subroutine to_struct_plasma

******************************************************************************

      subroutine prep_plasma_ext (nx, ny, ns, nCv, nCi, nFc, st_ext)

        implicit none
        integer, intent(in) :: nx, ny, ns, nCi, nCv, nFc
        type(B2StateExt), intent(inout) :: st_ext

        ! For now, only create a dummy external state.
        ! Functionality from b2mod_external to be carried over!
        call createB2StateExt( nCv, nFc, ns, st_ext )

      end subroutine prep_plasma_ext

******************************************************************************

      subroutine conv_regions (nnreg)
        implicit none
        integer, intent(in) :: nnreg(0:2)
        integer :: ii, iout
#ifndef F2003
        integer newunit
        external newunit
#endif

! csc Provide conversion of x/y directed regions
#ifdef F2003
        open(newunit=iout,file='b2us.regions.parameters')
#else
        iout = newunit()
        open(unit=iout,file='b2us.regions.parameters')
#endif
        write (iout,'(1x,a,i0)') 'number of face regions  = ',
     &   nnreg(1)+nnreg(2)

        write (iout,'(1x,a,1x,a)') 'x-directed regions','| fcReg'
        do ii = 1, nnreg(1)
          write (iout,'(1x,i0,1x,a,1x,i0)') ii, '|',ii
        end do

        write (iout,'(1x,a,1x,a)') 'y-directed regions','| fcReg'
        do ii = 1, nnreg(2)
          write (iout,'(1x,i0,1x,a,1x,i0)') ii,'|', ii+nnreg(1)
        end do

        close (iout)

      return
      end subroutine conv_regions

******************************************************************************

      subroutine conv_midplane(nx,ny)                                   !srv 16.11.22 {
        integer :: nx, ny
        real (kind=R8) :: ro(2), zo(2), rzomp(2,2)
        real (kind=R8) :: ri(2), zi(2), rzimp(2,2)
        integer :: ix, iy, iymin, iymax, iout, n_omp, n_imp
#ifndef F2003
        integer newunit
        external newunit
#endif

        ! This subroutine calculates arrays rzomp and rzimp with the R-Z coordinates
        ! of two points laying on the straight line going through all
        ! radial outer and inner midlane cells
        ! Writes the coordinates into b2.user.parameters

*   ..find the poloidal number of outer midplane cell as the cell with min Bz

        if(nnreg(0).eq.8) then
          n_omp = rightcut(2)
          iy = min(topcut(1),topcut(2))
          do ix = rightcut(2)+1, rightcut(1)-1
            if(bb(n_omp,iy,2).eq.0.0_R8) n_omp=ix
            if(abs(bb(ix,iy,2)).lt.abs(bb(n_omp,iy,2))) n_omp=ix
          enddo
        elseif(nnreg(0).ge.4) then
          n_omp = leftcut(1)
          iy = topcut(1)
          do ix = leftcut(1)+1, rightcut(1)-1
            if(.not.isRealCell(cflags(ix,iy,CELLFLAG_TYPE))) cycle
            if(bb(n_omp,iy,2).eq.0.0_R8) n_omp=ix
            if(abs(bb(ix,iy,2)).lt.abs(bb(n_omp,iy,2))) n_omp=ix
          enddo
        else
          write(*,*) ' conv_midplane: nnreg = ', nnreg(0)
          n_omp = (nx+2)/2
          iy = ny/2
          do ix = -1, nx
            if(.not.isRealCell(cflags(ix,iy,CELLFLAG_TYPE))) cycle
            if(bb(n_omp,iy,2).eq.0.0_R8) n_omp=ix
            if(abs(bb(ix,iy,2)).lt.abs(bb(n_omp,iy,2))) n_omp=ix
          enddo
        endif
        write(*,*) ' Found OMP location at ix = ' , n_omp
        iymax = iy
        do while (isRealCell(cflags(n_omp,iymax+1,CELLFLAG_TYPE)))
          iymax = iymax + 1
        end do
        iymin = iy
        do while (isRealCell(cflags(n_omp,iymin-1,CELLFLAG_TYPE)))
          iymin = iymin - 1
        end do
        write(*,*) 'Inner and outermost flux tube at OMP at '//
     &             'iymin, iymax  = ',iymin, iymax

        

        ro(1)=0.25_R8*(crx(n_omp,iymin,0)+crx(n_omp,iymin,1)+ ! R-coordinate of center of [n_omp,-1] cell
     &                 crx(n_omp,iymin,2)+crx(n_omp,iymin,3))
        zo(1)=0.25_R8*(cry(n_omp,iymin,0)+cry(n_omp,iymin,1)+ ! Z-coordinate of center of [n_omp,-1] cell
     &                 cry(n_omp,iymin,2)+cry(n_omp,iymin,3))
        ro(2)=0.25_R8*(crx(n_omp,iymax,0)+crx(n_omp,iymax,1)+ ! R-coordinate of center of [n_omp,ny] cell
     &                 crx(n_omp,iymax,2)+crx(n_omp,iymax,3))
        zo(2)=0.25_R8*(cry(n_omp,iymax,0)+cry(n_omp,iymax,1)+ ! Z-coordinate of center of [n_omp,ny] cell
     &                 cry(n_omp,iymax,2)+cry(n_omp,iymax,3))
        if(ro(2).gt.ro(1)) then
          rzomp(1,1)=ro(1)*0.9_R8                                         ! take a coordinate to the left of ro(1)
          rzomp(2,1)=zo(1)+                                               ! calculate, corresponding to rzomp(1,1), Z-coordinate of point
     &              (rzomp(1,1)-ro(1))*(zo(2)-zo(1))/(ro(2)-ro(1))        ! laying on the straight line going through the points ro(1),zo(1) and ro(2),zo(2)
          rzomp(1,2)=ro(2)*1.1_R8                                         ! take a coordinate to the right of ro(2)
          rzomp(2,2)=zo(1)+                                               ! calculate, corresponding to rzomp(1,2), Z-coordinate of point
     &              (rzomp(1,2)-ro(1))*(zo(2)-zo(1))/(ro(2)-ro(1))        ! laying on the straight line going through the points ro(1),zo(1) and ro(2),zo(2)
        else
          write(*,*) ' ro(1:2) :', ro(1), ro(2)
          write(*,*) ' zo(1:2) :', zo(1), zo(2)
          call xerrab('conv_midplane: unexpected coordinates for rzomp')
        endif

*   ..find the poloidal number of the inner midplane cell as the cell with max Bz

        if(nnreg(0).eq.8) then
          n_imp = leftcut(2)
          iy = min(topcut(1),topcut(2))
          do ix = leftcut(1)+1, leftcut(2)-1
            if(.not.isRealCell(cflags(ix,iy,CELLFLAG_TYPE))) cycle
            if(bb(n_imp,iy,2).eq.0.0_R8) n_imp=ix
            if(abs(bb(ix,iy,2)).gt.abs(bb(n_imp,iy,2))) n_imp=ix
          enddo
        elseif(nnreg(0).ge.4) then
          n_imp = leftcut(1)
          iy = topcut(1)
          do ix = leftcut(1)+1, rightcut(1)-1
            if(.not.isRealCell(cflags(ix,iy,CELLFLAG_TYPE))) cycle
            if(bb(n_imp,iy,2).eq.0.0_R8) n_imp=ix
            if(abs(bb(ix,iy,2)).gt.abs(bb(n_imp,iy,2))) n_imp=ix
          enddo
        else
          write(*,*) ' conv_midplane: nnreg = ', nnreg(0)
          n_imp = (nx+2)/4
          iy = ny/2
          do ix = -1, nx
            if(.not.isRealCell(cflags(ix,iy,CELLFLAG_TYPE))) cycle
            if(bb(n_imp,iy,2).eq.0.0_R8) n_imp=ix
            if(abs(bb(ix,iy,2)).gt.abs(bb(n_imp,iy,2))) n_imp=ix
          enddo
        endif
        write(*,*) ' Found IMP location at ix = ' , n_imp
        iymax = iy
        do while (isRealCell(cflags(n_imp,iymax+1,CELLFLAG_TYPE)))
          iymax = iymax + 1
        end do
        iymin = iy
        do while (isRealCell(cflags(n_imp,iymin-1,CELLFLAG_TYPE)))
          iymin = iymin - 1
        end do
        write(*,*) 'Inner and outermost flux tube at OMP at '//
     &             'iymin, iymax  = ', iymin, iymax
       
        ri(1)=0.25_R8*(crx(n_imp,iymax,0)+crx(n_imp,iymax,1)+ ! R-coordinate of center of [n_imp,ny] cell
     &                 crx(n_imp,iymax,2)+crx(n_imp,iymax,3))
        zi(1)=0.25_R8*(cry(n_imp,iymax,0)+cry(n_imp,iymax,1)+ ! Z-coordinate of center of [n_imp,ny] cell
     &                 cry(n_imp,iymax,2)+cry(n_imp,iymax,3))
        ri(2)=0.25_R8*(crx(n_imp,iymin,0)+crx(n_imp,iymin,1)+ ! R-coordinate of center of [n_imp,-1] cell
     &                 crx(n_imp,iymin,2)+crx(n_imp,iymin,3))
        zi(2)=0.25_R8*(cry(n_imp,iymin,0)+cry(n_imp,iymin,1)+ ! Z-coordinate of center of [n_imp,-1] cell
     &                 cry(n_imp,iymin,2)+cry(n_imp,iymin,3))

        if(ri(2).gt.ri(1)) then
          rzimp(1,1)=ri(1)*0.9_R8                                         ! take a coordinate to the left of ri(1)
          rzimp(2,1)=zi(1)+                                               ! calculate, corresponding to rzimp(1,1), Z-coordinate of point
     &              (rzimp(1,1)-ri(1))*(zi(2)-zi(1))/(ri(2)-ri(1))        ! laying on the straight line going through the points ri(1),zi(1) and ri(2),zi(2)
          rzimp(1,2)=ri(2)*1.1_R8                                         ! take a coordinate to the right of ri(2)
          rzimp(2,2)=zi(1)+                                               ! calculate, corresponding to rzimp(1,2), Z-coordinate of point
     &              (rzimp(1,2)-ri(1))*(zi(2)-zi(1))/(ri(2)-ri(1))        ! laying on the straight line going through the points ri(1),zi(1) and ri(2),zi(2)
        else
          write(*,*) ' ri(1:2) :', ri(1), ri(2)
          write(*,*) ' zi(1:2) :', zi(1), zi(2)
          call xerrab('conv_midplane: unexpected coordinates for rzimp')
        endif

#ifdef F2003
        open(newunit=iout,file='info.dat')
#else
        iout = newunit()
        open(unit=iout,file='info.dat')
#endif
        write(iout,'(a,i4)') ' Number of outer midplane cell =', n_omp
        write(iout,'(a,f6.2,a,f6.2)') ' Ro1 = ', ro(1),' Zo1 = ', zo(1)
        write(iout,'(a,f6.2,a,f6.2)') ' Ro2 = ', ro(2),' Zo2 = ', zo(2)
        write(iout,'(a,i4)') ' Number of inner midplane cell =', n_imp
        write(iout,'(a,f6.2,a,f6.2)') ' Ri1 = ', ri(1),' Zi1 = ', zi(1)
        write(iout,'(a,f6.2,a,f6.2)') ' Ri2 = ', ri(2),' Zi2 = ', zi(2)
        close(iout)
#ifdef F2003
        open(newunit=iout,file='b2us.user.parameters')
#else
        open(unit=iout,file='b2us.user.parameters')
#endif
        write(iout,'(a)') '&user'
        write(iout,'(a,f8.2,a)') 'rzomp(1,1) = ',rzomp(1,1),','
        write(iout,'(a,f8.2,a)') 'rzomp(2,1) = ',rzomp(2,1),','
        write(iout,'(a,f8.2,a)') 'rzomp(1,2) = ',rzomp(1,2),','
        write(iout,'(a,f8.2,a)') 'rzomp(2,2) = ',rzomp(2,2),','
        write(iout,'(a,f8.2,a)') 'rzimp(1,1) = ',rzimp(1,1),','
        write(iout,'(a,f8.2,a)') 'rzimp(2,1) = ',rzimp(2,1),','
        write(iout,'(a,f8.2,a)') 'rzimp(1,2) = ',rzimp(1,2),','
        write(iout,'(a,f8.2,a)') 'rzimp(2,2) = ',rzimp(2,2),','
        write(iout,"(a)") '/'
        close (iout)

      return
      end subroutine conv_midplane                                       !srv 16.11.22 }

******************************************************************************

      end module b2us_prep

