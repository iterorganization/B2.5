      subroutine conv_triangles (mpg, gm,
     .    nmxnodes, nmxelem, nnodes, nelem, npnodes, npelem,
     .    xco, yco, node_vert, elem, elem_cv, trifc, neigh,
     .    cntr, tria_mode)

      use b2mod_types
      use b2mod_indirect
      use b2us_map
      use b2us_geo
      use b2us_plasma
      use b2us_prep
#ifdef B25_EIRENE
      use eirmod_module_avltree
#endif
      implicit none

      integer, intent(in) :: nmxnodes, nmxelem, npnodes, npelem,
     .  tria_mode
      integer ::  nnodes, nelem
      type(mapping), intent(in) :: mpg
      type(geometry), intent(in) :: gm
      integer :: node_vert(nmxnodes), elem(3,nmxelem), elem_cv(nmxelem),
     .  trifc(3,nmxelem), neigh(3,2,nmxelem), cntr(3,nmxelem)
      real(R8) :: xco(nmxnodes), yco(nmxnodes)

#ifdef B25_EIRENE
      type(TAVLTree), pointer :: baum
      logical :: inserted
      integer :: nx, ny
      integer :: nco_old, nelem_conv, ic, ico, ivx, iel, iel2,
     .           i1, i2, i3, idum, i, nel2,
     .           is, is1, is2, ifc, icv, icut, ixt, iyt,
     .           ix, iy, isig, iused, ineir,
     .           ntrfrm, ismin, ismax, icmin, icmax, cfce(4,3),
     .           ixl, ixr, iv, icside

      real(R8), allocatable :: xc_old(:), yc_old(:)
      integer, allocatable :: im_old_new(:), im_new_old(:)
      integer, allocatable :: ixtri(:), iytri(:),
     .     ix_b2_eir(:), ix_eir_b25(:), iy_eir(:)
      integer :: artificial_slab = 0, overwrite_ntrfrm = 0,
     .           is_narrow_grid = 0, use_old_co = 0
      real(R8) :: xc, yc, dist = 1.0_R8
      real(R8) :: tol_align = 1.0e-4_R8
      logical :: ffound
      character(500) :: line

      type :: tri_elem
        integer :: itri
        type(tri_elem), pointer :: next_tri => null()
      end type tri_elem

      type :: tri_list
        type(tri_elem), pointer :: ptri => null()
      end type tri_list

      type(tri_list), allocatable :: trias(:,:)
      type(tri_elem), pointer :: tel
#ifndef F2003
      integer newunit
      external newunit
#endif

      call ipgeti ('artificial_slab', artificial_slab)
      call ipgeti ('overwrite_ntrfrm', overwrite_ntrfrm)
      call ipgeti ('is_narrow_grid', is_narrow_grid)
      call ipgetr ('conv_triangles_dist',dist)
      call ipgeti ('conv_triangles_old_co', use_old_co)
      call ipgetr ('conv_triangles_tol', tol_align)

!  if input.dat exists check which form of fort.33 is used
      inquire (file='input.dat',exist=ffound)
      if (ffound) then
#ifdef F2003
        open (newunit=ineir,file='input.dat')
#else
        ineir = newunit()
        open (unit=ineir,file='input.dat')
#endif
        line = repeat(' ',500)
        do while (index(line,'*** 14') == 0)
          read (ineir,'(a500)') line
        end do
        read (ineir,'(a500)') line
        read (ineir,'(12i6)') idum, ntrfrm
        close (ineir)
      else
        ntrfrm = 0
      endif

!  read coordinates from fort.33
      open (unit=33,access='SEQUENTIAL',form='FORMATTED')
      read (33,*) nco_old

      allocate (xc_old(nco_old))
      allocate (yc_old(nco_old))
      allocate (im_old_new(nco_old+mpg%nVx))
      allocate (im_new_old(nco_old+mpg%nVx))
      im_old_new = 0
      im_new_old = 0

      if (overwrite_ntrfrm.eq.1) then
        ntrfrm = 1
      endif
      if (ntrfrm == 0) then
        read (33,*) (xc_old(i),i=1,nco_old)
        read (33,*) (yc_old(i),i=1,nco_old)
      else
        do i = 1, nco_old
          read (33,*) idum, xc_old(i), yc_old(i)
        end do
      end if
      close (33)

!  now find new coordinate indices

!  first introduce coordinates from unstructured mesh into
!  self-balancing binary search tree
      baum => EIRENE_NewTree()
      call eirene_setcon
      do i = 1, npnodes
        IC=i
        inserted=.false.
        xc = xco(i)
        yc = yco(i)
        call EIRENE_insert (baum, xc, yc,
     .                      0.0_R8, dist, ic, inserted)
        if (.not.inserted) then
            write(*,*) '------------'
            write(*,*) 'problem with plasma node ', i
            write(*,*) 'x = ', xc
            write(*,*) 'y = ', yc
            write(*,*) 'not inserted in node list'
            write(*,*) '------------'
        endif
      end do


!  now check coordinates from old fort.33 file for new index
      nnodes = npnodes
      do ico = 1, nco_old
        IC=NNODES+1
        inserted=.false.
        call EIRENE_insert (baum, xc_old(ico), yc_old(ico),
     .                      0._R8, dist, ic, inserted)
        IF (INSERTED) THEN
          if (is_narrow_grid.eq.1) then
            !! if we have a narrow grid, this should not happen because every node of the Eirene triangle grid should be in the plasma grid.
            !! therefore this has to do with round-off error of faulty input files. Therefore these cells are printed.
            write(*,*) '------------'
            write(*,*) 'Narrow grid - faulty node ', ico
            write(*,*) 'x = ', xc_old(ico)
            write(*,*) 'y = ', yc_old(ico)
            write(*,*) '------------'
          else
            NNODES=NNODES+1
            if (NNODES .gt. (nmxnodes))
     .        call xerrab ('insufficient array space for new nodes')
            xco(NNODES) = xc_old(ico)
            yco(NNODES) = yc_old(ico)
            !write (*,*) 'Inserting node ', nnodes, ' for old node ', ico
            im_old_new(ico) = nnodes
            im_new_old(nnodes) = ico
          endif
        else
          if (im_new_old(ic) .ne.0) then
            write (*,*) 'Problem converting nodes:'
            write (*,*) 'Trying to link old node ', ico, ' to new '//
     .                  'node ', ic, 'but new node already linked '//
     .                  'to old node ', im_new_old(ic)
!            call xerrab ('node already used')
          endif
          !write (*,*) 'Linking node ', ic, ' to old node ', ico'
          if (use_old_co.eq.1) then
            xco(ic) = xc_old(ico)
            yco(ic) = yc_old(ico)
          endif
          im_old_new(ico) = ic
          im_new_old(ic) = ico
        END IF

      end do

      !! Some sanity checks
      !! 1) check whether all nodes from original fort.33 have been
      !!    included in conversion
      ico = 0
      do i = 1, nco_old
        iused = count(im_new_old == i)
        if (iused <= 0) then
          ico = ico + 1
          write (*,*) 'Node ', i, ' from original fort.33 not used'
          write (*,*) 'X, Y = ', xc_old(i), yc_old(i)
        elseif (iused > 1) then
          ico = ico + 1
          write (*,*) 'Node ', i, ' from original fort.33 linked to ' //
     .                'multiple new nodes after conversion'
          write (*,*) 'X, Y = ', xc_old(i), yc_old(i)
        endif
      enddo

      !! 2) check whether each new node corresponds to an old node
      do i = 1, nnodes
        iused = count(im_old_new == i)
        if (iused <= 0 .and.
     .      .not.(i.le.npnodes.and.tria_mode.eq.1)) then ! for tria_mode 1, only chech the nodes that are not plasma vertices
          ico = ico + 1
          write (*,*) 'Node ', i, ' in new fort.33 not present in '//
     .                'original file'
          write (*,*) 'X, Y = ', xco(i), yco(i)
        elseif (iused > 1) then
          ico = ico + 1
          write (*,*) 'Node ', i, ' in new fort.33 linked to ' //
     .                'multiple nodes in original file'
          write (*,*) 'X, Y = ', xco(i), yco(i)
        endif
      enddo

      !! 3) check whether plasma vertices are correctly linked to a node
      if (nnodes /= nco_old.and.tria_mode.eq.0) then
        write (*,*) 'Created ', nnodes, ' nodes, ',
     .              'but expected ', nco_old, ' nodes'
        !! check wether each plasma vertex has been identified as a
        !! node in fort.33
        do i = 1, NPNODES
          if (im_new_old(i) == 0) then
            ivx = node_vert(i)
            ico = ico + 1
            write (*,*)
     .       ' Vertex ivx, vxX, vxY = ', ivx, gm%vxX(ivx), gm%vxY(ivx),
     .       ' not found in original fort.33 file'
          end if
        end do
      end if
      if (ico > 0) then
        write (*,*) 'Problem(s) converting coordinates in fort.33'
        write (*,*) 'Check whether triangle grid and plasma grid match,'
        write (*,*) 'with sufficient numerical precision.'
!        call xerrab ('conv_triangles -- problem with node conversion')
      endif

!  conversion of coordinates finished

!  now read element information from old fort.34
      open (unit=34,access='SEQUENTIAL',form='FORMATTED')
      read (34,*) nelem_conv

!  assume npelem elements may already have been present (plasma)
!  => append the additional nelem_conv ones from fort.34 at the end
      nelem = npelem + nelem_conv
      do iel = npelem + 1, nelem
        read (34,*) idum, i1, i2, i3
        elem(1,iel) = im_old_new(i1)
        elem(2,iel) = im_old_new(i2)
        elem(3,iel) = im_old_new(i3)
      end do
      close(34)

!  conversion of elements finished

!  now read neighbor information from old fort.35
      open (unit=35,access='SEQUENTIAL',form='FORMATTED')
      read (35,*) nel2

      if (nel2 /= nelem_conv)
     .  call xerrab ('conv_triangles -- '//
     .               'wrong no. of triangles in fort.35')

      if (tria_mode.eq.1) then ! assume new format fort.35; void triangles only, already unstructured format

        do i = npelem + 1, nelem
          read (35,*) idum,(neigh(is,1,i),neigh(is,2,i),cntr(is,i),
     .                is = 1,3),
     .                elem_cv(i), trifc(1,i), trifc(2,i), trifc(3,i)
        end do
        close(35)

        ! consistency check
        ico = 0
        do iel = npelem + 1, nelem
          if (elem_cv(iel).ne.-1 .or. trifc(1,iel).ne.0 .or.
     .        trifc(2,iel).ne.0  .or. trifc(3,iel).ne.0 ) then
            ico = ico + 1
            write (*,*) 'Problem with void triangle ',iel
            write (*,*) 'iCv, iFc1, iFc2, iFc3 = ', elem_cv(iel),
     .       trifc(1:3,iel)
          endif
        enddo
        call xertst(ico.eq.0, 'Problem with void triangles')

      elseif (tria_mode.eq.0) then ! assume structured format fort.35; conversion of ix,iy indexing needed

        ! set size of structured grid to convert plasma triangles
        if (mpg%isClassicalGrid == 1) then
          nx = mpg%nx
          ny = mpg%ny
        else
          call xerrab 
     .     ('Need structured grid info to convert plasma triangles.')
        endif

        allocate (ixtri(nelem_conv))
        allocate (iytri(nelem_conv))

        do i = npelem + 1, nelem
          read (35,*) idum,(neigh(is,1,i),neigh(is,2,i),cntr(is,i),
     .                is = 1,3), ixtri(i), iytri(i)
        end do
        close(35)

!  setup conversion arrays

        allocate (ix_b2_eir(1:nx))
        allocate (ix_eir_b25(1:nx+3*nncut))
        ix_eir_b25 = -2

        ix_b2_eir = [ (ix, ix = 1, nx) ]
        do icut = 1, nncut
          ix_b2_eir(leftcut(icut)+1:nx)=ix_b2_eir(leftcut(icut)+1:nx)+1
          ix_b2_eir(rightcut(icut)+1:nx)=
     .       ix_b2_eir(rightcut(icut)+1:nx)+1
        end do
        do icut = 1, mpg%nctcnt
          ix_b2_eir(mpg%ictcntrl(icut)+1:nx)=
     .       ix_b2_eir(mpg%ictcntrl(icut)+1:nx)+1
        end do

        do ix = 1, nx
          ixt = ix_b2_eir(ix)
          ix_eir_b25(ixt) = ix - 1
        end do

        allocate (iy_eir(ny))
        iy_eir = 0
        do iy = 1, ny
          iy_eir(iy) = iy - 1
        end do

!  convert plasma cell index from (ixtri,iytri) to icv
        if (artificial_slab > 0) then
          allocate(trias(-1:nx,-1:ny))
        end if

        do iel = npelem + 1, nelem
          if (ixtri(iel) < 0) then
            elem_cv(iel) = -1
          else
!  ixtri, iytri are in Eirene cell counting
!  find cell number in B2.5 counting
            ixt = ix_eir_b25(ixtri(iel))
            iyt = iy_eir(iytri(iel))
            if (ixt == -2) then
              elem_cv(iel) = -1
            else
              icv = mpg%imapCv(ixt,iyt)
              elem_cv(iel) = icv
            end if
            if  (artificial_slab .eq. 1) then
              if ((ixt > -2) .and. (ixt <= nx) .and.
     .            (iyt > -2) .and. (iyt <= ny)) then
                allocate (tel)
                tel%itri = iel
                tel%next_tri => trias(ixt,iyt)%ptri
                trias(ixt,iyt)%ptri => tel
              end if
            end if
          end if
        end do
      endif

      ! shift neighbor info for presence of (plasma) triangles
      do iel = npelem + 1, nelem
        do is = 1,3
          if (neigh(is,1,iel).ne.0)
     .      neigh(is,1,iel) = neigh(is,1,iel) + npelem
        enddo
      enddo

!  find face numbers for triangle sides
      do iel = npelem + 1, nelem
        if (elem_cv(iel) < 0) then
!  triangle belongs to void region
!  check if it is neighboring a plasma cell
          do is = 1, 3
            if (triFc(is,iel) == 0) then
              is1 = is + 1
              if (is1 > 3) is1 = 1
              if (neigh(is,1,iel) == 0) then
                ! Triangle w/o neighbor => edge of triangle grid.
                ! For stitching void grid to plasma grid: check if there is a plasma triangle on other side
                call find_plasma_neighbor (iCv, iFc, iel2, is2,
     .            elem(is, iel), elem(is1, iel), nnodes, node_vert,
     .            npelem, nelem, elem, mpg)
                ! If yes, adapt neighbor information
                if (iCv > 0) then
                  neigh(is,1,iel) = iel2
                  neigh(is,2,iel) = is2
                  neigh(is2,1,iel2) = iel
                  neigh(is2,2,iel2) = is
                  isig = 1
                  if (iCv == mpg%fcCv(iFc,2)) isig = -1
                  triFc(is,iel) = -isig * iFc
                  triFc(is2,iel2) = isig * iFc
                endif
              else
                ! Triangle with neighbor => check cell of neighbor
                iCv = elem_cv(neigh(is,1,iel))
                if (iCv > 0) then
                  call find_face_number (iCv, elem(is,iel),
     .              elem(is1,iel), nnodes, npnodes, node_vert,
     .              xco, yco, ifc, mpg, gm, tol_align)
                  isig = 1
                  if (iFc > 0) then
                    if (iCv == mpg%fcCv(iFc,2)) isig = -1
                  end if
                  trifc(is,iel) = -isig * iFc
                  trifc(neigh(is,2,iel),neigh(is,1,iel)) = isig * iFc
                end if
              endif
            end if
          end do
        else
!  triangle is in plasma region
          iCv = elem_cv(iel)
          do is = 1, 3
            if (triFc(is,iel) == 0) then
              is1 = is + 1
              if (is1 > 3) is1 = 1
              call find_face_number (iCv, elem(is,iel),
     .          elem(is1,iel), nnodes, npnodes, node_vert,
     .          xco, yco, ifc, mpg, gm, tol_align)
              isig = 1
              if (iFc > 0) then
                if (iCv == mpg%fcCv(iFc,2)) isig = -1
              end if
              triFc(is,iel) = isig * iFc
              if (neigh(is,1,iel) /= 0)
     .          triFc(neigh(is,2,iel),neigh(is,1,iel)) = -isig * iFc
            end if
          end do
          if (count(trifc(1:3,iel) /= 0) == 0) then
            write (*,*) ' triangle ',iel,' matched no faces '
          end if
        end if
      end do

      if (artificial_slab > 0) then
        write (*,*) ' artificial slab is used'
        do icut = 1, nncut
          do iy = bottomcut(icut)+1, topcut(icut)-1
           do icside = 1, 2
             if (icside == 1) then
               ixr = leftcut(icut)
               ixl = leftcut(icut)-1
               icv = mpg%imapcv(ixl,iy)
! bottom face
               cfce(1,1:3)=[ mpg%imapVx(ixl,iy), mpg%imapVx(ixr,iy),
     .                       mpg%imapFCy(ixl,iy) ]
! right face
               cfce(2,1:3)=[ mpg%imapVx(ixr,iy), mpg%imapVx(ixr,iy+1),
     .                       mpg%imapFcx(rightcut(icut),iy) ]
! top face
               cfce(3,1:3)=[ mpg%imapVx(ixl,iy+1), mpg%imapVx(ixr,iy+1),
     .                       mpg%imapFCy(ixr,iy+1) ]
! left face
               cfce(4,1:3)=[ mpg%imapVx(ixl,iy), mpg%imapVx(ixl,iy+1),
     .                       mpg%imapFcx(ixl,iy) ]
               tel => trias(ixl,iy)%ptri
             else
               ixr = rightcut(icut)
               ixl = rightcut(icut)-1
               icv = mpg%imapcv(ixl,iy)
! bottom face
               cfce(1,1:3)=[ mpg%imapVx(ixl,iy), mpg%imapVx(ixr,iy),
     .                         mpg%imapFCy(ixl,iy) ]
! right face
               cfce(2,1:3)=[ mpg%imapVx(ixr,iy), mpg%imapVx(ixr,iy+1),
     .                         mpg%imapFcx(leftcut(icut),iy) ]
! top face
               cfce(3,1:3)=[ mpg%imapVx(ixl,iy+1), mpg%imapVx(ixr,iy+1),
     .                         mpg%imapFCy(ixr,iy+1) ]
! left face
               cfce(4,1:3)=[ mpg%imapVx(ixl,iy), mpg%imapVx(ixl,iy+1),
     .                         mpg%imapFcx(ixl,iy) ]
               tel => trias(ixl,iy)%ptri
             end if
            do while (associated(tel))
              iel = tel%itri
              write (*,*) ' triangle ',iel,' treated '
              do is = 1, 3
                is1 = is + 1
                if (is1 > 3) is1 = 1
                if (trifc(is,iel) == 0) then
                  ismin = min(node_vert(elem(is,iel)),
     .                        node_vert(elem(is1,iel)))
                  ismax = max(node_vert(elem(is,iel)),
     .                        node_vert(elem(is1,iel)))
                  do iv = 1, 4
                    icmin = min(cfce(iv,1),cfce(iv,2))
                    icmax = max(cfce(iv,1),cfce(iv,2))
                    if ((ismin == icmin) .and. (ismax == icmax)) then
                      ifc = cfce(iv,3)
                      if (ifc > 0) then
                        if (icv == mpg%fcCv(ifc,2)) isig = -1
                      end if
                      trifc(is,iel) = isig*ifc
                      exit
                    end if
                  end do
                end if
              end do
              if (count(trifc(1:3,iel) /= 0) == 0) then
                write (*,*) ' triangle at cut ',iel,' matched no faces '
              end if
              tel => tel%next_tri
            end do
           end do
          end do
        end do
      end if

      call EIRENE_DestroyTree(baum)
      if (allocated(xc_old)) deallocate (xc_old)
      if (allocated(yc_old)) deallocate (yc_old)
      if (allocated(im_old_new)) deallocate (im_old_new)
      if (allocated(im_new_old)) deallocate (im_new_old)
      if (allocated(ixtri)) deallocate (ixtri)
      if (allocated(iytri)) deallocate (iytri)
      if (allocated(ix_b2_eir)) deallocate (ix_b2_eir)
      if (allocated(ix_eir_b25)) deallocate (ix_eir_b25)
      if (allocated(iy_eir)) deallocate (iy_eir)
#endif

      return
      end subroutine conv_triangles

      subroutine find_face_number (iCv, in1, in2,
     .                             nnodes, npnodes, node_vert,
     .                             xco, yco, ifc, mpg, gm,
     .                             tol_align)

      use b2mod_types
      use b2us_map
      use b2us_geo

      implicit none
      integer, intent(in) :: iCv, in1, in2, nnodes, npnodes,
     .                       node_vert(nnodes)
      real(R8), intent(in) :: xco(nnodes), yco(nnodes)
      type(mapping), intent(in) :: mpg
      type(geometry), intent(in) :: gm
      integer, intent(out) :: ifc

      integer, allocatable :: jvx(:)
      integer :: nv, nf, j, jFc, jv1, jv2
      real(R8) :: px1, py1, px2, py2, vpx, vpy,
     .            fx1, fy1, fx2, fy2, vfx, vfy, d, det, dq, lq, vvp, vvf
      real(R8) :: tol_align
      logical :: found

      nv = mpg%cvVxP(iCv,2)
      allocate (jvx(nv))
      jvx = mpg%cvVx(mpg%cvVxp(iCv,1) : mpg%cvVxp(iCv,1)+nv-1)

      nf = mpg%cvFcP(iCv,2)
      iFc = 0
      found=.false.

      if ((in1 <= npnodes) .and. (in2 <= npnodes)) then
!  both in1 and in2 are nodes in the unstructured plasma mesh
!  check faces of cell icv
        do j = 1, nf
          jFc = mpg%cvFc(mpg%cvFcP(iCv,1)+j-1)
          jv1 = mpg%fcVx(jFc,1)
          jv2 = mpg%fcVx(jFc,2)
          if ((min(node_vert(in1),node_vert(in2)) == min(jv1,jv2)) .and.
     .        (max(node_vert(in1),node_vert(in2)) == max(jv1,jv2))) then
!  face identified
            iFc = jFc
            found=.true.
            return
          end if
        end do
        if (.not.found) then
!          write (*,*) 'No face found for nodes ', in1, in2
        endif
      else
!  at least one vertex was introduced in triangle creation
        px1 = xco(in1)
        py1 = yco(in1)
        px2 = xco(in2)
        py2 = yco(in2)
        vpx = px2 - px1
        vpy = py2 - py1
        vvp = sqrt(vpx**2 + vpy**2)
        vpx = vpx/vvp
        vpy = vpy/vvp
        do j = 1, nf
          jfc = mpg%cvFc(mpg%cvFcP(icv,1)+j-1)
          jv1 = mpg%fcVx(jfc,1)
          jv2 = mpg%fcVx(jfc,2)
          fx1 = gm%vxX(jv1)*100._R8
          fy1 = gm%vxY(jv1)*100._R8
          fx2 = gm%vxX(jv2)*100._R8
          fy2 = gm%vxY(jv2)*100._R8
          vfx = fx2 - fx1
          vfy = fy2 - fy1
          vvf = sqrt(vfx**2 + vfy**2)
          vfx = vfx/vvf
          vfy = vfy/vvf
!  check if p and f are parallel
          det = abs(vfx*vpy - vfy*vpx)
          if (det > tol_align) cycle
!  check that the "new" point lies on interval (x1,y1) -- (x2,y2)
          if ((node_vert(in1) /= jv1) .and. (node_vert(in1) /= jv2))then
            d = vvf*vfy*px1 - vvf*vfx*py1 + (fx2*fy1-fx1*fy2)
          else
            d = vvf*vfy*px2 - vvf*vfx*py2 + (fx2*fy1-fx1*fy2)
          end if
          dq = d*d
          lq = vvf*vvf
          if (dq < 1.D-6*LQ) then
            ifc = jfc
            return
          end if
        end do
!        write (*,*) ' No face found for nodes ', in1, in2
!        call xerrab ('find_face_number -- no face found ')
      end if

      deallocate (jvx)
      return
      end subroutine find_face_number

      subroutine find_plasma_neighbor (iCv, iFc, iel2, is2, in1, in2,
     .  nnodes, node_vert, npelem, nelem, elem, mpg)

      use b2mod_types
      use b2us_map
      use b2us_geo
      implicit none

      type(mapping), intent(in) :: mpg

      integer :: icv, ifc, icv1, ifc1, iel, iel2, is, is1, is2,
     .  in1, in2, iv1, iv2, nnodes, npelem, nelem
      integer :: node_vert(nnodes), elem(3,nelem)
      logical :: found

      iCv = -1 ! default: assume void neighbor
      iFc = -1
      iel2 = -1
      is2  = -1

      ! In the new unstructured code, we assume that plasma boundary faces
      ! will not have been split for triangulation
      ! => we simply need to detect whether nodes of the triangle face correspond to
      !    those of a plasma boundary face
      do iCv1 = mpg%nCi+1,mpg%nCv  ! loop over guard cells
        iFc1 = mpg%cvFc(mpg%cvFcP(iCv1,1)) ! number of the guard cell face
        iv1 = mpg%fcVx(iFc1,1)
        iv2 = mpg%fcVx(iFc1,2)
        if ((min(node_vert(in1),node_vert(in2)) == min(iv1,iv2)) .and.
     .      (max(node_vert(in1),node_vert(in2)) == max(iv1,iv2))) then
          ! Cell and face number identified
          iCv = min(mpg%fcCv(iFc1,1),mpg%fcCv(iFc1,2))
          iFc = iFc1
          ! Now, find the plasma triangle that has this face in common
          found = .false.
          do iel = 1, npelem
            do is = 1, 3
              is1 = is + 1
              if (is1 > 3) is1 = 1
              iv1 = elem(is,iel)
              iv2 = elem(is1,iel)
              if ((min(in1,in2) == min(iv1,iv2)) .and.
     .            (max(in1,in2) == max(iv1,iv2))) then
                if (found) call xerrab ('neibhbor already found')
                found = .true.
                iel2 = iel
                is2 = is
              endif
            enddo
          enddo
          if (.not.found) then
            write (*,*) 'Did not find triangle for iCv, iFc =  ',iCv,iFc
            call xerrab ('Problem stitching plasma & void grids.')
          endif
        endif
      enddo

      return
      end subroutine find_plasma_neighbor
