      subroutine make_plasma_triangles (mpg, gm, nx, ny,
     .  nmxnodes, nmxelem, nelem,
     .  vert_node, elem, elem_cv, trifc, neighbor)

      use b2mod_types
      !use b2mod_boundary_namelist
      !use b2mod_neutrals_namelist
      use b2us_map
      use b2us_geo

      implicit none

      integer, intent(in) :: nx, ny, nmxnodes, nmxelem
      integer :: nelem
      type(mapping), intent(in) :: mpg
      type(geometry), intent(in) :: gm
      integer :: vert_node(mpg%nVx), elem(3,nmxelem), elem_cv(nmxelem),
     .  trifc(3,nmxelem), neighbor(3,2,nmxelem)

      ! local types to build the triangles
      type :: elem_queue
        integer :: iv1, iv2, iv3, ic, ifc1, ifc2, ifc3
        type(elem_queue), pointer :: next
      end type 

      type :: tri
        integer :: notri
        type(tri), pointer :: nxt_tri
      end type tri

      type :: tri_list
        integer :: co
        type(tri), pointer :: p
      end type tri_list

      type(elem_queue), pointer :: qhead, qlast, qp
      type(tri_list), allocatable :: ctri(:)
      type(tri), pointer :: ct

      ! local variables
      integer :: ic, i1, i2, i3, mxCvFc, ielem, is, ifc, i, jfc, 
     .           j, is1, js, js1, imin, imax, jmin, jmax, jelem, 
     .           ifcmin, ifcmax, it, 
     .           unused, ifc1, ifc2, ifc3,
     .           ileft, last(1), isig(3)
      integer, allocatable :: ifcused(:)
      logical, allocatable :: fcused(:)
      logical :: ffound
      intrinsic pack

      nullify(qhead)
      nullify(qlast)

      mxCvFc = maxval(mpg%cvFcP(:,2))
      allocate (fcused(mxCvFc))
      allocate (ifcused(mxCvFc))

      nelem = 0

!  cut each plasma cell into triangles and store the vertices in elem_queue
      do ic = 1, mpg%nCi
        
        if (mpg%cvVxP(ic,2) == 3) then

!  cell is a triangle, store vertex numbers
          i1 = mpg%cvVx(mpg%cvVxP(ic,1))
          i2 = mpg%cvVx(mpg%cvVxP(ic,1)+1)
          i3 = mpg%cvVx(mpg%cvVxP(ic,1)+2)
          ifc1 = find_face(i1,i2,ic)
          ifc2 = find_face(i2,i3,ic)
          ifc3 = find_face(i3,i1,ic)
          call set_element(i1,i2,i3,ic,ifc1,ifc2,ifc3)
          nelem = nelem + 1

        else

!  cell has more then 3 vertices and needs to be split in triangles
          fcused = .false.
          do ifc = 1, mpg%cvFcP(ic,2)
! face numbers for cell is
            ifcused(ifc) = mpg%cvFc(mpg%cvFcP(ic,1)+ifc-1) 
          end do
!  find faces connected with first vertex of cell
          i1 = mpg%cvVx(mpg%cvVxP(ic,1))  ! first vertex of cell ic
          i2 = 0
          i3 = 0
!  we need to identify the two faces connected to vertex i1 of cell ic
!  providing us with vertices i2 and i3 of the first triangle
!  taken from cell ic 
         do i = 1, mpg%vxFcP(i1,2)   ! for all faces connected to vertex i1
!           ifc = mpg%vxFc(mpg%vxFcP(i1,1)+i-1)
!           if ((mpg%fcCv(ifc,1) == ic) .or.   ! cell ic is connected
!     .         (mpg%fcCv(ifc,2) == ic)) then  ! to face ifc
!  look for faces of cell ic which are connected to vertex i1 
              jfc = 0
              do j = 1, mpg%cvFcP(ic,2) ! for all faces belonging to cell ic
                ifc = ifcused(j)
!               if (.not.fcused(j) .and. (ifcused(j) == ifc)) then
                if (.not.fcused(j)) then
!  check if first vertex i1 of cell ic belongs to face ifc
                  if ((mpg%fcVx(ifc,1) == i1) .or.
     .                (mpg%fcVx(ifc,2) == i1)) then
                    jfc = j   ! face starting at vertex i1 found
                    exit
                  end if
                end if
              end do
              if ((jfc > 0) .and. .not.fcused(jfc)) then
                fcused(jfc) = .true.  ! mark face as used
!  identified the second vertex of face ifc
!  one vertex of face ifc is equal to i1, thus summing up the 
!  vertex numbers and subtracting i1 will provide the vertex number we require
                if (i2 == 0) then
                  i2 = mpg%fcVx(ifc,1) + mpg%fcVx(ifc,2) - i1 
                  ifc1 = ifc
                else if (i3 == 0) then
                  i3 = mpg%fcVx(ifc,1) + mpg%fcVx(ifc,2) - i1
                  ifc3 = ifc
                  exit
                else
                  call xerrab('write_triangles -- too many faces')
                end if
              else
                call xerrab('write_triangles -- wrong face')
              end if
!           end if
         end do

          ifc2 = 0  ! new face inside cell ic 
          call set_element(i1,i2,i3,ic,ifc1,ifc2,ifc3)
          nelem = nelem + 1

!  first triangle finished, now find all others

          do 
            ileft = count(.not.fcused)
            if (ileft == 1) exit
!  use first vertex of side 2 of previously found triangle as starting point
            i1 = qlast%iv2
            i2 = qlast%iv3
            ifc1 = qlast%ifc2
            i3 = 0
            do jfc = 1, mpg%cvFcP(ic,2)
              if (fcused(jfc)) cycle
              ifc = ifcused(jfc)
              if ((mpg%fcVx(ifc,1) == i1) .or. 
     .            (mpg%fcVx(ifc,2) == i1)) then
!  face needed to form the next triangle found
                i3 = mpg%fcVx(ifc,1) + mpg%fcVx(ifc,2) - i1 
                ifc3 = ifc
                fcused(jfc) = .true.
!  only two usable faces left => both belong to triangle
                if (ileft == 2) then
!  identify face number of last triangle side
                  last = pack(ifcused,.not.fcused) 
                  ifc2 = last(1)
                end if
                exit
              end if
            end do
            if (i3 == 0) 
     .        call xerrab ('write_triangles -- ')

            call set_element(i1,i2,i3,ic,ifc1,ifc2,ifc3)
            nelem = nelem + 1
            
          end do

!  check that the last unused face of cell ic is part of the
!  last constructed triangle
          unused = 0
          do jfc = 1, mpg%cvFcP(ic,2)
            if (.not.fcused(jfc)) unused = jfc
          end do
          ifcmin = min(mpg%fcVx(ifcused(unused),1),
     .                 mpg%fcVx(ifcused(unused),2))
          ifcmax = max(mpg%fcVx(ifcused(unused),1),
     .                 mpg%fcVx(ifcused(unused),2))
          if (.not.(((ifcmin == min(i1,i2)) .or. 
     .               (ifcmin == min(i2,i3)) .or.
     .               (ifcmin == min(i1,i3))) .and.
     .              ((ifcmax == max(i1,i2)) .or. 
     .               (ifcmax == max(i2,i3)) .or.
     .               (ifcmax == max(i1,i3)))))
     .      call xerrab ('write_triangles -- wrong face')                              

        end if
        
      end do

!  fetch element information from queue and store it in arrays for easier access

      allocate (ctri(mpg%nVx))
      do i = 1, mpg%nVx
        ctri(i)%co = 0
        nullify(ctri(i)%p)
      end do

      ielem = 0
      do while (associated(qhead))
        ielem = ielem + 1
        if (ielem > nelem) 
     .    call xerrab('write_triangles -- element number')
        qp => qhead
        qhead => qhead%next
        if (vert_node(qp%iv1).eq.0 .or. vert_node(qp%iv2).eq.0 .or.
     .      vert_node(qp%iv3).eq.0 ) then
          call xerrab('make_plasma_triangles -- problem node def.')
        endif
        elem(1:3,ielem) = (/ vert_node(qp%iv1), vert_node(qp%iv2), 
     .    vert_node(qp%iv3) /)
        elem_cv(ielem) = qp%ic
        trifc(1:3,ielem) = (/ qp%ifc1, qp%ifc2, qp%ifc3 /) 
        isig = 1
        if (trifc(1,ielem) > 0) then
          if (elem_cv(ielem) == mpg%fcCv(trifc(1,ielem),2)) isig(1) = -1
        end if
        if (trifc(2,ielem) > 0) then
          if (elem_cv(ielem) == mpg%fcCv(trifc(2,ielem),2)) isig(2) = -1
        end if
        if (trifc(3,ielem) > 0) then
          if (elem_cv(ielem) == mpg%fcCv(trifc(3,ielem),2)) isig(3) = -1
        end if
        trifc(1:3,ielem) = trifc(1:3,ielem) * isig
        
        call add_tri(vert_node(qp%iv1), ielem)
        call add_tri(vert_node(qp%iv2), ielem)
        call add_tri(vert_node(qp%iv3), ielem)
      end do
      if (ielem /= nelem) 
     .    call xerrab('write_triangles -- total element number')

!  find neighbors
      neighbor = 0
      
      do ielem = 1, nelem
        do is = 1, 3
          if (neighbor(is,1,ielem) == 0) then
            is1 = is + 1
            if (is1 > 3) is1 = 1
            imin = min(elem(is,ielem),elem(is1,ielem))
            imax = max(elem(is,ielem),elem(is1,ielem))
            ct => ctri(elem(is,ielem))%p
!            write (*,*) ' ielem, is ',ielem, is,
!     .        ' no of triangles ',ctri(elem(is,ielem))%co
            jloop: do while (associated(ct))
              jelem = ct%notri
              if (jelem /= ielem) then  ! do not check triangle ielem
                do js = 1, 3
                  if (neighbor(js,1,jelem) == 0) then
                    js1 = js + 1
                    if (js1 > 3) js1 = 1
                    jmin = min(elem(js,jelem),elem(js1,jelem))
                    jmax = max(elem(js,jelem),elem(js1,jelem))
                    if ((imin == jmin) .and. (imax == jmax)) then ! neighbor found
                      neighbor(is,1,ielem) = jelem
                      neighbor(is,2,ielem) = js
                      neighbor(js,1,jelem) = ielem
                      neighbor(js,2,jelem) = is
                      exit jloop
                    end if
                  end if
                end do ! js
              end if
              ct => ct%nxt_tri
            end do jloop
          end if
        end do
      end do

      contains 

      subroutine add_tri (iv, ic)
      integer, intent(in) :: iv, ic
        allocate (ct)
        ct%notri = ic
        nullify(ct%nxt_tri)
        if (associated(ctri(iv)%p)) then
          ct%nxt_tri => ctri(iv)%p
        end if
        ctri(iv)%p => ct
        ctri(iv)%co = ctri(iv)%co + 1
      end subroutine add_tri

      function areat(x1, y1, x2, y2, x3, y3)
      real (kind=R8) :: areat
      real (kind=R8), intent(in) :: x1, y1, x2, y2, x3, y3
        areat = 0.5_R8*(x1*(y2-y3)+x2*(y3-y1)+x3*(y1-y2))
      end function areat

      subroutine set_element(i1, i2, i3, ic, ifc1, ifc2, ifc3)
      integer, intent(in) :: i1, i2, i3, ic, ifc1, ifc2, ifc3
      real(kind=R8) :: at
        allocate(qp)
        nullify(qp%next)
        qp%ic = ic
        qp%iv1 = i1
        qp%iv2 = i2
        qp%iv3 = i3
        qp%ifc1 = ifc1
        qp%ifc2 = ifc2
        qp%ifc3 = ifc3
        at = areat(gm%vxX(qp%iv1),gm%vxY(qp%iv1),
     .             gm%vxX(qp%iv2),gm%vxY(qp%iv2),
     .             gm%vxX(qp%iv3),gm%vxY(qp%iv3))
        if (at < 0._R8) call switch_points
        call push_queue
      end subroutine set_element

      subroutine switch_points
      integer :: ih, ifh
        ih = qp%iv2
        qp%iv2 = qp%iv3
        qp%iv3 = ih
        ifh = qp%ifc1
        qp%ifc1 = qp%ifc3
        qp%ifc3 = ifh
      end subroutine switch_points
      
      subroutine push_queue       
        if (.not.associated(qhead)) then
          qhead => qp
          qlast => qp
        else
          qlast%next => qp
          qlast => qp
        end if
      end subroutine push_queue

      function find_face(i1, i2, ic) result(ifc)
      integer, intent(in) :: i1, i2, ic
      integer :: ifc, j1, j2, i
        do i = 1, mpg%cvFcP(ic,2)
          ifc = mpg%cvFc(mpg%cvFcP(ic,1)+i-1)
          j1 = mpg%fcVx(ifc,1) 
          j2 = mpg%fcVx(ifc,2) 
          if (((i1 == j1) .and. (i2 == j2)) .or.
     .        ((i1 == j2) .and. (i2 == j1))) return
        end do
        ifc = 0
        return
      end function find_face
       
      end subroutine make_plasma_triangles