      subroutine write_triangles (nout, mpg, gm)

      use b2mod_types
      use b2mod_boundary_namelist
      use b2mod_neutrals_namelist
      use b2us_map
      use b2us_geo

      implicit none

      integer, intent(in) :: nout
      type(mapping), intent(in) :: mpg
      type(geometry), intent(in) :: gm

      type :: elem_queue
        integer :: iv1, iv2, iv3, ic, iFc1, iFc2, iFc3
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

      integer :: ic, i1, i2, i3, mxCvFc, nelem, ielem, is, iFc, i, jfc,
     .           j, is1, js, js1, imin, imax, jmin, jmax, jelem,
     .           iFcmin, iFcmax, it,
     .           unused, nopen, iFc1, iFc2, iFc3,
     .           ileft, last(1), isig(3)
      integer, allocatable :: iFcused(:), elem(:,:), elem_cv(:),
     .                        neighbor(:,:,:), bndmdl(:,:),
     .                        triFc(:,:), node_vert(:)
      logical, allocatable :: fcused(:)
      logical :: ffound
      intrinsic pack

      nullify(qhead)
      nullify(qlast)

      mxCvFc = maxval(mpg%cvFcP(:,2))
      allocate (fcused(mxCvFc))
      allocate (iFcused(mxCvFc))

      nelem = 0

!  cut each plasma cell into triangles and store the vertices in elem_queue
      do ic = 1, mpg%nCi

        if (mpg%cvVxP(ic,2) == 3) then

!  cell is a triangle, store vertex numbers
          i1 = mpg%cvVx(mpg%cvVxP(ic,1))
          i2 = mpg%cvVx(mpg%cvVxP(ic,1)+1)
          i3 = mpg%cvVx(mpg%cvVxP(ic,1)+2)
          iFc1 = find_face(i1,i2,ic)
          iFc2 = find_face(i2,i3,ic)
          iFc3 = find_face(i3,i1,ic)
          call set_element(i1,i2,i3,ic,iFc1,iFc2,iFc3)
          nelem = nelem + 1

        else

!  cell has more then 3 vertices and needs to be split in triangles
          fcused = .false.
          do iFc = 1, mpg%cvFcP(ic,2)
! face numbers for cell is
            iFcused(iFc) = mpg%cvFc(mpg%cvFcP(ic,1)+iFc-1)
          end do
!  find faces connected with first vertex of cell
          i1 = mpg%cvVx(mpg%cvVxP(ic,1))  ! first vertex of cell ic
          i2 = 0
          i3 = 0
!  we need to identify the two faces connected to vertex i1 of cell ic
!  providing us with vertices i2 and i3 of the first triangle
!  taken from cell ic
         do i = 1, mpg%vxFcP(i1,2)   ! for all faces connected to vertex i1
!           iFc = mpg%vxFc(mpg%vxFcP(i1,1)+i-1)
!           if ((mpg%fcCv(iFc,1) == ic) .or.   ! cell ic is connected
!     .         (mpg%fcCv(iFc,2) == ic)) then  ! to face iFc
!  look for faces of cell ic which are connected to vertex i1
              jfc = 0
              do j = 1, mpg%cvFcP(ic,2) ! for all faces belonging to cell ic
                iFc = iFcused(j)
!               if (.not.fcused(j) .and. (iFcused(j) == iFc)) then
                if (.not.fcused(j)) then
!  check if first vertex i1 of cell ic belongs to face iFc
                  if ((mpg%fcVx(iFc,1) == i1) .or.
     .                (mpg%fcVx(iFc,2) == i1)) then
                    jfc = j   ! face starting at vertex i1 found
                    exit
                  end if
                end if
              end do
              if ((jfc > 0) .and. .not.fcused(jfc)) then
                fcused(jfc) = .true.  ! mark face as used
!  identified the second vertex of face iFc
!  one vertex of face iFc is equal to i1, thus summing up the
!  vertex numbers and subtracting i1 will provide the vertex number we require
                if (i2 == 0) then
                  i2 = mpg%fcVx(iFc,1) + mpg%fcVx(iFc,2) - i1
                  iFc1 = iFc
                else if (i3 == 0) then
                  i3 = mpg%fcVx(iFc,1) + mpg%fcVx(iFc,2) - i1
                  iFc3 = iFc
                  exit
                else
                  call xerrab('write_triangles -- too many faces')
                end if
              else
                call xerrab('write_triangles -- wrong face')
              end if
!           end if
         end do

          iFc2 = 0  ! new face inside cell ic
          call set_element(i1,i2,i3,ic,iFc1,iFc2,iFc3)
          nelem = nelem + 1

!  first triangle finished, now find all others

          do
            ileft = count(.not.fcused)
            if (ileft == 1) exit
!  use first vertex of side 2 of previously found triangle as starting point
            i1 = qlast%iv2
            i2 = qlast%iv3
            iFc1 = qlast%iFc2
            i3 = 0
            do jfc = 1, mpg%cvFcP(ic,2)
              if (fcused(jfc)) cycle
              iFc = iFcused(jfc)
              if ((mpg%fcVx(iFc,1) == i1) .or.
     .            (mpg%fcVx(iFc,2) == i1)) then
!  face needed to form the next triangle found
                i3 = mpg%fcVx(iFc,1) + mpg%fcVx(iFc,2) - i1
                iFc3 = iFc
                fcused(jfc) = .true.
!  only two usable faces left => both belong to triangle
                if (ileft == 2) then
!  identify face number of last triangle side
                  last = pack(iFcused,.not.fcused)
                  iFc2 = last(1)
                end if
                exit
              end if
            end do
            if (i3 == 0)
     .        call xerrab ('write_triangles -- ')

            call set_element(i1,i2,i3,ic,iFc1,iFc2,iFc3)
            nelem = nelem + 1

          end do

!  check that the last unused face of cell ic is part of the
!  last constructed triangle
          unused = 0
          do jfc = 1, mpg%cvFcP(ic,2)
            if (.not.fcused(jfc)) unused = jfc
          end do
          iFcmin = min(mpg%fcVx(iFcused(unused),1),
     .                 mpg%fcVx(iFcused(unused),2))
          iFcmax = max(mpg%fcVx(iFcused(unused),1),
     .                 mpg%fcVx(iFcused(unused),2))
          if (.not.(((iFcmin == min(i1,i2)) .or.
     .               (iFcmin == min(i2,i3)) .or.
     .               (iFcmin == min(i1,i3))) .and.
     .              ((iFcmax == max(i1,i2)) .or.
     .               (iFcmax == max(i2,i3)) .or.
     .               (iFcmax == max(i1,i3)))))
     .      call xerrab ('write_triangles -- wrong face')

        end if

      end do

!  fetch element information from queue and store it in arrays for easier access

      allocate (elem(3,nelem))
      allocate (elem_cv(nelem))
      allocate (neighbor(3,2,nelem))
      allocate (bndmdl(3,nelem))
      allocate (triFc(3,nelem))
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
        elem(1:3,ielem) = (/ qp%iv1, qp%iv2, qp%iv3 /)
        elem_cv(ielem) = qp%ic
        triFc(1:3,ielem) = (/ qp%ifc1, qp%ifc2, qp%ifc3 /)
        isig = 1
        if (triFc(1,ielem) > 0) then
          if (elem_cv(ielem) == mpg%fcCv(triFc(1,ielem),2)) isig(1) = -1
        end if
        if (triFc(2,ielem) > 0) then
          if (elem_cv(ielem) == mpg%fcCv(triFc(2,ielem),2)) isig(2) = -1
        end if
        if (triFc(3,ielem) > 0) then
          if (elem_cv(ielem) == mpg%fcCv(triFc(3,ielem),2)) isig(3) = -1
        end if
        triFc(1:3,ielem) = trifc(1:3,ielem) * isig

        call add_tri(qp%iv1, ielem)
        call add_tri(qp%iv2, ielem)
        call add_tri(qp%iv3, ielem)
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

!  set recycling boundary model indicators
      inquire (file='input.dat',exist=ffound)
      allocate (node_vert(mpg%nVx))
      node_vert = [(i, i = 1, mpg%nVx)]
      if (ffound) then
        call conv_surface_models (mpg,
     .     bndmdl, triFc, elem_cv, nelem, node_vert, mpg%nVx)
      else
        call find_nondef_surfaces (nout, mpg, elem, bndmdl,
     .     elem_cv,  nelem, node_vert, mpg%nVx)
      end if
      deallocate (node_vert)


!  check for holes in the triangle mesh
      nopen = 0
      do it = 1, nelem
        do is = 1, 3
          if (neighbor(is,1,it)+abs(bndmdl(is,it)) == 0) then
            write (*,*) ' is, it ',is,it,' open'
            nopen = nopen + 1
          end if
        end do
      end do
      if (nopen > 0)
     .  call xerrab ('write_triangles -- open mesh')

      call write_npco (mpg%nVx, gm%vxX, gm%vxY, 100._R8)

      call write_elemente (nelem, elem)

      call write_neighbors (nelem, neighbor, bndmdl, elem_cv, triFc)

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

      subroutine set_element(i1, i2, i3, ic, iFc1, ifc2, ifc3)
      integer, intent(in) :: i1, i2, i3, ic, iFc1, ifc2, ifc3
      real(kind=R8) :: at
        allocate(qp)
        nullify(qp%next)
        qp%ic = ic
        qp%iv1 = i1
        qp%iv2 = i2
        qp%iv3 = i3
        qp%iFc1 = ifc1
        qp%iFc2 = ifc2
        qp%iFc3 = ifc3
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
        ifh = qp%iFc1
        qp%iFc1 = qp%ifc3
        qp%iFc3 = ifh
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

      function find_face(i1, i2, ic) result(iFc)
      integer, intent(in) :: i1, i2, ic
      integer :: iFc, j1, j2, i
        do i = 1, mpg%cvFcP(ic,2)
          iFc = mpg%cvFc(mpg%cvFcP(ic,1)+i-1)
          j1 = mpg%fcVx(iFc,1)
          j2 = mpg%fcVx(iFc,2)
          if (((i1 == j1) .and. (i2 == j2)) .or.
     .        ((i1 == j2) .and. (i2 == j1))) return
        end do
        iFc = 0
        return
      end function find_face

      end subroutine write_triangles


      subroutine write_npco (nco, xco, yco, sclfac)

      use b2mod_types

      implicit none
      integer, intent(in) :: nco
      real(R8), intent(in) :: xco(nco), yco(nco), sclfac
      integer :: iun, i
#ifndef F2003
      integer newunit
      external newunit
#endif

#ifdef F2003
      open (newunit=iun,file='triangle.coordinates',
     .           access='SEQUENTIAL',form='FORMATTED')
#else
      iun = newunit()
      open (unit=iun,file='triangle.coordinates',
     .           access='SEQUENTIAL',form='FORMATTED')
#endif
      write (iun,'(i9)') nco
      do i = 1, nco
        write (iun,'(i9,2es24.16)') i, xco(i)*sclfac,
     .                                 yco(i)*sclfac
      end do
      close (iun)

      return
      end subroutine write_npco


      subroutine write_elemente (nel, elem)

      use b2mod_types

      implicit none
      integer, intent(in) :: nel, elem(3,nel)
      integer :: iun, i, is
#ifndef F2003
      integer newunit
      external newunit
#endif

#ifdef F2003
      open (newunit=iun,file='triangle.elements',
     .           access='SEQUENTIAL',form='FORMATTED')
#else
      iun = newunit()
      open (unit=iun,file='triangle.elements',
     .           access='SEQUENTIAL',form='FORMATTED')
#endif
      write (iun,'(i9)') nel
      do i = 1, nel
        write (iun,'(4i9)') i, (elem(is,i), is = 1,3)
      end do
      close (iun)

      return
      end subroutine write_elemente


      subroutine write_neighbors (nel, neighbor, bndmdl, elem_cv, triFc)

      use b2mod_types

      implicit none
      integer, intent(in) :: nel, neighbor(3,2,nel),bndmdl(3,nel),
     .                       elem_cv(nel), triFc(3,nel)
      integer :: iun, i, is, js
#ifndef F2003
      integer newunit
      external newunit
#endif

#ifdef F2003
      open (newunit=iun,file='triangle.neighbors',
     .           access='SEQUENTIAL',form='FORMATTED')
#else
      iun = newunit()
      open (unit=iun,file='triangle.neighbors',
     .           access='SEQUENTIAL',form='FORMATTED')
#endif
      write (iun,'(i9)') nel
      do i = 1, nel
        write (iun,'(i9,3x,3(3i9,3x),4i9)')
     .     i, (neighbor(is,1,i), neighbor(is,2,i),
     .         bndmdl(is,i), is = 1,3), elem_cv(i),
     .        (triFc(js,i), js = 1,3)
      end do
      close (iun)

      return
      end subroutine write_neighbors
