      subroutine prep_triangles (nout, mpg, gm, tria_mode)

      use b2mod_types
      use b2us_map
      use b2us_geo
#ifdef B25_EIRENE
      use eirmod_module_avltree
#endif
      implicit none

      integer, intent(in) :: nout, tria_mode
      type(mapping), intent(in) :: mpg
      type(geometry), intent(in) :: gm

      ! data for the (new/converted) triangles
      integer :: nnodes, npnodes, nelem, npelem
      integer, allocatable :: elem(:,:), neighbor(:,:,:), bndmdl(:,:),
     .    elem_cv(:), trifc(:,:), node_vert(:), vert_node(:), cntr(:,:)
      real(R8), allocatable :: xco(:), yco(:)

      ! aux variables
      integer :: ivx, is, it, iel, iel2, nmxnodes, nmxelem
      integer :: ico, io_33, io_34, io_35, iused, nopen
      logical linq_33, linq_34, linq_35, lcheck, ex

      ! check for presence of consistent fort.3[3-5] files
      ! for tria_mode = 0: triangle files to be converted (void + plasma)
      ! for tria_mode = 1: assume void triangles only; plasma triangles to be added
      linq_33 = .false.
      linq_34 = .false.
      linq_35 = .false.
      ico  = 0
      iel  = 0
      iel2 = 0
      inquire(unit=33,exist=ex)
      if (ex) then
        open (unit=33)
        read (33,*,iostat=io_33) ico
        linq_33 = (io_33 == 0)
        close (33)
      endif
      inquire(unit=34,exist=ex)
      if (ex) then
        open (unit=34)
        read (34,*,iostat=io_34) iel
        linq_34 = (io_34 == 0)
        close (34)
      endif
      inquire(unit=35,exist=ex)
      if (ex) then
        open (unit=35)
        read (35,*,iostat=io_35) iel2
        linq_35 = (io_35 == 0)
        close (35)
      endif
      lcheck = (ico > 0) .and. (iel > 0) .and. (iel2 > 0) .and.
     .         (iel == iel2)
      if (.not.lcheck) then
        if (tria_mode.eq.0) write (*,*) 
     .    'No valid fort.3[3-5] files found for conversion. '//
     .    'Creating new triangles instead (plasma domain only).'
        if (tria_mode.eq.1) write (*,*) 
     .    'No fort.3[3-5] files with void triangles found. '//
     .    'Only creating triangles in plasma domain.'
        ico  = 0
        iel  = 0
        iel2 = 0
      endif
      ! Pre-allocate the required arrays, with some margin
      nmxnodes = ico + mpg%nVx
      nmxelem  = iel + 2*mpg%nCv
      allocate (xco(nmxnodes))
      allocate (yco(nmxnodes))
      allocate (node_vert(nmxnodes))
      allocate (vert_node(mpg%nVx))
      allocate (elem(3,nmxelem))
      allocate (neighbor(3,2,nmxelem))
      allocate (bndmdl(3,nmxelem))
      allocate (elem_cv(nmxelem))
      allocate (trifc(3,nmxelem))
      allocate (cntr(3,nmxelem))
      xco = 0.0_R8
      yco = 0.0_R8
      node_vert = 0
      vert_node = 0
      elem = 0
      neighbor = 0
      bndmdl = 0
      elem_cv = 0
      trifc = 0
      cntr = 0

      ! Start building node list by adding plasma vertices
      npnodes = 0
      npelem  = 0
      do ivx = 1, mpg%nVx
        !! Check whether vertex is actually used in the plasma grid
        !! => remove X-point duplicates
        iused = count(mpg%cvVx == ivx)
        if (iused > 0) then
          npnodes = npnodes + 1
          xco(npnodes) = gm%vxX(ivx)*100._R8    ! convert to cm
          yco(npnodes) = gm%vxY(ivx)*100._R8    ! convert to cm
          node_vert(npnodes) = ivx              ! store plasma vertex number corresponding to the node
          vert_node(ivx) = npnodes
        else
          write (*,'(a,i6,a,1p,1e15.8,a,1e15.8,a/a)')
     .     'Vertex ', ivx ,' with coordinates (x, y) = (',
     .      gm%vxX(ivx),',',gm%vxY(ivx),') not used. ',
     .     'Likely X-point duplicate; not included among '//
     .     'triangle nodes.'
        endif
      end do
      write (*,*) 'Found ', npnodes, ' unique plasma vertices/nodes.'
      nnodes = npnodes
      nelem  = npelem

      if (tria_mode.eq.1 .or. .not.lcheck) then
        ! Make triangles in the plasma domain
        write (*,*) 'Building triangles for plasma region'
        call make_plasma_triangles(mpg, gm, nmxelem, nelem,
     .    vert_node, elem, elem_cv, trifc, neighbor)
        npelem = nelem
      endif

      if (lcheck) then
        ! Read and convert triangles from existing files
	! tria_mode 0: plasma + void
        ! tria_mode 1: void only
        if (tria_mode.eq.0) 
     .    write (*,*) 'Converting triangles'
        if (tria_mode.eq.1) 
     .    write (*,*) 'Converting triangles for void region'
        call conv_triangles(mpg, gm,
     .    nmxnodes, nmxelem, nnodes, nelem, npnodes, npelem,
     .    xco, yco, node_vert, elem, elem_cv, trifc, neighbor,
     .    cntr, tria_mode)
      endif

      ! Set recycling boundary model indicators
      inquire (file='input.dat',exist=ex)
      if (ex) then
        write (*,*) 'Converting surface models from input.dat'
        call conv_surface_models (mpg,
     .     bndmdl, trifc, elem_cv, nelem, node_vert, nnodes)
      else
        write (*,*) 'Identifying nondef. standard surfaces'
        call find_nondef_surfaces (nout, mpg,
     .     elem, bndmdl, elem_cv, nelem, node_vert, nnodes, neighbor)
      end if

      do it = 1, nelem
        do is = 1, 3
          write (*,*) 'Tria, side, neigh, bndmdl, cntr '
          write (*,*) it,is,neighbor(is,1,it),bndmdl(is,it),cntr(is,it)
        enddo
      enddo

      ! Check for holes in the triangle mesh
      nopen = 0
      do it = 1, nelem
        do is = 1, 3
          if (neighbor(is,1,it) + abs(bndmdl(is,it))
     .                          + abs(cntr(is,it)) == 0) then
            write (*,*) 'Triangle side is, it ',is,it,' open. '//
     .                  'Node coordinates: '
            write (*,*) xco(elem(1,it)),xco(elem(2,it)),xco(elem(3,it))
            write (*,*) yco(elem(1,it)),yco(elem(2,it)),yco(elem(3,it))
            write (*,*) 'Neighbor, bndmdl, cntr = ',neighbor(is,1,it),
     .        bndmdl(is,it), cntr(is,it)
            nopen = nopen + 1
          end if
        end do
      end do
      if (nopen > 0)
     .  call xerrab ('write_triangles -- open mesh')

      ! write the output files
      call write_triangles (nnodes, nelem, xco, yco, elem,
     .    neighbor, bndmdl, elem_cv, trifc)


      deallocate (xco)
      deallocate (yco)
      deallocate (node_vert)
      deallocate (vert_node)
      deallocate (elem)
      deallocate (neighbor)
      deallocate (bndmdl)
      deallocate (elem_cv)
      deallocate (trifc)
      deallocate (cntr)

      end subroutine prep_triangles
