      module b2us_prep

      use b2mod_types
      use b2mod_geo
      use b2mod_geo_corner
      use b2mod_indirect
      use b2mod_plasma
      use b2mod_transport
      use b2mod_anomalous_transport
      use b2mod_work
      use b2mod_residuals
      use b2mod_input_profile
      use b2mod_rates
      use b2mod_diag
      use b2us_map
      use b2us_geo
      use b2us_plasma

      implicit none

      private

      public :: read_mapping, calc_mapping, write_mapping,
     .          read_geometry, prep_geometry, write_geometry,
     .          calc_mapping_bc, calc_mapping_neut,
     .          prep_plasma, 
     .          prep_plasma_ext

cwdk  eventually: take this outside of a module! (or: make sure this module is
cwdk  only used at b2mod_driver level)
      type(mapping), public, save :: mpg
      type(geometry), public, save :: gm
      type(B2State), public, save :: state
      type(B2StateExt), public, save :: state_ext

      contains

***********************************************************************

      subroutine read_mapping (fname)

      implicit none
      character(*),intent(in) :: fname
      integer :: iun, io, idum(0:9)

      open (newunit=iun,file=fname,status='old',iostat=io)

      if (io /= 0) then
        write (*,*) ' file ', fname, ' not readable'
        mpg%nCv = 0
        mpg%nFc = 0
        mpg%nVx = 0
      else
        call cfruin (iun, 1, idum, 'nCi')
        mpg%nCi = idum(0) 
        call cfruin (iun, 1, idum, 'nCg')
        mpg%nCg = idum(0) 
        call cfruin (iun, 1, idum, 'nCv')
        mpg%nCv = idum(0)
        call cfruin (iun, 1, idum, 'nFc')
        mpg%nFc = idum(0)
        call cfruin (iun, 1, idum, 'nVx')
        mpg%nVx = idum(0)
        call cfruin (iun, 1, idum, 'nFs')
        mpg%nFs = idum(0)
        call cfruin (iun, 1, idum, 'nCmxVx')
        mpg%nCmxVx = idum(0)
        call cfruin (iun, 1, idum, 'nCmxFc')
        mpg%nCmxFc = idum(0)
        call cfruin (iun, 1, idum, 'mxnFs')
        mpg%mxnFs = idum(0)

        call alloc_mapping (mpg)
        
        call cfruin (iun, mpg%nCv*2,          mpg%cvFcP, 'cvFcP')
        call cfruin (iun, mpg%nCmxFc*mpg%nCv, mpg%cvFc,  'cvFc')
        call cfruin (iun, mpg%nFc*2,          mpg%fcCv,  'fcCv')
        call cfruin (iun, mpg%nFc*2,          mpg%fcVx,  'fcVx')
        call cfruin (iun, mpg%nCv*2,          mpg%cvVxP, 'cvVxP')
        call cfruin (iun, mpg%nCmxVx*mpg%nCv, mpg%cvVx,  'cvVx')
        call cfruin (iun, mpg%nVx*2,          mpg%vxFcP, 'vxFcP')
        call cfruin (iun, mpg%nCmxFc*mpg%nVx, mpg%vxFc,  'vxFc')
        call cfruin (iun, mpg%nVx*2,          mpg%vxCvP, 'vxCvP')
        call cfruin (iun, mpg%nCmxVx*mpg%nCv, mpg%vxCv,  'vxCv')
        call cfruin (iun, mpg%mxnFs*2,        mpg%cvFsP, 'cvFsP')
        call cfruin (iun, mpg%nCv,            mpg%cvFs,  'cvFs')
        call cfruin (iun, mpg%mxnFs*2,        mpg%fcFsP, 'fcFsP')
        call cfruin (iun, mpg%nCv,            mpg%fcFs,  'fcFs')
        call cfruin (iun, mpg%nFc,            mpg%fcReg,  'fcReg')
        call cfruin (iun, mpg%nCv,            mpg%cvReg,  'cvReg')
        call cfruin (iun, mpg%nCv,            mpg%cvOnClosedSurface,
     .                                            'cvOnClosedSurface')
      end if

      close  (iun)

      return
      end subroutine read_mapping 	

***********************************************************************

      subroutine write_mapping (fname)

      implicit none
      character(*),intent(in) :: fname
      integer :: nout, io, idum(0:9)

      open (newunit=nout,file=fname,iostat=io)

      if (io /= 0) then
        write (*,*) ' file ', fname, ' cannot be written'
      else

        idum(0) = mpg%nCi
        call cfwuin (nout, 1, idum, 'nCi')
        idum(0) = mpg%nCg
        call cfwuin (nout, 1, idum, 'nCg')
        idum(0) = mpg%nCv
        call cfwuin (nout, 1, idum, 'nCv')
        idum(0) = mpg%nFc
        call cfwuin (nout, 1, idum, 'nFc')
        idum(0) = mpg%nVx
        call cfwuin (nout, 1, idum, 'nVx')
        idum(0) = mpg%nFs
        call cfwuin (nout, 1, idum, 'nFs')
        idum(0) = mpg%nCmxVx
        call cfwuin (nout, 1, idum, 'nCmxVx')
        idum(0) = mpg%nCmxFc
        call cfwuin (nout, 1, idum, 'nCmxFc')
        idum(0) = mpg%mxnFs
        call cfwuin (nout, 1, idum, 'mxnFs')

        call cfwuin (nout, mpg%nCv*2,          mpg%cvFcP, 'cvFcP')
        call cfwuin (nout, mpg%nCmxFc*mpg%nCv, mpg%cvFc,  'cvFc')
        call cfwuin (nout, mpg%nFc*2,          mpg%fcCv,  'fcCv')
        call cfwuin (nout, mpg%nFc*2,          mpg%fcVx,  'fcVx')
        call cfwuin (nout, mpg%nCv*2,          mpg%cvVxP, 'cvVxP')
        call cfwuin (nout, mpg%nCmxVx*mpg%nCv, mpg%cvVx,  'cvVx')
        call cfwuin (nout, mpg%nVx*2,          mpg%vxFcP, 'vxFcP')
        call cfwuin (nout, mpg%nCmxFc*mpg%nVx, mpg%vxFc,  'vxFc')
        call cfwuin (nout, mpg%nVx*2,          mpg%vxCvP, 'vxCvP')
        call cfwuin (nout, mpg%nCmxVx*mpg%nCv, mpg%vxCv,  'vxCv')
        call cfwuin (nout, mpg%mxnFs*2,        mpg%cvFsP, 'cvFsP')
        call cfwuin (nout, mpg%nCv,            mpg%cvFs,  'cvFs')
        call cfwuin (nout, mpg%mxnFs*2,        mpg%fcFsP, 'fcFsP')
        call cfwuin (nout, mpg%nCv,            mpg%fcFs,  'fcFs')
        call cfwuin (nout, mpg%nFc,            mpg%fcReg,  'fcReg')
        call cfwuin (nout, mpg%nCv,            mpg%cvReg,  'cvReg')
        call cfwuin (nout, mpg%nCv,            mpg%cvOnClosedSurface,
     .                                            'cvOnClosedSurface')

      end if

      close  (nout)

      return
      end subroutine write_mapping 	

***********************************************************************

      subroutine calc_mapping(nx,ny,
     .           leftix,leftiy,rightix,rightiy,
     .           topix,topiy,bottomix,bottomiy,
     .           nncut,leftcut,rightcut,topcut,bottomcut)

      implicit none
      integer, intent(in) :: nx, ny, nncut
      integer, intent(in) :: leftix(-1:nx,-1:ny),
     .                       leftiy(-1:nx,-1:ny),
     .                       rightix(-1:nx,-1:ny),
     .                       rightiy(-1:nx,-1:ny),
     .                       topix(-1:nx,-1:ny),
     .                       topiy(-1:nx,-1:ny),
     .                       bottomix(-1:nx,-1:ny),
     .                       bottomiy(-1:nx,-1:ny)
      integer, intent(in) :: leftcut(nncut), rightcut(nncut),
     .                       topcut(nncut), bottomcut(nncut)
      integer :: ix, iy, i, if1, if2, if3, if4, ic, ifc, 
     .           vxFc_end, vxCv_end, icut, ix1, ix2, iy1, iy2,
     .           ivx, ivx1, ivx2, nxx, mxreg, ifs, icv, ireg, istart
      integer, allocatable :: indFc(:), indCv(:), fh(:), ch(:)

!  assume a structured quadrilateral grid

!  no. of cells
      mpg%nCi = nx*ny                         ! internal cells
      mpg%nCg = 2*nx + 2*ny                   ! guard cells
      mpg%nCv = mpg%nCi + mpg%nCg             ! total number of cells

!  no. of cell faces
      mpg%nFc = ny*(nx+1) + nx*(ny+1)

! no. of vertices
      mpg%nVx = (nx+1) * (ny+1)

! max. no. of vertices per cell
      mpg%nCmxVx = 4

! max. no. of faces per cell
      mpg%nCmxFc = 4

! max no. of flux surfaces
      mpg%mxnFs = 8*nncut*ny

      call alloc_mapping (mpg)

! no. of boundaries to be considered for this cell
      mpg%cvFcP(1:mpg%nCi,2) = 4
      mpg%cvFcP(mpg%nCi+1:mpg%nCv,2) = 1

! no. of vertices to be considered for this cell
      mpg%cvVxP(1:mpg%nCi,2) = 4
      mpg%cvVxP(mpg%nCi+1:mpg%nCv,2) = 2
      
! starting index in the list of cell faces
! starting index for every cell in the cvVx list and its number of vertices
      mpg%cvFcP(1,1) = 1
      mpg%cvVxP(1,1) = 1
      do i=2, mpg%nCv
        mpg%cvFcP(i,1) = mpg%cvFcP(i-1,1) + mpg%cvFcP(i-1,2)
        mpg%cvVxP(i,1) = mpg%cvVxP(i-1,1) + mpg%cvVxP(i-1,2)
      end do

! list of cell faces and vertices for inner cells
! interpolation coefficients for interpolation to cell centers
      nxx = (nx+1)*ny
      do iy = 1, ny
        do ix = 1, nx
          ic = ix + (iy-1)*nx
          if1 = mpg%cvFcP(ic,1)
          if2 = if1 + 1
          if3 = if2 + 1
          if4 = if3 + 1
!          mpg%cvFc(if1) = ix + (iy-1)*(nx+1)
          mpg%cvFc(if1) = leftix(ix-1,iy-1)+2 + (iy-1)*(nx+1)
          mpg%cvFc(if2) = ix+1 + (iy-1)*(nx+1)
          mpg%cvFc(if3) = nxx + iy + (ix-1)*(ny+1)
          mpg%cvFc(if4) = nxx + iy+1 + (ix-1)*(ny+1)
          mpg%intcellP(if1) = 1.0_R8
          mpg%intcellP(if2) = 1.0_R8
          mpg%intcellP(if3) = 0.0_R8
          mpg%intcellP(if4) = 0.0_R8
          mpg%intcellR(if1) = 0.0_R8
          mpg%intcellR(if2) = 0.0_R8
          mpg%intcellR(if3) = 1.0_R8
          mpg%intcellR(if4) = 1.0_R8
!          mpg%cvVx(if1) = ix + (iy-1)*(nx+1)
          mpg%cvVx(if1) = leftix(ix-1,iy-1)+2 + (iy-1)*(nx+1)
          mpg%cvVx(if2) = ix+1 + (iy-1)*(nx+1)
!          mpg%cvVx(if3) = ix + iy*(nx+1)
          mpg%cvVx(if3) = leftix(ix-1,iy-1)+2 + iy*(nx+1)
          mpg%cvVx(if4) = ix+1 + iy*(nx+1)
        end do
      end do

! list of cell faces and vertices for guard cells
      do ix = 1, nx
        ic = mpg%nCi + ix
        mpg%cvFc(mpg%cvFcP(ic,1)) = nxx + 1 + (ix-1)*(ny+1)
        mpg%intcellP(mpg%cvFcP(ic,1)) = 1.0_R8
        mpg%intcellR(mpg%cvFcP(ic,1)) = 1.0_R8
        mpg%cvVx(mpg%cvVxP(ic,1)) = ix 
        mpg%cvVx(mpg%cvVxP(ic,1)+1) = ix+1 
        ic = mpg%nCi + nx + ix
        mpg%cvFc(mpg%cvFcP(ic,1)) = nxx + ny+1 + (ix-1)*(ny+1)
        mpg%intcellP(mpg%cvFcP(ic,1)) = 1.0_R8
        mpg%intcellR(mpg%cvFcP(ic,1)) = 1.0_R8
        mpg%cvVx(mpg%cvVxP(ic,1)) = ix + ny*(nx+1)
        mpg%cvVx(mpg%cvVxP(ic,1)+1) = ix+1 + ny*(nx+1)
      end do

      do iy = 1, ny
        ic = mpg%nCi + 2*nx + iy
        mpg%cvFc(mpg%cvFcP(ic,1)) = 1 + (iy-1)*(nx+1)
        mpg%intcellP(mpg%cvFcP(ic,1)) = 1.0_R8
        mpg%intcellR(mpg%cvFcP(ic,1)) = 1.0_R8
        mpg%cvVx(mpg%cvVxP(ic,1)) = 1 + (iy-1)*(nx+1) 
        mpg%cvVx(mpg%cvVxP(ic,1)+1) = 1 + iy*(nx+1) 
        ic = mpg%nCi + 2*nx + ny + iy
        mpg%cvFc(mpg%cvFcP(ic,1)) = nx+1 + (iy-1)*(nx+1)
        mpg%intcellP(mpg%cvFcP(ic,1)) = 1.0_R8
        mpg%intcellR(mpg%cvFcP(ic,1)) = 1.0_R8
        mpg%cvVx(mpg%cvVxP(ic,1)) = nx+1 + (iy-1)*(nx+1) 
        mpg%cvVx(mpg%cvVxP(ic,1)+1) = nx+1 + iy*(nx+1) 
      end do

! listing for every face the two corresponding cells
! and the corresponding vertices
! x faces
       do iy = 1, ny
         do ix = 2, nx
           ic = ix + (iy-1)*nx
           ifc = ix + (iy-1)*(nx+1)
           mpg%fcCv(ifc,1) = ic - 1
           mpg%fcCv(ifc,2) = rightix(ix-2,iy-1)+1 + (iy-1)*nx
           mpg%fcNv(ifc,1) = 6
           mpg%fcNv(ifc,2) = 5
           mpg%fcVx(ifc,1) = ix + (iy-1)*(nx+1)
           mpg%fcVx(ifc,2) = ix + iy*(nx+1)
         end do

! left guard cells
         ic = 1 + (iy-1)*nx
         ifc = 1 + (iy-1)*(nx+1)
         mpg%fcCv(ifc,1) = mpg%nCi + 2*nx + iy
         mpg%fcCv(ifc,2) = ic
         mpg%fcNv(ifc,1) = 4
         mpg%fcNv(ifc,2) = 5
         mpg%fcVx(ifc,1) = 1 + (iy-1)*(nx+1)
         mpg%fcVx(ifc,2) = 1 + iy*(nx+1)

! right guard cells
         ic = nx+1 + (iy-1)*nx
         ifc = nx+1 + (iy-1)*(nx+1)
         mpg%fcCv(ifc,1) = ic - 1
         mpg%fcCv(ifc,2) = mpg%nCi + 2*nx + ny + iy
         mpg%fcNv(ifc,1) = 6
         mpg%fcNv(ifc,2) = 4
         mpg%fcVx(ifc,1) = nx+1 + (iy-1)*(nx+1)
         mpg%fcVx(ifc,2) = nx+1 + iy*(nx+1)
       end do

! y faces
       do ix = 1, nx
         do iy = 2, ny
           ic = ix + (iy-1)*nx
           ifc = nxx + iy + (ix-1)*(ny+1)
           mpg%fcCv(ifc,2) = ic
           mpg%fcCv(ifc,1) = ic - nx
           mpg%fcNv(ifc,1) = 8
           mpg%fcNv(ifc,2) = 3
           mpg%fcVx(ifc,2) = ix + (iy-1)*(nx+1)
           mpg%fcVx(ifc,1) = ix+1 + (iy-1)*(nx+1)
         end do

! lower guard cells
         ic = ix 
         ifc = nxx + ix
         mpg%fcCv(ifc,2) = ic
         mpg%fcCv(ifc,1) = mpg%nCi + ix
         mpg%fcNv(ifc,1) = 5
         mpg%fcNv(ifc,2) = 3
         mpg%fcVx(ifc,2) = ix 
         mpg%fcVx(ifc,1) = ix+1 

! upper guard cells
         ic = ix + (ny-1)*(nx+1)
         ifc = nxx + nx + ix
         mpg%fcCv(ifc,2) = mpg%nCi + nx + ix
         mpg%fcCv(ifc,1) = ic
         mpg%fcNv(ifc,1) = 8
         mpg%fcNv(ifc,2) = 3
         mpg%fcVx(ifc,2) = ix + ny*(nx+1)
         mpg%fcVx(ifc,1) = ix+1 + ny*(nx+1)
       end do

! merging x-point vertex
       do icut = 1, nncut
         ix1 = leftcut(icut)
         iy1 = topcut(icut)
         ix2 = rightcut(icut)
         iy2 = topcut(icut)
         ivx1 = ix1+1 + iy1*(nx+1)
         ivx2 = ix2+1 + iy2*(nx+1)

         where (mpg%cvVx == ivx2)
           mpg%cvVx = ivx1
         end where

         where (mpg%fcVx == ivx2)
           mpg%fcVx = ivx1
         end where

       end do

! listing for every vertex the cell and face numbers
! and the starting points in the list
       allocate (indFc(mpg%nFc))
       indFc = (/ (i, i=1,mpg%nFc) /)

       allocate (indCv(mpg%nCv))
       indCv = (/ (i, i=1,mpg%nCv) /)

       allocate (fh(mpg%nFc))
       fh(1:mpg%nFc) = mpg%fcVx(1:mpg%nFc,1)

       allocate (ch(mpg%ncV))
       ch(1:mpg%ncV) = mpg%cvVx(1:mpg%ncV)
       
       mpg%vxFcP(1,1) = 1
       mpg%vxCvP(1,1) = 1
       do ivx = 1, (nx+1)*(ny+1)-1
         mpg%vxFcP(ivx,2) = count(fh == ivx)
         mpg%vxFcP(ivx+1,1) = mpg%vxFcP(ivx,1) + mpg%vxFcP(ivx,2)  
         if (mpg%vxFcP(ivx,2) > 0) then
           mpg%vxFc(mpg%vxFcP(ivx,1):mpg%vxFcP(ivx+1,1)-1) =
     .          pack( indFc, fh == ivx )
         end if
         mpg%vxCvP(ivx,2) = count(ch == ivx)
         mpg%vxCvP(ivx+1,1) = mpg%vxCvP(ivx,1) + mpg%vxCvP(ivx,2)  
         if (mpg%vxCvP(ivx,2) > 0) then
           mpg%vxCv(mpg%vxCvP(ivx,1):mpg%vxCvP(ivx+1,1)-1) =
     .          pack( indCv, ch == ivx )
         end if
       end do

       ivx = (nx+1)*(ny+1)
       mpg%vxFcP(ivx,2) = count(fh == ivx)
       vxFc_end = mpg%vxFcP(ivx,1) + mpg%vxFcP(ivx,2)
       if (mpg%vxFcP(ivx,2) > 0) then
         mpg%vxFc(mpg%vxFcP(ivx,1):vxFc_end) =
     .          pack( indFc, fh == ivx )
       end if
       mpg%vxCvP(ivx,2) = count(ch == ivx)
       vxCv_end = mpg%vxCvP(ivx,1) + mpg%vxCvP(ivx,2)
       if (mpg%vxCvP(ivx,2) > 0) then
         mpg%vxCv(mpg%vxCvP(ivx,1):vxCv_end) =
     .          pack( indCv, ch == ivx )
       end if

       deallocate (indFc)
       deallocate (indCv)
       deallocate (fh)
       deallocate (ch)

! regions (cell volumes)
      do iy = 0, ny-1
        do ix = 0, nx-1
          ic = (ix+1) + iy*nx
          mpg%cvReg(ic) = region(ix,iy,0)
          mpg%cvOnClosedSurface(ic) = on_closed_surface(ix,iy)
	end do
      end do	

      do ix = 0, nx-1
        ic = mpg%nCi + ix+1
        mpg%cvReg(ic) = region(ix,-1,0)
        mpg%cvOnClosedSurface(ic) = on_closed_surface(ix,-1)

        ic = mpg%nCi + nx + ix+1
        mpg%cvReg(ic) = region(ix,ny,0)
        mpg%cvOnClosedSurface(ic) = on_closed_surface(ix,ny)
      end do
 
      do iy = 0, ny-1
        ic = mpg%nCi + 2*nx + iy+1
        mpg%cvReg(ic) = region(-1,iy,0)
        mpg%cvOnClosedSurface(ic) = on_closed_surface(-1,iy)

        ic = mpg%nCi + 2*nx + ny + iy+1
        mpg%cvReg(ic) = region(nx,iy,0)
        mpg%cvOnClosedSurface(ic) = on_closed_surface(nx,iy)
      end do

! regions (faces)
 ! x faces
       do iy = 0, ny-1
         do ix = 1, nx-1
           ifc = ix+1 + iy*(nx+1)
           mpg%fcReg(ifc) = region(ix,iy,1)
         end do

! left boundary face (to guard cells)
         ix = 0
         ifc = ix+1 + iy*(nx+1)
         mpg%fcReg(ifc) = region(ix,iy,1)

! right guard cells
         ix = nx
         ifc = ix+1 + iy*(nx+1)
         mpg%fcReg(ifc) = region(ix,iy,1)
        end do

! y faces
       do ix = 0, nx-1
         do iy = 1, ny-1
           ifc = nxx + iy+1 + ix*(ny+1)
           mpg%fcReg(ifc) = region(ix,iy,2)
         end do

! lower guard cells
         iy = 0
         ifc = nxx + ix+1
         mpg%fcReg(ifc) = region(ix,iy,2)

! upper guard cells
         iy = ny
         ifc = nxx + nx + ix+1
         mpg%fcReg(ifc) = region(ix,iy,2)
      end do

      mxReg = maxval(mpg%cvReg(1:mpg%nCv))

      ifs = 0
      icv = 0
      do ireg= 1, mxReg
	do iy = 0, ny-1
	  istart = -1
	  do ix = 0, nx-1
            ic = (ix+1) + iy*nx
            ifc = ix + iy*(nx+1)
            if (mpg%cvReg(ic) == ireg) then
	      if (istart == -1) then
                ifs = ifs + 1
	        icv = icv + 1
                mpg%cvFs(icv) = ic
	        mpg%cvFsP(ifs,1) = icv
	        mpg%fcFs(icv) = ifc
                mpg%fcFsP(ifs,1) = icv
	        istart = 1
              else
	        icv = icv + 1
		mpg%cvFs(icv) = ic
	        mpg%fcFs(icv) = ifc
              end if
	    else if (istart /= -1) then
	      mpg%cvFsP(ifs,2) = icv - mpg%cvFsP(ifs,1) + 1
	      mpg%fcFsP(ifs,2) = icv - mpg%fcFsP(ifs,1) + 1
            end if
          end do
        end do	
      end do
      mpg%nFs = ifs

      end subroutine calc_mapping

***********************************************************************

      subroutine calc_mapping_bc (nx, ny, ns, m)
      use b2mod_boundary_namelist
      implicit none
      integer :: nx, ny, ns, ib, ibc, is, inbc
      integer :: bcs(nbc), bce(nbc), bcl(nbc)
      type (mapping), intent (inout) :: m



      m%nBc   = nbc     ! number of boundary conditions
      m%mxnBc = 0       ! number of cells with imposed BC
      do ib = 1, nbc
        bcs(ib) = bcstart(ib)
        bce(ib) = bcend(ib)
        ! exclude corner cells if present
        if (bcchar(ib).eq.'N' .or. bcchar(ib).eq.'S') then
          if (leftix (bcstart(ib),bcpos(ib)) .eq. -2)  
     >      bcs(ib) = bcs(ib) + 1
          if (rightix(bcend(ib),  bcpos(ib)) .eq. nx+1) 
     >      bce(ib) = bce(ib) - 1
        endif
        if (bcchar(ib).eq.'E' .or. bcchar(ib).eq.'W') then
          if (bottomiy(bcpos(ib),bcstart(ib)) .eq. -2)
     >      bcs(ib) = bcs(ib) + 1
          if (topiy   (bcpos(ib),bcend(ib)  ) .eq. ny+1)
     >      bce(ib) = bce(ib) - 1
        endif
        bcl(ib) = bce(ib) - bcs(ib) + 1
        m%mxnBc = m%mxnBc + bcl(ib)
      enddo
      if (m%mxnBc .ne. m%nCg ) 
     >  call xerrab ('calc_mapping_bc -- wrong mxnBc')


      call alloc_mapping_bc (m)

      inbc = 0 ! index in bcCv
      do ib = 1, nbc
        mpg%bcCvP(ib,1) = inbc + 1
        mpg%bcCvP(ib,2) = bcl(ib)
        if (bcchar(ib).eq.'S') then
          do ibc = 1, bcl(ib)
            inbc = inbc + 1
            mpg%bcCv(inbc,1) = mpg%nCi + bcs(ib) + ibc ! number of the guard cell
            mpg%bcCv(inbc,2) = bcs(ib) + ibc           ! number of the corresponding domain cell
            mpg%bcFcOr(inbc) = -1.0_R8
          enddo
        elseif (bcchar(ib).eq.'N') then
          do ibc = 1, bcl(ib)
            inbc = inbc + 1
            mpg%bcCv(inbc,1) = mpg%nCi + nx + bcs(ib) + ibc
            mpg%bcCv(inbc,2) = nx*(ny-1) + bcs(ib) + ibc 
            mpg%bcFcOr(inbc) = 1.0_R8
          enddo
        elseif (bcchar(ib).eq.'W') then
          do ibc = 1, bcl(ib)
            inbc = inbc + 1
            mpg%bcCv(inbc,1) = mpg%nCi + 2*nx + bcs(ib) + ibc
            mpg%bcCv(inbc,2) = nx*bcs(ib) + 1
            mpg%bcFcOr(inbc) = -1.0_R8
          enddo
        elseif (bcchar(ib).eq.'E') then 
          do ibc = 1, bcl(ib)
            inbc = inbc + 1
            mpg%bcCv(inbc,1) = mpg%nCi + 2*nx + ny + bcs(ib) + ibc
            mpg%bcCv(inbc,2) = nx*(bcs(ib) + 1)
            mpg%bcFcOr(inbc) = 1.0_R8
          enddo
        else
          call xerrab ('calc_mapping_bc -- wrong BCCHAR')
        endif
      enddo     

      if (inbc .ne. m%nCg ) 
     >  call xerrab ('calc_mapping_bc -- problem setting up bc arrays')

      end subroutine calc_mapping_bc


***********************************************************************



      subroutine calc_mapping_neut (nx, ny, ns, m)
      use b2mod_neutrals_namelist
      implicit none
      integer :: nx, ny, ns
      type (mapping), intent (inout) :: m

      call xerrab('to be written!')

      end subroutine calc_mapping_neut


***********************************************************************


      subroutine read_geometry (fname)

      implicit none
      character(*),intent(in) :: fname
      integer :: iun, io, idum(0:9)

      open (newunit=iun,file=fname,iostat=io)

      if (io /= 0) then
        write (*,*) ' file ', fname, ' cannot be written'
      else

        call alloc_geometry (gm, mpg%nCv, mpg%nFc, mpg%nVx, mpg%nFs)

! cell volumes 
        call cfrure (iun, mpg%nCv*4, gm%cvBb,   'cvBb')
        call cfrure (iun, mpg%nCv,   gm%cvX,    'cvX')
        call cfrure (iun, mpg%nCv,   gm%cvY,    'cvY')
        call cfrure (iun, mpg%nCv,   gm%cvS,    'cvS')
        call cfrure (iun, mpg%nCv,   gm%cvH,    'cvH')
        call cfrure (iun, mpg%nCv*2, gm%cvQgam, 'cvQgam')
        call cfrure (iun, mpg%nCv,   gm%cvVol,  'cvVol')
        call cfrure (iun, mpg%nCv,   gm%cvOnedBsq,'cvOnedBsq')
        call cfrure (iun, mpg%nCv,   gm%cvBzb,  'cvBzb')

! face quantities
        call cfrure (iun, mpg%nFc*4, gm%fcBb,   'fcBb')
        call cfrure (iun, mpg%nFc,   gm%fcS,    'fcS')
        call cfrure (iun, mpg%nFc*2, gm%fcHc,   'fcHc')
        call cfrure (iun, mpg%nFc,   gm%fcHt,   'fcHt')
        call cfrure (iun, mpg%nFc*2, gm%fcQgam, 'fcQgam')
        call cfrure (iun, mpg%nFc*2, gm%fcQalf, 'fcQalf')
        call cfrure (iun, mpg%nFc*2, gm%fcQbet, 'fcQbet')
        call cfrure (iun, mpg%nFc,   gm%fcPbs,  'fcPbs')
        call cfrure (iun, mpg%nFc,   gm%fcPbshz,'fcPbshz')

! vertice quantities
        call cfrure (iun, mpg%nVx,   gm%vxX,    'vxX')
        call cfrure (iun, mpg%nVx,   gm%vxY,    'vxY')
        call cfrure (iun, mpg%nVx,   gm%vxFfbz, 'vxFfbz')
        call cfrure (iun, mpg%nVx,   gm%vxFpsi, 'vxFpsi')

        call cfrure (iun, mpg%nFs,   gm%fsConn, 'fsConn')

      end if

      close  (iun)

      return
      end subroutine read_geometry

***********************************************************************

      subroutine write_geometry (fname)

      implicit none
      character(*),intent(in) :: fname
      integer :: nout, io, idum(0:9)

      open (newunit=nout,file=fname,iostat=io)

      if (io /= 0) then
        write (*,*) ' file ', fname, ' cannot be written'
      else
! cell volumes 
        call cfwure (nout, mpg%nCv*4, gm%cvBb,   'cvBb')
        call cfwure (nout, mpg%nCv,   gm%cvX,    'cvX')
        call cfwure (nout, mpg%nCv,   gm%cvY,    'cvY')
        call cfwure (nout, mpg%nCv,   gm%cvS,    'cvS')
        call cfwure (nout, mpg%nCv,   gm%cvH,    'cvH')
        call cfwure (nout, mpg%nCv*2, gm%cvQgam, 'cvQgam')
        call cfwure (nout, mpg%nCv,   gm%cvVol,  'cvVol')
        call cfwure (nout, mpg%nCv,   gm%cvOnedBsq,'cvOnedBsq')
        call cfwure (nout, mpg%nCv,   gm%cvBzb,  'cvBzb')

! face quantities
        call cfwure (nout, mpg%nFc*4, gm%fcBb,   'fcBb')
        call cfwure (nout, mpg%nFc,   gm%fcS,    'fcS')
        call cfwure (nout, mpg%nFc*2, gm%fcHc,   'fcHc')
        call cfwure (nout, mpg%nFc,   gm%fcHt,   'fcHt')
        call cfwure (nout, mpg%nFc*2, gm%fcQgam, 'fcQgam')
        call cfwure (nout, mpg%nFc*2, gm%fcQalf, 'fcQalf')
        call cfwure (nout, mpg%nFc*2, gm%fcQbet, 'fcQbet')
        call cfwure (nout, mpg%nFc,   gm%fcPbs,  'fcPbs')
        call cfwure (nout, mpg%nFc,   gm%fcPbshz,'fcPbshz')

! vertice quantities
        call cfwure (nout, mpg%nVx,   gm%vxX,    'vxX')
        call cfwure (nout, mpg%nVx,   gm%vxY,    'vxY')
        call cfwure (nout, mpg%nVx,   gm%vxFfbz, 'vxFfbz')
        call cfwure (nout, mpg%nVx,   gm%vxFpsi, 'vxFpsi')

        call cfwure (nout, mpg%nFs,   gm%fsConn, 'fsConn')

      end if

      close  (nout)

      return
      end subroutine write_geometry

***********************************************************************

      subroutine prep_geometry (nx, ny, nCv, nCi, nFc, nVx, nFs)

      implicit none
      integer, intent(in) :: nx, ny, nCv, nCi, nFc, nVx, nFs
      integer :: ix, iy, ic, ifc, nxx, ivx, i, j
      real (kind=R8) :: hy1(-1:nx,-1:ny)

      hy1 = hy*qz(:,:,1)

      nxx = (nx+1)*ny

      call alloc_geometry (gm, nCv, nFc, nVx, nFs)

! cell volumes 
      do iy = 0, ny-1
        do ix = 0, nx-1
          ic = (ix+1) + iy*nx
! b field
          gm%cvBb(ic,0:3) = bb(ix,iy,0:3)
! cell center
          gm%cvX(ic) = 0.25_R8*sum(crx(ix,iy,0:3))
          gm%cvY(ic) = 0.25_R8*sum(cry(ix,iy,0:3))
          gm%cvS(ic) = gs(ix,iy,2)
          gm%cvH(ic) = hz(ix,iy)
          gm%cvHx(ic) = hx(ix,iy)
          gm%cvQgam(ic,0) = qz(ix,iy,1)    ! cos
          gm%cvQgam(ic,1) = -qz(ix,iy,0)   ! sin
	  gm%cvVol(ic) = vol(ix,iy)
          gm%cvOnedBsq(ic) = OnedBsq(ix,iy)
          gm%cvBzb(ic) = bzb(ix,iy)
        end do
      end do

      do ix = 0, nx-1
        ic = nCi + ix+1
        gm%cvBb(ic,0:3) = bb(ix,-1,0:3)
        gm%cvX(ic) = 0.25_R8*sum(crx(ix,-1,0:3))
        gm%cvY(ic) = 0.25_R8*sum(cry(ix,-1,0:3))
        gm%cvS(ic) = gs(ix,-1,2)
        gm%cvH(ic) = hz(ix,-1)
        gm%cvHx(ic) = hx(ix,-1)
        gm%cvQgam(ic,0) = qz(ix,-1,1)
        gm%cvQgam(ic,1) = -qz(ix,-1,0)
        gm%cvVol(ic) = vol(ix,-1)
        gm%cvOnedBsq(ic) = OnedBsq(ix,-1)
        gm%cvBzb(ic) = bzb(ix,-1)

        ic = nCi + nx + ix+1
        gm%cvBb(ic,0:3) = bb(ix,ny,0:3)
        gm%cvX(ic) = 0.25_R8*sum(crx(ix,ny,0:3))
        gm%cvY(ic) = 0.25_R8*sum(cry(ix,ny,0:3))
        gm%cvS(ic) = gs(ix,ny,2)
        gm%cvH(ic) = hz(ix,ny)
        gm%cvHx(ic) = hx(ix,ny)
        gm%cvQgam(ic,0) = qz(ix,ny,1)
        gm%cvQgam(ic,1) = -qz(ix,ny,0)
        gm%cvVol(ic) = vol(ix,ny)
        gm%cvOnedBsq(ic) = OnedBsq(ix,ny)
        gm%cvBzb(ic) = bzb(ix,ny)
      end do
 
      do iy = 0, ny-1
        ic = nCi + 2*nx + iy+1
        gm%cvBb(ic,0:3) = bb(-1,iy,0:3)
        gm%cvX(ic) = 0.25_R8*sum(crx(-1,iy,0:3))
        gm%cvY(ic) = 0.25_R8*sum(cry(-1,iy,0:3))
        gm%cvS(ic) = gs(-1,iy,2)
        gm%cvH(ic) = hz(-1,iy)
        gm%cvHx(ic) = hx(-1,iy)
        gm%cvQgam(ic,0) = qz(-1,iy,1)
        gm%cvQgam(ic,1) = -qz(-1,iy,0)
	gm%cvVol(ic) = vol(-1,iy)
        gm%cvOnedBsq(ic) = OnedBsq(-1,iy)
        gm%cvBzb(ic) = bzb(-1,iy)

        ic = nCi + 2*nx + ny + iy+1
        gm%cvBb(ic,0:3) = bb(nx,iy,0:3)
        gm%cvX(ic) = 0.25_R8*sum(crx(nx,iy,0:3))
        gm%cvY(ic) = 0.25_R8*sum(cry(nx,iy,0:3))
        gm%cvS(ic) = gs(nx,iy,2)
        gm%cvH(ic) = hz(nx,iy)
        gm%cvHx(ic) = hx(nx,iy)
        gm%cvQgam(ic,0) = qz(nx,iy,1)
        gm%cvQgam(ic,1) = -qz(nx,iy,0)
	gm%cvVol(ic) = vol(nx,iy)
        gm%cvOnedBsq(ic) = OnedBsq(nx,iy)
        gm%cvBzb(ic) = bzb(nx,iy)
      end do

! x faces
       do iy = 0, ny-1
         do ix = 1, nx-1
           ifc = ix+1 + iy*(nx+1)
           gm%fcBb(ifc,0:3) = wbbl(ix,iy,0:3)
           gm%fcS(ifc) = gs(ix,iy,0)
           gm%fcHc(ifc,1) = 0.5_R8*hx(leftix(ix,iy),leftiy(ix,iy))
           gm%fcHc(ifc,2) = 0.5_R8*hx(ix,iy)
           gm%fcHt(ifc)   = 
     .       0.5_R8*(hy(leftix(ix,iy),leftiy(ix,iy))+hy(ix,iy))             !cwdk: computation could be made exact
           gm%fcVol(ifc,1) = 0.5_R8*vol(leftix(ix,iy),leftiy(ix,iy))
           gm%fcVol(ifc,2) = 0.5_R8*vol(ix,iy)
           gm%fcQgam(ifc,0) = qc(ix,iy,0)
           gm%fcQgam(ifc,1) = -qs(ix,iy,0)
           gm%fcQalf(ifc,0) = gm%fcQgam(ifc,0)
           gm%fcQalf(ifc,1) = gm%fcQgam(ifc,1)
           gm%fcQbet(ifc,0) = 1.0_R8
           gm%fcQbet(ifc,1) = 0.0_R8
           gm%fcPbs(ifc) = pbs(ix,iy,0)
           gm%fcPbshz(ifc) = pbshz(ix,iy,0)
         end do

! left boundary face (to guard cells)
         ix = 0
         ifc = ix+1 + iy*(nx+1)
         gm%fcBb(ifc,0:3) = wbbl(ix,iy,0:3)
         gm%fcS(ifc) = gs(ix,iy,0)
         gm%fcHc(ifc,1) = 0.5_R8*hx(leftix(ix,iy),leftiy(ix,iy))
         gm%fcHc(ifc,2) = 0.5_R8*hx(ix,iy)
         gm%fcHt(ifc)   = hy(leftix(ix,iy),leftiy(ix,iy))                   !cwdk: computation could be made exact
         gm%fcVol(ifc,1) = 0.5_R8*vol(leftix(ix,iy),leftiy(ix,iy))
         gm%fcVol(ifc,2) = 0.5_R8*vol(ix,iy)
         gm%fcQgam(ifc,0) = qc(ix,iy,0)
         gm%fcQgam(ifc,1) = -qs(ix,iy,0)
         gm%fcQalf(ifc,0) = gm%fcQgam(ifc,0)
         gm%fcQalf(ifc,1) = gm%fcQgam(ifc,1)
         gm%fcQbet(ifc,0) = 1.0_R8
         gm%fcQbet(ifc,1) = 0.0_R8
         gm%fcPbs(ifc) = pbs(ix,iy,0)
         gm%fcPbshz(ifc) = pbshz(ix,iy,0)

! right guard cells
         ix = nx
         ifc = ix+1 + iy*(nx+1)
         gm%fcBb(ifc,0:3) = wbbl(ix,iy,0:3)
	 gm%fcS(ifc) = gs(ix,iy,0)
         gm%fcHc(ifc,1) = 0.5_R8*hx(leftix(ix,iy),leftiy(ix,iy))
         gm%fcHc(ifc,2) = 0.5_R8*hx(ix,iy)
         gm%fcHt(ifc)   = hy(ix,iy)                                         !cwdk: computation could be made exact
         gm%fcVol(ifc,1) = 0.5_R8*vol(leftix(ix,iy),leftiy(ix,iy))
         gm%fcVol(ifc,2) = 0.5_R8*vol(ix,iy)
         gm%fcQgam(ifc,0) = qc(ix,iy,0)
         gm%fcQgam(ifc,1) = -qs(ix,iy,0)
         gm%fcQalf(ifc,0) = gm%fcQgam(ifc,0)
         gm%fcQalf(ifc,1) = gm%fcQgam(ifc,1)
         gm%fcQbet(ifc,0) = 1.0_R8
         gm%fcQbet(ifc,1) = 0.0_R8
         gm%fcPbs(ifc) = pbs(ix,iy,0)
         gm%fcPbshz(ifc) = pbshz(ix,iy,0)
       end do

! y faces
!       allocate (qs(-1:nx,-1:ny,0:1))
!       call intfacev(nx, ny, bottomix, bottomiy, vol, qz(-1,-1,0),
!     .               qs(-1,-1,1))
       do ix = 0, nx-1
         do iy = 1, ny-1
           ifc = nxx + iy+1 + ix*(ny+1)
           gm%fcBb(ifc,0:3) = wbbv(ix,iy,0:3)
           gm%fcS(ifc) = gs(ix,iy,1)
           gm%fcHc(ifc,1) = 0.5_R8*hy(bottomix(ix,iy),bottomiy(ix,iy))       !cwdk: change orientation ??
           gm%fcHc(ifc,2) = 0.5_R8*hy(ix,iy)
           gm%fcHt(ifc)   = 
     .       0.5_R8*(hx(bottomix(ix,iy),bottomiy(ix,iy)) + hx(ix,iy))       !cwdk: computation could be made exact
           gm%fcVol(ifc,1) = 0.5_R8*vol(bottomix(ix,iy),bottomiy(ix,iy))
           gm%fcVol(ifc,2) = 0.5_R8*vol(ix,iy)
           gm%fcQgam(ifc,0) = qc(ix,iy,1)
           gm%fcQgam(ifc,1) = -qs(ix,iy,1)
           gm%fcQalf(ifc,0) = 0.0_R8
           gm%fcQalf(ifc,1) = 1.0_R8
           gm%fcQbet(ifc,0) = gm%fcQgam(ifc,1)
           gm%fcQbet(ifc,1) = -gm%fcQgam(ifc,0)
           gm%fcPbs(ifc) = pbs(ix,iy,1)
           gm%fcPbshz(ifc) = pbshz(ix,iy,1)
         end do

! lower guard cells
         iy = 0
         ifc = nxx + ix+1
         gm%fcBb(ifc,0:3) = wbbv(ix,iy,0:3)
         gm%fcS(ifc) = gs(ix,iy,1)
         gm%fcHc(ifc,1) = 0.5_R8*hy(bottomix(ix,iy),bottomiy(ix,iy))
         gm%fcHc(ifc,2) = 0.5_R8*hy(ix,iy)
         gm%fcHt(ifc)   = hx(bottomix(ix,iy),bottomiy(ix,iy))               !cwdk: computation could be made exact
         gm%fcVol(ifc,1) = 0.5_R8*vol(bottomix(ix,iy),bottomiy(ix,iy))
         gm%fcVol(ifc,2) = 0.5_R8*vol(ix,iy)
         gm%fcQgam(ifc,0) = qc(ix,iy,1)
         gm%fcQgam(ifc,1) = -qs(ix,iy,1)
         gm%fcQalf(ifc,0) = 0.0_R8
         gm%fcQalf(ifc,1) = 1.0_R8
         gm%fcQbet(ifc,0) = gm%fcQgam(ifc,1)
         gm%fcQbet(ifc,1) = -gm%fcQgam(ifc,0)
         gm%fcPbs(ifc) = pbs(ix,iy,1)
         gm%fcPbshz(ifc) = pbshz(ix,iy,1)

! upper guard cells
         iy = ny
         ifc = nxx + nx + ix+1
         gm%fcBb(ifc,0:3) = wbbv(ix,iy,0:3)
	 gm%fcS(ifc) = gs(ix,iy,1)
         gm%fcHc(ifc,1) = 0.5_R8*hy(bottomix(ix,iy),bottomiy(ix,iy))
         gm%fcHc(ifc,2) = 0.5_R8*hy(ix,iy)
         gm%fcHt(ifc)   = hx(ix,iy)                                         !cwdk: computation could be made exact
         gm%fcVol(ifc,1) = 0.5_R8*vol(bottomix(ix,iy),bottomiy(ix,iy))
         gm%fcVol(ifc,2) = 0.5_R8*vol(ix,iy)
         gm%fcQgam(ifc,0) = qc(ix,iy,1)
         gm%fcQgam(ifc,1) = -qs(ix,iy,1)
         gm%fcQalf(ifc,0) = 0.0_R8
         gm%fcQalf(ifc,1) = 1.0_R8
         gm%fcQbet(ifc,0) = gm%fcQgam(ifc,1)
         gm%fcQbet(ifc,1) = -gm%fcQgam(ifc,0)
         gm%fcPbs(ifc) = pbs(ix,iy,1)
         gm%fcPbshz(ifc) = pbshz(ix,iy,1)
       end do      

! vertices
       do iy = 0, ny
         do ix = 0, nx
           ivx = ix+1 + iy*(nx+1)
           gm%vxX(ivx) = crx(ix,iy,1)
           gm%vxY(ivx) = cry(ix,iy,1)
           gm%vxFfbz(ivx) = ffbz(ix,iy,1)
           gm%vxFpsi(ivx) = fpsi(ix,iy,1)
         end do
       end do

       deallocate(qs)
 
! connection length for core region

      gm%fsConn = 0.0_R8
      do i = 1, nFs 
	ic = mpg%cvFS(mpg%cvFsP(i,1))
	if (mod(mpg%cvReg(ic),4) == 1) then
	  do j=mpg%fcFsP(i,1), mpg%fcFsP(i,1)+mpg%fcFsp(i,2)-1
	    ifc = mpg%fcFs(j)
	    gm%fsconn(i) = gm%fsconn(i) + (gm%fcHc(ifc,1)+gm%fcHc(ifc,2))*
     .                     gm%fcBb(ifc,3)/abs(gm%fcBb(ifc,0))
	  end do
        end if
      end do      

      end subroutine prep_geometry

**********************************************************************
      
      subroutine prep_plasma (nx, ny, ns, nCv, nCi, nFc, st)
      
      implicit none
      integer, intent(in) :: nx, ny, ns, nCi, nCv, nFc
      type(B2State), intent(inout) :: st
      integer :: ix, iy, ic, ifc

      call createB2Plasma( nCv, ns, st % pl )
      call createB2Coeff( nCv, nFC, ns, st % co )
      call createB2Derivatives( nCv, nFc, ns, st % dv )
      call createB2Source( nCv, ns, st % sr )   
      call createB2SourceWork( nCv, ns, st % srw )   
      call createB2Rates( nCv, ns, st % rt )
      call createB2RatesWork( nCv, ns, st % rtw )
      call createB2PlasmaSnapshot( nCv, nFc, ns, st % psnl )
      call createB2PlasmaSnapshot( nCv, nFc, ns, st % psnc )
      call createB2Diagnostic( nCv, ns, st % diag )
      call createB2Update( nCv, ns, st % update )

      st % ns = ns
      st % ismain = 1

! cell volumes 
      do iy = 0, ny-1
        do ix = 0, nx-1
          ic = (ix+1) + iy*nx
! B2plasma
          st % pl % na(ic,0:ns-1) = na(ix,iy,0:ns-1)
          st % pl % ua(ic,0:ns-1) = ua(ix,iy,0:ns-1)
          st % pl % po(ic) = po(ix,iy)
          st % pl % te(ic) = te(ix,iy)
          st % pl % ti(ic) = ti(ix,iy)

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

        end do
      end do
      
! gueard cells iy = -1 and iy = ny
      do ix = 0, nx-1
        ic = nCi + ix+1
        iy = -1
! B2plasma
        st % pl % na(ic,0:ns-1) = na(ix,iy,0:ns-1)
        st % pl % ua(ic,0:ns-1) = ua(ix,iy,0:ns-1)
        st % pl % po(ic) = po(ix,iy)
        st % pl % te(ic) = te(ix,iy)
        st % pl % ti(ic) = ti(ix,iy)

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

        ic = nCi + nx + ix+1
        iy = ny
! B2plasma
        st % pl % na(ic,0:ns-1) = na(ix,iy,0:ns-1)
        st % pl % ua(ic,0:ns-1) = ua(ix,iy,0:ns-1)
        st % pl % po(ic) = po(ix,iy)
        st % pl % te(ic) = te(ix,iy)
        st % pl % ti(ic) = ti(ix,iy)
        
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
      end do

! guard cells ix = -1 and ix = nx
      do iy = 0, ny-1
        ic = nCi + 2*nx + iy+1
        ix = -1
! B2plasma
        st % pl % na(ic,0:ns-1) = na(ix,iy,0:ns-1)
        st % pl % ua(ic,0:ns-1) = ua(ix,iy,0:ns-1)
        st % pl % po(ic) = po(ix,iy)
        st % pl % te(ic) = te(ix,iy)
        st % pl % ti(ic) = ti(ix,iy)
        
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


        ic = nCi + 2*nx + ny + iy+1
        ix = nx
! B2plasma
        st % pl % na(ic,0:ns-1) = na(ix,iy,0:ns-1)
        st % pl % ua(ic,0:ns-1) = ua(ix,iy,0:ns-1)
        st % pl % po(ic) = po(ix,iy)
        st % pl % te(ic) = te(ix,iy)
        st % pl % ti(ic) = ti(ix,iy)
        
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
      end do



! x faces
      do iy = 0, ny-1
        do ix = 1, nx-1
          ifc = ix+1 + iy*(nx+1)

! B2coeff
          st % co % csig_an(ifc,0:1) = csig_an(ix,iy,0,0:1)
          st % co % calf(ifc,0:1) = calf(ix,iy,0,0:1)

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
          st % co % cthe(ifc,0:ns-1) = cthe(ix,iy,0:ns-1)
          st % co % cthi(ifc,0:ns-1) = cthi(ix,iy,0:ns-1)
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
          st % dv % fna_he(ifc,0:1,0:ns-1) = fna_he(ix,iy,0,0:1,0:ns-1)
          st % dv % fna_32(ifc,0:1,0:ns-1) = fna_32(ix,iy,0,0:1,0:ns-1)
          st % dv % fna_52(ifc,0:1,0:ns-1) = fna_52(ix,iy,0,0:1,0:ns-1)
          st % dv % fna_32nd(ifc,0:1,0:ns-1) = 
     .              fna_32nd(ix,iy,0,0:1,0:ns-1)
          st % dv % fna_52nd(ifc,0:1,0:ns-1) = 
     .              fna_52nd(ix,iy,0,0:1,0:ns-1)
          st % dv % fna_nodrift(ifc,0:1,0:ns-1) = 
     .              fna_nodrift(ix,iy,0,0:1,0:ns-1)

          st % dv % fhe(ifc,0:1) = fhe(ix,iy,0,0:1)
          st % dv % fhe_mdf(ifc,0:1) = fhe_mdf(ix,iy,0,0:1)
          st % dv % fhet(ifc,0:1) = fhet(ix,iy,0,0:1)
          st % dv % fhe_eir(ifc,0:1) = fhe_eir(ix,iy,0,0:1)
          st % dv % fhi(ifc,0:1) = fhi(ix,iy,0,0:1)
          st % dv % fhi_mdf(ifc,0:1) = fhi_mdf(ix,iy,0,0:1)
          st % dv % fhit(ifc,0:1) = fhit(ix,iy,0,0:1)
          st % dv % fhi_eir(ifc,0:1) = fhi_eir(ix,iy,0,0:1)

! derived quantities
          st % dv % uadia(ifc,0:1,0:ns-1) = 
     .              uadia(ix,iy,0,0:1,0:ns-1)
          st % dv % vadia(ifc,0:1,0:ns-1) = 
     .              vadia(ix,iy,0,0:1,0:ns-1)
          st % dv % wadia(ifc,0:1,0:ns-1) = 
     .              wadia(ix,iy,0,0:1,0:ns-1)
          st % dv % vaecrb(ifc,0:1,0:ns-1) = 
     .              vaecrb(ix,iy,0,0:1,0:ns-1)
          st % dv % wedia(ifc,0:1) = 
     .              wedia(ix,iy,0,0:1)
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

!B2PlasmaSnapshot
          st % psnl % fch(ifc,0:1) = fch(ix,iy,0,0:1)
          st % psnl % fna(ifc,0:1,0:ns-1) = fna(ix,iy,0,0:1,0:ns-1)
          st % psnl % fhe(ifc,0:1) = fhe(ix,iy,0,0:1)
          st % psnl % fhi(ifc,0:1) = fhi(ix,iy,0,0:1)
          st % psnc % fch(ifc,0:1) = fch(ix,iy,0,0:1)
          st % psnc % fna(ifc,0:1,0:ns-1) = fna(ix,iy,0,0:1,0:ns-1)
          st % psnc % fhe(ifc,0:1) = fhe(ix,iy,0,0:1)
          st % psnc % fhi(ifc,0:1) = fhi(ix,iy,0,0:1)

        end do
      end do



! y faces

      do ix = 0, nx-1
        do iy = 1, ny-1
          ifc = nxx + iy+1 + ix*(ny+1)

! B2coeff
          st % co % csig_an(ifc,0:1) = csig_an(ix,iy,1,0:1)
          st % co % calf(ifc,0:1) = calf(ix,iy,1,0:1)
        
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
          st % co % cthe(ifc,0:ns-1) = cthe(ix,iy,0:ns-1)
          st % co % cthi(ifc,0:ns-1) = cthi(ix,iy,0:ns-1)
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
          st % dv % fna_nodrift(ifc,0:1,0:ns-1) = 
     .              fna_nodrift(ix,iy,1,0:1,0:ns-1)
        
          st % dv % fhe(ifc,0:1) = fhe(ix,iy,1,0:1)
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
          
! B2PlasmaSnapshot
          st % psnl % fch(ifc,0:1) = fch(ix,iy,1,0:1)
          st % psnl % fna(ifc,0:1,0:ns-1) = fna(ix,iy,1,0:1,0:ns-1)
          st % psnl % fhe(ifc,0:1) = fhe(ix,iy,1,0:1)
          st % psnl % fhi(ifc,0:1) = fhi(ix,iy,1,0:1)
          st % psnc % fch(ifc,0:1) = fch(ix,iy,1,0:1)
          st % psnc % fna(ifc,0:1,0:ns-1) = fna(ix,iy,1,0:1,0:ns-1)
          st % psnc % fhe(ifc,0:1) = fhe(ix,iy,1,0:1)
          st % psnc % fhi(ifc,0:1) = fhi(ix,iy,1,0:1)

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
      
******************************************************************************


      subroutine prep_plasma_ext (nx, ny, ns, nCv, nCi, nFc, st_ext)
      
        implicit none
        integer, intent(in) :: nx, ny, ns, nCi, nCv, nFc
        type(B2StateExt), intent(inout) :: st_ext
        integer :: ix, iy, ic, ifc

        ! For now, only create a dummy external state.
        ! Functionality from b2mod_external to be carried over!
        call createB2StateExt( nCv, nFc, ns, st_ext )

      end subroutine prep_plasma_ext

      end module b2us_prep
