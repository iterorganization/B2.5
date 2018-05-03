      module b2us_prep

      use b2mod_types
      use b2mod_geo
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

      public :: read_mapping, calc_mapping, prep_geometry, prep_plasma

      type(mapping), public, save :: mpg
      type(geometry), public, save :: gm
      type(B2State), public, save :: state

      contains

***********************************************************************

      subroutine read_mapping (fname)

      implicit none
      character(*),intent(in) :: fname
      integer :: iun, io

      open (newunit=iun,file=fname,status='old',iostat=io)

      if (io /= 0) then
        write (*,*) ' file ', fname, ' not readable'
        mpg%nCv = 0
        mpg%nFc = 0
        mpg%nVx = 0
      else
      end if

      close  (iun)

      return
      end subroutine read_mapping 	

***********************************************************************

      subroutine calc_mapping(nxp,nyp,
     .           leftix,leftiy,rightix,rightiy,
     .           topix,topiy,bottomix,bottomiy,
     .           nncut,leftcut,rightcut,topcut,bottomcut)

      implicit none
      integer, intent(in) :: nxp, nyp, nncut
      integer, intent(in) :: leftix(-1:nxp,-1:nyp),
     .                       leftiy(-1:nxp,-1:nyp),
     .                       rightix(-1:nxp,-1:nyp),
     .                       rightiy(-1:nxp,-1:nyp),
     .                       topix(-1:nxp,-1:nyp),
     .                       topiy(-1:nxp,-1:nyp),
     .                       bottomix(-1:nxp,-1:nyp),
     .                       bottomiy(-1:nxp,-1:nyp)
      integer, intent(in) :: leftcut(nncut), rightcut(nncut),
     .                       topcut(nncut), bottomcut(nncut)
      integer :: ix, iy, i, if1, if2, if3, if4, ic, ifc, 
     .           vxFc_end, vxCv_end, nx, ny, icut, ix1, ix2, iy1, iy2,
     .           ivx, ivx1, ivx2, nxx, mxreg, ifs, icv, ireg, istart
      integer, allocatable :: indFc(:), indCv(:), fh(:), ch(:)

!  assume a structured quadrilateral grid

      nx = nxp - 1
      ny = nyp - 1

!  no. of cells
      mpg%nCi = nx*ny 
      mpg%nCg = 2*nx + 2*ny   
      mpg%nCv = mpg%nCi + mpg%nCg  

!  no. of cell faces
      mpg%nFc = ny*(nx+1) + nx*(ny+1)

! no. of vertices
      mpg%nVx = (nx+1) * (ny+1)

! max. no. of vertices per cell
      mpg%nCmxVx = 4

! max. no. of faces per cell
      mpg%nCmxFc = 4

! max no. of flux surfaces
      mpg%nFs = nncut*ny

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
        mpg%cvVx(mpg%cvVxP(ic,1)) = ix 
        mpg%cvVx(mpg%cvVxP(ic,1)+1) = ix+1 
        ic = mpg%nCi + nx + ix
        mpg%cvFc(mpg%cvFcP(ic,1)) = nxx + ny+1 + (ix-1)*(ny+1)
        mpg%cvVx(mpg%cvVxP(ic,1)) = ix + ny*(nx+1)
        mpg%cvVx(mpg%cvVxP(ic,1)+1) = ix+1 + ny*(nx+1)
      end do

      do iy = 1, ny
        ic = mpg%nCi + 2*nx + iy
        mpg%cvFc(mpg%cvFcP(ic,1)) = 1 + (iy-1)*(nx+1)
        mpg%cvVx(mpg%cvVxP(ic,1)) = 1 + (iy-1)*(nx+1) 
        mpg%cvVx(mpg%cvVxP(ic,1)+1) = 1 + iy*(nx+1) 
        ic = mpg%nCi + 2*nx + ny + iy
        mpg%cvFc(mpg%cvFcP(ic,1)) = nx+1 + (iy-1)*(nx+1)
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
           mpg%fcVx(ifc,1) = ix + (iy-1)*(nx+1)
           mpg%fcVx(ifc,2) = ix + iy*(nx+1)
         end do

! left guard cells
         ic = 1 + (iy-1)*nx
         ifc = 1 + (iy-1)*(nx+1)
         mpg%fcCv(ifc,1) = mpg%nCi + 2*nx + iy
         mpg%fcCv(ifc,2) = ic
         mpg%fcVx(ifc,1) = 1 + (iy-1)*(nx+1)
         mpg%fcVx(ifc,2) = 1 + iy*(nx+1)

! right guard cells
         ic = nx+1 + (iy-1)*nx
         ifc = nx+1 + (iy-1)*(nx+1)
         mpg%fcCv(ifc,1) = ic - 1
         mpg%fcCv(ifc,2) = mpg%nCi + 2*nx + ny + iy
         mpg%fcVx(ifc,1) = nx+1 + (iy-1)*(nx+1)
         mpg%fcVx(ifc,2) = nx+1 + iy*(nx+1)
       end do

! y faces
       do ix = 1, nx
         do iy = 2, ny
           ic = ix + (iy-1)*nx
           ifc = nxx + iy + (ix-1)*(ny+1)
           mpg%fcCv(ifc,1) = ic
           mpg%fcCv(ifc,2) = ic - nx
           mpg%fcVx(ifc,1) = ix + (iy-1)*(nx+1)
           mpg%fcVx(ifc,2) = ix+1 + (iy-1)*(nx+1)
         end do

! lower guard cells
         ic = ix 
         ifc = nxx + ix
         mpg%fcCv(ifc,1) = ic
         mpg%fcCv(ifc,2) = mpg%nCi + ix
         mpg%fcVx(ifc,1) = ix 
         mpg%fcVx(ifc,2) = ix+1 

! upper guard cells
         ic = ix + (ny-1)*(nx+1)
         ifc = nxx + nx + ix
         mpg%fcCv(ifc,1) = mpg%nCi + nx + ix
         mpg%fcCv(ifc,2) = ic
         mpg%fcVx(ifc,1) = ix + ny*(nx+1)
         mpg%fcVx(ifc,2) = ix+1 + ny*(nx+1)
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
       do ivx = 1, (nx+1)*(ny+1)-1
         mpg%vxFcP(ivx,2) = count(fh == ivx)
         mpg%vxFcP(ivx+1,1) = mpg%vxFcP(ivx,1) + mpg%vxFcP(ivx+1,2)  
         if (mpg%vxFcP(ivx,2) > 0) then
           mpg%vxFcP(mpg%vxFcP(ivx,1):mpg%vxFcP(ivx+1,1)-1,1) =
     .          pack( indFc, fh == ivx )
         end if
         mpg%vxCvP(ivx,2) = count(ch == ivx)
         mpg%vxCvP(ivx+1,1) = mpg%vxCvP(ivx,1) + mpg%vxCvP(ivx+1,2)  
         if (mpg%vxCvP(ivx,2) > 0) then
           mpg%vxCvP(mpg%vxCvP(ivx,1):mpg%vxCvP(ivx+1,1)-1,1) =
     .          pack( indCv, ch == ivx )
         end if
       end do

       ivx = (nx+1)*(ny+1)
       mpg%vxFcP(ivx,2) = count(fh == ivx)
       vxFc_end = mpg%vxFcP(ivx,1) + mpg%vxFcP(ivx,2)
       if (mpg%vxFcP(ivx,2) > 0) then
         mpg%vxFcP(mpg%vxFcP(ivx,1):vxFc_end,1) =
     .          pack( indFc, fh == ivx )
       end if
       mpg%vxCvP(ivx,2) = count(ch == ivx)
       vxCv_end = mpg%vxCvP(ivx,1) + mpg%vxCvP(ivx,2)
       if (mpg%vxCvP(ivx,2) > 0) then
         mpg%vxCvP(mpg%vxCvP(ivx,1):vxCv_end,1) =
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
	end do
      end do	

      do ix = 0, nx-1
        ic = mpg%nCi + ix+1
        mpg%cvReg(ic) = region(ix,-1,0)

        ic = mpg%nCi + nx + ix+1
        mpg%cvReg(ic) = region(ix,ny,0)
      end do
 
      do iy = 0, ny-1
        ic = mpg%nCi + 2*nx + iy+1
        mpg%cvReg(ic) = region(-1,iy,0)

        ic = mpg%nCi + 2*nx + ny + iy+1
        mpg%cvReg(ic) = region(nx,iy,0)
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

      subroutine prep_geometry (nxp, nyp, nCv, nCi, nFc, nVx, nFs)

      implicit none
      integer, intent(in) :: nxp, nyp, nCv, nCi, nFc, nVx, nFs
      integer :: nx, ny, ix, iy, ic, ifc, nxx, ivx, i, j
      real(r8), allocatable :: qs(:,:,:)

      nx = nxp - 1
      ny = nyp - 1
      nxx = (nx+1)*ny

      call alloc_geometry (gm, nCv, nFc, nVx, nFs)
! cell volumes 
      do iy = 0, ny-1
        do ix = 0, nx-1
          ic = (ix+1) + iy*nx
! b field
          gm%cvBb(ic,0:3) = bb(ix,iy,0:3)
! cell center
          gm%cvX(ic) = 0.25d0*sum(crx(ix,iy,0:3))
          gm%cvY(ic) = 0.25d0*sum(cry(ix,iy,0:3))
          gm%cvS(ic) = gs(ix,iy,2)
          gm%cvH(ic) = hz(ix,iy)
          gm%cvQgam(ic,0) = qz(ix,iy,1)    ! cos
          gm%cvQgam(ic,1) = -qz(ix,iy,0)   ! sin
	  gm%cvVol(ic) = vol(ix,iy)
        end do
      end do

      do ix = 0, nx-1
        ic = nCi + ix+1
        gm%cvBb(ic,0:3) = bb(ix,-1,0:3)
        gm%cvX(ic) = 0.25d0*sum(crx(ix,-1,0:3))
        gm%cvY(ic) = 0.25d0*sum(cry(ix,-1,0:3))
        gm%cvS(ic) = gs(ix,-1,2)
        gm%cvH(ic) = hz(ix,-1)
        gm%cvQgam(ic,0) = qz(ix,-1,1)
        gm%cvQgam(ic,1) = -qz(ix,-1,0)
        gm%cvVol(ic) = vol(ix,-1)

        ic = nCi + nx + ix+1
        gm%cvBb(ic,0:3) = bb(ix,ny,0:3)
        gm%cvX(ic) = 0.25d0*sum(crx(ix,ny,0:3))
        gm%cvY(ic) = 0.25d0*sum(cry(ix,ny,0:3))
        gm%cvS(ic) = gs(ix,ny,2)
        gm%cvH(ic) = hz(ix,ny)
        gm%cvQgam(ic,0) = qz(ix,ny,1)
        gm%cvQgam(ic,1) = -qz(ix,ny,0)
        gm%cvVol(ic) = vol(ix,ny)
      end do
 
      do iy = 0, ny-1
        ic = nCi + 2*nx + iy+1
        gm%cvBb(ic,0:3) = bb(-1,iy,0:3)
        gm%cvX(ic) = 0.25d0*sum(crx(-1,iy,0:3))
        gm%cvY(ic) = 0.25d0*sum(cry(-1,iy,0:3))
        gm%cvS(ic) = gs(-1,iy,2)
        gm%cvH(ic) = hz(-1,iy)
        gm%cvQgam(ic,0) = qz(-1,iy,1)
        gm%cvQgam(ic,1) = -qz(-1,iy,0)
	gm%cvVol(ic) = vol(-1,iy)

        ic = nCi + 2*nx + ny + iy+1
        gm%cvBb(ic,0:3) = bb(nx,iy,0:3)
        gm%cvX(ic) = 0.25d0*sum(crx(nx,iy,0:3))
        gm%cvY(ic) = 0.25d0*sum(cry(nx,iy,0:3))
        gm%cvS(ic) = gs(nx,iy,2)
        gm%cvH(ic) = hz(nx,iy)
        gm%cvQgam(ic,0) = qz(nx,iy,1)
        gm%cvQgam(ic,1) = -qz(nx,iy,0)
	gm%cvVol(ic) = vol(nx,iy)
      end do

! x faces
       do iy = 0, ny-1
         do ix = 1, nx-1
           ifc = ix+1 + iy*(nx+1)
           gm%fcBb(ifc,0:3) = wbbl(ix,iy,0:3)
           gm%fcS(ifc) = gs(ix,iy,0)
           gm%fcHc(ifc,0) = 0.5d0*hx(leftix(ix,iy),leftiy(ix,iy))
           gm%fcHc(ifc,1) = 0.5d0*hx(ix,iy)
           gm%fcHt(ifc,0) = hy(leftix(ix,iy),leftiy(ix,iy))
	   gm%fcHt(ifc,1) = hy(ix,iy)
           gm%fcQgam(ifc,0) = qc(ix,iy,0)
           gm%fcQgam(ifc,1) = sqrt(1.d0-qc(ix,iy,0)**2)
           gm%fcQalf(ifc,0) = gm%fcQgam(ifc,0)
           gm%fcQalf(ifc,1) = gm%fcQgam(ifc,1)
           gm%fcQbet(ifc,0) = 1.d0
           gm%fcQbet(ifc,1) = 0.d0
           gm%fcPbs(ifc) = pbs(ix,iy,0)
           gm%fcPbshz(ifc) = pbshz(ix,iy,0)
         end do

! left boundary face (to guard cells)
         ix = 0
         ifc = ix+1 + iy*(nx+1)
         gm%fcBb(ifc,0:3) = wbbl(ix,iy,0:3)
         gm%fcS(ifc) = gs(ix,iy,0)
         gm%fcHc(ifc,0) = 0.5d0*hx(leftix(ix,iy),leftiy(ix,iy))
         gm%fcHc(ifc,1) = 0.5d0*hx(ix,iy)
         gm%fcHt(ifc,0) = hy(leftix(ix,iy),leftiy(ix,iy))
         gm%fcHt(ifc,1) = hy(ix,iy)
         gm%fcQgam(ifc,0) = qc(ix,iy,0)
         gm%fcQgam(ifc,1) = sqrt(1.d0-qc(ix,iy,0)**2)
         gm%fcQalf(ifc,0) = gm%fcQgam(ifc,0)
         gm%fcQalf(ifc,1) = gm%fcQgam(ifc,1)
         gm%fcQbet(ifc,0) = 1.d0
         gm%fcQbet(ifc,1) = 0.d0
         gm%fcPbs(ifc) = pbs(ix,iy,0)
         gm%fcPbshz(ifc) = pbshz(ix,iy,0)

! right guard cells
         ix = nx
         ifc = ix+1 + iy*(nx+1)
         gm%fcBb(ifc,0:3) = wbbl(ix,iy,0:3)
	 gm%fcS(ifc) = gs(ix,iy,0)
         gm%fcHc(ifc,0) = 0.5d0*hx(leftix(ix,iy),leftiy(ix,iy))
         gm%fcHc(ifc,1) = 0.5d0*hx(ix,iy)
         gm%fcHt(ifc,0) = hy(leftix(ix,iy),leftiy(ix,iy))
         gm%fcHt(ifc,1) = hy(ix,iy)
         gm%fcQgam(ifc,0) = qc(ix,iy,0)
         gm%fcQgam(ifc,1) = sqrt(1.d0-qc(ix,iy,0)**2)
         gm%fcQalf(ifc,0) = gm%fcQgam(ifc,0)
         gm%fcQalf(ifc,1) = gm%fcQgam(ifc,1)
         gm%fcQbet(ifc,0) = 1.d0
         gm%fcQbet(ifc,1) = 0.d0
         gm%fcPbs(ifc) = pbs(ix,iy,0)
         gm%fcPbshz(ifc) = pbshz(ix,iy,0)
       end do

! y faces
       allocate (qs(-1:nx,-1:ny,0:1))
       call intfacev(nx, ny, bottomix, bottomiy, vol, qz(-1,-1,0),
     .               qs(-1,-1,1))
       do ix = 0, nx-1
         do iy = 1, ny-1
           ifc = nxx + iy+1 + ix*(ny+1)
           gm%fcBb(ifc,0:3) = wbbv(ix,iy,0:3)
           gm%fcS(ifc) = gs(ix,iy,1)
           gm%fcHc(ifc,0) = 0.5d0*hy(bottomix(ix,iy),bottomiy(ix,iy))
           gm%fcHc(ifc,1) = 0.5d0*hy(ix,iy)
           gm%fcHt(ifc,0) = hx(bottomix(ix,iy),bottomiy(ix,iy))
           gm%fcHt(ifc,1) = hx(ix,iy)
           gm%fcQgam(ifc,0) = sqrt(1.d0 - qs(ix,iy,1)**2)
           gm%fcQgam(ifc,1) = -qs(ix,iy,1)
           gm%fcQalf(ifc,0) = 1.d0
           gm%fcQalf(ifc,1) = 0.d0
           gm%fcQbet(ifc,0) = 0.d0
           gm%fcQbet(ifc,1) = 1.d0
           gm%fcPbs(ifc) = pbs(ix,iy,1)
           gm%fcPbshz(ifc) = pbshz(ix,iy,1)
         end do

! lower guard cells
         iy = 0
         ifc = nxx + ix+1
         gm%fcBb(ifc,0:3) = wbbv(ix,iy,0:3)
         gm%fcS(ifc) = gs(ix,iy,1)
         gm%fcHc(ifc,0) = 0.5d0*hy(bottomix(ix,iy),bottomiy(ix,iy))
         gm%fcHc(ifc,1) = 0.5d0*hy(ix,iy)
         gm%fcHt(ifc,0) = hx(bottomix(ix,iy),bottomiy(ix,iy))
         gm%fcHt(ifc,1) = hx(ix,iy)
         gm%fcQgam(ifc,0) = sqrt(1.d0 - qs(ix,iy,1)**2)
         gm%fcQgam(ifc,1) = -qs(ix,iy,1)
         gm%fcQalf(ifc,0) = 1.d0
         gm%fcQalf(ifc,1) = 0.d0
         gm%fcQbet(ifc,0) = 0.d0
         gm%fcQbet(ifc,1) = 1.d0
         gm%fcPbs(ifc) = pbs(ix,iy,1)
         gm%fcPbshz(ifc) = pbshz(ix,iy,1)

! upper guard cells
         iy = ny
         ifc = nxx + nx + ix+1
         gm%fcBb(ifc,0:3) = wbbv(ix,iy,0:3)
	 gm%fcS(ifc) = gs(ix,iy,1)
         gm%fcHc(ifc,0) = 0.5d0*hy(bottomix(ix,iy),bottomiy(ix,iy))
         gm%fcHc(ifc,1) = 0.5d0*hy(ix,iy)
         gm%fcHt(ifc,0) = hx(bottomix(ix,iy),bottomiy(ix,iy))
         gm%fcHt(ifc,1) = hx(ix,iy)
         gm%fcQgam(ifc,0) = sqrt(1.d0 - qs(ix,iy,1)**2)
         gm%fcQgam(ifc,1) = -qs(ix,iy,1)
         gm%fcQalf(ifc,0) = 1.d0
         gm%fcQalf(ifc,1) = 0.d0
         gm%fcQbet(ifc,0) = 0.d0
         gm%fcQbet(ifc,1) = 1.d0
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

      gm%fsConn = 0.d0
      do i = 1, nFs 
	ic = mpg%cvFS(mpg%cvFsP(i,1))
	if (mod(mpg%cvReg(ic),4) == 1) then
	  do j=mpg%fcFsP(i,1), mpg%fcFsP(i,1)+mpg%fcFsp(i,2)-1
	    ifc = mpg%fcFs(j)
	    gm%fsconn(i) = gm%fsconn(i) + (gm%fcHc(ifc,0)+gm%fcHc(ifc,1))*
     .                     gm%fcBb(ifc,3)/abs(gm%fcBb(ifc,0))
	  end do
        end if
      end do      

      end subroutine prep_geometry

**********************************************************************
      
      subroutine prep_plasma (nxp, nyp, ns, nCv, nCi, nFc, st)
      
      implicit none
      integer, intent(in) :: nxp, nyp, ns, nCi, nCv, nFc
      type(B2State), intent(inout) :: st
      integer :: nx, ny, ix, iy, ic, ifc

      nx = nxp - 1
      ny = nyp - 1

      call createB2Plasma( nCv, ns, st % pl )
      call createB2Coeff( nCv, nFC, ns, st % co )
      call createB2Derivatives( nCv, nFc, ns, st % dv )
      call createB2Source( nCv, ns, st % sr )   
      call createB2SourceWork( nCv, ns, st % srw )   
      call createB2Rates( nCv, ns, st % rt )
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
          st % co % uadia(ic,0:1,0:1,0:ns-1) = 
     .              uadia(ix,iy,0:1,0:1,0:ns-1)
          st % co % wadia(ic,0:1,0:1,0:ns-1) = 
     .              wadia(ix,iy,0:1,0:1,0:ns-1)
          st % co % vaecrb(ic,0:1,0:1,0:ns-1) = 
     .              vaecrb(ix,iy,0:1,0:1,0:ns-1)

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
          st % dv % pcca(ic,0:1,0:1,0:ns-1) = pcca(ix,iy,0:1,0:1,0:ns-1)
          st % dv % pccm(ic,0:1,0:1) = pccm(ix,iy,0:1,0:1)

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
        st % co % uadia(ic,0:1,0:1,0:ns-1) = 
     &            uadia(ix,iy,0:1,0:1,0:ns-1)
        st % co % wadia(ic,0:1,0:1,0:ns-1) = 
     &            wadia(ix,iy,0:1,0:1,0:ns-1)
        st % co % vaecrb(ic,0:1,0:1,0:ns-1) = 
     &            vaecrb(ix,iy,0:1,0:1,0:ns-1)

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
        st % dv % pcca(ic,0:1,0:1,0:ns-1) = pcca(ix,iy,0:1,0:1,0:ns-1)
        st % dv % pccm(ic,0:1,0:1) = pccm(ix,iy,0:1,0:1)
        
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
        st % co % uadia(ic,0:1,0:1,0:ns-1) = 
     &            uadia(ix,iy,0:1,0:1,0:ns-1)
        st % co % wadia(ic,0:1,0:1,0:ns-1) = 
     &            wadia(ix,iy,0:1,0:1,0:ns-1)
        st % co % vaecrb(ic,0:1,0:1,0:ns-1) = 
     &            vaecrb(ix,iy,0:1,0:1,0:ns-1)
        
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
        st % dv % pcca(ic,0:1,0:1,0:ns-1) = pcca(ix,iy,0:1,0:1,0:ns-1)
        st % dv % pccm(ic,0:1,0:1) = pccm(ix,iy,0:1,0:1)

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
        st % co % uadia(ic,0:1,0:1,0:ns-1) = 
     &            uadia(ix,iy,0:1,0:1,0:ns-1)
        st % co % wadia(ic,0:1,0:1,0:ns-1) = 
     &            wadia(ix,iy,0:1,0:1,0:ns-1)
        st % co % vaecrb(ic,0:1,0:1,0:ns-1) = 
     &            vaecrb(ix,iy,0:1,0:1,0:ns-1)
        
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
        st % dv % pcca(ic,0:1,0:1,0:ns-1) = pcca(ix,iy,0:1,0:1,0:ns-1)
        st % dv % pccm(ic,0:1,0:1) = pccm(ix,iy,0:1,0:1)

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
        st % co % uadia(ic,0:1,0:1,0:ns-1) = 
     &            uadia(ix,iy,0:1,0:1,0:ns-1)
        st % co % wadia(ic,0:1,0:1,0:ns-1) = 
     &            wadia(ix,iy,0:1,0:1,0:ns-1)
        st % co % vaecrb(ic,0:1,0:1,0:ns-1) = 
     &            vaecrb(ix,iy,0:1,0:1,0:ns-1)
        
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
        st % dv % pcca(ic,0:1,0:1,0:ns-1) = pcca(ix,iy,0:1,0:1,0:ns-1)
        st % dv % pccm(ic,0:1,0:1) = pccm(ix,iy,0:1,0:1)

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
          st % co % csig_an(ifc,0,0) = csig_an(ix,iy,0,0)
          st % co % calf(ifc,0,0) = calf(ix,iy,0,0)

          st % co % csig(ifc,0,0) = csig(ix,iy,0,0)
          st % co % csigin(ifc,0,0,0:ns-1,0:ns-1) = 
     .              csigin(ix,iy,0,0,0:ns-1,0:ns-1)

          st % co % chve(ifc,0,0) = chve(ix,iy,0,0)
          st % co % chce(ifc,0,0) = chce(ix,iy,0,0)
          st % co % chvemx(ifc,0) = chvemx(ix,iy)

          st % co % chvi(ifc,0,0) = chvi(ix,iy,0,0)
          st % co % chci(ifc,0,0) = chci(ix,iy,0,0)
          st % co % chvimx(ifc,0) = chvimx(ix,iy)

          st % co % cvla(ifc,0,0,0:ns-1) = cvla(ix,iy,0,0,0:ns-1)
          st % co % cdna(ifc,0,0,0:ns-1) = cdna(ix,iy,0,0,0:ns-1)
          st % co % cdpa(ifc,0,0,0:ns-1) = cdpa(ix,iy,0,0,0:ns-1)

          st % co % cvsa(ifc,0,0,0:ns-1) = cvsa(ix,iy,0,0,0:ns-1)
          st % co % cthe(ifc,0:ns-1) = cthe(ix,iy,0:ns-1)
          st % co % cthi(ifc,0:ns-1) = cthi(ix,iy,0:ns-1)
          st % co % cvsa_cl(ifc,0,0,0:ns-1) = cvsa_cl(ix,iy,0,0,0:ns-1)

          st % co % fllim0fhi(ifc,0,0,0:ns-1) = 
     &              fllim0fhi(ix,iy,0,0,0:ns-1)
          st % co % fllimvisc(ifc,0:ns-1) = fllimvisc(ix,iy,0:ns-1)
      
! fluxes
          st % dv % fch(ifc,0,0) = fch(ix,iy,0,0)
          st % dv % fch_32(ifc,0,0) = fch_32(ix,iy,0,0)
          st % dv % fch_52(ifc,0,0) = fch_52(ix,iy,0,0)
          st % dv % fch_p(ifc,0,0) = fch_p(ix,iy,0,0)

          st % dv % fne(ifc,0,0) = fne(ix,iy,0,0)
          st % dv % fne_32(ifc,0,0) = fne_32(ix,iy,0,0)
          st % dv % fne_52(ifc,0,0) = fne_52(ix,iy,0,0)

          st % dv % fni(ifc,0,0) = fni(ix,iy,0,0)
          st % dv % fni_32(ifc,0,0) = fni_32(ix,iy,0,0)
          st % dv % fni_52(ifc,0,0) = fni_52(ix,iy,0,0)

          st % dv % fna(ifc,0,0,0:ns-1) = fna(ix,iy,0,0,0:ns-1)
          st % dv % fna_32(ifc,0,0,0:ns-1) = fna_32(ix,iy,0,0,0:ns-1)
          st % dv % fna_52(ifc,0,0,0:ns-1) = fna_52(ix,iy,0,0,0:ns-1)
          st % dv % fna_nodrift(ifc,0,0,0:ns-1) = 
     .              fna_nodrift(ix,iy,0,0,0:ns-1)

          st % dv % fhe(ifc,0,0) = fhe(ix,iy,0,0)
          st % dv % fhet(ifc,0,0) = fhet(ix,iy,0,0)
          st % dv % fhi(ifc,0,0) = fhi(ix,iy,0,0)
          st % dv % fhit(ifc,0,0) = fhit(ix,iy,0,0)
      
! numerical coefficients
          st % dv % conc(ifc,0,0) = conc(ix,iy,0,0)
          st % dv % fllime(ifc) = fllime(ix,iy)
          st % dv % floe(ifc,0,0) = floe(ix,iy,0,0)
          st % dv % cone(ifc,0,0,:) = cone(ix,iy,0,0,:)
          st % dv % fllimi(ifc) = fllimi(ix,iy)
          st % dv % floi(ifc,0,0) = floi(ix,iy,0,0)
          st % dv % coni(ifc,0,0,:) = coni(ix,iy,0,0,:)

!B2PlasmaSnapshot
          st % psnl % fch(ifc,0,0) = fch(ix,iy,0,0)
          st % psnl % fna(ifc,0,0,0:ns-1) = fna(ix,iy,0,0,0:ns-1)
          st % psnl % fhe(ifc,0,0) = fhe(ix,iy,0,0)
          st % psnl % fhi(ifc,0,0) = fhi(ix,iy,0,0)
          st % psnc % fch(ifc,0,0) = fch(ix,iy,0,0)
          st % psnc % fna(ifc,0,0,0:ns-1) = fna(ix,iy,0,0,0:ns-1)
          st % psnc % fhe(ifc,0,0) = fhe(ix,iy,0,0)
          st % psnc % fhi(ifc,0,0) = fhi(ix,iy,0,0)

        end do



! left boundary face (to guard cells)
        ix = 0
        ifc = ix+1 + iy*(nx+1)

! B2coeff
        st % co % csig_an(ifc,0,0) = csig_an(ix,iy,0,0)
        st % co % calf(ifc,0,0) = calf(ix,iy,0,0)
        
        st % co % csig(ifc,0,0) = csig(ix,iy,0,0)
        st % co % csigin(ifc,0,0,0:ns-1,0:ns-1) = 
     .            csigin(ix,iy,0,0,0:ns-1,0:ns-1)

        st % co % chve(ifc,0,0) = chve(ix,iy,0,0)
        st % co % chce(ifc,0,0) = chce(ix,iy,0,0)
        st % co % chvemx(ifc,0) = chvemx(ix,iy)
        
        st % co % chvi(ifc,0,0) = chvi(ix,iy,0,0)
        st % co % chci(ifc,0,0) = chci(ix,iy,0,0)
        st % co % chvimx(ifc,0) = chvimx(ix,iy)
        
        st % co % cvla(ifc,0,0,0:ns-1) = cvla(ix,iy,0,0,0:ns-1)
        st % co % cdna(ifc,0,0,0:ns-1) = cdna(ix,iy,0,0,0:ns-1)
        st % co % cdpa(ifc,0,0,0:ns-1) = cdpa(ix,iy,0,0,0:ns-1)
        
        st % co % cvsa(ifc,0,0,0:ns-1) = cvsa(ix,iy,0,0,0:ns-1)
        st % co % cthe(ifc,0:ns-1) = cthe(ix,iy,0:ns-1)
        st % co % cthi(ifc,0:ns-1) = cthi(ix,iy,0:ns-1)
        st % co % cvsa_cl(ifc,0,0,0:ns-1) = cvsa_cl(ix,iy,0,0,0:ns-1)
        
        st % co % fllim0fhi(ifc,0,0,0:ns-1) = 
     &            fllim0fhi(ix,iy,0,0,0:ns-1)
        st % co % fllimvisc(ifc,0:ns-1) = fllimvisc(ix,iy,0:ns-1)
      
! fluxes
        st % dv % fch(ifc,0,0) = fch(ix,iy,0,0)
        st % dv % fch_32(ifc,0,0) = fch_32(ix,iy,0,0)
        st % dv % fch_52(ifc,0,0) = fch_52(ix,iy,0,0)
        st % dv % fch_p(ifc,0,0) = fch_p(ix,iy,0,0)
        
        st % dv % fne(ifc,0,0) = fne(ix,iy,0,0)
        st % dv % fne_32(ifc,0,0) = fne_32(ix,iy,0,0)
        st % dv % fne_52(ifc,0,0) = fne_52(ix,iy,0,0)
        
        st % dv % fni(ifc,0,0) = fni(ix,iy,0,0)
        st % dv % fni_32(ifc,0,0) = fni_32(ix,iy,0,0)
        st % dv % fni_52(ifc,0,0) = fni_52(ix,iy,0,0)
        
        st % dv % fna(ifc,0,0,0:ns-1) = fna(ix,iy,0,0,0:ns-1)
        st % dv % fna_32(ifc,0,0,0:ns-1) = fna_32(ix,iy,0,0,0:ns-1)
        st % dv % fna_52(ifc,0,0,0:ns-1) = fna_52(ix,iy,0,0,0:ns-1)
        st % dv % fna_nodrift(ifc,0,0,0:ns-1) = 
     .            fna_nodrift(ix,iy,0,0,0:ns-1)

        st % dv % fhe(ifc,0,0) = fhe(ix,iy,0,0)
        st % dv % fhet(ifc,0,0) = fhet(ix,iy,0,0)
        st % dv % fhi(ifc,0,0) = fhi(ix,iy,0,0)
        st % dv % fhit(ifc,0,0) = fhit(ix,iy,0,0)
        
! numerical coefficients
        st % dv % conc(ifc,0,0) = conc(ix,iy,0,0)
        st % dv % fllime(ifc) = fllime(ix,iy)
        st % dv % floe(ifc,0,0) = floe(ix,iy,0,0)
        st % dv % cone(ifc,0,0,:) = cone(ix,iy,0,0,:)
        st % dv % fllimi(ifc) = fllimi(ix,iy)
        st % dv % floi(ifc,0,0) = floi(ix,iy,0,0)
        st % dv % coni(ifc,0,0,:) = coni(ix,iy,0,0,:)
        
! B2PlasmaSnapshot
        st % psnl % fch(ifc,0,0) = fch(ix,iy,0,0)
        st % psnl % fna(ifc,0,0,0:ns-1) = fna(ix,iy,0,0,0:ns-1)
        st % psnl % fhe(ifc,0,0) = fhe(ix,iy,0,0)
        st % psnl % fhi(ifc,0,0) = fhi(ix,iy,0,0)
        st % psnc % fch(ifc,0,0) = fch(ix,iy,0,0)
        st % psnc % fna(ifc,0,0,0:ns-1) = fna(ix,iy,0,0,0:ns-1)
        st % psnc % fhe(ifc,0,0) = fhe(ix,iy,0,0)
        st % psnc % fhi(ifc,0,0) = fhi(ix,iy,0,0)
  
      
! right guard cells
        ix = nx
        ifc = ix+1 + iy*(nx+1)

! B2coeff
        st % co % csig_an(ifc,0,0) = csig_an(ix,iy,0,0)
        st % co % calf(ifc,0,0) = calf(ix,iy,0,0)
        
        st % co % csig(ifc,0,0) = csig(ix,iy,0,0)
        st % co % csigin(ifc,0,0,0:ns-1,0:ns-1) = 
     .            csigin(ix,iy,0,0,0:ns-1,0:ns-1)

        st % co % chve(ifc,0,0) = chve(ix,iy,0,0)
        st % co % chce(ifc,0,0) = chce(ix,iy,0,0)
        st % co % chvemx(ifc,0) = chvemx(ix,iy)
        
        st % co % chvi(ifc,0,0) = chvi(ix,iy,0,0)
        st % co % chci(ifc,0,0) = chci(ix,iy,0,0)
        st % co % chvimx(ifc,0) = chvimx(ix,iy)
        
        st % co % cvla(ifc,0,0,0:ns-1) = cvla(ix,iy,0,0,0:ns-1)
        st % co % cdna(ifc,0,0,0:ns-1) = cdna(ix,iy,0,0,0:ns-1)
        st % co % cdpa(ifc,0,0,0:ns-1) = cdpa(ix,iy,0,0,0:ns-1)
        
        st % co % cvsa(ifc,0,0,0:ns-1) = cvsa(ix,iy,0,0,0:ns-1)
        st % co % cthe(ifc,0:ns-1) = cthe(ix,iy,0:ns-1)
        st % co % cthi(ifc,0:ns-1) = cthi(ix,iy,0:ns-1)
        st % co % cvsa_cl(ifc,0,0,0:ns-1) = cvsa_cl(ix,iy,0,0,0:ns-1)
        
        st % co % fllim0fhi(ifc,0,0,0:ns-1) = 
     &            fllim0fhi(ix,iy,0,0,0:ns-1)
        st % co % fllimvisc(ifc,0:ns-1) = fllimvisc(ix,iy,0:ns-1)
        
! fluxes
        st % dv % fch(ifc,0,0) = fch(ix,iy,0,0)
        st % dv % fch_32(ifc,0,0) = fch_32(ix,iy,0,0)
        st % dv % fch_52(ifc,0,0) = fch_52(ix,iy,0,0)
        st % dv % fch_p(ifc,0,0) = fch_p(ix,iy,0,0)
        
        st % dv % fne(ifc,0,0) = fne(ix,iy,0,0)
        st % dv % fne_32(ifc,0,0) = fne_32(ix,iy,0,0)
        st % dv % fne_52(ifc,0,0) = fne_52(ix,iy,0,0)
        
        st % dv % fni(ifc,0,0) = fni(ix,iy,0,0)
        st % dv % fni_32(ifc,0,0) = fni_32(ix,iy,0,0)
        st % dv % fni_52(ifc,0,0) = fni_52(ix,iy,0,0)
        
        st % dv % fna(ifc,0,0,0:ns-1) = fna(ix,iy,0,0,0:ns-1)
        st % dv % fna_32(ifc,0,0,0:ns-1) = fna_32(ix,iy,0,0,0:ns-1)
        st % dv % fna_52(ifc,0,0,0:ns-1) = fna_52(ix,iy,0,0,0:ns-1)
        st % dv % fna_nodrift(ifc,0,0,0:ns-1) = 
     .            fna_nodrift(ix,iy,0,0,0:ns-1)
        
        st % dv % fhe(ifc,0,0) = fhe(ix,iy,0,0)
        st % dv % fhet(ifc,0,0) = fhet(ix,iy,0,0)
        st % dv % fhi(ifc,0,0) = fhi(ix,iy,0,0)
        st % dv % fhit(ifc,0,0) = fhit(ix,iy,0,0)
        
! numerical coefficients
        st % dv % conc(ifc,0,0) = conc(ix,iy,0,0)
        st % dv % fllime(ifc) = fllime(ix,iy)
        st % dv % floe(ifc,0,0) = floe(ix,iy,0,0)
        st % dv % cone(ifc,0,0,:) = cone(ix,iy,0,0,:)
        st % dv % fllimi(ifc) = fllimi(ix,iy)
        st % dv % floi(ifc,0,0) = floi(ix,iy,0,0)
        st % dv % coni(ifc,0,0,:) = coni(ix,iy,0,0,:)
        
!B2PlasmaSnapshot
        st % psnl % fch(ifc,0,0) = fch(ix,iy,0,0)
        st % psnl % fna(ifc,0,0,0:ns-1) = fna(ix,iy,0,0,0:ns-1)
        st % psnl % fhe(ifc,0,0) = fhe(ix,iy,0,0)
        st % psnl % fhi(ifc,0,0) = fhi(ix,iy,0,0)
        st % psnc % fch(ifc,0,0) = fch(ix,iy,0,0)
        st % psnc % fna(ifc,0,0,0:ns-1) = fna(ix,iy,0,0,0:ns-1)
        st % psnc % fhe(ifc,0,0) = fhe(ix,iy,0,0)
        st % psnc % fhi(ifc,0,0) = fhi(ix,iy,0,0)
      end do



! y faces

      do ix = 0, nx-1
        do iy = 1, ny-1
          ifc = nxx + iy+1 + ix*(ny+1)

! B2coeff
          st % co % csig_an(ifc,1,1) = csig_an(ix,iy,1,1)
          st % co % calf(ifc,1,1) = calf(ix,iy,1,1)
        
          st % co % csig(ifc,1,1) = csig(ix,iy,1,1)
          st % co % csigin(ifc,1,1,0:ns-1,0:ns-1) = 
     .              csigin(ix,iy,1,1,0:ns-1,0:ns-1)

          st % co % chve(ifc,1,1) = chve(ix,iy,1,1)
          st % co % chce(ifc,1,1) = chce(ix,iy,1,1)
          st % co % chvemx(ifc,1) = chvemx(ix,iy)
        
          st % co % chvi(ifc,1,1) = chvi(ix,iy,1,1)
          st % co % chci(ifc,1,1) = chci(ix,iy,1,1)
          st % co % chvimx(ifc,1) = chvimx(ix,iy)
        
          st % co % cvla(ifc,1,1,0:ns-1) = cvla(ix,iy,1,1,0:ns-1)
          st % co % cdna(ifc,1,1,0:ns-1) = cdna(ix,iy,1,1,0:ns-1)
          st % co % cdpa(ifc,1,1,0:ns-1) = cdpa(ix,iy,1,1,0:ns-1)
          
          st % co % cvsa(ifc,1,1,0:ns-1) = cvsa(ix,iy,1,1,0:ns-1)
          st % co % cthe(ifc,0:ns-1) = cthe(ix,iy,0:ns-1)
          st % co % cthi(ifc,0:ns-1) = cthi(ix,iy,0:ns-1)
          st % co % cvsa_cl(ifc,1,1,0:ns-1) = cvsa_cl(ix,iy,1,1,0:ns-1)
          
          st % co % fllim0fhi(ifc,1,1,0:ns-1) = 
     &              fllim0fhi(ix,iy,1,1,0:ns-1)
          st % co % fllimvisc(ifc,0:ns-1) = fllimvisc(ix,iy,0:ns-1)
          
! fluxes
          st % dv % fch(ifc,1,1) = fch(ix,iy,1,1)
          st % dv % fch_32(ifc,1,1) = fch_32(ix,iy,1,1)
          st % dv % fch_52(ifc,1,1) = fch_52(ix,iy,1,1)
          st % dv % fch_p(ifc,1,1) = fch_p(ix,iy,1,1)
          
          st % dv % fne(ifc,1,1) = fne(ix,iy,1,1)
          st % dv % fne_32(ifc,1,1) = fne_32(ix,iy,1,1)
          st % dv % fne_52(ifc,1,1) = fne_52(ix,iy,1,1)
          
          st % dv % fni(ifc,1,1) = fni(ix,iy,1,1)
          st % dv % fni_32(ifc,1,1) = fni_32(ix,iy,1,1)
          st % dv % fni_52(ifc,1,1) = fni_52(ix,iy,1,1)
        
          st % dv % fna(ifc,1,1,0:ns-1) = fna(ix,iy,1,1,0:ns-1)
          st % dv % fna_32(ifc,1,1,0:ns-1) = fna_32(ix,iy,1,1,0:ns-1)
          st % dv % fna_52(ifc,1,1,0:ns-1) = fna_52(ix,iy,1,1,0:ns-1)
          st % dv % fna_nodrift(ifc,1,1,0:ns-1) = 
     .              fna_nodrift(ix,iy,1,1,0:ns-1)
        
          st % dv % fhe(ifc,1,1) = fhe(ix,iy,1,1)
          st % dv % fhet(ifc,1,1) = fhet(ix,iy,1,1)
          st % dv % fhi(ifc,1,1) = fhi(ix,iy,1,1)
          st % dv % fhit(ifc,1,1) = fhit(ix,iy,1,1)
        
! numerical coefficients
          st % dv % conc(ifc,1,1) = conc(ix,iy,1,1)
          st % dv % fllime(ifc) = fllime(ix,iy)
          st % dv % floe(ifc,1,1) = floe(ix,iy,1,1)
          st % dv % cone(ifc,1,1,:) = cone(ix,iy,1,1,:)
          st % dv % fllimi(ifc) = fllimi(ix,iy)
          st % dv % floi(ifc,1,1) = floi(ix,iy,1,1)
          st % dv % coni(ifc,1,1,:) = coni(ix,iy,1,1,:)
          
! B2PlasmaSnapshot
          st % psnl % fch(ifc,1,1) = fch(ix,iy,1,1)
          st % psnl % fna(ifc,1,1,0:ns-1) = fna(ix,iy,1,1,0:ns-1)
          st % psnl % fhe(ifc,1,1) = fhe(ix,iy,1,1)
          st % psnl % fhi(ifc,1,1) = fhi(ix,iy,1,1)
          st % psnc % fch(ifc,1,1) = fch(ix,iy,1,1)
          st % psnc % fna(ifc,1,1,0:ns-1) = fna(ix,iy,1,1,0:ns-1)
          st % psnc % fhe(ifc,1,1) = fhe(ix,iy,1,1)
          st % psnc % fhi(ifc,1,1) = fhi(ix,iy,1,1)

        end do

! lower guard cells
        iy = 0
        ifc = nxx + ix+1

! B2coeff
        st % co % csig_an(ifc,1,1) = csig_an(ix,iy,1,1)
        st % co % calf(ifc,1,1) = calf(ix,iy,1,1)
        
        st % co % csig(ifc,1,1) = csig(ix,iy,1,1)
        st % co % csigin(ifc,1,1,0:ns-1,0:ns-1) = 
     .            csigin(ix,iy,1,1,0:ns-1,0:ns-1)

        st % co % chve(ifc,1,1) = chve(ix,iy,1,1)
        st % co % chce(ifc,1,1) = chce(ix,iy,1,1)
        st % co % chvemx(ifc,1) = chvemx(ix,iy)
        
        st % co % chvi(ifc,1,1) = chvi(ix,iy,1,1)
        st % co % chci(ifc,1,1) = chci(ix,iy,1,1)
        st % co % chvimx(ifc,1) = chvimx(ix,iy)
        
        st % co % cvla(ifc,1,1,0:ns-1) = cvla(ix,iy,1,1,0:ns-1)
        st % co % cdna(ifc,1,1,0:ns-1) = cdna(ix,iy,1,1,0:ns-1)
        st % co % cdpa(ifc,1,1,0:ns-1) = cdpa(ix,iy,1,1,0:ns-1)
          
        st % co % cvsa(ifc,1,1,0:ns-1) = cvsa(ix,iy,1,1,0:ns-1)
        st % co % cthe(ifc,0:ns-1) = cthe(ix,iy,0:ns-1)
        st % co % cthi(ifc,0:ns-1) = cthi(ix,iy,0:ns-1)
        st % co % cvsa_cl(ifc,1,1,0:ns-1) = cvsa_cl(ix,iy,1,1,0:ns-1)
          
        st % co % fllim0fhi(ifc,1,1,0:ns-1) = 
     &            fllim0fhi(ix,iy,1,1,0:ns-1)
        st % co % fllimvisc(ifc,0:ns-1) = fllimvisc(ix,iy,0:ns-1)
          
! fluxes
        st % dv % fch(ifc,1,1) = fch(ix,iy,1,1)
        st % dv % fch_32(ifc,1,1) = fch_32(ix,iy,1,1)
        st % dv % fch_52(ifc,1,1) = fch_52(ix,iy,1,1)
        st % dv % fch_p(ifc,1,1) = fch_p(ix,iy,1,1)
          
        st % dv % fne(ifc,1,1) = fne(ix,iy,1,1)
        st % dv % fne_32(ifc,1,1) = fne_32(ix,iy,1,1)
        st % dv % fne_52(ifc,1,1) = fne_52(ix,iy,1,1)
          
        st % dv % fni(ifc,1,1) = fni(ix,iy,1,1)
        st % dv % fni_32(ifc,1,1) = fni_32(ix,iy,1,1)
        st % dv % fni_52(ifc,1,1) = fni_52(ix,iy,1,1)
        
        st % dv % fna(ifc,1,1,0:ns-1) = fna(ix,iy,1,1,0:ns-1)
        st % dv % fna_32(ifc,1,1,0:ns-1) = fna_32(ix,iy,1,1,0:ns-1)
        st % dv % fna_52(ifc,1,1,0:ns-1) = fna_52(ix,iy,1,1,0:ns-1)
        st % dv % fna_nodrift(ifc,1,1,0:ns-1) = 
     .            fna_nodrift(ix,iy,1,1,0:ns-1)
        
        st % dv % fhe(ifc,1,1) = fhe(ix,iy,1,1)
        st % dv % fhet(ifc,1,1) = fhet(ix,iy,1,1)
        st % dv % fhi(ifc,1,1) = fhi(ix,iy,1,1)
        st % dv % fhit(ifc,1,1) = fhit(ix,iy,1,1)
        
! numerical coefficients
        st % dv % conc(ifc,1,1) = conc(ix,iy,1,1)
        st % dv % fllime(ifc) = fllime(ix,iy)
        st % dv % floe(ifc,1,1) = floe(ix,iy,1,1)
        st % dv % cone(ifc,1,1,:) = cone(ix,iy,1,1,:)
        st % dv % fllimi(ifc) = fllimi(ix,iy)
        st % dv % floi(ifc,1,1) = floi(ix,iy,1,1)
        st % dv % coni(ifc,1,1,:) = coni(ix,iy,1,1,:)
          
! B2PlasmaSnapshot
        st % psnl % fch(ifc,1,1) = fch(ix,iy,1,1)
        st % psnl % fna(ifc,1,1,0:ns-1) = fna(ix,iy,1,1,0:ns-1)
        st % psnl % fhe(ifc,1,1) = fhe(ix,iy,1,1)
        st % psnl % fhi(ifc,1,1) = fhi(ix,iy,1,1)
        st % psnc % fch(ifc,1,1) = fch(ix,iy,1,1)
        st % psnc % fna(ifc,1,1,0:ns-1) = fna(ix,iy,1,1,0:ns-1)
        st % psnc % fhe(ifc,1,1) = fhe(ix,iy,1,1)
        st % psnc % fhi(ifc,1,1) = fhi(ix,iy,1,1)

! upper guard cells
        iy = ny
        ifc = nxx + nx + ix+1

! B2coeff
        st % co % csig_an(ifc,1,1) = csig_an(ix,iy,1,1)
        st % co % calf(ifc,1,1) = calf(ix,iy,1,1)
        
        st % co % csig(ifc,1,1) = csig(ix,iy,1,1)
        st % co % csigin(ifc,1,1,0:ns-1,0:ns-1) = 
     .            csigin(ix,iy,1,1,0:ns-1,0:ns-1)

        st % co % chve(ifc,1,1) = chve(ix,iy,1,1)
        st % co % chce(ifc,1,1) = chce(ix,iy,1,1)
        st % co % chvemx(ifc,1) = chvemx(ix,iy)
        
        st % co % chvi(ifc,1,1) = chvi(ix,iy,1,1)
        st % co % chci(ifc,1,1) = chci(ix,iy,1,1)
        st % co % chvimx(ifc,1) = chvimx(ix,iy)
        
        st % co % cvla(ifc,1,1,0:ns-1) = cvla(ix,iy,1,1,0:ns-1)
        st % co % cdna(ifc,1,1,0:ns-1) = cdna(ix,iy,1,1,0:ns-1)
        st % co % cdpa(ifc,1,1,0:ns-1) = cdpa(ix,iy,1,1,0:ns-1)
          
        st % co % cvsa(ifc,1,1,0:ns-1) = cvsa(ix,iy,1,1,0:ns-1)
        st % co % cthe(ifc,0:ns-1) = cthe(ix,iy,0:ns-1)
        st % co % cthi(ifc,0:ns-1) = cthi(ix,iy,0:ns-1)
        st % co % cvsa_cl(ifc,1,1,0:ns-1) = cvsa_cl(ix,iy,1,1,0:ns-1)
          
        st % co % fllim0fhi(ifc,1,1,0:ns-1) = 
     &            fllim0fhi(ix,iy,1,1,0:ns-1)
        st % co % fllimvisc(ifc,0:ns-1) = fllimvisc(ix,iy,0:ns-1)
        
! fluxes
        st % dv % fch(ifc,1,1) = fch(ix,iy,1,1)
        st % dv % fch_32(ifc,1,1) = fch_32(ix,iy,1,1)
        st % dv % fch_52(ifc,1,1) = fch_52(ix,iy,1,1)
        st % dv % fch_p(ifc,1,1) = fch_p(ix,iy,1,1)
          
        st % dv % fne(ifc,1,1) = fne(ix,iy,1,1)
        st % dv % fne_32(ifc,1,1) = fne_32(ix,iy,1,1)
        st % dv % fne_52(ifc,1,1) = fne_52(ix,iy,1,1)
          
        st % dv % fni(ifc,1,1) = fni(ix,iy,1,1)
        st % dv % fni_32(ifc,1,1) = fni_32(ix,iy,1,1)
        st % dv % fni_52(ifc,1,1) = fni_52(ix,iy,1,1)
        
        st % dv % fna(ifc,1,1,0:ns-1) = fna(ix,iy,1,1,0:ns-1)
        st % dv % fna_32(ifc,1,1,0:ns-1) = fna_32(ix,iy,1,1,0:ns-1)
        st % dv % fna_52(ifc,1,1,0:ns-1) = fna_52(ix,iy,1,1,0:ns-1)
        st % dv % fna_nodrift(ifc,1,1,0:ns-1) = 
     .            fna_nodrift(ix,iy,1,1,0:ns-1)
        
        st % dv % fhe(ifc,1,1) = fhe(ix,iy,1,1)
        st % dv % fhet(ifc,1,1) = fhet(ix,iy,1,1)
        st % dv % fhi(ifc,1,1) = fhi(ix,iy,1,1)
        st % dv % fhit(ifc,1,1) = fhit(ix,iy,1,1)
        
! numerical coefficients
        st % dv % conc(ifc,1,1) = conc(ix,iy,1,1)
        st % dv % fllime(ifc) = fllime(ix,iy)
        st % dv % floe(ifc,1,1) = floe(ix,iy,1,1)
        st % dv % cone(ifc,1,1,:) = cone(ix,iy,1,1,:)
        st % dv % fllimi(ifc) = fllimi(ix,iy)
        st % dv % floi(ifc,1,1) = floi(ix,iy,1,1)
        st % dv % coni(ifc,1,1,:) = coni(ix,iy,1,1,:)
          
! B2PlasmaSnapshot
        st % psnl % fch(ifc,1,1) = fch(ix,iy,1,1)
        st % psnl % fna(ifc,1,1,0:ns-1) = fna(ix,iy,1,1,0:ns-1)
        st % psnl % fhe(ifc,1,1) = fhe(ix,iy,1,1)
        st % psnl % fhi(ifc,1,1) = fhi(ix,iy,1,1)
        st % psnc % fch(ifc,1,1) = fch(ix,iy,1,1)
        st % psnc % fna(ifc,1,1,0:ns-1) = fna(ix,iy,1,1,0:ns-1)
        st % psnc % fhe(ifc,1,1) = fhe(ix,iy,1,1)
        st % psnc % fhi(ifc,1,1) = fhi(ix,iy,1,1)
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



      end module b2us_prep
