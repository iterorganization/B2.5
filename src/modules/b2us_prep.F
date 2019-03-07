      module b2us_prep

      use b2mod_types
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

      implicit none

      private

      public :: calc_mapping,
     .          prep_geometry,
     .          prep_plasma, prep_plasma_ext,
     .          to_struct_plasma, read_map_geo, write_map_geo

cwdk  eventually: take this outside of a module! (or: make sure this module is
cwdk  only used at b2mod_driver level)
cwdk      type(mapping), public, save :: mpg
cwdk      type(geometry), public, save :: gm
cwdk      type(B2State), public, save :: state
cwdk      type(B2StateExt), public, save :: state_ext

      contains

***********************************************************************

      subroutine read_map_geo (ninp, mpg, gm)
      use b2mod_types
      use b2mod_version
      use b2mod_b2cmfs

      implicit none
      integer, intent(in) :: ninp
      type (mapping), intent (inout) :: mpg
      type (geometry), intent (inout) :: gm

* ..read the mapping and geometry

      call read_mapping (ninp, mpg)
      call read_geometry (ninp, mpg%nCv, mpg%nFc, mpg%nVx, mpg%nFs, 
     .  mpg%nVmxCv, gm)
      
      end subroutine read_map_geo

***********************************************************************

      subroutine write_map_geo (nout, mpg, gm)
      use b2mod_types
      use b2mod_version
      use b2mod_b2cmfs

      implicit none
      integer :: nout
      type (mapping), intent (in) :: mpg
      type (geometry), intent (in) :: gm

* ..write the geometry

      call write_mapping (nout, mpg)
      call write_geometry (nout, mpg%nCv, mpg%nFc, mpg%nVx, mpg%nFs, 
     .  mpg%nVmxCv, gm)
      
      end subroutine write_map_geo

***********************************************************************

      subroutine calc_mapping(nx, ny, m)
#ifdef WG_DEBUG
      use b2us_debug
#endif
      implicit none
      type (mapping), intent (inout) :: m
      integer, intent(in) :: nx, ny
      integer :: ix, iy, i, if1, if2, if3, if4, ic, ifc, ic1, ic2, 
     .           vxFc_end, vxCv_end, icut, ix1, ix2, iy1, iy2,
     .           ivx, ivx1, ivx2, nfx, mxreg, ifs, icv, ireg, istart,
     .           inv, inv1, inv2, invmx, ifcc, icount, n2, nghosti,
     .           nghostcvi, ncorneri, nghostfcxi, nghostfcyi, 
     .           nghostvxi, ici, icg
      integer, allocatable :: indFc(:), indCv(:), fh(:), ch(:), 
     .           ighost(:,:), icorner(:,:), ighostcv(:,:),
     .           ixghostcv(:), iyghostcv(:), ixcorner(:), iycorner(:),
     .           ighostfcx(:,:), ixghostfcx(:), iyghostfcx(:), 
     .           ighostfcy(:,:), ixghostfcy(:), iyghostfcy(:), 
     .           ighostvx(:,:), ixghostvx(:), iyghostvx(:)
      logical lfound

!  assume a structured quadrilateral grid

!  store information for easy back-transformation to structured format
      m%isClassicalGrid = 1
      m%nx    = nx
      m%ny    = ny
      m%nncut = nncut

      write (*,*) 'Creating mapping ...'

#ifdef WG_DEBUG
      nx_st = nx
      ny_st = ny
#endif

!  no. of cells
! Count internal ghost cells and "corner" points
! Remark in this mapping routine the old B2-numbering is used; do loop over internal cells run from 1->nx, etc. 
! CTB TO BE DONE ALLOCATE THE NEW ARRAYS

! Initialize the mapping matrices
      n2 = (nx+2)*(ny+2)
      allocate (ighost(-1:nx,-1:ny))
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

! Count internal ghost cells
      icount = 0
      do iy = 1, ny
        do ix = 1, nx 
          if((leftix(ix-1,iy-1)==-2).or.(rightix(ix-1,iy-1)==nx+1).or.   
     .       (bottomiy(ix-1,iy-1)==-2).or.(topiy(ix-1,iy-1)==ny+1)) then  
            ighostcv(ix-1,iy-1)=1
            icount=icount+1
            ixghostcv(icount) = ix-1
            iyghostcv(icount) = iy-1
		  endif	
        enddo
      enddo
      nghosti = icount              ! internal ghost cells
      icount = 0
      do ix = 0, nx+1
        ighost(ix-1,-1)=1
        ighost(ix-1,ny)=1
      enddo
      do iy = 1, ny
        ighost(-1,iy-1)=1
        ighost(nx,iy-1)=1
      enddo
      nghostcvi = icount 
      
! Count internal "corner" points, by checking whether it is surrounded by ghost cells
      icount = 0      
      do iy = 1, ny
        do ix = 1, nx 
          if((ighostcv(ix-1,iy-1)==1).and.
     .       (ighostcv(rightix(ix-1,iy-1),rightiy(ix-1,iy-1))==1).and.
     .       (ighostcv(leftix(ix-1,iy-1),leftiy(ix-1,iy-1))==1).and.   
     .       (ighostcv(bottomix(ix-1,iy-1),bottomiy(ix-1,iy-1))==1).and.
     .       (ighostcv(topix(ix-1,iy-1),topiy(ix-1,iy-1))==1)) then  
            icorner(ix-1,iy-1)=1
            icount=icount+1
            ixcorner(icount) = ix-1
            iycorner(icount) = iy-1
          endif	
        enddo
      enddo
      ncorneri = icount

      m%nCi = nx*ny-nghostcvi             ! internal cells
      m%nCg = 2*nx + 2*ny + nghostcvi -ncorneri               ! guard cells
      m%nCv = m%nCi + m%nCg             ! total number of cells

      write (*,*) 'nCv, nCi, nCg', m%nCv, m%nCi, m%nCg

!  no. of cell faces
! Count internal faces in between ghost cells (first x-faces then y-faces)
      icount = 0      
      do iy = 1, ny
        do ix = 1, nx+1
          if((ighostcv(ix-2,iy-1)==1).and.(ighostcv(ix-1,iy-1)==1))then
            ighostfcx(ix-1,iy-1)=1
            icount=icount+1
            ixghostfcx(icount) = ix-1
            iyghostfcx(icount) = iy-1
		  endif	
        enddo
      enddo
      nghostfcxi = icount

      icount = 0      
      do ix = 1, nx
        do iy = 1, ny+1
          if((ighostcv(ix-1,iy-2)==1).and.(ighostcv(ix-1,iy-1)==1)) then
            ighostfcy(ix-1,iy-1)=1
            icount=icount+1
            ixghostfcy(icount) = ix-1
            iyghostfcy(icount) = iy-1
		  endif	
        enddo
      enddo
      nghostfcyi = icount

      m%nFc = ny*(nx+1) + nx*(ny+1)-nghostfcxi-nghostfcyi

! no. of vertices
! Count internal vertices surrrounded by ghost cells
      icount = 0      
      do iy = 1, ny+1
        do ix = 1, nx+1
          if((ighostcv(ix-2,iy-1)==1).and.(ighostcv(ix-1,iy-1)==1).and.
     .       (ighostcv(ix-2,iy-2)==1).and.(ighostcv(ix-1,iy-2)==1)) then           
            ighostvx(ix-1,iy-1)=1
            icount=icount+1
            ixghostvx(icount) = ix-1
            iyghostvx(icount) = iy-1
		  endif	
        enddo
      enddo
      nghostvxi = icount
      m%nVx = (nx+1) * (ny+1) - nghostvxi

! no. of flux surfaces (conservative estimate)
      m%nFs = max(1,(1 + 2*nncut)*ny)

! max. no. of entries in cvVx (vertex list of each cell)
      m%nCmxVx = 4 * m%nCi + 2*m%nCg

! max. no. of entries in cvFc (face list of each cell)
      m%nCmxFc = 4 * m%nCi + m%nCg

! max. no. of entries in vxCv (cell list of each vertex)
! corner vertices have only 3; 
! each X-point vertex has 8; duplicate X-point vertices for slabs
      m%nVmxCv = 4 * m%nVx - 4*nncut + 2*nncut*4 

! max. no. of entries is vxFc (face list of each vertex)
! corner vertices have only 2;
! other vertices on boundaries have only 3; 
! each X-point vertex has 8; duplicate X-point vertices for slabs
      m%nVmxFc = 4 * m%nVx - 4*2*nncut - 2*((nx-1) + (ny-1)) + 2*nncut*4

! max no. of entries in matrix
      m%nCmxNv = 9*m%nCi - 4*nncut + ! internal cells (corners only 8)
     .           6*m%nCg - 8*nncut + ! guard cells (at corners only 5)
     .           2*8*nncut*4         ! correction for cells around X-points

      write (*,*) 'Mapping defined, with'
      write (*,*) 'nCv, nCi, nCg  = ', m%nCv, m%nCi, m%nCg
      write (*,*) 'nFc, nVx       = ', m%nFc, m%nVx
      write (*,*) 'nFs            = ', m%nFs
      write (*,*) 'nCmxFc, nCmxVx = ', m%nCmxFc, m%nCmxVx
      write (*,*) 'nVmxFc, nVmxCv = ', m%nVmxFc, m%nVmxCv
      write (*,*) 'nCmxNv         = ', m%nCmxNv

      call alloc_mapping (m)

! Compose the mapping matrix from structured to unstructured index

! for Cv's
      m%imapCv= 0 ! 0 means not in use in unstructured
      ici =0
      icg = m%nCi
      do iy = 1, ny
        do ix = 1, nx
          if(ighostcv(ix-1,iy-1)==0) then
            ici = ici + 1          
            m%imapCv(ix-1,iy-1) = ici
          elseif(icorner(ix-1,iy-1)==0) then
            icg = icg+1
            m%imapCv(ix-1,iy-1) = icg
          endif
        enddo
      enddo

      iy = 0
      do ix = 1, nx
        if(icorner(ix-1,iy-1)==0) then
          icg = icg+1
          m%imapCv(ix-1,iy-1) = icg
        endif
      enddo
      iy=ny+1
      do ix = 1, nx
        if(icorner(ix-1,iy-1)==0) then
          icg = icg+1
          m%imapCv(ix-1,iy-1) = icg
        endif
      enddo
      ix = 0
      do iy = 1, ny
        if(icorner(ix-1,iy-1)==0) then
          icg = icg+1
          m%imapCv(ix-1,iy-1) = icg
        endif
      enddo     
      ix = nx+1
      do iy = 1, ny
        if(icorner(ix-1,iy-1)==0) then
          icg = icg+1
          m%imapCv(ix-1,iy-1) = icg
        endif
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
      ici =0
      do iy = 1, ny
        do ix = 1, nx+1
          if(ighostfcx(ix-1,iy-1)==0) then
            ici = ici + 1          
            m%imapFcx(ix-1,iy-1) = ici
           endif
        enddo
      enddo

      m%imapFcy= 0 ! 0 means not in use in unstructured
      do ix = 1, nx
        do iy = 1, ny+1
          if(ighostfcy(ix-1,iy-1)==0) then
            ici = ici + 1          
            m%imapFcy(ix-1,iy-1) = ici
           endif
        enddo
      enddo
   
      if (ici /= m%nFc) then
        call xerrab 
     .  ('Inconsistency between counting and mapping in nFc!')
      endif
 
! for Vx
      m%imapVx= 0 ! 0 means not in use in unstructured
      ici =0
      do iy = 1, ny+1
        do ix = 1, nx+1
          if(ighostfcx(ix-1,iy-1)==0) then
            ici = ici + 1          
            m%imapVx(ix-1,iy-1) = ici
           endif
        enddo
      enddo

      if (ici /= m%nVx) then
        call xerrab 
     .    ('Inconsistency between counting and mapping in nVx!')
      endif

! no. of faces to be considered for this cell
      m%cvFcP(1:m%nCi,2) = 4
      m%cvFcP(m%nCi+1:m%nCv,2) = 1

! no. of vertices to be considered for this cell
      m%cvVxP(1:m%nCi,2) = 4
      m%cvVxP(m%nCi+1:m%nCv,2) = 2
      
! starting index in the list of cell faces
! starting index for every cell in the cvVx list and its number of vertices
      m%cvFcP(1,1) = 1
      m%cvVxP(1,1) = 1
      do i=2, m%nCv
        m%cvFcP(i,1) = m%cvFcP(i-1,1) + m%cvFcP(i-1,2)
        m%cvVxP(i,1) = m%cvVxP(i-1,1) + m%cvVxP(i-1,2)
      end do

! list of cell faces and vertices for inner cells
! interpolation coefficients for interpolation to cell centers
      nfx = (nx+1)*ny
      do iy = 1, ny
        do ix = 1, nx
          if(m%imapCv(ix-1,iy-1) /= 0) then
            ic = m%imapCv(ix-1,iy-1)!! seems not te be needed in this loop??
!          ic = ix + (iy-1)*nx
            if1 = m%cvFcP(ic,1)
            if2 = if1 + 1
            if3 = if2 + 1
            if4 = if3 + 1
            m%cvFc(if1) = m%imapFcx(ix-1,iy-1)    
!          m%cvFc(if1) = ix + (iy-1)*(nx+1)
!!          m%cvFc(if1) = leftix(ix-1,iy-1)+2 + (iy-1)*(nx+1)
!!          m%cvFc(if2) = ix+1 + (iy-1)*(nx+1)
            m%cvFc(if2) = m%imapFcx(rightix(ix-1,iy-1),
     .                              rightiy(ix-1,iy-1))
!          m%cvFc(if2) = rightix(ix-1,iy-1)+1 + rightiy(ix-1,iy-1)*(nx+1)
            m%cvFc(if3) = m%imapFcy(ix-1,iy-1)   
!          m%cvFc(if3) = nfx + iy + (ix-1)*(ny+1)
!!          m%cvFc(if4) = nfx + iy+1 + (ix-1)*(ny+1)
            m%cvFc(if4) = m%imapFcy(topix(ix-1,iy-1),topiy(ix-1,iy-1))
!          m%cvFc(if4) = nfx + topiy(ix-1,iy-1)+1+topix(ix-1,iy-1)*(ny+1)
            m%intcellP(if1) = 1.0_R8
            m%intcellP(if2) = 1.0_R8
            m%intcellP(if3) = 0.0_R8
            m%intcellP(if4) = 0.0_R8
            m%intcellR(if1) = 0.0_R8
            m%intcellR(if2) = 0.0_R8
            m%intcellR(if3) = 1.0_R8
            m%intcellR(if4) = 1.0_R8
            m%cvVx(if1) = m%imapVx(ix-1,iy-1)
!          m%cvVx(if1) = ix + (iy-1)*(nx+1)
!!          m%cvVx(if1) = leftix(ix-1,iy-1)+2 + (iy-1)*(nx+1)
!!          m%cvVx(if2) = ix+1 + (iy-1)*(nx+1)
            m%cvVx(if2) = m%imapVx(rightix(ix-1,iy-1),
     .                             rightiy(ix-1,iy-1))
!          m%cvVx(if2) = rightix(ix-1,iy-1)+1 + rightiy(ix-1,iy-1)*(nx+1)
            m%cvVx(if3) = m%imapVx(topix(ix-1,iy-1),topiy(ix-1,iy-1))
!          m%cvVx(if3) = topix(ix-1,iy-1)+1 + topiy(ix-1,iy-1)*(nx+1)
!!          m%cvVx(if3) = leftix(ix-1,iy-1)+2 + iy*(nx+1)
!!          m%cvVx(if4) = ix+1 + iy*(nx+1)
            m%cvVx(if4) = m%imapVx(rightix(topix(ix-1,iy-1),
     .                             topiy(ix-1,iy-1)), 
     .                 rightiy(topix(ix-1,iy-1),topiy(ix-1,iy-1)))
!          m%cvVx(if4) = rightix(topix(ix-1,iy-1),topiy(ix-1,iy-1))+1 + 
!     .                 rightiy(topix(ix-1,iy-1),topiy(ix-1,iy-1))*(nx+1)
          endif
        end do
      end do



! list of cell faces and vertices for guard cells
      do ix = 1, nx
        if(m%imapCv(ix-1,-1) /= 0 ) then
          ic = m%imapCv(ix-1,-1)
          m%cvFc(m%cvFcP(ic,1)) = m%imapFcy(topix(ix-1,-1),
     .                                      topiy(ix-1,-1))
          m%intcellP(m%cvFcP(ic,1)) = 1.0_R8
          m%intcellR(m%cvFcP(ic,1)) = 1.0_R8
          m%cvVx(m%cvVxP(ic,1)) = m%imapVx(topix(ix-1,-1),
     .                                     topiy(ix-1,-1))
          m%cvVx(m%cvVxP(ic,1)+1)= m%imapVx(
     .               rightix(topix(ix-1,-1),topiy(ix-1,-1)), 
     .               rightiy(topix(ix-1,-1),topiy(ix-1,-1))) 
!        ic = m%nCi + ix
!        m%cvFc(m%cvFcP(ic,1)) = nfx + 
!     .               topiy(ix-1,-1)+1 + topix(ix-1,-1)*(ny+1)
!        m%intcellP(m%cvFcP(ic,1)) = 1.0_R8
!        m%intcellR(m%cvFcP(ic,1)) = 1.0_R8
!        m%cvVx(m%cvVxP(ic,1)) = 
!     .               topix(ix-1,-1)+1 + topiy(ix-1,-1)*(nx+1)
!        m%cvVx(m%cvVxP(ic,1)+1) = 
!     .               rightix(topix(ix-1,-1),topiy(ix-1,-1))+1 + 
!     .               rightiy(topix(ix-1,-1),topiy(ix-1,-1))*(nx+1) 
        endif
        if(m%imapCv(ix-1,ny)/=0) then
          ic = m%imapCv(ix-1,ny)
          m%cvFc(m%cvFcP(ic,1)) = m%imapFcy(ix-1,ny)
          m%intcellP(m%cvFcP(ic,1)) = 1.0_R8
          m%intcellR(m%cvFcP(ic,1)) = 1.0_R8
          m%cvVx(m%cvVxP(ic,1)) = m%imapVx(ix-1,ny)
          m%cvVx(m%cvVxP(ic,1)+1)= m%imapVx(rightix(ix-1,ny),
     .                                     rightiy(ix-1,ny))         
!        ic = m%nCi + nx + ix
!     .               rightiy(topix(ix-1,-1),topiy(ix-1,-1))*(nx+1) 
!        m%cvFc(m%cvFcP(ic,1)) = nfx + ny+1 + (ix-1)*(ny+1)
!        m%intcellP(m%cvFcP(ib2us_prep.Fc,1)) = 1.0_R8
!        m%intcellR(m%cvFcP(ic,1)) = 1.0_R8
!        m%cvVx(m%cvVxP(ic,1)) = ix + ny*(nx+1)
!        m%cvVx(m%cvVxP(ic,1)+1) = 
!     .               rightix(ix-1,ny)+1 + rightiy(ix-1,ny)*(nx+1)
         endif      
      end do

      do iy = 1, ny
        if(m%imapCv(-1,iy-1)/=0) then
          ic = m%imapCv(-1,iy-1)
          m%cvFc(m%cvFcP(ic,1)) = m%imapFcx(rightix(-1,iy-1),
     .                                      rightiy(-1,iy-1))
          m%intcellP(m%cvFcP(ic,1)) = 1.0_R8
          m%intcellR(m%cvFcP(ic,1)) = 1.0_R8
          m%cvVx(m%cvVxP(ic,1)) = m%imapVx(rightix(-1,iy-1),
     .                                     rightiy(-1,iy-1))
          m%cvVx(m%cvVxP(ic,1)+1)= m%imapVx(rightix(topix(-1,iy-1),
     .                                              topiy(-1,iy-1)),
     .                                      rightiy(topix(-1,iy-1),
     .                                              topiy(-1,iy-1)))
        endif 
!        ic = m%nCi + 2*nx + iy
!        m%cvFc(m%cvFcP(ic,1)) = 
!     .               rightix(-1,iy-1)+1 + rightiy(-1,iy-1)*(nx+1)
!        m%intcellP(m%cvFcP(ic,1)) = 1.0_R8
!        m%intcellR(m%cvFcP(ic,1)) = 1.0_R8
!        m%cvVx(m%cvVxP(ic,1)) = 
!     .               rightix(-1,iy-1)+1 + rightiy(-1,iy-1)*(nx+1)
!        m%cvVx(m%cvVxP(ic,1)+1) = 
!     .               rightix(topix(-1,iy-1),topiy(-1,iy-1))+1 + 
!     .               rightiy(topix(-1,iy-1),topiy(-1,iy-1))*(nx+1)
        if(m%imapCv(nx,iy-1)/=0) then        
          ic = m%imapCv(nx,iy-1)
          m%cvFc(m%cvFcP(ic,1)) = m%imapFcx(nx,iy-1)
          m%intcellP(m%cvFcP(ic,1)) = 1.0_R8
          m%intcellR(m%cvFcP(ic,1)) = 1.0_R8
          m%cvVx(m%cvVxP(ic,1)) = m%imapVx(nx,iy-1)
          m%cvVx(m%cvVxP(ic,1)+1)= m%imapVx(topix(nx,iy-1),
     .                                      topiy(nx,iy-1))
        endif
!        ic = m%nCi + 2*nx + ny + iy
!        m%cvFc(m%cvFcP(ic,1)) = nx+1 + (iy-1)*(nx+1)
!        m%intcellP(m%cvFcP(ic,1)) = 1.0_R8
!        m%intcellR(m%cvFcP(ic,1)) = 1.0_R8
!        m%cvVx(m%cvVxP(ic,1)) = nx+1 + (iy-1)*(nx+1) 
!        m%cvVx(m%cvVxP(ic,1)+1) = 
!     .               topix(nx,iy-1)+1 + topiy(nx,iy-1)*(nx+1)
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
             m%fcVx(ifc,2) = m%imapCv(topix(ix-1,iy-1),topiy(ix-1,iy-1))
           endif  
         end do
       end do

!       do iy = 1, ny
!         do ix = 2, nx
!           ic = ix + (iy-1)*nx
!           ifc = ix + (iy-1)*(nx+1)
!           m%fcCv(ifc,1) = leftix(ix-1,iy-1)+1 + (iy-1)*nx
!           m%fcCv(ifc,2) = ic
!           m%fcVx(ifc,1) = ix + (iy-1)*(nx+1)
!           m%fcVx(ifc,2) = ix + iy*(nx+1)
!         end do

!! left guard cells
!         ic = 1 + (iy-1)*nx
!         ifc = 1 + (iy-1)*(nx+1)
!         m%fcCv(ifc,1) = m%nCi + 2*nx + iy
!         m%fcCv(ifc,2) = ic
!         m%fcVx(ifc,1) = 1 + (iy-1)*(nx+1)
!         m%fcVx(ifc,2) = 1 + iy*(nx+1)

!! right guard cells
!         ic = nx+1 + (iy-1)*nx
!         ifc = nx+1 + (iy-1)*(nx+1)
!         m%fcCv(ifc,1) = ic - 1
!         m%fcCv(ifc,2) = m%nCi + 2*nx + ny + iy
!         m%fcVx(ifc,1) = nx+1 + (iy-1)*(nx+1)
!         m%fcVx(ifc,2) = nx+1 + iy*(nx+1)
!       end do

! y faces

       do ix = 1, nx 
         do iy = 1, ny+1
           if(m%imapFcy(ix-1,iy-1)/=0) then
             ifc =m%imapFcy(ix-1,iy-1)
             m%fcCv(ifc,2) = m%imapCv(ix-1,iy-1)
             m%fcCv(ifc,1) = m%imapCv(bottomix(ix-1,iy-1),
     .                                bottomiy(ix-1,iy-1))
             m%fcVx(ifc,2) = m%imapVx(ix-1,iy-1)
             m%fcVx(ifc,1) = m%imapCv(rightix(ix-1,iy-1),
     .                                rightiy(ix-1,iy-1))
           endif
         end do
       end do

!       do ix = 1, nx
!         do iy = 2, ny
!           ic = ix + (iy-1)*nx
!           ifc = nfx + iy + (ix-1)*(ny+1)
!           m%fcCv(ifc,2) = ic
!           m%fcCv(ifc,1) = ic - nx
!           m%fcVx(ifc,2) = ix + (iy-1)*(nx+1)
!           m%fcVx(ifc,1) = rightix(ix-1,iy-1)+1 + (iy-1)*(nx+1)
!         end do

!! lower guard cells
!         ic = ix
!         ifc = nfx + 1 + (ix-1)*(ny+1)
!         m%fcCv(ifc,2) = ic
!         m%fcCv(ifc,1) = m%nCi + ix
!         m%fcVx(ifc,2) = ix 
!         m%fcVx(ifc,1) = rightix(ix-1,0)+1 

! upper guard cells
!         ic = ix + (ny-1)*nx
!         ifc = nfx + ix*(ny+1)
!         m%fcCv(ifc,2) = m%nCi + nx + ix
!         m%fcCv(ifc,1) = ic
!         m%fcVx(ifc,2) = ix + ny*(nx+1)
!         m%fcVx(ifc,1) = rightix(ix-1,ny)+1 + ny*(nx+1)
!       end do

! merging x-point vertex !!! TB not changed... seems not to be necessary
       do icut = 1, nncut
         ix1 = leftcut(icut)
         iy1 = topcut(icut)
         ix2 = rightcut(icut)
         iy2 = topcut(icut)
         ivx1 = m%imapVx(ix1,iy1)
         ivx2 = m%imapVx(ix2,iy2)

         where (m%imapVx == ivx2)
           m%imapVx = ivx1
         end where

         where (m%cvVx == ivx2)
           m%cvVx = ivx1
         end where

         where (m%fcVx == ivx2)
           m%fcVx = ivx1
         end where

       end do

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
       do ivx = 1, (nx+1)*(ny+1)-1
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

       ivx = (nx+1)*(ny+1)
       m%vxFcP(ivx,2) = count(fh == ivx)
       vxFc_end = m%vxFcP(ivx,1) + m%vxFcP(ivx,2)
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
      do iy = 0, ny+1
        do ix = 0, nx+1
          if(m%imapCv(ix-1,iy-1)/=0) then          
            ic = m%imapCv(ix-1,iy-1)
!          ic = (ix+1) + iy*nx
            m%cvReg(ic) = region(ix-1,iy-1,0)
            m%cvOnClosedSurface(ic) = on_closed_surface(ix-1,iy-1)
          endif
        end do
      end do	

!      do ix = 1, nx
!        ic = m%nCi + ix+1
!        m%cvReg(ic) = region(ix,-1,0)
!        m%cvOnClosedSurface(ic) = on_closed_surface(ix,-1)
!        ic = m%nCi + nx + ix+1
!        m%cvReg(ic) = region(ix,ny,0)
!        m%cvOnClosedSurface(ic) = on_closed_surface(ix,ny)
!      end do
 
!      do iy = 0, ny-1
!        ic = m%nCi + 2*nx + iy+1
!        m%cvReg(ic) = region(-1,iy,0)
!        m%cvOnClosedSurface(ic) = on_closed_surface(-1,iy)
!
!        ic = m%nCi + 2*nx + ny + iy+1
!        m%cvReg(ic) = region(nx,iy,0)
!        m%cvOnClosedSurface(ic) = on_closed_surface(nx,iy)
!      end do

! regions (faces)
! x faces
       do iy = 1, ny
         do ix = 1, nx+1
           if(m%imapFcx(ix-1,iy-1)/=0) then
             ifc = m%imapFcx(ix-1,iy-1)
             m%fcReg(ifc) = region(ix-1,iy-1,1)
           endif
         end do

!       do iy = 0, ny-1
!         do ix = 1, nx-1
!           ifc = ix+1 + iy*(nx+1)
!           m%fcReg(ifc) = region(ix,iy,1)
!         end do

!! left boundary face (to guard cells)
!         ix = 0
!         ifc = ix+1 + iy*(nx+1)
!         m%fcReg(ifc) = region(ix,iy,1)

! right guard cells
!         ix = nx
!         ifc = ix+1 + iy*(nx+1)
!         m%fcReg(ifc) = region(ix,iy,1)
        end do

! y faces
       do ix = 1, nx
         do iy = 1, ny+1
           if(m%imapFcy(ix-1,iy-1)/=0) then
             ifc = m%imapFcy(ix-1,iy-1)
             m%fcReg(ifc) = region(ix-1,iy-1,2)
           end if
         end do

!       do ix = 0, nx-1
!         do iy = 1, ny-1
!           ifc = nfx + iy+1 + ix*(ny+1)
!           m%fcReg(ifc) = region(ix,iy,2)
!         end do

!! lower guard cells
!         iy = 0
!         ifc = nfx + 1 + ix*(ny+1)
!         m%fcReg(ifc) = region(ix,iy,2)

!! upper guard cells
!         iy = ny
!         ifc = nfx + (ix+1)*(ny+1)
!         m%fcReg(ifc) = region(ix,iy,2)
      end do

      mxReg = maxval(m%cvReg(1:m%nCv))

! identify the various flux surfaces
      m%fsCv = 0 ! init to 0
      m%fsFc = 0 ! init to 0

      ifs = 0
      icv = 0
      ifcc = 0
      do ireg= 1, mxReg
    	do iy = 1, ny
	      istart = -1
	      do ix = 1, nx
            if(m%imapCv(ix-1,iy-1)/=0) then
              ic = m%imapCv(ix-1,iy-1)
              ifc = m%imapFcx(ix-1,iy-1)
!	do iy = 0, ny-1
!	  istart = -1
!	  do ix = 0, nx-1
!            ic = (ix+1) + iy*nx
!            ifc = ix+1 + iy*(nx+1)
              if (m%cvReg(ic) == ireg) then
	            if (istart == -1) then
                  ifs = ifs + 1
                  icv = icv + 1
                  ifcc = ifcc + 1
                  m%fsCv(icv) = ic
	              m%fsCvP(ifs,1) = icv
	              m%fsFc(ifcc) = ifc
                  m%fsFcP(ifs,1) = ifcc
     	          istart = 1
                else
     	          icv = icv + 1
                  ifcc = ifcc + 1
        		  m%fsCv(icv) = ic
     	          m%fsFc(ifcc) = ifc
                  m%fsCvP(ifs,2) = icv - m%fsCvP(ifs,1) + 1
                  m%fsFcP(ifs,2) = ifcc - m%fsFcP(ifs,1) + 1
                  if (ix == (nx-1)) then
                    ifcc = ifcc + 1
                    m%fsFc(ifcc) = ifc + 1
                    m%fsFcP(ifs,2) = ifcc - m%fsFcP(ifs,1) + 1
                  end if
                end if
              end if
            end if 
          end do
        end do	
      end do
      m%nFs = ifs


#ifdef WG_TODO
! moved to init_mapping!

! setting up connectivity information for the matrix stencil
! consider writing separate routine/module for matrix setup; is generic

      allocate (indCv(4*4+4+1)) ! max 4 vertices per cell, max 8 cells for each vertex (very conservative)
      inv = 1
      do icv = 1, m%nCv
        m%cvNvP(icv,1) = inv
        ! init: each cell appears in its own stencil
        indCv = icv
        invmx = 1
        ! list all cells surrounding cell icv, by checking the
        ! cells connecting to all vertices of cell icv
        do ivx = 1, m%cvVxP(icv,2)
          ivx1 = m%cvVx(m%cvVxP(icv,1)+ivx-1)     ! the number in the vertex list
          indCv(invmx+1:invmx+m%vxCvP(ivx1,2)) =  ! all cells connecting to this vertex
     .      m%vxCv(m%vxCvP(ivx1,1):m%vxCvP(ivx1,1)+m%vxCvP(ivx1,2)-1)
          invmx = invmx + m%vxCvP(ivx1,2)
        enddo
        ! remove duplicates
        m%cvNv(m%cvNvP(icv,1)) = indCv(1)
        i = 1
        outer: do inv1 = 2, invmx
          do inv2 = 1, i 
            if (m%cvNv(m%cvNvP(icv,1)+inv2-1).eq.indCv(inv1)) then
              ! duplicate element
              cycle outer
            endif
          enddo
          ! new element found
          i = i + 1
          m%cvNv(m%cvNvP(icv,1)+i-1) = indCv(inv1)
        end do outer
        ! TODO: sort by cell number, keeping first cell fixed (main diagonal)
        m%cvNvP(icv,2) = i
        inv = inv + i
      enddo

      deallocate(indCv)
#endif
      return
      end subroutine calc_mapping

***********************************************************************

      subroutine prep_geometry(nx, ny, nCv, nCi, nFc, nVx, nFs, nVmxCv,
     .  gm, mpg)

      implicit none
      integer, intent(in) :: nx, ny, nCv, nCi, nFc, nVx, nFs, nVmxCv
      integer :: ix, iy, ic, ifc, nfx, ivx, i, j
      real (kind=R8) :: hy1(-1:nx,-1:ny), du, dux, duy
      type (mapping) , intent (in) :: mpg
      type (geometry), intent (out) :: gm

      hy1 = hy*qz(:,:,1)

      nfx = (nx+1)*ny

      call alloc_geometry (gm, nCv, nFc, nVx, nFs, nVmxCv)

! cell volumes 
!      do iy = 0, ny-1
!        do ix = 0, nx-1
      do iy = -1, ny
        do ix = -1, nx
           if(mpg%imapCv(ix,iy)/=0) then           
             ic = mpg%imapCv(ix,iy)
!          ic = (ix+1) + iy*nx
! b field
            gm%cvBb(ic,0:3) = bb(ix,iy,0:3)
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
            gm%cvBzb(ic) = bzb(ix,iy)
          end if 
        end do
      end do

!      do ix = 0, nx-1
!        ic = nCi + ix+1
!        gm%cvBb(ic,0:3) = bb(ix,-1,0:3)
!        gm%cvX(ic) = 0.25_R8*sum(crx(ix,-1,0:3))
!        gm%cvY(ic) = 0.25_R8*sum(cry(ix,-1,0:3))
!        gm%cvSz(ic) = gs(ix,-1,2)
!        gm%cvHz(ic) = hz(ix,-1)
!        gm%cvHx(ic) = hx(ix,-1)
!        gm%cvQgam(ic,0) = qz(ix,-1,1)
!        gm%cvQgam(ic,1) = -qz(ix,-1,0)
!        gm%cvVol(ic) = vol(ix,-1)
!        gm%cvOnedBsq(ic) = OnedBsq(ix,-1)
!        gm%cvBzb(ic) = bzb(ix,-1)

!        ic = nCi + nx + ix+1
!        gm%cvBb(ic,0:3) = bb(ix,ny,0:3)
!        gm%cvX(ic) = 0.25_R8*sum(crx(ix,ny,0:3))
!        gm%cvY(ic) = 0.25_R8*sum(cry(ix,ny,0:3))
!        gm%cvSz(ic) = gs(ix,ny,2)
!        gm%cvHz(ic) = hz(ix,ny)
!        gm%cvHx(ic) = hx(ix,ny)
!        gm%cvQgam(ic,0) = qz(ix,ny,1)
!        gm%cvQgam(ic,1) = -qz(ix,ny,0)
!        gm%cvVol(ic) = vol(ix,ny)
!        gm%cvOnedBsq(ic) = OnedBsq(ix,ny)
!        gm%cvBzb(ic) = bzb(ix,ny)
!      end do
 
!      do iy = 0, ny-1
!        ic = nCi + 2*nx + iy+1
!        gm%cvBb(ic,0:3) = bb(-1,iy,0:3)
!        gm%cvX(ic) = 0.25_R8*sum(crx(-1,iy,0:3))
!        gm%cvY(ic) = 0.25_R8*sum(cry(-1,iy,0:3))
!        gm%cvSz(ic) = gs(-1,iy,2)
!        gm%cvHz(ic) = hz(-1,iy)
!        gm%cvHx(ic) = hx(-1,iy)
!        gm%cvQgam(ic,0) = qz(-1,iy,1)
!        gm%cvQgam(ic,1) = -qz(-1,iy,0)
!        gm%cvVol(ic) = vol(-1,iy)
!        gm%cvOnedBsq(ic) = OnedBsq(-1,iy)
!        gm%cvBzb(ic) = bzb(-1,iy)

!        ic = nCi + 2*nx + ny + iy+1
!        gm%cvBb(ic,0:3) = bb(nx,iy,0:3)
!        gm%cvX(ic) = 0.25_R8*sum(crx(nx,iy,0:3))
!        gm%cvY(ic) = 0.25_R8*sum(cry(nx,iy,0:3))
!        gm%cvSz(ic) = gs(nx,iy,2)
!        gm%cvHz(ic) = hz(nx,iy)
!        gm%cvHx(ic) = hx(nx,iy)
!        gm%cvQgam(ic,0) = qz(nx,iy,1)
!        gm%cvQgam(ic,1) = -qz(nx,iy,0)
!        gm%cvVol(ic) = vol(nx,iy)
!        gm%cvOnedBsq(ic) = OnedBsq(nx,iy)
!        gm%cvBzb(ic) = bzb(nx,iy)
!      end do

! x faces
       do iy = 0, ny-1
         do ix = 1, nx-1
!           ifc = ix+1 + iy*(nx+1)
           if(mpg%imapFcx(ix,iy)/=0) then
             ifc = mpg%imapFcx(ix,iy)
             gm%fcBb(ifc,0:3) = wbbl(ix,iy,0:3)
             gm%fcS(ifc) = gs(ix,iy,0)
             gm%fcHc(ifc,1) = 0.5_R8*hx(leftix(ix,iy),leftiy(ix,iy))
             gm%fcHc(ifc,2) = 0.5_R8*hx(ix,iy)
             gm%fcHt(ifc)   = 
     .         0.5_R8*(hy(leftix(ix,iy),leftiy(ix,iy))+hy(ix,iy))             !cwdk: computation could be made exact
             gm%fcVol(ifc,1) = 0.5_R8*vol(leftix(ix,iy),leftiy(ix,iy))
             gm%fcVol(ifc,2) = 0.5_R8*vol(ix,iy)
             gm%fcQgam(ifc,0) = qc(ix,iy,0)
             gm%fcQgam(ifc,1) = -qs(ix,iy,0)
           gm%fcQalf(ifc,0) = gm%fcQgam(ifc,0)*
     .       gm%fcBb(ifc,0)/abs(gm%fcBb(ifc,0))
           gm%fcQalf(ifc,1) = gm%fcQgam(ifc,1)*
     .       gm%fcBb(ifc,0)/abs(gm%fcBb(ifc,0))
           gm%fcQbet(ifc,0) = 1.0_R8*
     .       gm%fcBb(ifc,0)/abs(gm%fcBb(ifc,0))
             gm%fcQbet(ifc,1) = 0.0_R8
             gm%fcPbs(ifc) = pbs(ix,iy,0)
             gm%fcPbshz(ifc) = pbshz(ix,iy,0)
           end if
         end do

! left boundary face (to guard cells)
         ix = 0
!         ifc = ix+1 + iy*(nx+1)
         if(mpg%imapFcx(ix,iy)/=0) then
           ifc = mpg%imapFcx(ix,iy)
           gm%fcBb(ifc,0:3) = wbbl(ix,iy,0:3)
           gm%fcS(ifc) = gs(ix,iy,0)
           gm%fcHc(ifc,1) = 0.5_R8*hx(leftix(ix,iy),leftiy(ix,iy))
           gm%fcHc(ifc,2) = 0.5_R8*hx(ix,iy)
           gm%fcHt(ifc)   = hy(leftix(ix,iy),leftiy(ix,iy))                   !cwdk: computation could be made exact
           gm%fcVol(ifc,1) = 0.5_R8*vol(leftix(ix,iy),leftiy(ix,iy))
           gm%fcVol(ifc,2) = 0.5_R8*vol(ix,iy)
           gm%fcQgam(ifc,0) = qc(ix,iy,0)
           gm%fcQgam(ifc,1) = -qs(ix,iy,0)
         gm%fcQalf(ifc,0) = gm%fcQgam(ifc,0)*
     .       gm%fcBb(ifc,0)/abs(gm%fcBb(ifc,0))
         gm%fcQalf(ifc,1) = gm%fcQgam(ifc,1)*
     .       gm%fcBb(ifc,0)/abs(gm%fcBb(ifc,0))
         gm%fcQbet(ifc,0) = 1.0_R8*
     .       gm%fcBb(ifc,0)/abs(gm%fcBb(ifc,0))
           gm%fcQbet(ifc,1) = 0.0_R8
           gm%fcPbs(ifc) = pbs(ix,iy,0)
           gm%fcPbshz(ifc) = pbshz(ix,iy,0)
         end if
! right guard cells
         ix = nx
!         ifc = ix+1 + iy*(nx+1)
         if(mpg%imapFcx(ix,iy)/=0) then 
           ifc = mpg%imapFcx(ix,iy)
           gm%fcBb(ifc,0:3) = wbbl(ix,iy,0:3)
           gm%fcS(ifc) = gs(ix,iy,0)
           gm%fcHc(ifc,1) = 0.5_R8*hx(leftix(ix,iy),leftiy(ix,iy))
           gm%fcHc(ifc,2) = 0.5_R8*hx(ix,iy)
           gm%fcHt(ifc)   = hy(ix,iy)                                         !cwdk: computation could be made exact
           gm%fcVol(ifc,1) = 0.5_R8*vol(leftix(ix,iy),leftiy(ix,iy))
           gm%fcVol(ifc,2) = 0.5_R8*vol(ix,iy)
           gm%fcQgam(ifc,0) = qc(ix,iy,0)
           gm%fcQgam(ifc,1) = -qs(ix,iy,0)
         gm%fcQalf(ifc,0) = gm%fcQgam(ifc,0)*
     .       gm%fcBb(ifc,0)/abs(gm%fcBb(ifc,0))
         gm%fcQalf(ifc,1) = gm%fcQgam(ifc,1)*
     .       gm%fcBb(ifc,0)/abs(gm%fcBb(ifc,0))
         gm%fcQbet(ifc,0) = 1.0_R8*
     .       gm%fcBb(ifc,0)/abs(gm%fcBb(ifc,0))
           gm%fcQbet(ifc,1) = 0.0_R8
           gm%fcPbs(ifc) = pbs(ix,iy,0)
           gm%fcPbshz(ifc) = pbshz(ix,iy,0)
         endif
       end do

! y faces
!       allocate (qs(-1:nx,-1:ny,0:1))
!       call intfacev(nx, ny, bottomix, bottomiy, vol, qz(-1,-1,0),
!     .               qs(-1,-1,1))
       do ix = 0, nx-1
         do iy = 1, ny-1
!           ifc = nfx + iy+1 + ix*(ny+1)
           if(mpg%imapFcy(ix,iy)/=0) then
             ifc = mpg%imapFcy(ix,iy)
             gm%fcBb(ifc,0:3) = wbbv(ix,iy,0:3)
             gm%fcS(ifc) = gs(ix,iy,1)
             gm%fcHc(ifc,1) = 0.5_R8*hy(bottomix(ix,iy),bottomiy(ix,iy))       !cwdk: change orientation ??
             gm%fcHc(ifc,2) = 0.5_R8*hy(ix,iy)
             gm%fcHt(ifc)   = 
     .         0.5_R8*(hx(bottomix(ix,iy),bottomiy(ix,iy)) + hx(ix,iy))       !cwdk: computation could be made exact
             gm%fcVol(ifc,1) = 0.5_R8*vol(bottomix(ix,iy),
     .                                    bottomiy(ix,iy))
             gm%fcVol(ifc,2) = 0.5_R8*vol(ix,iy)
             gm%fcQgam(ifc,0) = qc(ix,iy,1)
             gm%fcQgam(ifc,1) = qs(ix,iy,1)
             gm%fcQalf(ifc,0) = 0.0_R8
           gm%fcQalf(ifc,1) = 1.0_R8*
     .       gm%fcBb(ifc,0)/abs(gm%fcBb(ifc,0))
           gm%fcQbet(ifc,0) = gm%fcQgam(ifc,1)*
     .       gm%fcBb(ifc,0)/abs(gm%fcBb(ifc,0))
           gm%fcQbet(ifc,1) = -gm%fcQgam(ifc,0)*
     .       gm%fcBb(ifc,0)/abs(gm%fcBb(ifc,0))
             gm%fcPbs(ifc) = pbs(ix,iy,1)
             gm%fcPbshz(ifc) = pbshz(ix,iy,1)
           end if
         end do

! lower guard cells
         iy = 0
!         ifc = nfx + 1 + ix*(ny+1)
         if(mpg%imapFcy(ix,iy)/=0) then
           ifc = mpg%imapFcy(ix,iy)  
           gm%fcBb(ifc,0:3) = wbbv(ix,iy,0:3)
           gm%fcS(ifc) = gs(ix,iy,1)
           gm%fcHc(ifc,1) = 0.5_R8*hy(bottomix(ix,iy),bottomiy(ix,iy))
           gm%fcHc(ifc,2) = 0.5_R8*hy(ix,iy)
           gm%fcHt(ifc)   = hx(bottomix(ix,iy),bottomiy(ix,iy))               !cwdk: computation could be made exact
           gm%fcVol(ifc,1) = 0.5_R8*vol(bottomix(ix,iy),bottomiy(ix,iy))
           gm%fcVol(ifc,2) = 0.5_R8*vol(ix,iy)
           gm%fcQgam(ifc,0) = qc(ix,iy,1)
           gm%fcQgam(ifc,1) = qs(ix,iy,1)
           gm%fcQalf(ifc,0) = 0.0_R8
         gm%fcQalf(ifc,1) = 1.0_R8*
     .     gm%fcBb(ifc,0)/abs(gm%fcBb(ifc,0))
         gm%fcQbet(ifc,0) = gm%fcQgam(ifc,1)*
     .     gm%fcBb(ifc,0)/abs(gm%fcBb(ifc,0))
         gm%fcQbet(ifc,1) = -gm%fcQgam(ifc,0)*
     .     gm%fcBb(ifc,0)/abs(gm%fcBb(ifc,0))
           gm%fcPbs(ifc) = pbs(ix,iy,1)
           gm%fcPbshz(ifc) = pbshz(ix,iy,1)
         endif
! upper guard cells
         iy = ny
!         ifc = nfx + (ix+1)*(ny+1)
         if(mpg%imapFcy(ix,iy)/=0) then
           ifc = mpg%imapFcy(ix,iy)
           gm%fcBb(ifc,0:3) = wbbv(ix,iy,0:3)
           gm%fcS(ifc) = gs(ix,iy,1)
           gm%fcHc(ifc,1) = 0.5_R8*hy(bottomix(ix,iy),bottomiy(ix,iy))
           gm%fcHc(ifc,2) = 0.5_R8*hy(ix,iy)
           gm%fcHt(ifc)   = hx(ix,iy)                                         !cwdk: computation could be made exact
           gm%fcVol(ifc,1) = 0.5_R8*vol(bottomix(ix,iy),bottomiy(ix,iy))
           gm%fcVol(ifc,2) = 0.5_R8*vol(ix,iy)
           gm%fcQgam(ifc,0) = qc(ix,iy,1)
           gm%fcQgam(ifc,1) = qs(ix,iy,1)
           gm%fcQalf(ifc,0) = 0.0_R8
         gm%fcQalf(ifc,1) = 1.0_R8*
     .     gm%fcBb(ifc,0)/abs(gm%fcBb(ifc,0))
         gm%fcQbet(ifc,0) = gm%fcQgam(ifc,1)*
     .     gm%fcBb(ifc,0)/abs(gm%fcBb(ifc,0))
         gm%fcQbet(ifc,1) = -gm%fcQgam(ifc,0)*
     .     gm%fcBb(ifc,0)/abs(gm%fcBb(ifc,0))
           gm%fcPbs(ifc) = pbs(ix,iy,1)
           gm%fcPbshz(ifc) = pbshz(ix,iy,1)
         end if
       end do      

! vertices
       do iy = 0, ny
         do ix = 0, nx
!           ivx = ix+1 + iy*(nx+1)
           if(mpg%imapVx(ix,iy)/=0) then
             ivx = mpg%imapVx(ix,iy)
             gm%vxX(ivx) = crx(ix,iy,0)
             gm%vxY(ivx) = cry(ix,iy,0)
             gm%vxFfbz(ivx) = ffbz(ix,iy,1)
             gm%vxFpsi(ivx) = fpsi(ix,iy,1)
             do i = 1, mpg%vxCvP(ivx,2)
               ic = mpg%vxCv(mpg%vxCvP(ivx,1) + i - 1)
               gm%vxVol(mpg%vxCvP(ivx,1) + i - 1) = 0.25*gm%cvVol(ic)
             end do
           end if
         end do
       end do

 
! connection length for core region

      gm%fsConn = 0.0_R8
      do i = 1, nFs 
       ic = mpg%fsCv(mpg%fsCvP(i,1))
       if (mpg%cvOnClosedSurface(ic)) then        
        do j=mpg%fsFcP(i,1), mpg%fsFcP(i,1)+mpg%fsFcP(i,2)-1
          ifc = mpg%fsFc(j)
cwdk  there is an implicit assumption of field-alignment of fcHc(ifc,:) here!
          gm%fsconn(i) = gm%fsconn(i) + 
     .                   gm%fcHc(ifc,1)*gm%cvBb(mpg%fcCv(ifc,1),3)/
     .                              abs(gm%cvBb(mpg%fcCv(ifc,1),0)) + 
     .                   gm%fcHc(ifc,2)*gm%cvBb(mpg%fcCv(ifc,2),3)/
     .                              abs(gm%cvBb(mpg%fcCv(ifc,2),0))
        end do
       end if
      end do      


*  .. unit vector of the magnetic field in X,Y,Z coorinates

cwdk  cell centers
cwdk  assumption: poloidal direction is direction of line between 
cwdk              poloidal face centers (E/W)
      do iy = 0, ny-1
        do ix = 0, nx-1
          ic = (ix+1) + iy*nx

          ! Construct unit vector in poloidal direction:
          dux = (crx(ix,iy,1)+crx(ix,iy,3))/2.0_R8-
     &          (crx(ix,iy,0)+crx(ix,iy,2))/2.0_R8
          duy = (cry(ix,iy,1)+cry(ix,iy,3))/2.0_R8-
     &          (cry(ix,iy,0)+cry(ix,iy,2))/2.0_R8

          ! length should be hx; recomputed for safety
          du  = sqrt(dux**2+duy**2)

          ! implicit assumption: gm%cvBb(ic,1) = 0.0
          gm%cvEb(ic,0) = gm%cvBb(ic,0)/gm%cvBb(ic,3)*dux/du
          gm%cvEb(ic,1) = gm%cvBb(ic,0)/gm%cvBb(ic,3)*duy/du
          gm%cvEb(ic,2) = gm%cvBb(ic,2)/gm%cvBb(ic,3)

        enddo
      enddo

      do ix = 0, nx-1
        ! bottom guard cells
        ic = nCi + ix+1
        iy = -1
        
        ! Construct unit vector in poloidal direction:
        dux = (crx(ix,iy,1)+crx(ix,iy,3))/2.0_R8-
     &        (crx(ix,iy,0)+crx(ix,iy,2))/2.0_R8
        duy = (cry(ix,iy,1)+cry(ix,iy,3))/2.0_R8-
     &        (cry(ix,iy,0)+cry(ix,iy,2))/2.0_R8

        ! length should be hx; recomputed for safety
        du  = sqrt(dux**2+duy**2)

        ! implicit assumption: gm%cvBb(ic,1) = 0.0
        gm%cvEb(ic,0) = gm%cvBb(ic,0)/gm%cvBb(ic,3)*dux/du
        gm%cvEb(ic,1) = gm%cvBb(ic,0)/gm%cvBb(ic,3)*duy/du
        gm%cvEb(ic,2) = gm%cvBb(ic,2)/gm%cvBb(ic,3)

        ! top guard cells
        ic = nCi + nx + ix+1
        iy = ny

        ! Construct unit vector in poloidal direction:
        dux = (crx(ix,iy,1)+crx(ix,iy,3))/2.0_R8-
     &        (crx(ix,iy,0)+crx(ix,iy,2))/2.0_R8
        duy = (cry(ix,iy,1)+cry(ix,iy,3))/2.0_R8-
     &        (cry(ix,iy,0)+cry(ix,iy,2))/2.0_R8

        ! length should be hx; recomputed for safety
        du  = sqrt(dux**2+duy**2)

        ! implicit assumption: gm%cvBb(ic,1) = 0.0
        gm%cvEb(ic,0) = gm%cvBb(ic,0)/gm%cvBb(ic,3)*dux/du
        gm%cvEb(ic,1) = gm%cvBb(ic,0)/gm%cvBb(ic,3)*duy/du
        gm%cvEb(ic,2) = gm%cvBb(ic,2)/gm%cvBb(ic,3)

      end do
 
      do iy = 0, ny-1
        ic = nCi + 2*nx + iy+1
        ix = -1

        ! Construct unit vector in poloidal direction:
        dux = (crx(ix,iy,1)+crx(ix,iy,3))/2.0_R8-
     &        (crx(ix,iy,0)+crx(ix,iy,2))/2.0_R8
        duy = (cry(ix,iy,1)+cry(ix,iy,3))/2.0_R8-
     &        (cry(ix,iy,0)+cry(ix,iy,2))/2.0_R8

        ! length should be hx; recomputed for safety
        du  = sqrt(dux**2+duy**2)

        ! implicit assumption: gm%cvBb(ic,1) = 0.0
        gm%cvEb(ic,0) = gm%cvBb(ic,0)/gm%cvBb(ic,3)*dux/du
        gm%cvEb(ic,1) = gm%cvBb(ic,0)/gm%cvBb(ic,3)*duy/du
        gm%cvEb(ic,2) = gm%cvBb(ic,2)/gm%cvBb(ic,3)

        ic = nCi + 2*nx + ny + iy+1
        ix = nx

        ! Construct unit vector in poloidal direction:
        dux = (crx(ix,iy,1)+crx(ix,iy,3))/2.0_R8-
     &        (crx(ix,iy,0)+crx(ix,iy,2))/2.0_R8
        duy = (cry(ix,iy,1)+cry(ix,iy,3))/2.0_R8-
     &        (cry(ix,iy,0)+cry(ix,iy,2))/2.0_R8

        ! length should be hx; recomputed for safety
        du  = sqrt(dux**2+duy**2)

        ! implicit assumption: gm%cvBb(ic,1) = 0.0
        gm%cvEb(ic,0) = gm%cvBb(ic,0)/gm%cvBb(ic,3)*dux/du
        gm%cvEb(ic,1) = gm%cvBb(ic,0)/gm%cvBb(ic,3)*duy/du
        gm%cvEb(ic,2) = gm%cvBb(ic,2)/gm%cvBb(ic,3)

      end do

      ! normalize
      gm%cvEb(:,0) = gm%cvEb(:,0)/
     &  sqrt(gm%cvEb(:,0)**2 + gm%cvEb(:,1)**2 + gm%cvEb(:,2)**2)
      gm%cvEb(:,1) = gm%cvEb(:,1)/
     &  sqrt(gm%cvEb(:,0)**2 + gm%cvEb(:,1)**2 + gm%cvEb(:,2)**2)
      gm%cvEb(:,2) = gm%cvEb(:,2)/
     &  sqrt(gm%cvEb(:,0)**2 + gm%cvEb(:,1)**2 + gm%cvEb(:,2)**2)

      end subroutine prep_geometry

**********************************************************************
      
      subroutine prep_plasma (nx, ny, ns, nCv, nCi, nFc, mpg, st)
      
      implicit none
      integer, intent(in) :: nx, ny, ns, nCi, nCv, nFc
      type(mapping), intent(in) :: mpg
      type(B2State), intent(inout) :: st
      integer :: ix, iy, ic, ifc, nfx

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
            st % dv % fna_he(ifc,0:1,0:ns-1)=fna_he(ix,iy,0,0:1,0:ns-1)
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
          st % dv % fne_eir(ifc,0:1) = 
     .              fne_eir(ix,iy,1,0:1)
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
      
      subroutine to_struct_plasma (ns, mpg, st)
      
      implicit none
      integer, intent(in) :: ns
      type(mapping), intent(in) :: mpg
      type(B2State), intent(in) :: st
      integer :: nx, ny, nCi, nCv, nFc
      integer :: ix, iy, ic, ifc, nfx
      integer :: nout_st, io, idum(0:9)
      character lblst*120, lblps*120


      nCv = mpg%nCv
      nCi = mpg%nCi
      nFc = mpg%nFc

      if (mpg.isClassicalGrid == 1) then
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

      nfx = (nx+1)*ny

! cell volumes 
      do iy = 0, ny-1
        do ix = 0, nx-1
          ic = (ix+1) + iy*nx
! B2plasma
          na(ix,iy,0:ns-1) = st % pl % na(ic,0:ns-1)
          ua(ix,iy,0:ns-1) = st % pl % ua(ic,0:ns-1)
          po(ix,iy)        = st % pl % po(ic)
          te(ix,iy)        = st % pl % te(ic) 
          ti(ix,iy)        = st % pl % ti(ic) 

! B2coeff
          ceqp(ix,iy) = st % co % ceqp(ic)

! from b2coeffwork
          dna0(ix,iy,0:ns-1)     = st % co % dna0(ic,0:ns-1)
          dpa0(ix,iy,0:ns-1)     = st % co % dpa0(ic,0:ns-1)
          vla0(ix,iy,0:1,0:ns-1) = st % co % vla0(ic,0:1,0:ns-1)
          vsa0(ix,iy,0:ns-1)     = st % co % vsa0(ic,0:ns-1)
          hci0(ix,iy)            = st % co % hci0(ic)
          hcib(ix,iy,0:ns-1)     = st % co % hcib(ic,0:ns-1)
          hce0(ix,iy)            = st % co % hce0(ic) 
          sig0(ix,iy)            = st % co % sig0(ic)
          alf0(ix,iy)            = st % co % alf0(ic)
      
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
!         resht(ix,iy)        = st % dv % resht(ic)
          resmt(ix,iy)        = st % dv % resmt(ic)
      
! corrections
          corua(ix,iy,0:ns-1) = st % dv % corua(ic,0:ns-1)
          corpa(ix,iy,0:ns-1) = st % dv % corpa(ic,0:ns-1)
          corut(ix,iy)        = st % dv % corut(ic)
          corpo(ix,iy)        = st % dv % corpo(ic)
!         cortt(ix,iy)        = st % dv % cortt(ic)
          corte(ix,iy)        = st % dv % corte(ic)
          corti(ix,iy)        = st % dv % corti(ic)
      
! coupling coefficients
          pcca(ix,iy,0,0,0:ns-1) = st % dv % pcca(ic,0,0:ns-1)
          pccm(ix,iy,0,0)        = st % dv % pccm(ic,0)
          pcca(ix,iy,1,1,0:ns-1) = st % dv % pcca(ic,1,0:ns-1)
          pccm(ix,iy,1,1)        = st % dv % pccm(ic,1)

! B2Source
          sch(ix,iy,0:3)        = st % sr % sch(ic,0:3)
          she(ix,iy,0:3)        = st % sr % she(ic,0:3)
          shi(ix,iy,0:3)        = st % sr % shi(ic,0:3)
          smo(ix,iy,0:3,0:ns-1) = st % sr % smo(ic,0:3,0:ns-1)
          smq(ix,iy,0:3,0:ns-1) = st % sr % smq(ic,0:3,0:ns-1)
          sna(ix,iy,0:1,0:ns-1) = st % sr % sna(ic,0:1,0:ns-1)

! B2rates
          rlqa(ix,iy,0:1,0:ns-1) = st % rt % rlqa(ic,0:1,0:ns-1)
          rlrd(ix,iy,0:1,0:ns-1) = st % rt % rlrd(ic,0:1,0:ns-1)
          rlqa(ix,iy,0:1,0:ns-1) = st % rt % rlqa(ic,0:1,0:ns-1)
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

        end do
      end do
      
! gueard cells iy = -1 and iy = ny
      do ix = 0, nx-1
        ic = nCi + ix+1
        iy = -1
! B2plasma
        na(ix,iy,0:ns-1) = st % pl % na(ic,0:ns-1)
        ua(ix,iy,0:ns-1) = st % pl % ua(ic,0:ns-1)
        po(ix,iy)        = st % pl % po(ic)
        te(ix,iy)        = st % pl % te(ic)
        ti(ix,iy)        = st % pl % ti(ic)

! B2coeff
        ceqp(ix,iy) = st % co % ceqp(ic)

! from b2coeffwork
        dna0(ix,iy,0:ns-1)     = st % co % dna0(ic,0:ns-1)
        dpa0(ix,iy,0:ns-1)     = st % co % dpa0(ic,0:ns-1)
        vla0(ix,iy,0:1,0:ns-1) = st % co % vla0(ic,0:1,0:ns-1)
        vsa0(ix,iy,0:ns-1)     = st % co % vsa0(ic,0:ns-1)
        hci0(ix,iy)            = st % co % hci0(ic)
        hcib(ix,iy,0:ns-1)     = st % co % hcib(ic,0:ns-1)
        hce0(ix,iy)            = st % co % hce0(ic)
        sig0(ix,iy)            = st % co % sig0(ic)
        alf0(ix,iy)            = st % co % alf0(ic)
        
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
!       resht(ix,iy)        = st % dv % resht(ic)
        resmt(ix,iy)        = st % dv % resmt(ic)
      
! corrections
        corua(ix,iy,0:ns-1) = st % dv % corua(ic,0:ns-1)
        corpa(ix,iy,0:ns-1) = st % dv % corpa(ic,0:ns-1)
        corut(ix,iy)        = st % dv % corut(ic)
        corpo(ix,iy)        = st % dv % corpo(ic)
!       cortt(ix,iy)        = st % dv % cortt(ic)
        corte(ix,iy)        = st % dv % corte(ic)
        corti(ix,iy)        = st % dv % corti(ic)
      
! coupling coefficients
        pcca(ix,iy,0,0,0:ns-1) = st % dv % pcca(ic,0,0:ns-1)
        pccm(ix,iy,0,0)        = st % dv % pccm(ic,0)
        pcca(ix,iy,1,1,0:ns-1) = st % dv % pcca(ic,1,0:ns-1)
        pccm(ix,iy,1,1)        = st % dv % pccm(ic,1)
        
! B2SourceWork
        sch(ix,iy,0:3)        = st % sr % sch(ic,0:3)
        she(ix,iy,0:3)        = st % sr % she(ic,0:3)
        shi(ix,iy,0:3)        = st % sr % shi(ic,0:3)
        smo(ix,iy,0:3,0:ns-1) = st % sr % smo(ic,0:3,0:ns-1)
        smq(ix,iy,0:3,0:ns-1) = st % sr % smq(ic,0:3,0:ns-1)
        sna(ix,iy,0:1,0:ns-1) = st % sr % sna(ic,0:1,0:ns-1)
       
! B2rates
        rlqa(ix,iy,0:1,0:ns-1) = st % rt % rlqa(ic,0:1,0:ns-1)
        rlrd(ix,iy,0:1,0:ns-1) = st % rt % rlrd(ic,0:1,0:ns-1)
        rlqa(ix,iy,0:1,0:ns-1) = st % rt % rlqa(ic,0:1,0:ns-1)
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

        ic = nCi + nx + ix+1
        iy = ny
! B2plasma
        na(ix,iy,0:ns-1) = st % pl % na(ic,0:ns-1)
        ua(ix,iy,0:ns-1) = st % pl % ua(ic,0:ns-1)
        po(ix,iy)        = st % pl % po(ic)
        te(ix,iy)        = st % pl % te(ic)
        ti(ix,iy)        = st % pl % ti(ic)
        
! B2coeff
        ceqp(ix,iy) = st % co % ceqp(ic)
        
! from b2coeffwork
        dna0(ix,iy,0:ns-1)     = st % co % dna0(ic,0:ns-1)
        dpa0(ix,iy,0:ns-1)     = st % co % dpa0(ic,0:ns-1)
        vla0(ix,iy,0:1,0:ns-1) = st % co % vla0(ic,0:1,0:ns-1)
        vsa0(ix,iy,0:ns-1)     = st % co % vsa0(ic,0:ns-1)
        hci0(ix,iy)            = st % co % hci0(ic)
        hcib(ix,iy,0:ns-1)     = st % co % hcib(ic,0:ns-1)
        hce0(ix,iy)            = st % co % hce0(ic)
        sig0(ix,iy)            = st % co % sig0(ic)
        alf0(ix,iy)            = st % co % alf0(ic)
        
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
!       resht(ix,iy)        = st % dv % resht(ic)
        resmt(ix,iy)        = st % dv % resmt(ic)
      
! corrections
        corua(ix,iy,0:ns-1) = st % dv % corua(ic,0:ns-1)
        corpa(ix,iy,0:ns-1) = st % dv % corpa(ic,0:ns-1)
        corut(ix,iy)        = st % dv % corut(ic)
        corpo(ix,iy)        = st % dv % corpo(ic)
!       cortt(ix,iy)        = st % dv % cortt(ic)
        corte(ix,iy)        = st % dv % corte(ic)
        corti(ix,iy)        = st % dv % corti(ic)
      
! coupling coefficients
        pcca(ix,iy,0,0,0:ns-1) = st % dv % pcca(ic,0,0:ns-1)
        pccm(ix,iy,0,0)        = st % dv % pccm(ic,0)
        pcca(ix,iy,1,1,0:ns-1) = st % dv % pcca(ic,1,0:ns-1)
        pccm(ix,iy,1,1)        = st % dv % pccm(ic,1)

! B2SourceWork
        sch(ix,iy,0:3)        = st % sr % sch(ic,0:3)
        she(ix,iy,0:3)        = st % sr % she(ic,0:3)
        shi(ix,iy,0:3)        = st % sr % shi(ic,0:3)
        smo(ix,iy,0:3,0:ns-1) = st % sr % smo(ic,0:3,0:ns-1)
        smq(ix,iy,0:3,0:ns-1) = st % sr % smq(ic,0:3,0:ns-1)
        sna(ix,iy,0:1,0:ns-1) = st % sr % sna(ic,0:3,0:ns-1)

! B2rates
        rlqa(ix,iy,0:1,0:ns-1) = st % rt % rlqa(ic,0:1,0:ns-1)
        rlrd(ix,iy,0:1,0:ns-1) = st % rt % rlrd(ic,0:1,0:ns-1)
        rlqa(ix,iy,0:1,0:ns-1) = st % rt % rlqa(ic,0:1,0:ns-1)
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
      end do

! guard cells ix = -1 and ix = nx
      do iy = 0, ny-1
        ic = nCi + 2*nx + iy+1
        ix = -1
! B2plasma
        na(ix,iy,0:ns-1) = st % pl % na(ic,0:ns-1)
        ua(ix,iy,0:ns-1) = st % pl % ua(ic,0:ns-1) 
        po(ix,iy)        = st % pl % po(ic)
        te(ix,iy)        = st % pl % te(ic)
        ti(ix,iy)        = st % pl % ti(ic)
        
! B2coeff
        ceqp(ix,iy) = st % co % ceqp(ic)
        
! from b2coeffwork
        dna0(ix,iy,0:ns-1)     = st % co % dna0(ic,0:ns-1)
        dpa0(ix,iy,0:ns-1)     = st % co % dpa0(ic,0:ns-1)
        vla0(ix,iy,0:1,0:ns-1) = st % co % vla0(ic,0:1,0:ns-1)
        vsa0(ix,iy,0:ns-1)     = st % co % vsa0(ic,0:ns-1)
        hci0(ix,iy)            = st % co % hci0(ic)
        hcib(ix,iy,0:ns-1)     = st % co % hcib(ic,0:ns-1)
        hce0(ix,iy)            = st % co % hce0(ic)
        sig0(ix,iy)            = st % co % sig0(ic)
        alf0(ix,iy)            = st % co % alf0(ic)
        
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
!       resht(ix,iy)        = st % dv % resht(ic)
        resmt(ix,iy)        = st % dv % resmt(ic)
      
! corrections
        corua(ix,iy,0:ns-1) = st % dv % corua(ic,0:ns-1)
        corpa(ix,iy,0:ns-1) = st % dv % corpa(ic,0:ns-1)
        corut(ix,iy)        = st % dv % corut(ic)
        corpo(ix,iy)        = st % dv % corpo(ic)
!       cortt(ix,iy)        = st % dv % cortt(ic)
        corte(ix,iy)        = st % dv % corte(ic)
        corti(ix,iy)        = st % dv % corti(ic)
      
! coupling coefficients
        pcca(ix,iy,0,0,0:ns-1) = st % dv % pcca(ic,0,0:ns-1)
        pccm(ix,iy,0,0)        = st % dv % pccm(ic,0)
        pcca(ix,iy,1,1,0:ns-1) = st % dv % pcca(ic,1,0:ns-1)
        pccm(ix,iy,1,1)        = st % dv % pccm(ic,1)

! B2SourceWork
        sch(ix,iy,0:3)        = st % sr % sch(ic,0:3)
        she(ix,iy,0:3)        = st % sr % she(ic,0:3)
        shi(ix,iy,0:3)        = st % sr % shi(ic,0:3)
        smo(ix,iy,0:3,0:ns-1) = st % sr % smo(ic,0:3,0:ns-1)
        smq(ix,iy,0:3,0:ns-1) = st % sr % smq(ic,0:3,0:ns-1)
        sna(ix,iy,0:1,0:ns-1) = st % sr % sna(ic,0:1,0:ns-1)

! B2rates
        rlqa(ix,iy,0:1,0:ns-1) = st % rt % rlqa(ic,0:1,0:ns-1)
        rlrd(ix,iy,0:1,0:ns-1) = st % rt % rlrd(ic,0:1,0:ns-1)
        rlqa(ix,iy,0:1,0:ns-1) = st % rt % rlqa(ic,0:1,0:ns-1)
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


        ic = nCi + 2*nx + ny + iy+1
        ix = nx
! B2plasma
        na(ix,iy,0:ns-1) = st % pl % na(ic,0:ns-1)
        ua(ix,iy,0:ns-1) = st % pl % ua(ic,0:ns-1)
        po(ix,iy)        = st % pl % po(ic)
        te(ix,iy)        = st % pl % te(ic)
        ti(ix,iy)        = st % pl % ti(ic)
        
! B2coeff
        ceqp(ix,iy) = st % co % ceqp(ic)
        
! from b2coeffwork
        dna0(ix,iy,0:ns-1)     = st % co % dna0(ic,0:ns-1)
        dpa0(ix,iy,0:ns-1)     = st % co % dpa0(ic,0:ns-1)
        vla0(ix,iy,0:1,0:ns-1) = st % co % vla0(ic,0:1,0:ns-1)
        vsa0(ix,iy,0:ns-1)     = st % co % vsa0(ic,0:ns-1)
        hci0(ix,iy)            = st % co % hci0(ic)
        hcib(ix,iy,0:ns-1)     = st % co % hcib(ic,0:ns-1)
        hce0(ix,iy)            = st % co % hce0(ic)
        sig0(ix,iy)            = st % co % sig0(ic)
        alf0(ix,iy)            = st % co % alf0(ic)
        
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
!       resht(ix,iy)        = st % dv % resht(ic)
        resmt(ix,iy)        = st % dv % resmt(ic)
      
! corrections
        corua(ix,iy,0:ns-1) = st % dv % corua(ic,0:ns-1)
        corpa(ix,iy,0:ns-1) = st % dv % corpa(ic,0:ns-1)
        corut(ix,iy)        = st % dv % corut(ic)
        corpo(ix,iy)        = st % dv % corpo(ic)
!       cortt(ix,iy)        = st % dv % cortt(ic)
        corte(ix,iy)        = st % dv % corte(ic)
        corti(ix,iy)        = st % dv % corti(ic)
      
! coupling coefficients
        pcca(ix,iy,0,0,0:ns-1) = st % dv % pcca(ic,0,0:ns-1)
        pccm(ix,iy,0,0)        = st % dv % pccm(ic,0)
        pcca(ix,iy,1,1,0:ns-1) = st % dv % pcca(ic,1,0:ns-1)
        pccm(ix,iy,1,1)        = st % dv % pccm(ic,1)

! B2SourceWork
        sch(ix,iy,0:3)        = st % sr % sch(ic,0:3)
        she(ix,iy,0:3)        = st % sr % she(ic,0:3)
        shi(ix,iy,0:3)        = st % sr % shi(ic,0:3)
        smo(ix,iy,0:3,0:ns-1) = st % sr % smo(ic,0:3,0:ns-1)
        smq(ix,iy,0:3,0:ns-1) = st % sr % smq(ic,0:3,0:ns-1)
        sna(ix,iy,0:1,0:ns-1) = st % sr % sna(ic,0:3,0:ns-1)

! B2rates
        rlqa(ix,iy,0:1,0:ns-1) = st % rt % rlqa(ic,0:1,0:ns-1)
        rlrd(ix,iy,0:1,0:ns-1) = st % rt % rlrd(ic,0:1,0:ns-1)
        rlqa(ix,iy,0:1,0:ns-1) = st % rt % rlqa(ic,0:1,0:ns-1)
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
          ifc = ix+1 + iy*(nx+1)

! B2coeff
          csig_an(ix,iy,0,0:1) = st % co % csig_an(ifc,0:1)
          calf(ix,iy,0,0:1)    = st % co % calf(ifc,0:1)

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
          cthe(ix,iy,0:ns-1) = st % co % cthe(ifc,0:ns-1)
          cthi(ix,iy,0:ns-1) = st % co % cthi(ifc,0:ns-1)
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
          fhe_mdf(ix,iy,0,0:1) = st % dv % fhe_mdf(ifc,0:1)
          fhet(ix,iy,0,0:1)    = st % dv % fhet(ifc,0:1)
          fhe_eir(ix,iy,0,0:1) = st % dv % fhe_eir(ifc,0:1)
          fhi(ix,iy,0,0:1)     = st % dv % fhi(ifc,0:1)
          fhi_mdf(ix,iy,0,0:1) = st % dv % fhi_mdf(ifc,0:1)
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

        end do
      end do



! y faces

      do ix = 0, nx-1
        do iy = 0, ny
          ifc = nfx + iy+1 + ix*(ny+1)

! B2coeff
          csig_an(ix,iy,1,0:1) = st % co % csig_an(ifc,0:1)
          calf(ix,iy,1,0:1)    = st % co % calf(ifc,0:1)
        
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
          cthe(ix,iy,0:ns-1)          = st % co % cthe(ifc,0:ns-1)
          cthi(ix,iy,0:ns-1)          = st % co % cthi(ifc,0:ns-1)
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
          fhe_mdf(ix,iy,1,0:1) = st % dv % fhe_mdf(ifc,0:1)
          fhePSch(ix,iy,1,0:1) = st % dv % fhePSch(ifc,0:1)
          fhet(ix,iy,1,0:1)    = st % dv % fhet(ifc,0:1)
          fhe_eir(ix,iy,1,0:1) = st % dv % fhe_eir(ifc,0:1)
          fhi(ix,iy,1,0:1)     = st % dv % fhi(ifc,0:1)
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
      open (newunit=nout_st,file='b2fstate_st',iostat=io)
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
     &  ti, po, fna, fhe, fhi, fch, fch_32, fch_52, fch_p, kinrgy,
     &  fna_mdf, fhe_mdf, fhi_mdf, fna_fcor, fna_nodrift, fna_he,
     &  fnaPSch, fhePSch, fhiPSch, fna_eir, fne_eir, fhe_eir, fhi_eir,
     &  fna_32, fna_52, fni_32, fni_52, fne_32, fne_52,
     &  fchdia, fchin, fchvispar, fchvisper, fchvisq, fchinert,
     &  vaecrb, vadia, wadia, veecrb, vedia, floe_noc, floi_noc)

      close(nout_st)

*   ..write b2fplasmf_st in structured format
      open (newunit=nout_st,file='b2fplasmf_st',iostat=io)
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
        integer :: ix, iy, ic, ifc

        ! For now, only create a dummy external state.
        ! Functionality from b2mod_external to be carried over!
        call createB2StateExt( nCv, nFc, ns, st_ext )

      end subroutine prep_plasma_ext




      end module b2us_prep

