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
      use b2mod_time
      use b2mod_version
      use b2mod_b2cmpa
      use b2us_map
      use b2us_geo
      use b2us_plasma
      use b2mod_b2cmfs

      implicit none

      private

      public :: read_mapping, calc_mapping, write_mapping,
     .          read_geometry, prep_geometry, write_geometry,
     .          calc_mapping_bc, calc_mapping_rc,
     .          prep_plasma, prep_plasma_ext,
     .          to_struct_plasma, read_map_geo, write_map_geo,
     .          write_plasma_state, read_plasma_state,
     .          conv_bcparms, write_bcparms

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
      type (mapping), intent (out) :: mpg
      type (geometry), intent (out) :: gm
      integer :: idum(0:9)
      character lblgm*120, idout*80

* ..read the geometry
* ..obtain version number

      call read_mapping (ninp, mpg)
      call read_geometry (ninp, mpg, gm)

!      close (ninp)
      
      end subroutine read_map_geo

***********************************************************************

      subroutine write_map_geo (fname, mpg, gm)
      use b2mod_types
      use b2mod_version
      use b2mod_b2cmfs

      implicit none
      character(*),intent(in) :: fname
      type (mapping), intent (in) :: mpg
      type (geometry), intent (in) :: gm
      integer :: nout, io, idum(0:9)
      character lblgm*120, idout*80

      open (newunit=nout,file=fname,iostat=io)

      if (io /= 0) then
        write (*,*) ' file ', fname, ' not writable'
        call xerrab ('write_map_geo -- no geometry file ')
      end if
* ..set new version number
      b2fgmtry_version = '03.002.000'

* ..read the geometry
* ..obtain version number
      call cfverw (nout, b2fgmtry_version)
      lblgm = 'Test for b2fgmtry_us'
      call cfwuch (nout, 120, lblgm, 'label')
      idum(0) = isymm
      call cfwuin (nout, 1, idum, 'isymm')

      call write_mapping (nout, mpg)
      call write_geometry (nout, mpg, gm)

      close (nout)
      
      end subroutine write_map_geo

***********************************************************************


      subroutine read_mapping (iun, mpg)

      implicit none
      integer, intent(in) :: iun
      integer :: idum(0:9)
      type (mapping), intent (out) :: mpg

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
      call cfruin (iun, 1, idum, 'nVmxCv')
      mpg%nVmxCv = idum(0)
      call cfruin (iun, 1, idum, 'nVmxFc')
      mpg%nVmxFc = idum(0)

      call alloc_mapping (mpg)
        
      call cfruin (iun, mpg%nCv*2,   mpg%cvFcP, 'cvFcP')
      call cfruin (iun, mpg%nCmxFc,  mpg%cvFc,  'cvFc')
      call cfruin (iun, mpg%nFc*2,   mpg%fcCv,  'fcCv')
      call cfruin (iun, mpg%nFc*2,   mpg%fcVx,  'fcVx')
      call cfruin (iun, mpg%nCv*2,   mpg%cvVxP, 'cvVxP')
      call cfruin (iun, mpg%nCmxVx,  mpg%cvVx,  'cvVx')
      call cfruin (iun, mpg%nVx*2,   mpg%vxFcP, 'vxFcP')
      call cfruin (iun, mpg%nVmxFc,  mpg%vxFc,  'vxFc')
      call cfruin (iun, mpg%nVx*2,   mpg%vxCvP, 'vxCvP')
      call cfruin (iun, mpg%nVmxCv,  mpg%vxCv,  'vxCv')
      call cfruin (iun, mpg%nFs*2,   mpg%fsCvP, 'fsCvP')
      call cfruin (iun, mpg%nCv,     mpg%fsCv,  'fsCv')
      call cfruin (iun, mpg%nFs*2,   mpg%fsFcP, 'fsFcP')
      call cfruin (iun, mpg%nFc,     mpg%fsFc,  'fsFc')
      call cfruin (iun, mpg%nFc,     mpg%fcReg,  'fcReg')
      call cfruin (iun, mpg%nCv,     mpg%cvReg,  'cvReg')
!      call cfruin (iun, mpg%nCv,     mpg%cvOnClosedSurface,
!     .                                  'cvOnClosedSurface')

      nnreg(0) = maxval(mpg%cvReg)
      nnreg(1) = maxval(mpg%fcReg)

      return
      end subroutine read_mapping

***********************************************************************

      subroutine write_mapping (nout, mpg)

      implicit none
      integer, intent(in) :: nout
      integer :: idum(0:9)
      type (mapping), intent (in) :: mpg

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
      idum(0) = mpg%nVmxCv
      call cfwuin (nout, 1, idum, 'nVmxCv')
      idum(0) = mpg%nVmxFc
      call cfwuin (nout, 1, idum, 'nVmxFc')

      call cfwuin (nout, mpg%nCv*2,   mpg%cvFcP, 'cvFcP')
      call cfwuin (nout, mpg%nCmxFc,  mpg%cvFc,  'cvFc')
      call cfwuin (nout, mpg%nFc*2,   mpg%fcCv,  'fcCv')
      call cfwuin (nout, mpg%nFc*2,   mpg%fcVx,  'fcVx')
      call cfwuin (nout, mpg%nCv*2,   mpg%cvVxP, 'cvVxP')
      call cfwuin (nout, mpg%nCmxVx,  mpg%cvVx,  'cvVx')
      call cfwuin (nout, mpg%nVx*2,   mpg%vxFcP, 'vxFcP')
      call cfwuin (nout, mpg%nVmxFc,  mpg%vxFc,  'vxFc')
      call cfwuin (nout, mpg%nVx*2,   mpg%vxCvP, 'vxCvP')
      call cfwuin (nout, mpg%nVmxCv,  mpg%vxCv,  'vxCv')
      call cfwuin (nout, mpg%nFs*2,   mpg%fsCvP, 'fsCvP')
      call cfwuin (nout, mpg%nCv,     mpg%fsCv,  'fsCv')
      call cfwuin (nout, mpg%nFs*2,   mpg%fsFcP, 'fsFcP')
      call cfwuin (nout, mpg%nFc,     mpg%fsFc,  'fsFc')
      call cfwuin (nout, mpg%nFc,     mpg%fcReg,  'fcReg')
      call cfwuin (nout, mpg%nCv,     mpg%cvReg,  'cvReg')
!      call cfwuin (nout, mpg%nCv,     mpg%cvOnClosedSurface,
!     .                                   'cvOnClosedSurface')

      return
      end subroutine write_mapping

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
     .           inv, inv1, inv2, invmx, ifcc
      integer, allocatable :: indFc(:), indCv(:), fh(:), ch(:)
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
      m%nCi = nx*ny                         ! internal cells
      m%nCg = 2*nx + 2*ny                   ! guard cells
      m%nCv = m%nCi + m%nCg             ! total number of cells

      write (*,*) 'nCv, nCi, nCg', m%nCv, m%nCi, m%nCg

!  no. of cell faces
      m%nFc = ny*(nx+1) + nx*(ny+1)

! no. of vertices
      m%nVx = (nx+1) * (ny+1)

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
          ic = ix + (iy-1)*nx
          if1 = m%cvFcP(ic,1)
          if2 = if1 + 1
          if3 = if2 + 1
          if4 = if3 + 1
          m%cvFc(if1) = ix + (iy-1)*(nx+1)
!          m%cvFc(if1) = leftix(ix-1,iy-1)+2 + (iy-1)*(nx+1)
!          m%cvFc(if2) = ix+1 + (iy-1)*(nx+1)
          m%cvFc(if2) = rightix(ix-1,iy-1)+1 + rightiy(ix-1,iy-1)*(nx+1)
          m%cvFc(if3) = nfx + iy + (ix-1)*(ny+1)
!          m%cvFc(if4) = nfx + iy+1 + (ix-1)*(ny+1)
          m%cvFc(if4) = nfx + topiy(ix-1,iy-1)+1+topix(ix-1,iy-1)*(ny+1)
          m%intcellP(if1) = 1.0_R8
          m%intcellP(if2) = 1.0_R8
          m%intcellP(if3) = 0.0_R8
          m%intcellP(if4) = 0.0_R8
          m%intcellR(if1) = 0.0_R8
          m%intcellR(if2) = 0.0_R8
          m%intcellR(if3) = 1.0_R8
          m%intcellR(if4) = 1.0_R8
          m%cvVx(if1) = ix + (iy-1)*(nx+1)
!          m%cvVx(if1) = leftix(ix-1,iy-1)+2 + (iy-1)*(nx+1)
!          m%cvVx(if2) = ix+1 + (iy-1)*(nx+1)
          m%cvVx(if2) = rightix(ix-1,iy-1)+1 + rightiy(ix-1,iy-1)*(nx+1)
          m%cvVx(if3) = topix(ix-1,iy-1)+1 + topiy(ix-1,iy-1)*(nx+1)
!          m%cvVx(if3) = leftix(ix-1,iy-1)+2 + iy*(nx+1)
!          m%cvVx(if4) = ix+1 + iy*(nx+1)
          m%cvVx(if4) = rightix(topix(ix-1,iy-1),topiy(ix-1,iy-1))+1 + 
     .                 rightiy(topix(ix-1,iy-1),topiy(ix-1,iy-1))*(nx+1)
        end do
      end do

! list of cell faces and vertices for guard cells
      do ix = 1, nx
        ic = m%nCi + ix
        m%cvFc(m%cvFcP(ic,1)) = nfx + 
     .               topiy(ix-1,-1)+1 + topix(ix-1,-1)*(ny+1)
        m%intcellP(m%cvFcP(ic,1)) = 1.0_R8
        m%intcellR(m%cvFcP(ic,1)) = 1.0_R8
        m%cvVx(m%cvVxP(ic,1)) = 
     .               topix(ix-1,-1)+1 + topiy(ix-1,-1)*(nx+1)
        m%cvVx(m%cvVxP(ic,1)+1) = 
     .               rightix(topix(ix-1,-1),topiy(ix-1,-1))+1 + 
     .               rightiy(topix(ix-1,-1),topiy(ix-1,-1))*(nx+1) 
        ic = m%nCi + nx + ix
        m%cvFc(m%cvFcP(ic,1)) = nfx + ny+1 + (ix-1)*(ny+1)
        m%intcellP(m%cvFcP(ic,1)) = 1.0_R8
        m%intcellR(m%cvFcP(ic,1)) = 1.0_R8
        m%cvVx(m%cvVxP(ic,1)) = ix + ny*(nx+1)
        m%cvVx(m%cvVxP(ic,1)+1) = 
     .               rightix(ix-1,ny)+1 + rightiy(ix-1,ny)*(nx+1)
      end do

      do iy = 1, ny
        ic = m%nCi + 2*nx + iy
        m%cvFc(m%cvFcP(ic,1)) = 
     .               rightix(-1,iy-1)+1 + rightiy(-1,iy-1)*(nx+1)
        m%intcellP(m%cvFcP(ic,1)) = 1.0_R8
        m%intcellR(m%cvFcP(ic,1)) = 1.0_R8
        m%cvVx(m%cvVxP(ic,1)) = 
     .               rightix(-1,iy-1)+1 + rightiy(-1,iy-1)*(nx+1)
        m%cvVx(m%cvVxP(ic,1)+1) = 
     .               rightix(topix(-1,iy-1),topiy(-1,iy-1))+1 + 
     .               rightiy(topix(-1,iy-1),topiy(-1,iy-1))*(nx+1)
        ic = m%nCi + 2*nx + ny + iy
        m%cvFc(m%cvFcP(ic,1)) = nx+1 + (iy-1)*(nx+1)
        m%intcellP(m%cvFcP(ic,1)) = 1.0_R8
        m%intcellR(m%cvFcP(ic,1)) = 1.0_R8
        m%cvVx(m%cvVxP(ic,1)) = nx+1 + (iy-1)*(nx+1) 
        m%cvVx(m%cvVxP(ic,1)+1) = 
     .               topix(nx,iy-1)+1 + topiy(nx,iy-1)*(nx+1)
      end do

! listing for every face the two corresponding cells
! and the corresponding vertices
! x faces
       do iy = 1, ny
         do ix = 2, nx
           ic = ix + (iy-1)*nx
           ifc = ix + (iy-1)*(nx+1)
           m%fcCv(ifc,1) = leftix(ix-1,iy-1)+1 + (iy-1)*nx
           m%fcCv(ifc,2) = ic
           m%fcVx(ifc,1) = ix + (iy-1)*(nx+1)
           m%fcVx(ifc,2) = ix + iy*(nx+1)
         end do

! left guard cells
         ic = 1 + (iy-1)*nx
         ifc = 1 + (iy-1)*(nx+1)
         m%fcCv(ifc,1) = m%nCi + 2*nx + iy
         m%fcCv(ifc,2) = ic
         m%fcVx(ifc,1) = 1 + (iy-1)*(nx+1)
         m%fcVx(ifc,2) = 1 + iy*(nx+1)

! right guard cells
         ic = nx+1 + (iy-1)*nx
         ifc = nx+1 + (iy-1)*(nx+1)
         m%fcCv(ifc,1) = ic - 1
         m%fcCv(ifc,2) = m%nCi + 2*nx + ny + iy
         m%fcVx(ifc,1) = nx+1 + (iy-1)*(nx+1)
         m%fcVx(ifc,2) = nx+1 + iy*(nx+1)
       end do

! y faces
       do ix = 1, nx
         do iy = 2, ny
           ic = ix + (iy-1)*nx
           ifc = nfx + iy + (ix-1)*(ny+1)
           m%fcCv(ifc,2) = ic
           m%fcCv(ifc,1) = ic - nx
           m%fcVx(ifc,2) = ix + (iy-1)*(nx+1)
           m%fcVx(ifc,1) = rightix(ix-1,iy-1)+1 + (iy-1)*(nx+1)
         end do

! lower guard cells
         ic = ix
         ifc = nfx + 1 + (ix-1)*(ny+1)
         m%fcCv(ifc,2) = ic
         m%fcCv(ifc,1) = m%nCi + ix
         m%fcVx(ifc,2) = ix 
         m%fcVx(ifc,1) = rightix(ix-1,0)+1 

! upper guard cells
         ic = ix + (ny-1)*nx
         ifc = nfx + ix*(ny+1)
         m%fcCv(ifc,2) = m%nCi + nx + ix
         m%fcCv(ifc,1) = ic
         m%fcVx(ifc,2) = ix + ny*(nx+1)
         m%fcVx(ifc,1) = rightix(ix-1,ny)+1 + ny*(nx+1)
       end do

! merging x-point vertex
       do icut = 1, nncut
         ix1 = leftcut(icut)
         iy1 = topcut(icut)
         ix2 = rightcut(icut)
         iy2 = topcut(icut)
         ivx1 = ix1+1 + iy1*(nx+1)
         ivx2 = ix2+1 + iy2*(nx+1)

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
      do iy = 0, ny-1
        do ix = 0, nx-1
          ic = (ix+1) + iy*nx
          m%cvReg(ic) = region(ix,iy,0)
          m%cvOnClosedSurface(ic) = on_closed_surface(ix,iy)
	end do
      end do	

      do ix = 0, nx-1
        ic = m%nCi + ix+1
        m%cvReg(ic) = region(ix,-1,0)
        m%cvOnClosedSurface(ic) = on_closed_surface(ix,-1)

        ic = m%nCi + nx + ix+1
        m%cvReg(ic) = region(ix,ny,0)
        m%cvOnClosedSurface(ic) = on_closed_surface(ix,ny)
      end do
 
      do iy = 0, ny-1
        ic = m%nCi + 2*nx + iy+1
        m%cvReg(ic) = region(-1,iy,0)
        m%cvOnClosedSurface(ic) = on_closed_surface(-1,iy)

        ic = m%nCi + 2*nx + ny + iy+1
        m%cvReg(ic) = region(nx,iy,0)
        m%cvOnClosedSurface(ic) = on_closed_surface(nx,iy)
      end do

! regions (faces)
 ! x faces
       do iy = 0, ny-1
         do ix = 1, nx-1
           ifc = ix+1 + iy*(nx+1)
           m%fcReg(ifc) = region(ix,iy,1)
         end do

! left boundary face (to guard cells)
         ix = 0
         ifc = ix+1 + iy*(nx+1)
         m%fcReg(ifc) = region(ix,iy,1)

! right guard cells
         ix = nx
         ifc = ix+1 + iy*(nx+1)
         m%fcReg(ifc) = region(ix,iy,1)
        end do

! y faces
       do ix = 0, nx-1
         do iy = 1, ny-1
           ifc = nfx + iy+1 + ix*(ny+1)
           m%fcReg(ifc) = region(ix,iy,2)
         end do

! lower guard cells
         iy = 0
         ifc = nfx + 1 + ix*(ny+1)
         m%fcReg(ifc) = region(ix,iy,2)

! upper guard cells
         iy = ny
         ifc = nfx + (ix+1)*(ny+1)
         m%fcReg(ifc) = region(ix,iy,2)
      end do

      mxReg = maxval(m%cvReg(1:m%nCv))

      ifs = 0
      icv = 0
      ifcc = 0
      do ireg= 1, mxReg
	do iy = 0, ny-1
	  istart = -1
	  do ix = 0, nx-1
            ic = (ix+1) + iy*nx
            ifc = ix+1 + iy*(nx+1)
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
          end do
        end do	
      end do
      m%nFs = ifs

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

      return
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
        elseif (bcchar(ib).eq.'E' .or. bcchar(ib).eq.'W') then
          if (bottomiy(bcpos(ib),bcstart(ib)) .eq. -2)
     >      bcs(ib) = bcs(ib) + 1
          if (topiy   (bcpos(ib),bcend(ib)  ) .eq. ny+1)
     >      bce(ib) = bce(ib) - 1
        else
          call xerrab ('bchar not yet available for WG')
        endif
        bcl(ib) = bce(ib) - bcs(ib) + 1
        m%mxnBc = m%mxnBc + bcl(ib)
      enddo
      if (m%mxnBc .ne. m%nCg ) 
     >  call xerrab ('calc_mapping_bc -- wrong mxnBc')


      call alloc_mapping_bc (m)

      inbc = 0 ! index in bcCv
      do ib = 1, nbc
        m%bcCvP(ib,1) = inbc + 1
        m%bcCvP(ib,2) = bcl(ib)
        if (bcchar(ib).eq.'S') then
          do ibc = 1, bcl(ib)
            inbc = inbc + 1
            m%bcCv(inbc,1) = m%nCi + bcs(ib) + ibc ! number of the guard cell
            m%bcCv(inbc,2) = bcs(ib) + ibc           ! number of the corresponding domain cell
            m%bcFcOr(inbc) = -1.0_R8
          enddo
        elseif (bcchar(ib).eq.'N') then
          do ibc = 1, bcl(ib)
            inbc = inbc + 1
            m%bcCv(inbc,1) = m%nCi + nx + bcs(ib) + ibc
            m%bcCv(inbc,2) = nx*(ny-1) + bcs(ib) + ibc 
            m%bcFcOr(inbc) = 1.0_R8
          enddo
        elseif (bcchar(ib).eq.'W') then
          do ibc = 1, bcl(ib)
            inbc = inbc + 1
            m%bcCv(inbc,1) = m%nCi + 2*nx + bcs(ib) + ibc
            m%bcCv(inbc,2) = nx*(bcs(ib) + ibc - 1) + 1
            m%bcFcOr(inbc) = -1.0_R8
          enddo
        elseif (bcchar(ib).eq.'E') then 
          do ibc = 1, bcl(ib)
            inbc = inbc + 1
            m%bcCv(inbc,1) = m%nCi + 2*nx + ny + bcs(ib) + ibc
            m%bcCv(inbc,2) = nx*(bcs(ib) +ibc)
            m%bcFcOr(inbc) = 1.0_R8
          enddo
        else
          call xerrab ('calc_mapping_bc -- wrong BCCHAR')
        endif
      enddo     

      if (inbc .ne. m%nCg ) 
     >  call xerrab ('calc_mapping_bc -- problem setting up bc arrays')

      return
      end subroutine calc_mapping_bc


***********************************************************************



      subroutine calc_mapping_rc (nx, ny, ns, m)
      use b2mod_neutrals_namelist
      implicit none
      integer :: nx, ny, ns, istra, istrac, instra
      integer :: rcs(nstrat), rce(nstrat), rcl(nstrat)
      type (mapping), intent (inout) :: m

      m%nRc   = nstrat  ! number of recycling boundaries
      m%mxnRc = 0       ! number of cells with imposed recycling
      do istra = 1, nstrat
        rcs(istra) = rcstart(istra)
        rce(istra) = rcend(istra)
        ! exclude corner cells if present
        if (crcstra(istra).eq.'N' .or. crcstra(istra).eq.'S') then
          if (leftix (rcstart(istra),rcpos(istra)) .eq. -2)  
     >      rcs(istra) = rcs(istra) + 1
          if (rightix(rcend(istra),  rcpos(istra)) .eq. nx+1) 
     >      rce(istra) = rce(istra) - 1
          rcl(istra) = rce(istra) - rcs(istra) + 1
          m%mxnRc = m%mxnRc + rcl(istra)
        elseif (crcstra(istra).eq.'E' .or. crcstra(istra).eq.'W') then
          if (bottomiy(rcpos(istra),rcstart(istra)) .eq. -2)
     >      rcs(istra) = rcs(istra) + 1
          if (topiy   (rcpos(istra),rcend(istra)  ) .eq. ny+1)
     >      rce(istra) = rce(istra) - 1
          rcl(istra) = rce(istra) - rcs(istra) + 1
          m%mxnRc = m%mxnRc + rcl(istra)
        elseif (crcstra(istra).eq.'V' .or. crcstra(istra).eq.'T') then
          rcl(istra) = 0
        else
          call xerrab ('Wrong CRCSTRA or not yet available for WG')
        endif
      enddo

      call alloc_mapping_rc (m)

      instra = 0 ! index in rcCv
      do istra = 1, nstrat
        m%rcCvP(istra,1) = instra + 1
        m%rcCvP(istra,2) = rcl(istra)
        if (crcstra(istra).eq.'S') then
          do istrac = 1, rcl(istra)
            instra = instra + 1
            m%rcCv(instra,1) = m%nCi + rcs(istra) + istrac ! number of the guard cell
            m%rcCv(instra,2) = rcs(istra) + istrac           ! number of the corresponding domain cell
            m%rcFcOr(instra) = -1.0_R8
          enddo
        elseif (crcstra(istra).eq.'N') then
          do istrac = 1, rcl(istra)
            instra = instra + 1
            m%rcCv(instra,1) = m%nCi + nx + rcs(istra) + istrac
            m%rcCv(instra,2) = nx*(ny-1) + rcs(istra) + istrac
            m%rcFcOr(instra) = 1.0_R8
          enddo
        elseif (crcstra(istra).eq.'W') then
          do istrac = 1, rcl(istra)
            instra = instra + 1
            m%rcCv(instra,1) = m%nCi + 2*nx + rcs(istra) + istrac
            m%rcCv(instra,2) = nx*(rcs(istra) + istrac - 1) + 1
            m%rcFcOr(instra) = -1.0_R8
          enddo
        elseif (crcstra(istra).eq.'E') then 
          do istrac = 1, rcl(istra)
            instra = instra + 1
            m%rcCv(instra,1) = m%nCi + 2*nx + ny + rcs(istra) + 
     .        istrac
            m%rcCv(instra,2) = nx*(rcs(istra) + istrac)
            m%rcFcOr(instra) = 1.0_R8
          enddo
        endif
      enddo 
    
      return
      end subroutine calc_mapping_rc


***********************************************************************


      subroutine read_geometry (iun, mpg, gm)

      implicit none
      integer, intent(in) :: iun
      integer :: idum(0:9)
      type (mapping) , intent (in) :: mpg
      type (geometry), intent (out) :: gm

      call alloc_geometry (gm, mpg%nCv, mpg%nFc, mpg%nVx, mpg%nFs,
     .    mpg%nVmxCv)

! cell volumes 
      call cfrure (iun, mpg%nCv*4, gm%cvBb,   'cvBb')
      call cfrure (iun, mpg%nCv,   gm%cvX,    'cvX')
      call cfrure (iun, mpg%nCv,   gm%cvY,    'cvY')
      call cfrure (iun, mpg%nCv,   gm%cvS,    'cvS')
      call cfrure (iun, mpg%nCv,   gm%cvHz,   'cvHz')
      call cfrure (iun, mpg%nCv*2, gm%cvQgam, 'cvQgam')
      call cfrure (iun, mpg%nCv,   gm%cvVol,  'cvVol')
!      call cfrure (iun, mpg%nCv,   gm%cvOnedBsq,'cvOnedBsq')
!      call cfrure (iun, mpg%nCv,   gm%cvBzb,  'cvBzb')

! face quantities
!      call cfrure (iun, mpg%nFc*4, gm%fcBb,   'fcBb')
      call cfrure (iun, mpg%nFc,   gm%fcS,    'fcS')
      call cfrure (iun, mpg%nFc*2, gm%fcHc,   'fcHc')
      call cfrure (iun, mpg%nFc,   gm%fcHt,   'fcHt')
      call cfrure (iun, mpg%nFc*2, gm%fcVol,  'fcVol')
      call cfrure (iun, mpg%nFc*2, gm%fcQgam, 'fcQgam')
      call cfrure (iun, mpg%nFc*2, gm%fcQalf, 'fcQalf')
      call cfrure (iun, mpg%nFc*2, gm%fcQbet, 'fcQbet')
      call cfrure (iun, mpg%nFc,   gm%fcPbs,  'fcPbs')
!      call cfrure (iun, mpg%nFc,   gm%fcPbshz,'fcPbshz')

! vertice quantities
      call cfrure (iun, mpg%nVx,   gm%vxX,    'vxX')
      call cfrure (iun, mpg%nVx,   gm%vxY,    'vxY')
      call cfrure (iun, mpg%nVmxCv,gm%vxVol,  'vxVol')
      call cfrure (iun, mpg%nVx,   gm%vxFfbz, 'vxFfbz')
      call cfrure (iun, mpg%nVx,   gm%vxFpsi, 'vxFpsi')

      call cfrure (iun, mpg%nFs,   gm%fsConn, 'fsConn')

      return
      end subroutine read_geometry

***********************************************************************

      subroutine write_geometry (nout, mpg, gm)

      implicit none
      integer, intent(in) :: nout
      integer :: idum(0:9)
      type (mapping) , intent (in) :: mpg
      type (geometry), intent (in) :: gm

! cell volumes 
      call cfwure (nout, mpg%nCv*4, gm%cvBb,   'cvBb')
      call cfwure (nout, mpg%nCv,   gm%cvX,    'cvX')
      call cfwure (nout, mpg%nCv,   gm%cvY,    'cvY')
      call cfwure (nout, mpg%nCv,   gm%cvS,    'cvS')
      call cfwure (nout, mpg%nCv,   gm%cvHz,   'cvHz')
      call cfwure (nout, mpg%nCv*2, gm%cvQgam, 'cvQgam')
      call cfwure (nout, mpg%nCv,   gm%cvVol,  'cvVol')
!      call cfwure (nout, mpg%nCv,   gm%cvOnedBsq,'cvOnedBsq')
!      call cfwure (nout, mpg%nCv,   gm%cvBzb,  'cvBzb')

! face quantities
!      call cfwure (nout, mpg%nFc*4, gm%fcBb,   'fcBb')
      call cfwure (nout, mpg%nFc,   gm%fcS,    'fcS')
      call cfwure (nout, mpg%nFc*2, gm%fcHc,   'fcHc')
      call cfwure (nout, mpg%nFc,   gm%fcHt,   'fcHt')
      call cfwure (nout, mpg%nFc*2, gm%fcVol,  'fcVol')
      call cfwure (nout, mpg%nFc*2, gm%fcQgam, 'fcQgam')
      call cfwure (nout, mpg%nFc*2, gm%fcQalf, 'fcQalf')
      call cfwure (nout, mpg%nFc*2, gm%fcQbet, 'fcQbet')
      call cfwure (nout, mpg%nFc,   gm%fcPbs,  'fcPbs')
!      call cfwure (nout, mpg%nFc,   gm%fcPbshz,'fcPbshz')

! vertice quantities
      call cfwure (nout, mpg%nVx,   gm%vxX,    'vxX')
      call cfwure (nout, mpg%nVx,   gm%vxY,    'vxY')
      call cfwure (nout, mpg%nVmxCv,gm%vxVol,  'vxVol')
      call cfwure (nout, mpg%nVx,   gm%vxFfbz, 'vxFfbz')
      call cfwure (nout, mpg%nVx,   gm%vxFpsi, 'vxFpsi')

      call cfwure (nout, mpg%nFs,   gm%fsConn, 'fsConn')

      return
      end subroutine write_geometry

***********************************************************************

      subroutine prep_geometry(nx, ny, nCv, nCi, nFc, nVx, nFs, nVmxCv,
     .  gm, mpg)

      implicit none
      integer, intent(in) :: nx, ny, nCv, nCi, nFc, nVx, nFs, nVmxCv
      integer :: ix, iy, ic, ifc, nfx, ivx, i, j
      real (kind=R8) :: hy1(-1:nx,-1:ny)
      type (mapping) , intent (in) :: mpg
      type (geometry), intent (out) :: gm

      hy1 = hy*qz(:,:,1)

      nfx = (nx+1)*ny

      call alloc_geometry (gm, nCv, nFc, nVx, nFs, nVmxCv)

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
          gm%cvHz(ic) = hz(ix,iy)
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
        gm%cvHz(ic) = hz(ix,-1)
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
        gm%cvHz(ic) = hz(ix,ny)
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
        gm%cvHz(ic) = hz(-1,iy)
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
        gm%cvHz(ic) = hz(nx,iy)
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
           ifc = nfx + iy+1 + ix*(ny+1)
           gm%fcBb(ifc,0:3) = wbbv(ix,iy,0:3)
           gm%fcS(ifc) = gs(ix,iy,1)
           gm%fcHc(ifc,1) = 0.5_R8*hy(bottomix(ix,iy),bottomiy(ix,iy))       !cwdk: change orientation ??
           gm%fcHc(ifc,2) = 0.5_R8*hy(ix,iy)
           gm%fcHt(ifc)   = 
     .       0.5_R8*(hx(bottomix(ix,iy),bottomiy(ix,iy)) + hx(ix,iy))       !cwdk: computation could be made exact
           gm%fcVol(ifc,1) = 0.5_R8*vol(bottomix(ix,iy),bottomiy(ix,iy))
           gm%fcVol(ifc,2) = 0.5_R8*vol(ix,iy)
           gm%fcQgam(ifc,0) = qc(ix,iy,1)
           gm%fcQgam(ifc,1) = qs(ix,iy,1)
           gm%fcQalf(ifc,0) = 0.0_R8
           gm%fcQalf(ifc,1) = 1.0_R8
           gm%fcQbet(ifc,0) = gm%fcQgam(ifc,1)
           gm%fcQbet(ifc,1) = -gm%fcQgam(ifc,0)
           gm%fcPbs(ifc) = pbs(ix,iy,1)
           gm%fcPbshz(ifc) = pbshz(ix,iy,1)
         end do

! lower guard cells
         iy = 0
         ifc = nfx + 1 + ix*(ny+1)
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
         gm%fcQalf(ifc,1) = 1.0_R8
         gm%fcQbet(ifc,0) = gm%fcQgam(ifc,1)
         gm%fcQbet(ifc,1) = -gm%fcQgam(ifc,0)
         gm%fcPbs(ifc) = pbs(ix,iy,1)
         gm%fcPbshz(ifc) = pbshz(ix,iy,1)

! upper guard cells
         iy = ny
         ifc = nfx + (ix+1)*(ny+1)
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
           gm%vxX(ivx) = crx(ix,iy,0)
           gm%vxY(ivx) = cry(ix,iy,0)
           gm%vxFfbz(ivx) = ffbz(ix,iy,1)
           gm%vxFpsi(ivx) = fpsi(ix,iy,1)
           do i = 1, mpg%vxCvP(ivx,2)
             ic = mpg%vxCv(mpg%vxCvP(ivx,1) + i - 1)
             gm%vxVol(mpg%vxCvP(ivx,1) + i - 1) = 0.25*gm%cvVol(ic)
           end do
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

      end subroutine prep_geometry

**********************************************************************
      
      subroutine prep_plasma (nx, ny, ns, nCv, nCi, nFc, st)
      
      implicit none
      integer, intent(in) :: nx, ny, ns, nCi, nCv, nFc
      type(B2State), intent(inout) :: st
      integer :: ix, iy, ic, ifc, nfx

      call createB2State (nCv, nFc, ns, st)

      nfx = (nx+1)*ny

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
      
! guard cells iy = -1 and iy = ny
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
        do ix = 0, nx
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
     .              vadia(ix,iy,0,0:1,0:ns-1)
          st % dv % wadia(ifc,0:1,0:ns-1) = 
     .              wadia(ix,iy,0,0:1,0:ns-1)
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

        end do
      end do



! y faces

      do ix = 0, nx-1
        do iy = 0, ny
          ifc = nfx + iy+1 + ix*(ny+1)

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
      
      subroutine to_struct_plasma (nx, ny, ns, nCv, nCi, nFc, st)
      
      implicit none
      integer, intent(in) :: nx, ny, ns, nCi, nCv, nFc
      type(B2State), intent(in) :: st
      integer :: ix, iy, ic, ifc, nfx

      if (.not.allocated(fch)) call alloc_b2mod_plasma(nx,ny,ns)
      if (.not.allocated(ceqp)) call alloc_b2mod_transport(nx,ny,ns)
      if (.not.allocated(dna0)) 
     .    call alloc_b2mod_anomalous_transport(nx,ny,ns)
      if (.not.allocated(pcca)) call alloc_b2mod_work(nx,ny,ns)
      if (.not.allocated(pcca)) call alloc_b2mod_residuals(nx,ny,ns)

!      ns = st % ns

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
!          ue(ix,iy)     = st % dv % ue(ic)
      
! residuals
          resmo(ix,iy,0:ns-1) = st % dv % resmo(ic,0:ns-1)
          resco(ix,iy,0:ns-1) = st % dv % resco(ic,0:ns-1)
          respo(ix,iy)        = st % dv % respo(ic)
          reshe(ix,iy)        = st % dv % reshe(ic)
          reshi(ix,iy)        = st % dv % reshi(ic)
!          resht(ix,iy)        = st % dv % resht(ic)
          resmt(ix,iy)        = st % dv % resmt(ic)
      
! corrections
          corua(ix,iy,0:ns-1) = st % dv % corua(ic,0:ns-1)
          corpa(ix,iy,0:ns-1) = st % dv % corpa(ic,0:ns-1)
          corut(ix,iy)        = st % dv % corut(ic)
          corpo(ix,iy)        = st % dv % corpo(ic)
!          cortt(ix,iy)        = st % dv % cortt(ic)
          corte(ix,iy)        = st % dv % corte(ic)
          corti(ix,iy)        = st % dv % corti(ic)
      
! coupling coefficients
          pcca(ix,iy,0,0,0:ns-1) = st % dv % pcca(ic,0,0:ns-1)
          pccm(ix,iy,0,0)        = st % dv % pccm(ic,0)
          pcca(ix,iy,1,1,0:ns-1) = st % dv % pcca(ic,1,0:ns-1)
          pccm(ix,iy,1,1)        = st % dv % pccm(ic,1)

! B2SourceWork
!         sch0(ix,iy,0:3)        = st % srw % sch0(ic,0:3)
!         she0(ix,iy,0:3)        = st % srw % she0(ic,0:3)
!         shi0(ix,iy,0:3)        = st % srw % shi0(ic,0:3)
!         smo0(ix,iy,0:3,0:ns-1) = st % srw % smo0(ic,0:3,0:ns-1)
!         smq0(ix,iy,0:3,0:ns-1) = st % srw % smq0(ic,0:3,0:ns-1)
!         sna0(ix,iy,0:3,0:ns-1) = st % srw % sna0(ic,0:3,0:ns-1)

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
!        ue(ix,iy)     = st % dv % ue(ic)
        
! residuals
        resmo(ix,iy,0:ns-1) = st % dv % resmo(ic,0:ns-1)
        resco(ix,iy,0:ns-1) = st % dv % resco(ic,0:ns-1)
        respo(ix,iy)        = st % dv % respo(ic)
        reshe(ix,iy)        = st % dv % reshe(ic)
        reshi(ix,iy)        = st % dv % reshi(ic)
!        resht(ix,iy)        = st % dv % resht(ic)
        resmt(ix,iy)        = st % dv % resmt(ic)
      
! corrections
        corua(ix,iy,0:ns-1) = st % dv % corua(ic,0:ns-1)
        corpa(ix,iy,0:ns-1) = st % dv % corpa(ic,0:ns-1)
        corut(ix,iy)        = st % dv % corut(ic)
        corpo(ix,iy)        = st % dv % corpo(ic)
!        cortt(ix,iy)        = st % dv % cortt(ic)
        corte(ix,iy)        = st % dv % corte(ic)
        corti(ix,iy)        = st % dv % corti(ic)
      
! coupling coefficients
        pcca(ix,iy,0,0,0:ns-1) = st % dv % pcca(ic,0,0:ns-1)
        pccm(ix,iy,0,0)        = st % dv % pccm(ic,0)
        pcca(ix,iy,1,1,0:ns-1) = st % dv % pcca(ic,1,0:ns-1)
        pccm(ix,iy,1,1)        = st % dv % pccm(ic,1)
        
! B2SourceWork
!       sch0(ix,iy,0:3) = st % srw % sch0(ic,0:3)
!       she0(ix,iy,0:3) = st % srw % she0(ic,0:3)
!       shi0(ix,iy,0:3) = st % srw % shi0(ic,0:3)
!       smo0(ix,iy,0:3,0:ns-1) = st % srw % smo0(ic,0:3,0:ns-1)
!       smq0(ix,iy,0:3,0:ns-1) = st % srw % smq0(ic,0:3,0:ns-1)
!       sna0(ix,iy,0:3,0:ns-1) = st % srw % sna0(ic,0:3,0:ns-1)
        
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
!        ue(ix,iy)     = st % dv % ue(ic)
        
! residuals
        resmo(ix,iy,0:ns-1) = st % dv % resmo(ic,0:ns-1)
        resco(ix,iy,0:ns-1) = st % dv % resco(ic,0:ns-1)
        respo(ix,iy)        = st % dv % respo(ic)
        reshe(ix,iy)        = st % dv % reshe(ic)
        reshi(ix,iy)        = st % dv % reshi(ic)
!        resht(ix,iy)        = st % dv % resht(ic)
        resmt(ix,iy)        = st % dv % resmt(ic)
      
! corrections
        corua(ix,iy,0:ns-1) = st % dv % corua(ic,0:ns-1)
        corpa(ix,iy,0:ns-1) = st % dv % corpa(ic,0:ns-1)
        corut(ix,iy)        = st % dv % corut(ic)
        corpo(ix,iy)        = st % dv % corpo(ic)
!        cortt(ix,iy)        = st % dv % cortt(ic)
        corte(ix,iy)        = st % dv % corte(ic)
        corti(ix,iy)        = st % dv % corti(ic)
      
! coupling coefficients
        pcca(ix,iy,0,0,0:ns-1) = st % dv % pcca(ic,0,0:ns-1)
        pccm(ix,iy,0,0)        = st % dv % pccm(ic,0)
        pcca(ix,iy,1,1,0:ns-1) = st % dv % pcca(ic,1,0:ns-1)
        pccm(ix,iy,1,1)        = st % dv % pccm(ic,1)

! B2SourceWork
!       sch0(ix,iy,0:3)        = st % srw % sch0(ic,0:3)
!       she0(ix,iy,0:3)        = st % srw % she0(ic,0:3)
!       shi0(ix,iy,0:3)        = st % srw % shi0(ic,0:3)
!       smo0(ix,iy,0:3,0:ns-1) = st % srw % smo0(ic,0:3,0:ns-1)
!       smq0(ix,iy,0:3,0:ns-1) = st % srw % smq0(ic,0:3,0:ns-1)
!       sna0(ix,iy,0:3,0:ns-1) = st % srw % sna0(ic,0:3,0:ns-1)

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
!        ue(ix,iy)     = st % dv % ue(ic)
        
! residuals
        resmo(ix,iy,0:ns-1) = st % dv % resmo(ic,0:ns-1)
        resco(ix,iy,0:ns-1) = st % dv % resco(ic,0:ns-1)
        respo(ix,iy)        = st % dv % respo(ic)
        reshe(ix,iy)        = st % dv % reshe(ic)
        reshi(ix,iy)        = st % dv % reshi(ic)
!        resht(ix,iy)        = st % dv % resht(ic)
        resmt(ix,iy)        = st % dv % resmt(ic)
      
! corrections
        corua(ix,iy,0:ns-1) = st % dv % corua(ic,0:ns-1)
        corpa(ix,iy,0:ns-1) = st % dv % corpa(ic,0:ns-1)
        corut(ix,iy)        = st % dv % corut(ic)
        corpo(ix,iy)        = st % dv % corpo(ic)
!        cortt(ix,iy)        = st % dv % cortt(ic)
        corte(ix,iy)        = st % dv % corte(ic)
        corti(ix,iy)        = st % dv % corti(ic)
      
! coupling coefficients
        pcca(ix,iy,0,0,0:ns-1) = st % dv % pcca(ic,0,0:ns-1)
        pccm(ix,iy,0,0)        = st % dv % pccm(ic,0)
        pcca(ix,iy,1,1,0:ns-1) = st % dv % pcca(ic,1,0:ns-1)
        pccm(ix,iy,1,1)        = st % dv % pccm(ic,1)

! B2SourceWork
!       sch0(ix,iy,0:3)        = st % srw % sch0(ic,0:3)
!       she0(ix,iy,0:3)        = st % srw % she0(ic,0:3)
!       shi0(ix,iy,0:3)        = st % srw % shi0(ic,0:3)
!       smo0(ix,iy,0:3,0:ns-1) = st % srw % smo0(ic,0:3,0:ns-1)
!       smq0(ix,iy,0:3,0:ns-1) = st % srw % smq0(ic,0:3,0:ns-1)
!       sna0(ix,iy,0:3,0:ns-1) = st % srw % sna0(ic,0:3,0:ns-1)

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
!        ue(ix,iy)     = st % dv % ue(ic)
        
! residuals
        resmo(ix,iy,0:ns-1) = st % dv % resmo(ic,0:ns-1)
        resco(ix,iy,0:ns-1) = st % dv % resco(ic,0:ns-1)
        respo(ix,iy)        = st % dv % respo(ic)
        reshe(ix,iy)        = st % dv % reshe(ic)
        reshi(ix,iy)        = st % dv % reshi(ic)
!        resht(ix,iy)        = st % dv % resht(ic)
        resmt(ix,iy)        = st % dv % resmt(ic)
      
! corrections
        corua(ix,iy,0:ns-1) = st % dv % corua(ic,0:ns-1)
        corpa(ix,iy,0:ns-1) = st % dv % corpa(ic,0:ns-1)
        corut(ix,iy)        = st % dv % corut(ic)
        corpo(ix,iy)        = st % dv % corpo(ic)
!        cortt(ix,iy)        = st % dv % cortt(ic)
        corte(ix,iy)        = st % dv % corte(ic)
        corti(ix,iy)        = st % dv % corti(ic)
      
! coupling coefficients
        pcca(ix,iy,0,0,0:ns-1) = st % dv % pcca(ic,0,0:ns-1)
        pccm(ix,iy,0,0)        = st % dv % pccm(ic,0)
        pcca(ix,iy,1,1,0:ns-1) = st % dv % pcca(ic,1,0:ns-1)
        pccm(ix,iy,1,1)        = st % dv % pccm(ic,1)

! B2SourceWork
!       sch0(ix,iy,0:3)        = st % srw % sch0(ic,0:3)
!       she0(ix,iy,0:3)        = st % srw % she0(ic,0:3)
!       shi0(ix,iy,0:3)        = st % srw % shi0(ic,0:3)
!       smo0(ix,iy,0:3,0:ns-1) = st % srw % smo0(ic,0:3,0:ns-1)
!       smq0(ix,iy,0:3,0:ns-1) = st % srw % smq0(ic,0:3,0:ns-1)
!       sna0(ix,iy,0:3,0:ns-1) = st % srw % sna0(ic,0:3,0:ns-1)

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

! Corner cells: add-hoc interpolation
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


!B2PlasmaSnapshot
!          fch(ix,iy,0,0:1)        = st % psnl % fch(ifc,0:1)
!          fna(ix,iy,0,0:1,0:ns-1) = st % psnl % fna(ifc,0:1,0:ns-1)
!          fhe(ix,iy,0,0:1)        = st % psnl % fhe(ifc,0:1)
!          fhi(ix,iy,0,0:1)        = st % psnl % fhi(ifc,0:1)
!          fch(ix,iy,0,0:1)        = st % psnc % fch(ifc,0:1)
!          fna(ix,iy,0,0:1,0:ns-1) = st % psnc % fna(ifc,0:1,0:ns-1)
!          fhe(ix,iy,0,0:1)        = st % psnc % fhe(ifc,0:1)
!          fhi(ix,iy,0,0:1)        = st % psnc % fhi(ifc,0:1)

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
          
! B2PlasmaSnapshot
!          fch(ix,iy,1,0:1)        = st % psnl % fch(ifc,0:1)
!          fna(ix,iy,1,0:1,0:ns-1) = st % psnl % fna(ifc,0:1,0:ns-1)
!          fhe(ix,iy,1,0:1)        = st % psnl % fhe(ifc,0:1)
!          fhi(ix,iy,1,0:1)        = st % psnl % fhi(ifc,0:1)
!          fch(ix,iy,1,0:1)        = st % psnc % fch(ifc,0:1)
!          fna(ix,iy,1,0:1,0:ns-1) = st % psnc % fna(ifc,0:1,0:ns-1)
!          fhe(ix,iy,1,0:1)        = st % psnc % fhe(ifc,0:1)
!          fhi(ix,iy,1,0:1)        = st % psnc % fhi(ifc,0:1)

        end do
      end do     

      aresco(0:ns-1) = st % diag % aresco(0:ns-1)
      aresmo(0:ns-1) = st % diag % aresmo(0:ns-1)

!     acorpa(0:ns-1) = st % diag % acorpa(0:ns-1)
!     acorua(0:ns-1) = st % diag % acorua(0:ns-1)

!     rescoreg(0:nnreg(0),0:ns-1) =
!    .         st % diag % rescoreg(0:nnreg(0),0:ns-1)
!     resmoreg(0:nnreg(0),0:ns-1)
!    .         st % diag % resmoreg(0:nnreg(0),0:ns-1)
!     reshereg(0:nnreg(0)) = st % diag % reshereg(0:nnreg(0))
!     reshireg(0:nnreg(0)) = st % diag % reshireg(0:nnreg(0))

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


***********************************************************************


      subroutine write_plasma_state (fname, ns, nCv, nFc, st)
      
      implicit none
      character(*), intent(in) :: fname
      integer, intent(in) :: ns, nCv, nFc
      type(B2State), intent(inout) :: st
      integer ::  nput, io, idum(0:2) 
      real (kind=R8) :: time(1)
      character lble*120
      character*10, save ::
     *   b2fstate_version =  '03.002.000'
      real (kind=R8) ::
     *   zpdum(0:nsdecl-1,4)

*-----------------------------------------------------------------------
*.computation

* ..preliminaries
*   ..subprogram start-up calls
      call subini ('write_plasma_state')
*   ..test ns
      call xertst (1.le.ns, 'faulty argument ns')
      
* ..save state and fluxes
*   ..save plasma state

      open (newunit=nput,file=fname,iostat=io)

*   ..get label
      call b2xxid ('b2mn', lble(1:60))
      lble(61:120) = 
     *  '                                                            '

*   ..write header
      idum(0) = ns
      idum(1) = nCv
      idum(2) = nFc
      call cfverw (nput, b2fstate_version)
      call cfwuin (nput, 3, idum, 'ns,nCv,nFc')
      call cfwuch (nput, 120, lble, 'label')
      call b2wuzd (nput, b2fstate_version, ns, zamin, zamax, zn, am)

      call cfwure (nput, nCv*ns,   st % pl % na, 'na')
      call cfwure (nput, nCv,      st % dv % ne, 'ne')
      call cfwure (nput, nCv*ns,   st % pl % ua, 'ua')
      call cfwure (nput, nFc*2*ns, st % dv % uadia, 'uadia') 
      call cfwure (nput, nCv,      st % pl % te, 'te')
      call cfwure (nput, nCv,      st % pl % ti, 'ti')
      call cfwure (nput, nCv,      st % pl % po, 'po')

*   ..save fluxes
      call cfwure (nput, nFc*2*ns, st % dv % fna, 'fna')
      call cfwure (nput, nFc*2,    st % dv % fhe, 'fhe')
      call cfwure (nput, nFc*2,    st % dv % fhi, 'fhi')
      call cfwure (nput, nFc*2,    st % dv % fch, 'fch')
      call cfwure (nput, nFc*2,    st % dv % fch_32, 'fch_32')
      call cfwure (nput, nFc*2,    st % dv % fch_52, 'fch_52')
      call cfwure (nput, nCv*ns,   st % dv % kinrgy, 'kinrgy')

*   ..save time
      time(1)=tim
      call cfwure (nput, 1, time, 'time')

*   ..save parallel current
      call cfwure (nput, nFc*2,    st % dv % fch_p, 'fch_p')

*   ..saving various fluxes and velocities needed by advanced SOLPS-ITER features
      call cfwure (nput, nFc*2*ns, st % dv % fna_mdf, 'fna_mdf')
      call cfwure (nput, nFc*2,    st % dv % fhe_mdf, 'fhe_mdf')
      call cfwure (nput, nFc*2,    st % dv % fhi_mdf, 'fhi_mdf')
      call cfwure (nput, nFc*2*ns, st % dv % fna_fcor, 'fna_fcor')
      call cfwure (nput, nFc*2*ns, st % dv % fna_nodrift, 'fna_nodrift')
      call cfwure (nput, nFc*2*ns, st % dv % fna_he, 'fna_he')
      call cfwure (nput, nFc*2*ns, st % dv % fnaPSch, 'fnaPSch')
      call cfwure (nput, nFc*2,    st % dv % fhePSch, 'fhePSch')
      call cfwure (nput, nFc*2,    st % dv % fhiPSch, 'fhiPSch')
      call cfwure (nput, nFc*2*ns, st % dv % fna_eir, 'fna_eir')
      call cfwure (nput, nFc*2,    st % dv % fne_eir, 'fne_eir')
      call cfwure (nput, nFc*2,    st % dv % fhe_eir, 'fhe_eir')
      call cfwure (nput, nFc*2,    st % dv % fhi_eir, 'fhi_eir')
      call cfwure (nput, nFc*2*ns, st % dv % fna_32, 'fna_32')
      call cfwure (nput, nFc*2*ns, st % dv % fna_52, 'fna_52')
      call cfwure (nput, nFc*2,    st % dv % fni_32, 'fni_32')
      call cfwure (nput, nFc*2,    st % dv % fni_52, 'fni_52')
      call cfwure (nput, nFc*2,    st % dv % fne_32, 'fne_32')
      call cfwure (nput, nFc*2,    st % dv % fne_52, 'fne_52')
      call cfwure (nput, nFc*2,    st % dv % fchdia, 'fchdia')
      call cfwure (nput, nFc*2,    st % dv % fchin, 'fchin')
      call cfwure (nput, nFc*2,    st % dv % fchvispar, 'fchvispar')
      call cfwure (nput, nFc*2,    st % dv % fchvisper, 'fchvisper')
      call cfwure (nput, nFc*2,    st % dv % fchvisq, 'fchvisq')
      call cfwure (nput, nFc*2,    st % dv % fchinert, 'fchinert')
      call cfwure (nput, nFc*2*ns, st % dv % vaecrb, 'vaecrb')
      call cfwure (nput, nFc*2*ns, st % dv % vadia, 'vadia')
      call cfwure (nput, nFc*2*ns, st % dv % wadia, 'wadia')
      call cfwure (nput, nFc*2,    st % dv % veecrb, 'veecrb')
      call cfwure (nput, nFc*2,    st % dv % vedia, 'vedia')
      call cfwure (nput, nFc*2,    st % dv % floe_noc, 'floe_noc')
      call cfwure (nput, nFc*2,    st % dv % floi_noc, 'floi_noc')

* ..return
      call subend ()
      return

      end subroutine write_plasma_state


***********************************************************************


      subroutine read_plasma_state (nput, ns, nCv, nFc, st)
      
      implicit none
      integer, intent(in) :: nput, ns, nCv, nFc
      type(B2State), intent(inout) :: st
      integer :: ix, iy, ic, ifc, nfx, idum(3)
      real (kind=R8) :: time(1)
      character lbli*120
      real (kind=R8) ::
     *   zpdum(0:nsdecl-1,4)
      character*10, save :: b2fstate_version 

*-----------------------------------------------------------------------
*.computation

* ..preliminaries
*   ..subprogram start-up calls
      call subini ('read_plasma_state')
*   ..test nput
      call xertst (1.le.nput, 'faulty argument nput')
*   ..test ns
      call xertst (1.le.ns, 'faulty argument ns')
 
      call createB2State (nCv, nFc, ns, st)
     
* ..read state and fluxes
*   ..read plasma state
      call cfverr (nput, b2fstati_version)
      call cfruin (nput, 3, idum, 'ns,nCv,nFc')
      call cfruch (nput, 120, lbli, 'label')
      call b2ruzd (nput, b2fstati_version, ns, 
     &          zpdum(0,1), zpdum(0,2), zpdum(0,3), zpdum(0,4), .false.)

      call cfrure (nput, nCv*ns,   st % pl % na, 'na')
      call cfrure (nput, nCv,      st % dv % ne, 'ne')
      call cfrure (nput, nCv*ns,   st % pl % ua, 'ua')
      call cfrure (nput, nFc*2*ns, st % dv % uadia, 'uadia') 
      call cfrure (nput, nCv,      st % pl % te, 'te')
      call cfrure (nput, nCv,      st % pl % ti, 'ti')
      call cfrure (nput, nCv,      st % pl % po, 'po')

*   ..read fluxes
      call cfrure (nput, nFc*2*ns, st % dv % fna, 'fna')
      call cfrure (nput, nFc*2,    st % dv % fhe, 'fhe')
      call cfrure (nput, nFc*2,    st % dv % fhi, 'fhi')
      call cfrure (nput, nFc*2,    st % dv % fch, 'fch')
      call cfrure (nput, nFc*2,    st % dv % fch_32, 'fch_32')
      call cfrure (nput, nFc*2,    st % dv % fch_52, 'fch_52')
      call cfrure (nput, nCv*ns,   st % dv % kinrgy, 'kinrgy')

*   ..read time
      call cfrure (nput, 1, time, 'time')
      tim = time(1)

*   ..read parallel current
      call cfrure (nput, nFc*2,    st % dv % fch_p, 'fch_p')

*   ..read various fluxes and velocities needed by advanced SOLPS-ITER features
      call cfrure (nput, nFc*2*ns, st % dv % fna_mdf, 'fna_mdf')
      call cfrure (nput, nFc*2,    st % dv % fhe_mdf, 'fhe_mdf')
      call cfrure (nput, nFc*2,    st % dv % fhi_mdf, 'fhi_mdf')
      call cfrure (nput, nFc*2*ns, st % dv % fna_fcor, 'fna_fcor')
      call cfrure (nput, nFc*2*ns, st % dv % fna_nodrift, 'fna_nodrift')
      call cfrure (nput, nFc*2*ns, st % dv % fna_he, 'fna_he')
      call cfrure (nput, nFc*2*ns, st % dv % fnaPSch, 'fnaPSch')
      call cfrure (nput, nFc*2,    st % dv % fhePSch, 'fhePSch')
      call cfrure (nput, nFc*2,    st % dv % fhiPSch, 'fhiPSch')
      call cfrure (nput, nFc*2*ns, st % dv % fna_eir, 'fna_eir')
      call cfrure (nput, nFc*2,    st % dv % fne_eir, 'fne_eir')
      call cfrure (nput, nFc*2,    st % dv % fhe_eir, 'fhe_eir')
      call cfrure (nput, nFc*2,    st % dv % fhi_eir, 'fhi_eir')
      call cfrure (nput, nFc*2*ns, st % dv % fna_32, 'fna_32')
      call cfrure (nput, nFc*2*ns, st % dv % fna_52, 'fna_52')
      call cfrure (nput, nFc*2,    st % dv % fni_32, 'fni_32')
      call cfrure (nput, nFc*2,    st % dv % fni_52, 'fni_52')
      call cfrure (nput, nFc*2,    st % dv % fne_32, 'fne_32')
      call cfrure (nput, nFc*2,    st % dv % fne_52, 'fne_52')
      call cfrure (nput, nFc*2,    st % dv % fchdia, 'fchdia')
      call cfrure (nput, nFc*2,    st % dv % fchin, 'fchin')
      call cfrure (nput, nFc*2,    st % dv % fchvispar, 'fchvispar')
      call cfrure (nput, nFc*2,    st % dv % fchvisper, 'fchvisper')
      call cfrure (nput, nFc*2,    st % dv % fchvisq, 'fchvisq')
      call cfrure (nput, nFc*2,    st % dv % fchinert, 'fchinert')
      call cfrure (nput, nFc*2*ns, st % dv % vaecrb, 'vaecrb')
      call cfrure (nput, nFc*2*ns, st % dv % vadia, 'vadia')
      call cfrure (nput, nFc*2*ns, st % dv % wadia, 'wadia')
      call cfrure (nput, nFc*2,    st % dv % veecrb, 'veecrb')
      call cfrure (nput, nFc*2,    st % dv % vedia, 'vedia')
      call cfrure (nput, nFc*2,    st % dv % floe_noc, 'floe_noc')
      call cfrure (nput, nFc*2,    st % dv % floi_noc, 'floi_noc')

* ..return
      call subend ()
      return

      end subroutine read_plasma_state

***********************************************************************

      subroutine conv_bcparms (nx, ny, ns, m)
      use b2mod_boundary_namelist
      implicit none
      integer, intent(in) :: nx, ny, ns
      type (mapping), intent (inout) :: m
      integer :: ib, ibc, is, inbc, nfx
      integer :: bcs(nbc), bce(nbc), bcl(nbc)

!      if (any(bc_list_size > 0))
!     >   call xerrab ('conv_bcparms -- bc_list options not yet ready')

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
        elseif (bcchar(ib).eq.'E' .or. bcchar(ib).eq.'W') then
          if (bottomiy(bcpos(ib),bcstart(ib)) .eq. -2)
     >      bcs(ib) = bcs(ib) + 1
          if (topiy   (bcpos(ib),bcend(ib)  ) .eq. ny+1)
     >      bce(ib) = bce(ib) - 1
        else
          call xerrab ('bchar not yet available for WG')
        endif
        bcl(ib) = bce(ib) - bcs(ib) + 1
        m%mxnBc = m%mxnBc + bcl(ib)
      enddo
      if (m%mxnBc .ne. m%nCg ) 
     >  call xerrab ('conv_bcparms -- wrong mxnBc')


      call alloc_mapping_bc (m)

      nfx = (nx+1)*ny     ! number of x-faces
      inbc = 0 ! index in bcCv
      do ib = 1, nbc
        m%bcFcP(ib,1) = inbc + 1
        m%bcFcP(ib,2) = bcl(ib)
        if (bcchar(ib).eq.'S') then
          do ibc = 1, bcl(ib)
            inbc = inbc + 1
            m%bcFc(inbc) = nfx + (bcs(ib)+ibc-1)*(ny+1) + 1    ! number of the corresponding domain face
            m%bcFcOr(inbc) = -1.0_R8
          enddo
        elseif (bcchar(ib).eq.'N') then
          do ibc = 1, bcl(ib)
            inbc = inbc + 1
            m%bcFc(inbc) = nfx + (bcs(ib)+ibc-1)*(ny+1) + ny+1
            m%bcFcOr(inbc) = 1.0_R8
          enddo
        elseif (bcchar(ib).eq.'W') then
          do ibc = 1, bcl(ib)
            inbc = inbc + 1
            m%bcFc(inbc) = (bcs(ib)+ibc-1)*(nx+1) + 1
            m%bcFcOr(inbc) = -1.0_R8
          enddo
        elseif (bcchar(ib).eq.'E') then 
          do ibc = 1, bcl(ib)
            inbc = inbc + 1
            m%bcFc(inbc) = (bcs(ib)+ibc-1)*(nx+1) + nx+1
            m%bcFcOr(inbc) = 1.0_R8
          enddo
        else
          call xerrab ('conv_bcparms -- wrong BCCHAR')
        endif
      enddo     

      if (inbc .ne. m%nCg ) 
     >  call xerrab ('conv_bcparms -- problem setting up bc arrays')

      return
      end subroutine conv_bcparms


***********************************************************************

      subroutine write_bcparms (ns, m)
      use b2mod_boundary_namelist
      implicit none
      integer, intent(in) :: ns
      type (mapping), intent (inout) :: m
      integer :: ib, iout, ian, ien, is, i
      character(1) :: ap =''''
      character(10) :: new_boundary_version = '03.002.000'

      open(newunit=iout,file='b2us.boundary.parameters')

      call cfverw(iout,new_boundary_version)

      write (iout,'(1x,a)') '&boundary'
      write (iout,'(1x,a,i0)') 'NBC     = ',m%nBc
      write (iout,'(1x,a,(10(2x,4a1)/))') 
     .       'BCCHAR  =',(ap,bcchar(ib),ap,',',ib=1,m%nBc)

      write (iout,'(1x,a,(10(1x,i4,a1)/))') 
     .       'NFACES  =',(m%bcFcP(ib,2),',',ib=1,m%nBc)
      do ib = 1, m%nBc
        ian = m%bcFcP(ib,1)
        ien = ian + m%bcFcP(ib,2) - 1
        write (iout,'(1x,a,i2,a,(1x,T15,10(i5,a1,1x)))') 
     .    'FACES(1,',ib,')= ',(m%bcFc(i),',',i=ian,ien)
      end do

      do ib = 1, m%nBc
        ian = m%bcFcP(ib,1)
        ien = ian + m%bcFcP(ib,2) - 1
        write (iout,'(1x,a,i2,a,(1x,T18,10(f5.2,a1,1x)))') 
     .    'FACE_ORI(1,',ib,')= ',(m%bcFcOr(i),',',i=ian,ien)
      end do

      write (iout,'(1x,a,(10(1x,i9,a1)/))') 
     .       'BCENE   =   ',(bcene(ib),',',ib=1,m%nBc)
      write (iout,'(1x,a,(10(es10.2,a1)/))') 
     .       'ENEPAR(1,1)=',(enepar(ib,1),',',ib=1,m%nBc)
      write (iout,'(1x,a,(10(es10.2,a1)/))') 
     .       'ENEPAR(1,2)=',(enepar(ib,2),',',ib=1,m%nBc)

      write (iout,'(1x,a,(10(1x,i9,a1)/))') 
     .       'BCENI   =   ',(bceni(ib),',',ib=1,m%nBc)
      write (iout,'(1x,a,(10(es10.2,a1)/))') 
     .       'ENIPAR(1,1)=',(enipar(ib,1),',',ib=1,m%nBc)
      write (iout,'(1x,a,(10(es10.2,a1)/))') 
     .       'ENIPAR(1,2)=',(enipar(ib,2),',',ib=1,m%nBc)

      write (iout,'(1x,a,(10(1x,i9,a1)/))') 
     .      'BCPOT   =   ',(bcpot(ib),',',ib=1,m%nBc)
      write (iout,'(1x,a,(10(es10.2,a1)/))') 
     .       'POTPAR(1,1)=',(potpar(ib,1),',',ib=1,m%nBc)
      write (iout,'(1x,a,(10(es10.2,a1)/))') 
     .       'POTPAR(1,2)=',(potpar(ib,2),',',ib=1,m%nBc)

      do ib = 1, m%nBc
        write (iout,'(1x,a,i2,a,(10(1x,i9,a1)/))') 
     .       'BCCON(0,',ib,') =  ',(bccon(is,ib),',',is=0,ns-1)
        write (iout,'(1x,a,i2,a,(10(es10.2,a1)/))') 
     .       'CONPAR(0,',ib,',1)=',(conpar(is,ib,1),',',is=0,ns-1)
        write (iout,'(1x,a,i2,a,(10(es10.2,a1)/))') 
     .       'CONPAR(0,',ib,',2)=',(conpar(is,ib,2),',',is=0,ns-1)
      end do

      do ib = 1, m%nBc
        write (iout,'(1x,a,i2,a,(10(1x,i9,a1)/))') 
     .       'BCMOM(0,',ib,') =  ',(bcmom(is,ib),',',is=0,ns-1)
        write (iout,'(1x,a,i2,a,(10(es10.2,a1)/))') 
     .       'MOMPAR(0,',ib,',1)=',(mompar(is,ib,1),',',is=0,ns-1)
        write (iout,'(1x,a,i2,a,(10(es10.2,a1)/))') 
     .       'MOMPAR(0,',ib,',2)=',(mompar(is,ib,2),',',is=0,ns-1)
      end do

      write (iout,'(1x,a,f10.7,a1)') 'GAMMAI=',gammai,','
      write (iout,'(1x,a,f10.7,a1)') 'GAMMAE=',gammae,','
      write (iout,'(1x,a,l1,a1)') 'LBNDUSR=',lbndusr,','
      write (iout,'(1x,a,l1,a1)') 'LFEEDBACK=',lfeedback,','
      
      write (iout,'(1x,a,i0,a1)') 'NNISO = ',nniso,','
      
      if (nniso > 0) 
     .   call xerrab ('write_bcparms -- nniso > 0')        

      write (iout,'(1x,a,es10.2,a1)') 'TEISO = ',teiso,','
      write (iout,'(1x,a,es10.2,a1)') 'TIISO = ',tiiso,','
      write (iout,'(1x,a,es10.2,a1)') 'PHIISO = ',phiiso,','
      write (iout,'(1x,a,(10(es10.2,a1)/))') 
     .       'NIISO(0)=',(niiso(is),',',is=0,ns)

      write (iout,'(1x,a,(10(1x,i5,a1)/))') 
     .      'lcbs=',(lcbs(ib),',',ib=1,ncbs)

      write (iout,'(1x,5a)') 'BOUNDARY_FILENAME=',ap,
     .                       trim(boundary_filename),ap,','
      write (iout,'(1x,a,es10.2,a1)') 
     .       'BOUNDARY_TIME_MOD = ',boundary_time_mod,','
      write (iout,'(1x,a,es10.2,a1)') 
     .       'BOUNDARY_TIME_SWITCH = ',boundary_time_switch,','
      write (iout,'(1x,a,l1,a1)') 'WRITE_NML_BND=',write_nml_bnd,','
    
      write (iout,'(1x,a)') '/'
      close (iout)

      end subroutine write_bcparms

      end module b2us_prep

