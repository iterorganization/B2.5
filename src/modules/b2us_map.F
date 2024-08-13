      module b2us_map

      use b2mod_types
      use b2mod_version
      implicit none

      private
      public :: alloc_mapping, dealloc_mapping,
     .          alloc_mapping_bc, dealloc_mapping_bc,
     .          alloc_mapping_rc, dealloc_mapping_rc,
     .          alloc_mapping_cf, dealloc_mapping_cf,
     .          read_mapping, write_mapping,
     .          init_mapping, get_version,
     .          calc_add_connectivity

      type, public :: mapping
        integer :: nCv, nFc, nVx, nCg, nCi, nCmxVx, nCmxFc, nVmxCv,
     &             nVmxFc, nFs, nFt, nBc, mxnBc, nRc, mxnRc, nCmxNv,
     &             nCf, mxnCf, mxStencil, nFsVxmx, nXpt, nStr, nTgc
        integer :: nnreg(0:1), periodic_bc, iFssep, iFssep2

        integer, allocatable :: cvFcP(:,:)             ! (nCv,2) pointing for every cell to the first index number in the cvFc list and its number of faces
        integer, allocatable :: cvFc(:)                !         listing for every cell the corresponding face numbers
        integer, allocatable :: fcCv(:,:)              ! (nFc,2) listing for every face the two corresponding cells
        integer, allocatable :: fcVx(:,:)              ! (nFc,2) listing for every face the two corresponding vertices
        integer, allocatable :: cvVxP(:,:)             ! (nCv,2) pointing for every cell to the first index number in the cvVx list and its number of vertices
        integer, allocatable :: cvVx(:)                !         listing for every cell the corresponding vertex numbers
        integer, allocatable :: vxFcP(:,:)             ! (nVx,2) pointing for every vertex to the first index number in the vxFc list and its number of faces
        integer, allocatable :: vxFc(:)                !         listing for every vertex the corresponding face numbers
        integer, allocatable :: vxCvP(:,:)             ! (nVx,2) pointing for every cell to the first index number in the vxCv list and its number of cells
        integer, allocatable :: vxCv(:)                !         listing for every vertex the corresponding cell numbers
        integer, allocatable :: ftCvP(:,:)             ! (nFt,2) pointing for every flux tube to the first index number in the ftCv list and its number of cells
        integer, allocatable :: ftCv(:)                !         listing for every flux tube the corresponding cell numbers
        integer, allocatable :: ftFcP(:,:)             ! (nFt,2) pointing for every flux tube to the first index number in the ftFc list and its number of cell faces
        integer, allocatable :: ftFc(:)                !         listing for every flux tube the corresponding cell face numbers
        integer, allocatable :: cvFt(:)                !         listing for every cell the corresponding flux tube
        integer, allocatable :: fsFcP(:,:)             ! (nFs,2) pointing for every flux surface to the first index number in the fsFc list and its number of cell faces
        integer, allocatable :: fsFc(:)                !         listing for every flux surface the corresponding cell face numbers
        integer, allocatable :: fcFs(:)                ! for every cell face the corresponding flux surface number, if any

        integer, allocatable :: fcReg(:)               ! Region numbers of the cell faces
        integer, allocatable :: cvReg(:)               ! Region numbers of the cells
        integer, allocatable :: ftReg(:)               ! Region numbers of the flux tubes
        integer, allocatable :: fcAligned(:)           ! Integer to indicate whether a face is aligned with the magnetic field
        logical, allocatable :: cvOnClosedSurface(:)   ! Logical to indicate whether a cell is on a closed flux surface

        integer, allocatable :: bcCvP(:,:)             ! (nBc,2) pointing for each boundary condition to the first index in the bcCv list, and its number of guard cells
        integer, allocatable :: bcCv(:,:)              !         listing for each boundary condition the guard cells, and the corresponding domain cells
        integer, allocatable :: bcFcP(:,:)             ! (nBc,2) pointing for each boundary condition to the first index in the bcFc list, and its number of boundary faces
        integer, allocatable :: bcFc(:)                !         listing for each boundary condition the boundary faces
        real (kind=R8), allocatable :: bcFcOr(:)       ! Multiplier to indicate whether outward normal coincides with the face normal (1.0_R8) or not (-1.0_R8)

        integer, allocatable :: rcCvP(:,:)             ! (nRc,2) pointing for each stratum to the first index in the rcCv list, and its number of guard cells
        integer, allocatable :: rcCv(:,:)              !         listing for each stratum the guard cells, and the corresponding domain cells
        integer, allocatable :: rcFcP(:,:)             ! (nRc,2) pointing for each stratum to the first index in the rcFc list, and its number of boundary faces
        integer, allocatable :: rcFc(:)                !         listing for each stratum the boundary faces
        real (kind=R8), allocatable :: rcFcOr(:)       ! Multiplier to indicate whether outward normal coincides with the face normal (1.0_R8) or not (-1.0_R8)

        integer, allocatable :: cfRegP(:,:)            ! (nCf,2) pointing for each cost function to the first index in the cfReg list, and its number of items
        integer, allocatable :: cfReg(:)               !         listing for each cost function the corresponding domain cells or faces
        logical, allocatable :: cfOnCv(:)              ! (nCf) Logical to indicate whether a cf domain is on Cv (true) or Fc (false)
        real (kind=R8), allocatable :: cfFcOr (:)      ! Multiplier to indicate whether outward normal coincides with the face normal (1.0_R8) or not (-1.0_R8)

!        real (kind=R8), allocatable :: intfaceV(:,:)   ! volume factors for interpolation from centers to faces; currently: geo%fcVol
!        real (kind=R8), allocatable :: intvertexV(:)   ! volume factors for interpolation from centers to vertices; currently: geo%vxVol
        real (kind=R8), allocatable :: intcellP(:)     ! factors for interpolation from faces to centers in poloidal direction
        real (kind=R8), allocatable :: intcellR(:)     ! factors for interpolation from faces to centers in radial direction

        integer, allocatable :: cvNvP(:,:)             !(nCv,2): pointing for each cell to the start in cvNv array, and the number of points in its stencil
        integer, allocatable :: cvNv(:)                !         listing for each cell the numbers of the cells in its stencil

        integer, allocatable :: fsVxP(:,:)             !(nFs,2): for each flux surface, starting point of its list of vertices in fsVx array [(,1) component], and the number of vertices on the flux surface [(,2) component]
        integer, allocatable :: fsVx(:)                ! list of vertices belonging to the flux surfaces
        integer, allocatable :: vxFs(:)                ! for each vertex, the flux surface it belongs to

        integer, allocatable :: Xpt(:)                 ! vertex numbers for the X-points
        integer, allocatable :: strVxP(:,:)            !(nXpt,2): for each Xpt, the list of associated strike point vertices, and their number
        integer, allocatable :: strVx(:)               ! list of vertex numbers for the strike points
        integer, allocatable :: tgVx(:)                ! list of vertex numbers for the tangency points
        integer, allocatable :: strDiv(:)              ! For every strike point, the index of the divertor target on which it is located
                                                       ! The divertor targets are labelled from 1 to 4 starting travelling from West to East in the computational domain
        integer, allocatable :: divFcP(:,:)            ! For each divertor target, the list of associated divertor faces, and their number
        integer, allocatable :: divFc(:)               ! ORDERED list of faces forming the divertor targets
        real (kind=R8), allocatable :: divFcOr(:)      ! Multiplier to indicate whether outward normal coincides with the face normal (1.0_R8) or not (-1.0_R8)
        integer, allocatable :: ifdiv(:)               ! List index of the face immediately outside the strike point
        integer, allocatable :: ivdiv(:)               ! Index of the active strike-point vertex on the divertor target

        integer :: isClassicalGrid                     ! indicates whether grid is derived from an "original", structured grid
        integer :: nx, ny, nncut, nctcnt               ! some parameters related to the original grid
        integer, allocatable :: imapCv(:,:)            ! back projection of cv indices to structured grid
        integer, allocatable :: imapFcx(:,:)           ! back projection of fc indices for x-faces to structured grid
        integer, allocatable :: imapFcy(:,:)           ! back projection of fc indices for y-faces to structured grid
        integer, allocatable :: imapVx(:,:)            ! back projection of Vx indices to structured grid

        integer, allocatable :: icornVx(:)             ! corner indices; to be used in converting target information
        integer, allocatable :: ictcntrl(:)            ! indices of central cuts; to be used in converting triangle information

        integer, allocatable :: fcLbl(:)               ! (nFc) indicates for each face the corresponding face label
        integer, allocatable :: ftLbl(:)               ! (nFt) indicates for each flux tube the corresponding label
        integer, allocatable :: cvLbl(:)               ! (nCv) indicates for each control volume the corresponding label
        end type mapping

      contains

**********************************************************************

      subroutine alloc_mapping (m)
#ifdef WG_DEBUG
      use b2us_debug
#endif

      implicit none
      type(mapping), intent(inout) :: m

      if (allocated(m%cvFcP)) return

      allocate (m%cvFcP(m%nCv,2))
      allocate (m%cvFc(m%nCmxFc))
      allocate (m%fcCv(m%nFc,2))
      allocate (m%fcVx(m%nFc,2))
      allocate (m%cvVxP(m%nCv,2))
      allocate (m%cvVx(m%nCmxVx))
      allocate (m%vxFcP(m%nVx,2))
      allocate (m%vxFc(m%nVmxFc))
      allocate (m%vxCvP(m%nVx,2))
      allocate (m%vxCv(m%nVmxCv))
      allocate (m%ftCv(m%nCv))
      allocate (m%ftCvP(m%nFt,2))
      allocate (m%ftFc(m%nFc))
      allocate (m%ftFcP(m%nFt,2))
      allocate (m%cvFt(m%nCv))
      allocate (m%fsFc(m%nFc))
      allocate (m%fsFcP(m%nFs,2))
      allocate (m%fcFs(m%nFc))
      allocate (m%fcReg(m%nFc))
      allocate (m%cvReg(m%nCv))
      allocate (m%ftReg(m%nFt))
      allocate (m%fcAligned(m%nFc))
      allocate (m%cvOnClosedSurface(m%nCv))

!      allocate (m%intfaceV(m%nFc,2))
!      allocate (m%intvertexV(m%nCmxVx * m%nCv))
      allocate (m%intcellP(m%nCmxFc))
      allocate (m%intcellR(m%nCmxFc))

!      allocate (m%cvNvP(m%nCv,2))
!      allocate (m%cvNv(m%nCmxNv))

      allocate (m%fsVxP(m%nFs,2))
      allocate (m%vxFs(m%nVx))

      if (m%isClassicalGrid.eq.1) then
        allocate (m%imapCv(-1:m%nx,-1:m%ny))
        allocate (m%imapFcx(-1:m%nx,-1:m%ny))
        allocate (m%imapFcy(-1:m%nx,-1:m%ny))
        allocate (m%imapVx(-1:m%nx,-1:m%ny))
        allocate (m%icornVx(1:m%nx))
#ifdef WG_DEBUG
        allocate (imapCv(-1:m%nx,-1:m%ny))
        allocate (imapFcx(-1:m%nx,-1:m%ny))
        allocate (imapFcy(-1:m%nx,-1:m%ny))
        allocate (imapVx(-1:m%nx,-1:m%ny))
#endif
      else
        allocate (m%imapCv(0,0))
        allocate (m%imapFcx(0,0))
        allocate (m%imapFcy(0,0))
        allocate (m%imapVx(0,0))
        allocate (m%icornVx(0))
#ifdef WG_DEBUG
        allocate (imapCv(0,0))
        allocate (imapFcx(0,0))
        allocate (imapFcy(0,0))
        allocate (imapVx(0,0))
#endif
      endif

      allocate (m%fcLbl(m%nFc))
      allocate (m%ftLbl(m%nFt))
      allocate (m%cvLbl(m%nCv))

      call initialize_mapping (m)

      write (*,*) '=================================='
      write (*,*) 'Mapping allocated'
      write (*,*) 'nCv, nCi, nCg ', m%nCv, m%nCi, m%nCg
      write (*,*) 'nFc, nVx, nFt ', m%nFc, m%nVx, m%nFt
      write (*,*) 'nFs           ', m%nFs
      write (*,*) 'nCmxVx, nCmxFc', m%nCmxVx, m%nCmxFc
      write (*,*) 'nVmxCv, nVmxFc', m%nVmxCv, m%nVmxFc
      write (*,*) '=================================='

      return
      end subroutine alloc_mapping

**********************************************************************

      subroutine initialize_mapping (m)
#ifdef WG_DEBUG
      use b2us_debug
#endif
      implicit none
      type(mapping), intent(inout) :: m

      m%cvFcP = 0
      m%cvFc = 0
      m%fcCv = 0
      m%fcVx = 0
      m%cvVxP = 0
      m%cvVx = 0
      m%vxFcP = 0
      m%vxFc = 0
      m%vxCvP = 0
      m%vxCv = 0
      m%ftCv = 0
      m%ftCvP = 0
      m%ftFc = 0
      m%ftFcP = 0
      m%cvFt = 0
      m%fsFc = 0
      m%fsFcP = 0
      m%fcFs = 0
      m%fcReg = 0
      m%cvReg = 0
      m%ftReg = 0
      m%fcAligned = 0
      m%cvOnClosedSurface = .false.

!      m%intfaceV = 0._R8
!      m%intvertexV = 0._R8
      m%intcellP = 0._R8
      m%intcellR = 0_R8

!      m%cvNvP = 0
!      m%cvNv = 0

      m%fsVxP = 0
      m%vxFs = 0

      if (m%isClassicalGrid.eq.1) then
        m%imapCv = 0
        m%imapFcx = 0
        m%imapFcy = 0
        m%imapVx = 0
        m%icornVx = 0
#ifdef WG_DEBUG
        imapCv = 0
        imapFcx = 0
        imapFcy = 0
        imapVx = 0
#endif
      endif

      m%fcLbl = 0
      m%ftLbl = 0
      m%cvLbl = 0

      m%iFssep = 0
      m%iFssep2 = 0

      write (*,*) '=================================='
      write (*,*) 'Mapping initialized'
      write (*,*) '=================================='

      return
      end subroutine initialize_mapping

**********************************************************************

      subroutine dealloc_mapping (m)
#ifdef WG_DEBUG
      use b2us_debug
#endif

      implicit none
      type(mapping), intent(in out) :: m

      call dealloc_mapping_bc(m)
      call dealloc_mapping_rc(m)
      call dealloc_mapping_cf(m)

      if (.not.allocated(m%cvFcP)) return

      deallocate (m%cvFcP)
      deallocate (m%cvFc)
      deallocate (m%fcCv)
      deallocate (m%fcVx)
      deallocate (m%cvVxP)
      deallocate (m%cvVx)
      deallocate (m%vxFcP)
      deallocate (m%vxFc)
      deallocate (m%vxCvP)
      deallocate (m%vxCv)
      deallocate (m%ftCv)
      deallocate (m%ftCvP)
      deallocate (m%ftFc)
      deallocate (m%ftFcP)
      deallocate (m%cvFt)
      deallocate (m%fsFc)
      deallocate (m%fsFcP)
      deallocate (m%fcFs)
      deallocate (m%fcReg)
      deallocate (m%cvReg)
      deallocate (m%ftReg)
      deallocate (m%fcAligned)
      deallocate (m%cvOnClosedSurface)

!      deallocate (m%intfaceV)
!      deallocate (m%intvertexV)
      deallocate (m%intcellP)
      deallocate (m%intcellR)

      if (allocated(m%cvNvP)) deallocate (m%cvNvP)
      if (allocated(m%cvNv))  deallocate (m%cvNv)

      deallocate (m%fsVxP)
      if (allocated(m%fsVx)) deallocate (m%fsVx)
      deallocate (m%vxFs)

      if (allocated(m%Xpt)) deallocate (m%Xpt)
      if (allocated(m%strVxP)) deallocate (m%strVxP)
      if (allocated(m%strVx)) deallocate (m%strVx)
      if (allocated(m%tgVx)) deallocate (m%tgVx)
      if (allocated(m%strDiv)) deallocate (m%strDiv)
      if (allocated(m%divFcP)) deallocate (m%divFcP)
      if (allocated(m%divFc)) deallocate (m%divFc)
      if (allocated(m%divFcOr)) deallocate (m%divFcOr)
      if (allocated(m%ifdiv)) deallocate (m%ifdiv)
      if (allocated(m%ivdiv)) deallocate (m%ivdiv)

      deallocate (m%imapCv)
      deallocate (m%imapFcx)
      deallocate (m%imapFcy)
      deallocate (m%imapVx)
      deallocate (m%icornVx)
      if (allocated(m%ictcntrl)) deallocate (m%ictcntrl)
#ifdef WG_DEBUG
      deallocate (imapCv)
      deallocate (imapFcx)
      deallocate (imapFcy)
      deallocate (imapVx)
#endif

      deallocate (m%fcLbl)
      deallocate (m%ftLbl)
      deallocate (m%cvLbl)

      return
      end subroutine dealloc_mapping

**********************************************************************

      subroutine alloc_mapping_bc (m)

      implicit none
      type(mapping), intent(in out) :: m

      if (allocated(m%bcCvP)) return

      allocate (m%bcCvP  ( m % nBc,   2 ) )
      allocate (m%bcCv   ( m % mxnBc, 2 ) )
      allocate (m%bcFcP  ( m % nBc,   2 ) )
      allocate (m%bcFc   ( m % mxnBc    ) )
      allocate (m%bcFcOr ( m % mxnBc    ) )

      m%bcCvP = 0
      m%bcCv = 0
      m%bcFcP = 0
      m%bcFc = 0
      m%bcFcOr = 0._R8

      return
      end subroutine alloc_mapping_bc

**********************************************************************

      subroutine dealloc_mapping_bc (m)

      implicit none
      type(mapping), intent(in out) :: m

      if (.not.allocated(m%bcCvP)) return

      deallocate ( m % bcCvP  )
      deallocate ( m % bcCv   )
      deallocate ( m % bcFcP  )
      deallocate ( m % bcFc   )
      deallocate ( m % bcFcOr )

      return
      end subroutine dealloc_mapping_bc

**********************************************************************

      subroutine alloc_mapping_rc (m)

      implicit none
      type(mapping), intent(in out) :: m

      if (allocated(m%rcCvP)) return

      allocate (m%rcCvP  ( m % nRc,   2 ) )
      allocate (m%rcCv   ( m % mxnRc, 2 ) )
      allocate (m%rcFcP  ( m % nRc,   2 ) )
      allocate (m%rcFc   ( m % mxnRc    ) )
      allocate (m%rcFcOr ( m % mxnRc    ) )

      m%rcCvP = 0
      m%rcCv = 0
      m%rcFcP = 0
      m%rcFc = 0
      m%rcFcOr = 0._R8

      return
      end subroutine alloc_mapping_rc

**********************************************************************

      subroutine dealloc_mapping_rc (m)

      implicit none
      type(mapping), intent(in out) :: m

      if (.not.allocated(m%rcCvP)) return

      deallocate ( m % rcCvP  )
      deallocate ( m % rcCv   )
      deallocate ( m % rcFcP  )
      deallocate ( m % rcFc   )
      deallocate ( m % rcFcOr )

      return
      end subroutine dealloc_mapping_rc

**********************************************************************

      subroutine alloc_mapping_cf (m)

      implicit none
      type(mapping), intent(in out) :: m

      if (allocated(m%cfRegP)) return

      allocate (m%cfRegP  ( m % nCf,   2 ) )
      allocate (m%cfReg   ( m % mxnCf    ) )
      allocate (m%cfOnCv  ( m % nCf      ) )
      allocate (m%cfFcOr  ( m % mxnCf    ) )

      m%cfRegP = 0
      m%cfReg = 0
      m%cfOnCv = .true.
      m%cfFcOr = 0._R8

      return
      end subroutine alloc_mapping_cf

**********************************************************************

      subroutine dealloc_mapping_cf (m)

      implicit none
      type(mapping), intent(in out) :: m

      if (.not.allocated(m%cfRegP)) return

      deallocate ( m % cfRegP  )
      deallocate ( m % cfReg   )
      deallocate ( m % cfOnCv )
      deallocate ( m % cfFcOr   )

      return
      end subroutine dealloc_mapping_cf

**********************************************************************

      subroutine read_mapping (iun, mpg)
#ifdef WG_DEBUG
      use b2us_debug
#endif
      implicit none
      integer, intent(in) :: iun
      integer :: n2, nx
      type (mapping), intent (inout) :: mpg
      external cfruin, cfrure

#ifdef WG_DEBUG
      nx_st = mpg%nx
      ny_st = mpg%ny
      isClassicalGrid = mpg%isClassicalGrid
#endif

      call alloc_mapping (mpg)

! write the mapping information
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
      call cfruin (iun, mpg%nFt*2,   mpg%ftCvP, 'ftCvP')
      call cfruin (iun, mpg%nCv,     mpg%ftCv,  'ftCv')
      call cfruin (iun, mpg%nFt*2,   mpg%ftFcP, 'ftFcP')
      call cfruin (iun, mpg%nFc,     mpg%ftFc,  'ftFc')
      call cfruin (iun, mpg%nCv,     mpg%cvFt,  'cvFt')
      call cfruin (iun, mpg%nFs*2,   mpg%fsFcP, 'fsFcP')
      call cfruin (iun, mpg%nFc,     mpg%fsFc,  'fsFc')
      call cfruin (iun, mpg%nFc,     mpg%fcReg, 'fcReg')
      call cfruin (iun, mpg%nCv,     mpg%cvReg, 'cvReg')
      call cfruin (iun, mpg%nFt,     mpg%ftReg, 'ftReg')
      call cfrure (iun, mpg%nCmxFc,  mpg%intcellP, 'intcellP')
      call cfrure (iun, mpg%nCmxFc,  mpg%intcellR, 'intcellR')

      n2 = 0
      nx = 0
      if (mpg%isClassicalGrid == 1) n2 = (mpg%nx+2)*(mpg%ny+2)
      if (mpg%isClassicalGrid == 1) nx =  mpg%nx
      call cfruin (iun, n2, mpg%imapCv,  'imapCv')
      call cfruin (iun, n2, mpg%imapFcx, 'imapFcx')
      call cfruin (iun, n2, mpg%imapFcy, 'imapFcy')
      call cfruin (iun, n2, mpg%imapVx,  'imapVx')
      call cfruin (iun, nx, mpg%icornVX, 'icornVx')

      call cfruin (iun, mpg%nFc, mpg%fcLbl, 'fcLbl')
      call cfruin (iun, mpg%nCv, mpg%cvLbl, 'cvLbl')
      call cfruin (iun, mpg%nFt, mpg%ftLbl, 'ftLbl')

#ifdef WG_DEBUG
      imapCv  = mpg%imapCv
      imapFcx = mpg%imapFcx
      imapFcy = mpg%imapFcy
      imapVx  = mpg%imapVx
#endif

      call init_mapping (mpg)

      return
      end subroutine read_mapping

***********************************************************************

      subroutine write_mapping (nout, mpg)

      implicit none
      integer, intent(in) :: nout
      integer :: idum(0:9), n2
      type (mapping), intent (in) :: mpg
      external cfwuin, cfwure

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
      call cfwuin (nout, mpg%nFt*2,   mpg%ftCvP, 'ftCvP')
      call cfwuin (nout, mpg%nCv,     mpg%ftCv,  'ftCv')
      call cfwuin (nout, mpg%nFt*2,   mpg%ftFcP, 'ftFcP')
      call cfwuin (nout, mpg%nFc,     mpg%ftFc,  'ftFc')
      call cfwuin (nout, mpg%nCv,     mpg%cvFt,  'cvFt')
      call cfwuin (nout, mpg%nFs*2,   mpg%fsFcP, 'fsFcP')
      call cfwuin (nout, mpg%nFc,     mpg%fsFc,  'fsFc')
      call cfwuin (nout, mpg%nFc,     mpg%fcReg, 'fcReg')
      call cfwuin (nout, mpg%nCv,     mpg%cvReg, 'cvReg')
      call cfwuin (nout, mpg%nFt,     mpg%ftReg, 'ftReg')
      call cfwure (nout, mpg%nCmxFc,  mpg%intcellP, 'intcellP')
      call cfwure (nout, mpg%nCmxFc,  mpg%intcellR, 'intcellR')

      if (mpg%isClassicalGrid.eq.1) then
        n2 = (mpg%nx+2)*(mpg%ny+2)
        call cfwuin (nout, n2, mpg%imapCv,  'imapCv')
        call cfwuin (nout, n2, mpg%imapFcx, 'imapFcx')
        call cfwuin (nout, n2, mpg%imapFcy, 'imapFcy')
        call cfwuin (nout, n2, mpg%imapVx,  'imapVx')
        call cfwuin (nout, mpg%nx, mpg%icornVx, 'icornVx')
      else
        call cfwuin (nout, 0, idum, 'imapCv')
        call cfwuin (nout, 0, idum, 'imapFcx')
        call cfwuin (nout, 0, idum, 'imapFcy')
        call cfwuin (nout, 0, idum, 'imapVx')
        call cfwuin (nout, 0, idum, 'icornVx')
      endif

      call cfwuin (nout, mpg%nFc, mpg%fcLbl, 'fcLbl')
      call cfwuin (nout, mpg%nCv, mpg%cvLbl, 'cvLbl')
      call cfwuin (nout, mpg%nFt, mpg%ftLbl, 'ftLbl')

      return
      end subroutine write_mapping

**********************************************************************

      subroutine init_mapping(m)
      implicit none
      type (mapping), intent (inout) :: m
      integer :: i, j, k, l, indd
      integer :: iCv, iCv1, iCv2, iFc, iFcc, iFs, iVa, iVb,
     *           iNv, iNv1, iNv2, iVx, iVx1, iVx2, iNvmx, inew
      integer :: totalcount, fscount, iactive, icount, fccount
      integer, allocatable :: indCv(:), indFsVx(:), indXpt(:),
     *                        indFc(:), cvNvLoc(:), freg(:),
     *                        old_face_list(:), new_face_list(:)
      logical :: duplicate, left_match_found, right_match_found
      external ipgeti, xerrab
cwdk  This initialization routine precomputes a number of arrays
cwdk  in the mapping that are not stored in the b2fgmtry file, but
cwdk  can be derived directly from data in that file

*   ..regions
      m%nnreg(0) = maxval(m%cvReg)
      m%nnreg(1) = maxval(m%fcReg)

*   ..set up connectivity information for the matrix stencil
      allocate (m%cvNvP(m%nCv,2))
      allocate (cvNvLoc(10*m%nCv)) ! initial overestimate of size cvNv array
      allocate (indCv(100))        ! (over)estimate of max number number of neighbors an individual cell could have

      iNv = 1
      do iCv = 1, m%nCv
        m%cvNvP(iCv,1) = iNv
        ! init: each cell appears in its own stencil
        indCv = iCv
        iNvmx = 1
        ! list all cells surrounding cell iCv, by listing the
        ! cells connecting to all vertices of cell iCv
        do iVx = 1, m%cvVxP(iCv,2)
          iVx1 = m%cvVx(m%cvVxP(iCv,1)+iVx-1)     ! the number in the vertex list
          indCv(iNvmx+1:iNvmx+m%vxCvP(iVx1,2)) =  ! all cells connecting to this vertex
     .      m%vxCv(m%vxCvP(iVx1,1):m%vxCvP(iVx1,1)+m%vxCvP(iVx1,2)-1)
          iNvmx = iNvmx + m%vxCvP(iVx1,2)
        enddo
        ! remove duplicates
        cvNvLoc(m%cvNvP(iCv,1)) = indCv(1)
        i = 1
        outer: do iNv1 = 2, iNvmx
          do iNv2 = 1, i
            if (cvNvLoc(m%cvNvP(iCv,1)+iNv2-1).eq.indCv(iNv1)) then
              ! duplicate element
              cycle outer
            endif
          enddo
          ! new element found
          i = i + 1
          cvNvLoc(m%cvNvP(iCv,1)+i-1) = indCv(iNv1)
        end do outer
        ! TODO: sort by cell number, keeping first cell fixed (main diagonal)
        m%cvNvP(iCv,2) = i
        iNv = iNv + i
      enddo

      ! allocate and fill m%cvNv
      m%nCmxNv = iNv - 1
      allocate (m%cvNv(m%nCmxNv))
      m%cvNv = cvNvLoc(1:m%nCmxNv)

!csc this is only needed for adjoint AD in b2usco, b2usmo, etc.
      m%mxStencil = maxval(m%cvNvP(1:m%nCv,2))

      ! clean up
      deallocate (cvNvLoc)
      deallocate (indCv)

      allocate (indFsVx(m%nVx))
      m%fcFs = 0
      indFsVx = 0
      totalcount = 0
      do iFs = 1, m%nFs
        fscount = 0
        do i = m%fsFcP(iFs,1), m%fsFcP(iFs,1) + m%fsFcP(iFs,2) - 1
          iFc = m%fsFc(i)
          m%fcFs(iFc) = iFs
          iVx1 = m%fcVx(iFc,1)
          iVx2 = m%fcVx(iFc,2)
          if (m%vxFs(iVx1).eq.0) then
            m%vxFs(iVx1) = iFs
            m%nFsVxmx = m%nFsVxmx + 1
            fscount = fscount + 1
            totalcount = totalcount + 1
            indFsVx(totalcount) = iVx1
          end if
          if (m%vxFs(iVx2).eq.0) then
            m%vxFs(iVx2) = iFs
            m%nFsVxmx = m%nFsVxmx + 1
            fscount = fscount + 1
            totalcount = totalcount + 1
            indFsVx(totalcount) = iVx2
          end if
        end do
        if (fscount.gt.0) then
          m%fsVxP(iFs,1) = totalcount - fscount + 1
          m%fsVxP(iFs,2) = fscount
        end if
      end do
      m%nFsVxmx = totalcount
      allocate(m%fsVx(m%nFsVxmx))
      m%fsVx(1:m%nFsVxmx) = indFsVx(1:m%nFsVxmx)
      deallocate(indFsVx)

      allocate(indXpt(m%nVx))
      indXpt = 0
      m%nXpt = 0
      iactive = 0
      do iVx = 1, m%nVx
        if (m%vxFcP(iVx,2).ge.4) then
          fccount = 0
          do i = m%vxFcP(iVx,1), m%vxFcP(iVx,1) + m%vxFcP(iVx,2) - 1
            iFc = m%vxFc(i)
            if (m%fcFs(iFc).eq.0) cycle
            fccount = fccount + 1
          end do
          if (fccount.eq.4) then
            m%nXpt = m%nXpt + 1
            indXpt(m%nXpt) = iVx
            if (iactive.eq.0) then
              do i = m%vxCvP(iVx,1), m%vxCvP(iVx,1) + m%vxCvP(iVx,2) - 1
                if (iactive.ne.0) cycle
                iCv = m%vxCv(i)
                if (m%cvReg(iCv).eq.1) iactive = m%nXpt
              end do
            end if
          end if
        end if
      end do

      m%nStr = 0
      m%nTgc = 0
      if (m%nXpt.gt.0) then
        allocate(m%Xpt(m%nXpt))
        m%Xpt(1:m%nXpt) = indXpt(1:m%nXpt)
        if (iactive.ne.0) then
          do i = 1, m%vxFcP(m%Xpt(iactive),2)
            iFs = m%fcFs(m%vxFc(m%vxFcP(m%Xpt(iactive),1)+i-1))
            if (iFs.eq.0) cycle
            if (m%iFssep.eq.0) m%iFssep = iFs
            if (m%iFssep2.eq.0.and.iFs.ne.m%iFssep) m%iFssep2 = iFs
          end do
        end if
        indXpt = 0
        icount = 0
        allocate(m%strVxP(m%nXpt,2))
        m%strVxP = 0
        m%strVxP(1,1) = 1
        do i = 1, m%nXpt
          if (i.gt.1) m%strVxP(i,1) = m%strVxP(i-1,1) + m%strVxP(i-1,2)
          iFs = m%iFssep
          do j = m%fsFcP(iFs,1), m%fsFcP(iFs,1) + m%fsFcP(iFs,2) - 1
            iFc = m%fsFc(j)
            do k = 1, 2
              iVx = m%fcVx(iFc,k)
              fccount = 0
              do l = m%vxFcP(iVx,1), m%vxFcP(iVx,1) + m%vxFcP(iVx,2) - 1
                if (m%fcFs(m%vxFc(l)).eq.iFs) fccount = fccount + 1
              end do
              if (fccount.eq.1) then
                icount = icount + 1
                m%strVxP(i,2) = m%strVxP(i,2) + 1
                indXpt(icount) = iVx
              end if
            end do
          end do
          if (m%iFssep2.gt.0) then
            iFs = m%iFssep2
            do j = m%fsFcP(iFs,1), m%fsFcP(iFs,1) + m%fsFcP(iFs,2) - 1
              iFc = m%fsFc(j)
              do k = 1, 2
                iVx = m%fcVx(iFc,k)
                fccount = 0
                do l = m%vxFcP(iVx,1),
     &                 m%vxFcP(iVx,1) + m%vxFcP(iVx,2) - 1
                  if (m%fcFs(m%vxFc(l)).eq.iFs) fccount = fccount + 1
                end do
                if (fccount.eq.1) then
                  icount = icount + 1
                  m%strVxP(i,2) = m%strVxP(i,2) + 1
                  indXpt(icount) = iVx
                end if
              end do
            end do
          end if
        end do
        m%nStr = icount
        allocate(m%strVx(m%nStr))
        m%strVx(1:m%nStr) = indXpt(1:m%nStr)
        allocate(m%strDiv(m%nStr))
        allocate(freg(m%nStr))
        m%strDiv = 0
        freg = 0
        do i = 1, m%nStr
          iVx = m%strVx(i)
          do j = m%vxFcP(iVx,1), m%vxFcP(iVx,1) + m%vxFcP(iVx,2) - 1
            if (m%strDiv(i).ne.0) cycle
            iFc = m%vxFc(j)
            if (m%fcReg(iFc).eq.1) then
              m%strDiv(i) = 1
              freg(i) = 1
            else if (m%fcReg(iFc).eq.4) then
              m%strDiv(i) = 2
              freg(i) = 4
            else if (m%fcReg(iFc).eq.5) then
              if (m%nXpt.eq.1) then
                m%strDiv(i) = 1
              else
                m%strDiv(i) = 3
              endif
              freg(i) = 5
            else if (m%fcReg(iFc).eq.8) then
              if (m%nXpt.eq.1) then
                m%strDiv(i) = 2
              else
                m%strDiv(i) = 4
              endif
              freg(i) = 8
            end if
          end do
        end do
      else if (m%nnreg(0).eq.2) then
        allocate(freg(2))
        iFc = 1
        icount = 0
        do while (m%iFssep.eq.0.and.iFc.le.m%nFc)
          iCv1 = m%fcCv(iFc,1)
          iCv2 = m%fcCv(iFc,2)
          if (m%cvReg(iCv1).ne.m%cvReg(iCv2).and.
     &        m%cvReg(iCv1).ne.0.and.m%cvReg(iCv2).ne.0) then
            m%iFssep = m%fcFs(iFc)
          end if
          if (m%iFssep.eq.0) iFc = iFc + 1
        end do
        do i = m%fsFcP(m%iFssep,1),
     &         m%fsFcP(m%iFssep,1) + m%fsFcP(m%iFssep,2) - 1
          iFc = m%fsFc(i)
          do k = 1, 2
            iVx = m%fcVx(iFc,k)
            do j = m%vxFcP(iVx,1), m%vxFcP(iVx,1) + m%vxFcP(iVx,2) - 1
              iFcc = m%vxFc(j)
              if (m%fcLbl(iFcc).ne.0) then
                duplicate = .false.
                do l = 1, icount
                  if (indXpt(l).eq.iVx) duplicate = .true.
                end do
                if (.not.duplicate) then
                  icount = icount + 1
                  indXpt(icount) = iVx
                end if
              end if
            end do
          end do
        end do
        m%nTgc = icount
        allocate(m%tgVx(m%nTgc))
        m%tgVx(1:m%nTgc) = indXpt(1:m%nTgc)
        allocate(m%strDiv(2))
        m%strDiv(1) = 1
        m%strDiv(2) = 2
        freg(1) = 1
        freg(2) = 2
      else if (m%nnreg(0).eq.1) then
        call ipgeti ('b2mwti_jsep', m%iFssep)
        allocate(freg(2))
        freg(1) = 1
        freg(2) = 2
        if (m%iFssep.gt.0) then
          indXpt = 0
          icount = 0
          allocate(m%strVxP(1,2))
          m%strVxP(1,1) = 1
          m%strVxP(1,2) = 2
          do j = m%fsFcP(m%iFssep,1),
     &           m%fsFcP(m%iFssep,1) + m%fsFcP(m%iFssep,2) - 1
            iFc = m%fsFc(j)
            do k = 1, 2
              iVx = m%fcVx(iFc,k)
              fccount = 0
              do l = m%vxFcP(iVx,1), m%vxFcP(iVx,1) + m%vxFcP(iVx,2) - 1
                if (m%fcFs(m%vxFc(l)).eq.m%iFssep) fccount = fccount + 1
              end do
              if (fccount.eq.1) then
                icount = icount + 1
                indXpt(icount) = iVx
              end if
            end do
          end do
          m%nStr = icount
          if (m%nStr.gt.0) then
            allocate(m%strVx(m%nStr))
            m%strVx(1:m%nStr) = indXpt(1:m%nStr)
            allocate(m%strDiv(m%nStr))
            m%strDiv = 0
            do i = 1, m%nStr
              iVx = m%strVx(i)
              do j = m%vxFcP(iVx,1), m%vxFcP(iVx,1) + m%vxFcP(iVx,2) - 1
                if (m%strDiv(i).ne.0) cycle
                iFc = m%vxFc(j)
                if (m%fcReg(iFc).eq.1) then
                  m%strDiv(i) = 1
                else if (m%fcReg(iFc).eq.2) then
                  m%strDiv(i) = 2
                end if
              end do
            end do
          end if
        end if
      end if
      deallocate(indXpt)

      if (allocated(m%strDiv)) then
       indd = maxval(m%strDiv)
       if (indd.gt.0) then
        allocate(m%ifdiv(maxval(m%strDiv)))
        allocate(m%ivdiv(maxval(m%strDiv)))
        allocate(m%divFcP(maxval(m%strDiv),2))
        allocate(indFc(m%nFc))
        icount = 0
        m%ifdiv = 0
        m%ivdiv = 0
        m%divFcP(1,1) = 1
        m%divFcP(:,2) = 0
        do i = 1, indd  ! first assemble all the faces
          if (i.gt.1) m%divFcP(i,1) = m%divFcP(i-1,1) + m%divFcP(i-1,2)
          do iFc = 1, m%nFc
            if (m%fcReg(iFc).eq.freg(i)) then
              icount = icount + 1
              m%divFcP(i,2) = m%divFcP(i,2) + 1
              indFc(icount) = iFc
            end if
          end do
        end do
        allocate(m%divFc(icount))
        allocate(m%divFcOr(icount))
        m%divFcOr = 0.0_R8
        m%divFc(1:icount) = indFc(1:icount)
        do i = 1, indd ! now re-order them
          allocate(old_face_list(m%divFcP(i,2)))
          allocate(new_face_list(m%divFcP(i,2)))
          old_face_list(1:m%divFcP(i,2))=
     &     m%divFc(m%divFcP(i,1):m%divFcP(i,1)+m%divFcP(i,2)-1)
          new_face_list(:) = 0
          new_face_list(1) = m%divFc(m%divFcP(i,1))
          old_face_list(1) = 0
          iVa = m%fcVx(m%divFc(m%divFcP(i,1)),1)
          iVb = m%fcVx(m%divFc(m%divFcP(i,1)),2)
          inew = 1
          do while (inew.lt.m%divFcP(i,2))
            left_match_found = .false.
            do j = 1, m%divFcP(i,2)
              if (left_match_found) cycle
              if (old_face_list(j).eq.0) cycle
              if (iVa.eq.m%fcVx(old_face_list(j),1)) then
                left_match_found = .true.
                iVa = m%fcVx(old_face_list(j),2)
              else if (iVa.eq.m%fcVx(old_face_list(j),2)) then
                left_match_found = .true.
                iVa = m%fcVx(old_face_list(j),1)
              end if
              if (left_match_found) then
                indFc(1:inew) = new_face_list(1:inew)
                new_face_list(1) = old_face_list(j)
                new_face_list(2:inew+1) = indFc(1:inew)
                old_face_list(j) = 0
                inew = inew + 1
              end if
            end do
            right_match_found = .false.
            do j = 1, m%divFcP(i,2)
              if (right_match_found) cycle
              if (old_face_list(j).eq.0) cycle
              if (iVb.eq.m%fcVx(old_face_list(j),1)) then
                right_match_found = .true.
                iVb = m%fcVx(old_face_list(j),2)
              else if (iVb.eq.m%fcVx(old_face_list(j),2)) then
                right_match_found = .true.
                iVb = m%fcVx(old_face_list(j),1)
              end if
              if (right_match_found) then
                new_face_list(inew+1) = old_face_list(j)
                old_face_list(j) = 0
                inew = inew + 1
              end if
            end do
            if (.not.(left_match_found.or.right_match_found).and.
     &          inew.lt.m%divFcP(i,2)) then
              do j = 1, inew
                iFc = new_face_list(j)
                iVa = m%fcVx(iFc,1)
                iVb = m%fcVx(iFc,2)
                write(6,'(a,4i6)')
     &        'Chain of faces and vertices in new list ',j,iFc,iVa,iVb
              end do
              do j = 1, m%divFcP(i,2)
                if (old_face_list(j).eq.0) cycle
                iFc = m%divFc(m%divFcP(i,1)+j-1)
                iVa = m%fcVx(iFc,1)
                iVb = m%fcVx(iFc,2)
                write(6,'(a,4i6)')
     &        'Left over faces and vertices in old list',j,iFc,iVa,iVb
              end do
              call xerrab ('Broken face chain along target')
            end if
          end do
          m%divFc(m%divFcP(i,1):m%divFcP(i,1)+m%divFcP(i,2)-1)=
     &     new_face_list(1:m%divFcP(i,2))
          deallocate(old_face_list)
          deallocate(new_face_list)
        end do
        deallocate(indFc)
       end if
      end if
      deallocate(freg)

      write (*,*) 'nXpt, nStr, nTgc ', m%nXpt, m%nStr, m%nTgc
      write (*,*) 'nFsVxmx          ', m%nFsVxmx

      return
      end subroutine init_mapping

**********************************************************************

      subroutine get_version(boundary_filename,version)
      implicit none
      logical :: file_ok
      character*260 boundary_filename
      character*10 version
      character*260 filename
      external find_file, xerrab

      filename=trim(boundary_filename)
      call find_file(filename,file_ok)
      if(file_ok) then
        open(99,file=filename)
        call cfverr(99,version)
        close(99)
      else
        call xerrab ('No '//trim(filename))
      endif
      return
      end subroutine get_version

**********************************************************************

!> Computes additional connectivity data based
!! on the connectivity data of the traduitoutb2us file
      subroutine calc_add_connectivity(m)
      implicit none

      !!.. arguments
      type(mapping) :: m

      !!.. local variables
      integer :: i, ic, ivx, ifc, n, vxFc_end, vxCv_end
      integer, allocatable :: indFc(:), indCv(:), fh(:), ch(:), fch(:),
     .                    indcvfc(:)
      intrinsic pack, count
      external xerrab

      !!Read out fcVx, cvFc, cvFcP, cvVx and cvVxP
      !!Construct vxCv, vxCvP, vxFcP, fcCv !(because of ghostcell, every face has two cells)

      !Construct vxFcP
      allocate (indFc(2*m%nFc))
      indFc(1:m%nFc)          = (/ (i, i=1,m%nFc) /)
      indFc(m%nFc+1:2*m%nFc)    = (/ (i, i=1,m%nFc) /)

      allocate (fh(2*m%nFc))
      fh(1:m%nFc)             = m%fcVx(1:m%nFc,1)
      fh(m%nFc+1:2*m%nFc)     = m%fcVx(1:m%nFc,2)

      allocate (indCv(m%nCmxVx))
      do ic = 1,m%nCv
         indCv(m%cvVxP(ic,1):m%cvVxP(ic,1)+m%cvVxP(ic,2)-1) = ic
      enddo

      allocate (ch(m%nCmxVx))
      ch(1:m%nCmxVx) = m%cvVx(1:m%nCmxVx)

      m%vxFcP(1,1) = 1
      m%vxCvP(1,1) = 1
      do ivx = 1, m%nVx-1
         m%vxFcP(ivx,2) = count (fh == ivx)
         m%vxFcP(ivx+1,1) = m%vxFcP(ivx,1) + m%vxFcP(ivx,2)
         if (m%vxFcP(ivx,2) > 0) then
             m%vxFc(m%vxFcP(ivx,1):m%vxFcP(ivx+1,1)-1) =
     .        pack(indFc,fh == ivx)
         endif

         m%vxCvP(ivx,2) = count(ch == ivx)
         m%vxCvP(ivx+1,1) = m%vxCvP(ivx,1) + m%vxCvP(ivx,2)
         if (m%vxCvP(ivx,2) > 0 ) then
            m%vxCv(m%vxCvP(ivx,1):m%vxCvP(ivx+1,1)-1) =
     .         pack( indCv, ch == ivx )
         endif
      enddo

      ivx = m%nVx
      m%vxFcP(ivx,2) = count(fh == ivx)
      vxFc_end = m%vxFcP(ivx,1) + m%vxFcP(ivx,2) -1
      if (m%vxFcP(ivx,2) > 0) then
          m%vxFc(m%vxFcP(ivx,1):vxFc_end) =
     .      pack( indFc, fh == ivx )
      end if
      m%vxCvP(ivx,2) = count(ch == ivx)
      vxCv_end = m%vxCvP(ivx,1) + m%vxCvP(ivx,2) - 1
      if (m%vxCvP(ivx,2) > 0) then
         m%vxCv(m%vxCvP(ivx,1):vxCv_end) =
     .     pack( indCv, ch == ivx )
      endif

      !Try analog for fcCv(nFc,2), similar to vxCv - check this in debugger

      allocate (indcvfc(m%nCmxFc))
      do ic = 1,m%nCv
         indcvfc(m%cvFcP(ic,1):m%cvFcP(ic,1)+m%cvFcP(ic,2)-1) = ic
      enddo

      allocate (fch(m%nCmxFc))
      fch(1:m%nCmxFc) = m%cvFc(1:m%nCmxFc)

      !here loop over the faces

      do ifc = 1,m%nFc
           n = count(fch == ifc)
           if (n .eq. 2) then
             m%fcCv(ifc,1:2) = pack(indcvfc, fch == ifc)
           else
             call xerrab
     .        ('calc_add_connectivity: Wrong # of cells per face')
           endif
      enddo

      deallocate (indFc)
      deallocate (indCv)
      deallocate (indcvfc)
      deallocate (fh)
      deallocate (ch)
      deallocate (fch)

      return
      end subroutine calc_add_connectivity

**********************************************************************

      end module b2us_map
