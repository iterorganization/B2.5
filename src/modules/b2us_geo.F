      module b2us_geo

      use b2mod_types
      implicit none
      
      private
      public :: alloc_geometry, dealloc_geometry,
     .          read_geometry, write_geometry,
     .          init_geometry

      type, public :: geometry
        
!  cell quantites
        real (kind=R8), allocatable  :: cvBb(:,:)      ! Magnetic field: (nCv,0:3) Bx, By, Bz, |B|
        real (kind=R8), allocatable  :: cvX(:)         ! X-co of cell center
        real (kind=R8), allocatable  :: cvY(:)         ! Y-co of cell center
        real (kind=R8), allocatable  :: cvSz(:)        ! Area perp. to third (z) direction
        real (kind=R8), allocatable  :: cvHz(:)        ! Length in third (z) direction
        real (kind=R8), allocatable  :: cvHx(:)        ! Temp: characteristic poloidal length; remove?
        real (kind=R8), allocatable  :: cvQgam(:,:)    ! Temp: complementary angle between poloidal and radial directions; remove?
        real (kind=R8), allocatable  :: cvVol(:)       ! Cell volume
        real (kind=R8), allocatable  :: cvOnedBsq(:)   ! 1/B**2 in cell centers
        real (kind=R8), allocatable  :: cvBzb(:)       ! Factor for P-S/mdf flows
        real (kind=R8), allocatable  :: cvEb(:,:)      ! 

!  face based quantities
        real (kind=R8), allocatable  :: fcBb(:,:)      ! Magnetic field: (nFc,0:3) to be recomputed at cell faces (Bx,By,Bz,|B|)
        real (kind=R8), allocatable  :: fcS(:)         ! Total cell face area
        real (kind=R8), allocatable  :: fcHc(:,:)      ! Lengths between cell face center and neighboring cell centers ("connector" direction)
        real (kind=R8), allocatable  :: fcHt(:)        ! Lengths along cell face ("tangential" direction)
        real (kind=R8), allocatable  :: fcVol(:,:)     ! Part of volume of neighboring cell centers associated with the face; for interpolation
        real (kind=R8), allocatable  :: fcQgam(:,:)    ! Cosine (:,0) and sine (:,1) of gamma
        real (kind=R8), allocatable  :: fcQalf(:,:)    ! Cosine (:,0) and sine (:,1) of alpha
        real (kind=R8), allocatable  :: fcQbet(:,:)    ! Cosine (:,0) and sine (:,1) of beta
        real (kind=R8), allocatable  :: fcPbs(:)       ! Magnetic contact area
        real (kind=R8), allocatable  :: fcPbshz(:)     ! Magnetic contact area * hz (at faces)
        real (kind=R8), allocatable  :: fcEb(:,:)      ! 
        real (kind=R8), allocatable  :: vxBb(:,:)      ! Magnetic field: (nVx,0:3) to be recomputed at vertices (Bx,By,Bz,|B|)
        real (kind=R8), allocatable  :: vxX(:)         ! X-co of vertex
        real (kind=R8), allocatable  :: vxY(:)         ! Y-co of vertex
        real (kind=R8), allocatable  :: vxVol(:)       ! Part of volume of surrounding cell centers associated with the vertex; for interpolation
        real (kind=R8), allocatable  :: vxFfbz(:)      !
        real (kind=R8), allocatable  :: vxFpsi(:)      ! 
        real (kind=R8), allocatable  :: vxOnedBsq(:)   ! 1/B**2 at vertices

!  flux surface quantities
        real (kind=R8), allocatable  :: fsConn(:)      ! Connection length

      end type geometry

      contains

***********************************************************************

      subroutine alloc_geometry (g, nCv, nFc, nVx, nFs, nVmxCv)

      implicit none
      integer, intent(in) :: nCv, nFc, nVx, nFs, nVmxCv
      type(geometry), intent(in out) :: g

      if (allocated(g%cvBb)) return

!  cell quantites
      allocate (g%cvBb(nCv,0:3)) 
      allocate (g%cvX(nCv))
      allocate (g%cvY(nCv))
      allocate (g%cvSz(nCv))
      allocate (g%cvHz(nCv))
      allocate (g%cvHx(nCv))
      allocate (g%cvQgam(nCv,0:1))
      allocate (g%cvVol(nCv))
      allocate (g%cvOnedBsq(nCv))
      allocate (g%cvBzb(nCv))
      allocate (g%cvEb(nCv,0:2)) 

!  face quantites
      allocate (g%fcBb(nFc,0:3)) 
      allocate (g%fcS(nFc))
      allocate (g%fcHc(nFc,1:2))
      allocate (g%fcHt(nFc))
      allocate (g%fcVol(nFc,1:2))      
      allocate (g%fcQgam(nFc,0:1))
      allocate (g%fcQalf(nFc,0:1))
      allocate (g%fcQbet(nFc,0:1))
      allocate (g%fcPbs(nFc))
      allocate (g%fcPbshz(nFc))
      allocate (g%fcEb(nFc,0:2)) 

!  vertex quantites
      allocate (g%vxBb(nVx,0:3))
      allocate (g%vxX(nVx))
      allocate (g%vxY(nVx))
      allocate (g%vxVol(nVmxCv))
      allocate (g%vxFfbz(nVx))
      allocate (g%vxFpsi(nVx))
      allocate (g%vxOnedBsq(nVx))

!  flux surface quantites
      allocate (g%fsConn(nFs))
      
      return
      end subroutine alloc_geometry

***********************************************************************

      subroutine dealloc_geometry (g)

      implicit none
      type(geometry), intent(in out) :: g

      if (.not.allocated(g%cvBb)) return

!  cell quantites
      deallocate (g%cvBb) 
      deallocate (g%cvX)
      deallocate (g%cvY)
      deallocate (g%cvSz)
      deallocate (g%cvHz)
      deallocate (g%cvHx)
      deallocate (g%cvQgam)
      deallocate (g%cvVol)
      deallocate (g%cvOnedBsq)
      deallocate (g%cvBzb)
      deallocate (g%cvEb) 
      
!  face quantites
      deallocate (g%fcBb) 
      deallocate (g%fcS)
      deallocate (g%fcHc)
      deallocate (g%fcHt)
      deallocate (g%fcVol)      
      deallocate (g%fcQgam)
      deallocate (g%fcQalf)
      deallocate (g%fcQbet)
      deallocate (g%fcPbs)
      deallocate (g%fcPbshz)
      deallocate (g%fcEb) 
      
!  vertex quantites
      deallocate (g%vxBb)
      deallocate (g%vxX)
      deallocate (g%vxY)
      deallocate (g%vxVol)
      deallocate (g%vxFfbz)
      deallocate (g%vxFpsi)
      deallocate (g%vxOnedBsq)

!  flux surface quantites
      deallocate (g%fsConn)
      
      return
      end subroutine dealloc_geometry

***********************************************************************


      subroutine read_geometry (iun, nCv, nFc, nVx, nFs, nVmxCv, gm)

      implicit none
      integer, intent(in) :: iun, nCv, nFc, nVx, nFs, nVmxCv
      integer :: idum(0:9)
      type (geometry), intent (inout) :: gm

      call alloc_geometry (gm, nCv, nFc, nVx, nFs, nVmxCv)

! cell volumes 
      call cfrure (iun, nCv*4, gm%cvBb,   'cvBb')
      call cfrure (iun, nCv,   gm%cvX,    'cvX')
      call cfrure (iun, nCv,   gm%cvY,    'cvY')
      call cfrure (iun, nCv,   gm%cvSz,   'cvSz')
      call cfrure (iun, nCv,   gm%cvHz,   'cvHz')
      call cfrure (iun, nCv,   gm%cvHx,   'cvHx') !WG temp!
      call cfrure (iun, nCv*2, gm%cvQgam, 'cvQgam')
      call cfrure (iun, nCv,   gm%cvVol,  'cvVol')

! face quantities
      call cfrure (iun, nFc,   gm%fcS,    'fcS')
      call cfrure (iun, nFc*2, gm%fcHc,   'fcHc')
      call cfrure (iun, nFc,   gm%fcHt,   'fcHt')
      call cfrure (iun, nFc*2, gm%fcVol,  'fcVol')
      call cfrure (iun, nFc*2, gm%fcQgam, 'fcQgam')
      call cfrure (iun, nFc*2, gm%fcQalf, 'fcQalf')
      call cfrure (iun, nFc*2, gm%fcQbet, 'fcQbet')
      call cfrure (iun, nFc,   gm%fcPbs,  'fcPbs')

! vertex quantities
      call cfrure (iun, nVx,   gm%vxX,    'vxX')
      call cfrure (iun, nVx,   gm%vxY,    'vxY')
      call cfrure (iun, nVmxCv,gm%vxVol,  'vxVol')
      call cfrure (iun, nVx,   gm%vxFfbz, 'vxFfbz')
      call cfrure (iun, nVx,   gm%vxFpsi, 'vxFpsi')

      call cfrure (iun, nFs,   gm%fsConn, 'fsConn')

      call cfrure (iun, nCv*3, gm%cvEb,   'cvEb')
      call cfrure (iun, nFc*3, gm%fcEb,   'fcEb')

      return
      end subroutine read_geometry

***********************************************************************

      subroutine write_geometry (nout, nCv, nFc, nVx, nFs, nVmxCv, gm)

      implicit none
      integer, intent(in) :: nout, nCv, nFc, nVx, nFs, nVmxCv
      integer :: idum(0:9)
      type (geometry), intent (in) :: gm

! cell volumes 
      call cfwure (nout, nCv*4, gm%cvBb,   'cvBb')
      call cfwure (nout, nCv,   gm%cvX,    'cvX')
      call cfwure (nout, nCv,   gm%cvY,    'cvY')
      call cfwure (nout, nCv,   gm%cvSz,   'cvSz')
      call cfwure (nout, nCv,   gm%cvHz,   'cvHz')
      call cfwure (nout, nCv,   gm%cvHx,   'cvHx') !WG temp!
      call cfwure (nout, nCv*2, gm%cvQgam, 'cvQgam')
      call cfwure (nout, nCv,   gm%cvVol,  'cvVol')

! face quantities
      call cfwure (nout, nFc,   gm%fcS,    'fcS')
      call cfwure (nout, nFc*2, gm%fcHc,   'fcHc')
      call cfwure (nout, nFc,   gm%fcHt,   'fcHt')
      call cfwure (nout, nFc*2, gm%fcVol,  'fcVol')
      call cfwure (nout, nFc*2, gm%fcQgam, 'fcQgam')
      call cfwure (nout, nFc*2, gm%fcQalf, 'fcQalf')
      call cfwure (nout, nFc*2, gm%fcQbet, 'fcQbet')
      call cfwure (nout, nFc,   gm%fcPbs,  'fcPbs')

! vertex quantities
      call cfwure (nout, nVx,   gm%vxX,    'vxX')
      call cfwure (nout, nVx,   gm%vxY,    'vxY')
      call cfwure (nout, nVmxCv,gm%vxVol,  'vxVol')
      call cfwure (nout, nVx,   gm%vxFfbz, 'vxFfbz')
      call cfwure (nout, nVx,   gm%vxFpsi, 'vxFpsi')

      call cfwure (nout, nFs,   gm%fsConn, 'fsConn')

      call cfwure (nout, nCv*3, gm%cvEb,   'cvEb')
      call cfwure (nout, nFc*3, gm%fcEb,   'fcEb')

      return
      end subroutine write_geometry

***********************************************************************

      subroutine init_geometry (mpg, gm)
      use b2mod_constants
      use b2us_map
      implicit none
      type (mapping), intent (in) :: mpg
      type (geometry), intent (inout) :: gm
      integer :: nCv, nFc, nVx, nFx
      integer :: use_9pt_stencil = 1
      real (kind=R8) :: hzconst, fac_hz = 0.0_R8

      call ipgetr ('b2mndr_hz', fac_hz)
      call xertst (0.0.le.fac_hz .and. 1.0.ge.fac_hz,
     &  'faulty parameter b2mndr_hz')
      call ipgeti ('b2mndr_use_9pt_stencil',use_9pt_stencil)

      nCv = mpg%nCv
      nFc = mpg%nFc
      nVx = mpg%nVx

*   ..interpolate magnetic field to cell faces and vertices
      call intface1(nCv,nFc,2,mpg%fcCv,gm%fcVol,gm%cvBb,gm%fcBb)
      call intvertex1(nCv,nVx,2,mpg,gm%vxVol,gm%cvBb,gm%vxBb)

*   ..compute fcBb(:,3) and vxBb(:,3) components
      ! guarantee that bb(:,3)=sqrt(bb(:,0:2)**2)
      gm%fcBb(:,3) = 
     .  sqrt(gm%fcBb(:,0)**2 + gm%fcBb(:,1)**2 + gm%fcBb(:,2)**2)
      gm%vxBb(:,3) = 
     .  sqrt(gm%vxBb(:,0)**2 + gm%vxBb(:,1)**2 + gm%vxBb(:,2)**2)

*  .. calculate 1/B/B
      gm%cvOnedBsq = 1.0_R8/gm%cvBb(:,3)**2

*  .. interpolate OnedBsq for corners of cell
      call intvertex(nCv,nVx,mpg,gm%vxVol,gm%cvOnedBsq,gm%vxOnedBsq) !WG: interpolate 1/cvBb**2, or compute as 1/vxBb**2 ???

*  .. for the time being
      gm%cvEb(:,0) = 0.0_R8
      gm%cvEb(:,1) = 0.0_R8
      gm%cvEb(:,2) = 1.0_R8

      gm%fcEb(:,0) = 0.0_R8
      gm%fcEb(:,1) = 0.0_R8
      gm%fcEb(:,2) = 1.0_R8

*   ..compute cvHz
      hzconst = 2.0_R8*pi*6.0_R8   ! for ITER for testing
      hzconst = 2.0_R8*pi*1.6_R8   ! for ASDEX-UPGRADE for testing
      hzconst = 1.0_R8
      gm%cvHz = hzconst*(1.0e0_R8-fac_hz)+fac_hz*gm%cvVol/gm%cvSz

*   ..compute fcPbshz (magnetic field pitch * cell face area * hz)
      call intface(nCv,nFc,mpg%fcCv,gm%fcVol,gm%cvHz,gm%fcPbshz)
      gm%fcPbshz = gm%fcPbs * gm%fcPbshz

*   ..compute Bz/B**2 - Bz/<B**2>
      call b2xbzb (nCv, mpg, gm%cvBb, gm%cvVol, gm%cvBzb)

*   ..modify geometry to recover 5-point stencil behavior
      ! only supported for cases derived from original, structured grid
      if (use_9pt_stencil.eq.0) then
        if (mpg%isClassicalGrid.eq.1) then
          nfx = (mpg%nx + 1)*mpg%ny
          ! 'poloidal' faces
          gm%fcQgam(1:nfx,1) = 0.0_R8
          gm%fcQalf(1:nfx,1) = 0.0_R8
          ! 'radial' faces
          gm%fcQgam(nfx+1:mpg%nFc,1) = 0.0_R8
          gm%fcQbet(nfx+1:mpg%nFc,0) = 0.0_R8
        else
          call xerrab ('5-point stencil option not supported for UG')
        endif
      endif

      return
      end subroutine init_geometry


***********************************************************************

      end module b2us_geo
