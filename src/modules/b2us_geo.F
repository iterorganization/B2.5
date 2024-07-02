      module b2us_geo

      use b2mod_types
      implicit none

      private
      public :: alloc_geometry, dealloc_geometry,
     .          read_geometry, write_geometry,
     .          init_geometry, calc_conn, check_geometry,
     .          divide_by_area, split_face_flux, pitch

      type, public :: geometry

!  cell quantites
        real (kind=R8), allocatable  :: cvBb(:,:)      ! Magnetic field: (nCv,0:3) Bx, By, Bz, |B|
        real (kind=R8), allocatable  :: cvX(:)         ! X-co of cell center
        real (kind=R8), allocatable  :: cvY(:)         ! Y-co of cell center
        real (kind=R8), allocatable  :: cvSz(:)        ! Area perp. to third (z) direction
        real (kind=R8), allocatable  :: cvHz(:)        ! Length in third (z) direction
        real (kind=R8), allocatable  :: cvHx(:)        ! Temp: characteristic poloidal length; eventually remove, because meaning unclear in WG context
        real (kind=R8), allocatable  :: cvQgam(:,:)    ! Temp: complementary angle between poloidal and radial directions; remove?
        real (kind=R8), allocatable  :: cvVol(:)       ! Cell volume
        real (kind=R8), allocatable  :: cvOnedBsq(:)   ! 1/B**2 in cell centers
        real (kind=R8), allocatable  :: cvBzb(:)       ! Factor for P-S/mdf flows
        real (kind=R8), allocatable  :: cvEb(:,:)      ! Unit vector pointing along magnetic field (Cartesian X,Y,Z coordinates)
        real (kind=R8), allocatable  :: cvFpsi(:)      ! Psi in each cell

!  face based quantities
        real (kind=R8), allocatable  :: fcBb(:,:)      ! Magnetic field: (nFc,0:3) to be recomputed at cell faces (Bx,By,Bz,|B|)
        real (kind=R8), allocatable  :: fcS(:)         ! Total cell face area
        real (kind=R8), allocatable  :: fcHc(:,:)      ! Lengths between cell face center and neighboring cell centers ("connector" direction)
        real (kind=R8), allocatable  :: fcHt(:)        ! Lengths along cell face ("tangential" direction)
        real (kind=R8), allocatable  :: fcHz(:)        ! Length in third (z) direction
        real (kind=R8), allocatable  :: fcVol(:,:)     ! Part of volume of neighboring cell centers associated with the face; for interpolation
        real (kind=R8), allocatable  :: fcQgam(:,:)    ! Cosine (:,0) and sine (:,1) of gamma
        real (kind=R8), allocatable  :: fcQalf(:,:)    ! Cosine (:,0) and sine (:,1) of alpha
        real (kind=R8), allocatable  :: fcQbet(:,:)    ! Cosine (:,0) and sine (:,1) of beta
        real (kind=R8), allocatable  :: fcPbs(:)       ! Magnetic contact area
        real (kind=R8), allocatable  :: fcPbshz(:)     ! Magnetic contact area * hz (at faces)
        real (kind=R8), allocatable  :: fcBzb(:)       ! Factor for P-S/mdf flows
        real (kind=R8), allocatable  :: fcEb(:,:)      ! Unit vector pointing along magnetic field (Cartesian X,Y,Z coordinates)
        real (kind=R8), allocatable  :: fcFpsi(:)      ! Psi in each face center

!  vertex quantities
        real (kind=R8), allocatable  :: vxBb(:,:)      ! Magnetic field: (nVx,0:3) to be recomputed at vertices (Bx,By,Bz,|B|)
        real (kind=R8), allocatable  :: vxX(:)         ! X-co of vertex
        real (kind=R8), allocatable  :: vxY(:)         ! Y-co of vertex
        real (kind=R8), allocatable  :: vxHz(:)        ! Length in third (z) direction
        real (kind=R8), allocatable  :: vxVol(:)       ! Part of volume of surrounding cell centers associated with the vertex; for interpolation
        real (kind=R8), allocatable  :: vxFfbz(:)      !
        real (kind=R8), allocatable  :: vxFpsi(:)      !
        real (kind=R8), allocatable  :: vxOnedBsq(:)   ! 1/B**2 at vertices
        real (kind=R8), allocatable  :: vxBzb(:)       ! Factor for P-S/mdf flows
        real (kind=R8), allocatable  :: vxEb(:,:)      ! Unit vector pointing along magnetic field (Cartesian X,Y,Z coordinates)

!  flux tube/surface quantities
        real (kind=R8), allocatable  :: cvConn(:)      ! (nCv) Connection length for each CV
        real (kind=R8), allocatable  :: vxConn(:)      ! (nVx) Connection length for each vertex
        real (kind=R8), allocatable  :: ftConn(:)      ! (nFt) Connection length for each flux tube
        real (kind=R8), allocatable  :: fsConn(:)      ! (nFs) Connection length for each flux surface
        real (kind=R8), allocatable  :: ftEps(:)       ! (nFt) Inverse aspect ratio
        real (kind=R8), allocatable  :: ftBbav2(:)     ! (nFt) Average B**2 on the flux tube

        real (kind=R8), allocatable  :: ftFpsi(:)      ! (nFt) Psi of a flux tube; for a cell with cvFt(iCv)=iFt, we should have ftFpsi(iFt) = cvFpsi(cvFt(iCv))
        real (kind=R8), allocatable  :: fsPsi(:)       ! (nFs) Psi on each flux surface; for  face with fcFs(iFc)=iFs, we should have fsPsi(iFs) = fcFpsi(fcFs(iFc))

        real (kind=R8), allocatable  :: ftDs(:)        ! (nFt) A flux tube coordinate based on the radial distance at the OMP
        real (kind=R8), allocatable  :: fsDs(:)        ! (nFs) A flux surface coordinate based on the radial distance at the OMP

        logical :: LSN ! Identifies topologies where the western target is
                       ! the inner one
        logical :: psi_increasing ! .true. if the psi value increases from the O-point to the active X-point
      end type geometry

      contains

***********************************************************************

      subroutine alloc_geometry (g, nCv, nFc, nVx, nFs, nFt, nVmxCv)

      implicit none
      integer, intent(in) :: nCv, nFc, nVx, nFs, nFt, nVmxCv
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
      allocate (g%cvFpsi(nCv))

!  face quantites
      allocate (g%fcBb(nFc,0:3))
      allocate (g%fcS(nFc))
      allocate (g%fcHc(nFc,1:2))
      allocate (g%fcHt(nFc))
      allocate (g%fcHz(nFc))
      allocate (g%fcVol(nFc,1:2))
      allocate (g%fcQgam(nFc,0:1))
      allocate (g%fcQalf(nFc,0:1))
      allocate (g%fcQbet(nFc,0:1))
      allocate (g%fcPbs(nFc))
      allocate (g%fcPbshz(nFc))
      allocate (g%fcBzb(nFc))
      allocate (g%fcEb(nFc,0:2))
      allocate (g%fcFpsi(nFc))

!  vertex quantites
      allocate (g%vxBb(nVx,0:3))
      allocate (g%vxX(nVx))
      allocate (g%vxY(nVx))
      allocate (g%vxHz(nVx))
      allocate (g%vxVol(nVmxCv))
      allocate (g%vxFfbz(nVx))
      allocate (g%vxFpsi(nVx))
      allocate (g%vxOnedBsq(nVx))
      allocate (g%vxBzb(nVx))
      allocate (g%vxEb(nVx,0:2))

!  flux tube quantites
      allocate (g%cvConn(nCv))
      allocate (g%ftConn(nFt))
      allocate (g%ftEps(nFt))
      allocate (g%ftBbav2(nFt))

!  flux surface quantities
      allocate (g%vxConn(nVx))
      allocate (g%fsConn(nFs))
      allocate (g%fsPsi(nFs))
      allocate (g%ftFpsi(nFt))
      allocate (g%fsDs(nFs))
      allocate (g%ftDs(nFt))

      call initialize_geometry(g)

      return
      end subroutine alloc_geometry

***********************************************************************

      subroutine initialize_geometry (g)

      implicit none
      type(geometry), intent(in out) :: g

!  cell quantites
      g%cvBb = 0._R8
      g%cvX = 0._R8
      g%cvY = 0._R8
      g%cvSz = 0._R8
      g%cvHz = 0._R8
      g%cvHx = 0._R8
      g%cvQgam = 0._R8
      g%cvVol = 0._R8
      g%cvOnedBsq = 0._R8
      g%cvBzb = 0._R8
      g%cvEb = 0._R8
      g%cvFpsi = 0._R8

!  face quantites
      g%fcBb = 0._R8
      g%fcS = 0._R8
      g%fcHc = 0._R8
      g%fcHt = 0._R8
      g%fcHz = 0._R8
      g%fcVol = 0._R8
      g%fcQgam = 0._R8
      g%fcQalf = 0._R8
      g%fcQbet = 0._R8
      g%fcPbs = 0._R8
      g%fcPbshz = 0._R8
      g%fcBzb = 0._R8
      g%fcEb = 0._R8
      g%fcFpsi = 0._R8

!  vertex quantites
      g%vxBb = 0._R8
      g%vxX = 0._R8
      g%vxY = 0._R8
      g%vxHz = 0._R8
      g%vxVol = 0._R8
      g%vxFfbz = 0._R8
      g%vxFpsi = 0._R8
      g%vxOnedBsq = 0._R8
      g%vxBzb = 0._R8
      g%vxEb = 0._R8

!  flux tube quantites
      g%cvConn = 0._R8
      g%ftConn = 0._R8
      g%ftEps = 0._R8
      g%ftBbav2 = 0._R8

!  flux surface quantities
      g%vxConn = 0._R8
      g%fsConn = 0._R8
      g%fsPsi = 0._R8
      g%ftFpsi = 0._R8

      g%fsDs = 0._R8
      g%ftDs = 0._R8

      g%LSN = .true.
      g%psi_increasing = .true.

      return
      end subroutine initialize_geometry

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
      deallocate (g%cvFpsi)

!  face quantites
      deallocate (g%fcBb)
      deallocate (g%fcS)
      deallocate (g%fcHc)
      deallocate (g%fcHt)
      deallocate (g%fcHz)
      deallocate (g%fcVol)
      deallocate (g%fcQgam)
      deallocate (g%fcQalf)
      deallocate (g%fcQbet)
      deallocate (g%fcPbs)
      deallocate (g%fcPbshz)
      deallocate (g%fcBzb)
      deallocate (g%fcEb)
      deallocate (g%fcFpsi)

!  vertex quantites
      deallocate (g%vxBb)
      deallocate (g%vxX)
      deallocate (g%vxY)
      deallocate (g%vxHz)
      deallocate (g%vxVol)
      deallocate (g%vxFfbz)
      deallocate (g%vxFpsi)
      deallocate (g%vxOnedBsq)
      deallocate (g%vxBzb)
      deallocate (g%vxEb)

!  flux tube quantites
      deallocate (g%cvConn)
      deallocate (g%ftConn)
      deallocate (g%ftEps)
      deallocate (g%ftBbav2)

!  flux surface quantites
      deallocate (g%vxConn)
      deallocate (g%fsConn)
      deallocate (g%fsPsi)
      deallocate (g%ftFpsi)
      deallocate (g%fsDs)
      deallocate (g%ftDs)

      return
      end subroutine dealloc_geometry

***********************************************************************


      subroutine read_geometry (iun, nCv, nFc, nVx, nFs, nFt, nVmxCv,
     .  gm)

      implicit none
      integer, intent(in) :: iun, nCv, nFc, nVx, nFs, nFt, nVmxCv
      type (geometry), intent (inout) :: gm
      external cfrure

      call alloc_geometry (gm, nCv, nFc, nVx, nFs, nFt, nVmxCv)

! cell volumes
      call cfrure (iun, nCv*4, gm%cvBb,   'cvBb')
      call cfrure (iun, nCv*3, gm%cvEb,   'cvEb')
      call cfrure (iun, nCv,   gm%cvX,    'cvX')
      call cfrure (iun, nCv,   gm%cvY,    'cvY')
      call cfrure (iun, nCv,   gm%cvSz,   'cvSz')
      call cfrure (iun, nCv,   gm%cvHz,   'cvHz')
      call cfrure (iun, nCv,   gm%cvHx,   'cvHx') !WG temp!
      call cfrure (iun, nCv*2, gm%cvQgam, 'cvQgam')
      call cfrure (iun, nCv,   gm%cvVol,  'cvVol')

! face quantities
      call cfrure (iun, nFc*4, gm%fcBb,   'fcBb')
      call cfrure (iun, nFc,   gm%fcS,    'fcS')
      call cfrure (iun, nFc*2, gm%fcHc,   'fcHc')
      call cfrure (iun, nFc,   gm%fcHt,   'fcHt')
      call cfrure (iun, nFc*2, gm%fcQgam, 'fcQgam')
      call cfrure (iun, nFc*2, gm%fcQalf, 'fcQalf')
      call cfrure (iun, nFc*2, gm%fcQbet, 'fcQbet')
      call cfrure (iun, nFc,   gm%fcPbs,  'fcPbs')

! vertex quantities
      call cfrure (iun, nVx*4, gm%vxBb,   'vxBb')
      call cfrure (iun, nVx,   gm%vxX,    'vxX')
      call cfrure (iun, nVx,   gm%vxY,    'vxY')
      call cfrure (iun, nVx,   gm%vxFfbz, 'vxFfbz')
      call cfrure (iun, nVx,   gm%vxFpsi, 'vxFpsi')

! flux surface quantities
      call cfrure (iun, nCv,   gm%cvConn, 'cvConn')
      call cfrure (iun, nFs,   gm%fsPsi,  'fsPsi')

! check consistency of input
      call check_geometry(nCv, nFc, nVx, nFs, nFt, gm)

      return
      end subroutine read_geometry

***********************************************************************

      subroutine write_geometry (nout, nCv, nFc, nVx, nFs, nFt,
     .                           nVmxCv, gm)

      implicit none
      integer, intent(in) :: nout, nCv, nFc, nVx, nFs, nFt, nVmxCv
      type (geometry), intent (in) :: gm
      external cfwure

! cell volumes
      call cfwure (nout, nCv*4, gm%cvBb,   'cvBb')
      call cfwure (nout, nCv*3, gm%cvEb,   'cvEb')
      call cfwure (nout, nCv,   gm%cvX,    'cvX')
      call cfwure (nout, nCv,   gm%cvY,    'cvY')
      call cfwure (nout, nCv,   gm%cvSz,   'cvSz')
      call cfwure (nout, nCv,   gm%cvHz,   'cvHz')
      call cfwure (nout, nCv,   gm%cvHx,   'cvHx') !WG temp!
      call cfwure (nout, nCv*2, gm%cvQgam, 'cvQgam')
      call cfwure (nout, nCv,   gm%cvVol,  'cvVol')

! face quantities
      call cfwure (nout, nFc*4, gm%fcBb,   'fcBb')
      call cfwure (nout, nFc,   gm%fcS,    'fcS')
      call cfwure (nout, nFc*2, gm%fcHc,   'fcHc')
      call cfwure (nout, nFc,   gm%fcHt,   'fcHt')
      call cfwure (nout, nFc*2, gm%fcQgam, 'fcQgam')
      call cfwure (nout, nFc*2, gm%fcQalf, 'fcQalf')
      call cfwure (nout, nFc*2, gm%fcQbet, 'fcQbet')
      call cfwure (nout, nFc,   gm%fcPbs,  'fcPbs')

! vertex quantities
      call cfwure (nout, nVx*4, gm%vxBb,   'vxBb')
      call cfwure (nout, nVx,   gm%vxX,    'vxX')
      call cfwure (nout, nVx,   gm%vxY,    'vxY')
      call cfwure (nout, nVx,   gm%vxFfbz, 'vxFfbz')
      call cfwure (nout, nVx,   gm%vxFpsi, 'vxFpsi')

! flux surface quantities
      call cfwure (nout, nCv,   gm%cvConn, 'cvConn')
      call cfwure (nout, nFs,   gm%fsPsi,  'fsPsi')

      return
      end subroutine write_geometry

***********************************************************************

      subroutine check_geometry (nCv, nFc, nVx, nFs, nFt, gm)
      implicit none
      integer, intent(in) :: nCv, nFc, nVx, nFs, nFt
      type (geometry), intent (in) :: gm
      external b2xvsg, xertst

      call xertst (nCv.gt.0, 'faulty argument nCv !')
      call xertst (nFc.gt.0, 'faulty argument nFc !')
      call xertst (nVx.gt.0, 'faulty argument nVx !')
      call xertst (nFs.gt.0, 'faulty argument nFs !')
      call xertst (nFt.gt.0, 'faulty argument nFt !')

*   ..check sign of magnetic field components
      call b2xvsg (nCv, gm % cvBb(:,0), 1, 'cvBb0', '.ge.')
      call b2xvsg (nCv, gm % cvBb(:,1), 1, 'cvBb1', '.ge.')
      call b2xvsg (nCv, gm % cvBb(:,2), 1, 'cvBb2', '.ge.')
      call b2xvsg (nCv, gm % cvBb(:,3), 1, 'cvBb3', '.ge.')

      call b2xvsg (nFc, gm % fcBb(:,0), 1, 'fcBb0', '.ge.')
      call b2xvsg (nFc, gm % fcBb(:,1), 1, 'fcBb1', '.ge.')
      call b2xvsg (nFc, gm % fcBb(:,2), 1, 'fcBb2', '.ge.')
      call b2xvsg (nFc, gm % fcBb(:,3), 1, 'fcBb3', '.ge.')

      call b2xvsg (nVx, gm % vxBb(:,0), 1, 'vxBb0', '.ge.')
      call b2xvsg (nVx, gm % vxBb(:,1), 1, 'vxBb1', '.ge.')
      call b2xvsg (nVx, gm % vxBb(:,2), 1, 'vxBb2', '.ge.')
      call b2xvsg (nVx, gm % vxBb(:,3), 1, 'vxBb3', '.ge.')

      return
      end subroutine check_geometry

***********************************************************************

      subroutine init_geometry (mpg, gm, switch)
      use b2mod_constants
      use b2mod_b2cmfs
     &  , only : isymm
      use b2mod_switches
      use b2us_map
      implicit none
      type (switches), intent (in) :: switch
      type (mapping), intent (inout) :: mpg
      type (geometry), intent (inout) :: gm
      integer :: i, j, k, l, nCv, nFc, nVx, nFx,
     &           iCv, iFc, iVx, inCv, iXpt, iFc1, iFc2, iVx1, iVx2,
     &           i1, i2, inv_dist(mpg%nVx) , vv(2)
      integer, allocatable :: old_face_list(:), verts(:)
      real (kind=R8) :: hzconst, dux, duy, du, sbf, Psi1, Psi2,
     .  count_up, count_down, count_eq, dpsi(mpg%nCv),
     .  dpsi_max_down, dpsi_max_up
      logical active, match_found
      intrinsic maxval, maxloc, minloc, abs, sqrt
      external b2xbzb, intface, intvertex, xerrab, xertst

      nCv = mpg%nCv
      nFc = mpg%nFc
      nVx = mpg%nVx

*   ..identify periodicity
      call is_periodic( mpg, gm )

*   ..identify cells on closed surfaces
      if (mpg%nnreg(0).eq.1) then
        mpg%cvOnClosedSurface = mpg%periodic_bc.eq.1
      else if (mpg%nnreg(0).eq.2 .or. mpg%nnreg(0).eq.4 .or.
     .         mpg%nnreg(0).eq.5 .or. mpg%nnreg(0).eq.8) then
        mpg%cvOnClosedSurface = mod(mpg%cvReg,4).eq.1
      else if (mpg%nnreg(0).eq.7) then
        mpg%cvOnClosedSurface = mpg%cvReg.eq.1
      else
        mpg%cvOnClosedSurface = mpg%cvReg.eq.1
        write (*,*) 'Warning: unknown number of regions'
        write (*,*) 'Assuming core is region 1'
        !call xerrab ('Unknown number of regions!')
      end if

*   ..set weights for interpolation to faces
      if (switch%fcVol_style.eq.0) then
       ! weights for inverse volume interpolation
       do iFc = 1, nFc
         gm%fcVol(iFc,1) = 0.5_R8*gm%cvVol(mpg%fcCv(iFc,1))
         gm%fcVol(iFc,2) = 0.5_R8*gm%cvVol(mpg%fcCv(iFc,2))
       enddo
      elseif (switch%fcVol_style.eq.1) then
       ! weights for simple averaging
       do iFc = 1, nFc
         gm%fcVol(iFc,1) = 0.5_R8
         gm%fcVol(iFc,2) = 0.5_R8
       enddo
      else
        call xerrab ('wrong value of b2mndr_fcVol_style')
      endif

*   ..precompute psi in each cell according to vxVol_style=2
      do iCv = 1, nCv
        allocate(verts(1:mpg%cvVxP(iCv,2)))
        verts = mpg%cvVx(mpg%cvVxP(iCv,1):mpg%cvVxP(iCv,1)+
     &          mpg%cvVxP(iCv,2)-1)
        gm%cvFpsi(iCv) = 0.5_R8*(
     &    maxval(gm%vxFpsi(verts(1:mpg%cvVxP(iCv,2)))) +
     &    minval(gm%vxFpsi(verts(1:mpg%cvVxP(iCv,2)))) )
        dpsi(iCv) = abs(
     &    maxval(gm%vxFpsi(verts(1:mpg%cvVxP(iCv,2)))) -
     &    minval(gm%vxFpsi(verts(1:mpg%cvVxP(iCv,2)))) )
       deallocate(verts)
      enddo

*   ..precompute psi in each face center
      do iFc = 1, nFc
        gm%fcFpsi(iFc) = 0.5_R8*(gm%vxFpsi(mpg%fcVx(iFc,1)) +
     &    gm%vxFpsi(mpg%fcVx(iFc,2)) )
      enddo

*   ..precompute psi for flux tube
      do iFt = 1, nFt
        gm%ftFpsi(iFt) = gm%cvFpsi(mpg%ftCv(mpg%ftCvP(iFt,1)))
        do i = 1, mpg%ftCvP(iFt,2)
          iCv = mpg%ftCv(mpg%ftCvP(iFt,1)+i-1)
          if (gm%cvFpsi(iCv).ne.gm%ftFpsi(iFt))
     &      call xerrab('Inconsistent definition of cvFpsi')
        enddo
      enddo

*   ..set weights for interpolation to vertices
      if (switch%vxVol_style.eq.0) then
       ! weights for inverse volume interpolation
       ! idea: each cell divides its volume equally
       !       between its vertices
       do iVx = 1, nVx
        do inCv = 1, mpg%vxCvP(iVx,2)
         iCv = mpg%vxCv(mpg%vxCvP(iVx,1) + inCv - 1)
         gm%vxVol(mpg%vxCvP(iVx,1) + inCv - 1) =
     &     (1._R8/mpg%cvVxP(iCv,2))*gm%cvVol(iCv)
        enddo
       enddo
      elseif (switch%vxVol_style.eq.1) then
       ! weights based on relative volume fraction
       ! idea: each cell assigns equal volume fraction
       !       to each of its vertices
       do iVx = 1, nVx
        do inCv = 1, mpg%vxCvP(iVx,2)
         iCv = mpg%vxCv(mpg%vxCvP(iVx,1) + inCv - 1)
         if (iCv.le.mpg%nCi) then
          gm%vxVol(mpg%vxCvP(iVx,1) + inCv - 1) = 1._R8/mpg%cvVxP(iCv,2)
         else
          ! Make sure guard cells dominate in the interpolation
          ! to vertices => use volume
          gm%vxVol(mpg%vxCvP(iVx,1) + inCv - 1) =
     &      (1._R8/mpg%cvVxP(iCv,2))*gm%cvVol(iCv)
         endif
        enddo
       enddo
      elseif (switch%vxVol_style.eq.2) then
       ! Assumption: all internal 'radial' faces aligned (boundaries faces possibly non-aligned)
       ! Idea: assign equal total weight to all cells on either side of the flux surface (inverse volume weighting!)
       !       (expected: in case of zero gradient along poloidal direction, value at vertex will be
       !                  average between flux tube values)
       ! Count cells per vertex with higher, lower Psi
       ! Exceptional case: equal Psi (pentagon splitting);
       ! will also catch case where Psi = 0 everywhere
       do iVx = 1, nVx
         count_down    = 0.0_R8
         count_up      = 0.0_R8
         count_eq      = 0.0_R8
         do inCv = 1, mpg%vxCvP(iVx,2)
           iCv = mpg%vxCv(mpg%vxCvP(iVx,1) + inCv - 1)
           if (gm%cvFpsi(iCv) .gt. gm%vxFpsi(iVx)) then
             count_up = count_up + 1.0e0_R8
           elseif (gm%cvFpsi(iCv) .lt. gm%vxFpsi(iVx)) then
             count_down = count_down + 1.0e0_R8
           else
             count_eq = count_eq + 1.0e0_R8
           endif 
         enddo
         do inCv = 1, mpg%vxCvP(iVx,2)
           iCv = mpg%vxCv(mpg%vxCvP(iVx,1) + inCv - 1)
           if (iCv.le.mpg%nCi) then
            if ( gm%cvFpsi(iCv) .gt. gm%vxFpsi(iVx)) then
              gm%vxVol(mpg%vxCvP(iVx,1) + inCv - 1) = 
     &          2.0_R8*count_up
            elseif ( gm%cvFpsi(iCv) .lt. gm%vxFpsi(iVx)) then
              gm%vxVol(mpg%vxCvP(iVx,1) + inCv - 1) = 
     &          2.0_R8*count_down
            else
              gm%vxVol(mpg%vxCvP(iVx,1) + inCv - 1) = 
     &          count_eq
            endif
           else
            ! Make sure guard cells dominate in the interpolation
            ! to vertices => use small number
            gm%vxVol(mpg%vxCvP(iVx,1) + inCv - 1) = 0.25e-3_R8          
           endif 
         enddo
       enddo
      elseif (switch%vxVol_style.eq.3) then
       ! equal weight for each cell around the vertex
       gm%vxVol = 1.0_R8
      elseif (switch%vxVol_style.eq.4) then
       ! try-out: 
       ! assign 1/4, except guard cells
       do iVx = 1, nVx
        do inCv = 1, mpg%vxCvP(iVx,2)
         iCv = mpg%vxCv(mpg%vxCvP(iVx,1) + inCv - 1)
         if (iCv.le.mpg%nCi) then
          gm%vxVol(mpg%vxCvP(iVx,1) + inCv - 1) = 0.25_R8
         else
          ! Make sure guard cells dominate in the interpolation
          ! to vertices => use volume
          gm%vxVol(mpg%vxCvP(iVx,1) + inCv - 1) = 0.25e-3_R8
         endif
        enddo
        if (mpg%vxCvP(iVx,2).eq.3) then
          do inCv = 1, mpg%vxCvP(iVx,2)
           iCv = mpg%vxCv(mpg%vxCvP(iVx,1) + inCv - 1)
           if (iCv.le.mpg%nCi) then
            if (mpg%cvVxP(iCv,2) .eq. 5) then
             gm%vxVol(mpg%vxCvP(iVx,1) + inCv - 1) = 0.125_R8
            else
             gm%vxVol(mpg%vxCvP(iVx,1) + inCv - 1) = 0.25_R8            
            endif
           else
            ! Make sure guard cells dominate in the interpolation
            ! to vertices => use volume
            gm%vxVol(mpg%vxCvP(iVx,1) + inCv - 1) = 0.25e-3_R8
           endif
        enddo        
        endif
       enddo
      elseif (switch%vxVol_style.eq.5) then     
       do iVx = 1, nVx
        do inCv = 1, mpg%vxCvP(iVx,2)
         iCv = mpg%vxCv(mpg%vxCvP(iVx,1) + inCv - 1)
         gm%vxVol(mpg%vxCvP(iVx,1) + inCv - 1) = 
     &         1.0_R8/gm%cvVol(iCv)
        enddo
       enddo   
      elseif (switch%vxVol_style.eq.6) then 
       ! try-out: 
       ! assign 1/4, except guard cells
       ! 1/8 for pent at common vertex
       ! 1/2 for pent at flux surface verts
       do iVx = 1, nVx

        do inCv = 1, mpg%vxCvP(iVx,2)
         iCv = mpg%vxCv(mpg%vxCvP(iVx,1) + inCv - 1)
         if (iCv.le.mpg%nCi) then
          gm%vxVol(mpg%vxCvP(iVx,1) + inCv - 1) = 0.25_R8
         else
          ! Make sure guard cells dominate in the interpolation
          ! to vertices => use volume
          gm%vxVol(mpg%vxCvP(iVx,1) + inCv - 1) = 0.25e-3_R8
         endif
        enddo

        if (mpg%vxCvP(iVx,2).eq.3) then
          do inCv = 1, mpg%vxCvP(iVx,2)
           iCv = mpg%vxCv(mpg%vxCvP(iVx,1) + inCv - 1)
           if (iCv.le.mpg%nCi) then
            if (mpg%cvVxP(iCv,2) .eq. 5) then
             gm%vxVol(mpg%vxCvP(iVx,1) + inCv - 1) = 0.125_R8
            else
             gm%vxVol(mpg%vxCvP(iVx,1) + inCv - 1) = 0.25_R8            
            endif
           else
            ! Make sure guard cells dominate in the interpolation
            ! to vertices => use volume
            gm%vxVol(mpg%vxCvP(iVx,1) + inCv - 1) = 0.25e-3_R8
           endif
          enddo

        elseif (mpg%vxCvP(iVx,2).eq.4) then
         do inCv = 1, mpg%vxCvP(iVx,2)   
           iCv = mpg%vxCv(mpg%vxCvP(iVx,1) + inCv - 1)
           if (iCv.le.mpg%nCi) then
             gm%vxVol(mpg%vxCvP(iVx,1) + inCv - 1) = 
     &            gm%cvVol(iCv)
           else
            ! Make sure guard cells dominate in the interpolation
            ! to vertices => use volume
            gm%vxVol(mpg%vxCvP(iVx,1) + inCv - 1) = 
     &           gm%cvVol(iCv) / 1e2_R8       
           endif
         enddo
                       
        endif

       enddo

      elseif (switch%vxVol_style .eq. 7) then 
         !recalculate cvFpsi
         do iCv = 1, nCv
            allocate(verts(1:mpg%cvVxP(iCv,2)))
            verts = mpg%cvVx(mpg%cvVxP(iCv,1):mpg%cvVxP(iCv,1)+
     &         mpg%cvVxP(iCv,2))
            gm%cvFpsi(iCv) = sum(gm%vxFpsi(verts))/mpg%cvVxP(iCv,2)
            i1 = minloc(gm%vxFpsi(verts),1) -1
            i2 = maxloc(gm%vxFpsi(verts),1) -1
            dpsi(iCv) = abs(gm%vxFpsi(verts(i2)) - 
     1                gm%vxFpsi(verts(i1)) )
            deallocate(verts)
         enddo
         !count cells per vertex with higher and lower Psi
         do iVx = 1, nVx
            count_down    = 0.0_R8
            count_up      = 0.0_R8
            do inCv = 1, mpg%vxCvP(iVx,2)
               iCv = mpg%vxCv(mpg%vxCvP(iVx,1)+ inCv - 1)
               if ( gm%cvFpsi(iCv) .gt. gm%vxFpsi(iVx)) then
                  count_up = count_up + 1.0e0_R8
               elseif ( gm%cvFpsi(iCv) .lt. gm%vxFpsi(iVx)) then
                  count_down = count_down + 1.0e0_R8
               endif 
            enddo
            do inCv = 1, mpg%vxCvP(iVx,2)
              iCv = mpg%vxCv(mpg%vxCvP(iVx,1)+ inCv - 1)
              if (iCv.le.mpg%nCi) then
               if ( gm%cvFpsi(iCv) .gt. gm%vxFpsi(iVx)) then
                   gm%vxVol(mpg%vxCvP(iVx,1)+inCv - 1) = 
     &               count_up / mpg%vxCvP(iVx,2)               
               elseif ( gm%cvFpsi(iCv) .lt. gm%vxFpsi(iVx)) then
                    gm%vxVol(mpg%vxCvP(iVx,1)+inCv - 1) = 
     &                count_down / mpg%vxCvP(iVx,2)                  
               endif
              else 
                 ! Make guards cells dominate and account for 
                 ! length difference via distance to vertex
                 gm%vxVol(mpg%vxCvP(iVx,1)+inCv - 1) = 0.25e-3_R8 *
     &             sqrt( (gm%cvX(iCv) - gm%vxX(iVx))**2  + 
     &                   (gm%cvY(iCv) - gm%vxY(iVx))**2  )            
              endif 
            enddo
            !Special case for pentagon splitting
            if (mpg%vxCvP(iVx,2).eq.3) then
              do inCv = 1, mpg%vxCvP(iVx,2)
               iCv = mpg%vxCv(mpg%vxCvP(iVx,1)+ inCv - 1)
               if (iCv .le. mpg%nCi) then
                  if (mpg%cvVxP(iCv,2).eq.5) then
                    gm%vxVol(mpg%vxCvP(iVx,1)+inCv - 1) = 0.125_R8
                  else
                    gm%vxVol(mpg%vxCvP(iVx,1)+inCv - 1) = 0.25_R8  
                  endif 
               else
                 gm%vxVol(mpg%vxCvP(iVx,1)+inCv - 1) = 0.25e-3_R8 *
     &             sqrt( (gm%cvX(iCv) - gm%vxX(iVx))**2  + 
     &                   (gm%cvY(iCv) - gm%vxY(iVx))**2  )                 
               endif

              enddo

            endif

         enddo
      elseif (switch%vxVol_style .eq. 8) then
      !Same as style 2 but with weighting for boundary vertices
       ! assign 1/4, except guard cells
       do iVx = 1, nVx
        do inCv = 1, mpg%vxCvP(iVx,2)
         iCv = mpg%vxCv(mpg%vxCvP(iVx,1) + inCv - 1)
         if (iCv.le.mpg%nCi) then
          gm%vxVol(mpg%vxCvP(iVx,1) + inCv - 1) = 0.25_R8
         else
          ! Make sure guard cells dominate in the interpolation
          ! to vertices => use volume
          gm%vxVol(mpg%vxCvP(iVx,1) + inCv - 1) = 0.25e-3_R8 *
     &             sqrt( (gm%cvX(iCv) - gm%vxX(iVx))**2  + 
     &                   (gm%cvY(iCv) - gm%vxY(iVx))**2  ) 
         endif
        enddo
       enddo 

      elseif (switch%vxVol_style .eq. 9) then 
         !recalculate cvFpsi
         do iCv = 1, nCv
            allocate(verts(1:mpg%cvVxP(iCv,2)))
            verts = mpg%cvVx(mpg%cvVxP(iCv,1):mpg%cvVxP(iCv,1)+
     &         mpg%cvVxP(iCv,2))
            gm%cvFpsi(iCv) = sum(gm%vxFpsi(verts))/mpg%cvVxP(iCv,2)
            i1 = minloc(gm%vxFpsi(verts),1) -1
            i2 = maxloc(gm%vxFpsi(verts),1) -1
            dpsi(iCv) = abs(gm%vxFpsi(verts(i2)) - 
     1                gm%vxFpsi(verts(i1)) )
            deallocate(verts)
         enddo
         !count cells per vertex with higher and lower Psi
         do iVx = 1, nVx
            count_down    = 0.0_R8
            count_up      = 0.0_R8
            dpsi_max_down = 0.0_R8
            dpsi_max_up   = 0.0_R8
            do inCv = 1, mpg%vxCvP(iVx,2)
               iCv = mpg%vxCv(mpg%vxCvP(iVx,1)+ inCv - 1)
               if ( gm%cvFpsi(iCv) .gt. gm%vxFpsi(iVx)) then
                  count_up = count_up + 1.0e0_R8
                  if ( dpsi(iCv) .gt. dpsi_max_up) then
                     dpsi_max_up = dpsi(iCv)
                  endif
               elseif ( gm%cvFpsi(iCv) .lt. gm%vxFpsi(iVx)) then
                  count_down = count_down + 1.0e0_R8
                  if ( dpsi(iCv) .gt. dpsi_max_down) then
                     dpsi_max_down = dpsi(iCv)
                  endif
               endif 
            enddo
            do inCv = 1, mpg%vxCvP(iVx,2)
              iCv = mpg%vxCv(mpg%vxCvP(iVx,1)+ inCv - 1)
              if (iCv.le.mpg%nCi) then
               if ( gm%cvFpsi(iCv) .gt. gm%vxFpsi(iVx)) then
                 if (mpg%cvVxP(iCv,2) .eq. 3) then
                   gm%vxVol(mpg%vxCvP(iVx,1)+inCv - 1) = 
     &               count_up / (mpg%vxCvP(iVx,2)*sqrt(2.0_R8)* 
     &                sqrt(dpsi_max_up/dpsi(iCv)))               
                 else
                   gm%vxVol(mpg%vxCvP(iVx,1)+inCv - 1) = 
     &               count_up / mpg%vxCvP(iVx,2)
                 endif 
               elseif ( gm%cvFpsi(iCv) .lt. gm%vxFpsi(iVx)) then
                  if (mpg%cvVxP(iCv,2) .eq. 3) then
                    gm%vxVol(mpg%vxCvP(iVx,1)+inCv - 1) = 
     &                count_down / (mpg%vxCvP(iVx,2)*sqrt(2.0_R8)*
     &                   sqrt(dpsi_max_up/dpsi(iCv)))                  
                  else
                    gm%vxVol(mpg%vxCvP(iVx,1)+inCv - 1) = 
     &                count_down / mpg%vxCvP(iVx,2)
                  endif
               endif
              else 
                 ! Make guards cells dominate and account for 
                 ! length difference via distance to vertex
                 gm%vxVol(mpg%vxCvP(iVx,1)+inCv - 1) = 0.25e-3_R8 *
     &             sqrt( (gm%cvX(iCv) - gm%vxX(iVx))**2  + 
     &                   (gm%cvY(iCv) - gm%vxY(iVx))**2  )            
              endif 
            enddo
            !Special case for pentagon splitting
            if (mpg%vxCvP(iVx,2).eq.3) then
              do inCv = 1, mpg%vxCvP(iVx,2)
               iCv = mpg%vxCv(mpg%vxCvP(iVx,1)+ inCv - 1)
               if (iCv .le. mpg%nCi) then
                  if (mpg%cvVxP(iCv,2).eq.5) then
                    gm%vxVol(mpg%vxCvP(iVx,1)+inCv - 1) = 0.125_R8
                  else
                    gm%vxVol(mpg%vxCvP(iVx,1)+inCv - 1) = 0.25_R8  
                  endif 
               else
                 gm%vxVol(mpg%vxCvP(iVx,1)+inCv - 1) = 0.25e-3_R8 *
     &             sqrt( (gm%cvX(iCv) - gm%vxX(iVx))**2  + 
     &                   (gm%cvY(iCv) - gm%vxY(iVx))**2  )                 
               endif

              enddo

            endif

         enddo

      elseif (switch%vxVol_style .eq. 10) then
       do iVx = 1, nVx
        do inCv = 1, mpg%vxCvP(iVx,2)
         iCv = mpg%vxCv(mpg%vxCvP(iVx,1) + inCv - 1)
         if (iCv.le.mpg%nCi) then
          gm%vxVol(mpg%vxCvP(iVx,1) + inCv - 1) = 
     &             sqrt( (gm%cvX(iCv) - gm%vxX(iVx))**2  + 
     &                   (gm%cvY(iCv) - gm%vxY(iVx))**2  )     
         else
          ! Make sure guard cells dominate in the interpolation
          ! to vertices => use volume
          gm%vxVol(mpg%vxCvP(iVx,1) + inCv - 1) = 0.25e-3_R8 *
     &             sqrt( (gm%cvX(iCv) - gm%vxX(iVx))**2  + 
     &                   (gm%cvY(iCv) - gm%vxY(iVx))**2  ) 
         endif
        enddo
       enddo
      elseif (switch%vxVol_style .eq. 11) then
       ! Check which type to interpolation type to use 
       inv_dist(:) = 0
       do iVx = 1,nVx
        do inCv = 1, mpg%vxCvP(iVx,2)
         iCv = mpg%vxCv(mpg%vxCvP(iVx,1) + inCv - 1)
         if ((mpg%cvVxP(iCv,2).eq.3) .or.
     &      (mpg%cvVxP(iCv,2).gt.4)) then
           inv_dist(iVx) = 1
         endif
        enddo         
       enddo
       
       ! Compute the interpolation weights
       do iVx = 1,nVx
        if (inv_dist(iVx).eq.1) then
         do inCv = 1, mpg%vxCvP(iVx,2)
          iCv = mpg%vxCv(mpg%vxCvP(iVx,1) + inCv - 1)
           gm%vxVol(mpg%vxCvP(iVx,1) + inCv - 1) = 
     &             sqrt( (gm%cvX(iCv) - gm%vxX(iVx))**2  + 
     &                   (gm%cvY(iCv) - gm%vxY(iVx))**2  )     
          
         enddo   
        else
         do inCv = 1, mpg%vxCvP(iVx,2)
          iCv = mpg%vxCv(mpg%vxCvP(iVx,1) + inCv - 1)
          if (iCv.le.mpg%nCi) then
           gm%vxVol(mpg%vxCvP(iVx,1) + inCv - 1) = 0.25_R8
          else
          ! Make sure guard cells dominate in the interpolation
          ! to vertices => use volume
           gm%vxVol(mpg%vxCvP(iVx,1) + inCv - 1) = 0.25e-3_R8
          endif
         enddo         
        endif    
       enddo
     
      else
        call xerrab ('wrong value of b2mndr_vxVol_style')
      endif

*  .. calculate 1/B/B
      gm%cvOnedBsq = 1.0_R8/gm%cvBb(:,3)**2
      gm%vxOnedBsq = 1.0_R8/(gm%vxBb(:,3)**2 + 1.0e-30_R8)

*  .. unit vector of magnetic field in X,Y,Z coordinates at faces
      do iFc = 1, mpg%nFc

        ! Construct unit vector normal to face
        dux =   gm%vxY(mpg%fcVx(iFc,2))-gm%vxY(mpg%fcVx(iFc,1))
        duy = -(gm%vxX(mpg%fcVx(iFc,2))-gm%vxX(mpg%fcVx(iFc,1)))

        ! Normalize
        du  = sqrt(dux**2+duy**2)
        dux = dux/du
        duy = duy/du

        ! Sign of toroidal field (1: out of page; -1: into page)
        sbf = sign(1.0_R8,gm%fcQalf(iFc,0)*gm%fcQbet(iFc,0) -
     &                    gm%fcQalf(iFc,1)*gm%fcQbet(iFc,1))

        ! Rotate over angle -alpha; scale with poloidal pitch
        gm%fcEb(iFc,0) = gm%fcBb(iFc,0)/gm%fcBb(iFc,3)*
     &    (     gm%fcQalf(iFc,0)*dux + sbf*gm%fcQalf(iFc,1)*duy)
        gm%fcEb(iFc,1) = gm%fcBb(iFc,0)/gm%fcBb(iFc,3)*
     &    (-sbf*gm%fcQalf(iFc,1)*dux +     gm%fcQalf(iFc,0)*duy)
        gm%fcEb(iFc,2) = sbf*gm%fcBb(iFc,2)/gm%fcBb(iFc,3)

      end do

      ! normalize
      gm%fcEb(:,0) = gm%fcEb(:,0)/
     &  sqrt(gm%fcEb(:,0)**2 + gm%fcEb(:,1)**2 + gm%fcEb(:,2)**2)
      gm%fcEb(:,1) = gm%fcEb(:,1)/
     &  sqrt(gm%fcEb(:,0)**2 + gm%fcEb(:,1)**2 + gm%fcEb(:,2)**2)
      gm%fcEb(:,2) = gm%fcEb(:,2)/
     &  sqrt(gm%fcEb(:,0)**2 + gm%fcEb(:,1)**2 + gm%fcEb(:,2)**2)

!WG_TODO: Need to provide a calculation for gm%vxEb

*   ..compute cvHz
      hzconst = 2.0_R8*pi*6.0_R8   ! for ITER for testing
      hzconst = 2.0_R8*pi*1.6_R8   ! for ASDEX-UPGRADE for testing
      hzconst = 1.0_R8
      gm%cvHz = hzconst*(1.0e0_R8-switch%b2mndr_hz)+switch%b2mndr_hz*
     & gm%cvVol/gm%cvSz

*   ..compute fcHz and fcPbshz
      if (switch%fcHz_style.eq.0) then
        call intface(nCv,nFc,mpg%fcCv,gm%fcVol,gm%cvHz,gm%fcHz)
      else
        if (isymm.eq.0) then
          gm%fcHz = hzconst*(1.0e0_R8-switch%b2mndr_hz)+
     &      switch%b2mndr_hz*switch%b2agmt_width
        elseif (isymm.eq.1 .or. isymm.eq.2) then
          do iFc = 1, mpg%nFc
            gm%fcHz(iFc) = hzconst*(1.0e0_R8-switch%b2mndr_hz)+
     &        switch%b2mndr_hz*
     &        pi*(gm%vxX(mpg%fcVx(iFc,1))+gm%vxX(mpg%fcVx(iFc,2)))
          enddo
        else if (isymm.eq.3 .or. isymm.eq.4) then
          do iFc = 1, mpg%nFc
            gm%fcHz(iFc) = hzconst*(1.0e0_R8-switch%b2mndr_hz)+
     &        switch%b2mndr_hz*
     &        pi*(gm%vxY(mpg%fcVx(iFc,1))+gm%vxY(mpg%fcVx(iFc,2)))
          enddo
        endif
      endif
      gm%fcPbshz = gm%fcPbs * gm%fcHz

*   ..compute vxHz
      if (switch%vxHz_style.eq.0) then
        call intvertex(nCv,nVx,mpg,gm%vxVol,gm%cvHz,gm%vxHz)
      else
        if (isymm.eq.0) then
          gm%vxHz = hzconst*(1.0e0_R8-switch%b2mndr_hz)+
     &      switch%b2mndr_hz*switch%b2agmt_width
        elseif (isymm.eq.1 .or. isymm.eq.2) then
          gm%vxHz = hzconst*(1.0e0_R8-switch%b2mndr_hz)+
     &      switch%b2mndr_hz*2.0_R8*pi*gm%vxX
        else if (isymm.eq.3 .or. isymm.eq.4) then
          gm%vxHz = hzconst*(1.0e0_R8-switch%b2mndr_hz)+
     &      switch%b2mndr_hz*2.0_R8*pi*gm%vxY
        endif
      endif

*   ..compute Bz/B**2 - Bz/<B**2>
      call b2xbzb (nCv, nFc, nVx, mpg, gm%cvBb, gm%fcBb, gm%vxBb,
     &  gm%cvVol, gm%cvBzb, gm%fcBzb, gm%vxBzb)

*   ..modify geometry to recover 5-point stencil behavior
      ! only supported for cases derived from original, structured grid
      if (switch%use_9pt_stencil.eq.0) then
        if (mpg%isClassicalGrid.eq.1) then
          nfx = maxval(mpg%imapFcx)
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

*   ..calculate connection length
      call calc_conn(mpg, gm)

*   ..determine topology and sign of psi derivative
!! LSN = .false. if the tangency or X-point is above or
!!               at the same level and to the left of the O-point
      iXpt = 0
      if (mpg%nnreg(0).eq.1) then
!! TODO: Find a way to ascertain, if flux surfaces are defined,
!!       whether the psi value is increasing or decreasing radially
      elseif (mpg%nnreg(0).eq.2) then
        do i = mpg%vxFcP(mpg%tgVx(1),1),
     &         mpg%vxFcP(mpg%tgVx(1),1)+mpg%vxFcP(mpg%tgVx(1),2)-1
          iFc = mpg%vxFc(i)
          if ( mpg%cvReg(mpg%fcCv(iFc,1)).eq.1 .and.
     &         mpg%cvReg(mpg%fcCv(iFc,2)).eq.1) then
            if ( mpg%fcVx(iFc,1).eq.mpg%tgVx(1)) then
              iVx = mpg%fcVx(iFc,2)
            else
              iVx = mpg%fcVx(iFc,1)
            end if
            if ( gm%vxY(mpg%tgVx(1)).gt.gm%vxY(iVx) .or.
     &          (gm%vxY(mpg%tgVx(1)).eq.gm%vxY(iVx) .and.
     &           gm%vxX(mpg%tgVx(1)).gt.gm%vxX(iVx)) ) then
              gm%LSN = .false.
            end if
            gm%psi_increasing = gm%fsPsi(mpg%vxFs(iVx)).lt.
     &                          gm%fsPsi(mpg%vxFs(mpg%tgVx(1))).or.
     &                          maxval(abs(gm%fsPsi)).eq.0.0_R8
          end if
        end do
      else if (mpg%nnreg(0).eq.4) then
        do i = mpg%vxFcP(mpg%Xpt(1),1),
     &         mpg%vxFcP(mpg%Xpt(1),1)+mpg%vxFcP(mpg%Xpt(1),2)-1
          iFc = mpg%vxFc(i)
          if ( mpg%cvReg(mpg%fcCv(iFc,1)).eq.1 .and.
     &         mpg%cvReg(mpg%fcCv(iFc,2)).eq.1) then
            if ( mpg%fcVx(iFc,1).eq.mpg%Xpt(1)) then
              iVx = mpg%fcVx(iFc,2)
            else
              iVx = mpg%fcVx(iFc,1)
            end if
            if ( gm%vxY(mpg%Xpt(1)).gt.gm%vxY(iVx) .or.
     &          (gm%vxY(mpg%Xpt(1)).eq.gm%vxY(iVx) .and.
     &           gm%vxX(mpg%Xpt(1)).gt.gm%vxX(iVx)) ) then
              gm%LSN = .false.
            end if
          end if
        end do
        iXpt = 1
      else if (mpg%nnreg(0).eq.5) then
        iXpt = 1
      else if (mpg%nnreg(0).ge.7) then
        active = .false.
        do i = mpg%vxCvP(mpg%Xpt(1),1), mpg%vxCvP(mpg%Xpt(1),1) +
     &                                  mpg%vxCvP(mpg%Xpt(1),2) - 1
          iCv = mpg%vxCv(i)
          if (mpg%cvReg(iCv).eq.1) active = .true.
        end do
        if (active) then
          iXpt = 1
        else
          iXpt = 2
        end if
      end if
      if (iXpt.ne.0) then
       match_found = .false.
       do i = mpg%vxCvP(mpg%Xpt(iXpt),1), mpg%vxCvP(mpg%Xpt(iXpt),1) +
     &                                    mpg%vxCvP(mpg%Xpt(iXpt),2) - 1
        if (match_found) cycle
        iCv = mpg%vxCv(i)
        if (mpg%cvReg(iCv).eq.1) then
          do j = mpg%cvVxP(iCv,1),
     &           mpg%cvVxP(iCv,1) + mpg%cvVxP(iCv,2) - 1
            if (match_found) cycle
            if (mpg%vxFs(mpg%cvVx(j)).eq.mpg%vxFs(mpg%Xpt(iXpt))) cycle
            gm%psi_increasing = gm%fsPsi(mpg%vxFs(mpg%cvVx(j))).lt.
     &                          gm%fsPsi(mpg%vxFs(mpg%Xpt(iXpt))).or.
     &                          maxval(abs(gm%fsPsi)).eq.0.0_R8
            match_found = .true.
          end do
        end if
       end do
      end if

!! Invert the divertor face list if not in correct psi order
!! and identify strike point face indices
      do i = 1, maxval(mpg%strDiv)
        match_found = .false.
        k = 1
        do while (.not.match_found)
          iFc1 = mpg%divFc(mpg%divFcP(i,1)+k-1)
!nh check if the vertices belong to a flux surface
          if ((mpg%vxFs(mpg%fcVx(iFc1,1)).ne.0).and.    
     &      (mpg%vxFs(mpg%fcVx(iFc1,2)).ne.0)) then ! linear interpolation
            Psi1 = 0.5_R8 * (gm%fsPsi(mpg%vxFs(mpg%fcVx(iFc1,1)))+
     &                     gm%fsPsi(mpg%vxFs(mpg%fcVx(iFc1,2))))
            match_found = .true.
          elseif (mpg%vxFs(mpg%fcVx(iFc1,1)).ne.0) then ! take value from vertex 1
            Psi1 = gm%fsPsi(mpg%vxFs(mpg%fcVx(iFc1,1)))
            match_found = .true.
          elseif (mpg%vxFs(mpg%fcVx(iFc1,2)).ne.0) then ! take value from vertex 2
            Psi1 = gm%fsPsi(mpg%vxFs(mpg%fcVx(iFc1,2)))
            match_found = .true.
          else
            if (k.eq.mpg%divFcP(i,2)) then ! last face
              Psi1 = 0.0_R8
              write (*,*) 'Warning: no flux surface assigned to any ' //
     &         'of the grid vertices coinciding with '
              write (*,*) 'target ', i 
              write (*,*) 'Assume Psi = 0'
              match_found = .true.
            else
              k = k+1 ! check next face
            end if
          end if 
        end do
        match_found = .false.
        k = mpg%divFcP(i,2)
        do while (.not.match_found)
          iFc2 = mpg%divFc(mpg%divFcP(i,1)+k-1)
!nh check if the vertices belong to a flux surface
          if ((mpg%vxFs(mpg%fcVx(iFc2,1)).ne.0).and.    
     &      (mpg%vxFs(mpg%fcVx(iFc2,2)).ne.0)) then ! linear interpolation
            Psi2 = 0.5_R8 * (gm%fsPsi(mpg%vxFs(mpg%fcVx(iFc2,1)))+
     &                     gm%fsPsi(mpg%vxFs(mpg%fcVx(iFc2,2))))
            match_found = .true.
          elseif (mpg%vxFs(mpg%fcVx(iFc2,1)).ne.0) then ! take value from vertex 1
            Psi2 = gm%fsPsi(mpg%vxFs(mpg%fcVx(iFc2,1)))
            match_found = .true.
          elseif (mpg%vxFs(mpg%fcVx(iFc2,2)).ne.0) then ! take value from vertex 2
            Psi2 = gm%fsPsi(mpg%vxFs(mpg%fcVx(iFc2,2)))
            match_found = .true.
          else
            if (k.eq.1) then ! last face
              Psi2 = 0.0_R8
              match_found = .true.
            else
              k = k-1 ! check next face
            end if
          end if 
        end do
        if ((gm%psi_increasing .and. Psi1.gt.Psi2) .or.
     &      (.not.gm%psi_increasing .and. Psi2.gt.Psi1)) then
          allocate(old_face_list(1:mpg%divFcP(i,2)))
          old_face_list(1:mpg%divFcP(i,2)) =
     &     mpg%divFc(mpg%divFcP(i,1):mpg%divFcP(i,1)+mpg%divFcP(i,2)-1)
          do j = 1, mpg%divFcP(i,2)
            mpg%divFc(mpg%divFcP(i,1)+j-1) =
     &       old_face_list(mpg%divFcP(i,2)-j+1)
          end do
          deallocate(old_face_list)
        end if
        match_found = .false.
        do j = mpg%divFcP(i,1), mpg%divFcP(i,1)+mpg%divFcP(i,2)-1
          iVx1 = mpg%fcVx(mpg%divFc(j),1)
          iVx2 = mpg%fcVx(mpg%divFc(j),2)
          do k = 1, mpg%nXpt
            do l = mpg%strVxP(k,1), mpg%strVxP(k,1)+mpg%strVxP(k,2)-1
              if (iVx1.eq.mpg%strVx(l).and.
     &            mpg%vxFs(iVx1).eq.mpg%vxFs(mpg%Xpt(k))) then
                if((gm%psi_increasing .and.
     &              gm%vxFpsi(iVx1).le.gm%vxFpsi(iVx2)).or.
     &             (.not.gm%psi_increasing .and.
     &              gm%vxFpsi(iVx1).gt.gm%vxFpsi(iVx2))) then
                  if (mpg%ifdiv(i).eq.0) then
                    match_found = .true.
                    mpg%ifdiv(i) = j - mpg%divFcP(i,1) + 1
                    mpg%ivdiv(i) = iVx1
                  else ! Two strike points on the same target
                       ! We label the active one only
                    if (mpg%vxFs(iVx1).eq.mpg%iFssep) then
                      mpg%ifdiv(i) = j - mpg%divFcP(i,1) + 1
                      mpg%ivdiv(i) = iVx1
                    end if
                  end if
                end if
              else if (iVx2.eq.mpg%strVx(l).and.
     &                 mpg%vxFs(iVx2).eq.mpg%vxFs(mpg%Xpt(k))) then
                if((gm%psi_increasing .and.
     &              gm%vxFpsi(iVx2).le.gm%vxFpsi(iVx1)).or.
     &             (.not.gm%psi_increasing .and.
     &              gm%vxFpsi(iVx2).gt.gm%vxFpsi(iVx1))) then
                  if (mpg%ifdiv(i).eq.0) then
                    match_found = .true.
                    mpg%ifdiv(i) = j - mpg%divFcP(i,1) + 1
                    mpg%ivdiv(i) = iVx2
                  else ! Two strike points on the same target
                       ! We label the active one only
                    if (mpg%vxFs(iVx2).eq.mpg%iFssep) then
                      mpg%ifdiv(i) = j - mpg%divFcP(i,1) + 1
                      mpg%ivdiv(i) = iVx2
                    end if
                  end if
                end if
              end if
            end do
          end do
        end do
        call xertst (match_found .or. mpg%nXpt.eq.0,
     &   'No matching face found for strike point !')
      end do

      return
      end subroutine init_geometry


***********************************************************************

      subroutine calc_conn (mpg, gm)
      use b2us_map
      implicit none
      type (mapping), intent (in) :: mpg
      type (geometry), intent (inout) :: gm
      integer :: nfx, iCv, iCv1, iCv2, ic
      integer :: iFc, iFs, iFt, iVx, i, j
      integer, save :: conn_style = 0, iout = 0
      real (kind=r8) :: wrk0(mpg%nCv), wrk1(mpg%nCv), L1, L2, LL,
     .  tt0(mpg%nFt), tt1(mpg%nFt), tt2(mpg%nFt), tt3(mpg%nFt)
      external ipgeti, my_out_us, xerrab

      call ipgeti ('b2mndr_conn_style', conn_style)
      call ipgeti ('b2mndr_conn_iout', iout)

      ! build poloidal length of cells into wrk1
      wrk1 = 0.0_R8
      do iFt = 1, mpg%nFt 
        if (mpg%ftFcP(iFt,2) .gt. 0) then
          do j = mpg%ftFcP(iFt,1), mpg%ftFcP(iFt,1)+mpg%ftFcP(iFt,2)-1
            iFc = mpg%ftFc(j)
            iCv1 = mpg%fcCv(iFc,1)
            iCv2 = mpg%fcCv(iFc,2)
            wrk1(iCv1) = wrk1(iCv1) + gm%fcHc(iFc,1)
            wrk1(iCv2) = wrk1(iCv2) + gm%fcHc(iFc,2)
          enddo
        elseif (mpg%ftFcP(iFt,2) .eq. 0) then
          ! Special case of flux tubes consisting of only guard cells
          ! along plasma-void interfaces => these do not have faces
          do j = mpg%ftCvP(iFt,1), mpg%ftCvP(iFt,1)+mpg%ftCvP(iFt,2)-1
            iCv = mpg%ftCv(j)
            if (iCv .le. mpg%nCi) then
              call xerrab ('expected guard cell')
            endif
            iFc = mpg%cvFc(mpg%cvFcP(iCv,1))
            wrk1(iCv) = gm%fcHt(iFc)
          enddo
        else
          call xerrab ('Problem with flux tube definition')
        endif
      enddo

      ! construct wrk arrays and connection length of flux tubes
      wrk0 = 0.0_R8
      tt0 = 0.0_R8
      tt1 = 0.0_R8
      tt2 = 0.0_R8
      tt3 = 0.0_R8
      gm%ftConn = 0.0_R8
      do iFt = 1, mpg%nFt 
        do j = mpg%ftCvP(iFt,1), mpg%ftCvP(iFt,1)+mpg%ftCvP(iFt,2)-1
          iCv = mpg%ftCv(j)
          wrk0(iCv) = wrk1(iCv)*gm%cvBb(iCv,3)/gm%cvBb(iCv,0)
          tt0(iFt) = tt0(iFt) + wrk1(iCv)
          tt1(iFt) = tt1(iFt) + wrk1(iCv)*gm%cvVol(iCv)/gm%cvSz(iCv)
          tt2(iFt) = tt2(iFt) + gm%cvBb(iCv,3)*gm%cvBb(iCv,3)*
     .        gm%cvVol(iCv)
          tt3(iFt) = tt3(iFt) + gm%cvVol(iCv)
          gm%ftConn(iFt) = gm%ftConn(iFt) + wrk0(iCv)
        enddo
      enddo

      ! construct ftEps and ftBbav2
      gm%ftEps = 0.0_R8
      gm%ftBbav2 = 0.0_R8
      do iFt = 1, mpg%nFt 
        if (mpg%cvOnClosedSurface(mpg%ftCv(mpg%ftCvP(iFt,1)))) then
          gm%ftEps(iFt) = tt0(iFt)**2/tt1(iFt)
          gm%ftBbav2(iFt) = tt2(iFt)/tt3(iFt)
        endif
      enddo

      ! 'effective' connection length for use in k-eps model
      if (conn_style.ne.0.or.maxval(gm%cvConn).eq.0.0_R8) then ! if 0, assumed to be read from b2fgmtry
       gm%cvConn = 0.0_R8
       do iFt = 1, mpg%nFt 
        iCv = mpg%ftCv(mpg%ftCvP(iFt,1))
        if (mpg%cvOnClosedSurface(iCv)) then
         !core region: actual connection length (same as ftConn)
         gm%cvConn(mpg%ftCv(mpg%ftCvP(iFt,1):mpg%ftCvP(iFt,1)+
     &                        mpg%ftCvP(iFt,2)-1)) = gm%ftConn(iFt)
        else
         ! SOL/PF: effective distance to targets/boundaries
         ! implicit assumption that cells in ftCv are listed starting
         ! at one  of the targets
         do ic = mpg%ftCvP(iFt,1),mpg%ftCvP(iFt,1)+mpg%ftCvP(iFt,2)-1
          iCv = mpg%ftCv(ic)
          L1 = 0.5_R8*wrk0(iCv)
          do j = ic+1, mpg%ftCvP(iFt,1) + mpg%ftCvP(iFt,2)-1
           L1 = L1 + wrk0(mpg%ftCv(j))
          end do
          L2 = 0.5_R8*wrk0(iCv)
          do j = mpg%ftCvP(iFt,1), ic-1
           L2 = L2 + wrk0(mpg%ftCv(j))
          end do
          if (conn_style.ne.2) then
           gm%cvConn(iCv) = min(L1,L2)
          elseif (conn_style.eq.2) then
           call xerrab ('Use of conn_style.eq.2 is deprecated')
          endif
         end do
        endif
       end do
      endif

      ! connection length evaluated for flux surfaces and vertices
      gm%fsConn = 0.0_R8
      gm%vxConn = 0.0_R8
      do iFs = 1, mpg%nFs
        do i = mpg%fsFcP(iFs,1), mpg%fsFcP(iFs,1)+mpg%fsFcP(iFs,2)-1
          iFc = mpg%fsFc(i)
          gm%fsConn(iFs) = gm%fsConn(iFs) +
     &                     gm%fcHt(iFc)*gm%fcBb(iFc,3)/gm%fcBb(iFc,0)
        end do
        iFc = mpg%fsFc(mpg%fsFcP(iFs,1))
        iCv1 = mpg%fcCv(iFc,1)
        iCv2 = mpg%fcCv(iFc,2)
        if (((mpg%cvOnClosedSurface(iCv1).or.iCv1.gt.mpg%nCi) .and.
     &       (mpg%cvOnClosedSurface(iCv2).or.iCv2.gt.mpg%nCi)) .or.
     &        iFs.eq.mpg%iFssep ) then
          do j = mpg%fsVxP(iFs,1), mpg%fsVxP(iFs,1)+mpg%fsVxP(iFs,2)-1
            iVx = mpg%fsVx(j)
            gm%vxConn(iVx) = gm%fsConn(iFs)
          end do
        else
          do i = mpg%fsFcP(iFs,1), mpg%fsFcP(iFs,1)+mpg%fsFcP(iFs,2)-1
            iFc = mpg%fsFc(i)

            ! First follow first vertex
            iVx = mpg%fcVx(iFc,1)
            L1 = 0.0_R8
            call follow_flux_surface(mpg, gm, ifs, iFc, iVx, L1)

            ! Second follow second vertex
            iVx = mpg%fcVx(iFc,2)
            L2 = 0.0_R8
            call follow_flux_surface(mpg, gm, ifs, iFc, iVx, L2)

            ! connection length of current face
            LL = gm%fcHt(iFc)*gm%fcBb(iFc,3)/gm%fcBb(iFc,0)

            ! first vertex --> L2 must be increased by LL
            iVx = mpg%fcVx(iFc,1)
            gm%vxConn(iVx) = min(L2+LL,L1)
            ! second vertex --> L1 must be increased by LL
            iVx = mpg%fcVx(iFc,2)
            gm%vxConn(iVx) = min(L2,L1+LL)
          end do
        end if
      end do
      if (iout.ne.0) then
       call my_out_us(70,mpg%nCv,0,gm%cvConn,'calc_conn_cvConn')
       call my_out_us(70,mpg%nFt,0,gm%ftConn,'calc_conn_ftConn')
       call my_out_us(70,mpg%nFt,0,gm%ftEps,'calc_conn_ftEps')
       call my_out_us(70,mpg%nFt,0,gm%ftBbav2,'calc_conn_ftBbav2')
      endif

      return
      end subroutine calc_conn

**********************************************************************

      subroutine is_periodic ( mpg, gm )
      use b2us_map
      implicit none
      type(mapping), intent(inout) :: mpg
      type(geometry), intent(in) :: gm
      integer :: iCvm, iCvn, iCvp, iCv, iFc, iFt, icount
      logical :: closed, end_found, all_closed

      mpg%periodic_bc = 0
      if (mpg%nnreg(0).eq.2 .or. mpg%nnreg(0).eq.5) then
        mpg%periodic_bc = 1
      else if (mpg%nnreg(0).eq.1) then
        all_closed = mpg%nFt.gt.0
        do iFt = 1, mpg%nFt
          closed = .false.
          iCv = mpg%ftCv(mpg%ftCvP(iFt,1))
          if (mpg%ftCvP(iFt,2).gt.1) then
            iCvm = iCv
            iCvn = 0
            do iFc = mpg%cvFcP(iCv,1),
     &               mpg%cvFcP(iCv,1) + mpg%cvFcP(iCv,2) - 1
              if (gm%fcQalf(mpg%cvFc(iFc),0).eq.0.0_R8) cycle
              if (mpg%fcCv(mpg%cvFc(iFc),1).eq.iCv) then
                iCvn = mpg%fcCv(mpg%cvFc(iFc),2)
              else
                iCvn = mpg%fcCv(mpg%cvFc(iFc),1)
              end if
            end do
            closed = iCvn.eq.0
            end_found = iCvn.eq.0
            do while (.not.end_found)
              iCvp = 0
              do iFc = mpg%cvFcP(iCvn,1),
     &                 mpg%cvFcP(iCvn,1) + mpg%cvFcP(iCvn,2) - 1
                if (gm%fcQalf(mpg%cvFc(iFc),0).eq.0.0_R8) cycle
                if (mpg%fcCv(mpg%cvFc(iFc),1).eq.iCvm) cycle
                if (mpg%fcCv(mpg%cvFc(iFc),2).eq.iCvm) cycle
                if (mpg%fcCv(mpg%cvFc(iFc),1).eq.iCvn) then
                  iCvp = mpg%fcCv(mpg%cvFc(iFc),2)
                else
                  iCvp = mpg%fcCv(mpg%cvFc(iFc),1)
                end if
              end do
              iCvm = iCvn
              iCvn = iCvp
              closed = iCvp.eq.iCv
              end_found = closed .or. iCvp.eq.0
            end do
          else if (mpg%ftCvP(iFt,2).eq.1) then
            icount = 0
            iCv = mpg%ftCv(mpg%ftCvP(iFt,1))
            do iFc = mpg%cvFcP(iCv,1), mpg%cvFcP(iCv,1) +
     &                                 mpg%cvFcP(iCv,2) - 1
              if (gm%fcQalf(mpg%cvFc(iFc),0).eq.0.0_R8) cycle
              if (mpg%fcCv(mpg%cvFc(iFc),1).eq.iCv .and.
     &            mpg%fcCv(mpg%cvFc(iFc),2).eq.iCv) icount = icount + 1
            end do
            closed = icount.eq.2
          else if (mpg%ftCvP(iFt,2).eq.0) then
            closed = .true.
          end if
          all_closed = all_closed .and. closed
        end do
        if (all_closed) mpg%periodic_bc = 1
      end if

      return
      end subroutine is_periodic

***********************************************************************

      subroutine divide_by_area( nFc, geo, flux, flux_density )
      implicit none
      integer, intent(in) :: nFc
      type(geometry), intent(in) :: geo
      real(kind=R8), intent(in) :: flux( nFc, 0:1 )
      real(kind=R8), intent(out) :: flux_density( nFc, 0:1 )
      integer iFc

      do iFc = 1, nFc
        if ( geo%fcQalf(iFc,0).ne.0.0_R8 ) then
          flux_density(iFc,0) = flux(iFc,0)/
     &      (geo%fcS(iFc)*geo%fcQalf(iFc,0))
        else
          flux_density(iFc,0) = 0.0_R8
        end if
        if ( geo%fcQalf(iFc,1).ne.0.0_R8 ) then
          flux_density(iFc,1) = flux(iFc,1)/
     &      (geo%fcS(iFc)*geo%fcQalf(iFc,1))
        else
          flux_density(iFc,1) = 0.0_R8
        end if
      end do

      return
      end subroutine divide_by_area

***********************************************************************

      subroutine split_face_flux( mpg, geo, flxFace, tmpFace )
      use b2us_map
      implicit none
      type(mapping), intent(in) :: mpg
      type(geometry), intent(in) :: geo
      real(kind=R8), intent(in) :: flxFace( mpg%nFc, 0:1 )
      real(kind=R8), intent(out) :: tmpFace( mpg%nFc )
      integer iFc, nFx

cxpb  This routine is a hack to assign the poloidal value to 'poloidal' faces
cxpb  and the radial value to 'radial' faces.

      if (mpg%isClassicalGrid.eq.1) then
        nFx = maxval(mpg%imapFcx)
        do iFc = 1, nFx
          tmpFace(iFc) = flxFace(iFc,0)
        end do
        do iFc = nFx+1, mpg%nFc
          tmpFace(iFc) = flxFace(iFc,1)
        end do
      else
        do iFc = 1, mpg%nFc
          if ( geo%fcQalf(iFc,0).eq.0.0_R8 ) then
            tmpFace(iFc) = flxFace(iFc,1)
          else
            tmpFace(iFc) = flxFace(iFc,0)
          end if
        end do
      end if

      return
      end subroutine split_face_flux

***********************************************************************

      subroutine follow_flux_surface( mpg, geo, ifs, ifc, ivx, LL )
      use b2us_map
      implicit none
      type(mapping), intent(in) :: mpg
      type(geometry), intent(in) :: geo
      integer, intent(in) :: ifs
      integer, intent(in) :: ifc, ivx
      real(kind=R8), intent (inout) :: LL
      integer fc_current, vx_current, fc_next, vx1, vx2, j
      logical found_next

csc   This routine follows a flux surface starting from a vertex of a face
csc   until a target, and calculates its connection length.
csc   The idea is: find a face that is on same flux tube but is different
csc   than previous face, then identify next vertex
      fc_current = ifc
      vx_current = ivx
      found_next = .true.
      do while (found_next)
         ! we have a new face, loop over all faces connected to the current vertex
         found_next = 0
         do j = mpg%vxFcP(vx_current,1), 
     &          mpg%vxFcP(vx_current,1)+mpg%vxFcP(vx_current,2)-1
           fc_next = mpg%vxFc(j)
           if (fc_next==fc_current) cycle  ! we ingore the current face
           ! check if face belongs to same flux tube
           if (mpg%fcFs(fc_next)==ifs) then
             ! face belongs to flux tube and is not the previous face
             found_next = .true.
             LL = LL + geo%fcHt(fc_next)*geo%fcBb(fc_next,3)/
     &                 geo%fcBb(fc_next,0)
             vx1 = mpg%fcVx(fc_next,1)
             vx2 = mpg%fcVx(fc_next,2)
             ! identify vertex of next face which is not the current vertex
             if (vx1==vx_current) then
               vx_current = vx2
             else
               vx_current = vx1
             endif
             ! update the current face number
             fc_current = fc_next
           endif
           ! stop the for loop if we found the next face
           if (found_next) exit
         enddo
      end do! while


      return
      end subroutine follow_flux_surface

***********************************************************************

      function pitch(g,iCv)
      implicit none
      integer, intent(in) :: iCv
      type(geometry), intent(in) :: g
      real(kind=R8) :: pitch
      pitch = sqrt(g%cvBb(iCv,0)**2+g%cvBb(iCv,1)**2)/g%cvBb(iCv,3)
      return
      end function pitch

***********************************************************************

      end module b2us_geo
