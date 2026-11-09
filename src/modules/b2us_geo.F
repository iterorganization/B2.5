      module b2us_geo

      use b2mod_types
      implicit none

      private
      public :: alloc_geometry, dealloc_geometry,
     .          read_geometry, write_geometry,
     .          init_geometry, init_walls_and_targets,
     .          calc_conn, check_geometry,
     .          divide_by_area, split_face_flux, pitch

      type, public :: geometry

!  cell quantites
        real (kind=R8), allocatable  :: cvBb(:,:)      ! Magnetic field: (nCv,0:3) Bx, By, Bz, |B|
        real (kind=R8), allocatable  :: cvX(:)         ! X-coordinate of cell center
        real (kind=R8), allocatable  :: cvY(:)         ! Y-coordinate of cell center
        real (kind=R8), allocatable  :: cvSz(:)        ! Area perp. to third (z) direction
        real (kind=R8), allocatable  :: cvHz(:)        ! Length in third (z) direction
        real (kind=R8), allocatable  :: cvHx(:)        ! Temp: characteristic poloidal length; eventually remove, because meaning unclear in WG context
        real (kind=R8), allocatable  :: cvHy(:)        ! Estimate of the radial width of a cell
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
        real (kind=R8), allocatable  :: fcQalf(:,:)    ! Cosine (:,0) and sine (:,1) of alpha for positive Bz;
                                                       ! Cosine (:,0) and (-1)*sine (:,1) of alpha for negative Bz
        real (kind=R8), allocatable  :: fcQbet(:,:)    ! Cosine (:,0) and sine (:,1) of beta for positive Bz;
                                                       ! (-1)*Cosine (:,0) and sine (:,1) of beta for negative Bz
        real (kind=R8), allocatable  :: fcPbs(:)       ! Magnetic contact area
        real (kind=R8), allocatable  :: fcPbshz(:)     ! Magnetic contact area * hz (at faces)
        real (kind=R8), allocatable  :: fcBzb(:)       ! Factor for P-S/mdf flows
        real (kind=R8), allocatable  :: fcEb(:,:)      ! Unit vector pointing along magnetic field (Cartesian X,Y,Z coordinates)
        real (kind=R8), allocatable  :: fcFpsi(:)      ! Psi in each face center

!  vertex quantities
        real (kind=R8), allocatable  :: vxBb(:,:)      ! Magnetic field: (nVx,0:3) to be recomputed at vertices (Bx,By,Bz,|B|)
        real (kind=R8), allocatable  :: vxX(:)         ! X-coordinate of vertex
        real (kind=R8), allocatable  :: vxY(:)         ! Y-coordinate of vertex
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

        real (kind=R8) :: signmf  ! Direction of poloidal field
                                  ! +1: Going clockwise around plasma
                                  ! -1: Going counterclockwise around plasma

        real (kind=R8) :: Qalfmin, Qalfmax
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
      allocate (g%cvHy(nCv))
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
      g%cvHy = 0._R8
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

      g%Qalfmax = 1.0e-3_R8
      g%Qalfmin = 1.0e-3_R8

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
      deallocate (g%cvHy)
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
      real(kind=R8) :: rdum(1)
      external cfrure, cfrure_opt, ipgetr

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

      rdum(1) = 1.0e-3_R8
      call ipgetr('b2us_prep_Qalfmin', rdum(1))
      call cfrure_opt (iun,   1,   rdum, 'Qalfmin')
      gm%Qalfmin = rdum(1)
      rdum(1) = 1.0e-3_R8
      call ipgetr('b2us_prep_Qalfmax', rdum(1))
      call cfrure_opt (iun,   1,   rdum, 'Qalfmax')
      gm%Qalfmax = rdum(1)

! check consistency of input
      call check_geometry(nCv, nFc, nVx, nFs, nFt, gm)

      return
      end subroutine read_geometry

***********************************************************************

      subroutine write_geometry (nout, nCv, nFc, nVx, nFs, gm)

      implicit none
      integer, intent(in) :: nout, nCv, nFc, nVx, nFs
      type (geometry), intent (in) :: gm
      real(kind=R8) :: rdum(1)
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

      rdum(1) = gm%Qalfmin
      call cfwure (nout,   1,   rdum,  'Qalfmin')
      rdum(1) = gm%Qalfmax
      call cfwure (nout,   1,   rdum,  'Qalfmax')

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
     &           iCv, iFc, iVx, inCv, iXpt, iFc1, iFc2, iVx1, iVx2, iFt,
     &           inv_dist(mpg%nVx)
      integer :: count_up, count_down, count_eq
      integer, allocatable :: old_face_list(:), verts(:)
      integer, save :: temp_ignore_geo = 0
      real (kind=R8) :: hzconst, r0, z0, t0, dux, duy, du, sbf,
     .  Psi1, Psi2, dpsi(mpg%nCv), psi_active,
     .  dpsi_max_down, dpsi_max_up
      logical match_found
      intrinsic maxval, abs, sqrt
      external b2xbzb, intface, intvertex, ipgeti, xerrab, xertst

      call ipgeti ('b2mndr_temp_ignore_geo', temp_ignore_geo)

      nCv = mpg%nCv
      nFc = mpg%nFc
      nVx = mpg%nVx

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
      do iFt = 1, mpg%nFt
        gm%ftFpsi(iFt) = gm%cvFpsi(mpg%ftCv(mpg%ftCvP(iFt,1)))
        do i = 1, mpg%ftCvP(iFt,2)
          iCv = mpg%ftCv(mpg%ftCvP(iFt,1)+i-1)
          if (gm%cvFpsi(iCv).ne.gm%ftFpsi(iFt)) then
            write (*,'(a,i6,a,i4)')
     &       'Warning: cvFpsi value in cell ', iCv,
     &       ' not equal to value in its flux tube ', iFt
            write (*,*) gm%cvFpsi(iCv), gm%ftFpsi(iFt)
          endif
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
     &       (1._R8/mpg%cvVxP(iCv,2))*gm%cvVol(iCv)
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
              gm%vxVol(mpg%vxCvP(iVx,1) + inCv - 1) =
     &         1._R8/mpg%cvVxP(iCv,2)
            else
              ! Make sure guard cells dominate in the interpolation
              ! to vertices => use volume
              gm%vxVol(mpg%vxCvP(iVx,1) + inCv - 1) =
     &         (1._R8/mpg%cvVxP(iCv,2))*gm%cvVol(iCv)
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
          count_down = 0
          count_up   = 0
          count_eq   = 0
          do inCv = 1, mpg%vxCvP(iVx,2)
            iCv = mpg%vxCv(mpg%vxCvP(iVx,1) + inCv - 1)
            if (gm%cvFpsi(iCv) .gt. gm%vxFpsi(iVx)) then
              count_up = count_up + 1
            elseif (gm%cvFpsi(iCv) .lt. gm%vxFpsi(iVx)) then
              count_down = count_down + 1
            else
              count_eq = count_eq + 1
            endif
          enddo
          do inCv = 1, mpg%vxCvP(iVx,2)
            iCv = mpg%vxCv(mpg%vxCvP(iVx,1) + inCv - 1)
            if (iCv.le.mpg%nCi) then
              if ( gm%cvFpsi(iCv) .gt. gm%vxFpsi(iVx)) then
                gm%vxVol(mpg%vxCvP(iVx,1) + inCv - 1) =
     &           2.0_R8*count_up
              elseif ( gm%cvFpsi(iCv) .lt. gm%vxFpsi(iVx)) then
                gm%vxVol(mpg%vxCvP(iVx,1) + inCv - 1) =
     &           2.0_R8*count_down
              else
                gm%vxVol(mpg%vxCvP(iVx,1) + inCv - 1) =
     &           count_eq
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
     &           gm%cvVol(iCv) / 1.0e2_R8
              endif
            enddo
          endif
        enddo
      elseif (switch%vxVol_style .eq. 7 .or.
     &        switch%vxVol_style .eq. 8 ) then
        ! wdk Style 7 and 8 have become basically equal to vxVol_style 2, except for treatment of guard cells
        ! wdk need to remove obsolete options eventually

        ! Assumption: all internal 'radial' faces aligned (boundaries faces possibly non-aligned)
        ! Idea: assign equal total weight to all cells on either side of the flux surface (inverse volume weighting!)
        !       (expected: in case of zero gradient along poloidal direction, value at vertex will be
        !                  average between flux tube values)
        ! Count cells per vertex with higher, lower Psi
        ! Exceptional case: equal Psi (pentagon splitting);
        ! will also catch case where Psi = 0 everywhere
        do iVx = 1, nVx
          count_down = 0
          count_up   = 0
          count_eq   = 0
          do inCv = 1, mpg%vxCvP(iVx,2)
            iCv = mpg%vxCv(mpg%vxCvP(iVx,1) + inCv - 1)
            if ( gm%cvFpsi(iCv) .gt. gm%vxFpsi(iVx) ) then
              count_up = count_up + 1
            elseif ( gm%cvFpsi(iCv) .lt. gm%vxFpsi(iVx) ) then
              count_down = count_down + 1
            else
              count_eq = count_eq + 1
            endif
          enddo
          do inCv = 1, mpg%vxCvP(iVx,2)
            iCv = mpg%vxCv(mpg%vxCvP(iVx,1) + inCv - 1)
            if (iCv.le.mpg%nCi) then
              if ( gm%cvFpsi(iCv) .gt. gm%vxFpsi(iVx)) then
                gm%vxVol(mpg%vxCvP(iVx,1) + inCv - 1) =
     &            2.0_R8*count_up
              elseif ( gm%cvFpsi(iCv) .lt. gm%vxFpsi(iVx)) then
                gm%vxVol(mpg%vxCvP(iVx,1) + inCv - 1) =
     &            2.0_R8*count_down
              else
                gm%vxVol(mpg%vxCvP(iVx,1) + inCv - 1) =
     &            real(count_eq,R8)
              endif
            else
              ! Make guards cells dominate and account for
              ! length difference via distance to vertex
              gm%vxVol(mpg%vxCvP(iVx,1) + inCv - 1) = 0.25e-3_R8 *
     &             sqrt( (gm%cvX(iCv) - gm%vxX(iVx))**2  +
     &                   (gm%cvY(iCv) - gm%vxY(iVx))**2  )
            endif
          enddo
        enddo
      elseif (switch%vxVol_style .eq. 9) then
        !count cells per vertex with higher and lower Psi
        do iVx = 1, nVx
          count_down = 0
          count_up   = 0
          dpsi_max_down = 0.0_R8
          dpsi_max_up   = 0.0_R8
          do inCv = 1, mpg%vxCvP(iVx,2)
            iCv = mpg%vxCv(mpg%vxCvP(iVx,1)+ inCv - 1)
            if ( gm%cvFpsi(iCv) .gt. gm%vxFpsi(iVx) ) then
              count_up = count_up + 1
              if ( dpsi(iCv) .gt. dpsi_max_up) then
                dpsi_max_up = dpsi(iCv)
              endif
            elseif ( gm%cvFpsi(iCv) .lt. gm%vxFpsi(iVx) ) then
              count_down = count_down + 1
              if ( dpsi(iCv) .gt. dpsi_max_down) then
                dpsi_max_down = dpsi(iCv)
              endif
            endif
          enddo
          do inCv = 1, mpg%vxCvP(iVx,2)
            iCv = mpg%vxCv(mpg%vxCvP(iVx,1)+ inCv - 1)
            if (iCv .le. mpg%nCi) then
              if ( gm%cvFpsi(iCv) .gt. gm%vxFpsi(iVx) ) then
                if (mpg%cvVxP(iCv,2) .eq. 3) then
                  gm%vxVol(mpg%vxCvP(iVx,1)+inCv - 1) =
     &               real(count_up,R8) / (mpg%vxCvP(iVx,2)*sqrt(2.0_R8)*
     &               sqrt(dpsi_max_up/dpsi(iCv)))
                else
                  gm%vxVol(mpg%vxCvP(iVx,1)+inCv - 1) =
     &               real(count_up,R8) / mpg%vxCvP(iVx,2)
                endif
              elseif ( gm%cvFpsi(iCv) .lt. gm%vxFpsi(iVx) ) then
                if (mpg%cvVxP(iCv,2) .eq. 3) then
                  gm%vxVol(mpg%vxCvP(iVx,1)+inCv - 1) =
     &             real(count_down,R8) / (mpg%vxCvP(iVx,2)*sqrt(2.0_R8)*
     &                   sqrt(dpsi_max_up/dpsi(iCv)))
                else
                  gm%vxVol(mpg%vxCvP(iVx,1)+inCv - 1) =
     &             real(count_down,R8) / mpg%vxCvP(iVx,2)
                endif
              endif
            else
              ! Make guards cells dominate and account for
              ! length difference via distance to vertex
              gm%vxVol(mpg%vxCvP(iVx,1)+inCv - 1) = 0.25e-3_R8 *
     &             sqrt( (gm%cvX(iCv) - gm%vxX(iVx))**2 +
     &                   (gm%cvY(iCv) - gm%vxY(iVx))**2 )
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
     &          (mpg%cvVxP(iCv,2).gt.4)) then
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

        ! Normalize
        du = sqrt(gm%fcEb(iFc,0)**2+gm%fcEb(iFc,1)**2+gm%fcEb(iFc,2)**2)
        gm%fcEb(iFc,0) = gm%fcEb(iFc,0)/du
        gm%fcEb(iFc,1) = gm%fcEb(iFc,1)/du
        gm%fcEb(iFc,2) = gm%fcEb(iFc,2)/du

      end do

!WG_TODO: Need to provide a calculation for gm%vxEb

*   ..estimate the radial width of cells: cvHy
      do iCv = 1, nCv
        do i = 1, mpg%cvFcP(iCv,2)
          iFc = mpg%cvFc(mpg%cvFcP(iCv,1)+i-1)
          gm%cvHy(iCv) = gm%cvHy(iCv)+
     .                   gm%fcHt(iFc)*abs(gm%fcQalf(iFc,0))
        enddo
      enddo
      do iCv = 1, mpg%nCi
        gm%cvHy(iCv) = 0.5_R8*gm%cvHy(iCv)
      enddo

*   ..compute cvHz
      hzconst = 2.0_R8*pi*6.0_R8   ! for ITER for testing
      hzconst = 2.0_R8*pi*1.6_R8   ! for ASDEX-UPGRADE for testing
      hzconst = 1.0_R8
      gm%cvHz = hzconst*(1.0_R8-switch%b2mndr_hz)+switch%b2mndr_hz*
     & gm%cvVol/gm%cvSz

*   ..compute fcHz and fcPbshz
      if (switch%fcHz_style.eq.0) then
        call intface(nCv,nFc,mpg%fcCv,gm%fcVol,gm%cvHz,gm%fcHz)
      else
        if (isymm.eq.0) then
          gm%fcHz = hzconst*(1.0_R8-switch%b2mndr_hz)+
     &      switch%b2mndr_hz*switch%b2agmt_width
        elseif (isymm.eq.1 .or. isymm.eq.2) then
          do iFc = 1, mpg%nFc
            gm%fcHz(iFc) = hzconst*(1.0_R8-switch%b2mndr_hz)+
     &        switch%b2mndr_hz*
     &        pi*(gm%vxX(mpg%fcVx(iFc,1))+gm%vxX(mpg%fcVx(iFc,2)))
          enddo
        else if (isymm.eq.3 .or. isymm.eq.4) then
          do iFc = 1, mpg%nFc
            gm%fcHz(iFc) = hzconst*(1.0_R8-switch%b2mndr_hz)+
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
          gm%vxHz = hzconst*(1.0_R8-switch%b2mndr_hz)+
     &      switch%b2mndr_hz*switch%b2agmt_width
        elseif (isymm.eq.1 .or. isymm.eq.2) then
          gm%vxHz = hzconst*(1.0_R8-switch%b2mndr_hz)+
     &      switch%b2mndr_hz*2.0_R8*pi*gm%vxX
        else if (isymm.eq.3 .or. isymm.eq.4) then
          gm%vxHz = hzconst*(1.0_R8-switch%b2mndr_hz)+
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
          nFx = maxval(mpg%imapFcx)
          ! 'poloidal' faces
          gm%fcQgam(1:nFx,1) = 0.0_R8
          gm%fcQalf(1:nFx,1) = 0.0_R8
          ! 'radial' faces
          gm%fcQgam(nFx+1:mpg%nFc,1) = 0.0_R8
          gm%fcQbet(nFx+1:mpg%nFc,0) = 0.0_R8
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
      if ((mpg%nnreg(0).eq.1.and.mpg%isClassicalGrid.eq.1).or.
     &    .not.allocated(mpg%Xpt)) then
!! TODO: Find a way to ascertain, if flux surfaces are defined,
!!       whether the psi value is increasing or decreasing radially
      elseif (mpg%nXpt.eq.0) then
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
     &                          maxval(dabs(gm%fsPsi)).eq.0.0_R8
          end if
        end do
      else if (mpg%nnreg(0).eq.5 .and. mpg%isClassicalGrid.eq.1) then
        iXpt = 1
      else if (mpg%nXpt.eq.1) then
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
      else if (mpg%nXpt.ge.2) then
        iXpt = mpg%iactive
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

      if (mpg%iactive.gt.0 .and. mpg%iinactive.eq.0 .and. mpg%nXpt.ge.2)
     > then
!! look for the inactive X-point as the one closest in psi-space to the
!! active X-point
        psi_active = gm%fsPsi(mpg%vxFs(mpg%Xpt(mpg%iactive)))
        j = 0
        do i = 1, mpg%nXpt
          if (i.eq.mpg%iactive) cycle
          if (j.eq.0) then
            j = i
          else
            if (abs( gm%fsPsi(mpg%vxFs(mpg%Xpt(i)))-psi_active ).lt.
     .          abs( gm%fsPsi(mpg%vxFs(mpg%Xpt(j)))-psi_active ))
     >        j = i
          end if
        end do
        mpg%iinactive = j
      end if

      ! Compute signmf
      r0=0.0_R8
      z0=0.0_R8
      t0=0.0_R8
      iFt=0
      gm%signmf=0.0_R8
      if (any(mpg%cvOnClosedSurface(1:mpg%nCi))) then
        do while (t0.eq.0.0_R8.and.iFt.lt.mpg%nFt)
          iFt=iFt+1
          iCv=mpg%ftCv(mpg%ftCvP(iFt,1))
          if (mpg%cvOnClosedSurface(iCv)) then
            do i=1,mpg%ftCvP(iFt,2)
              r0=r0+gm%cvX(mpg%ftCv(mpg%ftCvP(iFt,1)+i-1))*
     .              gm%cvHx(mpg%ftCv(mpg%ftCvP(iFt,1)+i-1))
              z0=z0+gm%cvY(mpg%ftCv(mpg%ftCvP(iFt,1)+i-1))*
     .              gm%cvHx(mpg%ftCv(mpg%ftCvP(iFt,1)+i-1))
              t0=t0+
     .              gm%cvHx(mpg%ftCv(mpg%ftCvP(iFt,1)+i-1))
            end do
          end if
        end do
      end if
      if (t0.gt.0.0_R8) then
        r0=r0/t0
        z0=z0/t0
        iFt=0
        do while (gm%signmf.eq.0.0_R8.and.iFt.lt.mpg%nFt)
          iFt=iFt+1
          iCv=mpg%ftCv(mpg%ftCvP(iFt,1))
          if (mpg%cvOnClosedSurface(iCv)) then
            do i=1,mpg%ftCvP(iFt,2)
              gm%signmf=gm%signmf
     .         -((gm%cvX(mpg%ftCv(mpg%ftCvP(iFt,1)+i-1))-r0)*
     .            gm%cvEb(mpg%ftCv(mpg%ftCvP(iFt,1)+i-1),1))
     .         +((gm%cvY(mpg%ftCv(mpg%ftCvP(iFt,1)+i-1))-z0)*
     .            gm%cvEb(mpg%ftCv(mpg%ftCvP(iFt,1)+i-1),0))
            end do
          end if
        end do
! signmf is positive if the poloidal field travels clockwise around the O-point
        gm%signmf=sign(1.0_R8,gm%signmf)
      else
        gm%signmf=sign(1.0_R8,gm%cvEb(1,0))
      end if

      if (temp_ignore_geo.eq.0) then !! This is a temporary solution to be able to run extended cases, while more general solutions are being developed
!! Invert the divertor face list if not in correct psi order
!! and identify strike point face indices
      if (allocated(mpg%strDiv)) then
       do i = 1, maxval(mpg%strDiv)
        match_found = .false.
        k = 1
        do while (.not.match_found)
          iFc1 = mpg%divFc(mpg%divFcP(i,1)+k-1)
!nh check if the vertices belong to a flux surface
          if ((mpg%vxFs(mpg%fcVx(iFc1,1)).ne.0).and.
     &        (mpg%vxFs(mpg%fcVx(iFc1,2)).ne.0)) then ! linear interpolation
            Psi1 = 0.5_R8 * (gm%fsPsi(mpg%vxFs(mpg%fcVx(iFc1,1)))+
     &                       gm%fsPsi(mpg%vxFs(mpg%fcVx(iFc1,2))))
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
     &        (mpg%vxFs(mpg%fcVx(iFc2,2)).ne.0)) then ! linear interpolation
            Psi2 = 0.5_R8 * (gm%fsPsi(mpg%vxFs(mpg%fcVx(iFc2,1)))+
     &                       gm%fsPsi(mpg%vxFs(mpg%fcVx(iFc2,2))))
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
              if (iVx1.eq.mpg%strVx(l)) then
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
                    if (mpg%vxFs(iVx1).eq.mpg%iFssep.or.
     &                 (mpg%vxFs(iVx1).eq.mpg%iFssep2.and.
     &                  mpg%iFssep2.gt.0)) then
                      mpg%ifdiv(i) = j - mpg%divFcP(i,1) + 1
                      mpg%ivdiv(i) = iVx1
                    end if
                  end if
                end if
              else if (iVx2.eq.mpg%strVx(l)) then
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
                    if (mpg%vxFs(iVx2).eq.mpg%iFssep.or.
     &                 (mpg%vxFs(iVx2).eq.mpg%iFssep2.and.
     &                  mpg%iFssep2.gt.0)) then
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
      end if

      else
        write (*,*) 'Warning: temp_ignore_geo switch used'
        write (*,*) 'This is a temporary solution to run unstructured',
     & ' geometries that do not fit in existing case categorizations.'
        write (*,*) 'This option limits post-processing options.'
      endif

      if (allocated(mpg%fcs_wall)) then
      ! Count nFci
        mpg%nFci = 0
        do j = 2, mpg%nCg
          if (maxval(mpg%cvFt(mpg%fcCv(mpg%fcs_wall(j), 1:2))) /=
     &        maxval(mpg%cvFt(mpg%fcCv(mpg%fcs_wall(j-1), 1:2)))) then
            mpg%nFci = mpg%nFci + 1
          end if
        end do
        mpg%nFci = mpg%nFci + 1

        allocate(mpg%fcs_wall_ind(mpg%nFci, 2),
     &           mpg%fcs_wall_type(mpg%nFci))

        ! Populate mapping array
        k = 1
        i = 1
        do j = 2, mpg%nCg
          if (maxval(mpg%cvFt(mpg%fcCv(mpg%fcs_wall(j), 1:2))) /=
     &        maxval(mpg%cvFt(mpg%fcCv(mpg%fcs_wall(j-1), 1:2)))) then
            mpg%fcs_wall_ind(k, 1) = i
            mpg%fcs_wall_ind(k, 2) = j - i
            if (mpg%cvOnClosedSurface
     &          (minval(mpg%fcCv(mpg%fcs_wall(i), 1:2)))) then
              mpg%fcs_wall_type(k) = 1
            else
              mpg%fcs_wall_type(k) = 2
            end if
            k = k + 1
            i = j
          end if
        end do

        ! Handle last entry
        mpg%fcs_wall_ind(mpg%nFci, 1) = mpg%fcs_wall_ind(mpg%nFci-1, 1)
     &    + mpg%fcs_wall_ind(mpg%nFci-1, 2)
        mpg%fcs_wall_ind(mpg%nFci, 2) = mpg%nCg -
     &    mpg%fcs_wall_ind(mpg%nFci, 1) + 1
        if (mpg%cvOnClosedSurface
     &      (minval(mpg%fcCv(mpg%fcs_wall(i), 1:2)))) then
          mpg%fcs_wall_type(k) = 1
        else
          mpg%fcs_wall_type(k) = 2
        end if
      endif

      write(*,*) "list of boundary faces in one flux tube and type:"
      do j=1,mpg%nFci
        write(*,'(4i4)') j, mpg%fcs_wall_ind(j,1),
     &                      mpg%fcs_wall_ind(j,2), mpg%fcs_wall_type(j)
        write(*,*) mpg%fcs_wall(mpg%fcs_wall_ind(j,1):
     &                          mpg%fcs_wall_ind(j,1)+
     &                          mpg%fcs_wall_ind(j,2)-1)
      enddo
      return
      end subroutine init_geometry

***********************************************************************

      subroutine calc_conn (mpg, gm)
      use b2us_map
      implicit none
      type (mapping), intent (in) :: mpg
      type (geometry), intent (inout) :: gm
      integer :: iCv, iCv1, iCv2, ic
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
         ! at one of the targets
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
     &        iFs.eq.mpg%iFssep .or.
     &       (iFs.eq.mpg%iFssep2.and.mpg%iFssep2.gt.0)) then
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

      subroutine init_walls_and_targets( mpg, gm )
      use b2us_map
      implicit none
      type (mapping), intent (inout) :: mpg
      type (geometry), intent (in) :: gm
      integer :: i, j, k, l, indd, inew,
     *           iCv, iCv1, iCv2, iCn, iFc, iFcc, iFs, iVx, iVa, iVb
      integer :: icount, fccount, artificial_slab = 0, found = 0
      integer, allocatable :: indXpt(:), indFc(:), regcount(:), freg(:),
     *                        old_face_list(:), new_face_list(:)
      integer :: fcsb(maxval(mpg%fcLbl),mpg%nCi),
     * fcsbcounter(maxval(mpg%fcLbl)), cvsb(maxval(mpg%fcLbl),mpg%nCi)
      logical :: duplicate, left_match_found, right_match_found
      external ipgeti, xerrab

*   ..identify periodicity
      call is_periodic( mpg, gm )

*   ..identify cells on closed surfaces
*   ..for Goat grids, was already done in is_periodic
      if (mpg%isClassicalGrid.eq.1) then
        if (mpg%nnreg(0).eq.1) then
          mpg%cvOnClosedSurface(1:mpg%nCi) = mpg%periodic_bc.eq.1
        else if (mpg%nnreg(0).eq.2 .or. mpg%nnreg(0).eq.4 .or.
     .           mpg%nnreg(0).eq.5 .or. mpg%nnreg(0).eq.8) then
          mpg%cvOnClosedSurface(1:mpg%nCi) =
     .      mod(mpg%cvReg(1:mpg%nCi),4).eq.1
        else if (mpg%nnreg(0).eq.7) then
          mpg%cvOnClosedSurface(1:mpg%nCi) = mpg%cvReg(1:mpg%nCi).eq.1
        else
          mpg%cvOnClosedSurface(1:mpg%nCi) = mpg%cvReg(1:mpg%nCi).eq.1
          write (*,*) 'Warning: unknown number of regions'
          write (*,*) 'Assuming core is region 1'
          !call xerrab ('Unknown number of regions!')
        end if
      end if
      mpg%cvOnClosedSurface(mpg%nCi+1:mpg%nCv) = .false.

!! standard approach to find X-points: find vertices where 4 of the surrounding faces lie on a flux surface.
      allocate(indXpt(mpg%nVx))
      indXpt = 0
      mpg%nXpt = 0
      mpg%iactive = 0
      mpg%iinactive = 0
      do iVx = 1, mpg%nVx
        if (mpg%vxFcP(iVx,2).ge.4) then
          fccount = 0
          do i = mpg%vxFcP(iVx,1),
     &           mpg%vxFcP(iVx,1) + mpg%vxFcP(iVx,2) - 1
            iFc = mpg%vxFc(i)
            if (mpg%fcFs(iFc).eq.0) cycle
            fccount = fccount + 1
          end do
          if (fccount.eq.4) then
            mpg%nXpt = mpg%nXpt + 1
            indXpt(mpg%nXpt) = iVx
            do i = mpg%vxCvP(iVx,1),
     &             mpg%vxCvP(iVx,1) + mpg%vxCvP(iVx,2) - 1
              iCv = mpg%vxCv(i)
              if (mpg%cvReg(iCv).eq.1) then
                if (mpg%iactive.eq.0) then
                  mpg%iactive = mpg%nXpt
                else if (mpg%iinactive.eq.0) then
                  mpg%iinactive = mpg%nXpt
                end if
              end if
            end do
          end if
        end if
      end do

!! use a different approach for slab cases with cuts,
!! because in this case the two X-point vertices are each surrounded by 2 faces that lie on a flux surface, just like all other internal vertices.
!! An X-point vertex in this case is surrounded by cvReg values 1, 2, and 3 or 4
      call ipgeti ('artificial_slab', artificial_slab)
      if (mpg%nXpt.eq.0.and.artificial_slab.eq.1) then
        do iVx = 1, mpg%nVx
          if (found.eq.1) cycle
          if (mpg%vxCvP(iVx,2).ge.4) then
            allocate(regcount(mpg%vxCvP(iVx,2)))
            do i = mpg%vxCvP(iVx,1),
     &             mpg%vxCvP(iVx,1) + mpg%vxCvP(iVx,2) - 1
              iCv = mpg%vxCv(i)
              regcount(i-mpg%vxCvP(iVx,1)+1) = mpg%cvReg(iCv)
            end do
            if (ANY(regcount==1).and.ANY(regcount==2).and.
     &         (ANY(regcount==3).or.ANY(regcount==4))) then
              mpg%nXpt = mpg%nXpt + 1
              indXpt(mpg%nXpt) = iVx
              found = 1
              if (mpg%iactive.eq.0) then
                do i = mpg%vxCvP(iVx,1),
     &                 mpg%vxCvP(iVx,1) + mpg%vxCvP(iVx,2) - 1
                  if (mpg%iactive.ne.0) cycle
                  iCv = mpg%vxCv(i)
                  if (mpg%cvReg(iCv).eq.1) mpg%iactive = mpg%nXpt
                end do
              end if
            end if
            deallocate(regcount)
          end if
        end do
      endif

      mpg%nStr = 0
      mpg%nTgc = 0
      if (mpg%nXpt.gt.0) then
        allocate(mpg%Xpt(mpg%nXpt))
        mpg%Xpt(1:mpg%nXpt) = indXpt(1:mpg%nXpt)
        if (mpg%iactive.ne.0) then
          do i = 1, mpg%vxFcP(mpg%Xpt(mpg%iactive),2)
            iFs =
     .       mpg%fcFs(mpg%vxFc(mpg%vxFcP(mpg%Xpt(mpg%iactive),1)+i-1))
            if (iFs.eq.0) cycle
            if (mpg%iFssep.eq.0) mpg%iFssep = iFs
            if (mpg%iFssep2.eq.0.and.iFs.ne.mpg%iFssep)
     >       mpg%iFssep2 = iFs
          end do
        end if
        indXpt = 0
        icount = 0
        allocate(mpg%strVxP(mpg%nXpt,2))
        mpg%strVxP = 0
        mpg%strVxP(1,1) = 1
        do i = 1, mpg%nXpt
          if (i.gt.1)
     >     mpg%strVxP(i,1) = mpg%strVxP(i-1,1) + mpg%strVxP(i-1,2)
          if (i.eq.mpg%iactive) then
            iFs = mpg%iFssep
          else
            iFs = mpg%vxFs(mpg%Xpt(i))
          end if
          do j = mpg%fsFcP(iFs,1),
     >           mpg%fsFcP(iFs,1) + mpg%fsFcP(iFs,2) - 1
            iFc = mpg%fsFc(j)
            do k = 1, 2
              iVx = mpg%fcVx(iFc,k)
              fccount = 0
              do l = mpg%vxFcP(iVx,1),
     >               mpg%vxFcP(iVx,1) + mpg%vxFcP(iVx,2) - 1
                if (mpg%fcFs(mpg%vxFc(l)).eq.iFs) fccount = fccount + 1
              end do
              if (fccount.eq.1) then
                icount = icount + 1
                mpg%strVxP(i,2) = mpg%strVxP(i,2) + 1
                indXpt(icount) = iVx
              end if
            end do
          end do
          if (mpg%iFssep2.gt.0) then
            if (i.eq.mpg%iinactive .or. i.eq.mpg%iactive) then
              iFs = mpg%iFssep2
              do j = mpg%fsFcP(iFs,1),
     &               mpg%fsFcP(iFs,1) + mpg%fsFcP(iFs,2) - 1
                iFc = mpg%fsFc(j)
                do k = 1, 2
                  iVx = mpg%fcVx(iFc,k)
                  fccount = 0
                  do l = mpg%vxFcP(iVx,1),
     &                   mpg%vxFcP(iVx,1) + mpg%vxFcP(iVx,2) - 1
                    if (mpg%fcFs(mpg%vxFc(l)).eq.iFs)
     >               fccount = fccount + 1
                  end do
                  if (fccount.eq.1) then
                    icount = icount + 1
                    mpg%strVxP(i,2) = mpg%strVxP(i,2) + 1
                    indXpt(icount) = iVx
                  end if
                end do
              end do
            end if
          end if
        end do
        mpg%nStr = icount
        allocate(mpg%strVx(mpg%nStr))
        mpg%strVx(1:mpg%nStr) = indXpt(1:mpg%nStr)
        allocate(mpg%strDiv(mpg%nStr))
        allocate(freg(mpg%nStr))
        mpg%strDiv = 0
        freg = 0
        do i = 1, mpg%nXpt
          do j = mpg%strVxP(i,1),
     &           mpg%strVxP(i,1) + mpg%strVxP(i,2) - 1
            iVx = mpg%strVx(j)
            do k = mpg%vxFcP(iVx,1),
     &             mpg%vxFcP(iVx,1) + mpg%vxFcP(iVx,2) - 1
              if (mpg%strDiv(j).ne.0) cycle
              iFc = mpg%vxFc(k)
              if (mpg%fcReg(iFc).eq.1) then
                mpg%strDiv(j) = 1
                freg(j) = 1
              else if (mpg%fcReg(iFc).eq.4) then
                mpg%strDiv(j) = 2
                freg(j) = 4
              else if (mpg%fcReg(iFc).eq.5) then
                if (mpg%nXpt.eq.1) then
                  mpg%strDiv(j) = 1
                  freg(j) = 5
                else if (i.eq.mpg%iactive .or. i.eq.mpg%iinactive) then
                  mpg%strDiv(j) = 3
                  freg(j) = 5
                endif
              else if (mpg%fcReg(iFc).eq.8) then
                if (mpg%nXpt.eq.1) then
                  mpg%strDiv(j) = 2
                  freg(j) = 8
                else if (i.eq.mpg%iactive .or. i.eq.mpg%iinactive) then
                  mpg%strDiv(j) = 4
                  freg(j) = 8
                endif
              end if
            end do
          end do
        end do
        allocate(mpg%tgVx(mpg%nTgc)) !csc needed for adjoint AD
      else
        allocate(freg(2))
        iFc = 1
        icount = 0
        do while (mpg%iFssep.eq.0.and.iFc.le.mpg%nFc)
          iCv1 = mpg%fcCv(iFc,1)
          iCv2 = mpg%fcCv(iFc,2)
          if (mpg%cvReg(iCv1).ne.mpg%cvReg(iCv2).and.
     &        mpg%cvReg(iCv1).ne.0.and.mpg%cvReg(iCv2).ne.0) then
            mpg%iFssep = mpg%fcFs(iFc)
          end if
          if (mpg%iFssep.eq.0) iFc = iFc + 1
        end do
        if (mpg%iFssep.gt.0) then
          do i = mpg%fsFcP(mpg%iFssep,1),
     &           mpg%fsFcP(mpg%iFssep,1) + mpg%fsFcP(mpg%iFssep,2) - 1
            iFc = mpg%fsFc(i)
            do k = 1, 2
              iVx = mpg%fcVx(iFc,k)
              do j = mpg%vxFcP(iVx,1),
     &               mpg%vxFcP(iVx,1) + mpg%vxFcP(iVx,2) - 1
                iFcc = mpg%vxFc(j)
                if (mpg%fcLbl(iFcc).ne.0) then
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
          mpg%nTgc = icount
          allocate(mpg%tgVx(mpg%nTgc))
          mpg%tgVx(1:mpg%nTgc) = indXpt(1:mpg%nTgc)
          allocate(mpg%strDiv(2))
          mpg%strDiv(1) = 1
          mpg%strDiv(2) = 2
          freg(1) = 1
          freg(2) = 2
        end if
      end if

      if (mpg%nnreg(0).eq.1 .and. mpg%iFssep.eq.0) then
        do iCv = 1, mpg%nCi
          if (mpg%iFssep.gt.0) cycle
          do j = mpg%cvFcP(iCv,1),
     &           mpg%cvFcP(iCv,1) + mpg%cvFcP(iCv,2) - 1
            iFc = mpg%cvFc(j)
            if (mpg%fcCv(iFc,1).eq.iCv) iCn = mpg%fcCv(iFc,2)
            if (mpg%fcCv(iFc,2).eq.iCv) iCn = mpg%fcCv(iFc,1)
            if (iCn.gt.mpg%nCi) cycle
            if
     .     (mpg%cvOnClosedSurface(iCv).neqv.mpg%cvOnClosedSurface(iCn))
     >       mpg%iFssep = mpg%fcFs(iFc)
          end do
        end do
        if (mpg%iFssep.eq.0) call ipgeti ('b2mwti_jsep', mpg%iFssep)
        if (.not.allocated(freg)) allocate(freg(2))
        freg(1) = 1
        freg(2) = 2
        if (mpg%iFssep.gt.0) then
          indXpt = 0
          icount = 0
          allocate(mpg%strVxP(1,2))
          mpg%strVxP(1,1) = 1
          mpg%strVxP(1,2) = 2
          do j = mpg%fsFcP(mpg%iFssep,1),
     &           mpg%fsFcP(mpg%iFssep,1) + mpg%fsFcP(mpg%iFssep,2) - 1
            iFc = mpg%fsFc(j)
            do k = 1, 2
              iVx = mpg%fcVx(iFc,k)
              fccount = 0
              do l = mpg%vxFcP(iVx,1),
     &               mpg%vxFcP(iVx,1) + mpg%vxFcP(iVx,2) - 1
                if (mpg%fcFs(mpg%vxFc(l)).eq.mpg%iFssep)
     >           fccount = fccount + 1
              end do
              if (fccount.eq.1) then
                icount = icount + 1
                indXpt(icount) = iVx
              end if
            end do
          end do
          mpg%nStr = icount
          if (mpg%nStr.gt.0) then
            allocate(mpg%strVx(mpg%nStr))
            mpg%strVx(1:mpg%nStr) = indXpt(1:mpg%nStr)
            allocate(mpg%strDiv(mpg%nStr))
            mpg%strDiv = 0
            do i = 1, mpg%nStr
              iVx = mpg%strVx(i)
              do j = mpg%vxFcP(iVx,1),
     &               mpg%vxFcP(iVx,1) + mpg%vxFcP(iVx,2) - 1
                if (mpg%strDiv(i).ne.0) cycle
                iFc = mpg%vxFc(j)
                if (mpg%fcReg(iFc).eq.1) then
                  mpg%strDiv(i) = 1
                else if (mpg%fcReg(iFc).eq.2) then
                  mpg%strDiv(i) = 2
                end if
              end do
            end do
          end if
        end if
      end if
      deallocate(indXpt)

      if (allocated(mpg%strDiv)) then
       indd = size(mpg%strDiv)
       if (indd.gt.0 .and. maxval(mpg%strDiv).gt.0) then
        allocate(mpg%ifdiv(maxval(mpg%strDiv)))
        allocate(mpg%ivdiv(maxval(mpg%strDiv)))
        allocate(mpg%divFcP(maxval(mpg%strDiv),2))
        allocate(indFc(mpg%nFc))
        indFc = 0
        icount = 0
        mpg%ifdiv = 0
        mpg%ivdiv = 0
        mpg%divFcP(:,2) = 0
        do i = 1, indd  ! first assemble all the faces
          if (mpg%strDiv(i).eq.0 .or. freg(i).eq.0) cycle
          if (mpg%divFcP(mpg%strDiv(i),2).ne.0) cycle ! already done, skip
          do iFc = 1, mpg%nFc
            if (mpg%fcReg(iFc).eq.freg(i)) then
              icount = icount + 1
              mpg%divFcP(mpg%strDiv(i),2) =
     &         mpg%divFcP(mpg%strDiv(i),2) + 1
              indFc(icount) = iFc
            end if
          end do
        end do
        mpg%divFcP(1,1) = 1
        do i = 2, maxval(mpg%strDiv)
          mpg%divFcP(i,1) = mpg%divFcP(i-1,1) + mpg%divFcP(i-1,2)
        end do
        icount = count(indFc(:).ne.0)
        allocate(mpg%divFc(icount))
        allocate(mpg%divFcOr(icount))
        mpg%divFc = 0
        mpg%divFcOr = 0.0_R8
        icount = 1
        do i = 1, indd
          if (mpg%strDiv(i).eq.0) cycle
          if (mpg%divFc(mpg%divFcP(mpg%strDiv(i),1)).ne.0) cycle ! already done
          mpg%divFc(mpg%divFcP(mpg%strDiv(i),1):
     .              mpg%divFcP(mpg%strDiv(i),1)+
     .              mpg%divFcP(mpg%strDiv(i),2)-1) =
     .      indFc(icount:icount+mpg%divFcP(mpg%strDiv(i),2)-1)
          icount = icount + mpg%divFcP(mpg%strDiv(i),2)
        end do
        do i = 1, maxval(mpg%strDiv) ! now re-order them
          allocate(old_face_list(mpg%divFcP(i,2)))
          allocate(new_face_list(mpg%divFcP(i,2)))
          old_face_list(1:mpg%divFcP(i,2))=
     &     mpg%divFc(mpg%divFcP(i,1):mpg%divFcP(i,1)+mpg%divFcP(i,2)-1)
          new_face_list(:) = 0
          new_face_list(1) = mpg%divFc(mpg%divFcP(i,1))
          old_face_list(1) = 0
          iVa = mpg%fcVx(mpg%divFc(mpg%divFcP(i,1)),1)
          iVb = mpg%fcVx(mpg%divFc(mpg%divFcP(i,1)),2)
          inew = 1
          do while (inew.lt.mpg%divFcP(i,2))
            left_match_found = .false.
            do j = 1, mpg%divFcP(i,2)
              if (left_match_found) cycle
              if (old_face_list(j).eq.0) cycle
              if (iVa.eq.mpg%fcVx(old_face_list(j),1)) then
                left_match_found = .true.
                iVa = mpg%fcVx(old_face_list(j),2)
              else if (iVa.eq.mpg%fcVx(old_face_list(j),2)) then
                left_match_found = .true.
                iVa = mpg%fcVx(old_face_list(j),1)
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
            do j = 1, mpg%divFcP(i,2)
              if (right_match_found) cycle
              if (old_face_list(j).eq.0) cycle
              if (iVb.eq.mpg%fcVx(old_face_list(j),1)) then
                right_match_found = .true.
                iVb = mpg%fcVx(old_face_list(j),2)
              else if (iVb.eq.mpg%fcVx(old_face_list(j),2)) then
                right_match_found = .true.
                iVb = mpg%fcVx(old_face_list(j),1)
              end if
              if (right_match_found) then
                new_face_list(inew+1) = old_face_list(j)
                old_face_list(j) = 0
                inew = inew + 1
              end if
            end do
            if (.not.(left_match_found.or.right_match_found).and.
     &          inew.lt.mpg%divFcP(i,2)) then
              do j = 1, inew
                iFc = new_face_list(j)
                iVa = mpg%fcVx(iFc,1)
                iVb = mpg%fcVx(iFc,2)
                write(6,'(a,4i6)')
     &        'Chain of faces and vertices in new list ',j,iFc,iVa,iVb
              end do
              do j = 1, mpg%divFcP(i,2)
                if (old_face_list(j).eq.0) cycle
                iFc = mpg%divFc(mpg%divFcP(i,1)+j-1)
                iVa = mpg%fcVx(iFc,1)
                iVb = mpg%fcVx(iFc,2)
                write(6,'(a,4i6)')
     &        'Left over faces and vertices in old list',j,iFc,iVa,iVb
              end do
              call xerrab ('Broken face chain along target')
            end if
          end do
          mpg%divFc(mpg%divFcP(i,1):mpg%divFcP(i,1)+mpg%divFcP(i,2)-1)=
     &     new_face_list(1:mpg%divFcP(i,2))
          deallocate(old_face_list)
          deallocate(new_face_list)
        end do
        deallocate(indFc)
       end if
      end if

      if (allocated(freg)) then
        deallocate(freg)
      else
        write(*,*) 'Targetless topology detected :'
        write(*,*) 'Strike point data not defined'
      end if

      allocate(mpg%cvAlongBoundary(mpg%nCv))
      mpg%cvAlongBoundary = .false.
      fcsbcounter = 0
      do iCv = mpg%nCi+1, mpg%nCv
        iFc = mpg%cvFc(mpg%cvFcP(iCv,1))
        if (mpg%fcLbl(iFc).le.0) cycle
        fcsbcounter(mpg%fcLbl(iFc)) = fcsbcounter(mpg%fcLbl(iFc))+1
        fcsb(mpg%fcLbl(iFc),fcsbcounter(mpg%fcLbl(iFc))) = iFc
        iCv1 = mpg%fcCv(iFc,1)
        iCv2 = mpg%fcCv(iFc,2)
        cvsb(mpg%fcLbl(iFc),fcsbcounter(mpg%fcLbl(iFc))) =
     .   min(iCv1,iCv2)
      enddo
      do k=1,maxval(mpg%fcLbl)
        do j=1,fcsbcounter(k)
          iCv2 = cvsb(k,j)
          iFc = fcsb(k,j)
          iCv1 = maxval(mpg%fcCv(iFc,1:2))
          if (j.lt.fcsbcounter(k)) then
            if (mpg%cvFt(iCv2).eq.mpg%cvFt(cvsb(k,j+1))) then
              mpg%cvAlongBoundary(iCv2) = .true.
            endif
          endif
          if (j.gt.1) then
            if (mpg%cvFt(iCv2).eq.mpg%cvFt(cvsb(k,j-1))) then
              mpg%cvAlongBoundary(iCv2) = .true.
            endif
          endif
        enddo
      enddo

      allocate(mpg%fcs_wall(mpg%nCg))

      k = 1
      do l = mpg%nCi+1, mpg%nCv
        mpg%fcs_wall(k) = mpg%cvFc(mpg%cvFcP(l,1))
        k = k+1
      enddo

      do k = 1, mpg%nCg
        do l = k+1, mpg%nCg
          if
     &     (mpg%fcVx(mpg%fcs_wall(k),1) == mpg%fcVx(mpg%fcs_wall(l),1)
     &      .or.
     &      mpg%fcVx(mpg%fcs_wall(k),1) == mpg%fcVx(mpg%fcs_wall(l),2)
     &      .or.
     &      mpg%fcVx(mpg%fcs_wall(k),2) == mpg%fcVx(mpg%fcs_wall(l),1)
     &      .or.
     &      mpg%fcVx(mpg%fcs_wall(k),2) == mpg%fcVx(mpg%fcs_wall(l),2))
     >     then
            inew = mpg%fcs_wall(k+1)
            mpg%fcs_wall(k+1) = mpg%fcs_wall(l)
            mpg%fcs_wall(l) = inew
            exit
          endif
        enddo
      enddo

      write (*,*) 'nXpt, nStr, nTgc ', mpg%nXpt, mpg%nStr, mpg%nTgc

      return
      end subroutine init_walls_and_targets

**********************************************************************

      subroutine is_periodic ( mpg, gm )
      use b2us_map
      implicit none
      type(mapping), intent(inout) :: mpg
      type(geometry), intent(in) :: gm
      integer :: iCvm, iCvn, iCvp, iCv, iFc, iFt, icount, keep_count
      integer :: i, iCvtmp, iVc1, iVc2
      logical :: closed, end_found, all_closed, pinch

      mpg%periodic_bc = 0
      if (mpg%isClassicalGrid.eq.1 .and.
     &   (mpg%nnreg(0).eq.2 .or. mpg%nnreg(0).eq.5)) then
        mpg%periodic_bc = 1
      else if (mpg%nnreg(0).eq.1 .or. mpg%isClassicalGrid.eq.0) then
        all_closed = mpg%nFt.gt.0
        do iFt = 1, mpg%nFt
          closed = .false.
          iCv = mpg%ftCv(mpg%ftCvP(iFt,1))
          if (mpg%ftCvP(iFt,2).gt.1) then
            iCvm = iCv
            iCvn = 0
            keep_count = 0
            do iFc = mpg%cvFcP(iCv,1),
     &               mpg%cvFcP(iCv,1) + mpg%cvFcP(iCv,2) - 1
              if (gm%fcQalf(mpg%cvFc(iFc),0).eq.0.0_R8) cycle
              if (mpg%fcCv(mpg%cvFc(iFc),1).eq.iCv) then
                iCvtmp = mpg%fcCv(mpg%cvFc(iFc),2)
              else
                iCvtmp = mpg%fcCv(mpg%cvFc(iFc),1)
              end if
              if (iCvtmp.gt.mpg%nCi) cycle
              keep_count = keep_count + 1
              iCvn = iCvtmp
            end do
            closed = iCvn.eq.0
            end_found = iCvn.eq.0
            do while (.not.end_found)
              iCvp = 0
              keep_count = 0
              do iFc = mpg%cvFcP(iCvn,1),
     &                 mpg%cvFcP(iCvn,1) + mpg%cvFcP(iCvn,2) - 1
                if (gm%fcQalf(mpg%cvFc(iFc),0).eq.0.0_R8) cycle
                if (mpg%fcCv(mpg%cvFc(iFc),1).eq.iCvm) cycle
                if (mpg%fcCv(mpg%cvFc(iFc),2).eq.iCvm) cycle
                if (mpg%fcCv(mpg%cvFc(iFc),1).eq.iCvn) then
                  iCvtmp = mpg%fcCv(mpg%cvFc(iFc),2)
                else
                  iCvtmp = mpg%fcCv(mpg%cvFc(iFc),1)
                end if
                if (iCvtmp.gt.mpg%nCi) cycle
                keep_count = keep_count + 1
                iCvp = iCvtmp
              end do
              if (iCvp.eq.0) then
                ! genuine dead end: check whether the cell we are
                ! stuck at shares a vertex with the flux tube
                ! original seed cell -- a valid closed-loop pinch-point
                ! (e.g. limiter-tangent) with no direct face link
                pinch = .false.
                do iVc1 = mpg%cvVxP(iCvn,1),
     &                    mpg%cvVxP(iCvn,1)+mpg%cvVxP(iCvn,2)-1
                  do iVc2 = mpg%cvVxP(iCv,1),
     &                      mpg%cvVxP(iCv,1)+mpg%cvVxP(iCv,2)-1
                    if (mpg%cvVx(iVc1).eq.mpg%cvVx(iVc2))
     &                pinch = .true.
                  end do
                end do
                closed = pinch
                end_found = .true.
              else
                iCvm = iCvn
                iCvn = iCvp
                closed = iCvp.eq.iCv
                end_found = closed
              end if
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
          if (closed .and. mpg%isClassicalGrid.eq.0) then
            do i = 1, mpg%ftCvP(iFt,2)
              iCv = mpg%ftCv(mpg%ftCvP(iFt,1)+i-1)
              mpg%cvOnClosedSurface(iCv) = .true.
            end do
          end if
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
         found_next = .false.
         do j = mpg%vxFcP(vx_current,1),
     &          mpg%vxFcP(vx_current,1)+mpg%vxFcP(vx_current,2)-1
           fc_next = mpg%vxFc(j)
           if (fc_next==fc_current) cycle  ! we ignore the current face
           ! check if face belongs to same flux tube
           if (mpg%fcFs(fc_next)==ifs) then
             ! if we are back at the starting face, this
             ! flux surface is closed: stop here instead
             ! of looping forever
             if (fc_next==ifc) then
               found_next = .false.
               exit
             end if
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
