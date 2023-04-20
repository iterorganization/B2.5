      module b2us_geo

      use b2mod_types
      implicit none
      
      private
      public :: alloc_geometry, dealloc_geometry,
     .          read_geometry, write_geometry,
     .          init_geometry, calc_conn, check_geometry

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

!  flux tube/surface quantities
        real (kind=R8), allocatable  :: cvConn(:)      ! (nCv) Connection length for each CV
        real (kind=R8), allocatable  :: ftConn(:)      ! (nFt) Connection length for each flux tube
        real (kind=R8), allocatable  :: ftEps(:)       ! (nFt) Inverse aspect ratio
        real (kind=R8), allocatable  :: ftBbav2(:)     ! (nFt) Average B**2 on the flux tube

        real (kind=R8), allocatable  :: fsPsi(:)       ! (nFs) Psi on each flux surface

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

!  flux tube quantites
      allocate (g%cvConn(nCv))
      allocate (g%ftConn(nFt))
      allocate (g%ftEps(nFt))
      allocate (g%ftBbav2(nFt))

!  flux surface quantities
      allocate (g%fsPsi(nFs))

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

!  flux surface quantites
      g%cvConn = 0._R8
      g%ftConn = 0._R8
      g%ftEps = 0._R8
      g%ftBbav2 = 0._R8

      g%fsPsi = 0._R8
      
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

!  flux surface quantites
      deallocate (g%cvConn)
      deallocate (g%ftConn)
      deallocate (g%ftEps)
      deallocate (g%ftBbav2)
      deallocate (g%fsPsi)
      
      return
      end subroutine dealloc_geometry

***********************************************************************


      subroutine read_geometry (iun, nCv, nFc, nVx, nFs, nFt, nVmxCv,
     .  gm)

      implicit none
      integer, intent(in) :: iun, nCv, nFc, nVx, nFs, nFt, nVmxCv
      integer :: idum(0:9)
      type (geometry), intent (inout) :: gm

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
      integer :: idum(0:9)
      type (geometry), intent (in) :: gm

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
      type (mapping), intent (in) :: mpg
      type (geometry), intent (inout) :: gm
      integer :: nCv, nFc, nVx, nFx, iCv, iFc, iVx, inCv, i1, i2,
     . inv_dist(mpg%nVx) 
      real (kind=R8) :: hzconst, fac_hz = 0.0_R8, dux, duy, du, width,
     .  sbf, cvPsi(mpg%nCv), count_up, count_down, dpsi(mpg%nCv),
     .  dpsi_max_down, dpsi_max_up
      integer, allocatable :: verts(:)
      intrinsic maxval, maxloc, minloc, abs, sqrt
      data width /1.0_R8/
      save width

      call ipgetr ('b2mndr_hz', fac_hz)
      call xertst (0.0.le.fac_hz .and. 1.0.ge.fac_hz,
     &  'faulty parameter b2mndr_hz')
      call ipgetr ('b2agmt_1d_width', width)

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
         !calculate cvPsi
         do iCv = 1, nCv
            allocate(verts(1:mpg%cvVxP(iCv,2)))
            verts = mpg%cvVx(mpg%cvVxP(iCv,1):mpg%cvVxP(iCv,1)+
     &         mpg%cvVxP(iCv,2))
            cvPsi(iCv) = sum(gm%vxFpsi(verts))/mpg%cvVxP(iCv,2)
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
               if ( cvPsi(iCv) .gt. gm%vxFpsi(iVx)) then
                  count_up = count_up + 1.0e0_R8
               elseif ( cvPsi(iCv) .lt. gm%vxFpsi(iVx)) then
                  count_down = count_down + 1.0e0_R8
               endif 
            enddo
            do inCv = 1, mpg%vxCvP(iVx,2)
              iCv = mpg%vxCv(mpg%vxCvP(iVx,1)+ inCv - 1)
              if (iCv.le.mpg%nCi) then
               if ( cvPsi(iCv) .gt. gm%vxFpsi(iVx)) then
                   gm%vxVol(mpg%vxCvP(iVx,1)+inCv - 1) = 
     &               count_up / mpg%vxCvP(iVx,2)               
               elseif ( cvPsi(iCv) .lt. gm%vxFpsi(iVx)) then
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
         !calculate cvPsi
         do iCv = 1, nCv
            allocate(verts(1:mpg%cvVxP(iCv,2)))
            verts = mpg%cvVx(mpg%cvVxP(iCv,1):mpg%cvVxP(iCv,1)+
     &         mpg%cvVxP(iCv,2))
            cvPsi(iCv) = sum(gm%vxFpsi(verts))/mpg%cvVxP(iCv,2)
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
               if ( cvPsi(iCv) .gt. gm%vxFpsi(iVx)) then
                  count_up = count_up + 1.0e0_R8
                  if ( dpsi(iCv) .gt. dpsi_max_up) then
                     dpsi_max_up = dpsi(iCv)
                  endif
               elseif ( cvPsi(iCv) .lt. gm%vxFpsi(iVx)) then
                  count_down = count_down + 1.0e0_R8
                  if ( dpsi(iCv) .gt. dpsi_max_down) then
                     dpsi_max_down = dpsi(iCv)
                  endif
               endif 
            enddo
            do inCv = 1, mpg%vxCvP(iVx,2)
              iCv = mpg%vxCv(mpg%vxCvP(iVx,1)+ inCv - 1)
              if (iCv.le.mpg%nCi) then
               if ( cvPsi(iCv) .gt. gm%vxFpsi(iVx)) then
                 if (mpg%cvVxP(iCv,2) .eq. 3) then
                   gm%vxVol(mpg%vxCvP(iVx,1)+inCv - 1) = 
     &               count_up / (mpg%vxCvP(iVx,2)*sqrt(2.0_R8)* 
     &                sqrt(dpsi_max_up/dpsi(iCv)))               
                 else
                   gm%vxVol(mpg%vxCvP(iVx,1)+inCv - 1) = 
     &               count_up / mpg%vxCvP(iVx,2)
                 endif 
               elseif ( cvPsi(iCv) .lt. gm%vxFpsi(iVx)) then
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
      gm%vxOnedBsq = 1.0_R8/gm%vxBb(:,3)**2

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


*   ..compute cvHz
      hzconst = 2.0_R8*pi*6.0_R8   ! for ITER for testing
      hzconst = 2.0_R8*pi*1.6_R8   ! for ASDEX-UPGRADE for testing
      hzconst = 1.0_R8
      gm%cvHz = hzconst*(1.0e0_R8-fac_hz)+fac_hz*gm%cvVol/gm%cvSz

*   ..compute fcHz and fcPbshz
      if (switch%fcHz_style.eq.0) then
        call intface(nCv,nFc,mpg%fcCv,gm%fcVol,gm%cvHz,gm%fcHz)
      else
        if (isymm.eq.0) then
          gm%fcHz = hzconst*(1.0e0_R8-fac_hz)+fac_hz*width
        elseif (isymm.eq.1 .or. isymm.eq.2) then
          do iFc = 1, mpg%nFc
            gm%fcHz(iFc) = hzconst*(1.0e0_R8-fac_hz)+fac_hz*
     &        pi*(gm%vxX(mpg%fcVx(iFc,1))+gm%vxX(mpg%fcVx(iFc,2)))
          enddo
        else if (isymm.eq.3 .or. isymm.eq.4) then
          do iFc = 1, mpg%nFc
            gm%fcHz(iFc) = hzconst*(1.0e0_R8-fac_hz)+fac_hz*
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
          gm%vxHz = hzconst*(1.0e0_R8-fac_hz)+fac_hz*width
        elseif (isymm.eq.1 .or. isymm.eq.2) then
          gm%vxHz = hzconst*(1.0e0_R8-fac_hz)+
     &      fac_hz*2.0_R8*pi*gm%vxX
        else if (isymm.eq.3 .or. isymm.eq.4) then
          gm%vxHz = hzconst*(1.0e0_R8-fac_hz)+
     &      fac_hz*2.0_R8*pi*gm%vxY
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

      return
      end subroutine init_geometry


***********************************************************************

      subroutine calc_conn (mpg, gm)
      use b2us_map
      implicit none
      type (mapping), intent (in) :: mpg
      type (geometry), intent (inout) :: gm
      integer :: nfx, iCv, iCv1, iCv2, iFt, i, j, ic, iFc
      integer, save :: conn_style = 0, iout = 0
      real (kind=r8) :: wrk0(mpg%nCv), wrk1(mpg%nCv), Lr, Ll,
     .  tt0(mpg%nFt), tt1(mpg%nFt), tt2(mpg%nFt), tt3(mpg%nFt)

      call ipgeti ('b2mndr_conn_style',conn_style)
      call ipgeti ('b2mndr_conn_iout',iout)

      if (mpg%isClassicalGrid.eq.1) then
       ! reconstruct cvHx into wrk1
       nfx = maxval(mpg%imapFcx)
       wrk1 = 0.0_R8
       do ifc = 1, nfx
         icv1 = mpg%fcCv(ifc,1)
         icv2 = mpg%fcCv(ifc,2)
         wrk1(icv1) = wrk1(icv1) + gm%fcHc(ifc,1)
         wrk1(icv2) = wrk1(icv2) + gm%fcHc(ifc,2)
       enddo
       ! 'N/S guard cells': take length of face
       do ifc = nfx+1, mpg%nFc
         icv = max(mpg%fcCv(iFc,1), mpg%fcCv(iFc,2))
         if (icv.gt.mpg%nCi) wrk1(icv) = gm%fcHt(ifc)
       enddo
       ! construct wrk arrays
       wrk0 = 0.0_R8
       tt0 = 0.0_R8
       tt1 = 0.0_R8
       tt2 = 0.0_R8
       tt3 = 0.0_R8
       do iFt = 1, mpg%nFt 
         do j = mpg%ftCvP(iFt,1), mpg%ftCvP(iFt,1)+mpg%ftCvP(iFt,2)-1
           icv = mpg%ftCv(j)
           wrk0(icv) = wrk1(icv)*gm%cvBb(icv,3)/gm%cvBb(icv,0)
           tt0(iFt) = tt0(iFt) + wrk1(icv)
           tt1(iFt) = tt1(iFt) + wrk1(icv)*gm%cvVol(icv)/gm%cvSz(icv)
           tt2(iFt) = tt2(iFt) + gm%cvBb(iCv,3)*gm%cvBb(iCv,3)*
     .         gm%cvVol(iCv)
           tt3(iFt) = tt3(iFt) + gm%cvVol(iCv)
         enddo
       enddo
       ! full connection length of each flux tube
       gm%ftConn = 0.0_R8
       gm%ftEps = 0.0_R8
       gm%ftBbav2 = 0.0_R8
       do iFt = 1, mpg%nFt 
         Ll = 0.0_R8
         do j=mpg%ftCvP(iFt,1), mpg%ftCvP(iFt,1)+mpg%ftCvP(iFt,2)-1
          Ll = Ll + wrk0(mpg%ftCv(j))
         enddo
         gm%ftConn(iFt) = Ll
         if (mpg%cvOnClosedSurface(mpg%ftCv(mpg%ftCvP(iFt,1)))) then
           gm%ftEps(iFt) = tt0(iFt)**2/tt1(iFt)
           gm%ftBbav2(iFt) = tt2(iFt)/tt3(iFt)
         endif
       enddo
       ! 'effective' connection length for use in k-eps model
       if (conn_style.ne.0) then ! if 0, assumed to be read from b2fgmtry
        gm%cvConn = 0.0_R8
        do iFt = 1, mpg%nFt 
         iCv = mpg%ftCv(mpg%ftCvP(iFt,1))
         if (mpg%cvOnClosedSurface(iCv)) then
          !core region: actual connection length (same as ftConn)
          gm%cvConn(mpg%ftCv(mpg%ftCvP(iFt,1):mpg%ftCvP(iFt,1)+
     &                         mpg%ftCvP(iFt,2)-1)) = gm%ftConn(iFt)
         else
          ! SOL/PF: effective distance to targets/boundaries
          do ic = mpg%ftCvP(iFt,1),mpg%ftCvP(iFt,1)+mpg%ftCvP(iFt,2)-1
           iCv = mpg%ftCv(ic)
           Lr = 0.5_R8*wrk0(iCv)
           do j = ic+1, mpg%ftCvP(iFt,1) + mpg%ftCvP(iFt,2)-1
            Lr = Lr + wrk0(mpg%ftCv(j))
           end do
           Ll = 0.5_R8*wrk0(iCv)
           do j = mpg%ftCvP(iFt,1), ic-1
            Ll = Ll + wrk0(mpg%ftCv(j))
           end do
           if (conn_style.eq.1) then
            gm%cvConn(iCv) = min(Lr,Ll)
           elseif (conn_style.eq.2) then
            gm%cvConn(iCv) = Lr*Ll/(Lr+Ll)
           else
            write (*,*) 'b2mndr_conn_style=',conn_style
            call xerrab ('parameter b2mndr_conn_style not coded')
           endif
          end do
         endif
        end do
       endif
      else
       write(*,*) 
     &   'WARNING: connection length not yet available for US grids'
       ! reconstruct cvHx into wrk1
       do ifc = 1, mpg%nFc
         if (gm%fcQalf(ifc,0) .ne. 0) then ! non-aligned face
           icv1 = mpg%fcCv(ifc,1)
           icv2 = mpg%fcCv(ifc,2)
           wrk1(icv1) = wrk1(icv1) + gm%fcHc(ifc,1)
           wrk1(icv2) = wrk1(icv2) + gm%fcHc(ifc,2)           
         else ! 'N/S guard cells': take length of face
           icv = max(mpg%fcCv(iFc,1), mpg%fcCv(iFc,2))
           if (icv.gt.mpg%nCi) wrk1(icv) = gm%fcHt(ifc)          
         endif
       enddo
       ! construct wrk arrays
       wrk0 = 0.0_R8
       tt0 = 0.0_R8
       tt1 = 0.0_R8
       tt2 = 0.0_R8
       tt3 = 0.0_R8
       do iFt = 1, mpg%nFt 
         do j = mpg%ftCvP(iFt,1), mpg%ftCvP(iFt,1)+mpg%ftCvP(iFt,2)-1
           icv = mpg%ftCv(j)
           wrk0(icv) = wrk1(icv)*gm%cvBb(icv,3)/gm%cvBb(icv,0)
           tt0(iFt) = tt0(iFt) + wrk1(icv)
           tt1(iFt) = tt1(iFt) + wrk1(icv)*gm%cvVol(icv)/gm%cvSz(icv)
           tt2(iFt) = tt2(iFt) + gm%cvBb(iCv,3)*gm%cvBb(iCv,3)*
     .         gm%cvVol(iCv)
           tt3(iFt) = tt3(iFt) + gm%cvVol(iCv)
         enddo
       enddo
       ! full connection length of considered flux tubes
       gm%ftConn = 0.0_R8
       gm%ftEps = 0.0_R8
       gm%ftBbav2 = 0.0_R8
       do iFt = 1, mpg%nFt 
         Ll = 0.0_R8
         do j=mpg%ftCvP(iFt,1), mpg%ftCvP(iFt,1)+mpg%ftCvP(iFt,2)-1
          Ll = Ll + wrk0(mpg%ftCv(j))
         enddo
         gm%ftConn(iFt) = Ll
         if (mpg%cvOnClosedSurface(mpg%ftCv(mpg%ftCvP(iFt,1)))) then
           gm%ftEps(iFt) = tt0(iFt)**2/tt1(iFt)
           gm%ftBbav2(iFt) = tt2(iFt)/tt3(iFt)
         endif
       enddo
      endif  

      if (iout.ne.0) then
       call my_out_us(70,mpg%nCv,0,gm%cvConn,'calc_conn_cvConn')
      endif

      return
      end subroutine calc_conn


***********************************************************************

      end module b2us_geo
