      module b2us_geo

      use b2mod_types
      implicit none
      
      private
      public :: alloc_geometry, dealloc_geometry

      type, public :: geometry
        
!  cell quantiites
        real (kind=R8), allocatable  :: cvBb(:,:) ! (nCv,4) Bx,By,Bz, |B|
        real (kind=R8), allocatable  :: cvX(:)
        real (kind=R8), allocatable  :: cvY(:)
        real (kind=R8), allocatable  :: cvS(:)
        real (kind=R8), allocatable  :: cvH(:)
        real (kind=R8), allocatable  :: cvHx(:)
        real (kind=R8), allocatable  :: cvQgam(:,:)
        real (kind=R8), allocatable  :: cvVol(:)
        real (kind=R8), allocatable  :: cvOnedBsq(:)
        real (kind=R8), allocatable  :: cvBzb(:)

!  face based quantities
        real (kind=R8), allocatable  :: fcBb(:,:) ! (4,nFc,4) to be recomputed at cell faces (Bx,By,Bz,|B|)        real (kind=R8), allocatable  :: crx(:,:,:)
        real (kind=R8), allocatable  :: fcS(:)
        real (kind=R8), allocatable  :: fcHc(:,:)
        real (kind=R8), allocatable  :: fcHt(:)
        real (kind=R8), allocatable  :: fcQgam(:,:)
        real (kind=R8), allocatable  :: fcQalf(:,:)
        real (kind=R8), allocatable  :: fcQbet(:,:)
        real (kind=R8), allocatable  :: fcPbs(:)
        real (kind=R8), allocatable  :: fcPbshz(:)

        real (kind=R8), allocatable  :: fcVol(:,:)

!  vertex quantities
        real (kind=R8), allocatable  :: vxX(:)
        real (kind=R8), allocatable  :: vxY(:)
        real (kind=R8), allocatable  :: vxFfbz(:)
        real (kind=R8), allocatable  :: vxFpsi(:)

        real (kind=R8), allocatable  :: fsConn(:)
      end type geometry

      contains

***********************************************************************

      subroutine alloc_geometry (g, nCv, nFc, nVx, nFs)

      implicit none
      integer, intent(in) :: nCv, nFc, nVx, nFs
      type(geometry), intent(in out) :: g

      allocate (g%cvBb(nCv,0:3)) 
      allocate (g%cvX(nCv))
      allocate (g%cvY(nCv))
      allocate (g%cvS(nCv))
      allocate (g%cvH(nCv))      ! hz
      allocate (g%cvHx(nCv))      ! hz
      allocate (g%cvQgam(nCv,0:1))
      allocate (g%cvVol(nCv))
      allocate (g%cvOnedBsq(nCv))
      allocate (g%cvBzb(nCv))

      allocate (g%fcBb(nFc,0:3)) 
      allocate (g%fcS(nFc))      !gs
      allocate (g%fcHc(nFc,1:2)) ! different definition
      allocate (g%fcHt(nFc)) ! different definition
      allocate (g%fcVol(nFc,1:2))      
      allocate (g%fcQgam(nFc,0:1))
      allocate (g%fcQalf(nFc,0:1))
      allocate (g%fcQbet(nFc,0:1))
      allocate (g%fcPbs(nFc))
      allocate (g%fcPbshz(nFc))

      allocate (g%vxX(nVx))
      allocate (g%vxY(nVx))
      allocate (g%vxFfbz(nVx))
      allocate (g%vxFpsi(nVx))
      
      allocate (g%fsConn(nFs))
      
      return
      end subroutine alloc_geometry

**********************************************************************

      subroutine dealloc_geometry (g)

      implicit none
      type(geometry), intent(in out) :: g

      deallocate (g%cvBb) 
      deallocate (g%cvX)
      deallocate (g%cvY)
      deallocate (g%cvS)
      deallocate (g%cvH)    ! hz
      deallocate (g%cvHx)      ! hz
      deallocate (g%cvQgam)
      deallocate (g%cvVol)
      deallocate (g%cvOnedBsq)
      deallocate (g%cvBzb)
      
      deallocate (g%fcBb) 
      deallocate (g%fcS)     !gs
      deallocate (g%fcHc) ! different definition
      deallocate (g%fcHt) ! different definition
      deallocate (g%fcVol)      
      deallocate (g%fcQgam)
      deallocate (g%fcQalf)
      deallocate (g%fcQbet)
      deallocate (g%fcPbs)
      deallocate (g%fcPbshz)
      
      deallocate (g%vxX)
      deallocate (g%vxY)
      deallocate (g%vxFfbz)
      deallocate (g%vxFpsi)
      
      deallocate (g%fsConn)
      
      return
      end subroutine dealloc_geometry

**********************************************************************

      end module b2us_geo
