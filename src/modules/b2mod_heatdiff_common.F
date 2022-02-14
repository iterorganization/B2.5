!> common blocks from heatdiff1D moved here to make it simpler
!> to declare them as threadprivate
      module b2mod_heatdiff_common
      use b2mod_types
      use b2mod_wall
      implicit none

      !common /zzline/ dx,x0,xn,x
      real (kind=R8) :: dx(ndepth-1), x0, xn
      real (kind=R8) :: x(ndepth)
      !$OMP THREADPRIVATE(dx,x0,xn,x)

      ! common /phyparams/ Qtot,Tcbc,active
      real (kind=R8) :: Qtot, Tcbc
      logical :: active
      !$OMP THREADPRIVATE(Qtot,Tcbc,active)

      ! common /location/ iwall, ns
      integer :: iwall, ns
      !$OMP THREADPRIVATE(iwall, ns)

      end module b2mod_heatdiff_common
