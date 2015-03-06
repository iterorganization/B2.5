      module b2mod_b2plot_ank
      use b2mod_types
      implicit none

      real (kind=R8), allocatable :: data(:)
      integer ndata

      contains

      subroutine alloc_b2mod_b2plot_ank(nwall)
      implicit none
      integer nwall

      allocate(data(nwall))

      return

      end subroutine alloc_b2mod_b2plot_ank

      end module b2mod_b2plot_ank
