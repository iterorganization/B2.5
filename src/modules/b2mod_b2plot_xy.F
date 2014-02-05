      module b2mod_b2plot_xy
      use b2mod_types

      integer, save :: LNG,nnx,nny

      real (kind=R4), allocatable, save :: 
     *  X(:),
     *  Y(:)

      contains

      subroutine alloc_b2mod_b2plot_xy(lng_in)
      integer lng_in

      LNG=lng_in
      allocate(X(1:LNG))
      allocate(Y(1:LNG))

      end subroutine alloc_b2mod_b2plot_xy

      subroutine read_b2mod_b2plot_xy(nread)
      integer nread

      end subroutine read_b2mod_b2plot_xy

      subroutine write_b2mod_b2plot_xy(nwrite)
      integer nwrite

      end subroutine write_b2mod_b2plot_xy

      end module b2mod_b2plot_xy
