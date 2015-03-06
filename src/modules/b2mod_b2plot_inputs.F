      module b2mod_b2plot_inputs
      use b2mod_types
      implicit none

c  VERSION : 20.03.97 13:08
      integer nwldprtx, nwldprt, npumplis, lppp
      parameter (nwldprtx=12)
      integer lwldstrt(nwldprtx), nwldlis(nwldprtx)
      integer, allocatable :: lwldlis(:,:), lpumplis(:)
      real (kind=R4), allocatable :: ppp1(:,:), ppp2(:,:)

      contains

      subroutine alloc_b2plot_inputs(nlim)
      implicit none
      integer nlim

      allocate(lwldlis(nlim,nwldprtx))
      allocate(lpumplis(nlim))
      allocate(ppp1(3,nlim),ppp2(3,nlim))

      return
      end subroutine alloc_b2plot_inputs

      end module b2mod_b2plot_inputs
