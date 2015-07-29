      module b2mod_b2plot_wall_loading
      use b2mod_types
      implicit none
#include <DIMENSIONS.F>
c======================================================================
c*** Header file for the routines calculating the wall loading
c*** !!! If used, must be preceded by DIMENSIONS.F !!!
c======================================================================
c
c  version : 06.09.97 19:59
c
      integer nwldprtx,nadsrcm
      parameter (nwldprtx=4+12, nadsrcm=400)
      integer lppp, nwldprt, lwldstrt(nwldprtx), nwldlis(nwldprtx),
     ,	jedge(4),
     ,	nadsrc, nwldrad, npumplis 
      integer, allocatable :: indwll(:), luwrk(:), lwldrad(:),
     ,  lwldlis(:,:), lpumplis(:)
      logical l2rdldw, l2rdldt, l2pwinp, lwlddir(nwldprtx)
      real (kind=R4) :: ppp1(3,DEF_NLIM), ppp2(3,DEF_NLIM)
      real (kind=R4) :: addsrc_xy(2,nadsrcm), add_rad
      real (kind=R8), allocatable :: segm_area(:), segm_len(:),
     ,	wldrfdst(:,:), wldrdld(:), co_data(:)
c======================================================================
c*** nwldprtx : max. number of the "plot zones" + 1 for shadowing
c*** nsgmx    : max. number of the wall segments + refined targets
c*** nadsrcm  : max. number of additional radiation sources
c*** lppp     : actual number of the wall segments read from DG
c*** nwldprt  : actual number of the "plot zones" + 1 for shadowing
c*** lwldstrt : list of starting segments for the "plot zones"
c*** nwldlis  : list of lengths of the "plot zones"
c*** lwlddir  : orientation of the "plot zones" (true = CCW)
c*** lwldlis  : list of segments belonging to the "plot zones"
c*** nwldrad  : number of segments where rad. load is to be calculated
c*** lwldrad  : list of segments where rad. load is to be calculated
c*** jedgi1   : inner edge of inner target => jedge(1)
c*** jedgi2   : inner edge of outer target => jedge(2)
c*** jedgo1   : outer edge of inner target => jedge(3)
c*** jedgo2   : outer edge of outer target => jedge(4)
c*** nadsrc   : actual number of additional radiation sources
c*** indwll   : index of the wall segments to be plotted
c*** ppp1     : coordinates of the starting points of the segments
c*** ppp2     : coordinates of the ending points of the segments
c*** addsrc_xy: coordinates of the additional radiation sources
c*** add_rad  : total power radiated from the additional sources
c*** segm_area: areas of the wall segments
c*** segm_len : lengths of the wall segments
c*** wldrfdst : reference distances used as x-ordinate for 1-D plots
c*** wldrdld  : radiation load calculated in b2pwrld
c*** l2pwinp  : geometry input block was called
c*** l2rdldw  : radiation load was calculated for the walls
c*** l2rdldt  : radiation load was calculated for the targets
c*** luwrk    : working array
c*** co_data  : working array (originally, for average energy)
c======================================================================

      contains

      subroutine alloc_b2plot_wall_loading(nlim,nsgmx)
      implicit none
      integer nlim, nsgmx
      
      allocate(indwll(nsgmx), luwrk(nsgmx), lwldrad(nsgmx))
      allocate(lwldlis(nsgmx,nwldprtx))
      allocate(lpumplis(nlim))
      allocate(segm_area(nsgmx), segm_len(nsgmx),
     ,	wldrfdst(nsgmx,nwldprtx), wldrdld(nsgmx), co_data(nsgmx))

      return
      end subroutine alloc_b2plot_wall_loading

      end module b2mod_b2plot_wall_loading
