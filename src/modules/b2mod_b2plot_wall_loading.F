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
      integer nwldprtx,nadsrcm,ntrgsx,mskipel,lnwlmtrl,mlblprtx
      integer nsptmx
      parameter (ntrgsx=4,nwldprtx=ntrgsx+12, nadsrcm=400)
      parameter (mskipel=DEF_NLIM,lnwlmtrl=8, mlblprtx=8)
      parameter (nsptmx=3)
      integer ildparm, irdparm, lppp, nwldprt,
     ,  lwldstrt(nwldprtx), nwldlis(nwldprtx),
     ,  jedge(2*ntrgsx), ltrgpos(ntrgsx),
     ,  lseppos(ntrgsx), lpfrpos(ntrgsx),
     ,  ntrgts, ntrgs,
     ,  nadsrc, nwldrad, npumplis, ntrgsa, lnumboff, nskipel, lytg, nytg
      integer, allocatable :: indwll(:), luwrk(:), lwldrad(:),
     ,  lwldlis(:,:), leirlis(:,:), lpumplis(:)
      integer lelmind(mskipel), lskipel(mskipel), lcspwll(mskipel),
     ,        lpspwll(mskipel), lpsptrg(ntrgsx), lcsptrg(ntrgsx),
     ,        ltrgedge(2,ntrgsx), mtrgedge(2,ntrgsx)
      logical lvdg20
      logical l2rdldw, l2rdldt, l2pwinp, lrdlsp,
     ,  lwlddir(nwldprtx), lskip(DEF_NLIM)
      real (kind=R4) :: ppp1(3,DEF_NLIM), ppp2(3,DEF_NLIM)
      real (kind=R8) :: addsrc_xy(2,nadsrcm), add_rad
      real (kind=R8), allocatable :: segm_area(:), segm_len(:),
     ,  wldrfdst(:,:), wldrdld(:), co_data(:), wldrdlsp(:,:)
      character*8 :: pltzlbl(nwldprtx)
      character mtrtrgsg(ntrgsx)*(lnwlmtrl),
     ,          lblwldpt(nwldprtx)*(mlblprtx)
      character, allocatable :: mtrwldsg(:)*(lnwlmtrl)
c======================================================================
c*** nwldprtx : max. number of the "plot zones" + 1 for shadowing
c*** nsgmx    : max. number of the wall segments + refined targets
c*** ntrgsx   : max. number of targets to be treated
c*** mlblprtx : max. length of a plot zone label
c*** nsptmx   : max. number of different sputtered materials
c*** nadsrcm  : max. number of additional radiation sources
c*** lppp     : actual number of the wall segments read from DG
c*** nwldprt  : actual number of the "plot zones" + 1 for shadowing
c*** lwldstrt : list of starting segments for the "plot zones"
c*** pltzlbl  : list of labels for the "plot zones"
c*** nwldlis  : list of lengths of the "plot zones"
c*** ntrgs    : number of targets to be treated
c*** ntrgts   : number of targets in the grid (from cuts)
c*** lwlddir  : orientation of the "plot zones" (true = CCW)
c*** lwldlis  : list of segments belonging to the "plot zones" (DG numbering)
c*** leirlis  : list of segments belonging to the "plot zones" (Eirene numbering)
c*** lnumboff : numbering offset for the wall segments in Eirene
c*** nskipel  : number of segments to be ignored by Eirene
c*** lskipel  : list of segments to be ignored by Eirene
c*** lelmind  : segment index in Eirene notation (lelmind(i) is the
c***          : Eirene segment number for i-th element in DG)
c*** lcsptrg  : chemical sputtering flag for targets
c*** lpsptrg  : physical sputtering flag for targets
c*** lcspwll  : chemical sputtering flag for wall nsegments
c*** lpspwll  : physical sputtering flag for wall nsegments
c*** lwlddir  : orientation of the "plot zones" (true = CCW)
c*** lvdg20   : file created with DG v. 2.0 or later
c*** lblwldpt : plot zone labels to appear in the output file names
c*** lnwlmtrl : maximum length of surface material name
c*** mtrwldsg : material of the wall segment
c*** mtrtrgsg : material of the target
c*** ildparm  : parameter passed to load calculation routines
c*** lskip    : list of DG segments not for Eirene
c*** nwldrad  : number of segments where rad. load is to be calculated
c*** lwldrad  : list of segments where rad. load is to be calculated
c*** ltrgedge : list of segments to be linked to the target edges
c***            (1,j) means PFR edge of target j
c***            if not specified, jedg[io][12] are used
c*** mtrgedge : list of segments to be linked to the target edges
c***            if not specified, ltrgedge or jedg[io][12] are used
c*** jedgi1   : inner edge of inner target => jedge(1)
c*** jedgi2   : inner edge of outer target => jedge(2)
c*** jedgo1   : outer edge of inner target => jedge(3) (SN) or jedge(5) (DN)
c*** jedgo2   : outer edge of outer target => jedge(4) (SN) or jedge(6) (DN)
c*** jedgi3   : inner edge of upper inner target => jedge(3)
c*** jedgi4   : inner edge of upper outer target => jedge(4)
c*** jedgo3   : outer edge of upper inner target => jedge(7)
c*** jedgo4   : outer edge of upper outer target => jedge(8)
c*** nytg     : number of radial cells in the target ( < ny-1 for limiter)
c*** lytg     : starting radial cell for targets ( > 0 for limiter)
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
c*** wldrdlsp : B2 radiation load split by species
c*** l2pwinp  : geometry input block was called
c*** l2rdldw  : radiation load was calculated for the walls
c*** l2rdldt  : radiation load was calculated for the targets
c*** lrdlsp   : radiation load split was calculated
c*** luwrk    : working array
c*** co_data  : working array (originally, for average energy)
c======================================================================

      contains

      subroutine alloc_b2plot_wall_loading(nlim,nsgmx)
      implicit none
      integer nlim, nsgmx
      
      allocate(indwll(nsgmx), luwrk(nsgmx), lwldrad(nsgmx))
      allocate(lwldlis(nsgmx,nwldprtx))
      allocate(leirlis(nsgmx,nwldprtx))
      allocate(lpumplis(nlim))
      allocate(wldrdlsp(nsgmx,DEF_NATM))
      allocate(segm_area(nsgmx), segm_len(nsgmx),
     ,  wldrfdst(nsgmx,nwldprtx), wldrdld(nsgmx), co_data(nsgmx))
      allocate(mtrwldsg(nsgmx))
      lskip = .false.

      return
      end subroutine alloc_b2plot_wall_loading

      subroutine dealloc_b2plot_wall_loading
      implicit none
      
      deallocate(indwll, luwrk, lwldrad)
      deallocate(lwldlis)
      deallocate(leirlis)
      deallocate(lpumplis)
      deallocate(wldrdlsp)
      deallocate(segm_area, segm_len, wldrfdst, wldrdld, co_data)
      deallocate(mtrwldsg)

      return
      end subroutine dealloc_b2plot_wall_loading

      end module b2mod_b2plot_wall_loading
