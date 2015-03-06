      subroutine b2pwinm(unm)
      use b2mod_types
      use b2mod_b2plot_wall_loading
c
c  version : 11.11.97 21:36
c
c======================================================================
c*** Reads the data from dg output file (called from imprt)
c*** MUST be compiled in single precision!!!
c======================================================================
#ifdef USE_EIRENE
      use eirmod_parmmod
      use eirmod_extrab25
#endif
      implicit none
#ifndef USE_EIRENE
#include "COUPLE/KOPPLDIM.F"
#endif
      integer nnns
      parameter (nnns=14)
      character*8 unm,unms(nnns),u
      integer ll, luswch
      integer i
      data unms /
     1        'starting','segmlist','plotzone','shdw_str','p1      ',
     2        'p2      ','rad_srcs','radpower','sources ','jedgi1  ',
     3        'jedgi2  ','jedgo1  ','jedgo2  ','elemlist'/
      data ll / 0 /, luswch / 0 /
c======================================================================
c
c      print *,'b2pwinm >> ',unm
      call locase(unm,u,8)
      do i=1,nnns
        if(u.eq.unms(i)) go to (10,20,30,40,50,60,70,80,90,100,
     ,                                               110,120,130,140),i
      end do
      ll=0
      call scipit
      return
c----------------------------------------------------------------------
c*** Starting segment
 10   if(ll.eq.0) then
        call scipit
      else if(ll.lt.0 .or. ll.gt.nwldprtx) then
        write(*,'(1x,a,i4,a,i4)') 'Wrong wall segment number:',ll,
     ,                            '.  Must be positive and =< ',nwldprtx
        if(ll.gt.nwldprtx) 
     >   write (*,*) 'Check nwldprtx in b2mod_b2plot_wall_loading.F'
      else
        call resimi(lwldstrt(ll))
      end if
      return
c----------------------------------------------------------------------
c*** Target list
 20   if(ll.eq.0) then
        call scipit
      else if(ll.le.0 .or. ll.gt.nwldprtx) then
        write(*,'(1x,a,i4,a,i4)') 'Wrong wall segment number:',ll,
     ,                            '.  Must be positive and =< ',nwldprtx
        if(ll.gt.nwldprtx) 
     >   write (*,*) 'Check nwldprtx in b2mod_b2plot_wall_loading.F'
      else
        call rearri(lwldlis(1,ll),nlim,nwldlis(ll))
      end if
      return
c----------------------------------------------------------------------
c*** Plotting zone
 30   call resimi(ll)
      ll=ll+1
      nwldprt=max(nwldprt,ll)
      return
c----------------------------------------------------------------------
c*** Shadowing structure
 40   ll=1
      nwldprt=max(nwldprt,ll)
      call scipit
      return
c----------------------------------------------------------------------
c*** Starting points of the segments
 50   call rearre(ppp1,3*nlim,ll)
      lppp=ll/3
      return
c----------------------------------------------------------------------
c*** Ending points of the segments
 60   call rearre(ppp2,3*nlim,ll)
      lppp=ll/3
      return
c----------------------------------------------------------------------
c*** Additional source group (not implemented yet)
 70   call scipit
      return
c----------------------------------------------------------------------
c*** Power radiated from the group of additional sources (single group)
 80   call resime(add_rad)
      return
c----------------------------------------------------------------------
c*** Coordinates of the additional source points
 90   call rearre(addsrc_xy,2*nadsrcm,nadsrc)
      nadsrc=nadsrc/2
      do i=1,nadsrc
        addsrc_xy(1,i)=0.001*addsrc_xy(1,i)
        addsrc_xy(2,i)=0.001*addsrc_xy(2,i)
      end do
      return
c----------------------------------------------------------------------
c*** Inner edge of inner target
 100  call resimi(i)
      luswch=1
      return
c----------------------------------------------------------------------
c*** Inner edge of outer target
 110  call resimi(i)
      luswch=2
      return
c----------------------------------------------------------------------
c*** Outer edge of inner target
 120  call resimi(i)
      luswch=3
      return
c----------------------------------------------------------------------
c*** Outer edge of outer target
 130  call resimi(i)
      luswch=4
      return
c----------------------------------------------------------------------
c*** Values for the target edges
 140  if(luswch.ge.1 .and. luswch.le.4) then
        call rearri(jedge(luswch),1,i)
      else
        call scipit
      end if
      luswch=0
      return
c======================================================================
      end subroutine b2pwinm
