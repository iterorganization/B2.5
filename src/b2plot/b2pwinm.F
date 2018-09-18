      subroutine b2pwinm(unm)
      use b2mod_types
      use b2mod_b2plot
      use b2mod_b2plot_wall_loading
c
c  version : 11.11.97 21:36
c
c======================================================================
c*** Reads the data from dg output file (called from imprt)
c*** MUST be compiled in single precision!!!
c======================================================================
#ifdef B25_EIRENE
      use eirmod_parmmod
      use eirmod_extrab25
#endif
      implicit none
#include <DIMENSIONS.F>
#ifndef B25_EIRENE
#include "COUPLE/KOPPLDIM.F"
#endif
#include "impcon.inc"
      integer nnns
      real (kind=R4) :: add_rdi, add_srci(2, nadsrcm)
      parameter (nnns=31)
      character*8 unm,unms(nnns),u
      integer ll, lutrgt, luswch
      integer i, j
      logical streql
      external streql
      data unms /
     1        'starting','segmlist','plotzone','shdw_str','p1      ',
     2        'p2      ','rad_srcs','radpower','sources ','jedgi1  ',
     3        'jedgi2  ','jedgo1  ','jedgo2  ','elemlist','edge_trg',
     4        'edge_sol','edge_pfr','sol_edge','pfr_edge','trg_spcf',
     5        'skipelem','misselem','laddsrfs','gleirdat','pltz_lbl',
     6        'material','made_of' ,'tgspcsw' ,'tgsppsw' ,'jlsptp'  ,
     7        'jlsptc'  /
      data ll / 0 /, luswch / 0 /, lutrgt / 0 /
c======================================================================
c
c      print *,'b2pwinm >> ',unm
      call locase(unm,u,8)
      do i=1,nnns
        if(streql(u,unms(i))) go to
     ,     ( 10, 20, 30, 40, 50, 60, 70, 80, 90,100,
     ,      110,120,130,140,150,160,170,180,190,200,
     ,      210,210,230,240,250,260,270,280,290,300,
     ,      310),i
      end do
      ll=0
      call scipit
      return
c----------------------------------------------------------------------
c*** Starting segment
   10 if(ll.eq.0) then
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
   20 if(ll.eq.0) then
        call scipit
      else if(ll.le.0 .or. ll.gt.nwldprtx) then
        write(*,'(1x,a,i4,a,i4)') 'Wrong wall segment number:',ll,
     ,                            '.  Must be positive and =< ',nwldprtx
        if(ll.gt.nwldprtx) 
     >   write (*,*) 'Check nwldprtx in b2mod_b2plot_wall_loading.F'
      else
        call rearri(lwldlis(1,ll),nsgmx,nwldlis(ll))
      end if
      return
c----------------------------------------------------------------------
c*** Plotting zone
   30 call resimi(ll)
      ll=ll+1
      nwldprt=max(nwldprt,ll)
      return
c----------------------------------------------------------------------
c*** Shadowing structure
   40 ll=1
      nwldprt=max(nwldprt,ll)
      call scipit
      return
c----------------------------------------------------------------------
c*** Starting points of the segments
   50 call rearre(ppp1,3*DEF_NLIM,ll)
      lppp=ll/3
      do i=1,lppp !{
        pp1(1,i)=0.001_R8*ppp1(1,i)
        pp1(2,i)=0.001_R8*ppp1(2,i)
      end do !}
      return
c----------------------------------------------------------------------
c*** Ending points of the segments
   60 call rearre(ppp2,3*DEF_NLIM,ll)
      lppp=ll/3
      do i=1,lppp !{
        pp2(1,i)=0.001_R8*ppp2(1,i)
        pp2(2,i)=0.001_R8*ppp2(2,i)
      end do !}
      return
c----------------------------------------------------------------------
c*** Additional source group (not implemented yet)
   70 call scipit
      return
c----------------------------------------------------------------------
c*** Power radiated from the group of additional sources (single group)
   80 call resime(add_rdi)
      add_rad=add_rdi
      return
c----------------------------------------------------------------------
c*** Coordinates of the additional source points
   90 call rearre(add_srci,2*nadsrcm,nadsrc)
      nadsrc=nadsrc/2
      do i=1,nadsrc
        addsrc_xy(1,i)=0.001_R8*add_srci(1,i)
        addsrc_xy(2,i)=0.001_R8*add_srci(2,i)
      end do
      return
c----------------------------------------------------------------------
c*** Inner edge of inner target
  100 call resimi(i)
      luswch=1
      goto 140
c----------------------------------------------------------------------
c*** Inner edge of outer target
  110 call resimi(i)
      luswch=2
      goto 140
c----------------------------------------------------------------------
c*** Outer edge of inner target
  120 call resimi(i)
      luswch=3
      goto 140
c----------------------------------------------------------------------
c*** Outer edge of outer target
  130 call resimi(i)
      luswch=4
      goto 140
c----------------------------------------------------------------------
c*** Values for the target edges
  140 if(luswch.ge.1 .and. luswch.le.4) then
        call rearri(jedge(luswch),1,i)
      else
        call scipit
      end if
      luswch=0
      return
c-----------------------------------------------------------------------
c*** Target edges
  150 call resimi(ll)
      if(ll.lt.0 .or. ll.gt.ntrgsx) then !{
        write(*,'(1x,a,i3,a,i3)') 'Wrong target number:',ll,
     ,                              '.  Must be positive =< ',ntrgsx
      else !}{
        ntrgs=max(ntrgs,ll)
      end if !}
      return
c-----------------------------------------------------------------------
c*** SOL side edge
  160 if(ll.eq.0) then          !{
        call scipit
      else if(ll.lt.0 .or. ll.gt.ntrgsx) then !{
        write(*,'(1x,a,i3,a,i3)') 'Wrong target number:',ll,
     ,                              '.  Must be positive =< ',ntrgsx
      else !}{
        call resimi(ltrgedge(2,ll))
      end if !}
      return
c-----------------------------------------------------------------------
c*** PFR side edge
  170 if(ll.eq.0) then          !{
        call scipit
      else if(ll.lt.0 .or. ll.gt.ntrgsx) then !}{
        write(*,'(1x,a,i3,a,i3)') 'Wrong target number:',ll,
     ,                              '.  Must be positive =< ',ntrgsx
      else !}{
        call resimi(ltrgedge(1,ll))
      end if !}
      return
c-----------------------------------------------------------------------
c*** SOL side edge
  180 if(lutrgt.le.0 .or. lutrgt.gt.ntrgsx) then !{
        call scipit
      else !}{
        if(mtrgedge(2,lutrgt).gt.0) then !{
          write(*,*) 'Target',lutrgt,': repeated specification of ',
     ,                                                          unms(i)
          write (*,*) '==> corrupted param.dg?'
        end if !}
        call resimi(mtrgedge(2,lutrgt))
      end if !}
      return
c-----------------------------------------------------------------------
c*** PFR side edge
  190 if(lutrgt.le.0 .or. lutrgt.gt.ntrgsx) then !{
        call scipit
      else !}{
        if(mtrgedge(1,lutrgt).gt.0) then !{
          write(*,*) 'Target',lutrgt,': repeated specification of ',
     ,                                                          unms(i)
          write (*,*) '==> corrupted param.dg?'
        end if !}
        call resimi(mtrgedge(1,lutrgt))
      end if !}
      return
c----------------------------------------------------------------------
c*** Target specification
  200 call resimi(lutrgt)
      if(lutrgt.lt.0 .or. lutrgt.gt.ntrgsx) then !{
        write(*,'(1x,a,i4,a,i4)') 'Wrong target number:',lutrgt,
     ,                            '.  Must be positive and =< ', ntrgsx
        write (*,*) '==> corrupted param.dg?'
      else !}{
        ntrgsa=max(ntrgsa,lutrgt)
      end if !}
      return
c-----------------------------------------------------------------------
c*** List of segments to be skipped
  210 call rearri(lskipel(nskipel+1),mskipel-nskipel,i)
      do j = nskipel+1, nskipel+i
        lskip(lskipel(j)) = .true.
      end do
      nskipel=nskipel+i
      return
c-----------------------------------------------------------------------
c*** Numbering offset
  230 call resimi(lnumboff)
      return
c-----------------------------------------------------------------------
c*** DG v. 2.0 or later
  240 call resimi(i)
      lvdg20=.true.
      return
c-----------------------------------------------------------------------
c*** Plot zone label
  250 if(ll.eq.0) then          !{
        call scipit
      else if(ll.le.0 .or. ll.gt.nwldprtx) then !}{
        write(*,'(1x,a,i3,a,i3)') 'Wrong wall segment number:',ll,
     ,                              '.  Must be positive =< ',nwldprtx
        if(ll.gt.nwldprtx) write (*,*) '==> check nwldprtx in b2pinpc.h'
      else !}{
        i=lena
        lena=mlblprtx
        lblwldpt(ll)=' '
        call resima(lblwldpt(ll))
        lena=i
      end if !}
      return
c-----------------------------------------------------------------------
c*** Wall material
  260 i=lena
      lena=lnwlmtrl
      call rearra(mtrwldsg,nsgmx,j)
      lena=i
      return
c-----------------------------------------------------------------------
c*** Target material
  270 if(lutrgt.le.0 .or. lutrgt.gt.ntrgsx) then !{
        call scipit
      else !}{
        i=lena
        lena=lnwlmtrl
        call resima(mtrtrgsg(lutrgt))
        lena=i
      end if !}
      return
c-----------------------------------------------------------------------
c*** Target sputtering - chemical
  280 if(lutrgt.le.0 .or. lutrgt.gt.ntrgsx) then !{
        call scipit
      else !}{
        call resimi(lcsptrg(lutrgt))
      end if !}
      return
c-----------------------------------------------------------------------
c*** Target sputtering - physical
  290 if(lutrgt.le.0 .or. lutrgt.gt.ntrgsx) then !{
        call scipit
      else !}{
        call resimi(lpsptrg(lutrgt))
      end if !}
      return
c-----------------------------------------------------------------------
c*** Wall sputtering - physical
  300 call rearri(lpspwll,mskipel,i)
      return
c-----------------------------------------------------------------------
c*** Wall sputtering - chemical
  310 call rearri(lcspwll,mskipel,i)
      return
c======================================================================
      end subroutine b2pwinm
