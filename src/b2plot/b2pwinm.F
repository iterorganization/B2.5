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
      use b2mod_dimensions
      implicit none
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
        if(streql(u,unms(i))) exit
      end do
      select case (i)
c----------------------------------------------------------------------
c*** Starting segment
      case (1)
       if(ll.eq.0) then
        call scipit
       else if(ll.lt.0 .or. ll.gt.nwldprtx) then
        write(*,'(1x,a,i4,a,i4)') 'Wrong wall segment number:',ll,
     ,                            '.  Must be positive and =< ',nwldprtx
        if(ll.gt.nwldprtx) 
     >   write (*,*) 'Check nwldprtx in b2mod_b2plot_wall_loading.F'
       else
        call resimi(lwldstrt(ll))
       end if
c----------------------------------------------------------------------
c*** Target list
      case (2)
       if(ll.eq.0) then
        call scipit
       else if(ll.le.0 .or. ll.gt.nwldprtx) then
        write(*,'(1x,a,i4,a,i4)') 'Wrong wall segment number:',ll,
     ,                            '.  Must be positive and =< ',nwldprtx
        if(ll.gt.nwldprtx) 
     >   write (*,*) 'Check nwldprtx in b2mod_b2plot_wall_loading.F'
       else
        call rearri(lwldlis(1,ll),nsgmx,nwldlis(ll))
       end if
c----------------------------------------------------------------------
c*** Plotting zone
      case (3)
       call resimi(ll)
       ll=ll+1
       nwldprt=max(nwldprt,ll)
c----------------------------------------------------------------------
c*** Shadowing structure
      case (4)
       ll=1
       nwldprt=max(nwldprt,ll)
       pltzlbl(ll)='Sh. str.'
       call scipit
c----------------------------------------------------------------------
c*** Starting points of the segments
      case (5)
       call rearre(ppp1,3*DEF_NLIM,ll)
       lppp=ll/3
       if(lppp.gt.DEF_NLIM) then !{
         write(0,*) 'Too many wall elements: maximum ',DEF_NLIM
         write(0,*) 'Current number of requested elements: ',lppp
         call xerrab
     .    (' ==> Increase DEF_NLIM in b2mod_dimensions.F')
       end if !}
       do i=1,lppp !{
        pp1(1,i)=0.001_R8*ppp1(1,i)
        pp1(2,i)=0.001_R8*ppp1(2,i)
       end do !}
c----------------------------------------------------------------------
c*** Ending points of the segments
      case (6)
       call rearre(ppp2,3*DEF_NLIM,ll)
       lppp=ll/3
       if(lppp.gt.DEF_NLIM) then !{
         write(0,*) 'Too many wall elements: maximum ',DEF_NLIM
         write(0,*) 'Current number of requested elements: ',lppp
         call xerrab
     .    (' ==> Increase DEF_NLIM in b2mod_dimensions.F')
       end if !}
       do i=1,lppp !{
        pp2(1,i)=0.001_R8*ppp2(1,i)
        pp2(2,i)=0.001_R8*ppp2(2,i)
       end do !}
c----------------------------------------------------------------------
c*** Additional source group (not implemented yet)
      case (7)
       call scipit
c----------------------------------------------------------------------
c*** Power radiated from the group of additional sources (single group)
      case (8)
       call resime(add_rdi)
       add_rad=add_rdi
c----------------------------------------------------------------------
c*** Coordinates of the additional source points
      case (9)
       call rearre(add_srci,2*nadsrcm,nadsrc)
       nadsrc=nadsrc/2
       do i=1,nadsrc
        addsrc_xy(1,i)=0.001_R8*add_srci(1,i)
        addsrc_xy(2,i)=0.001_R8*add_srci(2,i)
       end do
c----------------------------------------------------------------------
c*** 10: Inner edge of Western target
c*** 11: Inner edge of Eastern target
c*** 12: Outer edge of Western target
c*** 13: Outer edge of Eastern target
c*** 14: Values for the target edges
      case (10:14)
       if(i.le.13) then
         call resimi(i)
         luswch=i-9
       end if
       if(luswch.ge.1 .and. luswch.le.4) then
        call rearri(jedge(luswch),1,i)
       else
        call scipit
       end if
       luswch=0
c-----------------------------------------------------------------------
c*** Target edges
      case (15)
       call resimi(ll)
       if(ll.lt.0 .or. ll.gt.ntrgsx) then !{
        write(*,'(1x,a,i3,a,i3)') 'Wrong target number:',ll,
     ,                              '.  Must be positive =< ',ntrgsx
       else !}{
        ntrgs=max(ntrgs,ll)
       end if !}
c-----------------------------------------------------------------------
c*** SOL side edge
      case (16)
       if(ll.eq.0) then          !{
        call scipit
       else if(ll.lt.0 .or. ll.gt.ntrgsx) then !{
        write(*,'(1x,a,i3,a,i3)') 'Wrong target number:',ll,
     ,                              '.  Must be positive =< ',ntrgsx
       else !}{
        call resimi(ltrgedge(2,ll))
       end if !}
c-----------------------------------------------------------------------
c*** PFR side edge
      case (17)
       if(ll.eq.0) then          !{
        call scipit
       else if(ll.lt.0 .or. ll.gt.ntrgsx) then !}{
        write(*,'(1x,a,i3,a,i3)') 'Wrong target number:',ll,
     ,                              '.  Must be positive =< ',ntrgsx
       else !}{
        call resimi(ltrgedge(1,ll))
       end if !}
c-----------------------------------------------------------------------
c*** SOL side edge
      case (18)
       if(lutrgt.le.0 .or. lutrgt.gt.ntrgsx) then !{
        call scipit
       else !}{
        if(mtrgedge(2,lutrgt).gt.0) then !{
          write(*,*) 'Target',lutrgt,': repeated specification of ',
     ,                                                          unms(i)
          write (*,*) '==> corrupted param.dg?'
        end if !}
        call resimi(mtrgedge(2,lutrgt))
       end if !}
c-----------------------------------------------------------------------
c*** PFR side edge
      case (19)
       if(lutrgt.le.0 .or. lutrgt.gt.ntrgsx) then !{
        call scipit
       else !}{
        if(mtrgedge(1,lutrgt).gt.0) then !{
          write(*,*) 'Target',lutrgt,': repeated specification of ',
     ,                                                          unms(i)
          write (*,*) '==> corrupted param.dg?'
        end if !}
        call resimi(mtrgedge(1,lutrgt))
       end if !}
c----------------------------------------------------------------------
c*** Target specification
      case (20)
       call resimi(lutrgt)
       if(lutrgt.lt.0 .or. lutrgt.gt.ntrgsx) then !{
        write(*,'(1x,a,i4,a,i4)') 'Wrong target number:',lutrgt,
     ,                            '.  Must be positive and =< ', ntrgsx
        write (*,*) '==> corrupted param.dg?'
       else !}{
        ntrgsa=max(ntrgsa,lutrgt)
       end if !}
c-----------------------------------------------------------------------
c*** List of segments to be skipped
      case (21:22)
       call rearri(lskipel(nskipel+1),mskipel-nskipel,i)
       do j = nskipel+1, nskipel+i
        if(lskipel(j).gt.DEF_NLIM) then !{
         write(0,*) 'Too many wall elements: maximum ',DEF_NLIM
         write(0,*) 'Current number of requested elements: ',lskipel(j)
         call xerrab
     .    (' ==> Increase DEF_NLIM in b2mod_dimensions.F')
        end if !}
        lskip(lskipel(j)) = .true.
       end do
       nskipel=nskipel+i
c-----------------------------------------------------------------------
c*** Numbering offset
      case (23)
       call resimi(lnumboff)
c-----------------------------------------------------------------------
c*** DG v. 2.0 or later
      case (24)
       call resimi(i)
       lvdg20=.true.
c-----------------------------------------------------------------------
c*** Plot zone label
      case (25)
       if(ll.eq.0) then          !{
        call scipit
       else if(ll.le.0 .or. ll.gt.nwldprtx) then !}{
        write(*,'(1x,a,i3,a,i3)') 'Wrong wall segment number:',ll,
     ,                              '.  Must be positive =< ',nwldprtx
        if(ll.gt.nwldprtx) write (*,*) '==> check nwldprtx in b2pinpc.h'
       else !}{
        i=lena
        lena=mlblprtx
        pltzlbl(ll)=repeat(' ',8)
        call resima(pltzlbl(ll))
        lena=i
       end if !}
c-----------------------------------------------------------------------
c*** Wall material
      case (26)
       i=lena
       lena=lnwlmtrl
       call rearra(mtrwldsg,nsgmx,j)
       lena=i
c-----------------------------------------------------------------------
c*** Target material
      case (27)
       if(lutrgt.le.0 .or. lutrgt.gt.ntrgsx) then !{
        call scipit
       else !}{
        i=lena
        lena=lnwlmtrl
        call resima(mtrtrgsg(lutrgt))
        lena=i
       end if !}
c-----------------------------------------------------------------------
c*** Target sputtering - chemical
      case (28)
       if(lutrgt.le.0 .or. lutrgt.gt.ntrgsx) then !{
        call scipit
       else !}{
        call resimi(lcsptrg(lutrgt))
       end if !}
c-----------------------------------------------------------------------
c*** Target sputtering - physical
      case (29)
       if(lutrgt.le.0 .or. lutrgt.gt.ntrgsx) then !{
        call scipit
       else !}{
        call resimi(lpsptrg(lutrgt))
       end if !}
c-----------------------------------------------------------------------
c*** Wall sputtering - physical
      case (30)
       call rearri(lpspwll,mskipel,i)
c-----------------------------------------------------------------------
c*** Wall sputtering - chemical
      case (31)
       call rearri(lcspwll,mskipel,i)
c-----------------------------------------------------------------------
      case default
       ll=0
       call scipit
      end select
c======================================================================
      return
      end subroutine b2pwinm
