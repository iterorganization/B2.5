      subroutine b2pinp(iret)
      use b2mod_types
      use b2mod_dimensions
      implicit none
      integer iret
c
c  version : 11.11.97 21:29
c
c----------------------------------------------------------------------
c*** Data input from dg files.
c*** MUST be compiled in single precision!!!
c----------------------------------------------------------------------
#include "impcon.inc"
      integer lun
      parameter (lun=96)
      logical open_file, file_ok
      character*24 ufile
      external b2prdn, open_file
c----------------------------------------------------------------------
c*** Set up the input routine
c
      iret=0
      lena=8
      jmess=2
      lfin=72
      keylen=8
      immed=.false.
      mimp=lun
      fixkey=.false.
      finish = 'finish'
      promp(1) = repeat(' ',8)
      promp(2) = repeat(' ',8)
      ufile='param.dg'
      call find_file(ufile,file_ok)
c      open(lun,file=ufile,status='old',err=998)
      if(open_file(lun,ufile,'modules/Uinp')) then
c*** ... and read the data from dg
        call import(b2prdn)
        close(lun)
      else
        goto 998
      endif
      return
c----------------------------------------------------------------------
  998 print *,'Error opening the file ',ufile
      iret=1
      return
c----------------------------------------------------------------------
      end subroutine b2pinp
c
c
c======================================================================
      subroutine b2prdn(unm)
      use b2mod_types
      use b2mod_b2plot_wall_loading
#ifdef B25_EIRENE
      use eirmod_parmmod
      use eirmod_extrab25
#endif
      use b2mod_dimensions
      implicit none
c
c*** Reads the data from dg output file
c----------------------------------------------------------------------
#ifndef B25_EIRENE
#include "COUPLE/KOPPLDIM.F"
#endif
#include "impcon.inc"
      integer nnns
      real (kind=R4) :: add_rdi, add_srci(2, nadsrcm)
      parameter (nnns=25)
      character*8 unm,unms(nnns),u
      data unms /
     1        'starting','partlist','wlldpart','duct    ','p1      ',
     2        'p2      ','segmlist','plotzone','jedgi1  ','jedgi2  ',
     3        'jedgo1  ','jedgo2  ','elemlist','sol_edge','pfr_edge',
     4        'trg_spcf','pltz_lbl','misselem','skipelem','shdw_str',
     5        'radpower','sources ','gleirdat','fclbl'   ,'grd_mode'/
      integer i, j
      integer ll, lt, luswch
      logical streql
      external streql
      data ll / 0 /, luswch / 0 /, lt / 0 /
      save
c======================================================================
c
c      print *,'b2prdn >> ',unm
      call locase(unm,u,8)
      i = 0
      do i = 1, nnns
        if (streql(u,unms(i))) exit
      end do
      select case (i)
      case (1)
        if(ll.le.0 .or. ll.gt.nwldprtx) then
          print '(1x,a,i4,a,i4)','Wrong wall segment number:',ll,
     ,                           '.  Must be positive and =< ',nwldprtx
          if(ll.gt.nwldprtx) 
     >     print *,'Check nwldprtx in b2mod_b2plot_wall_loading.F'
        else
          call resimi(lwldstrt(ll))
        end if
      case (2,7)
        if(ll.le.0 .or. ll.gt.nwldprtx) then
          print '(1x,a,i4,a,i4)','Wrong wall segment number:',ll,
     ,                           '.  Must be positive and =< ',nwldprtx
          if(ll.gt.nwldprtx) 
     >     print *,'Check nwldprtx in b2mod_b2plot_wall_loading.F'
        else
          call rearri(lwldlis(1,ll),nsgmx,nwldlis(ll))
        end if
      case (3,8)
        call resimi(ll)
        ll=ll+1
        nwldprt=max(nwldprt,ll)
      case (4)
        call rearri(lpumplis,nlim,npumplis)
      case (5)
        call rearre(ppp1,3*DEF_NLIM,ll)
        lppp=ll/3
        if(lppp.gt.DEF_NLIM) then !{
          write(0,*) 'Too many wall elements: maximum ',DEF_NLIM
          write(0,*) 'Current number of requested elements: ',lppp
          call xerrab
     .     (' ==> Increase DEF_NLIM in b2mod_dimensions.F')
        end if !}
      case (6)
        call rearre(ppp2,3*DEF_NLIM,ll)
        lppp=ll/3
        if(lppp.gt.DEF_NLIM) then !{
          write(0,*) 'Too many wall elements: maximum ',DEF_NLIM
          write(0,*) 'Current number of requested elements: ',lppp
          call xerrab
     .     (' ==> Increase DEF_NLIM in b2mod_dimensions.F')
        end if !}
      case (9:13)
        if (streql(u,unms(9))) then
c*** Inner edge of Western target (inner for LSN, outer for USN)
          call resimi(i)
          luswch=1
        else if (streql(u,unms(10))) then
c*** Inner edge of Eastern target (outer for LSN, inner for USN)
          call resimi(i)
          luswch=2
        else if (streql(u,unms(11))) then
c*** Outer edge of Western target (inner for LSN, outer for USN)
          call resimi(i)
          luswch=3
        else if (streql(u,unms(12))) then
c*** Outer edge of Eastern target (outer for LSN, inner for USN)
          call resimi(i)
          luswch=4
        end if
c*** Values for the target edges (old style, only worked for SN topologies)
        if(luswch.ge.1 .and. luswch.le.4) then
          call rearri(jedge(luswch),1,i)
        else
          call scipit
        end if
        luswch=0
      case (14)
c*** New style values for the target SOL edges
        if (lt.eq.0) then
          call scipit
        else if(lt.lt.0 .or. lt.gt.ntrgsx) then
          write(*,'(1x,a,i4,a,i4)') 'Wrong target number:',lt,
     ,                            '.  Must be positive and =< ', ntrgsx
        else if (lt.ge.1 .and. lt.le.ntrgsx) then
          call resimi(jedge(ntrgsx+lt))
        end if
      case (15)
c*** New style values for the target PFR edges
        if (lt.eq.0) then
          call scipit
        else if(lt.lt.0 .or. lt.gt.ntrgsx) then
          write(*,'(1x,a,i4,a,i4)') 'Wrong target number:',lt,
     ,                            '.  Must be positive and =< ', ntrgsx
        else if (lt.ge.1 .and. lt.le.ntrgsx) then
          call resimi(jedge(lt))
        end if
      case (16)
c*** Target specifications
        call resimi(lt)
        if(lt.lt.0 .or. lt.gt.ntrgsx)
     >   write(*,'(1x,a,i4,a,i4)') 'Wrong target number:',lt,
     ,                            '.  Must be positive and =< ', ntrgsx
      case (17)
c*** Plotzone label
        i=lena
        lena=8
        call resima(pltzlbl(ll))
        lena=i
      case (18:19)
c*** Elements not for Eirene
        call rearri(lskipel(nskipel+1),mskipel-nskipel,i)
        do j = nskipel+1, nskipel+i
          if(lskipel(j).gt.DEF_NLIM) then !{
            write(0,*) 'Too many wall elements: maximum ',DEF_NLIM
            write(0,*) 'Current number of requested elements: ',
     .        lskipel(j)
            call xerrab
     .       (' ==> Increase DEF_NLIM in b2mod_dimensions.F')
          end if !}
          lskip(lskipel(j)) = .true.
        end do
        nskipel=nskipel+i
      case (20)
c*** Shadowing structure
        ll=1
        nwldprt=max(nwldprt,ll)
        pltzlbl(ll)='Sh. Str.'
        call scipit
      case (21)
c*** Power radiated from the group of additional sources (single group)
        call resime(add_rdi)
        add_rad=add_rdi
      case (22)
c*** Coordinates of the additional source points
        call rearre(add_srci,2*nadsrcm,nadsrc)
        nadsrc=nadsrc/2
        do i=1,nadsrc
          addsrc_xy(1,i)=0.001_R8*add_srci(1,i)
          addsrc_xy(2,i)=0.001_R8*add_srci(2,i)
        end do
c*** DG v. 2.0 or later
      case (23)
        call resimi(i)
        lvdg20=.true.
c*** Face labels
      case (24)
        call rearri(lfclbl,DEF_NLIM,ll)
        if(ll.gt.DEF_NLIM) then !{
          write(0,*) 'Too many wall elements: maximum ',DEF_NLIM
          write(0,*) 'Current number of requested elements: ',ll
          call xerrab
     .     (' ==> Increase DEF_NLIM in b2mod_dimensions.F')
        end if !}
c*** Grid mode
      case (25)
        call resimi(grid_mode)
      case default
        call scipit
      end select
      return
c======================================================================
      end subroutine b2prdn
