      subroutine b2pinp(*)
      use b2mod_types
      implicit none
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
      character*16 ufile
      external b2prdn
c----------------------------------------------------------------------
c*** Set up the input routine
c
      immed=.false.
      mimp=lun
      ufile='param.dg'
      open(lun,file=ufile,status='old',err=998)
c*** ... and read the data from dg
      call import(b2prdn)
      close(lun)
      return
c----------------------------------------------------------------------
 998  print *,'Error opening the file ',ufile
      return 1
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
      implicit none
c
c*** Reads the data from dg output file
c----------------------------------------------------------------------
#ifndef B25_EIRENE
#include "COUPLE/KOPPLDIM.F"
#endif
      integer nnns
      parameter (nnns=6)
      character*8 unm,unms(nnns),u
      data unms /
     1        'starting','partlist','wlldpart','duct    ','p1      ',
     2        'p2      '/
      integer i
      integer ll
      data ll / 0 /
c-----------------------------------------------------
c
c      print *,'b2prdn >> ',unm
      call locase(unm,u,8)
      do 1000 i=1,nnns
          if(u.eq.unms(i)) go to (10,20,30,40,50,60),i
 1000 continue
      call scipit
      return
c----------------------------------------------------------------------
 10   if(ll.le.0 .or. ll.gt.nwldprtx) then
          print '(1x,a,i4,a,i4)','Wrong wall segment number:',ll,
     ,                           '.  Must be positive and =< ',nwldprtx
          if(ll.gt.nwldprtx) 
     >     print *,'Check nwldprtx in b2mod_b2plot_wall_loading.F'
          return
      end if
      call resimi(lwldstrt(ll))
      return
c----------------------------------------------------------------------
 20   if(ll.le.0 .or. ll.gt.nwldprtx) then
          print '(1x,a,i4,a,i4)','Wrong wall segment number:',ll,
     ,                           '.  Must be positive and =< ',nwldprtx
          if(ll.gt.nwldprtx) 
     >     print *,'Check nwldprtx in b2mod_b2plot_wall_loading.F'
          return
      end if
      call rearri(lwldlis(1,ll),nlim,nwldlis(ll))
      return
c----------------------------------------------------------------------
 30   call resimi(ll)
      nwldprt=max(nwldprt,ll)
      return
c----------------------------------------------------------------------
 40   call rearri(lpumplis,nlim,npumplis)
      return
c----------------------------------------------------------------------
 50   call rearre(ppp1,3*nlim,ll)
      lppp=ll/3
      return
c----------------------------------------------------------------------
 60   call rearre(ppp2,3*nlim,ll)
      lppp=ll/3
      return
c======================================================================
      end subroutine b2prdn
