      subroutine b2pnwl(*)
c
c  version : 11.11.97 20:52
c
c======================================================================
c*** b2pnwl: neutral wall load (temporary solution)
c***
c======================================================================
      use b2mod_types
      use b2mod_b2plot
      use b2mod_constants
csw 02feb2011
      use eirmod_extrab25
csw
      implicit none
csw 07feb2011
#     include "COUPLE/KOPPLDIM.F"
csw#     include "COUPLE/EIRDIAG.F"
#     include "b2pinpc.h"
c*** Local variables
c
      integer lun
      parameter (lun=96)
c
c*** flxctff (flux cutoff) represents zero for the particle fluxes
c*** to the wall in order to allow one to plot them on logarithmic
c*** scale. enrctff and pwrctff do the same for the particle energy
c*** and for the power flux.
      real (kind=R8) ::
     *   flxctff,enrctff,pwrctff
      parameter (flxctff=0.1e0_R8, enrctff=1.e-8_R8, pwrctff=1.e-8_R8)
      integer il,is,i,ii,ia,k,ic
      real (kind=R8) ::
     *   x,dist,h,a,aa,hlpx,hlpw,wwt,wwk,wwp,wwyn,wwyg,
     ,     wwj(natm+nmol),wwe(natm+nmol)
      character*8 text_an(natm),text_mn(natm)
      intrinsic max
c----------------------------------------------------------------------
      print *,'Starting b2pnwl...',nlimi,nnlimi
c*** Set the default values
      npumplis=0
      nwldprt=0
c*** Read the data from dg files
      call b2pinp(*999)
c
c*** Convert the input wall co-ordinates into double precision
c*** and scale them to meters.
c*** We have now two sets of the wall segments, one saved by Eirene
c*** and (another) one read from dg file. Only the segments which are
c*** identical in both sets are to be used for plotting, the rest
c*** will be filled with the value representing zero in the list-
c*** sorting routine b2pwsrt. The reference distance along the wall
c*** will however be calculated from the geometry read from dg file.
c*** This should allow one to have the same reference frame for the
c*** radiation and neutral wall load.
c
      do i=1,lppp
        pp1(1,i)=0.001e0_R8*ppp1(1,i)
        pp1(2,i)=0.001e0_R8*ppp1(2,i)
        pp2(1,i)=0.001e0_R8*ppp2(1,i)
        pp2(2,i)=0.001e0_R8*ppp2(2,i)
      end do
cc%%%
c      print *,'After b2pinp',nwldprt
c      print '(a,(20i5))','lwldstrt:',(lwldstrt(i),i=1,nwldprt)
c      print '(a,(20i5))','nwldlis :',(nwldlis(i),i=1,nwldprt)
c      do il=1,nwldprt
c        print '(a,(20i5))','lwldlis :',il,(lwldlis(i,il),i=1,nwldprt)
c      end do
cc%%%
c*** Then sort the target lists
      call b2pwsrt(*999)
cc<<<
c      print *,'After b2pwsrt',nwldprt
c      print '(a,(20i5))','lwldstrt:',(lwldstrt(i),i=1,nwldprt)
c      print '(a,(20i5))','nwldlis :',(nwldlis(i),i=1,nwldprt)
c      do il=1,nwldprt
c        print '(a,(20i5))','lwldlis :',il,
c     ,                                  (lwldlis(i,il),i=1,nwldlis(il))
c      end do
cc>>>
c*** ... and produce the output file
      do i=1,natmi
        text_an(i)=textan(i)(2:8)
      end do
      do i=1,nmoli
        text_mn(i)=textmn(i)(2:8)
      end do
      open (lun,file='wlld.dat')
      write(lun,
     ,      '(2x,a4,5x,a3,6(5x,a2,3x),4(5x,a4,3x),16(3x,a4,a4,1x))')
     -     'zone','elm',
     ,     'R1','Z1','R2','Z2',' X','Xr',
     ,     'area','Wtot','Wkin','Wpot',
     ,     ('flx_',text_an(i),'Eav_',text_an(i),i=1,natmi),
     ,     ('flx_',text_mn(i),'Eav_',text_mn(i),i=1,nmoli),
     ,      'net_','ersn','grs_','ersn'
c<<<
c      print *,'nwldprt,natmi,nmoli',nwldprt,natmi,nmoli
c      do il=1,nwldprt
c        print '(1x,20i4)',il,nwldlis(il)
c        print '(1x,20i4)',(lwldlis(i,il),i=1,nwldlis(il))
c      end do
c>>>
c
c*** Check whether carbon atoms are present in the calculation
c
      ic=0
      do i=1,natmi
        if(text_an(i) .eq. 'C') ic=i
      end do
c
      do il=1,nwldprt
        if(nwldlis(il) .gt. 0) then
          write(lun,*)
c*** calculate the length of the target segment
          hlpx=0.e0_R8
          do i=1,nwldlis(il)
            ia=lwldlis(i,il)
            is=iabs(ia)
            dist=sqrt((p1(1,is)-p2(1,is))**2+(p1(2,is)-p2(2,is))**2)
            hlpx=hlpx+dist
          end do
c
          x=0.e0_R8
          dist=0.e0_R8
          hlpw=0.e0_R8
          do i=1,nwldlis(il)
            ia=lwldlis(i,il)
            is=iabs(ia)
            h=dist
            dist=sqrt((pp1(1,is)-pp2(1,is))**2+(pp1(2,is)-pp2(2,is))**2)
            a=pi*(pp1(1,is)+pp2(1,is))*dist
            dist=0.5e0_R8*dist
            x=x+h+dist
            if(ia .gt. 0) then
              aa=1.e-6_R8/a
              hlpw=hlpw+wldnek(is,0)+wldnep(is,0)
              wwk=wldnek(is,0)*aa
              wwp=wldnep(is,0)*aa
              wwt=max(wwk+wwp,pwrctff)
c              wwk=max(wwk,pwrctff)
c              wwp=max(wwp,pwrctff)
              k=0
              do ii=1,natmi
                k=k+1
                wwj(k)=max(wldna(is,ii,0),flxctff)/a
                wwe(k)=max(ewlda(is,ii,0),enrctff)
              end do
              do ii=1,nmoli
                k=k+1
                wwj(k)=max(wldnm(is,ii,0),flxctff)/a
                wwe(k)=max(ewldm(is,ii,0),enrctff)
              end do
c*** Calculate the net and gross erosion
c*** (temporary solution, no molecules involved)
              if(ic .gt. 0) then
                wwyg=wldspt(is,0)/a
                wwyn=(wldra(is,ic,0)-wldna(is,ic,0))/a
              else
                wwyn=flxctff
                wwyg=flxctff
              end if
            else
              wwyn=flxctff
              wwyg=flxctff
              wwk=pwrctff
              wwp=pwrctff
              wwt=pwrctff
              k=natmi+nmoli
              do ii=1,k
                wwj(ii)=0.1e0_R8*flxctff
                wwe(ii)=enrctff
              end do
            end if
            write (lun,'(1x,a5,i2.2,i5,5f10.4,1p,24e12.4)') 'zone#',
     ,                il,is,pp1(1,is),pp1(2,is),pp2(1,is),pp2(2,is),
     ,                x,x/hlpx,a,wwt,wwk,wwp,(wwj(ii),wwe(ii),ii=1,k),
     ,                wwyn,wwyg
          end do
          write(*,*) 'Neutral power to wall: ',il,hlpw
        end if
      end do
      close(lun)
c
c      if(natmi .gt. 1  .and.  npumplis .gt. 0) then
cc*** calculate the helium concentration at the pump
c        a=0.e0_R8
c        aa=0.e0_R8
c        do il=1,npumplis
c          is=lpumplis(il)
c          if(is .le. 0  .or.  is .gt. nlimi) then
c            print *,'=== Wrong segment number for the duct ',is
c            go to 999
c          end if
c          a=a+wldna(is,2,0)
c          do i=1,natmi
c            aa=aa+wldna(is,i,0)
c          end do
c          do i=1,nmoli
c            aa=aa+wldnm(is,i,0)*2.e0_R8
c          end do
c        end do
c        a=a/(a+aa)
c        print '(a,1p,e11.3)','He fraction to the pump: ',a
c        open(lun,file='fraction.he')
c        write(lun,'(1p,e16.3)') a
c        close(lun)
c      end if
      write(*,*) '*** Target loading with neutrals ***'
      write(*,700) 'Total power released',
     ,                               wldnek(nlim+7,0)+wldnep(nlim+7,0),
     ,                               wldnek(nlim+8,0)+wldnep(nlim+8,0)
      write(*,700) 'incl. potential energy',wldnep(nlim+7,0),
     ,                                                 wldnep(nlim+8,0)
      do i=1,natmi
        write(*,710) text_an(i),' atom flux',
     ,                              wldna(nlim+7,i,0),wldna(nlim+8,i,0)
        write(*,710) text_an(i),' average energy',
     ,                              ewlda(nlim+7,i,0),ewlda(nlim+8,i,0)
      end do
      if(ic .gt. 0) then
        wwyn=wldra(nlim+7,ic,0)-wldna(nlim+7,ic,0)
        wwyg=wldra(nlim+8,ic,0)-wldna(nlim+8,ic,0)
        write(*,700) 'Net   C erosion',wwyn,wwyg
        write(*,700) 'Gross C erosion',wldspt(nlim+7,0),wldspt(nlim+8,0)
      end if
c      print 700,'Hydrogen atom flux',wldna(nlim+7,1,0),wldna(nlim+8,1,0
c      print 700,'Average energy',ewlda(nlim+7,1,0),ewlda(nlim+8,1,0)
c      print 700,'Helium   atom flux',wldna(nlim+7,2,0),wldna(nlim+8,2,0
c      print 700,'Average energy',ewlda(nlim+7,2,0),ewlda(nlim+8,2,0)
c      if(ns .gt. 2) then
c          print 700,'Neon     atom flux',wldna(nlim+7,3,0),
c     ,                                                 wldna(nlim+8,3,0
c          print 700,'Average energy',ewlda(nlim+7,3,0),ewlda(nlim+8,3,0
c      end if
  700 format(1x,a22,':',1p,2e14.3)
  710 format(1x,a2,a20,':',1p,2e14.3)
      return
c----------------------------------------------------------------------
 999  return 1
c----------------------------------------------------------------------
      end subroutine b2pnwl
c
c
c======================================================================
      subroutine b2pwsrt(*)
c
c*** sort the lists of target sections to simply connected chains
c----------------------------------------------------------------------
c*** Global variables
c***
c***  nwldprt : the number of lists to be sorted
c***  lwldstrt: the initial segments of the chains
c***  nwldlis : the list lengths
c***  lwldlis : the list themselvs
c
      use b2mod_types
      use b2mod_b2plot
      implicit none
#     include "COUPLE/KOPPLDIM.F"
#     include "COUPLE/EIRDIAG.F"
#     include "b2pinpc.h"
c
c*** Local variables
c***
c*** goup  : the end of the last segment of the chain is still free
c*** placed: the segment is already in the chain
c*** lwsrtwrk: a temporary array for sorting
c*** nwsrtwrk: a counter for this array
c*** wsrttol : a tolerance parameter for defining adjacent segments
c
      logical goup,placed(nlim)
      integer lwsrtwrk(nlim),il,is,ii,jj,js,nwsrtwrk,k,ik,j
      real (kind=R8) ::
     *   wsrttol
      parameter (wsrttol=2.0e-4_R8)
c----------------------------------------------------------------------
c
cc%%%
c      print *,'b2pwsrt',nwldprt
c      do il=1,nwldprt
c        print '(1x,20i4)',il,nwldlis(il),lwldstrt(il)
c        print '(1x,20i4)',(lwldlis(i,il),i=1,nwldlis(il))
c        do ii=1,nwldlis(il)
c          is=lwldlis(ii,il)
c          print '(4f12.4)',pp1(1,is),pp1(2,is),pp2(1,is),pp2(2,is)
c        end do
c      end do
cc%%%
      do il=1,nwldprt
        is=lwldstrt(il)
        if(is .gt. nlim  .or.  is .le. 0) then
          write(*,*) 'Starting segment outside the limits',is,il
          nwldlis(il)=0.e0_R8
        else
          goup=.false.
          do ii=1,nwldlis(il)
            placed(ii)=.false.
            if(is .eq. lwldlis(ii,il)) then
              placed(ii)=.true.
              goup=.true.
            end if
          end do
          if(.not.goup) then
            write(*,*) 'Starting segment outside the list',is,il
            nwldlis(il)=0.e0_R8
          else
c*** initialise the chain counter and the chain list
            nwsrtwrk=1
            lwsrtwrk(nwsrtwrk)=is
c*** both ends of the first element will be checked (goup=.true. here)
            k=2
            do 100 ii=2,nwldlis(il)
c              print *,'do 100:',ii,is,goup
              do ik=1,k
c                print *,'do ik:',ik,is,goup
                do jj=1,nwldlis(il)
c                  print *,'do jj:',jj,is,goup,placed(jj)
                  if(.not.placed(jj)) then
                    js=lwldlis(jj,il)
c                    print *,'js:',js,is,goup
                    if(goup) then
                      if(abs(pp2(1,is)-pp1(1,js)) .le. wsrttol  .and. 
     .                 abs(pp2(2,is)-pp1(2,js)) .le. wsrttol) then
                        goup=.true.
                        is=js
                        nwsrtwrk=nwsrtwrk+1
                        lwsrtwrk(nwsrtwrk)=is
                        placed(jj)=.true.
c                        print *,'2->1:',is,nwsrtwrk
                        go to 100
                      end if
                      if(abs(pp2(1,is)-pp2(1,js)) .le. wsrttol  .and. 
     .                 abs(pp2(2,is)-pp2(2,js)) .le. wsrttol) then
                        goup=.false.
                        is=js
                        nwsrtwrk=nwsrtwrk+1
                        lwsrtwrk(nwsrtwrk)=is
                        placed(jj)=.true.
c                        print *,'2->2:',is,nwsrtwrk
                        go to 100
                      end if
                    else
                      if(abs(pp1(1,is)-pp2(1,js)) .le. wsrttol  .and. 
     .                 abs(pp1(2,is)-pp2(2,js)) .le. wsrttol) then
                        goup=.false.
                        is=js
                        nwsrtwrk=nwsrtwrk+1
                        lwsrtwrk(nwsrtwrk)=is
                        placed(jj)=.true.
c                        print *,'1->2:',is,nwsrtwrk
                        go to 100
                      end if
                      if(abs(pp1(1,is)-pp1(1,js)) .le. wsrttol  .and. 
     .                 abs(pp1(2,is)-pp1(2,js)) .le. wsrttol) then
                        goup=.true.
                        is=js
                        nwsrtwrk=nwsrtwrk+1
                        lwsrtwrk(nwsrtwrk)=is
                        placed(jj)=.true.
c                        print *,'1->1:',is,nwsrtwrk
                        go to 100
                      end if
                    end if
                  end if
                end do
c*** This only works for the first element.
c*** Look in the opposite direction
                goup=.false.
              end do
c*** No connection found
              go to 50
 100        k=1
 50         if(nwsrtwrk .ne. nwldlis(il)) then
              write(*,*) 'The chain is broken: ',il
              nwldlis(il)=nwsrtwrk
            end if
            do ii=1,nwldlis(il)
              j=lwsrtwrk(ii)
              if(abs(p1(1,j)-pp1(1,j)) .le. wsrttol  .and. 
     .         abs(p1(2,j)-pp1(2,j)) .le. wsrttol  .and. 
     .         abs(p2(1,j)-pp2(1,j)) .le. wsrttol  .and. 
     .         abs(p2(2,j)-pp2(2,j)) .le. wsrttol) then
                lwldlis(ii,il)=j
              else
                lwldlis(ii,il)=-j
              end if
            end do
c            print '(a,i2/(1x,20i4))','il =',il,
c     ,                                (lwldlis(ii,il),ii=1,nwldlis(il))
          end if
        end if
      end do
c----------------------------------------------------------------------
      return
      end subroutine b2pwsrt
