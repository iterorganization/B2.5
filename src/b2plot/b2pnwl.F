      subroutine b2pnwl(iret)
c
c  version : 11.11.97 20:52
c
c======================================================================
c*** b2pnwl: neutral wall load (temporary solution)
c***
c======================================================================
      use b2mod_types
      use b2mod_b2plot
      use b2mod_b2cmfs
      use b2mod_indirect
      use b2mod_constants
      use b2mod_neutrals_namelist
      use b2mod_b2plot_wall_loading
#ifdef B25_EIRENE
      use eirmod_extrab25
#endif
      implicit none
#include <DIMENSIONS.F>
      integer iret
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
      integer il,is,i,ii,ia,ib,k,ic
      character*74 hlp_frm
      real (kind=R8) ::
     *   x,dist,h,a,aa,hlpx,hlpw,wwt,wwk,wwp,wwyn,wwyg,
     ,     wwj(natm+nmol),wwe(natm+nmol)
      real (kind=R8) :: width
      data width /1.0_R8/
      save width
      character*8 text_an(natm),text_mn(nmol)
      logical streql
      external streql, ipgetr
c----------------------------------------------------------------------
      print *,'Starting b2pnwl...',nlim,DEF_NLIM
c*** Input the DG data and prepare the geometry
c
      call ipgetr ('b2agmt_1d_width', width)
c*** Set the default values
      npumplis=0
      nwldprt=0
c*** Read the data from dg files
      iret=0
      call b2pinp(iret)
      if (iret.eq.1) goto 999
c*** Re-order jedge data if less than ntrgsx targets
      if (ntrgsa.lt.ntrgsx) then
        do i=1,ntrgsa
          jedge(ntrgsa+i)=jedge(ntrgsx+i)
          jedge(ntrgsx+i)=0
        end do
      end if
c
c*** Convert the input wall coordinates into double precision
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
c      write (hlp_frm,'(a,i2,a)') '(a,',nwldprt,'i5)'
c      write (6,hlp_frm) 'lwldstrt:',(lwldstrt(i),i=1,nwldprt)
c      write (6,hlp_frm) 'nwldlis :',(nwldlis(i),i=1,nwldprt)
c      do il=1,nwldprt
c        print '(a,i5/(5x,15i5))','lwldlis :',il,
c     ,                                  (lwldlis(i,il),i=1,nwldlis(il))
c      end do
cc%%%
c*** Then sort the target lists
c*** The shadowing structure must be closed, therefore it starts from
c*** an arbitrary point.

      lwldstrt(1)=iabs(lwldlis(1,1))
      call b2pwsrt(iret)
      if (iret.eq.1) goto 999
cc<<<
c      print *,'After b2pwsrt',nwldprt
c      write (hlp_frm,'(a,i2,a)') '(a,',nwldprt,'i5)'
c      write (6,hlp_frm) 'lwldstrt:',(lwldstrt(i),i=1,nwldprt)
c      write (6,hlp_frm) 'nwldlis :',(nwldlis(i),i=1,nwldprt)
c      do il=1,nwldprt
c        print '(a,i5/(5x,15i5))','lwldlis :',il,
c     ,                                  (lwldlis(i,il),i=1,nwldlis(il))
c        print '(a,i5/(5x,15i5))','leirlis :',il,
c     ,                                  (leirlis(i,il),i=1,nwldlis(il))
c      end do
cc>>>
c*** ... and produce the output file
      do i=1,natmi
        text_an(i)=textan(i-1)(2:8)
      end do
      do i=1,nmoli
        text_mn(i)=textmn(i-1)(2:8)
      end do
      open (lun,file='wlld.dat')
      write (hlp_frm,'(a,i2,a)')
     .      '(2x,a4,5x,a3,1x,a4,6(5x,a2,3x),1x,4(5x,a4,3x),3x,',
     .      2*natmi+2*nmoli,
     .      '(a4,a8),2(1x,a8,3x))'
      if (nmoli.gt.0) then
        write(lun,hlp_frm)
     -     'zone','DG#','EIR#',
     ,     'R1','Z1','R2','Z2',' X','Xr',
     ,     'area','Wtot','Wkin','Wpot',
     ,     ('flx_',text_an(i),'Eav_',text_an(i),i=1,natmi),
     ,     ('flx_',text_mn(i),'Eav_',text_mn(i),i=1,nmoli),
     ,      'net_ersn','grs_ersn'
      else
        write(lun,hlp_frm)
     -     'zone','DG#','EIR#',
     ,     'R1','Z1','R2','Z2',' X','Xr',
     ,     'area','Wtot','Wkin','Wpot',
     ,     ('flx_',text_an(i),'Eav_',text_an(i),i=1,natmi),
     ,      'net_ersn','grs_ersn'
      end if
#ifdef DBG
      print *,'nwldprt,natmi,nmoli',nwldprt,natmi,nmoli
      do il=1,nwldprt
        print '(1x,2i4,1x,a)',il,nwldlis(il),pltzlbl(il)
        write (6,'(1x,20i4)') (lwldlis(i,il),i=1,nwldlis(il))
      end do
#endif
c
c*** Check whether carbon atoms are present in the calculation
c
      ic=0
      do i=1,natmi
        if(streql(text_an(i),'C')) ic=i
      end do
c
      do il=2,nwldprt
        if(nwldlis(il) .gt. 0) then
          write(lun,*)
c*** calculate the length of the target segment
          hlpx=0.e0_R8
          do i=1,nwldlis(il)
            ia=lwldlis(i,il)
            is=leirlis(i,il)
            if (is.gt.0) then
              dist=sqrt((p1(1,is)-p2(1,is))**2+(p1(2,is)-p2(2,is))**2)
              hlpx=hlpx+dist
            end if
          end do
c
          x=0.e0_R8
          dist=0.e0_R8
          hlpw=0.e0_R8
          do i=1,nwldlis(il)
            ia=lwldlis(i,il)
            ib=iabs(ia)
            is=leirlis(i,il)
            h=dist
            dist=sqrt((pp1(1,ib)-pp2(1,ib))**2+
     .                (pp1(2,ib)-pp2(2,ib))**2)
            if (isymm.eq.0) then
              a=width*dist
            else if (isymm.eq.1 .or. isymm.eq.2) then
              a=pi*(pp1(1,ib)+pp2(1,ib))*dist
            else if (isymm.eq.3 .or. isymm.eq.4) then
              a=pi*(pp1(2,ib)+pp2(2,ib))*dist
            end if
            aa=1.e-6_R8/a
            dist=0.5e0_R8*dist
            x=x+h+dist
            if(is.gt.0) then
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
              if(ic .gt. 0) then
                wwyg=wldspt(is,0)/a
                wwyn=(wldra(is,ic,0)-wldna(is,ic,0)
     .                              -wldpa(is,ic,0))/a
                do ii=1,nmoli
                  if (mlcmp(ic,ii).gt.0)
     >             wwyn=wwyn+mlcmp(ic,ii)*(wldrm(is,ii,0)-wldnm(is,ii,0)
     .                                    -wldpm(is,ii,0))/a
                end do
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
            write (hlp_frm,'(a,i2,a)')
     ,                 '(1x,a8,2i5,1p,5f10.4,',2*k+12,'e12.4)'
            write (lun,hlp_frm) pltzlbl(il),
     ,                ib,is,pp1(1,ib),pp1(2,ib),pp2(1,ib),pp2(2,ib),
     ,                x,x/hlpx,a,wwt,wwk,wwp,(wwj(ii),wwe(ii),ii=1,k),
     ,                wwyn,wwyg
          end do
          write(*,*) 'Neutral power to wall: ',pltzlbl(il),il,hlpw
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
c          a=a+wldna(is,j_he_at,0)
c          do i=1,natmi
c            aa=aa+wldna(is,i,0)
c          end do
c          do i=1,nmoli
c            do j=1,natmi
c              aa=aa+wldnm(is,i,0)*mlcmp(j,i)
c            end do
c          end do
c        end do
c        a=a/(a+aa)
c        print '(a,1p,e11.3)','He fraction to the pump: ',a
c        open(lun,file='fraction.he')
c        write(lun,'(1p,e16.3)') a
c        close(lun)
c      end if
      write(*,*) '*** Target loading with neutrals ***'
      if (ltns(2).eq.0) then
       write(*,700) 'Total power released',
     ,              wldnek(nlim+ltns(1),0)+wldnep(nlim+ltns(1),0)
       write(*,700) 'incl. potential energy',wldnep(nlim+ltns(1),0)
       do i=1,natmi
        write(*,710) text_an(i),' atom flux',wldna(nlim+ltns(1),i,0)
        write(*,710) text_an(i),' average energy',
     ,                                       ewlda(nlim+ltns(1),i,0)
       end do
       if(ic .gt. 0) then
        wwyn=wldra(nlim+ltns(1),ic,0)+wldpa(nlim+ltns(1),ic,0)
     .      -wldna(nlim+ltns(1),ic,0)
        do i=1,nmoli
          if (mlcmp(ic,i).gt.0) then
            wwyn=wwyn+mlcmp(ic,i)*(wldrm(nlim+ltns(1),i,0)
     .                            +wldpm(nlim+ltns(1),i,0)
     .                            -wldnm(nlim+ltns(1),i,0))
          end if
        end do
        write(*,700) 'Net   C erosion',wwyn
        write(*,700) 'Gross C erosion',wldspt(nlim+ltns(1),0)
       end if
      else
       write(*,700) 'Total power released',
     ,              wldnek(nlim+ltns(1),0)+wldnep(nlim+ltns(1),0),
     ,              wldnek(nlim+ltns(2),0)+wldnep(nlim+ltns(2),0)
       write(*,700) 'incl. potential energy',wldnep(nlim+ltns(1),0),
     ,                                       wldnep(nlim+ltns(2),0)
       do i=1,natmi
        write(*,710) text_an(i),' atom flux',
     ,               wldna(nlim+ltns(1),i,0),wldna(nlim+ltns(2),i,0)
        write(*,710) text_an(i),' average energy',
     ,               ewlda(nlim+ltns(1),i,0),ewlda(nlim+ltns(2),i,0)
       end do
       if(ic .gt. 0) then
        wwyn=wldra(nlim+ltns(1),ic,0)+wldpa(nlim+ltns(1),ic,0)
     .      -wldna(nlim+ltns(1),ic,0)
        wwyg=wldra(nlim+ltns(2),ic,0)+wldpa(nlim+ltns(2),ic,0)
     .      -wldna(nlim+ltns(2),ic,0)
        do i=1,nmoli
          if (mlcmp(ic,i).gt.0) then
            wwyn=wwyn+mlcmp(ic,i)*(wldrm(nlim+ltns(1),i,0)
     .                            +wldpm(nlim+ltns(1),i,0)
     .                            -wldnm(nlim+ltns(1),i,0))
            wwyg=wwyg+mlcmp(ic,i)*(wldrm(nlim+ltns(2),i,0)
     .                            +wldpm(nlim+ltns(2),i,0)
     .                            -wldnm(nlim+ltns(2),i,0))
          end if
        end do
        write(*,700) 'Net   C erosion',wwyn,wwyg
        write(*,700) 'Gross C erosion',wldspt(nlim+ltns(1),0),
     .                                 wldspt(nlim+ltns(2),0)
       end if
      end if
c     if (j_h_at(1).ne.0) then
c      print 700,'Hydrogen  atom flux',wldna(nlim+ltns(1),j_h_at(1),0),
c                                      wldna(nlim+ltns(2),j_h_at(1),0)
c      print 700,'Average energy',ewlda(nlim+ltns(1),j_h_at(1),0),
c                                 ewlda(nlim+ltns(2),j_h_at(1),0)
c     end if
c     if (j_h_at(2).ne.0) then
c      print 700,'Deuterium atom flux',wldna(nlim+ltns(1),j_h_at(2),0),
c                                      wldna(nlim+ltns(2),j_h_at(2),0)
c      print 700,'Average energy',ewlda(nlim+ltns(1),j_h_at(2),0),
c                                 ewlda(nlim+ltns(2),j_h_at(2),0)
c     end if
c     if (j_h_at(3).ne.0) then
c      print 700,'Tritium   atom flux',wldna(nlim+ltns(1),j_h_at(3),0),
c                                      wldna(nlim+ltns(2),j_h_at(3),0)
c      print 700,'Average energy',ewlda(nlim+ltns(1),j_h_at(3),0),
c                                 ewlda(nlim+ltns(2),j_h_at(3),0)
c     end if
c     if (j_he_at.ne.0) then
c      print 700,'Helium    atom flux',wldna(nlim+ltns(1),j_he_at,0),
c                                      wldna(nlim+ltns(2),j_he_at,0)
c      print 700,'Average energy',ewlda(nlim+ltns(1),j_he_at,0),
c                                 ewlda(nlim+ltns(2),j_he_at,0)
c     end if
c     if(j_ne_at.ne.0) then
c      print 700,'Neon      atom flux',wldna(nlim+ltns(1),j_ne_at,0),
c     ,                                wldna(nlim+ltns(2),j_ne_at,0)
c      print 700,'Average energy',ewlda(nlim+ltns(1),j_ne_at,0),
c                                 ewlda(nlim+ltns(2),j_ne_at,0)
c     end if
  700 format(1x,a22,':',1p,2e14.3)
  710 format(1x,a2,a20,':',1p,2e14.3)
      return
c----------------------------------------------------------------------
  999 return
c----------------------------------------------------------------------
      end subroutine b2pnwl
c
c
c======================================================================
      subroutine b2pwsrt(iret)
c
c*** sort the lists of target sections to simply connected chains
c----------------------------------------------------------------------
c*** Global variables
c***
c***  nwldprt : the number of lists to be sorted
c***  lwldstrt: the initial segments of the chains
c***  nwldlis : the list lengths
c***  lwldlis : the lists themselves (DG numbering)
c***  leirlis : the lists themselves (Eirene numbering)
c
      use b2mod_types
      use b2mod_b2plot
      use b2mod_indirect
      use b2mod_b2plot_wall_loading
#ifdef B25_EIRENE
      use eirmod_parmmod
      use eirmod_extrab25
#else
      use b2mod_eirdiag
#endif
      implicit none
      integer iret
#include <DIMENSIONS.F>
#ifndef B25_EIRENE
#     include "COUPLE/KOPPLDIM.F"
#endif
c
c*** Local variables
c***
c*** goup  : the end of the last segment of the chain is still free
c*** placed: the segment is already in the chain
c*** lwsrtwrk: a temporary array for sorting
c*** nwsrtwrk: a counter for this array
c*** wsrttol : a tolerance parameter for defining adjacent segments
c
      logical goup,placed(DEF_NLIM),found,direction_same,no_fort44
      integer lwsrtwrk(DEF_NLIM),il,is,ii,jj,js,nwsrtwrk,k,ik,j
      real (kind=R8) ::
     *   wsrttol, d2etol
      parameter (wsrttol=2.0e-4_R8, d2etol=1.0e-2_R8)
c----------------------------------------------------------------------
c
cc%%%
c      print *,'b2pwsrt',nwldprt
c      if (nnreg(0).lt.8) then
c        print '(a,4i6)','jedgi1, jedgi2, jedgo1, jedgo2 =',jedge(1:4)
c      else
c        print '(a,4i6)','jedgi1, jedgi2, jedgi3, jedgi4 =',jedge(1:4)
c        print '(a,4i6)','jedgo1, jedgo2, jedgo3, jedgo4 =',jedge(5:8)
c      end if
c      do il=1,nwldprt
c        print '(1x,20i4)',il,nwldlis(il),lwldstrt(il)
c        print '(1x,20i4)',(lwldlis(ii,il),ii=1,nwldlis(il))
c        do ii=1,nwldlis(il)
c          is=lwldlis(ii,il)
c          print '(4f12.4)',pp1(1,is),pp1(2,is),pp2(1,is),pp2(2,is)
c        end do
c      end do
cc%%%
      no_fort44=maxval(p1).eq.0.0_R8 .and. minval(p1).eq.0.0_R8 .and.
     &          maxval(p2).eq.0.0_R8 .and. minval(p2).eq.0.0_R8
      iret=0
      do il=1,nwldprt
        is=lwldstrt(il)
        if(is .gt. lppp  .or.  is .le. 0) then
          write(*,*) 'Starting segment outside the limits',is,il
          nwldlis(il)=0
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
            nwldlis(il)=0
          else
c*** initialise the chain counter and the chain list
            nwsrtwrk=1
            lwsrtwrk(nwsrtwrk)=is
c*** both ends of the first element will be checked (goup=.true. here)
            k=2
            do ii=2,nwldlis(il)
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
     .                   abs(pp2(2,is)-pp1(2,js)) .le. wsrttol) then
                        goup=.true.
                        is=js
                        nwsrtwrk=nwsrtwrk+1
                        lwsrtwrk(nwsrtwrk)=is
                        placed(jj)=.true.
c                        print *,'2->1:',is,nwsrtwrk
                        go to 100
                      end if
                      if(abs(pp2(1,is)-pp2(1,js)) .le. wsrttol  .and.
     .                   abs(pp2(2,is)-pp2(2,js)) .le. wsrttol) then
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
     .                   abs(pp1(2,is)-pp2(2,js)) .le. wsrttol) then
                        goup=.false.
                        is=js
                        nwsrtwrk=nwsrtwrk+1
                        lwsrtwrk(nwsrtwrk)=is
                        placed(jj)=.true.
c                        print *,'1->2:',is,nwsrtwrk
                        go to 100
                      end if
                      if(abs(pp1(1,is)-pp1(1,js)) .le. wsrttol  .and.
     .                   abs(pp1(2,is)-pp1(2,js)) .le. wsrttol) then
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
  100         k=1
            end do
   50       if(nwsrtwrk .ne. nwldlis(il)) then
              write(*,*) 'The chain is broken: ',il
              nwldlis(il)=nwsrtwrk
            end if
            do ii=1,nwldlis(il)
              j = lwsrtwrk(ii)
              if (lskip(j)) then
                lwldlis(ii,il)=j
                leirlis(ii,il)=0
              else if (no_fort44) then
                lwldlis(ii,il)=j
                leirlis(ii,il)=j
              else if (j.eq.jedge(1).or.j.eq.jedge(2).or.
     .                 j.eq.jedge(3).or.j.eq.jedge(4).or.
     .               ((j.eq.jedge(5).or.j.eq.jedge(6).or.
     .                 j.eq.jedge(7).or.j.eq.jedge(8)).and.
     .                 nnreg(0).eq.8)) then
                k = 1
                found = .false.
                do while (.not.found .and. k.le.nlim)
                  if(abs(p1(1,k)-pp1(1,j)) .le. wsrttol  .and.
     .               abs(p1(2,k)-pp1(2,j)) .le. wsrttol) then
                    direction_same = abs(
     .               (((p2(1,k)-p1(1,k))*(pp2(1,j)-pp1(1,j))+
     .                 (p2(2,k)-p1(2,k))*(pp2(2,j)-pp1(2,j)))/
     .                (sqrt(( p2(1,k)- p1(1,k))**2+
     .                      ( p2(2,k)- p1(2,k))**2)*
     .                 sqrt((pp2(1,j)-pp1(1,j))**2+
     .                      (pp2(2,j)-pp1(2,j))**2)))-1.0_R8).le.d2etol
                  else
     .            if(abs(p1(1,k)-pp2(1,j)) .le. wsrttol  .and.
     .               abs(p1(2,k)-pp2(2,j)) .le. wsrttol) then
                    direction_same = abs(
     .               (((p2(1,k)-p1(1,k))*(pp1(1,j)-pp2(1,j))+
     .                 (p2(2,k)-p1(2,k))*(pp1(2,j)-pp2(2,j)))/
     .                (sqrt(( p2(1,k)- p1(1,k))**2+
     .                      ( p2(2,k)- p1(2,k))**2)*
     .                 sqrt((pp2(1,j)-pp1(1,j))**2+
     .                      (pp2(2,j)-pp1(2,j))**2)))-1.0_R8).le.d2etol
                  else
     .            if(abs(p2(1,k)-pp1(1,j)) .le. wsrttol  .and.
     .               abs(p2(2,k)-pp1(2,j)) .le. wsrttol) then
                    direction_same = abs(
     .               (((p1(1,k)-p2(1,k))*(pp2(1,j)-pp1(1,j))+
     .                 (p1(2,k)-p2(2,k))*(pp2(2,j)-pp1(2,j)))/
     .                (sqrt(( p2(1,k)- p1(1,k))**2+
     .                      ( p2(2,k)- p1(2,k))**2)*
     .                 sqrt((pp2(1,j)-pp1(1,j))**2+
     .                      (pp2(2,j)-pp1(2,j))**2)))-1.0_R8).le.d2etol
                  else
     .            if(abs(p2(1,k)-pp2(1,j)) .le. wsrttol  .and.
     .               abs(p2(2,k)-pp2(2,j)) .le. wsrttol) then
                    direction_same = abs(
     .               (((p1(1,k)-p2(1,k))*(pp1(1,j)-pp2(1,j))+
     .                 (p1(2,k)-p2(2,k))*(pp1(2,j)-pp2(2,j)))/
     .                (sqrt(( p2(1,k)- p1(1,k))**2+
     .                      ( p2(2,k)- p1(2,k))**2)*
     .                 sqrt((pp2(1,j)-pp1(1,j))**2+
     .                      (pp2(2,j)-pp1(2,j))**2)))-1.0_R8).le.d2etol
                  else
                    direction_same = .false.
                  end if
                  if (direction_same) then
                    found = .true.
                    lwldlis(ii,il)=j
                    leirlis(ii,il)=k
                  else
                    k = k + 1
                  end if
                end do
                if (.not.found) write(*,'(a,i4,a)')
     .            'Found no correspondence (1A) between DG wall # ', j ,
     .            ' and Eirene wall list'
              else
                k = 1
                found = .false.
                do while (.not.found .and. k.le.nlim)
                  if((abs(p1(1,k)-pp1(1,j)) .le. wsrttol  .and.
     .                abs(p1(2,k)-pp1(2,j)) .le. wsrttol  .and.
     .                abs(p2(1,k)-pp2(1,j)) .le. wsrttol  .and.
     .                abs(p2(2,k)-pp2(2,j)) .le. wsrttol) .or.
     .               (abs(p1(1,k)-pp2(1,j)) .le. wsrttol  .and.
     .                abs(p1(2,k)-pp2(2,j)) .le. wsrttol  .and.
     .                abs(p2(1,k)-pp1(1,j)) .le. wsrttol  .and.
     .                abs(p2(2,k)-pp1(2,j)) .le. wsrttol)) then
                    lwldlis(ii,il)=j
                    leirlis(ii,il)=k
                    found = .true.
                  end if
                  k = k + 1
                end do
                if (.not.found) write(*,'(a,i4,a)')
     .           'Found no correspondence (2A) between DG wall # ', j ,
     .           ' and Eirene wall list'
              end if
            end do
c            print '(a,i2/(1x,20i4))','il =',il,
c     ,                                (lwldlis(ii,il),ii=1,nwldlis(il))
c            print '(1x,20i4)', (leirlis(ii,il),ii=1,nwldlis(il))
          end if
        end if
      end do
c----------------------------------------------------------------------
      return
      end subroutine b2pwsrt
