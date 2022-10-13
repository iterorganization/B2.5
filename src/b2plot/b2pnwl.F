      subroutine b2pnwl(iret)
c
c  version : 11.11.97 20:52
c
c======================================================================
c*** b2pnwl: neutral wall load (temporary solution)
c***
c======================================================================
      use b2mod_types
      use b2mod_geo
      use b2mod_b2plot
      use b2mod_b2cmfs
      use b2mod_indirect
      use b2mod_constants
      use b2mod_neutrals_namelist
      use b2mod_b2plot_wall_loading
#ifdef B25_EIRENE
      use eirmod_wneutrals
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
      integer il,is,i,ii,ia,ib,k,ic,ias,ix,iy,iss
      character*74 hlp_frm
      real (kind=R8) ::
     *   x,distance,h,a,aa,hlpx,hlpw,wwt,wwk,wwp,wwyn,wwyg,wwy(4),
     ,     wwj(natm+nmol),wwe(natm+nmol)
      real (kind=R8) :: hlp_w(4),hlp_p(4),hlp_s(4)
      real (kind=R8), allocatable :: hlp_a(:,:),hlp_e(:,:),
     ,                               hlp_m(:,:),hlp_f(:,:)
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
      call b2pwlprp
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
            if (is.gt.0.and.is.le.nlim) then
              distance=sqrt((p1(1,is)-p2(1,is))**2+
     .                      (p1(2,is)-p2(2,is))**2)
              hlpx=hlpx+distance
            else if (is.gt.nlim) then
              ias=is-nlim
              do iss=1,nsts !{
                if (eirdiag_nds_typ(iss).ne.2) cycle
                if (ias.lt.eirdiag_nds_ind(iss)) cycle
                if (ias.gt.eirdiag_nds_ind(iss+1)) cycle
                ix=eirdiag_nds_srf(iss)
                iy=ias+(eirdiag_nds_start(iss)-1)-eirdiag_nds_ind(iss)
                hlpx=hlpx+distc(ix,iy,0,ix,iy,2)
              end do !}
            end if
          end do
c
          x=0.e0_R8
          distance=0.e0_R8
          hlpw=0.e0_R8
          do i=1,nwldlis(il)
            ia=lwldlis(i,il)
            ib=iabs(ia)
            is=leirlis(i,il)
            h=distance
            distance=sqrt((pp1(1,ib)-pp2(1,ib))**2+
     .                    (pp1(2,ib)-pp2(2,ib))**2)
            if (isymm.eq.0) then
              a=width*distance
            else if (isymm.eq.1 .or. isymm.eq.2) then
              a=pi*(pp1(1,ib)+pp2(1,ib))*distance
            else if (isymm.eq.3 .or. isymm.eq.4) then
              a=pi*(pp1(2,ib)+pp2(2,ib))*distance
            end if
            aa=1.e-6_R8/a
            distance=0.5e0_R8*distance
            x=x+h+distance
            if(is.gt.0 .and. is.le.nlim) then
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
            else if (is.gt.nlim) then
              ias=is-nlim
              wwk=0.0_R8
              wwp=0.0_R8
              k=0
              do ii=1,natmi
                k=k+1
                wwk=wwk+(ewlda_res(ii,ias)-ewldea_res(ii,ias))
                wwj(k)=max(wldna_res(ii,ias),flxctff)/a
                wwe(k)=max(ewlda_res(ii,ias),enrctff)
              end do
              do ii=1,nmoli
                k=k+1
                wwk=wwk+(ewldm_res(ii,ias)-ewldem_res(ii,ias))
                wwp=wwp+ewldmr_res(ii,ias)
                wwj(k)=max(wldnm_res(ii,ias),flxctff)/a
                wwe(k)=max(ewldm_res(ii,ias),enrctff)
              end do
              if(ic.gt.0) then
                wwyg=wldspta_res(ic,ias)/a
                wwyn=(wldspta_res(ic,ias)-wldna_res(ic,ias)-
     >                                    wlpump_res(ic,ias))/a
                do ii=1,nmoli
                  if (mlcmp(ic,ii).gt.0)
     >             wwyn=wwyn+mlcmp(ic,ii)*
     .              (wldsptm_res(ii,ias)-wldnm_res(ii,ias)-
     .                                   wlpump(ii+natmi,ias))/a
                end do
              else
                wwyn=flxctff
                wwyg=flxctff
              end if
              hlpw=hlpw+wwk+wwp
              wwk=wwk*aa
              wwp=wwp*aa
              wwt=max(wwk+wwp,pwrctff)
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
      allocate(hlp_a(ntrgsx,natmi))
      allocate(hlp_e(ntrgsx,natmi))
      allocate(hlp_m(ntrgsx,nmoli))
      allocate(hlp_f(ntrgsx,nmoli))
      hlp_a=0.0_R8
      hlp_m=0.0_R8
      hlp_e=0.0_R8
      hlp_f=0.0_R8
      hlp_w=0.0_R8
      hlp_p=0.0_R8
      do k=1,ntns
        if (ktns(k).eq.0) cycle
        if (itrg(ktns(k)).eq.0) cycle
        hlp_w(itrg(ktns(k))) = hlp_w(itrg(ktns(k)))+
     .                   wldnek(nlim+ltns(k),0)+wldnep(nlim+ltns(k),0)
        hlp_p(itrg(ktns(k))) = hlp_p(itrg(ktns(k)))+
     .                   wldnep(nlim+ltns(k),0)
        do i=1,natmi
          hlp_a(itrg(ktns(k)),i) = hlp_a(itrg(ktns(k)),i)+
     .                   wldna(nlim+ltns(k),i,0)
          hlp_e(itrg(ktns(k)),i) = hlp_e(itrg(ktns(k)),i)+
     .                   ewlda(nlim+ltns(k),i,0)*wldna(nlim+ltns(k),i,0)
        end do
        do i=1,nmoli
          hlp_m(itrg(ktns(k)),i) = hlp_m(itrg(ktns(k)),i)+
     .                   wldnm(nlim+ltns(k),i,0)
          hlp_f(itrg(ktns(k)),i) = hlp_f(itrg(ktns(k)),i)+
     .                   ewldm(nlim+ltns(k),i,0)*wldnm(nlim+ltns(k),i,0)
        end do
      end do
      do k=1,ntrgts
        if (itrg(k).eq.0) cycle
        do i=1,natmi
          if (hlp_a(itrg(k),i).ne.0.0_R8)
     >     hlp_e(itrg(k),i) = hlp_e(itrg(k),i)/hlp_a(itrg(k),i)
        end do
        do i=1,nmoli
          if (hlp_m(itrg(k),i).ne.0.0_R8)
     >     hlp_f(itrg(k),i) = hlp_f(itrg(k),i)/hlp_m(itrg(k),i)
        end do
      end do
      write(*,700) 'Total power released',(hlp_w(itrg(i)),i=1,ntrgts)
      write(*,700) 'incl. potential energy',
     ,                                    (hlp_p(itrg(i)),i=1,ntrgts)
      do i=1,natmi
        write(*,710) text_an(i),' atom flux',
     ,                                    (hlp_a(itrg(k),i),k=1,ntrgts)
        write(*,710) text_an(i),' average energy',
     ,                                    (hlp_e(itrg(k),i),k=1,ntrgts)
      end do
      do i=1,nmoli
        write(*,710) text_mn(i),' molecule flux',
     ,                                    (hlp_m(itrg(k),i),k=1,ntrgts)
        write(*,710) text_mn(i),' average energy',
     ,                                    (hlp_f(itrg(k),i),k=1,ntrgts)
      end do
      if(ic .gt. 0) then
        wwy=0.0_R8
        hlp_s=0.0_R8
        do k=1,ntns
          if (ktns(k).eq.0) cycle
          if (itrg(ktns(k)).eq.0) cycle
          wwy(itrg(ktns(k)))=wwy(itrg(ktns(k)))+
     .                 wldra(nlim+ltns(k),ic,0)+wldpa(nlim+ltns(k),ic,0)
     .                -wldna(nlim+ltns(k),ic,0)
          do i=1,nmoli
            if (mlcmp(ic,i).gt.0) then
              wwy(itrg(ktns(k)))=wwy(itrg(ktns(k)))+mlcmp(ic,i)*
     .                            (wldrm(nlim+ltns(k),i,0)
     .                            +wldpm(nlim+ltns(k),i,0)
     .                            -wldnm(nlim+ltns(k),i,0))
            end if
          end do
          hlp_s(itrg(ktns(k))) = hlp_s(itrg(ktns(k)))+
     .                          wldspt(nlim+ltns(k),0)
        end do
        write(*,700) 'Net   C erosion',(wwy(itrg(i)),i=1,ntrgts)
        write(*,700) 'Gross C erosion',(hlp_s(itrg(i)),i=1,ntrgts)
      end if
      deallocate(hlp_a,hlp_e,hlp_m,hlp_f)
c
c*** Calculate the load
c
      call b2pwrdld

  700 format(1x,a22,':',1p,4e14.3)
  710 format(1x,a2,a20,':',1p,4e14.3)
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
      logical goup,placed(DEF_NLIM),found,direction_same
      integer lwsrtwrk(DEF_NLIM),il,is,ii,jj,js,nwsrtwrk,k,ik,j
      real (kind=R8) ::
     *   wsrttol, d2etol
      parameter (wsrttol=2.0e-4_R8, d2etol=1.0e-2_R8)
c----------------------------------------------------------------------
c
cc%%%
c      print *,'b2pwsrt',nwldprt
c      if (nnreg(0).lt.7) then
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
              else if (nlimi.eq.0) then
                lwldlis(ii,il)=j
                leirlis(ii,il)=j
              else if (j.eq.jedge(1).or.j.eq.jedge(2).or.
     .                 j.eq.jedge(3).or.j.eq.jedge(4).or.
     .               ((j.eq.jedge(5).or.j.eq.jedge(6).or.
     .                 j.eq.jedge(7).or.j.eq.jedge(8)).and.
     .                 nnreg(0).ge.7)) then
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
