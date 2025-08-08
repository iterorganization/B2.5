      subroutine b2pnwl(iret,from_b2plot)
c
c  version : 11.11.97 20:52
c
c======================================================================
c*** b2pnwl: neutral wall load (temporary solution)
c***
c======================================================================
      use b2us_data
      use b2mod_types
      use b2mod_b2plot
      use b2mod_b2cmfs
      use b2mod_b2cmpa
      use b2mod_constants
      use b2mod_neutrals_namelist
      use b2mod_b2plot_wall_loading
#ifdef B25_EIRENE
      use eirmod_wneutrals
#endif
      use b2mod_dimensions
      implicit none
      integer iret
      logical from_b2plot
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
      parameter (flxctff=0.1_R8,enrctff=1.0e-8_R8,pwrctff=1.0e-8_R8)
      integer il,is,i,ii,ia,ib,j,js,k,ic,ias,iend,istart,iFc,iCv
      character*74 hlp_frm
      real (kind=R8) ::
     *   x,distance,h,a,aa,u,hlpx,hlpw,wwt,wwk,wwp,wwyn,wwyg,wwy(4),
     ,     wwj(DEF_NATM+DEF_NMOL),wwe(DEF_NATM+DEF_NMOL)
      real (kind=R8) :: hlp_w(4),hlp_p(4),hlp_s(4)
      real (kind=R8), allocatable :: hlp_a(:,:),hlp_e(:,:),
     ,                               hlp_m(:,:),hlp_f(:,:)
      real (kind=R8) :: width
      data width /1.0_R8/
      save width
      character*8 text_an(DEF_NATM),text_mn(DEF_NMOL)
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
      call b2pwlprp
      if (iret.eq.1) goto 999
c*** ... and produce the output file
      do i=1,natmi
        text_an(i)=textan(i-1)(2:8)
      end do
      do i=1,nmoli
        text_mn(i)=textmn(i-1)(2:8)
      end do
c*** Check whether carbon atoms are present in the calculation
      ic=0
      do i=1,natmi
        if(streql(text_an(i),'C')) ic=i
      end do
c
      if (from_b2plot) then
        open (lun,file='wlld.dat')
        if (ic.eq.0) then
          write (hlp_frm,'(a,i2,a)')
     .      '(2x,a4,5x,a3,1x,a4,6(5x,a2,3x),1x,4(5x,a4,3x),3x,',
     .      2*natmi+2*nmoli,
     .      '(a4,a8))'
        else
          write (hlp_frm,'(a,i2,a)')
     .      '(2x,a4,5x,a3,1x,a4,6(5x,a2,3x),1x,4(5x,a4,3x),3x,',
     .      2*natmi+2*nmoli,
     .      '(a4,a8),2(1x,a8,3x))'
        end if
        if (ic.eq.0) then
          if (nmoli.gt.0) then
            write(lun,hlp_frm)
     -       'zone','DG#','EIR#',
     ,       'R1','Z1','R2','Z2',' X','Xr',
     ,       'area','Wtot','Wkin','Wpot',
     ,       ('flx_',text_an(i),'Eav_',text_an(i),i=1,natmi),
     ,       ('flx_',text_mn(i),'Eav_',text_mn(i),i=1,nmoli)
          else
            write(lun,hlp_frm)
     -       'zone','DG#','EIR#',
     ,       'R1','Z1','R2','Z2',' X','Xr',
     ,       'area','Wtot','Wkin','Wpot',
     ,       ('flx_',text_an(i),'Eav_',text_an(i),i=1,natmi)
          end if
        else
          if (nmoli.gt.0) then
            write(lun,hlp_frm)
     -       'zone','DG#','EIR#',
     ,       'R1','Z1','R2','Z2',' X','Xr',
     ,       'area','Wtot','Wkin','Wpot',
     ,       ('flx_',text_an(i),'Eav_',text_an(i),i=1,natmi),
     ,       ('flx_',text_mn(i),'Eav_',text_mn(i),i=1,nmoli),
     ,        'net_ersn','grs_ersn'
          else
            write(lun,hlp_frm)
     -       'zone','DG#','EIR#',
     ,       'R1','Z1','R2','Z2',' X','Xr',
     ,       'area','Wtot','Wkin','Wpot',
     ,       ('flx_',text_an(i),'Eav_',text_an(i),i=1,natmi),
     ,        'net_ersn','grs_ersn'
          end if
        end if
      end if
#ifdef DBG
      print *,'nwldprt,natmi,nmoli',nwldprt,natmi,nmoli
      do il=1,nwldprt
        print '(1x,2i4,1x,a)',il,nwldlis(il),pltzlbl(il)
        write (6,'(1x,20i5)') (lwldlis(i,il),i=1,nwldlis(il))
      end do
#endif
c
      do il=1,nwldprt
        if(nwldlis(il) .gt. 0) then
          if (from_b2plot) write(lun,*)
c*** calculate the length of the target segment
          hlpx=0.0_R8
          do i=1,nwldlis(il)
            ib=iabs(lwldlis(i,il))
            hlpx=hlpx+sqrt((pp1(1,ib)-pp2(1,ib))**2+
     .                     (pp1(2,ib)-pp2(2,ib))**2)
          end do
c
          x=0.0_R8
          distance=0.0_R8
          hlpw=0.0_R8
          do i=1,nwldlis(il)
            ia=lwldlis(i,il)
            ib=iabs(ia)
            if (switch%use_eirene.ne.0) then
              is=leirlis(i,il)
            else
              is=lb25lis(i,il)
            end if
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
            aa=1.0e-6_R8/a
            distance=0.5_R8*distance
            x=x+h+distance
            if (switch%use_eirene.ne.0) then
              if(is.gt.0 .and. is.le.nlimi) then
                hlpw=hlpw+wldnek(is,0)+wldnep(is,0)
                wwk=wldnek(is,0)*aa
                wwp=wldnep(is,0)*aa
                wwt=max(wwk+wwp,pwrctff)
c                wwk=max(wwk,pwrctff)
c                wwp=max(wwp,pwrctff)
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
     .                                -wldpa(is,ic,0))/a
                  do ii=1,nmoli
                    if (mlcmp(ic,ii).gt.0)
     >               wwyn=wwyn+
     .                mlcmp(ic,ii)*(wldrm(is,ii,0)-wldnm(is,ii,0)
     .                                            -wldpm(is,ii,0))/a
                  end do
                else
                  wwyn=flxctff
                  wwyg=flxctff
                end if
              else if (ib.gt.lppp) then
                iend = lppp
                if (grid_mode.lt.2) then
                  j = 0
                  js = 0
                  do ii = 1, ntrgts
                    if (itrg(ii).eq.0) cycle
                    if (j.gt.0.and.js.gt.0) cycle
                    istart = iend+1
                    iend = iend+mpg%divFcP(itrg(ii),2)
                    if (ib.ge.istart.and.ib.le.iend) then
                      j = ii
                      js = ib-istart+1
                    end if
                  end do
                  iFc = mpg%divFc(mpg%divFcP(itrg(j),1)+js-1)
                  ias = b2_fnmti(iFc)
                else
                  ias = ib-lppp
                  do ii = 1, ncns
                    ias = ias+eirdiag_nds_ind(ii+1)-eirdiag_nds_ind(ii)
                  end do
                  ias = ias + 1
                end if
                wwk=0.0_R8
                wwp=0.0_R8
                k=0
                do ii=1,natmi
                  k=k+1
                  if (ias.eq.0) cycle
                  wwk=wwk+(ewlda_res(ii,ias)-ewldea_res(ii,ias))
                  wwj(k)=max(wldna_res(ii,ias),flxctff)/a
                  wwe(k)=max(ewlda_res(ii,ias),enrctff)
                end do
                do ii=1,nmoli
                  k=k+1
                  if (ias.eq.0) cycle
                  wwk=wwk+(ewldm_res(ii,ias)-ewldem_res(ii,ias))
                  wwp=wwp+ewldmr_res(ii,ias)
                  wwj(k)=max(wldnm_res(ii,ias),flxctff)/a
                  wwe(k)=max(ewldm_res(ii,ias),enrctff)
                end do
                if(ic.gt.0.and.ias.gt.0) then
                  wwyg=wldspta_res(ic,ias)/a
                  wwyn=(wldspta_res(ic,ias)-wldna_res(ic,ias)-
     >                                      wlpump_res(ic,ias))/a
                  do ii=1,nmoli
                    if (mlcmp(ic,ii).gt.0)
     >               wwyn=wwyn+mlcmp(ic,ii)*
     .                (wldsptm_res(ii,ias)-wldnm_res(ii,ias)-
     .                                     wlpump(ii+natmi,ias))/a
                  end do
                else
                  wwyn=flxctff
                  wwyg=flxctff
                end if
                hlpw=hlpw+wwk+wwp
                wwk=max(wwk*aa,pwrctff)
                wwp=max(wwp*aa,pwrctff)
                wwt=max(wwk+wwp,pwrctff)
              else
                wwyn=flxctff
                wwyg=flxctff
                wwk=pwrctff
                wwp=pwrctff
                wwt=pwrctff
                k=natmi+nmoli
                do ii=1,k
                  wwj(ii)=0.1_R8*flxctff
                  wwe(ii)=enrctff
                end do
              end if
              if (from_b2plot) then
                write (hlp_frm,'(a,i2,a)')
     ,                 '(1x,a8,2i5,1p,5f10.4,',2*k+12,'e12.4)'
                if (ic.eq.0) then
                  write (lun,hlp_frm) pltzlbl(il),
     ,                ib,is,pp1(1,ib),pp1(2,ib),pp2(1,ib),pp2(2,ib),
     ,                x,x/hlpx,a,wwt,wwk,wwp,(wwj(ii),wwe(ii),ii=1,k)
                else
                  write (lun,hlp_frm) pltzlbl(il),
     ,                ib,is,pp1(1,ib),pp1(2,ib),pp2(1,ib),pp2(2,ib),
     ,                x,x/hlpx,a,wwt,wwk,wwp,(wwj(ii),wwe(ii),ii=1,k),
     ,                wwyn,wwyg
                end if
              end if
            else if (grid_mode.eq.2) then
              u=mpg%bcFcOr(mpg%fcBc(is,1))
              iCv=min(mpg%fcCv(is,1),mpg%fcCv(is,2))
              k=0
              wwk=0.0_R8
              wwp=0.0_R8
              do ii=0,state%ns-1
                if (.not.is_neutral(ii)) cycle
                k=k+1
                hlpw=hlpw+
     .               u*(state%dv%fhm(is,0,ii)+state%dv%fhm(is,1,ii)+
     .                  state%dv%fhp(is,0,ii)+state%dv%fhp(is,1,ii))
                wwk=wwk+
     .               u*(state%dv%fhm(is,0,ii)+state%dv%fhm(is,1,ii))*aa
                wwp=wwp+
     .               u*(state%dv%fhp(is,0,ii)+state%dv%fhp(is,1,ii))*aa
                wwj(k)=max(flxctff,
     .               u*(state%dv%fna(is,0,ii)+state%dv%fna(is,1,ii)))/a
                if (nint(zn(ii)).eq.1) then
                  wwe(k)=max(state%pl%tn(iCv),enrctff)
                else
                  wwe(k)=max(state%pl%ti(iCv),enrctff)
                end if
              enddo
              wwt=max(wwk+wwp,pwrctff)
              wwyn=flxctff
              wwyg=flxctff
              if (from_b2plot) then
                write (hlp_frm,'(a,i2,a)')
     ,                 '(1x,a8,2i5,1p,5f10.4,',2*k+12,'e12.4)'
                write (lun,hlp_frm) pltzlbl(il),
     ,                ib,is,pp1(1,ib),pp1(2,ib),pp2(1,ib),pp2(2,ib),
     ,                x,x/hlpx,a,wwt,wwk,wwp,(wwj(ii),wwe(ii),ii=1,k)
              end if
            end if
          end do
          write(*,*) 'Neutral power to wall: ',pltzlbl(il),il,hlpw
        end if
      end do
      if (from_b2plot) close(lun)
c
c      if(natmi .gt. 1  .and.  npumplis .gt. 0) then
cc*** calculate the helium concentration at the pump
c        a=0.0_R8
c        aa=0.0_R8
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
        if (switch%use_eirene.ne.0) then
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
        else
          do ii=1,mpg%divFcP(itrg(ktns(k)),2)
            iFc=mpg%divFc(mpg%divFcP(itrg(ktns(k)),1)+ii-1)
            u=mpg%divFcOr(mpg%divFcP(itrg(ktns(k)),1)+ii-1)
            iCv=min(mpg%fcCv(iFc,1),mpg%fcCv(iFc,2))
            do i=1,nspecies
              hlp_w(itrg(ktns(k))) = hlp_w(itrg(ktns(k)))+u*
     .         (state%dv%fhm(iFc,0,b2espcr(i))+
     .          state%dv%fhm(iFc,1,b2espcr(i))+
     .          state%dv%fhp(iFc,0,b2espcr(i))+
     .          state%dv%fhp(iFc,1,b2espcr(i)))
              hlp_p(itrg(ktns(k))) = hlp_p(itrg(ktns(k)))+u*
     .         (state%dv%fhp(iFc,0,b2espcr(i))+
     .          state%dv%fhp(iFc,1,b2espcr(i)))
              hlp_a(itrg(ktns(k)),i) = hlp_a(itrg(ktns(k)),i)+u*
     .         (state%dv%fna(iFc,0,b2espcr(i))+
     .          state%dv%fna(iFc,1,b2espcr(i)))
              if (nint(zn(b2espcr(i))).eq.1) then
                hlp_e(itrg(ktns(k)),i) = hlp_a(itrg(ktns(k)),i)+u*
     .          ((state%dv%fna_32(iFc,0,b2espcr(i))+
     .            state%dv%fna_32(iFc,1,b2espcr(i)))*1.5_R8+
     .           (state%dv%fna_52(iFc,0,b2espcr(i))+
     .            state%dv%fna_52(iFc,1,b2espcr(i)))*2.5_R8)*
     .            state%pl%tn(iCv)
              else
                hlp_e(itrg(ktns(k)),i) = hlp_a(itrg(ktns(k)),i)+u*
     .          ((state%dv%fna_32(iFc,0,b2espcr(i))+
     .            state%dv%fna_32(iFc,1,b2espcr(i)))*1.5_R8+
     .           (state%dv%fna_52(iFc,0,b2espcr(i))+
     .            state%dv%fna_52(iFc,1,b2espcr(i)))*2.5_R8)*
     .            state%pl%ti(iCv)
              end if
            end do
          end do
        end if
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
