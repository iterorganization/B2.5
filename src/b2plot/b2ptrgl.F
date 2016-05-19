      subroutine b2ptrgl(ismain,ismain0,BoRiS)
c
c  version : 07.04.98 13:12
c
c======================================================================
c*** Calculate and write down the power loading onto the targets,
c*** related profiles (ne, Te, and so on), and the peak power loads.
c***  Power in W/m**2 (previously, MW/m**2) !!!
c======================================================================
      use b2mod_types
      use b2mod_geo
      use b2mod_diag
      use b2mod_time
      use b2mod_rates
      use b2mod_plasma
      use b2mod_b2plot
      use b2mod_b2cmpa
      use b2mod_indirect
      use b2mod_b2cmpa
      use b2mod_b2cmfs
      use b2mod_constants
      use b2mod_external
      use b2mod_b2plot_wall_loading
#ifdef B25_EIRENE
      use eirmod_extrab25
#endif
      implicit none
#include <DIMENSIONS.F>
      integer ismain, ismain0
      real (kind=R8) :: BoRiS
      integer lun
      parameter (lun=96)
c
      integer i,j,jj,k,l,m,n,iss,iy,js,nclm_ld,nclm_er,nx_hlp
      integer i1,i2,i3,ia,ib,ic,is
      integer jcf,jci,jhei,jhef
      integer m_hlp(ntrgsx)
      real (kind=R8) ::
     *  h,t,u,a,b,d,p,wo,wt,mion,mmol,
     ,  z_hlp(ntrgsx),xtg_hlp(-1:ny-1,ntrgsx),xmp_hlp(0:ny-1,ntrgsx),
     ,  xmc_hlp(0:ny-1,ntrgsx),s(0:ny-1,ntrgsx),wi(ntrgsx),
     ,  zctrg(0:ny-1,ntrgsx),rctrg(0:ny-1,ntrgsx),xctrg(0:ny-1,ntrgsx),
     ,  xtrg(0:ny-1,ntrgsx),netrg(0:ny-1,ntrgsx),
     ,  tetrg(0:ny-1,ntrgsx),wpttrg(0:ny-1,ntrgsx),
     ,  wpltrg(0:ny-1,ntrgsx),wnutrg(0:ny-1,ntrgsx),
     ,  wrdtrg(0:ny-1,ntrgsx,0:DEF_NATM),wtttrg(0:ny-1,ntrgsx),
     ,  wpk_tt(0:ntrgsx),wpk_rd(0:ntrgsx),
     ,  wpk_pl(0:ntrgsx),wpk_nu(0:ntrgsx),
     ,  wlc_tt(ntrgsx),wlc_rd(ntrgsx),wlc_pl(ntrgsx),wlc_nu(ntrgsx),
     ,  pwr_tt(ntrgsx),pwr_rd(ntrgsx),fwht_tt(ntrgsx),fwxt_tt(ntrgsx),
     ,  hwht_tt(ntrgsx),hwhm_tt(ntrgsx),hwhx_tt(ntrgsx),hwxt_tt(ntrgsx),
     ,  hwxm_tt(ntrgsx)
c
c*** Arrays for particle fluxes and erosion
c***  sjca  impinging flux of C atoms
c***  sjcb  "net" sputtering of C (negative means deposition)
c***  sjcg  "gross" erosion rate
c***  sjcn  "net" erosion rate
c***  sjcp  impinging flux of C ions
c***  sjcs  "gross" sputtering of C
c***  sjct  total flux of C particles to the target
c***  sjha  impinging flux of D atoms
c***  sjhp  impinging flux of D ions
c***  sjht  total flux of C particles to the target
c***    suffix "_t" denotes totals from distributions,
c***    and "_w" the totals from Eirene.
c***  sgro  is the mass density of the target material [g/cm**3]
c***  ba    is the atomic mass of the target material (carbon)
c
      real (kind=R8) ::
     *  sjca(0:ny-1,ntrgsx),  sjcb(0:ny-1,ntrgsx), sjcg(0:ny-1,ntrgsx),
     ,  sjcn(0:ny-1,ntrgsx),  sjcp(0:ny-1,ntrgsx), sjcs(0:ny-1,ntrgsx),
     ,  sjct(0:ny-1,ntrgsx),
     * sjhea(0:ny-1,ntrgsx), sjhep(0:ny-1,ntrgsx), sjhet(0:ny-1,ntrgsx),
     * sphea(0:ny-1,ntrgsx),
     *  sjha(0:ny-1,ntrgsx), sjhp(0:ny-1,ntrgsx), sjht(0:ny-1,ntrgsx),
     *  spha(0:ny-1,ntrgsx), thea(0:ny-1,ntrgsx)
      real (kind=R8) ::
     *   sjca_t(ntrgsx),  sjha_t(ntrgsx),  sjcs_t(ntrgsx), 
     ,   sjcp_t(ntrgsx),  sjhp_t(ntrgsx),
     ,  sjhea_t(ntrgsx), sjhep_t(ntrgsx),
     ,   sjca_w(0:ntrgsx),  sjha_w(0:ntrgsx), sjcs_w(0:ntrgsx),
     ,   sjcp_w(0:ntrgsx),  sjhp_w(0:ntrgsx), sjht_w(0:ntrgsx),
     ,   sjct_w(0:ntrgsx),  sjcb_w(0:ntrgsx),
     ,  sjhea_w(0:ntrgsx), sjhep_w(0:ntrgsx), sjhet_w(0:ntrgsx)
      real (kind=R8) ::
     *   sgro, ba, ers_factor, nsecyr
      parameter (ba=12.0_R8, sgro=2.0_R8)
      parameter (nsecyr=24.0_R8*3600.0_R8*365.25_R8)
      parameter (ers_factor=nsecyr*1.0e-3_R8*ba/(avogr*sgro))
      real (kind=R8) ::
     *   ss(ntrgsx),ss_m,ersn_a(ntrgsx),ersn_am,ersn_p(ntrgsx),ersn_pm,
     *   ersn_l(ntrgsx)
      intrinsic sqrt, max
c
c*** Text variables:
c***  ftrg    filenames for the profiles along the targets
c***  fsumm   filenames for the summaries
      integer n_intgrs
      parameter (n_intgrs=8)
      character*12 ftrg(ntrgsx,2),fsumm(4),fn_hlp
      character fmt_peak*80,tab*1,tl*1,
     ,  trg_label(ntrgsx)*9,trg_short(ntrgsx)*2,hlp_trg(ntrgsx)*2,
     ,  pkl_label(2)*9,hdr_2(n_intgrs)*8
      data ftrg /
     ,  'ld_tg_i.dat', 'ld_tg_iu.dat', 'ld_tg_ou.dat', 'ld_tg_o.dat',
     ,  'fp_tg_i.dat', 'fp_tg_iu.dat', 'fp_tg_ou.dat', 'fp_tg_o.dat'/
      data fsumm /
     ,  'peak_power','sput_trg.sum','peak_power_u','peak_power_w'/
      data trg_label /'inner','inner_top','outer_top','outer'/
      data trg_short / 'il'  ,    'iu'   ,    'ou'   , 'ol'  /
      data pkl_label /'     ','- top'/
      data hdr_2 / 'int_rad','int_tot','fwht_tot','fwxt_tot',
     ,  'hwht_tot','hwhm_tot','hwxt_tot','hwxm_tot' /
c======================================================================
c*** At present, only a crude estimate is available for the neutrals:
c*** power flux proportional to n*T**(3/2),
c*** particle flux proportional to n*T**(1/2),
c*** sputtering proportional to the D particle flux (chem. dominant),
c*** re-deposition proportional to the C particle flux.
c*** Here n and T are the density and temperature of the hydrogenic
c*** atoms near the target.
c*** The proportionality factor is adjusted to match the integrals
c*** which are calculated in Eirene.
c
c      print *,'Starting b2ptrgl...'
c      print *,'ers_factor= ',ers_factor
      tab=char(9)

c----------------------------------------------------------------------
c*** Power loading:
c*** Radiation load onto the targets
c
      call b2ptrdl(wrdtrg)
c      print *,'After b2ptrdl: ntrgts,lytg,nytg,lrdlsp = ',   !###
c     ,         ntrgts,lytg,nytg,lrdlsp                       !###
c-----------------------------------------------------------------------
c*** Target geometry

      do j=1,ntrgts !{
        z_hlp(j)=0.0_R8
        k=ltrgpos(j)
        do iy=1,lseppos(j) !{
          if (mod(j,2).eq.1) then ! West targets
            z_hlp(j)=z_hlp(j) - sqrt((crx(k,iy,3)-crx(k,iy,1))**2+
     +                               (cry(k,iy,3)-cry(k,iy,1))**2)
          else ! East targets
            z_hlp(j)=z_hlp(j) - sqrt((crx(k,iy,2)-crx(k,iy,0))**2+
     +                               (cry(k,iy,2)-cry(k,iy,0))**2)
          endif
        end do !}
        xtg_hlp(-1,j)=z_hlp(j)
        do iy=0,ny-1 !{
          if (mod(j,2).eq.1) then ! West targets
            rctrg(iy,j)=0.5_R8*(crx(k,iy,1)+crx(k,iy,3))
            zctrg(iy,j)=0.5_R8*(cry(k,iy,1)+cry(k,iy,3))
            h=sqrt((crx(k,iy,3)-crx(k,iy,1))**2+
     +             (cry(k,iy,3)-cry(k,iy,1))**2)
            s(iy,j)=gs(k,iy,0)
          else ! East targets
            rctrg(iy,j)=0.5_R8*(crx(k,iy,0)+crx(k,iy,2))
            zctrg(iy,j)=0.5_R8*(cry(k,iy,0)+cry(k,iy,2))
            h=sqrt((crx(k,iy,2)-crx(k,iy,0))**2+
     +             (cry(k,iy,2)-cry(k,iy,0))**2)
            s(iy,j)=gs(k,iy+1,0)
          endif
          z_hlp(j)=z_hlp(j)+0.5_R8*h
          xtrg(iy,j)=z_hlp(j)
          z_hlp(j)=z_hlp(j)+0.5_R8*h
          xtg_hlp(iy,j)=z_hlp(j)
        end do !}
      end do !}
c!###{
c      print *
c      print *,'ny,xtg_hlp:',ny
c      do j=1,ntrgts !{
c        print '(1x,i3,10f9.5/(4x,10f9.6))',lseppos(j),
c     ,                                            (xtg_hlp(i,j),i=1,ny)
c      end do !}
c!###}

c*** Mid-plane co-ordinates

      l=jxi
      m=0
      xmc_hlp=0.0_R8
      do j=1,ntrgts !{
        u=0.0_R8
        k=l
        n=2*m
        do i=0,lseppos(j) !{
          h=0.5_R8*sqrt((crx(k,i,n+1)-crx(k,i,n))**2+
     +                  (cry(k,i,n+1)-cry(k,i,n))**2)
          u=u+h
          xmp_hlp(i,j)=u
          u=u+h
        end do !}
        do i=0,lseppos(j) !{
          xmp_hlp(i,j)=xmp_hlp(i,j)-u
        end do !}
        u=0.0_R8
        xmc_hlp(lseppos(j),j)=0.0_R8
        do i=lseppos(j)+1,ny-1 !{
          h=0.5_R8*sqrt((crx(k,i,n+1)-crx(k,i,n))**2+
     +                  (cry(k,i,n+1)-cry(k,i,n))**2)
          u=u+h
          xmp_hlp(i,j)=u
          u=u+h
          xmc_hlp(i,j)=u
        end do !}
        if(j.eq.ntrgts/2) l=jxa
        m=1-m
      end do !}

c!###{
c      print *
c      print *,'ny,xmc_hlp:',ny
c      do j=1,ntrgts !{
c        print '(1x,i3,10f9.5/(4x,10f9.6))',lseppos(j),
c     ,                                            (xmc_hlp(i,j),i=0,ny-1)
c      end do !}
c!###}
c-----------------------------------------------------------------------
c*** Power loading:

c
c*** Plasma and neutrals
c

      l=1
      u=-1.0_R8
      do j=1,ntrgts !{
        k=ltrgpos(j)
        do iy=lytg,ny-1 !{
          wo=fhe(k,iy,0)+fhi(k,iy,0)+fhi_ext(k,iy,0)
          do is=0,ns-1 !{
            wo=wo+fna(k,iy,0,is)*rpt(k,iy,is)*ev
     .           +(1.0_R8-BoRiS)*kinrgy(k,iy,is)
          end do !}
          do is=0,ns_ext-1 !{
            wo=wo+fa_ext(k,iy,0,is)*pt_ext(k,iy,is)*ev
     .           +(1.0_R8-BoRiS)*(0.5_R8*am_ext(is)*mp)*
     .             ua_ext(k,iy,is)**2
          end do !}
c        wo=-1.0e-6_R8*wo/s(iy,1)
          wo= wo/s(iy,j)
          wpltrg(iy,j)=wo*u
          netrg(iy,j)=ne(k+l,iy)
          tetrg(iy,j)=te(k+l,iy)/ev
          wnutrg(iy,j)=0.0_R8
          do is=1,natmi
            t=tab2(k+l+1,iy+1,is,1)/ev
            wnutrg(iy,j)=wnutrg(iy,1)+
     .       s(iy,j)*dab2(k+l+1,iy+1,is,1)*tab2(k+l+1,iy+1,is,1)*
     *       sqrt(max(t/am(eb2spcr(is)),0.0_R8))*0.25_R8*2.0e4_R8
          enddo
          do is=1,nmoli
            mmol=0.0_R8
            do iss = 1,natmi
              mmol=mmol+mlcmp(iss,is)*am(eb2spcr(iss))
            enddo
            t=tmb2(k+l+1,iy+1,is,1)/ev
            wnutrg(iy,j)=wnutrg(iy,j)+
     .       s(iy,j)*dmb2(k+l+1,iy+1,is,1)*tmb2(k+l+1,iy+1,is,1)*
     *       sqrt(max(t/mmol,0.0_R8))*0.25_R8*2.0e4_R8
          enddo
          do is=1,nioni
            mion=0.0_R8
            do iss = 1,natmi
              mion=mion+mlcmp(iss,is)*am(eb2spcr(iss))
            enddo
            t=tib2(k+l+1,iy+1,is,1)/ev
cxpb  The expression below needs to be corrected and to include recombination energy
            wnutrg(iy,j)=wnutrg(iy,j)+
     .       s(iy,j)*dib2(k+l+1,iy+1,is,1)*tib2(k+l+1,iy+1,is,1)*
     *       sqrt(max(t/mion,0.0_R8))*0.25_R8*2.0e4_R8
          enddo
        end do !}
        l=1-l
        u=-u
      end do !}
  
c
c*** Re-scale the neutral target loading and calculate the sums
c
      do j=1,ntrgts !{
        wi(j)=0.0_R8
        do i=lytg,ny-1 !{
          wi(j)=wi(j)+wnutrg(i,j)
        end do !}
      end do !}
c      print *,'ntrgts, wi: ',ntrgts,wi
      do j=1,ntrgts !{
        if(ntns .eq. ntrgts) then !{
          js=ltns(j)+nlim
          wo=wldnek(js,0)+wldnep(js,0)
c          print *,'j,js,wo = ',j,js,wo
          if(wi(j) .gt. 0.0_R8) then !{
            wi(j)=wo/wi(j)
          else !}{
            wi(j)=0.0_R8
          end if !}
        else !}{
          write(*,*) 'b2ptrgl: ntns.ne.ntrgts (',ntns,' vs.',ntrgts,
     ,               ' ) - no neutral load onto the targets calculated'
          write(*,*) '  ==> check specification of ltns in ',
     ,                               'b2.neutrals.parameters'
          wi(j)=0.0_R8
        end if !}
        do iy=lytg,ny-1 !{
c          wnutrg(iy,j)=1.0e-6_R8*wnutrg(iy,j)*wi(j)/s(iy,j)
          wnutrg(iy,j)=wnutrg(iy,j)*wi(j)/s(iy,j)
          wpttrg(iy,j)=wnutrg(iy,j)+wpltrg(iy,j)
          wtttrg(iy,j)=wpttrg(iy,j)+wrdtrg(iy,j,0)
        end do !}
      end do !}
      write(*,'(1x,a,1p,16e12.3)')
     -                   'Scaling factors for Wn : ',(wi(i),i=1,ntrgts)
c
      nclm_ld=11
      if(lrdlsp) nclm_ld=nclm_ld+natmi
      do j=1,ntrgts !{
        jj=j
        if(j.eq.ntrgts) jj=ntrgsx
        open(lun,file=ftrg(jj,1))
        write(fmt_peak,710) nclm_ld,tab
        write(lun,740) 'Target loading - ',trg_label(jj),label
        if(lrdlsp) then !{
          write(lun,fmt_peak) 'x','ne','Te','Wtot','Wpart',
     ,            'Wrad','Wpls','Wneut','xMP','Rclfc','Zclfc',
     ,            ('rad_'//textan(is)(2:3),is=0,natmi-1)
        else !}{
          write(lun,fmt_peak) 'x','ne','Te','Wtot','Wpart',
     ,            'Wrad','Wpls','Wneut','xMP','Rclfc','Zclfc'
        end if !}
        write(fmt_peak,720) nclm_ld,tab
        if(ntrgts.eq.4) then !{
          k=j
        else !}{
          k=ntrgts
        end if !}
        do i=lytg,ny-1 !{
          if(lrdlsp) then !{
            write(lun,fmt_peak) xtrg(i,j),netrg(i,j),tetrg(i,j),
     ,              wtttrg(i,j),wpttrg(i,j),wrdtrg(i,j,0),wpltrg(i,j),
     ,              wnutrg(i,j),xmp_hlp(i,k),rctrg(i,j),zctrg(i,j),
     ,              (wrdtrg(i,j,is),is=1,natmi)
          else !}{
            write(lun,fmt_peak) xtrg(i,j),netrg(i,j),tetrg(i,j),
     ,              wtttrg(i,j),wpttrg(i,j),wrdtrg(i,j,0),wpltrg(i,j),
     ,              wnutrg(i,j),xmp_hlp(i,k),rctrg(i,j),zctrg(i,j)
          end if !}
        end do !}
        close (lun)
      end do !}
c----------------------------------------------------------------------
c*** Calculate the peak power loads onto the targets and write them
c*** down to the summary file
c
      wpk_tt(0)=0
      wpk_rd(0)=0
      wpk_pl(0)=0
      wpk_nu(0)=0
      do j=1,ntrgts !{
        wpk_tt(j)=wtttrg(lytg,j)
        wpk_rd(j)=wrdtrg(lytg,j,0)
        wpk_pl(j)=wpltrg(lytg,j)
        wpk_nu(j)=wnutrg(lytg,j)
        wlc_tt(j)=xtrg(lytg,j)
        wlc_rd(j)=xtrg(lytg,j)
        wlc_pl(j)=xtrg(lytg,j)
        wlc_nu(j)=xtrg(lytg,j)
        m_hlp(j)=lytg
        do i=lytg+1,ny-1 !{
          if(wtttrg(i,j) .gt. wpk_tt(j)) then !{
            wpk_tt(j)=wtttrg(i,j)
            wlc_tt(j)=xtrg(i,j)
            m_hlp(j)=i
          end if !}
          if(wrdtrg(i,j,0) .gt. wpk_rd(j)) then !{
            wpk_rd(j)=wrdtrg(i,j,0)
            wlc_rd(j)=xtrg(i,j)
          end if !}
          if(wpltrg(i,j) .gt. wpk_pl(j)) then !{
            wpk_pl(j)=wpltrg(i,j)
            wlc_pl(j)=xtrg(i,j)
          end if !}
          if(wnutrg(i,j) .gt. wpk_nu(j)) then !{
            wpk_nu(j)=wnutrg(i,j)
            wlc_nu(j)=xtrg(i,j)
          end if !}
        end do !}
      end do !}

      write(fmt_peak,700) 20,tab
  700 format('(1p,e12.5,',i3,'(1h',a1,',e12.4))')
      i1=1
      i2=ntrgts
      fn_hlp=fsumm(1)
      tl=' '
      do j=1,2 !{
        if(j.eq.1 .or. ntrgts.eq.4) then !{
          open(lun,file=fn_hlp)
          write(lun,740) 'Peak powers ',pkl_label(j),label
          i3=i2-i1
          wpk_tt(0)=max(wpk_tt(i1),wpk_tt(i2))
          wpk_rd(0)=max(wpk_rd(i1),wpk_rd(i2))
          wpk_pl(0)=max(wpk_pl(i1),wpk_pl(i2))
          wpk_nu(0)=max(wpk_nu(i1),wpk_nu(i2))
          write(lun,'(a5,20(a1,2a))') '#time',
     -     tab,'peak_t_i',tl,tab,'peak_t_o',tl,tab,'peak_t_x',tl,
     ,     tab, 'loc_t_i',tl,tab, 'loc_t_o',tl,
     ,     tab,'peak_r_i',tl,tab,'peak_r_o',tl,tab,'peak_r_x',tl,
     ,     tab, 'loc_r_i',tl,tab, 'loc_r_o',tl,
     ,     tab,'peak_p_i',tl,tab,'peak_p_o',tl,tab,'peak_p_x',tl,
     ,     tab, 'loc_p_i',tl,tab, 'loc_p_o',tl,
     ,     tab,'peak_n_i',tl,tab,'peak_n_o',tl,tab,'peak_n_x',tl,
     ,     tab, 'loc_n_i',tl,tab, 'loc_n_o',tl
          write(lun,fmt_peak) tim,
     ,     (wpk_tt(i),i=i1,i2,i3), wpk_tt(0), (wlc_tt(i),i=1,2),
     ,     (wpk_rd(i),i=i1,i2,i3), wpk_rd(0), (wlc_rd(i),i=1,2),
     ,     (wpk_pl(i),i=i1,i2,i3), wpk_pl(0), (wlc_pl(i),i=1,2),
     ,     (wpk_nu(i),i=i1,i2,i3), wpk_nu(0), (wlc_nu(i),i=1,2)
          close (lun)
        end if !}
        i1=2
        i2=3
        tl='u'
        fn_hlp=fsumm(3)
      end do !}
c-----------------------------------------------------------------------
c*** Calculate the integral loads onto the targets and power load widths
c*** and write them down to a summary file

      do j=1,ntrgts !{
        hlp_trg(j)=trg_short(j)
        pwr_tt(j)=0.0_R8
        pwr_rd(j)=0.0_R8
        fwht_tt(j)=0.0_R8
        fwxt_tt(j)=0.0_R8
        hwht_tt(j)=0.0_R8
        hwhm_tt(j)=0.0_R8
        hwxt_tt(j)=0.0_R8
        hwxm_tt(j)=0.0_R8
        do i=lytg+1,ny-1 !{
          pwr_tt(j)=pwr_tt(j)+wtttrg(i,j)*s(i,j)
          pwr_rd(j)=pwr_rd(j)+wrdtrg(i,j,0)*s(i,j)
        end do !}
        p=0.5_R8*wpk_tt(j)
        a=xtg_hlp(-1,j)
        if(wtttrg(0,j).lt.p) then !{
          do i=m_hlp(j)-1,0,-1 !{
            if(wtttrg(i,j).lt.p) then !{
              h=(wtttrg(i+1,j)-wtttrg(i,j))
              if(h.gt.0.0_R8) then !{
                h=(p-wtttrg(i,j))/h
              else !}{
                h=1.0_R8
              end if !}
              a=xtg_hlp(i,j)*(1.0_R8-h)+xtg_hlp(i+1,j)*h
              ia=i  !###
              exit
            end if !}
          end do !}
        end if !}
        b=xtg_hlp(ny-1,j)
        d=xmc_hlp(ny-1,j)
        if(wtttrg(ny-1,j).lt.p) then !{
          do i=m_hlp(j)+1,ny-1 !{
            if(wtttrg(i,j).lt.p) then !{
              h=(wtttrg(i-1,j)-wtttrg(i,j))
              if(h.gt.0.0_R8) then !{
                h=(p-wtttrg(i,j))/h
              else !}{
                h=1.0_R8
              end if !}
              b=xtg_hlp(i+1,j)*(1.0_R8-h)+xtg_hlp(i,j)*h
              d=xmc_hlp(i+1,j)*(1.0_R8-h)+xmc_hlp(i,j)*h
              ib=i  !###
              exit
            end if !}
          end do !}
        end if !}
        fwht_tt(j)=b-a
        hwht_tt(j)=b-xtrg(m_hlp(j),j)
        hwhm_tt(j)=d-xmp_hlp(m_hlp(j),j)
c*** Total power outside the strike point position and power half-width
        wt=0.0_R8
        do i=lseppos(j)+1,ny-1 !{
          wt=wt+wtttrg(i,j)*s(i,j)
        end do !}
        fwxt_tt(j)=pwr_tt(j)/(wpk_tt(j)*2.0_R8*pi*rctrg(m_hlp(j),j))
        hwxt_tt(j)=wt/(wpk_tt(j)*2.0_R8*pi*rctrg(m_hlp(j),j))
        hwxm_tt(j)=xmc_hlp(ny-1,j)
        wo=0.0_R8
        do i=lseppos(j)+1,ny-1 !{
          if(xtg_hlp(i,j).ge.hwxt_tt(j)) then !{
            h=(hwxt_tt(j)-h)/(xtg_hlp(i,j)-xtg_hlp(i-1,j))
            hwxm_tt(j)=wo*(1.0_R8-h)+xmc_hlp(i,j)*h
            ic=i  !###
            exit
          end if !}
          h=xtg_hlp(i,j)
          wo=xmc_hlp(i,j)
        end do !}
!###{
        if(fwht_tt(j).lt.0.0_R8 .or. hwht_tt(j).lt.0.0_R8 .or.
     .                               hwhm_tt(j).lt.0.0_R8) then !{
          print *,'!!!! NEGATIVE POWER WIDTH - target ',j,' time ',tim
          print *,'fwht,hwht,hwhm = ',fwht_tt(j),hwht_tt(j),hwhm_tt(j)
          print *,'xtg_hlp:'
          print '(1x,i3,10f10.5/(4x,10f10.5))',lseppos(j)+1,
     ,                                          (xtg_hlp(i,j),i=0,ny-1)
          print *,'xmc_hlp:'
          print '(1x,i3,10f10.5/(4x,10f10.5))',lseppos(j)+1,
     ,                                          (xmc_hlp(i,j),i=0,ny-1)
          print *,'m_hlp,wtttrg:',m_hlp(j)
          print '(1x,i3,1p,10e10.2/(4x,10e10.2))',lseppos(j)+1,
     ,                                          (wtttrg(i,j),i=0,ny-1)
          print '(5a4,a7,2a12)','j','imx','ia','ib','ic','a','b','d'
          print '(5i4,1p,3e12.4)',j,m_hlp(j),ia,ib,ic,a,b,d
        end if !}
!###}
      end do !}
      hlp_trg(ntrgts)=trg_short(ntrgsx)
      fn_hlp=fsumm(4)
      tl=' '
      open(lun,file=fn_hlp)
      write(lun,740) 'Power integrals ','and widths',label
      write(lun,750) 'Divertors are marked "il", "ol" ,"iu", "ou"'
      write(lun,750) hdr_2(1),' : integrated radiation power to the ',
     ,  'target'
      write(lun,750) hdr_2(2),' : integrated total power to the target'
      write(lun,750) hdr_2(3),' : FWHM of total power at the target'
      write(lun,750) hdr_2(4),' : double-sided exponential width of ',
     ,  'total power at the target'
      write(lun,750) hdr_2(5),' : HWHM of total power at the target'
      write(lun,750) hdr_2(6),' : the same mapped to the OMP'
      write(lun,750) hdr_2(7),' : exponential width of the SOL side ',
     ,  'of total power at the target'
      write(lun,750) hdr_2(8),' : the same mapped to the OMP'
      write(fmt_peak,'(a,i3,a)') '(a6,',n_intgrs*ntrgs,'(a1,a))'
      write(lun,fmt_peak) '# time',
     -  ((tab,trim(hdr_2(i))//'_'//hlp_trg(j),i=1,n_intgrs),j=1,ntrgs)
      write(fmt_peak,705) n_intgrs*ntrgs,tab
  705 format('(1p,e12.5,',i3,'(1h',a1,',e12.4))')
      write(lun,fmt_peak) tim,
     ,        (pwr_rd(j),pwr_tt(j),fwht_tt(j),fwxt_tt(j),
     ,         hwht_tt(j),hwhm_tt(j),hwxt_tt(j),hwxm_tt(j),j=1,ntrgs)
      close (lun)
c----------------------------------------------------------------------
c*** Calculate the particle fluxes and erosion - if present
c*** At present, only for carbon - so determine where it is
c
      jcf=0
      jci=-1
      jhef=0
      jhei=-1
      if(nspecies .gt. 1) then
        do is=1,nspecies
          if (is_neutral(jci+1)) jci=jci+1
          if (is_neutral(jhei+1)) jhei=jhei+1
          if (nint(zn(jci+1)) .eq. 6) then
            jcf=is
          else
            jci=jci+nfluids(is)
          end if
          if (nint(zn(jhei+1)) .eq. 2) then
            jhef=is
          else
            jhei=jhei+nfluids(is)
          end if
        end do
      end if
c      print *,'jcf,jci = ',jcf,jci
c      print *,'jhef,jhei = ',jhef,jhei
c
c*** Now jcf > 0 shows the isonuclear sequence index for C and jhef > 0 for He
c*** - assuming the same particle order in B2 and Eirene -
c*** and jci and jhei mark the first ion species for C and He

c*** Calculate the particle fluxes

c      print *,'Particle fluxes...' !###

      if(ntns.eq.ntrgts) then !{

        do j=1,ntrgts !{
          do i=0,ny-1 !{
            sjca(i,j)=0.0_R8
            sjcb(i,j)=0.0_R8
            sjcg(i,j)=0.0_R8
            sjcn(i,j)=0.0_R8
            sjcp(i,j)=0.0_R8
            sjcs(i,j)=0.0_R8
            sjct(i,j)=0.0_R8
            sjha(i,j)=0.0_R8
            sjhp(i,j)=0.0_R8
            sjht(i,j)=0.0_R8
            sjhea(i,j)=0.0_R8
            sjhep(i,j)=0.0_R8
            sjhet(i,j)=0.0_R8
            sphea(i,j)=0.0_R8
            spha(i,j)=0.0_R8
            thea(i,j)=0.0_R8
          end do !}
          sjca_w(j)=0.0_R8
          sjcb_w(j)=0.0_R8
          sjcp_w(j)=0.0_R8
          sjcs_w(j)=0.0_R8
          sjct_w(j)=0.0_R8
          sjha_w(j)=0.0_R8
          sjhp_w(j)=0.0_R8
          sjht_w(j)=0.0_R8
          sjhea_w(j)=0.0_R8
          sjhep_w(j)=0.0_R8
          sjhet_w(j)=0.0_R8
          ersn_a(j)=0.0_R8
          ersn_p(j)=0.0_R8
          ersn_l(j)=0.0_R8
        end do !}

        u=-1.0_R8
        l=1
        do j=1,ntrgts !{
          k=ltrgpos(j)
c          print *,'j,ltrgpos = ',j,k  !###
          sjca_t(j)=0.0_R8
          sjcs_t(j)=0.0_R8
          sjha_t(j)=0.0_R8
          sjhp_t(j)=0.0_R8
          sjhea_t(j)=0.0_R8
          sjhep_t(j)=0.0_R8
          sjcp_t(j)=0.0_R8
          do iy=lytg,ny-1 !{
            sjhp(iy,j)=u*fna(k,iy,0,ismain)
            t=tab2(k+l+1,iy+1,b2espcr(ismain0),1)/ev
            sjha(iy,j)=s(iy,2)*dab2(k+l+1,iy+1,b2espcr(ismain0),1)*
     *                 sqrt(max(t/am(ismain),0.0_R8))*0.25_R8*1.0e4_R8
            spha(iy,j)=dab2(k+l+1,iy+1,b2espcr(ismain0),1)*
     *                 tab2(k+l+1,iy+1,b2espcr(ismain0),1)
            sjhp_t(j)=sjhp_t(j)+sjhp(iy,j)
            sjha_t(j)=sjha_t(j)+sjha(iy,j)
            sjcp(iy,j)=0.0_R8
            sjca(iy,j)=0.0_R8
            if(jcf.gt.0) then !{
              do is=1,6 !{
                sjcp(iy,j)=sjcp(iy,j)+fna(k,iy,0,jci-1+is)
              end do !}
              sjcp(iy,j)=u*sjcp(iy,j)
              t=tab2(k+l+1,iy+1,jcf,1)/ev
              sjca(iy,j)=s(iy,j)*dab2(k+l+1,iy+1,jcf,1)*
     *                   sqrt(max(t/am(jci),0.0_R8))*0.25_R8*1.0e4_R8
              sjcp_t(j)=sjcp_t(j)+sjcp(iy,j)
              sjca_t(j)=sjca_t(j)+sjca(iy,j)
            end if !}
            if(jhef.gt.0) then !{
              do is=1,2 !{
                sjhep(iy,j)=sjhep(iy,j)+fna(k,iy,0,jhei-1+is)
              end do !}
              sjhep(iy,j)=u*sjhep(iy,j)
              t=tab2(k+l+1,iy+1,jhef,1)/ev
              sjhea(iy,j)=s(iy,j)*dab2(k+l+1,iy+1,jhef,1)*
     *                    sqrt(max(t/am(jhei),0.0_R8))*0.25_R8*1.0e4_R8
              sjhep_t(j)=sjhep_t(j)+sjhep(iy,j)
              sjhea_t(j)=sjhea_t(j)+sjhea(iy,j)
              sphea(iy,j)=dab2(k+l+1,iy+1,jhef,1)*
     *                    tab2(k+l+1,iy+1,jhef,1)
              thea(iy,j)=tab2(k+l+1,iy+1,jhef,1)/ev
            end if !}
          end do !}
          l=1-l
          u=-u
        end do !}
c!###{
c        print '(6(3x,a7,2x)/1p(6e12.4))','sjha_t','sjhea_t',
c     ,    'sjca_t','sjhp_t','sjhep_t','sjcp_t',
c     ,    (sjha_t(i),sjhea_t(i),sjca_t(i),sjhp_t(i),
c     ,                          sjhep_t(i),sjcp_t(i),i=1,ntrgts)
c          print '(2a5,6(3x,a7,2x))','j','js',                        !###
c     ,      'sjha_w','sjhea_w','sjca_w','sjhp_w','sjhep_w','sjcp_w'
c!###}
c
        do j=1,ntrgts !{
          js=ltns(j)+nlim
          sjha_w(j)=wldna(js,1,0)
          sjhp_w(j)=wldpp(js,1,0)
          sjht_w(j)=sjha_w(j)+sjhp_w(j)
          if(jcf.gt.0) then !{
            sjca_w(j)=wldna(js,jcf,0)
            sjcp_w(j)=0
            do is=1,6 !{
              sjcp_w(j)=sjcp_w(j)+wldpp(js,jci-1+is,0)
            end do !}
            sjct_w(j)=sjca_w(j)+sjcp_w(j)
          end if !}
          if(jhef.gt.0) then !{
            sjhea_w(j)=wldna(js,jhef,0)
            sjhep_w(j)=0
            do is=1,2 !{
              sjhep_w(j)=sjhep_w(j)+wldpp(js,jhei-1+is,0)
            end do !}
            sjhet_w(j)=sjhea_w(j)+sjhep_w(j)
          end if !}
c!###{
c          print '(2i5,1p,6e12.4)',j,js,sjha_w(j),sjhea_w(j),
c     ,      sjca_w(j),sjhp_w(j),sjhep_w(j),sjcp_w(j)
c!###}
          if(sjha_t(j).gt.0.0_R8) then !{
            sjha_t(j)=sjha_w(j)/sjha_t(j)
          else !}{
            sjha_t(j)=0.0_R8
          end if !}
          if(sjhp_t(j).gt.0.0_R8) then !{
            sjhp_t(j)=sjhp_w(j)/sjhp_t(j)
          else !}{
            sjhp_t(j)=0.0_R8
          end if !}
          if(jcf.gt.0) then !{
            if(sjca_t(j).gt.0.0_R8) then !{
              sjca_t(j)=sjca_w(j)/sjca_t(j)
            else !}{
              sjca_t(j)=0.0_R8
            end if !}
            if(sjcp_t(j).gt.0.0_R8) then !{
              sjcp_t(j)=sjcp_w(j)/sjcp_t(j)
            else !}{
              sjcp_t(j)=0.0_R8
            end if !}
          end if !}
          if(jhef.gt.0) then !{
            if(sjhea_t(j).gt.0.0_R8) then !{
              sjhea_t(j)=sjhea_w(j)/sjhea_t(j)
            else !}{
              sjhea_t(j)=0.
            end if !}
            if(sjhep_t(j).gt.0.0_R8) then !{
              sjhep_t(j)=sjhep_w(j)/sjhep_t(j)
            else !}{
              sjhep_t(j)=0.0_R8
            end if !}
          end if !}
          do iy=lytg,ny-1 !{
            sjha(iy,j)=sjha(iy,j)*sjha_t(j)
            sjhp(iy,j)=sjhp(iy,j)*sjhp_t(j)
            sjht(iy,j)=sjha(iy,j)+sjhp(iy,j)
            sjca(iy,j)=sjca(iy,j)*sjca_t(j)
            sjcp(iy,j)=sjcp(iy,j)*sjcp_t(j)
            sjct(iy,j)=sjca(iy,j)+sjcp(iy,j)
            sjhea(iy,j)=sjhea(iy,j)*sjhea_t(j)
            sjhep(iy,j)=sjhep(iy,j)*sjhep_t(j)
            sjhet(iy,j)=sjhea(iy,j)+sjhep(iy,j)
          end do !}
        end do !}
        write(*,'(1x,a,1p,16e12.3)') 'Scaling factors for jD : ',
     ,                                           (sjha_t(i),i=1,ntrgts)
        write(*,'(1x,a,1p,16e12.3)') 'Scaling factors for jD+: ',
     ,                                           (sjhp_t(i),i=1,ntrgts)
        write(*,'(1x,a,1p,16e12.3)') 'Scaling factors for jC : ',
     ,                                           (sjca_t(i),i=1,ntrgts)
        write(*,'(1x,a,1p,16e12.3)') 'Scaling factors for jC+: ',
     ,                                           (sjcp_t(i),i=1,ntrgts)
        write(*,'(1x,a,1p,16e12.3)') 'Scaling factors for jHe : ',
     ,                                           (sjhea_t(i),i=1,ntrgts)
        write(*,'(1x,a,1p,16e12.3)') 'Scaling factors for jHe+: ',
     ,                                           (sjhep_t(i),i=1,ntrgts)

c*** Calculate sputtering

c      print *,'Sputtering...'
        if(jcf .gt. 0) then !{
          do j=1,ntrgts !{
            do iy=lytg,ny-1 !{
              sjcs(iy,j)=sjht(iy,j)
              sjcs_t(j)=sjcs_t(j)+sjcs(iy,j)
            end do !}
            if(ntns.eq.2 .or. ntns.eq.4) then !{
              js=ltns(j)+nlim
              sjcs_w(j)=wldspt(js,0)
              sjcb_w(j)=sjcs_w(j)-sjct_w(j)
c              print *,'j,js,sjcs_w,sjcb_w = ',        !###
c     ,                      j,js,sjcs_w(j),sjcb_w(j)  !###
              if(sjcs_t(j).gt.0.0_R8) then !{
                sjcs_t(j)=sjcs_w(j)/sjcs_t(j)
              else !}{
                sjcs_t(j)=0.0_R8
              end if !}
            else !}{
              sjcs_t(j)=0.0_R8
            end if !}
            do iy=lytg,ny !{
              sjcs(iy,j)=sjcs_t(j)*sjcs(iy,j)
              sjcb(iy,j)=(sjcs(iy,j)-sjct(iy,j))/s(iy,j)
              sjcn(iy,j)=ers_factor*sjcb(iy,j)
              sjcs(iy,j)=sjcs(iy,j)/s(iy,j)
              sjcg(iy,j)=ers_factor*sjcs(iy,j)
              sjca(iy,j)=sjca(iy,j)/s(iy,j)
              sjha(iy,j)=sjha(iy,j)/s(iy,j)
              sjcp(iy,j)=sjcp(iy,j)/s(iy,j)
              sjhp(iy,j)=sjhp(iy,j)/s(iy,j)
              sjct(iy,j)=sjct(iy,j)/s(iy,j)
              sjht(iy,j)=sjht(iy,j)/s(iy,j)
              sjhea(iy,j)=sjhea(iy,j)/s(iy,j)
              sjhep(iy,j)=sjhep(iy,j)/s(iy,j)
              sjhet(iy,j)=sjhet(iy,j)/s(iy,j)
            end do !}
          end do !}
          write(*,'(1x,a,1p,16e12.3)') 'Scaling factors for Sp : ',
     ,                                           (sjcs_t(i),i=1,ntrgts)
        end if !}

        nclm_er=16
        do j=1,ntrgts !{
          nx_hlp=ltrgpos(j)
          jj=j
          if(j.eq.ntrgts) jj=ntrgsx
          open(lun,file=ftrg(jj,2))
          write(fmt_peak,710) nclm_er,tab
          write(lun,740) 'Target erosion - ',trg_label(jj),label
          write(lun,fmt_peak) 'x','Sput_net','Sput_gross',
     ,            'Ersn_net','Ersn_gross','D_flux_atm','D_flux_ion',
     ,            'C_flux_atm','C_flux_ion','T_atm_D','He_flx_atm',
     ,            'He_flx_ion','T_atm_He','pr_atm_He','pr_atm_D','xMP'
          write(fmt_peak,720) nclm_er,tab
          if(ntrgts.eq.4) then !{
            k=j
          else !}{
            k=ntrgts
          end if !}
          do i=lytg,ny-1 !{
            write(lun,fmt_peak) xtrg(i,j),sjcb(i,j),sjcs(i,j),
     ,            sjcn(i,j),sjcg(i,j),sjha(i,j),sjhp(i,j),
     ,            sjca(i,j),sjcp(i,j),
     ,            tab2(nx_hlp+1,i+1,b2espcr(ismain0),1)/ev,
     ,            sjhea(i,j),sjhep(i,j),thea(i,j),sphea(i,j),spha(i,j),
     ,            xmp_hlp(i,k)
          end do !}
          close (lun)
        end do !}
c----------------------------------------------------------------------
c*** Produce a summary file for target sputtering
c

        sjhp_w(0)=0.0_R8
        sjha_w(0)=0.0_R8
        sjht_w(0)=0.0_R8
        sjcp_w(0)=0.0_R8
        sjca_w(0)=0.0_R8
        sjct_w(0)=0.0_R8
        sjcs_w(0)=0.0_R8
        sjcb_w(0)=0.0_R8
        do j=1,ntrgts !{
          sjhp_w(0)=sjhp_w(0)+sjhp_w(j)
          sjha_w(0)=sjha_w(0)+sjha_w(j)
          sjht_w(0)=sjht_w(0)+sjht_w(j)
          sjcp_w(0)=sjcp_w(0)+sjcp_w(j)
          sjca_w(0)=sjca_w(0)+sjca_w(j)
          sjct_w(0)=sjct_w(0)+sjct_w(j)
          sjcs_w(0)=sjcs_w(0)+sjcs_w(j)
          sjcb_w(0)=sjcb_w(0)+sjcb_w(j)
        end do !}
        k=1
        if(ntrgts.eq.2) k=ntrgsx-1
        open(lun,file=fsumm(2))
        write(lun,*) 'Target erosion summary: ',label
        write(lun,*)
        write(lun,'(12x,3(a1,a12))')
     -                   (tab,trg_label(i),i=1,ntrgsx,k),tab,'total   '
        write(fmt_peak,730) ntrgts+1,tab
  730   format('(1p,1x,a11,',i3,'(1h',a1,',e12.4))')
        write(lun,*) 'Particle fluxes [1/s]'
        write(lun,fmt_peak) 'D ions on',
     ,                               (sjhp_w(i),i=1,ntrgts),sjhp_w(0)
        write(lun,fmt_peak) 'D atoms on',
     ,                               (sjha_w(i),i=1,ntrgts),sjha_w(0)
        write(lun,fmt_peak) 'D total on',
     ,                               (sjht_w(i),i=1,ntrgts),sjht_w(0)
        write(lun,fmt_peak) 'He ions on',
     ,                               (sjhep_w(i),i=1,ntrgts),sjhep_w(0)
        write(lun,fmt_peak) 'He atoms on',
     ,                               (sjhea_w(i),i=1,ntrgts),sjhea_w(0)
        write(lun,fmt_peak) 'He total on',
     ,                               (sjhet_w(i),i=1,ntrgts),sjhet_w(0)
        write(lun,fmt_peak) 'C ions on',
     ,                               (sjcp_w(i),i=1,ntrgts),sjcp_w(0)
        write(lun,fmt_peak) 'C atoms on',
     ,                               (sjca_w(i),i=1,ntrgts),sjca_w(0)
        write(lun,fmt_peak) 'C total on',
     ,                               (sjct_w(i),i=1,ntrgts),sjct_w(0)
        write(lun,fmt_peak) 'C sputtered',
     ,                               (sjcs_w(i),i=1,ntrgts),sjcs_w(0)
        write(lun,fmt_peak) 'C removed',
     ,                               (sjcb_w(i),i=1,ntrgts),sjcb_w(0)
        write(lun,*)
        write(lun,*) 'Erosion [mm/y]'
        ersn_pm=0.0_R8
        ersn_am=0.0_R8
        ss_m=0.0_R8
        do j=1,ntrgts !{
          ss(j)=0.0_R8
          ersn_p(j)=0.0_R8
          do iy=lytg,ny-1 !{
            ss(j)=ss(j)+s(iy,j)
            if(sjcn(iy,j) .gt. ersn_p(j)) then !{
              ersn_p(j)=sjcn(iy,j)
              ersn_l(j)=xtrg(iy,j)
            end if !}
          end do !}
          ersn_a(j)=sjcb_w(j)*ers_factor/ss(j)
          ersn_pm=max(ersn_p(j),ersn_pm)
          ss_m=ss_m+ss(j)
          ersn_am=ersn_am+sjcb_w(j)
        end do !}
        ersn_am=ersn_am*ers_factor/ss_m
        write(lun,fmt_peak) 'average',(ersn_a(i),i=1,ntrgts),ersn_am
        write(lun,fmt_peak) 'peak',(ersn_p(i),i=1,ntrgts),ersn_pm
        write(lun,fmt_peak) 'peak pos.',(ersn_l(i),i=1,ntrgts)

        close (lun)
      else !}{
        write(*,*) 'b2ptrgl: ntns.ne.ntrgts (',ntns,' vs.',ntrgts,
     ,                                 ' ) - no sputtering calculated.'
        write(*,*) '  ==> check specification of ltns in ',
     ,                               'b2.parameters (section NEUTRALS)'
        do j=1,ntrgts !{
          sjha_t(j)=0.0_R8
          sjca_t(j)=0.0_R8
          do iy=0,ny-1 !{
            sjha(iy,j)=0.0_R8
            sjhp(iy,j)=0.0_R8
            sjht(iy,j)=0.0_R8
            sjca(iy,j)=0.0_R8
            sjcp(iy,j)=0.0_R8
            sjct(iy,j)=0.0_R8
          end do !}
        end do !}
      end if !}
      return
c======================================================================
  710 format('(2h# ,',i3,'(1x,a11,1h',a1,'))')
  720 format('(1p,',i3,'(e12.4,1h',a1,'))')
  740 format('# ',2a,' : ',a)
  750 format('# ',4a)
c=======================================================================
      end subroutine b2ptrgl
