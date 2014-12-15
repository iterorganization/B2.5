      subroutine b2ptrgl(BoRiS)
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
      use b2mod_indirect
      use b2mod_b2cmpa
      use b2mod_constants
      use eirmod_extrab25
      implicit none
      real (kind=R8) :: BoRiS
      integer lun
      parameter (lun=96)
c
      integer i,j,is,iy,js,jcf,jci,nclm_ld,nclm_er,nx_hlp
      real (kind=R8) ::
     *  h,t,u,wo,zwest,zeast,s(0:ny-1,2),wi(2),
     ,  xtrg(0:ny-1,2),netrg(0:ny-1,2),tetrg(0:ny-1,2),wpttrg(0:ny-1,2),
     ,  wpltrg(0:ny-1,2),wnutrg(0:ny-1,2),
     ,  wrdtrg(0:ny-1,2),wtttrg(0:ny-1,2),
     ,  wpk_tt(3),wpk_rd(3),wpk_pl(3),wpk_nu(3),
     ,  wlc_tt(2),wlc_rd(2),wlc_pl(2),wlc_nu(2)
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
     *   sjca(ny,2), sjcb(ny,2), sjcg(ny,2), sjcn(ny,2),
     ,   sjcp(ny,2), sjcs(ny,2), sjct(ny,2),
     ,   sjha(ny,2), sjhp(ny,2), sjht(ny,2)
      real (kind=R8) ::
     *   sjca_t(2), sjha_t(2), sjcs_t(2), sjcp_t(2), sjhp_t(2),
     ,   sjca_w(2), sjha_w(2), sjcs_w(2), sjcp_w(2),
     ,   sjhp_w(2), sjht_w(2), sjct_w(2), sjcb_w(2)
      real (kind=R8) ::
     *   sgro, ba, ers_factor, nsecyr
      parameter (ba=12.0_R8, sgro=2.0_R8)
      parameter (nsecyr=24.0_R8*3600.0_R8*365.25_R8)
      parameter (ers_factor=nsecyr*1.0e-3_R8*ba/(avogr*sgro))
      real (kind=R8) ::
     *   ss(2),ss_m,ersn_a(2),ersn_am,ersn_p(2),ersn_pm,ersn_l(2)
      intrinsic sqrt, max
c
c*** Text variables:
c***  ftrg    filenames for the profiles along the targets
c***  fsumm   filenames for the summaries
      character*12 ftrg(2,2),fsumm(2),fmt_peak*80,tab*1,trg_label(2)*5
      data ftrg /'ld_tg_i.dat', 'ld_tg_o.dat',
     ,           'fp_tg_i.dat', 'fp_tg_o.dat'/
      data fsumm /'peak_power','sput_trg.sum'/
      data trg_label /'inner','outer'/
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
c
c*** Target geometry
c
      zwest=0.0_R8
      zeast=0.0_R8
      do iy=0,topcut(1)
        zwest=zwest - sqrt((crx(0,iy,2)-crx(0,iy,0))**2+
     +                     (cry(0,iy,2)-cry(0,iy,0))**2)
        zeast=zeast - sqrt((crx(nx-1,iy,3)-crx(nx-1,iy,1))**2+
     +                     (cry(nx-1,iy,3)-cry(nx-1,iy,1))**2)
      end do
      do iy=0,ny-1
        u=0.5_R8*(crx(0,iy,2)+crx(0,iy,0))*2.0_R8*pi
        h=sqrt((crx(0,iy,2)-crx(0,iy,0))**2+
     +         (cry(0,iy,2)-cry(0,iy,0))**2)
        s(iy,1)=u*h
        zwest=zwest+0.5_R8*h
        xtrg(iy,1)=zwest
        u=0.5_R8*(crx(nx-1,iy,3)+crx(nx-1,iy,1))*2.0_R8*pi
        h=sqrt((crx(nx-1,iy,3)-crx(nx-1,iy,1))**2+
     +         (cry(nx-1,iy,3)-cry(nx-1,iy,1))**2)
        s(iy,2)=u*h
        zeast=zeast+0.5_R8*h
        xtrg(iy,2)=zeast
        zwest=zwest+0.5_R8*h
        zeast=zeast+0.5_R8*h
      end do
c----------------------------------------------------------------------
c*** Power loading:
c*** Radiation load onto the targets
c
      call b2ptrdl(wrdtrg)
c
c*** Plasma and neutrals
c
      do iy=0,ny-1
        wo=fhe(0,iy,0)+fhi(0,iy,0)
        do is=0,ns-1 !{
          wo=wo+fna(0,iy,0,is)*(
     .     (rpt(0,iy,is)*hx(-1,iy)+rpt(-1,iy,is)*hx(0,iy))/
     .     (hx(-1,iy)+hx(0,iy))*ev+(1.0_R8-BoRiS)*
     .     (kinrgy(0,iy,is)*hx(-1,iy)+kinrgy(-1,iy,is)*hx(0,iy))/
     .     (hx(-1,iy)+hx(0,iy)))
        end do !}
c        wo=-1.0e-6_R8*wo/s(iy,1)
        wo=-wo/s(iy,1)
        wpltrg(iy,1)=wo
        netrg(iy,1)=ne(0,iy)
        tetrg(iy,1)=te(0,iy)/ev
        t=tab2(1,iy+1,1,1)/ev
        wnutrg(iy,1)=s(iy,1)*dab2(1,iy+1,1,1)*tab2(1,iy+1,1,1)*
     *   sqrt(max(t/am(1),0.0_R8))*0.25_R8*2.0e4_R8
c
        wo= fhe(nx,iy,0)+fhi(nx,iy,0)
        do is=0,ns-1 !{
          wo=wo+fna(nx,iy,0,is)*(
     .     (rpt(nx,iy,is)*hx(nx-1,iy)+rpt(nx-1,iy,is)*hx(nx,iy))/
     .     (hx(nx-1,iy)+hx(nx,iy))*ev+(1.0_R8-BoRiS)*
     .     (kinrgy(nx,iy,is)*hx(nx-1,iy)+kinrgy(nx-1,iy,is)*hx(nx,iy))/
     .     (hx(nx-1,iy)+hx(nx,iy)))
        end do !}
c        wo= 1.0e-6_R8*wo/s(iy,2)
        wo= wo/s(iy,2)
        wpltrg(iy,2)=wo
        netrg(iy,2)=ne(nx-1,iy)
        tetrg(iy,2)=te(nx-1,iy)/ev
        t=tab2(nx,iy+1,1,1)/ev
        wnutrg(iy,2)=s(iy,2)*dab2(nx,iy+1,1,1)*tab2(nx,iy+1,1,1)*
     *   sqrt(max(t/am(1),0.0_R8))*0.25_R8*2.0e4_R8
      end do
c
c*** Re-scale the neutral target loading and calculate the sums
c
      do j=1,2
        wi(j)=0.0_R8
        do i=0,ny-1
          wi(j)=wi(j)+wnutrg(i,j)
        end do
      end do
c      print *,'ntns, wi: ',ntns,wi
      do j=1,2
        if(ntns .eq. 2) then
          js=ltns(j)+nlim
          wo=wldnek(js,0)+wldnep(js,0)
c          print *,'j,js,wo = ',j,js,wo
          if(wi(j) .gt. 0.0_R8) then
            wi(j)=wo/wi(j)
          else
            wi(j)=0.0_R8
          end if
        else
          wi(j)=0.0_R8
        end if
        do iy=0,ny-1
c          wnutrg(iy,j)=1.0e-6_R8*wnutrg(iy,j)*wi(j)/s(iy,j)
          wnutrg(iy,j)=wnutrg(iy,j)*wi(j)/s(iy,j)
          wpttrg(iy,j)=wnutrg(iy,j)+wpltrg(iy,j)
          wtttrg(iy,j)=wpttrg(iy,j)+wrdtrg(iy,j)
        end do
      end do
      print *,'Scaling factors for Wn: ',wi
c
      nclm_ld=8
      do j=1,2
        open(lun,file=ftrg(j,1))
        write(fmt_peak,710) nclm_ld,tab
        write(lun,740) 'Target loading - ',trg_label(j),label
        write(lun,fmt_peak) 'x','ne','Te',
     ,                             'Wtot','Wpart','Wrad','Wpls','Wneut'
        write(fmt_peak,720) nclm_ld,tab
        do i=0,ny-1
          write(lun,fmt_peak) xtrg(i,j),netrg(i,j),tetrg(i,j),
     ,      wtttrg(i,j),wpttrg(i,j),wrdtrg(i,j),wpltrg(i,j),wnutrg(i,j)
        end do
        close (lun)
      end do
c----------------------------------------------------------------------
c*** Calculate the peak power loads onto the targets and write them
c*** down to the summary file
c
      do j=1,2
        wpk_tt(j)=wtttrg(0,j)
        wpk_rd(j)=wrdtrg(0,j)
        wpk_pl(j)=wpltrg(0,j)
        wpk_nu(j)=wnutrg(0,j)
        wlc_tt(j)=xtrg(0,j)
        wlc_rd(j)=xtrg(0,j)
        wlc_pl(j)=xtrg(0,j)
        wlc_nu(j)=xtrg(0,j)
        do i=1,ny-1
          if(wtttrg(i,j) .gt. wpk_tt(j)) then
            wpk_tt(j)=wtttrg(i,j)
            wlc_tt(j)=xtrg(i,j)
          end if
          if(wrdtrg(i,j) .gt. wpk_rd(j)) then
            wpk_rd(j)=wrdtrg(i,j)
            wlc_rd(j)=xtrg(i,j)
          end if
          if(wpltrg(i,j) .gt. wpk_pl(j)) then
            wpk_pl(j)=wpltrg(i,j)
            wlc_pl(j)=xtrg(i,j)
          end if
          if(wnutrg(i,j) .gt. wpk_nu(j)) then
            wpk_nu(j)=wnutrg(i,j)
            wlc_nu(j)=xtrg(i,j)
          end if
        end do
      end do
      wpk_tt(3)=max(wpk_tt(1),wpk_tt(2))
      wpk_rd(3)=max(wpk_rd(1),wpk_rd(2))
      wpk_pl(3)=max(wpk_pl(1),wpk_pl(2))
      wpk_nu(3)=max(wpk_nu(1),wpk_nu(2))
c
      write(fmt_peak,700) 20,tab
  700 format('(1p,e12.5,',i3,'(1h',a1,',e12.4))')
      open(lun,file=fsumm(1))
      write(lun,'(a5,20(a1,a))') '#time',
     -  tab,'peak_t_i',tab,'peak_t_o',tab,'peak_t_x',
     ,  tab,'loc_t_i',tab,'loc_t_o',
     ,  tab,'peak_r_i',tab,'peak_r_o',tab,'peak_r_x',
     ,  tab,'loc_r_i',tab,'loc_r_o',
     ,  tab,'peak_p_i',tab,'peak_p_o',tab,'peak_p_x',
     ,  tab,'loc_p_i',tab,'loc_p_o',
     ,  tab,'peak_n_i',tab,'peak_n_o',tab,'peak_n_x',
     ,  tab,'loc_n_i',tab,'loc_n_o'
      write(lun,fmt_peak) tim,
     ,  (wpk_tt(i),i=1,3), (wlc_tt(i),i=1,2),
     ,  (wpk_rd(i),i=1,3), (wlc_rd(i),i=1,2),
     ,  (wpk_pl(i),i=1,3), (wlc_pl(i),i=1,2),
     ,  (wpk_nu(i),i=1,3), (wlc_nu(i),i=1,2)
c
        close (lun)
c----------------------------------------------------------------------
c*** Calculate the particle fluxes and erosion - if present
c*** At present, only for carbon - so determine where it is
c
      jcf=0
      if(nspecies .gt. 1) then
        do is=0,nspecies-1
          if(nfluids(is) .eq. 6) jcf=is
        end do
      end if
      jci=0
      do j=2,jcf
        jci=jci+nfluids(j-1)
      end do
c      print *,'jcf,jci = ',jcf,jci
c
c*** Now jcf > 0 shows the particle index for C
c*** - assuming the same particle order in B2 and Eirene -
c*** and jci is the offset for the carbon ions in the whole ion set
c
      do j=1,2
        sjca_t(j)=0.0_R8
        sjcs_t(j)=0.0_R8
        sjha_t(j)=0.0_R8
        sjhp_t(j)=0.0_R8
        sjcp_t(j)=0.0_R8
      end do
c
c*** Calculate the particle fluxes
c
c      print *,'Particle fluxes...'
      do iy=0,ny-1
        sjhp(iy,1)=-fna(0,iy,1,0)
        sjhp(iy,2)= fna(nx,iy,1,0)
        sjcp(iy,1)=0
        sjcp(iy,2)=0
        do is=1,6
          sjcp(iy,1)=sjcp(iy,1)-fna(0,iy,jci+is,0)
          sjcp(iy,2)=sjcp(iy,2)+fna(nx,iy,jci+is,0)
        end do
        t=tab2(1,iy+1,1,1)/ev
        sjha(iy,1)=s(iy,1)*dab2(1,iy+1,1,1)*
     *   sqrt(max(t/am(1),0.0_R8))*0.25_R8*1.0e4_R8
        t=tab2(nx,iy+1,1,1)/ev
        sjha(iy,2)=s(iy,2)*dab2(nx,iy+1,1,1)*
     *   sqrt(max(t/am(1),0.0_R8))*0.25_R8*1.0e4_R8
        if(jcf .gt. 0) then
          t=tab2(1,iy+1,jcf,1)/ev
          sjca(iy,1)=s(iy,1)*dab2(1,iy+1,jcf,1)*
     *     sqrt(max(t/am(jci+1),0.0_R8))*0.25_R8*1.0e4_R8
          t=tab2(nx,iy+1,jcf,1)/ev
          sjca(iy,2)=s(iy,2)*dab2(nx,iy+1,jcf,1)*
     *     sqrt(max(t/am(jci+1),0.0_R8))*0.25_R8*1.0e4_R8
        end if
      end do
      do j=1,2
        do i=1,ny
          sjca_t(j)=sjca_t(j)+sjca(i,j)
          sjha_t(j)=sjha_t(j)+sjha(i,j)
          sjhp_t(j)=sjhp_t(j)+sjhp(i,j)
          sjcp_t(j)=sjcp_t(j)+sjcp(i,j)
        end do
      end do
c      print '(a,1p/(4e12.3))','sjha_t,sjca_t,sjhp_t,sjcp_t : ',
c     ,                  (sjha_t(i),sjca_t(i),sjhp_t(i),sjcp_t(i),i=1,2)
      do j=1,2
        if(ntns .eq. 2) then
          js=ltns(j)+nlim
          sjha_w(j)=wldna(js,1,0)
          sjhp_w(j)=wldpp(js,1,0)
          sjht_w(j)=sjha_w(j)+sjhp_w(j)
          if(jcf .gt. 0) then
            sjca_w(j)=wldna(js,jcf,0)
            sjcp_w(j)=0
            do is=1,6
              sjcp_w(j)=sjcp_w(j)+wldpp(js,jci+is,0)
            end do
          sjct_w(j)=sjca_w(j)+sjcp_w(j)
          end if
c          print '(a,2i5,1p,4e12.3)',
c     ,              'j,js,sjha_w,sjca_w,sjhp_w,sjcp_w = ',
c     ,               j,js,sjha_w(j),sjca_w(j),sjhp_w(j),sjcp_w(j)
          if(sjha_t(j) .gt. 0.0_R8) then
            sjha_t(j)=sjha_w(j)/sjha_t(j)
          else
            sjha_t(j)=0.0_R8
          end if
          if(sjhp_t(j) .gt. 0.0_R8) then
            sjhp_t(j)=sjhp_w(j)/sjhp_t(j)
          else
            sjhp_t(j)=0.0_R8
          end if
          if(jcf .gt. 0) then
            if(sjca_t(j) .gt. 0.0_R8) then
              sjca_t(j)=sjca_w(j)/sjca_t(j)
            else
              sjca_t(j)=0.0_R8
            end if
            if(sjcp_t(j) .gt. 0.0_R8) then
              sjcp_t(j)=sjcp_w(j)/sjcp_t(j)
            else
              sjcp_t(j)=0.0_R8
            end if
          end if
        else
          sjha_t(j)=0.0_R8
          sjca_t(j)=0.0_R8
        end if
        do iy=0,ny-1
          sjha(iy,j)=sjha(iy,j)*sjha_t(j)
          sjhp(iy,j)=sjhp(iy,j)*sjhp_t(j)
          sjht(iy,j)=sjha(iy,j)+sjhp(iy,j)
          sjca(iy,j)=sjca(iy,j)*sjca_t(j)
          sjcp(iy,j)=sjcp(iy,j)*sjcp_t(j)
          sjct(iy,j)=sjca(iy,j)+sjcp(iy,j)
        end do
      end do
      print *,'Scaling factors for jD: ',sjha_t
      print *,'Scaling factors for jD+: ',sjhp_t
      print *,'Scaling factors for jC: ',sjca_t
      print *,'Scaling factors for jC+: ',sjcp_t
c
c*** Calculate sputtering
c
c      print *,'Sputtering...'
      if(jcf .gt. 0) then
        do iy=0,ny-1
          sjcs(iy,1)=sjht(iy,1)
          sjcs(iy,2)=sjht(iy,2)
          sjcs_t(1)=sjcs_t(1)+sjcs(iy,1)
          sjcs_t(2)=sjcs_t(2)+sjcs(iy,2)
        end do
      end if
      do j=1,2
        if(ntns .eq. 2) then
          js=ltns(j)+nlim
          sjcs_w(j)=wldspt(js,0)
          sjcb_w(j)=sjcs_w(j)-sjct_w(j)
c          print *,'j,js,sjcs_w,sjcb_w = ',j,js,sjcs_w(j),sjcb_w(j)
          if(sjcs_t(j) .gt. 0.0_R8) then
            sjcs_t(j)=sjcs_w(j)/sjcs_t(j)
          else
            sjcs_t(j)=0.0_R8
          end if
        else
          sjcs_t(j)=0.0_R8
        end if
        do iy=0,ny-1
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
        end do
      end do
      print *,'Scaling factors for Sp: ',sjcs_t
      nclm_er=10
      nx_hlp=1
      do j=1,2
        open(lun,file=ftrg(j,2))
        write(fmt_peak,710) nclm_er,tab
        write(lun,740) 'Target erosion - ',trg_label(j),label
        write(lun,fmt_peak) 'x','Sput_net','Sput_gross',
     ,            'Ersn_net','Ersn_gross','D_flux_atm','D_flux_ion',
     ,            'C_flux_atm','C_flux_ion','T_neut_D'
        write(fmt_peak,720) nclm_er,tab
        do i=0,ny-1
          write(lun,fmt_peak) xtrg(i,j),sjcb(i,j),sjcs(i,j),
     ,              sjcn(i,j),sjcg(i,j),sjha(i,j),sjhp(i,j),
     ,              sjca(i,j),sjcp(i,j),tab2(nx_hlp,i,1,1)/ev
        end do
        close (lun)
        nx_hlp=nx
      end do
c----------------------------------------------------------------------
c*** Produce a summary file for target sputtering
c
      open(lun,file=fsumm(2))
      write(lun,*) 'Target erosion summary: ',label
      write(lun,*)
      write(lun,'(12x,3(a1,a12))')
     -                     tab,'inner   ',tab,'outer   ',tab,'total   '
      write(fmt_peak,730) 3,tab
  730 format('(1p,1x,a11,',i3,'(1h',a1,',e12.4))')
      write(lun,*) 'Particle fluxes [1/s]'
      write(lun,fmt_peak) 'D ions on',sjhp_w,sjhp_w(1)+sjhp_w(2)
      write(lun,fmt_peak) 'D atoms on',sjha_w,sjha_w(1)+sjha_w(2)
      write(lun,fmt_peak) 'D total on',sjht_w,sjht_w(1)+sjht_w(2)
      write(lun,fmt_peak) 'C ions on',sjcp_w,sjcp_w(1)+sjcp_w(2)
      write(lun,fmt_peak) 'C atoms on',sjca_w,sjca_w(1)+sjca_w(2)
      write(lun,fmt_peak) 'C total on',sjct_w,sjct_w(1)+sjct_w(2)
      write(lun,fmt_peak) 'C sputtered',sjcs_w,sjcs_w(1)+sjcs_w(2)
      write(lun,fmt_peak) 'C removed',sjcb_w,sjcb_w(1)+sjcb_w(2)
      write(lun,*)
      write(lun,*) 'Erosion [mm/y]'
      ersn_pm=0.0_R8
      ersn_am=0.0_R8
      ss_m=0.0_R8
      do j=1,2
        ss(j)=0.0_R8
        ersn_p(j)=0.0_R8
        do iy=0,ny-1
          ss(j)=ss(j)+s(iy,j)
          if(sjcn(iy,j) .gt. ersn_p(j)) then
            ersn_p(j)=sjcn(iy,j)
            ersn_l(j)=xtrg(iy,j)
          end if
        end do
        ersn_a(j)=sjcb_w(j)*ers_factor/ss(j)
        ersn_pm=max(ersn_p(j),ersn_pm)
        ss_m=ss_m+ss(j)
        ersn_am=ersn_am+sjcb_w(j)
      end do
      ersn_am=ersn_am*ers_factor/ss_m
      write(lun,fmt_peak) 'average',ersn_a,ersn_am
      write(lun,fmt_peak) 'peak',ersn_p,ersn_pm
      write(lun,fmt_peak) 'peak pos.',ersn_l
c
      close (lun)
      return
c======================================================================
  710 format('(2h# ,',i3,'(1x,a10,1h',a1,'))')
  720 format('(1p,',i3,'(e12.4,1h',a1,'))')
  740 format('# ',2a,' : ',a)
      end subroutine b2ptrgl
