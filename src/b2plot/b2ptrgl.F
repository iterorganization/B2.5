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
      use eirmod_ccona
#endif
      implicit none
#include <DIMENSIONS.F>
      integer ismain, ismain0
      real (kind=R8) :: BoRiS
      integer lun
      parameter (lun=96)
#ifndef B25_EIRENE
      real (kind=R8) :: EionH2
      parameter (EionH2=15.4_R8)
#endif
c
      integer i,j,jj,k,l,m,n,iss,imol,iy,js,nclm_ld,nclm_er,nx_hlp
      integer i1,i2,i3,ia,ib,ic,is,ix,ik,imin,imax
      integer ixstart(ntrgsx,2),ixend(ntrgsx,2)
      integer m_hlp(ntrgsx)
      logical match_found
      real (kind=R8) ::
     *  h,t,u,a,b,d,p,wo,wt,mion,mmol,t1(2),
     ,  z_hlp(ntrgsx),xtg_hlp(-1:ny-1,ntrgsx),xmp_hlp(0:ny-1,ntrgsx),
     ,  xmc_hlp(0:ny-1,ntrgsx),s(0:ny-1,ntrgsx),wi(ntrgsx),
     ,  zctrg(0:ny-1,ntrgsx),rctrg(0:ny-1,ntrgsx),xctrg(0:ny-1,ntrgsx),
     ,  xtrg(0:ny-1,ntrgsx),netrg(0:ny-1,ntrgsx),tetrg(0:ny-1,ntrgsx),
     ,  titrg(0:ny-1,ntrgsx),wpttrg(0:ny-1,ntrgsx),
     ,  wpltrg(0:ny-1,ntrgsx),wnutrg(0:ny-1,ntrgsx),
     ,  wrdtrg(0:ny-1,ntrgsx,0:DEF_NATM),wtttrg(0:ny-1,ntrgsx),
     ,  wpk_tt(0:ntrgsx),wpk_rd(0:ntrgsx),
     ,  wpk_pl(0:ntrgsx),wpk_nu(0:ntrgsx),
     ,  wlc_tt(ntrgsx),wlc_rd(ntrgsx),wlc_pl(ntrgsx),wlc_nu(ntrgsx),
     ,  pwr_tt(ntrgsx),pwr_rd(ntrgsx),fwht_tt(ntrgsx),fwxt_tt(ntrgsx),
     ,  hwht_tt(ntrgsx),hwhm_tt(ntrgsx),hwhx_tt(ntrgsx),hwxt_tt(ntrgsx),
     ,  hwxm_tt(ntrgsx),lcnnt(0:ny-1,ntrgsx),lcnnx(0:ny-1,ntrgsx)
c
c*** Arrays for particle fluxes and erosion
c***  sja  impinging flux of atoms
c***  sjp  impinging flux of ions
c***  sjt  total flux of particles to the target
c***  sjb  "net" sputtering of target material (negative means deposition) [g/y]
c***  sjs  "gross" sputtering of target material [g/y]
c***  sjg  "gross" erosion rate [mm/y]
c***  sjn  "net" erosion rate [mm/y]
c***    suffix "_t" denotes totals from distributions,
c***    and "_w" the totals from Eirene.
c
      real (kind=R8) ::
     *  sja(0:ny-1,ntrgsx,natm), sjp(0:ny-1,ntrgsx,natm), 
     *  sjt(0:ny-1,ntrgsx,natm), spa(0:ny-1,ntrgsx,natm),
     *  sjb(0:ny-1,ntrgsx,ntargsp), sjg(0:ny-1,ntrgsx,ntargsp),
     *  sjn(0:ny-1,ntrgsx,ntargsp), sjs(0:ny-1,ntrgsx,ntargsp)
      real (kind=R8) ::
     *   sja_t(0:ntrgsx,natm), sjp_t(0:ntrgsx,natm),
     ,   sja_w(0:ntrgsx,natm), sjp_w(0:ntrgsx,natm), 
     ,   sjt_w(0:ntrgsx,natm), spt_out(4*natm),
     ,   sjs_t(ntrgsx,ntargsp), sjs_w(0:ntrgsx,ntargsp),
     ,   sjb_w(0:ntrgsx,ntargsp),sb_t,ss_t,sn_t,sg_t
      integer sumlis(natm,2)
      character*11 spt_str(4,natm)
      real (kind=R8) ::
     *   ers_factor, nsecyr
      parameter (nsecyr=24.0_R8*3600.0_R8*365.25_R8*0.9925_R8)
      parameter (ers_factor=nsecyr*1.0e-3_R8/avogr)
      real (kind=R8) ::
     *   ss(ntrgsx),ersn_a(ntrgsx),ss_m,ersn_am,ersn_p(ntrgsx),ersn_pm,
     *   ersn_l(ntrgsx)
      logical ft44
      logical streql
      external streql
      intrinsic sqrt
c*** Data for erosion calculations
c***  n_srfs_m  : number of possible target materials
c***  srfs_m    : list of possible target materials
c***  srfs_x    : corresponding species (say, Fe for SS)
c***  srfs_d    : specific masses of possible target materials [g/cm**3]
c***  srfs_a    : atomic weight of possible target materials
      integer n_srfs_m,spt_mat,spt_atm
      parameter (n_srfs_m=8)
      real (kind=R8) :: srfs_d(n_srfs_m),srfs_a(n_srfs_m)
      character srfs_m(n_srfs_m)*4,srfs_x(n_srfs_m)*2
      data srfs_m /'Be','C','CFC','Fe','SS','ss','Mo' ,'W'/
      data srfs_x /'Be','C','C'  ,'Fe','Fe','Fe','Mo' ,'W'/
      data srfs_a /  9.,12., 12. , 56., 56., 56., 96. , 184./
      data srfs_d /1.8 ,2.0, 2.0 , 7.8, 7.8, 7.8, 10.2, 19.3/
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
      inquire(file='fort.44',exist=ft44)

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
        do iy=0,lseppos(j) !{
          z_hlp(j)=z_hlp(j) - sqrt((crx(k,iy,3)-crx(k,iy,1))**2+
     +                             (cry(k,iy,3)-cry(k,iy,1))**2)
        end do !}
        xtg_hlp(-1,j)=z_hlp(j)
        do iy=0,ny-1 !{
          rctrg(iy,j)=0.5_R8*(crx(k,iy,1)+crx(k,iy,3))
          zctrg(iy,j)=0.5_R8*(cry(k,iy,1)+cry(k,iy,3))
          h=sqrt((crx(k,iy,3)-crx(k,iy,1))**2+
     +           (cry(k,iy,3)-cry(k,iy,1))**2)
          s(iy,j)=gs(k+1,iy,0)
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
c        print '(1x,i3,10f9.5/(4x,10f9.6))',lseppos(j)+1,
c     ,                                         (xtg_hlp(i,j),i=0,ny-1)
c      end do !}
c!###}

c*** Mid-plane co-ordinates

      l=jxi
      m=0
      xmc_hlp=0.0_R8
      do j=1,ntrgts !{
        u=0.0_R8
        k=lpfrpos(j)
        n=1-m
        do i=0,lseppos(j) !{
          h=0.5_R8*sqrt((crx(k,i,n+2)-crx(k,i,n))**2+
     +                  (cry(k,i,n+2)-cry(k,i,n))**2)
          u=u+h
          xmp_hlp(i,j)=u
          u=u+h
        end do !}
        do i=0,lseppos(j) !{
          xmp_hlp(i,j)=xmp_hlp(i,j)-u
        end do !}
        u=0.0_R8
        k=l
        xmc_hlp(lseppos(j),j)=0.0_R8
        do i=lseppos(j)+1,ny-1 !{
          h=0.5_R8*sqrt((crx(k,i,3)-crx(k,i,1))**2+
     +                  (cry(k,i,3)-cry(k,i,1))**2)
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
c      print *,'ny,xmp_hlp:',ny,jxi,jxa
c      do j=1,ntrgts !{
c        print '(1x,2i3,10f9.5/(4x,10f9.6))',lseppos(j)+1,lpfrpos(j),
c     ,                                         (xmp_hlp(i,j),i=0,ny-1)
c      end do !}
c!###}
c-----------------------------------------------------------------------
c*** Power loading:

c
c*** Plasma and neutrals
c

      l=1
      u=-1.0_R8
      ixstart(1,1)=jxi
      ixend(1,1)=0
      ixstart(ntrgts,1)=jxa
      ixend(ntrgts,1)=nx
      ixstart(1,2)=leftcut(1)-1
      ixend(1,2)=0
      ixstart(ntrgts,2)=rightcut(1)
      ixend(ntrgts,2)=nx
      if (ntrgts.gt.2) then
        ixstart(2,1)=jxi
        ixend(2,1)=nxtl
        ixstart(3,1)=jxa
        ixend(3,1)=nxtr
        ixstart(2,2)=leftcut(2)
        ixend(2,2)=nxtl
        ixstart(2,2)=rightcut(2)-1
        ixend(3,2)=nxtr
      endif
      do j=1,ntrgts !{
        k=ltrgpos(j)
        do iy=lytg,ny-1 !{
          wnutrg(iy,j)=0.0_R8
          wo=fhe(k+1-l,iy,0)+fhi(k+1-l,iy,0)+fhi_ext(k+1-l,iy,0)
          do is=0,ns-1 !{
            if (ft44) then
              wo=wo+fna(k+1-l,iy,0,is)*(rpt(k+l,iy,is)*ev
     .           +(1.0_R8-BoRiS)*(kinrgy(k+l,iy,is)+ti(k+l,iy)))
c carlis 31jan2017  added ti as suggested by sw in b2blne.F
            else
              if (is_neutral(is)) then
                wnutrg(iy,j)=wnutrg(iy,j)+fna(k+1-l,iy,0,is)
     .           *(1.0_R8-BoRiS)*(kinrgy(k+l,iy,is)+ti(k+l,iy))
              else
                wo=wo+fna(k+1-l,iy,0,is)*(rpt(k+l,iy,is)*ev
     .           +(1.0_R8-BoRiS)*(kinrgy(k+l,iy,is)+ti(k+l,iy)))
c carlis 31jan2017  added ti as suggested by sw in b2blne.F
              endif
            endif
          end do !}
          do is=0,ns_ext-1 !{
            if (ft44) then
              wo=wo+fa_ext(k+1-l,iy,0,is)*(pt_ext(k+l,iy,is)*ev
     .           +(1.0_R8-BoRiS)*(0.5_R8*am_ext(is)*mp*
     .             ua_ext(k+l,iy,is)**2+ta_ext(k+l,iy,is)))
c carlis 31jan2017  added ta as suggested by sw in b2blne.F
            else
              if (is_ext_neutral(is)) then
                wnutrg(iy,j)=wnutrg(iy,j)+fa_ext(k+1-l,iy,0,is)
     .           *(1.0_R8-BoRiS)*(0.5_R8*am_ext(is)*mp*
     .             ua_ext(k+l,iy,is)**2+ta_ext(k+l,iy,is))
              else
                wo=wo+fa_ext(k+1-l,iy,0,is)*(pt_ext(k+l,iy,is)*ev
     .           +(1.0_R8-BoRiS)*(0.5_R8*am_ext(is)*mp*
     .             ua_ext(k+l,iy,is)**2+ta_ext(k+l,iy,is)))
c carlis 31jan2017  added ta as suggested by sw in b2blne.F
              endif
            end if
          end do !}
          wo=wo+(1.0_R8-BoRiS)*fne(k+1-l,iy,0)*te(k+l,iy)
c carlis 31jan2017  added electron contribution as suggested by sw in b2blne.
c        wo=-1.0e-6_R8*wo/s(iy,1)
          wo= wo/s(iy,j)
          wpltrg(iy,j)=wo*u
          if(.not.ft44) wnutrg(iy,j)=wnutrg(iy,j)*u/s(iy,j)
          netrg(iy,j)=ne(k+l,iy)
          tetrg(iy,j)=te(k+l,iy)/ev
		  titrg(iy,j)=ti(k+l,iy)/ev
          if (ft44) then
            do is=1,natmi
              t=tab2(k+l+1,iy+1,is,1)/ev
              wnutrg(iy,j)=wnutrg(iy,j)+
     .         s(iy,j)*dab2(k+l+1,iy+1,is,1)*tab2(k+l+1,iy+1,is,1)*
     *         sqrt(max(t/am(eb2spcr(is)),0.0_R8))*0.25_R8*2.0e4_R8
            enddo
            do is=1,nmoli
              mmol=0.0_R8
              do iss = 1,natmi
                mmol=mmol+mlcmp(iss,is)*am(eb2spcr(iss))
              enddo
              t=tmb2(k+l+1,iy+1,is,1)/ev
              wnutrg(iy,j)=wnutrg(iy,j)+
     .         s(iy,j)*dmb2(k+l+1,iy+1,is,1)*tmb2(k+l+1,iy+1,is,1)*
     *         sqrt(max(t/mmol,0.0_R8))*0.25_R8*2.0e4_R8
            enddo
            do is=1,nioni
              mion=0.0_R8
              imol=1
              match_found=.false.
              do while (.not.match_found.and.imol.le.nmoli)
                if (streql(textin(is-1),textmn(imol-1)//'+')) then
                  match_found=.true.
                  do iss = 1,natmi
                    mion=mion+mlcmp(iss,imol)*am(eb2spcr(iss))
                  enddo
                endif
                imol=imol+1
              enddo
              t=tib2(k+l+1,iy+1,is,1)/ev
              if (mion.gt.0.0_R8)
     >         wnutrg(iy,j)=wnutrg(iy,j)+
     .          s(iy,j)*dib2(k+l+1,iy+1,is,1)*
     *          sqrt(max(t/mion,0.0_R8))*0.25_R8*1.0e4_R8*
     .          (2.0_R8*tib2(k+l+1,iy+1,is,1)+EionH2*ev)
            enddo
          endif
c carlis  21 feb 2017 added calculation of connection length from midplane to targets
c         and from x-point to targets
          do ik=1,2
            ix=ixstart(j,ik)
            t1(ik)=0.0_R8
            if (ik.eq.1.and.iy.gt.jsep) then
c carlis If starting @midplane in SOL then remove half of starting cell
              t1(ik)=t1(ik)-hx(ix,iy)/2.0_R8/
     1             (sqrt(bb(ix,iy,0)**2+bb(ix,iy,1)**2)/bb(ix,iy,3))            
            end if
            if (l.eq.1) then
              do ix=ixend(j,ik),ixstart(j,ik)
                t1(ik)=t1(ik)+hx(ix,iy)/
     1                (sqrt(bb(ix,iy,0)**2+bb(ix,iy,1)**2)/bb(ix,iy,3))
              end do	  
            else
              do ix=ixstart(j,ik),ixend(j,ik)
                t1(ik)=t1(ik)+hx(ix,iy)/
     1                (sqrt(bb(ix,iy,0)**2+bb(ix,iy,1)**2)/bb(ix,iy,3))
              end do			  
            endif
          end do
          lcnnt(iy,j)=t1(1)
          lcnnx(iy,j)=t1(2)
        end do !}
        l=1-l
        u=-u
      end do !}
      do iy=lytg,jsep !{
c carlis  If in PFR region then connection length from midplane to target will
c         be equal to connection length from target to target
        lcnnt(iy,1)=lcnnx(iy,1)+lcnnx(iy,2)
        lcnnt(iy,ntrgts)=lcnnx(iy,ntrgts)+lcnnx(iy,ntrgts-1)
        if (ntrgts.gt.2) then
          lcnnt(iy,2)=lcnnt(iy,1)
          lcnnt(iy,3)=lcnnt(iy,ntrgts)
        end if
      end do

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
        if (ft44) then
          if(ntns .eq. ntrgts) then !{
            js=ltns(j)+nlim
            wo=wldnek(js,0)+wldnep(js,0)
c          print *,'j,js,wldnek,wldnep = ',j,js,wldnek(js,0),wldnep(js,0)
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
        endif
        do iy=lytg,ny-1 !{
c          if (ft44) wnutrg(iy,j)=1.0e-6_R8*wnutrg(iy,j)*wi(j)/s(iy,j)
          if (ft44) wnutrg(iy,j)=wnutrg(iy,j)*wi(j)/s(iy,j)
          wpttrg(iy,j)=wnutrg(iy,j)+wpltrg(iy,j)
          wtttrg(iy,j)=wpttrg(iy,j)+wrdtrg(iy,j,0)
        end do !}
      end do !}
c      write(*,'(1x,a,1p,16e12.3)')
c     -                   'Scaling factors for Wn : ',(wi(i),i=1,ntrgts)
c
      nclm_ld=14
      if(lrdlsp) nclm_ld=nclm_ld+natmi
      do j=1,ntrgts !{
        jj=j
        if(j.eq.ntrgts) jj=ntrgsx
        open(lun,file=ftrg(jj,1))
        write(fmt_peak,710) nclm_ld,tab
        write(lun,740) 'Target loading - ',trg_label(jj),label
        write(lun,740) 'x',tab,'co-ordinate along target [m]. <0 in PFR'
        write(lun,740) 'ne',tab,'electron density [1/m^3]'
        write(lun,740) 'Te',tab,'electron temperature [eV]'
        write(lun,740) 'Ti',tab,'ion temperature [eV]'
        write(lun,740) 'Wtot',tab,'total power flux [W/m^2]'
        write(lun,740) 'Wpart',tab,'power flux with particles [W/m^2]'
        write(lun,740) 'Wrad',tab,'power flux with radiation [W/m^2]'
        write(lun,740) 'Wpls',tab,'power flux with charged particles
     ,[W/m^2]'
        write(lun,740) 'Wneut',tab,'power flux with neutral particles 
     ,[W/m^2]'
        write(lun,740) 'xMP',tab,'radial co-ordinate in mid-plane, 
     ,(IMP for inner target(s), OMP for outer(s)) [m]'
        write(lun,740) 'Rclfc',tab,'radial co-ordinate on target [m]'
        write(lun,740) 'Zclfc',tab,'vertical co-ordinate on target [m]'
c        write(lun,740) 'Wpar_xpt',tab,'parallel power flux around 
c     ,x-point [W/m^2]'
c        write(lun,740) 'WWpar',tab,'would-be parallel power flux 
c     ,without dissipation in divertor [W/m^2]'
c        write(lun,740) 'WWtrg',tab,'would-be target power flux 
c     ,without dissipation in divertor [W/m^2]'
        write(lun,740) 'Lcnnt',tab,'connection length from mid-plane 
     ,to target [m]'
        write(lun,740) 'Lcnnx',tab,'connection length from x-point 
     ,to target [m]'
        if(lrdlsp) then !{
          write(lun,fmt_peak) 'x','ne','Te','Ti','Wtot','Wpart',
     ,            'Wrad','Wpls','Wneut','xMP','Rclfc','Zclfc',
     ,            'Lcnnt','Lcnnx',
     ,            ('rad_'//textan(is)(2:3),is=0,natmi-1)
        else !}{
          write(lun,fmt_peak) 'x','ne','Te','Ti','Wtot','Wpart',
     ,            'Wrad','Wpls','Wneut','xMP','Rclfc','Zclfc',
     ,            'Lcnnt','Lcnnx'
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
     ,              titrg(i,j),wtttrg(i,j),wpttrg(i,j),wrdtrg(i,j,0),
     ,              wpltrg(i,j),wnutrg(i,j),xmp_hlp(i,k),rctrg(i,j),
     ,              zctrg(i,j),lcnnt(i,j),lcnnx(i,j),
     ,              (wrdtrg(i,j,is),is=1,natmi)
          else !}{
            write(lun,fmt_peak) xtrg(i,j),netrg(i,j),tetrg(i,j),
     ,              titrg(i,j),wtttrg(i,j),wpttrg(i,j),wrdtrg(i,j,0),
     ,              wpltrg(i,j),wnutrg(i,j),xmp_hlp(i,k),rctrg(i,j),
     ,              zctrg(i,j),lcnnt(i,j),lcnnx(i,j)
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
c     carlis 22 feb 2017 calculating startiing/ending indexes of each plasma species      
      do ia=1,natm
        imin=1000
        imax=0
        do is=0,ns-1
          if (b2espcr(is).eq.ia) then
            imin=min(imin,is);
            imax=max(imax,is);
          end if
        end do
        sumlis(ia,1)=imin
        sumlis(ia,2)=imax
      end do

c*** Calculate the particle fluxes

c      print *,'Particle fluxes...' !###
      if(ntns.eq.ntrgts) then !{
        do j=1,ntrgts !{
          do i=0,ny-1 !{
            do ia=1,natm
               sja(i,j,ia)=0.0_R8
               sjp(i,j,ia)=0.0_R8
               sjt(i,j,ia)=0.0_R8
            end do
          end do !}
          ersn_a(j)=0.0_R8
          ersn_p(j)=0.0_R8
          ersn_l(j)=0.0_R8   
          do is=1,ntargsp
            sjs_w(j,is)=0.0_R8
            sjb_w(j,is)=0.0_R8
            sjs_t(j,is)=0.0_R8
          end do    
          do ia=1,natm
            sja_w(j,ia)=0.0_R8
            sjp_w(j,ia)=0.0_R8
            sjt_w(j,ia)=0.0_R8
            sja_t(j,ia)=0.0_R8
            sjp_t(j,ia)=0.0_R8
          end do      
        end do !}

c        do is=0,ns-1
c          print *,'is,am(is)',is,am(is)
c        end do
c        js=ltns(1)+nlim
c        do is=0,ns+1
c          print *,'is,wldna(js,is,0)',is,wldna(js,is,0)
c          print *,'is,wldpp(js,is,0)',is,wldpp(js,is,0)
c        end do
c        print *,'ismain,ismain0',ismain,ismain0
        u=-1.0_R8
        l=1
        do j=1,ntrgts !{
          k=ltrgpos(j)
c          print *,'j,ltrgpos = ',j,k  !###
          do iy=lytg,ny-1 !{
            do ia=1,natm
              do is=sumlis(ia,1)+1,sumlis(ia,2)
                sjp(iy,j,ia)=sjp(iy,j,ia)+u*fna(k+1-l,iy,0,is)
              end do
              if (ft44) then
                t=tab2(k+l+1,iy+1,ia,1)/ev
                sja(iy,j,ia)=s(iy,2)*dab2(k+l+1,iy+1,ia,1)*
     *                       sqrt(max(t/am(sumlis(ia,1)),0.0_R8))*
     *                       0.25_R8*1.0e4_R8
                spa(iy,j,ia)=dab2(k+l+1,iy+1,ia,1)*
     *                       tab2(k+l+1,iy+1,ia,1)
              else
                sja(iy,j,ia)=u*fna(k+1-l,iy,0,sumlis(ia,1))
                spa(iy,j,ia)=na(k+l,iy,sumlis(ia,1))*ti(k+l,iy)
              endif
              sjp_t(j,ia)=sjp_t(j,ia)+sjp(iy,j,ia)
              sja_t(j,ia)=sja_t(j,ia)+sja(iy,j,ia)
            end do
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
c        print *,'size wldpp(,,)',size(wldpp,1),size(wldpp,2),size(wldpp,3)
        if (ft44) then
          do j=1,ntrgts !{
            js=ltns(j)+nlim
            do ia=1,natm
              sja_w(j,ia)=wldna(js,ia,0)
              do is=sumlis(ia,1)+2-ia,sumlis(ia,2)+1-ia
                sjp_w(j,ia)=sjp_w(j,ia)+wldpp(js,is,0)
              end do
              sjt_w(j,ia)=sja_w(j,ia)+sjp_w(j,ia)
              if(sja_t(j,ia).gt.0.0_R8) then !{
                sja_t(j,ia)=sja_w(j,ia)/sja_t(j,ia)
              else !}{
                sja_t(j,ia)=0.0_R8
              end if !}
              if(sjp_t(j,ia).gt.0.0_R8) then !{
                sjp_t(j,ia)=sjp_w(j,ia)/sjp_t(j,ia)
              else !}{
                sjp_t(j,ia)=0.0_R8
              end if !}
              do iy=lytg,ny-1 !{
                sja(iy,j,ia)=sja(iy,j,ia)*sja_t(j,ia)
                sjp(iy,j,ia)=sjp(iy,j,ia)*sjp_t(j,ia)
                sjt(iy,j,ia)=sja(iy,j,ia)+sjp(iy,j,ia)
              end do !}
            end do
          end do !}
        end if
c        write(*,'(1x,a,1p,16e12.3)') 'Scaling factors for jD : ',
c     ,                                           (sjha_t(i),i=1,ntrgts)
c        write(*,'(1x,a,1p,16e12.3)') 'Scaling factors for jD+: ',
c     ,                                           (sjhp_t(i),i=1,ntrgts)
c        write(*,'(1x,a,1p,16e12.3)') 'Scaling factors for jC : ',
c     ,                                           (sjca_t(i),i=1,ntrgts)
c        write(*,'(1x,a,1p,16e12.3)') 'Scaling factors for jC+: ',
c     ,                                           (sjcp_t(i),i=1,ntrgts)
c        write(*,'(1x,a,1p,16e12.3)') 'Scaling factors for jHe : ',
c     ,                                           (sjhea_t(i),i=1,ntrgts)
c        write(*,'(1x,a,1p,16e12.3)') 'Scaling factors for jHe+: ',
c     ,                                           (sjhep_t(i),i=1,ntrgts)

c*** Calculate sputtering

        if (ntargsp.gt.0.and.ft44) then !sputtering happening
          do is=1,ntargsp
c           Check if sputtering is happening and which material is being sputtered
c           and its atomic mass and density
            print *,'sputtered target material = ', track_species(is)
            do ia=1,n_srfs_m
             if (streql(track_species(is),srfs_m(ia))) then
               spt_mat=ia
             end if
            end do
            spt_atm=0
            do ia=1,natm
              if (streql(track_species(is),atoms_list(ia))) then
                spt_atm=ia
              end if
            end do
            do j=1,ntrgts
              do iy=lytg,ny-1
                sjs(iy,j,is)=sjt(iy,j,spt_atm)
                sjs_t(j,is)=sjs_t(j,is)+sjs(iy,j,is)
              end do
              if(ntns.eq.2 .or. ntns.eq.4) then
                js=ltns(j)+nlim
                sjs_w(j,is)=wldspt(js,0)
                sjb_w(j,is)=sjs_w(j,is)-sjt_w(j,spt_atm)
                if(sjs_t(j,is).gt.0.0_R8) then
                  sjs_t(j,is)=sjs_w(j,is)/sjs_t(j,is)
                else
                  sjs_t(j,is)=0.0_R8
                end if
              else
                sjs_t(j,is)=0.0_R8
              end if
              do iy=lytg,ny-1!{
                sjs(iy,j,is)=sjs_t(j,is)*sjs(iy,j,is)
                sjb(iy,j,is)=(sjs(iy,j,is)-sjt(iy,j,spt_atm))/s(iy,j)
                sjn(iy,j,is)=ers_factor*srfs_a(spt_mat)*sjb(iy,j,is)/
     .                     srfs_a(spt_mat)
                sjs(iy,j,is)=sjs(iy,j,is)/s(iy,j)
                sjg(iy,j,is)=ers_factor*srfs_a(spt_mat)*sjs(iy,j,is)/
     .                     srfs_a(spt_mat)
                do ia=1,natm
                  sjp(iy,j,ia)=sjp(iy,j,ia)/s(iy,j)
                  sjt(iy,j,ia)=sjt(iy,j,ia)/s(iy,j) 
                  sja(iy,j,ia)=sja(iy,j,ia)/s(iy,j)
                end do
              end do
            end do
            write(*,'(1x,a,1p,16e12.3)') 
     .    'Scaling factors for '//track_species(is)//' Sp : ',
     ,    (sjs_t(j,is),j=1,ntrgts)
          end do
        end if
        
        nclm_er=6+4*natm
        l=1
	do ia=1,natm
          spt_str(1,ia)=trim(atoms_list(ia))//'_flux_atm'
          spt_str(2,ia)=trim(atoms_list(ia))//'_flux_ion'
          spt_str(3,ia)=trim(atoms_list(ia))//'_temp_atm'
          spt_str(4,ia)=trim(atoms_list(ia))//'_pres_atm'
        end do
        do j=1,ntrgts !{
          nx_hlp=ltrgpos(j)+l
          jj=j
          if(j.eq.ntrgts) jj=ntrgsx
           open(lun,file=ftrg(jj,2))
           write(fmt_peak,710) nclm_er,tab
           write(lun,740) 'Target erosion - ',trg_label(jj),label
           write(lun,740) 'x',tab,'co-ordinate along target [m].
     , <0 in PFR'
           write(lun,740) 'Sput_net',tab,'net sputtering [g/y] 
     ,(negative means deposition)'
           write(lun,740) 'Sput_gross',tab,'gross sputtering [g/y]'
           write(lun,740) 'Ersn_net',tab,'net erosion [mm/y] 
     ,(negative means deposition)'
           write(lun,740) 'Ersn_gross',tab,'gross erosion [mm/y]'
           write(lun,740) 'XX_flux_atm',tab,'flux of neutral atoms 
     ,of species XX [1/m^2/s]'
           write(lun,740) 'XX_flux_ion',tab,'flux of ions of species 
     ,XX [1/m^2/s]'
           write(lun,740) 'XX_temp_atm',tab,'temperature of neutral 
     ,atoms of species XX [eV]'
           write(lun,740) 'XX_pres_atm',tab,'pressure of neutral 
     ,atoms of species XX [Pa]'
           write(lun,740) 'xMP',tab,'radial co-ordinate in mid-plane, 
     ,(IMP for inner target(s), OMP for outer(s)) [m]'
           write(lun,fmt_peak) 'x','Sput_net','Sput_gross',
     ,            'Ersn_net','Ersn_gross',spt_str,'xMP'
           write(fmt_peak,720) nclm_er,tab
           if(ntrgts.eq.4) then !{
             k=j
           else !}{
             k=ntrgts
           end if !}
           do i=lytg,ny-1 !{
             do ia=1,natm
               spt_out(1+4*(ia-1))=sja(i,j,ia)
               spt_out(2+4*(ia-1))=sjp(i,j,ia)
               spt_out(3+4*(ia-1))=tab2(nx_hlp+1,i+1,ia,1)/ev
               spt_out(4+4*(ia-1))=spa(i,j,ia)
             end do
             sb_t=0.0_R8
             ss_t=0.0_R8
             sn_t=0.0_R8
             sg_t=0.0_R8
             if (ntargsp.gt.0.and.ft44) then
               do is=1,ntargsp
                 sb_t=sb_t+sjb(i,j,is)
                 ss_t=ss_t+sjs(i,j,is)
                 sn_t=sn_t+sjn(i,j,is)
                 sg_t=sg_t+sjg(i,j,is)
               end do
             end if
             write(lun,fmt_peak) xtrg(i,j),sb_t,ss_t,sn_t,sg_t,
     ,                           spt_out,xmp_hlp(i,k)
           end do !}
           close (lun)
           l=1-l
         end do !}
c----------------------------------------------------------------------
c*** Produce a summary file for target sputtering
c
        do ia=1,natm
          sjp_w(0,ia)=0.0_R8
          sja_w(0,ia)=0.0_R8
          sjt_w(0,ia)=0.0_R8
        end do
        do is=1,ntargsp
          sjs_w(0,is)=0.0_R8
          sjb_w(0,is)=0.0_R8
        end do
        do j=1,ntrgts
          do ia=1,natm
            sjp_w(0,ia)=sjp_w(0,ia)+sjp_w(j,ia)
            sja_w(0,ia)=sja_w(0,ia)+sja_w(j,ia)
            sjt_w(0,ia)=sjt_w(0,ia)+sjt_w(j,ia)
          end do
        end do !}
        k=1     
        if(ntrgts.eq.2) k=ntrgsx-1
          open(lun,file=fsumm(2))
          write(lun,*) 'Target erosion summary: ',label
          write(lun,*)
          write(lun,'(12x,3(a1,a12))')
     -                   (tab,trg_label(i),i=1,ntrgsx,k),tab,'total   '
          write(fmt_peak,730) ntrgts+1,tab
  730     format('(1p,1x,a12,',i3,'(1h',a1,',e12.4))')
          write(lun,*) 'Particle fluxes [1/s]'
          do ia=1,natm
            write(lun,fmt_peak) trim(atoms_list(ia))//' ions on',
     ,                          (sjp_w(i,ia),i=1,ntrgts),sjp_w(0,ia)
            write(lun,fmt_peak) trim(atoms_list(ia))//' atoms on',
     ,                          (sja_w(i,ia),i=1,ntrgts),sja_w(0,ia)
            write(lun,fmt_peak) trim(atoms_list(ia))//' total on',
     ,                          (sjt_w(i,ia),i=1,ntrgts),sjt_w(0,ia)
          end do
          if (ntargsp.gt.0.and.ft44) then
            do is=1,ntargsp
              do j=1,ntrgts
                sjs_w(0,is)=sjs_w(0,is)+sjs_w(j,is)
                sjb_w(0,is)=sjb_w(0,is)+sjb_w(j,is)
              end do
            end do
            do is=1,ntargsp
              write(lun,fmt_peak) trim(track_species(is))//' sputtered',
     ,                          (sjs_w(i,is),i=1,ntrgts),sjs_w(0,is)
              write(lun,fmt_peak) trim(track_species(is))//' removed',
     ,                          (sjb_w(i,is),i=1,ntrgts),sjb_w(0,is)
            end do
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
                sn_t=0.0_R8
                do is=1,ntargsp
                  sn_t=sn_t+sjn(i,j,is)
                end do
                if(sn_t .gt. ersn_p(j)) then !{
                  ersn_p(j)=sn_t
                  ersn_l(j)=xtrg(iy,j)
                end if !}
              end do !}
              ersn_a(j)=0.0_R8
              do is=1,ntargsp
                ersn_a(j)=ersn_a(j)+sjb_w(j,is)*ers_factor*
     .           srfs_a(spt_mat)/srfs_d(spt_mat)/ss(j)
              end do
              ersn_pm=max(ersn_p(j),ersn_pm)
              ss_m=ss_m+ss(j)
              ersn_am=ersn_am+ersn_a(j)*ss(j)
            end do !}
            ersn_am=ersn_am/ss_m
            write(lun,fmt_peak) 'average',(ersn_a(i),i=1,ntrgts),ersn_am
            write(lun,fmt_peak) 'peak',(ersn_p(i),i=1,ntrgts),ersn_pm
            write(lun,fmt_peak) 'peak pos.',(ersn_l(i),i=1,ntrgts)
            close (lun)
          end if
        else
          write(*,*) 'b2ptrgl: ntns.ne.ntrgts (',ntns,' vs.',ntrgts,
     ,                                 ' ) - no sputtering calculated.'
          write(*,*) '  ==> check specification of ltns in ',
     ,                               'b2.parameters (section NEUTRALS)'
        end if    

      return
c======================================================================
  710 format('(2h# ,',i3,'(1x,a11,1h',a1,'))')
  720 format('(1p,',i3,'(e12.4,1h',a1,'))')
  740 format('# ',2a,' : ',a)
  750 format('# ',4a)
c=======================================================================
      end subroutine b2ptrgl
