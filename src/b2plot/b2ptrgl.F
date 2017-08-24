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
      character*8 TEXTA(DEF_NATM)
#endif
c
      integer i,j,jj,k,l,m,n,iss,imol,iy,js,nclm_ld,nclm_er,nx_hlp
      integer i1,i2,i3,ia,ib,ic,is,ix,ik,imin,imax,ind,ims,ias
      integer ixstart(ntrgsx,2),ixend(ntrgsx,2)
      integer m_hlp(ntrgsx)
      logical match_found,do_sput(ntrgsx)
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
c***  sjm  impinging flux of molecules
c***  sjs  sputtering of target material [g/y]
c***  sje  erosion rate [mm/y]
c***    suffix "_t" denotes totals from distributions,
c***    and "_w" the totals from Eirene.
c
c*** Recap of fluxes from eirene
c***    wldnek  :   heat transferred with neutrals. Evaluated as EOTAT+EOTML+EOTIO
c***                   -erfaat-erfmat-erfiat-erfaml-erfmml-erfiml
c***    wldnep  :   potential energy released by neutrals. Evaluated as diss_pot_H2*prfaml
c***    wldna   :   flux of atoms impinging onto the surface. Equal to POTAT
c***    ewlda   :   their average energy. Equal to EOTAT/POTAT*wldna
c***    wldnm   :   flux of molecules impinging onto the surface. Equal to POTML
c***    ewldm   :   their average energy. Equal to EOTML/POTML*wldnm
c***    wldra   :   flux of reflected atoms. Equal to PRFAAT+PRFMAT+PRFIAT
c***    wldrm   :   flux of reflected molecules. Equal to PRFAML+PRFMML+PRFIML
c***    wldpp   :   flux of plasma ions impinging onto the surface. Equal to POTPL
c***    wldpa   :   flux of emitted atoms from BI. Equal to PRFPAT
c***    wldpm   :   flux of emitted molecules from BI. Equal to PRFPML
c***    wldpeb  :   power carried away by these atoms and molecules. Equal to ERFPAT+ERFPML+ERFPIO
c***    wldspt  :   flux of sputtered wall material
c***    wldspta :   sputtered flux by emitted atom
c***    wldsptm :   sputtered flux by emitted molecule
c***    wlarea  :   areas of the surface segments from Eirene
c***    wlabsrp :   absorption at the surfaces (1-recyct from Eirene)
c***    ewldea_res  :   emitted energy flux of atoms. Equal to ERFAAT+ERFMAT+ERFIAT+ERFPHAT
c***    ewldem_res  :   emitted energy flux due of molecules. Equal to ERFAML+ERFMML+
c***                      ERFIML+ERFPHML
c***    ewldrp_res  :   kinetic energy of reflected atoms+molecules originated from ions. ERFPAT+ERFPML
c***    ewldmr_res  :   energy due to recombination of atoms into molecules. ()PRFAML+PRFPML)*H_EDISS

      real (kind=R8) ::
     *  sja(0:ny-1,ntrgsx,DEF_NATM), sjp(0:ny-1,ntrgsx,DEF_NATM), 
     *  sjt(0:ny-1,ntrgsx,DEF_NATM), spa(0:ny-1,ntrgsx,DEF_NATM),
     *  sje(0:ny-1,ntrgsx), sjs(0:ny-1,ntrgsx),
     *  sjm(0:ny-1,ntrgsx,DEF_NMOL),spm(0:ny-1,ntrgsx,DEF_NMOL)
      real (kind=R8) ::
     *   sja_t(0:ntrgsx,DEF_NATM), sjp_t(0:ntrgsx,DEF_NATM),
     ,   sja_w(0:ntrgsx,DEF_NATM), sjp_w(0:ntrgsx,DEF_NATM),
     ,   sjm_t(0:ntrgsx,DEF_NMOL), sjm_w(0:ntrgsx,DEF_NMOL), 
     ,   sjt_w(0:ntrgsx,DEF_NATM), fpt_out(4*DEF_NATM),
     ,   sjs_t(ntrgsx), sjs_w(0:ntrgsx),sptmass,
     ,   sb_t,ss_t,sn_t,sg_t, spt_out(2),mol_out(2*DEF_NMOL)
      integer eb2atm(natm)
      character*11 spt_str(2,ntrgsx),fpt_str(4,DEF_NATM)
      character*15 mol_str(2,DEF_NMOL)
      character*10 ssmol,ssion
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
      integer n_srfs_m,spt_mat(ntrgsx),sptlis(ntrgsx)
      parameter (n_srfs_m=7)
      real (kind=R8) :: srfs_d(n_srfs_m),srfs_a(n_srfs_m)
      character srfs_m(n_srfs_m)*4,srfs_x(n_srfs_m)*2
      data srfs_m /'Be','C','CFC','Fe','SS','Mo' ,'W'/
      data srfs_x /'Be','C','C'  ,'Fe','Fe','Mo' ,'W'/
      data srfs_a /  9.,12., 12. , 56., 56., 96. , 184./
      data srfs_d /1.8 ,2.0, 2.0 , 7.8, 7.8, 10.2, 19.3/
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

c*** Midplane coordinates

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
        ISS=ltns(j)
        IND=eirdiag_nds_ind(ISS)
        IAS=IND+1
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
c             carlis Using values directly calculated from Eirene
              wnutrg(iy,j)=wnutrg(iy,j)+(ewlda_res(is,iy+IAS)
     .                                 -ewldea_res(is,iy+IAS))*
     .                                  sarea_res(iy+IAS)
            enddo
            do is=1,nmoli
c             carlis Using values directly calculated from Eirene
              wnutrg(iy,j)=wnutrg(iy,j)+(ewldm_res(is,iy+IAS)
     .                                 -ewldem_res(is,iy+IAS)
     .                                 +ewldmr_res(is,iy+IAS))*
     .                                  sarea_res(iy+IAS)
            enddo
            do is=1,nioni
              mion=0.0_R8
              imol=1
              ssion=textin(is-1)
              ssion=trim(adjustl(ssion))
              match_found=.false.
              do while (.not.match_found.and.imol.le.nmoli)
                ssmol=textmn(imol-1)
                ssmol=trim(adjustl(ssmol))//'+'
                if (streql(ssion,ssmol)) then
                  match_found=.true.
                  do ims = 1,natmi
                    mion=mion+mlcmp(ims,imol)*am(eb2spcr(ims))
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
        write(lun,740) 'x',tab,'coordinate along target [m]. <0 in PFR'
        write(lun,740) 'ne',tab,'electron density [1/m^3]'
        write(lun,740) 'Te',tab,'electron temperature [eV]'
        write(lun,740) 'Ti',tab,'ion temperature [eV]'
        write(lun,740) 'Wtot',tab,'total power flux [W/m^2]'
        write(lun,740) 'Wpart',tab,'power flux with particles [W/m^2]'
        write(lun,740) 'Wrad',tab,'power flux with radiation [W/m^2]'
        write(lun,740) 'Wpls',tab,
     .   'power flux with charged particles [W/m^2]'
        write(lun,740) 'Wneut',tab,
     .   'power flux with neutral particles [W/m^2]'
        write(lun,740) 'xMP',tab,'radial coordinate in midplane '//
     .   '(IMP for inner target(s), OMP for outer(s)) [m]'
        write(lun,740) 'Rclfc',tab,'radial coordinate on target [m]'
        write(lun,740) 'Zclfc',tab,'vertical coordinate on target [m]'
c        write(lun,740) 'Wpar_xpt',tab,'parallel power flux around 
c     ,x-point [W/m^2]'
c        write(lun,740) 'WWpar',tab,'would-be parallel power flux 
c     ,without dissipation in divertor [W/m^2]'
c        write(lun,740) 'WWtrg',tab,'would-be target power flux 
c     ,without dissipation in divertor [W/m^2]'
        write(lun,740) 'Lcnnt',tab,
     .   'connection length from midplane to target [m]'
        write(lun,740) 'Lcnnx',tab,
     .   'connection length from x-point to target [m]'
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
c     carlis eb2atm contains the indexes of the B2 atoms mapped to eirene
c            so if in Eirene is [D He Be Ne] and in B2 [D Ne He Be] then
c            eb2atm will be [1 4 2 3]
#ifndef B25_EIRENE      
      do ia=1,natmi
        TEXTA(ia)=adjustl(TEXTAN(ia-1))
      end do
#endif
      do ia=1,natmi !initialize
        eb2atm(ia)=ia
      end do
      do is=1,natmi !species over B2
        do ia=1,natmi !species over Eirene
          if (streql(TEXTA(ia),atoms_list(is))) then
           eb2atm(is)=ia
          end if
        end do
      end do
      print *, 'eb2spcr',eb2spcr
      print *, 'b2ptrgl: list of atoms according to B2: ',atoms_list
      print *, 'list of atoms according to Eirene: ',TEXTA   
      print *, 'ordering of B2 atoms according to Eirene: ',
     ,          (atoms_list(eb2atm(i)),i=1,natmi)
c*** Calculate the particle fluxes

c      print *,'Particle fluxes...' !###
      if(ntns.eq.ntrgts) then !{
        do j=1,ntrgts !{
          do iy=0,ny-1 !{
            do ia=1,natmi
              sja(iy,j,ia)=0.0_R8
              sjp(iy,j,ia)=0.0_R8
              sjt(iy,j,ia)=0.0_R8
            end do
            do imol=1,nmoli
              sjm(iy,j,imol)=0.0_R8
            end do
          end do !}
          ersn_a(j)=0.0_R8
          ersn_p(j)=0.0_R8
          ersn_l(j)=0.0_R8   
          sjs_w(j)=0.0_R8
          sjs_t(j)=0.0_R8
          do ia=1,natmi
            sja_w(j,ia)=0.0_R8
            sjp_w(j,ia)=0.0_R8
            sjt_w(j,ia)=0.0_R8
            sja_t(j,ia)=0.0_R8
            sjp_t(j,ia)=0.0_R8
          end do    
          do imol=1,nmoli
            sjm_w(j,imol)=0.0_R8
            sjm_t(j,imol)=0.0_R8
          end do  
        end do !}

c        do is=0,ns-1
c          print *,'is,am(is)',is,am(is)
c        end do
c        js=ltns(1)+nlim
c        do is=1,natmi
c          print *,'is,wldna(js,is,0)',is,wldna(js,is,0)
c        end do
c        do is=1,ns-nspecies
c          print *,'is,wldpp(js,is,0)',is,wldpp(js,is,0)
c        end do
c        print *,'ismain,ismain0',ismain,ismain0
        u=-1.0_R8
        l=1
        do j=1,ntrgts !{
          k=ltrgpos(j)
c          print *,'j,ltrgpos = ',j,k  !###
          ISS=ltns(j)
          IND=eirdiag_nds_ind(ISS)
          IAS=IND+1
          do iy=lytg,ny-1 !{
            do ia=1,natmi !loop on atoms according to Eirene
              i1=eb2spcr(ia)+1
              i2=eb2spcr(ia)+nfluids(ia)
              do is=i1,i2
                sjp(iy,j,ia)=sjp(iy,j,ia)+u*fna(k+1-l,iy,0,is)
              end do
              if (ft44) then
                sja(iy,j,ia)=wldna_res(ia,iy+IAS)*
     .                       sarea_res(iy+IAS)
                spa(iy,j,ia)=dab2(k+l+1,iy+1,eb2atm(ia),1)*
     *                       tab2(k+l+1,iy+1,eb2atm(ia),1)
              else
                sja(iy,j,ia)=u*fna(k+1-l,iy,0,eb2spcr(ia))
                spa(iy,j,ia)=na(k+l,iy,eb2spcr(ia))*ti(k+l,iy)
              endif
              sjp_t(j,ia)=sjp_t(j,ia)+sjp(iy,j,ia)
              sja_t(j,ia)=sja_t(j,ia)+sja(iy,j,ia)
            end do
            do imol=1,nmoli
              if (ft44) then
                sjm(iy,j,imol)=wldnm_res(imol,iy+IAS)*
     .                       sarea_res(iy+IAS)
                spm(iy,j,imol)=dmb2(k+l+1,iy+1,imol,1)*
     *                       tmb2(k+l+1,iy+1,imol,1)
              else
                sjm(iy,j,imol)=0.0
                spm(iy,j,imol)=0.0
              endif
              sjm_t(j,imol)=sjm_t(j,imol)+sjm(iy,j,imol)
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
            i1=1
            do ia=1,natmi !loop over eirene atoms
              sja_w(j,ia)=wldna(js,ia,0)
              i2=i1+nfluids(ia)-1
              do is=i1,i2
                sjp_w(j,ia)=sjp_w(j,ia)+wldpp(js,ia,0)
              end do
              i1=i2+1
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
            do imol=1,nmoli
              sjm_w(j,imol)=wldnm(js,imol,0)
              if(sjm_t(j,imol).gt.0.0_R8) then !{
                sjm_t(j,imol)=sjm_w(j,imol)/sjm_t(j,imol)
              else !}{
                sjm_t(j,imol)=0.0_R8
              end if !}
              do iy=lytg,ny-1 !{
                sjm(iy,j,imol)=sjm(iy,j,imol)*sjm_t(j,imol)
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
        do_sput=.false.
        if (ntargsp.gt.0.and.ft44) then !sputtering happening
          do j=1,ntrgts
c           Check if sputtering is happening and which material is being sputtered
c           and its atomic mass and density
c carlis It is needed a variable that tells for each target or surface which is the sputtered atom.
c        Temporary solution for this routine is to use the bulk track_species
            do ia=1,n_srfs_m
             if (streql(mtrtrgsg(j),srfs_m(ia))) then
               spt_mat(j)=ia
             end if
            end do
            ISS=ltns(j)
            IND=eirdiag_nds_ind(ISS)
            IAS=IND+1
            sptlis(j)=0
            do iy=lytg,ny-1
              do ia=1,natmi
                if (wldspta_res(ia,iy+IAS).gt.0.0_R8) then
                  sptlis(j)=ia
                  sjs(iy,j)=wldspta_res(ia,iy+IAS)*
     .                    sarea_res(iy+IAS)
                  sjs_t(j)=sjs_t(j)+sjs(iy,j)
                end if
              end do
              do imol=1,nmoli
                if (wldsptm_res(imol,iy+IAS).gt.0.0_R8) then
                  sptlis(j)=imol+natmi
                  sjs(iy,j)=wldsptm_res(imol,iy+IAS)*
     .                    sarea_res(iy+IAS)
                  sjs_t(j)=sjs_t(j)+sjs(iy,j)
                end if
              end do
            end do

            if(ntns.eq.2 .or. ntns.eq.4 .and.
     .         sptlis(j).ne.0) then
              js=ltns(j)+nlim
              if (sptlis(j).le.natmi) then
                sjs_w(j)=wldspta(js,sptlis(j),0)
              else
                sjs_w(j)=wldsptm(js,sptlis(j)-natmi,0)
              end if
              if(sjs_t(j).gt.0.0_R8) then
                sjs_t(j)=sjs_w(j)/sjs_t(j)
              else
                sjs_t(j)=0.0_R8
              end if
            else
              sjs_t(j)=0.0_R8
            end if

            if (sjs_t(j).gt.0.0_R8) then
              do_sput(j)=.true.
              do iy=lytg,ny-1!{
                sjs(iy,j)=sjs_t(j)*sjs(iy,j)
c carlis The erosion of the wall is evaluated using the wall material and not the emitted atom/molecule
c carlis sje [mm/y] = ers_factor[s/y*mol/g*m2/cm2*mm/cm]*srfs_a[g/mol]*sjs[1/s]/srfs_d[g/cm3]/s[m2]
                sje(iy,j)=ers_factor*srfs_a(spt_mat(j))*sjs(iy,j)/
     .                   srfs_d(spt_mat(j))/s(iy,j)
                if (sptlis(j).le.natmi) then
                  sptmass=am(eb2spcr(sptlis(j)))
                else
                  imol = sptlis(j)-natmi
                  sptmass=0.0_R8
                  do ia = 1, natmi
                    if (mlcmp(ia,imol).gt.0) then
                      sptmass=sptmass+mlcmp(ia,imol)*am(eb2spcr(ia))
                    end if
                  end do
                end if
c carlis The emitted amount in g/y is evaluated using the sputtered atom or molecule
c carlis sjs (only in output) [g/y]=sjs[1/s]*srfs_a[g/mol]*nsecyr[s/y]/avogr[1/mol]
                sjs(iy,j)=sjs(iy,j)*sptmass*nsecyr/avogr
                sjs_w(j)=sjs_w(j)+sjs(iy,j)
              end do
              if (sptlis(j).le.natmi) then
                ssmol=texta(sptlis(j))
              else
                ssmol=textmn(sptlis(j)-natmi-1)
                ssmol=adjustl(ssmol)
              end if
              print *,'target material = ', mtrtrgsg(j)
              print *,'sputtered atom from target = ',
     ,                ssmol
              write(*,'(1x,a,1p,16e12.3)') 
     .         'Scaling factors for '//mtrtrgsg(j)//' Sp : ',
     ,          (sjs_t(i),i=1,ntrgts)
            else ! no sputtering
              print *,
     ,'b2ptrgl, zero sputtering flux for target j=',j
              sjs_w(j)=0.0
            end if
          end do !loop on targets
        else
          print *, 'b2ptrgl, no sputtering calculated:'
          print *, '         track_species=0 or no ft44 available'
        end if 
        
        if (do_sput(j)) then
          nclm_er=2+4*natmi+2*ntrgts+2*nmoli
        else
          nclm_er=2+4*natmi+2*nmoli
        endif
        l=1
	do ia=1,natmi
          fpt_str(1,ia)=trim(TEXTA(ia))//'_flux_atm'
          fpt_str(2,ia)=trim(TEXTA(ia))//'_flux_ion'
          fpt_str(3,ia)=trim(TEXTA(ia))//'_temp_atm'
          fpt_str(4,ia)=trim(TEXTA(ia))//'_pres_atm'
        end do
        do k=1,ntrgts
          if (do_sput(k)) then
            if (sptlis(j).le.natmi) then
              spt_str(1,k)=trim(TEXTA(sptlis(j)))//'_sput'
            else
              ssmol=textmn(sptlis(j)-natmi-1)
              ssmol=adjustl(ssmol)
              spt_str(1,k)=trim(ssmol)//'_sput'
            end if
            spt_str(2,k)=trim(mtrtrgsg(k))//'_ersn'
          end if
        end do
        do imol=1,nmoli
          ssmol=textmn(imol-1)
          ssmol=adjustl(ssmol)
          mol_str(1,imol)=trim(ssmol)//'_flux_mol'
          mol_str(2,imol)=trim(ssmol)//'_pres_mol'
        end do
        do j=1,ntrgts !{
          nx_hlp=ltrgpos(j)+l
          jj=j
          if(j.eq.ntrgts) jj=ntrgsx
          open(lun,file=ftrg(jj,2))
          write(fmt_peak,710) nclm_er,tab
          write(lun,740) 'Target erosion - ',trg_label(jj),label
          write(lun,740) 'x',tab,
     .     'coordinate along target [m]. <0 in PFR'
          write(lun,740) 'XX_flux_atm',tab,
     .     'incident flux of neutral atoms of species XX [1/m^2/s]'
          write(lun,740) 'XX_flux_ion',tab,
     .     'incident flux of ions of species XX [1/m^2/s]'
          write(lun,740) 'XX_temp_atm',tab,
     .     'temperature of neutral atoms of species XX [eV]'
          write(lun,740) 'XX_pres_atm',tab,
     .     'pressure of neutral atoms of species XX [Pa]'
          write(lun,740) 'XX_flux_mol',tab,
     .     'incident flux of molecules XX [1/m^2/s]'
          write(lun,740) 'XX_pres_mol',tab,
     .     'pressure of molecules XX [Pa]'
          write(lun,740) 'XX_sput',tab,
     .     'gross sputtering of material XX [g/y]'
          write(lun,740) 'XX_ersn',tab,
     .     'gross erosion of material XX [mm/y]'
          write(lun,740) 'xMP',tab,'radial coordinate in midplane '//
     .     '(IMP for inner target(s), OMP for outer(s)) [m]'
          if (do_sput(j)) then
            write(lun,fmt_peak) 'x',
     ,           (fpt_str(:,ia),ia=1,natmi),
     ,           (mol_str(:,ia),ia=1,nmoli),spt_str(:,j),'xMP'
            write(fmt_peak,720) nclm_er,tab
          else 
            nclm_er=2+4*natmi+2*nmoli
            write(lun,fmt_peak) 'x',
     ,           (fpt_str(:,ia),ia=1,natmi),
     ,           (mol_str(:,ia),ia=1,nmoli),'xMP'
            write(fmt_peak,720) nclm_er,tab
          end if
          if(ntrgts.eq.4) then !{
            k=j
          else !}{
            k=ntrgts
          end if !}
          do iy=lytg,ny-1 !{
            do ia=1,natmi
              fpt_out(1+4*(ia-1))=sja(iy,j,ia)/s(iy,j)
              fpt_out(2+4*(ia-1))=sjp(iy,j,ia)/s(iy,j)
              fpt_out(3+4*(ia-1))=tab2(nx_hlp+1,iy+1,eb2atm(ia),1)/ev
              fpt_out(4+4*(ia-1))=spa(iy,j,ia)
            end do
            do imol=1,nmoli
              mol_out(1+2*(imol-1))=sjm(iy,j,imol)/s(iy,j)
              mol_out(2+2*(imol-1))=spm(iy,j,imol)
            end do
            if (do_sput(j)) then
              spt_out(1)=sjs(iy,j)
              spt_out(2)=sje(iy,j)
              write(lun,fmt_peak) xtrg(iy,j),
     ,         (fpt_out(ia),ia=1,4*natmi),
     ,         (mol_out(ia),ia=1,2*nmoli),
     ,                          spt_out,xmp_hlp(iy,k)
            else 
              write(lun,fmt_peak) xtrg(iy,j),
     ,         (fpt_out(ia),ia=1,4*natmi),
     ,         (mol_out(ia),ia=1,2*nmoli),
     ,                          xmp_hlp(iy,k)
            end if
          end do !}
          close (lun)
          l=1-l
        end do !}
c----------------------------------------------------------------------
c*** Produce a summary file for target sputtering
c
        do ia=1,natmi
          sjp_w(0,ia)=0.0_R8
          sja_w(0,ia)=0.0_R8
          sjt_w(0,ia)=0.0_R8
        end do
        do imol=1,nmoli
          sjm_w(0,imol)=0.0_R8
        end do
        sjs_w(0)=0.0_R8
        do j=1,ntrgts
          do ia=1,natmi
            sjp_w(0,ia)=sjp_w(0,ia)+sjp_w(j,ia)
            sja_w(0,ia)=sja_w(0,ia)+sja_w(j,ia)
            sjt_w(0,ia)=sjt_w(0,ia)+sjt_w(j,ia)
          end do
          do imol=1,nmoli
            sjm_w(0,imol)=sjm_w(0,imol)+sjm_w(j,imol)
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
          do ia=1,natmi
            write(lun,fmt_peak) trim(TEXTA(ia))//' ions on',
     ,                          (sjp_w(i,ia),i=1,ntrgts),sjp_w(0,ia)
            write(lun,fmt_peak) trim(TEXTA(ia))//' atoms on',
     ,                          (sja_w(i,ia),i=1,ntrgts),sja_w(0,ia)
            write(lun,fmt_peak) trim(TEXTA(ia))//' total on',
     ,                          (sjt_w(i,ia),i=1,ntrgts),sjt_w(0,ia)
          end do
          write(lun,*)
          write(lun,*) 'Molecule fluxes [1/s]'
          do imol=1,nmoli
            ssmol=textmn(imol-1)
            ssmol=adjustl(ssmol)
            write(lun,fmt_peak) trim(ssmol)//' molecules on',
     ,                          (sjm_w(i,imol),i=1,ntrgts),sjm_w(0,imol)
          end do
          if (ntargsp.gt.0.and.ft44) then
            do j=1,ntrgts
              sjs_w(0)=sjs_w(0)+sjs_w(j)
            end do
            do j=1,ntrgts
              write(lun,*)
              write(lun,*) 'Total sputtering [g/y]'
              write(lun,fmt_peak) trim(mtrtrgsg(j)),
     ,             (sjs_w(i),i=1,ntrgts),sjs_w(0)
            end do
            write(lun,*)
            write(lun,*) 'Erosion [mm/y]'
            ersn_pm=0.0_R8
            ersn_am=0.0_R8
            ss_m=0.0_R8
            do j=1,ntrgts !{
              if (do_sput(j)) then
                ss(j)=0.0_R8
                ersn_p(j)=0.0_R8
                do iy=lytg,ny-1 !{
                  sn_t=0.0_R8
                  ss(j)=ss(j)+s(iy,j)
                  sn_t=sn_t+sje(i,j)
                  if(sn_t .gt. ersn_p(j)) then !{
                    ersn_p(j)=sn_t
                    ersn_l(j)=xtrg(iy,j)
                  end if !}
                end do !}
                ersn_a(j)=0.0_R8
                ersn_a(j)=ersn_a(j)+sjs_w(j)/ss(j)
                ersn_pm=max(ersn_p(j),ersn_pm)
                ss_m=ss_m+ss(j)
                ersn_am=ersn_am+ersn_a(j)*ss(j)
              else
                ersn_a(j)=0.0_R8     
                ersn_p(j)=0.0_R8  
                ersn_l(j)=0.0_R8       
              end if
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
  720 format('(2x,1p,',i3,'(e12.4,1h',a1,'))')
  740 format('# ',2a,' : ',a)
  750 format('# ',4a)
c=======================================================================
      end subroutine b2ptrgl
