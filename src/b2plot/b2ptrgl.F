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
      use b2mod_b2cmpa
      use b2mod_indirect
      use b2mod_b2cmpa
      use b2mod_b2cmfs
      use b2mod_constants
      use b2mod_external
      use b2mod_b2plot_wall_loading
#ifdef B25_EIRENE
      use eirmod_wneutrals
      use eirmod_ccona
      use eirmod_cinit
#else
      use b2mod_eirbra
#endif
      implicit none
#include <DIMENSIONS.F>
      real (kind=R8) :: BoRiS
      integer lun
      parameter (lun=96)
      integer, save :: use_eirene = 0
#ifdef B25_EIRENE
      integer eb2atm(nspecies)
#else
      real (kind=R8) :: EionH2
      parameter (EionH2=15.4_R8)
      character*8 TEXTA(DEF_NATM)
#endif
c
      integer i,j,jj,k,l,m,n,iss,imol,iy,js,nn,nclm_ld,nclm_er,nx_hlp
      integer i1,i2,i3,ia,ib,ic,is,ix,ik,ind,ims,ias,nclm_na,nclm_nm
      integer kval, kflx
      integer ixstart(ntrgsx,2),ixend(ntrgsx,2)
      integer m_hlp(ntrgsx),npar_hlp(ntrgsx),ntrg_hlp(ntrgsx)
      logical match_found,do_sput(ntrgsx)
c*  s:        area of target segments cut by the grid
c*  wid_par:  ideal parallel power flux at the targets
c*  wid_trg:  ideal power on the targets
c*  xwhm_tt:  exp. decay length on target, calc. from hw & mapped to OMP
      real (kind=R8) ::
     *  h,ha,hb,hc,t,u,a,b,d,p,wo,wt,mion,t1(2),
     ,  z_hlp(ntrgsx),xtg_hlp(-1:ny-1,ntrgsx),xmp_hlp(0:ny-1,ntrgsx),
     ,  xmc_hlp(0:ny-1,ntrgsx),s(0:ny-1,ntrgsx),wi(ntrgsx),
     ,  zctrg(0:ny-1,ntrgsx),rctrg(0:ny-1,ntrgsx),
     ,  xtrg(0:ny-1,ntrgsx),netrg(0:ny-1,ntrgsx),tetrg(0:ny-1,ntrgsx),
     ,  titrg(0:ny-1,ntrgsx),wpttrg(0:ny-1,ntrgsx),
     ,  wpltrg(0:ny-1,ntrgsx),wnutrg(0:ny-1,ntrgsx),
     ,  wrdtrg(0:ny-1,ntrgsx,0:DEF_NATM),wtttrg(0:ny-1,ntrgsx),
     ,  wpk_rd(0:ntrgsx),wpk_te(ntrgsx),wpk_ne(ntrgsx),
!!!     wpk_tt(0:ntrgsx) now in module b2mod_b2plot_all_loading  iyv 11.04.22
     ,  wpk_pl(0:ntrgsx),wpk_nu(0:ntrgsx),wp_xpt(0:ny-1,ntrgsx),
     ,  wlc_tt(ntrgsx),wlc_rd(ntrgsx),wlc_pl(ntrgsx),wlc_nu(ntrgsx),
     ,  wid_par(0:ny-1,ntrgsx),wid_trg(0:ny-1,ntrgsx),
     ,  whtpl(0:ny-1,ntrgsx),wptpl(0:ny-1,ntrgsx),wptnu(0:ny-1,ntrgsx),
     ,  whtnu(0:ny-1,ntrgsx),
c* Integrals of total and radiated power to the targets, power widths
     ,  pwr_tt(ntrgsx),pwr_rd(ntrgsx),fwht_tt(ntrgsx),fwxt_tt(ntrgsx),
     ,  pwr_cr(ntrgsx),frc_cr(ntrgsx),
     ,  hwht_tt(ntrgsx),hwhm_tt(ntrgsx),hwxt_tt(ntrgsx),
     ,  hwxm_tt(ntrgsx),lcnnt(0:ny-1,ntrgsx),lcnnx(0:ny-1,ntrgsx),
     ,  hwxm_idp(ntrgsx),hwxm_idt(ntrgsx),
     ,  hwhm_idp(ntrgsx),hwhm_idt(ntrgsx),
     ,  wpar_hlp(ntrgsx),wtrg_hlp(ntrgsx),
     ,  xwhm_idp(ntrgsx),xwhm_idt(ntrgsx),xwhm_tt(ntrgsx),
cvk kinetic energy of reflected neutrals originating from ions
     ,  wnutpr(0:ny-1,ntrgsx)
c
c*** Arrays for particle fluxes and erosion
c***  sja  impinging flux of atoms
c***  sjp  impinging flux of ions
c***  sjt  total flux of particles to the target
c***  sjm  impinging flux of molecules
c***  sjs  sputtering of target material [g/y]
c***  sje  erosion rate [mm/y]
c***  spa  pressure of atoms along the target
c***    suffix "_t" denotes totals from distributions,
c***    and "_w" the totals from Eirene.
c
c*** Recap of fluxes from Eirene
c***    wldnek  :   net kinergy energy deposited by neutrals. Evaluated as
c***                    EOTAT+EOTML+EOTIO-ERFAAT-ERFMAT-ERFIAT-ERFAML-ERFMML-ERFIML
c***    wldnep  :   energy released by recombination into molecules. Evaluated as
c***                    diss_pot_H2*(PRFAML+PRFPML)
c***    wldna   :   flux of atoms impinging onto the surface. Equal to POTAT
c***    ewlda   :   their average energy. Equal to EOTAT/POTAT
c***    wldnm   :   flux of molecules impinging onto the surface. Equal to POTML
c***    ewldm   :   their average energy. Equal to EOTML/POTML
c***    wldra   :   flux of reflected atoms. Equal to PRFAAT+PRFMAT+PRFIAT
c***    wldrm   :   flux of reflected molecules. Equal to PRFAML+PRFMML+PRFIML
c***    wldpp   :   flux of plasma ions impinging onto the surface. Equal to POTPL
c***    wldpa   :   flux of emitted atoms from BI. Equal to PRFPAT
c***    wldpm   :   flux of emitted molecules from BI. Equal to PRFPML
c***    wldpeb  :   power carried away by these atoms and molecules. Equal to
c***                  ERFPAT+ERFPML+ERFPIO
c***    wldspt  :   flux of sputtered wall material
c***    wldspta :   sputtered flux by emitted atom
c***    wldsptm :   sputtered flux by emitted molecule
c***    wlarea  :   areas of the surface segments from Eirene
c***    wlabsrp :   absorption at the surfaces (1-RECYCT from Eirene)
c***    ewldea_res  :   emitted energy flux of atoms. Equal to
c***                      ERFAAT+ERFMAT+ERFIAT+ERFPHAT
c***    ewldem_res  :   emitted energy flux due to molecules. Equal to
c***                      ERFAML+ERFMML+ERFIML+ERFPHML
c***    ewldrp_res  :   kinetic energy of reflected atoms+molecules originated from ions.
c***                      ERFPAT+ERFPML
c***    ewldmr_res  :   energy due to recombination of atoms and atomic ions into molecules.
c***                     (PRFAML+PRFPML)*H_EDISS

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
     ,   sn_t,spt_out(2),mol_out(2*DEF_NMOL)
      character*11 spt_str(2,ntrgsx),fpt_str(4*DEF_NATM)
      character*15 mol_str(2,DEF_NMOL)
      character*16 ssmol,ssion
      real (kind=R8) ::
     *   ers_factor, nsecyr
      parameter (nsecyr=24.0_R8*3600.0_R8*365.25_R8*0.9925_R8)
      parameter (ers_factor=nsecyr*1.0e-3_R8/avogr)
      real (kind=R8) ::
     *   ss(ntrgsx),ersn_a(ntrgsx),ss_m,ersn_am,ersn_p(ntrgsx),ersn_pm,
     *   ersn_l(ntrgsx)
      logical ft44
      logical streql
      integer, save :: ncall=0   ! iyv 19.07.22
      external streql
      intrinsic sqrt
c*** Data for erosion calculations
c***  n_srfs_m  : number of possible target materials
c***  srfs_m    : list of possible target materials
c***  srfs_d    : specific masses of possible target materials [g/cm**3]
c***  srfs_a    : atomic weight of possible target materials
      integer n_srfs_m,spt_mat(ntrgsx),sptlis(ntrgsx)
      parameter (n_srfs_m=7)
      real (kind=R8) :: srfs_d(n_srfs_m),srfs_a(n_srfs_m)
      real (kind=R8) :: tti, ptext, faext, uaext, taext
      character srfs_m(n_srfs_m)*4
      data srfs_m /'Be','C','CFC','Fe','SS','Mo' ,'W'/
      data srfs_a /  9.0_R8, 12.0_R8, 12.0_R8, 56.0_R8,
     .              56.0_R8, 96.0_R8,184.0_R8/
      data srfs_d /  1.8_R8,  2.0_R8,  2.0_R8,  7.8_R8,
     .               7.8_R8, 10.2_R8, 19.3_R8/
c
c*** Text variables:
c***  ftrg    filenames for the profiles along the targets
c***  fsumm   filenames for the summaries
      integer n_intgrs,n_lrdlcr,n_trgprj
      parameter (n_intgrs=8, n_lrdlcr=2, n_trgprj=7)
      character*16 hlp_frm
      character*256 filename
      character*12 ftrg(ntrgsx,2),fsumm(4),fn_hlp
      character fmt_peak*80,tab*1,tl*1,
     ,  trg_label(ntrgsx)*9,trg_short(ntrgsx)*2,hlp_trg(ntrgsx)*2,
     ,  pkl_label(2)*9,
     ,  hdr_2(n_intgrs)*8,hdr_3(0:n_lrdlcr)*8,hdr_4(n_trgprj)*8
      data ftrg /
     ,  'ld_tg_i.dat', 'ld_tg_iu.dat', 'ld_tg_ou.dat', 'ld_tg_o.dat',
     ,  'fp_tg_i.dat', 'fp_tg_iu.dat', 'fp_tg_ou.dat', 'fp_tg_o.dat'/
      data fsumm /
     ,  'peak_power','sput_trg.sum','peak_power_u','peak_power_w'/
      data trg_label /'inner','inner_top','outer_top','outer'/
      data trg_short / 'il'  ,    'iu'   ,    'ou'   , 'ol'  /
      data pkl_label /'     ','- top'/
      data hdr_2 / ' int_rad',' int_tot','fwht_tot','fwxt_tot',
     ,  'hwht_tot','hwhm_tot','hwxt_tot','hwxm_tot' /
      data hdr_3 / 'core_pwr','core_trg','frcr_trg'/
      data hdr_4 / 'hwxm_par','hwxm_trg','hwhm_par','hwhm_trg',
     ,  'xwhm_par','xwhm_trg','xwhm_tot'/
c======================================================================
c*** If no data on neutral fluxes to the targets is present in the
c*** Eirene output, then a crude estimate is used for the neutrals:
c*** power flux proportional to n*T**(3/2),
c*** particle flux proportional to n*T**(1/2),
c*** sputtering proportional to the D particle flux (chem. dominant),
c*** re-deposition proportional to the C particle flux.
c*** Here n and T are the density and temperature of the hydrogenic
c*** atoms near the target.
c*** The proportionality factor is adjusted to match the integrals
c*** which are calculated in Eirene.
#ifdef DBG
      print *,'Starting b2ptrgl...'
      print *,'ers_factor,lrdlsp= ',ers_factor,lrdlsp
#endif
      tab=char(9)
      if (ncall.eq.0)       ! iyv 19.07.22
     , call ipgeti ('b2mndr_eirene', use_eirene)
      if (use_eirene.ne.0) then
        filename=fort_lc//'44'
        call find_file(filename,ft44)
      else
        ft44 = .false.
      end if

c----------------------------------------------------------------------
c*** Power loading:
c*** Radiation load onto the targets

      call b2ptrdl(wrdtrg)
#ifdef DBG
      print '(2a,3i5,2(1x,l1))','After b2ptrdl: ',
     ,        'ntrgts,lytg,nytg,lrdlsp,lrdlcr = ',   !###
     ,         ntrgts,lytg,nytg,lrdlsp,lrdlcr        !###
#endif
c-----------------------------------------------------------------------
c*** Target geometry

      do j=1,ntrgts !{
        i=itrg(j)
        z_hlp(i)=0.0_R8
        k=ltrgpos(i)
        do iy=0,lseppos(i) !{
          z_hlp(i)=z_hlp(i) - sqrt((crx(k,iy,3)-crx(k,iy,1))**2+
     +                             (cry(k,iy,3)-cry(k,iy,1))**2)
        end do !}
        xtg_hlp(-1,i)=z_hlp(i)
        do iy=0,ny-1 !{
          rctrg(iy,i)=0.5_R8*(crx(k,iy,1)+crx(k,iy,3))
          zctrg(iy,i)=0.5_R8*(cry(k,iy,1)+cry(k,iy,3))
          h=sqrt((crx(k,iy,3)-crx(k,iy,1))**2+
     +           (cry(k,iy,3)-cry(k,iy,1))**2)
          s(iy,i)=gs(k+1,iy,0)
          z_hlp(i)=z_hlp(i)+0.5_R8*h
          xtrg(iy,i)=z_hlp(i)
          z_hlp(i)=z_hlp(i)+0.5_R8*h
          xtg_hlp(iy,i)=z_hlp(i)
        end do !}
      end do !}
#ifdef DBG
      print *
      print *,'ny,xtg_hlp:',ny
      do j=1,ntrgts !{
        k=itrg(j)
        print '(1x,i3,10(1x,f9.6)/(4x,10(1x,f9.6)))',lseppos(k)+1,
     ,                                         (xtg_hlp(i,k),i=-1,ny-1)
      end do !}
#endif

c*** Midplane coordinates

      l=jxi
      m=0
      xmc_hlp=0.0_R8
      do j=1,2*max(1,nncut) !{ Do not use ntrgts here to allow for sign flips
        u=0.0_R8
        if (itrg(j).ne.0) then !{
          k=lpfrpos(itrg(j))
          n=1-m
          do i=0,lseppos(itrg(j)) !{
            h=0.5_R8*sqrt((crx(k,i,n+2)-crx(k,i,n))**2+
     +                    (cry(k,i,n+2)-cry(k,i,n))**2)
            u=u+h
            xmp_hlp(i,itrg(j))=u
            u=u+h
          end do !}
          do i=0,lseppos(itrg(j)) !{
            xmp_hlp(i,itrg(j))=xmp_hlp(i,itrg(j))-u
          end do !}
          u=0.0_R8
          k=l
          xmc_hlp(lseppos(itrg(j)),itrg(j))=0.0_R8
          do i=lseppos(itrg(j))+1,ny-1 !{
            h=0.5_R8*sqrt((crx(k,i,3)-crx(k,i,1))**2+
     +                    (cry(k,i,3)-cry(k,i,1))**2)
            u=u+h
            xmp_hlp(i,itrg(j))=u
            u=u+h
            xmc_hlp(i,itrg(j))=u
          end do !}
        end if !}
        if(j.eq.ntrgts/2) l=jxa
        m=1-m
      end do !}

c!###{
c      print *
c      print *,'ny,xmp_hlp:',ny,jxi,jxa
c      do j=1,ntrgts !{
c        print '(1x,2i3,10f9.5/(4x,10f9.6))',lseppos(itrg(j))+1,
c                                            lpfrpos(itrg(j)),
c     ,                                     (xmp_hlp(i,itrg(j)),i=0,ny-1)
c      end do !}
c!###}
c-----------------------------------------------------------------------
c*** Projected (ideal) power flux on targets

      wid_par=0.0_R8
      wid_trg=0.0_R8
      l=-1
      u=-1.0_R8
      do j=1,2*max(1,nncut) !{ Do not use ntrgts here to allow for sign flips
        if (itrg(j).gt.0) then !{
          m=ltrgpos(itrg(j))
          k=lpfrpos(itrg(j))-l
          i=lseppos(itrg(j))+1
          wo=0.0_R8
          do iy=i,ny-1 !{
            wo=wo+fht(k,iy,0)-fhj(k,iy,0)
          end do !}
          wo=wo*u
          if (wo.le.0.0_R8) then !{
            write(*,'(a,3i5,1p,2(1x,1e12.5))')
     .       'Wrong sign of qpar. j,k,i,u,wo=',j,k,i,u,wo
          else !}{
            n=k
            nn=lpfrpos(itrg(j-l))+l
            do ix=k-l,nn,-l !{
              wt=wo
              wo=0.0_R8
              do iy=i,ny-1 !{
                wo=wo+fht(ix,iy,0)-fhj(ix,iy,0)
              end do !}
              wo=wo*u
              if(wo.lt.wt) exit
              n=ix
            end do !}
c!###{
c          print *,'  j  nn   n   m   k   wo          wt:'
c          print '(5i4,1p,2e12.5)',j,nn,n,m,k,wo,wt
c          print '(a4,a9,5a12)','iy','fe+','sx ','pbs','Rk ','Rm ','Stg'
c!###}
            k=n
            do iy=i,ny-1 !{
              wo=fht(k,iy,0)-fhj(k,iy,0)
c!###{
c            print '(1p,i4,6e12.4)',iy,wo,
c     ,        sx(k,iy),pbs(k,iy,0),cr(k,iy),cr(m,iy),s(iy,itrg(j))
c!###}
              wp_xpt(iy,itrg(j))=wo*u/abs(pbs(k,iy,0))
              wid_par(iy,itrg(j))=wp_xpt(iy,itrg(j))*cr(k,iy)/cr(m,iy)
              wid_trg(iy,itrg(j))=wo*u/s(iy,itrg(j))
            end do !}
          end if !}
        end if !}
        l=-l
        u=-u
      end do !}

c-----------------------------------------------------------------------
c*** Power loading:

      ixstart(1,1)=jxi
      ixend(1,1)=0
      ixstart(1,2)=leftcut(1)-1
      ixend(1,2)=0
      if (ntrgts.gt.0) then
        ixstart(ntrgts,1)=jxa
        ixend(ntrgts,1)=nx-1
        ixstart(ntrgts,2)=rightcut(1)
        ixend(ntrgts,2)=nx-1
      end if
      if (ntrgts.gt.2) then
        ixstart(2,1)=jxi
        ixend(2,1)=nxtl-1
        ixstart(3,1)=jxa
        ixend(3,1)=nxtr+1
        ixstart(2,2)=leftcut(2)
        ixend(2,2)=nxtl-1
        ixstart(3,2)=rightcut(2)-1
        ixend(3,2)=nxtr+1
      endif
c
c*** Plasma and neutrals
c
      l=1
      u=-1.0_R8
      do j=1,2*max(1,nncut) !{ Do not use ntrgts here to allow for sign flips
        if (itrg(j).gt.0) then !{
          k=ltrgpos(itrg(j))
          kval=k+l
          kflx=k+1-l
          if( .not. from_b2plot) kflx=k+1
          do iy=lytg,ny-1 !{
            wo=0.0_R8
            wnutrg(iy,itrg(j))=0.0_R8
            whtpl(iy,itrg(j))=fhe(kflx,iy,0)+fhi(kflx,iy,0)+
     .                    fhi_ext(kflx,iy,0)
            do is=0,ns-1 !{
              tti=ti(kval,iy)
              if( .not. from_b2plot )
     *            tti=0.5_R8*(ti(kflx,iy)+ti(k,iy))
              if (.not.is_neutral(is).or.ft44) then
                wo=wo+fhp(kflx,iy,0,is)
c carlis 31jan2017  added ti as suggested by sw in b2blne.F
              else
                wnutrg(iy,itrg(j))=wnutrg(iy,itrg(j))+(1.0_R8-BoRiS)*
     *           (fhm(kflx,iy,0,is)+fna(kflx,iy,0,is)*tti)
              endif
              whtpl(iy,itrg(j))=whtpl(iy,itrg(j))+(1.0_R8-BoRiS)*
     *         (fhm(kflx,iy,0,is)+fna(kflx,iy,0,is)*tti)
            end do !}
            do is=0,ns_ext-1 !{
              ptext=pt_ext(kval,iy,is)
              uaext=ua_ext(kval,iy,is)
              faext=fa_ext(kflx,iy,0,is)
              if( .not. from_b2plot) then
                ptext=0.5_R8*(pt_ext(kflx,iy,is)+pt_ext(k,iy,is))
                uaext=0.5_R8*(ua_ext(kflx,iy,is)+ua_ext(k,iy,is))
                faext=0.5_R8*(fa_ext(kflx,iy,0,is)+fa_ext(k,iy,0,is))
                taext=0.5_R8*(ta_ext(kflx,iy,is)+ta_ext(k,iy,is))
              end if
              if (.not.is_ext_neutral(is).or.ft44) then
                wo=wo+faext*ptext*ev
c carlis 31jan2017  added ta as suggested by sw in b2blne.F
              else
                wnutrg(iy,itrg(j))=wnutrg(iy,itrg(j))+
     .                 faext
     .               *(1.0_R8-BoRiS)*(0.5_R8*am_ext(is)*mp*
     .                 uaext**2+taext)
              end if
              whtpl(iy,itrg(j))=whtpl(iy,itrg(j))+
     .               faext
     .             *(1.0_R8-BoRiS)*(taext+0.5_R8*
     .               am_ext(is)*mp*uaext**2)
            end do !}
            whtpl(iy,itrg(j))=whtpl(iy,itrg(j))+fne(kflx,iy,0)
     .           *(1.0_R8-BoRiS)*te(kval,iy)
c carlis 31jan2017  added electron contribution as suggested by sw in b2blne.
c           wo=-1.0e-6_R8*wo/s(iy,itrg(j))
            wptpl(iy,itrg(j))=wo/s(iy,itrg(j))*u
            whtpl(iy,itrg(j))=whtpl(iy,itrg(j))/s(iy,itrg(j))*u
            wpltrg(iy,itrg(j))=wptpl(iy,itrg(j))+whtpl(iy,itrg(j))
            if(.not.ft44) wnutrg(iy,itrg(j))=
     .                    wnutrg(iy,itrg(j))*u/s(iy,itrg(j))
            netrg(iy,itrg(j))=ne(kval,iy)
            tetrg(iy,itrg(j))=te(kval,iy)/ev
            titrg(iy,itrg(j))=ti(kval,iy)/ev
            if (ft44) then
              do iss=1,nsts !{
                if (eirdiag_nds_typ(iss).ne.2) cycle
                if (eirdiag_nds_srf(iss).ne.ltrgpos(itrg(j))+1) cycle
                if (eirdiag_nds_start(iss).gt.(iy+1)) cycle
                if (eirdiag_nds_end(iss).lt.(iy+1)) cycle
                IND=eirdiag_nds_ind(ISS)
                IAS=IND+1-(eirdiag_nds_start(iss)-1)
                do is=1,natmi
c             carlis Using values directly calculated from Eirene
                  wnutrg(iy,itrg(j))=wnutrg(iy,itrg(j))+
     .                                   (ewlda_res(is,iy+IAS)
     .                                   -ewldea_res(is,iy+IAS))*
     .                                    sarea_res(iy+IAS)
                enddo
                do is=1,nmoli
c             carlis Using values directly calculated from Eirene
                  wnutrg(iy,itrg(j))=wnutrg(iy,itrg(j))+
     .                                   (ewldm_res(is,iy+IAS)
     .                                   -ewldem_res(is,iy+IAS)
     .                                   +ewldmr_res(is,iy+IAS))*
     .                                    sarea_res(iy+IAS)
                enddo
c             kinetic energy of recycled and sputtered neutrals
                wnutpr(iy,itrg(j))=ewldrp_res(iy+IAS)
              enddo !}
              whtpl(iy,itrg(j))=whtpl(iy,itrg(j))-wnutpr(iy,itrg(j))
              wpltrg(iy,itrg(j))=wpltrg(iy,itrg(j))-wnutpr(iy,itrg(j))
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
                      mion=mion+mlcmp(ims,imol)*am(eb2atcr(ims))
                    enddo
                  endif
                  imol=imol+1
                enddo
                t=tib2(kval+1,iy+1,is,1)/ev
                if (mion.gt.0.0_R8)
     >           wnutrg(iy,itrg(j))=wnutrg(iy,itrg(j))+
     .            s(iy,itrg(j))*dib2(kval+1,iy+1,is,1)*
     *            sqrt(max(t/mion,0.0_R8))*0.25_R8*1.0e4_R8*
     .           (2.0_R8*tib2(kval+1,iy+1,is,1)+EionH2*ev)
              enddo
            endif
c carlis  21 feb 2017 added calculation of connection length from midplane to targets
c         and from x-point to targets
            do ik=1,2
              ix=ixstart(itrg(j),ik)
              t1(ik)=0.0_R8
              if (ix.eq.-2) cycle
              if (ik.eq.1.and.iy.gt.jsep) then
c carlis If starting @midplane in SOL then remove half of starting cell
                t1(ik)=t1(ik)-hx(ix,iy)/2.0_R8/
     1             (sqrt(bb(ix,iy,0)**2+bb(ix,iy,1)**2)/bb(ix,iy,3))
              end if
              if (l.eq.1) then
                do ix=ixend(itrg(j),ik),ixstart(itrg(j),ik)
                  t1(ik)=t1(ik)+hx(ix,iy)/
     1                (sqrt(bb(ix,iy,0)**2+bb(ix,iy,1)**2)/bb(ix,iy,3))
                end do
              else
                do ix=ixstart(itrg(j),ik),ixend(itrg(j),ik)
                  t1(ik)=t1(ik)+hx(ix,iy)/
     1                (sqrt(bb(ix,iy,0)**2+bb(ix,iy,1)**2)/bb(ix,iy,3))
                end do
              endif
            end do
            lcnnt(iy,itrg(j))=t1(1)
            lcnnx(iy,itrg(j))=t1(2)
          end do !}
        end if !}
        l=1-l
        u=-u
      end do !}
      if (ntrgts.ge.2) then
        do iy=lytg,jsep !{
c carlis  If in PFR region then connection length from midplane to target will
c         be equal to connection length from target to target
          lcnnt(iy,itrg(1))=lcnnx(iy,itrg(1))+lcnnx(iy,itrg(2))
          lcnnt(iy,itrg(ntrgts))=lcnnx(iy,itrg(ntrgts))+
     .                           lcnnx(iy,itrg(ntrgts-1))
          if (ntrgts.gt.2) then
            lcnnt(iy,itrg(2))=lcnnt(iy,itrg(1))
            lcnnt(iy,itrg(3))=lcnnt(iy,itrg(ntrgts))
          end if
        end do
      end if

c
c*** Re-scale the neutral target loading and calculate the sums
c
      do j=1,ntrgts !{
        wi(itrg(j))=0.0_R8
        do i=lytg,ny-1 !{
          wi(itrg(j))=wi(itrg(j))+wnutrg(i,itrg(j))
        end do !}
      end do !}
c      print *,'ntrgts, wi: ',ntrgts,(wi(itrg(i)),i=1,ntrgts)  !####
      do j=1,ntrgts !{
        if(itrg(j).eq.0) cycle
        if (ft44) then
          if(count(ktns(:)==itrg(j)).gt.0) then !{
            wo=0.0_R8
            do iss=1,ntns !{
              if (ktns(iss).ne.itrg(j)) cycle
              js=ltns(iss)+nlim
              wo=wo+wldnek(js,0)+wldnep(js,0)
c             print *,'j,js,wldnek,wldnep = ',
c     ,          itrg(j),js,wldnek(js,0),wldnep(js,0)
            end do !}
            if(wi(itrg(j)) .gt. 0.0_R8) then !{
              wi(itrg(j))=wo/wi(itrg(j))
            else !}{
              wi(itrg(j))=0.0_R8
            end if !}
          else !}{
            write(*,'(a,i2,a,i2,a)')
     ,               'b2ptrgl: no Eirene surface for target ',j,
     ,               ' - no neutral load onto the target calculated'
            write(*,'(2a)') '  ==> check specification of ltns in ',
     ,                               'b2.neutrals.parameters'
            wi(itrg(j))=0.0_R8
          end if !}
        endif
        do iy=lytg,ny-1 !{
c          if (ft44) wnutrg(iy,itrg(j))=
c     ,     1.0e-6_R8*wnutrg(iy,itrg(j))*wi(itrg(j))/s(iy,itrg(j))
          if (ft44) then
            wnutrg(iy,itrg(j))=
     .       wnutrg(iy,itrg(j))*wi(itrg(j))/s(iy,itrg(j))
            if(abs(wo).gt.0.0_R8) then !{
              wptnu(iy,itrg(j))=wnutrg(iy,itrg(j))*wldnep(js,0)/wo
              whtnu(iy,itrg(j))=wnutrg(iy,itrg(j))-wptnu(iy,itrg(j))
            else !}{
              wptnu(iy,itrg(j))=0.0_R8
              whtnu(iy,itrg(j))=0.0_R8
            end if !}
          else !}{
            wptnu(iy,itrg(j))=0.0_R8
            whtnu(iy,itrg(j))=0.0_R8
          endif
          wpttrg(iy,itrg(j))=wnutrg(iy,itrg(j))+wpltrg(iy,itrg(j))
          wtttrg(iy,itrg(j))=wpttrg(iy,itrg(j))+wrdtrg(iy,itrg(j),0)
        end do !}
      end do !}
      write(*,'(1x,a,1p,16e12.3)')
     -           'Scaling factors for Wn : ',(wi(itrg(i)),i=1,ntrgts)

      nclm_ld=21
#ifdef B25_EIRENE
      nclm_na=natmi
      nclm_nm=nmoli
#else
      nclm_na=nspecies
      nclm_nm=0
#endif
      if(lrdlsp) nclm_ld=nclm_ld+nclm_na
      if(lrdlcr) nclm_ld=nclm_ld+1

      if (from_b2plot) then    ! iyv 11.04.22

      do j=1,ntrgts !{
        jj=itrg(j)
        if(j.eq.ntrgts.and..not.(ntrgts.eq.1.and.itrg(1).eq.1))
     >   jj=ntrgsx
        if(.not.LSN) jj=ntrgsx+1-jj
        open(lun,file=ftrg(jj,1))
        write(fmt_peak,710) nclm_ld
        write(lun,740) 'Target loading - ',trg_label(jj),trim(label)
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
        write(lun,740) 'Wheat',tab,
     .   'power flux due to particles impinging on wall [W/m^2]'
        write(lun,740) 'Wpot',tab,
     .   'power flux due to recombination [W/m^2]'
        write(lun,740) 'Whtpl',tab,
     .   'power flux due to plasma particles impinging on wall [W/m^2]'
        write(lun,740) 'Wptpl',tab,
     .   'power flux due to plasma recombination [W/m^2]'
        write(lun,740) 'xMP',tab,'radial coordinate in midplane '//
     .   '(IMP for inner target(s), OMP for outer(s)) [m]'
        write(lun,740) 'Rclfc',tab,'radial coordinate on target [m]'
        write(lun,740) 'Zclfc',tab,'vertical coordinate on target [m]'
        write(lun,740) 'Wpar_xpt',tab,'parallel power flux around',
     ,   ' x-point [W/m^2]'
        write(lun,740) 'WWpar',tab,'would-be parallel power flux',
     ,   ' without dissipation in divertor [W/m^2]'
        write(lun,740) 'WWtrg',tab,'would-be target power flux',
     ,   ' without dissipation in divertor [W/m^2]'
        write(lun,740) 'Lcnnt',tab,
     .   'connection length from midplane to target [m]'
        write(lun,740) 'Lcnnx',tab,
     .   'connection length from x-point to target [m]'
        if(lrdlsp) then !{
          write(lun,740) 'rad_XX',tab,
     .     'radiation power flux from species XX [W/m^2]'
          write(lun,fmt_peak) '#','x','ne','Te','Ti','Wtot','Wpart',
     ,            'Wrad','Wpls','Wneut','Wheat','Wpot','Whtpl',
     ,            'Wptpl','xMP','Rclfc','Zclfc','Wpar_xpt',
     ,            'WWpar','WWtrg','Lcnnt','Lcnnx',
     ,            ('rad_'//textan(is)(2:3),is=0,nclm_na-1)
        else if(lrdlcr) then !}{
          write(lun,740) 'rd_cr',tab,
     ,     'radiation power flux from core [W/m^2]'
          write(lun,fmt_peak) '#','x','ne','Te','Ti','Wtot','Wpart',
     ,            'Wrad','Wpls','Wneut','Wheat','Wpot','Whtpl',
     ,            'Wptpl','xMP','Rclfc','Zclfc','Wpar_xpt',
     ,            'WWpar','WWtrg','Lcnnt','Lcnnx','rd_cr'
        else !}{
          write(lun,fmt_peak) '#','x','ne','Te','Ti','Wtot','Wpart',
     ,            'Wrad','Wpls','Wneut','Wheat','Wpot','Whtpl',
     ,            'Wptpl','xMP','Rclfc','Zclfc','Wpar_xpt',
     ,            'WWpar','WWtrg','Lcnnt','Lcnnx'
        end if !}
        write(fmt_peak,720) nclm_ld
        if(ntrgts.eq.4) then !{
          k=itrg(j)
        else !}{
          k=itrg(ntrgts)
        end if !}
        do i=lytg,ny-1 !{
          if(lrdlsp) then !{
            write(lun,fmt_peak) xtrg(i,itrg(j)),
     ,              netrg(i,itrg(j)),tetrg(i,itrg(j)),titrg(i,itrg(j)),
     ,              wtttrg(i,itrg(j)),wpttrg(i,itrg(j)),
     ,              wrdtrg(i,itrg(j),0),wpltrg(i,itrg(j)),
     ,              wnutrg(i,itrg(j)),whtnu(i,itrg(j))+whtpl(i,itrg(j)),
     ,              wptnu(i,itrg(j))+wptpl(i,itrg(j)),whtpl(i,itrg(j)),
     ,              wptpl(i,itrg(j)),
     ,              xmp_hlp(i,k),rctrg(i,itrg(j)),zctrg(i,itrg(j)),
     ,              wp_xpt(i,itrg(j)),
     ,              wid_par(i,itrg(j)),wid_trg(i,itrg(j)),
     ,              lcnnt(i,itrg(j)),lcnnx(i,itrg(j)),
     ,             (wrdtrg(i,itrg(j),is),is=1,nclm_na)
          else if(lrdlcr) then !}{
            is=1
            write(lun,fmt_peak) xtrg(i,itrg(j)),
     ,              netrg(i,itrg(j)),tetrg(i,itrg(j)),titrg(i,itrg(j)),
     ,              wtttrg(i,itrg(j)),wpttrg(i,itrg(j)),
     ,              wrdtrg(i,itrg(j),0),wpltrg(i,itrg(j)),
     ,              wnutrg(i,itrg(j)),whtnu(i,itrg(j))+whtpl(i,itrg(j)),
     ,              wptnu(i,itrg(j))+wptpl(i,itrg(j)),
     ,              whtpl(i,itrg(j)),wptpl(i,itrg(j)),
     ,              xmp_hlp(i,k),rctrg(i,itrg(j)),zctrg(i,itrg(j)),
     ,              wp_xpt(i,itrg(j)),
     ,              wid_par(i,itrg(j)),wid_trg(i,itrg(j)),
     ,              lcnnt(i,itrg(j)),lcnnx(i,itrg(j)),
     ,              wrdtrg(i,itrg(j),is)
          else !}{
            write(lun,fmt_peak) xtrg(i,itrg(j)),
     ,              netrg(i,itrg(j)),tetrg(i,itrg(j)),titrg(i,itrg(j)),
     ,              wtttrg(i,itrg(j)),wpttrg(i,itrg(j)),
     ,              wrdtrg(i,itrg(j),0),wpltrg(i,itrg(j)),
     ,              wnutrg(i,itrg(j)),whtnu(i,itrg(j))+whtpl(i,itrg(j)),
     ,              wptnu(i,itrg(j))+wptpl(i,itrg(j)),whtpl(i,itrg(j)),
     ,              wptpl(i,itrg(j)),
     ,              xmp_hlp(i,k),rctrg(i,itrg(j)),zctrg(i,itrg(j)),
     ,              wp_xpt(i,itrg(j)),
     ,              wid_par(i,itrg(j)),wid_trg(i,itrg(j)),
     ,              lcnnt(i,itrg(j)),lcnnx(i,itrg(j))
          end if !}
        end do !}
        close (lun)
      end do !}
      end if                 ! iyv 11.04.22
c----------------------------------------------------------------------
c*** Calculate the peak power loads onto the targets and write them
c*** down to the summary file
c
      l=1
      wpk_tt(0)=0
      wpk_rd(0)=0
      wpk_pl(0)=0
      wpk_nu(0)=0
      do j=1,2*max(nncut,1) !{ Do not use ntrgts here to allow for sign flips
        if (itrg(j).gt.0) then !{
          k=ltrgpos(itrg(j))
          wpk_tt(itrg(j))=wtttrg(lytg,itrg(j))
          wpk_rd(itrg(j))=wrdtrg(lytg,itrg(j),0)
          wpk_pl(itrg(j))=wpltrg(lytg,itrg(j))
          wpk_nu(itrg(j))=wnutrg(lytg,itrg(j))
          wpk_te(itrg(j))=te(k+l,lytg)/ev
          wpk_ne(itrg(j))=ne(k+l,lytg)
          wlc_tt(itrg(j))=xtrg(lytg,itrg(j))
          wlc_rd(itrg(j))=xtrg(lytg,itrg(j))
          wlc_pl(itrg(j))=xtrg(lytg,itrg(j))
          wlc_nu(itrg(j))=xtrg(lytg,itrg(j))
          m_hlp(itrg(j))=lytg
          do i=lytg+1,ny-1 !{
            if(wtttrg(i,itrg(j)) .gt. wpk_tt(itrg(j))) then !{
              wpk_tt(itrg(j))=wtttrg(i,itrg(j))
              wpk_te(itrg(j))=te(k+l,i)/ev
              wpk_ne(itrg(j))=ne(k+l,i)
              wlc_tt(itrg(j))=xtrg(i,itrg(j))
              m_hlp(itrg(j))=i
            end if !}

            if(from_b2plot) then           !iyv 11.04.22

            if(wrdtrg(i,itrg(j),0) .gt. wpk_rd(itrg(j))) then !{
              wpk_rd(itrg(j))=wrdtrg(i,itrg(j),0)
              wlc_rd(itrg(j))=xtrg(i,itrg(j))
            end if !}
            if(wpltrg(i,itrg(j)) .gt. wpk_pl(itrg(j))) then !{
              wpk_pl(itrg(j))=wpltrg(i,itrg(j))
              wlc_pl(itrg(j))=xtrg(i,itrg(j))
            end if !}
            if(wnutrg(i,itrg(j)) .gt. wpk_nu(itrg(j))) then !{
              wpk_nu(itrg(j))=wnutrg(i,itrg(j))
              wlc_nu(itrg(j))=xtrg(i,itrg(j))
            end if !}
            end if ! if(from_b2plot)   !iyv 11.04.22
          end do !}
        end if !}
        l=1-l
      end do !}


      if ( .not. from_b2plot) then
        ncall = ncall + 1
        return       !iyv 11.04.22
      endif

      write(fmt_peak,700) 24
  700 format('(1p,e12.5,',i3,'(a1,e12.5))')
      i1=1
      i2=ntrgts
      fn_hlp=fsumm(1)
      tl=' '
      do j=1,2 !{
        if((j.eq.1 .or. ntrgts.eq.4).and.ntrgts.gt.0) then !{
          open(lun,file=fn_hlp)
          write(lun,740) 'Peak powers ',pkl_label(j),trim(label)
          if(.not.LSN) then
            i3=i2
            i2=i1
            i1=i3
          endif
          i3=i2-i1
          if (i3.eq.0 .and. ntrgts.eq.1) i3=1

          wpk_tt(0)=max(wpk_tt(i1),wpk_tt(i2))
          wpk_rd(0)=max(wpk_rd(i1),wpk_rd(i2))
          wpk_pl(0)=max(wpk_pl(i1),wpk_pl(i2))
          wpk_nu(0)=max(wpk_nu(i1),wpk_nu(i2))
          if (i1.ne.i2) then
           write(lun,'(a6,a1,24(a1,2a))') '# time',tab,
     -      tab,'peak_t_i',tl,tab,'peak_t_o',tl,tab,'peak_t_x',tl,
     ,      tab, 'loc_t_i',tl,tab, 'loc_t_o',tl,
     ,      tab,'peak_r_i',tl,tab,'peak_r_o',tl,tab,'peak_r_x',tl,
     ,      tab, 'loc_r_i',tl,tab, 'loc_r_o',tl,
     ,      tab,'peak_p_i',tl,tab,'peak_p_o',tl,tab,'peak_p_x',tl,
     ,      tab, 'loc_p_i',tl,tab, 'loc_p_o',tl,
     ,      tab,'peak_n_i',tl,tab,'peak_n_o',tl,tab,'peak_n_x',tl,
     ,      tab, 'loc_n_i',tl,tab, 'loc_n_o',tl,
     ,      tab,'peak_tei',tl,tab,'peak_teo',tl,
     ,      tab,'peak_nei',tl,tab,'peak_neo',tl
          else
           write(lun,'(a6,a1,14(a1,2a))') '# time',tab,
     -      tab,'peak_t_i',tl,tab,'peak_t_x',tl,
     ,      tab, 'loc_t_i',tl,
     ,      tab,'peak_r_i',tl,tab,'peak_r_x',tl,
     ,      tab, 'loc_r_i',tl,
     ,      tab,'peak_p_i',tl,tab,'peak_p_x',tl,
     ,      tab, 'loc_p_i',tl,
     ,      tab,'peak_n_i',tl,tab,'peak_n_x',tl,
     ,      tab, 'loc_n_i',tl,
     ,      tab,'peak_tei',tl,tab,'peak_nei',tl
          end if
          write(lun,fmt_peak) tim,
     ,     (tab, wpk_tt(i),i=i1,i2,i3), tab, wpk_tt(0),
     ,     (tab, wlc_tt(i),i=i1,i2,i3),
     ,     (tab, wpk_rd(i),i=i1,i2,i3), tab, wpk_rd(0),
     ,     (tab, wlc_rd(i),i=i1,i2,i3),
     ,     (tab, wpk_pl(i),i=i1,i2,i3), tab, wpk_pl(0),
     ,     (tab, wlc_pl(i),i=i1,i2,i3),
     ,     (tab, wpk_nu(i),i=i1,i2,i3), tab, wpk_nu(0),
     ,     (tab, wlc_nu(i),i=i1,i2,i3),
     ,     (tab, wpk_te(i),i=i1,i2,i3),
     ,     (tab, wpk_ne(i),i=i1,i2,i3)
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

      wpar_hlp=maxval(wid_par,1)
      npar_hlp=maxloc(wid_par,1)-1
      wtrg_hlp=maxval(wid_trg,1)
      ntrg_hlp=maxloc(wid_trg,1)-1
c!###{
c      print *,'wid_par:'
c      do j=1,ntrgts !{
c        print '(1p,10e12.4)',(wid_par(i,itrg(j)),i=0,ny-1)
c      end do !}
c      print *,'wid_trg:'
c      do j=1,ntrgts !{
c        print '(1p,10e12.4)',(wid_trg(i,itrg(j)),i=0,ny-1)
c      end do !}
c      print '(a,1p,8e12.4)','wpar_hlp: ',(wpar_hlp(itrg(j)),j=1,ntrgts)
c      print '(a,i8,7i12)','npar_hlp: ',(npar_hlp(itrg(j)),j=1,ntrgts)
c      print '(a,1p,8e12.4)','wtrg_hlp: ',(wtrg_hlp(itrg(j)),j=1,ntrgts)
c      print '(a,i8,7i12)','ntrg_hlp: ',(ntrg_hlp(itrg(j)),j=1,ntrgts)
c!###}
      ia=-2
      ib=-2
      ic=-2
      ha=0.0_R8
      hb=0.0_R8
      hc=0.0_R8
      do j=1,ntrgts !{
        hlp_trg(itrg(j))=trg_short(itrg(j))
        pwr_tt(itrg(j))=0.0_R8
        pwr_rd(itrg(j))=0.0_R8
        pwr_cr(itrg(j))=0.0_R8
        frc_cr(itrg(j))=0.0_R8
        fwht_tt(itrg(j))=0.0_R8
        fwxt_tt(itrg(j))=0.0_R8
        hwht_tt(itrg(j))=0.0_R8
        hwhm_tt(itrg(j))=0.0_R8
        hwxt_tt(itrg(j))=0.0_R8
        hwxm_tt(itrg(j))=0.0_R8
        do i=lytg+1,ny-1 !{
          pwr_tt(itrg(j))=pwr_tt(itrg(j))+wtttrg(i,itrg(j))*s(i,itrg(j))
          pwr_rd(itrg(j))=pwr_rd(itrg(j))+
     .                                  wrdtrg(i,itrg(j),0)*s(i,itrg(j))
          if(lrdlcr) pwr_cr(itrg(j))=pwr_cr(itrg(j))+
     .                                  wrdtrg(i,itrg(j),1)*s(i,itrg(j))
        end do !}
        if(add_rad.gt.0.0_R8)
     >   frc_cr(itrg(j))=1.0e-6_R8*pwr_cr(itrg(j))/add_rad
        p=0.5_R8*wpk_tt(itrg(j))
        a=xtg_hlp(-1,itrg(j))
        if(wtttrg(0,itrg(j)).lt.p) then !{
          do i=m_hlp(itrg(j))-1,0,-1 !{
            if(wtttrg(i,itrg(j)).lt.p) then !{
              h=(wtttrg(i+1,itrg(j))-wtttrg(i,itrg(j)))
              if(h.gt.0.0_R8) then !{
                h=(p-wtttrg(i,itrg(j)))/h
              else !}{
                h=1.0_R8
              end if !}
              a=xtg_hlp(i,itrg(j))*(1.0_R8-h)+xtg_hlp(i+1,itrg(j))*h
              ia=i  !###
              ha=h  !###
              exit
            end if !}
          end do !}
        end if !}
        b=xtg_hlp(ny-1,itrg(j))
        d=xmc_hlp(ny-1,itrg(j))
        if(wtttrg(ny-1,itrg(j)).lt.p) then !{
          do i=m_hlp(itrg(j))+1,ny-1 !{
            if(wtttrg(i,itrg(j)).lt.p) then !{
              h=(wtttrg(i-1,itrg(j))-wtttrg(i,itrg(j)))
              if(h.gt.0.0_R8) then !{
                h=(p-wtttrg(i,itrg(j)))/h
              else !}{
                h=1.0_R8
              end if !}
              b=xtg_hlp(i+1,itrg(j))*(1.0_R8-h)+xtg_hlp(i,itrg(j))*h
              d=xmc_hlp(i+1,itrg(j))*(1.0_R8-h)+xmc_hlp(i,itrg(j))*h
              ib=i  !###
              hb=h  !###
              exit
            end if !}
          end do !}
        end if !}
        fwht_tt(itrg(j))=b-a
        hwht_tt(itrg(j))=b-xtrg(m_hlp(itrg(j)),itrg(j))
        hwhm_tt(itrg(j))=d-xmp_hlp(m_hlp(itrg(j)),itrg(j))
        xwhm_tt(itrg(j))=hwhm_tt(itrg(j))/log(2.0_R8)
c*** Total power outside the strike point position and power half-width
        wt=0.0_R8
        do i=lseppos(itrg(j))+1,ny-1 !{
          wt=wt+wtttrg(i,itrg(j))*s(i,itrg(j))
        end do !}
        fwxt_tt(itrg(j))=pwr_tt(itrg(j))/
     .   (wpk_tt(itrg(j))*2.0_R8*pi*rctrg(m_hlp(itrg(j)),itrg(j)))
        hwxt_tt(itrg(j))=wt/
     .   (wpk_tt(itrg(j))*2.0_R8*pi*rctrg(m_hlp(itrg(j)),itrg(j)))
        hwxm_tt(itrg(j))=xmc_hlp(ny-1,itrg(j))
        wo=0.0_R8
        do i=lseppos(itrg(j))+1,ny-1 !{
          if(xtg_hlp(i,itrg(j)).ge.hwxt_tt(itrg(j))) then !{
            h=(hwxt_tt(itrg(j))-xtg_hlp(i-1,itrg(j)))/
     .       (xtg_hlp(i,itrg(j))-xtg_hlp(i-1,itrg(j)))
            hwxm_tt(itrg(j))=
     .        xmc_hlp(i-1,itrg(j))*(1.0_R8-h)+xmc_hlp(i,itrg(j))*h
            ic=i  !###
            hc=h  !###
            exit
          end if !}
c*          h=xtg_hlp(i,itrg(j))
c*          wo=xmc_hlp(i,itrg(j))
        end do !}
        wo=wpar_hlp(itrg(j))
        m=npar_hlp(itrg(j))
        p=wtrg_hlp(itrg(j))
        n=ntrg_hlp(itrg(j))
c        print '(1x,a,2i4,1p,2e12.3)','m,n,wo,p=',m,n,wo,p !###
        wt=0.0_R8
        t=0.0_R8
        do i=lseppos(itrg(j))+1,ny-1 !{
          a=abs(xmc_hlp(i,itrg(j))-xmc_hlp(i-1,itrg(j)))
          wt=wt+wid_par(i,itrg(j))*a
          t=t+wid_trg(i,itrg(j))*a
        end do !}
c        print '(a,i4,1p,4e12.4)','j,wt,wo,t,p: ',j,wt,wo,t,p  !###
        if(wo.gt.0.0_R8) hwxm_idp(itrg(j))=wt/wo
        if(p.gt.0.0_R8) hwxm_idt(itrg(j))=t/p
c*** xwhm: exponential decay length calculated from half-width
c        print '(a,i4,1p,2e12.3)','m,wo,wid_par(m,j):',
c                                  m,wo,wid_par(m,itrg(j)) !###
        do i=m,ny-1 !{
          if(wid_par(i,itrg(j)).lt.0.5_R8*wo) then !{
            d=(0.5_R8*wo-wid_par(i,itrg(j)))/
     .                  (wid_par(i-1,itrg(j))-wid_par(i,itrg(j)))
            hwhm_idp(itrg(j))=d*xmp_hlp(i-1,itrg(j))+
     .        (1-d)*xmp_hlp(i,itrg(j))-xmp_hlp(m,itrg(j))
            xwhm_idp(itrg(j))=hwhm_idp(itrg(j))/log(2.0_R8)
c!###{
c            print '(a,i4,1p,5e12.3)',
c     ,        'i,d,wid(i-1),wid(i),xmp(i-1),xmp(i):',i,d,
c     ,        wid_par(i-1,itrg(j)),wid_par(i,itrg(j)),
c     ,        xmp_hlp(i-1,itrg(j)),xmp_hlp(i,itrg(j))
c!###}
            exit
          end if !}
        end do !}
c        print '(a,i4,1p,2e12.3)','n,p,wid_trg(n,j):',
c                                  n,p,wid_trg(n,itrg(j)) !###
        do i=n,ny-1 !{
          if(wid_trg(i,itrg(j)).lt.0.5_R8*p) then !{
            d=(0.5_R8*p-wid_trg(i,itrg(j)))/
     .        (wid_trg(i-1,itrg(j))-wid_trg(i,itrg(j)))
            hwhm_idt(itrg(j))=d*xmp_hlp(i-1,itrg(j))+
     .        (1.0_R8-d)*xmp_hlp(i,itrg(j))-xmp_hlp(n,itrg(j))
            xwhm_idt(itrg(j))=hwhm_idt(itrg(j))/log(2.0_R8)
c!###{
c            print '(a,i4,1p,5e12.3)',
c     ,        'i,d,wid(i-1),wid(i),xmp(i-1),xmp(i):',i,d,
c     ,        wid_trg(i-1,itrg(j)),wid_trg(i,itrg(j)),
c     ,        xmp_hlp(i-1,itrg(j)),xmp_hlp(i,itrg(j))
c!###}
            exit
          end if !}
        end do !}
!###{
        if(fwht_tt(itrg(j)).lt.0.0_R8 .or.
     .     hwht_tt(itrg(j)).lt.0.0_R8 .or.
     .     hwhm_tt(itrg(j)).lt.0.0_R8) then !{
          print *,'!!!! NEGATIVE POWER WIDTH - target ',itrg(j),
     ,            ' time ',tim
          print '(1x,a,3(f10.5,4x))','fwht,hwht,hwhm = ',
     ,              fwht_tt(itrg(j)),hwht_tt(itrg(j)),hwhm_tt(itrg(j))
          print *,'xtrg:'
          print '(1x,i3,10f10.5/(4x,10f10.5))',lseppos(itrg(j))+1,
     ,                                        (xtrg(i,itrg(j)),i=0,ny-1)
          print *,'xtg_hlp:',xtg_hlp(-1,itrg(j))
          print '(1x,i3,10f10.5/(4x,10f10.5))',lseppos(itrg(j))+1,
     ,                                     (xtg_hlp(i,itrg(j)),i=0,ny-1)
          print *,'xmc_hlp:'
          print '(1x,i3,10f10.5/(4x,10f10.5))',lseppos(itrg(j))+1,
     ,                                     (xmc_hlp(i,itrg(j)),i=0,ny-1)
          print *,'m_hlp,wtttrg:',m_hlp(itrg(j))
          print '(1x,i3,1p,10e10.2/(4x,10e10.2))',lseppos(itrg(j))+1,
     ,                                     (wtttrg(i,itrg(j)),i=0,ny-1)
          print '(5a4,a7,9a12)','j','imx','ia','ib','ic','a','b','d',
     ,      'wpk','pwr','p','ha','hb','hc'
          print '(5i4,1p,9e12.4)',itrg(j),m_hlp(itrg(j)),ia,ib,ic,a,b,d,
     ,      wpk_tt(itrg(j)),pwr_tt(itrg(j)),p,ha,hb,hc
        end if !}
!###}
      end do !}
      if(ntrgts.gt.0) hlp_trg(ntrgts)=trg_short(ntrgsx)
      if(.not.LSN) then
        hlp_trg(1)=trg_short(ntrgsx)
        if(ntrgts.gt.0) hlp_trg(ntrgts)=trg_short(1)
      end if
      fn_hlp=fsumm(4)
      tl=' '
      open(lun,file=fn_hlp)
      write(lun,740) 'Power integrals ','and widths',trim(label)
      write(lun,750) 'Divertors are marked "il", "ol" ,"iu", "ou"'
      write(lun,750) trim(hdr_2(1)),' : integrated radiation power ',
     ,  'to the target'
      write(lun,750) trim(hdr_2(2)),' : integrated total power ',
     ,  'to the target'
      write(lun,750) trim(hdr_2(3)),' : FWHM of total power ',
     ,  'at the target'
      write(lun,750) trim(hdr_2(4)),' : double-sided exponential ',
     ,  'width of total power at the target'
      write(lun,750) trim(hdr_2(5)),' : HWHM of total power ',
     ,  'at the target'
      write(lun,750) trim(hdr_2(6)),' : the same mapped to the OMP'
      write(lun,750) trim(hdr_2(7)),' : exponential width of the SOL ',
     ,  'side of total power at the target'
      write(lun,750) trim(hdr_2(8)),' : the same mapped to the OMP'
      if(lrdlcr) then !{
        write(lun,750) trim(hdr_3(0)),' : extra power radiated from ',
     ,    'the core sources'
        write(lun,750) trim(hdr_3(1)),' : integrated radiation power ',
     ,    'to the target from the core sources'
        write(lun,750) trim(hdr_3(2)),' : fraction of the ',
     ,    'core radiation that reaches the target'
      end if !}
      write(lun,750) trim(hdr_4(1)),' : exponential width of the max. ',
     ,  'parallel flux mapped to the OMP'
      write(lun,750) trim(hdr_4(2)),' : exponential width of the max. ',
     ,  'target flux mapped to the OMP'
      write(lun,750) trim(hdr_4(3)),' : HWHM of the max. ',
     ,  'parallel flux mapped to the OMP'
      write(lun,750) trim(hdr_4(4)),' : HWHM of the max. ',
     ,  'target flux mapped to the OMP'
      write(lun,750) trim(hdr_4(5)),' : exponential decay length of ',
     ,  'the max. parallel flux, ',
     ,  'calculated from HM and mapped to the OMP'
      write(lun,750) trim(hdr_4(6)),' : exponential decay length of ',
     ,  'the max. target flux, calculated from HM and mapped to the OMP'
      write(lun,750) trim(hdr_4(7)),' : exponential decay length of ',
     ,  'the real target flux, ',
     ,  'calculated from HM and mapped to the OMP'
      n=(n_intgrs+n_trgprj)*ntrgs
      if(lrdlcr) n=n+n_lrdlcr*ntrgs+1
      write(fmt_peak,'(a,i3,a)') '(a6,a1,',n,'(a1,a))'
      if(lrdlcr) then !{
        write(lun,fmt_peak) '# time',tab,
     -   ((tab,hdr_2(i)//'_'//hlp_trg(itrg(j)),i=1,n_intgrs),j=1,ntrgs),
     ,   tab,hdr_3(0),
     ,   ((tab,hdr_3(i)//'_'//hlp_trg(itrg(j)),i=1,n_lrdlcr),j=1,ntrgs),
     ,   ((tab,hdr_4(i)//'_'//hlp_trg(itrg(j)),i=1,n_trgprj),j=1,ntrgs)
        write(fmt_peak,705) n
        write(lun,fmt_peak) tim,
     ,    (tab, pwr_rd(itrg(j)),tab, pwr_tt(itrg(j)),
     ,     tab,fwht_tt(itrg(j)),tab,fwxt_tt(itrg(j)),
     ,     tab,hwht_tt(itrg(j)),tab,hwhm_tt(itrg(j)),
     ,     tab,hwxt_tt(itrg(j)),tab,hwxm_tt(itrg(j)),j=1,ntrgs),
     ,     tab,1.e6*add_rad,
     ,    (tab,  pwr_cr(itrg(j)),tab,  frc_cr(itrg(j)),j=1,ntrgs),
     ,    (tab,hwxm_idp(itrg(j)),tab,hwxm_idt(itrg(j)),
     ,     tab,hwhm_idp(itrg(j)),tab,hwhm_idt(itrg(j)),
     ,     tab,xwhm_idp(itrg(j)),tab,xwhm_idt(itrg(j)),
     ,     tab,xwhm_tt(itrg(j)),j=1,ntrgs)
      else !}{
        write(lun,fmt_peak) '# time',tab,
     -   ((tab,hdr_2(i)//'_'//hlp_trg(itrg(j)),i=1,n_intgrs),j=1,ntrgs),
     ,   ((tab,hdr_4(i)//'_'//hlp_trg(itrg(j)),i=1,n_trgprj),j=1,ntrgs)
        write(fmt_peak,705) n
  705   format('(1p,e12.5,',i3,'(a1,e12.5))')
        write(lun,fmt_peak) tim,
     ,    (tab, pwr_rd(itrg(j)),tab, pwr_tt(itrg(j)),
     ,     tab,fwht_tt(itrg(j)),tab,fwxt_tt(itrg(j)),
     ,     tab,hwht_tt(itrg(j)),tab,hwhm_tt(itrg(j)),
     ,     tab,hwxt_tt(itrg(j)),tab,hwxm_tt(itrg(j)),j=1,ntrgs),
     ,    (tab,hwxm_idp(itrg(j)),tab,hwxm_idt(itrg(j)),
     ,     tab,hwhm_idp(itrg(j)),tab,hwhm_idt(itrg(j)),
     ,     tab,xwhm_idp(itrg(j)),tab,xwhm_idt(itrg(j)),
     ,     tab,xwhm_tt(itrg(j)),j=1,ntrgs)
      end if !}
      close (lun)
c----------------------------------------------------------------------
c*** Calculate the particle fluxes and erosion - if present
c     carlis eb2atm contains the indices of the B2 atoms mapped to Eirene
c            so if in Eirene is [D He Be Ne] and in B2 [D Ne He Be] then
c            eb2atm will be [1 4 2 3]
#ifndef B25_EIRENE
      do ia=1,nspecies
        TEXTA(ia)=species_list(ia)
      end do
#else
      do ia=1,nspecies !initialize
        eb2atm(ia)=ia
      end do
      do is=1,nspecies !species over B2
        do ia=1,natm !species over Eirene
          if (streql(TEXTA(ia),species_list(is))) then
           eb2atm(is)=ia
          end if
        end do
      end do
#endif
      write(hlp_frm,'(a,i2,a)') '(a,',nspecies,'(1x,i3))'
      write(6,hlp_frm) 'eb2spcr',(eb2spcr(is),is=1,nspecies)
      write(hlp_frm,'(a,i2,a)') '(a,',nspecies,'(1x,a2))'
      write(6,hlp_frm) 'b2ptrgl: list of species according to B2: ',
     ,                 (species_list(is),is=1,nspecies)
#ifdef B25_EIRENE
      write(hlp_frm,'(a,i2,a)') '(a,',natmi,'(1x,a8))'
      write(6,hlp_frm) 'List of atoms according to Eirene: ',
     .                 (TEXTA(is),is=1,natmi)
      write(hlp_frm,'(a,i2,a)') '(a,',nspecies,'(1x,a2))'
      write(6,hlp_frm) 'Ordering of B2 atoms according to Eirene: ',
     ,                 (atoms_list(eb2atm(i)),i=1,nspecies)
#endif
c*** Calculate the particle fluxes

c      print *,'Particle fluxes...' !###
      do j=1,ntrgsx !{
        do iy=0,ny-1 !{
          do ia=1,nclm_na
            sja(iy,j,ia)=0.0_R8
          end do
          do is=1,nspecies
            sjp(iy,j,is)=0.0_R8
            sjt(iy,j,is)=0.0_R8
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
        do ia=1,nclm_na
          sja_w(j,ia)=0.0_R8
          sja_t(j,ia)=0.0_R8
        end do
        do is=1,nspecies
          sjp_w(j,is)=0.0_R8
          sjt_w(j,is)=0.0_R8
          sjp_t(j,is)=0.0_R8
        end do
        do imol=1,nmoli
          sjm_w(j,imol)=0.0_R8
          sjm_t(j,imol)=0.0_R8
        end do
      end do !}

c      do is=0,ns-1
c        print *,'is,am(is)',is,am(is)
c      end do
c      do iss=1,ntns
c        if(ktns(iss).ne.1) cycle
c        js=ltns(iss)+nlim
c        do is=1,nclm_na
c          print *,'is,js,wldna(js,is,0)',is,js,wldna(js,is,0)
c        end do
c      end do
c      do is=1,ns-nspecies
c        print *,'is,wldpp(js,is,0)',is,wldpp(js,is,0)
c      end do
      u=-1.0_R8
      l=1
      do j=1,2*max(1,nncut) !{ Do not use ntrgts here to allow for sign flips
        if (itrg(j).ne.0) then !{
          k=ltrgpos(itrg(j))
c          print *,'j,ltrgpos = ',itrg(j),k  !###
          do iy=lytg,ny-1 !{
            do ia=1,nspecies !loop on species according to B2.5
              i1=eb2spcr(ia)+1
              i2=eb2spcr(ia)+nfluids(ia)
              do is=i1,i2
                sjp(iy,itrg(j),ia)=sjp(iy,itrg(j),ia)+
     .                             u*fna(k+1-l,iy,0,is)
              end do
            end do
            if (ft44) then
              do ia=1,natmi !loop on atoms according to Eirene
                spa(iy,itrg(j),ia)=
     .           dab2(k+l+1,iy+1,ia,1)*tab2(k+l+1,iy+1,ia,1)
              end do
              do imol=1,nmoli
                spm(iy,itrg(j),imol)=dmb2(k+l+1,iy+1,imol,1)*
     *                               tmb2(k+l+1,iy+1,imol,1)
              end do
              do iss=1,nsts !{
                if (eirdiag_nds_typ(iss).ne.2) cycle
                if (eirdiag_nds_srf(iss).ne.ltrgpos(itrg(j))+1) cycle
                if (eirdiag_nds_start(iss).gt.(iy+1)) cycle
                if (eirdiag_nds_end(iss).lt.(iy+1)) cycle
                IND=eirdiag_nds_ind(ISS)
                IAS=IND+1-(eirdiag_nds_start(iss)-1)
                do ia=1,natmi !loop on atoms according to Eirene
                  sja(iy,itrg(j),ia)=
     .              wldna_res(ia,iy+IAS)*sarea_res(iy+IAS)
                end do
                do imol=1,nmoli
                  sjm(iy,itrg(j),imol)=wldnm_res(imol,iy+IAS)*
     .                         sarea_res(iy+IAS)
                end do
              end do !}
            else
              do is=1,nspecies !loop on species according to B2.5
                sja(iy,itrg(j),is)=u*fna(k+1-l,iy,0,eb2spcr(is))
                spa(iy,itrg(j),is)=na(k+l,iy,eb2spcr(is))*ti(k+l,iy)
              end do
              do imol=1,nmoli
                sjm(iy,itrg(j),imol)=0.0_R8
                spm(iy,itrg(j),imol)=0.0_R8
              end do
            endif
            do is=1,nspecies
              sjp_t(itrg(j),is)=sjp_t(itrg(j),is)+sjp(iy,itrg(j),is)
            end do
            do ia=1,nclm_na
              sja_t(itrg(j),ia)=sja_t(itrg(j),ia)+sja(iy,itrg(j),ia)
            end do
            do imol=1,nmoli
              sjm_t(itrg(j),imol)=
     .         sjm_t(itrg(j),imol)+sjm(iy,itrg(j),imol)
            end do
          end do !}
        end if !}
        l=1-l
        u=-u
      end do !}
c!###{
c      print '(6(3x,a7,2x)/1p(6e12.4))','sjha_t','sjhea_t',
c     ,    'sjca_t','sjhp_t','sjhep_t','sjcp_t',
c     ,    ( sjha_t(itrg(i)),sjhea_t(itrg(i)),
c     ,      sjca_t(itrg(i)), sjhp_t(itrg(i)),
c     ,     sjhep_t(itrg(i)), sjcp_t(itrg(i)),i=1,ntrgts)
c      print '(2a5,6(3x,a7,2x))','j','js',                        !###
c     ,      'sjha_w','sjhea_w','sjca_w','sjhp_w','sjhep_w','sjcp_w'
c!###}
c
c      print *,'size wldpp(,,)',size(wldpp,1),size(wldpp,2),size(wldpp,3)
      if (ft44) then
        do j=1,ntrgts !{
          do iss=1,ntns
            if (ktns(iss).ne.itrg(j)) cycle
            js=ltns(iss)+nlim
            do ia=1,natmi !loop over Eirene atoms
              sja_w(itrg(j),ia)=sja_w(itrg(j),ia)+wldna(js,ia,0)
            end do
            do is=1,ns-nspecies !loop over Eirene background ions
              ia=latmscl(lkindp(is))
              if (ia.eq.0) cycle
              sjp_w(itrg(j),ia)=sjp_w(itrg(j),ia)+wldpp(js,is,0)
            end do
            do imol=1,nmoli
              sjm_w(itrg(j),imol)=sjm_w(itrg(j),imol)+wldnm(js,imol,0)
            end do
          end do
          do ia=1,natmi !loop over Eirene atoms
            if(sja_t(itrg(j),ia).gt.0.0_R8) then !{
              sja_t(itrg(j),ia)=sja_w(itrg(j),ia)/sja_t(itrg(j),ia)
            else !}{
              sja_t(itrg(j),ia)=0.0_R8
            end if !}
            do iy=lytg,ny-1 !{
              sja(iy,itrg(j),ia)=sja(iy,itrg(j),ia)*sja_t(itrg(j),ia)
              sjt(iy,itrg(j),latmscl(ia))=
     .                 sjt(iy,itrg(j),latmscl(ia))+sjp(iy,itrg(j),ia)
            end do !}
          end do
          do is=1,nspecies !loop over B2.5 species
            sjt_w(itrg(j),is)=sja_w(itrg(j),is)+sjp_w(itrg(j),is)
            if(sjp_t(itrg(j),is).gt.0.0_R8) then !{
              sjp_t(itrg(j),is)=sjp_w(itrg(j),is)/sjp_t(itrg(j),is)
            else !}{
              sjp_t(itrg(j),is)=0.0_R8
            end if !}
            do iy=lytg,ny-1 !{
              sjp(iy,itrg(j),is)=sjp(iy,itrg(j),is)*sjp_t(itrg(j),is)
              sjt(iy,itrg(j),is)=sjt(iy,itrg(j),is)+sjp(iy,itrg(j),is)
            end do !}
          end do
          do imol=1,nmoli
            if(sjm_t(itrg(j),imol).gt.0.0_R8) then !{
              sjm_t(itrg(j),imol)=
     .         sjm_w(itrg(j),imol)/sjm_t(itrg(j),imol)
            else !}{
              sjm_t(itrg(j),imol)=0.0_R8
            end if !}
            do iy=lytg,ny-1 !{
              sjm(iy,itrg(j),imol)=
     .         sjm(iy,itrg(j),imol)*sjm_t(itrg(j),imol)
            end do !}
          end do
        end do !}
      end if
c      write(*,'(1x,a,1p,16e12.3)') 'Scaling factors for atoms : ',
c     ,           ((sja_t(itrg(i),is),i=1,ntrgts),is=1,natmi)
c      write(*,'(1x,a,1p,16e12.3)') 'Scaling factors for species: ',
c     ,           ((sjp_t(itrg(i),is),i=1,ntrgts),is=1,nspecies)

c*** Calculate sputtering
      l=1
      do_sput=.false.
      if (ntargsp.gt.0.and.ft44) then !sputtering happening
        do j=1,ntrgts
c           Check if sputtering is happening and which material is being sputtered
c           and its atomic mass and density
c carlis A variable is needed that tells for each target or surface which is the sputtered atom.
c        Temporary solution for this routine is to use the bulk track_species
          do ia=1,n_srfs_m
            if (streql(mtrtrgsg(itrg(j)),srfs_m(ia))) then
              spt_mat(itrg(j))=ia
            end if
          end do
          sptlis(itrg(j))=0
          do iy=lytg,ny-1
            do iss=1,nsts
              if (eirdiag_nds_typ(iss).ne.2) cycle
              if (eirdiag_nds_srf(iss).ne.ltrgpos(itrg(j))+1) cycle
              if (eirdiag_nds_start(iss).gt.(iy+1)) cycle
              if (eirdiag_nds_end(iss).lt.(iy+1)) cycle
              IND=eirdiag_nds_ind(ISS)
              IAS=IND+1-(eirdiag_nds_start(iss)-1)
              do ia=1,natmi
                if (wldspta_res(ia,iy+IAS).gt.0.0_R8) then
                  sptlis(itrg(j))=ia
                  sjs(iy,itrg(j))=wldspta_res(ia,iy+IAS)*
     .                      sarea_res(iy+IAS)
                  sjs_t(itrg(j))=sjs_t(itrg(j))+sjs(iy,itrg(j))
                end if
              end do
              do imol=1,nmoli
                if (wldsptm_res(imol,iy+IAS).gt.0.0_R8) then
                  sptlis(itrg(j))=imol+natmi
                  sjs(iy,itrg(j))=wldsptm_res(imol,iy+IAS)*
     .                      sarea_res(iy+IAS)
                  sjs_t(itrg(j))=sjs_t(itrg(j))+sjs(iy,itrg(j))
                end if
              end do
            end do
          end do

          if(sptlis(itrg(j)).ne.0) then
            do iss=1,ntns
              if (ktns(iss).ne.itrg(j)) cycle
              js=ltns(iss)+nlim
              if (sptlis(itrg(j)).le.natmi) then
                sjs_w(itrg(j))=sjs_w(itrg(j))+
     .                         wldspta(js,sptlis(itrg(j)),0)
              else
                sjs_w(itrg(j))=sjs_w(itrg(j))+
     .                         wldsptm(js,sptlis(itrg(j))-natmi,0)
              end if
            end do
            if(sjs_t(itrg(j)).gt.0.0_R8) then
              sjs_t(itrg(j))=sjs_w(itrg(j))/sjs_t(itrg(j))
            else
              sjs_t(itrg(j))=0.0_R8
            end if
          else
            sjs_t(itrg(j))=0.0_R8
          end if

          if (sjs_t(itrg(j)).gt.0.0_R8) then
            do_sput(itrg(j))=.true.
            do iy=lytg,ny-1!{
              sjs(iy,itrg(j))=sjs_t(itrg(j))*sjs(iy,itrg(j))
c carlis The erosion of the wall is evaluated using the wall material and not the emitted atom/molecule
c carlis sje [mm/y] = ers_factor[s/y*mol/g*m2/cm2*mm/cm]*srfs_a[g/mol]*sjs[1/s]/srfs_d[g/cm3]/s[m2]
              sje(iy,itrg(j))=
     .          ers_factor*srfs_a(spt_mat(itrg(j)))*sjs(iy,itrg(j))/
     .                     srfs_d(spt_mat(itrg(j)))/s(iy,itrg(j))
              if (sptlis(itrg(j)).le.natmi) then
                sptmass=am(eb2atcr(sptlis(itrg(j))))
              else
                imol = sptlis(itrg(j))-natmi
                sptmass=0.0_R8
                do ia = 1, natmi
                  if (mlcmp(ia,imol).gt.0) then
                    sptmass=sptmass+mlcmp(ia,imol)*am(eb2atcr(ia))
                  end if
                end do
              end if
c carlis The emitted amount in g/y is evaluated using the sputtered atom or molecule
c carlis sjs (only in output) [g/y]=sjs[1/s]*srfs_a[g/mol]*nsecyr[s/y]/avogr[1/mol]
              sjs(iy,itrg(j))=sjs(iy,itrg(j))*sptmass*nsecyr/avogr
              sjs_w(itrg(j))=sjs_w(itrg(j))+sjs(iy,itrg(j))
            end do
            if (sptlis(itrg(j)).le.natmi) then
              ssmol=texta(sptlis(itrg(j)))
            else
              ssmol=textmn(sptlis(itrg(j))-natmi-1)
              ssmol=adjustl(ssmol)
            end if
            print *,'target material = ', mtrtrgsg(itrg(j))
            print *,'sputtered atom from target = ',
     ,               ssmol
            write(*,'(1x,a,1p,16e12.3)')
     .         'Scaling factors for '//mtrtrgsg(itrg(j))//' Sp : ',
     ,          (sjs_t(itrg(i)),i=1,ntrgts)
          else ! no sputtering
            print *,'b2ptrgl, zero sputtering flux for target j=',
     ,                                                   itrg(j)
            sjs_w(itrg(j))=0.0_R8
          end if
        end do !loop on targets
      else
        print '(4a)', 'b2ptrgl, no sputtering calculated:',
     ,     ' track_species=0 or no ',fort_lc,'44 file available'
      end if

      l=1
      do ia=1,nclm_na
        fpt_str(1+3*(ia-1))=trim(TEXTA(ia))//'_flux_atm'
        fpt_str(2+3*(ia-1))=trim(TEXTA(ia))//'_temp_atm'
        fpt_str(3+3*(ia-1))=trim(TEXTA(ia))//'_pres_atm'
      end do
      do is=1,nspecies
        fpt_str(is+3*nclm_na)=trim(species_list(is))//'_flux_ion'
      end do
      do k=1,ntrgts
        j=itrg(k)
        if (do_sput(j)) then
          if (sptlis(j).le.natmi) then
            spt_str(1,j)=trim(TEXTA(sptlis(j)))//'_sput'
          else
            ssmol=textmn(sptlis(j)-natmi-1)
            ssmol=adjustl(ssmol)
            spt_str(1,j)=trim(ssmol)//'_sput'
          end if
          spt_str(2,j)=trim(mtrtrgsg(j))//'_ersn'
        end if
      end do
      do imol=1,nmoli
        ssmol=textmn(imol-1)
        ssmol=adjustl(ssmol)
        mol_str(1,imol)=trim(ssmol)//'_flux_mol'
        mol_str(2,imol)=trim(ssmol)//'_pres_mol'
      end do
      do j=1,2*max(nncut,1) !{ Do not use ntrgts here to allow for sign flips
        if (itrg(j).gt.0) then !{
          nx_hlp=ltrgpos(itrg(j))+l
          jj=itrg(j)
          if(j.eq.ntrgts.and..not.(ntrgts.eq.1.and.itrg(1).eq.1))
     >     jj=ntrgsx
          if(.not.LSN) jj=ntrgsx+1-jj
          open(lun,file=ftrg(jj,2))
          if (do_sput(itrg(j))) then
            nclm_er=2+3*nclm_na+nspecies+2*ntrgts+2*nclm_nm
          else
            nclm_er=2+3*nclm_na+nspecies+2*nclm_nm
          endif
          write(fmt_peak,710) nclm_er
          write(lun,740) 'Target erosion - ',trg_label(jj),trim(label)
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
#ifdef B25_EIRENE
          write(lun,740) 'XX_flux_mol',tab,
     .     'incident flux of molecules XX [1/m^2/s]'
          write(lun,740) 'XX_pres_mol',tab,
     .     'pressure of molecules XX [Pa]'
#endif
          write(lun,740) 'XX_sput',tab,
     .     'gross sputtering of material XX [g/y]'
          write(lun,740) 'XX_ersn',tab,
     .     'gross erosion of material XX [mm/y]'
          write(lun,740) 'xMP',tab,'radial coordinate in midplane '//
     .     '(IMP for inner target(s), OMP for outer(s)) [m]'
          if (do_sput(itrg(j))) then
            write(lun,fmt_peak) '#','x',
     ,           (fpt_str(ia),ia=1,3*nclm_na+nspecies),
     ,           (mol_str(:,ia),ia=1,nclm_nm),spt_str(:,itrg(j)),'xMP'
            write(fmt_peak,720) nclm_er
          else
            write(lun,fmt_peak) '#','x',
     ,           (fpt_str(ia),ia=1,3*nclm_na+nspecies),
     ,           (mol_str(:,ia),ia=1,nclm_nm),'xMP'
            write(fmt_peak,720) nclm_er
          end if
          if(ntrgts.eq.4) then !{
            k=itrg(j)
          else !}{
            k=itrg(ntrgts)
          end if !}
          do iy=lytg,ny-1 !{
            do ia=1,nclm_na
              fpt_out(1+3*(ia-1))=sja(iy,itrg(j),ia)/s(iy,itrg(j))
#ifdef B25_EIRENE
              fpt_out(2+3*(ia-1))=tab2(nx_hlp+1,iy+1,ia,1)/ev
#else
              fpt_out(2+3*(ia-1))=ti(nx_hlp,iy)/ev
#endif
              fpt_out(3+3*(ia-1))=spa(iy,itrg(j),ia)
            end do
            do is=1,nspecies
              fpt_out(is+3*nclm_na)=sjp(iy,itrg(j),is)/s(iy,itrg(j))
            end do
            do imol=1,nmoli
              mol_out(1+2*(imol-1))=sjm(iy,itrg(j),imol)/s(iy,itrg(j))
              mol_out(2+2*(imol-1))=spm(iy,itrg(j),imol)
            end do
            if (do_sput(itrg(j))) then
              spt_out(1)=sjs(iy,itrg(j))
              spt_out(2)=sje(iy,itrg(j))
              write(lun,fmt_peak) xtrg(iy,itrg(j)),
     ,         (fpt_out(ia),ia=1,3*nclm_na+nspecies),
     ,         (mol_out(ia),ia=1,2*nclm_nm),
     ,          spt_out,xmp_hlp(iy,k)
            else
              write(lun,fmt_peak) xtrg(iy,itrg(j)),
     ,         (fpt_out(ia),ia=1,3*nclm_na+nspecies),
     ,         (mol_out(ia),ia=1,2*nclm_nm),
     ,          xmp_hlp(iy,k)
            end if
          end do !}
          close (lun)
        end if !}
        l=1-l
      end do !}
c----------------------------------------------------------------------
c*** Produce a summary file for target sputtering
c
      do ia=1,nclm_na
        sja_w(0,ia)=0.0_R8
      end do
      do is=1,nspecies
        sjp_w(0,is)=0.0_R8
        sjt_w(0,is)=0.0_R8
      end do
      do imol=1,nmoli
        sjm_w(0,imol)=0.0_R8
      end do
      sjs_w(0)=0.0_R8
      do j=1,ntrgts
        do ia=1,nclm_na
          sja_w(0,ia)=sja_w(0,ia)+sja_w(itrg(j),ia)
        end do
        do is=1,nspecies
          sjp_w(0,is)=sjp_w(0,is)+sjp_w(itrg(j),is)
          sjt_w(0,is)=sjt_w(0,is)+sjt_w(itrg(j),is)
        end do
        do imol=1,nmoli
          sjm_w(0,imol)=sjm_w(0,imol)+sjm_w(itrg(j),imol)
        end do
      end do !}
      k=1
      if(ntrgts.eq.1) k=ntrgsx
      if(ntrgts.eq.2) k=ntrgsx-1
      open(lun,file=fsumm(2))
      write(lun,'(2a)') 'Target erosion summary: ',trim(label)
      write(lun,*)
      if (LSN .and. itrg(1).eq.1) then
        write(lun,'(12x,3(a1,a12))')
     -                   (tab,trg_label(i),i=1,ntrgsx,k),tab,'total   '
      else
        write(lun,'(12x,3(a1,a12))')
     -                   (tab,trg_label(i),i=ntrgsx,1,-k),tab,'total   '
      endif
      write(fmt_peak,730) ntrgts+1
  730 format('(1p,1x,a12,',i3,'(a1,e12.4))')
      write(lun,*) 'Particle fluxes [1/s]'
      do is=1,nspecies
        write(lun,fmt_peak) trim(species_list(is))//' ions on',
     ,                       (tab,sjp_w(itrg(i),is),i=1,ntrgts),
     ,                        tab,sjp_w(0,is)
      end do
      do ia=1,nclm_na
        write(lun,fmt_peak) trim(TEXTA(ia))//' atoms on',
     ,                       (tab,sja_w(itrg(i),ia),i=1,ntrgts),
     ,                        tab,sja_w(0,ia)
      end do
      do is=1,nspecies
        write(lun,fmt_peak) trim(species_list(is))//' total on',
     ,                       (tab,sjt_w(itrg(i),is),i=1,ntrgts),
     ,                        tab,sjt_w(0,is)
      end do
      if (use_eirene.ne.0.and.nmoli.gt.0) then
        write(lun,*)
        write(lun,*) 'Molecule fluxes [1/s]'
        do imol=1,nmoli
          ssmol=textmn(imol-1)
          ssmol=adjustl(ssmol)
          write(lun,fmt_peak) trim(ssmol)//' mols on',
     ,               (tab,sjm_w(itrg(i),imol),i=1,ntrgts),
     ,                tab,sjm_w(0,imol)
        end do
      end if
      if (ntargsp.gt.0.and.ft44) then
        k=0
        do j=1,ntrgts
          sjs_w(0)=sjs_w(0)+sjs_w(itrg(j))
          if(k.eq.0.and..not.streql(mtrtrgsg(itrg(j)),' ')) k=itrg(j)
        end do
        if (k.ne.0) then
          write(lun,*)
          write(lun,*) 'Total sputtering [g/y]'
          write(lun,fmt_peak) trim(mtrtrgsg(k)),
     ,                       (tab,sjs_w(itrg(i)),i=1,ntrgts),
     ,                        tab,sjs_w(0)
        end if
        write(lun,*)
        write(lun,*) 'Erosion [mm/y]'
        ersn_pm=0.0_R8
        ersn_am=0.0_R8
        ss_m=0.0_R8
        do j=1,ntrgts !{
          if (do_sput(itrg(j))) then
            ss(itrg(j))=0.0_R8
            ersn_p(itrg(j))=0.0_R8
            do iy=lytg,ny-1 !{
              sn_t=0.0_R8
              ss(itrg(j))=ss(itrg(j))+s(iy,itrg(j))
              sn_t=sn_t+sje(iy,itrg(j))
              if(sn_t .gt. ersn_p(itrg(j))) then !{
                ersn_p(itrg(j))=sn_t
                ersn_l(itrg(j))=xtrg(iy,itrg(j))
              end if !}
            end do !}
            ersn_a(itrg(j))=0.0_R8
            ersn_a(itrg(j))=ersn_a(itrg(j))+sjs_w(itrg(j))/ss(itrg(j))
            ersn_pm=max(ersn_p(itrg(j)),ersn_pm)
            ss_m=ss_m+ss(itrg(j))
            ersn_am=ersn_am+ersn_a(itrg(j))*ss(itrg(j))
          else
            ersn_a(itrg(j))=0.0_R8
            ersn_p(itrg(j))=0.0_R8
            ersn_l(itrg(j))=0.0_R8
          end if
        end do !}
        ersn_am=ersn_am/ss_m
        write(lun,fmt_peak) 'average',(tab,ersn_a(itrg(i)),i=1,ntrgts),
     .                                 tab,ersn_am
        write(lun,fmt_peak) 'peak',(tab,ersn_p(itrg(i)),i=1,ntrgts),
     .                              tab,ersn_pm
        write(lun,fmt_peak) 'peak pos.',(tab,ersn_l(itrg(i)),i=1,ntrgts)
        close (lun)
      end if

      ncall = ncall+1   !iyv 19.07.22
      return

c======================================================================
  710 format('(a1,',i3,'(1x,a11,2x))')
  720 format('(1x,1p,',i3,'(e12.5,2x))')
  740 format('# ',2a,' : ',2a)
  750 format('# ',4a)
c=======================================================================
      end subroutine b2ptrgl
